#!/bin/bash
# run_p2_pipeline.sh
#
# Driver for the full P2 (BASKET) cosmic-bench QA pipeline (stages 02..12) on a
# single long run registered in p2_qa_config.py. It reproduces the exact
# stage / spark-veto matrix that build_final_pdf.py expects: every veto-capable
# stage is run twice (once spark-vetoed, once not) so the report can show both.
#
# Usage:
#   ./run_p2_pipeline.sh [run_key]                # long run, stages 02..12
#   ./run_p2_pipeline.sh [run_key] --scan drift   # drift scan (stage 16) + PDF
#   ./run_p2_pipeline.sh [run_key] --scan mesh    # mesh scan (11+16) + PDF
#
# Default (no --scan) is a *long run at a single working point*: stages 02..12,
# then build_final_pdf.py. The --scan modes instead run the voltage-scan stage
# for that run key and build the matching scan report:
#   drift -> 16_drift_scan_efficiency.py --scan drift -> build_drift_scan_pdf.py
#   mesh  -> 11_hv_scan_efficiency.py + 16 --scan mesh -> build_hv_scan_pdf.py
#
# Stages 09 (pedestal QA) needs the per-FEU hits_root/ tree; if it is absent the
# stage is skipped with a warning rather than failing the whole pipeline.

set -e

# Re-exec the whole pipeline inside a cgroup memory cap: a runaway stage is
# then OOM-killed by the kernel instead of freezing the machine (14 GB RAM,
# no swap). MemoryHigh throttles first; MemoryMax kills.
if [[ -z "${P2_MEMCAPPED:-}" ]] && command -v systemd-run >/dev/null 2>&1; then
    exec env P2_MEMCAPPED=1 systemd-run --user --scope --quiet \
        -p MemoryHigh=7G -p MemoryMax=8G -- "$0" "$@"
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# -- args: [run_key] [--scan drift|mesh] ------------------------------------- #
RUN_KEY=""
SCAN=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --scan)   SCAN="$2"; shift 2 ;;
        --scan=*) SCAN="${1#*=}"; shift ;;
        -*)       echo "unknown option: $1" >&2; exit 2 ;;
        *)        RUN_KEY="$1"; shift ;;
    esac
done
RUN_KEY="${RUN_KEY:-det1_long2}"
if [[ -n "${SCAN}" && "${SCAN}" != "drift" && "${SCAN}" != "mesh" ]]; then
    echo "invalid --scan '${SCAN}' (use drift or mesh)" >&2; exit 2
fi

if [[ -f "${SCRIPT_DIR}/../.venv/bin/activate" ]]; then
    # shellcheck disable=SC1091
    source "${SCRIPT_DIR}/../.venv/bin/activate"
fi

cd "${SCRIPT_DIR}"

echo "=================================================================="
echo " P2 QA pipeline  |  run_key = ${RUN_KEY}${SCAN:+  |  ${SCAN} scan}"
echo "=================================================================="

run () {   # run <stage.py> <args...>
    echo ""
    echo ">>> python3 $*"
    python3 "$@"
}

# -- voltage-scan pipelines: scan stage + matching PDF, then done ------------- #
if [[ "${SCAN}" == "drift" ]]; then
    run 16_drift_scan_efficiency.py "${RUN_KEY}" --scan drift
    run build_drift_scan_pdf.py "${RUN_KEY}"
    echo ""
    echo "=================================================================="
    echo " Drift scan done for ${RUN_KEY}. Report: p2_<det>_drift_scan.pdf"
    echo "=================================================================="
    exit 0
elif [[ "${SCAN}" == "mesh" ]]; then
    run 11_hv_scan_efficiency.py "${RUN_KEY}"
    run 16_drift_scan_efficiency.py "${RUN_KEY}" --scan mesh
    run build_hv_scan_pdf.py "${RUN_KEY}"
    echo ""
    echo "=================================================================="
    echo " Mesh scan done for ${RUN_KEY}. Report: p2_<det>_hv_scan.pdf"
    echo "=================================================================="
    exit 0
fi

# -- 02 map validation: non-vetoed (+ strategy comparison) and vetoed --------- #
run 02_map_validation.py "${RUN_KEY}" --no-veto-sparks --compare-strategies
run 02_map_validation.py "${RUN_KEY}"

# -- 03 M3 alignment: non-vetoed and vetoed ---------------------------------- #
run 03_m3_alignment.py "${RUN_KEY}" --no-veto-sparks
run 03_m3_alignment.py "${RUN_KEY}"

# -- 04 M3 reference QA: reference tracker only (no veto variants) ------------ #
run 04_m3_reference_qa.py "${RUN_KEY}"

# -- 05 detector deep QA: non-vetoed and vetoed ------------------------------ #
run 05_detector_deep_qa.py "${RUN_KEY}" --no-veto-sparks
run 05_detector_deep_qa.py "${RUN_KEY}"

# -- 06 efficiency maps: non-vetoed and vetoed (vetoed saves pad_footprint) --- #
run 06_efficiency_maps.py "${RUN_KEY}" --no-veto-sparks
run 06_efficiency_maps.py "${RUN_KEY}"

# -- 07 HV spark QA: the spark analysis itself (single, dead-suffix only) ----- #
run 07_hv_spark_qa.py "${RUN_KEY}"

# -- 08 per-pad spark QA: residual sparks on veto-cleaned hits (vetoed only) -- #
run 08_pad_spark_qa.py "${RUN_KEY}"

# -- 09 pedestal QA: needs per-FEU hits_root/ (skip if absent) ---------------- #
HITS_ROOT="$(python3 - "${RUN_KEY}" <<'PY'
import sys, p2_qa_config as qa
print(qa.get_config(sys.argv[1]).hits_root_dir)
PY
)"
if [[ -d "${HITS_ROOT}" ]] && ls "${HITS_ROOT}"/*.root >/dev/null 2>&1; then
    run 09_pedestal_qa.py "${RUN_KEY}" --no-veto-sparks
else
    echo ""
    echo ">>> SKIP 09 pedestal QA: no hits_root/ files at ${HITS_ROOT}"
fi

# -- 10 sliding-window efficiency map (vetoed, writes into 06_efficiency/) ---- #
# Fine 5 mm kernel: resolves the ~6 mm big insulation pillars. --min 10 keeps
# the small kernel populated (~15-20 rays/kernel at cosmic statistics).
run 10_efficiency_map_sliding.py "${RUN_KEY}" --kernel 5 --grid 240 --min 10

# -- 12 structural-bias validation (vetoed) ---------------------------------- #
run 12_validation.py "${RUN_KEY}"

# Surface timing map: per-pad peak time + resolution, and (--continuous) a
# sliding map in the M3 reference frame. Must follow 06, whose ray_hit_miss
# list the continuous mode joins against.
run 19_timing_surface.py "${RUN_KEY}" --continuous

echo ""
echo "=================================================================="
echo " Pipeline done for ${RUN_KEY}. Build the report with:"
echo "   python3 build_final_pdf.py ${RUN_KEY}"
echo "=================================================================="
