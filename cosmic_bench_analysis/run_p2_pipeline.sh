#!/bin/bash
# run_p2_pipeline.sh
#
# Driver for the full P2 (BASKET) cosmic-bench QA pipeline (stages 02..12) on a
# single long run registered in p2_qa_config.py. It reproduces the exact
# stage / spark-veto matrix that build_final_pdf.py expects: every veto-capable
# stage is run twice (once spark-vetoed, once not) so the report can show both.
#
# Usage:
#   ./run_p2_pipeline.sh [run_key]      # default: det1_long2
#
# This is for a *long run at a single working point*, NOT an HV scan. The HV
# scan uses stage 11 + build_hv_scan_pdf.py instead (see build_hv_scan_pdf.py).
#
# Stages 09 (pedestal QA) needs the per-FEU hits_root/ tree; if it is absent the
# stage is skipped with a warning rather than failing the whole pipeline.

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
RUN_KEY="${1:-det1_long2}"

if [[ -f "${SCRIPT_DIR}/../.venv/bin/activate" ]]; then
    # shellcheck disable=SC1091
    source "${SCRIPT_DIR}/../.venv/bin/activate"
fi

cd "${SCRIPT_DIR}"

echo "=================================================================="
echo " P2 QA pipeline  |  run_key = ${RUN_KEY}"
echo "=================================================================="

run () {   # run <stage.py> <args...>
    echo ""
    echo ">>> python3 $*"
    python3 "$@"
}

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
run 10_efficiency_map_sliding.py "${RUN_KEY}"

# -- 12 structural-bias validation (vetoed) ---------------------------------- #
run 12_validation.py "${RUN_KEY}"

echo ""
echo "=================================================================="
echo " Pipeline done for ${RUN_KEY}. Build the report with:"
echo "   python3 build_final_pdf.py ${RUN_KEY}"
echo "=================================================================="
