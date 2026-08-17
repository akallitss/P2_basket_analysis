#!/bin/bash
# run_stage19_pass.sh
#
# Waits for a running re-analysis batch to finish, wires stage 19 into
# run_p2_pipeline.sh so future runs include it, then runs a stage-19 pass over
# every long-run key so all detectors have a timing map on equal footing.
#
#   ./run_stage19_pass.sh [--wait-for BATCH.log]
#
# Why it waits: the batch invokes run_p2_pipeline.sh fresh per key, so editing
# that script mid-flight would give some runs stage 19 and others not -- the
# same silent inconsistency that produced the scan-vs-long-run mismatch.
#
# Scans (hv_scan / drift_scan / mesh_scan keys) are NOT included: their
# SUB_RUN is a pseudo-name covering many physical sub_runs, so stage 19 has no
# single combined_hits dir to read. Covering them needs a per-sub_run loop.

set -u
export P2_MINIMAL_CUTS=1
export P2_MINIMAL_EXCEPT=HOT_PAD_RATIO
export P2_HOT_PAD_RATIO=3.0

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${HERE}"

WAIT_LOG=""
[[ "${1:-}" == "--wait-for" ]] && WAIT_LOG="${2:-}"

if [[ -n "${WAIT_LOG}" ]]; then
    echo "waiting for batch to finish: ${WAIT_LOG}"
    while ! grep -q "^ done 2026" "${WAIT_LOG}" 2>/dev/null; do sleep 60; done
    echo "batch finished $(date '+%F %T')"
fi

# -- 1. wire stage 19 into the long-run pipeline (idempotent) ---------------- #
if grep -q "19_timing_surface" run_p2_pipeline.sh; then
    echo "stage 19 already wired into run_p2_pipeline.sh"
else
    python3 - <<'EOF'
p = 'run_p2_pipeline.sh'
s = open(p).read()
anchor = 'run 12_validation.py "${RUN_KEY}"'
add = anchor + '''

# Surface timing map: per-pad peak time + resolution, and (--continuous) a
# sliding map in the M3 reference frame. Must follow 06, whose ray_hit_miss
# list the continuous mode joins against.
run 19_timing_surface.py "${RUN_KEY}" --continuous'''
assert s.count(anchor) == 1, s.count(anchor)
open(p, 'w').write(s.replace(anchor, add, 1))
print('  wired stage 19 into run_p2_pipeline.sh')
EOF
fi

# -- 2. stage-19 pass over every long-run key -------------------------------- #
KEYS=(det1_long det1_long3 det1_long4 det1_initial1 det1_long5
      det2_long1 det3_initial1 det3_final1 det4_initial1 det4_long2)

LOGDIR=$(python3 -c 'import p2_qa_config as C; print(C.ANALYSIS_ROOT)')/_stage19_logs
mkdir -p "${LOGDIR}"
SUM="${LOGDIR}/summary.tsv"
: > "${SUM}"

echo "=================================================================="
echo " stage 19 pass — ${#KEYS[@]} long-run keys — $(date '+%F %T')"
echo "=================================================================="
ok=0; fail=0
for k in "${KEYS[@]}"; do
    t0=$(date +%s)
    printf '>>> %-16s ' "${k}"
    if python3 19_timing_surface.py "${k}" --continuous \
            > "${LOGDIR}/${k}.log" 2>&1; then
        st=OK; ok=$((ok+1))
    else
        st=FAILED; fail=$((fail+1))
    fi
    dt=$(( $(date +%s) - t0 ))
    echo "${st} in ${dt}s"
    [[ "${st}" == FAILED ]] && tail -3 "${LOGDIR}/${k}.log" | sed 's/^/      | /'
    # pull the headline numbers straight out of the log
    line=$(grep -E "genuine pad-to-pad|per-pad peak-time sigma" "${LOGDIR}/${k}.log" 2>/dev/null | tr '\n' ' ')
    printf '%s\t%s\t%ds\t%s\n' "${k}" "${st}" "${dt}" "${line}" >> "${SUM}"
done

echo ""
echo "=================================================================="
echo " stage 19 pass done $(date '+%F %T')   ok=${ok} failed=${fail}"
echo " logs: ${LOGDIR}"
echo "=================================================================="
cut -f1-3 "${SUM}" | column -t
