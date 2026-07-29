#!/usr/bin/env bash
#
# Build the single merged ROOT file of all sub-runs, for uRWELL + P2 tracking.
#
# Wraps merge_subruns.py with the environment it needs.  See that script's
# docstring for why hadd cannot be used here and what geventId is for.
#
# Usage:
#   ./build_merged_root.sh                        # drift_mesh_scan_1 + highstat_eff_1
#   ./build_merged_root.sh -n                     # dry run, list the sub-runs
#   ./build_merged_root.sh -r drift_mesh_scan_1   # one run only
#   ./build_merged_root.sh -f 1                   # uRWELLs only (FEU 1)
#   ./build_merged_root.sh -o /path/out.root
#
# Long job: expect roughly an hour per 10 GB of input.  Run it under nohup or
# in a tmux/screen session, not in a shell you will close.
#
set -u -o pipefail

unset PYTHONPATH                          # handoff §2 - not optional
PY=/local/home/banco/DAQ_Control_Dream_Beam/.venv/bin/python

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
OUT=/local/home/banco/P2_data/TB_July2026_H4/merged/all_subruns_hits.root

RUNS=()
EXTRA=()

while getopts ":r:o:f:nAah" opt; do
    case "$opt" in
        r) RUNS+=("$OPTARG") ;;
        o) OUT="$OPTARG" ;;
        f) EXTRA+=(--feu "$OPTARG") ;;
        n) EXTRA+=(--dry-run) ;;
        A) EXTRA+=(--include-incomplete) ;;
        a) EXTRA+=(--all-branches) ;;
        h) sed -n '2,17p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        \?) echo "unknown option -$OPTARG" >&2; exit 1 ;;
    esac
done
[ ${#RUNS[@]} -eq 0 ] && RUNS=(drift_mesh_scan_1 highstat_eff_1)

[ -x "$PY" ] || { echo "interpreter not found: $PY" >&2; exit 1; }

# refuse to start if the destination filesystem cannot hold the result
mkdir -p "$(dirname "$OUT")"
avail_gb=$(df -BG --output=avail "$(dirname "$OUT")" | tail -1 | tr -dc '0-9')
if [ "${avail_gb:-0}" -lt 40 ]; then
    echo "only ${avail_gb}G free on $(dirname "$OUT") - the merge needs room" >&2
    exit 1
fi

echo "runs   : ${RUNS[*]}"
echo "output : $OUT"
echo "free   : ${avail_gb}G"
echo

exec "$PY" "$HERE/merge_subruns.py" \
    --runs "${RUNS[@]}" \
    --out "$OUT" \
    "${EXTRA[@]}"
