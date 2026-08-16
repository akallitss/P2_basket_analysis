#!/bin/bash
# ---------------------------------------------------------------------------
# sweep_efficiency.sh -- uRWELL-track-referenced VMM efficiency, whole campaign.
#
#   ./sweep_efficiency.sh [--jobs N] [--dry-run] [--only run_32,run_61] [--force]
#
# Runs urw_vmm_efficiency.py on every sub_run whose streams are matched (a
# match summary on EOS with a usable lock), and stores the result next to the
# matches:
#
#   <eos>/vmm/matching/efficiency/vmm_eff_<run>_<sub>.json
#
# Same conventions as sweep_matching.sh: work in $TMPDIR, push with xrdcp,
# nothing of size on the AFS home, resumable.
# ---------------------------------------------------------------------------
set -uo pipefail

# NB: not BASE -- the LCG setup.sh defines BASE, and sourcing it
# after this line silently repointed URW_DET_DIR into /cvmfs.
EOSBASE=/eos/experiment/ntof/data/x17/p2_sps_july
EOS_URL=root://eospublic.cern.ch
MATCH=$EOSBASE/vmm/matching
OUT=$MATCH/efficiency
LCG=/cvmfs/sft.cern.ch/lcg/views/LCG_110/x86_64-el9-gcc13-opt/setup.sh
JOBS=4
DRY=0
FORCE=0
ONLY=""
MIN_FRAC=0.5
EXTRA=""

while [ $# -gt 0 ]; do
    case "$1" in
        --jobs) JOBS=$2; shift 2 ;;
        --dry-run) DRY=1; shift ;;
        --force) FORCE=1; shift ;;
        --only) ONLY=$2; shift 2 ;;
        --min-frac) MIN_FRAC=$2; shift 2 ;;
        --out-dir) OUT=$2; shift 2 ;;
        --extra) EXTRA=$2; shift 2 ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
done

CODE=$(cd "$(dirname "$0")" && pwd)
WORK=${TMPDIR:-/tmp}/p2_eff_sweep.$USER
mkdir -p "$WORK"

set +u
# shellcheck disable=SC1090
source "$LCG"
set -u
export PYTHONPATH=$CODE:${PYTHONPATH:-}
export SACLAY_MM_DIR=${SACLAY_MM_DIR:-$CODE}
export URW_DET_DIR=${URW_DET_DIR:-$EOSBASE/config/detectors/}
xrdfs "$EOS_URL" mkdir -p "$OUT" 2>/dev/null

# --- work list: matched sub_runs -------------------------------------------
LIST=$WORK/joblist.txt
python3 - "$MATCH" "$OUT" "$ONLY" "$FORCE" "$MIN_FRAC" > "$LIST" <<'PY'
import glob, json, os, sys
match, out, only, force, min_frac = sys.argv[1:6]
only = [x for x in only.split(',') if x]
for fn in sorted(glob.glob(f'{match}/match_*.json')):
    try:
        d = json.load(open(fn))
    except Exception:
        continue
    if only and d['run'] not in only:
        continue
    f = d.get('match_frac_covered')
    f = d.get('match_frac_dream', 0.0) if f is None else f
    if f < float(min_frac):
        continue
    if force == '0' and os.path.exists(f"{out}/vmm_eff_{d['run']}_{d['sub']}.json"):
        continue
    print(d['run'], d['sub'])
PY
n=$(wc -l < "$LIST")
echo "$n sub_run(s) to measure (jobs=$JOBS, force=$FORCE, min_frac=$MIN_FRAC)"
[ "$DRY" = 1 ] && { cat "$LIST"; exit 0; }
[ "$n" = 0 ] && exit 0

cat > "$WORK/one.sh" <<EOS
#!/bin/bash
set -uo pipefail
run=\$1; sub=\$2
tag=\${run}_\${sub}
d=$WORK/\$tag
mkdir -p "\$d"
python3 $CODE/urw_vmm_efficiency.py "\$run" "\$sub" --out "\$d" $EXTRA \
    > "\$d/log.txt" 2>&1
rc=\$?
f=\$d/vmm_eff_\$tag.json
if [ \$rc -ne 0 ] || [ ! -f "\$f" ]; then
    echo "FAIL(\$rc) \$tag  \$(tail -2 "\$d/log.txt" | tr '\n' ' ')"
    rm -rf "\$d"; exit \$rc
fi
xrdcp -f -s "\$f" "$EOS_URL/$OUT/\$(basename \$f)" || {
    echo "UPLOAD-FAIL \$tag"; exit 3; }
echo "OK \$tag  \$(python3 -c "
import json;d=json.load(open('\$f'))
print(' '.join('%s=%.3f'%(s['station'],s['efficiency']['value'])
               for s in d['stations'] if 'efficiency' in s))")"
rm -rf "\$d"
EOS
chmod +x "$WORK/one.sh"

xargs -a "$LIST" -n 2 -P "$JOBS" "$WORK/one.sh"
echo "EFF_SWEEP_DONE $(date -u +%H:%M:%S)"
