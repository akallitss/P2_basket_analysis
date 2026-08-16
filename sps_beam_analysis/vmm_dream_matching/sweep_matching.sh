#!/bin/bash
# ---------------------------------------------------------------------------
# sweep_matching.sh -- match every VMM sub_run of the campaign to DREAM, once.
#
#   ./sweep_matching.sh [--jobs N] [--dry-run] [--only run_57,run_58] [--force]
#
# Run it on lxplus. For each (run, sub_run) that has BOTH a DREAM
# combined_hits_root and a VMM capture directory it produces, on EOS under
# vmm/matching/:
#
#   vmm_triggers_<run>_<sub>.npz   all VMM trigger-channel (VMM 0 ch 44) times
#   match_<run>_<sub>.json         per-spill fit + match summary
#   match_<run>_<sub>.npz          one row per DREAM event: matched?, t_vmm,
#                                  residual, index into the trigger array
#
# Everything is written to $TMPDIR first and pushed with xrdcp: the AFS home is
# 5 GB and a single match npz is tens of MB, so nothing of size may land there.
# The sweep is resumable -- a sub_run whose match json is already on EOS is
# skipped unless --force.
# ---------------------------------------------------------------------------
set -uo pipefail

BASE=/eos/experiment/ntof/data/x17/p2_sps_july
EOS_URL=root://eospublic.cern.ch
OUT=$BASE/vmm/matching
LCG=/cvmfs/sft.cern.ch/lcg/views/LCG_110/x86_64-el9-gcc13-opt/setup.sh
JOBS=4
DRY=0
FORCE=0
ONLY=""
# vmm_decode.py lives in the online-analysis repo, not this one; the pcapng
# path needs it on PYTHONPATH.
DECODE_DIR=${VMM_DECODE_DIR:-$HOME/p2_match}

while [ $# -gt 0 ]; do
    case "$1" in
        --jobs) JOBS=$2; shift 2 ;;
        --dry-run) DRY=1; shift ;;
        --force) FORCE=1; shift ;;
        --only) ONLY=$2; shift 2 ;;
        --decode-dir) DECODE_DIR=$2; shift 2 ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
done

CODE=$(cd "$(dirname "$0")" && pwd)
WORK=${TMPDIR:-/tmp}/p2_match_sweep.$USER
mkdir -p "$WORK"

# --- the work list: DREAM sub_runs that also have VMM data ------------------
LIST=$WORK/joblist.txt
: > "$LIST"
for d in "$BASE"/runs/run_*/*/combined_hits_root; do
    sub=$(basename "$(dirname "$d")")
    run=$(basename "$(dirname "$(dirname "$d")")")
    [ -n "$ONLY" ] && [[ ",$ONLY," != *",$run,"* ]] && continue
    v=$BASE/vmm/runs/$run/$sub
    # VMM side: reduced columns if the online pass kept them, raw captures
    # otherwise -- extract_vmm_triggers picks, and the two agree bit for bit.
    if [ -d "$v/hits_store" ] || [ -d "$v/raw_daq_data" ]; then
        if [ "$FORCE" = 0 ] && [ -f "$OUT/match_${run}_${sub}.json" ]; then
            continue
        fi
        echo "$run $sub" >> "$LIST"
    fi
done
n=$(wc -l < "$LIST")
echo "$n sub_run(s) to match (jobs=$JOBS, force=$FORCE)"
[ "$DRY" = 1 ] && { cat "$LIST"; exit 0; }
[ "$n" = 0 ] && exit 0

# --- one sub_run ------------------------------------------------------------
cat > "$WORK/one.sh" <<EOS
#!/bin/bash
set -uo pipefail
run=\$1; sub=\$2
tag=\${run}_\${sub}
d=$WORK/\$tag
mkdir -p "\$d"
log=\$d/log.txt
trg=\$d/vmm_triggers_\$tag.npz
json=\$d/match_\$tag.json
frac_of() { python3 -c "import json;print(json.load(open('\$1'))['match_frac_dream'])" 2>/dev/null || echo 0; }
# a subshell, so a failing stage leaves only the subshell and this script can
# still report which sub_run died and why
(
  echo "=== \$tag  \$(date -u +%H:%M:%S)"
  # the usual path: the trigger fan-out on VMM 0 ch 44
  if python3 $CODE/extract_vmm_triggers.py "\$run" "\$sub" --out "\$trg"; then
      python3 $CODE/match_streams.py "\$run" "\$sub" --out "\$d" --tol-us 2 \
              --vmm-npz "\$trg" || exit 2
  fi
  # ... but the cabling is not the same all campaign long (run_25 has no hits
  # on that channel at all), so when the default channel gives no usable lock,
  # go and find which channel the trigger is really on and match again.
  f=\$(frac_of "\$json")
  if python3 -c "import sys;sys.exit(0 if \$f < 0.5 else 1)"; then
      echo "--- default channel gave frac=\$f, searching for the trigger channel"
      python3 $CODE/find_trigger_channel.py "\$run" "\$sub" --out "\$trg" || {
          # no channel coincides with DREAM at all: keep the failed match as
          # the record of that, rather than leaving the sub_run to be retried
          # by every future sweep
          echo "--- no trigger channel found"
          [ -f "\$json" ] && exit 0 || exit 4; }
      python3 $CODE/match_streams.py "\$run" "\$sub" --out "\$d" --tol-us 2 \
              --vmm-npz "\$trg" || exit 2
  fi
  [ -f "\$json" ] || exit 1
) > "\$log" 2>&1
rc=\$?
if [ \$rc -ne 0 ]; then
    echo "FAIL(\$rc) \$tag  \$(tail -2 "\$log" | tr '\n' ' ')"
    rm -rf "\$d"
    exit \$rc
fi
for f in "\$d"/vmm_triggers_\$tag.npz "\$d"/match_\$tag.json "\$d"/match_\$tag.npz; do
    xrdcp -f -s "\$f" "$EOS_URL/$OUT/\$(basename \$f)" || {
        echo "UPLOAD-FAIL \$tag \$(basename \$f)"; exit 3; }
done
frac=\$(python3 -c "import json;print('%.3f'%json.load(open('\$d/match_\$tag.json'))['match_frac_dream'])")
rms=\$(python3 -c "import json;d=json.load(open('\$d/match_\$tag.json'));print('%.1f'%(d['residual_rms_ns'] or -1))")
echo "OK \$tag  frac=\$frac  rms=\${rms}ns"
rm -rf "\$d"
EOS
chmod +x "$WORK/one.sh"

# shellcheck disable=SC1090
# the LCG setup.sh reads variables it never sets, so it dies instantly under
# `set -u`; drop the flag across the source and put it back
set +u
# shellcheck disable=SC1090
source "$LCG"
set -u
[ -f "$DECODE_DIR/vmm_decode.py" ] || echo \
    "WARNING: no vmm_decode.py in $DECODE_DIR — sub_runs without reduced hit" \
    "columns will fail (pass --decode-dir)"
export PYTHONPATH=$CODE:$DECODE_DIR:${PYTHONPATH:-}
xrdfs "$EOS_URL" mkdir -p "$OUT" 2>/dev/null

xargs -a "$LIST" -n 2 -P "$JOBS" "$WORK/one.sh"
echo "SWEEP_DONE $(date -u +%H:%M:%S)"
