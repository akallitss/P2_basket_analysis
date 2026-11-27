#!/bin/bash
# HTCondor job (lxplus): muon_funnel.py on every weekend run except the pulser runs 270 / 283, one process per run,
# 4 in parallel (~1.5 GB each). Submit from an AFS dir with muon_funnel.sub; outputs $W/analysis/muon_funnel/run*.json
# OUTDIR (default analysis/muon_funnel) and EXTRA (extra muon_funnel.py options, e.g. --pointing FILE) from the
# environment (muon_funnel_pointing.sub)
W=${W:-/eos/experiment/ntof/data/x17/p2_sps_july/vmm_cosmic_run_oct26/analysis_lxplus/weekend_20261119}
# S = the area holding analysis/muon_funnel/*.py, paddles and pointing.json (default: W itself)
S=${S:-$W}
cd $W || exit 1
export OUTDIR=${OUTDIR:-analysis/muon_funnel} EXTRA=${EXTRA:-}
source pipeline/env.sh
export PYTHONPATH=$W/pipeline/cva:$PYTHONPATH
RUNS=$(python3 -c "
import json; d=json.load(open('pseudo/chunks.json'))
print(' '.join(str(r) for r in sorted({v['run'] for v in d.values() if '_tp' not in v['cfg']})))")
echo "start $(date '+%F %T') on $(hostname): runs $RUNS"
one() {
    r=$1
    # already done (delete the json, or FORCE=1, to redo); runs not dumped yet are skipped by having no chunks
    if [ -s $OUTDIR/run$r.json ] && [ -z "$FORCE" ]; then echo "run $r: exists"; return; fi
    J=$(python3 -c "
import json, os; d=json.load(open('pseudo/chunks.json'))
f = [f'tags/eff_run{k}.json' for k, v in sorted(d.items()) if v['run'] == $r]
print(' '.join(f) if f and all(os.path.exists(x) for x in f) else '')")
    if [ -z "$J" ]; then echo "run $r: not fully dumped yet, skipped"; return; fi
    python3 $S/analysis/muon_funnel/muon_funnel.py $J --paddles $S/analysis/muon_funnel/paddles_fit_20261104.json $EXTRA \
        --out $OUTDIR/run$r.json > $OUTDIR/run$r.log 2>&1
    echo "run $r rc=$? $(date '+%T')"
}
export -f one
echo $RUNS | tr ' ' '\n' | xargs -P 4 -I{} bash -c 'one {}'
echo "end $(date '+%F %T')"
