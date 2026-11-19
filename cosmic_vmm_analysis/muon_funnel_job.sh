#!/bin/bash
# HTCondor job (lxplus): muon_funnel.py on every weekend run except the pulser runs 270 / 283, one process per run,
# 4 in parallel (~1.5 GB each). Submit from an AFS dir with muon_funnel.sub; outputs $W/analysis/muon_funnel/run*.json
W=/eos/experiment/ntof/data/x17/p2_sps_july/vmm_cosmic_run_oct26/analysis_lxplus/weekend_20261119
cd $W || exit 1
source pipeline/env.sh
export PYTHONPATH=$W/pipeline/cva:$PYTHONPATH
RUNS=$(python3 -c "
import json; d=json.load(open('pseudo/chunks.json'))
print(' '.join(str(r) for r in sorted({v['run'] for v in d.values()}) if r not in (270, 283)))")
echo "start $(date '+%F %T') on $(hostname): runs $RUNS"
one() {
    r=$1
    J=$(python3 -c "
import json; d=json.load(open('pseudo/chunks.json'))
print(' '.join(f'tags/eff_run{k}.json' for k, v in sorted(d.items()) if v['run'] == $r))")
    python3 analysis/muon_funnel/muon_funnel.py $J --out analysis/muon_funnel/run$r.json > analysis/muon_funnel/run$r.log 2>&1
    echo "run $r rc=$? $(date '+%T')"
}
export -f one
echo $RUNS | tr ' ' '\n' | xargs -P 4 -I{} bash -c 'one {}'
echo "end $(date '+%F %T')"
