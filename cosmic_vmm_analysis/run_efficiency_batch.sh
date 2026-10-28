#!/bin/bash
# Per-run efficiency and straight-track rate, one run at a time (never pooled), same map and veto for all runs given.
# Runs on lxplus from ~/cva (LCG view + VMM_QA_DIR set by the caller).
#
# Usage: run_efficiency_batch.sh OUT_DIR CONN_SHIFT "VETO" "ORDER1;ORDER2;..." RUN [RUN ...]
#   e.g. run_efficiency_batch.sh ~/tpold 0 2:11,11:35 "P2_IN:4=bot:rev|top:lin;P2_IN:6=bot:rev|top:lin" $(seq 76 111)
# Candidates (dump_candidates.py output) are read from ~/tp/cand/runN.npz.
set -u
out=$1; shift
shift_c=$1; shift
veto=$1; shift
orders=$1; shift
runs_dir=/eos/experiment/ntof/data/x17/p2_sps_july/vmm_cosmic_run_oct26/runs
mkdir -p "$out"
oargs=()
IFS=';' read -ra ol <<< "$orders"
for o in "${ol[@]}"; do [ -n "$o" ] && oargs+=(--order "$o"); done
for r in "$@"; do
    d=$(ls -d $runs_dir/run_$r/cfg_* | head -1)
    files=$(ls $d/raw_daq_data/*.pcapng | sort | paste -sd,)
    python3 efficiency_tp.py "$files" "$out/eff_run$r.json" --conn-shift="$shift_c" --veto "$veto" "${oargs[@]}" \
        --title "run_$r $(basename $d | sed 's/cfg_gain4.5_peaktime200_opt_cosmic_//')" > "$out/eff_run$r.log" 2>&1
    python3 straight_track_rate.py "$out/str_run$r.json" "run_$r:$HOME/tp/cand/run$r.npz" --conn-shift="$shift_c" \
        --veto "$veto" "${oargs[@]}" > "$out/str_run$r.log" 2>&1
    echo "run_$r: $(grep -c probe_det "$out/eff_run$r.json" 2>/dev/null) probes; $(tail -1 "$out/str_run$r.log")"
done
