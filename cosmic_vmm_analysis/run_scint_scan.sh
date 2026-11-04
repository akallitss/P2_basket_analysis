#!/bin/bash
# Scintillator-referenced efficiency for every finished run of a threshold scan, one run at a time (never pooled).
# Runs on lxplus from ~/cva2 (LCG view + VMM_QA_DIR set by the caller). A run is done when its capture is closed
# (.capture_done) and it has at least MIN_FILES capture files; its outputs are skipped if already there, so the
# script can be re-run while the scan goes on.
#
# Sparks: when the run has an hv_monitor.csv (runs >= 170) spark_veto.py writes OUT_DIR/veto_run_N.json and both
# scintillator efficiencies drop the triggers in its intervals (--bad-times); NOVETO=1 turns this off.
#
# Usage: run_scint_scan.sh OUT_DIR LABEL_PATTERN MIN_FILES
#   e.g. run_scint_scan.sh ~/wk1 adwk1 100     (the 3-hour points of the weekend scan, map of runs >= 112)
#        run_scint_scan.sh ~/wk2 adwk2 100     (the November scan, HV-monitor spark veto applied)
set -u
out=$1; pat=$2; minf=$3
runs_dir=/eos/experiment/ntof/data/x17/p2_sps_july/vmm_cosmic_run_oct26/runs
order=(--conn-shift=-1 --order "P2_IN:3=bot:rev|top:lin" --order "P2_IN:5=bot:rev|top:lin")
mkdir -p "$out"
for d in $(ls -d $runs_dir/run_*/cfg_*"$pat"* 2>/dev/null | sort -t_ -k2 -n); do
    r=$(basename "$(dirname "$d")" | sed 's/run_//')
    n=$(ls "$d"/raw_daq_data/*.pcapng 2>/dev/null | wc -l)
    [ "$n" -lt "$minf" ] && continue
    [ -f "$d/raw_daq_data/.capture_done" ] || continue
    [ -f "$out/eff_run$r.json" ] && continue
    files=$(ls "$d"/raw_daq_data/*.pcapng | sort | paste -sd,)
    lab=$(basename "$d" | sed 's/cfg_gain4.5_peaktime200_opt_cosmic_//')
    bt=()
    if [ -z "${NOVETO:-}" ] && [ -f "$d/hv_monitor.csv" ]; then
        python3 spark_veto.py "$(dirname "$d")" --out "$out" > "$out/veto_run$r.log" 2>&1 && bt=(--bad-times "$out/veto_run_$r.json")
    fi
    python3 scint_efficiency.py "$files" "$out/scint3_run$r.json" --trig 0:40 "${order[@]}" "${bt[@]}" --title "run_$r $lab" \
        > "$out/scint3_run$r.log" 2>&1
    python3 scint_efficiency.py "$files" "$out/scintor_run$r.json" --top 0:41,0:45 --bottom 0:44 "${order[@]}" "${bt[@]}" \
        --title "run_$r $lab" > "$out/scintor_run$r.log" 2>&1
    python3 efficiency_tp.py "$files" "$out/eff_run$r.json" "${order[@]}" --dead-vmm 0,1 --title "run_$r $lab" \
        > "$out/eff_run$r.log" 2>&1
    echo "run_$r $lab: $n files done"
done
