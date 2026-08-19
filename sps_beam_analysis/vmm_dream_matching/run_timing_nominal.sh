#!/usr/bin/env bash
# Redo the timing budget at the nominal working point, and fit the dt peak at
# each detector's best-timing point per gas.
#
# Needs a Kerberos ticket:   kinit akallits@CERN.CH
# Run it on lxplus, or locally once the sub_runs below are staged.
#
# The four working points come from vmm_timing_by_subrun.csv, taking the
# smallest sigma among sub_runs that still have a real signal (efficiency
# > 0.30).  Without that cut the "best" points are dead detectors: mesh 350 V
# gives sigma 4.8 ns at zero efficiency, which is below the 6.5 ns
# quantisation floor and is not a coincidence at all.
set -euo pipefail

EOS=/eos/experiment/ntof/data/x17/p2_sps_july/vmm/runs
WORK=${TMPDIR:-/tmp}/vmm_timing
OUT=$(cd "$(dirname "$0")" && pwd)
FIGS=$OUT/../../mpgd2026/figs
mkdir -p "$WORK"

# run/sub_run                              tag           label
POINTS=(
  "run_48/cfg_gain4.5_peaktime200_opt|midA|P2_MID best, gas A - mesh 450 / drift 750, gain 4.5, peak 200 ns"
  "run_66/cfg_gain4.5_peaktime200_opt|midB|P2_MID + P2_OUT best, gas B - mesh 450 / drift 750, gain 4.5, peak 200 ns"
  "run_57/driftscan_gap400V|outA|P2_OUT best, gas A - mesh 450 / drift 850, gain 3.0, peak 200 ns"
  "run_36/operating_00|nomA|nominal working point, gas A - mesh 450 / drift 750, gain 3.0, peak 200 ns"
)

for p in "${POINTS[@]}"; do
  IFS='|' read -r sub tag label <<<"$p"
  dst=$WORK/$tag
  if [ ! -d "$dst" ]; then
    echo "staging $sub -> $dst"
    mkdir -p "$dst"
    # hits_store holds the reduced column store; fall back to the raw captures
    xrdcp -r --silent "root://eosuser.cern.ch/$EOS/$sub/hits_store" "$dst" \
      || cp -r "$EOS/$sub/hits_store/." "$dst/"
  fi
  echo "== $tag: $label"
  python3 "$OUT/vmm_timing_peaks.py" --store "$dst" --tag "$tag" \
      --captures 24 --label "$label"
done

echo
echo "== timing budget at the nominal working point =="
python3 "$OUT/vmm_timing_budget.py" --store "$WORK/nomA" --captures 12 \
    | tee "$OUT/TIMING_BUDGET_nominal.txt"

echo
python3 "$OUT/../../mpgd2026/make_timing_peaks_fig.py" \
    "$OUT/timing_peaks_midA.npz" "$OUT/timing_peaks_midB.npz" \
    --col-labels "gas A  Ar/CO2/iC4H10 93/5/2" "gas B  Ar/CF4/iC4H10 88/10/2" \
    --stations P2_MID P2_OUT -o "$FIGS" --name vmm_timing_peaks
echo "figures in $FIGS"
