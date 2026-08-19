#!/usr/bin/env bash
# Timing budget at the nominal working point + the dt peak fits per detector
# per gas.  RUN THIS ON LXPLUS -- there is no xrdcp and no EOS mount on the
# laptop, and the CERN principal is `akallits`, not the local username.
#
#   ssh lxplus            # the ssh config alias; keeps a 1-day control master
#   ./run_timing_nominal.sh
#
# Two things this has to work around, both learned the hard way:
#
#  * lxplus's default python3 has no pandas, which vmm_decode needs, so an LCG
#    view is sourced below.
#  * MOST sub_runs were reduced with --drop-columns and keep only
#    counts/meta/scalars, and the ENTIRE gas-B period (run_61 onward) has no
#    hits_store at all.  Only these keep full columns: run_32/33/34/35/36,
#    38, 40-43, 49-52, 56.  Everything else has to be decoded from pcapng,
#    which is what decode_to_store.py is for (~12 s per capture).
#
# Working points, chosen from vmm_timing_by_subrun.csv by smallest sigma among
# sub_runs that still have a real signal (efficiency > 0.30).  Without that cut
# the "best" points are dead detectors: mesh 350 V gives sigma 4.8 ns at zero
# efficiency, below the 6.5 ns quantisation floor, so not a coincidence at all.
set -euo pipefail

source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc13-opt/setup.sh

E=/eos/experiment/ntof/data/x17/p2_sps_july/vmm/runs
HERE=$(cd "$(dirname "$0")" && pwd)
WORK=${TMPDIR:-/tmp}/$(whoami)/vmmtiming
mkdir -p "$WORK"
cd "$WORK"

# --- gas B has no column store: rebuild it from the raw captures ----------- #
if [ ! -f "$WORK/store_gasB/.done" ]; then
  python3 "$HERE/decode_to_store.py" \
      "$E/run_66/cfg_gain4.5_peaktime200_opt/raw_daq_data" \
      "$WORK/store_gasB" --captures 12
  touch "$WORK/store_gasB/.done"
fi

# --- the dt peaks ---------------------------------------------------------- #
python3 "$HERE/vmm_timing_peaks.py" --store "$E/run_36/operating_00/hits_store" \
    --tag nomA --captures 24 \
    --label "gas A - run_36/operating_00, mesh 450 / drift 750, gain 3.0, peak 200 ns"

python3 "$HERE/vmm_timing_peaks.py" --store "$E/run_35/driftscan_gap400V/hits_store" \
    --tag outA --captures 16 \
    --label "gas A - run_35/driftscan_gap400V, mesh 450 / drift 850, gain 3.0, peak 200 ns"

python3 "$HERE/vmm_timing_peaks.py" --store "$WORK/store_gasB" \
    --tag gasB --captures 12 \
    --label "gas B - run_66/cfg_gain4.5_peaktime200_opt, mesh 450 / drift 750, gain 4.5, peak 200 ns"

# --- the budget, at the nominal point, both estimators --------------------- #
python3 "$HERE/vmm_timing_budget.py" --store "$E/run_36/operating_00/hits_store" \
    --captures 10 | tee "$WORK/TIMING_BUDGET_nominal.txt"

echo
echo "npz + budget in $WORK ; copy them back and run"
echo "  python3 mpgd2026/make_timing_peaks_fig.py timing_peaks_nomA.npz timing_peaks_gasB.npz \\"
echo "      --col-labels 'gas A  Ar/CO2/iC4H10 93/5/2' 'gas B  Ar/CF4/iC4H10 88/10/2' \\"
echo "      --stations P2_MID P2_OUT -o mpgd2026/figs --name vmm_timing_peaks"
