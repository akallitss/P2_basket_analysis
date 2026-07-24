#!/usr/bin/env bash
# run_beam_qa.sh -- the beam-time "is this run good?" pass, on banco.
#
#   ./run_beam_qa.sh <run_name> [extra args for 25_commissioning_qa.py]
#
# Runs, for every sub_run of <run_name> that already has combined-hits ROOT:
#   25_commissioning_qa.py   rates / occupancy / latency / HV / verdicts
#   24_event_sync_qa.py      FEU-to-FEU trigger alignment
# and writes under $SPS_ANALYSIS_ROOT, which is where the DAQ GUI's Analysis
# tab reads from and what the EOS backup watcher syncs.
#
# Safe to re-run while the DAQ is still writing: sub_runs the processor has not
# combined yet are simply skipped, so re-run it as the run fills in.
set -euo pipefail

RUN=${1:?usage: run_beam_qa.sh <run_name> [extra args]}
shift || true

BASE=${P2_BASE:-/local/home/banco/P2_data/TB_July2026_H4}
export SPS_DATA_ROOT=${SPS_DATA_ROOT:-$BASE/runs}
export SPS_ANALYSIS_ROOT=${SPS_ANALYSIS_ROOT:-$BASE/analysis}
export SPS_RUN=$RUN
PY=${SPS_PYTHON:-/local/home/banco/DAQ_Control_Dream_Beam/.venv/bin/python}

unset PYTHONPATH || true          # banco's ISEG SDK PYTHONPATH shadows ROOT/uproot
cd "$(dirname "$0")"

echo "=== $RUN   data=$SPS_DATA_ROOT   analysis=$SPS_ANALYSIS_ROOT"
"$PY" 25_commissioning_qa.py live "$@"
echo
echo "=== event sync"
"$PY" 24_event_sync_qa.py live
echo
echo "=== HV spark QA"
"$PY" 26_hv_spark_qa.py live
echo
echo "=== pedestal QA"
"$PY" 27_pedestal_qa.py live --with-rate
echo
echo "=== timing QA"
"$PY" 28_timing_qa.py live
