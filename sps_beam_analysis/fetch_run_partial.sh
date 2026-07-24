#!/usr/bin/env bash
# fetch_run_partial.sh -- pull only the files the analysis needs for one run,
# from banco to the local LaCie mirror, skipping the huge raw .fdf files and the
# decoded_root / hits_root intermediates. Safe to run against a LIVE run: it only
# fetches sub_runs that carry a .subrun_complete marker, so the sub_run the DAQ
# is currently writing (whose combined-hits file may be partial) is left alone.
# Re-run it as the run fills in; rsync is incremental.
#
#   ./fetch_run_partial.sh <run_name>
#
# What it fetches per completed sub_run:
#   combined_hits_root/*.root   (what every stage reads)
#   hv_monitor.csv              (HV / spark veto)
#   raw_daq_data/run_time.txt   (durations)
# plus, once per run: run_config.json and config/ + dream_config/ if present.
set -uo pipefail

RUN=${1:?usage: fetch_run_partial.sh <run_name>}
REMOTE=${REMOTE:-banco_cern}
RBASE=${RBASE:-/local/home/banco/P2_data/TB_July2026_H4}
LBASE=${LBASE:-/media/ak271430/LaCie/Extras/Physics/Post-Doc-Saclay/data/SPS_beam_test/TB_July2026_H4}
RSYNC="rsync -rt --no-perms --no-owner --no-group --modify-window=2 --partial --timeout=120"

echo "== $RUN : $REMOTE:$RBASE/runs/$RUN  ->  $LBASE/runs/$RUN"

# 1. Which sub_runs are complete? (marker file written by the DAQ at sub_run end)
mapfile -t DONE < <(ssh -o ConnectTimeout=30 -o ServerAliveInterval=15 "$REMOTE" \
  "for d in $RBASE/runs/$RUN/*/; do [ -f \$d/.subrun_complete ] && basename \$d; done" 2>/dev/null)
if [ ${#DONE[@]} -eq 0 ]; then
  echo "  no completed sub_runs yet (or banco unreachable) -- nothing to fetch"
  exit 0
fi
echo "  completed sub_runs: ${DONE[*]}"

# 2. Run-level metadata (small, always refresh).
mkdir -p "$LBASE/runs/$RUN"
$RSYNC "$REMOTE:$RBASE/runs/$RUN/run_config.json" "$LBASE/runs/$RUN/" 2>/dev/null
$RSYNC "$REMOTE:$RBASE/config"        "$LBASE/" 2>/dev/null
$RSYNC "$REMOTE:$RBASE/dream_config"  "$LBASE/" 2>/dev/null

# 3. Per completed sub_run, only the analysis inputs.
for s in "${DONE[@]}"; do
  d="$LBASE/runs/$RUN/$s"
  mkdir -p "$d/combined_hits_root" "$d/raw_daq_data"
  $RSYNC "$REMOTE:$RBASE/runs/$RUN/$s/combined_hits_root/" "$d/combined_hits_root/" 2>/dev/null
  $RSYNC "$REMOTE:$RBASE/runs/$RUN/$s/hv_monitor.csv"      "$d/" 2>/dev/null
  $RSYNC "$REMOTE:$RBASE/runs/$RUN/$s/raw_daq_data/run_time.txt" "$d/raw_daq_data/" 2>/dev/null
  # per-FEU recorded-trigger sets (extract_recorded_events.py on the DAQ host)
  $RSYNC "$REMOTE:$RBASE/runs/$RUN/$s/recorded_events.npz" "$d/" 2>/dev/null
  n=$(ls "$d/combined_hits_root/"*.root 2>/dev/null | wc -l)
  echo "  $s : $n combined-hits file(s)"
done
echo "== done: $(du -sh "$LBASE/runs/$RUN" 2>/dev/null | cut -f1) in $LBASE/runs/$RUN"
