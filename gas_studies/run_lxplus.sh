#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# run_lxplus.sh -- build p2_gas_scan against Garfield++ on lxplus and run the
# Magboltz transport scan for one or more gases, fetching the CSVs into results/.
#
#   ./run_lxplus.sh                      # scan every gas in gases.py
#   ./run_lxplus.sh ar_co2_iso_95_3_2    # scan one gas by key
#
# Env overrides:
#   LXPLUS_HOST      ssh target                (default: lxplus)
#   LCG_VIEW         LCG view with Garfield++  (default: LCG_105 el9/gcc13)
#   REMOTE_BASE      AFS staging dir on lxplus (default: p2_gas_scan, rel to $HOME)
#   EXTRA_SSH_OPTS   extra ssh/scp options, e.g. reuse a live control socket:
#                    EXTRA_SSH_OPTS="-o ControlPath=$HOME/.ssh/master-...@lxplus..:22"
#
# Resilience: the Magboltz scan can take tens of minutes, so it is launched
# DETACHED (nohup) on the worker node and writes an AFS ".DONE"/".FAIL" marker
# when finished -- a dropped ssh connection or expired Kerberos ticket no longer
# kills it. We then poll the (shared) AFS marker and fetch the CSVs. If polling
# is interrupted you can re-fetch later; nothing is recomputed.
#
# Notes: lxplus is load-balanced, so /tmp is node-local and not shared between
# ssh calls -- sources and CSVs live under AFS ($REMOTE_BASE, shared across
# nodes) but the build + run happen in a node-local /tmp dir (AFS quota is tight).
# -----------------------------------------------------------------------------
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LXPLUS_HOST="${LXPLUS_HOST:-lxplus}"
LCG_VIEW="${LCG_VIEW:-/cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc13-opt}"
REMOTE_BASE="${REMOTE_BASE:-p2_gas_scan}"   # relative to remote $HOME (AFS)
EXTRA_SSH_OPTS="${EXTRA_SSH_OPTS:-}"

# keepalive so a brief network hiccup does not tear the connection down
SSH="ssh -o BatchMode=yes -o ServerAliveInterval=30 -o ServerAliveCountMax=6 $EXTRA_SSH_OPTS"
SCP="scp -o BatchMode=yes -o ServerAliveInterval=30 -o ServerAliveCountMax=6 $EXTRA_SSH_OPTS"

mkdir -p "$HERE/results"

# --fetch : pull whatever CSVs already exist remotely and exit (no re-launch).
# The scanner flushes every row as it is computed, so a partial fetch mid-run
# gives you the drift/timing rows while the slow Townsend scan is still going.
FETCH_ONLY=0
if [ "${1:-}" = "--fetch" ]; then FETCH_ONLY=1; shift; fi

# 1. Build the gas list (single source of truth = gases.py) ------------------
GASLIST="$HERE/results/.gaslist.txt"
python3 "$HERE/gases.py" --emit-args "$@" > "$GASLIST"

if [ "$FETCH_ONLY" = "1" ]; then
  echo "== fetching current remote CSVs (partial runs included) =="
  while read -r key _; do
    [ -z "$key" ] && continue
    if $SCP "$LXPLUS_HOST:~/$REMOTE_BASE/out/$key.csv" "$HERE/results/$key.csv" 2>/dev/null; then
      echo "  results/$key.csv  ($(( $(wc -l < "$HERE/results/$key.csv") - 1 )) rows)"
    else
      echo "  (no $key.csv yet)"
    fi
  done < "$GASLIST"
  exit 0
fi
echo "== gases to scan =="
cat "$GASLIST"

# 2. Stage sources on AFS (shared across all lxplus nodes) -------------------
echo "== staging sources on $LXPLUS_HOST:~/$REMOTE_BASE =="
$SSH "$LXPLUS_HOST" "mkdir -p ~/$REMOTE_BASE/src ~/$REMOTE_BASE/out"
$SCP "$HERE/p2_gas_scan.cpp" "$HERE/CMakeLists.txt" "$GASLIST" \
     "$LXPLUS_HOST:~/$REMOTE_BASE/src/"
$SSH "$LXPLUS_HOST" "mv ~/$REMOTE_BASE/src/.gaslist.txt ~/$REMOTE_BASE/src/gaslist.txt"

# 3. Build, then launch the scan DETACHED (survives connection drops) --------
echo "== building on $LXPLUS_HOST (LCG: $LCG_VIEW), then launching detached scan =="
$SSH "$LXPLUS_HOST" \
    "LCG_VIEW='$LCG_VIEW' REMOTE_BASE='$REMOTE_BASE' bash -s" <<'REMOTE'
set -eo pipefail
OUT="$HOME/$REMOTE_BASE/out"
rm -f "$OUT/.DONE" "$OUT/.FAIL" "$OUT"/*.DONE "$OUT"/*.FAIL
echo "-- sourcing LCG view"
set +u                       # LCG setup.sh references unbound vars
source "$LCG_VIEW/setup.sh"
set -u
WORK="$(mktemp -d /tmp/p2gas.XXXXXX)"
cp "$HOME/$REMOTE_BASE/src/p2_gas_scan.cpp" "$HOME/$REMOTE_BASE/src/CMakeLists.txt" "$WORK/"
cd "$WORK"
echo "-- cmake configure"; cmake -DCMAKE_BUILD_TYPE=Release -S . -B build >/dev/null
echo "-- make";           cmake --build build -j4 >/dev/null
echo "-- launching detached scan on $(hostname)"
# One process PER GAS, all in parallel: Magboltz is single-threaded and the
# gases are independent, so N gases cost the wall time of the slowest one
# instead of the sum. Each gas drops its own <key>.DONE/<key>.FAIL marker and
# its own log, so a single failing mixture no longer hides the others.
setsid nohup bash -c '
  while read -r key rest; do
    [ -z "$key" ] && continue
    echo ">>> launching $key : $rest"
    (
      "'"$WORK"'/build/p2_gas_scan" "'"$OUT"'/$key.csv" $rest \
        && touch "'"$OUT"'/$key.DONE" || touch "'"$OUT"'/$key.FAIL"
    ) > "'"$OUT"'/$key.log" 2>&1 &
  done < "'"$HOME"'/'"$REMOTE_BASE"'/src/gaslist.txt"
  wait
  touch "'"$OUT"'/.DONE"
  rm -rf "'"$WORK"'"
' > "$OUT/.scan.log" 2>&1 < /dev/null &
disown || true
echo "-- detached scan launched (pid $!). Poll marker: $OUT/.DONE"
REMOTE

# 4. Poll the AFS marker until the scan finishes -----------------------------
echo "== waiting for remote scan (polling AFS marker; safe to interrupt) =="
done=0
for i in $(seq 1 300); do            # up to ~100 min
  # one round trip: finished marker + per-gas row counts (progress readout)
  st=$($SSH "$LXPLUS_HOST" "cd ~/$REMOTE_BASE/out 2>/dev/null || exit 0
       test -f .DONE && echo ALLDONE
       for f in *.csv; do [ -e \"\$f\" ] && echo \"  \${f%.csv}: \$(( \$(wc -l < \"\$f\") - 1 )) rows\"; done
       for f in *.FAIL; do [ -e \"\$f\" ] && echo \"  !! \${f%.FAIL} FAILED\"; done" 2>/dev/null || true)
  echo "$st" | grep -v ALLDONE | grep . || true
  if echo "$st" | grep -q ALLDONE; then echo "-- remote scan DONE"; done=1; break; fi
  sleep 20
done
if [ "$done" != "1" ]; then
  echo "!! poll timed out; the scan may still be running remotely."
  echo "   Re-run this script (it will re-fetch once ~/$REMOTE_BASE/out/.DONE exists),"
  echo "   or fetch manually from ~/$REMOTE_BASE/out/."
  exit 3
fi

# 5. Fetch CSVs back ---------------------------------------------------------
echo "== fetching results =="
while read -r key _; do
  [ -z "$key" ] && continue
  if $SCP "$LXPLUS_HOST:~/$REMOTE_BASE/out/$key.csv" "$HERE/results/$key.csv" 2>/dev/null; then
    echo "  results/$key.csv"
  else
    echo "  !! $key.csv missing -- see ~/$REMOTE_BASE/out/$key.log"
  fi
done < "$GASLIST"

echo "== done. Now run:  python3 analyze.py"
