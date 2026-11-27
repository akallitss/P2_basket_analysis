#!/usr/bin/env bash
# run_lxplus.sh -- build water_scan on lxplus (Garfield++ from an LCG view) and run
# one detached Magboltz process per gas key of gaslist.txt (or of the keys given).
#
#   ./run_lxplus.sh launch [key ...]   # stage, build, launch detached (max NPAR at once)
#   SKIP_BUILD=1 GASLIST=gaslist_combos.txt NPAR=4 ./run_lxplus.sh launch   # extra batch, same binary
#   ./run_lxplus.sh status             # row counts per gas
#   ./run_lxplus.sh fetch              # copy CSVs + logs into results/magboltz/
#
# Env: SSH_OPTS (default: reuse the open ControlMaster socket), LCG_VIEW, NCOLL (x1e7),
#      NPAR (parallel processes, default 6), GASLIST (alternative list file).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HOST="${HOST:-akallits@lxplus.cern.ch}"
SSH_OPTS="${SSH_OPTS:--o BatchMode=yes -o ControlPath=$HOME/.ssh/master-%r@%h:%p}"
LCG_VIEW="${LCG_VIEW:-/cvmfs/sft.cern.ch/lcg/views/LCG_109/x86_64-el9-gcc13-opt}"
NCOLL="${NCOLL:-1}"
NPAR="${NPAR:-6}"
GASLIST="${GASLIST:-$HERE/gaslist.txt}"
RB=water_scan            # remote dir under AFS $HOME
SSH="ssh $SSH_OPTS $HOST"
SCP="scp $SSH_OPTS"
cmd="${1:-status}"; shift || true

case "$cmd" in
launch)
  sel=$(mktemp)
  grep -v '^#' "$GASLIST" | grep . > "$sel.all"
  if [ $# -gt 0 ]; then for k in "$@"; do grep -E "^$k[[:space:]]" "$sel.all"; done > "$sel"; else cp "$sel.all" "$sel"; fi
  $SSH "mkdir -p ~/$RB/src ~/$RB/out"
  [ "${SKIP_BUILD:-0}" = 1 ] || $SCP "$HERE/water_scan.cpp" "$HERE/CMakeLists.txt" "$HOST:~/$RB/src/"
  TODO="todo_$(date +%s).txt"
  $SCP "$sel" "$HOST:~/$RB/src/$TODO"
  $SSH "LCG_VIEW=$LCG_VIEW NCOLL=$NCOLL NPAR=$NPAR RB=$RB TODO=$TODO SKIP_BUILD=${SKIP_BUILD:-0} bash -s" <<'REMOTE'
set -eo pipefail
set +u; source "$LCG_VIEW/setup.sh"; set -u
cd ~/$RB/src
if [ "${SKIP_BUILD:-0}" != 1 ]; then cmake -DCMAKE_BUILD_TYPE=Release -S . -B build >/dev/null && cmake --build build -j4 >/dev/null; fi
echo "built on $(hostname) with $LCG_VIEW"
OUT=~/$RB/out
setsid nohup bash -c "
  cat ~/$RB/src/$TODO | xargs -P $NPAR -L 1 bash -c '
    k=\$0; shift 0
    ~/$RB/src/build/water_scan $OUT/\$k.csv $NCOLL 293.15 760 \"\$@\" > $OUT/\$k.log 2>&1 && touch $OUT/\$k.DONE || touch $OUT/\$k.FAIL'
" > $OUT/.launch_$TODO.log 2>&1 < /dev/null &
disown || true
echo launched
REMOTE
  ;;
status)
  $SSH "cd ~/$RB/out && for f in *.csv; do echo \"\$f \$((\$(wc -l < \$f)-1)) \$(ls \${f%.csv}.DONE \${f%.csv}.FAIL 2>/dev/null)\"; done; uptime" ;;
fetch)
  mkdir -p "$HERE/results/magboltz"
  $SCP "$HOST:~/$RB/out/*.csv" "$HOST:~/$RB/out/*.log" "$HERE/results/magboltz/" ;;
*) echo "unknown command $cmd"; exit 1 ;;
esac
