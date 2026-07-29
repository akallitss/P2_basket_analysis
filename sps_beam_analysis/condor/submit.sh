#!/usr/bin/env bash
# submit.sh -- drive an HTCondor analysis sweep on lxplus from banco.
#
#   ./submit.sh [--group rec|hits|wave] [--source ntof|salsachip|user|all]
#              [--run RUN[,RUN...]] [--dry-run]
#
# Packs the analysis code, builds the joblist from what is actually on EOS,
# ships both to lxplus and submits. The code is packed FRESH every time and
# travels in the input sandbox, so a submission is pinned to the code as it
# stands right now -- editing the repo afterwards cannot retroactively change
# what a running sweep is doing.
#
# Kerberos: every ssh/scp here needs GSSAPIDelegateCredentials. Without it the
# session reaches lxplus but has no AFS token, and cannot so much as read its
# own home directory -- an error that looks nothing like a Kerberos problem.
set -euo pipefail

GROUP=hits
RUNS=""
SOURCE=ntof
DRY=0
while [ $# -gt 0 ]; do
    case "$1" in
        --group)   GROUP=$2; shift 2 ;;
        --run)     RUNS=$2; shift 2 ;;
        --source)  SOURCE=$2; shift 2 ;;
        --dry-run) DRY=1; shift ;;
        -h|--help) sed -n '2,16p' "$0"; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done

HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$HERE/../.." && pwd)          # .../P2_basket_analysis
LX=${LX:-akallits@lxplus.cern.ch}
LXDIR=${LXDIR:-p2_condor}
SSH="ssh -o GSSAPIAuthentication=yes -o GSSAPIDelegateCredentials=yes -o ConnectTimeout=30"
SCP="scp -o GSSAPIAuthentication=yes -o GSSAPIDelegateCredentials=yes -o ConnectTimeout=30"
export KRB5_CONFIG=${KRB5_CONFIG:-/local/home/banco/DAQ_Control_Dream_Beam/config/krb5_cern.conf}

klist -s 2>/dev/null || { echo "No Kerberos ticket -- run: kinit akallits@CERN.CH" >&2; exit 1; }

# --- 1. pack the code ------------------------------------------------------
echo "== packing code from $REPO"
tar --exclude=__pycache__ --exclude='*.pyc' --exclude=logbook \
    --exclude=patches --exclude=condor \
    -czf "$HERE/p2code.tgz" \
    -C "$REPO" sps_beam_analysis cosmic_bench_analysis \
    Detector_Mapping/P2_BASKET/P2_BASKET_mapping.csv
echo "   p2code.tgz  $(du -h "$HERE/p2code.tgz" | cut -f1)"

# --- 2. build the joblist from EOS -----------------------------------------
echo "== building joblist (group=$GROUP source=$SOURCE)"
python3 "$HERE/make_joblist.py" --group "$GROUP" --source "$SOURCE" \
    ${RUNS:+--run "$RUNS"} -o "$HERE/joblist.txt"
N=$(wc -l < "$HERE/joblist.txt")
[ "$N" -gt 0 ] || { echo "joblist is empty -- nothing to submit" >&2; exit 1; }

if [ "$DRY" = "1" ]; then
    echo "== dry run: $N job(s) would be submitted; joblist at $HERE/joblist.txt"
    exit 0
fi

# --- 3. ship and submit ----------------------------------------------------
echo "== shipping to $LX:~/$LXDIR"
$SSH "$LX" "mkdir -p ~/$LXDIR/logs ~/$LXDIR/products"
$SCP "$HERE/job.sh" "$HERE/analysis.sub" "$HERE/p2code.tgz" \
     "$HERE/joblist.txt" "$LX:~/$LXDIR/"

echo "== submitting $N job(s)"
$SSH "$LX" "cd ~/$LXDIR && chmod +x job.sh && condor_submit analysis.sub"

cat <<EOF

Submitted. To watch:
  ssh -o GSSAPIDelegateCredentials=yes $LX 'condor_q'
When it is done, merge the scan-level products and pull them to the GUI:
  $HERE/merge_and_pull.sh --group $GROUP
EOF
