#!/bin/bash
# ---------------------------------------------------------------------------
# job.sh -- one HTCondor job = one (run, sub_run) of the SPS P2 analysis.
#
#   ./job.sh <run> <sub_run> <stage_group>
#
# stage_group:
#   probe   diagnostics + 21 + 28 only (Step 0 shakedown, see README_CONDOR.md)
#   hits    21 -> 22 -> 20 -> 23 -> 26 -> 28   (input: combined_hits_root)
#   wave    29                                 (input: decoded_root, ~370 MB)
#
# The job is self-contained: it stages its inputs from EOS into the condor
# scratch dir, rebuilds the DAQ-side directory layout the stages expect, runs
# them with SPS_DATA_ROOT/SPS_ANALYSIS_ROOT pointed at that scratch, and pushes
# the (small) products back to EOS. Nothing is written to AFS.
#
# Why stage in rather than read over xrootd directly: the stages use plain
# os.path/glob on a local tree (sps_config.RunConfig), and 21 -> 22 within a
# sub-run re-read the same hits file several times. One sequential 400 MB pull
# beats many random remote reads.
# ---------------------------------------------------------------------------
set -uo pipefail

RUN=${1:?usage: job.sh <run> <sub_run> <stage_group> [prod_sub]}
SUB=${2:?usage: job.sh <run> <sub_run> <stage_group> [prod_sub]}
GROUP=${3:-hits}
# Where stages 20/22 file their per-point products; make_joblist.py sets it to
# 'scan' for a multi-point run. '-' or empty means "let the stage decide".
PROD_SUB=${4:-}
[ "$PROD_SUB" = "-" ] && PROD_SUB=""
# Which EOS holds this run's DATA. The campaign was backed up to three
# different places while quotas were being sorted out, and the early runs only
# ever reached salsachip.
SRC=${5:-ntof}

case "$SRC" in
    ntof)
        EOS_URL=root://eospublic.cern.ch
        EOS_BASE=/eos/experiment/ntof/data/x17/p2_sps_july ;;
    salsachip)
        EOS_URL=root://eosproject.cern.ch
        EOS_BASE=/eos/project/s/salsachip/Data/T2_tests/P2_SPS_Dream_Data ;;
    user)
        EOS_URL=root://eosuser.cern.ch
        EOS_BASE=/eos/user/a/akallits/P2_SPS_backup_temp ;;
    *)  echo "FATAL: unknown source '$SRC'" >&2; exit 1 ;;
esac

# Output ALWAYS goes to nTOF, whatever the input source. salsachip has been
# over quota since 2026-07-25 (every write there returns [3021]) and the user
# EOS space was only ever a stopgap, so those are read-only as far as this
# pipeline is concerned. Keeping one output location also means the products of
# all three sources land in a single tree that merge_and_pull.sh can pull in
# one pass.
OUT_URL=root://eospublic.cern.ch
OUT_BASE=/eos/experiment/ntof/data/x17/p2_sps_july
LCG_VIEW=${LCG_VIEW:-/cvmfs/sft.cern.ch/lcg/views/LCG_110}

SCRATCH=${_CONDOR_SCRATCH_DIR:-$PWD}
SPS_ANALYSIS_ROOT="$SCRATCH/analysis"
T0=$(date +%s)

say() { echo "[$(date -u +%H:%M:%S)] $*"; }
fail() { echo "FATAL: $*" >&2; exit 1; }

# --- 0. products.tgz must exist from the first second of the job -------------
# transfer_output_files names it unconditionally, and a condor job whose named
# output file is missing goes ON HOLD rather than returning its logs -- which
# is exactly backwards when the job died early and the logs are the whole
# point. So: create it empty now, and repack whatever exists on the way out,
# however we exit.
# UPLOAD_OK is set to 1 once every product has reached EOS. In that case the
# sandbox copy is packed EMPTY (~45 bytes): the submit file still gets the file
# it asked for, but the AFS home -- which has well under 1 GB free -- does not
# collect a redundant ~2 MB per sub_run on top of what is already on EOS.
UPLOAD_OK=0
pack_products() {
    if [ "$UPLOAD_OK" = "1" ]; then
        tar -czf "$SCRATCH/products.tgz" -T /dev/null 2>/dev/null
        return
    fi
    tar -czf "$SCRATCH/products.tgz" -C "$SPS_ANALYSIS_ROOT" . 2>/dev/null \
        || tar -czf "$SCRATCH/products.tgz" -T /dev/null 2>/dev/null
}
mkdir -p "$SPS_ANALYSIS_ROOT"
pack_products
trap pack_products EXIT

say "=== job.sh  run=$RUN  sub_run=$SUB  group=$GROUP  src=$SRC"
say "host=$(hostname -f)  user=$(whoami)  scratch=$SCRATCH"
df -h "$SCRATCH" | tail -1

# --- 1. Kerberos ------------------------------------------------------------
# HTCondor with MY.SendCredential=true drops a krb5 ccache in $_CONDOR_CREDS.
# Reads from the nTOF space are world-readable so they work unauthenticated,
# but the write-back leg needs us to be akallits.
if [ -n "${_CONDOR_CREDS:-}" ] && [ -d "${_CONDOR_CREDS}" ]; then
    say "condor creds dir: $_CONDOR_CREDS"
    ls -la "$_CONDOR_CREDS" || true
    cc=$(ls "$_CONDOR_CREDS"/*.cc 2>/dev/null | head -1)
    [ -n "$cc" ] && export KRB5CCNAME="FILE:$cc"
fi
say "KRB5CCNAME=${KRB5CCNAME:-<unset>}"
klist 2>&1 | head -8 || say "klist: no ticket"

# --- 2. Software environment ------------------------------------------------
# LCG views carry numpy/pandas/matplotlib/scipy/uproot; nothing to build.
PLAT=""
for p in x86_64-el9-gcc13-opt x86_64-el9-gcc14-opt x86_64-el9-gcc12-opt \
         x86_64-centos9-gcc12-opt x86_64-el8-gcc13-opt; do
    if [ -f "$LCG_VIEW/$p/setup.sh" ]; then PLAT=$p; break; fi
done
[ -n "$PLAT" ] || fail "no LCG platform found under $LCG_VIEW"
say "LCG view: $LCG_VIEW/$PLAT"
# The LCG setup.sh reads variables it never sets (COMPILER, ...), so it dies
# instantly under `set -u`. Drop the flag across the source and put it back.
set +u
# shellcheck disable=SC1090
source "$LCG_VIEW/$PLAT/setup.sh" || fail "LCG setup.sh failed"
set -u

python3 - <<'PY' || fail "python dependency check failed"
import sys
mods = ['numpy', 'pandas', 'matplotlib', 'scipy', 'uproot', 'awkward']
print('python', sys.version.split()[0])
bad = []
for m in mods:
    try:
        mod = __import__(m)
        print(f'  {m:12s} {getattr(mod, "__version__", "?")}')
    except Exception as e:
        print(f'  {m:12s} MISSING ({e})'); bad.append(m)
sys.exit(1 if bad else 0)
PY

# matplotlib must not try to open a display or write to a read-only HOME.
export MPLBACKEND=Agg
export MPLCONFIGDIR="$SCRATCH/.mpl"
export HOME="$SCRATCH"          # keep any stray dotfile writes off AFS
mkdir -p "$MPLCONFIGDIR"

# --- 3. Unpack the analysis code -------------------------------------------
CODE="$SCRATCH/code"
mkdir -p "$CODE"
tar -xzf p2code.tgz -C "$CODE" || fail "could not unpack p2code.tgz"
[ -f "$CODE/sps_beam_analysis/sps_config.py" ] || fail "code payload incomplete"
say "code unpacked: $(du -sh "$CODE" | cut -f1)"

# --- 3b. 'rec' group: recorded-trigger sets, read straight off EOS ----------
# 22_tag_probe_efficiency needs <sub_run>/recorded_events.npz to correct for
# DAQ overlap; without it a trigger the probe FEU never recorded is
# indistinguishable from one where it recorded nothing, and the efficiency is
# biased low. That npz is made from decoded_root -- but ONLY from its eventId
# branch, so there is no reason to move the ~370 MB: uproot reads that one
# column over xrootd in well under a second (measured 2.4 M eventIds in 0.5 s).
# This group therefore stages nothing and writes its result back beside the
# data, where a later `hits` sweep picks it up automatically.
if [ "$GROUP" = "rec" ]; then
    DEC="$EOS_BASE/runs/$RUN/$SUB/decoded_root"
    mapfile -t URLS < <(xrdfs "$EOS_URL" ls "$DEC" 2>/dev/null \
                        | grep '\.root$' | grep -v '/\.' | sed "s#^#$EOS_URL/#")
    [ ${#URLS[@]} -gt 0 ] || fail "no decoded .root files in $DEC"
    say "extracting recorded triggers from ${#URLS[@]} chunk(s), no staging"
    NPZ="$SCRATCH/recorded_events.npz"
    cd "$CODE/sps_beam_analysis" || fail "cannot cd into the code dir"
    python3 extract_recorded_events.py --files "${URLS[@]}" --out "$NPZ" --force \
        || fail "extract_recorded_events failed"
    [ -s "$NPZ" ] || fail "extraction produced no output"
    # Written to the nTOF tree, NOT beside the source data: salsachip is over
    # quota and cannot be written to at all. `hits` looks here first.
    DEST="$OUT_BASE/runs/$RUN/$SUB/recorded_events.npz"
    xrdfs "$OUT_URL" mkdir -p "$(dirname "$DEST")" 2>/dev/null
    xrdcp -f -s "$NPZ" "$OUT_URL/$DEST" || fail "upload failed: $DEST"
    say "uploaded $(du -h "$NPZ" | cut -f1) -> $DEST"
    say "=== done  rc=0  total $(( $(date +%s) - T0 ))s"
    exit 0
fi

# --- 4. Stage inputs from EOS ----------------------------------------------
# Rebuild <DATA_ROOT>/<run>/<sub_run>/... exactly as the DAQ writes it, so the
# stages need no path knowledge of EOS at all.
DATA="$SCRATCH/data"
SRC_URL="$EOS_URL/$EOS_BASE/runs/$RUN"
DST="$DATA/$RUN/$SUB"
mkdir -p "$DST/raw_daq_data"

xrd_get() {  # xrd_get <remote-rel-path> <local-path>   (missing file is OK)
    xrdcp -f -s "$SRC_URL/$1" "$2" 2>/dev/null && return 0
    say "  (absent: $1)"; return 1
}

say "staging inputs from $SRC_URL"
xrd_get "run_config.json" "$DATA/$RUN/run_config.json" \
    || fail "run_config.json missing -- cannot build a RunConfig"
xrd_get "$SUB/hv_monitor.csv"            "$DST/hv_monitor.csv"
xrd_get "$SUB/raw_daq_data/run_time.txt" "$DST/raw_daq_data/run_time.txt"
# recorded_events.npz is written by the `rec` group into the nTOF tree
# regardless of where the data came from, so look there FIRST and fall back to
# beside-the-data for anything produced before that split existed.
xrdcp -f -s "$OUT_URL/$OUT_BASE/runs/$RUN/$SUB/recorded_events.npz" \
      "$DST/recorded_events.npz" 2>/dev/null \
    || xrd_get "$SUB/recorded_events.npz" "$DST/recorded_events.npz"

case "$GROUP" in
    wave) NEED=decoded_root ;;
    *)    NEED=combined_hits_root ;;
esac
# `raweff` is the one RUN-level group: stage 30 walks every sub_run of the run
# and reads a capped slice of ONE decoded chunk per station, so staging is both
# unnecessary and impossible -- a whole run's decoded_root is tens of GB against
# an 8 GB sandbox. It reads over xrootd instead and needs only run_config.json,
# which is already down. Skip straight to the stage.
if [ "$GROUP" != "raweff" ]; then
mkdir -p "$DST/$NEED"
say "staging $NEED ..."
# `grep -v '/\.'` drops EOS's atomic-version placeholders, which are named
# .sys.v#.<realfile>.root and so match *.root. They are not copyable files:
# xrdcp fails on them, and because staging treats a failed copy as fatal, one
# placeholder killed the entire job. This cost 29 of 266 sub_runs on the first
# full sweep (2026-07-28) -- all on salsachip, which has versioning enabled.
mapfile -t FILES < <(xrdfs "$EOS_URL" ls "$EOS_BASE/runs/$RUN/$SUB/$NEED" 2>/dev/null \
                     | grep '\.root$' | grep -v '/\.')
[ ${#FILES[@]} -gt 0 ] || fail "no .root files in $NEED for $RUN/$SUB"
for f in "${FILES[@]}"; do
    xrdcp -f -s "$EOS_URL/$f" "$DST/$NEED/$(basename "$f")" \
        || fail "xrdcp failed for $f"
done
say "staged: $(du -sh "$DST/$NEED" | cut -f1) in ${#FILES[@]} file(s), $(( $(date +%s) - T0 ))s"
fi

# --- 5. Run the stages ------------------------------------------------------
export SPS_DATA_ROOT="$DATA"
export SPS_ANALYSIS_ROOT      # set at the top, so the EXIT trap can pack it
export SPS_RUN="$RUN"
cd "$CODE/sps_beam_analysis" || fail "cannot cd into the code dir"

case "$GROUP" in
    probe) STAGES=(21_telescope_align.py 28_timing_qa.py) ;;
    hits)  STAGES=(21_telescope_align.py 22_tag_probe_efficiency.py
                   20_beam_spectra.py 23_beam_profile.py
                   24_event_sync_qa.py
                   26_hv_spark_qa.py 28_timing_qa.py) ;;
    wave)  STAGES=(29_waveform_timing.py) ;;
    raweff) STAGES=(30_raw_stream_efficiency.py) ;;
    *)     fail "unknown stage group '$GROUP'" ;;
esac

# Exit codes, which analysis.sub's retry policy depends on:
#   0  every stage succeeded
#   1  infrastructure failure (staging, LCG, missing input) -- RETRY helps
#   2  the job ran but >=1 stage failed -- deterministic, retrying just burns
#      another 10 minutes to reach the same place. `fail()` above exits 1.
RC=0
failed_stages=()
for s in "${STAGES[@]}"; do
    # Only 20 and 22 write per-point products that belong together under a
    # scan directory; the rest already file per sub_run.
    extra=()
    case "$s" in
        20_beam_spectra.py|22_tag_probe_efficiency.py)
            [ -n "$PROD_SUB" ] && extra=(--prod-sub "$PROD_SUB") ;;
        30_raw_stream_efficiency.py)
            # run-level: every sub_run at once, read straight off EOS
            extra=(--eos-url "$EOS_URL" --eos-base "$EOS_BASE/runs") ;;
    esac
    say "--- $s ${extra[*]:-}"
    t=$(date +%s)
    if [ "$GROUP" = "raweff" ]; then
        python3 "$s" live ${extra[@]+"${extra[@]}"}
    else
        python3 "$s" live --sub-run "$SUB" ${extra[@]+"${extra[@]}"}
    fi
    r=$?
    say "--- $s exit=$r  $(( $(date +%s) - t ))s"
    if [ $r -ne 0 ]; then
        RC=2
        failed_stages+=("$s($r)")
    fi
done
# One greppable line per job, so a 169-job sweep can be triaged with
#   grep -h 'STAGE FAILURES' logs/*.out
if [ ${#failed_stages[@]} -gt 0 ]; then
    say "STAGE FAILURES: $RUN/$SUB: ${failed_stages[*]}"
fi

# --- 6. Push products back to EOS ------------------------------------------
# Products are PNG/CSV/JSON (a few MB per sub_run). The layout mirrors the DAQ
# GUI's Analysis tab: <analysis>/<det_tag>/<run>/<sub_run>/<stage>/...
OUT_EOS="$OUT_BASE/analysis"
n_files=$(find "$SPS_ANALYSIS_ROOT" -type f 2>/dev/null | wc -l)
n_up=0; n_fail=0

# ONE recursive xrdcp per top-level station dir, not one per file: every xrdcp
# invocation opens a fresh connection and pays a Kerberos handshake, so this is
# 3-4 round trips instead of ~50. Not a huge win in absolute terms -- the
# per-file fallback measured 50 files in 6 s on 2026-07-28 -- but it is the
# difference that matters when a sub_run produces per-pad maps by the hundred.
upload_recursive() {
    local ok=1
    xrdfs "$OUT_URL" mkdir -p "$OUT_EOS" 2>/dev/null
    for d in "$SPS_ANALYSIS_ROOT"/*/; do
        [ -d "$d" ] || continue
        # Skip station dirs that contain no files. out_dir() mkdir -p's a tree
        # for every detector before it knows whether there is anything to put
        # in it, so the uRWELL references (no pad map, no products) leave empty
        # trees behind -- and xrdcp -r on an empty tree fails, which used to
        # drag the whole upload into the per-file fallback path.
        [ -n "$(find "$d" -type f -print -quit 2>/dev/null)" ] || continue
        local tag t
        tag=$(basename "${d%/}"); t=$(date +%s)
        if xrdcp -r -f -s "${d%/}" "$OUT_URL/$OUT_EOS/" 2>/dev/null; then
            say "  uploaded $tag/ ($(( $(date +%s) - t ))s)"
        else
            say "  recursive upload failed for $tag/"
            ok=0
        fi
    done
    return $((1 - ok))
}

# Per-file fallback: slow, but it is the path that is already proven to work,
# so a quirk in recursive xrdcp cannot cost us a whole sweep's results.
upload_per_file() {
    while IFS= read -r f; do
        local rel=${f#"$SPS_ANALYSIS_ROOT"/}
        xrdfs "$OUT_URL" mkdir -p "$OUT_EOS/$(dirname "$rel")" 2>/dev/null
        if xrdcp -f -s "$f" "$OUT_URL/$OUT_EOS/$rel" 2>/dev/null; then
            n_up=$((n_up + 1))
        else
            n_fail=$((n_fail + 1)); say "  UPLOAD FAILED: $rel"
        fi
    done < <(find "$SPS_ANALYSIS_ROOT" -type f)
}

if [ "$n_files" -eq 0 ]; then
    say "no products to upload"
    n_fail=1
elif upload_recursive; then
    n_up=$n_files
else
    say "falling back to per-file upload"
    upload_per_file
fi
say "uploaded $n_up/$n_files file(s) to $OUT_EOS, $n_fail failure(s)"
[ "$n_fail" -eq 0 ] && [ "$n_up" -gt 0 ] && UPLOAD_OK=1

# The EXIT trap repacks products.tgz on the way out, so the products also come
# home through condor's own output sandbox whatever happens above.
[ "$n_fail" -gt 0 ] && say "NOTE: EOS write-back partially failed -- rely on products.tgz"

say "=== done  rc=$RC  total $(( $(date +%s) - T0 ))s"
exit $RC
