#!/usr/bin/env bash
# merge_and_pull.sh -- after a condor sweep: bring the products home, build the
# scan-level curves, and leave everything where the DAQ GUI's Analysis tab
# reads from.
#
#   ./merge_and_pull.sh [--run RUN[,RUN...]] [--no-merge] [--no-push]
#
# Three steps, in this order for a reason:
#
#   1. PULL   EOS -> banco. The per-sub_run jobs wrote their products (and a
#             scan_row.json per point) straight to EOS from the worker nodes.
#   2. MERGE  run stages 20/22/26/28 with --scan-only, ON BANCO. Each reads the
#             scan_row.json files -- a few kB -- and rebuilds the scan-level
#             curve + CSV that a serial whole-run pass would have produced.
#             This needs NO hits files, which is the point: the raw data can
#             stay on EOS, or be pruned entirely.
#   3. PUSH   the merged scan products back to EOS, so the EOS copy is complete
#             too and a future pull gets everything.
#
# The pull writes into <analysis>/<det_tag>/... and <analysis>/telescope/...,
# never into <analysis>/<run>/..., which is qa_watcher's live territory. The
# two trees are disjoint by construction, so this cannot overwrite live QA.
set -euo pipefail

RUNS=""
DO_MERGE=1
DO_PUSH=1
while [ $# -gt 0 ]; do
    case "$1" in
        --run)      RUNS=$2; shift 2 ;;
        --no-merge) DO_MERGE=0; shift ;;
        --no-push)  DO_PUSH=0; shift ;;
        -h|--help)  sed -n '2,25p' "$0"; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done

HERE=$(cd "$(dirname "$0")" && pwd)
STAGES_DIR=$(cd "$HERE/.." && pwd)
EOS_URL=${EOS_URL:-root://eospublic.cern.ch}
EOS_BASE=${EOS_BASE:-/eos/experiment/ntof/data/x17/p2_sps_july}
BASE=${P2_BASE:-/local/home/banco/P2_data/TB_July2026_H4}
ANALYSIS=${SPS_ANALYSIS_ROOT:-$BASE/analysis}
PY=${SPS_PYTHON:-/local/home/banco/DAQ_Control_Dream_Beam/.venv/bin/python}
export KRB5_CONFIG=${KRB5_CONFIG:-/local/home/banco/DAQ_Control_Dream_Beam/config/krb5_cern.conf}

klist -s 2>/dev/null || { echo "No Kerberos ticket -- run: kinit akallits@CERN.CH" >&2; exit 1; }
mkdir -p "$ANALYSIS"

# --- 1. pull ---------------------------------------------------------------
# One recursive xrdcp per top-level station dir: a single process reusing one
# connection. Per-FILE xrdcp would pay a fresh Kerberos handshake each time and
# take minutes over a tree this size.
echo "== pulling products from EOS -> $ANALYSIS"
mapfile -t TAGS < <(xrdfs "$EOS_URL" ls "$EOS_BASE/analysis" 2>/dev/null \
                    | xargs -rn1 basename)
if [ ${#TAGS[@]} -eq 0 ]; then
    echo "   nothing under $EOS_BASE/analysis yet -- has a sweep run?" >&2
else
    for tag in "${TAGS[@]}"; do
        echo "   $tag"
        xrdcp -r -f -s "$EOS_URL/$EOS_BASE/analysis/$tag" "$ANALYSIS/" \
            || echo "   WARNING: pull failed for $tag" >&2
    done
fi

# --- 2. merge --------------------------------------------------------------
# Which runs to merge: those named, else every run that has scan rows.
if [ -n "$RUNS" ]; then
    IFS=',' read -r -a RUN_LIST <<< "$RUNS"
else
    # Layout is <tag>/<run>/<sub_run>/<stage>/scan_row.json, so the run name is
    # four fields up from the file.
    mapfile -t RUN_LIST < <(find "$ANALYSIS" -name scan_row.json 2>/dev/null \
        | awk -F/ '{print $(NF-3)}' | sort -u)
fi

if [ "$DO_MERGE" = "1" ] && [ ${#RUN_LIST[@]} -gt 0 ]; then
    # --- run_config.json mirror --------------------------------------------
    # --scan-only needs no hits files, but it DOES need run_config.json: every
    # stage calls detectors() / daq_info() / scan_axis() through it. That file
    # lives with the raw data, so for any run pruned from local disk (i.e. most
    # of the campaign) the merge died with FileNotFoundError -- on 2026-07-28
    # it worked for exactly the 6 runs still held locally and failed for the
    # other 16. Mirror it from EOS: a few kB per run, and it makes the merge
    # genuinely independent of where the data lives.
    META="$BASE/runs_meta"
    mkdir -p "$META"
    echo "== mirroring run_config.json for ${#RUN_LIST[@]} run(s) -> $META"
    for run in "${RUN_LIST[@]}"; do
        dest="$META/$run/run_config.json"
        [ -s "$dest" ] && continue
        mkdir -p "$META/$run"
        # Prefer the local copy when the run is still on disk: it is the same
        # file, and any local edit is the authoritative one.
        if [ -s "$BASE/runs/$run/run_config.json" ]; then
            cp "$BASE/runs/$run/run_config.json" "$dest"
            continue
        fi
        for spec in \
            "root://eospublic.cern.ch|/eos/experiment/ntof/data/x17/p2_sps_july" \
            "root://eosproject.cern.ch|/eos/project/s/salsachip/Data/T2_tests/P2_SPS_Dream_Data" \
            "root://eosuser.cern.ch|/eos/user/a/akallits/P2_SPS_backup_temp"; do
            u=${spec%%|*}; b=${spec#*|}
            if xrdcp -f -s "$u/$b/runs/$run/run_config.json" "$dest" 2>/dev/null; then
                break
            fi
        done
        [ -s "$dest" ] || echo "   WARNING: no run_config.json anywhere for $run" >&2
    done

    echo "== merging scan-level products for ${#RUN_LIST[@]} run(s)"
    # Point at the mirror, not $BASE/runs: the mirror has a run_config.json for
    # EVERY run, local or not.
    export SPS_DATA_ROOT="$META"
    export SPS_ANALYSIS_ROOT=$ANALYSIS
    unset PYTHONPATH || true     # banco's ISEG SDK PYTHONPATH shadows uproot
    cd "$STAGES_DIR"
    for run in "${RUN_LIST[@]}"; do
        echo "-- $run"
        export SPS_RUN=$run
        # 21/23/24/25 have no scan-level aggregation; the rest do. 29 is last
        # because it merges from its own summary JSONs rather than scan_row.json,
        # and it exits non-zero for a run the wave group has not covered yet --
        # a WARNING there means "no stage 29 products", not a failure.
        for stage in 22_tag_probe_efficiency 28_timing_qa 26_hv_spark_qa \
                     20_beam_spectra 29_waveform_timing; do
            echo "   $stage --scan-only"
            "$PY" "$stage.py" live --scan-only \
                || echo "   WARNING: $stage --scan-only failed for $run" >&2
        done
    done
else
    echo "== merge skipped"
fi

# --- 3. push the merged scan products back ---------------------------------
if [ "$DO_PUSH" = "1" ] && [ ${#TAGS[@]} -gt 0 ]; then
    echo "== pushing merged products back to EOS"
    for tag in "${TAGS[@]}"; do
        [ -d "$ANALYSIS/$tag" ] || continue
        xrdcp -r -f -s "$ANALYSIS/$tag" "$EOS_URL/$EOS_BASE/analysis/" \
            || echo "   WARNING: push failed for $tag" >&2
    done
fi

echo
echo "Done. The GUI Analysis tab (port 5001) reads $ANALYSIS directly —"
echo "no restart needed, the products are simply there now."
