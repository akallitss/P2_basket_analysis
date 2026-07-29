#!/usr/bin/env bash
#
# Run align_and_track.py over every sub-run of one or more runs.
#
# Products go to $SPS_ANALYSIS_ROOT/urw_alignment/<run>/, NOT into the repo -
# one log, plot and alignment json per sub-run, plus <run>/summary.tsv.  The
# combined <outdir>/summary.tsv is rebuilt from all runs present.  -o overrides.
#
# The uRWELLs sat at a fixed operating point for the whole campaign (drift 600,
# resistive 420), so the alignment constants should come out the same in every
# sub-run of a scan - the scans vary P2 voltages, not uRWELL ones.  A sub-run
# whose dx/dy differs from the others is telling you something is wrong with
# that sub-run, not that the detector moved.
#
# Usage:
#   ./loop_subruns.sh                                  # drift_mesh_scan_1, all sub-runs
#   ./loop_subruns.sh -r highstat_eff_1                # another run
#   ./loop_subruns.sh -r drift_mesh_scan_1 -r highstat_eff_1
#   ./loop_subruns.sh -o /tmp/scan -m 5e5              # quick pass, fewer hits
#   ./loop_subruns.sh -t                               # also save per-event tracks
#   ./loop_subruns.sh -n                               # dry run, just list
#
set -u -o pipefail

# --- environment. Both lines matter (handoff §2): the banco login shell puts an
# --- ISEG HV SDK on PYTHONPATH whose partial module copies shadow ROOT/uproot.
unset PYTHONPATH
PY=/local/home/banco/DAQ_Control_Dream_Beam/.venv/bin/python

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
RUNS_DIR=/local/home/banco/P2_data/TB_July2026_H4/runs

RUNS=()
# results belong under the analysis root, never in the git working tree
ANALYSIS_ROOT="${SPS_ANALYSIS_ROOT:-/local/home/banco/P2_data/TB_July2026_H4/analysis}"
OUTDIR="$ANALYSIS_ROOT/urw_alignment"
MIN_AMP=0            # handoff §7: an amplitude cut throws away 97% of the
                     # four-fold coincidences and makes the correlation worse
MAX_HITS=3e6
DRY=0
SKIP_INCOMPLETE=1
PLOTS=1
TRACKS=0            # per-event tracks are ~100 MB a sub-run; off unless asked

usage() { sed -n '2,21p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }

while getopts ":r:o:m:M:ntAPh" opt; do
    case "$opt" in
        r) RUNS+=("$OPTARG") ;;
        o) OUTDIR="$OPTARG" ;;
        m) MAX_HITS="$OPTARG" ;;
        M) MIN_AMP="$OPTARG" ;;
        n) DRY=1 ;;
        t) TRACKS=1 ;;
        A) SKIP_INCOMPLETE=0 ;;
        P) PLOTS=0 ;;
        h) usage 0 ;;
        \?) echo "unknown option -$OPTARG" >&2; usage 1 ;;
    esac
done
[ ${#RUNS[@]} -eq 0 ] && RUNS=(drift_mesh_scan_1)

[ -x "$PY" ] || { echo "interpreter not found: $PY" >&2; exit 1; }

HEADER='run\tsub_run\tn_matched\tslope_x\tdx_mm\tspread_x\tslope_y\tdy_mm\tspread_y\tangle_x_mrad\tangle_y_mrad\tseconds\tstatus\n'

# The summary is written per run, and the combined <outdir>/summary.tsv is
# rebuilt from those at the end.  Doing it the other way round meant that
# running one run silently truncated another run's rows out of the file.
# A dry run must leave nothing behind, so it creates no directory at all.
if [ "$DRY" -eq 1 ]; then
    SUMMARY=/dev/null
else
    mkdir -p "$OUTDIR"
    SUMMARY="$OUTDIR/summary.tsv"
fi

n_ok=0; n_fail=0; n_skip=0
t_all=$(date +%s)

for run in "${RUNS[@]}"; do
    run_dir="$RUNS_DIR/$run"
    if [ ! -d "$run_dir" ]; then
        echo "no such run: $run_dir" >&2; n_fail=$((n_fail+1)); continue
    fi

    # <analysis root>/urw_alignment/<run>/ - the layout the rest of the
    # analysis tree uses.  Its summary starts fresh on every re-run of THIS run.
    dest="$OUTDIR/$run"
    RUN_TSV=/dev/null
    if [ "$DRY" -eq 0 ]; then
        mkdir -p "$dest"
        RUN_TSV="$dest/summary.tsv"
        printf "$HEADER" > "$RUN_TSV"
    fi

    for sub_dir in "$run_dir"/*/; do
        sub=$(basename "$sub_dir")

        # a sub-run is usable only if the DAQ finished it and the hits exist
        if [ ! -d "$sub_dir/combined_hits_root" ] || \
           [ -z "$(ls -A "$sub_dir/combined_hits_root" 2>/dev/null)" ]; then
            echo "-- skip $run/$sub: no combined_hits_root"
            n_skip=$((n_skip+1)); continue
        fi
        if [ "$SKIP_INCOMPLETE" -eq 1 ] && [ ! -f "$sub_dir/.subrun_complete" ]; then
            echo "-- skip $run/$sub: no .subrun_complete marker (-A to force)"
            n_skip=$((n_skip+1)); continue
        fi

        if [ "$DRY" -eq 1 ]; then
            nch=$(ls -1 "$sub_dir"/combined_hits_root/*.root 2>/dev/null | wc -l)
            echo "-- would run $run/$sub ($nch chunks)"
            continue
        fi

        log="$dest/${sub}.log"
        extra=()
        [ "$PLOTS" -eq 1 ] && extra+=(--plot "$dest/${sub}.png")
        # .npz, not .parquet: neither pyarrow nor fastparquet is installed in
        # the venv, so a parquet name fails after the whole sub-run is done
        [ "$TRACKS" -eq 1 ] && extra+=(--save-tracks "$dest/${sub}_tracks.npz")

        echo "== $run/$sub"
        t0=$(date +%s)
        "$PY" "$HERE/align_and_track.py" \
            --run-dir "$run_dir" --sub-run "$sub" \
            --min-amp "$MIN_AMP" --max-hits "$MAX_HITS" \
            "${extra[@]}" \
            > "$log" 2>&1
        rc=$?
        t1=$(date +%s); dt=$((t1-t0))

        if [ $rc -ne 0 ]; then
            echo "   FAILED (rc=$rc) - see $log"
            tail -3 "$log" | sed 's/^/     /'
            printf '%s\t%s\t\t\t\t\t\t\t\t\t\t%d\tFAILED\n' "$run" "$sub" "$dt" >> "$RUN_TSV"
            n_fail=$((n_fail+1)); continue
        fi

        # align_and_track.py writes this next to itself; move it out of the
        # working tree into the analysis root with the rest of the products
        json="$HERE/alignment_${run}_${sub}.json"
        if [ -f "$json" ]; then
            mv -f "$json" "$dest/alignment_${sub}.json"
            json="$dest/alignment_${sub}.json"
            "$PY" - "$run" "$sub" "$json" "$dt" <<'EOF' >> "$RUN_TSV"
import json, sys
run, sub, path, dt = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
d = json.load(open(path))
a, g = d['align'], d['angles']
print('\t'.join([run, sub, str(d['n_matched']),
                 f"{a['x']['slope']:.4f}", f"{a['x']['offset']:.3f}", f"{a['x']['rms']:.2f}",
                 f"{a['y']['slope']:.4f}", f"{a['y']['offset']:.3f}", f"{a['y']['rms']:.2f}",
                 f"{g['x']['median']:.2f}", f"{g['y']['median']:.2f}", dt, 'ok']))
EOF
            grep -E 'matched on eventId|^x:|^y:' "$log" | sed 's/^/   /'
        else
            echo "   ran but no alignment json appeared - see $log"
            printf '%s\t%s\t\t\t\t\t\t\t\t\t\t%d\tNO_JSON\n' "$run" "$sub" "$dt" >> "$RUN_TSV"
        fi
        echo "   ${dt}s"
        n_ok=$((n_ok+1))
    done
done

echo
echo "=================================================================="
echo "ok=$n_ok  failed=$n_fail  skipped=$n_skip   total $(( $(date +%s) - t_all ))s"
if [ "$DRY" -eq 0 ]; then
    # rebuild the combined table from every run present under $OUTDIR, so that
    # re-running one run does not drop the others out of it
    { printf "$HEADER"
      for f in "$OUTDIR"/*/summary.tsv; do
          [ -e "$f" ] && tail -n +2 "$f"
      done
    } > "$SUMMARY"
    echo "outputs in $OUTDIR"
    echo
    column -t -s $'\t' "$SUMMARY"
fi
[ "$n_fail" -eq 0 ]
