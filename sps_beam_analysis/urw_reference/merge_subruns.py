#!/usr/bin/env python3
"""
Merge the per-sub-run `combined_hits_root` files of TB_July2026_H4 into one
ROOT file, for uRWELL + P2 track reconstruction across all sub-runs at once.

Why this script exists rather than `hadd`
-----------------------------------------
1. `hadd` is installed here as a snap and refuses to touch /local/home
   ("Sorry, home directories outside of /home needs configuration"), so it
   cannot read the data or write the output.  We merge with uproot instead.

2. `eventId` restarts at 1 in every sub-run.  Verified on disk:
       drift_mesh_scan_1/nominal_00   eventId 1 .. 2795239
       drift_mesh_scan_1/drift_450    eventId 1 .. 2792733
   so a plain concatenation would collide on the one branch everything is
   joined by.  This script adds two columns:

       subrun_id  int32   index into the `subruns` lookup tree written here
       geventId   uint64  subrun_id * GEVENT_STRIDE + eventId

   **Join on `geventId`, never on `eventId` alone.**  Within one sub-run
   geventId is monotonic in eventId, so it also sorts sensibly.

3. The source files carry several *cycles* of the `hits` tree from
   reprocessing - `hits;21`, `hits;20` (handoff §5.3).  `uproot.open(f)['hits']`
   resolves to the highest cycle, which is the newest and correct one.  Cycles
   are never summed; `hadd` would have had to be told about this explicitly.

What is safe here and how it was checked
----------------------------------------
Chunk files concatenate without duplication.  Per FEU the eventId ranges of
consecutive chunks are strictly disjoint - on highstat_eff_1/beam_commissioning_00,
FEU 4 ends chunk 002 at 7283632 and begins chunk 003 at 7283633.  The
whole-file ranges *look* like they overlap (002 spans ..7283648 and 003 spans
7283633..) only because the four FEUs stop at different events.  Nothing is
double counted and no event is split across two chunk files.

The output tree keeps the source's FEU-major ordering within each sub-run
(handoff §5.4 trap: it is ordered by FEU, not by event).  Do not truncate a
read of this file and conclude a FEU is empty.

Usage
-----
    unset PYTHONPATH
    PY=/local/home/banco/DAQ_Control_Dream_Beam/.venv/bin/python
    $PY merge_subruns.py --runs drift_mesh_scan_1 highstat_eff_1 \\
        --out /local/home/banco/P2_data/TB_July2026_H4/merged/all_subruns_hits.root

    $PY merge_subruns.py --runs drift_mesh_scan_1 --dry-run     # just list
"""
import os
import sys
import glob
import json
import time
import argparse

import numpy as np

RUNS_DIR = '/local/home/banco/P2_data/TB_July2026_H4/runs'

# subrun_id * GEVENT_STRIDE + eventId.  The largest eventId seen in this
# campaign is 8_384_316 (highstat_eff_1/beam_commissioning_00), four orders of
# magnitude below the stride, and the script asserts the margin per sub-run.
GEVENT_STRIDE = 10 ** 10

# The branches urw_lib actually reads, plus `feu` to tell the detectors apart.
# The rest of the analyser's fit outputs are not used for tracking (handoff §5.3).
TRACKING_BRANCHES = ['eventId', 'trigger_timestamp_ns', 'channel', 'amplitude',
                     'time_of_max', 'significance', 'saturated', 'feu']

FEU_DETECTOR = {1: 'EIC_uRWELL_front + EIC_uRWELL_back',
                3: 'P2_IN', 4: 'P2_MID', 5: 'P2_OUT'}


def discover_subruns(runs, runs_dir, include_incomplete=False):
    """Every sub-run of `runs` that has combined_hits_root files, in order.

    Sub-runs without a `.subrun_complete` marker are skipped unless asked for:
    the marker is written only when the DAQ finished the sub-run cleanly.
    """
    found = []
    for run in runs:
        run_dir = os.path.join(runs_dir, run)
        if not os.path.isdir(run_dir):
            sys.exit(f'no such run directory: {run_dir}')
        for sub in sorted(os.listdir(run_dir)):
            sub_dir = os.path.join(run_dir, sub)
            if not os.path.isdir(sub_dir):
                continue
            files = sorted(glob.glob(os.path.join(
                sub_dir, 'combined_hits_root', '*_feu-combined_hits.root')))
            complete = os.path.exists(os.path.join(sub_dir, '.subrun_complete'))
            if not files:
                print(f'  skip {run}/{sub}: no combined_hits_root files')
                continue
            if not complete and not include_incomplete:
                print(f'  skip {run}/{sub}: no .subrun_complete marker')
                continue
            found.append(dict(run=run, sub_run=sub, sub_dir=sub_dir,
                              files=files, complete=complete))
    return found


def missing_events(sub_dir):
    """Per-FEU count of events the readout dropped, from recorded_events.npz.

    Handoff §5.5: check this before trusting a sub-run.  It is empty for
    nominal_00 but FEU 1 dropped ~632 000 events in drift_450.  Missing events
    are a readout problem, not a physics one, but they skew any efficiency.
    """
    path = os.path.join(sub_dir, 'recorded_events.npz')
    if not os.path.exists(path):
        return {}
    try:
        with np.load(path, allow_pickle=True) as d:
            return {k: int(np.asarray(d[k]).size)
                    for k in d.files if k.endswith('_missing')}
    except Exception as exc:                                   # noqa: BLE001
        print(f'    warning: cannot read {path}: {exc}')
        return {}


def merge(subruns, out_path, branches, step, feu_filter, compression_level):
    """Stream every chunk file into one `hits` TTree, adding the sub-run columns.

    Note on uproot: `fout['hits'] = {...}` writes an **RNTuple**, not a TTree,
    in uproot 5.7.5.  RNTuple is ROOT's new experimental format and old ROOT
    and C++ TTree-based code cannot read it, so the tree is created with
    `mktree`, which does produce a genuine TTree (Model_TTree_v20 on readback).
    """
    import uproot

    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    summary = []
    n_total = 0
    started = time.time()

    # dtypes must be pinned up front: mktree fixes the branch types and every
    # later extend() has to match them exactly
    src = uproot.open(subruns[0]['files'][0])['hits']
    types = {b.name: src[b.name].interpretation.numpy_dtype
             for b in src.branches if b.name in branches}
    types['subrun_id'] = np.dtype(np.int32)
    types['geventId'] = np.dtype(np.uint64)

    with uproot.recreate(out_path,
                         compression=uproot.ZLIB(compression_level)) as fout:
        hits = fout.mktree('hits', types)
        for sid, entry in enumerate(subruns):
            tag = f'{entry["run"]}/{entry["sub_run"]}'
            n_sub = 0
            ev_min, ev_max = None, None
            per_feu = {}
            t0 = time.time()
            print(f'[{sid + 1}/{len(subruns)}] {tag}  '
                  f'({len(entry["files"])} chunks)', flush=True)

            for path in entry['files']:
                tree = uproot.open(path)['hits']       # highest cycle, never summed
                for chunk in tree.iterate(branches, library='np',
                                          step_size=step):
                    if feu_filter is not None:
                        sel = np.isin(chunk['feu'], feu_filter)
                        chunk = {k: v[sel] for k, v in chunk.items()}
                    ev = chunk['eventId']
                    if not len(ev):
                        continue

                    lo, hi = int(ev.min()), int(ev.max())
                    if hi >= GEVENT_STRIDE:
                        sys.exit(f'eventId {hi} in {tag} exceeds the '
                                 f'geventId stride {GEVENT_STRIDE}')
                    ev_min = lo if ev_min is None else min(ev_min, lo)
                    ev_max = hi if ev_max is None else max(ev_max, hi)

                    out = {k: chunk[k].astype(types[k], copy=False)
                           for k in branches}
                    out['subrun_id'] = np.full(len(ev), sid, dtype=np.int32)
                    out['geventId'] = (np.uint64(sid) * np.uint64(GEVENT_STRIDE)
                                       + ev.astype(np.uint64))
                    hits.extend(out)

                    n_sub += len(ev)
                    for f, c in zip(*np.unique(chunk['feu'],
                                               return_counts=True)):
                        per_feu[int(f)] = per_feu.get(int(f), 0) + int(c)

            n_total += n_sub
            dt = time.time() - t0
            print(f'      {n_sub:,} hits  eventId {ev_min}..{ev_max}  '
                  f'{dt:.0f}s  ({n_sub / max(dt, 1e-9) / 1e6:.1f} Mhit/s)',
                  flush=True)

            summary.append(dict(
                subrun_id=sid, run=entry['run'], sub_run=entry['sub_run'],
                n_chunks=len(entry['files']), n_hits=n_sub,
                event_id_min=ev_min, event_id_max=ev_max,
                complete=entry['complete'],
                hits_per_feu=per_feu,
                missing_events=missing_events(entry['sub_dir'])))

        if n_total == 0:
            sys.exit('nothing was written - no hits survived the selection')

        # lookup tree: subrun_id -> which run and sub-run the hits came from
        import awkward as ak
        lookup = fout.mktree('subruns', {
            'subrun_id': np.int32, 'run': str, 'sub_run': str,
            'n_chunks': np.int32, 'n_hits': np.int64,
            'event_id_min': np.uint64, 'event_id_max': np.uint64,
            'n_missing_feu1': np.int64})
        lookup.extend({
            'subrun_id': np.array([s['subrun_id'] for s in summary], np.int32),
            'run': ak.Array([s['run'] for s in summary]),
            'sub_run': ak.Array([s['sub_run'] for s in summary]),
            'n_chunks': np.array([s['n_chunks'] for s in summary], np.int32),
            'n_hits': np.array([s['n_hits'] for s in summary], np.int64),
            'event_id_min': np.array([s['event_id_min'] for s in summary],
                                     np.uint64),
            'event_id_max': np.array([s['event_id_max'] for s in summary],
                                     np.uint64),
            'n_missing_feu1': np.array(
                [s['missing_events'].get('feu1_missing', -1) for s in summary],
                np.int64),
        })

        meta = dict(created=time.strftime('%Y-%m-%d %H:%M:%S'),
                    gevent_stride=GEVENT_STRIDE,
                    branches=list(branches) + ['subrun_id', 'geventId'],
                    feu_filter=sorted(feu_filter) if feu_filter else 'all',
                    feu_detector=FEU_DETECTOR,
                    n_hits=n_total, subruns=summary)
        fout['meta'] = json.dumps(meta, indent=2)

    with open(os.path.splitext(out_path)[0] + '_index.json', 'w') as fh:
        json.dump(meta, fh, indent=2)

    dt = time.time() - started
    size = os.path.getsize(out_path) / 1e9
    print(f'\nwrote {out_path}')
    print(f'  {n_total:,} hits from {len(subruns)} sub-runs, '
          f'{size:.2f} GB, {dt / 60:.1f} min')
    print(f'  index: {os.path.splitext(out_path)[0] + "_index.json"}')
    return n_total


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--runs', nargs='+',
                    default=['drift_mesh_scan_1', 'highstat_eff_1'],
                    help='run directories under --runs-dir to merge')
    ap.add_argument('--runs-dir', default=RUNS_DIR)
    ap.add_argument('--out', default=os.path.join(
        os.path.dirname(RUNS_DIR), 'merged', 'all_subruns_hits.root'))
    ap.add_argument('--all-branches', action='store_true',
                    help='keep every branch of the source tree, not just the '
                         'ones tracking needs')
    ap.add_argument('--feu', nargs='+', type=int, default=None,
                    help='keep only these FEUs (1=both uRWELLs, 3/4/5=P2 '
                         'IN/MID/OUT). Default: all.')
    ap.add_argument('--step', type=int, default=2_000_000,
                    help='entries read and written per basket')
    ap.add_argument('--compression', type=int, default=1,
                    help='ZLIB level. lz4 and zstd are not installed here.')
    ap.add_argument('--include-incomplete', action='store_true',
                    help='also merge sub-runs with no .subrun_complete marker')
    ap.add_argument('--dry-run', action='store_true',
                    help='list what would be merged and stop')
    args = ap.parse_args()

    print(f'runs: {", ".join(args.runs)}')
    subruns = discover_subruns(args.runs, args.runs_dir, args.include_incomplete)
    if not subruns:
        sys.exit('no sub-runs found')

    n_files = sum(len(s['files']) for s in subruns)
    n_bytes = sum(os.path.getsize(f) for s in subruns for f in s['files'])
    print(f'\n{len(subruns)} sub-runs, {n_files} chunk files, '
          f'{n_bytes / 1e9:.1f} GB of input')
    for i, s in enumerate(subruns):
        miss = missing_events(s['sub_dir'])
        flag = ''
        if miss.get('feu1_missing'):
            flag = f"   !! feu1 dropped {miss['feu1_missing']:,} events"
        print(f'  {i:3d}  {s["run"]}/{s["sub_run"]:<24} '
              f'{len(s["files"])} chunks{flag}')

    if args.dry_run:
        return

    branches = None if args.all_branches else TRACKING_BRANCHES
    if branches is None:
        import uproot
        branches = [b.name for b in uproot.open(subruns[0]['files'][0])['hits'].branches]
    if args.feu is not None and 'feu' not in branches:
        branches = branches + ['feu']
    print(f'\nbranches kept: {", ".join(branches)} (+ subrun_id, geventId)')

    merge(subruns, args.out, branches, args.step,
          set(args.feu) if args.feu else None, args.compression)


if __name__ == '__main__':
    main()
