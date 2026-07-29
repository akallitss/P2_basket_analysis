#!/usr/bin/env python3
"""
Check a merged file produced by merge_subruns.py against its source sub-runs.

Checks:
  1. the file holds a real TTree, not an RNTuple (uproot's dict-assignment
     writes RNTuple, which ROOT/C++ TTree code cannot read - merge_subruns.py
     uses mktree to avoid that, and this asserts it stayed that way);
  2. `subruns` lookup tree and `meta` agree with the tree length;
  3. geventId == subrun_id * 10^10 + eventId on every entry, and sub-runs
     occupy contiguous, non-overlapping entry ranges;
  4. every branch matches the source chunk files entry for entry, in order.

Note geventId is NOT unique per entry and is not meant to be: one entry is one
strip that fired, so all strips of one event share it (handoff §5.3).  It is
unique per (sub-run, event), which is what a join needs.

Usage:
    unset PYTHONPATH
    $PY verify_merged.py /path/to/merged.root [--full]

Sub-runs are stored contiguously, so the source comparison reads the merged
tree by entry range rather than with a `cut`, which would scan all 740M
entries per sub-run.  Without --full only the first and last sub-run are
compared entry by entry; --full does all of them.
"""
import os
import sys
import glob
import json
import argparse

import numpy as np
import uproot

RUNS_DIR = '/local/home/banco/P2_data/TB_July2026_H4/runs'
STRIDE = 10 ** 10


def compare_source(tree, names, sub_dir, start, stop, step):
    """True when merged entries [start, stop) equal the sub-run's source files."""
    src_names = [n for n in names if n not in ('subrun_id', 'geventId')]
    files = sorted(glob.glob(os.path.join(
        sub_dir, 'combined_hits_root', '*_feu-combined_hits.root')))
    if not files:
        return None, 'sources not found'

    src = uproot.iterate([f + ':hits' for f in files], src_names,
                         library='np', step_size=step)
    out = tree.iterate(names, library='np', step_size=step,
                       entry_start=start, entry_stop=stop)

    buf_s = {k: np.empty(0, dtype=None) for k in src_names}
    n_s = n_o = 0
    pending = None
    for chunk_o in out:
        need = len(chunk_o['eventId'])
        n_o += need
        # pull enough source entries to cover this merged chunk
        while pending is None or len(pending['eventId']) < need:
            try:
                nxt = next(src)
            except StopIteration:
                break
            n_s += len(nxt['eventId'])
            pending = nxt if pending is None else {
                k: np.concatenate([pending[k], nxt[k]]) for k in nxt}
        if pending is None or len(pending['eventId']) < need:
            return False, 'source ran out after {} entries'.format(n_s)
        for k in src_names:
            if not np.array_equal(pending[k][:need], chunk_o[k]):
                return False, 'branch {} differs'.format(k)
        pending = {k: v[need:] for k, v in pending.items()}

    left = len(pending['eventId']) if pending is not None else 0
    for rest in src:
        left += len(rest['eventId'])
    if left:
        return False, 'source has {} extra entries'.format(left)
    return True, '{:,} entries'.format(n_o)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('merged')
    ap.add_argument('--runs-dir', default=RUNS_DIR)
    ap.add_argument('--step', type=int, default=5_000_000)
    ap.add_argument('--full', action='store_true',
                    help='compare every sub-run against its sources, not just '
                         'the first and last')
    args = ap.parse_args()

    r = uproot.open(args.merged)
    tree = r['hits']
    fail = 0

    # ---- 1. real TTree?
    kind = type(tree).__name__
    is_ttree = kind.startswith('Model_TTree')
    print('object type : {}  {}'.format(kind, 'ok' if is_ttree else
                                        'FAIL - not a TTree'))
    fail += 0 if is_ttree else 1

    names = [b.name for b in tree.branches]
    print('branches    : {}'.format(', '.join(names)))
    print('entries     : {:,}'.format(tree.num_entries))
    print('size on disk: {:.2f} GB'.format(os.path.getsize(args.merged) / 1e9))

    sr = r['subruns'].arrays(library='np')
    meta = json.loads(r['meta'])
    n_sub = len(sr['subrun_id'])
    print('sub-runs    : {}'.format(n_sub))

    # ---- 2. bookkeeping consistency
    if int(sr['n_hits'].sum()) != tree.num_entries:
        print('FAIL: subruns n_hits sums to {:,}, tree has {:,}'.format(
            int(sr['n_hits'].sum()), tree.num_entries))
        fail += 1
    if meta['n_hits'] != tree.num_entries:
        print('FAIL: meta n_hits {:,} != tree {:,}'.format(
            meta['n_hits'], tree.num_entries))
        fail += 1

    # ---- 3. streaming scan: geventId identity + contiguous sub-run blocks
    print('\nscanning {:,} entries...'.format(tree.num_entries))
    counts = np.zeros(n_sub, dtype=np.int64)
    ev_lo = np.full(n_sub, np.iinfo(np.int64).max, dtype=np.int64)
    ev_hi = np.full(n_sub, -1, dtype=np.int64)
    prev_max = -1
    contiguous = True
    for chunk in tree.iterate(['subrun_id', 'eventId', 'geventId'],
                              library='np', step_size=args.step):
        sid, ev, g = chunk['subrun_id'], chunk['eventId'], chunk['geventId']
        if not np.array_equal(g // STRIDE, sid.astype(np.uint64)):
            print('FAIL: geventId // stride != subrun_id'); fail += 1
        if not np.array_equal(g % STRIDE, ev):
            print('FAIL: geventId % stride != eventId'); fail += 1
        if sid.min() < prev_max:
            contiguous = False
        prev_max = int(sid.max())
        u, c = np.unique(sid, return_counts=True)
        counts[u] += c
        for s in u:
            m = sid == s
            ev_lo[s] = min(ev_lo[s], int(ev[m].min()))
            ev_hi[s] = max(ev_hi[s], int(ev[m].max()))
    print('  geventId decodes to (subrun_id, eventId): ok')
    print('  sub-runs stored contiguously            : {}'.format(contiguous))
    if not contiguous:
        fail += 1

    # ---- 4. per sub-run counts vs the lookup tree
    print('\nper sub-run (counts vs the subruns tree):')
    bad = 0
    for i in range(n_sub):
        ok = (counts[i] == int(sr['n_hits'][i])
              and ev_lo[i] == int(sr['event_id_min'][i])
              and ev_hi[i] == int(sr['event_id_max'][i]))
        bad += 0 if ok else 1
        print('  {:3d} {:<24}/{:<24} {:>12,}  ev {}..{}  {}'.format(
            int(sr['subrun_id'][i]), sr['run'][i], sr['sub_run'][i],
            int(counts[i]), ev_lo[i], ev_hi[i], 'ok' if ok else 'FAIL'))
    fail += bad

    # ---- 5. entry-by-entry against the sources
    offsets = np.concatenate([[0], np.cumsum(sr['n_hits'])])
    which = range(n_sub) if args.full else sorted({0, n_sub - 1})
    print('\nentry-by-entry against sources ({}):'.format(
        'all sub-runs' if args.full else 'first and last only, use --full for all'))
    for i in which:
        sub_dir = os.path.join(args.runs_dir, sr['run'][i], sr['sub_run'][i])
        ok, note = compare_source(tree, names, sub_dir,
                                  int(offsets[i]), int(offsets[i + 1]),
                                  args.step)
        if ok is False:
            fail += 1
        print('  {:<24}/{:<24} identical={}  ({})'.format(
            sr['run'][i], sr['sub_run'][i], ok, note))

    print('\n' + ('ALL CHECKS PASSED' if fail == 0
                  else '{} CHECK(S) FAILED'.format(fail)))
    return 1 if fail else 0


if __name__ == '__main__':
    sys.exit(main())
