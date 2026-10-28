#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Channel order per connector from straight tracks: which wiring option makes det2 sit on the line?

For one connector at a time, every option of orientation_check.OPTIONS (channel order of the bottom and top VMM,
swapped halves, interleaved) is tried while the other connectors keep the SPS order. Score = number of groups whose
best triplet has det2 within --tol mm of the line (after the median det2 offset), minus the same count with the det2
hits taken from another group (chance, --n-shuffle times). Groups are split in two halves (even / odd group index):
the option is chosen on one half and its score is reported on the other, so a free choice among ~40 options cannot
fit noise.

Usage: half_order_scan.py cand.npz[,more.npz] OUT.json [--conn-shift -1] [--tol 20] [--max-hits 12] [--connectors P2_MID:3,...]
       [--veto 2:27,6:44]
"""
import argparse
import json
import sys

import numpy as np

import orientation_check as O
from straight_track_rate import groups, parse_veto
from event_display import best_triplets, shift_of, best_per_group


def score(h, tol, rng, n_shuffle):
    t, evs = best_triplets(h)
    s = shift_of(t)
    n = int((best_per_group(t, s)['r'] < tol).sum())
    ch = [int((best_per_group(best_triplets(h, rng.permutation(len(evs)))[0], s)['r'] < tol).sum())
          for _ in range(n_shuffle)]
    return n - float(np.mean(ch)), n, float(np.mean(ch))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('cand')
    ap.add_argument('out_json')
    ap.add_argument('--conn-shift', type=int, default=-1)
    ap.add_argument('--tol', type=float, default=20)
    ap.add_argument('--max-hits', type=int, default=12)
    ap.add_argument('--n-shuffle', type=int, default=5)
    ap.add_argument('--connectors', default=None, help='e.g. P2_MID:3,P2_OUT:4 (default: all)')
    ap.add_argument('--veto', default='', help='vmm:ch,... whose groups are dropped, e.g. 2:27,6:44')
    a = ap.parse_args()
    veto = parse_veto(a.veto)
    rng = np.random.default_rng(11)
    mp = O.Mapper(a.conn_shift)
    keys = sorted(mp.conn)
    if a.connectors:
        want = {(c.split(':')[0], int(c.split(':')[1])) for c in a.connectors.split(',')}
        keys = [k for k in keys if k in want]

    X, Y = mp.xy_table()
    h, _ = groups(a.cand.split(','), a.max_hits, X, Y, veto)
    evs = np.unique(h['ev'].to_numpy())
    half = {e: i % 2 for i, e in enumerate(evs)}
    base = {p: score(h[h['ev'].map(half) == p], a.tol, rng, a.n_shuffle) for p in (0, 1)}
    print(f'SPS order everywhere: signal half A {base[0][0]:.0f}, half B {base[1][0]:.0f}')
    out = {'params': vars(a), 'base': base, 'connectors': {}}
    for key in keys:
        res = {}
        for name in O.OPTIONS:
            Xc, Yc = mp.xy_table({key: name})
            hc, _ = groups(a.cand.split(','), a.max_hits, Xc, Yc, veto)
            hp = hc['ev'].map(half)
            res[name] = [score(hc[hp == p], a.tol, rng, a.n_shuffle)[0] for p in (0, 1)]
        best_a = max(res, key=lambda k: res[k][0])
        best_b = max(res, key=lambda k: res[k][1])
        top = sorted(res, key=lambda k: -(res[k][0] + res[k][1]))[:4]
        out['connectors'][f'{key[0]}:c{key[1]}'] = {'scores': res, 'best_on_A': best_a, 'best_on_B': best_b}
        print(f'{key[0]} c{key[1]} (VMM {mp.conn[key]["bot"]}/{mp.conn[key]["top"]}): '
              f'SPS {res[O.SPS][0]:.0f}/{res[O.SPS][1]:.0f} | chosen on A: {best_a} -> held-out B {res[best_a][1]:.0f} '
              f'| chosen on B: {best_b} -> held-out A {res[best_b][0]:.0f} | top: '
              + ', '.join(f'{k} {res[k][0]:.0f}/{res[k][1]:.0f}' for k in top), flush=True)
    json.dump(out, open(a.out_json, 'w'), indent=1)


if __name__ == '__main__':
    sys.exit(main())
