#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Is the fixed 140-310 ns 'any hit' window of the muon funnel the same for the three chambers? (20 Nov 2026)

For every scintillator coincidence where BOTH tag chambers of a probe have a clean cluster 140-310 ns after the
scintillator (the funnel's tag condition), histogram, in 10 ns bins from -1000 to +3000 ns relative to the scintillator:
  * t_first: the earliest hit of the probe chamber (any channel) in that range;
  * t_all: every probe hit (the time spectrum).
The same is filled --offset ns later (random hits). Prints, per probe, the fraction of tagged muons with a hit in the
funnel window, in a wide window (-200..+1500), the median / 1st / 99th percentile of t_first, and how many of the
window-silent muons have their first hit later in the wide window (on minus off).

Usage: time_window_check.py EFF_JSON [EFF_JSON ...] --out OUT.json   (efficiency_tp.py jsons, as muon_funnel.py)
"""
import argparse
import json
import sys

import numpy as np

import orientation_check as O
from efficiency_tp import load, split_groups
from muon_funnel import Chamber, NAME, TAGS
from straight_track_rate import parse_veto

E = np.arange(-1000, 3010, 10)


def one_run(path, a, acc):
    J = json.load(open(path))
    p = J['params']
    X, Y = O.Mapper(p['conn_shift']).xy_table(O.parse_orders(p['order']))
    V, T, C, A = load(p['files'])
    counts = np.zeros((32, 64), int)
    np.add.at(counts, (V, C), 1)
    live = np.isfinite(X) & (counts >= p['min_hits'])
    for v in {int(v) for v in str(p['dead_vmm']).split(',') if v}:
        live[v] = False
    vmask = np.zeros((32, 64), bool)
    for v, c in parse_veto(p['veto']):
        live[v, c] = False
        vmask[v, c] = True
    g = split_groups(T, p['gap'])
    burst_hit = np.bincount(g)[g] > p['burst']
    CH = {d: Chamber(d, V, T, C, A, X, Y, live, vmask[V, C], burst_hit, p) for d in (3, 2, 1)}
    sv, sc = map(int, p['scint_ch'].split(':'))
    ts_all = np.sort(T[(V == sv) & (C == sc)])
    ts_all = ts_all[(ts_all > T[0] + 1e4) & (ts_all < T[-1] - 1e4)]
    S = {d: [CH[d].state(ts, a.lo, a.hi)[0] for ts in ts_all] for d in (3, 2, 1)}
    for pr, (ta, tb) in TAGS.items():
        r = acc.setdefault(NAME[pr], {w: {'first': np.zeros(len(E) - 1), 'all': np.zeros(len(E) - 1), 'n': 0}
                                      for w in ('on', 'off')})
        t = CH[pr].t
        for i, ts in enumerate(ts_all):
            if S[ta][i] != 'clean' or S[tb][i] != 'clean':
                continue
            for w, sh in (('on', 0.0), ('off', a.offset)):
                s0, e0 = np.searchsorted(t, ts + sh + E[0]), np.searchsorted(t, ts + sh + E[-1])
                dt = t[s0:e0] - ts - sh
                r[w]['n'] += 1
                if len(dt):
                    r[w]['all'] += np.histogram(dt, E)[0]
                    r[w]['first'][min(np.searchsorted(E, dt.min(), side='right') - 1, len(E) - 2)] += 1


def frac(h, lo, hi):
    return h[(E[:-1] >= lo) & (E[:-1] < hi)].sum()


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('jsons', nargs='+')
    ap.add_argument('--out', required=True)
    ap.add_argument('--lo', type=float, default=140)
    ap.add_argument('--hi', type=float, default=310)
    ap.add_argument('--offset', type=float, default=20000)
    a = ap.parse_args()
    acc = {}
    for j in a.jsons:
        one_run(j, a, acc)
        print(j, flush=True)
    for k, r in acc.items():
        n = r['on']['n']
        f = r['on']['first'] - r['off']['first'] * n / max(r['off']['n'], 1)
        h = r['on']['all'] - r['off']['all'] * n / max(r['off']['n'], 1)
        c = np.cumsum(h) / max(h.sum(), 1)
        q = {x: float(E[np.searchsorted(c, x)]) for x in (0.01, 0.05, 0.5, 0.95, 0.99)}
        print(f'{k}: tagged {n}  hit in {a.lo:.0f}-{a.hi:.0f} ns {frac(r["on"]["first"], a.lo, a.hi) / n:.3f}  '
              f'first hit in -200..1500 {frac(r["on"]["first"], -200, 1500) / n:.3f} (off {frac(r["off"]["first"], -200, 1500) / max(r["off"]["n"], 1):.3f})'
              f'  all-hit time quantiles (on-off) {q}')
        print(f'   first hit before {a.lo:.0f}: {frac(f, -1000, a.lo) / n:.4f}  after {a.hi:.0f} (to 1500): '
              f'{frac(f, a.hi, 1500) / n:.4f}  (random-subtracted, per tagged muon)')
        for w in ('on', 'off'):
            r[w] = {kk: (vv.tolist() if hasattr(vv, 'tolist') else vv) for kk, vv in r[w].items()}
    json.dump({'edges': E.tolist(), 'lo': a.lo, 'hi': a.hi, 'probes': acc}, open(a.out, 'w'))


if __name__ == '__main__':
    sys.exit(main())
