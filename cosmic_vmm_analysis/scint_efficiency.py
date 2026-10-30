#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Chamber efficiency with the external scintillator coincidence as the reference (lab runs with paddles on H0).

The coincidence of the paddles (--trig vmm:ch on hybrid 0) marks a particle crossing the stack, independently of the
chambers. Per trigger and chamber: the clean clusters (as in efficiency_tp.py: 1..--max-tag live pads within
--cluster-r mm) and any live hit in [t + lo, t + hi], the window around the chamber latency (measured, median of the
first hit after the trigger); random hits from the same window --sideband ns before the trigger.

  1. any hit: fraction of triggers with a hit in the chamber (corrected for random hits);
  2. scintillator-confirmed tag-and-probe: of the triggers where the two OTHER chambers have a clean cluster, the
     fraction where the probe has a hit anywhere / within --probe-r mm of the line through the two clusters (+ the
     median probe offset); random hits as in 1.

Map as in efficiency_tp.py (--conn-shift, --order, --veto, --dead-vmm). Output OUT.json.

Usage: scint_efficiency.py f1.pcapng[,...] OUT.json --trig 0:40 --conn-shift -1 --order ... [--title run_130]
"""
import argparse
import json
import sys

import numpy as np

import orientation_check as O
from coincidence_check import read_times
from straight_track_rate import parse_veto

Z = {3: 320.0, 2: 630.0, 1: 940.0}
NAME = {3: 'det3 P2_IN', 2: 'det2 P2_MID', 1: 'det1 P2_OUT'}
DET = np.zeros(32, int)
DET[2:8], DET[8:14], DET[14:20] = 3, 2, 1


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('files')
    ap.add_argument('out_json')
    ap.add_argument('--trig', default='0:40')
    ap.add_argument('--conn-shift', type=int, default=-1)
    ap.add_argument('--order', action='append')
    ap.add_argument('--veto', default='')
    ap.add_argument('--dead-vmm', default='')
    ap.add_argument('--half-win', type=float, default=150)
    ap.add_argument('--sideband', type=float, default=3000)
    ap.add_argument('--max-tag', type=int, default=6)
    ap.add_argument('--cluster-r', type=float, default=15)
    ap.add_argument('--probe-r', type=float, default=15)
    ap.add_argument('--title', default='')
    a = ap.parse_args()

    X, Y = O.Mapper(a.conn_shift).xy_table(O.parse_orders(a.order))
    state = {}
    V, T, C = [], [], []
    for p in a.files.split(','):
        for v, t, c in read_times(p, state):
            V.append(v); T.append(t); C.append(c)
    V, T, C = (np.concatenate(x) for x in (V, T, C))
    o = np.argsort(T, kind='stable')
    V, T, C = V[o].astype(int), T[o], C[o].astype(int)
    span = (T[-1] - T[0]) * 1e-9
    tv, tc = (int(x) for x in a.trig.split(':'))
    trig = T[(V == tv) & (C == tc)]
    counts = np.zeros((32, 64), int)
    np.add.at(counts, (V, C), 1)
    live = np.isfinite(X) & (counts > 0)
    for v in (int(x) for x in a.dead_vmm.split(',') if x):
        live[v] = False
    for v, c in parse_veto(a.veto):
        live[v, c] = False

    res = {'title': a.title, 'span_s': span, 'n_trig': int(len(trig)), 'trig_rate_hz': len(trig) / span,
           'params': vars(a), 'chambers': {}}
    print(f'{a.title}: {len(trig)} triggers on {a.trig} ({len(trig) / span:.3f} Hz) in {span / 60:.1f} min')
    ev = {}
    for d in (3, 2, 1):
        m = (DET[V] == d) & live[V, C]
        ht, hx, hy = T[m], X[V[m], C[m]], Y[V[m], C[m]]
        # latency: first live hit within 0..1000 ns after the trigger
        s = np.searchsorted(ht, trig)
        e = np.searchsorted(ht, trig + 1000)
        first = np.array([ht[x] - t for x, y, t in zip(s, e, trig) if y > x])
        lat = float(np.median(first)) if len(first) else 0.0
        lo, hi = lat - a.half_win, lat + a.half_win
        s = np.searchsorted(ht, trig + lo)
        e = np.searchsorted(ht, trig + hi, side='right')
        s2 = np.searchsorted(ht, trig + lo - a.sideband)
        e2 = np.searchsorted(ht, trig + hi - a.sideband, side='right')
        has, acc = e > s, float(np.mean(e2 > s2))
        cx, cy, clean = np.full(len(trig), np.nan), np.full(len(trig), np.nan), np.zeros(len(trig), bool)
        for k in np.nonzero(has)[0]:
            x, y = hx[s[k]:e[k]], hy[s[k]:e[k]]
            mx, my = x.mean(), y.mean()
            if len(x) <= a.max_tag and np.max(np.hypot(x - mx, y - my)) <= a.cluster_r:
                cx[k], cy[k], clean[k] = mx, my, True
        ev[d] = {'has': has, 'acc': acc, 's': s, 'e': e, 'ht': ht, 'hx': hx, 'hy': hy, 'cx': cx, 'cy': cy, 'clean': clean}
        eff = (has.mean() - acc) / (1 - acc)
        res['chambers'][NAME[d]] = {'latency_ns': lat, 'any_hit': float(eff), 'any_hit_err': float(np.sqrt(eff * (1 - eff) / len(trig))),
                                    'random': acc, 'clean_cluster_frac': float(clean.mean())}
        print(f'{NAME[d]}: latency {lat:.0f} ns; hit anywhere {eff:.3f} +- {np.sqrt(eff * (1 - eff) / len(trig)):.3f} '
              f'(random {acc:.4f}); one clean cluster {clean.mean():.3f}')
    allthree = ev[3]['has'] & ev[2]['has'] & ev[1]['has']
    res['all_three'] = float(allthree.mean())
    print(f'all three chambers hit: {allthree.mean():.3f}')
    for p, (ta, tb) in ((2, (3, 1)), (1, (3, 2)), (3, (2, 1))):
        A, B, P = ev[ta], ev[tb], ev[p]
        sel = A['clean'] & B['clean']
        f = (Z[p] - Z[ta]) / (Z[tb] - Z[ta])
        px = A['cx'] + f * (B['cx'] - A['cx'])
        py = A['cy'] + f * (B['cy'] - A['cy'])
        dmin = np.full(len(trig), np.inf)
        for k in np.nonzero(sel & P['has'])[0]:
            dmin[k] = np.min(np.hypot(P['hx'][P['s'][k]:P['e'][k]] - px[k], P['hy'][P['s'][k]:P['e'][k]] - py[k]))
        near = sel & (dmin < 60)
        ox = float(np.median((P['cx'] - px)[near & P['clean']])) if (near & P['clean']).sum() > 10 else 0.0
        oy = float(np.median((P['cy'] - py)[near & P['clean']])) if (near & P['clean']).sum() > 10 else 0.0
        for k in np.nonzero(sel & P['has'])[0]:
            dmin[k] = np.min(np.hypot(P['hx'][P['s'][k]:P['e'][k]] - px[k] - ox, P['hy'][P['s'][k]:P['e'][k]] - py[k] - oy))
        n = int(sel.sum())
        acc = P['acc']
        r = {'tags': [NAME[ta], NAME[tb]], 'n_tagged': n, 'offset_mm': [ox, oy]}
        for lab, ok in (('any_position', P['has']), (f'within_{a.probe_r:.0f}mm', dmin < a.probe_r), ('within_50mm', dmin < 50)):
            k = int((sel & ok).sum())
            e = (k / n - acc) / (1 - acc) if n else float('nan')
            r[lab] = {'k': k, 'eff': e, 'err': float(np.sqrt(max(e * (1 - e), 0) / n)) if n else None}
        res[f'probe_{NAME[p]}'] = r
        print(f'scintillator-confirmed tag-and-probe, probe {NAME[p]} (clean clusters in {NAME[ta]} and {NAME[tb]}: {n}): '
              + ', '.join(f'{lab} {v["eff"]:.3f} +- {v["err"]:.3f}' for lab, v in r.items() if isinstance(v, dict) and 'eff' in v))
    json.dump(res, open(a.out_json, 'w'), indent=1, default=float)


if __name__ == '__main__':
    sys.exit(main())
