#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Paddle efficiency with the CHAMBERS as the reference: three-chamber straight tracks found in the raw hit stream
(no paddle in the selection), then for each paddle channel on hybrid 0 whether a pulse sits at the chamber time
minus the chamber latency, within +- --win ns.

Track: a det2 live hit with det3 and det1 live hits within +- --half-win ns, at most --max-tag hits per chamber all
within --cluster-r mm of their mean (the clean cluster of scint_efficiency.py), det2 within --tol mm of the midpoint
of det3 and det1 (z = 320 / 630 / 940). One track per 300 ns. Chance level: the same with the paddle pulses looked
up --sideband ns earlier.

Output OUT.npz: per track the three cluster positions, the det2 time, and per paddle the pulse time offset (NaN if
none) in the signal and in the sideband window; OUT.json summary.

Usage: paddle_efficiency.py f1.pcapng[,...] OUT.npz --conn-shift -1 --order ... [--paddles 0:40,0:41,0:44,0:45]
       [--latency 242] [--win 120] [--bad-times veto.json]
"""
import argparse
import json
import sys

import numpy as np

import orientation_check as O
from coincidence_check import read_times
from straight_track_rate import parse_veto
from spark_veto import WallClock, load_veto, in_intervals

DET = np.zeros(32, int)
DET[2:8], DET[8:14], DET[14:20] = 3, 2, 1
Z = {3: 320.0, 2: 630.0, 1: 940.0}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('files')
    ap.add_argument('out')
    ap.add_argument('--paddles', default='0:40,0:41,0:44,0:45')
    ap.add_argument('--latency', type=float, default=242.0, help='chamber time - paddle time [ns]')
    ap.add_argument('--win', type=float, default=120.0)
    ap.add_argument('--sideband', type=float, default=5000.0)
    ap.add_argument('--conn-shift', type=int, default=-1)
    ap.add_argument('--order', action='append')
    ap.add_argument('--veto', default='')
    ap.add_argument('--dead-vmm', default='')
    ap.add_argument('--half-win', type=float, default=150)
    ap.add_argument('--max-tag', type=int, default=6)
    ap.add_argument('--cluster-r', type=float, default=15)
    ap.add_argument('--tol', type=float, default=25)
    ap.add_argument('--title', default='')
    ap.add_argument('--bad-times', default='')
    a = ap.parse_args()

    X, Y = O.Mapper(a.conn_shift).xy_table(O.parse_orders(a.order))
    state = {}
    V, T, C = [], [], []
    clock = WallClock()
    for p in a.files.split(','):
        recs = read_times(p, state)
        ts = [t for _, t, _ in recs if len(t)]
        clock.add_file(p, min(t.min() for t in ts) if ts else None)
        for v, t, c in recs:
            V.append(v); T.append(t); C.append(c)
    V, T, C = (np.concatenate(x) for x in (V, T, C))
    o = np.argsort(T, kind='stable')
    V, T, C = V[o].astype(int), T[o], C[o].astype(int)
    span = (T[-1] - T[0]) * 1e-9
    counts = np.zeros((32, 64), int)
    np.add.at(counts, (V, C), 1)
    live = np.isfinite(X) & (counts > 0)
    for v in (int(x) for x in a.dead_vmm.split(',') if x):
        live[v] = False
    for v, c in parse_veto(a.veto):
        live[v, c] = False

    H = {}
    for d in (3, 2, 1):
        m = (DET[V] == d) & live[V, C]
        H[d] = (T[m], X[V[m], C[m]], Y[V[m], C[m]])
        print(f'det{d}: {m.sum()} live hits ({m.sum() / span:.1f} Hz)')
    t2 = H[2][0]
    # seeds: det2 hits with >= 1 det3 and >= 1 det1 hit within +- half_win
    def near(d):
        t = H[d][0]
        s = np.searchsorted(t, t2 - a.half_win)
        e = np.searchsorted(t, t2 + a.half_win, side='right')
        return s, e
    s3, e3 = near(3)
    s1, e1 = near(1)
    s2, e2 = near(2)
    cand = np.nonzero((e3 > s3) & (e1 > s1))[0]
    # one track per 300 ns
    keep = np.concatenate([[True], np.diff(t2[cand]) > 300]) if len(cand) else np.zeros(0, bool)
    cand = cand[keep]
    print(f'{len(cand)} three-chamber time groups')
    rows = []
    for k in cand:
        cl = {}
        ok = True
        for d, s, e in ((3, s3, e3), (2, s2, e2), (1, s1, e1)):
            x, y = H[d][1][s[k]:e[k]], H[d][2][s[k]:e[k]]
            if len(x) > a.max_tag:
                ok = False
                break
            mx, my = x.mean(), y.mean()
            if np.max(np.hypot(x - mx, y - my)) > a.cluster_r:
                ok = False
                break
            cl[d] = (mx, my, len(x))
        if not ok:
            continue
        r = np.hypot(cl[2][0] - (cl[3][0] + cl[1][0]) / 2, cl[2][1] - (cl[3][1] + cl[1][1]) / 2)
        rows.append((t2[k], cl[3][0], cl[3][1], cl[2][0], cl[2][1], cl[1][0], cl[1][1], cl[3][2], cl[2][2], cl[1][2], r))
    R = np.array(rows)
    print(f'{len(R)} clean triple-cluster groups, {int((R[:, -1] < a.tol).sum())} straight (det2 within {a.tol:g} mm)')
    out = {'t2': R[:, 0], 'x3': R[:, 1], 'y3': R[:, 2], 'x2': R[:, 3], 'y2': R[:, 4], 'x1': R[:, 5], 'y1': R[:, 6],
           'n3': R[:, 7], 'n2': R[:, 8], 'n1': R[:, 9], 'r2': R[:, 10], 'wall': clock(R[:, 0])}
    if a.bad_times:
        out['bad'] = in_intervals(out['wall'], load_veto(a.bad_times))
    else:
        out['bad'] = np.zeros(len(R), bool)
    t0 = R[:, 0] - a.latency
    summ = {'title': a.title, 'span_s': span, 'n_groups': int(len(cand)), 'n_clean': int(len(R)),
            'n_straight': int((R[:, -1] < a.tol).sum()), 'params': vars(a), 'paddles': {}}
    for spec in a.paddles.split(','):
        v, c = (int(s) for s in spec.split(':'))
        tp = T[(V == v) & (C == c)]
        key = 'p' + spec.replace(':', '_')
        for lab, tt in (('', t0), ('_sb', t0 - a.sideband)):
            i = np.clip(np.searchsorted(tp, tt), 1, max(len(tp) - 1, 1))
            dts = np.c_[tp[i - 1] - tt, tp[np.minimum(i, len(tp) - 1)] - tt] if len(tp) else np.full((len(tt), 2), np.nan)
            j = np.argmin(np.abs(dts), 1)
            dt = dts[np.arange(len(tt)), j]
            dt[np.abs(dt) > a.win] = np.nan
            out[key + lab] = dt
        st = (R[:, -1] < a.tol) & ~out['bad']
        e = float(np.isfinite(out[key][st]).mean())
        sb = float(np.isfinite(out[key + '_sb'][st]).mean())
        summ['paddles'][spec] = {'eff_straight': e, 'chance': sb, 'singles_hz': len(tp) / span}
        print(f'paddle {spec}: {len(tp) / span:.1f} Hz singles; on straight tracks: {e:.3f} (chance {sb:.4f})')
    np.savez_compressed(a.out, **out)
    json.dump(summ, open(a.out.replace('.npz', '') + '.json', 'w'), indent=1, default=float)


if __name__ == '__main__':
    sys.exit(main())
