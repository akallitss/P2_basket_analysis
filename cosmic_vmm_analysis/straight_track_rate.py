#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Rate of straight three-chamber tracks above the chance level, per run -- the muon signal of the stand.

Input: .npz of dump_candidates.py. Groups that touch all three chambers with at most --max-hits hits; in each,
every triplet (one hit of det3, det2, det1 on mapped pads) is tested: STRAIGHT when the det2 pad is within --tol mm
of the midpoint of the det3 and det1 pads (the line through the outer hits, z = 320 / 630 / 940 mm, crosses det2
half way), after removing the median det2 offset of the run. A group counts once if any triplet is straight.
Chance level: the same count with the det2 hits taken from another group (groups shuffled, --n-shuffle times).
signal = straight - chance, per hour; it is proportional to the muon flux times the three efficiencies (and the
acceptance), so across threshold settings its maximum marks the best setting. The pad map is the SPS one
(orientation_check.py found no wiring that does better on held-out events).

A run spec may end with @f when only a fraction f of its capture files was decoded (the span from the first to
the last hit is then multiplied by f).

From run_112 on the chambers are read on connectors c3-c5 instead of c4-c6: --conn-shift -1.

Usage: straight_track_rate.py OUT.json run_76:cand76.npz[,more.npz] run_95:cand95.npz@0.25 ... [--tol 30] [--max-hits 12]
       [--conn-shift -1]
"""
import argparse
import json
import sys

import numpy as np
import pandas as pd

import orientation_check as O


def parse_veto(s):
    """'2:27,6:44' -> {(2, 27), (6, 44)}"""
    return {tuple(int(x) for x in p.split(':')) for p in s.split(',')} if s else set()


def groups(files, max_hits, X, Y, veto=()):
    """veto: (vmm, ch) whose groups are dropped whole (a discharging channel and the pickup it induces elsewhere)."""
    out = []
    span = 0.0
    off = 0
    for f in files:
        z = np.load(f)
        span += float(z['span_s'])
        h = pd.DataFrame({k: z[k] for k in ('vmm', 'ch', 'adc', 't', 'ev', 'det')})
        h['ev'] += off
        off = int(h['ev'].max()) + 1
        h['vmm'] = h['vmm'].astype(int)
        h['ch'] = h['ch'].astype(int)
        g = h.groupby('ev')
        h = h[(g['det'].transform('nunique') == 3) & (g['vmm'].transform('size') <= max_hits)]
        if veto:
            bad = pd.Series([(v, c) in veto for v, c in zip(h['vmm'], h['ch'])], index=h.index)
            h = h[~bad.groupby(h['ev']).transform('any')]
        h = h.assign(x=X[h['vmm'], h['ch']], y=Y[h['vmm'], h['ch']]).dropna(subset=['x'])
        out.append(h)
    return pd.concat(out, ignore_index=True), span


def count_straight(h, tol, perm=None):
    """number of groups with a straight triplet; perm maps group index -> group whose det2 hits are used."""
    evs = np.unique(h['ev'].to_numpy())
    idx = {e: i for i, e in enumerate(evs)}
    parts = {d: h[h['det'] == d] for d in (1, 2, 3)}
    gi = {d: parts[d]['ev'].map(idx).to_numpy() for d in (1, 2, 3)}
    if perm is not None:
        gi[2] = perm[gi[2]]
    # outer pairs per group
    a = pd.DataFrame({'g': gi[3], 'x3': parts[3]['x'].to_numpy(), 'y3': parts[3]['y'].to_numpy()})
    b = pd.DataFrame({'g': gi[1], 'x1': parts[1]['x'].to_numpy(), 'y1': parts[1]['y'].to_numpy()})
    c = pd.DataFrame({'g': gi[2], 'x2': parts[2]['x'].to_numpy(), 'y2': parts[2]['y'].to_numpy()})
    t = a.merge(b, on='g').merge(c, on='g')
    rx = t['x2'] - 0.5 * (t['x3'] + t['x1'])
    ry = t['y2'] - 0.5 * (t['y3'] + t['y1'])
    return t['g'].to_numpy(), rx.to_numpy(), ry.to_numpy(), len(evs)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('out_json')
    ap.add_argument('runs', nargs='+')
    ap.add_argument('--tol', type=float, default=30)
    ap.add_argument('--max-hits', type=int, default=12)
    ap.add_argument('--n-shuffle', type=int, default=10)
    ap.add_argument('--conn-shift', type=int, default=0, help='-1 for runs >= 112 (read-out on c3-c5)')
    ap.add_argument('--order', action='append', help="channel order of one connector, e.g. 'P2_IN:3=bot:rev|top:lin'")
    ap.add_argument('--veto', default='', help='vmm:ch,... whose groups are dropped, e.g. 2:27,6:44')
    a = ap.parse_args()
    X, Y = O.Mapper(a.conn_shift).xy_table(O.parse_orders(a.order))
    veto = parse_veto(a.veto)
    rng = np.random.default_rng(5)
    res = {'params': vars(a), 'runs': {}}
    print(f'{"run":>8} {"span h":>7} {"groups/h":>9} {"straight/h":>11} {"chance/h":>9} {"signal/h":>9} {"+-":>5}')
    for spec in a.runs:
        name, files = spec.split(':', 1)
        files, _, frac = files.partition('@')
        h, span = groups(files.split(','), a.max_hits, X, Y, veto)
        span *= float(frac) if frac else 1.0
        hours = span / 3600
        g, rx, ry, ng = count_straight(h, a.tol)
        near = np.hypot(rx, ry) < 80
        sx, sy = (np.median(rx[near]), np.median(ry[near])) if near.sum() > 20 else (0.0, 0.0)
        ok = np.hypot(rx - sx, ry - sy) < a.tol
        n_str = len(np.unique(g[ok]))
        ch = []
        for _ in range(a.n_shuffle):
            gs, sxr, syr, _ = count_straight(h, a.tol, rng.permutation(ng))
            ch.append(len(np.unique(gs[np.hypot(sxr - sx, syr - sy) < a.tol])))
        chance = float(np.mean(ch))
        # stricter null: each det2 hit keeps its chip (and so any chip-level noise correlation between chambers)
        # but gets the channel of another det2 hit of the same chip
        cc = []
        is2 = (h['det'] == 2).to_numpy()
        for _ in range(a.n_shuffle):
            hc = h.copy()
            c2 = hc.loc[is2, 'ch'].to_numpy().copy()
            v2 = hc.loc[is2, 'vmm'].to_numpy()
            for v in np.unique(v2):
                m = v2 == v
                c2[m] = rng.permutation(c2[m])
            hc.loc[is2, 'ch'] = c2
            hc.loc[is2, 'x'] = X[v2, c2]
            hc.loc[is2, 'y'] = Y[v2, c2]
            gs, sxr, syr, _ = count_straight(hc, a.tol)
            cc.append(len(np.unique(gs[np.hypot(sxr - sx, syr - sy) < a.tol])))
        chance_chip = float(np.mean(cc))
        sig = n_str - chance
        err = float(np.sqrt(n_str + np.var(ch) + chance / a.n_shuffle))
        sig_c = n_str - chance_chip
        err_c = float(np.sqrt(n_str + np.var(cc) + chance_chip / a.n_shuffle))
        res['runs'][name] = {'span_s': span, 'groups': ng, 'straight': n_str, 'chance': chance, 'signal': sig,
                             'signal_err': err, 'signal_per_h': sig / hours, 'signal_per_h_err': err / hours,
                             'chance_chip': chance_chip, 'signal_chip': sig_c, 'signal_chip_err': err_c,
                             'signal_chip_per_h': sig_c / hours, 'signal_chip_per_h_err': err_c / hours,
                             'det2_shift_mm': [float(sx), float(sy)]}
        print(f'{name:>8} {hours:7.2f} {ng / hours:9.0f} {n_str / hours:11.1f} {chance / hours:9.1f} {sig / hours:9.1f} {err / hours:5.1f}'
              f'   chip-null {chance_chip / hours:7.1f} -> signal {sig_c / hours:7.1f} +- {err_c / hours:5.1f}', flush=True)
    json.dump(res, open(a.out_json, 'w'), indent=1)


if __name__ == '__main__':
    sys.exit(main())
