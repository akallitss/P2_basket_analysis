#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Per-trigger dump of the scintillator-referenced lab runs (the event-level companion of scint_efficiency.py).

Same trigger (offline (top OR) AND bottom, --bad-times veto), same latency window and clean-cluster definition as
scint_efficiency.py, but instead of the tag-and-probe statistics it writes every trigger with
  * the paddle pulses on hybrid 0 (--paddles 0:40,0:41,0:44,0:45): ADC and time relative to the trigger of the
    pulse nearest the trigger within +- --coinc-win ns (NaN if none);
  * per chamber: number of live hits in the latency window, the clean flag, the cluster centroid;
  * the live hits of the window: trigger index, chamber, vmm, channel, adc, t - t_trig, pad x, y.
Output OUT.npz (arrays) + OUT.json (run summary, map options). Analyse with handoff_scratch / trigger_dump_*.py.

Usage: scint_dump.py f1.pcapng[,...] OUT.npz --top 0:41,0:45 --bottom 0:44 --conn-shift -1 --order ... [--bad-times veto.json]
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


def nearest(ts, adcs, t0, win):
    """For each t0: adc and dt of the pulse in ts nearest to t0 within +- win (NaN if none)."""
    i = np.searchsorted(ts, t0)
    best_dt = np.full(len(t0), np.nan)
    best_adc = np.full(len(t0), np.nan)
    for k in (i - 1, i):
        k = np.clip(k, 0, max(len(ts) - 1, 0))
        if len(ts) == 0:
            break
        dt = ts[k] - t0
        ok = np.abs(dt) <= win
        better = ok & (np.isnan(best_dt) | (np.abs(dt) < np.abs(best_dt)))
        best_dt[better] = dt[better]
        best_adc[better] = adcs[k][better]
    return best_adc, best_dt


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('files')
    ap.add_argument('out')
    ap.add_argument('--top', default='0:41,0:45')
    ap.add_argument('--bottom', default='0:44')
    ap.add_argument('--paddles', default='0:40,0:41,0:44,0:45')
    ap.add_argument('--coinc-win', type=float, default=60)
    ap.add_argument('--conn-shift', type=int, default=-1)
    ap.add_argument('--order', action='append')
    ap.add_argument('--veto', default='')
    ap.add_argument('--dead-vmm', default='')
    ap.add_argument('--half-win', type=float, default=150)
    ap.add_argument('--max-tag', type=int, default=6)
    ap.add_argument('--cluster-r', type=float, default=15)
    ap.add_argument('--title', default='')
    ap.add_argument('--bad-times', default='')
    a = ap.parse_args()

    X, Y = O.Mapper(a.conn_shift).xy_table(O.parse_orders(a.order))
    state = {}
    V, T, C, A = [], [], [], []
    clock = WallClock()
    for p in a.files.split(','):
        recs = read_times(p, state, with_adc=True)
        ts = [t for _, t, _, _ in recs if len(t)]
        clock.add_file(p, min(t.min() for t in ts) if ts else None)
        for v, t, c, ad in recs:
            V.append(v); T.append(t); C.append(c); A.append(ad)
    V, T, C, A = (np.concatenate(x) for x in (V, T, C, A))
    o = np.argsort(T, kind='stable')
    V, T, C, A = V[o].astype(int), T[o], C[o].astype(int), A[o].astype(int)
    span = (T[-1] - T[0]) * 1e-9

    def chan(spec):
        v, c = (int(s) for s in spec.split(':'))
        m = (V == v) & (C == c)
        return T[m], A[m]

    def times(spec):
        return np.sort(np.concatenate([chan(s)[0] for s in spec.split(',')]))

    top, bot = times(a.top), times(a.bottom)
    s = np.searchsorted(bot, top - a.coinc_win)
    e = np.searchsorted(bot, top + a.coinc_win)
    trig = top[e > s]
    trig = trig[np.concatenate([[True], np.diff(trig) > 300])] if len(trig) else trig
    n_before = int(len(trig))
    veto_info = None
    if a.bad_times:
        iv = load_veto(a.bad_times)
        bad = in_intervals(clock(trig), iv)
        veto_info = {'file': a.bad_times, 'n_trig_before': n_before, 'n_trig_vetoed': int(bad.sum())}
        trig = trig[~bad]
    counts = np.zeros((32, 64), int)
    np.add.at(counts, (V, C), 1)
    live = np.isfinite(X) & (counts > 0)
    for v in (int(x) for x in a.dead_vmm.split(',') if x):
        live[v] = False
    for v, c in parse_veto(a.veto):
        live[v, c] = False

    out = {'t_trig': trig, 'wall': clock(trig)}
    # paddle pulses nearest the trigger
    for spec in a.paddles.split(','):
        ts, ad = chan(spec)
        adc, dt = nearest(ts, ad, trig, a.coinc_win)
        key = 'p' + spec.replace(':', '_')
        out[key + '_adc'], out[key + '_dt'] = adc, dt
    # chamber hits in the latency window
    hits = {k: [] for k in ('itrig', 'det', 'vmm', 'ch', 'adc', 'dt', 'x', 'y')}
    summary = {'title': a.title, 'span_s': span, 'n_trig': int(len(trig)), 'n_trig_before_veto': n_before,
               'bad_times': veto_info, 'params': vars(a), 'chambers': {}}
    for d in (3, 2, 1):
        m = (DET[V] == d) & live[V, C]
        ht, hv, hc, ha = T[m], V[m], C[m], A[m]
        hx, hy = X[hv, hc], Y[hv, hc]
        s = np.searchsorted(ht, trig)
        e = np.searchsorted(ht, trig + 1000)
        first = np.array([ht[x] - t for x, y, t in zip(s, e, trig) if y > x])
        lat = float(np.median(first)) if len(first) else 0.0
        lo, hi = lat - a.half_win, lat + a.half_win
        s = np.searchsorted(ht, trig + lo)
        e = np.searchsorted(ht, trig + hi, side='right')
        n = e - s
        cx, cy, clean = np.full(len(trig), np.nan), np.full(len(trig), np.nan), np.zeros(len(trig), bool)
        for k in np.nonzero(n > 0)[0]:
            x, y = hx[s[k]:e[k]], hy[s[k]:e[k]]
            mx, my = x.mean(), y.mean()
            cx[k], cy[k] = mx, my
            if len(x) <= a.max_tag and np.max(np.hypot(x - mx, y - my)) <= a.cluster_r:
                clean[k] = True
            hits['itrig'].append(np.full(n[k], k))
            hits['det'].append(np.full(n[k], d))
            hits['vmm'].append(hv[s[k]:e[k]])
            hits['ch'].append(hc[s[k]:e[k]])
            hits['adc'].append(ha[s[k]:e[k]])
            hits['dt'].append(ht[s[k]:e[k]] - trig[k])
            hits['x'].append(x)
            hits['y'].append(y)
        out[f'd{d}_n'], out[f'd{d}_clean'], out[f'd{d}_cx'], out[f'd{d}_cy'] = n, clean, cx, cy
        summary['chambers'][d] = {'latency_ns': lat, 'any_hit': float((n > 0).mean()), 'clean': float(clean.mean())}
        print(f'det{d}: latency {lat:.0f} ns, hit {np.mean(n > 0):.3f}, clean {clean.mean():.3f}')
    for k, v in hits.items():
        out['h_' + k] = np.concatenate(v) if v else np.array([])
    out['live'] = live
    out['padX'], out['padY'] = X, Y
    np.savez_compressed(a.out, **out)
    json.dump(summary, open(a.out.replace('.npz', '') + '.json', 'w'), indent=1, default=float)
    print(f'{a.title}: {len(trig)} triggers written to {a.out}')


if __name__ == '__main__':
    sys.exit(main())
