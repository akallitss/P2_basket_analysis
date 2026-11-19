#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Follow every scintillator coincidence through the whole tag-and-probe chain, one class per step (20 Nov 2026).

The reference is the scintillator coincidence (--scint-ch, H0 VMM 0 ch 40): a muon went through the paddles. For each
probe chamber P (tags A, B = the other two, as in efficiency_tp.py: MID <- IN+OUT, OUT <- IN+MID, IN <- MID+OUT) every
coincidence ends in exactly one class:

  T0  no_tag_A / no_tag_B / no_tag_both   a tag chamber has no hit 140-310 ns after the scintillator (the muon missed
                                          its read-out area, or that chamber was inefficient)
  T1  tag_rejected                        both tag chambers fired, but at least one hit group is not a clean cluster
                                          (> --max-tag hits, spread > --cluster-r mm, dead/masked pad, stack burst)
  F   outside_area / dead_pad             tags clean, but the predicted point (straight line + the run's median probe
                                          offset) is > --fid-r mm from any probe pad, or on a dead/masked pad
  P0  probe_silent                        tags clean, prediction on a live pad, the probe has NO hit at all 140-310 ns
  P1  probe_far                           the probe fired, but no live pad within --probe-r mm of the prediction
  P2  probe_near_late                     a live pad within --probe-r mm fired 140-310 ns after the scintillator, but
                                          not within +-(--win) ns of tag time + latency (the tag-and-probe window)
  OK  found                               a live pad within --probe-r mm in the tag-and-probe window

Random hits: the probe classes are re-evaluated with the probe window shifted by --offset ns (same tags), and every
class count gets the shifted-window version subtracted in the summary.

Extras written for the explanation page: distance of the nearest probe hit (P1 and OK), the predicted chip (per-chip
table), for rejected tag groups the spatial pieces (single linkage within --cluster-r mm) with their separation, chip
and channel difference from the main (largest-ADC) piece, and up to --examples example events per class with all hits
of the three chambers within -200..+600 ns of the scintillator.

Usage: muon_funnel.py EFF_JSON [EFF_JSON ...] --out OUT.json   (efficiency_tp.py jsons: files + cuts from 'params')
"""
import argparse
import json
import sys
from collections import Counter, defaultdict

import numpy as np

import orientation_check as O
from efficiency_tp import DET, load, split_groups, Z
from straight_track_rate import parse_veto

NAME = {3: 'P2_IN', 2: 'P2_MID', 1: 'P2_OUT'}
TAGS = {2: (3, 1), 1: (3, 2), 3: (2, 1)}
CLASSES = ['no_tag_both', 'no_tag_A', 'no_tag_B', 'tag_rejected', 'outside_area', 'dead_pad', 'probe_silent',
           'probe_far', 'probe_near_late', 'found']


def pieces(x, y, r):
    """single-linkage labels of points within r mm."""
    n = len(x)
    lab = -np.ones(n, int)
    k = 0
    for i in range(n):
        if lab[i] >= 0:
            continue
        stack, lab[i] = [i], k
        while stack:
            j = stack.pop()
            nb = np.nonzero((lab < 0) & (np.hypot(x - x[j], y - y[j]) <= r))[0]
            lab[nb] = k
            stack += nb.tolist()
        k += 1
    return lab


class Chamber:
    """time groups of one chamber with the efficiency_tp.chamber_clusters cuts."""

    def __init__(self, d, V, T, C, A, X, Y, live, bad_hit, burst_hit, p):
        m = DET[V] == d
        self.t, self.v, self.c, self.a = T[m], V[m], C[m], A[m]
        self.x, self.y = X[self.v, self.c], Y[self.v, self.c]
        self.ok = live[self.v, self.c]
        g = split_groups(self.t, p['gap'])
        ng = int(g[-1]) + 1 if len(g) else 0
        n_all = np.bincount(g, minlength=ng)
        n_live = np.bincount(g, weights=self.ok, minlength=ng)
        sx = np.bincount(g, weights=np.where(self.ok, self.x, 0), minlength=ng)
        sy = np.bincount(g, weights=np.where(self.ok, self.y, 0), minlength=ng)
        with np.errstate(invalid='ignore', divide='ignore'):
            self.mx, self.my = sx / n_live, sy / n_live
        r = np.hypot(self.x - self.mx[g], self.y - self.my[g])
        spread = np.zeros(ng)
        np.maximum.at(spread, g, np.where(self.ok, r, 0))
        bad = np.bincount(g, weights=bad_hit[m] | burst_hit[m], minlength=ng) > 0
        self.t0 = np.full(ng, np.inf); np.minimum.at(self.t0, g, self.t)
        self.t1 = np.full(ng, -np.inf); np.maximum.at(self.t1, g, self.t)
        self.n_all = n_all
        self.reason = {'burst': bad, 'dead': n_live < n_all, 'npads': n_all > p['max_tag'],
                       'wide': spread > p['cluster_r']}
        self.clean = (n_live >= 1) & (n_all <= p['max_tag']) & (n_live == n_all) & (spread <= p['cluster_r']) & ~bad
        self.starts = np.searchsorted(g, np.arange(ng))
        self.ends = np.r_[self.starts[1:], len(g)]
        self.lx, self.ly, self.lt = self.x[self.ok], self.y[self.ok], self.t[self.ok]   # live hits (time-ordered)

    def state(self, ts, lo, hi):
        """('none'|'clean'|'rejected', group index)."""
        a = np.searchsorted(self.t1, ts + lo)
        b = np.searchsorted(self.t0, ts + hi, side='right')
        if b <= a:
            return 'none', -1
        gs = np.arange(a, b)
        cl = gs[self.clean[gs]]
        if len(cl):
            return 'clean', int(cl[np.argmin(np.abs(self.t0[cl] - ts - 0.5 * (lo + hi)))])
        return 'rejected', int(gs[np.argmax(self.n_all[gs])])

    def any_hit(self, t_lo, t_hi):
        return np.searchsorted(self.t, t_hi, side='right') > np.searchsorted(self.t, t_lo)

    def nearest_live(self, px, py, t_lo, t_hi):
        s, e = np.searchsorted(self.lt, t_lo), np.searchsorted(self.lt, t_hi, side='right')
        if e <= s:
            return np.inf
        return float(np.min(np.hypot(self.lx[s:e] - px, self.ly[s:e] - py)))

    def hits(self, t_lo, t_hi, ts):
        s, e = np.searchsorted(self.t, t_lo), np.searchsorted(self.t, t_hi, side='right')
        return [[round(float(self.x[k]), 1) if np.isfinite(self.x[k]) else None,
                 round(float(self.y[k]), 1) if np.isfinite(self.y[k]) else None,
                 round(float(self.t[k] - ts)), int(self.a[k]), int(self.v[k]), int(self.c[k]), bool(self.ok[k])]
                for k in range(s, e)]


def piece_info(ch, gidx, r):
    """spatial pieces of a rejected group: per extra piece (distance to main, n hits, max ADC, same chip, |dch|, dt)."""
    s0, s1 = ch.starts[gidx], ch.ends[gidx]
    ok = ch.ok[s0:s1]
    x, y, a, v, c, t = (q[s0:s1][ok] for q in (ch.x, ch.y, ch.a, ch.v, ch.c, ch.t))
    if not len(x):
        return 0, []
    lab = pieces(x, y, r)
    k = lab.max() + 1
    adc = np.array([a[lab == i].sum() for i in range(k)])
    main = int(np.argmax(adc))
    mx, my = x[lab == main].mean(), y[lab == main].mean()
    out = []
    for i in range(k):
        if i == main:
            continue
        q = lab == i
        same = np.isin(v[q], v[lab == main])
        dch = min((abs(int(cc) - int(c0)) for vv, cc in zip(v[q], c[q]) for v0, c0 in zip(v[lab == main], c[lab == main])
                   if vv == v0), default=-1)
        out.append([round(float(np.hypot(x[q].mean() - mx, y[q].mean() - my)), 1), int(q.sum()), int(a[q].max()),
                    bool(same.any()), int(dch), round(float(t[q].min() - t[lab == main].min()))])
    return k, out


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
    # pads of every chamber (for the fiducial and the drawings)
    pads = {}
    for d in (3, 2, 1):
        vv, cc = np.nonzero(np.isfinite(X))
        k = DET[vv] == d
        pads[d] = (X[vv[k], cc[k]], Y[vv[k], cc[k]], vv[k], cc[k], live[vv[k], cc[k]])
        acc.setdefault('pads', {}).setdefault(NAME[d], [[round(float(x), 1), round(float(y), 1), int(v), int(c)]
                                                        for x, y, v, c in zip(*pads[d][:4])])
    sv, sc = map(int, p['scint_ch'].split(':'))
    ts_all = np.sort(T[(V == sv) & (C == sc)])
    ts_all = ts_all[(ts_all > T[0] + 1e4) & (ts_all < T[-1] - 1e4)]
    lo, hi = a.lo, a.hi
    S = {d: [CH[d].state(ts, lo, hi) for ts in ts_all] for d in (3, 2, 1)}
    for d in (3, 2, 1):
        r_ = acc['rejects'].setdefault(NAME[d], {'states': Counter(), 'reasons': Counter(), 'pieces': Counter(),
                                                 'extra': []})
        r_.setdefault('states_off', Counter())
        for ts in ts_all:
            r_['states_off'][CH[d].state(ts + a.offset, lo, hi)[0]] += 1
        for s_, g_ in S[d]:
            r_['states'][s_] += 1
            if s_ == 'rejected':
                r_['reasons']['+'.join(r for r in ('burst', 'dead', 'npads', 'wide') if CH[d].reason[r][g_])] += 1
                k, info = piece_info(CH[d], g_, p['cluster_r'])
                r_['pieces'][str(min(k, 5))] += 1
                r_['extra'] += info
    for pr, (ta, tb) in TAGS.items():
        key = NAME[pr]
        res = acc['probes'].setdefault(key, {'on': Counter(), 'off': Counter(), 'chip': defaultdict(Counter),
                                             'd_far': Counter(), 'd_found': Counter()})
        lat = J[f'probe_det{pr}']['latency_ns']
        off_xy = J[f'probe_det{pr}']['offset_mm']
        cx, cy, cv, cc, cl = pads[pr]
        for i_ts, ts in enumerate(ts_all):
            sa, ga = S[ta][i_ts]
            sb, gb = S[tb][i_ts]
            px = py = None
            extra = {}
            if sa == 'none' and sb == 'none':
                cls = 'no_tag_both'
            elif sa == 'none':
                cls = 'no_tag_A'
            elif sb == 'none':
                cls = 'no_tag_B'
            elif sa == 'rejected' or sb == 'rejected':
                cls = 'tag_rejected'
                for d_, s_, g_ in ((ta, sa, ga), (tb, sb, gb)):
                    if s_ == 'rejected':
                        extra[NAME[d_]] = '+'.join(r for r in ('burst', 'dead', 'npads', 'wide') if CH[d_].reason[r][g_])
            else:
                f = (Z[pr] - Z[ta]) / (Z[tb] - Z[ta])
                xa, ya, xb, yb = CH[ta].mx[ga], CH[ta].my[ga], CH[tb].mx[gb], CH[tb].my[gb]
                px = xa + f * (xb - xa) + off_xy[0]
                py = ya + f * (yb - ya) + off_xy[1]
                d2 = np.hypot(cx - px, cy - py)
                q = int(np.argmin(d2))
                t_tag = CH[ta].t0[ga]
                chip = int(cv[q])
                if d2[q] >= p['fid_r']:
                    cls = 'outside_area'
                elif not cl[q]:
                    cls = 'dead_pad'
                else:
                    for w, sh in (('on', 0.0), ('off', a.offset)):
                        if not CH[pr].any_hit(ts + lo + sh, ts + hi + sh):
                            c_ = 'probe_silent'
                        else:
                            dn = CH[pr].nearest_live(px, py, ts + lo + sh, ts + hi + sh)
                            if dn > p['probe_r']:
                                c_ = 'probe_far'
                            elif CH[pr].nearest_live(px, py, t_tag + lat - p['win'] + sh,
                                                     t_tag + lat + p['win'] + sh) <= p['probe_r']:
                                c_ = 'found'
                            else:
                                c_ = 'probe_near_late'
                            if w == 'on':
                                b = 'inf' if not np.isfinite(dn) else str(int(min(dn, 300) // 5 * 5))
                                res['d_far' if c_ == 'probe_far' else 'd_found'][b] += 1
                        if w == 'on':
                            cls = c_
                            res['chip'][str(chip)][c_] += 1
                        else:
                            res['off'][c_] += 1
            res['on'][cls] += 1
            ex = acc['examples'].setdefault(key, {}).setdefault(cls, [])
            if len(ex) < a.examples and (cls != 'found' or len(ex) < 2):
                ex.append({'run': J['params']['title'], 'class': cls, 'extra': extra,
                           'pred': None if px is None else [round(float(px), 1), round(float(py), 1)],
                           'hits': {NAME[d]: CH[d].hits(ts - 200, ts + 600, ts) for d in (3, 2, 1)}})
    acc['n_scint'] += len(ts_all)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('jsons', nargs='+')
    ap.add_argument('--out', required=True)
    ap.add_argument('--lo', type=float, default=140)
    ap.add_argument('--hi', type=float, default=310)
    ap.add_argument('--offset', type=float, default=4000)
    ap.add_argument('--examples', type=int, default=4)
    a = ap.parse_args()
    acc = {'n_scint': 0, 'probes': {}, 'rejects': {}, 'examples': {}}
    for j in a.jsons:
        one_run(j, a, acc)
        print(j, acc['n_scint'], flush=True)
    for k, r in acc['probes'].items():
        r['chip'] = {c: dict(v) for c, v in r['chip'].items()}
    json.dump(acc, open(a.out, 'w'), default=lambda o: dict(o) if isinstance(o, Counter) else str(o))


if __name__ == '__main__':
    sys.exit(main())
