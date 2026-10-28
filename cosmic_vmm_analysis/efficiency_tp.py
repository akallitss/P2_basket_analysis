#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tag-and-probe efficiency of each chamber for ONE run, with the fake tags subtracted from the data.

Without an external trigger most coincidences of a clean cluster in two chambers are noise, so a plain tag-and-probe
counts every fake tag as a missed probe. Here the tag pairs are built from the time difference of the two clusters:

  * clusters per chamber: hits of that chamber grouped in time (gap --gap ns); a cluster is CLEAN when it has 1..--max-tag
    mapped live pads, all within --cluster-r mm of their mean, and its time is not inside a whole-stack burst (a group of
    all chambers with more than --burst hits); position = mean pad centre, time = earliest hit;
  * tag pairs (chambers a, b) with dt = t_b - t_a: ON-time |dt - dt0| < --on ns around the muon peak dt0 (median of
    the pairs within +-200 ns), OFF-time --off-lo < |dt - dt0| < --off-hi ns: random pairs, scaled by
    r = 2 --on / (2 (--off-hi - --off-lo));
  * prediction in the probe chamber: straight line through the two tag positions (z 320 / 630 / 940 mm), plus the
    median probe offset measured on on-time pairs (probe hit within 60 mm);
  * fiducial: the nearest pad to the prediction is a LIVE pad of the probe, within --fid-r mm (live = mapped, not in
    --dead-vmm, and with at least --min-hits hits in the run, which removes masked and dead channels);
  * HIT: a live probe pad within --probe-r mm of the prediction, hit within --win ns of the probe latency (median
    probe-minus-tag time on on-time pairs); random probe hits: the same search --sideband ns later (acc);
  * eff = (k_on - r k_off) / (n_on - r n_off), then corrected for random probe hits: (eff - acc) / (1 - acc).

Map: read-out connectors shifted by --conn-shift, per-connector channel order from --order (orientation_check.Mapper).
--veto vmm:ch,...: hits of these channels are dropped and so are the clusters of their chamber group (a discharging
channel and what it induces).

Output OUT.json: per probe chamber the counts, eff with error, acc, the dependence on --probe-r and on the inclination
tan(theta) of the tag pair, and on the pulse height of the tags (the smaller of the two clusters' largest ADC): noise
firing two chambers at once passes the off-time subtraction and has a cosmic-like angle distribution (run_76 shows both),
but sits near threshold, so the efficiency rising with the tag ADC cut and levelling off measures the real one; per chip
and per pad.

Usage: efficiency_tp.py f1.pcapng[,f2,...] OUT.json --conn-shift -1 --order 'P2_IN:3=bot:rev|top:lin' ...
       [--veto 2:27,6:44] [--dead-vmm 4,7,9] [--title run_124]
"""
import argparse
import json
import sys

import numpy as np

import orientation_check as O
from coincidence_check import read_times, CHAMBER
from straight_track_rate import parse_veto

Z = {3: 320.0, 2: 630.0, 1: 940.0}
NAME = {3: 'det3 P2_IN', 2: 'det2 P2_MID', 1: 'det1 P2_OUT'}
DET = np.zeros(32, int)
for _v, _d in CHAMBER.items():
    DET[_v] = _d


def load(files):
    state = {}
    V, T, C, A = [], [], [], []
    for p in files.split(','):
        for v, t, ch, adc in read_times(p, state, with_adc=True):
            V.append(v); T.append(t); C.append(ch); A.append(adc)
    V, T, C, A = np.concatenate(V), np.concatenate(T), np.concatenate(C), np.concatenate(A)
    o = np.argsort(T, kind='stable')
    return V[o].astype(int), T[o], C[o].astype(int), A[o].astype(int)


def split_groups(t, gap):
    return np.concatenate([[0], np.cumsum(np.diff(t) > gap)]) if len(t) else np.zeros(0, int)


def chamber_clusters(d, V, T, C, A, X, Y, live, bad_hit, burst_hit, a):
    """clean clusters of chamber d: arrays x, y, t, largest ADC."""
    m = DET[V] == d
    v, t, c, x, y, adc = V[m], T[m], C[m], X[m], Y[m], A[m]
    ok_hit = live[v, c]
    g = split_groups(t, a.gap)
    ng = int(g[-1]) + 1 if len(g) else 0
    n_live = np.bincount(g, weights=ok_hit, minlength=ng)
    n_all = np.bincount(g, minlength=ng)
    bad = np.bincount(g, weights=bad_hit[m] | burst_hit[m], minlength=ng) > 0
    sx = np.bincount(g, weights=np.where(ok_hit, x, 0), minlength=ng)
    sy = np.bincount(g, weights=np.where(ok_hit, y, 0), minlength=ng)
    with np.errstate(invalid='ignore', divide='ignore'):
        mx, my = sx / n_live, sy / n_live
    r = np.hypot(x - mx[g], y - my[g])
    spread = np.zeros(ng)
    np.maximum.at(spread, g, np.where(ok_hit, r, 0))
    t0 = np.full(ng, np.inf)
    np.minimum.at(t0, g, t)
    amax = np.zeros(ng, int)
    np.maximum.at(amax, g, adc)
    clean = (n_live >= 1) & (n_all <= a.max_tag) & (n_live == n_all) & (spread <= a.cluster_r) & ~bad
    return mx[clean], my[clean], t0[clean], amax[clean]


def pairs(ta, tb, lo, hi):
    """all index pairs (i, j) with lo <= tb[j] - ta[i] <= hi (both sorted)."""
    s = np.searchsorted(tb, ta + lo)
    e = np.searchsorted(tb, ta + hi, side='right')
    n = e - s
    i = np.repeat(np.arange(len(ta)), n)
    j = np.concatenate([np.arange(a0, b0) for a0, b0 in zip(s, e)]) if n.sum() else np.zeros(0, int)
    return i, j


def nearest_hit(px, py, pt, hx, hy, ht, lo, hi):
    """distance from (px, py) to the nearest probe hit with time in [pt + lo, pt + hi] (inf if none)."""
    s = np.searchsorted(ht, pt + lo)
    e = np.searchsorted(ht, pt + hi, side='right')
    out = np.full(len(px), np.inf)
    for k in np.nonzero(e > s)[0]:
        out[k] = np.min(np.hypot(hx[s[k]:e[k]] - px[k], hy[s[k]:e[k]] - py[k]))
    return out


def probe(p, ta_d, tb_d, cl, hits, live_pads, a):
    (xa, ya, ta, aa), (xb, yb, tb, ab) = cl[ta_d], cl[tb_d]
    i, j = pairs(ta, tb, -a.off_hi - 300, a.off_hi + 300)
    dt = tb[j] - ta[i]
    near0 = np.abs(dt) < 200
    dt0 = float(np.median(dt[near0])) if near0.sum() > 20 else 0.0
    on = np.abs(dt - dt0) < a.on
    off = (np.abs(dt - dt0) > a.off_lo) & (np.abs(dt - dt0) < a.off_hi)
    r_scale = (2 * a.on) / (2 * (a.off_hi - a.off_lo))
    f = (Z[p] - Z[ta_d]) / (Z[tb_d] - Z[ta_d])
    px = xa[i] + f * (xb[j] - xa[i])
    py = ya[i] + f * (yb[j] - ya[i])
    pt = ta[i]
    hx, hy, ht = hits[p]
    # latency and offset from on-time pairs
    sel = np.nonzero(on)[0]
    s = np.searchsorted(ht, pt[sel] - 300)
    e = np.searchsorted(ht, pt[sel] + 300, side='right')
    lat_l, ox, oy = [], [], []
    for k, s0, e0 in zip(sel, s, e):
        if e0 > s0:
            dd = np.hypot(hx[s0:e0] - px[k], hy[s0:e0] - py[k])
            q = int(np.argmin(dd))
            if dd[q] < 60:
                lat_l.append(ht[s0 + q] - pt[k]); ox.append(hx[s0 + q] - px[k]); oy.append(hy[s0 + q] - py[k])
    lat = float(np.median(lat_l)) if lat_l else 0.0
    off_xy = (float(np.median(ox)), float(np.median(oy))) if len(ox) > 20 else (0.0, 0.0)
    px, py = px + off_xy[0], py + off_xy[1]
    # fiducial: nearest pad of the probe is live and close
    cx, cy, cv, cc, clive = live_pads[p]
    near = np.empty(len(px), int)
    dpad = np.empty(len(px))
    for s0 in range(0, len(px), 4000):
        d2 = (px[s0:s0 + 4000, None] - cx[None]) ** 2 + (py[s0:s0 + 4000, None] - cy[None]) ** 2
        near[s0:s0 + 4000] = np.argmin(d2, 1)
        dpad[s0:s0 + 4000] = np.sqrt(d2.min(1))
    fid = (dpad < a.fid_r) & clive[near]
    dmin = nearest_hit(px, py, pt, hx, hy, ht, lat - a.win, lat + a.win)
    smin = nearest_hit(px, py, pt, hx, hy, ht, lat - a.win + a.sideband, lat + a.win + a.sideband)

    def est(mask_extra=None, probe_r=a.probe_r):
        m = fid if mask_extra is None else fid & mask_extra
        n_on, n_off = int((m & on).sum()), int((m & off).sum())
        k_on, k_off = int((m & on & (dmin < probe_r)).sum()), int((m & off & (dmin < probe_r)).sum())
        acc_n = int((m & on & (smin < probe_r)).sum())
        n_sig = n_on - r_scale * n_off
        k_sig = k_on - r_scale * k_off
        acc = acc_n / n_on if n_on else float('nan')
        if n_sig <= 0:
            return {'n_on': n_on, 'n_off': n_off, 'k_on': k_on, 'k_off': k_off, 'n_sig': n_sig, 'eff': None}
        e_raw = k_sig / n_sig
        if not acc < 1:
            return {'n_on': n_on, 'n_off': n_off, 'k_on': k_on, 'k_off': k_off, 'n_sig': n_sig, 'acc': acc, 'eff': None}
        e = (e_raw - acc) / (1 - acc)
        # error: binomial on the signal plus the Poisson spread of the subtraction
        var = (e_raw * (1 - e_raw) / n_sig if 0 <= e_raw <= 1 else 0) + \
              (r_scale ** 2 * (k_off + e_raw ** 2 * n_off)) / n_sig ** 2
        return {'n_on': n_on, 'n_off': n_off, 'k_on': k_on, 'k_off': k_off, 'n_sig': round(n_sig, 1),
                'purity': round(n_sig / n_on, 3) if n_on else None, 'eff_raw_on': round(k_on / n_on, 4) if n_on else None,
                'acc': round(acc, 4), 'eff': round(e, 4), 'err': round(float(np.sqrt(var)) / (1 - acc), 4)}

    out = {'tags': [NAME[ta_d], NAME[tb_d]], 'dt0_ns': dt0, 'latency_ns': lat, 'offset_mm': off_xy,
           'r_scale': r_scale, 'dt_hist_ns': np.histogram(dt - dt0, bins=np.arange(-2000, 2001, 20))[0].tolist(),
           'all': est(), 'vs_probe_r': {str(r): est(probe_r=r) for r in (5, 10, 15, 20, 30, 50)}}
    slope = np.hypot(xb[j] - xa[i], yb[j] - ya[i]) / abs(Z[tb_d] - Z[ta_d])
    edges = (0, 0.1, 0.2, 0.3, 0.5, 0.8, 1.5, 10)
    out['vs_tan_theta'] = {f'{lo}-{hi}': est((slope >= lo) & (slope < hi)) for lo, hi in zip(edges[:-1], edges[1:])}
    tag_adc = np.minimum(aa[i], ab[j])
    out['vs_tag_adc'] = {str(c): est(tag_adc >= c) for c in (0, 100, 150, 200, 300, 400, 600, 1000)}
    out['per_chip'] = {str(v): est(cv[near] == v) for v in sorted(set(cv.tolist()))}
    pads = []
    for q in np.unique(near[fid]):
        r = est(near == q)
        pads.append([int(cv[q]), int(cc[q]), round(float(cx[q]), 1), round(float(cy[q]), 1), r['n_on'], r['n_off'],
                     r['k_on'], r['k_off']])
    out['per_pad'] = pads
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('files')
    ap.add_argument('out_json')
    ap.add_argument('--conn-shift', type=int, default=-1)
    ap.add_argument('--order', action='append')
    ap.add_argument('--veto', default='')
    ap.add_argument('--dead-vmm', default='', help='chips not taking data (e.g. parked at sdt 1023): 4,7,9')
    ap.add_argument('--min-hits', type=int, default=1)
    ap.add_argument('--gap', type=float, default=100)
    ap.add_argument('--burst', type=int, default=36)
    ap.add_argument('--max-tag', type=int, default=6)
    ap.add_argument('--cluster-r', type=float, default=15)
    ap.add_argument('--probe-r', type=float, default=15)
    ap.add_argument('--fid-r', type=float, default=9)
    ap.add_argument('--win', type=float, default=100)
    ap.add_argument('--on', type=float, default=60)
    ap.add_argument('--off-lo', type=float, default=300)
    ap.add_argument('--off-hi', type=float, default=1900)
    ap.add_argument('--sideband', type=float, default=4000)
    ap.add_argument('--title', default='')
    a = ap.parse_args()

    X, Y = O.Mapper(a.conn_shift).xy_table(O.parse_orders(a.order))
    veto = parse_veto(a.veto)
    dead = {int(v) for v in a.dead_vmm.split(',') if v}
    V, T, C, A = load(a.files)
    span = float(T[-1] - T[0]) * 1e-9
    counts = np.zeros((32, 64), int)
    np.add.at(counts, (V, C), 1)
    live = np.isfinite(X) & (counts >= a.min_hits)
    for v in dead:
        live[v] = False
    for v, c in veto:
        live[v, c] = False
    vmask = np.zeros((32, 64), bool)
    for v, c in veto:
        vmask[v, c] = True
    bad_hit = vmask[V, C]
    g = split_groups(T, a.gap)
    burst_hit = (np.bincount(g)[g] > a.burst)
    HX, HY = X[V, C], Y[V, C]
    cl = {d: chamber_clusters(d, V, T, C, A, HX, HY, live, bad_hit, burst_hit, a) for d in (3, 2, 1)}
    hits, live_pads = {}, {}
    for d in (3, 2, 1):
        m = (DET[V] == d) & live[V, C]
        hits[d] = (HX[m], HY[m], T[m])
        vv, cc = np.nonzero(np.isfinite(X))
        k = DET[vv] == d
        live_pads[d] = (X[vv[k], cc[k]], Y[vv[k], cc[k]], vv[k], cc[k], live[vv[k], cc[k]])
    res = {'title': a.title, 'files': len(a.files.split(',')), 'n_hits': int(len(T)), 'span_s': round(span, 1),
           'params': vars(a), 'clean_clusters': {NAME[d]: int(len(cl[d][2])) for d in cl},
           'live_channels': {NAME[d]: int(live_pads[d][4].sum()) for d in live_pads},
           'dead_or_masked': {NAME[d]: [f'{int(v)}:{int(c)}' for v, c, l in zip(live_pads[d][2], live_pads[d][3],
                                                                                 live_pads[d][4]) if not l]
                              for d in live_pads}}
    print(f'{a.title}: {len(T)} hits in {span / 60:.1f} min; clean clusters', res['clean_clusters'],
          '; live channels', res['live_channels'])
    for p, tags in ((2, (3, 1)), (1, (3, 2)), (3, (2, 1))):
        r = probe(p, tags[0], tags[1], cl, hits, live_pads, a)
        res[f'probe_det{p}'] = r
        s = r['all']
        print(f'{NAME[p]} (tags {r["tags"][0]} + {r["tags"][1]}): on {s["n_on"]} off {s["n_off"]} -> tags {s["n_sig"]} '
              f'(purity {s.get("purity")}), eff {s["eff"]} +- {s.get("err")} (raw on-time {s.get("eff_raw_on")}, acc '
              f'{s.get("acc")}); dt0 {r["dt0_ns"]:.0f} ns, lat {r["latency_ns"]:.0f} ns, offset '
              f'{np.round(r["offset_mm"], 1)} mm', flush=True)
        print('   vs probe-r:', {k: v['eff'] for k, v in r['vs_probe_r'].items()})
        print('   vs smaller tag ADC >= cut: eff (tags):', {k: (v['eff'], v['n_sig']) for k, v in r['vs_tag_adc'].items()})
        print('   vs tan(theta) of the tag pair: eff (tags):',
              {k: (v['eff'], v['n_sig']) for k, v in r['vs_tan_theta'].items()})
        print('   per chip:', {k: (v['eff'], v.get('err')) for k, v in r['per_chip'].items() if v['eff'] is not None})
    json.dump(res, open(a.out_json, 'w'), indent=1, default=float)


if __name__ == '__main__':
    sys.exit(main())
