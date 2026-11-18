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

OWN-PAD efficiency (per_pad_own, all_own): the same tags, but a HIT only when the predicted pad's OWN channel fired
within --win ns of the latency (its neighbours do not count); random hits of that channel from the --sideband window.
This is the efficiency of one channel, to compare with its own threshold. per_pad_own rows:
[vmm, ch, x, y, n_on, n_off, k_on, k_off, acc_on] (acc_on = on-time tags with a sideband hit on the pad).

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
from spark_veto import WallClock

Z = {3: 320.0, 2: 630.0, 1: 940.0}
NAME = {3: 'det3 P2_IN', 2: 'det2 P2_MID', 1: 'det1 P2_OUT'}
DET = np.zeros(32, int)
for _v, _d in CHAMBER.items():
    DET[_v] = _d


def load(files, clock=None):
    """hits of the capture files, time-ordered; with clock (a spark_veto.WallClock) every file is also anchored to wall
    time by its first packet, so clock(t_ns) gives the DAQ-box epoch seconds of a hit."""
    state = {}
    V, T, C, A = [], [], [], []
    for p in files.split(','):
        recs = read_times(p, state, with_adc=True)
        if clock is not None:
            ts = [r[1] for r in recs if len(r[1])]
            clock.add_file(p, min(x.min() for x in ts) if ts else None)
        for v, t, ch, adc in recs:
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


def time_since_previous(times_by_key, keys, t):
    """for each (key, t): t minus the last time of that key strictly before t (inf if none)."""
    out = np.full(len(t), np.inf)
    for k in set(keys.tolist()):
        tt = times_by_key.get(k)
        if tt is None or not len(tt):
            continue
        m = np.nonzero(keys == k)[0]
        idx = np.searchsorted(tt, t[m]) - 1
        ok = idx >= 0
        out[m[ok]] = t[m[ok]] - tt[idx[ok]]
    return out


def probe(p, ta_d, tb_d, cl, hits, live_pads, a, chip_times=None, chan_times=None):
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

    def est(mask_extra=None, probe_r=a.probe_r, dmin=dmin, smin=smin):
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
    # probe criteria: time window, and a hit anywhere in the probe chamber (position ignored)
    out['vs_win_ns'] = {str(w): est(dmin=nearest_hit(px, py, pt, hx, hy, ht, lat - w, lat + w)) for w in (50, 100, 300, 1000)}
    out['any_position'] = est(probe_r=np.inf)
    # time of the probe hits within --probe-r of the prediction, relative to tag + latency (on minus scaled off)
    edges_t = np.arange(-2000, 2001, 25)
    h_on, h_off = np.zeros(len(edges_t) - 1), np.zeros(len(edges_t) - 1)
    for k in np.nonzero(fid & (on | off))[0]:
        s0, e0 = np.searchsorted(ht, pt[k] + lat - 2000), np.searchsorted(ht, pt[k] + lat + 2000)
        q = np.hypot(hx[s0:e0] - px[k], hy[s0:e0] - py[k]) < a.probe_r
        hh = np.histogram(ht[s0:e0][q] - pt[k] - lat, bins=edges_t)[0]
        if on[k]:
            h_on += hh
        else:
            h_off += hh
    out['probe_dt_hist'] = {'edges_ns': edges_t.tolist(), 'on': h_on.tolist(), 'off_scaled': (r_scale * h_off).tolist()}
    # dead time: time since the previous hit (any ADC, also noise) on the chip / the channel of the predicted pad
    if chip_times is not None:
        t_ref = pt + lat - a.win
        g_chip = time_since_previous(chip_times, cv[near], t_ref)
        g_chan = time_since_previous(chan_times, cv[near] * 64 + cc[near], t_ref)
        gb = (0, 1e3, 1e4, 1e5, 1e6, 1e7, np.inf)
        lab = ('<1us', '1-10us', '10-100us', '0.1-1ms', '1-10ms', '>10ms')
        out['vs_prev_chip_hit'] = {l: est((g_chip >= lo) & (g_chip < hi)) for l, lo, hi in zip(lab, gb[:-1], gb[1:])}
        out['vs_prev_channel_hit'] = {l: est((g_chan >= lo) & (g_chan < hi)) for l, lo, hi in zip(lab, gb[:-1], gb[1:])}
    out['per_chip'] = {str(v): est(cv[near] == v) for v in sorted(set(cv.tolist()))}
    pads = []
    for q in np.unique(near[fid]):
        r = est(near == q)
        pads.append([int(cv[q]), int(cc[q]), round(float(cx[q]), 1), round(float(cy[q]), 1), r['n_on'], r['n_off'],
                     r['k_on'], r['k_off']])
    out['per_pad'] = pads
    if chan_times is not None:
        # own-pad hit: the predicted pad's own channel fired in the window (0 = hit, inf = none, as a "distance")
        key = cv[near] * 64 + cc[near]
        own = np.full(len(px), np.inf)
        own_s = np.full(len(px), np.inf)
        for kk in np.unique(key[fid]):
            tt = chan_times.get(int(kk))
            if tt is None or not len(tt):
                continue
            m = np.nonzero(fid & (key == kk))[0]
            for dst, lo in ((own, lat - a.win), (own_s, lat - a.win + a.sideband)):
                s0 = np.searchsorted(tt, pt[m] + lo)
                e0 = np.searchsorted(tt, pt[m] + lo + 2 * a.win, side='right')
                dst[m[e0 > s0]] = 0.0
        out['all_own'] = est(dmin=own, smin=own_s)
        pads = []
        for q in np.unique(near[fid]):
            mq = fid & (near == q)
            r = est(near == q, dmin=own, smin=own_s)
            pads.append([int(cv[q]), int(cc[q]), round(float(cx[q]), 1), round(float(cy[q]), 1), r['n_on'], r['n_off'],
                         r['k_on'], r['k_off'], int((mq & on & (own_s < np.inf)).sum())])
        out['per_pad_own'] = pads
        if a.dump_tags:
            # per-tag record for the resolution/efficiency unfolding: every fiducial on/off tag, its predicted point,
            # the predicted pad, own-pad hit (window / sideband), the tag-pair slope, every live probe pad that
            # fired within --dump-r mm of the prediction in the window (fired_*) and in the sideband window (side_*),
            # as pad centres flattened with the tag index; dt_scint for the scintillator-confirmed selection
            sel = np.nonzero(fid & (on | off))[0]
            def fired_lists(t0):
                """flattened (tag index, pad centre) of the live probe pads fired within --dump-r mm of the prediction
                in the window [t0 - win, t0 + win] relative to the tag time."""
                s0 = np.searchsorted(ht, pt[sel] + t0 - a.win)
                e0 = np.searchsorted(ht, pt[sel] + t0 + a.win, side='right')
                fi, fx, fy = [], [], []
                for n_, (k, s1, e1) in enumerate(zip(sel, s0, e0)):
                    if e1 > s1:
                        q = np.hypot(hx[s1:e1] - px[k], hy[s1:e1] - py[k]) < a.dump_r
                        xy = np.unique(np.c_[hx[s1:e1][q], hy[s1:e1][q]], axis=0)
                        fi += [n_] * len(xy); fx += xy[:, 0].tolist(); fy += xy[:, 1].tolist()
                return np.array(fi, int), np.array(fx), np.array(fy)
            fi, fx, fy = fired_lists(lat)
            si, sx, sy = fired_lists(lat + a.sideband)   # random hits: same tags, window shifted by --sideband
            slope = np.hypot(xb[j] - xa[i], yb[j] - ya[i]) / abs(Z[tb_d] - Z[ta_d])
            # signed time from the nearest scintillator coincidence (H0 VMM 0 ch 40) to the tag (+ ~210 ns for muons)
            sc = a._scint
            ks = np.clip(np.searchsorted(sc, pt[sel]), 1, max(len(sc) - 1, 1))
            dts = (np.where(np.abs(pt[sel] - sc[ks - 1]) < np.abs(pt[sel] - sc[ks]), pt[sel] - sc[ks - 1], pt[sel] - sc[ks])
                   if len(sc) > 1 else np.full(len(sel), np.inf))
            a._dump[f'p{p}'] = dict(px=px[sel], py=py[sel], on=on[sel], vmm=cv[near[sel]], ch=cc[near[sel]], dt_scint=dts,
                                    cx=cx[near[sel]], cy=cy[near[sel]], own=np.isfinite(own[sel]),
                                    own_s=np.isfinite(own_s[sel]), slope=slope[sel],
                                    tag_adc=np.minimum(aa[i], ab[j])[sel], fired_i=fi, fired_x=fx, fired_y=fy,
                                    side_i=si, side_x=sx, side_y=sy, r_scale=np.array(r_scale),
                                    dump_r=np.array(a.dump_r))
            if a.dump_diag:
                a._dump[f'p{p}'].update(diag_fields(p, sel, px, py, pt, lat, hx, hy, ht, g_chip, g_chan, a))
    return out


def diag_fields(p, sel, px, py, pt, lat, hx, hy, ht, g_chip, g_chan, a):
    """--dump-diag: per tag (same order as the other per-tag arrays) its time, wall time, the probe-chamber activity
    around it and the dead-time gaps; flattened (fired_nr_*) every live probe hit within --diag-r mm of the prediction
    and within +-(--diag-wide) ns of tag + latency: its time from tag + latency, ADC, distance and vmm*64+ch.
    The fired_ prefix makes channel_efficiency.select_tags drop the flattened arrays instead of mis-indexing them."""
    ha, hk = a._hits_adc[p], a._hits_key[p]
    W, n = a.diag_wide, len(sel)
    pts = pt[sel] + lat
    s0 = np.searchsorted(ht, pts - W)
    e0 = np.searchsorted(ht, pts + W, side='right')
    n_any = np.zeros(n, int)                    # live probe hits anywhere in the chamber within +-win
    d_win = np.full(n, np.inf)                  # nearest of them to the prediction (mm)
    fi, fdt, fadc, fd, fk = [], [], [], [], []
    for n_, (k, s1, e1) in enumerate(zip(sel, s0, e0)):
        if e1 <= s1:
            continue
        dt_ = ht[s1:e1] - pts[n_]
        dd = np.hypot(hx[s1:e1] - px[k], hy[s1:e1] - py[k])
        inw = np.abs(dt_) <= a.win
        n_any[n_] = int(inw.sum())
        if n_any[n_]:
            d_win[n_] = float(dd[inw].min())
        q = dd < a.diag_r
        if q.any():
            fi += [n_] * int(q.sum()); fdt += dt_[q].tolist(); fadc += ha[s1:e1][q].tolist()
            fd += dd[q].tolist(); fk += hk[s1:e1][q].tolist()
    # chamber activity around the tag: all live hits and saturated hits (ADC >= --sat-adc) within +-100 us
    sat_t = ht[ha >= a.sat_adc]
    w = 1e5
    n100 = np.searchsorted(ht, pts + w, side='right') - np.searchsorted(ht, pts - w)
    s100 = np.searchsorted(sat_t, pts + w, side='right') - np.searchsorted(sat_t, pts - w)
    return dict(t_tag=pt[sel], wall=a._clock(pt[sel]) if a._clock.anchors else np.full(n, np.nan),
                n_any=n_any, d_win=d_win, n100us=n100, sat100us=s100,
                prev_chip=g_chip[sel], prev_chan=g_chan[sel], latency=np.array(lat), win=np.array(a.win),
                fired_nr_i=np.array(fi, int), fired_nr_dt=np.array(fdt), fired_nr_adc=np.array(fadc, int),
                fired_nr_d=np.array(fd), fired_nr_key=np.array(fk, int))


def second_rates(V, T, A, clock, sat_adc):
    """per wall-clock second and chamber: live+dead hits and saturated hits (for the spark-flag check)."""
    ws = clock(T)
    t0 = np.floor(ws.min())
    edges = np.arange(t0, np.ceil(ws.max()) + 1)
    out = {'t0': np.array(t0)}
    for d in (3, 2, 1):
        m = DET[V] == d
        out[f'hits_p{d}'] = np.histogram(ws[m], edges)[0].astype(np.int32)
        out[f'sat_p{d}'] = np.histogram(ws[m & (A >= sat_adc)], edges)[0].astype(np.int32)
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
    ap.add_argument('--dump-tags', default='', help='write the per-tag records of every probe to this .npz')
    ap.add_argument('--dump-r', type=float, default=80, help='fired pads kept within this distance (mm) in the dump')
    ap.add_argument('--scint-ch', default='0:40', help='vmm:ch of the scintillator coincidence (dt_scint in the dump)')
    ap.add_argument('--dump-diag', action='store_true', help='also dump the diagnostic fields (diag_fields) and the '
                    'per-second chamber rates (rates_*): wall time of every tag, probe hits near the prediction with '
                    'time and ADC, chamber activity, dead-time gaps')
    ap.add_argument('--diag-wide', type=float, default=2000, help='+-ns around tag + latency for the fired_nr_* hits')
    ap.add_argument('--diag-r', type=float, default=15, help='mm around the prediction for the fired_nr_* hits')
    ap.add_argument('--sat-adc', type=int, default=1000, help='ADC from which a hit counts as saturated')
    a = ap.parse_args()
    a._dump = {}
    a._clock = WallClock()

    X, Y = O.Mapper(a.conn_shift).xy_table(O.parse_orders(a.order))
    veto = parse_veto(a.veto)
    dead = {int(v) for v in a.dead_vmm.split(',') if v}
    V, T, C, A = load(a.files, a._clock if a.dump_diag else None)
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
    a._hits_adc, a._hits_key = {}, {}
    for d in (3, 2, 1):
        m = (DET[V] == d) & live[V, C]
        hits[d] = (HX[m], HY[m], T[m])
        if a.dump_diag:
            a._hits_adc[d], a._hits_key[d] = A[m], V[m] * 64 + C[m]
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
    key = V * 64 + C
    o = np.argsort(key, kind='stable')
    ks, starts = np.unique(key[o], return_index=True)
    To = T[o]
    chan_times = {int(k): To[s0:e0] for k, s0, e0 in zip(ks, starts, list(starts[1:]) + [len(o)])}
    chip_times = {int(v): T[V == v] for v in np.unique(V)}
    sv, sc_ = map(int, a.scint_ch.split(':'))
    a._scint = np.sort(T[(V == sv) & (C == sc_)])
    for p, tags in ((2, (3, 1)), (1, (3, 2)), (3, (2, 1))):
        r = probe(p, tags[0], tags[1], cl, hits, live_pads, a, chip_times, chan_times)
        res[f'probe_det{p}'] = r
        s = r['all']
        print(f'{NAME[p]} (tags {r["tags"][0]} + {r["tags"][1]}): on {s["n_on"]} off {s["n_off"]} -> tags {s["n_sig"]} '
              f'(purity {s.get("purity")}), eff {s["eff"]} +- {s.get("err")} (raw on-time {s.get("eff_raw_on")}, acc '
              f'{s.get("acc")}); dt0 {r["dt0_ns"]:.0f} ns, lat {r["latency_ns"]:.0f} ns, offset '
              f'{np.round(r["offset_mm"], 1)} mm', flush=True)
        if 'all_own' in r:
            print(f'   OWN pad only: eff {r["all_own"]["eff"]} +- {r["all_own"].get("err")} (acc {r["all_own"].get("acc")})')
        print('   vs probe-r:', {k: v['eff'] for k, v in r['vs_probe_r'].items()})
        print('   vs window +-ns:', {k: v['eff'] for k, v in r['vs_win_ns'].items()},
              '| hit anywhere in the chamber:', r['any_position']['eff'])
        print('   vs time since previous CHIP hit: eff (tags):',
              {k: (v['eff'], v['n_sig']) for k, v in r['vs_prev_chip_hit'].items()})
        print('   vs time since previous CHANNEL hit: eff (tags):',
              {k: (v['eff'], v['n_sig']) for k, v in r['vs_prev_channel_hit'].items()})
        print('   vs smaller tag ADC >= cut: eff (tags):', {k: (v['eff'], v['n_sig']) for k, v in r['vs_tag_adc'].items()})
        print('   vs tan(theta) of the tag pair: eff (tags):',
              {k: (v['eff'], v['n_sig']) for k, v in r['vs_tan_theta'].items()})
        print('   per chip:', {k: (v['eff'], v.get('err')) for k, v in r['per_chip'].items() if v['eff'] is not None})
    vars(a).pop('_scint')
    for k in ('_hits_adc', '_hits_key'):
        vars(a).pop(k)
    clock = vars(a).pop('_clock')
    if a.dump_diag and a.dump_tags and clock.anchors:
        a._dump['rates'] = second_rates(V, T, A, clock, a.sat_adc)
    dump = vars(a).pop('_dump')  # res['params'] is vars(a): keep the arrays out of the json
    if a.dump_tags:
        np.savez_compressed(a.dump_tags, **{f'{p}_{k}': v for p, d in dump.items() for k, v in d.items()})
    json.dump(res, open(a.out_json, 'w'), indent=1, default=float)


if __name__ == '__main__':
    sys.exit(main())
