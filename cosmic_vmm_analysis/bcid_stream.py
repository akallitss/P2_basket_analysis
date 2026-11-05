#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Where do the hits of a reconstructed muon sit in the read-out stream, and what are their BCIDs?

The FEC sends UDP frames of 6-byte words: SRS markers (one per VMM, carrying the 42-bit FEC timestamp in 25 ns
ticks) interleaved with hit words (channel, ADC, TDC, 12-bit Gray BCID and the 5-bit offset = number of 4096-BC
overflows since that VMM's last marker). The hit time is, as in coincidence_check.read_times and the QA,
    t = srs * TICK + (offset*4096 + bcid + 1) * 22.5 ns - tdc * 60/255 ns
so  t - srs*TICK = the hit time relative to its SRS marker ("chip time"), and offset*4096 + bcid the same in BC.
TICK (--srs-tick) is 22.5 ns: the marker timestamp advances by 2.6700e9 ticks per 60.08 s of pcap wall time (runs
121 / 130 / 196), and the markers come every 65536 ticks = 16 x 4096 BC, i.e. the FEC timestamp counts the BC clock.
coincidence_check.read_times and the QA use 25 ns, which stretches absolute times by 11 % (rates and the wall-clock
mapping of the spark veto), but not time differences inside one marker period.

This script keeps, for every hit word, its place in the stream: global hit number (seq, all hit words incl. the
neighbour hits, flag 0), UDP frame number (pkt), word position in the frame (w), number of words since the marker
of its VMM (wsm), plus bcid / offset / tdc / srs. It then finds the chamber-only straight tracks of
paddle_efficiency.py (det2 seed, det3 and det1 within +- --half-win ns, clean clusters, det2 within --tol mm of the
det3-det1 midpoint) and, per track and chamber, writes where the cluster hits are in the stream and their BCIDs:

  per chamber (order det3, det2, det1): first (earliest) cluster hit -> seq, pkt, w, wsm, vmm, bcid, offset, srs,
      chip time; min / max seq of the cluster hits
  per track: seq span of all cluster hits of the three chambers, frame span, number of OTHER hits inside that span
      (over-threshold chamber hits, neighbour hits, hybrid-0 paddle hits), the stream order of the three chambers
  chance reference: the same span for the three chambers' first hits --sideband ns later (unrelated hits)
  paddles (--paddles): nearest pulse to t2 - --latency within --win ns: dt, seq - seq(det2), bcid, srs

Global (all hits): word-position, words-since-marker and chip-time histograms of the chamber over-threshold hits;
"lateness" of every chamber hit in the stream (max time of the hits already read - its own time, i.e. how far the
stream is out of time order); marker period and marker-timestamp agreement between VMMs (first file).

Output OUT.npz (per track arrays + histograms) and OUT.json (summary). Figures: bcid_stream_plots.py.

Usage: bcid_stream.py f1.pcapng[,...] OUT.npz --conn-shift -1 --order 'P2_IN:3=bot:rev|top:lin' --order ...
       [--bad-times veto.json] [--max-files N] [--title run_196]
"""
import argparse
import json
import struct
import sys

import numpy as np

import orientation_check as O
from coincidence_check import gray
from straight_track_rate import parse_veto
from spark_veto import WallClock, load_veto, in_intervals

DET = np.zeros(32, int)
DET[2:8], DET[8:14], DET[14:20] = 3, 2, 1
ORDER = (3, 2, 1)
BC = 22.5


def read_stream(path, state, cnt, mk_log=None):
    """per frame arrays of every hit word with its stream position. state[vmm] = (srs, global word index) of the
    last marker; cnt = [frames, words] running counters (updated). mk_log: list to append (vmm, srs, gw) markers."""
    out = []
    with open(path, 'rb', buffering=16 << 20) as f:
        head = f.read(12)
        if len(head) < 12 or struct.unpack('<I', head[8:12])[0] != 0x1A2B3C4D:
            raise SystemExit(f'{path}: not a little-endian pcapng')
        f.seek(0)
        while True:
            hdr = f.read(8)
            if len(hdr) < 8:
                break
            bt, bl = struct.unpack('<II', hdr)
            body = f.read(bl - 8)
            if len(body) < bl - 8:
                break
            if bt != 6:
                continue
            caplen = struct.unpack('<I', body[12:16])[0]
            data = body[20:20 + caplen]
            off = 14
            if data[12:14] == b'\x81\x00':
                off = 18
            elif data[12:14] != b'\x08\x00':
                continue
            if len(data) < off + 28 or data[off + 9] != 17:
                continue
            p = data[off + (data[off] & 0x0F) * 4 + 8:]
            if len(p) < 22 or p[4:7] != b'VM3':
                continue
            nw = (len(p) - 22 + 5) // 6
            b = np.frombuffer(p, dtype=np.uint8, count=6 * nw, offset=16).reshape(-1, 6).astype(np.int64)
            d1 = ((b[:, 0] << 24) | (b[:, 1] << 16) | (b[:, 2] << 8) | b[:, 3])
            d2 = ((b[:, 4] << 8) | b[:, 5])
            ishit = (d2 & 0x8000) != 0
            hv = (d1 >> 22) & 0x1F
            mid = (d2 >> 10) & 0x1F
            msrs = (d1 << 10) | (d2 & 0x3FF)
            gw0 = cnt[1]
            pos_h = np.nonzero(ishit)[0]
            srs_h = np.zeros(len(pos_h), np.int64)
            mgw_h = np.full(len(pos_h), -1, np.int64)
            mpos_all = np.nonzero(~ishit)[0]
            if mk_log is not None and len(mpos_all):
                mk_log.append(np.c_[mid[mpos_all], msrs[mpos_all], gw0 + mpos_all])
            for vv in np.unique(np.concatenate([hv[pos_h], mid[mpos_all]])):
                mpos = mpos_all[mid[mpos_all] == vv]
                hmask = hv[pos_h] == vv
                base = state.get(int(vv), (0, -1))
                if hmask.any():
                    k = np.searchsorted(mpos, pos_h[hmask]) - 1
                    if len(mpos):
                        kk = np.maximum(k, 0)
                        srs_h[hmask] = np.where(k >= 0, msrs[mpos][kk], base[0])
                        mgw_h[hmask] = np.where(k >= 0, gw0 + mpos[kk], base[1])
                    else:
                        srs_h[hmask], mgw_h[hmask] = base
                if len(mpos):
                    state[int(vv)] = (int(msrs[mpos[-1]]), int(gw0 + mpos[-1]))
            if len(pos_h):
                h1, h2 = d1[pos_h], d2[pos_h]
                offs = (h1 >> 27) & 0x1F
                offs = np.where(offs < 16, offs, offs - 32)
                bcid = gray(h1 & 0xFFF)
                tdc = h2 & 0xFF
                gw = gw0 + pos_h
                out.append({'vmm': hv[pos_h].astype(np.int8), 'ch': ((h2 >> 8) & 0x3F).astype(np.int8),
                            'flag': ((h2 >> 14) & 1).astype(bool), 'adc': ((h1 >> 12) & 0x3FF).astype(np.int16),
                            'bcid': bcid.astype(np.int16), 'offs': offs.astype(np.int8), 'tdc': tdc.astype(np.uint8),
                            'srs': srs_h, 'pkt': np.full(len(pos_h), cnt[0], np.int32), 'w': pos_h.astype(np.int16),
                            'wsm': np.where(mgw_h >= 0, gw - mgw_h, -1).astype(np.int32)})
            cnt[0] += 1
            cnt[1] += nw
    return out


def hist(x, edges):
    return np.histogram(x, bins=edges)[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('files')
    ap.add_argument('out')
    ap.add_argument('--paddles', default='0:44,0:45,0:41')
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
    ap.add_argument('--max-files', type=int, default=0)
    ap.add_argument('--title', default='')
    ap.add_argument('--bad-times', default='')
    ap.add_argument('--srs-tick', type=float, default=22.5, help='SRS marker timestamp tick [ns] (22.5 measured; the old decoder used 25)')
    ap.add_argument('--examples', type=int, default=12, help='stream excerpts of this many straight tracks in the json')
    a = ap.parse_args()

    X, Y = O.Mapper(a.conn_shift).xy_table(O.parse_orders(a.order))
    files = a.files.split(',')
    if a.max_files:
        files = files[:a.max_files]
    state, cnt, recs, mk_log = {}, [0, 0], [], []
    clock = WallClock()
    for i, p in enumerate(files):
        r = read_stream(p, state, cnt, mk_log if i == 0 else None)
        recs += r
        # wall-clock anchor: first valid hit time of the file
        t0 = [x['srs'][x['srs'] > 0] * a.srs_tick for x in r]
        t0 = [x.min() for x in t0 if len(x)]
        clock.add_file(p, min(t0) if t0 else None)
    S = {k: np.concatenate([r[k] for r in recs]) for k in recs[0]}
    del recs
    n_all = len(S['vmm'])
    vmm, ch = S['vmm'].astype(int), S['ch'].astype(int)
    chip = (S['offs'].astype(np.float64) * 4096 + S['bcid'] + 1) * BC - S['tdc'] * 60.0 / 255.0
    T = S['srs'] * a.srs_tick + chip
    valid = S['srs'] > 0
    det = DET[vmm]
    print(f'{len(files)} files, {cnt[0]} frames, {cnt[1]} words, {n_all} hit words '
          f'({S["flag"].sum()} over threshold, {(~S["flag"]).sum()} neighbour)')

    # ---- global stream properties of the chamber over-threshold hits
    ot = S['flag'] & valid & (det > 0)
    # lateness: how much earlier than the latest hit already in the stream (all valid hits, stream order)
    tv = np.where(valid, T, -np.inf)
    runmax = np.maximum.accumulate(tv)
    late = np.r_[0.0, runmax[:-1]] - T          # > 0: this hit is older than one already read
    late = np.where(valid, late, np.nan)
    # per VMM: is the stream time-ordered within one VMM?
    inv_vmm = {}
    for v in np.unique(vmm[ot]):
        m = np.nonzero(S['flag'] & valid & (vmm == v))[0]
        dt = np.diff(T[m])
        inv_vmm[int(v)] = {'n': int(len(m)), 'frac_backwards': float(np.mean(dt < 0)) if len(dt) else None,
                           'p99_backwards_ns': float(np.percentile(-dt[dt < 0], 99)) if (dt < 0).any() else 0.0}
    E = {'w': np.arange(0, 1500, 10), 'wsm': np.arange(0, 3001, 10), 'chip': np.arange(-16 * 4096 * BC, 16 * 4096 * BC + 1, 4096 * BC / 8),
         'bcid': np.arange(0, 4097, 64), 'offs': np.arange(-16.5, 16), 'late': np.r_[-1e9, np.arange(-1000, 100001, 500), 1e12]}
    G = {}
    for lab, m in (('ot', ot), ('nb', (~S['flag']) & valid & (det > 0)), ('pad', S['flag'] & valid & (vmm <= 1))):
        G[f'h_w_{lab}'] = hist(S['w'][m], E['w'])
        G[f'h_wsm_{lab}'] = hist(S['wsm'][m], E['wsm'])
        G[f'h_chip_{lab}'] = hist(chip[m], E['chip'])
        G[f'h_bcid_{lab}'] = hist(S['bcid'][m], E['bcid'])
        G[f'h_offs_{lab}'] = hist(S['offs'][m], E['offs'])
        G[f'h_late_{lab}'] = hist(late[m], E['late'])
    # markers of the first file: period per VMM and agreement between VMMs of one marker round
    M = np.concatenate(mk_log) if mk_log else np.zeros((0, 3), np.int64)
    mk = {}
    for v in np.unique(M[:, 0]):
        s = M[M[:, 0] == v, 1]
        d = np.diff(s)
        mk[int(v)] = {'n': int(len(s)), 'period_ticks_median': float(np.median(d)) if len(d) else None,
                      'period_ticks_p05_p95': np.percentile(d, [5, 95]).tolist() if len(d) else None}
    # marker rounds: consecutive markers in the stream with the same timestamp
    if len(M):
        o = np.argsort(M[:, 2])
        Ms = M[o]
        new = np.r_[True, np.diff(Ms[:, 1]) != 0]
        rid = np.cumsum(new) - 1
        nper = np.bincount(rid)
        mk_rounds = {'n_rounds': int(len(nper)), 'vmms_per_round_median': float(np.median(nper)),
                     'frac_rounds_all_vmms': float(np.mean(nper >= len(mk))),
                     'round_step_ticks_median': float(np.median(np.diff(Ms[new, 1])))}
    else:
        mk_rounds = {}

    # ---- tracks (chamber-only straight tracks, as paddle_efficiency.py)
    o = np.argsort(np.where(valid, T, np.inf), kind='stable')
    o = o[valid[o]]
    counts = np.zeros((32, 64), int)
    np.add.at(counts, (vmm[ot], ch[ot]), 1)
    live = np.isfinite(X) & (counts > 0)
    for v in (int(x) for x in a.dead_vmm.split(',') if x):
        live[v] = False
    for v, c in parse_veto(a.veto):
        live[v, c] = False
    H = {}
    for d in ORDER:
        j = o[(det[o] == d) & S['flag'][o] & live[vmm[o], ch[o]]]
        H[d] = (T[j], X[vmm[j], ch[j]], Y[vmm[j], ch[j]], j)
    span = (T[o[-1]] - T[o[0]]) * 1e-9
    t2 = H[2][0]

    def near(d, tt):
        t = H[d][0]
        return np.searchsorted(t, tt - a.half_win), np.searchsorted(t, tt + a.half_win, side='right')
    se = {d: near(d, t2) for d in ORDER}
    cand = np.nonzero((se[3][1] > se[3][0]) & (se[1][1] > se[1][0]))[0]
    keep = np.concatenate([[True], np.diff(t2[cand]) > 300]) if len(cand) else np.zeros(0, bool)
    cand = cand[keep]
    rows, members = [], []
    for k in cand:
        cl, ok = {}, True
        for d in ORDER:
            s, e = se[d][0][k], se[d][1][k]
            x, y = H[d][1][s:e], H[d][2][s:e]
            if len(x) > a.max_tag or np.max(np.hypot(x - x.mean(), y - y.mean())) > a.cluster_r:
                ok = False
                break
            cl[d] = (x.mean(), y.mean(), H[d][3][s:e])
        if not ok:
            continue
        r = np.hypot(cl[2][0] - (cl[3][0] + cl[1][0]) / 2, cl[2][1] - (cl[3][1] + cl[1][1]) / 2)
        rows.append((t2[k], r, *(v for d in ORDER for v in cl[d][:2])))
        members.append([cl[d][2] for d in ORDER])
    R = np.array(rows)
    n = len(R)
    print(f'{len(cand)} three-chamber time groups, {n} clean, {int((R[:, 1] < a.tol).sum())} straight')

    P = {'t2': R[:, 0], 'r2': R[:, 1], 'xy': R[:, 2:].reshape(n, 3, 2)}
    P['wall'] = clock(P['t2'])
    P['bad'] = in_intervals(P['wall'], load_veto(a.bad_times)) if a.bad_times else np.zeros(n, bool)
    f3 = lambda dt: np.zeros((n, 3), dt)  # noqa: E731
    for key, dt in (('n', int), ('seq', np.int64), ('seq_min', np.int64), ('seq_max', np.int64), ('pkt', np.int64),
                    ('w', int), ('wsm', int), ('vmm', int), ('bcid', int), ('offs', int), ('srs', np.int64),
                    ('chip', float), ('t', float), ('nvmm', int)):
        P[key] = f3(dt)
    P['span'], P['pkt_span'] = np.zeros(n, np.int64), np.zeros(n, np.int64)
    P['foreign_ot'], P['foreign_nb'], P['foreign_pad'], P['foreign_ot_track_vmm'] = (np.zeros(n, np.int64) for _ in range(4))
    P['own_nb'] = np.zeros(n, np.int64)
    flag, isch = S['flag'], det > 0
    for i, mem in enumerate(members):
        allm = np.concatenate(mem)
        for c, idx in enumerate(mem):
            f = idx[np.argmin(T[idx])]
            P['n'][i, c] = len(idx)
            P['seq'][i, c], P['seq_min'][i, c], P['seq_max'][i, c] = f, idx.min(), idx.max()
            P['pkt'][i, c], P['w'][i, c], P['wsm'][i, c], P['vmm'][i, c] = S['pkt'][f], S['w'][f], S['wsm'][f], vmm[f]
            P['bcid'][i, c], P['offs'][i, c], P['srs'][i, c] = S['bcid'][f], S['offs'][f], S['srs'][f]
            P['chip'][i, c], P['t'][i, c] = chip[f], T[f]
            P['nvmm'][i, c] = len(np.unique(vmm[idx]))
        lo, hi = allm.min(), allm.max()
        P['span'][i], P['pkt_span'][i] = hi - lo, S['pkt'][hi] - S['pkt'][lo]
        inside = np.ones(hi - lo + 1, bool)
        inside[allm - lo] = False
        sl = slice(lo, hi + 1)
        fo = inside & flag[sl] & isch[sl]
        P['foreign_ot'][i] = fo.sum()
        P['foreign_ot_track_vmm'][i] = (fo & np.isin(vmm[sl], vmm[allm])).sum()
        # neighbour hits: of the track's own pads' VMMs within the window, or unrelated
        nb = inside & ~flag[sl] & isch[sl]
        tw = np.abs(T[sl] - P['t2'][i]) <= a.half_win + 200
        P['own_nb'][i] = (nb & tw & np.isin(vmm[sl], vmm[allm])).sum()
        P['foreign_nb'][i] = nb.sum() - P['own_nb'][i]
        P['foreign_pad'][i] = (inside & (vmm[sl] <= 1)).sum()
    P['order'] = np.argsort(P['seq'], axis=1)      # chamber columns in stream order (0 = det3, 1 = det2, 2 = det1)
    # chance reference: first hit of each chamber after t2 + sideband
    sb = np.zeros((n, 3), np.int64)
    sbt = np.zeros((n, 3))
    for c, d in enumerate(ORDER):
        t, j = H[d][0], H[d][3]
        q = np.clip(np.searchsorted(t, P['t2'] + a.sideband), 0, len(t) - 1)
        sb[:, c], sbt[:, c] = j[q], t[q]
    P['sb_seq'], P['sb_t'] = sb, sbt
    P['sb_span'] = sb.max(1) - sb.min(1)
    # paddles
    for spec in a.paddles.split(','):
        v, c = (int(s) for s in spec.split(':'))
        jp = o[(vmm[o] == v) & (ch[o] == c) & flag[o]]
        tp = T[jp]
        key = 'p' + spec.replace(':', '_')
        tt = P['t2'] - a.latency
        q = np.clip(np.searchsorted(tp, tt), 1, max(len(tp) - 1, 1))
        cand2 = np.c_[q - 1, np.minimum(q, len(tp) - 1)]
        dts = tp[cand2] - tt[:, None]
        pick = cand2[np.arange(n), np.argmin(np.abs(dts), 1)]
        dt = tp[pick] - tt
        hitp = np.abs(dt) <= a.win
        P[key + '_dt'] = np.where(hitp, dt, np.nan)
        P[key + '_dseq'] = np.where(hitp, jp[pick] - P['seq'][:, 1], 0)
        P[key + '_bcid'] = np.where(hitp, S['bcid'][jp[pick]], -1)
        P[key + '_srs'] = np.where(hitp, S['srs'][jp[pick]], -1)
        P[key + '_chip'] = np.where(hitp, chip[jp[pick]], np.nan)
        P[key + '_hit'] = hitp

    # ---- summary
    st = (P['r2'] < a.tol) & ~P['bad']
    def q(x):
        x = np.asarray(x)[st]
        return {'median': float(np.median(x)), 'p16_p84': np.percentile(x, [16, 84]).tolist(),
                'p05_p95': np.percentile(x, [5, 95]).tolist(), 'max': float(np.max(x))} if len(x) else None
    pairs = ((1, 0), (2, 0), (2, 1))            # det2-det3, det1-det3, det1-det2 (column indices)
    pn = {(1, 0): 'det2-det3', (2, 0): 'det1-det3', (2, 1): 'det1-det2'}
    summ = {'title': a.title, 'files': len(files), 'frames': cnt[0], 'words': cnt[1], 'hit_words': int(n_all),
            'hits_over_threshold': int(S['flag'].sum()), 'span_s': span, 'n_groups': int(len(cand)), 'n_clean': n,
            'n_straight_good': int(st.sum()), 'params': vars(a), 'markers_first_file': mk, 'marker_rounds': mk_rounds,
            'vmm_time_order': inv_vmm,
            'chamber_hits_lateness_ns': {'frac_late_gt_0': float(np.nanmean(late[ot] > 0)),
                                         'p50_p90_p99_p999': np.nanpercentile(late[ot], [50, 90, 99, 99.9]).tolist()},
            'tracks': {}}
    tr = summ['tracks']
    for c, d in enumerate(ORDER):
        tr[f'det{d}'] = {'word_in_frame': q(P['w'][:, c]), 'words_since_marker': q(P['wsm'][:, c]),
                         'offset': q(P['offs'][:, c]), 'bcid': q(P['bcid'][:, c]), 'chip_time_ns': q(P['chip'][:, c]),
                         'n_hits': q(P['n'][:, c])}
    for (i, j) in pairs:
        tr[pn[(i, j)]] = {'dt_ns': q(P['t'][:, i] - P['t'][:, j]),
                          'dbcid': q(((P['bcid'][:, i] - P['bcid'][:, j] + 2048) % 4096) - 2048),
                          'dbc_offs_bcid': q((P['offs'][:, i] * 4096 + P['bcid'][:, i]) - (P['offs'][:, j] * 4096 + P['bcid'][:, j])),
                          'dsrs_ticks': q(P['srs'][:, i] - P['srs'][:, j]),
                          'frac_same_srs': float(np.mean((P['srs'][:, i] == P['srs'][:, j])[st])),
                          'dseq': q(P['seq'][:, i] - P['seq'][:, j]), 'abs_dseq': q(np.abs(P['seq'][:, i] - P['seq'][:, j])),
                          'dpkt': q(P['pkt'][:, i] - P['pkt'][:, j]),
                          'frac_same_frame': float(np.mean((P['pkt'][:, i] == P['pkt'][:, j])[st]))}
    tr['span'] = q(P['span'])
    tr['frames_spanned'] = q(P['pkt_span'])
    tr['frac_all_in_one_frame'] = float(np.mean((P['pkt_span'] == 0)[st]))
    tr['foreign_ot'] = q(P['foreign_ot'])
    tr['foreign_neighbour'] = q(P['foreign_nb'])
    tr['own_neighbour'] = q(P['own_nb'])
    tr['foreign_paddle'] = q(P['foreign_pad'])
    tr['chance_span'] = q(P['sb_span'])
    oc = {}
    for row in P['order'][st]:
        key = ''.join(str(ORDER[x]) for x in row)
        oc[key] = oc.get(key, 0) + 1
    tr['stream_order_of_chambers'] = dict(sorted(oc.items(), key=lambda kv: -kv[1]))
    for spec in a.paddles.split(','):
        key = 'p' + spec.replace(':', '_')
        m = st & P[key + '_hit']
        tr[f'paddle {spec}'] = {'n': int(m.sum()),
                                'dseq_vs_det2': {'median': float(np.median(P[key + '_dseq'][m])) if m.any() else None,
                                                 'p05_p95': np.percentile(P[key + '_dseq'][m], [5, 95]).tolist() if m.any() else None},
                                'frac_same_srs_as_det2': float(np.mean((P[key + '_srs'] == P['srs'][:, 1])[m])) if m.any() else None,
                                'dt_det2_minus_paddle_ns': float(np.nanmedian(P['t'][m, 1] - (P['t2'][m] - a.latency + P[key + '_dt'][m]))) if m.any() else None}
    # stream excerpts: the hit words around a few straight tracks (markers in between counted by wsm)
    ex = []
    for i in np.nonzero(st)[0][:a.examples]:
        mem = set(np.concatenate(members[i]).tolist())
        lo, hi = int(P['seq_min'][i].min()) - 6, int(P['seq_max'][i].max()) + 6
        rows_ = []
        for j in range(max(lo, 0), min(hi, n_all - 1) + 1):
            rows_.append({'seq': j, 'frame': int(S['pkt'][j]), 'word': int(S['w'][j]), 'words_since_marker': int(S['wsm'][j]),
                          'vmm': int(vmm[j]), 'det': int(det[j]), 'ch': int(ch[j]), 'over_thr': bool(flag[j]),
                          'offset': int(S['offs'][j]), 'bcid': int(S['bcid'][j]), 'tdc': int(S['tdc'][j]),
                          'srs': int(S['srs'][j]), 't_minus_t2_ns': round(float(T[j] - P['t2'][i]), 1),
                          'track': j in mem})
        ex.append({'track': int(i), 'r2_mm': float(P['r2'][i]), 'hits': rows_})
    summ['examples'] = ex
    np.savez_compressed(a.out, **P, **G, **{f'edges_{k}': v for k, v in E.items()})
    json.dump(summ, open(a.out.replace('.npz', '') + '.json', 'w'), indent=1, default=float)
    print(json.dumps(tr, indent=1, default=float))


if __name__ == '__main__':
    sys.exit(main())
