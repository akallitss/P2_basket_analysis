#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""bcid_stream.py for the SPS H4 VMM data (July/August 2026): BCIDs and read-out-stream positions of the
trigger-referenced three-station tracks.

Same FEC, same firmware and the same cabling as the cosmic stand (vmm_stations.py): P2_IN = VMM 2-7, P2_MID = VMM
8-13, P2_OUT = VMM 14-19, the beam trigger (scintillators) on VMM 0 channel 44. The `offset` field counts 0..15 here
as in the lab, the markers come every 65536 ticks in rounds of all 20 VMMs, so the hit time is the same as in
bcid_stream.py:  t = srs * 22.5 ns + (offset*4096 + bcid + 1) * 22.5 ns - tdc * 60/255 ns.

Track: a trigger hit with no other trigger within --iso ns, and in each of the three stations a clean cluster
(<= --max-tag over-threshold hits on live mapped pads, all within --cluster-r mm of their mean) within +- --half-win
ns of trigger + that station's latency (the peak of the station - trigger time difference, fitted on the first
capture). Straight: P2_MID within --tol mm of the P2_IN-P2_OUT midpoint (the stations are equally spaced in z,
320 / 630 / 940 mm), after subtracting the median offset of that residual (relative alignment). Hot channels
(> --hot-ratio x the larger of the median channel of their VMM and of their station, first capture) are not used for clusters but stay in the stream.

The output has the keys of bcid_stream.py so that bcid_stream_plots.py --setup sps reads it: chamber columns
0 / 1 / 2 = P2_IN / P2_MID / P2_OUT (the slots of det3 / det2 / det1), the trigger is "paddle" p0_44 (dt = 0 by
construction), t2 = trigger time. Captures are read one at a time (about 33 M hit words each); seq and pkt run on
across captures, the stream state (last marker per VMM) is carried over.

Usage (lxplus, ~/cva3 with bcid_stream.py and its imports, VMM_QA_DIR = the vmm_qa dir with vmm_stations.py):
  bcid_stream_sps.py RAW_DIR OUT.npz [--max-files N] [--title run_46]
"""
import argparse
import glob
import json
import os
import sys
import time

import numpy as np

from bcid_stream import read_stream, hist

sys.path.insert(0, os.environ.get('VMM_QA_DIR', os.path.expanduser('~/tpc/vmm_qa')))
import vmm_stations as VS  # noqa: E402

BC = 22.5
STATION = np.zeros(32, int)
STATION[2:8], STATION[8:14], STATION[14:20] = 1, 2, 3        # column + 1: P2_IN, P2_MID, P2_OUT
NAMES = ('P2_IN', 'P2_MID', 'P2_OUT')


def pad_xy():
    tab = VS.build_pad_table()
    X, Y = np.full((32, 64), np.nan), np.full((32, 64), np.nan)
    m = tab['mapped'].values
    X[tab['vmm'].values[m], tab['ch'].values[m]] = tab['pad_cx'].values[m]
    Y[tab['vmm'].values[m], tab['ch'].values[m]] = tab['pad_cy'].values[m]
    return X, Y


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('raw')
    ap.add_argument('out')
    ap.add_argument('--trigger', default='0:44')
    ap.add_argument('--iso', type=float, default=1000.0, help='no other trigger within +- this [ns]')
    ap.add_argument('--half-win', type=float, default=150.0)
    ap.add_argument('--max-tag', type=int, default=6)
    ap.add_argument('--cluster-r', type=float, default=15.0)
    ap.add_argument('--tol', type=float, default=25.0)
    ap.add_argument('--hot-ratio', type=float, default=8.0)
    ap.add_argument('--sideband', type=float, default=5000.0)
    ap.add_argument('--max-files', type=int, default=0)
    ap.add_argument('--examples', type=int, default=12)
    ap.add_argument('--title', default='')
    a = ap.parse_args()

    files = sorted(glob.glob(os.path.join(a.raw, '*.pcapng')))
    if a.max_files:
        files = files[:a.max_files]
    tv, tc = (int(s) for s in a.trigger.split(':'))
    X, Y = pad_xy()
    E = {'w': np.arange(0, 1500, 10), 'wsm': np.arange(0, 3001, 10),
         'chip': np.arange(-16 * 4096 * BC, 16 * 4096 * BC + 1, 4096 * BC / 8),
         'bcid': np.arange(0, 4097, 64), 'offs': np.arange(-16.5, 16), 'late': np.r_[-1e9, np.arange(-1000, 100001, 500), 1e12]}
    G = {}
    state, cnt = {}, [0, 0]
    seq0 = 0
    rows = []                      # per track dicts of scalars / 3-vectors
    ex = []
    lat, live, hot = None, None, {}
    lat_hist = {}
    t_first = t_last = None
    late_ot_all = []
    inv = {}
    n_trig = n_iso = n_cand = 0
    tstart = time.time()
    for fi, p in enumerate(files):
        r = read_stream(p, state, cnt)
        S = {k: np.concatenate([x[k] for x in r]) for k in r[0]}
        del r
        n = len(S['vmm'])
        vmm, ch, flag = S['vmm'].astype(int), S['ch'].astype(int), S['flag']
        chip = (S['offs'].astype(np.float64) * 4096 + S['bcid'] + 1) * BC - S['tdc'] * 60.0 / 255.0
        T = S['srs'] * BC + chip
        valid = S['srs'] > 0
        st_of = STATION[vmm]
        ot = flag & valid & (st_of > 0)
        trg = flag & valid & (vmm == tv) & (ch == tc)
        tvv = T[valid]
        if len(tvv):
            t_first = tvv.min() if t_first is None else min(t_first, tvv.min())
            t_last = tvv.max() if t_last is None else max(t_last, tvv.max())
        # global stream properties
        tvm = np.where(valid, T, -np.inf)
        runmax = np.maximum.accumulate(tvm)
        late = np.where(valid, np.r_[0.0, runmax[:-1]] - T, np.nan)
        del tvm, runmax
        for lab, m in (('ot', ot), ('nb', (~flag) & valid & (st_of > 0)), ('pad', flag & valid & (vmm <= 1))):
            for key, x in (('w', S['w']), ('wsm', S['wsm']), ('chip', chip), ('bcid', S['bcid']), ('offs', S['offs']), ('late', late)):
                G[f'h_{key}_{lab}'] = G.get(f'h_{key}_{lab}', 0) + hist(x[m], E[key])
        if fi < 4:
            late_ot_all.append(late[ot][::7])
            for v in np.unique(vmm[ot]):
                mm = np.nonzero(ot & (vmm == v))[0]
                d = np.diff(T[mm])
                q_ = inv.setdefault(int(v), [0, 0, []])
                q_[0] += len(d)
                q_[1] += int((d < 0).sum())
                q_[2].append(-d[d < 0][:20000])
        # first capture: hot channels, live pads, latency per station
        if live is None:
            counts = np.zeros((32, 64), int)
            np.add.at(counts, (vmm[ot], ch[ot]), 1)
            live = np.isfinite(X) & (counts > 0)
            st_med = {}
            for s_ in (1, 2, 3):
                c = counts[STATION == s_].ravel()
                st_med[s_] = np.median(c[c > 0]) if (c > 0).any() else 0
            for v in range(2, 20):
                c = counts[v]
                # reference: the larger of the VMM's and the station's median channel (a VMM at the edge of the
                # beam has a low median, and its beam-spot pads must not look hot)
                med = max(np.median(c[c > 0]) if (c > 0).any() else 0, st_med[STATION[v]])
                h_ = np.nonzero(c > a.hot_ratio * max(med, 1))[0]
                if len(h_):
                    hot[v] = h_.tolist()
                    live[v, h_] = False
        o = np.argsort(np.where(valid, T, np.inf), kind='stable')
        o = o[valid[o]]
        jt = o[trg[o]]
        tt = T[jt]
        n_trig += len(jt)
        dprev = np.r_[np.inf, np.diff(tt)]
        dnext = np.r_[np.diff(tt), np.inf]
        iso = (dprev > a.iso) & (dnext > a.iso)
        n_iso += int(iso.sum())
        H = {}
        for c in range(3):
            j = o[(st_of[o] == c + 1) & flag[o] & live[vmm[o], ch[o]]]
            H[c] = (T[j], X[vmm[j], ch[j]], Y[vmm[j], ch[j]], j)
        if lat is None:
            lat = []
            for c in range(3):
                t = H[c][0]
                q = np.clip(np.searchsorted(tt, t), 1, len(tt) - 1)
                d = np.where(np.abs(t - tt[q - 1]) < np.abs(t - tt[q]), t - tt[q - 1], t - tt[q])
                h_, e_ = np.histogram(d, bins=np.arange(-2000, 3001, 10))
                lat_hist[NAMES[c]] = h_.tolist()
                k = np.argmax(h_)
                lat.append(float(e_[k] + 5))
            lat_hist['edges'] = e_.tolist()
            print('latency per station [ns]', lat, 'hot', hot, flush=True)
        se = []
        for c in range(3):
            t = H[c][0]
            se.append((np.searchsorted(t, tt + lat[c] - a.half_win), np.searchsorted(t, tt + lat[c] + a.half_win, side='right')))
        cand = np.nonzero(iso & np.all([s[1] > s[0] for s in se], axis=0))[0]
        n_cand += len(cand)
        isch = st_of > 0
        it_all = jt
        sb_idx = []
        for c in range(3):
            t = H[c][0]
            sb_idx.append(np.clip(np.searchsorted(t, tt + a.sideband), 0, len(t) - 1))
        for k in cand:
            cl, ok = [], True
            for c in range(3):
                s, e = se[c][0][k], se[c][1][k]
                x, y = H[c][1][s:e], H[c][2][s:e]
                if len(x) > a.max_tag or np.max(np.hypot(x - x.mean(), y - y.mean())) > a.cluster_r:
                    ok = False
                    break
                cl.append((x.mean(), y.mean(), H[c][3][s:e]))
            if not ok:
                continue
            mem = [z[2] for z in cl]
            allm = np.concatenate(mem)
            R = {'t2': tt[k], 'd2': np.array([cl[1][0] - (cl[0][0] + cl[2][0]) / 2, cl[1][1] - (cl[0][1] + cl[2][1]) / 2]),
                 'xy': np.array([[z[0], z[1]] for z in cl])}
            for key in ('n', 'seq', 'seq_min', 'seq_max', 'pkt', 'w', 'wsm', 'vmm', 'bcid', 'offs', 'srs', 'chip', 't', 'nvmm'):
                R[key] = []
            for idx in mem:
                f = idx[np.argmin(T[idx])]
                R['n'].append(len(idx)); R['seq'].append(seq0 + f); R['seq_min'].append(seq0 + idx.min())
                R['seq_max'].append(seq0 + idx.max()); R['pkt'].append(S['pkt'][f]); R['w'].append(S['w'][f])
                R['wsm'].append(S['wsm'][f]); R['vmm'].append(vmm[f]); R['bcid'].append(S['bcid'][f])
                R['offs'].append(S['offs'][f]); R['srs'].append(S['srs'][f]); R['chip'].append(chip[f])
                R['t'].append(T[f]); R['nvmm'].append(len(np.unique(vmm[idx])))
            lo, hi = allm.min(), allm.max()
            R['span'], R['pkt_span'] = hi - lo, S['pkt'][hi] - S['pkt'][lo]
            inside = np.ones(hi - lo + 1, bool)
            inside[allm - lo] = False
            sl = slice(lo, hi + 1)
            fo = inside & flag[sl] & isch[sl]
            R['foreign_ot'] = int(fo.sum())
            R['foreign_ot_track_vmm'] = int((fo & np.isin(vmm[sl], vmm[allm])).sum())
            nb = inside & ~flag[sl] & isch[sl]
            tw = np.abs(T[sl] - tt[k] - lat[1]) <= a.half_win + 200
            R['own_nb'] = int((nb & tw & np.isin(vmm[sl], vmm[allm])).sum())
            R['foreign_nb'] = int(nb.sum()) - R['own_nb']
            R['foreign_pad'] = int((inside & (vmm[sl] <= 1)).sum())
            R['foreign_hot'] = int((inside & isch[sl] & ~live[vmm[sl], ch[sl]]).sum())
            R['sb_seq'] = [seq0 + H[c][3][sb_idx[c][k]] for c in range(3)]
            R['sb_t'] = [H[c][0][sb_idx[c][k]] for c in range(3)]
            f_t = it_all[k]
            R['p0_44_dt'], R['p0_44_dseq'] = 0.0, int(f_t - (R['seq'][1] - seq0))
            R['p0_44_dseq_first'] = int(f_t - (min(R['seq']) - seq0))
            R['p0_44_bcid'], R['p0_44_srs'], R['p0_44_chip'], R['p0_44_hit'] = int(S['bcid'][f_t]), int(S['srs'][f_t]), float(chip[f_t]), True
            R['p0_44_seq'] = seq0 + f_t
            rows.append(R)
            if len(ex) < a.examples and fi == 0 and len(rows) % 50 == 1:
                memset = set(allm.tolist()) | {int(f_t)}
                lo2, hi2 = int(min(lo, f_t)) - 6, int(max(hi, f_t)) + 6
                ex.append({'track': len(rows) - 1, 'hits': [
                    {'seq': int(seq0 + j), 'frame': int(S['pkt'][j]), 'word': int(S['w'][j]), 'words_since_marker': int(S['wsm'][j]),
                     'vmm': int(vmm[j]), 'station': NAMES[st_of[j] - 1] if st_of[j] else ('trigger' if (vmm[j] == tv and ch[j] == tc) else f'vmm{vmm[j]}'),
                     'ch': int(ch[j]), 'over_thr': bool(flag[j]), 'offset': int(S['offs'][j]), 'bcid': int(S['bcid'][j]),
                     'tdc': int(S['tdc'][j]), 'srs': int(S['srs'][j]), 't_minus_trigger_ns': round(float(T[j] - tt[k]), 1),
                     'track': j in memset} for j in range(max(lo2, 0), min(hi2, n - 1) + 1)]})
        seq0 += n
        print(f'{fi + 1}/{len(files)} {os.path.basename(p)}: {n} hit words, {len(jt)} triggers, {int(iso.sum())} isolated, '
              f'{len(cand)} three-station, {len(rows)} clean tracks so far, {time.time() - tstart:.0f} s', flush=True)
        del S, T, chip, late, o, H
    # assemble
    nt = len(rows)
    P = {}
    for key in rows[0]:
        P[key] = np.array([r_[key] for r_ in rows])
    off = np.median(P['d2'], axis=0)
    P['r2'] = np.hypot(*(P['d2'] - off).T)
    P['r2_raw'] = np.hypot(*P['d2'].T)
    P['bad'] = np.zeros(nt, bool)
    P['order'] = np.argsort(P['seq'], axis=1)
    P['sb_span'] = P['sb_seq'].max(1) - P['sb_seq'].min(1)
    st = P['r2'] < a.tol
    span = (t_last - t_first) * 1e-9
    lo_all = np.concatenate(late_ot_all)

    def q(x):
        x = np.asarray(x, float)[st]
        return {'median': float(np.median(x)), 'p16_p84': np.percentile(x, [16, 84]).tolist(),
                'p05_p95': np.percentile(x, [5, 95]).tolist(), 'max': float(np.max(x))} if len(x) else None
    summ = {'title': a.title, 'setup': 'sps', 'files': len(files), 'frames': cnt[0], 'words': cnt[1], 'hit_words': int(seq0),
            'span_s': span, 'n_triggers': n_trig, 'n_isolated_triggers': n_iso, 'n_groups': n_cand, 'n_clean': nt,
            'n_straight_good': int(st.sum()), 'params': vars(a), 'latency_ns': dict(zip(NAMES, lat)), 'latency_hist': lat_hist,
            'hot_channels': hot, 'alignment_offset_mm': off.tolist(),
            'r2_p50_p90': np.percentile(P['r2'], [50, 90]).tolist(),
            'vmm_time_order': {v: {'n': x[0], 'frac_backwards': x[1] / max(x[0], 1),
                                   'p99_backwards_ns': float(np.percentile(np.concatenate(x[2]), 99)) if x[1] else 0.0}
                               for v, x in inv.items()},
            'chamber_hits_lateness_ns': {'frac_late_gt_0': float(np.nanmean(lo_all > 0)),
                                         'p50_p90_p99_p999': np.nanpercentile(lo_all, [50, 90, 99, 99.9]).tolist()},
            'marker_rounds': {}, 'tracks': {}}
    tr = summ['tracks']
    tr['span'] = q(P['span'])
    tr['frac_all_in_one_frame'] = float(np.mean((P['pkt_span'] == 0)[st]))
    tr['foreign_hot'] = q(P['foreign_hot'])
    summ['examples'] = ex
    np.savez_compressed(a.out, **P, **G, **{f'edges_{k}': v for k, v in E.items()})
    json.dump(summ, open(a.out.replace('.npz', '') + '.json', 'w'), indent=1, default=float)
    print(json.dumps({k: summ[k] for k in ('n_triggers', 'n_isolated_triggers', 'n_groups', 'n_clean', 'n_straight_good',
                                           'latency_ns', 'alignment_offset_mm', 'r2_p50_p90', 'tracks')}, indent=1, default=float))


if __name__ == '__main__':
    sys.exit(main())
