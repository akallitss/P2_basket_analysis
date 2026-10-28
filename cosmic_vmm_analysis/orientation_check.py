#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Channel order and chamber orientation of the cosmic stand, without an external trigger.

Input: the .npz files of dump_candidates.py (small multi-chamber time groups). The CLEAN sample is the groups that
touch all three chambers with at most --max-hits hits.

Step 1 -- channel order inside each connector (as cosmic_bench_analysis/02_map_validation.py, which picks the
ordering whose occupancy looks like a real illumination). Per connector the candidate wirings are:
  * each VMM on a 64-strip half (as at SPS, or the two halves swapped) with an order inside the half among
    linear / reverse / pairswap / reverse+pairswap (16 x 2 options), or
  * the two VMMs on the odd / even strips of the whole connector, each linear or reverse (2 x 4 options).
Two scores, computed separately so they can be compared:
  A. adjacency: same-chip hit pairs of the clean sample (|dch| > 4, i.e. no electronic neighbours, ADC > --adc)
     whose pads are closer than 20 mm -- charge shared between neighbouring pads;
  B. smoothness: hits of the clean sample per pad compared with the median of the 6 nearest pads (the speckle
     metric of vmm_stations.map_roughness), minimised connector by connector.
Step 2 -- orientation between chambers (as cosmic_bench_analysis/03_m3_alignment.py): for events with one compact
cluster per chamber, rotate (and mirror) chamber B's positions about its own pad-cloud centre and scan the angle;
the figure of merit is r_xx + r_yy of the position correlation with chamber A. Lying flat, the right answer is a
rotation about the vertical axis by some angle, possibly mirrored (chamber turned upside down).

Usage: orientation_check.py OUT.json cand1.npz [cand2.npz ...] [--adc 150] [--max-hits 12] [--use A|B|sps]
"""
import argparse
import itertools
import json
import sys

import numpy as np
import pandas as pd

import vmm_qa_path                                          # noqa: F401
import vmm_stations as VS

DET_STATION = {3: 'P2_IN', 2: 'P2_MID', 1: 'P2_OUT'}
W = np.arange(64)
ORDERS = {'lin': W, 'rev': 63 - W, 'pswap': W ^ 1, 'rev-pswap': (63 - W) ^ 1}


def connector_options():
    """{name: (strips for the bottom VMM, strips for the top VMM)}, strips 0..127 within the connector."""
    opt = {}
    for (ob, sb), (ot, st) in itertools.product(ORDERS.items(), ORDERS.items()):
        opt[f'bot:{ob}|top:{ot}'] = (sb, 64 + st)
        opt[f'SWAP bot:{ob}|top:{ot}'] = (64 + sb, st)
    for kb in (0, 1):
        for rb, rt in itertools.product((False, True), repeat=2):
            b = 2 * W + kb
            t = 2 * W + (1 - kb)
            opt[f'IL bot:{"odd" if kb else "even"}{"-rev" if rb else ""}|top:{"rev" if rt else "fwd"}'] = (
                127 - b if rb else b, 127 - t if rt else t)
    return opt


OPTIONS = connector_options()
SPS = 'bot:rev|top:rev'


def parse_orders(specs):
    """['P2_IN:3=bot:rev|top:lin', ...] -> {('P2_IN', 3): 'bot:rev|top:lin'} (connector numbers as read out)."""
    out = {}
    for s in specs or []:
        k, name = s.split('=', 1)
        st, n = k.split(':')
        if name not in OPTIONS:
            raise SystemExit(f'unknown order {name!r}')
        out[(st, int(n))] = name
    return out


def load(files, max_hits):
    out = []
    off = 0
    for f in files:
        z = np.load(f)
        ev = z['ev'] + off
        off = int(ev.max()) + 1
        out.append(pd.DataFrame({'ev': ev, 'vmm': z['vmm'].astype(int), 'ch': z['ch'].astype(int),
                                 'adc': z['adc'].astype(int), 't': z['t'], 'det': z['det'].astype(int)}))
    h = pd.concat(out, ignore_index=True)
    g = h.groupby('ev')
    n = g['vmm'].transform('size')
    nch = g['det'].transform('nunique')
    return h[(nch == 3) & (n <= max_hits)].reset_index(drop=True)


class Mapper:
    def __init__(self, conn_shift=0):
        """conn_shift: added to the SPS connector number of every hybrid (-1 from run_112 on: c4-c6 -> c3-c5)."""
        self.pad = pd.read_csv(VS.MAP_CSV).set_index('channel_id').sort_index()
        self.halves = VS.vmm_half_map('bot')
        # connectors: (station, n) -> (bottom vmm, top vmm)
        self.conn = {}
        for v, (st, n, half) in self.halves.items():
            self.conn.setdefault((st, n + conn_shift), {})[half] = v
        self.choice = {k: SPS for k in self.conn}

    def xy_table(self, choice=None):
        """(32, 64) arrays of pad x, y for the current (or given) choice."""
        ch = {**self.choice, **(choice or {})}
        X = np.full((32, 64), np.nan)
        Y = np.full((32, 64), np.nan)
        for (st, n), vm in self.conn.items():
            sb, stp = OPTIONS[ch[(st, n)]]
            for v, strips in ((vm['bot'], sb), (vm['top'], stp)):
                r = self.pad.reindex((n - 1) * 128 + strips)
                X[v] = r['pad_cx'].to_numpy()
                Y[v] = r['pad_cy'].to_numpy()
        return X, Y


def same_chip_pairs(h, adc):
    s = h[h['adc'] > adc]
    pairs = []
    for (_, v), grp in s.groupby(['ev', 'vmm']):
        c = grp['ch'].to_numpy()
        for i in range(len(c)):
            for j in range(i + 1, len(c)):
                if abs(c[i] - c[j]) > 4:
                    pairs.append((v, c[i], c[j]))
    return np.array(pairs, dtype=int).reshape(-1, 3)


def adjacency(pairs, X, Y, vmms):
    m = np.isin(pairs[:, 0], vmms)
    p = pairs[m]
    if not len(p):
        return float('nan'), 0
    d = np.hypot(X[p[:, 0], p[:, 1]] - X[p[:, 0], p[:, 2]], Y[p[:, 0], p[:, 1]] - Y[p[:, 0], p[:, 2]])
    return float(np.nanmean(d < 20)), int(len(p))


def roughness(counts, X, Y, vmms):
    """speckle of the occupancy over the pads of the given vmms (same formula as vmm_stations.map_roughness)."""
    x = X[vmms].ravel()
    y = Y[vmms].ravel()
    n = counts[vmms].ravel().astype(float)
    ok = np.isfinite(x)
    x, y, n = x[ok], y[ok], n[ok]
    d = np.hypot(x[:, None] - x[None], y[:, None] - y[None])
    np.fill_diagonal(d, np.inf)
    nb = np.argsort(d, 1)[:, :6]
    med = np.median(n[nb], 1)
    live = med > 0.05 * med.max()
    return float(np.median(np.abs(n[live] - med[live]) / (n[live] + med[live] + 1)))


def step1(h, mp, adc):
    pairs = same_chip_pairs(h, adc)
    counts = np.zeros((32, 64))
    np.add.at(counts, (h['vmm'].to_numpy(), h['ch'].to_numpy()), 1)
    res = {'n_pairs': int(len(pairs)), 'connectors': {}}
    rng = np.random.default_rng(1)
    best_A, best_B = {}, {}
    for key, vm in sorted(mp.conn.items()):
        vmms = [vm['bot'], vm['top']]
        scA = {}
        for name in OPTIONS:
            X, Y = mp.xy_table({key: name})
            scA[name] = adjacency(pairs, X, Y, vmms)[0]
        X, Y = mp.xy_table()
        p = pairs[np.isin(pairs[:, 0], vmms)].copy()
        rnd = []
        for _ in range(20):
            q = p.copy()
            q[:, 2] = rng.permutation(q[:, 2])
            rnd.append(adjacency(q, X, Y, vmms)[0])
        best_A[key] = max(scA, key=lambda k: -1 if np.isnan(scA[k]) else scA[k])
        res['connectors'][f'{key[0]} c{key[1]}'] = {
            'vmms': vmms, 'n_pairs': adjacency(pairs, X, Y, vmms)[1], 'random_adjacency': float(np.nanmean(rnd)),
            'adjacency_sps': scA[SPS], 'adjacency_best': [best_A[key], scA[best_A[key]]],
            'adjacency_top5': sorted(((k, round(v, 3)) for k, v in scA.items() if not np.isnan(v)), key=lambda kv: -kv[1])[:5]}
    # smoothness: coordinate descent per station, starting from the SPS wiring
    for st in ('P2_IN', 'P2_MID', 'P2_OUT'):
        keys = [k for k in mp.conn if k[0] == st]
        vmms = sorted(v for k in keys for v in mp.conn[k].values())
        cur = {k: SPS for k in keys}
        X, Y = mp.xy_table(cur)
        r_sps = roughness(counts, X, Y, vmms)
        for _ in range(3):
            changed = False
            for k in keys:
                sc = {}
                for name in OPTIONS:
                    X, Y = mp.xy_table({**cur, k: name})
                    sc[name] = roughness(counts, X, Y, vmms)
                b = min(sc, key=sc.get)
                if b != cur[k] and sc[b] < sc[cur[k]] - 1e-4:
                    cur[k] = b
                    changed = True
            if not changed:
                break
        X, Y = mp.xy_table(cur)
        # random reference: shuffle the counts among the station's channels
        cnt = counts.copy()
        flat = cnt[vmms].ravel()
        rr = []
        for _ in range(10):
            cnt[vmms] = rng.permutation(flat).reshape(len(vmms), 64)
            rr.append(roughness(cnt, X, Y, vmms))
        res[f'smoothness_{st}'] = {'roughness_sps': r_sps, 'roughness_best': roughness(counts, X, Y, vmms),
                                   'roughness_random': float(np.mean(rr)),
                                   'best': {f'c{k[1]}': cur[k] for k in keys}}
        best_B.update(cur)
    return res, best_A, best_B


def clusters(h, X, Y, max_spread=30.0):
    """one position per (event, chamber) when the mapped pads are within max_spread mm of their mean."""
    d = h.assign(x=X[h['vmm'], h['ch']], y=Y[h['vmm'], h['ch']]).dropna(subset=['x'])
    g = d.groupby(['ev', 'det'])
    d = d.assign(mx=g['x'].transform('mean'), my=g['y'].transform('mean'))
    d = d.assign(far=np.hypot(d['x'] - d['mx'], d['y'] - d['my']) > max_spread)
    c = d.groupby(['ev', 'det']).agg(x=('mx', 'first'), y=('my', 'first'), far=('far', 'any')).reset_index()
    c = c[~c['far']]
    wide = c.pivot(index='ev', columns='det', values=['x', 'y']).dropna()
    return wide


def rotation_scan(xa, ya, xb, yb, dtheta=2.0):
    xb0, yb0 = xb - xb.mean(), yb - yb.mean()
    out = {}
    best = None
    for mirror in (1, -1):
        th = np.arange(0, 360, dtheta)
        fom = []
        for t in np.radians(th):
            yy = mirror * yb0
            xr = np.cos(t) * xb0 - np.sin(t) * yy
            yr = np.sin(t) * xb0 + np.cos(t) * yy
            rxx = np.corrcoef(xa, xr)[0, 1]
            ryy = np.corrcoef(ya, yr)[0, 1]
            fom.append((rxx, ryy, rxx + ryy))
        fom = np.array(fom)
        out[mirror] = (th, fom)
        j = int(np.argmax(fom[:, 2]))
        if best is None or fom[j, 2] > best[3]:
            best = (mirror, float(th[j]), float(fom[j, 0]), float(fom[j, 2]), float(fom[j, 1]))
    return out, best


def step2(h, X, Y, n_shuffle=20):
    w = clusters(h, X, Y)
    res = {'n_events': int(len(w))}
    rng = np.random.default_rng(2)
    for a, b in ((3, 2), (2, 1), (3, 1)):
        xa, ya, xb, yb = (w[('x', a)].to_numpy(), w[('y', a)].to_numpy(), w[('x', b)].to_numpy(), w[('y', b)].to_numpy())
        scan, best = rotation_scan(xa, ya, xb, yb)
        null = []
        for _ in range(n_shuffle):
            p = rng.permutation(len(xb))
            null.append(rotation_scan(xa, ya, xb[p], yb[p], dtheta=6.0)[1][3])
        th1, f1 = scan[1]
        res[f'det{a}-det{b}'] = {'best_mirror': best[0], 'best_theta_deg': best[1], 'best_fom': best[3],
                                 'fom_at_identity': float(f1[0, 2]), 'null_fom_mean': float(np.mean(null)),
                                 'null_fom_max': float(np.max(null)),
                                 'fom_curve_mirror+1': f1[:, 2].round(3).tolist(), 'fom_curve_mirror-1': scan[-1][1][:, 2].round(3).tolist()}
        print(f'det{a}-det{b}: {len(xa)} events, best mirror {best[0]:+d} theta {best[1]:.0f} deg FOM {best[3]:.3f} '
              f'(identity {f1[0, 2]:.3f}; shuffled max {np.max(null):.3f})')
    return res


def single_hit_events(h):
    """events with exactly one hit per chamber: arrays (n, 3) of vmm and ch, columns det3, det2, det1."""
    g = h.groupby(['ev', 'det'])['vmm'].transform('size')
    s = h[g == 1]
    full = s.groupby('ev')['det'].transform('size') == 3
    s = s[full].sort_values(['ev', 'det'], ascending=[True, False])
    return s['vmm'].to_numpy().reshape(-1, 3), s['ch'].to_numpy().reshape(-1, 3), s['ev'].to_numpy().reshape(-1, 3)[:, 0]


def midpoint_score(V, C, X, Y, tol):
    """fraction of events whose det2 pad is within tol mm of the midpoint of the det3 and det1 pads, after the
    median offset of det2 (chamber shift) is removed."""
    x, y = X[V, C], Y[V, C]
    ok = np.isfinite(x).all(1)
    rx = x[ok, 1] - 0.5 * (x[ok, 0] + x[ok, 2])
    ry = y[ok, 1] - 0.5 * (y[ok, 0] + y[ok, 2])
    near = np.hypot(rx, ry) < 80
    sx, sy = (np.median(rx[near]), np.median(ry[near])) if near.sum() > 20 else (0.0, 0.0)
    return float(np.mean(np.hypot(rx - sx, ry - sy) < tol)), (float(sx), float(sy)), int(ok.sum())


def step1_tracks(h, mp, tol=30.0, sweeps=4):
    """step 1c: choose the wiring of all nine connectors that makes the most events straight lines
    (det2 at the midpoint of det3 and det1), by coordinate descent from the SPS wiring; validated on the other
    half of the events (fit on even, test on odd event ids) and against shuffled events."""
    V, C, E = single_hit_events(h)
    fit, test = (E % 2 == 0), (E % 2 == 1)
    keys = sorted(mp.conn)
    cur = {k: SPS for k in keys}
    hist = []
    for sweep in range(sweeps):
        changed = False
        for k in keys:
            sc = {}
            for name in OPTIONS:
                X, Y = mp.xy_table({**cur, k: name})
                sc[name] = midpoint_score(V[fit], C[fit], X, Y, tol)[0]
            b = max(sc, key=sc.get)
            if b != cur[k] and sc[b] > sc[cur[k]] + 1e-4:
                cur[k] = b
                changed = True
        X, Y = mp.xy_table(cur)
        hist.append(midpoint_score(V[fit], C[fit], X, Y, tol)[0])
        if not changed:
            break
    rng = np.random.default_rng(3)
    out = {'n_events': int(len(V)), 'tol_mm': tol, 'sweeps': hist}
    for label, ch in (('sps', {k: SPS for k in keys}), ('fit', cur)):
        X, Y = mp.xy_table(ch)
        s_fit = midpoint_score(V[fit], C[fit], X, Y, tol)
        s_test = midpoint_score(V[test], C[test], X, Y, tol)
        sh = []
        for _ in range(10):
            Vs, Cs = V[test].copy(), C[test].copy()
            p = rng.permutation(len(Vs))
            Vs[:, 1], Cs[:, 1] = Vs[p, 1], Cs[p, 1]
            sh.append(midpoint_score(Vs, Cs, X, Y, tol)[0])
        out[label] = {'score_fit_half': s_fit[0], 'score_test_half': s_test[0], 'det2_shift_mm': s_test[1],
                      'score_shuffled_test': float(np.mean(sh)), 'n_test': s_test[2]}
    out['wiring'] = {f'{k[0]} c{k[1]}': v for k, v in cur.items()}
    return out, cur


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('out_json')
    ap.add_argument('files', nargs='+')
    ap.add_argument('--adc', type=int, default=150)
    ap.add_argument('--max-hits', type=int, default=12)
    ap.add_argument('--use', default='B', choices=['A', 'B', 'sps'])
    a = ap.parse_args()
    h = load(a.files, a.max_hits)
    print(f'clean sample: {h["ev"].nunique()} events, {len(h)} hits')
    mp = Mapper()
    res1, bA, bB = step1(h, mp, a.adc)
    for k, r in res1['connectors'].items():
        print(f'{k:10s} VMM {r["vmms"]}: pairs {r["n_pairs"]:4d}  random {r["random_adjacency"]:.2f}  SPS {r["adjacency_sps"]:.2f}  '
              f'best {r["adjacency_best"][0]} {r["adjacency_best"][1]:.2f}')
    for st in ('P2_IN', 'P2_MID', 'P2_OUT'):
        r = res1[f'smoothness_{st}']
        print(f'{st}: roughness SPS {r["roughness_sps"]:.3f}, best {r["roughness_best"]:.3f}, random {r["roughness_random"]:.3f}: {r["best"]}')
    agree = {f'{k[0]} c{k[1]}': bA[k] == bB[k] for k in bA}
    print('A and B agree on:', [k for k, v in agree.items() if v])
    out = {'params': vars(a), 'step1': res1, 'choice_A': {f'{k[0]} c{k[1]}': v for k, v in bA.items()},
           'choice_B': {f'{k[0]} c{k[1]}': v for k, v in bB.items()}}
    r1c, bC = step1_tracks(h, mp)
    out['step1_tracks'] = r1c
    print(f'\nstep 1c (straight tracks): {r1c["n_events"]} single-hit events, sweeps {np.round(r1c["sweeps"], 3)}')
    for label in ('sps', 'fit'):
        r = r1c[label]
        print(f'  {label}: fit half {r["score_fit_half"]:.3f}, test half {r["score_test_half"]:.3f}, '
              f'shuffled {r["score_shuffled_test"]:.3f}, det2 shift {np.round(r["det2_shift_mm"], 1)} mm')
    print('  wiring:', r1c['wiring'])
    for label, choice in (('sps', {}), ('A', bA), ('B', bB), ('tracks', bC)):
        print(f'\nstep 2 with the {label} wiring:')
        X, Y = mp.xy_table(choice)
        out[f'step2_{label}'] = step2(h, X, Y)
    json.dump(out, open(a.out_json, 'w'), indent=1, default=str)


if __name__ == '__main__':
    sys.exit(main())
