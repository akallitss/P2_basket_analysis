#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Absolute efficiency of single VMM channels on the cosmic bench, without a model of the tracking resolution.

THE MEASUREMENT (details, checks and how-to: CHANNEL_EFFICIENCY.md)
  Tags: a muon seen by two chambers predicts a point in the third (the probe); only tags with a scintillator coincidence
  are used (self-made tags contain correlated noise). The prediction is wrong by ~3 mm (probe in the middle) to ~8 mm
  (probe outside), comparable to the 12 mm pads, so "did the predicted pad fire" measures the tracking, not the pad.
  Counting instead: the error moves muons, it does not create or remove them. For pad j

      fired_j = eps_j x (muons that crossed pad j) + other hits,   muons that crossed = density x area_j,

  the density counted in pad cells around the pad. K(R) = firings of the pads from tags predicted within R (random
  hits from an off-time window removed) = eps E + b R^2 once R >= R_lo > the prediction error; b from the increments
  beyond R_lo, eps = (K(R_lo) - b R_lo^2) / E. No resolution model (the predictions sit on the pad lattice).

SUBCOMMANDS (run each with -h for all options)
  dump       decode the raw captures of runs: efficiency_tp.py --dump-tags per run (per-tag records + json)
  sdt-table  global threshold (sdt) of every chip from the ext files of the runs (run on the DAQ box)
  thr-table  threshold and baseline of every channel in every run, in mV (monitor map + sdt calibration)
  measure    the counting estimator: per run, per (chamber, threshold - baseline) bin, per channel; figures
  toy        closure test on simulated muons: is the true efficiency recovered for a given resolution?

Example chain (wk2 scan, every path real):
  python3 channel_efficiency.py dump --runs-dir /eos/experiment/ntof/data/x17/p2_sps_july/vmm_cosmic_run_oct26/runs \\
      --out ~/wk2tags --jobs 3 181 184 187
  python3 channel_efficiency.py thr-table TAGS_DIR --map map_Km18_i5_unmasked.json --sdtcal sdtcal.json \\
      --run-sdt run_sdt.json -o thr.csv
  python3 channel_efficiency.py measure TAGS_DIR thr.csv OUT_DIR --scint 140:310
"""
import argparse
import glob
import json
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DETV = {3: range(2, 8), 2: range(8, 14), 1: range(14, 20)}
NAME = {3: 'P2_IN', 2: 'P2_MID', 1: 'P2_OUT'}
PNUM = {v: k for k, v in NAME.items()}
DEFAULT_ORDERS = ['P2_IN:3=bot:rev|top:lin', 'P2_IN:5=bot:rev|top:lin']


# ----------------------------------------------------------------------------------------------------------------- core
RADII = np.arange(4, 81, 4.0)     # counting radii (mm); the dump keeps fired pads within 80 mm
R_DENS = 40.0                     # radius (mm) of the circle the local muon density is counted in
PAD_WH = (12.05, 11.6)            # median P2 pad size (mm): tangential x radial


def pos_key(x, y):
    """integer key of a pad-centre position (0.1 mm), to match fired pads to pads."""
    return np.round(np.asarray(x) * 10).astype(np.int64) * 100000 + np.round(np.asarray(y) * 10).astype(np.int64)


def _fired_sets(fi, fx, fy):
    out = {}
    for i, kk in zip(np.asarray(fi).tolist(), pos_key(fx, fy).tolist()):
        out.setdefault(kk, set()).add(i)
    return out


def count_curves(tags, pad_x, pad_y, radii=RADII, r_dens=R_DENS, live=None):
    """For every pad j: n[j] = local tag density x pi r_dens^2, and for every radius R: k[j, R] = tags predicted within R
    for which pad j fired in the muon window, s[j, R] the same in the sideband (random) window, v[j, R] = variance of
    k - s. Tags are weighted +1 (on-time) / -r_scale (off-time pairs).
    Density, in pad cells: every tag belongs to the cell of the live pad nearest to its predicted point; the density
    around pad j = tags in the cells whose centres are within r_dens of pad j / the area of those cells (live = (x, y,
    area) of every live pad of the chamber). The cells tile the plane, so this is exact whether the predicted points
    are spread continuously or sit on the pad lattice (single-pad tags), and it follows the edges and masked pads.
    Without live: tags within r_dens / (pi r_dens^2)."""
    from scipy.spatial import cKDTree
    px, py = tags['px'], tags['py']
    w = np.where(tags['on'], 1.0, -float(tags['r_scale']))
    tree = cKDTree(np.c_[px, py])
    F = _fired_sets(tags['fired_i'], tags['fired_x'], tags['fired_y'])
    S = _fired_sets(tags['side_i'], tags['side_x'], tags['side_y']) if 'side_i' in tags else {}
    pk = pos_key(pad_x, pad_y)
    rmax = max(radii[-1], r_dens)
    if live is not None:
        lx, ly, la = (np.asarray(q, float) for q in live)
        lt = cKDTree(np.c_[lx, ly])
        cell = lt.query(np.c_[px, py])[1]                         # each tag's cell (nearest live pad)
        wcell = np.bincount(cell, w, len(lx))                     # weighted tags per cell
        nbr = [np.asarray(lt.query_ball_point([x, y], r_dens), int) for x, y in zip(pad_x, pad_y)]
        n_cells = np.array([wcell[q].sum() for q in nbr])
        dens_area = np.array([la[q].sum() for q in nbr])
    else:
        n_cells, dens_area = None, np.full(len(pad_x), np.pi * r_dens ** 2)
    n = np.zeros(len(pad_x))
    k, s, v = (np.zeros((len(pad_x), len(radii))) for _ in range(3))
    for j, (x, y) in enumerate(zip(pad_x, pad_y)):
        idx = np.asarray(tree.query_ball_point([x, y], rmax), int)
        if not len(idx):
            continue
        d = np.hypot(px[idx] - x, py[idx] - y)
        n[j] = (n_cells[j] if n_cells is not None else w[idx[d < r_dens]].sum()) * np.pi * r_dens ** 2 / dens_area[j]
        f = np.isin(idx, list(F.get(int(pk[j]), ())))
        g = np.isin(idx, list(S.get(int(pk[j]), ())))
        for ir, r in enumerate(radii):
            q = d < r
            k[j, ir] = (w[idx[q]] * f[q]).sum()
            s[j, ir] = (w[idx[q]] * g[q]).sum()
            v[j, ir] = (w[idx[q]] ** 2 * (f[q] | g[q])).sum()
    return n, k, s, v


RLO_SCAN = (24.0, 32.0, 40.0, 48.0, 56.0)


def fit_group(n, k, s, v, area, radii=RADII, r_dens=R_DENS, rlo=40.0):
    """eps of a group of pads, no resolution model.
    E = sum_j n_j area_j / (pi r_dens^2) = expected muon crossings of the pads (density from the tags within r_dens).
    K(R) = sum_j (k - s): pad firings from tags predicted within R, random (sideband) hits removed.
    Once rlo is larger than the prediction error, all the pads' own muons are inside and K(R) = eps E + b R^2 for
    R >= rlo, b R^2 being firings correlated with muons that crossed elsewhere. b from the INCREMENTS of K beyond rlo
    (independent counts), then eps = (K(rlo) - b rlo^2) / E; errors from the counts (Poisson, off-time weights)."""
    E = float((n * np.asarray(area) / (np.pi * r_dens ** 2)).sum())
    K = (k - s).sum(0)
    V = v.sum(0)
    out = dict(E=E, K=K, Kerr=np.sqrt(np.maximum(V, 0)), eps=np.nan, err=np.nan, b=np.nan, chi2ndf=np.nan)
    if E <= 0:
        return out
    i0 = int(np.argmin(np.abs(radii - rlo)))
    dK = np.diff(K[i0:])
    dV = np.maximum(np.diff(V[i0:]), 1.0)
    dR2 = np.diff(radii[i0:] ** 2)
    if len(dK):
        wts = 1 / dV
        b = float((wts * dK * dR2).sum() / (wts * dR2 ** 2).sum())       # weighted LS through the origin
        vb = float(1 / (wts * dR2 ** 2).sum())
        chi2 = float((wts * (dK - b * dR2) ** 2).sum())
        vb *= max(1.0, chi2 / max(len(dK) - 1, 1))
    else:
        b, vb, chi2 = 0.0, 0.0, np.nan
    R0 = radii[i0]
    eps = (K[i0] - b * R0 ** 2) / E
    err = np.sqrt(max(V[i0], 1.0) + (R0 ** 2) ** 2 * vb) / E
    out.update(eps=float(eps), err=float(err), b=b / E, chi2ndf=chi2 / max(len(dK) - 1, 1), rlo=float(R0))
    # stability against the choice of rlo (must be flat once rlo exceeds the prediction error)
    for r in RLO_SCAN:
        if r <= radii[-2]:
            j = int(np.argmin(np.abs(radii - r)))
            out[f'eps_rlo{int(r)}'] = float((K[j] - b * radii[j] ** 2) / E)
    return out


def select_tags(tags, m):
    """keep the tags in mask m; the fired-pad lists (fired_*, side_*) are re-indexed."""
    new = np.full(len(m), -1)
    new[m] = np.arange(m.sum())
    out = {k: (v[m] if np.ndim(v) and len(v) == len(m) else v) for k, v in tags.items()
           if not k.startswith(('fired_', 'side_'))}
    for pre in ('fired', 'side'):
        if f'{pre}_i' in tags:
            keep = m[tags[f'{pre}_i']]
            out[f'{pre}_i'] = new[tags[f'{pre}_i'][keep]]
            out[f'{pre}_x'], out[f'{pre}_y'] = tags[f'{pre}_x'][keep], tags[f'{pre}_y'][keep]
    return out


def load_tags(npz, p):
    Z = np.load(npz)
    return {k.split('_', 1)[1]: Z[k] for k in Z.files if k.startswith(f'p{p}_')}


def run_info(tags_dir, run):
    """title, ext-file tag, K and dead/masked pads of a run, from its efficiency_tp.py json."""
    E = json.load(open(f'{tags_dir}/eff_run{run}.json'))
    tag = E['title'].split(' ', 1)[1] if ' ' in E['title'] else ''
    m = re.search(r'_K([pm]?\d+)', tag)
    K = int(m[1].replace('p', '+').replace('m', '-')) if m else None
    dead = {tuple(map(int, s.split(':'))) for lst in E['dead_or_masked'].values() for s in lst}
    return dict(title=E['title'], tag=tag, K=K, dead=dead, params=E['params'])


def pad_geometry(conn_shift=-1, orders=None):
    """(32, 64) arrays x, y, area of the pad of each read-out channel, with the map of the efficiency runs."""
    sys.path.insert(0, HERE)
    import orientation_check as O
    import vmm_stations as VS
    from scipy.spatial import cKDTree
    X, Y = O.Mapper(conn_shift).xy_table(O.parse_orders(orders or DEFAULT_ORDERS))
    p = pd.read_csv(VS.MAP_CSV)
    area = np.full((32, 64), np.nan)
    ok = np.isfinite(X)
    d, i = cKDTree(np.c_[p.pad_cx, p.pad_cy]).query(np.c_[X[ok], Y[ok]])
    assert d.max() < 0.01, f'pad table mismatch {d.max()}'
    area[ok] = p.pad_area.to_numpy()[i]
    return X, Y, area


# ------------------------------------------------------------------------------------------------------------- commands
def cmd_dump(a):
    """efficiency_tp.py --dump-tags on every run; skips runs already done."""
    os.makedirs(a.out, exist_ok=True)

    def one(run):
        if os.path.exists(f'{a.out}/tags_run{run}.npz') and not a.force:
            return f'run_{run}: exists, skipped'
        d = sorted(glob.glob(f'{a.runs_dir}/run_{run}/cfg_*'))
        if not d:
            return f'run_{run}: no cfg_* directory under {a.runs_dir}/run_{run}'
        files = sorted(glob.glob(f'{d[0]}/raw_daq_data/*.pcapng'))
        lab = os.path.basename(d[0]).replace('cfg_gain4.5_peaktime200_opt_cosmic_', '')
        cmd = [sys.executable, os.path.join(HERE, 'efficiency_tp.py'), ','.join(files), f'{a.out}/eff_run{run}.json',
               f'--conn-shift={a.conn_shift}', '--dead-vmm', a.dead_vmm, '--dump-tags', f'{a.out}/tags_run{run}.npz',
               '--scint-ch', a.scint_ch, '--title', f'run_{run} {lab}']
        if a.diag:
            cmd.append('--dump-diag')
        for o in a.order or DEFAULT_ORDERS:
            cmd += ['--order', o]
        if a.mem:
            cmd = ['systemd-run', '--user', '--scope', '-q', '-p', f'MemoryMax={a.mem}'] + cmd
        with open(f'{a.out}/eff_run{run}.log', 'w') as log:
            rc = subprocess.call(cmd, stdout=log, stderr=subprocess.STDOUT)
        return f'run_{run}: rc={rc} ({len(files)} files)'
    with ThreadPoolExecutor(a.jobs) as ex:
        for msg in ex.map(one, a.runs):
            print(msg, flush=True)


def cmd_sdt_table(a):
    """global threshold of every chip (sdt = 64 * bank1spi2 + bank1spi1) and masked channels of each ext file."""
    out = {}
    for f in a.files:
        lo, hi, sm = {}, {}, {}
        for line in open(os.path.join(a.ext_dir, f)):
            m = re.match(r'\s*Sys \d+ Bec \d+ Hyb (\d+) Vmm (\d+) Reg (\S+)\s+Field (\S+)\s+(\d+)', line)
            if not m:
                continue
            chip, reg, fld, val = 2 * int(m[1]) + int(m[2]), m[3], m[4], int(m[5])
            if fld == 'sdt' and reg == 'bank1spi1':
                lo[chip] = val
            if fld == 'sdt' and reg == 'bank1spi2':
                hi[chip] = val
            if fld == 'sm' and val == 1:
                sm.setdefault(chip, []).append(int(re.search(r'(\d+)', reg)[1]))
        out[f] = {'sdt': {c: 64 * hi.get(c, 0) + lo.get(c, 0) for c in sorted(set(lo) | set(hi))}, 'masked': sm}
    json.dump(out, open(a.o, 'w'), indent=1) if a.o != '-' else print(json.dumps(out))


def cmd_thr_table(a):
    """one row per (run, vmm, ch): threshold, baseline, distance in mV."""
    runmap = dict(x.split('=', 1) for x in a.run_map)
    base = json.load(open(a.map)) if a.map else None
    cal = json.load(open(a.sdtcal))['fit'] if a.sdtcal else None
    rsdt = json.load(open(a.run_sdt)) if a.run_sdt else {}
    rows = []
    for fz in sorted(glob.glob(f'{a.tags_dir}/tags_run*.npz')):
        run = re.search(r'run(\d+)', fz)[1]
        info = run_info(a.tags_dir, run)
        if run in runmap:      # a map measured with this run's own configuration (e.g. trimmed): taken as is
            mp, shift = json.load(open(runmap[run])), None
        else:                  # the reference map, shifted to this run's global threshold
            mp = base
            sdt = rsdt[f'{a.ext_prefix}{info["tag"]}.txt']['sdt']
            shift = {int(v): (sdt[str(v)] - mp['sdt'][str(v)]) * cal[str(v)]['mV_per_sdt'] for v in mp['sdt']}
        for v, chans in mp['baseline'].items():
            for c, bl in chans.items():
                thr = mp['threshold'][v][c] + (shift[int(v)] if shift else 0.0)
                rows.append((int(run), info['K'], int(v), int(c), thr, bl, thr - bl))
    T = pd.DataFrame(rows, columns=['run', 'K', 'vmm', 'ch', 'threshold_mV', 'baseline_mV', 'dist_mV'])
    T.to_csv(a.o, index=False)
    print(f'{len(T)} rows, {T.run.nunique()} runs -> {a.o}')


def real_run(pid):
    """real run number of a pseudo-run id (run x 100 + chunk, or run x 1000 + chunk for >= 100 chunks); a real run
    number (< 1000) is returned as is."""
    return pid if pid < 1000 else (pid // 1000 if pid >= 100000 else pid // 100)


def bad_time_mask(tags, run, bad_times):
    """True for the tags inside the HV/LV bad-time intervals of spark_veto.py (BAD_TIMES/veto_run_<run>.json);
    needs the wall times of a --dump-diag dump. Missing veto file or wall times: an error, never silently no veto."""
    from spark_veto import load_veto, in_intervals
    f = f'{bad_times}/veto_run_{real_run(run)}.json'
    if not os.path.exists(f):
        raise SystemExit(f'--bad-times: no {f}')
    if 'wall' not in tags:
        raise SystemExit(f'--bad-times: run {run} has no wall times (dump it with --dump-diag)')
    return in_intervals(tags['wall'], load_veto(f))


def hull_edge_distance(pts, q):
    """distance (mm) of the points q to the boundary of the convex hull of the pad centres pts (>= 0 inside)."""
    from scipy.spatial import ConvexHull
    hv = pts[ConvexHull(pts).vertices]
    a, b = hv, np.roll(hv, -1, axis=0)
    ab = b - a
    t = np.clip(((q[:, None, :] - a[None]) * ab[None]).sum(-1) / (ab ** 2).sum(-1)[None], 0, 1)
    d = np.hypot(*(q[:, None, :] - (a[None] + t[..., None] * ab[None])).transpose(2, 0, 1))
    return d.min(1)


def parse_fiducial(spec):
    """'P2_IN:20,P2_OUT:20' -> {'P2_IN': 20.0, 'P2_OUT': 20.0}"""
    out = {}
    for item in filter(None, (spec or '').split(',')):
        det, mm = item.split(':')
        if det not in NAME.values():
            raise SystemExit(f'--fiducial: unknown chamber {det} (one of {sorted(NAME.values())})')
        out[det] = float(mm)
    return out


def measure_pads(tags_dir, thr, scint=None, clean=40.0, radii=RADII, conn_shift=-1, orders=None, verbose=True,
                 bad_times=None, veto_log=None, fiducial=None):
    """per (run, chamber, live pad): the count curves (n, k, s, v) and the pad's threshold - baseline distance;
    plus the tag - scintillator times of all on-time tags (to check the window). bad_times: directory of
    spark_veto.py json files; tags inside those intervals are dropped (counts per run and chamber in veto_log).
    fiducial: {chamber: mm}; a pad whose centre is closer than mm to the outline (convex hull of all its chamber's pad
    centres) counts as not clean, so it leaves every result (measure, regroup.py, checks.py) like a pad near a dead one;
    d_edge (mm) is kept per pad either way."""
    sys.path.insert(0, HERE)
    from scipy.spatial import cKDTree
    X, Y, area = pad_geometry(conn_shift, orders)
    thr = thr.set_index(['run', 'vmm', 'ch'])
    out, dts = [], []
    for fz in sorted(glob.glob(f'{tags_dir}/tags_run*.npz')):
        run = int(re.search(r'run(\d+)', fz)[1])
        info = run_info(tags_dir, run)
        for p in (3, 2, 1):
            tags = load_tags(fz, p)
            dts.append(tags['dt_scint'][tags['on']])
            if bad_times:
                bad = bad_time_mask(tags, run, bad_times)
                if veto_log is not None:
                    sc = (tags['dt_scint'] >= scint[0]) & (tags['dt_scint'] <= scint[1]) if scint else True
                    veto_log.append(dict(run=run, det=NAME[p], tags=int(len(bad)), vetoed=int(bad.sum()),
                                         on_scint=int((tags['on'] & sc).sum()),
                                         on_scint_vetoed=int((tags['on'] & sc & bad).sum())))
                tags = select_tags(tags, ~bad)
            if scint:
                tags = select_tags(tags, (tags['dt_scint'] >= scint[0]) & (tags['dt_scint'] <= scint[1]))
            allp = [(v, c) for v in DETV[p] for c in range(64) if np.isfinite(X[v, c])]
            live = [(v, c) for v, c in allp if (v, c) not in info['dead']]
            deadxy = np.array([(X[v, c], Y[v, c]) for v, c in allp if (v, c) in info['dead']]).reshape(-1, 2)
            px, py = np.array([X[v, c] for v, c in live]), np.array([Y[v, c] for v, c in live])
            ok_clean = (cKDTree(deadxy).query(np.c_[px, py])[0] > clean) if (clean and len(deadxy)) else \
                np.ones(len(live), bool)
            allxy = np.array([(X[v, c], Y[v, c]) for v, c in allp])
            d_edge = hull_edge_distance(allxy, np.c_[px, py])
            ok_fid = d_edge >= (fiducial or {}).get(NAME[p], -1.0)
            n, k, s, v = count_curves(tags, px, py, radii, live=(px, py, [area[q] for q in live]))
            for j, (vv, cc) in enumerate(live):
                d = thr.loc[(run, vv, cc)] if (run, vv, cc) in thr.index else None
                out.append(dict(run=run, K=info['K'], det=NAME[p], vmm=vv, ch=cc, x=px[j], y=py[j], area=area[vv, cc],
                                clean=bool(ok_clean[j] and ok_fid[j]), d_edge=float(d_edge[j]), fid=bool(ok_fid[j]),
                                dist_mV=float(d.dist_mV) if d is not None else np.nan,
                                n=n[j], k=k[j], s=s[j], v=v[j]))
        if verbose:
            print(f'run_{run} counted', flush=True)
    return out, np.concatenate(dts)


def group_eps(rows, radii=RADII, rlo=40.0):
    a = {key: np.array([r[key] for r in rows]) for key in ('n', 'k', 's', 'v', 'area')}
    return fit_group(a['n'], a['k'], a['s'], a['v'], a['area'], radii, rlo=rlo)


def cmd_measure(a):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    sys.path.insert(0, HERE)
    from event_display import style, INK, MUTED, COL, MRK
    os.makedirs(a.out_dir, exist_ok=True)
    thr = pd.read_csv(a.thr_table)
    scint = tuple(map(float, a.scint.split(':'))) if a.scint else None
    vlog = []
    fid = parse_fiducial(a.fiducial)
    rows, dts = measure_pads(a.tags_dir, thr, scint, a.clean, bad_times=a.bad_times or None, veto_log=vlog,
                             fiducial=fid)
    if fid:
        for det, mm in fid.items():
            rr = [r for r in rows if r['det'] == det]
            print(f'--fiducial {det}: pads closer than {mm:g} mm to the outline dropped: '
                  f'{sum(not r["fid"] for r in rr)} of {len(rr)} (run x pad)')
    if a.bad_times:
        V = pd.DataFrame(vlog)
        V.to_csv(f'{a.out_dir}/bad_times_vetoed.csv', index=False)
        t = V.groupby('det')[['on_scint', 'on_scint_vetoed']].sum()
        print('--bad-times', a.bad_times, ': scintillator-confirmed on-time tags vetoed per chamber:',
              {d: f'{r.on_scint_vetoed}/{r.on_scint} ({100 * r.on_scint_vetoed / max(r.on_scint, 1):.1f} %)'
               for d, r in t.iterrows()})
    h, e = np.histogram(dts[np.isfinite(dts)], np.arange(-1000, 1001, 25))
    print('tag time - scintillator time (on-time tags), 25 ns bins, non-empty:',
          {int(e[i]): int(h[i]) for i in np.nonzero(h)[0]})
    use = [r for r in rows if r['clean'] and r['dist_mV'] >= 0]
    blind = [r for r in rows if r['clean'] and r['dist_mV'] < 0]
    res, curves = [], {}

    def add(kind, key, sel, keep_curve=False):
        if not len(sel):
            return
        g = group_eps(sel, rlo=a.rlo)
        rec = dict(kind=kind, **key, dist_mV=float(np.mean([r['dist_mV'] for r in sel])), pads=len(sel),
                   **{k: g.get(k, np.nan) for k in ['eps', 'err', 'b', 'E', 'chi2ndf'] +
                      [f'eps_rlo{int(r)}' for r in RLO_SCAN]})
        res.append(rec)
        if keep_curve:
            curves[(kind, key.get('det'))] = g
    for (run, det) in sorted({(r['run'], r['det']) for r in use}):
        add('run', dict(run=run, K=next(r['K'] for r in use if r['run'] == run), det=det),
            [r for r in use if r['run'] == run and r['det'] == det])
    bins = np.array(a.bins)
    for det in NAME.values():
        for lo, hi in zip(bins[:-1], bins[1:]):
            add('distance', dict(det=det, bin=f'{lo:g}-{hi:g}'),
                [r for r in use if r['det'] == det and lo <= r['dist_mV'] < hi])
        add('plateau', dict(det=det), [r for r in use if r['det'] == det and a.plateau[0] <= r['dist_mV'] < a.plateau[1]],
            keep_curve=True)
    add('blind', dict(det='all'), blind)
    if a.per_channel:
        for (vv, cc) in sorted({(r['vmm'], r['ch']) for r in use}):
            rr = [r for r in use if r['vmm'] == vv and r['ch'] == cc]
            for lo, hi in zip(bins[:-1], bins[1:]):
                add('channel', dict(det=rr[0]['det'], vmm=vv, ch=cc, bin=f'{lo:g}-{hi:g}'),
                    [r for r in rr if lo <= r['dist_mV'] < hi])
    R = pd.DataFrame(res)
    R.to_csv(f'{a.out_dir}/efficiency.csv', index=False)
    pd.DataFrame([{k: v for k, v in r.items() if k not in ('n', 'k', 's', 'v')} for r in rows]).to_csv(
        f'{a.out_dir}/pads.csv', index=False)
    np.savez_compressed(f'{a.out_dir}/curves.npz', radii=RADII, r_dens=R_DENS,
                        **{key: np.array([r[key] for r in rows]) for key in ('n', 'k', 's', 'v', 'area', 'run', 'vmm',
                                                                             'ch', 'dist_mV', 'clean')},
                        det=np.array([r['det'] for r in rows]))
    pd.set_option('display.width', 220)
    print(f'\nplateau ({a.plateau[0]:g}-{a.plateau[1]:g} mV), R_lo = {a.rlo:g} mm; eps_rloX = the same with R_lo = X '
          f'(must agree once X > the prediction error):')
    print(R[R.kind == 'plateau'][['det', 'eps', 'err', 'b', 'chi2ndf', 'pads'] +
                                 [f'eps_rlo{int(r)}' for r in RLO_SCAN]].round(4).to_string(index=False))
    print('blind pads (baseline above threshold):',
          R[R.kind == 'blind'][['eps', 'err', 'pads']].round(4).to_dict('records'))
    dd = R[R.kind == 'distance']
    print(dd.pivot(index='bin', columns='det', values='eps').round(3))
    print(dd.pivot(index='bin', columns='det', values='err').round(3))

    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for det, p in PNUM.items():
        g = curves.get(('plateau', det))
        if g is None or not np.isfinite(g['eps']):
            continue
        ax[0].errorbar(RADII, g['K'] / g['E'], g['Kerr'] / g['E'], marker=MRK[p], color=COL[p], ms=4, lw=1, elinewidth=1,
                       label=f'{det}: eps {g["eps"]:.2f} +- {g["err"]:.2f}')
        rr = np.linspace(0, RADII[-1], 50)
        ax[0].plot(rr, g['eps'] + g['b'] * rr ** 2, color=COL[p], lw=0.9, ls='--')
        d = dd[(dd.det == det) & dd.eps.notna()]
        ax[1].errorbar(d.dist_mV, d.eps, d.err, marker=MRK[p], color=COL[p], ms=4, lw=1.5, label=det)
    ax[0].set_xlabel('R (mm): tags counted within R of the pad centre', fontsize=8, color=INK)
    ax[0].set_ylabel('pad firings (random removed) / expected crossings', fontsize=8, color=INK)
    ax[0].axvline(a.rlo, color=MUTED, lw=0.8, ls=':')
    ax[0].set_title(f'plateau pads ({a.plateau[0]:g}-{a.plateau[1]:g} mV); dashed: eps + b R$^2$ from R >= {a.rlo:g} mm',
                    fontsize=9, color=INK, loc='left')
    ax[1].set_ylim(0, 1.1)
    ax[1].set_xlabel('threshold - baseline of the channel (mV)', fontsize=8, color=INK)
    ax[1].set_ylabel('single-channel efficiency', fontsize=8, color=INK)
    ax[1].set_title('absolute efficiency' + (' (scintillator-confirmed tags)' if scint else ''), fontsize=9, color=INK,
                    loc='left')
    for x in ax:
        x.legend(frameon=False, fontsize=7)
        style(x)
    fig.tight_layout()
    fig.savefig(f'{a.out_dir}/efficiency.png', dpi=150)
    print(f'\nwrote {a.out_dir}/efficiency.csv, pads.csv, curves.npz, efficiency.png')


def toy(sigma, eps=0.85, n_mu=40000, pitch=(12.0, 11.6), half=12, noise=0.0, tail=(0.0, 20.0), impure=0.0, seed=1,
        radii=RADII, rlo=40.0, quantize=False):
    """closure test. Muons uniform over a (2 half + 1)^2 grid of pads; a crossed pad fires with probability eps;
    prediction = truth + Gauss(sigma) per coordinate (a fraction tail[0] with tail[1] instead); noise = probability
    per tag that any pad fires at random (in the muon AND the sideband window, independently); impure = fraction of
    fake tags (random position, nothing of theirs fires). The central pads are measured as the real ones are.
    Returns the fitted eps (and sigma) and the naive 'did the nearest pad fire' efficiency."""
    rng = np.random.default_rng(seed)
    W, H = pitch
    gx, gy = np.meshgrid(np.arange(-half, half + 1) * W, np.arange(-half, half + 1) * H, indexing='ij')
    pad_x, pad_y = gx.ravel(), gy.ravel()
    L = (half + 0.5) * np.array([W, H])
    tx, ty = rng.uniform(-L[0], L[0], n_mu), rng.uniform(-L[1], L[1], n_mu)
    crossed = (np.round(tx / W).astype(int) + half) * (2 * half + 1) + np.round(ty / H).astype(int) + half
    s = np.where(rng.random(n_mu) < tail[0], tail[1], sigma)
    px, py = tx + rng.normal(0, 1, n_mu) * s, ty + rng.normal(0, 1, n_mu) * s
    if quantize:   # predictions on the pad lattice, as when the tag clusters are single pads (pad centres)
        px, py = np.round(px / W) * W, np.round(py / H) * H
    hit = rng.random(n_mu) < eps
    fi = list(np.nonzero(hit)[0])
    fx, fy = list(pad_x[crossed[hit]]), list(pad_y[crossed[hit]])
    si, sx, sy = [], [], []
    if noise:
        for lst in ((fi, fx, fy), (si, sx, sy)):
            t, j = np.nonzero(rng.random((n_mu, len(pad_x))) < noise)
            lst[0].extend(t); lst[1].extend(pad_x[j]); lst[2].extend(pad_y[j])
    n_fake = int(impure * n_mu / (1 - impure)) if impure else 0
    tags = dict(px=np.r_[px, rng.uniform(-L[0], L[0], n_fake)], py=np.r_[py, rng.uniform(-L[1], L[1], n_fake)],
                on=np.ones(n_mu + n_fake, bool), r_scale=np.array(0.0), fired_i=np.array(fi, int), fired_x=np.array(fx),
                fired_y=np.array(fy), side_i=np.array(si, int), side_x=np.array(sx), side_y=np.array(sy))
    central = (np.abs(pad_x) <= (half - 7) * W) & (np.abs(pad_y) <= (half - 7) * H)
    n, k, sd, v = count_curves(tags, pad_x[central], pad_y[central], radii,
                               live=(pad_x, pad_y, np.full(len(pad_x), W * H)))
    g = fit_group(n, k, sd, v, np.full(central.sum(), W * H), radii, rlo=rlo)
    nx, ny = np.round(px / W).astype(int) + half, np.round(py / H).astype(int) + half
    ok = (nx >= 0) & (nx <= 2 * half) & (ny >= 0) & (ny <= 2 * half)
    naive = float(np.mean(hit[ok] & (nx[ok] * (2 * half + 1) + ny[ok] == crossed[ok])))
    return dict(sigma=sigma, eps_true=eps, eps_est=g['eps'], err=g['err'], b=g['b'], chi2ndf=g['chi2ndf'], naive=naive,
                **{k: v for k, v in g.items() if k.startswith('eps_rlo')})


def cmd_toy(a):
    print(f'{"sigma":>6} {"tail":>9} {"noise":>7} {"impure":>6}   true   estimate (R_lo {a.rlo:g})       R_lo scan '
          f'{RLO_SCAN}        naive own-pad')
    for sg in a.sigmas:
        r = toy(sg, a.eps, a.muons, noise=a.noise, tail=(a.tail_frac, a.tail_sigma), impure=a.impure, seed=a.seed,
                rlo=a.rlo, quantize=a.quantize)
        pull = (r['eps_est'] - a.eps) / r['err']
        print(f'{sg:6.1f} {a.tail_frac:4.2f}x{a.tail_sigma:<4.0f} {a.noise:7.4f} {a.impure:6.2f}   {a.eps:.3f}  '
              f'{r["eps_est"]:.3f} +- {r["err"]:.3f} (pull {pull:+.1f})   '
              + ' '.join(f'{r.get(f"eps_rlo{int(x)}", np.nan):.3f}' for x in RLO_SCAN) + f'   {r["naive"]:.3f}')


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0], formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog=__doc__.split('\n', 1)[1])
    sp = ap.add_subparsers(dest='cmd', required=True)
    d = sp.add_parser('dump', help='decode runs into per-tag records (efficiency_tp.py --dump-tags)')
    d.add_argument('runs', nargs='+', help='run numbers')
    d.add_argument('--runs-dir', required=True, help='directory holding run_N/cfg_*/raw_daq_data/*.pcapng')
    d.add_argument('--out', required=True)
    d.add_argument('--jobs', type=int, default=1, help='runs decoded in parallel (each needs a few GB of RAM)')
    d.add_argument('--mem', default='10G', help="memory cap per run via systemd-run ('' = none)")
    d.add_argument('--conn-shift', type=int, default=-1)
    d.add_argument('--order', action='append', help=f'connector orders (default {DEFAULT_ORDERS})')
    d.add_argument('--dead-vmm', default='0,1')
    d.add_argument('--scint-ch', default='0:40')
    d.add_argument('--force', action='store_true')
    d.add_argument('--diag', action='store_true', help='efficiency_tp.py --dump-diag: wall time of every tag (needed by '
                   'measure --bad-times), probe hits near the prediction with time and ADC, chamber activity')
    d.set_defaults(f=cmd_dump)
    s = sp.add_parser('sdt-table', help='chip thresholds (sdt) and masks of ext files')
    s.add_argument('ext_dir')
    s.add_argument('files', nargs='+')
    s.add_argument('-o', default='-', help='output json (- = stdout)')
    s.set_defaults(f=cmd_sdt_table)
    t = sp.add_parser('thr-table', help='threshold, baseline, distance (mV) per run and channel')
    t.add_argument('tags_dir')
    t.add_argument('--map', help='monitor map json (baseline + threshold per channel) of the reference ext file')
    t.add_argument('--sdtcal', help='sdtcal.json (mV per sdt count per chip)')
    t.add_argument('--run-sdt', help='sdt-table output for the runs')
    t.add_argument('--run-map', action='append', default=[], help='RUN=map.json: a map measured with that run\'s '
                   'own configuration (used as is)')
    t.add_argument('--ext-prefix', default='p2b-config-cern-ext_gain4.5_peaktime200_opt_cosmic_')
    t.add_argument('-o', required=True)
    t.set_defaults(f=cmd_thr_table)
    m = sp.add_parser('measure', help='the counting estimator')
    m.add_argument('tags_dir')
    m.add_argument('thr_table')
    m.add_argument('out_dir')
    m.add_argument('--scint', default='140:310', help="LO:HI ns window of (tag - scintillator) time; '' = all tags")
    m.add_argument('--rlo', type=float, default=40.0, help='radius (mm) from which all of a pad\'s own muons are '
                   'counted (> 4 x the prediction error); check the printed eps_rloX columns agree')
    m.add_argument('--per-channel', action='store_true', help='also one fit per channel and distance bin (slow, '
                   'low statistics)')
    m.add_argument('--clean', type=float, default=40, help='drop pads with a dead/masked pad within this distance')
    m.add_argument('--fiducial', default='', help="CHAMBER:MM[,...], e.g. 'P2_IN:20,P2_OUT:20': also drop the pads of "
                   'that chamber whose centre is closer than MM to the outline of its pads (default: none)')
    m.add_argument('--bins', type=float, nargs='+', default=[0, 14, 18, 22, 26, 30, 34, 40, 46, 52, 58, 64, 70, 100])
    m.add_argument('--plateau', type=float, nargs=2, default=[18, 40], help='distance range of the plateau summary')
    m.add_argument('--bad-times', default='', help='directory of spark_veto.py veto_run_<run>.json: drop the tags in '
                   'those HV/LV bad-time intervals (needs a --dump-diag dump, which has the wall time of every tag)')
    m.set_defaults(f=cmd_measure)
    o = sp.add_parser('toy', help='closure test on simulated muons')
    o.add_argument('--sigmas', type=float, nargs='+', default=[3, 8, 12])
    o.add_argument('--eps', type=float, default=0.85)
    o.add_argument('--muons', type=int, default=40000)
    o.add_argument('--noise', type=float, default=0.0)
    o.add_argument('--tail-frac', type=float, default=0.0)
    o.add_argument('--tail-sigma', type=float, default=20.0)
    o.add_argument('--impure', type=float, default=0.0)
    o.add_argument('--seed', type=int, default=1)
    o.add_argument('--rlo', type=float, default=40.0)
    o.add_argument('--quantize', action='store_true', help='predictions snapped to the pad lattice (single-pad tags)')
    o.set_defaults(f=cmd_toy)
    a = ap.parse_args()
    a.f(a)


if __name__ == '__main__':
    main()
