#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
11_hv_scan_efficiency.py

P2 detector efficiency (and resolution) vs mesh HV, adapted from
nTof_x17/mx_june_cosmic_qa/10_hv_scan_efficiency.py.

The mesh-HV scan run (p2_det1_mesh_hv_scan_7-2-26) has one sub_run per HV point,
named mesh_<NNN>V_drift_<MMM>V (mesh 345..420 V in 5 V steps, drift = mesh+180).
Each sub_run holds combined_hits_root / m3_tracking_root / hv_monitor.csv but NO
run_config -> the det1 long-run wiring is borrowed (cfg.run_config_ref).

Method (mirrors 06_efficiency_maps.py, one detector geometry for all HV points):
  1. Build the channel table once (borrowed run_config + pad map, dead conn dropped).
  2. Load each sub_run's M3 rays + P2 pad centroids (optional per-sub_run HV
     spark veto using that sub_run's hv_monitor.csv).
  3. Fit the pad->M3 rigid transform ONCE on the pooled matched events across all
     HV points (the mounting geometry does not change with HV), and freeze the
     transformed pad footprint -> a single active area used for every HV point,
     so the efficiency denominator region is identical across the scan.
  4. Per HV point: reco = transform(pad centroid); for every clean M3 single
     track in the active area ask whether a reco is within R mm (numerator) and
     whether P2 fired any pad (has_any). A track with no P2 event is a genuine
     miss kept in the denominator.
  5. Efficiency vs HV and core spatial resolution (robust sigma of the aligned
     residual) vs HV.

Products (<Analysis>/<detN>/<scan_run>/hv_scan/11_hv_scan_efficiency/):
  efficiency_vs_hv<suffix>.png / resolution_vs_hv<suffix>.png
  efficiency_vs_hv<suffix>.csv

Usage:
  python3 11_hv_scan_efficiency.py [det1_hvscan] [--r 20] [--active-r 30]
        [--z 246] [--fit-fiducial 300] [--min-valid 50] [--no-veto-sparks]
"""

import os
import re
import glob
import argparse
import matplotlib
matplotlib.use('Agg')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree

import p2_qa_config as qa
import p2_mapping as pmap
import p2_align as pa
import p2_sparks as ps


def find_subruns(cfg):
    """Discover mesh_<NNN>V sub_runs under the scan run, ascending in mesh HV."""
    run_dir = cfg.run_dir
    out = []
    for name in sorted(os.listdir(run_dir)):
        m = re.search(r'mesh_(\d+)V', name)
        if not m:
            continue
        hits = glob.glob(os.path.join(run_dir, name, 'combined_hits_root', '*.root'))
        m3 = glob.glob(os.path.join(run_dir, name, 'm3_tracking_root', '*.root'))
        if hits and m3:
            out.append((name, int(m.group(1))))
        elif os.path.isdir(os.path.join(run_dir, name)):
            print(f'  [skip] {name}: no hits/m3 root files (empty sub_run)')
    return sorted(out, key=lambda x: x[1])


def load_subrun(cfg, ct, subrun, z, chi2_cut, veto_sparks):
    """Return (m3, p2, hit_events) for one HV sub_run, optionally spark-vetoed."""
    sub = cfg.subrun_dir(subrun)
    hits_dir = os.path.join(sub, 'combined_hits_root')
    m3_dir = os.path.join(sub, 'm3_tracking_root')
    m3 = pa.load_m3_positions(m3_dir, z, chi2_cut)
    p2, hit_events = pa.load_p2_centroids(hits_dir, ct)
    n_veto = 0
    if veto_sparks:
        hv_csv = os.path.join(sub, 'hv_monitor.csv')
        if os.path.isfile(hv_csv):
            sv = ps.SparkVeto.from_csv(hv_csv, cfg)
            bad = sv.vetoed_ids_from_hits(hits_dir, ct.attrs['feus'])
            n_veto = len(bad)
            m3 = m3[~m3['eventId'].isin(bad)].copy()
            p2 = p2[~p2['eventId'].isin(bad)].copy()
            hit_events = hit_events - bad
    return m3, p2, hit_events, n_veto


def _robust_sigma(v):
    v = v[np.isfinite(v)]
    if len(v) < 10:
        return np.nan
    return float(1.4826 * np.median(np.abs(v - np.median(v))))


def main():
    ap = argparse.ArgumentParser(description='P2 efficiency vs mesh HV.')
    ap.add_argument('run_key', nargs='?', default='det1_hvscan')
    ap.add_argument('--strategy', default='reverse',
                    choices=['linear', 'reverse', 'pairswap'])
    ap.add_argument('--r', type=float, default=None,
                    help='match radius [mm]; default = run-config MATCH_R.')
    ap.add_argument('--active-r', type=float, default=30.0,
                    help='ray is in the active area if within this of a pad [mm].')
    ap.add_argument('--z', type=float, default=None,
                    help='M3 projection plane z [mm]; default = run-config '
                         'det_plane_z (measured PLANE_Z wins over nominal).')
    ap.add_argument('--chi2-cut', type=float, default=qa.M3_CHI2_CUT)
    ap.add_argument('--fit-fiducial', type=float, default=300.0,
                    help='|x_m3|,|y_m3| window used only to fit the transform.')
    ap.add_argument('--min-valid', type=int, default=50,
                    help='min matched events for an HV point to be plotted.')
    ap.add_argument('--veto-sparks', action=argparse.BooleanOptionalAction,
                    default=True, help='per-sub_run HV spark veto (default on).')
    args = ap.parse_args()

    cfg = qa.get_config(args.run_key)
    print(cfg)
    if args.r is None:
        args.r = cfg.MATCH_R
    if args.z is None:
        args.z = cfg.det_plane_z()
    print(f'match radius R = {args.r:g} mm, projection z = {args.z:g} mm')
    out_dir = cfg.out_dir('11_hv_scan_efficiency')
    suffix = cfg.product_suffix(args.veto_sparks)

    ct = pmap.build_channel_table(cfg.run_config_path, cfg.MAP_CSV_PATH,
                                  det_type=cfg.DET_TYPE, det_name=cfg.DET_NAME,
                                  strategy=args.strategy,
                                  drop_connectors=cfg.DEAD_CONNECTORS)

    subruns = find_subruns(cfg)
    if not subruns:
        print(f'No mesh_<NNN>V sub_runs under {cfg.run_dir}')
        return
    print(f'HV-scan points ({len(subruns)}): ' +
          ', '.join(f'{hv}V' for _, hv in subruns))

    # --- load every HV point, pool matched events for one transform fit ---
    data = []
    pooled = []
    for subrun, hv in subruns:
        m3, p2, hit_events, n_veto = load_subrun(cfg, ct, subrun, args.z,
                                                 args.chi2_cut, args.veto_sparks)
        matched = m3.merge(p2, on='eventId', how='inner')
        print(f'  mesh {hv}V: M3 {len(m3):,} | P2 {len(p2):,} | matched '
              f'{len(matched):,}' + (f' | spark-vetoed {n_veto}' if n_veto else ''))
        data.append(dict(subrun=subrun, hv=hv, m3=m3, p2=p2,
                         hit_events=hit_events, matched=matched))
        pooled.append(matched)
    pool = pd.concat(pooled, ignore_index=True)
    fitm = pool
    if args.fit_fiducial > 0:
        fitm = pool[(pool['x_m3'].abs() < args.fit_fiducial) &
                    (pool['y_m3'].abs() < args.fit_fiducial)]
    T = pa.fit_transform(fitm['x_m3'], fitm['y_m3'], fitm['x_pad'], fitm['y_pad'])
    print(f'\nPooled transform (all HV): rotation {T.rotation_deg:.2f} deg, '
          f'reflection {T.reflection}, scale {T.s:.3f}, RMSE {T.rmse:.1f} mm '
          f'(N={len(fitm):,})')

    # frozen active area from the transformed pad footprint (same for all HV)
    padc = ct.drop_duplicates('channel_id')
    pcx, pcy = T.apply(padc['pad_cx'].to_numpy(), padc['pad_cy'].to_numpy())
    tree = cKDTree(np.column_stack([pcx, pcy]))

    # --- per-HV efficiency inside the frozen active area ---
    rows = []
    for a in data:
        m3, p2 = a['m3'], a['p2']
        reco = {int(e): (xx, yy) for e, xx, yy in zip(
            p2['eventId'], *T.apply(p2['x_pad'].to_numpy(), p2['y_pad'].to_numpy()))}
        d = m3.rename(columns={'x_m3': 'x', 'y_m3': 'y'}).copy()
        if len(d) == 0:
            continue
        nn, _ = tree.query(np.column_stack([d['x'], d['y']]))
        d['in_active'] = nn <= args.active_r
        d['has_any'] = d['eventId'].isin(a['hit_events'])
        within = np.zeros(len(d), dtype=bool)
        dxres = np.full(len(d), np.nan); dyres = np.full(len(d), np.nan)
        ev = d['eventId'].to_numpy(); xs = d['x'].to_numpy(); ys = d['y'].to_numpy()
        for i in range(len(d)):
            rc = reco.get(int(ev[i]))
            if rc is not None:
                dxres[i] = xs[i] - rc[0]; dyres[i] = ys[i] - rc[1]
                within[i] = (dxres[i] ** 2 + dyres[i] ** 2) ** 0.5 <= args.r
        d['within'] = within

        da = d[d['in_active']]
        n = len(da)
        if n < args.min_valid:
            print(f'  [skip] mesh {a["hv"]}V: only {n} active-area tracks')
            continue
        eff = da['within'].mean(); err = np.sqrt(eff * (1 - eff) / n)
        eff_any = da['has_any'].mean()
        sx = _robust_sigma(dxres[d['in_active'].to_numpy()])
        sy = _robust_sigma(dyres[d['in_active'].to_numpy()])
        rows.append(dict(hv=a['hv'], subrun=a['subrun'], n_active=n,
                         n_within=int(da['within'].sum()),
                         eff_reco=eff, eff_reco_err=err, eff_anyhit=eff_any,
                         sigma_x_mm=sx, sigma_y_mm=sy))
        print(f'  mesh {a["hv"]}V: eff(reco<{args.r:g}mm)={eff:.3f}+-{err:.3f}  '
              f'eff(any)={eff_any:.3f}  sigma=({sx:.1f},{sy:.1f})mm  '
              f'({int(da["within"].sum())}/{n})')

    if not rows:
        print('No HV points had enough active-area tracks.')
        return
    df = pd.DataFrame(rows).sort_values('hv').reset_index(drop=True)
    df.to_csv(os.path.join(out_dir, f'efficiency_vs_hv{suffix}.csv'), index=False)

    tag = f'{cfg.DET_TAG} {cfg.RUN}'
    # efficiency vs HV
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.errorbar(df['hv'], df['eff_reco'], yerr=df['eff_reco_err'], fmt='o-',
                color='steelblue', capsize=4, lw=2, ms=7,
                label=f'reco within {args.r:g} mm')
    ax.plot(df['hv'], df['eff_anyhit'], 's--', color='darkorange', ms=6, alpha=0.8,
            label='any pad fired')
    ax.set_xlabel('mesh HV [V]'); ax.set_ylabel('efficiency (fixed active area)')
    ax.set_ylim(0, 1.02); ax.grid(True, alpha=0.3); ax.legend()
    ax.set_title(f'{cfg.DET_NAME} efficiency vs mesh HV — {tag}\n'
                 f'(r<{args.r:g} mm, frozen active area, drift = mesh+180 V)')
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, f'efficiency_vs_hv{suffix}.png'),
                dpi=200, bbox_inches='tight'); plt.close(fig)

    # resolution vs HV
    fig2, ax2 = plt.subplots(figsize=(8, 5))
    ax2.plot(df['hv'], df['sigma_x_mm'], 'o-', color='steelblue', lw=2, ms=7, label='σ_x')
    ax2.plot(df['hv'], df['sigma_y_mm'], 's--', color='darkorange', lw=2, ms=6, label='σ_y')
    ax2.set_xlabel('mesh HV [V]'); ax2.set_ylabel('core residual σ [mm]')
    ax2.set_ylim(0, None); ax2.grid(True, alpha=0.3); ax2.legend()
    ax2.set_title(f'{cfg.DET_NAME} residual width vs mesh HV — {tag}\n'
                  f'(robust σ of aligned P2−M3 residual, active area)')
    fig2.tight_layout()
    fig2.savefig(os.path.join(out_dir, f'resolution_vs_hv{suffix}.png'),
                 dpi=200, bbox_inches='tight'); plt.close(fig2)

    print(f'\n{"HV[V]":>6}  {"eff_reco":>9}  {"+-err":>6}  {"eff_any":>8}  {"tracks":>7}')
    for _, r in df.iterrows():
        print(f'{r.hv:>6.0f}  {r.eff_reco:>9.3f}  {r.eff_reco_err:>6.3f}  '
              f'{r.eff_anyhit:>8.3f}  {r.n_active:>7.0f}')
    print(f'\nWritten to: {out_dir}')


if __name__ == '__main__':
    main()
