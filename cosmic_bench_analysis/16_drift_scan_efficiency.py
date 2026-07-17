#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
16_drift_scan_efficiency.py

P2 detector efficiency (and resolution) vs DRIFT voltage at fixed mesh HV,
following 11_hv_scan_efficiency.py (itself adapted from
nTof_x17/mx_june_cosmic_qa/10_hv_scan_efficiency.py).

The drift-scan runs take BOTH detectors simultaneously: one sub_run per
voltage point, named like

    drift_scan_det4_<mesh>_<drift>_det3_<mesh>_<drift>

so the same run is analysed twice, once per run_key (det3_driftscan1 /
det4_driftscan1); the sub_run name is parsed for THAT detector's (mesh, drift).

Method (identical to stage 11):
  1. Build the channel table once from the run-level run_config.json.
  2. Per sub_run: load M3 rays + P2 pad centroids (cfg.MIN_AMP kills the
     stale-pedestal noise floor), optional per-sub_run HV spark veto (burst
     counting also above cfg.MIN_AMP).
  3. Fit the pad->M3 transform ONCE on the pooled matched events (geometry
     does not change with HV) and freeze the transformed pad footprint ->
     one active area for every point.
  4. Per point: eff(reco within R), eff(any pad fired), robust residual sigma.

Products (<Analysis>/<detN>/<scan_run>/drift_scan/16_drift_scan_efficiency/):
  efficiency_vs_drift<suffix>.png / resolution_vs_drift<suffix>.png
  efficiency_vs_drift<suffix>.csv

Usage:
  python3 16_drift_scan_efficiency.py det4_driftscan1 [--r 20] [--active-r 30]
        [--z 232] [--fit-fiducial 300] [--min-valid 50] [--no-veto-sparks]
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
    """Discover drift-scan sub_runs and this detector's (mesh, drift) from the
    sub_run name (…_<det_tag>_<mesh>_<drift>…), ascending in drift V."""
    pat = re.compile(rf'{cfg.DET_TAG}_(\d+)_(\d+)')
    out = []
    for name in sorted(os.listdir(cfg.run_dir)):
        if not name.startswith('drift_scan'):
            continue
        m = pat.search(name)
        if not m:
            continue
        hits = glob.glob(os.path.join(cfg.run_dir, name,
                                      'combined_hits_root', '*.root'))
        m3 = glob.glob(os.path.join(cfg.run_dir, name,
                                    'm3_tracking_root', '*.root'))
        if hits and m3:
            out.append((name, int(m.group(1)), int(m.group(2))))
        else:
            print(f'  [skip] {name}: no hits/m3 root files yet')
    return sorted(out, key=lambda x: x[2])


def load_subrun(cfg, ct, subrun, z, chi2_cut, veto_sparks):
    """Return (m3, p2, hit_events, n_veto) for one drift point."""
    sub = cfg.subrun_dir(subrun)
    hits_dir = os.path.join(sub, 'combined_hits_root')
    m3_dir = os.path.join(sub, 'm3_tracking_root')
    m3 = pa.load_m3_positions(m3_dir, z, chi2_cut)
    p2, hit_events = pa.load_p2_centroids(hits_dir, ct, min_amp=cfg.MIN_AMP,
                                          drop_pads=cfg.NOISY_PADS)
    n_veto = 0
    if veto_sparks:
        hv_csv = os.path.join(sub, 'hv_monitor.csv')
        if os.path.isfile(hv_csv):
            sv = ps.SparkVeto.from_csv(hv_csv, cfg)
            bad = sv.vetoed_ids_from_hits(hits_dir, ct.attrs['feus'],
                                          min_amp=cfg.MIN_AMP)
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
    ap = argparse.ArgumentParser(description='P2 efficiency vs drift HV.')
    ap.add_argument('run_key', nargs='?', default='det4_driftscan1')
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
                    help='min active-area tracks for a point to be plotted.')
    ap.add_argument('--min-fit', type=int, default=100,
                    help='min pooled matched events to trust the transform fit.')
    ap.add_argument('--veto-sparks', action=argparse.BooleanOptionalAction,
                    default=True, help='per-sub_run HV spark veto (default on).')
    args = ap.parse_args()

    cfg = qa.get_config(args.run_key)
    print(cfg)
    if args.r is None:
        args.r = cfg.MATCH_R
    if args.z is None:
        args.z = cfg.det_plane_z()
    print(f'match radius R = {args.r:g} mm, projection z = {args.z:g} mm, '
          f'min_amp = {cfg.MIN_AMP:g} ADC')
    out_dir = cfg.out_dir('16_drift_scan_efficiency')
    suffix = cfg.product_suffix(args.veto_sparks)

    ct = pmap.build_channel_table(cfg.run_config_path, cfg.MAP_CSV_PATH,
                                  det_type=cfg.DET_TYPE, det_name=cfg.DET_NAME,
                                  strategy=args.strategy,
                                  drop_connectors=cfg.DEAD_CONNECTORS)
    if cfg.DEAD_CONNECTORS:
        print(f'  dropped dead connectors: {list(cfg.DEAD_CONNECTORS)}')

    subruns = find_subruns(cfg)
    if not subruns:
        print(f'No drift_scan sub_runs for {cfg.DET_TAG} under {cfg.run_dir}')
        return
    mesh_v = subruns[0][1]
    print(f'Drift-scan points ({len(subruns)}, mesh {mesh_v} V): ' +
          ', '.join(f'{dv}V' for _, _, dv in subruns))

    # --- load every point, pool matched events for one transform fit -------- #
    data = []
    pooled = []
    for subrun, mv, dv in subruns:
        m3, p2, hit_events, n_veto = load_subrun(cfg, ct, subrun, args.z,
                                                 args.chi2_cut, args.veto_sparks)
        matched = m3.merge(p2, on='eventId', how='inner')
        print(f'  drift {dv}V: M3 {len(m3):,} | P2 {len(p2):,} | matched '
              f'{len(matched):,}' + (f' | spark-vetoed {n_veto}' if n_veto else ''))
        data.append(dict(subrun=subrun, drift=dv, m3=m3, p2=p2,
                         hit_events=hit_events, matched=matched))
        pooled.append(matched)
    pool = pd.concat(pooled, ignore_index=True)
    fitm = pool
    if args.fit_fiducial > 0:
        fitm = pool[(pool['x_m3'].abs() < args.fit_fiducial) &
                    (pool['y_m3'].abs() < args.fit_fiducial)]
    if len(fitm) < args.min_fit:
        print(f'\nOnly {len(fitm)} pooled matched events (< {args.min_fit}) — '
              'not enough signal to fit the pad->M3 transform. The detector '
              'is not (yet) responding above the noise floor; no efficiency '
              'can be extracted from this scan.')
        return
    T = pa.fit_transform(fitm['x_m3'], fitm['y_m3'], fitm['x_pad'], fitm['y_pad'])
    print(f'\nPooled transform (all points): rotation {T.rotation_deg:.2f} deg, '
          f'reflection {T.reflection}, scale {T.s:.3f}, RMSE {T.rmse:.1f} mm '
          f'(N={len(fitm):,})')
    if not (0.7 < T.s < 1.4):
        print('  WARNING: fitted scale far from 1 — the "signal" sample is '
              'probably noise/discharges, not track-correlated. Treat the '
              'efficiency below as an upper limit on real response.')

    # frozen active area from the transformed pad footprint (same for all pts)
    padc = ct[ct['mapped']].drop_duplicates('channel_id')
    pcx, pcy = T.apply(padc['pad_cx'].to_numpy(), padc['pad_cy'].to_numpy())
    tree = cKDTree(np.column_stack([pcx, pcy]))

    # --- per-point efficiency inside the frozen active area ----------------- #
    rows = []
    for a in data:
        m3, p2 = a['m3'], a['p2']
        if len(m3) == 0:
            continue
        reco = {int(e): (xx, yy) for e, xx, yy in zip(
            p2['eventId'], *T.apply(p2['x_pad'].to_numpy(),
                                    p2['y_pad'].to_numpy()))}
        d = m3.rename(columns={'x_m3': 'x', 'y_m3': 'y'}).copy()
        nn, _ = tree.query(np.column_stack([d['x'], d['y']]))
        d['in_active'] = nn <= args.active_r
        d['has_any'] = d['eventId'].isin(a['hit_events'])
        within = np.zeros(len(d), dtype=bool)
        dxres = np.full(len(d), np.nan)
        dyres = np.full(len(d), np.nan)
        ev = d['eventId'].to_numpy()
        xs = d['x'].to_numpy()
        ys = d['y'].to_numpy()
        for i in range(len(d)):
            rc = reco.get(int(ev[i]))
            if rc is not None:
                dxres[i] = xs[i] - rc[0]
                dyres[i] = ys[i] - rc[1]
                within[i] = (dxres[i] ** 2 + dyres[i] ** 2) ** 0.5 <= args.r
        d['within'] = within

        da = d[d['in_active']]
        n = len(da)
        if n < args.min_valid:
            print(f'  [skip] drift {a["drift"]}V: only {n} active-area tracks')
            continue
        eff = da['within'].mean()
        err = np.sqrt(max(eff * (1 - eff), 1e-12) / n)
        eff_any = da['has_any'].mean()
        sx = _robust_sigma(dxres[d['in_active'].to_numpy()])
        sy = _robust_sigma(dyres[d['in_active'].to_numpy()])
        rows.append(dict(drift=a['drift'], mesh=mesh_v, subrun=a['subrun'],
                         n_active=n, n_within=int(da['within'].sum()),
                         n_p2_events=len(p2),
                         eff_reco=eff, eff_reco_err=err, eff_anyhit=eff_any,
                         sigma_x_mm=sx, sigma_y_mm=sy))
        print(f'  drift {a["drift"]}V: eff(reco<{args.r:g}mm)={eff:.3f}+-{err:.3f}  '
              f'eff(any)={eff_any:.3f}  sigma=({sx:.1f},{sy:.1f})mm  '
              f'({int(da["within"].sum())}/{n})')

    if not rows:
        print('No drift points had enough active-area tracks.')
        return
    df = pd.DataFrame(rows).sort_values('drift').reset_index(drop=True)
    df.to_csv(os.path.join(out_dir, f'efficiency_vs_drift{suffix}.csv'),
              index=False)

    tag = f'{cfg.DET_TAG} {cfg.RUN}'
    # efficiency vs drift V
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.errorbar(df['drift'], df['eff_reco'], yerr=df['eff_reco_err'], fmt='o-',
                color='steelblue', capsize=4, lw=2, ms=7,
                label=f'reco within {args.r:g} mm')
    ax.plot(df['drift'], df['eff_anyhit'], 's--', color='darkorange', ms=6,
            alpha=0.8, label='any pad fired')
    ax.set_xlabel('drift HV [V]')
    ax.set_ylabel('efficiency (fixed active area)')
    ax.set_ylim(0, 1.02)
    ax.grid(True, alpha=0.3)
    ax.legend()
    ax.set_title(f'{cfg.DET_NAME} efficiency vs drift HV — {tag}\n'
                 f'(mesh {mesh_v} V fixed, r<{args.r:g} mm, '
                 f'min_amp {cfg.MIN_AMP:g} ADC, frozen active area)')
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, f'efficiency_vs_drift{suffix}.png'),
                dpi=200, bbox_inches='tight')
    plt.close(fig)

    # resolution vs drift V
    fig2, ax2 = plt.subplots(figsize=(8, 5))
    ax2.plot(df['drift'], df['sigma_x_mm'], 'o-', color='steelblue', lw=2,
             ms=7, label='σ_x')
    ax2.plot(df['drift'], df['sigma_y_mm'], 's--', color='darkorange', lw=2,
             ms=6, label='σ_y')
    ax2.set_xlabel('drift HV [V]')
    ax2.set_ylabel('core residual σ [mm]')
    ax2.set_ylim(0, None)
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    ax2.set_title(f'{cfg.DET_NAME} residual width vs drift HV — {tag}\n'
                  f'(robust σ of aligned P2−M3 residual, active area)')
    fig2.tight_layout()
    fig2.savefig(os.path.join(out_dir, f'resolution_vs_drift{suffix}.png'),
                 dpi=200, bbox_inches='tight')
    plt.close(fig2)

    print(f'\n{"drift[V]":>8}  {"eff_reco":>9}  {"+-err":>6}  {"eff_any":>8}  '
          f'{"n_p2":>6}  {"tracks":>7}')
    for _, r in df.iterrows():
        print(f'{r.drift:>8.0f}  {r.eff_reco:>9.3f}  {r.eff_reco_err:>6.3f}  '
              f'{r.eff_anyhit:>8.3f}  {r.n_p2_events:>6.0f}  {r.n_active:>7.0f}')
    print(f'\nWritten to: {out_dir}')


if __name__ == '__main__':
    main()
