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
import p2_io as p2io

_HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_GAS_TABLE = os.path.join(_HERE, 'garfield_inputs', 'gas_table.csv')


def load_drift_velocity(gas_csv):
    """Magboltz drift-velocity table for the fill gas: return (v_d(E), E_min,
    E_max) where v_d maps E [V/cm] -> v_d [um/ns] (linear interp). Outside the
    tabulated range it returns NaN below E_min (v_d -> 0 as E -> 0, but that
    region is NOT tabulated so we do not fake a flat value) and clips flat above
    E_max. gas_table.csv has E_Vcm and vz [cm/ns] (see 14_timing_simulation.py)."""
    if not os.path.isfile(gas_csv):
        return None, None, None
    g = pd.read_csv(gas_csv).sort_values('E_Vcm')
    e = g['E_Vcm'].to_numpy(float)
    v = g['vz'].to_numpy(float) * 1e4          # cm/ns -> um/ns
    return (lambda E: np.interp(E, e, v, left=np.nan, right=v[-1]),
            float(e[0]), float(e[-1]))


def find_subruns(cfg, scan='drift'):
    """Discover <scan>_scan sub_runs and this detector's (mesh, drift) from the
    sub_run name (…<det_tag>_<mesh>_<drift>…), ascending in the scanned V.

    scan='drift': drift_scan_* names, x = drift (mesh fixed).
    scan='mesh' : mesh_scan_*  names, x = mesh (drift moves in tandem)."""
    pat = re.compile(rf'{cfg.DET_TAG}_(\d+)_(\d+)')
    out = []
    for name in sorted(os.listdir(cfg.run_dir)):
        if not name.startswith(f'{scan}_scan'):
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
    return sorted(out, key=lambda x: x[1] if scan == 'mesh' else x[2])


def load_subrun(cfg, ct, subrun, z, chi2_cut, veto_sparks, sig_amp=300.0):
    """Return (m3, p2, hit_events, n_veto, amp, tspread) for one scan point.

    amp is the per-point pad-amplitude summary (gain proxy); tspread is the
    signal-hit peak-time spread (drift-time proxy). Both are over the same hits
    that feed the efficiency (hot pads dropped, spark-vetoed events excluded)."""
    sub = cfg.subrun_dir(subrun)
    hits_dir = os.path.join(sub, 'combined_hits_root')
    m3_dir = os.path.join(sub, 'm3_tracking_root')
    drop = p2io.drop_pads_for(cfg, ct, hits_dir=hits_dir)
    m3 = pa.load_m3_positions(m3_dir, z, chi2_cut)
    p2, hit_events = pa.load_p2_centroids(hits_dir, ct, min_amp=cfg.MIN_AMP,
                                          drop_pads=drop)
    n_veto = 0
    bad = None
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
    amp = p2io.pad_amp_stats(hits_dir, ct, min_amp=cfg.MIN_AMP,
                             drop_pads=drop, exclude_events=bad)
    tspread = p2io.pad_time_spread(hits_dir, ct, sig_amp=sig_amp,
                                   min_amp=cfg.MIN_AMP, drop_pads=drop,
                                   exclude_events=bad)
    tres = p2io.event_time_resolution(hits_dir, ct, sig_amp=sig_amp,
                                      min_amp=cfg.MIN_AMP, drop_pads=drop,
                                      exclude_events=bad)
    return m3, p2, hit_events, n_veto, amp, tspread, tres


def _robust_sigma(v):
    v = v[np.isfinite(v)]
    if len(v) < 10:
        return np.nan
    return float(1.4826 * np.median(np.abs(v - np.median(v))))


def _best_timing_window(d, col='time_res_ns', frac=1.10):
    """(vlo, vhi) of the efficient points whose timing metric is within `frac`
    of its minimum -- the operable best-timing band. None if no efficient pts."""
    good = d[(d['eff_reco'] > 0.5) & d[col].notna()]
    if not len(good):
        return None
    opt = good[good[col] <= frac * good[col].min()]
    return opt['x'].min(), opt['x'].max()


def _drift_velocity_plot(df, cfg, tag, out_dir, suffix, vd_fn, args):
    """Two-panel drift figure: (left) efficiency + nominal Magboltz v_d(E),
    (right) the measured peak-time spread (a drift-time proxy). The Magboltz
    curve is drawn ONLY over the tabulated field range (below it v_d is not
    tabulated and physically -> 0 as E -> 0; no flat extrapolation)."""
    d = df[df['eff_reco'].notna()].copy()
    gap_cm = args.drift_gap_mm / 10.0
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5.2))

    # -- left: efficiency + Magboltz v_d ---------------------------------- #
    axL.errorbar(d['x'], d['eff_reco'], yerr=d['eff_reco_err'], fmt='o-',
                 color='steelblue', capsize=4, lw=2, ms=7,
                 label=f'efficiency (reco <{args.r:g} mm)')
    axL.set_xlabel('drift HV [V]')
    axL.set_ylabel('efficiency (fixed active area)', color='steelblue')
    axL.tick_params(axis='y', labelcolor='steelblue')
    axL.set_ylim(0, 1.02); axL.grid(True, alpha=0.3)
    axv = axL.twinx()
    if vd_fn is not None:
        mesh = d['mesh'].iloc[0]
        vv = np.linspace(mesh + 1, d['drift'].max(), 400)
        vd = vd_fn((vv - mesh) / gap_cm)          # NaN below the table -> gap
        axv.plot(vv, vd, '-', color='seagreen', lw=2,
                 label='Magboltz $v_d$ (nominal gas)')
        if np.isfinite(vd).any():
            ipk = int(np.nanargmax(vd))
            axv.annotate(f'nominal $v_d$ peak\n{(vv[ipk]-mesh)/gap_cm:.0f} V/cm',
                         xy=(vv[ipk], vd[ipk]),
                         xytext=(vv[ipk] + 55, vd[ipk] * 0.70),
                         fontsize=7.5, color='seagreen', ha='left', va='top',
                         arrowprops=dict(arrowstyle='->', color='seagreen',
                                         lw=0.8))
            # honest low-field note: table stops, v_d -> 0 as E -> 0
            vlo = float(vv[np.isfinite(vd)][0])
            axv.annotate('table stops here\n($v_d\\!\\to\\!0$ as $E\\!\\to\\!0$)',
                         xy=(vlo, vd[np.isfinite(vd)][0]),
                         xytext=(vlo + 8, vd[np.isfinite(vd)][0] * 0.55),
                         fontsize=6.8, color='dimgrey', ha='left', va='top',
                         arrowprops=dict(arrowstyle='->', color='dimgrey',
                                         lw=0.7))
    axv.set_ylabel('drift velocity $v_d$ [µm/ns]', color='seagreen')
    axv.tick_params(axis='y', labelcolor='seagreen')
    axv.set_ylim(0, None)
    h1, l1 = axL.get_legend_handles_labels()
    h2, l2 = axv.get_legend_handles_labels()
    axL.legend(h1 + h2, l1 + l2, fontsize=8, loc='lower right')
    axL.set_title('Efficiency + nominal drift velocity', fontsize=11)

    # -- right: measured peak-time spread --------------------------------- #
    axR.plot(d['x'], d['tom_spread_ns'], 'o-', color='crimson', lw=2, ms=7,
             label='peak-time spread p90−p10 (data)')
    axR.set_xlabel('drift HV [V]')
    axR.set_ylabel('peak-time spread [ns]  (smaller = better timing)',
                   color='crimson')
    axR.tick_params(axis='y', labelcolor='crimson')
    axR.set_ylim(0, None); axR.grid(True, alpha=0.3)
    win = _best_timing_window(d, col='tom_spread_ns')
    if win:
        axR.axvspan(win[0] - 10, win[1] + 10, color='gold', alpha=0.20,
                    label='best-timing window (efficient)')
    axR.legend(fontsize=8, loc='upper right', framealpha=0.92)
    axR.set_title('Measured drift-time spread (all signal hits)', fontsize=11)

    fig.suptitle(
        f'{cfg.DET_NAME} drift velocity & timing vs drift HV — {tag}\n'
        f'(mesh {d["mesh"].iloc[0]:.0f} V fixed, {args.drift_gap_mm:g} mm gap; '
        'nominal-gas Magboltz vs measured peak-time spread)')
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, f'drift_velocity_vs_drift{suffix}.png'),
                dpi=200, bbox_inches='tight')
    plt.close(fig)


def _time_resolution_plot(df, cfg, tag, out_dir, suffix, args):
    """Estimated detector time resolution sigma_t (per-event leading-pad
    time_of_max, referenced to the scintillator trigger / readout start) vs
    drift HV, ZOOMED to the minimum region (the near-zero-field point, whose
    drift time blows up, is dropped). The efficient best-timing band is shaded.
    """
    all_d = df[df['time_res_ns'].notna() & (df['e_drift_Vcm'] > 50)].copy()
    if len(all_d) < 2:
        return
    smin = all_d['time_res_ns'].min()
    # zoom on the minimum: show only points within 2.5x of the best sigma_t
    # (drops the low-field turn-on knee that would otherwise squash the valley)
    d = all_d[all_d['time_res_ns'] <= 2.5 * smin].copy()
    n_hidden = len(all_d) - len(d)
    fig, ax = plt.subplots(figsize=(8.5, 5.4))
    ax.errorbar(d['x'], d['time_res_ns'], fmt='o-', color='indigo', lw=2, ms=8,
                capsize=4, label='σ$_t$ (leading-pad, robust)')
    ax.plot(d['x'], d['time_res_mad_ns'], '^:', color='mediumpurple', ms=6,
            alpha=0.7, label='σ$_t$ (MAD cross-check)')
    best = d.loc[d['time_res_ns'].idxmin()]
    span = d['time_res_ns'].max() - d['time_res_ns'].min() + 1
    ax.annotate(f'min σ$_t$ = {best.time_res_ns:.1f} ns\n'
                f'@ drift {best.x:.0f} V ({best.e_drift_Vcm:.0f} V/cm), '
                f'eff {best.eff_reco:.0%}',
                xy=(best.x, best.time_res_ns),
                xytext=(best.x, best.time_res_ns + 0.42 * span),
                fontsize=8.5, ha='center', color='indigo',
                arrowprops=dict(arrowstyle='->', color='indigo'))
    win = _best_timing_window(all_d, col='time_res_ns')
    if win:
        ax.axvspan(win[0] - 12, win[1] + 12, color='gold', alpha=0.22,
                   label='best-timing window (within 10% of min, efficient)')
    lo = d['time_res_ns'].min(); hi = d['time_res_ns'].max()
    ax.set_ylim(max(0, lo - 0.20 * (hi - lo) - 0.5), hi + 0.55 * (hi - lo) + 0.5)
    if n_hidden:
        ax.text(0.99, 0.02, f'{n_hidden} low-field turn-on point(s) '
                'above the zoom', transform=ax.transAxes, ha='right',
                va='bottom', fontsize=7, color='grey', style='italic')
    ax.set_xlabel('drift HV [V]')
    ax.set_ylabel('estimated time resolution σ$_t$ [ns]')
    ax.grid(True, alpha=0.3); ax.legend(fontsize=9, loc='upper left')
    ax.set_title(
        f'{cfg.DET_NAME} estimated time resolution vs drift HV — {tag}\n'
        f'(mesh {d["mesh"].iloc[0]:.0f} V fixed, per-event leading-pad '
        'time_of_max vs trigger; ~0 V drift point dropped)')
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, f'time_resolution_vs_drift{suffix}.png'),
                dpi=200, bbox_inches='tight')
    plt.close(fig)


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
    ap.add_argument('--min-amp', type=float, default=None,
                    help='override cfg.MIN_AMP [ADC] (e.g. 60: still above '
                         'the det3 noise-floor q99~48 but keeps more of the '
                         'low-gain signal Landau than the default 100).')
    ap.add_argument('--fit-top', type=int, default=3,
                    help='fit the pad->M3 transform on the N points with the '
                         'most matched events (0 = pool all points). Keeps '
                         'noise-matches from near-dead points out of the fit.')
    ap.add_argument('--scan', default='drift', choices=['drift', 'mesh'],
                    help='scanned voltage: drift (drift_scan_* sub_runs, mesh '
                         'fixed) or mesh (mesh_scan_* sub_runs, drift in '
                         'tandem). Default drift.')
    ap.add_argument('--drift-gap-mm', type=float, default=3.0,
                    help='conversion (drift) gap thickness [mm], for the drift '
                         'field E=(V_drift-V_mesh)/gap and the drift-velocity '
                         'overlay. Default 3.0 (P2 Micromegas).')
    ap.add_argument('--amp-gap-um', type=float, default=150.0,
                    help='amplification gap [um], only for the reported '
                         'avalanche field. Default 150 (P2 Micromegas).')
    ap.add_argument('--gas-table', default=DEFAULT_GAS_TABLE,
                    help='Magboltz drift-velocity table (E_Vcm, vz[cm/ns]) for '
                         'the fill gas; overlaid on the drift scan.')
    ap.add_argument('--sig-amp', type=float, default=300.0,
                    help='signal amplitude threshold [ADC] for the data '
                         'arrival-time (drift-time) spread measurement.')
    args = ap.parse_args()

    cfg = qa.get_config(args.run_key)
    print(cfg)
    if args.min_amp is not None:
        cfg.MIN_AMP = float(args.min_amp)
    if args.r is None:
        args.r = cfg.MATCH_R
    if args.z is None:
        args.z = cfg.det_plane_z()
    print(f'match radius R = {args.r:g} mm, projection z = {args.z:g} mm, '
          f'min_amp = {cfg.MIN_AMP:g} ADC')
    out_dir = cfg.out_dir(f'16_{args.scan}_scan_efficiency')
    suffix = cfg.product_suffix(args.veto_sparks)
    if args.min_amp is not None:
        suffix += f'_minamp{args.min_amp:g}'   # never overwrite the defaults

    # drift-field geometry + Magboltz drift velocity (drift scan only)
    gap_cm = args.drift_gap_mm / 10.0
    vd_fn = vd_emin = vd_emax = None
    if args.scan == 'drift':
        vd_fn, vd_emin, vd_emax = load_drift_velocity(args.gas_table)
        print(f'drift gap {args.drift_gap_mm:g} mm, amp gap {args.amp_gap_um:g} '
              f'um; drift field E=(V_drift-V_mesh)/{gap_cm:g}cm')
        if vd_fn is None:
            print(f'  [warn] no gas table at {args.gas_table} — vd overlay off')
        else:
            print(f'  Magboltz vd table covers E={vd_emin:.0f}-{vd_emax:.0f} '
                  f'V/cm (below E_min vd is NOT tabulated: vd->0 as E->0)')

    ct = pmap.build_channel_table(cfg.run_config_path, cfg.MAP_CSV_PATH,
                                  det_type=cfg.DET_TYPE, det_name=cfg.DET_NAME,
                                  strategy=args.strategy,
                                  drop_connectors=cfg.DEAD_CONNECTORS)
    if cfg.DEAD_CONNECTORS:
        print(f'  dropped dead connectors: {list(cfg.DEAD_CONNECTORS)}')

    subruns = find_subruns(cfg, args.scan)
    if not subruns:
        print(f'No {args.scan}_scan sub_runs for {cfg.DET_TAG} under {cfg.run_dir}')
        return
    if args.scan == 'mesh':
        fixed_lbl = (f'drift in tandem '
                     f'{subruns[0][2]}-{subruns[-1][2]} V')
        print(f'Mesh-scan points ({len(subruns)}, {fixed_lbl}): ' +
              ', '.join(f'{mv}V' for _, mv, _ in subruns))
    else:
        fixed_lbl = f'mesh {subruns[0][1]} V fixed'
        print(f'Drift-scan points ({len(subruns)}, {fixed_lbl}): ' +
              ', '.join(f'{dv}V' for _, _, dv in subruns))

    # --- load every point, pool matched events for one transform fit -------- #
    data = []
    pooled = []
    for subrun, mv, dv in subruns:
        m3, p2, hit_events, n_veto, amp, tspread, tres = load_subrun(
            cfg, ct, subrun, args.z, args.chi2_cut, args.veto_sparks,
            sig_amp=args.sig_amp)
        matched = m3.merge(p2, on='eventId', how='inner')
        xv = mv if args.scan == 'mesh' else dv
        print(f'  {args.scan} {xv}V: M3 {len(m3):,} | P2 {len(p2):,} | matched '
              f'{len(matched):,}' + (f' | spark-vetoed {n_veto}' if n_veto else '')
              + f' | pad amp mean {amp["mean"]:.0f} ADC'
              + f' | t_max spread {tspread["spread"]:.0f} ns'
              + f' | sigma_t {tres["sigma_ns"]:.1f} ns')
        data.append(dict(subrun=subrun, mesh=mv, drift=dv, x=xv, m3=m3, p2=p2,
                         hit_events=hit_events, matched=matched, amp=amp,
                         tspread=tspread, tres=tres))
        pooled.append(matched)
    # transform fit sample: only the points with real response — pooling the
    # near-dead points too mixes noise-matches into the fit and biases the
    # frozen active area (seen on the det3 mesh scan: scale 0.85/RMSE 73 mm
    # pooled vs 0.96/40 mm clean).
    if args.fit_top and len(data) > args.fit_top:
        best = sorted(data, key=lambda a: len(a['matched']),
                      reverse=True)[:args.fit_top]
        print(f'transform fit restricted to the {args.fit_top} points with '
              f'most matches: ' + ', '.join(a['subrun'] for a in best))
        pool = pd.concat([a['matched'] for a in best], ignore_index=True)
    else:
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
            print(f'  [skip] {args.scan} {a["x"]}V: only {n} active-area tracks')
            continue
        eff = da['within'].mean()
        err = np.sqrt(max(eff * (1 - eff), 1e-12) / n)
        eff_any = da['has_any'].mean()
        sx = _robust_sigma(dxres[d['in_active'].to_numpy()])
        sy = _robust_sigma(dyres[d['in_active'].to_numpy()])
        amp = a['amp']
        ts = a['tspread']
        tr = a['tres']
        # drift field + drift velocity (drift scan only). The measured p90-p10
        # window of the peak-time spans ~ 0.8 * d_gap / v_drift for a uniformly
        # illuminated gap, so an "apparent" v_d follows from the spread (carries
        # a constant diffusion+trigger-phase floor; the trend is the signal).
        e_drift = (a['drift'] - a['mesh']) / gap_cm if args.scan == 'drift' else np.nan
        vd_mag = float(vd_fn(e_drift)) if (vd_fn is not None and e_drift > 0) else np.nan
        vd_app = (0.8 * args.drift_gap_mm * 1e3 / ts['spread']
                  if (args.scan == 'drift' and ts['spread'] and ts['spread'] > 0)
                  else np.nan)
        rows.append(dict(x=a['x'], drift=a['drift'], mesh=a['mesh'],
                         subrun=a['subrun'],
                         n_active=n, n_within=int(da['within'].sum()),
                         n_p2_events=len(p2),
                         eff_reco=eff, eff_reco_err=err, eff_anyhit=eff_any,
                         sigma_x_mm=sx, sigma_y_mm=sy,
                         amp_mean=amp['mean'], amp_mean_err=amp['sem'],
                         amp_median=amp['median'], n_hits_amp=amp['n'],
                         e_drift_Vcm=e_drift, vd_magboltz_um_ns=vd_mag,
                         tom_spread_ns=ts['spread'], tom_iqr_ns=ts['iqr'],
                         vd_apparent_um_ns=vd_app,
                         time_res_ns=tr['sigma_ns'],
                         time_res_mad_ns=tr['mad_ns'], n_time_ev=tr['n_events']))
        print(f'  {args.scan} {a["x"]}V: eff(reco<{args.r:g}mm)={eff:.3f}+-{err:.3f}  '
              f'eff(any)={eff_any:.3f}  sigma=({sx:.1f},{sy:.1f})mm  '
              f'({int(da["within"].sum())}/{n})')

    if not rows:
        print(f'No {args.scan} points had enough active-area tracks.')
        return
    df = pd.DataFrame(rows).sort_values('x').reset_index(drop=True)
    df.to_csv(os.path.join(out_dir, f'efficiency_vs_{args.scan}{suffix}.csv'),
              index=False)

    tag = f'{cfg.DET_TAG} {cfg.RUN}'
    xlbl = f'{args.scan} HV [V]'
    # efficiency vs scanned V
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.errorbar(df['x'], df['eff_reco'], yerr=df['eff_reco_err'], fmt='o-',
                color='steelblue', capsize=4, lw=2, ms=7,
                label=f'reco within {args.r:g} mm')
    ax.plot(df['x'], df['eff_anyhit'], 's--', color='darkorange', ms=6,
            alpha=0.8, label='any pad fired')
    ax.set_xlabel(xlbl)
    ax.set_ylabel('efficiency (fixed active area)')
    ax.set_ylim(0, 1.02)
    ax.grid(True, alpha=0.3)
    ax.legend()
    ax.set_title(f'{cfg.DET_NAME} efficiency vs {args.scan} HV — {tag}\n'
                 f'({fixed_lbl}, r<{args.r:g} mm, '
                 f'min_amp {cfg.MIN_AMP:g} ADC, frozen active area)')
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, f'efficiency_vs_{args.scan}{suffix}.png'),
                dpi=200, bbox_inches='tight')
    plt.close(fig)

    # mean pad amplitude vs scanned V -- linear + log panels. In a MESH scan
    # amplitude tracks the (exponential) gas gain; in a DRIFT scan the gain is
    # fixed, so amplitude instead tracks charge collection / mesh transparency
    # through the drift gap -- keep the log-panel caption scan-aware.
    fig3, (ax3, ax3b) = plt.subplots(1, 2, figsize=(12, 5))
    for axi in (ax3, ax3b):
        axi.errorbar(df['x'], df['amp_mean'], yerr=df['amp_mean_err'], fmt='o-',
                     color='steelblue', capsize=4, lw=2, ms=7,
                     label='mean pad amp')
        axi.plot(df['x'], df['amp_median'], 's--', color='darkorange', ms=6,
                 alpha=0.8, label='median pad amp')
        axi.set_xlabel(xlbl); axi.set_ylabel('pad amplitude [ADC]')
        axi.grid(True, alpha=0.3); axi.legend()
    ax3.set_ylim(0, None); ax3.set_title('linear')
    ax3b.set_yscale('log')
    ax3b.set_title('log y — exponential gas gain is a straight line'
                   if args.scan == 'mesh'
                   else 'log y (gain fixed: amplitude = charge collection)')
    amp_note = ('mapped pad hits, hot pads + spark events removed'
                if args.scan == 'mesh' else
                'gain fixed at the mesh working point; amplitude tracks '
                'drift-gap charge collection')
    fig3.suptitle(f'{cfg.DET_NAME} mean pad amplitude vs {args.scan} HV — {tag}\n'
                  f'({fixed_lbl}, {amp_note})')
    fig3.tight_layout()
    fig3.savefig(os.path.join(out_dir, f'amplitude_vs_{args.scan}{suffix}.png'),
                 dpi=200, bbox_inches='tight')
    plt.close(fig3)

    # resolution vs scanned V
    fig2, ax2 = plt.subplots(figsize=(8, 5))
    ax2.plot(df['x'], df['sigma_x_mm'], 'o-', color='steelblue', lw=2,
             ms=7, label='σ_x')
    ax2.plot(df['x'], df['sigma_y_mm'], 's--', color='darkorange', lw=2,
             ms=6, label='σ_y')
    ax2.set_xlabel(xlbl)
    ax2.set_ylabel('core residual σ [mm]')
    ax2.set_ylim(0, None)
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    ax2.set_title(f'{cfg.DET_NAME} residual width vs {args.scan} HV — {tag}\n'
                  f'(robust σ of aligned P2−M3 residual, active area)')
    fig2.tight_layout()
    fig2.savefig(os.path.join(out_dir, f'resolution_vs_{args.scan}{suffix}.png'),
                 dpi=200, bbox_inches='tight')
    plt.close(fig2)

    # drift velocity & timing (drift scan only): the fill-gas Magboltz v_d(E)
    # mapped onto each drift point, plus the measured peak-time spread (a data
    # proxy for the drift-time = d_gap/v_d) -> which drift HV gives BEST timing.
    if args.scan == 'drift' and df['e_drift_Vcm'].notna().any():
        _drift_velocity_plot(df, cfg, tag, out_dir, suffix, vd_fn, args)
        _time_resolution_plot(df, cfg, tag, out_dir, suffix, args)

    print(f'\n{args.scan + "[V]":>8}  {"eff_reco":>9}  {"+-err":>6}  {"eff_any":>8}  '
          f'{"n_p2":>6}  {"tracks":>7}  {"amp_mean":>9}  {"amp_med":>8}')
    for _, r in df.iterrows():
        print(f'{r.x:>8.0f}  {r.eff_reco:>9.3f}  {r.eff_reco_err:>6.3f}  '
              f'{r.eff_anyhit:>8.3f}  {r.n_p2_events:>6.0f}  {r.n_active:>7.0f}  '
              f'{r.amp_mean:>9.0f}  {r.amp_median:>8.0f}')
    if args.scan == 'drift' and df['e_drift_Vcm'].notna().any():
        print(f'\n{"drift[V]":>8}  {"E[V/cm]":>8}  {"vd_mag":>7}  {"t_spread":>9}  '
              f'{"sigma_t":>8}   (vd um/ns, spread & sigma_t ns)')
        for _, r in df.iterrows():
            print(f'{r.x:>8.0f}  {r.e_drift_Vcm:>8.0f}  {r.vd_magboltz_um_ns:>7.1f}  '
                  f'{r.tom_spread_ns:>9.0f}  {r.time_res_ns:>8.1f}')
        good = df[(df['eff_reco'] > 0.5) & df['time_res_ns'].notna()]
        if len(good):
            best = good.loc[good['time_res_ns'].idxmin()]
            print(f'\n  BEST TIMING (min sigma_t) among efficient points: '
                  f'drift {best.x:.0f} V (E={best.e_drift_Vcm:.0f} V/cm), '
                  f'sigma_t {best.time_res_ns:.1f} ns, eff {best.eff_reco:.2f}')
    print(f'\nWritten to: {out_dir}')


if __name__ == '__main__':
    main()
