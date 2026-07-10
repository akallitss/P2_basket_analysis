#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
12_validation.py

Structural validation of the P2 efficiency chain. The goal is NOT to re-plot the
efficiency but to prove the number is *depicted in the data* and is not being
pushed low (or high) by an algorithmic choice. Every test has a predetermined
pass condition; a failing test points at a specific bias.

It deliberately reuses the SAME functions as the real pipeline (p2_align.fit_transform,
load_m3_positions/endpoints, load_p2_centroids; p2_mapping.build_channel_table) so
the validator cannot diverge from what it validates.

Tests (this first cut — the highest-value, fastest trio + cheap redundancy):
  1. KNOB PLATEAU     eff vs match-radius R and vs active-area radius. An honest
                      cut sits on a plateau. The R plot also overlays the ideal
                      Gaussian containment 1-exp(-R^2/2 sigma^2): the part of the
                      "loss" that is pure geometry of the cut, not the detector.
  2. z PLATEAU        eff vs M3 projection plane z (re-project + re-fit the
                      transform at each z). Shows whether the shipped z sits on
                      the plateau or is leaving efficiency on the table.
  3. ALIGNMENT        free-scale fit vs scale-fixed-to-1 (pads are physical truth)
                      vs robust outlier-rejected fit. Reports core sigma and the
                      efficiency each gives. If cleaning the fit RAISES eff, the
                      shipped number was an alignment artefact.
  4. MAPPING vs DEAD  per-pad firing fraction vs per-pad local efficiency. A real
                      dead region is dark in BOTH; a mapping bug is dark in eff
                      but bright in firing (pads fire, reco lands elsewhere).
  5. ESTIMATOR        charge-weighted centroid-within-R vs nearest-fired-pad-within-R.
                      Quantifies how much the centroiding choice itself costs.
  6. STABILITY        odd/even-event split — must agree within statistics.

Products (<Analysis>/<detN>/<run>/<sub_run>/12_validation/):
  validation_<test>.png, validation_report.pdf, validation_summary.txt

Usage:
  python3 12_validation.py [run_key] [--z Z] [--strategy reverse] [--r 20] [--active-r 30]
"""
import os
import json
import argparse
import datetime

import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.spatial import cKDTree

import p2_qa_config as qa
import p2_mapping as pmap
import p2_align as pa
import p2_sparks as ps


# --------------------------------------------------------------------------- #
# helpers (mirror the pipeline)
# --------------------------------------------------------------------------- #
def det_plane_z(cfg):
    """Same z the pipeline stages use: measured PLANE_Z (03 z-scan) wins over
    the run_config det_center z."""
    return cfg.det_plane_z()


def robust_sigma(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if len(v) < 10:
        return np.nan
    return float(1.4826 * np.median(np.abs(v - np.median(v))))


def fix_scale_one(T0):
    """Same rotation as the free fit but scale forced to 1 (pads = physical truth)."""
    T = pa.Transform(T0.R, 1.0, T0.pad_mean, T0.m3_mean)
    return T


def robust_fit(x_m3, y_m3, x_pad, y_pad, n_iter=3, k=3.0):
    """Iterative outlier-rejected rigid+scale fit: down-weight rays whose residual
    is > k*MAD from the median, refit. Guards the core against tail contamination."""
    x_m3, y_m3 = np.asarray(x_m3, float), np.asarray(y_m3, float)
    x_pad, y_pad = np.asarray(x_pad, float), np.asarray(y_pad, float)
    mask = np.ones(len(x_m3), dtype=bool)
    T = pa.fit_transform(x_m3, y_m3, x_pad, y_pad)
    for _ in range(n_iter):
        rx, ry = T.apply(x_pad, y_pad)
        r = np.hypot(x_m3 - rx, y_m3 - ry)
        med, mad = np.median(r[mask]), np.median(np.abs(r[mask] - np.median(r[mask])))
        thr = med + k * 1.4826 * (mad if mad > 0 else 1.0)
        newmask = r <= thr
        if newmask.sum() < 100 or newmask.sum() == mask.sum():
            mask = newmask
            break
        mask = newmask
        T = pa.fit_transform(x_m3[mask], y_m3[mask], x_pad[mask], y_pad[mask])
    return T, mask


def reco_and_tree(T, p2, ct):
    """P2 reco positions (event->xy, M3 frame) + transformed-pad footprint KD-tree."""
    rx, ry = T.apply(p2['x_pad'].to_numpy(), p2['y_pad'].to_numpy())
    recodf = pd.DataFrame({'eventId': p2['eventId'].to_numpy(),
                           'rx': rx, 'ry': ry})
    padc = ct.drop_duplicates('channel_id')
    pcx, pcy = T.apply(padc['pad_cx'].to_numpy(), padc['pad_cy'].to_numpy())
    tree = cKDTree(np.column_stack([pcx, pcy]))
    return recodf, tree, np.column_stack([pcx, pcy]), padc


def ray_table(m3, recodf, tree, hit_events):
    """Per-ray table: nn_dist to nearest pad, |r| residual (+dx,dy), has_any.
    R / active_r cuts are applied downstream so scans are pure thresholding."""
    d = m3.rename(columns={'x_m3': 'x', 'y_m3': 'y'}).merge(recodf, on='eventId', how='left')
    nn_dist, nn_idx = tree.query(np.column_stack([d['x'], d['y']]))
    d['nn_dist'] = nn_dist
    d['nn_idx'] = nn_idx
    d['dx'] = d['x'] - d['rx']
    d['dy'] = d['y'] - d['ry']
    d['resid'] = np.hypot(d['dx'], d['dy'])
    d['has_any'] = d['eventId'].isin(hit_events)
    return d


def eff_of(d, R, active_r):
    sub = d[d['nn_dist'] <= active_r]
    n = len(sub)
    if n == 0:
        return float('nan'), float('nan'), 0
    within = (sub['resid'] <= R).mean() * 100
    any_ = sub['has_any'].mean() * 100
    return within, any_, n


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description='Structural validation of P2 efficiency.')
    ap.add_argument('run_key', nargs='?', default='det1_long')
    ap.add_argument('--strategy', default='reverse',
                    choices=['linear', 'reverse', 'pairswap'])
    ap.add_argument('--r', type=float, default=None,
                    help='match radius [mm]; default = run-config MATCH_R.')
    ap.add_argument('--active-r', type=float, default=30.0)
    ap.add_argument('--chi2-cut', type=float, default=qa.M3_CHI2_CUT)
    ap.add_argument('--fit-fiducial', type=float, default=300.0)
    ap.add_argument('--z', type=float, default=None,
                    help='M3 projection z; default = 06 lookup (run_config / 232).')
    ap.add_argument('--veto-sparks', action=argparse.BooleanOptionalAction, default=True)
    args = ap.parse_args()

    cfg = qa.get_config(args.run_key)
    print(cfg)
    out = cfg.out_dir('12_validation')
    sfx = cfg.product_suffix(args.veto_sparks)
    z0 = args.z if args.z is not None else det_plane_z(cfg)
    R0 = args.r if args.r is not None else cfg.MATCH_R
    AR0 = args.active_r
    verdicts = []   # (test, status, message)

    ct = pmap.build_channel_table(cfg.run_config_path, cfg.MAP_CSV_PATH,
                                  det_type=cfg.DET_TYPE, det_name=cfg.DET_NAME,
                                  strategy=args.strategy,
                                  drop_connectors=cfg.DEAD_CONNECTORS)
    hits_dir, m3_dir = cfg.combined_hits_dir, cfg.m3_tracking_dir

    print(f'Loading M3 (z={z0:.0f}) + P2 centroids ...')
    m3 = pa.load_m3_positions(m3_dir, z0, args.chi2_cut)
    p2, hit_events = pa.load_p2_centroids(hits_dir, ct)
    if args.veto_sparks:
        sv = ps.SparkVeto.from_cfg(cfg)
        bad = sv.vetoed_ids_from_hits(hits_dir, ct.attrs['feus'])
        m3 = m3[~m3['eventId'].isin(bad)].copy()
        p2 = p2[~p2['eventId'].isin(bad)].copy()
        hit_events = hit_events - set(int(b) for b in bad)
        print(f'  spark veto: {100*(1-sv.live_fraction()):.2f}% deadtime removed')

    matched = m3.merge(p2, on='eventId', how='inner')
    fitm = matched[(matched['x_m3'].abs() < args.fit_fiducial) &
                   (matched['y_m3'].abs() < args.fit_fiducial)]
    T0 = pa.fit_transform(fitm['x_m3'], fitm['y_m3'], fitm['x_pad'], fitm['y_pad'])
    print(f'nominal transform: rot {T0.rotation_deg:.2f} deg, scale {T0.s:.3f}, '
          f'RMSE {T0.rmse:.1f} mm')

    recodf, tree, pad_xy, padc = reco_and_tree(T0, p2, ct)
    d = ray_table(m3, recodf, tree, hit_events)
    eff0, any0, n0 = eff_of(d, R0, AR0)
    core_sx = robust_sigma(d.loc[d['nn_dist'] <= AR0, 'dx'])
    core_sy = robust_sigma(d.loc[d['nn_dist'] <= AR0, 'dy'])
    print(f'SHIPPED point: eff(reco<={R0:g})={eff0:.1f}%  any={any0:.1f}%  '
          f'core sigma=({core_sx:.1f},{core_sy:.1f}) mm  (N={n0:,})')

    # ---------------------------------------------------------------- Test 1 #
    R_scan = [8, 10, 12, 14, 16, 18, 20, 25, 30, 35, 40, 50, 60]
    eR = [eff_of(d, R, AR0)[0] for R in R_scan]
    AR_scan = [15, 20, 25, 30, 35, 40, 50, 60]
    eAR = [eff_of(d, R0, ar)[0] for ar in AR_scan]

    # ideal Gaussian containment scaled to the R-plateau (uses core sigma)
    sig = float(np.hypot(core_sx, core_sy) / np.sqrt(2)) if np.isfinite(core_sx) else 10.0
    plateau = eR[R_scan.index(50)]
    contain = [plateau * (1 - np.exp(-(R ** 2) / (2 * sig ** 2))) for R in R_scan]

    eff_plateau = eR[R_scan.index(40)]
    r_gap = eff_plateau - eff0
    contain0 = 1 - np.exp(-(R0 ** 2) / (2 * sig ** 2))   # fraction of a Gaussian inside R0
    verdicts.append(('R-cut plateau',
                     'FLAG' if r_gap > 3 else 'PASS',
                     f'eff rises {eff0:.1f}%→{eff_plateau:.1f}% from R={R0:g}→40 mm '
                     f'(+{r_gap:.1f} pt). The R={R0:g} mm cut geometry alone discards '
                     f'~{100*(1-contain0):.0f}% of genuine hits into the tail (σ≈{sig:.1f} mm).'))
    ar_var = max(eAR[AR_scan.index(20):AR_scan.index(40)+1]) - \
        min(eAR[AR_scan.index(20):AR_scan.index(40)+1])
    verdicts.append(('active-area plateau',
                     'FLAG' if ar_var > 2 else 'PASS',
                     f'eff varies {ar_var:.1f} pt over active_r 20–40 mm '
                     f'(shipped {AR0:g} mm → {eff0:.1f}%).'))

    fig, axs = plt.subplots(1, 2, figsize=(13, 5))
    axs[0].plot(R_scan, eR, 'o-', color='steelblue', lw=2, label='measured eff(reco)')
    axs[0].plot(R_scan, contain, '--', color='grey',
                label=f'ideal Gaussian containment (σ≈{sig:.1f} mm)')
    axs[0].axvline(R0, color='crimson', ls=':', label=f'shipped R={R0:g} mm')
    axs[0].axhline(any0, color='darkorange', ls='--', alpha=0.7,
                   label=f'any-pad ceiling = {any0:.1f}%')
    axs[0].set_xlabel('match radius R [mm]'); axs[0].set_ylabel('efficiency [%]')
    axs[0].set_title('eff vs match radius — is R on the plateau?')
    axs[0].grid(True, alpha=0.3); axs[0].legend(fontsize=8)
    axs[1].plot(AR_scan, eAR, 's-', color='seagreen', lw=2)
    axs[1].axvline(AR0, color='crimson', ls=':', label=f'shipped active_r={AR0:g} mm')
    axs[1].set_xlabel('active-area radius [mm]'); axs[1].set_ylabel('efficiency [%]')
    axs[1].set_title('eff vs active-area radius (denominator)')
    axs[1].grid(True, alpha=0.3); axs[1].legend(fontsize=8)
    fig.suptitle(f'{cfg.DET_NAME} — knob plateau tests  ({cfg.RUN})')
    fig.tight_layout()
    fig.savefig(f'{out}/validation_knob_plateau{sfx}.png', dpi=150, bbox_inches='tight')
    plt.close(fig)

    # ---------------------------------------------------------------- Test 2 #
    print('z-scan (re-project + re-fit at each z) ...')
    ep = pa.load_m3_endpoints(m3_dir, args.chi2_cut)
    keep = set(int(e) for e in m3['eventId'])
    ep = ep[ep['eventId'].isin(keep)].copy()
    z_scan = list(range(226, 259, 4))
    if int(round(z0)) not in z_scan:
        z_scan = sorted(z_scan + [int(round(z0))])
    eZ, sZ = [], []
    for z in z_scan:
        zx, zy = pa.project_to_z(ep, z)
        mz = pd.DataFrame({'eventId': ep['eventId'].to_numpy(), 'x_m3': zx, 'y_m3': zy})
        mz = mz[np.isfinite(mz['x_m3']) & np.isfinite(mz['y_m3'])]
        mm = mz.merge(p2, on='eventId', how='inner')
        fm = mm[(mm['x_m3'].abs() < args.fit_fiducial) & (mm['y_m3'].abs() < args.fit_fiducial)]
        Tz = pa.fit_transform(fm['x_m3'], fm['y_m3'], fm['x_pad'], fm['y_pad'])
        rc, tr, _, _ = reco_and_tree(Tz, p2, ct)
        dz = ray_table(mz, rc, tr, hit_events)
        e, _, _ = eff_of(dz, R0, AR0)
        eZ.append(e)
        sZ.append(robust_sigma(dz.loc[dz['nn_dist'] <= AR0, 'dx']))
    z_best = z_scan[int(np.nanargmax(eZ))]
    e_best = np.nanmax(eZ)
    verdicts.append(('z plateau',
                     'FLAG' if (e_best - eff0) > 3 else 'PASS',
                     f'shipped z={z0:.0f} → {eff0:.1f}%; eff peaks at z={z_best} → '
                     f'{e_best:.1f}% (+{e_best-eff0:.1f} pt).'))
    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax.plot(z_scan, eZ, 'o-', color='purple', lw=2)
    ax.axvline(z0, color='crimson', ls=':', label=f'shipped z={z0:.0f} mm')
    ax.axvline(z_best, color='green', ls='--', label=f'best z={z_best} mm')
    ax.set_xlabel('M3 projection plane z [mm]'); ax.set_ylabel('efficiency (reco) [%]')
    ax.set_title(f'{cfg.DET_NAME} — eff vs projection z (re-fit at each z)\n{cfg.RUN}')
    ax.grid(True, alpha=0.3); ax.legend()
    fig.tight_layout()
    fig.savefig(f'{out}/validation_z_plateau{sfx}.png', dpi=150, bbox_inches='tight')
    plt.close(fig)

    # ---------------------------------------------------------------- Test 3 #
    print('alignment variants (free scale / scale=1 / robust) ...')
    variants = {'free scale': T0, 'scale = 1': fix_scale_one(T0)}
    Trob, rmask = robust_fit(fitm['x_m3'].to_numpy(), fitm['y_m3'].to_numpy(),
                             fitm['x_pad'].to_numpy(), fitm['y_pad'].to_numpy())
    variants['robust fit'] = Trob
    align_rows = []
    for name, T in variants.items():
        rc, tr, _, _ = reco_and_tree(T, p2, ct)
        dd = ray_table(m3, rc, tr, hit_events)
        e, a, n = eff_of(dd, R0, AR0)
        sx = robust_sigma(dd.loc[dd['nn_dist'] <= AR0, 'dx'])
        sy = robust_sigma(dd.loc[dd['nn_dist'] <= AR0, 'dy'])
        align_rows.append(dict(name=name, scale=T.s, rmse=getattr(T, 'rmse', np.nan),
                               eff=e, any=a, sx=sx, sy=sy))
    best_align = max(align_rows, key=lambda r: r['eff'])
    verdicts.append(('alignment scale/robust',
                     'FLAG' if (best_align['eff'] - eff0) > 3 else 'PASS',
                     f'free-scale (s={T0.s:.3f}) → {eff0:.1f}%; best variant '
                     f'"{best_align["name"]}" → {best_align["eff"]:.1f}% '
                     f'(+{best_align["eff"]-eff0:.1f} pt, core σx {best_align["sx"]:.1f} mm).'))
    fig, ax = plt.subplots(figsize=(8, 4.6)); ax.axis('off')
    tbl = [['variant', 'scale', 'fit RMSE', 'eff(reco)', 'any', 'σx', 'σy']]
    for r in align_rows:
        tbl.append([r['name'], f'{r["scale"]:.3f}', f'{r["rmse"]:.1f}',
                    f'{r["eff"]:.1f}%', f'{r["any"]:.1f}%',
                    f'{r["sx"]:.1f}', f'{r["sy"]:.1f}'])
    t = ax.table(cellText=tbl, loc='center', cellLoc='center')
    t.auto_set_font_size(False); t.set_fontsize(10); t.scale(1, 1.6)
    for j in range(len(tbl[0])):
        t[0, j].set_facecolor('#1f3a5f'); t[0, j].set_text_props(color='w', weight='bold')
    ax.set_title(f'{cfg.DET_NAME} — alignment variants: does cleaning the fit raise eff?\n'
                 f'{cfg.RUN}', fontsize=11)
    fig.tight_layout()
    fig.savefig(f'{out}/validation_alignment{sfx}.png', dpi=150, bbox_inches='tight')
    plt.close(fig)

    # ---------------------------------------------------------------- Test 4 #
    print('mapping-vs-dead: per-pad firing vs per-pad efficiency ...')
    # per-pad hit counts from the raw hits (attach pads, count by channel_id)
    import glob
    import uproot
    feu_set = set(ct.attrs['feus'])
    parts = []
    for fp in sorted(glob.glob(os.path.join(hits_dir, '*.root'))):
        a = uproot.open(f'{fp}:hits').arrays(['eventId', 'channel', 'amplitude', 'feu'],
                                             library='pd')
        parts.append(a[a['feu'].isin(feu_set)].copy())
    hits = pmap.attach_pads_to_hits(pd.concat(parts, ignore_index=True), ct)
    hits = hits[hits['mapped'] & hits['pad_cx'].notna()]
    if args.veto_sparks:
        hits = hits[~hits['eventId'].isin(bad)]
    n_events = m3['eventId'].nunique()
    fire = hits.groupby('channel_id').size().rename('n_fire')
    padc = padc.merge(fire, on='channel_id', how='left')
    padc['n_fire'] = padc['n_fire'].fillna(0)
    padc['fire_frac'] = padc['n_fire'] / n_events
    # per-pad local efficiency: assign each active ray to its nearest pad (nn_idx)
    da = d[d['nn_dist'] <= AR0]
    padeff = da.groupby('nn_idx')['resid'].apply(lambda s: (s <= R0).mean())
    nrays = da.groupby('nn_idx').size()
    padc = padc.reset_index(drop=True)
    padc['loc_eff'] = padc.index.map(padeff).astype(float)
    padc['n_rays'] = padc.index.map(nrays).astype(float)
    # dead vs mapping discriminator: do the low-efficiency pads still FIRE?
    #   dead region      -> low eff AND low firing (both suppressed)
    #   mapping/seam bug  -> low eff BUT normal firing (pads fire, reco lands elsewhere)
    seen = padc[padc['n_rays'] >= 10]
    lowmask = seen['loc_eff'] < 0.15
    good_fire = seen.loc[seen['loc_eff'] > 0.5, 'fire_frac'].median()
    n_low = int(lowmask.sum())
    low_fire = seen.loc[lowmask, 'fire_frac'].median() if n_low else 0.0
    ratio = (low_fire / good_fire) if good_fire and np.isfinite(good_fire) else 0.0
    mapping = ratio > 0.5 and n_low >= 10
    verdicts.append(('mapping vs dead',
                     'FLAG' if mapping else 'PASS',
                     f'{n_low} pads have local eff<15%; they fire at {100*ratio:.0f}% of the '
                     f'good-pad rate → '
                     f'{"MAPPING/seam suspect (pads fire, reco lands elsewhere)" if mapping else "consistent with genuinely dead (firing also suppressed)"}.'))
    fig, axs = plt.subplots(1, 2, figsize=(14, 6))
    for ax, col, ttl, cm, vlim in [
            (axs[0], 'fire_frac', 'per-pad firing fraction (raw hits)', 'viridis', None),
            (axs[1], 'loc_eff', 'per-pad local efficiency (nearest-ray)', 'RdYlGn', (0, 1))]:
        sc = ax.scatter(padc['pad_cx'], padc['pad_cy'], c=padc[col], s=14, cmap=cm,
                        vmin=(vlim[0] if vlim else None), vmax=(vlim[1] if vlim else None))
        plt.colorbar(sc, ax=ax, label=col)
        ax.set_aspect('equal'); ax.set_xlabel('pad_cx [mm]'); ax.set_ylabel('pad_cy [mm]')
        ax.set_title(ttl, fontsize=10)
    # mark the dead-looking-but-firing pads on the efficiency panel
    hot = seen[lowmask & (seen['fire_frac'] >= 0.5 * good_fire)]
    axs[1].scatter(hot['pad_cx'], hot['pad_cy'], s=60, facecolors='none',
                   edgecolors='blue', lw=1.2, label='low-eff but firing')
    if len(hot):
        axs[1].legend(fontsize=8, loc='upper right')
    fig.suptitle(f'{cfg.DET_NAME} — dead region or mapping bug?  '
                 f'(dark in BOTH = dead; dark-eff/bright-fire = mapping)\n{cfg.RUN}')
    fig.tight_layout()
    fig.savefig(f'{out}/validation_mapping_vs_dead{sfx}.png', dpi=150, bbox_inches='tight')
    plt.close(fig)

    # ---------------------------------------------------------------- Test 5 #
    print('estimator: centroid-within-R vs nearest-fired-pad-within-R ...')
    hx, hy = T0.apply(hits['pad_cx'].to_numpy(), hits['pad_cy'].to_numpy())
    hit_xy = {}
    for ev, x, y in zip(hits['eventId'].to_numpy(), hx, hy):
        hit_xy.setdefault(int(ev), []).append((x, y))
    hit_xy = {k: np.asarray(v) for k, v in hit_xy.items()}
    within_nn = np.zeros(len(da), dtype=bool)
    for i, (ev, x, y) in enumerate(zip(da['eventId'].to_numpy(),
                                       da['x'].to_numpy(), da['y'].to_numpy())):
        pts = hit_xy.get(int(ev))
        if pts is not None:
            within_nn[i] = np.min(np.hypot(pts[:, 0] - x, pts[:, 1] - y)) <= R0
    eff_nn = 100 * within_nn.mean()
    verdicts.append(('estimator (centroid vs nearest pad)',
                     'FLAG' if (eff_nn - eff0) > 3 else 'PASS',
                     f'centroid-within-R {eff0:.1f}% vs nearest-fired-pad-within-R '
                     f'{eff_nn:.1f}% (+{eff_nn-eff0:.1f} pt): centroiding '
                     f'{"costs efficiency" if eff_nn-eff0>3 else "is not the limitation"}.'))

    # ---------------------------------------------------------------- Test 6 #
    parity = (d['eventId'] % 2).astype(int)
    eO = eff_of(d[parity == 1], R0, AR0)
    eE = eff_of(d[parity == 0], R0, AR0)
    diff = abs(eO[0] - eE[0])
    err = np.hypot(np.sqrt(eO[0]*(100-eO[0])/max(eO[2], 1)),
                   np.sqrt(eE[0]*(100-eE[0])/max(eE[2], 1)))
    verdicts.append(('odd/even stability',
                     'FLAG' if diff > 3 * err else 'PASS',
                     f'odd {eO[0]:.1f}% vs even {eE[0]:.1f}% (Δ={diff:.1f} pt, '
                     f'{diff/err:.1f}σ).'))

    # ---------------------------------------------------------------- report #
    lines = [
        f'P2 EFFICIENCY VALIDATION — {cfg.DET_TAG} {cfg.RUN}/{cfg.SUB_RUN} [{args.strategy}]',
        f'  shipped point: z={z0:.0f} mm, R={R0:g} mm, active_r={AR0:g} mm',
        f'  eff(reco) = {eff0:.1f}%   any-pad = {any0:.1f}%   '
        f'core σ = ({core_sx:.1f}, {core_sy:.1f}) mm   N_active = {n0:,}',
        f'  nominal transform: rotation {T0.rotation_deg:.2f}°, scale {T0.s:.3f}, '
        f'fit RMSE {T0.rmse:.1f} mm',
        '',
        f'  {"TEST":<34} {"VERDICT":<6}  DETAIL',
        '  ' + '-' * 96,
    ]
    for name, status, msg in verdicts:
        lines.append(f'  {name:<34} {status:<6}  {msg}')
    flags = [v for v in verdicts if v[1] == 'FLAG']
    lines += ['', f'  {len(flags)}/{len(verdicts)} tests FLAG a possible bias.']
    if flags:
        upper = max([eff0] + [best_align['eff'], e_best, eff_nn, eff_plateau])
        lines.append(f'  Best honest estimate of the position-reco efficiency lies between the '
                     f'shipped {eff0:.1f}% and ~{upper:.1f}% (any-pad ceiling {any0:.1f}%).')
    report = '\n'.join(lines)
    print('\n' + report)
    with open(f'{out}/validation_summary{sfx}.txt', 'w') as f:
        f.write(report + '\n')

    # one-page PDF: verdict table + the four diagnostic figures
    pdf_path = os.path.join(qa.DATA_ROOT, 'Analysis', cfg.DET_TAG,
                            f'p2_{cfg.DET_TAG}_validation.pdf')
    with PdfPages(pdf_path) as pdf:
        fig = plt.figure(figsize=(8.27, 11.69))
        fig.text(0.06, 0.965, 'P2 BASKET — Detector 1', fontsize=24, fontweight='bold',
                 color='#1f3a5f', va='top')
        fig.text(0.06, 0.928, 'Efficiency validation — structural bias audit', fontsize=13,
                 color='#555', va='top')
        fig.text(0.06, 0.905, f'{cfg.RUN}   shipped: z={z0:.0f} mm · R={R0:g} mm · '
                 f'eff {eff0:.1f}% (any {any0:.1f}%)', fontsize=9, family='monospace',
                 color='#333', va='top')
        ax = fig.add_axes([0.06, 0.66, 0.88, 0.20]); ax.axis('off')
        tb = [['Test', 'Verdict']]
        for name, status, msg in verdicts:
            tb.append([name, status])
        tt = ax.table(cellText=tb, loc='center', cellLoc='left', colWidths=[0.72, 0.28])
        tt.auto_set_font_size(False); tt.set_fontsize(8.5); tt.scale(1, 1.4)
        for j in range(2):
            tt[0, j].set_facecolor('#1f3a5f'); tt[0, j].set_text_props(color='w', weight='bold')
        for i, (_, status, _) in enumerate(verdicts, start=1):
            col = '#f7d7d5' if status == 'FLAG' else '#d9efd9'
            for j in range(2):
                tt[i, j].set_facecolor(col)
        # details as wrapped bullets under the table
        det = '\n'.join(f'• {name}: {msg}' for name, status, msg in verdicts)
        fig.text(0.06, 0.635, det, fontsize=8.2, va='top', ha='left', wrap=True,
                 linespacing=1.5, color='#1a1a1a')
        n_flag = sum(1 for v in verdicts if v[1] == 'FLAG')
        fig.text(0.06, 0.15,
                 f'How to read this: every test has a predetermined pass condition, so a FLAG points at a\n'
                 f'specific way the pipeline could be depressing (or inflating) the efficiency — not a matter\n'
                 f'of opinion. PASS = the number is stable against that choice.  '
                 f'{n_flag}/{len(verdicts)} tests flag a bias;\n'
                 f'the honest position-reco efficiency lies between the shipped {eff0:.1f}% and ~'
                 f'{eff_plateau:.1f}% (any-pad ceiling {any0:.1f}%).',
                 fontsize=8.6, va='top', color='#222',
                 bbox=dict(boxstyle='round,pad=0.6', fc='#eef3f8', ec='#9db8d2'))
        fig.text(0.96, 0.008, datetime.date.today().isoformat(), ha='right',
                 fontsize=6, color='grey')
        pdf.savefig(fig, dpi=200); plt.close(fig)
        for png in ['validation_knob_plateau', 'validation_z_plateau',
                    'validation_alignment', 'validation_mapping_vs_dead']:
            p = f'{out}/{png}{sfx}.png'
            if os.path.isfile(p):
                fig = plt.figure(figsize=(11.69, 8.27))
                a = fig.add_axes([0, 0, 1, 1]); a.axis('off')
                a.imshow(plt.imread(p))
                pdf.savefig(fig, dpi=150); plt.close(fig)
    print(f'\nPlots + summary -> {out}\nReport PDF -> {pdf_path}')


if __name__ == '__main__':
    main()
