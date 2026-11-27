#!/usr/bin/env python3
"""50_deform_scan_microtpc.py -- drift time vs drift field from pad PAIRS, at hit level.

The micro-TPC of stage 48 (an inclined track that crosses two pads: the pad holding the upper part of the
path fires later by about half the drift time) needs only one time per pad, so it can be done on every
point of a drift scan, where no waveforms were kept. That gives the drift time T(dV) with NO pulse-shape
response and NO zero-time reference: the two nuisance parameters that let the peak-time scan (stage 44 +
gas_studies/water_scan) trade the gap against the water content of the gas.

Per sub-run: adjacent pad pairs under an M3 track (both hits > 120 ADC, not saturated), slope s of the track
from the larger-pulse pad to the other, dt = time_of_max(B) - time_of_max(A);
D = [median dt (s > S_MIN) - median dt (s < -S_MIN)] / 2,  D = kappa * T / 2  (kappa = dilution by transverse
diffusion and thresholds; forward model for |s| > 0.15 with the slope mix of the data: 0.67 with dry-gas
diffusion, 0.76 with the diffusion of 0.5-0.75 % water, x 0.85).

Then T(dV) = d / v_f(dV / d) is fitted with the Magboltz water family of gas_studies/water_scan (f = water
fraction, d = gap of the zone), for kappa fixed at several values and free.

Zones: pads of the nominal-gap class of stage 49 (3.6-4.4 mm on the dry-gas scale), and all pads.
Out (<scan run>/50_deform_scan_microtpc/): pairs.csv.gz, D_vs_dV.csv, fit.csv, scan_microtpc.png, summary.txt
Usage: python3 50_deform_scan_microtpc.py [scan_key] [ref_key]
"""
import glob
import importlib
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'gas_studies', 'water_scan'))
import p2_qa_config as qa       # noqa: E402
import p2_align as pa           # noqa: E402
import p2_io as p2io            # noqa: E402
import p2_sparks as ps          # noqa: E402
import p2_mapping as pmap       # noqa: E402
import p2_study_cache as sc     # noqa: E402
s42 = importlib.import_module('42_pillar_shadow_vs_field')

AMP_B, S_MIN, NEAR = 120.0, 0.15, 22.0
RNG = np.random.default_rng(20261128)


def pairs_of(cfg, ct, d, T, z):
    hits_dir = os.path.join(d, 'combined_hits_root'); m3_dir = os.path.join(d, 'm3_tracking_root')
    drop = set(p2io.drop_pads_for(cfg, ct, hits_dir=hits_dir))
    ep = pa.load_m3_endpoints(m3_dir, qa.M3_CHI2_CUT, qa.M3_MIN_NCLUS)
    bad = set()
    hv = os.path.join(d, 'hv_monitor.csv')
    if os.path.isfile(hv):
        bad = set(ps.SparkVeto.from_csv(hv, cfg).vetoed_ids_from_hits(hits_dir, ct.attrs['feus'], min_amp=cfg.MIN_AMP))
    ep = ep[~ep.eventId.isin(bad)].drop_duplicates('eventId')
    x, y = pa.project_to_z(ep, z)
    pxy = s42.inverse(T, np.column_stack([x, y]))
    dz = ep['Z_Up'].to_numpy() - ep['Z_Down'].to_numpy()
    tm = np.column_stack([(ep['X_Up'].to_numpy() - ep['X_Down'].to_numpy()) / dz,
                          (ep['Y_Up'].to_numpy() - ep['Y_Down'].to_numpy()) / dz])
    tp = (tm / T.s) @ T.R
    trk = pd.DataFrame({'eventId': ep.eventId.to_numpy(np.int64), 'xp': pxy[:, 0], 'yp': pxy[:, 1],
                        'txp': tp[:, 0], 'typ': tp[:, 1]})
    parts = []
    for df in p2io.iter_hits(hits_dir, ['eventId', 'channel', 'amplitude', 'feu', 'time_of_max', 'time', 'saturated'],
                             ct.attrs['feus'], min_amp=cfg.MIN_AMP, progress=False):
        h = pmap.attach_pads_to_hits(df, ct)
        h = h[h['mapped'] & h['pad_cx'].notna() & ~h['channel_id'].isin(drop)]
        parts.append(h[['eventId', 'channel_id', 'pad_cx', 'pad_cy', 'amplitude', 'time_of_max', 'time', 'saturated']])
    if not parts:
        return pd.DataFrame(), len(ep)
    h = pd.concat(parts, ignore_index=True)
    h['eventId'] = h.eventId.astype(np.int64)
    h = h.merge(trk, on='eventId')
    h = h[(np.hypot(h.pad_cx - h.xp, h.pad_cy - h.yp) < NEAR) & (h.amplitude > AMP_B) & ~h.saturated.astype(bool)]
    h = h.sort_values(['eventId', 'amplitude'], ascending=[True, False]).drop_duplicates(['eventId', 'channel_id'])
    h['rk'] = h.groupby('eventId').cumcount()
    a = h[h.rk == 0].set_index('eventId'); b = h[h.rk == 1].set_index('eventId')
    p = a.join(b, how='inner', lsuffix='_a', rsuffix='_b')
    p['dist'] = np.hypot(p.pad_cx_a - p.pad_cx_b, p.pad_cy_a - p.pad_cy_b)
    p = p[p.dist < 14.5].copy()
    nx = (p.pad_cx_b - p.pad_cx_a) / p.dist; ny = (p.pad_cy_b - p.pad_cy_a) / p.dist
    p['s'] = p.txp_a * nx + p.typ_a * ny
    p['u'] = (p.xp_a - (p.pad_cx_a + p.pad_cx_b) / 2) * nx + (p.yp_a - (p.pad_cy_a + p.pad_cy_b) / 2) * ny
    p['dt'] = p.time_of_max_b - p.time_of_max_a
    p['dt30'] = p.time_b - p.time_a
    p = p[p.u.abs() < 3.0]
    return p[['channel_id_a', 'channel_id_b', 's', 'u', 'amplitude_a', 'amplitude_b', 'dt', 'dt30']].reset_index(), len(ep)


def dval(p, col='dt', nboot=300):
    a = p[p.s > S_MIN][col].to_numpy(); b = p[p.s < -S_MIN][col].to_numpy()
    if len(a) < 8 or len(b) < 8:
        return np.nan, np.nan, len(a) + len(b)
    D = 0.5 * (np.median(a) - np.median(b))
    bs = [0.5 * (np.median(a[RNG.integers(0, len(a), len(a))]) - np.median(b[RNG.integers(0, len(b), len(b))]))
          for _ in range(nboot)]
    return D, float(np.std(bs)), len(a) + len(b)


def main():
    scan_key = sys.argv[1] if len(sys.argv) > 1 else 'det2_initial2'
    ref_key = sys.argv[2] if len(sys.argv) > 2 else 'det2_long_repos'
    cfg = qa.get_config(scan_key); cfg_r = qa.get_config(ref_key)
    out = os.path.join(cfg.run_dir.replace(cfg.DATA_ROOT, os.path.dirname(os.path.dirname(os.path.dirname(cfg.OUT_BASE)))), '')
    out = os.path.abspath(os.path.join(cfg.OUT_BASE, '..', '50_deform_scan_microtpc')); os.makedirs(out, exist_ok=True)
    ct = sc.channel_table(cfg)
    T, z, nfit = s42.align_series(cfg, ct)
    L = [f'{scan_key}: alignment of the series from its alignment run (trimmed fit, z = {z:.1f} mm, scale {T.s:.4f}); '
         f'pairs: both hits > {AMP_B:g} ADC, adjacent pads, track within 3 mm of the shared edge; D from |slope| > {S_MIN}']
    gp = pd.read_csv(cfg_r.out_dir('49_deform_decomposition') + '/pads.csv').set_index('channel_id').gap
    cache = f'{out}/pairs.csv.gz'
    if os.path.isfile(cache) and '--rebuild' not in sys.argv:
        P = pd.read_csv(cache)
    else:
        subs = [('alignment', 430, 730, cfg.sub_run_dir)]
        for d in sorted(glob.glob(os.path.join(cfg.run_dir, '*_scan_*'))):
            nm = os.path.basename(d); a_, b_ = nm.split('_')[-2:]
            subs.append(('drift' if nm.startswith('drift') else 'mesh', int(a_), int(b_), d))
        parts = []
        for kind, mesh, drift, d in subs:
            if not glob.glob(os.path.join(d, 'combined_hits_root', '*.root')):
                continue
            p, ntr = pairs_of(cfg, ct, d, T, z)
            if len(p):
                p['kind'], p['mesh'], p['drift'], p['ntrk'] = kind, mesh, drift, ntr
                parts.append(p)
            print(f'{kind:9s} {mesh}/{drift}: {ntr:6d} tracks, {len(p):4d} pairs', flush=True)
        P = pd.concat(parts, ignore_index=True)
        P.to_csv(cache, index=False)
    P['gap'] = 0.5 * (gp.reindex(P.channel_id_a.to_numpy()).to_numpy() + gp.reindex(P.channel_id_b.to_numpy()).to_numpy())
    P['dV'] = P.drift - P.mesh
    zones = {'nominal': (P.gap >= 3.6) & (P.gap < 4.4), 'all': P.gap > 0, 'ridge': P.gap >= 5.5, 'thin': P.gap < 3.6}
    rows = []
    for zn, m in zones.items():
        q = P[m & (P.mesh == 430)]
        # the alignment run and the drift-scan point at the same voltages are pooled
        for dV, g in q.groupby('dV'):
            D, e, n = dval(g); D30, e30, _ = dval(g, 'dt30')
            rows.append(dict(zone=zn, dV=dV, E_nom=dV / 0.4, D=D, e=e, D30=D30, e30=e30, n_pairs=n, n_all=len(g)))
    R = pd.DataFrame(rows); R.to_csv(f'{out}/D_vs_dV.csv', index=False)
    L.append('D(dV) [ns], nominal-gap pads (dry-scale gap 3.6-4.4 mm of stage 49):')
    for r in R[R.zone == 'nominal'].itertuples():
        L.append(f'   dV {r.dV:5.0f} V  E_nom {r.E_nom:6.0f} V/cm  D = {r.D:7.1f} +- {r.e:5.1f} (peak)  {r.D30:7.1f} +- {r.e30:5.1f} (30 % edge)  pairs {r.n_pairs}')

    # ---- fit T(dV) = d / v_f(dV/d), D = kappa T / 2 ---------------------------------------------
    import analyze_water as aw
    gas = aw.load_gases()
    tabs = sorted([(f, gas[k]) for k, f in aw.FAMILIES['h2o'] if k in gas], key=lambda a: a[0])
    Fg = np.arange(0.0, 2.01, 0.05); Dg = np.arange(3.0, 6.51, 0.1)
    fits = []
    for zn in ('nominal', 'all'):
        q = R[(R.zone == zn) & R.D.notna() & (R.dV > 20)]
        dV, y, ye = q.dV.to_numpy(), q.D.to_numpy(), np.hypot(q.e.to_numpy(), 2.0)
        grid = []
        for f in Fg:
            gm = aw.GasModel.mix(tabs, f)
            for d in Dg:
                Th = 0.5 * (d / 10.0) / gm.v(dV / (d / 10.0))
                kfree = float((y * Th / ye ** 2).sum() / (Th * Th / ye ** 2).sum())
                row = dict(zone=zn, f=f, d=d, kappa_free=kfree, chi2_free=float((((y - kfree * Th) / ye) ** 2).sum()))
                for k in (0.67, 0.76, 0.85):
                    row[f'chi2_k{k}'] = float((((y - k * Th) / ye) ** 2).sum())
                grid.append(row)
        G = pd.DataFrame(grid); fits.append(G)
        L.append(f'fit, zone {zn}: {len(q)} points (dV > 20 V), errors + 2 ns in quadrature')
        for col, nm in (('chi2_k0.67', 'kappa = 0.67 (dry diffusion)'), ('chi2_k0.76', 'kappa = 0.76 (wet-gas diffusion)'), ('chi2_k0.85', 'kappa = 0.85'), ('chi2_free', 'kappa free')):
            b = G.loc[G[col].idxmin()]; ok = G[G[col] < G[col].min() + 2.3]
            L.append(f'   {nm:13s}: water {b.f:.2f} % ({ok.f.min():.2f}-{ok.f.max():.2f}), gap {b.d:.1f} mm ({ok.d.min():.1f}-{ok.d.max():.1f}), '
                     f'chi2 {b[col]:.1f}/{len(q) - (3 if col == "chi2_free" else 2)}' + (f', kappa {b.kappa_free:.2f}' if col == 'chi2_free' else ''))
        d0 = G[np.isclose(G.f, 0.0)]
        L.append(f'   dry gas, any gap, kappa free: best chi2 {d0.chi2_free.min():.0f}')
    Gall = pd.concat(fits); Gall.to_csv(f'{out}/fit.csv', index=False)

    fig, ax = plt.subplots(1, 3, figsize=(21, 5.6))
    q = R[R.zone == 'nominal']
    ax[0].errorbar(q.E_nom, q.D, q.e, fmt='o', color='#257de0', label='nominal-gap pads')
    q2 = R[R.zone == 'ridge']
    ax[0].errorbar(q2.E_nom, q2.D, q2.e, fmt='s', color='#0e3d7e', mfc='none', label='ridge pads')
    G = fits[0]; b = G.loc[G['chi2_k0.76'].idxmin()]
    xx = np.linspace(25, 525, 200)
    for (f, d, ls, lab) in ((b.f, b.d, '-', f'best: {b.f:.2f} % water, {b.d:.1f} mm'), (0.75, 4.0, '--', '0.75 % water, 4.0 mm'), (0.0, 4.0, ':', 'dry, 4.0 mm')):
        gm = aw.GasModel.mix(tabs, f)
        ax[0].plot(xx / 0.4, 0.76 * 0.5 * (d / 10.0) / gm.v(xx / (d / 10.0)), ls, color='#5c5c5c', label=lab + ' (kappa 0.76)')
    ax[0].set(yscale='log', xlabel='nominal drift field dV / 4 mm [V/cm]', ylabel='pad-pair time difference D [ns]',
              title='Half the drift time, straight from pad pairs'); ax[0].legend(frameon=False, fontsize=8)
    for a_, col, ttl in ((ax[1], 'chi2_k0.76', 'kappa = 0.76 (forward model, wet-gas diffusion)'), (ax[2], 'chi2_free', 'kappa free')):
        Z = G.pivot(index='d', columns='f', values=col); zz = Z.to_numpy() - Z.to_numpy().min()
        a_.contourf(Z.columns, Z.index, zz, levels=[0, 2.3, 6.2, 11.8], colors=['#257de0', '#84b9f3', '#cce3fd'])
        a_.set(xlabel='water fraction [%]', ylabel='gap of the nominal zone [mm]', title=f'68 / 95 / 99.7 % regions, {ttl}')
        a_.axhline(4.0, color='#9c9c9c', lw=1)
    for a_ in ax:
        a_.grid(color='#e6e6e6')
    fig.tight_layout(); fig.savefig(f'{out}/scan_microtpc.png', dpi=110); plt.close(fig)
    with open(f'{out}/summary.txt', 'w') as fh:
        fh.write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
