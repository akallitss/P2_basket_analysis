#!/usr/bin/env python3
"""42_pillar_shadow_vs_field.py -- the small-pillar shadow (on the MEASURED
lattice of stage 41) across a chamber's drift-field and mesh-voltage scans.

The lattice is fixed to the board, so stage 41's lattice_pad.npz (pad frame)
of a reference run is reused as is on every sub-run of the scan series; the
origin is cross-checked by a Fourier phase on the pooled vertical tracks.

Per scan point and pooled field / mesh groups, for vertical (< 6 deg), 6-12
and >= 12 deg tracks:
  no-hit vs distance to the nearest lattice point (0-0.5 / 0.5-1 / >= 1 mm),
  integrated shadow = no-hit(all) - no-hit(>= 1 mm),
  r_eff: hard dead disc (+) Gaussian pointing fitted to the 3-bin vertical
         profile with the far bin as baseline (1 sigma = chi2 + 1).

det2 (11-24 series, audit 2026-11-27): r_eff falls from 0.70 mm (< 200 V/cm)
to 0.25 mm at 800-1200 V/cm and rises again above 1200 V/cm (mesh
transparency); on the mesh scan it falls from 0.70 mm (340 V) to 0.3 mm
(420 V). The shadow depth follows the signal margin, not the geometry alone.

Usage: python3 42_pillar_shadow_vs_field.py [scan_key] [lattice_key]
   defaults: det2_initial2 (11-24 series)  det2_long_repos (stage 41 lattice)
Out: <lattice run>/../42_pillar_shadow_vs_field_<scan>/{shadow_vs_field.csv,
     shadow_vs_field.png, summary.txt}
"""
import glob
import os
import pickle
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2_qa_config as qa       # noqa: E402
import p2_align as pa           # noqa: E402
import p2_io as p2io            # noqa: E402
import p2_sparks as ps          # noqa: E402
import p2_study_cache as sc     # noqa: E402
import p2_study_common as cm    # noqa: E402
import p2_mapping as pmap       # noqa: E402

GRID = '#e6e6e6'
EDGES = np.array([0, 0.5, 1.0, 1.75])
TH = [(0, 6), (6, 12), (12, 90)]
FIELD_GROUPS = [(-1, 200), (200, 400), (400, 600), (600, 800), (800, 1000), (1000, 1200), (1200, 1400)]
RNG = np.random.default_rng(42)


def inverse(T, xy):
    return ((np.asarray(xy) - T.m3_mean) / T.s) @ T.R + T.pad_mean


def lattice_points(la, lo, hi, margin=5.0):
    n = 250
    ii = np.arange(-n, n + 1)
    I, J = np.meshgrid(ii, ii)
    L = la['origin'] + I.ravel()[:, None] * la['a1'] + J.ravel()[:, None] * la['a2']
    return L[(L > lo - margin).all(1) & (L < hi + margin).all(1)]


class DiscLookup:
    """3-bin profile of a hard disc of radius r on the lattice, pointing sigma."""

    def __init__(self, L, sigma, nmc=300000):
        c = L[len(L) // 2]
        Lc = L[np.hypot(*(L - c).T) < 15]
        t = c + RNG.uniform(-5, 5, (nmc, 2))
        self.dt = cKDTree(Lc).query(t)[0]
        dm = cKDTree(Lc).query(t + sigma * RNG.normal(0, 1, (nmc, 2)))[0]
        self.n, _ = np.histogram(dm, EDGES)
        self.rs = np.arange(0.0, 0.71, 0.01)
        self.M0 = {r: np.histogram(dm, EDGES, weights=(self.dt < r).astype(float))[0] / self.n for r in self.rs}

    def fit(self, Y, N):
        b = Y[2]
        w = N / np.clip(Y * (1 - Y), 1e-4, None)
        chis = {r: float((w[:2] * (Y[:2] - (b + (1 - b) * self.M0[r][:2])) ** 2).sum()) for r in self.rs}
        r0 = min(chis, key=chis.get)
        ok = [r for r in self.rs if chis[r] < chis[r0] + 1]
        return r0, min(ok), max(ok)


def align_series(cfg, ct):
    """Trimmed alignment + slope-fitted plane z on the series' alignment run."""
    ali = cfg.sub_run_dir
    ep = pa.load_m3_endpoints(os.path.join(ali, 'm3_tracking_root'), qa.M3_CHI2_CUT, qa.M3_MIN_NCLUS)
    drop = p2io.drop_pads_for(cfg, ct, hits_dir=os.path.join(ali, 'combined_hits_root'))
    cen, _ = pa.load_p2_centroids(os.path.join(ali, 'combined_hits_root'), ct, drop_pads=drop)
    tx, ty, _ = cm.ray_dirs(ep)
    z = float(cfg.det_plane_z())
    for _ in range(4):
        x, y = pa.project_to_z(ep, z)
        m = pd.DataFrame({'eventId': ep.eventId, 'x': x, 'y': y, 'tx': tx, 'ty': ty}).merge(cen, on='eventId')
        m = m[(m.x.abs() < 300) & (m.y.abs() < 300)]
        T = pa.fit_transform(m.x, m.y, m.x_pad, m.y_pad, trim_mm=15.0)
        px, py = T.apply(m.x_pad, m.y_pad)
        rx, ry = px - m.x, py - m.y
        k = np.hypot(rx, ry) < 15
        tt = np.r_[m.tx[k], m.ty[k]]; rr = np.r_[rx[k], ry[k]]
        z += float((tt * rr).sum() / (tt * tt).sum())
    return T, z, len(m)


def point(cfg, ct, d, T, z, Lfull, bigp):
    hits_dir = os.path.join(d, 'combined_hits_root'); m3_dir = os.path.join(d, 'm3_tracking_root')
    drop = p2io.drop_pads_for(cfg, ct, hits_dir=hits_dir)
    ep = pa.load_m3_endpoints(m3_dir, qa.M3_CHI2_CUT, qa.M3_MIN_NCLUS)
    _, hit_ev = pa.load_p2_centroids(hits_dir, ct, min_amp=cfg.MIN_AMP, drop_pads=drop)
    hv = os.path.join(d, 'hv_monitor.csv')
    if os.path.isfile(hv):
        sv = ps.SparkVeto.from_csv(hv, cfg)
        bad = sv.vetoed_ids_from_hits(hits_dir, ct.attrs['feus'], min_amp=cfg.MIN_AMP)
        ep = ep[~ep.eventId.isin(bad)]
    x, y = pa.project_to_z(ep, z); xy = np.column_stack([x, y])
    padc = ct[ct['mapped']].drop_duplicates('channel_id')
    pxy = np.column_stack(T.apply(padc.pad_cx.to_numpy(), padc.pad_cy.to_numpy()))
    nn, _ = cKDTree(pxy).query(xy)
    act = pa.mask_dead_cells(nn, xy, ct, T) <= cfg.ACTIVE_R
    _, _, th = cm.ray_dirs(ep)
    nohit = ~ep.eventId.isin(hit_ev).to_numpy()
    xp = inverse(T, xy)
    db = cKDTree(bigp).query(xp)[0] if len(bigp) else np.full(len(xp), np.inf)
    dl = cKDTree(Lfull).query(xp)[0]
    k = act & (db > 10)
    return dict(xp=xp[k], nohit=nohit[k], th=th[k], dl=dl[k])


def profile(sets, lo, hi):
    N = np.zeros(3); K = np.zeros(3)
    for s in sets:
        k = (s['th'] >= lo) & (s['th'] < hi)
        n, _ = np.histogram(s['dl'][k], EDGES); kk, _ = np.histogram(s['dl'][k], EDGES, weights=s['nohit'][k].astype(float))
        N += n; K += kk
    return N, K / np.maximum(N, 1), K.sum() / max(N.sum(), 1)


def main():
    scan_key = sys.argv[1] if len(sys.argv) > 1 else 'det2_initial2'
    lat_key = sys.argv[2] if len(sys.argv) > 2 else 'det2_long_repos'
    cfg = qa.get_config(scan_key)
    cfg_l = qa.get_config(lat_key)
    lat = os.path.join(cfg_l.out_dir('41_pillar_lattice'), 'lattice_pad.npz')
    if not os.path.isfile(lat):
        sys.exit(f'run 41_pillar_lattice.py {lat_key} first ({lat} missing)')
    la = dict(np.load(lat))
    out = os.path.abspath(os.path.join(cfg_l.OUT_BASE, '..', f'42_pillar_shadow_vs_field_{scan_key}'))
    os.makedirs(out, exist_ok=True)
    ct = sc.channel_table(cfg)
    ct.attrs['dead_guard'] = 2.0; ct.attrs['edge_guard'] = 2.0
    A = np.load(os.path.join(cfg_l.out_dir('34_clustering_matching'), 'alignment.npz'))
    sig_pt = float(A['sigma_point'])
    L_ = []
    T, z, nfit = align_series(cfg, ct)
    L_.append(f'{scan_key}: alignment run z = {z:.1f} mm, scale {T.s:.4f}, rot {T.rotation_deg:.2f} deg, '
              f'rmse {T.rmse:.2f} mm (n {nfit}); lattice from {lat_key} stage 41 (pitch {float(la["pitch"]):.3f} mm)')
    pil = pmap.load_pillars(cfg.MASK_GBR_PATH)
    bigp = pil[pil.big][['x', 'y']].to_numpy() if len(pil) else np.empty((0, 2))
    padc = ct[ct['mapped']].drop_duplicates('channel_id')
    lo = padc[['pad_cx', 'pad_cy']].to_numpy().min(0); hi = padc[['pad_cx', 'pad_cy']].to_numpy().max(0)
    Lfull = lattice_points(la, lo, hi, margin=15)
    look = DiscLookup(Lfull, sig_pt)

    subs = [('alignment', cfg.mesh_v if hasattr(cfg, 'mesh_v') else 430, cfg.drift_v if hasattr(cfg, 'drift_v') else 730,
             cfg.sub_run_dir)]
    for d in sorted(glob.glob(os.path.join(cfg.run_dir, '*_scan_*'))):
        nm = os.path.basename(d)
        try:
            a, b = nm.split('_')[-2:]
            subs.append(('drift' if nm.startswith('drift') else 'mesh', int(a), int(b), d))
        except ValueError:
            continue
    store = {}
    for kind, mesh, drift, d in subs:
        if not glob.glob(os.path.join(d, 'combined_hits_root', '*.root')):
            continue
        s = point(cfg, ct, d, T, z, Lfull, bigp)
        s.update(kind=kind, mesh=mesh, drift=drift, field=(drift - mesh) / (cfg.DRIFT_GAP_MM / 10))
        store[(kind, mesh, drift)] = s
        kv = s['th'] < 6
        print(f'{kind:9s} {mesh}/{drift}: interior {len(kv):5d}, vertical {kv.sum():4d}: no-hit d<0.5 '
              f'{100*s["nohit"][kv & (s["dl"] < 0.5)].mean():5.1f} %  d>=1 {100*s["nohit"][kv & (s["dl"] >= 1)].mean():5.1f} %',
              flush=True)
    with open(f'{out}/store.pkl', 'wb') as f:
        pickle.dump(store, f)

    # origin cross-check: Fourier phase at the lattice's reciprocal vectors on pooled vertical tracks
    K = 2 * np.pi * np.linalg.inv(np.c_[la['a1'], la['a2']])      # rows k_i, k_i . a_j = 2 pi delta_ij
    pooled = [s for s in store.values() if s['kind'] == 'alignment' or (s['kind'] == 'drift' and 600 < s['field'] < 1200)]
    P = np.vstack([s['xp'][s['th'] < 8] for s in pooled]) - la['origin']
    nh = np.concatenate([s['nohit'][s['th'] < 8] for s in pooled])
    p0 = nh.mean(); sig = np.sqrt(len(P) * p0 * (1 - p0))
    for nm, k in [('k1', K[0]), ('k2', K[1]), ('k1-k2', K[0] - K[1])]:
        ph = np.exp(1j * (P @ k)); F = ph[nh].sum() - p0 * ph.sum()
        L_.append(f'origin check ({len(P)} pooled vertical tracks): {nm} |F|/sigma {abs(F)/sig:.1f}, phase '
                  f'{np.degrees(np.angle(F)):+.0f} deg (0 = same lattice origin as {lat_key})')

    rows = []

    def one(label, kind, xval, sets):
        for lo_, hi_ in TH:
            N, Y, allp = profile(sets, lo_, hi_)
            r, rlo, rhi = look.fit(Y, N) if lo_ == 0 else (np.nan, np.nan, np.nan)
            rows.append(dict(set=label, kind=kind, x=xval, th_lo=lo_, th_hi=hi_, n=int(N.sum()),
                             nohit_all=100 * allp, nohit_far=100 * Y[2], shadow_pts=100 * (allp - Y[2]),
                             near=100 * Y[0], mid=100 * Y[1], n_near=int(N[0]), r_eff=r, r_lo=rlo, r_hi=rhi))
    for lo_, hi_ in FIELD_GROUPS:
        sets = [s for s in store.values() if s['kind'] == 'drift' and lo_ < s['field'] <= hi_]
        if sets:
            one(f'drift {max(lo_, 0)}-{hi_} V/cm', 'drift', 0.5 * (max(lo_, 0) + hi_), sets)
    for m in sorted({s['mesh'] for s in store.values() if s['kind'] == 'mesh'}):
        one(f'mesh {m} V', 'mesh', m, [s for s in store.values() if s['kind'] == 'mesh' and s['mesh'] == m])
    for key, s in store.items():
        if s['kind'] == 'alignment':
            one(f'alignment {key[1]}/{key[2]}', 'alignment', s['field'], [s])
    R = pd.DataFrame(rows); R.to_csv(f'{out}/shadow_vs_field.csv', index=False)
    L_.append('vertical tracks: no-hit all / far (>= 1 mm) / near (< 0.5 mm) [%], integrated shadow [pts], r_eff [mm]')
    for r in R[R.th_lo == 0].itertuples():
        L_.append(f'   {r.set:22s} n {r.n:5d}: all {r.nohit_all:5.2f}  far {r.nohit_far:5.2f}  near {r.near:5.1f} ({r.n_near})  '
                  f'shadow {r.shadow_pts:+5.2f}  r_eff {r.r_eff:.2f} ({r.r_lo:.2f}-{r.r_hi:.2f})')
    L_.append('>= 12 deg tracks: integrated shadow [pts]: ' + ', '.join(
        f'{r.set}: {r.shadow_pts:+.2f}' for r in R[R.th_lo == 12].itertuples()))

    fig, axs = plt.subplots(1, 3, figsize=(18, 4.8))
    v = R[(R.th_lo == 0) & (R.kind == 'drift')]
    axs[0].errorbar(v.x, v.r_eff, [v.r_eff - v.r_lo, v.r_hi - v.r_eff], fmt='o-', color='#ec4647', label='drift scan')
    a_ = R[(R.th_lo == 0) & (R.kind == 'alignment')]
    axs[0].errorbar(a_.x, a_.r_eff, [a_.r_eff - a_.r_lo, a_.r_hi - a_.r_eff], fmt='s', color='#257de0', label='alignment run')
    axs[0].axhline(0.4, color='#9c9c9c', ls=':', label='drawn pillar radius 0.40 mm')
    axs[0].set(xlabel='drift field [V/cm]', ylabel='effective dead radius [mm]', ylim=(0, 0.75),
               title='Small-pillar shadow vs drift field (vertical tracks)')
    axs[0].grid(color=GRID); axs[0].legend(frameon=False, fontsize=8)
    m = R[(R.th_lo == 0) & (R.kind == 'mesh')]
    axs[1].errorbar(m.x, m.r_eff, [m.r_eff - m.r_lo, m.r_hi - m.r_eff], fmt='o-', color='#ec4647')
    axs[1].axhline(0.4, color='#9c9c9c', ls=':')
    axs[1].set(xlabel='mesh voltage [V]', ylabel='effective dead radius [mm]', ylim=(0, 0.75),
               title='Small-pillar shadow vs mesh voltage (gain)')
    axs[1].grid(color=GRID)
    for (lo_, hi_), col in zip(TH, ['#ec4647', '#f2a600', '#257de0']):
        q = R[(R.th_lo == lo_) & (R.kind == 'drift')]
        axs[2].plot(q.x, q.shadow_pts, 'o-', color=col, label=f'{lo_}-{hi_} deg')
    axs[2].axhline(0, color='#9c9c9c', lw=1)
    axs[2].set(xlabel='drift field [V/cm]', ylabel='integrated shadow [points of efficiency]',
               title='Tracks lost to the small pillars vs drift field')
    axs[2].grid(color=GRID); axs[2].legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(f'{out}/shadow_vs_field.png', dpi=130); plt.close(fig)
    with open(f'{out}/summary.txt', 'w') as f:
        f.write('\n'.join(L_) + '\n')
    print('\n'.join(L_)); print('->', out)


if __name__ == '__main__':
    main()
