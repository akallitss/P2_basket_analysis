#!/usr/bin/env python3
"""41_pillar_lattice.py -- find the REAL pillar lattice of a chamber from the
no-hit tracks themselves (no Gerber assumed), then measure its shadow.

WHY. On det2_long_repos the Gerber (Saclay V1/V2: 0.8 mm pillars on a 4.0 mm
triangular lattice) shows no shadow at all, while the no-hit rate depends
strongly on the track angle (5 % vertical -> 1 % above 16 deg). A model-free
2D Fourier search of the vertical-track no-hit positions finds a triangular
lattice of 3.00 mm pitch (row spacing 2.60 mm) at ~10 sigma in each of its
three directions, in the SAME orientation as the Gerber lattice but NOT at its
pitch. The Gerber and the bulk-lab sheets (Mask_M2_V2) both say 4 mm: the
physical pillar film differs from the drawings. Folding the tracks on the
fitted 3 mm lattice gives a ~20 % no-hit rate within 0.5 mm of a lattice point
for vertical tracks, falling to the ~1 % floor by 1.25 mm, and the shadow
shrinks with angle exactly as the angle trend does. Audit 2026-11-27.

WHAT IT DOES (on the 34 cache + 34 alignment, like 36):
  1. periodicity search: |F(k)|/sigma with F(k) = sum_nohit e^{ik.x} - p0 sum_all e^{ik.x}
     over wavelengths 1.5-6 mm and all directions, vertical interior tracks;
     reports the strongest peaks and whether they form a 60-deg triplet.
  2. lattice: refine the three peaks, solve the origin from their phases, build
     the real-space basis (a1, a2) and the lattice over the footprint; the same
     lattice is written in the PAD frame (lattice_pad.npz) so that other runs of
     the same chamber (stage 42, drift/mesh scans) can reuse it directly: the
     pillars are fixed to the board.
  3. registration check: lattice origin refitted per quadrant (fixed k).
  4. shadow: no-hit vs distance to the nearest lattice point per angle class,
     folded 2D maps, a hard-disc (+) Gaussian-pointing fit of the vertical
     profile (r_eff, sigma_point, baseline) and the matched-cluster charge vs
     distance (is the edge partial?).
  5. big pillars: per-pillar disc fit (centre, radius, edge sigma) -> mask
     registration offset and edge sharpness, independent of the small lattice.

Out: <run>/41_pillar_lattice/{summary.txt, peaks.csv, lattice_m3.npy,
     lattice_pad.npz, profile_by_angle.csv, folded_nohit.png, profile_nohit.png,
     periodicity.png, big_pillars.csv}
Usage: python3 41_pillar_lattice.py [run_key]
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.spatial import cKDTree
from scipy.stats import rice

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2_study_cache as sc     # noqa: E402
import p2_study_common as cm    # noqa: E402
import p2_mapping as pmap       # noqa: E402

GRID = '#e6e6e6'
C1, C2, C3, C4 = '#9c9c9c', '#257de0', '#ec4647', '#1baf7a'
RNG = np.random.default_rng(20261128)
D_EDGES = np.arange(0, 1.76, 0.25)
CLASSES = [(0, 4), (4, 8), (8, 12), (12, 16), (16, 20), (20, 90)]


def inverse(T, xy):
    """M3 -> pad frame (inverse of p2_align.Transform.apply)."""
    return ((np.asarray(xy) - T.m3_mean) / T.s) @ T.R + T.pad_mean


def fourier(P, nh):
    p0 = nh.mean()
    sig = np.sqrt(len(P) * p0 * (1 - p0))

    def F(kx, ky):
        ph = np.exp(1j * (kx * P[:, 0] + ky * P[:, 1]))
        return ph[nh].sum() - p0 * ph.sum()
    return F, sig


def periodicity_scan(P, nh, lam_rng=(1.5, 6.0), dlam=0.02, dang=1.0):
    F, sig = fourier(P, nh)
    rows = []
    for lam in np.arange(lam_rng[0], lam_rng[1] + 1e-9, dlam):
        k = 2 * np.pi / lam
        for a in np.arange(0, 180, dang):
            rows.append((lam, a, abs(F(k * np.cos(np.radians(a)), k * np.sin(np.radians(a)))) / sig))
    return pd.DataFrame(rows, columns=['lam', 'ang', 'z'])


def refine_peak(P, nh, lam0, ang0, dl=0.05, da=2.0):
    F, sig = fourier(P, nh)
    best = (0, lam0, ang0, 0j)
    for lam in np.arange(lam0 - dl, lam0 + dl + 1e-9, 0.002):
        k = 2 * np.pi / lam
        for a in np.arange(ang0 - da, ang0 + da + 1e-9, 0.1):
            f = F(k * np.cos(np.radians(a)), k * np.sin(np.radians(a)))
            if abs(f) > best[0]:
                best = (abs(f), lam, a, f)
    return dict(z=best[0] / sig, lam=best[1], ang=best[2], phase=np.angle(best[3]), F=best[3])


def kvec(pk):
    return 2 * np.pi / pk['lam'] * np.array([np.cos(np.radians(pk['ang'])), np.sin(np.radians(pk['ang']))])


def build_lattice(pk1, pk2, lo, hi, margin=5.0):
    """Origin from the two phases (F ~ e^{i k.x0} for no-hits concentrated at
    x0 + lattice), basis a_j with k_i.a_j = 2 pi delta_ij, points in [lo, hi]."""
    K = np.vstack([kvec(pk1), kvec(pk2)])
    x0 = np.linalg.solve(K, np.array([pk1['phase'], pk2['phase']]))
    Am = 2 * np.pi * np.linalg.inv(K)
    a1, a2 = Am[:, 0], Am[:, 1]
    n = int(np.ceil((hi - lo).max() / min(np.linalg.norm(a1), np.linalg.norm(a2)))) + 5
    ii = np.arange(-n, n + 1)
    I, J = np.meshgrid(ii, ii)
    L = x0 + I.ravel()[:, None] * a1 + J.ravel()[:, None] * a2
    L = L[(L > lo - margin).all(1) & (L < hi + margin).all(1)]
    return L, x0, a1, a2


def binned(x, y, edges):
    c = np.digitize(x, edges) - 1
    ok = (c >= 0) & (c < len(edges) - 1)
    n = np.bincount(c[ok], minlength=len(edges) - 1)
    k = np.bincount(c[ok], weights=y[ok], minlength=len(edges) - 1)
    return n, np.divide(k, n, out=np.full(len(n), np.nan), where=n > 0)


class DiscModel:
    """Hard dead disc of radius r_eff on the lattice, true position smeared by
    a Gaussian pointing of sigma, folded like the data (distance to the nearest
    lattice point of the MEASURED position)."""

    def __init__(self, L, nmc=300000):
        c = L[len(L) // 2]
        self.Lc = L[np.hypot(*(L - c).T) < 15]
        self.t = c + RNG.uniform(-5, 5, (nmc, 2))
        self.dt = cKDTree(self.Lc).query(self.t)[0]
        self.u = RNG.normal(0, 1, (nmc, 2))

    def profile(self, r_eff, sig, b, edges):
        dm = cKDTree(self.Lc).query(self.t + sig * self.u)[0]
        nh = np.where(self.dt < r_eff, 1.0, b)
        n, _ = np.histogram(dm, edges)
        k, _ = np.histogram(dm, edges, weights=nh)
        return k / np.maximum(n, 1)

    def fit(self, Y, N, edges, sig_fixed=None):
        """Grid in (r_eff, sigma); the baseline b enters linearly
        (P = b + (1 - b) M0, M0 the b = 0 profile) and is solved by weighted LS."""
        var = np.clip(Y * (1 - Y), 1e-4, None) / np.maximum(N, 1)
        w = 1 / var
        sigs = [sig_fixed] if sig_fixed else np.arange(0.4, 1.31, 0.1)
        best = (np.inf, 0, 0, 0)
        for sg in sigs:
            dm = cKDTree(self.Lc).query(self.t + sg * self.u)[0]
            n, _ = np.histogram(dm, edges)
            for r in np.arange(0.10, 0.71, 0.02):
                k, _ = np.histogram(dm, edges, weights=(self.dt < r).astype(float))
                M0 = k / np.maximum(n, 1)
                # Y = M0 + b (1 - M0)  ->  b = sum w (Y - M0)(1 - M0) / sum w (1 - M0)^2
                b = float(((w * (Y - M0) * (1 - M0)).sum()) / ((w * (1 - M0) ** 2).sum()))
                b = max(b, 0.0)
                chi = float((w * (Y - (b + (1 - b) * M0)) ** 2).sum())
                if chi < best[0]:
                    best = (chi, r, sg, b)
        chi, r, sg, b = best
        return float(r), float(sg), float(b), float(chi)


def main():
    run_key = sys.argv[1] if len(sys.argv) > 1 else 'det2_long_repos'
    hits, rays, g, cfg = sc.load(run_key)
    out = cfg.out_dir('41_pillar_lattice')
    d34 = cfg.out_dir('34_clustering_matching')
    A = np.load(f'{d34}/alignment.npz')
    T = cm.make_T(g, R=A['R'], s=A['s'], pad_mean=A['pad_mean'], m3_mean=A['m3_mean'])
    ra = pd.read_pickle(f'{d34}/rays_aligned.pkl')
    xy = ra[['x', 'y']].to_numpy(); th = ra.theta.to_numpy()
    nohit = ~ra.eventId.isin(hits.eventId).to_numpy()
    pil = pmap.load_pillars(cfg.MASK_GBR_PATH)
    px, py = T.apply(pil.x.to_numpy(), pil.y.to_numpy())
    big = np.c_[px, py][pil.big.to_numpy()]
    sm = np.c_[px, py][~pil.big.to_numpy()]
    dbig = cKDTree(big).query(xy)[0] if len(big) else np.full(len(xy), np.inf)
    act = (ra.sd.to_numpy() >= 2.0)
    interior = act & (dbig > 10)
    L_ = []

    # ---- 1. periodicity search (vertical interior tracks) --------------------------
    sel = interior & (th < 8)
    P, nh = xy[sel], nohit[sel]
    S = periodicity_scan(P, nh)
    S.to_csv(f'{out}/peaks.csv', index=False)
    n_trials = len(S)
    top = S.sort_values('z', ascending=False)
    L_.append(f'{run_key}: periodicity search on {sel.sum():,} vertical (<8 deg) interior tracks, '
              f'{nh.sum():,} no-hit ({100*nh.mean():.2f} %); {n_trials} (wavelength, direction) trials, '
              f'noise max ~ {np.sqrt(2*np.log(n_trials)):.1f} sigma')
    L_.append('strongest peaks (wavelength [mm], direction [deg], |F|/sigma):')
    # cluster the top entries (one line per distinct peak)
    peaks = []
    for _, r in top.iterrows():
        if any(abs(r.lam - p['lam']) < 0.15 and min(abs(r.ang - p['ang']), 180 - abs(r.ang - p['ang'])) < 6
               for p in peaks):
            continue
        peaks.append(dict(lam=r.lam, ang=r.ang, z=r.z))
        if len(peaks) >= 8:
            break
    for p in peaks:
        L_.append(f'   {p["lam"]:.2f} mm  {p["ang"]:6.1f} deg  {p["z"]:5.1f} sigma')
    # Gerber lattice for reference: its own reciprocal peaks
    if len(sm) > 10:
        d_, j_ = cKDTree(sm).query(sm, k=2)
        v = sm[j_[:, 1]] - sm
        L_.append(f'Gerber ({os.path.basename(cfg.MASK_GBR_PATH)}): nearest-neighbour pitch '
                  f'{np.median(d_[:, 1]):.2f} mm, row spacing {np.median(d_[:, 1])*np.sqrt(3)/2:.2f} mm, '
                  f'row directions {sorted(set(np.round(np.degrees(np.arctan2(v[:, 1], v[:, 0])) % 180 // 5 * 5)))[:6]}')

    # ---- 2. the lattice: a 60-degree triplet among the strong peaks -----------------
    strong = [p for p in peaks if p['z'] > np.sqrt(2 * np.log(n_trials)) + 1.5]
    trip = None
    for i in range(len(strong)):
        for j in range(i + 1, len(strong)):
            pi, pj = strong[i], strong[j]
            if abs(pi['lam'] - pj['lam']) < 0.1 and abs(min(abs(pi['ang'] - pj['ang']),
                                                            180 - abs(pi['ang'] - pj['ang'])) - 60) < 4:
                trip = (pi, pj)
                break
        if trip:
            break
    if trip is None:
        L_.append('NO 60-degree peak triplet above the noise: no periodic pillar shadow detected.')
        with open(f'{out}/summary.txt', 'w') as f:
            f.write('\n'.join(L_) + '\n')
        print('\n'.join(L_))
        return
    pk1 = refine_peak(P, nh, trip[0]['lam'], trip[0]['ang'])
    pk2 = refine_peak(P, nh, trip[1]['lam'], trip[1]['ang'])
    lam = 0.5 * (pk1['lam'] + pk2['lam'])
    pitch = lam * 2 / np.sqrt(3)
    L, x0, a1, a2 = build_lattice(pk1, pk2, xy[interior].min(0), xy[interior].max(0))
    # third direction consistency
    k3 = kvec(pk1) - kvec(pk2)
    F, sig = fourier(P, nh)
    f3 = F(*k3)
    ph3_pred = ((k3 @ x0 + np.pi) % (2 * np.pi)) - np.pi
    L_.append(f'lattice: triangular, row spacing {lam:.3f} mm -> pitch {pitch:.3f} mm; peaks {pk1["z"]:.1f} / '
              f'{pk2["z"]:.1f} sigma at {pk1["ang"]:.1f} / {pk2["ang"]:.1f} deg; third direction {abs(f3)/sig:.1f} '
              f'sigma, phase {np.degrees(np.angle(f3)):.0f} deg (predicted {np.degrees(ph3_pred):.0f})')
    L_.append(f'   basis a1 {np.round(a1, 3)} a2 {np.round(a2, 3)} (|a| {np.linalg.norm(a1):.3f} / '
              f'{np.linalg.norm(a2):.3f} mm, angle {np.degrees(np.arccos(a1 @ a2 / np.linalg.norm(a1) / np.linalg.norm(a2))):.1f} deg); '
              f'origin (M3 frame) {np.round(x0, 3)}')
    np.save(f'{out}/lattice_m3.npy', L)
    Lp = inverse(T, L)
    o_p = inverse(T, x0[None, :])[0]
    a1_p = inverse(T, (x0 + a1)[None, :])[0] - o_p
    a2_p = inverse(T, (x0 + a2)[None, :])[0] - o_p
    np.savez(f'{out}/lattice_pad.npz', origin=o_p, a1=a1_p, a2=a2_p, pitch=pitch, run_key=run_key)
    L_.append(f'   in the PAD frame: origin {np.round(o_p, 3)}, a1 {np.round(a1_p, 3)}, a2 {np.round(a2_p, 3)} '
              f'(lattice_pad.npz; fixed to the board, reusable on other runs of this chamber)')

    # ---- 3. registration per quadrant --------------------------------------------
    K = np.vstack([kvec(pk1), kvec(pk2)])
    xm, ym = np.median(P, 0)
    L_.append('lattice origin refitted per quadrant (fixed k), shift relative to the global origin:')
    for nm, q in [('x- y-', (P[:, 0] < xm) & (P[:, 1] < ym)), ('x- y+', (P[:, 0] < xm) & (P[:, 1] >= ym)),
                  ('x+ y-', (P[:, 0] >= xm) & (P[:, 1] < ym)), ('x+ y+', (P[:, 0] >= xm) & (P[:, 1] >= ym))]:
        Fq, sq = fourier(P[q], nh[q])
        fs = [Fq(*k) for k in K]
        xq = np.linalg.solve(K, np.array([np.angle(f) for f in fs]))
        c = np.linalg.solve(np.c_[a1, a2], xq - x0); c -= np.round(c)
        dq = c[0] * a1 + c[1] * a2
        L_.append(f'   {nm}: n {q.sum():5d}  |F|/sigma {abs(fs[0])/sq:.1f} / {abs(fs[1])/sq:.1f}  '
                  f'shift ({dq[0]:+.2f}, {dq[1]:+.2f}) mm')

    # ---- 4. shadow profiles ----------------------------------------------------------
    d, j = cKDTree(L).query(xy)
    off = xy - L[j]
    rows = []
    L_.append('no-hit [%] vs distance to the nearest lattice point, interior tracks, per angle class:')
    hdr = ' '.join(f'{lo:.2f}-{hi:.2f}' for lo, hi in zip(D_EDGES[:-1], D_EDGES[1:]))
    L_.append(f'   class     {hdr}')
    for lo, hi in CLASSES:
        k = interior & (th >= lo) & (th < hi)
        n, y = binned(d[k], nohit[k].astype(float), D_EDGES)
        rows += [dict(th_lo=lo, th_hi=hi, d_lo=D_EDGES[i], d_hi=D_EDGES[i + 1], n=int(n[i]), p=100 * y[i])
                 for i in range(len(n))]
        L_.append(f'   {lo:2d}-{hi:2d} deg  ' + ' '.join(f'{100*v:9.2f}' for v in y) + f'   (n {k.sum()})')
    pd.DataFrame(rows).to_csv(f'{out}/profile_by_angle.csv', index=False)
    cell = np.sqrt(3) / 2 * pitch ** 2
    vert = interior & (th < 6)
    N, Y = binned(d[vert], nohit[vert].astype(float), D_EDGES)
    dm = DiscModel(L)
    r, s, b, chi = dm.fit(Y, N, D_EDGES)
    r2, _, b2, chi2_ = dm.fit(Y, N, D_EDGES, sig_fixed=float(A['sigma_point']))
    far = interior & (d >= 1.25)
    L_.append(f'hard-disc fit, vertical (<6 deg): r_eff {r:.3f} mm, pointing sigma {s:.2f} mm, baseline '
              f'{100*b:.2f} %, chi2 {chi:.1f}/{len(Y)-3}; with sigma fixed {float(A["sigma_point"]):.2f}: '
              f'r_eff {r2:.3f} mm, baseline {100*b2:.2f} %, chi2 {chi2_:.1f}/{len(Y)-2}')
    L_.append(f'   geometric dead fraction pi r_eff^2 / cell: {100*np.pi*r**2/cell:.2f} % (cell {cell:.2f} mm^2); '
              f'vertical no-hit {100*nohit[vert].mean():.2f} % = baseline + {100*(nohit[vert].mean()-b):.2f} pts')
    L_.append('integrated shadow per angle class = no-hit(all) - no-hit(>= 1.25 mm from a lattice point):')
    for lo, hi in CLASSES:
        k = interior & (th >= lo) & (th < hi)
        L_.append(f'   {lo:2d}-{hi:2d} deg: all {100*nohit[k].mean():5.2f} %, far {100*nohit[k & far].mean():5.2f} %, '
                  f'shadow {100*(nohit[k].mean()-nohit[k & far].mean()):+5.2f} pts (n {k.sum()})')
    # matched-cluster charge vs distance
    cc = f'{cfg.out_dir("34_study_cache")}/clusters_cc20.pkl'
    if os.path.isfile(cc):
        clu = pd.read_pickle(cc)
        cxy = cm.apply_T(T, clu[['x_pad', 'y_pad']].to_numpy())
        dc, jc = cm.nearest(ra.eventId.to_numpy(), xy, cxy, clu.eventId.to_numpy())
        q = np.where(jc >= 0, clu.q.to_numpy()[np.clip(jc, 0, None)], np.nan)
        L_.append('matched-cluster charge (hit tracks, centroid < 20 mm) vs distance to lattice point, '
                  'vertical: median / fraction < 300 ADC')
        for lo, hi in zip(D_EDGES[:-1], D_EDGES[1:]):
            k = vert & (dc < 20) & (d >= lo) & (d < hi)
            L_.append(f'   {lo:.2f}-{hi:.2f}: {np.nanmedian(q[k]):6.0f} ADC   {100*np.mean(q[k] < 300):5.1f} %  (n {k.sum()})')

    # ---- 5. big pillars --------------------------------------------------------------
    brows = []
    if len(big):
        L_.append('big pillars, per-pillar disc fit (theta < 12 deg, within 9 mm): centre offset from the Gerber, '
                  'R_eff, edge sigma (pointing), baseline')
        for i, bc in enumerate(big):
            dd = np.hypot(*(xy - bc).T)
            k = act & (dd < 9) & (th < 12)
            if k.sum() < 30:
                L_.append(f'   big {i} at ({bc[0]:7.1f}, {bc[1]:7.1f}): n {k.sum()} (outside the live area)')
                continue
            Pk, nk = xy[k], nohit[k].astype(float)

            def nll(p):
                cx, cy, R, sg, bb = p
                if R < 0.5 or sg < 0.25 or sg > 4 or bb < 0.001 or bb > 0.3:
                    return 1e9
                rr = np.hypot(Pk[:, 0] - cx, Pk[:, 1] - cy)
                pr = np.clip(bb + (1 - bb) * rice.cdf(R / sg, b=rr / sg), 1e-6, 1 - 1e-6)
                return -(nk * np.log(pr) + (1 - nk) * np.log(1 - pr)).sum()
            best = min((minimize(nll, [bc[0] + dx, bc[1] + dy, 3.0, 0.8, 0.03], method='Nelder-Mead',
                                 options=dict(maxiter=2000, xatol=1e-3, fatol=1e-3))
                        for dx in (-1.5, 0, 1.5) for dy in (-1.5, 0, 1.5)), key=lambda r_: r_.fun)
            cx, cy, R, sg, bb = best.x
            brows.append(dict(i=i, x_gerber=bc[0], y_gerber=bc[1], dx=cx - bc[0], dy=cy - bc[1], R=R, sigma=sg,
                              base=bb, n=int(k.sum())))
            L_.append(f'   big {i} at ({bc[0]:7.1f}, {bc[1]:7.1f}): n {k.sum():3d}, offset ({cx-bc[0]:+.2f}, '
                      f'{cy-bc[1]:+.2f}) mm, R_eff {R:.2f} mm (Gerber {float(pil[pil.big].r.median()):.2f}), '
                      f'edge sigma {sg:.2f} mm, baseline {100*bb:.1f} %')
        if brows:
            B = pd.DataFrame(brows); B.to_csv(f'{out}/big_pillars.csv', index=False)
            L_.append(f'   mean offset ({B.dx.mean():+.2f} +- {B.dx.std()/np.sqrt(len(B)):.2f}, '
                      f'{B.dy.mean():+.2f} +- {B.dy.std()/np.sqrt(len(B)):.2f}) mm: the mask film is shifted '
                      f'relative to the pads; mean R_eff {B.R.mean():.2f} mm')

    # ---- figures ---------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(9, 4.5))
    piv = S.pivot(index='lam', columns='ang', values='z')
    im = ax.imshow(piv.to_numpy(), origin='lower', aspect='auto', cmap='magma',
                   extent=[piv.columns.min(), piv.columns.max(), piv.index.min(), piv.index.max()])
    plt.colorbar(im, ax=ax, label='|F| / sigma')
    ax.set(xlabel='direction of k [deg]', ylabel='wavelength [mm]',
           title=f'{run_key}: periodicity of the vertical-track no-hits (peaks: {lam:.2f} mm rows -> {pitch:.2f} mm pitch)')
    fig.tight_layout(); fig.savefig(f'{out}/periodicity.png', dpi=130); plt.close(fig)

    fig, axs = plt.subplots(1, 3, figsize=(16, 5))
    for ax, (lab, s_) in zip(axs, [('theta < 6 deg', th < 6), ('6-12 deg', (th >= 6) & (th < 12)),
                                   ('theta >= 16 deg', th >= 16)]):
        s_ = s_ & interior
        rng_ = [[-1.8, 1.8], [-1.8, 1.8]]
        H, _, _ = np.histogram2d(off[s_, 0], off[s_, 1], bins=24, range=rng_, weights=nohit[s_].astype(float))
        Nn, _, _ = np.histogram2d(off[s_, 0], off[s_, 1], bins=24, range=rng_)
        im = ax.imshow(100 * (H / np.maximum(Nn, 1)).T, origin='lower', extent=[-1.8, 1.8, -1.8, 1.8],
                       cmap='Reds', vmin=0, vmax=25)
        plt.colorbar(im, ax=ax, label='no hit [%]')
        ax.add_patch(plt.Circle((0, 0), r, fill=False, color='k', lw=1, ls='--'))
        ax.set(title=f'{lab} (n={s_.sum():,})', xlabel='x - nearest lattice point [mm]', ylabel='y [mm]')
    fig.suptitle(f'{run_key}: no-hit rate around the nearest point of the fitted {pitch:.2f} mm lattice '
                 f'(dashed: fitted r_eff {r:.2f} mm)')
    fig.tight_layout(); fig.savefig(f'{out}/folded_nohit.png', dpi=130); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.5, 5))
    xc = 0.5 * (D_EDGES[1:] + D_EDGES[:-1])
    for (lo, hi), col in zip([(0, 6), (6, 12), (12, 16), (16, 90)], [C3, '#f2a600', C4, C2]):
        k = interior & (th >= lo) & (th < hi)
        n_, y_ = binned(d[k], nohit[k].astype(float), D_EDGES)
        ax.errorbar(xc, 100 * y_, 100 * np.sqrt(np.clip(y_ * (1 - y_), 1e-4, None) / np.maximum(n_, 1)),
                    fmt='o-', color=col, label=f'{lo}-{hi} deg')
    ax.plot(xc, 100 * dm.profile(r, s, b, D_EDGES), 'k--',
            label=f'hard disc r_eff {r:.2f} mm, pointing {s:.2f} mm, base {100*b:.1f} %')
    ax.set(xlabel='distance to the nearest lattice point [mm]', ylabel='no hit [%]',
           title=f'{run_key}: pillar shadow on the fitted {pitch:.2f} mm lattice')
    ax.grid(color=GRID); ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(f'{out}/profile_nohit.png', dpi=130); plt.close(fig)

    with open(f'{out}/summary.txt', 'w') as f:
        f.write('\n'.join(L_) + '\n')
    print('\n'.join(L_))
    print('->', out)


if __name__ == '__main__':
    main()
