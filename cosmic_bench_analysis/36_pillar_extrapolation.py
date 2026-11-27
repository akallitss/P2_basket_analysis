#!/usr/bin/env python3
"""36_pillar_extrapolation.py -- measure the small-pillar shadow directly and
extrapolate the loss to smaller pillar diameters.

Needs 34_clustering_matching (robust alignment, 0.7 mm pointing, fiducial).

DATA
  * folded-cell map: every interior ray folded into one 4 mm pillar cell
    (position modulo the lattice) -> no-hit map. With 0.7 mm pointing the
    shadow of a 0.8 mm pillar is now visible directly (stage 31 could only see
    it through the angle dependence, with the old 3 mm pointing).
  * no-hit vs distance to the nearest small-pillar centre, by angle class.
  * no-hit vs zenith angle (as 31).

MODEL (electron level, fitted to the three data views jointly)
  primary clusters  Heed mu- 4 GeV in Ar/iC4H10 95/5 (garfield_inputs):
                    Poisson along the path, cluster sizes from Heed
  drift             transverse diffusion sigma = D_T * sqrt(z), D_T from the
                    Magboltz table at the run's drift field
  pillar            an electron landing within r_eff = r_geo + delta of a pillar
                    centre is lost (delta = effective dead ring: field
                    distortion near the insulator, charging)
  detection         signal = sum of Polya(theta=0.5) avalanche gains of the
                    collected electrons; hit if >= thr (mean-gain units)
  pointing          measured position = true + Gauss(sigma_point)
Fit parameters: delta, thr (grid), b = angle/position-independent no-hit (linear).

EXTRAPOLATION: loss vs pillar diameter (0.3-1.0 mm), for the fitted dead ring
held constant AND scaled with the radius, for this run's cosmic angular
distribution and for normal incidence (beam).

Usage: python3 36_pillar_extrapolation.py [run_key]
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2_study_cache as sc     # noqa: E402
import p2_study_common as cm    # noqa: E402
import p2_mapping as pmap       # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
# garfield_inputs is not tracked (ignored) -> in a worktree, read the main checkout's copy
_GI = os.path.join(HERE, 'garfield_inputs')
if not os.path.isfile(os.path.join(_GI, 'gas_table.csv')) and '/.claude/worktrees/' in HERE:
    _GI = os.path.join(HERE.split('/.claude/worktrees/')[0], 'cosmic_bench_analysis', 'garfield_inputs')
INK2, GRID = '#5c5c5c', '#e6e6e6'
C1, C2, C3, C4 = '#9c9c9c', '#257de0', '#ec4647', '#1baf7a'
POLYA_THETA = 0.5
RNG = np.random.default_rng(20261127)


def lattice_frame(pil_m3):
    """Small-pillar lattice in the M3 frame: origin (a pillar), unit vectors, pitch."""
    sm = pil_m3[~pil_m3.big][['x', 'y']].to_numpy()
    from scipy.spatial import cKDTree
    d, j = cKDTree(sm).query(sm, k=2)
    v = sm[j[:, 1]] - sm
    ang = np.degrees(np.arctan2(v[:, 1], v[:, 0])) % 90
    phi = np.radians(np.median(ang))
    pitch = float(np.median(d[:, 1]))
    u1 = np.array([np.cos(phi), np.sin(phi)]); u2 = np.array([-np.sin(phi), np.cos(phi)])
    return sm[0], u1, u2, pitch


def fold(xy, origin, u1, u2, a):
    p = xy - origin
    c1 = (p @ u1 + a / 2) % a - a / 2
    c2 = (p @ u2 + a / 2) % a - a / 2
    return c1, c2


def gas_dt(e_vcm):
    t = pd.read_csv(os.path.join(_GI, 'gas_table.csv'))
    return float(np.interp(e_vcm, t.E_Vcm, t.dt))       # sqrt(cm)


def heed_sizes():
    h = pd.read_csv(os.path.join(_GI, 'heed_clusters.csv.gz'))
    h0 = h[h.angle_deg == 0]
    lam = len(h0) / h0['track'].nunique() / 0.3          # clusters per cm (3 mm Heed slab)
    return lam, h0['ne'].to_numpy()


class ElectronMC:
    """Pre-generated tracks; pillar radius / threshold are applied afterwards."""

    def __init__(self, theta, phi, gap_mm, a, dt_sqrtcm, n_per=1, sigma_point=0.7):
        lam, sizes = heed_sizes()
        th = np.repeat(theta, n_per); ph = np.repeat(phi, n_per)
        n = len(th)
        g = gap_mm
        path_cm = g / 10 / np.cos(np.radians(th))
        ncl = RNG.poisson(lam * path_cm)
        tid = np.repeat(np.arange(n), ncl)
        ne = RNG.choice(sizes, size=len(tid))
        eid = np.repeat(tid, ne)
        z = RNG.uniform(0, g, len(tid))                       # height above the mesh [mm]
        z_e = np.repeat(z, ne)
        x0 = RNG.uniform(-a / 2, a / 2, n); y0 = RNG.uniform(-a / 2, a / 2, n)   # true pos at mid-gap
        t_ = np.tan(np.radians(th))
        dx = (z_e - g / 2) * t_[eid] * np.cos(np.radians(ph[eid]))
        dy = (z_e - g / 2) * t_[eid] * np.sin(np.radians(ph[eid]))
        s = dt_sqrtcm * np.sqrt(z_e / 10) * 10                # mm
        ex = x0[eid] + dx + RNG.normal(0, 1, len(eid)) * s
        ey = y0[eid] + dy + RNG.normal(0, 1, len(eid)) * s
        ex = (ex + a / 2) % a - a / 2; ey = (ey + a / 2) % a - a / 2
        self.r_e = np.hypot(ex, ey)                           # distance to the pillar centre
        self.gain = RNG.gamma(1 + POLYA_THETA, 1 / (1 + POLYA_THETA), len(eid))
        self.eid, self.n = eid, n
        self.theta = th
        mx = x0 + RNG.normal(0, sigma_point, n); my = y0 + RNG.normal(0, sigma_point, n)
        self.mx, self.my = (mx + a / 2) % a - a / 2, (my + a / 2) % a - a / 2
        self.d_meas = np.hypot(self.mx, self.my)
        self.ne_tot = np.bincount(eid, minlength=n)

    def nohit(self, r_eff, thr):
        coll = self.r_e > r_eff
        amp = np.bincount(self.eid[coll], weights=self.gain[coll], minlength=self.n)
        return amp < thr


def binned(x, y, edges):
    c = np.digitize(x, edges) - 1
    n = np.bincount(c[(c >= 0) & (c < len(edges) - 1)], minlength=len(edges) - 1)
    k = np.bincount(c[(c >= 0) & (c < len(edges) - 1)], weights=y[(c >= 0) & (c < len(edges) - 1)],
                    minlength=len(edges) - 1)
    return n, np.divide(k, n, out=np.full(len(n), np.nan), where=n > 0)


def main():
    run_key = sys.argv[1] if len(sys.argv) > 1 else 'det2_long_repos'
    hits, rays, g, cfg = sc.load(run_key)
    out = cfg.out_dir('36_pillar_extrapolation')
    d34 = cfg.out_dir('34_clustering_matching')
    A = np.load(f'{d34}/alignment.npz')
    T = cm.make_T(g, R=A['R'], s=A['s'], pad_mean=A['pad_mean'], m3_mean=A['m3_mean'])
    sig_pt = float(A['sigma_point'])
    ra = pd.read_pickle(f'{d34}/rays_aligned.pkl')
    L = []

    pil = pmap.load_pillars(cfg.MASK_GBR_PATH)
    px, py = T.apply(pil.x.to_numpy(), pil.y.to_numpy())
    pil_m3 = pil.assign(x=px, y=py, r=pil.r * T.s)
    r_geo = float(pil_m3[~pil_m3.big].r.median()); R_big = float(pil_m3[pil_m3.big].r.median())
    origin, u1, u2, a = lattice_frame(pil_m3)
    L.append(f'small pillars: r_geo {r_geo:.3f} mm on a {a:.3f} mm lattice; {int(pil.big.sum())} big '
             f'(R {R_big:.2f} mm); pointing sigma {sig_pt:.2f} mm')

    # interior rays: fiducial, >= 4 mm inside the live area, >= 10 mm from a big pillar
    xy = ra[['x', 'y']].to_numpy()
    big = pil_m3[pil_m3.big][['x', 'y']].to_numpy()
    dbig = np.min(np.hypot(xy[:, None, 0] - big[None, :, 0], xy[:, None, 1] - big[None, :, 1]), 1)
    nohit = ~ra.eventId.isin(hits.eventId).to_numpy()
    inter = (ra.sd.to_numpy() >= 4.0) & (dbig > 10)
    th = ra.theta.to_numpy(); phi = np.degrees(np.arctan2(ra.ty, ra.tx)).to_numpy()
    # offset to the ACTUAL nearest small pillar (the mask is not one global lattice:
    # rows are offset pad by pad, so folding modulo the pitch would smear the shadow)
    from scipy.spatial import cKDTree
    smxy = pil_m3[~pil_m3.big][['x', 'y']].to_numpy()
    dpil, jp = cKDTree(smxy).query(xy)
    off = xy - smxy[jp]
    c1, c2 = off @ u1, off @ u2
    L.append(f'interior rays: {inter.sum():,}; no-hit {100*nohit[inter].mean():.2f} %')

    # ---- data views ----------------------------------------------------------------
    # The no-hit rate depends on angle EVERYWHERE (also >1.5 mm from any pillar), so the
    # pillar term is fitted only to the LOCAL excess near pillars: per angle class a free
    # baseline b_c, plus the electron-level pillar term M(d, theta; delta, thr).
    th_edges = np.array([0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 23, 26, 30])
    d_edges = np.arange(0, 2.01, 0.25)
    cls = [(0, 4), (4, 8), (8, 12), (12, 16), (16, 90)]
    m = inter
    far = m & (dpil >= 1.4)
    n_t, y_t = binned(th[m], nohit[m].astype(float), th_edges)
    n_tf, y_tf = binned(th[far], nohit[far].astype(float), th_edges)
    prof_d = []
    for lo, hi in cls:
        k = m & (th >= lo) & (th < hi)
        prof_d.append(binned(dpil[k], nohit[k].astype(float), d_edges))
    vert = m & (th < 6)
    n_dv, y_dv = binned(dpil[vert], nohit[vert].astype(float), d_edges)
    incl = m & (th >= 12)
    n_di, y_di = binned(dpil[incl], nohit[incl].astype(float), d_edges)
    near_v = vert & (dpil < 0.4); far_v = vert & (dpil >= 1.2)
    L.append(f'vertical (<6 deg) no-hit: within 0.4 mm of a pillar {100*nohit[near_v].mean():.2f} % '
             f'(n={near_v.sum()}), >= 1.2 mm away {100*nohit[far_v].mean():.2f} % (n={far_v.sum()})')
    L.append('no-hit vs angle, all interior / >= 1.4 mm from any pillar [%]: ' + ', '.join(
        f'{lo}-{hi}: {100*y_t[i]:.1f}/{100*y_tf[i]:.1f}' for i, (lo, hi) in
        enumerate(zip(th_edges[:-1], th_edges[1:])) if n_t[i] > 100))

    # ---- model fit -------------------------------------------------------------------
    e_drift = (730 - 430) / (float(g['gap']) / 10)                      # V/cm
    dt = gas_dt(e_drift)
    L.append(f'drift field {e_drift:.0f} V/cm -> D_T {dt:.4f} sqrt(cm): sigma_T over the gap '
             f'{10*dt*np.sqrt(float(g["gap"])/10):.3f} mm')
    mc = ElectronMC(th[m], phi[m], float(g['gap']), a, dt, n_per=4, sigma_point=sig_pt)
    deltas = np.round(np.arange(-0.30, 0.41, 0.02), 3)
    thrs = [1, 2, 3, 5]
    grid = []
    bs = np.linspace(0, 0.10, 401)
    for de in deltas:
        for tr in thrs:
            nh = mc.nohit(max(r_geo + de, 0), tr).astype(float)
            chi, nd, bfit = 0.0, 0, []
            for (lo, hi), (N, Y) in zip(cls, prof_d):
                km = (mc.theta >= lo) & (mc.theta < hi)
                _, M = binned(mc.d_meas[km], nh[km], d_edges)
                ok = (N > 30) & np.isfinite(M) & np.isfinite(Y)
                var = np.clip(Y * (1 - Y), 2e-3, None) / np.maximum(N, 1)
                c_ = [(((Y - (1 - (1 - b_) * (1 - M))) ** 2 / var)[ok]).sum() for b_ in bs]
                j = int(np.argmin(c_)); chi += c_[j]; nd += int(ok.sum()) - 1; bfit.append(bs[j])
            grid.append(dict(delta=de, thr=tr, chi2=chi, ndf=nd - 1,
                             b=';'.join(f'{100*b_:.2f}' for b_ in bfit),
                             loss_mc=float(nh.mean() - mc.nohit(0, tr).mean())))
    G = pd.DataFrame(grid); G.to_csv(f'{out}/fit_grid.csv', index=False)
    bb = G.loc[G.chi2.idxmin()]
    de, tr = float(bb.delta), float(bb.thr)
    b_cls = [float(v) / 100 for v in bb.b.split(';')]
    prof = G.groupby('delta').chi2.min()
    ok_d = prof[prof <= prof.min() + 4].index   # 95% CL band (delta chi2 < 4)
    d_lo, d_hi = float(ok_d.min()), float(ok_d.max())
    G0 = G[G.delta == round(-r_geo, 3)] if (G.delta == round(-r_geo, 3)).any() else G[G.delta == G.delta.min()]
    L.append(f'local fit: dead ring delta = {de:+.2f} mm (95% CL: {d_lo:+.2f} .. {d_hi:+.2f}) -> r_eff '
             f'{r_geo+de:.2f} mm (Gerber {r_geo:.2f}); thr {tr:g} x mean single-e gain; chi2/ndf '
             f'{bb.chi2:.1f}/{int(bb.ndf)}  (no-pillar hypothesis chi2 {G0.chi2.min():.1f}); '
             f'baselines per class {", ".join(f"{100*b_:.2f}" for b_ in b_cls)} %')
    nh_best = mc.nohit(max(r_geo + de, 0), tr).astype(float)
    loss_run = float(nh_best.mean() - mc.nohit(0, tr).mean())
    loss_lo = float(mc.nohit(max(r_geo + d_lo, 0), tr).mean() - mc.nohit(0, tr).mean())
    loss_hi = float(mc.nohit(max(r_geo + d_hi, 0), tr).mean() - mc.nohit(0, tr).mean())
    L.append(f'small-pillar loss for this run\'s angles (model): {100*loss_run:.2f} % '
             f'({100*loss_lo:.2f} .. {100*loss_hi:.2f})')

    # ---- extrapolation ---------------------------------------------------------------
    diam = np.round(np.arange(0.30, 1.001, 0.05), 3)
    mc_beam = ElectronMC(np.zeros(60000), np.zeros(60000), float(g['gap']), a, dt, sigma_point=sig_pt)
    th_c = RNG.choice(th[m], 60000); ph_c = RNG.choice(phi[m], 60000)
    mc_cos = ElectronMC(th_c, ph_c, float(g['gap']), a, dt, sigma_point=sig_pt)
    rows = []
    for D in diam:
        r = D / 2
        for ring, r_e in [('dead ring constant', r + de), ('dead ring scaled', r * (r_geo + de) / r_geo),
                          ('no dead ring', r), ('constant, lo', r + d_lo), ('constant, hi', r + d_hi)]:
            for nm, mm_ in [('cosmic (this run)', mc_cos), ('normal incidence', mc_beam)]:
                base = mm_.nohit(0, tr).mean()
                rows.append(dict(diameter=D, ring=ring, angles=nm,
                                 loss=100 * (mm_.nohit(max(r_e, 0), tr).mean() - base),
                                 geo=100 * np.pi * max(r_e, 0) ** 2 / a ** 2))
    X = pd.DataFrame(rows); X.to_csv(f'{out}/loss_vs_diameter.csv', index=False)
    L.append(f'\nsmall-pillar loss [%] vs diameter (pitch {a:.1f} mm), model with fitted delta/thr:')
    piv = X.pivot_table(index='diameter', columns=['angles', 'ring'], values='loss')
    L.append(piv.round(2).to_string())
    # big pillars: blurred-disc with the measured pointing
    from scipy.stats import rice
    dbig_m = dbig[ra.fid.to_numpy()]
    nb = ~ra.eventId.isin(hits.eventId).to_numpy()[ra.fid.to_numpy()]
    be = np.arange(0, 10.01, 0.5)
    n_b, y_b = binned(dbig_m, nb.astype(float), be)
    xb = 0.5 * (be[1:] + be[:-1])
    from scipy.optimize import curve_fit
    okb = n_b > 20
    bl = float(np.nanmean(y_b[(xb > 7) & okb]))
    f_ = lambda d, R: bl + (1 - bl) * rice.cdf(R / sig_pt, b=d / sig_pt)   # noqa: E731
    (Rb,), cvb = curve_fit(f_, xb[okb], y_b[okb], p0=[R_big + 0.5],
                           sigma=np.sqrt(np.clip(y_b[okb] * (1 - y_b[okb]), 1e-3, None) / n_b[okb]))
    area = ra.fid.sum() / max(1, 1)
    dens = ((dbig_m < 10).sum()) / (len(big) * np.pi * 100)
    loss_big = 100 * len(big) * np.pi * Rb ** 2 * dens / ra.fid.sum()
    L.append(f'\nbig pillars: R_eff {Rb:.2f} +- {np.sqrt(cvb[0,0]):.2f} mm (Gerber {R_big:.2f}) with sigma fixed '
             f'{sig_pt:.2f}; loss {loss_big:.2f} % of fiducial rays')

    # ---- figures ---------------------------------------------------------------------
    fig, axs = plt.subplots(2, 2, figsize=(13, 10))
    ax = axs[0, 0]
    sel = vert
    H, xe, ye = np.histogram2d(c1[sel], c2[sel], bins=20, range=[[-a/2, a/2], [-a/2, a/2]],
                               weights=nohit[sel].astype(float))
    N, _, _ = np.histogram2d(c1[sel], c2[sel], bins=20, range=[[-a/2, a/2], [-a/2, a/2]])
    im = ax.imshow(100 * (H / np.maximum(N, 1)).T, origin='lower', extent=[-a/2, a/2, -a/2, a/2],
                   cmap='Reds', vmin=0)
    ax.add_patch(plt.Circle((0, 0), r_geo, fill=False, color=C2, lw=1.5))
    ax.add_patch(plt.Circle((0, 0), r_geo + de, fill=False, color=C2, lw=1.5, ls='--'))
    plt.colorbar(im, ax=ax, label='no hit [%]')
    ax.set(title=f'Track offset from its nearest pillar, < 6 deg (n={sel.sum():,})', xlabel='lattice u [mm]',
           ylabel='lattice v [mm]')
    ax = axs[0, 1]
    xd = 0.5 * (d_edges[1:] + d_edges[:-1])
    for (lo, hi), (N, Y), b_, col in zip([cls[0], cls[3]], [prof_d[0], prof_d[3]],
                                         [b_cls[0], b_cls[3]], [C3, C2]):
        ax.errorbar(xd, 100 * Y, 100 * np.sqrt(np.clip(Y * (1 - Y), 1e-4, None) / np.maximum(N, 1)),
                    fmt='o', color=col, label=f'data {lo}-{hi} deg')
        km = (mc.theta >= lo) & (mc.theta < hi)
        _, mv = binned(mc.d_meas[km], nh_best[km], d_edges)
        ax.plot(xd, 100 * (1 - (1 - b_) * (1 - mv)), color=col, lw=2,
                label=f'baseline {100*b_:.1f}% + pillar model')
    ax.axvline(r_geo, color=INK2, ls=':', lw=1)
    ax.set(xlabel='measured distance to the nearest small-pillar centre [mm]', ylabel='no hit [%]',
           title=f'Local pillar shadow (pointing {sig_pt:.2f} mm)')
    ax.grid(color=GRID); ax.legend(frameon=False, fontsize=8)
    ax = axs[1, 0]
    xt = 0.5 * (th_edges[1:] + th_edges[:-1])
    okt = n_t > 100
    ax.errorbar(xt[okt], 100 * y_t[okt], 100 * np.sqrt(np.clip(y_t * (1 - y_t), 1e-4, None) / n_t)[okt],
                fmt='o', color=C3, label='data: all interior rays')
    ax.errorbar(xt[okt] + .3, 100 * y_tf[okt],
                100 * np.sqrt(np.clip(y_tf * (1 - y_tf), 1e-4, None) / np.maximum(n_tf, 1))[okt],
                fmt='s', color=C1, label='data: >= 1.4 mm from any small pillar')
    _, mt = binned(mc.theta, nh_best, th_edges)
    ax.plot(xt[okt], 100 * mt[okt], color=C2, lw=2, label='pillar term alone (model)')
    ax.set(xlabel='track zenith angle [deg]', ylabel='no hit [%]',
           title='No hit vs angle: the angle trend is there far from pillars too', ylim=(0, None))
    ax.grid(color=GRID); ax.legend(frameon=False, fontsize=8)
    ax = axs[1, 1]
    for (ring, ls), (ang, col) in [((r_, l_), (a_, c_)) for r_, l_ in
                                   [('dead ring constant', '-'), ('dead ring scaled', '--')]
                                   for a_, c_ in [('cosmic (this run)', C2), ('normal incidence', C3)]]:
        s_ = X[(X.ring == ring) & (X.angles == ang)]
        ax.plot(s_.diameter, s_.loss, ls, color=col, lw=2, label=f'{ang}, {ring}')
        if ring == 'dead ring constant':
            lo_ = X[(X.ring == 'constant, lo') & (X.angles == ang)]
            hi_ = X[(X.ring == 'constant, hi') & (X.angles == ang)]
            ax.fill_between(lo_.diameter, lo_.loss, hi_.loss.to_numpy(), color=col, alpha=.12, lw=0)
    ax.axvline(2 * r_geo, color=INK2, ls=':', lw=1); ax.text(2 * r_geo + .01, ax.get_ylim()[1] * .05,
                                                             'current 0.8 mm', color=INK2, fontsize=8)
    ax.axvline(0.5, color=C4, ls=':', lw=1.5); ax.text(0.51, ax.get_ylim()[1] * .9, '0.5 mm min.',
                                                       color=C4, fontsize=8)
    ax.set(xlabel='small-pillar diameter [mm] (4 mm pitch)', ylabel='tracks lost [%]',
           title='Extrapolated small-pillar loss', ylim=(0, None))
    ax.grid(color=GRID); ax.legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(f'{out}/pillar_extrapolation.png', dpi=140); plt.close(fig)

    with open(f'{out}/summary.txt', 'w') as f:
        f.write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
