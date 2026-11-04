#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Geometric expectation of the scintillator-referenced tag statistics: where should the tags of scint_eff_study.py
fall in each probe chamber if all three chambers were 100 % efficient?

Stand (lab, from 30 Oct 2026; z downwards, P2 pad frame in x-y):
  top paddles      two 200 x 600 x 50 mm paddles exactly on top of each other, z = -100 .. 0 (OR in the trigger)
  det3 P2_IN       z = 320
  det2 P2_MID      z = 630
  det1 P2_OUT      z = 940
  bottom paddle    200 x 600 x 50 mm, z = --z-bottom +- 25 (default 1250, ~30 cm below det1)
The paddle position and long-axis angle in x-y are not recorded: they are FITTED to the measured tags (pooled
scint_efficiency.py jsons, tags_xy of the three probes), starting from "aligned with the read-out band" (user).

Trigger (as the scintor jsons): (top 0:41 OR top 0:45) AND bottom 0:44. A paddle is crossed when the straight track
is inside its rectangle on its top OR bottom face. Muons: cos^n(theta) (default n = 2), flux through a horizontal
plane J (default 1 /cm2/min), thrown uniformly on z = 0 over a box larger than any top-paddle position tried.

Tag for probe P (as scint_efficiency.py): trigger AND a read-out pad hit in each of the two other chambers (the
"clean cluster"; here a 100 % efficient single pad). Predicted point at z_P = the line through the CENTRES of the two
hit pads (the data use cluster centroids), so the pad quantisation of the prediction is reproduced.

The muon pool is computed once (pad hit per chamber, crossing points on the four paddle faces); a paddle placement
is then only a mask on the pool, so the fit (binned Poisson shape likelihood of the three probes' tag positions,
MC normalised to the data per probe) is fast.

Products in OUT_DIR:
  fit.json                    fitted geometry, likelihood, absolute rates (trigger, tags at eff = 1) vs data
  geometry_3d.json            geometry + a sample of triggered muons for the 3D page
  density_<label>.png         data vs MC (eff = 1) tag density per probe, and data/MC
  kernel_<label>.png          the kernel-radius map of scint_eff_study.py for the data and for the MC scaled to the
                              same number of tags (what the statistics plot looks like from geometry alone)
  projections_<label>.png     tag distribution along / across the read-out band, data vs MC (fitted and aligned)

Usage: trigger_acceptance_mc.py OUT_DIR scintor_run*.json [--n 4e7] [--z-bottom 1250] [--maps-runs 187,196]
"""
import argparse
import json
import os
import re

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.spatial import cKDTree

import vmm_qa_path                                          # noqa: F401
import vmm_stations as VS
from cosmic_acceptance_mc import Chamber, CM2_MIN
from scint_eff_study import inside_depth, HALF_DIAG

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

from event_display import style, INK, MUTED, COL  # noqa: E402

DETS = ((3, 'det3 P2_IN'), (2, 'det2 P2_MID'), (1, 'det1 P2_OUT'))
NAMES = [nm for _, nm in DETS]
PAD_L, PAD_W, PAD_T = 600.0, 200.0, 50.0
BIN = 20.0
MC_COL = '#555555'


# ---------------------------------------------------------------- data
def load_data(paths):
    runs = []
    for f in paths:
        E = json.load(open(f))
        run = int(re.findall(r'run_?(\d+)', E['title'] + ' ' + os.path.basename(f))[0])
        bt = E.get('bad_times') or {}
        r = {'run': run, 'title': E['title'], 'n_trig': E['n_trig'],
             'live_s': E['span_s'] - bt.get('vetoed_s', 0.0), 'all_three': E.get('all_three'),
             'clean': {nm: E['chambers'][nm]['clean_cluster_frac'] for nm in NAMES}}
        for nm in NAMES:
            P = E[f'probe_{nm}']
            r[nm] = (np.array(P['tags_xy']['x'], float), np.array(P['tags_xy']['y'], float))
            r[nm + ':d'] = np.array(P['tags_xy']['d'], float)
            r[nm + ':pads'] = np.c_[P['pads']['x'], P['pads']['y']]
        runs.append(r)
    return sorted(runs, key=lambda r: r['run'])


def local_efficiency(runs, nm, pads, R, k=200):
    """Measured local efficiency of chamber nm (as a probe): mean of (probe hit within R) over the k nearest
    pooled tags predicted on its read-out footprint; returns f(xy) -> eff."""
    x = np.concatenate([r[nm][0] for r in runs])
    y = np.concatenate([r[nm][1] for r in runs])
    w = np.concatenate([r[nm + ':d'] for r in runs]) <= R
    on = cKDTree(pads).query(np.c_[x, y])[0] <= HALF_DIAG
    tt, ww = cKDTree(np.c_[x[on], y[on]]), w[on].astype(float)

    def f(xy):
        out = np.empty(len(xy))
        for i0 in range(0, len(xy), 500_000):
            _, idx = tt.query(xy[i0:i0 + 500_000], k=k)
            out[i0:i0 + 500_000] = ww[idx].mean(1)
        return out
    return f


def pooled(runs, nm, which=None):
    rs = [r for r in runs if which is None or r['run'] in which]
    return (np.concatenate([r[nm][0] for r in rs]), np.concatenate([r[nm][1] for r in rs]))


# ---------------------------------------------------------------- geometry
def band_axis(pads):
    """Centre and long-axis angle (deg) of the read-out band (PCA of the pad centres)."""
    c = pads.mean(0)
    w, v = np.linalg.eigh(np.cov((pads - c).T))
    u = v[:, np.argmax(w)]
    if u[0] < 0:
        u = -u
    return c, np.degrees(np.arctan2(u[1], u[0]))


def in_rect(x, y, cx, cy, ang):
    a = np.radians(ang)
    dx, dy = x - cx, y - cy
    lu = dx * np.cos(a) + dy * np.sin(a)
    lv = -dx * np.sin(a) + dy * np.cos(a)
    return (np.abs(lu) <= PAD_L / 2) & (np.abs(lv) <= PAD_W / 2)


def rect_corners(cx, cy, ang):
    a = np.radians(ang)
    u, v = np.array([np.cos(a), np.sin(a)]), np.array([-np.sin(a), np.cos(a)])
    return [list(np.array([cx, cy]) + su * PAD_L / 2 * u + sv * PAD_W / 2 * v)
            for su, sv in ((-1, -1), (1, -1), (1, 1), (-1, 1))]


class Pool:
    """Muons thrown once; everything that does not depend on the paddle placement is stored."""

    def __init__(self, rng, n, box, zs, zfaces, n_exp, chambers, centres, smear, chunk=2_000_000):
        self.box, self.n_thrown, self.zfaces = box, int(n), zfaces
        keep = {k: [] for k in ('x0', 'y0', 'tx', 'ty', 'th') + tuple(f'pad{i}' for i in range(3))}
        left = int(n)
        while left > 0:
            m = min(chunk, left)
            left -= m
            x = rng.uniform(box[0], box[1], m)
            y = rng.uniform(box[2], box[3], m)
            cth = rng.random(m) ** (1.0 / (n_exp + 2))      # flux through a horizontal plane: cos^(n+1) dOmega
            tth = np.sqrt(1 - cth ** 2) / cth
            ph = rng.uniform(0, 2 * np.pi, m)
            tx, ty = tth * np.cos(ph), tth * np.sin(ph)
            pads = [chambers[nm].pad_of(np.c_[x + zs[nm] * tx, y + zs[nm] * ty]) for nm in NAMES]
            ok = (np.vstack(pads) >= 0).sum(0) >= 2          # only these can ever make a tag
            for k, v in (('x0', x), ('y0', y), ('tx', tx), ('ty', ty), ('th', np.degrees(np.arccos(cth)))):
                keep[k].append(v[ok])
            for i in range(3):
                keep[f'pad{i}'].append(pads[i][ok])
        for k, v in keep.items():
            setattr(self, k, np.concatenate(v))
        self.pad = np.vstack([self.pad0, self.pad1, self.pad2])
        self.hit = self.pad >= 0
        self.true = {nm: np.c_[self.x0 + zs[nm] * self.tx, self.y0 + zs[nm] * self.ty] for nm in NAMES}
        # measured point per chamber: the hit pad centre (smear < 0) or the true point + Gaussian smear (cluster
        # centroid resolution; pad pitch / sqrt 12 = 3.4 mm by default)
        meas = {}
        for i, nm in enumerate(NAMES):
            if smear < 0:
                meas[nm] = centres[nm][np.maximum(self.pad[i], 0)]
            else:
                meas[nm] = self.true[nm] + rng.normal(0, smear, self.true[nm].shape)
        # predicted point of each probe from the two other chambers' measured points
        self.pred = {}
        self.weight = {nm: np.ones(len(self.x0)) for nm in NAMES}
        for i, nm in enumerate(NAMES):
            j, k = [t for t in range(3) if t != i]
            a, b = NAMES[j], NAMES[k]
            ok = self.hit[j] & self.hit[k]
            pa = np.where(ok[:, None], meas[a], np.nan)
            pb = np.where(ok[:, None], meas[b], np.nan)
            f = (zs[nm] - zs[a]) / (zs[b] - zs[a])
            self.pred[nm] = (ok, pa + f * (pb - pa))
        # crossing points on the paddle faces
        self.face = {z: (self.x0 + z * self.tx, self.y0 + z * self.ty) for z in zfaces}

    def set_tag_efficiency(self, effmap):
        """Weight of each pool muon as a tag of probe P = local efficiency of the two OTHER chambers at its true
        crossing points (effmap[nm](xy) -> eff); None = all chambers 100 % efficient."""
        for i, nm in enumerate(NAMES):
            if effmap is None:
                self.weight[nm] = np.ones(len(self.x0))
                continue
            j, k = [t for t in range(3) if t != i]
            ok = self.pred[nm][0]
            w = np.zeros(len(self.x0))
            w[ok] = effmap[NAMES[j]](self.true[NAMES[j]][ok]) * effmap[NAMES[k]](self.true[NAMES[k]][ok])
            self.weight[nm] = w

    def area_m2(self):
        return (self.box[1] - self.box[0]) * (self.box[3] - self.box[2]) * 1e-6

    def trigger(self, p):
        zt0, zt1, zb0, zb1 = self.zfaces
        top = in_rect(*self.face[zt0], p[0], p[1], p[2]) | in_rect(*self.face[zt1], p[0], p[1], p[2])
        bot = in_rect(*self.face[zb0], p[3], p[4], p[5]) | in_rect(*self.face[zb1], p[3], p[4], p[5])
        return top & bot


def trigger_rate(rng, p, zfaces, n_exp, J, n=4_000_000):
    """Trigger rate [Hz] of the paddle telescope alone (100 % efficient paddles), thrown over the top paddle box."""
    zt0, zt1, zb0, zb1 = zfaces
    cn = np.array(rect_corners(p[0], p[1], p[2]))
    lo, hi = cn.min(0) - 120, cn.max(0) + 120
    x = rng.uniform(lo[0], hi[0], n)
    y = rng.uniform(lo[1], hi[1], n)
    cth = rng.random(n) ** (1.0 / (n_exp + 2))
    tth = np.sqrt(1 - cth ** 2) / cth
    ph = rng.uniform(0, 2 * np.pi, n)
    tx, ty = tth * np.cos(ph), tth * np.sin(ph)
    top = in_rect(x + zt0 * tx, y + zt0 * ty, p[0], p[1], p[2]) | in_rect(x + zt1 * tx, y + zt1 * ty, p[0], p[1], p[2])
    bot = in_rect(x + zb0 * tx, y + zb0 * ty, p[3], p[4], p[5]) | in_rect(x + zb1 * tx, y + zb1 * ty, p[3], p[4], p[5])
    return J * np.prod(hi - lo) * 1e-6 * np.mean(top & bot)


# ---------------------------------------------------------------- likelihood
class ShapeFit:
    """Binned Poisson shape likelihood of the tag points ON THE PROBE'S READ-OUT FOOTPRINT (the area of the maps; the
    data tails far outside are mis-predicted tags the MC does not model)."""

    def __init__(self, pool, data_xy, edges, pads):
        self.pool, self.edges = pool, edges
        self.D, self.on = {}, {}
        for nm in NAMES:
            tr = cKDTree(pads[nm])
            x, y = data_xy[nm]
            k = tr.query(np.c_[x, y])[0] <= HALF_DIAG
            self.D[nm] = np.histogram2d(x[k], y[k], bins=edges)[0]
            ok, xy = pool.pred[nm]
            on = np.zeros(len(ok), bool)
            on[ok] = tr.query(xy[ok])[0] <= HALF_DIAG
            self.on[nm] = on
        self.calls = 0

    def mc_hist(self, p, mask=None):
        trig = self.pool.trigger(p) if mask is None else mask
        out = {}
        for nm in NAMES:
            ok, xy = self.pool.pred[nm]
            s = trig & ok & self.on[nm]
            out[nm] = np.histogram2d(xy[s, 0], xy[s, 1], bins=self.edges, weights=self.pool.weight[nm][s])[0]
        return out

    def nll(self, p):
        self.calls += 1
        M = self.mc_hist(p)
        tot = 0.0
        for nm in NAMES:
            d, m = self.D[nm], M[nm]
            if m.sum() == 0:
                return 1e12
            mu = d.sum() * (m + 0.5) / (m + 0.5).sum()       # +0.5: no empty MC bin under a data tag
            tot += float((mu - d * np.log(mu)).sum())
        return tot

    def chi2(self, p):
        M = self.mc_hist(p)
        out = {}
        for nm in NAMES:
            d, m = self.D[nm], M[nm]
            mu = d.sum() * m / max(m.sum(), 1)
            s = (d > 0) | (mu > 0)
            out[nm] = {'chi2_neyman_like': float((((d - mu) ** 2) / np.maximum(mu, 1))[s].sum()), 'bins': int(s.sum())}
        return out


# ---------------------------------------------------------------- figures
def grid_hist(x, y, edges):
    return np.histogram2d(x, y, bins=edges)[0]


def density_fig(label, data_xy, mc_xy, pads, edges, path):
    fig, axs = plt.subplots(3, 3, figsize=(16, 10), constrained_layout=True)
    ext = [edges[0][0], edges[0][-1], edges[1][0], edges[1][-1]]
    allp = np.vstack(list(pads.values()))
    view = (allp[:, 0].min() - 80, allp[:, 0].max() + 80, allp[:, 1].min() - 80, allp[:, 1].max() + 80)
    for col, (dd, nm) in enumerate(DETS):
        D = grid_hist(*data_xy[nm], edges)
        M = grid_hist(*mc_xy[nm], edges)
        Ms = M * D.sum() / max(M.sum(), 1)
        vmax = max(D.max(), Ms.max())
        for row, (Z, t, cm, vmin, vm) in enumerate((
                (D, f'{nm}: data, {int(D.sum())} tags', 'magma_r', 0, vmax),
                (Ms, f'{nm}: MC at eff = 1, scaled to the same total', 'magma_r', 0, vmax),
                (np.where(Ms > 3, D / np.maximum(Ms, 1e-9), np.nan), f'{nm}: data / MC (bins with MC > 3)', 'RdBu_r', 0, 2))):
            ax = axs[row, col]
            style(ax)
            ax.grid(False)
            Zp = np.where(Z > 0, Z, np.nan) if row < 2 else Z
            im = ax.imshow(Zp.T, origin='lower', extent=ext, cmap=cm, vmin=vmin, vmax=vm, aspect='equal',
                           interpolation='nearest')
            ax.scatter(pads[nm][:, 0], pads[nm][:, 1], s=0.8, color=MUTED, lw=0, zorder=3)
            ax.set_xlim(view[0], view[1])
            ax.set_ylim(view[2], view[3])
            ax.set_title(t, fontsize=9, color=INK)
            ax.set_xlabel('x [mm]', fontsize=8)
            ax.set_ylabel('y [mm]', fontsize=8)
            if col == 2:
                fig.colorbar(im, ax=axs[row, :], shrink=0.8,
                             label=f'tags per {BIN:g} x {BIN:g} mm' if row < 2 else 'data / MC')
    fig.suptitle(f'{label}: predicted tag point in each probe, data against the geometric expectation '
                 '(trigger x read-out pads of the two tag chambers, all chambers 100 % efficient)', fontsize=11)
    fig.savefig(path, dpi=105)
    plt.close(fig)


def kernel_map(x, y, pads, target, max_kernel=60.0):
    tree = cKDTree(pads)
    gx = np.arange(pads[:, 0].min() - 10, pads[:, 0].max() + 10, 4.0)
    gy = np.arange(pads[:, 1].min() - 10, pads[:, 1].max() + 10, 4.0)
    GX, GY = np.meshgrid(gx, gy, indexing='ij')
    g = np.column_stack([GX.ravel(), GY.ravel()])
    on = tree.query(g)[0] <= HALF_DIAG
    R = np.full(len(g), np.nan)
    if len(x) >= target:
        R[on] = cKDTree(np.column_stack([x, y])).query(g[on], k=target)[0][:, -1]
    return R.reshape(GX.shape), [gx[0] - 2, gx[-1] + 2, gy[0] - 2, gy[-1] + 2]


def kernel_fig(label, data_xy, mcs, pads, margin, target, path, n_full):
    """Kernel radius (radius holding `target` tags) as in scint_eff_study.maps_fig. Row 0: data. Next rows: each MC
    with the data's number of read-out-area tags (shape only). Last row: geometry with the number of tags 100 %
    efficient chambers would give in the same live time (at the measured trigger rate)."""
    nrow = 2 + len(mcs)
    fig, axs = plt.subplots(nrow, 3, figsize=(17, 4.1 * nrow), constrained_layout=True)
    for col, (dd, nm) in enumerate(DETS):
        tree = cKDTree(pads[nm])
        dx, dy = data_xy[nm]
        din = inside_depth(tree, dx, dy) >= margin
        nd = int(din.sum())
        rows = [(dx[din], dy[din], target, f'{nm}: DATA, {nd} tags in the read-out area')]
        for lab, mc in mcs:
            mx, my = map(np.asarray, mc[nm])
            mi = inside_depth(tree, mx, my) >= margin
            rows.append((mx[mi], my[mi], max(int(round(target * mi.sum() / nd)), 1),
                         f'{nm}: MC {lab}, scaled to {nd} tags'))
        mx, my = map(np.asarray, mcs[0][1][nm])
        mi = inside_depth(tree, mx, my) >= margin
        rows.append((mx[mi], my[mi], max(int(round(target * mi.sum() / max(n_full[nm], 1))), 1),
                     f'{nm}: MC {mcs[0][0]}, {n_full[nm]:.0f} tags (all chambers 100 % efficient)'))
        for row, (x, y, t, ttl) in enumerate(rows):
            R, ext = kernel_map(x, y, pads[nm], t)
            ax = axs[row, col]
            style(ax)
            ax.grid(False)
            im = ax.imshow(R.T, origin='lower', extent=ext, cmap=plt.get_cmap('magma_r'), vmin=0, vmax=60,
                           aspect='equal', interpolation='nearest')
            ax.scatter(pads[nm][:, 0], pads[nm][:, 1], s=0.8, color=MUTED, lw=0, zorder=3)
            ax.set_title(ttl, fontsize=9, color=INK)
            ax.set_xlabel('x [mm]', fontsize=8)
            ax.set_ylabel('y [mm]', fontsize=8)
    fig.colorbar(im, ax=axs, shrink=0.4, label=f'kernel radius [mm] (radius holding {target} tags)')
    fig.suptitle(f'{label}: kernel radius of the efficiency maps (scint_eff_study.py, margin {margin:g} mm), '
                 'data against geometry', fontsize=11)
    fig.savefig(path, dpi=95)
    plt.close(fig)


def projections_fig(label, data_xy, curves, centre, ang, path):
    a = np.radians(ang)
    u, v = np.array([np.cos(a), np.sin(a)]), np.array([-np.sin(a), np.cos(a)])
    fig, axs = plt.subplots(2, 3, figsize=(16, 8), constrained_layout=True)
    for col, (dd, nm) in enumerate(DETS):
        for row, (ax_vec, lab, bins) in enumerate(((u, f'along the read-out band ({ang:.0f} deg to x) [mm]',
                                                    np.arange(-360, 361, 15)),
                                                   (v, 'across the band [mm]', np.arange(-200, 201, 10)))):
            ax = axs[row, col]
            style(ax)
            dx, dy = data_xy[nm]
            sd = (np.c_[dx, dy] - centre) @ ax_vec
            hd, _ = np.histogram(sd, bins)
            c = 0.5 * (bins[1:] + bins[:-1])
            ax.errorbar(c, hd, yerr=np.sqrt(hd), fmt='o', ms=3, color=COL[dd], label=f'data ({len(sd)} tags)')
            for lb, mc, colr, st in curves:
                mx, my = mc[nm]
                sm = (np.c_[mx, my] - centre) @ ax_vec
                hm, _ = np.histogram(sm, bins)
                ax.step(c, hm * len(sd) / max(len(sm), 1), where='mid', color=colr, ls=st, lw=1.3, label=lb)
            ax.set_xlabel(lab, fontsize=8)
            ax.set_ylabel('tags', fontsize=8)
            ax.set_title(nm, fontsize=9, color=INK)
            if row == 0:
                ax.legend(fontsize=7, frameon=False, loc='upper left')
    fig.suptitle(f'{label}: tag positions projected on the band axes (all tags, inside and outside the read-out '
                 'area; MC normalised to the data)', fontsize=11)
    fig.savefig(path, dpi=105)
    plt.close(fig)


# ---------------------------------------------------------------- main
NAMES_P = ['top_x', 'top_y', 'top_angle', 'bottom_x', 'bottom_y', 'bottom_angle']


def do_fit(fit, p_aligned, ang, mode):
    def full(q):
        if mode == 'both':
            return q
        if mode == 'common':
            return np.array([q[0], q[1], q[4], q[2], q[3], q[4]])
        return np.array([q[0], q[1], ang, q[2], q[3], ang])
    q0 = {'both': p_aligned, 'common': np.r_[p_aligned[[0, 1, 3, 4]], ang], 'none': p_aligned[[0, 1, 3, 4]]}[mode]
    it = {'both': (0, 1, 3, 4), 'common': (0, 1, 2, 3), 'none': (0, 1, 2, 3)}[mode]
    # coarse grid of the top / bottom displacement ALONG the band (the direction the data constrain least well and
    # where Nelder-Mead stalls on the step-function MC likelihood); the best two grid points seed the simplex
    u = np.array([np.cos(np.radians(ang)), np.sin(np.radians(ang))])
    grid = []
    for st in range(-100, 151, 50):
        for sb in range(-350, 101, 50):
            pp = p_aligned.copy()
            pp[0:2] += st * u
            pp[3:5] += sb * u
            grid.append((fit.nll(pp), st, sb))
    grid.sort()
    print('    grid best (nll, top along, bottom along):', [tuple(round(v, 1) for v in g) for g in grid[:4]])
    best = None
    for _, st, sb in grid[:2]:
        q = q0.copy()
        q[it[0]] += st * u[0]
        q[it[1]] += st * u[1]
        q[it[2]] += sb * u[0]
        q[it[3]] += sb * u[1]
        sh = (st, sb)
        steps = np.full(len(q), 40.0)
        if mode == 'both':
            steps[[2, 5]] = 8.0
        elif mode == 'common':
            steps[4] = 8.0
        simplex = np.vstack([q] + [q + np.eye(len(q))[k] * steps[k] for k in range(len(q))])
        r = minimize(lambda t: fit.nll(full(t)), q, method='Nelder-Mead',
                     options={'initial_simplex': simplex, 'xatol': 1.0, 'fatol': 0.5, 'maxfev': 1500})
        print(f'    start {sh}: nll {r.fun:.1f} at {np.round(full(r.x), 1)} ({r.nfev} calls)')
        if best is None or r.fun < best.fun:
            best = r
    p = full(best.x)
    nll = fit.nll(p)
    scans = {}
    for k, nmk in enumerate(NAMES_P):
        if 'angle' in nmk and mode == 'none':
            continue
        h = 3.0 if 'angle' in nmk else 15.0
        grid = np.arange(-4, 5) * h
        vals = []
        for dlt in grid:
            pp = p.copy()
            pp[k] += dlt
            if mode == 'common' and 'angle' in nmk:
                pp[2] = pp[5] = p[k] + dlt
            vals.append(fit.nll(pp) - nll)
        c2 = np.polyfit(grid, vals, 2)
        scans[nmk] = {'grid': grid.tolist(), 'dnll': list(map(float, vals)),
                      'sigma_parabola': float(1 / np.sqrt(2 * c2[0])) if c2[0] > 0 else float('nan')}
    return p, nll, scans


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('out_dir')
    ap.add_argument('jsons', nargs='+')
    ap.add_argument('--n', type=float, default=4e7, help='muons thrown on z = 0')
    ap.add_argument('--J', type=float, default=1.0, help='flux through a horizontal plane, /cm2/min')
    ap.add_argument('--n-exp', type=float, default=2.0, help='cos^n zenith exponent')
    ap.add_argument('--z', default='320,630,940', help='z of det3, det2, det1 [mm]')
    ap.add_argument('--z-top', type=float, default=0.0, help='lower face of the top paddle pair [mm] (pair = 100 mm)')
    ap.add_argument('--z-bottom', type=float, default=1250.0, help='centre of the bottom paddle [mm]')
    ap.add_argument('--fit-angles', default='both', choices=('both', 'common', 'none'),
                    help='paddle angles: fitted separately, one common angle, or fixed on the band axis')
    ap.add_argument('--smear', type=float, default=12.0 / np.sqrt(12),
                    help='cluster position resolution [mm]; < 0 = use the hit pad centre')
    ap.add_argument('--r', type=float, default=30.0, help='match radius of the measured local efficiency [mm]')
    ap.add_argument('--maps-runs', default='187,196', help='runs pooled for the K = 0 figures')
    ap.add_argument('--margin', type=float, default=10.0)
    ap.add_argument('--target', type=int, default=150)
    ap.add_argument('--seed', type=int, default=1)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    J = a.J * CM2_MIN
    zs = dict(zip(NAMES, (float(v) for v in a.z.split(','))))
    zfaces = (a.z_top, a.z_top - 2 * PAD_T, a.z_bottom - PAD_T / 2, a.z_bottom + PAD_T / 2)

    runs = load_data(a.jsons)
    pads = {nm: runs[0][nm + ':pads'] for nm in NAMES}
    allpads = pd.read_csv(VS.MAP_CSV).set_index('channel_id').sort_index()
    tree = cKDTree(allpads[['pad_cx', 'pad_cy']].to_numpy())
    chambers, centres = {}, {}
    for nm in NAMES:
        d, i = tree.query(pads[nm])
        assert d.max() < 0.05, f'{nm}: json pads not on the pad map'
        read = allpads.iloc[np.sort(i)]
        chambers[nm] = Chamber(read, np.ones(len(read), bool))
        centres[nm] = read[['pad_cx', 'pad_cy']].to_numpy()
    centre, ang = band_axis(np.vstack([pads[nm] for nm in NAMES]))
    print(f'read-out band: centre ({centre[0]:.0f}, {centre[1]:.0f}) mm, long axis {ang:.1f} deg to x')

    box = (centre[0] - 750, centre[0] + 750, centre[1] - 700, centre[1] + 700)
    pool = Pool(rng, a.n, box, zs, zfaces, a.n_exp, chambers, centres, a.smear)
    print(f'pool: {a.n:.2g} thrown on {pool.area_m2():.2f} m2, {len(pool.x0)} hit >= 2 read-out areas')
    effmap = {nm: local_efficiency(runs, nm, pads[nm], a.r) for nm in NAMES}

    data_all = {nm: pooled(runs, nm) for nm in NAMES}
    lo = np.floor((centre - 520) / BIN) * BIN
    edges = (np.arange(lo[0], lo[0] + 1041, BIN), np.arange(lo[1] - 100, lo[1] + 1041, BIN))
    fit = ShapeFit(pool, data_all, edges, pads)
    p_aligned = np.array([centre[0], centre[1], ang, centre[0], centre[1], ang])

    res = {'params': vars(a), 'zfaces_mm': zfaces, 'z_chambers_mm': zs, 'band_centre_mm': centre.tolist(),
           'band_angle_deg': ang, 'paddle_mm': [PAD_L, PAD_W, PAD_T],
           'aligned': dict(zip(NAMES_P, map(float, p_aligned)))}
    models = {}
    for model, em in (('geometry', None), ('geometry_x_tag_eff', effmap)):
        pool.set_tag_efficiency(em)
        print(f'fit, model {model}:')
        p, nll, scans = do_fit(fit, p_aligned, ang, a.fit_angles)
        models[model] = p
        res[model] = {'fit': dict(zip(NAMES_P, map(float, p))), 'nll': nll, 'nll_aligned': fit.nll(p_aligned),
                      'chi2': fit.chi2(p), 'chi2_aligned': fit.chi2(p_aligned), 'scans': scans,
                      'offset_bottom_minus_top_mm': [float(p[3] - p[0]), float(p[4] - p[1])]}
    pool.set_tag_efficiency(None)

    # ---- absolute numbers (all chambers 100 % efficient, paddles of the geometry fit)
    p_fit = models['geometry']
    trig_fit = trigger_rate(rng, p_fit, zfaces, a.n_exp, J)
    mask = pool.trigger(p_fit)
    live = sum(r['live_s'] for r in runs)
    ntr = sum(r['n_trig'] for r in runs)
    res['trigger_rate_hz_eff1'] = trig_fit
    res['trigger_rate_hz_aligned_eff1'] = trigger_rate(rng, p_aligned, zfaces, a.n_exp, J)
    res['data'] = {'runs': [r['run'] for r in runs], 'live_s': live, 'n_trig': ntr, 'trigger_rate_hz': ntr / live}
    sc = J * pool.area_m2() / pool.n_thrown          # Hz per pool muon
    per = {}
    for i, nm in enumerate(NAMES):
        ok = pool.pred[nm][0]
        tag_hz = int((mask & ok).sum()) * sc
        cross = int((mask & pool.hit[i]).sum()) * sc
        nd = sum(len(r[nm][0]) for r in runs)
        cl = float(np.average([r['clean'][nm] for r in runs], weights=[r['n_trig'] for r in runs]))
        per[nm] = {'tags_per_trigger_mc': tag_hz / trig_fit, 'tags_per_trigger_data': nd / ntr,
                   'implied_tag_pair_eff_product': (nd / ntr) / (tag_hz / trig_fit),
                   'readout_crossing_per_trigger_mc': cross / trig_fit, 'clean_cluster_per_trigger_data': cl,
                   'implied_clean_cluster_eff': cl / (cross / trig_fit)}
    res['per_probe'] = per
    res['all_three_per_trigger_mc'] = int((mask & pool.hit.all(0)).sum()) * sc / trig_fit
    res['all_three_per_trigger_data'] = float(np.average([r['all_three'] for r in runs],
                                                         weights=[r['n_trig'] for r in runs]))
    th = pool.th[mask & (pool.hit.sum(0) >= 2)]
    res['zenith_deg_tag_muons'] = {'median': float(np.median(th)), 'p90': float(np.percentile(th, 90))}

    # ---- figures
    def mc_xy(p, em):
        """Tag points of the MC; with em, each tag kept with probability = its tag-chamber efficiency weight."""
        pool.set_tag_efficiency(em)
        m = pool.trigger(p)
        out = {}
        for nm in NAMES:
            ok, xy = pool.pred[nm]
            s = m & ok
            if em is not None:
                s &= rng.random(len(s)) < pool.weight[nm]
            out[nm] = tuple(xy[s].T)
        pool.set_tag_efficiency(None)
        return out
    mgeo = mc_xy(models['geometry'], None)
    mal = mc_xy(p_aligned, None)
    meff = mc_xy(models['geometry_x_tag_eff'], effmap)
    which = {int(t) for t in a.maps_runs.split(',')}
    data_k0 = {nm: pooled(runs, nm, which) for nm in NAMES}
    lab_all = f'runs {runs[0]["run"]}-{runs[-1]["run"]} pooled'
    lab_k0 = 'K+0 (run ' + '+'.join(str(r) for r in sorted(which)) + ')'
    density_fig(lab_all + ', MC = geometry only', data_all, mgeo, pads, edges,
                os.path.join(a.out_dir, 'density_geometry.png'))
    density_fig(lab_all + ', MC = geometry x measured tag-chamber efficiency', data_all, meff, pads, edges,
                os.path.join(a.out_dir, 'density_geometry_x_tag_eff.png'))
    curves = (('MC: geometry, fitted paddles', mgeo, INK, '-'),
              ('MC: geometry, paddles centred on the band', mal, MUTED, '--'),
              ('MC: geometry x tag-chamber efficiency', meff, '#c0392b', '-'))
    projections_fig(lab_all, data_all, curves, centre, ang, os.path.join(a.out_dir, 'projections_all.png'))
    projections_fig(lab_k0, data_k0, curves, centre, ang, os.path.join(a.out_dir, 'projections_K+0.png'))
    live_k0 = sum(r['live_s'] for r in runs if r['run'] in which)
    n_full = {}
    for nm in NAMES:
        mx, my = mgeo[nm]
        fin = np.mean(inside_depth(cKDTree(pads[nm]), np.asarray(mx), np.asarray(my)) >= a.margin)
        n_full[nm] = per[nm]['tags_per_trigger_mc'] * trig_fit * live_k0 * fin
    ntr_k0 = sum(r['n_trig'] for r in runs if r['run'] in which)
    res['K0'] = {'live_s': live_k0, 'n_trig': ntr_k0, 'n_readout_tags_eff1_nominal_flux': n_full,
                 'n_readout_tags_eff1_at_measured_triggers': {nm: v * (ntr_k0 / live_k0) / trig_fit
                                                              for nm, v in n_full.items()}}
    n_meas = res['K0']['n_readout_tags_eff1_at_measured_triggers']
    kernel_fig(lab_k0, data_k0, [('geometry only (eff = 1)', mgeo), ('geometry x tag-chamber efficiency', meff)],
               pads, a.margin, a.target, os.path.join(a.out_dir, 'kernel_K+0.png'), n_meas)

    # ---- geometry + triggered tracks for the 3D page
    sel = np.flatnonzero(mask & pool.hit.all(0))
    sel2 = np.flatnonzero(mask & (pool.hit.sum(0) == 2))
    idx = np.r_[rng.choice(sel, min(250, len(sel)), replace=False),
                rng.choice(sel2, min(80, len(sel2)), replace=False)]

    def padlist(df):
        return {'cx': np.round(df['pad_cx'].to_numpy(), 1).tolist(), 'cy': np.round(df['pad_cy'].to_numpy(), 1).tolist(),
                'w': np.round(df['pad_w'].to_numpy(), 2).tolist(), 'h': np.round(df['pad_h'].to_numpy(), 2).tolist(),
                'a': np.round(df['pad_angle'].to_numpy(), 2).tolist()}
    g3 = {'paddles': {'top': {'corners': rect_corners(*p_fit[:3]), 'z': [zfaces[1], zfaces[0]]},
                      'bottom': {'corners': rect_corners(*p_fit[3:]), 'z': [zfaces[2], zfaces[3]]}},
          'chambers': {nm: {'z': zs[nm]} for nm in NAMES},
          'readout_pads': {nm: padlist(chambers[nm].pads) for nm in NAMES},
          'all_pads': padlist(allpads),
          'tracks': [{'x0': round(float(pool.x0[i]), 1), 'y0': round(float(pool.y0[i]), 1),
                      'tx': round(float(pool.tx[i]), 4), 'ty': round(float(pool.ty[i]), 4),
                      'hits': [bool(h) for h in pool.hit[:, i]]} for i in idx]}
    json.dump(g3, open(os.path.join(a.out_dir, 'geometry_3d.json'), 'w'))
    json.dump(res, open(os.path.join(a.out_dir, 'fit.json'), 'w'), indent=1, default=float)
    short = {k: res[k] for k in ('trigger_rate_hz_eff1', 'trigger_rate_hz_aligned_eff1', 'data',
                                 'all_three_per_trigger_mc', 'all_three_per_trigger_data', 'per_probe', 'K0',
                                 'zenith_deg_tag_muons')}
    for m in models:
        short[m] = {k: v for k, v in res[m].items() if k != 'scans'}
        short[m]['sigma'] = {k: round(v['sigma_parabola'], 1) for k, v in res[m]['scans'].items()}
    print(json.dumps(short, indent=1, default=float))


if __name__ == '__main__':
    main()
