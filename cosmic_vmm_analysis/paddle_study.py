#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Why the tag distributions of the scintillator-referenced runs are not the paddle geometry: the paddles'
position-dependent efficiency, measured two ways, and the MC closure.

Inputs
  dump_run*.npz   scint_dump.py: the triggered sample (per trigger: paddle pulses on hybrid 0, clean clusters)
  peff_run*.npz   paddle_efficiency.py: chamber-referenced straight tracks (no paddle in the selection) with the
                  paddle pulses found at the chamber time minus the latency
  scintor_run*.json  the tag positions of scint_efficiency.py (data of the closure figure)

Products in OUT_DIR
  trigger_composition.json  which top paddle fired, slopes and tag medians per category
  paddle_images.png         triggered sample: tracks extrapolated to the paddle planes, fraction with 41 / 45
  paddle_efficiency.png     chamber-referenced efficiency maps of 45, 41, 44 at their planes + along profiles
  paddle_efficiency.json    the along / across profiles, mean efficiencies, propagation-delay slopes
  closure.png / closure.json  tag distributions per probe: data, MC with the measured paddle maps, MC with uniform
                            paddles centred on the band

Band coordinates: along = unit vector at the band angle (21.5 deg to x) from the band centre (pad centroid),
across = its normal; the along axis points to the upper-right (wide) end of the read-out wedge.

Usage: paddle_study.py OUT_DIR --dumps DIR --jsons 'scintor_run*.json' [--n 3e7] [--z-top -50] [--z-bottom 1250]
"""
import argparse
import glob
import json
import os

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

import vmm_qa_path                                          # noqa: F401
import trigger_acceptance_mc as T

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

Z = np.array([320., 630., 940.])
PADDLES = ('p0_45', 'p0_41', 'p0_44')
PLANE = {'p0_45': 'top', 'p0_41': 'top', 'p0_44': 'bottom'}


def cat(zs, k):
    return np.concatenate([z[k] for z in zs])


def band_frame(jsons):
    runs = T.load_data(jsons)
    pads = {nm: runs[0][nm + ':pads'] for nm in T.NAMES}
    centre, ang = T.band_axis(np.vstack([pads[nm] for nm in T.NAMES]))
    a = np.radians(ang)
    R = np.c_[[np.cos(a), np.sin(a)], [-np.sin(a), np.cos(a)]]       # columns: along, across
    return runs, pads, centre, ang, R


def fit_lines(S, cl):
    """Least-squares line per track through the clean clusters S (N, 3, 2) in band coordinates; (N, 2 [intercept,
    slope], 2 [along, across])."""
    N = len(S)
    out = np.full((N, 2, 2), np.nan)
    for i in np.nonzero(cl.sum(1) >= 2)[0]:
        m = cl[i]
        A = np.c_[np.ones(m.sum()), Z[m]]
        for j in range(2):
            out[i, :, j] = np.linalg.lstsq(A, S[i, m, j], rcond=None)[0]
    return out


def profile(x, hit, edges, nmin=10):
    h = np.histogram(x, edges)[0]
    k = np.histogram(x[hit], edges)[0]
    eff = np.where(h >= nmin, k / np.maximum(h, 1), np.nan)
    err = np.where(h >= nmin, np.sqrt(np.maximum(k * (h - k), 0) / np.maximum(h, 1) ** 3), np.nan)
    return h, eff, err


# ---------------------------------------------------------------- A + B: the triggered sample
def triggered(out, dumps, centre, R, zt, zb):
    zs = [np.load(f) for f in dumps]
    a41 = np.isfinite(cat(zs, 'p0_41_adc'))
    a45 = np.isfinite(cat(zs, 'p0_45_adc'))
    a40 = np.isfinite(cat(zs, 'p0_40_adc'))
    cl = np.c_[cat(zs, 'd3_clean'), cat(zs, 'd2_clean'), cat(zs, 'd1_clean')]
    P = np.stack([np.c_[cat(zs, f'd{d}_cx'), cat(zs, f'd{d}_cy')] for d in (3, 2, 1)], 1)
    S = (P - centre) @ R
    F = fit_lines(S, cl)
    top, bot, slope = F[:, 0] + F[:, 1] * zt, F[:, 0] + F[:, 1] * zb, F[:, 1]
    sel = cl.sum(1) >= 2
    cats = {'41 only': a41 & ~a45, '45 only': a45 & ~a41, 'both': a41 & a45}
    res = {'n_trig': int(len(a41)), 'runs': [os.path.basename(f) for f in dumps],
           'fraction': {k: float(m.mean()) for k, m in cats.items()},
           'hardware_triple_ch40': {'fraction': float(a40.mean()), 'equals_both': float((a40 == (a41 & a45)).mean())},
           'paddle_adc_median': {k: float(np.nanmedian(cat(zs, k + '_adc'))) for k in ('p0_40', 'p0_41', 'p0_44', 'p0_45')},
           'tracks_2clean': int(sel.sum()), 'tracks_3clean': int((cl.sum(1) == 3).sum()), 'per_category': {}}
    # tag medians per probe and category (the scint_efficiency.py tag = clean in the two other chambers)
    for k, m in cats.items():
        s = sel & m
        r = {'n_tracks': int(s.sum()),
             'slope_along_mean': float(np.nanmean(slope[s, 0])), 'slope_along_rms': float(np.nanstd(slope[s, 0])),
             'slope_across_mean': float(np.nanmean(slope[s, 1])),
             'top_plane_along_mean': float(np.nanmean(top[s, 0])), 'bottom_plane_along_mean': float(np.nanmean(bot[s, 0]))}
        for i, nm in enumerate(T.NAMES):
            j, kk = [t for t in range(3) if t != i]
            tag = cl[:, j] & cl[:, kk] & m
            f = (Z[i] - Z[j]) / (Z[kk] - Z[j])
            pred = S[:, j] + f * (S[:, kk] - S[:, j])
            r[f'tag_median_along_{nm.split()[0]}'] = float(np.median(pred[tag, 0]))
            r[f'n_tags_{nm.split()[0]}'] = int(tag.sum())
        res['per_category'][k] = r
    json.dump(res, open(os.path.join(out, 'trigger_composition.json'), 'w'), indent=1)
    print('trigger composition:', json.dumps(res['fraction']), 'ch40 =', res['hardware_triple_ch40'])
    for k, r in res['per_category'].items():
        print(f'  {k:8s} n {r["n_tracks"]:5d} slope along {r["slope_along_mean"]:+.3f} | tag medians along det3/det2/det1: '
              f'{r["tag_median_along_det3"]:6.0f} {r["tag_median_along_det2"]:6.0f} {r["tag_median_along_det1"]:6.0f}')
    # figure: paddle images
    eb_al, eb_ac = np.arange(-450, 451, 25), np.arange(-200, 201, 25)
    ext = [eb_al[0], eb_al[-1], eb_ac[0], eb_ac[-1]]

    def H(xy, m):
        return np.histogram2d(xy[m, 0], xy[m, 1], [eb_al, eb_ac])[0]
    Ht = H(top, sel)
    panels = ((Ht, f'top plane z = {zt:g}: all triggered tracks (>= 2 clean chambers)', 'magma_r', None),
              (np.where(Ht >= 8, H(top, sel & a41) / np.maximum(Ht, 1), np.nan), 'fraction with paddle 41', 'viridis', 1),
              (np.where(Ht >= 8, H(top, sel & a45) / np.maximum(Ht, 1), np.nan), 'fraction with paddle 45', 'viridis', 1),
              (H(bot, sel), f'bottom plane z = {zb:g}: all triggered tracks (shaped by paddle 44)', 'magma_r', None),
              (H(top, sel & a41 & ~a45), 'top plane: 41 only', 'magma_r', None),
              (H(top, sel & a45 & ~a41), 'top plane: 45 only', 'magma_r', None))
    fig, axs = plt.subplots(2, 3, figsize=(18, 9), constrained_layout=True)
    for ax, (Zh, t, cm, vm) in zip(axs.ravel(), panels):
        im = ax.imshow(Zh.T, origin='lower', extent=ext, aspect='equal', cmap=cm, vmin=0, vmax=vm, interpolation='nearest')
        ax.set_title(t, fontsize=9)
        ax.set_xlabel('along the band [mm]')
        ax.set_ylabel('across [mm]')
        fig.colorbar(im, ax=ax, shrink=0.7)
        al = np.linspace(-288, 180, 50)
        w = 15 + (al + 288) * (75 / 468)
        ax.plot(al, w, 'c-', lw=0.6)
        ax.plot(al, -w, 'c-', lw=0.6)
    fig.suptitle('triggered sample (scint_dump.py, runs pooled): chamber tracks extrapolated to the paddle planes; '
                 'cyan = read-out wedge', fontsize=11)
    fig.savefig(os.path.join(out, 'paddle_images.png'), dpi=90)
    plt.close(fig)


# ---------------------------------------------------------------- C: chamber-referenced paddle efficiency
def chamber_referenced(out, peffs, centre, R, zt, zb, tol):
    zz = [np.load(f) for f in peffs]
    st = (cat(zz, 'r2') < tol) & ~cat(zz, 'bad')
    S3 = (np.c_[cat(zz, 'x3'), cat(zz, 'y3')] - centre) @ R
    S1 = (np.c_[cat(zz, 'x1'), cat(zz, 'y1')] - centre) @ R
    sl = (S1 - S3) / (Z[2] - Z[0])
    plane = {'top': S3 + sl * (zt - Z[0]), 'bottom': S3 + sl * (zb - Z[0])}
    dt = {k: cat(zz, k) for k in PADDLES + ('p0_40',)}
    hit = {k: np.isfinite(v) for k, v in dt.items()}
    al_e, ac_e = np.arange(-450, 451, 30), np.arange(-200, 201, 20)
    res = {'runs': [os.path.basename(f) for f in peffs], 'n_straight': int(st.sum()),
           'eff_on_straight_tracks': {k: float(h[st].mean()) for k, h in hit.items()},
           'eff_trigger_41or45_and_44': float(((hit['p0_41'] | hit['p0_45']) & hit['p0_44'])[st].mean()),
           'along_edges_mm': al_e.tolist(), 'across_edges_mm': ac_e.tolist(), 'paddles': {}}
    maps = {}
    fig, axs = plt.subplots(3, 3, figsize=(18, 12), constrained_layout=True)
    for col, k in enumerate(PADDLES):
        pl = plane[PLANE[k]]
        inb = np.abs(pl[:, 1]) < 60
        h, eff, err = profile(pl[st & inb, 0], hit[k][st & inb], al_e)
        inb2 = np.abs(pl[:, 0]) < 200
        h2, eff2, err2 = profile(pl[st & inb2, 1], hit[k][st & inb2], ac_e)
        Hn = np.histogram2d(pl[st, 0], pl[st, 1], [al_e, ac_e])[0]
        Hk = np.histogram2d(pl[st & hit[k], 0], pl[st & hit[k], 1], [al_e, ac_e])[0]
        E = np.where(Hn >= 5, Hk / np.maximum(Hn, 1), 0.0)
        maps[k] = (al_e, ac_e, E)
        # propagation delay: median pulse time vs along (both the paddle and 45 fired, 45 as the reference when k != 45)
        c_al = 0.5 * (al_e[1:] + al_e[:-1])
        ref = dt['p0_45'] if k != 'p0_45' else 0.0
        m = st & hit[k] & (np.abs(pl[:, 1]) < 60) & (hit['p0_45'] if k != 'p0_45' else True)
        tmed = np.array([np.median((dt[k] - ref)[m & (pl[:, 0] >= lo) & (pl[:, 0] < hi)])
                         if (m & (pl[:, 0] >= lo) & (pl[:, 0] < hi)).sum() >= 25 else np.nan
                         for lo, hi in zip(al_e[:-1], al_e[1:])])
        ok = np.isfinite(tmed) & (np.abs(c_al) <= 250)
        pfit = np.polyfit(c_al[ok], tmed[ok], 1) if ok.sum() > 3 else [np.nan, np.nan]
        # geometric edges: where the profile crosses half its plateau
        res['paddles'][k] = {'plane': PLANE[k], 'along_profile': [None if np.isnan(x) else round(float(x), 3) for x in eff],
                             'along_n': h.tolist(), 'across_profile': [None if np.isnan(x) else round(float(x), 3) for x in eff2],
                             'across_n': h2.tolist(),
                             'eff_lower_left_half': float(np.nanmean(eff[(c_al > -300) & (c_al < -100)])),
                             'eff_upper_right_half': float(np.nanmean(eff[(c_al > 50) & (c_al < 250)])),
                             'delay_ns_per_m_along': float(pfit[0] * 1000), 'delay_reference': '45' if k != 'p0_45' else 'det2 - latency'}
        ax = axs[0, col]
        im = ax.imshow(np.where(Hn >= 6, E, np.nan).T, origin='lower', extent=[al_e[0], al_e[-1], ac_e[0], ac_e[-1]],
                       aspect='equal', cmap='viridis', vmin=0, vmax=1, interpolation='nearest')
        ax.set_title(f'paddle {k[3:]} efficiency on its plane ({PLANE[k]}, z = {zt if PLANE[k] == "top" else zb:g}), '
                     'chamber straight tracks', fontsize=9)
        ax.set_xlabel('along the band [mm]')
        ax.set_ylabel('across [mm]')
        fig.colorbar(im, ax=ax, shrink=0.7)
        ax = axs[1, col]
        ax.errorbar(c_al, eff, yerr=err, fmt='o', ms=3)
        ax.set_ylim(0, 1.05)
        ax.set_xlabel('along the band at the paddle plane [mm]')
        ax.set_ylabel('efficiency')
        ax.set_title(f'{k[3:]}: along profile, |across| < 60 mm', fontsize=9)
        ax.grid(alpha=.3)
        ax = axs[2, col]
        ax.plot(c_al, tmed, 'o', ms=3)
        ax.plot(c_al[ok], np.polyval(pfit, c_al[ok]), '-', lw=1, color='grey')
        ax.set_xlabel('along the band at the paddle plane [mm]')
        ax.set_ylabel(f'median pulse time - {"t(45)" if k != "p0_45" else "(t(det2) - latency)"} [ns]')
        ax.set_title(f'{k[3:]}: pulse time against position, {pfit[0] * 1000:+.1f} ns/m', fontsize=9)
        ax.grid(alpha=.3)
    fig.suptitle(f'paddle efficiency from {int(st.sum())} chamber-only straight tracks (paddle_efficiency.py, '
                 f'runs {", ".join(os.path.basename(f)[8:-4] for f in peffs)})', fontsize=11)
    fig.savefig(os.path.join(out, 'paddle_efficiency.png'), dpi=90)
    plt.close(fig)
    json.dump(res, open(os.path.join(out, 'paddle_efficiency.json'), 'w'), indent=1)
    print('chamber-referenced paddle efficiency on straight tracks:', json.dumps(res['eff_on_straight_tracks']),
          'trigger', round(res['eff_trigger_41or45_and_44'], 3))
    for k in PADDLES:
        r = res['paddles'][k]
        print(f'  {k}: lower-left half {r["eff_lower_left_half"]:.2f}, upper-right half {r["eff_upper_right_half"]:.2f}, '
              f'delay {r["delay_ns_per_m_along"]:+.1f} ns/m (ref {r["delay_reference"]})')
    return maps


# ---------------------------------------------------------------- D: MC closure
def closure(out, runs, pads, centre, ang, R, maps, n, zt, zb, seed):
    allp = pd.read_csv(T.VS.MAP_CSV).set_index('channel_id').sort_index()
    tree = cKDTree(allp[['pad_cx', 'pad_cy']].to_numpy())
    zs = dict(zip(T.NAMES, Z))
    ch, ce = {}, {}
    for nm in T.NAMES:
        d, i = tree.query(pads[nm])
        read = allp.iloc[np.sort(i)]
        ch[nm] = T.Chamber(read, np.ones(len(read), bool))
        ce[nm] = read[['pad_cx', 'pad_cy']].to_numpy()
    box = (centre[0] - 750, centre[0] + 750, centre[1] - 700, centre[1] + 700)
    rng = np.random.default_rng(seed)
    P = T.Pool(rng, n, box, zs, (0., -100., zb - 25, zb + 25), 2.0, ch, ce, 3.5)

    def plane(z):
        return (np.c_[P.x0 + z * P.tx, P.y0 + z * P.ty] - centre) @ R
    tp, bp = plane(zt), plane(zb)

    def eff2d(k, pl):
        al_e, ac_e, E = maps[k]
        i = np.clip(np.searchsorted(al_e, pl[:, 0]) - 1, 0, len(al_e) - 2)
        j = np.clip(np.searchsorted(ac_e, pl[:, 1]) - 1, 0, len(ac_e) - 2)
        o = E[i, j]
        o[(pl[:, 0] < al_e[0]) | (pl[:, 0] >= al_e[-1]) | (pl[:, 1] < ac_e[0]) | (pl[:, 1] >= ac_e[-1])] = 0
        return o
    e41, e45, e44 = eff2d('p0_41', tp), eff2d('p0_45', tp), eff2d('p0_44', bp)
    w_meas = (1 - (1 - e41) * (1 - e45)) * e44
    w_unif = ((np.abs(tp[:, 0]) < 300) & (np.abs(tp[:, 1]) < 100) & (np.abs(bp[:, 0]) < 300) & (np.abs(bp[:, 1]) < 100)).astype(float)
    data = {nm: T.pooled(runs, nm) for nm in T.NAMES}
    bins_al, bins_ac = np.arange(-360, 361, 15), np.arange(-200, 201, 10)
    fig, axs = plt.subplots(2, 3, figsize=(16, 8), constrained_layout=True)
    res = {'n_thrown': n, 'pool': int(len(P.x0)), 'mean_trigger_prob_measured': float(w_meas.mean()), 'per_probe': {}}
    for col, nm in enumerate(T.NAMES):
        ok, xy = P.pred[nm]
        s = (xy[ok] - centre) @ R
        dx, dy = data[nm]
        sd = (np.c_[dx, dy] - centre) @ R
        r = {'n_data': int(len(sd)), 'data_mean_along': float(sd[:, 0].mean()), 'data_median_along': float(np.median(sd[:, 0]))}
        for lab, w in (('measured', w_meas[ok]), ('uniform', w_unif[ok])):
            r[f'mc_{lab}_mean_along'] = float(np.average(s[:, 0], weights=w))
            for j, bins, ax_lab in ((0, bins_al, 'along'), (1, bins_ac, 'across')):
                hd = np.histogram(sd[:, j], bins)[0]
                hm = np.histogram(s[:, j], bins, weights=w)[0]
                hm2 = np.histogram(s[:, j], bins, weights=w ** 2)[0]
                mu = hm * hd.sum() / hm.sum()
                var_mc = hm2 * (hd.sum() / hm.sum()) ** 2
                k = (mu > 2) | (hd > 2)
                r[f'chi2_{ax_lab}_{lab}'] = float((((hd - mu) ** 2) / np.maximum(mu + var_mc, 1))[k].sum())
                r[f'bins_{ax_lab}_{lab}'] = int(k.sum())
        res['per_probe'][nm] = r
        print(f'  {nm}: data mean along {r["data_mean_along"]:6.1f} | MC measured maps {r["mc_measured_mean_along"]:6.1f} '
              f'(chi2 along {r["chi2_along_measured"]:.0f}/{r["bins_along_measured"]}) | uniform paddles '
              f'{r["mc_uniform_mean_along"]:6.1f} (chi2 {r["chi2_along_uniform"]:.0f}/{r["bins_along_uniform"]})')
        for row, (j, bins, lab) in enumerate(((0, bins_al, f'along the band ({ang:.0f} deg to x) [mm]'), (1, bins_ac, 'across the band [mm]'))):
            ax = axs[row, col]
            hd = np.histogram(sd[:, j], bins)[0]
            c = 0.5 * (bins[1:] + bins[:-1])
            ax.errorbar(c, hd, yerr=np.sqrt(hd), fmt='o', ms=3, label=f'data ({len(sd)} tags)')
            for w, colr, ls, lb in ((w_meas[ok], 'k', '-', 'MC: measured paddle efficiency maps, chambers 100 % efficient'),
                                    (w_unif[ok], 'grey', '--', 'MC: uniform 600 x 200 paddles centred on the band')):
                hm = np.histogram(s[:, j], bins, weights=w)[0]
                ax.step(c, hm * hd.sum() / hm.sum(), where='mid', color=colr, ls=ls, lw=1.3, label=lb)
            ax.set_xlabel(lab)
            ax.set_ylabel('tags')
            ax.set_title(nm)
            if row == 0 and col == 0:
                ax.legend(fontsize=7, frameon=False)
    fig.suptitle('tag positions per probe (scint_efficiency.py, runs pooled): data against the geometric MC with the '
                 'chamber-measured paddle efficiencies (no free parameter)', fontsize=11)
    fig.savefig(os.path.join(out, 'closure.png'), dpi=100)
    plt.close(fig)
    json.dump(res, open(os.path.join(out, 'closure.json'), 'w'), indent=1)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('out_dir')
    ap.add_argument('--dumps', required=True, help='directory with dump_run*.npz and peff_run*.npz')
    ap.add_argument('--jsons', required=True, help='glob of the scintor_run*.json')
    ap.add_argument('--n', type=float, default=3e7)
    ap.add_argument('--z-top', type=float, default=-50.0, help='mid-plane of the top paddle pair [mm]')
    ap.add_argument('--z-bottom', type=float, default=1250.0)
    ap.add_argument('--tol', type=float, default=25.0, help='straightness of the chamber-referenced tracks [mm]')
    ap.add_argument('--seed', type=int, default=1)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    runs, pads, centre, ang, R = band_frame(sorted(glob.glob(a.jsons)))
    print(f'band centre ({centre[0]:.1f}, {centre[1]:.1f}), angle {ang:.2f} deg')
    dumps = sorted(glob.glob(os.path.join(a.dumps, 'dump_run*.npz')))
    peffs = sorted(glob.glob(os.path.join(a.dumps, 'peff_run*.npz')))
    triggered(a.out_dir, dumps, centre, R, a.z_top, a.z_bottom)
    maps = chamber_referenced(a.out_dir, peffs, centre, R, a.z_top, a.z_bottom, a.tol)
    closure(a.out_dir, runs, pads, centre, ang, R, maps, a.n, a.z_top, a.z_bottom, a.seed)


if __name__ == '__main__':
    main()
