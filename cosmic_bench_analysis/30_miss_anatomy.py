#!/usr/bin/env python3
"""30_miss_anatomy.py -- where and why the remaining misses happen.

Usage:  python3 30_miss_anatomy.py <run_key> [--no-veto-sparks] [--wf-chunks N]

Works on stage 06's active-area ray list (same selection as the efficiency).
1. Distance map from the Gerber pad tiles (pad frame, 0.25 mm raster): every
   ray gets its distance to the nearest NON-LIVE region and that region's kind
   (dead/shorted pad tile vs outside the instrumented array). Miss rates vs
   that distance show how far the dead-pad / edge effect reaches.
2. Pointing resolution of the reference at the P2 plane: the core residual
   sigma minus the pad-size term (pitch/sqrt(12)) in quadrature.
3. Interior misses (far from any boundary): vs M3 track inclination, M3 chi2,
   time in run, per-pad clustering (Poisson test), and -- from decoded_root --
   whether a sub-threshold pulse sits on the pad under the track.
Outputs: 30_miss_anatomy/miss_anatomy<sfx>.png + .csv tables
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy import ndimage
from scipy.spatial import cKDTree
import awkward as ak
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path as MPath

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2_qa_config as qa      # noqa: E402
import p2_mapping as pmap      # noqa: E402
import p2_align as pa          # noqa: E402
import p2_waveforms as pw      # noqa: E402
from M3RefTracking import M3RefTracking  # noqa: E402

C_DEAD, C_EDGE, C_RECO, C_ANY = '#ec4647', '#6a4ce2', '#257de0', '#f47230'
INK, INK2, GRID = '#1f1f1f', '#5c5c5c', '#e6e6e6'
RES = 0.25   # raster [mm]


def rate_err(k, n):
    p = k / np.maximum(n, 1)
    return p, np.sqrt(np.clip(p * (1 - p), 0, None) / np.maximum(n, 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run_key', nargs='?', default=qa.DEFAULT_RUN)
    ap.add_argument('--veto-sparks', action=argparse.BooleanOptionalAction, default=True)
    ap.add_argument('--wf-chunks', type=int, default=2, help='decoded chunks for the waveform check')
    ap.add_argument('--interior', type=float, default=6.0,
                    help='interior = farther than this [mm] from any non-live region')
    args = ap.parse_args()
    cfg = qa.get_config(args.run_key)
    out_dir = cfg.out_dir('30_miss_anatomy')
    sfx = cfg.product_suffix(args.veto_sparks)

    # ---------- geometry: live / dead tiles in the pad frame ----------
    full = pmap.build_channel_table(cfg.run_config_path, cfg.MAP_CSV_PATH, det_type=cfg.DET_TYPE,
                                    det_name=cfg.DET_NAME, strategy_overrides=cfg.STRATEGY_OVERRIDES)
    live = pmap.build_channel_table(cfg.run_config_path, cfg.MAP_CSV_PATH, det_type=cfg.DET_TYPE,
                                    det_name=cfg.DET_NAME, drop_connectors=cfg.DEAD_CONNECTORS,
                                    drop_channels=cfg.DEAD_CHANNELS, cluster_link=cfg.CLUSTER_LINK_MM,
                                    strategy_overrides=cfg.STRATEGY_OVERRIDES)
    # all 1280 pads of the map are geometry; live ones are those in `live`
    mp = pd.read_csv(cfg.MAP_CSV_PATH)
    tab_all = mp[['channel_id', 'pad_cx', 'pad_cy', 'pad_w', 'pad_h', 'pad_angle']].copy()
    tab_all['mapped'] = True
    _, verts_all = pmap.pad_tiles(tab_all)
    verts_all = np.asarray(verts_all, float)
    live_ids = set(live.loc[live['mapped'], 'channel_id'])
    is_live = tab_all['channel_id'].isin(live_ids).to_numpy()
    x0, y0 = verts_all[:, :, 0].min() - 20, verts_all[:, :, 1].min() - 20
    x1, y1 = verts_all[:, :, 0].max() + 20, verts_all[:, :, 1].max() + 20
    nx, ny = int((x1 - x0) / RES) + 1, int((y1 - y0) / RES) + 1
    label = np.zeros((ny, nx), np.int8)            # 0 outside, 1 live, 2 dead
    for v, lv in zip(verts_all, is_live):
        ix0 = max(int((v[:, 0].min() - x0) / RES) - 2, 0); ix1 = min(int((v[:, 0].max() - x0) / RES) + 3, nx)
        iy0 = max(int((v[:, 1].min() - y0) / RES) - 2, 0); iy1 = min(int((v[:, 1].max() - y0) / RES) + 3, ny)
        gx, gy = np.meshgrid(x0 + RES * np.arange(ix0, ix1), y0 + RES * np.arange(iy0, iy1))
        inside = MPath(v).contains_points(np.c_[gx.ravel(), gy.ravel()], radius=0.35).reshape(gx.shape)
        sub = label[iy0:iy1, ix0:ix1]
        sub[inside] = 1 if lv else 2
    # close the tiny inter-pad gaps (median 0.07 mm) so they count as live
    live_m = ndimage.binary_closing(label == 1, iterations=2)
    dead_m = (label == 2) & ~live_m
    nonlive = ~live_m
    dist, (iy, ix) = ndimage.distance_transform_edt(nonlive == 0, sampling=RES, return_indices=True)
    nearest_kind = np.where(dead_m[iy, ix], 2, 0)    # kind of the nearest non-live pixel

    # ---------- rays -> pad frame ----------
    padc = live[live['mapped']].drop_duplicates('channel_id').reset_index(drop=True)
    fp = pd.read_csv(f'{cfg.out_dir("06_efficiency")}/pad_footprint{sfx}.csv')
    P = np.c_[padc.pad_cx, padc.pad_cy, np.ones(len(padc))]
    coef, *_ = np.linalg.lstsq(P, fp[['x', 'y']].to_numpy(), rcond=None)
    M, b0 = coef[:2].T, coef[2]
    Minv = np.linalg.inv(M)
    rays = pd.read_csv(f'{cfg.out_dir("06_efficiency")}/ray_hit_miss_list{sfx}.csv')
    for c in ('in_active', 'has_any', 'within'):
        rays[c] = rays[c].astype(str).str.lower().isin(('true', '1'))
    ra = rays[rays['in_active']].copy().reset_index(drop=True)
    q = (ra[['x', 'y']].to_numpy() - b0) @ Minv.T
    ra['px'], ra['py'] = q[:, 0], q[:, 1]
    gi = np.clip(((ra.px - x0) / RES).round().astype(int), 0, nx - 1)
    gj = np.clip(((ra.py - y0) / RES).round().astype(int), 0, ny - 1)
    ra['d_bound'] = dist[gj, gi]
    ra['bound'] = np.where(nearest_kind[gj, gi] == 2, 'dead pad', 'outer edge')
    ra['miss'] = ~ra['within']
    ra['nohit'] = ~ra['has_any']
    N = len(ra)

    # ---------- pointing resolution ----------
    core = ra[ra.within & ra.resid.notna()]
    pitch = np.median(cKDTree(padc[['pad_cx', 'pad_cy']].to_numpy()).query(
        padc[['pad_cx', 'pad_cy']].to_numpy(), k=2)[0][:, 1])
    # 2-D radial residual -> per-axis sigma via Rayleigh median: med = sigma*sqrt(2 ln2)
    s_tot = np.median(core.resid) / np.sqrt(2 * np.log(2))
    s_pad = pitch / np.sqrt(12)
    s_point = np.sqrt(max(s_tot ** 2 - s_pad ** 2, 0))
    print(f'core residual sigma/axis {s_tot:.2f} mm; pad term {s_pad:.2f} mm (pitch {pitch:.2f}); '
          f'=> reference+alignment pointing ~{s_point:.2f} mm per axis')

    # ---------- miss rate vs distance to boundary ----------
    edges = np.array([0, 1, 2, 3, 4, 5, 6, 8, 10, 15, 100])
    rows = []
    for kind in ('dead pad', 'outer edge'):
        s = ra[ra.bound == kind]
        for lo, hi in zip(edges[:-1], edges[1:]):
            t = s[(s.d_bound >= lo) & (s.d_bound < hi)]
            pm, em = rate_err(t.miss.sum(), len(t))
            pn, en = rate_err(t.nohit.sum(), len(t))
            rows.append(dict(kind=kind, lo=lo, hi=hi, n=len(t), miss=pm, miss_err=em, nohit=pn, nohit_err=en))
    prof = pd.DataFrame(rows)
    prof.to_csv(f'{out_dir}/miss_vs_boundary{sfx}.csv', index=False)
    print(prof.round(4).to_string(index=False))

    interior = ra[ra.d_bound >= args.interior].copy()
    base_miss = interior.miss.mean(); base_nohit = interior.nohit.mean()
    excess = {}
    for kind in ('dead pad', 'outer edge'):
        s = ra[(ra.bound == kind) & (ra.d_bound < args.interior)]
        excess[kind] = (len(s), (s.miss.sum() - base_miss * len(s)) / N, (s.nohit.sum() - base_nohit * len(s)) / N)
    print(f'\ninterior (>= {args.interior:g} mm from any non-live region): n {len(interior)} '
          f'({100*len(interior)/N:.1f}% of rays)  miss {100*base_miss:.2f}%  no-hit {100*base_nohit:.2f}%')
    for k, (n, em, en) in excess.items():
        print(f'  band < {args.interior:g} mm of {k}: n {n}, EXCESS miss {100*em:.2f}% of all rays '
              f'(no-hit part {100*en:.2f}%)')
    for guard in (2, 3, 4, 6):
        k_ = ra[~((ra.bound == 'dead pad') & (ra.d_bound < guard))]
        print(f'  efficiency with a {guard} mm guard around dead pads: {100*k_.within.mean():.2f}% '
              f'(keeps {100*len(k_)/N:.1f}% of rays)')

    # ---------- interior misses: angle, chi2, time, pad clustering ----------
    m3 = M3RefTracking(os.path.join(cfg.m3_tracking_dir, ''), single_track=True,
                       chi2_cut=qa.M3_CHI2_CUT, min_nclus=qa.M3_MIN_NCLUS)
    rd = m3.ray_data

    def col(n):
        return np.asarray(ak.to_numpy(ak.ravel(rd[n])), dtype=float)
    ang = pd.DataFrame({'eventId': np.asarray(ak.to_numpy(ak.ravel(rd['evn'])), dtype=np.int64),
                        'theta': np.degrees(np.arctan2(np.hypot(col('X_Up') - col('X_Down'),
                                                                col('Y_Up') - col('Y_Down')),
                                                       np.abs(col('Z_Up') - col('Z_Down')))),
                        'chi2': np.maximum(col('Chi2X'), col('Chi2Y'))}).drop_duplicates('eventId')
    interior = interior.merge(ang, on='eventId', how='left')
    tb = []
    for name, colname, bins in (('theta [deg]', 'theta', [0, 10, 20, 30, 40, 50, 90]),
                                ('M3 chi2', 'chi2', [0, 0.5, 1, 2, 3, 5]),
                                ('time [h]', 't_h', None)):
        if colname == 't_h':
            interior['t_h'] = interior['t_sec'] / 3600.0
            bins = list(np.linspace(0, np.ceil(interior.t_h.max()), 9))
        cut = pd.cut(interior[colname], bins)
        g = interior.groupby(cut, observed=False).agg(n=('miss', 'size'), miss=('miss', 'mean'), nohit=('nohit', 'mean'))
        g['var'] = name
        tb.append(g.reset_index().rename(columns={colname: 'bin'}).assign(bin=lambda t: t['bin'].astype(str)))
        print(f'\ninterior misses vs {name}:'); print(g[['n', 'miss', 'nohit']].round(4).to_string())
    pd.concat(tb).to_csv(f'{out_dir}/interior_vs_vars{sfx}.csv', index=False)

    # per-pad clustering of interior no-hits: dispersion vs binomial expectation
    tree = cKDTree(fp[['x', 'y']].to_numpy())
    _, ipad = tree.query(interior[['x', 'y']].to_numpy())
    interior['pad'] = ipad
    pp = interior.groupby('pad').agg(n=('nohit', 'size'), k=('nohit', 'sum'))
    pp = pp[pp.n >= 20]
    p0 = pp.k.sum() / pp.n.sum()
    chi2_disp = (((pp.k - pp.n * p0) ** 2) / (pp.n * p0 * (1 - p0))).sum()
    dof = len(pp) - 1
    worst = pp.assign(rate=pp.k / pp.n).sort_values('rate', ascending=False).head(10)
    print(f'\nper-pad interior no-hit: mean {100*p0:.2f}%, dispersion chi2/dof = {chi2_disp:.0f}/{dof} '
          f'= {chi2_disp/dof:.2f} (1 = uniform/Poisson)')
    padinfo = padc.iloc[worst.index][['feu', 'channel', 'connector_N', 'half']].reset_index(drop=True)
    print(pd.concat([worst.reset_index(drop=True), padinfo], axis=1).round(3).to_string(index=False))

    # ---------- waveform check on interior no-hit events ----------
    wf_rows = []
    dec = pw.list_decoded_files(cfg.decoded_root_dir)
    dec = dec[dec.chunk < args.wf_chunks]
    nh = interior[interior.nohit].copy()
    ok_hit = interior[interior.within].sample(min(3000, int(interior.within.sum())), random_state=1)
    padc_feu, padc_ch = padc.feu.to_numpy(), padc.channel.to_numpy()
    nh['feu'], nh['ch'] = padc_feu[nh.pad], padc_ch[nh.pad]
    ok_hit = ok_hit.assign(pad=tree.query(ok_hit[['x', 'y']].to_numpy())[1])
    ok_hit['feu'], ok_hit['ch'] = padc_feu[ok_hit.pad], padc_ch[ok_hit.pad]
    for _, fr in dec.iterrows():
        feu = int(fr.feu)
        ped = pw.pedestal_from_data(fr.path); rms = pw.noise_rms_from_data(fr.path, ped)
        want = {**{int(e): ('nohit', int(c)) for e, c in zip(nh.eventId[nh.feu == feu], nh.ch[nh.feu == feu])}}
        ctrl = {int(e): ('hit', int(c)) for e, c in zip(ok_hit.eventId[ok_hit.feu == feu], ok_hit.ch[ok_hit.feu == feu])}
        for ev, ts, ft, wf in pw.iter_decoded_events(fr.path, step=2000):
            sel = np.isin(ev, list(want) + list(ctrl))
            if not sel.any():
                continue
            w = pw.subtract_common_noise(wf[sel] - ped[None, :, None])
            for k, e in enumerate(ev[sel]):
                kind, c = want.get(int(e)) or ctrl.get(int(e))
                wmax = w[k, c].max()
                wf_rows.append(dict(kind=kind, feu=feu, ch=c, peak=wmax, sig=wmax / max(rms[c], 1.0),
                                    # random control: a channel 20 away in the same 64-block
                                    ctrl_sig=w[k, (c // 64) * 64 + (c % 64 + 20) % 64].max()
                                    / max(rms[(c // 64) * 64 + (c % 64 + 20) % 64], 1.0)))
    wfd = pd.DataFrame(wf_rows)
    wfd.to_csv(f'{out_dir}/nohit_waveforms{sfx}.csv', index=False)
    if len(wfd):
        a = wfd[wfd.kind == 'nohit']
        print(f'\nwaveforms on the pad under the track, interior no-hit events (n={len(a)}): '
              f'peak > 3 sigma {100*(a.sig > 3).mean():.1f}%, > 4 sigma {100*(a.sig > 4).mean():.1f}%; '
              f'control channel > 3 sigma {100*(a.ctrl_sig > 3).mean():.1f}%')

    # ---------- figure ----------
    plt.rcParams.update({'font.size': 9.5, 'axes.edgecolor': INK2, 'xtick.color': INK2,
                         'ytick.color': INK2, 'axes.labelcolor': INK})
    fig, axs = plt.subplots(2, 3, figsize=(17, 10))
    ax = axs[0, 0]
    misses = ra[ra.miss]
    ax.scatter(ra.x[~ra.miss], ra.y[~ra.miss], s=1, color='#c9c9c9', rasterized=True)
    col_m = np.where(misses.d_bound < args.interior,
                     np.where(misses.bound == 'dead pad', C_DEAD, C_EDGE), INK)
    ax.scatter(misses.x, misses.y, s=4, c=col_m, rasterized=True)
    for c_, lab in ((C_DEAD, f'miss < {args.interior:g} mm from a dead pad'),
                    (C_EDGE, f'miss < {args.interior:g} mm from the outer edge'), (INK, 'interior miss')):
        ax.scatter([], [], s=12, color=c_, label=lab)
    ax.set_aspect('equal'); ax.legend(fontsize=8, frameon=False, loc='lower left')
    ax.set_title('A  Misses by location', loc='left', fontweight='bold', color=INK)
    ax.set_xlabel('M3 X [mm]'); ax.set_ylabel('M3 Y [mm]')

    ax = axs[0, 1]
    for kind, c_ in (('dead pad', C_DEAD), ('outer edge', C_EDGE)):
        s = prof[(prof.kind == kind) & (prof.n > 30)]
        xm = (s.lo + np.minimum(s.hi, 20)) / 2
        ax.errorbar(xm, 100 * s.miss, 100 * s.miss_err, fmt='o-', color=c_, lw=2, ms=5,
                    label=f'next to: {kind}')
    ax.axhline(100 * base_miss, color=INK2, ls='--', lw=1)
    ax.text(19.5, 100 * base_miss, ' interior', va='bottom', ha='right', fontsize=8, color=INK2)
    ax.axvspan(0, 2 * s_point, color='#f2d7d7', alpha=0.5, lw=0)
    ax.text(2 * s_point, ax.get_ylim()[1] * 0.95 if ax.get_ylim()[1] else 1, f' 2 x pointing sigma\n ({2*s_point:.1f} mm)',
            fontsize=8, color=INK2, va='top')
    ax.set_xlim(0, 20); ax.set_xlabel('distance from the track to the nearest non-live region [mm]')
    ax.set_ylabel('miss rate (not reco within R) [%]')
    ax.set_title('B  Miss rate vs distance to a dead pad / the edge', loc='left', fontweight='bold', color=INK)
    ax.legend(frameon=False, fontsize=8.5)

    ax = axs[0, 2]
    ax.axis('off')
    lines = [f'Remaining misses: {100*ra.miss.mean():.2f}% of {N:,} active-area tracks',
             f'  band < {args.interior:g} mm of a dead pad : +{100*excess["dead pad"][1]:.2f}% (excess)',
             f'  band < {args.interior:g} mm of the edge   : +{100*excess["outer edge"][1]:.2f}% (excess)',
             f'  interior rate                     : {100*base_miss:.2f}%  '
             f'(no-hit {100*base_nohit:.2f}%)',
             '',
             f'Pointing (M3 + alignment) at P2: ~{s_point:.1f} mm per axis',
             f'  (core residual {s_tot:.1f} mm, pad term {s_pad:.1f} mm)',
             '',
             f'Interior no-hit per pad: chi2/dof = {chi2_disp/dof:.2f}',
             '  (1 = uniform; >1 = some pads worse)']
    if len(wfd):
        a = wfd[wfd.kind == 'nohit']
        lines += ['', f'Interior no-hit waveforms (n={len(a)}):',
                  f'  pad under track > 3 sigma: {100*(a.sig > 3).mean():.1f}%',
                  f'  control channel > 3 sigma: {100*(a.ctrl_sig > 3).mean():.1f}%']
    for i, l in enumerate(lines):
        ax.text(0, 1 - i * 0.062, l, fontsize=10 if i == 0 else 9.5, color=INK if i == 0 else INK2,
                fontweight='bold' if i == 0 else None, va='top', family='DejaVu Sans Mono')
    ax.set_title('C  Summary', loc='left', fontweight='bold', color=INK)

    for ax, var, xl in ((axs[1, 0], 'theta [deg]', 'M3 track zenith angle [deg]'),
                        (axs[1, 1], 'time [h]', 'time in run [h]')):
        g = [t for t in tb if t['var'].iloc[0] == var][0]
        g = g[g.n > 30]
        xx = np.arange(len(g))
        pm, em = 100 * g.miss, 100 * np.sqrt(g.miss * (1 - g.miss) / g.n)
        pn, en = 100 * g.nohit, 100 * np.sqrt(g.nohit * (1 - g.nohit) / g.n)
        ax.errorbar(xx, pm, em, fmt='o-', color=C_RECO, lw=2, ms=5, label='miss (not reco within R)')
        ax.errorbar(xx, pn, en, fmt='s--', color=C_ANY, lw=2, ms=5, label='no hit at all')
        ax.set_xticks(xx); ax.set_xticklabels(g['bin'], fontsize=8, rotation=20)
        ax.set_xlabel(xl); ax.set_ylabel('interior rate [%]'); ax.legend(frameon=False, fontsize=8.5)
        ax.set_ylim(0, max(pm.max() + em.max(), 1) * 1.3)
        ax.set_title(f'{"D" if "theta" in var else "E"}  Interior misses vs {xl.split(" [")[0].lower()}',
                     loc='left', fontweight='bold', color=INK)

    ax = axs[1, 2]
    if len(wfd):
        bins = np.linspace(-2, 15, 52)
        for kind, c_, lab in (('nohit', C_DEAD, 'no-hit: pad under track'),
                              ('hit', C_RECO, 'reconstructed: pad under track')):
            s = wfd[wfd.kind == kind]
            ax.hist(np.clip(s.sig, -2, 15), bins, density=True, histtype='step', lw=2, color=c_, label=f'{lab} (n={len(s)})')
        s = wfd[wfd.kind == 'nohit']
        ax.hist(np.clip(s.ctrl_sig, -2, 15), bins, density=True, histtype='step', lw=1.5, ls='--', color=INK2,
                label='no-hit: control channel')
        ax.set_xlabel('max waveform sample / channel noise RMS (clipped at 15)')
        ax.set_ylabel('density'); ax.legend(frameon=False, fontsize=8)
        ax.set_title('F  Is there a sub-threshold pulse under the track?', loc='left', fontweight='bold', color=INK)
    for a_ in axs.ravel():
        for s_ in ('top', 'right'):
            a_.spines[s_].set_visible(False)
        if a_ is not axs[0, 2] and a_ is not axs[0, 0]:
            a_.grid(True, color=GRID, lw=0.8); a_.set_axisbelow(True)
    fig.suptitle(f'{cfg.DET_NAME} anatomy of the remaining misses — {cfg.RUN}/{cfg.SUB_RUN}',
                 x=0.01, ha='left', color=INK, fontsize=13)
    fig.tight_layout()
    fig.savefig(f'{out_dir}/miss_anatomy{sfx}.png', dpi=140)
    plt.close(fig)
    print(f'  saved {out_dir}/miss_anatomy{sfx}.png')


if __name__ == '__main__':
    main()
