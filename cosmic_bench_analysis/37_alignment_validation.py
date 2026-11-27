#!/usr/bin/env python3
"""37_alignment_validation.py -- what each alignment / active-area change does,
one step at a time, on ONE run (validation before any pipeline-wide switch).

Steps (cumulative unless marked), all with the stage-06 efficiency definition
(closest CC20 cluster within MATCH_R; active = ACTIVE_R + dead-pad Voronoi cells
+ guards), computed with the SAME p2_align code the pipeline uses:

  A  stage 06 as it was   plain fit, z 240, dead guard 6 mm, ACTIVE_R edge
  B  + trimmed fit        (scale 0.978 -> 1.000)
  C  + plane z 246        (slope fit)
  D  + dead guard 2 mm
  E  + outer-edge guard 2 mm   (= the new det2_long_repos config)
  A2 control: guards 2/2 with the OLD alignment (what the guard change alone does)

Per step: active rays, efficiency (within R / any hit), the accidental rate
(event mixing) and accidental-corrected efficiency, the robust residual
sigma/axis, the transform scale; plus the dead-border and outer-edge miss
profiles and the residual maps for A vs E.

Usage: python3 37_alignment_validation.py [run_key]
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
import p2_align as pa           # noqa: E402

INK2, GRID = '#5c5c5c', '#e6e6e6'
C_OLD, C_NEW, C_CTL = '#9c9c9c', '#257de0', '#ec4647'


def main():
    run_key = sys.argv[1] if len(sys.argv) > 1 else 'det2_long_repos'
    hits, rays, g, cfg = sc.load(run_key)
    out = cfg.out_dir('37_alignment_validation')
    cpath = cfg.out_dir('34_study_cache') + '/clusters_cc20.pkl'
    if os.path.isfile(cpath):
        clu = pd.read_pickle(cpath)
    else:
        clu = cm.cc_clusters(hits, 20.0); clu.to_pickle(cpath)
    ct = sc.channel_table(cfg)
    ev = rays['eventId'].to_numpy()
    cpad = clu[['x_pad', 'y_pad']].to_numpy()
    R = float(g['match_r']); AR = float(g['active_r'])
    hit_any = np.isin(ev, hits['eventId'].to_numpy())
    padc = ct[ct['mapped']].drop_duplicates('channel_id')
    z_old = 240.0
    z_new = float(cfg.det_plane_z())
    trim = cfg.FIT_TRIM_MM or 15.0

    def align(z, trimmed):
        x, y = pa.project_to_z(rays, z); xy = np.column_stack([x, y])
        # the 06 fit: (track, max-charge cluster) pairs in the fiducial window
        cen = clu.loc[clu.groupby('eventId')['q'].idxmax()]
        m = pd.DataFrame({'eventId': ev, 'x': x, 'y': y}).merge(cen, on='eventId')
        m = m[(m.x.abs() < 300) & (m.y.abs() < 300)]
        T = pa.fit_transform(m.x, m.y, m.x_pad, m.y_pad, trim_mm=trim if trimmed else None)
        return xy, T

    def evaluate(xy, T, dead_guard, edge_guard):
        ct.attrs['dead_guard'] = dead_guard
        ct.attrs['edge_guard'] = edge_guard
        pxy = np.column_stack(T.apply(padc['pad_cx'].to_numpy(), padc['pad_cy'].to_numpy()))
        from scipy.spatial import cKDTree
        nn, _ = cKDTree(pxy).query(xy)
        nn = pa.mask_dead_cells(nn, xy, ct, T)
        act = nn <= AR
        cxy = cm.apply_T(T, cpad)
        d, j = cm.nearest(ev, xy, cxy, clu['eventId'].to_numpy())
        hit = d <= R
        # accidentals: pair each ray with the clusters of an event 1000 triggers later
        ue = np.unique(ev); mp = dict(zip(ue, np.roll(ue, -1000)))
        dm, _ = cm.nearest(np.array([mp[e] for e in ev]), xy, cxy, clu['eventId'].to_numpy())
        acc = (dm <= R)[act].mean()
        e_ = hit[act].mean()
        k = act & (d < 15)
        r = cxy[j[k]] - xy[k]
        sig = 1.4826 * np.median(np.abs(r - np.median(r, 0)), 0)
        return dict(n=int(act.sum()), eff=100 * e_, err=cm.binom_err(hit[act].sum(), act.sum()),
                    any=100 * hit_any[act].mean(), acc=100 * acc, eff_corr=100 * (e_ - acc) / (1 - acc),
                    sig_x=sig[0], sig_y=sig[1], scale=T.s), act, hit

    steps = [('A', 'stage 06 as it was (plain fit, z 240, guard 6, ACTIVE_R edge)', z_old, False, 6.0, None),
             ('B', '+ trimmed alignment fit', z_old, True, 6.0, None),
             ('C', f'+ plane z {z_new:g}', z_new, True, 6.0, None),
             ('D', '+ dead-pad guard 6 -> 2 mm', z_new, True, 2.0, None),
             ('E', '+ outer-edge guard 2 mm (new config)', z_new, True, 2.0, 2.0),
             ('A2', 'control: old alignment with guards 2/2', z_old, False, 2.0, 2.0)]
    rows, keep = [], {}
    for key, lab, z, tr, dg, eg in steps:
        xy, T = align(z, tr)
        res, act, hit = evaluate(xy, T, dg, eg)
        rows.append(dict(step=key, change=lab, **res))
        keep[key] = (xy, T, act, hit)
        print(f'{key:3s} {lab:62s} N {res["n"]:6d}  eff {res["eff"]:.2f} +- {res["err"]:.2f}  '
              f'acc {res["acc"]:.2f}  corr {res["eff_corr"]:.2f}  any {res["any"]:.2f}  '
              f'sigma {res["sig_x"]:.2f}/{res["sig_y"]:.2f}  s {res["scale"]:.4f}')
    tab = pd.DataFrame(rows)
    tab.to_csv(f'{out}/steps.csv', index=False)

    # ---- profiles A vs E: miss vs distance to dead tiles, and to the outer edge ----
    fig, axs = plt.subplots(1, 3, figsize=(18, 4.8))
    for key, col in [('A', C_OLD), ('E', C_NEW)]:
        xy, T, act, hit = keep[key]
        ct.attrs['dead_guard'] = 0.0; ct.attrs['edge_guard'] = None
        dd = pa.dist_to_dead_tiles(xy, ct, T)
        ed = pa.outer_edge_distance(xy, ct, T)
        from scipy.spatial import cKDTree
        pxy = np.column_stack(T.apply(padc['pad_cx'].to_numpy(), padc['pad_cy'].to_numpy()))
        nn, _ = cKDTree(pxy).query(xy)
        cell = pa.mask_dead_cells(nn, xy, ct, T) <= AR      # live Voronoi cell, no guards
        interior = 1 - hit[cell & (dd > 10) & (ed > 10)].mean()
        b = np.arange(0, 8.01, 0.5); xb = b[:-1] + 0.25
        sel = cell & (ed >= 3)
        mr = [100 * (1 - hit[sel & (dd >= lo) & (dd < lo + .5)].mean()) for lo in b[:-1]]
        axs[0].plot(xb, mr, '-o', color=col, lw=2, ms=4,
                    label=f'{key}: {"stage 06 alignment" if key == "A" else "corrected alignment"}')
        axs[0].axhline(100 * interior, color=col, ls=':', lw=1)
        be = np.arange(-4, 8.01, 0.5); xe = be[:-1] + 0.25
        sel2 = (dd > 10) & (nn <= 12)
        me = [100 * (1 - hit[sel2 & (ed >= lo) & (ed < lo + .5)].mean())
              if (sel2 & (ed >= lo) & (ed < lo + .5)).sum() > 30 else np.nan for lo in be[:-1]]
        axs[1].plot(xe, me, '-o', color=col, lw=2, ms=4, label=key)
    axs[0].axvline(2, color=C_NEW, ls='--', lw=1); axs[0].axvline(6, color=C_OLD, ls='--', lw=1)
    axs[0].set(xlabel='distance from the track to the nearest dead/shorted/unplugged pad tile [mm]',
               ylabel='miss rate [%]', title='Dead-pad border (live cells, >= 3 mm from the outer edge)')
    axs[1].axvline(0, color=INK2, lw=1); axs[1].axvline(2, color=C_NEW, ls='--', lw=1)
    axs[1].set(xlabel='signed distance to the OUTER edge of the pad array [mm] (+ inside)',
               ylabel='miss rate [%]', title='Outer edge (> 10 mm from dead pads)', yscale='log')
    for ax in axs[:2]:
        ax.grid(color=GRID); ax.legend(frameon=False, fontsize=8)
    # waterfall-ish: eff and N per step
    ax = axs[2]
    t2 = tab[tab.step != 'A2']
    ax.errorbar(range(len(t2)), t2.eff, t2.err, fmt='-o', color=C_NEW, lw=2, label='efficiency within R')
    ax.plot(range(len(t2)), t2.eff_corr, 's--', color=C_NEW, alpha=.6, label='accidental-corrected')
    c2 = tab[tab.step == 'A2']
    ax.errorbar([len(t2) - 0.15], c2.eff, c2.err, fmt='D', color=C_CTL, label='A2: guards 2/2, OLD alignment')
    for i, r in enumerate(t2.itertuples()):
        ax.annotate(f'N={r.n:,}', (i, r.eff), textcoords='offset points', xytext=(0, 9), ha='center',
                    fontsize=8, color=INK2)
    ax.set_xticks(range(len(t2)), t2.step)
    ax.set(ylabel='efficiency [%]', title='Step by step (see steps.csv)')
    ax.grid(color=GRID); ax.legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(f'{out}/alignment_validation.png', dpi=140); plt.close(fig)

    # residual maps A vs E
    fig, axs = plt.subplots(1, 2, figsize=(12, 5.2))
    for ax, key in zip(axs, ['A', 'E']):
        xy, T, act, hit = keep[key]
        cxy = cm.apply_T(T, cpad)
        d, j = cm.nearest(ev, xy, cxy, clu['eventId'].to_numpy())
        k = act & (d < 15)
        r = cxy[j[k]] - xy[k]
        df = pd.DataFrame({'a': np.floor(xy[k, 0] / 40), 'b': np.floor(xy[k, 1] / 40), 'rx': r[:, 0], 'ry': r[:, 1]})
        gb = df.groupby(['a', 'b']).agg(rx=('rx', 'median'), ry=('ry', 'median'), n=('rx', 'size'))
        gb = gb[gb.n > 100].reset_index()
        q = ax.quiver((gb.a + .5) * 40, (gb.b + .5) * 40, gb.rx, gb.ry, np.hypot(gb.rx, gb.ry),
                      cmap='viridis', clim=(0, 2.5), angles='xy', scale_units='xy', scale=0.06)
        plt.colorbar(q, ax=ax, label='|median cluster - track| [mm]')
        r_ = tab[tab.step == key].iloc[0]
        ax.set(title=f'{key}: scale {r_.scale:.3f}, residual sigma {r_.sig_x:.1f}/{r_.sig_y:.1f} mm (arrows x16.7)',
               xlabel='x M3 [mm]', ylabel='y M3 [mm]', aspect='equal')
    fig.tight_layout(); fig.savefig(f'{out}/residual_maps_A_vs_E.png', dpi=140); plt.close(fig)
    print(f'-> {out}')


if __name__ == '__main__':
    main()
