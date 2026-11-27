#!/usr/bin/env python3
"""40_nohit_angle_vs_field.py -- does the near-vertical no-hit excess depend on
the drift field (collection / transparency / attachment) or on the mesh voltage
(gain)? A purely geometric cause would not move with either.

Uses the 11-24 det2 series (same detector, before it was repositioned):
alignment run (430/730), drift scan (mesh 430, drift 430..955) and mesh scan
(gap 300 V, mesh 340..420). One alignment for the whole series: trimmed fit +
slope-fitted plane z on the alignment run (the detector did not move during it).
Active area as stage 06 with the det2_long_repos rules (ACTIVE_R + dead-pad
cells + 2 mm dead and outer-edge guards). No-hit = no mapped P2 hit in the event.
Spark-vetoed events are dropped per scan point (stage-16 veto).

Out: 40_nohit_vs_field/{nohit_vs_field.csv,.png}, summary.txt
Usage: python3 40_nohit_angle_vs_field.py
"""
import os
import sys
import glob

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

TH = [(0, 6), (6, 12), (12, 90)]
COL = {(0, 6): '#ec4647', (6, 12): '#f2a600', (12, 90): '#257de0'}
GRID = '#e6e6e6'


def point(cfg, ct, subdir, T, z, drop=None):
    hits_dir = os.path.join(subdir, 'combined_hits_root')
    m3_dir = os.path.join(subdir, 'm3_tracking_root')
    drop = p2io.drop_pads_for(cfg, ct, hits_dir=hits_dir)
    ep = pa.load_m3_endpoints(m3_dir, qa.M3_CHI2_CUT, qa.M3_MIN_NCLUS)
    _, hit_ev = pa.load_p2_centroids(hits_dir, ct, min_amp=cfg.MIN_AMP, drop_pads=drop)
    hv = os.path.join(subdir, 'hv_monitor.csv')
    if os.path.isfile(hv):
        sv = ps.SparkVeto.from_csv(hv, cfg)
        bad = sv.vetoed_ids_from_hits(hits_dir, ct.attrs['feus'], min_amp=cfg.MIN_AMP)
        ep = ep[~ep.eventId.isin(bad)]
    x, y = pa.project_to_z(ep, z); xy = np.column_stack([x, y])
    padc = ct[ct['mapped']].drop_duplicates('channel_id')
    pxy = np.column_stack(T.apply(padc.pad_cx.to_numpy(), padc.pad_cy.to_numpy()))
    nn, _ = cKDTree(pxy).query(xy)
    act = pa.mask_dead_cells(nn, xy, ct, T) <= cfg.ACTIVE_R
    tx, ty, th = cm.ray_dirs(ep)
    nohit = ~ep.eventId.isin(hit_ev).to_numpy()
    return act, nohit, th


def main():
    cfg = qa.get_config('det2_initial2')
    out = os.path.join(qa.get_config('det2_long_repos').OUT_BASE, '..', '40_nohit_vs_field_11-24')
    out = os.path.abspath(out); os.makedirs(out, exist_ok=True)
    ct = sc.channel_table(cfg)
    ct.attrs['dead_guard'] = 2.0
    ct.attrs['edge_guard'] = 2.0
    L = []
    # ---- alignment on the alignment run: trimmed fit, z from residual-vs-slope ----
    ali = cfg.sub_run_dir
    ep = pa.load_m3_endpoints(os.path.join(ali, 'm3_tracking_root'), qa.M3_CHI2_CUT, qa.M3_MIN_NCLUS)
    drop = p2io.drop_pads_for(cfg, ct, hits_dir=os.path.join(ali, 'combined_hits_root'))
    cen, _ = pa.load_p2_centroids(os.path.join(ali, 'combined_hits_root'), ct, drop_pads=drop)
    tx, ty, _ = cm.ray_dirs(ep)
    z = 240.0
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
    L.append(f'11-24 alignment run: z = {z:.1f} mm, scale {T.s:.4f}, rot {T.rotation_deg:.2f}, '
             f'rmse {T.rmse:.2f} mm (n {len(m)})')

    rows = []
    subs = [('alignment', 430, 730, ali)]
    for d in sorted(glob.glob(os.path.join(cfg.run_dir, '*_scan_det2_*'))):
        nm = os.path.basename(d); a, b = nm.split('_')[-2:]
        subs.append(('drift' if nm.startswith('drift') else 'mesh', int(a), int(b), d))
    for kind, mesh, drift, d in subs:
        if not glob.glob(os.path.join(d, 'combined_hits_root', '*.root')):
            continue
        act, nohit, th = point(cfg, ct, d, T, z)
        for lo, hi in TH + [(0, 90)]:
            k = act & (th >= lo) & (th < hi)
            rows.append(dict(scan=kind, mesh=mesh, drift=drift, field=(drift - mesh) / 0.4,
                             th_lo=lo, th_hi=hi, n=int(k.sum()), k=int(nohit[k].sum())))
        print(f'{kind:9s} {mesh}/{drift}: active {act.sum():5d}  no-hit {100*nohit[act].mean():5.2f} %  '
              f'<6deg {100*nohit[act & (th < 6)].mean():5.2f} %  >=12deg {100*nohit[act & (th >= 12)].mean():5.2f} %')
    D = pd.DataFrame(rows)
    D['p'] = 100 * D.k / D.n.clip(lower=1)
    D['err'] = 100 * np.sqrt(np.clip(D.k, 1, None)) / D.n.clip(lower=1)
    D.to_csv(f'{out}/nohit_vs_field.csv', index=False)

    # ---- plots: no-hit vs drift field (drift scan) and vs mesh V (mesh scan), per angle class ----
    fig, axs = plt.subplots(1, 3, figsize=(18, 4.8))
    for ax, sc_, xcol, xl in [(axs[0], 'drift', 'field', 'drift field [V/cm] (mesh 430 V)'),
                              (axs[1], 'mesh', 'mesh', 'mesh voltage [V] (drift gap 300 V)')]:
        s = D[D.scan == sc_]
        for (lo, hi) in TH:
            q = s[(s.th_lo == lo) & (s.th_hi == hi)].sort_values(xcol)
            ax.errorbar(q[xcol], q.p, q.err, fmt='o-', color=COL[(lo, hi)], lw=1.8, ms=4, label=f'{lo}-{hi} deg')
        ax.set(xlabel=xl, ylabel='no P2 hit at all [%]', yscale='log', title=f'{sc_} scan (11-24)')
        ax.grid(color=GRID); ax.legend(frameon=False, fontsize=8)
    # excess ratio vertical / inclined vs drift field (grouped in 4 bins for stats)
    s = D[(D.scan == 'drift')]
    v = s[s.th_hi == 6].set_index('drift'); i = s[s.th_lo == 12].set_index('drift')
    ax = axs[2]
    ax.errorbar(v.field, v.p - i.p, np.hypot(v.err, i.err), fmt='o-', color='#5c5c5c', lw=1.8,
                label='no-hit(<6 deg) - no-hit(>=12 deg)')
    ax.axhline(0, color='#9c9c9c', lw=1)
    ax.set(xlabel='drift field [V/cm]', ylabel='vertical excess [points]', title='Vertical-track excess vs drift field')
    ax.grid(color=GRID); ax.legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(f'{out}/nohit_vs_field.png', dpi=140); plt.close(fig)

    # pooled numbers: drift field groups
    D3 = D[~((D.th_lo == 0) & (D.th_hi == 90))]          # angle classes only (no all-angle rows)
    s = D3[D3.scan == 'drift'].copy()
    s['grp'] = pd.cut(s.field, [0, 200, 400, 600, 800, 1000, 1200, 1400])
    g_ = s.groupby(['grp', 'th_lo'], observed=True)[['n', 'k']].sum().reset_index()
    g_['p'] = 100 * g_.k / g_.n
    L.append('\ndrift scan, pooled field groups [V/cm]: no-hit % (n) per angle class')
    L.append(g_.pivot(index='grp', columns='th_lo', values='p').round(2).to_string())
    L.append(g_.pivot(index='grp', columns='th_lo', values='n').to_string())
    s = D3[D3.scan == 'mesh'].copy()
    s['grp'] = pd.cut(s.mesh, [335, 355, 375, 395, 415, 425])
    g_ = s.groupby(['grp', 'th_lo'], observed=True)[['n', 'k']].sum().reset_index()
    g_['p'] = 100 * g_.k / g_.n
    L.append('\nmesh scan, pooled mesh-voltage groups: no-hit % per angle class')
    L.append(g_.pivot(index='grp', columns='th_lo', values='p').round(2).to_string())
    a = D[D.scan == 'alignment']
    L.append('\nalignment run 430/730: ' + ', '.join(f'{r.th_lo}-{r.th_hi}: {r.p:.2f} % (n {r.n})' for r in a.itertuples()))
    with open(f'{out}/summary.txt', 'w') as f:
        f.write('\n'.join(L) + '\n')
    print('\n'.join(L))
    print('->', out)


if __name__ == '__main__':
    main()
