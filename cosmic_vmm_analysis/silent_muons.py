#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The "chamber silent" muons: scintillator-confirmed muons for which the probe chamber recorded no live hit at all
within the time window (n_any == 0 in a --dump-diag dump). They are ~10 % of the muons (24-26 % on chips 2 and 14,
FINDINGS_diag_20261118.md). Is it acceptance (the track crosses a region that is not read out) or read-out (data lost)?

Every fraction is on ON-time minus r x OFF-time tags (muons), scintillator window --scint, as channel_efficiency measure.

  map       per chamber: muons and silent muons in 10 mm bins of the predicted position (px, py), the silent fraction,
            with the pad centres coloured by chip                                   -> silent_map_<det>.png, silent_map.csv
  edges     silent fraction vs distance of the prediction to the chamber outline (convex hull of its pads) and to the
            nearest pad of another chip (chip boundary)                              -> silent_vs_edge.csv/.png
  chips     silent fraction per predicted chip (and hybrid), per configuration       -> silent_per_chip.csv
  time      per run and per chunk; clustering: P(silent | the previous muon in the same chamber was silent) vs P(silent |
            not silent), by the time between the two; silent fraction vs the time since the previous hit on the
            predicted chip (prev_chip) and vs the chamber hit rate of that second    -> silent_time*.csv/.png
Clustering and the rate dependence use ON-time tags only (OFF-time tags are fakes, mostly silent; their share is
printed).

Usage: silent_muons.py TAGS_DIR PSEUDO_DIR VETO_DIR OUT_DIR [--exclude-runs 233] [--scint 140:310] [--configs REGEX]
(lxplus, LCG view; reads the dumps only.)
"""
import argparse
import os
import re
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import probe_diagnostics as PD  # noqa: E402

DETS = ['P2_IN', 'P2_MID', 'P2_OUT']
DETV = {'P2_IN': range(2, 8), 'P2_MID': range(8, 14), 'P2_OUT': range(14, 20)}


def frac(s, w):
    W = np.sum(w)
    return (float(np.sum(w * s) / W) if W > 0 else np.nan), float(W)


def seg_dist(P, A, B):
    """distance of points P (n,2) to segments A->B (m,2): (n,m)."""
    AB = B - A
    t = np.clip(((P[:, None, :] - A[None]) * AB[None]).sum(-1) / (AB ** 2).sum(-1)[None], 0, 1)
    Q = A[None] + t[..., None] * AB[None]
    return np.hypot(*(P[:, None, :] - Q).transpose(2, 0, 1))


def geometry():
    from channel_efficiency import pad_geometry
    from scipy.spatial import ConvexHull, cKDTree, Delaunay
    X, Y, _ = pad_geometry()
    geo = {}
    for det, chips in DETV.items():
        pts, chip = [], []
        for c in chips:
            ok = np.isfinite(X[c])
            pts.append(np.c_[X[c][ok], Y[c][ok]])
            chip += [c] * int(ok.sum())
        pts, chip = np.vstack(pts), np.array(chip)
        h = ConvexHull(pts)
        hv = pts[h.vertices]
        geo[det] = dict(pts=pts, chip=chip, hull=hv, tri=Delaunay(hv),
                        trees={c: cKDTree(pts[chip == c]) for c in chips})
    return geo


def edge_distances(t, geo):
    """d_edge: distance to the chamber hull (negative outside); d_chip: distance to the nearest pad of another chip
    than the predicted pad's chip."""
    d_edge = np.full(len(t), np.nan)
    d_chip = np.full(len(t), np.nan)
    for det, g in geo.items():
        m = (t.det == det).to_numpy()
        if not m.any():
            continue
        P = np.c_[t.px.to_numpy()[m], t.py.to_numpy()[m]]
        A, B = g['hull'], np.roll(g['hull'], -1, axis=0)
        d = seg_dist(P, A, B).min(1)
        inside = g['tri'].find_simplex(P) >= 0
        d_edge[m] = np.where(inside, d, -d)
        own = t.vmm.to_numpy()[m]
        best = np.full(m.sum(), np.inf)
        for c, tr in g['trees'].items():
            dd = tr.query(P)[0]
            best = np.where(own != c, np.minimum(best, dd), best)
        d_chip[m] = best
    return d_edge, d_chip


def plot_maps(t, geo, out, bin_mm=10):
    import matplotlib.pyplot as plt
    rows = []
    for det in DETS:
        s = t[t.det == det]
        if s.empty:
            continue
        g = geo[det]
        x0, y0 = g['pts'].min(0) - 30
        x1, y1 = g['pts'].max(0) + 30
        bx, by = np.arange(x0, x1 + bin_mm, bin_mm), np.arange(y0, y1 + bin_mm, bin_mm)
        H = np.histogram2d(s.px, s.py, [bx, by], weights=s.w)[0]
        S = np.histogram2d(s.px, s.py, [bx, by], weights=s.w * s.silent)[0]
        with np.errstate(invalid='ignore', divide='ignore'):
            F = np.where(H >= 15, S / H, np.nan)
        for i in range(len(bx) - 1):
            for j in range(len(by) - 1):
                if H[i, j] > 0:
                    rows.append(dict(det=det, x=bx[i] + bin_mm / 2, y=by[j] + bin_mm / 2, muons=H[i, j], silent=S[i, j]))
        fig, ax = plt.subplots(1, 3, figsize=(19, 6.2))
        for a_, M, ttl, kw in ((ax[0], H, 'muons (on - r off)', dict(cmap='viridis')),
                               (ax[1], S, 'silent muons', dict(cmap='magma')),
                               (ax[2], F, 'silent fraction (bins with >= 15 muons)', dict(cmap='coolwarm', vmin=0, vmax=0.5))):
            im = a_.pcolormesh(bx, by, M.T, shading='flat', **kw)
            fig.colorbar(im, ax=a_, fraction=0.046)
            for c in DETV[det]:
                q = g['pts'][g['chip'] == c]
                a_.scatter(q[:, 0], q[:, 1], s=1.5, c='w' if a_ is not ax[2] else 'k', alpha=0.35)
                a_.text(*q.mean(0), str(c), ha='center', va='center', fontsize=11, weight='bold',
                        color='w' if a_ is not ax[2] else 'k')
            hv = np.vstack([g['hull'], g['hull'][:1]])
            a_.plot(hv[:, 0], hv[:, 1], 'k-' if a_ is ax[2] else 'w-', lw=1)
            a_.set_title(f'{det}: {ttl}')
            a_.set_xlabel('x (mm)'); a_.set_ylabel('y (mm)'); a_.set_aspect('equal')
        fig.tight_layout(); fig.savefig(f'{out}/silent_map_{det}.png', dpi=110); plt.close(fig)
    pd.DataFrame(rows).to_csv(f'{out}/silent_map.csv', index=False)


def binned(t, col, edges, by=('det',)):
    rows = []
    lab = [f'{lo:g}..{hi:g}' for lo, hi in zip(edges[:-1], edges[1:])]
    t = t.assign(_b=pd.cut(t[col], edges, labels=lab, right=False))
    for k, g in t.groupby(list(by) + ['_b'], observed=True):
        f, W = frac(g.silent.to_numpy(), g.w.to_numpy())
        rows.append(dict(zip(list(by) + [col], k if isinstance(k, tuple) else (k,)), silent_frac=round(f, 4),
                         muons=round(W, 1)))
    return pd.DataFrame(rows)


def clustering(t):
    """ON tags, per (run, det) in wall order: P(silent | previous silent) vs P(silent | previous not silent),
    by the gap to the previous ON tag."""
    rows = []
    gaps = (0, 1, 10, 60, 600, np.inf)
    o = t[t.on].sort_values(['run', 'det', 'wall'])
    prev_s = o.groupby(['run', 'det']).silent.shift(1)
    gap = o.wall - o.groupby(['run', 'det']).wall.shift(1)
    o = o.assign(prev_s=prev_s, gap=gap).dropna(subset=['prev_s'])
    for det, g in o.groupby('det'):
        for lo, hi in zip(gaps[:-1], gaps[1:]):
            q = g[(g.gap >= lo) & (g.gap < hi)]
            a = q[q.prev_s == 1]
            b = q[q.prev_s == 0]
            rows.append(dict(det=det, gap_s=f'{lo:g}..{hi:g}', n_prev_silent=len(a),
                             p_silent_after_silent=round(a.silent.mean(), 4) if len(a) else np.nan,
                             n_prev_ok=len(b), p_silent_after_ok=round(b.silent.mean(), 4) if len(b) else np.nan))
    return pd.DataFrame(rows)


def dispersion(t, bin_s=300):
    """ON tags: silent counts in bin_s windows per (run, det) against binomial: index = var / (n p (1-p))."""
    rows = []
    o = t[t.on]
    for (run, det), g in o.groupby(['run', 'det']):
        b = ((g.wall - g.wall.min()) // bin_s).astype(int)
        c = g.groupby(b).silent.agg(['sum', 'count'])
        c = c[c['count'] >= 5]
        if len(c) < 4:
            continue
        p = c['sum'].sum() / c['count'].sum()
        exp_var = np.mean(c['count'] * p * (1 - p))
        rows.append(dict(run=run, det=det, bins=len(c), p_silent_on=round(p, 4),
                         dispersion=round(float(np.var(c['sum'] - c['count'] * p) / exp_var), 3) if exp_var > 0 else np.nan,
                         max_bin_frac=round(float((c['sum'] / c['count']).max()), 3)))
    return pd.DataFrame(rows)


def chamber_rate(t, RATES):
    """hits/s of the probe chamber in the tag's second (from the dump's per-second rates)."""
    r = np.full(len(t), np.nan)
    key = {'P2_IN': 3, 'P2_MID': 2, 'P2_OUT': 1}
    byp = {R['pid']: R for R in RATES}
    for pid, idx in t.groupby('pid').indices.items():
        R = byp.get(pid)
        if R is None:
            continue
        s = t.iloc[idx]
        for det, j in s.groupby('det').indices.items():
            arr = np.asarray(R[f'hits_{key[det]}'])
            k = np.floor(s.wall.to_numpy()[j] - R['t0']).astype(int)
            ok = (k >= 0) & (k < len(arr))
            v = np.full(len(k), np.nan)
            v[ok] = arr[k[ok]]
            r[idx[j]] = v
    return r


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('tags'); ap.add_argument('pseudo'); ap.add_argument('veto'); ap.add_argument('out')
    ap.add_argument('--exclude-runs', default='233')
    ap.add_argument('--scint', default='140:310')
    ap.add_argument('--configs', default='', help='regex on the configuration label: keep only these')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    log = open(f'{a.out}/silent_summary.txt', 'w')

    def say(*x):
        s = ' '.join(str(v) for v in x)
        print(s, flush=True)
        log.write(s + '\n')

    tags, near, RATES, vetoes, win = PD.load_all(a)
    if a.configs:
        tags = tags[tags.cfg.str.contains(a.configs, regex=True)]
    tags = tags.copy()
    tags['silent'] = (tags.n_any == 0).astype(float)
    # px/py are not kept by load_all: read them back per chunk (same selection order)
    tags = add_positions(tags, a)
    geo = geometry()
    tags['d_edge'], tags['d_chip'] = edge_distances(tags, geo)
    tags['hyb'] = tags.vmm // 2
    tags['rate'] = chamber_rate(tags, RATES)
    say(f'{len(tags)} tags ({tags.on.sum()} on), runs {tags.run.min()}-{tags.run.max()}, win {win} ns')
    for det in DETS:
        s = tags[tags.det == det]
        f, W = frac(s.silent.to_numpy(), s.w.to_numpy())
        off_share = float(-s.w[~s.on].sum() / s.w[s.on].sum()) if s.on.any() else np.nan
        say(f'{det}: muons {W:.0f}, silent fraction {f:.4f}; off-time share of the on tags {off_share:.3f}')
    plot_maps(tags, geo, a.out)

    e1 = binned(tags, 'd_edge', [-50, 0, 5, 10, 15, 20, 30, 40, 60, 100, 1000])
    e2 = binned(tags, 'd_chip', [0, 2, 4, 6, 8, 10, 15, 20, 30, 1000])
    e1.to_csv(f'{a.out}/silent_vs_edge.csv', index=False)
    e2.to_csv(f'{a.out}/silent_vs_chipboundary.csv', index=False)
    say('\n== silent fraction vs distance to the chamber outline (mm; negative = outside the pad hull)')
    say(e1.to_string(index=False))
    say('\n== silent fraction vs distance to the nearest pad of another chip (mm)')
    say(e2.to_string(index=False))

    rows = []
    for (det, cfg, v), g in tags.groupby(['det', 'cfg', 'vmm']):
        f, W = frac(g.silent.to_numpy(), g.w.to_numpy())
        rows.append(dict(det=det, cfg=cfg, chip=v, hyb=v // 2, silent_frac=round(f, 4), muons=round(W, 1)))
    pc = pd.DataFrame(rows)
    pc.to_csv(f'{a.out}/silent_per_chip_config.csv', index=False)
    rows = []
    for (det, v), g in tags.groupby(['det', 'vmm']):
        f, W = frac(g.silent.to_numpy(), g.w.to_numpy())
        fi, Wi = frac(g.silent.to_numpy()[g.d_edge.to_numpy() > 20], g.w.to_numpy()[g.d_edge.to_numpy() > 20])
        fc, Wc = frac(g.silent.to_numpy()[g.d_chip.to_numpy() > 10], g.w.to_numpy()[g.d_chip.to_numpy() > 10])
        rows.append(dict(det=det, chip=v, hyb=v // 2, silent_frac=round(f, 4), muons=round(W, 1),
                         silent_frac_edge_gt20=round(fi, 4), muons_edge_gt20=round(Wi, 1),
                         silent_frac_chipb_gt10=round(fc, 4), muons_chipb_gt10=round(Wc, 1)))
    pchip = pd.DataFrame(rows)
    pchip.to_csv(f'{a.out}/silent_per_chip.csv', index=False)
    say('\n== per predicted chip (all configs); also > 20 mm inside the outline, > 10 mm from another chip')
    say(pchip.to_string(index=False))

    rows = []
    for (det, run, cfg), g in tags.groupby(['det', 'run', 'cfg']):
        f, W = frac(g.silent.to_numpy(), g.w.to_numpy())
        rows.append(dict(det=det, run=run, cfg=cfg, silent_frac=round(f, 4), muons=round(W, 1),
                         t_start=float(g.wall.min())))
    pr = pd.DataFrame(rows)
    pr.to_csv(f'{a.out}/silent_per_run.csv', index=False)
    rows = []
    for (det, pid), g in tags.groupby(['det', 'pid']):
        f, W = frac(g.silent.to_numpy(), g.w.to_numpy())
        rows.append(dict(det=det, pid=pid, run=g.run.iloc[0], silent_frac=round(f, 4), muons=round(W, 1),
                         t_mid=float(g.wall.median())))
    pk = pd.DataFrame(rows)
    pk.to_csv(f'{a.out}/silent_per_chunk.csv', index=False)

    cl = clustering(tags)
    cl.to_csv(f'{a.out}/silent_clustering.csv', index=False)
    say('\n== clustering (ON tags): P(silent) after a silent / a non-silent muon of the same chamber, by the gap')
    say(cl.to_string(index=False))
    dp = dispersion(tags)
    dp.to_csv(f'{a.out}/silent_dispersion.csv', index=False)
    say('\n== dispersion of silent counts in 300 s bins per run (1 = binomial, >1 = bursts); per chamber median / max')
    say(dp.groupby('det').dispersion.describe()[['count', '50%', 'max']].to_string())

    pcb = binned(tags[tags.on], 'prev_chip', [0, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, np.inf])
    pcb = pcb.rename(columns={'silent_frac': 'silent_frac_on', 'muons': 'on_tags'})
    pcb.to_csv(f'{a.out}/silent_vs_prev_chip.csv', index=False)
    say('\n== ON tags: silent fraction vs time since the previous hit on the predicted chip (ns)')
    say(pcb.to_string(index=False))
    rq = tags[tags.on & tags.rate.notna()].copy()
    rq['rate_rel'] = rq.rate / rq.groupby(['run', 'det']).rate.transform('median')
    rb = binned(rq, 'rate_rel', [0, 0.5, 0.8, 0.9, 1.0, 1.1, 1.25, 1.5, 2, 5, np.inf])
    rb = rb.rename(columns={'silent_frac': 'silent_frac_on', 'muons': 'on_tags'})
    rb.to_csv(f'{a.out}/silent_vs_rate.csv', index=False)
    say('\n== ON tags: silent fraction vs the probe chamber hit rate in that second, relative to the median of the '
        'run and chamber')
    say(rb.to_string(index=False))

    plot_time(pk, pr, e1, e2, a.out)
    tags[['run', 'pid', 'cfg', 'det', 'on', 'w', 'vmm', 'ch', 'px', 'py', 'wall', 'silent', 'd_edge', 'd_chip',
          'prev_chip', 'rate', 'n100us', 'sat100us']].to_parquet(f'{a.out}/silent_tags.parquet') \
        if _has_parquet() else None
    say(f'\noutputs in {a.out}')


def _has_parquet():
    try:
        import pyarrow  # noqa: F401
        return True
    except ImportError:
        return False


def add_positions(tags, a):
    """load_all keeps vmm/ch but not the predicted position: read px, py back with the same selection."""
    import glob
    lo, hi = map(float, a.scint.split(':'))
    pos = []
    key = {'P2_IN': 3, 'P2_MID': 2, 'P2_OUT': 1}
    for pid in tags.pid.unique():
        fz = f'{a.tags}/tags_run{pid}.npz'
        Z = np.load(fz)
        for det, p in key.items():
            dts = Z[f'p{p}_dt_scint']
            idx = np.nonzero((dts >= lo) & (dts <= hi))[0]
            if not len(idx):
                continue
            pos.append(pd.DataFrame(dict(tid=idx + p * 10_000_000 + pid * 100_000_000,
                                         px=Z[f'p{p}_px'][idx], py=Z[f'p{p}_py'][idx])))
    return tags.merge(pd.concat(pos, ignore_index=True), on='tid', how='left')


def plot_time(pk, pr, e1, e2, out):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 2, figsize=(16, 9))
    for det in DETS:
        s = pk[(pk.det == det) & (pk.muons > 5)]
        ax[0, 0].plot(pd.to_datetime(s.t_mid, unit='s'), s.silent_frac, '.', ms=3, label=det)
        s = pr[(pr.det == det)]
        ax[0, 1].errorbar(s.run, s.silent_frac, yerr=np.sqrt(np.clip(s.silent_frac * (1 - s.silent_frac), 0, None)
                                                            / np.maximum(s.muons, 1)), fmt='o', ms=3, label=det)
        s = e1[e1.det == det]
        ax[1, 0].plot(range(len(s)), s.silent_frac, 'o-', label=det)
        ax[1, 0].set_xticks(range(len(s))); ax[1, 0].set_xticklabels(s.d_edge, rotation=45)
        s = e2[e2.det == det]
        ax[1, 1].plot(range(len(s)), s.silent_frac, 'o-', label=det)
        ax[1, 1].set_xticks(range(len(s))); ax[1, 1].set_xticklabels(s.d_chip, rotation=45)
    ax[0, 0].set_title('silent fraction per chunk (~7 min)'); ax[0, 0].set_ylim(-0.05, 0.6)
    ax[0, 1].set_title('silent fraction per run'); ax[0, 1].set_xlabel('run')
    ax[1, 0].set_title('vs distance to the chamber outline (mm)')
    ax[1, 1].set_title('vs distance to the nearest pad of another chip (mm)')
    for x in ax.flat:
        x.legend(); x.grid(alpha=0.3); x.set_ylabel('silent fraction')
    fig.tight_layout(); fig.savefig(f'{out}/silent_time_edges.png', dpi=110); plt.close(fig)


if __name__ == '__main__':
    main()
