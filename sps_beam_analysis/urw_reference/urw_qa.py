#!/usr/bin/env python3
"""
Stand-alone uRWELL QA for TB_July2026_H4 -- route 1 of handoff §6.

Per (run, sub_run, detector) it produces, from the strip maps + hits_root:
  01_occupancy   hits and charge vs strip position, per view, dead-strip list
  02_spectra     amplitude / significance / cluster-charge spectra, Landau MPV
  03_cluster     cluster size and multiplicity, broken out per pitch zone
  04_beamspot    2D map of leading-cluster (x, y) + projections with core fits
  05_timing      time_of_max distribution and rate vs time in the sub_run
  summary.csv    one row per detector/view/zone with the numbers above

Usage
  PY=/local/home/banco/DAQ_Control_Dream_Beam/.venv/bin/python
  unset PYTHONPATH
  $PY urw_qa.py --sub-run nominal_00
  $PY urw_qa.py --sub-run all --max-hits 4000000     # whole drift/mesh scan
"""

import os
import re
import sys
import json
import argparse

import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from urw_lib import (UrwGeometry, feu_hit_files, iter_hits, cluster_hits,  # noqa: E402
                     leading_points, URW_DETS, DET_DIR_DEFAULT)

RUN_DIR_DEFAULT = '/local/home/banco/P2_data/TB_July2026_H4/runs/drift_mesh_scan_1'
# products go under the analysis root like every other stage, not into the repo
OUT_DEFAULT = os.path.join(
    os.environ.get('SPS_ANALYSIS_ROOT',
                   '/local/home/banco/P2_data/TB_July2026_H4/analysis'),
    'urw_referenced_efficiency', 'qa')


# ------------------------------------------------------------------ fits ----

def _gauss(x, a, mu, sig):
    return a * np.exp(-0.5 * ((x - mu) / sig) ** 2)


def fit_core(vals, nbins=100, n_iter=3):
    """Gaussian core fit (refit inside +-2 sigma) -- same recipe as 23_beam_profile."""
    vals = np.asarray(vals, float)
    vals = vals[np.isfinite(vals)]
    out = dict(mu=np.nan, sigma=np.nan, ok=False)
    if len(vals) < 50:
        return out
    counts, edges = np.histogram(vals, bins=nbins)
    ctr = 0.5 * (edges[:-1] + edges[1:])
    mu, sig = float(np.median(vals)), float(vals.std())
    for _ in range(n_iter):
        win = np.abs(ctr - mu) < 2 * sig
        if win.sum() < 5:
            break
        try:
            p, _ = curve_fit(_gauss, ctr[win], counts[win],
                             p0=[counts[win].max(), mu, sig], maxfev=8000)
        except (RuntimeError, ValueError):
            return out
        mu, sig = float(p[1]), abs(float(p[2]))
    out.update(mu=mu, sigma=sig, ok=True)
    return out


def landau_mpv(vals, nbins=120, qlo=1, qhi=99):
    """Most-probable value from the mode of the charge histogram (robust, fit-free)."""
    vals = np.asarray(vals, float)
    vals = vals[np.isfinite(vals) & (vals > 0)]
    if len(vals) < 50:
        return np.nan
    lo, hi = np.percentile(vals, [qlo, qhi])
    counts, edges = np.histogram(vals, bins=nbins, range=(lo, hi))
    # smooth over 5 bins so the mode is not a statistical spike
    k = np.ones(5) / 5
    sm = np.convolve(counts, k, mode='same')
    return float(0.5 * (edges[np.argmax(sm)] + edges[np.argmax(sm) + 1]))


# ------------------------------------------------------------ accumulation --

class DetAccumulator:
    def __init__(self, geo):
        self.geo = geo
        self.n_hits_raw = 0
        self.n_hits_sel = 0
        self.hit_counts = np.zeros(512, np.int64)
        self.hit_charge = np.zeros(512, np.float64)
        self.amp_hist = np.zeros(200, np.int64)
        self.amp_edges = np.linspace(0, 2000, 201)
        self.sig_hist = np.zeros(200, np.int64)
        self.sig_edges = np.linspace(0, 100, 201)
        self.clusters = []
        self.points = []
        self.event_ids = set()

    def add(self, chunk, args):
        ch = chunk['channel'].astype(np.int32)
        mine = self.geo.mapped[ch]
        self.n_hits_raw += int(mine.sum())
        if mine.any():
            amp = chunk['amplitude'][mine].astype(np.float64)
            sig = chunk['significance'][mine].astype(np.float64)
            self.amp_hist += np.histogram(amp, bins=self.amp_edges)[0]
            self.sig_hist += np.histogram(sig, bins=self.sig_edges)[0]

        sel = mine.copy()
        if args.min_amp > 0:
            sel &= chunk['amplitude'] >= args.min_amp
        if args.min_signif > 0:
            sel &= chunk['significance'] >= args.min_signif
        self.n_hits_sel += int(sel.sum())
        if sel.any():
            np.add.at(self.hit_counts, ch[sel], 1)
            np.add.at(self.hit_charge, ch[sel], chunk['amplitude'][sel].astype(np.float64))

        cl = cluster_hits(chunk, self.geo, min_amp=args.min_amp,
                          min_signif=args.min_signif, gap_factor=args.gap_factor)
        if len(cl):
            self.clusters.append(cl)
            self.points.append(leading_points(cl))

    def finish(self):
        self.clusters = (pd.concat(self.clusters, ignore_index=True)
                         if self.clusters else cluster_hits({}, self.geo))
        self.points = (pd.concat(self.points, ignore_index=True)
                       if len(self.points) else pd.DataFrame())


# --------------------------------------------------------------- plotting ---

def plot_occupancy(acc, out_dir, tag):
    geo, counts = acc.geo, acc.hit_counts
    fig, axes = plt.subplots(2, 2, figsize=(13, 7))
    for j, view in enumerate(('x', 'y')):
        m = geo.mapped & (geo.view == view)
        idx = np.flatnonzero(m)
        order = np.argsort(geo.pos[idx])
        idx = idx[order]
        p, c = geo.pos[idx], counts[idx]
        colors = plt.cm.tab10(geo.zone[idx] % 10)

        ax = axes[0, j]
        ax.bar(p, c, width=geo.pitch[idx] * 0.9, color=colors)
        ax.set_xlabel(f'local {view} [mm]')
        ax.set_ylabel('hits')
        ax.set_title(f'{view} view -- hits vs strip position')

        ax = axes[1, j]
        mean_amp = np.divide(acc.hit_charge[idx], np.maximum(c, 1))
        ax.bar(p, mean_amp, width=geo.pitch[idx] * 0.9, color=colors)
        ax.set_xlabel(f'local {view} [mm]')
        ax.set_ylabel('mean amplitude [ADC]')
        ax.set_title(f'{view} view -- mean hit amplitude')
    fig.suptitle(f'{tag} -- occupancy (colour = pitch zone)')
    fig.tight_layout()
    fig.savefig(f'{out_dir}/01_occupancy.png', dpi=110)
    plt.close(fig)

    dead = np.flatnonzero(geo.mapped & (counts == 0))
    hot = np.flatnonzero(geo.mapped & (counts > 5 * np.median(counts[geo.mapped])))
    return dead, hot


def plot_spectra(acc, out_dir, tag):
    cl = acc.clusters
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    ctr = 0.5 * (acc.amp_edges[:-1] + acc.amp_edges[1:])
    axes[0].step(ctr, acc.amp_hist, where='mid')
    axes[0].set_yscale('log')
    axes[0].set_xlabel('hit amplitude [ADC]')
    axes[0].set_ylabel('hits')
    axes[0].set_title('raw hit amplitude')

    ctr = 0.5 * (acc.sig_edges[:-1] + acc.sig_edges[1:])
    axes[1].step(ctr, acc.sig_hist, where='mid')
    axes[1].set_yscale('log')
    axes[1].set_xlabel('significance')
    axes[1].set_title('hit significance')

    mpvs = {}
    for view in ('x', 'y'):
        q = cl.loc[cl['view'] == view, 'charge'].values
        if len(q) < 50:
            continue
        hi = np.percentile(q, 99.5)
        axes[2].hist(q, bins=120, range=(0, hi), histtype='step', label=f'{view} view')
        mpvs[view] = landau_mpv(q)
    axes[2].set_xlabel('cluster charge [sum ADC]')
    axes[2].set_title('cluster charge  MPV ' +
                      ', '.join(f'{v}={m:.0f}' for v, m in mpvs.items()))
    axes[2].legend()
    fig.suptitle(f'{tag} -- spectra')
    fig.tight_layout()
    fig.savefig(f'{out_dir}/02_spectra.png', dpi=110)
    plt.close(fig)
    return mpvs


def plot_clusters(acc, out_dir, tag):
    cl, geo = acc.clusters, acc.geo
    zt = geo.zone_table()
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    for _, z in zt.iterrows():
        s = cl.loc[(cl['zone'] == z['zone']) & (~cl['zone_mixed']), 'size'].values
        if len(s) < 50:
            continue
        axes[0].hist(s, bins=np.arange(0.5, 15.5), histtype='step', density=True,
                     label=f"{z['view']} p={z['pitch']}mm  <n>={s.mean():.2f}")
    axes[0].set_xlabel('cluster size [strips]')
    axes[0].set_ylabel('fraction')
    axes[0].set_title('cluster size per pitch zone (pure clusters)')
    axes[0].legend(fontsize=7)

    for _, z in zt.iterrows():
        q = cl.loc[(cl['zone'] == z['zone']) & (~cl['zone_mixed']), 'charge'].values
        if len(q) < 50:
            continue
        axes[1].hist(q, bins=100, range=(0, np.percentile(q, 99.5)), histtype='step',
                     density=True, label=f"{z['view']} p={z['pitch']}mm  MPV={landau_mpv(q):.0f}")
    axes[1].set_xlabel('cluster charge [sum ADC]')
    axes[1].set_title('cluster charge per pitch zone')
    axes[1].legend(fontsize=7)

    for view in ('x', 'y'):
        n = cl.loc[cl['view'] == view].groupby('eventId').size()
        if not len(n):
            continue
        axes[2].hist(n, bins=np.arange(0.5, 12.5), histtype='step', density=True,
                     label=f'{view}  <N>={n.mean():.2f}')
    axes[2].set_xlabel('clusters / event')
    axes[2].set_title('cluster multiplicity')
    axes[2].legend()
    fig.suptitle(f'{tag} -- clustering')
    fig.tight_layout()
    fig.savefig(f'{out_dir}/03_cluster.png', dpi=110)
    plt.close(fig)


def plot_beamspot(acc, out_dir, tag):
    pts, geo = acc.points, acc.geo
    fits = {}
    if not len(pts):
        return fits
    both = pts.dropna(subset=['x', 'y'])
    fig = plt.figure(figsize=(13, 5))
    ax = fig.add_subplot(1, 3, 1)
    if len(both):
        ax.hist2d(both['x'], both['y'], bins=[128, 128],
                  range=[[0, geo.active_size[0]], [0, geo.active_size[1]]], cmin=1)
    ax.set_xlabel('local x [mm]')
    ax.set_ylabel('local y [mm]')
    ax.set_title(f'leading-cluster xy  (n={len(both)})')
    ax.set_aspect('equal')

    for i, view in enumerate(('x', 'y')):
        ax = fig.add_subplot(1, 3, 2 + i)
        v = pts[view].dropna().values
        if len(v) < 50:
            continue
        size = geo.active_size[i]
        ax.hist(v, bins=128, range=(0, size), histtype='step')
        f = fit_core(v, nbins=128)
        fits[view] = f
        if f['ok']:
            xs = np.linspace(0, size, 400)
            h, e = np.histogram(v, bins=128, range=(0, size))
            amp = h.max()
            ax.plot(xs, amp * np.exp(-0.5 * ((xs - f['mu']) / f['sigma']) ** 2), 'r-',
                    label=f"core mu={f['mu']:.2f} sig={f['sigma']:.2f} mm")
            ax.legend(fontsize=8)
        ax.set_xlabel(f'local {view} [mm]')
        ax.set_title(f'{view} profile')
    fig.suptitle(f'{tag} -- beam profile')
    fig.tight_layout()
    fig.savefig(f'{out_dir}/04_beamspot.png', dpi=110)
    plt.close(fig)
    return fits


def plot_timing(acc, out_dir, tag):
    cl = acc.clusters
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for view in ('x', 'y'):
        t = cl.loc[cl['view'] == view, 'time'].values
        t = t[np.isfinite(t)]
        if len(t) < 50:
            continue
        axes[0].hist(t, bins=100, histtype='step',
                     label=f'{view}  median={np.median(t):.1f}')
    axes[0].set_xlabel('charge-weighted time_of_max [samples]')
    axes[0].set_title('cluster time')
    axes[0].legend()

    t = cl['t_ns'].values / 1e9
    if len(t):
        t = t - t.min()
        axes[1].hist(t, bins=200, histtype='step')
        axes[1].set_xlabel('time in sub-run [s]')
        axes[1].set_ylabel('clusters / bin')
        axes[1].set_title('rate vs time (spill structure)')
    fig.suptitle(f'{tag} -- timing')
    fig.tight_layout()
    fig.savefig(f'{out_dir}/05_timing.png', dpi=110)
    plt.close(fig)


# ------------------------------------------------------------------- main ---

def analyse(geo, sub_run_dir, args, out_dir, tag):
    files = feu_hit_files(sub_run_dir, feu=geo.feu_num)
    if not files:
        print(f'  !! no FEU-{geo.feu_num} hits files in {sub_run_dir}')
        return None
    acc = DetAccumulator(geo)
    # feu= matters when feu_hit_files falls back to the combined files
    for chunk in iter_hits(files, max_hits=args.max_hits, progress=args.verbose,
                           feu=geo.feu_num):
        acc.add(chunk, args)
    acc.finish()
    os.makedirs(out_dir, exist_ok=True)

    dead, hot = plot_occupancy(acc, out_dir, tag)
    mpvs = plot_spectra(acc, out_dir, tag)
    plot_clusters(acc, out_dir, tag)
    fits = plot_beamspot(acc, out_dir, tag)
    plot_timing(acc, out_dir, tag)

    cl = acc.clusters
    n_ev = cl['eventId'].nunique() if len(cl) else 0
    rows = []
    zt = geo.zone_table()
    for _, z in zt.iterrows():
        m = (cl['zone'] == z['zone']) & (~cl['zone_mixed'])
        sub = cl[m]
        rows.append(dict(
            detector=geo.name, sub_run=os.path.basename(sub_run_dir), view=z['view'],
            zone=int(z['zone']), pitch=z['pitch'], interpitch=z['interpitch'],
            n_strips=int(z['n_strips']), n_clusters=int(m.sum()),
            mean_size=float(sub['size'].mean()) if len(sub) else np.nan,
            mpv_charge=landau_mpv(sub['charge'].values) if len(sub) else np.nan,
            mean_charge=float(sub['charge'].mean()) if len(sub) else np.nan,
            median_time=float(sub['time'].median()) if len(sub) else np.nan))
    summary = pd.DataFrame(rows)
    summary.to_csv(f'{out_dir}/summary_zones.csv', index=False)

    pts = acc.points
    meta = dict(
        detector=geo.name, det_type=geo.det_type, sub_run=os.path.basename(sub_run_dir),
        feu=geo.feu_num, connectors=geo.feu_connectors,
        active_size=[float(v) for v in geo.active_size[:2]], z_mm=float(geo.center[2]),
        n_hits_mapped=acc.n_hits_raw, n_hits_selected=acc.n_hits_sel,
        n_events_with_cluster=int(n_ev), n_clusters=int(len(cl)),
        frac_events_xy=float(len(pts.dropna(subset=['x', 'y'])) / n_ev) if n_ev else np.nan,
        mixed_zone_cluster_frac=float(cl['zone_mixed'].mean()) if len(cl) else np.nan,
        dead_channels=dead.tolist(), n_dead=int(len(dead)), hot_channels=hot.tolist(),
        mpv_cluster_charge={k: float(v) for k, v in mpvs.items()},
        core_fit={v: {kk: (float(vv) if isinstance(vv, (int, float, np.floating)) else bool(vv))
                      for kk, vv in f.items()} for v, f in fits.items()},
        view_mode=geo.view_mode, axis_flip=geo.axis_flip,
        cuts=dict(min_amp=args.min_amp, min_signif=args.min_signif,
                  gap_factor=args.gap_factor, max_hits=args.max_hits))
    with open(f'{out_dir}/summary.json', 'w') as fh:
        json.dump(meta, fh, indent=2)
    if args.save_clusters:
        cl.to_parquet(f'{out_dir}/clusters.parquet')
        pts.to_parquet(f'{out_dir}/points.parquet')
    return dict(meta=meta, zones=summary, points=pts, clusters=cl)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run-dir', default=RUN_DIR_DEFAULT)
    ap.add_argument('--det-dir', default=DET_DIR_DEFAULT)
    ap.add_argument('--sub-run', default='nominal_00',
                    help="sub-run name, a regex, or 'all'")
    ap.add_argument('--det', default='both', choices=['both', 'front', 'back'])
    ap.add_argument('--out', default=OUT_DEFAULT)
    ap.add_argument('--min-amp', type=float, default=0.0)
    ap.add_argument('--min-signif', type=float, default=5.0)
    ap.add_argument('--gap-factor', type=float, default=1.05)
    ap.add_argument('--max-hits', type=int, default=3_000_000,
                    help='per detector per sub-run; 0 = all')
    ap.add_argument('--raw-view-mode', action='store_true',
                    help="take the strip map's connector order and channel order "
                         "as-is ('AB' on every view), instead of the measured "
                         "wiring in urw_lib.VIEW_MODE_DEFAULT -- for comparison "
                         'only; the raw order is wrong on all four views')
    ap.add_argument('--save-clusters', action='store_true')
    ap.add_argument('--verbose', action='store_true')
    args = ap.parse_args()
    if args.max_hits <= 0:
        args.max_hits = None

    run_json = os.path.join(args.run_dir, 'run_config.json')
    run_name = os.path.basename(args.run_dir.rstrip('/'))
    sub_runs = sorted(d for d in os.listdir(args.run_dir)
                      if os.path.isdir(os.path.join(args.run_dir, d)))
    if args.sub_run != 'all':
        pat = re.compile(args.sub_run)
        sub_runs = [s for s in sub_runs if pat.fullmatch(s) or s == args.sub_run]
    if not sub_runs:
        sys.exit(f'no sub-runs matching {args.sub_run!r} in {args.run_dir}')

    dets = URW_DETS if args.det == 'both' else \
        [d for d in URW_DETS if d.endswith(args.det)]

    all_zones, all_meta = [], []
    for sr in sub_runs:
        sr_dir = os.path.join(args.run_dir, sr)
        for name in dets:
            mode = {'x': 'AB', 'y': 'AB'} if args.raw_view_mode else None
            geo = UrwGeometry(name, run_json, args.det_dir, sub_run_name=sr,
                              view_mode=mode)
            tag = f'{name}  {run_name}/{sr}'
            print(f'== {tag}', flush=True)
            out_dir = os.path.join(args.out, name, run_name, sr)
            res = analyse(geo, sr_dir, args, out_dir, tag)
            if res is None:
                continue
            all_zones.append(res['zones'])
            all_meta.append(res['meta'])
            m = res['meta']
            print(f"   hits {m['n_hits_mapped']} -> sel {m['n_hits_selected']}, "
                  f"{m['n_clusters']} clusters in {m['n_events_with_cluster']} events, "
                  f"dead {m['n_dead']}, xy frac {m['frac_events_xy']:.3f}", flush=True)

    if all_zones:
        os.makedirs(args.out, exist_ok=True)
        z = pd.concat(all_zones, ignore_index=True)
        z.to_csv(os.path.join(args.out, f'summary_zones_{run_name}.csv'), index=False)
        pd.DataFrame(all_meta).to_json(
            os.path.join(args.out, f'summary_{run_name}.json'), orient='records', indent=2)
        print(f'\nwrote {args.out}/summary_zones_{run_name}.csv')


if __name__ == '__main__':
    main()
