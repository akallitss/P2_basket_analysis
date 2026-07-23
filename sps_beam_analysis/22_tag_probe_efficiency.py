#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
22_tag_probe_efficiency.py

Tag-and-probe efficiency of the P2 telescope stations (replaces the bench's
M3-referenced 06/11 efficiency). BANCO/ALPIDE provides no tracking, so this IS
the efficiency method: the other telescope planes tag a beam particle and the
plane under test is probed for a matching hit.

Definitions
-----------
For a PROBE plane, an event is TAGGED when at least --min-tag of the OTHER
planes carry a clean single cluster (sps_cluster: leading pad + pads within
--cluster-r, no pad outside). When >=3 planes exist a MAJORITY is required
(ceil(n_other/2) tags, and at least --min-tag). Each tagging plane's cluster is
mapped into the probe plane's frame through the alignment JSON written by
21_telescope_align (plane->reference->probe rigid transforms); the tag position
is the mean of those predictions. The probe is HIT when the probe plane has a
cluster within --probe-r mm of the tag position.

    efficiency = N(tagged AND probe hit) / N(tagged)

Clopper-Pearson (exact binomial) 68.27% intervals are used everywhere.

Honest caveat (baked into every plot subtitle)
----------------------------------------------
This is efficiency RELATIVE TO THE TAG SELECTION and only over the geometric
overlap acceptance of the tagging planes with the probe -- NOT an absolute
detector efficiency (which would need an external tracker). Regions the tag
never illuminates are simply absent from the denominator.

Products (<Analysis>/<probe_tag>/<run>/<sub_run-or-scan>/22_tag_probe_efficiency/):
  eff_map_<PROBE>_<pt><suffix>.png   per-pad efficiency over the overlap (+csv)
  tag_probe_efficiency<suffix>.png   efficiency vs HV / sub_run (all probes)
  tag_probe_efficiency<suffix>.csv   per-point, per-probe table

Usage:
  python3 22_tag_probe_efficiency.py [run_key] [--sub-run NAME] [--probe-r 8]
        [--cluster-r 15] [--min-tag 1] [--no-veto-sparks]
"""

import os
import json
import argparse

import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import beta

import sps_config as sc
import p2_mapping as pmap
import p2_io as p2io
import p2_sparks as ps
import sps_cluster as scl


def clopper_pearson(k, n, cl=0.6827):
    """Exact binomial (Clopper-Pearson) interval. Returns (lo, hi)."""
    if n == 0:
        return 0.0, 1.0
    a = 1.0 - cl
    lo = beta.ppf(a / 2, k, n - k + 1) if k > 0 else 0.0
    hi = beta.ppf(1 - a / 2, k + 1, n - k) if k < n else 1.0
    return float(lo), float(hi)


def _spark_veto(cfg, det, hv_csv):
    class _Shim:
        SPARK_CHANNEL = det.spark_channel
        SPARK_IMON_THR = cfg.SPARK_IMON_THR
        SPARK_GUARD_BEFORE = cfg.SPARK_GUARD_BEFORE
        SPARK_GUARD_AFTER = cfg.SPARK_GUARD_AFTER
        BURST_NPADS = cfg.BURST_NPADS
    return ps.SparkVeto.from_csv(hv_csv, _Shim)


def load_clusters(cfg, det, sub_run, cluster_r, veto_sparks):
    hits_dir = cfg.combined_hits_dir(sub_run)
    chunk0 = p2io.hit_files(hits_dir)[0]
    hv_csv = cfg.hv_monitor_csv(sub_run)
    t_min = scl.settle_t_min(hv_csv, det.spark_channel, chunk0)
    veto = (_spark_veto(cfg, det, hv_csv)
            if veto_sparks and os.path.isfile(hv_csv) and det.spark_channel
            else None)
    ct = cfg.channel_table(det)
    return scl.stream_event_clusters(hits_dir, ct, cluster_r,
                                     min_amp=cfg.MIN_AMP, veto=veto,
                                     t_min=t_min)


def load_transforms(cfg, sub_run):
    """{det_name: RigidTransform plane->reference} from the alignment JSON.
    Falls back to any alignment.json in the run if this sub_run has none."""
    path = cfg.alignment_json(sub_run)
    if not os.path.isfile(path):
        import glob
        cands = glob.glob(os.path.join(cfg.ANALYSIS_ROOT, sc.TELESCOPE_TAG,
                                       cfg.RUN, '*', '21_telescope_align',
                                       'alignment.json'))
        if not cands:
            return None, None
        path = sorted(cands)[0]
    with open(path) as f:
        al = json.load(f)
    tf = {name: scl.RigidTransform.from_dict(p)
          for name, p in al['planes'].items()}
    return tf, path


def predict_in_probe(tag_positions, probe_name, tf):
    """Map a list of (name, x, y) tag clusters into the probe plane frame and
    average. tf[name] maps plane->reference; probe mapping is ref->probe."""
    inv_probe = tf[probe_name].inverse()
    xs, ys = [], []
    for name, x, y in tag_positions:
        rx, ry = tf[name].apply(x, y)          # plane -> reference
        px, py = inv_probe.apply(rx, ry)       # reference -> probe
        xs.append(float(px)); ys.append(float(py))
    return float(np.mean(xs)), float(np.mean(ys))


def eval_probe(cfg, probe, others, clusters, tf, probe_r, min_tag):
    """Tag-and-probe for one probe plane. Returns (df_tag, n_tag, n_hit) where
    df_tag has columns pred_x, pred_y, hit (bool) for each tagged event."""
    # single clean clusters per plane, indexed by eventId
    single = {name: c[c['single']].set_index('eventId')
              for name, c in clusters.items()}
    probe_all = clusters[probe].set_index('eventId')  # any cluster on probe

    n_other = len(others)
    need = max(min_tag, (n_other // 2 + 1) if n_other >= 3 else min_tag)

    # union of eventIds seen by any tagging plane's single clusters
    ev_union = set()
    for name in others:
        ev_union |= set(single[name].index.tolist())

    rows = []
    for ev in ev_union:
        tags = []
        for name in others:
            if ev in single[name].index:
                r = single[name].loc[ev]
                tags.append((name, float(r['x']), float(r['y'])))
        if len(tags) < need:
            continue
        px, py = predict_in_probe(tags, probe, tf)
        hit = False
        if ev in probe_all.index:
            pr = probe_all.loc[ev]
            # a probe event can appear once (Series); guard for duplicates
            pxv = np.atleast_1d(pr['x']); pyv = np.atleast_1d(pr['y'])
            d = np.hypot(pxv - px, pyv - py)
            hit = bool(np.nanmin(d) <= probe_r)
        rows.append((px, py, hit))
    df = pd.DataFrame(rows, columns=['pred_x', 'pred_y', 'hit'])
    return df, len(df), int(df['hit'].sum()) if len(df) else 0


def plot_eff_map(df, ct, probe, lbl, sub, caveat, out_png, out_csv,
                 min_pad_tags=5):
    """Per-pad efficiency over the overlap: assign each tag to the nearest
    probe pad, efficiency per pad (Clopper-Pearson counts stored in the CSV)."""
    from matplotlib.collections import PolyCollection
    pads = ct.drop_duplicates('channel_id')[['channel_id', 'pad_cx',
                                              'pad_cy']].reset_index(drop=True)
    if not len(df):
        pd.DataFrame(columns=['channel_id', 'n_tag', 'n_hit', 'eff']).to_csv(
            out_csv, index=False)
        return
    px = pads['pad_cx'].to_numpy(); py = pads['pad_cy'].to_numpy()
    # nearest pad for each tag prediction. A dense (n_tag x n_pad) distance
    # matrix is O(n_tag*n_pad) memory -- ~19 GB at 4.7M tags x 512 pads -- so
    # use a KDTree over the (few hundred) pads: O(n_tag log n_pad), MB not GB.
    tx = df['pred_x'].to_numpy(); ty = df['pred_y'].to_numpy()
    from scipy.spatial import cKDTree
    _, idx = cKDTree(np.column_stack([px, py])).query(
        np.column_stack([tx, ty]), k=1)
    pads['n_tag'] = np.bincount(idx, minlength=len(pads))
    pads['n_hit'] = np.bincount(idx, weights=df['hit'].to_numpy(float),
                                minlength=len(pads)).astype(int)
    pads['eff'] = np.where(pads['n_tag'] > 0, pads['n_hit'] / pads['n_tag'],
                           np.nan)
    pads['eff'] = pads['eff'].where(pads['n_tag'] >= min_pad_tags)
    pads.to_csv(out_csv, index=False)

    fig, ax = plt.subplots(figsize=(7.6, 6.2))
    if pmap.has_tile_geometry(ct):
        tpads, verts = pmap.pad_tiles(ct)
        e = (pads.set_index('channel_id')['eff']
             .reindex(tpads['channel_id']).to_numpy(float))
        good = np.isfinite(e)
        ax.add_collection(PolyCollection(verts[~good], facecolors='0.92',
                                         edgecolors='0.7', linewidths=0.3))
        if good.any():
            pc = PolyCollection(verts[good], array=e[good], cmap='viridis',
                                edgecolors='face', linewidths=0.2,
                                clim=(0, 1))
            ax.add_collection(pc)
            fig.colorbar(pc, ax=ax, label='tag-probe efficiency')
        ax.autoscale_view(); ax.set_aspect('equal')
    else:
        s = ax.scatter(pads['pad_cx'], pads['pad_cy'], c=pads['eff'],
                       cmap='viridis', vmin=0, vmax=1, s=14)
        fig.colorbar(s, ax=ax, label='tag-probe efficiency')
        ax.set_aspect('equal')
    ax.set_xlabel('pad_cx [mm]'); ax.set_ylabel('pad_cy [mm]')
    ax.set_title(f'{probe} tag-probe efficiency map — {lbl} ({sub})\n{caveat}',
                 fontsize=9)
    fig.tight_layout()
    fig.savefig(out_png, dpi=160, bbox_inches='tight')
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description='Tag-and-probe telescope '
                                             'efficiency.')
    ap.add_argument('run_key', nargs='?', default=sc.DEFAULT_RUN)
    ap.add_argument('--sub-run', default=None,
                    help='sub_run (default: all discovered sub_runs).')
    ap.add_argument('--probe-r', type=float, default=8.0,
                    help='probe-hit search radius around the tag position [mm].')
    ap.add_argument('--cluster-r', type=float, default=15.0)
    ap.add_argument('--min-tag', type=int, default=None,
                    help='min OTHER planes with a clean cluster to tag '
                         '(default: run-config MIN_TAG).')
    ap.add_argument('--min-pad-tags', type=int, default=5,
                    help='min tags in a pad for it to get a map entry.')
    ap.add_argument('--veto-sparks', action=argparse.BooleanOptionalAction,
                    default=True)
    args = ap.parse_args()

    cfg = sc.get_config(args.run_key)
    print(cfg)
    # Tag-and-probe runs on the planes that have a pad map (the uRWELL
    # references have none); alignment (21) already ran on the same set.
    dets = cfg.mappable_detectors()
    skipped = [d.name for d in cfg.detectors() if d not in dets]
    if skipped:
        print(f'  (no pad map, skipped: {", ".join(skipped)})')
    if len(dets) < 2:
        print('Need >=2 mappable stations for tag-and-probe.')
        return
    min_tag = args.min_tag if args.min_tag is not None else cfg.MIN_TAG
    subruns = [args.sub_run] if args.sub_run else cfg.find_subruns()
    if not subruns:
        print('No sub_runs with combined hits.')
        return
    prod_sub = 'scan' if len(subruns) > 1 else subruns[0]
    suffix = cfg.product_suffix(args.veto_sparks)
    caveat = ('efficiency relative to the tag selection / overlap acceptance '
              'only (not absolute)')

    rows = []
    for sub in subruns:
        print(f'  sub_run {sub}:')
        tf, al_path = load_transforms(cfg, sub)
        if tf is None:
            print('    no alignment.json (run 21_telescope_align first) — '
                  'skipping')
            continue
        print(f'    alignment: {os.path.relpath(al_path, cfg.ANALYSIS_ROOT)}')
        clusters = {d.name: load_clusters(cfg, d, sub, args.cluster_r,
                                          args.veto_sparks) for d in dets}
        for probe in dets:
            others = [d.name for d in dets if d.name != probe.name]
            if probe.name not in tf or any(o not in tf for o in others):
                print(f'    {probe.name}: missing transform, skipping')
                continue
            df, n_tag, n_hit = eval_probe(cfg, probe.name, others, clusters, tf,
                                          args.probe_r, min_tag)
            hv = cfg.subrun_mesh_hv(sub, probe)
            pt = hv if hv is not None else subruns.index(sub)
            lbl = f'{hv}V' if hv is not None else sub
            eff = n_hit / n_tag if n_tag else np.nan
            lo, hi = clopper_pearson(n_hit, n_tag)
            print(f'    probe {probe.name}: {n_tag} tagged, {n_hit} found, '
                  f'eff = '
                  + (f'{eff:.3f} [+{hi-eff:.3f}/-{eff-lo:.3f}]'
                     if n_tag else 'n/a (no tags)'))
            out_dir = cfg.out_dir(probe.det_tag, prod_sub,
                                  '22_tag_probe_efficiency')
            plot_eff_map(df, cfg.channel_table(probe), probe.name, lbl, sub,
                         caveat,
                         os.path.join(out_dir,
                                      f'eff_map_{probe.det_tag}_{lbl}{suffix}.png'),
                         os.path.join(out_dir,
                                      f'eff_map_{probe.det_tag}_{lbl}{suffix}.csv'),
                         min_pad_tags=args.min_pad_tags)
            rows.append(dict(sub_run=sub, point=pt, hv=hv,
                             probe=probe.name, probe_tag=probe.det_tag,
                             n_tag=n_tag, n_hit=n_hit, eff=eff,
                             eff_lo=lo, eff_hi=hi))

    if not rows:
        print('No efficiency points produced.')
        return
    df = pd.DataFrame(rows)
    # write one CSV per probe under its own tree, plus a combined print
    for probe_tag, sub_df in df.groupby('probe_tag'):
        out_dir = cfg.out_dir(probe_tag, prod_sub, '22_tag_probe_efficiency')
        sub_df.sort_values('point').to_csv(
            os.path.join(out_dir, f'tag_probe_efficiency{suffix}.csv'),
            index=False)

    # efficiency vs HV / point, all probes
    has_hv = df['hv'].notna().any()
    xcol = 'hv' if has_hv else 'point'
    xlab = 'mesh HV [V]' if has_hv else 'sub_run index'
    fig, ax = plt.subplots(figsize=(8, 5))
    for probe, s in df.groupby('probe'):
        s = s.sort_values('point')
        yerr = np.vstack([s['eff'] - s['eff_lo'], s['eff_hi'] - s['eff']])
        ax.errorbar(s[xcol], s['eff'], yerr=yerr, fmt='o-', capsize=4, lw=2,
                    ms=7, label=f'probe {probe}')
    ax.set_xlabel(xlab); ax.set_ylabel('tag-probe efficiency')
    ax.set_ylim(0, 1.02); ax.grid(True, alpha=0.3); ax.legend(fontsize=9)
    ax.set_title(f'Tag-probe efficiency vs {xlab} — {cfg.RUN}\n{caveat}',
                 fontsize=10)
    fig.tight_layout()
    for probe_tag in df['probe_tag'].unique():
        out_dir = cfg.out_dir(probe_tag, prod_sub, '22_tag_probe_efficiency')
        fig.savefig(os.path.join(out_dir,
                                 f'tag_probe_efficiency{suffix}.png'),
                    dpi=170, bbox_inches='tight')
    plt.close(fig)
    print(f'\nWrote tag-probe products for probes: '
          f'{sorted(df["probe"].unique())}')


if __name__ == '__main__':
    main()
