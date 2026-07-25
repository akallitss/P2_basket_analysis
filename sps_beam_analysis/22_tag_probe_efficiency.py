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
import re
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
    # per-plane post-fit residual [mm]; the reference has 0 by construction.
    rmse = {name: (float(p.get('rmse_post') or 0.0))
            for name, p in al['planes'].items()}
    return tf, path, rmse


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


def recorded_mask(event_ids, rec_entry):
    """Boolean mask: which eventIds the FEU actually recorded. `rec_entry` is
    the (lo, hi, missing_set) triple from RunConfig.recorded_events."""
    lo, hi, missing = rec_entry
    ev = np.asarray(event_ids, dtype=np.int64)
    ok = (ev >= lo) & (ev <= hi)
    if missing:
        ok &= ~np.fromiter((int(e) in missing for e in ev), bool, len(ev))
    return ok


def eval_probe(cfg, probe, others, clusters, tf, probe_r, min_tag,
               fiducial=None, probe_rec=None):
    """Tag-and-probe for one probe plane. Returns (df_tag, n_tag, n_hit) where
    df_tag has columns pred_x, pred_y, hit (bool) for each tagged event.

    `others` is already restricted to the tag-ELIGIBLE planes (well-aligned
    ones). `fiducial`, if given as (xmin, xmax, ymin, ymax), keeps only events
    whose predicted impact point falls inside the probe's active area, so the
    efficiency is intrinsic (a track that misses the probe's pads is not counted
    as an inefficiency). `probe_rec`, if given as a (lo, hi, missing) recorded-
    events triple for the probe's FEU, adds a `recorded` column: tags in
    triggers the probe FEU never recorded are DAQ losses, not detector
    inefficiency, and are excluded from the corrected denominator.
    """
    # Vectorised over events (the per-event .loc loop was O(1e6) Python-level
    # and dominated the runtime). Each tag plane's single clusters are mapped
    # into the probe frame in one shot; predictions are averaged per event.
    inv_probe = tf[probe].inverse()
    parts = []
    for name in others:
        c = clusters[name]
        s = c[c['single']]
        if not len(s):
            continue
        rx, ry = tf[name].apply(s['x'].to_numpy(), s['y'].to_numpy())  # ->ref
        px, py = inv_probe.apply(rx, ry)                               # ->probe
        parts.append(pd.DataFrame({'eventId': s['eventId'].to_numpy(),
                                   'px': px, 'py': py}))
    empty = (pd.DataFrame(columns=['pred_x', 'pred_y', 'hit']), 0, 0, 0)
    if not parts:
        return empty
    allpred = pd.concat(parts, ignore_index=True)
    g = allpred.groupby('eventId').agg(px=('px', 'mean'), py=('py', 'mean'),
                                       n=('px', 'size'))
    need = max(min_tag, (len(others) // 2 + 1) if len(others) >= 3 else min_tag)
    g = g[g['n'] >= need]
    if not len(g):
        return empty

    n_outside = 0
    if fiducial is not None:
        x0, x1, y0, y1 = fiducial
        inside = ((g['px'] >= x0) & (g['px'] <= x1) &
                  (g['py'] >= y0) & (g['py'] <= y1))
        n_outside = int((~inside).sum())
        g = g[inside]
        if not len(g):
            return (pd.DataFrame(columns=['pred_x', 'pred_y', 'hit']),
                    0, 0, n_outside)

    # probe-hit test: nearest probe cluster (any) to the prediction, per event
    pc = clusters[probe][['eventId', 'x', 'y']]
    m = g.reset_index()[['eventId', 'px', 'py']].merge(pc, on='eventId',
                                                       how='left')
    m['d'] = np.hypot(m['x'] - m['px'], m['y'] - m['py'])
    dmin = m.groupby('eventId')['d'].min()
    hit = (dmin.reindex(g.index) <= probe_r).fillna(False).to_numpy()
    df = pd.DataFrame({'pred_x': g['px'].to_numpy(),
                       'pred_y': g['py'].to_numpy(), 'hit': hit})
    df['recorded'] = (recorded_mask(g.index.to_numpy(), probe_rec)
                      if probe_rec is not None else True)
    return df, len(df), int(df['hit'].sum()), n_outside


def fiducial_box(ct, margin):
    """(xmin, xmax, ymin, ymax) of the probe's pad centres, shrunk by `margin`
    mm on every side so a prediction near the edge still has room for its true
    pad to be inside the array."""
    px = ct['pad_cx'].to_numpy(); py = ct['pad_cy'].to_numpy()
    return (float(px.min()) + margin, float(px.max()) - margin,
            float(py.min()) + margin, float(py.max()) - margin)


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
    ap.add_argument('--subruns-glob', default=None,
                    help='comma-separated fnmatch patterns selecting sub_runs, '
                         "e.g. 'drift_*' or 'meshscan_*,nominal_*' -- for runs "
                         'that mix a drift arm and a mesh arm.')
    ap.add_argument('--probe-r', type=float, default=None,
                    help='probe-hit search radius [mm]. Default: '
                         '--probe-nsigma x the tag planes\' alignment residual.')
    ap.add_argument('--probe-nsigma', type=float, default=3.0,
                    help='probe radius in units of the alignment residual when '
                         '--probe-r is not given.')
    ap.add_argument('--cluster-r', type=float, default=15.0)
    ap.add_argument('--min-tag', type=int, default=None,
                    help='min OTHER planes with a clean cluster to tag '
                         '(default: run-config MIN_TAG).')
    ap.add_argument('--tag-max-rmse', type=float, default=40.0,
                    help='a plane may TAG only if its alignment post-fit '
                         'residual is below this [mm]; excludes planes that '
                         'do not align (e.g. a detector not tracking the beam).')
    ap.add_argument('--no-fiducial', action='store_true',
                    help='do not require the predicted impact inside the probe '
                         'active area (then efficiency is acceptance-relative).')
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
    prod_sub = 'scan'
    if args.subruns_glob:
        import fnmatch
        pats = [p.strip() for p in args.subruns_glob.split(',') if p.strip()]
        subruns = [s for s in subruns
                   if any(fnmatch.fnmatch(s, p) for p in pats)]
        # separate product dir per arm so drift and mesh curves don't overwrite
        prod_sub = 'scan_' + re.sub(r'[^A-Za-z0-9]+', '_',
                                    args.subruns_glob).strip('_')
    if not subruns:
        print('No sub_runs with combined hits (after --subruns-glob filter).')
        return
    if len(subruns) == 1 and not args.subruns_glob:
        prod_sub = subruns[0]
    # Which electrode does this run scan? Decided PER PROBE: a probe held at
    # fixed HV while another plane is scanned (P2_IN during the MID/OUT drift
    # leg) is a control curve — its x is the program's scanned value, and its
    # products are labelled by sub_run so points don't overwrite each other.
    # A run where nothing varies (stability run) plots vs sub_run index.
    own_axis = {d.name: cfg.scan_axis(subruns, d)[0] for d in dets}
    axis_det = next((d for d in dets if own_axis[d.name]), None)
    scan_ax, scan_label = (cfg.scan_axis(subruns, axis_det) if axis_det
                           else (None, 'sub_run index'))
    if scan_ax:
        print(f'  scan axis: {scan_label} (from {axis_det.name})')
    suffix = cfg.product_suffix(args.veto_sparks)
    caveat = ('intrinsic within the probe fiducial area; tag planes gated by '
              f'alignment residual < {args.tag_max_rmse:.0f} mm'
              if not args.no_fiducial else
              'acceptance-relative (no fiducial cut)')

    rows = []
    for sub in subruns:
        print(f'  sub_run {sub}:')
        tf, al_path, rmse = load_transforms(cfg, sub)
        if tf is None:
            print('    no alignment.json (run 21_telescope_align first) — '
                  'skipping')
            continue
        print(f'    alignment: {os.path.relpath(al_path, cfg.ANALYSIS_ROOT)}')
        # A plane may TAG only if it is itself well aligned; a plane that does
        # not track the beam gives garbage predictions and must not tag.
        eligible = {n for n, r in rmse.items() if r <= args.tag_max_rmse}
        excluded = [n for n in tf if n not in eligible]
        if excluded:
            print(f'    tag-ineligible (residual > {args.tag_max_rmse:.0f} mm): '
                  f'{", ".join(f"{n} ({rmse[n]:.0f} mm)" for n in excluded)}')
        clusters = {d.name: load_clusters(cfg, d, sub, args.cluster_r,
                                          args.veto_sparks) for d in dets}
        # per-FEU recorded-trigger sets (DAQ-overlap correction); None -> raw
        rec = cfg.recorded_events(sub)
        if rec is None:
            print('    (no recorded_events.npz -- efficiencies are not '
                  'DAQ-overlap corrected; run extract_recorded_events.py '
                  'on the DAQ host)')
        for probe in dets:
            others = [d.name for d in dets
                      if d.name != probe.name and d.name in eligible]
            if probe.name not in tf or not others:
                print(f'    {probe.name}: no eligible tag plane, skipping')
                continue
            probe_rec = rec.get(probe.feus[0]) if rec else None
            # probe radius: explicit, else N-sigma of the prediction residual.
            # The scale is the pairwise plane-to-plane residual, taken over the
            # WELL-ALIGNED planes in the pair (the probe's own residual is used
            # only if the probe itself aligned -- otherwise a broken probe like
            # P2_IN would blow the radius up with its garbage 132 mm fit).
            scales = [rmse[o] for o in others]
            if rmse.get(probe.name, 1e9) <= args.tag_max_rmse:
                scales.append(rmse[probe.name])
            res = max(scales) if scales else 0.0
            probe_r = (args.probe_r if args.probe_r is not None
                       else max(20.0, args.probe_nsigma * res))
            fid = (None if args.no_fiducial
                   else fiducial_box(cfg.channel_table(probe), probe_r))
            df, n_tag, n_hit, n_out = eval_probe(
                cfg, probe.name, others, clusters, tf, probe_r, min_tag,
                fiducial=fid, probe_rec=probe_rec)
            ax_p = own_axis[probe.name]
            hv = cfg.subrun_scan_hv(sub, probe, ax_p) if ax_p else None
            if hv is not None:
                pt, lbl = hv, f'{hv}V'
            elif scan_ax:
                pt = cfg.subrun_scan_hv(sub, axis_det, scan_ax)
                pt = pt if pt is not None else subruns.index(sub)
                lbl = sub
            else:
                pt, lbl = subruns.index(sub), sub
            eff = n_hit / n_tag if n_tag else np.nan
            lo, hi = clopper_pearson(n_hit, n_tag)
            # DAQ-overlap corrected: drop tags in triggers the probe FEU never
            # recorded (a hit implies recorded, so the numerator is unchanged).
            n_tag_rec = int(df['recorded'].sum()) if len(df) else 0
            eff_c = n_hit / n_tag_rec if n_tag_rec else np.nan
            lo_c, hi_c = clopper_pearson(n_hit, n_tag_rec)
            print(f'    probe {probe.name} (tags {"+".join(others)}, '
                  f'r={probe_r:.0f}mm): {n_tag} in-fiducial tagged '
                  f'({n_out} outside), {n_hit} found, eff = '
                  + (f'{eff:.3f}' if n_tag else 'n/a (no tags)')
                  + (f'  | recorded {n_tag_rec} '
                     f'({100 * (1 - n_tag_rec / n_tag):.1f}% DAQ loss) '
                     f'-> eff_corr = {eff_c:.3f} '
                     f'[+{hi_c - eff_c:.3f}/-{eff_c - lo_c:.3f}]'
                     if probe_rec is not None and n_tag else ''))
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
                             tags='+'.join(others), probe_r_mm=probe_r,
                             n_tag=n_tag, n_out_fiducial=n_out,
                             n_hit=n_hit, eff=eff,
                             eff_lo=lo, eff_hi=hi,
                             n_tag_recorded=(n_tag_rec if probe_rec is not None
                                             else None),
                             eff_corr=(eff_c if probe_rec is not None
                                       else None),
                             eff_corr_lo=(lo_c if probe_rec is not None
                                          else None),
                             eff_corr_hi=(hi_c if probe_rec is not None
                                          else None)))

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

    # efficiency vs scan point (per-probe own HV, or the program's scanned
    # value for fixed-HV control probes; sub_run index for stability runs)
    xcol, xlab = 'point', scan_label
    has_corr = df['eff_corr'].notna().any()
    fig, ax = plt.subplots(figsize=(8, 5))
    for probe, s in df.groupby('probe'):
        s = s.sort_values('point')
        if has_corr and s['eff_corr'].notna().any():
            yerr = np.vstack([s['eff_corr'] - s['eff_corr_lo'],
                              s['eff_corr_hi'] - s['eff_corr']])
            line = ax.errorbar(s[xcol], s['eff_corr'], yerr=yerr, fmt='o-',
                               capsize=4, lw=2, ms=7,
                               label=f'probe {probe} (DAQ-corr.)')
            ax.plot(s[xcol], s['eff'], 'o--', ms=4, lw=1, alpha=0.4,
                    color=line[0].get_color())
        else:
            yerr = np.vstack([s['eff'] - s['eff_lo'], s['eff_hi'] - s['eff']])
            ax.errorbar(s[xcol], s['eff'], yerr=yerr, fmt='o-', capsize=4,
                        lw=2, ms=7, label=f'probe {probe}')
    ax.set_xlabel(xlab)
    ax.set_ylabel('tag-probe efficiency'
                  + (' (solid = DAQ-overlap corrected, faint = raw)'
                     if has_corr else ''))
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
