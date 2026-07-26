#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
urw_p2_efficiency.py -- P2 residuals and efficiency referenced to uRWELL tracks.

The P2 telescope's own stage 22 (tag-and-probe) can only measure a P2 station
against the OTHER P2 stations, so it is limited by their 12 mm pads and it is
honest only about efficiency *relative to the tag selection*.  The two EIC
uRWELL references bracket the three stations (z = 0 and 1370 mm, with P2_IN /
MID / OUT at 320 / 630 / 940 mm) and are strip detectors with ~1 mm pitch, so
once they are mapped correctly they give an external track that points at every
P2 plane to about 1 mm -- three times better than a P2 pad.  That turns the
measurement into a real tracking-referenced efficiency.

What this stage does, per sub_run:

  1. Builds uRWELL points: leading cluster per view per plane, from the FEU-1
     hits files (urw_lib).  Front and back share an eventId because both sit on
     FEU 1, so no cross-FEU synchronisation is needed.
  2. Aligns back to front (linear, per axis) and cuts on their agreement.  That
     cut is the track chi^2: with the correct strip mapping the core of
     back - front is ~0.9 mm wide, and the ~10 % outside it is pile-up and
     mis-assigned leading clusters.
  3. Interpolates to each P2 plane's survey z.
  4. Streams that station's combined hits into a leading cluster per event
     (plus the leading pad of any SECOND cluster, so a second particle in the
     event cannot fake a miss), matches on eventId, and fits the uRWELL frame
     onto the P2 pad frame with a free affine (see fit_frame: the affine is
     applied, the rigid fit is reported as the check that it stayed orthogonal).
  5. Residuals, and efficiency = N(track predicted in the pads AND a P2 cluster
     within --probe-r) / N(track predicted in the pads), with Clopper-Pearson
     68.27 % intervals, as a number, a profile and a 2D map.

Why the affine is a real check, not decoration
----------------------------------------------
The uRWELL -> P2 map has to be a proper rotation: the two detectors are viewed
from the same side along the same beam, so no reflection is allowed and no
scale is allowed (both frames are in mm).  Fitting a free 2x2 matrix and
finding it orthogonal with det = +1 is therefore a genuine test that both strip
maps and the pad map are right.  Measured on highstat_eff_1: -60 deg,
singular values 0.99-1.02, det +1.00 to +1.02 at all three stations, and the
1-2 % departure is beam divergence, not a map error (fit_frame).

Note on the front->back scale removed in build_tracks: because the beam
diverges, back = 1.027 * front in y is physical, and dividing it out puts the
track projection in a "parallel beam" frame.  The per-station affine then puts
the magnification back, so the two steps are consistent - and the front-back
agreement cut stays a pure noise cut rather than one that grows with position.

Denominator hygiene
-------------------
An efficiency is only as good as the events it does NOT count.  Dropped from
the denominator, and each reported:
  * triggers the probe's FEU never recorded (recorded_events.npz), which are
    not the detector's fault;
  * events inside an HV spark window for that station, and events before its
    mesh had settled -- the same p2_sparks veto stage 21/22 use, but applied by
    TIME to the track list so the veto removes numerator and denominator
    together rather than only killing the P2 hits;
  * tracks whose prediction is not inside the pad plane (nearest pad centre
     further than --fid-r), so acceptance edges are not counted as inefficiency.
Dead or dropped connectors are deliberately NOT removed: they stay in the
denominator and show up as holes in the efficiency map.

Usage
  unset PYTHONPATH
  /local/home/banco/DAQ_Control_Dream_Beam/.venv/bin/python urw_p2_efficiency.py \
      --run highstat_eff_1 --sub-run beam_commissioning_00 --out out_eff
"""
import os
import sys
import json
import glob
import argparse

import numpy as np
import pandas as pd

# this package, then sps_beam_analysis for sps_config, whose setup_paths()
# adds the shared core in cosmic_bench_analysis (p2_io, p2_mapping, ...)
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.dirname(_HERE))
import sps_config as _sc  # noqa: E402
_sc.setup_paths()

import urw_lib as U            # noqa: E402
import sps_config as sc        # noqa: E402
import sps_cluster as scl      # noqa: E402
import p2_mapping as pmap      # noqa: E402
import p2_sparks as ps         # noqa: E402

DATA_ROOT = '/local/home/banco/P2_data/TB_July2026_H4/runs'
FRONT, BACK = 'EIC_uRWELL_front', 'EIC_uRWELL_back'


# --------------------------------------------------------------------------- #
# uRWELL tracks
# --------------------------------------------------------------------------- #
def urw_points(det, run_json, sub_dir, max_chunks=None):
    geo = U.UrwGeometry(det, run_json, sub_run_name=os.path.basename(sub_dir))
    files = U.feu_hit_files(sub_dir, feu=geo.feu_num)
    if max_chunks:
        files = files[:max_chunks]
    parts = []
    # feu= is not optional: where the per-FEU hits_root files have been deleted
    # (all 23 sub-runs of drift_mesh_scan_1) feu_hit_files falls back to the
    # combined files, which carry all four FEUs, and without the filter the P2
    # channels would be silently mapped onto uRWELL strips.
    for chunk in U.iter_hits(files, progress=False, feu=geo.feu_num):
        cl = U.cluster_hits(chunk, geo)
        if len(cl):
            parts.append(U.leading_points(cl))
    pts = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    return geo, pts


def robust_line(x, y, nsig=2.5, iters=6):
    keep = np.isfinite(x) & np.isfinite(y)
    slope = offset = np.nan
    for _ in range(iters):
        if keep.sum() < 50:
            break
        slope, offset = np.polyfit(x[keep], y[keep], 1)
        r = y - (slope * x + offset)
        new = np.isfinite(x) & np.isfinite(y) & (np.abs(r) < nsig * r[keep].std())
        if new.sum() == keep.sum():
            break
        keep = new
    return float(slope), float(offset)


def build_tracks(run_json, sub_dir, max_chunks, track_cut):
    """Two-point uRWELL tracks with a front-back agreement cut.

    Returns (tracks, info).  `tracks` carries eventId, the front and
    back-in-front-frame positions, the slope, and t_ns for the spark veto.
    """
    geo_f, pf = urw_points(FRONT, run_json, sub_dir, max_chunks)
    geo_b, pb = urw_points(BACK, run_json, sub_dir, max_chunks)
    dz = float(geo_b.center[2] - geo_f.center[2])

    ev = pf[['eventId', 'x', 'y', 'x_size', 'y_size', 't_ns']].rename(
        columns={'x': 'xf', 'y': 'yf', 'x_size': 'nsxf', 'y_size': 'nsyf'}).merge(
        pb[['eventId', 'x', 'y', 'x_size', 'y_size']].rename(
            columns={'x': 'xb', 'y': 'yb', 'x_size': 'nsxb', 'y_size': 'nsyb'}),
        on='eventId')
    n_pair = len(ev)
    ev = ev.dropna(subset=['xf', 'yf', 'xb', 'yb']).reset_index(drop=True)

    info = {'n_front': int(len(pf)), 'n_back': int(len(pb)),
            'n_matched': int(n_pair), 'n_4coord': int(len(ev)), 'dz_mm': dz,
            'front_back': {}}
    for ax in ('x', 'y'):
        s, o = robust_line(ev[f'{ax}f'].to_numpy(), ev[f'{ax}b'].to_numpy())
        ev[f'{ax}b_c'] = (ev[f'{ax}b'] - o) / s
        ev[f'd{ax}'] = ev[f'{ax}b_c'] - ev[f'{ax}f']
        q = np.percentile(ev[f'd{ax}'], [25, 75])
        info['front_back'][ax] = {
            'slope': s, 'offset_mm': o,
            'sigma_iqr_mm': float((q[1] - q[0]) / 1.349),
            'rms_mm': float(ev[f'd{ax}'].std())}
    good = (ev['dx'].abs() < track_cut) & (ev['dy'].abs() < track_cut)
    info['track_cut_mm'] = float(track_cut)
    info['n_good_track'] = int(good.sum())
    info['good_track_frac'] = float(good.mean())
    ev = ev[good].reset_index(drop=True)
    for ax in ('x', 'y'):
        ev[f'slope_{ax}'] = ev[f'd{ax}'] / dz
    return ev, info, dz


def project(tracks, z, dz):
    t = z / dz
    return np.column_stack([
        tracks['xf'] + t * (tracks['xb_c'] - tracks['xf']),
        tracks['yf'] + t * (tracks['yb_c'] - tracks['yf'])])


# --------------------------------------------------------------------------- #
# P2 clusters
# --------------------------------------------------------------------------- #
def p2_clusters(files, channel_table, cluster_r, min_amp=0.0):
    """Leading cluster per event, plus the leading pad of any second cluster.

    Same reduction as sps_cluster.stream_event_clusters (leading pad + mapped
    pads within cluster_r, charge-weighted centroid, `single` = nothing
    outside), over an explicit file list so --max-chunks can bound a first
    pass, and with (x2, y2) added: without a second candidate a second particle
    in the same trigger pulls the leading cluster away from the track and fakes
    a miss.
    """
    import uproot
    branches = ['eventId', 'channel', 'amplitude', 'feu', 'trigger_timestamp_ns']
    feus = set(channel_table.attrs['feus'])
    parts = []
    for fp in files:
        with uproot.open(fp) as f:
            arrs = f['hits'].arrays(branches, library='np')
        df = pd.DataFrame(arrs)
        del arrs
        df = df[df['feu'].isin(feus)]
        if min_amp > 0:
            df = df[df['amplitude'] >= min_amp]
        h = pmap.attach_pads_to_hits(df, channel_table)
        h = h[h['mapped'] & h['pad_cx'].notna()]
        del df
        if not len(h):
            continue
        lead = h.loc[h.groupby('eventId')['amplitude'].idxmax(),
                     ['eventId', 'amplitude', 'pad_cx', 'pad_cy', 'channel_id',
                      'trigger_timestamp_ns']].rename(
            columns={'amplitude': 'a_lead', 'pad_cx': 'lx', 'pad_cy': 'ly',
                     'channel_id': 'lead_pad'})
        npad = h.groupby('eventId').size().rename('n_pad')
        h = h.merge(lead[['eventId', 'lx', 'ly']], on='eventId')
        near = ((h['pad_cx'] - h['lx']) ** 2 +
                (h['pad_cy'] - h['ly']) ** 2) <= cluster_r ** 2
        hc, ho = h[near], h[~near]
        w = hc['amplitude'].clip(lower=0).astype(np.float64)
        g = pd.DataFrame({'eventId': hc['eventId'], '_wx': w * hc['pad_cx'],
                          '_wy': w * hc['pad_cy'], '_w': w})
        agg = g.groupby('eventId').agg(_wx=('_wx', 'sum'), _wy=('_wy', 'sum'),
                                       _w=('_w', 'sum'), n_clus=('_w', 'size'))
        e = pd.DataFrame({'x': agg['_wx'] / agg['_w'], 'y': agg['_wy'] / agg['_w'],
                          'q': agg['_w'], 'n_clus': agg['n_clus']})
        e = e.join(npad).join(lead.set_index('eventId')[
            ['a_lead', 'lead_pad', 'trigger_timestamp_ns']])
        if len(ho):
            second = ho.loc[ho.groupby('eventId')['amplitude'].idxmax(),
                            ['eventId', 'pad_cx', 'pad_cy']].rename(
                columns={'pad_cx': 'x2', 'pad_cy': 'y2'}).set_index('eventId')
            e = e.join(second)
        else:
            e['x2'] = np.nan
            e['y2'] = np.nan
        parts.append(e.reset_index())
        del h, hc, ho
    if not parts:
        return pd.DataFrame(columns=['eventId', 'x', 'y', 'q', 'n_clus', 'n_pad',
                                     'a_lead', 'lead_pad', 'x2', 'y2', 't'])
    ev = pd.concat(parts, ignore_index=True)
    # an event straddling a chunk boundary appears twice: keep the larger half
    ev = (ev.sort_values('n_pad').drop_duplicates('eventId', keep='last')
          .reset_index(drop=True))
    ev['single'] = ev['n_pad'] == ev['n_clus']
    ev['t'] = ev['trigger_timestamp_ns'].astype(np.int64) / 1e9
    ev['eventId'] = ev['eventId'].astype(np.int64)
    return ev.drop(columns='trigger_timestamp_ns')


# --------------------------------------------------------------------------- #
# uRWELL frame -> P2 pad frame
# --------------------------------------------------------------------------- #
class Affine:
    """X = A @ x + t, with A a free 2x2."""

    def __init__(self, A, t):
        self.A = np.asarray(A, float)
        self.t = np.asarray(t, float)

    def apply(self, x, y):
        p = np.column_stack([np.asarray(x, float), np.asarray(y, float)])
        q = p @ self.A.T + self.t
        return q[:, 0], q[:, 1]

    def to_dict(self):
        return {'A': self.A.tolist(), 't': self.t.tolist()}


def fit_frame(src, dst, win):
    """Fit src (uRWELL frame) -> dst (P2 pad frame).

    The transform APPLIED downstream is the free affine.  The rigid fit is also
    done and reported, as the check that the geometry makes sense: two planes
    seen from the same side along the same beam can only differ by a rotation
    and a translation, so the affine must come out orthogonal.

    It comes out orthogonal to 1-2 %, and that 1-2 % is understood: it is an
    anisotropic magnification growing linearly with z (explore6_divergence.py).
    The stretch tensor is diagonal in the uRWELL basis to better than 0.2 % -
    no shear - and its y term extrapolates to the independently measured
    front->back slope at dz = 1370 mm to within 0.16 %, giving a virtual source
    41.7 m upstream in y and 224 m in x.  That is a divergent beam, focused in
    one plane and not the other, not a detector distortion.  Using the affine
    keeps that known optics term out of the residuals; the rigid rmse and the
    singular values below are what tell you it stayed small.
    """
    A, t = np.eye(2), np.zeros(2)
    keep = np.ones(len(src), bool)
    for i in range(4):
        M = np.column_stack([src[keep], np.ones(keep.sum())])
        sol, *_ = np.linalg.lstsq(M, dst[keep], rcond=None)
        A, t = sol[:2].T, sol[2]
        if i == 0:
            continue          # first pass starts from identity: window is blind
        r = dst - (src @ A.T + t)
        keep = (np.abs(r[:, 0]) < win) & (np.abs(r[:, 1]) < win)
    u, s, vt = np.linalg.svd(A)
    rot = u @ vt
    stretch = vt.T @ np.diag(s) @ vt      # symmetric, in the uRWELL basis
    tf, rmse = scl.fit_rigid(src[keep, 0], src[keep, 1], dst[keep, 0], dst[keep, 1])
    return dict(
        affine_A=A.tolist(), affine_t=t.tolist(),
        singular_values=[float(v) for v in s], det=float(np.linalg.det(A)),
        affine_rotation_deg=float(np.degrees(np.arctan2(rot[1, 0], rot[0, 0]))),
        stretch_xx=float(stretch[0, 0]), stretch_yy=float(stretch[1, 1]),
        stretch_xy=float(stretch[0, 1]),
        rigid=tf.to_dict(), rigid_rmse_mm=float(rmse),
        n_fit=int(keep.sum()), frac_in_window=float(keep.mean())), Affine(A, t), keep


# --------------------------------------------------------------------------- #
# efficiency
# --------------------------------------------------------------------------- #
def clopper_pearson(k, n, cl=0.6827):
    from scipy.stats import beta
    if n == 0:
        return 0.0, 1.0
    a = 1.0 - cl
    lo = beta.ppf(a / 2, k, n - k + 1) if k > 0 else 0.0
    hi = beta.ppf(1 - a / 2, k + 1, n - k) if k < n else 1.0
    return float(lo), float(hi)


def recorded_mask(event_ids, rec_entry):
    lo, hi, missing = rec_entry
    ev = np.asarray(event_ids, dtype=np.int64)
    ok = (ev >= lo) & (ev <= hi)
    if missing:
        ok &= ~np.isin(ev, np.fromiter(missing, np.int64, len(missing)))
    return ok


def run_station(cfg, det, sub_run, tracks, dz, files, args, veto, t_min):
    """Everything for one probe station.  Returns (summary dict, DataFrame)."""
    from scipy.spatial import cKDTree
    ct = cfg.channel_table(det)
    pads = ct[ct['mapped']].drop_duplicates('channel_id')
    tree = cKDTree(pads[['pad_cx', 'pad_cy']].to_numpy())

    p2 = p2_clusters(files, ct, args.cluster_r, min_amp=cfg.MIN_AMP)
    print(f'  {det.name}: {len(p2)} events with a P2 cluster')

    # -- denominator hygiene, in the order the cuts are reported -------------- #
    tr = tracks.copy()
    n0 = len(tr)
    rec = (cfg.recorded_events(sub_run) or {}).get(det.feus[0])
    if rec is not None:
        tr = tr[recorded_mask(tr['eventId'], rec)]
    n_rec = len(tr)
    if t_min > 0:
        tr = tr[tr['t_ns'] / 1e9 >= t_min]
    n_settle = len(tr)
    if veto is not None and veto.intervals:
        tr = tr[veto.event_mask(tr['t_ns'].to_numpy())]
    n_veto = len(tr)

    src = project(tr, det.z, dz)

    # -- frame fit, on the matched subset ------------------------------------ #
    j = pd.DataFrame({'eventId': tr['eventId'].to_numpy(),
                      'px': src[:, 0], 'py': src[:, 1]}).merge(
        p2[['eventId', 'x', 'y', 'x2', 'y2', 'n_clus', 'n_pad', 'q', 'single']],
        on='eventId', how='inner')
    if len(j) < args.min_events:
        print(f'    only {len(j)} tracks match a P2 cluster '
              f'(need {args.min_events}) - skipping this station')
        return None, None, None
    fit, tf, keep = fit_frame(j[['px', 'py']].to_numpy(),
                              j[['x', 'y']].to_numpy(), args.match_win)
    print(f'    frame: rot {fit["affine_rotation_deg"]:+.3f} deg  '
          f'det {fit["det"]:+.4f}  sing {fit["singular_values"][0]:.4f}/'
          f'{fit["singular_values"][1]:.4f}  rigid rmse {fit["rigid_rmse_mm"]:.2f} mm'
          f'  stretch {fit["stretch_xx"]:.4f}/{fit["stretch_yy"]:.4f} '
          f'(shear {fit["stretch_xy"]:+.4f})')

    # -- residuals ----------------------------------------------------------- #
    fx, fy = tf.apply(j['px'].to_numpy(), j['py'].to_numpy())
    j['ex'], j['ey'] = fx, fy                       # track in the P2 pad frame
    j['rx'], j['ry'] = j['x'] - fx, j['y'] - fy
    core = (j['rx'].abs() < args.match_win) & (j['ry'].abs() < args.match_win)
    resid = {}
    for ax in ('x', 'y'):
        d = j.loc[core, f'r{ax}']
        q = np.percentile(d, [25, 75])
        resid[ax] = {'median_mm': float(d.median()), 'rms_mm': float(d.std()),
                     'sigma_iqr_mm': float((q[1] - q[0]) / 1.349)}
    print(f'    residual rms  dx {resid["x"]["rms_mm"]:.2f}  '
          f'dy {resid["y"]["rms_mm"]:.2f} mm  (n={int(core.sum())})')

    # -- efficiency ---------------------------------------------------------- #
    ex, ey = tf.apply(src[:, 0], src[:, 1])
    d_pad, _ = tree.query(np.column_stack([ex, ey]))
    fid = d_pad < args.fid_r
    n_fid = int(fid.sum())

    hit = pd.DataFrame({'eventId': tr['eventId'].to_numpy()[fid],
                        'ex': ex[fid], 'ey': ey[fid],
                        # the same tracks in the uRWELL's own frame, so a
                        # feature can be attributed to the reference or to the
                        # probe: a uRWELL artefact sits at a fixed (ux, uy) for
                        # every station, a P2 defect at a fixed (ex, ey) for one
                        'ux': src[fid, 0], 'uy': src[fid, 1]}).merge(
        p2[['eventId', 'x', 'y', 'x2', 'y2', 'n_clus', 'q', 'a_lead']],
        on='eventId', how='left')
    d1 = np.hypot(hit['x'] - hit['ex'], hit['y'] - hit['ey'])
    d2 = np.hypot(hit['x2'] - hit['ex'], hit['y2'] - hit['ey'])
    hit['dmin'] = np.fmin(d1.fillna(np.inf), d2.fillna(np.inf))
    hit['found'] = hit['dmin'] < args.probe_r
    k, n = int(hit['found'].sum()), len(hit)
    lo, hi = clopper_pearson(k, n)
    # a miss is either "the station recorded nothing at all for this trigger"
    # (a real inefficiency) or "it fired somewhere else" (a real inefficiency
    # only if the track is right; otherwise a reference artefact)
    n_empty = int(hit['x'].isna().sum())
    n_far = n - k - n_empty
    print(f'    efficiency  {k}/{n} = {k / max(n, 1):.4f}  '
          f'[{lo:.4f}, {hi:.4f}]  (68.27% CP)')
    print(f'      misses: {n_empty} with no P2 hit at all, {n_far} with a '
          f'cluster further than {args.probe_r:g} mm')

    # Efficiency versus the matching radius: the systematic on the number above.
    # A real hit lands within about a pad; beyond that the curve only picks up
    # accidental matches, so the slope of the plateau IS the accidental rate.
    dm = hit['dmin'].to_numpy()
    scan = [(float(r), float(np.mean(dm < r))) for r in
            (5, 8, 10, 12, 15, 18, 20, 25, 30, 40)]
    i15 = [r for r, _ in scan].index(15.0)
    slope = ((scan[-1][1] - scan[i15][1]) / (scan[-1][0] - scan[i15][0])
             if len(scan) > i15 + 1 else 0.0)
    print('      eff vs probe r: ' +
          '  '.join(f'{r:g}:{e:.4f}' for r, e in scan) +
          f'   -> accidental slope {slope * 10:.4f} per 10 mm')

    summary = dict(
        station=det.name, z_mm=det.z, feu=det.feus[0],
        n_tracks=n0, n_after_recorded=n_rec, n_after_settle=n_settle,
        n_after_spark_veto=n_veto, n_in_fiducial=n_fid,
        fid_r_mm=args.fid_r, probe_r_mm=args.probe_r,
        cluster_r_mm=args.cluster_r, match_win_mm=args.match_win,
        frame=fit, residual=resid,
        n_matched_for_frame=int(len(j)),
        mean_pads_per_cluster=float(j['n_clus'].mean()),
        efficiency=dict(k=k, n=n, value=k / max(n, 1), lo=lo, hi=hi,
                        n_miss_no_hit=n_empty, n_miss_far=n_far,
                        vs_probe_r=[{'r_mm': r, 'eff': e} for r, e in scan],
                        accidental_per_10mm=float(slope * 10)),
        settle_t_min_s=float(t_min),
        spark_intervals=int(len(veto.intervals)) if veto is not None else 0)
    return summary, j[core], hit


# --------------------------------------------------------------------------- #
# plots
# --------------------------------------------------------------------------- #
def plot_station(det, j, hit, summary, args, out_png):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    def eff_map(a, xcol, ycol, xlabel, ylabel, title):
        bs = args.bin
        xe = np.arange(hit[xcol].min(), hit[xcol].max() + bs, bs)
        ye = np.arange(hit[ycol].min(), hit[ycol].max() + bs, bs)
        tot, _, _ = np.histogram2d(hit[xcol], hit[ycol], bins=[xe, ye])
        ok, _, _ = np.histogram2d(hit.loc[hit['found'], xcol],
                                  hit.loc[hit['found'], ycol], bins=[xe, ye])
        with np.errstate(invalid='ignore', divide='ignore'):
            e = np.where(tot >= args.min_bin, ok / tot, np.nan)
        im = a.pcolormesh(xe, ye, e.T, cmap='viridis', vmin=args.eff_vmin, vmax=1.0)
        fig.colorbar(im, ax=a, label='efficiency')
        a.set(xlabel=xlabel, ylabel=ylabel, title=title)
        a.set_aspect('equal')

    fig, ax = plt.subplots(2, 4, figsize=(23, 10))
    w = args.match_win
    b = np.linspace(-w, w, 100)
    a = ax[0, 0]
    for c, lbl in (('rx', 'dx'), ('ry', 'dy')):
        a.hist(j[c], bins=b, histtype='step', lw=1.6,
               label=f'{lbl}  rms {j[c].std():.2f}  '
                     f'sigma(IQR) {summary["residual"][lbl[1]]["sigma_iqr_mm"]:.2f} mm')
    a.set(xlabel='P2 cluster - uRWELL track [mm]', ylabel='events',
          title=f'{det.name} residuals')
    a.legend(fontsize=8)
    a.grid(alpha=0.3)

    a = ax[0, 1]
    hh = a.hist2d(j['rx'], j['ry'], bins=[b, b], cmin=1)
    fig.colorbar(hh[3], ax=a)
    a.set(xlabel='dx [mm]', ylabel='dy [mm]', title='2D residual')
    a.set_aspect('equal')

    for col, axis, lbl in ((2, 'ex', 'x'), (3, 'ey', 'y')):
        a = ax[0, col]
        for c, cl in (('rx', 'dx'), ('ry', 'dy')):
            g = j.groupby(pd.cut(j[axis], 40), observed=True)[c]
            m, e = g.median(), g.std() / np.sqrt(g.size().clip(lower=1))
            xs = [iv.mid for iv in m.index]
            a.errorbar(xs, m.values, yerr=e.values, fmt='o', ms=3, label=cl)
        a.axhline(0, color='k', lw=0.8)
        a.set(xlabel=f'track {lbl} in P2 frame [mm]',
              ylabel='median residual [mm]',
              title=f'residual vs {lbl} (leftover rotation/scale)')
        a.legend(fontsize=8)
        a.grid(alpha=0.3)

    a = ax[1, 0]
    hh = a.hist2d(hit['ux'], hit['uy'], bins=130, cmin=1)
    fig.colorbar(hh[3], ax=a)
    a.set(xlabel='track x in the uRWELL frame [mm]', ylabel='track y [mm]',
          title=f'track density, uRWELL frame ({len(hit)} in fiducial)')
    a.set_aspect('equal')

    eff_map(ax[1, 1], 'ux', 'uy', 'track x in the uRWELL frame [mm]',
            'track y [mm]', 'efficiency in the uRWELL frame\n'
            '(features here belong to the REFERENCE)')
    eff_map(ax[1, 2], 'ex', 'ey', 'track x in the P2 frame [mm]',
            'track y [mm]', f'efficiency in the P2 pad frame\n'
            f'({args.bin:g} mm bins, >={args.min_bin} tracks)')

    a = ax[1, 3]
    d = np.linspace(0, 40, 81)
    a.hist(hit['dmin'].replace(np.inf, np.nan).dropna(), bins=d, histtype='step',
           lw=1.6)
    a.axvline(args.probe_r, color='r', ls='--', label=f'probe r = {args.probe_r:g} mm')
    e = summary['efficiency']
    a.set(xlabel='distance track -> nearest P2 cluster [mm]', ylabel='tracks',
          yscale='log',
          title=f'matching distance\n{e["n_miss_no_hit"]} misses with no P2 hit, '
                f'{e["n_miss_far"]} out of range')
    a.legend(fontsize=8)
    a.grid(alpha=0.3)

    fig.suptitle(
        f'{det.name}  z = {det.z:.0f} mm   uRWELL-referenced   '
        f'efficiency {e["value"]:.4f} (+{e["hi"] - e["value"]:.4f}/'
        f'-{e["value"] - e["lo"]:.4f}) on {e["n"]} tracks   '
        f'residual rms {summary["residual"]["x"]["rms_mm"]:.2f}/'
        f'{summary["residual"]["y"]["rms_mm"]:.2f} mm', fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_png, dpi=95)
    plt.close(fig)
    print(f'    wrote {out_png}')


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', default='highstat_eff_1')
    ap.add_argument('--data-root', default=DATA_ROOT)
    ap.add_argument('--sub-run', action='append', default=[],
                    help='repeatable; default = every sub_run with combined hits')
    ap.add_argument('--max-chunks', type=int, default=0,
                    help='0 = all chunk files')
    ap.add_argument('--track-cut', type=float, default=3.0,
                    help='|back - front| cut per axis [mm]')
    ap.add_argument('--cluster-r', type=float, default=15.0,
                    help='P2 pad clustering radius [mm]')
    ap.add_argument('--probe-r', type=float, default=15.0,
                    help='track <-> P2 cluster matching radius [mm]')
    ap.add_argument('--fid-r', type=float, default=9.0,
                    help='fiducial: nearest pad centre within this [mm]')
    ap.add_argument('--match-win', type=float, default=25.0,
                    help='frame-fit / residual window [mm]')
    ap.add_argument('--bin', type=float, default=4.0,
                    help='efficiency map bin [mm]')
    ap.add_argument('--min-bin', type=int, default=25,
                    help='minimum tracks in an efficiency-map bin')
    ap.add_argument('--eff-vmin', type=float, default=0.5)
    ap.add_argument('--min-events', type=int, default=500,
                    help='skip a station with fewer matched tracks than this')
    ap.add_argument('--no-veto-sparks', action='store_true')
    ap.add_argument('--out', default='out_eff')
    args = ap.parse_args()

    run_dir = os.path.join(args.data_root, args.run)
    run_json = os.path.join(run_dir, 'run_config.json')
    cfg = sc.RunConfig(key='urw_p2', run=args.run, data_root=args.data_root)
    os.makedirs(args.out, exist_ok=True)

    sub_runs = args.sub_run or cfg.find_subruns()
    dets = cfg.mappable_detectors()
    print(f'{args.run}: {len(sub_runs)} sub_run(s), probes '
          f'{[d.name for d in dets]}')

    all_summaries = []
    for sub_run in sub_runs:
        sub_dir = os.path.join(run_dir, sub_run)
        print(f'\n=== {sub_run}')
        tracks, tinfo, dz = build_tracks(run_json, sub_dir,
                                         args.max_chunks or None, args.track_cut)
        print(f'  uRWELL: {tinfo["n_4coord"]} four-coordinate events, '
              f'{tinfo["n_good_track"]} pass |back-front| < {args.track_cut:g} mm '
              f'({tinfo["good_track_frac"]:.1%})')
        for ax in ('x', 'y'):
            fb = tinfo['front_back'][ax]
            print(f'    {ax}: back = {fb["slope"]:+.4f} * front {fb["offset_mm"]:+.2f} mm'
                  f'   core sigma {fb["sigma_iqr_mm"]:.2f} mm')

        hits_dir = cfg.combined_hits_dir(sub_run)
        files = sorted(glob.glob(os.path.join(hits_dir, '*.root')))
        if args.max_chunks:
            files = files[:args.max_chunks]
        hv_csv = cfg.hv_monitor_csv(sub_run)

        t_span = float(tracks['t_ns'].max()) / 1e9 if len(tracks) else 0.0
        for det in dets:
            # The HV-settle cut is derived by matching hv_monitor.csv timestamps
            # against a DAQ start parsed out of the chunk filename.  On
            # drift_mesh_scan_1 that lands ~9.6 h off and the cut would silently
            # discard every event, so refuse a cut that does not fit inside the
            # sub_run rather than return an empty efficiency.
            t_min = scl.settle_t_min(hv_csv, det.spark_channel, files[0])
            if t_min > 0 and t_min >= t_span:
                print(f'  [{det.name}: settle cut {t_min:.0f} s exceeds the '
                      f'{t_span:.0f} s of trigger timestamps in this sub_run - '
                      f'hv_monitor.csv does not line up, cut disabled]')
                t_min = 0.0
            veto = None
            if not args.no_veto_sparks and os.path.isfile(hv_csv) and det.spark_channel:
                class _Shim:
                    SPARK_CHANNEL = det.spark_channel
                    SPARK_IMON_THR = cfg.SPARK_IMON_THR
                    SPARK_GUARD_BEFORE = cfg.SPARK_GUARD_BEFORE
                    SPARK_GUARD_AFTER = cfg.SPARK_GUARD_AFTER
                    BURST_NPADS = cfg.BURST_NPADS
                veto = ps.SparkVeto.from_csv(hv_csv, _Shim)
            s, j, hit = run_station(cfg, det, sub_run, tracks, dz, files,
                                    args, veto, t_min)
            if s is None:
                continue
            # which electrode this run scans for THIS station (mesh, drift or
            # neither) - drift_mesh_scan_1 scans P2_MID and P2_OUT while P2_IN
            # sits fixed as the control, so it is per station, not per run
            axis, label = cfg.scan_axis(sub_runs, det)
            s.update(run=args.run, sub_run=sub_run, urwell=tinfo,
                     scan_axis=axis, scan_label=label,
                     scan_hv=cfg.subrun_scan_hv(sub_run, det, axis),
                     mesh_hv=cfg.subrun_mesh_hv(sub_run, det),
                     drift_hv=cfg.subrun_drift_hv(sub_run, det))
            all_summaries.append(s)
            plot_station(det, j, hit, s, args,
                         os.path.join(args.out, f'{det.name}_{sub_run}.png'))
        # flush after every sub_run: a 23-point HV scan is hours of work and a
        # crash at the last point must not cost the first twenty-two
        _write_tables(all_summaries, args)

    df = _write_tables(all_summaries, args)
    print(df.to_string(index=False))

    if df['sub_run'].nunique() > 1:
        plot_summary(df, args, os.path.join(args.out, f'summary_{args.run}.png'))


def _write_tables(all_summaries, args):
    out_json = os.path.join(args.out, f'urw_p2_efficiency_{args.run}.json')
    with open(out_json, 'w') as fh:
        json.dump(all_summaries, fh, indent=2)
    rows = [{'run': s['run'], 'sub_run': s['sub_run'],
             'station': s['station'], 'z_mm': s['z_mm'],
             'n': s['efficiency']['n'], 'eff': s['efficiency']['value'],
             'lo': s['efficiency']['lo'], 'hi': s['efficiency']['hi'],
             'res_x_mm': s['residual']['x']['rms_mm'],
             'res_y_mm': s['residual']['y']['rms_mm'],
             'rot_deg': s['frame']['affine_rotation_deg'],
             'det': s['frame']['det'],
             'scan_axis': s['scan_axis'], 'scan_label': s['scan_label'],
             'scan_hv': s['scan_hv'], 'mesh_hv': s['mesh_hv'],
             'drift_hv': s['drift_hv'],
             'accidental_per_10mm': s['efficiency']['accidental_per_10mm']}
            for s in all_summaries]
    df = pd.DataFrame(rows)
    out_csv = os.path.join(args.out, f'urw_p2_efficiency_{args.run}.csv')
    df.to_csv(out_csv, index=False)
    print(f'  [wrote {out_json} and .csv: {len(df)} rows]', flush=True)
    return df


def plot_summary(df, args, out_png):
    """Efficiency and residual width per station across the sub_runs."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    subs = sorted(df['sub_run'].unique())
    xi = {s: i for i, s in enumerate(subs)}
    # A station whose own electrode is being scanned gets its HV on the x axis;
    # one held fixed (the control in a scan run, or every station in a run at a
    # single working point) falls back to the sub_run index.
    scanned = df[df['scan_hv'].notna()]
    use_hv = bool(len(scanned)) and \
        scanned.groupby('station')['scan_hv'].nunique().max() > 1
    # ...but only if the scanned electrode actually identifies the working
    # point.  drift_mesh_scan_1 scans BOTH: its ten drift_* sub_runs all sit at
    # mesh 450 V, so plotting against the mesh alone stacks them in a vertical
    # spike.  When that happens fall back to the sub_run index here - the
    # dedicated efficiency_vs_* figures from plot_hv_curves.py are the ones that
    # separate the two electrodes properly.
    if use_hv:
        for _, g in scanned.groupby('station'):
            partner = 'drift_hv' if (g['scan_axis'] == 'mesh').any() else 'mesh_hv'
            if partner in g and g.groupby('scan_hv')[partner].nunique().max() > 1:
                use_hv = False
                break
    xlabel = (' / '.join(sorted(scanned['scan_label'].dropna().unique()))
              if use_hv else 'sub_run')

    fig, ax = plt.subplots(1, 4, figsize=(22, 5))
    for station, g in df.groupby('station'):
        g = g.sort_values('scan_hv' if use_hv and g['scan_hv'].notna().all()
                          else 'sub_run')
        if use_hv and g['scan_hv'].notna().all() and g['scan_hv'].nunique() > 1:
            x, tag = g['scan_hv'].to_numpy(float), station
        elif use_hv:
            x, tag = np.full(len(g), np.nan), f'{station} (fixed)'
        else:
            x, tag = np.array([xi[s] for s in g['sub_run']], float), station
        ax[0].errorbar(x, g['eff'],
                       yerr=[g['eff'] - g['lo'], g['hi'] - g['eff']],
                       fmt='o-', ms=5, capsize=3, label=tag)
        ax[1].plot(x, g['res_x_mm'], 'o-', ms=5, label=f'{station} dx')
        ax[1].plot(x, g['res_y_mm'], 's--', ms=5, label=f'{station} dy')
        ax[2].plot(x, g['rot_deg'], 'o-', ms=5, label=station)
        ax[3].plot(x, g['accidental_per_10mm'], 'o-', ms=5, label=station)
    for a, ylab, ttl in (
            (ax[0], 'efficiency', 'uRWELL-referenced efficiency'),
            (ax[1], 'residual rms [mm]', 'residual width (pad-limited)'),
            (ax[2], 'uRWELL -> P2 rotation [deg]', 'frame stability'),
            (ax[3], 'd(eff)/dr per 10 mm', 'accidental-match rate')):
        if not use_hv:
            a.set_xticks(range(len(subs)))
            a.set_xticklabels(subs, rotation=45, ha='right',
                              fontsize=7 if len(subs) > 10 else 8)
        a.set_xlabel(xlabel)
        a.set_ylabel(ylab)
        a.set_title(ttl)
        a.grid(alpha=0.3)
        a.legend(fontsize=7)
    ax[0].axhline(1.0, color='k', lw=0.8, ls=':')
    fig.suptitle(f'{df["run"].iloc[0]}: '
                 f'P2 stations referenced to uRWELL tracks '
                 f'(probe r = {args.probe_r:g} mm, error bars 68.27% Clopper-Pearson)')
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out_png, dpi=100)
    plt.close(fig)
    print(f'wrote {out_png}')


if __name__ == '__main__':
    main()
