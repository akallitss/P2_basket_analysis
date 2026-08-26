#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""p2_selftrack.py -- what the three P2 stations know about a track, judged
against the uRWELL reference.

`urw_p2_efficiency.py` asks one question of one station at a time: given a
reference track pointing at it, did it fire?  That measures *detection*.  In
the experiment there is no reference telescope: the three P2 stations are the
tracker, and what matters is where their own three-point track says the
particle went.  This measures *that*, with the uRWELL standing in for truth.

Everything about the reference, the alignment, the fiducial and the probe
radius is imported from `urw_p2_efficiency`, unchanged -- only the question is
new.  Three things come out:

  1. **Per-station pointing.**  The residual of each station's cluster to the
     reference track, split by cluster size, so the pitch/sqrt(12) floor and
     whatever charge sharing buys can be told apart.
  2. **Self-track pointing.**  A straight line through the three P2 clusters,
     compared with the reference track: offset and angle, and the error of
     extrapolating that line to a plane downstream of the basket -- which is
     the number the experiment actually cares about.
  3. **What self-referencing costs.**  P2_MID's residual to the P2_IN-P2_OUT
     line next to its residual to the reference.  Stage 22's tag-and-probe has
     to use the former; the gap between them is the systematic that buys.

The alignment is the reference's: each station's uRWELL -> pad affine comes
from `fit_frame`, so the *mean* of every delta below is zero by construction.
The *widths* are the measurement.  An experiment would align against something
external too, so this is the honest configuration, not a cheat -- but it does
mean nothing here bounds a global scale or rotation error.

Usage (lxplus, LCG_110, see run_selftrack.sh)
    python3 p2_selftrack.py --run eff_nominal_1 --out <dir>
"""
import os
import sys
import json
import glob
import argparse

import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import urw_p2_efficiency as E     # noqa: E402  (brings sps_config paths with it)
import sps_config as sc           # noqa: E402
import sps_cluster as scl         # noqa: E402
import p2_sparks as ps            # noqa: E402

STATIONS = ('P2_IN', 'P2_MID', 'P2_OUT')

# Histogram binning.  Fixed here rather than derived, so the arrays from two
# runs can be added without re-deriving edges.
RES_EDGES = np.arange(-25.0, 25.0 + 1e-9, 0.25)      # station residual [mm]
DMIN_EDGES = np.arange(0.0, 40.0 + 1e-9, 0.25)       # track <-> cluster [mm]
DPOS_EDGES = np.arange(-25.0, 25.0 + 1e-9, 0.1)      # self-track offset [mm]
DANG_EDGES = np.arange(-20.0, 20.0 + 1e-9, 0.1)      # self-track angle [mrad]
MAP_BIN = 4.0                                        # pad-frame map bin [mm]
# The face of a single pad, folded: every track's offset from the centre of the
# pad it points at.  A pad detector's efficiency is not flat across its own
# face -- near an edge the charge splits and the leading pad may drop under
# threshold -- and that is invisible in any map binned in absolute mm.
PAD_EDGES = np.arange(-7.0, 7.0 + 1e-9, 0.4)

# Where the self-track gets extrapolated to.  z is the beamline coordinate of
# the handoff (uRWELL front = 0); the stations are at 320 / 630 / 940 and the
# back uRWELL at 1370, so 1370 is a measured cross-check and the rest is the
# lever arm an experiment downstream of the basket would be working with.
Z_EXTRAP = (940.0, 1140.0, 1370.0, 1740.0, 2340.0, 2940.0)


def robust_sigma(d):
    """IQR-based width.  The residual of a pad detector is a box, not a
    Gaussian, so the rms is the honest summary and this is the companion that
    is not moved by the 1-in-10^3 tail; both are reported."""
    if len(d) < 10:
        return float('nan')
    q = np.percentile(d, [25, 75])
    return float((q[1] - q[0]) / 1.349)


def summarize(d, name=''):
    d = np.asarray(d, float)
    d = d[np.isfinite(d)]
    if not len(d):
        return {'n': 0}
    return {'n': int(len(d)), 'mean_mm': float(d.mean()),
            'median_mm': float(np.median(d)), 'rms_mm': float(d.std()),
            'sigma_iqr_mm': robust_sigma(d),
            'p95_abs_mm': float(np.percentile(np.abs(d), 95))}


# --------------------------------------------------------------------------- #
# one station: clusters, alignment, and the per-track table
# --------------------------------------------------------------------------- #
def station_table(cfg, det, sub_run, files, tracks, dz, args):
    """Per reference track, what this station did about it.

    Returns (table, frame_fit, pad_extent).  The table is one row per
    reference track -- NOT per cluster -- so a station that recorded nothing
    still occupies its row, which is what makes the three joinable.
    """
    from scipy.spatial import cKDTree
    ct = cfg.channel_table(det)
    pads = ct[ct['mapped']].drop_duplicates('channel_id')
    tree = cKDTree(pads[['pad_cx', 'pad_cy']].to_numpy())

    p2 = E.p2_clusters(files, ct, args.cluster_r, min_amp=cfg.MIN_AMP)
    print(f'  {det.name}: {len(p2)} events with a cluster')
    sys.stdout.flush()
    if len(p2) < args.min_events:
        return None, None, None

    src = E.project(tracks, det.z, dz)          # reference, uRWELL frame, at z

    # -- alignment, exactly as the efficiency stage does it ------------------ #
    j = pd.DataFrame({'eventId': tracks['eventId'].to_numpy(),
                      'px': src[:, 0], 'py': src[:, 1]}).merge(
        p2[['eventId', 'x', 'y']], on='eventId', how='inner')
    if len(j) < args.min_events:
        return None, None, None
    fit, tf, _ = E.fit_frame(j[['px', 'py']].to_numpy(),
                             j[['x', 'y']].to_numpy(), args.match_win)
    print(f'    frame rot {fit["affine_rotation_deg"]:+.3f} deg  '
          f'det {fit["det"]:+.4f}')

    # -- the reference track, in this station's pad frame -------------------- #
    ex, ey = tf.apply(src[:, 0], src[:, 1])
    d_pad, i_pad = tree.query(np.column_stack([ex, ey]))
    pc = pads[['pad_cx', 'pad_cy']].to_numpy()

    t = pd.DataFrame({'eventId': tracks['eventId'].to_numpy(),
                      'ex': ex, 'ey': ey, 'fid': d_pad < args.fid_r,
                      # where on the face of its own pad the track landed
                      'qx': ex - pc[i_pad, 0], 'qy': ey - pc[i_pad, 1],
                      # and which pad that was -- the handle a dead-pad mask
                      # needs.  Unused by this stage's own products.
                      'pid': pads['channel_id'].to_numpy()[i_pad],
                      'pa': pads['pad_angle'].to_numpy()[i_pad]})
    t = t.merge(p2[['eventId', 'x', 'y', 'x2', 'y2', 'n_clus', 'n_pad',
                    'a_lead', 'q']], on='eventId', how='left')
    t['dx'] = t['x'] - t['ex']
    t['dy'] = t['y'] - t['ey']
    # ... and the same residual along the pad's OWN axes.  The pads are a fan,
    # so a residual taken along the board axes is the projection of a rotated
    # rectangle and its edge is smeared by however much the pad angle varies
    # across the beam spot -- a few tenths of a millimetre, which is the size
    # of the thing the edge is used to measure.
    ca = np.cos(np.deg2rad(t['pa'].to_numpy()))
    sa = np.sin(np.deg2rad(t['pa'].to_numpy()))
    t['rw'] = t['dx'] * ca + t['dy'] * sa
    t['rh'] = -t['dx'] * sa + t['dy'] * ca
    d1 = np.hypot(t['dx'], t['dy'])
    d2 = np.hypot(t['x2'] - t['ex'], t['y2'] - t['ey'])
    t['dmin'] = np.fmin(d1.fillna(np.inf), d2.fillna(np.inf))
    t['found'] = t['dmin'] < args.probe_r

    # -- the cluster, back in the uRWELL frame ------------------------------- #
    # This is what makes a P2-only track possible: three stations, one frame.
    Ai = np.linalg.inv(np.asarray(fit['affine_A'], float))
    tv = np.asarray(fit['affine_t'], float)
    p = np.column_stack([t['x'].to_numpy(), t['y'].to_numpy()]) - tv
    u = p @ Ai.T
    t['ux'], t['uy'] = u[:, 0], u[:, 1]

    ext = (float(pads['pad_cx'].min()), float(pads['pad_cx'].max()),
           float(pads['pad_cy'].min()), float(pads['pad_cy'].max()))
    return t, fit, ext


# --------------------------------------------------------------------------- #
# the P2-only track
# --------------------------------------------------------------------------- #
def line_through(z, u):
    """Least-squares line u(z) for a stack of tracks.

    `z` is (nstation,), `u` is (ntrack, nstation).  Returns slope and the
    intercept at z = 0.  With three equally-weighted points this is the same
    fit the experiment would do, and with two it is exact.
    """
    z = np.asarray(z, float)
    zbar = z.mean()
    dz = z - zbar
    denom = float((dz ** 2).sum())
    ubar = u.mean(axis=1)
    slope = (u * dz).sum(axis=1) / denom
    return slope, ubar - slope * zbar


def self_track_block(m, zs, args):
    """Compare the P2-only track with the reference, for the events where all
    three stations found the particle."""
    out = {}
    all3 = np.ones(len(m), bool)
    for s in STATIONS:
        all3 &= m[f'found_{s}'].to_numpy() & m[f'fid_{s}'].to_numpy()
    out['n_all3'] = int(all3.sum())
    if all3.sum() < args.min_events:
        return out, None

    d = m[all3]
    hists = {}
    for ax in ('x', 'y'):
        u = np.column_stack([d[f'u{ax}_{s}'].to_numpy() for s in STATIONS])
        sl, ic = line_through(zs, u)
        # the reference line, in the same frame
        sl_ref = d[f'slope_{ax}'].to_numpy()
        ic_ref = d[f'{ax}f'].to_numpy()          # uRWELL front is z = 0

        dang = (sl - sl_ref) * 1e3               # mrad
        out[f'angle_{ax}_mrad'] = summarize(dang)
        out[f'angle_{ax}_mrad']['sigma_iqr_mrad'] = \
            out[f'angle_{ax}_mrad'].pop('sigma_iqr_mm')
        hists[f'dang_{ax}'] = np.histogram(dang, DANG_EDGES)[0]

        out[f'position_{ax}'] = {}
        for z in Z_EXTRAP:
            dpos = (ic + sl * z) - (ic_ref + sl_ref * z)
            out[f'position_{ax}'][f'z{z:.0f}'] = summarize(dpos)
            if abs(z - 940.0) < 1e-6:
                hists[f'dpos_{ax}'] = np.histogram(dpos, DPOS_EDGES)[0]

        # -- what self-referencing costs -------------------------------------
        # P2_MID predicted from the P2_IN <-> P2_OUT line (stage 22's geometry)
        # against P2_MID predicted from the reference.  Same events, same
        # clusters; only the ruler changes.
        z_io = np.array([zs[0], zs[2]])
        u_io = np.column_stack([d[f'u{ax}_P2_IN'].to_numpy(),
                                d[f'u{ax}_P2_OUT'].to_numpy()])
        sl_io, ic_io = line_through(z_io, u_io)
        pred_self = ic_io + sl_io * zs[1]
        pred_ref = ic_ref + sl_ref * zs[1]
        seen = d[f'u{ax}_P2_MID'].to_numpy()
        out[f'mid_vs_self_{ax}'] = summarize(seen - pred_self)
        out[f'mid_vs_ref_{ax}'] = summarize(seen - pred_ref)
        hists[f'mid_self_{ax}'] = np.histogram(seen - pred_self, RES_EDGES)[0]
        hists[f'mid_ref_{ax}'] = np.histogram(seen - pred_ref, RES_EDGES)[0]

    return out, hists


# --------------------------------------------------------------------------- #
def build_joined(cfg, dets, sub_run, args, run_dir, run_json):
    """Everything up to the joined per-track table, shared with p2_pillars.py.

    Split out of `run_subrun` unchanged so a second stage can ask a different
    question of exactly the same tracks, alignment and vetoes -- the moment two
    stages build their own denominators they stop being comparable.
    """
    sub_dir = os.path.join(run_dir, sub_run)
    print(f'\n=== {sub_run}')
    sys.stdout.flush()
    tracks, tinfo, dz = E.build_tracks(run_json, sub_dir,
                                       args.max_chunks or None, args.track_cut)
    print(f'  uRWELL: {tinfo["n_good_track"]} good tracks '
          f'({tinfo["good_track_frac"]:.1%} of 4-coordinate events)')
    if len(tracks) < args.min_events:
        print('  too few tracks, skipping')
        return None

    files = sorted(glob.glob(os.path.join(cfg.combined_hits_dir(sub_run),
                                          '*.root')))
    if args.max_chunks:
        files = files[:args.max_chunks]
    if not files:
        print('  no combined hits, skipping')
        return None
    hv_csv = cfg.hv_monitor_csv(sub_run)
    t_span = float(tracks['t_ns'].max()) / 1e9 if len(tracks) else 0.0

    # -- denominator hygiene, per station, then intersected ------------------ #
    # A three-station track needs every station to have been live for the
    # trigger, so the masks are ANDed rather than applied per station.  Doing
    # it any other way would let a station that was vetoed out contribute a
    # cluster to somebody else's track.
    keep = np.ones(len(tracks), bool)
    ev = tracks['eventId'].to_numpy()
    tns = tracks['t_ns'].to_numpy()
    veto_n = 0
    for det in dets:
        rec = (cfg.recorded_events(sub_run) or {}).get(det.feus[0])
        if rec is not None:
            keep &= E.recorded_mask(ev, rec)
        t_min = scl.settle_t_min(hv_csv, det.spark_channel, files[0])
        if t_min > 0 and t_min < t_span:
            keep &= (tns / 1e9) >= t_min
        if not args.no_veto_sparks and os.path.isfile(hv_csv) \
                and det.spark_channel:
            class _Shim:
                SPARK_CHANNEL = det.spark_channel
                SPARK_IMON_THR = cfg.SPARK_IMON_THR
                SPARK_GUARD_BEFORE = cfg.SPARK_GUARD_BEFORE
                SPARK_GUARD_AFTER = cfg.SPARK_GUARD_AFTER
                BURST_NPADS = cfg.BURST_NPADS
            v = ps.SparkVeto.from_csv(hv_csv, _Shim)
            if v is not None and v.intervals:
                keep &= v.event_mask(tns)
                veto_n += len(v.intervals)
    tracks = tracks[keep].reset_index(drop=True)
    print(f'  {len(tracks)} tracks survive recorded/settle/spark '
          f'({veto_n} spark intervals over the three stations)')
    sys.stdout.flush()

    # -- per station --------------------------------------------------------- #
    tabs, fits, exts, zs = {}, {}, {}, []
    for det in dets:
        t, fit, ext = station_table(cfg, det, sub_run, files, tracks, dz, args)
        if t is None:
            print(f'  {det.name}: not enough matched tracks, skipping sub_run')
            return None
        tabs[det.name], fits[det.name], exts[det.name] = t, fit, ext
        zs.append(det.z)
    zs = np.array(zs, float)

    # -- join, keeping one row per reference track --------------------------- #
    m = pd.DataFrame({'eventId': tracks['eventId'].to_numpy(),
                      'xf': tracks['xf'].to_numpy(),
                      'yf': tracks['yf'].to_numpy(),
                      'slope_x': tracks['slope_x'].to_numpy(),
                      'slope_y': tracks['slope_y'].to_numpy()})
    for s in STATIONS:
        cols = ['ex', 'ey', 'qx', 'qy', 'dx', 'dy', 'dmin', 'found', 'fid',
                'ux', 'uy', 'n_clus', 'a_lead', 'pid', 'rw', 'rh']
        m = m.join(tabs[s][cols].add_suffix(f'_{s}'))
    return m, tracks, tinfo, zs, fits, exts


def run_subrun(cfg, dets, sub_run, args, run_dir, run_json):
    got = build_joined(cfg, dets, sub_run, args, run_dir, run_json)
    if got is None:
        return None
    m, tracks, tinfo, zs, fits, exts = got

    res = {'run': args.run, 'sub_run': sub_run, 'z_mm': zs.tolist(),
           'n_tracks': int(len(tracks)), 'urwell': tinfo,
           'probe_r_mm': args.probe_r, 'fid_r_mm': args.fid_r,
           'frame': {s: fits[s] for s in STATIONS},
           'stations': {}}
    hists = {}

    # -- 1. per-station pointing --------------------------------------------- #
    for s in STATIONS:
        f = m[f'fid_{s}'].to_numpy() & m[f'found_{s}'].to_numpy()
        blk = {'n_fid': int(m[f'fid_{s}'].sum()),
               'n_found': int(f.sum()),
               'efficiency': float(f.sum() / max(int(m[f'fid_{s}'].sum()), 1)),
               'mean_pads_per_cluster': float(
                   m.loc[f, f'n_clus_{s}'].mean())}
        nc = m.loc[f, f'n_clus_{s}'].to_numpy()
        for ax in ('x', 'y'):
            d = m.loc[f, f'd{ax}_{s}'].to_numpy()
            core = np.abs(d) < args.match_win
            blk[f'resid_{ax}'] = summarize(d[core])
            # split by cluster size: a one-pad cluster can only say "this pad",
            # so it is the pitch/sqrt(12) box; anything better has to come from
            # a centroid over two or more.
            for tag, sel in (('single', nc <= 1), ('multi', nc >= 2)):
                blk[f'resid_{ax}_{tag}'] = summarize(d[core & sel])
                hists[f'res_{ax}_{s}_{tag}'] = np.histogram(
                    d[core & sel], RES_EDGES)[0]
        blk['frac_multi_pad'] = float((nc >= 2).mean()) if f.sum() else 0.0
        hists[f'dmin_{s}'] = np.histogram(
            m.loc[m[f'fid_{s}'], f'dmin_{s}'].to_numpy(), DMIN_EDGES)[0]
        res['stations'][s] = blk
        print(f'  {s}: eff {blk["efficiency"]:.4f}  '
              f'resid rms {blk["resid_x"]["rms_mm"]:.2f}/'
              f'{blk["resid_y"]["rms_mm"]:.2f} mm  '
              f'multi-pad {blk["frac_multi_pad"]:.1%}')

    # -- 2. how many stations fired at all ----------------------------------- #
    fid_all = np.ones(len(m), bool)
    nfound = np.zeros(len(m), int)
    for s in STATIONS:
        fid_all &= m[f'fid_{s}'].to_numpy()
        nfound += (m[f'found_{s}'].to_numpy() & m[f'fid_{s}'].to_numpy())
    res['multiplicity'] = {
        'n_fid_all3': int(fid_all.sum()),
        'n_station_hist': [int(((nfound == k) & fid_all).sum())
                           for k in range(4)],
        'frac_3of3': float(((nfound == 3) & fid_all).sum()
                           / max(int(fid_all.sum()), 1)),
        'frac_ge2of3': float(((nfound >= 2) & fid_all).sum()
                             / max(int(fid_all.sum()), 1))}
    print(f'  3-station tracks: {res["multiplicity"]["frac_3of3"]:.4f} '
          f'of the {res["multiplicity"]["n_fid_all3"]} tracks inside all three')

    # -- 3. the self-track ---------------------------------------------------- #
    st, sthist = self_track_block(m, zs, args)
    res['selftrack'] = st
    if sthist:
        hists.update({f'st_{k}': v for k, v in sthist.items()})
        print(f'  self-track vs reference at z=940: '
              f'dx rms {st["position_x"]["z940"]["rms_mm"]:.2f} mm, '
              f'angle sigma {st["angle_x_mrad"]["sigma_iqr_mrad"]:.2f} mrad')

    # -- 4. pad-frame maps ---------------------------------------------------- #
    # Counts, not ratios, so sub-runs can be added.  Efficiency and mean
    # residual are formed offline from these.
    maps = {}
    for s in STATIONS:
        x0, x1, y0, y1 = exts[s]
        xe = np.arange(x0 - MAP_BIN, x1 + 2 * MAP_BIN, MAP_BIN)
        ye = np.arange(y0 - MAP_BIN, y1 + 2 * MAP_BIN, MAP_BIN)
        fid = m[f'fid_{s}'].to_numpy()
        ex, ey = m[f'ex_{s}'].to_numpy(), m[f'ey_{s}'].to_numpy()
        fnd = fid & m[f'found_{s}'].to_numpy()
        maps[f'map_edges_x_{s}'] = xe
        maps[f'map_edges_y_{s}'] = ye
        maps[f'map_n_{s}'] = np.histogram2d(ex[fid], ey[fid], [xe, ye])[0]
        maps[f'map_k_{s}'] = np.histogram2d(ex[fnd], ey[fnd], [xe, ye])[0]
        r = np.hypot(m[f'dx_{s}'].to_numpy(), m[f'dy_{s}'].to_numpy())
        w = fnd & (r < args.match_win)
        maps[f'map_rsum_{s}'] = np.histogram2d(ex[w], ey[w], [xe, ye],
                                               weights=r[w])[0]
        maps[f'map_rn_{s}'] = np.histogram2d(ex[w], ey[w], [xe, ye])[0]

        # -- the same thing folded onto one pad ------------------------------
        qx, qy = m[f'qx_{s}'].to_numpy(), m[f'qy_{s}'].to_numpy()
        pe = [PAD_EDGES, PAD_EDGES]
        maps[f'padmap_n_{s}'] = np.histogram2d(qx[fid], qy[fid], pe)[0]
        maps[f'padmap_k_{s}'] = np.histogram2d(qx[fnd], qy[fnd], pe)[0]
        a = m[f'a_lead_{s}'].to_numpy()
        wa = fnd & np.isfinite(a)
        maps[f'padmap_asum_{s}'] = np.histogram2d(qx[wa], qy[wa], pe,
                                                  weights=a[wa])[0]
        maps[f'padmap_an_{s}'] = np.histogram2d(qx[wa], qy[wa], pe)[0]
        nc = m[f'n_clus_{s}'].to_numpy()
        wc = fnd & np.isfinite(nc)
        maps[f'padmap_csum_{s}'] = np.histogram2d(qx[wc], qy[wc], pe,
                                                  weights=nc[wc])[0]
        # 1D: efficiency against distance from the pad centre, which is the
        # projection that actually gets quoted
        rq = np.hypot(qx, qy)
        maps[f'padr_n_{s}'] = np.histogram(rq[fid], PAD_EDGES[PAD_EDGES >= 0])[0]
        maps[f'padr_k_{s}'] = np.histogram(rq[fnd], PAD_EDGES[PAD_EDGES >= 0])[0]
    hists.update(maps)
    hists['pad_edges'] = PAD_EDGES

    # -- 5. a sample of the joined table, for anything not foreseen ---------- #
    keepcols = ['xf', 'yf', 'slope_x', 'slope_y'] + [
        f'{c}_{s}' for s in STATIONS
        for c in ('ex', 'ey', 'qx', 'qy', 'dx', 'dy', 'dmin', 'found', 'fid',
                  'ux', 'uy', 'n_clus', 'a_lead')]
    sub = m[fid_all]
    if args.sample and len(sub) > args.sample:
        sub = sub.iloc[::int(np.ceil(len(sub) / args.sample))]
    hists['sample_cols'] = np.array(keepcols)
    hists['sample'] = sub[keepcols].to_numpy(np.float32)

    return res, hists


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', default='eff_nominal_1')
    ap.add_argument('--data-root', default=E.DATA_ROOT)
    ap.add_argument('--sub-run', action='append', default=[])
    ap.add_argument('--max-chunks', type=int, default=0)
    ap.add_argument('--track-cut', type=float, default=3.0)
    ap.add_argument('--cluster-r', type=float, default=15.0)
    ap.add_argument('--probe-r', type=float, default=15.0)
    ap.add_argument('--fid-r', type=float, default=9.0)
    ap.add_argument('--match-win', type=float, default=25.0)
    ap.add_argument('--min-events', type=int, default=500)
    ap.add_argument('--no-veto-sparks', action='store_true')
    ap.add_argument('--sample', type=int, default=40000,
                    help='rows of the joined table to keep per sub_run')
    ap.add_argument('--out', default='out_selftrack')
    args = ap.parse_args()

    run_dir = os.path.join(args.data_root, args.run)
    run_json = os.path.join(run_dir, 'run_config.json')
    cfg = sc.RunConfig(key='urw_p2', run=args.run, data_root=args.data_root)
    os.makedirs(args.out, exist_ok=True)

    dets = [d for d in cfg.mappable_detectors() if d.name in STATIONS]
    dets.sort(key=lambda d: d.z)
    if len(dets) != 3:
        sys.exit(f'need all three stations, found {[d.name for d in dets]}')
    sub_runs = args.sub_run or cfg.find_subruns()
    print(f'{args.run}: {len(sub_runs)} sub_run(s), stations '
          f'{[(d.name, d.z) for d in dets]}')

    summaries, acc = [], {}
    for sub_run in sub_runs:
        try:
            got = run_subrun(cfg, dets, sub_run, args, run_dir, run_json)
        except Exception as exc:                      # noqa: BLE001
            print(f'  !! {sub_run} failed: {exc!r}')
            continue
        if got is None:
            continue
        res, hists = got
        summaries.append(res)
        for k, v in hists.items():
            if k == 'sample_cols':
                acc[k] = v
            elif k.startswith('map_edges') or k == 'pad_edges':
                acc[k] = v
            elif k == 'sample':
                acc.setdefault(k, []).append(v)
            else:
                acc[k] = acc[k] + v if k in acc else v
        # flush after every sub_run: this is hours of work
        _write(summaries, acc, args)
    _write(summaries, acc, args)
    print(f'\ndone: {len(summaries)} sub_run(s)')


def _write(summaries, acc, args):
    base = os.path.join(args.out, f'p2_selftrack_{args.run}')
    with open(base + '.json', 'w') as fh:
        json.dump(summaries, fh, indent=1)
    out = dict(acc)
    if isinstance(out.get('sample'), list):
        out['sample'] = np.concatenate(out['sample'])
    out['res_edges'] = RES_EDGES
    out['dmin_edges'] = DMIN_EDGES
    out['dpos_edges'] = DPOS_EDGES
    out['dang_edges'] = DANG_EDGES
    np.savez_compressed(base + '.npz', **out)
    print(f'  wrote {base}.json / .npz')


if __name__ == '__main__':
    main()
