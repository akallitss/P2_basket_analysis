#!/usr/bin/env python3
"""Exploratory: how does the uRWELL local frame map onto the P2 pad frame?

Reads ONE chunk of highstat_eff_1/beam_commissioning_00, builds uRWELL
two-point tracks (FEU 1) and P2 leading clusters (FEU 3/4/5), matches on
eventId and fits a general 2D affine

    (X, Y)_P2 = A @ (x, y)_track + t

by robust least squares.  A general affine (not a rigid transform) is used on
purpose: the uRWELL axis parity is not known a priori - the connector-order
work of 2026-07-25 fixed the ordering only up to an overall axis reflection -
so the fit has to be free to come out with det(A) < 0 or with the axes
interchanged.  Whether A is close to orthogonal is then the check that the
mapping is physical.
"""
import os
import sys

import numpy as np
import pandas as pd

# this package, then sps_beam_analysis for sps_config, whose setup_paths()
# adds the shared core in cosmic_bench_analysis (p2_io, p2_mapping, ...)
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.dirname(_HERE))
import sps_config as _sc  # noqa: E402
_sc.setup_paths()

import urw_lib as U           # noqa: E402
import sps_config as sc       # noqa: E402
import sps_cluster as scl     # noqa: E402
import p2_mapping as pmap     # noqa: E402
import uproot                 # noqa: E402

RUN_DIR = '/local/home/banco/P2_data/TB_July2026_H4/runs/highstat_eff_1'
SUB = 'beam_commissioning_00'
N_CHUNK = 1


def urw_tracks(run_json, sub_dir, n_chunk):
    """eventId, x/y at both uRWELLs plus a slope, from the FEU-1 hits files."""
    pts = {}
    for det in ('EIC_uRWELL_front', 'EIC_uRWELL_back'):
        geo = U.UrwGeometry(det, run_json, sub_run_name=os.path.basename(sub_dir))
        files = U.feu_hit_files(sub_dir, feu=geo.feu_num)[:n_chunk]
        parts = []
        for chunk in U.iter_hits(files, progress=False, feu=geo.feu_num):
            cl = U.cluster_hits(chunk, geo)
            if len(cl):
                parts.append(U.leading_points(cl))
        pts[det] = (geo, pd.concat(parts, ignore_index=True))
        print(f'  {det}: {len(pts[det][1])} events with a leading cluster')
    (gf, pf), (gb, pb) = pts['EIC_uRWELL_front'], pts['EIC_uRWELL_back']
    dz = float(gb.center[2] - gf.center[2])
    ev = pf[['eventId', 'x', 'y', 'n_x', 'n_y']].rename(
        columns={'x': 'xf', 'y': 'yf', 'n_x': 'nxf', 'n_y': 'nyf'}).merge(
        pb[['eventId', 'x', 'y', 'n_x', 'n_y']].rename(
            columns={'x': 'xb', 'y': 'yb', 'n_x': 'nxb', 'n_y': 'nyb'}),
        on='eventId')
    ev = ev.dropna(subset=['xf', 'yf', 'xb', 'yb']).reset_index(drop=True)
    return ev, dz


def robust_affine(src, dst, nsig=3.0, iters=6):
    """dst = A @ src + t by least squares with nsig-sigma clipping.

    src, dst: (N, 2).  Returns (A, t, resid_rms, keep).
    """
    keep = np.ones(len(src), bool)
    A = np.eye(2)
    t = np.zeros(2)
    for _ in range(iters):
        M = np.column_stack([src[keep], np.ones(keep.sum())])
        sol, *_ = np.linalg.lstsq(M, dst[keep], rcond=None)
        A, t = sol[:2].T, sol[2]
        r = dst - (src @ A.T + t)
        d = np.hypot(r[:, 0], r[:, 1])
        new = d < nsig * d[keep].std()
        if new.sum() < 50 or new.sum() == keep.sum():
            keep = new
            break
        keep = new
    r = dst - (src @ A.T + t)
    return A, t, r, keep


def main():
    run_json = os.path.join(RUN_DIR, 'run_config.json')
    sub_dir = os.path.join(RUN_DIR, SUB)
    print(f'uRWELL tracks from {N_CHUNK} chunk(s) of {SUB}')
    ev, dz = urw_tracks(run_json, sub_dir, N_CHUNK)
    print(f'  matched front+back with all 4 coords: {len(ev)}   dz = {dz:.0f} mm')

    # front->back scale+offset, so the back can be carried into the front frame
    for ax in ('x', 'y'):
        s, o = np.polyfit(ev[f'{ax}f'], ev[f'{ax}b'], 1)
        for _ in range(4):
            r = ev[f'{ax}b'] - (s * ev[f'{ax}f'] + o)
            m = np.abs(r) < 2.5 * r.std()
            s, o = np.polyfit(ev.loc[m, f'{ax}f'], ev.loc[m, f'{ax}b'], 1)
        ev[f'{ax}b_c'] = (ev[f'{ax}b'] - o) / s
        ev[f'slope_{ax}'] = (ev[f'{ax}b_c'] - ev[f'{ax}f']) / dz
        print(f'  {ax}: back = {s:+.4f} * front {o:+.3f} mm')

    os.environ.setdefault('SPS_RUN', 'highstat_eff_1')
    cfg = sc.RunConfig(key='x', run='highstat_eff_1',
                       data_root='/local/home/banco/P2_data/TB_July2026_H4/runs')
    hits_dir = cfg.combined_hits_dir(SUB)
    files = sorted(os.listdir(hits_dir))[:N_CHUNK]
    print(f'\nP2 clusters from {files}')

    for det in cfg.mappable_detectors():
        ct = cfg.channel_table(det)
        # stream only the first chunk: point stream_event_clusters at a dir
        # listing we control by monkeypatching hit_files is overkill, so read
        # the chunk directly here.
        import p2_io as p2io
        allf = p2io.hit_files(hits_dir)
        p2 = _clusters_from_files(allf[:N_CHUNK], ct, cfg.CLUSTER_R)
        z = det.z
        m = ev.copy()
        m['px'] = m['xf'] + m['slope_x'] * z
        m['py'] = m['yf'] + m['slope_y'] * z
        j = m.merge(p2[['eventId', 'x', 'y', 'single', 'n_pad', 'q']],
                    on='eventId', how='inner')
        j = j[j['single']]
        print(f'\n=== {det.name}  z={z:.0f}  matched single clusters: {len(j)}')
        if len(j) < 500:
            print('   too few, skipping')
            continue
        src = j[['px', 'py']].to_numpy()
        dst = j[['x', 'y']].to_numpy()
        A, t, r, keep = robust_affine(src, dst)
        u, s, vt = np.linalg.svd(A)
        print(f'   A = [[{A[0,0]:+.4f} {A[0,1]:+.4f}], [{A[1,0]:+.4f} {A[1,1]:+.4f}]]')
        print(f'   t = ({t[0]:+.2f}, {t[1]:+.2f}) mm')
        print(f'   singular values {s[0]:.4f} {s[1]:.4f}   det = {np.linalg.det(A):+.4f}')
        rot = u @ vt
        print(f'   nearest rotation angle = {np.degrees(np.arctan2(rot[1,0], rot[0,0])):+.2f} deg')
        print(f'   residual rms  x {r[keep,0].std():.2f}  y {r[keep,1].std():.2f} mm'
              f'   ({keep.sum()}/{len(j)} kept)')
        # raw correlations, to see the parity by eye
        for a, b in (('px', 'x'), ('px', 'y'), ('py', 'x'), ('py', 'y')):
            print(f'     corr({a},{b}) = {np.corrcoef(j[a], j[b])[0,1]:+.3f}')


def _clusters_from_files(files, ct, cluster_r):
    """sps_cluster.stream_event_clusters over an explicit file list."""
    import p2_io as p2io
    branches = ['eventId', 'channel', 'amplitude', 'feu', 'trigger_timestamp_ns']
    feus = set(ct.attrs['feus'])
    parts = []
    for fp in files:
        with uproot.open(fp) as f:
            arrs = f['hits'].arrays(branches, library='np')
        df = pd.DataFrame(arrs)
        df = df[df['feu'].isin(feus)]
        h = pmap.attach_pads_to_hits(df, ct)
        h = h[h['mapped'] & h['pad_cx'].notna()]
        del df, arrs
        lead = h.loc[h.groupby('eventId')['amplitude'].idxmax(),
                     ['eventId', 'amplitude', 'pad_cx', 'pad_cy', 'channel_id']
                     ].rename(columns={'amplitude': 'a_lead', 'pad_cx': 'lx',
                                       'pad_cy': 'ly', 'channel_id': 'lead_pad'})
        npad = h.groupby('eventId').size().rename('n_pad')
        h = h.merge(lead[['eventId', 'lx', 'ly']], on='eventId')
        near = ((h['pad_cx'] - h['lx']) ** 2 +
                (h['pad_cy'] - h['ly']) ** 2) <= cluster_r ** 2
        hc = h[near]
        w = hc['amplitude'].clip(lower=0).astype(np.float64)
        g = pd.DataFrame({'eventId': hc['eventId'], '_wx': w * hc['pad_cx'],
                          '_wy': w * hc['pad_cy'], '_w': w})
        agg = g.groupby('eventId').agg(_wx=('_wx', 'sum'), _wy=('_wy', 'sum'),
                                       _w=('_w', 'sum'), n_clus=('_w', 'size'))
        e = pd.DataFrame({'x': agg['_wx'] / agg['_w'], 'y': agg['_wy'] / agg['_w'],
                          'q': agg['_w'], 'n_clus': agg['n_clus']})
        e = e.join(npad).join(lead.set_index('eventId')[['a_lead', 'lead_pad']])
        parts.append(e.reset_index())
    ev = pd.concat(parts, ignore_index=True)
    ev = (ev.sort_values('n_pad').drop_duplicates('eventId', keep='last')
          .reset_index(drop=True))
    ev['single'] = ev['n_pad'] == ev['n_clus']
    ev['eventId'] = ev['eventId'].astype(np.int64)
    return ev


if __name__ == '__main__':
    main()
