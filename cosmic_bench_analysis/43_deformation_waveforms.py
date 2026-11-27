#!/usr/bin/env python3
"""43_deformation_waveforms.py -- the waveform dataset of the deformation study.

One pass over the decoded (non-zero-suppressed) DREAM files of a run. Keeps the
32-sample waveform (pedestal- and common-noise-subtracted, 60 ns samples) of

  * every mapped pad with a pulse (peak > 5 sigma and >= 3 samples > 3 sigma,
    the same rule as p2_waveforms.extract_hits), in EVERY triggered event, and
  * every mapped pad whose centre lies within --near mm of an M3 track (stage-34
    trimmed alignment), whatever its amplitude,

so that all later deformation stages (pulse shape per pad, drift-time spread,
pad-pair micro-TPC, sub-threshold sharing) run in seconds on a ~100 MB table
instead of re-reading 6.6 GB.

Products (<run>/43_deformation/):
  wf_table.pkl   one row per kept (event, pad): eventId, feu, channel,
                 channel_id, connector_N, half, pad_cx, pad_cy, rms, ftst,
                 timestamp, t_h (hours since run start; NaN if the event has no M3 record), is_hit, near_track,
                 saturated, vetoed (spark / burst event, stage-34 list), peak,
                 and for tracked events xp, yp (track at the plane, PAD frame),
                 theta, tx_p, ty_p (track slopes rotated into the pad frame), sd
                 (signed distance to the live-area edge), d_trk (track - pad centre)
  wf.npy         (n_rows, 32) float32, same order
  regions.csv    pad -> region partition (k-means, K = 8 / 24 / 48) used by every
                 deformation stage; written once, never overwritten
  summary.txt    counts and the assumptions (alignment, cuts)

Assumptions written to summary.txt: alignment = 34_clustering_matching/alignment.npz
(trimmed fit), plane z = cfg.det_plane_z(), no amplitude cut, no time cut applied
here (t_h and vetoed are stored; cut downstream).

Usage: python3 43_deformation_waveforms.py [run_key] [--near 22]
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2_study_cache as sc     # noqa: E402
import p2_study_common as cm    # noqa: E402
import p2_align as pa           # noqa: E402
import p2_sparks as ps          # noqa: E402
import p2_waveforms as pw       # noqa: E402

ADC_MAX = 4095.0


def inverse(T, xy):
    """M3 -> pad frame (inverse of p2_align.Transform.apply)."""
    return ((np.asarray(xy) - T.m3_mean) / T.s) @ T.R + T.pad_mean


def write_regions(ctm, out):
    """Pad -> region partition shared by the deformation stages: k-means on the
    live pad centres (pad frame), K = 8, 24, 48; labels ordered by centroid x."""
    from scipy.cluster.vq import kmeans2
    c = ctm.drop_duplicates('channel_id')
    xy = c[['pad_cx', 'pad_cy']].to_numpy(float)
    df = pd.DataFrame({'channel_id': c.channel_id.to_numpy(int), 'pad_cx': xy[:, 0], 'pad_cy': xy[:, 1],
                       'connector_N': c.connector_N.to_numpy(int), 'feu': c.feu.to_numpy(int),
                       'channel': c.channel.to_numpy(int)})
    for K in (8, 24, 48):
        cen, lab = kmeans2(xy, K, minit='++', seed=20261128, iter=50)
        order = np.lexsort((cen[:, 1], cen[:, 0])); rel = np.empty(K, int); rel[order] = np.arange(K)
        df[f'reg{K}'] = rel[lab]
    df.to_csv(f'{out}/regions.csv', index=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run_key', nargs='?', default='det2_long_repos')
    ap.add_argument('--near', type=float, default=22.0, help='keep pads within this of a track [mm]')
    args = ap.parse_args()
    hits, rays, g, cfg = sc.load(args.run_key)
    out = cfg.out_dir('43_deformation')
    d34 = cfg.out_dir('34_clustering_matching')
    A = np.load(f'{d34}/alignment.npz')
    T = cm.make_T(g, R=A['R'], s=A['s'], pad_mean=A['pad_mean'], m3_mean=A['m3_mean'])
    ra = pd.read_pickle(f'{d34}/rays_aligned.pkl').drop_duplicates('eventId').reset_index(drop=True)
    pxy = inverse(T, ra[['x', 'y']].to_numpy())
    # slopes in the pad frame: direction vectors transform with the rotation only
    tp = (ra[['tx', 'ty']].to_numpy() / T.s) @ T.R
    ra['xp'], ra['yp'], ra['tx_p'], ra['ty_p'] = pxy[:, 0], pxy[:, 1], tp[:, 0], tp[:, 1]
    trk = ra.set_index('eventId')[['xp', 'yp', 'theta', 'tx_p', 'ty_p', 'sd']]

    ct = sc.channel_table(cfg)
    ctm = ct[ct['mapped']].drop_duplicates(['feu', 'channel'])
    if not os.path.isfile(f'{out}/regions.csv'):
        write_regions(ctm, out)
    sv = ps.SparkVeto.from_cfg(cfg)
    bad = set(int(e) for e in sv.vetoed_ids_from_hits(cfg.combined_hits_dir, ct.attrs['feus'],
                                                       min_amp=cfg.MIN_AMP))
    et = pa.load_event_times(cfg.m3_tracking_dir)

    dec = pw.list_decoded_files(cfg.decoded_root_dir, feus=ct.attrs['feus'])
    tabs, waves = [], []
    n_ev = 0
    for _, fr in dec.iterrows():
        feu = int(fr.feu)
        c = ctm[ctm.feu == feu]
        chans = c.channel.to_numpy(int)
        cxy = c[['pad_cx', 'pad_cy']].to_numpy(float)
        tree = cKDTree(cxy)
        is_mapped = np.zeros(pw.N_CH, bool); is_mapped[chans] = True
        ped = pw.pedestal_from_data(fr.path)
        rms = pw.noise_rms_from_data(fr.path, ped)
        for evs, ts, ft, wf in pw.iter_decoded_events(fr.path, step=2000):
            if feu == dec.feu.min():
                n_ev += len(evs)
            w = pw.subtract_common_noise(wf - ped[None, :, None])
            peak = w.max(axis=2)
            nab = (w > (3.0 * rms)[None, :, None]).sum(axis=2)
            hit = (peak > 5.0 * rms[None, :]) & (nab >= 3) & is_mapped[None, :]
            keep = hit.copy()
            near = np.zeros_like(hit)
            has_t = np.isin(evs, trk.index.to_numpy())
            if has_t.any():
                k_idx = np.flatnonzero(has_t)
                txy = trk.loc[evs[k_idx], ['xp', 'yp']].to_numpy()
                for k, cand in zip(k_idx, tree.query_ball_point(txy, r=args.near)):
                    if len(cand):
                        near[k, chans[cand]] = True
                keep |= near
            iev, ich = np.nonzero(keep)
            if not len(iev):
                continue
            tabs.append(pd.DataFrame({
                'eventId': evs[iev].astype(np.int64), 'feu': feu, 'channel': ich,
                'rms': rms[ich], 'ftst': ft[iev].astype(np.int16), 'timestamp': ts[iev],
                'is_hit': hit[iev, ich], 'near_track': near[iev, ich],
                'saturated': wf[iev, ich, :].max(axis=1) >= ADC_MAX - 4,
                'peak': peak[iev, ich]}))
            waves.append(w[iev, ich, :].astype(np.float32))
        print(f'  {os.path.basename(fr.path)}: {sum(len(t) for t in tabs):,} rows so far', flush=True)

    R = pd.concat(tabs, ignore_index=True)
    WF = np.vstack(waves)
    R = R.merge(ctm[['feu', 'channel', 'channel_id', 'connector_N', 'half', 'pad_cx', 'pad_cy']],
                on=['feu', 'channel'], how='left')
    R = R.join(trk, on='eventId')
    R['d_trk'] = np.hypot(R.pad_cx - R.xp, R.pad_cy - R.yp)
    R['vetoed'] = R.eventId.isin(bad)
    R = R.merge(et[['eventId', 't_sec']], on='eventId', how='left')
    R['t_h'] = R['t_sec'] / 3600.0          # hours since run start (trigger clock)
    R = R.drop(columns='t_sec')
    R.to_pickle(f'{out}/wf_table.pkl')
    np.save(f'{out}/wf.npy', WF)
    L = [f'{args.run_key}: {len(dec)} decoded files, {n_ev:,} events on FEU {dec.feu.min()}',
         f'rows kept {len(R):,}: hit pads {int(R.is_hit.sum()):,}, track-near pads {int(R.near_track.sum()):,} '
         f'(within {args.near:g} mm of a track), both {int((R.is_hit & R.near_track).sum()):,}',
         f'events with a row {R.eventId.nunique():,}; tracked events {int(R.xp.notna().groupby(R.eventId).any().sum()):,}; '
         f'vetoed rows {int(R.vetoed.sum()):,}; saturated hit rows {int((R.saturated & R.is_hit).sum()):,}',
         f'assumptions: alignment = 34_clustering_matching/alignment.npz (trimmed fit, scale {float(T.s):.4f}), '
         f'plane z {cfg.det_plane_z():g} mm, hit rule peak > 5 sigma and >= 3 samples > 3 sigma, '
         f'pedestal = median of the first 400 events per file, CNS per 64-channel block; no amplitude or time cut '
         f'(cfg t_max_h = {cfg.T_MAX_H}: cut on t_h downstream)']
    with open(f'{out}/summary.txt', 'w') as f:
        f.write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
