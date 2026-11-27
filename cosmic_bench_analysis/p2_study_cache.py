#!/usr/bin/env python3
"""p2_study_cache.py -- one-time per-run cache for the 34-36 studies.

Builds (and reloads) everything the clustering / 3D-alignment / pillar studies
need, so each study does not re-stream the combined hits:

  hits.pkl   one row per mapped, kept P2 pad hit (spark-vetoed events and the
                 cfg noisy/hot pads removed, cfg time window applied):
                 eventId, channel_id, x_pad, y_pad, amp, time [ns], tot
  rays.pkl   every clean M3 single track (chi2 cut, time window, spark veto):
                 eventId, X_Up, Y_Up, Z_Up, X_Down, Y_Down, Z_Down
  geom.npz       live pad centres + tiles (pad frame), dead pad centres + tiles,
                 the stage-06 transform (R, s, pad_mean, m3_mean) at PLANE_Z.

Usage:  python3 p2_study_cache.py <run_key> [--rebuild]
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2_qa_config as qa      # noqa: E402
import p2_mapping as pmap      # noqa: E402
import p2_io as p2io           # noqa: E402
import p2_align as pa          # noqa: E402
import p2_sparks as ps         # noqa: E402


def channel_table(cfg):
    return pmap.build_channel_table(cfg.run_config_path, cfg.MAP_CSV_PATH,
                                    det_type=cfg.DET_TYPE, det_name=cfg.DET_NAME,
                                    strategy='reverse',
                                    drop_connectors=cfg.DEAD_CONNECTORS,
                                    drop_channels=cfg.DEAD_CHANNELS,
                                    cluster_link=cfg.CLUSTER_LINK_MM,
                                    dead_guard=cfg.DEAD_GUARD_MM,
                                    strategy_overrides=cfg.STRATEGY_OVERRIDES)


def build(run_key):
    cfg = qa.get_config(run_key)
    out = cfg.out_dir('34_study_cache')
    ct = channel_table(cfg)
    drop = set(p2io.drop_pads_for(cfg, ct))
    sv = ps.SparkVeto.from_cfg(cfg)
    bad = sv.vetoed_ids_from_hits(cfg.combined_hits_dir, ct.attrs['feus'], min_amp=cfg.MIN_AMP)

    parts = []
    for df in p2io.iter_hits(cfg.combined_hits_dir,
                             ['eventId', 'channel', 'amplitude', 'feu', 'time',
                              'time_over_threshold'],
                             ct.attrs['feus'], t_max_h=cfg.T_MAX_H, t_min_h=cfg.T_MIN_H,
                             min_amp=cfg.MIN_AMP):
        h = pmap.attach_pads_to_hits(df, ct)
        h = h[h['mapped'] & h['pad_cx'].notna() & ~h['channel_id'].isin(drop)
              & ~h['eventId'].isin(bad)]
        parts.append(pd.DataFrame({
            'eventId': h['eventId'].to_numpy(np.int64),
            'channel_id': h['channel_id'].to_numpy(),
            'x_pad': h['pad_cx'].to_numpy(float), 'y_pad': h['pad_cy'].to_numpy(float),
            'amp': h['amplitude'].clip(lower=0).to_numpy(float),
            'time': h['time'].to_numpy(float),
            'tot': h['time_over_threshold'].to_numpy(float)}))
    hits = pd.concat(parts, ignore_index=True)
    # one row per pad per event (duplicated channel rows -> keep the largest pulse)
    hits = hits.sort_values('amp', ascending=False).drop_duplicates(['eventId', 'channel_id'])
    hits.to_pickle(f'{out}/hits.pkl')

    ep = pa.load_m3_endpoints(cfg.m3_tracking_dir, qa.M3_CHI2_CUT, qa.M3_MIN_NCLUS)
    ep = pa.filter_events_by_time(ep, cfg.m3_tracking_dir, cfg.T_MAX_H, t_min_h=cfg.T_MIN_H)
    ep = ep[~ep['eventId'].isin(bad)].reset_index(drop=True)
    ep.to_pickle(f'{out}/rays.pkl')

    # stage-06 transform: rigid+scale fit of the clustered centroid at PLANE_Z
    z = cfg.det_plane_z()
    cen, _ = pa.load_p2_centroids(cfg.combined_hits_dir, ct, min_amp=cfg.MIN_AMP,
                                  t_max_h=cfg.T_MAX_H, t_min_h=cfg.T_MIN_H,
                                  drop_pads=tuple(drop))
    x, y = pa.project_to_z(ep, z)
    m = pd.DataFrame({'eventId': ep['eventId'], 'x_m3': x, 'y_m3': y}).merge(cen, on='eventId')
    m = m[(m['x_m3'].abs() < 300) & (m['y_m3'].abs() < 300)]
    T = pa.fit_transform(m['x_m3'], m['y_m3'], m['x_pad'], m['y_pad'])

    padc = ct[ct['mapped']].drop_duplicates('channel_id')
    tpads, verts = pmap.pad_tiles(ct)
    dxy = ct.attrs.get('dead_pad_xy'); dv = ct.attrs.get('dead_pad_verts')
    np.savez(f'{out}/geom.npz',
             pad_id=padc['channel_id'].to_numpy(), pad_xy=padc[['pad_cx', 'pad_cy']].to_numpy(float),
             tile_verts=np.asarray(verts, float),
             dead_xy=np.asarray(dxy if dxy is not None else np.empty((0, 2)), float),
             dead_verts=np.asarray(dv if dv is not None else np.empty((0, 4, 2)), float),
             R=T.R, s=T.s, pad_mean=T.pad_mean, m3_mean=T.m3_mean, plane_z=z,
             match_r=cfg.MATCH_R, active_r=cfg.ACTIVE_R, dead_guard=cfg.DEAD_GUARD_MM,
             gap=cfg.DRIFT_GAP_MM)
    print(f'cache: {len(hits):,} pad hits, {hits.eventId.nunique():,} events, {len(ep):,} rays; '
          f'T rot {T.rotation_deg:.2f} s {T.s:.4f} rmse {T.rmse:.2f} at z={z}')
    return out


def load(run_key, rebuild=False):
    cfg = qa.get_config(run_key)
    out = cfg.out_dir('34_study_cache')
    if rebuild or not os.path.isfile(f'{out}/geom.npz'):
        build(run_key)
    g = dict(np.load(f'{out}/geom.npz', allow_pickle=True))
    return (pd.read_pickle(f'{out}/hits.pkl'), pd.read_pickle(f'{out}/rays.pkl'), g, cfg)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('run_key', nargs='?', default='det2_long_repos')
    ap.add_argument('--rebuild', action='store_true')
    a = ap.parse_args()
    load(a.run_key, a.rebuild)
