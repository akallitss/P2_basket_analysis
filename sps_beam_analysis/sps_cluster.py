#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sps_cluster.py

Shared beam-telescope infrastructure for the SPS stages, kept in one place so
the alignment (21), tag-and-probe (22) and profile (23) stages reuse exactly
the same per-event cluster extraction and the same plane-to-plane rigid
transform (no copies).

Two pieces:

1. `stream_event_clusters` -- memory-safe (p2_io chunk streaming) reduction of a
   station's combined hits into ONE leading cluster per event: leading pad plus
   every mapped pad within `cluster_r` of it, charge-weighted centroid, and a
   `single` flag that is True when the event has no pad outside that radius
   (i.e. one clean cluster). This is the "clean single cluster" used to tag and
   to align.

2. `RigidTransform` / `fit_rigid` -- the beam is ~parallel to z with negligible
   multiple scattering, so the map between two P2 planes (same PCB geometry) is
   a rigid 2D transform: translation + rotation about z, NO scale. `fit_rigid`
   is the scale-free Kabsch/Procrustes solution; `RigidTransform` carries the
   (dx, dy, theta) that the alignment JSON stores and applies/inverts them so a
   tag position can be carried from any plane into any other.
"""

import os
import re

import numpy as np
import pandas as pd

import p2_mapping as pmap
import p2_io as p2io


# --------------------------------------------------------------------------- #
# HV-settle cut (mesh still ramping at the start of a point) -- adapted from
# cosmic_bench_analysis/18_fe55_spectra.settle_t_min, keyed on a station's mesh.
# --------------------------------------------------------------------------- #
def settle_t_min(hv_csv, mesh_channel, chunk0, margin=5.0):
    """Event time [s since DAQ start] before which the station mesh was still
    ramping to its set point. 0 when there is no ramp / no monitor file."""
    if not mesh_channel or not os.path.isfile(hv_csv):
        return 0.0
    df = pd.read_csv(hv_csv)
    vcol, scol = f'{mesh_channel} vmon', f'{mesh_channel} v0'
    if vcol not in df.columns or scol not in df.columns:
        return 0.0
    ts = pd.to_datetime(df['timestamp'])
    ok = (df[vcol].astype(float) - df[scol].astype(float)).abs() < 2.0
    if not ok.any():
        return 0.0
    i_ok = int(np.argmax(ok.to_numpy()))
    if (ts.iloc[i_ok] - ts.iloc[0]).total_seconds() < 10.0:
        return 0.0
    t_settle = ts.iloc[i_ok]
    m = re.search(r'(\d{6})_(\d{2})H(\d{2})', os.path.basename(chunk0))
    if not m:
        return 0.0
    daq_start = pd.Timestamp(f'20{m.group(1)[:2]}-{m.group(1)[2:4]}-'
                             f'{m.group(1)[4:]} {m.group(2)}:{m.group(3)}:00')
    return max(0.0, (t_settle - daq_start).total_seconds() + margin)


# --------------------------------------------------------------------------- #
# Per-event leading cluster (streamed)
# --------------------------------------------------------------------------- #
def stream_event_clusters(hits_dir, channel_table, cluster_r, min_amp=0.0,
                          veto=None, t_min=0.0, drop_pads=(), progress=False):
    """Per-event leading cluster for one station, streamed chunk by chunk.

    A cluster = the leading (max-amplitude) pad plus every mapped pad within
    `cluster_r` mm of it. Returns a DataFrame indexed 0..N-1 with columns:
        eventId, x, y      charge-weighted cluster centroid [mm]
        q                  cluster charge [ADC]
        n_clus             pads in the cluster
        n_pad              total mapped pads in the event
        a_lead             leading-pad amplitude [ADC]
        lead_pad           leading pad channel_id
        t                  leading-pad trigger time [s since run start]
        single             True when n_pad == n_clus (one clean cluster)

    veto: optional p2_sparks.SparkVeto (spark + burst events removed).
    t_min: HV-settle cut in seconds (drop earlier events).
    """
    branches = ['eventId', 'channel', 'amplitude', 'feu', 'trigger_timestamp_ns']
    drop = set(drop_pads)
    parts = []
    for df in p2io.iter_hits(hits_dir, branches, channel_table.attrs['feus'],
                             min_amp=min_amp, progress=progress):
        if t_min > 0:
            df = df[df['trigger_timestamp_ns'].astype(np.int64) / 1e9 >= t_min]
        if veto is not None and len(df):
            df, _ = veto.apply(df)
        if not len(df):
            continue
        h = pmap.attach_pads_to_hits(df, channel_table)
        h = h[h['mapped'] & h['pad_cx'].notna()]
        if drop:
            h = h[~h['channel_id'].isin(drop)]
        del df
        if not len(h):
            continue
        # leading pad per event
        lead = h.loc[h.groupby('eventId')['amplitude'].idxmax(),
                     ['eventId', 'amplitude', 'pad_cx', 'pad_cy', 'channel_id',
                      'trigger_timestamp_ns']].rename(
            columns={'amplitude': 'a_lead', 'pad_cx': 'lx', 'pad_cy': 'ly',
                     'channel_id': 'lead_pad'})
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
        ev = pd.DataFrame({'x': agg['_wx'] / agg['_w'],
                           'y': agg['_wy'] / agg['_w'],
                           'q': agg['_w'], 'n_clus': agg['n_clus']})
        ev = ev.join(npad).join(
            lead.set_index('eventId')[['a_lead', 'lead_pad',
                                       'trigger_timestamp_ns']])
        ev['t'] = ev['trigger_timestamp_ns'].astype(np.int64) / 1e9
        parts.append(ev.drop(columns='trigger_timestamp_ns').reset_index())

    if not parts:
        return pd.DataFrame(columns=['eventId', 'x', 'y', 'q', 'n_clus',
                                     'n_pad', 'a_lead', 'lead_pad', 't',
                                     'single'])
    ev = pd.concat(parts, ignore_index=True)
    # an event straddling a chunk boundary appears twice: keep the larger half
    ev = (ev.sort_values('n_pad').drop_duplicates('eventId', keep='last')
          .reset_index(drop=True))
    ev['single'] = (ev['n_pad'] == ev['n_clus'])
    ev['eventId'] = ev['eventId'].astype(np.int64)
    return ev


# --------------------------------------------------------------------------- #
# Channel-space reductions (no pad geometry available)
# --------------------------------------------------------------------------- #
# Used for large detectors whose (feu, channel) -> pad wiring / pad map is not
# yet known (e.g. the Nov-2025 5-FEU P2_SPS25). Everything is done in the
# electronics address space (feu, channel) -- no x/y, no physical clustering
# across FEUs -- so the products are occupancy and an amplitude/charge Landau,
# honestly labelled as channel-space.
def stream_event_charge_chanspace(hits_dir, feus, cluster_chan=3, min_amp=0.0,
                                  veto=None, t_min=0.0, progress=True):
    """Per-event leading charge cluster in CHANNEL space (no geometry).

    For each event: the leading (max-amplitude) hit sets the reference (its FEU
    and channel); the cluster charge sums every hit on the SAME FEU within
    +-cluster_chan channels of the leading channel. This is a genuine
    leading-pad + neighbours MIP estimator in electronics space (adjacency is
    channel index, not physical pad -- documented). Returns a DataFrame keyed
    0..N-1 with eventId, lfeu, lchan, a_lead, q, n_chan, n_hit, t.
    """
    branches = ['eventId', 'channel', 'amplitude', 'feu', 'trigger_timestamp_ns']
    parts = []
    for df in p2io.iter_hits(hits_dir, branches, feus, min_amp=min_amp,
                             progress=progress):
        if t_min > 0:
            df = df[df['trigger_timestamp_ns'].astype(np.int64) / 1e9 >= t_min]
        if veto is not None and len(df):
            df, _ = veto.apply(df)
        if not len(df):
            continue
        lead = df.loc[df.groupby('eventId')['amplitude'].idxmax(),
                      ['eventId', 'feu', 'channel', 'amplitude',
                       'trigger_timestamp_ns']].rename(
            columns={'feu': 'lfeu', 'channel': 'lchan', 'amplitude': 'a_lead'})
        n_hit = df.groupby('eventId').size().rename('n_hit')
        d = df.merge(lead[['eventId', 'lfeu', 'lchan']], on='eventId')
        near = (d['feu'] == d['lfeu']) & \
            ((d['channel'].astype(np.int32) - d['lchan'].astype(np.int32))
             .abs() <= cluster_chan)
        dc = d[near]
        agg = dc.groupby('eventId').agg(q=('amplitude', 'sum'),
                                        n_chan=('amplitude', 'size'))
        ev = agg.join(lead.set_index('eventId')[['lfeu', 'lchan', 'a_lead',
                                                 'trigger_timestamp_ns']])
        ev = ev.join(n_hit)
        ev['t'] = ev['trigger_timestamp_ns'].astype(np.int64) / 1e9
        parts.append(ev.drop(columns='trigger_timestamp_ns').reset_index())
    if not parts:
        return pd.DataFrame(columns=['eventId', 'lfeu', 'lchan', 'a_lead', 'q',
                                     'n_chan', 'n_hit', 't'])
    ev = pd.concat(parts, ignore_index=True)
    ev = (ev.sort_values('n_chan').drop_duplicates('eventId', keep='last')
          .reset_index(drop=True))
    ev['eventId'] = ev['eventId'].astype(np.int64)
    return ev


def stream_chan_occupancy(hits_dir, feus, min_amp=0.0, veto=None, t_min=0.0,
                          progress=True):
    """Streamed per-(feu, channel) hit counts + per-event (n_hit, t).
    Returns (counts, events): counts is a Series indexed by (feu, channel);
    events a DataFrame indexed by eventId with n_hit and t [s]."""
    counts = None
    ev_parts = []
    for df in p2io.iter_hits(hits_dir, ['eventId', 'channel', 'amplitude',
                                        'feu', 'trigger_timestamp_ns'], feus,
                             min_amp=min_amp, progress=progress):
        if t_min > 0:
            df = df[df['trigger_timestamp_ns'].astype(np.int64) / 1e9 >= t_min]
        if veto is not None and len(df):
            df, _ = veto.apply(df)
        if not len(df):
            continue
        c = df.groupby(['feu', 'channel']).size()
        counts = c if counts is None else counts.add(c, fill_value=0)
        g = df.groupby('eventId')
        ev_parts.append(pd.DataFrame(
            {'n_hit': g.size(),
             't': g['trigger_timestamp_ns'].first().astype(np.int64) / 1e9}))
    if counts is None:
        return pd.Series(dtype=np.int64), pd.DataFrame(columns=['n_hit', 't'])
    counts = counts.astype(np.int64)
    ev = pd.concat(ev_parts).groupby(level=0).agg(n_hit=('n_hit', 'sum'),
                                                  t=('t', 'min'))
    return counts, ev


# --------------------------------------------------------------------------- #
# Rigid (translation + rotation, no scale) transform between planes
# --------------------------------------------------------------------------- #
class RigidTransform:
    """2D map  X = R(theta) @ x + (dx, dy).  Stored as (dx, dy, theta_deg)."""

    def __init__(self, dx=0.0, dy=0.0, theta_deg=0.0):
        self.dx = float(dx)
        self.dy = float(dy)
        self.theta_deg = float(theta_deg)

    @property
    def R(self):
        a = np.radians(self.theta_deg)
        c, s = np.cos(a), np.sin(a)
        return np.array([[c, -s], [s, c]])

    def apply(self, x, y):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        R = self.R
        return (R[0, 0] * x + R[0, 1] * y + self.dx,
                R[1, 0] * x + R[1, 1] * y + self.dy)

    def inverse(self):
        R = self.R.T
        t = -R @ np.array([self.dx, self.dy])
        return RigidTransform(t[0], t[1],
                              np.degrees(np.arctan2(R[1, 0], R[0, 0])))

    def to_dict(self):
        return {'dx': self.dx, 'dy': self.dy, 'theta_deg': self.theta_deg}

    @classmethod
    def from_dict(cls, d):
        return cls(d.get('dx', 0.0), d.get('dy', 0.0), d.get('theta_deg', 0.0))

    def __repr__(self):
        return (f'<Rigid dx={self.dx:.2f} dy={self.dy:.2f} '
                f'theta={self.theta_deg:.3f} deg>')


def fit_rigid(x_src, y_src, x_dst, y_dst):
    """Best-fit rigid transform (rotation+translation, NO scale) src -> dst,
    via scale-free Kabsch. Returns (RigidTransform, rmse [mm])."""
    src = np.column_stack([np.asarray(x_src, float), np.asarray(y_src, float)])
    dst = np.column_stack([np.asarray(x_dst, float), np.asarray(y_dst, float)])
    mu_s, mu_d = src.mean(0), dst.mean(0)
    H = (src - mu_s).T @ (dst - mu_d)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, d]) @ U.T
    t = mu_d - R @ mu_s
    tf = RigidTransform(t[0], t[1], np.degrees(np.arctan2(R[1, 0], R[0, 0])))
    fx, fy = tf.apply(src[:, 0], src[:, 1])
    rmse = float(np.sqrt(np.mean((fx - dst[:, 0]) ** 2 +
                                 (fy - dst[:, 1]) ** 2)))
    return tf, rmse
