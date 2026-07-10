#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p2_align.py

Shared P2<->M3 alignment helpers: load M3 single-track impacts and P2 pad
centroids, and fit the rigid (rotation + reflection + isotropic scale +
translation) transform that maps a P2 pad position into the M3 bench frame.

The transform is the one validated in 03_m3_alignment.py (reverse ordering,
~89 deg rotation, scale ~1.0). Both 03 and 06 use these helpers so the
efficiency stage reuses exactly the alignment logic, not a copy.
"""

import os
import sys
import glob

import numpy as np
import pandas as pd
import uproot
import awkward as ak

import p2_mapping as pmap
from p2_qa_config import M3_CHI2_CUT, M3_MIN_NCLUS

_M3_PKG = os.path.expanduser(
    '~/Documents/PostDocSaclay/nTof_x17/cosmic_bench_analysis')
if _M3_PKG not in sys.path:
    sys.path.insert(0, _M3_PKG)
from M3RefTracking import M3RefTracking  # noqa: E402

_HIT_BRANCHES = ['eventId', 'channel', 'amplitude', 'feu']


# --------------------------------------------------------------------------- #
def load_m3_positions(m3_dir, z, chi2_cut=M3_CHI2_CUT, single_track=True,
                      min_nclus=M3_MIN_NCLUS):
    """Single-track M3 impact positions at plane z, keyed by eventId.

    The good-track recipe is `Chi2X,Chi2Y < chi2_cut` AND (for tracking-v2 rays)
    `NClusX,NClusY >= min_nclus`; see p2_qa_config.M3_CHI2_CUT / M3_MIN_NCLUS."""
    m3 = M3RefTracking(os.path.join(m3_dir, ''), single_track=single_track,
                       chi2_cut=chi2_cut, min_nclus=min_nclus)
    x, y, evn = m3.get_xy_positions(z)
    df = pd.DataFrame({'eventId': np.asarray(evn, dtype=np.int64),
                       'x_m3': np.asarray(x, dtype=float),
                       'y_m3': np.asarray(y, dtype=float)})
    df = df[np.isfinite(df['x_m3']) & np.isfinite(df['y_m3'])]
    return df.drop_duplicates('eventId')


def load_event_times(m3_dir):
    """Per-event elapsed time [s] since run start, from the m3 tracking files.
    `evttime` is trigger_timestamp_ns/10, so *1e-8 gives seconds."""
    files = sorted(glob.glob(os.path.join(m3_dir, '*.root')))
    parts = []
    for fp in files:
        a = uproot.open(f'{fp}:T').arrays(['evn', 'evttime'], library='np')
        parts.append(pd.DataFrame({'eventId': a['evn'].astype(np.int64),
                                   't_sec': a['evttime'] * 1e-8}))
    return pd.concat(parts, ignore_index=True).drop_duplicates('eventId')


def load_m3_endpoints(m3_dir, chi2_cut=M3_CHI2_CUT, min_nclus=M3_MIN_NCLUS):
    """Single-track M3 line endpoints (X/Y at Z_Up and Z_Down) per eventId, so a
    track can be reprojected to ANY z cheaply for the z alignment scan."""
    m3 = M3RefTracking(os.path.join(m3_dir, ''), single_track=True,
                       chi2_cut=chi2_cut, min_nclus=min_nclus)
    rd = m3.ray_data

    def col(name):
        return np.asarray(ak.to_numpy(ak.ravel(rd[name])), dtype=float)

    df = pd.DataFrame({
        'eventId': np.asarray(ak.to_numpy(ak.ravel(rd['evn'])), dtype=np.int64),
        'X_Up': col('X_Up'), 'Y_Up': col('Y_Up'),
        'X_Down': col('X_Down'), 'Y_Down': col('Y_Down'),
        'Z_Up': col('Z_Up'), 'Z_Down': col('Z_Down'),
    })
    return df.drop_duplicates('eventId')


def project_to_z(ep, z):
    """Project M3 endpoint frame `ep` (from load_m3_endpoints) to plane z.
    Returns (x, y) numpy arrays."""
    f = (z - ep['Z_Down']) / (ep['Z_Up'] - ep['Z_Down'])
    x = ep['X_Down'] + (ep['X_Up'] - ep['X_Down']) * f
    y = ep['Y_Down'] + (ep['Y_Up'] - ep['Y_Down']) * f
    return np.asarray(x), np.asarray(y)


def load_p2_centroids(hits_dir, channel_table, min_amp=0.0, leading_pad=False):
    """Per-event P2 pad centroid (charge-weighted or leading-pad), keyed by eventId.
    Also returns the set of eventIds with any mapped P2 hit."""
    files = sorted(glob.glob(os.path.join(hits_dir, '*.root')))
    feu_set = set(channel_table.attrs['feus'])
    parts = []
    for fp in files:
        arr = uproot.open(f'{fp}:hits').arrays(_HIT_BRANCHES, library='pd')
        parts.append(arr[arr['feu'].isin(feu_set)].copy())
    hits = pd.concat(parts, ignore_index=True)
    if min_amp > 0:
        hits = hits[hits['amplitude'] >= min_amp]
    hits = pmap.attach_pads_to_hits(hits, channel_table)
    hits = hits[hits['mapped'] & hits['pad_cx'].notna()]
    hit_events = set(int(e) for e in hits['eventId'].unique())

    if leading_pad:
        idx = hits.groupby('eventId')['amplitude'].idxmax()
        cen = hits.loc[idx, ['eventId', 'pad_cx', 'pad_cy']].rename(
            columns={'pad_cx': 'x_pad', 'pad_cy': 'y_pad'})
        cen['n_pad'] = 1
    else:
        w = hits['amplitude'].clip(lower=0)
        hits = hits.assign(_wx=w * hits['pad_cx'], _wy=w * hits['pad_cy'], _w=w)
        g = hits.groupby('eventId')
        cen = pd.DataFrame({'x_pad': g['_wx'].sum() / g['_w'].sum(),
                            'y_pad': g['_wy'].sum() / g['_w'].sum(),
                            'n_pad': g.size()}).reset_index()
    cen = cen[np.isfinite(cen['x_pad']) & np.isfinite(cen['y_pad'])]
    return cen, hit_events


# --------------------------------------------------------------------------- #
class Transform:
    """Rigid+scale map  pad(x,y) -> M3(x,y):  m3 = m3_mean + s*(pad - pad_mean)@R^T."""

    def __init__(self, R, s, pad_mean, m3_mean):
        self.R = R
        self.s = s
        self.pad_mean = pad_mean
        self.m3_mean = m3_mean
        self.rotation_deg = np.degrees(np.arctan2(R[1, 0], R[0, 0])) % 360.0
        self.reflection = bool(np.linalg.det(R) < 0)

    def apply(self, x_pad, y_pad):
        P = np.column_stack([np.asarray(x_pad) - self.pad_mean[0],
                             np.asarray(y_pad) - self.pad_mean[1]])
        M = self.s * (P @ self.R.T) + self.m3_mean
        return M[:, 0], M[:, 1]


def fit_transform(x_m3, y_m3, x_pad, y_pad):
    """Best-fit rigid+reflection+scale transform pad->m3 (Procrustes/SVD)."""
    x_m3, y_m3 = np.asarray(x_m3), np.asarray(y_m3)
    x_pad, y_pad = np.asarray(x_pad), np.asarray(y_pad)
    pad_mean = np.array([x_pad.mean(), y_pad.mean()])
    m3_mean = np.array([x_m3.mean(), y_m3.mean()])
    P = np.column_stack([x_pad - pad_mean[0], y_pad - pad_mean[1]])
    Q = np.column_stack([x_m3 - m3_mean[0], y_m3 - m3_mean[1]])
    H = P.T @ Q
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    s = S.sum() / (P ** 2).sum()
    t = Transform(R, s, pad_mean, m3_mean)
    fit = s * (P @ R.T)
    t.rmse = float(np.sqrt(np.mean(np.sum((fit - Q) ** 2, axis=1))))
    return t
