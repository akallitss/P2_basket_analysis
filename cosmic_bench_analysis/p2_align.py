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
import p2_io as p2io
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


def load_p2_centroids(hits_dir, channel_table, min_amp=0.0, leading_pad=False,
                      t_max_h=None, drop_pads=()):
    """Per-event P2 pad centroid (charge-weighted or leading-pad), keyed by eventId.
    Also returns the set of eventIds with any mapped P2 hit.

    Streams the combined-hits chunks (p2_io) so memory stays bounded no matter
    how many chunks the run has. min_amp / t_max_h / drop_pads are the per-run
    data-quality cuts (cfg.MIN_AMP / cfg.T_MAX_H / cfg.NOISY_PADS)."""
    cen, hit_ev, _ = p2io.event_centroids(hits_dir, channel_table,
                                          min_amp=min_amp,
                                          leading_pad=leading_pad,
                                          t_max_h=t_max_h,
                                          drop_pads=drop_pads)
    return cen, set(int(e) for e in hit_ev)


def filter_events_by_time(df, m3_dir, t_max_h, id_col='eventId'):
    """Drop rows whose event happened later than t_max_h hours after run
    start (event times from the m3 tracking files, same clock as the hits'
    trigger_timestamp_ns). Used so a detector trip mid-run does not leave
    dead-time rays in the efficiency denominator."""
    if t_max_h is None:
        return df
    # t_sec (= trigger_timestamp_ns/1e9) starts at ~0 at run start, so the cut
    # is on absolute run time. Do NOT reference the first m3 event: the
    # on-the-fly ray processor can start hours into a run (det4 7-15: first
    # ray file begins 2.56 h in), which would shift the window.
    et = load_event_times(m3_dir)
    good = set(et.loc[et['t_sec'] <= t_max_h * 3600.0, 'eventId'])
    kept = df[df[id_col].isin(good)]
    print(f'  t_max_h={t_max_h:g} h: kept {len(kept):,}/{len(df):,} '
          f'{id_col} rows.')
    return kept


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
