#!/usr/bin/env python3
"""p2_deform.py -- shared pieces of the deformation study (stages 43-5x).

load(run_key)      the stage-43 waveform dataset with per-pulse shape features
features(WF, ...)  shape features of (n, 32) pedestal/CNS-subtracted waveforms

All times are in ns relative to the trigger: t = sample * 60 + ftst * 10 (the
p2_waveforms convention, validated in stage 13). Shape features (differences of
two times of the same pulse) do not depend on that convention.

Feature definitions (leading edge = walking left from the maximum, trailing
edge = walking right; crossings linearly interpolated between 60 ns samples):
  amp        parabola-interpolated maximum [ADC]
  t_pk       time of that maximum
  t_l10..50  leading-edge crossing of 10/20/30/50 % of amp
  t_t50,t_t30 trailing-edge crossing of 50/30 % of amp (NaN if beyond the window)
  rise       t_pk - t_l30         (the quantity of stage 24)
  r2050      t_l50 - t_l20        (leading-edge slope only)
  fwhm       t_t50 - t_l50
  fall       t_t50 - t_pk
  area       sum of samples from (l10 sample - 1) to the end of the window [ADC*sample]
  wid        area / amp * 60      (effective width, ns; includes the ion tail in the window)
  t_cen      amplitude-weighted mean time of the samples above 20 % of amp
  base       mean of the first 3 samples (pre-pulse baseline after CNS)
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2_qa_config as qa       # noqa: E402

TS = 60.0      # ns per sample
TF = 10.0      # ns per ftst unit


def _cross(w, level, i0, step):
    """Crossing of `level` walking from sample i0 in direction step (-1/+1)."""
    n, m = w.shape
    i = i0.copy()
    t = np.full(n, np.nan)
    done = np.zeros(n, bool)
    rows = np.arange(n)
    for _ in range(m):
        j = i + step
        ok = (~done) & (j >= 0) & (j < m)
        if not ok.any():
            break
        yj = np.where(ok, w[rows, np.clip(j, 0, m - 1)], np.nan)
        yi = w[rows, i]
        hit = ok & (yj <= level)
        den = np.where(np.abs(yi - yj) > 1e-9, yi - yj, np.nan)
        t[hit] = i[hit] + step * (yi[hit] - level[hit]) / den[hit]
        done |= hit | ~ok
        i = np.where(ok & ~hit, j, i)
    return t


def features(WF, ftst=None):
    """Shape features of waveforms WF (n, 32). Returns a DataFrame (times in ns)."""
    w = np.asarray(WF, np.float64)
    n, m = w.shape
    rows = np.arange(n)
    imax = w.argmax(1)
    inner = (imax > 0) & (imax < m - 1)
    y1 = w[rows, np.clip(imax - 1, 0, m - 1)]; y2 = w[rows, imax]; y3 = w[rows, np.clip(imax + 1, 0, m - 1)]
    den = y1 - 2 * y2 + y3
    off = np.where(inner & (np.abs(den) > 1e-9), 0.5 * (y1 - y3) / np.where(np.abs(den) > 1e-9, den, 1), 0.0)
    off = np.clip(off, -1, 1)
    amp = np.where(inner, y2 - 0.25 * (y1 - y3) * off, y2)
    tpk = imax + off
    out = {'amp': amp, 'imax': imax}
    lead = {}
    for f in (10, 20, 30, 50):
        lead[f] = _cross(w, f / 100.0 * amp, imax, -1)
    trail = {f: _cross(w, f / 100.0 * amp, imax, +1) for f in (50, 30)}
    t0 = (ftst.astype(float) * TF if ftst is not None else 0.0)
    out['t_pk'] = tpk * TS + t0
    for f in lead:
        out[f't_l{f}'] = lead[f] * TS + t0
    for f in trail:
        out[f't_t{f}'] = trail[f] * TS + t0
    out['rise'] = (tpk - lead[30]) * TS
    out['r2050'] = (lead[50] - lead[20]) * TS
    out['r1050'] = (lead[50] - lead[10]) * TS
    out['fwhm'] = (trail[50] - lead[50]) * TS
    out['fall'] = (trail[50] - tpk) * TS
    out['fall30'] = (trail[30] - tpk) * TS
    # area from one sample before the 10 % crossing to the end of the window
    i0 = np.clip(np.floor(np.nan_to_num(lead[10], nan=0.0)).astype(int) - 1, 0, m - 1)
    idx = np.arange(m)[None, :]
    msk = idx >= i0[:, None]
    area = (w * msk).sum(1)
    out['area'] = area
    out['wid'] = area / np.where(amp > 0, amp, np.nan) * TS
    above = w > 0.2 * amp[:, None]
    # contiguous region around the maximum only
    lft = np.floor(np.nan_to_num(lead[20], nan=0.0)).astype(int)
    rgt = np.ceil(np.nan_to_num(trail[30], nan=m - 1)).astype(int)
    reg = (idx >= lft[:, None]) & (idx <= rgt[:, None]) & above
    s = (w * reg).sum(1)
    out['t_cen'] = np.where(s > 0, (w * reg * idx).sum(1) / np.where(s > 0, s, 1), np.nan) * TS + t0
    out['base'] = w[:, :3].mean(1)
    out['last'] = w[:, -2:].mean(1)
    return pd.DataFrame(out)


def load(run_key='det2_long_repos', with_features=True):
    """Stage-43 dataset: (R, WF, cfg). R gets the feature columns and
    'good' = not vetoed, inside cfg time window."""
    cfg = qa.get_config(run_key)
    out = cfg.out_dir('43_deformation')
    R = pd.read_pickle(f'{out}/wf_table.pkl')
    WF = np.load(f'{out}/wf.npy')
    ok = ~R.vetoed.to_numpy()
    if cfg.T_MAX_H is not None:
        ok &= ~(R.t_h.to_numpy() > cfg.T_MAX_H)
    R['good'] = ok
    if with_features:
        F = features(WF, R.ftst.to_numpy())
        for c in F.columns:
            R[c] = F[c].to_numpy()
    reg = pd.read_csv(f'{out}/regions.csv')
    R = R.merge(reg[['channel_id', 'reg8', 'reg24', 'reg48']], on='channel_id', how='left')
    return R, WF, cfg


def s68(v):
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    if len(v) < 5:
        return np.nan
    return 0.5 * (np.percentile(v, 84.13) - np.percentile(v, 15.87))


def boot(v, fn=np.median, n=200, rng=None):
    """Bootstrap error of fn(v)."""
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    if len(v) < 5:
        return np.nan
    rng = rng or np.random.default_rng(1)
    return np.std([fn(v[rng.integers(0, len(v), len(v))]) for _ in range(n)])
