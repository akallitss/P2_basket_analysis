#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p2_io.py

Memory-safe streaming access to the combined-hits ROOT files.

The det4 long run carries ~17.7 M hits per 400 MB chunk (~680 hits/event at
mesh 480 V), so the old pattern -- load every chunk into one pandas DataFrame,
then reduce -- needs many GB of RAM and OOM-kills the analysis on a 14 GB
machine with no swap. All stages now iterate chunk by chunk with `iter_hits`
and reduce each chunk (per-pad / per-event aggregation) before the next one is
read, so peak memory stays at a single chunk (<~1 GB) no matter how many
chunks the run has.

`event_centroids` is the shared streaming reduction to per-event pad
centroids used by 03 (via its spark-vetoed wrapper) and by p2_align (06/11/12).
"""

import os
import glob

import numpy as np
import pandas as pd
import uproot

import p2_mapping as pmap

# The native combined-hits dtypes are wasteful (uint64 eventId, int32 feu...);
# downcast on read so a transient chunk costs about half the naive frame.
_DTYPES = {
    'eventId': np.int64,           # continuous DAQ counter; keep wide
    'trigger_timestamp_ns': np.uint64,
    'channel': np.uint16,
    'amplitude': np.float32,
    'saturated': np.bool_,
    'feu': np.int16,
    'time': np.float32,
}


def hit_files(hits_dir):
    files = sorted(glob.glob(os.path.join(hits_dir, '*.root')))
    if not files:
        raise FileNotFoundError(f'No combined-hits ROOT files in {hits_dir}')
    return files


def iter_hits(hits_dir, branches, feus=None, progress=True,
              t_max_h=None, min_amp=0.0):
    """Yield one downcast, FEU-filtered DataFrame per combined-hits chunk.

    Reduce each yielded frame before requesting the next one; never
    concatenate the raw frames (that is exactly the OOM this module removes).

    Per-run data-quality cuts (cfg.T_MAX_H / cfg.MIN_AMP):
      t_max_h  drop hits later than this many hours after the first trigger
               (e.g. the detector tripped mid-run); chunks are time-ordered so
               iteration stops early once the window is exhausted.
      min_amp  drop hits below this amplitude [ADC] (stale-pedestal noise
               floor). Applied BEFORE any spark/burst veto the caller runs, so
               burst pad counts see only signal-like hits.
    """
    files = hit_files(hits_dir)
    feu_set = set(feus) if feus is not None else None
    read = list(branches)
    if t_max_h is not None and 'trigger_timestamp_ns' not in read:
        read.append('trigger_timestamp_ns')
    if min_amp and 'amplitude' not in read:
        read.append('amplitude')
    t_end = None
    for i, fp in enumerate(files):
        if progress:
            print(f'    chunk {i + 1}/{len(files)}: {os.path.basename(fp)}',
                  flush=True)
        with uproot.open(fp) as f:
            arrs = f['hits'].arrays(read, library='np')
        df = pd.DataFrame({b: arrs[b].astype(_DTYPES.get(b, arrs[b].dtype),
                                              copy=False)
                           for b in read})
        del arrs
        if t_max_h is not None:
            ts = df['trigger_timestamp_ns'].astype(np.int64)
            if t_end is None:
                t_end = int(ts.iloc[0]) + int(t_max_h * 3.6e12)
            if ts.iloc[0] > t_end:   # chunks are time-ordered: nothing left
                if progress:
                    print(f'    t_max_h={t_max_h:g} h reached — skipping the '
                          f'remaining {len(files) - i} chunk(s).', flush=True)
                break
            if ts.iloc[-1] > t_end:
                df = df[ts <= t_end]
        m = np.ones(len(df), dtype=bool)
        if feu_set is not None and 'feu' in df.columns:
            m &= df['feu'].isin(feu_set).to_numpy()
        if min_amp:
            m &= (df['amplitude'].to_numpy() >= min_amp)
        if not m.all():
            df = df[m]
        extra = [c for c in read if c not in branches]
        if extra:
            df = df.drop(columns=extra)
        yield df


def _empty_cen():
    return pd.DataFrame({'eventId': pd.Series(dtype=np.int64),
                         'x_pad': pd.Series(dtype=np.float64),
                         'y_pad': pd.Series(dtype=np.float64),
                         'n_pad': pd.Series(dtype=np.int64)})


def event_centroids(hits_dir, channel_table, min_amp=0.0, leading_pad=False,
                    spark_veto=None, t_max_h=None, drop_pads=()):
    """Per-event P2 pad centroid (charge-weighted or leading-pad), streamed.

    Returns (cen, hit_events, veto_stats):
      cen        DataFrame [eventId, x_pad, y_pad, n_pad]
      hit_events sorted int64 ndarray of eventIds with any mapped P2 hit
      veto_stats {'n_rm': rows dropped, 'n_burst': burst events} (zeros if
                 spark_veto is None)

    Per-event partial sums are combined across chunks with a groupby-sum, so
    an event straddling a chunk boundary is still reduced exactly once.
    min_amp / t_max_h are applied at read time (see iter_hits), i.e. before
    the spark veto's burst counting.
    """
    branches = ['eventId', 'channel', 'amplitude', 'feu']
    if spark_veto is not None:
        branches.append('trigger_timestamp_ns')
    parts, ev_parts = [], []
    n_rm = n_burst = 0
    for df in iter_hits(hits_dir, branches, channel_table.attrs['feus'],
                        t_max_h=t_max_h, min_amp=min_amp):
        if spark_veto is not None:
            df, rm = spark_veto.apply(df)
            n_rm += rm
            n_burst += spark_veto.last_burst_events
        h = pmap.attach_pads_to_hits(df, channel_table)
        h = h[h['mapped'] & h['pad_cx'].notna()]
        if drop_pads:   # known-noisy channels (cfg.NOISY_PADS)
            h = h[~h['channel_id'].isin(set(drop_pads))]
        del df
        if not len(h):
            continue
        ev_parts.append(h['eventId'].unique())
        if leading_pad:
            idx = h.groupby('eventId')['amplitude'].idxmax()
            parts.append(h.loc[idx, ['eventId', 'amplitude',
                                     'pad_cx', 'pad_cy']].copy())
        else:
            w = h['amplitude'].clip(lower=0).astype(np.float64)
            part = pd.DataFrame({'eventId': h['eventId'],
                                 '_wx': w * h['pad_cx'],
                                 '_wy': w * h['pad_cy'],
                                 '_w': w})
            parts.append(part.groupby('eventId').agg(
                _wx=('_wx', 'sum'), _wy=('_wy', 'sum'),
                _w=('_w', 'sum'), _n=('_w', 'size')))

    hit_events = (np.unique(np.concatenate(ev_parts)) if ev_parts
                  else np.array([], dtype=np.int64))
    if leading_pad:
        if parts:
            allp = pd.concat(parts, ignore_index=True)
            idx = allp.groupby('eventId')['amplitude'].idxmax()
            cen = allp.loc[idx, ['eventId', 'pad_cx', 'pad_cy']].rename(
                columns={'pad_cx': 'x_pad', 'pad_cy': 'y_pad'})
        else:  # e.g. the veto removed every event: keep float dtypes so the
            # isfinite cut below still works on the empty frame
            cen = _empty_cen()[['eventId', 'x_pad', 'y_pad']]
        cen['n_pad'] = 1
    else:
        if parts:
            tot = pd.concat(parts).groupby(level=0).sum()
            cen = pd.DataFrame({'x_pad': tot['_wx'] / tot['_w'],
                                'y_pad': tot['_wy'] / tot['_w'],
                                'n_pad': tot['_n'].astype(np.int64)}
                               ).reset_index()
        else:
            cen = _empty_cen()
    cen = cen[np.isfinite(cen['x_pad']) & np.isfinite(cen['y_pad'])]
    return cen, hit_events, {'n_rm': n_rm, 'n_burst': n_burst}
