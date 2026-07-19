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


def hot_pad_mask(counts, ratio=5.0, min_n=30):
    """Boolean mask of 'hot' pads from per-pad hit counts.

    A pad is hot when it carries more than `ratio` times the median occupancy
    of the pads that fired at all AND at least `min_n` hits (so a low-stats
    sub_run cannot flag Poisson fluctuations). Cosmic illumination is smooth,
    so genuine pads stay within ~2-3x of the median; constantly-firing pads
    (oscillating channels, HV pickup) sit an order of magnitude above it
    without ever tripping the spark/burst vetoes."""
    counts = np.asarray(counts, dtype=float)
    fired = counts > 0
    if not fired.any():
        return np.zeros(len(counts), dtype=bool)
    med = float(np.median(counts[fired]))
    return (counts > ratio * med) & (counts >= min_n)


_HOT_CACHE = {}


def auto_hot_pads(hits_dir, channel_table, min_amp=0.0, t_max_h=None,
                  ratio=5.0, min_n=30, verbose=True):
    """Auto-detect constantly-firing pads in a sub_run (see hot_pad_mask).

    Streamed occupancy pre-pass: only per-pad counts are held (one row per
    pad), so memory stays bounded like every other reduction here. The result
    is cached per (hits_dir, cuts) so stages that load the same sub_run more
    than once in a process scan it only once.

    Returns a sorted tuple of hot channel_ids (empty if ratio is falsy)."""
    if not ratio:
        return ()
    key = (hits_dir, float(min_amp), t_max_h, float(ratio), int(min_n))
    if key in _HOT_CACHE:
        return _HOT_CACHE[key]
    counts = None
    for df in iter_hits(hits_dir, ['channel', 'feu'],
                        channel_table.attrs['feus'], progress=False,
                        t_max_h=t_max_h, min_amp=min_amp):
        h = pmap.attach_pads_to_hits(df, channel_table)
        h = h[h['mapped'] & h['pad_cx'].notna()]
        del df
        if not len(h):
            continue
        c = h.groupby('channel_id').size()
        counts = c if counts is None else counts.add(c, fill_value=0)
    if counts is None:
        _HOT_CACHE[key] = ()
        return ()
    hot = hot_pad_mask(counts.to_numpy(), ratio=ratio, min_n=min_n)
    ids = tuple(sorted(int(i) for i in counts.index[hot]))
    if verbose and ids:
        med = float(np.median(counts[counts > 0]))
        det = ', '.join(f'{i} ({int(counts[i]):,} hits, '
                        f'{counts[i] / med:.1f}x median)' for i in ids)
        print(f'Hot-pad cut (>{ratio:g}x median occupancy, >={min_n} hits): '
              f'masking {len(ids)} pad(s): {det}', flush=True)
    _HOT_CACHE[key] = ids
    return ids


def drop_pads_for(cfg, channel_table, hits_dir=None):
    """Channels to drop at load time: manual cfg.NOISY_PADS plus the pads
    auto-flagged by the occupancy cut (cfg.HOT_PAD_RATIO; 0 disables).
    Pass hits_dir for stages that loop per-sub_run directories (11/16)."""
    drop = set(cfg.NOISY_PADS)
    ratio = getattr(cfg, 'HOT_PAD_RATIO', 0.0)
    if ratio:
        drop |= set(auto_hot_pads(hits_dir or cfg.combined_hits_dir,
                                  channel_table, min_amp=cfg.MIN_AMP,
                                  t_max_h=cfg.T_MAX_H, ratio=ratio))
    return tuple(sorted(drop))


def pad_amp_stats(hits_dir, channel_table, min_amp=0.0, t_max_h=None,
                  drop_pads=(), exclude_events=None, amax=20000.0, dbin=2.0):
    """Streamed pad-amplitude summary for one sub_run (gain proxy vs HV).

    Mean and its standard error come from exact running sums; the median comes
    from a fine ADC histogram (0..amax, dbin-wide bins) so no per-hit array is
    ever held -- memory stays bounded like every other reduction here (the det4
    chunks carry ~18M hits). Amplitudes above amax are clipped into the top bin
    (P2 saturates well below 20k ADC, so this only guards pathological values).

    drop_pads (hot/noisy pads) and exclude_events (e.g. spark-vetoed eventIds,
    so the amplitude matches the efficiency sample) are applied before the
    reduction. Returns dict(n, mean, sem, median).
    """
    nb = int(round(amax / dbin))
    edges = np.linspace(0.0, amax, nb + 1)
    hist = np.zeros(nb, dtype=np.int64)
    excl = set(exclude_events) if exclude_events is not None else None
    n = 0
    s = 0.0
    s2 = 0.0
    for df in iter_hits(hits_dir, ['eventId', 'channel', 'amplitude', 'feu'],
                        channel_table.attrs['feus'], progress=False,
                        t_max_h=t_max_h, min_amp=min_amp):
        h = pmap.attach_pads_to_hits(df, channel_table)
        h = h[h['mapped'] & h['pad_cx'].notna()]
        if drop_pads:
            h = h[~h['channel_id'].isin(set(drop_pads))]
        if excl is not None and len(h):
            h = h[~h['eventId'].isin(excl)]
        del df
        if not len(h):
            continue
        a = h['amplitude'].to_numpy(dtype=np.float64)
        n += len(a)
        s += float(a.sum())
        s2 += float((a * a).sum())
        hist += np.histogram(np.clip(a, 0.0, amax), bins=edges)[0]
    if n == 0:
        return dict(n=0, mean=np.nan, sem=np.nan, median=np.nan)
    mean = s / n
    var = max(s2 / n - mean * mean, 0.0)
    sem = float(np.sqrt(var / n))
    cum = np.cumsum(hist)
    k = int(np.searchsorted(cum, 0.5 * n))
    k = min(k, nb - 1)
    median = 0.5 * (edges[k] + edges[k + 1])
    return dict(n=n, mean=mean, sem=sem, median=float(median))


def pad_time_spread(hits_dir, channel_table, sig_amp=300.0, min_amp=0.0,
                    t_max_h=None, drop_pads=(), exclude_events=None,
                    tmax=2000.0, dbin=2.0):
    """Streamed spread of the hit peak-time (time_of_max) for signal-band pad
    hits (amplitude >= sig_amp), a data proxy for the drift-time spread.

    For cosmic tracks the ionisation is deposited ~uniformly across the drift
    gap, so electrons arrive over a window ~ d_gap / v_drift: the WIDTH of the
    time_of_max distribution measures that drift-time spread. The DREAM shaping
    offset (~140 ns rise) and the per-event trigger phase are common to all
    hits in a point, so they cancel in the spread (they shift, not widen it),
    leaving diffusion + trigger-phase jitter as a roughly constant floor across
    the scan -- so the TREND and the location of the minimum are robust even
    though the absolute value carries that floor.

    Percentiles come from a fine time_of_max histogram (0..tmax, dbin-wide) so
    no per-hit array is held. Returns dict(n, p10, p50, p90, spread, iqr) [ns].
    """
    nb = int(round(tmax / dbin))
    edges = np.linspace(0.0, tmax, nb + 1)
    hist = np.zeros(nb, dtype=np.int64)
    excl = set(exclude_events) if exclude_events is not None else None
    thr = max(float(min_amp), float(sig_amp))
    for df in iter_hits(hits_dir, ['eventId', 'channel', 'amplitude', 'feu',
                                   'time_of_max'],
                        channel_table.attrs['feus'], progress=False,
                        t_max_h=t_max_h, min_amp=min_amp):
        h = pmap.attach_pads_to_hits(df, channel_table)
        h = h[h['mapped'] & h['pad_cx'].notna() & (h['amplitude'] >= thr)]
        if drop_pads:
            h = h[~h['channel_id'].isin(set(drop_pads))]
        if excl is not None and len(h):
            h = h[~h['eventId'].isin(excl)]
        del df
        if not len(h):
            continue
        tom = h['time_of_max'].to_numpy(dtype=np.float64)
        hist += np.histogram(np.clip(tom, 0.0, tmax), bins=edges)[0]
    n = int(hist.sum())
    if n == 0:
        return dict(n=0, p10=np.nan, p50=np.nan, p90=np.nan,
                    spread=np.nan, iqr=np.nan)
    cum = np.cumsum(hist)
    ctr = 0.5 * (edges[:-1] + edges[1:])

    def pct(q):
        return float(ctr[min(int(np.searchsorted(cum, q * n)), nb - 1)])
    p10, p25, p50, p75, p90 = (pct(q) for q in (0.10, 0.25, 0.50, 0.75, 0.90))
    return dict(n=n, p10=p10, p50=p50, p90=p90,
                spread=p90 - p10, iqr=p75 - p25)


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
