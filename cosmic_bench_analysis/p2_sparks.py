#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p2_sparks.py

Flag HV sparks on the P2 mesh from hv_monitor.csv and turn them into an
event-level veto that every analysis stage can apply.

Signal
------
The P2 amplification mesh is a single CAEN HV channel (cfg.SPARK_CHANNEL, e.g.
'1:0' at 420 V). A spark is a micro-discharge across the gap: it shows up as a
brief spike in the monitored current `imon` (and, for a hard discharge, a dip in
`vmon`), while the crate keeps the channel powered. Baseline imon on this mesh is
~0.07 uA; sparks reach several uA up to the ISet compliance (~10 uA). We flag a
sample as spark-like when imon >= cfg.SPARK_IMON_THR (default 2 uA, well above
the ~0.14 uA 90th-percentile baseline), then group consecutive spark samples
into spark *intervals* (one physical discharge / burst).

Time alignment
--------------
hv_monitor.csv carries wall-clock timestamps sampled every ~2 s. DAQ hits carry
`trigger_timestamp_ns`, which runs continuously across the combined-hits files
from 0 at run start. The HV monitor also starts at run start, so both share a
common t=0: HV time = (timestamp - timestamp[0]) in seconds equals event time =
trigger_timestamp_ns / 1e9. (Verified: the two spans agree to ~0.3 %.)

Veto
----
Each spark interval [t_start, t_end] is padded by cfg.SPARK_GUARD_BEFORE (clock
skew) and cfg.SPARK_GUARD_AFTER (post-spark recovery, during which the detector
response is degraded); overlapping padded intervals are merged. `apply_veto`
drops every event whose trigger time falls inside a vetoed interval, so sparks
are removed identically from every stage. `live_fraction` reports the surviving
live-time for efficiency normalisation.

Typical use in a stage
----------------------
    sp = p2_sparks.SparkVeto.from_cfg(cfg)          # detect + build intervals
    hits, n_removed = sp.apply(hits)                # drop spark-time events
    print(sp.summary())                             # text block for the log
"""

import numpy as np
import pandas as pd

# numpy 2.x renamed trapz -> trapezoid
_trapz = getattr(np, 'trapezoid', getattr(np, 'trapz', None))


# --------------------------------------------------------------------------- #
# Loading
# --------------------------------------------------------------------------- #
def load_hv(hv_csv, channel):
    """Load hv_monitor.csv for one channel. Returns a DataFrame with columns
    t (seconds since the first sample), imon [uA], vmon [V], power (0/1)."""
    df = pd.read_csv(hv_csv)
    ts = pd.to_datetime(df['timestamp'])
    out = pd.DataFrame({
        't': (ts - ts.iloc[0]).dt.total_seconds().to_numpy(),
        'imon': df[f'{channel} imon'].to_numpy(dtype=float),
        'vmon': df[f'{channel} vmon'].to_numpy(dtype=float),
        'power': df[f'{channel} power'].to_numpy(dtype=float),
    })
    return out


# --------------------------------------------------------------------------- #
# Spark detection
# --------------------------------------------------------------------------- #
def detect_sparks(hv, i_thr, merge_gap=6.0):
    """Group consecutive over-threshold mesh samples into spark intervals.

    hv        : DataFrame from load_hv (needs t, imon).
    i_thr     : imon threshold [uA].
    merge_gap : samples separated by <= this [s] are treated as one spark burst
                (default 6 s ~ up to two missed 2 s samples).

    Returns a DataFrame, one row per spark, with columns:
        t_start, t_end, t_peak, dur, peak_imon, charge (uA*s = uC), n_samples.
    """
    t = hv['t'].to_numpy()
    im = hv['imon'].to_numpy()
    spark = im >= i_thr
    if not spark.any():
        return pd.DataFrame(columns=['t_start', 't_end', 't_peak', 'dur',
                                     'peak_imon', 'charge', 'n_samples'])

    # index runs of consecutive over-threshold samples
    idx = np.flatnonzero(spark)
    breaks = np.flatnonzero(np.diff(idx) > 1)
    starts = np.concatenate([[idx[0]], idx[breaks + 1]])
    ends = np.concatenate([idx[breaks], [idx[-1]]])           # inclusive

    # merge runs whose time gap is within merge_gap (one flickering discharge)
    groups = [[starts[0], ends[0]]]
    for s, e in zip(starts[1:], ends[1:]):
        if t[s] - t[groups[-1][1]] <= merge_gap:
            groups[-1][1] = e
        else:
            groups.append([s, e])

    rows = []
    for s, e in groups:
        tt, ii = t[s:e + 1], im[s:e + 1]
        pk = int(np.argmax(ii))
        charge = float(_trapz(ii, tt)) if len(tt) > 1 else float(ii[0] * 2.0)
        rows.append({
            't_start': float(tt[0]), 't_end': float(tt[-1]),
            't_peak': float(tt[pk]), 'dur': float(tt[-1] - tt[0]),
            'peak_imon': float(ii[pk]), 'charge': charge,
            'n_samples': int(len(tt)),
        })
    return pd.DataFrame(rows, columns=['t_start', 't_end', 't_peak', 'dur',
                                       'peak_imon', 'charge', 'n_samples'])


# --------------------------------------------------------------------------- #
# Veto intervals
# --------------------------------------------------------------------------- #
def veto_intervals(sparks, guard_before, guard_after):
    """Pad each spark by (guard_before, guard_after) and merge overlaps.
    Returns a sorted list of (lo, hi) run-time intervals [s]."""
    if len(sparks) == 0:
        return []
    lo = sparks['t_start'].to_numpy() - guard_before
    hi = sparks['t_end'].to_numpy() + guard_after
    order = np.argsort(lo)
    lo, hi = lo[order], hi[order]
    merged = [[lo[0], hi[0]]]
    for a, b in zip(lo[1:], hi[1:]):
        if a <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    return [(float(a), float(b)) for a, b in merged]


def _in_intervals(t, intervals):
    """Boolean mask: is each t [s] inside any (lo, hi) interval?"""
    mask = np.zeros(len(t), dtype=bool)
    for lo, hi in intervals:
        mask |= (t >= lo) & (t <= hi)
    return mask


# --------------------------------------------------------------------------- #
# Convenience wrapper
# --------------------------------------------------------------------------- #
class SparkVeto:
    """Detected sparks + veto intervals for one run, ready to apply to hits."""

    def __init__(self, hv, sparks, intervals, channel, i_thr,
                 guard_before, guard_after):
        self.hv = hv
        self.sparks = sparks
        self.intervals = intervals
        self.channel = channel
        self.i_thr = i_thr
        self.guard_before = guard_before
        self.guard_after = guard_after
        self.t_run = float(hv['t'].max()) if len(hv) else 0.0

    @classmethod
    def from_cfg(cls, cfg):
        hv = load_hv(cfg.hv_monitor_csv, cfg.SPARK_CHANNEL)
        sparks = detect_sparks(hv, cfg.SPARK_IMON_THR)
        intervals = veto_intervals(sparks, cfg.SPARK_GUARD_BEFORE,
                                   cfg.SPARK_GUARD_AFTER)
        return cls(hv, sparks, intervals, cfg.SPARK_CHANNEL, cfg.SPARK_IMON_THR,
                   cfg.SPARK_GUARD_BEFORE, cfg.SPARK_GUARD_AFTER)

    # -- veto application --------------------------------------------------- #
    @property
    def vetoed_seconds(self):
        return float(sum(hi - lo for lo, hi in self.intervals))

    def live_fraction(self):
        if self.t_run <= 0:
            return 1.0
        return max(0.0, 1.0 - self.vetoed_seconds / self.t_run)

    def event_mask(self, trigger_ns):
        """Boolean 'keep' mask for an array of trigger_timestamp_ns."""
        t = np.asarray(trigger_ns, dtype=float) / 1e9
        return ~_in_intervals(t, self.intervals)

    def apply(self, df, ts_col='trigger_timestamp_ns'):
        """Return (kept_df, n_removed_rows), dropping rows in a spark window."""
        if not self.intervals or ts_col not in df.columns:
            return df, 0
        keep = self.event_mask(df[ts_col].to_numpy())
        return df[keep].copy(), int((~keep).sum())

    def vetoed_event_ids(self, df, ts_col='trigger_timestamp_ns',
                         id_col='eventId'):
        """Set of eventIds that fall inside a spark window (for id-based cuts)."""
        if not self.intervals or ts_col not in df.columns:
            return set()
        bad = ~self.event_mask(df[ts_col].to_numpy())
        return set(df.loc[bad, id_col].astype(int))

    def vetoed_ids_from_hits(self, hits_dir, feus):
        """Read event times straight from the combined-hits ROOT files and
        return the set of eventIds inside a spark window. Used by stages that
        cut by eventId (e.g. efficiency, which must drop the same events from
        both the M3 ray list and the P2 centroids)."""
        import glob
        import os
        import uproot
        if not self.intervals:
            return set()
        feu_set = set(feus)
        ids = []
        for fp in sorted(glob.glob(os.path.join(hits_dir, '*.root'))):
            a = uproot.open(f'{fp}:hits').arrays(
                ['eventId', 'trigger_timestamp_ns', 'feu'], library='pd')
            a = a[a['feu'].isin(feu_set)]
            bad = ~self.event_mask(a['trigger_timestamp_ns'].to_numpy())
            ids.extend(a.loc[bad, 'eventId'].astype(int).tolist())
        return set(ids)

    # -- reporting ---------------------------------------------------------- #
    def summary(self):
        n = len(self.sparks)
        rate = n / (self.t_run / 3600.0) if self.t_run > 0 else 0.0
        lines = [
            f'HV spark veto  (channel {self.channel}, imon>={self.i_thr:g} uA)',
            f'  run duration        : {self.t_run/3600:.2f} h ({self.t_run:.0f} s)',
            f'  sparks detected     : {n}  ({rate:.1f} / h)',
            f'  guard window        : -{self.guard_before:g} / +{self.guard_after:g} s',
            f'  vetoed live-time    : {self.vetoed_seconds:.0f} s '
            f'({100*(1-self.live_fraction()):.2f}% of run)',
        ]
        if n:
            lines.append(f'  peak imon (max)     : {self.sparks["peak_imon"].max():.2f} uA')
            lines.append(f'  total spark charge  : {self.sparks["charge"].sum():.1f} uC')
        return '\n'.join(lines)
