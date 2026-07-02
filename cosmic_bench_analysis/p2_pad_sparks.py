#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p2_pad_sparks.py

Hit-level (per-pad) spark flagging for the P2 pad detector.

Motivation
----------
The HV-monitor veto (p2_sparks) removes whole-detector discharges: events taken
while the *total* mesh current spikes. But individual pads can micro-spark on
their own -- localised discharges too small to move the integrated mesh current
above threshold, and shorter than the 2 s HV logging interval, so they are
invisible to the power-supply record and survive the HV veto. On P2 det1 they
show up as specific "hot" pads (a cluster on connector 6) with:

  * an abnormal fraction of saturated / very-high-amplitude hits
    (discharge charge, not minimum-ionising cosmics), and
  * an elevated, time-clustered hit rate,

while ordinary large outer pads merely have a higher rate (bigger pad area) but
a *normal* amplitude/saturation profile.

So the discriminator is the discharge signature -- saturation rate and
high-amplitude fraction -- NOT the raw rate (which is confounded by pad area).
A pad is flagged when both are robust outliers above the detector-wide median.

This module only *computes and flags*; 08_pad_spark_qa.py renders the review so
the flag can be judged before deciding whether to mask these pads in the
analysis (like the dead-connector drop).
"""

import numpy as np
import pandas as pd

HI_ADC = 1000.0        # a hit above this amplitude is "high" (discharge-like)
NOISE_ADC = 150.0      # a hit below this amplitude is near-pedestal / threshold


def robust_z(x):
    """Median/MAD z-score (1.4826*MAD ~ sigma for a Gaussian)."""
    x = np.asarray(x, dtype=float)
    m = np.median(x)
    mad = np.median(np.abs(x - m)) * 1.4826
    return (x - m) / (mad if mad > 0 else 1.0)


def pad_spark_metrics(hits, channel_table, n_events, hi_adc=HI_ADC,
                      noise_adc=NOISE_ADC):
    """Per-pad amplitude/rate metrics from mapped hits.

    hits must already carry channel_id / pad geometry (attach_pads_to_hits) and
    should be the hits you actually analyse (i.e. after the HV-monitor veto).

    Returns a DataFrame indexed by channel_id with:
      n, fire_frac, sat_rate, hi_frac, lo_frac, mean_amp, med_amp,
      connector_N, radius, pad_cx, pad_cy
    where hi_frac = fraction of hits > hi_adc (discharge-like) and
    lo_frac = fraction of hits < noise_adc (near-pedestal / noise-like).
    """
    g = hits.groupby('channel_id')
    m = pd.DataFrame({
        'n': g.size(),
        'fire_frac': g['eventId'].nunique() / max(n_events, 1),
        'sat_rate': g['saturated'].mean(),
        'hi_frac': g['amplitude'].apply(lambda a: float((a > hi_adc).mean())),
        'lo_frac': g['amplitude'].apply(lambda a: float((a < noise_adc).mean())),
        'mean_amp': g['amplitude'].mean(),
        'med_amp': g['amplitude'].median(),
    })
    geom = channel_table.drop_duplicates('channel_id').set_index('channel_id')[
        ['connector_N', 'radius', 'pad_cx', 'pad_cy']]
    return m.join(geom)


def flag_noise_pads(metrics, noise_amp=NOISE_ADC, z_rate=3.0, min_hits=20):
    """Flag noisy channels: pads that fire at an abnormally high rate but with a
    near-pedestal median amplitude (electronic noise, not a real Landau signal).

    A pad is noise if its median amplitude < noise_amp AND its firing rate is a
    robust outlier (fire_z > z_rate). This is the amplitude-opposite of a spark
    (low amp, not saturating), so masking noise first cleanly separates the two
    populations. Adds 'fire_z' and 'noise_pad' columns.
    """
    m = metrics.copy()
    ok = m['n'] >= min_hits
    m['fire_z'] = np.nan
    if ok.any():
        m.loc[ok, 'fire_z'] = robust_z(m.loc[ok, 'fire_frac'])
    m['noise_pad'] = ok & (m['med_amp'] < noise_amp) & (m['fire_z'] > z_rate)
    return m


def flag_spark_pads(metrics, z_sat=4.0, z_hi=3.0, min_hits=20):
    """Flag pads whose saturation rate AND high-amplitude fraction are both
    robust outliers. Pads with < min_hits are never flagged (too few stats).

    Adds columns sat_z, hi_z, spark_pad (bool) and returns the DataFrame.
    Robust z-scores are computed only over pads with >= min_hits so a few
    near-silent pads do not distort the median/MAD.
    """
    m = metrics.copy()
    ok = m['n'] >= min_hits
    m['sat_z'] = np.nan
    m['hi_z'] = np.nan
    if ok.any():
        m.loc[ok, 'sat_z'] = robust_z(m.loc[ok, 'sat_rate'])
        m.loc[ok, 'hi_z'] = robust_z(m.loc[ok, 'hi_frac'])
    m['spark_pad'] = ok & (m['sat_z'] > z_sat) & (m['hi_z'] > z_hi)
    return m


class PadSparkMask:
    """Per-pad noise + spark flags for one run, ready to review or apply.

    Two populations are flagged in order:
      1. noise pads  -- high rate, near-pedestal amplitude (electronic noise);
      2. spark pads  -- saturated / high-amplitude discharge signature, judged
                        on the pads that are NOT noise.
    Masking noise first keeps the (opposite-amplitude) noise channels from
    contaminating the discharge statistics.
    """

    def __init__(self, metrics, z_sat, z_hi, min_hits, hi_adc, noise_amp, z_rate):
        self.metrics = metrics
        self.z_sat = z_sat
        self.z_hi = z_hi
        self.min_hits = min_hits
        self.hi_adc = hi_adc
        self.noise_amp = noise_amp
        self.z_rate = z_rate

    @classmethod
    def from_hits(cls, hits, channel_table, n_events, z_sat=4.0, z_hi=3.0,
                  min_hits=20, hi_adc=HI_ADC, noise_amp=NOISE_ADC, z_rate=3.0):
        met = pad_spark_metrics(hits, channel_table, n_events, hi_adc, noise_amp)
        # 1) noise first
        met = flag_noise_pads(met, noise_amp, z_rate, min_hits)
        # 2) discharge on the non-noise pads (recompute robust z excluding noise)
        keep = ~met['noise_pad']
        sub = flag_spark_pads(met[keep], z_sat, z_hi, min_hits)
        met['sat_z'] = np.nan
        met['hi_z'] = np.nan
        met['spark_pad'] = False
        met.loc[keep, ['sat_z', 'hi_z', 'spark_pad']] = sub[
            ['sat_z', 'hi_z', 'spark_pad']]
        return cls(met, z_sat, z_hi, min_hits, hi_adc, noise_amp, z_rate)

    @property
    def noise(self):
        return self.metrics[self.metrics['noise_pad']].sort_values(
            'fire_frac', ascending=False)

    @property
    def flagged(self):
        """Discharge/spark-flagged pads, worst (highest sat_rate) first."""
        return self.metrics[self.metrics['spark_pad']].sort_values(
            'sat_rate', ascending=False)

    @property
    def noise_ids(self):
        return set(int(c) for c in self.noise.index)

    @property
    def flagged_ids(self):
        return set(int(c) for c in self.flagged.index)

    @property
    def all_masked_ids(self):
        return self.noise_ids | self.flagged_ids

    def apply_hits(self, hits, which='all'):
        """Drop hits on masked pads. which in {'all','noise','spark'}.
        Returns (kept, n_removed)."""
        ids = {'all': self.all_masked_ids, 'noise': self.noise_ids,
               'spark': self.flagged_ids}[which]
        if not ids:
            return hits, 0
        keep = ~hits['channel_id'].isin(ids)
        return hits[keep].copy(), int((~keep).sum())

    def summary(self):
        n_active = int((self.metrics['n'] >= self.min_hits).sum())
        nz, f = self.noise, self.flagged
        lines = [
            f'Per-pad flags over {n_active} active pads (min_hits={self.min_hits}):',
            f'  NOISE pads (med_amp<{self.noise_amp:g} ADC & fire_z>{self.z_rate:g}) '
            f': {len(nz)}',
        ]
        if len(nz):
            lines.append('    per connector: ' + ', '.join(
                f'c{int(c)}={int(v)}' for c, v in nz.groupby('connector_N').size().items()))
        lines.append(f'  SPARK pads (sat_z>{self.z_sat:g} & hi_z>{self.z_hi:g}, '
                     f'hi_adc={self.hi_adc:g}) : {len(f)}')
        if len(f):
            lines.append('    per connector: ' + ', '.join(
                f'c{int(c)}={int(v)}' for c, v in f.groupby('connector_N').size().items()))
            lines.append(f'    worst pad   : ch_id {int(f.index[0])} '
                         f'(conn {int(f.iloc[0].connector_N)}, '
                         f'sat_rate {f.iloc[0].sat_rate:.3f}, '
                         f'hi_frac {f.iloc[0].hi_frac:.3f})')
        return '\n'.join(lines)
