#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
25_commissioning_qa.py

DAQ commissioning / validation QA for an SPS beam run: "is the whole telescope
actually taking good data?", answered from the combined-hits ROOT files, the HV
monitor CSV and the DAQ's own run bookkeeping, with an explicit PASS/WARN/FAIL
verdict per check.

This is the stage you run FIRST at the beam, before any physics stage. Where
20/21/22/23 measure the detectors, this one validates the *data taking*: every
detector reading out, every connector plugged, the trigger arriving at the
expected rate with the SPS spill structure, the signal sitting inside the
sampling window (latency), and the HV actually at its setpoint while the data
were written.

It is deliberately geometry-free -- it works in raw (FEU, channel) space -- so
it runs on the P2 stations and on the EIC uRWELL references alike, and it does
not depend on the pad map being correct. That makes it valid on day one of a
beam test, when the mapping/orientation is exactly what you are still checking.

Per detector, a detector is one CONTIGUOUS CHANNEL SLICE of one FEU, taken from
the run_config `dream_feus` wiring (Dream connector c -> channels
(c-1)*64 ... c*64-1). So the two uRWELLs sharing FEU 1 are separated correctly
(front = connectors 1-4, back = 5-8) without any hardcoding.

Checks (each -> PASS / WARN / FAIL, thresholds are CLI-tunable)
--------------------------------------------------------------
  readout        detector present in the data at all
  trigger_share  fraction of the run's triggers in which this detector has >=1
                 hit -- the per-plane "hit efficiency per trigger". A plane far
                 below its neighbours is unplugged, mis-thresholded, out of the
                 beam spot, or has the wrong latency.
  connectors     every Dream connector alive (a fully dead 64-channel block is
                 an unplugged/broken cable)
  dead_channels  fraction of channels far below the detector's median occupancy
  hot_channels   fraction of channels far above it (noise / stale pedestal)
  latency        the max-sample peak sits inside the sampling window and away
                 from both edges -- the signal is not clipped by a wrong Dream
                 latency (see 26 / the latency scan for the actual optimum)
  saturation     fraction of saturated hits
  rate_stability the trigger rate does not collapse mid-run (a FEU dropping out
                 or the DAQ backpressuring shows up here)
  hv_setpoint    fraction of the sub_run's monitored time with every channel at
                 its setpoint (catches data taken during the HV ramp)
  hv_current     imon excursions (sparking) during the sub_run

Products (<Analysis>/telescope/<run>/<sub_run>/25_commissioning_qa/):
  rates_<sub_run>.png          trigger rate vs time (spill structure), per-
                               detector event counts, hits/event, cumulative
                               events (FEU drop-out)
  occupancy_<sub_run>.png      per-channel occupancy per detector with the
                               Dream connector blocks marked, dead/hot flagged
  signal_<sub_run>.png         amplitude spectra, max-sample (latency) and
                               hit-time distributions per detector
  hv_<sub_run>.png             vmon / imon vs time for every HV channel
  commissioning_qa_<sub_run>.json   every number + every verdict
  commissioning_qa_<sub_run>.csv    the per-detector metric table

Run level (<Analysis>/telescope/<run>/25_commissioning_qa/):
  trend_<run>.png              every metric vs sub_run -- and when the DAQ
                               recorded a scan variable in the sub_runs list
                               (e.g. `latency`), vs THAT, which turns this
                               panel into the latency-scan curve
  commissioning_summary_<run>.json  roll-up + the worst verdict per check

Usage:
  SPS_DATA_ROOT=/local/home/banco/P2_data/TB_July2026_H4/runs \
  SPS_ANALYSIS_ROOT=/local/home/banco/P2_data/TB_July2026_H4/analysis \
  SPS_RUN=beam_commissioning_1 python3 25_commissioning_qa.py live
  # options: [--sub-run NAME] [--max-chunks N] [--min-amp ADC] [--no-hv]
"""

import os
import re
import csv
import json
import argparse
import datetime as dt

import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import sps_config as sc
import p2_io as p2io

CH_PER_CONNECTOR = 64          # Dream connector granularity on a FEU
N_CH_PER_FEU = 512

# Verdict thresholds (overridable on the command line).
DEFAULTS = dict(
    share_warn=0.50,           # trigger_share below this vs the best plane
    share_fail=0.05,
    dead_frac_ch=0.10,         # channel counted dead below 10% of the median
    dead_warn=0.10,            # ... and WARN when >10% of channels are dead
    dead_fail=0.30,
    hot_factor=10.0,           # hot above 10x the median
    hot_warn=0.02,
    hot_fail=0.10,
    edge_samples=2,            # max-sample peak this close to an edge -> bad
    sat_warn=0.02,
    sat_fail=0.10,
    gap_frac_warn=0.25,        # longest trigger-free gap / DAQ run time
    livetime_warn=0.80,        # trigger-timestamp span / DAQ run time
    livetime_fail=0.40,
    empty_warn=0.30,           # fraction of triggers with no hit anywhere
    empty_fail=0.70,
    hv_setpoint_tol=0.02,      # |vmon - v0| / v0
    hv_at_warn=0.95,           # fraction of the sub_run at setpoint
    hv_at_fail=0.50,
    imon_warn=2.0,             # uA
    imon_fail=10.0,
)

PASS, WARN, FAIL, NA = 'PASS', 'WARN', 'FAIL', 'N/A'
_RANK = {PASS: 0, NA: 1, WARN: 2, FAIL: 3}
_COLOR = {PASS: '#2e7d32', WARN: '#ef6c00', FAIL: '#c62828', NA: '#757575'}


def worst(*verdicts):
    vs = [v for v in verdicts if v]
    return max(vs, key=lambda v: _RANK[v]) if vs else NA


# --------------------------------------------------------------------------- #
# Detector = one contiguous channel slice of one FEU, from the DAQ wiring
# --------------------------------------------------------------------------- #
class Slice:
    """A detector's (feu, channel-range) footprint plus its connector labels."""

    def __init__(self, det, dream_feus):
        self.name = det.name
        self.det_type = det.det_type
        self.z = det.z
        # {dream connector index -> wiring key}, e.g. {1: 'c_4_bot'} / {5: 'x1'}
        self.conn_label = {}
        feus = set()
        for key, (feu, conn) in ((k, (int(v[0]), int(v[1])))
                                 for k, v in dream_feus.items()):
            feus.add(feu)
            self.conn_label[conn] = key
        if len(feus) != 1:
            # A detector spread over several FEUs cannot be one slice; fall back
            # to the whole of its lowest FEU and say so.
            print(f'  [warn] {self.name}: wiring spans FEUs {sorted(feus)} -- '
                  f'QA slice uses FEU {min(feus)} only')
        self.feu = min(feus)
        self.connectors = sorted(self.conn_label)
        self.ch_lo = (min(self.connectors) - 1) * CH_PER_CONNECTOR
        self.ch_hi = max(self.connectors) * CH_PER_CONNECTOR - 1   # inclusive
        self.n_ch = self.ch_hi - self.ch_lo + 1

    @property
    def tag(self):
        return f'{self.name} (FEU {self.feu} ch {self.ch_lo}-{self.ch_hi})'

    def mask(self, feu_arr, ch_arr):
        return ((feu_arr == self.feu) & (ch_arr >= self.ch_lo) &
                (ch_arr <= self.ch_hi))

    def __repr__(self):
        return f'<Slice {self.tag}>'


def build_slices(cfg):
    return [Slice(d, d.dream_feus) for d in cfg.detectors() if d.dream_feus]


# --------------------------------------------------------------------------- #
# Streamed reduction of one sub_run
# --------------------------------------------------------------------------- #
BRANCHES = ['eventId', 'trigger_timestamp_ns', 'feu', 'channel', 'amplitude',
            'max_sample', 'time', 'saturated']

AMP_BINS = np.linspace(0, 4096, 129)
N_MS_BINS = 64                 # max_sample is a float: sub-sample resolution


def reduce_subrun(cfg, sub_run, slices, n_samples, sample_period,
                  max_chunks=None, min_amp=0.0):
    """One pass over the combined hits -> per-detector accumulators.

    Everything is reduced chunk by chunk (p2_io.iter_hits) so peak memory stays
    at one chunk however long the run is.
    """
    hits_dir = cfg.combined_hits_dir(sub_run)
    acc = {}
    for s in slices:
        acc[s.name] = dict(
            slice=s,
            n_hits=0, n_sat=0,
            occ=np.zeros(s.n_ch, dtype=np.int64),
            occ_amp=np.zeros(s.n_ch, dtype=np.float64),
            amp_hist=np.zeros(len(AMP_BINS) - 1, dtype=np.int64),
            ms_hist=np.zeros(N_MS_BINS, dtype=np.int64),
            t_hist=np.zeros(N_MS_BINS, dtype=np.int64),
            ev_parts=[],
        )
    ms_edges = np.linspace(0, max(n_samples, 1), N_MS_BINS + 1)
    t_edges = np.linspace(0, max(n_samples, 1) * sample_period, N_MS_BINS + 1)

    n_chunks = 0
    for df in p2io.iter_hits(hits_dir, BRANCHES, progress=True,
                             min_amp=min_amp):
        n_chunks += 1
        if not len(df):
            continue
        feu = df['feu'].to_numpy()
        ch = df['channel'].to_numpy()
        for s in slices:
            m = s.mask(feu, ch)
            if not m.any():
                continue
            a = acc[s.name]
            sub = df[m]
            rel = ch[m] - s.ch_lo
            amp = sub['amplitude'].to_numpy(np.float64)
            a['n_hits'] += int(m.sum())
            a['n_sat'] += int(sub['saturated'].to_numpy().sum())
            a['occ'] += np.bincount(rel, minlength=s.n_ch)
            a['occ_amp'] += np.bincount(rel, weights=amp, minlength=s.n_ch)
            a['amp_hist'] += np.histogram(amp, bins=AMP_BINS)[0]
            a['ms_hist'] += np.histogram(sub['max_sample'].to_numpy(),
                                         bins=ms_edges)[0]
            a['t_hist'] += np.histogram(sub['time'].to_numpy(),
                                        bins=t_edges)[0]
            g = sub.groupby('eventId')
            a['ev_parts'].append(pd.DataFrame({
                'n_hits': g.size(),
                'ts': g['trigger_timestamp_ns'].min()}))
        if max_chunks and n_chunks >= max_chunks:
            print(f'    [--max-chunks {max_chunks}] stopping early', flush=True)
            break

    for a in acc.values():
        parts = a.pop('ev_parts')
        if parts:
            ev = pd.concat(parts)
            ev = ev.groupby(level=0).agg(n_hits=('n_hits', 'sum'),
                                         ts=('ts', 'min')).sort_index()
        else:
            ev = pd.DataFrame(columns=['n_hits', 'ts'])
        a['events'] = ev
    return acc, ms_edges, t_edges, n_chunks


# --------------------------------------------------------------------------- #
# Per-detector metrics + verdicts
# --------------------------------------------------------------------------- #
def detector_metrics(a, thr, n_triggers, ms_edges, t_edges, n_samples,
                     duration_s):
    s = a['slice']
    ev = a['events']
    occ = a['occ']
    m = {'detector': s.name, 'det_type': s.det_type, 'feu': s.feu,
         'z_mm': s.z, 'ch_lo': s.ch_lo, 'ch_hi': s.ch_hi,
         'n_hits': int(a['n_hits']), 'n_events': int(len(ev))}

    m['trigger_share'] = float(len(ev) / n_triggers) if n_triggers else None
    m['hit_rate_hz'] = float(a['n_hits'] / duration_s) if duration_s else None
    m['event_rate_hz'] = float(len(ev) / duration_s) if duration_s else None
    m['mean_hits_per_event'] = float(ev['n_hits'].mean()) if len(ev) else 0.0
    m['median_hits_per_event'] = float(ev['n_hits'].median()) if len(ev) else 0.0
    m['sat_frac'] = float(a['n_sat'] / a['n_hits']) if a['n_hits'] else None

    # --- occupancy: dead / quiet / hot channels, per-connector liveness ----- #
    # A beam spot makes occupancy intrinsically non-uniform, so the three are
    # deliberately defined differently:
    #   dead  -- LITERALLY zero hits: no beam-spot argument can explain it, so
    #            this is the one that means "channel not reading out".
    #   quiet -- far below the median: informational (usually just outside the
    #            beam spot on a plane the beam clips).
    #   hot   -- far above the 95th percentile, NOT the median: inside the spot
    #            a channel can legitimately sit 10x over the median, so the
    #            median would flag the whole beam spot as noise.
    med = float(np.median(occ[occ > 0])) if (occ > 0).any() else 0.0
    p95 = float(np.percentile(occ, 95)) if (occ > 0).any() else 0.0
    m['median_channel_occupancy'] = med
    m['p95_channel_occupancy'] = p95
    dead = occ == 0
    quiet = (occ < thr['dead_frac_ch'] * med) & ~dead if med else \
        np.zeros_like(occ, bool)
    hot = occ > thr['hot_factor'] * p95 if p95 else np.zeros_like(occ, bool)
    m['n_dead_channels'] = int(dead.sum())
    m['dead_frac'] = float(dead.mean())
    m['n_quiet_channels'] = int(quiet.sum())
    m['quiet_frac'] = float(quiet.mean())
    m['n_hot_channels'] = int(hot.sum())
    m['hot_frac'] = float(hot.mean())
    m['hot_channels'] = (np.flatnonzero(hot) + s.ch_lo).tolist()[:32]

    conn_hits, dead_conns = {}, []
    for c in s.connectors:
        lo = (c - 1) * CH_PER_CONNECTOR - s.ch_lo
        blk = occ[lo:lo + CH_PER_CONNECTOR]
        label = s.conn_label.get(c, str(c))
        conn_hits[f'{c}:{label}'] = int(blk.sum())
        if (blk == 0).mean() > 0.9:      # a whole 64-ch block silent = cabling
            dead_conns.append(f'{c}:{label}')
    m['connector_hits'] = conn_hits
    m['dead_connectors'] = dead_conns

    # --- amplitude ---------------------------------------------------------- #
    ah = a['amp_hist']
    centres = 0.5 * (AMP_BINS[:-1] + AMP_BINS[1:])
    if ah.sum():
        m['amp_peak_adc'] = float(centres[int(np.argmax(ah))])
        cum = np.cumsum(ah) / ah.sum()
        m['amp_median_adc'] = float(np.interp(0.5, cum, centres))
        m['amp_p90_adc'] = float(np.interp(0.9, cum, centres))
    else:
        m['amp_peak_adc'] = m['amp_median_adc'] = m['amp_p90_adc'] = None

    # --- latency: where in the sampling window does the pulse peak? --------- #
    mh = a['ms_hist']
    ms_c = 0.5 * (ms_edges[:-1] + ms_edges[1:])
    if mh.sum():
        m['max_sample_peak'] = float(ms_c[int(np.argmax(mh))])
        m['max_sample_mean'] = float((ms_c * mh).sum() / mh.sum())
        cum = np.cumsum(mh) / mh.sum()
        m['max_sample_median'] = float(np.interp(0.5, cum, ms_c))
        edge = thr['edge_samples']
        m['frac_first_edge'] = float(mh[ms_c < edge].sum() / mh.sum())
        m['frac_last_edge'] = float(
            mh[ms_c > n_samples - edge].sum() / mh.sum())
    else:
        for k in ('max_sample_peak', 'max_sample_mean', 'max_sample_median',
                  'frac_first_edge', 'frac_last_edge'):
            m[k] = None
    th = a['t_hist']
    t_c = 0.5 * (t_edges[:-1] + t_edges[1:])
    m['hit_time_peak_ns'] = float(t_c[int(np.argmax(th))]) if th.sum() else None

    # --- rate stability over the sub_run ------------------------------------ #
    # NOT a rate-vs-time-bin metric: SPS slow extraction means most 1 s bins are
    # legitimately empty between spills, so any "min bin / median bin" measure
    # flags a perfectly healthy run. What a dropout actually looks like is a
    # GAP much longer than the SPS cycle, so that is what is measured -- plus
    # where the detector's first and last trigger sit inside the sub_run.
    m['max_gap_s'] = m['first_hit_s'] = m['last_hit_s'] = m['gap_frac'] = None
    if len(ev) > 50:
        ts = np.sort(ev['ts'].to_numpy(np.float64)) / 1e9
        t = ts - ts[0]
        span = float(t[-1])
        if 0 < span < 24 * 3600:
            m['max_gap_s'] = float(np.max(np.diff(t))) if len(t) > 1 else 0.0
            m['first_hit_s'] = 0.0
            m['last_hit_s'] = span
            if duration_s:
                m['gap_frac'] = float(m['max_gap_s'] / duration_s)
        m['ts_span_s'] = span
    return m


def detector_verdicts(m, thr, best_share):
    v = {}
    v['readout'] = FAIL if not m['n_hits'] else PASS

    sh = m['trigger_share']
    if sh is None or not best_share:
        v['trigger_share'] = NA
    else:
        r = sh / best_share
        v['trigger_share'] = (FAIL if r < thr['share_fail'] else
                              WARN if r < thr['share_warn'] else PASS)

    v['connectors'] = FAIL if m['dead_connectors'] else PASS
    v['dead_channels'] = (FAIL if m['dead_frac'] > thr['dead_fail'] else
                          WARN if m['dead_frac'] > thr['dead_warn'] else PASS)
    if m['det_type'] != 'P2' and v['dead_channels'] == FAIL:
        # A P2 uses all 64 channels of every connector it occupies, so a dead
        # channel there is a fault. A strip detector (the uRWELL references)
        # generally has fewer strips per view than the connector has channels,
        # so its unconnected channels are silent BY CONSTRUCTION -- never let
        # that FAIL the run; it is reported, not judged.
        v['dead_channels'] = WARN
        m['dead_note'] = ('non-P2 detector: silent channels are most likely '
                          'unconnected strips, not a readout fault')
    v['hot_channels'] = (FAIL if m['hot_frac'] > thr['hot_fail'] else
                         WARN if m['hot_frac'] > thr['hot_warn'] else PASS)

    if m['max_sample_peak'] is None:
        v['latency'] = NA
    else:
        edge = max(m['frac_first_edge'], m['frac_last_edge'])
        v['latency'] = (FAIL if edge > 0.30 else WARN if edge > 0.10 else PASS)

    sf = m['sat_frac']
    v['saturation'] = (NA if sf is None else
                       FAIL if sf > thr['sat_fail'] else
                       WARN if sf > thr['sat_warn'] else PASS)

    gf = m['gap_frac']
    v['rate_stability'] = (NA if gf is None else
                           FAIL if gf > 2 * thr['gap_frac_warn'] else
                           WARN if gf > thr['gap_frac_warn'] else PASS)
    return v


# --------------------------------------------------------------------------- #
# HV monitor
# --------------------------------------------------------------------------- #
def read_hv(path):
    """hv_monitor.csv -> (DataFrame indexed by elapsed seconds, {chan: {...}})."""
    if not os.path.isfile(path):
        return None, {}
    df = pd.read_csv(path)
    if not len(df):
        return None, {}
    ts = pd.to_datetime(df['timestamp'])
    df['t_s'] = (ts - ts.iloc[0]).dt.total_seconds()
    chans = sorted({c.rsplit(' ', 1)[0] for c in df.columns if ' vmon' in c},
                   key=lambda c: tuple(int(x) for x in c.split(':')))
    return df, {c: dict(v0=f'{c} v0', vmon=f'{c} vmon', imon=f'{c} imon')
                for c in chans}


def hv_metrics(df, chans, thr):
    out = {}
    for c, cols in chans.items():
        v0 = df[cols['v0']].to_numpy(float)
        vm = df[cols['vmon']].to_numpy(float)
        im = df[cols['imon']].to_numpy(float)
        with np.errstate(divide='ignore', invalid='ignore'):
            rel = np.abs(vm - v0) / np.where(v0 == 0, np.nan, v0)
        at = np.isfinite(rel) & (rel <= thr['hv_setpoint_tol'])
        out[c] = dict(v0=float(v0[-1]),
                      vmon_mean=float(np.mean(vm)),
                      vmon_last=float(vm[-1]),
                      at_setpoint_frac=float(at.mean()),
                      imon_mean=float(np.mean(im)),
                      imon_max=float(np.max(im)),
                      n_imon_over_warn=int((im > thr['imon_warn']).sum()))
    return out


def hv_verdicts(hvm, thr):
    if not hvm:
        return {'hv_setpoint': NA, 'hv_current': NA}, None, None
    at = min(v['at_setpoint_frac'] for v in hvm.values())
    imx = max(v['imon_max'] for v in hvm.values())
    return ({'hv_setpoint': (FAIL if at < thr['hv_at_fail'] else
                             WARN if at < thr['hv_at_warn'] else PASS),
             'hv_current': (FAIL if imx > thr['imon_fail'] else
                            WARN if imx > thr['imon_warn'] else PASS)},
            at, imx)


# --------------------------------------------------------------------------- #
# Plots
# --------------------------------------------------------------------------- #
def _det_colors(slices):
    cmap = plt.get_cmap('tab10')
    return {s.name: cmap(i % 10) for i, s in enumerate(slices)}


def plot_rates(acc, slices, colors, sub_run, duration_s, out_png, title_extra):
    fig, axes = plt.subplots(2, 2, figsize=(14, 8.5))
    ax = axes[0, 0]
    # Trigger rate vs time from the union of all detectors' event timestamps.
    ref = max(acc.values(), key=lambda a: len(a['events']))
    ev = ref['events']
    if len(ev) > 10:
        t = (ev['ts'].to_numpy(np.float64) - float(ev['ts'].min())) / 1e9
        span = max(t.max(), 1e-3)
        for bw, style in ((0.2, dict(lw=0.8, alpha=0.55, color='#1565c0')),
                          (1.0, dict(lw=1.8, color='#0d47a1'))):
            nb = max(2, int(span / bw))
            cnt, edges = np.histogram(t, bins=nb, range=(0, span))
            ax.step(edges[:-1], cnt / (edges[1] - edges[0]), where='post',
                    label=f'{bw:g} s bins', **style)
        if duration_s:
            ax.axvline(duration_s, color='k', ls='--', lw=1,
                       label=f'DAQ run time {duration_s:.0f} s')
        ax.set_xlabel('time in sub_run [s]  (FEU trigger timestamp)')
        ax.set_ylabel('trigger rate [Hz]')
        ax.legend(fontsize=8)
    ax.set_title(f'trigger rate vs time -- SPS spill structure\n'
                 f'reference {ref["slice"].name}', fontsize=9)
    ax.grid(alpha=.3)

    ax = axes[0, 1]
    names = [s.name for s in slices]
    n_ev = [len(acc[n]['events']) for n in names]
    n_hit = [acc[n]['n_hits'] for n in names]
    x = np.arange(len(names))
    ax.bar(x - 0.2, n_ev, 0.4, label='events with >=1 hit',
           color=[colors[n] for n in names])
    ax.bar(x + 0.2, n_hit, 0.4, label='hits', alpha=.45,
           color=[colors[n] for n in names])
    ax.set_xticks(x)
    ax.set_xticklabels([f'{n}\nFEU {acc[n]["slice"].feu}' for n in names],
                       fontsize=7)
    ax.set_yscale('log')
    ax.set_ylabel('count')
    ax.legend(fontsize=8)
    ax.set_title('per-detector statistics', fontsize=9)
    ax.grid(alpha=.3, axis='y')

    ax = axes[1, 0]
    for n in names:
        ev = acc[n]['events']
        if not len(ev):
            continue
        v = ev['n_hits'].to_numpy()
        ax.hist(v, bins=np.arange(0, min(v.max(), 200) + 2) - 0.5,
                histtype='step', lw=1.5, label=n, color=colors[n], log=True)
    ax.set_xlabel('hits per triggered event')
    ax.set_ylabel('events')
    ax.legend(fontsize=8)
    ax.set_title('hit multiplicity', fontsize=9)
    ax.grid(alpha=.3)

    ax = axes[1, 1]
    for n in names:
        ev = acc[n]['events']
        if len(ev) < 10:
            continue
        t = (ev['ts'].to_numpy(np.float64) - float(ev['ts'].min())) / 1e9
        ax.plot(np.sort(t), np.arange(1, len(t) + 1), lw=1.5, label=n,
                color=colors[n])
    ax.set_xlabel('time in sub_run [s]')
    ax.set_ylabel('cumulative events')
    ax.legend(fontsize=8)
    ax.set_title('cumulative events -- a flat stretch = that FEU dropped out',
                 fontsize=9)
    ax.grid(alpha=.3)

    fig.suptitle(f'{sub_run} -- rates & trigger stream   {title_extra}',
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


def plot_occupancy(acc, slices, colors, metrics, sub_run, out_png, thr):
    n = len(slices)
    fig, axes = plt.subplots(n, 1, figsize=(14, 2.3 * n + 1), squeeze=False)
    for ax, s in zip(axes[:, 0], slices):
        a, m = acc[s.name], metrics[s.name]
        chans = np.arange(s.ch_lo, s.ch_hi + 1)
        ax.step(chans, a['occ'], where='mid', lw=0.8, color=colors[s.name])
        med = m['median_channel_occupancy']
        if med:
            ax.axhline(med, ls='--', lw=0.8, color='k', alpha=.5)
            ax.axhline(thr['dead_frac_ch'] * med, ls=':', lw=0.8, color='r',
                       alpha=.6)
            dead = a['occ'] < thr['dead_frac_ch'] * med
            if dead.any():
                ax.plot(chans[dead], np.full(dead.sum(), 0.5), '|', ms=6,
                        color='r', alpha=.7)
        for c in s.connectors:
            x0 = (c - 1) * CH_PER_CONNECTOR
            ax.axvline(x0, color='gray', lw=0.6, alpha=.5)
            ax.text(x0 + CH_PER_CONNECTOR / 2, ax.get_ylim()[1],
                    s.conn_label.get(c, str(c)), fontsize=6, ha='center',
                    va='top', color='gray')
        ax.set_yscale('log')
        ax.set_ylabel('hits')
        ax.set_xlim(s.ch_lo - 2, s.ch_hi + 2)
        ax.set_title(
            f'{s.tag}   dead {m["n_dead_channels"]}/{s.n_ch} '
            f'({100 * m["dead_frac"]:.1f}%)   hot {m["n_hot_channels"]}   '
            f'dead connectors: {", ".join(m["dead_connectors"]) or "none"}',
            fontsize=8)
        ax.grid(alpha=.25)
    axes[-1, 0].set_xlabel('FEU channel')
    fig.suptitle(f'{sub_run} -- channel occupancy (dashed = median, '
                 f'dotted = dead threshold)', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


def plot_signal(acc, slices, colors, metrics, sub_run, ms_edges, t_edges,
                n_samples, sample_period, latency, out_png):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
    ms_c = 0.5 * (ms_edges[:-1] + ms_edges[1:])
    t_c = 0.5 * (t_edges[:-1] + t_edges[1:])
    amp_c = 0.5 * (AMP_BINS[:-1] + AMP_BINS[1:])

    ax = axes[0]
    for s in slices:
        h = acc[s.name]['amp_hist']
        if h.sum():
            med = metrics[s.name]['amp_median_adc']
            ax.step(amp_c, h, where='mid', lw=1.4, color=colors[s.name],
                    label=f'{s.name}  med {med:.0f}' if med is not None
                    else s.name)
    ax.set_xlabel('hit amplitude [ADC]')
    ax.set_ylabel('hits')
    ax.set_yscale('log')
    ax.legend(fontsize=7)
    ax.set_title('amplitude spectra', fontsize=9)
    ax.grid(alpha=.3)

    ax = axes[1]
    for s in slices:
        h = acc[s.name]['ms_hist']
        if h.sum():
            ax.step(ms_c, h / h.sum(), where='mid', lw=1.4,
                    color=colors[s.name], label=s.name)
    ax.axvspan(0, 2, color='r', alpha=.08)
    ax.axvspan(n_samples - 2, n_samples, color='r', alpha=.08)
    ax.set_xlim(0, n_samples)
    ax.set_xlabel('sample of waveform maximum')
    ax.set_ylabel('fraction of hits')
    ax.legend(fontsize=7)
    ax.set_title(f'LATENCY check -- window = {n_samples} x {sample_period} ns'
                 + (f', Dream latency {latency}' if latency is not None else ''),
                 fontsize=9)
    ax.grid(alpha=.3)

    ax = axes[2]
    for s in slices:
        h = acc[s.name]['t_hist']
        if h.sum():
            ax.step(t_c, h / h.sum(), where='mid', lw=1.4,
                    color=colors[s.name], label=s.name)
    ax.set_xlabel('hit time [ns]')
    ax.set_ylabel('fraction of hits')
    ax.legend(fontsize=7)
    ax.set_title('hit-time distribution', fontsize=9)
    ax.grid(alpha=.3)

    fig.suptitle(f'{sub_run} -- signal shape & timing', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


def plot_hv(df, chans, hvm, sub_run, out_png):
    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    t = df['t_s'].to_numpy()
    cmap = plt.get_cmap('tab20')
    for i, (c, cols) in enumerate(chans.items()):
        col = cmap(i % 20)
        axes[0].plot(t, df[cols['vmon']], lw=1.3, color=col,
                     label=f'{c} (set {hvm[c]["v0"]:.0f} V, '
                           f'{100 * hvm[c]["at_setpoint_frac"]:.0f}% at setpoint)')
        axes[0].plot(t, df[cols['v0']], lw=0.7, ls='--', color=col, alpha=.5)
        axes[1].plot(t, df[cols['imon']], lw=1.1, color=col, label=c)
    axes[0].set_ylabel('vmon [V]  (dashed = setpoint)')
    axes[0].legend(fontsize=6, ncol=2)
    axes[0].grid(alpha=.3)
    axes[1].set_ylabel('imon [uA]')
    axes[1].set_xlabel('time in sub_run [s]')
    axes[1].set_yscale('symlog', linthresh=0.1)
    axes[1].legend(fontsize=6, ncol=4)
    axes[1].grid(alpha=.3)
    fig.suptitle(f'{sub_run} -- HV during the sub_run', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


def plot_trend(results, slices, colors, run, scan_var, out_png):
    """Every metric vs sub_run -- vs the DAQ's scan variable when there is one
    (e.g. `latency`), which makes this the latency-scan summary."""
    names = [s.name for s in slices]
    subs = [r['sub_run'] for r in results]
    xs = [r.get('scan_value') for r in results]
    if scan_var and all(v is not None for v in xs):
        x = np.array(xs, float)
        xlabel = scan_var
    else:
        x = np.arange(len(subs), dtype=float)
        xlabel = 'sub_run'

    panels = [('trigger_share', 'fraction of triggers with a hit', True),
              ('event_rate_hz', 'events with a hit [Hz]', False),
              ('max_sample_peak', 'peak sample of waveform max', False),
              ('amp_median_adc', 'median hit amplitude [ADC]', False),
              ('mean_hits_per_event', 'mean hits / triggered event', False),
              ('dead_frac', 'fraction of dead channels', True)]
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    for ax, (key, ylab, zero_base) in zip(axes.ravel(), panels):
        for n in names:
            y = [r['detectors'].get(n, {}).get(key) for r in results]
            ok = [i for i, v in enumerate(y) if v is not None]
            if not ok:
                continue
            ax.plot(x[ok], [y[i] for i in ok], 'o-', lw=1.5, ms=4,
                    color=colors[n], label=n)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylab, fontsize=8)
        if zero_base:
            ax.set_ylim(bottom=0)
        if xlabel == 'sub_run':
            ax.set_xticks(x)
            ax.set_xticklabels(subs, rotation=45, ha='right', fontsize=6)
        ax.grid(alpha=.3)
    axes[0, 0].legend(fontsize=7)
    ttl = f'{run} -- commissioning trend'
    if xlabel != 'sub_run':
        ttl += (f'   (scan variable: {scan_var} -- the peak of '
                f'"fraction of triggers with a hit" is the optimum)')
    fig.suptitle(ttl, fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Reporting
# --------------------------------------------------------------------------- #
CHECKS = ['readout', 'trigger_share', 'connectors', 'dead_channels',
          'hot_channels', 'latency', 'saturation', 'rate_stability']


def print_table(res):
    ef = res['empty_trigger_frac']
    print(f'\n  sub_run {res["sub_run"]}   triggers {res["n_triggers"]} '
          f'({res["n_triggers_with_hits"]} with hits, '
          f'{100 * ef:.0f}% empty [{res["run_checks"]["empty_triggers"]}])')
    print(f'    DAQ run time {res["duration_s"]} s, trigger-timestamp livetime '
          f'{res["livetime_s"]} s ({res["livetime_frac"]} '
          f'[{res["run_checks"]["livetime"]}])   '
          f'rate {res["trigger_rate_hz"]} Hz '
          f'({res["trigger_rate_live_hz"]} Hz live)')
    hdr = (f'    {"detector":<18}{"FEU":>4}{"events":>9}{"share":>7}'
           f'{"hits/ev":>9}{"dead":>6}{"hot":>5}{"peak smp":>9}'
           f'{"amp med":>9}   verdict')
    print(hdr)
    print('    ' + '-' * (len(hdr) - 4))
    for n, m in res['detectors'].items():
        v = res['verdicts'][n]
        bad = [k for k in CHECKS if v.get(k) in (WARN, FAIL)]
        print(f'    {n:<18}{m["feu"]:>4}{m["n_events"]:>9}'
              f'{(m["trigger_share"] or 0):>7.2f}'
              f'{m["mean_hits_per_event"]:>9.1f}'
              f'{m["n_dead_channels"]:>6}{m["n_hot_channels"]:>5}'
              f'{(m["max_sample_peak"] or 0):>9.1f}'
              f'{(m["amp_median_adc"] or 0):>9.0f}   '
              f'{worst(*v.values())}'
              + (f'  <- {", ".join(bad)}' if bad else ''))
    hv = res['hv']
    if hv.get('channels'):
        print(f'    HV: min at-setpoint {100 * hv["min_at_setpoint_frac"]:.0f}%'
              f'   max imon {hv["max_imon_uA"]:.2f} uA'
              f'   [{hv["verdicts"]["hv_setpoint"]}/'
              f'{hv["verdicts"]["hv_current"]}]')
    print(f'    ==> {res["sub_run"]} verdict: {res["verdict"]}')


def write_csv(res, path):
    rows = []
    for n, m in res['detectors'].items():
        row = {k: v for k, v in m.items()
               if not isinstance(v, (dict, list))}
        row.update({f'v_{k}': v for k, v in res['verdicts'][n].items()})
        row['sub_run'] = res['sub_run']
        rows.append(row)
    if not rows:
        return
    keys = list(rows[0])
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def _jsonable(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(type(o))


# --------------------------------------------------------------------------- #
def analyse_subrun(cfg, sub_run, slices, colors, thr, args):
    daq = cfg.daq_info()
    n_samples = int(daq.get('n_samples_per_waveform') or 16)
    sample_period = float(daq.get('sample_period') or 60)
    meta = cfg.subrun_meta(sub_run)
    latency = meta.get('latency', daq.get('latency'))

    print(f'\n== {sub_run}', flush=True)
    acc, ms_edges, t_edges, n_chunks = reduce_subrun(
        cfg, sub_run, slices, n_samples, sample_period,
        max_chunks=args.max_chunks, min_amp=args.min_amp)

    # The trigger stream. With zero suppression a trigger that fired no channel
    # anywhere leaves NO row in the hits file, so the number of eventIds that
    # carry hits undercounts the triggers. The DAQ's eventId is a contiguous
    # per-run trigger counter, so its SPAN is the true trigger count and the
    # difference is the fraction of triggers that were completely empty --
    # itself a headline commissioning number (a noisy/accidental trigger).
    all_ids = set()
    for a in acc.values():
        all_ids.update(a['events'].index.tolist())
    n_with_hits = len(all_ids)
    n_triggers = (max(all_ids) - min(all_ids) + 1) if all_ids else 0
    empty_frac = (1.0 - n_with_hits / n_triggers) if n_triggers else None

    duration_s, start_ts = cfg.run_time(sub_run)
    # Live time = the span the FEU trigger timestamps actually cover. A live
    # time well short of the DAQ run time means the trigger stopped (beam off,
    # spill structure at the end of the run, or the DAQ stalled).
    spans = [((float(a['events']['ts'].max()) -
               float(a['events']['ts'].min())) / 1e9)
             for a in acc.values() if len(a['events']) > 1]
    livetime_s = float(max(spans)) if spans else None
    if not duration_s:
        duration_s = livetime_s

    metrics, verdicts = {}, {}
    shares = [len(a['events']) / n_triggers for a in acc.values()] \
        if n_triggers else [0]
    best_share = max(shares) if shares else 0
    for s in slices:
        m = detector_metrics(acc[s.name], thr, n_triggers, ms_edges, t_edges,
                             n_samples, duration_s)
        metrics[s.name] = m
        verdicts[s.name] = detector_verdicts(m, thr, best_share)

    # HV
    hv = {'channels': {}, 'verdicts': {'hv_setpoint': NA, 'hv_current': NA}}
    hv_df, hv_chans = (None, {}) if args.no_hv else \
        read_hv(cfg.hv_monitor_csv(sub_run))
    if hv_df is not None and hv_chans:
        hvm = hv_metrics(hv_df, hv_chans, thr)
        hvv, at, imx = hv_verdicts(hvm, thr)
        hv = {'channels': hvm, 'verdicts': hvv,
              'min_at_setpoint_frac': at, 'max_imon_uA': imx,
              'n_samples': int(len(hv_df))}

    res = {
        'run': cfg.RUN, 'sub_run': sub_run,
        'analysed_at': dt.datetime.now().isoformat(timespec='seconds'),
        'n_chunks': n_chunks,
        'n_triggers': n_triggers,
        'n_triggers_with_hits': n_with_hits,
        'empty_trigger_frac': empty_frac,
        'duration_s': round(duration_s, 2) if duration_s else None,
        'livetime_s': round(livetime_s, 2) if livetime_s else None,
        'livetime_frac': (round(livetime_s / duration_s, 3)
                          if livetime_s and duration_s else None),
        'trigger_rate_hz': round(n_triggers / duration_s, 2)
                           if duration_s else None,
        'trigger_rate_live_hz': round(n_triggers / livetime_s, 2)
                                if livetime_s else None,
        'start_unix_ts': start_ts,
        'n_samples_per_waveform': n_samples, 'sample_period_ns': sample_period,
        'latency': latency,
        'scan_value': meta.get('latency'),
        'zero_suppress': daq.get('zero_suppress'),
        'detectors': metrics, 'verdicts': verdicts, 'hv': hv,
    }
    lf, ef = res['livetime_frac'], empty_frac
    res['run_checks'] = {
        'livetime': (NA if lf is None else
                     FAIL if lf < thr['livetime_fail'] else
                     WARN if lf < thr['livetime_warn'] else PASS),
        'empty_triggers': (NA if ef is None else
                           FAIL if ef > thr['empty_fail'] else
                           WARN if ef > thr['empty_warn'] else PASS),
    }
    res['verdict'] = worst(*[worst(*v.values()) for v in verdicts.values()],
                           *hv['verdicts'].values(),
                           *res['run_checks'].values())

    out = cfg.out_dir(sc.TELESCOPE_TAG, sub_run, '25_commissioning_qa')
    title = (f'{n_triggers} triggers, {res["trigger_rate_hz"]} Hz, '
             f'{res["duration_s"]} s')
    plot_rates(acc, slices, colors, sub_run, duration_s,
               os.path.join(out, f'rates_{sub_run}.png'), title)
    plot_occupancy(acc, slices, colors, metrics, sub_run,
                   os.path.join(out, f'occupancy_{sub_run}.png'), thr)
    plot_signal(acc, slices, colors, metrics, sub_run, ms_edges, t_edges,
                n_samples, sample_period, latency,
                os.path.join(out, f'signal_{sub_run}.png'))
    if hv_df is not None and hv_chans:
        plot_hv(hv_df, hv_chans, hv['channels'], sub_run,
                os.path.join(out, f'hv_{sub_run}.png'))
    with open(os.path.join(out, f'commissioning_qa_{sub_run}.json'), 'w') as fh:
        json.dump(res, fh, indent=2, default=_jsonable)
    write_csv(res, os.path.join(out, f'commissioning_qa_{sub_run}.csv'))
    print_table(res)
    print(f'    -> {out}')
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('run_key', nargs='?', default=None)
    ap.add_argument('--sub-run', default=None,
                    help='analyse only this sub_run (default: all on disk)')
    ap.add_argument('--max-chunks', type=int, default=None,
                    help='stop after N combined-hits files (quick look)')
    ap.add_argument('--min-amp', type=float, default=0.0,
                    help='drop hits below this amplitude [ADC]')
    ap.add_argument('--no-hv', action='store_true')
    for k, v in DEFAULTS.items():
        ap.add_argument(f'--{k.replace("_", "-")}', type=type(v), default=v)
    args = ap.parse_args()
    thr = {k: getattr(args, k) for k in DEFAULTS}

    cfg = sc.get_config(args.run_key)
    print(cfg)
    slices = build_slices(cfg)
    if not slices:
        raise SystemExit('No detector wiring in the run_config -- nothing to do')
    for s in slices:
        print('  ', s.tag, '  connectors',
              ', '.join(f'{c}:{s.conn_label[c]}' for c in s.connectors))
    colors = _det_colors(slices)

    subs = [args.sub_run] if args.sub_run else cfg.find_subruns()
    if not subs:
        raise SystemExit(f'No sub_run with combined-hits ROOT under '
                         f'{cfg.run_dir} (still decoding?)')
    print(f'  sub_runs: {subs}')

    results = [analyse_subrun(cfg, s, slices, colors, thr, args) for s in subs]

    # -- run-level roll-up --------------------------------------------------- #
    scan_var = 'latency' if any(r['scan_value'] is not None for r in results) \
        else None
    out = os.path.join(cfg.ANALYSIS_ROOT, sc.TELESCOPE_TAG, cfg.RUN,
                       '25_commissioning_qa')
    os.makedirs(out, exist_ok=True)
    plot_trend(results, slices, colors, cfg.RUN, scan_var,
               os.path.join(out, f'trend_{cfg.RUN}.png'))

    per_check = {}
    for c in CHECKS:
        per_check[c] = worst(*[r['verdicts'][n].get(c)
                               for r in results for n in r['verdicts']])
    for c in ('hv_setpoint', 'hv_current'):
        per_check[c] = worst(*[r['hv']['verdicts'].get(c) for r in results])
    for c in ('livetime', 'empty_triggers'):
        per_check[c] = worst(*[r['run_checks'].get(c) for r in results])
    summary = {
        'run': cfg.RUN,
        'analysed_at': dt.datetime.now().isoformat(timespec='seconds'),
        'sub_runs': [r['sub_run'] for r in results],
        'detectors': [s.tag for s in slices],
        'scan_variable': scan_var,
        'per_check_verdict': per_check,
        'verdict': worst(*[r['verdict'] for r in results]),
        'per_sub_run': {r['sub_run']: {
            'verdict': r['verdict'], 'n_triggers': r['n_triggers'],
            'trigger_rate_hz': r['trigger_rate_hz'],
            'empty_trigger_frac': r['empty_trigger_frac'],
            'livetime_frac': r['livetime_frac'],
            'scan_value': r['scan_value'],
            'trigger_share': {n: m['trigger_share']
                              for n, m in r['detectors'].items()},
            'max_sample_peak': {n: m['max_sample_peak']
                                for n, m in r['detectors'].items()},
        } for r in results},
    }
    with open(os.path.join(out, f'commissioning_summary_{cfg.RUN}.json'),
              'w') as fh:
        json.dump(summary, fh, indent=2, default=_jsonable)

    print(f'\n{"=" * 78}\n  RUN {cfg.RUN}: {summary["verdict"]}')
    for c, v in per_check.items():
        print(f'    {c:<16} {v}')
    if scan_var:
        print(f'\n  scan variable "{scan_var}" -- per sub_run:')
        for r in results:
            shares = ' '.join(
                f'{n.split("_")[-1]}={m["trigger_share"]:.2f}'
                for n, m in r['detectors'].items()
                if m['trigger_share'] is not None)
            print(f'    {r["sub_run"]:<18} {scan_var}={r["scan_value"]}  '
                  f'rate={r["trigger_rate_hz"]} Hz  share: {shares}')
    print(f'  -> {out}')


if __name__ == '__main__':
    main()
