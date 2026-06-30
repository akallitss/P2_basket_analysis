#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p2_channel_qa.py

Raw channel-level QA for a P2 (BASKET) detector under test in the cosmic
bench, built from the *_feu-combined_hits.root files. No spatial mapping is
applied here -- everything is per-FEU / per-channel, which is what you want for
a first sanity look at a freshly bulked detector.

The FEUs (and their instrumented connectors) that belong to the P2 detector are
read straight out of the run's run_config.json, so the same script works for any
bench run without editing FEU lists by hand. A DREAM FEU connector maps to a
contiguous block of 64 channels:  connector c -> channels [(c-1)*64, c*64).

Products (written to <out_dir>/):
  Overall P2 summary
    p2_summary.png            pulse-height dist, hit multiplicity,
                              saturated hits per FEU, hit rate vs time
  Per FEU
    p2_feu<N>_channels.png    occupancy (hits vs channel) +
                              pulse-height per channel (2D + mean profile)

Run headless (Agg); everything is saved to disk.

Usage
  python3 p2_channel_qa.py \
      --hits-dir /path/to/.../combined_hits_root \
      --run-config /path/to/.../run_config.json

If --run-config is omitted the script looks for run_config.json two levels above
the hits dir (the usual bench layout: <run>/<sub_run>/combined_hits_root/).
"""

import matplotlib
matplotlib.use('Agg')  # headless: save figures, never block on a window

import os
import glob
import json
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import uproot

# Branches we actually need -- read only these so the (huge, noisy) reference
# FEUs in the same files don't blow up memory before we filter them out.
_BRANCHES = ['eventId', 'trigger_timestamp_ns', 'channel', 'amplitude',
             'saturated', 'feu', 'time']
# Minimal branches for the (large) reference detector: only what the timing
# coincidence needs.
_REF_BRANCHES = ['eventId', 'feu', 'time']

CH_PER_CONNECTOR = 64      # DREAM FEU connector -> 64 contiguous channels
N_CH_PER_FEU = 512         # 8 connectors per FEU
AMP_LABEL = 'Pulse height (amplitude) [ADC]'


# ---------------------------------------------------------------------------
# Configuration from run_config.json
# ---------------------------------------------------------------------------

def parse_p2_feus(run_config_path, det_type='P2', det_name=None):
    """
    From run_config.json, return for the P2 detector:
      feu_connectors : {feu_number: sorted list of instrumented connector ids}
      feu_channels   : {feu_number: np.array of instrumented channel numbers}
      det_name       : the detector's name as written in the config

    The detector is selected by det_type (default 'P2'), or by det_name if given.
    """
    with open(run_config_path) as fh:
        cfg = json.load(fh)

    detectors = cfg.get('detectors', [])
    chosen = None
    for d in detectors:
        if det_name is not None and d.get('name') == det_name:
            chosen = d
            break
        if det_name is None and d.get('det_type') == det_type:
            chosen = d
            break
    if chosen is None:
        raise ValueError(
            f'No detector with '
            f'{"name="+det_name if det_name else "det_type="+det_type} '
            f'in {run_config_path}. Available: '
            f'{[(d.get("name"), d.get("det_type")) for d in detectors]}')

    feu_connectors = {}
    for _conn_name, (feu, connector) in chosen['dream_feus'].items():
        feu_connectors.setdefault(int(feu), set()).add(int(connector))

    feu_connectors = {f: sorted(c) for f, c in sorted(feu_connectors.items())}
    feu_channels = {}
    for feu, conns in feu_connectors.items():
        chans = np.concatenate([
            np.arange((c - 1) * CH_PER_CONNECTOR, c * CH_PER_CONNECTOR)
            for c in conns
        ])
        feu_channels[feu] = np.sort(chans)

    return feu_connectors, feu_channels, chosen.get('name', det_type)


def parse_reference_feus(run_config_path, ref_det_type='mx17', exclude_feus=()):
    """Return (sorted FEU list, name) of the reference tracker used to tag
    cosmic tracks. Picks the detector whose det_type == ref_det_type; any FEUs
    in exclude_feus (e.g. the P2 FEUs) are dropped. Returns ([], None) if none.
    """
    with open(run_config_path) as fh:
        detectors = json.load(fh).get('detectors', [])
    for d in detectors:
        if d.get('det_type') == ref_det_type:
            feus = sorted({int(f) for f, _c in d['dream_feus'].values()}
                          - set(exclude_feus))
            return feus, d.get('name', ref_det_type)
    return [], None


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_p2_hits(hits_dir, p2_feus):
    """Load only the P2 FEU hits from every combined-hits ROOT file in hits_dir.

    Filtering per file keeps the reference-detector FEUs (which can carry >1M
    noisy hits) out of memory. eventId and trigger_timestamp_ns are already
    continuous across files, so a plain concat preserves the run time axis.
    """
    files = sorted(f for f in glob.glob(os.path.join(hits_dir, '*.root'))
                   if '_datrun_' in os.path.basename(f)
                   or 'combined_hits' in os.path.basename(f))
    if not files:
        raise FileNotFoundError(f'No combined-hits ROOT files in {hits_dir}')

    p2_set = set(p2_feus)
    parts = []
    for fp in files:
        arr = uproot.open(f'{fp}:hits').arrays(_BRANCHES, library='pd')
        parts.append(arr[arr['feu'].isin(p2_set)].copy())
    df = pd.concat(parts, ignore_index=True)

    print(f'Loaded {len(df):,} P2 hits (FEUs {sorted(p2_set)}) from '
          f'{len(files)} files, over {df["eventId"].nunique():,} P2-firing events.')
    return df, files


def load_reference_hits(files, ref_feus):
    """Load minimal (eventId, time) reference hits for the timing coincidence.

    Only the reference FEU rows and only the branches the coincidence needs are
    kept, so the >1M-hit noisy reference stays manageable.
    """
    ref_set = set(ref_feus)
    parts = []
    for fp in files:
        arr = uproot.open(f'{fp}:hits').arrays(_REF_BRANCHES, library='pd')
        parts.append(arr[arr['feu'].isin(ref_set)][['eventId', 'time']].copy())
    ref = pd.concat(parts, ignore_index=True)
    print(f'Loaded {len(ref):,} reference hits (FEUs {sorted(ref_set)}) for '
          f'timing coincidence.')
    return ref


def load_m3_tracks(m3_dir):
    """Load reconstructed m3 telescope tracks, one row per event.

    The m3_tracking_root tree 'T' has evn (== combined-hits eventId), rayN
    (number of reconstructed cosmic tracks) and evttime. A few evn values repeat
    across files; we keep the max rayN per event. Returns a DataFrame indexed by
    eventId with columns: n_ray, has_track.
    """
    files = sorted(glob.glob(os.path.join(m3_dir, '*.root')))
    if not files:
        raise FileNotFoundError(f'No m3 tracking ROOT files in {m3_dir}')
    parts = []
    for fp in files:
        parts.append(uproot.open(f'{fp}:T').arrays(['evn', 'rayN'], library='pd'))
    m3 = pd.concat(parts, ignore_index=True)
    m3 = (m3.groupby('evn')['rayN'].max()
          .rename('n_ray').to_frame())
    m3.index.name = 'eventId'
    m3['has_track'] = m3['n_ray'] >= 1
    print(f'Loaded m3 tracks for {len(m3):,} events '
          f'({int(m3["has_track"].sum()):,} with >=1 reconstructed ray) '
          f'from {len(files)} files.')
    return m3


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------

def _save_fig(fig, out_dir, name, dpi=150):
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, name)
    fig.savefig(path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f'  wrote {path}')


def _feu_colors(feus):
    cmap = plt.get_cmap('tab10')
    return {f: cmap(i % 10) for i, f in enumerate(feus)}


# ---------------------------------------------------------------------------
# Overall P2 summary
# ---------------------------------------------------------------------------

def plot_summary(df, title, out_dir):
    """2x2 P2 summary: pulse height, multiplicity, saturation/FEU, rate vs time."""
    feus = sorted(df['feu'].unique())
    colors = _feu_colors(feus)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # (1) Pulse-height distribution -------------------------------------------
    ax = axes[0, 0]
    amp = df['amplitude'].values
    hi = np.percentile(amp, 99.9)
    bins = np.linspace(0, max(hi, 1.0), 120)
    ax.hist(amp, bins=bins, color='0.7', edgecolor='none', label='all P2')
    for f in feus:
        ax.hist(df.loc[df['feu'] == f, 'amplitude'], bins=bins,
                histtype='step', lw=1.5, color=colors[f], label=f'FEU {f}')
    ax.set_yscale('log')
    ax.set_xlabel(AMP_LABEL)
    ax.set_ylabel('Hits / bin')
    ax.set_title(f'Pulse-height distribution (x-axis clipped at 99.9% = '
                 f'{hi:.0f} ADC)')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # (2) Hit multiplicity (P2 hits per event) --------------------------------
    ax = axes[0, 1]
    mult = df.groupby('eventId').size().values
    m_hi = int(np.percentile(mult, 99))
    ax.hist(mult, bins=np.arange(0.5, max(m_hi, 2) + 1.5, 1),
            color='steelblue', edgecolor='none')
    ax.set_xlabel('P2 hits per event')
    ax.set_ylabel('Events / bin')
    ax.set_title(f'Hit multiplicity  (mean = {mult.mean():.1f}, '
                 f'median = {int(np.median(mult))})')
    ax.grid(True, axis='y', alpha=0.3)

    # (3) Saturated hits per FEU ----------------------------------------------
    ax = axes[1, 0]
    sat_counts, tot_counts, fracs = [], [], []
    for f in feus:
        sub = df[df['feu'] == f]
        s = int(sub['saturated'].sum())
        t = len(sub)
        sat_counts.append(s)
        tot_counts.append(t)
        fracs.append(100.0 * s / t if t else 0.0)
    xpos = np.arange(len(feus))
    bars = ax.bar(xpos, sat_counts, color=[colors[f] for f in feus],
                  edgecolor='black', lw=0.5)
    for x, s, frac in zip(xpos, sat_counts, fracs):
        ax.text(x, s, f'{frac:.1f}%', ha='center', va='bottom', fontsize=9)
    ax.set_xticks(xpos)
    ax.set_xticklabels([f'FEU {f}' for f in feus])
    ax.set_ylabel('Saturated hits')
    overall = 100.0 * sum(sat_counts) / sum(tot_counts) if sum(tot_counts) else 0.0
    ax.set_title(f'Saturated hits per FEU  (overall {overall:.1f}%)')
    ax.grid(True, axis='y', alpha=0.3)

    # (4) Hit rate vs time ----------------------------------------------------
    ax = axes[1, 1]
    ts = df['trigger_timestamp_ns'].values
    t_sec = (ts - ts.min()) / 1e9
    n_bins = 120
    counts, edges = np.histogram(t_sec, bins=n_bins)
    width = edges[1] - edges[0]
    centres = 0.5 * (edges[:-1] + edges[1:])
    ax.plot(centres, counts / width, color='steelblue', lw=1.3)
    mean_rate = len(t_sec) / (t_sec.max() - t_sec.min())
    ax.axhline(mean_rate, color='crimson', ls='--', lw=1,
               label=f'mean {mean_rate:.1f} Hz')
    ax.set_xlabel('Time since run start [s]')
    ax.set_ylabel('P2 hit rate [Hz]')
    ax.set_title(f'Hit rate vs time  (run = {t_sec.max() / 3600:.2f} h)')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    fig.suptitle(title, fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    _save_fig(fig, out_dir, 'p2_summary.png')


# ---------------------------------------------------------------------------
# Per-FEU channel QA
# ---------------------------------------------------------------------------

def plot_feu_channels(df, feu, instrumented_ch, title, out_dir):
    """Per-FEU figure: hit occupancy + pulse height per channel.

    instrumented_ch : channels wired for this FEU (from run_config connectors).
                      Channels outside this set are shaded -- a hit there is a
                      mapping/noise red flag.
    """
    sub = df[df['feu'] == feu]
    ch = sub['channel'].values
    amp = sub['amplitude'].values

    inst = np.asarray(instrumented_ch)
    ch_lo, ch_hi = int(inst.min()), int(inst.max())
    edges = np.arange(ch_lo - 0.5, ch_hi + 1.5, 1)
    centres = np.arange(ch_lo, ch_hi + 1)

    # occupancy + dead-channel count
    occ, _ = np.histogram(ch, bins=edges)
    dead = int(np.sum(occ[np.isin(centres, inst)] == 0))
    n_inst = len(inst)

    fig, axes = plt.subplots(3, 1, figsize=(13, 11), sharex=True)

    # connector boundaries every 64 channels
    def _draw_connectors(ax):
        for b in range(0, N_CH_PER_FEU + 1, CH_PER_CONNECTOR):
            if ch_lo - 1 <= b <= ch_hi + 1:
                ax.axvline(b - 0.5, color='0.8', lw=0.7, zorder=0)

    # shade channels that are NOT instrumented for this FEU
    def _shade_uninstrumented(ax):
        not_inst = np.setdiff1d(centres, inst)
        for c in not_inst:
            ax.axvspan(c - 0.5, c + 0.5, color='mistyrose', zorder=-1, lw=0)

    # (1) Occupancy -----------------------------------------------------------
    ax = axes[0]
    _shade_uninstrumented(ax)
    _draw_connectors(ax)
    ax.bar(centres, occ, width=1.0, color='steelblue', edgecolor='none')
    ax.set_ylabel('Hits')
    ax.set_title(f'Occupancy — hits vs channel   '
                 f'({n_inst} instrumented ch, {dead} with zero hits, '
                 f'pink = not wired)')
    ax.grid(True, axis='y', alpha=0.3)

    # (2) Mean pulse height per channel (profile) -----------------------------
    ax = axes[1]
    _shade_uninstrumented(ax)
    _draw_connectors(ax)
    sum_amp, _ = np.histogram(ch, bins=edges, weights=amp)
    with np.errstate(invalid='ignore', divide='ignore'):
        mean_amp = np.where(occ > 0, sum_amp / occ, np.nan)
    ax.plot(centres, mean_amp, color='darkorange', lw=1.0, marker='.', ms=2)
    ax.set_ylabel('Mean pulse height [ADC]')
    ax.set_title('Mean pulse height per channel')
    ax.grid(True, alpha=0.3)

    # (3) Pulse height per channel — 2D density -------------------------------
    ax = axes[2]
    amp_hi = np.percentile(amp, 99.5) if len(amp) else 1.0
    amp_bins = np.linspace(0, max(amp_hi, 1.0), 100)
    h = ax.hist2d(ch, amp, bins=[edges, amp_bins], cmap='viridis',
                  cmin=1)
    fig.colorbar(h[3], ax=ax, label='Hits / bin')
    ax.set_xlabel('Channel number')
    ax.set_ylabel(AMP_LABEL)
    ax.set_title(f'Pulse height per channel (2D, y clipped at 99.5% = '
                 f'{amp_hi:.0f} ADC)')

    fig.suptitle(f'{title}\nFEU {feu} — {len(sub):,} hits', fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    _save_fig(fig, out_dir, f'p2_feu{feu}_channels.png')

    return dict(feu=feu, n_hits=len(sub), n_instrumented=n_inst,
                n_dead=dead, sat_frac=100.0 * sub['saturated'].mean()
                if len(sub) else 0.0)


# ---------------------------------------------------------------------------
# Reference-coincidence (cosmic-track-tagged) QA
# ---------------------------------------------------------------------------

def tag_reference_coincidence(df_p2, df_ref):
    """For every P2 hit, find the nearest-in-time reference hit in the SAME
    event and return its signed dt = t_P2 - t_ref [ns].

    Vectorised with merge_asof (nearest match, grouped by eventId). P2 hits in
    events with no reference hit get dt = NaN. Returns a copy of df_p2 with an
    added 'dt_ref_ns' column. merge_asof requires distinct left/right key names,
    so the times are renamed before matching.
    """
    out = df_p2.copy()
    left = (df_p2[['eventId', 'time']].rename(columns={'time': 't_p2'})
            .reset_index()                       # 'index' maps the match back
            .sort_values('t_p2', kind='mergesort'))
    right = (df_ref[['eventId', 'time']].rename(columns={'time': 't_ref'})
             .sort_values('t_ref', kind='mergesort'))
    m = pd.merge_asof(left, right, left_on='t_p2', right_on='t_ref',
                      by='eventId', direction='nearest')
    m['dt'] = m['t_p2'] - m['t_ref']
    out['dt_ref_ns'] = m.set_index('index')['dt']
    return out


def plot_correlated(df_p2, df_ref, feu_channels, title, out_dir,
                    window_ns=100.0, ref_name='reference'):
    """Cosmic-track-tagged QA via timing coincidence with the reference tracker.

    Four panels:
      (1) dt = t_P2 - t_ref distribution -> coincidence peak over flat accidentals
      (2) pulse-height spectrum: all P2 hits vs reference-coincident hits
      (3) per-FEU occupancy: all vs coincident (are channel patterns the same?)
      (4) P2 self-coincidence: fraction of P2 hits that are reference-tagged,
          per FEU, plus the multiplicity of P2 FEUs firing together.
    """
    feus = sorted(df_p2['feu'].unique())
    colors = _feu_colors(feus)
    tagged = df_p2[df_p2['dt_ref_ns'].abs() <= window_ns]
    frac = len(tagged) / len(df_p2) if len(df_p2) else 0.0

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # (1) dt distribution -----------------------------------------------------
    ax = axes[0, 0]
    dt = df_p2['dt_ref_ns'].dropna().values
    rng = 6 * window_ns
    ax.hist(dt, bins=200, range=(-rng, rng), color='steelblue', edgecolor='none')
    ax.axvspan(-window_ns, window_ns, color='orange', alpha=0.25,
               label=f'±{window_ns:.0f} ns window')
    ax.set_yscale('log')
    ax.set_xlabel(r'$\Delta t = t_{P2} - t_{ref}$ [ns]')
    ax.set_ylabel('P2 hits / bin')
    ax.set_title(f'Timing coincidence with {ref_name}\n'
                 f'{frac * 100:.0f}% of P2 hits in window')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # (2) pulse-height: all vs tagged -----------------------------------------
    ax = axes[0, 1]
    amp_all = df_p2['amplitude'].values
    hi = np.percentile(amp_all, 99.9)
    bins = np.linspace(0, max(hi, 1.0), 120)
    ax.hist(amp_all, bins=bins, color='0.75', edgecolor='none',
            label=f'all P2 ({len(df_p2):,})')
    ax.hist(tagged['amplitude'], bins=bins, histtype='step', lw=1.8,
            color='crimson', label=f'{ref_name}-coincident ({len(tagged):,})')
    ax.set_yscale('log')
    ax.set_xlabel(AMP_LABEL)
    ax.set_ylabel('Hits / bin')
    ax.set_title('Pulse height — all vs cosmic-tagged')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # (3) per-FEU occupancy: all vs tagged ------------------------------------
    ax = axes[1, 0]
    width = 0.38
    xpos = np.arange(len(feus))
    all_counts = [int((df_p2['feu'] == f).sum()) for f in feus]
    tag_counts = [int((tagged['feu'] == f).sum()) for f in feus]
    ax.bar(xpos - width / 2, all_counts, width, color='0.7', label='all P2')
    ax.bar(xpos + width / 2, tag_counts, width,
           color=[colors[f] for f in feus], edgecolor='black', lw=0.5,
           label='cosmic-tagged')
    for x, a, t in zip(xpos, all_counts, tag_counts):
        ax.text(x + width / 2, t, f'{100 * t / a:.0f}%' if a else '–',
                ha='center', va='bottom', fontsize=9)
    ax.set_xticks(xpos)
    ax.set_xticklabels([f'FEU {f}' for f in feus])
    ax.set_ylabel('Hits')
    ax.set_title('Hits per FEU — all vs cosmic-tagged (label = tagged fraction)')
    ax.legend(fontsize=9)
    ax.grid(True, axis='y', alpha=0.3)

    # (4) P2 FEU multiplicity (self-coincidence) ------------------------------
    ax = axes[1, 1]
    feus_per_event = df_p2.groupby('eventId')['feu'].nunique()
    vals, cnts = np.unique(feus_per_event.values, return_counts=True)
    ax.bar(vals, cnts, color='seagreen', edgecolor='black', lw=0.5)
    for v, c in zip(vals, cnts):
        ax.text(v, c, f'{100 * c / cnts.sum():.0f}%', ha='center',
                va='bottom', fontsize=9)
    ax.set_xticks(vals)
    ax.set_xlabel('Number of P2 FEUs firing in the event')
    ax.set_ylabel('Events')
    ax.set_title('P2 self-coincidence (multi-FEU = likely real track)')
    ax.grid(True, axis='y', alpha=0.3)

    fig.suptitle(f'{title}\nCosmic-tagged QA — coincidence window ±{window_ns:.0f} ns',
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    _save_fig(fig, out_dir, 'p2_correlated.png')

    return frac


# ---------------------------------------------------------------------------
# m3-track-tagged QA (reconstructed-track reference, event-level match)
# ---------------------------------------------------------------------------

def tag_m3_tracks(df_p2, m3):
    """Attach m3 'n_ray'/'has_track' to each P2 hit by exact eventId match.

    The m3 telescope reconstructs cosmic tracks per event; evn matches the
    combined-hits eventId exactly, so this is an event-level tag (no timing
    window). Hits in events m3 never saw get has_track = False, n_ray = 0.
    """
    out = df_p2.copy()
    out['n_ray'] = out['eventId'].map(m3['n_ray']).fillna(0).astype(int)
    out['has_track'] = out['n_ray'] >= 1
    return out


def plot_correlated_m3(df_p2, m3, title, out_dir, ref_name='m3'):
    """Cosmic-track-tagged QA using reconstructed m3 telescope tracks.

    Four panels:
      (1) pulse-height spectrum: all P2 hits vs hits in m3-tracked events
      (2) per-FEU hits: all vs m3-tracked (label = tracked fraction)
      (3) P2 hit multiplicity for events WITH vs WITHOUT a reconstructed track
      (4) reference context: m3 track multiplicity (rayN) + event-count summary
    """
    feus = sorted(df_p2['feu'].unique())
    colors = _feu_colors(feus)
    tagged = df_p2[df_p2['has_track']]
    frac = len(tagged) / len(df_p2) if len(df_p2) else 0.0

    # event-level counts
    p2_events = df_p2['eventId'].unique()
    n_p2 = len(p2_events)
    n_p2_track = int(pd.Series(p2_events).map(m3['has_track']).fillna(False).sum())

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # (1) pulse-height: all vs m3-tagged -------------------------------------
    ax = axes[0, 0]
    amp_all = df_p2['amplitude'].values
    hi = np.percentile(amp_all, 99.9)
    bins = np.linspace(0, max(hi, 1.0), 120)
    ax.hist(amp_all, bins=bins, color='0.75', edgecolor='none',
            label=f'all P2 ({len(df_p2):,})')
    ax.hist(tagged['amplitude'], bins=bins, histtype='step', lw=1.8,
            color='crimson', label=f'{ref_name}-track tagged ({len(tagged):,})')
    ax.set_yscale('log')
    ax.set_xlabel(AMP_LABEL)
    ax.set_ylabel('Hits / bin')
    ax.set_title('Pulse height — all vs cosmic-track-tagged')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # (2) per-FEU hits: all vs tagged ----------------------------------------
    ax = axes[0, 1]
    width = 0.38
    xpos = np.arange(len(feus))
    all_counts = [int((df_p2['feu'] == f).sum()) for f in feus]
    tag_counts = [int((tagged['feu'] == f).sum()) for f in feus]
    ax.bar(xpos - width / 2, all_counts, width, color='0.7', label='all P2')
    ax.bar(xpos + width / 2, tag_counts, width,
           color=[colors[f] for f in feus], edgecolor='black', lw=0.5,
           label='track-tagged')
    for x, a, t in zip(xpos, all_counts, tag_counts):
        ax.text(x + width / 2, t, f'{100 * t / a:.0f}%' if a else '–',
                ha='center', va='bottom', fontsize=9)
    ax.set_xticks(xpos)
    ax.set_xticklabels([f'FEU {f}' for f in feus])
    ax.set_ylabel('Hits')
    ax.set_title('Hits per FEU — all vs tagged (label = tracked fraction)')
    ax.legend(fontsize=9)
    ax.grid(True, axis='y', alpha=0.3)

    # (3) P2 multiplicity: track vs no-track events --------------------------
    ax = axes[1, 0]
    mult = df_p2.groupby('eventId').size()
    ev_has_track = (pd.Series(mult.index).map(m3['has_track'])
                    .fillna(False).astype(bool).values)
    m_hi = int(np.percentile(mult.values, 99))
    mbins = np.arange(0.5, max(m_hi, 2) + 1.5, 1)
    ax.hist(mult.values[ev_has_track], bins=mbins, histtype='stepfilled',
            color='crimson', alpha=0.55, label='event has m3 track')
    ax.hist(mult.values[~ev_has_track], bins=mbins, histtype='step', lw=1.8,
            color='navy', label='no m3 track')
    ax.set_yscale('log')
    ax.set_xlabel('P2 hits per event')
    ax.set_ylabel('Events / bin')
    ax.set_title('P2 multiplicity — track vs no-track events')
    ax.legend(fontsize=9)
    ax.grid(True, axis='y', alpha=0.3)

    # (4) m3 reference context: rayN + counts --------------------------------
    ax = axes[1, 1]
    ray_vals, ray_cnts = np.unique(m3['n_ray'].values, return_counts=True)
    ax.bar(ray_vals, ray_cnts, color='seagreen', edgecolor='black', lw=0.5)
    ax.set_yscale('log')
    ax.set_xticks(ray_vals)
    ax.set_xlabel('m3 reconstructed tracks per event (rayN)')
    ax.set_ylabel('Events')
    ax.set_title('m3 reference — track multiplicity')
    ax.grid(True, axis='y', alpha=0.3)
    txt = (f'P2-firing events: {n_p2:,}\n'
           f'  with m3 track:  {n_p2_track:,} ({100 * n_p2_track / n_p2:.0f}%)\n'
           f'P2 hits tagged:   {frac * 100:.0f}%')
    ax.text(0.97, 0.95, txt, transform=ax.transAxes, ha='right', va='top',
            fontsize=10, family='monospace',
            bbox=dict(boxstyle='round', fc='white', ec='0.6'))

    fig.suptitle(f'{title}\nCosmic-track-tagged QA — reference: {ref_name} '
                 f'reconstructed tracks', fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    _save_fig(fig, out_dir, 'p2_correlated_m3.png')

    return frac, n_p2_track / n_p2 if n_p2 else 0.0


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description='P2 raw channel QA from combined-hits ROOT files.')
    ap.add_argument('--hits-dir', required=True,
                    help='directory with *_feu-combined_hits.root files')
    ap.add_argument('--run-config', default=None,
                    help='run_config.json (default: two dirs above --hits-dir)')
    ap.add_argument('--det-type', default='P2',
                    help='detector det_type to select in the config (default P2)')
    ap.add_argument('--det-name', default=None,
                    help='detector name to select (overrides --det-type)')
    ap.add_argument('--out-dir', default=None,
                    help='output directory (default: ./plots/p2_qa/<run>/<subrun>)')
    ap.add_argument('--ref-mode', default='auto',
                    choices=['auto', 'm3', 'mx17', 'none'],
                    help='cosmic-track tag reference: "m3" = reconstructed m3 '
                         'telescope tracks (event-level, cleanest); "mx17" = '
                         'timing coincidence with mx17 strip hits; "auto" = m3 '
                         'if m3_tracking_root is present, else mx17; "none" = '
                         'skip correlated QA.')
    ap.add_argument('--m3-dir', default=None,
                    help='m3_tracking_root dir (default: sibling of --hits-dir)')
    ap.add_argument('--ref-det-type', default='mx17',
                    help='det_type of the reference strip tracker for mx17-mode '
                         'timing coincidence (default mx17).')
    ap.add_argument('--coinc-window-ns', type=float, default=100.0,
                    help='|t_P2 - t_ref| window for the mx17 timing tag [ns]')
    args = ap.parse_args()

    hits_dir = os.path.abspath(args.hits_dir)

    # default run_config.json: <run>/<sub_run>/combined_hits_root -> <run>/run_config.json
    run_config = args.run_config
    if run_config is None:
        run_root = os.path.dirname(os.path.dirname(hits_dir))
        run_config = os.path.join(run_root, 'run_config.json')
    if not os.path.isfile(run_config):
        raise FileNotFoundError(f'run_config.json not found at {run_config} '
                                f'(pass --run-config explicitly)')

    feu_connectors, feu_channels, det_name = parse_p2_feus(
        run_config, det_type=args.det_type, det_name=args.det_name)
    p2_feus = sorted(feu_connectors)

    with open(run_config) as fh:
        run_name = json.load(fh).get('run_name', os.path.basename(os.path.dirname(run_config)))
    sub_run = os.path.basename(os.path.dirname(hits_dir))

    print(f'Detector  : {det_name}  ({args.det_type})')
    print(f'Run       : {run_name} / {sub_run}')
    for f in p2_feus:
        print(f'  FEU {f}: connectors {feu_connectors[f]} '
              f'-> channels {feu_channels[f].min()}-{feu_channels[f].max()} '
              f'({len(feu_channels[f])} ch)')

    out_dir = args.out_dir or os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        'plots', 'p2_qa', run_name, sub_run)

    df, files = load_p2_hits(hits_dir, p2_feus)

    title = f'P2 channel QA — {det_name} — {run_name} / {sub_run}'
    plot_summary(df, title, out_dir)

    print('\nPer-FEU channel QA:')
    stats = []
    for f in p2_feus:
        stats.append(plot_feu_channels(df, f, feu_channels[f], title, out_dir))

    # --- cosmic-track-tagged QA ---
    # Pick the reference: m3 reconstructed tracks (cleanest, event-level match)
    # if available, else mx17 strip-hit timing coincidence.
    m3_dir = args.m3_dir or os.path.join(os.path.dirname(hits_dir), 'm3_tracking_root')
    m3_available = os.path.isdir(m3_dir) and bool(glob.glob(os.path.join(m3_dir, '*.root')))
    mode = args.ref_mode
    if mode == 'auto':
        mode = 'm3' if m3_available else 'mx17'

    if mode == 'm3':
        if not m3_available:
            print(f'\nm3 mode requested but no m3 tracks at {m3_dir} — '
                  f'skipping correlated QA.')
        else:
            print('\nCosmic-track-tagged QA vs m3 reconstructed tracks:')
            m3 = load_m3_tracks(m3_dir)
            df = tag_m3_tracks(df, m3)
            frac, ev_frac = plot_correlated_m3(df, m3, title, out_dir, ref_name='m3')
            print(f'  {ev_frac * 100:.1f}% of P2 events have an m3 track; '
                  f'{frac * 100:.1f}% of P2 hits tagged.')
    elif mode == 'mx17':
        ref_feus, ref_name = parse_reference_feus(
            run_config, ref_det_type=args.ref_det_type, exclude_feus=p2_feus)
        if not ref_feus:
            print(f'\nNo reference detector (det_type={args.ref_det_type!r}) '
                  f'in config — skipping correlated QA.')
        else:
            print(f'\nCosmic-tagged QA vs reference {ref_name} (FEUs {ref_feus}):')
            df_ref = load_reference_hits(files, ref_feus)
            df = tag_reference_coincidence(df, df_ref)
            frac = plot_correlated(df, df_ref, feu_channels, title, out_dir,
                                   window_ns=args.coinc_window_ns,
                                   ref_name=ref_name)
            print(f'  {frac * 100:.1f}% of P2 hits coincident within '
                  f'±{args.coinc_window_ns:.0f} ns')

    # text summary
    print('\n' + '=' * 64)
    print(f'{"FEU":>4} {"hits":>10} {"instr.ch":>9} {"dead":>6} {"sat%":>7}')
    for s in stats:
        print(f'{s["feu"]:>4} {s["n_hits"]:>10,} {s["n_instrumented"]:>9} '
              f'{s["n_dead"]:>6} {s["sat_frac"]:>6.1f}%')
    print('=' * 64)
    print(f'\nAll plots written under: {out_dir}')


if __name__ == '__main__':
    main()
