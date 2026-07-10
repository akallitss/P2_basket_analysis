#!/usr/bin/env python3
"""
m3_comparison.py

Compares M3 reference tracker reconstruction between a non-zero-suppressed
(non-ZS) and a zero-suppressed (ZS) run on the cosmic bench.

The M3 tracker is the reference instrument — if ZS degrades its reconstruction
quality, any downstream efficiency or signal comparison is unreliable.

Metrics
-------
  - Track multiplicity per event (rayN)
  - Chi2 distributions in X and Y (before and after cut)
  - Single-good-track yield at each stage of the quality cuts
  - Beam profile at the top M3 plane
  - Track angle distributions (theta_x, theta_y)
  - Cosmic rate vs time (confirms both runs sampled the same flux)
"""

import sys
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import awkward as ak

# M3RefTracking lives in the nTof_x17 analysis repo
sys.path.insert(0, os.path.expanduser(
    '~/Documents/PostDocSaclay/nTof_x17/cosmic_bench_analysis'))
from M3RefTracking import M3RefTracking, get_xy_angles

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

DATA_BASE = os.path.expanduser(
    '~/Documents/PostDocSaclay/data/Cosmic_Bench/'
    'mx17_non_and_ZS_Test_6-17-26'
)

RUNS = {
    'Non-ZS': f'{DATA_BASE}/run_non_zs/m3_tracking_root/',
    'ZS':     f'{DATA_BASE}/run_zs/m3_tracking_root/',
}

COLORS   = {'Non-ZS': 'steelblue', 'ZS': 'tomato'}
# Tracking-v2 recommended recipe: Chi2X,Chi2Y < CHI2_CUT AND NClusX,NClusY >= MIN_NCLUS.
# (Matches p2_qa_config.M3_CHI2_CUT / M3_MIN_NCLUS; kept local so this standalone
# ZS/non-ZS diagnostic has no pipeline dependency.) MIN_NCLUS drops the exact
# 2-point-per-coordinate fits whose denormal-tiny chi2 slips through a chi2-only cut.
CHI2_CUT = 5.0
MIN_NCLUS = 3
EVTTIME_TO_S = 1e-8   # evttime is in 10 ns ticks → seconds

PLOTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'plots')

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def to_np(arr):
    """Convert awkward or numpy array to a flat numpy array."""
    if isinstance(arr, ak.highlevel.Array):
        return ak.to_numpy(arr)
    return np.asarray(arr)


def flatten_np(arr):
    """Flatten a ragged awkward array to a 1-D numpy array."""
    if isinstance(arr, ak.highlevel.Array):
        return ak.to_numpy(ak.flatten(arr))
    return np.asarray(arr).ravel()

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------

def load_runs(chi2_cut=CHI2_CUT):
    """
    Load each run twice:
      raw  — no quality cuts (single_track=False): gives access to full rayN
              and chi2 distributions before any selection.
      good — standard quality cuts (single_track=True): events with exactly
              one track passing chi2 and detector-size cuts.
    """
    raw, good = {}, {}
    for label, ray_dir in RUNS.items():
        print(f'\n--- Loading {label} (no cuts) ---')
        raw[label]  = M3RefTracking(ray_dir, single_track=False, chi2_cut=chi2_cut,
                                    min_nclus=MIN_NCLUS)
        print(f'--- Loading {label} (with quality cuts) ---')
        good[label] = M3RefTracking(ray_dir, single_track=True,  chi2_cut=chi2_cut,
                                    min_nclus=MIN_NCLUS)
    return raw, good

# ---------------------------------------------------------------------------
# Per-run summary
# ---------------------------------------------------------------------------

def compute_summary(raw, good):
    rows = {}
    for label in RUNS:
        rd   = raw[label].ray_data
        gd   = good[label].ray_data

        rayn        = to_np(rd['rayN'])
        n_total     = len(rayn)
        n_has_track = int(np.sum(rayn > 0))
        n_good      = len(to_np(gd['X_Up']))

        chi2x = flatten_np(rd['Chi2X'])
        chi2y = flatten_np(rd['Chi2Y'])

        t = to_np(rd['evttime']) * EVTTIME_TO_S
        duration_s = float(np.max(t) - np.min(t)) if len(t) > 1 else 1.0

        rows[label] = {
            'n_total':               n_total,
            'n_has_track':           n_has_track,
            'track_reco_rate_%':     100 * n_has_track / n_total,
            'n_single_good_track':   n_good,
            'good_track_rate_%':     100 * n_good / n_total,
            'median_chi2x':          float(np.median(chi2x)),
            'median_chi2y':          float(np.median(chi2y)),
            'duration_s':            duration_s,
            'good_tracks_per_s':     n_good / duration_s if duration_s > 0 else 0,
        }
    return rows


def print_summary(rows):
    labels = list(RUNS.keys())
    metrics = [
        ('Total events',             'n_total',             '{:,}'),
        ('Events with ≥1 track',     'n_has_track',         '{:,}'),
        ('Track reco rate [%]',      'track_reco_rate_%',   '{:.1f}'),
        ('Single good-track events', 'n_single_good_track', '{:,}'),
        ('Good-track rate [%]',      'good_track_rate_%',   '{:.1f}'),
        ('Median Chi2 X',            'median_chi2x',        '{:.3f}'),
        ('Median Chi2 Y',            'median_chi2y',        '{:.3f}'),
        ('Run duration [s]',         'duration_s',          '{:.1f}'),
        ('Good tracks / s',          'good_tracks_per_s',   '{:.2f}'),
    ]
    header = f'{"Metric":<32}' + ''.join(f'{lb:>14}' for lb in labels)
    print(f'\n{"="*60}')
    print('M3 Tracker Quality Summary')
    print('='*60)
    print(header)
    print('-'*60)
    for name, key, fmt in metrics:
        row = f'{name:<32}' + ''.join(fmt.format(rows[lb][key]) + ' ' * (13 - len(fmt.format(rows[lb][key]))) for lb in labels)
        print(row)

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def plot_multiplicity(raw, ax):
    bins = np.arange(-0.5, 8.5)
    for label, m3 in raw.items():
        rayn = to_np(m3.ray_data['rayN'])
        ax.hist(rayn, bins=bins, histtype='step', lw=2,
                color=COLORS[label],
                label=f'{label}  (N={len(rayn):,})',
                density=True)
    ax.set_xticks(range(8))
    ax.set_xlabel('Tracks per event (rayN)')
    ax.set_ylabel('Fraction of events')
    ax.set_title('Track multiplicity')
    ax.legend()


def plot_chi2(raw, ax, coord):
    """coord: 'X' or 'Y'"""
    bins = np.linspace(0, 15, 75)
    for label, m3 in raw.items():
        chi2 = flatten_np(m3.ray_data[f'Chi2{coord}'])
        ax.hist(chi2, bins=bins, histtype='step', lw=2,
                color=COLORS[label],
                label=f'{label}  (N={len(chi2):,})',
                density=True)
    ax.axvline(CHI2_CUT, color='k', ls='--', lw=1.2, label=f'cut = {CHI2_CUT}')
    ax.set_yscale('log')
    ax.set_xlabel('χ²')
    ax.set_ylabel('Density')
    ax.set_title(f'Chi2 {coord}  (all tracks, before cut)')
    ax.legend(fontsize=8)


def plot_yield_bar(rows, ax):
    labels_list = list(RUNS.keys())
    stages = ['n_total', 'n_has_track', 'n_single_good_track']
    stage_labels = ['All events', '≥1 track\n(raw)', 'Single good\ntrack (cut)']
    x = np.arange(len(stages))
    width = 0.35

    for i, label in enumerate(labels_list):
        vals = [rows[label][s] for s in stages]
        offset = (i - 0.5) * width
        bars = ax.bar(x + offset, vals, width,
                      label=label, color=COLORS[label], alpha=0.85)
        for bar, val in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() * 1.02,
                    f'{val:,}', ha='center', va='bottom', fontsize=8)

    ax.set_xticks(x)
    ax.set_xticklabels(stage_labels, fontsize=9)
    ax.set_ylabel('Event count')
    ax.set_title('Track yield at each cut stage')
    ax.legend()


def plot_beam_profile(good, axes):
    for ax, label in zip(axes, RUNS.keys()):
        m3 = good[label]
        x = to_np(m3.ray_data['X_Up'])
        y = to_np(m3.ray_data['Y_Up'])
        h = ax.hist2d(x, y, bins=40,
                      range=[[-250, 250], [-250, 250]],
                      cmap='viridis')
        plt.colorbar(h[3], ax=ax, label='Tracks')
        ax.set_xlabel('X at top M3 [mm]')
        ax.set_ylabel('Y at top M3 [mm]')
        ax.set_title(f'Beam profile — {label}\n({len(x):,} single good tracks)')


def plot_angles(good, ax_x, ax_y):
    bins = np.linspace(-20, 20, 60)
    for label, m3 in good.items():
        tx, ty, _ = get_xy_angles(m3.ray_data)
        tx_deg = np.degrees(np.asarray(tx))
        ty_deg = np.degrees(np.asarray(ty))
        kw = dict(bins=bins, histtype='step', lw=2,
                  color=COLORS[label], label=label, density=True)
        ax_x.hist(tx_deg, **kw)
        ax_y.hist(ty_deg, **kw)
    for ax, coord in [(ax_x, 'X'), (ax_y, 'Y')]:
        ax.set_xlabel(f'θ_{coord.lower()} [deg]')
        ax.set_ylabel('Density')
        ax.set_title(f'Track angle θ_{coord.lower()}')
        ax.legend()


def plot_rate_vs_time(raw, ax, n_bins=20):
    """Track count per time bin — confirms both runs saw the same cosmic rate."""
    for label, m3 in raw.items():
        evttime = to_np(m3.ray_data['evttime']) * EVTTIME_TO_S
        rayn    = to_np(m3.ray_data['rayN'])
        has_track = rayn > 0
        t_min, t_max = evttime.min(), evttime.max()
        bins = np.linspace(t_min, t_max, n_bins + 1)
        counts, edges = np.histogram(evttime[has_track], bins=bins)
        dt = np.diff(edges)
        centres = 0.5 * (edges[:-1] + edges[1:])
        ax.step(centres - t_min, counts / dt,
                where='mid', lw=2, color=COLORS[label], label=label)
    ax.set_xlabel('Time from run start [s]')
    ax.set_ylabel('Reconstructed tracks / s')
    ax.set_title('Cosmic track rate vs time')
    ax.legend()

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    os.makedirs(PLOTS_DIR, exist_ok=True)

    raw, good = load_runs()
    rows = compute_summary(raw, good)
    print_summary(rows)

    # --- Figure 1: track quality (2 × 3 grid) ---
    fig1, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig1.suptitle('M3 Reference Tracker — Non-ZS vs ZS', fontsize=14, y=1.01)

    plot_multiplicity(raw, axes[0, 0])
    plot_chi2(raw, axes[0, 1], 'X')
    plot_chi2(raw, axes[0, 2], 'Y')
    plot_yield_bar(rows, axes[1, 0])
    plot_beam_profile(good, [axes[1, 1], axes[1, 2]])

    fig1.tight_layout()
    out1 = os.path.join(PLOTS_DIR, 'm3_quality_comparison.png')
    fig1.savefig(out1, dpi=150, bbox_inches='tight')
    print(f'\nSaved → {out1}')

    # --- Figure 2: angles + rate (1 × 3 grid) ---
    fig2, axes2 = plt.subplots(1, 3, figsize=(18, 5))
    fig2.suptitle('M3 Track Angles and Rate — Non-ZS vs ZS', fontsize=14)

    plot_angles(good, axes2[0], axes2[1])
    plot_rate_vs_time(raw, axes2[2])

    fig2.tight_layout()
    out2 = os.path.join(PLOTS_DIR, 'm3_angle_rate_comparison.png')
    fig2.savefig(out2, dpi=150, bbox_inches='tight')
    print(f'Saved → {out2}')

    # plt.show()  # disabled: use Agg backend for headless runs


if __name__ == '__main__':
    main()
