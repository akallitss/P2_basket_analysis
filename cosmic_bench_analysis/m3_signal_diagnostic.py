#!/usr/bin/env python3
"""
m3_signal_diagnostic.py

Reads the cluster-level analyse trees produced by run_m3_pipeline.sh for both
the non-ZS and ZS runs, and answers:

  Q1. Do all 8 M3 detector planes produce clusters in ZS?
      (A blind plane → no track possible)
  Q2. Are ZS cluster sizes below the ClusSizeCut_Min=3 threshold?
      (Losing border strips due to ZS → cluster rejected)
  Q3. Is the cluster peak sample (MaxSample) within the cut window [7, 26]?
      (Timing shift in ZS data → cluster rejected)
  Q4. Is cluster TOT too low for the TOTCut=3 threshold?
  Q5. Are cluster amplitudes reasonable (positive) in ZS?
      (Double pedestal subtraction → negative amplitudes → garbage chi2)

Usage:
    python3 m3_signal_diagnostic.py

Output:
    plots/m3_cluster_occupancy.png  — clusters per plane per event
    plots/m3_cluster_properties.png — size, MaxSample, TOT, amplitude
"""

import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import uproot

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
TRACKING_DIR = os.path.expanduser(
    '~/Documents/PostDocSaclay/cosmic_bench_m3_tracking')

ANALYSE_FILES = {
    'Non-ZS': os.path.join(TRACKING_DIR, 'root_files/analyse_non_zs.root'),
    'ZS':     os.path.join(TRACKING_DIR, 'root_files/analyse_zs.root'),
}

PLOTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'plots')
COLORS = {'Non-ZS': 'steelblue', 'ZS': 'tomato'}

# ---------------------------------------------------------------------------
# M3 detector layout (from config_ref_2022.json, mg_n order 0-7)
# ---------------------------------------------------------------------------
N_DET    = 8
MAX_CLUS = 300   # MGv2_Detector::MaxNClus

DET_LABELS = [
    'mg0 Y  z=1302', 'mg1 X  z=1302',
    'mg2 Y  z=1185', 'mg3 X  z=1185',
    'mg4 Y  z=144',  'mg5 X  z=144',
    'mg6 X  z=24',   'mg7 Y  z=24',
]

# Config cuts (config_ref_2022.json / per-detector values)
CUT_SIZE_MIN   = {'mg0': 3, 'mg1': 2, 'mg2': 3, 'mg3': 2,
                  'mg4': 3, 'mg5': 2, 'mg6': 3, 'mg7': 2}
CUT_SAMPLE_MIN = 7
CUT_SAMPLE_MAX = 26
CUT_TOT_MIN    = 3   # global TOTCut

# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load(path):
    """
    Load the analyse ROOT tree.
    Returns a dict with numpy arrays for each cluster property.
    Branch shapes in ROOT: MGv2_NClus[8], MGv2_ClusSize[8][300], etc.
    Uproot reads them as (n_events, 8) and (n_events, 2400) respectively.
    """
    if not os.path.exists(path):
        print(f'  [ERROR] File not found: {path}')
        print('  → Run run_m3_pipeline.sh first.')
        return None

    f = uproot.open(path)
    t = f['T']
    print(f'  Branches: {[k for k in t.keys() if "MGv2" in k]}')

    keys = ['MGv2_NClus', 'MGv2_ClusSize', 'MGv2_ClusMaxSample',
            'MGv2_ClusTOT', 'MGv2_ClusAmpl', 'MGv2_ClusPos']
    raw = t.arrays(keys, library='np')

    n_ev = len(raw['MGv2_NClus'])

    d = {}
    d['n_events'] = n_ev
    # NClus: (n_events, 8)
    nclus = raw['MGv2_NClus']
    if nclus.ndim == 1:
        nclus = nclus.reshape(n_ev, N_DET)
    d['NClus'] = nclus

    # Cluster arrays: (n_events, 8*300) → reshape to (n_events, 8, 300)
    for key in ['ClusSize', 'ClusMaxSample', 'ClusTOT', 'ClusAmpl', 'ClusPos']:
        arr = raw[f'MGv2_{key}']
        if arr.ndim == 1:
            arr = arr.reshape(n_ev, N_DET * MAX_CLUS)
        d[key] = arr.reshape(n_ev, N_DET, MAX_CLUS)

    return d


def valid_clusters(d, det, key):
    """Return flat array of valid cluster values for detector index `det`."""
    nclus  = d['NClus'][:, det]           # (n_events,)
    values = d[key][:, det, :]            # (n_events, MAX_CLUS)
    idx = np.arange(MAX_CLUS)[np.newaxis, :]
    mask = idx < nclus[:, np.newaxis]
    return values[mask]


# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------

def print_summary(data):
    print(f'\n{"="*68}')
    print(f'{"Cluster yield per M3 detector plane":^68}')
    print(f'{"="*68}')
    header = f'{"Detector":<22}' + ''.join(f'{lb:>16}' for lb in data)
    print(header)
    print('-'*68)
    for det in range(N_DET):
        row = f'{DET_LABELS[det]:<22}'
        for label, d in data.items():
            n_with = int(np.sum(d['NClus'][:, det] > 0))
            frac   = 100 * n_with / d['n_events']
            row += f'  {n_with:>5,} / {d["n_events"]:>5,} ({frac:>4.1f}%)'
        print(row)

    print('-'*68)
    # Events with ≥1 cluster in ALL 8 planes
    row = f'{"All 8 planes":<22}'
    for label, d in data.items():
        has_all = np.all(d['NClus'] > 0, axis=1)
        n = int(np.sum(has_all))
        row += f'  {n:>5,} / {d["n_events"]:>5,} ({100*n/d["n_events"]:>4.1f}%)'
    print(row)


# ---------------------------------------------------------------------------
# Plot 1: Cluster occupancy per detector plane
# ---------------------------------------------------------------------------

def plot_occupancy(data, out_dir):
    fig, axes = plt.subplots(2, 4, figsize=(18, 8), sharey=False)
    fig.suptitle('M3 Cluster Occupancy per Detector Plane — Non-ZS vs ZS',
                 fontsize=13)

    for det, ax in enumerate(axes.ravel()):
        bins = np.arange(-0.5, 12.5)
        for label, d in data.items():
            nclus = d['NClus'][:, det]
            ax.hist(nclus, bins=bins, histtype='step', lw=2,
                    color=COLORS[label], density=True,
                    label=f'{label}  (≥1: {100*np.mean(nclus>0):.0f}%)')
        ax.set_title(DET_LABELS[det], fontsize=9)
        ax.set_xlabel('Clusters per event')
        ax.set_ylabel('Fraction')
        ax.legend(fontsize=7)
        ax.set_xticks(range(0, 12, 2))

    fig.tight_layout()
    out = os.path.join(out_dir, 'm3_cluster_occupancy.png')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    print(f'Saved → {out}')


# ---------------------------------------------------------------------------
# Plot 2: Cluster properties (all detectors combined)
# ---------------------------------------------------------------------------

def plot_properties(data, out_dir):
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('M3 Cluster Properties — Non-ZS vs ZS  (all planes combined)',
                 fontsize=13)

    # --- ClusSize ---
    ax = axes[0, 0]
    bins = np.arange(0.5, 25.5)
    for label, d in data.items():
        vals = np.concatenate([valid_clusters(d, det, 'ClusSize')
                               for det in range(N_DET)])
        ax.hist(vals, bins=bins, histtype='step', lw=2,
                color=COLORS[label], density=True,
                label=f'{label}  (N={len(vals):,})')
    # Show the per-detector ClusSizeCut_Min values (range 2-3)
    ax.axvline(2, color='k', ls=':', lw=1, label='cut_min = 2 (X planes)')
    ax.axvline(3, color='gray', ls='--', lw=1, label='cut_min = 3 (Y planes)')
    ax.set_xlabel('Cluster size [strips]')
    ax.set_ylabel('Density')
    ax.set_title('Cluster size')
    ax.legend(fontsize=8)

    # --- ClusMaxSample ---
    ax = axes[0, 1]
    bins = np.linspace(0, 32, 33)
    for label, d in data.items():
        vals = np.concatenate([valid_clusters(d, det, 'ClusMaxSample')
                               for det in range(N_DET)])
        ax.hist(vals, bins=bins, histtype='step', lw=2,
                color=COLORS[label], density=True,
                label=f'{label}  (N={len(vals):,})')
    ax.axvline(CUT_SAMPLE_MIN, color='k', ls='--', lw=1.2,
               label=f'cut [{CUT_SAMPLE_MIN}, {CUT_SAMPLE_MAX}]')
    ax.axvline(CUT_SAMPLE_MAX, color='k', ls='--', lw=1.2)
    ax.set_xlabel('Peak sample index')
    ax.set_ylabel('Density')
    ax.set_title('ClusMaxSample (timing)')
    ax.legend(fontsize=8)

    # --- ClusTOT ---
    ax = axes[1, 0]
    bins = np.arange(-0.5, 32.5)
    for label, d in data.items():
        vals = np.concatenate([valid_clusters(d, det, 'ClusTOT')
                               for det in range(N_DET)])
        ax.hist(vals, bins=bins, histtype='step', lw=2,
                color=COLORS[label], density=True,
                label=f'{label}  (N={len(vals):,})')
    ax.axvline(CUT_TOT_MIN, color='k', ls='--', lw=1.2,
               label=f'TOTCut = {CUT_TOT_MIN}')
    ax.set_xlabel('Cluster TOT [samples]')
    ax.set_ylabel('Density')
    ax.set_title('Cluster TOT')
    ax.legend(fontsize=8)

    # --- ClusAmpl ---
    ax = axes[1, 1]
    for label, d in data.items():
        vals = np.concatenate([valid_clusters(d, det, 'ClusAmpl')
                               for det in range(N_DET)])
        # Use log-spaced bins to show the full range including negatives
        finite = vals[np.isfinite(vals)]
        bins = np.linspace(np.percentile(finite, 1),
                           np.percentile(finite, 99), 80)
        ax.hist(vals, bins=bins, histtype='step', lw=2,
                color=COLORS[label], density=True,
                label=f'{label}  median={np.median(finite):.0f}')
    ax.axvline(0, color='k', ls='-', lw=0.8, alpha=0.5)
    ax.set_xlabel('Cluster amplitude [ADC]')
    ax.set_ylabel('Density')
    ax.set_title('Cluster amplitude\n(negative = double pedestal subtraction)')
    ax.legend(fontsize=8)

    fig.tight_layout()
    out = os.path.join(out_dir, 'm3_cluster_properties.png')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    print(f'Saved → {out}')


# ---------------------------------------------------------------------------
# Plot 3: Per-plane cut survival rates
# ---------------------------------------------------------------------------

def plot_cut_survival(data, out_dir):
    """
    For each detector plane, show what fraction of clusters survive each cut.
    Directly identifies which cut is the bottleneck in ZS.
    """
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    fig.suptitle('Fraction of clusters surviving each cut — per detector plane',
                 fontsize=13)

    cut_labels = ['Size ≥ cut_min', f'MaxSample ∈ [{CUT_SAMPLE_MIN},{CUT_SAMPLE_MAX}]',
                  f'TOT > {CUT_TOT_MIN}', 'All cuts']

    for ax, (label, d) in zip(axes, data.items()):
        x = np.arange(N_DET)
        bar_data = {c: [] for c in cut_labels}

        for det in range(N_DET):
            size    = valid_clusters(d, det, 'ClusSize')
            sample  = valid_clusters(d, det, 'ClusMaxSample')
            tot     = valid_clusters(d, det, 'ClusTOT')
            n_total = len(size)
            if n_total == 0:
                for c in cut_labels:
                    bar_data[c].append(0.0)
                continue

            # Per-detector ClusSizeCut_Min
            size_cut = list(CUT_SIZE_MIN.values())[det]
            m_size   = size >= size_cut
            m_sample = (sample >= CUT_SAMPLE_MIN) & (sample <= CUT_SAMPLE_MAX)
            m_tot    = tot > CUT_TOT_MIN

            bar_data['Size ≥ cut_min'].append(np.mean(m_size) * 100)
            bar_data[f'MaxSample ∈ [{CUT_SAMPLE_MIN},{CUT_SAMPLE_MAX}]'].append(
                np.mean(m_sample) * 100)
            bar_data[f'TOT > {CUT_TOT_MIN}'].append(np.mean(m_tot) * 100)
            bar_data['All cuts'].append(np.mean(m_size & m_sample & m_tot) * 100)

        width = 0.2
        cut_colors = ['#4daf4a', '#ff7f00', '#984ea3', '#e41a1c']
        for i, (cut, vals) in enumerate(bar_data.items()):
            ax.bar(x + i * width, vals, width, label=cut,
                   color=cut_colors[i], alpha=0.85)

        ax.set_xticks(x + 1.5 * width)
        ax.set_xticklabels([f'mg{i}' for i in range(N_DET)], fontsize=8)
        ax.set_ylabel('% clusters surviving cut')
        ax.set_ylim(0, 110)
        ax.set_title(label, fontsize=11)
        ax.legend(fontsize=8)
        ax.axhline(100, color='k', ls=':', lw=0.8, alpha=0.4)

    fig.tight_layout()
    out = os.path.join(out_dir, 'm3_cut_survival.png')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    print(f'Saved → {out}')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    os.makedirs(PLOTS_DIR, exist_ok=True)

    data = {}
    for label, path in ANALYSE_FILES.items():
        print(f'\nLoading {label}: {path}')
        d = load(path)
        if d is None:
            sys.exit(1)
        print(f'  {d["n_events"]:,} events loaded')
        data[label] = d

    print_summary(data)
    plot_occupancy(data, PLOTS_DIR)
    plot_properties(data, PLOTS_DIR)
    plot_cut_survival(data, PLOTS_DIR)

    print('\nDone. Check plots/ for the diagnostic figures.')


if __name__ == '__main__':
    main()
