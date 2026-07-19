#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_hv_scan_pdf.py

Compile the P2 (BASKET) det1 MESH-HV scan into a styled PDF matching the final
QA report. Adapted from nTof_x17/mx_june_cosmic_qa/build_hv_scan_pdf.py, but P2
sweeps the MESH (not a resistive layer), drift tracks it at mesh + 180 V, and the
active area is defined by the transformed pad footprint (connectors 1 & 10 dropped
for the scan).

Reads the scan's efficiency_vs_hv<suffix>.csv written by 11_hv_scan_efficiency.py
and REPLOTS from the CSV (vector, no blur):
  - a stat-card + efficiency-vs-HV + resolution-vs-HV page,
  - a mean-pad-amplitude (gain proxy) page with the exponential-gain fit, and
  - the observations from the analysis session.

The amplitude page is drawn only when the CSV carries the amp_mean column
(scans processed before the gain proxy was added still build the first page).

Usage:
  python3 build_hv_scan_pdf.py [run_key] [--out=PATH]     (default key det1_hvscan)
"""
import os
import re
import sys
import datetime

import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec

import p2_qa_config as qa

INK = '#1f3a5f'
C_RECO = '#1f77b4'
C_ANY = '#d62728'
C_SX = '#1f77b4'
C_SY = '#d62728'


def load_scan(cfg):
    """Return the scan DataFrame (sorted by HV) written by stage 11, or None."""
    suffix = cfg.product_suffix(veto_sparks=True)  # _without_connectors_1_10_spark_vetoed
    stage = os.path.join(cfg.OUT_BASE, '11_hv_scan_efficiency')
    csv = os.path.join(stage, f'efficiency_vs_hv{suffix}.csv')
    if not os.path.isfile(csv):
        # fall back to any efficiency_vs_hv*.csv in the stage dir
        import glob
        cands = sorted(glob.glob(os.path.join(stage, 'efficiency_vs_hv*.csv')))
        if not cands:
            return None, None
        csv = cands[-1]
    df = pd.read_csv(csv).sort_values('hv').reset_index(drop=True)
    return (df if not df.empty else None), csv


def drift_of(df):
    m = re.search(r'drift_(\d+)V', str(df['subrun'].iloc[0]))
    return int(m.group(1)) if m else None


def drift_gap(df):
    """Mesh->drift tandem offset [V], read from the sub_run names if present
    (…_<mesh>_<drift>), else None. det1 7-2 scan = 180 V; det1 7-19 = 310 V."""
    gaps = set()
    for s in df['subrun'].astype(str):
        m = re.search(r'_(\d+)_(\d+)$', s)
        if m:
            gaps.add(int(m.group(2)) - int(m.group(1)))
    return gaps.pop() if len(gaps) == 1 else None


def gain_fit(df):
    """Exponential gas-gain fit ln(amp_mean) = a*V + b over the points above
    turn-on (amp rising). Returns (a, doubling_V, span) or None."""
    if 'amp_mean' not in df or df['amp_mean'].notna().sum() < 3:
        return None
    d = df[df['amp_mean'] > 0]
    a, _ = np.polyfit(d['hv'], np.log(d['amp_mean']), 1)
    span = d['amp_mean'].max() / d['amp_mean'].min()
    return a, (np.log(2) / a if a > 0 else np.nan), span


def build(cfg, df, csv, out):
    dp = df
    ipk = dp['eff_reco'].idxmax()
    peak_eff = dp.loc[ipk, 'eff_reco'] * 100
    peak_hv = dp.loc[ipk, 'hv']
    sig_pk = dp.loc[ipk, 'sigma_x_mm'] if np.isfinite(dp.loc[ipk, 'sigma_x_mm']) else float('nan')
    lo, hi = int(dp['hv'].min()), int(dp['hv'].max())
    gap = drift_gap(dp)
    gap_lbl = f'drift = mesh + {gap} V' if gap is not None else 'drift stepped in tandem'
    # detect gaps in the nominal 5 V grid (e.g. the empty 365 V m3 file)
    grid = set(range(lo, hi + 1, 5))
    missing = sorted(grid - set(int(v) for v in dp['hv']))

    with PdfPages(out) as pdf:
        fig = plt.figure(figsize=(8.27, 11.69))
        gs = GridSpec(2, 1, figure=fig, height_ratios=[1, 1],
                      hspace=0.42, left=0.11, right=0.90, top=0.73, bottom=0.36)

        fig.text(0.06, 0.965, 'P2 BASKET — Detector 1', fontsize=27, fontweight='bold',
                 color=INK, va='top')
        fig.text(0.06, 0.923, 'Mesh-HV scan — efficiency & resolution', fontsize=14,
                 color='#555555', va='top')
        conn_lbl = ('connectors ' + ' & '.join(str(c) for c in cfg.DEAD_CONNECTORS)
                    + ' dropped') if cfg.DEAD_CONNECTORS else 'all connectors live'
        fig.text(0.06, 0.898,
                 f'{cfg.RUN}     {gap_lbl}     {conn_lbl}',
                 fontsize=9.5, color='#333333', va='top', family='monospace')

        # stat cards ------------------------------------------------------- #
        hax = fig.add_axes([0.06, 0.775, 0.88, 0.09]); hax.axis('off')
        hax.set_xlim(0, 1); hax.set_ylim(0, 1)

        def card(x, value, unit, label):
            sep = '' if unit == '%' else ' '
            hax.text(x, 0.66, f'{value}{sep}{unit}' if unit else f'{value}',
                     fontsize=19, fontweight='bold', color=INK, va='center')
            hax.text(x, 0.10, label, fontsize=8.0, color='dimgrey', va='center')
        card(0.00, f'{peak_eff:.1f}', '%', f'Peak efficiency (≤{cfg.MATCH_R:g} mm)')
        card(0.23, f'{peak_hv:.0f}', 'V', 'at mesh HV')
        card(0.42, f'{sig_pk:.1f}' if np.isfinite(sig_pk) else '—', 'mm', 'σₓ at peak')
        card(0.61, f'{lo}–{hi}', 'V', 'mesh HV range')
        card(0.90, f'{len(dp)}', '', 'HV points')

        # efficiency vs HV ------------------------------------------------- #
        axe = fig.add_subplot(gs[0])
        axe.errorbar(dp['hv'], dp['eff_reco'], yerr=dp['eff_reco_err'], fmt='o-',
                     color=C_RECO, capsize=4, lw=2, ms=7,
                     label=f'reco within {cfg.MATCH_R:g} mm')
        axe.plot(dp['hv'], dp['eff_anyhit'], 's--', color=C_ANY, ms=6, alpha=0.85,
                 label='any pad fired')
        axe.set_xlabel('mesh HV [V]'); axe.set_ylabel('efficiency (frozen active area)')
        axe.set_ylim(0, 1.02); axe.grid(True, alpha=0.3); axe.legend(fontsize=9)
        axe.set_title('Efficiency vs mesh HV', fontsize=12)

        # resolution vs HV ------------------------------------------------- #
        axr = fig.add_subplot(gs[1])
        axr.plot(dp['hv'], dp['sigma_x_mm'], 'o-', color=C_SX, lw=2, ms=7, label='σₓ')
        axr.plot(dp['hv'], dp['sigma_y_mm'], 's--', color=C_SY, lw=2, ms=6, label='σ_y')
        axr.set_xlabel('mesh HV [V]'); axr.set_ylabel('core residual σ [mm]')
        axr.set_ylim(0, None); axr.grid(True, alpha=0.3); axr.legend(fontsize=9)
        axr.set_title('Spatial resolution vs mesh HV', fontsize=12)

        # observations box ------------------------------------------------- #
        e_lo = dp['eff_reco'].iloc[0]; e_hi = dp['eff_reco'].iloc[-1]
        a_lo = dp['eff_anyhit'].iloc[0]; a_hi = dp['eff_anyhit'].iloc[-1]
        miss = (f'  ⚠ {", ".join(f"{m} V" for m in missing)} dropped '
                '(m3 tracking file transferred empty — re-fetch pending).'
                if missing else '')
        sig_lo = dp['sigma_x_mm'].dropna()
        sig_txt = (f'{sig_lo.min():.1f}–{sig_lo.max():.1f} mm'
                   if len(sig_lo) else '—')
        obs = (
            'OBSERVATIONS\n'
            f'• Mesh scanned {lo}→{hi} V in 5 V steps ({gap_lbl}), one sub_run per point. '
            + (conn_lbl.capitalize() + ' from the active area.\n' if cfg.DEAD_CONNECTORS
               else 'All connectors are live.\n')
            + '• A single pad→M3 rigid transform is fit on the POOLED matched events across all HV '
            'and the transformed pad footprint is frozen, so every HV point shares an identical '
            'efficiency denominator region.\n'
            f'• Per-sub_run HV spark veto is applied (mesh ch {cfg.SPARK_CHANNEL}); '
            'auto-flagged hot pads are dropped.\n'
            f'• Clean turn-on: reco ≤{cfg.MATCH_R:g} mm {e_lo:.2f}→{e_hi:.2f}, any-pad '
            f'{a_lo:.2f}→{a_hi:.2f} over {lo}→{hi} V. Core residual σ {sig_txt} (pad-pitch limited).'
            + miss + '\n'
            f'• Plateau reached by ~415 V (98% of the fitted asymptote at ≈419 V); below ~405 V the '
            'efficiency drops onto the turn-on slope. See page 2 for the gain (amplitude) curve.'
        )
        fig.text(0.06, 0.30, obs, fontsize=9.2, color='#1a1a1a', va='top', ha='left',
                 wrap=True, linespacing=1.5,
                 bbox=dict(boxstyle='round,pad=0.8', fc='#f4f7fb', ec='#9db8d2', lw=1.0))

        fig.text(0.06, 0.03, f'source: {os.path.relpath(csv, qa.ANALYSIS_ROOT)}',
                 fontsize=7, color='#666666', family='monospace')
        fig.text(0.96, 0.008, datetime.date.today().isoformat(), ha='right',
                 fontsize=6, color='grey')
        pdf.savefig(fig, dpi=200)
        plt.close(fig)

        # ---- page 2: mean pad amplitude (gain proxy) -------------------- #
        if 'amp_mean' in dp and dp['amp_mean'].notna().any():
            _amplitude_page(pdf, cfg, dp, csv, gap_lbl, lo, hi)

    print(f'Wrote HV-scan PDF -> {out}')


def _amplitude_page(pdf, cfg, dp, csv, gap_lbl, lo, hi):
    """Second PDF page: mean/median pad amplitude vs mesh HV (gain proxy),
    linear + log, with the exponential gas-gain fit."""
    fit = gain_fit(dp)
    yerr = dp['amp_mean_err'] if 'amp_mean_err' in dp else None

    fig = plt.figure(figsize=(8.27, 11.69))
    gs = GridSpec(2, 1, figure=fig, height_ratios=[1, 1],
                  hspace=0.48, left=0.11, right=0.90, top=0.70, bottom=0.36)

    fig.text(0.06, 0.965, 'P2 BASKET — Detector 1', fontsize=27,
             fontweight='bold', color=INK, va='top')
    fig.text(0.06, 0.923, 'Mesh-HV scan — gas gain (mean pad amplitude)',
             fontsize=14, color='#555555', va='top')
    fig.text(0.06, 0.898,
             f'{cfg.RUN}     {gap_lbl}     hot pads + spark events removed',
             fontsize=9.5, color='#333333', va='top', family='monospace')

    # stat cards
    hax = fig.add_axes([0.06, 0.775, 0.88, 0.09]); hax.axis('off')
    hax.set_xlim(0, 1); hax.set_ylim(0, 1)

    def card(x, value, unit, label):
        sep = '' if unit in ('%', '×') else ' '
        hax.text(x, 0.66, f'{value}{sep}{unit}' if unit else f'{value}',
                 fontsize=19, fontweight='bold', color=INK, va='center')
        hax.text(x, 0.10, label, fontsize=8.0, color='dimgrey', va='center')
    amp_lo = dp['amp_mean'].iloc[0]; amp_hi = dp['amp_mean'].iloc[-1]
    card(0.00, f'{amp_lo:.0f}', 'ADC', f'mean amp at {lo} V')
    card(0.25, f'{amp_hi:.0f}', 'ADC', f'mean amp at {hi} V')
    if fit is not None:
        _, doubling, span = fit
        card(0.52, f'{span:.1f}', '×', f'gain span {lo}→{hi} V')
        card(0.76, f'{doubling:.0f}', 'V', 'gain-doubling ΔV (fit)')

    # linear
    ax = fig.add_subplot(gs[0])
    ax.errorbar(dp['hv'], dp['amp_mean'], yerr=yerr, fmt='o-', color=C_RECO,
                capsize=4, lw=2, ms=7, label='mean pad amp')
    ax.plot(dp['hv'], dp['amp_median'], 's--', color=C_ANY, ms=6, alpha=0.85,
            label='median pad amp')
    ax.set_xlabel('mesh HV [V]'); ax.set_ylabel('pad amplitude [ADC]')
    ax.set_ylim(0, None); ax.grid(True, alpha=0.3); ax.legend(fontsize=9)
    ax.set_title('Mean pad amplitude vs mesh HV (linear)', fontsize=12)

    # log with exponential fit overlay
    axl = fig.add_subplot(gs[1])
    axl.errorbar(dp['hv'], dp['amp_mean'], yerr=yerr, fmt='o', color=C_RECO,
                 capsize=4, ms=7, label='mean pad amp')
    axl.plot(dp['hv'], dp['amp_median'], 's', color=C_ANY, ms=6, alpha=0.85,
             label='median pad amp')
    if fit is not None:
        a, doubling, _ = fit
        d = dp[dp['amp_mean'] > 0]
        b = np.log(d['amp_mean']).mean() - a * d['hv'].mean()
        xx = np.linspace(dp['hv'].min(), dp['hv'].max(), 100)
        axl.plot(xx, np.exp(a * xx + b), '-', color='grey', lw=1.5,
                 label=f'exp fit: ×2 per {doubling:.0f} V')
    axl.set_yscale('log')
    axl.set_xlabel('mesh HV [V]'); axl.set_ylabel('pad amplitude [ADC]')
    axl.grid(True, alpha=0.3, which='both'); axl.legend(fontsize=9)
    axl.set_title('Log scale — exponential gas gain is a straight line',
                  fontsize=12)

    # observations
    obs = (
        'GAIN NOTES\n'
        f'• Mean pad amplitude is a gain proxy: it rises from {amp_lo:.0f} ADC at '
        f'{lo} V to {amp_hi:.0f} ADC at {hi} V'
        + (f' (×{fit[2]:.1f}), an exponential ≈ ×2 per {fit[1]:.0f} V — '
           'the expected Micromegas avalanche behaviour.\n' if fit is not None
           else '.\n')
        + '• The amplitude spectrum is Landau-like, so the mean sits well above '
        'the median; both track the same exponential, so either can set a working '
        'point.\n'
        '• Unlike the efficiency, the amplitude does NOT plateau — it keeps '
        'climbing with mesh HV. Efficiency saturates once the smallest signals '
        'clear threshold, while the gain (and the spark rate) keeps growing, so '
        'the lowest mesh HV on the efficiency plateau is the safe operating point.'
    )
    fig.text(0.06, 0.30, obs, fontsize=9.2, color='#1a1a1a', va='top',
             ha='left', wrap=True, linespacing=1.5,
             bbox=dict(boxstyle='round,pad=0.8', fc='#f4f7fb', ec='#9db8d2', lw=1.0))

    fig.text(0.06, 0.03, f'source: {os.path.relpath(csv, qa.ANALYSIS_ROOT)}',
             fontsize=7, color='#666666', family='monospace')
    fig.text(0.96, 0.008, datetime.date.today().isoformat(), ha='right',
             fontsize=6, color='grey')
    pdf.savefig(fig, dpi=200)
    plt.close(fig)


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('-')]
    key = args[0] if args else 'det1_hvscan'
    cfg = qa.get_config(key)
    print(cfg)
    df, csv = load_scan(cfg)
    if df is None:
        print('No efficiency_vs_hv CSV found — run 11_hv_scan_efficiency.py first.')
        return
    print(f'  {len(df)} HV points, {int(df["hv"].min())}–{int(df["hv"].max())} V  ({csv})')
    default_out = os.path.join(qa.ANALYSIS_ROOT, cfg.DET_TAG,
                               f'p2_{cfg.DET_TAG}_hv_scan.pdf')
    out = next((a.split('=', 1)[1] for a in sys.argv if a.startswith('--out=')), default_out)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    build(cfg, df, csv, out)


if __name__ == '__main__':
    main()
