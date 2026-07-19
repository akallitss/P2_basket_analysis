#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_drift_scan_pdf.py

Compile the P2 (BASKET) det1 DRIFT-field scan into a styled PDF matching the
mesh-scan report (build_hv_scan_pdf.py). The drift scan holds the mesh at a
fixed working point and steps the drift (cathode) HV, so it probes the drift
region -- charge collection and, crucially, the TIME RESOLUTION -- at fixed
gain.

Reads efficiency_vs_drift<suffix>.csv written by 16_drift_scan_efficiency.py
(--scan drift) and REPLOTS from the CSV (vector, no blur):
  page 1  efficiency + estimated time resolution sigma_t vs drift HV
  page 2  charge collection (mean pad amplitude) + drift velocity
          (Magboltz nominal) with the measured peak-time spread

The time-resolution and drift-velocity panels are drawn only when the CSV
carries those columns (older scans still build what they have).

Usage:
  python3 build_drift_scan_pdf.py [run_key] [--out=PATH]   (default det1_driftscan2)
"""
import os
import sys
import glob
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
C_TRES = '#4b0082'        # indigo, matches stage-16 time-resolution plot
C_VD = '#2e8b57'         # seagreen, matches the Magboltz overlay
C_SPREAD = '#d62728'
PANEL_BG = '#f4f7fb'
PANEL_EC = '#9db8d2'


def load_drift_scan(cfg):
    """Return (df sorted by drift V, csv path) from 16_drift_scan_efficiency,
    or (None, None)."""
    suffix = cfg.product_suffix(veto_sparks=True)
    stage = os.path.join(cfg.OUT_BASE, '16_drift_scan_efficiency')
    csv = os.path.join(stage, f'efficiency_vs_drift{suffix}.csv')
    if not os.path.isfile(csv):
        cands = sorted(glob.glob(os.path.join(stage, 'efficiency_vs_drift*.csv')))
        if not cands:
            return None, None
        csv = cands[-1]
    df = pd.read_csv(csv).sort_values('x').reset_index(drop=True)
    return (df if not df.empty else None), csv


def efficient(df):
    """Rows with a real efficiency plateau value (drops discharge/no-track
    sub_runs like the det1 drift 765 V point at 9%)."""
    return df[df['eff_reco'] > 0.5]


def best_timing(df):
    """Row with the smallest sigma_t among efficient points, or None."""
    d = efficient(df)
    d = d[d['time_res_ns'].notna()] if 'time_res_ns' in d else d.iloc[0:0]
    return d.loc[d['time_res_ns'].idxmin()] if len(d) else None


def timing_window(df, col='time_res_ns', frac=1.10):
    """(vlo, vhi) of efficient points within `frac` of the min of `col`."""
    d = efficient(df)
    d = d[d[col].notna()]
    if not len(d):
        return None
    opt = d[d[col] <= frac * d[col].min()]
    return opt['x'].min(), opt['x'].max()


def _header(fig, subtitle, cfg, mesh_v):
    fig.text(0.06, 0.965, 'P2 BASKET — Detector 1', fontsize=27,
             fontweight='bold', color=INK, va='top')
    fig.text(0.06, 0.923, subtitle, fontsize=14, color='#555555', va='top')
    conn = ('connectors ' + ' & '.join(str(c) for c in cfg.DEAD_CONNECTORS)
            + ' dropped') if cfg.DEAD_CONNECTORS else 'all connectors live'
    mesh_lbl = f'mesh {mesh_v:.0f} V fixed' if mesh_v is not None else 'mesh fixed'
    fig.text(0.06, 0.898, f'{cfg.RUN}     {mesh_lbl}     {conn}',
             fontsize=9.5, color='#333333', va='top', family='monospace')


def _cards(fig, entries):
    hax = fig.add_axes([0.06, 0.775, 0.88, 0.09]); hax.axis('off')
    hax.set_xlim(0, 1); hax.set_ylim(0, 1)
    for x, value, unit, label in entries:
        sep = '' if unit in ('%', '×', '') else ' '
        hax.text(x, 0.66, f'{value}{sep}{unit}', fontsize=18.5,
                 fontweight='bold', color=INK, va='center')
        hax.text(x, 0.10, label, fontsize=7.8, color='dimgrey', va='center')


def _obs_box(fig, text):
    fig.text(0.06, 0.30, text, fontsize=9.2, color='#1a1a1a', va='top',
             ha='left', wrap=True, linespacing=1.5,
             bbox=dict(boxstyle='round,pad=0.8', fc=PANEL_BG, ec=PANEL_EC, lw=1.0))


def _footer(fig, csv):
    fig.text(0.06, 0.03, f'source: {os.path.relpath(csv, qa.ANALYSIS_ROOT)}',
             fontsize=7, color='#666666', family='monospace')
    fig.text(0.96, 0.008, datetime.date.today().isoformat(), ha='right',
             fontsize=6, color='grey')


def _page_eff_timing(pdf, cfg, df, csv, mesh_v):
    dp = df
    eff_pts = efficient(dp)
    peak_eff = eff_pts['eff_reco'].max() * 100
    lo, hi = int(dp['x'].min()), int(dp['x'].max())
    has_tres = 'time_res_ns' in dp and dp['time_res_ns'].notna().any()
    best = best_timing(dp) if has_tres else None

    fig = plt.figure(figsize=(8.27, 11.69))
    gs = GridSpec(2, 1, figure=fig, height_ratios=[1, 1],
                  hspace=0.42, left=0.11, right=0.94, top=0.73, bottom=0.36)
    _header(fig, 'Drift-field scan — efficiency & time resolution', cfg, mesh_v)

    if best is not None:
        cards = [(0.00, f'{peak_eff:.1f}', '%', 'Plateau efficiency'),
                 (0.22, f'{best.time_res_ns:.1f}', 'ns', 'Best σ$_t$'),
                 (0.42, f'{best.x:.0f}', 'V', 'at drift HV'),
                 (0.62, f'{lo}–{hi}', 'V', 'drift HV range'),
                 (0.90, f'{len(dp)}', '', 'drift points')]
    else:
        cards = [(0.00, f'{peak_eff:.1f}', '%', 'Plateau efficiency'),
                 (0.34, f'{lo}–{hi}', 'V', 'drift HV range'),
                 (0.90, f'{len(dp)}', '', 'drift points')]
    _cards(fig, cards)

    # efficiency vs drift
    axe = fig.add_subplot(gs[0])
    axe.errorbar(dp['x'], dp['eff_reco'], yerr=dp['eff_reco_err'], fmt='o-',
                 color=C_RECO, capsize=4, lw=2, ms=7,
                 label=f'reco within {cfg.MATCH_R:g} mm')
    axe.plot(dp['x'], dp['eff_anyhit'], 's--', color=C_ANY, ms=6, alpha=0.85,
             label='any pad fired')
    axe.set_xlabel('drift HV [V]')
    axe.set_ylabel('efficiency (frozen active area)')
    axe.set_ylim(0, 1.02); axe.grid(True, alpha=0.3); axe.legend(fontsize=9)
    axe.set_title('Efficiency vs drift HV', fontsize=12)

    # time resolution vs drift (zoomed to the minimum)
    axt = fig.add_subplot(gs[1])
    if has_tres:
        d = dp[dp['time_res_ns'].notna() & (dp['e_drift_Vcm'] > 50)]
        smin = d['time_res_ns'].min()
        z = d[d['time_res_ns'] <= 2.5 * smin]        # zoom to the valley
        axt.plot(z['x'], z['time_res_ns'], 'o-', color=C_TRES, lw=2, ms=7,
                 label='σ$_t$ (leading-pad, robust)')
        win = timing_window(dp, 'time_res_ns')
        if win:
            axt.axvspan(win[0] - 12, win[1] + 12, color='gold', alpha=0.22,
                        label='best-timing window (efficient)')
        if best is not None:
            axt.annotate(f'{best.time_res_ns:.1f} ns',
                         xy=(best.x, best.time_res_ns),
                         xytext=(best.x, best.time_res_ns + 0.45 *
                                 (z['time_res_ns'].max() - smin + 1)),
                         ha='center', fontsize=9, color=C_TRES,
                         arrowprops=dict(arrowstyle='->', color=C_TRES))
        loy, hiy = z['time_res_ns'].min(), z['time_res_ns'].max()
        axt.set_ylim(max(0, loy - 0.2 * (hiy - loy) - 0.5),
                     hiy + 0.55 * (hiy - loy) + 0.5)
        axt.legend(fontsize=9, loc='upper left')
        axt.set_title('Estimated time resolution vs drift HV '
                      '(zoom; ~0 V point dropped)', fontsize=12)
    else:
        axt.text(0.5, 0.5, 'time resolution not in this CSV\n'
                 '(re-run 16_drift_scan_efficiency.py)', ha='center',
                 va='center', transform=axt.transAxes, color='grey')
        axt.set_title('Estimated time resolution vs drift HV', fontsize=12)
    axt.set_xlabel('drift HV [V]')
    axt.set_ylabel('time resolution σ$_t$ [ns]')
    axt.grid(True, alpha=0.3)

    # observations
    e_lo = dp['eff_reco'].iloc[0]
    win = timing_window(dp, 'time_res_ns') if has_tres else None
    win_txt = (f'{win[0]:.0f}–{win[1]:.0f} V' if win else '—')
    tline = (
        f'• Estimated time resolution (per-event leading-pad time_of_max vs the '
        f'scintillator trigger, fine-time-step corrected) bottoms out at '
        f'σ$_t$ = {best.time_res_ns:.1f} ns at drift {best.x:.0f} V '
        f'(E ≈ {best.e_drift_Vcm:.0f} V/cm), at {best.eff_reco:.0%} efficiency. '
        f'The efficient best-timing window is {win_txt}.\n'
        if best is not None else '')
    obs = (
        'OBSERVATIONS\n'
        f'• Mesh held at {mesh_v:.0f} V (fixed gain); the drift/cathode HV is '
        f'stepped {lo}→{hi} V, so this scan probes the DRIFT region only.\n'
        f'• Efficiency turns on from {e_lo:.0%} (near-zero drift field) to the '
        f'{peak_eff:.0f}% plateau, reached by drift ≈ 565 V (drift gap ≈ 150 V).\n'
        + tline +
        '• σ$_t$ rises on BOTH sides of the optimum: low field is limited by '
        'diffusion and poor charge collection, high field by the slower drift '
        'velocity. See page 2 for charge collection and the drift velocity.'
    )
    _obs_box(fig, obs)
    _footer(fig, csv)
    pdf.savefig(fig, dpi=200)
    plt.close(fig)


def _page_transport(pdf, cfg, df, csv, mesh_v):
    dp = df
    has_amp = 'amp_mean' in dp and dp['amp_mean'].notna().any()
    has_vd = 'vd_magboltz_um_ns' in dp and dp['vd_magboltz_um_ns'].notna().any()
    has_spread = 'tom_spread_ns' in dp and dp['tom_spread_ns'].notna().any()
    if not (has_amp or has_vd or has_spread):
        return

    fig = plt.figure(figsize=(8.27, 11.69))
    gs = GridSpec(2, 1, figure=fig, height_ratios=[1, 1],
                  hspace=0.42, left=0.11, right=0.90, top=0.73, bottom=0.37)
    _header(fig, 'Drift-field scan — charge collection & drift velocity',
            cfg, mesh_v)

    best = best_timing(dp)
    cards = []
    if has_amp:
        amp_plateau = efficient(dp)['amp_mean'].median()
        cards.append((0.00, f'{amp_plateau:.0f}', 'ADC',
                      'plateau mean amp'))
    if has_vd and best is not None and np.isfinite(best.vd_magboltz_um_ns):
        cards.append((0.30, f'{best.vd_magboltz_um_ns:.0f}', 'µm/ns',
                      'Magboltz v$_d$ at optimum'))
    if has_spread:
        smin = efficient(dp)['tom_spread_ns'].min()
        cards.append((0.66, f'{smin:.0f}', 'ns', 'min peak-time spread'))
    if cards:
        _cards(fig, cards)

    # top: charge collection (mean pad amplitude, gain fixed)
    axa = fig.add_subplot(gs[0])
    if has_amp:
        yerr = dp['amp_mean_err'] if 'amp_mean_err' in dp else None
        axa.errorbar(dp['x'], dp['amp_mean'], yerr=yerr, fmt='o-', color=C_RECO,
                     capsize=4, lw=2, ms=7, label='mean pad amp')
        if 'amp_median' in dp:
            axa.plot(dp['x'], dp['amp_median'], 's--', color=C_ANY, ms=6,
                     alpha=0.85, label='median pad amp')
        axa.legend(fontsize=9)
    axa.set_xlabel('drift HV [V]'); axa.set_ylabel('pad amplitude [ADC]')
    axa.set_ylim(0, None); axa.grid(True, alpha=0.3)
    axa.set_title('Charge collection (mean pad amplitude, gain fixed)',
                  fontsize=12)

    # bottom: Magboltz drift velocity + measured peak-time spread
    axv = fig.add_subplot(gs[1])
    gap_cm = 0.3
    if has_vd:
        mesh = dp['mesh'].iloc[0]
        # smooth Magboltz curve over the tabulated field range only
        vv = np.linspace(dp['x'].min(), dp['x'].max(), 400)
        vd = np.interp((vv - mesh) / gap_cm,
                       (dp['e_drift_Vcm']).to_numpy(),
                       dp['vd_magboltz_um_ns'].to_numpy(),
                       left=np.nan, right=np.nan)
        axv.plot(vv, vd, '-', color=C_VD, lw=2, label='Magboltz v$_d$ (nominal)')
        axv.plot(dp['x'], dp['vd_magboltz_um_ns'], 'o', color=C_VD, ms=5)
    axv.set_xlabel('drift HV [V]')
    axv.set_ylabel('drift velocity v$_d$ [µm/ns]', color=C_VD)
    axv.tick_params(axis='y', labelcolor=C_VD); axv.set_ylim(0, None)
    axv.grid(True, alpha=0.3)
    axs = axv.twinx()
    if has_spread:
        axs.plot(dp['x'], dp['tom_spread_ns'], 's--', color=C_SPREAD, ms=6,
                 alpha=0.85, label='peak-time spread (data)')
        axs.set_ylabel('peak-time spread [ns]', color=C_SPREAD)
        axs.tick_params(axis='y', labelcolor=C_SPREAD); axs.set_ylim(0, None)
    h1, l1 = axv.get_legend_handles_labels()
    h2, l2 = axs.get_legend_handles_labels()
    axv.legend(h1 + h2, l1 + l2, fontsize=8, loc='center right')
    axv.set_title('Drift velocity (nominal Magboltz) vs measured spread',
                  fontsize=12)

    obs = (
        'CHARGE-TRANSPORT NOTES\n'
        '• Gain is fixed (mesh constant), so the mean pad amplitude tracks '
        'charge COLLECTION through the drift gap, not gas gain: it rises out of '
        'the low-field region and is roughly flat once the field fully collects '
        'the ionisation.\n'
        '• The nominal-gas Magboltz v$_d$(E) is shown only over its tabulated '
        'range (v$_d$ → 0 as E → 0 below it, not drawn). It falls monotonically '
        'above ≈ 200 V/cm, so the high-drift timing degradation on page 1 is '
        'consistent with the slower nominal drift velocity.\n'
        '• The peak-time spread (all signal hits) is a drift-time proxy; unlike '
        'the per-event σ$_t$ on page 1 it also carries the intra-track depth '
        'spread. A clean absolute v$_d$ measurement needs the full timing '
        'chain (trigger latency + fine-time-step), not this proxy.'
    )
    _obs_box(fig, obs)
    _footer(fig, csv)
    pdf.savefig(fig, dpi=200)
    plt.close(fig)


def build(cfg, df, csv, out):
    mesh_v = float(df['mesh'].iloc[0]) if 'mesh' in df else None
    with PdfPages(out) as pdf:
        _page_eff_timing(pdf, cfg, df, csv, mesh_v)
        _page_transport(pdf, cfg, df, csv, mesh_v)
    print(f'Wrote drift-scan PDF -> {out}')


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('-')]
    key = args[0] if args else 'det1_driftscan2'
    cfg = qa.get_config(key)
    print(cfg)
    df, csv = load_drift_scan(cfg)
    if df is None:
        print('No efficiency_vs_drift CSV found — run '
              '16_drift_scan_efficiency.py --scan drift first.')
        return
    print(f'  {len(df)} drift points, {int(df["x"].min())}–{int(df["x"].max())} V'
          f'  ({csv})')
    default_out = os.path.join(qa.ANALYSIS_ROOT, cfg.DET_TAG,
                               f'p2_{cfg.DET_TAG}_drift_scan.pdf')
    out = next((a.split('=', 1)[1] for a in sys.argv if a.startswith('--out=')),
               default_out)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    build(cfg, df, csv, out)


if __name__ == '__main__':
    main()
