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
  - a stat-card + efficiency-vs-HV + resolution-vs-HV page, plus
  - the observations from the analysis session.

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


def build(cfg, df, csv, out):
    dp = df
    ipk = dp['eff_reco'].idxmax()
    peak_eff = dp.loc[ipk, 'eff_reco'] * 100
    peak_hv = dp.loc[ipk, 'hv']
    sig_pk = dp.loc[ipk, 'sigma_x_mm'] if np.isfinite(dp.loc[ipk, 'sigma_x_mm']) else float('nan')
    lo, hi = int(dp['hv'].min()), int(dp['hv'].max())
    # detect gaps in the nominal 5 V grid (e.g. the empty 365 V m3 file)
    grid = set(range(lo, hi + 1, 5))
    missing = sorted(grid - set(int(v) for v in dp['hv']))

    with PdfPages(out) as pdf:
        fig = plt.figure(figsize=(8.27, 11.69))
        gs = GridSpec(2, 1, figure=fig, height_ratios=[1, 1],
                      hspace=0.40, left=0.11, right=0.90, top=0.73, bottom=0.30)

        fig.text(0.06, 0.965, 'P2 BASKET — Detector 1', fontsize=27, fontweight='bold',
                 color=INK, va='top')
        fig.text(0.06, 0.923, 'Mesh-HV scan — efficiency & resolution', fontsize=14,
                 color='#555555', va='top')
        fig.text(0.06, 0.898,
                 f'{cfg.RUN}     drift = mesh + 180 V     connectors 1 & 10 dropped',
                 fontsize=9.5, color='#333333', va='top', family='monospace')

        # stat cards ------------------------------------------------------- #
        hax = fig.add_axes([0.06, 0.775, 0.88, 0.09]); hax.axis('off')
        hax.set_xlim(0, 1); hax.set_ylim(0, 1)

        def card(x, value, unit, label):
            sep = '' if unit == '%' else ' '
            hax.text(x, 0.66, f'{value}{sep}{unit}' if unit else f'{value}',
                     fontsize=19, fontweight='bold', color=INK, va='center')
            hax.text(x, 0.10, label, fontsize=8.0, color='dimgrey', va='center')
        card(0.00, f'{peak_eff:.1f}', '%', 'Peak efficiency (≤20 mm)')
        card(0.23, f'{peak_hv:.0f}', 'V', 'at mesh HV')
        card(0.42, f'{sig_pk:.1f}' if np.isfinite(sig_pk) else '—', 'mm', 'σₓ at peak')
        card(0.61, f'{lo}–{hi}', 'V', 'mesh HV range')
        card(0.90, f'{len(dp)}', '', 'HV points')

        # efficiency vs HV ------------------------------------------------- #
        axe = fig.add_subplot(gs[0])
        axe.errorbar(dp['hv'], dp['eff_reco'], yerr=dp['eff_reco_err'], fmt='o-',
                     color=C_RECO, capsize=4, lw=2, ms=7, label='reco within 20 mm')
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
        obs = (
            'OBSERVATIONS\n'
            f'• Mesh scanned {lo}→{hi} V in 5 V steps (drift = mesh + 180 V), one sub_run per point. '
            'Connectors 1 & 10 were disconnected during the scan and are dropped from the active area.\n'
            '• A single pad→M3 rigid transform is fit on the POOLED matched events across all HV '
            '(≈88.4°, scale 0.925 — consistent with the long run) and the transformed pad footprint is '
            'frozen, so every HV point shares an identical efficiency denominator region.\n'
            '• Per-sub_run HV spark veto is applied (mesh ch 1:0).\n'
            f'• Clean turn-on: reco ≤20 mm {e_lo:.2f}→{e_hi:.2f}, any-pad {a_lo:.2f}→{a_hi:.2f} over '
            f'{lo}→{hi} V, still rising at {hi} V. Core residual σ ~ 10–11 mm, flat (pad-pitch limited).'
            + miss + '\n'
            f'• Run-to-run caveat: the scan {hi} V point ({e_hi:.2f}) reads higher than the 6-30 long-run '
            '420 V (0.50) — same nominal voltage two days apart, i.e. conditioning / active-area '
            'definition, not chased here.'
        )
        fig.text(0.06, 0.265, obs, fontsize=9.2, color='#1a1a1a', va='top', ha='left',
                 wrap=True, linespacing=1.5,
                 bbox=dict(boxstyle='round,pad=0.8', fc='#f4f7fb', ec='#9db8d2', lw=1.0))

        fig.text(0.06, 0.03, f'source: {os.path.relpath(csv, qa.ANALYSIS_ROOT)}',
                 fontsize=7, color='#666666', family='monospace')
        fig.text(0.96, 0.008, datetime.date.today().isoformat(), ha='right',
                 fontsize=6, color='grey')
        pdf.savefig(fig, dpi=200)
        plt.close(fig)

    print(f'Wrote HV-scan PDF -> {out}')


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
