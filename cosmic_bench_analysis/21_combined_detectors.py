#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
21_combined_detectors.py

Conference figures: all four detectors on one set of axes.

  combined_efficiency.png   efficiency per detector (run average AND the
                            charged-up late window, side by side)
  combined_mesh_scan.png    efficiency vs mesh HV, every detector that has one
  combined_drift_scan.png   efficiency + time resolution vs drift field
  combined_summary.csv      the numbers behind all three

Two things this figure set is careful about, because both would mislead:

  * A run-integrated efficiency is NOT the detector's efficiency when the run
    contains charging-up. det4 reads 88.3 % integrated and 92.3 % charged-up.
    Both are plotted, so the difference is visible rather than hidden by a
    choice.
  * A detector whose efficiency FALLS during the run (det3) has no meaningful
    single number at all. It is drawn hatched and labelled, not quietly
    averaged in with the others.

Usage:
  python3 21_combined_detectors.py [-o OUTDIR]
"""

import argparse
import glob
import os
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import p2_qa_config as qa

# Standard run per detector (see the per-detector DET*_EFFICIENCY.md).
STD = {
    'det1': ('det1_long5', 'det1_meshscan1', 'det1_driftscan2'),
    'det2': ('det2_long1', 'det2_hvscan1', None),
    'det3': ('det3_initial1', None, 'det3_driftscan1'),
    'det4': ('det4_long2', None, 'det4_driftscan1'),
}
MIN_EV_TIMING = 500      # below this a sigma_t is statistics, not resolution
COL = {'det1': '#1f77b4', 'det2': '#2ca02c', 'det3': '#d62728', 'det4': '#ff7f0e'}


def _find(cfg, stage, stem):
    pats = sorted(glob.glob(os.path.join(cfg.OUT_BASE, stage,
                                         f'{stem}*spark_vetoed.csv')))
    pats = [p for p in pats if 'last' not in os.path.basename(p)]
    return pats[-1] if pats else None


def charging(cfg):
    """(full, late, verdict) from stage 20's summary, if it has been run."""
    f = sorted(glob.glob(os.path.join(cfg.OUT_BASE, '20_charging_up',
                                      'charging_summary*.txt')))
    if not f:
        return np.nan, np.nan, ''
    txt = open(f[-1]).read()
    def g(pat):
        m = re.search(pat, txt)
        return float(m.group(1)) if m else np.nan
    full = g(r'FULL RUN \(as quoted\)\s*: eff ([\d.]+)')
    late = g(r'late\s+\([^)]*\)\s*: eff ([\d.]+)')
    ver = 'DEGRADING' if 'DEGRADING' in txt else (
        'RISING' if 'STILL RISING' in txt else 'PLATEAUED')
    return full, late, ver


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-o', '--out', default=None)
    args = ap.parse_args()
    out = args.out or os.path.join(qa.ANALYSIS_ROOT, 'combined')
    os.makedirs(out, exist_ok=True)

    rows = []
    for det, (lk, mk, dk) in STD.items():
        cfg = qa.get_config(lk)
        full, late, ver = charging(cfg)
        rows.append(dict(det=det, run=cfg.RUN, sub_run=cfg.SUB_RUN,
                         eff_full=full, eff_late=late, verdict=ver))
    summ = pd.DataFrame(rows)
    summ.to_csv(os.path.join(out, 'combined_summary.csv'), index=False)
    print(summ.to_string(index=False))

    # ---------------------------------------------------------- efficiency #
    fig, ax = plt.subplots(figsize=(9, 5.5))
    x = np.arange(len(summ))
    w = 0.38
    for i, r in summ.iterrows():
        deg = r['verdict'] == 'DEGRADING'
        ax.bar(x[i] - w / 2, r['eff_full'], w, color=COL[r['det']], alpha=.45,
               hatch='//' if deg else None,
               label='run average' if i == 0 else None)
        ax.bar(x[i] + w / 2, r['eff_late'], w, color=COL[r['det']],
               hatch='//' if deg else None,
               label='charged-up (late window)' if i == 0 else None)
        for xx, v in ((x[i] - w / 2, r['eff_full']), (x[i] + w / 2, r['eff_late'])):
            if np.isfinite(v):
                ax.text(xx, v + 0.8, f'{v:.1f}', ha='center', fontsize=8.5)
        if deg:
            ax.text(x[i], 40, 'DEGRADING\nno single value\nis meaningful',
                    ha='center', fontsize=8, color='crimson', weight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels([f"{r['det']}\n{r['sub_run'][:22]}" for _, r in summ.iterrows()],
                       fontsize=8)
    ax.set_ylabel('efficiency (reco within R) [%]')
    ax.set_ylim(0, 105); ax.grid(axis='y', alpha=.3); ax.legend(fontsize=9)
    ax.set_title('P2 cosmic-bench efficiency — all detectors\n'
                 'hatched = efficiency fell during the run, so neither bar is '
                 '"the" efficiency', fontsize=10)
    fig.tight_layout()
    p = os.path.join(out, 'combined_efficiency.png')
    fig.savefig(p, dpi=200, bbox_inches='tight'); plt.close(fig)
    print('saved', p)

    # ----------------------------------------------------------- mesh scan #
    fig, ax = plt.subplots(figsize=(9, 5.5))
    n = 0
    for det, (lk, mk, dk) in STD.items():
        if not mk:
            continue
        cfg = qa.get_config(mk)
        f = _find(cfg, '11_hv_scan_efficiency', 'efficiency_vs_hv')
        if not f:
            continue
        d = pd.read_csv(f).sort_values('hv')
        ax.errorbar(d['hv'], 100 * d['eff_reco'],
                    yerr=100 * d.get('eff_reco_err', 0), fmt='o-',
                    color=COL[det], label=f'{det} ({mk})', capsize=2)
        n += 1
    ax.set_xlabel('mesh HV [V]'); ax.set_ylabel('efficiency [%]')
    ax.grid(alpha=.3); ax.legend(fontsize=9)
    ax.set_title('Mesh (gain) scans — drift stepped in tandem\n'
                 'det3 excluded: its scan ran while the detector was collapsing',
                 fontsize=10)
    fig.tight_layout()
    p = os.path.join(out, 'combined_mesh_scan.png')
    fig.savefig(p, dpi=200, bbox_inches='tight'); plt.close(fig)
    print(f'saved {p} ({n} detectors)')

    # ---------------------------------------------------------- drift scan #
    fig, axs = plt.subplots(1, 2, figsize=(14, 5.4))
    n = 0
    for det, (lk, mk, dk) in STD.items():
        if not dk:
            continue
        cfg = qa.get_config(dk)
        f = _find(cfg, '16_drift_scan_efficiency', 'efficiency_vs_drift')
        if not f:
            continue
        d = pd.read_csv(f).sort_values('e_drift_Vcm')
        axs[0].errorbar(d['e_drift_Vcm'], 100 * d['eff_reco'],
                        yerr=100 * d.get('eff_reco_err', 0), fmt='o-',
                        color=COL[det], label=f'{det} ({dk})', capsize=2)
        if 'time_res_ns' in d.columns:
            # A sigma_t from a handful of events is not a resolution. det4's
            # drift scan yields ~90 events per point (no drift field, so almost
            # nothing is detected) and returns 8-11 ns -- which would read as
            # the BEST timing on the plot. Draw sparse points hollow and say so.
            nlow = d['n_active'] < MIN_EV_TIMING if 'n_active' in d.columns \
                else pd.Series(False, index=d.index)
            axs[1].plot(d.loc[~nlow, 'e_drift_Vcm'], d.loc[~nlow, 'time_res_ns'],
                        'o-', color=COL[det], label=f'{det}')
            if nlow.any():
                axs[1].plot(d.loc[nlow, 'e_drift_Vcm'], d.loc[nlow, 'time_res_ns'],
                            'o:', color=COL[det], mfc='none', alpha=.55,
                            label=f'{det} (<{MIN_EV_TIMING} evt — not a resolution)')
        n += 1
    axs[0].set_ylabel('efficiency [%]')
    axs[0].set_title('Efficiency vs drift field\n'
                     'zero field is a null test: it must collapse', fontsize=10)
    axs[1].set_ylabel(r'time resolution $\sigma_t$ [ns]')
    axs[1].set_yscale('log')
    axs[1].set_title('Time resolution vs drift field\n'
                     'leading-pad, $\\sigma=(p_{84.1}-p_{15.9})/2$; hollow = too '
                     'few events to mean anything', fontsize=10)
    for a in axs:
        a.set_xlabel('drift field [V/cm]'); a.grid(alpha=.3); a.legend(fontsize=9)
    fig.tight_layout()
    p = os.path.join(out, 'combined_drift_scan.png')
    fig.savefig(p, dpi=200, bbox_inches='tight'); plt.close(fig)
    print(f'saved {p} ({n} detectors)')


if __name__ == '__main__':
    main()
