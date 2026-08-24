#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_stability_fig.py -- the error budget of a single-run VMM efficiency.

Backup-slide material, and the reason no VMM number in the talk carries a
Clopper-Pearson error bar.  Every point below is the SAME chip configuration
file (gain 3.0 mV/fC, peaking 200 ns), the SAME HV (mesh 450 / drift 750) and
the SAME gas (Ar/CO2/iC4H10 93/5/2), read off each run's own run_config.json.

  within one run, minutes apart      0.2 - 0.3 points
  across 24 h                        P2_OUT 5.6, P2_MID 10.3, P2_IN 3.6 points

So the uncertainty on any single-run VMM efficiency is a several-point
run-to-run systematic, not the +-0.04 % binomial interval.  What it is NOT:
the per-chip ADC percentiles and the live-channel counts are flat across the
whole window, so it is neither a gain droop nor a threshold change.

Usage:  python3 make_stability_fig.py [-o OUTDIR]
"""

import argparse
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

TABLE = ('/local/home/ak271430/Documents/PostDocSaclay/P2_basket_analysis/'
         'sps_beam_analysis/vmm_dream_matching/efficiency_conditions.csv')
CHIP = 'p2b-config-cern-ext_gain3.0_peaktime200.txt'
C = {'P2_IN': '#1f77b4', 'P2_MID': '#ff7f0e', 'P2_OUT': '#2ca02c'}

plt.rcParams.update({
    'font.size': 14, 'axes.titlesize': 14, 'axes.labelsize': 17,
    'xtick.labelsize': 16, 'ytick.labelsize': 16,
    'legend.fontsize': 10, 'figure.facecolor': 'white',
    'axes.grid': True, 'grid.alpha': 0.3, 'axes.axisbelow': True,
})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-o', '--out', default='figs')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    t = pd.read_csv(TABLE)
    t['start'] = pd.to_datetime(t['start'])
    r = t[(t.chip == CHIP) & (t.mesh == 450) & (t.drift == 750)
          & (t.gas == 'Ar/CO2/Iso 93/5/2')].sort_values('start')

    fig, ax = plt.subplots(figsize=(11.0, 5.8))
    for st, col in (('P2_OUT', 'eff_OUT'), ('P2_MID', 'eff_MID'),
                    ('P2_IN', 'eff_IN')):
        ax.plot(r['start'], r[col], marker='o', ms=9, lw=2.0, color=C[st],
                label=f'{st}   ({100 * (r[col].max() - r[col].min()):.1f} '
                      f'points end to end)')
        first, last = r.iloc[0], r.iloc[-1]
        for row, dx, ha in ((first, -8, 'right'), (last, 8, 'left')):
            ax.annotate(f"{row[col]:.3f}", (row['start'], row[col]),
                        textcoords='offset points', xytext=(dx, 0),
                        ha=ha, va='center', fontsize=10, color=C[st],
                        fontweight='bold')

    for _, row in r.drop_duplicates('run').iterrows():
        ax.annotate(row['run'], (row['start'], 0.795), rotation=90,
                    fontsize=8.5, color='#666', ha='center', va='top')

    ax.set_ylim(0.08, 0.82)
    ax.set_ylabel('absolute efficiency (uRWELL-referenced)')
    ax.set_xlabel('time')
    ax.set_title('Identical chip config, identical HV, identical gas, 24 h '
                 'apart\n'
                 'within a run it repeats to 0.2-0.3 points; across the day '
                 'it does not', fontsize=11.5)
    ax.legend(loc='center left', framealpha=0.96)
    fig.autofmt_xdate()
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(a.out, f'vmm_stability.{ext}'), dpi=170,
                    bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print('  [vmm_stability]')


if __name__ == '__main__':
    main()
