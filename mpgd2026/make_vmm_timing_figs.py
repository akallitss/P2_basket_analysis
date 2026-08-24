#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_vmm_timing_figs.py -- the VMM timing arm on two panels.

(a) Coincidence width vs drift voltage, gas A against gas B, all three
    stations, from the trigger-referenced measurement that ran on every
    capture of the campaign (vmm_dream_matching/vmm_trigger_timing.py).
    Gas B is better everywhere and the advantage is largest at low field,
    which is where a faster drift velocity should show.

(b) Where that width comes from, and what would remove it
    (vmm_dream_matching/vmm_timing_budget.py).  Not the clock: the TDC fine
    time is already applied and toggling it changes nothing, and BCID
    quantisation alone is only 22.5/sqrt(12) = 6.5 ns.  The trigger channel --
    a scintillator read through a VMM discriminator -- is the largest single
    removable term.

Usage:  python3 make_vmm_timing_figs.py [-o OUTDIR]
"""

import argparse
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BY_HV = ('/local/home/ak271430/Documents/PostDocSaclay/P2_basket_analysis/'
         'sps_beam_analysis/vmm_dream_matching/vmm_trigger_timing_by_hv.csv')
GAS_A, GAS_B = 'Ar/CO2/Iso 93/5/2', 'Ar/CF4/Iso 88/10/2'
C = {'P2_IN': '#1f77b4', 'P2_MID': '#ff7f0e', 'P2_OUT': '#2ca02c'}

# vmm_timing_budget.py at the NOMINAL working point (run_36/operating_00,
# mesh 450 / drift 750, gas A), solved with the Gaussian core-fit estimator --
# the same one vmm_efficiency.py fits, so these are directly comparable with
# every sigma quoted elsewhere. The rms-inside-+-120 ns estimator gives the
# same ordering but sits ~30 % higher; the two must never be mixed.
BUDGET = {'measured P2_MID vs trigger': 24.2,
          'trigger channel': 18.4,
          'P2_MID intrinsic': 15.6,
          'per-channel t0 spread': 10.8,
          'BCID quantisation': 22.5 / np.sqrt(12)}
# same solve on gas B (run_66): trigger 16.4, P2_MID 12.2, P2_OUT 24.7 ns.
PROJECTED = np.sqrt(15.6 ** 2 + 10.4 ** 2)   # DREAM timestamp in place of it

plt.rcParams.update({
    'font.size': 14, 'axes.titlesize': 14, 'axes.labelsize': 17,
    'xtick.labelsize': 16, 'ytick.labelsize': 16,
    'legend.fontsize': 9.5, 'figure.facecolor': 'white',
    'axes.grid': True, 'grid.alpha': 0.3, 'axes.axisbelow': True,
})
NL = chr(10)


def panel_gas(ax):
    d = pd.read_csv(BY_HV)
    d = d[(d.mesh == 450) & (d.n_captures >= 5)]
    for st in ('P2_MID', 'P2_OUT'):
        for gas, ls, mk, lbl in ((GAS_A, '--', 'o', 'gas A  Ar/CO2/iC4H10 93/5/2'),
                                 (GAS_B, '-', 's', 'gas B  Ar/CF4/iC4H10 88/10/2')):
            s = d[(d.station == st) & (d.gas == gas)].sort_values('drift')
            if len(s) < 2:
                continue
            ax.plot(s['drift'], s['sigma_ns'], ls=ls, marker=mk, ms=8,
                    lw=2.0, color=C[st],
                    label=f'{st}, {lbl.split()[1]}' if st == 'P2_MID'
                    else None)
    for st, dy in (('P2_MID', -13), ('P2_OUT', 11)):
        s = d[(d.station == st) & (d.gas == GAS_A)].sort_values('drift')
        if len(s):
            ax.annotate(st, (s['drift'].iloc[-1], s['sigma_ns'].iloc[-1]),
                        textcoords='offset points', xytext=(9, dy),
                        color=C[st], fontsize=10.5, fontweight='bold',
                        va='center')
    ax.axhline(20, color='#d62728', ls=':', lw=1.8)
    ax.text(602, 21.5, 'the 20 ns P2 goal', color='#d62728', fontsize=9.5)
    ax.set_xlim(590, 890)
    ax.set_ylim(0, 62)
    ax.set_xlabel('drift voltage [V]   (mesh fixed at 450 V)')
    ax.set_ylabel('coincidence width sigma [ns]')
    ax.set_title('(a) trigger-referenced width vs drift field' + NL
                 + 'dashed = gas A, solid = gas B (CF4): better everywhere, '
                   'most at low field' + NL
                 + '(P2_IN omitted -- raised thresholds, few captures, fits '
                   'unstable)', fontsize=10.8)
    h = [plt.Line2D([], [], color='#555', ls='--', marker='o',
                    label='gas A  Ar/CO2/iC4H10 93/5/2'),
         plt.Line2D([], [], color='#555', ls='-', marker='s',
                    label='gas B  Ar/CF4/iC4H10 88/10/2')]
    ax.legend(handles=h, loc='upper right', framealpha=0.96)


def panel_budget(ax):
    order = ['measured\nP2_MID vs trigger', 'trigger\nchannel',
             'P2_MID\nintrinsic', 'per-channel\nt0 spread',
             'BCID\nquantisation']
    cols = ['#333333', '#d62728', '#ff7f0e', '#9467bd', '#7f9f7f']
    y = np.arange(len(order))[::-1]
    vals = [BUDGET[k.replace(chr(10), ' ')] for k in order]
    ax.barh(y, vals, color=cols, height=0.62)
    for yi, v in zip(y, vals):
        ax.text(v + 0.5, yi, f'{v:.1f} ns', va='center', fontsize=10.5,
                fontweight='bold')
    ax.set_yticks(y)
    ax.set_yticklabels(order)
    ax.axvline(PROJECTED, color='#1f77b4', ls='--', lw=2.0)
    ax.text(PROJECTED + 0.5, 2.5,
            f'with a DREAM-referenced timestamp{NL}'
            f'instead: {PROJECTED:.1f} ns',
            color='#1f77b4', fontsize=9.2, va='center')
    ax.set_xlim(0, 30)
    ax.set_xlabel('rms contribution [ns]')
    ax.set_title('(b) what the width is made of' + NL
                 + 'the clock is not the limit; the trigger channel is larger '
                   'than the chamber', fontsize=11.3)
    ax.grid(axis='y', alpha=0)
    ax.text(0.97, 0.06, 'Gaussian core fit, run_36/operating_00,'
            + NL + 'mesh 450 / drift 750, gas A',
            transform=ax.transAxes, ha='right', va='bottom',
            fontsize=8.2, color='#666', style='italic',
            bbox=dict(fc='white', ec='#ccc', boxstyle='round,pad=0.3'))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-o', '--out', default='figs')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14.6, 5.6),
                                   layout='constrained')
    panel_gas(ax1)
    panel_budget(ax2)
    fig.suptitle('VMM3a timing: the gas works, and the clock was never the '
                 'limit', fontsize=13.5, y=1.02)
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(a.out, f'vmm_timing.{ext}'), dpi=170,
                    bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print('  [vmm_timing]')


if __name__ == '__main__':
    main()
