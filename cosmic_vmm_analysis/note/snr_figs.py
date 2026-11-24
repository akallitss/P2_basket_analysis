#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Section 8 of "Where do the missed muons go?" in signal-to-noise terms (24 Nov 2026).

S = median muon pulse height (largest pad of the found muons, minus the 40 ADC baseline), in mV with 1.07 ADC/mV
    (vmmsc Offline_ADC, Nov 2025: 3.2 ADC/fC at 3 mV/fC; the VMM ADC reads voltage, so this is gain-independent).
T = threshold above baseline, in mV:
    T_need : median over the chamber of the per-channel thresholds that give the noise rate R (bscan_target
             reports), i.e. what the noise asks for (effect 1 only);
    T_set  : median over the chamber of the threshold actually set (dist_mV of the efficiency pipeline, in
             gain_scan_eff.csv), i.e. after the trim (adds effect 2, the clipped channels).
At a fixed noise rate the threshold is a fixed number of noise sigmas (rate ~ exp(-T^2 / 2 sigma^2)), so S/T is
S/N divided by a constant: its change with gain is the change of S/N.

Usage: snr_figs.py GS_OUT_DIR REPORT_DIR OUT_DIR   (GS_OUT_DIR = $W/analysis/gain_scan/out)
Outputs: c_snr.png, snr_numbers.csv
"""
import json
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

GS, REP, OUT = sys.argv[1:4]
os.makedirs(OUT, exist_ok=True)
ADC_PER_MV, BASE = 1.07, 40.0
DETS = ['P2_IN', 'P2_MID', 'P2_OUT']
COL = {'P2_IN': '#2a78d6', 'P2_MID': '#eb6834', 'P2_OUT': '#1baf7a'}   # the note's chamber colours
GREY, DARK = '#6B6B6B', '#262626'
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Liberation Sans', 'DejaVu Sans'],
                     'font.size': 13, 'axes.edgecolor': GREY, 'axes.labelcolor': DARK, 'xtick.color': GREY,
                     'ytick.color': GREY, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': True,
                     'grid.color': '#E6E6E6', 'legend.frameon': False, 'savefig.dpi': 200, 'savefig.bbox': 'tight',
                     'figure.facecolor': 'white'})

MA = pd.read_csv(os.path.join(GS, 'muon_adc_gain.csv'))
GE = pd.read_csv(os.path.join(GS, 'gain_scan_eff.csv'))
GC = pd.read_csv(os.path.join(GS, 'gain_scan_chip.csv'))
chipdet = GC.groupby('chip').det.first().to_dict()
rows = []
for g in (3.0, 6.0, 9.0):
    for R in (2, 10, 50):
        rep = json.load(open(os.path.join(REP, f'hvon_g{g}_R{R}_report.json')))
        for det in DETS:
            need, ncl = [], 0
            for c, v in rep.items():
                if chipdet.get(int(c)) != det or 'targets' not in v:
                    continue
                med = set(int(x) for x in (v.get('chip_median_used') or []))
                need += [x for k, x in v['targets'].items() if int(k) not in med]
                ncl += len(v.get('clipped') or [])
            m = MA[(MA.gain == g) & (MA.R == R) & (MA.det == det) & MA.label.str.contains('gscan')].iloc[0]
            e = GE[(GE.gain == g) & (GE.R == R) & (GE.det == det) & (GE.family == 'gscan')].iloc[0]
            S = (m.adc_median - BASE) / ADC_PER_MV
            Tn, Ts = float(np.median(need)), float(e.thr_mV)
            rows.append(dict(gain=g, R=R, det=det, S_mV=S, S_fC=S / g, T_need_mV=Tn, T_set_mV=Ts,
                             T_need_fC=Tn / g, ST_need=S / Tn, ST_set=S / Ts, n_chan=len(need), n_clipped=ncl,
                             silent=e.silent))
T = pd.DataFrame(rows)
T.to_csv(os.path.join(OUT, 'snr_numbers.csv'), index=False, float_format='%.3f')
pd.set_option('display.width', 200)
print(T[T.R == 10].round(2).to_string(index=False))

# ---------------------------------------------------------------- S/T vs gain, R = 10 Hz
f, ax = plt.subplots(figsize=(6.8, 4.4))
t = T[T.R == 10]
for d in DETS:
    s = t[t.det == d].sort_values('gain')
    ax.plot(s.gain, s.ST_need, color=COL[d], lw=2, marker='o', ms=8, mec='white', mew=1.5, zorder=3)
    ax.plot(s.gain, s.ST_set, color=COL[d], lw=2, ls='--', marker='o', ms=8, mfc='white', mew=2, zorder=3)
    ax.fill_between(s.gain, s.ST_set, s.ST_need, color=COL[d], alpha=0.10, lw=0)
    ax.text(9.25, s.ST_need.iloc[-1], d, color=DARK, va='center', fontsize=12)
ax.plot([], [], color=GREY, lw=2, marker='o', ms=8, mec='white', label='threshold the noise asks for (effect 1)')
ax.plot([], [], color=GREY, lw=2, ls='--', marker='o', ms=8, mfc='white', mew=2,
        label='threshold actually set, after trim (+ effect 2)')
ax.set_xlim(2.5, 10.3); ax.set_ylim(0, 13)
ax.set_xticks([3, 6, 9])
ax.set_xlabel('VMM gain setting (mV/fC)')
ax.set_ylabel('signal / threshold\n(median muon / threshold for 10 Hz noise)')
ax.legend(fontsize=11, loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=1)
f.savefig(os.path.join(OUT, 'c_snr.png'))
print('written', os.path.join(OUT, 'c_snr.png'))
