#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Is the VMM gain what the configuration says, and where does the threshold sit relative to the muon signal?
Conversion-free check of the electronics gain scan (20 Nov 2026), from the weekend tag dumps.

Per configuration and probe chamber:
  muon ADC   for scintillator-confirmed on-time tags (dt_scint 140-310 ns): the largest ADC of the live probe hits within
             15 mm of the prediction and +-150 ns of tag + latency (the muon's pad), its quartiles and the saturated
             share (ADC >= 1000). Minus the ADC baseline (--baseline, ~40 ADC) it is proportional to the REAL gain.
  edge ADC   the per-channel lower edge of the all-hits ADC spectrum (1st percentile of the 'adchist' histogram,
             16-ADC bins), median over the chamber's channels: the threshold as the ADC sees it.
  ratio      (edge - baseline) / (median muon ADC - baseline): the threshold relative to the muon signal. It does
             not use the nominal mV/fC. If raising the electronics gain helped, this ratio would fall.

Real-gain check: (median muon ADC - baseline) at gain g divided by the same at gain 3, against g / 3. The median is
biased low where many muons saturate (reported).

Usage (from $W, env.sh): muon_adc_gain.py tags pseudo/chunks.json OUT.csv [--baseline 40]
"""
import argparse
import glob
import json
import re

import numpy as np
import pandas as pd

NAME = {3: 'P2_IN', 2: 'P2_MID', 1: 'P2_OUT'}
CHIPS = {3: range(2, 8), 2: range(8, 14), 1: range(14, 20)}


def label(cfg):
    g = re.search(r'gain([0-9.]+)', cfg)[1]
    return f"g{g} " + re.sub(r'^cfg_gain[0-9.]+_peaktime200_opt_cosmic_', '', cfg)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('tags'); ap.add_argument('chunks'); ap.add_argument('out')
    ap.add_argument('--baseline', type=float, default=40.0)
    a = ap.parse_args()
    chunks = json.load(open(a.chunks))
    adc = {}
    hist = {}
    for f in sorted(glob.glob(f'{a.tags}/tags_run*.npz')):
        pid = re.search(r'run(\d+)\.npz', f)[1]
        info = chunks.get(pid)
        if info is None or '_tp' in info['cfg']:
            continue
        lab = label(info['cfg'])
        Z = np.load(f)
        for p in (3, 2, 1):
            k = f'p{p}_'
            if k + 'fired_nr_adc' not in Z.files:
                continue
            n = len(Z[k + 'on'])
            dts = Z[k + 'dt_scint']
            good = Z[k + 'on'] & (dts >= 140) & (dts <= 310)
            fi, fa, fd, fdt = (Z[k + x] for x in ('fired_nr_i', 'fired_nr_adc', 'fired_nr_d', 'fired_nr_dt'))
            m = (fd < 15) & (np.abs(fdt) <= 150)
            best = np.full(n, -1)
            np.maximum.at(best, fi[m], fa[m])
            sel = good & (best >= 0)
            adc.setdefault((lab, p), []).append(best[sel])
        if 'adchist_h' in Z.files:
            H = Z['adchist_h']
            hist[lab] = hist.get(lab, 0) + H
    rows = []
    for (lab, p), lst in adc.items():
        x = np.concatenate(lst)
        if not len(x):
            continue
        H = hist.get(lab)
        edges = []
        if H is not None:
            for v in CHIPS[p]:
                for c in range(64):
                    h = H[v, c]
                    if h.sum() > 200:
                        cdf = np.cumsum(h) / h.sum()
                        edges.append(16 * int(np.searchsorted(cdf, 0.01)) + 8)
        edge = float(np.median(edges)) if edges else np.nan
        q1, med, q3 = np.percentile(x, [25, 50, 75])
        rows.append(dict(label=lab, gain=float(re.search(r'g([0-9.]+)', lab)[1]), det=NAME[p], muons=len(x),
                         adc_q1=q1, adc_median=med, adc_q3=q3, saturated=float((x >= 1000).mean()),
                         edge_adc=edge, ratio=(edge - a.baseline) / (med - a.baseline) if med > a.baseline else np.nan))
    T = pd.DataFrame(rows)
    # real gain vs gain 3, same chamber, gscan configurations of the same R
    T['R'] = T.label.str.extract(r'hvon(\d+)').astype(float)
    ref = T[T.gain == 3.0].set_index(['det', 'R']).adc_median
    T['gain_ratio_vs_g3'] = [(r.adc_median - a.baseline) / (ref.get((r.det, r.R), np.nan) - a.baseline)
                             if 'gscan' in r.label else np.nan for r in T.itertuples()]
    T['nominal_ratio_vs_g3'] = T.gain / 3.0
    T = T.sort_values(['det', 'gain', 'R'])
    T.to_csv(a.out, index=False, float_format='%.3f')
    pd.set_option('display.width', 220)
    print(T.drop(columns=['R']).round(3).to_string(index=False))


if __name__ == '__main__':
    main()
