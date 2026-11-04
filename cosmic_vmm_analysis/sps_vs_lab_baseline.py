#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Were the high-baseline blocks of VMM 4, 7 and 9 already there at the SPS (July 2026)?

The lab measures every channel's baseline directly with neighbour read-out (sng 1, neighbor_baseline.py).
The SPS VMM runs were taken with sng 0 (no sub-threshold read-outs), so there the comparison uses what both
data sets have: the ADC of each channel's OWN (over-threshold) hits. In the lab the high-baseline channels
have own-hit ADCs far above the rest of their chip (their pulses sit on a raised level). If a block had
existed at the SPS, the same channels would stand out the same way there. VMM 16 ch 33 (lab baseline ~700)
is the positive control: it also stands out in the SPS data.

Inputs:
  --lab   neighbor_baseline.csv of the HV-on sng 1 runs (analysis_hvon/sng/), columns from neighbor_baseline.py
  --sps   one or more SPS hits_store sub-run directories (vmm.npy, ch.npy, adc.npy), e.g.
          .../TB_July26_H4_vmm/runs/run_33/driftscan_gap150V/hits_store
A channel is "high" when its lab baseline is more than --margin ADC above the median of the ordinary
channels of its chip (median taken over channels below 80 ADC, so a big block does not pull it up).

Writes <out>/sps_vs_lab_baseline.json (per chip and channel) and <out>/sps_vs_lab_baseline.png.
Usage:  sps_vs_lab_baseline.py --lab <csv> --sps <hits_store>[,<hits_store>...] --out <dir> [--margin 30]
"""

import argparse
import csv
import glob
import json
import os

import numpy as np


def load_lab(path, min_nb=50):
    lab = {}
    for r in csv.DictReader(open(path)):
        if int(r['n_nb']) < min_nb:
            continue
        own = r['adc_over_median']
        lab[(int(r['vmm']), int(r['ch']))] = (float(r['adc_nb_median']),
                                              float(own) if own not in ('', 'nan') else np.nan)
    return lab


def load_sps(stores):
    v, c, a = [], [], []
    for s in stores:
        for d in sorted(glob.glob(os.path.join(s, '*'))):
            if not os.path.exists(os.path.join(d, 'vmm.npy')):
                continue
            v.append(np.load(os.path.join(d, 'vmm.npy')))
            c.append(np.load(os.path.join(d, 'ch.npy')))
            a.append(np.load(os.path.join(d, 'adc.npy')))
    return np.concatenate(v), np.concatenate(c), np.concatenate(a)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--lab', required=True)
    ap.add_argument('--sps', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--margin', type=float, default=30.0)
    ap.add_argument('--min-sps', type=int, default=20, help='min SPS own hits for an ADC median')
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    lab = load_lab(args.lab)
    v, c, a = load_sps(args.sps.split(','))
    key = v.astype(np.int64) * 64 + c
    order = np.argsort(key, kind='stable')
    key, a = key[order], a[order]
    edges = np.searchsorted(key, np.arange(20 * 64 + 1))

    res = {}
    print(f'{"VMM":>3} {"high ch":>7} | {"lab BL":>13} | {"lab own ADC":>13} | {"SPS own ADC":>13} | {"SPS hits":>13} | ratio high/normal: lab own, SPS own')
    for vm in range(2, 20):
        bl = {ch: lab[(vm, ch)][0] for ch in range(64) if (vm, ch) in lab}
        if not bl:
            continue
        med = float(np.median([b for b in bl.values() if b < 80]))
        rows = []
        for ch in range(64):
            seg = a[edges[vm * 64 + ch]:edges[vm * 64 + ch + 1]]
            n = len(seg)
            lb, lo = lab.get((vm, ch), (np.nan, np.nan))
            rows.append(dict(ch=ch, lab_bl=lb, lab_own=lo, sps_n=int(n),
                             sps_own=float(np.median(seg)) if n >= args.min_sps else np.nan,
                             sps_p10=float(np.percentile(seg, 10)) if n >= args.min_sps else np.nan,
                             high=bool(ch in bl and bl[ch] > med + args.margin)))
        hi = [r for r in rows if r['high']]
        no = [r for r in rows if not r['high'] and r['ch'] in bl]

        def m(rs, k):
            return float(np.nanmedian([r[k] for r in rs])) if rs else np.nan
        summ = {g: {k: m(rs, k) for k in ('lab_bl', 'lab_own', 'sps_own', 'sps_p10', 'sps_n')}
                for g, rs in (('normal', no), ('high', hi))}
        res[vm] = dict(chip_median=med, high=[r['ch'] for r in hi], summary=summ, channels=rows)
        if hi:
            s, h = summ['normal'], summ['high']
            print(f'{vm:>3} {len(hi):>7} | {s["lab_bl"]:5.0f} -> {h["lab_bl"]:5.0f} | {s["lab_own"]:5.0f} -> {h["lab_own"]:5.0f} | '
                  f'{s["sps_own"]:5.0f} -> {h["sps_own"]:5.0f} | {s["sps_n"]:5.0f} -> {h["sps_n"]:5.0f} | '
                  f'{h["lab_own"] / s["lab_own"]:.2f}, {h["sps_own"] / s["sps_own"]:.2f}   ch {[r["ch"] for r in hi]}')
    json.dump(res, open(os.path.join(args.out, 'sps_vs_lab_baseline.json'), 'w'), indent=1)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    show = [vm for vm in (4, 7, 9, 16) if vm in res] + [3, 5]
    fig, axes = plt.subplots(len(show), 1, figsize=(13, 2.6 * len(show)), sharex=True)
    for ax, vm in zip(axes, show):
        R = res[vm]['channels']
        ch = np.array([r['ch'] for r in R])
        hi = np.array([r['high'] for r in R])
        for k, col, lab_ in (('lab_bl', '#1f6fb4', 'lab baseline (neighbour read-out, HV on, Oct)'),
                             ('lab_own', '#e8743b', 'lab own-hit ADC (median)'),
                             ('sps_own', '#222222', 'SPS own-hit ADC (median, July run_33)')):
            y = np.array([r[k] for r in R], float)
            ax.plot(ch, y, 'o', ms=4, color=col, label=lab_, mfc=col if k != 'sps_own' else 'none')
        for x in ch[hi]:
            ax.axvspan(x - 0.5, x + 0.5, color='#e8743b', alpha=0.10, lw=0)
        ax.set_yscale('log')
        ax.set_ylim(25, 1100)
        ax.set_ylabel(f'VMM {vm}\nADC')
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=8, ncol=3, loc='upper left')
    axes[-1].set_xlabel('channel (shaded: lab baseline > chip median + %g ADC)' % args.margin)
    fig.suptitle('High-baseline blocks: lab (Oct 2026) against SPS (July 2026), same chips and channels')
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, 'sps_vs_lab_baseline.png'), dpi=110)
    print('wrote', args.out)


if __name__ == '__main__':
    main()
