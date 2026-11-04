#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Baselines of the VMM channels from every neighbour read-out of the campaign, and the chip threshold zero point
from the SPS and from the lab, compared.

Inputs:
  * neighbour read-out groups (neighbor_baseline.py outputs, one directory each): per channel the median ADC of the
    neighbour read-outs = the channel's baseline; a channel counts when it has at least --min-nb read-outs;
  * lab edges (edge_vs_sdt.py): per run and chip the threshold edge (1 % of the ADC spectrum) and the sdt;
  * SPS edges: the per-chip ADC spectra of the July runs (vmm_adc_all.npz) and their sdt (runcfg/sdt.json), as in
    sps_beam_analysis/vmm_threshold_tradeoff/note_figs4.py;
  * the July per-chip "baseline in sdt" of that analysis (BASE: the edge slid down at 0.9 ADC per sdt step to an
    assumed 38 ADC).

The threshold zero point of a chip is the sdt at which its threshold voltage sits at its baseline:
    lab: a straight line edge = e0 + slope * sdt through the chip's lab points with the edge above its baseline;
         zero point = (baseline - e0) / slope;
    SPS: zero point = sdt - (edge - baseline) / slope, for each July run, with the chip's lab slope (0.9 if none).
Each is given with the measured baseline (the chip median of group --ref) and with 38 ADC (which reproduces July).

Usage: baseline_compare.py OUT_STEM --groups A:dirA,B:dirB,... --ref C --lab-edges edges.json --sps-npz npz
       --sps-sdt sdt.json [--min-nb 50]
"""
import argparse
import csv
import json
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from event_display import style, INK, MUTED, COL  # noqa: E402

JULY = {2: 208.1, 3: 182.2, 4: 197.8, 5: 211.7, 7: 200.6, 8: 196.7, 9: 191.1, 10: 210.0, 11: 204.1, 12: 210.6, 13: 188,
        14: 201, 15: 180, 16: 199, 17: 187, 18: 182, 19: 204}
SPS_RUNS = ('run_48', 'run_52', 'run_66')            # the July "opt" runs used for JULY (run_46 left out, as there)
GCOL = {'A': MUTED, 'B': '#7a3e9d', 'C': COL[2], 'D': COL[1]}
CHIPS = range(2, 20)


def load_group(d, min_nb):
    out = {}
    for r in csv.DictReader(open(d + '/neighbor_baseline.csv')):
        # a masked channel is still read out as a neighbour, with ADC 0: not a baseline
        if int(r['n_nb']) >= min_nb and float(r['adc_nb_median'] or 0) >= 20:
            out[(int(r['vmm']), int(r['ch']))] = float(r['adc_nb_median'])
    return out


def edge(h, q=0.01):
    c = np.cumsum(h) / h.sum()
    return int(np.searchsorted(c, q))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('out_stem')
    ap.add_argument('--groups', required=True)
    ap.add_argument('--ref', default='C')
    ap.add_argument('--lab-edges', required=True)
    ap.add_argument('--sps-npz', required=True)
    ap.add_argument('--sps-sdt', required=True)
    ap.add_argument('--min-nb', type=int, default=50)
    ap.add_argument('--block', type=float, default=90, help='channels above this baseline are left out of the chip reference')
    a = ap.parse_args()

    G = {g.split(':')[0]: load_group(g.split(':')[1], a.min_nb) for g in a.groups.split(',')}
    med = {g: {v: float(np.median([b for (vv, c), b in G[g].items() if vv == v])) if any(vv == v for vv, _ in G[g]) else None
               for v in CHIPS} for g in G}
    cov = {g: {v: sum(1 for vv, _ in G[g] if vv == v) for v in CHIPS} for g in G}
    print('channels with a baseline per group:', {g: len(G[g]) for g in G})
    pairs = {}
    for g1, g2 in (('A', 'C'), ('C', 'D'), ('A', 'D')):
        if g1 in G and g2 in G:
            common = [k for k in G[g1] if k in G[g2] and k[0] not in (4, 7, 9)]
            d = np.array([G[g2][k] - G[g1][k] for k in common])
            pairs[f'{g1}->{g2}'] = d
            print(f'{g1} -> {g2}: {len(common)} channels (normal chips): median diff {np.median(d):+.1f} ADC, '
                  f'|diff| <= 3: {np.mean(np.abs(d) <= 3):.2f}, <= 5: {np.mean(np.abs(d) <= 5):.2f}')

    E = json.load(open(a.lab_edges))
    # reference baseline of a chip: median of its ordinary channels (the high-baseline blocks of VMM 4, 7, 9, above
    # --block ADC, do not respond like the rest and would pull the chip median up)
    ref = {}
    for v in CHIPS:
        vals = [b for (vv, _), b in G[a.ref].items() if vv == v and b < a.block]
        ref[v] = float(np.median(vals)) if vals else None
    lab = {}
    for v in CHIPS:
        pts = [(r['sdt'][str(v)], r['edge'][str(v)], int(run)) for run, r in E.items()
               if r['edge'][str(v)] is not None and r['sdt'][str(v)] not in (None, 1023)]
        b = ref.get(v)
        use = [p for p in pts if b is not None and p[1] > b + 5]
        fit = None
        if len({p[0] for p in use}) >= 3:
            x = np.array([p[0] for p in use], float)
            y = np.array([p[1] for p in use], float)
            keep = np.ones(len(x), bool)
            for _ in range(3):                      # robust: drop points far from the line (noisy periods of a chip)
                s, e0 = np.polyfit(x[keep], y[keep], 1)
                res = y - (e0 + s * x)
                mad = 1.4826 * np.median(np.abs(res[keep])) + 1.0
                keep = np.abs(res) < 3 * mad
                if keep.sum() < 3 or len(set(x[keep])) < 3:
                    break
            fit = (float(s), float(e0)) if 0.5 < s < 1.3 else None
        lab[v] = {'pts': pts, 'used': len(use), 'fit': fit}

    z = np.load(a.sps_npz)
    SDT = json.load(open(a.sps_sdt))
    sps = {v: [] for v in CHIPS}
    for run in SPS_RUNS:
        keys = [k for k in z.files if k.startswith(run + '__') and k.endswith('__adc') and 'gain4.5' in k]
        if not keys:
            continue
        h = z[keys[0]].astype(float)
        for v in CHIPS:
            if h[v].sum() > 500:
                sps[v].append((SDT[run][str(v)], edge(h[v]), run))

    rows = []
    print(f'{"VMM":>3} {"bl A":>5} {"bl C":>5} {"bl D":>5} {"slope":>6} {"July(38)":>8} {"SPS(38)":>8} {"SPS(meas)":>9} {"lab(meas)":>9} {"lab(38)":>8}')
    for v in CHIPS:
        b = ref.get(v)
        s = lab[v]['fit'][0] if lab[v]['fit'] else 0.9
        zl = (b - lab[v]['fit'][1]) / s if (lab[v]['fit'] and b is not None) else None
        zl38 = (38 - lab[v]['fit'][1]) / s if lab[v]['fit'] else None
        zs = [sd - (e - b) / s for sd, e, _ in sps[v]] if b is not None else []
        zs38 = [sd - (e - 38) / 0.9 for sd, e, _ in sps[v]]
        row = {'vmm': v, 'baseline': {g: med[g].get(v) for g in G}, 'coverage': {g: cov[g][v] for g in G}, 'slope_lab': s if lab[v]['fit'] else None,
               'n_lab_points': lab[v]['used'], 'zero_july_38': JULY.get(v), 'zero_sps_38': float(np.mean(zs38)) if zs38 else None,
               'zero_sps_meas': float(np.mean(zs)) if zs else None, 'zero_lab_meas': zl, 'zero_lab_38': zl38}
        rows.append(row)

        def f(x, w=8):
            return f'{x:{w}.1f}' if x is not None else f'{"-":>{w}}'
        print(f'{v:>3} {f(med.get("A", {}).get(v), 5)} {f(med.get("C", {}).get(v), 5)} {f(med.get("D", {}).get(v), 5)} '
              f'{f(row["slope_lab"], 6)} {f(row["zero_july_38"])} {f(row["zero_sps_38"])} {f(row["zero_sps_meas"], 9)} '
              f'{f(row["zero_lab_meas"], 9)} {f(row["zero_lab_38"])}')
    json.dump({'chips': rows, 'pairs': {k: {'n': int(len(d)), 'median': float(np.median(d)),
                                            'within3': float(np.mean(np.abs(d) <= 3))} for k, d in pairs.items()}},
              open(a.out_stem + '_baselines.json', 'w'), indent=1)

    # figure 1: baselines per chip and group
    fig, axs = plt.subplots(1, 2, figsize=(15, 4.8), constrained_layout=True, gridspec_kw={'width_ratios': [2.2, 1]})
    ax = axs[0]
    style(ax)
    gs = [g for g in ('A', 'B', 'C', 'D') if g in G]
    for i, g in enumerate(gs):
        for v in CHIPS:
            vals = [b for (vv, _), b in G[g].items() if vv == v]
            if not vals:
                continue
            x = v + (i - (len(gs) - 1) / 2) * 0.2
            q = np.percentile(vals, [10, 50, 90])
            ax.plot([x, x], [q[0], q[2]], color=GCOL[g], lw=2)
            ax.plot(x, q[1], 'o', color=GCOL[g], ms=4, label=g if v == min(vv for vv, _ in G[g]) else None)
    ax.axhline(38, color=INK, ls=':', lw=1)
    ax.text(19.6, 39, '38 (assumed)', fontsize=8, color=INK, ha='right')
    ax.set_ylim(25, 110)
    ax.set_xticks(list(CHIPS))
    ax.set_xlabel('VMM', fontsize=8)
    ax.set_ylabel('channel baseline [ADC]: median and 10-90 % of the channels', fontsize=8)
    ax.set_title('baselines from every neighbour read-out (A: 23 Oct HV off, B: VMM 4/7/9 HV off, C: 23 Oct HV on, '
                 'D: 31 Oct, new connectors)', fontsize=9, color=INK)
    ax.legend(fontsize=8, frameon=False, ncol=4)
    ax = axs[1]
    style(ax)
    for k, d in pairs.items():
        ax.hist(np.clip(d, -15, 15), bins=np.arange(-15.5, 16, 1), histtype='step', lw=2, label=f'{k} ({len(d)} ch)')
    ax.set_xlabel('difference of the same channel between two measurements [ADC]', fontsize=8)
    ax.set_title('is a channel\'s baseline stable?', fontsize=9, color=INK)
    ax.legend(fontsize=8, frameon=False)
    fig.savefig(a.out_stem + '_baselines.png', dpi=120)

    # figure 2: edge against sdt per chip
    fig, axs = plt.subplots(3, 6, figsize=(18, 9), constrained_layout=True, sharey=True)
    for v in CHIPS:
        ax = axs.flat[v - 2]
        style(ax)
        pts = lab[v]['pts']
        if pts:
            x = np.array([p[0] for p in pts])
            y = np.array([p[1] for p in pts])
            new = np.array([p[2] >= 112 for p in pts])
            ax.plot(x[~new], y[~new], 'o', ms=3.5, color=COL[2], label='lab, runs < 112')
            ax.plot(x[new], y[new], 's', ms=3.5, mfc='white', color=COL[1], label='lab, runs >= 112')
        if sps[v]:
            ax.plot([p[0] for p in sps[v]], [p[1] for p in sps[v]], 'D', ms=5, color=INK, label='SPS July')
        b = ref.get(v)
        if b is not None:
            ax.axhline(b, color=COL[3], lw=1.2, ls='--')
        ax.axhline(38, color=MUTED, lw=0.8, ls=':')
        if lab[v]['fit']:
            s, e0 = lab[v]['fit']
            xx = np.array([150, 320])
            ax.plot(xx, e0 + s * xx, color=INK, lw=1)
        for key, mk, col in (('zero_lab_meas', '^', COL[2]), ('zero_sps_meas', 'v', INK), ('zero_july_38', 'x', MUTED)):
            zz = rows[v - 2][key]
            if zz is not None and b is not None:
                ax.plot(zz, b if key != 'zero_july_38' else 38, mk, color=col, ms=8)
        ax.set_xlim(150, 320)
        ax.set_ylim(20, 160)
        ax.set_title(f'VMM {v}' + (f'  slope {lab[v]["fit"][0]:.2f}' if lab[v]['fit'] else ''), fontsize=8.5, color=INK)
        if v == 2:
            ax.legend(fontsize=6.5, frameon=False, loc='upper left')
    fig.supxlabel('threshold register sdt', fontsize=9)
    fig.supylabel('threshold edge [ADC] (dashed: measured baseline, dotted: 38)', fontsize=9)
    fig.suptitle('Threshold edge against sdt: where it meets the baseline is the chip\'s threshold zero point '
                 '(triangle up: lab, down: SPS, x: July with 38)', fontsize=11)
    fig.savefig(a.out_stem + '_edges.png', dpi=110)


if __name__ == '__main__':
    sys.exit(main())
