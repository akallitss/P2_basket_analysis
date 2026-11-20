#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Efficiency vs the VMM electronics settings of the weekend grid (20 Nov 2026), from the muon funnel.

Per configuration (electronics gain 3 / 4.5 / 6 / 9 mV/fC x noise-rate target R = 2 / 10 / 50 Hz, the gain-4.5
threshold scan "bscan", and REF) and per probe chamber: the funnel classes of the muons expected on a live pad
(found / silent / far / late, with and without the paddle cut), next to the threshold of that configuration
expressed in charge:

  thr_fC = (threshold - baseline) [mV] / gain [mV/fC]

with (threshold - baseline) = dist_mV of the efficiency pipeline (results/current/by_config.csv, by_chip_config.csv:
the threshold DAC setting of each channel against its measured baseline, median over the chamber / chip). It is a
nominal charge: the gain is the VMM nominal setting, not a measured transfer function.

Per chip as well (the chip of the predicted pad): found / silent vs that chip's thr_fC, which spreads the lever arm.

Usage (from $W, env.sh): gain_scan_eff.py analysis/muon_funnel/summary.json results/current OUT_DIR
Outputs: gain_scan_eff.csv (per configuration and chamber), gain_scan_chip.csv (per configuration and chip),
         gain_scan_eff.json (for the note), printed tables.
"""
import argparse
import json
import os
import re

import numpy as np
import pandas as pd

K = ['found', 'probe_silent', 'probe_far', 'probe_near_late']
INSIDE = ['inside >60', 'inside 30-60', 'inside 15-30', 'inside 0-15']


def key(lab):
    """funnel label 'g4.5 bscan_hvon2_v1_hvm' -> (pipeline config name, gain, R, family)."""
    g, name = lab.split(' ', 1)
    gain = float(g[1:])
    name = re.sub(r'(_v1)?_hvm$', '', name)
    cfg = f'gain{g[1:]}_{name}' if name.startswith('gscan') else name
    m = re.search(r'hvon(\d+)', name)
    fam = 'REF' if name.startswith('trimmed') else ('bscan' if name.startswith('bscan') else 'gscan')
    return cfg, gain, int(m[1]) if m else None, fam


def frac(c):
    n = sum(c.get(x, 0) for x in K)
    return n, {x: c.get(x, 0) / n if n else np.nan for x in K}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('summary'); ap.add_argument('resdir'); ap.add_argument('out')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    G = json.load(open(a.summary))['groups']
    B = pd.read_csv(os.path.join(a.resdir, 'by_config.csv'))
    B = B[B.config != 'config'].astype({'dist_mV': float})
    BC = pd.read_csv(os.path.join(a.resdir, 'by_chip_config.csv')).astype({'dist_mV': float})
    rows, crows = [], []
    for lab, g in G.items():
        if lab in ('REF', 'ALL'):
            continue
        cfg, gain, R, fam = key(lab)
        for p, P in g['probes'].items():
            n, f = frac(P['on'])
            ins = {}
            for b in INSIDE:
                for x, v in P.get('paddle', {}).get(b, {}).items():
                    ins[x] = ins.get(x, 0) + v
            ni, fi = frac(ins)
            b = B[(B.config == cfg) & (B.det == p)]
            d = float(b.dist_mV.iloc[0]) if len(b) else np.nan
            rows.append(dict(label=lab, config=cfg, family=fam, gain=gain, R=R, det=p, runs=len(set(g['runs'])),
                             thr_mV=d, thr_fC=d / gain, n=n, found=f['found'], silent=f['probe_silent'],
                             far=f['probe_far'], late=f['probe_near_late'], n_in=ni, found_in=fi['found'],
                             silent_in=fi['probe_silent']))
            for chip, c in P['chip'].items():
                nc, fc = frac(c)
                bc = BC[(BC.config == cfg) & (BC.chip == int(chip))]
                dc = float(bc.dist_mV.iloc[0]) if len(bc) else np.nan
                crows.append(dict(label=lab, config=cfg, family=fam, gain=gain, R=R, det=p, chip=int(chip),
                                  thr_mV=dc, thr_fC=dc / gain, n=nc, found=fc['found'], silent=fc['probe_silent']))
    T = pd.DataFrame(rows).sort_values(['det', 'family', 'gain', 'R'])
    C = pd.DataFrame(crows).sort_values(['det', 'chip', 'gain', 'R'])
    T.to_csv(os.path.join(a.out, 'gain_scan_eff.csv'), index=False, float_format='%.4f')
    C.to_csv(os.path.join(a.out, 'gain_scan_chip.csv'), index=False, float_format='%.4f')
    # simple straight-line fits of silent and found vs thr_fC (config level, binomial weights), per chamber
    fits = {}

    def line(t, y, tag):
        t = t[np.isfinite(t.thr_fC) & (t.n > 100)]
        if len(t) < 3:
            return
        w = t.n / np.clip(t[y] * (1 - t[y]), 1e-3, None)
        p_, cov = np.polyfit(t.thr_fC, t[y], 1, w=np.sqrt(w), cov='unscaled')
        chi2 = float((w * (t[y] - np.polyval(p_, t.thr_fC)) ** 2).sum())
        fits[tag] = {'slope_per_fC': float(p_[0]), 'slope_err': float(np.sqrt(cov[0, 0])), 'at_0fC': float(p_[1]),
                     'at_0fC_err': float(np.sqrt(cov[1, 1])), 'chi2_ndf': chi2 / max(len(t) - 2, 1),
                     'range_fC': [float(t.thr_fC.min()), float(t.thr_fC.max())], 'points': int(len(t))}

    # silent is the position-free class (valid also where the position mapping is suspect); 'found' only for gscan/REF
    for d, t in T.groupby('det'):
        line(t, 'silent', f'{d}:silent:all')
        line(t[t.thr_fC < 5], 'silent', f'{d}:silent:below5fC')
        line(t[t.family != 'bscan'], 'found', f'{d}:found:no_bscan')
    for d, c in C.groupby('det'):
        line(c, 'silent', f'{d}:silent:chips_all')
        line(c[c.thr_fC < 5], 'silent', f'{d}:silent:chips_below5fC')
    json.dump({'configs': T.replace({np.nan: None}).to_dict('records'),
               'chips': C.replace({np.nan: None}).to_dict('records'), 'fits': fits},
              open(os.path.join(a.out, 'gain_scan_eff.json'), 'w'), indent=1)
    pd.set_option('display.width', 220)
    cols = ['label', 'runs', 'gain', 'R', 'thr_mV', 'thr_fC', 'n', 'found', 'silent', 'far', 'found_in', 'silent_in']
    for d, t in T.groupby('det'):
        print(f'== {d}')
        print(t[cols].sort_values('thr_fC').round(3).to_string(index=False))
    print('== straight lines vs thr_fC (config level)')
    for k, v in fits.items():
        print(f'  {k}: {v["slope_per_fC"] * 100:+.1f} +- {v["slope_err"] * 100:.1f} points per fC; at 0 fC '
              f'{v["at_0fC"]:.3f} +- {v["at_0fC_err"]:.3f}; chi2/ndf {v["chi2_ndf"]:.1f}; {v["points"]} points '
              f'{v["range_fC"][0]:.1f}-{v["range_fC"][1]:.1f} fC')


if __name__ == '__main__':
    main()
