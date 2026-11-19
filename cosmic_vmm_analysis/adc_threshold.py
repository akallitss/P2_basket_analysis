#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Per-channel threshold in ADC and the muon pulse-height spectrum above it (analysis task 2, 19 Nov 2026).

Question (user): at high gain, does the low tail of the muon ADC distribution sit above the threshold while the
efficiency does not change? Then the gain (the charge lost under threshold) is not what limits the efficiency.

Inputs: --dump-diag dumps made with the 19 Nov fields (adchist_h = per-channel ADC histogram of ALL hits, 64 bins of
16 ADC aligned on the ADC DNL period; fired_nr_* = live probe hits near the prediction), pseudo/chunks.json.

  edges   per configuration and channel (summed over its chunks), the threshold in ADC from the lower edge of the
          all-hits spectrum (noise piles up just above threshold), two definitions (user: compute both):
            p01    ADC below which 1 % of the channel's hits lie (linear inside the 16-ADC bin)
            half   50 % point of the rising edge: first bin from below reaching half the spectrum's maximum,
                   interpolated between bin centres
          channels with fewer than --min-hits hits get no edge                        -> adc_edges.csv
  muons   every muon (ON - r OFF tags, --scint window) whose predicted pad's own channel fired within +-win: its
          ADC (largest own-channel hit) minus the channel's edge, both definitions     -> adc_own_minus_edge.csv
  loss    per configuration x chamber (and x chip): distribution of ADC - edge; fit of a Moyal (Landau-like) shape
          to the part above --fit-lo ADC above the edge (truncated likelihood), integral of the fitted shape below 0
          = expected fraction of own-pad signals below threshold                      -> adc_loss.csv, adc_spectra_*.png

Usage: adc_threshold.py TAGS_DIR PSEUDO_DIR OUT_DIR [--exclude-runs 233] [--scint 140:310] [--min-hits 200]
(lxplus, LCG view; reads the dumps only.)
"""
import argparse
import glob
import json
import os
import re
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
NAME = {3: 'P2_IN', 2: 'P2_MID', 1: 'P2_OUT'}
BIN = 16


def label(cfg):
    try:
        from common import label as L
        return L(cfg)
    except ImportError:
        return cfg


def edge_p01(h):
    n = h.sum()
    if n <= 0:
        return np.nan
    c = np.cumsum(h) / n
    i = int(np.searchsorted(c, 0.01))
    prev = c[i - 1] if i > 0 else 0.0
    f = (0.01 - prev) / max(c[i] - prev, 1e-12)
    return BIN * (i + f)


def edge_half(h):
    if h.sum() <= 0:
        return np.nan
    m = h.max()
    i = int(np.argmax(h >= 0.5 * m))
    if i == 0:
        return BIN * 0.5
    x0, x1 = BIN * (i - 0.5), BIN * (i + 0.5)
    y0, y1 = h[i - 1], h[i]
    return x0 + (0.5 * m - y0) / max(y1 - y0, 1e-12) * (x1 - x0)


def gain_of(cfg):
    m = re.search(r'gain(\d+\.\d)', cfg)
    return float(m[1]) if m else np.nan


def fit_loss(x, w, lo):
    """Moyal fit to x >= lo (weights w, truncated likelihood); returns (loc, scale, fraction of the full shape below 0,
    fraction below lo)."""
    from scipy import optimize, stats
    q = x >= lo
    x, w = x[q], w[q]
    if w.sum() < 50:
        return (np.nan,) * 4
    wp = np.clip(w, 0, None)

    def nll(p):
        loc, ls = p
        s = np.exp(ls)
        tail = stats.moyal.sf(lo, loc, s)
        if tail <= 0:
            return 1e30
        return -np.sum(wp * (stats.moyal.logpdf(x, loc, s) - np.log(tail)))
    med = np.average(x, weights=wp)
    r = optimize.minimize(nll, [med * 0.7, np.log(max(np.std(x) * 0.5, 5))], method='Nelder-Mead')
    loc, s = r.x[0], np.exp(r.x[1])
    return loc, s, float(stats.moyal.cdf(0, loc, s)), float(stats.moyal.cdf(lo, loc, s))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('tags'); ap.add_argument('pseudo'); ap.add_argument('out')
    ap.add_argument('--exclude-runs', default='233')
    ap.add_argument('--scint', default='140:310')
    ap.add_argument('--min-hits', type=int, default=200)
    ap.add_argument('--fit-lo', type=float, default=16.0, help='fit the spectrum from this many ADC above the edge')
    ap.add_argument('--pulsed', default='', help='pulsed_channels.json: those channels and their two neighbours get no '
                    'edge in their pulser configurations (the test pulses dominate their spectrum)')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    chunks = {int(k): v for k, v in json.load(open(f'{a.pseudo}/chunks.json')).items()}
    excl = {int(x) for x in a.exclude_runs.split(',') if x}
    lo_s, hi_s = map(float, a.scint.split(':'))

    # ---- 1. edges per configuration and channel
    H = {}
    files = []
    for fz in sorted(glob.glob(f'{a.tags}/tags_run*.npz')):
        pid = int(re.search(r'run(\d+)\.npz', fz)[1])
        if pid not in chunks or chunks[pid]['run'] in excl:
            continue
        Z = np.load(fz)
        if 'adchist_h' not in Z.files:
            continue
        cfg = label(chunks[pid]['cfg'])
        H[cfg] = H.get(cfg, 0) + Z['adchist_h'].astype(np.int64)
        files.append((pid, cfg, fz))
    print(f'{len(files)} chunks with adchist_h, {len(H)} configurations', flush=True)
    skip = set()
    if a.pulsed:
        for k, lst in json.load(open(a.pulsed)).items():
            if not k.startswith('_'):
                skip |= {(label('cfg_' + k), v, c + dc) for v, c in lst for dc in (-1, 0, 1)}
    rows = []
    for cfg, h in H.items():
        for v in range(2, 20):
            for c in range(64):
                hh = h[v, c]
                n = int(hh.sum())
                e1, e2 = (edge_p01(hh), edge_half(hh)) if n >= a.min_hits and (cfg, v, c) not in skip else (np.nan, np.nan)
                rows.append(dict(cfg=cfg, gain=gain_of(cfg), vmm=v, ch=c, hits=n, edge_p01=e1, edge_half=e2))
    E = pd.DataFrame(rows)
    E.to_csv(f'{a.out}/adc_edges.csv', index=False)
    print('\n== edge per configuration (median over channels, ADC)')
    print(E.groupby('cfg')[['edge_p01', 'edge_half', 'hits']].median().round(1).to_string())
    ed = {(r.cfg, r.vmm, r.ch): (r.edge_p01, r.edge_half) for r in E.itertuples()}

    # ---- 2. own-pad ADC of every muon, minus the channel edge
    M = []
    for pid, cfg, fz in files:
        Z = np.load(fz)
        for p in (3, 2, 1):
            g = lambda k: Z[f'p{p}_{k}']
            if f'p{p}_fired_nr_i' not in Z.files:
                continue
            dts = g('dt_scint')
            sel = (dts >= lo_s) & (dts <= hi_s)
            if not sel.any():
                continue
            r = float(g('r_scale'))
            win = float(g('win'))
            key = g('vmm') * 64 + g('ch')
            fi, fdt, fadc, fk = g('fired_nr_i'), g('fired_nr_dt'), g('fired_nr_adc'), g('fired_nr_key')
            own = (fk == key[fi]) & (np.abs(fdt) <= win) & sel[fi]
            if not own.any():
                continue
            df = pd.DataFrame(dict(i=fi[own], adc=fadc[own])).groupby('i').adc.max()
            idx = df.index.to_numpy()
            w = np.where(g('on')[idx], 1.0, -r)
            v, c = g('vmm')[idx], g('ch')[idx]
            e = np.array([ed.get((cfg, int(vv), int(cc)), (np.nan, np.nan)) for vv, cc in zip(v, c)])
            M.append(pd.DataFrame(dict(pid=pid, run=chunks[pid]['run'], cfg=cfg, gain=gain_of(cfg), det=NAME[p],
                                       vmm=v, ch=c, w=w, adc=df.to_numpy(), edge_p01=e[:, 0], edge_half=e[:, 1])))
    M = pd.concat(M, ignore_index=True)
    M['d_p01'] = M.adc - M.edge_p01
    M['d_half'] = M.adc - M.edge_half
    M.to_csv(f'{a.out}/adc_own_minus_edge.csv.gz', index=False)

    # ---- 3. loss estimate
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows = []
    for (cfg, det), s in M.groupby(['cfg', 'det']):
        row = dict(cfg=cfg, gain=gain_of(cfg), det=det, muons=round(float(s.w.sum()), 1),
                   adc_median=float(np.median(s.adc[s.w > 0])))
        for k in ('p01', 'half'):
            q = s[f'd_{k}'].notna()
            x, w = s[f'd_{k}'][q].to_numpy(), s.w[q].to_numpy()
            loc, sc, f0, flo = fit_loss(x, w, a.fit_lo)
            row.update({f'{k}_edge_med': float(np.nanmedian(s[f'edge_{k}'])), f'{k}_d_p05': float(np.nanpercentile(x, 5)),
                        f'{k}_d_median': float(np.nanmedian(x)), f'{k}_fit_loc': loc, f'{k}_fit_scale': sc,
                        f'{k}_loss_below_0': round(f0, 4) if np.isfinite(f0) else np.nan,
                        f'{k}_frac_obs_below_fitlo': round(float(w[x < a.fit_lo].sum() / w.sum()), 4)})
        rows.append(row)
    L = pd.DataFrame(rows).sort_values(['det', 'gain', 'cfg'])
    L.to_csv(f'{a.out}/adc_loss.csv', index=False)
    print('\n== own-pad ADC minus threshold edge; Moyal fit above +%g ADC, its integral below 0 = loss under threshold'
          % a.fit_lo)
    print(L[['cfg', 'det', 'muons', 'adc_median', 'p01_edge_med', 'p01_d_median', 'p01_loss_below_0',
             'half_edge_med', 'half_d_median', 'half_loss_below_0']].to_string(index=False))
    rows = []
    for (cfg, det, v), s in M.groupby(['cfg', 'det', 'vmm']):
        for k in ('p01', 'half'):
            q = s[f'd_{k}'].notna()
            loc, sc, f0, _ = fit_loss(s[f'd_{k}'][q].to_numpy(), s.w[q].to_numpy(), a.fit_lo)
            rows.append(dict(cfg=cfg, det=det, chip=v, edge=k, muons=round(float(s.w.sum()), 1),
                             d_median=float(np.nanmedian(s[f'd_{k}'])), loss_below_0=round(f0, 4) if np.isfinite(f0) else np.nan))
    pd.DataFrame(rows).to_csv(f'{a.out}/adc_loss_per_chip.csv', index=False)

    for det in ('P2_IN', 'P2_MID', 'P2_OUT'):
        s0 = M[M.det == det]
        if s0.empty:
            continue
        cfgs = sorted(s0.cfg.unique(), key=lambda c: (gain_of(c), c))
        fig, ax = plt.subplots(1, 2, figsize=(16, 6))
        bins = np.arange(-160, 1000, 16)
        for k, x_ in zip(('p01', 'half'), ax):
            for cfg in cfgs:
                s = s0[s0.cfg == cfg]
                h = np.histogram(s[f'd_{k}'], bins, weights=s.w)[0]
                x_.step(bins[:-1], h / max(h.sum(), 1), where='post', label=f'{cfg} ({gain_of(cfg):g})')
            x_.axvline(0, color='k', lw=1)
            x_.set_title(f'{det}: own-pad ADC minus threshold edge ({k})')
            x_.set_xlabel('ADC - edge'); x_.set_ylabel('fraction of muons / 16 ADC'); x_.grid(alpha=0.3)
            x_.legend(fontsize=7)
        fig.tight_layout(); fig.savefig(f'{a.out}/adc_spectra_{det}.png', dpi=110); plt.close(fig)
    print(f'\noutputs in {a.out}')


if __name__ == '__main__':
    main()
