#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Two tests of where the silent muons come from (20 Nov 2026), from the weekend tag dumps.

For every scintillator-confirmed on-time tag (dt_scint 140-310 ns) whose probe FIRED near the prediction (live pads
within 15 mm and +-150 ns of tag + latency, the 'found' muons):
  max ADC   the largest pad ADC (the per-pad trigger sees this one)
  sum ADC   the sum over those pads, baseline-subtracted per pad (--baseline)
  npads     how many pads

A. Spectrum near threshold (full ADC resolution). The VMM ADC has a period-16 differential non-linearity (some codes
   collect ~2x); the code-width correction is measured from the data themselves: counts per (code mod 16) in the smooth
   region 200-700 ADC, normalised to 1. The corrected max-ADC spectrum, as a function of (ADC - baseline) / (median -
   baseline), shows how dense the muons are just above the threshold edge. Two simple extrapolations of the density
   below the edge down to zero signal (flat = upper estimate, falling linearly to 0 = lower estimate) give the share of
   muons that would be below the edge, to compare with the measured silent share.

B. Charge sharing. max ADC and sum ADC by npads (1, 2, 3+). If sharing between pads is what keeps muons below the
   per-pad threshold, the multi-pad clusters have a much lower max-pad signal than their sum, and the muons close to the
   edge are mostly multi-pad ones.

Groups: REF (gain 4.5 trimmed), gain 3 (all R), gain 6 (R 10 and 50), gain 9 (all R).
Usage (from $W, env.sh): threshold_spectrum.py tags pseudo/chunks.json OUT.json [--baseline 40]
"""
import argparse
import glob
import json
import re

import numpy as np

NAME = {3: 'P2_IN', 2: 'P2_MID', 1: 'P2_OUT'}


def group(cfg):
    if '_tp' in cfg:
        return None
    if 'trimmed_Km18' in cfg:
        return 'REF g4.5'
    m = re.search(r'gain([0-9.]+)_peaktime200_opt_cosmic_gscan_hvon(\d+)', cfg)
    if not m:
        return None
    if m[1] == '6.0' and m[2] == '2':
        return None
    return f'gain {m[1]}'


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('tags'); ap.add_argument('chunks'); ap.add_argument('out')
    ap.add_argument('--baseline', type=float, default=40.0)
    a = ap.parse_args()
    chunks = json.load(open(a.chunks))
    acc = {}
    for f in sorted(glob.glob(f'{a.tags}/tags_run*.npz')):
        pid = re.search(r'run(\d+)\.npz', f)[1]
        info = chunks.get(pid)
        g = group(info['cfg']) if info else None
        if g is None:
            continue
        Z = np.load(f)
        for p in (3, 2, 1):
            k = f'p{p}_'
            if k + 'fired_nr_adc' not in Z.files:
                continue
            n = len(Z[k + 'on'])
            dts = Z[k + 'dt_scint']
            good = Z[k + 'on'] & (dts >= 140) & (dts <= 310)
            fi, fa, fd, fdt, fk = (Z[k + x] for x in ('fired_nr_i', 'fired_nr_adc', 'fired_nr_d', 'fired_nr_dt',
                                                       'fired_nr_key'))
            m = (fd < 15) & (np.abs(fdt) <= 150) & good[fi]
            fi, fa, fk = fi[m], fa[m], fk[m]
            # one entry per (tag, pad): largest ADC of that pad in the window
            o = np.lexsort((-fa, fk, fi))
            fi, fa, fk = fi[o], fa[o], fk[o]
            first = np.r_[True, (np.diff(fi) != 0) | (np.diff(fk) != 0)]
            fi, fa = fi[first], fa[first]
            mx = np.full(n, -1); np.maximum.at(mx, fi, fa)
            sm = np.zeros(n); np.add.at(sm, fi, np.clip(fa - a.baseline, 0, None))
            npd = np.bincount(fi, minlength=n)
            sel = mx >= 0
            r = acc.setdefault((g, NAME[p]), {'max': [], 'sum': [], 'np': [], 'tags': 0})
            r['max'].append(mx[sel]); r['sum'].append(sm[sel]); r['np'].append(npd[sel])
            r['tags'] += int(good.sum())
    out = {}
    for (g, d), r in sorted(acc.items()):
        mx, sm, npd = (np.concatenate(r[x]) for x in ('max', 'sum', 'np'))
        if len(mx) < 200:
            continue
        codes = np.bincount(np.clip(mx, 0, 1023), minlength=1024).astype(float)
        smooth = codes[200:700].reshape(-1, 50)  # not used directly; DNL below
        res = np.array([codes[200:700][(np.arange(200, 700) % 16) == k].sum() for k in range(16)])
        dnl = res / res.mean() if res.sum() else np.ones(16)
        width = dnl[np.arange(1024) % 16]
        dens = codes / np.where(width > 0, width, 1)
        med = float(np.median(mx))
        base = a.baseline
        edge = float(np.percentile(mx, 1))
        # density just above the edge (first 10 % of (median - baseline)) vs the peak density, 8-ADC smoothing
        k8 = np.convolve(dens, np.ones(8) / 8, mode='same')
        span = med - base
        lo_e, hi_e = int(edge), int(edge + 0.1 * span)
        d_edge = float(k8[lo_e:hi_e + 1].mean())
        n_found = len(mx)
        # extrapolations below the edge, in muons: flat to 0 signal (upper), linear to 0 at baseline (lower)
        below_flat = d_edge * max(edge - base, 0)
        below_lin = 0.5 * below_flat
        rel = (np.arange(1024) - base) / span
        hist_rel = np.histogram(rel, bins=np.linspace(0, 3, 61), weights=dens)[0]
        out[f'{g}|{d}'] = {
            'group': g, 'det': d, 'tags': r['tags'], 'found': n_found, 'median_adc': med, 'edge_adc_p1': edge,
            'edge_over_signal': (edge - base) / span,
            'dnl_mod16': dnl.round(3).tolist(),
            'density_rel_edges': np.linspace(0, 3, 61).round(3).tolist(), 'density_rel': hist_rel.round(2).tolist(),
            'below_edge_upper': below_flat / (n_found + below_flat), 'below_edge_lower': below_lin / (n_found + below_lin),
            'npads_share': {str(k): float((np.minimum(npd, 3) == k).mean()) for k in (1, 2, 3)},
            'max_median_by_np': {str(k): float(np.median(mx[np.minimum(npd, 3) == k])) for k in (1, 2, 3)
                                 if (np.minimum(npd, 3) == k).sum() > 20},
            'sum_median_by_np': {str(k): float(np.median(sm[np.minimum(npd, 3) == k])) for k in (1, 2, 3)
                                 if (np.minimum(npd, 3) == k).sum() > 20},
            'near_edge_share_by_np': {str(k): float(((mx - edge) < 0.2 * span)[np.minimum(npd, 3) == k].mean())
                                      for k in (1, 2, 3) if (np.minimum(npd, 3) == k).sum() > 20},
            'max_hist_edges': list(range(0, 1040, 16)),
            'max_hist_by_np': {str(k): np.histogram(mx[np.minimum(npd, 3) == k], bins=range(0, 1040, 16))[0].tolist()
                               for k in (1, 2, 3)},
            'sum_hist_edges': list(range(0, 2080, 32)),
            'sum_hist': np.histogram(sm, bins=range(0, 2080, 32))[0].tolist(),
            'max_hist_all': np.histogram(mx - base, bins=range(0, 2080, 32))[0].tolist(),
        }
        o_ = out[f'{g}|{d}']
        print(f'{g:9s} {d:6s} found {n_found:5d} median {med:5.0f} edge(p1) {edge:4.0f} edge/signal '
              f'{o_["edge_over_signal"]:.2f} | below-edge extrapolation {100 * o_["below_edge_lower"]:.1f}-'
              f'{100 * o_["below_edge_upper"]:.1f} % | npads 1/2/3+ '
              + '/'.join(f'{100 * o_["npads_share"][str(k)]:.0f}' for k in (1, 2, 3)) + ' % | max median by np '
              + str({k: round(v) for k, v in o_['max_median_by_np'].items()}) + ' | near-edge share by np '
              + str({k: round(100 * v) for k, v in o_['near_edge_share_by_np'].items()}))
    print('DNL (code mod 16) example:', next(iter(out.values()))['dnl_mod16'])
    json.dump(out, open(a.out, 'w'))


if __name__ == '__main__':
    main()
