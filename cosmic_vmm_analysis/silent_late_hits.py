#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Are 'silent' muons at high electronics gain really silent, or late? (20 Nov 2026)

For scintillator-confirmed on-time tags (dt_scint 140-310 ns) with no live probe hit within +-win of tag + latency
(n_any == 0, 'silent'): is there a live probe hit within 15 mm of the prediction later (or earlier), in bins of time
from tag + latency (the dumps keep every hit within 15 mm and +-2 us: fired_nr_*)? Random rate: the same for the
non-silent tags is not a control, so the rate is compared with the -2 us .. -1 us side band of the same tags.
Also the time distribution of the nearest hit for the FOUND tags (within +-2 us), per group: does the muon hit
time move with the gain?

Usage (from $W, env.sh): silent_late_hits.py tags pseudo/chunks.json
"""
import glob
import json
import re
import sys

import numpy as np

NAME = {3: 'P2_IN', 2: 'P2_MID', 1: 'P2_OUT'}
BINS = [-2000, -1000, -100, 100, 200, 300, 500, 1000, 2000]


def group(cfg):
    if '_tp' in cfg:
        return None
    if 'trimmed_Km18' in cfg:
        return 'REF g4.5'
    m = re.search(r'gain([0-9.]+)_peaktime200_opt_cosmic_gscan_hvon(\d+)', cfg)
    return f'gain {m[1]}' if m else None


def main():
    tags, chunks = sys.argv[1], json.load(open(sys.argv[2]))
    acc = {}
    for f in sorted(glob.glob(f'{tags}/tags_run*.npz')):
        pid = re.search(r'run(\d+)\.npz', f)[1]
        info = chunks.get(pid)
        g = group(info['cfg']) if info else None
        if g is None:
            continue
        Z = np.load(f)
        for p in (3, 2, 1):
            k = f'p{p}_'
            if k + 'n_any' not in Z.files:
                continue
            n = len(Z[k + 'on'])
            dts = Z[k + 'dt_scint']
            good = Z[k + 'on'] & (dts >= 140) & (dts <= 310)
            sil = good & (Z[k + 'n_any'] == 0)
            fnd = good & (Z[k + 'n_any'] > 0)
            fi, fdt = Z[k + 'fired_nr_i'], Z[k + 'fired_nr_dt']
            a = acc.setdefault((g, NAME[p]), {'n_sil': 0, 'sil_bins': np.zeros(len(BINS) - 1),
                                               'found_dt': []})
            a['n_sil'] += int(sil.sum())
            for lo, hi, j in zip(BINS[:-1], BINS[1:], range(len(BINS) - 1)):
                m = (fdt >= lo) & (fdt < hi)
                has = np.zeros(n, bool)
                has[fi[m]] = True
                a['sil_bins'][j] += (has & sil).sum()
            # nearest-in-time hit of the found tags
            o = np.argsort(np.abs(fdt))
            first = np.full(n, np.nan)
            seen = np.zeros(n, bool)
            for i_, t_ in zip(fi[o], fdt[o]):
                if not seen[i_]:
                    seen[i_] = True
                    first[i_] = t_
            a['found_dt'].append(first[fnd & np.isfinite(first)])
    for (g, d), a in sorted(acc.items()):
        ns = a['n_sil']
        side = a['sil_bins'][0] / 1000.0          # per ns, -2 .. -1 us
        cells = []
        for lo, hi, c in zip(BINS[:-1], BINS[1:], a['sil_bins']):
            exp = side * (hi - lo)
            cells.append(f'{lo}..{hi}: {100 * c / ns:.1f} % (rand {100 * exp / ns:.1f})')
        fd = np.concatenate(a['found_dt'])
        q = np.percentile(fd, [5, 50, 95]) if len(fd) else [np.nan] * 3
        print(f'{g:9s} {d:6s} silent {ns:5d}: near hit (<=15 mm) at ns from tag+lat: ' + ' | '.join(cells))
        print(f'{"":16s} found: hit time from tag+lat 5/50/95 % = {q[0]:.0f} / {q[1]:.0f} / {q[2]:.0f} ns')


if __name__ == '__main__':
    main()
