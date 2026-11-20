#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Are the extra silent muons at high electronics gain the ones on CLIPPED channels? (20 Nov 2026)

bscan_target.py could not trim some channels down to their noise-rate target (trim > 31: the channel stays higher).
For every scintillator-confirmed on-time tag (dt_scint 140-310 ns) of a gscan configuration, the predicted pad's
channel (vmm, ch) is looked up in the report of that gain and R: clipped, without a measured crossing (chip median
used), or normal. Silent = no live probe hit at all within +-win of tag + latency (n_any == 0).

Usage (from $W, env.sh): silent_vs_clipped.py tags pseudo/chunks.json REPORT_DIR
"""
import glob
import json
import os
import re
import sys

import numpy as np


def main():
    tags, chunks, rdir = sys.argv[1], json.load(open(sys.argv[2])), sys.argv[3]
    reports = {}
    acc = {}
    for f in sorted(glob.glob(f'{tags}/tags_run*.npz')):
        pid = re.search(r'run(\d+)\.npz', f)[1]
        info = chunks.get(pid)
        if not info or '_tp' in info['cfg']:
            continue
        m = re.search(r'gain([0-9.]+)_peaktime200_opt_cosmic_gscan_hvon(\d+)_hvm$', info['cfg'])
        if not m:
            continue
        key = (m[1], int(m[2]))
        if key not in reports:
            p = os.path.join(rdir, f'hvon_g{m[1]}_R{m[2]}_report.json')
            if not os.path.exists(p):
                continue
            R = json.load(open(p))
            reports[key] = {
                'clipped': {(int(c), int(ch)) for c, v in R.items() for ch in (v['clipped'] or [])},
                'median': {(int(c), int(ch)) for c, v in R.items() for ch in (v.get('chip_median_used') or [])}}
        rep = reports[key]
        Z = np.load(f)
        for p in (3, 2, 1):
            k = f'p{p}_'
            if k + 'n_any' not in Z.files:
                continue
            dts = Z[k + 'dt_scint']
            sel = Z[k + 'on'] & (dts >= 140) & (dts <= 310)
            vm, ch, sil = Z[k + 'vmm'][sel], Z[k + 'ch'][sel], (Z[k + 'n_any'] == 0)[sel]
            for v, c, s in zip(vm, ch, sil):
                cls = 'clipped' if (int(v), int(c)) in rep['clipped'] else (
                    'no crossing (chip median)' if (int(v), int(c)) in rep['median'] else 'normal')
                a = acc.setdefault((f'gain {key[0]} R={key[1]}', cls), [0, 0])
                a[0] += 1
                a[1] += int(s)
    for (g, cls), (n, s) in sorted(acc.items()):
        e = np.sqrt(s / n * (1 - s / n) / n) if n else 0
        print(f'{g:16s} {cls:26s} muons {n:5d}  silent {100 * s / n:5.1f} +- {100 * e:.1f} %')


if __name__ == '__main__':
    main()
