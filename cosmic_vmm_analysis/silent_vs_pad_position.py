#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Are the silent muons the ones that cross near a pad edge? (20 Nov 2026)

The VMM triggers per pad. A muon crossing near a pad boundary shares its charge between 2-4 pads, each of which can
stay below threshold although the total is large. If that is what makes the chamber silent, the silent share must
grow with the distance of the predicted point from the centre of its pad.

Per configuration group and probe chamber, scintillator-confirmed on-time tags (dt_scint 140-310 ns) from the weekend
dumps: silent = no live probe hit at all within +-win of tag + latency (n_any == 0); position in the pad = the
prediction (px, py) minus the predicted pad's centre (cx, cy), as |dx|, |dy| and r. Bins of r (mm) and of the
'edge distance' min(6 - |dx|, 5.8 - |dy|) (how close to the nearest pad edge, pads ~12 x 11.6 mm; the pad axes of the
fan are rotated, so this is approximate).

The prediction has a resolution of a few mm, which smears the pattern; the trend, not the absolute, is the test.

Usage (from $W, env.sh): silent_vs_pad_position.py tags pseudo/chunks.json
"""
import glob
import json
import re
import sys

import numpy as np

NAME = {3: 'P2_IN', 2: 'P2_MID', 1: 'P2_OUT'}


def group(cfg):
    if '_tp' in cfg:
        return None
    if 'trimmed_Km18' in cfg:
        return 'REF (gain 4.5)'
    m = re.search(r'gain([0-9.]+)_peaktime200_opt_cosmic_gscan_hvon(\d+)', cfg)
    if m and m[1] in ('3.0', '9.0'):
        return f'gain {m[1]} (all R)'
    return None


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
            dts = Z[k + 'dt_scint']
            sel = Z[k + 'on'] & (dts >= 140) & (dts <= 310)
            dx = np.abs(Z[k + 'px'] - Z[k + 'cx'])[sel]
            dy = np.abs(Z[k + 'py'] - Z[k + 'cy'])[sel]
            sil = (Z[k + 'n_any'] == 0)[sel]
            a = acc.setdefault((g, NAME[p]), {'dx': [], 'dy': [], 's': []})
            a['dx'].append(dx); a['dy'].append(dy); a['s'].append(sil)
    rb = [0, 1.5, 3, 4.5, 6, 7.5, 9]
    eb = [-9, 0, 1, 2, 3, 4, 7]
    for (g, d), a in sorted(acc.items()):
        dx, dy, s = (np.concatenate(a[x]) for x in ('dx', 'dy', 's'))
        r = np.hypot(dx, dy)
        e = np.minimum(6 - dx, 5.8 - dy)
        print(f'== {g}, {d}: {len(s)} muons, silent {100 * s.mean():.1f} %')
        print('   r from pad centre (mm): ' + ' | '.join(
            f'{lo}-{hi}: {100 * s[(r >= lo) & (r < hi)].mean():.1f} % (n {((r >= lo) & (r < hi)).sum()})'
            for lo, hi in zip(rb[:-1], rb[1:]) if ((r >= lo) & (r < hi)).sum() > 30))
        print('   distance to the nearest pad edge (mm; < 0 = predicted outside the pad): ' + ' | '.join(
            f'{lo}..{hi}: {100 * s[(e >= lo) & (e < hi)].mean():.1f} % (n {((e >= lo) & (e < hi)).sum()})'
            for lo, hi in zip(eb[:-1], eb[1:]) if ((e >= lo) & (e < hi)).sum() > 30))


if __name__ == '__main__':
    main()
