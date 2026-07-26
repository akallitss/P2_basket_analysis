#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
write_analysis_readme.py -- describe an urw_referenced_efficiency output tree.

Generated, not hand-written, so the numbers in the README cannot drift away
from the CSV they describe.  Run after urw_p2_efficiency.py + plot_hv_curves.py.
"""
import os
import glob
import argparse
from datetime import datetime

import numpy as np
import pandas as pd


def repeat_check(df):
    """Same HV, different time: what the efficiency did in between.

    drift_700 (01:57) and nominal_00 (02:50) are the identical HV point, so the
    difference between them is a direct handle on drift with time/dose that the
    HV curves would otherwise absorb.
    """
    have = set(df['sub_run'])
    if not {'drift_700', 'nominal_00'} <= have:
        return None
    out = []
    for st, g in df.groupby('station'):
        a = g[g['sub_run'] == 'drift_700']
        b = g[g['sub_run'] == 'nominal_00']
        if len(a) and len(b):
            out.append((st, float(a['eff'].iloc[0]), float(b['eff'].iloc[0])))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', required=True)
    ap.add_argument('--run', required=True)
    args = ap.parse_args()

    csv = os.path.join(args.dir, f'urw_p2_efficiency_{args.run}.csv')
    df = pd.read_csv(csv)
    pngs = sorted(os.path.basename(p) for p in glob.glob(os.path.join(args.dir, '*.png')))
    curves = [p for p in pngs if p.startswith('efficiency_vs_')]
    summary = [p for p in pngs if p.startswith('summary_')]
    # the per-station PNGs live in per_subrun/ once the tree is tidied, but are
    # still at the top level straight out of the stage
    per_sub = ([os.path.basename(p) for p in
                glob.glob(os.path.join(args.dir, 'per_subrun', '*.png'))]
               or [p for p in pngs if p not in curves + summary])

    L = []
    w = L.append
    w(f'# P2 efficiency referenced to uRWELL tracks — `{args.run}`\n')
    w(f'Generated {datetime.now().isoformat(timespec="seconds")} by '
      f'`~/P2_basket_analysis/sps_beam_analysis/urw_reference/urw_p2_efficiency.py` + `plot_hv_curves.py`.\n')
    w('The three P2 BASKET stations are measured against tracks from the two EIC')
    w('uRWELL references that bracket them (z = 0 and 1370 mm). The uRWELL points')
    w('at each P2 plane to about 1 mm, three times better than a 12 mm P2 pad, so')
    w('this is a tracking-referenced efficiency rather than a tag-and-probe one.')
    w('Method, systematics and traps: `~/P2_basket_analysis/sps_beam_analysis/urw_reference/URW_TRACKING_HANDOFF_2026-07-25.md`')
    w('§13. The mapping and alignment it rests on are frozen in')
    w('`../MAPPING_AND_ALIGNMENT.md`.\n')

    w('## Efficiency curves\n')
    for p in curves:
        w(f'* `{p}`')
    if not curves:
        w('* (none — no scan family had enough points)')
    w('')
    w('**How the points are grouped.** A mesh scan here does not hold the drift')
    w('fixed: the drift tracks the mesh so the drift *field* stays constant, so the')
    w('invariant that defines a scan family is **drift − mesh**, not drift. The')
    w('`drift_*` sub-runs instead hold the mesh and walk the drift, and get their')
    w('own figure. Repeated HV settings are combined by adding counts.\n')

    w('## Files\n')
    w('| file | what |')
    w('|---|---|')
    w(f'| `urw_p2_efficiency_{args.run}.csv` | one row per (sub_run, station): '
      'efficiency + interval, residual widths, frame fit, HV |')
    w(f'| `urw_p2_efficiency_{args.run}.json` | the same plus per-station cut '
      'flow, efficiency vs probe radius, full affine matrices |')
    w(f'| `efficiency_table_{args.run}.csv` | efficiency pivoted by '
      '(sub_run, mesh, drift) × station |')
    for p in summary:
        w(f'| `{p}` | efficiency / residual width / frame stability per sub_run |')
    w(f'| `per_subrun/` | {len(per_sub)} eight-panel PNGs, one per station per '
      'sub_run: residuals, efficiency maps in both frames, matching distance |')
    w('')

    w('## Results\n')
    piv = df.pivot_table(index=['sub_run', 'mesh_hv', 'drift_hv'],
                         columns='station', values='eff')
    w('```')
    w(piv.round(4).to_string())
    w('```\n')

    w('## Reading these numbers\n')
    zero = df[(df['drift_hv'] - df['mesh_hv']) <= 0]
    if len(zero):
        subs = sorted(set(zero['sub_run']))
        w(f'* **Zero drift field.** In `{", ".join(subs)}` the drift and mesh are at')
        w('  the same voltage for the scanned stations, so there is no field to')
        w(f'  collect ionisation. The ~{zero["eff"].mean():.0%} there is real, not a bug.')
    rc = repeat_check(df)
    if rc:
        w('* **Time drift.** `drift_700` and `nominal_00` are the *same* HV point ~53')
        w('  minutes apart, so the difference is drift with time/dose, not voltage:')
        for st, a, b in rc:
            w(f'  {st} {a:.4f} → {b:.4f} ({(b - a) * 100:+.2f} points).')
        w('  The mesh scan was taken with the voltage *decreasing* over ~2 hours, so')
        w('  any such drift adds coherently to the HV dependence. It is small next to')
        w('  the HV effect but it is not zero.')
    w('* **Systematics.** The statistical errors are Clopper-Pearson at 68.27 % and')
    w('  are tiny at these statistics. The real uncertainty is ~0.2 % from the')
    w('  matching radius and ~0.3 % from the track-quality cut (handoff §13.3).')
    w('* **Acceptance.** Efficiency is measured only where the uRWELL illuminates')
    w('  the pads — a ~130 × 130 mm patch — not over the whole P2.')

    path = os.path.join(args.dir, 'README.md')
    with open(path, 'w') as fh:
        fh.write('\n'.join(L) + '\n')
    print(f'wrote {path}')


if __name__ == '__main__':
    main()
