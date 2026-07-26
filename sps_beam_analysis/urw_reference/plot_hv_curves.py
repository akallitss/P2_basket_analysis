#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_hv_curves.py -- P2 efficiency versus HV, from urw_p2_efficiency.py output.

Reads the CSV that `urw_p2_efficiency.py` writes for a scan run and turns it
into efficiency curves, in the style of
`analysis/<det>/<run>/scan/22_tag_probe_efficiency/tag_probe_efficiency*.png`:
one panel, one coloured line per probe station, Clopper-Pearson error bars.

Grouping the points is the only subtle part
-------------------------------------------
A "mesh scan" does not hold the drift electrode fixed.  In `drift_mesh_scan_1`
the drift tracks the mesh so that the drift FIELD stays constant -- P2_MID and
P2_OUT run at drift = mesh + 250 V, P2_IN at mesh + 200 V -- while the separate
`drift_*` sub-runs hold the mesh at 450 V and walk the drift from 450 to 900 V.
So the invariant that defines a scan family is **drift - mesh**, not drift.

Each station's points are therefore grouped by `drift - mesh`; a group holding
more than one mesh value is a mesh-scan family.  Families are then matched
across stations by rank (each station's lowest drift-field family together,
etc.), which is what produces one figure per drift setting for a genuine 2D
scan and a single figure for a 1D one.  Whatever ends up in a figure is printed
and written into its caption, so the grouping is auditable rather than implied.

The drift scan gets its own figure the same way, grouped by mesh.

Usage
  unset PYTHONPATH
  $PY plot_hv_curves.py --csv <analysis>/urw_p2_efficiency_drift_mesh_scan_1.csv \
      --out <analysis>
"""
import os
import argparse
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

COLOURS = {'P2_IN': '#1f77b4', 'P2_MID': '#ff7f0e', 'P2_OUT': '#2ca02c'}
ORDER = ['P2_IN', 'P2_MID', 'P2_OUT']


def families(df, scan_col, fixed_col, min_points=3):
    """{station: [(fixed_value, sub-frame sorted by scan_col), ...]}.

    A family is a set of points at one value of `drift - mesh` (when scanning
    the mesh) or one mesh (when scanning the drift), holding at least
    `min_points` distinct values of the scanned electrode.
    """
    out = {}
    for st, g in df.groupby('station'):
        fams = []
        for key, gg in g.groupby(fixed_col):
            gg = gg.sort_values(scan_col)
            if gg[scan_col].nunique() >= min_points:
                fams.append((key, gg))
        if fams:
            out[st] = sorted(fams, key=lambda kv: kv[0])
    return out


def by_key(keys, fmt):
    """'200 V (P2_IN), 250 V (P2_MID, P2_OUT)' from {station: key}.

    Stations sharing a setting are listed together, which keeps the caption
    short enough to read when two of the three run at the same value.
    """
    inv = defaultdict(list)
    for st, k in keys.items():
        inv[k].append(st)
    return ', '.join(
        f'{fmt(k)} ({", ".join(sorted(v, key=lambda s: ORDER.index(s) if s in ORDER else 99))})'
        for k, v in sorted(inv.items()))


def merge_duplicates(g, scan_col):
    """One point per HV setting, counts combined, interval recomputed.

    A setting can appear in several sub_runs - P2_IN holds mesh 430 V through
    the whole drift scan, and nominal_00 repeats drift_700 exactly - and those
    are genuine repeats of the same point, so they are added rather than drawn
    on top of each other.  `n_sub_runs` records how many went in.
    """
    from scipy.stats import beta

    def cp(k, n, cl=0.6827):
        if n == 0:
            return 0.0, 1.0
        a = 1.0 - cl
        return (float(beta.ppf(a / 2, k, n - k + 1)) if k > 0 else 0.0,
                float(beta.ppf(1 - a / 2, k + 1, n - k)) if k < n else 1.0)

    rows = []
    for x, gg in g.groupby(scan_col):
        n = int(gg['n'].sum())
        k = int(round((gg['eff'] * gg['n']).sum()))
        lo, hi = cp(k, n)
        rows.append({scan_col: x, 'n': n, 'eff': k / max(n, 1), 'lo': lo,
                     'hi': hi, 'n_sub_runs': len(gg)})
    return pd.DataFrame(rows).sort_values(scan_col)


def plot_family(groups, scan_col, xlabel, title, subtitle, out_png):
    fig, ax = plt.subplots(figsize=(9.5, 6))
    n_tot = 0
    for st in ORDER:
        if st not in groups:
            continue
        key, g = groups[st]
        m = merge_duplicates(g, scan_col)
        x = m[scan_col].to_numpy(float)
        y = m['eff'].to_numpy(float)
        lo = y - m['lo'].to_numpy(float)
        hi = m['hi'].to_numpy(float) - y
        rep = m['n_sub_runs'].max()
        lbl = f'probe {st}  ({key})'
        if rep > 1:
            lbl += f'  [up to {rep} sub_runs merged]'
        ax.errorbar(x, y, yerr=[lo, hi], fmt='o-', ms=6, lw=1.8, capsize=3,
                    color=COLOURS.get(st), label=lbl)
        n_tot += len(m)
    ax.set_xlabel(xlabel)
    ax.set_ylabel('efficiency (uRWELL-track referenced)')
    ax.set_ylim(0.0, 1.02)
    ax.grid(alpha=0.3)
    ax.legend(loc='lower right', fontsize=9)
    import textwrap
    ax.set_title(f'{title}\n' + '\n'.join(textwrap.wrap(subtitle, 88)),
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(out_png, dpi=130)
    plt.close(fig)
    print(f'wrote {out_png}  ({n_tot} points)')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--min-points', type=int, default=3)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    df = pd.read_csv(args.csv)
    df = df[df['mesh_hv'].notna() & df['drift_hv'].notna()].copy()
    df['dfield'] = (df['drift_hv'] - df['mesh_hv']).astype(int)
    run = df['run'].iloc[0]
    print(f'{run}: {len(df)} station-points, '
          f'{df["sub_run"].nunique()} sub_runs')

    made = []

    # ---- mesh scans: one figure per drift-field setting -------------------- #
    fam = families(df, 'mesh_hv', 'dfield', args.min_points)
    n_fig = max((len(v) for v in fam.values()), default=0)
    for i in range(n_fig):
        groups = {st: v[i] for st, v in fam.items() if len(v) > i}
        if not groups:
            continue
        keys = {st: g[0] for st, g in groups.items()}
        sub = 'drift = mesh + ' + by_key(keys, lambda k: f'{k} V')
        tag = '_'.join(str(k) for k in sorted(set(keys.values())))
        png = os.path.join(args.out, f'efficiency_vs_mesh_dfield{tag}.png')
        plot_family(groups, 'mesh_hv', 'mesh HV [V]',
                    f'P2 efficiency vs mesh voltage — {run}',
                    f'referenced to uRWELL tracks; {sub}', png)
        made.append(png)
        for st, (k, g) in groups.items():
            m = merge_duplicates(g, 'mesh_hv')
            print(f'   {st} (drift-mesh = {k} V): ' +
                  ', '.join(f'{int(v)}V:{e:.4f}'
                            for v, e in zip(m['mesh_hv'], m['eff'])))

    # ---- drift scans: one figure per mesh setting -------------------------- #
    fam = families(df, 'drift_hv', 'mesh_hv', args.min_points)
    n_fig = max((len(v) for v in fam.values()), default=0)
    for i in range(n_fig):
        groups = {st: v[i] for st, v in fam.items() if len(v) > i}
        if not groups:
            continue
        keys = {st: int(g[0]) for st, g in groups.items()}
        sub = 'mesh held at ' + by_key(keys, lambda k: f'{k} V')
        tag = '_'.join(str(k) for k in sorted(set(keys.values())))
        png = os.path.join(args.out, f'efficiency_vs_drift_mesh{tag}.png')
        plot_family(groups, 'drift_hv', 'drift HV [V]',
                    f'P2 efficiency vs drift voltage — {run}',
                    f'referenced to uRWELL tracks; {sub}', png)
        made.append(png)

    if not made:
        print('no scan family had enough points; nothing plotted')

    # ---- the table behind the curves --------------------------------------- #
    piv = df.pivot_table(index=['sub_run', 'mesh_hv', 'drift_hv'],
                         columns='station', values='eff')
    tab = os.path.join(args.out, f'efficiency_table_{run}.csv')
    piv.round(4).to_csv(tab)
    print(f'wrote {tab}')
    print(piv.round(4).to_string())


if __name__ == '__main__':
    main()
