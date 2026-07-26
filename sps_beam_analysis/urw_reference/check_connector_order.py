#!/usr/bin/env python3
"""
Determine, from data, whether each uRWELL view's two Dream connectors are
assigned to the strip map in the right order.

Method (no free parameters):
  The combined_hits_root file is event-synchronised across all four FEUs, so for
  every trigger we can take the leading (highest-amplitude) uRWELL channel in a
  view and the leading P2_MID pad, whose position comes from the independently
  validated P2_BASKET pad map.  Plotting the median P2 pad coordinate against
  the uRWELL channel index must give a single continuous monotone relation:
  the beam is one spot, and both detectors see the same tracks.

  If the two 64-channel Dream connectors of a view are interchanged with respect
  to the map, that relation breaks into two disjoint branches with a jump at the
  connector boundary.  Interchanging them, (idx + 64) % 128, repairs it.  The
  test is which of the two orderings has the higher |Spearman rho|.

Result on drift_mesh_scan_1/nominal_00 (2026-07-25):
  front_x swap, front_y keep, back_x swap, back_y swap

SUPERSEDED 2026-07-26 -- like check_view_transform.py, this test is blind to a
reversal of channel order INSIDE a connector, and that is what the back has.
Its back rows are the wrong member of a degenerate pair.  Use
explore4_back_map.py; the live table is urw_lib.VIEW_MODE_DEFAULT.

Usage:
  PY=/local/home/banco/DAQ_Control_Dream_Beam/.venv/bin/python
  unset PYTHONPATH
  $PY check_connector_order.py [--sub-run nominal_00] [--plot out.png]
"""

import os
import sys
import glob
import argparse

import numpy as np
import pandas as pd
import uproot

# this package, then sps_beam_analysis for sps_config, whose setup_paths()
# adds the shared core in cosmic_bench_analysis (p2_io, p2_mapping, ...)
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.dirname(_HERE))
import sps_config as _sc  # noqa: E402
_sc.setup_paths()
import p2_mapping as pmap  # noqa: E402
import sps_config as sc  # noqa: E402

RUN_DIR = '/local/home/banco/P2_data/TB_July2026_H4/runs/drift_mesh_scan_1'
VIEWS = {'front_x': (0, 128), 'front_y': (128, 256),
         'back_x': (256, 384), 'back_y': (384, 512)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run-dir', default=RUN_DIR)
    ap.add_argument('--sub-run', default='nominal_00')
    ap.add_argument('--anchor', default='P2_MID')
    ap.add_argument('--min-amp', type=float, default=150.0)
    ap.add_argument('--plot', default='')
    args = ap.parse_args()

    run_json = os.path.join(args.run_dir, 'run_config.json')
    combined = sorted(glob.glob(os.path.join(
        args.run_dir, args.sub_run, 'combined_hits_root', '*_feu-combined_hits.root')))
    if not combined:
        sys.exit(f'no combined_hits_root file in {args.run_dir}/{args.sub_run}')
    ct = pmap.build_channel_table(run_json, sc.MAP_CSV_PATH, det_type='P2',
                                  det_name=args.anchor)
    anchor_feu = ct.attrs['feus'][0]

    parts = []
    for chunk in uproot.open(combined[0])['hits'].iterate(
            ['eventId', 'channel', 'amplitude', 'feu'], library='np', step_size=4_000_000):
        d = pd.DataFrame(chunk)
        parts.append(d[d['amplitude'] > args.min_amp])
    df = pd.concat(parts, ignore_index=True)
    print(f'{os.path.basename(combined[0])}: {len(df)} hits with amp > {args.min_amp}')

    tab = ct[['channel', 'pad_cx', 'pad_cy']].set_index('channel')
    s = df[df['feu'] == anchor_feu].copy()
    s[['pad_cx', 'pad_cy']] = tab.reindex(s['channel']).values
    p2 = (s.dropna(subset=['pad_cx']).sort_values('amplitude', ascending=False)
          .groupby('eventId').first()[['pad_cx', 'pad_cy']])
    print(f'{args.anchor} (feu {anchor_feu}): {len(p2)} events with a mapped pad')

    u = df[df['feu'] == 1]
    rows, panels = [], []
    for vn, (lo, hi) in VIEWS.items():
        sv = u[(u['channel'] >= lo) & (u['channel'] < hi)]
        idx = (sv.sort_values('amplitude', ascending=False)
               .groupby('eventId').first()['channel'] - lo).rename('idx')
        j = pd.concat([idx, p2], axis=1, join='inner').dropna()
        best = max(['pad_cx', 'pad_cy'],
                   key=lambda c: abs(j['idx'].corr(j[c], method='spearman')))
        k_nat = j['idx'].values
        k_swap = (k_nat + 64) % 128
        rho_nat = pd.Series(k_nat).corr(pd.Series(j[best].values), method='spearman')
        rho_swap = pd.Series(k_swap).corr(pd.Series(j[best].values), method='spearman')
        verdict = 'SWAP' if abs(rho_swap) > abs(rho_nat) else 'keep'
        rows.append(dict(view=vn, anchor_axis=best, n=len(j),
                         rho_natural=rho_nat, rho_swapped=rho_swap, verdict=verdict))
        panels.append((vn, best, k_nat, k_swap, j[best].values))

    out = pd.DataFrame(rows)
    print()
    print(out.to_string(index=False,
                        formatters={'rho_natural': '{:+.3f}'.format,
                                    'rho_swapped': '{:+.3f}'.format}))
    print('\nCONNECTOR_SWAP =', {vn: (v == 'SWAP') for vn, v in
                                 zip(out['view'], out['verdict'])})

    if args.plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(2, 4, figsize=(19, 8))
        for col, (vn, best, k_nat, k_swap, pv) in enumerate(panels):
            for row, (k, lab) in enumerate([(k_nat, 'natural'), (k_swap, 'conn-swapped')]):
                ax = axes[row, col]
                g = pd.DataFrame({'k': k, 'p': pv}).groupby('k')['p']
                med, cnt = g.median(), g.size()
                m = cnt >= 50
                ax.plot(med.index[m], med[m], '.', ms=4)
                ax.axvline(64, color='r', ls='--', lw=0.8)
                rho = pd.Series(k).corr(pd.Series(pv), method='spearman')
                ax.set_title(f'{vn} [{lab}] vs {best}\nspearman={rho:+.3f}', fontsize=9)
                ax.set_xlabel('channel idx')
                ax.set_ylabel(f'median {args.anchor} {best} [mm]')
        fig.suptitle('uRWELL channel -> position continuity test against the P2 pad map')
        fig.tight_layout()
        fig.savefig(args.plot, dpi=100)
        print(f'wrote {args.plot}')


if __name__ == '__main__':
    main()
