#!/usr/bin/env python3
"""Merge the per-sub_run condor products of a run into one scan CSV.

The condor pipeline (`p2_condor_urw/urwtb_job.sh`) runs `urw_p2_efficiency.py`
one sub_run at a time and uploads to `analysis/urw_timebins/<run>/<sub_run>/`.
Every HV-scan figure, though, wants the whole ladder in one table, laid out the
way `urw_referenced_efficiency/<run>/urw_p2_efficiency_<run>.csv` already is
for the runs that were processed serially.  This stitches the former into the
latter, so both processing routes end in the same file.

    python3 merge_timebin_products.py --run low_mesh_scan_1 \
        [--src DIR] [--out DIR]
"""

import argparse
import glob
import os

import pandas as pd

DEF_SRC = os.path.expanduser(
    '~/Documents/PostDocSaclay/data/SPS_Beam_Test/mpgd26_workspace/'
    'products/urw_timebins_x/urw_timebins')
DEF_OUT = ('/media/ak271430/LaCie/Extras/Physics/Post-Doc-Saclay/data/'
           'SPS_Beam_Test/TB_July2026_H4/analysis/urw_referenced_efficiency')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', required=True)
    ap.add_argument('--src', default=DEF_SRC)
    ap.add_argument('--out', default=DEF_OUT)
    a = ap.parse_args()

    files = sorted(glob.glob(f'{a.src}/{a.run}/*/urw_p2_efficiency_*.csv'))
    if not files:
        raise SystemExit(f'no per-sub_run CSVs under {a.src}/{a.run}')

    frames = []
    for f in files:
        df = pd.read_csv(f)
        if 'sub_run' not in df.columns:      # older products carried it only
            df['sub_run'] = os.path.basename(os.path.dirname(f))
        frames.append(df)
    d = pd.concat(frames, ignore_index=True)
    d = d.sort_values(['station', 'mesh_hv', 'drift_hv', 'sub_run'])

    out_dir = os.path.join(a.out, a.run)
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, f'urw_p2_efficiency_{a.run}.csv')
    d.to_csv(out, index=False)
    print(f'{len(files)} sub_runs -> {len(d)} rows  {out}')
    for st, g in d.groupby('station'):
        print(f'  {st:<8} mesh {g.mesh_hv.min():.0f}-{g.mesh_hv.max():.0f} V   '
              f'eff {g.eff.min():.3f}-{g.eff.max():.3f}')


if __name__ == '__main__':
    main()
