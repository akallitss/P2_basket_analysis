#!/usr/bin/env python3
"""
Align the two uRWELL references to each other and build two-point tracks.

Both uRWELLs sit on FEU 1, so front and back hits are already in the same
hits_root file and share an eventId - no cross-FEU synchronisation is needed
for uRWELL-only tracking.

Model: translation only.  The detectors are mounted square to the beam, so we
fit a straight line of back-vs-front position per axis, check the slope is
compatible with +1 (that is the "no rotation, no scale" check), and take the
offset as the alignment constant.

Run:
  unset PYTHONPATH
  /local/home/banco/DAQ_Control_Dream_Beam/.venv/bin/python align_and_track.py \
      --run-dir /local/home/banco/P2_data/TB_July2026_H4/runs/drift_mesh_scan_1 \
      --sub-run nominal_00 --plot align.png
"""
import os
import sys
import json
import argparse

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import urw_lib as U  # noqa: E402

RUN_DIR = '/local/home/banco/P2_data/TB_July2026_H4/runs/drift_mesh_scan_1'
FRONT, BACK = 'EIC_uRWELL_front', 'EIC_uRWELL_back'
P2_Z = {'P2_IN': 320.0, 'P2_MID': 630.0, 'P2_OUT': 940.0}


def robust_line(x, y, nsig=2.5, iters=6):
    """Least squares slope/offset with iterative nsig-sigma clipping."""
    keep = np.isfinite(x) & np.isfinite(y)
    slope = offset = np.nan
    for _ in range(iters):
        if keep.sum() < 20:
            break
        slope, offset = np.polyfit(x[keep], y[keep], 1)
        resid = y - (slope * x + offset)
        sigma = resid[keep].std()
        new = np.isfinite(x) & np.isfinite(y) & (np.abs(resid) < nsig * sigma)
        if new.sum() == keep.sum():
            keep = new
            break
        keep = new
    resid = y - (slope * x + offset)
    return dict(slope=float(slope), offset=float(offset),
                rms=float(resid[keep].std()), n_used=int(keep.sum()),
                frac_used=float(keep.mean())), keep


def leading_per_event(det, sub_run_dir, run_json, min_amp, max_hits):
    geo = U.UrwGeometry(det, run_json,
                        sub_run_name=os.path.basename(sub_run_dir))
    parts = []
    for chunk in U.iter_hits(U.feu_hit_files(sub_run_dir, feu=geo.feu_num),
                             max_hits=max_hits, progress=False,
                             feu=geo.feu_num):
        clusters = U.cluster_hits(chunk, geo, min_amp=min_amp)
        if len(clusters):
            parts.append(U.leading_points(clusters))
    pts = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    return geo, pts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run-dir', default=RUN_DIR)
    ap.add_argument('--sub-run', default='nominal_00')
    ap.add_argument('--min-amp', type=float, default=100.0)
    ap.add_argument('--max-hits', type=float, default=3e6)
    ap.add_argument('--plot', default='')
    ap.add_argument('--save-tracks', default='')
    args = ap.parse_args()

    run_json = os.path.join(args.run_dir, 'run_config.json')
    sub_dir = os.path.join(args.run_dir, args.sub_run)

    geo_f, pf = leading_per_event(FRONT, sub_dir, run_json,
                                  args.min_amp, int(args.max_hits))
    geo_b, pb = leading_per_event(BACK, sub_dir, run_json,
                                  args.min_amp, int(args.max_hits))
    dz = float(geo_b.center[2] - geo_f.center[2])
    print(f'front: {len(pf)} events with a leading cluster')
    print(f'back : {len(pb)} events with a leading cluster')
    print(f'lever arm dz = {dz:.1f} mm')

    ev = pf[['eventId', 'x', 'y']].rename(columns={'x': 'xf', 'y': 'yf'}).merge(
        pb[['eventId', 'x', 'y']].rename(columns={'x': 'xb', 'y': 'yb'}),
        on='eventId')
    both = ev.dropna(subset=['xf', 'xb', 'yf', 'yb']).reset_index(drop=True)
    print(f'matched on eventId: {len(ev)}, with all four coordinates: {len(both)}')

    align = {}
    for axis in ('x', 'y'):
        fit, keep = robust_line(both[f'{axis}f'].values, both[f'{axis}b'].values)
        align[axis] = fit
        print(f'\n{axis}: back = {fit["slope"]:+.4f} * front {fit["offset"]:+.3f} mm'
              f'   spread {fit["rms"]:.2f} mm'
              f'   ({fit["n_used"]} of {len(both)} used, {fit["frac_used"]:.1%})')
        print(f'   slope - 1 = {fit["slope"] - 1:+.4f}  '
              f'-> {"consistent with no rotation/scale" if abs(fit["slope"] - 1) < 0.05 else "CHECK: not 1"}')

    # translation-only alignment: shift the back plane onto the front
    for axis in ('x', 'y'):
        both[f'{axis}b_al'] = both[f'{axis}b'] - align[axis]['offset']
        both[f'd{axis}'] = both[f'{axis}b_al'] - both[f'{axis}f']

    # two-point track: position and slope, extrapolate to the P2 planes
    for axis in ('x', 'y'):
        both[f'slope_{axis}'] = (both[f'{axis}b_al'] - both[f'{axis}f']) / dz
    for name, z in P2_Z.items():
        for axis in ('x', 'y'):
            both[f'{name}_{axis}'] = both[f'{axis}f'] + both[f'slope_{axis}'] * z

    print('\ntrack angles (mrad), after alignment:')
    summary = {'run': os.path.basename(args.run_dir), 'sub_run': args.sub_run,
               'min_amp': args.min_amp, 'dz_mm': dz,
               'n_matched': int(len(both)), 'align': align, 'angles': {}}
    for axis in ('x', 'y'):
        a = both[f'slope_{axis}'].values * 1000.0
        med, iqr = np.median(a), np.percentile(a, 75) - np.percentile(a, 25)
        summary['angles'][axis] = dict(median=float(med), rms=float(a.std()),
                                       sigma_from_iqr=float(iqr / 1.349))
        print(f'  {axis}: median {med:+.2f}  rms {a.std():.2f}  '
              f'sigma(IQR) {iqr / 1.349:.2f}')

    print('\nextrapolated beam position at the P2 planes (mm, front local frame):')
    for name, z in P2_Z.items():
        print(f'  {name} (z={z:.0f}): '
              f'x {both[f"{name}_x"].mean():6.1f} +- {both[f"{name}_x"].std():5.1f}   '
              f'y {both[f"{name}_y"].mean():6.1f} +- {both[f"{name}_y"].std():5.1f}')

    out_json = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            f'alignment_{os.path.basename(args.run_dir)}_{args.sub_run}.json')
    with open(out_json, 'w') as fh:
        json.dump(summary, fh, indent=2)
    print(f'\nwrote {out_json}')

    if args.save_tracks:
        # Format follows the extension.  Parquet needs pyarrow or fastparquet
        # and NEITHER is installed in the venv on this machine, so a .parquet
        # name fails after all the work is done - .npz is the safe default.
        ext = os.path.splitext(args.save_tracks)[1].lower()
        if ext == '.parquet':
            both.to_parquet(args.save_tracks)
        elif ext == '.csv':
            both.to_csv(args.save_tracks, index=False)
        else:
            np.savez_compressed(args.save_tracks,
                                **{c: both[c].values for c in both.columns})
        print(f'wrote {args.save_tracks} ({len(both)} tracks)')

    if args.plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2, 3, figsize=(16, 9))
        for row, axis in enumerate(('x', 'y')):
            ax[row, 0].hist2d(both[f'{axis}f'], both[f'{axis}b'], bins=120, cmin=1)
            ax[row, 0].set(xlabel=f'front {axis} [mm]', ylabel=f'back {axis} [mm]',
                           title=f'{axis}: back vs front (slope {align[axis]["slope"]:+.3f})')
            ax[row, 1].hist(both[f'd{axis}'], bins=200, range=(-40, 40))
            ax[row, 1].set(xlabel=f'back - front, aligned [mm]',
                           title=f'{axis} spread rms {align[axis]["rms"]:.2f} mm')
            ax[row, 2].hist(both[f'slope_{axis}'] * 1000, bins=200, range=(-40, 40))
            ax[row, 2].set(xlabel=f'track angle {axis} [mrad]', title=f'{axis} angle')
        fig.suptitle(f'{os.path.basename(args.run_dir)}/{args.sub_run}  '
                     f'{len(both)} matched events, amplitude > {args.min_amp:.0f}')
        fig.tight_layout()
        fig.savefig(args.plot, dpi=100)
        print(f'wrote {args.plot}')


if __name__ == '__main__':
    main()
