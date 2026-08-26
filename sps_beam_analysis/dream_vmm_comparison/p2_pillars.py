#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""p2_pillars.py -- the fine-grained pass: dead-area masking and the bulk
pillar lattice.

`p2_selftrack.py` measured what the three stations do over the whole beam
spot.  Two questions it cannot answer come out of the same tracks:

  1. **How does the stack do where it is alive?**  Part of the inefficiency is
     hardware -- a weak connector on P2_IN, the five 6.15 mm bulk pillars, pads
     that never respond.  Those are known, localised, and would be QC'd out of a
     production detector, so the stack's performance with them masked is a
     different (and more transferable) number from its performance with them in.
     This stage keeps a per-pad ledger so any mask can be applied offline, and
     writes the *whole* joined table rather than a sample so the masked widths
     are measured, not extrapolated.

  2. **Can the reference see the bulk pillars?**  The CERN bulk mask
     (`P2_Mask2.gbr`) puts 41 366 pillars of diameter 0.500 mm on an exact
     2.000 mm square lattice, at even-integer millimetre coordinates of the
     board frame -- 4.9 % of the amplification area, with *zero* free
     parameters in the predicted phase.  The pad frame this stage bins in is
     that board frame, so the pillar signal is a modulation at a known
     wavevector.  Whether it is visible is a question about the reference
     telescope's pointing resolution, and it is answered by a number: the
     amplitude of the modulation at k = (pi, 0) and (0, pi) rad/mm.

Both need a map far finer than the pad, so this stage carries a 0.20 mm map --
ten bins per pillar period -- over a fixed window around the beam spot.  It is
counts, so sub-runs add; everything else is done offline in `pillar_stats.py`.

Nothing about the reference, the alignment, the fiducial or the vetoes is
re-derived: `build_joined` is imported from `p2_selftrack`, so the two stages
share one denominator.

Usage (lxplus, LCG_110, see run_pillars.sh)
    python3 p2_pillars.py --run eff_nominal_1 --out <dir>
"""
import os
import sys
import json
import argparse

import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import p2_selftrack as ST        # noqa: E402
import sps_config as sc          # noqa: E402

STATIONS = ST.STATIONS

# The fine map.  0.20 mm is ten bins per pillar period, so the lattice
# frequency lands exactly on a bin of the map's own FFT and the bin-width
# correction is an analytic sinc, not a fudge.  +-64 mm covers the beam spot
# (sigma ~ 22 mm) and is exactly 64 pillar periods across, so a window edge is
# never half a period.
FINE_BIN = 0.20
FINE_HALF = 64.0
PILLAR_PITCH = 2.0               # mm, exact, from P2_Mask2.gbr
PILLAR_DIA = 0.500               # mm, exact, from P2_Mask2.gbr

# The residual to the reference, taken along the pad's own axes and binned
# finely enough to resolve the EDGE of the box.  A one-pad cluster reports the
# pad centre, so this residual is (pad box) convolved with the reference's
# pointing error and nothing else: the edge width IS the pointing resolution,
# measured on the same tracks, in the same frame, with no model of the
# telescope in between.  0.02 mm bins over +-12 mm.
EDGE_EDGES = np.arange(-12.0, 12.0 + 1e-9, 0.02)


def fine_edges(centre):
    n = int(round(2 * FINE_HALF / FINE_BIN))
    return centre - FINE_HALF + FINE_BIN * np.arange(n + 1)


def snap(v):
    """Centre the fine window on a pillar node, so the window spans a whole
    number of periods and the fold has no edge phase."""
    return float(PILLAR_PITCH * np.round(v / PILLAR_PITCH))


def fine_maps(m, s, centre):
    """n / found / amplitude / cluster-size, binned at 0.20 mm."""
    xe, ye = fine_edges(centre[0]), fine_edges(centre[1])
    fid = m[f'fid_{s}'].to_numpy()
    fnd = fid & m[f'found_{s}'].to_numpy()
    ex, ey = m[f'ex_{s}'].to_numpy(), m[f'ey_{s}'].to_numpy()
    out = {f'fine_n_{s}': np.histogram2d(ex[fid], ey[fid], [xe, ye])[0],
           f'fine_k_{s}': np.histogram2d(ex[fnd], ey[fnd], [xe, ye])[0]}
    a = m[f'a_lead_{s}'].to_numpy()
    wa = fnd & np.isfinite(a)
    out[f'fine_asum_{s}'] = np.histogram2d(ex[wa], ey[wa], [xe, ye],
                                           weights=a[wa])[0]
    out[f'fine_an_{s}'] = np.histogram2d(ex[wa], ey[wa], [xe, ye])[0]
    nc = m[f'n_clus_{s}'].to_numpy()
    wc = fnd & np.isfinite(nc)
    out[f'fine_csum_{s}'] = np.histogram2d(ex[wc], ey[wc], [xe, ye],
                                           weights=nc[wc])[0]
    return out


def edge_hists(m, s):
    """The pad-box residual, split by cluster size.

    Only the one-pad clusters carry the clean box: a two-pad centroid sits
    somewhere between two centres and rounds the edge for a reason that has
    nothing to do with the reference.
    """
    out = {}
    fnd = m[f'fid_{s}'].to_numpy() & m[f'found_{s}'].to_numpy()
    nc = m[f'n_clus_{s}'].to_numpy()
    for ax in ('rw', 'rh'):
        d = m[f'{ax}_{s}'].to_numpy()
        for tag, sel in (('single', nc <= 1), ('multi', nc >= 2)):
            w = fnd & sel & np.isfinite(d)
            out[f'edge_{ax}_{s}_{tag}'] = np.histogram(d[w], EDGE_EDGES)[0]
    return out


def pad_ledger(m, s, npad):
    """Per pad: tracks pointed at it, tracks it found, summed lead amplitude.

    Indexed by channel_id, so a mask is a list of ids and nothing has to be
    matched by position.
    """
    fid = m[f'fid_{s}'].to_numpy()
    fnd = fid & m[f'found_{s}'].to_numpy()
    pid = np.nan_to_num(m[f'pid_{s}'].to_numpy(), nan=-1).astype(int)
    ok = fid & (pid >= 0) & (pid < npad)
    a = m[f'a_lead_{s}'].to_numpy()
    wa = fnd & np.isfinite(a) & (pid >= 0) & (pid < npad)
    return {f'pad_n_{s}': np.bincount(pid[ok], minlength=npad),
            f'pad_k_{s}': np.bincount(pid[fnd & ok], minlength=npad),
            f'pad_an_{s}': np.bincount(pid[wa], minlength=npad),
            f'pad_asum_{s}': np.bincount(pid[wa], a[wa], npad)}


SAMPLE_COLS = ['xf', 'yf', 'slope_x', 'slope_y'] + [
    f'{c}_{s}' for s in STATIONS
    for c in ('ex', 'ey', 'found', 'fid', 'ux', 'uy', 'a_lead', 'n_clus',
              'pid', 'qx', 'qy', 'dx', 'dy', 'rw', 'rh')]


def run_subrun(cfg, dets, sub_run, args, run_dir, run_json, centres):
    got = ST.build_joined(cfg, dets, sub_run, args, run_dir, run_json)
    if got is None:
        return None
    m, tracks, tinfo, zs, fits, exts = got
    npad = args.npad

    out = {}
    res = {'run': args.run, 'sub_run': sub_run, 'z_mm': zs.tolist(),
           'n_tracks': int(len(tracks)), 'urwell': tinfo,
           'frame': {s: fits[s] for s in STATIONS}, 'stations': {}}

    for s in STATIONS:
        if s not in centres:
            # freeze the window on the first sub_run that has tracks, so every
            # sub_run adds into the same bins
            centres[s] = (snap(np.median(m[f'ex_{s}'])),
                          snap(np.median(m[f'ey_{s}'])))
            print(f'  {s}: fine window centred on {centres[s]}')
        out.update(fine_maps(m, s, centres[s]))
        out.update(pad_ledger(m, s, npad))
        out.update(edge_hists(m, s))
        fid = m[f'fid_{s}'].to_numpy()
        fnd = fid & m[f'found_{s}'].to_numpy()
        res['stations'][s] = {'n_fid': int(fid.sum()), 'n_found': int(fnd.sum()),
                              'efficiency': float(fnd.sum() / max(fid.sum(), 1)),
                              'centre': list(centres[s])}
        print(f'  {s}: eff {res["stations"][s]["efficiency"]:.4f}  '
              f'({int(fnd.sum())}/{int(fid.sum())})')

    fid_all = np.ones(len(m), bool)
    for s in STATIONS:
        fid_all &= m[f'fid_{s}'].to_numpy()
    sub = m[fid_all]
    if args.sample and len(sub) > args.sample:
        sub = sub.iloc[::int(np.ceil(len(sub) / args.sample))]
    out['sample'] = sub[SAMPLE_COLS].to_numpy(np.float32)
    res['n_fid_all3'] = int(fid_all.sum())
    res['n_sample'] = int(len(sub))
    sys.stdout.flush()
    return res, out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', default='eff_nominal_1')
    ap.add_argument('--data-root', default=ST.E.DATA_ROOT)
    ap.add_argument('--sub-run', action='append', default=[])
    ap.add_argument('--max-chunks', type=int, default=0)
    ap.add_argument('--track-cut', type=float, default=3.0)
    ap.add_argument('--cluster-r', type=float, default=15.0)
    ap.add_argument('--probe-r', type=float, default=15.0)
    ap.add_argument('--fid-r', type=float, default=9.0)
    ap.add_argument('--match-win', type=float, default=25.0)
    ap.add_argument('--min-events', type=int, default=500)
    ap.add_argument('--no-veto-sparks', action='store_true')
    ap.add_argument('--npad', type=int, default=1400,
                    help='length of the per-pad ledger (channel_id space)')
    ap.add_argument('--sample', type=int, default=0,
                    help='0 = keep every track fiducial in all three')
    ap.add_argument('--out', default='out_pillars')
    args = ap.parse_args()

    run_dir = os.path.join(args.data_root, args.run)
    run_json = os.path.join(run_dir, 'run_config.json')
    cfg = sc.RunConfig(key='urw_p2', run=args.run, data_root=args.data_root)
    os.makedirs(args.out, exist_ok=True)

    dets = [d for d in cfg.mappable_detectors() if d.name in STATIONS]
    dets.sort(key=lambda d: d.z)
    if len(dets) != 3:
        sys.exit(f'need all three stations, found {[d.name for d in dets]}')
    sub_runs = args.sub_run or cfg.find_subruns()
    print(f'{args.run}: {len(sub_runs)} sub_run(s); fine bin {FINE_BIN} mm, '
          f'pillar pitch {PILLAR_PITCH} mm dia {PILLAR_DIA} mm')

    summaries, acc, centres = [], {}, {}
    for sub_run in sub_runs:
        try:
            got = run_subrun(cfg, dets, sub_run, args, run_dir, run_json,
                             centres)
        except Exception as exc:                      # noqa: BLE001
            print(f'  !! {sub_run} failed: {exc!r}')
            continue
        if got is None:
            continue
        res, out = got
        summaries.append(res)
        for k, v in out.items():
            if k == 'sample':
                acc.setdefault(k, []).append(v)
            else:
                acc[k] = acc[k] + v if k in acc else v
        _write(summaries, acc, args, centres)
    _write(summaries, acc, args, centres)
    print(f'\ndone: {len(summaries)} sub_run(s)')


def _write(summaries, acc, args, centres):
    base = os.path.join(args.out, f'p2_pillars_{args.run}')
    with open(base + '.json', 'w') as fh:
        json.dump(summaries, fh, indent=1)
    out = {}
    for k, v in acc.items():
        if k == 'sample':
            out[k] = np.concatenate(v)
        elif k.startswith('fine_n_') or k.startswith('fine_k_') \
                or k.startswith('fine_an_') or k.startswith('pad_') \
                or k.startswith('edge_'):
            out[k] = np.asarray(v).astype(np.int32)
        else:
            out[k] = np.asarray(v).astype(np.float32)
    out['sample_cols'] = np.array(SAMPLE_COLS)
    out['fine_bin'] = np.array([FINE_BIN])
    out['fine_half'] = np.array([FINE_HALF])
    out['pillar'] = np.array([PILLAR_PITCH, PILLAR_DIA])
    out['edge_edges'] = EDGE_EDGES
    for s in STATIONS:
        if s in centres:
            out[f'centre_{s}'] = np.array(centres[s])
    np.savez_compressed(base + '.npz', **out)
    print(f'  wrote {base}.json / .npz '
          f'({os.path.getsize(base + ".npz") / 1e6:.1f} MB)')


if __name__ == '__main__':
    main()
