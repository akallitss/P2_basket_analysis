#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""urw_p2_padadc.py -- per-pad pulse height and efficiency on the DREAM readout.

The VMM-side twin of this is `pad_pulse_height()` in `eff_autopsy_report.py`:
for every pad, the efficiency of the uRWELL tracks that point at it and the
pulse-height distribution of the ones it actually recorded.  This runs the same
measurement on the DREAM readout of the same three detectors, so the two
per-pad gain maps can be put side by side.

Nothing about the reference or the selection is re-invented here: tracks, the
frame fit, the fiducial cut and the probe radius all come from
`urw_p2_efficiency.run_station`'s recipe with its defaults, because the point
of the comparison is that only the electronics differ.

What is added is the pad INDEX of the track prediction (`tree.query` returns it
and the parent drops it), which turns a station-level efficiency into a per-pad
one, and a per-pad histogram of the leading-pad amplitude `a_lead` of the
tracks the station DID record -- the DREAM analogue of the VMM's `win_adc`.

Usage (lxplus, LCG_110, see run_padadc.sh)
    python3 urw_p2_padadc.py --run eff_nominal_1 --station P2_OUT \
        --sub-run eff_nominal_00 --out $HOME/p2_eff_dn/out/padadc
"""
import os
import sys
import json
import glob
import argparse

import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import urw_p2_efficiency as E     # noqa: E402  (brings sps_config paths with it)
import sps_config as sc           # noqa: E402
import sps_cluster as scl         # noqa: E402
import p2_sparks as ps            # noqa: E402

# a_lead histogram: DREAM amplitudes run to a few thousand ADC (P2_OUT cluster
# MPV is ~274), so 8-ADC bins over 0..4096 resolve the peak and still hold the
# tail.  Same shape as the VMM side's adc_vs_ch, one row per pad.
AMP_BIN = 8.0
AMP_NBIN = 512


def pad_table(cfg, det):
    ct = cfg.channel_table(det)
    pads = ct[ct['mapped']].drop_duplicates('channel_id').reset_index(drop=True)
    return pads


def run_subrun(cfg, det, sub_run, args, acc):
    """Accumulate one sub_run into `acc` (per-pad counters + histograms)."""
    from scipy.spatial import cKDTree
    run_dir = os.path.join(args.data_root, args.run)
    run_json = os.path.join(run_dir, 'run_config.json')
    sub_dir = os.path.join(run_dir, sub_run)

    tracks, tinfo, dz = E.build_tracks(run_json, sub_dir,
                                       args.max_chunks or None, args.track_cut)
    hits_dir = cfg.combined_hits_dir(sub_run)
    files = sorted(glob.glob(os.path.join(hits_dir, '*.root')))
    if args.max_chunks:
        files = files[:args.max_chunks]
    if not files:
        print(f'  {sub_run}: no combined hits'), sys.stdout.flush()
        return None
    hv_csv = cfg.hv_monitor_csv(sub_run)

    # same denominator hygiene as run_station, in the same order
    t_span = float(tracks['t_ns'].max()) / 1e9 if len(tracks) else 0.0
    t_min = scl.settle_t_min(hv_csv, det.spark_channel, files[0])
    if t_min > 0 and t_min >= t_span:
        t_min = 0.0
    veto = None
    if not args.no_veto_sparks and os.path.isfile(hv_csv) and det.spark_channel:
        class _Shim:
            SPARK_CHANNEL = det.spark_channel
            SPARK_IMON_THR = cfg.SPARK_IMON_THR
            SPARK_GUARD_BEFORE = cfg.SPARK_GUARD_BEFORE
            SPARK_GUARD_AFTER = cfg.SPARK_GUARD_AFTER
            BURST_NPADS = cfg.BURST_NPADS
        veto = ps.SparkVeto.from_csv(hv_csv, _Shim)

    pads = acc['pads']
    pads_id = pads['channel_id'].to_numpy()
    tree = cKDTree(pads[['pad_cx', 'pad_cy']].to_numpy())
    ct = cfg.channel_table(det)
    p2 = E.p2_clusters(files, ct, args.cluster_r, min_amp=cfg.MIN_AMP)

    tr = tracks.copy()
    rec = (cfg.recorded_events(sub_run) or {}).get(det.feus[0])
    if rec is not None:
        tr = tr[E.recorded_mask(tr['eventId'], rec)]
    if t_min > 0:
        tr = tr[tr['t_ns'] / 1e9 >= t_min]
    if veto is not None and veto.intervals:
        tr = tr[veto.event_mask(tr['t_ns'].to_numpy())]

    src = E.project(tr, det.z, dz)

    j = pd.DataFrame({'eventId': tr['eventId'].to_numpy(),
                      'px': src[:, 0], 'py': src[:, 1]}).merge(
        p2[['eventId', 'x', 'y', 'x2', 'y2', 'n_clus', 'n_pad', 'q', 'single']],
        on='eventId', how='inner')
    if len(j) < args.min_events:
        print(f'  {sub_run}: only {len(j)} matched tracks - skipped')
        return None
    fit, tf, _ = E.fit_frame(j[['px', 'py']].to_numpy(),
                             j[['x', 'y']].to_numpy(), args.match_win)

    ex, ey = tf.apply(src[:, 0], src[:, 1])
    d_pad, i_pad = tree.query(np.column_stack([ex, ey]))
    fid = d_pad < args.fid_r

    hit = pd.DataFrame({'eventId': tr['eventId'].to_numpy()[fid],
                        'ex': ex[fid], 'ey': ey[fid],
                        'ipad': i_pad[fid]}).merge(
        p2[['eventId', 'x', 'y', 'x2', 'y2', 'a_lead', 'lead_pad']],
        on='eventId', how='left')
    d1 = np.hypot(hit['x'] - hit['ex'], hit['y'] - hit['ey'])
    d2 = np.hypot(hit['x2'] - hit['ex'], hit['y2'] - hit['ey'])
    hit['dmin'] = np.fmin(d1.fillna(np.inf), d2.fillna(np.inf))
    hit['found'] = hit['dmin'] < args.probe_r

    ip = hit['ipad'].to_numpy()
    fnd = hit['found'].to_numpy()
    np.add.at(acc['n'], ip, 1)
    np.add.at(acc['k'], ip[fnd], 1)

    a = hit['a_lead'].to_numpy()
    m = fnd & np.isfinite(a) & (a > 0)
    b = np.clip((a[m] / AMP_BIN).astype(np.int64), 0, AMP_NBIN - 1)
    np.add.at(acc['hist'], (ip[m], b), 1)
    acc['amp_sum'] += float(np.sum(a[m]))
    acc['amp_n'] += int(m.sum())

    # The VMM twin takes the ADC of the hit NEAREST the prediction, not the
    # largest in the event, so `a_lead` is only the same quantity when the pad
    # that led is the pad the track pointed at.  Keep that subset separately:
    # it is the strictly like-for-like histogram, and the count of it against
    # the count above says how often the difference can matter at all.
    same = m & (hit['lead_pad'].to_numpy() == pads_id[ip])
    bs = np.clip((a[same] / AMP_BIN).astype(np.int64), 0, AMP_NBIN - 1)
    np.add.at(acc['hist_own'], (ip[same], bs), 1)

    k, n = int(fnd.sum()), len(hit)
    print(f'  {sub_run}: {n} fiducial tracks, eff {k / max(n, 1):.4f}, '
          f'rot {fit["affine_rotation_deg"]:+.3f} deg, '
          f'{int(m.sum())} amplitudes ({int(same.sum())} on the predicted pad)')
    sys.stdout.flush()
    return dict(sub_run=sub_run, n=n, k=k, eff=k / max(n, 1),
                rot_deg=fit['affine_rotation_deg'],
                det=fit['det'], n_amp=int(m.sum()))


def quantile_from_hist(h, q):
    """Interpolated quantile of a histogram row, in ADC.

    The same estimator the VMM side uses for its DNL-clean median, so the two
    per-pad medians are read off the data the same way -- here it is only a
    convenience (DREAM has no comb), there it is a correction.
    """
    tot = h.sum()
    if tot == 0:
        return np.nan
    cum = np.cumsum(h) / tot
    i = int(np.searchsorted(cum, q, side='left'))
    below = cum[i - 1] if i > 0 else 0.0
    frac = (q - below) / max(cum[i] - below, 1e-12)
    return float((i + frac) * AMP_BIN)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', default='eff_nominal_1')
    ap.add_argument('--data-root', default=E.DATA_ROOT)
    ap.add_argument('--station', default='P2_OUT')
    ap.add_argument('--sub-run', action='append', default=[])
    ap.add_argument('--max-chunks', type=int, default=0)
    ap.add_argument('--track-cut', type=float, default=3.0)
    ap.add_argument('--cluster-r', type=float, default=15.0)
    ap.add_argument('--probe-r', type=float, default=15.0)
    ap.add_argument('--fid-r', type=float, default=9.0)
    ap.add_argument('--match-win', type=float, default=25.0)
    ap.add_argument('--min-events', type=int, default=500)
    ap.add_argument('--no-veto-sparks', action='store_true')
    ap.add_argument('--out', default='out_padadc')
    args = ap.parse_args()

    cfg = sc.RunConfig(key='urw_p2', run=args.run, data_root=args.data_root)
    os.makedirs(args.out, exist_ok=True)
    dets = [d for d in cfg.mappable_detectors() if d.name == args.station]
    if not dets:
        sys.exit(f'{args.station} is not a mappable detector of {args.run}')
    det = dets[0]

    sub_runs = args.sub_run or cfg.find_subruns()
    pads = pad_table(cfg, det)
    npad = len(pads)
    print(f'{args.run} / {det.name}: {npad} mapped pads, '
          f'{len(sub_runs)} sub_run(s)')

    acc = dict(pads=pads, n=np.zeros(npad, np.int64), k=np.zeros(npad, np.int64),
               hist=np.zeros((npad, AMP_NBIN), np.int64),
               hist_own=np.zeros((npad, AMP_NBIN), np.int64),
               amp_sum=0.0, amp_n=0)
    rows = []
    for sub_run in sub_runs:
        r = run_subrun(cfg, det, sub_run, args, acc)
        if r:
            rows.append(r)

    h, ho = acc['hist'], acc['hist_own']
    med = np.array([quantile_from_hist(h[i], 0.5) for i in range(npad)])
    med_own = np.array([quantile_from_hist(ho[i], 0.5) for i in range(npad)])
    p25 = np.array([quantile_from_hist(h[i], 0.25) for i in range(npad)])
    p75 = np.array([quantile_from_hist(h[i], 0.75) for i in range(npad)])
    tab = pd.DataFrame({
        'pad_id': pads['channel_id'].to_numpy(),
        'pad_cx': pads['pad_cx'].to_numpy(), 'pad_cy': pads['pad_cy'].to_numpy(),
        'n_track': acc['n'], 'k_hit': acc['k'],
        'eff': np.where(acc['n'] > 0, acc['k'] / np.maximum(acc['n'], 1), np.nan),
        'amp_n': h.sum(1), 'amp_med': med, 'amp_p25': p25, 'amp_p75': p75,
        'amp_n_own': ho.sum(1), 'amp_med_own': med_own,
        'amp_mean': np.where(h.sum(1) > 0,
                             (h * (np.arange(AMP_NBIN) + 0.5) * AMP_BIN).sum(1)
                             / np.maximum(h.sum(1), 1), np.nan)})
    for c in ('pad_area', 'pad_r', 'pad_row', 'pad_col'):
        if c in pads.columns:
            tab[c] = pads[c].to_numpy()
    base = os.path.join(args.out, f'dream_padadc_{args.run}_{det.name}')
    tab.to_csv(base + '.csv', index=False)
    np.savez_compressed(base + '_hist.npz', hist=h, hist_own=ho,
                        amp_bin=AMP_BIN, pad_id=tab['pad_id'].to_numpy())
    with open(base + '.json', 'w') as fh:
        json.dump(dict(run=args.run, station=det.name, sub_runs=rows,
                       amp_bin=AMP_BIN, amp_nbin=AMP_NBIN,
                       n_pads=npad, args=vars(args)), fh, indent=2, default=str)
    ok = acc['n'] >= 500
    print(f'\n{int(ok.sum())} pads with >=500 tracks; '
          f'total {int(acc["n"].sum())} tracks, '
          f'{int(acc["k"].sum())} hits, '
          f'eff {acc["k"].sum() / max(acc["n"].sum(), 1):.4f}')
    print(f'wrote {base}.csv / _hist.npz / .json')


if __name__ == '__main__':
    main()
