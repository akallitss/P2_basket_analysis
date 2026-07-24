#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
27_pedestal_qa.py

Per-channel pedestal QA for the P2 telescope stations at the SPS, the beam-arm
counterpart of the cosmic bench 09_pedestal_qa.py.

The DREAM zero suppression keeps a pulse only if amplitude >= nsigma *
pedestalRMS_channel, so each channel's effective hit threshold is set by its own
pedestal RMS. A channel with an inflated RMS has a raised threshold (it goes
quiet / dead-looking); a channel with a collapsed RMS fires on noise. This stage
maps the pedestal mean (baseline) and RMS (noise) on the pad plane and checks
whether pedestal structure explains any per-pad hit-rate structure.

Unlike the bench stage, which reads the `pedestals` TTree embedded in the
per-FEU hits_root, this computes the pedestal DIRECTLY from the dedicated
pedestal run's raw samples (the `nt` tree: per-event vectors of channel /
amplitude over the sampling window). Common-noise subtraction is off for this
campaign, so the raw per-channel mean / RMS is the applied pedestal, and no
hits_root fetch is needed -- the pedestal run is small and already local.

Products (<Analysis>/<det_tag>/<run>/pedestals/27_pedestal_qa/):
  pedestal_rms_map_<det>.png     RMS on the pad plane
  pedestal_mean_map_<det>.png    baseline (mean) on the pad plane
  pedestal_rms_dist_<det>.png    RMS distribution + threshold + outliers
  pedestal_rms_vs_rate_<det>.png RMS vs per-pad firing fraction (--with-rate)
  pedestal_qa_<det>.csv          per-channel table + outlier flag

Usage:
  SPS_DATA_ROOT=.../runs SPS_ANALYSIS_ROOT=.../analysis SPS_RUN=<run> \
      python3 27_pedestal_qa.py live [--pedestal-dir DIR] [--nsigma 5]
      [--with-rate] [--rate-sub-run NAME]
"""

import os
import csv
import glob
import argparse

import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import uproot

import sps_config as sc
import p2_mapping as pmap

N_CH_PER_FEU = 512


def find_pedestal_dir(cfg, override=None):
    """Locate the pedestal run's per-FEU .root files. Override wins; else look
    for a `pedestals/<run>/pedestals` tree beside the DATA_ROOT `runs` dir and
    take the newest."""
    if override:
        return override
    parent = os.path.dirname(cfg.DATA_ROOT.rstrip('/'))
    cands = sorted(glob.glob(os.path.join(parent, 'pedestals', '*', 'pedestals')))
    if not cands:
        cands = sorted(glob.glob(os.path.join(parent, 'pedestals', 'pedestals')))
    if not cands:
        raise SystemExit(f'No pedestal run found under {parent}/pedestals/. '
                         'Pass --pedestal-dir.')
    return cands[-1]


def pedestal_for_feu(ped_dir, feu):
    """Per-channel (mean, rms) for one FEU from the pedestal run's `nt` tree,
    accumulated over all samples and events. Returns a DataFrame keyed on
    channel (0..511)."""
    # The per-FEU decoded pedestal run: '..._000_<FEU>.root'. The .root glob
    # already excludes the _ped.prg / _thr.prg RunCtrl outputs.
    files = sorted(glob.glob(os.path.join(ped_dir, f'*_{feu:02d}.root')))
    if not files:
        raise FileNotFoundError(f'No pedestal .root for FEU {feu} in {ped_dir}')
    n = N_CH_PER_FEU
    s = np.zeros(n); s2 = np.zeros(n); cnt = np.zeros(n)
    with uproot.open(files[0]) as f:
        t = f['nt']
        for ch, amp in zip(*[t[b].array(library='np') for b in
                             ('channel', 'amplitude')]):
            ch = np.asarray(ch, np.int64)
            amp = np.asarray(amp, np.float64)
            ok = (ch >= 0) & (ch < n)
            ch, amp = ch[ok], amp[ok]
            s += np.bincount(ch, weights=amp, minlength=n)
            s2 += np.bincount(ch, weights=amp * amp, minlength=n)
            cnt += np.bincount(ch, minlength=n)
    with np.errstate(invalid='ignore', divide='ignore'):
        mean = s / cnt
        var = s2 / cnt - mean ** 2
    rms = np.sqrt(np.clip(var, 0, None))
    return pd.DataFrame({'feu': feu, 'channel': np.arange(n),
                         'ped_mean': mean, 'ped_rms': rms, 'n_samp': cnt})


def load_pedestals(ped_dir, ct):
    """Pedestal (mean, rms) for a station's FEUs, mapped to its connected pads."""
    feus = ct.attrs['feus']
    peds = pd.concat([pedestal_for_feu(ped_dir, f) for f in feus],
                     ignore_index=True)
    m = ct.merge(peds, on=['feu', 'channel'], how='left')
    return m[m['mapped'] & m['ped_rms'].notna() & (m['n_samp'] > 0)].copy()


def load_firing(cfg, ct, sub_run, veto_sparks=True):
    """Per-pad firing fraction in one sub_run, streamed."""
    import p2_io as p2io
    import p2_sparks as ps
    hits_dir = cfg.combined_hits_dir(sub_run)
    br = ['eventId', 'channel', 'feu', 'trigger_timestamp_ns']
    sv = None
    if veto_sparks:
        hv = cfg.hv_monitor_csv(sub_run)
        if os.path.isfile(hv):
            class _S:
                SPARK_CHANNEL = None
                SPARK_IMON_THR = cfg.SPARK_IMON_THR
                SPARK_GUARD_BEFORE = cfg.SPARK_GUARD_BEFORE
                SPARK_GUARD_AFTER = cfg.SPARK_GUARD_AFTER
                BURST_NPADS = cfg.BURST_NPADS
            # station-agnostic here; per-station spark veto handled in 26.
    fire = None
    ev = []
    for a in p2io.iter_hits(hits_dir, br, ct.attrs['feus'], progress=False):
        h = pmap.attach_pads_to_hits(a, ct)
        h = h[h['mapped'] & h['pad_cx'].notna()]
        if not len(h):
            continue
        ev.append(h['eventId'].to_numpy())
        c = h.groupby('channel_id')['eventId'].nunique()
        fire = c if fire is None else fire.add(c, fill_value=0)
    nev = len(np.unique(np.concatenate(ev))) if ev else 0
    return (fire / max(nev, 1)).rename('fire') if fire is not None \
        else pd.Series(dtype=float), nev


def _pad_map(ax, df, val, label, cmap='viridis'):
    s = ax.scatter(df['pad_cx'], df['pad_cy'], c=df[val], s=16, marker='s',
                   cmap=cmap, linewidths=0)
    ax.set_aspect('equal'); ax.set_xlabel('pad_cx [mm]')
    ax.set_ylabel('pad_cy [mm]')
    import matplotlib.pyplot as plt
    plt.colorbar(s, ax=ax, label=label, fraction=0.046, pad=0.04)


def robust_z(x):
    med = np.median(x)
    mad = np.median(np.abs(x - med)) or 1.0
    return 0.6745 * (x - med) / mad


def main():
    import matplotlib.pyplot as plt
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('run_key', nargs='?', default=None)
    ap.add_argument('--pedestal-dir', default=None,
                    help='pedestal run dir with per-FEU .root (default: auto).')
    ap.add_argument('--nsigma', type=float, default=5.0,
                    help='ZS threshold = nsigma * pedestal RMS (default 5).')
    ap.add_argument('--z-thr', type=float, default=5.0,
                    help='robust-z on RMS above which a channel is an outlier.')
    ap.add_argument('--with-rate', action='store_true',
                    help='also load combined hits for the RMS-vs-firing panel.')
    ap.add_argument('--rate-sub-run', default=None,
                    help='sub_run for the firing rate (default: first on disk).')
    args = ap.parse_args()

    cfg = sc.get_config(args.run_key)
    print(cfg)
    ped_dir = find_pedestal_dir(cfg, args.pedestal_dir)
    print(f'  pedestal run: {ped_dir}')
    dets = cfg.mappable_detectors()
    skipped = [d.name for d in cfg.detectors() if d not in dets]
    if skipped:
        print(f'  (no pad map, skipped: {", ".join(skipped)})')

    rate_sub = args.rate_sub_run
    if args.with_rate and rate_sub is None:
        subs = cfg.find_subruns()
        rate_sub = subs[0] if subs else None

    for det in dets:
        print(f'\n== {det.name}  FEUs {det.feus}')
        ct = cfg.channel_table(det)
        try:
            ped = load_pedestals(ped_dir, ct)
        except FileNotFoundError as e:
            print(f'  {e}; skipping'); continue
        thr = args.nsigma * ped['ped_rms']
        ped['threshold_adc'] = thr
        ped['rms_z'] = robust_z(ped['ped_rms'].to_numpy())
        ped['outlier'] = np.abs(ped['rms_z']) > args.z_thr
        med_rms = float(ped['ped_rms'].median())
        print(f'  {len(ped)} mapped channels  median RMS {med_rms:.1f} ADC  '
              f'-> threshold {args.nsigma:g}sigma = {args.nsigma*med_rms:.0f} ADC'
              f'  outliers {int(ped["outlier"].sum())}')

        out = cfg.out_dir(det.det_tag, 'pedestals', '27_pedestal_qa')

        has_geom = pmap.has_tile_geometry(ct)
        if has_geom:
            fig, ax = plt.subplots(figsize=(7.4, 6))
            _pad_map(ax, ped, 'ped_rms', 'pedestal RMS [ADC]')
            ax.set_title(f'{det.name} pedestal RMS (noise) on pads')
            fig.tight_layout()
            fig.savefig(os.path.join(out, f'pedestal_rms_map_{det.det_tag}.png'),
                        dpi=150, bbox_inches='tight'); plt.close(fig)

            fig, ax = plt.subplots(figsize=(7.4, 6))
            _pad_map(ax, ped, 'ped_mean', 'pedestal mean [ADC]', cmap='cividis')
            ax.set_title(f'{det.name} pedestal baseline (mean) on pads')
            fig.tight_layout()
            fig.savefig(os.path.join(out, f'pedestal_mean_map_{det.det_tag}.png'),
                        dpi=150, bbox_inches='tight'); plt.close(fig)

        fig, ax = plt.subplots(figsize=(8.5, 5))
        ax.hist(ped['ped_rms'], bins=60, color='steelblue')
        ax.axvline(med_rms, color='k', ls='--',
                   label=f'median {med_rms:.1f} ADC')
        if ped['outlier'].any():
            for r in ped.loc[ped['outlier'], 'ped_rms']:
                ax.axvline(r, color='crimson', lw=0.6, alpha=0.6)
        ax.set_xlabel('pedestal RMS [ADC]'); ax.set_ylabel('channels')
        ax.set_title(f'{det.name} pedestal RMS distribution '
                     f'({int(ped["outlier"].sum())} outliers)')
        ax.legend(); ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(os.path.join(out, f'pedestal_rms_dist_{det.det_tag}.png'),
                    dpi=150, bbox_inches='tight'); plt.close(fig)

        if args.with_rate and rate_sub:
            fire, nev = load_firing(cfg, ct, rate_sub)
            ped = ped.merge(fire.rename('fire'), left_on='channel_id',
                            right_index=True, how='left')
            ped['fire'] = ped['fire'].fillna(0.0)
            fig, ax = plt.subplots(figsize=(8.5, 5.5))
            s = ax.scatter(ped['ped_rms'], ped['fire'], s=14,
                           c=ped['outlier'].map({True: 'crimson',
                                                 False: 'steelblue'}),
                           linewidths=0, alpha=0.6)
            ax.set_xlabel('pedestal RMS [ADC]')
            ax.set_ylabel(f'firing fraction ({rate_sub}, {nev} events)')
            ax.set_title(f'{det.name} pedestal RMS vs firing rate '
                         f'(red = RMS outlier)')
            ax.grid(True, alpha=0.3)
            fig.tight_layout()
            fig.savefig(os.path.join(out,
                                     f'pedestal_rms_vs_rate_{det.det_tag}.png'),
                        dpi=150, bbox_inches='tight'); plt.close(fig)

        cols = [c for c in ['channel_id', 'feu', 'channel', 'pad_cx', 'pad_cy',
                            'ped_mean', 'ped_rms', 'threshold_adc', 'rms_z',
                            'outlier', 'fire'] if c in ped.columns]
        ped[cols].to_csv(os.path.join(out, f'pedestal_qa_{det.det_tag}.csv'),
                         index=False)
        print(f'  -> {out}')


if __name__ == '__main__':
    main()
