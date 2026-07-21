#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
04_m3_reference_qa.py

M3 reference-tracker QA, independent of the P2 detector under test. Mirrors
nTof_x17/mx_june_cosmic_qa/02_m3_reference_qa.py (same M3RefTracking reader).

Products (written to <Analysis>/<detN>/<run>/<sub_run>/04_m3_reference_qa/):
  m3_chi2_distributions.png  Chi2X / Chi2Y (raw, pre single-track cut)
  m3_track_multiplicity.png  good tracks per event
  m3_angles.png              theta_x / theta_y of clean single tracks
  m3_beam_profile_detz.png   projected (x,y) at the P2 plane z (run_config)
  m3_up_down_positions.png    hit maps at the up + down tracker stations

Usage: python3 04_m3_reference_qa.py [run_key] [--chi2-cut C]

With a non-default --chi2-cut every product gets a '_chi2' filename tag
(p2_qa_config.chi2_tag) so cut-variant plots never overwrite the standard ones.
"""

import os
import sys
import glob
import json
import argparse
import datetime
import matplotlib
matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import uproot

import p2_qa_config as qa

CHI2_BIN_MIN = 10.0   # chi2-acceptance-vs-time bin width [min]

_M3_PKG = os.path.expanduser(
    '~/Documents/PostDocSaclay/nTof_x17/cosmic_bench_analysis')
if _M3_PKG not in sys.path:
    sys.path.insert(0, _M3_PKG)
import awkward as ak  # noqa: E402
from M3RefTracking import M3RefTracking, get_ray_data, get_xy_angles  # noqa: E402

# Tracking-v2 recommended recipe: Chi2X,Chi2Y < chi2_cut AND NClusX,NClusY >= MIN_NCLUS.
# The old naive chi2<20 let 2-point-per-coordinate fits (exact -> denormal-tiny chi2)
# slip through; MIN_NCLUS drops them. Single source of truth in p2_qa_config;
# the chi2 cut can be overridden per invocation with --chi2-cut.
MIN_NCLUS = qa.M3_MIN_NCLUS


def _det_plane_z(cfg):
    with open(cfg.run_config_path) as f:
        det = json.load(f)['detectors']
    for d in det:
        if d.get('name') == cfg.DET_NAME:
            return float(d['det_center_coords']['z'])
    return 232.0


def _sym_range(*arrays, pad=1.05):
    """Symmetric [-L, L] range covering the 0.5-99.5 percentiles of the data."""
    v = np.concatenate([np.abs(a[np.isfinite(a)]) for a in arrays])
    L = np.percentile(v, 99.5) * pad
    return [[-L, L], [-L, L]]


def plot_chi2(cfg, m3_dir, out_dir, chi2_cut, sfx):
    raw = get_ray_data(m3_dir)
    chi2x = ak.to_numpy(ak.flatten(raw['Chi2X']))
    chi2y = ak.to_numpy(ak.flatten(raw['Chi2Y']))
    n_tracks = ak.to_numpy(ak.num(raw['Chi2X'], axis=1))

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    for ax, data, lbl in [(axes[0], chi2x, 'Chi2X'), (axes[1], chi2y, 'Chi2Y')]:
        finite = data[np.isfinite(data)]
        ax.hist(finite[finite < 50], bins=100, color='indianred', edgecolor='none')
        ax.axvline(chi2_cut, color='k', ls='--', lw=1, label=f'cut = {chi2_cut:g}')
        ax.set_xlabel(lbl); ax.set_ylabel('Tracks'); ax.set_yscale('log')
        ax.legend(); ax.grid(True, alpha=0.3)
    fig.suptitle(f'M3 raw track chi2 — {cfg.RUN}')
    fig.tight_layout()
    fig.savefig(f'{out_dir}/m3_chi2_distributions{sfx}.png', dpi=150,
                bbox_inches='tight')
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(n_tracks, bins=range(0, int(n_tracks.max()) + 2), align='left',
            color='steelblue', edgecolor='white')
    ax.set_xlabel('Reconstructed tracks per event'); ax.set_ylabel('Events')
    ax.set_yscale('log')
    ax.set_title(f'M3 track multiplicity — {cfg.RUN}\n'
                 f'{len(n_tracks):,} events, mean {n_tracks.mean():.2f} tracks/event')
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(f'{out_dir}/m3_track_multiplicity{sfx}.png', dpi=150,
                bbox_inches='tight')
    plt.close(fig)
    return len(n_tracks)


def plot_angles_and_positions(cfg, m3_dir, out_dir, det_z, chi2_cut, sfx):
    rays = M3RefTracking(m3_dir, chi2_cut=chi2_cut, min_nclus=MIN_NCLUS)
    n_clean = len(ak.to_numpy(rays.ray_data['X_Up']))

    x_ang, y_ang, _ = get_xy_angles(rays.ray_data)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    for ax, data, lbl in [(axes[0], np.degrees(x_ang), r'$\theta_x$ [deg]'),
                          (axes[1], np.degrees(y_ang), r'$\theta_y$ [deg]')]:
        ax.hist(data, bins=120, range=(-40, 40), color='seagreen', edgecolor='none')
        ax.set_xlabel(lbl); ax.set_ylabel('Tracks')
        ax.set_title(f'{lbl}  (median {np.median(data):.2f}, std {np.std(data):.2f})')
        ax.grid(True, alpha=0.3)
    fig.suptitle(f'M3 single-track angles — {cfg.RUN}  ({n_clean:,} tracks)')
    fig.tight_layout()
    fig.savefig(f'{out_dir}/m3_angles{sfx}.png', dpi=150, bbox_inches='tight')
    plt.close(fig)

    xs, ys, _ = rays.get_xy_positions(det_z)
    xs, ys = np.asarray(xs), np.asarray(ys)
    rng = _sym_range(xs, ys)
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    h = ax.hist2d(xs, ys, bins=80, range=rng, cmap='viridis')
    fig.colorbar(h[3], ax=ax, label='Tracks')
    ax.set_xlabel('X [mm]'); ax.set_ylabel('Y [mm]'); ax.set_aspect('equal')
    ax.set_title(f'M3 track projection at P2 plane z = {det_z:.0f} mm\n{cfg.RUN}')
    fig.tight_layout()
    fig.savefig(f'{out_dir}/m3_beam_profile_detz{sfx}.png', dpi=150,
                bbox_inches='tight')
    plt.close(fig)

    rd = rays.ray_data
    xu, yu = ak.to_numpy(rd['X_Up']), ak.to_numpy(rd['Y_Up'])
    xd, yd = ak.to_numpy(rd['X_Down']), ak.to_numpy(rd['Y_Down'])
    zu = float(np.mean(ak.to_numpy(rd['Z_Up'])))
    zd = float(np.mean(ak.to_numpy(rd['Z_Down'])))
    rng2 = _sym_range(xu, yu, xd, yd)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    for ax, xx, yy, z in [(axes[0], xu, yu, zu), (axes[1], xd, yd, zd)]:
        h = ax.hist2d(xx, yy, bins=80, range=rng2, cmap='viridis')
        fig.colorbar(h[3], ax=ax, label='Tracks')
        ax.set_xlabel('X [mm]'); ax.set_ylabel('Y [mm]'); ax.set_aspect('equal')
        ax.set_title(f'Station at z = {z:.0f} mm')
    fig.suptitle(f'M3 station hit maps — {cfg.RUN}')
    fig.tight_layout()
    fig.savefig(f'{out_dir}/m3_up_down_positions{sfx}.png', dpi=150,
                bbox_inches='tight')
    plt.close(fig)
    return n_clean


def plot_chi2_acceptance_vs_time(cfg, out_dir, chi2_cut, sfx, summary):
    """Fraction of M3 single tracks passing the Chi2X,Chi2Y < chi2_cut cut, in
    time bins across the run. The M3 telescope is the efficiency REFERENCE, so a
    transient collapse of its track quality (chi2 tail blowing up) silently
    corrupts the efficiency denominator without touching the P2 detector -- it
    reads as a fake efficiency dip. This trace exposes such reference glitches
    and time-stamps them; a glitch is flagged when the acceptance drops below
    50% of its own run median while the raw single-track rate stays up.

    Reads evttime + rayN + Chi2X/Chi2Y straight from the m3 tracking files
    (Chi2X/Chi2Y are per-track std::vector<double>; single tracks -> length 1)."""
    files = sorted(glob.glob(os.path.join(cfg.m3_tracking_dir, '*.root')))
    if not files:
        summary.append('  chi2 acceptance vs time: no m3 tracking files')
        return
    ev, c2x, c2y, ok1 = [], [], [], []
    for f in files:
        with uproot.open(f) as fh:
            a = fh[fh.keys()[0]].arrays(['evttime', 'rayN', 'Chi2X', 'Chi2Y'],
                                        library='np')
        one = a['rayN'] == 1                       # single-track events
        if not one.any():
            continue
        ev.append(a['evttime'][one].astype(np.float64))
        c2x.append(np.array([x[0] if len(x) else np.nan for x in a['Chi2X'][one]]))
        c2y.append(np.array([y[0] if len(y) else np.nan for y in a['Chi2Y'][one]]))
    if not ev:
        summary.append('  chi2 acceptance vs time: no single tracks')
        return
    ev = np.concatenate(ev)
    c2x = np.concatenate(c2x); c2y = np.concatenate(c2y)
    good = (c2x < chi2_cut) & (c2y < chi2_cut)
    t_h = (ev * 10 - ev.min() * 10) / 1e9 / 3600.0   # evttime = trig_ns/10
    bw = CHI2_BIN_MIN / 60.0
    edges = np.arange(0.0, t_h.max() + bw, bw)
    n_all = np.histogram(t_h, bins=edges)[0]
    n_good = np.histogram(t_h[good], bins=edges)[0]
    enough = n_all >= 50
    acc = np.divide(n_good, n_all, out=np.full(len(n_all), np.nan),
                    where=enough)
    ctr_h = edges[:-1] + bw / 2

    try:
        t0 = datetime.datetime.fromisoformat(
            json.load(open(cfg.run_config_path))['start_time'])
        x = np.array([t0 + datetime.timedelta(hours=float(h)) for h in ctr_h])
        use_wall = True
    except Exception:
        x, use_wall = ctr_h, False

    base = np.nanmedian(acc)                          # run-median acceptance
    glitch = enough & (acc < 0.5 * base)
    fig, ax = plt.subplots(figsize=(11, 4.2))
    ax.plot(x, 100 * acc, 'o-', color='crimson', ms=4, lw=1.4)
    ax.axhline(100 * base, ls=':', color='grey',
               label=f'run median {100*base:.0f}%')
    win = None
    if glitch.any():
        gx = x[glitch]
        lo, hi = gx.min(), gx.max()
        pad = (datetime.timedelta(minutes=CHI2_BIN_MIN) if use_wall
               else CHI2_BIN_MIN / 60.0)
        ax.axvspan(lo - pad, hi + pad, color='gold', alpha=0.25,
                   label='reference glitch')
        win = (lo, hi)
    if use_wall:
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        ax.set_xlabel(f'wall clock (run start {t0:%Y-%m-%d %H:%M})')
    else:
        ax.set_xlabel('time since run start [h]')
    ax.set_ylim(0, 100); ax.grid(True, alpha=0.3); ax.legend(fontsize=8)
    ttl = (f'{cfg.DET_NAME} M3 reference: fraction of tracks passing '
           f'chi2<{chi2_cut:g} vs time — {cfg.RUN}')
    if win is not None:
        wl = (f'{win[0]:%m-%d %H:%M}-{win[1]:%H:%M}' if use_wall
              else f'{win[0]:.1f}-{win[1]:.1f} h')
        ttl += f'\n⚠ REFERENCE GLITCH: acceptance collapsed at {wl}'
    ax.set_title(ttl, fontsize=10)
    fig.tight_layout()
    fig.savefig(f'{out_dir}/m3_chi2_acceptance_vs_time{sfx}.png', dpi=150,
                bbox_inches='tight')
    plt.close(fig)

    if win is not None:
        wl = (f'{win[0]:%m-%d %H:%M}-{win[1]:%H:%M}' if use_wall
              else f'{win[0]:.1f}-{win[1]:.1f} h')
        summary.append(f'  ⚠ M3 REFERENCE GLITCH: chi2<{chi2_cut:g} acceptance '
                       f'collapsed to <{50}% of the {100*base:.0f}% median at '
                       f'~{wl} — efficiency in that window is unreliable '
                       f'(reference, not the P2 detector).')
    else:
        summary.append(f'  M3 chi2<{chi2_cut:g} acceptance steady at '
                       f'~{100*base:.0f}% (no reference glitch).')


def main():
    ap = argparse.ArgumentParser(description='M3 reference-tracker QA.')
    ap.add_argument('run_key', nargs='?', default=qa.DEFAULT_RUN)
    ap.add_argument('--chi2-cut', type=float, default=qa.M3_CHI2_CUT)
    args = ap.parse_args()

    cfg = qa.get_config(args.run_key)
    print(cfg)
    sfx = qa.chi2_tag(args.chi2_cut)
    out_dir = cfg.out_dir('04_m3_reference_qa')
    m3_dir = os.path.join(cfg.m3_tracking_dir, '')  # trailing sep for the reader
    det_z = _det_plane_z(cfg)
    print(f'P2 plane z = {det_z:.1f} mm | chi2 cut = {args.chi2_cut:g}')
    n_events = plot_chi2(cfg, m3_dir, out_dir, args.chi2_cut, sfx)
    n_clean = plot_angles_and_positions(cfg, m3_dir, out_dir, det_z,
                                        args.chi2_cut, sfx)
    summary = [f'M3 reference QA — {cfg.DET_NAME}  {cfg.RUN}/{cfg.SUB_RUN}',
               f'  M3 events: {n_events:,} | clean single tracks: {n_clean:,} '
               f'({100 * n_clean / max(n_events, 1):.1f}%)']
    plot_chi2_acceptance_vs_time(cfg, out_dir, args.chi2_cut, sfx, summary)
    txt = '\n'.join(summary)
    print(txt)
    with open(f'{out_dir}/m3_reference_summary{sfx}.txt', 'w') as f:
        f.write(txt + '\n')
    print(f'M3 reference QA written to: {out_dir}')


if __name__ == '__main__':
    main()
