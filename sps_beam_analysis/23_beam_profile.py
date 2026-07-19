#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
23_beam_profile.py

Beam-spot and time-structure characterisation per station / sub_run at the SPS.

Products per point (sub_run), for each station:
  1. Beam-spot hit map (per-pad counts on the real pad tiles) + a 2D histogram
     of the per-event cluster centroids.
  2. x and y projections of the cluster centroids with Gaussian CORE fits
     (iteratively fit the central +-2 sigma so the beam halo does not bias the
     core width).
  3. Rate vs wall-clock time within the sub_run -- the SPS slow-extraction spill
     structure appears as rate modulation. With --beam-csv the DAQ beam_monitor
     intensity is overlaid on the same time axis.
  4. Pile-up indicator: mean pads/event vs instantaneous trigger rate (rate in a
     sliding time window around each event).

--beam-csv PATH overlays the per-day DAQ beam-monitor CSV
(columns: timestamp, unix_ts, intensity_e10). The run start wall-clock is read
from the chunk-0 filename (…_YYMMDD_HHhMM_…) so event time [s] maps to
wall-clock; without the CSV the rate panel still shows the DAQ trigger rate.

Products (<Analysis>/<det_tag>/<run>/<sub_run>/23_beam_profile/):
  beam_spot_<pt><suffix>.png       hit map + centroid 2D histogram
  profiles_<pt><suffix>.png        x/y projections + Gaussian core fits
  timing_<pt><suffix>.png          rate vs wall-clock (+ beam CSV), pile-up
  beam_profile_<pt><suffix>.csv    core centres/sigmas, rate, pile-up slope

Usage:
  python3 23_beam_profile.py [run_key] [--det P2_OUT] [--sub-run NAME]
        [--beam-csv PATH] [--no-veto-sparks]
"""

import os
import re
import argparse

import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import sps_config as sc
import p2_mapping as pmap
import p2_io as p2io
import p2_sparks as ps
import sps_cluster as scl


def _spark_veto(cfg, det, hv_csv):
    class _Shim:
        SPARK_CHANNEL = det.spark_channel
        SPARK_IMON_THR = cfg.SPARK_IMON_THR
        SPARK_GUARD_BEFORE = cfg.SPARK_GUARD_BEFORE
        SPARK_GUARD_AFTER = cfg.SPARK_GUARD_AFTER
        BURST_NPADS = cfg.BURST_NPADS
    return ps.SparkVeto.from_csv(hv_csv, _Shim)


def pad_hit_counts(hits_dir, ct, min_amp, t_min, veto):
    counts = None
    for df in p2io.iter_hits(hits_dir, ['eventId', 'channel', 'amplitude',
                                        'feu', 'trigger_timestamp_ns'],
                             ct.attrs['feus'], min_amp=min_amp, progress=False):
        if t_min > 0:
            df = df[df['trigger_timestamp_ns'].astype(np.int64) / 1e9 >= t_min]
        if veto is not None and len(df):
            df, _ = veto.apply(df)
        if not len(df):
            continue
        h = pmap.attach_pads_to_hits(df, ct)
        h = h[h['mapped'] & h['pad_cx'].notna()]
        if not len(h):
            continue
        c = h.groupby('channel_id').size()
        counts = c if counts is None else counts.add(c, fill_value=0)
    return counts.astype(np.int64) if counts is not None else pd.Series(dtype=np.int64)


def _gauss(x, a, mu, sig):
    return a * np.exp(-0.5 * ((x - mu) / sig) ** 2)


def fit_core(vals, nbins=80, n_iter=3):
    """Gaussian CORE fit: fit, then refit within +-2 sigma. Returns
    dict(mu, sigma, ok, edges, counts)."""
    vals = np.asarray(vals, float)
    vals = vals[np.isfinite(vals)]
    counts, edges = np.histogram(vals, bins=nbins)
    ctr = 0.5 * (edges[:-1] + edges[1:])
    out = dict(mu=np.nan, sigma=np.nan, ok=False, edges=edges, counts=counts)
    if len(vals) < 30:
        return out
    mu, sig = float(np.median(vals)), float(vals.std())
    for _ in range(n_iter):
        win = np.abs(ctr - mu) < 2 * sig
        if win.sum() < 5:
            break
        try:
            p, _ = curve_fit(_gauss, ctr[win], counts[win],
                             p0=[counts[win].max(), mu, sig], maxfev=8000)
        except (RuntimeError, ValueError):
            return out
        mu, sig = float(p[1]), abs(float(p[2]))
    out.update(mu=mu, sigma=sig, ok=True)
    return out


def run_start_wallclock(chunk0):
    m = re.search(r'(\d{6})_(\d{2})H(\d{2})', os.path.basename(chunk0))
    if not m:
        return None
    return pd.Timestamp(f'20{m.group(1)[:2]}-{m.group(1)[2:4]}-'
                        f'{m.group(1)[4:]} {m.group(2)}:{m.group(3)}:00')


def load_beam_csv(path):
    if not path or not os.path.isfile(path):
        return None
    df = pd.read_csv(path)
    if 'unix_ts' not in df.columns or 'intensity_e10' not in df.columns:
        print(f'    [beam-csv] {path} missing unix_ts/intensity_e10; ignoring')
        return None
    df['wall'] = pd.to_datetime(df['unix_ts'], unit='s')
    return df


def plot_spot(pad_hits, ev, ct, det, lbl, sub, out_png):
    from matplotlib.collections import PolyCollection
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.6))
    ax = axes[0]
    if pmap.has_tile_geometry(ct):
        pads, verts = pmap.pad_tiles(ct)
        counts = (pad_hits.reindex(pads['channel_id']).fillna(0).to_numpy(float))
        fired = counts > 0
        ax.add_collection(PolyCollection(verts[~fired], facecolors='0.92',
                                         edgecolors='0.7', linewidths=0.3))
        if fired.any():
            pc = PolyCollection(verts[fired], array=counts[fired],
                                cmap='inferno',
                                norm=matplotlib.colors.LogNorm(vmin=1),
                                edgecolors='face', linewidths=0.2)
            ax.add_collection(pc)
            fig.colorbar(pc, ax=ax, label='pad hits (log)')
        ax.autoscale_view()
    ax.set_aspect('equal'); ax.set_xlabel('pad_cx [mm]')
    ax.set_ylabel('pad_cy [mm]'); ax.set_title('per-pad hit map')

    ax = axes[1]
    if len(ev):
        hb = ax.hexbin(ev['x'], ev['y'], gridsize=50, cmap='viridis', mincnt=1)
        fig.colorbar(hb, ax=ax, label='events / bin')
    ax.set_aspect('equal'); ax.set_xlabel('cluster x [mm]')
    ax.set_ylabel('cluster y [mm]'); ax.set_title('cluster-centroid density')
    fig.suptitle(f'{det.name} beam spot — {lbl} ({sub})', fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_png, dpi=150, bbox_inches='tight'); plt.close(fig)


def plot_profiles(ev, fx, fy, det, lbl, sub, out_png):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for ax, f, col, name in [(axes[0], fx, 'x', 'x'), (axes[1], fy, 'y', 'y')]:
        ctr = 0.5 * (f['edges'][:-1] + f['edges'][1:])
        ax.step(ctr, f['counts'], where='mid', color='steelblue', lw=1.3,
                label='cluster centroids')
        if f['ok']:
            xs = np.linspace(f['mu'] - 4 * f['sigma'], f['mu'] + 4 * f['sigma'],
                             300)
            amp = f['counts'][np.argmin(np.abs(ctr - f['mu']))]
            ax.plot(xs, _gauss(xs, amp, f['mu'], f['sigma']), '-',
                    color='crimson', lw=2,
                    label=f'core {name}0={f["mu"]:.1f}, '
                          rf'$\sigma$={f["sigma"]:.1f} mm')
        ax.set_xlabel(f'cluster {name} [mm]'); ax.set_ylabel('events / bin')
        ax.grid(True, alpha=0.3); ax.legend(fontsize=9)
    fig.suptitle(f'{det.name} beam profile projections — {lbl} ({sub})',
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_png, dpi=150, bbox_inches='tight'); plt.close(fig)


def plot_timing(ev, det, lbl, sub, t0_wall, beam_df, out_png):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    t = ev['t'].to_numpy()
    t = t - t.min()
    # rate vs time
    ax = axes[0]
    dt = 0.5
    nb = max(4, int((t.max() - t.min()) / dt) + 1)
    counts, edges = np.histogram(t, bins=nb)
    ctr = 0.5 * (edges[:-1] + edges[1:])
    rate = counts / (edges[1] - edges[0])
    ax.plot(ctr, rate, '-', color='steelblue', lw=1.3, label='DAQ trigger rate')
    ax.set_xlabel('time since first event [s]'); ax.set_ylabel('rate [Hz]',
                                                              color='steelblue')
    ax.grid(True, alpha=0.3)
    if beam_df is not None and t0_wall is not None:
        w0 = (t0_wall + pd.to_timedelta(ev['t'].min(), unit='s'))
        sel = ((beam_df['wall'] >= w0 - pd.Timedelta(seconds=2)) &
               (beam_df['wall'] <= w0 + pd.to_timedelta(t.max() + 2, unit='s')))
        b = beam_df[sel]
        if len(b):
            ax2 = ax.twinx()
            bt = (b['wall'] - w0).dt.total_seconds()
            ax2.plot(bt, b['intensity_e10'], 's--', color='darkorange', ms=3,
                     alpha=0.8, label='beam intensity')
            ax2.set_ylabel(r'intensity [$10^{10}$]', color='darkorange')
    ax.set_title('rate vs time (spill structure)')

    # pile-up: pads/event vs instantaneous rate
    ax = axes[1]
    if len(ev) > 20:
        order = np.argsort(ev['t'].to_numpy())
        ts = ev['t'].to_numpy()[order]
        npad = ev['n_pad'].to_numpy()[order]
        win = 0.2  # s
        lo = np.searchsorted(ts, ts - win)
        hi = np.searchsorted(ts, ts + win)
        inst_rate = (hi - lo) / (2 * win)
        # profile: mean npad in rate bins
        rb = np.linspace(inst_rate.min(), np.percentile(inst_rate, 99), 15)
        idx = np.clip(np.digitize(inst_rate, rb), 1, len(rb) - 1)
        mean_np = [npad[idx == k].mean() if (idx == k).any() else np.nan
                   for k in range(1, len(rb))]
        ax.plot(0.5 * (rb[:-1] + rb[1:]), mean_np, 'o-', color='purple', lw=1.5)
        ax.set_xlabel('instantaneous rate [Hz]')
        ax.set_ylabel('mean pads / event')
        ax.grid(True, alpha=0.3)
    ax.set_title('pile-up indicator')
    fig.suptitle(f'{det.name} beam timing — {lbl} ({sub})', fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_png, dpi=150, bbox_inches='tight'); plt.close(fig)


def process_det(cfg, det, args, subruns, beam_df):
    ct = cfg.channel_table(det)
    print(f'  {det.name}: FEUs {det.feus}')
    rows = []
    for sub in subruns:
        hv = cfg.subrun_mesh_hv(sub, det)
        lbl = f'{hv}V' if hv is not None else sub
        out_dir = cfg.out_dir(det.det_tag, sub, '23_beam_profile')
        hits_dir = cfg.combined_hits_dir(sub)
        chunk0 = p2io.hit_files(hits_dir)[0]
        hv_csv = cfg.hv_monitor_csv(sub)
        t_min = scl.settle_t_min(hv_csv, det.spark_channel, chunk0)
        veto = (_spark_veto(cfg, det, hv_csv)
                if args.veto_sparks and os.path.isfile(hv_csv)
                and det.spark_channel else None)
        ev = scl.stream_event_clusters(hits_dir, ct, args.cluster_r,
                                       min_amp=cfg.MIN_AMP, veto=veto,
                                       t_min=t_min)
        if not len(ev):
            print(f'    {lbl}: no events')
            continue
        pad_hits = pad_hit_counts(hits_dir, ct, cfg.MIN_AMP, t_min, veto)
        fx = fit_core(ev['x']); fy = fit_core(ev['y'])
        t0 = run_start_wallclock(chunk0)
        plot_spot(pad_hits, ev, ct, det, lbl, sub,
                  os.path.join(out_dir, f'beam_spot_{lbl}{args.suffix}.png'))
        plot_profiles(ev, fx, fy, det, lbl, sub,
                      os.path.join(out_dir, f'profiles_{lbl}{args.suffix}.png'))
        plot_timing(ev, det, lbl, sub, t0, beam_df,
                    os.path.join(out_dir, f'timing_{lbl}{args.suffix}.png'))
        live = float(ev['t'].max() - ev['t'].min())
        row = dict(sub_run=sub, hv=hv, n_events=len(ev),
                   rate_hz=len(ev) / live if live > 0 else np.nan,
                   x0_mm=fx['mu'], sigma_x_mm=fx['sigma'],
                   y0_mm=fy['mu'], sigma_y_mm=fy['sigma'],
                   mean_pads=float(ev['n_pad'].mean()))
        rows.append(row)
        pd.DataFrame([row]).to_csv(
            os.path.join(out_dir, f'beam_profile_{lbl}{args.suffix}.csv'),
            index=False)
        print(f'    {lbl}: {len(ev):,} ev, spot sigma_x={fx["sigma"]:.1f} '
              f'sigma_y={fy["sigma"]:.1f} mm, mean pads/ev {row["mean_pads"]:.1f}')
    return rows


def process_det_chanspace(cfg, det, args, subruns, beam_df):
    """Channel-space beam profiling when no pad geometry is available: per-(feu,
    channel) occupancy, per-FEU hit fractions, hits/event, rate vs wall-clock
    (spill structure) and a pile-up indicator. No x/y beam spot (needs the pad
    map). Products under <det_tag>/<run>/<sub>/23_beam_profile_chanspace/."""
    feus = det.feus
    print(f'  {det.name}: CHANNEL-SPACE mode (no pad map), FEUs {feus}')
    for sub in subruns:
        out_dir = cfg.out_dir(det.det_tag, sub, '23_beam_profile_chanspace')
        hits_dir = cfg.combined_hits_dir(sub)
        try:
            chunk0 = p2io.hit_files(hits_dir)[0]
            counts, ev = scl.stream_chan_occupancy(hits_dir, feus,
                                                   min_amp=cfg.MIN_AMP,
                                                   progress=True)
        except Exception as e:  # noqa: BLE001
            print(f'    {sub}: unreadable (incomplete transfer?), skipping — {e}')
            continue
        if not len(ev):
            print(f'    {sub}: no events')
            continue
        t0 = run_start_wallclock(chunk0)
        # occupancy grid: FEU rows x channel cols
        occ = np.zeros((len(feus), 512), dtype=float)
        feu_idx = {f: i for i, f in enumerate(feus)}
        for (f, ch), n in counts.items():
            if f in feu_idx and 0 <= ch < 512:
                occ[feu_idx[f], int(ch)] = n
        per_feu = occ.sum(axis=1)
        pd.DataFrame(occ, index=[f'FEU{f}' for f in feus]).to_csv(
            os.path.join(out_dir, f'occupancy_{sub}.csv'))

        fig, axes = plt.subplots(2, 2, figsize=(15, 9))
        # (0,0) occupancy heatmap
        ax = axes[0, 0]
        im = ax.imshow(occ, aspect='auto', origin='lower', cmap='inferno',
                       norm=matplotlib.colors.LogNorm(
                           vmin=max(1, occ[occ > 0].min()) if (occ > 0).any()
                           else 1),
                       extent=[0, 512, -0.5, len(feus) - 0.5])
        fig.colorbar(im, ax=ax, label='hits (log)')
        ax.set_yticks(range(len(feus)))
        ax.set_yticklabels([f'FEU{f}' for f in feus])
        ax.set_xlabel('channel (0-511)'); ax.set_title('per-(FEU, channel) '
                                                       'occupancy')
        # (0,1) per-FEU hit fraction
        ax = axes[0, 1]
        frac = 100 * per_feu / per_feu.sum()
        ax.bar([f'FEU{f}' for f in feus], frac, color='steelblue')
        for i, v in enumerate(frac):
            ax.text(i, v, f'{v:.1f}%', ha='center', va='bottom', fontsize=9)
        ax.set_ylabel('hit fraction [%]'); ax.set_title('beam illumination '
                                                       'across FEUs')
        ax.grid(True, axis='y', alpha=0.3)
        # (1,0) rate vs time
        ax = axes[1, 0]
        t = ev['t'].to_numpy(); t = t - t.min()
        nb = max(4, int(t.max() / 0.5) + 1)
        h, e = np.histogram(t, bins=nb)
        cr = 0.5 * (e[:-1] + e[1:])
        ax.plot(cr, h / (e[1] - e[0]), '-', color='steelblue', lw=1.2,
                label='DAQ trigger rate')
        ax.set_xlabel('time since first event [s]')
        ax.set_ylabel('rate [Hz]', color='steelblue'); ax.grid(True, alpha=0.3)
        if beam_df is not None and t0 is not None:
            w0 = t0 + pd.to_timedelta(ev['t'].min(), unit='s')
            sel = ((beam_df['wall'] >= w0 - pd.Timedelta(seconds=2)) &
                   (beam_df['wall'] <= w0 + pd.to_timedelta(t.max() + 2,
                                                            unit='s')))
            b = beam_df[sel]
            if len(b):
                ax2 = ax.twinx()
                ax2.plot((b['wall'] - w0).dt.total_seconds(),
                         b['intensity_e10'], 's--', color='darkorange', ms=3,
                         alpha=0.8)
                ax2.set_ylabel(r'intensity [$10^{10}$]', color='darkorange')
        ax.set_title('rate vs time (spill structure)')
        # (1,1) pile-up: mean hits/event vs instantaneous rate
        ax = axes[1, 1]
        order = np.argsort(ev['t'].to_numpy())
        ts = ev['t'].to_numpy()[order]; nh = ev['n_hit'].to_numpy()[order]
        win = 0.2
        lo = np.searchsorted(ts, ts - win); hi = np.searchsorted(ts, ts + win)
        inst = (hi - lo) / (2 * win)
        rb = np.linspace(inst.min(), np.percentile(inst, 99), 15)
        idx = np.clip(np.digitize(inst, rb), 1, len(rb) - 1)
        mean_nh = [nh[idx == k].mean() if (idx == k).any() else np.nan
                   for k in range(1, len(rb))]
        ax.plot(0.5 * (rb[:-1] + rb[1:]), mean_nh, 'o-', color='purple', lw=1.5)
        ax.set_xlabel('instantaneous rate [Hz]')
        ax.set_ylabel('mean hits / event'); ax.grid(True, alpha=0.3)
        ax.set_title('pile-up indicator')

        fig.suptitle(f'{det.name} beam profile (channel space) — {sub}   '
                     f'[{len(ev):,} events]', fontsize=13)
        fig.tight_layout(rect=(0, 0, 1, 0.97))
        fig.savefig(os.path.join(out_dir, f'profile_chanspace_{sub}.png'),
                    dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f'    {sub}: {len(ev):,} ev, per-FEU % = '
              f'{dict(zip(["FEU%d" % f for f in feus], frac.round(1)))}, '
              f'median hits/ev {ev["n_hit"].median():.0f}')
        print(f'    -> {out_dir}')


def main():
    ap = argparse.ArgumentParser(description='Beam-spot + time-structure '
                                             'profiling.')
    ap.add_argument('run_key', nargs='?', default=sc.DEFAULT_RUN)
    ap.add_argument('--det', default=None,
                    help='station name/det_tag (default: all stations).')
    ap.add_argument('--sub-run', default=None,
                    help='sub_run (default: all discovered sub_runs).')
    ap.add_argument('--cluster-r', type=float, default=15.0)
    ap.add_argument('--beam-csv', default=None,
                    help='DAQ beam_monitor CSV (timestamp, unix_ts, '
                         'intensity_e10) to overlay on the rate-vs-time panel.')
    ap.add_argument('--veto-sparks', action=argparse.BooleanOptionalAction,
                    default=True)
    args = ap.parse_args()

    cfg = sc.get_config(args.run_key)
    print(cfg)
    args.suffix = cfg.product_suffix(args.veto_sparks)
    subruns = [args.sub_run] if args.sub_run else cfg.find_subruns()
    if not subruns:
        print('No sub_runs with combined hits.')
        return
    beam_df = load_beam_csv(args.beam_csv)
    if args.beam_csv and beam_df is None:
        print('  [beam-csv] not usable; continuing without the overlay.')

    dets = cfg.detectors()
    if args.det:
        dets = [d for d in dets if args.det in (d.name, d.det_tag)]
    for det in dets:
        if cfg.HAS_GEOMETRY:
            process_det(cfg, det, args, subruns, beam_df)
        else:
            process_det_chanspace(cfg, det, args, subruns, beam_df)


if __name__ == '__main__':
    main()
