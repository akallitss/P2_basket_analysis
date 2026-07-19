#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
20_beam_spectra.py

Per-sub_run cluster-charge spectra and Landau MPV vs HV / sub_run for the P2
telescope stations at the SPS beam test. Adapted from cosmic_bench_analysis/
18_fe55_spectra.py: the streaming/clustering, HV-settle cut, HV spark veto and
per-pad gain-map machinery are reused unchanged, but the signal is a minimum-
ionising beam particle (Landau-distributed cluster charge) rather than the Fe55
5.9 keV photopeak, so the per-point fit is a LANDAU MPV, not a Gaussian mean.

Landau approximation
--------------------
scipy has no Landau distribution; `scipy.stats.moyal` is the standard closed-
form approximation to the Landau (same asymmetric shape, exponential low side +
power-law high tail). Its location parameter `loc` IS the distribution mode, so
the fitted `loc` is taken as the MPV and `scale` as the width (Landau FWHM ~=
3.59 * scale). This choice is documented here so downstream readers know the
"MPV" is a moyal-loc, adequate for tracking gas gain / working points; it is
not a full Landau*Gauss convolution.

Method per point (sub_run)
--------------------------
  1. Channel table per station from the run-level run_config wiring.
  2. Stream the station's combined hits (its FEU only): HV-settle cut (mesh
     still ramping), optional HV spark veto (cfg.BURST_NPADS = 0 by default so
     high-multiplicity beam events are kept). One clean leading cluster per
     event (sps_cluster): leading pad + pads within --cluster-r.
  3. Cluster charge histogram -> Landau (moyal) MPV fit above --fit-min.
  4. Per-pad gain map (leading-pad median cluster charge) + per-pad hit map.

Products (<Analysis>/<det_tag>/<run>/<sub_run-or-scan>/20_beam_spectra/):
  spectra/spectrum_<pt><suffix>.png     per-point spectrum + Landau fit
  gain_maps/gain_map_<pt><suffix>.png   per-pad MPV proxy on the pad tiles (+csv)
  hit_maps/hit_map_<pt><suffix>.png     per-pad hit counts (+csv)
  beam_spectra_overlay<suffix>.png      all points, log-y
  beam_mpv_vs_hv<suffix>.png            Landau MPV vs mesh HV (or point index)
  beam_width_vs_hv<suffix>.png          Landau width / MPV vs HV
  beam_rate_vs_hv<suffix>.png           event rate + saturated fraction
  beam_mpv_vs_hv<suffix>.csv            machine-readable per-point table

Usage:
  python3 20_beam_spectra.py [run_key] [--det P2_OUT] [--cluster-r 15]
        [--fit-min 200] [--no-veto-sparks]
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
from scipy.stats import moyal

import sps_config as sc
import p2_mapping as pmap
import p2_io as p2io
import p2_sparks as ps
import sps_cluster as scl

_FWHM_MOYAL = 3.591   # Landau/moyal FWHM = _FWHM_MOYAL * scale


def _spark_veto(cfg, det, hv_csv):
    class _Shim:
        SPARK_CHANNEL = det.spark_channel
        SPARK_IMON_THR = cfg.SPARK_IMON_THR
        SPARK_GUARD_BEFORE = cfg.SPARK_GUARD_BEFORE
        SPARK_GUARD_AFTER = cfg.SPARK_GUARD_AFTER
        BURST_NPADS = cfg.BURST_NPADS
    return ps.SparkVeto.from_csv(hv_csv, _Shim)


def pad_hit_counts(hits_dir, ct, min_amp, t_min, veto):
    """Per-pad TOTAL hit counts (Series indexed by channel_id), streamed,
    after the same settle cut / spark veto used for the spectrum."""
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
    return (counts.astype(np.int64) if counts is not None
            else pd.Series(dtype=np.int64))


def _landau(x, A, loc, scale):
    return A * moyal.pdf(x, loc=loc, scale=scale)


def fit_landau(q, fit_min, nbins=200):
    """Landau (moyal) MPV fit of cluster charges above fit_min. Returns
    dict(mpv, mpv_err, scale, fwhm, ok, edges, counts)."""
    q = np.asarray(q, float)
    q = q[np.isfinite(q)]
    hi = np.quantile(q, 0.999) if len(q) else 1.0
    counts, edges = np.histogram(q, bins=nbins,
                                 range=(0.0, max(hi, fit_min * 2)))
    out = dict(mpv=np.nan, mpv_err=np.nan, scale=np.nan, fwhm=np.nan,
               ok=False, edges=edges, counts=counts)
    # fit on a finer histogram capped near the 99th pct so the long Landau
    # tail does not swamp the peak binning
    hi_fit = np.quantile(q, 0.99) if len(q) else 1.0
    fc, fe = np.histogram(q, bins=nbins, range=(0.0, max(hi_fit, fit_min * 2)))
    ctr = 0.5 * (fe[:-1] + fe[1:])
    sel = ctr >= fit_min
    if sel.sum() < 8 or fc[sel].sum() < 200:
        return out
    smooth = np.convolve(fc.astype(float), np.ones(5) / 5, mode='same')
    # seed the MPV at the smoothed peak (the rising edge above fit_min locates
    # the Landau maximum, not the low-charge noise floor)
    mpv = float(ctr[sel][np.argmax(smooth[sel])])
    cov = None
    # iterate a moyal fit in a window AROUND the peak: the Landau tail is
    # heavier than a single moyal, so fitting the whole [fit_min, end] range
    # would drag `loc` off the mode -- a peak window keeps loc == MPV.
    for _ in range(3):
        win = (ctr >= max(fit_min, 0.35 * mpv)) & (ctr <= 3.0 * mpv)
        if win.sum() < 6:
            return out
        try:
            p, cov = curve_fit(_landau, ctr[win], fc[win],
                               p0=[fc[win].max() * 0.3 * mpv, mpv, 0.15 * mpv],
                               maxfev=10000,
                               bounds=([0, 0, 1e-3],
                                       [np.inf, edges[-1], edges[-1]]))
        except (RuntimeError, ValueError):
            return out
        mpv, scale = float(p[1]), abs(float(p[2]))
        if not (0 < mpv < edges[-1]):
            return out
    out.update(mpv=mpv, mpv_err=float(np.sqrt(max(cov[1][1], 0.0))),
               scale=scale, fwhm=_FWHM_MOYAL * scale, ok=True)
    return out


def _windowed_median(q, fit_min, n_iter=3):
    q = q[q >= fit_min]
    if not len(q):
        return np.nan
    m = float(np.median(q))
    for _ in range(n_iter):
        w = q[(q > 0.5 * m) & (q < 2.0 * m)]
        if not len(w):
            break
        m = float(np.median(w))
    return m


def plot_map(values, ct, title, cbar, out_png, log=False, grey_below=None):
    """Generic per-pad map on the real pad tiles (values indexed by
    channel_id). grey_below draws pads with value < it as grey."""
    from matplotlib.collections import PolyCollection
    fig, ax = plt.subplots(figsize=(7.5, 6))
    if pmap.has_tile_geometry(ct):
        pads, verts = pmap.pad_tiles(ct)
        v = values.reindex(pads['channel_id']).to_numpy(dtype=float)
        good = np.isfinite(v) & (v > (grey_below if grey_below is not None
                                      else -np.inf))
        ax.add_collection(PolyCollection(verts[~good], facecolors='0.92',
                                         edgecolors='0.7', linewidths=0.3))
        if good.any():
            norm = matplotlib.colors.LogNorm(vmin=max(1, np.nanmin(v[good]))) \
                if log else None
            pc = PolyCollection(verts[good], array=v[good], cmap='viridis',
                                norm=norm, edgecolors='face', linewidths=0.2)
            ax.add_collection(pc)
            fig.colorbar(pc, ax=ax, label=cbar)
        ax.autoscale_view(); ax.set_aspect('equal')
    else:
        ax.scatter(ct['pad_cx'], ct['pad_cy'], c='0.7', s=6)
    ax.set_xlabel('pad_cx [mm]'); ax.set_ylabel('pad_cy [mm]')
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out_png, dpi=160, bbox_inches='tight')
    plt.close(fig)


def process_det(cfg, det, args):
    subruns = cfg.find_subruns()
    if not subruns:
        print('  no sub_runs with hits')
        return
    # products dir tag: the run may be a scan (many points) -> use 'scan',
    # else the single sub_run name
    prod_sub = 'scan' if len(subruns) > 1 else subruns[0]
    suffix = cfg.product_suffix(args.veto_sparks)
    base = cfg.out_dir(det.det_tag, prod_sub, '20_beam_spectra')
    spec_dir = cfg.out_dir(det.det_tag, prod_sub, '20_beam_spectra', 'spectra')
    map_dir = cfg.out_dir(det.det_tag, prod_sub, '20_beam_spectra', 'gain_maps')
    hit_dir = cfg.out_dir(det.det_tag, prod_sub, '20_beam_spectra', 'hit_maps')
    ct = cfg.channel_table(det)
    print(f'  {det.name}: FEUs {det.feus}, {ct["channel_id"].nunique()} pads')

    rows, spectra = [], []
    for i, sub in enumerate(subruns):
        hv = cfg.subrun_mesh_hv(sub, det)
        pt = hv if hv is not None else i
        lbl = f'{hv}V' if hv is not None else sub
        hits_dir = cfg.combined_hits_dir(sub)
        chunk0 = p2io.hit_files(hits_dir)[0]
        hv_csv = cfg.hv_monitor_csv(sub)
        t_min = scl.settle_t_min(hv_csv, det.spark_channel, chunk0)
        veto = (_spark_veto(cfg, det, hv_csv)
                if args.veto_sparks and os.path.isfile(hv_csv)
                and det.spark_channel else None)
        ev = scl.stream_event_clusters(hits_dir, ct, args.cluster_r,
                                       min_amp=args.min_amp, veto=veto,
                                       t_min=t_min)
        if not len(ev):
            print(f'    {lbl}: no events after cuts')
            continue
        live_s = float(ev['t'].max() - ev['t'].min())
        if veto is not None:
            live_s = max(0.0, live_s - veto.vetoed_seconds)
        rate = len(ev) / live_s if live_s > 0 else np.nan
        fit = fit_landau(ev['q'], args.fit_min)
        # per-pad gain proxy + hit map
        gp = ev[ev['q'] >= args.fit_min].groupby('lead_pad')['q']
        per_pad = gp.apply(lambda s: _windowed_median(s.to_numpy(),
                                                      args.fit_min))
        nev_pad = gp.size()
        per_pad = per_pad.where(nev_pad >= args.map_min_events)
        padc = ct.drop_duplicates('channel_id').set_index('channel_id')
        gaincsv = padc[['connector_N', 'pad_cx', 'pad_cy']].copy()
        gaincsv['n_events'] = nev_pad.reindex(gaincsv.index).fillna(0).astype(int)
        gaincsv['mpv_proxy_adc'] = per_pad.reindex(gaincsv.index)
        gaincsv.to_csv(os.path.join(map_dir, f'gain_map_{lbl}{suffix}.csv'))
        plot_map(per_pad, ct,
                 f'{det.name} per-pad MPV proxy — {lbl} ({sub})',
                 'median cluster charge [ADC]',
                 os.path.join(map_dir, f'gain_map_{lbl}{suffix}.png'))
        pad_hits = pad_hit_counts(hits_dir, ct, args.min_amp, t_min, veto)
        pd.DataFrame({'n_hits': pad_hits}).to_csv(
            os.path.join(hit_dir, f'hit_map_{lbl}{suffix}.csv'))
        plot_map(pad_hits.astype(float), ct,
                 f'{det.name} hit map — {lbl} ({sub}), '
                 f'{int(pad_hits.sum()):,} hits',
                 'pad hits (log)',
                 os.path.join(hit_dir, f'hit_map_{lbl}{suffix}.png'), log=True,
                 grey_below=0)

        rows.append(dict(point=pt, hv=hv, subrun=sub, n_events=len(ev),
                         live_s=live_s, rate_hz=rate,
                         mpv_adc=fit['mpv'], mpv_err=fit['mpv_err'],
                         width_adc=fit['scale'], fwhm_adc=fit['fwhm'],
                         fwhm_over_mpv=(fit['fwhm'] / fit['mpv']
                                       if fit['ok'] else np.nan),
                         median_nclus=float(ev['n_clus'].median()),
                         sat_frac=np.nan))
        spectra.append(dict(pt=pt, lbl=lbl, fit=fit, hv=hv))
        print(f'    {lbl}: {len(ev):,} ev, {rate:.0f} Hz, MPV = '
              + (f'{fit["mpv"]:.0f}+-{fit["mpv_err"]:.0f} ADC, '
                 f'FWHM/MPV = {fit["fwhm"]/fit["mpv"]:.2f}'
                 if fit['ok'] else 'FIT FAILED'))

        # per-point spectrum
        fig, ax = plt.subplots(figsize=(7, 4.5))
        ctr = 0.5 * (fit['edges'][:-1] + fit['edges'][1:])
        ax.step(ctr, fit['counts'], where='mid', color='steelblue', lw=1.2,
                label=f'cluster charge (N={len(ev):,})')
        if fit['ok']:
            xs = np.linspace(args.fit_min, fit['edges'][-1], 400)
            model = _landau(xs, 1.0, fit['mpv'], fit['scale'])
            model = model / model.max() * fit['counts'][
                np.argmin(np.abs(ctr - fit['mpv']))]
            ax.plot(xs, model, '-', color='crimson', lw=2,
                    label=f'Landau MPV {fit["mpv"]:.0f} ADC, '
                          f'FWHM/MPV {fit["fwhm"]/fit["mpv"]:.2f}')
        ax.axvline(args.fit_min, color='gray', ls=':', lw=1,
                   label=f'fit-min {args.fit_min:g}')
        ax.set_xlabel('cluster charge [ADC]'); ax.set_ylabel('events / bin')
        ax.set_yscale('log'); ax.grid(True, alpha=0.3); ax.legend(fontsize=8)
        ax.set_title(f'{det.name} beam spectrum — {lbl} ({sub})')
        fig.tight_layout()
        fig.savefig(os.path.join(spec_dir, f'spectrum_{lbl}{suffix}.png'),
                    dpi=160, bbox_inches='tight')
        plt.close(fig)

    if not rows:
        print('  no usable points')
        return
    df = pd.DataFrame(rows).sort_values('point').reset_index(drop=True)
    df.to_csv(os.path.join(base, f'beam_mpv_vs_hv{suffix}.csv'), index=False)
    has_hv = df['hv'].notna().any()
    xcol = 'hv' if has_hv else 'point'
    xlab = 'mesh HV [V]' if has_hv else 'sub_run index'
    tag = f'{det.name}  {cfg.RUN}'
    cmap = plt.get_cmap('viridis')

    # overlay
    fig, ax = plt.subplots(figsize=(8, 5))
    pts = [s['pt'] for s in spectra]
    lo, hi = min(pts), max(pts)
    for s in spectra:
        f = s['fit']
        ctr = 0.5 * (f['edges'][:-1] + f['edges'][1:])
        c = cmap((s['pt'] - lo) / max(1, hi - lo))
        ax.step(ctr, f['counts'], where='mid', color=c, lw=1.1, label=s['lbl'])
    ax.set_xlabel('cluster charge [ADC]'); ax.set_ylabel('events / bin')
    ax.set_yscale('log'); ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7, ncol=2, title=xlab)
    ax.set_title(f'{det.name} beam spectra — {tag}')
    fig.tight_layout()
    fig.savefig(os.path.join(base, f'beam_spectra_overlay{suffix}.png'),
                dpi=180, bbox_inches='tight'); plt.close(fig)

    ok = df[df['mpv_adc'].notna()]
    # MPV vs HV
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.errorbar(ok[xcol], ok['mpv_adc'], yerr=ok['mpv_err'], fmt='o-',
                color='steelblue', capsize=4, lw=2, ms=7, label='Landau MPV')
    if has_hv and len(ok) >= 3:
        k, b = np.polyfit(ok['hv'], np.log(ok['mpv_adc']), 1)
        xs = np.linspace(ok['hv'].min(), ok['hv'].max(), 100)
        ax.plot(xs, np.exp(b + k * xs), '--', color='crimson', lw=1.5,
                label=f'exp fit: x2 every {np.log(2)/k:.1f} V' if k > 0
                else 'exp fit')
        ax.set_yscale('log')
    ax.set_xlabel(xlab); ax.set_ylabel('Landau MPV [ADC]')
    ax.grid(True, alpha=0.3, which='both'); ax.legend()
    ax.set_title(f'{det.name} Landau MPV (moyal) vs {xlab} — {tag}')
    fig.tight_layout()
    fig.savefig(os.path.join(base, f'beam_mpv_vs_hv{suffix}.png'), dpi=180,
                bbox_inches='tight'); plt.close(fig)

    # width vs HV
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(ok[xcol], ok['fwhm_over_mpv'], 'o-', color='steelblue', lw=2, ms=7)
    ax.set_xlabel(xlab); ax.set_ylabel('Landau FWHM / MPV')
    ax.set_ylim(0, None); ax.grid(True, alpha=0.3)
    ax.set_title(f'{det.name} Landau width/MPV vs {xlab} — {tag}')
    fig.tight_layout()
    fig.savefig(os.path.join(base, f'beam_width_vs_hv{suffix}.png'), dpi=180,
                bbox_inches='tight'); plt.close(fig)

    # rate vs HV
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(df[xcol], df['rate_hz'], 'o-', color='steelblue', lw=2, ms=7)
    ax.set_xlabel(xlab); ax.set_ylabel('event rate [Hz]')
    ax.grid(True, alpha=0.3)
    ax.set_title(f'{det.name} event rate vs {xlab} — {tag}')
    fig.tight_layout()
    fig.savefig(os.path.join(base, f'beam_rate_vs_hv{suffix}.png'), dpi=180,
                bbox_inches='tight'); plt.close(fig)

    print(f'  -> {base}')


def _drift_of(sub):
    m = re.search(r'(-?\d+)', sub)
    return int(m.group(1)) if m else None


def process_det_chanspace(cfg, det, args):
    """Channel-space beam spectra when no pad geometry is available: a genuine
    MIP Landau of the leading-channel + neighbours cluster charge in (feu,
    channel) space (no x/y). Products go under <det_tag>/<run>/<sub>/
    20_beam_spectra_chanspace/."""
    subruns = cfg.find_subruns()
    if not subruns:
        print('  no sub_runs with hits')
        return
    print(f'  {det.name}: CHANNEL-SPACE mode (no pad map), FEUs {det.feus}')
    rows, spectra = [], []
    for sub in subruns:
        out_dir = cfg.out_dir(det.det_tag, sub, '20_beam_spectra_chanspace')
        hits_dir = cfg.combined_hits_dir(sub)
        try:
            ev = scl.stream_event_charge_chanspace(
                hits_dir, det.feus, cluster_chan=args.cluster_chan,
                min_amp=args.min_amp, progress=True)
        except Exception as e:  # noqa: BLE001 (incomplete rsync = no 'hits')
            print(f'    {sub}: unreadable (incomplete transfer?), skipping — {e}')
            continue
        if not len(ev):
            print(f'    {sub}: no events')
            continue
        live = float(ev['t'].max() - ev['t'].min())
        rate = len(ev) / live if live > 0 else np.nan
        fit = fit_landau(ev['q'], args.fit_min)
        drift = _drift_of(sub)
        rows.append(dict(sub_run=sub, drift=drift, n_events=len(ev),
                         rate_hz=rate, mpv_adc=fit['mpv'],
                         mpv_err=fit['mpv_err'], width_adc=fit['scale'],
                         fwhm_adc=fit['fwhm'],
                         fwhm_over_mpv=(fit['fwhm'] / fit['mpv']
                                       if fit['ok'] else np.nan),
                         median_nchan=float(ev['n_chan'].median()),
                         median_nhit=float(ev['n_hit'].median())))
        spectra.append((sub, fit))
        print(f'    {sub}: {len(ev):,} ev, {rate:.0f} Hz, channel-space Landau '
              f'MPV = '
              + (f'{fit["mpv"]:.0f}+-{fit["mpv_err"]:.0f} ADC, '
                 f'FWHM/MPV = {fit["fwhm"]/fit["mpv"]:.2f}'
                 if fit['ok'] else 'FIT FAILED'))
        # per-point spectrum
        fig, ax = plt.subplots(figsize=(7.5, 4.8))
        ctr = 0.5 * (fit['edges'][:-1] + fit['edges'][1:])
        ax.step(ctr, fit['counts'], where='mid', color='steelblue', lw=1.2,
                label=f'cluster charge (N={len(ev):,})')
        if fit['ok']:
            xs = np.linspace(args.fit_min, fit['edges'][-1], 400)
            model = _landau(xs, 1.0, fit['mpv'], fit['scale'])
            model = model / model.max() * fit['counts'][
                np.argmin(np.abs(ctr - fit['mpv']))]
            ax.plot(xs, model, '-', color='crimson', lw=2,
                    label=f'Landau MPV {fit["mpv"]:.0f} ADC, '
                          f'FWHM/MPV {fit["fwhm"]/fit["mpv"]:.2f}')
        ax.axvline(args.fit_min, color='gray', ls=':', lw=1,
                   label=f'fit-min {args.fit_min:g}')
        ax.set_xlabel('leading-cluster charge [ADC] (channel space, '
                      f'+-{args.cluster_chan} ch)')
        ax.set_ylabel('events / bin'); ax.set_yscale('log')
        ax.grid(True, alpha=0.3); ax.legend(fontsize=8)
        ax.set_title(f'{det.name} beam MIP spectrum (channel space) — {sub}')
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir,
                                 f'spectrum_chanspace_{sub}.png'),
                    dpi=160, bbox_inches='tight')
        plt.close(fig)

    if not rows:
        print('  no usable points')
        return
    base = cfg.out_dir(det.det_tag, 'scan' if len(subruns) > 1 else subruns[0],
                       '20_beam_spectra_chanspace')
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(base, 'beam_mpv_vs_drift_chanspace.csv'),
              index=False)
    if len(df) > 1 and df['drift'].notna().all():
        fig, ax = plt.subplots(figsize=(7.5, 4.8))
        ax.errorbar(df['drift'], df['mpv_adc'], yerr=df['mpv_err'], fmt='o-',
                    color='steelblue', capsize=4, lw=2, ms=7)
        ax.set_xlabel('drift HV [V]'); ax.set_ylabel('channel-space Landau MPV '
                                                    '[ADC]')
        ax.grid(True, alpha=0.3)
        ax.set_title(f'{det.name} MIP Landau MPV vs drift HV (channel space)')
        fig.tight_layout()
        fig.savefig(os.path.join(base, 'beam_mpv_vs_drift_chanspace.png'),
                    dpi=170, bbox_inches='tight')
        plt.close(fig)
    print(f'  -> {base}')


def main():
    ap = argparse.ArgumentParser(description='Beam cluster spectra + Landau MPV '
                                             'vs HV / sub_run.')
    ap.add_argument('run_key', nargs='?', default=sc.DEFAULT_RUN)
    ap.add_argument('--det', default=None,
                    help='station name or det_tag (default: all stations).')
    ap.add_argument('--strategy', default='reverse',
                    choices=['linear', 'reverse', 'pairswap'])
    ap.add_argument('--cluster-r', type=float, default=15.0)
    ap.add_argument('--cluster-chan', type=int, default=3,
                    help='channel-space mode: +-this many channels around the '
                         'leading channel (same FEU) form the cluster.')
    ap.add_argument('--min-amp', type=float, default=None)
    ap.add_argument('--fit-min', type=float, default=200.0,
                    help='Landau fit starts above this cluster charge [ADC].')
    ap.add_argument('--map-min-events', type=int, default=30)
    ap.add_argument('--veto-sparks', action=argparse.BooleanOptionalAction,
                    default=True)
    args = ap.parse_args()

    cfg = sc.get_config(args.run_key)
    cfg.STRATEGY = args.strategy
    print(cfg)
    if args.min_amp is None:
        args.min_amp = cfg.MIN_AMP

    dets = cfg.detectors()
    if args.det:
        dets = [d for d in dets if args.det in (d.name, d.det_tag)]
        if not dets:
            print(f'No station matching {args.det!r}')
            return
    for det in dets:
        if cfg.HAS_GEOMETRY:
            process_det(cfg, det, args)
        else:
            process_det_chanspace(cfg, det, args)


if __name__ == '__main__':
    main()
