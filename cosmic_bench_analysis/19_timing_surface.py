#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
19_timing_surface.py

Surface timing map: per-pad peak time (time_of_max) across the detector plane.

Stage 13 characterises timing ALGORITHMS and resolution but is spatially blind.
This stage asks the orthogonal question -- does the timing depend on WHERE the
muon crossed?

What structure means
--------------------
In a parallel-plate geometry the drift time depends on the DEPTH z at which the
primary cluster forms, not on (x, y). So the physics prediction for this map is
FLAT: a Garfield/HEED simulation of the drift would predict no (x, y) structure
and cannot be used to interpret one. Anything seen here is therefore
instrumental, and that is the point of the map:

  * per-connector / per-FEU steps  -> cable-length and front-end delays
  * smooth gradients               -> a tilted detector plane, or a residual
                                      trigger-time correlation with position
  * isolated pads                  -> a channel with a bad baseline or shaping

CAVEAT: combined_hits carries `time_of_max` WITHOUT the DREAM fine-timestamp
(ftst) correction -- ftst lives in decoded_root and is applied in stage 13.
ftst is a per-EVENT phase, so it shifts hits coherently and largely averages
out in a per-pad mean over many events; it is not a per-channel offset. Treat
absolute values as relative, not calibrated.

Memory: per-pad time histograms (n_pads x n_bins), so the pass is streamed and
bounded like every other reduction here -- no per-hit array is held.

Products (<Analysis>/<detN>/<run>/<sub_run>/19_timing_surface/):
  timing_surface<sfx>.png      median peak time + spread per pad, on the PCB
  timing_per_connector<sfx>.png  offsets grouped by connector, the instrumental view
  timing_surface<sfx>.csv      per-pad n / median / p10 / p90 / iqr

Usage:
  python3 19_timing_surface.py [run_key] [--sig-amp 300] [--bins 1] [--no-veto-sparks]
"""

import argparse
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import p2_qa_config as qa
qa.setup_paths()
import p2_io as p2io
import p2_mapping as pmap
import p2_sparks as ps
from scipy.spatial import cKDTree

# Histogram window for the peak time [ns]. The signal band sits around
# 420-480 ns on det1 (DREAM 32-sample window), but time_of_max has a long tail
# and occasional large negative outliers, so the default spans the whole
# plausible range and the stage reports what falls outside it -- a clipped
# window silently biases the quantiles (a 0-400 ns window put every pad in the
# top bin and produced a spurious 0.01 ns "offset spread").
T_LO, T_HI = 0.0, 2000.0


def quantile_from_hist(counts, edges, q):
    """Linear-interpolated quantile from a histogram row."""
    tot = counts.sum()
    if tot <= 0:
        return np.nan
    c = np.cumsum(counts)
    target = q * tot
    i = int(np.searchsorted(c, target))
    if i >= len(counts):
        return edges[-1]
    c0 = c[i - 1] if i else 0.0
    w = counts[i]
    frac = (target - c0) / w if w > 0 else 0.0
    return edges[i] + frac * (edges[i + 1] - edges[i])


def sliding_timing_map(ev_xy, ev_t, x_grid, y_grid, kernel, min_ev):
    """Median offset and Gaussian-core sigma of the per-event pad time over the
    events inside `kernel` mm of each grid point (overlapping -> smooth)."""
    tree = cKDTree(ev_xy)
    gx, gy = np.meshgrid(x_grid, y_grid, indexing='ij')
    pts = np.column_stack([gx.ravel(), gy.ravel()])
    med = np.full(len(pts), np.nan)
    sig = np.full(len(pts), np.nan)
    cnt = np.full(len(pts), np.nan)   # NaN outside the acceptance -> blank
    for i, idx in enumerate(tree.query_ball_point(pts, kernel)):
        if len(idx):
            cnt[i] = len(idx)      # record occupancy even below min_ev
        if len(idx) < min_ev:
            continue
        t = ev_t[idx]
        p16, p50, p84 = np.percentile(t, [15.865, 50.0, 84.135])
        med[i] = p50
        sig[i] = 0.5 * (p84 - p16)
    sh = (len(x_grid), len(y_grid))
    return med.reshape(sh), sig.reshape(sh), cnt.reshape(sh)


def build_continuous(cfg, args, sfx, out, drop, ct, bad):
    """Sliding-window timing map on the M3 reference frame.

    Uses the per-ray list written by stage 06 (same x,y frame as the efficiency
    map) joined to a per-event pad time, so the two maps overlay exactly.

    RESOLUTION LIMIT: each event contributes the time of ONE pad, so the map
    cannot resolve structure finer than the 11.8 mm pad pitch, and the M3 track
    position carries ~5 mm of its own. A kernel below ~5 mm therefore only
    oversamples -- it looks smoother without adding information.
    """
    import glob as _glob
    cand = os.path.join(cfg.OUT_BASE, '06_efficiency',
                        f'ray_hit_miss_list{sfx}.csv')
    if not os.path.isfile(cand):
        alt = sorted(_glob.glob(os.path.join(cfg.OUT_BASE, '06_efficiency',
                                             'ray_hit_miss_list*.csv')))
        if not alt:
            print('  !! no ray_hit_miss_list from stage 06 -- run 06 first')
            return
        cand = alt[-1]
        print(f'  !! exact suffix not found, using {os.path.basename(cand)}')
    rays = pd.read_csv(cand)
    rays = rays[rays['in_active']] if 'in_active' in rays.columns else rays
    print(f'  continuous: {len(rays):,} rays from {os.path.basename(cand)}')

    # per-event leading-pad time (highest-amplitude signal-band pad)
    best_t, best_a = {}, {}
    for df in p2io.iter_hits(cfg.combined_hits_dir,
                             ['eventId', 'channel', 'feu', 'amplitude',
                              'time_of_max'],
                             ct.attrs['feus'], progress=False,
                             min_amp=cfg.MIN_AMP,
                             t_max_h=cfg.T_MAX_H, t_min_h=cfg.T_MIN_H):
        h = pmap.attach_pads_to_hits(df, ct)
        h = h[h['mapped'] & h['pad_cx'].notna() & (h['amplitude'] >= args.sig_amp)]
        if drop:
            h = h[~h['channel_id'].isin(drop)]
        if bad:
            h = h[~h['eventId'].isin(bad)]
        if not len(h):
            continue
        sel = h.loc[h.groupby('eventId')['amplitude'].idxmax(),
                    ['eventId', 'amplitude', 'time_of_max']]
        for ev, a, t in sel.itertuples(index=False):
            if ev not in best_a or a > best_a[ev]:
                best_a[ev] = float(a); best_t[ev] = float(t)
    tt = pd.DataFrame({'eventId': list(best_t.keys()),
                       't_pad': list(best_t.values())})
    m = rays.merge(tt, on='eventId', how='inner')
    if len(m) < 200:
        print(f'  !! only {len(m)} events with both a ray and a pad time')
        return
    ref = float(np.median(m['t_pad']))
    print(f'  continuous: {len(m):,} matched events, reference {ref:.1f} ns, '
          f'kernel {args.kernel:g} mm, grid {args.grid}')

    xy = m[['x', 'y']].to_numpy()
    t = m['t_pad'].to_numpy() - ref
    xg = np.linspace(xy[:, 0].min(), xy[:, 0].max(), args.grid)
    yg = np.linspace(xy[:, 1].min(), xy[:, 1].max(), args.grid)

    # Auto-widen the kernel for low-statistics runs. A fixed 6 mm kernel needs
    # ~80k matched events to fill; det3_initial1 has 4.9k and filled 0/25600
    # grid points, silently saving a blank figure. Grow the kernel until a
    # usable fraction fills, and put the kernel actually used on the plot so a
    # heavily smoothed map is never mistaken for a finely resolved one.
    kernel = args.kernel
    med, sig, cnt = sliding_timing_map(xy, t, xg, yg, kernel, args.min_ev)
    filled = int(np.isfinite(med).sum())
    while filled < 0.05 * med.size and kernel < args.max_kernel:
        kernel = min(kernel * 1.6, args.max_kernel)
        print(f'  continuous: only {filled}/{med.size} filled -> widening '
              f'kernel to {kernel:.1f} mm')
        med, sig, cnt = sliding_timing_map(xy, t, xg, yg, kernel, args.min_ev)
        filled = int(np.isfinite(med).sum())
    print(f'  continuous: {filled}/{med.size} grid points filled '
          f'(kernel {kernel:.1f} mm)')
    if filled == 0:
        print('  !! no grid point reached --min-ev; NOT writing an empty map. '
              f'Only {len(m):,} matched events -- lower --min-ev or raise '
              '--max-kernel if a coarse map is still wanted.')
        return

    ext = [xg[0], xg[-1], yg[0], yg[-1]]
    # NaN must render as blank in every panel. Left as-is, a 0-valued count
    # paints the whole out-of-acceptance region black and buries the edge.
    cmaps = {}
    for nm in ('coolwarm', 'viridis', 'magma'):
        c = matplotlib.colormaps[nm].copy()
        c.set_bad('white')
        cmaps[nm] = c
    fig, axs = plt.subplots(1, 3, figsize=(19, 5.6))
    lim = float(np.nanpercentile(np.abs(med), 95)) or 1.0
    im0 = axs[0].imshow(med.T, origin='lower', extent=ext, cmap=cmaps['coolwarm'],
                        vmin=-lim, vmax=lim, interpolation='bilinear')
    axs[0].set_title(f'Peak-time offset (median − {ref:.1f} ns)')
    plt.colorbar(im0, ax=axs[0], label='offset [ns]')
    im1 = axs[1].imshow(sig.T, origin='lower', extent=ext, cmap=cmaps['viridis'],
                        vmin=float(np.nanpercentile(sig, 2)),
                        vmax=float(np.nanpercentile(sig, 98)),
                        interpolation='bilinear')
    axs[1].set_title('Time resolution  σ=(p84.1−p15.9)/2')
    plt.colorbar(im1, ax=axs[1], label='σ [ns]')
    im2 = axs[2].imshow(cnt.T, origin='lower', extent=ext, cmap=cmaps['magma'],
                        vmin=0, vmax=float(np.nanpercentile(cnt, 99)),
                        interpolation='nearest')
    axs[2].set_title(f'events per kernel (r={kernel:g} mm)\n'
                     f'blank = no events; filled needs >= {args.min_ev}')
    plt.colorbar(im2, ax=axs[2], label='events')
    for a in axs:
        a.set_xlabel('reference X [mm]'); a.set_ylabel('reference Y [mm]')
        a.set_aspect('equal')
    fig.suptitle(f'{cfg.DET_NAME} sliding timing map — {cfg.RUN}/{cfg.SUB_RUN}\n'
                 f'M3 reference frame (same as the stage-10 efficiency map); '
                 f'kernel {kernel:.1f} mm on {len(m):,} events; structure finer '
                 f'than the 11.8 mm pad pitch is not resolved',
                 fontsize=10)
    fig.tight_layout()
    f = os.path.join(out, f'timing_map_sliding{sfx}.png')
    fig.savefig(f, dpi=150, bbox_inches='tight'); plt.close(fig)
    print('  saved', f)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run_key', nargs='?', default=qa.DEFAULT_RUN)
    ap.add_argument('--sig-amp', type=float, default=300.0,
                    help='use only signal-band hits above this amplitude [ADC]')
    ap.add_argument('--bins', type=float, default=1.0, help='histogram bin [ns]')
    ap.add_argument('--t-lo', type=float, default=T_LO, help='peak-time window low [ns]')
    ap.add_argument('--t-hi', type=float, default=T_HI, help='peak-time window high [ns]')
    ap.add_argument('--min-hits', type=int, default=20,
                    help='pads with fewer hits are left blank')
    ap.add_argument('--veto-sparks', action=argparse.BooleanOptionalAction,
                    default=True)
    ap.add_argument('--continuous', action='store_true',
                    help='also build a sliding-window map in the M3 reference '
                         'frame, directly comparable with the stage-10 '
                         'efficiency map')
    ap.add_argument('--kernel', type=float, default=6.0,
                    help='sliding kernel radius [mm] for --continuous')
    ap.add_argument('--grid', type=int, default=160, help='grid points per axis')
    ap.add_argument('--max-kernel', type=float, default=30.0,
                    help='cap when auto-widening the kernel on sparse runs')
    ap.add_argument('--min-ev', type=int, default=25,
                    help='minimum events in a kernel to fill a grid point')
    args = ap.parse_args()

    cfg = qa.get_config(args.run_key)
    print(cfg)
    out = cfg.out_dir('19_timing_surface')
    sfx = cfg.product_suffix(args.veto_sparks)

    ct = pmap.build_channel_table(cfg.run_config_path, cfg.MAP_CSV_PATH,
                                  det_type=cfg.DET_TYPE, det_name=cfg.DET_NAME,
                                  drop_connectors=cfg.DEAD_CONNECTORS)
    drop = set(p2io.drop_pads_for(cfg, ct))
    print(f'  FEUs {ct.attrs["feus"]}  dead connectors {list(cfg.DEAD_CONNECTORS)}  '
          f'dropped pads {len(drop)}')

    sv = ps.SparkVeto.from_cfg(cfg) if args.veto_sparks else None
    bad = set()
    if sv is not None:
        bad = set(int(b) for b in
                  sv.vetoed_ids_from_hits(cfg.combined_hits_dir,
                                          ct.attrs['feus'], min_amp=cfg.MIN_AMP))
        print(f'  spark veto: {len(bad):,} events removed')

    edges = np.arange(args.t_lo, args.t_hi + args.bins, args.bins)
    nb = len(edges) - 1
    pads = ct.drop_duplicates('channel_id')
    pads = pads[pads['mapped'] & pads['pad_cx'].notna()].set_index('channel_id')
    pad_ids = list(pads.index)
    idx = {p: i for i, p in enumerate(pad_ids)}
    hist = np.zeros((len(pad_ids), nb), dtype=np.int64)
    # Split-half (even/odd eventId) copies: two statistically independent
    # estimates of the same per-pad offset. Their difference measures the
    # STATISTICAL error directly, so the map can be tested for significance
    # instead of asserted -- a per-pad median from ~200 hits of a ~21 ns
    # distribution already carries ~2 ns of noise on its own.
    hist_a = np.zeros_like(hist)
    hist_b = np.zeros_like(hist)
    n_out, n_tot = [0], [0]      # hits falling outside the histogram window

    for df in p2io.iter_hits(cfg.combined_hits_dir,
                             ['channel', 'feu', 'amplitude', 'time_of_max',
                              'eventId'],
                             ct.attrs['feus'], progress=True,
                             min_amp=cfg.MIN_AMP,
                             t_max_h=cfg.T_MAX_H, t_min_h=cfg.T_MIN_H):
        h = pmap.attach_pads_to_hits(df, ct)
        h = h[h['mapped'] & h['pad_cx'].notna()]
        if bad:
            h = h[~h['eventId'].isin(bad)]
        if drop:
            h = h[~h['channel_id'].isin(drop)]
        h = h[h['amplitude'] >= args.sig_amp]
        if not len(h):
            continue
        rows = h['channel_id'].map(idx).to_numpy()
        good = ~pd.isna(rows)
        rows = rows[good].astype(int)
        tom = h['time_of_max'].to_numpy()[good].astype(float)
        inside = (tom >= args.t_lo) & (tom < args.t_hi)
        n_out[0] += int((~inside).sum()); n_tot[0] += int(inside.size)
        rows = rows[inside]
        cols = ((tom[inside] - args.t_lo) / args.bins).astype(int)
        np.add.at(hist, (rows, cols), 1)
        ev = h['eventId'].to_numpy()[good][inside]
        half = (ev % 2) == 0
        np.add.at(hist_a, (rows[half], cols[half]), 1)
        np.add.at(hist_b, (rows[~half], cols[~half]), 1)

    if n_tot[0]:
        frac = 100.0 * n_out[0] / n_tot[0]
        print(f'  peak-time window {args.t_lo:g}-{args.t_hi:g} ns: '
              f'{n_out[0]:,}/{n_tot[0]:,} hits outside ({frac:.2f} %)')
        if frac > 5:
            print('  !! >5 % outside the window -- widen --t-lo/--t-hi, '
                  'the quantiles below are biased')
    n = hist.sum(axis=1)
    med = np.array([quantile_from_hist(hist[i], edges, 0.5) for i in range(len(pad_ids))])
    p10 = np.array([quantile_from_hist(hist[i], edges, 0.10) for i in range(len(pad_ids))])
    p90 = np.array([quantile_from_hist(hist[i], edges, 0.90) for i in range(len(pad_ids))])
    p25 = np.array([quantile_from_hist(hist[i], edges, 0.25) for i in range(len(pad_ids))])
    p75 = np.array([quantile_from_hist(hist[i], edges, 0.75) for i in range(len(pad_ids))])
    # Gaussian-core sigma, identical definition to p2_io.event_time_resolution
    # (stage 16's time_res_ns), so the per-pad map is directly comparable with
    # the ~21 ns quoted by the drift scan. p90-p10 = 2.563 sigma for a Gaussian,
    # which is why the raw spread reads ~55 ns for the same distribution.
    p16 = np.array([quantile_from_hist(hist[i], edges, 0.15865) for i in range(len(pad_ids))])
    p84 = np.array([quantile_from_hist(hist[i], edges, 0.84135) for i in range(len(pad_ids))])
    sigma = 0.5 * (p84 - p16)
    med_a = np.array([quantile_from_hist(hist_a[i], edges, 0.5) for i in range(len(pad_ids))])
    med_b = np.array([quantile_from_hist(hist_b[i], edges, 0.5) for i in range(len(pad_ids))])
    ok = n >= args.min_hits
    med[~ok] = np.nan; p10[~ok] = np.nan; p90[~ok] = np.nan
    p25[~ok] = np.nan; p75[~ok] = np.nan
    sigma[~ok] = np.nan; med_a[~ok] = np.nan; med_b[~ok] = np.nan

    d = pd.DataFrame({'channel_id': pad_ids, 'n': n, 'median_ns': med,
                      'p10_ns': p10, 'p90_ns': p90,
                      'spread_ns': p90 - p10, 'iqr_ns': p75 - p25,
                      'sigma_ns': sigma, 'median_even_ns': med_a,
                      'median_odd_ns': med_b})
    d = d.join(pads[['pad_cx', 'pad_cy', 'connector_N']], on='channel_id')
    csv = os.path.join(out, f'timing_surface{sfx}.csv')
    d.to_csv(csv, index=False)

    live = d[d['n'] >= args.min_hits]
    if not len(live):
        print('!! no pad passed --min-hits; nothing to plot')
        return
    ref = float(np.nanmedian(live['median_ns']))
    print(f'  pads with >= {args.min_hits} hits: {len(live)}/{len(d)}   '
          f'detector median peak time {ref:.1f} ns')
    off = live['median_ns'] - ref
    print(f'  per-pad peak-time sigma (p84.1-p15.9)/2: '
          f'median {np.nanmedian(live["sigma_ns"]):.1f} ns  '
          f'[stage 16 event-level time_res is the same definition]')
    print(f'  per-pad offset spread: sigma {np.nanstd(off):.2f} ns, '
          f'p10..p90 {np.nanpercentile(off,10):.2f}..{np.nanpercentile(off,90):.2f} ns')

    # --- is the offset structure real, or just per-pad statistics? --------- #
    # Two independent halves (even/odd eventId) of the same pads:
    #   var(diff) = 2 * stat^2                  -> stat error per pad
    #   var(mean) = real^2 + stat^2/2           -> genuine pad-to-pad spread
    sh = live.dropna(subset=['median_even_ns', 'median_odd_ns'])
    real = stat = corr = float('nan')
    if len(sh) > 20:
        diff = (sh['median_even_ns'] - sh['median_odd_ns']).to_numpy()
        mean = (0.5 * (sh['median_even_ns'] + sh['median_odd_ns'])).to_numpy()
        stat = float(np.nanstd(diff) / np.sqrt(2.0))
        var_real = float(np.nanvar(mean) - 0.5 * stat ** 2)
        real = float(np.sqrt(var_real)) if var_real > 0 else 0.0
        corr = float(np.corrcoef(sh['median_even_ns'], sh['median_odd_ns'])[0, 1])
        print(f'  SPLIT-HALF TEST on {len(sh)} pads:')
        print(f'    statistical error per pad : {stat:.2f} ns')
        print(f'    genuine pad-to-pad spread : {real:.2f} ns')
        print(f'    even/odd correlation      : {corr:.3f} '
              f'({"structure reproduces" if corr > 0.3 else "NOT reproducible -- consistent with noise"})')

    # ---- surface map -------------------------------------------------- #
    fig, axs = plt.subplots(1, 2, figsize=(14.5, 6))
    lim = float(np.nanpercentile(np.abs(off), 95)) or 1.0
    s0 = axs[0].scatter(live['pad_cx'], live['pad_cy'], c=off, s=16,
                        marker='s', cmap='coolwarm', vmin=-lim, vmax=lim)
    axs[0].set_title(f'Median peak time − detector median ({ref:.1f} ns)\n'
                     f'genuine spread {real:.1f} ns vs {stat:.1f} ns statistical '
                     f'(split-half r={corr:.2f})')
    plt.colorbar(s0, ax=axs[0], label='offset [ns]')
    s1 = axs[1].scatter(live['pad_cx'], live['pad_cy'], c=live['sigma_ns'],
                        s=16, marker='s', cmap='viridis',
                        vmin=float(np.nanpercentile(live['sigma_ns'], 2)),
                        vmax=float(np.nanpercentile(live['sigma_ns'], 98)))
    axs[1].set_title(f'Per-pad time resolution  σ=(p84.1−p15.9)/2\n'
                     f'median {np.nanmedian(live["sigma_ns"]):.1f} ns '
                     f'(same definition as stage 16)')
    plt.colorbar(s1, ax=axs[1], label='σ [ns]')
    for a in axs:
        a.set_xlabel('pad_cx [mm]'); a.set_ylabel('pad_cy [mm]'); a.set_aspect('equal')
    fig.suptitle(f'{cfg.DET_NAME} pad timing surface — {cfg.RUN}/{cfg.SUB_RUN}\n'
                 f'signal band > {args.sig_amp:g} ADC; drift time depends on depth z, '
                 f'not on (x,y) — structure here is instrumental', fontsize=10)
    fig.tight_layout()
    f1 = os.path.join(out, f'timing_surface{sfx}.png')
    fig.savefig(f1, dpi=150, bbox_inches='tight'); plt.close(fig)

    # ---- per-connector view ------------------------------------------- #
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.6))
    g = live.groupby('connector_N')['median_ns']
    conns = sorted(g.groups)
    axs[0].errorbar(conns, [g.get_group(c).median() - ref for c in conns],
                    yerr=[g.get_group(c).std() for c in conns],
                    fmt='o-', color='steelblue', capsize=3)
    axs[0].axhline(0, color='grey', ls=':')
    axs[0].set_xlabel('physical connector'); axs[0].set_ylabel('median offset [ns]')
    axs[0].set_title('Timing offset per connector (cable / front-end delays)')
    axs[0].grid(alpha=0.3)
    axs[1].hist(off.dropna(), bins=40, color='steelblue', alpha=0.85)
    axs[1].set_xlabel('per-pad offset from detector median [ns]')
    axs[1].set_ylabel('pads'); axs[1].grid(alpha=0.3)
    axs[1].set_title(f'σ = {np.nanstd(off):.2f} ns over {len(live)} pads')
    fig.suptitle(f'{cfg.DET_NAME} timing offsets — {cfg.RUN}/{cfg.SUB_RUN}', fontsize=10)
    fig.tight_layout()
    f2 = os.path.join(out, f'timing_per_connector{sfx}.png')
    fig.savefig(f2, dpi=150, bbox_inches='tight'); plt.close(fig)

    print('  saved', f1)
    print('  saved', f2)
    print('  saved', csv)

    if args.continuous:
        build_continuous(cfg, args, sfx, out, drop, ct, bad)


if __name__ == '__main__':
    main()
