#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
18_fe55_spectra.py

Fe55 source spectra and gain vs mesh HV for the self-triggered banco runs
(run keys det2_fe55scan1 / det3_fe55scan1). There is NO reference tracker
here: the DREAM DAQ self-triggers on the Fe55 source (TCM multiplicity), so
instead of efficiency the scan observable is the Fe55 photopeak position
(~gain), the energy resolution and the trigger rate, per mesh-HV point.

Method per HV point (sub_run fe55_<NN>_mesh_out<V>_mid<V>):
  1. Channel table from the run-level run_config (P2_OUT / P2_MID wiring),
     built once.
  2. Stream the sub_run's combined hits (this detector's FEU only), with
     - an HV-settle cut: events taken before the mesh vmon first reaches its
       set point (hv_monitor.csv) are dropped (matters for the first point,
       which ramps 200 V -> setpoint);
     - an optional per-sub_run HV spark veto (mesh imon spikes, as stage 11);
       the cosmic burst veto is disabled via cfg.BURST_NPADS = 0.
  3. Per-event clustering: leading pad + every pad within --cluster-r mm of
     it -> cluster charge q_clus (Fe55 estimator); the leading-pad amplitude
     alone is kept as a cross-check spectrum. Because both detectors share
     the trigger, events where only the OTHER detector converted a photon
     leave just noise here -- they populate a low-charge bulge that the
     photopeak fit must ignore (--fit-min).
  4. Photopeak: mode of the smoothed q_clus histogram above --fit-min, then
     an iterative Gaussian fit in a +-30% window around the running mean.

Products (<Analysis>/<detN>/<run>/fe55_scan/18_fe55_spectra/):
  spectra/spectrum_<HV>V<suffix>.png   per-point spectrum + photopeak fit
  gain_maps/gain_map_<HV>V<suffix>.png per-pad Fe55 peak on the pad tiles
  gain_maps/gain_map_<HV>V<suffix>.csv (+ per-pad leading-event counts)
  hit_maps/hit_map_<HV>V<suffix>.png   per-pad hit counts, linear + log
  hit_maps/hit_map_<HV>V<suffix>.csv
  fe55_spectra_overlay<suffix>.png     all HV points, log-y
  fe55_gain_vs_hv<suffix>.png          peak vs mesh HV, semilog-y + exp fit
  fe55_resolution_vs_hv<suffix>.png    FWHM/peak vs mesh HV
  fe55_rate_vs_hv<suffix>.png          trigger rate + saturated fraction
  fe55_gain_vs_hv<suffix>.csv

The gain map assigns each event to its LEADING pad and estimates that pad's
photopeak with an iterative windowed median of the cluster charge (stable at
the ~50-100 events/pad a 5-min point gives; a per-pad Gaussian fit is not).
Pads with fewer than --map-min-events events are drawn grey.

Usage:
  python3 18_fe55_spectra.py [det2_fe55scan1] [--cluster-r 15] [--min-amp 0]
        [--fit-min 200] [--no-veto-sparks]
"""

import os
import re
import glob
import argparse
import matplotlib
matplotlib.use('Agg')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import p2_qa_config as qa
import p2_mapping as pmap
import p2_io as p2io
import p2_sparks as ps

# which mesh voltage in fe55_<NN>_mesh_out<V>_mid<V> belongs to which detector
_MESH_GROUP = {'P2_OUT': 1, 'P2_MID': 2}
_SUBRUN_RE = re.compile(r'fe55_\d+_mesh_out(\d+)_mid(\d+)')


def find_subruns(cfg):
    """Discover Fe55-scan sub_runs with hits, ascending in THIS detector's
    mesh HV. Returns [(sub_run_name, mesh_hv)]."""
    grp = _MESH_GROUP[cfg.DET_NAME]
    out = []
    for name in sorted(os.listdir(cfg.run_dir)):
        m = _SUBRUN_RE.search(name)
        if not m:
            continue
        hits = glob.glob(os.path.join(cfg.subrun_dir(name),
                                      'combined_hits_root', '*.root'))
        if hits:
            out.append((name, int(m.group(grp))))
        else:
            print(f'  [skip] {name}: no combined hits yet (still processing?)')
    return sorted(out, key=lambda x: x[1])


def settle_t_min(cfg, sub_dir, chunk0, margin=5.0):
    """Event-time [s since DAQ start] before which the mesh HV was still
    ramping. hv_monitor.csv starts at HV-set time; the data chunks start
    later (pedestal run first). The chunk-0 filename carries the DAQ start
    to the minute -- assume the start of that minute (conservative: cuts a
    little extra data rather than keeping ramp data)."""
    hv_csv = os.path.join(sub_dir, 'hv_monitor.csv')
    if not os.path.isfile(hv_csv):
        return 0.0
    df = pd.read_csv(hv_csv)
    vcol, scol = f'{cfg.SPARK_CHANNEL} vmon', f'{cfg.SPARK_CHANNEL} v0'
    if vcol not in df.columns:
        return 0.0
    ts = pd.to_datetime(df['timestamp'])
    ok = (df[vcol].astype(float) - df[scol].astype(float)).abs() < 2.0
    if not ok.any():
        return 0.0
    i_ok = int(np.argmax(ok.to_numpy()))
    # a 5 V step settles in 2-3 s, well inside the pedestal run that precedes
    # the data chunks -- only a genuine ramp (first point, trips) needs a cut
    if (ts.iloc[i_ok] - ts.iloc[0]).total_seconds() < 10.0:
        return 0.0
    t_settle = ts.iloc[i_ok]
    m = re.search(r'(\d{6})_(\d{2})H(\d{2})', os.path.basename(chunk0))
    if not m:
        return 0.0
    daq_start = pd.Timestamp(f'20{m.group(1)[:2]}-{m.group(1)[2:4]}-'
                             f'{m.group(1)[4:]} {m.group(2)}:{m.group(3)}:00')
    return max(0.0, (t_settle - daq_start).total_seconds() + margin)


def load_events(cfg, ct, subrun, cluster_r, min_amp, veto_sparks):
    """Per-event Fe55 observables for one sub_run, streamed chunk by chunk.

    Returns (events, pad_hits, live_s, n_veto):
      events   DataFrame [eventId, t [s], q_clus [ADC], n_clus, a_lead [ADC],
               lead_pad, sat]
      pad_hits per-pad TOTAL hit counts (Series indexed by channel_id, after
               the same cuts/veto/mapping) for the hit map
    """
    sub = cfg.subrun_dir(subrun)
    hits_dir = os.path.join(sub, 'combined_hits_root')
    chunk0 = p2io.hit_files(hits_dir)[0]
    t_min = settle_t_min(cfg, sub, chunk0)
    if t_min > 0:
        print(f'    HV-settle cut: dropping events before t = {t_min:.0f} s')

    veto = None
    if veto_sparks:
        hv_csv = os.path.join(sub, 'hv_monitor.csv')
        if os.path.isfile(hv_csv):
            veto = ps.SparkVeto.from_csv(hv_csv, cfg)

    branches = ['eventId', 'channel', 'amplitude', 'feu',
                'trigger_timestamp_ns', 'saturated']
    parts, hit_parts = [], []
    n_veto = 0
    for df in p2io.iter_hits(hits_dir, branches, ct.attrs['feus'],
                             min_amp=min_amp):
        if t_min > 0:
            df = df[df['trigger_timestamp_ns'].astype(np.int64) / 1e9 >= t_min]
        if veto is not None and len(df):
            n0 = df['eventId'].nunique()
            df, _ = veto.apply(df)
            n_veto += n0 - df['eventId'].nunique()
        if not len(df):
            continue
        h = pmap.attach_pads_to_hits(df, ct)
        h = h[h['mapped'] & h['pad_cx'].notna()]
        if cfg.NOISY_PADS:
            h = h[~h['channel_id'].isin(set(cfg.NOISY_PADS))]
        del df
        if not len(h):
            continue
        hit_parts.append(h.groupby('channel_id').size())
        lead = h.loc[h.groupby('eventId')['amplitude'].idxmax(),
                     ['eventId', 'amplitude', 'pad_cx', 'pad_cy', 'channel_id',
                      'trigger_timestamp_ns']].rename(
            columns={'amplitude': 'a_lead', 'pad_cx': 'lx', 'pad_cy': 'ly',
                     'channel_id': 'lead_pad'})
        h = h.merge(lead[['eventId', 'a_lead', 'lx', 'ly']], on='eventId')
        near = ((h['pad_cx'] - h['lx']) ** 2 +
                (h['pad_cy'] - h['ly']) ** 2) <= cluster_r ** 2
        g = h[near].groupby('eventId')
        ev = pd.DataFrame({'q_clus': g['amplitude'].sum(),
                           'n_clus': g.size(),
                           'sat': g['saturated'].any()})
        ev = ev.join(lead.set_index('eventId')[['a_lead', 'lead_pad',
                                                'trigger_timestamp_ns']])
        ev['t'] = ev['trigger_timestamp_ns'].astype(np.int64) / 1e9
        parts.append(ev.drop(columns='trigger_timestamp_ns').reset_index())

    if not parts:
        return pd.DataFrame(), pd.Series(dtype=np.int64), 0.0, n_veto
    pad_hits = (pd.concat(hit_parts, axis=1).fillna(0).sum(axis=1)
                .astype(np.int64) if hit_parts
                else pd.Series(dtype=np.int64))
    events = pd.concat(parts, ignore_index=True)
    # an event straddling a chunk boundary appears twice: keep the bigger half
    events = (events.sort_values('n_clus').drop_duplicates('eventId',
                                                           keep='last')
              .sort_values('t').reset_index(drop=True))
    live_s = float(events['t'].max() - events['t'].min())
    if veto is not None:
        live_s = max(0.0, live_s - veto.vetoed_seconds)
    return events, pad_hits, live_s, n_veto


def _gauss(x, a, mu, sig):
    return a * np.exp(-0.5 * ((x - mu) / sig) ** 2)


def fit_photopeak(q, fit_min, nbins=160):
    """Iterative Gaussian fit of the Fe55 photopeak.

    q       : per-event cluster charges [ADC]
    fit_min : ignore the low-charge noise bulge below this

    Returns dict(mu, mu_err, sigma, ok, edges, counts) -- edges/counts are the
    full-range histogram for plotting; mu is NaN when the fit fails.
    """
    q = np.asarray(q, dtype=float)
    q = q[np.isfinite(q)]
    hi = np.quantile(q, 0.999) if len(q) else 1.0
    counts, edges = np.histogram(q, bins=nbins, range=(0.0, max(hi, fit_min * 2)))
    out = dict(mu=np.nan, mu_err=np.nan, sigma=np.nan, ok=False,
               near_floor=False, edges=edges, counts=counts)
    # fit on a FINE histogram capped at the 98th percentile: the full-range
    # (99.9%) binning above is for plotting only -- its bin width is set by
    # the long cosmic/pile-up tail and can leave the +-30% window around a
    # low-lying photopeak with too few bins to fit.
    hi_fit = np.quantile(q, 0.98) if len(q) else 1.0
    fcounts, fedges = np.histogram(q, bins=200,
                                   range=(0.0, max(hi_fit, fit_min * 2)))
    ctr = 0.5 * (fedges[:-1] + fedges[1:])
    counts = fcounts
    sel = ctr >= fit_min
    if sel.sum() < 8 or counts[sel].sum() < 200:
        return out
    smooth = np.convolve(counts.astype(float), np.ones(5) / 5, mode='same')
    mu = float(ctr[sel][np.argmax(smooth[sel])])
    for _ in range(3):
        win = (ctr > 0.70 * mu) & (ctr < 1.30 * mu) & (ctr >= fit_min)
        if win.sum() < 5:
            return out
        try:
            p, cov = curve_fit(_gauss, ctr[win], counts[win],
                               p0=[counts[win].max(), mu, 0.15 * mu],
                               maxfev=5000)
        except (RuntimeError, ValueError):
            return out
        mu, sig = float(p[1]), abs(float(p[2]))
        if not (fit_min < mu < edges[-1]):
            return out
    out.update(mu=mu, mu_err=float(np.sqrt(max(cov[1][1], 0.0))),
               sigma=sig, ok=True,
               # a "peak" hugging the fit floor is usually the noise bulge
               # (e.g. det3 fits 270 ADC at every HV): keep it in the CSV but
               # exclude it from the gain curve's exponential fit
               near_floor=bool(mu < 1.5 * fit_min))
    return out


def _windowed_median_peak(q, fit_min, n_iter=3):
    """Robust per-pad photopeak: median of q above fit_min, then re-median
    inside a +-40% window around the running estimate. Tracks the photopeak
    even with an escape-peak/noise tail, at far lower statistics than a fit."""
    q = q[q >= fit_min]
    if not len(q):
        return np.nan
    m = float(np.median(q))
    for _ in range(n_iter):
        w = q[(q > 0.6 * m) & (q < 1.4 * m)]
        if not len(w):
            break
        m = float(np.median(w))
    return m


def plot_gain_map(ev, ct, cfg, hv, out_png, out_csv, fit_min, min_events,
                  subrun):
    """Per-pad Fe55 peak (events assigned to their leading pad) on the real
    pad tiles. Grey = live pad with < min_events entries. Returns the per-pad
    table (channel_id, connector_N, n_events, peak_adc)."""
    from matplotlib.collections import PolyCollection

    sig = ev[ev['q_clus'] >= fit_min]
    g = sig.groupby('lead_pad')['q_clus']
    per_pad = pd.DataFrame({
        'n_events': g.size(),
        'peak_adc': g.apply(lambda q: _windowed_median_peak(q.to_numpy(),
                                                            fit_min)),
    })
    per_pad.loc[per_pad['n_events'] < min_events, 'peak_adc'] = np.nan
    padc = (ct.drop_duplicates('channel_id')
            [['channel_id', 'connector_N', 'pad_cx', 'pad_cy']]
            .merge(per_pad, left_on='channel_id', right_index=True,
                   how='left'))
    padc['n_events'] = padc['n_events'].fillna(0).astype(int)
    padc.to_csv(out_csv, index=False)

    ok = padc['peak_adc'].notna()
    spread = (float(padc.loc[ok, 'peak_adc'].std() /
                    padc.loc[ok, 'peak_adc'].mean()) if ok.sum() > 1 else np.nan)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    if pmap.has_tile_geometry(ct):
        pads, verts = pmap.pad_tiles(ct)
        idx = padc.set_index('channel_id').reindex(pads['channel_id'])
        peak = idx['peak_adc'].to_numpy()
        nev = idx['n_events'].fillna(0).to_numpy()
        good = np.isfinite(peak)
        # left: gain map
        ax = axes[0]
        ax.add_collection(PolyCollection(verts[~good], facecolors='0.92',
                                         edgecolors='0.7', linewidths=0.3))
        if good.any():
            pc = PolyCollection(verts[good], array=peak[good], cmap='viridis',
                                edgecolors='face', linewidths=0.2)
            ax.add_collection(pc)
            fig.colorbar(pc, ax=ax, label='Fe55 peak [ADC]')
        ax.set_title(f'per-pad Fe55 peak — grey = < {min_events} events\n'
                     f'{int(good.sum())}/{len(good)} pads, '
                     f'spread std/mean = {spread:.2f}')
        # right: per-pad statistics
        ax = axes[1]
        fired = nev > 0
        ax.add_collection(PolyCollection(verts[~fired], facecolors='0.92',
                                         edgecolors='0.7', linewidths=0.3))
        if fired.any():
            pc = PolyCollection(verts[fired], array=nev[fired], cmap='inferno',
                                norm=matplotlib.colors.LogNorm(vmin=1),
                                edgecolors='face', linewidths=0.2)
            ax.add_collection(pc)
            fig.colorbar(pc, ax=ax, label='leading-pad events (log)')
        ax.set_title('per-pad event count (source illumination)')
        for ax in axes:
            ax.set_xlabel('pad_cx [mm]'); ax.set_ylabel('pad_cy [mm]')
            ax.autoscale_view(); ax.set_aspect('equal')
    else:
        for ax, col, cm in [(axes[0], 'peak_adc', 'viridis'),
                            (axes[1], 'n_events', 'inferno')]:
            s = ax.scatter(padc['pad_cx'], padc['pad_cy'], c=padc[col],
                           cmap=cm, s=18)
            fig.colorbar(s, ax=ax, label=col)
            ax.set_aspect('equal')
    fig.suptitle(f'{cfg.DET_NAME} Fe55 gain map — mesh {hv} V ({subrun})')
    fig.tight_layout()
    fig.savefig(out_png, dpi=170, bbox_inches='tight')
    plt.close(fig)
    return padc, spread


def plot_hit_map(pad_hits, ct, cfg, hv, out_png, out_csv, subrun):
    """Per-pad TOTAL hit counts on the real pad tiles, linear + log panels
    (05_detector_deep_qa surface-hitmap style). Grey = live pad, no hits."""
    from matplotlib.collections import PolyCollection

    padc = (ct.drop_duplicates('channel_id')
            [['channel_id', 'connector_N', 'pad_cx', 'pad_cy']].copy())
    padc['n_hits'] = (pad_hits.reindex(padc['channel_id']).fillna(0)
                      .astype(int).to_numpy())
    padc.to_csv(out_csv, index=False)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    if pmap.has_tile_geometry(ct):
        pads, verts = pmap.pad_tiles(ct)
        counts = (padc.set_index('channel_id')['n_hits']
                  .reindex(pads['channel_id']).fillna(0).to_numpy())
        fired = counts > 0
        for ax, norm, tag in [(axes[0], None, 'linear'),
                              (axes[1],
                               matplotlib.colors.LogNorm(vmin=1), 'log')]:
            ax.add_collection(PolyCollection(verts[~fired], facecolors='0.92',
                                             edgecolors='0.7',
                                             linewidths=0.3))
            if fired.any():
                pc = PolyCollection(verts[fired], array=counts[fired],
                                    cmap='inferno', norm=norm,
                                    edgecolors='face', linewidths=0.2)
                ax.add_collection(pc)
                fig.colorbar(pc, ax=ax, label=f'pad hits ({tag})')
            ax.set_xlabel('pad_cx [mm]'); ax.set_ylabel('pad_cy [mm]')
            ax.autoscale_view(); ax.set_aspect('equal')
            ax.set_title(f'per-pad hits ({tag}) — grey = never fired')
    else:
        for ax, norm, tag in [(axes[0], None, 'linear'),
                              (axes[1],
                               matplotlib.colors.LogNorm(vmin=1), 'log')]:
            s = ax.scatter(padc['pad_cx'], padc['pad_cy'], c=padc['n_hits'],
                           cmap='inferno', norm=norm, s=18)
            fig.colorbar(s, ax=ax, label=f'pad hits ({tag})')
            ax.set_aspect('equal')
    fig.suptitle(f'{cfg.DET_NAME} Fe55 hit map — mesh {hv} V ({subrun}), '
                 f'{int(padc["n_hits"].sum()):,} hits')
    fig.tight_layout()
    fig.savefig(out_png, dpi=170, bbox_inches='tight')
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description='Fe55 spectra + gain vs mesh HV.')
    ap.add_argument('run_key', nargs='?', default='det2_fe55scan1')
    ap.add_argument('--strategy', default='reverse',
                    choices=['linear', 'reverse', 'pairswap'])
    ap.add_argument('--cluster-r', type=float, default=15.0,
                    help='pads within this of the leading pad join the '
                         'cluster charge [mm].')
    ap.add_argument('--min-amp', type=float, default=None,
                    help='per-hit amplitude floor [ADC]; default = run-config '
                         'MIN_AMP.')
    ap.add_argument('--fit-min', type=float, default=200.0,
                    help='photopeak search starts above this cluster charge '
                         '[ADC] (skips the other-detector-trigger noise bulge).')
    ap.add_argument('--veto-sparks', action=argparse.BooleanOptionalAction,
                    default=True, help='per-sub_run HV spark veto (default on).')
    ap.add_argument('--map-min-events', type=int, default=30,
                    help='min leading-pad events for a pad to get a gain-map '
                         'entry.')
    args = ap.parse_args()

    cfg = qa.get_config(args.run_key)
    print(cfg)
    if args.min_amp is None:
        args.min_amp = cfg.MIN_AMP
    out_dir = cfg.out_dir('18_fe55_spectra')
    spec_dir = cfg.out_dir('18_fe55_spectra', 'spectra')
    map_dir = cfg.out_dir('18_fe55_spectra', 'gain_maps')
    hit_dir = cfg.out_dir('18_fe55_spectra', 'hit_maps')
    suffix = cfg.product_suffix(args.veto_sparks)

    ct = pmap.build_channel_table(cfg.run_config_path, cfg.MAP_CSV_PATH,
                                  det_type=cfg.DET_TYPE, det_name=cfg.DET_NAME,
                                  strategy=args.strategy,
                                  drop_connectors=cfg.DEAD_CONNECTORS)
    conns = sorted(ct['connector_N'].unique())
    print(f'wired connectors (run_config dream_feus): {conns}, '
          f'{ct["channel_id"].nunique()} pads on FEUs {ct.attrs["feus"]}')

    subruns = find_subruns(cfg)
    if not subruns:
        print(f'No fe55_* sub_runs with hits under {cfg.run_dir}')
        return
    print(f'Fe55 scan points ({len(subruns)}): '
          + ', '.join(f'{hv}V' for _, hv in subruns))

    rows, spectra = [], []
    for subrun, hv in subruns:
        print(f'  mesh {hv}V ({subrun}):')
        ev, pad_hits, live_s, n_veto = load_events(
            cfg, ct, subrun, args.cluster_r, args.min_amp, args.veto_sparks)
        if not len(ev):
            print('    no events after cuts, skipping')
            continue
        plot_hit_map(pad_hits, ct, cfg, hv,
                     os.path.join(hit_dir, f'hit_map_{hv}V{suffix}.png'),
                     os.path.join(hit_dir, f'hit_map_{hv}V{suffix}.csv'),
                     subrun)
        fit = fit_photopeak(ev['q_clus'], args.fit_min)
        rate = len(ev) / live_s if live_s > 0 else np.nan
        res = 2.355 * fit['sigma'] / fit['mu'] if fit['ok'] else np.nan
        _, map_spread = plot_gain_map(
            ev, ct, cfg, hv,
            os.path.join(map_dir, f'gain_map_{hv}V{suffix}.png'),
            os.path.join(map_dir, f'gain_map_{hv}V{suffix}.csv'),
            args.fit_min, args.map_min_events, subrun)
        rows.append(dict(hv=hv, subrun=subrun, n_events=len(ev),
                         live_s=live_s, rate_hz=rate, n_spark_veto=n_veto,
                         peak_adc=fit['mu'], peak_err=fit['mu_err'],
                         sigma_adc=fit['sigma'], fwhm_over_peak=res,
                         sat_frac=float(ev['sat'].mean()),
                         median_nclus=float(ev['n_clus'].median()),
                         gain_map_spread=map_spread,
                         near_floor=fit['near_floor']))
        spectra.append(dict(hv=hv, fit=fit, q=ev['q_clus'].to_numpy()))
        print(f'    {len(ev):,} events, {rate:.0f} Hz, peak = '
              + (f'{fit["mu"]:.0f} +- {fit["mu_err"]:.0f} ADC, '
                 f'FWHM/peak = {res:.2f}' if fit['ok'] else 'FIT FAILED')
              + f', sat {ev["sat"].mean():.1%}'
              + (f', spark-vetoed {n_veto}' if n_veto else ''))

        # per-point spectrum
        fig, ax = plt.subplots(figsize=(7, 4.5))
        ctr = 0.5 * (fit['edges'][:-1] + fit['edges'][1:])
        ax.step(ctr, fit['counts'], where='mid', color='steelblue', lw=1.2,
                label=f'cluster charge (N={len(ev):,})')
        if fit['ok']:
            xs = np.linspace(0.6 * fit['mu'], 1.4 * fit['mu'], 200)
            amp = fit['counts'][np.argmin(np.abs(ctr - fit['mu']))]
            ax.plot(xs, _gauss(xs, amp, fit['mu'], fit['sigma']), '-',
                    color='crimson', lw=2,
                    label=(f'photopeak {fit["mu"]:.0f} ADC, '
                           f'FWHM/peak {res:.2f}'))
        ax.axvline(args.fit_min, color='gray', ls=':', lw=1,
                   label=f'fit-min {args.fit_min:g}')
        ax.set_xlabel('cluster charge [ADC]'); ax.set_ylabel('events / bin')
        ax.set_yscale('log'); ax.grid(True, alpha=0.3); ax.legend(fontsize=8)
        ax.set_title(f'{cfg.DET_NAME} Fe55 spectrum — mesh {hv} V ({subrun})')
        fig.tight_layout()
        fig.savefig(os.path.join(spec_dir, f'spectrum_{hv}V{suffix}.png'),
                    dpi=170, bbox_inches='tight'); plt.close(fig)

    if not rows:
        print('No usable HV points.')
        return
    df = pd.DataFrame(rows).sort_values('hv').reset_index(drop=True)
    df.to_csv(os.path.join(out_dir, f'fe55_gain_vs_hv{suffix}.csv'),
              index=False)

    tag = f'{cfg.DET_TAG} {cfg.RUN}'
    cmap = plt.get_cmap('viridis')

    # spectra overlay
    fig, ax = plt.subplots(figsize=(8, 5))
    hvs = [s['hv'] for s in spectra]
    for s in spectra:
        f = s['fit']
        ctr = 0.5 * (f['edges'][:-1] + f['edges'][1:])
        c = cmap((s['hv'] - min(hvs)) / max(1, max(hvs) - min(hvs)))
        ax.step(ctr, f['counts'], where='mid', color=c, lw=1.2,
                label=f'{s["hv"]} V')
    ax.set_xlabel('cluster charge [ADC]'); ax.set_ylabel('events / bin')
    ax.set_yscale('log'); ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, ncol=2, title='mesh HV')
    ax.set_title(f'{cfg.DET_NAME} Fe55 spectra vs mesh HV — {tag}')
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, f'fe55_spectra_overlay{suffix}.png'),
                dpi=200, bbox_inches='tight'); plt.close(fig)

    # gain vs HV (+ exponential fit when >= 3 good, non-floor points)
    ok = df[df['peak_adc'].notna() & ~df['near_floor']]
    floor = df[df['peak_adc'].notna() & df['near_floor']]
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.errorbar(ok['hv'], ok['peak_adc'], yerr=ok['peak_err'], fmt='o-',
                color='steelblue', capsize=4, lw=2, ms=7, label='Fe55 photopeak')
    if len(floor):
        ax.plot(floor['hv'], floor['peak_adc'], 'o', mfc='none',
                color='gray', ms=7, label='near fit floor (noise bulge?)')
    if len(ok) >= 3:
        k, b = np.polyfit(ok['hv'], np.log(ok['peak_adc']), 1)
        xs = np.linspace(ok['hv'].min(), ok['hv'].max(), 100)
        ax.plot(xs, np.exp(b + k * xs), '--', color='crimson', lw=1.5,
                label=f'exp fit: x2 every {np.log(2) / k:.1f} V')
    ax.set_yscale('log')
    ax.set_xlabel('mesh HV [V]'); ax.set_ylabel('Fe55 peak [ADC]')
    ax.grid(True, alpha=0.3, which='both'); ax.legend()
    ax.set_title(f'{cfg.DET_NAME} gain (Fe55 peak) vs mesh HV — {tag}')
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, f'fe55_gain_vs_hv{suffix}.png'),
                dpi=200, bbox_inches='tight'); plt.close(fig)

    # resolution vs HV
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(ok['hv'], ok['fwhm_over_peak'], 'o-', color='steelblue', lw=2, ms=7)
    ax.set_xlabel('mesh HV [V]'); ax.set_ylabel('FWHM / peak')
    ax.set_ylim(0, None); ax.grid(True, alpha=0.3)
    ax.set_title(f'{cfg.DET_NAME} Fe55 energy resolution vs mesh HV — {tag}')
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, f'fe55_resolution_vs_hv{suffix}.png'),
                dpi=200, bbox_inches='tight'); plt.close(fig)

    # rate + saturation vs HV
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(df['hv'], df['rate_hz'], 'o-', color='steelblue', lw=2, ms=7,
            label='trigger rate')
    ax.set_xlabel('mesh HV [V]'); ax.set_ylabel('event rate [Hz]',
                                                color='steelblue')
    ax.grid(True, alpha=0.3)
    ax2 = ax.twinx()
    ax2.plot(df['hv'], 100 * df['sat_frac'], 's--', color='darkorange', ms=6,
             label='saturated events')
    ax2.set_ylabel('saturated events [%]', color='darkorange')
    ax2.set_ylim(0, None)
    ax.set_title(f'{cfg.DET_NAME} rate / saturation vs mesh HV — {tag}')
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, f'fe55_rate_vs_hv{suffix}.png'),
                dpi=200, bbox_inches='tight'); plt.close(fig)

    print(f'\n{"HV[V]":>6}  {"peak":>7}  {"+-":>5}  {"FWHM/pk":>8}  '
          f'{"rate":>7}  {"sat%":>5}  {"events":>7}')
    for _, r in df.iterrows():
        print(f'{r.hv:>6.0f}  {r.peak_adc:>7.0f}  {r.peak_err:>5.0f}  '
              f'{r.fwhm_over_peak:>8.2f}  {r.rate_hz:>7.0f}  '
              f'{100 * r.sat_frac:>5.1f}  {r.n_events:>7.0f}')
    print(f'\nWritten to: {out_dir}')


if __name__ == '__main__':
    main()
