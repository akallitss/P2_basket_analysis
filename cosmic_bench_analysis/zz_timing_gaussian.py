#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zz_timing_gaussian.py  (ad-hoc, not part of the numbered pipeline)

Show the underlying per-event timestamp distribution that produces the
`time_res_ns` value plotted by 16_drift_scan_efficiency.py, and compare the
two per-event estimators.

For every drift point we reproduce EXACTLY the selection used in
p2_io.event_time_resolution (signal-band pad hits, amplitude >= sig_amp, the
same min_amp / dropped-pad / spark-veto cuts as the scan) but in ONE streaming
pass we build BOTH per-event estimators and keep the arrays:
    'leading'  = time_of_max of the highest-amplitude pad   (pipeline default)
    'earliest' = earliest time_of_max among signal pads      (leading edge)
sigma_t = (p84.135 - p15.865)/2 (robust Gaussian core) -- the exact CSV number.

Products (-> 16_drift_scan_efficiency/):
  time_resolution_distribution_615V_gaussian.png   (3 points, leading)
  time_resolution_distribution_grid_all.png        (all 12 points, leading)
  time_resolution_estimator_comparison.png         (leading vs earliest)
"""

import os
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import p2_qa_config as qa
import p2_mapping as pmap
import p2_sparks as ps
import p2_io as p2io

MESH = 415.0
GAP_CM = 0.3
DRIFTS = [415, 465, 515, 565, 615, 665, 715, 765, 815, 865, 915, 965]


def per_event_times(cfg, ct, subrun, sig_amp=300.0, veto_sparks=True):
    """One streaming pass -> dict(leading=array, earliest=array) of per-event
    time_of_max [ns], same selection as p2_io.event_time_resolution."""
    sub = cfg.subrun_dir(subrun)
    hits_dir = os.path.join(sub, 'combined_hits_root')
    drop = set(p2io.drop_pads_for(cfg, ct, hits_dir=hits_dir))
    bad = None
    if veto_sparks:
        hv_csv = os.path.join(sub, 'hv_monitor.csv')
        if os.path.isfile(hv_csv):
            sv = ps.SparkVeto.from_csv(hv_csv, cfg)
            bad = sv.vetoed_ids_from_hits(hits_dir, ct.attrs['feus'],
                                          min_amp=cfg.MIN_AMP)
    excl = set(bad) if bad is not None else None
    thr = max(float(cfg.MIN_AMP), float(sig_amp))
    lead_a, lead_t, early_t = {}, {}, {}
    for df in p2io.iter_hits(hits_dir,
                             ['eventId', 'channel', 'amplitude', 'feu',
                              'time_of_max'],
                             ct.attrs['feus'], progress=False,
                             min_amp=cfg.MIN_AMP):
        h = pmap.attach_pads_to_hits(df, ct)
        h = h[h['mapped'] & h['pad_cx'].notna() & (h['amplitude'] >= thr)]
        if drop:
            h = h[~h['channel_id'].isin(drop)]
        if excl is not None and len(h):
            h = h[~h['eventId'].isin(excl)]
        del df
        if not len(h):
            continue
        # leading = max-amplitude pad
        idx = h.groupby('eventId')['amplitude'].idxmax()
        for ev, a, t in h.loc[idx, ['eventId', 'amplitude',
                                    'time_of_max']].itertuples(index=False):
            if ev not in lead_a or a > lead_a[ev]:
                lead_a[ev] = float(a)
                lead_t[ev] = float(t)
        # earliest = min time_of_max
        g = h.groupby('eventId')['time_of_max'].min()
        for ev, t in g.items():
            if ev not in early_t or t < early_t[ev]:
                early_t[ev] = float(t)
    return dict(leading=np.fromiter(lead_t.values(), dtype=np.float64),
                earliest=np.fromiter(early_t.values(), dtype=np.float64))


def robust(t):
    if len(t) < 20:
        return np.nan, np.nan, np.nan, np.nan
    p16, p50, p84 = np.percentile(t, [15.865, 50.0, 84.135])
    return p50, 0.5 * (p84 - p16), p16, p84


def _hist_gauss(ax, t, color, label_data, nsig=4.5, annotate=True):
    p50, sig, p16, p84 = robust(t)
    if not np.isfinite(sig) or sig <= 0:
        ax.text(0.5, 0.5, f'{label_data}\nN={len(t)} (too few)',
                transform=ax.transAxes, ha='center', va='center', fontsize=9)
        return p50, sig
    lo, hi = p50 - nsig * sig, p50 + nsig * sig
    sel = t[(t >= lo) & (t <= hi)]
    bw = max(2.0, sig / 6.0)
    bins = np.arange(lo, hi + bw, bw)
    ax.hist(sel, bins=bins, color=color, alpha=0.35, density=True,
            label=f'{label_data} (N={len(t):,})')
    xs = np.linspace(lo, hi, 400)
    g = np.exp(-0.5 * ((xs - p50) / sig) ** 2) / (sig * np.sqrt(2 * np.pi))
    ax.plot(xs, g, '-', color=color, lw=2.0,
            label=fr'Gauss $\sigma_t$={sig:.1f} ns')
    if annotate:
        ax.axvline(p16, color=color, ls=':', lw=1.2, alpha=0.9)
        ax.axvline(p84, color=color, ls=':', lw=1.2, alpha=0.9)
    ax.grid(True, alpha=0.3)
    return p50, sig


def fig_three(cfg, cache, out_dir):
    """Original 3-point figure (turn-on / best / high field), leading."""
    pts = [(465, 'turn-on', 'darkorange'),
           (615, r'BEST (min $\sigma_t$)', 'indigo'),
           (965, 'high field', 'crimson')]
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.2))
    for (dv, lbl, col), ax in zip(pts, axes):
        t = cache[dv]['leading']
        p50, sig = _hist_gauss(ax, t, col, 'data')
        e = (dv - MESH) / GAP_CM
        ax.set_title(f'{lbl}: drift {dv} V ({e:.0f} V/cm) '
                     fr'$\to\ \sigma_t$={sig:.1f} ns', fontsize=10.5)
        ax.set_xlabel('leading-pad time_of_max [ns] (ref. trigger)')
        ax.set_ylabel('probability density')
        ax.legend(fontsize=8.5, loc='upper right')
    fig.suptitle(f'{cfg.DET_NAME} — per-event timestamp distribution behind the '
                 f'time-resolution estimate  ({cfg.RUN})\n'
                 r'leading-pad time_of_max, $\sigma_t=(p_{84.1}-p_{15.9})/2$',
                 fontsize=11.5)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    p = os.path.join(out_dir, 'time_resolution_distribution_615V_gaussian.png')
    fig.savefig(p, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print('Written:', p)


def fig_grid(cfg, cache, out_dir):
    """4x3 grid: every drift point, leading estimator, Gaussian overlay."""
    fig, axes = plt.subplots(4, 3, figsize=(15, 16))
    for dv, ax in zip(DRIFTS, axes.ravel()):
        t = cache[dv]['leading']
        col = 'indigo' if dv == 615 else 'steelblue'
        p50, sig = _hist_gauss(ax, t, col, 'data')
        e = (dv - MESH) / GAP_CM
        star = '  ★' if dv == 615 else ''
        ax.set_title(f'drift {dv} V  ({e:.0f} V/cm)  '
                     fr'$\sigma_t$={sig:.1f} ns{star}', fontsize=10)
        ax.set_xlabel('leading-pad time_of_max [ns]', fontsize=8)
        ax.legend(fontsize=7.5, loc='upper right')
    fig.suptitle(f'{cfg.DET_NAME} — per-event timestamp distribution, ALL drift '
                 f'points  ({cfg.RUN})\n'
                 'leading-pad time_of_max ref. trigger; mesh 415 V fixed, 3 mm '
                 r'gap; $\sigma_t=(p_{84.1}-p_{15.9})/2$.  415 V = near-zero '
                 'field (few events, huge spread); 765 V timing normal though '
                 'its tracking failed', fontsize=11.5)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    p = os.path.join(out_dir, 'time_resolution_distribution_grid_all.png')
    fig.savefig(p, dpi=170, bbox_inches='tight')
    plt.close(fig)
    print('Written:', p)


def fig_estimator(cfg, cache, out_dir):
    """Leading-amplitude vs earliest-hit: sigma_t(drift) + core overlay @615V."""
    e_ok = [d for d in DRIFTS if (d - MESH) / GAP_CM > 50]
    rows = []
    for dv in e_ok:
        _, sl, _, _ = robust(cache[dv]['leading'])
        _, se, _, _ = robust(cache[dv]['earliest'])
        rows.append((dv, (dv - MESH) / GAP_CM, sl, se))
    x = [r[0] for r in rows]
    sl = [r[2] for r in rows]
    se = [r[3] for r in rows]

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(15, 5.6))
    # left: sigma_t vs drift, both estimators
    axL.plot(x, sl, 'o-', color='indigo', lw=2, ms=8,
             label='leading pad (max amplitude) — pipeline')
    axL.plot(x, se, 's--', color='teal', lw=2, ms=7,
             label='earliest pad (leading edge)')
    for xi, a, b in zip(x, sl, se):
        axL.annotate(f'{b-a:+.1f}', xy=(xi, min(a, b)), xytext=(0, -13),
                     textcoords='offset points', ha='center', fontsize=7,
                     color='teal')
    iL = int(np.argmin(sl))
    iE = int(np.argmin(se))
    axL.annotate(f'min {sl[iL]:.1f} ns @ {x[iL]} V', xy=(x[iL], sl[iL]),
                 xytext=(x[iL], sl[iL] + 4), fontsize=8, color='indigo',
                 ha='center', arrowprops=dict(arrowstyle='->', color='indigo'))
    axL.set_xlabel('drift HV [V]')
    axL.set_ylabel(r'estimated time resolution $\sigma_t$ [ns]')
    axL.set_ylim(0, None)
    axL.grid(True, alpha=0.3)
    axL.legend(fontsize=9, loc='upper left')
    axL.set_title('Estimator comparison across the drift scan\n'
                  '(teal numbers = earliest − leading, ns)', fontsize=10.5)
    # right: the two cores at the best point
    for est, col in (('leading', 'indigo'), ('earliest', 'teal')):
        t = cache[615][est]
        t = t - np.median(t)          # centre both at 0 to compare widths
        _hist_gauss(axR, t, col, est, annotate=False)
    axR.set_xlabel('time_of_max − median [ns]  (ref. trigger)')
    axR.set_ylabel('probability density')
    axR.legend(fontsize=9, loc='upper right')
    axR.set_title('Core shape at the best point (drift 615 V, 667 V/cm)\n'
                  'both centred at 0 to compare widths', fontsize=10.5)
    fig.suptitle(f'{cfg.DET_NAME} — leading-amplitude vs earliest-hit timing '
                 f'estimator  ({cfg.RUN})', fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    p = os.path.join(out_dir, 'time_resolution_estimator_comparison.png')
    fig.savefig(p, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print('Written:', p)
    print(f'\n{"drift":>6} {"E[V/cm]":>8} {"lead":>7} {"early":>7} {"diff":>6}')
    for dv, e, a, b in rows:
        print(f'{dv:>6} {e:>8.0f} {a:>7.1f} {b:>7.1f} {b-a:>+6.1f}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run-key', default='det1_driftscan2')
    ap.add_argument('--sig-amp', type=float, default=300.0)
    args = ap.parse_args()
    cfg = qa.get_config(args.run_key)
    ct = pmap.build_channel_table(cfg.run_config_path, cfg.MAP_CSV_PATH,
                                  det_type=cfg.DET_TYPE, det_name=cfg.DET_NAME,
                                  strategy='reverse',
                                  drop_connectors=cfg.DEAD_CONNECTORS)
    out_dir = cfg.out_dir('16_drift_scan_efficiency')

    cache = {}
    for dv in DRIFTS:
        subrun = f'drift_scan_det1_{int(MESH)}_{dv}'
        cache[dv] = per_event_times(cfg, ct, subrun, sig_amp=args.sig_amp)
        pl = robust(cache[dv]['leading'])[1]
        pe = robust(cache[dv]['earliest'])[1]
        print(f'drift {dv:>3} V  E={(dv-MESH)/GAP_CM:>5.0f} V/cm  '
              f'N_lead={len(cache[dv]["leading"]):>4}  '
              f'sig_lead={pl:>6.2f}  sig_early={pe:>6.2f} ns')

    fig_three(cfg, cache, out_dir)
    fig_grid(cfg, cache, out_dir)
    fig_estimator(cfg, cache, out_dir)


if __name__ == '__main__':
    main()
