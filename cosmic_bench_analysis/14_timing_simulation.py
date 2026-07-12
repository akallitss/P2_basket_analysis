#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
14_timing_simulation.py

Physics limit of the P2 (Micromegas, Ar/iC4H10 95:5) timing resolution vs the
applied drift field, for a 3 mm (+- 0.5 mm) conversion gap — is the ~20 ns
goal reachable, and are we algorithm-limited or physics-limited?

Inputs — computed with Garfield++ on lxplus (p2_gas_heed.cpp):
  gas_table.csv      Magboltz: drift velocity + longitudinal diffusion vs E
                     for Ar/iC4H10 95:5 (vz in cm/ns, dl in cm/sqrt(cm))
  heed_clusters.csv  HEED primary-ionisation clusters for 4 GeV cosmic muons
                     crossing a 3 mm slab at zenith 0..50 deg (z = depth above
                     the mesh [cm], ne = electrons in the cluster)

Monte Carlo (per event)
-----------------------
  1. Draw a zenith angle (cos^2 weighting over the HEED angle grid, or
     --vertical) and a HEED track (cluster depths + sizes) at that angle.
  2. Gap tolerance: for gap < 3 mm clusters beyond the gap are dropped; for
     gap > 3 mm the track is extended with a second HEED track shifted by
     3 mm (cluster statistics stay HEED-exact per unit length).
  3. Each electron drifts to the mesh: t = z/v_d(E) + N(0, sigma_L(E)
     sqrt(z)/v_d), gain ~ exponential (conservative avalanche fluctuation).

Timing estimators on the (gain-weighted) arrival series
-------------------------------------------------------
  first_e        arrival of the first electron (ideal threshold limit)
  q05/q10/q20/q50  time at 5/10/20/50% of the total avalanche charge
                 (ideal charge-fraction discriminator = CFD physics limit)
  dream_thr      REALISTIC chain: arrivals convolved with the measured DREAM
                 pulse shape (rise ~140 ns), sampled at 60 ns with measured
                 noise, 5-sigma threshold + ftst-style phase correction —
                 connects the simulation to the stage-13 data numbers.

Products (<Analysis>/<det>/<run>/<sub_run>/14_timing_sim/):
  sim_gas_transport.png        Magboltz v_d(E), sigma_L(E)
  sim_resolution_vs_field.png  sigma_t vs drift field per estimator, gap band
  sim_working_point.png        arrival-time anatomy at the current field
  sim_summary.txt / sim_resolution_vs_field.csv

Usage:
  python3 14_timing_simulation.py [run_key] [--n-events 4000] [--vertical]
          [--gas-table gas_table.csv] [--heed heed_clusters.csv]
"""

import os
import argparse
import matplotlib
matplotlib.use('Agg')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import p2_qa_config as qa

RNG = np.random.RandomState(12345)
# Garfield++ outputs shipped with the repo (generated on lxplus by
# garfield_inputs/p2_gas_heed.cpp; pandas reads the .gz transparently).
GARFIELD_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            'garfield_inputs')


# --------------------------------------------------------------------------- #
# Garfield inputs
# --------------------------------------------------------------------------- #
def load_gas_table(path):
    """Magboltz scan -> (E [V/cm], vd [cm/us], dl [um/sqrt(cm)]).
    Garfield native units: vz in cm/ns, dl in cm/sqrt(cm)."""
    g = pd.read_csv(path).sort_values('E_Vcm')
    E = g['E_Vcm'].to_numpy(float)
    vd = g['vz'].to_numpy(float) * 1e3          # cm/ns -> cm/us
    dl = g['dl'].to_numpy(float) * 1e4          # cm/sqrt(cm) -> um/sqrt(cm)
    return E, vd, dl


def load_heed_tracks(path):
    """HEED cluster dump -> {angle_deg: list of (z_cm array, ne array)}."""
    d = pd.read_csv(path)
    lib = {}
    for (ang, trk), g in d.groupby(['angle_deg', 'track']):
        lib.setdefault(float(ang), []).append(
            (g['z_cm'].to_numpy(float), g['ne'].to_numpy(int)))
    return lib


def sample_track(lib, angle, gap_cm):
    """One HEED track at `angle`, adapted to gap_cm (HEED slab was 0.3 cm)."""
    tracks = lib[angle]
    z, ne = tracks[RNG.randint(len(tracks))]
    if gap_cm <= 0.3:
        m = z < gap_cm
        z, ne = z[m], ne[m]
    else:
        z2, ne2 = tracks[RNG.randint(len(tracks))]
        m = z2 < (gap_cm - 0.3)
        z = np.concatenate([z, z2[m] + 0.3])
        ne = np.concatenate([ne, ne2[m]])
    return z, ne


# --------------------------------------------------------------------------- #
# event -> electron arrival series -> estimators
# --------------------------------------------------------------------------- #
def arrival_series(z, ne, vd_cm_us, dl_um_scm):
    """Drift the cluster electrons to the mesh: sorted times [ns] + gains."""
    zz = np.repeat(z, ne)
    if not len(zz):
        return None, None
    vd = vd_cm_us * 1e-3                       # cm/ns
    sig = (dl_um_scm * 1e-4) * np.sqrt(zz) / vd
    t = zz / vd + RNG.normal(0, 1, len(zz)) * sig
    gain = RNG.exponential(1.0, len(zz))
    order = np.argsort(t)
    return t[order], gain[order]


def charge_fraction_time(t, w, frac):
    c = np.cumsum(w)
    return float(np.interp(frac * c[-1], c, t))


def dream_kernel(dt_ns=5.0, rise_ns=140.0, fall_ns=350.0):
    t = np.arange(0, 2000, dt_ns)
    k = (1 - np.exp(-t / rise_ns)) * np.exp(-t / fall_ns)
    return k / k.max()


_KER = dream_kernel()


def dream_threshold_time(t, w, amp_adc=300.0, noise_adc=5.6, tps=60.0,
                         n_sig=5.0, ftst_ns=10.0):
    """Realistic chain: 5 ns binning, DREAM shaping, 60 ns sampling with a
    random (ftst-quantised) clock phase, 5-sigma threshold + interpolation,
    phase undone (the ftst correction)."""
    dt, tmax = 5.0, 2000.0
    edges = np.arange(0, tmax + dt, dt)
    h, _ = np.histogram(t, bins=edges, weights=w)
    sig = np.convolve(h, _KER)[:len(edges) - 1]
    if sig.max() <= 0:
        return np.nan
    sig = sig / sig.max() * amp_adc
    phase = RNG.randint(0, int(tps / ftst_ns)) * ftst_ns
    samp_t = np.arange(phase, tmax, tps)
    idx = np.clip((samp_t / dt).astype(int), 0, len(sig) - 1)
    wf = sig[idx] + RNG.normal(0, noise_adc, len(idx))
    thr = n_sig * noise_adc
    above = np.nonzero(wf > thr)[0]
    if not len(above) or above[0] == 0:
        return np.nan
    i = above[0]
    y0, y1 = wf[i - 1], wf[i]
    return samp_t[i - 1] + (thr - y0) / (y1 - y0) * tps - phase


def robust_sigma(x, clip=2.5, iters=4):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if len(x) < 10:
        return np.nan
    mu, sig = np.median(x), 1.4826 * np.median(np.abs(x - np.median(x)))
    for _ in range(iters):
        m = np.abs(x - mu) < clip * max(sig, 1e-3)
        if m.sum() < 10:
            break
        mu, sig = x[m].mean(), x[m].std()
    return float(sig)


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description='P2 timing physics simulation.')
    ap.add_argument('run_key', nargs='?', default='det2_long1')
    ap.add_argument('--gas-table',
                    default=os.path.join(GARFIELD_DIR, 'gas_table.csv'))
    ap.add_argument('--heed',
                    default=os.path.join(GARFIELD_DIR, 'heed_clusters.csv.gz'))
    ap.add_argument('--n-events', type=int, default=4000)
    ap.add_argument('--gap-mm', type=float, default=3.0)
    ap.add_argument('--gap-tol-mm', type=float, default=0.5)
    ap.add_argument('--vertical', action='store_true',
                    help='vertical tracks only (default: cos^2 cosmics).')
    args = ap.parse_args()

    cfg = qa.get_config(args.run_key)
    out_dir = cfg.out_dir('14_timing_sim')
    print(cfg)

    E_tab, vd_tab, dl_tab = load_gas_table(args.gas_table)
    lib = load_heed_tracks(args.heed)
    angles = sorted(lib.keys())
    print(f'  Magboltz table: {len(E_tab)} points, E {E_tab.min():.0f}-'
          f'{E_tab.max():.0f} V/cm')
    print(f'  HEED library : angles {angles}, '
          f'{sum(len(v) for v in lib.values()):,} tracks, '
          f'{np.mean([len(t[0]) for t in lib[0.0]]):.1f} clusters/track '
          f'(vertical, 3 mm)')

    def vd_of(E):
        return np.interp(E, E_tab, vd_tab)

    def dl_of(E):
        return np.interp(E, E_tab, dl_tab)

    v_gap = 170.0                                # V across the gap
    e_now = v_gap / (args.gap_mm / 10.0)         # V/cm
    print(f'  working point: {v_gap:.0f} V / {args.gap_mm:g} mm = '
          f'{e_now:.0f} V/cm (vd {vd_of(e_now):.2f} cm/us, '
          f'DL {dl_of(e_now):.0f} um/sqrt(cm))')

    # cos^2-weighted angle sampling probabilities on the HEED grid
    if args.vertical:
        ang_p = np.array([1.0] + [0.0] * (len(angles) - 1))
    else:
        w = np.array([np.cos(np.radians(a)) ** 2 * np.sin(np.radians(max(a, 2)))
                      for a in angles])
        ang_p = w / w.sum()

    e_scan = E_tab.copy()
    gaps_cm = [args.gap_mm / 10.0,
               (args.gap_mm - args.gap_tol_mm) / 10.0,
               (args.gap_mm + args.gap_tol_mm) / 10.0]
    estimators = ['first_e', 'q05', 'q10', 'q20', 'q50', 'dream_thr']
    fracs = {'q05': 0.05, 'q10': 0.10, 'q20': 0.20, 'q50': 0.50}

    def run_point(E, gap_cm, n_ev):
        res = {k: [] for k in estimators}
        for _ in range(n_ev):
            ang = angles[RNG.choice(len(angles), p=ang_p)]
            z, ne = sample_track(lib, ang, gap_cm)
            t, w = arrival_series(z, ne, vd_of(E), dl_of(E))
            if t is None:
                continue
            res['first_e'].append(t[0])
            for k, f in fracs.items():
                res[k].append(charge_fraction_time(t, w, f))
            res['dream_thr'].append(dream_threshold_time(t, w))
        return {k: robust_sigma(res[k]) for k in estimators}

    rows = []
    for E in e_scan:
        for gi, g in enumerate(gaps_cm):
            n_ev = args.n_events if gi == 0 else args.n_events // 2
            r = run_point(E, g, n_ev)
            rows.append(dict(E=E, gap_mm=g * 10, **r))
            if gi == 0:
                print(f'  E {E:6.0f} V/cm: ' +
                      '  '.join(f'{k} {r[k]:5.1f}' for k in estimators) + ' ns')
    df = pd.DataFrame(rows)
    df.to_csv(f'{out_dir}/sim_resolution_vs_field.csv', index=False)

    labels = {'first_e': 'first electron (ideal)',
              'q05': '5% charge', 'q10': '10% charge', 'q20': '20% charge',
              'q50': '50% charge (CFD limit)',
              'dream_thr': 'DREAM 60 ns chain, 5$\\sigma$ thr (realistic)'}
    colors = {'first_e': 'tab:green', 'q05': 'tab:olive', 'q10': 'tab:blue',
              'q20': 'tab:cyan', 'q50': 'tab:purple', 'dream_thr': 'tab:red'}

    # gas transport
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.5))
    ee = np.linspace(E_tab.min(), E_tab.max(), 400)
    axs[0].plot(ee, vd_of(ee), color='steelblue', lw=2)
    axs[0].plot(E_tab, vd_tab, 'o', color='steelblue', ms=4)
    axs[0].axvline(e_now, color='crimson', ls=':',
                   label=f'working point {e_now:.0f} V/cm')
    axs[0].set_xlabel('drift field [V/cm]')
    axs[0].set_ylabel('v$_d$ [cm/$\\mu$s]')
    axs[0].set_title('drift velocity — Ar/iC$_4$H$_{10}$ 95:5 (Magboltz)')
    axs[0].grid(True, alpha=0.3); axs[0].legend()
    axs[1].plot(ee, dl_of(ee), color='seagreen', lw=2)
    axs[1].plot(E_tab, dl_tab, 'o', color='seagreen', ms=4)
    axs[1].axvline(e_now, color='crimson', ls=':')
    axs[1].set_xlabel('drift field [V/cm]')
    axs[1].set_ylabel('$\\sigma_L$ [$\\mu$m/$\\sqrt{cm}$]')
    axs[1].set_title('longitudinal diffusion (Magboltz)')
    axs[1].grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(f'{out_dir}/sim_gas_transport.png', dpi=150, bbox_inches='tight')
    plt.close(fig)

    # resolution vs field
    fig, ax = plt.subplots(figsize=(9.5, 6))
    nom = df[df['gap_mm'] == args.gap_mm]
    lo = df[df['gap_mm'] == args.gap_mm - args.gap_tol_mm]
    hi = df[df['gap_mm'] == args.gap_mm + args.gap_tol_mm]
    for k in estimators:
        ax.plot(nom['E'], nom[k], 'o-', color=colors[k], lw=2, label=labels[k])
        ax.fill_between(nom['E'],
                        np.interp(nom['E'], lo['E'], lo[k]),
                        np.interp(nom['E'], hi['E'], hi[k]),
                        color=colors[k], alpha=0.12)
    ax.axhline(20, color='k', ls='--', lw=1.2, label='20 ns goal')
    ax.axvline(e_now, color='crimson', ls=':', lw=1.5,
               label=f'current: {v_gap:.0f} V / {args.gap_mm:g} mm = {e_now:.0f} V/cm')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel('drift field [V/cm]')
    ax.set_ylabel('timing resolution $\\sigma_t$ [ns]')
    tracks = 'vertical tracks' if args.vertical else 'cos$^2$ cosmic zenith'
    ax.set_title(f'Simulated timing resolution vs drift field — Ar/iso 95:5 '
                 f'(Magboltz+HEED), gap {args.gap_mm:g}$\\pm${args.gap_tol_mm:g} mm '
                 f'(band), {tracks}')
    ax.grid(True, alpha=0.3, which='both'); ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(f'{out_dir}/sim_resolution_vs_field.png', dpi=150,
                bbox_inches='tight')
    plt.close(fig)

    # anatomy at the working point (vertical tracks)
    t_first, t_q10, t_q50 = [], [], []
    for _ in range(3000):
        z, ne = sample_track(lib, 0.0, args.gap_mm / 10.0)
        t, w = arrival_series(z, ne, vd_of(e_now), dl_of(e_now))
        if t is None:
            continue
        t_first.append(t[0])
        t_q10.append(charge_fraction_time(t, w, 0.1))
        t_q50.append(charge_fraction_time(t, w, 0.5))
    fig, ax = plt.subplots(figsize=(9, 5))
    tmaxg = args.gap_mm / 10 / vd_of(e_now) * 1e3
    bins = np.arange(0, max(80, tmaxg * 1.3), 1.0)
    for x, lab, c in [(t_first, 'first electron', 'tab:green'),
                      (t_q10, '10% charge', 'tab:blue'),
                      (t_q50, '50% charge', 'tab:purple')]:
        x = np.asarray(x)
        ax.hist(x, bins=bins, histtype='step', lw=2, color=c,
                label=f'{lab}  ($\\sigma$ {robust_sigma(x):.1f} ns)')
    ax.set_xlabel('arrival time after crossing [ns]'); ax.set_ylabel('events')
    ax.set_title(f'Arrival-time anatomy at the working point ({e_now:.0f} V/cm, '
                 f'vertical)\nv$_d$ {vd_of(e_now):.2f} cm/$\\mu$s — full-gap '
                 f'drift {tmaxg:.0f} ns')
    ax.grid(True, alpha=0.3); ax.legend()
    fig.tight_layout()
    fig.savefig(f'{out_dir}/sim_working_point.png', dpi=150, bbox_inches='tight')
    plt.close(fig)

    # summary
    nom_now = run_point(e_now, args.gap_mm / 10.0, args.n_events)
    best_e = {k: float(nom.loc[nom[k].idxmin(), 'E']) for k in estimators}
    best_s = {k: float(nom[k].min()) for k in estimators}
    lines = [
        f'P2 timing physics simulation — Ar/iC4H10 95:5 (Magboltz+HEED via '
        f'Garfield++ on lxplus), gap {args.gap_mm:g} +- {args.gap_tol_mm:g} mm, '
        f'{"vertical" if args.vertical else "cos^2 cosmic"} tracks',
        f'  HEED: {np.mean([len(t[0]) for t in lib[0.0]]):.1f} clusters / 3 mm '
        f'vertical ({np.mean([len(t[0]) for t in lib[0.0]])/0.3:.0f} /cm), '
        f'mean first-cluster depth '
        f'{0.3e4/np.mean([len(t[0]) for t in lib[0.0]]):.0f} um',
        f'  working point: {v_gap:.0f} V / {args.gap_mm:g} mm = {e_now:.0f} V/cm, '
        f'vd {vd_of(e_now):.2f} cm/us, full-gap drift '
        f'{args.gap_mm/10/vd_of(e_now)*1e3:.0f} ns',
        '',
        '  sigma_t at the CURRENT working point:']
    for k in estimators:
        lines.append(f'    {labels[k]:44s}: {nom_now[k]:5.1f} ns')
    lines += ['', '  best over the field scan (optimal E):']
    for k in estimators:
        lines.append(f'    {labels[k]:44s}: {best_s[k]:5.1f} ns at '
                     f'{best_e[k]:.0f} V/cm')
    txt = '\n'.join(lines)
    print('\n' + txt)
    with open(f'{out_dir}/sim_summary.txt', 'w') as f:
        f.write(txt + '\n')
    print(f'\nWritten to: {out_dir}')


if __name__ == '__main__':
    main()
