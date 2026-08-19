#!/usr/bin/env python3
"""Roadmap 5b.5 -- track-level performance vs beam rate, from muons alone.

No pion run was taken, so there is no high-density point.  What the campaign
does have is the time structure *inside* each sub-run: `urw_p2_efficiency.py
--time-bins` records, per bin, how many reference tracks arrived and what
fraction the station matched.  Dividing the two gives an instantaneous track
rate and an efficiency measured under exactly the same conditions -- same
chamber, same HV, same gas, same alignment, minutes apart.  That is a real
rate dependence, and it is free of the trap that killed the first attempt:
sub-run-average rates are dominated by dead time and by which run it was, not
by the beam.

Second panel: the accidental (fake) match probability per 10 mm probe radius,
which is what a 3-plane coincidence has to reject.  It is a per-sub-run number,
so it is plotted against the sub-run's mean track rate.

P2_IN is drawn per chamber, never as one series -- three different chambers sat
in that station during the campaign (see ../chamber_history.py).
"""

import glob
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import p2style as st  # noqa: E402
import chamber_history as ch  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

from paths import S, RD, OUT  # noqa: E402

DETS = ['P2_MID', 'P2_OUT']       # det1 and det3: unchanged all campaign
MIN_TRACKS = 3000                 # per bin, so the error bar stays small


def load_bins():
    rows = []
    for f in glob.glob(f'{S}/products/urw_timebins_x/urw_timebins/*/*/'
                       'urw_p2_efficiency_*.json'):
        try:
            blocks = json.load(open(f))
        except (OSError, ValueError):
            continue
        for b in blocks:
            for t in (b.get('efficiency_vs_time') or []):
                dt = t['t_hi_s'] - t['t_lo_s']
                if dt <= 0 or not t['n']:
                    continue
                rows.append(dict(
                    run=b['run'], sub_run=b['sub_run'], station=b['station'],
                    mesh_hv=b.get('mesh_hv'), drift_hv=b.get('drift_hv'),
                    t_mid_s=t['t_mid_s'], n=t['n'], eff=t['eff'],
                    lo=t['lo'], hi=t['hi'], rate_hz=t['n'] / dt))
    return pd.DataFrame(rows)


def main():
    b = load_bins()
    b = b[(b.n >= MIN_TRACKS) & (b.mesh_hv == 450) &
          b.drift_hv.between(690, 760) & b.station.isin(DETS)]
    print(f'{len(b)} time bins at the working point, '
          f'{b.sub_run.nunique()} sub-runs, {b.run.nunique()} runs')

    per = pd.concat(
        [pd.read_csv(f) for f in glob.glob(
            f'{S}/products/urw_timebins_x/urw_timebins/*/*/'
            'urw_p2_efficiency_*.csv')], ignore_index=True)
    mean_rate = b.groupby(['run', 'sub_run'], as_index=False).rate_hz.mean()
    acc = (per[(per.mesh_hv == 450) & per.drift_hv.between(690, 760)]
           .merge(mean_rate, on=['run', 'sub_run'], how='inner'))

    # Fit WITHIN each sub-run.  Across sub-runs the efficiency has genuine
    # offsets -- P2_MID sits ~4 points higher on 25 Jul than on 27-28 Jul at
    # the same nominal working point -- and those offsets happen to correlate
    # with rate, so a single fit through everything returns +12 points/kHz,
    # which is the run step and not the beam.  Inside one sub-run nothing
    # changes but the load.
    slopes = []
    for (run, sub, det), g in b.groupby(['run', 'sub_run', 'station']):
        if len(g) < 5 or g.rate_hz.max() / g.rate_hz.min() < 1.15:
            continue
        a, b0 = np.polyfit(g.rate_hz, g.eff, 1)
        slopes.append(dict(run=run, sub_run=sub, station=det,
                           slope_pts_per_khz=a * 1e3 * 100,
                           rate_lo=g.rate_hz.min(), rate_hi=g.rate_hz.max(),
                           a=a, b=b0, n=len(g)))
    sl = pd.DataFrame(slopes)

    fig, axes = plt.subplots(1, 3, figsize=(16.4, 5.0))
    ax, axs, ax2 = axes

    runs = sorted(b.run.unique())
    rcol = dict(zip(runs, [st.C_VIO, st.C_YEL, st.C_MAG, st.C_GRN][:len(runs)]))
    for det, mk in zip(DETS, ('o', 's')):
        g = b[b.station == det]
        for run, gg in g.groupby('run'):
            ax.plot(gg.rate_hz, gg.eff, mk, ms=3.4, alpha=0.45,
                    color=rcol[run], mec='none')
    for _, r in sl.iterrows():
        xs = np.array([r.rate_lo, r.rate_hi])
        ax.plot(xs, r.a * xs + r.b, lw=1.1, alpha=0.85,
                color=st.DET_COLOR[r.station])
    for run in runs:
        ax.plot([], [], 'o', color=rcol[run], label=run)
    ax.plot([], [], '-', color=st.DET_COLOR['P2_MID'], label='P2_MID fits')
    ax.plot([], [], '-', color=st.DET_COLOR['P2_OUT'], label='P2_OUT fits')
    ax.set_xlabel('uRWELL track rate in the time bin [Hz]')
    ax.set_ylabel('absolute efficiency')
    ax.set_title('Each line is one sub-run\n(the offsets between runs are '
                 'not a rate effect)', fontsize=11)
    ax.legend(loc='lower right', fontsize=8)

    for det in DETS:
        s = sl[sl.station == det].slope_pts_per_khz
        if not len(s):
            continue
        axs.hist(s, bins=np.linspace(-12, 12, 25), histtype='step', lw=2,
                 color=st.DET_COLOR[det],
                 label=f'{det}: {s.mean():+.1f} $\\pm$ '
                       f'{s.std(ddof=1) / np.sqrt(len(s)):.1f} pts/kHz '
                       f'({len(s)} sub-runs)')
        axs.axvline(s.mean(), color=st.DET_COLOR[det], ls='--', lw=1.4)
        print(f'  {det}: within-sub_run slope '
              f'{s.mean():+.2f} +- {s.std(ddof=1) / np.sqrt(len(s)):.2f} '
              f'points per kHz  ({len(s)} sub-runs)')
    axs.axvline(0, color=st.TEXT2, lw=1.2)
    axs.set_xlabel('d(efficiency) / d(rate)  [points per kHz of tracks]')
    axs.set_ylabel('sub-runs')
    axs.set_title('Within-sub-run rate dependence\nconsistent with zero',
                  fontsize=11)
    axs.legend(loc='upper left', fontsize=8.5)

    for det in DETS + ['P2_IN']:
        g = acc[acc.station == det]
        if not len(g):
            continue
        if det == 'P2_IN':
            for run, gg in g.groupby('run'):
                ax2.plot(gg.rate_hz, gg.accidental_per_10mm * 100, '^', ms=5,
                         alpha=0.8,
                         label=f'{ch.label("P2_IN", run)}  [{run}]')
        else:
            ax2.plot(g.rate_hz, g.accidental_per_10mm * 100, 'o', ms=5.5,
                     color=st.DET_COLOR[det],
                     label=ch.label_any(det))
    ax2.set_xlabel('mean uRWELL track rate in the sub-run [Hz]')
    ax2.set_ylabel('accidental match per 10 mm probe [%]')
    ax2.set_title('Fake-match background\n0.05–0.11 % over the whole campaign',
                  fontsize=11)
    ax2.set_ylim(0, 0.13)
    ax2.legend(loc='lower left', fontsize=8)

    tr = pd.read_csv(f'{RD}/dream_raw_stream.csv')
    ntr = tr.groupby(['run', 'sub_run'], as_index=False).n_triggers.median()
    ntk = (per[per.station == 'P2_MID'][['run', 'sub_run', 'n']]
           .merge(ntr, on=['run', 'sub_run'], how='inner'))
    ratio = float((ntk.n_triggers / ntk.n).median()) if len(ntk) else np.nan
    if np.isfinite(ratio):
        hi = b.rate_hz.max() * ratio / 1e3
        ax.text(0.03, 0.55, f'{ratio:.1f} triggers per reference track,\n'
                            f'so the top of this axis is\n'
                            f'~{hi:.1f} kHz of beam trigger',
                transform=ax.transAxes, fontsize=9, color='#444',
                bbox=dict(fc='white', ec='#ddd', boxstyle='round,pad=0.3'))
        print(f'  {ratio:.2f} triggers per track -> top of axis '
              f'~{hi:.1f} kHz trigger rate')

    fig.suptitle('Track-level performance vs beam load — roadmap 5b.5 '
                 '(muons; no pion run was taken)\n'
                 'working point only (mesh 450 V, drift 700–750 V); each point '
                 'is one time bin inside a sub-run, so chamber, HV and gas are '
                 'held fixed', fontweight='bold')
    st.finish(fig, f'{OUT}/rate_performance.png')


if __name__ == '__main__':
    main()
