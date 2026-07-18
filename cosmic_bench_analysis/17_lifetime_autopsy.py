#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
17_lifetime_autopsy.py

Lifetime autopsy of det3 (P2_3) across its whole data-taking life
(2026-07-16 12:10 -> 2026-07-17 15:57): stitch EVERY sub_run of the drift-scan
run and the mesh-scan/initial/final run onto one absolute wall-clock axis and
track, in fine time bins:

  - M3 reference health: trigger rate, fraction of events with a good single
    track (chi2<5, NClus>=3)  -> is the REFERENCE stable?
  - det3 efficiency (within MATCH_R of a frozen, healthy-period pad->M3
    transform) and any-pad response
  - det3 GAIN proxies: matched-cluster total charge (sum amp) and cluster size
  - det4 any-pad response + charge while it was powered (shared-gas control,
    HV off from 21:45 on 7-16)

Also produces:
  - amplitude spectra of matched clusters per period (gain evolution)
  - per-connector spatial breakdown per period (uniform vs localized death)
  - a mesh-scan amplitude calibration amp(meshV) and the equivalent effective
    mesh voltage V_eff(t) implied by the measured cluster charge vs time

Products -> <ANALYSIS_ROOT>/det3/lifetime_autopsy/
Usage: python3 17_lifetime_autopsy.py [--bin-min 15]
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
import matplotlib.dates as mdates
import uproot
import awkward as ak
from scipy.spatial import cKDTree

import p2_qa_config as qa
import p2_mapping as pmap
import p2_align as pa

DET3_Z = 702.0
MATCH_R = 40.0
ACTIVE_R = 30.0
NOISY_PADS = (510,)

MESH_RUN = 'p2_det3_mesh_scan_det4_initial_7-16-26'
DRIFT_RUN = 'p2_det3_det4_drift_scan_7-16-26'


# --------------------------------------------------------------------------- #
def wall_t0(sub_dir):
    """Absolute start time parsed from the combined-hits datrun name."""
    fs = glob.glob(os.path.join(sub_dir, 'combined_hits_root', '*.root'))
    m = re.search(r'_datrun_(\d{2})(\d{2})(\d{2})_(\d{2})H(\d{2})_', os.path.basename(fs[0]))
    yy, mm, dd, H, M = (int(g) for g in m.groups())
    return pd.Timestamp(2000 + yy, mm, dd, H, M)


def m3_health(m3_dir):
    """Per-event M3 health frame: t_sec (run-relative), n_rays, good (>=1 ray
    with chi2<5 & NClus>=3 in both coordinates)."""
    parts = []
    for fp in sorted(glob.glob(os.path.join(m3_dir, '*.root'))):
        a = uproot.open(fp)['T'].arrays(
            ['evn', 'evttime', 'rayN', 'Chi2X', 'Chi2Y', 'NClusX', 'NClusY'])
        good = ak.any((a['Chi2X'] < qa.M3_CHI2_CUT) & (a['Chi2Y'] < qa.M3_CHI2_CUT) &
                      (a['NClusX'] >= qa.M3_MIN_NCLUS) & (a['NClusY'] >= qa.M3_MIN_NCLUS),
                      axis=1)
        parts.append(pd.DataFrame({
            'eventId': ak.to_numpy(a['evn']).astype(np.int64),
            't_sec': ak.to_numpy(a['evttime']) * 1e-8,
            'n_rays': ak.to_numpy(a['rayN']),
            'good': ak.to_numpy(good)}))
    return pd.concat(parts, ignore_index=True).drop_duplicates('eventId')


def det_events(hits_dir, ct, feus, drop_pads=()):
    """Per-event detector frame from the (clean) combined hits: t_sec, sumamp,
    maxamp, npad, charge-weighted centroid."""
    parts = []
    for fp in sorted(glob.glob(os.path.join(hits_dir, '*.root'))):
        a = uproot.open(f'{fp}:hits').arrays(
            ['eventId', 'trigger_timestamp_ns', 'channel', 'amplitude', 'feu'],
            library='pd')
        a = a[a['feu'].isin(set(feus))]
        h = pmap.attach_pads_to_hits(a, ct)
        h = h[h['mapped'] & h['pad_cx'].notna()]
        if drop_pads:
            h = h[~h['channel_id'].isin(set(drop_pads))]
        if not len(h):
            continue
        w = h['amplitude'].clip(lower=0).astype(float)
        t = h.assign(_wx=w * h['pad_cx'], _wy=w * h['pad_cy'], _w=w)
        g = t.groupby('eventId')
        parts.append(pd.DataFrame({
            't_sec': g['trigger_timestamp_ns'].first() / 1e9,
            'sumamp': g['amplitude'].sum(),
            'maxamp': g['amplitude'].max(),
            'npad': g['channel_id'].nunique(),
            'x_pad': g['_wx'].sum() / g['_w'].sum(),
            'y_pad': g['_wy'].sum() / g['_w'].sum()}).reset_index())
    if not parts:
        return pd.DataFrame(columns=['eventId', 't_sec', 'sumamp', 'maxamp',
                                     'npad', 'x_pad', 'y_pad'])
    return pd.concat(parts, ignore_index=True)


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description='det3 lifetime autopsy.')
    ap.add_argument('--bin-min', type=float, default=15.0)
    args = ap.parse_args()
    out = os.path.join(qa.ANALYSIS_ROOT, 'det3', 'lifetime_autopsy')
    os.makedirs(out, exist_ok=True)

    mesh_dir = os.path.join(qa.DATA_ROOT, MESH_RUN)
    drift_dir = os.path.join(qa.DATA_ROOT, DRIFT_RUN)

    # channel tables (per run_config; same physical wiring)
    ct3 = {}
    ct4 = {}
    for run_dir in (mesh_dir, drift_dir):
        rc = os.path.join(run_dir, 'run_config.json')
        ct3[run_dir] = pmap.build_channel_table(rc, qa.MAP_CSV_PATH, det_type='P2',
                                                det_name='P2_3', strategy='reverse',
                                                drop_connectors=(1, 8, 9, 10))
        ct4[run_dir] = pmap.build_channel_table(rc, qa.MAP_CSV_PATH, det_type='P2',
                                                det_name='P2_4', strategy='reverse',
                                                drop_connectors=(1, 10))

    # --- frozen det3 transform from the (healthy) initial run --------------- #
    init = os.path.join(mesh_dir, 'initial_run_det3_420_820_det4_430_830')
    m3i = pa.load_m3_positions(os.path.join(init, 'm3_tracking_root'), DET3_Z)
    p2i = det_events(os.path.join(init, 'combined_hits_root'),
                     ct3[mesh_dir], ct3[mesh_dir].attrs['feus'], NOISY_PADS)
    mi = m3i.merge(p2i, on='eventId')
    T = pa.fit_transform(mi['x_m3'], mi['y_m3'], mi['x_pad'], mi['y_pad'])
    print(f'frozen transform (initial run): rot {T.rotation_deg:.2f} deg, '
          f'scale {T.s:.3f}, RMSE {T.rmse:.1f} mm (N={len(mi):,})')
    padc = ct3[mesh_dir][ct3[mesh_dir]['mapped']].drop_duplicates('channel_id')
    pcx, pcy = T.apply(padc['pad_cx'].to_numpy(), padc['pad_cy'].to_numpy())
    tree = cKDTree(np.column_stack([pcx, pcy]))

    # --- enumerate sub_runs on the wall clock ------------------------------- #
    subs = ([(drift_dir, os.path.basename(s)) for s in
             sorted(glob.glob(os.path.join(drift_dir, 'drift_scan_*')))]
            + [(mesh_dir, 'initial_run_det3_420_820_det4_430_830')]
            + [(mesh_dir, os.path.basename(s)) for s in
               sorted(glob.glob(os.path.join(mesh_dir, 'mesh_scan_det3_*')))]
            + [(mesh_dir, 'final_run_det3_420_820')])

    rows = []          # per-ray rows for efficiency
    m3rows = []        # per-event M3 health rows
    d4rows = []        # det4 per-event rows
    periods = []       # (name, t0, t1, det3 HV) annotations
    for run_dir, sub in subs:
        sdir = os.path.join(run_dir, sub)
        if not glob.glob(os.path.join(sdir, 'combined_hits_root', '*.root')):
            continue
        t0 = wall_t0(sdir)
        mh = m3_health(os.path.join(sdir, 'm3_tracking_root'))
        mh['wall'] = t0 + pd.to_timedelta(mh['t_sec'], unit='s')
        mh['sub'] = sub
        m3rows.append(mh)

        m3p = pa.load_m3_positions(os.path.join(sdir, 'm3_tracking_root'), DET3_Z)
        p3 = det_events(os.path.join(sdir, 'combined_hits_root'),
                        ct3[run_dir], ct3[run_dir].attrs['feus'], NOISY_PADS)
        p4 = det_events(os.path.join(sdir, 'combined_hits_root'),
                        ct4[run_dir], ct4[run_dir].attrs['feus'])
        if len(p4):
            p4 = p4.assign(wall=t0 + pd.to_timedelta(p4['t_sec'], unit='s'), sub=sub)
            d4rows.append(p4)

        # per-ray table with frozen transform / active area
        d = m3p.merge(mh[['eventId', 't_sec']], on='eventId', how='left')
        nn, _ = tree.query(np.column_stack([d['x_m3'], d['y_m3']]))
        d = d[nn <= ACTIVE_R].copy()
        rec = p3.set_index('eventId')
        rx, ry = T.apply(rec['x_pad'].to_numpy(), rec['y_pad'].to_numpy())
        recx = pd.Series(rx, index=rec.index)
        recy = pd.Series(ry, index=rec.index)
        d['has_any'] = d['eventId'].isin(rec.index)
        dx = d['x_m3'] - d['eventId'].map(recx)
        dy = d['y_m3'] - d['eventId'].map(recy)
        d['within'] = np.hypot(dx, dy) <= MATCH_R
        d['sumamp'] = d['eventId'].map(rec['sumamp'])
        d['npad'] = d['eventId'].map(rec['npad'])
        d['wall'] = t0 + pd.to_timedelta(d['t_sec'], unit='s')
        d['sub'] = sub
        rows.append(d)

        mesh_v = (re.search(r'det3_(\d+)_', sub).group(1)
                  if 'det3_' in sub else '420')
        periods.append((sub, t0, mh['wall'].max(), mesh_v))
        print(f'{sub}: rays {len(d):,}  within {d["within"].mean()*100:5.1f}%  '
              f'any {d["has_any"].mean()*100:5.1f}%  ({t0:%d %H:%M})')

    rays = pd.concat(rows, ignore_index=True)
    m3h = pd.concat(m3rows, ignore_index=True)
    det4 = pd.concat(d4rows, ignore_index=True)

    # ----------------------------------------------------------------------- #
    # time-binned series
    # ----------------------------------------------------------------------- #
    b = f'{int(args.bin_min)}min'
    rb = rays.set_index('wall').resample(b)
    eff = rb.agg(eff=('within', 'mean'), any=('has_any', 'mean'),
                 n=('within', 'size'))
    matched = rays[rays['within']].set_index('wall')
    amp = matched.resample(b).agg(amp_med=('sumamp', 'median'),
                                  amp_q25=('sumamp', lambda s: s.quantile(.25)),
                                  amp_q75=('sumamp', lambda s: s.quantile(.75)),
                                  npad=('npad', 'mean'), n=('sumamp', 'size'))
    m3b = m3h.set_index('wall').resample(b).agg(
        rate=('eventId', lambda s: len(s) / (args.bin_min * 60)),
        frac_good=('good', 'mean'))
    d4b = det4.set_index('wall').resample(b).agg(
        n=('sumamp', 'size'), amp_med=('sumamp', 'median'))

    eff.to_csv(os.path.join(out, 'det3_eff_timeline.csv'))
    amp.to_csv(os.path.join(out, 'det3_amp_timeline.csv'))
    m3b.to_csv(os.path.join(out, 'm3_health_timeline.csv'))
    d4b.to_csv(os.path.join(out, 'det4_control_timeline.csv'))

    # ----------------------------------------------------------------------- #
    # master figure
    # ----------------------------------------------------------------------- #
    fig, axes = plt.subplots(4, 1, figsize=(15, 13), sharex=True,
                             gridspec_kw=dict(height_ratios=[1.4, 1, 1, 1]))
    ok = eff['n'] >= 100
    axes[0].plot(eff.index[ok], 100 * eff['eff'][ok], 'o-', ms=4,
                 color='steelblue', label=f'within {MATCH_R:.0f} mm')
    axes[0].plot(eff.index[ok], 100 * eff['any'][ok], 's--', ms=3, alpha=0.7,
                 color='darkorange', label='any pad')
    axes[0].set_ylabel('det3 efficiency [%]')
    axes[0].set_ylim(0, 100)
    axes[0].legend(loc='upper right', fontsize=9)

    okm = amp['n'] >= 20
    axes[1].plot(amp.index[okm], amp['amp_med'][okm], 'o-', ms=4, color='seagreen',
                 label='matched cluster charge (median)')
    axes[1].fill_between(amp.index[okm], amp['amp_q25'][okm], amp['amp_q75'][okm],
                         alpha=0.25, color='seagreen', label='25-75%')
    axes[1].set_ylabel('cluster Σamp [ADC]')
    axes[1].legend(loc='upper right', fontsize=9)
    ax1b = axes[1].twinx()
    ax1b.plot(amp.index[okm], amp['npad'][okm], ':', color='grey', alpha=0.8)
    ax1b.set_ylabel('pads/cluster', color='grey')

    axes[2].plot(m3b.index, m3b['rate'], '-', color='slategrey', label='M3 trigger rate')
    axes[2].set_ylabel('M3 rate [Hz]')
    ax2b = axes[2].twinx()
    ax2b.plot(m3b.index, 100 * m3b['frac_good'], '-', color='purple', alpha=0.7)
    ax2b.set_ylabel('good single track [%]', color='purple')
    ax2b.set_ylim(0, 60)
    axes[2].legend(loc='upper left', fontsize=9)

    axes[3].plot(d4b.index, d4b['n'] / (args.bin_min * 60), 'o-', ms=3,
                 color='firebrick', label='det4 signal-event rate (430/830 until 21:45)')
    axes[3].set_ylabel('det4 rate [Hz]')
    ax3b = axes[3].twinx()
    ax3b.plot(d4b.index, d4b['amp_med'], 's--', ms=3, color='chocolate', alpha=0.7)
    ax3b.set_ylabel('det4 Σamp med [ADC]', color='chocolate')
    axes[3].legend(loc='upper right', fontsize=9)
    axes[3].set_xlabel('wall clock')

    # annotations: HV values along the top + notable moments
    for sub, ta, tb, meshv in periods:
        axes[0].axvspan(ta, tb, color='k', alpha=0.03)
        if 'mesh_scan' in sub:
            axes[0].text(ta, 96, meshv, fontsize=6, rotation=90, va='top')
    for t, lbl, c in [
            (pd.Timestamp(2026, 7, 16, 12, 0), 'gas 6→8 L/h (mixer)', 'teal'),
            (pd.Timestamp(2026, 7, 16, 17, 47), 'HV off (run end)', 'red'),
            (pd.Timestamp(2026, 7, 16, 19, 42), 'HV on', 'green'),
            (pd.Timestamp(2026, 7, 16, 21, 45), 'det4 HV→0 / mesh scan starts', 'firebrick'),
            (pd.Timestamp(2026, 7, 17, 5, 56), 'final run starts (420/820)', 'navy')]:
        for ax in axes:
            ax.axvline(t, color=c, ls='--', lw=1, alpha=0.6)
        axes[0].text(t, 103, lbl, fontsize=7, color=c, rotation=20, ha='left')
    axes[0].set_title('det3 (P2_3) lifetime autopsy — 2026-07-16/17 — frozen transform, '
                      'clean (re-pedestalled) data, min_amp=0', pad=30)
    for ax in axes:
        ax.grid(True, alpha=0.3)
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%d/%m %H:%M'))
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(os.path.join(out, 'lifetime_master.png'), dpi=170,
                bbox_inches='tight')
    plt.close(fig)

    # ----------------------------------------------------------------------- #
    # amplitude spectra + spatial breakdown per period
    # ----------------------------------------------------------------------- #
    period_sel = [
        ('drift scan 820V (16th 17:15)', rays['sub'] == 'drift_scan_det4_450_850_det3_420_820'),
        ('initial run (16th 19:43)', rays['sub'] == 'initial_run_det3_420_820_det4_430_830'),
        ('mesh scan 420V (16th 21:45)', rays['sub'] == 'mesh_scan_det3_420_820'),
        ('final run (17th 05:56+)', rays['sub'] == 'final_run_det3_420_820'),
    ]
    fig2, ax2 = plt.subplots(figsize=(9, 5.5))
    bins = np.geomspace(20, 20000, 60)
    for name, sel in period_sel:
        s = rays.loc[sel & rays['within'], 'sumamp'].dropna()
        if len(s) < 10:
            ax2.plot([], [], label=f'{name} (n={len(s)})')
            continue
        ax2.hist(s, bins=bins, histtype='step', lw=1.8, density=True,
                 label=f'{name}  med={s.median():.0f} (n={len(s):,})')
    ax2.set_xscale('log')
    ax2.set_xlabel('matched cluster Σamp [ADC]')
    ax2.set_ylabel('normalised')
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)
    ax2.set_title('det3 matched-cluster charge spectra per period (gain evolution)')
    fig2.tight_layout()
    fig2.savefig(os.path.join(out, 'amp_spectra_periods.png'), dpi=170,
                 bbox_inches='tight')
    plt.close(fig2)

    # spatial: efficiency per pad-plane cell per period (frozen frame)
    fig3, ax3s = plt.subplots(2, 4, figsize=(19, 8.5))
    for i, (name, sel) in enumerate(period_sel):
        d = rays[sel]
        for j, col in enumerate(('within', 'has_any')):
            ax = ax3s[j][i]
            if len(d) < 50:
                ax.set_axis_off()
                continue
            hxy, xe, ye = np.histogram2d(d['x_m3'], d['y_m3'], bins=18)
            hin, _, _ = np.histogram2d(d.loc[d[col], 'x_m3'], d.loc[d[col], 'y_m3'],
                                       bins=[xe, ye])
            with np.errstate(invalid='ignore'):
                emap = np.where(hxy >= 5, hin / hxy, np.nan)
            im = ax.imshow(emap.T, origin='lower', vmin=0, vmax=1, cmap='RdYlGn',
                           extent=[xe[0], xe[-1], ye[0], ye[-1]], aspect='equal')
            ax.set_title(f'{name}\n{col}', fontsize=9)
            plt.colorbar(im, ax=ax, fraction=0.046)
    fig3.suptitle('det3 spatial efficiency per period (M3 frame, ≥5 rays/cell) — '
                  'uniform fade vs localized death')
    fig3.tight_layout()
    fig3.savefig(os.path.join(out, 'spatial_periods.png'), dpi=160,
                 bbox_inches='tight')
    plt.close(fig3)

    # ----------------------------------------------------------------------- #
    # mesh-scan amplitude calibration -> V_eff(t)
    # ----------------------------------------------------------------------- #
    cal = (rays[rays['sub'].str.startswith('mesh_scan') & rays['within']]
           .assign(v=lambda d: d['sub'].str.extract(r'det3_(\d+)_')[0].astype(float))
           .groupby('v')['sumamp'].agg(['median', 'size']))
    cal = cal[cal['size'] >= 15]
    print('\namp(meshV) calibration (matched clusters):')
    print(cal.round(1).to_string())
    if len(cal) >= 3:
        lv = np.log(cal['median'].to_numpy())
        vv = cal.index.to_numpy()
        pfit = np.polyfit(vv, lv, 1)   # ln(amp) = a*V + b (exponential gain)
        slope_dec = np.log(10) / pfit[0]
        print(f'exponential fit: amp ~ exp({pfit[0]:.4f} * V)  '
              f'(x10 every {slope_dec:.1f} V)')
        ok2 = amp['n'] >= 20
        veff = (np.log(amp['amp_med'][ok2]) - pfit[1]) / pfit[0]
        fig4, ax4 = plt.subplots(figsize=(13, 4.5))
        ax4.plot(amp.index[ok2], veff, 'o-', ms=4, color='crimson')
        ax4.axhline(420, color='k', ls='--', alpha=0.6)
        ax4.text(amp.index[ok2][0], 421, 'applied mesh V (420)', fontsize=8)
        ax4.set_ylabel('effective mesh V implied by gain [V]')
        ax4.set_xlabel('wall clock')
        ax4.grid(True, alpha=0.3)
        ax4.xaxis.set_major_formatter(mdates.DateFormatter('%d/%m %H:%M'))
        ax4.set_title('det3 effective mesh voltage from matched-cluster charge '
                      '(mesh-scan amplitude calibration)')
        fig4.autofmt_xdate()
        fig4.tight_layout()
        fig4.savefig(os.path.join(out, 'veff_timeline.png'), dpi=170,
                     bbox_inches='tight')
        plt.close(fig4)

    print(f'\nWritten to: {out}')


if __name__ == '__main__':
    main()
