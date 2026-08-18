#!/usr/bin/env python3
"""Stability figures: campaign HV timeline, sparks, efficiency drift."""
import os, sys, glob, json
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2style as st
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from paths import S, A, URW, RD, OUT  # noqa: E402
DETS = ['P2_IN', 'P2_MID', 'P2_OUT']

# ------------------------------------------------ S1: campaign HV timeline
hv = pd.read_csv(f'{RD}/hv_campaign_30s.csv', parse_dates=['time'])
hv = hv[hv.time > '2026-07-23']
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13.5, 7.2), sharex=True,
                               height_ratios=[2, 1])
for det in DETS:
    g = hv[(hv.det == det) & (hv.electrode == 'mesh')].sort_values('time')
    # gaps between runs: break the line when >30 min silent
    t = g.time.to_numpy()
    v = g.vmon.to_numpy()
    brk = np.where(np.diff(t) > np.timedelta64(30, 'm'))[0]
    v2 = v.astype(float).copy()
    tt = t.copy()
    ax1.plot(np.insert(tt, brk + 1, tt[brk] + np.timedelta64(1, 'm')),
             np.insert(v2, brk + 1, np.nan),
             color=st.DET_COLOR[det], lw=1.2, label=f'{det} mesh')
ax1.set_ylabel('mesh voltage readback [V]')
ax1.set_ylim(0, 480)
ax1.legend(loc='lower left', fontsize=9, ncol=3)
# gas periods
gas_swap0 = pd.Timestamp('2026-08-01 21:19')
gas_swap1 = pd.Timestamp('2026-08-02 11:30')
ax1.axvspan(gas_swap0, gas_swap1, color=st.GRID, alpha=0.5, lw=0)
ax1.annotate('gas exchange\n(P2 HV off)', xy=(gas_swap0, 460), fontsize=8.5,
             color=st.TEXT2, ha='left', va='top')
ax1.annotate('Ar/CO$_2$/iso 93/5/2', xy=(pd.Timestamp('2026-07-26'), 466),
             fontsize=10, color=st.TEXT2, fontweight='bold')
ax1.annotate('Ar/CF$_4$/iso 88/10/2', xy=(pd.Timestamp('2026-08-02 12:00'), 466),
             fontsize=10, color=st.TEXT2, fontweight='bold')
ax1.annotate('P2_IN chamber swap', xy=(pd.Timestamp('2026-07-28 09:00'), 30),
             fontsize=8.5, color=st.C_IN, rotation=90, va='bottom')

for det in DETS:
    g = hv[(hv.det == det) & (hv.electrode == 'mesh')].sort_values('time')
    ax2.plot(g.time, g.imon_max.clip(upper=2.2), color=st.DET_COLOR[det],
             lw=0.9, alpha=0.85)
ax2.set_ylabel('mesh current, 30 s max [µA]')
ax2.set_yscale('symlog', linthresh=0.1)
ax2.set_ylim(0, 2.4)
ax2.xaxis.set_major_formatter(mdates.DateFormatter('%d %b'))
ax2.set_xlabel('2026')
fig.suptitle('Campaign HV timeline — mesh voltage and current, all stations',
             fontweight='bold')
st.finish(fig, f'{OUT}/hv_campaign_timeline.png')

# ------------------------------------------------- S2: spark summary per det
sp = pd.read_csv(f'{RD}/dream_sparks.csv')
sp = sp[~sp.run.str.startswith('run_') | (sp.run == 'run_1')]  # July P2 runs
tot = (sp.groupby('det')
         .agg(n_sparks=('n_sparks', 'sum'), hours=('duration_s', lambda x: x.sum()/3600),
              worst_live=('live_fraction', 'min')))
tot['sparks_per_h'] = tot.n_sparks / tot.hours
fig, (axa, axb) = plt.subplots(1, 2, figsize=(12.6, 4.8))
dets = [d for d in DETS if d in tot.index]
y = np.arange(len(dets))
axa.barh(y, [tot.loc[d, 'sparks_per_h'] for d in dets],
         color=[st.DET_COLOR[d] for d in dets], height=0.55)
axa.set_yticks(y)
axa.set_yticklabels(dets)
axa.invert_yaxis()
for i, d in enumerate(dets):
    axa.text(tot.loc[d, 'sparks_per_h'] + 0.02, i,
             f"{tot.loc[d, 'n_sparks']:.0f} sparks / {tot.loc[d, 'hours']:.0f} h",
             va='center', fontsize=9, color=st.TEXT2)
axa.set_xlabel('spark rate [sparks / h] — whole July campaign')
axa.set_title('Sparks per station (imon > 2 µA excursions)')

# spark rate vs drift voltage for P2_OUT (the discharging station)
g = sp[(sp.det == 'P2_OUT')].copy()
if 'mesh_v' in g.columns:
    gg = g.groupby('mesh_v').agg(rate=('spark_rate_per_min', 'mean'),
                                 n=('n_sparks', 'sum')).reset_index()
d2 = sp[sp.det == 'P2_OUT']
axb.scatter(d2.mesh_v, d2.spark_rate_per_min * 60, s=22, alpha=0.6,
            color=st.C_OUT, edgecolors='none')
axb.set_xlabel('mesh voltage [V]')
axb.set_ylabel('spark rate [sparks / h]')
axb.set_title('P2_OUT spark rate vs mesh HV (per sub-run)')
fig.suptitle('HV stability — spark activity (Ar/CO$_2$/iso period)',
             fontweight='bold')
st.finish(fig, f'{OUT}/spark_summary.png')

# ----------------------------- S3: efficiency vs time (charging), if present
tb = []
for f in glob.glob(f'{S}/products/urw_timebins/*/*/urw_p2_efficiency_*.json'):
    try:
        d = json.load(open(f))
    except Exception:
        continue
    run = f.split('/')[-3]
    for sub, stations in (d.get('sub_runs') or {}).items():
        for det, summ in (stations or {}).items():
            for b in summ.get('efficiency_vs_time', []):
                tb.append(dict(run=run, sub_run=sub, det=det, **b))
if tb:
    tbd = pd.DataFrame(tb)
    tbd.to_csv(f'{RD}/eff_vs_time_bins.csv', index=False)
    print(f'eff_vs_time_bins.csv {len(tbd)} rows')
else:
    print('no time-bin JSONs yet (condor fleet still running)')

# ------------------------------ S4: eff drift across highstat (from urw csv)
u = pd.read_csv(f'{S}/products/urw_local/urw_referenced_efficiency/'
                'highstat_eff_1/urw_p2_efficiency_highstat_eff_1.csv')
fig, ax = plt.subplots(figsize=(8.6, 5.0))
order = sorted(u.sub_run.unique())
x = np.arange(len(order))
for det in DETS:
    g = u[u.station == det].set_index('sub_run').reindex(order)
    ax.errorbar(x, g.eff, yerr=[g.eff - g.lo, g.hi - g.eff],
                color=st.DET_COLOR[det], marker='o', lw=2, capsize=2.5,
                label=det)
    st.direct_label(ax, x[-1], g.eff.iloc[-1], det, st.DET_COLOR[det])
ax.axvspan(4.5, 5.5, color=st.GRID, alpha=0.4, lw=0)
ax.annotate('after 2.5 h pause,\nrate 4.6→3.2 kHz', xy=(5, 0.9445),
            fontsize=8.5, color=st.TEXT2, ha='center')
ax.set_xticks(x)
ax.set_xticklabels([s.replace('beam_commissioning_', 'sub ') for s in order])
ax.set_ylabel('absolute efficiency')
ax.set_xlabel('highstat_eff_1 sub-runs (~45 min each, 4.6 kHz)')
ax.set_title('Efficiency drift over 4.5 h at fixed HV — charging or rate?')
ax.legend(loc='lower left')
st.finish(fig, f'{OUT}/eff_drift_highstat.png')
print('done')
