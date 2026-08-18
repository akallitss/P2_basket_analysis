#!/usr/bin/env python3
"""Efficiency vs time WITHIN sub-runs (charging-up study, taskplan T1)."""
import os, sys, glob, json
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2style as st
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from paths import S, A, URW, RD, OUT  # noqa: E402
DETS = ['P2_IN', 'P2_MID', 'P2_OUT']

rows = []
for f in glob.glob(f'{S}/products/urw_timebins_x/urw_timebins/*/*/urw_p2_efficiency_*.json'):
    try:
        lst = json.load(open(f))
    except Exception:
        continue
    for e in lst:
        for b in e.get('efficiency_vs_time', []):
            rows.append(dict(run=e['run'], sub_run=e['sub_run'],
                             det=e['station'], eff_sub=e['efficiency']['value'],
                             **b))
tb = pd.DataFrame(rows)
tb.to_csv(f'{RD}/eff_vs_time_bins.csv', index=False)
print(f'{len(tb)} time bins from '
      f'{tb.groupby(["run","sub_run"]).ngroups} sub_runs')

# sub_run wall-clock start times from the HV monitor
hv = pd.read_csv(f'{RD}/hv_campaign_30s.csv', parse_dates=['time'])
t0 = hv.groupby(['run', 'sub_run'])['time'].min().to_dict()

for run in sorted(tb.run.unique()):
    g0 = tb[tb.run == run]
    fig, ax = plt.subplots(figsize=(11.5, 5.0))
    for det in DETS:
        g = g0[g0.det == det]
        if not len(g):
            continue
        xs, ys, los, his = [], [], [], []
        for sub in sorted(g.sub_run.unique()):
            gg = g[g.sub_run == sub].sort_values('t_mid_s')
            start = t0.get((run, sub))
            if start is None:
                continue
            rel = gg.t_mid_s - gg.t_lo_s.min()
            x = start + pd.to_timedelta(rel, unit='s')
            ax.errorbar(x, gg.eff, yerr=[gg.eff - gg.lo, gg.hi - gg.eff],
                        color=st.DET_COLOR[det], marker='o', ms=3.5, lw=1.4,
                        capsize=0, alpha=0.9)
        # one legend entry per det
        ax.plot([], [], color=st.DET_COLOR[det], marker='o', lw=1.4, label=det)
    ax.legend(loc='lower left', ncol=3)
    ax.set_ylabel('absolute efficiency (uRWELL-referenced)')
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    day = t0.get((run, sorted(g0.sub_run.unique())[0]))
    ax.set_xlabel(f'time of day, {day:%d %b %Y}' if day is not None else 'time')
    ax.set_title(f'Efficiency vs time within sub-runs — {run} '
                 f'(2-min bins, constant HV)')
    st.finish(fig, f'{OUT}/eff_vs_time_{run}.png')
print('done')
