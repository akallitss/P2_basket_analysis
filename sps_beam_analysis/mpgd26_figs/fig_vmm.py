#!/usr/bin/env python3
"""VMM figures: gas A/B turn-on, drift latency, config scan, DREAM comparison."""
import os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2style as st
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import chamber_history as ch

from paths import S, A, URW, RD, OUT  # noqa: E402
DETS = ['P2_IN', 'P2_MID', 'P2_OUT']
DY = {'P2_IN': 0, 'P2_MID': 9, 'P2_OUT': -9}

v = pd.read_csv(f'{RD}/vmm_subruns.csv')
v = v[v.n_triggers > 5e4]                     # drop starved points

# ---------------------------- V1: efficiency turn-on, gas A vs gas B, 3 panels
# gas A meshscans: run_31+32 (full, mesh 350-450); gas B: run_61.
mesh = v[v.sub_run.str.startswith('meshscan')].copy()
GASSETS = [('run_32', st.GAS_A, '-', 'o', 'run_32 (31 Jul)'),
           ('run_61', st.GAS_B, '--', '^', 'run_61 (2 Aug)')]
fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.8), sharey=True)
for ax, det in zip(axes, DETS):
    for run, gas, ls, mk, lab in GASSETS:
        g = mesh[(mesh.run == run) & (mesh.det == det)].sort_values('mesh_v')
        if not len(g):
            continue
        ax.errorbar(g.mesh_v, g.eff, yerr=g.eff_err, color=st.DET_COLOR[det],
                    ls=ls, marker=mk, capsize=2.5, lw=2,
                    alpha=1.0 if gas == st.GAS_A else 0.75,
                    label=f'{st.GAS_SHORT[gas]}  [{lab}]')
    ax.set_title(ch.label(det, 'run_32'), color=st.DET_COLOR[det])
    ax.set_xlabel('mesh voltage [V]')
    ax.set_ylim(-0.02, 0.75)
axes[0].set_ylabel('trigger-referenced efficiency\n(whole station, c4–c6 instrumented)')
axes[0].legend(loc='upper left', fontsize=9)
fig.suptitle('VMM efficiency turn-on: old vs new gas (same detectors, same '
             'electronics config)', fontweight='bold')
st.finish(fig, f'{OUT}/vmm_turnon_gasAB.png')

# ------------------------- V2: coincidence latency + sigma vs drift, gas A/B
dr = v[v.sub_run.str.startswith('driftscan') & (v.det == 'P2_MID')].copy()
SETS = [('run_26', st.GAS_A, st.C_MID, '-', 'o', 'run_26 (30 Jul)'),
        ('run_57', st.GAS_A, st.TEXT2, '-', 's', 'run_57 (1 Aug)'),
        ('run_61', st.GAS_B, st.C_VIO, '--', '^', 'run_61 (2 Aug)'),
        ('run_62', st.GAS_B, st.C_MAG, '--', 'v', 'run_62 (2 Aug)')]
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.2, 5.0))
for run, gas, col, ls, mk, lab in SETS:
    g = dr[dr.run == run].sort_values('drift_v')
    if not len(g):
        continue
    ax1.plot(g.drift_v, g.mu_ns, color=col, ls=ls, marker=mk, lw=2,
             label=f'{lab} — {st.GAS_SHORT[gas]}')
    ax2.plot(g.drift_v, g.sigma_ns, color=col, ls=ls, marker=mk, lw=2,
             label=f'{lab} — {st.GAS_SHORT[gas]}')
ax1.set_xlabel('drift voltage [V] (mesh 450 V)')
ax1.set_ylabel('trigger coincidence latency $\\mu$ [ns]')
ax1.set_title('Hit latency vs drift HV — the drift-velocity handle')
ax1.legend(fontsize=8.5)
ax2.set_xlabel('drift voltage [V] (mesh 450 V)')
ax2.set_ylabel('coincidence width $\\sigma$ [ns]')
ax2.set_title('Coincidence width vs drift HV (BCID-quantised readout)')
ax2.legend(fontsize=8.5)
fig.suptitle('P2_MID on VMM readout: gas comparison in the drift scan',
             fontweight='bold')
st.finish(fig, f'{OUT}/vmm_drift_gasAB.png')

# ------------------------------------------------- V3: config scan matrix
cfg = v[v.sub_run.str.startswith('cfg_')].copy()
cfg = cfg[cfg.gas == st.GAS_A]                      # keep one gas for clarity
cfg['gain'] = cfg.sub_run.str.extract(r'gain([\d.]+)').astype(float)
cfg['pt'] = cfg.sub_run.str.extract(r'peaktime(\d+)').astype(float)
cfg['variant'] = np.where(cfg.sub_run.str.contains('deflt'), 'deflt',
                  np.where(cfg.sub_run.str.contains('opt'), 'opt', 'base'))
agg = (cfg.groupby(['det', 'gain', 'pt', 'variant'], as_index=False)
          .agg(eff=('eff', 'median'), sigma_ns=('sigma_ns', 'median'),
               n=('n_triggers', 'sum')))
fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.0))
for ax, (col, lab) in zip(axes, [('eff', 'trigger-referenced efficiency'),
                                 ('sigma_ns', 'coincidence $\\sigma$ [ns]')]):
    for det in DETS:
        for gain, mk, ls in [(3.0, 'o', '-'), (4.5, '^', '--')]:
            g = agg[(agg.det == det) & (agg.gain == gain)
                    & (agg.variant == 'base')].sort_values('pt')
            if not len(g):
                continue
            ax.plot(g.pt, g[col], color=st.DET_COLOR[det], marker=mk, ls=ls,
                    lw=1.8, label=f'{det}  gain {gain} mV/fC')
    ax.set_xscale('log')
    ax.set_xticks([25, 50, 100, 200])
    ax.set_xticklabels(['25', '50', '100', '200'])
    ax.set_xlabel('VMM peaking time [ns]')
    ax.set_ylabel(lab)
axes[0].legend(fontsize=8, ncol=2)
fig.suptitle('VMM config scan (Ar/CO$_2$/iso runs, 30 Jul – 1 Aug): '
             'efficiency and timing vs shaping', fontweight='bold')
st.finish(fig, f'{OUT}/vmm_config_scan.png')

# also dump the full config table for the report
agg.sort_values(['det', 'gain', 'pt', 'variant']).to_csv(
    f'{RD}/vmm_config_scan_table.csv', index=False)

# ------------------------------------------ V4: DREAM vs VMM on one detector
fig, ax = plt.subplots(figsize=(8.0, 5.2))
u = pd.read_csv(f'{S}/products/urw_local/urw_referenced_efficiency/'
                'drift_mesh_scan_1/urw_p2_efficiency_drift_mesh_scan_1.csv')
for det in DETS:
    g = u[(u.sub_run.str.startswith('meshscan')) & (u.station == det)].sort_values('mesh_hv')
    ax.plot(g.mesh_hv, g.eff, color=st.DET_COLOR[det], marker='o', lw=2,
            label=f'{det} — DREAM, absolute (uRWELL tracks)')
    gv = mesh[(mesh.run == 'run_32') & (mesh.det == det)].sort_values('mesh_v')
    ax.plot(gv.mesh_v, gv.eff, color=st.DET_COLOR[det], marker='s', lw=1.6,
            ls='--', alpha=0.65,
            label=f'{det} — VMM, whole-station vs trigger')
ax.set_xlabel('mesh voltage [V]')
ax.set_ylabel('efficiency')
ax.set_ylim(0, 1.02)
ax.legend(fontsize=8, loc='upper left')
ax.set_title('DREAM vs VMM readout on the same detectors — NOT the same\n'
             'denominator: VMM counts every trigger, only c4–c6 instrumented')
st.finish(fig, f'{OUT}/vmm_vs_dream_caveat.png')
print('done')
