#!/usr/bin/env python3
"""Timing figures: sigma vs mesh/drift, Magboltz overlay, 2D surface, ladder."""
import os, sys, json, glob
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2style as st
import matplotlib.pyplot as plt

from paths import S, A, URW, RD, OUT  # noqa: E402
DETS = ['P2_IN', 'P2_MID', 'P2_OUT']
DY = {'P2_IN': 0, 'P2_MID': 9, 'P2_OUT': -9}

tm = pd.read_csv(f'{RD}/dream_timing_scans.csv')

# ------------------------------------------------- T1: sigma vs mesh voltage
fig, ax = plt.subplots(figsize=(7.6, 5.2))
g0 = tm[(tm.run == 'drift_mesh_scan_1') & (tm.axis == 'mesh')]
for det in DETS:
    mv = g0[f'mesh_v_{det}'].astype(float)
    sg = g0[f'{det}_sigma'].astype(float)
    ok = np.isfinite(mv) & np.isfinite(sg)
    o = np.argsort(mv[ok].to_numpy())
    x, y = mv[ok].to_numpy()[o], sg[ok].to_numpy()[o]
    ax.plot(x, y, color=st.DET_COLOR[det], marker='o', lw=2, label=det)
    st.direct_label(ax, x[-1], y[-1], det, st.DET_COLOR[det], dy=DY[det])
ax.axhline(20, color=st.C_RED, lw=1.2, ls=':')
ax.annotate('P2 goal: 20 ns', xy=(0.02, 20), xycoords=('axes fraction', 'data'),
            xytext=(0, 5), textcoords='offset points', fontsize=9,
            color=st.C_RED)
ax.set_xlabel('mesh voltage [V]')
ax.set_ylabel('single-station time resolution $\\sigma$ [ns]')
ax.set_title('Timing vs mesh HV (waveform TOA, walk-corrected) — '
             'drift_mesh_scan_1')
ax.legend(loc='upper right')
st.finish(fig, f'{OUT}/timing_vs_mesh.png')

# ------------------------- T2: sigma vs drift voltage + Magboltz gas overlay
gm = pd.read_csv(f'{RD}/gas_timing_model.csv')
g0 = tm[(tm.run == 'drift_mesh_scan_1') & (tm.axis == 'drift')]
fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13.2, 5.2))
for det in ['P2_MID', 'P2_OUT']:
    mv = g0[f'mesh_v_{det}'].astype(float)
    dv = g0[f'drift_v_{det}'].astype(float)
    sg = g0[f'{det}_sigma'].astype(float)
    ok = np.isfinite(dv) & np.isfinite(sg) & (dv > mv)   # drop zero-field point
    o = np.argsort(dv[ok].to_numpy())
    x, y = dv[ok].to_numpy()[o], sg[ok].to_numpy()[o]
    ax.plot(x, y, color=st.DET_COLOR[det], marker='o', lw=2, label=det)
    st.direct_label(ax, x[-1], y[-1], det, st.DET_COLOR[det], dy=DY[det])
ax.axhline(20, color=st.C_RED, lw=1.2, ls=':')
ax.set_xlabel('drift voltage [V] (mesh 450 V)')
ax.set_ylabel('single-station $\\sigma$ [ns]')
ax.set_yscale('log')
ax.set_yticks([10, 20, 50, 100, 150])
ax.set_yticklabels(['10', '20', '50', '100', '150'])
ax.set_title('Measured: $\\sigma$ vs drift HV (Ar/CO$_2$/iso)')
ax.legend()

# right: measured MID (gap = drift-450) with gas-model curves, sigma_0 fitted
det = 'P2_MID'
mv = g0[f'mesh_v_{det}'].astype(float)
dv = g0[f'drift_v_{det}'].astype(float)
sg = g0[f'{det}_sigma'].astype(float)
ok = np.isfinite(dv) & np.isfinite(sg) & (dv > mv)
gap = (dv[ok] - mv[ok]).to_numpy()
sig = sg[ok].to_numpy()
mA = gm[gm.key == 'ar_co2_iso_93_5_2']
mB = gm[gm.key == 'ar_cf4_iso_88_10_2']
sgasA = np.interp(gap, mA.dV, mA.sigma_gas_ns)
# electronics/walk floor from the plateau points (gap >= 250)
plateau = gap >= 250
s0 = float(np.sqrt(np.clip(np.median(sig[plateau] ** 2 - sgasA[plateau] ** 2),
                           0, None)))
gx = np.arange(60, 460, 5.0)
predA = np.sqrt(s0 ** 2 + np.interp(gx, mA.dV, mA.sigma_gas_ns) ** 2)
predB = np.sqrt(s0 ** 2 + np.interp(gx, mB.dV, mB.sigma_gas_ns) ** 2)
ax2.plot(gap, sig, 'o', color=st.C_MID, ms=7, label='measured (P2_MID)')
ax2.plot(gx, predA, '-', color=st.C_MID, lw=2,
         label=f'Magboltz Ar/CO$_2$/iso $\\oplus$ {s0:.1f} ns floor')
ax2.plot(gx, predB, '--', color=st.C_VIO, lw=2,
         label='same floor, Ar/CF$_4$/iso 88/10/2 (prediction)')
ax2.plot(gx, np.interp(gx, mB.dV, mB.sigma_gas_ns), ':', color=st.C_VIO, lw=1.5,
         label='Ar/CF$_4$/iso gas-only floor')
ax2.axhline(20, color=st.C_RED, lw=1.2, ls=':')
ax2.set_xlabel('drift $-$ mesh voltage [V] (3 mm gap)')
ax2.set_ylabel('single-station $\\sigma$ [ns]')
ax2.set_ylim(0, 45)
ax2.set_title('Magboltz expectation and the CF$_4$ gas gain')
ax2.legend(fontsize=9)
st.finish(fig, f'{OUT}/timing_vs_drift_magboltz.png')
print(f'  fitted non-gas floor sigma_0 = {s0:.2f} ns')

# ------------------------------------------ T3: 2D sigma(mesh,gap) surface
ts = pd.read_csv(f'{RD}/dream_timing_persubrun.csv')
hv = pd.read_csv(f'{S}/eos_inventory/hv_setpoints.csv').rename(
    columns={'mesh_or_resist': 'mesh_v', 'drift': 'drift_v'})
run2d = 'drift_mesh_2d_2'
g = ts[(ts.run == run2d) & (ts.kind == 'station')]
fig, axes = plt.subplots(1, 3, figsize=(14, 4.6))
im = None
for ax, det in zip(axes, DETS):
    gg = g[g.det == det].merge(
        hv[hv.det == det][['run', 'sub_run', 'mesh_v', 'drift_v']],
        on=['run', 'sub_run'], how='left').dropna(subset=['mesh_v', 'drift_v'])
    gg['mesh_v'] = gg.mesh_v.astype(float)
    gg['gap'] = gg.drift_v.astype(float) - gg.mesh_v
    gg = gg[np.isfinite(gg.sigma_ns) & (gg.sigma_ns < 200)]
    piv = gg.pivot_table(index='gap', columns='mesh_v', values='sigma_ns',
                         aggfunc='median')
    im = ax.pcolormesh(piv.columns, piv.index, piv.values, cmap=st.SEQ_CMAP_T,
                       vmin=10, vmax=60, edgecolors=st.SURFACE, linewidth=1.5)
    for yy in piv.index:
        for xx in piv.columns:
            v = piv.loc[yy, xx]
            if np.isfinite(v):
                ax.text(xx, yy, f'{v:.0f}', ha='center', va='center',
                        fontsize=7.2,
                        color='white' if v > 40 else st.TEXT)
    ax.set_title(det, color=st.DET_COLOR[det])
    ax.set_xlabel('mesh voltage [V]')
    ax.grid(False)
axes[0].set_ylabel('drift $-$ mesh voltage [V]')
fig.colorbar(im, ax=axes, label='$\\sigma$ [ns] (walk-corrected)',
             fraction=0.02)
fig.suptitle(f'Time resolution over the (mesh, drift) plane — {run2d}',
             fontweight='bold')
st.finish(fig, f'{OUT}/timing_2d_{run2d}.png', tight=False)
print('wrote', f'{OUT}/timing_2d_{run2d}.png')

# ------------------------- T3b: the same surface as curves, not boxes
fig, axes = plt.subplots(1, 3, figsize=(14, 4.8), sharey=True)
for ax, det in zip(axes, DETS):
    gg = g[g.det == det].merge(
        hv[hv.det == det][['run', 'sub_run', 'mesh_v', 'drift_v']],
        on=['run', 'sub_run'], how='left').dropna(subset=['mesh_v', 'drift_v'])
    gg['mesh_v'] = gg.mesh_v.astype(float)
    gg['gap'] = gg.drift_v.astype(float) - gg.mesh_v
    gg = gg[np.isfinite(gg.sigma_ns) & (gg.sigma_ns < 200)]
    gaps = sorted(gg.gap.unique())
    cmap = plt.get_cmap('viridis')
    for i, gap in enumerate(gaps):
        s = gg[gg.gap == gap].groupby('mesh_v', as_index=False).sigma_ns.median()
        s = s.sort_values('mesh_v')
        col = cmap(i / max(len(gaps) - 1, 1))
        ax.plot(s.mesh_v, s.sigma_ns, marker='o', ms=4.5, lw=1.8, color=col,
                label=f'{gap:.0f} V')
    ax.axhline(20, color=st.C_RED, lw=1.2, ls=':')
    ax.set_ylim(12, 70)          # P2_IN below 370 V runs to 180 ns, off scale
    ax.set_title(det, color=st.DET_COLOR[det])
    ax.set_xlabel('mesh voltage [V]')
axes[0].set_ylabel('single-station $\\sigma$ [ns]')
axes[0].annotate('P2 goal: 20 ns', xy=(0.03, 20), xycoords=('axes fraction', 'data'),
                 xytext=(0, 4), textcoords='offset points', fontsize=9,
                 color=st.C_RED)
axes[-1].legend(title='drift $-$ mesh [V]', fontsize=8, title_fontsize=8,
                loc='upper right', ncol=2)
fig.suptitle('Time resolution over the (mesh, drift) plane, as curves — '
             f'{run2d}\nmore gain always helps, and a 150 V gap costs ~15 ns; '
             'above a 300 V gap the curves converge\n'
             '(P2_IN below 370 V reaches 180 ns and is off scale)',
             fontweight='bold')
st.finish(fig, f'{OUT}/timing_2d_curves_{run2d}.png')

# ------------------- T3c: the same surface transposed -- sigma vs drift gap
# T3b asks "what does more gain buy"; this one asks "what does the drift field
# buy", which is the axis the gas change acts on.  Mesh 330 and 350 V are left
# out: at that gain P2_IN runs to 180 ns and the walk correction is not
# meaningful, so they would only stretch the y axis.
MESH_MIN = 370
fig, axes = plt.subplots(1, 3, figsize=(14, 4.8), sharey=True)
for ax, det in zip(axes, DETS):
    gg = g[g.det == det].merge(
        hv[hv.det == det][['run', 'sub_run', 'mesh_v', 'drift_v']],
        on=['run', 'sub_run'], how='left').dropna(subset=['mesh_v', 'drift_v'])
    gg['mesh_v'] = gg.mesh_v.astype(float)
    gg['gap'] = gg.drift_v.astype(float) - gg.mesh_v
    gg = gg[np.isfinite(gg.sigma_ns) & (gg.sigma_ns < 200)
            & (gg.mesh_v >= MESH_MIN)]
    meshes = sorted(gg.mesh_v.unique())
    cmap = plt.get_cmap('plasma')
    for i, mv in enumerate(meshes):
        s = gg[gg.mesh_v == mv].groupby('gap', as_index=False).sigma_ns.median()
        s = s.sort_values('gap')
        col = cmap(0.08 + 0.78 * i / max(len(meshes) - 1, 1))
        ax.plot(s.gap, s.sigma_ns, marker='o', ms=4.5, lw=1.8, color=col,
                label=f'{mv:.0f} V')
    ax.axhline(20, color=st.C_RED, lw=1.2, ls=':')
    ax.axvline(250, color=st.GRID, lw=1.4, ls='--')
    ax.set_ylim(12, 46)
    ax.set_title(det, color=st.DET_COLOR[det])
    ax.set_xlabel('drift $-$ mesh voltage [V]')
    # the physical axis: 3 mm gap, so 150 V = 500 V/cm
    top = ax.secondary_xaxis('top', functions=(lambda v: v / 0.3,
                                               lambda e: e * 0.3))
    top.set_xlabel('drift field [V/cm]', fontsize=9.5)
    top.tick_params(labelsize=8.5)
axes[0].set_ylabel('single-station $\\sigma$ [ns]')
axes[0].annotate('P2 goal: 20 ns', xy=(0.03, 20), xycoords=('axes fraction', 'data'),
                 xytext=(0, 4), textcoords='offset points', fontsize=9,
                 color=st.C_RED)
axes[0].annotate('working point', xy=(250, 44), fontsize=8.5, color='#666',
                 ha='center', va='top')
axes[-1].legend(title='mesh [V]', fontsize=8, title_fontsize=8,
                loc='upper right', ncol=2)
fig.suptitle('Time resolution vs drift voltage, one curve per mesh setting — '
             f'{run2d}\nat mesh 450 V a 150 V gap costs 8–13 ns; going from the '
             '250 V working point to 350–450 V would still buy 1.0 ns on MID, '
             '3.5 on OUT, 5.5 on IN\n'
             f'(mesh $\\geq$ {MESH_MIN} V only — below that the stations are '
             'amplitude-limited and run off scale)',
             fontweight='bold')
st.finish(fig, f'{OUT}/timing_2d_curves_vs_drift_{run2d}.png')

# ------------- T3d: what the drift choice is worth, as a table figure
# The numbers a reviewer writes down off T3c, at the working-point gain
# (mesh 450 V): what a too-small gap costs, where we ran, and what the scan's
# best point would have given.
rows = []
for det in DETS:
    gg = g[g.det == det].merge(
        hv[hv.det == det][['run', 'sub_run', 'mesh_v', 'drift_v']],
        on=['run', 'sub_run'], how='left').dropna(subset=['mesh_v', 'drift_v'])
    gg['mesh_v'] = gg.mesh_v.astype(float)
    gg['gap'] = gg.drift_v.astype(float) - gg.mesh_v
    gg = gg[np.isfinite(gg.sigma_ns) & (gg.sigma_ns < 200)
            & (gg.mesh_v == 450)]
    s = gg.groupby('gap', as_index=False).sigma_ns.median().sort_values('gap')
    if not len(s):
        continue
    best = s.loc[s.sigma_ns.idxmin()]

    def at(gap):
        r = s[s.gap == gap]
        return r.sigma_ns.iloc[0] if len(r) else np.nan

    rows.append(dict(det=det, s150=at(150), s250=at(250),
                     sbest=best.sigma_ns, gbest=best.gap))

fig, ax = plt.subplots(figsize=(9.6, 2.6))
ax.axis('off')
cols = ['station', '$\\sigma$ at 150 V gap\n(500 V/cm)',
        '$\\sigma$ at 250 V gap\n(833 V/cm, working point)',
        'best in scan']
cells, colours = [], []
for r in rows:
    cells.append([r['det'], f"{r['s150']:.1f} ns", f"{r['s250']:.1f} ns",
                  f"{r['sbest']:.1f} ns  at {r['gbest']:.0f} V"])
    colours.append([st.DET_COLOR[r['det']]] + ['#000'] * 3)
tab = ax.table(cellText=cells, colLabels=cols, cellLoc='center',
               loc='center', colWidths=[0.17, 0.26, 0.33, 0.26])
tab.auto_set_font_size(False)
tab.set_fontsize(11)
tab.scale(1, 2.15)
for (r, c), cell in tab.get_celld().items():
    cell.set_edgecolor(st.GRID)
    if r == 0:
        cell.set_facecolor('#f1f0ec')
        cell.set_text_props(fontweight='bold', fontsize=9.8)
    else:
        row = rows[r - 1]
        if c == 0:
            cell.set_text_props(color=st.DET_COLOR[row['det']],
                                fontweight='bold')
        # highlight the column that beats the 20 ns goal
        val = [None, row['s150'], row['s250'], row['sbest']][c]
        if c and val is not None and val < 20:
            cell.set_facecolor('#e6f2e2')
ax.set_title('What the drift field is worth at the working-point gain — '
             f'{run2d}, mesh 450 V\n'
             'walk-corrected waveform TOA vs the trigger; the 20 ns P2 goal is '
             'shaded green\n'
             'the campaign ran at a 250 V gap: 1.0 ns from optimal on MID, '
             '3.5 on OUT, 5.5 on IN', fontsize=11.5, pad=12)
st.finish(fig, f'{OUT}/timing_drift_choice_table_{run2d}.png', tight=False)

# -------------------------------- T4: correction ladder at the working point
f = f'{A}/telescope/highstat_eff_1/beam_commissioning_00/29_waveform_timing/waveform_timing_summary.json'
d = json.load(open(f))
fig, ax = plt.subplots(figsize=(8.2, 5.0))
steps = ['sigma_raw_ns', 'sigma_ftst_ns', 'sigma_walk_ns']
labels = ['raw TOA', '+ ftst clock-phase', '+ time-walk']
xpos = np.arange(len(steps))
for st_ in d['stations']:
    det = st_['detector']
    alg = st_['algorithms'][st_['best_algorithm']]
    y = [alg[k] for k in steps]
    ax.plot(xpos, y, color=st.DET_COLOR[det], marker='o', lw=2,
            label=f"{det} ({st_['best_algorithm']})")
    st.direct_label(ax, xpos[-1], y[-1], f"{y[-1]:.1f} ns",
                    st.DET_COLOR[det], dy=DY[det])
pair = {p['pair']: p for p in d['pairs']}
ax.set_xticks(xpos)
ax.set_xticklabels(labels)
ax.set_ylabel('$\\sigma$ vs trigger [ns]')
ax.set_title('Timing correction ladder at the working point '
             '(highstat_eff_1, mesh 450 / drift 700)')
ax.legend(loc='upper right')
st.finish(fig, f'{OUT}/timing_ladder.png')

# pair-based single-station numbers for the report
print('pair-derived single-station sigma (ns):',
      {k: v['sigma_single_ns'] for k, v in pair.items()})
print('done')
