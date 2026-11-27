#!/usr/bin/env python3
"""joint_gap_water.py -- can the drift scan tell the water fraction from the gap?

The drift scan alone (median peak time - zero-field time vs dV) constrains a COMBINATION of the gap d and
the water fraction f. The drift time at the operating point, measured independently of the scan and of the
pulse-shape response (pad-pair micro-TPC of stage 48: 118 +- 3 ns with the transverse diffusion of wet gas,
129 +- 3 ns with dry-gas diffusion; mean-pulse fit of stage 47: 125 ns), is
T(dV = 300 V) = d / v_d(f, 300 V / d): a second combination with the opposite slope. This script profiles
both on a (f, d) grid.

  scan chi2(f, d): response scale s and time offset dt0 free (dt0 with its 5 ns prior), 2 ns per-point syst.
  T constraint   : ((T(300; f, d) - T_OP) / E_OP)^2, T_OP = 120 ns, E_OP = 6 ns (covers the estimates)

Out: results/joint_gap_water.csv, results/joint_gap_water.png, printed summary.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import analyze_water as aw

T_OP, E_OP = 120.0, 6.0
gas = aw.load_gases(); resp = aw.Response(); data = aw.load_data()
tabs = sorted([(f, gas[k]) for k, f in aw.FAMILIES['h2o'] if k in gas], key=lambda a: a[0])
F = np.arange(0.0, 2.01, 0.05); Dg = np.arange(3.0, 6.01, 0.1)
rows = []
for name in ('chamber', 'reg24_good'):
    for f in F:
        gm = aw.GasModel.mix(tabs, f)
        for d in Dg:
            r = aw.fit_one(gm, resp, data[name], d_free=False, s_free=True, sys_ns=2.0, d0=d, dvmin=20.0)
            T = (d / 10.0) / gm.v(300.0 / (d / 10.0))
            rows.append(dict(data=name, f=f, d=d, chi2_scan=r['chi2'], s=r['s'], dt0=r['dt0'], T_op=float(T),
                             chi2_T=((T - T_OP) / E_OP) ** 2))
R = pd.DataFrame(rows); R['chi2_joint'] = R.chi2_scan + R.chi2_T
R.to_csv(os.path.join(aw.RES, 'joint_gap_water.csv'), index=False)
fig, axs = plt.subplots(1, 2, figsize=(14, 5.6))
for ax, name in zip(axs, ('chamber', 'reg24_good')):
    q = R[R.data == name]
    for col, c, lab in (('chi2_scan', '#257de0', 'drift scan'), ('chi2_T', '#f47230', f'T(300 V) = {T_OP:g} +- {E_OP:g} ns'),
                        ('chi2_joint', 'k', 'both')):
        Z = q.pivot(index='d', columns='f', values=col)
        z = Z.to_numpy() - Z.to_numpy().min()
        ax.contour(Z.columns, Z.index, z, levels=[2.3, 6.2], colors=c, linestyles=['-', '--'])
        ax.plot([], [], color=c, label=lab)
    b = q.loc[q.chi2_joint.idxmin()]
    ax.plot(b.f, b.d, 'k*', ms=12)
    ax.set(xlabel='water fraction [%]', ylabel='drift gap of the nominal zone [mm]', title=f'{name}: 68 % (solid) and 95 % (dashed) regions')
    ax.axhline(4.0, color='#9c9c9c', lw=1); ax.grid(color='#e6e6e6'); ax.legend(frameon=False)
    s_ = q.loc[q.chi2_scan.idxmin()]
    print(f'{name}: scan alone best f {s_.f:.2f} % d {s_.d:.1f} mm chi2 {s_.chi2_scan:.1f} (s {s_.s:.2f}, dt0 {s_.dt0:+.1f}); '
          f'joint best f {b.f:.2f} % d {b.d:.1f} mm (scan chi2 {b.chi2_scan:.1f}, T_op {b.T_op:.0f} ns, s {b.s:.2f}, dt0 {b.dt0:+.1f})')
    ok = q[q.chi2_joint < q.chi2_joint.min() + 2.3]
    print(f'   joint 68 % region: f {ok.f.min():.2f}-{ok.f.max():.2f} %, d {ok.d.min():.1f}-{ok.d.max():.1f} mm')
    oks = q[q.chi2_scan < q.chi2_scan.min() + 2.3]
    print(f'   scan-alone 68 % region: f {oks.f.min():.2f}-{oks.f.max():.2f} %, d {oks.d.min():.1f}-{oks.d.max():.1f} mm')
    for d0 in (3.5, 4.0, 4.5, 5.0):
        qq = q[np.isclose(q.d, d0)]; m_ = qq.loc[qq.chi2_scan.idxmin()]
        print(f'   d = {d0:.1f} mm: scan prefers f = {m_.f:.2f} % (chi2 {m_.chi2_scan:.1f}, s {m_.s:.2f}); there T(300 V) = {m_.T_op:.0f} ns')
    q0 = q[np.isclose(q.f, 0.0)]; m0 = q0.loc[q0.chi2_scan.idxmin()]
    print(f'   dry gas, any gap in 3-6 mm: best chi2 {m0.chi2_scan:.0f} at d = {m0.d:.1f} mm')
fig.tight_layout(); fig.savefig(os.path.join(aw.RES, 'joint_gap_water.png'), dpi=110)
