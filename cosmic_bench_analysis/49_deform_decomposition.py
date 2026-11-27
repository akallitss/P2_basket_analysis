#!/usr/bin/env python3
"""49_deform_decomposition.py -- the deformation of a P2 chamber, put together.

Combines the deformation stages into one per-pad table and a few figures:

  DRIFT GAP    T_d per pad = mean-pulse fit (stage 47, fine structure) with its
               large-scale level tied to the pulse-to-pulse spread (stage 47,
               single-pad events), which the ion tail and the electronics cannot
               move:  T_d = T_fit - <T_fit - T_spread>_(R = --s-radius, 45 mm).
               If the gain varies on a shorter scale than that (det1), the
               mean-pulse fit follows the ion tail, not the gap: use --s-radius 20.
               gap = solution of  gap / v_d((V_drift - V_mesh) / gap) = T_d.
  ION TAIL     S = <T_fit - T_spread>_(45 mm): the part of the mean-pulse width
               that is NOT drift time (single-electron response: ion transit in
               the amplification gap, or electronics shaping).
  GAIN         track-matched amplitude (pad under the track, central 4 mm,
               x cos theta) divided by the three things the drift gap does to it:
               primary charge (proportional to the gap), the field dependence of
               the collection/transparency (measured: chamber amplitude vs drift
               voltage, stage 44) at the LOCAL field, and the ballistic deficit
               of a longer arrival time (forward model). What is left is the
               intrinsic amplification gain G.
  AMP. GAP     if all of the variation of G is the amplification gap:
               d(gap) = - dlnG / k,  k = -(alpha - E dalpha/dE) from Magboltz at the
               mesh voltage (printed). The test of that 'if' is the ion tail: a
               larger gap must also give a longer ion transit, so S must fall
               with ln G; the slope is reported.
  SURFACES     with the track-geometry height of stage 45 (z_eff = mesh height +
               half the gap, M3 frame): mesh z = z_eff - gap/2, window z = z_eff +
               gap/2, per region.
  EFFECTS      per drift-gap class: time resolution, amplitude, gain, hit
               fraction, no-hit rate by angle, pillar-shadow depth (stage-41
               lattice).

Cross-checks printed: gap per region against the drift-scan estimates of stage
44 (timing lambda and amplitude-curve stretch, another day) and the micro-TPC
drift times of stage 48.

Out (<run>/49_deform_decomposition/): pads.csv, regions.csv, regions_all.csv, classes.csv,
deformation_maps.png, consistency.png, effects.png, ridge_profile.png (+ csv), summary.txt
Usage: python3 49_deform_decomposition.py [run_key] [--scan-dir DIR] [--gas CSV]
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2_qa_config as qa        # noqa: E402
import p2_deform as pdf          # noqa: E402
import p2_deform_model as dm     # noqa: E402
import p2_study_cache as sc      # noqa: E402
import p2_study_common as cm     # noqa: E402

RNG = np.random.default_rng(20261128)
GAP0 = 4.0
CLASS_EDGES = [0, 3.6, 4.4, 5.5, 30]
CLASS_NAMES = ['< 3.6', '3.6-4.4', '4.4-5.5', '> 5.5']
SER_PAR = dict(T_ion=278.0, t_peak=174.0, f_e=0.6)      # model SER fitted to the reference mean pulse


def make_vd(gas_csv=None):
    """v_d(E) [mm/ns]; default the dry Magboltz table of p2_deform_model."""
    if gas_csv is None:
        return dm.vdrift, 'Magboltz dry Ar/iC4H10 95/5 (garfield_inputs/gas_table.csv)'
    g = pd.read_csv(gas_csv)
    ecol = [c for c in g.columns if c.lower().startswith('e')][0]
    vcol = [c for c in g.columns if c.lower().startswith('v')][0]
    x = np.r_[0.0, g[ecol].to_numpy(float)]; y = np.r_[0.0, g[vcol].to_numpy(float)]
    if np.nanmax(y) < 0.01:           # cm/ns -> mm/ns
        y = y * 10.0
    return (lambda E: np.interp(E, x, y)), os.path.basename(gas_csv)


def gap_of_T(T, dV, vd):
    if not np.isfinite(T):
        return np.nan
    f = lambda d: d / float(vd(dV / (d / 10.0))) - T
    try:
        return brentq(f, 1.5, 12.0)
    except ValueError:
        return np.nan


def ballistic(Ts):
    """Peak of box(T) (x) SER per unit charge, relative to T = 125 ns (model SER)."""
    ts, s = dm.ser(dt=2.0, **SER_PAR)
    out = []
    for T in np.atleast_1d(Ts):
        n = max(int(round(T / 2.0)), 1)
        out.append(np.convolve(s, np.ones(n) / n)[:len(s)].max())
    out = np.array(out)
    n0 = int(round(125 / 2.0))
    return out / np.convolve(s, np.ones(n0) / n0)[:len(s)].max()


def slide(P, col, wcol, R):
    """Weighted mean of col over the pads within R mm."""
    xy = P[['pad_cx', 'pad_cy']].to_numpy()
    v = P[col].to_numpy(float); w = P[wcol].fillna(0).to_numpy(float) * np.isfinite(v)
    v = np.nan_to_num(v)
    out = np.full(len(P), np.nan)
    for i, nb in enumerate(cKDTree(xy).query_ball_point(xy, R)):
        ws = w[nb].sum()
        if ws > 0:
            out[i] = (w[nb] * v[nb]).sum() / ws
    return out


def wls(x, y, e):
    w = 1 / np.asarray(e, float) ** 2
    X = np.c_[np.ones(len(x)), x]
    C = np.linalg.inv(X.T @ (w[:, None] * X)); b = C @ X.T @ (w * y)
    chi = float((w * (y - X @ b) ** 2).sum())
    return b, np.sqrt(np.diag(C)), chi


def smap(ax, P, c, title, cmap='coolwarm', lo=None, hi=None, sym=None, label=''):
    v = P[c].to_numpy(float)
    if lo is None:
        lo, hi = np.nanpercentile(v, [2, 98])
    if sym is not None:
        a = max(abs(lo - sym), abs(hi - sym)); lo, hi = sym - a, sym + a
    sc_ = ax.scatter(P.pad_cx, P.pad_cy, c=v, s=20, marker='s', cmap=cmap, vmin=lo, vmax=hi)
    ax.set_aspect('equal'); ax.set_title(title, fontsize=11)
    ax.set_xlabel('pad x [mm]'); ax.set_ylabel('pad y [mm]')
    plt.colorbar(sc_, ax=ax, fraction=0.04, label=label)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run_key', nargs='?', default='det2_long_repos')
    ap.add_argument('--scan-dir', default=None, help='stage-44 output folder (drift scan of the same chamber)')
    ap.add_argument('--gas', default=None, help='CSV with E [V/cm] and v_d columns; default dry Magboltz')
    ap.add_argument('--dv', type=float, default=None, help='V_drift - V_mesh [V]; default from the sub-run name')
    ap.add_argument('--s-radius', type=float, default=45.0,
                    help='smoothing radius [mm] of the ion-tail term S; use ~20 when the gain varies on a short scale (det1)')
    a = ap.parse_args()
    cfg = qa.get_config(a.run_key)
    out = cfg.out_dir('49_deform_decomposition')
    d47 = cfg.out_dir('47_deform_pulse_maps')
    L = []
    vd, gas_name = make_vd(a.gas)
    if a.dv is None:
        nums = [int(x) for x in cfg.SUB_RUN.split('_') if x.isdigit()]
        a.dv = float(abs(nums[-1] - nums[-2]))
    E0 = a.dv / (GAP0 / 10.0)
    L.append(f'{a.run_key} ({cfg.DET_NAME}): V_drift - V_mesh = {a.dv:g} V (nominal {E0:.0f} V/cm for {GAP0:g} mm); '
             f'drift-velocity table: {gas_name}; T_d({GAP0:g} mm) = {GAP0/float(vd(E0)):.1f} ns at the nominal field')

    P = pd.read_csv(f'{d47}/pads.csv')
    # ---- A. drift time and ion-tail term ------------------------------------------------
    P['dser'] = P.sl_T_fit - P.sl_T_s
    P['S'] = slide(P, 'dser', 'sl_n1', a.s_radius)
    P['T_d'] = P.sl_T_fit - P.S
    P['gap_dry'] = [gap_of_T(t, a.dv, vd) for t in P.T_d]
    # ---- drift velocity measured by the chamber itself (drift scan), relative to the nominal field ----
    # In the regions whose gap is nominal, (median peak time - zero-drift time) is proportional to
    # T_d(E) = GAP0 / v_d(E), so its ratio to the value at the nominal field gives v_d(E) / v_d(E0)
    # with no gas model. Only the anchor v_d(E0) is taken from the table (checked by the micro-TPC).
    if a.scan_dir and os.path.isfile(f'{a.scan_dir}/region_stats.csv') and a.gas is None:
        rs = pd.read_csv(f'{a.scan_dir}/region_stats.csv')
        g24 = P.groupby('reg24').gap_dry.mean()
        good = g24[(g24 - GAP0).abs() < 0.25].index.tolist()
        d_ = rs[(rs.scheme == 'reg24') & (rs.kind == 'drift') & rs.region.isin(good)]
        w_ = 1 / d_.tom_med_e_300 ** 2
        cur = (d_.assign(wt=w_ * d_.tom_med_300, w=w_).groupby('dV')[['wt', 'w']].sum())
        cur['tom'] = cur.wt / cur.w
        cur['E'] = cur.index / (GAP0 / 10.0)
        t0 = float(cur.tom[cur.index < 5].iloc[0])
        cur = cur[cur.index > 90]
        tn = float(np.interp(a.dv, cur.index, cur.tom)) - t0
        cur['T_rel'] = (cur.tom - t0) / tn
        v0 = float(vd(E0))
        ex, vy = cur.E.to_numpy(), v0 / cur.T_rel.to_numpy()
        vd_dry = vd
        vd = lambda E: np.interp(E, ex, vy)
        gas_name = (f'measured in the chamber: v_d(E)/v_d({E0:.0f} V/cm) from the drift scan in the nominal-gap regions '
                    f'{good} (peak time minus the zero-drift time {t0:.0f} ns), anchored to the dry-gas table at {E0:.0f} V/cm')
        L.append('drift velocity used for the gap: ' + gas_name)
        L.append('   E [V/cm]: ' + ' '.join(f'{e:6.0f}' for e in ex))
        L.append('   measured v_d [um/ns]: ' + ' '.join(f'{1e3*v:6.1f}' for v in vy))
        L.append('   dry Magboltz [um/ns]: ' + ' '.join(f'{1e3*float(vd_dry(e)):6.1f}' for e in ex))
        pd.DataFrame({'E_Vcm': ex, 'vd_measured_mm_ns': vy, 'vd_dry_magboltz_mm_ns': [float(vd_dry(e)) for e in ex]}).to_csv(
            f'{out}/drift_velocity_measured.csv', index=False)
    P['gap'] = [gap_of_T(t, a.dv, vd) for t in P.T_d]
    L.append(f'   gap with the dry-gas table instead: median {P.gap_dry.median():.2f}, 95 % {P.gap_dry.quantile(.95):.2f}, max {P.gap_dry.max():.2f} mm')
    P['gap_fit_only'] = [gap_of_T(t, a.dv, vd) for t in P.sl_T_fit]
    P['gap_spread_only'] = [gap_of_T(t, a.dv, vd) for t in P.sl_T_s]
    P['E_loc'] = a.dv / (P.gap / 10.0)
    P['gclass'] = pd.cut(P.gap, CLASS_EDGES, labels=False)
    ok = P.gap.notna()
    L.append(f'drift gap per pad (sliding 18 mm, {ok.sum()} pads): median {P.gap.median():.2f} mm, '
             f'min {P.gap.min():.2f}, 5 % {P.gap.quantile(.05):.2f}, 95 % {P.gap.quantile(.95):.2f}, max {P.gap.max():.2f} mm; '
             f'pads per class {dict(zip(CLASS_NAMES, P.gclass.value_counts().sort_index().tolist()))}')
    L.append(f'   same from the mean-pulse fit alone: 5-95 % {P.gap_fit_only.quantile(.05):.2f}-{P.gap_fit_only.quantile(.95):.2f}, '
             f'max {P.gap_fit_only.max():.2f}; from the spread alone: {P.gap_spread_only.quantile(.05):.2f}-'
             f'{P.gap_spread_only.quantile(.95):.2f}, max {P.gap_spread_only.max():.2f} mm')
    L.append(f'ion-tail term S = <T_fit - T_spread> over {a.s_radius:g} mm: median {P.S.median():+.1f} ns, 5-95 % {P.S.quantile(.05):+.1f} .. '
             f'{P.S.quantile(.95):+.1f} ns; fitted extra delay dt0 (sliding): 5-95 % {P.sl_dt0.quantile(.05):+.1f} .. '
             f'{P.sl_dt0.quantile(.95):+.1f} ns')

    # ---- C. gain ----------------------------------------------------------------------------
    R, WF, _ = pdf.load(a.run_key, with_features=True)
    tk = R[R.good & R.near_track & (R.d_trk < 4.0) & (R.sd >= 2)].copy()
    tk['ampc'] = tk.amp * np.cos(np.radians(tk.theta))
    hit = tk[tk.is_hit]
    by = {int(c): g.ampc.to_numpy() for c, g in hit.groupby('channel_id')}
    xy = P[['pad_cx', 'pad_cy']].to_numpy(); ids = P.channel_id.to_numpy()
    amp_sl, n_sl = np.full(len(P), np.nan), np.zeros(len(P), int)
    for i, nb in enumerate(cKDTree(xy).query_ball_point(xy, 18.0)):
        v = [by[int(c)] for c in ids[nb] if int(c) in by]
        if v:
            v = np.concatenate(v); n_sl[i] = len(v)
            if len(v) >= 30:
                amp_sl[i] = np.median(v)
    P['amp_sl'], P['n_amp_sl'] = amp_sl, n_sl
    # field dependence of the collection: chamber median amplitude vs drift field (stage 44)
    F_rel = lambda E: np.ones_like(np.asarray(E, float))
    f_note = 'no drift scan given: field dependence of the collection NOT corrected'
    if a.scan_dir and os.path.isfile(f'{a.scan_dir}/region_stats.csv'):
        rs = pd.read_csv(f'{a.scan_dir}/region_stats.csv')
        c = rs[(rs.scheme == 'all') & (rs.kind == 'drift') & (rs.dV > 60)].sort_values('E')
        pc = np.polyfit(np.log(c.E), np.log(c.amp_med), 4, w=c.amp_med / c.amp_med_e)
        f0 = np.polyval(pc, np.log(E0))
        F_rel = lambda E: np.exp(np.polyval(pc, np.log(np.clip(E, c.E.min(), c.E.max()))) - f0)
        f_note = (f'collection vs field from the chamber amplitude curve of {os.path.basename(a.scan_dir)} '
                  f'(4th-order polynomial in ln E, {len(c)} points): F(450)/F({E0:.0f}) = {float(F_rel(450.)):.3f}, '
                  f'F(900) = {float(F_rel(900.)):.3f}, F(1200) = {float(F_rel(1200.)):.3f}')
    L.append(f_note)
    P['f_prim'] = P.gap / GAP0
    P['f_field'] = F_rel(P.E_loc.to_numpy())
    # ballistic deficit beyond what the chamber curve already contains at that field for the nominal gap
    T_nom_at_E = GAP0 / np.asarray(vd(P.E_loc.to_numpy()), float)
    P['f_ball'] = ballistic(P.T_d.fillna(125).to_numpy()) / ballistic(np.nan_to_num(T_nom_at_E, nan=125.0))
    P['lnG'] = np.log(P.amp_sl / (P.f_prim * P.f_field * P.f_ball))
    P['lnG'] -= np.nanmedian(P.lnG)
    P['lnA'] = np.log(P.amp_sl) - np.nanmedian(np.log(P.amp_sl))
    # Magboltz: d ln G / d gap at the mesh voltage
    k_note = 'no amplification table'
    kg = np.nan
    amp_csv = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'gas_studies', 'results', 'ar_iso_95_5.csv')
    nums = [int(x) for x in cfg.SUB_RUN.split('_') if x.isdigit()]
    v_mesh = float(min(nums[-2:]))
    if os.path.isfile(amp_csv):
        g = pd.read_csv(amp_csv); g = g[g.region == 'amp'].sort_values('E_Vcm')
        g_um = 128.0
        Ea = v_mesh / (g_um * 1e-4)
        al = np.interp(Ea, g.E_Vcm, g.alpha_cm)
        dal = np.gradient(g.alpha_cm.to_numpy(), g.E_Vcm.to_numpy())
        kg = -(al - Ea * np.interp(Ea, g.E_Vcm, dal)) * 1e-4          # per um
        k_note = (f'Magboltz (dry, no Penning): alpha = {al:.0f} /cm at {Ea/1e3:.1f} kV/cm ({v_mesh:g} V over {g_um:g} um), '
                  f'd ln G / d gap = {-kg:+.4f} per um -> 1 e-fold of gain = {1/kg:.0f} um of amplification gap')
    L.append(k_note)
    P['dgap_amp_um'] = -P.lnG / kg
    L.append(f'track-matched amplitude (sliding 18 mm): spread of ln A {pdf.s68(P.lnA):.3f} (robust sd); after removing '
             f'primary charge, field and ballistic terms the intrinsic gain has spread of ln G {pdf.s68(P.lnG):.3f} '
             f'-> amplification-gap rms {pdf.s68(P.dgap_amp_um):.1f} um if it is all gap '
             f'(5-95 % {P.dgap_amp_um.quantile(.05):+.0f} .. {P.dgap_amp_um.quantile(.95):+.0f} um)')
    L.append(f'   correlation ln A vs ln(gap): {P[["lnA"]].assign(g=np.log(P.gap)).corr().iloc[0, 1]:+.2f}; '
             f'ln G vs ln(gap): {P[["lnG"]].assign(g=np.log(P.gap)).corr().iloc[0, 1]:+.2f}')

    # ---- E. ion tail against gain -----------------------------------------------------------
    reg = pd.read_csv(f'{d47}/regions.csv')
    rows = []
    for K in (24, 48):
        for r, gq in P.groupby(f'reg{K}'):
            hh = hit[hit[f'reg{K}'] == r]
            rr = reg[(reg.K == K) & (reg.region == r)].iloc[0]
            gapr = gap_of_T(rr.T_s, a.dv, vd)
            T_nom = GAP0 / float(vd(a.dv / (gapr / 10.0)))
            fb = float(ballistic([rr.T_s])[0] / ballistic([T_nom])[0])
            lnG = np.log(hh.ampc.median() / ((gapr / GAP0) * float(F_rel(a.dv / (gapr / 10.0))) * fb))
            rows.append(dict(K=K, region=int(r), x=gq.pad_cx.mean(), y=gq.pad_cy.mean(), n=int(rr.n), n_trk=len(hh),
                             T_fit=rr.T_fit, e_T_fit=rr.e_T_fit, T_s=rr.T_s, e_T_s=rr.e_T_s, dt0=rr.dt0, e_dt0=rr.e_dt0,
                             S=rr.T_fit - rr.T_s, e_S=np.hypot(rr.e_T_fit, rr.e_T_s), gap_spread=gapr,
                             gap=float(np.nanmean(gq.gap)), amp=hh.ampc.median(),
                             e_lnA=1.25 * pdf.s68(np.log(hh.ampc)) / np.sqrt(max(len(hh), 1)), lnG=lnG,
                             conn=int(gq.connector_N.mode().iloc[0])))
    RG = pd.DataFrame(rows)
    for K in (24, 48):
        q = RG[RG.K == K].copy(); q['lnG'] -= q.lnG.median()
        RG.loc[RG.K == K, 'lnG'] = q.lnG
        b, e, chi = wls(q.lnG.to_numpy(), q.S.to_numpy(), q.e_S.to_numpy())
        r_ = np.corrcoef(q.lnG, q.S)[0, 1]
        sl_boot = [np.polyfit(q.lnG.to_numpy()[i], q.S.to_numpy()[i], 1)[0]
                   for i in (RNG.integers(0, len(q), len(q)) for _ in range(400))]
        L.append(f'ion tail vs gain, {K} regions: S = T_fit - T_spread against ln G: r = {r_:+.2f}, slope '
                 f'{np.polyfit(q.lnG, q.S, 1)[0]:+.1f} +- {np.std(sl_boot):.1f} ns per e-fold of gain (weighted {b[1]:+.1f}); '
                 f'S against ln A (no correction): r = {np.corrcoef(np.log(q.amp), q.S)[0, 1]:+.2f}')
        # within-connector (electronics constants removed)
        qs = q.S - q.groupby('conn').S.transform('mean'); qg = q.lnG - q.groupby('conn').lnG.transform('mean')
        if (qg.abs() > 0).sum() > 6:
            L.append(f'      within connectors: r = {np.corrcoef(qg, qs)[0, 1]:+.2f}, slope {np.polyfit(qg, qs, 1)[0]:+.1f} ns per e-fold')
    Pk = P.dropna(subset=['S', 'lnG'])
    L.append(f'   pad level (sliding values, {len(Pk)} pads): r(S, ln G) = {np.corrcoef(Pk.S, Pk.lnG)[0, 1]:+.2f}')
    ch = P.groupby(['connector_N']).agg(S=('S', 'median'), dt0=('sl_dt0', 'median'), lnG=('lnG', 'median'),
                                        gap=('gap', 'median'), rise=('sl_med_rise', 'median'), fall=('sl_med_fall', 'median'))
    L.append('per connector medians (electronics-sized terms):  conn   S[ns]  dt0[ns]  lnG   gap[mm]  rise  fall')
    for c_, r_ in ch.iterrows():
        L.append(f'      {int(c_):2d}  {r_.S:+6.1f}  {r_.dt0:+5.1f}  {r_.lnG:+.2f}  {r_.gap:5.2f}  {r_.rise:6.1f} {r_.fall:6.1f}')

    # ---- cross-checks against the drift scan (other day) and the micro-TPC ----------------
    q24 = RG[RG.K == 24].set_index('region')
    if a.scan_dir and os.path.isfile(f'{a.scan_dir}/gap_from_lambda_highfield.csv'):
        gl = pd.read_csv(f'{a.scan_dir}/gap_from_lambda_highfield.csv'); gl = gl[gl.scheme == 'reg24'].set_index('region')
        ac = pd.read_csv(f'{a.scan_dir}/amplitude_collapse.csv')
        ac = ac[(ac.scheme == 'reg24') & (ac.col == 'amp_med')]; ac = ac[ac.dVmin == ac.dVmin.min()].set_index('region')
        q24 = q24.join(gl[['d_mm', 'd_e']].rename(columns={'d_mm': 'gap_scan', 'd_e': 'e_gap_scan'}))
        q24 = q24.join((ac.s * GAP0).rename('gap_ampcurve')).join((ac.s_e_scaled * GAP0).rename('e_gap_ampcurve'))
        q24['lnG_scan'] = np.log(ac.g / ac.s); q24['lnG_scan'] -= q24.lnG_scan.median()
        for col, nm in (('gap_scan', 'drift-scan timing'), ('gap_ampcurve', 'amplitude-curve stretch')):
            r_ = np.corrcoef(q24.gap, q24[col])[0, 1]
            p_ = np.polyfit(q24.gap, q24[col], 1)
            L.append(f'gap per region (this stage, 24 regions) vs {nm} (stage 44, other day): r = {r_:+.2f}, '
                     f'{col} = {p_[0]:.2f} * gap {p_[1]:+.2f}, rms of the difference {np.std(q24[col] - q24.gap):.2f} mm')
        L.append(f'intrinsic gain per region vs the drift-scan gain (stage 44, g/s): r = {np.corrcoef(q24.lnG, q24.lnG_scan)[0, 1]:+.2f}; '
                 f'S vs scan gain: r = {np.corrcoef(q24.S, q24.lnG_scan)[0, 1]:+.2f}, slope '
                 f'{np.polyfit(q24.lnG_scan, q24.S, 1)[0]:+.1f} ns per e-fold')
    d48 = cfg.out_dir('48_deform_microtpc')
    if os.path.isfile(f'{d48}/drift_time_by_zone.csv') and os.path.isfile(f'{d48}/pairs.csv'):
        mt = pd.read_csv(f'{d48}/drift_time_by_zone.csv'); mt = mt[mt.est == 't_pk']
        pr = pd.read_csv(f'{d48}/pairs.csv')
        gmap = P.set_index('channel_id')
        pr['gapf'] = 0.5 * (gmap.gap.reindex(pr.channel_id_a).to_numpy() + gmap.gap.reindex(pr.channel_id_b).to_numpy())
        pr['Tf'] = 0.5 * (gmap.T_d.reindex(pr.channel_id_a).to_numpy() + gmap.T_d.reindex(pr.channel_id_b).to_numpy())
        L.append('micro-TPC (stage 48, zones by the PRELIMINARY gap) against this stage\'s T_d averaged over the same pairs:')
        for z, (lo, hi) in zip(('gap < 3.6 mm', 'gap 3.6-4.4 mm', 'gap 4.4-5.5 mm', 'gap > 5.5 mm'),
                               ((0, 3.6), (3.6, 4.4), (4.4, 5.5), (5.5, 30))):
            m_ = (pr.gap >= lo) & (pr.gap < hi) & (pr.s.abs() > 0.10)
            t_ = mt[mt.zone == z]
            if len(t_):
                L.append(f'   {z:15s}: micro-TPC T_d = {t_.T_d.iloc[0]:6.1f} +- {t_.e_T_d.iloc[0]:4.1f} ns;  pulse shape '
                         f'{pr[m_].Tf.mean():6.1f} ns  (gap {pr[m_].gapf.mean():.2f} mm, {int(m_.sum())} pairs)')

    # ---- H. surfaces -------------------------------------------------------------------------
    d45 = cfg.out_dir('45_deform_local_z')
    if os.path.isfile(f'{d45}/dz_reg24.csv'):
        z = pd.read_csv(f'{d45}/dz_reg24.csv').set_index('region')
        q24 = q24.join(z[['dz_c_corr', 'err_c_corr', 'dz_b', 'err_b', 'dz_b1', 'err_b1']])
        for zc, ec, nm in (('dz_c_corr', 'err_c_corr', 'centroid'), ('dz_b', 'err_b', 'boundary'),
                           ('dz_b1', 'err_b1', 'boundary, single pad')):
            k = q24[zc].notna()
            b, e, chi = wls((q24.gap[k] - GAP0).to_numpy(), q24.loc[k, zc].to_numpy(), q24.loc[k, ec].to_numpy())
            chis = {h: float((((q24.loc[k, zc] - np.average(q24.loc[k, zc] - h * (q24.gap[k] - GAP0), weights=1 / q24.loc[k, ec] ** 2)
                                - h * (q24.gap[k] - GAP0)) / q24.loc[k, ec]) ** 2).sum()) for h in (-0.5, 0.0, 0.5)}
            L.append(f'track height z_eff ({nm}) against the gap, 24 regions: slope {b[1]:+.2f} +- {e[1]:.2f} mm per mm '
                     f'(chi2 {chi:.1f}); chi2 for a fixed slope -0.5 (board moves) / 0 (both, symmetric) / +0.5 (window moves): '
                     f'{chis[-0.5]:.1f} / {chis[0.0]:.1f} / {chis[0.5]:.1f}')
        q24['z_mesh'] = q24.dz_c_corr - (q24.gap - GAP0) / 2
        q24['z_window'] = q24.dz_c_corr + (q24.gap - GAP0) / 2
    q24.reset_index().to_csv(f'{out}/regions.csv', index=False)
    RG.to_csv(f'{out}/regions_all.csv', index=False)       # 24 and 48 regions: T_fit, T_s, S, lnG, gap

    # ---- ridge geometry (if there is a thick zone) --------------------------------------------
    thick = P[P.gap > P.gap.median() + 2.0]
    prof_t = prof_l = None
    if len(thick) >= 20:
        c0 = thick[['pad_cx', 'pad_cy']].mean().to_numpy()
        _, sv, vt = np.linalg.svd(thick[['pad_cx', 'pad_cy']].to_numpy() - c0)
        axv = vt[0]; nrm = np.array([-axv[1], axv[0]])
        P['along'] = (P[['pad_cx', 'pad_cy']].to_numpy() - c0) @ axv
        P['perp'] = (P[['pad_cx', 'pad_cy']].to_numpy() - c0) @ nrm
        allp = pd.read_csv(cfg.MAP_CSV_PATH)[['pad_cx', 'pad_cy']].dropna().to_numpy(float)
        ap_, aa_ = (allp - c0) @ nrm, (allp - c0) @ axv
        side = 1.0 if abs(ap_[np.abs(aa_) < 60].max()) < abs(ap_[np.abs(aa_) < 60].min()) else -1.0
        edge = float((side * ap_[np.abs(aa_) < 60]).max())           # nearest edge of the full pad array
        P['perp'] *= side
        qt = P[np.abs(P.along) < 70]
        prof_t = qt.groupby(pd.cut(qt.perp, np.arange(-72, edge + 12, 12)), observed=True).agg(
            n=('gap', 'size'), perp=('perp', 'mean'), gap=('gap', 'mean'), T_fit_pad=('T_fit', 'mean'),
            T_s_pad=('T_s', 'mean')).dropna()
        ql = P[np.abs(P.perp) < 12]
        prof_l = ql.groupby(pd.cut(ql.along, np.arange(-240, 270, 30)), observed=True).agg(
            n=('gap', 'size'), along=('along', 'mean'), gap=('gap', 'mean'), x=('pad_cx', 'mean'), y=('pad_cy', 'mean')).dropna()
        ipk = prof_t.gap.idxmax()
        half = prof_l[prof_l.gap > 0.5 * (prof_l.gap.max() + P.gap.median())]
        L.append(f'thick zone ({len(thick)} pads more than 2 mm above the median): a ridge centred at pad ({c0[0]:.0f}, {c0[1]:.0f}) mm, '
                 f'axis at {np.degrees(np.arctan2(axv[1], axv[0])) % 180:.0f} deg; peak gap {P.gap.max():.1f} mm; across the ridge the gap '
                 f'peaks {edge - prof_t.perp[ipk]:.0f} mm inside the edge of the full pad array (live pads reach {edge - qt.perp.max():.0f} mm '
                 f'from it) and is {prof_t.gap.iloc[0]:.1f} mm at {edge - prof_t.perp.iloc[0]:.0f} mm; along the ridge it is above '
                 f'half height over {half.along.max() - half.along.min() + 30:.0f} mm and back to {prof_l.gap.iloc[0]:.1f} / {prof_l.gap.iloc[-1]:.1f} mm at the two ends')
        prof_t.to_csv(f'{out}/ridge_profile_across.csv'); prof_l.to_csv(f'{out}/ridge_profile_along.csv')
    thin = P[P.gap < P.gap.median() - 0.5]
    if len(thin) >= 10:
        L.append(f'thin zone ({len(thin)} pads more than 0.5 mm below the median): centred at pad ({thin.pad_cx.mean():.0f}, '
                 f'{thin.pad_cy.mean():.0f}) mm, minimum {P.gap.min():.2f} mm')
    GC = None
    if os.path.isfile(f'{d45}/gap_classes.csv'):
        GC = pd.read_csv(f'{d45}/gap_classes.csv')
        gf = pd.read_csv(f'{d45}/gap_fits.csv')
        L.append('track height against the gap, pad level (stage 45; slope in mm of mid-gap height per mm of gap: -0.5 = only the '
                 'board moves, +0.5 = only the window moves, 0 = both by the same amount):')
        for r_ in gf[gf.fit == 'no tilt'].itertuples():
            L.append(f'      {r_.est:32s} b = {r_.b:+.3f} +- {r_.eb:.3f}   delta-chi2 for b = -0.5 / 0 / +0.5: '
                     f'{r_.dchi2_m05:.0f} / {r_.dchi2_0:.1f} / {r_.dchi2_p05:.0f}')
        lo_, hi_ = GC.iloc[(GC.gap_mean_c - GAP0).abs().argmin()], GC.iloc[GC.gap_mean_c.argmax()]
        dz1, e1 = (hi_.c_1p - lo_.c_1p) / 0.54, np.hypot(hi_.e_c1p, lo_.e_c1p) / 0.54     # single-pad centroid: response 0.54 (toy)
        dz2, e2 = hi_.b_1f - lo_.b_1f, np.hypot(hi_.e_1f, lo_.e_1f)
        w1, w2 = 1 / e1 ** 2, 1 / e2 ** 2
        dzm, edz = (w1 * dz1 + w2 * dz2) / (w1 + w2), 1 / np.sqrt(w1 + w2)
        dg = hi_.gap_mean_c - lo_.gap_mean_c
        L.append(f'   thickest class ({hi_.gap_mean_c:.2f} mm) against the nominal class ({lo_.gap_mean_c:.2f} mm): mid-gap height changes by '
                 f'{dzm:+.2f} +- {edz:.2f} mm for a gap change of {dg:+.2f} mm -> window {dg/2 + dzm:+.2f} +- {edz:.2f} mm, '
                 f'board (mesh) {-dg/2 + dzm:+.2f} +- {edz:.2f} mm (M3 z, towards the window positive)')

    # ---- G. effects per gap class ------------------------------------------------------------
    hits, rays, g, _ = sc.load(a.run_key)
    d34 = cfg.out_dir('34_clustering_matching')
    A = np.load(f'{d34}/alignment.npz')
    T = cm.make_T(g, R=A['R'], s=A['s'], pad_mean=A['pad_mean'], m3_mean=A['m3_mean'])
    ra = pd.read_pickle(f'{d34}/rays_aligned.pkl')
    pxy = ((ra[['x', 'y']].to_numpy() - T.m3_mean) / T.s) @ T.R + T.pad_mean
    dd, ii = cKDTree(xy).query(pxy)
    ra['gap'] = P.gap.to_numpy()[ii]; ra['dpad'] = dd; ra['lnG'] = P.lnG.to_numpy()[ii]
    ra['nohit'] = ~ra.eventId.isin(hits.eventId)
    lat = cfg.out_dir('41_pillar_lattice') + '/lattice_pad.npz'
    if os.path.isfile(lat):
        Lp = np.load(lat)
        ij = np.array([(i, j) for i in range(-300, 301) for j in range(-300, 301)])
        pts = Lp['origin'] + ij[:, :1] * Lp['a1'] + ij[:, 1:] * Lp['a2']
        pts = pts[(pts[:, 0] > xy[:, 0].min() - 20) & (pts[:, 0] < xy[:, 0].max() + 20)
                  & (pts[:, 1] > xy[:, 1].min() - 20) & (pts[:, 1] < xy[:, 1].max() + 20)]
        ra['dl'] = cKDTree(pts).query(pxy)[0]
    qr = ra[(ra.sd >= 2) & (ra.dpad < 9)].copy()
    qr['gc'] = pd.cut(qr.gap, CLASS_EDGES, labels=False)
    tk['gap'] = P.set_index('channel_id').gap.reindex(tk.channel_id.to_numpy()).to_numpy()
    tk['gc'] = pd.cut(tk.gap, CLASS_EDGES, labels=False)
    lead = R[R.is_hit & R.good & ~R.saturated & (R.amp > 300) & (R.amp < 3500) & (R.imax >= 4) & (R.imax <= 12)]
    lead = lead.sort_values('amp', ascending=False).drop_duplicates('eventId')
    o_ = R.loc[R.is_hit & R.good & (R.amp > 100), ['eventId', 'pad_cx', 'pad_cy']].merge(
        lead[['eventId', 'pad_cx', 'pad_cy']], on='eventId', suffixes=('', '_L'))
    o_ = o_[np.hypot(o_.pad_cx - o_.pad_cx_L, o_.pad_cy - o_.pad_cy_L) < 25.0]
    lead = lead.join(o_.groupby('eventId').size().rename('nh'), on='eventId')
    lead = lead[lead.nh.fillna(1) == 1].copy()
    lead['gap'] = P.set_index('channel_id').gap.reindex(lead.channel_id.to_numpy()).to_numpy()
    lead['gc'] = pd.cut(lead.gap, CLASS_EDGES, labels=False)
    lead['t_r'] = lead.t_l20 - lead.groupby('channel_id').t_l20.transform('median')
    lead['tp_r'] = lead.t_pk - lead.groupby('channel_id').t_pk.transform('median')
    crow = []
    pe = lambda p, n: 100 * np.sqrt(max(p * (1 - p), 1e-9) / max(n, 1))
    for c_, nm in enumerate(CLASS_NAMES):
        pp = P[P.gclass == c_]; t_ = tk[tk.gc == c_]; hh = t_[t_.is_hit]; rr = qr[qr.gc == c_]; ll = lead[lead.gc == c_]
        o = dict(gap_class=nm, n_pads=len(pp), gap=pp.gap.mean(), T_d=pp.T_d.mean(), E_loc=pp.E_loc.mean(),
                 sigma_t_l20=pdf.s68(ll.t_r), sigma_t_pk=pdf.s68(ll.tp_r), n_single=len(ll),
                 amp_med=hh.ampc.median(), amp_mpv=dm.mpv(hh.ampc.to_numpy()), landau_lnwidth=pdf.s68(np.log(hh.ampc)),
                 frac_below_300=(hh.ampc < 300).mean(), lnG=pp.lnG.median(), S=pp.S.median(), n_tracks=len(rr))
        for (lo, hi) in ((0, 8), (8, 14), (14, 40)):
            s_ = rr[(rr.theta >= lo) & (rr.theta < hi)]
            o[f'nohit_{lo}_{hi}'] = 100 * s_.nohit.mean(); o[f'e_nohit_{lo}_{hi}'] = pe(s_.nohit.mean(), len(s_))
        o['nohit_all'] = 100 * rr.nohit.mean(); o['e_nohit_all'] = pe(rr.nohit.mean(), len(rr))
        if 'dl' in rr:
            for (lo, hi) in ((0, 8), (8, 14)):
                s_ = rr[(rr.theta >= lo) & (rr.theta < hi)]
                nr, fr = s_[s_.dl < 0.5], s_[s_.dl > 1.25]
                o[f'near_{lo}_{hi}'] = 100 * nr.nohit.mean(); o[f'e_near_{lo}_{hi}'] = pe(nr.nohit.mean(), len(nr))
                o[f'far_{lo}_{hi}'] = 100 * fr.nohit.mean(); o[f'e_far_{lo}_{hi}'] = pe(fr.nohit.mean(), len(fr))
        crow.append(o)
    CL = pd.DataFrame(crow)
    CL.to_csv(f'{out}/classes.csv', index=False)
    L.append('effects per drift-gap class [mm]:')
    L.append('   class    pads  gap   T_d   E_loc  sigma_t(l20/pk)  amp(med) ln-width  <300ADC   lnG    no-hit % (0-8 / 8-14 / >14 deg)      shadow: near / far, vertical')
    for r_ in CL.itertuples():
        sh = (f'{r_.near_0_8:5.1f}+-{r_.e_near_0_8:3.1f} / {r_.far_0_8:4.2f}+-{r_.e_far_0_8:4.2f}' if 'near_0_8' in CL else '')
        L.append(f'   {r_.gap_class:8s} {r_.n_pads:4d} {r_.gap:5.2f} {r_.T_d:5.0f} {r_.E_loc:6.0f}   {r_.sigma_t_l20:4.1f} / {r_.sigma_t_pk:4.1f}     '
                 f'{r_.amp_med:6.0f}   {r_.landau_lnwidth:.3f}   {100*r_.frac_below_300:4.1f} %  {r_.lnG:+.2f}   '
                 f'{r_.nohit_0_8:4.2f}+-{r_.e_nohit_0_8:.2f} / {r_.nohit_8_14:4.2f}+-{r_.e_nohit_8_14:.2f} / '
                 f'{r_.nohit_14_40:4.2f}+-{r_.e_nohit_14_40:.2f}   {sh}')
    # ---- lost muons against local gap, local gain and angle (track-level logistic fit) ----------
    ql = qr[qr.gap.notna() & qr.lnG.notna()].copy()
    ql['tan'] = np.tan(np.radians(ql.theta))

    def logit(cols, y):
        from scipy.optimize import minimize
        X = np.column_stack([np.ones(len(y))] + cols)
        nll = lambda b_: np.sum(np.logaddexp(0, X @ b_) - y * (X @ b_))
        jac = lambda b_: X.T @ (1 / (1 + np.exp(-(X @ b_))) - y)
        r_ = minimize(nll, np.zeros(X.shape[1]), jac=jac, method='BFGS')
        pr = 1 / (1 + np.exp(-(X @ r_.x)))
        return r_.x, np.sqrt(np.diag(np.linalg.inv(X.T @ (X * (pr * (1 - pr))[:, None])))), r_.fun
    yl = ql.nohit.to_numpy(float)
    b3, e3, n3 = logit([np.log(ql.gap / GAP0).to_numpy(), ql.lnG.to_numpy(), ql.tan.to_numpy()], yl)
    _, _, n1 = logit([ql.tan.to_numpy()], yl)
    L.append(f'lost muons, logistic fit on {len(ql):,} tracks (no-hit {100*yl.mean():.2f} %): log-odds = a '
             f'{b3[1]:+.2f} +- {e3[1]:.2f} * ln(gap / {GAP0:g} mm) {b3[2]:+.2f} +- {e3[2]:.2f} * ln(G / median) '
             f'{b3[3]:+.2f} +- {e3[3]:.2f} * tan(theta); -lnL {n3:.1f} against {n1:.1f} with the angle only. '
             f'One e-fold of intrinsic gain multiplies the odds of losing a muon by {np.exp(b3[2]):.2f}, '
             f'a gap 1.5 x thicker by {np.exp(b3[1]*np.log(1.5)):.2f}.')
    if 'dl' in ql:
        for lab, m_ in (('within 0.6 mm of a pillar', ql.dl < 0.6), ('more than 1.2 mm from a pillar', ql.dl > 1.2)):
            bb, ee, _ = logit([np.log(ql.gap[m_] / GAP0).to_numpy(), ql.lnG[m_].to_numpy(), ql.tan[m_].to_numpy()],
                              yl[m_.to_numpy()])
            L.append(f'   {lab} ({int(m_.sum()):,} tracks, no-hit {100*yl[m_.to_numpy()].mean():.2f} %): ln(gap) {bb[1]:+.2f} +- {ee[1]:.2f}, '
                     f'ln G {bb[2]:+.2f} +- {ee[2]:.2f}, tan(theta) {bb[3]:+.2f} +- {ee[3]:.2f}')
    nom = ql[(ql.gap > 3.6) & (ql.gap < 4.4)]
    if len(nom) > 3000:
        nom = nom.assign(gt=pd.qcut(nom.lnG, 3, labels=False))
        L.append('   nominal-gap class by intrinsic-gain tercile (ln G mean: no-hit all angles / below 8 deg): '
                 + '; '.join(f'{gq.lnG.mean():+.2f}: {100*gq.nohit.mean():.2f} / {100*gq.nohit[gq.theta < 8].mean():.2f} %'
                             for _, gq in nom.groupby('gt')))
    L.append('assumptions: reference drift time = the stage-47 argument (see its summary; the micro-TPC of stage 48 gives the absolute scale), '
             f'gap = T_d inverted with {gas_name}; lattice of stage 41 (measured, pad frame) for the shadow; '
             'active tracks = at least 2 mm inside the live area (stage-34 trimmed alignment); no-hit = no P2 hit in the event; '
             'gain correction assumes primary charge proportional to the gap; amplification-gap conversion assumes a 128 um gap '
             'and dry-gas Townsend coefficients (no Penning transfer), i.e. only the logarithmic derivative is used.')

    P.to_csv(f'{out}/pads.csv', index=False)

    # ---- figures -----------------------------------------------------------------------------
    fig, ax = plt.subplots(2, 4, figsize=(26, 11.5))
    smap(ax[0, 0], P, 'gap', 'drift gap [mm]', 'viridis', label='mm')
    smap(ax[0, 1], P, 'T_d', 'drift time across the gap [ns]', 'viridis', label='ns')
    smap(ax[0, 2], P, 'sl_s_l20', 'time spread of single-pad pulses, s68 [ns]', 'viridis', label='ns')
    smap(ax[0, 3], P, 'amp_sl', 'track-matched amplitude x cos(theta), median [ADC]', 'magma', label='ADC')
    smap(ax[1, 0], P, 'lnG', 'intrinsic gain, ln(G / median)  (gap effects removed)', sym=0.0)
    smap(ax[1, 1], P, 'dgap_amp_um', 'amplification-gap deviation if G is all gap [um]', sym=0.0, label='um')
    smap(ax[1, 2], P, 'S', 'ion-tail term S = mean-pulse width not due to drift [ns]', sym=float(np.nanmedian(P.S)), label='ns')
    smap(ax[1, 3], P, 'sl_dt0', 'fitted extra delay dt0 [ns]', sym=0.0, label='ns')
    fig.suptitle(f'{cfg.DET_NAME} {a.run_key}: deformation maps (per pad, sliding 18 mm)', fontsize=14)
    fig.tight_layout(); fig.savefig(f'{out}/deformation_maps.png', dpi=110); plt.close(fig)

    fig, ax = plt.subplots(1, 4, figsize=(25, 5.8))
    q = q24.reset_index()
    if 'gap_scan' in q:
        ax[0].errorbar(q.gap, q.gap_scan, q.e_gap_scan, fmt='o', color='#257de0', label='drift-scan timing (24 Nov)')
        ax[0].errorbar(q.gap, q.gap_ampcurve, q.e_gap_ampcurve, fmt='s', color='#ec4647', mfc='none',
                       label='amplitude-curve stretch (24 Nov)')
    lim = [2.8, max(7.5, float(q.gap.max()) + 0.3)]
    ax[0].plot(lim, lim, color='#9c9c9c', lw=1)
    ax[0].set(xlabel='drift gap from the pulse shape, this run [mm]', ylabel='drift gap from the drift scan [mm]',
              title='Gap per region: two days, two methods', xlim=lim, ylim=lim)
    ax[0].legend(frameon=False, fontsize=9)
    ax[1].errorbar(q.T_s, q.T_fit, xerr=q.e_T_s, yerr=q.e_T_fit, fmt='o', color='#1baf7a')
    tl = [95, 215]; ax[1].plot(tl, tl, color='#9c9c9c', lw=1)
    ax[1].set(xlabel='T_d from the pulse-to-pulse spread [ns]', ylabel='T_d from the mean-pulse fit [ns]',
              title='Two pulse-shape estimators, 24 regions')
    q48 = RG[RG.K == 48]
    ax[2].errorbar(q48.lnG, q48.S, q48.e_S, fmt='o', color='#7a4fd1')
    pf = np.polyfit(q48.lnG, q48.S, 1); xx = np.array([q48.lnG.min(), q48.lnG.max()])
    ax[2].plot(xx, np.polyval(pf, xx), color='#7a4fd1', lw=1)
    ax[2].set(xlabel='intrinsic gain, ln(G / median)', ylabel='S = T_fit - T_spread [ns]',
              title=f'Ion tail vs gain, 48 regions: {pf[0]:+.0f} ns per e-fold')
    if GC is not None:
        z0 = float(np.average(GC.b_1f, weights=1 / GC.e_1f ** 2)); z1 = float(np.average(GC.c_1p, weights=1 / GC.e_c1p ** 2))
        ax[3].errorbar(GC.gap_mean_b, GC.b_1f - z0, GC.e_1f, fmt='s', color='#f2a600', label='pad-boundary step, single-pad events')
        ax[3].errorbar(GC.gap_mean_c + 0.04, (GC.c_1p - z1) / 0.54, GC.e_c1p / 0.54, fmt='o', color='#257de0',
                       label='centroid of single-pad clusters (/ 0.54 response)')
        gx = np.array([3.0, float(GC.gap_mean_c.max()) + 0.3])
        for h, ls, nm in ((-0.5, '--', 'only the board moves'), (0.5, ':', 'only the window moves'), (0.0, '-', 'both, equally')):
            ax[3].plot(gx, h * (gx - GAP0), ls, color='#5c5c5c', lw=1.2, label=nm)
        ax[3].set(xlabel='drift gap class [mm]', ylabel='mid-gap height from the M3 tracks, relative [mm]',
                  title='Which surface is displaced? (stage 45, by gap class)')
        ax[3].legend(frameon=False, fontsize=8)
    for a_ in ax:
        a_.grid(color='#e6e6e6')
    fig.tight_layout(); fig.savefig(f'{out}/consistency.png', dpi=120); plt.close(fig)

    fig, ax = plt.subplots(1, 4, figsize=(25, 5.6))
    xg = CL.gap.to_numpy()
    for (lo, hi), col in zip(((0, 8), (8, 14), (14, 40)), ('#ec4647', '#f2a600', '#257de0')):
        ax[0].errorbar(xg, CL[f'nohit_{lo}_{hi}'], CL[f'e_nohit_{lo}_{hi}'], fmt='o-', color=col, label=f'{lo}-{hi} deg')
    ax[0].set(xlabel='drift gap [mm]', ylabel='tracks without any hit [%]', title='Lost muons vs drift gap', ylim=(0, None))
    ax[0].legend(frameon=False)
    if 'near_0_8' in CL:
        ax[1].errorbar(xg, CL.near_0_8, CL.e_near_0_8, fmt='o-', color='#ec4647', label='< 0.5 mm from a pillar, 0-8 deg')
        ax[1].errorbar(xg, CL.near_8_14, CL.e_near_8_14, fmt='s-', color='#f2a600', label='< 0.5 mm from a pillar, 8-14 deg')
        ax[1].errorbar(xg, CL.far_0_8, CL.e_far_0_8, fmt='o--', color='#9c9c9c', label='> 1.25 mm away, 0-8 deg')
        ax[1].set(xlabel='drift gap [mm]', ylabel='tracks without any hit [%]', title='Pillar shadow depth vs drift gap', ylim=(0, None))
        ax[1].legend(frameon=False, fontsize=9)
    ax[2].plot(xg, CL.amp_med / CL.amp_med.iloc[1], 'o-', color='#7a4fd1', label='measured amplitude')
    ax[2].plot(xg, xg / xg[1], '--', color='#9c9c9c', label='primary charge (gap)')
    ax[2].plot(xg, np.exp(CL.lnG - CL.lnG.iloc[1]), 's-', color='#1baf7a', mfc='none', label='intrinsic gain (after correction)')
    ax[2].set(xlabel='drift gap [mm]', ylabel='relative to the 3.6-4.4 mm class', title='Signal size vs drift gap')
    ax[2].legend(frameon=False, fontsize=9)
    ax[3].plot(CL.T_d, CL.sigma_t_l20, 'o-', color='#257de0', label='20 % leading edge')
    ax[3].plot(CL.T_d, CL.sigma_t_pk, 's-', color='#ec4647', mfc='none', label='peak')
    ax[3].plot([0, CL.T_d.max()], [0, CL.sigma_t_l20.iloc[1] / CL.T_d.iloc[1] * CL.T_d.max()], color='#9c9c9c', lw=1)
    ax[3].set(xlabel='drift time across the gap [ns]', ylabel='time spread of single-pad pulses, s68 [ns]',
              title='Time resolution vs drift time', xlim=(0, None), ylim=(0, None))
    ax[3].legend(frameon=False)
    for a_ in ax:
        a_.grid(color='#e6e6e6')
    fig.tight_layout(); fig.savefig(f'{out}/effects.png', dpi=120); plt.close(fig)

    if prof_t is not None:
        fig, ax = plt.subplots(1, 2, figsize=(14, 5.2))
        ax[0].plot(edge - prof_t.perp, prof_t.gap, 'o-', color='#257de0')
        ax[0].axvline(edge - qt.perp.max(), color='#9c9c9c', ls=':'); ax[0].text(edge - qt.perp.max() + 1, P.gap.median(), 'last live pads', fontsize=8, rotation=90)
        ax[0].set(xlabel='distance from the edge of the pad array [mm] (across the ridge, central 140 mm)', ylabel='drift gap [mm]',
                  title='Across the thick zone', xlim=(0, None))
        ax[1].plot(prof_l.along, prof_l.gap, 'o-', color='#ec4647')
        ax[1].set(xlabel='position along the ridge [mm]', ylabel='drift gap [mm]', title='Along the thick zone (within 12 mm of its axis)')
        for a_ in ax:
            a_.axhline(P.gap.median(), color='#9c9c9c', lw=1); a_.grid(color='#e6e6e6')
        fig.tight_layout(); fig.savefig(f'{out}/ridge_profile.png', dpi=120); plt.close(fig)

    with open(f'{out}/summary.txt', 'w') as f:
        f.write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
