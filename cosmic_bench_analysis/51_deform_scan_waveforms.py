#!/usr/bin/env python3
"""51_deform_scan_waveforms.py -- drift time vs drift field from the WAVEFORMS of a drift scan.

RESULT ON det2 (27 Nov 2026): the method closes on simulated pulses, but on the data the absolute value is
driven by the single lowest-field point and trades against a +17 ns delay of that pulse (153 ns without
delays, flat between 80 and 150 ns with them), so it does NOT fix the absolute scale with 30 min per point.
What the scan waveforms do give: drift times relative to the reference point, the time-spread law
s68 = k T at every field (k constant to 4 % from 104 to 334 ns), and a pulse shape that hardly changes with
the mesh voltage (ion transit short).


On one set of pads the single-electron response (ion tail (x) shaper) and the delays are the same at every
drift voltage; only the arrival box changes:  M_k(t) = box(T_k) (x) SER.  Hence for two scan points

        box(T_j) (x) M_k  =  box(T_k) (x) M_j        (same SER, same delay)

and a fit of this relation over all points gives every T_k in ABSOLUTE terms: no pulse-response model, no
zero-time reference, no diffusion/threshold correction (the limits of stages 44, 48, 50). The lever arm is
large (T from ~110 ns at high field to several hundred ns at low field), so the reference drift time itself
is well determined. T(dV) = d / v_d(dV / d) then gives the gap d and the gas (water fraction f) separately.

Steps
  1. extract (once, cached per sub-run): hit and track-near waveforms of every sub-run that has decoded_root,
     same rules as stage 43, tracks with the series alignment (as stages 42/50).
  2. per zone (pads of a drift-gap class of stage 49 of the reference run) and scan point: mean pulse on a
     10 ns grid (sum of waveforms / sum of amplitudes, leading pad of each event), peak time, single-pad
     spread, amplitude.
  3. drift scan: profile of the reference drift time T_ref (the nominal-voltage point): for each T_ref, every
     other point k gets its best (T_k, scale, optionally a delay); the summed squared residual vs T_ref has
     its minimum at the true value. Bootstrap over pulses for the errors.
  4. gas fit: T_k = d / v_f(dV_k / d) on the Magboltz water family (gas_studies/water_scan), f and d free.
  5. mesh scan (constant drift field): the ion transit time scales as 1 / V_mesh, so the same relation with
     boxes I_ref * V_ref / V gives the ion transit time I_ref in absolute terms.

Out (<scan run>/51_deform_scan_waveforms/): <sub_run>/{wf_table.pkl, wf.npy}, points.csv, drift_times.csv,
gas_fit.csv, spread_vs_drift_time.csv, mean_pulses.npz, scan_waveforms.png, summary.txt
Usage: python3 51_deform_scan_waveforms.py [scan_key] [ref_key] [--rebuild]
"""
import glob
import importlib
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'gas_studies', 'water_scan'))
import p2_qa_config as qa       # noqa: E402
import p2_align as pa           # noqa: E402
import p2_sparks as ps          # noqa: E402
import p2_study_cache as sc     # noqa: E402
import p2_waveforms as pw       # noqa: E402
import p2_deform as pdf         # noqa: E402
s42 = importlib.import_module('42_pillar_shadow_vs_field')

TG = np.arange(100.0, 1900.0, 10.0)
FIT = (TG > 140) & (TG < 1850)
RNG = np.random.default_rng(20261128)
ADC_MAX, NEAR, AMP_MIN = 4095.0, 22.0, 150.0
ZONES = {'nominal': (3.6, 4.4), 'thin': (0.0, 3.6), 'ridge': (5.5, 30.0)}


# ------------------------------------------------------------------------------------------ extraction
def extract(cfg, ct, d, T, z, out_sub):
    ctm = ct[ct['mapped']].drop_duplicates(['feu', 'channel'])
    hits_dir = os.path.join(d, 'combined_hits_root'); m3_dir = os.path.join(d, 'm3_tracking_root')
    ep = pa.load_m3_endpoints(m3_dir, qa.M3_CHI2_CUT, qa.M3_MIN_NCLUS).drop_duplicates('eventId')
    bad = set()
    hv = os.path.join(d, 'hv_monitor.csv')
    if os.path.isfile(hv):
        bad = set(int(e) for e in ps.SparkVeto.from_csv(hv, cfg).vetoed_ids_from_hits(
            hits_dir, ct.attrs['feus'], min_amp=cfg.MIN_AMP))
    x, y = pa.project_to_z(ep, z)
    pxy = s42.inverse(T, np.column_stack([x, y]))
    dz = ep['Z_Up'].to_numpy() - ep['Z_Down'].to_numpy()
    tm = np.column_stack([(ep['X_Up'].to_numpy() - ep['X_Down'].to_numpy()) / dz,
                          (ep['Y_Up'].to_numpy() - ep['Y_Down'].to_numpy()) / dz])
    tp = (tm / T.s) @ T.R
    trk = pd.DataFrame({'xp': pxy[:, 0], 'yp': pxy[:, 1], 'tx_p': tp[:, 0], 'ty_p': tp[:, 1]},
                       index=ep.eventId.to_numpy(np.int64))
    dec = pw.list_decoded_files(os.path.join(d, 'decoded_root'), feus=ct.attrs['feus'])
    tabs, waves = [], []
    for _, fr in dec.iterrows():
        feu = int(fr.feu)
        c = ctm[ctm.feu == feu]
        chans = c.channel.to_numpy(int)
        tree = cKDTree(c[['pad_cx', 'pad_cy']].to_numpy(float))
        is_mapped = np.zeros(pw.N_CH, bool); is_mapped[chans] = True
        ped = pw.pedestal_from_data(fr.path)
        rms = pw.noise_rms_from_data(fr.path, ped)
        for evs, ts, ft, wf in pw.iter_decoded_events(fr.path, step=2000):
            w = pw.subtract_common_noise(wf - ped[None, :, None])
            peak = w.max(axis=2)
            nab = (w > (3.0 * rms)[None, :, None]).sum(axis=2)
            hit = (peak > 5.0 * rms[None, :]) & (nab >= 3) & is_mapped[None, :]
            keep = hit.copy(); near = np.zeros_like(hit)
            has_t = np.isin(evs, trk.index.to_numpy())
            if has_t.any():
                k_idx = np.flatnonzero(has_t)
                txy = trk.loc[evs[k_idx], ['xp', 'yp']].to_numpy()
                for k, cand in zip(k_idx, tree.query_ball_point(txy, r=NEAR)):
                    if len(cand):
                        near[k, chans[cand]] = True
                keep |= near
            iev, ich = np.nonzero(keep)
            if not len(iev):
                continue
            tabs.append(pd.DataFrame({
                'eventId': evs[iev].astype(np.int64), 'feu': feu, 'channel': ich, 'rms': rms[ich],
                'ftst': ft[iev].astype(np.int16), 'is_hit': hit[iev, ich], 'near_track': near[iev, ich],
                'saturated': wf[iev, ich, :].max(axis=1) >= ADC_MAX - 4, 'peak': peak[iev, ich]}))
            waves.append(w[iev, ich, :].astype(np.float32))
    if not tabs:
        return None
    R = pd.concat(tabs, ignore_index=True)
    R = R.merge(ctm[['feu', 'channel', 'channel_id', 'connector_N', 'pad_cx', 'pad_cy']], on=['feu', 'channel'], how='left')
    R = R.join(trk, on='eventId')
    R['d_trk'] = np.hypot(R.pad_cx - R.xp, R.pad_cy - R.yp)
    R['good'] = ~R.eventId.isin(bad)
    os.makedirs(out_sub, exist_ok=True)
    R.to_pickle(f'{out_sub}/wf_table.pkl'); np.save(f'{out_sub}/wf.npy', np.vstack(waves))
    return len(R)


# ------------------------------------------------------------------------------------------ pulse tools
def mean_pulse(WF, rows, ftst, amp):
    w = WF[rows].astype(float)
    tb = ((np.arange(32)[None, :] * 60 + ftst[:, None] * 10) // 10).astype(int)
    nb = 32 * 6 + 6
    num = np.bincount(tb.ravel(), weights=w.ravel(), minlength=nb)
    den = np.bincount(tb.ravel(), weights=np.repeat(amp, 32), minlength=nb)
    t = np.arange(nb) * 10.0
    m = num / np.where(den > 0, den, np.nan)
    ok = np.isfinite(m)
    return np.interp(TG, t[ok], m[ok])


def boxconv(m, T, dt=10.0):
    c = np.concatenate([[0.0], np.cumsum(m) * dt])
    te = np.concatenate([[TG[0] - dt / 2], TG + dt / 2])
    return (np.interp(TG, te, c) - np.interp(TG - T, te, c, left=0.0)) / T


def pair_fit(m_k, m_ref, t_ref, free_dt=False, t0=None):
    """box(t_ref) (x) m_k(t - dt) * A = box(T_k) (x) m_ref  ->  (T_k, A, dt, sse)."""
    bref = None

    def resid(p):
        T, A, dt = p[0], p[1], (p[2] if free_dt else 0.0)
        mk = np.interp(TG - dt, TG, m_k, left=0.0, right=0.0) if free_dt else m_k
        return (boxconv(mk, t_ref) * A - boxconv(m_ref, T))[FIT]
    best = None
    for T0 in ([t0] if t0 else [t_ref, 2 * t_ref, 4 * t_ref]):
        p0 = [T0, 1.0] + ([0.0] if free_dt else [])
        lo = [15, 0.2] + ([-60] if free_dt else []); hi = [1500, 5] + ([60] if free_dt else [])
        r = least_squares(resid, p0, bounds=(lo, hi))
        if best is None or r.cost < best.cost:
            best = r
    return best.x[0], best.x[1], (best.x[2] if free_dt else 0.0), 2 * best.cost


def profile(M, ref, grid, free_dt=False):
    """Summed squared residual of all points against the reference, vs the trial reference drift time."""
    out = []
    for tr in grid:
        sse, Ts = 0.0, {}
        for k, m in M.items():
            if k == ref:
                continue
            T, A, dt, s_ = pair_fit(m, M[ref], tr, free_dt)
            sse += s_; Ts[k] = (T, A, dt)
        out.append((tr, sse, Ts))
    return out


def parab_min(x, y):
    i = int(np.argmin(y))
    if 0 < i < len(x) - 1:
        a, b, c = np.polyfit(x[i - 1:i + 2], y[i - 1:i + 2], 2)
        return -b / (2 * a)
    return x[i]


# ------------------------------------------------------------------------------------------ main
def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    scan_key = args[0] if args else 'det2_initial2'
    ref_key = args[1] if len(args) > 1 else 'det2_long_repos'
    cfg = qa.get_config(scan_key); cfg_r = qa.get_config(ref_key)
    out = os.path.abspath(os.path.join(cfg.OUT_BASE, '..', '51_deform_scan_waveforms')); os.makedirs(out, exist_ok=True)
    ct = sc.channel_table(cfg)
    T, z, _ = s42.align_series(cfg, ct)
    L = [f'{scan_key}: waveforms of the scan sub-runs; series alignment z = {z:.1f} mm; leading pad per event, '
         f'amplitude > {AMP_MIN:g} ADC, not saturated, peak sample 3-27; mean pulse = sum of waveforms / sum of amplitudes']
    subs = [('alignment', 430, 730, cfg.sub_run_dir)]
    for d in sorted(glob.glob(os.path.join(cfg.run_dir, '*_scan_*'))):
        nm = os.path.basename(d); a_, b_ = nm.split('_')[-2:]
        subs.append(('drift' if nm.startswith('drift') else 'mesh', int(a_), int(b_), d))
    gp = pd.read_csv(cfg_r.out_dir('49_deform_decomposition') + '/pads.csv').set_index('channel_id').gap
    data = {}
    for kind, mesh, drift, d in subs:
        if not glob.glob(os.path.join(d, 'decoded_root', '*.root')):
            continue
        osub = os.path.join(out, os.path.basename(d))
        if '--rebuild' in sys.argv or not os.path.isfile(f'{osub}/wf.npy'):
            n = extract(cfg, ct, d, T, z, osub)
            print(f'extracted {os.path.basename(d)}: {n} rows', flush=True)
        if not os.path.isfile(f'{osub}/wf.npy'):
            continue
        R = pd.read_pickle(f'{osub}/wf_table.pkl'); WF = np.load(f'{osub}/wf.npy')
        F = pdf.features(WF, R.ftst.to_numpy())
        for c in F.columns:
            R[c] = F[c].to_numpy()
        R['row'] = np.arange(len(R)); R['gap'] = gp.reindex(R.channel_id.to_numpy()).to_numpy()
        data[(kind, mesh, drift)] = (R, WF)
    if not data:
        sys.exit('no decoded_root in any sub-run of ' + cfg.run_dir)

    # ---- per point and zone: leading pulses, mean pulse, simple observables -------------------
    rows, MP = [], {}
    for key, (R, WF) in data.items():
        kind, mesh, drift = key
        h = R[R.is_hit & R.good & ~R.saturated]
        h = h.sort_values('amp', ascending=False).drop_duplicates('eventId')
        o = R.loc[R.is_hit & R.good & (R.amp > 100), ['eventId', 'pad_cx', 'pad_cy']].merge(
            h[['eventId', 'pad_cx', 'pad_cy']], on='eventId', suffixes=('', '_L'))
        o = o[np.hypot(o.pad_cx - o.pad_cx_L, o.pad_cy - o.pad_cy_L) < 25.0]
        h = h.join(o.groupby('eventId').size().rename('nh'), on='eventId')
        s = h[(h.amp > AMP_MIN) & (h.amp < 3500) & (h.imax >= 3) & (h.imax <= 27)]
        for zn, (lo, hi) in ZONES.items():
            q = s[(s.gap >= lo) & (s.gap < hi)]
            if len(q) < 150:
                continue
            MP[(zn,) + key] = (q.row.to_numpy(), q.ftst.to_numpy(), q.amp.to_numpy(), WF)
            q1 = q[q.nh.fillna(1) == 1]
            tr = q1.t_l20 - q1.groupby('channel_id').t_l20.transform('median')
            rows.append(dict(zone=zn, kind=kind, mesh=mesh, drift=drift, dV=drift - mesh, n=len(q), n_single=len(q1),
                             amp_med=q.amp.median(), t_pk=q.t_pk.median(), t_l20=q.t_l20.median(),
                             s68_l20=pdf.s68(tr), rise=q.rise.median(), fwhm=q.fwhm.median()))
    PT = pd.DataFrame(rows); PT.to_csv(f'{out}/points.csv', index=False)

    def pulses(zn, kind, pool_nominal=True, boot=False):
        """Mean pulse per dV (drift) or per mesh V; the alignment run is pooled with the same-voltage point."""
        acc = {}
        for (z_, k_, mesh, drift), (rw, ft, am, WF) in MP.items():
            if z_ != zn:
                continue
            if kind == 'drift' and not (k_ == 'drift' or (k_ == 'alignment' and pool_nominal)):
                continue
            if kind == 'mesh' and not (k_ == 'mesh' or k_ == 'alignment'):
                continue
            x = (drift - mesh) if kind == 'drift' else mesh
            if boot:
                i = RNG.integers(0, len(rw), len(rw)); rw, ft, am = rw[i], ft[i], am[i]
            acc.setdefault(x, []).append((WF, rw, ft, am))
        M = {}
        for x, lst in acc.items():
            num = sum(mean_pulse(WF, rw, ft, am) * am.sum() for WF, rw, ft, am in lst)
            M[x] = num / sum(am.sum() for _, _, _, am in lst)
        return M

    # ---- drift scan: absolute drift times ------------------------------------------------------
    import analyze_water as aw
    gas = aw.load_gases()
    tabs = sorted([(f, gas[k]) for k, f in aw.FAMILIES['h2o'] if k in gas], key=lambda a: a[0])
    grid = np.arange(60.0, 201.0, 5.0)
    res_rows, fit_rows, keep = [], [], {}
    for zn in ZONES:
        M = pulses(zn, 'drift')
        M = {k: v for k, v in M.items() if k >= 45}                 # dV = 25 V: longitudinal diffusion blurs the box
        if len(M) < 6 or 300 not in M:
            continue
        ref = 300
        pr = profile(M, ref, grid)
        t_ref = parab_min(np.array([p[0] for p in pr]), np.array([p[1] for p in pr]))
        prd = profile(M, ref, grid, free_dt=True)
        t_ref_dt = parab_min(np.array([p[0] for p in prd]), np.array([p[1] for p in prd]))
        best = {k: pair_fit(m, M[ref], t_ref) for k, m in M.items() if k != ref}
        bestd = {k: pair_fit(m, M[ref], t_ref, True) for k, m in M.items() if k != ref}
        # bootstrap over pulses
        bt, bT = [], {k: [] for k in best}
        for _ in range(12):
            Mb = {k: v for k, v in pulses(zn, 'drift', boot=True).items() if k in M}
            pb = profile(Mb, ref, grid[::2])
            tb = parab_min(np.array([p[0] for p in pb]), np.array([p[1] for p in pb]))
            bt.append(tb)
            for k in best:
                bT[k].append(pair_fit(Mb[k], Mb[ref], tb, t0=best[k][0])[0])
        e_ref = float(np.std(bt))
        L.append(f'zone {zn}: {len(M)} drift points; reference drift time (dV = {ref} V) from the profile over all points: '
                 f'T_ref = {t_ref:.1f} +- {e_ref:.1f} ns (bootstrap over pulses); with a free delay per point: {t_ref_dt:.1f} ns, '
                 f'fitted delays {np.median([v[2] for v in bestd.values()]):+.1f} ns median, {np.min([v[2] for v in bestd.values()]):+.1f} .. '
                 f'{np.max([v[2] for v in bestd.values()]):+.1f}')
        sse0 = {p_[0]: p_[1] for p_ in pr}; sse1 = {p_[0]: p_[1] for p_ in prd}
        lowk = min(M)
        L.append(f'   caution: the profile is driven by the lowest-field point (dV = {lowk} V). Summed squared residual x 1e3 at '
                 f'T_ref = 110 / 120 / 150 ns: no delay {1e3*sse0[110.0]:.1f} / {1e3*sse0[120.0]:.1f} / {1e3*sse0[150.0]:.1f}; free delays '
                 f'{1e3*sse1[110.0]:.2f} / {1e3*sse1[120.0]:.2f} / {1e3*sse1[150.0]:.2f} (flat). A delay of '
                 f'{pair_fit(M[lowk], M[ref], 115.0, True)[2]:+.0f} ns of the dV = {lowk} V pulse makes T_ref = 115 ns fit as well. '
                 'So this method does NOT fix the absolute drift time with 30 min per point; the no-delay value disagrees '
                 'with the pad pairs (stage 48: raw plateau alone gives T >= 105 ns and <= 136 ns for any diffusion).')
        res_rows.append(dict(zone=zn, dV=ref, T=t_ref, e_T=e_ref, scale=1.0, dt_free=0.0))
        for k in sorted(best):
            res_rows.append(dict(zone=zn, dV=k, T=best[k][0], e_T=float(np.std(bT[k])), scale=best[k][1], dt_free=bestd[k][2]))
        keep[zn] = (M, pr, t_ref)
        # gas fit: T_k = d / v_f(dV / d)
        q = pd.DataFrame([r for r in res_rows if r['zone'] == zn]).sort_values('dV')
        dV, y = q.dV.to_numpy(float), q['T'].to_numpy()
        ye = np.hypot(q.e_T.to_numpy(), 0.02 * y)            # 2 % per point for what the bootstrap does not see
        G = []
        for f in np.arange(0.0, 2.01, 0.05):
            gm = aw.GasModel.mix(tabs, f)
            for d in np.arange(3.0, 6.51, 0.05):
                Tm = (d / 10.0) / gm.v(dV / (d / 10.0))
                G.append(dict(zone=zn, f=f, d=d, chi2=float((((y - Tm) / ye) ** 2).sum())))
        G = pd.DataFrame(G); fit_rows.append(G)
        b = G.loc[G.chi2.idxmin()]; ok = G[G.chi2 < G.chi2.min() + 2.3]
        d0 = G[np.isclose(G.f, 0.0)]
        L.append(f'   gas fit (water family, no nuisance parameter): water {b.f:.2f} % ({ok.f.min():.2f}-{ok.f.max():.2f}), '
                 f'gap {b.d:.2f} mm ({ok.d.min():.2f}-{ok.d.max():.2f}), chi2 {b.chi2:.1f}/{len(q)-2}; dry gas, best gap: chi2 {d0.chi2.min():.0f}')
        L.append('   dV [V]: ' + ' '.join(f'{v:5.0f}' for v in dV))
        L.append('   T [ns]: ' + ' '.join(f'{v:5.0f}' for v in y))
        L.append('   +-    : ' + ' '.join(f'{v:5.1f}' for v in q.e_T))
    DT = pd.DataFrame(res_rows); DT.to_csv(f'{out}/drift_times.csv', index=False)
    if fit_rows:
        pd.concat(fit_rows).to_csv(f'{out}/gas_fit.csv', index=False)

    # ---- mesh scan: ion transit time ----------------------------------------------------------
    Mm = pulses('nominal', 'mesh')
    if len(Mm) >= 4 and 430 in Mm:
        ig = np.arange(20.0, 301.0, 10.0); sse = []
        for I in ig:
            tot = 0.0
            for V, m in Mm.items():
                if V == 430:
                    continue
                # box(I) (x) m_V * A = box(I * 430 / V) (x) m_430, scale A free
                lhs = boxconv(m, I)[FIT]; rhs = boxconv(Mm[430], I * 430.0 / V)[FIT]
                A = float((lhs * rhs).sum() / (lhs * lhs).sum())
                tot += float(((lhs * A - rhs) ** 2).sum())
            sse.append(tot)
        sse = np.array(sse); I_best = parab_min(ig, sse)
        L.append(f'mesh scan (nominal zone, {len(Mm)} voltages {min(Mm)}-{max(Mm)} V at constant drift field): ion transit time at 430 V from '
                 f'box(I) scaling as 1/V: I = {I_best:.0f} ns (summed squared residual {sse.min():.3e} at the minimum, '
                 f'{sse[0]:.3e} at 20 ns, {sse[-1]:.3e} at 300 ns)')
        mrow = PT[(PT.zone == 'nominal') & PT.kind.isin(['mesh', 'alignment'])].sort_values('mesh')
        L.append('   mesh V: ' + ' '.join(f'{v:6.0f}' for v in mrow.mesh))
        L.append('   t_pk  : ' + ' '.join(f'{v:6.1f}' for v in mrow.t_pk))
        L.append('   fwhm  : ' + ' '.join(f'{v:6.1f}' for v in mrow.fwhm))
        L.append('   amp   : ' + ' '.join(f'{v:6.0f}' for v in mrow.amp_med))
        L.append('   n     : ' + ' '.join(f'{v:6d}' for v in mrow.n))
    # ---- time-spread law at other fields: s68 / T with T relative to the pad-pair value at the reference -----
    d48 = cfg_r.out_dir('48_deform_microtpc_dt0.85') + '/drift_time_by_zone.csv'
    tref_ext = {'nominal': 118.0, 'thin': 102.0, 'ridge': 166.0}
    if os.path.isfile(d48):
        t48 = pd.read_csv(d48); t48 = t48[t48.est == 't_pk'].set_index('zone').T_d
        tref_ext = {'nominal': float(t48.get('gap 3.6-4.4 mm', 118.0)), 'thin': float(t48.get('gap < 3.6 mm', 102.0)),
                    'ridge': float(t48.get('gap > 5.5 mm', 166.0))}
    krows = []
    for zn, (M, _, _) in keep.items():
        q = PT[(PT.zone == zn) & PT.kind.isin(['drift', 'alignment'])]
        q = q.groupby('dV').apply(lambda g_: pd.Series({'s68': np.average(g_.s68_l20, weights=g_.n_single), 'n1': g_.n_single.sum()})).reset_index()
        for r_ in q.itertuples():
            k_ = int(r_.dV)
            if k_ not in M:
                continue
            Tk = tref_ext[zn] if k_ == 300 else pair_fit(M[k_], M[300], tref_ext[zn], True)[0]
            krows.append(dict(zone=zn, dV=k_, T=Tk, s68=r_.s68, k=r_.s68 / Tk, n_single=int(r_.n1)))
    KT = pd.DataFrame(krows); KT.to_csv(f'{out}/spread_vs_drift_time.csv', index=False)
    for zn, q in KT.groupby('zone', sort=False):
        hi_ = q[q.dV >= 100]
        L.append(f'time-spread law, zone {zn} (drift times relative to {tref_ext[zn]:.0f} ns at 300 V, pad pairs with wet-gas diffusion): '
                 f'single-pad s68 / T = {hi_.k.mean():.3f} +- {hi_.k.std():.3f} (rms over {len(hi_)} points, dV 100-525 V, '
                 f'T {hi_["T"].min():.0f}-{hi_["T"].max():.0f} ns); at dV 50 / 75 V (T {q["T"].max():.0f} ns): '
                 + ' / '.join(f'{v:.3f}' for v in q[q.dV < 100].k))
    # ---- absolute drift time from the SCALING of the time spread: if s68 = k T with k independent of the field,
    #      the low-field points (long lever arm) fix the reference drift time: s68_k / s68_ref = T_k(T_ref) / T_ref
    for zn, (M, _, _) in keep.items():
        q = PT[(PT.zone == zn) & PT.kind.isin(['drift', 'alignment'])]
        q = q.groupby('dV').apply(lambda g_: pd.Series({'s68': np.average(g_.s68_l20, weights=g_.n_single), 'n1': g_.n_single.sum()}))
        if 300 not in q.index:
            continue
        low = [k_ for k_ in sorted(M) if k_ <= 125 and k_ in q.index]
        tg_ = np.arange(90.0, 171.0, 5.0); c2 = []
        for tr_ in tg_:
            c_ = 0.0
            for k_ in low:
                Tk = pair_fit(M[k_], M[300], tr_, True)[0]
                pred = q.s68[300] * Tk / tr_; err = q.s68[k_] * np.hypot(1 / np.sqrt(q.n1[k_]), 1 / np.sqrt(q.n1[300]))
                c_ += ((q.s68[k_] - pred) / err) ** 2
            c2.append(c_)
        c2 = np.array(c2); tb_ = parab_min(tg_, c2); ok_ = tg_[c2 <= c2.min() + 1.0]
        L.append(f'absolute drift time from the time-spread scaling, zone {zn}: the spread at dV = {", ".join(str(k_) for k_ in low)} V '
                 f'relative to 300 V equals the drift-time ratio for T_ref = {tb_:.0f} ns ({ok_.min():.0f}-{ok_.max():.0f}; chi2 {c2.min():.1f} '
                 f'for {len(low)} points; chi2 at 150 ns: {c2[np.argmin(abs(tg_ - 150))]:.0f}). Assumes s68 / T_d independent of the '
                 'drift field (found the same within 10 % for the four gap classes at 475-890 V/cm); the no-delay profile value '
                 'would need s68 / T_d to be 30 % larger at 125 V/cm.')
        pd.DataFrame(dict(zone=zn, T_ref=tg_, chi2=c2)).to_csv(f'{out}/drift_time_from_spread_{zn}.csv', index=False)
    if len(Mm) >= 4:
        kk = (TG >= 150) & (TG <= 1000)
        cen = {V_: float((TG[kk] * m_[kk]).sum() / m_[kk].sum()) for V_, m_ in Mm.items()}
        vs_ = sorted(cen); sl = np.polyfit(np.log(np.array(vs_) / 430.0), [cen[v_] for v_ in vs_], 1)[0]
        rms_ = float(np.std([cen[v_] for v_ in vs_]))
        L.append('   mesh scan, centroid of the mean pulse [ns] (150-1000 ns): ' + ' '.join(f'{v_}: {cen[v_]:.1f}' for v_ in vs_))
        L.append(f'   slope d centroid / d ln V = {sl:+.1f} ns (rms of the points {rms_:.1f} ns). An ion box I scaling as 1/V moves the '
                 f'centroid by -I/2 per unit ln V, so I = {max(-2*sl, 0):.0f} ns, below ~{2*(abs(sl) + 2*rms_/0.235):.0f} ns at 2 sigma: '
                 'the pulse shape hardly changes with the mesh voltage.')
        # same at FIXED pulse height (250-700 ADC), so that an amplitude-dependent response cannot hide a change
        def _shape(m_):
            c_ = float((TG[kk] * m_[kk]).sum() / m_[kk].sum()); i_ = int(np.argmax(m_)); h_ = m_[i_] / 2
            j_ = i_
            while m_[j_] > h_:
                j_ -= 1
            up_ = TG[j_] + (h_ - m_[j_]) / (m_[j_ + 1] - m_[j_]) * 10
            j_ = i_
            while m_[j_] > h_:
                j_ += 1
            return c_, TG[j_ - 1] + (m_[j_ - 1] - h_) / (m_[j_ - 1] - m_[j_]) * 10 - up_
        fa = {}
        for (z_, k_, mesh_, drift_), (rw, ft, am, WF_) in MP.items():
            if z_ != 'nominal' or not (k_ in ('mesh', 'alignment') or (k_ == 'drift' and drift_ - mesh_ == 300)):
                continue
            sel_ = (am >= 250) & (am < 700)
            if sel_.sum() >= 150:
                fa.setdefault(mesh_, []).append((int(sel_.sum()),) + _shape(mean_pulse(WF_, rw[sel_], ft[sel_], am[sel_])))
        if len(fa) >= 4:
            vs_ = sorted(fa); n_ = np.array([sum(x[0] for x in fa[v_]) for v_ in vs_])
            c_ = np.array([np.average([x[1] for x in fa[v_]], weights=[x[0] for x in fa[v_]]) for v_ in vs_])
            f_ = np.array([np.average([x[2] for x in fa[v_]], weights=[x[0] for x in fa[v_]]) for v_ in vs_])
            lx = np.log(np.array(vs_) / 430.0); sc_ = np.polyfit(lx, c_, 1, w=np.sqrt(n_))[0]; sf_ = np.polyfit(lx, f_, 1, w=np.sqrt(n_))[0]
            pd.DataFrame(dict(mesh=vs_, n=n_, centroid=c_, fwhm=f_)).to_csv(f'{out}/mesh_scan_shape_fixed_amp.csv', index=False)
            L.append('   at fixed pulse height (250-700 ADC): mesh V / centroid / FWHM [ns]: '
                     + '  '.join(f'{v_}: {a_:.1f} / {b_:.1f}' for v_, a_, b_ in zip(vs_, c_, f_)))
            L.append(f'   slopes: d centroid / d ln V = {sc_:+.1f} ns, d FWHM / d ln V = {sf_:+.1f} ns. One e-fold of gain is '
                     f'{1/0.032:.0f} V = {1/0.032/430:.3f} in ln V, so changing the gain by a factor e WITH THE VOLTAGE changes the '
                     f'pulse width by {sf_/0.032/430:+.1f} ns: the -16 +- 5 ns per e-fold seen across the board (stage 49) is '
                     'not a generic "more gain, shorter pulse" effect. A gap change is the 7x stronger lever on the ion transit '
                     '(T_ion ~ gap^2 / V: +52 % per e-fold through the gap, -7.5 % through the voltage). An ion box of 240-350 ns '
                     '(mobility 1.5-2.2 cm2/Vs, 150 um, 430 V) scaling as 1/V would move the centroid by 27-40 ns over this scan; '
                     'it moves by less than 3 ns. Centroid: I < ~40 ns; FWHM narrowing: I ~ 100 ns; board slope: I ~ 60 ns.')
    np.savez(f'{out}/mean_pulses.npz', t=TG, **{f'{zn}_dV{k}': m for zn, (M, _, _) in keep.items() for k, m in M.items()},
             **{f'mesh_{k}': m for k, m in Mm.items()})

    # ---- figure ----------------------------------------------------------------------------------
    fig, ax = plt.subplots(1, 4, figsize=(27, 5.8))
    if 'nominal' in keep:
        M, pr, t_ref = keep['nominal']
        cm_ = plt.get_cmap('viridis')
        ks = sorted(M)
        for i, k in enumerate(ks):
            ax[0].plot(TG, M[k], color=cm_(i / max(len(ks) - 1, 1)), lw=1.3, label=f'{k / 0.4:.0f} V/cm' if i % 3 == 0 else None)
        ax[0].set(xlim=(100, 1500), xlabel='time after the trigger [ns]', ylabel='mean pulse (sum wf / sum amp)',
                  title='Mean pulse of the nominal-gap pads at every drift field'); ax[0].legend(frameon=False, fontsize=8)
        ax[1].plot([p[0] for p in pr], [p[1] for p in pr], 'o-', color='#257de0')
        ax[1].axvline(t_ref, color='#ec4647'); ax[1].set(xlabel='trial drift time of the 300 V point [ns]', ylabel='summed squared residual, all points',
                                                          title=f'Reference drift time: {t_ref:.0f} ns', yscale='log')
    for zn, col in (('nominal', '#257de0'), ('thin', '#84b9f3'), ('ridge', '#0e3d7e')):
        q = DT[DT.zone == zn].sort_values('dV')
        if len(q):
            ax[2].errorbar(q.dV / 0.4, q['T'], q.e_T, fmt='o', color=col, label=f'{zn} pads')
    if fit_rows:
        G = fit_rows[0]; b = G.loc[G.chi2.idxmin()]
        xx = np.linspace(45, 530, 200)
        for (f, d, ls, lab) in ((b.f, b.d, '-', f'best: {b.f:.2f} % water, {b.d:.2f} mm'), (0.0, 4.0, ':', 'dry gas, 4.0 mm')):
            gm = aw.GasModel.mix(tabs, f)
            ax[2].plot(xx / 0.4, (d / 10.0) / gm.v(xx / (d / 10.0)), ls, color='#5c5c5c', label=lab)
        Z = G.pivot(index='d', columns='f', values='chi2'); zz = Z.to_numpy() - Z.to_numpy().min()
        ax[3].contourf(Z.columns, Z.index, zz, levels=[0, 2.3, 6.2, 11.8], colors=['#257de0', '#84b9f3', '#cce3fd'])
        ax[3].axhline(4.0, color='#9c9c9c', lw=1)
        ax[3].set(xlabel='water fraction [%]', ylabel='gap of the nominal zone [mm]', title='68 / 95 / 99.7 % regions, no nuisance parameter')
    ax[2].set(yscale='log', xlabel='nominal drift field dV / 4 mm [V/cm]', ylabel='drift time across the gap [ns]',
              title='Absolute drift time from the pulse shapes'); ax[2].legend(frameon=False, fontsize=8)
    for a_ in ax:
        a_.grid(color='#e6e6e6')
    fig.tight_layout(); fig.savefig(f'{out}/scan_waveforms.png', dpi=110); plt.close(fig)
    with open(f'{out}/summary.txt', 'w') as fh:
        fh.write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
