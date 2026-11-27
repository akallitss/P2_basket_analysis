#!/usr/bin/env python3
"""52_deform_tradeoff.py -- time resolution against lost muons, vs track angle, for different drift gaps.

det2 is a natural experiment: its drift gap goes from 3.3 to 7 mm across the board (stage 49), and the
24 Nov scans vary the drift voltage (0-525 V) and the mesh voltage (340-420 V). Together they separate
what the GAP does from what the drift FIELD and the GAIN do, which the long run alone cannot (at one
voltage a thicker gap is also a lower field).

  MISSES   every reference track is a Bernoulli trial (no hit in the chamber). Two-component model fitted by
           maximum likelihood on all tracks of the long run and of every scan point, with the local gap and
           local intrinsic gain of the pad under the track (stage 49 map), the local drift field
           E = dV / gap, the mesh voltage and the angle:
             P_miss = 1 - (1 - P_thr)(1 - P_sh)
             P_thr  = sigma(a0 + a_q lnQ + a_x x + a_x2 x^2 + day)        signal below threshold
             P_sh   = sigma(c0 + c_q lnQ + c_x x + c_x2 x^2 + c_gap ln(gap/4) + day) * exp(-gap tan(theta) / L0)
                                                                         pillar shadow, fading with the
                                                                         sideways reach of the track
             lnQ = lnG + kV (V_mesh - 430) + ln(gap/4)   (gain x primary charge),   x = ln(E / 750 V/cm)
           ('day' = offsets of the scan day against the long-run day: same voltages, different losses.)
           A plain logistic model in the same variables fits worse by 95 units of log-likelihood.
  TIMING   pulse-to-pulse spread s68 of the 20 % leading-edge time of the leading pad (long-run waveforms),
           per gap class and angle, for single-pad events and for all events. It is proportional to the drift
           time T = gap / v_d(E) (stages 47/49), so  sigma(gap, E, angle) = k(angle) * gap / v_d(E), with
           k(angle) measured and v_d(E) the velocity measured in the chamber (stage 49). The spread includes
           the trigger (scintillator) jitter.
  TRADE-OFF  for a chamber built with a given gap and run at a chosen drift FIELD (the drift voltage scaled
           with the gap) and at the median gain of det2: lost muons and time spread per angle.

All gaps are on the dry-gas scale of stage 49 unless --gap-scale is given (best estimate of the true gap =
scale x dry-scale gap, section 8 of the audit); the scale only relabels the gap axis.

Out (<run>/52_deform_tradeoff/): tracks.pkl (cache), miss_model.csv, timing_by_class_angle.csv,
tradeoff.csv, validation.csv, tradeoff.png, summary.txt
Usage: python3 52_deform_tradeoff.py [ref_key] [scan_key] [--gap-scale 1.11] [--rebuild]
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
from scipy.optimize import minimize
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import p2_qa_config as qa       # noqa: E402
import p2_study_cache as sc     # noqa: E402
import p2_study_common as cm    # noqa: E402
import p2_mapping as pmap       # noqa: E402
import p2_deform as pdf         # noqa: E402
RNG = np.random.default_rng(20261128)
s42 = importlib.import_module('42_pillar_shadow_vs_field')

ANG = [(0, 5), (5, 10), (10, 15), (15, 20), (20, 30)]
GAPS = [3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5]
CLASS_EDGES = [0, 3.6, 4.4, 5.5, 30]
E0 = 750.0
NAMES = ['a0', 'a_q', 'a_x', 'a_x2', 'c0', 'c_q', 'c_x', 'c_x2', 'c_gap', 'L0', 'kV', 'day_thr', 'day_sh']
P0 = np.array([-5.5, -1.9, 0.1, 0.05, -2.7, -0.6, -1.2, 0.45, -1.7, 0.8, 0.032, -1.0, -0.3])
_sig = lambda z_: 1 / (1 + np.exp(-z_))


def miss_parts(p, lnG, dVm, gap, tan, E, day):
    """(P_threshold, P_shadow) of the two-component miss model."""
    a0, aq, ax, ax2, c0, cq, cx, cx2, cg, L0, kV, d1, d2 = p
    x = np.log(np.clip(E, 20, None) / E0)
    lnQ = lnG + kV * dVm + np.log(gap / 4.0)
    pthr = _sig(a0 + aq * lnQ + ax * x + ax2 * x * x + d1 * day)
    psh = _sig(c0 + cq * lnQ + cx * x + cx2 * x * x + cg * np.log(gap / 4.0) + d2 * day) * np.exp(-gap * tan / abs(L0))
    return pthr, psh


def pmiss(p, *a):
    pt, ps_ = miss_parts(p, *a)
    return 1 - (1 - pt) * (1 - ps_)


def fit_miss(arrs, y, start=P0):
    def nll(p):
        q = np.clip(pmiss(p, *arrs), 1e-9, 1 - 1e-9)
        return -(y * np.log(q) + (1 - y) * np.log(1 - q)).sum()
    r = minimize(nll, start, method='L-BFGS-B', options=dict(maxiter=3000, maxfun=60000))
    return r.x, r.fun


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    ref_key = args[0] if args else 'det2_long_repos'
    scan_key = args[1] if len(args) > 1 else 'det2_initial2'
    scale = float(sys.argv[sys.argv.index('--gap-scale') + 1]) if '--gap-scale' in sys.argv else 1.0
    cfg = qa.get_config(ref_key); cfs = qa.get_config(scan_key)
    out = cfg.out_dir('52_deform_tradeoff')
    P = pd.read_csv(cfg.out_dir('49_deform_decomposition') + '/pads.csv')
    pxy_pad = P[['pad_cx', 'pad_cy']].to_numpy(); ptree = cKDTree(pxy_pad)
    L = [f'{ref_key} + {scan_key}: gap and gain per track from the stage-49 map of {ref_key} (nearest pad, pad frame); '
         f'gap scale applied to the gap axis: x {scale:g}']

    # ---- tracks ------------------------------------------------------------------------------
    cache = f'{out}/tracks.pkl'
    if os.path.isfile(cache) and '--rebuild' not in sys.argv:
        TR = pd.read_pickle(cache)
    else:
        parts = []
        hits, rays, g, _ = sc.load(ref_key)
        d34 = cfg.out_dir('34_clustering_matching'); A = np.load(f'{d34}/alignment.npz')
        T = cm.make_T(g, R=A['R'], s=A['s'], pad_mean=A['pad_mean'], m3_mean=A['m3_mean'])
        ra = pd.read_pickle(f'{d34}/rays_aligned.pkl')
        xp = ((ra[['x', 'y']].to_numpy() - T.m3_mean) / T.s) @ T.R + T.pad_mean
        k = (ra.sd >= 2).to_numpy()
        parts.append(pd.DataFrame({'xp': xp[k, 0], 'yp': xp[k, 1], 'th': ra.theta.to_numpy()[k],
                                   'nohit': (~ra.eventId.isin(hits.eventId)).to_numpy()[k], 'mesh': 430, 'drift': 730,
                                   'kind': 'long', 'day': 0}))
        ct = sc.channel_table(cfs); ct.attrs['dead_guard'] = 2.0; ct.attrs['edge_guard'] = 2.0
        Ts, zs, _ = s42.align_series(cfs, ct)
        la = dict(np.load(os.path.join(cfg.out_dir('41_pillar_lattice'), 'lattice_pad.npz')))
        padc = ct[ct['mapped']].drop_duplicates('channel_id')
        lo = padc[['pad_cx', 'pad_cy']].to_numpy().min(0); hi = padc[['pad_cx', 'pad_cy']].to_numpy().max(0)
        Lfull = s42.lattice_points(la, lo, hi, margin=15)
        pil = pmap.load_pillars(cfs.MASK_GBR_PATH)
        bigp = pil[pil.big][['x', 'y']].to_numpy() if len(pil) else np.empty((0, 2))
        subs = [('alignment', 430, 730, cfs.sub_run_dir)]
        for d in sorted(glob.glob(os.path.join(cfs.run_dir, '*_scan_*'))):
            nm = os.path.basename(d); a_, b_ = nm.split('_')[-2:]
            subs.append(('drift' if nm.startswith('drift') else 'mesh', int(a_), int(b_), d))
        for kind, mesh, drift, d in subs:
            if not glob.glob(os.path.join(d, 'combined_hits_root', '*.root')):
                continue
            s_ = s42.point(cfs, ct, d, Ts, zs, Lfull, bigp)
            parts.append(pd.DataFrame({'xp': s_['xp'][:, 0], 'yp': s_['xp'][:, 1], 'th': s_['th'], 'nohit': s_['nohit'],
                                       'mesh': mesh, 'drift': drift, 'kind': kind, 'day': 1}))
            print(f'{kind} {mesh}/{drift}: {len(s_["th"])} tracks', flush=True)
        TR = pd.concat(parts, ignore_index=True)
        TR.to_pickle(cache)
    dd, ii = ptree.query(TR[['xp', 'yp']].to_numpy())
    TR['gap'] = P.gap.to_numpy()[ii]; TR['lnG'] = P.lnG.to_numpy()[ii]
    TR = TR[(dd < 9) & TR.gap.notna() & TR.lnG.notna() & (TR.th < 30)].copy()
    TR['tan'] = np.tan(np.radians(TR.th)); TR['dV'] = TR.drift - TR.mesh
    TR['E'] = TR.dV / (TR.gap / 10.0)
    TR['gc'] = pd.cut(TR.gap, CLASS_EDGES, labels=False)
    L.append(f'tracks: long run {int((TR.kind == "long").sum()):,}, alignment {int((TR.kind == "alignment").sum()):,}, '
             f'drift scan {int((TR.kind == "drift").sum()):,}, mesh scan {int((TR.kind == "mesh").sum()):,}')

    # ---- miss model --------------------------------------------------------------------------
    fitm = TR[(TR.dV >= 100)].copy()                 # below 250 V/cm nominal the chamber is in another regime
    arrs = (fitm.lnG.to_numpy(), (fitm.mesh - 430).to_numpy(float), fitm.gap.to_numpy(), fitm.tan.to_numpy(),
            fitm.E.to_numpy(), fitm.day.to_numpy(float))
    y = fitm.nohit.to_numpy(float)
    b, nll = fit_miss(arrs, y)
    BS = []
    for _ in range(25):                              # bootstrap over tracks: parameter and prediction errors
        i = RNG.integers(0, len(y), len(y))
        BS.append(fit_miss(tuple(a_[i] for a_ in arrs), y[i], start=b)[0])
    BS = np.array(BS); se = BS.std(0)
    pd.DataFrame({'term': NAMES, 'coef': b, 'err': se}).to_csv(f'{out}/miss_model.csv', index=False)
    L.append(f'miss model: {len(fitm):,} tracks with dV >= 100 V ({100*y.mean():.2f} % without hit), -lnL {nll:.1f}; '
             'errors from 25 bootstraps over tracks')
    for t, bb, ee in zip(NAMES, b, se):
        L.append(f'   {t:8s} {bb:+7.3f} +- {ee:.3f}')
    L.append(f'   kV = d ln G / d V_mesh implied by the misses: {b[10]:.4f} +- {se[10]:.4f} per V (the amplitudes give 0.033 per V); '
             f'L0 = sideways reach over which the pillar shadow falls by e: {abs(b[9]):.2f} +- {se[9]:.2f} mm')
    # validation: data vs model in cells
    val = []
    fitm['pm'] = pmiss(b, *arrs)
    for kind, grp in (('long', fitm[fitm.kind == 'long']), ('scan dV 250-350', fitm[(fitm.day == 1) & (fitm.mesh == 430) & fitm.dV.between(250, 350)]),
                      ('scan dV 100-225', fitm[(fitm.day == 1) & (fitm.mesh == 430) & fitm.dV.between(100, 225)]),
                      ('scan dV 375-525', fitm[(fitm.day == 1) & (fitm.mesh == 430) & fitm.dV.between(375, 525)]),
                      ('mesh 340-390', fitm[(fitm.kind == 'mesh') & (fitm.mesh <= 390)])):
        for c in range(4):
            for lo, hi in ANG:
                q = grp[(grp.gc == c) & (grp.th >= lo) & (grp.th < hi)]
                if len(q) < 150:
                    continue
                p_ = q.nohit.mean(); pm_ = q.pm.mean()
                val.append(dict(sample=kind, gap_class=c, th_lo=lo, th_hi=hi, n=len(q), data=100 * p_,
                                err=100 * np.sqrt(pm_ * (1 - pm_) / len(q)), model=100 * pm_))   # binomial error at the model value
    V = pd.DataFrame(val); V.to_csv(f'{out}/validation.csv', index=False)
    chi = (((V.data - V.model) / V.err) ** 2).sum()
    L.append(f'   validation in {len(V)} cells (sample x gap class x angle, >= 150 tracks): chi2 {chi:.0f} for {len(V)} cells')
    for smp, q in V.groupby('sample', sort=False):
        L.append(f'      {smp:16s} chi2 {(((q.data - q.model) / q.err) ** 2).sum():6.1f} / {len(q)} cells')

    # ---- timing from the long-run waveforms ---------------------------------------------------
    R, WF, _ = pdf.load(ref_key)
    gmap = P.set_index('channel_id')
    h = R[R.is_hit & R.good & ~R.saturated & (R.amp > 300) & (R.amp < 3500) & (R.imax >= 4) & (R.imax <= 12) & R.near_track]
    h = h.sort_values('amp', ascending=False).drop_duplicates('eventId')
    h = h[(h.d_trk < 9) & (h.sd >= 2)]
    o = R.loc[R.is_hit & R.good & (R.amp > 100), ['eventId', 'pad_cx', 'pad_cy']].merge(
        h[['eventId', 'pad_cx', 'pad_cy']], on='eventId', suffixes=('', '_L'))
    o = o[np.hypot(o.pad_cx - o.pad_cx_L, o.pad_cy - o.pad_cy_L) < 25.0]
    h = h.join(o.groupby('eventId').size().rename('nh'), on='eventId')
    h['gap'] = gmap.gap.reindex(h.channel_id.to_numpy()).to_numpy(); h['T_d'] = gmap.T_d.reindex(h.channel_id.to_numpy()).to_numpy()
    h['gc'] = pd.cut(h.gap, CLASS_EDGES, labels=False)
    h['tr'] = h.t_l20 - h.groupby('channel_id').t_l20.transform('median')
    trow = []
    for c in range(4):
        for lo, hi in ANG + [(0, 30)]:
            q = h[(h.gc == c) & (h.theta >= lo) & (h.theta < hi)]
            q1 = q[q.nh.fillna(1) == 1]
            if len(q1) < 80:
                continue
            trow.append(dict(gap_class=c, th_lo=lo, th_hi=hi, gap=q.gap.mean(), T_d=q.T_d.mean(), n_all=len(q), n_single=len(q1),
                             s_all=pdf.s68(q.tr), e_all=pdf.boot(q.tr, pdf.s68, 100), s_single=pdf.s68(q1.tr),
                             e_single=pdf.boot(q1.tr, pdf.s68, 100), frac_multi=1 - len(q1) / len(q)))
    TM = pd.DataFrame(trow); TM.to_csv(f'{out}/timing_by_class_angle.csv', index=False)
    TM['k_all'] = TM.s_all / TM.T_d; TM['k_single'] = TM.s_single / TM.T_d
    L.append('time spread s68 [ns] of the leading pad (20 % edge, pad medians removed), tracked muons, by gap class and angle; k = s68 / T_d:')
    L.append('   class  angle     gap    T_d   single-pad      all events     multi-pad   k_single  k_all')
    for r in TM.itertuples():
        L.append(f'   {r.gap_class}     {r.th_lo:2d}-{r.th_hi:<2d}  {r.gap:5.2f}  {r.T_d:5.0f}  {r.s_single:5.1f} +- {r.e_single:3.1f}   {r.s_all:5.1f} +- {r.e_all:3.1f}   '
                 f'{100*r.frac_multi:5.1f} %    {r.k_single:.3f}    {r.k_all:.3f}')
    # k(angle): weighted mean over the gap classes (proportionality to T_d), per angle bin
    kang = {}
    for lo, hi in ANG:
        q = TM[(TM.th_lo == lo) & (TM.th_hi == hi)]
        w1 = (q.s_single / q.e_single) ** 2; w2 = (q.s_all / q.e_all) ** 2
        kang[(lo, hi)] = (float((w1 * q.k_single).sum() / w1.sum()), float((w2 * q.k_all).sum() / w2.sum()),
                          float(np.sqrt(((q.k_single - (w1 * q.k_single).sum() / w1.sum()) ** 2 * w1).sum() / w1.sum())))
    L.append('   k(angle) = s68 / T_d averaged over the gap classes (single-pad / all events; rms between classes): '
             + '; '.join(f'{lo}-{hi} deg: {v[0]:.3f} / {v[1]:.3f} ({v[2]:.3f})' for (lo, hi), v in kang.items()))

    # ---- is the time spread proportional to the drift time at other fields? (stage 51, scan waveforms) -------
    f51 = os.path.abspath(os.path.join(cfs.OUT_BASE, '..', '51_deform_scan_waveforms', 'spread_vs_drift_time.csv'))
    k51 = 'not checked (run 51_deform_scan_waveforms.py)'
    if os.path.isfile(f51):
        K5 = pd.read_csv(f51); L.append('time-spread law at other drift fields (stage 51, single-pad events, all angles, scan day):')
        for zn, q in K5.groupby('zone', sort=False):
            h_ = q[q.dV >= 100]
            L.append(f'   zone {zn:8s} {len(h_)} points dV 100-525 V, T_d {h_["T"].min():.0f}-{h_["T"].max():.0f} ns: s68 / T_d = '
                     f'{h_.k.mean():.3f}, rms {h_.k.std():.3f} ({100*h_.k.std()/h_.k.mean():.0f} %); dV 50 / 75 V (T_d up to '
                     f'{q["T"].max():.0f} ns): ' + ' / '.join(f'{v:.3f}' for v in q[q.dV < 100].k))
        qn = K5[(K5.zone == 'nominal')]
        k51 = (f'checked on the scan waveforms: in the nominal zone s68 / T_d = {qn.k.mean():.3f} with an rms of '
               f'{100*qn.k.std()/qn.k.mean():.0f} % over 20 drift fields (T_d {qn["T"].min():.0f}-{qn["T"].max():.0f} ns)')

    # ---- drift velocity measured in the chamber (stage 49) -------------------------------------
    vt = pd.read_csv(cfg.out_dir('49_deform_decomposition') + '/drift_velocity_measured.csv')
    vd = lambda E: np.interp(E, vt.E_Vcm, vt.vd_measured_mm_ns)

    # ---- trade-off table -----------------------------------------------------------------------
    rows = []
    for mode, fields in (('fixed field', [500.0, 750.0, 1000.0]), ('fixed voltage', [300.0])):
        for fv in fields:
            for gp in GAPS:
                E = fv if mode == 'fixed field' else fv / (gp / 10.0)
                if E < vt.E_Vcm.min() or E > vt.E_Vcm.max():
                    continue
                Td = gp / float(vd(E))
                for lo, hi in ANG:
                    tn = np.tan(np.radians(0.5 * (lo + hi)))
                    arg = (0.0, 0.0, gp, tn, E, 0.0)
                    pt, psd = miss_parts(b, *arg)
                    pm = 100 * float(pmiss(b, *arg)); pe = 100 * float(np.std([pmiss(bb, *arg) for bb in BS]))
                    rows.append(dict(mode=mode, value=fv, gap_dry=gp, gap=gp * scale, E=E, dV=E * gp / 10.0, th_lo=lo, th_hi=hi, T_d=Td,
                                     miss=pm, e_miss=pe, miss_thr=100 * float(pt), miss_shadow=100 * float(psd),
                                     miss_day2=100 * float(pmiss(b, 0.0, 0.0, gp, tn, E, 1.0)),
                                     s_single=kang[(lo, hi)][0] * Td, s_all=kang[(lo, hi)][1] * Td))
    TO = pd.DataFrame(rows); TO.to_csv(f'{out}/tradeoff.csv', index=False)
    for mode, fv in (('fixed field', 750.0), ('fixed voltage', 300.0)):
        q = TO[(TO['mode'] == mode) & (TO.value == fv)]
        L.append(f'trade-off, {mode} {fv:g} {"V/cm" if mode == "fixed field" else "V"}, median gain of det2, long-run day: '
                 'lost muons [%] / time spread of the leading pad, all events [ns]')
        L.append('   gap [mm]  dV [V]  ' + '  '.join(f'{lo:2d}-{hi:<2d} deg   ' for lo, hi in ANG))
        for gp, qq in q.groupby('gap'):
            L.append(f'   {gp:6.2f}   {qq.dV.iloc[0]:5.0f}   ' + '  '.join(f'{r.miss:4.1f} / {r.s_all:4.1f}' for r in qq.itertuples()))
    L.append('assumptions: misses = no hit in the event for a track at least 2 mm inside the live area (trimmed alignment); '
             'gap, gain and field of the pad under the track from the stage-49 map (dry-gas gap scale inside the model; the '
             'local field is dV / gap on that scale); the model is fitted for dV >= 100 V; the time spread includes the '
             'trigger jitter and is taken proportional to the drift time, with the drift velocity measured in the chamber; '
             'k(angle) from the 430/730 V run; proportionality to the drift time ' + k51 + '; the day offset is the difference between '
             '24 and 25 Nov at identical voltages (gas / charging state) and is a floor on how well misses transfer between runs.')

    # ---- figure --------------------------------------------------------------------------------
    fig, ax = plt.subplots(1, 4, figsize=(27, 6))
    cmap = plt.get_cmap('viridis'); cg = {gp: cmap(i / (len(GAPS) - 1)) for i, gp in enumerate(GAPS)}
    q = TO[(TO['mode'] == 'fixed field') & (TO.value == 750.0)]
    xa = [0.5 * (lo + hi) for lo, hi in ANG]
    for gp in GAPS:
        qq = q[q.gap_dry == gp]
        ax[0].errorbar(xa, qq.miss, qq.e_miss, fmt='o-', color=cg[gp], label=f'{gp * scale:.1f} mm')
        ax[1].plot(xa, qq.s_all, 'o-', color=cg[gp]); ax[1].plot(xa, qq.s_single, 's--', color=cg[gp], mfc='none', ms=4)
    ax[0].set(xlabel='track angle from the normal [deg]', ylabel='muons without a hit [%]', title='Lost muons at 750 V/cm, median gain (model)', ylim=(0, None))
    ax[0].legend(frameon=False, title='drift gap', ncol=2, fontsize=8)
    ax[1].set(xlabel='track angle from the normal [deg]', ylabel='time spread s68 [ns]', title='Time spread at 750 V/cm (solid: all events, dashed: single-pad)', ylim=(0, None))
    for (lo, hi), mk in zip(ANG, ('o', 's', '^', 'D', 'v')):
        qq = q[(q.th_lo == lo)]
        ax[2].plot(qq.s_all, qq.miss, '-', color='#9c9c9c', lw=1)
        ax[2].scatter(qq.s_all, qq.miss, c=[cg[g_] for g_ in qq.gap_dry], marker=mk, s=45, label=f'{lo}-{hi} deg', zorder=3)
    ax[2].set(xlabel='time spread s68, all events [ns]', ylabel='muons without a hit [%]', title='Trade-off at 750 V/cm: each line is one angle, colour = gap', ylim=(0, None))
    ax[2].legend(frameon=False, fontsize=8)
    v_ = V[V['sample'] == 'long']
    for c, col in zip(range(4), ('#84b9f3', '#338ced', '#1b64be', '#0e3d7e')):
        qq = v_[v_.gap_class == c]
        ax[3].errorbar(0.5 * (qq.th_lo + qq.th_hi), qq.data, qq.err, fmt='o', color=col)
        ax[3].plot(0.5 * (qq.th_lo + qq.th_hi), qq.model, '-', color=col, label=f'gap class {c}')
    ax[3].set(xlabel='track angle [deg]', ylabel='muons without a hit [%]', title='Check: long run, data (points) and model (lines) by gap class', ylim=(0, None))
    ax[3].legend(frameon=False, fontsize=8)
    for a_ in ax:
        a_.grid(color='#e6e6e6')
    fig.tight_layout(); fig.savefig(f'{out}/tradeoff.png', dpi=110); plt.close(fig)
    with open(f'{out}/summary.txt', 'w') as fh:
        fh.write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
