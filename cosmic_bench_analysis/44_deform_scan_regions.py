#!/usr/bin/env python3
"""44_deform_scan_regions.py -- hit-level, region-by-region view of the det2
(P2_2) drift-field and mesh-voltage scans of 24 Nov 2026, as input to the
mechanical-deformation audit (does the drift gap differ across the chamber?).

Per sub-run (alignment 430/730, drift scan 430/430..955, mesh scan 340..420 /
mesh+300) the LEADING pad of every event is kept (max amplitude among mapped,
kept pads; spark/burst-vetoed events and hot/dead/shorted pads removed as in
stages 16/42). Timing/shape quantities use leading pads that are neither
saturated nor truncated; amplitude statistics use every leading pad.

Hit-tree definitions (mm_dream_reconstruction/waveform_analysis, checked in
the source and on this data):
  time_of_max  parabola-interpolated peak position of the RAW waveform x 60 ns,
               + ftst x 10 ns (DREAM fine timestamp, applied in the analyzer;
               verified: left_sample of left-truncated pulses = ftst/6 exactly)
  time         30 % constant-fraction point of the leading edge (raw waveform,
               linear interpolation, relative to the local baseline) x 60 ns,
               + ftst x 10 ns
  time_over_threshold  [ns] interpolated crossings of the 5-sigma gate on the
               MATCHED-FILTER gate waveform (5-sample = 300 ns boxcar), x 60 ns
  left/right_sample   those crossings in SAMPLE units, + ftst/6
  integral     sum of (raw sample - local baseline) over the full pulse span
               (out to the return to baseline), ADC x samples
  amplitude    parabola-fitted peak - local baseline (median of 4 samples
               ending 2 samples before the pulse); saturated pulses get a
               linear edge extrapolation (can exceed 4095)

Regions: the 43_deformation partition (reg8/reg24, pad frame) of the 11-25
run, valid here (same chamber, same cabling).

Usage: python3 44_deform_scan_regions.py [--rebuild]
Out:   <Analysis>/det2/<11-24 run>/44_deform_scan_regions/
"""
import argparse
import glob
import os
import re
import sys
import warnings

import numpy as np
import pandas as pd
import uproot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from scipy.interpolate import CubicSpline
from scipy import stats as sst

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2_qa_config as qa        # noqa: E402
import p2_io as p2io             # noqa: E402
import p2_mapping as pmap        # noqa: E402
import p2_sparks as ps           # noqa: E402
import p2_study_cache as sc      # noqa: E402

warnings.filterwarnings('ignore', category=RuntimeWarning)
HERE = os.path.dirname(os.path.abspath(__file__))
REGIONS = ('/local/home/ak271430/Documents/PostDocSaclay/data/Cosmic_Bench/Analysis/det2/'
           'p2_det2_long_run_430_repositioned_11-25-26/long_run_det2_430_730/43_deformation/regions.csv')
GAS1 = os.path.join(HERE, 'garfield_inputs', 'gas_table.csv')
GAS2 = os.path.join(HERE, '..', 'gas_studies', 'results', 'ar_iso_95_5.csv')
GAP_CM = 0.40
CUTS = (300.0, 150.0)
BANDS = [(150, 300), (300, 600), (600, 1200), (1200, 2400), (2400, 4000)]
B = 300
RNG = np.random.default_rng(44)
BR = ['eventId', 'trigger_timestamp_ns', 'channel', 'feu', 'amplitude', 'time', 'time_of_max',
      'time_over_threshold', 'integral', 'saturated', 'trunc_left', 'trunc_right', 'significance']
GRID = '#e6e6e6'


# ----------------------------------------------------------------------------- gas
def vdrift_func():
    g1 = pd.read_csv(GAS1)[['E_Vcm', 'vz']]
    g2 = pd.read_csv(GAS2)
    g2 = g2[g2.region == 'drift'][['E_Vcm', 'vz_cm_ns']].rename(columns={'vz_cm_ns': 'vz'})
    g = pd.concat([pd.DataFrame({'E_Vcm': [0.0], 'vz': [0.0]}), g2[g2.E_Vcm < 100], g1]).sort_values('E_Vcm')
    E = g.E_Vcm.to_numpy(float); v = g.vz.to_numpy(float)

    def vd(e):  # cm/ns; linear (mobility) between 0 and 50 V/cm
        return np.interp(np.clip(e, 1.0, E[-1]), E, v)
    return vd, g


# ----------------------------------------------------------------------------- stats
def boot_q(x, qs=(10, 16, 50, 84, 90)):
    """percentiles + bootstrap errors of median and sigma68"""
    x = np.asarray(x, float); n = len(x)
    if n < 5:
        return dict(n=n, med=np.nan, med_e=np.nan, s68=np.nan, s68_e=np.nan, p10=np.nan, p90=np.nan)
    p = np.percentile(x, qs)
    idx = RNG.integers(0, n, (B, n))
    bp = np.percentile(x[idx], [16, 50, 84], axis=1)
    return dict(n=n, med=p[2], med_e=float(np.std(bp[1])), s68=0.5 * (p[3] - p[1]),
                s68_e=float(np.std(0.5 * (bp[2] - bp[0]))), p10=p[0], p90=p[4])


def mpv_logkde(a, nb=400, bw=0.06):
    """robust MPV: mode of a Gaussian-smoothed histogram of log10(amplitude)"""
    la = np.log10(a[a > 0])
    if len(la) < 20:
        return np.nan
    lo, hi = np.percentile(la, [0.5, 99.5])
    h, e = np.histogram(la, nb, (lo - 0.3, hi + 0.3))
    c = 0.5 * (e[1:] + e[:-1]); dx = c[1] - c[0]
    k = np.exp(-0.5 * (np.arange(-4 * bw, 4 * bw + dx, dx) / bw) ** 2)
    s = np.convolve(h, k, 'same')
    return 10 ** c[np.argmax(s)]


def amp_stats(a):
    a = np.asarray(a, float); n = len(a)
    if n < 20:
        return dict(n_all=n, amp_med=np.nan, amp_med_e=np.nan, amp_mpv=np.nan, amp_mpv_e=np.nan)
    idx = RNG.integers(0, n, (100, n))
    meds = np.median(a[idx], axis=1)
    mpvs = np.array([mpv_logkde(a[i]) for i in idx])
    return dict(n_all=n, amp_med=float(np.median(a)), amp_med_e=float(np.std(meds)),
                amp_mpv=mpv_logkde(a), amp_mpv_e=float(np.std(mpvs)))


def cell_stats(g, ntrig, live):
    """all quantities of one (sub-run, region) cell; g = leading pads of the cell"""
    out = dict(amp_stats(g.amplitude.to_numpy()))
    out['sat_frac'] = float(g.saturated.mean()) if len(g) else np.nan
    out['ev_per_trig'] = len(g) / (ntrig * live) if ntrig else np.nan
    out['nfired_med'] = float(g.nfired.median()) if len(g) else np.nan
    clean = g[~g.saturated & ~g.trunc_left & ~g.trunc_right]
    for A in CUTS:
        c = clean[clean.amplitude >= A]
        t = int(A)
        for nm, x in [('tom', c.tom), ('time', c.time), ('dtm', c.tom - c.time), ('tot', c.tot),
                      ('ioa', c.integral / c.amplitude)]:
            s = boot_q(x.to_numpy())
            if nm == 'tom':
                out[f'n_{t}'] = s['n']
            keys = ['med', 'med_e'] if nm in ('dtm', 'tot', 'ioa') else ['med', 'med_e', 's68', 's68_e', 'p10', 'p90']
            for k in keys:
                out[f'{nm}_{k}_{t}'] = s[k]
    return out


# ----------------------------------------------------------------------------- loading
def subruns(run_dir):
    out = []
    for d in sorted(glob.glob(os.path.join(run_dir, '*_det2_*'))):
        nm = os.path.basename(d)
        m = re.match(r'(alignment|drift_scan|mesh_scan)_det2_(\d+)_(\d+)$', nm)
        if not m or not glob.glob(os.path.join(d, 'combined_hits_root', '*.root')):
            continue
        f = os.path.basename(glob.glob(os.path.join(d, 'combined_hits_root', '*.root'))[0])
        st = re.search(r'datrun_(\d{6})_(\d\d)H(\d\d)', f)
        out.append(dict(subrun=nm, kind=m.group(1).split('_')[0], mesh_set=int(m.group(2)),
                        drift_set=int(m.group(3)), start=f'{st.group(1)} {st.group(2)}:{st.group(3)}', dir=d))
    return sorted(out, key=lambda r: r['start'])


def hv_info(d, sv, t_lo, t_hi):
    h = pd.read_csv(os.path.join(d, 'hv_monitor.csv'))
    t = (pd.to_datetime(h.timestamp) - pd.to_datetime(h.timestamp.iloc[0])).dt.total_seconds().to_numpy()
    w = (t >= t_lo) & (t <= t_hi)
    if sv.intervals:
        w &= ~ps._in_intervals(t, sv.intervals)
    o = {}
    for nm, c in [('mesh', '1:0'), ('drift', '1:1')]:
        v = h[f'{c} vmon'].to_numpy(float); v0 = h[f'{c} v0'].to_numpy(float)
        o[f'v{nm}'] = float(v[w].mean()); o[f'v{nm}_rms'] = float(v[w].std())
        dev = np.abs(v - v0) > 3
        o[f'{nm}_dev_s'] = float(0.5 * (dev & (t >= t_lo) & (t <= t_hi)).sum())   # 2 Hz monitor
        o[f'{nm}_pwr0'] = int((h[f'{c} power'] == 0)[(t >= t_lo) & (t <= t_hi)].sum())
    o['n_sparks'] = len(sv.sparks); o['live'] = sv.live_fraction()
    return o


def build(out):
    cfgs = {k: qa.get_config(v) for k, v in
            [('alignment', 'det2_initial2'), ('drift', 'det2_driftscan2'), ('mesh', 'det2_meshscan2')]}
    cfg0 = cfgs['alignment']
    ct = sc.channel_table(cfg0)
    reg = pd.read_csv(REGIONS)
    lead_all, runs, log = [], [], []
    for s in subruns(cfg0.run_dir):
        cfg = cfgs[s['kind']]; d = s['dir']
        hd = os.path.join(d, 'combined_hits_root')
        drop = set(p2io.drop_pads_for(cfg, ct, hits_dir=hd))
        sv = ps.SparkVeto.from_csv(os.path.join(d, 'hv_monitor.csv'), cfg)
        bad = sv.vetoed_ids_from_hits(hd, ct.attrs['feus'], min_amp=cfg.MIN_AMP)
        parts = []
        for df in p2io.iter_hits(hd, BR, ct.attrs['feus'], progress=False, min_amp=cfg.MIN_AMP):
            h = pmap.attach_pads_to_hits(df, ct)
            h = h[h['mapped'] & h['pad_cx'].notna()]
            ev_all = h.eventId.unique()
            h = h[~h.channel_id.isin(drop) & ~h.eventId.isin(bad)]
            parts.append(h)
        h = pd.concat(parts, ignore_index=True)
        ts = h.trigger_timestamp_ns.to_numpy(np.int64) / 1e9
        nf = h.groupby('eventId').channel_id.nunique().rename('nfired')
        lead = h.sort_values('amplitude', ascending=False).drop_duplicates('eventId')
        lead = lead.merge(nf, left_on='eventId', right_index=True)
        lead = lead.merge(reg[['channel_id', 'reg8', 'reg24', 'reg48']], on='channel_id', how='left')
        nmiss = int(lead.reg24.isna().sum())
        lead = lead[lead.reg24.notna()]
        m3 = glob.glob(os.path.join(d, 'm3_tracking_root', '*.root'))
        ntrig = uproot.open(m3[0])['T'].num_entries if m3 else np.nan
        hv = hv_info(d, sv, float(ts.min()), float(ts.max()))
        dV = hv['vdrift'] - hv['vmesh']
        r = dict(s, **hv, dV=dV, E=dV / GAP_CM, ntrig=ntrig, n_lead=len(lead), n_vetoed=len(bad),
                 n_dropped_pads=len(drop), n_lead_noregion=nmiss, evid_max=int(h.eventId.max()))
        r.pop('dir')
        runs.append(r)
        lead = pd.DataFrame(dict(subrun=s['subrun'], kind=s['kind'], eventId=lead.eventId.to_numpy(),
                                 channel_id=lead.channel_id.astype(int).to_numpy(),
                                 reg8=lead.reg8.astype(int).to_numpy(), reg24=lead.reg24.astype(int).to_numpy(),
                                 pad_cx=lead.pad_cx.to_numpy(), pad_cy=lead.pad_cy.to_numpy(),
                                 amplitude=lead.amplitude.to_numpy(), time=lead.time.to_numpy(),
                                 tom=lead.time_of_max.to_numpy(), tot=lead.time_over_threshold.to_numpy(),
                                 integral=lead.integral.to_numpy(), saturated=lead.saturated.to_numpy(bool),
                                 trunc_left=lead.trunc_left.to_numpy(bool), trunc_right=lead.trunc_right.to_numpy(bool),
                                 significance=lead.significance.to_numpy(), nfired=lead.nfired.to_numpy()))
        lead_all.append(lead)
        msg = (f"{s['subrun']:26s} {s['start']} dV {dV:7.2f} V  E {dV / GAP_CM:7.1f} V/cm  trig {ntrig}  "
               f"lead {len(lead)}  vetoed ev {len(bad)}  sparks {hv['n_sparks']}  live {hv['live']:.3f}  "
               f"dropped pads {len(drop)}")
        print(msg, flush=True); log.append(msg)
    L = pd.concat(lead_all, ignore_index=True)
    R = pd.DataFrame(runs)
    L.to_csv(os.path.join(out, 'leading_pad_hits.csv.gz'), index=False)
    R.to_csv(os.path.join(out, 'subruns_hv.csv'), index=False)
    return L, R


def region_table(L, R):
    rows = []
    for r in R.itertuples():
        g0 = L[L.subrun == r.subrun]
        for scheme in ('all', 'reg8', 'reg24'):
            groups = [(-1, g0)] if scheme == 'all' else list(g0.groupby(scheme))
            for rid, g in groups:
                o = dict(subrun=r.subrun, kind=r.kind, mesh_set=r.mesh_set, drift_set=r.drift_set,
                         vmesh=r.vmesh, vdrift=r.vdrift, dV=r.dV, E=r.E, scheme=scheme, region=int(rid))
                o.update(cell_stats(g, r.ntrig, r.live))
                rows.append(o)
        print(f'  stats {r.subrun}', flush=True)
    return pd.DataFrame(rows)


def band_table(L, R, kinds, schemes=('all', 'reg8')):
    rows = []
    Rk = R[R.kind.isin(kinds)]
    for r in Rk.itertuples():
        g0 = L[(L.subrun == r.subrun) & ~L.saturated & ~L.trunc_left & ~L.trunc_right]
        for lo, hi in BANDS:
            gb = g0[(g0.amplitude >= lo) & (g0.amplitude < hi)]
            for scheme in schemes:
                groups = [(-1, gb)] if scheme == 'all' else list(gb.groupby(scheme))
                for rid, g in groups:
                    o = dict(subrun=r.subrun, kind=r.kind, vmesh=r.vmesh, dV=r.dV, E=r.E, band_lo=lo,
                             band_hi=hi, scheme=scheme, region=int(rid), n=len(g),
                             amp_med=float(g.amplitude.median()) if len(g) else np.nan)
                    for nm, x in [('tom', g.tom), ('time', g.time), ('dtm', g.tom - g.time), ('tot', g.tot),
                                  ('ioa', g.integral / g.amplitude)]:
                        s = boot_q(x.to_numpy())
                        o[f'{nm}_med'] = s['med']; o[f'{nm}_med_e'] = s['med_e']
                        if nm in ('tom', 'time'):
                            o[f'{nm}_s68'] = s['s68']; o[f'{nm}_s68_e'] = s['s68_e']
                    rows.append(o)
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------- fits
def wlin(x, y, e):
    """weighted y = a + b x; returns a, ea, b, eb, chi2, ndf (errors NOT rescaled)"""
    w = 1 / e ** 2
    X = np.c_[np.ones_like(x), x]
    C = np.linalg.inv(X.T @ (X * w[:, None]))
    p = C @ (X.T @ (w * y))
    chi2 = float((w * (y - X @ p) ** 2).sum())
    return p[0], np.sqrt(C[0, 0]), p[1], np.sqrt(C[1, 1]), chi2, len(x) - 2


def deviation_regression(T, scheme, q, proxy, cut, dVmin):
    """region deviation from the chamber value of quantity q vs a drift-time proxy"""
    D = T[(T.kind == 'drift') & (T.dV >= dVmin)]
    ch = D[D.scheme == 'all'].set_index('subrun')
    rr = D[D.scheme == scheme]
    out = []
    for rid, g in rr.groupby('region'):
        g = g.set_index('subrun')
        y = (g[f'{q}_{cut}'] - ch.loc[g.index, f'{q}_{cut}']).to_numpy()
        e = g[f'{q}_e_{cut}'].to_numpy()
        x = proxy.loc[g.index].to_numpy()
        k = np.isfinite(y) & np.isfinite(e) & (e > 0)
        if k.sum() < 4:
            continue
        a, ea, b, eb, chi2, ndf = wlin(x[k], y[k], e[k])
        c0 = np.average(y[k], weights=1 / e[k] ** 2)
        chi2c = float((((y[k] - c0) / e[k]) ** 2).sum())
        out.append(dict(region=rid, a=a, a_e=ea, b=b, b_e=eb, chi2=chi2, ndf=ndf, const=c0,
                        const_e=1 / np.sqrt((1 / e[k] ** 2).sum()), chi2_const=chi2c, ndf_const=int(k.sum()) - 1,
                        npts=int(k.sum())))
    return pd.DataFrame(out)


def collapse_fit(T, scheme, col, dVmin=20, nknots=8, s_free=True, scale_free=True):
    """A_r(dV) = g_r F(dV/s_r); F = cubic spline in dV (log amplitude); sum log g = sum log s = 0"""
    D = T[(T.kind == 'drift') & (T.scheme == scheme) & (T.dV >= dVmin)].dropna(subset=[col, col + '_e'])
    regs = np.sort(D.region.unique()); nr = len(regs)
    ri = np.searchsorted(regs, D.region.to_numpy())
    x = D.dV.to_numpy(); y = np.log(D[col].to_numpy()); e = D[col + '_e'].to_numpy() / D[col].to_numpy()
    kn = np.linspace(x.min() * 0.8, x.max() * 1.15, nknots)
    y0 = np.interp(kn, *zip(*sorted(zip(x, y))))

    def unpack(p):
        f = p[:nknots]; lg = np.r_[p[nknots:nknots + nr - 1], 0]; lg[-1] = -lg[:-1].sum()
        if s_free:
            ls = np.r_[p[nknots + nr - 1:], 0]; ls[-1] = -ls[:-1].sum()
        else:
            ls = np.zeros(nr)
        return f, lg, ls

    def res(p):
        f, lg, ls = unpack(p)
        u = x / np.exp(ls[ri])
        return (lg[ri] + CubicSpline(kn, f, bc_type='natural', extrapolate=True)(u) - y) / e
    p0 = np.r_[y0, np.zeros(nr - 1), np.zeros(nr - 1) if s_free else []]
    sol = least_squares(res, p0, x_scale='jac')
    chi2 = float((sol.fun ** 2).sum()); ndf = len(y) - len(p0)
    J = sol.jac
    try:
        cov = np.linalg.pinv(J.T @ J)
    except np.linalg.LinAlgError:
        cov = np.full((len(p0), len(p0)), np.nan)
    f, lg, ls = unpack(sol.x)
    # errors of all nr log s (last one is minus the sum)
    out = dict(regs=regs, chi2=chi2, ndf=ndf, g=np.exp(lg), s=np.exp(ls), knots=kn, f=f)
    if s_free:
        i0 = nknots + nr - 1
        A = np.vstack([np.eye(nr - 1), -np.ones(nr - 1)])
        cs = A @ cov[i0:, i0:] @ A.T
        out['s_e'] = np.exp(ls) * np.sqrt(np.diag(cs))
        out['cov_ls'] = cs
    return out


class PhysFit:
    """sigma68_r = sqrt(s0^2 + (k T_r)^2), median_r = t0_r + c T_r,
    T_r = d_r / v_d((dV - V0)/d_r)"""

    def __init__(self, T, scheme, vd, cut=300, dVmin=50, use=('s68', 'med')):
        D = T[(T.kind == 'drift') & (T.scheme == scheme) & (T.dV >= dVmin)]
        D = D.dropna(subset=[f'tom_med_{cut}', f'tom_s68_{cut}'])
        D = D[(D[f'tom_med_e_{cut}'] > 0) & (D[f'tom_s68_e_{cut}'] > 0)]
        self.regs = np.sort(D.region.unique()); self.nr = len(self.regs)
        self.ri = np.searchsorted(self.regs, D.region.to_numpy())
        self.dV = D.dV.to_numpy()
        self.m = D[f'tom_med_{cut}'].to_numpy(); self.me = D[f'tom_med_e_{cut}'].to_numpy()
        self.s = D[f'tom_s68_{cut}'].to_numpy(); self.se = D[f'tom_s68_e_{cut}'].to_numpy()
        self.vd = vd; self.use = use

    def Tr(self, d_mm, V0, ri=None, dV=None):
        ri = self.ri if ri is None else ri; dV = self.dV if dV is None else dV
        dcm = d_mm[ri] / 10
        return dcm / self.vd((dV - V0) / dcm)

    def fit(self, V0_free=False, d_common=False, d_fixed=None):
        nr = self.nr
        nd = 0 if d_fixed is not None else (1 if d_common else nr)
        names = ['k', 'c', 's0'] + (['V0'] if V0_free else []) + [f'd{i}' for i in range(nd)] + \
                [f't0_{i}' for i in range(nr)]

        def unpack(p):
            k, c, s0 = p[:3]; i = 3
            V0 = 0.0
            if V0_free:
                V0 = p[i]; i += 1
            if d_fixed is not None:
                d = np.full(nr, d_fixed)
            elif d_common:
                d = np.full(nr, p[i]); i += 1
            else:
                d = p[i:i + nr]; i += nr
            t0 = p[i:i + nr]
            return k, c, s0, V0, d, t0

        def res(p):
            k, c, s0, V0, d, t0 = unpack(p)
            T = self.Tr(d, V0)
            r = []
            if 's68' in self.use:
                r.append((np.sqrt(s0 ** 2 + (k * T) ** 2) - self.s) / self.se)
            if 'med' in self.use:
                r.append((t0[self.ri] + c * T - self.m) / self.me)
            return np.concatenate(r)
        t0g = np.array([self.m[self.ri == i].min() for i in range(nr)]) - 30
        p0 = [0.3, 0.5, 15.0] + ([0.0] if V0_free else []) + ([4.0] * nd) + list(t0g)
        lo = [0, 0, 0] + ([-60] if V0_free else []) + [1.0] * nd + [-np.inf] * nr
        hi = [5, 5, 200] + ([60] if V0_free else []) + [12.0] * nd + [np.inf] * nr
        sol = least_squares(res, p0, bounds=(lo, hi), x_scale='jac')
        chi2 = float((sol.fun ** 2).sum()); ndf = len(sol.fun) - len(p0)
        cov = np.linalg.pinv(sol.jac.T @ sol.jac)
        err = np.sqrt(np.diag(cov))
        sd = np.sqrt(np.outer(np.diag(cov), np.diag(cov)))
        corr = cov / np.where(sd > 0, sd, np.nan)
        k, c, s0, V0, d, t0 = unpack(sol.x)
        out = dict(names=names, p=sol.x, err=err, cov=cov, corr=corr, chi2=chi2, ndf=ndf, k=k, c=c, s0=s0,
                   V0=V0, d=d, t0=t0, n=len(sol.fun), active=sol.active_mask)
        if nd == nr:
            i0 = names.index('d0')
            w = np.full(nr, 1 / nr)
            Cd = cov[i0:i0 + nr, i0:i0 + nr]
            out['d_mean_e'] = float(np.sqrt(w @ Cd @ w))
            out['d_e'] = np.sqrt(np.diag(Cd))
            # relative d_r / d_mean errors (removes the common scale)
            A = np.eye(nr) - w
            out['drel_e'] = np.sqrt(np.diag(A @ Cd @ A.T)) / d.mean()
            # correlation of k and c with mean d
            vk = cov[0, i0:i0 + nr] @ w; vc = cov[1, i0:i0 + nr] @ w
            out['corr_k_dmean'] = vk / (err[0] * out['d_mean_e'])
            out['corr_c_dmean'] = vc / (err[1] * out['d_mean_e'])
        return out


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--rebuild', action='store_true')
    args = ap.parse_args()
    cfg0 = qa.get_config('det2_initial2')
    out = os.path.join(qa.ANALYSIS_ROOT, 'det2', cfg0.RUN, '44_deform_scan_regions')
    os.makedirs(out, exist_ok=True)
    fL = os.path.join(out, 'leading_pad_hits.csv.gz'); fR = os.path.join(out, 'subruns_hv.csv')
    if args.rebuild or not (os.path.isfile(fL) and os.path.isfile(fR)):
        L, R = build(out)
    else:
        L = pd.read_csv(fL); R = pd.read_csv(fR)
    fT = os.path.join(out, 'region_stats.csv')
    if args.rebuild or not os.path.isfile(fT):
        T = region_table(L, R); T.to_csv(fT, index=False)
    else:
        T = pd.read_csv(fT)
    fB = os.path.join(out, 'band_stats.csv')
    if args.rebuild or not os.path.isfile(fB):
        Bt = band_table(L, R, ('mesh', 'drift', 'alignment')); Bt.to_csv(fB, index=False)
    else:
        Bt = pd.read_csv(fB)
    vd, gas = vdrift_func()
    S = []   # summary lines

    def P(s=''):
        print(s, flush=True); S.append(s)

    # ---------------- sub-run / HV table
    P('=== Sub-runs (chronological), HV from hv_monitor.csv vmon (mesh 1:0, drift 1:1), mean over the trigger')
    P('    span outside spark-veto windows; dev = seconds with |vmon - v0| > 3 V inside the span')
    for r in R.itertuples():
        P(f'  {r.start} {r.subrun:26s} mesh {r.vmesh:7.2f} drift {r.vdrift:7.2f} dV {r.dV:7.2f} V E {r.E:7.1f} V/cm '
          f'| trig {r.ntrig} lead {r.n_lead} ({r.n_lead / r.ntrig:.3f}/trig) | sparks {r.n_sparks} live {r.live:.3f} '
          f'mesh dev {r.mesh_dev_s:.0f}s drift dev {r.drift_dev_s:.0f}s | vetoed ev {r.n_vetoed} hot/noisy pads {r.n_dropped_pads}')
    P(f'  leading pads outside regions.csv: {int(R.n_lead_noregion.sum())} in total')

    # ---------------- 5a chamber-wide drift scan
    ch = T[(T.scheme == 'all')].copy()
    ch['Tgas'] = GAP_CM / vd(np.clip(ch.dV, 1e-3, None) / GAP_CM)
    ch.loc[ch.dV < 1, 'Tgas'] = np.nan
    dr = ch[ch.kind == 'drift'].sort_values('dV')
    P('\n=== 5a  Chamber-wide drift scan (leading pad, clean, tom = time_of_max)')
    P('   E[V/cm]  T_gas[ns] | A>=300: n  med(tom)  s68(tom)  med(time) s68(time) med(tom-time) TOT  int/amp | A>=150: n med(tom) s68(tom) | MPV  ev/trig')
    for r in dr.itertuples():
        P(f'   {r.E:7.1f} {r.Tgas:7.1f} | {r.n_300:5d} {r.tom_med_300:6.1f}+-{r.tom_med_e_300:3.1f} {r.tom_s68_300:5.1f}+-{r.tom_s68_e_300:3.1f} '
          f'{r.time_med_300:6.1f}+-{r.time_med_e_300:3.1f} {r.time_s68_300:5.1f}+-{r.time_s68_e_300:3.1f} {r.dtm_med_300:6.1f}+-{r.dtm_med_e_300:3.1f} '
          f'{r.tot_med_300:5.0f} {r.ioa_med_300:4.2f} | {r.n_150:5d} {r.tom_med_150:6.1f}+-{r.tom_med_e_150:3.1f} {r.tom_s68_150:5.1f}+-{r.tom_s68_e_150:3.1f} '
          f'| {r.amp_mpv:5.0f}+-{r.amp_mpv_e:3.0f} {r.ev_per_trig:.3f}')
    al = ch[ch.kind == 'alignment'].iloc[0]
    P(f'   alignment 430/730 (1 h): n {al.n_300} med(tom) {al.tom_med_300:.1f}+-{al.tom_med_e_300:.1f} s68 {al.tom_s68_300:.1f}+-{al.tom_s68_e_300:.1f} '
      f'med(time) {al.time_med_300:.1f} MPV {al.amp_mpv:.0f}')
    fit_rows = []
    for cut in (300, 150):
        for dmin in (50, 100, 200):
            k = (dr.dV >= dmin)
            for q, qe in [('tom_med', 'tom_med_e'), ('time_med', 'time_med_e'), ('tom_s68', 'tom_s68_e'),
                          ('time_s68', 'time_s68_e'), ('dtm_med', 'dtm_med_e')]:
                y = dr[f'{q}_{cut}'][k].to_numpy(); e = dr[f'{qe}_{cut}'][k].to_numpy(); x = dr.Tgas[k].to_numpy()
                a, ea, b, eb, chi2, ndf = wlin(x, y, e)
                fit_rows.append(dict(cut=cut, dVmin=dmin, q=q, a=a, a_e=ea, b=b, b_e=eb, chi2=chi2, ndf=ndf))
    FR = pd.DataFrame(fit_rows); FR.to_csv(os.path.join(out, 'chamber_vs_Tgas_linear.csv'), index=False)
    P('   Linear fits q = a + b * T_gas(E) (d = 4 mm; errors not rescaled; chi2/ndf given):')
    for r in FR.itertuples():
        P(f'     A>={r.cut} dV>={r.dVmin:3d}V  {r.q:9s}: a {r.a:7.1f}+-{r.a_e:4.1f}  b {r.b:6.3f}+-{r.b_e:5.3f}  chi2/ndf {r.chi2:6.1f}/{r.ndf}')
    # quadrature model for s68 (chamber)
    for cut in (300, 150):
        k = dr.dV >= 50
        x = dr.Tgas[k].to_numpy(); y = dr[f'tom_s68_{cut}'][k].to_numpy(); e = dr[f'tom_s68_e_{cut}'][k].to_numpy()
        sol = least_squares(lambda p: (np.sqrt(p[0] ** 2 + (p[1] * x) ** 2) - y) / e, [15, 0.3])
        cv = np.linalg.pinv(sol.jac.T @ sol.jac)
        P(f'   s68(tom) = sqrt(s0^2 + (k T_gas)^2), A>={cut}, dV>=50: s0 {sol.x[0]:.1f}+-{np.sqrt(cv[0, 0]):.1f} ns  '
          f'k {sol.x[1]:.3f}+-{np.sqrt(cv[1, 1]):.3f}  chi2/ndf {(sol.fun ** 2).sum():.1f}/{len(x) - 2}  '
          f'(uniform-depth ionisation alone would give k = 0.34)')

    # time-walk check on the drift scan: fixed amplitude bands
    P('   Fixed amplitude band (chamber, clean leading pads): med(tom) [ns] per band vs E')
    bd = Bt[(Bt.kind == 'drift') & (Bt.scheme == 'all')]
    for (lo, hi), g in bd.groupby(['band_lo', 'band_hi']):
        g = g.sort_values('dV')
        P(f'     {lo:4d}-{hi:4d} ADC: ' + ' '.join(f'{e:.0f}:{m:.0f}({n})' for e, m, n in zip(g.E, g.tom_med, g.n)))

    # ---------------- figure 1
    fig, axs = plt.subplots(2, 3, figsize=(17, 9.5))
    Eg = np.linspace(30, 1400, 300); Tg = GAP_CM / vd(Eg)
    for cut, col in [(300, '#c0392b'), (150, '#2471a3')]:
        axs[0, 0].errorbar(dr.E, dr[f'tom_med_{cut}'], dr[f'tom_med_e_{cut}'], fmt='o-', color=col, ms=4, label=f'time_of_max, A>={cut}')
        axs[0, 0].errorbar(dr.E, dr[f'time_med_{cut}'], dr[f'time_med_e_{cut}'], fmt='s--', color=col, ms=4, mfc='none', label=f'time (30% CFD), A>={cut}')
        axs[0, 1].errorbar(dr.E, dr[f'tom_s68_{cut}'], dr[f'tom_s68_e_{cut}'], fmt='o-', color=col, ms=4, label=f'time_of_max, A>={cut}')
        axs[0, 1].errorbar(dr.E, dr[f'time_s68_{cut}'], dr[f'time_s68_e_{cut}'], fmt='s--', color=col, ms=4, mfc='none', label=f'time, A>={cut}')
        axs[0, 2].errorbar(dr.E, dr[f'dtm_med_{cut}'], dr[f'dtm_med_e_{cut}'], fmt='o-', color=col, ms=4, label=f'A>={cut}')
    r = FR[(FR.cut == 300) & (FR.dVmin == 50) & (FR.q == 'tom_med')].iloc[0]
    axs[0, 0].plot(Eg, r.a + r.b * Tg, ':', color='k', label=f'{r.a:.0f} + {r.b:.2f} d/v_d(E) (fit, E>=125)')
    axs[0, 0].set(xlabel='drift field (Vdrift - Vmesh)/0.4 cm [V/cm]', ylabel='median [ns]', title='Leading-pad median time vs drift field')
    axs[0, 1].plot(Eg, 0.34 * Tg, ':', color='k', label='0.34 d/v_d (uniform depth, no floor)')
    axs[0, 1].set(xlabel='drift field [V/cm]', ylabel='sigma68 [ns]', title='Leading-pad time spread vs drift field')
    axs[0, 2].set(xlabel='drift field [V/cm]', ylabel='median(time_of_max - time) [ns]', title='Peak minus 30% CFD time')
    axs[1, 0].plot(Eg, Tg, 'k-'); axs[1, 0].set(xlabel='drift field [V/cm]', ylabel='d / v_d [ns] (d = 4 mm)',
                                               title='Magboltz Ar/iC4H10 95/5 full drift time', ylim=(0, 250))
    for cut, col in [(300, '#c0392b'), (150, '#2471a3')]:
        k = dr.dV >= 25
        axs[1, 1].errorbar(dr.Tgas[k], dr[f'tom_med_{cut}'][k], dr[f'tom_med_e_{cut}'][k], fmt='o', color=col, ms=4, label=f'median tom A>={cut}')
        axs[1, 1].errorbar(dr.Tgas[k], dr[f'tom_s68_{cut}'][k] + 400, dr[f'tom_s68_e_{cut}'][k], fmt='^', color=col, ms=4, mfc='none', label=f'sigma68 tom + 400 ns, A>={cut}')
    axs[1, 1].set(xlabel='T_gas = d / v_d(E) [ns]', ylabel='[ns]', title='Median and width vs predicted drift time')
    axs[1, 2].errorbar(dr.E, dr.amp_mpv, dr.amp_mpv_e, fmt='o-', color='#c0392b', ms=4, label='MPV (log-KDE mode)')
    axs[1, 2].errorbar(dr.E, dr.amp_med, dr.amp_med_e, fmt='s--', color='#2471a3', ms=4, label='median')
    axs[1, 2].set(xlabel='drift field [V/cm]', ylabel='leading-pad amplitude [ADC]', title='Leading-pad amplitude vs drift field (mesh 430 V)')
    for a in axs.ravel():
        a.grid(color=GRID); a.legend(frameon=False, fontsize=7)
    fig.suptitle('det2 drift scan 24 Nov 2026 (mesh 430 V): chamber-wide leading-pad timing')
    fig.tight_layout(); fig.savefig(os.path.join(out, 'drift_chamber.png'), dpi=120, facecolor='white'); plt.close(fig)

    # ---------------- 5b region deviations
    P('\n=== 5b  Region-to-region pattern vs drift field (drift scan, A>=300 unless stated)')
    P('   deviation D_r(E) = q_r(E) - q_chamber(E), regressed vs a drift-time proxy x(E): D_r = a_r + b_r x')
    drs = dr.set_index('subrun')
    zc0 = L[(L.subrun == 'drift_scan_det2_430_430') & ~L.saturated & ~L.trunc_left & ~L.trunc_right].tom.median()
    proxies = {'Tgas': drs.Tgas, 'chamber_med_tom-t0(0V)': drs.tom_med_300 - zc0}
    P(f'   proxies: T_gas = 4 mm / v_d(E) (Magboltz); chamber median tom minus the chamber zero-drift t0 '
      f'({zc0:.1f} ns, 430/430 point) -> with the 2nd proxy the intercept a_r is the field-independent offset')
    DEV = {}
    for scheme in ('reg24', 'reg8'):
        for q in ('tom_med', 'tom_s68', 'time_med'):
            for pn, px in proxies.items():
                for cut in (300, 150):
                    for dmin in (50, 25):
                        dv = deviation_regression(T, scheme, q, px, cut, dmin)
                        if not len(dv):
                            continue
                        DEV[(scheme, q, pn, cut, dmin)] = dv
                        dof = len(dv) - 1
                        chi_b = float(((dv.b / dv.b_e) ** 2).sum())
                        chi_a = float(((dv.a / dv.a_e) ** 2).sum())
                        chi_c = float(((dv.const / dv.const_e) ** 2).sum())
                        fitq = dv.chi2.sum() / dv.ndf.sum()
                        if pn == 'Tgas' and q == 'tom_med' and cut == 300 and dmin == 50:
                            pass
                        if dmin == 50 or (q == 'tom_med' and cut == 300):
                            P(f'   {scheme} {q:8s} A>={cut} dV>={dmin:2d} proxy {pn:20s}: slopes b_r all 0? chi2 {chi_b:6.1f}/{dof} '
                              f'(p {sst.chi2.sf(chi_b, dof):.2g}); intercepts a_r all 0? chi2 {chi_a:6.1f}/{dof} '
                              f'(p {sst.chi2.sf(chi_a, dof):.2g}); constant offsets c_r all 0? chi2 {chi_c:6.1f}/{dof}; '
                              f'linear-fit chi2/ndf pooled {fitq:.2f}; constant-only chi2/ndf pooled '
                              f'{dv.chi2_const.sum() / dv.ndf_const.sum():.2f}')
    pd.concat([v.assign(scheme=k[0], q=k[1], proxy=k[2], cut=k[3], dVmin=k[4]) for k, v in DEV.items()]) \
        .to_csv(os.path.join(out, 'region_deviation_regression.csv'), index=False)
    # reproducibility: alignment run vs drift_scan 430_730 region pattern
    for scheme in ('reg24', 'reg8'):
        a1 = T[(T.kind == 'alignment') & (T.scheme == scheme)].set_index('region')
        a2 = T[(T.subrun == 'drift_scan_det2_430_730') & (T.scheme == scheme)].set_index('region')
        c1 = T[(T.kind == 'alignment') & (T.scheme == 'all')].iloc[0]
        c2 = T[(T.subrun == 'drift_scan_det2_430_730') & (T.scheme == 'all')].iloc[0]
        d1 = a1.tom_med_300 - c1.tom_med_300; d2 = a2.tom_med_300 - c2.tom_med_300
        e = np.sqrt(a1.tom_med_e_300 ** 2 + a2.tom_med_e_300 ** 2)
        chi = float((((d1 - d2) / e) ** 2).sum())
        rho = np.corrcoef(d1, d2)[0, 1]
        P(f'   repeatability {scheme}: alignment (1 h) vs drift_scan 430_730 (30 min) region offsets of med(tom): '
          f'diff chi2 {chi:.1f}/{len(d1)}, Pearson r {rho:.2f}; rms of offsets {d1.std():.1f} / {d2.std():.1f} ns')
    # pooled region offsets (weighted mean over E >= 300 V/cm and over E 62-250 V/cm)
    P('   reg24 offsets of med(tom) [ns] from the chamber, weighted mean in field bands (A>=300):')
    R24 = T[(T.kind == 'drift') & (T.scheme == 'reg24')].merge(dr[['subrun', 'tom_med_300']].rename(columns={'tom_med_300': 'chm'}), on='subrun')
    R24['dev'] = R24.tom_med_300 - R24.chm
    sig_rows = []
    for rid, g in R24.groupby('region'):
        o = dict(region=rid)
        for nm, lo, hi in [('low', 50, 100), ('mid', 100, 250), ('high', 250, 600)]:
            gg = g[(g.dV >= lo) & (g.dV < hi)]
            w = 1 / gg.tom_med_e_300 ** 2
            o[f'dev_{nm}'] = float((gg.dev * w).sum() / w.sum()); o[f'dev_{nm}_e'] = float(1 / np.sqrt(w.sum()))
        sig_rows.append(o)
    SR = pd.DataFrame(sig_rows)
    dvT = DEV[('reg24', 'tom_med', 'Tgas', 300, 50)].set_index('region')
    dvS = DEV[('reg24', 'tom_s68', 'Tgas', 300, 50)].set_index('region')
    SR = SR.set_index('region').join(dvT[['a', 'a_e', 'b', 'b_e']]).join(dvS[['a', 'a_e', 'b', 'b_e']], rsuffix='_s68')
    SR.to_csv(os.path.join(out, 'reg24_offsets_summary.csv'))
    P('   reg  dev(E 125-250) dev(250-625) dev(625-1500) | a_r [ns]  b_r (vs T_gas)  | s68: a_r  b_r')
    for rid, r in SR.iterrows():
        P(f'   {rid:3d}  {r.dev_low:+5.1f}+-{r.dev_low_e:3.1f}  {r.dev_mid:+5.1f}+-{r.dev_mid_e:3.1f}  {r.dev_high:+5.1f}+-{r.dev_high_e:3.1f} | '
          f'{r.a:+5.1f}+-{r.a_e:3.1f} {r.b:+6.3f}+-{r.b_e:5.3f} | {r.a_s68:+5.1f}+-{r.a_e_s68:3.1f} {r.b_s68:+6.3f}+-{r.b_e_s68:5.3f}')
    rho = np.corrcoef(SR.dev_low, SR.dev_high)[0, 1]
    P(f'   Pearson r between reg24 offsets at E 125-250 and at E 625-1500 V/cm: {rho:.2f}; '
      f'rms {SR.dev_low.std():.1f} vs {SR.dev_high.std():.1f} ns (mean stat err {SR.dev_low_e.mean():.1f}, {SR.dev_high_e.mean():.1f})')

    # zero-drift reference: the 430/430 point (dV ~ 0.5 V) only shows charge created at / below the
    # mesh (no drift), so its per-region median time_of_max is the region's electronics + trigger t0.
    P('\n   Zero-drift reference (drift_scan_det2_430_430, dV ~ 0.5 V): per-region t0 = median tom of clean leading pads')
    z0 = L[(L.subrun == 'drift_scan_det2_430_430') & ~L.saturated & ~L.trunc_left & ~L.trunc_right]
    ZR = []
    for amin in (0, 100):
        zz = z0[z0.amplitude >= amin]
        zc = boot_q(zz.tom.to_numpy())
        P(f'   A>={amin}: chamber n {zc["n"]} med {zc["med"]:.1f}+-{zc["med_e"]:.1f} s68 {zc["s68"]:.1f} ns, '
          f'median amplitude {zz.amplitude.median():.0f} ADC')
        for scheme in ('reg8', 'reg24'):
            for rid, g in zz.groupby(scheme):
                s = boot_q(g.tom.to_numpy())
                ZR.append(dict(amin=amin, scheme=scheme, region=rid, n=s['n'], t0=s['med'], t0_e=s['med_e'],
                               off0=s['med'] - zc['med'], ch0=zc['med'], ch0_e=zc['med_e']))
    ZR = pd.DataFrame(ZR)
    LAM = []
    for scheme in ('reg8', 'reg24'):
        ZRs = ZR[(ZR.amin == 0) & (ZR.scheme == scheme)].set_index('region')
        Rg = T[(T.kind == 'drift') & (T.scheme == scheme) & (T.dV >= 25)][['subrun', 'region', 'dV', 'E', 'tom_med_300', 'tom_med_e_300']]
        Rg = Rg.merge(dr[['subrun', 'tom_med_300', 'tom_med_e_300']].rename(columns={'tom_med_300': 'chm', 'tom_med_e_300': 'chm_e'}), on='subrun')
        Rg['t0'] = ZRs.t0.reindex(Rg.region).to_numpy(); Rg['t0_e'] = ZRs.t0_e.reindex(Rg.region).to_numpy()
        ch0 = ZRs.ch0.iloc[0]; ch0e = ZRs.ch0_e.iloc[0]
        Rg['Dr'] = Rg.tom_med_300 - Rg.t0; Rg['Dc'] = Rg.chm - ch0
        Rg['lam'] = Rg.Dr / Rg.Dc
        # error: region median + region t0 (dominant) ; chamber terms are small
        Rg['lam_e'] = np.abs(Rg.lam) * np.sqrt((Rg.tom_med_e_300 ** 2 + Rg.t0_e ** 2) / Rg.Dr ** 2 + (Rg.chm_e ** 2 + ch0e ** 2) / Rg.Dc ** 2)
        for rid, g in Rg.groupby('region'):
            o = dict(scheme=scheme, region=rid, n0=int(ZRs.n[rid]), off0=ZRs.off0[rid], off0_e=ZRs.t0_e[rid])
            for nm, lo, hi in [('low', 25, 100), ('mid', 100, 250), ('high', 250, 600)]:
                gg = g[(g.dV >= lo) & (g.dV < hi)]
                w = 1 / gg.lam_e ** 2
                o[f'lam_{nm}'] = float((gg.lam * w).sum() / w.sum()); o[f'lam_{nm}_e'] = float(1 / np.sqrt(w.sum()))
                o[f'lam_{nm}_chi2'] = float((((gg.lam - o[f'lam_{nm}']) / gg.lam_e) ** 2).sum()); o[f'lam_{nm}_n'] = len(gg)
                w2 = 1 / gg.tom_med_e_300 ** 2
                o[f'off_{nm}'] = float(((gg.tom_med_300 - gg.chm) * w2).sum() / w2.sum()); o[f'off_{nm}_e'] = float(1 / np.sqrt(w2.sum()))
            LAM.append(o)
    LAM = pd.DataFrame(LAM); LAM.to_csv(os.path.join(out, 'zero_drift_t0_and_lambda.csv'), index=False)
    ZR.to_csv(os.path.join(out, 'zero_drift_t0.csv'), index=False)
    P('   Per region: offset of t0 (0 V) from chamber t0 vs offset of med(tom) at E 625-1500 V/cm; drift-length ratio')
    P('   lambda_r = (med_r(E) - t0_r)/(med_ch(E) - t0_ch) averaged in dV 25-100 (E 62-250), 100-250 (250-625), 250-525 (625-1312 V/cm)')
    for scheme in ('reg8', 'reg24'):
        G = LAM[LAM.scheme == scheme].set_index('region')
        chi_same = float((((G.off_high - G.off0) / np.sqrt(G.off0_e ** 2 + G.off_high_e ** 2)) ** 2).sum())
        chi_l1 = float((((G.lam_high - 1) / G.lam_high_e) ** 2).sum())
        P(f'   {scheme}: high-field offsets == t0 offsets? chi2 {chi_same:.1f}/{len(G)} (p {sst.chi2.sf(chi_same, len(G)):.2g}); '
          f'Pearson r(off0, off_high) {np.corrcoef(G.off0, G.off_high)[0, 1]:.2f}; rms off0 {G.off0.std():.1f} ns, rms off_high {G.off_high.std():.1f} ns; '
          f'lambda_high all 1? chi2 {chi_l1:.1f}/{len(G) - 1} (p {sst.chi2.sf(chi_l1, len(G) - 1):.2g}); rms lambda high/mid/low '
          f'{G.lam_high.std():.3f}/{G.lam_mid.std():.3f}/{G.lam_low.std():.3f}; '
          f'Pearson r(lam_high, lam_low) {np.corrcoef(G.lam_high, G.lam_low)[0, 1]:.2f}')
        P('   reg  n0   t0 offset    off(E>=625)  | lam_low (E 62-250)  lam_mid (250-625)  lam_high (625-1312) [chi2/n within band]')
        for rid, r in G.iterrows():
            P(f'   {rid:3d} {r.n0:4d} {r.off0:+6.1f}+-{r.off0_e:3.1f} {r.off_high:+6.1f}+-{r.off_high_e:3.1f} | '
              f'{r.lam_low:5.2f}+-{r.lam_low_e:4.2f} [{r.lam_low_chi2:4.1f}/{r.lam_low_n}]  {r.lam_mid:5.2f}+-{r.lam_mid_e:4.2f} [{r.lam_mid_chi2:4.1f}/{r.lam_mid_n}]  '
              f'{r.lam_high:5.3f}+-{r.lam_high_e:5.3f} [{r.lam_high_chi2:4.1f}/{r.lam_high_n}]')
    # high-field gap estimate: at each E >= 625 V/cm solve d_r/v(dV/d_r) = lambda_r(E) * 4mm/v(dV/4mm)
    # (Magboltz only used where it describes the chamber data; chamber mean gap ASSUMED 4.0 mm)
    from scipy.optimize import brentq
    DR = []
    for scheme in ('reg8', 'reg24'):
        ZRs = ZR[(ZR.amin == 0) & (ZR.scheme == scheme)].set_index('region')
        Rg = T[(T.kind == 'drift') & (T.scheme == scheme) & (T.dV >= 250)][['subrun', 'region', 'dV', 'tom_med_300', 'tom_med_e_300']]
        Rg = Rg.merge(dr[['subrun', 'tom_med_300']].rename(columns={'tom_med_300': 'chm'}), on='subrun')
        ch0 = ZRs.ch0.iloc[0]
        for rid, g in Rg.groupby('region'):
            ds, es = [], []
            for r in g.itertuples():
                lam = (r.tom_med_300 - ZRs.t0[rid]) / (r.chm - ch0)
                lam_e = np.hypot(r.tom_med_e_300, ZRs.t0_e[rid]) / (r.chm - ch0)
                Tc = GAP_CM / vd(r.dV / GAP_CM)
                f = lambda d, lm: (d / vd(r.dV / d)) - lm * Tc  # noqa: E731
                try:
                    d0 = brentq(f, 0.1, 1.2, args=(lam,)); d1 = brentq(f, 0.1, 1.2, args=(lam + lam_e,))
                except ValueError:
                    continue
                ds.append(d0 * 10); es.append(abs(d1 - d0) * 10)
            ds = np.array(ds); es = np.array(es); w = 1 / es ** 2
            DR.append(dict(scheme=scheme, region=rid, d_mm=float((ds * w).sum() / w.sum()), d_e=float(1 / np.sqrt(w.sum())),
                           chi2=float((((ds - (ds * w).sum() / w.sum()) / es) ** 2).sum()), npts=len(ds)))
    DR = pd.DataFrame(DR); DR.to_csv(os.path.join(out, 'gap_from_lambda_highfield.csv'), index=False)
    P('   Gap estimate from lambda_r at E 625-1312 V/cm (Magboltz v_d, chamber-mean gap ASSUMED 4.0 mm; error = timing stat only):')
    for scheme in ('reg8', 'reg24'):
        G_ = DR[DR.scheme == scheme]
        P(f'   {scheme}: ' + '  '.join(f'{int(r.region)}:{r.d_mm:.2f}+-{r.d_e:.2f}' + (f'[chi2 {r.chi2:.0f}/{r.npts - 1}]' if r.chi2 > 2 * (r.npts - 1) else '')
                                    for r in G_.itertuples()) + f'   (rms {G_.d_mm.std():.2f} mm)')

    fig, axs = plt.subplots(1, 3, figsize=(18, 5.2))
    G = LAM[LAM.scheme == 'reg24'].set_index('region')
    axs[0].errorbar(G.off0, G.off_high, G.off_high_e, G.off0_e, fmt='o', color='#2471a3')
    for rid, r in G.iterrows():
        axs[0].annotate(str(rid), (r.off0, r.off_high), fontsize=7)
    lim = [min(G.off0.min(), G.off_high.min()) - 3, max(G.off0.max(), G.off_high.max()) + 3]
    axs[0].plot(lim, lim, 'k:', lw=1)
    axs[0].set(xlabel='t0 offset at zero drift field (430/430) [ns]', ylabel='med(tom) offset at E 625-1312 V/cm [ns]',
               title='reg24: electronics/trigger offset vs high-field offset')
    axs[1].errorbar(G.lam_high, G.lam_low, G.lam_low_e, G.lam_high_e, fmt='o', color='#c0392b')
    for rid, r in G.iterrows():
        axs[1].annotate(str(rid), (r.lam_high, r.lam_low), fontsize=7)
    axs[1].set(xlabel='lambda_r at E 625-1312 V/cm', ylabel='lambda_r at E 62-250 V/cm',
               title='reg24 drift-time ratio region/chamber (t0-subtracted)')
    for a in axs[:2]:
        a.grid(color=GRID)
    LAMFIG = (fig, axs, G)

    # figure 2: region deviation maps + slopes
    reg = pd.read_csv(REGIONS)
    fig, axs = plt.subplots(2, 3, figsize=(18, 10))

    def pmap_(ax, vals, title, cmap='coolwarm', sym=True):
        v = reg.reg24.map(vals)
        vm = np.nanmax(np.abs(v)) if sym else None
        sc_ = ax.scatter(reg.pad_cx, reg.pad_cy, c=v, s=14, marker='s', cmap=cmap,
                         vmin=-vm if sym else None, vmax=vm if sym else None)
        for rid, g in reg.groupby('reg24'):
            ax.text(g.pad_cx.mean(), g.pad_cy.mean(), str(rid), fontsize=7, ha='center', va='center')
        plt.colorbar(sc_, ax=ax); ax.set(title=title, xlabel='pad x [mm]', ylabel='pad y [mm]'); ax.set_aspect('equal')
    lf, laxs, LG = LAMFIG
    pmap_(laxs[2], LG.lam_high - 1, 'reg24 lambda_r - 1 at E 625-1312 V/cm', sym=True)
    lf.tight_layout(); lf.savefig(os.path.join(out, 'zero_drift_lambda.png'), dpi=120, facecolor='white'); plt.close(lf)
    pmap_(axs[0, 0], SR.dev_high, 'med(tom) offset, E 625-1500 V/cm [ns]')
    pmap_(axs[0, 1], SR.dev_low, 'med(tom) offset, E 125-250 V/cm [ns]')
    pmap_(axs[0, 2], SR.b, 'slope b_r of offset vs T_gas [ns/ns]')
    axs[1, 0].errorbar(SR.dev_high, SR.dev_low, SR.dev_low_e, SR.dev_high_e, fmt='o', color='#2471a3')
    for rid, r in SR.iterrows():
        axs[1, 0].annotate(str(rid), (r.dev_high, r.dev_low), fontsize=7)
    lim = [min(SR.dev_high.min(), SR.dev_low.min()) - 2, max(SR.dev_high.max(), SR.dev_low.max()) + 2]
    axs[1, 0].plot(lim, lim, 'k:', lw=1)
    axs[1, 0].set(xlabel='offset at E 625-1500 V/cm [ns]', ylabel='offset at E 125-250 V/cm [ns]', title='reg24 med(tom) offsets: low vs high field')
    axs[1, 1].errorbar(SR.index, SR.b, SR.b_e, fmt='o', color='#c0392b', label='median tom')
    axs[1, 1].errorbar(SR.index + 0.25, SR.b_s68, SR.b_e_s68, fmt='s', color='#2471a3', label='sigma68 tom')
    axs[1, 1].axhline(0, color='k', lw=1)
    axs[1, 1].set(xlabel='reg24 region', ylabel='slope b_r vs T_gas [ns per ns]', title='Field-dependent part of the region offsets')
    axs[1, 1].legend(frameon=False, fontsize=8)
    for rid in range(8):
        g = T[(T.kind == 'drift') & (T.scheme == 'reg8') & (T.region == rid) & (T.dV >= 25)].sort_values('dV')
        c_ = dr.set_index('subrun').loc[g.subrun, 'tom_med_300'].to_numpy()
        axs[1, 2].errorbar(g.E, g.tom_med_300 - c_, g.tom_med_e_300, fmt='o-', ms=3, label=f'reg8 {rid}')
    axs[1, 2].axhline(0, color='k', lw=1)
    axs[1, 2].set(xlabel='drift field [V/cm]', ylabel='region - chamber med(tom) [ns]', title='reg8 offsets vs drift field (A>=300)')
    axs[1, 2].legend(frameon=False, fontsize=7, ncol=2)
    for a in axs.ravel():
        a.grid(color=GRID)
    fig.suptitle('det2 drift scan: region-to-region timing pattern vs field')
    fig.tight_layout(); fig.savefig(os.path.join(out, 'drift_region_offsets.png'), dpi=120, facecolor='white'); plt.close(fig)

    # ---------------- 5c amplitude collapse
    P('\n=== 5c  Amplitude vs dV per region: collapse A_r(dV) = g_r F(dV / s_r)  (F = cubic spline in dV, log amplitude)')
    COL = {}
    for scheme in ('reg8', 'reg24'):
        for col in ('amp_mpv', 'amp_med'):
            for dmin in (20, 50):
                f1 = collapse_fit(T, scheme, col, dVmin=dmin, s_free=True)
                f0 = collapse_fit(T, scheme, col, dVmin=dmin, s_free=False)
                dchi = f0['chi2'] - f1['chi2']; ddof = f0['ndf'] - f1['ndf']
                COL[(scheme, col, dmin)] = (f1, f0)
                P(f'   {scheme} {col} dV>={dmin}: s free chi2/ndf {f1["chi2"]:.1f}/{f1["ndf"]}; s=1 chi2/ndf {f0["chi2"]:.1f}/{f0["ndf"]}; '
                  f'Delta chi2 {dchi:.1f} for {ddof} dof (p {sst.chi2.sf(dchi, ddof):.2g}); '
                  f'rescaled by chi2/ndf of the free fit: p {sst.chi2.sf(dchi / max(1, f1["chi2"] / f1["ndf"]), ddof):.2g}; '
                  f's_r spread (rms) {np.std(f1["s"]):.3f}')
    for scheme in ('reg8', 'reg24'):
        f1 = COL[(scheme, 'amp_med', 20)][0]
        sc_ = max(1, np.sqrt(f1['chi2'] / f1['ndf']))
        P(f'   {scheme} s_r (amp_med, dV>=20; errors x {sc_:.2f} = sqrt(chi2/ndf) where > 1) and g_r:')
        P('     ' + '  '.join(f'{r}:{s:.3f}+-{e * sc_:.3f}/g{g:.2f}' for r, s, e, g in zip(f1['regs'], f1['s'], f1['s_e'], f1['g'])))
        f1m = COL[(scheme, 'amp_mpv', 20)][0]
        sm = max(1, np.sqrt(f1m['chi2'] / f1m['ndf']))
        P(f'   {scheme} s_r (amp_mpv, dV>=20; errors x {sm:.2f}):')
        P('     ' + '  '.join(f'{r}:{s:.3f}+-{e * sm:.3f}' for r, s, e in zip(f1m['regs'], f1m['s'], f1m['s_e'])))
        rr = np.corrcoef(f1['s'], f1m['s'])[0, 1]
        P(f'     Pearson r between s_r from median and from MPV: {rr:.2f}')
    for scheme in ('reg8', 'reg24'):
        dd = DR[DR.scheme == scheme].set_index('region').d_mm
        for col in ('amp_med', 'amp_mpv'):
            f1 = COL[(scheme, col, 20)][0]
            s = pd.Series(f1['s'], index=f1['regs'])
            P(f'   cross-check {scheme}: Pearson r(s_r from {col} collapse, d_r/4mm from high-field timing) = '
              f'{np.corrcoef(s.reindex(dd.index), dd / 4.0)[0, 1]:.2f}')
        # timing-shape cross-checks at E >= 625 V/cm
        h = T[(T.kind == 'drift') & (T.dV >= 250)]
        gm = h.groupby(['scheme', 'region'])[['tom_s68_300', 'dtm_med_300', 'amp_med']].mean()
        c = gm.loc['all'].iloc[0]; x = gm.loc[scheme].reindex(dd.index)
        P(f'   cross-check {scheme} (E>=625): Pearson r(d_r, sigma68 offset) {np.corrcoef(dd, x.tom_s68_300 - c.tom_s68_300)[0, 1]:.2f}, '
          f'r(d_r, rise time tom-time offset) {np.corrcoef(dd, x.dtm_med_300 - c.dtm_med_300)[0, 1]:.2f}, '
          f'r(d_r, amplitude ratio) {np.corrcoef(dd, x.amp_med / c.amp_med)[0, 1]:.2f}')
    rows = []
    for (scheme, col, dmin), (f1, f0) in COL.items():
        sc_ = max(1, np.sqrt(f1['chi2'] / f1['ndf']))
        for r, s, e, g in zip(f1['regs'], f1['s'], f1['s_e'], f1['g']):
            rows.append(dict(scheme=scheme, col=col, dVmin=dmin, region=r, s=s, s_e=e, s_e_scaled=e * sc_, g=g,
                             chi2=f1['chi2'], ndf=f1['ndf'], chi2_s1=f0['chi2'], ndf_s1=f0['ndf']))
    pd.DataFrame(rows).to_csv(os.path.join(out, 'amplitude_collapse.csv'), index=False)
    # figure 3
    fig, axs = plt.subplots(1, 3, figsize=(18, 5.2))
    f1 = COL[('reg8', 'amp_med', 20)][0]
    for i, rid in enumerate(f1['regs']):
        g = T[(T.kind == 'drift') & (T.scheme == 'reg8') & (T.region == rid)].sort_values('dV')
        l, = axs[0].plot(g.dV, g.amp_med, 'o-', ms=3, label=f'reg8 {rid}')
        k = g.dV >= 20
        axs[1].plot(g.dV[k] / f1['s'][i], g.amp_med[k] / f1['g'][i], 'o', ms=3, color=l.get_color())
    u = np.linspace(f1['knots'][0], f1['knots'][-1], 200)
    axs[1].plot(u, np.exp(CubicSpline(f1['knots'], f1['f'], bc_type='natural')(u)), 'k-', lw=1, label='F')
    axs[0].set(xlabel='Vdrift - Vmesh [V]', ylabel='median leading amplitude [ADC]', title='Leading-pad amplitude vs drift voltage, reg8', yscale='log')
    axs[1].set(xlabel='(Vdrift - Vmesh) / s_r [V]', ylabel='amplitude / g_r', title=f'Collapsed (chi2/ndf {f1["chi2"]:.0f}/{f1["ndf"]})', yscale='log')
    f24 = COL[('reg24', 'amp_med', 20)][0]
    pmap_(axs[2], pd.Series(f24['s'] - 1, index=f24['regs']), 'reg24 s_r - 1 (amp median collapse)')
    for a in axs:
        a.grid(color=GRID)
    axs[0].legend(frameon=False, fontsize=7, ncol=2); axs[1].legend(frameon=False, fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(out, 'amplitude_collapse.png'), dpi=120, facecolor='white'); plt.close(fig)

    # ---------------- 6 physical fit
    P('\n=== 6  Physical fit (drift scan, time_of_max, A>=300): s68_r = sqrt(s0^2 + (k T_r)^2), med_r = t0_r + c T_r,')
    P('       T_r = d_r / v_d((dV - V0)/d_r), v_d Magboltz (gas_table.csv + 50 V/cm point, linear to 0 below 50 V/cm)')
    PF = []
    variants = [('chamber', 'all', 50, {}), ('chamber V0 free', 'all', 50, dict(V0_free=True)),
                ('chamber dV>=25', 'all', 25, {}), ('chamber dV>=25 V0 free', 'all', 25, dict(V0_free=True)),
                ('chamber dV>=100', 'all', 100, {}), ('chamber dV>=200', 'all', 200, {})]
    for scheme in ('reg8', 'reg24'):
        variants += [(f'{scheme} d_r free', scheme, 50, {}), (f'{scheme} d common', scheme, 50, dict(d_common=True)),
                     (f'{scheme} d_r free V0 free', scheme, 50, dict(V0_free=True)),
                     (f'{scheme} d_r free dV>=25', scheme, 25, {}), (f'{scheme} d common dV>=25', scheme, 25, dict(d_common=True)),
                     (f'{scheme} d_r free dV>=100', scheme, 100, {}), (f'{scheme} d common dV>=100', scheme, 100, dict(d_common=True)),
                     (f'{scheme} d_r free dV>=200', scheme, 200, {}), (f'{scheme} d common dV>=200', scheme, 200, dict(d_common=True))]
    FITS = {}
    for name, scheme, dmin, kw in variants:
        pf = PhysFit(T, scheme, vd, cut=300, dVmin=dmin)
        o = pf.fit(**kw)
        FITS[name] = (o, pf)
        ie = {n: i for i, n in enumerate(o['names'])}
        dstr = (f"d_mean {o['d'].mean():.2f}+-{o.get('d_mean_e', o['err'][ie['d0']] if 'd0' in ie else np.nan):.2f} mm "
                f"d_r rms {o['d'].std():.3f} mm")
        P(f'   {name:28s}: chi2/ndf {o["chi2"]:7.1f}/{o["ndf"]:3d}  k {o["k"]:.3f}+-{o["err"][0]:.3f}  c {o["c"]:.3f}+-{o["err"][1]:.3f}  '
          f's0 {o["s0"]:.1f}+-{o["err"][2]:.1f} ns  V0 {o["V0"]:+.1f}' + (f'+-{o["err"][3]:.1f}' if 'V0' in ie else ' (fixed)') +
          f' V  {dstr}' + (f'  corr(k,d_mean) {o["corr_k_dmean"]:+.2f} corr(c,d_mean) {o["corr_c_dmean"]:+.2f}' if 'corr_k_dmean' in o else
                         (f'  corr(k,d) {o["corr"][0, ie["d0"]]:+.2f} corr(c,d) {o["corr"][1, ie["d0"]]:+.2f}' if 'd0' in ie else '')) +
          (f'  at bound: {[o["names"][i] for i in np.flatnonzero(o["active"])]}' if o['active'].any() else ''))
        PF.append(dict(name=name, scheme=scheme, dVmin=dmin, chi2=o['chi2'], ndf=o['ndf'], k=o['k'], k_e=o['err'][0],
                       c=o['c'], c_e=o['err'][1], s0=o['s0'], s0_e=o['err'][2], V0=o['V0'],
                       V0_e=o['err'][3] if 'V0' in ie else np.nan, d_mean=o['d'].mean(),
                       d_mean_e=o.get('d_mean_e', o['err'][ie['d0']] if 'd0' in ie else np.nan), d_rms=o['d'].std()))
    pd.DataFrame(PF).to_csv(os.path.join(out, 'physfit_variants.csv'), index=False)
    for scheme in ('reg8', 'reg24'):
        for dmin in ('', ' dV>=25', ' dV>=100', ' dV>=200'):
            of, _ = FITS[f'{scheme} d_r free{dmin}']; oc, _ = FITS[f'{scheme} d common{dmin}']
            dchi = oc['chi2'] - of['chi2']; ddof = oc['ndf'] - of['ndf']
            scl = max(1, of['chi2'] / of['ndf'])
            P(f'   {scheme}{dmin or " dV>=50"}: d_r differ? Delta chi2 (common - free) {dchi:.1f} for {ddof} dof, p {sst.chi2.sf(dchi, ddof):.2g}; '
              f'after dividing by chi2/ndf {scl:.2f}: p {sst.chi2.sf(dchi / scl, ddof):.2g}')
    for scheme in ('reg8', 'reg24'):
        o, pf = FITS[f'{scheme} d_r free']
        scl = np.sqrt(max(1, o['chi2'] / o['ndf']))
        P(f'   {scheme} d_r [mm] (dV>=50, V0=0; abs. error / error of d_r/d_mean, both x{scl:.2f}):')
        P('     ' + '  '.join(f'{r}:{d:.2f}+-{e * scl:.2f}({de * scl * 100:.1f}%)' for r, d, e, de in
                             zip(pf.regs, o['d'], o['d_e'], o['drel_e'])))
        o25, pf25 = FITS[f'{scheme} d_r free dV>=100']
        rr = np.corrcoef(o['d'], o25['d'])[0, 1]
        P(f'     Pearson r of d_r (dV>=50) with d_r (dV>=100): {rr:.2f}')
        pd.DataFrame(dict(region=pf.regs, d_mm=o['d'], d_e=o['d_e'], d_rel_e=o['drel_e'], t0=o['t0'],
                          d_mm_dV100=o25['d'], scale=scl)).to_csv(os.path.join(out, f'physfit_d_{scheme}.csv'), index=False)
    # figure 4: chamber fit + reg24 d_r map
    fig, axs = plt.subplots(1, 3, figsize=(18, 5.2))
    o, pf = FITS['chamber']
    dvg = np.linspace(20, 530, 300)
    Tm = pf.Tr(o['d'], o['V0'], ri=np.zeros(len(dvg), int), dV=dvg)
    axs[0].errorbar(pf.dV, pf.m, pf.me, fmt='o', color='#c0392b', label='median tom (data)')
    axs[0].plot(dvg, o['t0'][0] + o['c'] * Tm, 'k-', label=f't0 + c T, c={o["c"]:.2f}, d={o["d"][0]:.2f} mm')
    axs[0].set(xlabel='Vdrift - Vmesh [V]', ylabel='[ns]', title='Chamber: median time_of_max (A>=300)', ylim=(400, 600))
    axs[1].errorbar(pf.dV, pf.s, pf.se, fmt='o', color='#2471a3', label='sigma68 tom (data)')
    axs[1].plot(dvg, np.sqrt(o['s0'] ** 2 + (o['k'] * Tm) ** 2), 'k-', label=f'k={o["k"]:.2f}, s0={o["s0"]:.1f} ns')
    axs[1].set(xlabel='Vdrift - Vmesh [V]', ylabel='[ns]', title='Chamber: sigma68 time_of_max (A>=300)', ylim=(0, 120))
    o24, pf24 = FITS['reg24 d_r free']
    pmap_(axs[2], pd.Series(o24['d'], index=pf24.regs), 'reg24 fitted d_r [mm] (V0=0, dV>=50)', cmap='viridis', sym=False)
    for a in axs[:2]:
        a.grid(color=GRID); a.legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(out, 'physfit.png'), dpi=120, facecolor='white'); plt.close(fig)

    # ---------------- 7 mesh scan
    P('\n=== 7  Mesh scan (drift = mesh + 300 V, 750 V/cm nominal)')
    ms = ch[ch.kind == 'mesh'].sort_values('vmesh')
    P('   mesh V | n300 n150 | med tom(300) s68 | med tom(150) s68 | tom-time(300) | TOT(300) | int/amp(300) | MPV  median  sat  ev/trig')
    for r in ms.itertuples():
        P(f'   {r.vmesh:6.1f} | {r.n_300:5d} {r.n_150:5d} | {r.tom_med_300:6.1f}+-{r.tom_med_e_300:3.1f} {r.tom_s68_300:5.1f} | '
          f'{r.tom_med_150:6.1f}+-{r.tom_med_e_150:3.1f} {r.tom_s68_150:5.1f} | {r.dtm_med_300:6.1f}+-{r.dtm_med_e_300:3.1f} | '
          f'{r.tot_med_300:5.0f}+-{r.tot_med_e_300:3.0f} | {r.ioa_med_300:5.2f}+-{r.ioa_med_e_300:4.2f} | {r.amp_mpv:5.0f} {r.amp_med:5.0f} '
          f'{r.sat_frac:.3f} {r.ev_per_trig:.3f}')
    P('   Fixed amplitude band, chamber: slope vs mesh V [per 10 V] of med(tom), med(time), med(tom-time), med TOT, med int/amp')
    mb = Bt[(Bt.kind == 'mesh') & (Bt.scheme == 'all')]
    MB = []
    for (lo, hi), g in mb.groupby(['band_lo', 'band_hi']):
        g = g[g.n >= 30].sort_values('vmesh')
        line = f'     {lo:4d}-{hi:4d} ADC (n {g.n.min() if len(g) else 0}-{g.n.max() if len(g) else 0}, {len(g)} pts): '
        o = dict(band_lo=lo, band_hi=hi, npts=len(g))
        for q in ('tom', 'time', 'dtm', 'tot', 'ioa'):
            if len(g) < 3:
                continue
            a, ea, b, eb, chi2, ndf = wlin(g.vmesh.to_numpy(), g[f'{q}_med'].to_numpy(), g[f'{q}_med_e'].to_numpy())
            line += f'{q} {10 * b:+.2f}+-{10 * eb:.2f} (chi2 {chi2:.0f}/{ndf})  '
            o[f'{q}_slope10'] = 10 * b; o[f'{q}_slope10_e'] = 10 * eb; o[f'{q}_chi2'] = chi2; o[f'{q}_ndf'] = ndf
        MB.append(o); P(line)
    pd.DataFrame(MB).to_csv(os.path.join(out, 'mesh_band_slopes.csv'), index=False)
    # amplitude dependence at fixed voltage (time walk): pooled alignment run, chamber
    ab = Bt[(Bt.kind == 'alignment') & (Bt.scheme == 'all')]
    P('   Amplitude dependence at fixed voltage (alignment run 430/730, chamber): band: med tom / time / tom-time / TOT / int/amp')
    for r in ab.itertuples():
        P(f'     {r.band_lo:4d}-{r.band_hi:4d} (n {r.n:5d}, med amp {r.amp_med:5.0f}): {r.tom_med:6.1f}+-{r.tom_med_e:3.1f}  '
          f'{r.time_med:6.1f}+-{r.time_med_e:3.1f}  {r.dtm_med:6.1f}+-{r.dtm_med_e:3.1f}  {r.tot_med:5.0f}+-{r.tot_med_e:3.0f}  {r.ioa_med:5.2f}+-{r.ioa_med_e:4.2f}')
    # region dependence in mesh scan: offsets of med(tom) per reg8, constant across mesh V?
    P('   Mesh scan reg8 offsets of med(tom) from chamber (A>=150): constant across mesh V? (chi2 of constant / ndf)')
    m8 = T[(T.kind == 'mesh') & (T.scheme == 'reg8')].merge(ms[['subrun', 'tom_med_150']].rename(columns={'tom_med_150': 'chm'}), on='subrun')
    m8['dev'] = m8.tom_med_150 - m8.chm
    for rid, g in m8.groupby('region'):
        e = g.tom_med_e_150; w = 1 / e ** 2
        c0 = (g.dev * w).sum() / w.sum(); chi = (((g.dev - c0) / e) ** 2).sum()
        a, ea, b, eb, chi2, ndf = wlin(g.vmesh.to_numpy(), g.dev.to_numpy(), e.to_numpy())
        P(f'     reg8 {rid}: mean offset {c0:+5.1f}+-{1 / np.sqrt(w.sum()):.1f} ns, chi2 const {chi:.1f}/{len(g) - 1}, slope {10 * b:+.2f}+-{10 * eb:.2f} ns/10V')
    # compare to drift-scan reg8 offsets at 750 V/cm
    d8 = T[(T.kind == 'drift') & (T.scheme == 'reg8') & (T.dV >= 250)].merge(dr[['subrun', 'tom_med_150']].rename(columns={'tom_med_150': 'chm'}), on='subrun')
    d8['dev'] = d8.tom_med_150 - d8.chm
    a8 = d8.groupby('region').apply(lambda g: np.average(g.dev, weights=1 / g.tom_med_e_150 ** 2))
    b8 = m8.groupby('region').apply(lambda g: np.average(g.dev, weights=1 / g.tom_med_e_150 ** 2))
    P(f'     reg8 offsets (A>=150) drift scan E>=625 vs mesh scan: ' + ', '.join(f'{i}: {a8[i]:+.1f}/{b8[i]:+.1f}' for i in a8.index)
      + f'  Pearson r {np.corrcoef(a8, b8)[0, 1]:.2f}')
    # figure 5
    fig, axs = plt.subplots(2, 3, figsize=(17, 9.5))
    cols = plt.cm.viridis(np.linspace(0, 0.9, len(BANDS)))
    for (lo, hi), c in zip(BANDS, cols):
        g = mb[(mb.band_lo == lo) & (mb.n >= 30)].sort_values('vmesh')
        for ax, q, yl in [(axs[0, 0], 'tom', 'median time_of_max [ns]'), (axs[0, 1], 'time', 'median time (30% CFD) [ns]'),
                          (axs[0, 2], 'dtm', 'median time_of_max - time [ns]'), (axs[1, 0], 'tot', 'median TOT [ns]'),
                          (axs[1, 1], 'ioa', 'median integral / amplitude [samples]')]:
            ax.errorbar(g.vmesh, g[f'{q}_med'], g[f'{q}_med_e'], fmt='o-', color=c, ms=4, label=f'{lo}-{hi} ADC')
            ax.set(xlabel='mesh voltage [V]', ylabel=yl, title=f'{yl.split(" [")[0]} at fixed amplitude')
    axs[1, 2].errorbar(ms.vmesh, ms.amp_mpv, ms.amp_mpv_e, fmt='o-', color='#c0392b', ms=4, label='MPV')
    axs[1, 2].errorbar(ms.vmesh, ms.amp_med, ms.amp_med_e, fmt='s--', color='#2471a3', ms=4, label='median')
    axs[1, 2].set(xlabel='mesh voltage [V]', ylabel='leading amplitude [ADC]', yscale='log', title='Leading-pad amplitude vs mesh voltage')
    for a in axs.ravel():
        a.grid(color=GRID); a.legend(frameon=False, fontsize=7)
    fig.suptitle('det2 mesh scan 24 Nov 2026 (drift = mesh + 300 V): pulse shape/timing at fixed amplitude')
    fig.tight_layout(); fig.savefig(os.path.join(out, 'mesh_scan.png'), dpi=120, facecolor='white'); plt.close(fig)

    with open(os.path.join(out, 'results_log.txt'), 'w') as f:
        f.write('\n'.join(S) + '\n')
    print('->', out)


if __name__ == '__main__':
    main()
