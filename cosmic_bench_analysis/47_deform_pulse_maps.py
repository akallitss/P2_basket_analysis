#!/usr/bin/env python3
"""47_deform_pulse_maps.py -- drift time per pad from the pulse shapes (stage-43 waveforms).

IDEA. A muon ionises uniformly along the drift gap, so on one pad the mean pulse is

    M(t) = box(0, T_d)  (x)  SER(t - t0),      T_d = gap / v_d

(SER = ion tail (x) shaper, t0 = electronics delay). For two places with the SAME
SER, box(T_ref) (x) M_pad = box(T_pad) (x) M_ref: the pad pulse is the reference
pulse with a longer or shorter box, and a pure delay shows up as a separate shift
dt0. Three estimators of T_d per pad, with different systematics:

  T_fit   cross-convolution fit of the pad's mean pulse to the reference mean
          pulse (free: T_pad, dt0, scale; T_ref fixed). Uses the whole shape;
          a change of the ion transit time would leak into it.
  T_k2    from the variance (2nd cumulant) of the mean pulse only:
          T^2 = T_ref^2 + 12 (k2_pad - k2_ref).  Blind to delays.
  T_s     from the EVENT-TO-EVENT spread of the leading-edge time,
          T = T_ref * s68_pad / s68_ref (primary-cluster statistics; the ion
          tail and the electronics are deterministic and cannot change it).
          Computed on SINGLE-PAD events only (no other pad above 100 ADC
          within 25 mm of the leading pad):
          when an inclined track shares its charge with a neighbour, the pad
          sees only part of the gap and its time shifts (the micro-TPC effect
          of stage 48), which inflates the spread, and more so in a thicker gap.
          On single-pad events the spread is proportional to the drift time
          (gap classes: ratios 0.93 / 1 / 1.16 / 1.48 for T ratios 0.89 / 1 /
          1.14 / 1.50).

The reference drift time T_ref is an argument (default: 4 mm / v_d of dry
Magboltz Ar/iC4H10 95/5 at the run voltages). For det2_long_repos 125 ns was
used: the value the free band-vs-reference fit returns; Magboltz gives 128 ns
and the micro-TPC of stage 48 measures 130 +- 3 ns, independently of the pulse
shape.

Sample: the leading (max amplitude) hit pad of every good event, 300-3500 ADC,
not saturated, peak sample 4-12. No track required (statistics); the track-matched
amplitude map comes from the tracked subsample (pad under the track, track within
4 mm of the pad centre, path-length corrected).

Out (<run>/47_deform_pulse_maps/): pads.csv (per pad and sliding R = 18 mm), pad_gap_prelim.csv,
regions.csv (reg24/reg48), mean_pulse_ref.csv, maps_*.png, summary.txt
Usage: python3 47_deform_pulse_maps.py [run_key] [--t-ref 125] [--ref-conn 4 5 6 7]
       (det2_long_repos was run with --t-ref 125)
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2_deform as pdf          # noqa: E402
import p2_deform_model as dm     # noqa: E402

TG = np.arange(100.0, 1300.0, 10.0)      # absolute-time grid of the mean pulses [ns]
FIT = (TG > 150) & (TG < 1100)
RNG = np.random.default_rng(20261128)


def mean_pulse(WF, rows, ftst, amp):
    """Sum of raw waveforms / sum of amplitudes on a 10 ns absolute-time grid."""
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
    """Convolution of the gridded pulse m with a unit-area box of length T [ns]."""
    c = np.concatenate([[0.0], np.cumsum(m) * dt])
    te = np.concatenate([[TG[0] - dt / 2], TG + dt / 2])
    return (np.interp(TG, te, c) - np.interp(TG - T, te, c, left=0.0)) / T


def xfit(m_pad, m_ref, t_ref, free_dt=True):
    """box(t_ref) (x) m_pad(t - dt0) * A = box(T) (x) m_ref  ->  (T, dt0, A, rms)."""
    rhs_cache = {}

    def resid(p):
        T, dt0, A = p
        if not free_dt:
            dt0 = 0.0
        lhs = boxconv(np.interp(TG - dt0, TG, m_pad, left=0.0, right=0.0), t_ref) * A
        return (lhs - boxconv(m_ref, T))[FIT]
    best = None
    for T0 in (t_ref * 0.8, t_ref * 1.5):
        r = least_squares(resid, [T0, 0.0, 1.0], bounds=([15, -80, 0.3], [600, 80, 3]))
        if best is None or r.cost < best.cost:
            best = r
    return best.x[0], (best.x[1] if free_dt else 0.0), best.x[2], float(np.sqrt(2 * best.cost / FIT.sum()))


def cumulants(m, lo=120.0, hi=1000.0):
    k = (TG >= lo) & (TG <= hi)
    a = m[k].sum()
    k1 = (TG[k] * m[k]).sum() / a
    k2 = ((TG[k] - k1) ** 2 * m[k]).sum() / a
    return k1, k2


def measure(s, WF, m_ref, k_ref, s_ref, t_ref, nboot=0):
    """All estimators for one group of pulses (DataFrame s with 'row')."""
    rows = s.row.to_numpy(); ft = s.ftst.to_numpy(); amp = s.amp.to_numpy()
    m = mean_pulse(WF, rows, ft, amp)
    T, dt0, A, rms = xfit(m, m_ref, t_ref)
    T0, _, _, rms0 = xfit(m, m_ref, t_ref, free_dt=False)
    k1, k2 = cumulants(m)
    T2sq = t_ref ** 2 + 12 * (k2 - k_ref[1])
    out = dict(n=len(s), T_fit=T, dt0=dt0, A=A, rms=rms, T_fit0=T0, rms0=rms0,
               dk1=k1 - k_ref[0], dk2=k2 - k_ref[1], T_k2=np.sign(T2sq) * np.sqrt(abs(T2sq)),
               amp_med=float(np.median(amp)))
    s1 = s[s.nh == 1]                       # spreads: single-pad events only (see docstring)
    out.update(n1=len(s1), s_l20=pdf.s68(s1.t_l20), s_l50=pdf.s68(s1.t_l50), s_pk=pdf.s68(s1.t_pk))
    out['T_s'] = t_ref * out['s_l20'] / s_ref
    for c in ('t_l10', 't_l20', 't_pk', 't_t50', 'rise', 'r2050', 'fwhm', 'fall', 'wid'):
        out[f'med_{c}'] = float(np.nanmedian(s[c]))
    if nboot:
        b = []
        for _ in range(nboot):
            i = RNG.integers(0, len(s), len(s))
            mb = mean_pulse(WF, rows[i], ft[i], amp[i])
            Tb, db, _, _ = xfit(mb, m_ref, t_ref)
            _, k2b = cumulants(mb)
            j = RNG.integers(0, len(s1), len(s1))
            b.append((Tb, db, pdf.s68(s1.t_l20.to_numpy()[j]), k2b))
        b = np.array(b)
        out['e_T_fit'], out['e_dt0'] = b[:, 0].std(), b[:, 1].std()
        out['e_T_s'] = t_ref * b[:, 2].std() / s_ref
        out['e_dk2'] = b[:, 3].std()
    return out


def smap(ax, P, c, title, cmap='coolwarm', lo=None, hi=None, sym=None):
    v = P[c].to_numpy(float)
    if lo is None:
        lo, hi = np.nanpercentile(v, [2, 98])
    if sym is not None:
        a = max(abs(lo - sym), abs(hi - sym)); lo, hi = sym - a, sym + a
    sc = ax.scatter(P.pad_cx, P.pad_cy, c=v, s=20, marker='s', cmap=cmap, vmin=lo, vmax=hi)
    ax.set_aspect('equal'); ax.set_title(title, fontsize=11)
    ax.set_xlabel('pad x [mm]'); ax.set_ylabel('pad y [mm]')
    plt.colorbar(sc, ax=ax, fraction=0.04)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run_key', nargs='?', default='det2_long_repos')
    ap.add_argument('--t-ref', type=float, default=None,
                    help='drift time of the reference zone [ns]; default 4 mm / v_d(dry Magboltz) at the run voltages')
    ap.add_argument('--ref-conn', type=int, nargs='+', default=[4, 5, 6, 7])
    ap.add_argument('--slide-r', type=float, default=18.0)
    ap.add_argument('--nboot', type=int, default=40)
    ap.add_argument('--amp-min', type=float, default=300.0,
                    help='lower amplitude cut of the timing sample [ADC]; a non-default value writes to 47_deform_pulse_maps_amp<N>')
    a = ap.parse_args()
    R, WF, cfg = pdf.load(a.run_key)
    out = cfg.out_dir('47_deform_pulse_maps' + ('' if a.amp_min == 300.0 else f'_amp{a.amp_min:g}'))
    nums = [int(x) for x in cfg.SUB_RUN.split('_') if x.isdigit()]
    dV = float(abs(nums[-1] - nums[-2]))
    if a.t_ref is None:
        a.t_ref = float(dm.t_drift(4.0, dV))
    R['row'] = np.arange(len(R))
    h = R[R.is_hit & R.good & ~R.saturated]
    h = h.sort_values('amp', ascending=False).drop_duplicates('eventId')
    # nh = pads above 100 ADC within 25 mm of the leading pad (itself included): 1 = single-pad event.
    # Local, so that a noisy connector elsewhere does not turn every event into a multi-pad one.
    o = R.loc[R.is_hit & R.good & (R.amp > 100), ['eventId', 'pad_cx', 'pad_cy']].merge(
        h[['eventId', 'pad_cx', 'pad_cy']], on='eventId', suffixes=('', '_L'))
    o = o[np.hypot(o.pad_cx - o.pad_cx_L, o.pad_cy - o.pad_cy_L) < 25.0]
    h = h.join(o.groupby('eventId').size().rename('nh'), on='eventId')
    h['nh'] = h.nh.fillna(1)
    s = h[(h.amp > a.amp_min) & (h.amp < 3500) & (h.imax >= 4) & (h.imax <= 12)].copy()
    L = [f'{a.run_key}: {len(s):,} leading pulses ({a.amp_min:g}-3500 ADC, not saturated, peak sample 4-12) on '
         f'{s.channel_id.nunique()} pads; reference = connectors {a.ref_conn}, T_ref = {a.t_ref:g} ns (argument)']

    ref = s[s.connector_N.isin(a.ref_conn)]
    m_ref = mean_pulse(WF, ref.row.to_numpy(), ref.ftst.to_numpy(), ref.amp.to_numpy())
    k_ref = cumulants(m_ref)
    # reference spread: median of the per-pad s68 (not the pooled one, which contains pad-to-pad offsets)
    r1 = ref[ref.nh == 1]
    s_ref = float(pdf.s68(r1.t_l20 - r1.groupby('channel_id').t_l20.transform('median')))
    pd.DataFrame({'t_ns': TG, 'm_ref': m_ref}).to_csv(f'{out}/mean_pulse_ref.csv', index=False)
    L.append(f'reference: {len(ref):,} pulses; per-pad s68(t_l20) median {s_ref:.2f} ns -> s68/T_ref = {s_ref/a.t_ref:.4f}; '
             f'mean-pulse centroid {k_ref[0]:.1f} ns, variance {k_ref[1]:.0f} ns^2')

    # ---- per pad (own pulses) and sliding (pulses of all pads within slide_r) ----
    reg = pd.read_csv(cfg.out_dir('43_deformation') + '/regions.csv')
    pads = reg.set_index('channel_id')
    by_pad = {int(c): g for c, g in s.groupby('channel_id')}
    tree = cKDTree(reg[['pad_cx', 'pad_cy']].to_numpy())
    rows_p, rows_s = [], []
    for i, pr in enumerate(reg.itertuples()):
        g = by_pad.get(int(pr.channel_id))
        if g is not None and len(g) >= 25:
            rows_p.append(dict(channel_id=pr.channel_id, **measure(g, WF, m_ref, k_ref, s_ref, a.t_ref)))
        nb = reg.channel_id.to_numpy()[tree.query_ball_point([pr.pad_cx, pr.pad_cy], a.slide_r)]
        gs = [by_pad[int(c)] for c in nb if int(c) in by_pad]
        if gs:
            gg = pd.concat(gs)
            if len(gg) >= 100:
                # the spread must not contain pad-to-pad offsets: remove each pad's median first
                gg = gg.copy()
                for c in ('t_l20', 't_l50', 't_pk'):
                    gg[c] = gg[c] - gg.groupby('channel_id')[c].transform('median')
                rows_s.append(dict(channel_id=pr.channel_id, npad=len(gs),
                                   **measure(gg, WF, m_ref, k_ref, s_ref, a.t_ref)))
    P = reg.merge(pd.DataFrame(rows_p), on='channel_id', how='left')
    S = pd.DataFrame(rows_s).add_prefix('sl_').rename(columns={'sl_channel_id': 'channel_id'})
    P = P.merge(S, on='channel_id', how='left')

    # track-matched amplitude (pad under the track, central 4 mm, cos(theta) corrected)
    tk = R[R.good & R.near_track & (R.d_trk < 4.0) & (R.sd >= 2)]
    tk = tk.assign(ampc=tk.amp * np.cos(np.radians(tk.theta)))
    A = tk[tk.is_hit].groupby('channel_id').agg(n_trk=('ampc', 'size'), ampc_med=('ampc', 'median'),
                                                ln_s68=('ampc', lambda v: pdf.s68(np.log(v))))
    A['hit_frac'] = tk.groupby('channel_id').is_hit.mean()
    P = P.merge(A, on='channel_id', how='left')
    P.to_csv(f'{out}/pads.csv', index=False)
    # preliminary gap per pad (dry-gas table), used by stages 45 and 48 to class pads
    from scipy.optimize import brentq
    def _gap(T):
        if not np.isfinite(T):
            return np.nan
        try:
            return brentq(lambda d: float(dm.t_drift(d, dV)) - T, 1.0, 14.0)
        except ValueError:
            return np.nan
    G = P[['channel_id', 'pad_cx', 'pad_cy', 'connector_N', 'sl_T_fit', 'sl_T_s']].copy()
    G['T_comb'] = G[['sl_T_fit', 'sl_T_s']].mean(axis=1)
    G['gap_prelim_mm'] = [_gap(t) for t in G.T_comb]
    G['gap_class'] = pd.cut(G.gap_prelim_mm, [0, 3.5, 3.9, 4.3, 5.0, 6.0, 20], labels=False)
    G.to_csv(f'{out}/pad_gap_prelim.csv', index=False)

    # ---- regions, with bootstrap errors ----
    rr = []
    for K in (24, 48):
        for r, g in s.groupby(f'reg{K}'):
            g = g.copy()
            for c in ('t_l20', 't_l50', 't_pk'):
                g[c] = g[c] - g.groupby('channel_id')[c].transform('median') + g[c].median()
            o = measure(g, WF, m_ref, k_ref, s_ref, a.t_ref, nboot=a.nboot)
            t_ = tk[tk[f'reg{K}'] == r]
            hh = t_[t_.is_hit]
            o.update(K=K, region=int(r), x=g.pad_cx.mean(), y=g.pad_cy.mean(), n_trk=len(t_),
                     hit_frac=t_.is_hit.mean(), ampc_med=hh.ampc.median(), ampc_mpv=dm.mpv(hh.ampc.to_numpy()),
                     ln_s68=pdf.s68(np.log(hh.ampc)), frac_lt300=(hh.ampc < 300).mean())
            rr.append(o)
    Rg = pd.DataFrame(rr)
    Rg.to_csv(f'{out}/regions.csv', index=False)

    # ---- summary numbers ----
    ok = P.T_fit.notna()
    L.append(f'per pad ({ok.sum()} pads with >= 25 pulses, median {P.n.median():.0f}): '
             f'T_fit median {P.T_fit.median():.1f}, 5-95 % {P.T_fit.quantile(.05):.0f}-{P.T_fit.quantile(.95):.0f} ns; '
             f'dt0 median {P.dt0.median():+.1f}, robust sd {pdf.s68(P.dt0):.1f} ns')
    q = P[ok]
    L.append('correlations between the three estimators (per pad / sliding):')
    for x_, y_ in (('T_fit', 'T_s'), ('T_fit', 'T_k2'), ('T_s', 'T_k2')):
        L.append(f'   {x_} vs {y_}: r = {q[[x_, y_]].corr().iloc[0, 1]:+.2f} / '
                 f'{P[["sl_" + x_, "sl_" + y_]].corr().iloc[0, 1]:+.2f}')
    L.append('regions (reg24): region  n   T_fit +- err   dt0 +- err   T_s +- err   T_k2   rms[1e-3]  amp(track)  ln-width  hit frac')
    for r in Rg[Rg.K == 24].sort_values('T_fit').itertuples():
        L.append(f'   {r.region:2d} {r.n:6d}  {r.T_fit:6.1f} +- {r.e_T_fit:4.1f}  {r.dt0:+5.1f} +- {r.e_dt0:3.1f}  '
                 f'{r.T_s:6.1f} +- {r.e_T_s:4.1f}  {r.T_k2:6.1f}  {1e3*r.rms:5.2f}  {r.ampc_med:6.0f}  {r.ln_s68:.3f}  {r.hit_frac:.3f}')
    L.append('assumptions: T_ref is an input (see stage 48 for its fit); sliding values pool the pads within '
             f'{a.slide_r:g} mm; spreads are computed after removing each pad\'s median; no track requirement for the '
             'timing sample; alignment only enters the track-matched amplitude (stage-34 trimmed fit).')

    # ---- maps ----
    fig, ax = plt.subplots(3, 4, figsize=(25, 17))
    smap(ax[0, 0], P, 'sl_T_fit', 'drift time T_d from the mean-pulse fit [ns] (sliding 18 mm)', 'viridis')
    smap(ax[0, 1], P, 'sl_T_s', 'T_d from the event-to-event spread [ns] (sliding)', 'viridis',
         *np.nanpercentile(P.sl_T_fit, [2, 98]))
    smap(ax[0, 2], P, 'sl_T_k2', 'T_d from the mean-pulse variance [ns] (sliding)', 'viridis',
         *np.nanpercentile(P.sl_T_fit, [2, 98]))
    smap(ax[0, 3], P, 'sl_dt0', 'fitted extra delay dt0 [ns] (sliding)', sym=0.0)
    smap(ax[1, 0], P, 'T_fit', 'T_d mean-pulse fit, single pads [ns]', 'viridis', *np.nanpercentile(P.sl_T_fit, [2, 98]))
    smap(ax[1, 1], P, 'T_s', 'T_d spread, single pads [ns]', 'viridis', *np.nanpercentile(P.sl_T_fit, [2, 98]))
    smap(ax[1, 2], P, 'sl_rms', 'fit residual rms (sliding)', 'magma')
    smap(ax[1, 3], P, 'ampc_med', 'track-matched median amplitude x cos(theta) [ADC]')
    smap(ax[2, 0], P, 'sl_med_rise', 'median rise (peak - 30 %) [ns] (sliding)')
    smap(ax[2, 1], P, 'sl_med_fwhm', 'median FWHM [ns] (sliding)')
    smap(ax[2, 2], P, 'sl_med_fall', 'median fall (50 % trailing - peak) [ns] (sliding)')
    smap(ax[2, 3], P, 'connector_N', 'connector', 'tab10', 1.5, 9.5)
    fig.suptitle(f'{cfg.DET_NAME} {a.run_key}: drift time per pad from the pulse shape (T_ref = {a.t_ref:g} ns)')
    fig.tight_layout(); fig.savefig(f'{out}/maps_drift_time.png', dpi=110); plt.close(fig)

    with open(f'{out}/summary.txt', 'w') as f:
        f.write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
