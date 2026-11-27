#!/usr/bin/env python3
"""48_deform_microtpc.py -- drift time from pad PAIRS under inclined tracks (micro-TPC).

An inclined track that crosses the boundary between two pads leaves the lower
part of its path (near the mesh) on one pad and the upper part (near the window)
on the other. The upper charge drifts longer, so the two pulses are separated in
time by about half the full drift time, with a sign given by the track direction
(M3). No pulse-shape model is involved: only the time difference of two pulses
of the same event, antisymmetrised in the slope, so channel delays cancel.

  s      track slope along the direction from pad A (larger pulse) to pad B
         (the adjacent pad with the next largest pulse), pad frame
  dt     t(B) - t(A)
  D(|s|) = [ median dt (s > +|s|) - median dt (s < -|s|) ] / 2

The sign of D tells which side the window is on (D > 0: the charge at larger M3
z arrives later, i.e. the window is at the top). The size of D at large |s| is
T_d / 2 reduced by transverse diffusion and by the pulse thresholds; that
reduction is computed with the forward model (p2_deform_model.simulate_pair)
with the SAME selection, for the nominal gap and for the two ways a place can
have a longer drift time (thicker gap, or slower drift in the same gap).

Out (<run>/48_deform_microtpc/): pairs.csv, curve_by_zone.csv, toy_curves.csv,
microtpc.png, summary.txt
Usage: python3 48_deform_microtpc.py [run_key] [--dt-scale=0.85]

The dilution factor depends on the transverse diffusion, and that on the gas: with 0.5-0.75 % water it is
0.83-0.88 of the dry-gas value (Magboltz), which lowers the fitted drift times by ~12 % (det2: nominal zone
129 -> 113 ns). The S-curve shape itself prefers the smaller diffusion.
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2_deform as pdf          # noqa: E402
import p2_deform_model as dm     # noqa: E402

SB = np.array([0.0, 0.05, 0.10, 0.15, 0.20, 0.26, 0.34, 0.50])
AMP_B = 120.0
EST = ('t_pk', 't_l20', 't_l50', 't_cen', 't_t50')
RNG = np.random.default_rng(20261128)


def pairs_from_features(FA, FB, extra=None):
    """Order two pads by amplitude -> (A = larger, B = smaller) feature frames."""
    sw = FB.amp.to_numpy() > FA.amp.to_numpy()
    A = FA.copy(); B = FB.copy()
    A[sw] = FB[sw].to_numpy(); B[sw] = FA[sw].to_numpy()
    return A, B, sw


def curve(p, est='t_pk', nboot=200):
    """Antisymmetrised D(|s|) with bootstrap errors."""
    rows = []
    for lo, hi in zip(SB[:-1], SB[1:]):
        a = p[(p.s > lo) & (p.s <= hi)]['d_' + est].to_numpy()
        b = p[(p.s < -lo) & (p.s >= -hi)]['d_' + est].to_numpy()
        if len(a) < 6 or len(b) < 6:
            rows.append((0.5 * (lo + hi), len(a), len(b), np.nan, np.nan)); continue
        D = 0.5 * (np.median(a) - np.median(b))
        bs = [0.5 * (np.median(a[RNG.integers(0, len(a), len(a))]) - np.median(b[RNG.integers(0, len(b), len(b))]))
              for _ in range(nboot)]
        rows.append((0.5 * (lo + hi), len(a), len(b), D, np.std(bs)))
    return pd.DataFrame(rows, columns=['s', 'n_pos', 'n_neg', 'D', 'e'])


def toy_curve(d_mm, v_scale, ser_par, amp_scale, n=6000, dV=300.0, dt_scale=1.0):
    """Model D(|s|) with the data selection (both pads hit, B > AMP_B, A 300-3500)."""
    H = dm.Heed()
    rows = []
    for lo, hi in zip(SB[:-1], SB[1:]):
        sl = 0.5 * (lo + hi)
        WA, WB, ph, info = dm.simulate_pair(n, sl, d_mm=d_mm, dV=dV, v_scale=v_scale, adc_per_e=amp_scale, noise=12.0, dt_scale=dt_scale,
                                            heed=H, rng=RNG, **ser_par)
        FA = pdf.features(WA); FB = pdf.features(WB)
        A, B, sw = pairs_from_features(FA, FB)
        ok = (B.amp > AMP_B).to_numpy() & (A.amp > AMP_B).to_numpy() & (A.amp < 3800).to_numpy()
        # dt = t(B) - t(A); geometric sign: the generated slope is positive from pad 'A_geo' to 'B_geo',
        # so after reordering by amplitude the slope along A->B is -sl where swapped
        sgn = np.where(sw, -1.0, 1.0)
        out = dict(s=sl, n=int(ok.sum()), frac_pair=float(ok.mean()), T_d=info['T_d'])
        for c in EST:
            dt = (B[c] - A[c]).to_numpy()[ok]; ss = sgn[ok]
            pos, neg = dt[ss > 0], dt[ss < 0]
            out['D_' + c] = 0.5 * (np.median(pos) - np.median(neg)) if len(pos) > 5 and len(neg) > 5 else np.nan
        rows.append(out)
    return pd.DataFrame(rows)


def main():
    args = [x for x in sys.argv[1:] if not x.startswith('--')]
    run_key = args[0] if args else 'det2_long_repos'
    # transverse diffusion of the model relative to dry Magboltz: 1.0 (dry) or ~0.85 (0.5-0.75 % water);
    # a non-default value writes to 48_deform_microtpc_dt<scale>
    dt_scale = float([x.split('=')[1] for x in sys.argv if x.startswith('--dt-scale=')][0]) if any(
        x.startswith('--dt-scale=') for x in sys.argv) else 1.0
    R, WF, cfg = pdf.load(run_key)
    out = cfg.out_dir('48_deform_microtpc' + ('' if dt_scale == 1.0 else f'_dt{dt_scale:g}'))
    L = []
    t = R[R.good & R.near_track & R.is_hit & ~R.saturated & (R.amp > AMP_B) & (R.imax >= 4) & (R.imax <= 12)].copy()
    t = t.sort_values(['eventId', 'amp'], ascending=[True, False])
    t['rk'] = t.groupby('eventId').cumcount()
    a = t[t.rk == 0].set_index('eventId'); b = t[t.rk == 1].set_index('eventId')
    p = a.join(b, how='inner', lsuffix='_a', rsuffix='_b')
    p['dist'] = np.hypot(p.pad_cx_a - p.pad_cx_b, p.pad_cy_a - p.pad_cy_b)
    p = p[p.dist < 14.5].copy()
    nx = (p.pad_cx_b - p.pad_cx_a) / p.dist; ny = (p.pad_cy_b - p.pad_cy_a) / p.dist
    p['s'] = p.tx_p_a * nx + p.ty_p_a * ny
    p['u'] = (p.xp_a - (p.pad_cx_a + p.pad_cx_b) / 2) * nx + (p.yp_a - (p.pad_cy_a + p.pad_cy_b) / 2) * ny
    for c in EST:
        p['d_' + c] = p[c + '_b'] - p[c + '_a']
    p = p[p.u.abs() < 3.0]
    gp = pd.read_csv(cfg.out_dir('47_deform_pulse_maps') + '/pad_gap_prelim.csv').set_index('channel_id')
    p['gap_a'] = gp.gap_prelim_mm.reindex(p.channel_id_a.to_numpy()).to_numpy()
    p['gap_b'] = gp.gap_prelim_mm.reindex(p.channel_id_b.to_numpy()).to_numpy()
    p['gap'] = 0.5 * (p.gap_a + p.gap_b)
    keep = ['channel_id_a', 'channel_id_b', 's', 'u', 'amp_a', 'amp_b', 'gap', 'theta_a'] + ['d_' + c for c in EST]
    p[keep].to_csv(f'{out}/pairs.csv')
    L.append(f'{run_key}: {len(p):,} adjacent pad pairs under a track (both pulses > {AMP_B:g} ADC, not saturated, '
             f'track within 3 mm of the shared edge)')

    zones = {'gap < 3.6 mm': p.gap < 3.6, 'gap 3.6-4.4 mm': (p.gap >= 3.6) & (p.gap < 4.4),
             'gap 4.4-5.5 mm': (p.gap >= 4.4) & (p.gap < 5.5), 'gap > 5.5 mm': p.gap >= 5.5, 'all': p.gap > 0}
    cur = []
    for z, m in zones.items():
        for est in EST:
            c = curve(p[m], est); c['zone'] = z; c['est'] = est; c['gap_mean'] = p[m].gap.mean(); cur.append(c)
    C = pd.concat(cur, ignore_index=True)
    C.to_csv(f'{out}/curve_by_zone.csv', index=False)

    # ---- forward model with the same selection --------------------------------------
    # SER parameters: fit of the model to the reference mean pulse (stage 47 exploration)
    ser_par = dict(T_ion=278.0, t_peak=174.0, f_e=0.6)
    # amplitude scale: median single-pad amplitude of a vertical muon = data (~1150 ADC)
    nums = [int(x) for x in cfg.SUB_RUN.split('_') if x.isdigit()]
    dV = float(abs(nums[-1] - nums[-2]))
    # amplitude scale of the model = the data: median amplitude x cos(theta) of the pad under a central track
    tkc = R[R.good & R.near_track & R.is_hit & (R.d_trk < 4.0)]
    amp_data = float(np.median(tkc.amp * np.cos(np.radians(tkc.theta))))
    W0, _, _, _ = dm.simulate(3000, d_mm=4.0, dV=dV, adc_per_e=1.0, noise=0.0, rng=RNG, **ser_par)
    a0 = W0.max(1)
    # the data median is over hits above the 5 sigma threshold (~60 ADC): match on the same truncated sample
    amp_scale = amp_data / np.median(a0)
    for _ in range(4):
        amp_scale = amp_data / np.median(a0[a0 * amp_scale > 60.0]) 
    toys = {'gap 4.0 mm (nominal)': (4.0, 1.0), 'gap 3.3 mm': (3.3, 1.0), 'gap 5.0 mm': (5.0, 1.0),
            'gap 7.0 mm': (7.0, 1.0), 'gap 4.0 mm, drift 1.5x slower': (4.0, 1 / 1.5)}
    TC = []
    for name, (d, vs) in toys.items():
        c = toy_curve(d, vs, ser_par, amp_scale, dV=dV, dt_scale=dt_scale); c['model'] = name; c['d_mm'] = d; c['v_scale'] = vs; TC.append(c)
    TC = pd.concat(TC, ignore_index=True)
    TC.to_csv(f'{out}/toy_curves.csv', index=False)

    L.append('sign: D > 0 at every slope -> the charge at larger M3 z arrives later -> the window (cathode) is on '
             'the upper side, the mesh/board below (M3 z up).')
    L.append('forward model, t_pk estimator: D / (T_d / 2) at |s| = 0.30 and 0.42, pair fraction at 0.30')
    for name, g in TC.groupby('model', sort=False):
        g = g.set_index('s')
        s3, s4 = g.index[-2], g.index[-1]
        L.append(f'   {name:32s} T_d {g.T_d.iloc[0]:6.1f} ns  D/(T_d/2) = {2*g.D_t_pk[s3]/g.T_d.iloc[0]:.3f}, '
                 f'{2*g.D_t_pk[s4]/g.T_d.iloc[0]:.3f}   D = {g.D_t_pk[s3]:.1f}, {g.D_t_pk[s4]:.1f} ns   pairs {100*g.frac_pair[s3]:.1f} %')

    # ---- T_d per zone: scale the model curve of the closest gap to the data ----------
    L.append('drift time per zone: least-squares scale of the model curve D(|s|) (thick-gap model of the nearest gap) '
             'to the data, |s| > 0.10, for three time estimators')
    res = []
    for z in zones:
        gm = C[C.zone == z].gap_mean.iloc[0]
        name = min(toys, key=lambda k: abs(toys[k][0] - gm) + (100 if toys[k][1] != 1 else 0))
        tm = TC[TC.model == name].set_index('s')
        line = f'   {z:16s} (mean gap estimate {gm:.2f} mm, model "{name}")'
        for est in ('t_pk', 't_l20', 't_cen'):
            dd = C[(C.zone == z) & (C.est == est)].set_index('s')
            k = dd.D.notna() & (dd.index > 0.10)
            y, e, m_ = dd.D[k].to_numpy(), dd.e[k].to_numpy(), tm['D_' + est].reindex(dd.index[k]).to_numpy()
            w = 1 / e ** 2
            sc = (w * y * m_).sum() / (w * m_ * m_).sum(); esc = 1 / np.sqrt((w * m_ * m_).sum())
            chi = (w * (y - sc * m_) ** 2).sum()
            Td = sc * tm.T_d.iloc[0]
            res.append(dict(zone=z, est=est, model=name, T_d=Td, e_T_d=esc * tm.T_d.iloc[0], chi2=chi, ndf=int(k.sum()) - 1,
                            n_pairs=int(dd.n_pos.sum() + dd.n_neg.sum())))
            line += f'  {est}: T_d = {Td:5.1f} +- {esc*tm.T_d.iloc[0]:4.1f} ns (chi2 {chi:.1f}/{int(k.sum())-1})'
        L.append(line)
    pd.DataFrame(res).to_csv(f'{out}/drift_time_by_zone.csv', index=False)
    L.append('assumptions: model SER = fit to the reference mean pulse (T_ion 278, peaking 174 ns, f_e 0.6; only the '
             f'threshold/dilution factor depends on it), Magboltz dry Ar/iC4H10 95/5 for v_d, its transverse diffusion x {dt_scale:g}, thresholds '
             f'as in data (both pads > {AMP_B:g} ADC), amplitude scale set to the data median of central tracks ({amp_data:.0f} ADC); zones from '
             'the stage-47 preliminary per-pad gap (mean of the two pads).')

    fig, axs = plt.subplots(1, 3, figsize=(20, 5.6))
    cols = {'gap < 3.6 mm': '#257de0', 'gap 3.6-4.4 mm': '#1baf7a', 'gap 4.4-5.5 mm': '#f2a600', 'gap > 5.5 mm': '#ec4647'}
    ax = axs[0]
    pp = p.assign(sb=pd.cut(p.s, np.r_[-SB[::-1][:-1], SB]))
    g = pp.groupby('sb', observed=True).agg(s=('s', 'mean'), d=('d_t_pk', 'median'), n=('s', 'size'),
                                            e=('d_t_pk', lambda v: 1.25 * v.std() / np.sqrt(len(v))))
    ax.errorbar(g.s, g.d, g.e, fmt='o', color='k')
    ax.axhline(0, color='#bbbbbb', lw=1); ax.axvline(0, color='#bbbbbb', lw=1)
    ax.set(xlabel='track slope along pad A -> pad B (dx/dz)', ylabel='t_peak(B) - t_peak(A) [ns]',
           title='Pad-pair time difference vs track slope (all pairs)')
    ax = axs[1]
    for z, col in cols.items():
        dd = C[(C.zone == z) & (C.est == 't_pk')]
        ax.errorbar(dd.s, dd.D, dd.e, fmt='o-', color=col, label=f'{z} ({int(dd.n_pos.sum()+dd.n_neg.sum())} pairs)')
    ax.set(xlabel='|slope|', ylabel='antisymmetrised time difference D [ns]', title='Data, by drift-gap class (stage 47)')
    ax.legend(frameon=False, fontsize=9)
    ax = axs[2]
    for (name, g), ls in zip(TC.groupby('model', sort=False), ('-', ':', '-.', '--', '-')):
        ax.plot(g.s, g.D_t_pk, ls, lw=2, label=f'{name}: T_d {g.T_d.iloc[0]:.0f} ns')
    for z, col in cols.items():
        dd = C[(C.zone == z) & (C.est == 't_pk')]
        ax.errorbar(dd.s, dd.D, dd.e, fmt='o', color=col, ms=4)
    ax.set(xlabel='|slope|', ylabel='D [ns]', title='Forward model with the data selection (lines) and data (points)')
    ax.legend(frameon=False, fontsize=8)
    for ax in axs:
        ax.grid(color='#e6e6e6')
    fig.tight_layout(); fig.savefig(f'{out}/microtpc.png', dpi=120); plt.close(fig)
    with open(f'{out}/summary.txt', 'w') as f:
        f.write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
