#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
24_gain_timing.py

Does the cosmic pulse-timing map follow the detector gain measured
independently with Fe55? (Test A1 of HANDOFF_gap_timing_2026-10-15.md.)

Stage 23 read the det1 timing map as time walk on a gain gradient. This stage
tests that pad by pad against an independent gain map:

  timing  stage 13 per-hit table (toa_hits<sfx>.csv, DREAM waveforms).
          Hits: amp >= --amp-min, not saturated, |t_parabola - median| <
          --t-win, pads with >= --min-hits hits. Per pad, medians of
            peak   = t_parabola
            frac30 = t_frac30
            rise   = t_parabola - t_frac30
  gain    Fe55 per-pad photopeak map (VMM readout, a different electronics
          chain): gainmap_<run>_gainmap.csv, valid pads only. Default run
          161639 = reduced VMM gain, so the high-gain pads are not clipped by
          saturation (in the original-gain run 150423, 296 of them are, which
          removes the lever arm). Other Fe55 runs are fitted as cross-checks.

Per quantity q: q_pad = a * ln(G_Fe55) + c_connector (least squares).
Reported: a in ns per e-fold of gain (bootstrap error over pads), the
variance explained within connectors, a label-shuffle null within connectors,
and the slope of q vs ln(amplitude) WITHIN a pad. Within one pad the gaps are
fixed, so that slope is the time-walk / amplitude-selection control.
Electronics control: rise residual vs DREAM channel within the connector.
Spatial control: neighbour correlation of the residuals after the gain term.
When stage 19 exists for the same config, its full-coverage peak-time map is
fitted too (peak only: stage 19 has no frac30).

Products (<Analysis>/<det>/<run>/<sub_run>/24_gain_timing/):
  offsets_vs_fe55_gain<sfx>.png    rise / peak / frac30 offsets vs ln G
  maps_gain_vs_rise<sfx>.png       Fe55 gain map beside the rise-offset map
  rise_residual_vs_channel<sfx>.png  electronics control
  slopes<sfx>.csv                  slope / R2 / null per Fe55 run and quantity
  pads_fe55_<run><sfx>.csv         per-pad table (main Fe55 run)
  summary<sfx>.txt                 all numbers
  inputs<sfx>.txt                  code version, input files, cuts

Usage:
  python3 24_gain_timing.py [run_key] [--fe55-run 161639]
          [--fe55-check 163851,150423] [--amp-min 300] [--no-veto-sparks]
"""

import argparse
import json
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import PolyCollection
from scipy.spatial import cKDTree

import p2_qa_config as qa
import p2_mapping as pmap
import p2_provenance as prov

# det1 Fe55 position scan, central source position (2026-09-09, after the
# July -> September remount), analysed by the VMM source-test scripts.
FE55_DIR = ('/media/ak271430/LaCie/Extras/Physics/Post-Doc-Saclay/data/'
            'LAB_Measurements/source_tests/Fe55_data/position_scan/'
            'central_position/analysis')

QUANT = {'rise': 'rise (peak − frac30)', 'peak': 'peak (parabola)',
         'f30': 'frac30'}
COL = {'rise': '#2a78d6', 'peak': '#eb6834', 'f30': '#1baf7a'}

plt.rcParams.update({'font.family': 'sans-serif',
                     'font.sans-serif': ['Arial', 'Liberation Sans', 'DejaVu Sans'],
                     'font.size': 11, 'figure.facecolor': 'white',
                     'axes.facecolor': 'white', 'savefig.facecolor': 'white',
                     'axes.spines.top': False, 'axes.spines.right': False})


def select(h, amp_min, t_win, t0):
    """t0: reference peak time, fixed once from the main amplitude cut -- a
    looser cut lets the late noise population (~1600 ns) drag the median."""
    s = h[h['amp_max'] >= amp_min]
    s = s[(s['t_parabola'] - t0).abs() < t_win].copy()
    s['rise'] = s['t_parabola'] - s['t_frac30']
    s['lnA'] = np.log(s['amp_max'])
    return s


def pad_table(hits, min_hits):
    n = hits.groupby('channel_id').size()
    hits = hits[hits['channel_id'].isin(n[n >= min_hits].index)]
    return hits.groupby('channel_id').agg(
        peak=('t_parabola', 'median'), f30=('t_frac30', 'median'),
        rise=('rise', 'median'), lnA=('lnA', 'median'),
        conn=('connector_N', 'first'), feu=('feu', 'first'),
        x=('pad_cx', 'first'), y=('pad_cy', 'first'), n_hits=('rise', 'size'))


def load_fe55(run):
    g = pd.read_csv(os.path.join(FE55_DIR, f'gainmap_{run}_gainmap.csv'))
    g = g[g['valid']].set_index('channel_id')
    return pd.DataFrame({'lnG': np.log(g['rel_gain']), 'fe_hits': g['n_hits']})


def fe55_conditions(run):
    try:
        with open(os.path.join(FE55_DIR, f'gainmap_{run}_gainmap.json')) as f:
            return json.load(f).get('conditions', '?')
    except OSError:
        return '?'


def fit(y, x, grp):
    """y = a*x + c_grp. Returns a, within-group R2, residuals."""
    X = pd.get_dummies(grp).astype(float)
    X['x'] = x
    X = X.to_numpy()
    beta, *_ = np.linalg.lstsq(X, y.to_numpy(), rcond=None)
    res = y.to_numpy() - X @ beta
    y0 = y - y.groupby(grp).transform('mean')
    return beta[-1], 1 - (res ** 2).sum() / (y0 ** 2).sum(), res


def within_pad_slope(hits, col, min_hits):
    n = hits.groupby('channel_id').size()
    hits = hits[hits['channel_id'].isin(n[n >= min_hits].index)]
    y = hits[col] - hits.groupby('channel_id')[col].transform('mean')
    x = hits['lnA'] - hits.groupby('channel_id')['lnA'].transform('mean')
    ok = y.notna()
    return (x[ok] * y[ok]).sum() / (x[ok] ** 2).sum()


def neighbour_corr(xy, r, radius=18.0, n_null=300):
    """Correlation of each pad's value with the mean of its neighbours
    (radius ~1.5 pad pitch), and the 95 % quantile of a label-shuffle null."""
    nb = cKDTree(xy).query_ball_point(xy, r=radius)
    nb = [[j for j in l if j != i] for i, l in enumerate(nb)]
    ok = np.array([len(l) > 0 for l in nb])

    def corr(v):
        nm = np.array([v[l].mean() if l else np.nan for l in nb])
        return np.corrcoef(v[ok], nm[ok])[0, 1]
    rng = np.random.default_rng(3)
    return corr(r), np.quantile([corr(rng.permutation(r)) for _ in range(n_null)], 0.95)


def draw_map(ax, verts, vals, lim, label, outline):
    ax.add_collection(PolyCollection(list(outline), facecolor='#eeeeee',
                                     edgecolor='white', lw=0.3))
    pc = PolyCollection(list(verts.loc[vals.index]), array=vals.to_numpy(),
                        cmap='RdBu_r', edgecolor='white', lw=0.3)
    pc.set_clim(-lim, lim)
    ax.add_collection(pc)
    ax.autoscale_view()
    ax.set_aspect('equal')
    ax.set_xlabel('x [mm]')
    ax.set_ylabel('y [mm]')
    plt.colorbar(pc, ax=ax, shrink=0.8, label=label)


def main():
    ap = argparse.ArgumentParser(description='Pad timing vs Fe55 gain.')
    ap.add_argument('run_key', nargs='?', default=qa.DEFAULT_RUN)
    ap.add_argument('--fe55-run', default='161639',
                    help='main Fe55 gain map (gainmap_<run>_gainmap.csv)')
    ap.add_argument('--fe55-check', default='163851,150423',
                    help='comma-separated Fe55 runs fitted as cross-checks')
    ap.add_argument('--amp-min', type=float, default=300.0)
    ap.add_argument('--t-win', type=float, default=150.0,
                    help='|t_parabola - median| window [ns]')
    ap.add_argument('--min-hits', type=int, default=40)
    ap.add_argument('--n-boot', type=int, default=1000)
    ap.add_argument('--veto-sparks', action=argparse.BooleanOptionalAction,
                    default=True, help='read the spark-vetoed stage 13 table')
    args = ap.parse_args()

    cfg = qa.get_config(args.run_key)
    print(cfg)
    out = cfg.out_dir('24_gain_timing')
    sfx = cfg.product_suffix(args.veto_sparks)
    toa_csv = os.path.join(cfg.OUT_BASE, '13_timing', f'toa_hits{sfx}.csv')
    if not os.path.isfile(toa_csv):
        raise SystemExit(f'no stage 13 table {toa_csv} -- run 13_timing_waveforms.py first')
    surf_csv = os.path.join(cfg.OUT_BASE, '19_timing_surface',
                            f'timing_surface{sfx}.csv')
    fe_runs = [args.fe55_run] + [r for r in args.fe55_check.split(',') if r]

    lines = []

    def log(s=''):
        print(s)
        lines.append(s)

    h = pd.read_csv(toa_csv)
    h = h[~h['saturated']]
    chunks = sorted(h['chunk'].unique().tolist())
    t0 = h.loc[h['amp_max'] >= args.amp_min, 't_parabola'].median()
    sel = select(h, args.amp_min, args.t_win, t0)
    pads = pad_table(sel, args.min_hits)
    log(f'{cfg.DET_TAG} {cfg.RUN}/{cfg.SUB_RUN}  (config {cfg.KEY}, chunks {chunks})')
    log(f'cosmic: {len(sel):,} hits pass the cuts, {len(pads)} pads with '
        f'>= {args.min_hits} hits ({pads.n_hits.sum():,} hits on them)')

    cols = {'rise': 'rise', 'peak': 't_parabola', 'f30': 't_frac30'}
    wp = {q: within_pad_slope(sel, c, args.min_hits) for q, c in cols.items()}
    log('within-pad slope vs ln(amplitude) [ns / e-fold]: '
        + ', '.join(f'{q} {wp[q]:+.2f}' for q in QUANT))

    rng = np.random.default_rng(1)
    rows, tabs = [], {}
    for run in fe_runs:
        m = pads.join(load_fe55(run), how='inner')
        tabs[run] = m
        b, r2, _ = fit(m['lnA'], m['lnG'], m['conn'])
        log(f'\nFe55 {run} ({fe55_conditions(run)}): {len(m)} matched pads; '
            f'cosmic ln(amp) vs ln G: slope {b:.2f}, R2 {r2:.2f}')
        n_null = 1000 if run == args.fe55_run else 100
        for q in QUANT:
            a, r2, _ = fit(m[q], m['lnG'], m['conn'])
            boot = []
            for _ in range(args.n_boot):
                i = rng.integers(0, len(m), len(m))
                boot.append(fit(m[q].iloc[i], m['lnG'].iloc[i], m['conn'].iloc[i])[0])
            shuf = [fit(m[q], m.groupby('conn')['lnG'].transform(
                lambda v: rng.permutation(v.to_numpy())), m['conn'])[1]
                for _ in range(n_null)]
            rows.append(dict(fe55_run=run, quantity=q, n_pads=len(m), slope_ns=a,
                             slope_err=np.std(boot), r2_within_conn=r2,
                             null_r2_95=np.quantile(shuf, 0.95),
                             within_pad_slope_ns=wp[q]))
            log(f'  {q:5s} slope {a:+6.2f} ± {np.std(boot):.2f} ns/e-fold   '
                f'R2 {r2:.2f}  (shuffle null 95 %: {np.quantile(shuf, .95):.3f})')
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(out, f'slopes{sfx}.csv'), index=False)
    main_res = res[res['fe55_run'] == args.fe55_run].set_index('quantity')

    m = tabs[args.fe55_run]
    fe = load_fe55(args.fe55_run)
    log(f'\nrobustness (Fe55 {args.fe55_run}):')
    hg = h.join(fe, on='channel_id')
    hg = hg[hg['lnG'].notna()]
    variants = {f'amp >= {args.amp_min / 2:g}': select(h, args.amp_min / 2, args.t_win, t0),
                f'amp >= {args.amp_min * 5 / 3:g}': select(h, args.amp_min * 5 / 3, args.t_win, t0)}
    eq = select(hg, args.amp_min / 2, args.t_win, t0)
    variants[f'amp >= {args.amp_min:g} x G_rel (equal-charge cut)'] = \
        eq[eq['amp_max'] >= args.amp_min * np.exp(eq['lnG'])]
    for name, v in variants.items():
        mv = pad_table(v, args.min_hits).join(fe, how='inner')
        log(f'  {name:38s} {len(mv):3d} pads: ' + ', '.join(
            f'{q} {fit(mv[q], mv["lnG"], mv["conn"])[0]:+.1f}' for q in QUANT))
    log('  cosmic ln(amp) as regressor (no Fe55), all pads: ' + ', '.join(
        f'{q} {fit(pads[q], pads["lnA"], pads["conn"])[0]:+.1f}' for q in QUANT))
    rb = []
    for _ in range(args.n_boot):
        mm = m.iloc[rng.integers(0, len(m), len(m))]
        rb.append(fit(mm['f30'], mm['lnG'], mm['conn'])[0]
                  / fit(mm['peak'], mm['lnG'], mm['conn'])[0])
    ratio = main_res.loc['f30', 'slope_ns'] / main_res.loc['peak', 'slope_ns']
    log(f'  frac30 / peak slope ratio: {ratio:.2f} ± {np.std(rb):.2f}')

    if os.path.isfile(surf_csv):
        s19 = pd.read_csv(surf_csv).set_index('channel_id')
        s19 = s19[s19['n'] >= args.min_hits].join(fe, how='inner')
        a, r2, _ = fit(s19['median_ns'], s19['lnG'], s19['connector_N'])
        log(f'  stage 19 peak time (all chunks, {len(s19)} pads): '
            f'{a:+.1f} ns/e-fold, R2 {r2:.2f}')
    else:
        log(f'  (no stage 19 table for this config: {os.path.basename(surf_csv)})')

    # electronics control
    ct = pmap.build_channel_table(cfg.run_config_path, cfg.MAP_CSV_PATH,
                                  det_type=cfg.DET_TYPE, det_name=cfg.DET_NAME,
                                  strategy='reverse',
                                  drop_connectors=cfg.DEAD_CONNECTORS,
                                  strategy_overrides=cfg.STRATEGY_OVERRIDES)
    tiles, verts = pmap.pad_tiles(ct)
    verts = pd.Series(list(verts), index=tiles['channel_id'].to_numpy())
    chan = ct[ct['mapped']].drop_duplicates('channel_id').set_index('channel_id')
    m = m.join(chan[['within']])
    m['res_rise'] = fit(m['rise'], m['lnG'], m['conn'])[2]
    log('\nelectronics control (rise):')
    log('  connector medians of rise, all pads [ns]: ' + ', '.join(
        f'{c}: {v:.1f}' for c, v in pads.groupby('conn')['rise'].median().items()))
    g8 = m.groupby(m['within'] // 8)['res_rise']
    log('  residual after gain fit by DREAM channel block (within // 8): '
        + ', '.join(f'{k}: {v:+.1f}' for k, v in g8.median().items()))
    log(f'  fraction of residual variance in the 8-channel blocks: '
        f'{np.var(g8.transform("mean")) / np.var(m["res_rise"]):.2f}')
    log(f'  residual std {m.res_rise.std():.2f} ns vs rise-offset std '
        f'{(m.rise - m.groupby("conn").rise.transform("mean")).std():.2f} ns')

    log('\nspatial structure (neighbour correlation, 18 mm; shuffle null 95 %):')
    xy = m[['x', 'y']].to_numpy()
    for q in QUANT:
        c, n95 = neighbour_corr(xy, fit(m[q], m['lnG'], m['conn'])[2])
        log(f'  {q:5s} residual: {c:+.2f} (null {n95:.2f})')
    c, n95 = neighbour_corr(xy, (m.lnG - m.groupby('conn').lnG.transform('mean')).to_numpy())
    log(f'  ln G itself : {c:+.2f} (null {n95:.2f})')

    q05, q95 = np.quantile(fe['lnG'], [0.05, 0.95])
    log(f'\nFe55 {args.fe55_run}, all {len(fe)} valid pads: std(ln G) '
        f'{fe.lnG.std():.3f}, 5-95 % range of ln G {q95 - q05:.3f}, std/mean(G) '
        f'{np.exp(fe.lnG).std() / np.exp(fe.lnG).mean():.3f}')
    m.to_csv(os.path.join(out, f'pads_fe55_{args.fe55_run}{sfx}.csv'))

    # ---------------- figures ----------------
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.2), sharex=True)
    xs = np.linspace(m.lnG.min(), m.lnG.max(), 2)
    x0 = m.lnG.mean()
    for ax, q in zip(axs, QUANT):
        r = main_res.loc[q]
        # connector constants removed, connector-mean gain kept on the x axis
        off = (m[q] - m.groupby('conn')[q].transform('mean')
               + (m.lnG.groupby(m.conn).transform('mean') - x0) * r.slope_ns)
        ax.scatter(m.lnG, off, s=14, color=COL[q], alpha=0.7, lw=0)
        ax.plot(xs, r.slope_ns * (xs - x0), color='black', lw=2,
                label=f'between pads: {r.slope_ns:+.1f} ± {r.slope_err:.1f} ns/e-fold')
        ax.plot(xs, wp[q] * (xs - x0), color='black', lw=2, ls='--',
                label=f'within a pad: {wp[q]:+.1f} ns/e-fold')
        ax.set_title(f'{QUANT[q]}   R² = {r.r2_within_conn:.2f}')
        ax.set_xlabel('ln(Fe55 relative gain)')
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), frameon=False,
                  fontsize=9)
    axs[0].set_ylabel('pad offset, connector removed [ns]')
    fig.suptitle(f'{cfg.DET_TAG}: cosmic pulse timing vs Fe55 gain '
                 f'({len(m)} pads, Fe55 run {args.fe55_run})')
    fig.tight_layout()
    fig.savefig(os.path.join(out, f'offsets_vs_fe55_gain{sfx}.png'), dpi=150,
                bbox_inches='tight')
    plt.close(fig)

    outline = np.stack(list(verts.to_numpy()))
    fig, axs = plt.subplots(1, 2, figsize=(13, 5.2))
    g = fe['lnG'][fe.index.isin(verts.index)]
    draw_map(axs[0], verts, g, np.quantile(np.abs(g), 0.98), 'ln(relative gain)',
             outline)
    axs[0].set_title(f'Fe55 gain (run {args.fe55_run}, VMM)')
    ro = pads.rise - pads.groupby('conn').rise.transform('median')
    draw_map(axs[1], verts, -ro, np.quantile(np.abs(ro), 0.98),
             '−(rise offset) [ns]', outline)
    axs[1].set_title(f'cosmic rise time, sign flipped (DREAM, chunks {chunks})')
    fig.suptitle(f'{cfg.DET_TAG}: Fe55 gain and cosmic rise time')
    fig.tight_layout()
    fig.savefig(os.path.join(out, f'maps_gain_vs_rise{sfx}.png'), dpi=150,
                bbox_inches='tight')
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.scatter(m['within'], m['res_rise'], s=12, color=COL['rise'], alpha=0.7, lw=0)
    med = g8.median()
    ax.step(np.r_[med.index * 8, 64], np.r_[med.values, med.values[-1]],
            where='post', color='black', lw=2)
    ax.axhline(0, color='#999999', lw=0.8)
    ax.set_xlabel('DREAM channel within connector')
    ax.set_ylabel('rise residual after gain fit [ns]')
    ax.set_title('electronics control: rise residual after the gain term')
    fig.tight_layout()
    fig.savefig(os.path.join(out, f'rise_residual_vs_channel{sfx}.png'), dpi=150)
    plt.close(fig)

    with open(os.path.join(out, f'summary{sfx}.txt'), 'w') as f:
        f.write('\n'.join(lines) + '\n')
    inputs = {'stage 13 per-hit table': toa_csv,
              'stage 13 chunks': chunks,
              'time window [h]': f'{cfg.T_MIN_H or 0:g} - {cfg.T_MAX_H or "end"}'}
    for run in fe_runs:
        inputs[f'Fe55 {run}'] = os.path.join(FE55_DIR, f'gainmap_{run}_gainmap.csv')
        inputs[f'Fe55 {run} conditions'] = fe55_conditions(run)
    if os.path.isfile(surf_csv):
        inputs['stage 19 surface (peak cross-check)'] = surf_csv
    prov.write(out, sfx, __file__, inputs,
               settings={'amp_min [ADC]': args.amp_min, 't_win [ns]': args.t_win,
                         'min_hits per pad': args.min_hits,
                         'main Fe55 run': args.fe55_run, 'bootstrap': args.n_boot,
                         'spark veto': args.veto_sparks})
    print(f'\nwrote {out}')


if __name__ == '__main__':
    main()
