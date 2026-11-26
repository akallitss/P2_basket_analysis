#!/usr/bin/env python3
"""32_misreco_anatomy.py -- why are fired events mis-reconstructed?

Usage:  python3 32_misreco_anatomy.py <run_key> [--no-veto-sparks]

For every active-area track (stage 06 list) the event's P2 hits are clustered
exactly like the pipeline (link CLUSTER_LINK_MM, max-charge cluster wins).
  A = the cluster closest to the track (if any pad within MATCH_R fired)
  B = the winning (reconstructed) cluster
Classes of mis-reco (reco > MATCH_R):
  wrong cluster won  : A exists, B != A        -> B: crosstalk copy of A?
                       (B pads are Dream neighbours ch+-1/+-2 of A pads),
                       noise, 2nd particle? B distance, amplitude, timing
  only distant hits  : nothing within MATCH_R  -> compare with the random
                       noise-hit rate measured on tracks that MISS the det
  resolution tail    : reco in (R, 2R]
Pointing check: residual core/tails vs track angle and chi2, and the mean
residual vector vs position (systematic misalignment).
Outputs: 32_misreco_anatomy/misreco_anatomy<sfx>.png + csv
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
import awkward as ak
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2_qa_config as qa      # noqa: E402
import p2_mapping as pmap      # noqa: E402
import p2_align as pa          # noqa: E402
import p2_io as p2io           # noqa: E402
from M3RefTracking import M3RefTracking  # noqa: E402

INK, INK2, GRID = '#1f1f1f', '#5c5c5c', '#e6e6e6'
C1, C2, C3, C4 = '#257de0', '#f47230', '#00b894', '#6a4ce2'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run_key', nargs='?', default=qa.DEFAULT_RUN)
    ap.add_argument('--veto-sparks', action=argparse.BooleanOptionalAction, default=True)
    args = ap.parse_args()
    cfg = qa.get_config(args.run_key)
    out_dir = cfg.out_dir('32_misreco_anatomy')
    sfx = cfg.product_suffix(args.veto_sparks)
    R, LINK = cfg.MATCH_R, cfg.CLUSTER_LINK_MM

    ct = pmap.build_channel_table(cfg.run_config_path, cfg.MAP_CSV_PATH, det_type=cfg.DET_TYPE,
                                  det_name=cfg.DET_NAME, drop_connectors=cfg.DEAD_CONNECTORS,
                                  drop_channels=cfg.DEAD_CHANNELS, cluster_link=cfg.CLUSTER_LINK_MM,
                                  dead_guard=cfg.DEAD_GUARD_MM, strategy_overrides=cfg.STRATEGY_OVERRIDES)
    padc = ct[ct['mapped']].drop_duplicates('channel_id').reset_index(drop=True)
    ed = cfg.out_dir('06_efficiency')
    fp = pd.read_csv(f'{ed}/pad_footprint{sfx}.csv')
    P = np.c_[padc.pad_cx, padc.pad_cy, np.ones(len(padc))]
    coef, *_ = np.linalg.lstsq(P, fp[['x', 'y']].to_numpy(), rcond=None)
    rays = pd.read_csv(f'{ed}/ray_hit_miss_list{sfx}.csv')
    for c in ('in_active', 'has_any', 'within'):
        rays[c] = rays[c].astype(str).str.lower().isin(('true', '1'))
    act = rays[rays['in_active']].set_index('eventId')
    N = len(act)

    # M3 angle / chi2 for every single-track event (incl. tracks that miss P2)
    m3 = M3RefTracking(os.path.join(cfg.m3_tracking_dir, ''), single_track=True,
                       chi2_cut=qa.M3_CHI2_CUT, min_nclus=qa.M3_MIN_NCLUS)
    rd = m3.ray_data

    def col(n):
        return np.asarray(ak.to_numpy(ak.ravel(rd[n])), dtype=float)
    ang = pd.DataFrame({'eventId': np.asarray(ak.to_numpy(ak.ravel(rd['evn'])), dtype=np.int64),
                        'theta': np.degrees(np.arctan2(np.hypot(col('X_Up') - col('X_Down'), col('Y_Up') - col('Y_Down')),
                                                       np.abs(col('Z_Up') - col('Z_Down')))),
                        'chi2': np.maximum(col('Chi2X'), col('Chi2Y'))}).drop_duplicates('eventId').set_index('eventId')

    # hits (same cuts as the pipeline)
    mp = ct[ct['mapped']].set_index(['feu', 'channel'])[['pad_cx', 'pad_cy']]
    parts = []
    for df in p2io.iter_hits(cfg.combined_hits_dir, ['eventId', 'feu', 'channel', 'amplitude', 'time_of_max', 'saturated'],
                             ct.attrs['feus'], progress=False, t_max_h=cfg.T_MAX_H, t_min_h=cfg.T_MIN_H,
                             min_amp=cfg.MIN_AMP):
        df = pd.DataFrame(df) if not isinstance(df, pd.DataFrame) else df
        idx = pd.MultiIndex.from_arrays([df['feu'].astype(int), df['channel'].astype(int)])
        ok = idx.isin(mp.index)
        df = df[ok].copy()
        xy = mp.reindex(idx[ok]).to_numpy()
        q = np.c_[xy, np.ones(len(xy))] @ coef
        df['x'], df['y'] = q[:, 0], q[:, 1]
        parts.append(df)
    h = pd.concat(parts, ignore_index=True)
    h['feu'] = h['feu'].astype(int); h['channel'] = h['channel'].astype(int)
    h = h.sort_values('amplitude', ascending=False).drop_duplicates(['eventId', 'feu', 'channel'])

    # ---- noise occupancy: tracks that cross NOTHING of the detector (> 40 mm out) ----
    allr = rays.copy()
    dfp, _ = cKDTree(fp[['x', 'y']].to_numpy()).query(allr[['x', 'y']].to_numpy())
    outside = allr[dfp > 40]
    ev_hit = set(h['eventId'])
    p_noise = outside['eventId'].isin(ev_hit).mean() if len(outside) else np.nan
    print(f'random-hit rate: tracks > 40 mm outside the pad array: {len(outside):,}, '
          f'P(any P2 hit) = {100*p_noise:.2f}%')

    # ---- per-event clustering for the mis-reco tracks + a sample of good ones ----
    mis = act[act['has_any'] & ~act['within']]
    good = act[act['within']].sample(min(6000, int(act['within'].sum())), random_state=2)
    sel_ev = set(mis.index) | set(good.index)
    hs = h[h['eventId'].isin(sel_ev)]
    rows = []
    for e, g in hs.groupby('eventId', sort=False):
        tx, ty = act.loc[e, 'x'], act.loc[e, 'y']
        x, y, a = g.x.to_numpy(), g.y.to_numpy(), g.amplitude.to_numpy()
        ch, feu, tm, sat = g.channel.to_numpy(), g.feu.to_numpy(), g.time_of_max.to_numpy(), g.saturated.to_numpy()
        n = len(x)
        if n > 1:
            pr = cKDTree(np.c_[x, y]).query_pairs(LINK, output_type='ndarray')
            G = coo_matrix((np.ones(len(pr)), (pr[:, 0], pr[:, 1])), shape=(n, n)) if len(pr) else coo_matrix((n, n))
            nc, lab = connected_components(G, directed=False)
        else:
            nc, lab = 1, np.zeros(1, int)
        csum = np.bincount(lab, weights=np.clip(a, 0, None), minlength=nc)
        cx = np.bincount(lab, weights=np.clip(a, 0, None) * x, minlength=nc) / np.maximum(csum, 1e-9)
        cy = np.bincount(lab, weights=np.clip(a, 0, None) * y, minlength=nc) / np.maximum(csum, 1e-9)
        cd = np.hypot(cx - tx, cy - ty)
        dpad = np.hypot(x - tx, y - ty)
        iB = int(np.argmax(csum))
        near = dpad <= R
        iA = int(lab[near][np.argmin(dpad[near])]) if near.any() else -1
        r = dict(eventId=e, kind='good' if e in good.index else 'mis', n_clu=nc, n_pads=n,
                 dB=cd[iB], sumB=csum[iB], npB=int((lab == iB).sum()), hasA=iA >= 0)
        if iA >= 0 and iA != iB:
            A_ch = set(zip(feu[lab == iA], ch[lab == iA]))
            B_ch = list(zip(feu[lab == iB], ch[lab == iB]))
            nb = sum(any((f_, c_ + o) in A_ch for o in (-2, -1, 1, 2)) for f_, c_ in B_ch)
            r.update(sumA=csum[iA], npA=int((lab == iA).sum()), dA=cd[iA],
                     B_neigh_frac=nb / len(B_ch),
                     dt=float(np.median(tm[lab == iB]) - np.median(tm[lab == iA])),
                     satA=bool(sat[lab == iA].any()), satB=bool(sat[lab == iB].any()),
                     maxA=float(a[lab == iA].max()), maxB=float(a[lab == iB].max()))
        elif iA >= 0 and len(np.unique(lab)) > 1:      # good event that also had another cluster
            others = [k for k in range(nc) if k != iA]
            k2 = max(others, key=lambda k: csum[k])
            A_ch = set(zip(feu[lab == iA], ch[lab == iA]))
            B_ch = list(zip(feu[lab == k2], ch[lab == k2]))
            nb = sum(any((f_, c_ + o) in A_ch for o in (-2, -1, 1, 2)) for f_, c_ in B_ch)
            r.update(sumA=csum[iA], sum2=csum[k2], B_neigh_frac=nb / len(B_ch),
                     satA=bool(sat[lab == iA].any()), maxA=float(a[lab == iA].max()))
        rows.append(r)
    t = pd.DataFrame(rows).set_index('eventId')
    t = t.join(act[['resid', 'x', 'y']]).join(ang)
    t.to_csv(f'{out_dir}/misreco_events{sfx}.csv')
    m = t[t.kind == 'mis']
    tail = m[m.resid <= 2 * R]
    wrong = m[(m.resid > 2 * R) & m.hasA]
    far = m[(m.resid > 2 * R) & ~m.hasA]
    gd = t[t.kind == 'good']
    print(f'\nmis-reco {len(m)} ({100*len(m)/N:.2f}%): tail {100*len(tail)/N:.2f}%, wrong cluster '
          f'{100*len(wrong)/N:.2f}%, only distant {100*len(far)/N:.2f}%')

    # ---- wrong cluster won: what is B? ----
    w = wrong
    xt = (w.B_neigh_frac > 0.5)
    print('\nWRONG CLUSTER WON (A = cluster at the track, B = winner):')
    print(f'  B made of Dream neighbours (ch+-1/+-2) of A pads: {100*xt.mean():.1f}%  '
          f'(good events with a 2nd cluster: {100*(gd.B_neigh_frac.dropna() > 0.5).mean():.1f}% of their 2nd clusters)')
    print(f'  A saturated: {100*w.satA.mean():.1f}%  (good events: {100*gd.satA.dropna().mean():.1f}%)')
    print(f'  B/A summed charge: median {np.median(w.sumB / w.sumA):.2f};  B pads median {w.npB.median():.0f}, '
          f'A pads median {w.npA.median():.0f}')
    print(f'  |t_B - t_A| median {np.median(np.abs(w.dt)):.0f} ns; within 60 ns: {100*(np.abs(w.dt) <= 60).mean():.1f}%')
    print(f'  B distance from track: median {w.dB.median():.0f} mm (quartiles {w.dB.quantile(.25):.0f}-{w.dB.quantile(.75):.0f})')
    # crosstalk-explained share and the cost
    print(f'  -> crosstalk-like (B = neighbours of A, in time): '
          f'{100*(xt & (np.abs(w.dt) <= 60)).mean():.1f}% of wrong-cluster events')
    # what if A had been chosen: was A within R? (by construction yes -> A would be reconstructed)
    good2 = gd[gd.sum2.notna()] if 'sum2' in gd else gd.iloc[0:0]
    print(f'  good events with a competing 2nd cluster: {100*len(good2)/max(len(gd),1):.1f}%; '
          f'2nd/1st charge median {np.median(good2.sum2 / good2.sumA) if len(good2) else np.nan:.2f}')

    # ---- only distant hits ----
    print('\nONLY DISTANT HITS (nothing fired within R of the track):')
    print(f'  measured {100*len(far)/N:.2f}% of tracks; random-noise expectation = P(no signal at track) x '
          f'P(noise hit) = {100*(1 - act.has_any.mean() + len(far)/N):.2f}% x {100*p_noise:.2f}% '
          f'= {100*(1 - act.has_any.mean() + len(far)/N) * p_noise:.3f}%')
    print(f'  nearest fired pad distance: median {far.dB.median():.0f} mm; B pads median {far.npB.median():.0f}; '
          f'B charge median {far.sumB.median():.0f} vs good-event cluster {gd.sumB.median():.0f}')

    # ---- pointing check ----
    ok = act[act['within']].join(ang)
    print('\nPOINTING: residual of reconstructed tracks (|r| <= R):')
    for nm, cc, bins in (('theta', 'theta', [0, 5, 10, 15, 20, 30]), ('chi2', 'chi2', [0, 0.5, 1, 2, 3, 5])):
        g = ok.groupby(pd.cut(ok[cc], bins), observed=False)['resid'].agg(['size', 'median'])
        g['sigma_axis'] = g['median'] / np.sqrt(2 * np.log(2))
        print(f'  vs {nm}:'); print(g.round(2).to_string())
    tail_frac = act[act['has_any']].pipe(lambda d: ((d.resid > R) & (d.resid <= 2 * R)).mean())
    print(f'  fired tracks with reco in ({R:g},{2*R:g}] mm: {100*tail_frac:.2f}%')
    # systematic residual vs position: mean dx, dy in 6x6 regions
    rr = act[act['within']].copy()
    reco = None

    # ---- figure ----
    plt.rcParams.update({'font.size': 9.5, 'axes.edgecolor': INK2, 'xtick.color': INK2, 'ytick.color': INK2})
    fig, axs = plt.subplots(2, 3, figsize=(17, 9.5))
    ax = axs[0, 0]
    bins = np.linspace(0, 450, 46)
    ax.hist(w.dB.clip(upper=449), bins, color=C2, alpha=.85, label=f'wrong cluster won: B (n={len(w)})')
    ax.hist(far.dB.clip(upper=449), bins, histtype='step', lw=2, color=C4, label=f'only distant hits (n={len(far)})')
    ax.set_xlabel('distance of the reconstructed cluster from the track [mm]'); ax.set_ylabel('events')
    ax.set_title('A  How far away is the winning cluster?', loc='left', fontweight='bold', color=INK)
    ax.legend(frameon=False, fontsize=8.5)
    ax = axs[0, 1]
    cats = ['B = Dream neighbours\nof A (crosstalk)', 'B in time with A\n(|dt| <= 60 ns)', 'A saturated']
    v_w = [100 * xt.mean(), 100 * (np.abs(w.dt) <= 60).mean(), 100 * w.satA.mean()]
    v_g = [100 * (gd.B_neigh_frac.dropna() > 0.5).mean(), np.nan, 100 * gd.satA.dropna().mean()]
    xx = np.arange(3)
    ax.bar(xx - 0.18, v_w, 0.34, color=C2, label='wrong-cluster events')
    ax.bar(xx + 0.18, v_g, 0.34, color=C1, label='reconstructed events (2nd cluster / A)')
    for i, (a_, b_) in enumerate(zip(v_w, v_g)):
        ax.text(i - 0.18, a_ + 1, f'{a_:.0f}%', ha='center', fontsize=8.5)
        if np.isfinite(b_):
            ax.text(i + 0.18, b_ + 1, f'{b_:.0f}%', ha='center', fontsize=8.5)
    ax.set_xticks(xx); ax.set_xticklabels(cats, fontsize=8.5); ax.set_ylim(0, 110)
    ax.set_ylabel('% of events'); ax.legend(frameon=False, fontsize=8.5)
    ax.set_title('B  What is the winning cluster B?', loc='left', fontweight='bold', color=INK)
    ax = axs[0, 2]
    lb = np.linspace(0, 4, 41)
    ax.hist(np.log10(w.sumA.clip(lower=1)), lb, histtype='step', lw=2, color=C3, label='A (at the track)')
    ax.hist(np.log10(w.sumB.clip(lower=1)), lb, histtype='step', lw=2, color=C2, label='B (winner)')
    ax.hist(np.log10(gd.sumB.clip(lower=1)), lb, histtype='step', lw=1.5, ls='--', color=C1, density=False,
            weights=np.full(len(gd), len(w) / max(len(gd), 1)), label='reconstructed cluster (scaled)')
    ax.set_xlabel('log10 cluster charge [ADC]'); ax.set_ylabel('events')
    ax.set_title('C  Cluster charges in wrong-cluster events', loc='left', fontweight='bold', color=INK)
    ax.legend(frameon=False, fontsize=8.5)
    for k, (ax, cc, bins, xl) in enumerate(((axs[1, 0], 'theta', [0, 5, 10, 15, 20, 30], 'track angle [deg]'),
                                            (axs[1, 1], 'chi2', [0, 0.5, 1, 2, 3, 5], 'M3 chi2_max'))):
        tt = act.join(ang)
        tt['wrong'] = tt.index.isin(wrong.index); tt['far'] = tt.index.isin(far.index)
        tt['tail'] = tt.index.isin(tail.index)
        g = tt.groupby(pd.cut(tt[cc], bins), observed=False)[['wrong', 'far', 'tail']].mean() * 100
        nn_ = tt.groupby(pd.cut(tt[cc], bins), observed=False).size()
        xx = np.arange(len(g))
        for colname, c_, lab_ in (('wrong', C2, 'wrong cluster won'), ('far', C4, 'only distant hits'),
                                  ('tail', C3, f'resolution tail ({R:g}-{2*R:g} mm)')):
            e_ = np.sqrt(g[colname] / 100 * (1 - g[colname] / 100) / nn_) * 100
            ax.errorbar(xx, g[colname], e_, fmt='o-', color=c_, lw=2, ms=5, label=lab_)
        ax.set_xticks(xx); ax.set_xticklabels([str(i) for i in g.index], fontsize=8, rotation=15)
        ax.set_xlabel(xl); ax.set_ylabel('% of tracks'); ax.legend(frameon=False, fontsize=8.5)
        ax.set_title(f'{"D" if k == 0 else "E"}  Mis-reco classes vs {xl}', loc='left', fontweight='bold', color=INK)
    ax = axs[1, 2]
    g = ok.groupby(pd.cut(ok['theta'], [0, 5, 10, 15, 20, 30]), observed=False)['resid'].median() / np.sqrt(2 * np.log(2))
    ax.plot(np.arange(len(g)), g.values, 'o-', color=C1, lw=2, label='residual sigma per axis (reco tracks)')
    ax.set_xticks(np.arange(len(g))); ax.set_xticklabels([str(i) for i in g.index], fontsize=8, rotation=15)
    ax.axhline(np.sqrt(np.median(ok.resid / np.sqrt(2 * np.log(2))) ** 2), color=INK2, ls=':', lw=1)
    ax.set_xlabel('track angle [deg]'); ax.set_ylabel('sigma [mm]'); ax.legend(frameon=False, fontsize=8.5)
    ax.set_title('F  Pointing: residual width vs angle', loc='left', fontweight='bold', color=INK)
    for a_ in axs.ravel():
        a_.grid(True, color=GRID, lw=0.8); a_.set_axisbelow(True)
        for s_ in ('top', 'right'):
            a_.spines[s_].set_visible(False)
    fig.suptitle(f'{cfg.DET_NAME} anatomy of the mis-reconstructions — {cfg.RUN}/{cfg.SUB_RUN}  '
                 f'(wrong cluster {100*len(wrong)/N:.2f}%, only distant {100*len(far)/N:.2f}%, '
                 f'tail {100*len(tail)/N:.2f}%)', x=0.01, ha='left', color=INK, fontsize=12)
    fig.tight_layout()
    fig.savefig(f'{out_dir}/misreco_anatomy{sfx}.png', dpi=140)
    plt.close(fig)
    print(f'  saved {out_dir}/misreco_anatomy{sfx}.png')


if __name__ == '__main__':
    main()
