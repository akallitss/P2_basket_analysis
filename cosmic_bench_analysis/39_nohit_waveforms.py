#!/usr/bin/env python3
"""39_nohit_waveforms.py -- what is on the pads under a track that left NO hit?

Hit finding (analyze_waveforms) needs a pulse > 5 sigma (gate noise) for >= 3
samples, so a 3-4 sigma pulse is dropped. Question: in the no-hit events (3.2 %
of active tracks, ~5 % for near-vertical tracks) is there a sub-threshold pulse
on the pad(s) under the track, or nothing at all?

Samples (active area exactly as stage 06 with the run config):
  NOHIT  every active track with no P2 hit in the event
  HIT    a random sample of active tracks reconstructed within MATCH_R
  NOISE  the NOHIT pads read in an UNRELATED event (eventId of the track +-K
         positions in the decoded file) -> the noise-only reference
Pads read per track: every live pad whose tile is within PAD_R mm of the track.

Per (track, pad), on the common-noise-subtracted waveform:
  peak   max sample in the signal window   sig = peak / channel noise RMS
  mf     matched-filter amplitude: least-squares scale of the average HIT pulse
         shape in the window; mf_sig = mf / its noise sigma (from NOISE)
Per track: the pad with the largest mf_sig. Plus stacked average waveforms
per class and angle bin -> mean collected charge in no-hit events.

Also checks that every NOHIT event exists in the decoded data of both P2 FEUs
(a missing event would be a DAQ loss, not physics).

Usage: python3 39_nohit_waveforms.py [run_key] [--n-hit 4000] [--pad-r 3]
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2_study_cache as sc     # noqa: E402
import p2_study_common as cm    # noqa: E402
import p2_align as pa           # noqa: E402
import p2_waveforms as pw       # noqa: E402

INK2, GRID = '#5c5c5c', '#e6e6e6'
C_NH, C_HIT, C_NOISE = '#ec4647', '#257de0', '#9c9c9c'
TH_BINS = [(0, 6), (6, 12), (12, 90)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run_key', nargs='?', default='det2_long_repos')
    ap.add_argument('--n-hit', type=int, default=4000)
    ap.add_argument('--pad-r', type=float, default=3.0)
    ap.add_argument('--shift', type=int, default=997, help='event offset for the NOISE control')
    args = ap.parse_args()
    hits, rays, g, cfg = sc.load(args.run_key)
    out = cfg.out_dir('39_nohit_waveforms')
    L = []

    # ---- tracks, alignment and active area exactly as stage 06 (run config) ----
    ct = sc.channel_table(cfg)
    ct.attrs['edge_guard'] = cfg.EDGE_GUARD_MM
    clu = pd.read_pickle(cfg.out_dir('34_study_cache') + '/clusters_cc20.pkl')
    ev = rays['eventId'].to_numpy()
    x, y = pa.project_to_z(rays, cfg.det_plane_z()); xy = np.column_stack([x, y])
    cen = clu.loc[clu.groupby('eventId')['q'].idxmax()]
    m = pd.DataFrame({'eventId': ev, 'x': x, 'y': y}).merge(cen, on='eventId')
    m = m[(m.x.abs() < 300) & (m.y.abs() < 300)]
    T = pa.fit_transform(m.x, m.y, m.x_pad, m.y_pad, trim_mm=cfg.FIT_TRIM_MM)
    padc = ct[ct['mapped']].drop_duplicates('channel_id').reset_index(drop=True)
    pxy = np.column_stack(T.apply(padc.pad_cx.to_numpy(), padc.pad_cy.to_numpy()))
    nn, _ = cKDTree(pxy).query(xy)
    act = pa.mask_dead_cells(nn, xy, ct, T) <= cfg.ACTIVE_R
    d, _ = cm.nearest(ev, xy, cm.apply_T(T, clu[['x_pad', 'y_pad']].to_numpy()), clu.eventId.to_numpy())
    within = d <= cfg.MATCH_R
    has_any = np.isin(ev, hits.eventId.to_numpy())
    tx, ty, theta = cm.ray_dirs(rays)
    L.append(f'active tracks {act.sum():,}: within {100*within[act].mean():.2f} %, any {100*has_any[act].mean():.2f} %, '
             f'no-hit {100*(~has_any[act]).mean():.2f} %')

    # pads (feu, channel) whose tile is within PAD_R of each track
    _, verts = __import__('p2_mapping').pad_tiles(ct)
    vx, vy = T.apply(verts[:, :, 0].ravel(), verts[:, :, 1].ravel())
    tiles = np.stack([vx.reshape(-1, 4), vy.reshape(-1, 4)], axis=2)
    tpads = __import__('p2_mapping').pad_tiles(ct)[0].reset_index(drop=True)
    tc = tiles.mean(1)
    rng = np.random.default_rng(7)
    i_nh = np.flatnonzero(act & ~has_any)
    i_hit = rng.choice(np.flatnonzero(act & within), min(args.n_hit, int((act & within).sum())), replace=False)
    want = []
    for kind, idx in (('nohit', i_nh), ('hit', i_hit)):
        cand = cKDTree(tc).query_ball_point(xy[idx], r=args.pad_r + 9.0)
        for i, cs in zip(idx, cand):
            cs = np.asarray(cs, int)
            if not len(cs):
                continue
            dd = np.array([pa.dist_to_quads(xy[i][None], tiles[c][None])[0] for c in cs])
            for c, dist in zip(cs[dd <= args.pad_r], dd[dd <= args.pad_r]):
                want.append((kind, int(ev[i]), i, int(tpads.feu[c]), int(tpads.channel[c]), float(dist)))
    W = pd.DataFrame(want, columns=['kind', 'eventId', 'ray', 'feu', 'ch', 'dtile'])
    L.append(f'tracks read: no-hit {len(i_nh):,}, hit {len(i_hit):,}; (track, pad) pairs {len(W):,} '
             f'(mean {W.groupby("ray").size().mean():.2f} pads per track within {args.pad_r:g} mm)')

    # ---- read the decoded waveforms --------------------------------------------
    dec = pw.list_decoded_files(cfg.decoded_root_dir)
    feus = sorted(W.feu.unique())
    present = {f: set() for f in feus}
    rows, stacks = [], []
    for _, fr in dec.iterrows():
        feu = int(fr.feu)
        if feu not in feus:
            continue
        ped = pw.pedestal_from_data(fr.path)
        rms = pw.noise_rms_from_data(fr.path, ped)
        wf_ = W[W.feu == feu]
        by_ev = {e: grp for e, grp in wf_.groupby('eventId')}
        # file-order list of events for the NOISE control (another event, same channels)
        for evs, ts, ft, wf in pw.iter_decoded_events(fr.path, step=2000):
            present[feu].update(int(e) for e in evs)
            sel = np.flatnonzero(np.isin(evs, list(by_ev)))
            if not len(sel):
                continue
            w = pw.subtract_common_noise(wf - ped[None, :, None])
            for k in sel:
                e = int(evs[k])
                for r in by_ev[e].itertuples():
                    rows.append((r.kind, e, r.ray, feu, r.ch, r.dtile, rms[r.ch]))
                    stacks.append(w[k, r.ch].copy())
                    if r.kind == 'nohit':          # same channel, unrelated event in this chunk
                        k2 = (k + args.shift) % len(evs)
                        rows.append(('noise', e, r.ray, feu, r.ch, r.dtile, rms[r.ch]))
                        stacks.append(w[k2, r.ch].copy())
    R = pd.DataFrame(rows, columns=['kind', 'eventId', 'ray', 'feu', 'ch', 'dtile', 'rms'])
    WF = np.asarray(stacks, dtype=np.float32)
    nh_ev = set(ev[i_nh].tolist())
    miss = {f: len(nh_ev - present[f]) for f in feus}
    L.append(f'no-hit events missing from the decoded data: ' + ', '.join(f'FEU {f}: {miss[f]}' for f in feus)
             + f' (of {len(nh_ev)})')

    # ---- signal window and matched filter from HIT pulses ----------------------
    big = (R.kind == 'hit').to_numpy() & (WF.max(1) > 20 * R.rms.to_numpy())
    pk = WF[big].argmax(1)
    s0, s1 = int(np.percentile(pk, 2)) - 2, int(np.percentile(pk, 98)) + 2
    s0, s1 = max(s0, 0), min(s1, WF.shape[1] - 1)
    tmpl = np.median(WF[big] / WF[big].max(1, keepdims=True), 0)
    tw = tmpl[s0:s1 + 1]
    L.append(f'signal window samples {s0}-{s1} (peak sample of big hit pulses: median {np.median(pk):.0f}); '
             f'template from {big.sum():,} hit pulses')
    win = WF[:, s0:s1 + 1]
    R['peak'] = win.max(1)
    R['sig'] = R.peak / R.rms.clip(lower=1)
    R['mf'] = win @ tw / (tw @ tw)
    noise_mf = R[R.kind == 'noise']
    mf_sigma = 1.4826 * np.median(np.abs(noise_mf.mf - noise_mf.mf.median()))
    R['mf_sig'] = R.mf / mf_sigma
    R['theta'] = theta[R.ray.to_numpy()]
    R.to_csv(f'{out}/pad_waveform_features.csv', index=False)
    L.append(f'matched-filter noise sigma {mf_sigma:.1f} ADC (from the NOISE control); '
             f'median channel RMS {R.rms.median():.1f} ADC')

    # per track: the pad with the largest matched-filter significance
    best = R.loc[R.groupby(['kind', 'ray']).mf_sig.idxmax()].reset_index(drop=True)
    L.append('\nper track, best pad within %g mm: fraction with matched-filter significance above k' % args.pad_r)
    L.append('class     theta     n      >2      >3      >4      >5   median mf [ADC]')
    tab = []
    for kind in ('noise', 'nohit', 'hit'):
        for lo, hi in TH_BINS + [(0, 90)]:
            b = best[(best.kind == kind) & (best.theta >= lo) & (best.theta < hi)]
            if not len(b):
                continue
            fr_ = [100 * (b.mf_sig > k).mean() for k in (2, 3, 4, 5)]
            tab.append(dict(kind=kind, th_lo=lo, th_hi=hi, n=len(b), f2=fr_[0], f3=fr_[1], f4=fr_[2], f5=fr_[3],
                            mf_med=b.mf.median()))
            L.append(f'{kind:7s} {lo:3d}-{hi:<3d} {len(b):6d} ' + ' '.join(f'{v:6.1f}%' for v in fr_)
                     + f'   {b.mf.median():8.1f}')
    pd.DataFrame(tab).to_csv(f'{out}/best_pad_significance.csv', index=False)

    # excess of no-hit over noise = tracks with a real sub-threshold pulse
    for lo, hi in TH_BINS + [(0, 90)]:
        a = best[(best.kind == 'nohit') & (best.theta >= lo) & (best.theta < hi)]
        b = best[(best.kind == 'noise') & (best.theta >= lo) & (best.theta < hi)]
        if len(a) and len(b):
            ex3 = (a.mf_sig > 3).mean() - (b.mf_sig > 3).mean()
            L.append(f'  theta {lo}-{hi}: no-hit tracks with a > 3 sigma pulse beyond noise: {100*ex3:.1f} % '
                     f'(n={len(a)})')

    # stacked averages (pad nearest the track = crossed pad, dtile == 0)
    L.append('\nstacked waveform on the CROSSED pad (dtile = 0): mean peak in the window [ADC]')
    fig, axs = plt.subplots(1, 3, figsize=(17, 4.6))
    hit_mean_peak = None
    xs = np.arange(WF.shape[1])
    for j, (lo, hi) in enumerate(TH_BINS):
        ax = axs[j]
        for kind, col in (('hit', C_HIT), ('nohit', C_NH), ('noise', C_NOISE)):
            k = (R.kind == kind).to_numpy() & (R.dtile == 0).to_numpy() & \
                (R.theta >= lo).to_numpy() & (R.theta < hi).to_numpy()
            if k.sum() < 5:
                continue
            mw = WF[k].mean(0); se = WF[k].std(0) / np.sqrt(k.sum())
            scale = 1.0
            lab = f'{kind} (n={k.sum():,})'
            if kind == 'hit':
                hit_mean_peak = mw[s0:s1 + 1].max(); scale = 0.05
                lab += ' x0.05'
            ax.plot(xs, mw * scale, color=col, lw=2, label=lab)
            if kind != 'hit':
                ax.fill_between(xs, (mw - se) * scale, (mw + se) * scale, color=col, alpha=.2, lw=0)
            pkv = mw[s0:s1 + 1] @ tw / (tw @ tw)
            L.append(f'  theta {lo}-{hi} {kind:6s}: n {k.sum():5d}  matched-filter mean {pkv:8.2f} ADC'
                     + (f'  = {100*pkv/(WF[(R.kind=="hit").to_numpy() & (R.dtile==0).to_numpy() & (R.theta>=lo).to_numpy() & (R.theta<hi).to_numpy()].mean(0)[s0:s1+1] @ tw/(tw@tw)):.2f} % of the hit mean'
                        if kind != 'hit' else ''))
        ax.axvspan(s0, s1, color=GRID, alpha=.5, lw=0)
        ax.set(title=f'crossed pad, tracks {lo}-{hi} deg', xlabel='sample (60 ns)', ylabel='mean ADC (CNS)')
        ax.grid(color=GRID); ax.legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(f'{out}/stacked_waveforms.png', dpi=140); plt.close(fig)

    # significance distributions
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.6))
    bins = np.arange(-6, 15.01, 0.5)
    for kind, col in (('noise', C_NOISE), ('nohit', C_NH), ('hit', C_HIT)):
        b = best[best.kind == kind]
        axs[0].hist(b.mf_sig.clip(-6, 15), bins=bins, histtype='step', lw=2, color=col, density=True,
                    label=f'{kind} (n={len(b):,})')
    axs[0].set(xlabel='best pad within %g mm: matched-filter significance' % args.pad_r, ylabel='density',
               yscale='log', title='Is there a pulse under the no-hit tracks?')
    axs[0].axvline(5, color=INK2, ls=':', lw=1)
    for (lo, hi), col in zip(TH_BINS, ('#ec4647', '#f2a600', '#257de0')):
        a = best[(best.kind == 'nohit') & (best.theta >= lo) & (best.theta < hi)]
        axs[1].hist(a.mf_sig.clip(-6, 15), bins=bins, histtype='step', lw=2, color=col, density=True,
                    label=f'no-hit {lo}-{hi} deg (n={len(a)})')
    b = best[best.kind == 'noise']
    axs[1].hist(b.mf_sig.clip(-6, 15), bins=bins, histtype='stepfilled', color=C_NOISE, alpha=.3, density=True,
                label='noise control')
    axs[1].set(xlabel='matched-filter significance', ylabel='density', yscale='log', title='No-hit tracks by angle')
    for ax in axs:
        ax.grid(color=GRID); ax.legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(f'{out}/nohit_significance.png', dpi=140); plt.close(fig)

    with open(f'{out}/summary.txt', 'w') as f:
        f.write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
