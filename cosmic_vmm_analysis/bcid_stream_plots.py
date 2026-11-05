#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figures and pooled table of bcid_stream.py: BCIDs and read-out-stream positions of the three-chamber straight
tracks, against all hits.

  bcid_tracks.png   between chambers: BCID difference, time difference, stream distance (signed, with the chance
                    reference of unrelated hits 5 us later), track span in the hit list, other hits inside the
                    span, stream order of the three chambers
  bcid_marker.png   relative to the SRS marker: time to marker, offset field, words since marker, word position in
                    the UDP frame (track hits against all over-threshold chamber hits), lateness of the chamber hits
                    in the stream, paddle pulses' stream distance to det2
  bcid_summary.json per run and pooled numbers

Usage: bcid_stream_plots.py OUT_DIR bcid_run181.npz [bcid_run184.npz ...]   (the .json next to each npz is read)
"""
import json
import os
import sys

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

INK, MUTED, GRID = '#1f1f1e', '#6b6a64', '#d8d7d2'      # as event_display.py (not imported: it needs the pad map)


def style(ax):
    for k in ('top', 'right'):
        ax.spines[k].set_visible(False)
    for k in ('left', 'bottom'):
        ax.spines[k].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=7)
    ax.grid(color=GRID, lw=0.5)
    ax.set_axisbelow(True)

DETS = (3, 2, 1)
PAIRS = (((1, 0), 'det2 - det3'), ((2, 0), 'det1 - det3'), ((2, 1), 'det1 - det2'))
COL = ('#2a6fb0', '#c2571a', '#3a8f4a')
PADCOL = {'p0_44': '#7a3fa0', 'p0_45': '#b0902a', 'p0_41': '#2a9aa0'}
PADNAME = {'p0_44': 'paddle 44 (bottom)', 'p0_45': 'paddle 45 (top)', 'p0_41': 'paddle 41 (top)'}


def load(paths):
    P, G, runs = {}, {}, []
    for p in paths:
        z = np.load(p)
        js = json.load(open(p.replace('.npz', '.json')))
        tol = js['params']['tol']
        st = (z['r2'] < tol) & ~z['bad']
        runs.append((js, int(st.sum())))
        for k in z.files:
            if k.startswith('h_'):
                G[k] = G.get(k, 0) + z[k]
            elif k.startswith('edges_'):
                G[k] = z[k]
            elif z[k].shape[:1] == st.shape:
                P.setdefault(k, []).append(z[k][st])
    return {k: np.concatenate(v) for k, v in P.items()}, G, runs


def dbcid(P, i, j):
    return ((P['bcid'][:, i] - P['bcid'][:, j] + 2048) % 4096) - 2048


def q(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return {'n': int(len(x)), 'median': float(np.median(x)), 'p16_p84': np.percentile(x, [16, 84]).tolist(),
            'p05_p95': np.percentile(x, [5, 95]).tolist(), 'p99': float(np.percentile(x, 99))} if len(x) else None


def fig_tracks(P, out, ntr):
    fig, axs = plt.subplots(2, 3, figsize=(16, 8.5), constrained_layout=True)
    ax = axs[0, 0]
    style(ax)
    b = np.arange(-8.5, 9)
    for c, ((i, j), lab) in zip(COL, PAIRS):
        ax.hist(dbcid(P, i, j), bins=b, histtype='step', lw=1.8, color=c, label=lab)
    ax.set_xlabel('BCID difference of the first cluster hits [BC = 22.5 ns]')
    ax.set_ylabel('tracks')
    ax.set_title(f'BCID difference between chambers ({ntr} straight tracks)', fontsize=9)
    ax.legend(fontsize=8, frameon=False)
    ax = axs[0, 1]
    style(ax)
    b = np.arange(-160, 161, 8)
    for c, ((i, j), lab) in zip(COL, PAIRS):
        d = P['t'][:, i] - P['t'][:, j]
        ax.hist(d, bins=b, histtype='step', lw=1.8, color=c, label=f'{lab}: sd {np.std(d):.0f} ns')
    ax.set_xlabel('time difference (BCID + TDC) [ns]')
    ax.set_title('time difference between chambers', fontsize=9)
    ax.legend(fontsize=8, frameon=False)
    ax = axs[0, 2]
    style(ax)
    b = np.arange(-12.5, 13)
    for c, ((i, j), lab) in zip(COL, PAIRS):
        ax.hist(P['seq'][:, i] - P['seq'][:, j], bins=b, histtype='step', lw=1.8, color=c, label=f'{lab} (track)')
        ax.hist(P['sb_seq'][:, i] - P['sb_seq'][:, j], bins=b, histtype='step', lw=1, ls=':', color=c)
    ax.plot([], [], ls=':', color=MUTED, label='unrelated hits (5 us later)')
    ax.set_yscale('log')
    ax.set_xlabel('stream distance of the first cluster hits [hit words]')
    ax.set_title('where in the hit list: distance between chambers', fontsize=9)
    ax.legend(fontsize=7.5, frameon=False)
    ax = axs[1, 0]
    style(ax)
    b = np.arange(-0.5, 41)
    ax.hist(np.clip(P['span'], 0, 40), bins=b, histtype='step', lw=2, color=INK, label='track: all cluster hits')
    ax.hist(np.clip(P['sb_span'], 0, 40), bins=b, histtype='step', lw=1.2, ls=':', color=MUTED,
            label='first hits of the three chambers 5 us later')
    ax.set_yscale('log')
    ax.set_xlabel('span in the hit list, last - first [hit words] (40 = overflow)')
    ax.set_title('how spread the track hits are in the list', fontsize=9)
    ax.legend(fontsize=8, frameon=False)
    ax = axs[1, 1]
    style(ax)
    b = np.arange(-0.5, 21)
    for k, lab, c in (('foreign_ot', 'other over-threshold chamber hits', INK), ('foreign_nb', 'neighbour hits (flag 0)', COL[0]),
                      ('foreign_pad', 'hybrid-0 (paddle) hits', COL[1])):
        x = P[k]
        ax.hist(np.clip(x, 0, 20), bins=b, histtype='step', lw=1.6, color=c, label=f'{lab}: {100 * np.mean(x > 0):.1f} % of tracks > 0')
    ax.set_yscale('log')
    ax.set_xlabel('hits inside the track span that are not track hits (20 = overflow)')
    ax.set_title('non-track hits between the track hits', fontsize=9)
    ax.legend(fontsize=8, frameon=False)
    ax = axs[1, 2]
    style(ax)
    keys, cnt = np.unique([''.join(str(DETS[x]) for x in r) for r in P['order']], return_counts=True)
    o = np.argsort(-cnt)
    ax.bar(range(len(o)), cnt[o] / cnt.sum(), color=INK)
    ax.set_xticks(range(len(o)))
    ax.set_xticklabels(['-'.join('det' + c for c in keys[x]) for x in o], rotation=30, fontsize=7.5)
    ax.set_ylabel('fraction of tracks')
    ax.set_title('order of the chambers in the stream (first hit per chamber)', fontsize=9)
    fig.savefig(os.path.join(out, 'bcid_tracks.png'), dpi=100)
    plt.close(fig)


def fig_marker(P, G, out):
    fig, axs = plt.subplots(2, 3, figsize=(16, 8.5), constrained_layout=True)

    def comp(ax, key, edges, xs, lab, scale=1.0, logy=False, rebin=1):
        style(ax)
        h = G[f'h_{key}_ot'].astype(float)
        if rebin > 1:
            k = len(h) // rebin * rebin
            h = h[:k].reshape(-1, rebin).sum(1)
            edges = edges[:k + 1:rebin]
        ctr = (edges[1:] + edges[:-1]) / 2 * scale
        ht = np.histogram(xs, bins=edges)[0].astype(float)
        ax.step(ctr, h / max(h.sum(), 1), where='mid', color=MUTED, lw=1.2, label=f'all over-threshold chamber hits ({int(h.sum())})')
        ax.step(ctr, ht / max(ht.sum(), 1), where='mid', color=INK, lw=1.8, label=f'track hits, first per chamber ({len(xs)})')
        ax.set_xlabel(lab)
        ax.set_ylabel('fraction')
        if logy:
            ax.set_yscale('log')
        ax.legend(fontsize=7.5, frameon=False)
    E = {k[6:]: G[k] for k in G if k.startswith('edges_')}
    comp(axs[0, 0], 'chip', E['chip'], P['chip'].ravel(), 'hit time - its SRS marker [us]', scale=1e-3, rebin=8)
    axs[0, 0].set_xlim(-200, 1600)
    axs[0, 0].set_title('time to the SRS marker (marker every 65536 x 22.5 ns = 1.47 ms)', fontsize=9)
    comp(axs[0, 1], 'offs', E['offs'], P['offs'].ravel(), 'offset field (4096-BC overflows since the marker)', logy=True)
    axs[0, 1].set_title('offset field', fontsize=9)
    comp(axs[0, 2], 'wsm', E['wsm'], P['wsm'].ravel(), 'words since the marker of the same VMM', logy=True)
    axs[0, 2].set_xlim(0, 400)
    axs[0, 2].set_title('position after the marker word in the stream', fontsize=9)
    comp(axs[1, 0], 'w', E['w'], P['w'].ravel(), 'word position in the UDP frame (1491 words)', rebin=5)
    axs[1, 0].set_title('position in the UDP frame', fontsize=9)
    ax = axs[1, 1]
    style(ax)
    e = E['late']
    for key, lab, c in (('ot', 'over-threshold chamber hits', INK), ('pad', 'hybrid-0 (paddle) hits', COL[1])):
        h = G[f'h_late_{key}'].astype(float)
        frac = h[::-1].cumsum()[::-1] / h.sum()
        ax.step(e[:-1] * 1e-3, frac, where='post', color=c, lw=1.6, label=lab)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(0.5, 100)
    ax.set_ylim(1e-4, 1.1)
    ax.set_xlabel('lateness = latest time already read - own time [us]')
    ax.set_ylabel('fraction of hits later than x')
    ax.set_title('how far out of time order the stream is', fontsize=9)
    ax.legend(fontsize=8, frameon=False)
    ax = axs[1, 2]
    style(ax)
    b = np.arange(-10.5, 11)
    for k in ('p0_44', 'p0_45', 'p0_41'):
        if k + '_hit' in P:
            m = P[k + '_hit'].astype(bool)
            ax.hist(P[k + '_dseq'][m], bins=b, histtype='step', lw=1.6, color=PADCOL[k], label=f'{PADNAME[k]} ({m.sum()})')
    ax.set_xlabel('stream distance paddle pulse - det2 first hit [hit words]')
    ax.set_title('the paddle pulses of the track in the list', fontsize=9)
    ax.legend(fontsize=8, frameon=False)
    fig.savefig(os.path.join(out, 'bcid_marker.png'), dpi=100)
    plt.close(fig)


def main():
    out, paths = sys.argv[1], sys.argv[2:]
    os.makedirs(out, exist_ok=True)
    P, G, runs = load(paths)
    ntr = len(P['t2'])
    fig_tracks(P, out, ntr)
    fig_marker(P, G, out)
    S = {'n_tracks': ntr, 'runs': [{'title': js['title'], 'n_straight_good': n, 'span_s': js['span_s'],
                                    'tracks_per_hour': n / js['span_s'] * 3600,
                                    'frac_all_in_one_frame': js['tracks']['frac_all_in_one_frame'],
                                    'span_median': js['tracks']['span']['median'],
                                    'lateness_p99_ns': js['chamber_hits_lateness_ns']['p50_p90_p99_p999'][2],
                                    'marker_rounds': js['marker_rounds']} for js, n in runs]}
    for (i, j), lab in PAIRS:
        S[lab] = {'dbcid': q(dbcid(P, i, j)), 'frac_abs_dbcid_le_1': float(np.mean(np.abs(dbcid(P, i, j)) <= 1)),
                  'frac_abs_dbcid_le_3': float(np.mean(np.abs(dbcid(P, i, j)) <= 3)),
                  'dt_ns': q(P['t'][:, i] - P['t'][:, j]), 'frac_same_marker': float(np.mean(P['srs'][:, i] == P['srs'][:, j])),
                  'dseq': q(P['seq'][:, i] - P['seq'][:, j]), 'dseq_chance': q(P['sb_seq'][:, i] - P['sb_seq'][:, j]),
                  'frac_same_frame': float(np.mean(P['pkt'][:, i] == P['pkt'][:, j]))}
    for c, d in enumerate(DETS):
        S[f'det{d}'] = {'chip_time_us': q(P['chip'][:, c] * 1e-3), 'offset': q(P['offs'][:, c]),
                        'words_since_marker': q(P['wsm'][:, c]), 'word_in_frame': q(P['w'][:, c]),
                        'n_hits': q(P['n'][:, c]), 'n_vmm': q(P['nvmm'][:, c]), 'frac_negative_offset': float(np.mean(P['offs'][:, c] < 0))}
    S['span'] = q(P['span'])
    S['span_chance'] = q(P['sb_span'])
    S['frac_span_le_track_hits'] = float(np.mean(P['span'] + 1 == P['n'].sum(1)))
    S['frac_all_in_one_frame'] = float(np.mean(P['pkt_span'] == 0))
    for k in ('foreign_ot', 'foreign_nb', 'foreign_pad', 'own_nb'):
        S[k] = {'frac_gt0': float(np.mean(P[k] > 0)), 'mean': float(P[k].mean()), 'p99': float(np.percentile(P[k], 99))}
    keys, cnt = np.unique([''.join(str(DETS[x]) for x in r) for r in P['order']], return_counts=True)
    S['stream_order'] = {k: float(v / cnt.sum()) for k, v in sorted(zip(keys, cnt), key=lambda kv: -kv[1])}
    for k in ('p0_44', 'p0_45', 'p0_41'):
        if k + '_hit' in P:
            m = P[k + '_hit'].astype(bool)
            S[PADNAME[k]] = {'n': int(m.sum()), 'dseq_vs_det2': q(P[k + '_dseq'][m]),
                             'dbcid_det2_minus_paddle': q(((P['bcid'][m, 1] - P[k + '_bcid'][m] + 2048) % 4096) - 2048),
                             'frac_same_marker_as_det2': float(np.mean(P[k + '_srs'][m] == P['srs'][m, 1]))}
    h = G['h_late_ot'].astype(float)
    e = G['edges_late']
    S['lateness_all_chamber_hits'] = {f'frac_later_than_{x}us': float(h[e[:-1] >= x * 1e3].sum() / h.sum()) for x in (1, 5, 20, 50, 90)}
    json.dump(S, open(os.path.join(out, 'bcid_summary.json'), 'w'), indent=1)
    print(json.dumps({k: v for k, v in S.items() if k != 'runs'}, indent=1)[:6000])


if __name__ == '__main__':
    sys.exit(main())
