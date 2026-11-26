#!/usr/bin/env python3
"""29_crosstalk_illustration.py -- where the shorted / dead / noisy channels sit:
on the DETECTOR (pads, far apart) vs on the CONNECTOR footprint (pins, adjacent).

Usage:  python3 29_crosstalk_illustration.py <run_key>

Uses cfg.SHORTED_CHANNELS (grouped by connectivity on the Dream channel number,
links n<->n+1 and n<->n+2), cfg.DEAD_CHANNELS and the Gerber pad map's
connector footprint (pin_x, pin_y). Waveform panels come from decoded_root;
the coupling panel from stage 26's xtalk_rows.csv when present.
Output: 29_crosstalk_illustration/crosstalk_illustration.png (+ pins csv)
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2_qa_config as qa      # noqa: E402
import p2_mapping as pmap      # noqa: E402
import p2_waveforms as pw      # noqa: E402

PAL = ['#257de0', '#f47230', '#00b894', '#f2a600', '#ee7fac', '#008a17', '#6a4ce2', '#ec4647']
INK, INK2, GRID, PADC = '#1f1f1f', '#5c5c5c', '#e6e6e6', '#d9d9d9'


def groups_of(shorted):
    """Connected groups of shorted (feu, ch): link channels 1 or 2 apart."""
    s = sorted(set(shorted))
    par = {x: x for x in s}

    def f(x):
        while par[x] != x:
            x = par[x]
        return x
    for a in s:
        for d in (1, 2):
            b = (a[0], a[1] + d)
            if b in par:
                par[f(a)] = f(b)
    out = {}
    for x in s:
        out.setdefault(f(x), []).append(x)
    return [sorted(v) for v in sorted(out.values())]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run_key', nargs='?', default=qa.DEFAULT_RUN)
    args = ap.parse_args()
    cfg = qa.get_config(args.run_key)
    out_dir = cfg.out_dir('29_crosstalk_illustration')

    # FULL table (nothing dropped except unplugged connectors) + pin geometry
    ct = pmap.build_channel_table(cfg.run_config_path, cfg.MAP_CSV_PATH, det_type=cfg.DET_TYPE,
                                  det_name=cfg.DET_NAME, drop_connectors=cfg.DEAD_CONNECTORS,
                                  strategy_overrides=cfg.STRATEGY_OVERRIDES)
    mp = pd.read_csv(cfg.MAP_CSV_PATH).set_index('channel_id')
    t = ct.join(mp[['connector_number', 'connector_pin', 'pin_name', 'pin_x', 'pin_y', 'mec8_pin']],
                on='channel_id')
    tk = t.set_index(['feu', 'channel'])
    groups = [list(g) for g in cfg.SHORTED_GROUPS] or groups_of(cfg.SHORTED_CHANNELS)
    shorted = set(cfg.SHORTED_CHANNELS)
    dead = [x for x in cfg.DEAD_CHANNELS if x not in shorted]

    # classify dead: dead-pin rows (long runs of one parity) vs isolated/noisy
    rows = []
    for gi, g in enumerate(groups):
        for x in g:
            r = tk.loc[x]
            rows.append(dict(kind=f'short {gi+1}', feu=x[0], ch=x[1], **{k: r[k] for k in
                        ('connector_N', 'half', 'strip', 'channel_id', 'connector_number',
                         'connector_pin', 'pin_name', 'mec8_pin', 'pin_x', 'pin_y', 'pad_cx', 'pad_cy')}))
    for x in dead:
        r = tk.loc[x]
        rows.append(dict(kind='dead', feu=x[0], ch=x[1], **{k: r[k] for k in
                    ('connector_N', 'half', 'strip', 'channel_id', 'connector_number',
                     'connector_pin', 'pin_name', 'mec8_pin', 'pin_x', 'pin_y', 'pad_cx', 'pad_cy')}))
    pins = pd.DataFrame(rows)
    pins.to_csv(f'{out_dir}/shorted_dead_pins.csv', index=False)

    plt.rcParams.update({'font.size': 9.5, 'axes.edgecolor': INK2, 'xtick.color': INK2,
                         'ytick.color': INK2, 'axes.labelcolor': INK})
    fig = plt.figure(figsize=(16, 15))
    gs = fig.add_gridspec(3, 4, height_ratios=[1.35, 1.0, 0.85], hspace=0.38, wspace=0.28)

    # ---------- A: detector pads ----------
    axA = fig.add_subplot(gs[0, :2])
    if pmap.has_tile_geometry(ct):
        tp, verts = pmap.pad_tiles(ct)
        axA.add_collection(PolyCollection(verts, facecolor=PADC, edgecolor='white', lw=0.3))
    full = ct[ct['mapped']].drop_duplicates('channel_id')
    axA.scatter(full.pad_cx, full.pad_cy, s=0)  # autoscale
    dd = pins[pins.kind == 'dead']
    axA.scatter(dd.pad_cx, dd.pad_cy, marker='x', s=28, color=INK, lw=1.3, label='dead / dead-but-noisy pad', zorder=3)
    for gi, g in enumerate(groups):
        s = pins[pins.kind == f'short {gi+1}']
        c = PAL[gi % len(PAL)]
        for i in range(len(s)):
            for j in range(i + 1, len(s)):
                axA.plot([s.pad_cx.iloc[i], s.pad_cx.iloc[j]], [s.pad_cy.iloc[i], s.pad_cy.iloc[j]],
                         color=c, lw=1.2, alpha=0.8, zorder=4)
        axA.scatter(s.pad_cx, s.pad_cy, s=55, color=c, edgecolor='white', lw=1, zorder=5,
                    label=f'short {gi+1}: FEU{g[0][0]} ch ' + '/'.join(str(x[1]) for x in g))
    axA.set_aspect('equal')
    axA.set_title('A  On the DETECTOR: pads of each shorted group are 46-354 mm apart',
                  loc='left', fontweight='bold', color=INK)
    axA.set_xlabel('pad frame x [mm]'); axA.set_ylabel('pad frame y [mm]')
    axA.legend(fontsize=7.5, frameon=False, loc='upper left', ncol=1)
    for s_ in ('top', 'right'):
        axA.spines[s_].set_visible(False)

    # ---------- B: connector footprints (zooms) ----------
    conns = sorted(pins.connector_number.unique())
    ncol_b = 3
    sub = gs[0, 2:].subgridspec((len(conns) + ncol_b - 1) // ncol_b, ncol_b, hspace=0.5, wspace=0.35)
    for k, cn in enumerate(conns):
        ax = fig.add_subplot(sub[k // ncol_b, k % ncol_b])
        allp = mp[mp.connector_number == cn]
        here = pins[pins.connector_number == cn]
        # zoom on the involved pins (+ margin)
        x0, x1 = here.pin_x.min() - 3, here.pin_x.max() + 3
        y0, y1 = here.pin_y.min() - 3, here.pin_y.max() + 3
        sel = allp[(allp.pin_x.between(x0, x1)) & (allp.pin_y.between(y0, y1))]
        ax.scatter(sel.pin_x, sel.pin_y, s=14, marker='s', color=PADC, zorder=1)
        d_ = here[here.kind == 'dead']
        ax.scatter(d_.pin_x, d_.pin_y, s=26, marker='x', color=INK, lw=1.2, zorder=3)
        for gi in range(len(groups)):
            s = here[here.kind == f'short {gi+1}']
            if len(s):
                ax.scatter(s.pin_x, s.pin_y, s=42, marker='s', color=PAL[gi % len(PAL)],
                           edgecolor='white', lw=0.8, zorder=4)
                for _, r in s.iterrows():
                    ax.annotate(f'{int(r.ch)}', (r.pin_x, r.pin_y), xytext=(3, 3),
                                textcoords='offset points', fontsize=7, color=INK)
        cN = int(here.connector_N.iloc[0]) if len(here) else cn + 1
        ax.set_title(f'B{k+1}  board connector {cN} (footprint)', loc='left', fontsize=9, color=INK)
        ax.set_aspect('equal'); ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
        ax.tick_params(labelsize=7)
        for s_ in ('top', 'right'):
            ax.spines[s_].set_visible(False)

    # ---------- C: waveforms ----------
    dec = pw.list_decoded_files(cfg.decoded_root_dir)
    axC1 = fig.add_subplot(gs[1, :2]); axC2 = fig.add_subplot(gs[1, 2:])
    tns = np.arange(32) * 60
    g0 = groups[0]
    feu0 = g0[0][0]
    f0 = dec[dec.feu == feu0].iloc[0].path
    ped = pw.pedestal_from_data(f0)
    rms = pw.noise_rms_from_data(f0, ped)
    shown_short = shown_norm = False
    skip = set(c for f_, c in cfg.DEAD_CHANNELS if f_ == feu0)
    for ev, ts, ft, wf in pw.iter_decoded_events(f0, step=2000):
        w = pw.subtract_common_noise(wf - ped[None, :, None])
        if not shown_short:
            a = w[:, g0[0][1]].max(1)
            k = np.where((a > 1200) & (a < 3500))[0]
            if len(k):
                e = k[0]
                for i, x in enumerate(g0):
                    axC1.plot(tns, w[e, x[1]], color=PAL[0] if i == 0 else PAL[(i + 2) % 8], lw=2.4 - 0.4 * i,
                              ls='-' if i == 0 else '--', label=f'FEU{x[0]} ch {x[1]}  (pad {int(tk.loc[x].channel_id)})')
                axC1.set_title(f'C1  Shorted group 1, one event: identical waveforms on pads up to '
                               f'{np.ptp(pins[pins.kind=="short 1"].pad_cx):.0f} mm apart',
                               loc='left', fontweight='bold', color=INK, fontsize=10)
                shown_short = True
        if not shown_norm:
            amax = w.max(2)
            for e in range(len(ev)):
                cand = [c for c in np.argsort(amax[e])[::-1][:3]
                        if c not in skip and c + 1 not in skip and c - 1 not in skip and 2500 < amax[e, c] < 3600
                        and (c + 17) // 64 == c // 64 and (c - 1) // 64 == c // 64
                        and c + 2 not in skip]
                if cand and amax[e, cand[0] + 1] > 5 * rms[cand[0] + 1]:
                    c = cand[0]
                    axC2.plot(tns, w[e, c] / 50, color=PAL[0], lw=2.4, label=f'source ch {c}  (scaled x 1/50)')
                    for o, col in ((-1, PAL[1]), (1, PAL[2]), (2, PAL[4]), (17, '#9c9c9c')):
                        axC2.plot(tns, w[e, c + o], color=col, lw=1.6, ls='-' if o != 17 else ':',
                                  label=f'ch {c+o} ({o:+d}{", control" if o == 17 else ""})')
                    axC2.axhline(5 * rms[c + 1], color=INK2, lw=0.8, ls='--')
                    axC2.text(tns[-1], 5 * rms[c + 1], ' 5 sigma', va='bottom', ha='right', fontsize=8, color=INK2)
                    axC2.set_title(f'C2  Ordinary crosstalk, one large pulse ({amax[e, c]:.0f} ADC): '
                                   f'same-shape ~1% copy on ch +-1/+-2', loc='left', fontweight='bold',
                                   color=INK, fontsize=10)
                    shown_norm = True
                    break
        if shown_short and shown_norm:
            break
    for ax in (axC1, axC2):
        ax.set_xlabel('time [ns]'); ax.set_ylabel('ADC (pedestal + common-noise subtracted)')
        ax.grid(True, color=GRID, lw=0.8); ax.legend(fontsize=8, frameon=False)
        for s_ in ('top', 'right'):
            ax.spines[s_].set_visible(False)

    # ---------- D: coupling vs channel offset ----------
    axD = fig.add_subplot(gs[2, :2])
    xr = os.path.join(cfg.out_dir('26_crosstalk'), 'xtalk_rows.csv')
    if os.path.isfile(xr):
        x = pd.read_csv(xr)
        sh = set(cfg.SHORTED_CHANNELS)
        x = x[[(f, s_) not in sh and (f, s_ + o) not in sh for f, s_, o in zip(x.feu, x.src, x.off)]]
        g = x.groupby(['kind', 'off'])['at_src_peak'].median().reset_index()
        g = g.sort_values('off')
        lab = [f'{o:+d}' for o in g.off]
        cols = [PAL[0] if k == 'nb' else '#9c9c9c' for k in g.kind]
        axD.bar(range(len(g)), 100 * g.at_src_peak, color=cols, width=0.62)
        for i, v in enumerate(100 * g.at_src_peak):
            axD.text(i, v + (0.03 if v >= 0 else -0.03), f'{v:.2f}%', ha='center',
                     va='bottom' if v >= 0 else 'top', fontsize=8.5, color=INK)
        axD.set_xticks(range(len(g))); axD.set_xticklabels(lab)
        axD.axhline(0, color=INK2, lw=0.8)
        axD.set_xlabel('Dream channel offset from the source (grey = same-chip controls)')
        axD.set_ylabel('neighbour / source amplitude [%]')
        axD.set_title('D  Coupling (median, shorted channels excluded): +-1 and +-2 equal, +-3 small',
                      loc='left', fontweight='bold', color=INK, fontsize=10)
        axD.grid(True, axis='y', color=GRID, lw=0.8); axD.set_axisbelow(True)
        for s_ in ('top', 'right'):
            axD.spines[s_].set_visible(False)

    # ---------- E: text box ----------
    axE = fig.add_subplot(gs[2, 2:]); axE.axis('off')
    lines = ['How to read this',
             '- Connector footprint (B): 2 pin rows ~8.3 mm apart, 0.5 mm pitch in a row;',
             '  even Dream channels in one row, odd in the other: ch n+2 sits NEXT TO n,',
             '  ch n+1 FACES it. Numbers = Dream channel, x = dead / noisy pin.',
             '- Shorted group = channels carrying the SAME waveform (amp ratio 1.00 +- 0.04).',
             '  Their pads are far apart on the detector (A, C1) but their connector pins',
             '  are neighbours (B): 0.5 mm in a row (ch n / n+2) -> solder bridge on the',
             '  board connector; facing pins (ch n / n+1, 8.3 mm on the board) are adjacent',
             '  on the electronics side (ribbon/adapter) instead.',
             '- Ordinary crosstalk (C2, D) is a ~0.85 % same-shape copy on ch +-1/+-2 only;',
             '  it crosses threshold only for very large pulses.',
             '- Dead pin row: every even ch 384-446 on FEU 6 = one whole row of one',
             '  half-connector: no signal and almost no noise (open input).',
             '- Dead-but-noisy pads: runs of consecutive pins in one row (FEU 6 odd',
             '  217-227, even 246-254) + FEU 7 254/281: they fire on noise only,',
             '  never on a muon over their pad (floating input).']
    for i, l in enumerate(lines):
        axE.text(0, 1 - i * 0.066, l, fontsize=9.5 if i else 11, fontweight='bold' if i == 0 else None,
                 color=INK if i == 0 else INK2, va='top', family='DejaVu Sans')

    fig.suptitle(f'{cfg.DET_NAME} crosstalk, shorts and dead channels — {cfg.RUN}/{cfg.SUB_RUN}',
                 x=0.01, ha='left', fontsize=13, color=INK)
    fig.savefig(f'{out_dir}/crosstalk_illustration.png', dpi=140, bbox_inches='tight')
    plt.close(fig)
    print(f'  saved {out_dir}/crosstalk_illustration.png')
    print(pins.groupby('kind').size().to_string())


if __name__ == '__main__':
    main()
