#!/usr/bin/env python3
"""A five-detector coincidence, drawn as the detector actually sees it.

The figure this replaces put the three stations side by side as three nearly
identical pad maps -- three panels that all said "one pad fired somewhere over
there" -- plus an x-vs-z panel in which every track was a flat horizontal line
because the axis spanned 200 mm for a 2 mm effect.  Nothing in it showed the
one thing an event display exists to show: that an independent reference track
and three independent pad hits describe the same muon.

So: one hero event, four views.

  (a) face-on   the whole pad plane, the three fired pads on it, and where the
                uRWELL reference says the muon crossed -- with a zoom, because
                at full scale a 12 mm pad is a speck
  (b) x vs z    the telescope from the side, reference track drawn across the
                full 1370 mm uRWELL lever arm, pad hits at their true z and
                true 12 mm width
  (c) y vs z    the same, other projection
  (d) all 25    the per-station reference-to-pad distance for every extracted
                coincidence, so the hero event is visibly typical

Input: nTof_x17/mpgd26/data/coincident_events.json, written on lxplus by
tools/extract_coincident_event.py over the merged hit file.  Same source as
the 3D renders in the deck, so the numbers agree by construction.
"""
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2style as st  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

from paths import OUT  # noqa: E402

EVENTS = os.environ.get(
    'COINCIDENT_EVENTS_JSON',
    os.path.expanduser('~/Documents/PostDocSaclay/nTof_x17/mpgd26/data/'
                       'coincident_events.json'))
DETS = ['P2_IN', 'P2_MID', 'P2_OUT']
PAD_MM = 12.0
ZOOM_MM = 32.0


def pad_polygons():
    """Every pad of the plane as its true Gerber rectangle (board frame)."""
    cba = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__)))), 'cosmic_bench_analysis')
    sys.path.insert(0, cba)
    import p2_qa_config as qa
    import p2_mapping as pm
    m = pm.load_pad_map(qa.MAP_CSV_PATH)
    a = np.radians(m['pad_angle'].to_numpy())
    w, h = m['pad_w'].to_numpy(), m['pad_h'].to_numpy()
    cx, cy = m['pad_cx'].to_numpy(), m['pad_cy'].to_numpy()
    ca, sa = np.cos(a), np.sin(a)
    return np.stack([np.stack([cx + dx * w * ca - dy * h * sa,
                               cy + dx * w * sa + dy * h * ca], axis=1)
                     for dx, dy in ((-.5, -.5), (.5, -.5), (.5, .5), (-.5, .5))],
                    axis=1)


def ref_line(ev, axis):
    """The reference track IN THE PAD FRAME, as (z, value) to plot.

    The raw `urwell` block is in the uRWELL's own frame (x ~ 70 mm where the
    pads sit at ~ 425 mm); plotting it against pad coordinates would draw two
    unrelated things on one axis.  What is directly comparable is
    `pred_c{x,y}_mm`, the same reference track already carried into each
    station's pad frame by the extractor -- and it is exactly what the quoted
    residual is measured against.  Three points, collinear by construction.
    """
    z = np.array([ev['stations'][d]['z_mm'] for d in DETS], float)
    v = np.array([ev['stations'][d][f'pred_c{axis}_mm'] for d in DETS], float)
    m, b = np.polyfit(z, v, 1)
    return z, v, m, b


def side_view(ax, ev, axis, legend=True):
    """One projection of the telescope, seen from the side.

    The y range is a few pads wide on purpose: the whole point is whether the
    reference line threads the three fired pads, and that is a millimetre
    question.  Drawn over 400 mm of pad plane it would be three dots on a
    straight line and would show nothing.
    """
    zp, vp, m, b = ref_line(ev, axis)
    zs = np.array([0.0, ev['urwell']['z_back_mm']])
    ax.plot(zs, m * zs + b, color=st.C_YEL, lw=3.0, zorder=2,
            label='uRWELL reference track, in the pad frame')
    ax.plot(zp, vp, 'o', color=st.C_YEL, ms=9, mec='white', mew=1.4, zorder=4,
            label='where it crosses each P2 plane')

    for det, c in zip(DETS, (st.C_IN, st.C_MID, st.C_OUT)):
        st_ = ev['stations'][det]
        z, v = st_['z_mm'], st_[f'pad_c{axis}_mm']
        ax.axvline(z, color=st.GRID, lw=1.2, zorder=0)
        ax.add_patch(Rectangle((z - 22, v - PAD_MM / 2), 44, PAD_MM,
                               facecolor=c, edgecolor=c, alpha=0.75, zorder=3))
        ax.annotate(f'{det}\n{st_["amplitude"]:.0f} ADC', xy=(z, v),
                    xytext=(0, 30), textcoords='offset points', ha='center',
                    color=c, fontweight='bold', fontsize=11.5)
    c0 = float(np.mean(vp))
    ax.set_ylim(c0 - 26, c0 + 34)
    ax.set_xlim(-60, ev['urwell']['z_back_mm'] + 60)
    ax.set_xlabel('z along the beam [mm]')
    ax.set_ylabel(f'pad {axis} [mm]')
    if legend:
        ax.legend(loc='lower right', fontsize=11)
    ax.set_title(f'Side view, {axis} vs z — pads at their true 12 mm width, '
                 'axis only 60 mm tall', fontsize=13.5)


def main():
    doc = json.load(open(EVENTS))
    evs = doc['events']
    # the hero: three pads, one per plane, smallest reference-to-pad distance
    hero = min((e for e in evs if e['total_pads'] == 3),
               key=lambda e: e['max_residual_mm'])
    u = hero['urwell']

    fig = plt.figure(figsize=(18.6, 10.6))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.15, 1, 1], hspace=0.30,
                          wspace=0.26)
    ax_face = fig.add_subplot(gs[:, 0])
    ax_x = fig.add_subplot(gs[0, 1:])
    ax_y = fig.add_subplot(gs[1, 1])
    ax_r = fig.add_subplot(gs[1, 2])

    # ---- (a) face-on ------------------------------------------------------ #
    verts = pad_polygons()
    ax_face.add_collection(PolyCollection(verts, facecolors='#f5f5f3',
                                          edgecolors='#e0dfdb', linewidths=0.3,
                                          zorder=1))
    px = hero['stations']['P2_MID']['pred_cx_mm']
    py = hero['stations']['P2_MID']['pred_cy_mm']
    for det, c in zip(DETS, (st.C_IN, st.C_MID, st.C_OUT)):
        s = hero['stations'][det]
        ax_face.add_patch(Rectangle(
            (s['pad_cx_mm'] - PAD_MM / 2, s['pad_cy_mm'] - PAD_MM / 2),
            PAD_MM, PAD_MM, facecolor=c, edgecolor=c, alpha=0.75, zorder=3,
            label=f'{det} pad fired'))
    ax_face.plot(px, py, marker='x', ms=16, mew=3.5, color=st.C_RED, zorder=5,
                 ls='none', label='uRWELL track where it crosses P2_MID')
    # at 600 mm across, a 12 mm pad is a speck -- ring it so the eye finds it
    ax_face.plot(px, py, marker='o', ms=34, mfc='none', mec=st.C_RED, mew=2.0,
                 zorder=4, ls='none')
    ax_face.add_patch(Rectangle((px - ZOOM_MM, py - ZOOM_MM), 2 * ZOOM_MM,
                                2 * ZOOM_MM, facecolor='none',
                                edgecolor=st.TEXT2, lw=1.6, ls='--', zorder=4))
    ax_face.set_aspect('equal')
    ax_face.set_xlabel('pad x [mm]')
    ax_face.set_ylabel('pad y [mm]')
    ax_face.grid(False)
    ax_face.legend(loc='lower left', fontsize=11)
    ax_face.set_title('Face-on: the whole 1280-pad plane,\n'
                      'and where this muon went through it', fontsize=13.5)

    # zoom inset -- at full scale a 12 mm pad is a speck
    ins = ax_face.inset_axes([0.55, 0.54, 0.43, 0.38])
    ins.add_collection(PolyCollection(verts, facecolors='none',
                                      edgecolors='#b9b8b3', linewidths=0.8))
    for det, c in zip(DETS, (st.C_IN, st.C_MID, st.C_OUT)):
        s = hero['stations'][det]
        ins.add_patch(Rectangle(
            (s['pad_cx_mm'] - PAD_MM / 2, s['pad_cy_mm'] - PAD_MM / 2),
            PAD_MM, PAD_MM, facecolor=c, edgecolor=c, alpha=0.55))
    ins.plot(px, py, marker='x', ms=13, mew=3.0, color=st.C_RED, ls='none')
    ins.set_xlim(px - ZOOM_MM, px + ZOOM_MM)
    ins.set_ylim(py - ZOOM_MM, py + ZOOM_MM)
    ins.set_aspect('equal'); ins.grid(False)
    ins.set_xticks([]); ins.set_yticks([])
    same = len({(hero['stations'][d]['pad_cx_mm'],
                 hero['stations'][d]['pad_cy_mm']) for d in DETS}) == 1
    ins.text(0.5, 0.03, 'all three planes fired the SAME pad'
             if same else 'the three fired pads', transform=ins.transAxes,
             ha='center', va='bottom', fontsize=11, color=st.TEXT2,
             fontweight='bold',
             bbox=dict(boxstyle='round,pad=0.25', fc=st.SURFACE, ec='none',
                       alpha=0.9))
    for sp in ins.spines.values():
        sp.set_edgecolor(st.TEXT2); sp.set_linewidth(1.6)
        sp.set_linestyle('--')

    # ---- (b, c) side views ------------------------------------------------ #
    side_view(ax_x, hero, 'x')
    side_view(ax_y, hero, 'y', legend=False)
    ax_y.set_title('Side view, y vs z', fontsize=13.5)

    # ---- (d) every extracted coincidence ---------------------------------- #
    for det, c in zip(DETS, (st.C_IN, st.C_MID, st.C_OUT)):
        r = np.array([e['stations'][det]['residual_mm'] for e in evs])
        ax_r.plot(np.arange(len(r)), r, 'o-', color=c, ms=5, lw=1.3,
                  alpha=0.85, label=f'{det}  (median {np.median(r):.1f} mm)')
    ax_r.axhline(PAD_MM / np.sqrt(12), color=st.C_RED, ls='--', lw=2.0)
    ax_r.text(len(evs) - 0.5, PAD_MM / np.sqrt(12) - 0.14,
              '12 mm/$\\sqrt{12}$ = 3.46 mm', color=st.C_RED, ha='right',
              va='top', fontsize=11.5, fontweight='bold')
    hero_i = evs.index(hero)
    ax_r.axvline(hero_i, color=st.TEXT2, lw=1.6, ls=':')
    ax_r.annotate('the event above', xy=(hero_i, 0.06),
                  xycoords=('data', 'axes fraction'), xytext=(6, 0),
                  textcoords='offset points', ha='left', va='bottom',
                  color=st.TEXT2, fontsize=11.5, fontweight='bold')
    ax_r.set_ylim(0, 5.0)
    ax_r.set_xlabel('extracted coincidence')
    ax_r.set_ylabel('reference $-$ pad centre [mm]')
    ax_r.legend(loc='upper left', fontsize=10.5, ncol=3,
                columnspacing=0.8, handletextpad=0.35)
    ax_r.set_title(f'All {len(evs)} of them — the hero event is typical',
                   fontsize=13.5)

    src = doc['source']
    fig.suptitle('One five-detector coincidence, four ways — '
                 f"{src['run']}/{src['sub_run']}, working point\n"
                 'two uRWELL planes give the reference track; three P2 planes '
                 'fire one pad each, all within '
                 f"{hero['max_residual_mm']:.1f} mm of it — no external trigger",
                 fontweight='bold', fontsize=16)
    st.finish(fig, f'{OUT}/event_display.png', tight=False)
    print(f"hero geventId {hero['geventId']}, max residual "
          f"{hero['max_residual_mm']:.2f} mm")


if __name__ == '__main__':
    main()
