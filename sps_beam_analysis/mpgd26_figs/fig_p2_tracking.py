#!/usr/bin/env python3
"""Roadmap 5b -- can the three P2 planes track on their own?

Input: `track_pairs.csv` from `nTof_x17/mpgd26/tools/extract_track_pairs.py`,
run on lxplus over the merged hit file.  Per five-detector coincidence it
carries the uRWELL reference track and a straight line fitted through the three
P2 cluster centroids, both expressed in the same per-station pad frames, so the
comparison involves no frame composition.

The answer, and it is a design statement rather than a defect: with 12 mm pads
and a 620 mm lever arm the P2 telescope **localises but does not collimate**.
One pad step across that lever arm is 12/620 = 19 mrad, while the H4 muon beam
diverges by only ~1.2 mrad -- so in 73 % of coincidences the same pad fires on
all three planes and the P2-only track is exactly parallel to the beam axis by
construction.  What the three planes do deliver is the crossing point, to
3.3 mm at the middle plane = 12 mm/sqrt(12), which is precisely what the
2-of-3 in-situ tagging scheme needs.
"""

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2style as st  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

from paths import OUT  # noqa: E402

CSV = os.environ.get(
    'TRACK_PAIRS_CSV',
    os.path.expanduser('~/Documents/PostDocSaclay/nTof_x17/mpgd26/data/'
                       'track_pairs.csv'))
PAD_MM = 12.0
LEVER_MM = 620.0          # P2_IN -> P2_OUT


def core(x):
    q = np.percentile(x, [16, 84])
    return (q[1] - q[0]) / 2



# --------------------------------------------------------------------------- #
# The same result, told as geometry instead of as histograms.
#
# The histogram version (above) is correct but it asks the audience to convert
# "core sigma 0.00 mrad" into a picture of the detector.  This one draws the
# picture: three pad rows at their real z, the beam's angular spread against
# the angle one pad step subtends, and the crossing-point residual against the
# pad square it is quantised by.  Same numbers, no conversion required.
# --------------------------------------------------------------------------- #
Z = (320.0, 630.0, 940.0)              # P2_IN / P2_MID / P2_OUT, mm
Z_LABEL = ('P2_IN', 'P2_MID', 'P2_OUT')


def fig_explained(d, single):
    fig, axes = plt.subplots(1, 3, figsize=(19.0, 6.2))
    step = PAD_MM / LEVER_MM * 1e3                      # 19.4 mrad
    div = core(d.slope_x_urw) * 1e3                     # 1.3 mrad

    # ---- A: the geometry, to scale --------------------------------------- #
    ax = axes[0]
    zz = np.array([Z[0] - 60, Z[2] + 60])
    for z, lab, c in zip(Z, Z_LABEL, (st.C_IN, st.C_MID, st.C_OUT)):
        ax.axvline(z, color=c, lw=3.0, alpha=0.85, zorder=2)
        ax.text(z, 20.5, lab, color=c, ha='center', va='bottom',
                fontweight='bold', fontsize=13)
    # the pad columns themselves -- boundaries run the length of the telescope
    for k in range(-2, 3):
        ax.axhline((k + .5) * PAD_MM, color=st.TEXT2, lw=1.1, alpha=0.5,
                   zorder=1)
    ax.axhspan(-PAD_MM / 2, PAD_MM / 2, color=st.C_MID, alpha=0.10, zorder=0)
    ax.text(zz[1], PAD_MM / 2 + 0.6, 'one pad, 12 mm', color=st.C_MID,
            ha='right', va='bottom', fontweight='bold', fontsize=13)

    # the beam: every track within the measured divergence
    zc = np.linspace(*zz, 2)
    for k, sl in enumerate((-div, div)):
        ax.plot(zc, (zc - Z[1]) * sl * 1e-3, color=st.C_YEL, lw=2.0)
    ax.fill_between(zc, (zc - Z[1]) * -div * 1e-3, (zc - Z[1]) * div * 1e-3,
                    color=st.C_YEL, alpha=0.55, zorder=1,
                    label=f'the H4 beam: $\\pm${div:.1f} mrad divergence\n'
                          f'= $\\pm${div * LEVER_MM * 1e-3:.1f} mm over the '
                          'full 620 mm lever arm')
    # the smallest angle P2 can even represent
    for sgn in (-1, 1):
        ax.plot(zc, (zc - Z[1]) * sgn * step * 1e-3, color=st.C_RED, lw=2.2,
                ls='--', zorder=4,
                label=(f'the smallest angle P2 can resolve: one pad step,\n'
                       f'{step:.0f} mrad — {step / div:.0f}$\\times$ the '
                       'beam divergence' if sgn > 0 else None))
    ax.set_xlim(*zz); ax.set_ylim(-21, 25)
    ax.set_xlabel('z along the beam [mm]')
    ax.set_ylabel('x across the pads [mm]')
    ax.legend(loc='lower left', fontsize=11.5)
    ax.set_title('Why the P2-only track has no angle:\n'
                 'the whole beam fits inside one pad column', fontsize=13.5)

    # ---- B: the quantisation, seen directly ------------------------------ #
    ax = axes[1]
    ax.scatter(d.slope_x_urw * 1e3, d.slope_x_p2 * 1e3, s=26, alpha=0.55,
               color=st.C_VIO, edgecolors='none', zorder=3)
    for k in range(-1, 2):
        ax.axhline(k * step, color=st.C_RED, ls=':', lw=1.3, zorder=1)
        if k:
            ax.text(34, k * step, f'{k:+d} pad', color=st.C_RED, fontsize=11,
                    va='bottom', ha='right')
    lim = 36
    ax.plot([-lim, lim], [-lim, lim], color=st.TEXT2, lw=1.6, ls='-',
            alpha=0.7, zorder=2, label='a perfect tracker would sit here')
    ax.axhspan(-step / 2, step / 2, color=st.C_MID, alpha=0.10, zorder=0)
    zero = float((d.slope_x_p2.abs() < 1e-9).mean())
    ax.text(lim - 1, 1.2, f'{100 * zero:.0f} % land '
            'exactly on zero', color=st.C_MID, ha='right', va='bottom',
            fontweight='bold', fontsize=12.5)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel('uRWELL track slope dx/dz [mrad]')
    ax.set_ylabel('P2-only track slope dx/dz [mrad]')
    ax.legend(loc='upper left', fontsize=11.5)
    ax.set_title('The same coincidence, both ways:\n'
                 'P2 answers in pad steps, the beam moves between them',
                 fontsize=13.5)

    # ---- C: what P2 DOES deliver ----------------------------------------- #
    ax = axes[2]
    h = ax.hist2d(single.dx_at_mid_mm, single.dy_at_mid_mm,
                  bins=np.linspace(-13, 13, 27), cmap=st.SEQ_CMAP, zorder=1)
    from matplotlib.patches import Rectangle, Circle
    ax.add_patch(Rectangle((-PAD_MM / 2, -PAD_MM / 2), PAD_MM, PAD_MM,
                           facecolor='none', edgecolor=st.C_RED, lw=2.6,
                           zorder=4, label='one 12 mm pad'))
    sig = PAD_MM / np.sqrt(12)
    ax.add_patch(Circle((0, 0), sig, facecolor='none', edgecolor=st.C_YEL,
                        lw=2.6, ls='--', zorder=4,
                        label=f'12 mm/$\\sqrt{{12}}$ = {sig:.2f} mm — '
                              'the quantisation limit'))
    ax.set_aspect('equal')
    ax.set_xlabel('P2 $-$ uRWELL, x at the middle plane [mm]')
    ax.set_ylabel('P2 $-$ uRWELL, y at the middle plane [mm]')
    ax.legend(loc='lower left', fontsize=10.5, framealpha=0.9,
              frameon=True, facecolor=st.SURFACE)
    ax.set_title('What P2 does deliver: the crossing point,\n'
                 f'to {core(single.dx_at_mid_mm):.1f} mm — the pad box, '
                 'nothing worse', fontsize=13.5)
    fig.colorbar(h[3], ax=ax, fraction=0.046, pad=0.03,
                 label='coincidences')

    fig.suptitle('The P2 telescope as a standalone tracker, told as geometry — '
                 f'{len(d)} five-detector coincidences\n'
                 'localises to the pad, cannot collimate: this is the 12 mm '
                 'pad pitch talking, not a defect',
                 fontweight='bold', fontsize=16)
    st.finish(fig, f'{OUT}/p2_standalone_tracking_explained.png')


def main():
    d = pd.read_csv(CSV)
    d['nmulti'] = ((d.n_pads_in > 1).astype(int) + (d.n_pads_mid > 1).astype(int)
                   + (d.n_pads_out > 1).astype(int))
    single = d[d.nmulti == 0]

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13.6, 5.3))

    # ---- angle: the pad step swamps the beam divergence ------------------ #
    bins = np.linspace(-25, 25, 61)
    ax.hist(d.slope_x_urw * 1e3, bins=bins, histtype='stepfilled',
            color=st.C_YEL, alpha=0.55,
            label=f'uRWELL track slope — the beam itself\n'
                  f'core $\\sigma$ = {core(d.slope_x_urw) * 1e3:.2f} mrad')
    ax.hist(d.slope_x_p2 * 1e3, bins=bins, histtype='step', lw=2.2,
            color=st.C_VIO,
            label=f'P2-only track slope, all {len(d)} coincidences\n'
                  f'core $\\sigma$ = {core(d.slope_x_p2) * 1e3:.2f} mrad '
                  f'({100 * (d.nmulti == 0).mean():.0f} % are exactly zero)')
    step = PAD_MM / LEVER_MM * 1e3
    for k in (-1, 1):
        ax.axvline(k * step, color=st.C_RED, ls='--', lw=1.4)
    ax.text(step - 1.0, 0.62, f'one pad step over the 620 mm\nlever arm = '
            f'{step:.0f} mrad', transform=ax.get_xaxis_transform(),
            rotation=90, ha='right', va='center', fontsize=11,
            color=st.C_RED)
    ax.set_xlabel('track slope dx/dz [mrad]')
    ax.set_ylabel('coincidences')
    ax.legend(loc='upper left', fontsize=11)
    ax.set_title('Angle: the pad quantisation is 16$\\times$ the beam '
                 'divergence,\nso the P2-only track cannot measure direction',
                 fontsize=13.5)

    # ---- position: exactly the pad quantisation -------------------------- #
    b2 = np.linspace(-12, 12, 49)
    for v, c, lab in ((single.dx_at_mid_mm, st.C_IN, 'x'),
                      (single.dy_at_mid_mm, st.C_MID, 'y')):
        ax2.hist(v, bins=b2, histtype='step', lw=2.2, color=c,
                 label=f'{lab}: core $\\sigma$ = {core(v):.2f} mm')
    q = PAD_MM / np.sqrt(12)
    for k in (-1, 1):
        ax2.axvline(k * q, color=st.C_RED, ls='--', lw=1.4)
    ax2.annotate(f'12 mm/$\\sqrt{{12}}$ = {q:.2f} mm', xy=(q, 0),
                 xytext=(6, 8), textcoords='offset points', fontsize=11.5,
                 color=st.C_RED)
    ax2.set_xlabel('P2-only crossing point $-$ uRWELL track, '
                   'at the middle plane [mm]')
    ax2.set_ylabel('coincidences')
    ax2.legend(loc='upper left', fontsize=11.5)
    ax2.set_title('Position: the three planes localise the track to the pad\n'
                  'quantisation — which is what 2-of-3 tagging needs',
                  fontsize=13.5)

    fig.suptitle('The P2 telescope as a standalone tracker — '
                 f'{len(d)} five-detector coincidences, '
                 'highstat_eff_1/beam_commissioning_00\n'
                 'P2-only line through three cluster centroids vs the uRWELL '
                 'reference track, same pad frames, no external trigger',
                 fontweight='bold')
    st.finish(fig, f'{OUT}/p2_standalone_tracking.png')

    fig_explained(d, single)

    print(f'{len(d)} coincidences, {100 * (d.nmulti == 0).mean():.1f} % '
          f'single-pad on all three planes')
    for a, lab in ((single.dx_at_mid_mm, 'dx@mid'), (single.dy_at_mid_mm, 'dy@mid')):
        print(f'  {lab}: median {np.median(a):+.2f} mm, core sigma {core(a):.2f} mm')
    print(f'  beam divergence: {core(d.slope_x_urw) * 1e3:.2f} / '
          f'{core(d.slope_y_urw) * 1e3:.2f} mrad (x / y)')


if __name__ == '__main__':
    main()
