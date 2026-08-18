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
            rotation=90, ha='right', va='center', fontsize=8.5,
            color=st.C_RED)
    ax.set_xlabel('track slope dx/dz [mrad]')
    ax.set_ylabel('coincidences')
    ax.legend(loc='upper left', fontsize=8.5)
    ax.set_title('Angle: the pad quantisation is 16$\\times$ the beam '
                 'divergence,\nso the P2-only track cannot measure direction',
                 fontsize=11)

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
                 xytext=(6, 8), textcoords='offset points', fontsize=9,
                 color=st.C_RED)
    ax2.set_xlabel('P2-only crossing point $-$ uRWELL track, '
                   'at the middle plane [mm]')
    ax2.set_ylabel('coincidences')
    ax2.legend(loc='upper left', fontsize=9)
    ax2.set_title('Position: the three planes localise the track to the pad\n'
                  'quantisation — which is what 2-of-3 tagging needs',
                  fontsize=11)

    fig.suptitle('The P2 telescope as a standalone tracker — '
                 f'{len(d)} five-detector coincidences, '
                 'highstat_eff_1/beam_commissioning_00\n'
                 'P2-only line through three cluster centroids vs the uRWELL '
                 'reference track, same pad frames, no external trigger',
                 fontweight='bold')
    st.finish(fig, f'{OUT}/p2_standalone_tracking.png')

    print(f'{len(d)} coincidences, {100 * (d.nmulti == 0).mean():.1f} % '
          f'single-pad on all three planes')
    for a, lab in ((single.dx_at_mid_mm, 'dx@mid'), (single.dy_at_mid_mm, 'dy@mid')):
        print(f'  {lab}: median {np.median(a):+.2f} mm, core sigma {core(a):.2f} mm')
    print(f'  beam divergence: {core(d.slope_x_urw) * 1e3:.2f} / '
          f'{core(d.slope_y_urw) * 1e3:.2f} mrad (x / y)')


if __name__ == '__main__':
    main()
