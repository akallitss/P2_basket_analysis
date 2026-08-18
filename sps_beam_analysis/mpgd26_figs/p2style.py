"""Shared matplotlib style for the MPGD26 SPS beam report figures.

Light-committed figures (PNG) on the validated default palette:
station identity keeps a fixed hue everywhere; gas is encoded by
linestyle+marker (secondary encoding, never color-alone vs the station hue).
"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SURFACE = '#fcfcfb'
TEXT = '#0b0b0b'
TEXT2 = '#52514e'
GRID = '#d9d8d4'

# categorical slots (validated order; first three validate all-pairs)
C_IN = '#2a78d6'    # P2_IN  (blue, slot 1)
C_MID = '#eb6834'   # P2_MID (orange, slot 2)
C_OUT = '#1baf7a'   # P2_OUT (aqua, slot 3)
C_YEL = '#eda100'   # slot 4 (uRW front / misc)
C_MAG = '#e87ba4'   # slot 5 (uRW back / misc)
C_GRN = '#008300'   # slot 6
C_VIO = '#4a3aa7'   # slot 7 (mx17 / model curves)
C_RED = '#e34948'   # slot 8 (reserved-ish; avoid for series)

DET_COLOR = {'P2_IN': C_IN, 'P2_MID': C_MID, 'P2_OUT': C_OUT,
             'uRW_front': C_YEL, 'uRW_back': C_MAG, 'mx17_E': C_VIO}

GAS_A = 'Ar/CO2/Iso 93/5/2'
GAS_B = 'Ar/CF4/Iso 88/10/2'
GAS_STYLE = {GAS_A: dict(ls='-', marker='o'),
             GAS_B: dict(ls='--', marker='^')}
GAS_SHORT = {GAS_A: 'Ar/CO₂/iso 93/5/2', GAS_B: 'Ar/CF₄/iso 88/10/2'}

SEQ_CMAP = 'Blues'        # sequential magnitude (efficiency maps)
SEQ_CMAP_T = 'Oranges'    # sequential magnitude (timing maps)

plt.rcParams.update({
    'figure.facecolor': SURFACE, 'axes.facecolor': SURFACE,
    'savefig.facecolor': SURFACE, 'savefig.dpi': 150,
    'text.color': TEXT, 'axes.labelcolor': TEXT2,
    'xtick.color': TEXT2, 'ytick.color': TEXT2,
    'axes.edgecolor': GRID, 'axes.linewidth': 0.8,
    'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6,
    'grid.alpha': 0.6, 'axes.axisbelow': True,
    'axes.spines.top': False, 'axes.spines.right': False,
    'font.size': 11, 'axes.titlesize': 12, 'axes.titleweight': 'bold',
    'legend.frameon': False, 'legend.fontsize': 9.5,
    'lines.linewidth': 2.0, 'lines.markersize': 5.5,
})


def direct_label(ax, x, y, text, color, dx=6, dy=0, fontsize=9.5):
    ax.annotate(text, xy=(x, y), xytext=(dx, dy), textcoords='offset points',
                fontsize=fontsize, color=color, va='center', fontweight='bold')


def finish(fig, path, title=None, tight=True):
    """Write the figure as PNG *and* PDF -- slides want the vector one."""
    if tight:
        fig.tight_layout()
    for p in (path, os.path.splitext(path)[0] + '.pdf'):
        fig.savefig(p, bbox_inches='tight')
    plt.close(fig)
    print('wrote', path, '(+pdf)')
