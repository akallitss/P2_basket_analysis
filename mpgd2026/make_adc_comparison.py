#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_adc_comparison.py -- the figure behind the framing line:

    "we ran the VMM with one global threshold per chip and never used the
     per-channel trim DACs -- here is what that costs, measured"

Two panels, both from measured quantities on the SAME chamber (P2_OUT) at the
SAME working point (mesh 450 V, drift 750 V, Ar/CO2/iC4H10 93/5/2):

  (a) what DREAM sees -- the full Landau cluster-charge spectrum of the
      highest-statistics run of the campaign (highstat_eff_1 /
      beam_commissioning_00, 25 Jul, 7.88 M clusters), with the window the VMM
      actually recorded marked on it.  The point of the panel is that the VMM's
      recording window opens where the DREAM spectrum is still climbing.

  (b) why a single global threshold costs more than one hard cut would --
      per-pad VMM efficiency against per-pad median pulse height (run_46,
      1 Aug, the best configuration of the campaign).  With no per-channel
      trim the effective threshold moves with every channel's baseline, so the
      pads spread over 0.59-0.94 while DREAM measures 0.96 uniformly.

Caveats carried on the figure itself, because a slide has no caption:
  * DREAM ADC is cluster charge (12-bit, cluster size 1.05-1.08 pads); VMM ADC
    is per-pad (10-bit).  The two are NOT the same unit, so each axis is
    normalised to that readout's own Landau MPV and they never share a raw
    scale.
  * The DREAM and VMM data are 7 days apart -- the pads were recabled between
    them, so the two readouts could never run at once.

Usage:  python3 make_adc_comparison.py [-o OUTDIR] [--vmm-hist FILE.npz]

--vmm-hist takes an npz with `counts` and `edges` of the track-matched VMM
per-pad ADC spectrum of P2_OUT in run_46.  It is not on this machine (it needs
the autopsy parquet from EOS), so without it panel (a) draws the VMM window
from the measured quantiles below and says so on the figure.
"""

import argparse
import json
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ANALYSIS = ('/local/home/ak271430/Documents/PostDocSaclay/data/SPS_Beam_Test/'
            'mpgd26_workspace/products/analysis')
DREAM_SPEC = (f'{ANALYSIS}/P2_OUT/highstat_eff_1/beam_commissioning_00/'
              '20_beam_spectra/scan_row.json')

# --- measured VMM per-pad ADC quantiles, P2_OUT, run_46 -------------------- #
# EFFICIENCY_AUTOPSY.md "What it is", item 1.  In VMM ADC counts.
VMM = dict(lowest=47.0, p05=81.0, mpv=128.0)

# --- measured per-pad efficiency vs per-pad median pulse height ------------ #
# EFFICIENCY_AUTOPSY.md item 7, quintiles over the 69 illuminated P2_OUT pads.
QUINT = [((90, 127), 0.752, 151914), ((127, 139), 0.816, 241640),
         ((139, 168), 0.892, 163012), ((168, 202), 0.920, 174252),
         ((202, 271), 0.921, 106159)]
PAD_BEST = [0.962, 0.957, 0.954]          # best individual pads, same run
PAD_P05, PAD_MED, PAD_P95 = 0.59, 0.894, 0.943

# DREAM P2_OUT at run_46's own working point (mesh 450 / drift 750), from the
# three independent uRWELL-referenced runs -- see V2 in COMPARISON_CLAIMS.md.
DREAM_BAND = (0.9546, 0.9655)
DREAM_C, VMM_C = '#111111', '#2ca02c'

plt.rcParams.update({
    'xtick.labelsize': 16, 'ytick.labelsize': 16,
    'font.size': 14, 'axes.titlesize': 14, 'axes.labelsize': 17,
    'legend.fontsize': 10, 'figure.facecolor': 'white',
    'axes.grid': True, 'grid.alpha': 0.3, 'axes.axisbelow': True,
})

NL = chr(10)


def dream_spectrum():
    d = json.load(open(DREAM_SPEC))
    c = np.asarray(d['spectrum']['counts'], float)
    e = np.asarray(d['spectrum']['edges'], float)
    return c, e, float(d['row']['mpv_adc']), int(d['row']['n_events'])


def panel_a(ax, vmm_hist, label='(a) '):
    c, e, mpv, n = dream_spectrum()
    ctr = 0.5 * (e[:-1] + e[1:])
    x = ctr / mpv                                   # normalise to own MPV
    y = c / c.sum() / (np.diff(e) / mpv)            # density per MPV unit
    ymax = y.max()

    vmpv_ref = VMM['mpv']
    if vmm_hist is not None:
        _m = np.load(vmm_hist)
        if 'meta' in _m.files:
            vmpv_ref = float(json.loads(str(_m['meta']))
                             .get('mpv_adc_dnl_robust') or VMM['mpv'])
    lo, p05 = VMM['lowest'] / vmpv_ref, VMM['p05'] / vmpv_ref

    # the part of the spectrum the discriminator removes, under everything else
    ax.axvspan(0, lo, color='#d62728', alpha=0.13, zorder=0)
    ax.fill_between(x, 0, y, where=(x < lo), step='mid',
                    color='#d62728', alpha=0.35, zorder=1)

    ax.step(x, y, where='mid', color=DREAM_C, lw=2.2, zorder=3,
            label=f'DREAM, cluster charge ({n / 1e6:.2f} M clusters)')

    # How much of the DREAM Landau sits below where the VMM starts recording.
    # This is the cut-off stated as a number rather than left to the eye.
    frac_lost = float(c[x < lo].sum() / c.sum())

    note = ''
    if vmm_hist is not None:
        z = np.load(vmm_hist)
        vc, ve = np.asarray(z['counts'], float), np.asarray(z['edges'], float)
        # The npz is stored at the VMM's native unit-ADC resolution, which is
        # right for the file and far too fine to draw: 1024 bins of 715k
        # entries is a hedge of spikes. Rebin for display only -- the MPV the
        # axis is normalised by stays the unit-bin one.
        vm = json.loads(str(z['meta'])) if 'meta' in z.files else {}
        # Rebin by the DNL period, not by a display-pretty factor: the VMM ADC
        # over-populates every 16th code by 2-3x, and any binning that is not a
        # multiple of 16 leaves that comb standing in the drawn spectrum.
        k = int(vm.get('dnl_period') or 16)
        n = (len(vc) // k) * k
        vc = vc[:n].reshape(-1, k).sum(axis=1)
        ve = ve[:n + 1:k]
        # and normalise by the DNL-robust peak -- the unit-bin mode is a comb
        # tooth (128 for P2_OUT), not the Landau peak (104).
        vmpv = float(vm.get('mpv_adc_dnl_robust') or VMM['mpv'])
        vx = 0.5 * (ve[:-1] + ve[1:]) / vmpv
        vy = vc / vc.sum() / (np.diff(ve) / vmpv)
        ax.step(vx, vy, where='mid', color=VMM_C, lw=2.2, zorder=3,
                label='VMM3a, per-pad ADC (run_46)')
        # both curves must fit: the VMM peak is taller than the DREAM one
        ymax = max(ymax, float(vy[(vx > 0.2) & (vx < 3.0)].max()))
    else:
        note = ('VMM drawn as its measured low edge only;' + NL
                + 'the full histogram needs the autopsy sample from EOS')

    ax.axvline(lo, color=VMM_C, lw=2.4, zorder=4)
    ax.axvline(p05, color=VMM_C, ls='--', lw=1.6, zorder=4)
    ax.axvline(1.0, color=DREAM_C, ls=':', lw=1.4, zorder=4)

    ax.annotate('lowest pulse the VMM' + NL + f'ever recorded: {lo:.2f} x MPV'
                + NL + f'({frac_lost:.0%} of DREAM lies below)',
                xy=(lo, ymax * 0.55), xytext=(1.30, ymax * 0.84),
                fontsize=10, color=VMM_C, ha='left',
                arrowprops=dict(arrowstyle='->', color=VMM_C, lw=1.4),
                bbox=dict(fc='#f4faf4', ec='#8bbf8b',
                          boxstyle='round,pad=0.35'))
    ax.annotate('VMM 5th percentile' + NL + f'{p05:.2f} x MPV',
                xy=(p05, ymax * 0.10), xytext=(2.18, ymax * 0.62),
                fontsize=9.5, color=VMM_C, ha='left',
                arrowprops=dict(arrowstyle='->', color=VMM_C, lw=1.2))
    # Short enough to sit inside the shaded band in the two-panel layout;
    # the number itself rides on the low-edge annotation, which has room.
    ax.text(lo / 2, ymax * 0.62, 'removed by the' + NL + 'discriminator',
            fontsize=9.5, color='#7a2020', ha='center', va='center',
            fontweight='bold')
    ax.text(1.03, ymax * 1.15, 'MPV', color=DREAM_C, fontsize=9.5, va='top')

    if note:
        ax.text(0.985, 0.05, note, transform=ax.transAxes, fontsize=8.4,
                ha='right', color='#666', style='italic',
                bbox=dict(fc='white', ec='#ccc', boxstyle='round,pad=0.35'))

    ax.set_xlim(0, 3.0)
    ax.set_ylim(0, ymax * 1.20)
    ax.set_xlabel("pulse height / that readout's own Landau MPV")
    ax.set_ylabel('normalised density')
    ax.set_title(label + 'P2_OUT, mesh 450 / drift 750, same gas' + NL
                 + 'the VMM window opens while the spectrum is still climbing',
                 fontsize=11.5)
    ax.legend(loc='upper right', framealpha=0.96)


def panel_b(ax):
    xs = [0.5 * (a + b) for (a, b), _, _ in QUINT]
    xe = [[0.5 * (a + b) - a for (a, b), _, _ in QUINT],
          [b - 0.5 * (a + b) for (a, b), _, _ in QUINT]]
    ys = [e for _, e, _ in QUINT]

    ax.axhspan(*DREAM_BAND, color=DREAM_C, alpha=0.12, zorder=0)
    ax.axhline(float(np.mean(DREAM_BAND)), color=DREAM_C, ls='--', lw=1.6,
               zorder=2)
    ax.text(287, 0.9705,
            f'DREAM, this chamber, this working point: '
            f'{DREAM_BAND[0]:.3f} - {DREAM_BAND[1]:.3f}  (3 runs, 25-28 Jul)',
            fontsize=9.4, color=DREAM_C, ha='right', va='bottom')

    for b in PAD_BEST:
        ax.axhline(b, color='#8bbf8b', lw=1.0, ls=':', zorder=1)
    ax.text(88, 0.9495,
            'best individual VMM pads: '
            + ' / '.join(f'{b:.3f}' for b in PAD_BEST),
            fontsize=9, color='#4a7a4a', ha='left', va='top')

    ax.errorbar(xs, ys, xerr=xe, fmt='s', ms=10, lw=2.0, color=VMM_C,
                capsize=4, zorder=4,
                label='VMM pads, grouped in quintiles of their' + NL
                      + 'own median pulse height (run_46)')
    for x, y in zip(xs, ys):
        ax.annotate(f'{y:.3f}', (x, y), textcoords='offset points',
                    xytext=(0, -22), ha='center', fontsize=10, color=VMM_C,
                    fontweight='bold')

    ax.set_xlim(85, 290)
    ax.set_ylim(0.70, 1.005)
    ax.set_xlabel('pad median pulse height [VMM ADC]')
    ax.set_ylabel('per-pad efficiency (uRWELL-referenced)')
    ax.set_title('(b) one global threshold per chip, no per-channel trim' + NL
                 + f'-> pads spread {PAD_P05:.2f}-{PAD_P95:.2f} under one '
                   'uniform beam', fontsize=11.5)
    ax.legend(loc='lower right', framealpha=0.96)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-o', '--out', default='figs')
    ap.add_argument('--vmm-hist', default=None)
    ap.add_argument('--overlay', action='store_true',
                    help='also write adc_overlay_dream_vmm.{png,pdf}: panel '
                         '(a) alone, at slide size')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14.6, 5.8))
    panel_a(ax1, a.vmm_hist)
    panel_b(ax2)
    fig.suptitle('What the discriminator threshold costs, measured on one '
                 'chamber', fontsize=13.5, y=1.02)
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(a.out, f'adc_threshold_cost.{ext}'), dpi=170,
                    bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'  [adc_threshold_cost] -> {a.out}')

    if a.overlay:
        fig2, ax = plt.subplots(figsize=(10.6, 6.6))
        panel_a(ax, a.vmm_hist, label='')
        fig2.suptitle('DREAM and VMM3a pulse height on the same chamber — '
                      'where the discriminator cuts', fontsize=13.5, y=0.98)
        fig2.tight_layout()
        for ext in ('png', 'pdf'):
            fig2.savefig(os.path.join(a.out, f'adc_overlay_dream_vmm.{ext}'),
                         dpi=170, bbox_inches='tight', facecolor='white')
        plt.close(fig2)
        print(f'  [adc_overlay_dream_vmm] -> {a.out}')


if __name__ == '__main__':
    main()
