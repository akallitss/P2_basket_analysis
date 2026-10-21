"""Figure 15 for P2_IN (VMM 2-7) and P2_MID (VMM 8-13): threshold vs masking trade-off per chip."""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from tradeoff import onsets, BASE, eff_d, w

SCR = os.path.dirname(os.path.abspath(__file__)) + '/'
OUT = os.path.expanduser('~/Documents/PostDocSaclay/cern-site/notes/unlisted/'
                         '2026-10-21-testing-the-vmm-threshold-theory-step-by-step_files/')
INK, MUTED, GRID, NOISE = "#17202a", "#5b6774", "#e3e7eb", "#8a6fbf"
PAL = ["#2a78d6", "#1baf7a", "#e0a31b", "#c43a7a", "#6b4fc9", "#eb6834"]
plt.rcParams.update({
    "font.size": 12, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.8, "axes.axisbelow": True, "legend.frameon": False,
    "figure.dpi": 130, "savefig.bbox": "tight", "savefig.facecolor": "white"})
JULY = {int(k): v for k, v in json.load(open(SCR + 'runcfg/sdt.json'))['run_52'].items()}
# chips without a clean edge at a low sdt: baseline extrapolated from sdt 285-300 (less certain); VMM 6 has none
BASE_X = {**BASE, 4: 199, 7: 200, 9: 191}
NOTE = {4: 'walking burst', 7: 'walking burst', 9: 'walking burst', 5: 'whole chip noisy', 12: 'whole chip noisy',
        6: 'baseline unknown (ADC edge ≈ 10)'}
SDT = np.arange(210, 305, 5)
on45, _ = onsets(4.5)

for name, chips, fn in [('P2_IN', range(2, 8), 'fig15b_tradeoff_P2_IN.png'),
                        ('P2_MID', range(8, 14), 'fig15c_tradeoff_P2_MID.png')]:
    fig, axs = plt.subplots(2, 3, figsize=(13.5, 7.4), sharex=True, sharey=True)
    for ax, v, c in zip(axs.flat, chips, PAL):
        nm = np.array([(on45[v] > s).sum() for s in SDT])
        ax.axvspan(195, 207.5, color=GRID, alpha=0.7, lw=0)
        ax.text(201, 0.42, 'not\nscanned', ha='center', fontsize=9, color=MUTED)
        title = f'VMM {v}   ' + (f'(baseline ≈ sdt {BASE_X[v]}{"?" if v in (4, 7, 9) else ""})' if v in BASE_X else '')
        ax.set_title(title, loc='left', fontsize=11.5, color=INK)
        if v in NOTE:
            ax.text(0.98, 0.97, NOTE[v], transform=ax.transAxes, ha='right', va='top', fontsize=10, color=NOISE, fontweight='bold')
        if v not in BASE_X:
            ax.plot(SDT, 1 - nm / 64, 'o-', color=c, lw=2, ms=5, clip_on=False)
            ax.text(0.5, 0.5, 'no efficiency model\n(line = fraction of channels unmasked)', transform=ax.transAxes,
                    ha='center', fontsize=9.5, color=MUTED)
            ax.set_xlim(195, 302); ax.set_ylim(0.3, 1.0)
            continue
        sig = np.array([np.average(eff_d(s - BASE_X[v], 4.5), weights=w) for s in SDT])
        tot = sig * (1 - nm / 64)
        ax.plot(SDT, sig, color=c, lw=1.5, ls='--', label='signal only')
        ax.plot(SDT, tot, 'o-', color=c, lw=2.2, ms=5, label='after masking noisy channels')
        prev = None
        for s_, n_, t_ in zip(SDT, nm, tot):
            show = n_ and n_ != prev and t_ - 0.035 > 0.32
            prev = n_
            if show:
                ax.text(s_, t_ - 0.035, f'{n_}', ha='center', va='top', fontsize=9, color=NOISE)
        i = int(np.argmax(tot))
        ax.plot([SDT[i]], [tot[i]], '*', ms=17, color=INK, zorder=4)
        ax.text(SDT[i] + 1.5, tot[i] + 0.02, f'best {SDT[i]}: {tot[i]*100:.0f} %', fontsize=10, color=INK)
        jv = JULY[v]; jt = tot[list(SDT).index(jv)]
        ax.plot([jv], [jt], 'o', ms=10, mfc='white', mec=INK, mew=1.5, zorder=4)
        ax.text(jv - 1.5, jt - 0.06, f'July {jv}: {jt*100:.0f} %', fontsize=9.5, color=MUTED, ha='right' if jv > 250 else 'left')
        ax.set_xlim(195, 302); ax.set_ylim(0.3, 1.0)
        print(name, v, 'best', SDT[i], round(tot[i], 3), 'july', jv, round(jt, 3), 'mask', nm[i])
    for ax in axs[1]:
        ax.set_xlabel('sdt  [DAC counts]')
    for ax in axs[:, 0]:
        ax.set_ylabel('efficiency (model)')
    axs[0, 0].legend(loc='lower left', fontsize=9.5, bbox_to_anchor=(0.12, 0.12))
    fig.suptitle(f'{name}: same trade-off, 4.5 mV/fC (purple = channels masked; signal shape borrowed from P2_OUT)',
                 x=0.06, ha='left', fontsize=12.5, color=INK)
    fig.tight_layout(); fig.savefig(OUT + fn); plt.close(fig)
