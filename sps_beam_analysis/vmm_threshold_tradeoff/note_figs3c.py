"""Summary of the threshold/masking trade-off for all 18 chips: best efficiency and masks at the best sdt."""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from tradeoff import onsets, eff_d, w
from note_figs3b import BASE_X, JULY, SDT, OUT, INK, MUTED, GRID, NOISE

ST = {'P2_IN': (range(2, 8), "#e0a31b"), 'P2_MID': (range(8, 14), "#1baf7a"), 'P2_OUT': (range(14, 20), "#2a78d6")}
BURST = {4, 7, 9}
GAIN = float(os.environ.get('GAIN', 4.5))          # bench scan exists at 3.0 and 4.5 mV/fC
SFX = '' if GAIN == 4.5 else f"_g{GAIN:.0f}"
on, on45 = onsets(GAIN)[0], onsets(4.5)[0]        # July reference is always the real July point (4.5)
rows = []
for st, (chips, c) in ST.items():
    for v in chips:
        nm = np.array([(on[v] > s).sum() for s in SDT])
        if v in BASE_X:
            tot = np.array([np.average(eff_d(s - BASE_X[v], GAIN), weights=w) for s in SDT]) * (1 - nm / 64)
            i = int(np.argmax(tot)); best, eb = int(SDT[i]), float(tot[i])
            ej = np.average(eff_d(JULY[v] - BASE_X[v], 4.5), weights=w) * (1 - (on45[v] > JULY[v]).sum() / 64)
        else:   # VMM 6: no baseline -> lowest sdt with no noisy channel, no efficiency
            i = int(np.argmax(nm == 0)); best, eb, ej = int(SDT[i]), np.nan, np.nan
        rows.append(dict(vmm=v, st=st, c=c, best=best, eff=eb, july=ej, nmask=int(nm[i]), july_sdt=JULY[v],
                         masked=[int(k) for k in np.where(on[v] > best)[0]], floor=best == 210))
json.dump(rows, open(os.path.dirname(os.path.abspath(__file__)) + f'/part3_allchips{SFX}.json', 'w'), indent=1)
x = np.arange(len(rows)); labels = [f"{r['vmm']}" for r in rows]


def station_axis(ax, ytxt):
    for st, (chips, c) in ST.items():
        i0 = [k for k, r in enumerate(rows) if r['st'] == st]
        ax.axvspan(i0[0] - 0.5, i0[-1] + 0.5, color=c, alpha=0.07, lw=0)
        ax.text((i0[0] + i0[-1]) / 2, ytxt, st, ha='center', va='bottom', fontsize=12, fontweight='bold', color=INK)
    ax.set_xticks(x, labels); ax.set_xlabel('VMM id'); ax.set_xlim(-0.6, len(rows) - 0.4)
    ax.grid(axis='x', visible=False)


# ---------------------------------------------- fig 17: best efficiency ---
fig, ax = plt.subplots(figsize=(12.5, 5.4))
for k, r in enumerate(rows):
    if np.isnan(r['eff']):
        ax.text(k, 0.47, 'no\nmodel', ha='center', fontsize=9, color=MUTED); continue
    ax.plot([k, k], [r['july'], r['eff']], color=GRID, lw=3, zorder=1, solid_capstyle='round')
    ax.plot(k, r['july'], 'o', ms=9, mfc='white', mec=MUTED, mew=1.5, zorder=3)
    hollow = r['vmm'] in BURST
    ax.plot(k, r['eff'], 'o', ms=12, color=r['c'], mfc='white' if hollow else r['c'], mew=2.2, zorder=4)
    below = r['eff'] < r['july']        # label under the point so it does not sit on the July marker
    ax.text(k, r['eff'] + (-0.05 if below else 0.018), f"{r['eff']*100:.0f}", ha='center', fontsize=10.5, color=INK, fontweight='bold')
    if r['floor']:
        ax.annotate('', (k + 0.3, r['eff'] + 0.035), (k + 0.3, r['eff'] - 0.005), arrowprops=dict(arrowstyle='->', color=MUTED, lw=1.3))
    ax.text(k, 0.365, f"sdt\n{r['best']}", ha='center', va='bottom', fontsize=8.5, color=MUTED)
ax.axhline(0.956, color="#eb6834", lw=1.3, ls='--'); ax.text(4, 0.962, 'DREAM on P2_OUT\n95.6 %', ha='center',
                                                            va='bottom', fontsize=10, color=INK)
ax.plot([], [], 'o', ms=12, color=INK, label='best reachable (global sdt + masks)')
ax.plot([], [], 'o', ms=9, mfc='white', mec=MUTED, mew=1.5, label='July setting (4.5 mV/fC, as run)')
ax.plot([], [], 'o', ms=12, mfc='white', mec=INK, mew=2.2, label='walking burst: upper bound')
ax.plot([], [], ls='', marker=r'$\uparrow$', ms=11, color=MUTED, label='at the scan floor (210): could go lower')
ax.set_ylim(0.36, 1.06); ax.set_ylabel('efficiency (model)')
station_axis(ax, 1.01)
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.14), ncol=4, fontsize=9.5)
ax.set_title(f'best efficiency each chip can reach at {GAIN:.1f} mV/fC (P2_IN / P2_MID use P2_OUT signal spectra)', loc='left',
             fontsize=11.5, color=INK, pad=22)
fig.savefig(OUT + f'fig17_best_eff_per_vmm{SFX}.png'); plt.close(fig)

# ------------------------------------------- fig 18: masks at the optimum ---
fig, ax = plt.subplots(figsize=(12.5, 4.6))
for k, r in enumerate(rows):
    ax.bar(k, r['nmask'], width=0.62, color=r['c'], hatch='//' if r['vmm'] in BURST else None,
           edgecolor='white' if r['vmm'] not in BURST else r['c'], fill=r['vmm'] not in BURST, lw=1.5)
    lab = ', '.join(map(str, r['masked'])) if 0 < r['nmask'] <= 3 else ''
    ax.text(k, r['nmask'] + 0.25, f"{r['nmask']}", ha='center', fontsize=11, fontweight='bold', color=INK)
    if lab:
        ax.text(k, r['nmask'] + 1.6, f'ch {lab}', ha='center', fontsize=8.5, color=MUTED, rotation=90, va='bottom')
    ax.text(k, -1.7, f"{r['best']}", ha='center', fontsize=8.5, color=MUTED)
ax.text(-0.55, -1.7, 'sdt', ha='right', fontsize=8.5, color=MUTED)
ax.set_ylim(-2.4, 17.5); ax.set_ylabel('channels masked (of 64)')
ax.axhline(0, color=MUTED, lw=0.8)
station_axis(ax, 16.2)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(facecolor='white', hatch='//', edgecolor=MUTED, label='walking burst: a fixed mask list will not hold')],
          loc='upper left', bbox_to_anchor=(0.0, 0.86), fontsize=9.5)
ax.set_yticks(range(0, 16, 2))
ax.set_title(f"{GAIN:.1f} mV/fC: channels to mask at each chip's best sdt — {sum(r['nmask'] for r in rows)} of {64*len(rows)} in total,"
             f" {sum(r['nmask'] for r in rows if r['vmm'] not in BURST)} without the walking-burst chips", loc='left', fontsize=11.5, color=INK)
fig.savefig(OUT + f'fig18_masks_per_vmm{SFX}.png'); plt.close(fig)
for r in rows:
    print(r['vmm'], r['st'], r['july_sdt'], '->', r['best'], f"{r['july']:.3f} -> {r['eff']:.3f}", r['nmask'], r['masked'])
