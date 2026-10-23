"""Part 4 figures: what the chip baseline is, how it is estimated from the VMM's own ADC edge, how sure it is,
and how it sets the cosmic threshold scan (figs 19-22)."""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from tradeoff import rows, BASE

SCR = os.path.dirname(os.path.abspath(__file__)) + '/'
OUT = os.path.expanduser('~/Documents/PostDocSaclay/cern-site/notes/unlisted/'
                         '2026-10-21-testing-the-vmm-threshold-theory-step-by-step_files/')
INK, MUTED, GRID, NOISE, THR = "#17202a", "#5b6774", "#e3e7eb", "#8a6fbf", "#eb6834"
STC = {'P2_IN': "#e0a31b", 'P2_MID': "#1baf7a", 'P2_OUT': "#2a78d6"}
plt.rcParams.update({
    "font.size": 12, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.8, "axes.axisbelow": True, "legend.frameon": False,
    "figure.dpi": 130, "savefig.bbox": "tight", "savefig.facecolor": "white"})
ADC_PER_DAC, PED_ADC, DAC_MV, DAC_OFF = 0.9, 38.0, 0.801, 25.18
EXTRAP = {4, 7, 9}
st = lambda v: 'P2_IN' if v < 8 else 'P2_MID' if v < 14 else 'P2_OUT'
z = np.load(SCR + 'vmm_adc_all.npz')
SDT = json.load(open(SCR + 'runcfg/sdt.json'))
RUNS = {'run_46': 'run_46__cfg_gain4.5_peaktime200', 'run_48': 'run_48__cfg_gain4.5_peaktime200_opt',
        'run_52': 'run_52__cfg_gain4.5_peaktime200_opt', 'run_66': 'run_66__cfg_gain4.5_peaktime200_opt'}


def edge(h, q=0.01):
    c = np.cumsum(h) / h.sum()
    return int(np.searchsorted(c, q))


# per chip, per run: (sdt, edge, baseline estimate)
pts = {v: [] for v in range(2, 20)}
for r, k in RUNS.items():
    a = z[k + '__adc'].astype(float)
    for v in range(2, 20):
        s, e = SDT[r][str(v)], edge(a[v])
        pts[v].append((r, s, e, s - (e - PED_ADC) / ADC_PER_DAC))
num = {}
base, err_run, err_sys = {}, {}, {}
for v in range(2, 20):
    if v == 6:
        continue
    b = np.array([p[3] for p in pts[v] if p[0] != 'run_46'])       # the three July "opt" runs, as in Part 3
    base[v] = float(b.mean()); err_run[v] = float(b.std())
    lever = np.mean([p[2] for p in pts[v] if p[0] != 'run_46']) - PED_ADC
    # slope 0.85-0.95 ADC/DAC and pedestal 38 +- 5 ADC, added in quadrature
    err_sys[v] = float(np.hypot(lever / ADC_PER_DAC * 0.05 / 0.9, 5 / ADC_PER_DAC))
num['baseline'] = {v: dict(base=round(base[v], 1), run_spread=round(err_run[v], 1), sys=round(err_sys[v], 1)) for v in base}

# ------------------------------------------------------------- fig 19: cartoon ---
def vmm_pulse(t, tp):
    a = tp / 1.5
    p0, p1 = 1.263 / a, (1.149 - 0.789j) / a
    K1, K2 = 1.584, -0.792 - 0.115j
    f = K1 * np.exp(-t * p0) + 2 * np.abs(K2) * np.exp(-t * p1.real) * np.cos(-t * p1.imag + np.angle(K2))
    f = np.where(t > 0, f, 0.0)
    return f / f.max()

fig, (ax, bx) = plt.subplots(1, 2, figsize=(14, 5.3), gridspec_kw=dict(width_ratios=[1.35, 1]))
t = np.linspace(-80, 760, 1500)
sdt = 230; thr = DAC_MV * sdt + DAC_OFF                       # 209.4 mV
peds = [(DAC_MV * 182 + DAC_OFF, 'VMM 18 (baseline ≈ sdt 182)', "#6b4fc9"), (DAC_MV * 201 + DAC_OFF, 'VMM 14 (baseline ≈ sdt 201)', "#2a78d6")]
rng = np.random.default_rng(2)
for (pb, lab, c), t0 in zip(peds, (60, 60)):
    noise = np.convolve(rng.normal(0, 1, t.size), np.ones(20) / 20 * 4.5, 'same')
    y = pb + 42 * vmm_pulse(t - t0, 200) + noise
    ax.plot(t, y, color=c, lw=1.8)
    ax.axhline(pb, color=c, lw=1, ls=':')
    ax.text(750, pb - 1.2, lab, color=c, fontsize=9.5, ha='right', va='top')
ax.axhline(thr, color=THR, lw=2)
ax.text(640, thr - 1.2, f'threshold, sdt {sdt}:\n0.801·{sdt} + 25.2 = {thr:.0f} mV\n(the same on both chips)', color=THR,
        fontsize=10, ha='right', va='top', bbox=dict(fc='white', ec='none', alpha=0.85, pad=1.5))
for (pb, _, c), xx in zip(peds, (655, 725)):
    ax.annotate('', (xx, thr), (xx, pb), arrowprops=dict(arrowstyle='<->', color=c, lw=1.6))
ax.text(745, 230.5, 'd = sdt − baseline:\nVMM 18: 230 − 182 = 48 DAC ≈ 38 mV ≈ 43 ADC\nVMM 14: 230 − 201 = 29 DAC ≈ 23 mV ≈ 26 ADC',
        fontsize=9.8, color=INK, ha='right', va='top', linespacing=1.4)
ax.set_xlim(-80, 760); ax.set_ylim(160, 232)
ax.set_xlabel('time [ns]'); ax.set_ylabel('shaper output [mV, absolute]')
ax.set_title('a   one sdt, two chips: what counts is the distance d from each chip’s baseline', loc='left', fontsize=11.5, color=INK)
# b: what the ADC sees
x = np.arange(0, 400)
mp, w = 120, 22
land = np.exp(-((x - mp) / w) - np.exp(-(x - mp) / w))
for pa, e, c, lab in ((PED_ADC, 81, "#6b4fc9", 'VMM 18'), (PED_ADC, 64, "#2a78d6", 'VMM 14')):
    h = np.where(x >= e, land, 0)
    bx.fill_between(x, h / h.max() * (0.95 if lab == 'VMM 18' else 0.75), step='mid', color=c, alpha=0.25, lw=0)
    bx.step(x, h / h.max() * (0.95 if lab == 'VMM 18' else 0.75), where='mid', color=c, lw=1.6, label=f'{lab}: edge at {e} ADC')
bx.axvline(PED_ADC, color=MUTED, lw=1.5, ls='--')
bx.text(PED_ADC + 3, 0.62, 'baseline\n≈ 38 ADC\n(never\nrecorded)', fontsize=9.5, color=MUTED, va='top')
bx.plot(x, land / land.max() * 0.95, color=MUTED, lw=1, ls=':', label='signal below threshold (lost)')
bx.set_xlim(0, 330); bx.set_ylim(0, 1.1); bx.set_yticks([])
bx.set_xlabel('VMM ADC (peak height, baseline included)'); bx.set_ylabel('hits')
bx.legend(loc='upper right', fontsize=9.5)
bx.set_title('b   the chip only records what is above threshold', loc='left', fontsize=11.5, color=INK)
fig.tight_layout(); fig.savefig(OUT + 'fig19_baseline_cartoon.png'); plt.close(fig)

# ----------------------------------------------- fig 20: edge -> baseline per chip ---
fig, axs = plt.subplots(1, 3, figsize=(15, 5.2), sharey=True)
for ax, (name, chips) in zip(axs, (('P2_IN', range(2, 8)), ('P2_MID', range(8, 14)), ('P2_OUT', range(14, 20)))):
    cols = plt.cm.tab10(np.linspace(0, 0.55, 6))
    for c, v in zip(cols, chips):
        if v == 6:
            s, e = pts[v][2][1], pts[v][2][2]
            ax.plot(s, e, 'x', color=c, ms=10, mew=2)
            ax.annotate('VMM 6: edge ≈ 10 ADC,\nbelow the baseline →\nno estimate', (s, e), xytext=(8, 30), textcoords='offset points',
                        fontsize=9, color=INK)
            continue
        ss = np.array([p[1] for p in pts[v]]); ee = np.array([p[2] for p in pts[v]])
        ax.plot(ss, ee, 'o', color=c, ms=7, zorder=3, mfc='white' if v in EXTRAP else c, mew=1.8)
        xs = np.linspace(base[v], max(ss), 50)
        ax.plot(xs, PED_ADC + ADC_PER_DAC * (xs - base[v]), color=c, lw=1.3, ls='--' if v in EXTRAP else '-', alpha=0.85)
        ax.plot(base[v], PED_ADC, 'v', color=c, ms=9, zorder=4)
        rank = sorted(b_ for b_ in [base[u] for u in chips if u in base]).index(base[v])
        ax.text(base[v], PED_ADC - 5 - 8 * (rank % 3), f'{v}', color=c, fontsize=9.5, ha='center', va='top', fontweight='bold')
    ax.axhline(PED_ADC, color=MUTED, lw=1)
    ax.text(302, PED_ADC + 2, 'baseline 38 ADC', fontsize=9, color=MUTED, ha='right')
    ax.set_xlim(172, 305); ax.set_ylim(0, 150)
    ax.set_xlabel('sdt  [DAC counts]')
    ax.set_title(name, loc='left', fontsize=12, color=INK, fontweight='bold')
axs[0].set_ylabel('edge = lowest ADC the chip records  [VMM ADC]')
axs[0].plot([], [], 'o', color=INK, label='measured edge (runs 46, 48, 52, 66)')
axs[0].plot([], [], 'o', mfc='white', mec=INK, mew=1.8, label='long extrapolation (only sdt 285–300)')
axs[0].plot([], [], 'v', color=INK, label='baseline estimate = where the line hits 38 ADC')
axs[0].legend(loc='upper left', fontsize=9.2)
fig.suptitle('slide each chip’s measured edge down its 0.9 ADC/DAC line until it reaches the baseline: the sdt there is the chip’s baseline',
             x=0.05, ha='left', fontsize=12, color=INK)
fig.tight_layout(); fig.savefig(OUT + 'fig20_edge_to_baseline.png'); plt.close(fig)

# ------------------------------------------------- fig 21: how sure, and a cross-check ---
fig, (ax, bx) = plt.subplots(1, 2, figsize=(14.5, 5.2), gridspec_kw=dict(width_ratios=[1.25, 1]))
vs = [v for v in range(2, 20)]
for k, v in enumerate(vs):
    c = STC[st(v)]
    if v not in base:
        ax.text(k, 205, 'no\nestimate', ha='center', fontsize=9, color=MUTED); continue
    ax.errorbar(k, base[v], yerr=err_sys[v], fmt='none', ecolor=c, elinewidth=7, alpha=0.25, capsize=0)
    ax.errorbar(k, base[v], yerr=max(err_run[v], 0.3), fmt='o', color=c, ms=8, mfc='white' if v in EXTRAP else c, mew=2, capsize=3)
    ax.text(k, base[v] + err_sys[v] + 1.5, f'{base[v]:.0f}', ha='center', fontsize=9.5, color=INK)
for name, c in STC.items():
    i0 = [k for k, v in enumerate(vs) if st(v) == name]
    ax.axvspan(i0[0] - 0.5, i0[-1] + 0.5, color=c, alpha=0.07, lw=0)
    ax.text((i0[0] + i0[-1]) / 2, 226, name, ha='center', fontsize=11.5, fontweight='bold', color=INK)
ax.set_xticks(range(len(vs)), [str(v) for v in vs]); ax.set_xlabel('VMM id'); ax.set_ylim(168, 230)
ax.set_xlim(-0.6, len(vs) - 0.4); ax.grid(axis='x', visible=False)
ax.set_ylabel('chip baseline  [sdt DAC counts]')
ax.legend(handles=[plt.Line2D([], [], marker='o', ls='', color=INK, label='spread over runs 48 / 52 / 66'),
                   Patch(color=MUTED, alpha=0.3, label='systematic: baseline 38 ± 5 ADC, slope 0.85–0.95')],
          loc='lower left', fontsize=9.3)
ax.set_title('a   chip baselines: run-to-run ±1 DAC, systematic ±6–8 DAC', loc='left', fontsize=11.5, color=INK)
# b: bench noise at sdt 210 vs distance from baseline
rr = sorted([r for r in rows if r['gain'] == 4.5], key=lambda r: r['sdt'])
r210 = [r for r in rr if r['sdt'] == 210][0]
for v in base:
    n = int((np.array(r210['per'][str(v)]['ch']) / r210['span_s'] > 10).sum())
    d = 210 - base[v]
    bx.errorbar(d, n, xerr=err_sys[v], fmt='o', color=STC[st(v)], ms=8, mfc='white' if v in EXTRAP else STC[st(v)],
                mew=2, elinewidth=1, alpha=0.9, capsize=0)
    bx.annotate(f'{v}', (d, n), xytext=(5, 4), textcoords='offset points', fontsize=9, color=INK)
bx.set_xlabel('d at sdt 210 = 210 − baseline  [DAC]'); bx.set_ylabel('channels above 10 Hz on the bench (of 64)')
bx.set_xlim(-12, 38); bx.set_ylim(-3, 66)
bx.axvspan(-12, 0, color=GRID, alpha=0.6, lw=0); bx.text(-6, 60, 'threshold\nbelow\nbaseline', ha='center', fontsize=9, color=MUTED)
bx.set_title('b   cross-check: chips close to their baseline are the noisy ones', loc='left', fontsize=11.5, color=INK)
fig.tight_layout(); fig.savefig(OUT + 'fig21_baseline_uncertainty.png'); plt.close(fig)

# -------------------------------------------------------- fig 22: the cosmic scan grid ---
D_NOISE = [50, 40, 34, 30, 27, 24, 22, 20, 18, 16]
D_COSMIC = [50, 40, 32, 27, 24, 22, 20, 18]
JULY = {int(k): v for k, v in SDT['run_52'].items()}
fig, ax = plt.subplots(figsize=(14, 5.6))
for k, v in enumerate(vs):
    c = STC[st(v)]
    if v not in base:
        for s in range(230, 199, -5):
            ax.plot(k, s, 's', color=MUTED, ms=5, mfc='white')
        ax.text(k, 186, 'start 230,\nsteps of 5', ha='center', fontsize=8.5, color=MUTED)
    else:
        b = base[v]
        ax.add_patch(plt.Rectangle((k - 0.38, b - err_sys[v]), 0.76, 2 * err_sys[v], color=GRID, alpha=0.9, lw=0, zorder=0))
        ax.plot([k - 0.38, k + 0.38], [b, b], color=MUTED, lw=1.5)
        ax.plot([k - 0.18] * len(D_NOISE), [b + d for d in D_NOISE], 's', color=NOISE, ms=4.5, alpha=0.9)
        ax.plot([k + 0.18] * len(D_COSMIC), [b + d for d in D_COSMIC], 'o', color=c, ms=5.5, mfc='white' if v in EXTRAP else c, mew=1.5)
    ax.plot(k, JULY[v], '_', color=INK, ms=22, mew=2.5, zorder=5)
for name, c in STC.items():
    i0 = [k for k, v in enumerate(vs) if st(v) == name]
    ax.axvspan(i0[0] - 0.5, i0[-1] + 0.5, color=c, alpha=0.06, lw=0)
    ax.text((i0[0] + i0[-1]) / 2, 306, name, ha='center', fontsize=11.5, fontweight='bold', color=INK)
ax.set_xticks(range(len(vs)), [str(v) for v in vs]); ax.set_xlabel('VMM id'); ax.set_xlim(-0.6, len(vs) - 0.4)
ax.set_ylim(165, 315); ax.set_ylabel('sdt  [DAC counts]'); ax.grid(axis='x', visible=False)
ax.legend(handles=[plt.Line2D([], [], marker='s', ls='', color=NOISE, ms=6, label='noise scan: baseline + 50 … 16'),
                   plt.Line2D([], [], marker='o', ls='', color=INK, ms=6, label='cosmic scan: baseline + 50 … 18'),
                   plt.Line2D([], [], marker='_', ls='', color=INK, ms=16, mew=2.5, label='July setting (run_52)'),
                   Patch(color=GRID, label='baseline ± systematic')],
          loc='upper center', bbox_to_anchor=(0.5, -0.13), ncol=4, fontsize=9.8)
ax.set_title('the cosmic scan steps each chip down from its own baseline (to be re-centred on the vmmsc pedestals)', loc='left',
             fontsize=11.5, color=INK)
fig.savefig(OUT + 'fig22_cosmic_scan_grid.png'); plt.close(fig)

json.dump(num, open(SCR + 'part4_numbers.json', 'w'), indent=1)
for v in base:
    print(v, f"{base[v]:.1f} ±{err_run[v]:.1f} (run) ±{err_sys[v]:.1f} (sys)  edges", [(p[1], p[2]) for p in pts[v]])
