"""Part 4, where the 38 ADC baseline comes from (figs 23-25): two rulers for one threshold, the fit and how robust
it is, and how a neighbour-logic run would measure the baseline directly."""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCR = os.path.dirname(os.path.abspath(__file__)) + '/'
OUT = os.path.expanduser('~/Documents/PostDocSaclay/cern-site/notes/unlisted/'
                         '2026-10-21-testing-the-vmm-threshold-theory-step-by-step_files/')
INK, MUTED, GRID, THR, DREAM, VMM = "#17202a", "#5b6774", "#e3e7eb", "#eb6834", "#eb6834", "#2a78d6"
CHIPC = {16: "#e0a31b", 17: "#c43a7a", 18: "#6b4fc9", 19: "#eb6834"}
plt.rcParams.update({
    "font.size": 12, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.8, "axes.axisbelow": True, "legend.frameon": False,
    "figure.dpi": 130, "savefig.bbox": "tight", "savefig.facecolor": "white"})
z = np.load(SCR + 'vmm_adc_all.npz')
nn = json.load(open(SCR + 'note_numbers.json'))


def edge(h, q=0.01):
    c = np.cumsum(h) / h.sum()
    return int(np.searchsorted(c, q))


RUN = {'46': 'run_46__cfg_gain4.5_peaktime200', '48': 'run_48__cfg_gain4.5_peaktime200_opt'}
P = []                                   # (chip, run, T, edge)
for c in nn['chips']:
    for r in ('46', '48'):
        P.append((c['chip'], r, c['T' + r], edge(z[RUN[r] + '__adc'][c['chip']])))
T = np.array([p[2] for p in P], float); E = np.array([p[3] for p in P], float)
fit = lambda m: np.polyfit(T[m], E[m], 1)
slope, icpt = fit(np.ones(len(P), bool))
num = dict(points=[dict(chip=p[0], run=p[1], T=p[2], edge=p[3]) for p in P], slope=float(slope), intercept=float(icpt))

# ----------------------------------------------------- fig 23: two rulers ---
fig, (ax, bx) = plt.subplots(1, 2, figsize=(14.5, 5.0), gridspec_kw=dict(width_ratios=[1.2, 1]))
ax.set_xlim(-30, 330); ax.set_ylim(-0.6, 3.2); ax.axis('off')
for y, lab, col, dy in ((2.4, 'threshold as charge, from the efficiency fit (Part 1)   [DREAM ADC]', DREAM, 0.42),
                        (0.6, 'threshold as seen by the VMM: its ADC edge   [VMM ADC]', VMM, -0.62)):
    ax.plot([0, 300], [y, y], color=col, lw=2.5)
    ax.text(0, y + dy, lab, color=col, fontsize=11, fontweight='bold')
for t in range(0, 301, 50):
    ax.plot([t, t], [2.4, 2.48], color=DREAM, lw=1.5); ax.text(t, 2.55, f'{t}', ha='center', fontsize=9.5, color=DREAM)
for e_ in range(40, 101, 10):
    t = (e_ - icpt) / slope
    if 0 <= t <= 300:
        ax.plot([t, t], [0.52, 0.6], color=VMM, lw=1.5); ax.text(t, 0.3, f'{e_}', ha='center', fontsize=9.5, color=VMM)
ax.plot([0, 0], [0.52, 0.6], color=VMM, lw=1.5); ax.text(0, 0.3, f'{icpt:.0f}', ha='center', fontsize=10.5, color=VMM, fontweight='bold')
for (v, r, t, e_) in P:
    ax.plot([t, t], [2.4, 0.6], color=CHIPC[v], lw=1.2, alpha=0.8)
    ax.plot([t], [2.4], 'o', color=CHIPC[v], ms=6); ax.plot([t], [0.6], 'o', color=CHIPC[v], ms=6)
ax.annotate('', (0, 0.62), (0, 2.38), arrowprops=dict(arrowstyle='->', color=INK, lw=2, ls='--'))
ax.text(8, 1.5, 'zero charge on the top ruler\nlands at ~38 on the bottom one:\nthat is the baseline', fontsize=10.5, color=INK, va='center')
ax.text(150, -0.45, 'each coloured line = one chip in one run (VMM 16–19, runs 46 and 48)', ha='center', fontsize=9.5, color=MUTED)
ax.set_title('a   one threshold, measured on two rulers', loc='left', fontsize=11.5, color=INK)
# b: why the intercept is the baseline
x = np.linspace(0, 1, 10)
base_v = 38
for i, (tt, lab) in enumerate(((250, 'high threshold'), (125, 'lower'), (0, 'threshold on the baseline'))):
    xc = i * 1.25; h = 0.196 * tt
    bx.add_patch(plt.Rectangle((xc, base_v), 0.8, h + 0.01, color=THR, alpha=0.35, lw=0))
    bx.plot([xc - 0.05, xc + 0.85], [base_v, base_v], color=MUTED, lw=1.5)
    bx.plot([xc - 0.05, xc + 0.85], [base_v + h, base_v + h], color=THR, lw=2.2)
    bx.text(xc + 0.4, base_v + h + 1.2, f'edge = 38 + {h:.0f} = {38 + h:.0f}', ha='center', fontsize=10, color=VMM)
    bx.text(xc + 0.4, base_v - 3.2, f'T = {tt} DREAM ADC', ha='center', fontsize=9.5, color=DREAM)
    bx.text(xc + 0.4, base_v - 6.2, lab, ha='center', fontsize=9.5, color=INK)
bx.text(1.65, 92, 'edge = baseline + 0.196 × T', ha='center', fontsize=11, color=INK, fontweight='bold')
bx.text(-0.05, base_v + 0.6, 'baseline', fontsize=9.5, color=MUTED, va='bottom')
bx.set_xlim(-0.2, 3.6); bx.set_ylim(28, 96); bx.axis('off')
bx.set_title('b   slide the threshold down: the edge slides to the baseline', loc='left', fontsize=11.5, color=INK)
fig.tight_layout(); fig.savefig(OUT + 'fig23_two_rulers.png'); plt.close(fig)

# --------------------------------------------- fig 24: the fit and its robustness ---
fig, (ax, bx) = plt.subplots(1, 2, figsize=(14.5, 5.3), gridspec_kw=dict(width_ratios=[1.2, 1]))
rng = np.random.default_rng(1); boot = []
xx = np.linspace(0, 300, 50)
for _ in range(3000):
    k = rng.integers(0, len(P), len(P))
    if len(set(T[k])) < 3:
        continue
    b = np.polyfit(T[k], E[k], 1); boot.append(b)
    if len(boot) <= 150:
        ax.plot(xx, np.polyval(b, xx), color=MUTED, lw=0.6, alpha=0.12)
boot = np.array(boot); lo, mid, hi = np.percentile(boot[:, 1], [16, 50, 84])
ax.plot(xx, np.polyval([slope, icpt], xx), color=INK, lw=2, label=f'fit: edge = {icpt:.1f} + {slope:.3f} · T')
for (v, r, t, e_) in P:
    sus = v == 19 and r == '46'
    ax.plot(t, e_, 'o', color=CHIPC[v], ms=10, mfc='white' if sus else CHIPC[v], mew=2.2, zorder=4)
    off = {(17, '48'): (-92, 6), (18, '48'): (8, 4), (17, '46'): (8, -12), (16, '48'): (8, -10)}.get((v, r), (8, -12 if r == '46' else 6))
    ax.annotate(f'VMM {v}, run_{r}', (t, e_), xytext=off, textcoords='offset points', fontsize=9, color=INK)
ax.annotate('noise-dominated run:\nedge pulled down?', (92, 58), xytext=(30, 75), fontsize=9.5, color=INK,
            arrowprops=dict(arrowstyle='->', color=MUTED))
ax.axvspan(-5, 88, color=GRID, alpha=0.5, lw=0); ax.text(42, 97, 'no data here:\nthe line is extended', ha='center', fontsize=9.5, color=MUTED)
ax.plot(0, icpt, 's', color=INK, ms=9, zorder=5)
ax.text(5, icpt - 5, f'T = 0  →  {icpt:.0f} ADC', fontsize=10.5, color=INK, fontweight='bold')
ax.set_xlim(-5, 300); ax.set_ylim(25, 100)
ax.set_xlabel('threshold fitted from efficiency, T  [DREAM ADC]'); ax.set_ylabel('edge in the VMM ADC  [VMM ADC]')
ax.plot([], [], color=MUTED, lw=1, alpha=0.5, label='same fit on resampled points')
ax.legend(loc='lower right', fontsize=9.8)
ax.set_title('a   the 8 points and the straight line through them', loc='left', fontsize=11.5, color=INK)
# b: how much the answer moves
var = []
allm = np.ones(len(P), bool)
var.append(('all 8 points (the value used)', icpt, None))
var.append(('resampled points (68 % range)', mid, (lo, hi)))
m = allm.copy(); m[[i for i, p in enumerate(P) if p[0] == 19 and p[1] == '46']] = False
var.append(('without VMM 19 run_46', fit(m)[1], None))
var.append(('run_46 points only', fit(np.array([p[1] == '46' for p in P]))[1], None))
var.append(('run_48 points only', fit(np.array([p[1] == '48' for p in P]))[1], None))
for v in (16, 17, 18, 19):
    m = np.array([p[0] == v for p in P]); var.append((f'VMM {v} alone (its 2 points)', fit(m)[1], None))
ys = np.arange(len(var))[::-1]
bx.axvspan(icpt - 5, icpt + 5, color=VMM, alpha=0.10, lw=0)
bx.axvspan(40, 50, color=DREAM, alpha=0.12, lw=0)
for y, (lab, val, rngv) in zip(ys, var):
    if rngv:
        bx.plot(rngv, [y, y], color=INK, lw=2.5)
    bx.plot(val, y, 'o', color=INK if y == ys[0] else MUTED, ms=9 if y == ys[0] else 7)
    bx.text(val + 0.6, y + 0.18, f'{val:.0f}', fontsize=9.5, color=INK)
bx.set_yticks(ys, [v[0] for v in var], fontsize=10); bx.set_xlim(25, 55)
bx.text(icpt, len(var) - 0.25, 'used: 38 ± 5', ha='center', fontsize=9.5, color=VMM)
bx.text(45, len(var) - 0.25, 'November:\nneighbour hits', ha='center', fontsize=9.5, color=DREAM, va='bottom')
bx.set_ylim(-0.7, len(var) + 0.6)
bx.set_xlabel('baseline from the intercept  [VMM ADC]'); bx.grid(axis='y', visible=False)
bx.set_title('b   how much the answer moves', loc='left', fontsize=11.5, color=INK)
fig.tight_layout(); fig.savefig(OUT + 'fig24_baseline_fit.png'); plt.close(fig)
num['variants'] = {lab: float(val) for lab, val, _ in var}; num['boot_68'] = [float(lo), float(hi)]

# -------------------------------------------- fig 25: neighbour logic, direct ---
fig, (ax, bx) = plt.subplots(1, 2, figsize=(14, 4.8), gridspec_kw=dict(width_ratios=[1, 1]))
chs = np.arange(5); peak = np.array([2, 9, 46, 11, 1]) + 0.0          # mV above baseline, one muon on the middle pad
thr = 25
cols = ['#c9d3dd', VMM, VMM, VMM, '#c9d3dd']
ax.bar(chs, peak, width=0.62, bottom=0, color=[VMM if p > thr else '#9db8d8' for p in peak], edgecolor='white')
ax.axhline(thr, color=THR, lw=2); ax.text(4.45, thr + 1, 'threshold', color=THR, ha='right', fontsize=10)
ax.axhline(0, color=MUTED, lw=1.5); ax.text(4.45, -4, 'baseline', color=MUTED, ha='right', fontsize=10)
ax.text(2, 48, 'fires', ha='center', fontsize=10.5, color=INK, fontweight='bold')
for c_ in (1, 3):
    ax.text(c_, peak[c_] + 1.5, 'recorded\nanyway', ha='center', fontsize=9.5, color=INK)
for c_ in (0, 4):
    ax.text(c_, 4, 'not\nrecorded', ha='center', fontsize=9, color=MUTED)
ax.set_xticks(chs, ['ch n−2', 'ch n−1', 'ch n', 'ch n+1', 'ch n+2']); ax.set_yticks([])
ax.set_ylim(-8, 58); ax.set_ylabel('pulse height above baseline')
ax.grid(False)
ax.set_title('a   neighbour logic on (sng = 1): the neighbours of a fired channel are read out too', loc='left', fontsize=11, color=INK)
xa = np.arange(0, 200)
rng = np.random.default_rng(4)
neigh = np.concatenate([rng.normal(44, 4, 60000), 44 + rng.exponential(18, 9000)])
sig = 44 + rng.gamma(4.5, 18, 40000); sig = sig[sig > 44 + 30]
bx.hist(neigh, bins=xa, color='#9db8d8', alpha=0.9, label='neighbour hits (below threshold)')
bx.hist(sig, bins=xa, color=VMM, alpha=0.75, label='fired channels (above threshold)')
bx.axvline(44, color=MUTED, lw=1.5, ls='--'); bx.text(47, bx.get_ylim()[1] * 0.92, 'peak = that channel’s\nbaseline, read directly', fontsize=10, color=INK, va='top')
bx.axvline(44 + 30, color=THR, lw=1.5); bx.text(44 + 32, bx.get_ylim()[1] * 0.45, 'threshold', fontsize=10, color=THR)
bx.set_xlim(0, 200); bx.set_yticks([]); bx.set_xlabel('VMM ADC'); bx.set_ylabel('hits (one channel)')
bx.legend(loc='upper right', fontsize=9.5)
bx.set_title('b   what one channel would show (illustration)', loc='left', fontsize=11, color=INK)
fig.tight_layout(); fig.savefig(OUT + 'fig25_neighbour_baseline.png'); plt.close(fig)

json.dump(num, open(SCR + 'part4b_numbers.json', 'w'), indent=1)
print(json.dumps(num, indent=1))
