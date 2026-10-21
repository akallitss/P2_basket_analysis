"""Figure 14b: bench noise vs run_46 out-of-spill hits, per channel; and what the out-of-spill hits look like."""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
SCR = os.path.dirname(os.path.abspath(__file__)) + '/'
sys.path.insert(0, SCR + 'oos')
from oos_load import load, SDT46
from tradeoff import rows as BENCH, BASE
INK, MUTED, GRID = "#17202a", "#5b6774", "#e3e7eb"
STC = {'P2_IN': "#e0a31b", 'P2_MID': "#1baf7a", 'P2_OUT': "#2a78d6"}
OUT = os.path.expanduser('~/Documents/PostDocSaclay/cern-site/notes/unlisted/'
                         '2026-10-21-testing-the-vmm-threshold-theory-step-by-step_files/')
plt.rcParams.update({
    "font.size": 12, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.8, "axes.axisbelow": True, "legend.frameon": False,
    "figure.dpi": 130, "savefig.bbox": "tight", "savefig.facecolor": "white"})
BASE_X = {**BASE, 4: 199, 7: 200, 9: 191}
d = load(); ro = d['r_out']
z = np.load(SCR + 'oos/run46_adc.npz'); H, C = z['H'], z['coinc']
rr = {x['sdt']: x for x in BENCH if x['gain'] == 4.5}
st = lambda v: 'P2_IN' if v < 8 else 'P2_MID' if v < 14 else 'P2_OUT'
edge = {v: 38 + 0.9 * (SDT46[v] - BASE_X[v]) for v in BASE_X}
# out-of-spill rate of "edge" hits only (within 32 ADC of the threshold edge) = the noise-like part
bench, beam, beam_edge, span = np.zeros((20, 64)), np.zeros((20, 64)), np.zeros((20, 64)), {}
for v in range(2, 20):
    s = min(rr, key=lambda k: abs(k - SDT46[v])); x = rr[s]; span[v] = x['span_s']
    bench[v] = np.array(x['per'][str(v)]['ch']) / x['span_s']
    b = int(edge.get(v, 50) // 16)
    frac = H[0, v, :, :b + 2].sum(1) / np.maximum(H[0, v].sum(1), 1)
    beam[v], beam_edge[v] = ro[v], ro[v] * frac
FB, FY = 0.3, 1.0 / d['t_out'] * 0.6      # where zeros are drawn
num = {}

fig, (ax, bx) = plt.subplots(1, 2, figsize=(14.5, 5.8), gridspec_kw=dict(width_ratios=[1.15, 1]))
for name, c in STC.items():
    vs = [v for v in range(2, 20) if st(v) == name]
    xb = np.concatenate([bench[v] for v in vs]); yb = np.concatenate([beam_edge[v] for v in vs])
    ax.scatter(np.maximum(xb, FB), np.maximum(yb, FY), s=22, color=c, alpha=0.75, edgecolor='white', lw=0.4, label=name)
lim = (FY * 0.7, 1e6)
ax.plot(lim, lim, color=MUTED, lw=1, ls='--'); ax.text(1.5e5, 2.5e4, 'beam = bench', fontsize=9.5, color=MUTED, ha='left')
ax.plot(lim, np.array(lim) * 100, color=MUTED, lw=1, ls=':'); ax.text(1.2e3, 5e5, 'beam = 100× bench', fontsize=9.5, color=MUTED, ha='right')
ax.axvspan(0.2, 0.42, color=GRID, alpha=0.7, lw=0); ax.text(0.47, 2e5, '← silent on bench', ha='left', fontsize=9, color=MUTED)
ax.axhspan(lim[0], FY * 1.5, color=GRID, alpha=0.7, lw=0)
lab = {(19, 20): (8, -2), (19, 1): (-8, 9), (19, 36): (8, -3), (14, 34): (-8, -4), (11, 10): (6, -12), (4, 40): (8, 0), (4, 39): (-8, 0),
       (4, 41): (8, -6), (10, 22): (8, 0), (10, 24): (8, -4)}
for (v, c), off in lab.items():
    ax.annotate(f'VMM {v} ch {c}', (max(bench[v, c], FB), max(beam_edge[v, c], FY)), xytext=off, textcoords='offset points',
                fontsize=8.5, color=INK, ha='right' if off[0] < 0 else 'left', va='center')
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(0.2, 2e6); ax.set_ylim(*lim)
ax.set_xlabel('bench, no beam (30 Jul), nearest sdt:  hits per second')
ax.set_ylabel('run_46 between spills, threshold-edge hits:  per second')
ax.legend(loc='lower right', fontsize=10, markerscale=1.5)
ax.set_title('a   same channels, but ~100× hotter in the beam hall', loc='left', fontsize=11.5, color=INK)

# b: ADC spectra of P2_OUT out-of-spill hits, hot channels vs the rest, vs in-spill
xa = np.arange(64) * 16 + 8
hot = [(v, c) for v in range(14, 20) for c in range(64) if beam_edge[v, c] > 100]
rest = np.ones((20, 64), bool); rest[:14] = False
for v, c in hot: rest[v, c] = False
def spec(j, sel):
    out = np.zeros(70)
    for v in range(14, 20):
        sh = int(round(edge[v] / 16)); h = H[j, v][sel[v]].sum(0)
        out[: 64 - sh] += h[sh:]
    return out
selh = {v: np.array([(v, c) in hot for c in range(64)]) for v in range(14, 20)}
selr = {v: rest[v] for v in range(14, 20)}
xr = np.arange(70) * 16
for lbl, sp, c, ls in (('between spills, the noisy channels (beam-edge rate > 100 Hz)', spec(0, selh), "#8a6fbf", '-'),
                        ('between spills, all other channels', spec(0, selr), STC['P2_OUT'], '-'),
                        ('in spill, all other channels (beam signal + the same floor)', spec(1, selr), MUTED, '--')):
    bx.step(xr, sp / sp.sum() / 16, where='post', color=c, lw=2, ls=ls, label=lbl)
bx.set_yscale('log'); bx.set_xlim(0, 900); bx.set_ylim(1e-6, 0.1)
bx.set_xlabel('VMM ADC above the threshold edge'); bx.set_ylabel('fraction of hits per ADC count')
bx.legend(loc='upper right', fontsize=9.2)
cp = C[0, 14:20][rest[14:20]].sum() / C[:, 14:20][:, rest[14:20]].sum()
ch_ = C[0, 14:20][~rest[14:20]].sum() / C[:, 14:20][:, ~rest[14:20]].sum()
bx.text(0.97, 0.52, f'hit on another P2_OUT channel within ±150 ns:\nnoisy channels {ch_*100:.0f} %,  other channels {cp*100:.0f} %',
        transform=bx.transAxes, ha='right', fontsize=9.5, color=INK)
bx.set_title('b   P2_OUT between spills: two kinds of hit', loc='left', fontsize=11.5, color=INK)
fig.tight_layout(); fig.savefig(OUT + 'fig14b_bench_vs_out_of_spill.png'); plt.close(fig)

num.update(t_out=float(d['t_out']), t_in=float(d['t_in']), partner_rest=float(cp), partner_hot=float(ch_),
           n_hot_out=len(hot),
           floor_hz={s: float(np.median([np.median(beam[v]) for v in range(2, 20) if st(v) == s])) for s in STC},
           bench_hot=int((bench[2:] > 10).sum()), bench_hot_beam_hot=int(((bench[2:] > 10) & (beam_edge[2:] > 10)).sum()),
           beam_hot=int((beam_edge[2:] > 100).sum()), beam_hot_bench_hot=int(((beam_edge[2:] > 100) & (bench[2:] > 10)).sum()))
json.dump(num, open(SCR + 'oos/oos_numbers.json', 'w'), indent=1)
print(json.dumps(num, indent=1))
