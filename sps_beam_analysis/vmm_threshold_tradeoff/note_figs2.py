"""Part 2 figures: what sdt is, sdt -> mV -> VMM ADC -> fC, noise floor, trim."""
import json, glob, os, re, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import uproot

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/'
sys.path.insert(0, REPO + 'dream_vmm_comparison')
import threshold_model as tm

SCR = os.path.dirname(os.path.abspath(__file__)) + '/'
OUT = os.path.expanduser('~/Documents/PostDocSaclay/cern-site/notes/unlisted/'
                         '2026-10-21-testing-the-vmm-threshold-theory-step-by-step_files/')
VMM, DREAM, PRED = "#2a78d6", "#eb6834", "#3b4652"
INK, MUTED, GRID, NOISE = "#17202a", "#5b6774", "#e3e7eb", "#8a6fbf"
CHIPC = {14: "#2a78d6", 15: "#1baf7a", 16: "#e0a31b", 17: "#c43a7a", 18: "#6b4fc9", 19: "#eb6834"}
plt.rcParams.update({
    "font.size": 12, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.8, "axes.axisbelow": True, "legend.frameon": False,
    "figure.dpi": 130, "savefig.bbox": "tight", "savefig.facecolor": "white"})

DAC_MV, DAC_OFF = 0.801, 25.182      # vmmsc ThresholdDAC_to_mV
ADC_PER_DAC = 0.9                    # measured (VMM19 0.85, chips 14-18 ~0.93)
BASE_ADC = 38.0                      # edge = 37.7 + 0.196 T
ADC_PER_MV = ADC_PER_DAC / DAC_MV    # ~1.12 ; Nov pulser gives 1.07
EDGE_PER_T = 0.196                   # VMM ADC per DREAM ADC at 4.5 mV/fC

z = np.load(SCR + 'vmm_adc_all.npz')


def edge(h, q=0.01):
    c = np.cumsum(h) / h.sum()
    return int(np.searchsorted(c, q))


# ------------------------------------------------ fig 6: the pulse cartoon --
# VMM shaper impulse response from the Iakovidis et al. note (3rd order, 1 real + 2 cc poles)
def vmm_pulse(t, tp):
    a = tp / 1.5
    p0, p1 = 1.263 / a, (1.149 - 0.789j) / a
    K1, K2 = 1.584, -0.792 - 0.115j
    f = K1 * np.exp(-t * p0) + 2 * np.abs(K2) * np.exp(-t * p1.real) * np.cos(-t * p1.imag + np.angle(K2))
    f = np.where(t > 0, f, 0.0)
    return f / f.max()


t = np.linspace(-100, 900, 2000)
base_mv = 170.0
thr_mv = DAC_MV * 230 + DAC_OFF            # sdt 230 -> 209 mV
fig, ax = plt.subplots(figsize=(10.5, 4.8))
rng = np.random.default_rng(3)
noise = np.convolve(rng.normal(0, 1.0, t.size), np.ones(25) / 25 * 5, 'same') * 1.2
ax.plot(t, base_mv + noise, color=MUTED, lw=1, alpha=0.7)
for t0, A, lab in [(80, 75, 'big pulse: fires'), (420, 30, 'small pulse: never fires')]:
    y = base_mv + A * vmm_pulse(t - t0, 200) + noise
    fired = y.max() > thr_mv
    ax.plot(t, y, color=VMM if fired else DREAM, lw=2.2)
    k = np.argmax(y)
    ax.annotate(lab, xy=(t[k], y[k]), xytext=((t[k] + 40, y[k] + 8) if fired else (190, 193)), color=VMM if fired else DREAM, fontsize=11,
                arrowprops=dict(arrowstyle='-', color=MUTED, lw=0.8))
    if fired:
        ax.plot([t[k]], [y[k]], 'o', color=VMM, ms=8)
        ax.annotate('peak is digitised\n→ VMM ADC value', xy=(t[k], y[k]), xytext=(t[k] + 150, y[k] - 5),
                    color=VMM, fontsize=10.5, arrowprops=dict(arrowstyle='->', color=VMM))
ax.axhline(thr_mv, color=INK, lw=1.8, ls='--')
ax.text(880, thr_mv + 1.5, f'threshold = 0.801 × sdt + 25.2 mV\n(sdt 230 → {thr_mv:.0f} mV)', ha='right', va='bottom', fontsize=10.5)
ax.axhline(base_mv, color=MUTED, lw=1)
ax.text(880, base_mv - 2, 'baseline ("pedestal") ≈ 170 mV, different on each channel', ha='right', va='top', fontsize=10.5, color=MUTED)
ax.annotate('', xy=(-60, thr_mv), xytext=(-60, base_mv), arrowprops=dict(arrowstyle='<->', color=INK, lw=1.4))
ax.text(-52, (thr_mv + base_mv) / 2, ' what matters:\n threshold − baseline', va='center', fontsize=10.5)
ax.set_xlabel('time  [ns]'); ax.set_ylabel('shaper output  [mV]')
ax.set_xlim(-100, 900); ax.set_ylim(155, 260)
ax.set_title('Inside one VMM channel (pulse shape: VMM shaper transfer function, 200 ns peaking)', loc='left', fontsize=12, color=INK)
fig.savefig(OUT + 'fig6_pulse_cartoon.png'); plt.close(fig)

# ------------------------------------------- fig 7: ADC spectra and edges --
runs = [('run_46__cfg_gain4.5_peaktime200', 'run_46 · 4.5 mV/fC · sdt 224', '#9bb8e0', '-'),
        ('run_52__cfg_gain4.5_peaktime200_opt', 'run_52 · 4.5 mV/fC · sdt 230 (VMM19 265)', VMM, '-'),
        ('run_50__cfg_gain3.0_peaktime200_opt', 'run_50 · 3.0 mV/fC · sdt 230', VMM, ':'),
        ('run_53__cfg_gain4.5_peaktime200_deflt', 'run_53 · 4.5 mV/fC · default (sdt 300)', INK, '-')]
fig, axs = plt.subplots(1, 3, figsize=(15, 4.4), sharey=True)
for ax, v in zip(axs, [16, 18, 19]):
    for k, lab, c, ls in runs:
        h = z[k + '__adc'][v].astype(float); h16 = h.reshape(-1, 8).sum(1)
        ax.step(np.arange(128) * 8 + 4, h16 / h16.sum(), where='mid', color=c, ls=ls, lw=1.8, label=lab)
        e = edge(h)
        ax.plot([e], [3e-5], 'v', color=c, ms=8)
    ax.axvline(BASE_ADC, color=MUTED, lw=1)
    ax.text(BASE_ADC + 2, 0.4, 'baseline\n≈ 38', color=MUTED, fontsize=9.5)
    ax.set_yscale('log'); ax.set_xlim(0, 400); ax.set_ylim(2e-5, 1)
    ax.set_title(f'VMM {v}', loc='left', fontsize=12, color=INK); ax.set_xlabel('VMM ADC')
axs[0].set_ylabel('fraction of hits (log)')
axs[0].legend(loc='upper center', bbox_to_anchor=(1.6, -0.18), ncol=2, fontsize=10)
fig.savefig(OUT + 'fig7_adc_edges.png'); plt.close(fig)

# ------------------------------- fig 8: edge vs sdt, edge vs fitted T ------
nn = json.load(open(SCR + 'note_numbers.json'))
T46 = {c['chip']: c['T46'] for c in nn['chips']}; T48 = {c['chip']: c['T48'] for c in nn['chips']}
E = {}
for k, s in [('run_46__cfg_gain4.5_peaktime200', {v: 224 for v in range(14, 20)}),
             ('run_48__cfg_gain4.5_peaktime200_opt', {**{v: 230 for v in range(14, 19)}, 19: 265}),
             ('run_50__cfg_gain3.0_peaktime200_opt', {v: 230 for v in range(14, 20)}),
             ('run_53__cfg_gain4.5_peaktime200_deflt', {v: 300 for v in range(14, 20)}),
             ('run_51__cfg_gain3.0_peaktime200_deflt', {v: 300 for v in range(14, 20)})]:
    for v in range(14, 20):
        E.setdefault(v, []).append((s[v], edge(z[k + '__adc'][v]), 'g3' if 'gain3.0' in k else 'g45'))
fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.6))
ax = axs[0]
for v in range(14, 20):
    pts = np.array([(a, b) for a, b, _ in E[v]], float)
    g3 = np.array([gg == 'g3' for _, _, gg in E[v]])
    ax.plot(pts[~g3, 0], pts[~g3, 1], 'o', color=CHIPC[v], ms=8, label=f'VMM {v}')
    ax.plot(pts[g3, 0], pts[g3, 1], 's', mfc='none', mec=CHIPC[v], ms=10, mew=1.5)
    p = np.polyfit(pts[:, 0], pts[:, 1], 1)
    xx = np.array([215, 305]); ax.plot(xx, np.polyval(p, xx), color=CHIPC[v], lw=1, alpha=0.6)
ax.axhline(BASE_ADC, color=MUTED, lw=1); ax.text(302, BASE_ADC + 2, 'baseline', ha='right', color=MUTED, fontsize=10)
ax.set_xlabel('sdt  [DAC counts]  (top axis: threshold in mV)'); ax.set_ylabel('measured edge  [VMM ADC]')
ax.set_title('edge moves ~0.9 ADC per DAC count; ○ = 4.5, □ = 3.0 mV/fC', loc='left', fontsize=11.5, color=INK)
ax.legend(ncol=3, fontsize=9.5, loc='upper left')
sec = ax.secondary_xaxis('top', functions=(lambda d: DAC_MV * d + DAC_OFF, lambda m: (m - DAC_OFF) / DAC_MV))
sec.set_xlabel('threshold  [mV, absolute]', color=MUTED)
ax = axs[1]
x, y, c = [], [], []
for v in [16, 17, 18, 19]:
    for T, run in [(T46[v], 'run_46__cfg_gain4.5_peaktime200'), (T48[v], 'run_48__cfg_gain4.5_peaktime200_opt')]:
        x.append(T); y.append(edge(z[run + '__adc'][v])); c.append(CHIPC[v])
p = np.polyfit(x, y, 1)
ax.scatter(x, y, c=c, s=70, edgecolor='white', zorder=3)
xx = np.array([0, 300]); ax.plot(xx, np.polyval(p, xx), color=PRED, lw=1.5)
ax.set_xlim(0, 300); ax.set_ylim(30, 100)
ax.set_xlabel('threshold fitted from efficiency (DREAM spectra)  [DREAM ADC]')
ax.set_ylabel('edge seen in the VMM\'s own ADC  [VMM ADC]')
ax.set_title(f'two independent measurements: edge = {p[1]:.0f} + {p[0]:.3f}·T,  r = {np.corrcoef(x, y)[0,1]:.3f}', loc='left', fontsize=11.5, color=INK)
ax.text(5, p[1] + 2, '← intercept = baseline', color=PRED, fontsize=10)
fig.tight_layout(); fig.savefig(OUT + 'fig8_edge_calibration.png'); plt.close(fig)

# --------------------------------------- fig 9: bench noise scan ------------
rows = json.load(open(SCR + 'thrscan_summary.json'))
fig, axs = plt.subplots(1, 2, figsize=(13, 4.5), sharey=True)
for ax, gsel in zip(axs, [3.0, 4.5]):
    rr = sorted([r for r in rows if r['gain'] == gsel], key=lambda r: r['sdt'])
    for v in range(14, 20):
        s = [r['sdt'] for r in rr]; rate = [max(r['per'][str(v)]['n'] / r['span_s'], 0.3) for r in rr]
        ax.plot(s, rate, 'o-', color=CHIPC[v], ms=5, lw=1.5, label=f'VMM {v}')
    ax.set_yscale('log'); ax.set_ylim(0.2, 1e5); ax.set_xlim(207, 262)
    ax.axhline(0.3, color=MUTED, lw=0.8, ls=':'); ax.text(260, 0.35, 'nothing seen', ha='right', fontsize=9, color=MUTED)
    ax.axvline(230, color=INK, lw=1, ls='--'); ax.text(231, 3e4, 'beam setting 230', fontsize=9.5)
    ax.set_xlabel('sdt  [DAC counts]')
    ax.set_title(f'bench, no beam, {gsel} mV/fC', loc='left', fontsize=12, color=INK)
axs[0].set_ylabel('noise hits per second, whole chip (64 ch)')
axs[1].legend(ncol=2, fontsize=9.5)
fig.tight_layout(); fig.savefig(OUT + 'fig9_bench_noise_scan.png'); plt.close(fig)

# ------------------------------------ fig 10: Nov trim calibration ----------
d = pd.read_csv(SCR + 'Threshold_093415.csv', skipinitialspace=True, header=None, skiprows=1,
                names=['hyb', 'fec', 'vmm', 'ch', 'sd', 'thr', 'sd_cal', 'thr_cal', 'best'])
d = d[d.hyb != 'hybrid_id']
for col in ['fec', 'vmm', 'ch', 'sd', 'thr', 'sd_cal', 'thr_cal', 'best']:
    d[col] = pd.to_numeric(d[col])
d = d[d.thr > 0]
q = d[(d.vmm == 14)]
fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.5))
ax = axs[0]
for ch, qq in q.groupby('ch'):
    ax.plot(qq.sd, qq.thr, color=VMM, lw=0.8, alpha=0.35)
    ax.plot([qq.sd_cal.iloc[0]], [qq.thr_cal.iloc[0]], 'o', color=INK, ms=3)
ax.axhline(q.best.iloc[0], color=DREAM, lw=1.8)
ax.text(31, q.best.iloc[0] + 1.2, f'common target {q.best.iloc[0]:.0f} mV', ha='right', color=DREAM)
ax.set_xlabel('trim  sd  (0 … 31)'); ax.set_ylabel('threshold above baseline  [mV]')
ax.set_title('one chip (Nov 2025 bench numbering VMM 14, 3 mV/fC): 64 channels', loc='left', fontsize=11.5, color=INK)
ax = axs[1]
first = d[d.sd == 0].groupby(['hyb', 'vmm', 'ch']).thr.first()
cal = d.groupby(['hyb', 'vmm', 'ch']).thr_cal.first()
first = first[first > 20]; cal = cal[cal > 20]
bins = np.arange(30, 80, 1)
ax.hist(first - first.groupby(level=[0, 1]).transform('mean'), bins=np.arange(-20, 21, 1), color=DREAM, alpha=0.6, label=f'no trim (sd = 0): σ = {first.groupby(level=[0,1]).std().mean():.1f} mV')
ax.hist(cal - cal.groupby(level=[0, 1]).transform('mean'), bins=np.arange(-20, 21, 1), color=VMM, alpha=0.8, label=f'after trim: σ = {cal.groupby(level=[0,1]).std().median():.1f} mV')
ax.set_xlabel('channel threshold − chip mean  [mV]'); ax.set_ylabel('channels (10 chips)')
ax.legend(fontsize=10, loc='upper left')
ax.set_title('trim removes the channel-to-channel spread', loc='left', fontsize=11.5, color=INK)
fig.tight_layout(); fig.savefig(OUT + 'fig10_trim.png'); plt.close(fig)

# ------------------------- fig 11: efficiency vs threshold above baseline ---
g, H, b = tm.load(); S = tm.Spectra(H, b); w = g.n_track_v.values.astype(float); ed = g.eff_d.values
aa = np.arange(4, 60, 0.5)
eff = np.array([np.average(ed * S.frac_above(a / EDGE_PER_T), weights=w) for a in aa])
p10 = np.array([np.percentile(ed * S.frac_above(a / EDGE_PER_T), 10) for a in aa])
now = {v: edge(z['run_52__cfg_gain4.5_peaktime200_opt__adc'][v]) - BASE_ADC for v in range(14, 20)}
fig, ax = plt.subplots(figsize=(11, 4.8))
ax.axvspan(0, 12, color=NOISE, alpha=0.12, lw=0)
ax.axvspan(12, 20, color=NOISE, alpha=0.06, lw=0)
ax.text(6, 0.62, 'noise floor\n(bench: chips\nfill with noise)', ha='center', color=NOISE, fontsize=10)
ax.text(16, 0.62, 'only a few\nnoisy channels', ha='center', color=NOISE, fontsize=10)
ax.plot(aa, eff, color=VMM, lw=2.5, label='P2_OUT station efficiency (model)')
ax.plot(aa, p10, color=VMM, lw=1.5, ls='--', label='worst 10 % of pads')
ax.axhline(0.956, color=DREAM, lw=1.5); ax.text(59, 0.962, 'DREAM 95.6 %', ha='right', color=DREAM)
for v, a in now.items():
    ax.plot([a], [np.interp(a, aa, eff)], 'o', color=CHIPC[v], ms=10, mec='white', zorder=4)
    ax.text(a, np.interp(a, aa, eff) + 0.025, f'{v}', ha='center', color=CHIPC[v], fontsize=10, fontweight='bold')
ax.set_xlim(0, 60); ax.set_ylim(0.5, 1.0)
ax.set_xlabel('threshold above baseline  [VMM ADC]   (≈ mV)')
ax.set_ylabel('efficiency')
sec = ax.secondary_xaxis('top', functions=(lambda a: a / (ADC_PER_MV * 4.5), lambda q_: q_ * ADC_PER_MV * 4.5))
sec.set_xlabel('threshold in charge at 4.5 mV/fC  [fC]', color=MUTED)
ax.legend(loc='upper right', bbox_to_anchor=(1, 0.93), fontsize=10)
ax.set_title('dots: where each P2_OUT chip sits at sdt 230 (VMM 19 at 265)', loc='left', fontsize=11.5, color=INK)
fig.savefig(OUT + 'fig11_eff_vs_threshold.png'); plt.close(fig)

# ---------------------------------------- fig 12: rate vs distance (halo) ---
gg = pd.read_csv(REPO + 'dream_vmm_comparison/data/compare_dream_vmm_P2_OUT.csv')
wt = gg.n_track_v.values; cx = (gg.x_v * wt).sum() / wt.sum(); cy = (gg.y_v * wt).sum() / wt.sum()
gg['r'] = np.hypot(gg.x_v - cx, gg.y_v - cy) / 10.0


def ntr(rs):
    d_ = json.load(open(glob.glob(SCR + f'eff/vmm_eff_{rs}*.json')[0]))
    return [x for x in d_['stations'] if x['station'] == 'P2_OUT'][0]['n_tracks_in_captures']


fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.4))
for rs, lab, c in [('run_52__cfg_gain4.5_peaktime200_opt', 'sdt 230', VMM), ('run_53__cfg_gain4.5_peaktime200_deflt', 'sdt 300 (default)', INK)]:
    avc = z[rs + '__avc'].astype(float); n = ntr(rs.replace('__', '_'))
    h = np.array([avc[int(v), int(c_)].sum() for v, c_ in zip(gg.vmm, gg.ch)]) / n * 1e3
    m = h > 0
    axs[0].scatter(gg.r[m], h[m], s=10, color=c, alpha=0.5, label=lab)
axs[0].set_yscale('log'); axs[0].set_xlabel('pad distance from beam centre  [cm]')
axs[0].set_ylabel('hits per 1000 reference tracks')
axs[0].axvspan(0, 5, color=DREAM, alpha=0.08); axs[0].text(2.5, 2, 'beam\nspot', ha='center', color=DREAM, fontsize=10)
axs[0].legend(fontsize=10); axs[0].set_title('rate falls smoothly with distance, same at sdt 230 and 300', loc='left', fontsize=11.5, color=INK)
far = gg[gg.r > 20]
for rs, lab, c in [('run_52__cfg_gain4.5_peaktime200_opt', 'far pads (>20 cm), sdt 230', VMM), ('run_53__cfg_gain4.5_peaktime200_deflt', 'far pads (>20 cm), sdt 300', INK)]:
    avc = z[rs + '__avc'].astype(float)
    hs = np.zeros(128)
    for v, c_ in zip(far.vmm, far.ch): hs += avc[int(v), int(c_)]
    axs[1].step(np.arange(128) * 8 + 4, hs / hs.sum(), where='mid', color=c, lw=1.8, label=lab)
axs[1].set_yscale('log'); axs[1].set_xlim(0, 800); axs[1].set_xlabel('VMM ADC'); axs[1].set_ylabel('fraction of hits')
axs[1].legend(fontsize=10); axs[1].set_title('far pads: a muon-like spectrum, cut lower at 230', loc='left', fontsize=11.5, color=INK)
fig.tight_layout(); fig.savefig(OUT + 'fig12_halo.png'); plt.close(fig)
print('now (ADC above baseline):', now)
