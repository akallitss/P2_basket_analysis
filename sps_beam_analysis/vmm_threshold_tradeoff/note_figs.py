"""Figures for the 'testing the VMM threshold theory' note."""
import json, glob, os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/'
sys.path.insert(0, REPO + 'dream_vmm_comparison')
import threshold_model as tm

SCR = os.path.dirname(os.path.abspath(__file__)) + '/'
OUT = os.path.expanduser('~/Documents/PostDocSaclay/cern-site/notes/unlisted/'
                         '2026-10-21-testing-the-vmm-threshold-theory-step-by-step_files/')
os.makedirs(OUT, exist_ok=True)
LACIE = ('/media/ak271430/LaCie/Extras/Physics/Post-Doc-Saclay/data/SPS_beam_test/'
         'TB_July2026_H4/analysis/det3/')

VMM, DREAM, PRED = "#2a78d6", "#eb6834", "#3b4652"
INK, MUTED, GRID = "#17202a", "#5b6774", "#e3e7eb"
plt.rcParams.update({
    "font.size": 12, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.8, "axes.axisbelow": True, "legend.frameon": False,
    "figure.dpi": 130, "savefig.bbox": "tight", "savefig.facecolor": "white"})

# ------------------------------------------------------------------ data ----
g, H, binw = tm.load()
S = tm.Spectra(H, binw)
full = pd.read_csv(tm.DATA + '/compare_dream_vmm_P2_OUT.csv')
idx = full.index[full['use']].to_numpy()
effd = g.eff_d.to_numpy()
chip = g.vmm.to_numpy()
Ts = np.arange(4.0, 3000.0, 1.0)
F = np.array([S.frac_above(T) for T in Ts])


def vmm(rs):
    f = glob.glob(SCR + f'eff/vmm_eff_{rs}*.json')[0]
    s = [x for x in json.load(open(f))['stations'] if x['station'] == 'P2_OUT'][0]['per_pad']
    return np.array(s['n'])[idx].astype(float), np.array(s['k'])[idx].astype(float)


def chi_curve(n, k, m):
    e = k / np.maximum(n, 1)
    return ((n * (effd * F - e) ** 2)[:, m]).sum(1) / n[m].sum()


def fitT(n, k, m):
    c = chi_curve(n, k, m)
    return Ts[np.argmin(c)], c


def wmean(x, n, m):
    return (x[m] * n[m]).sum() / n[m].sum()


# ---------------------------------------------------- fig 1: the model -----
x = (np.arange(H.shape[1]) + 0.5) * binw
order = np.argsort(g.amp_med_d.to_numpy())
pads = {"strong pad": order[-4], "weak pad": order[3]}
fig, axs = plt.subplots(2, 1, figsize=(8.5, 5.6), sharex=True)
for ax, (lbl, j) in zip(axs, pads.items()):
    h = H[j] / H[j].sum() / binw
    h = np.convolve(h, np.ones(3) / 3, 'same')
    ax.plot(x, h, color=DREAM, lw=2)
    ax.fill_between(x, h, where=x >= 162, color=DREAM, alpha=0.18, lw=0)
    ax.fill_between(x, h, where=x <= 162, color=DREAM, alpha=0.75, lw=0)
    ax.axvline(27, color=DREAM, ls='--', lw=1.5)
    ax.axvline(162, color=VMM, lw=2.2)
    f = S.frac_above(162)[j]
    ax.text(0.99, 0.88, f"{lbl}\nVMM keeps {f*100:.0f} % of what DREAM sees",
            transform=ax.transAxes, ha='right', va='top', color=INK)
    ax.set_xscale('log'); ax.set_xlim(15, 3600); ax.set_yticks([])
    ax.set_ylabel('pulses (a.u.)')
axs[0].text(27, axs[0].get_ylim()[1] * 0.97, ' DREAM threshold\n ~27 ADC (5σ, per channel)',
            color=DREAM, va='top', fontsize=10.5)
axs[1].annotate('VMM threshold: 162 DREAM ADC\n(one DAC value for the whole chip)', xy=(162, axs[1].get_ylim()[1] * 0.55),
            xytext=(700, axs[1].get_ylim()[1] * 0.55), color=VMM, va='center', fontsize=10.5,
            arrowprops=dict(arrowstyle='->', color=VMM))
axs[1].set_xlabel('pulse height on the pad, measured by DREAM  [DREAM ADC, log scale]')
fig.savefig(OUT + 'fig1_model.png'); plt.close(fig)

# -------------------------------------------- fig 2: the fit on run_46 -----
n46, k46 = vmm('run_46_cfg_gain4.5_peaktime200')
allm = np.ones(len(g), bool)
T46, c46 = fitT(n46, k46, allm)
p46 = effd * S.frac_above(T46)
e46 = k46 / n46
fig, axs = plt.subplots(1, 2, figsize=(10.5, 4.3))
ax = axs[0]
ax.plot(Ts, np.sqrt(c46) * 100, color=PRED, lw=2)
ax.axvline(T46, color=VMM, lw=1.5)
ax.set_xlim(0, 400); ax.set_ylim(0, 25)
ax.set_xlabel('trial threshold T  [DREAM ADC]'); ax.set_ylabel('rms(predicted − measured)  [points]')
ax.set_title('scan the one free number', loc='left', color=INK, fontsize=12)
ax.text(T46 + 8, 20, f'best T = {T46:.0f}', color=VMM)
ax = axs[1]
ax.plot([0.45, 1], [0.45, 1], color=MUTED, lw=1, ls=':')
ax.scatter(p46, e46, s=36, color=VMM, edgecolor='white', lw=0.6, zorder=3)
ax.set_xlim(0.45, 1); ax.set_ylim(0.45, 1)
ax.set_xlabel('predicted = DREAM eff × fraction above T')
ax.set_ylabel('measured VMM efficiency')
ax.set_title(f'53 pads, r = {np.corrcoef(p46, e46)[0,1]:.2f}', loc='left', color=INK, fontsize=12)
fig.tight_layout(); fig.savefig(OUT + 'fig2_fit_run46.png'); plt.close(fig)

# ----------------------------------------------- fig 3: amplifier gain ------
not19 = chip != 19
pairs = [('run_48_cfg_gain4.5_peaktime200_opt', 'run_43_cfg_gain3.0_peaktime200_opt', 'run_48 → run_43'),
         ('run_52_cfg_gain4.5_peaktime200_opt', 'run_50_cfg_gain3.0_peaktime200_opt', 'run_52 → run_50')]
gain_rows = []
fig, axs = plt.subplots(1, 3, figsize=(14, 4.3), gridspec_kw=dict(width_ratios=[1.25, 1, 1]))
for i, (hi, lo, lbl) in enumerate(pairs):
    nh, kh = vmm(hi); nl, kl = vmm(lo)
    Th, ch = fitT(nh, kh, not19); Tl, cl = fitT(nl, kl, not19)
    j = np.argmin(abs(Ts - 1.5 * Th)); pl = effd * F[j]; el = kl / nl
    gain_rows.append(dict(pair=lbl, T45=Th, T30=Tl, ratio=Tl / Th,
                          pred=wmean(pl, nl, not19), meas=kl[not19].sum() / nl[not19].sum(),
                          r=np.corrcoef(pl[not19], el[not19])[0, 1], meas45=kh[not19].sum() / nh[not19].sum()))
    if i == 0:
        ax = axs[0]
        ax.plot(Ts, np.sqrt(ch) * 100, color=VMM, lw=2, label='gain 4.5 mV/fC (run_48)')
        ax.plot(Ts, np.sqrt(cl) * 100, color=VMM, lw=2, ls='--', label='gain 3.0 mV/fC (run_43)')
        ax.axvline(1.5 * Th, color=PRED, lw=1.2, ls=':')
        ax.annotate('', xy=(1.5 * Th, 3), xytext=(Th, 3),
                    arrowprops=dict(arrowstyle='->', color=PRED, lw=1.5))
        ax.text(1.25 * Th, 4, '× 1.5 expected', ha='center', color=PRED)
        ax.text(Th, 6.5, f'{Th:.0f}', color=VMM, ha='center'); ax.text(Tl, 8.5, f'{Tl:.0f}', color=VMM, ha='center')
        ax.set_xlim(0, 500); ax.set_ylim(0, 25)
        ax.set_xlabel('trial threshold T  [DREAM ADC]'); ax.set_ylabel('rms residual  [points]')
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2, fontsize=10)
        ax.set_title('fit each run on its own', loc='left', fontsize=12, color=INK)
    ax = axs[1 + i]
    ax.plot([0.3, 1], [0.3, 1], color=MUTED, lw=1, ls=':')
    ax.scatter(pl[not19], el[not19], s=34, color=VMM, edgecolor='white', lw=0.6, zorder=3)
    ax.set_xlim(0.3, 1); ax.set_ylim(0.3, 1)
    ax.set_xlabel('PREDICTED at gain 3.0 (no fit)'); ax.set_ylabel('measured at gain 3.0')
    ax.set_title(f'{lbl}: {gain_rows[-1]["pred"]:.3f} predicted, {gain_rows[-1]["meas"]:.3f} measured',
                 loc='left', fontsize=11.5, color=INK)
fig.tight_layout(); fig.savefig(OUT + 'fig3_gain_step.png'); plt.close(fig)

# ----------------------------------------------- fig 4: VMM 19 alone --------
n48, k48 = vmm('run_48_cfg_gain4.5_peaktime200_opt')
chips = [16, 17, 18, 19]
sdt46 = {c: 224 for c in chips}; sdt48 = {16: 230, 17: 230, 18: 230, 19: 265}
chip_rows = []
for c in chips:
    m = chip == c
    chip_rows.append(dict(chip=c, npad=int(m.sum()),
                          e46=k46[m].sum() / n46[m].sum(), e48=k48[m].sum() / n48[m].sum(),
                          T46=fitT(n46, k46, m)[0], T48=fitT(n48, k48, m)[0]))
cr = pd.DataFrame(chip_rows)
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
xx = np.arange(len(chips)); w = 0.36
for ax, (a, b, yl) in zip(axs, [('e46', 'e48', 'VMM efficiency'), ('T46', 'T48', 'fitted threshold  [DREAM ADC]')]):
    ax.bar(xx - w / 2 - 0.01, cr[a], w, color=VMM, alpha=0.45, label='run_46', edgecolor='white')
    ax.bar(xx + w / 2 + 0.01, cr[b], w, color=VMM, label='run_48', edgecolor='white')
    for i, c in enumerate(chips):
        ax.text(xx[i] - w / 2, cr[a][i], f'{sdt46[c]}', ha='center', va='bottom', fontsize=9.5, color=MUTED)
        ax.text(xx[i] + w / 2, cr[b][i], f'{sdt48[c]}', ha='center', va='bottom', fontsize=9.5,
                color=INK if c == 19 else MUTED, fontweight='bold' if c == 19 else 'normal')
    ax.set_xticks(xx, [f'VMM {c}\n({n} pads)' for c, n in zip(chips, cr.npad)])
    ax.set_ylabel(yl); ax.grid(axis='x', visible=False)
axs[0].set_ylim(0.5, 1.0); axs[0].legend(loc='upper left', ncol=2)
axs[0].set_title('numbers on bars = sdt DAC setting', loc='left', fontsize=11, color=MUTED)
axs[1].set_ylim(0, 320)
fig.tight_layout(); fig.savefig(OUT + 'fig4_vmm19.png'); plt.close(fig)

# ----------------------------------------------- fig 5: mesh scan ----------
spec = {}
for f in glob.glob(LACIE + 'beam_nominal_meshscan_1/*/20_beam_spectra/scan_row.json'):
    if 'nominal_01' in f: continue
    s = json.load(open(f))['spectrum']
    spec[int(s['hv'])] = (np.array(s['edges']), np.array(s['counts'], float))
te = pd.read_csv(LACIE + 'beam_nominal_meshscan_1/scan/22_tag_probe_efficiency/tag_probe_efficiency_spark_vetoed.csv')
effd_hv = {int(r.hv): r.eff for r in te.itertuples() if r.sub_run != 'nominal_01'}


def fa(e, c, T):
    cum = np.concatenate([[0], np.cumsum(c)]) / c.sum()
    return 1 - np.interp(T, e, cum)


cond = pd.read_csv(REPO + 'vmm_dream_matching/efficiency_conditions.csv')
v32 = cond[cond.run.isin(['run_32']) & cond['sub'].str.startswith('mesh')].set_index('mesh').eff_OUT
hvs = sorted(spec)
e450, c450 = spec[450]
Tm = Ts[np.argmin(abs(np.array([fa(e450, c450, t) for t in Ts]) - v32[450] / effd_hv[450]))]
pred = np.array([effd_hv[h] * fa(*spec[h], Tm) for h in hvs])
mesh_rows = [dict(hv=h, vmm=v32.get(h, np.nan), pred=p, dream=effd_hv[h]) for h, p in zip(hvs, pred)]

fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.6), gridspec_kw=dict(width_ratios=[1, 1.15]))
ax = axs[0]
for h, a in [(450, 1.0), (420, 0.65), (390, 0.35)]:
    e, c = spec[h]; xc = 0.5 * (e[1:] + e[:-1]); d = c / c.sum() / np.diff(np.log(e.clip(1)))
    ax.plot(xc, d, color=DREAM, lw=2, alpha=a, label=f'mesh {h} V')
ax.axvline(Tm, color=VMM, lw=2); ax.axvline(27, color=DREAM, ls='--', lw=1.3)
ax.annotate(f'VMM level (fitted at\n450 V only: {Tm:.0f})', xy=(Tm, ax.get_ylim()[1] * 0.62), xytext=(1300, ax.get_ylim()[1] * 0.62), color=VMM, va='center', fontsize=10, arrowprops=dict(arrowstyle='->', color=VMM))
ax.set_xscale('log'); ax.set_xlim(20, 8000); ax.set_yticks([])
ax.text(28, ax.get_ylim()[1] * 0.5, ' DREAM\n threshold', color=DREAM, fontsize=10)
ax.set_xlabel('DREAM pulse height  [ADC, log]'); ax.set_ylabel('pulses per log-interval (a.u.)')
ax.legend(loc='upper right', fontsize=10)
ax.set_title('DREAM\'s own spectrum at each mesh voltage', loc='left', fontsize=12, color=INK)
ax = axs[1]
ax.plot(hvs, [effd_hv[h] for h in hvs], 'o-', color=DREAM, lw=2, ms=8, label='DREAM measured')
ax.plot(hvs, pred, '-', color=PRED, lw=2, label='VMM predicted (DREAM spectrum cut at one level)')
ax.plot(hvs, [v32.get(h) for h in hvs], 'o', color=VMM, ms=9, mec='white', mew=1, label='VMM measured (run_32)', zorder=4)
ax.set_yscale('log'); ax.set_ylim(0.04, 1.3)
ax.set_yticks([0.05, 0.1, 0.2, 0.5, 1.0], ['5 %', '10 %', '20 %', '50 %', '100 %'])
ax.set_xlabel('P2_OUT mesh voltage  [V]'); ax.set_ylabel('efficiency (log)')
ax.legend(loc='lower right', fontsize=10)
ax.set_title('one level, seven voltages', loc='left', fontsize=12, color=INK)
fig.tight_layout(); fig.savefig(OUT + 'fig5_mesh_scan.png'); plt.close(fig)

json.dump(dict(T46=T46, gain=gain_rows, chips=chip_rows, Tmesh=Tm, mesh=mesh_rows),
          open(SCR + 'note_numbers.json', 'w'), indent=1, default=float)
print(json.dumps(dict(T46=T46, gain=gain_rows, chips=chip_rows, Tmesh=Tm, mesh=mesh_rows), indent=1, default=float))
