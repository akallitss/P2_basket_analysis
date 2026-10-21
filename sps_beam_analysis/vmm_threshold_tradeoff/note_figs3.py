"""Part 3 figures: per-channel noise map from the 29-30 Jul bench scan and the threshold/masking trade-off."""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from tradeoff import onsets, BASE, rows, eff_d, g, w, vm, ch, TRIM_DAC

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
JULY = {14: 230, 15: 230, 16: 230, 17: 230, 18: 230, 19: 265}
SDT = np.arange(210, 305, 5)
num = {}


def rate_map(gain, v):
    rr = sorted([r for r in rows if r['gain'] == gain], key=lambda r: r['sdt'])
    return np.array([r['sdt'] for r in rr]), np.array([np.array(r['per'][str(v)]['ch']) / r['span_s'] for r in rr])


# ------------------------------------- fig 13: per-channel noise map --------
fig, axs = plt.subplots(2, 3, figsize=(13.5, 7.2), sharex=True, sharey=True)
for ax, v in zip(axs.flat, range(14, 20)):
    s, R = rate_map(4.5, v)
    m = s <= 260
    im = ax.pcolormesh(np.arange(65) - 0.5, np.append(s[m], s[m][-1] + 5) - 2.5, np.where(R[m] > 0, R[m], np.nan),
                       cmap='Purples', norm=LogNorm(1, 1e4), shading='flat')
    if JULY[v] <= 260:
        ax.axhline(JULY[v], color=INK, lw=1.2, ls='--')
        ax.text(63, JULY[v] + 1, f'July {JULY[v]}', ha='right', va='bottom', fontsize=9.5, color=INK)
    else:
        ax.text(63, 261, f'July {JULY[v]} ↑ (off scale)', ha='right', va='top', fontsize=9.5, color=INK)
    nz = [(R[i] > 10).sum() for i in range(m.sum())]
    for i, si in enumerate(s[m]):
        if si <= 235:
            ax.text(66, si, f'{nz[i]}', va='center', fontsize=9, color=NOISE if nz[i] else MUTED)
    ax.set_title(f'VMM {v}   (baseline ≈ sdt {BASE[v]})', loc='left', fontsize=11.5, color=INK)
    ax.set_ylim(207.5, 262.5); ax.set_xlim(-0.5, 63.5); ax.grid(False)
for ax in axs[1]:
    ax.set_xlabel('channel')
for ax in axs[:, 0]:
    ax.set_ylabel('sdt  [DAC counts]')
axs[0, 2].text(70, 262, 'channels\n> 10 Hz', fontsize=9, color=NOISE, va='top')
cb = fig.colorbar(im, ax=axs, shrink=0.8, pad=0.06)
cb.set_label('noise hits per second, per channel')
fig.suptitle('P2_OUT chips, 4.5 mV/fC, no beam (30 July): which channel is noisy at which threshold', x=0.06, ha='left',
             fontsize=12.5, color=INK)
fig.savefig(OUT + 'fig13_channel_noise_map.png'); plt.close(fig)

# --------------------------- fig 14: bench hot channels = beam hot channels -
z = np.load(SCR + 'vmm_adc_all.npz')
avc = z['run_46__cfg_gain4.5_peaktime200__avc'].astype(float)  # all chips at sdt 224 except 5/12 (256)
fig, ax = plt.subplots(figsize=(8.5, 5.2))
mk = {14: 'o', 16: 's', 19: 'D', 11: '^', 10: 'v'}
col = {14: CHIPC[14], 16: CHIPC[16], 19: CHIPC[19], 11: PRED, 10: MUTED}
for v in [14, 16, 19, 11, 10]:
    s, R = rate_map(4.5, v)
    b = R[list(s).index(225)]
    cnt = avc[v].sum(1); rel = cnt / np.median(cnt)
    ax.scatter(np.maximum(b, 0.25), rel, s=34, marker=mk[v], color=col[v], alpha=0.8, edgecolor='white', lw=0.5,
               label=f'VMM {v}')
    for c in np.where(rel > 8)[0]:
        off = {(19, 1): (-8, 8, 'right'), (19, 36): (8, -8, 'left'), (14, 34): (-8, -6, 'right')}.get((v, int(c)), (6, 0, 'left'))
        ax.annotate(f'VMM {v} ch {c}' if (v, int(c)) in [(14, 34), (11, 10)] else f'ch {c}', (max(b[c], 0.25), rel[c]),
                    xytext=off[:2], textcoords='offset points', fontsize=9, color=col[v], va='center', ha=off[2])
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(0.18, 3e3)
ax.axvspan(0.18, 0.4, color=GRID, alpha=0.6, lw=0); ax.text(0.27, 1.5e-3, 'silent\non bench', ha='center', fontsize=9, color=MUTED)
ax.set_xlabel('bench, no beam, sdt 225:  noise hits per second')
ax.set_ylabel('beam run_46, sdt 224:  hits ÷ chip median')
ax.legend(loc='upper left', ncol=1, fontsize=10)
ax.set_title('the channels that are noisy on the bench are the ones that flood in the beam', loc='left', fontsize=11.5, color=INK)
fig.savefig(OUT + 'fig14_bench_vs_beam.png'); plt.close(fig)

# --------------------------------- fig 15: per-chip trade-off, no trim ------
on45, _ = onsets(4.5)
fig, axs = plt.subplots(2, 3, figsize=(13.5, 7.4), sharex=True, sharey=True)
best = {}
for ax, v in zip(axs.flat, range(14, 20)):
    sig = np.array([np.average(eff_d(s - BASE[v], 4.5), weights=w) for s in SDT])
    nm = np.array([(on45[v] > s).sum() for s in SDT])
    tot = sig * (1 - nm / 64)
    ax.axvspan(195, 207.5, color=GRID, alpha=0.7, lw=0)
    ax.text(201, 0.42, 'not\nscanned', ha='center', fontsize=9, color=MUTED)
    ax.plot(SDT, sig, color=CHIPC[v], lw=1.5, ls='--', label='signal only')
    ax.plot(SDT, tot, 'o-', color=CHIPC[v], lw=2.2, ms=5, label='after masking noisy channels')
    for s_, n_, t_ in zip(SDT, nm, tot):
        if n_ and s_ <= 235:
            ax.text(s_, t_ - 0.035, f'{n_}', ha='center', va='top', fontsize=9, color=NOISE)
    i = int(np.argmax(tot)); best[v] = int(SDT[i])
    ax.plot([SDT[i]], [tot[i]], '*', ms=17, color=INK, zorder=4)
    ax.text(SDT[i] + 1.5, tot[i] + 0.02, f'best {SDT[i]}: {tot[i]*100:.0f} %', fontsize=10, color=INK)
    jv = JULY[v]; jt = tot[list(SDT).index(jv)]
    ax.plot([jv], [jt], 'o', ms=10, mfc='white', mec=INK, mew=1.5, zorder=4)
    ax.text(jv - 1.5, jt - 0.06, f'July {jv}: {jt*100:.0f} %', fontsize=9.5, color=MUTED, ha='right' if jv > 250 else 'left')
    ax.set_title(f'VMM {v}   (baseline ≈ sdt {BASE[v]})', loc='left', fontsize=11.5, color=INK)
    ax.set_xlim(195, 272); ax.set_ylim(0.3, 1.0)
    num[f'chip{v}'] = dict(best=best[v], eff_best=float(tot[i]), eff_july=float(jt), nmask_best=int(nm[i]))
for ax in axs[1]:
    ax.set_xlabel('sdt  [DAC counts]')
for ax in axs[:, 0]:
    ax.set_ylabel('efficiency (model)')
axs[0, 0].legend(loc='lower right', fontsize=9.5)
fig.suptitle('lower the threshold until masking costs more than it gains  (4.5 mV/fC; purple numbers = channels masked)',
             x=0.06, ha='left', fontsize=12.5, color=INK)
fig.tight_layout(); fig.savefig(OUT + 'fig15_tradeoff_per_chip.png'); plt.close(fig)


# ----------------------------- fig 16: station, with and without trim -------
def station_eff(gain, d, pessimistic, R=10.0):
    on, _ = onsets(gain, R)
    e = np.empty(len(g)); nraise = 0; nunk = 0; nmask = 0
    for v in range(14, 20):
        St = BASE[v] + d
        meas = on[v] > 210
        req = np.where(meas, on[v], max(St, 210) if pessimistic else St)
        req = np.maximum(req, St)
        mk = (req - St > TRIM_DAC) | (on[v] >= 999)
        nraise += int((meas & (on[v] > St) & ~mk).sum()); nmask += int(mk.sum())
        nunk += int((~meas).sum()) if St < 210 else 0
        sel = np.where(vm == v)[0]
        cache = {r: eff_d(r - BASE[v], gain) for r in np.unique(req)}
        for k in sel:
            e[k] = 0 if mk[ch[k]] else cache[req[ch[k]]][k]
    return np.average(e, weights=w), nraise, nunk, nmask


dd = np.arange(4, 46, 1)
fig, axs = plt.subplots(2, 1, figsize=(11, 8), sharex=True, gridspec_kw=dict(height_ratios=[2.3, 1]))
ax = axs[0]
cur = {}
for gain, c, lab in [(4.5, VMM, '4.5 mV/fC'), (3.0, PRED, '3.0 mV/fC')]:
    opt = np.array([station_eff(gain, d, False)[0] for d in dd])
    pes = np.array([station_eff(gain, d, True)[0] for d in dd])
    ax.fill_between(dd * 0.9, pes, opt, color=c, alpha=0.15, lw=0)
    ax.plot(dd * 0.9, opt, color=c, lw=2.2, label=f'{lab}, trimmed: unscanned channels stay quiet')
    ax.plot(dd * 0.9, pes, color=c, lw=1.4, ls='--', label=f'{lab}, trimmed: unscanned channels held at sdt 210')
    cur[gain] = (opt, pes)
ax.axhline(0.956, color=DREAM, lw=1.5); ax.text(40.5, 0.959, 'DREAM 95.6 %', ha='right', va='bottom', color=DREAM)
# July and no-trim points (common-d axis not exact for them: place at P2_OUT track-weighted mean d)
def mean_d(cfg):
    return np.average([cfg[v] - BASE[v] for v in vm], weights=w) * 0.9
cfgJ = JULY; dJ = np.array([cfgJ[v] - BASE[v] for v in vm], float)
eJ = np.average(np.concatenate([[0]]) if False else [eff_d(cfgJ[v] - BASE[v], 4.5)[k] for k, v in enumerate(vm)], weights=w)
on, _ = onsets(4.5)
eB = []
for k, v in enumerate(vm):
    eB.append(0 if on[v][ch[k]] > best[v] else eff_d(best[v] - BASE[v], 4.5)[k])
eB = np.average(eB, weights=w)
mJ = np.average(g.eff_v, weights=w)
ax.plot([mean_d(cfgJ)], [eJ], 'o', ms=11, mfc='white', mec=INK, mew=1.6, zorder=5)
ax.text(mean_d(cfgJ) + 0.6, eJ - 0.012, f'July config: model {eJ*100:.1f} %, measured {mJ*100:.1f} %', fontsize=10, va='top')
ax.plot([mean_d(best)], [eB], '*', ms=18, color=INK, zorder=5)
ax.text(mean_d(best) + 0.8, eB + 0.006, f'per-chip best sdt, no trim: {eB*100:.1f} %', fontsize=10)
ax.set_ylim(0.72, 0.975); ax.set_ylabel('P2_OUT efficiency (model)')
ax.legend(loc='lower left', fontsize=9.5)
ax.set_title('one common threshold above baseline for every channel, noisy channels raised by the trim', loc='left',
             fontsize=11.5, color=INK)
ax = axs[1]
nr = np.array([station_eff(4.5, d, False)[1] for d in dd]); nu = np.array([station_eff(4.5, d, False)[2] for d in dd])
ax.bar(dd * 0.9, nr, width=0.75, color=NOISE, label='noisy on the bench → raised by trim')
ax.bar(dd * 0.9, nu, bottom=nr, width=0.75, color=GRID, label='below the scan: status unknown')
ax.set_ylabel('P2_OUT channels\n(of 384)'); ax.set_xlabel('threshold above baseline  [VMM ADC]   (≈ mV)')
sec = axs[0].secondary_xaxis('top', functions=(lambda a: a / (1.12 * 4.5), lambda q: q * 1.12 * 4.5))
sec.set_xlabel('in charge at 4.5 mV/fC  [fC]', color=MUTED)
ax.legend(loc='upper right', fontsize=9.5); ax.set_xlim(3, 41)
fig.tight_layout(); fig.savefig(OUT + 'fig16_station_tradeoff.png'); plt.close(fig)

for d in [10, 15, 20, 25]:
    num[f'trim45_d{d}'] = [float(x) for x in station_eff(4.5, d, False)[:1] + station_eff(4.5, d, True)[:1]]
    num[f'trim30_d{d}'] = [float(x) for x in station_eff(3.0, d, False)[:1] + station_eff(3.0, d, True)[:1]]
    num[f'nraise45_d{d}'] = station_eff(4.5, d, False)[1]
num.update(eJ=float(eJ), mJ=float(mJ), eB=float(eB), best=best)
json.dump(num, open(SCR + 'part3_numbers.json', 'w'), indent=1)
print(json.dumps(num, indent=1))
