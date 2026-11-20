#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figures for the weekly-meeting deck update of 20 Nov 2026 ("Where the VMM inefficiency comes from"), in the
deck's palette. Inputs = the products of muon_funnel_summary.py, gain_scan_eff.py, muon_adc_gain.py,
threshold_spectrum.py, silent_vs_clipped.py, silent_late_hits.py, wrapup_losses.py and the bscan_target reports.

Usage (lxplus, env.sh): gain_deck_figs.py GS_OUT_DIR REPORT_DIR FUNNEL_SUMMARY OUT_DIR
"""
import json
import os
import re
import statistics as st
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

GS, REP, FUN, OUT = sys.argv[1:5]
os.makedirs(OUT, exist_ok=True)
RED, NAVY, TEAL, GREEN, GREY, DARK, AMBER = '#E50019', '#3E4A83', '#0093A7', '#2E8B3E', '#6B6B6B', '#262626', '#E8A33D'
GAINC = {3.0: NAVY, 4.5: TEAL, 6.0: GREEN, 9.0: RED}
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Liberation Sans', 'DejaVu Sans'],
                     'font.size': 13, 'axes.edgecolor': GREY, 'axes.labelcolor': DARK, 'xtick.color': GREY,
                     'ytick.color': GREY, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': True,
                     'grid.color': '#E6E6E6', 'legend.frameon': False, 'savefig.dpi': 200, 'savefig.bbox': 'tight',
                     'figure.facecolor': 'white'})
DETS = ['P2_IN', 'P2_MID', 'P2_OUT']
K4 = ['found', 'probe_silent', 'probe_far', 'probe_near_late']
KCOL = {'found': GREEN, 'probe_silent': RED, 'probe_far': AMBER, 'probe_near_late': NAVY}
KLAB = {'found': 'found', 'probe_silent': 'silent: no hit at all', 'probe_far': 'hit, but > 15 mm away',
        'probe_near_late': 'hit near, but late'}
FS = json.load(open(FUN))['groups']
GE = json.load(open(os.path.join(GS, 'gain_scan_eff.json')))
TS = json.load(open(os.path.join(GS, 'threshold_spectrum.json')))
MA = pd.read_csv(os.path.join(GS, 'muon_adc_gain.csv'))
nums = {}


def save(f, name):
    f.savefig(os.path.join(OUT, name))
    plt.close(f)


# ---------------------------------------------------------------- 1. the efficiency step: found / silent / far / late
f, ax = plt.subplots(figsize=(7.6, 3.4))
for i, d in enumerate(DETS[::-1]):
    on = FS['REF']['probes'][d]['on']
    n = sum(on.get(k, 0) for k in K4)
    x = 0
    for k in K4:
        w = on.get(k, 0) / n
        ax.barh(i, w * 100, left=x, color=KCOL[k], height=0.62, label=KLAB[k] if i == 0 else None)
        if w > 0.045:
            ax.text(x + w * 50, i, f'{100 * w:.0f} %', ha='center', va='center', color='white', fontsize=12,
                    fontweight='bold')
        x += w * 100
    nums[f'split_{d}'] = {k: on.get(k, 0) / n for k in K4} | {'n': n}
ax.set_yticks(range(3), DETS[::-1])
ax.set_xlim(0, 100); ax.set_xlabel('share of the muons expected on a live pad (%)')
ax.grid(axis='y', visible=False)
ax.legend(ncol=2, loc='upper center', bbox_to_anchor=(0.5, -0.28), fontsize=11)
save(f, 'c_split.png')

# ---------------------------------------------------------------- 2. budget: read-out, dead time, sparks vs silent
sil_mid = nums['split_P2_MID']['probe_silent'] * 100
f, ax = plt.subplots(figsize=(7.6, 3.2))
items = [('read-out (test pulses)', 0.14, NAVY), ('dead time', 0.29, NAVY), ('sparks', 0.12, NAVY),
         ('all three', 0.55, GREY), ('P2_MID silent', sil_mid, RED)]
for i, (lab, v, c) in enumerate(items[::-1]):
    ax.barh(i, v, color=c, height=0.6)
    ax.text(v + 0.2, i, f'{v:.2f}' if v < 1 else f'{v:.1f}', va='center', fontsize=12)
ax.set_yticks(range(len(items)), [x[0] for x in items[::-1]])
ax.set_xlabel('efficiency points'); ax.grid(axis='y', visible=False)
save(f, 'c_budget.png')

# ---------------------------------------------------------------- 3. the gain is real
f, ax = plt.subplots(figsize=(6.4, 4.4))
g3 = MA[MA.gain == 3.0].groupby('det').adc_median.mean()
for d, mk in zip(DETS, ('o', 's', '^')):
    t = MA[MA.det == d]
    rat = (t.adc_median - 40) / (g3[d] - 40)
    ax.scatter(t.gain, rat, marker=mk, s=60, color=[GAINC[g] for g in t.gain], edgecolor='white', zorder=3, label=d)
ax.plot([0, 10], [0, 10 / 3], color=GREY, ls='--', lw=1.5, label='if the gain is nominal')
ax.set_xlim(2, 10); ax.set_ylim(0, 3.6)
ax.set_xlabel('VMM gain setting (mV/fC)'); ax.set_ylabel('muon pulse height\nrelative to gain 3')
leg = ax.legend(fontsize=11, loc='upper left')
for h in leg.legend_handles[:3]:
    h.set_color(GREY)
save(f, 'c_gain_real.png')

# ---------------------------------------------------------------- 4. effect 1: same channel, more charge for the same noise
per = {}
for g in ('3.0', '6.0', '9.0'):
    R = json.load(open(os.path.join(REP, f'hvon_g{g}_R10_report.json')))
    for c, v in R.items():
        med = set(int(x) for x in (v.get('chip_median_used') or []))
        for k, x in (v.get('targets') or {}).items():
            if int(k) not in med:
                per.setdefault((int(c), int(k)), {})[float(g)] = x / float(g)
f, ax = plt.subplots(figsize=(6.6, 4.2))
bins = np.linspace(0.5, 3.0, 51)
for g, c in ((6.0, GREEN), (9.0, RED)):
    r = [v[g] / v[3.0] for v in per.values() if 3.0 in v and g in v]
    ax.hist(r, bins=bins, histtype='stepfilled', alpha=0.35, color=c)
    ax.hist(r, bins=bins, histtype='step', lw=2, color=c, label=f'gain {g:g} / gain 3: median {np.median(r):.2f}')
    nums[f'noise_ratio_g{g:g}'] = float(np.median(r))
ax.axvline(1, color=GREY, ls='--', lw=1.5)
ax.annotate('ratio 1 = same charge:\nmore gain would help', xy=(1.0, ax.get_ylim()[1] * 0.55), xytext=(2.05, ax.get_ylim()[1] * 0.55),
            color=GREY, fontsize=11, va='center', arrowprops=dict(arrowstyle='->', color=GREY))
ax.set_xlabel('threshold for the same noise rate (R = 10 Hz), in charge,\nsame channel, ratio to gain 3')
ax.set_ylabel('channels'); ax.legend(fontsize=11, loc='upper right')
save(f, 'c_noise_fC.png')

# ---------------------------------------------------------------- 5. effect 2: trim range per chip
f, ax = plt.subplots(figsize=(6.6, 4.2))
for i, g in enumerate(('3.0', '6.0', '9.0')):
    R = json.load(open(os.path.join(REP, f'hvon_g{g}_R10_report.json')))
    sp = []
    for c, v in R.items():
        med = set(int(x) for x in (v.get('chip_median_used') or []))
        m = [x for k, x in (v.get('targets') or {}).items() if int(k) not in med]
        if len(m) > 2:
            sp.append(max(m) - min(m))
    xs = i + (np.random.default_rng(1).random(len(sp)) - 0.5) * 0.35
    ax.scatter(xs, sp, s=45, color=GAINC[float(g)], edgecolor='white', zorder=3)
    nums[f'trim_spread_g{g}'] = float(np.median(sp))
ax.axhline(31 * 1.1, color=RED, lw=2, ls='--')
ax.text(-0.4, 31 * 1.1 + 2, 'what a trim can correct (31 x 1.1 mV)', color=RED, ha='left', fontsize=11)
ax.set_xticks(range(3), ['gain 3', 'gain 6', 'gain 9'])
ax.set_ylabel('spread of the needed thresholds\ninside one chip (mV)'); ax.grid(axis='x', visible=False)
save(f, 'c_trim_spread.png')

# ---------------------------------------------------------------- 6. clipped channels are not the cause
txt = open(os.path.join(GS, 'silent_vs_clipped.txt')).read().splitlines()
rows = []
for ln in txt:
    m = re.match(r'(gain \S+ R=\d+)\s+(clipped|no crossing \(chip median\)|normal)\s+muons\s+(\d+)\s+silent\s+([\d.]+) \+- ([\d.]+)', ln)
    if m:
        rows.append((m[1], m[2], int(m[3]), float(m[4]), float(m[5])))
C = pd.DataFrame(rows, columns=['set', 'cls', 'n', 's', 'e'])
f, ax = plt.subplots(figsize=(7.2, 4.2))
sets = ['gain 6.0 R=10', 'gain 9.0 R=2', 'gain 9.0 R=10', 'gain 9.0 R=50']
for j, (cls, c) in enumerate((('normal', NAVY), ('clipped', RED), ('no crossing (chip median)', GREY))):
    for i, s_ in enumerate(sets):
        r = C[(C.set == s_) & (C.cls == cls)]
        if len(r) and r.n.iloc[0] > 50:
            ax.bar(i + (j - 1) * 0.26, r.s.iloc[0], 0.24, yerr=r.e.iloc[0], color=c, capsize=3,
                   label=cls.replace(' (chip median)', ': chip median') if i == 1 else None)
ax.set_xticks(range(len(sets)), [s_.replace('.0', '').replace(' R=', '\nR = ') + ' Hz' for s_ in sets])
ax.set_ylabel('silent muons (%)'); ax.grid(axis='x', visible=False)
ax.legend(fontsize=11, loc='upper left')
save(f, 'c_clipped.png')

# ---------------------------------------------------------------- 7. together: silent vs threshold in charge
conf = pd.DataFrame(GE['configs'])
for d, w_, name in (('P2_MID', 8.4, 'c_silent_thr_mid.png'), ('ALL', 12.6, 'c_silent_thr_all.png')):
    dets = DETS if d == 'ALL' else [d]
    f, axs = plt.subplots(1, len(dets), figsize=(w_, 4.4), sharey=True, squeeze=False)
    for ax, dd in zip(axs[0], dets):
        t = conf[(conf.det == dd) & conf.thr_fC.notna()]
        ax.axvspan(3.0, 4.45, color=GREEN, alpha=0.08)
        ax.text(3.05, 29, 'reachable', color=GREEN, fontsize=11, va='top')
        for _, r in t.iterrows():
            e = np.sqrt(r.silent * (1 - r.silent) / r.n) * 100
            ax.errorbar(r.thr_fC, r.silent * 100, yerr=e, fmt='s' if r.family == 'REF' else 'o', ms=8,
                        color=GAINC[r.gain], mfc='white' if r.family == 'bscan' else GAINC[r.gain], mew=2,
                        capsize=0, zorder=3)
        ax.set_xlim(3.0, 8.6); ax.set_ylim(0, 30)
        ax.set_xlabel('threshold in charge (fC)')
        if len(dets) > 1:
            ax.set_title(dd, fontsize=13, color=DARK)
    axs[0][0].set_ylabel('silent muons (%)')
    for g, c in GAINC.items():
        axs[0][-1].plot([], [], 'o', color=c, label=f'gain {g:g} mV/fC')
    axs[0][-1].legend(fontsize=11, loc='lower right')
    save(f, name)

# ---------------------------------------------------------------- 8. efficiency over the gain x noise-rate grid (P2_MID)
f, ax = plt.subplots(figsize=(7.2, 4.2))
gains, Rs = [3.0, 4.5, 6.0, 9.0], [50, 10, 2]
M = np.full((3, 4), np.nan); A = [[''] * 4 for _ in range(3)]
for _, r in conf[conf.det == 'P2_MID'].iterrows():
    if r.R in Rs and r.gain in gains:
        i, j = Rs.index(int(r.R)), gains.index(r.gain)
        M[i, j] = 100 * (1 - r.silent)
        A[i][j] = f'{100 * (1 - r.silent):.1f} %\n{r.thr_fC:.1f} fC'
im = ax.imshow(M, cmap='RdYlGn', vmin=74, vmax=90, aspect='auto')
for i in range(3):
    for j in range(4):
        if A[i][j]:
            ax.text(j, i, A[i][j], ha='center', va='center', fontsize=11, color=DARK)
ax.set_xticks(range(4), [f'gain {g:g}' + ('\n(threshold scan)' if g == 4.5 else '') for g in gains])
ax.set_yticks(range(3), [f'R = {R} Hz' for R in Rs]); ax.grid(False)
cb = f.colorbar(im, ax=ax, fraction=0.04); cb.set_label('P2_MID: muons with any hit (%)')
save(f, 'c_grid.png')

# ---------------------------------------------------------------- 9. is the silent share a threshold tail?
f, ax = plt.subplots(figsize=(7.4, 4.4))
for g, c in (('gain 3.0', NAVY), ('REF g4.5', TEAL), ('gain 9.0', RED)):
    T = TS[f'{g}|P2_MID']
    e = np.array(T['density_rel_edges']); h = np.array(T['density_rel'], float)
    h = h / h.sum() / (e[1] - e[0])
    ax.step(e[:-1], h, where='post', color=c, lw=2, label=f'{g.replace("g4.5", "(gain 4.5)")}')
    ax.axvline(T['edge_over_signal'], color=c, ls=':', lw=1.5)
ax.set_xlim(0, 2.5)
ax.set_xlabel('largest pad signal / median muon signal (above baseline)')
ax.set_ylabel('muons (normalised)')
ax.legend(fontsize=11, title='dotted line = threshold edge', title_fontsize=10)
save(f, 'c_spectrum.png')

# ---------------------------------------------------------------- 10. late hits of the 'silent' muons
late = {}
cur = None
for ln in open(os.path.join(GS, 'silent_late_hits.txt')):
    m = re.match(r'(REF g4\.5|gain \S+)\s+(P2_\w+)\s+silent\s+(\d+): (.*)', ln)
    if m:
        cells = re.findall(r'(-?\d+)\.\.(-?\d+): ([\d.]+) % \(rand ([\d.]+)\)', m[4])
        late[(m[1], m[2])] = [(int(a), int(b), float(c) - float(r_)) for a, b, c, r_ in cells]
f, ax = plt.subplots(figsize=(7.4, 4.2))
labs = ['-1000..-100', '100..200', '200..300', '300..500', '500..1000']
for k, (g, d, c) in enumerate((('REF g4.5', 'P2_MID', TEAL), ('gain 9.0', 'P2_MID', RED), ('gain 9.0', 'P2_IN', '#F08A97'))):
    v = {f'{a}..{b}': x for a, b, x in late[(g, d)]}
    ax.bar(np.arange(5) + (k - 1) * 0.27, [v.get(l_, 0) for l_ in labs], 0.25, color=c,
           label=f'{d}, {g.replace("REF g4.5", "REF (gain 4.5)")}')
ax.set_xticks(range(5), ['earlier\n(-1 µs..-100 ns)', '100-200 ns\nlate', '200-300 ns', '300-500 ns', '0.5-1 µs'])
ax.set_ylabel('silent muons with a hit\nnear the prediction (%)'); ax.grid(axis='x', visible=False)
ax.legend(fontsize=11)
save(f, 'c_late.png')

# ---------------------------------------------------------------- 11. what the silent muons are made of (P2_MID)
grp = {'gain 3': ('gain 3.0', lambda r: r.gain == 3.0), 'REF (gain 4.5)': ('REF g4.5', lambda r: r.family == 'REF'),
       'gain 6': ('gain 6.0', lambda r: r.gain == 6.0 and r.R in (10, 50)), 'gain 9': ('gain 9.0', lambda r: r.gain == 9.0)}
parts = {}
for lab, (tk, sel) in grp.items():
    t = conf[(conf.det == 'P2_MID')]
    t = t[[bool(sel(r)) for _, r in t.iterrows()]]
    sil = 100 * (t.silent * t.n).sum() / t.n.sum()
    T = TS[f'{tk}|P2_MID']
    tail = 100 * 0.5 * (T['below_edge_lower'] + T['below_edge_upper'])
    tail_rng = (100 * T['below_edge_lower'], 100 * T['below_edge_upper'])
    lt = sum(x for a, b, x in late[(tk, 'P2_MID')] if (a, b) != (-2000, -1000) and a != 1000 and (a, b) != (-100, 100))
    lt_pts = sil * lt / 100
    parts[lab] = dict(silent=sil, tail=tail, tail_lo=tail_rng[0], tail_hi=tail_rng[1], late=lt_pts,
                      rest=max(sil - tail - lt_pts, 0))
nums['silent_parts'] = parts
f, ax = plt.subplots(figsize=(7.6, 4.4))
for i, (lab, p) in enumerate(parts.items()):
    ax.bar(i, p['tail'], 0.6, color=AMBER, label='signal below threshold (spectrum tail)' if i == 0 else None)
    ax.errorbar(i, p['tail'], yerr=[[p['tail'] - p['tail_lo']], [p['tail_hi'] - p['tail']]], color=DARK, capsize=4)
    ax.bar(i, p['late'], 0.6, bottom=p['tail'], color=NAVY, label='hit, but outside the time window' if i == 0 else None)
    ax.bar(i, p['rest'], 0.6, bottom=p['tail'] + p['late'], color=RED, label='no hit at all, not explained' if i == 0 else None)
    ax.text(i, p['silent'] + 0.5, f'{p["silent"]:.1f} %', ha='center', fontsize=12, fontweight='bold')
ax.set_xticks(range(len(parts)), list(parts))
ax.set_ylabel('P2_MID silent muons (%)'); ax.set_ylim(0, 30); ax.grid(axis='x', visible=False)
ax.legend(fontsize=11, loc='upper left')
save(f, 'c_silent_parts.png')

json.dump(nums, open(os.path.join(OUT, 'deck_numbers.json'), 'w'), indent=1)
print(json.dumps(nums, indent=1))
