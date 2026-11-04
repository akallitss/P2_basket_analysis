#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Masked channels per VMM against the threshold, one point per run of a threshold scan.

Read from the applied chip config of each sub-run (run_N/p2b-config-*<sub>.txt, the file inventory_runs.py reads):
per VMM the threshold sdt = bank1spi2*64 + bank1spi1 and the channels with mask bit sm 1. K (threshold offset from
each chip's safe value, DAC) comes from the sub-run label (_Kp64_, _K0_, _Km8_ ...), as in summarize_scint_scan.py.

Writes OUT_STEM.json / OUT_STEM.csv (run, K, per VMM sdt, number and list of masked channels),
OUT_STEM_per_vmm.png (one panel per VMM, rows = chambers: masked channels and sdt against K) and
OUT_STEM_summary.png (masked channels per chamber against K, and the VMM x run table as a heat map).

Usage: mask_vs_threshold.py RUNS_DIR LABEL_PATTERN OUT_STEM [--runs 179-200] [--min-files 100]
   e.g. mask_vs_threshold.py /eos/.../runs adwk2 ~/wk2/masks
"""
import argparse
import csv
import glob
import json
import os
import re

import numpy as np

from inventory_runs import LINE, VMMS
from coincidence_check import CHAMBER

DETS = ((3, 'det3 P2_IN'), (2, 'det2 P2_MID'), (1, 'det1 P2_OUT'))


def k_of(label):
    m = re.search(r'_K([pm]?)(\d+)', label)
    return None if not m else (-1 if m[1] == 'm' else 1) * int(m[2])


def read_config(path):
    spi, sm = {}, {v: {} for v in VMMS}
    for line in open(path):
        m = LINE.search(line.split('#')[0])
        if not m:
            continue
        v = 2 * int(m[1]) + int(m[2])
        if v not in sm:
            continue
        reg, field, val = m[3], m[4], int(m[5])
        if field == 'sdt':
            spi[(v, reg)] = val
        elif field == 'sm':
            ch = re.findall(r'\d+', reg)
            if ch:
                sm[v][int(ch[-1])] = val      # last setting of a channel wins
    out = {}
    for v in VMMS:
        sdt = None
        if (v, 'bank1spi1') in spi and (v, 'bank1spi2') in spi:
            sdt = spi[(v, 'bank1spi2')] * 64 + spi[(v, 'bank1spi1')]
        out[v] = {'sdt': sdt, 'masked': sorted(c for c, x in sm[v].items() if x)}
    return out


def parse_range(s):
    a, b = s.split('-')
    return set(range(int(a), int(b) + 1))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('runs_dir')
    ap.add_argument('pattern')
    ap.add_argument('out_stem')
    ap.add_argument('--runs', default=None)
    ap.add_argument('--min-files', type=int, default=100,
                    help='a run with at least this many capture files is a full point; K values without one are '
                         'shown with their latest short run (open markers)')
    a = ap.parse_args()
    keep = parse_range(a.runs) if a.runs else None

    rows = []
    for sdir in glob.glob(os.path.join(a.runs_dir, 'run_*', 'cfg_*' + a.pattern + '*')):
        rdir, sub = os.path.dirname(sdir), os.path.basename(sdir)
        n = int(rdir.split('_')[-1])
        if keep and n not in keep:
            continue
        cfgs = glob.glob(os.path.join(rdir, 'p2b-config-*' + sub[len('cfg'):] + '.txt'))
        if not cfgs:
            print('no config for', sdir)
            continue
        c = read_config(cfgs[0])
        rows.append({'run': n, 'label': sub, 'K': k_of(sub),
                     'files': len(glob.glob(os.path.join(sdir, 'raw_daq_data', '*.pcapng'))),
                     'vmm': {v: c[v] for v in VMMS}})
    rows.sort(key=lambda r: r['run'])
    if not rows:
        raise SystemExit('no sub-runs matching ' + a.pattern)

    print(f'{"run":>4} {"K":>4} {"files":>5} | ' + ' '.join(f'{v:>3}' for v in VMMS) + ' | ' +
          ' '.join(f'det{d:<2}' for d, _ in DETS))
    for r in rows:
        nm = [len(r['vmm'][v]['masked']) for v in VMMS]
        per = [sum(len(r['vmm'][v]['masked']) for v in VMMS if CHAMBER[v] == d) for d, _ in DETS]
        print(f'{r["run"]:>4} {r["K"]!s:>4} {r["files"]:>5} | ' + ' '.join(f'{x:>3}' for x in nm) + ' | ' +
              ' '.join(f'{x:>5}' for x in per))
    json.dump(rows, open(a.out_stem + '.json', 'w'), indent=1)
    with open(a.out_stem + '.csv', 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['run', 'K', 'label', 'vmm', 'chamber', 'sdt', 'n_masked', 'masked_channels'])
        for r in rows:
            for v in VMMS:
                x = r['vmm'][v]
                w.writerow([r['run'], r['K'], r['label'], v, CHAMBER[v], x['sdt'], len(x['masked']),
                            ' '.join(map(str, x['masked']))])

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from event_display import style, INK, MUTED, COL

    # one point per K: every full run (repeats kept), else the latest short run of that K
    pts = []
    for k in sorted({r['K'] for r in rows if r['K'] is not None}, reverse=True):
        rk = [r for r in rows if r['K'] == k]
        full = [r for r in rk if r['files'] >= a.min_files]
        pts += full if full else [dict(rk[-1], short=True)]
    pts.sort(key=lambda r: (-r['K'], r['run']))
    Ks = np.array([r['K'] for r in pts], float)
    short = np.array([bool(r.get('short')) for r in pts])
    print('points:', ', '.join(f"run {r['run']} K {r['K']}{' (short)' if r.get('short') else ''}" for r in pts))

    # per VMM: masked channels (left axis) against K, actual sdt written next to each point
    fig, axs = plt.subplots(3, 6, figsize=(18, 9), sharex=True, sharey=True, constrained_layout=True)
    for i, (d, nm) in enumerate(DETS):
        for j, v in enumerate(sorted(u for u in VMMS if CHAMBER[u] == d)):
            ax = axs[i, j]
            style(ax)
            y = np.array([len(r['vmm'][v]['masked']) for r in pts])
            s = [r['vmm'][v]['sdt'] for r in pts]
            o = np.argsort(-Ks, kind='stable'); o = o[~short[o]]   # line through the full runs only
            ax.plot(Ks[o], y[o], '-', color=COL[d], lw=1.2, alpha=0.5)
            ax.plot(Ks[~short], y[~short], 'o', color=COL[d], ms=6, mec='white', mew=1)
            ax.plot(Ks[short], y[short], 'o', color=COL[d], ms=6, mfc='white', mew=1.2)
            for k, yy, ss in zip(Ks, y, s):
                ax.annotate(f'{ss}', (k, yy), textcoords='offset points', xytext=(0, 6), ha='center',
                            fontsize=6, color=MUTED)
            ax.set_title(f'VMM {v}  ({nm.split()[1]})', fontsize=9, color=INK)
            if j == 0:
                ax.set_ylabel('masked channels (of 64)', fontsize=8)
            if i == 2:
                ax.set_xlabel('threshold: safe + K [DAC]', fontsize=8)
    ax.invert_xaxis()
    top = max(len(r['vmm'][v]['masked']) for r in pts for v in VMMS)
    ax.set_ylim(-1, max(8, top * 1.25 + 2))
    fig.suptitle(f'{a.pattern}: masked channels per VMM against the threshold, one point per run '
                 '(grey numbers: applied sdt in DAC; open: short run only)', fontsize=11)
    fig.savefig(a.out_stem + '_per_vmm.png', dpi=120)

    # summary: per chamber against K, and the VMM x run table
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(16, 6), gridspec_kw={'width_ratios': [1, 1.6]}, constrained_layout=True)
    style(ax)
    o = np.argsort(-Ks, kind='stable'); o = o[~short[o]]   # line through the full runs only
    for d, nm in DETS:
        y = np.array([sum(len(r['vmm'][v]['masked']) for v in VMMS if CHAMBER[v] == d) for r in pts])
        ax.plot(Ks[o], y[o], '-', color=COL[d], lw=2)
        ax.plot(Ks[~short], y[~short], 'o', color=COL[d], ms=7, mec='white', mew=1, label=f'{nm} (384 ch)')
        ax.plot(Ks[short], y[short], 'o', color=COL[d], ms=7, mfc='white', mew=1.5)
        ax.annotate(nm.split()[0], (Ks[o][-1], y[o][-1]), textcoords='offset points', xytext=(6, 0),
                    va='center', fontsize=8, color=INK)
    ax.invert_xaxis()
    ax.set_xlabel('threshold: safe + K [DAC]', fontsize=9)
    ax.set_ylabel('masked channels per chamber', fontsize=9)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=8, frameon=False, loc='upper left', title='open: short run only (not a final mask set)',
              title_fontsize=7)
    ax.set_title('per chamber', fontsize=10, color=INK)

    M = np.array([[len(r['vmm'][v]['masked']) for r in pts] for v in VMMS])
    im = bx.imshow(M, aspect='auto', cmap='Blues', vmin=0, vmax=max(1, M.max()))
    for (i, j), x in np.ndenumerate(M):
        if x:
            bx.text(j, i, str(x), ha='center', va='center', fontsize=7,
                    color='white' if x > 0.6 * M.max() else INK)
    bx.set_yticks(range(len(VMMS)))
    bx.set_yticklabels([f'VMM {v} (det{CHAMBER[v]})' for v in VMMS], fontsize=7)
    bx.set_xticks(range(len(pts)))
    bx.set_xticklabels([f'K {r["K"]:+d}\nrun {r["run"]}' + ('*' if r.get('short') else '') for r in pts], fontsize=7)
    for b in (5.5, 11.5):
        bx.axhline(b, color='white', lw=2)
    bx.set_title('masked channels, VMM x run (* short run)', fontsize=10, color=INK)
    fig.colorbar(im, ax=bx, shrink=0.8, label='masked channels (of 64)')
    fig.suptitle(f'{a.pattern}: channel masks against the threshold', fontsize=11)
    fig.savefig(a.out_stem + '_summary.png', dpi=120)
    print('wrote', a.out_stem + '{.json,.csv,_per_vmm.png,_summary.png}')


if __name__ == '__main__':
    main()
