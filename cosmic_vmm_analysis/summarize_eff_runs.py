#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Table and figure of the per-run efficiencies (run_efficiency_batch.sh output), one point per run, never pooled.

Reads DIR/eff_runN.json (efficiency_tp.py) and DIR/str_runN.json (straight_track_rate.py). The threshold setting is the
sub-run name (safehvon_v1 = safe, thrhvon_p8 = safe + 8 DAC, ...). Left: efficiency of each chamber against run number,
coloured by setting; right: against the setting (safe + K DAC), each run its own point.

Usage: summarize_eff_runs.py DIR OUT_STEM [--probe-r 15]
"""
import argparse
import glob
import json
import os
import re
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from event_display import style, INK, MUTED, COL  # noqa: E402

PROBES = (('probe_det3', 3, 'det3 P2_IN'), ('probe_det2', 2, 'det2 P2_MID'), ('probe_det1', 1, 'det1 P2_OUT'))
SET_COL = {-8: '#7a3e9d', -6: '#d55181', -4: '#eb6834', 0: INK, 8: '#2a78d6', 16: '#1baf7a'}


def offset(title):
    m = re.search(r'thrhvon_([pm])(\d+)', title)
    if m:
        return (1 if m[1] == 'p' else -1) * int(m[2])
    return 0 if 'safe' in title else None


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('dir')
    ap.add_argument('out_stem')
    ap.add_argument('--probe-r', default='15')
    a = ap.parse_args()
    rows = []
    for f in sorted(glob.glob(os.path.join(a.dir, 'eff_run*.json')), key=lambda p: int(re.findall(r'\d+', p)[-1])):
        E = json.load(open(f))
        n = int(re.findall(r'\d+', os.path.basename(f))[0])
        title = E.get('title', '')
        k = offset(title)
        s = os.path.join(a.dir, f'str_run{n}.json')
        st = json.load(open(s))['runs'].get(f'run_{n}', {}) if os.path.exists(s) else {}
        row = {'run': n, 'setting': title.split(' ', 1)[-1], 'K': k, 'minutes': E['span_s'] / 60,
               'straight_per_h': st.get('signal_chip_per_h'), 'straight_err': st.get('signal_chip_per_h_err')}
        for key, d, _ in PROBES:
            r = E[key]['vs_probe_r'][a.probe_r]
            row[f'eff{d}'], row[f'err{d}'], row[f'tags{d}'] = r.get('eff'), r.get('err'), r.get('n_sig')
        rows.append(row)
    print(f'{"run":>4} {"setting":24} {"K":>4} {"min":>5} {"straight/h":>12} ' +
          ' '.join(f'{"eff det" + str(d):>16}' for _, d, _ in PROBES))
    for r in rows:
        sp = f'{r["straight_per_h"]:.0f}+-{r["straight_err"]:.0f}' if r['straight_per_h'] is not None else '-'
        effs = ' '.join(f'{r[f"eff{d}"]:6.3f}+-{r[f"err{d}"]:.3f} ({r[f"tags{d}"]:4.0f})' if r[f'eff{d}'] is not None
                        else f'{"-":>16}' for _, d, _ in PROBES)
        print(f'{r["run"]:>4} {r["setting"][:24]:24} {r["K"]!s:>4} {r["minutes"]:5.0f} {sp:>12} {effs}')
    json.dump(rows, open(a.out_stem + '_per_run.json', 'w'), indent=1)

    fig, axs = plt.subplots(2, 3, figsize=(15, 7.5), constrained_layout=True)
    for col, (_, d, name) in enumerate(PROBES):
        for row_i, xkey in enumerate(('run', 'K')):
            ax = axs[row_i, col]
            style(ax)
            for K in sorted({r['K'] for r in rows if r['K'] is not None}):
                rr = [r for r in rows if r['K'] == K and r[f'eff{d}'] is not None]
                if not rr:
                    continue
                x = np.array([r[xkey] for r in rr], float)
                if xkey == 'K':
                    x = x + np.linspace(-0.6, 0.6, len(rr))
                ax.errorbar(x, [r[f'eff{d}'] for r in rr], yerr=[r[f'err{d}'] for r in rr], fmt='o', ms=4, lw=1,
                            color=SET_COL.get(K, MUTED), label=f'safe {K:+d} DAC' if K else 'safe')
            ax.set_ylim(-0.05, 1.0)
            ax.set_xlabel('run' if xkey == 'run' else 'threshold: safe + K [DAC] (runs spread sideways)', fontsize=8)
            ax.set_ylabel(f'efficiency (probe r < {a.probe_r} mm)', fontsize=8)
            ax.set_title(name, fontsize=9, color=INK)
            if col == 0 and row_i == 0:
                ax.legend(fontsize=7, frameon=False, ncol=2)
    fig.suptitle(f'{os.path.basename(os.path.normpath(a.dir))}: efficiency per run (one point per run)', fontsize=11)
    fig.savefig(a.out_stem + '_per_run.png', dpi=120)


if __name__ == '__main__':
    sys.exit(main())
