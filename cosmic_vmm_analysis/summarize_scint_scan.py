#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Table and figure of a scintillator-referenced threshold scan (run_scint_scan.sh output), one point per run.

Per run: the threshold offset K from the sub-run label (..._Kp64_..., _K0_, _Km8_ = safe + 64, 0, -8 DAC), the
fraction of each chamber's pads that are live (masks grow at low thresholds), and the efficiency of each chamber on
all tags and on the live area only (scint_efficiency.py), with scintillator-confirmed tags (coincidence + one clean cluster in each of the other two chambers): probe
hit anywhere and within 50 mm, for the hardware triple (scint3) and the offline (top OR) AND bottom trigger (scintor);
plus the chamber-only tag-and-probe (efficiency_tp.py, a lower bound). "veto" = fraction of the run removed by the
HV-monitor spark veto (spark_veto.py, when run_scint_scan.sh applied it).

Usage: summarize_scint_scan.py DIR OUT_STEM
"""
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

DETS = ((3, 'det3 P2_IN'), (2, 'det2 P2_MID'), (1, 'det1 P2_OUT'))


def k_of(title):
    m = re.search(r'_K([pm]?)(\d+)', title)
    return None if not m else (-1 if m[1] == 'm' else 1) * int(m[2])


def main():
    d, stem = sys.argv[1], sys.argv[2]
    rows = []
    for f in sorted(glob.glob(os.path.join(d, 'eff_run*.json')), key=lambda p: int(re.findall(r'\d+', os.path.basename(p))[0])):
        n = int(re.findall(r'\d+', os.path.basename(f))[0])
        E = json.load(open(f))
        row = {'run': n, 'title': E['title'], 'K': k_of(E['title']), 'hours': E['span_s'] / 3600}
        for tag in ('scint3', 'scintor'):
            p = os.path.join(d, f'{tag}_run{n}.json')
            if not os.path.exists(p):
                continue
            S = json.load(open(p))
            row[f'{tag}_trig'] = S['n_trig']
            if S.get('bad_times'):
                row['vetoed_frac'] = S['bad_times']['vetoed_frac']
            for dd, nm in DETS:
                P = S[f'probe_{nm}']
                row[f'{tag}_n{dd}'] = P['n_tagged']
                row[f'{tag}_any{dd}'] = (P['any_position']['eff'], P['any_position']['err'])
                row[f'{tag}_50_{dd}'] = (P['within_50mm']['eff'], P['within_50mm']['err'])
                row[f'live{dd}'] = S['chambers'][nm].get('live_pad_frac')
                if 'live_area' in P:
                    L = P['live_area']
                    row[f'{tag}_la_any{dd}'] = (L['any_position']['eff'], L['any_position']['err'])
                    row[f'{tag}_la_50_{dd}'] = (L['within_50mm']['eff'], L['within_50mm']['err'])
                    row[f'{tag}_la_frac{dd}'] = P['frac_tags_on_live_area']
        for dd, nm in DETS:
            r = E[f'probe_det{dd}']['vs_probe_r']['50']
            row[f'tp50_{dd}'] = (r.get('eff'), r.get('err'))
        rows.append(row)
    print(f'{"run":>4} {"K":>4} {"h":>4} {"veto":>5} {"trig3":>6} {"trigOR":>6}  ' +
          '  '.join(f'{nm:^44s}' for _, nm in DETS))
    print(' ' * 36 + '  '.join(f'{"live":>5} {"OR any":>12} {"OR 50mm":>12} {"chambers":>9}' for _ in DETS))
    for r in rows:
        cells = []
        for dd, _ in DETS:
            a = r.get(f'scintor_any{dd}', (np.nan, np.nan))
            b = r.get(f'scintor_50_{dd}', (np.nan, np.nan))
            cells.append(f'{r.get(f"live{dd}", np.nan):5.2f} {a[0]:6.3f}+-{a[1]:.3f} {b[0]:6.3f}+-{b[1]:.3f} {r[f"tp50_{dd}"][0] or np.nan:9.3f}')
        print(f'{r["run"]:>4} {r["K"]!s:>4} {r["hours"]:4.1f} {100 * r.get("vetoed_frac", 0):4.1f}% {r.get("scint3_trig", 0):6d} {r.get("scintor_trig", 0):6d}  ' + '  '.join(cells))
    print('\nlive area only: tags whose predicted point has every nearby probe pad live (masked / silent pads excluded)')
    print(f'{"run":>4} {"K":>4}  ' + '  '.join(f'{nm:^37s}' for _, nm in DETS))
    print(' ' * 11 + '  '.join(f'{"pads":>5} {"tags":>5} {"OR any":>12} {"OR 50mm":>12}' for _ in DETS))
    for r in rows:
        cells = []
        for dd, _ in DETS:
            a = r.get(f'scintor_la_any{dd}', (np.nan, np.nan))
            b = r.get(f'scintor_la_50_{dd}', (np.nan, np.nan))
            cells.append(f'{r.get(f"live{dd}", np.nan):5.2f} {r.get(f"scintor_la_frac{dd}", np.nan):5.2f} '
                         f'{a[0]:6.3f}+-{a[1]:.3f} {b[0]:6.3f}+-{b[1]:.3f}')
        print(f'{r["run"]:>4} {r["K"]!s:>4}  ' + '  '.join(cells))
    json.dump(rows, open(stem + '_scan.json', 'w'), indent=1, default=float)

    fig, axs = plt.subplots(1, 4, figsize=(17, 4.3), constrained_layout=True)
    for j, (dd, nm) in enumerate(DETS):
        ax = axs[j]
        style(ax)
        for tag, key, mk, lab in (('scintor', 'any', 'o', '(top OR) AND bottom, anywhere'), ('scintor', '50_', 's', '(top OR) AND bottom, 50 mm'),
                                  ('scintor', 'la_any', 'D', '(top OR) AND bottom, anywhere, live area only'),
                                  ('scint3', 'any', '^', 'triple, anywhere')):
            pts = [(r['K'], r[f'{tag}_{key}{dd}']) for r in rows if r['K'] is not None and f'{tag}_{key}{dd}' in r]
            if pts:
                x = np.array([p[0] for p in pts], float)
                ax.errorbar(x + {'o': 0, 's': 1.2, '^': -1.2, 'D': 0.6}[mk], [p[1][0] for p in pts], yerr=[p[1][1] for p in pts], fmt=mk,
                            color=MUTED if mk == '^' else COL[dd], mfc='white' if mk == 's' else None,
                            mec=INK if mk == 'D' else None, ms=6, lw=1.2, label=lab)
        pts = [(r['K'], r[f'tp50_{dd}'][0]) for r in rows if r['K'] is not None and r[f'tp50_{dd}'][0] is not None]
        if pts:
            ax.plot([p[0] for p in pts], [p[1] for p in pts], 'x', color=INK, ms=6, label='chambers alone (50 mm, lower bound)')
        ax.set_ylim(0, 1)
        ax.invert_xaxis()
        ax.set_xlabel('threshold: safe + K [DAC]', fontsize=8)
        ax.set_ylabel('efficiency', fontsize=8)
        ax.set_title(nm, fontsize=9, color=INK)
        if j == 0:
            ax.legend(fontsize=7, frameon=False, loc='lower left')
    ax = axs[3]
    style(ax)
    for dd, nm in DETS:
        pts = [(r['K'], r[f'live{dd}']) for r in rows if r['K'] is not None and r.get(f'live{dd}') is not None]
        if pts:
            ax.plot([p[0] for p in pts], [p[1] for p in pts], 'o-', color=COL[dd], label=nm)
    ax.set_ylim(0, 1.02)
    ax.invert_xaxis()
    ax.set_xlabel('threshold: safe + K [DAC]', fontsize=8)
    ax.set_ylabel('live pads / mapped pads of the chamber', fontsize=8)
    ax.set_title('live pads', fontsize=9, color=INK)
    ax.legend(fontsize=7, frameon=False)
    fig.suptitle(f'{os.path.basename(os.path.normpath(d))}: efficiency against the threshold, scintillator reference, one point per run', fontsize=11)
    fig.savefig(stem + '_scan.png', dpi=120)


if __name__ == '__main__':
    sys.exit(main())
