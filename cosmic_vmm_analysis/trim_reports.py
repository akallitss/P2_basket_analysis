#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""What the per-gain threshold calibration (bscan_target.py reports) actually achieved, per gain and noise rate R
(20 Nov 2026). For each report: channels, channels WITHOUT a measured noise-rate crossing (given the chip median
target: 'chip_median_used'), clipped channels (trim would have to exceed 31), the measured targets in charge
(target mV / nominal gain), the per-chip spread of the measured targets against the trim range (31 x ~1.1 mV), and
how much the per-chip global threshold is set by its single highest channel.

For the channels measured at several gains: the SAME channel's noise-rate target in charge at gain 3 / 6 / 9 (does the
noise in fC grow with gain, or is it only the trim range?).

Usage: trim_reports.py REPORT_DIR   (hvon_g{3.0,6.0,9.0}_R{2,10,50}_report.json, hvon{2,10,50}_v1_report.json = 4.5)
"""
import json
import os
import statistics as st
import sys

D = sys.argv[1]
TRIM_RANGE_MV = 31 * 1.1


def files():
    for R in (2, 10, 50):
        yield 4.5, R, f'hvon{R}_v1_report.json'
    for g in ('3.0', '6.0', '9.0'):
        for R in (2, 10, 50):
            yield float(g), R, f'hvon_g{g}_R{R}_report.json'


per_ch = {}
print('gain  R | channels | no crossing -> chip median | clipped | target fC: p10 / median / p90 (measured channels) | '
      'per-chip spread max-min (mV), median | chips spread > trim range | top channel above 2nd (mV), median')
for g, R, f in files():
    p = os.path.join(D, f)
    if not os.path.exists(p):
        continue
    r = json.load(open(p))
    ntot = nmed = ncl = 0
    spreads, tg, gaps = [], [], []
    for chip, v in r.items():
        med = set(int(x) for x in (v.get('chip_median_used') or []))
        T = {int(k): x for k, x in (v.get('targets') or {}).items()}
        meas = {k: x for k, x in T.items() if k not in med}
        ntot += len(T)
        nmed += len(med)
        ncl += len(v['clipped']) if isinstance(v['clipped'], list) else int(v['clipped'])
        if len(meas) >= 2:
            s = sorted(meas.values())
            spreads.append(s[-1] - s[0])
            gaps.append(s[-1] - s[-2])
        tg += [x / g for x in meas.values()]
        for k, x in meas.items():
            per_ch.setdefault((int(chip), k, R), {})[g] = x / g
    if not tg:
        continue
    q = st.quantiles(tg, n=10)
    print(f'{g:4} {R:2} | {ntot:5} | {nmed:4} ({100 * nmed / max(ntot, 1):3.0f} %) | {ncl:4} | {q[0]:.2f} / '
          f'{st.median(tg):.2f} / {q[-1]:.2f} | {st.median(spreads):5.1f} | {sum(s > TRIM_RANGE_MV for s in spreads)}/'
          f'{len(spreads)} | {st.median(gaps):.1f}')

print('\nsame channel, same R: noise-rate target in fC at gain B divided by the one at gain A (median over channels)')
for R in (2, 10, 50):
    for a, b in ((3.0, 6.0), (3.0, 9.0), (6.0, 9.0), (4.5, 9.0)):
        rat = [v[b] / v[a] for (c, k, r_), v in per_ch.items() if r_ == R and a in v and b in v and v[a] > 0]
        if rat:
            qq = st.quantiles(rat, n=4)
            print(f'  R={R:2}: gain {b:g} / gain {a:g}: median {st.median(rat):.2f} [{qq[0]:.2f}, {qq[2]:.2f}] '
                  f'({len(rat)} channels)')
