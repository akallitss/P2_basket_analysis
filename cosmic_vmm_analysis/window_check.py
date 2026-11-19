"""Are the 'silent' muons late/early hits outside the +-100 ns window? From the fixed-clock tag dumps: for every
scintillator-confirmed tag (on - r off), the nearest-in-time live probe hit within 15 mm of the prediction out to
+-2 us (fired_nr_*), and the efficiency (hit within 15 mm) as a function of the half-window.
Usage (from $W, env.sh): window_check.py   (reads tags/ and pseudo/chunks.json; no side-band subtraction yet).
"""
import glob, json, re, sys
import numpy as np
chunks = {int(k): v for k, v in json.load(open('pseudo/chunks.json')).items()}
W = [50, 100, 150, 200, 300, 500, 1000, 2000]
acc = {}
for f in sorted(glob.glob('tags/tags_run*.npz')):
    pid = int(re.search(r'run(\d+)\.npz', f)[1])
    run = chunks.get(pid, {}).get('run')
    if run is None or run == 270 or run < 268:
        continue
    Z = np.load(f)
    for p in (3, 2, 1):
        if f'p{p}_fired_nr_i' not in Z.files: continue
        dts = Z[f'p{p}_dt_scint']; n = len(dts)
        sel = (dts >= 140) & (dts <= 310)
        on = Z[f'p{p}_on']; r = float(Z[f'p{p}_r_scale'])
        w = np.where(on, 1.0, -r)
        fi, fdt, fd = Z[f'p{p}_fired_nr_i'], Z[f'p{p}_fired_nr_dt'], Z[f'p{p}_fired_nr_d']
        best = np.full(n, np.inf)               # smallest |dt| of a hit within 15 mm
        m = fd < 15
        np.minimum.at(best, fi[m], np.abs(fdt[m]))
        a = acc.setdefault(p, dict(W=0.0, k=np.zeros(len(W)), dt=[]))
        a['W'] += w[sel].sum()
        for j, x in enumerate(W):
            a['k'][j] += (w[sel] * (best[sel] <= x)).sum()
        q = sel & on & np.isfinite(best)
        a['dt'].append(best[q])
NAME = {3: 'P2_IN', 2: 'P2_MID', 1: 'P2_OUT'}
print('efficiency (live hit within 15 mm of the prediction) vs half-window around tag + latency, muons = on - r off')
print('half-window ns:', W)
for p, a in acc.items():
    e = a['k'] / a['W']
    print(f'{NAME[p]}: muons {a["W"]:.0f}  eff', np.round(e, 4).tolist(), f' gain 100 ns -> 2 us: {e[-1] - e[1]:+.4f}')
    d = np.concatenate(a['dt'])
    print(f'    ON tags with a hit within 15 mm in +-2 us: |dt| quantiles 50/90/95/99 %:',
          np.round(np.percentile(d, [50, 90, 95, 99]), 0).tolist(), ' fraction beyond 100 ns', round(float((d > 100).mean()), 4))
