"""Per read-out connector, the beam-style efficiency on the cosmic stand, from the fixed-clock tag dumps:
muon = scintillator-confirmed tag (on - r x off), predicted pad live and within 9 mm (the dump's fiducial);
HIT = a live probe pad within 15 mm of the prediction within +-150 ns of tag + latency (the beam window: -15..+285 ns
around a 135 ns latency); ANY = any live probe hit in the chamber within +-100 ns (n_any > 0).
Connector of the prediction from its chip (cosmic stand, conn-shift -1): IN 2,3->c3 4,5->c4 6,7->c5; MID 8-13; OUT 14-19.
Usage (from $W, env.sh): conn_eff.py   (reads tags/ and pseudo/chunks.json; prints the table per configuration).
"""
import glob, json, re
import numpy as np
chunks = {int(k): v for k, v in json.load(open('pseudo/chunks.json')).items()}
cfg_of = lambda c: re.sub(r'^cfg_gain[0-9.]+_peaktime200_opt_cosmic_', '', c)
NAME = {3: 'P2_IN', 2: 'P2_MID', 1: 'P2_OUT'}
FIRST = {3: 2, 2: 8, 1: 14}
acc = {}
for f in sorted(glob.glob('tags/tags_run*.npz')):
    pid = int(re.search(r'run(\d+)\.npz', f)[1])
    info = chunks.get(pid)
    if info is None or info['run'] == 270 or info['run'] < 268:
        continue
    gain = re.search(r'gain([0-9.]+)', info['cfg'])[1]
    key_cfg = f"g{gain} {cfg_of(info['cfg'])}"
    Z = np.load(f)
    for p in (3, 2, 1):
        dts = Z[f'p{p}_dt_scint']; n = len(dts)
        sel = np.nonzero((dts >= 140) & (dts <= 310))[0]
        w = np.where(Z[f'p{p}_on'], 1.0, -float(Z[f'p{p}_r_scale']))
        conn = 3 + (Z[f'p{p}_vmm'] - FIRST[p]) // 2
        fi, fdt, fd = Z[f'p{p}_fired_nr_i'], Z[f'p{p}_fired_nr_dt'], Z[f'p{p}_fired_nr_d']
        hit = np.zeros(n, bool)
        m = (fd < 15) & (np.abs(fdt) <= 150)
        hit[fi[m]] = True
        anyh = Z[f'p{p}_n_any'] > 0
        for c in (3, 4, 5):
            q = sel[conn[sel] == c]
            for k in (key_cfg, 'ALL'):
                a = acc.setdefault((k, NAME[p], c), np.zeros(3))
                a += [w[q].sum(), (w[q] * hit[q]).sum(), (w[q] * anyh[q]).sum()]
order = ['ALL'] + sorted({k[0] for k in acc if k[0] != 'ALL'})
print('per connector: muons / HIT (pad within 15 mm, +-150 ns) / ANY (any chamber hit) -- on - r off weights')
for k in order:
    cells = []
    for det in ('P2_IN', 'P2_MID', 'P2_OUT'):
        for c in (3, 4, 5):
            Wm, H, A = acc.get((k, det, c), np.zeros(3))
            cells.append(f'{det[3:]} c{c} {H / Wm:.3f}/{A / Wm:.3f} ({Wm:.0f})' if Wm > 20 else f'{det[3:]} c{c} -')
    print(f'{k:28s} ' + ' | '.join(cells))
