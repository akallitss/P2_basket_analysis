"""Shape likelihood on a grid of bottom-paddle displacements (along / across the band), top paddle on the band."""
import sys, glob
import numpy as np
sys.path.insert(0, sys.argv[1])
import trigger_acceptance_mc as T
import pandas as pd
from scipy.spatial import cKDTree

J = '/local/home/ak271430/Documents/PostDocSaclay/data/LAB_Measurements/VMM-cosmics/vmm_threshold_scan/analysis_wk2/eff_definitions/json/'
runs = T.load_data(sorted(glob.glob(J + 'scintor_run*.json')))
pads = {nm: runs[0][nm + ':pads'] for nm in T.NAMES}
allp = pd.read_csv(T.VS.MAP_CSV).set_index('channel_id').sort_index()
tree = cKDTree(allp[['pad_cx', 'pad_cy']].to_numpy())
centre, ang = T.band_axis(np.vstack([pads[nm] for nm in T.NAMES]))
a = np.radians(ang)
u = np.array([np.cos(a), np.sin(a)])
v = np.array([-np.sin(a), np.cos(a)])
zs = {'det3 P2_IN': 320., 'det2 P2_MID': 630., 'det1 P2_OUT': 940.}
zf = (0., -100., 1225., 1275.)
ch, ce = {}, {}
for nm in T.NAMES:
    d, i = tree.query(pads[nm])
    read = allp.iloc[np.sort(i)]
    ch[nm] = T.Chamber(read, np.ones(len(read), bool))
    ce[nm] = read[['pad_cx', 'pad_cy']].to_numpy()
box = (centre[0] - 750, centre[0] + 750, centre[1] - 700, centre[1] + 700)
P = T.Pool(np.random.default_rng(5), float(sys.argv[2]), box, zs, zf, 2.0, ch, ce, 3.5)
print(len(P.x0), 'pool')
data = {nm: T.pooled(runs, nm) for nm in T.NAMES}
lo = np.floor((centre - 520) / T.BIN) * T.BIN
edges = (np.arange(lo[0], lo[0] + 1041, T.BIN), np.arange(lo[1] - 100, lo[1] + 1041, T.BIN))
F = T.ShapeFit(P, data, edges)
base = np.array([centre[0], centre[1], ang, centre[0], centre[1], ang])
ref = F.nll(base)
print('along \\ across', [c for c in range(-60, 61, 30)])
for al in range(0, -401, -50):
    row = []
    for ac in range(-60, 61, 30):
        p = base.copy()
        p[3:5] += al * u + ac * v
        row.append(F.nll(p) - ref)
    print(f'{al:5d}', ' '.join(f'{r:8.0f}' for r in row))
