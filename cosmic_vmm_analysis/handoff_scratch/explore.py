"""Which change moves the MC tag centroids like the data (along-band medians det3 +8, det2 -29, det1 -88 mm)?
(a) bottom paddle displaced along the band; (b) a lateral shear of the chamber stack (det1 frame displaced)."""
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

for nm in T.NAMES:
    x, y = T.pooled(runs, nm)
    s = (np.c_[x, y] - centre)
    print(f'data {nm}: along median {np.median(s @ u):6.1f} mean {np.mean(s @ u):6.1f} | across median {np.median(s @ v):6.1f}')

zs = {'det3 P2_IN': 320., 'det2 P2_MID': 630., 'det1 P2_OUT': 940.}
zf = (0., -100., 1225., 1275.)


def build(offsets):
    ch, ce = {}, {}
    for nm in T.NAMES:
        d, i = tree.query(pads[nm])
        read = allp.iloc[np.sort(i)]
        ch[nm] = T.Chamber(read, np.ones(len(read), bool), offsets.get(nm, (0., 0.)))
        ce[nm] = read[['pad_cx', 'pad_cy']].to_numpy()
    box = (centre[0] - 750, centre[0] + 750, centre[1] - 700, centre[1] + 700)
    rng = np.random.default_rng(3)
    P = T.Pool(rng, 4e6, box, zs, zf, 2.0, ch, ce, 3.5)
    # the predicted points are computed from TRUE stand-frame positions; data are in each chamber's own frame:
    # convert: measured in own frame = stand - offset
    for i, nm in enumerate(T.NAMES):
        j, k = [t for t in range(3) if t != i]
        A, B = T.NAMES[j], T.NAMES[k]
        ok = P.hit[j] & P.hit[k]
        pa = P.true[A] - np.array(offsets.get(A, (0, 0)))
        pb = P.true[B] - np.array(offsets.get(B, (0, 0)))
        f = (zs[nm] - zs[A]) / (zs[B] - zs[A])
        P.pred[nm] = (ok, pa + f * (pb - pa))
    return P


def report(P, p, lab):
    m = P.trigger(p)
    out = []
    for nm in T.NAMES:
        ok, xy = P.pred[nm]
        s = (xy[m & ok] - centre)
        out.append(f'{np.median(s @ u):6.1f}')
    print(f'{lab:55s} along medians det3/det2/det1: ' + ' '.join(out))


base = np.array([centre[0], centre[1], ang, centre[0], centre[1], ang])
P0 = build({})
report(P0, base, 'aligned paddles, aligned chambers')
for d in (-100, -200, -300):
    p = base.copy()
    p[3:5] += d * u
    report(P0, p, f'bottom paddle {d} mm along band')
    p = base.copy()
    p[0:2] -= d * u
    report(P0, p, f'top paddle {-d} mm along band')
for d in (40, 80, 120):
    P = build({'det2 P2_MID': tuple(0.5 * d * u), 'det1 P2_OUT': tuple(d * u)})
    report(P, base, f'stack shear: det1 frame +{d} mm along band (det2 half)')
    P = build({'det1 P2_OUT': tuple(d * u)})
    report(P, base, f'det1 alone shifted +{d} mm along band')
