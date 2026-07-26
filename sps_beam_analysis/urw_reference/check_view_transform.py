#!/usr/bin/env python3
"""
Per-view connector-order test that needs NO external anchor.

Both strip maps are uniform: global index k = (map connector)*64 + channel runs
monotonically with the Gerber coordinate over the full 0..127 of every view
(verified on inter_map.txt and strip_map.txt).  So interchanging a view's two
Dream connectors is exactly a cyclic shift of the position array by 64 strips.

The beam spot straddles the middle of both planes.  A cyclic half-plane shift
therefore SPLITS it in two and throws the pieces to opposite edges, which blows
up the width of the profile.  The correct assignment is the one that gives one
compact blob.  Metric: std of the leading-cluster position, and the fraction of
entries in the largest contiguous run of populated bins.

Also does the front<->back position correlation for all four combinations of
(front swapped, back swapped) per axis, which is independent again: it uses two
different map files and no P2 information.

SUPERSEDED 2026-07-26 -- kept for the record, do not use it to set the mapping.
The premise in the second paragraph is false for the BACK: strip_map.txt has
three pitch zones per view (1.0 / 1.5 / 0.5 mm), so the back's position array is
NOT uniform and interchanging its connectors is not a cyclic shift.  This test
also cannot see a reversal of channel order inside a connector at all, and that
is what the back actually has.  explore4_back_map.py scores all four candidate
wirings per view against the front (which the P2 stations independently confirm
points to <1 mm) and finds the back is 'AB_rev' on both views, where this test
said "interchange the connectors".  See urw_lib.VIEW_MODE_DEFAULT.
"""
import os, sys, glob, argparse
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import urw_lib as U

RUN_DIR = '/local/home/banco/P2_data/TB_July2026_H4/runs/drift_mesh_scan_1'
# the map's own order, so this script's shift emulation is the only transform
RAW_MODE = {'x': 'AB', 'y': 'AB'}
NO_FLIP = {'x': False, 'y': False}


def profile_stats(p, bin_mm=2.0, frac=0.10):
    p = p[np.isfinite(p)]
    lo, hi = np.nanmin(p), np.nanmax(p)
    nb = max(4, int(np.ceil((hi - lo) / bin_mm)))
    h, edges = np.histogram(p, bins=nb, range=(lo, hi))
    on = h > frac * h.max()
    # largest contiguous run of populated bins
    best = cur = 0
    for v in on:
        cur = cur + 1 if v else 0
        best = max(best, cur)
    m = np.zeros_like(on)
    # entries inside the best run
    run_end = 0; cur = 0
    for i, v in enumerate(on):
        cur = cur + 1 if v else 0
        if cur == best: run_end = i
    m[run_end - best + 1:run_end + 1] = True
    return dict(n=len(p), mean=p.mean(), std=p.std(),
                frac_in_largest_blob=h[m].sum() / h.sum())


def robust_line(x, y, nsig=2.5, iters=5):
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    keep = np.ones(len(x), bool)
    a = b = np.nan
    for _ in range(iters):
        if keep.sum() < 20: break
        a, b = np.polyfit(x[keep], y[keep], 1)
        r = y - (a * x + b)
        s = np.std(r[keep])
        keep = np.abs(r) < nsig * s
    return a, b, np.std(y[keep] - (a * x[keep] + b)) if keep.sum() else np.nan, keep.mean()


def load(det, sub_run_dir, run_json, mode, min_amp, max_hits):
    geo = U.UrwGeometry(det, run_json, sub_run_name=os.path.basename(sub_run_dir),
                        view_mode=mode, axis_flip=NO_FLIP)
    parts = []
    for chunk in U.iter_hits(U.feu_hit_files(sub_run_dir, feu=geo.feu_num),
                             max_hits=max_hits, progress=False, feu=geo.feu_num):
        c = U.cluster_hits(chunk, geo, min_amp=min_amp)
        if len(c): parts.append(U.leading_points(c))
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run-dir', default=RUN_DIR)
    ap.add_argument('--sub-run', default='nominal_00')
    ap.add_argument('--min-amp', type=float, default=100.0)
    ap.add_argument('--max-hits', type=float, default=3e6)
    args = ap.parse_args()

    run_json = os.path.join(args.run_dir, 'run_config.json')
    sub_dir = os.path.join(args.run_dir, args.sub_run)
    pts = {}
    for det in U.URW_DETS:
        pts[det] = load(det, sub_dir, run_json, RAW_MODE, args.min_amp, int(args.max_hits))
        print(f'{det}: {len(pts[det])} events with a leading cluster', flush=True)

    # ---- test 1: profile compactness, each view on its own -------------------
    print('\n=== profile compactness (no swap vs connectors interchanged) ===')
    print('interchanging the two connectors = cyclic shift of the position array '
          'by 64 strips\n')
    rows = []
    for det in U.URW_DETS:
        d = pts[det]
        for v in ('x', 'y'):
            p = d[v].values
            # emulate the swap by shifting the position by +-half the plane
            lo, hi = np.nanmin(p), np.nanmax(p)
            span = hi - lo
            half = span / 2.0
            ps = np.where(p < lo + half, p + half, p - half)
            a = profile_stats(p); b = profile_stats(ps)
            rows.append(dict(det=det.replace('EIC_uRWELL_', ''), view=v,
                             std_asis=a['std'], blob_asis=a['frac_in_largest_blob'],
                             std_swapped=b['std'], blob_swapped=b['frac_in_largest_blob'],
                             verdict='SWAP' if b['std'] < a['std'] else 'keep'))
    t = pd.DataFrame(rows)
    print(t.to_string(index=False, float_format=lambda x: f'{x:7.3f}'))

    # ---- test 2: front-back correlation, all 4 combos per axis ---------------
    print('\n=== front vs back position correlation, all four combinations ===')
    rows = []
    for v in ('x', 'y'):
        for fs in (False, True):
            for bs in (False, True):
                f = pts['EIC_uRWELL_front'][['eventId', v]].rename(columns={v: 'f'})
                b = pts['EIC_uRWELL_back'][['eventId', v]].rename(columns={v: 'b'})
                j = f.merge(b, on='eventId').dropna()
                def shift(p):
                    lo, hi = np.nanmin(p), np.nanmax(p); h = (hi - lo) / 2
                    return np.where(p < lo + h, p + h, p - h)
                fv = shift(j['f'].values) if fs else j['f'].values
                bv = shift(j['b'].values) if bs else j['b'].values
                a, c, rms, keep = robust_line(fv, bv)
                rows.append(dict(axis=v, front='swap' if fs else 'keep',
                                 back='swap' if bs else 'keep', n=len(j),
                                 pearson=np.corrcoef(fv, bv)[0, 1], slope=a,
                                 offset=c, resid_rms=rms, on_ridge=keep))
    t2 = pd.DataFrame(rows)
    print(t2.to_string(index=False, float_format=lambda x: f'{x:8.3f}'))
    print('\nbest per axis (slope closest to +1 with the highest on-ridge fraction):')
    for v in ('x', 'y'):
        s = t2[t2['axis'] == v].copy()
        s['score'] = s['on_ridge'] / (1 + np.abs(s['slope'] - 1))
        w = s.sort_values('score', ascending=False).iloc[0]
        print(f"  {v}: front={w['front']} back={w['back']}  slope={w['slope']:+.3f} "
              f"rms={w['resid_rms']:.2f}mm on_ridge={w['on_ridge']:.2f}")


if __name__ == '__main__':
    main()
