"""p2_study_common.py -- shared pieces of the 34/35/36 studies (on the
p2_study_cache products): transform, active area, baseline clustering and the
stage-06 efficiency definition, so a study's numbers are comparable to 06."""
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

import p2_align as pa


def make_T(g, R=None, s=None, pad_mean=None, m3_mean=None):
    return pa.Transform(g['R'] if R is None else R, float(g['s'] if s is None else s),
                        g['pad_mean'] if pad_mean is None else pad_mean,
                        g['m3_mean'] if m3_mean is None else m3_mean)


def apply_T(T, xy):
    x, y = T.apply(xy[..., 0].ravel(), xy[..., 1].ravel())
    return np.stack([x, y], -1).reshape(xy.shape)


def ray_dirs(rays):
    """Per-ray slopes tx = dx/dz, ty = dy/dz and zenith angle [deg]."""
    dz = rays['Z_Up'].to_numpy() - rays['Z_Down'].to_numpy()
    tx = (rays['X_Up'].to_numpy() - rays['X_Down'].to_numpy()) / dz
    ty = (rays['Y_Up'].to_numpy() - rays['Y_Down'].to_numpy()) / dz
    return tx, ty, np.degrees(np.arctan(np.hypot(tx, ty)))


def active_mask(xy, g, T, active_r=None, guard=None):
    """Stage-06 denominator: within ACTIVE_R of a live pad centre, outside the
    dead-pad Voronoi cells and the dead-tile guard band."""
    active_r = float(g['active_r'] if active_r is None else active_r)
    guard = float(g['dead_guard'] if guard is None else guard)
    live = apply_T(T, g['pad_xy'])
    nn, idx = cKDTree(live).query(xy)
    ok = nn <= active_r
    if len(g['dead_xy']):
        dd, _ = cKDTree(apply_T(T, g['dead_xy'])).query(xy)
        ok &= dd >= nn
        if guard > 0:
            ok &= pa.dist_to_quads(xy, apply_T(T, g['dead_verts'])) >= guard
    return ok, nn, idx


def cc_clusters(hits, link=20.0, time_cut=None):
    """All connected-component clusters per event (pad frame).
    time_cut [ns]: additionally require |t_i - t_j| < time_cut for a link.
    Returns DataFrame [eventId, x_pad, y_pad, q, n_pad, t, qmax]."""
    out = []
    for e, gdf in hits.groupby('eventId', sort=False):
        x = gdf['x_pad'].to_numpy(); y = gdf['y_pad'].to_numpy(); w = gdf['amp'].to_numpy()
        t = gdf['time'].to_numpy()
        n = len(x)
        if n > 1:
            pr = cKDTree(np.column_stack([x, y])).query_pairs(link, output_type='ndarray')
            if time_cut is not None and len(pr):
                pr = pr[np.abs(t[pr[:, 0]] - t[pr[:, 1]]) < time_cut]
            G = coo_matrix((np.ones(len(pr)), (pr[:, 0], pr[:, 1])), shape=(n, n)) if len(pr) \
                else coo_matrix((n, n))
            nc, lab = connected_components(G, directed=False)
        else:
            nc, lab = 1, np.zeros(1, int)
        q = np.bincount(lab, weights=w, minlength=nc)
        qs = np.where(q > 0, q, 1)
        cx = np.bincount(lab, weights=w * x, minlength=nc) / qs
        cy = np.bincount(lab, weights=w * y, minlength=nc) / qs
        ct = np.bincount(lab, weights=w * t, minlength=nc) / qs
        qm = np.zeros(nc); np.maximum.at(qm, lab, w)
        nn = np.bincount(lab, minlength=nc)
        for k in range(nc):
            out.append((e, cx[k], cy[k], q[k], nn[k], ct[k], qm[k]))
    return pd.DataFrame(out, columns=['eventId', 'x_pad', 'y_pad', 'q', 'n_pad', 't', 'qmax'])


def nearest(ev, xy, clu_xy_m3, clu_ev):
    """Distance from each ray to its event's closest cluster (inf if none) and
    the row index of that cluster in the cluster table (-1 if none)."""
    t = pd.DataFrame({'eventId': ev, 'i': np.arange(len(ev)), 'tx': xy[:, 0], 'ty': xy[:, 1]})
    c = pd.DataFrame({'eventId': clu_ev, 'j': np.arange(len(clu_ev)),
                      'cx': clu_xy_m3[:, 0], 'cy': clu_xy_m3[:, 1]})
    m = t.merge(c, on='eventId')
    d = np.full(len(ev), np.inf); j = np.full(len(ev), -1)
    if len(m):
        m['d'] = np.hypot(m.cx - m.tx, m.cy - m.ty)
        b = m.loc[m.groupby('i')['d'].idxmin()]
        d[b.i.to_numpy()] = b.d.to_numpy(); j[b.i.to_numpy()] = b.j.to_numpy()
    return d, j


def binom_err(k, n):
    p = k / n
    return 100 * np.sqrt(p * (1 - p) / n)
