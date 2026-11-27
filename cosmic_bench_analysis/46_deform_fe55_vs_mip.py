#!/usr/bin/env python3
"""46_deform_fe55_vs_mip.py -- det2 pad by pad: Fe55 (July 2026) vs cosmic MIP (November 2026).

Why: an Fe55 conversion is a point charge, so its pulse height measures the
AMPLIFICATION gain only and its pulse shape has no drift-time spread; a muon
ionises the whole drift gap, so its charge also scales with the local DRIFT-gap
thickness and its pulse is stretched by the drift time.  Comparing the two on
the same pads separates amplification-gap from drift-gap effects.

Data
  Fe55   det2_fe55scan1: p2_fe55_det2_det3_mesh_scan_7-18-26, 12 x 5 min,
         mesh 420 -> 365 V (drift = mesh + 280 V), DREAM self-trigger, FEU 3,
         connectors 4-7 (512 pads).  The chamber on FEU 3 / CAEN 8:0-8:1 is
         det2 according to both the Fe55 run_config ("P2_OUT ... (det2)") and
         the SPS run_configs (station P2_IN = det2 on FEU 3 / 8:0-8:1); the
         station NAME differs between the two (see chamber_history.py).
  MIP    det2_long_repos: 25 Nov 2026, mesh 430 / drift 730 V, 7.6 h, FEUs 6/7,
         trimmed stage-34 alignment (34_clustering_matching/alignment.npz,
         rays_aligned.pkl), pads in the PAD frame via the inverse transform.

Hit-level quantities (mm_dream_reconstruction/waveform_analysis):
  time         30 % constant-fraction leading-edge time (linear interpolation
               between samples), ftst-corrected, ns
  time_of_max  parabola vertex through the 3 samples around the maximum
               (saturated: crossing of the two edge-slope lines), ftst-corr., ns
  rise         = time_of_max - time  [ns]: the same definition in both analyzer
               versions (ftst cancels; absolute self-trigger times never used)
  integral     sum of baseline-subtracted samples [ADC x sample].
               JULY files (analyzer before 2026-07-24, commit 1cad061): sum over
               the samples strictly inside the threshold crossings only.
               NOVEMBER files (8ec6769+): full span to the baseline return.
  width        = 60 ns x integral / amplitude  [ns]: an effective width, but
               NOT the same definition in July and November.
  TOT          JULY: integer number of samples above 5 sigma (x 60 ns, quantised)
               on the raw waveform; NOVEMBER: interpolated crossings on the
               boxcar matched-filter gate.  NOT the same definition either.
  => only `rise` is compared ACROSS periods; width and TOT are used within
     each period (slopes vs gain) and their difference map only after removing
     a constant (structure), with that caveat.

Products: <det2_long_repos>/46_deform_fe55_vs_mip/ (pad tables, region tables,
figures, summary.txt).  Fe55 / MIP per-event tables are cached there
(fe55_events.pkl, mip_tracks.pkl); --rebuild re-streams the hit files.

Usage: python3 46_deform_fe55_vs_mip.py [--edge 2.5] [--rebuild]
"""
import argparse
import importlib.util
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                       # noqa: E402
from matplotlib.collections import PolyCollection     # noqa: E402
from scipy import odr, stats                          # noqa: E402
from scipy.spatial import cKDTree                     # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p2_qa_config as qa        # noqa: E402
import p2_mapping as pmap        # noqa: E402
import p2_io as p2io             # noqa: E402
import p2_sparks as ps           # noqa: E402
import p2_study_cache as sc      # noqa: E402
import p2_study_common as cm     # noqa: E402

plt.rcParams.update({'figure.facecolor': 'white', 'axes.facecolor': 'white',
                     'savefig.facecolor': 'white'})

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location('s18', os.path.join(HERE, '18_fe55_spectra.py'))
s18 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s18)

REGIONS_CSV = ('/local/home/ak271430/Documents/PostDocSaclay/data/Cosmic_Bench/Analysis/det2/'
               'p2_det2_long_run_430_repositioned_11-25-26/long_run_det2_430_730/43_deformation/regions.csv')
TPS = 60.0           # ns per sample, both runs (run_config dream_daq_info.sample_period)
CLUS_R = 15.0        # pads within this of the reference pad centre form the cluster [mm]
FE_FLOOR = 200.0     # Fe55 cluster-charge floor (other-detector-trigger noise bulge), as stage 18
WIN = 0.25           # per-pad photopeak window +-25 %
NBOOT = 80
RNG = np.random.default_rng(46)


# ----------------------------------------------------------------------------- helpers
def inverse(T, xy):
    """M3 -> pad frame (inverse of p2_align.Transform.apply)."""
    return ((np.asarray(xy) - T.m3_mean) / T.s) @ T.R + T.pad_mean


def peak_iter(q, m0):
    """Per-pad photopeak: median in [0.6, 1.6] m0, then 4 x median in +-WIN."""
    q = q[(q > 0.6 * m0) & (q < 1.6 * m0)]
    if len(q) < 5:
        return np.nan
    m = np.median(q)
    for _ in range(4):
        w = q[(q > (1 - WIN) * m) & (q < (1 + WIN) * m)]
        if len(w) < 5:
            return np.nan
        m = np.median(w)
    return m


def peak_tmean(q, m0):
    """Alternative estimator: mean inside +-WIN of the iterated median."""
    m = peak_iter(q, m0)
    if not np.isfinite(m):
        return np.nan
    w = q[(q > (1 - WIN) * m) & (q < (1 + WIN) * m)]
    return w.mean()


def boot(fn, x, n=NBOOT, *a):
    x = np.asarray(x)
    if len(x) < 5:
        return np.nan
    v = [fn(x[RNG.integers(0, len(x), len(x))], *a) for _ in range(n)]
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    return float(np.std(v)) if len(v) > n // 2 else np.nan


def tmean_low(x, f=0.7):
    x = np.sort(np.asarray(x))
    k = max(1, int(round(f * len(x))))
    return float(x[:k].mean())


def tmean10(x):
    return float(stats.trim_mean(x, 0.1)) if len(x) else np.nan


def wmean(v, e):
    ok = np.isfinite(v) & np.isfinite(e) & (e > 0)
    if ok.sum() == 0:
        return np.nan, np.nan, np.nan, 0
    w = 1 / e[ok] ** 2
    m = np.sum(w * v[ok]) / w.sum()
    chi2 = np.sum(w * (v[ok] - m) ** 2)
    n = int(ok.sum())
    err = np.sqrt(1 / w.sum()) * (np.sqrt(chi2 / (n - 1)) if n > 1 and chi2 > n - 1 else 1.0)
    return m, err, chi2, n


def odr_fit(x, y, sx, sy):
    """Straight line with errors in both coordinates. Returns slope, err, intercept."""
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(sx) & np.isfinite(sy) & (sx > 0) & (sy > 0)
    if ok.sum() < 4:
        return np.nan, np.nan, np.nan
    b0 = np.polyfit(x[ok], y[ok], 1)
    o = odr.ODR(odr.RealData(x[ok], y[ok], sx=sx[ok], sy=sy[ok]),
                odr.Model(lambda B, X: B[0] * X + B[1]), beta0=b0).run()
    return float(o.beta[0]), float(o.sd_beta[0]), float(o.beta[1])


def boot_corr(x, y, n=2000):
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) < 4:
        return np.nan, np.nan, 0
    r = np.corrcoef(x, y)[0, 1]
    idx = RNG.integers(0, len(x), (n, len(x)))
    rb = np.array([np.corrcoef(x[i], y[i])[0, 1] for i in idx])
    return float(r), float(np.std(rb)), len(x)


def moran(v, xy, link=CLUS_R, nperm=2000):
    """Moran's I of v on the pad-adjacency graph (centres within `link` mm),
    with a permutation p-value (one-sided, positive autocorrelation)."""
    ok = np.isfinite(v)
    v, xy = v[ok], xy[ok]
    if len(v) < 10:
        return np.nan, np.nan
    pairs = cKDTree(xy).query_pairs(link, output_type='ndarray')
    if not len(pairs):
        return np.nan, np.nan

    def I(z):
        z = z - z.mean()
        return len(z) / (2 * len(pairs)) * 2 * np.sum(z[pairs[:, 0]] * z[pairs[:, 1]]) / np.sum(z * z)
    i0 = I(v)
    perm = np.array([I(RNG.permutation(v)) for _ in range(nperm)])
    return float(i0), float((np.sum(perm >= i0) + 1) / (nperm + 1))


def demean_by(v, grp):
    s = pd.Series(v)
    return (s - s.groupby(np.asarray(grp)).transform('median')).to_numpy()


def reliability(a, b):
    """Split-half correlation and Spearman-Brown full-sample reliability."""
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 5:
        return np.nan, np.nan
    r = np.corrcoef(a[ok], b[ok])[0, 1]
    return float(r), float(2 * r / (1 + r)) if r > -1 else np.nan


# ----------------------------------------------------------------------------- Fe55
def fe55_extract(cfgj, ct):
    subruns = s18.find_subruns(cfgj)
    br = ['eventId', 'channel', 'feu', 'amplitude', 'time', 'time_of_max', 'integral',
          'time_over_threshold', 'saturated', 'trigger_timestamp_ns']
    out, stat = [], []
    for sub, hv in subruns:
        sd = cfgj.subrun_dir(sub)
        hd = os.path.join(sd, 'combined_hits_root')
        t_min = s18.settle_t_min(cfgj, sd, p2io.hit_files(hd)[0])
        hv_csv = os.path.join(sd, 'hv_monitor.csv')
        veto = ps.SparkVeto.from_csv(hv_csv, cfgj) if os.path.isfile(hv_csv) else None
        parts = []
        for df in p2io.iter_hits(hd, br, ct.attrs['feus'], progress=False, min_amp=cfgj.MIN_AMP):
            if t_min > 0:
                df = df[df['trigger_timestamp_ns'].astype(np.int64) / 1e9 >= t_min]
            if veto is not None and len(df):
                df, _ = veto.apply(df)
            h = pmap.attach_pads_to_hits(df, ct)
            h = h[h['mapped'] & h['pad_cx'].notna()]
            parts.append(pd.DataFrame({
                'eventId': h.eventId.to_numpy(np.int64), 'channel_id': h.channel_id.to_numpy(np.int32),
                'x': h.pad_cx.to_numpy(np.float32), 'y': h.pad_cy.to_numpy(np.float32),
                'a': h.amplitude.to_numpy(np.float32), 't': h.time.to_numpy(np.float32),
                'tm': h.time_of_max.to_numpy(np.float32), 'itg': h.integral.to_numpy(np.float32),
                'tot': h.time_over_threshold.to_numpy(np.float32), 'sat': h.saturated.to_numpy(bool)}))
            del df, h
        H = pd.concat(parts, ignore_index=True)
        del parts
        H = H.sort_values('a', ascending=False).drop_duplicates(['eventId', 'channel_id'])
        L = H.drop_duplicates('eventId').set_index('eventId')
        lx = L.x.reindex(H.eventId).to_numpy(); ly = L.y.reindex(H.eventId).to_numpy()
        d = np.hypot(H.x.to_numpy() - lx, H.y.to_numpy() - ly)
        H['near'] = d <= CLUS_R
        H['nb'] = (d > 0.5) & (d <= CLUS_R)
        g = H[H.near].groupby('eventId')
        gn = H[H.nb].groupby('eventId')
        ev = pd.DataFrame({'q': g.a.sum(), 'sat_any': g.sat.any()})
        ev['n_nb'] = gn.size().reindex(ev.index).fillna(0).astype(np.int16)
        ev['a_nb'] = gn.a.max().reindex(ev.index).fillna(0).astype(np.float32)
        ev = ev.join(L[['channel_id', 'a', 't', 'tm', 'itg', 'tot', 'sat']])
        ev['hv'] = hv
        n_all = H.groupby('eventId').size()
        ev['n_hits'] = n_all.reindex(ev.index).astype(np.int32)
        out.append(ev.reset_index())
        stat.append(dict(hv=hv, subrun=sub, n_events=len(ev), t_min_cut=t_min,
                         n_vetoed_s=(veto.vetoed_seconds if veto is not None else 0.0)))
        print(f'  Fe55 {hv} V: {len(ev):,} events', flush=True)
        del H, L
    return pd.concat(out, ignore_index=True), pd.DataFrame(stat)


def fe55_pad_peaks(E, hvs, mask=None):
    """Per (pad, HV) photopeak of q (cluster) and a (lead pad), with bootstrap errors."""
    rows = []
    for hv in hvs:
        e = E[(E.hv == hv)] if mask is None else E[(E.hv == hv) & mask]
        fit = s18.fit_photopeak(e.q.to_numpy(), FE_FLOOR)
        m0 = fit['mu'] if fit['ok'] else np.nan
        fa = s18.fit_photopeak(e.a.to_numpy(), FE_FLOOR)
        m0a = fa['mu'] if fa['ok'] else np.nan
        if not np.isfinite(m0) or fit['near_floor']:
            continue
        for pid, grp in e[e.q > FE_FLOOR].groupby('channel_id'):
            q = grp.q.to_numpy(); a = grp.a.to_numpy()
            mq = peak_iter(q, m0)
            ma = peak_iter(a, m0a) if np.isfinite(m0a) else np.nan
            row = dict(channel_id=pid, hv=hv, n=len(q), m0=m0, peak_q=mq,
                       peak_q_tm=peak_tmean(q, m0), peak_a=ma)
            if mask is None:
                row['err_q'] = boot(peak_iter, q, NBOOT, m0)
                row['err_a'] = boot(peak_iter, a, NBOOT, m0a) if np.isfinite(m0a) else np.nan
                w = q[(q > (1 - WIN) * mq) & (q < (1 + WIN) * mq)] if np.isfinite(mq) else []
                row['n_pk'] = len(w)
            else:
                row['err_q'] = (1.2533 * 0.12 * mq / np.sqrt(max(len(q), 1))) if np.isfinite(mq) else np.nan
            rows.append(row)
    return pd.DataFrame(rows)


def fe55_gain_fit(P, col='peak_q', ecol='err_q', vref=420.0, vmin=385, min_pts=4):
    """Per pad: ln(peak) = a + b (V - vref) over the points with peak > 1.6 x floor."""
    rows = []
    for pid, g in P.groupby('channel_id'):
        g = g[(g.hv >= vmin) & np.isfinite(g[col]) & np.isfinite(g[ecol]) & (g[col] > 1.6 * FE_FLOOR)]
        if len(g) < min_pts:
            continue
        x = g.hv.to_numpy(float) - vref
        y = np.log(g[col].to_numpy()); sy = g[ecol].to_numpy() / g[col].to_numpy()
        sy = np.maximum(sy, 0.005)
        w = 1 / sy ** 2
        A = np.vstack([np.ones_like(x), x]).T
        C = np.linalg.inv(A.T @ (A * w[:, None]))
        p = C @ (A.T @ (w * y))
        chi2 = float(np.sum(w * (y - A @ p) ** 2)); ndf = len(x) - 2
        s = np.sqrt(max(chi2 / ndf, 1.0)) if ndf > 0 else 1.0
        rows.append(dict(channel_id=pid, lnG420=p[0], lnG420_err=np.sqrt(C[0, 0]) * s,
                         slope=p[1], slope_err=np.sqrt(C[1, 1]) * s, cov=C[0, 1] * s * s,
                         chi2=chi2, ndf=ndf, npts=len(x)))
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------- MIP
def tile_locate(xy, g, edge_min):
    """Pad tile containing each point + distance to its edge (convex quads)."""
    V = g['tile_verts']; pid = g['pad_id']
    tr = cKDTree(g['pad_xy'])
    _, i3 = tr.query(xy, k=3)
    out_pid = np.full(len(xy), -1, np.int64); out_d = np.full(len(xy), np.nan)
    for k in range(3):
        v = V[i3[:, k]]                                  # (n,4,2)
        a = v; b = np.roll(v, -1, axis=1)
        ab = b - a; ap = xy[:, None, :] - a
        cr = ab[..., 0] * ap[..., 1] - ab[..., 1] * ap[..., 0]
        inside = (cr >= 0).all(1) | (cr <= 0).all(1)
        t = np.clip((ap * ab).sum(-1) / (ab * ab).sum(-1), 0, 1)
        dd = np.hypot(*(ap - t[..., None] * ab).transpose(2, 0, 1)).min(1)
        new = inside & (out_pid < 0)
        out_pid[new] = pid[i3[new, k]]; out_d[new] = dd[new]
    return out_pid, out_d


def mip_extract(cfg, g, edge_min):
    d34 = cfg.out_dir('34_clustering_matching')
    A = np.load(f'{d34}/alignment.npz')
    T = cm.make_T(g, R=A['R'], s=A['s'], pad_mean=A['pad_mean'], m3_mean=A['m3_mean'])
    ra = pd.read_pickle(f'{d34}/rays_aligned.pkl').drop_duplicates('eventId').reset_index(drop=True)
    xy = inverse(T, ra[['x', 'y']].to_numpy())
    ra['xp'], ra['yp'] = xy[:, 0], xy[:, 1]
    ra['pid'], ra['de'] = tile_locate(xy, g, edge_min)
    sel = ra[(ra.pid >= 0) & (ra.de > edge_min) & (ra.sd > 2.0)].copy()
    ct = sc.channel_table(cfg)
    drop = set(p2io.drop_pads_for(cfg, ct))
    bad = set(int(e) for e in ps.SparkVeto.from_cfg(cfg).vetoed_ids_from_hits(
        cfg.combined_hits_dir, ct.attrs['feus'], min_amp=cfg.MIN_AMP))
    sel = sel[~sel.eventId.isin(bad)]
    evs = set(sel.eventId.tolist())
    br = ['eventId', 'channel', 'feu', 'amplitude', 'time', 'time_of_max', 'integral',
          'time_over_threshold', 'saturated', 'trunc_left', 'trunc_right']
    parts = []
    for df in p2io.iter_hits(cfg.combined_hits_dir, br, ct.attrs['feus'], progress=False,
                             t_max_h=cfg.T_MAX_H, t_min_h=cfg.T_MIN_H, min_amp=cfg.MIN_AMP):
        df = df[df.eventId.isin(evs)]
        h = pmap.attach_pads_to_hits(df, ct)
        h = h[h['mapped'] & h['pad_cx'].notna() & ~h['channel_id'].isin(drop)]
        parts.append(pd.DataFrame({
            'eventId': h.eventId.to_numpy(np.int64), 'channel_id': h.channel_id.to_numpy(np.int32),
            'x': h.pad_cx.to_numpy(np.float32), 'y': h.pad_cy.to_numpy(np.float32),
            'a': h.amplitude.clip(lower=0).to_numpy(np.float32), 't': h.time.to_numpy(np.float32),
            'tm': h.time_of_max.to_numpy(np.float32), 'itg': h.integral.to_numpy(np.float32),
            'tot': h.time_over_threshold.to_numpy(np.float32), 'sat': h.saturated.to_numpy(bool),
            'trunc': (h.trunc_left.astype(bool) | h.trunc_right.astype(bool)).to_numpy()}))
    H = pd.concat(parts, ignore_index=True)
    H = H.sort_values('a', ascending=False).drop_duplicates(['eventId', 'channel_id'])
    pxy = pd.DataFrame(g['pad_xy'], columns=['px', 'py'], index=g['pad_id'])
    S = sel[['eventId', 'pid', 'de', 'theta', 'sd']].join(pxy, on='pid')
    H = H.merge(S[['eventId', 'pid', 'px', 'py']], on='eventId')
    d = np.hypot(H.x - H.px, H.y - H.py)
    own = H[H.channel_id == H.pid].set_index('eventId')[['a', 't', 'tm', 'itg', 'tot', 'sat', 'trunc']]
    cl = H[d <= CLUS_R].groupby('eventId')
    nb = H[(d > 0.5) & (d <= CLUS_R)].groupby('eventId')
    S = S.set_index('eventId').join(own)
    S['q'] = cl.a.sum().reindex(S.index)
    S['n_nb'] = nb.size().reindex(S.index).fillna(0).astype(np.int16)
    S['a_nb'] = nb.a.max().reindex(S.index).fillna(0)
    S['n_ev_hits'] = H.groupby('eventId').size().reindex(S.index).fillna(0)
    return S.reset_index(), dict(n_rays=len(ra), n_in_tile=int((ra.pid >= 0).sum()),
                                 n_sel=len(sel), n_bad=len(bad), scale=float(T.s))


# ----------------------------------------------------------------------------- plotting
def draw_map(ax, verts_by_id, ids, vals, title, cmap='viridis', label='', vlim=None, fig=None):
    allv = np.array(list(verts_by_id.values()))
    ax.add_collection(PolyCollection(allv, facecolors='0.93', edgecolors='0.8', linewidths=0.2))
    ok = np.isfinite(vals)
    if ok.any():
        V = np.array([verts_by_id[i] for i in np.asarray(ids)[ok]])
        lo, hi = (np.nanpercentile(vals[ok], [2, 98]) if vlim is None else vlim)
        pc = PolyCollection(V, array=vals[ok], cmap=cmap, edgecolors='face', linewidths=0.2)
        pc.set_clim(lo, hi)
        ax.add_collection(pc)
        (fig or ax.figure).colorbar(pc, ax=ax, label=label, shrink=0.8)
    ax.autoscale_view(); ax.set_aspect('equal')
    ax.set_xlabel('pad x [mm]'); ax.set_ylabel('pad y [mm]'); ax.set_title(title, fontsize=9)


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--edge', type=float, default=2.5, help='min track distance inside the pad tile [mm]')
    ap.add_argument('--rebuild', action='store_true')
    args = ap.parse_args()

    cfgj = qa.get_config('det2_fe55scan1')
    hits_n, rays_n, g, cfg = sc.load('det2_long_repos')
    del hits_n, rays_n
    base = cfg.out_dir('46_deform_fe55_vs_mip')
    # non-default edge cut: same Fe55 cache, products in a sub-directory
    out = base if args.edge == 2.5 else cfg.out_dir('46_deform_fe55_vs_mip', f'edge{args.edge:g}')
    L = []                                  # summary lines

    def say(s=''):
        print(s, flush=True); L.append(s)

    # ---- channel tables: same channel_id -> pad in both periods
    ctj = pmap.build_channel_table(cfgj.run_config_path, cfgj.MAP_CSV_PATH, det_type=cfgj.DET_TYPE,
                                   det_name=cfgj.DET_NAME, strategy='reverse',
                                   drop_connectors=cfgj.DEAD_CONNECTORS, drop_channels=cfgj.DEAD_CHANNELS,
                                   cluster_link=cfgj.CLUSTER_LINK_MM, dead_guard=cfgj.DEAD_GUARD_MM,
                                   strategy_overrides=cfgj.STRATEGY_OVERRIDES)
    ctn = sc.channel_table(cfg)
    pj = ctj[ctj.mapped].drop_duplicates('channel_id')[['channel_id', 'pad_cx', 'pad_cy', 'connector_N', 'half', 'feu', 'channel']]
    pn = ctn[ctn.mapped].drop_duplicates('channel_id')[['channel_id', 'pad_cx', 'pad_cy', 'feu', 'channel']]
    dropn = set(p2io.drop_pads_for(cfg, ctn))
    pn = pn[~pn.channel_id.isin(dropn)]
    PM = pj.merge(pn, on='channel_id', suffixes=('', '_n'))
    dxy = float(np.max(np.hypot(PM.pad_cx - PM.pad_cx_n, PM.pad_cy - PM.pad_cy_n)))
    PM['block_j'] = PM.channel // 64; PM['block_n'] = PM.feu_n * 8 + PM.channel_n // 64
    reg = pd.read_csv(REGIONS_CSV)[['channel_id', 'reg8', 'reg24', 'reg48']]
    PM = PM.merge(reg, on='channel_id', how='left')
    say('46_deform_fe55_vs_mip -- det2 Fe55 (det2_fe55scan1, 18 Jul 2026) vs MIP (det2_long_repos, 25 Nov 2026)')
    say(f'channel tables: July {len(pj)} pads (FEU 3, conn 4-7), Nov {len(pn)} live pads (after cfg drops); '
        f'common {len(PM)}; max |pad centre difference| = {dxy:.3g} mm (channel_id -> pad identical)')
    say('electronics: July conn4-7 = FEU3 ch 0-511; Nov conn4 = FEU6 256-383, conn5 = FEU6 384-511, '
        'conn6 = FEU7 0-127, conn7 = FEU7 128-255 (same index modulo 128): every pad on a different channel/ASIC')

    # ---- Fe55 events
    fe_pkl = f'{base}/fe55_events.pkl'
    if args.rebuild or not os.path.isfile(fe_pkl):
        E, Est = fe55_extract(cfgj, ctj)
        E.to_pickle(fe_pkl); Est.to_csv(f'{base}/fe55_subruns.csv', index=False)
    E = pd.read_pickle(fe_pkl)
    hvs = sorted(E.hv.unique())
    say(f'Fe55 events: {len(E):,} over {len(hvs)} HV points {hvs[0]}-{hvs[-1]} V')

    # neighbour sharing in July, per HV
    say('\n[Fe55 charge on a 2nd pad] fraction of photopeak-band events (q in 0.75-1.25 x global peak) '
        'with >= 1 neighbour (centre within 15 mm) above threshold, and median a_nb/a_lead:')
    share_rows = []
    for hv in hvs:
        e = E[E.hv == hv]
        f = s18.fit_photopeak(e.q.to_numpy(), FE_FLOOR)
        if not f['ok'] or f['near_floor']:
            continue
        b = e[(e.q > 0.75 * f['mu']) & (e.q < 1.25 * f['mu'])]
        r = (b.a_nb / b.a)[b.n_nb > 0]
        share_rows.append(dict(hv=hv, peak_q=f['mu'], n=len(b), f_nb=(b.n_nb > 0).mean(),
                               f_nb1=(b.n_nb == 1).mean(), r_med=r.median(),
                               a_single=b.a[b.n_nb == 0].median(), q_shared=b.q[b.n_nb > 0].median()))
    SH = pd.DataFrame(share_rows)
    for _, r in SH.iterrows():
        say(f'  {r.hv:.0f} V: peak_q {r.peak_q:6.0f}  N {r.n:6.0f}  f(nb>0) {r.f_nb:.3f}  f(nb=1) {r.f_nb1:.3f}  '
            f'median a_nb/a_lead {r.r_med:.3f}')

    # ---- per-pad Fe55 peaks + gain-vs-V fits
    P = fe55_pad_peaks(E, hvs)
    P.to_csv(f'{out}/fe55_pad_peaks.csv', index=False)
    G = fe55_gain_fit(P, 'peak_q', 'err_q')
    Ga = fe55_gain_fit(P, 'peak_a', 'err_a').rename(columns={'lnG420': 'lnA420', 'lnG420_err': 'lnA420_err',
                                                             'slope': 'slope_a', 'slope_err': 'slope_a_err'})
    Gt = fe55_gain_fit(P.assign(peak_q=P.peak_q_tm), 'peak_q', 'err_q')[['channel_id', 'lnG420']].rename(
        columns={'lnG420': 'lnG420_tm'})
    # split halves (odd/even eventId) for the reliability of the map
    Po = fe55_pad_peaks(E, hvs, E.eventId % 2 == 1); Pe = fe55_pad_peaks(E, hvs, E.eventId % 2 == 0)
    Go = fe55_gain_fit(Po)[['channel_id', 'lnG420', 'slope']].rename(columns={'lnG420': 'lnG_odd', 'slope': 'b_odd'})
    Ge = fe55_gain_fit(Pe)[['channel_id', 'lnG420', 'slope']].rename(columns={'lnG420': 'lnG_even', 'slope': 'b_even'})
    p420 = P[P.hv == 420][['channel_id', 'n', 'n_pk', 'peak_q', 'err_q', 'peak_a']].rename(
        columns={'n': 'fe_n420', 'n_pk': 'fe_npk420', 'peak_q': 'fe_q420', 'err_q': 'fe_q420_err', 'peak_a': 'fe_a420'})
    FE = (G.merge(Ga[['channel_id', 'lnA420', 'lnA420_err', 'slope_a']], on='channel_id', how='left')
           .merge(Gt, on='channel_id', how='left').merge(Go, on='channel_id', how='left')
           .merge(Ge, on='channel_id', how='left').merge(p420, on='channel_id', how='left'))
    say(f'\n[Fe55 per pad] estimator: events assigned to their LEADING pad; cluster charge q = lead + pads '
        f'with centre within {CLUS_R:g} mm (the 2nd-pad charge above is kept); per (pad, HV) photopeak = median '
        f'in [0.6,1.6] x global peak, then 4 x median in +-{WIN:.0%} (excludes the Ar escape peak at ~0.5 x); '
        f'bootstrap error ({NBOOT} resamples). ln(peak) vs V fitted per pad over V >= 385 with peak > {1.6*FE_FLOOR:.0f} '
        f'ADC (threshold/trigger censoring), >= 4 points; lnG420 = fitted value at 420 V.')
    n420 = P[P.hv == 420]
    say(f'  420 V: {len(n420)} pads with an estimate; events per pad (q > {FE_FLOOR:.0f}) median {n420.n.median():.0f} '
        f'(10-90 %: {n420.n.quantile(.1):.0f}-{n420.n.quantile(.9):.0f}); per-pad peak error median '
        f'{np.nanmedian(n420.err_q / n420.peak_q):.3f} (relative)')
    say(f'  per-pad fit: {len(G)} pads; lnG420 error median {G.lnG420_err.median():.4f}; '
        f'map spread sd(lnG420) = {G.lnG420.std():.3f}; chi2/ndf median {np.nanmedian(G.chi2 / G.ndf.clip(lower=1)):.2f}')
    say(f'  slope d lnG/dV: median {G.slope.median():.5f} /V (x2 every {np.log(2) / G.slope.median():.1f} V); '
        f'per-pad error median {G.slope_err.median():.5f}; sd over pads {G.slope.std():.5f}; '
        f'expected sd from errors {np.sqrt(np.mean(G.slope_err ** 2)):.5f}')
    ro, rel = reliability(FE.lnG_odd.to_numpy(), FE.lnG_even.to_numpy())
    rbo, relb = reliability(FE.b_odd.to_numpy(), FE.b_even.to_numpy())
    say(f'  split-half (odd/even eventId) reliability: lnG420 r = {ro:.3f} (full-sample {rel:.3f}); '
        f'slope r = {rbo:.3f} (full {relb:.3f})')
    rr = np.corrcoef(*FE[['lnG420', 'lnA420']].dropna().to_numpy().T)[0, 1]
    rt = np.corrcoef(*FE[['lnG420', 'lnG420_tm']].dropna().to_numpy().T)[0, 1]
    say(f'  estimator comparison: corr(lnG cluster q, lnG lead-pad a) = {rr:.3f}; offset median '
        f'{np.nanmedian(FE.lnG420 - FE.lnA420):.3f}; corr(median, window-mean estimator) = {rt:.3f}')

    # Fe55 shape at 420 V, photopeak events of each pad (unsaturated)
    pk = FE.set_index('channel_id').fe_q420
    e = E[E.hv == 420].join(pk.rename('pk'), on='channel_id')
    e = e[np.isfinite(e.pk) & (e.q > (1 - WIN) * e.pk) & (e.q < (1 + WIN) * e.pk) & ~e.sat]
    e = e.assign(rise=e.tm - e.t, width=TPS * e.itg / e.a, totn=e.tot)
    e = e[(e.rise > 0) & (e.rise < 600)]
    # amplitude dependence of the Fe55 proxies: within pads, over 395-420 V photopeak events
    ee = E[E.hv >= 395].copy()
    pkh = P.set_index(['channel_id', 'hv']).peak_q
    ee['pk'] = pkh.reindex(pd.MultiIndex.from_arrays([ee.channel_id, ee.hv])).to_numpy()
    ee = ee[np.isfinite(ee.pk) & (ee.q > (1 - WIN) * ee.pk) & (ee.q < (1 + WIN) * ee.pk) & ~ee.sat]
    ee = ee.assign(rise=ee.tm - ee.t, width=TPS * ee.itg / ee.a, totn=ee.tot, la=np.log(ee.a))
    ee = ee[(ee.rise > 0) & (ee.rise < 600)]
    say('\n[Fe55 shape] lead-pad hit of photopeak events (q within +-25 % of the pad peak), unsaturated, '
        '0 < rise < 600 ns.  Amplitude dependence within pads (pad-demeaned, 395-420 V pooled):')
    sF = {}
    for c in ['rise', 'width', 'totn']:
        y = demean_by(ee[c].to_numpy(), ee.channel_id); x = demean_by(ee.la.to_numpy(), ee.channel_id)
        k = np.polyfit(x, y, 1)[0]
        kb = [np.polyfit(x[i], y[i], 1)[0] for i in RNG.integers(0, len(x), (40, len(x)))]
        sF[c] = k
        med = ee.groupby('hv')[c].median()
        say(f'  {c:5s}: d/d ln(amp) = {k:7.2f} +- {np.std(kb):.2f} ns; median by HV '
            + ', '.join(f'{h}:{v:.1f}' for h, v in med.items()))
    fsh = e.groupby('channel_id').agg(fe_nsh=('rise', 'size'), fe_rise=('rise', tmean10),
                                      fe_width=('width', tmean10), fe_tot=('totn', 'mean'),
                                      fe_ash=('a', 'median'))
    for c, cc in [('rise', 'fe_rise'), ('width', 'fe_width'), ('totn', 'fe_tot')]:
        fsh[cc + '_err'] = e.groupby('channel_id')[c].apply(
            lambda v: boot(tmean10 if c != 'totn' else np.mean, v.to_numpy(), 60))
    # single-pad vs shared lead-pad rise (sharing check)
    rs1 = e[e.n_nb == 0].rise.median(); rs2 = e[e.n_nb > 0].rise.median()
    say(f'  420 V per pad: photopeak hits median {fsh.fe_nsh.median():.0f}; rise trimmed mean: map median '
        f'{fsh.fe_rise.median():.1f} ns, per-pad error median {fsh.fe_rise_err.median():.2f} ns, sd over pads '
        f'{fsh.fe_rise.std():.2f} ns; lead rise single-pad events {rs1:.1f} vs shared {rs2:.1f} ns')
    ro_ = e[e.eventId % 2 == 1].groupby('channel_id').rise.apply(tmean10)
    re_ = e[e.eventId % 2 == 0].groupby('channel_id').rise.apply(tmean10)
    jj = pd.concat([ro_, re_], axis=1).dropna().to_numpy()
    r1, r2 = reliability(jj[:, 0], jj[:, 1])
    say(f'  Fe55 rise map split-half (odd/even) r = {r1:.3f} (full {r2:.3f})')
    fsh['fe_rise_odd'] = ro_; fsh['fe_rise_even'] = re_
    FE = FE.merge(fsh.reset_index(), on='channel_id', how='left')

    # ---- MIP
    mip_pkl = f'{base}/mip_tracks_edge{args.edge:g}.pkl'
    if args.rebuild or not os.path.isfile(mip_pkl):
        S, mst = mip_extract(cfg, g, args.edge)
        S.to_pickle(mip_pkl)
        pd.Series(mst).to_csv(f'{base}/mip_extract_edge{args.edge:g}.csv')
    S = pd.read_pickle(mip_pkl)
    mst = pd.read_csv(f'{base}/mip_extract_edge{args.edge:g}.csv', index_col=0).iloc[:, 0]
    S['cth'] = np.cos(np.radians(S.theta))
    S['hit'] = S.a.notna() & (S.a > 0)
    S['ac'] = S.a * S.cth; S['qc'] = S.q * S.cth
    S['rise'] = S.tm - S.t; S['width'] = TPS * S.itg / S.a; S['totn'] = S.tot
    half = S.eventId < S.eventId.median()
    say(f'\n[MIP selection] aligned rays {int(mst.n_rays):,}; in a live tile {int(mst.n_in_tile):,}; '
        f'> {args.edge:g} mm inside the tile and sd > 2 mm (live area): {int(mst.n_sel):,}; minus spark/burst-vetoed '
        f'events: {len(S):,}.  Own-pad hit (any amplitude): {S.hit.mean():.3f}.')
    say(f'  charge on a neighbour (centre within 15 mm) for own-pad hits: f(n_nb>0) = {(S[S.hit].n_nb > 0).mean():.4f}; '
        f'saturated own hits {S[S.hit].sat.mean():.3f}; truncated {S[S.hit].trunc.mean():.3f}; '
        f'median theta {S.theta.median():.1f} deg')
    Sh = S[S.hit]
    grp = Sh.groupby('pid')
    MI = pd.DataFrame({'mip_ntrk': S.groupby('pid').size(), 'mip_nhit': grp.size(),
                       'mip_med': grp.ac.median(), 'mip_tm70': grp.ac.apply(lambda v: tmean_low(v.to_numpy())),
                       'mip_qmed': grp.qc.median(),
                       'mip_med_h1': Sh[half[Sh.index]].groupby('pid').ac.median(),
                       'mip_med_h2': Sh[~half[Sh.index]].groupby('pid').ac.median(),
                       'mip_sat': grp.sat.mean()})
    MI['mip_eff'] = MI.mip_nhit / MI.mip_ntrk
    MI['mip_med_err'] = grp.ac.apply(lambda v: boot(np.median, v.to_numpy(), 200))
    MI['mip_tm70_err'] = grp.ac.apply(lambda v: boot(tmean_low, v.to_numpy(), 200))
    MI['mip_qmed_err'] = grp.qc.apply(lambda v: boot(np.median, v.to_numpy(), 200))
    MI.index.name = 'channel_id'
    rmh, relm = reliability(np.log(MI.mip_med_h1.to_numpy()), np.log(MI.mip_med_h2.to_numpy()))
    say(f'[MIP per pad] MPV estimators on amp x cos(theta) of the own-pad hit: median (primary), mean of the '
        f'lowest 70 %, and median of the summed cluster (pads within 15 mm); bootstrap errors (200).')
    say(f'  pads with tracks {len(MI)}; hits per pad median {MI.mip_nhit.median():.0f} (10-90 %: '
        f'{MI.mip_nhit.quantile(.1):.0f}-{MI.mip_nhit.quantile(.9):.0f}); relative error of the median: median '
        f'{np.nanmedian(MI.mip_med_err / MI.mip_med):.3f}; of the low-70 % mean {np.nanmedian(MI.mip_tm70_err / MI.mip_tm70):.3f}')
    say(f'  split-half (first/second half of the run) reliability of ln median: r = {rmh:.3f} (full {relm:.3f}); '
        f'global median first/second half: {Sh[half[Sh.index]].ac.median():.0f} / {Sh[~half[Sh.index]].ac.median():.0f} ADC')
    rmt = np.corrcoef(np.log(MI[['mip_med', 'mip_tm70']].dropna().to_numpy().T))[0, 1]
    rmq = np.corrcoef(np.log(MI[['mip_med', 'mip_qmed']].dropna().to_numpy().T))[0, 1]
    say(f'  estimator comparison: corr(ln median, ln low-70 % mean) = {rmt:.3f}; corr(ln median, ln cluster median) = {rmq:.3f}')

    # MIP shape: own hit, unsaturated, untruncated; amplitude dependence within pads
    Ms = Sh[~Sh.sat & ~Sh.trunc & (Sh.rise > 0) & (Sh.rise < 600) & (Sh.a > 100)].copy()
    Ms['la'] = np.log(Ms.a)
    sM = {}
    say('[MIP shape] own-pad hit, unsaturated, untruncated, amp > 100 ADC. Amplitude dependence within pads:')
    for c in ['rise', 'width', 'totn']:
        y = demean_by(Ms[c].to_numpy(), Ms.pid); x = demean_by(Ms.la.to_numpy(), Ms.pid)
        k = np.polyfit(x, y, 1)[0]
        kb = [np.polyfit(x[i], y[i], 1)[0] for i in RNG.integers(0, len(x), (40, len(x)))]
        sM[c] = k
        bins = pd.cut(Ms.a, [100, 300, 600, 1000, 1500, 2500, 4000])
        med = Ms.groupby(bins, observed=True)[c].median()
        say(f'  {c:5s}: d/d ln(amp) = {k:7.2f} +- {np.std(kb):.2f} ns; median by amp bin '
            + ', '.join(f'{iv.left:.0f}-{iv.right:.0f}:{v:.1f}' for iv, v in med.items()))
    say(f'  theta dependence of rise (median, theta < 8 / 8-15 / > 15 deg): '
        + ' / '.join(f'{Ms[m].rise.median():.1f}' for m in [Ms.theta < 8, (Ms.theta >= 8) & (Ms.theta < 15), Ms.theta >= 15]))
    # Fe55-matched amplitude band per reg24 (16-84 % of Fe55 photopeak lead amplitudes at 420 V)
    e24 = e.merge(PM[['channel_id', 'reg24']], on='channel_id')
    band = e24.groupby('reg24').a.quantile([0.16, 0.84]).unstack()
    Ms = Ms.merge(PM[['channel_id', 'reg24']], left_on='pid', right_on='channel_id', how='left')
    lo = band[0.16].reindex(Ms.reg24).to_numpy(); hi = band[0.84].reindex(Ms.reg24).to_numpy()
    Ms['inband'] = (Ms.a >= lo) & (Ms.a <= hi)
    # amplitude-corrected proxies (each period's own within-pad slope) at A_ref = 1000 ADC
    for c in ['rise', 'width', 'totn']:
        Ms[c + '_c'] = Ms[c] - sM[c] * np.log(Ms.a / 1000.0)
        e[c + '_c'] = e[c] - sF[c] * np.log(e.a / 1000.0)
    gm = Ms.groupby('pid'); gb = Ms[Ms.inband].groupby('pid')
    MS = pd.DataFrame({'mip_nsh': gm.size(), 'mip_rise': gm.rise.apply(tmean10),
                       'mip_rise_c': gm.rise_c.apply(tmean10), 'mip_width': gm.width.apply(tmean10),
                       'mip_width_c': gm.width_c.apply(tmean10), 'mip_tot': gm.totn.apply(tmean10),
                       'mip_nband': gb.size(), 'mip_rise_band': gb.rise.apply(tmean10)})
    MS['mip_rise_err'] = gm.rise.apply(lambda v: boot(tmean10, v.to_numpy(), 60))
    MS['mip_rise_c_err'] = gm.rise_c.apply(lambda v: boot(tmean10, v.to_numpy(), 60))
    MS['mip_width_err'] = gm.width.apply(lambda v: boot(tmean10, v.to_numpy(), 60))
    MS['mip_tot_err'] = gm.totn.apply(lambda v: boot(tmean10, v.to_numpy(), 60))
    MS['mip_rise_band_err'] = gb.rise.apply(lambda v: boot(tmean10, v.to_numpy(), 60))
    # alternative correction: the ELECTRONIC amplitude slope measured with Fe55 (point charges)
    Ms['rise_cF'] = Ms.rise - sF['rise'] * np.log(Ms.a / 1000.0)
    MS['mip_rise_cF'] = Ms.groupby('pid').rise_cF.apply(tmean10)
    MS['mip_rise_cF_err'] = MS.mip_rise_err
    h1 = Ms.eventId < S.eventId.median()
    MS['mip_rise_h1'] = Ms[h1].groupby('pid').rise.apply(tmean10)
    MS['mip_rise_h2'] = Ms[~h1].groupby('pid').rise.apply(tmean10)
    r1, r2 = reliability(MS.mip_rise_h1.to_numpy(), MS.mip_rise_h2.to_numpy())
    say(f'  MIP rise map split-half (first/second half of run) r = {r1:.3f} (full {r2:.3f})')
    MS.index.name = 'channel_id'
    fc = e.groupby('channel_id').agg(fe_rise_c=('rise_c', tmean10), fe_width_c=('width_c', tmean10))
    say(f'  per pad: shape hits median {MS.mip_nsh.median():.0f}, rise map median {MS.mip_rise.median():.1f} ns, '
        f'per-pad error median {MS.mip_rise_err.median():.2f} ns; in the Fe55-matched band: hits per pad median '
        f'{MS.mip_nband.median():.0f}; band edges per reg24 (ADC) median {band[0.16].median():.0f}-{band[0.84].median():.0f}')

    # ---- merge per pad
    D = (PM.merge(FE, on='channel_id', how='left').merge(MI.reset_index(), on='channel_id', how='left')
           .merge(MS.reset_index(), on='channel_id', how='left').merge(fc.reset_index(), on='channel_id', how='left'))
    D['lnM'] = np.log(D.mip_med); D['lnM_err'] = D.mip_med_err / D.mip_med
    D['lnM70'] = np.log(D.mip_tm70); D['lnM70_err'] = D.mip_tm70_err / D.mip_tm70
    D['lnMq'] = np.log(D.mip_qmed); D['lnMq_err'] = D.mip_qmed_err / D.mip_qmed
    D['lnG430'] = D.lnG420 + 10 * D.slope
    D['lnG430_err'] = np.sqrt(D.lnG420_err ** 2 + 100 * D.slope_err ** 2 + 20 * D['cov'])
    D['lnR'] = D.lnM - D.lnG420; D['lnR_err'] = np.hypot(D.lnM_err, D.lnG420_err)
    D['drise'] = D.mip_rise - D.fe_rise; D['drise_err'] = np.hypot(D.mip_rise_err, D.fe_rise_err)
    D['drise_c'] = D.mip_rise_c - D.fe_rise_c
    D['dwidth'] = D.mip_width - D.fe_width; D['dwidth_err'] = np.hypot(D.mip_width_err, D.fe_width_err)
    D['dtot'] = D.mip_tot - D.fe_tot; D['dtot_err'] = np.hypot(D.mip_tot_err, D.fe_tot_err)
    D.to_csv(f'{out}/pads_fe55_vs_mip.csv', index=False)
    ok = np.isfinite(D.lnM) & np.isfinite(D.lnG420)
    say(f'\n[common pads] {len(D)} common live pads; with Fe55 gain fit {int(np.isfinite(D.lnG420).sum())}; with MIP MPV '
        f'{int(np.isfinite(D.lnM).sum())}; both {int(ok.sum())}; MIP hits/pad on these median {D.mip_nhit[ok].median():.0f}')

    # ---- 4a: ln MIP vs ln Fe55
    say('\n[4a] ln(MIP MPV) vs ln(Fe55 gain), per pad (ODR with per-pad errors on both axes; Pearson r with '
        'bootstrap error; r_max = sqrt(rel_M rel_G) is the ceiling set by measurement noise):')
    rel_noise = {}
    for nm, v, ev in [('lnG420', D.lnG420, D.lnG420_err), ('lnM', D.lnM, D.lnM_err)]:
        vv = v[ok]; rel_noise[nm] = 1 - np.nanmean(ev[ok] ** 2) / np.nanvar(vv)
    rmax = np.sqrt(max(rel_noise['lnG420'], 0) * max(rel_noise['lnM'], 0))
    say(f'  error-based reliabilities: Fe55 {rel_noise["lnG420"]:.3f}, MIP {rel_noise["lnM"]:.3f} -> r_max = {rmax:.3f}; '
        f'sd(lnG420) = {D.lnG420[ok].std():.3f}, sd(lnM) = {D.lnM[ok].std():.3f}, median errors {D.lnG420_err[ok].median():.3f} / '
        f'{D.lnM_err[ok].median():.3f}')
    res4a = {}
    for lab, xg, xge, ym, yme in [('median vs lnG420', 'lnG420', 'lnG420_err', 'lnM', 'lnM_err'),
                                  ('median vs lnG430 (extrap.)', 'lnG430', 'lnG430_err', 'lnM', 'lnM_err'),
                                  ('low-70% mean vs lnG420', 'lnG420', 'lnG420_err', 'lnM70', 'lnM70_err'),
                                  ('cluster median vs lnG420', 'lnG420', 'lnG420_err', 'lnMq', 'lnMq_err'),
                                  ('median vs lead-pad lnA420', 'lnA420', 'lnA420_err', 'lnM', 'lnM_err')]:
        for dm in ['none', 'connector', 'half']:
            x = D[xg].to_numpy(); y = D[ym].to_numpy()
            if dm != 'none':
                key = D.connector_N if dm == 'connector' else D.connector_N * 2 + (D.half == 'top')
                x = demean_by(x, key); y = demean_by(y, key)
            k, ke, _ = odr_fit(x, y, D[xge].to_numpy(), D[yme].to_numpy())
            r, re_, n = boot_corr(x, y, 500)
            okk = np.isfinite(x) & np.isfinite(y)
            kols = np.polyfit(x[okk], y[okk], 1)[0]
            rs = stats.spearmanr(x[okk], y[okk]).correlation
            scat = np.std(y[okk] - (k * x[okk] + np.mean(y[okk] - k * x[okk])))
            res4a[(lab, dm)] = (r, re_, k, ke)
            say(f'  {lab:28s} demean {dm:9s}: N {n}  r = {r:.3f} +- {re_:.3f} (Spearman {rs:.3f})  ODR slope = '
                f'{k:.3f} +- {ke:.3f}  OLS slope {kols:.3f}  scatter {scat:.3f}')
    # region level
    say('  region level (weighted means of per-pad values within the region):')
    RG = {}
    for rk in ['reg24', 'reg48']:
        rows = []
        for rid, gq in D[ok].groupby(rk):
            mg, sg, _, n = wmean(gq.lnG420.to_numpy(), gq.lnG420_err.to_numpy())
            mm, sm, _, _ = wmean(gq.lnM.to_numpy(), gq.lnM_err.to_numpy())
            mr, sr, c2, _ = wmean(gq.lnR.to_numpy(), gq.lnR_err.to_numpy())
            md, sdd, c2d, nd = wmean(gq.drise.to_numpy(), gq.drise_err.to_numpy())
            mfr, sfr, _, _ = wmean(gq.fe_rise.to_numpy(), gq.fe_rise_err.to_numpy())
            mmr, smr, _, _ = wmean(gq.mip_rise.to_numpy(), gq.mip_rise_err.to_numpy())
            rows.append(dict(region=rid, n=n, cx=gq.pad_cx.mean(), cy=gq.pad_cy.mean(),
                             conn=int(gq.connector_N.mode().iloc[0]),
                             lnG=mg, lnG_err=sg, lnM=mm, lnM_err=sm, lnR=mr, lnR_err=sr, lnR_chi2_in=c2,
                             fe_rise=mfr, fe_rise_err=sfr, mip_rise=mmr, mip_rise_err=smr,
                             drise=md, drise_err=sdd, drise_chi2_in=c2d))
        R = pd.DataFrame(rows); RG[rk] = R
        R.to_csv(f'{out}/regions_{rk}.csv', index=False)
        k, ke, _ = odr_fit(R.lnG.to_numpy(), R.lnM.to_numpy(), R.lnG_err.to_numpy(), R.lnM_err.to_numpy())
        r, re_, n = boot_corr(R.lnG.to_numpy(), R.lnM.to_numpy(), 2000)
        xd = demean_by(R.lnG.to_numpy(), R.conn); yd = demean_by(R.lnM.to_numpy(), R.conn)
        kd, kde, _ = odr_fit(xd, yd, R.lnG_err.to_numpy(), R.lnM_err.to_numpy())
        rd, rde, _ = boot_corr(xd, yd, 2000)
        say(f'  {rk}: N {n} regions ({R.n.median():.0f} pads each), r = {r:.3f} +- {re_:.3f}, ODR slope {k:.3f} +- {ke:.3f}; '
            f'connector-demeaned r = {rd:.3f} +- {rde:.3f}, slope {kd:.3f} +- {kde:.3f}')

    # ---- 4b: ratio map
    say('\n[4b] ln(MIP MPV / Fe55 gain) = lnM - lnG420 (absolute value meaningless: different electronics, HV, '
        'ADC scale; the MAP is the observable):')
    v = D.lnR.to_numpy(); ev = D.lnR_err.to_numpy(); okr = np.isfinite(v) & np.isfinite(ev)
    m, s, c2, n = wmean(v, ev)
    say(f'  per pad: N {n}, weighted mean {m:.3f}, sd {np.nanstd(v[okr]):.3f}, median error {np.nanmedian(ev[okr]):.3f}; '
        f'chi2/ndf vs flat = {c2:.1f}/{n - 1} = {c2 / (n - 1):.2f} (p = {stats.chi2.sf(c2, n - 1):.2g})')
    for key, nm in [(D.connector_N, 'connector'), (D.connector_N * 2 + (D.half == 'top'), 'half (64-ch block)')]:
        c2c = 0; nn = 0; ngrp = 0
        for _, gq in D[okr].groupby(key[okr]):
            mm_, _, cc, k_ = wmean(gq.lnR.to_numpy(), gq.lnR_err.to_numpy()); c2c += cc; nn += k_; ngrp += 1
        say(f'  after per-{nm} constants: chi2/ndf = {c2c:.1f}/{nn - ngrp} = {c2c / (nn - ngrp):.2f}')
    for c in sorted(D.connector_N.unique()):
        gq = D[okr & (D.connector_N == c)]
        mm_, ss_, _, k_ = wmean(gq.lnR.to_numpy(), gq.lnR_err.to_numpy())
        say(f'    conn {c}: lnR = {mm_:.3f} +- {ss_:.3f} (N {k_})')
    I, pI = moran(v, D[['pad_cx', 'pad_cy']].to_numpy())
    Ic, pIc = moran(demean_by(v, D.connector_N), D[['pad_cx', 'pad_cy']].to_numpy())
    say(f'  spatial autocorrelation (Moran I, neighbours within 15 mm): raw I = {I:.3f} (p = {pI:.3g}); '
        f'connector-demeaned I = {Ic:.3f} (p = {pIc:.3g})')
    # smooth-structure fit: plane / quadratic in (x, y) after connector constants
    xy = D[['pad_cx', 'pad_cy']].to_numpy(); xc = xy - np.nanmean(xy[okr], 0)
    Xq = np.column_stack([xc[:, 0], xc[:, 1], xc[:, 0] ** 2, xc[:, 1] ** 2, xc[:, 0] * xc[:, 1]]) / 100.0
    Cn = pd.get_dummies(D.connector_N).to_numpy(float)
    for nm, X in [('plane', np.column_stack([Cn, Xq[:, :2]])), ('quadratic', np.column_stack([Cn, Xq]))]:
        w = 1 / ev[okr] ** 2
        Xw = X[okr] * np.sqrt(w)[:, None]; yw = v[okr] * np.sqrt(w)
        beta, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
        c2f = float(np.sum((yw - Xw @ beta) ** 2))
        c20 = sum(wmean(gq.lnR.to_numpy(), gq.lnR_err.to_numpy())[2] for _, gq in D[okr].groupby(D.connector_N[okr]))
        dof = X.shape[1] - Cn.shape[1]
        say(f'  connector consts + {nm} in (x, y): chi2 {c20:.1f} -> {c2f:.1f} for {dof} extra params '
            f'(delta chi2 p = {stats.chi2.sf(c20 - c2f, dof):.2g})')
    for rk in ['reg24', 'reg48']:
        R = RG[rk]
        mR, sR, c2R, nR = wmean(R.lnR.to_numpy(), R.lnR_err.to_numpy())
        c2Rc = sum(wmean(gq.lnR.to_numpy(), gq.lnR_err.to_numpy())[2] for _, gq in R.groupby('conn'))
        ngc = R.conn.nunique()
        say(f'  {rk}: chi2/ndf vs flat {c2R:.1f}/{nR - 1}; with per-connector constants {c2Rc:.1f}/{nR - ngc}; '
            f'region values (lnR +- err): ' + ', '.join(f'{int(a)}:{b:.3f}+-{c:.3f}' for a, b, c in
                                                      zip(R.region, R.lnR, R.lnR_err)))

    # power-law residual: is there structure beyond lnM = k lnG + c ?
    kpl, kple, _ = odr_fit(D.lnG420.to_numpy(), D.lnM.to_numpy(), D.lnG420_err.to_numpy(), D.lnM_err.to_numpy())
    D['lnRk'] = D.lnM - kpl * D.lnG420
    v = D.lnRk.to_numpy(); ev = D.lnR_err.to_numpy()
    c2k = sum(wmean(gq.lnRk.to_numpy(), gq.lnR_err.to_numpy())[2] for _, gq in D[okr].groupby(D.connector_N[okr]))
    Ik, pIk = moran(demean_by(v, D.connector_N), D[['pad_cx', 'pad_cy']].to_numpy())
    say(f'  residual of lnM after lnM = k lnG (k = {kpl:.3f} +- {kple:.3f}) and connector constants: chi2/ndf = '
        f'{c2k:.1f}/{int(okr.sum()) - 5} = {c2k / (okr.sum() - 5):.2f}; Moran I = {Ik:.3f} (p = {pIk:.3g})')

    # ---- 4c: shapes
    hkey = (D.connector_N * 2 + (D.half == 'top')).to_numpy()
    say('\n[4c] pulse shape vs Fe55 gain (ODR with errors on both axes; ns per e-fold of Fe55 gain). Demeaning by '
        'connector, and by connector half = 64-channel DREAM ASIC (the halves are the same channel blocks in July '
        'and November, so half-demeaning removes per-ASIC constants of BOTH periods):')
    shp = [('Fe55 rise', 'fe_rise', 'fe_rise_err'), ('Fe55 width', 'fe_width', 'fe_width_err'),
           ('Fe55 TOT', 'fe_tot', 'fe_tot_err'), ('MIP rise', 'mip_rise', 'mip_rise_err'),
           ('MIP rise amp-corr(MIP)', 'mip_rise_c', 'mip_rise_c_err'),
           ('MIP rise amp-corr(Fe)', 'mip_rise_cF', 'mip_rise_cF_err'),
           ('MIP rise Fe55-band', 'mip_rise_band', 'mip_rise_band_err'),
           ('MIP width', 'mip_width', 'mip_width_err'), ('MIP TOT', 'mip_tot', 'mip_tot_err'),
           ('MIP-Fe55 rise', 'drise', 'drise_err'), ('lnMPV - lnG', 'lnR', 'lnR_err')]
    for dm, key in [('connector', D.connector_N.to_numpy()), ('half', hkey)]:
        say(f'  demeaned by {dm}:')
        for lab, yc, yce in shp:
            x = demean_by(D.lnG420.to_numpy(), key); y = demean_by(D[yc].to_numpy(), key)
            k, ke, _ = odr_fit(x, y, D.lnG420_err.to_numpy(), D[yce].to_numpy())
            r, re_, n = boot_corr(x, y, 500)
            okk = np.isfinite(x) & np.isfinite(y)
            kols = np.polyfit(x[okk], y[okk], 1)[0] if okk.sum() > 3 else np.nan
            say(f'    {lab:24s}: N {n}  slope {k:8.3f} +- {ke:.3f} (OLS {kols:7.3f})  r = {r:.3f} +- {re_:.3f}')
    say('  MIP rise vs Fe55 rise per pad:')
    for dm, key in [('connector', D.connector_N.to_numpy()), ('half', hkey)]:
        x = demean_by(D.fe_rise.to_numpy(), key); y = demean_by(D.mip_rise.to_numpy(), key)
        k, ke, _ = odr_fit(x, y, D.fe_rise_err.to_numpy(), D.mip_rise_err.to_numpy())
        r, re_, n = boot_corr(x, y, 500)
        say(f'    {dm}-demeaned N {n}: r = {r:.3f} +- {re_:.3f}, ODR slope {k:.3f} +- {ke:.3f}')
    say('  difference maps (MIP - Fe55): chi2 vs flat / per-connector / per-half constants; Moran I after connector and '
        'after half demeaning:')
    for lab, c, ce in [('rise', 'drise', 'drise_err'), ('width (defs differ)', 'dwidth', 'dwidth_err'),
                       ('TOT (defs differ)', 'dtot', 'dtot_err')]:
        v = D[c].to_numpy(); ev = D[ce].to_numpy(); okd = np.isfinite(v) & np.isfinite(ev) & (ev > 0)
        m, s, c2, n = wmean(v, ev)
        c2c = sum(wmean(gq[c].to_numpy(), gq[ce].to_numpy())[2] for _, gq in D[okd].groupby(D.connector_N[okd]))
        c2h = sum(wmean(gq[c].to_numpy(), gq[ce].to_numpy())[2] for _, gq in D[okd].groupby(hkey[okd]))
        Ic, pIc = moran(demean_by(v, D.connector_N), D[['pad_cx', 'pad_cy']].to_numpy())
        Ih, pIh = moran(demean_by(v, hkey), D[['pad_cx', 'pad_cy']].to_numpy())
        say(f'    {lab:20s}: mean {m:.2f} +- {s:.2f} ns, sd {np.nanstd(v[okd]):.2f}, median err {np.nanmedian(ev[okd]):.2f}; '
            f'chi2/ndf flat {c2 / (n - 1):.2f}, conn {c2c / (n - 4):.2f}, half {c2h / (n - 8):.2f} (ndf ~{n - 8}); '
            f'Moran I conn {Ic:.3f} (p {pIc:.3g}), half {Ih:.3f} (p {pIh:.3g})')
    R = RG['reg24']
    say('    reg24 rise (weighted means): ' + ', '.join(
        f'{int(a)}: Fe {b:.1f}+-{c:.1f} MIP {d:.1f}+-{e_:.1f}' for a, b, c, d, e_ in
        zip(R.region, R.fe_rise, R.fe_rise_err, R.mip_rise, R.mip_rise_err)))
    for rk in ['reg24', 'reg48']:
        R = RG[rk]
        r, re_, n = boot_corr(R.fe_rise.to_numpy(), R.mip_rise.to_numpy(), 2000)
        r2_, re2, _ = boot_corr(R.lnR.to_numpy(), R.drise.to_numpy(), 2000)
        k, ke, _ = odr_fit(R.lnR.to_numpy(), R.drise.to_numpy(), R.lnR_err.to_numpy(), R.drise_err.to_numpy())
        say(f'    {rk}: corr(Fe55 rise, MIP rise) = {r:.3f} +- {re_:.3f}; corr(lnR, MIP-Fe55 rise) = {r2_:.3f} +- {re2:.3f}, '
            f'ODR d(drise)/d(lnR) = {k:.1f} +- {ke:.1f} ns (N {n})')

    # ---- drift-gap signature: lnR and (MIP - Fe55) rise move together?
    say('\n[drift-gap test] a thicker drift gap raises the MIP charge (lnR up by dg/4 mm) and stretches the MIP pulse '
        '(rise up by ~ dg / v_drift x O(0.5); Ar/iC4H10 95/5 at 750 V/cm, v ~ 4 cm/us -> O(10) ns per mm), i.e. a '
        'POSITIVE correlation, expected d(drise)/d(lnR) of order +40 ns (factor-2 uncertain).  Per pad:')
    for dm, key in [('connector', D.connector_N.to_numpy()), ('half', hkey)]:
        for yc, lab in [('drise', 'MIP-Fe55 rise'), ('mip_rise', 'MIP rise alone'), ('mip_rise_cF', 'MIP rise amp-corr(Fe)')]:
            x = demean_by(D.lnR.to_numpy(), key); y = demean_by(D[yc].to_numpy(), key)
            k, ke, _ = odr_fit(x, y, D.lnR_err.to_numpy(), D.drise_err.to_numpy() if yc == 'drise' else D.mip_rise_err.to_numpy())
            r, re_, n = boot_corr(x, y, 500)
            say(f'    {dm:9s} {lab:22s} vs lnR: N {n}  r = {r:.3f} +- {re_:.3f}  ODR slope {k:.1f} +- {ke:.1f} ns')
    say('  CAUTION: lnR and the MIP rise above come from the SAME MIP hits, and within a pad the MIP rise grows with '
        'amplitude, so shared sampling noise correlates them.  Independent test: lnR from one time half of the run, '
        'rise from the other (half-demeaned):')
    for a_, b_ in [('h1', 'h2'), ('h2', 'h1'), ('h1', 'h1')]:
        x = demean_by((np.log(D['mip_med_' + a_]) - D.lnG420).to_numpy(), hkey)
        y = demean_by((D['mip_rise_' + b_] - D.fe_rise).to_numpy(), hkey)
        r, re_, n = boot_corr(x, y, 500)
        okk = np.isfinite(x) & np.isfinite(y)
        say(f'    lnR({a_}) vs MIP-Fe55 rise({b_}): N {n}  r = {r:.3f} +- {re_:.3f}  OLS slope '
            f'{np.polyfit(x[okk], y[okk], 1)[0]:.1f} ns (attenuated by the half-sample lnR noise)'
            + ('   <- same half, biased' if a_ == b_ else ''))
    for rk in ['reg24', 'reg48']:
        rows = []
        for rid, gq in D.groupby(rk):
            x1 = np.nanmedian(np.log(gq.mip_med_h1) - gq.lnG420); y2 = np.nanmean(gq.mip_rise_h2 - gq.fe_rise)
            x2 = np.nanmedian(np.log(gq.mip_med_h2) - gq.lnG420); y1 = np.nanmean(gq.mip_rise_h1 - gq.fe_rise)
            rows.append((x1, y2, x2, y1))
        Q = np.array(rows)
        ra, rae, n = boot_corr(Q[:, 0], Q[:, 1], 2000); rb_, rbe, _ = boot_corr(Q[:, 2], Q[:, 3], 2000)
        oka = np.isfinite(Q[:, 0]) & np.isfinite(Q[:, 1])
        say(f'    {rk} cross-half: corr(lnR h1, drise h2) = {ra:.3f} +- {rae:.3f}; corr(lnR h2, drise h1) = {rb_:.3f} +- '
            f'{rbe:.3f}; OLS slope {np.polyfit(Q[oka, 0], Q[oka, 1], 1)[0]:.1f} ns (N {n})')
    # robustness: drop the pads at the edge of either illumination
    rob = (D.fe_n420 >= 50) & (D.mip_nhit >= 12)
    x = demean_by(D.lnG420.to_numpy(), hkey); y = demean_by(D.lnM.to_numpy(), hkey)
    rA, reA, nA = boot_corr(x[rob], y[rob], 500)
    kA, keA, _ = odr_fit(x[rob], y[rob], D.lnG420_err.to_numpy()[rob], D.lnM_err.to_numpy()[rob])
    x = demean_by(D.lnR.to_numpy(), hkey); y = demean_by(D.drise.to_numpy(), hkey)
    rB, reB, nB = boot_corr(x[rob], y[rob], 500)
    Ib, pIb = moran(np.where(rob, demean_by(D.lnR.to_numpy(), hkey), np.nan), D[['pad_cx', 'pad_cy']].to_numpy())
    say(f'  robustness, pads with >= 50 Fe55 events and >= 12 MIP hits (N {nA}), half-demeaned: lnM vs lnG r = {rA:.3f} +- '
        f'{reA:.3f}, slope {kA:.3f} +- {keA:.3f}; lnR vs drise r = {rB:.3f} +- {reB:.3f}; lnR Moran I {Ib:.3f} (p {pIb:.3g})')
    D.to_csv(f'{out}/pads_fe55_vs_mip.csv', index=False)

    # ---- 5: caveats
    say('\n[5] caveats, quantified')
    say('  per-connector medians (July FEU3 ASICs / Nov FEU6-7 ASICs):')
    cols = ['lnG420', 'slope', 'fe_rise', 'fe_width', 'fe_tot', 'lnM', 'mip_rise', 'mip_width', 'mip_tot', 'lnR', 'drise', 'mip_eff']
    cmed = D.groupby('connector_N')[cols].median()
    say('    ' + cmed.round(4).to_string().replace('\n', '\n    '))
    cmed.to_csv(f'{out}/connector_medians.csv')
    # block-of-64 pattern: one-way ANOVA on halves (within connector) for each quantity
    say('  64-channel block (DREAM ASIC = connector half) effect, Kruskal-Wallis over the 8 halves '
        '(H, p) and the same after removing connector medians:')
    hkey = D.connector_N * 2 + (D.half == 'top')
    for c in ['lnG420', 'fe_rise', 'lnM', 'mip_rise', 'lnR', 'drise']:
        grps = [gq[c].dropna().to_numpy() for _, gq in D.groupby(hkey)]
        H1, p1 = stats.kruskal(*grps)
        dd = pd.Series(demean_by(D[c].to_numpy(), D.connector_N))
        H2, p2 = stats.kruskal(*[dd[hkey == k].dropna().to_numpy() for k in sorted(hkey.unique())])
        say(f'    {c:9s}: H = {H1:6.1f} (p {p1:.2g}); within connectors H = {H2:6.1f} (p {p2:.2g})')
    # HV difference: per-pad slope spread x 10 V
    say(f'  HV 420 -> 430 V: sd over pads of 10 V x slope = {10 * D.slope.std():.4f} (in lnG), against sd(lnG420) '
        f'{D.lnG420.std():.3f}; part of it is noise (expected from errors {10 * np.sqrt(np.nanmean(D.slope_err ** 2)):.4f}). '
        f'Drift field 700 V/cm (July, 280 V / 4 mm) vs 750 V/cm (Nov).')
    # illumination edge
    rr_ = stats.spearmanr(D.fe_n420, D.lnR, nan_policy='omit').correlation
    say(f'  Fe55 illumination: events/pad at 420 V min {D.fe_n420.min():.0f}, median {D.fe_n420.median():.0f}, max {D.fe_n420.max():.0f}; '
        f'Spearman(fe_n420, lnR) = {rr_:.3f}')
    # saturation
    s420 = E[(E.hv == 420)].groupby('channel_id').sat.mean()
    say(f'  saturation: Fe55 420 V lead-pad saturated fraction per pad median {s420.median():.4f} (max {s420.max():.3f}); '
        f'MIP own-hit saturated fraction per pad median {D.mip_sat.median():.3f} (max {D.mip_sat.max():.3f}); '
        f'the medians are insensitive to the extrapolated saturated values (rank statistics).')
    # threshold censoring in MIP
    rs_ = stats.spearmanr(D.lnG420, D.mip_eff, nan_policy='omit').correlation
    say(f'  MIP censoring: own-pad hit fraction per pad median {D.mip_eff.median():.3f}; Spearman(lnG420, eff) = {rs_:.3f}; '
        f'own-hit amp x cos 5th percentile {np.nanpercentile(Sh.ac, 5):.0f} ADC vs 5-sigma thresholds of a few tens of ADC.')
    hdr = [
        'CUTS AND ASSUMPTIONS',
        f'  chamber identity: FEU 3 / CAEN 8:0-8:1 in the 18 Jul Fe55 run = det2 (run_config "P2_OUT ... (det2)"; SPS run '
        'configs: P2_IN = det2 on FEU 3 / 8:0-8:1). Station names disagree (chamber_history.py), the chamber does not.',
        '  gas Ar/iC4H10 95/5 in both runs; 60 ns sampling, 32 samples in both; DREAM peaking-time register of the July '
        'self-trigger cfg (P2SelfTrigger.cfg) is not archived locally -> not verified equal to CosmicTb_P2.cfg.',
        '  analyzer: July hit files are from the pre-2026-07-24 WaveformAnalyzer (no trunc/significance branches): '
        'baseline = min of pre-samples, integral between threshold crossings, TOT = integer samples on the raw waveform. '
        'November: matched-filter gate, median baseline, full-span integral, interpolated TOT. Only rise = time_of_max - '
        'time has the same definition.',
        f'  Fe55: HV-settle cut and HV-spark veto as stage 18; no burst veto; per event lead pad + pads with centre within '
        f'{CLUS_R:g} mm (cluster q); q > {FE_FLOOR:.0f} ADC; per-pad photopeak = iterated median in +-{WIN:.0%}; '
        'ln(peak) vs mesh V fitted per pad over 385-420 V (points with peak > 320 ADC, >= 4 points). Shape: lead-pad hit of '
        '420 V photopeak events, unsaturated, 0 < rise < 600 ns, 10 % trimmed mean per pad (TOT: mean).',
        f'  MIP: det2_long_repos stage-34 trimmed alignment (scale {float(mst.scale):.4f}), plane z {cfg.det_plane_z():g} mm, '
        f'rays inside a live pad tile by > {args.edge:g} mm and sd > 2 mm, spark/burst-vetoed events removed, t < '
        f'{cfg.T_MAX_H} h; own-pad hit amp x cos(theta); MPV = median (primary), low-70 % mean, cluster median; shape: own '
        'hit unsaturated, untruncated, amp > 100 ADC, 0 < rise < 600 ns.',
        '  pillar lattice: not used here (3.0 mm measured lattice, rule 1); pillar shadows only reduce the MIP charge of a '
        'fraction of tracks and enter the MPV as a common low tail; the > edge cut keeps tracks off pad borders.',
        '  common pads: connectors 4-7 minus the November drops (FEU6 even 384-446 dead pin row, FEU7 254).',
        '']
    with open(f'{out}/summary.txt', 'w') as f:
        f.write('\n'.join(hdr + L) + '\n')

    # ---- figures
    vid = dict(zip(g['pad_id'], g['tile_verts']))
    vid = {k: v for k, v in vid.items() if k in set(D.channel_id)}
    ids = D.channel_id.to_numpy()
    fig, axs = plt.subplots(2, 3, figsize=(17, 10))
    draw_map(axs[0, 0], vid, ids, D.lnG420.to_numpy(), 'Fe55 ln(photopeak) at 420 V (July, per-pad fit)', label='ln ADC')
    draw_map(axs[0, 1], vid, ids, D.lnM.to_numpy(), f'MIP ln(median amp cos theta) (Nov, >{args.edge:g} mm inside)', label='ln ADC')
    draw_map(axs[0, 2], vid, ids, demean_by(D.lnR.to_numpy(), D.connector_N), 'ln(MIP/Fe55), connector medians removed',
             cmap='coolwarm', label='ln ratio', vlim=(-0.4, 0.4))
    draw_map(axs[1, 0], vid, ids, D.slope.to_numpy() * 100, 'Fe55 d lnG/dV per pad', label='% per V')
    draw_map(axs[1, 1], vid, ids, D.fe_n420.to_numpy(), 'Fe55 events per pad at 420 V (illumination)', cmap='inferno', label='events')
    draw_map(axs[1, 2], vid, ids, D.mip_nhit.to_numpy(), 'MIP own-pad hits per pad', cmap='inferno', label='hits')
    fig.suptitle('det2 connectors 4-7: Fe55 gain (18 Jul 2026) vs cosmic MIP (25 Nov 2026)')
    fig.tight_layout(); fig.savefig(f'{out}/maps_gain_mip_ratio.png', dpi=140); plt.close(fig)

    fig, axs = plt.subplots(2, 3, figsize=(17, 10))
    draw_map(axs[0, 0], vid, ids, D.fe_rise.to_numpy(), 'Fe55 rise (time_of_max - time), 420 V photopeak', label='ns')
    draw_map(axs[0, 1], vid, ids, D.mip_rise.to_numpy(), 'MIP rise (own-pad hit)', label='ns')
    draw_map(axs[0, 2], vid, ids, demean_by(D.drise.to_numpy(), D.connector_N), 'MIP - Fe55 rise, connector medians removed',
             cmap='coolwarm', label='ns', vlim=(-10, 10))
    draw_map(axs[1, 0], vid, ids, D.fe_width.to_numpy(), 'Fe55 60 ns x integral/amp (July definition)', label='ns')
    draw_map(axs[1, 1], vid, ids, D.mip_width.to_numpy(), 'MIP 60 ns x integral/amp (Nov definition)', label='ns')
    draw_map(axs[1, 2], vid, ids, D.fe_rise_err.to_numpy(), 'Fe55 rise per-pad error', cmap='inferno', label='ns')
    fig.suptitle('det2 pulse-shape proxies, connectors 4-7')
    fig.tight_layout(); fig.savefig(f'{out}/maps_shape.png', dpi=140); plt.close(fig)

    fig, axs = plt.subplots(1, 3, figsize=(17, 5.2))
    ax = axs[0]
    for c in sorted(D.connector_N.unique()):
        q = D[D.connector_N == c]
        ax.errorbar(q.lnG420, q.lnM, xerr=q.lnG420_err, yerr=q.lnM_err, fmt='o', ms=2.5, alpha=0.5, lw=0.5, label=f'conn {c}')
    ax.set_xlabel('ln Fe55 photopeak at 420 V [ln ADC]'); ax.set_ylabel('ln MIP median amp cos theta [ln ADC]')
    r_, re_, k_, ke_ = res4a[('median vs lnG420', 'none')]
    ax.set_title(f'per pad: r = {r_:.2f} +- {re_:.2f}, ODR slope {k_:.2f} +- {ke_:.2f}'); ax.legend(fontsize=7)
    ax = axs[1]
    for rk, mk in [('reg24', 's'), ('reg48', 'o')]:
        R = RG[rk]
        ax.errorbar(R.lnG, R.lnM, xerr=R.lnG_err, yerr=R.lnM_err, fmt=mk, ms=5, alpha=0.8, label=rk)
    ax.set_xlabel('ln Fe55 gain (region mean)'); ax.set_ylabel('ln MIP MPV (region mean)')
    ax.set_title('region level'); ax.legend()
    ax = axs[2]
    x = demean_by(D.lnG420.to_numpy(), D.connector_N)
    ax.errorbar(x, demean_by(D.fe_rise.to_numpy(), D.connector_N), yerr=D.fe_rise_err, fmt='o', ms=2.5, alpha=0.5,
                lw=0.5, label='Fe55 rise')
    ax.errorbar(x, demean_by(D.mip_rise.to_numpy(), D.connector_N), yerr=D.mip_rise_err, fmt='o', ms=2.5, alpha=0.4,
                lw=0.5, label='MIP rise')
    ax.set_xlabel('ln Fe55 gain, connector-demeaned'); ax.set_ylabel('rise, connector-demeaned [ns]')
    ax.set_title('rise vs Fe55 gain within connectors'); ax.legend()
    fig.tight_layout(); fig.savefig(f'{out}/scatter_mip_vs_fe55.png', dpi=140); plt.close(fig)

    fig, axs = plt.subplots(1, 3, figsize=(17, 5))
    ax = axs[0]
    for hv in [420, 405, 390]:
        q = E[(E.hv == hv) & (E.q > 0)]
        ax.hist(q.q, bins=np.linspace(0, 2500, 126), histtype='step', label=f'{hv} V cluster q')
        ax.hist(q.a, bins=np.linspace(0, 2500, 126), histtype='step', ls='--', label=f'{hv} V lead pad')
    ax.set_yscale('log'); ax.set_xlabel('ADC'); ax.set_ylabel('events / 20 ADC'); ax.legend(fontsize=7)
    ax.set_title('Fe55 spectra (all pads)')
    ax = axs[1]
    ax.plot(SH.hv, SH.f_nb, 'o-', label='f(>=1 neighbour fired), Fe55 photopeak')
    ax.plot(SH.hv, SH.r_med, 's-', label='median a_nb / a_lead')
    ax.axhline((S[S.hit].n_nb > 0).mean(), color='crimson', ls='--', label='MIP Nov: f(>=1 neighbour)')
    ax.set_xlabel('mesh HV [V] (July)'); ax.set_ylabel('fraction / ratio'); ax.set_ylim(0, 1); ax.legend(fontsize=7)
    ax.set_title('second-pad charge: Fe55 July vs MIP November')
    ax = axs[2]
    ax.errorbar(D.slope * 100, D.lnR - np.nanmedian(D.lnR), xerr=D.slope_err * 100, yerr=D.lnR_err, fmt='o', ms=2, alpha=0.4, lw=0.4)
    ax.set_xlabel('Fe55 d lnG/dV [% per V]'); ax.set_ylabel('ln(MIP/Fe55) - median')
    ax.set_title('ratio vs per-pad gain slope')
    fig.tight_layout(); fig.savefig(f'{out}/fe55_spectra_sharing.png', dpi=140); plt.close(fig)

    # amplitude dependence figure
    fig, axs = plt.subplots(1, 3, figsize=(17, 5))
    for ax, c, lab in zip(axs, ['rise', 'width', 'totn'], ['rise [ns]', '60 ns x integral/amp [ns]', 'TOT [ns]']):
        bins = np.geomspace(150, 4000, 16)
        cen = np.sqrt(bins[1:] * bins[:-1])
        mm = Ms.groupby(pd.cut(Ms.a, bins), observed=False)[c].median().to_numpy()
        fm = ee.groupby(pd.cut(ee.a, bins), observed=False)[c].median().to_numpy()
        ax.plot(cen, mm, 'o-', label='MIP Nov (own hit)'); ax.plot(cen, fm, 's-', label='Fe55 July (photopeak, 395-420 V)')
        ax.set_xscale('log'); ax.set_xlabel('hit amplitude [ADC]'); ax.set_ylabel(lab); ax.legend(fontsize=8)
        ax.set_title(f'{lab} vs amplitude')
    fig.tight_layout(); fig.savefig(f'{out}/shape_vs_amplitude.png', dpi=140); plt.close(fig)
    print(f'products in {out}')


if __name__ == '__main__':
    main()
