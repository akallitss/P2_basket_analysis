#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
23_drift_gap_bow.py

Is the structure in the stage-19 timing surface a MECHANICAL bow of the
mesh/readout assembly, or front-end electronics?

Stage 19 draws the per-pad peak-time map and shows (split-half) that its
structure is real, but it cannot say what causes it. This stage decomposes it
and tests the mechanical hypothesis three independent ways:

  1. MODEL       per-pad offset = per-connector constant (cable / front-end
                 delay, 8 discrete dof) + one smooth quadratic surface in
                 (pad_cx, pad_cy). If the residual falls to the split-half
                 statistical error, nothing else is present -- no isolated bad
                 channels, no unmodelled structure.
  2. AMPLITUDE   the same model fitted to per-pad median pulse height. A longer
                 drift path means more primary ionisation AND more drift time,
                 so a real gap variation moves both together; a cable delay
                 moves only the time. (Time-walk would make big pulses EARLIER,
                 so it cannot fake a positive correlation.)
  3. 1/v_d       across a drift scan the drift velocity changes by ~10x. A
                 geometric gap variation gives a timing swing proportional to
                 1/v_d, so swing * v_d is CONSTANT and equals the gap variation
                 in um. An electronic offset is constant in ns instead.

The sign convention is stage 19's: offset = median time_of_max - detector
median, so LATER (red) = longer drift = LARGER drift gap at that pad. A dome
(quadratic maximum in the middle) therefore means the gap is largest in the
middle, i.e. the mesh/readout assembly bows AWAY from the drift electrode --
whereas a drift foil sagging toward the mesh under gravity would give a bowl.

The gap variation in mm uses cfg.DRIFT_GAP_MM and the drift velocity; the
time_of_max-to-gap factor is not exactly 1, so timing and amplitude give a
RANGE, not a number. Quote the range.

Memory: the per-pad amplitude pass streams chunk by chunk through
p2_io.iter_hits and keeps only a per-pad histogram (n_pads x n_bins ints), so
nothing scales with the hit count -- same bound as every other reduction here.

Products (<Analysis>/<detN>/<run>/<sub_run>/23_drift_gap_bow/):
  bow_decomposition<sfx>.png   observed / model / residual maps + variance table
  bow_timing_vs_amplitude<sfx>.png  the two shape terms, maps and correlation
  bow_profile<sfx>.png         offset vs distance from the fitted apex
  bow_surface_3d<sfx>.png      the sag in 3D: flat cathode + bowed assembly
  bow_raw_binned<sfx>.png      the bow with NO surface model: raw pads, hard
                               spatial bins, and a sliding (running) average
  bow_residual_structure<sfx>.png  sliding map of what the model does NOT
                               explain, with a label-shuffle significance null
  bow_bin_scan<sfx>.png        precision vs bin size AND vs sliding radius;
                               why anything below the pad pitch buys nothing
  bow_vd_scaling<sfx>.png      swing and swing*v_d vs drift HV (needs --scan-key)
  bow_fit<sfx>.csv             fit parameters, apex, implied gap variation

Usage:
  python3 23_drift_gap_bow.py det5_long1 [--scan-key det5_driftscan1]
        [--sig-amp 300] [--min-hits 20] [--no-veto-sparks]
"""

import argparse
import glob
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

import p2_qa_config as qa
qa.setup_paths()
import p2_io as p2io
import p2_mapping as pmap
import p2_sparks as ps

# Per-pad amplitude histogram: 0..5000 ADC. P2 pad pulses sit well inside this
# (median ~500 ADC on det5); anything above is clipped into the top bin, which
# only affects the tail, never the median we take from it.
AMP_MAX = 5000.0
AMP_BIN = 25.0


def channel_table(cfg, strategy='reverse'):
    return pmap.build_channel_table(
        cfg.run_config_path, cfg.MAP_CSV_PATH, det_type=cfg.DET_TYPE,
        det_name=cfg.DET_NAME, strategy=strategy,
        drop_connectors=cfg.DEAD_CONNECTORS,
        strategy_overrides=cfg.STRATEGY_OVERRIDES)


def pad_median_amplitude(hits_dir, ct, sig_amp, exclude_events=None,
                         min_amp=0.0, drop_pads=()):
    """Per-pad median amplitude over signal hits, streamed.

    Returns a DataFrame(channel_id, amp_median, n_amp). Only hits above
    sig_amp count, matching the band stage 19 uses for the peak time, so the
    two per-pad quantities are built from the same sample -- which is also why
    drop_pads must be the SAME hot/noisy set stage 19 removed, or the amplitude
    and the timing would be describing different detectors.
    """
    m = ct[ct['mapped']]
    if len(drop_pads):
        m = m[~m['channel_id'].isin(set(drop_pads))]
    lut = pd.Series(m['channel_id'].to_numpy(),
                    index=pd.MultiIndex.from_arrays(
                        [m['feu'].astype(int), m['channel'].astype(int)]))
    pads = np.sort(m['channel_id'].unique())
    row_of = {int(p): i for i, p in enumerate(pads)}
    nb = int(round(AMP_MAX / AMP_BIN))
    hist = np.zeros((len(pads), nb), dtype=np.int64)
    excl = set(exclude_events) if exclude_events is not None else None

    for h in p2io.iter_hits(hits_dir, ['eventId', 'channel', 'amplitude', 'feu'],
                            feus=ct.attrs['feus'], min_amp=min_amp):
        if excl:
            h = h[~h['eventId'].isin(excl)]
        h = h[h['amplitude'] > sig_amp]
        if not len(h):
            continue
        cid = lut.reindex(pd.MultiIndex.from_arrays(
            [h['feu'].astype(int), h['channel'].astype(int)])).to_numpy()
        ok = pd.notna(cid)
        if not ok.any():
            continue
        rows = np.array([row_of[int(c)] for c in cid[ok]])
        cols = np.clip((h['amplitude'].to_numpy()[ok] / AMP_BIN).astype(int),
                       0, nb - 1)
        np.add.at(hist, (rows, cols), 1)

    n = hist.sum(axis=1)
    centres = (np.arange(nb) + 0.5) * AMP_BIN
    med = np.full(len(pads), np.nan)
    nz = n > 0
    cum = np.cumsum(hist[nz], axis=1)
    half = (n[nz] / 2.0)[:, None]
    med[nz] = centres[(cum >= half).argmax(axis=1)]
    return pd.DataFrame({'channel_id': pads, 'amp_median': med, 'n_amp': n})


def fit_bow(d, value_col):
    """Least-squares fit of `value_col` = per-connector constant + quadratic
    surface in (pad_cx, pad_cy).

    Returns (shape, resid, info): `shape` is the quadratic part with its median
    removed (the connector constants absorb the overall level, so only the
    SHAPE is meaningful), `resid` the leftover per pad, `info` a dict with the
    stationary point, the Hessian eigenvalues (both negative = dome, both
    positive = bowl) and the robust 2.5-97.5 % peak-to-peak of the shape.
    """
    x = d['pad_cx'].to_numpy(float)
    y = d['pad_cy'].to_numpy(float)
    z = d[value_col].to_numpy(float)
    dummies = pd.get_dummies(d['connector_N']).to_numpy(float)
    nq = dummies.shape[1]
    A = np.column_stack([dummies, x, y, x * x, y * y, x * y])
    c, *_ = np.linalg.lstsq(A, z, rcond=None)
    shape = A[:, nq:] @ c[nq:]
    shape = shape - np.median(shape)
    resid = z - A @ c
    a, b, axx, ayy, axy = c[nq:]
    H = np.array([[2 * axx, axy], [axy, 2 * ayy]])
    try:
        x0, y0 = np.linalg.solve(H, [-a, -b])
    except np.linalg.LinAlgError:
        x0 = y0 = np.nan
    ev = np.linalg.eigvalsh(H)
    kind = ('dome' if ev.max() < 0 else 'bowl' if ev.min() > 0 else 'saddle')
    return shape, resid, dict(
        apex_x=float(x0), apex_y=float(y0), kind=kind,
        eig_lo=float(ev[0]), eig_hi=float(ev[1]),
        anisotropy=float(abs(ev[0] / ev[1])) if ev[1] else np.nan,
        pp=float(np.percentile(shape, 97.5) - np.percentile(shape, 2.5)),
        conn_const={int(k): float(v) for k, v in
                    zip(sorted(d['connector_N'].unique()), c[:nq])})


def bootstrap_fit(d, value_col, n=300, seed=0):
    """Pad-resampling bootstrap of the surface fit.

    Reports what is actually constrained. The CURVATURE SIGN (dome vs bowl) and
    the peak-to-peak over the instrumented area are well determined; the apex
    POSITION is not -- the active area is a fan-shaped triangle and the surface
    is strongly elliptical, so its stationary point is an extrapolation and can
    land far outside the pads. Do not read the apex as "the centre of the bow".
    """
    rng = np.random.default_rng(seed)
    xs, ys, pps, kinds = [], [], [], []
    for _ in range(n):
        s = d.iloc[rng.choice(len(d), len(d), replace=True)]
        if s['connector_N'].nunique() < d['connector_N'].nunique():
            continue
        try:
            _, _, info = fit_bow(s, value_col)
        except np.linalg.LinAlgError:
            continue
        kinds.append(info['kind'])
        pps.append(info['pp'])
        if np.isfinite(info['apex_x']):
            xs.append(info['apex_x'])
            ys.append(info['apex_y'])
    if len(kinds) < 20:
        return dict(apex_sx=np.nan, apex_sy=np.nan, pp_sd=np.nan,
                    dome_frac=np.nan, n=len(kinds))
    return dict(apex_sx=float(np.std(xs)) if xs else np.nan,
                apex_sy=float(np.std(ys)) if ys else np.nan,
                pp_sd=float(np.std(pps)),
                dome_frac=kinds.count('dome') / len(kinds), n=len(kinds))


def raw_profile(d, ucol, nbin=12):
    """Median RAW per-pad timing offset and amplitude in bins of the dome
    coordinate, with the per-connector constants removed first.

    This is the honest cross-check: correlating the two FITTED surfaces is
    close to circular (both are smooth quadratics over the same pads), so the
    comparison must be made on the measurements themselves.
    """
    g = d.copy()
    for col, out in (('off', 'off_c'), ('amp_median', 'amp_c')):
        g[out] = g[col] - g.groupby('connector_N')[col].transform('mean')
    g['q'] = pd.qcut(g[ucol], nbin, labels=False)
    return g.groupby('q').agg(u=(ucol, 'median'), t=('off_c', 'median'),
                              a=('amp_c', 'median'), n=('off_c', 'size'))


def _pad_map(ax, d, v, title, label, cmap='coolwarm', sym=True):
    lim = float(np.nanpercentile(np.abs(v), 97.5)) or 1.0
    kw = dict(vmin=-lim, vmax=lim) if sym else {}
    s = ax.scatter(d['pad_cx'], d['pad_cy'], c=v, s=14, marker='s',
                   cmap=cmap, **kw)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel('pad_cx [mm]')
    ax.set_ylabel('pad_cy [mm]')
    ax.set_aspect('equal')
    plt.colorbar(s, ax=ax, label=label)


def plot_decomposition(d, off, shape, resid, info, stat, cfg, out, sfx, args,
                       gen, frac):
    model = off - resid
    # Report how much of the REAL structure the smooth surface accounts for.
    # "Residual == noise" is the wrong bar: with enough hits it never quite is.
    explained = frac >= 0.6
    verdict = (f'{gen:.1f} ns of real structure left, but the surface '
               f'already accounts for {100*frac:.0f} %'
               if explained else
               f'{gen:.1f} ns left — the surface explains only '
               f'{100*frac:.0f} %: NOT a smooth bow')
    fig, axs = plt.subplots(1, 4, figsize=(21, 5))
    _pad_map(axs[0], d, off, 'Observed per-pad offset\n(stage 19: red = LATER '
             '= larger drift gap)', 'offset [ns]')
    _pad_map(axs[1], d, model, 'Model: 8 connector constants\n+ one quadratic '
             'surface', 'model [ns]')
    _pad_map(axs[2], d, resid, f'Residual — std {resid.std():.2f} ns vs '
             f'{stat:.2f} ns statistical\n({verdict})', 'residual [ns]')
    ax = axs[3]
    ax.hist(resid, bins=40, color='0.7', edgecolor='k', lw=0.4,
            density=True, label=f'residual (std {resid.std():.2f} ns)')
    xs = np.linspace(resid.min(), resid.max(), 200)
    ax.plot(xs, np.exp(-0.5 * (xs / stat) ** 2) / (stat * np.sqrt(2 * np.pi)),
            'r-', lw=2, label=f'split-half statistical\nerror {stat:.2f} ns')
    ax.set_xlabel('residual [ns]')
    ax.set_ylabel('density')
    ax.legend(fontsize=8)
    ax.set_title(f'Smooth surface accounts for {100*frac:.0f} % of the\n'
                 'genuine structure' if explained else
                 f'Only {100*frac:.0f} % explained: real non-quadratic\n'
                 'structure dominates the residual', fontsize=10)
    ax.grid(alpha=0.3)
    fig.suptitle(
        f'{cfg.DET_NAME} drift-gap bow — decomposition of the timing surface — '
        f'{cfg.DET_TAG} {cfg.RUN}\n'
        f'quadratic part is a {info["kind"].upper()} (anisotropy '
        f'{info["anisotropy"]:.1f})' +
        (' — a dome means the gap is LARGEST in the middle'
         if info['kind'] == 'dome' and explained else
         ' — NOT a clean bow: do not read a gap variation off this'))
    fig.tight_layout()
    fig.savefig(os.path.join(out, f'bow_decomposition{sfx}.png'), dpi=200,
                bbox_inches='tight')
    plt.close(fig)


def plot_timing_vs_amplitude(d, st, sa, it, ia, amed, prof, r_raw,
                             bt, ba, cfg, out, sfx, args):
    fig, axs = plt.subplots(1, 3, figsize=(17.5, 5))
    _pad_map(axs[0], d, st, f'TIMING — FITTED quadratic only (smooth by\n'
             f'construction; stage 19 has the raw map). {it["kind"]}, '
             f'{it["pp"]:.1f} ± {bt["pp_sd"]:.1f} ns',
             'ns  (red = later)')
    _pad_map(axs[1], d, sa, f'AMPLITUDE — FITTED quadratic only. {ia["kind"]}, '
             f'\n{ia["pp"]:.0f} ADC = {100*ia["pp"]/amed:.0f} % of the median',
             'ADC  (red = bigger)', cmap='RdBu_r')
    ax = axs[2]
    sc = ax.scatter(prof['t'], prof['a'], c=prof['u'], s=80, cmap='viridis',
                    zorder=3, edgecolor='k', lw=0.5)
    ax.set_xlabel('median peak-time offset [ns]  (later →)')
    ax.set_ylabel('median amplitude offset [ADC]  (bigger →)')
    ax.grid(alpha=0.3)
    plt.colorbar(sc, ax=ax, label='distance from apex [mm]')
    ax.set_title(f'RAW measurements, binned by distance from the apex\n'
                 f'(connector constants removed): r = {r_raw:+.3f}', fontsize=10)
    agree = (it['kind'] == 'dome' and ia['kind'] == 'dome' and r_raw > 0.5)
    verdict = (
        'both are DOMES and the raw binned measurements track each other '
        f'(r = {r_raw:+.3f}) — the drift-gap signature'
        if agree else
        f'timing is a {it["kind"]}, amplitude a {ia["kind"]}, and the raw '
        f'measurements are correlated at r = {r_raw:+.3f} — '
        + ('ANTI-correlated, i.e. NOT a gap variation (bigger pulses peaking '
           'earlier is time-walk on a GAIN gradient)' if r_raw < 0 else
           'they do not track each other'))
    fig.suptitle(
        f'{cfg.DET_NAME} — timing vs pulse height — {cfg.DET_TAG} {cfg.RUN}\n'
        + verdict, y=1.0)
    fig.tight_layout()
    fig.savefig(os.path.join(out, f'bow_timing_vs_amplitude{sfx}.png'), dpi=200,
                bbox_inches='tight')
    plt.close(fig)


def _fan_triangulation(x, y, max_edge_factor=2.5):
    """Delaunay triangulation of the pad centres with the long triangles
    removed, so the rendered sheet follows the fan outline instead of bridging
    across the empty regions of the board."""
    import matplotlib.tri as mtri
    tri = mtri.Triangulation(x, y)
    t = tri.triangles
    edge = np.maximum.reduce([
        np.hypot(x[t[:, i]] - x[t[:, j]], y[t[:, i]] - y[t[:, j]])
        for i, j in ((0, 1), (1, 2), (2, 0))])
    tri.set_mask(edge > max_edge_factor * np.median(edge))
    return tri


def plot_3d_bow(d, shape, info, cfg, out, sfx, vd, gap_mm):
    """The sag, drawn as it physically is: a flat drift cathode above, and the
    bowed mesh/readout assembly below it, separated by the local drift gap.

    z is the measured gap converted to mm (shape_ns * v_d), so where the peak
    time is LATER the assembly sits FURTHER from the cathode. The z axis is in
    mm while x/y span ~500 mm, so the vertical is exaggerated by a large factor
    -- it is stated on the figure, because at true scale a ~1 mm bow across a
    500 mm board is invisible.
    """
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers '3d')
    x = d['pad_cx'].to_numpy(float)
    y = d['pad_cy'].to_numpy(float)
    dz = shape * vd / 1000.0                       # mm of extra gap
    gap = gap_mm + dz                              # local drift gap [mm]
    tri = _fan_triangulation(x, y)
    exag = (max(x.max() - x.min(), y.max() - y.min())
            / max(gap.max() - gap.min(), 1e-9))

    fig = plt.figure(figsize=(16, 7.4))

    # -- left: the real thing, cathode + bowed assembly ------------------- #
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    pad = 0.03 * (x.max() - x.min())
    xg = np.array([x.min() - pad, x.max() + pad, x.max() + pad, x.min() - pad])
    yg = np.array([y.min() - pad, y.min() - pad, y.max() + pad, y.max() + pad])
    ax.plot_trisurf(xg, yg, np.full(4, gap_mm), color='steelblue', alpha=0.20,
                    linewidth=0, shade=False)
    ax.text(x.min(), y.max(), gap_mm * 1.05, 'drift cathode (assumed flat)',
            color='steelblue', fontsize=9)
    # The cathode sits at z = gap_mm, so the assembly is at z = -dz: where the
    # peak time is LATER (dz > 0) the surface drops away from the cathode.
    # plot_trisurf colours per TRIANGLE from z, so colour by z with a reversed
    # map (= colouring by the gap, which is gap_mm - z) rather than set_array,
    # whose per-vertex length would not match the triangle count.
    ax.plot_trisurf(tri, -dz, cmap='coolwarm_r', linewidth=0.1,
                    antialiased=True)
    sm = plt.cm.ScalarMappable(
        cmap='coolwarm', norm=plt.Normalize(gap.min(), gap.max()))
    sm.set_array([])
    ax.set_xlabel('pad_cx [mm]')
    ax.set_ylabel('pad_cy [mm]')
    ax.set_zlabel('height [mm]')
    ax.set_zlim(min(-dz.max() * 1.6, -0.2), gap_mm * 1.15)
    ax.set_box_aspect((1, 1, 0.62))
    plt.colorbar(sm, ax=ax, shrink=0.6, pad=0.10, label='local drift gap [mm]')
    ax.set_title(f'Drift gap: flat cathode above, bowed mesh/readout below\n'
                 f'nominal {gap_mm:g} mm; vertical exaggeration ≈ {exag:.0f}×',
                 fontsize=10)

    # -- right: the deflection alone, with the per-pad measurements ------- #
    ax2 = fig.add_subplot(1, 2, 2, projection='3d')
    m2 = ax2.plot_trisurf(tri, dz, cmap='coolwarm', linewidth=0.1,
                          antialiased=True, alpha=0.85)
    pad_dz = (d['off'].to_numpy()
              - d.groupby('connector_N')['off'].transform('mean').to_numpy()
              ) * vd / 1000.0
    ax2.scatter(x, y, pad_dz, s=4, color='0.25', alpha=0.35, depthshade=False)
    ax2.set_xlabel('pad_cx [mm]')
    ax2.set_ylabel('pad_cy [mm]')
    ax2.set_zlabel('gap deviation [mm]')
    ax2.set_box_aspect((1, 1, 0.45))
    plt.colorbar(m2, ax=ax2, shrink=0.6, pad=0.10,
                 label='gap deviation from the mean [mm]')
    ax2.set_title('Fitted deflection surface, with the per-pad\nmeasurements '
                  '(grey, connector constants removed)', fontsize=10)

    # same robust measure the CSV and the other figures quote
    pp_mm = float(np.percentile(dz, 97.5) - np.percentile(dz, 2.5))
    fig.suptitle(
        f'{cfg.DET_NAME} drift-gap bow in 3D — {cfg.DET_TAG} {cfg.RUN}\n'
        f'peak-to-peak {pp_mm:.2f} mm on a {gap_mm:g} mm gap '
        f'({100*pp_mm/gap_mm:.0f} %); v_d = {vd:.1f} µm/ns — the DEPTH is '
        'measured, the SHAPE is the fitted quadratic', y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(out, f'bow_surface_3d{sfx}.png'), dpi=200,
                bbox_inches='tight')
    plt.close(fig)


def bin_size_scan(d, stat, sizes=(10, 12, 15, 20, 25, 30, 40, 55, 70, 90, 120)):
    """Precision vs bin size. Averaging only starts once a bin holds more than
    one pad, so nothing below the pad pitch (~11.8 mm here) buys anything: the
    bin error stays at the single-pad error and the 'span' just grows with the
    noise. Returns a DataFrame for the figure and the printout."""
    rows = []
    for b in sizes:
        ax = (np.floor(d['pad_cx'] / b) + 0.5) * b
        ay = (np.floor(d['pad_cy'] / b) + 0.5) * b
        t = d.groupby([ax, ay]).agg(m=('off_c', 'median'),
                                    n=('off_c', 'size')).reset_index()
        err = stat / np.sqrt(t['n'])
        rows.append(dict(bin_mm=b, n_bins=len(t),
                         pads_per_bin=float(t['n'].median()),
                         err_ns=float(err.median()),
                         span_ns=float(t['m'].max() - t['m'].min()),
                         max_sig=float(np.abs(t['m'] / err).max())))
    return pd.DataFrame(rows)


def radius_scan(d, stat, radii=(5, 10, 11, 12, 13, 15, 18, 20, 25, 30, 40,
                                55, 70)):
    """Pads per sliding window vs window radius.

    On a regular pad lattice this is a STAIRCASE, not a smooth curve: the window
    picks up whole shells of neighbours at a time, so it holds 1 pad until the
    radius crosses the pitch, then jumps. Anything below the pitch averages
    nothing at all -- exactly like a sub-pitch hard bin.
    """
    xy = d[['pad_cx', 'pad_cy']].to_numpy()
    tree = cKDTree(xy)
    vals = d['off_c'].to_numpy()
    rows = []
    for r in radii:
        idx = tree.query_ball_point(xy, r=r)
        n = np.array([len(i) for i in idx])
        m = np.array([np.median(vals[i]) for i in idx])
        err = stat / np.sqrt(n)
        rows.append(dict(radius_mm=r, pads_per_window=float(np.median(n)),
                         err_ns=float(np.median(err)),
                         span_ns=float(m.max() - m.min()),
                         max_sig=float(np.abs(m / err).max())))
    return pd.DataFrame(rows)


def plot_bin_scan(scan, rscan, pitch, cfg, out, sfx):
    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(19, 5))
    a1.plot(scan['bin_mm'], scan['pads_per_bin'], 'o-', color='steelblue',
            lw=2, ms=7, label='pads per bin (median)')
    a1.plot(scan['bin_mm'], scan['err_ns'], 's--', color='crimson', lw=2, ms=6,
            label='error per bin [ns]')
    a1.axvline(pitch, color='k', ls=':', lw=1.5)
    a1.text(pitch * 1.05, a1.get_ylim()[1] * 0.6,
            f'pad pitch\n{pitch:.1f} mm', fontsize=8)
    a1.set_xscale('log'); a1.set_yscale('log')
    a1.set_xlabel('bin size [mm]'); a1.grid(alpha=0.3, which='both')
    a1.legend(fontsize=9)
    a1.set_title('Below the pad pitch a bin holds ONE pad:\nno averaging, error '
                 'stays at the single-pad value', fontsize=10)
    a2.plot(scan['bin_mm'], scan['max_sig'], 'o-', color='seagreen', lw=2, ms=7)
    a2.axvline(pitch, color='k', ls=':', lw=1.5)
    a2.set_xscale('log')
    a2.set_xlabel('bin size [mm]')
    a2.set_ylabel('largest |median| / error  [σ]')
    a2.grid(alpha=0.3, which='both')
    a2.set_title('Significance of the spatial structure\n(bigger bins = more '
                 'pads averaged = clearer)', fontsize=10)
    a3.step(rscan['radius_mm'], rscan['pads_per_window'], where='post', lw=2,
            color='darkorange', label='pads per window (median)')
    a3.plot(rscan['radius_mm'], rscan['err_ns'], 's--', color='crimson', lw=2,
            ms=6, label='error per window [ns]')
    a3.axvline(pitch, color='k', ls=':', lw=1.5)
    a3.text(pitch * 1.05, a3.get_ylim()[1] * 0.5,
            f'pad pitch\n{pitch:.1f} mm', fontsize=8)
    a3.set_xscale('log'); a3.set_yscale('log')
    a3.set_xlabel('sliding-window radius [mm]')
    a3.grid(alpha=0.3, which='both'); a3.legend(fontsize=9)
    a3.set_title('SLIDING window: a STAIRCASE, not a curve\n'
                 '(whole shells of neighbours enter at once)', fontsize=10)

    fig.suptitle(f'{cfg.DET_NAME} — choosing the bin size — {cfg.DET_TAG} '
                 f'{cfg.RUN}\nthe bow is a ~200 mm-scale feature, so bins of '
                 '40-90 mm average hard without smearing it')
    fig.tight_layout()
    fig.savefig(os.path.join(out, f'bow_bin_scan{sfx}.png'), dpi=200,
                bbox_inches='tight')
    plt.close(fig)


def plot_raw_binned(d, cfg, out, sfx, vd, stat, bin_mm=55.0,
                    smooth_r=30.0):
    """The bow WITHOUT any surface model: raw per-pad medians, then the same
    numbers averaged in coarse SPATIAL bins.

    A single pad's peak time carries a ~stat ns error, comparable to the whole
    effect, which is why the per-pad map looks like noise and why some form of
    averaging is unavoidable. Averaging over a spatial bin assumes NOTHING about
    the shape -- unlike the quadratic fit, which is only a compact summary. If
    the bow is real it must already be visible here, at sqrt(N) better precision
    per bin.

    The per-connector constants (cable/front-end delays) are removed first by
    plain per-connector means -- also not a fit.
    """
    g = d.copy()
    g['off_c'] = g['off'] - g.groupby('connector_N')['off'].transform('mean')
    g['ax'] = (np.floor(g['pad_cx'] / bin_mm) + 0.5) * bin_mm
    g['ay'] = (np.floor(g['pad_cy'] / bin_mm) + 0.5) * bin_mm
    b = g.groupby(['ax', 'ay']).agg(m=('off_c', 'median'), n=('off_c', 'size'))
    # Decide the occupancy cut from ALL bins, not from the survivors: at a bin
    # size below the pad pitch a handful of boundary bins can still hold 5 pads
    # while the typical bin holds one, and cutting at 5 would then draw a map of
    # three squares instead of showing that no averaging is happening.
    typical = float(b['n'].median())
    min_pads = 5 if typical >= 5 else 1
    if typical < 2:
        print(f'  [warn] {bin_mm:g} mm bins hold a median of {typical:.1f} '
              'pad(s) — at or below the pad pitch, so NO averaging happens and '
              'the binned map is just the raw per-pad map redrawn')
    b = b[b['n'] >= min_pads].reset_index()
    b['err'] = stat / np.sqrt(b['n'])

    # --- sliding (running) average ------------------------------------- #
    # Same averaging idea as the hard bins, but the window is centred on every
    # pad instead of on a fixed grid, so there are no arbitrary bin edges and a
    # feature is not split between two bins. Evaluated AT the pad positions, so
    # it never extrapolates into the empty parts of the fan.
    # CAVEAT: neighbouring points share pads, so the map is smoother than its
    # information content -- the points are NOT independent. The error shown is
    # the error of one window, not of the map.
    xy = g[['pad_cx', 'pad_cy']].to_numpy()
    tree = cKDTree(xy)
    nb_idx = tree.query_ball_point(xy, r=smooth_r)
    vals = g['off_c'].to_numpy()
    g['slide'] = [float(np.median(vals[i])) for i in nb_idx]
    g['slide_n'] = [len(i) for i in nb_idx]
    g['slide_err'] = stat / np.sqrt(g['slide_n'])

    fig, axs = plt.subplots(1, 4, figsize=(24, 5.2))
    lim = float(np.nanpercentile(np.abs(g['off_c']), 97.5))
    s0 = axs[0].scatter(g['pad_cx'], g['pad_cy'], c=g['off_c'], s=14,
                        marker='s', cmap='coolwarm', vmin=-lim, vmax=lim)
    axs[0].set_title(f'RAW per-pad, connector constants removed\n'
                     f'(each pad ±{stat:.1f} ns — noise dominates locally)',
                     fontsize=10)
    plt.colorbar(s0, ax=axs[0], label='offset [ns]')

    lim2 = float(np.nanpercentile(np.abs(b['m']), 100))
    # marker AREA scales with the bin area so the squares tile the map
    s_bin = max(8.0, 380.0 * (bin_mm / 55.0) ** 2)
    s1 = axs[1].scatter(b['ax'], b['ay'], c=b['m'], s=s_bin, marker='s',
                        cmap='coolwarm', vmin=-lim2, vmax=lim2)
    if len(b) <= 120:                      # unreadable beyond this
        for _, r in b.iterrows():
            axs[1].text(r['ax'], r['ay'], f'{r["m"]:+.0f}', ha='center',
                        va='center', fontsize=6.2,
                        color='k' if abs(r['m']) < 0.6 * lim2 else 'w')
    axs[1].set_title(
        f'SAME NUMBERS, averaged in {bin_mm:.0f} mm bins '
        f'({b["n"].median():.0f} pad/bin)\n' +
        (f'no model, no fit — median ±{b["err"].median():.1f} ns per bin'
         if b['n'].median() >= 2 else
         f'±{b["err"].median():.1f} ns — NO AVERAGING: bin ≤ pad pitch'),
        fontsize=10)
    plt.colorbar(s1, ax=axs[1], label='median offset [ns]')
    lim3 = float(np.nanpercentile(np.abs(g['slide']), 100))
    s2 = axs[2].scatter(g['pad_cx'], g['pad_cy'], c=g['slide'], s=14,
                        marker='s', cmap='coolwarm', vmin=-lim3, vmax=lim3)
    axs[2].set_title(
        f'SLIDING average, r = {smooth_r:.0f} mm '
        f'({g["slide_n"].median():.0f} pads/window)\n'
        f'no bin edges; ±{g["slide_err"].median():.1f} ns per window '
        '(windows OVERLAP)', fontsize=10)
    plt.colorbar(s2, ax=axs[2], label='running median offset [ns]')
    for a in axs[:3]:
        a.set_xlabel('pad_cx [mm]')
        a.set_ylabel('pad_cy [mm]')
        a.set_aspect('equal')

    # fit-free radial profile about the ACTIVE-AREA CENTROID (no fitted apex)
    cx, cy = g['pad_cx'].mean(), g['pad_cy'].mean()
    g['u'] = np.hypot(g['pad_cx'] - cx, g['pad_cy'] - cy)
    q = g.groupby(pd.qcut(g['u'], 10, labels=False)).agg(
        u=('u', 'median'), m=('off_c', 'median'), n=('off_c', 'size'))
    q['err'] = stat / np.sqrt(q['n'])
    axs[3].errorbar(q['u'], q['m'], yerr=q['err'], fmt='o-', color='crimson',
                    lw=2, ms=7, capsize=4)
    axs[3].axhline(0, color='k', lw=0.6)
    axs[3].set_xlabel(f'distance from the active-area centroid '
                      f'({cx:.0f}, {cy:.0f}) mm')
    axs[3].set_ylabel('median offset [ns]  (later = larger gap)')
    axs[3].grid(alpha=0.3)
    sig = float(np.nanmax(np.abs(q['m'] / q['err'])))
    axs[3].set_title(f'Fit-free radial profile (10 equal-count bins)\n'
                     f'largest bin is {sig:.0f}σ from the centre value',
                     fontsize=10)
    sec = axs[3].secondary_yaxis(
        'right', functions=(lambda t: t * vd / 1000.0,
                            lambda z: z * 1000.0 / vd))
    sec.set_ylabel(f'implied gap change [mm]  (v_d = {vd:.1f} µm/ns)')

    fig.suptitle(
        f'{cfg.DET_NAME} — the bow in the MEASUREMENTS ALONE, no surface model '
        f'— {cfg.DET_TAG} {cfg.RUN}\n'
        'the quadratic fit elsewhere is a summary of this, not the evidence for '
        'it')
    fig.tight_layout()
    fig.savefig(os.path.join(out, f'bow_raw_binned{sfx}.png'), dpi=200,
                bbox_inches='tight')
    plt.close(fig)
    return b, q


def plot_residual_structure(d, resid, cfg, out, sfx, stat, r_mm, vd,
                            n_shuffle=300, seed=0):
    """Sliding map of the RESIDUAL -- what the connector constants and the one
    quadratic surface do NOT account for.

    The residual is orthogonal to the fitted model by construction, so the
    large-scale bow is gone and this isolates short-scale structure. A small
    radius is right here: r ~ 15 mm holds the pad plus its first shell of
    neighbours (7 pads), which resolves ~30 mm features.

    Significance needs care. The windows OVERLAP, so the per-window sigmas are
    not independent trials and the largest of them cannot be read off a Gaussian
    tail. Calibrate instead by SHUFFLING the residual values between pads and
    re-running the identical smoothing: that preserves the pad geometry and the
    window overlap exactly while destroying any real spatial correlation, so the
    shuffled max |sigma| is the right null for the observed one.
    """
    g = d.copy()
    g['res'] = resid
    xy = g[['pad_cx', 'pad_cy']].to_numpy()
    tree = cKDTree(xy)
    idx = tree.query_ball_point(xy, r=r_mm)
    n_win = np.array([len(i) for i in idx])
    vals = g['res'].to_numpy()
    sm = np.array([np.median(vals[i]) for i in idx])
    err = stat / np.sqrt(n_win)
    sig = sm / err

    rng = np.random.default_rng(seed)
    null_max = []
    for _ in range(n_shuffle):
        v = rng.permutation(vals)
        m = np.array([np.median(v[i]) for i in idx])
        null_max.append(np.abs(m / err).max())
    null_max = np.array(null_max)
    obs = float(np.abs(sig).max())
    p_val = float((null_max >= obs).mean())

    fig, axs = plt.subplots(1, 3, figsize=(18.5, 5.2))
    _pad_map(axs[0], g, vals, f'RAW residual per pad\n(±{stat:.1f} ns each)',
             'residual [ns]')
    lim = float(np.nanpercentile(np.abs(sm), 100))
    s1 = axs[1].scatter(g['pad_cx'], g['pad_cy'], c=sm, s=14, marker='s',
                        cmap='coolwarm', vmin=-lim, vmax=lim)
    axs[1].set_title(f'SLIDING residual, r = {r_mm:.0f} mm '
                     f'({np.median(n_win):.0f} pads/window)\n'
                     f'±{np.median(err):.1f} ns per window', fontsize=10)
    plt.colorbar(s1, ax=axs[1], label='running median residual [ns]')
    lim2 = max(3.0, float(np.nanpercentile(np.abs(sig), 99.5)))
    s2 = axs[2].scatter(g['pad_cx'], g['pad_cy'], c=sig, s=14, marker='s',
                        cmap='RdBu_r', vmin=-lim2, vmax=lim2)
    axs[2].set_title(f'Significance — max |σ| = {obs:.1f}\n'
                     f'shuffled null gives {null_max.mean():.1f} ± '
                     f'{null_max.std():.1f}  (p = {p_val:.3f})', fontsize=10)
    plt.colorbar(s2, ax=axs[2], label='window median / error [σ]')
    for a in axs:
        a.set_xlabel('pad_cx [mm]')
        a.set_ylabel('pad_cy [mm]')
        a.set_aspect('equal')
    verdict = ('REAL localised structure beyond the bow'
               if p_val < 0.05 else
               'consistent with noise once window overlap is accounted for')
    fig.suptitle(
        f'{cfg.DET_NAME} — what the bow model does NOT explain — {cfg.DET_TAG} '
        f'{cfg.RUN}\n{verdict}   (residual of connector constants + quadratic '
        f'surface; {1000*np.median(err)*vd/1e3:.0f} µm per window at '
        f'v_d = {vd:.1f} µm/ns)')
    fig.tight_layout()
    fig.savefig(os.path.join(out, f'bow_residual_structure{sfx}.png'), dpi=200,
                bbox_inches='tight')
    plt.close(fig)

    g['res_slide'] = sm
    g['res_sig'] = sig
    def _gen(v):
        return np.sqrt(max(v.var() - stat ** 2, 0.0))
    per_conn = g.groupby('connector_N')['res'].apply(_gen)
    per_half = (g.groupby(['connector_N', 'half'])['res'].apply(_gen)
                if 'half' in g.columns else None)
    return dict(obs_max_sig=obs, null_mean=float(null_max.mean()),
                null_sd=float(null_max.std()), p_value=p_val,
                pads_per_window=float(np.median(n_win)),
                err_ns=float(np.median(err))), per_conn, per_half


def plot_profile(d, off, shape, info, cfg, out, sfx, args, vd):
    u = np.hypot(d['pad_cx'] - info['apex_x'], d['pad_cy'] - info['apex_y'])
    fig, ax = plt.subplots(figsize=(8.5, 5.4))
    ax.scatter(u, off, s=9, alpha=0.30, color='0.6',
               label='per-pad offset (connector constants NOT removed)')
    q = pd.qcut(u, 12, labels=False)
    prof = pd.DataFrame({'u': u, 's': shape, 'q': q}).groupby('q').agg(
        u=('u', 'median'), s=('s', 'median'))
    ax.plot(prof['u'], prof['s'], 'o-', color='crimson', lw=2, ms=7,
            label='fitted quadratic surface (median per bin)')
    ax.axhline(0, color='k', lw=0.6)
    ax.set_xlabel(f'distance from the fitted apex '
                  f'({info["apex_x"]:.0f}, {info["apex_y"]:.0f}) mm')
    ax.set_ylabel('peak-time offset [ns]   (later = larger gap)')
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8.5)
    sec = ax.secondary_yaxis(
        'right', functions=(lambda t: t * vd / 1000.0,
                            lambda z: z * 1000.0 / vd))
    sec.set_ylabel(f'implied gap change [mm]  (v_d = {vd:.1f} µm/ns)')
    ax.set_title(f'{cfg.DET_NAME} drift-gap bow profile — {cfg.DET_TAG} '
                 f'{cfg.RUN}\nfalling away from the apex = the gap closes '
                 'toward the edges')
    fig.tight_layout()
    fig.savefig(os.path.join(out, f'bow_profile{sfx}.png'), dpi=200,
                bbox_inches='tight')
    plt.close(fig)


def vd_scaling(cfg_scan, ct, apex, sig_amp, nbin=6, min_drift=None, sfx=''):
    """Timing swing across the bow, per drift-scan point, with that point's
    apparent v_d. Returns a DataFrame(drift, n, swing_ns, vd, dz_um).

    A scan directory can hold several stage-16 CSVs from different dead-connector
    sets; take the one whose suffix matches this analysis, never just the first
    alphabetically (that silently picks a stale product).
    """
    d16 = cfg_scan.out_dir('16_drift_scan_efficiency')
    exact = os.path.join(d16, f'efficiency_vs_drift{sfx}.csv')
    if os.path.isfile(exact):
        pick = exact
    else:
        csv = sorted(glob.glob(os.path.join(d16, 'efficiency_vs_drift*.csv')))
        if not csv:
            return None
        pick = csv[0]
        print(f'  [warn] no stage-16 CSV matching {sfx!r}; falling back to '
              f'{os.path.basename(pick)}')
    print(f'  1/v_d test uses {os.path.basename(pick)}')
    sc = pd.read_csv(pick)
    vd = dict(zip(sc['x'], sc['vd_apparent_um_ns']))
    eff = dict(zip(sc['x'], sc['eff_reco']))
    m = ct[ct['mapped']]
    lut = pd.Series(m['channel_id'].to_numpy(),
                    index=pd.MultiIndex.from_arrays(
                        [m['feu'].astype(int), m['channel'].astype(int)]))
    geo = m.drop_duplicates('channel_id').set_index('channel_id')
    u = np.hypot(geo['pad_cx'] - apex[0], geo['pad_cy'] - apex[1])

    rows = []
    for sub in sorted(os.listdir(cfg_scan.run_dir)):
        if not sub.startswith('drift_scan'):
            continue
        try:
            dv = int(sub.split('_')[-1])
        except ValueError:
            continue
        # only points where the detector actually responds: a dead or
        # turning-on point has no bow to measure, just noise
        if eff.get(dv, 0) < 0.5 or not np.isfinite(vd.get(dv, np.nan)):
            continue
        if min_drift is not None and dv < min_drift:
            continue
        hits_dir = os.path.join(cfg_scan.subrun_dir(sub), 'combined_hits_root')
        parts = []
        for h in p2io.iter_hits(hits_dir, ['channel', 'amplitude', 'feu',
                                           'time_of_max'],
                                feus=ct.attrs['feus'], progress=False):
            h = h[h['amplitude'] > sig_amp]
            if not len(h):
                continue
            cid = lut.reindex(pd.MultiIndex.from_arrays(
                [h['feu'].astype(int), h['channel'].astype(int)])).to_numpy()
            ok = pd.notna(cid)
            if ok.any():
                parts.append(pd.DataFrame({
                    'u': u.reindex(cid[ok]).to_numpy(),
                    't': h['time_of_max'].to_numpy()[ok]}))
        if not parts:
            continue
        p = pd.concat(parts, ignore_index=True).dropna()
        if len(p) < 400:
            continue
        prof = p.groupby(pd.qcut(p['u'], nbin, labels=False))['t'].median()
        swing = float(prof.max() - prof.min())
        rows.append(dict(drift=dv, n=len(p), swing_ns=swing, vd=vd[dv],
                         dz_um=swing * vd[dv]))
    return pd.DataFrame(rows) if rows else None


def plot_vd_scaling(r, cfg, out, sfx, gap_mm):
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5))
    axL.plot(r['drift'], r['swing_ns'], 'o-', color='crimson', lw=2, ms=8)
    axL.set_xlabel('drift HV [V]')
    axL.set_ylabel('timing swing across the bow [ns]')
    axL.set_yscale('log')
    axL.grid(alpha=0.3, which='both')
    axL.set_title(f'Raw swing changes {r["swing_ns"].max()/r["swing_ns"].min():.1f}×\n'
                  'an ELECTRONIC offset would be flat here', fontsize=10)
    mean, std = r['dz_um'].mean(), r['dz_um'].std()
    axR.plot(r['drift'], r['dz_um'], 'o-', color='seagreen', lw=2, ms=8)
    axR.axhline(mean, color='k', ls='--', lw=1.2,
                label=f'mean {mean:.0f} µm ({100*mean/1000/gap_mm:.0f} % of the '
                      f'{gap_mm:g} mm gap)')
    axR.axhspan(mean - std, mean + std, color='seagreen', alpha=0.15,
                label=f'± {std:.0f} µm ({100*std/mean:.0f} %)')
    axR.set_xlabel('drift HV [V]')
    axR.set_ylabel('swing × $v_d$  =  gap variation [µm]')
    axR.set_ylim(0, None)
    axR.grid(alpha=0.3)
    axR.legend(fontsize=8.5)
    axR.set_title('…but swing × $v_d$ is CONSTANT\n= a fixed mechanical gap '
                  'variation', fontsize=10)
    fig.suptitle(f'{cfg.DET_NAME} — the bow scales as 1/$v_d$, so it is '
                 f'geometry, not electronics — {cfg.DET_TAG}\n'
                 '(drift scan, efficient points only; $v_d$ from each point’s '
                 'own measured peak-time spread)')
    fig.tight_layout()
    fig.savefig(os.path.join(out, f'bow_vd_scaling{sfx}.png'), dpi=200,
                bbox_inches='tight')
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('run_key', nargs='?', default='det5_long1')
    ap.add_argument('--scan-key', default=None,
                    help='drift-scan run_key for the 1/v_d test (e.g. '
                         'det5_driftscan1). Needs stage 16 to have run there.')
    ap.add_argument('--strategy', default='reverse')
    ap.add_argument('--sig-amp', type=float, default=300.0,
                    help='signal band: only hits above this feed the per-pad '
                         'amplitude (stage 19 uses the same for the time)')
    ap.add_argument('--resid-r', type=float, default=15.0,
                    help='sliding-window radius [mm] for the RESIDUAL map. '
                         '15 mm = pad + first shell of neighbours (7 pads), the '
                         'finest window that averages anything.')
    ap.add_argument('--smooth-r', type=float, default=30.0,
                    help='radius [mm] of the sliding (running) average window, '
                         'centred on every pad. ~30 mm holds a similar number '
                         'of pads to a 55 mm hard bin but has no bin edges.')
    ap.add_argument('--bin-mm', type=float, default=55.0,
                    help='spatial bin size [mm] for the model-free map. Must '
                         'exceed the pad pitch (~11.8 mm) to average anything; '
                         '40-90 mm is the useful range for a ~200 mm bow.')
    ap.add_argument('--min-hits', type=int, default=20,
                    help='pads with fewer signal hits are dropped')
    ap.add_argument('--vd', type=float, default=None,
                    help='drift velocity [um/ns] for the ns->mm axis; default '
                         'is taken from the drift scan if --scan-key is given, '
                         'else 27.2')
    ap.add_argument('--veto-sparks', action=argparse.BooleanOptionalAction,
                    default=True)
    args = ap.parse_args()

    cfg = qa.RUNS[args.run_key]
    print(cfg)
    sfx = cfg.product_suffix(args.veto_sparks)
    out = cfg.out_dir('23_drift_gap_bow')
    gap_mm = cfg.DRIFT_GAP_MM
    print(f'  drift gap {gap_mm:g} mm (cfg.DRIFT_GAP_MM)')

    ct = channel_table(cfg, args.strategy)
    sub = cfg.subrun_dir(cfg.SUB_RUN)
    hits_dir = os.path.join(sub, 'combined_hits_root')

    # stage 19 must have run: it owns the per-pad peak time and the split-half
    # statistical error this stage tests the residual against.
    t19 = glob.glob(os.path.join(cfg.out_dir('19_timing_surface'),
                                 f'timing_surface{sfx}.csv'))
    if not t19:
        print(f'!! no stage-19 CSV at {cfg.out_dir("19_timing_surface")} '
              f'(expected timing_surface{sfx}.csv) — run 19_timing_surface.py '
              'first')
        return
    d19 = pd.read_csv(t19[0])

    bad = None
    if args.veto_sparks:
        sv = ps.SparkVeto.from_csv(os.path.join(sub, 'hv_monitor.csv'), cfg)
        bad = sv.vetoed_ids_from_hits(hits_dir, ct.attrs['feus'],
                                      min_amp=cfg.MIN_AMP)
        print(f'  spark veto: {len(bad)} events removed')

    drop = p2io.drop_pads_for(cfg, ct, hits_dir=hits_dir)
    if len(drop):
        print(f'  hot/noisy pads dropped (same set as stage 19): {len(drop)}')
    amp = pad_median_amplitude(hits_dir, ct, args.sig_amp, exclude_events=bad,
                               min_amp=cfg.MIN_AMP, drop_pads=drop)
    geo = ct[ct['mapped']].drop_duplicates('channel_id')[
        ['channel_id', 'pad_cx', 'pad_cy', 'connector_N', 'half', 'feu']]
    d = (d19[['channel_id', 'n', 'median_ns', 'median_even_ns',
              'median_odd_ns']]
         .merge(amp, on='channel_id', how='inner')
         .merge(geo, on='channel_id', how='inner'))
    d = d[(d['n'] >= args.min_hits) & (d['n_amp'] >= args.min_hits) &
          d['median_ns'].notna() & d['amp_median'].notna()].copy()
    if len(d) < 100:
        print(f'!! only {len(d)} usable pads — not enough to fit a surface')
        return
    d['off'] = d['median_ns'] - d['median_ns'].median()
    print(f'  {len(d)} pads with >= {args.min_hits} signal hits; '
          f'offset std {d["off"].std():.2f} ns')

    # Split-half error. CAREFUL: std(even-odd)/sqrt(2) is the error of a
    # HALF-sample median (that is what stage 19 prints and uses internally).
    # What we fit here is `median_ns`, the FULL-sample median, whose error is
    # smaller by another sqrt(2). Comparing the residual against the half-sample
    # number would understate the leftover structure and can make a model look
    # perfect when it is not.
    sh = d.dropna(subset=['median_even_ns', 'median_odd_ns'])
    sd = float(np.nanstd(sh['median_even_ns'] - sh['median_odd_ns']))
    stat_half = sd / np.sqrt(2.0)
    stat = sd / 2.0
    print(f'  per-pad statistical error: {stat:.2f} ns on the full-sample '
          f'median (stage 19 quotes {stat_half:.2f} ns for a half sample)')

    # --- 1. decomposition ---------------------------------------------- #
    st, rt, it = fit_bow(d, 'off')
    gen_tot = np.sqrt(max(d['off'].var() - stat ** 2, 0.0))
    gen = np.sqrt(max(rt.var() - stat ** 2, 0.0))
    # The honest figure of merit is not "is the residual pure noise" (with
    # enough hits it never quite is) but HOW MUCH of the real structure the one
    # smooth surface accounts for.
    frac = 1.0 - (gen / gen_tot) ** 2 if gen_tot > 0 else np.nan
    print(f'\n  MODEL connector constants + quadratic surface:')
    print(f'    genuine structure in the map : {gen_tot:.2f} ns')
    print(f'    residual std {rt.std():.2f} ns vs {stat:.2f} ns statistical '
          f'-> genuine residual {gen:.2f} ns')
    print(f'    => the smooth surface explains {100*frac:.0f} % of the genuine '
          'variance')
    print(f'    quadratic part: {it["kind"].upper()}, anisotropy '
          f'{it["anisotropy"]:.1f} (5+ = closer to a cylindrical bend than a '
          f'symmetric sag), peak-to-peak {it["pp"]:.1f} ns')
    print('    connector constants [ns]: ' +
          ', '.join(f'c{k}={v:+.1f}' for k, v in it['conn_const'].items()))

    # --- 2. amplitude cross-check --------------------------------------- #
    sa, ra, ia = fit_bow(d, 'amp_median')
    amed = float(d['amp_median'].median())
    bt = bootstrap_fit(d, 'off')
    ba = bootstrap_fit(d, 'amp_median')
    # Correlating the two FITTED surfaces would be close to circular (both are
    # smooth quadratics over the same pads), so do it on the raw measurements
    # binned by distance from the apex, with the connector constants removed.
    d['u_apex'] = np.hypot(d['pad_cx'] - it['apex_x'], d['pad_cy'] - it['apex_y'])
    prof = raw_profile(d, 'u_apex')
    r_raw = float(np.corrcoef(prof['t'], prof['a'])[0, 1])
    print(f'\n  AMPLITUDE cross-check:')
    print(f'    {ia["kind"].upper()}, apex ({ia["apex_x"]:.0f}, '
          f'{ia["apex_y"]:.0f}) mm — {np.hypot(it["apex_x"]-ia["apex_x"], it["apex_y"]-ia["apex_y"]):.0f} mm '
          'from the timing apex')
    print(f'    BOOTSTRAP — what is actually constrained:')
    print(f'      curvature sign  : dome in {100*bt["dome_frac"]:.0f} % of '
          f'timing resamples, {100*ba["dome_frac"]:.0f} % of amplitude ones')
    print(f'      peak-to-peak    : {it["pp"]:.1f} +- {bt["pp_sd"]:.1f} ns')
    print(f'      apex position   : +-({bt["apex_sx"]:.0f}, {bt["apex_sy"]:.0f}) mm '
          '— NOT constrained (fan-shaped area + elliptical surface); do not '
          'read it as the centre of the bow')
    print(f'    shape peak-to-peak {ia["pp"]:.0f} ADC = '
          f'{100*ia["pp"]/amed:.1f} % of the {amed:.0f} ADC median')
    print(f'    correlation of the RAW binned profiles (connector constants '
          f'removed, {len(prof)} bins): r = {r_raw:+.3f}')
    print(f'    [the two fitted surfaces correlate at r = '
          f'{float(np.corrcoef(st, sa)[0, 1]):+.3f}, but that is partly '
          'circular — quote the raw number]')

    # --- 3. 1/v_d scaling ------------------------------------------------ #
    scan = None
    vd_used = args.vd
    if args.scan_key:
        cfg_scan = qa.RUNS[args.scan_key]
        scan = vd_scaling(cfg_scan, channel_table(cfg_scan, args.strategy),
                          (it['apex_x'], it['apex_y']), args.sig_amp,
                          sfx=cfg_scan.product_suffix(args.veto_sparks))
        if scan is not None and len(scan) >= 3:
            print(f'\n  1/v_d SCALING across {len(scan)} drift points:')
            for _, s in scan.iterrows():
                print(f'    drift {s.drift:.0f} V: swing {s.swing_ns:6.1f} ns  '
                      f'v_d {s.vd:5.1f} um/ns  -> {s.dz_um:6.0f} um')
            print(f'    raw swing varies '
                  f'{scan["swing_ns"].max()/scan["swing_ns"].min():.1f}x; '
                  f'swing*v_d = {scan["dz_um"].mean():.0f} +- '
                  f'{scan["dz_um"].std():.0f} um '
                  f'({100*scan["dz_um"].std()/scan["dz_um"].mean():.0f} % scatter)')
            if vd_used is None:
                vd_used = float(scan['vd'].max())
        else:
            print('\n  [skip] 1/v_d scaling: not enough usable drift points')
            scan = None
    if vd_used is None:
        vd_used = 34.0          # det5 @ drift 820 V, 4 mm gap; see --vd
        print(f'  [warn] no --scan-key and no --vd: assuming v_d = '
              f'{vd_used:.1f} um/ns for the ns->mm axis only')

    # --- implied gap variation ------------------------------------------- #
    dz_t = it['pp'] * vd_used / 1000.0
    dz_a = ia['pp'] / amed * gap_mm
    if frac < 0.6 or it['kind'] != 'dome' or r_raw < 0.5:
        print(f'\n  *** THE BOW MODEL DOES NOT DESCRIBE THIS DETECTOR ***')
        print(f'    the smooth surface explains only {100*frac:.0f} % of the '
              f'genuine variance, the curvature is a {it["kind"]}, and the '
              f'amplitude correlation is {r_raw:+.3f}. The numbers below are '
              'printed for completeness but should NOT be quoted as a gap '
              'variation.')
    print(f'\n  IMPLIED DRIFT-GAP VARIATION on the {gap_mm:g} mm gap:')
    print(f'    from timing    (v_d {vd_used:.1f} um/ns): {dz_t:.2f} mm '
          f'= {100*dz_t/gap_mm:.0f} %')
    print(f'    from amplitude (primaries ~ path length): {dz_a:.2f} mm '
          f'= {100*dz_a/gap_mm:.0f} %')
    if scan is not None:
        print(f'    from 1/v_d scaling: {scan["dz_um"].mean()/1000:.2f} mm '
              f'= {100*scan["dz_um"].mean()/1000/gap_mm:.0f} %')
    print(f'    -> the {it["kind"]} means the gap is '
          f'{"LARGEST in the middle: the mesh/readout assembly bows AWAY from the drift electrode" if it["kind"] == "dome" else "SMALLEST in the middle"}')

    # --- products --------------------------------------------------------- #
    plot_decomposition(d, d['off'].to_numpy(), st, rt, it, stat, cfg, out, sfx,
                       args, gen, frac)
    plot_timing_vs_amplitude(d, st, sa, it, ia, amed, prof, r_raw,
                             bt, ba, cfg, out, sfx, args)
    plot_profile(d, d['off'].to_numpy(), st, it, cfg, out, sfx, args, vd_used)
    plot_3d_bow(d, st, it, cfg, out, sfx, vd_used, gap_mm)
    bbin, bprof = plot_raw_binned(d, cfg, out, sfx, vd_used, stat,
                                  bin_mm=args.bin_mm,
                                  smooth_r=args.smooth_r)
    _xy = d[['pad_cx', 'pad_cy']].to_numpy()
    pitch = float(np.median(cKDTree(_xy).query(_xy, k=2)[0][:, 1]))
    gg = d.copy()
    gg['off_c'] = gg['off'] - gg.groupby('connector_N')['off'].transform('mean')
    bscan = bin_size_scan(gg, stat)
    rscan = radius_scan(gg, stat)
    plot_bin_scan(bscan, rscan, pitch, cfg, out, sfx)
    bscan.to_csv(os.path.join(out, f'bow_bin_scan{sfx}.csv'), index=False)
    rscan.to_csv(os.path.join(out, f'bow_radius_scan{sfx}.csv'), index=False)
    print(f'\n  BIN-SIZE SCAN (pad pitch {pitch:.1f} mm):')
    print('    bin[mm]  bins  pads/bin  err[ns]  span[ns]  max_sig')
    for _, r in bscan.iterrows():
        flag = '  <- below pad pitch' if r.bin_mm < pitch else ''
        print(f'    {r.bin_mm:7.0f} {r.n_bins:5.0f} {r.pads_per_bin:9.1f} '
              f'{r.err_ns:8.2f} {r.span_ns:9.1f} {r.max_sig:8.1f}{flag}')
    print(f'\n  SLIDING-RADIUS SCAN:')
    print('     r[mm]  pads/window  err[ns]  max_sig')
    for _, r in rscan.iterrows():
        flag = ('  <- 1 pad: no averaging' if r.pads_per_window < 2 else '')
        print(f'    {r.radius_mm:6.0f} {r.pads_per_window:12.1f} '
              f'{r.err_ns:8.2f} {r.max_sig:8.1f}{flag}')

    rinfo, per_conn, per_half = plot_residual_structure(
        d, rt, cfg, out, sfx, stat, args.resid_r, vd_used)
    print(f'\n  RESIDUAL STRUCTURE (sliding r = {args.resid_r:g} mm, '
          f'{rinfo["pads_per_window"]:.0f} pads/window, '
          f'±{rinfo["err_ns"]:.2f} ns):')
    print(f'    max |sigma| observed {rinfo["obs_max_sig"]:.1f}; '
          f'label-shuffled null {rinfo["null_mean"]:.1f} ± '
          f'{rinfo["null_sd"]:.1f} -> p = {rinfo["p_value"]:.3f}')
    print('    genuine residual RMS per connector [ns]: ' +
          ', '.join(f'c{int(k)}={v:.1f}' for k, v in per_conn.items()))
    if per_half is not None:
        worst = per_half.sort_values(ascending=False).head(4)
        print('    worst connector-HALVES [ns]: ' +
              ', '.join(f'c{int(k[0])}-{k[1]}={v:.1f}' for k, v in worst.items()))

    sig = float(np.nanmax(np.abs(bprof['m'] / bprof['err'])))
    print(f'\n  MODEL-FREE CHECK (spatial bins, no surface fit):')
    print(f'    {len(bbin)} bins of >=5 pads, median error '
          f'{bbin["err"].median():.2f} ns/bin')
    print(f'    binned offsets span {bbin["m"].min():+.1f} .. '
          f'{bbin["m"].max():+.1f} ns')
    print(f'    fit-free radial profile: largest bin {sig:.0f} sigma from the '
          'centre value')
    if scan is not None:
        plot_vd_scaling(scan, cfg, out, sfx, gap_mm)
        scan.to_csv(os.path.join(out, f'bow_vd_scaling{sfx}.csv'), index=False)

    pd.DataFrame([dict(
        run_key=args.run_key, n_pads=len(d), drift_gap_mm=gap_mm,
        stat_err_ns=stat, stat_err_halfsample_ns=stat_half,
        resid_std_ns=float(rt.std()), genuine_total_ns=gen_tot,
        genuine_resid_ns=gen, frac_variance_explained=frac,
        kind=it['kind'], apex_x_mm=it['apex_x'], apex_y_mm=it['apex_y'],
        anisotropy=it['anisotropy'], timing_pp_ns=it['pp'],
        amp_apex_x_mm=ia['apex_x'], amp_apex_y_mm=ia['apex_y'],
        amp_pp_adc=ia['pp'], amp_median_adc=amed,
        amp_pp_frac=ia['pp'] / amed, corr_raw_binned=r_raw,
        corr_fitted_surfaces=float(np.corrcoef(st, sa)[0, 1]),
        apex_err_x_mm=bt['apex_sx'], apex_err_y_mm=bt['apex_sy'],
        pp_bootstrap_sd_ns=bt['pp_sd'], dome_frac_timing=bt['dome_frac'],
        dome_frac_amp=ba['dome_frac'],
        vd_um_per_ns=vd_used, gap_var_timing_mm=dz_t, gap_var_amp_mm=dz_a,
        gap_var_vdscaling_mm=(scan['dz_um'].mean() / 1000.0
                              if scan is not None else np.nan),
    )]).to_csv(os.path.join(out, f'bow_fit{sfx}.csv'), index=False)
    print(f'\nWritten to: {out}')


if __name__ == '__main__':
    main()
