#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Chamber efficiency of a scintillator-referenced threshold scan, defined as on the cosmic bench (06/10/12 stages).

Input: the scint_efficiency.py json of each run (per probe chamber: every scintillator-confirmed tag with its
predicted point, the distance to the nearest probe hit, whether every nearby pad is live; and the probe's pads).

Definitions (cosmic_bench_analysis/06_efficiency_maps.py, adapted):
  read-out area   the predicted point lies at least --margin mm inside the footprint of the MAPPED probe pads (a point
                  is on the footprint when it is within the pad half-diagonal of a pad centre; "margin mm inside" =
                  the whole circle of that radius around it is on the footprint). Tracks predicted outside the
                  connectors we read are not in the denominator (acceptance, not efficiency). The margin covers the
                  pointing error of the two-chamber prediction (bench: sub-mm M3 pointing, margin 0).
  operating eff   read-out-area tags with a probe hit within R of the prediction. Masked and silent pads stay in the
                  footprint: a muon lost on them is a loss, as NOISY_PADS on the bench. Headline number.
  live-pad eff    the same, without the tags whose nearby pads are not all live (scint_efficiency.py live area); a
                  diagnostic of the response of the pads still on, always quoted with the fraction of tags it keeps.
  breakdown       read-out-area tags: hit within R / fired elsewhere (hit, but further than R) / no hit.
R and the margin are checked for a plateau (validation figure per run, as 12_validation.py); the random-hit rate in
the window is < 0.1 % (json 'random') and not subtracted.

Products in OUT_DIR:
  validation_run<N>.png      per probe: eff vs R (+ hit anywhere), eff vs margin; residual distribution
  maps_<label>.png           per probe: adaptive-kernel map (each point averages the --target nearest tags, so every
                             point has the same binomial error; blank where the kernel is wider than --max-kernel)
                             and the kernel radius (= local resolution of the map); not-live pads marked
  per_vmm.png                operating efficiency per VMM (tags assigned to the chip of the nearest pad) against K
  scan.png, scan.txt, scan.json   operating / live-pad efficiency, breakdown and live-pad fraction against K

Runs with the same K can be pooled for the maps (--pool, default on): their tags are summed.

Usage: scint_eff_study.py OUT_DIR scintor_run181.json scintor_run184.json ... [--r 30] [--margin 10]
"""
import argparse
import json
import os
import re
import sys
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

from event_display import style, INK, MUTED, GRID, COL  # noqa: E402

DETS = ((3, 'det3 P2_IN'), (2, 'det2 P2_MID'), (1, 'det1 P2_OUT'))
NOTLIVE = '#c0392b'
HALF_DIAG = 0.5 * np.hypot(12.07, 11.61)   # P2 pad half-diagonal [mm], the bench ACTIVE_R
R_SCAN = np.array([5, 8, 10, 12, 15, 20, 25, 30, 40, 50, 60, 80, 100, 150])
M_SCAN = np.array([0, 5, 10, 15, 20, 25, 30, 40])
NCIRC = 24


def k_of(title):
    m = re.search(r'_K([pm]?)(\d+)', title)
    return None if not m else (-1 if m[1] == 'm' else 1) * int(m[2])


def binom(k, n):
    e = k / n if n else np.nan
    return e, (np.sqrt(max(e * (1 - e), 0) / n) if n else np.nan)


class Probe:
    """Tags and pads of one probe chamber in one run."""

    def __init__(self, P):
        t = P['tags_xy']
        if 'd' not in t:
            raise SystemExit('json has no per-tag distance: re-run scint_efficiency.py (version with tags_xy["d"])')
        self.x, self.y = np.array(t['x'], float), np.array(t['y'], float)
        self.d = np.array(t['d'], float)
        self.any = np.array(t['any'], bool)
        self.live = np.array(t['live'], bool)
        p = P['pads']
        self.px, self.py = np.array(p['x'], float), np.array(p['y'], float)
        self.plive = np.array(p['live'], bool)
        self.pvmm = np.array(p['vmm'], int)
        self.tree = cKDTree(np.column_stack([self.px, self.py]))
        _, self.nearest = self.tree.query(np.column_stack([self.x, self.y]))
        self._depth = None

    def depth(self):
        """How far inside the footprint each predicted point is: largest scanned margin whose circle is on it."""
        if self._depth is None:
            self._depth = inside_depth(self.tree, self.x, self.y)
        return self._depth


def on_footprint(tree, x, y):
    return tree.query(np.column_stack([x, y]))[0] <= HALF_DIAG


def inside_depth(tree, x, y):
    """Per point: the largest margin of M_SCAN for which the centre and a circle of that radius are on the footprint
    (-1 = the point itself is off the footprint)."""
    dep = np.where(on_footprint(tree, x, y), 0.0, -1.0)
    ang = np.linspace(0, 2 * np.pi, NCIRC, endpoint=False)
    ok = dep >= 0
    for m in M_SCAN[1:]:
        cx = (x[:, None] + m * np.cos(ang)[None, :]).ravel()
        cy = (y[:, None] + m * np.sin(ang)[None, :]).ravel()
        allon = on_footprint(tree, cx, cy).reshape(len(x), NCIRC).all(axis=1)
        ok = ok & allon
        dep[ok] = m
    return dep


def measure(pr, R, margin):
    """Operating / live-pad efficiency, breakdown and outside fraction of one probe."""
    ins = pr.depth() >= margin
    n = int(ins.sum())
    w = pr.d <= R
    out = {'n_tags': int(len(pr.x)), 'n_area': n, 'frac_outside': 1 - n / max(len(pr.x), 1)}
    out['operating'] = binom(int((ins & w).sum()), n)
    out['any'] = binom(int((ins & pr.any).sum()), n)
    nl = int((ins & pr.live).sum())
    out['live'] = binom(int((ins & pr.live & w).sum()), nl)
    out['live_kept'] = nl / n if n else np.nan
    out['breakdown'] = {'within_R': float((ins & w).sum() / max(n, 1)),
                        'fired_elsewhere': float((ins & pr.any & ~w).sum() / max(n, 1)),
                        'no_hit': float((ins & ~pr.any).sum() / max(n, 1))}
    out['live_pad_frac'] = float(pr.plive.mean())
    vm = {}
    for v in np.unique(pr.pvmm):
        s = ins & (pr.pvmm[pr.nearest] == v)
        e, er = binom(int((s & w).sum()), int(s.sum()))
        vm[int(v)] = {'n': int(s.sum()), 'eff': e, 'err': er,
                      'pads': int((pr.pvmm == v).sum()), 'pads_not_live': int(((pr.pvmm == v) & ~pr.plive).sum())}
    out['per_vmm'] = vm
    return out


def validation_fig(run, title, probes, R, margin, path):
    fig, axs = plt.subplots(3, 3, figsize=(15, 11), constrained_layout=True)
    for col, (dd, nm) in enumerate(DETS):
        pr = probes[nm]
        ins = pr.depth() >= margin
        n = ins.sum()
        # eff vs R
        ax = axs[0, col]
        style(ax)
        e = np.array([((pr.d <= r) & ins).sum() / n for r in R_SCAN])
        er = np.sqrt(e * (1 - e) / n)
        ax.errorbar(R_SCAN, e, yerr=er, fmt='o-', color=COL[dd], ms=4, lw=1.2, label=f'hit within R ({n} tags)')
        ea = (pr.any & ins).sum() / n
        ax.axhline(ea, color=MUTED, lw=1, ls='--', label=f'hit anywhere {ea:.3f}')
        ax.axvline(R, color=INK, lw=0.8, ls=':', label=f'R used = {R:g} mm')
        ax.axvline(HALF_DIAG, color=GRID, lw=0.8)
        ax.set_xscale('log')
        ax.set_xlabel('match radius R [mm]', fontsize=8)
        ax.set_ylabel(f'efficiency (margin {margin:g} mm)', fontsize=8)
        ax.set_title(f'{nm}: efficiency against R', fontsize=9, color=INK)
        ax.legend(fontsize=7, frameon=False, loc='lower right')
        # eff vs margin
        ax = axs[1, col]
        style(ax)
        rows = []
        for m in M_SCAN:
            s = pr.depth() >= m
            k = int(((pr.d <= R) & s).sum())
            rows.append((m, *binom(k, int(s.sum())), s.mean()))
        rows = np.array(rows)
        ax.errorbar(rows[:, 0], rows[:, 1], yerr=rows[:, 2], fmt='o-', color=COL[dd], ms=4, lw=1.2,
                    label=f'hit within {R:g} mm')
        ax.axvline(margin, color=INK, lw=0.8, ls=':', label=f'margin used = {margin:g} mm')
        ax2 = ax.twinx()
        ax2.plot(rows[:, 0], rows[:, 3], 's--', color=MUTED, ms=3, lw=0.8)
        ax2.set_ylim(0, 1.05)
        ax2.set_ylabel('fraction of all tags kept', fontsize=8, color=MUTED)
        ax2.tick_params(labelsize=7, colors=MUTED)
        ax.set_xlabel('margin inside the read-out footprint [mm]', fontsize=8)
        ax.set_ylabel('efficiency', fontsize=8)
        ax.set_title(f'{nm}: efficiency against the margin', fontsize=9, color=INK)
        ax.legend(fontsize=7, frameon=False, loc='lower right')
        # residual distribution
        ax = axs[2, col]
        style(ax)
        dd_ = pr.d[ins & pr.any]
        ax.hist(np.clip(dd_, 0, 150), bins=np.arange(0, 152, 2), color=COL[dd], alpha=0.8)
        q = np.percentile(dd_, [50, 68, 90]) if len(dd_) else [np.nan] * 3
        for v, lab in zip(q, ('50 %', '68 %', '90 %')):
            ax.axvline(v, color=INK, lw=0.7, ls=':')
            ax.text(v, ax.get_ylim()[1] * 0.92, f' {lab} {v:.0f}', fontsize=7, color=INK)
        ax.set_xlabel('distance from the prediction to the nearest probe hit [mm] (150 = overflow)', fontsize=8)
        ax.set_ylabel('tags with a probe hit', fontsize=8)
        ax.set_title(f'{nm}: pointing (read-out area, tags with a hit)', fontsize=9, color=INK)
    fig.suptitle(f'{title}: choice of the match radius R and of the margin (plateau check, as cosmic-bench stage 12)',
                 fontsize=11)
    fig.savefig(path, dpi=110)
    plt.close(fig)


def maps_fig(label, plist, R, margin, target, max_kernel, path):
    """plist: {probe name: [Probe, ...]} (runs pooled)."""
    fig, axs = plt.subplots(2, 3, figsize=(17, 10), constrained_layout=True)
    cmap = plt.get_cmap('viridis').copy()
    cmap.set_bad('white')
    for col, (dd, nm) in enumerate(DETS):
        ps = plist[nm]
        x = np.concatenate([p.x for p in ps])
        y = np.concatenate([p.y for p in ps])
        w = np.concatenate([p.d <= R for p in ps])
        ins = np.concatenate([p.depth() >= margin for p in ps])
        x, y, w = x[ins], y[ins], w[ins]
        p0 = ps[-1]
        plive = np.all([p.plive for p in ps], axis=0)
        gx = np.arange(p0.px.min() - 10, p0.px.max() + 10, 4.0)
        gy = np.arange(p0.py.min() - 10, p0.py.max() + 10, 4.0)
        GX, GY = np.meshgrid(gx, gy, indexing='ij')
        g = np.column_stack([GX.ravel(), GY.ravel()])
        on = p0.tree.query(g)[0] <= HALF_DIAG
        E = np.full(len(g), np.nan)
        RAD = np.full(len(g), np.nan)
        if len(x) >= target:
            tt = cKDTree(np.column_stack([x, y]))
            dist, idx = tt.query(g[on], k=target)
            rad = dist[:, -1]
            e = w[idx].mean(axis=1)
            e[rad > max_kernel] = np.nan
            E[on], RAD[on] = e, rad
        E, RAD = E.reshape(GX.shape), RAD.reshape(GX.shape)
        ext = [gx[0] - 2, gx[-1] + 2, gy[0] - 2, gy[-1] + 2]
        n, k = len(x), int(w.sum())
        e0, er0 = binom(k, n)
        for row, (Z, lab, cm, vmin, vmax, cl) in enumerate((
                (E, f'efficiency (hit within {R:g} mm)', cmap, 0, 1, 'efficiency, each point = mean of the '
                 f'{target} nearest tags (+-{np.sqrt(0.25 / target):.2f} at most)'),
                (RAD, 'kernel radius [mm]', plt.get_cmap('magma_r'), 0, max_kernel, 'kernel radius [mm] '
                 f'(radius holding {target} tags; > {max_kernel:g} mm blank above)'))):
            ax = axs[row, col]
            style(ax)
            ax.grid(False)
            im = ax.imshow(Z.T, origin='lower', extent=ext, cmap=cm, vmin=vmin, vmax=vmax, aspect='equal',
                           interpolation='nearest')
            ax.scatter(p0.px[plive], p0.py[plive], s=1.2, color=GRID if row == 0 else MUTED, lw=0, zorder=3)
            ax.scatter(p0.px[~plive], p0.py[~plive], s=10, marker='x', color=NOTLIVE, lw=0.8, zorder=4,
                       label=f'pad not live ({(~plive).sum()} of {len(plive)})')
            ax.set_xlim(ext[0], ext[1])
            ax.set_ylim(ext[2], ext[3])
            ax.set_xlabel('x [mm]', fontsize=8)
            ax.set_ylabel('y [mm]', fontsize=8)
            if row == 0:
                ax.set_title(f'{nm}: operating efficiency {e0:.3f} +- {er0:.3f} ({n} tags in the read-out area)',
                             fontsize=9, color=INK)
            else:
                ax.set_title(f'{nm}: {lab}', fontsize=9, color=INK)
            ax.legend(fontsize=7, frameon=False, loc='upper left')
            if col == 2:
                fig.colorbar(im, ax=axs[row, :], shrink=0.8, label=cl)
    fig.suptitle(f'{label}: adaptive-kernel efficiency map at the predicted point (as cosmic-bench stage 10); '
                 f'read-out area with a {margin:g} mm margin', fontsize=11)
    fig.savefig(path, dpi=110)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('out_dir')
    ap.add_argument('jsons', nargs='+')
    ap.add_argument('--r', type=float, default=30, help='match radius [mm]')
    ap.add_argument('--margin', type=float, default=10, help='margin inside the read-out footprint [mm] (one of %s)'
                    % ','.join(str(m) for m in M_SCAN))
    ap.add_argument('--target', type=int, default=150, help='tags per kernel of the maps')
    ap.add_argument('--max-kernel', type=float, default=60)
    ap.add_argument('--no-validation', action='store_true')
    a = ap.parse_args()
    if a.margin not in M_SCAN:
        raise SystemExit(f'--margin must be one of {list(M_SCAN)}')
    if not os.path.isdir(a.out_dir):
        os.makedirs(a.out_dir)

    rows = []
    bylab = defaultdict(lambda: defaultdict(list))
    for f in a.jsons:
        E = json.load(open(f))
        run = int(re.findall(r'run_?(\d+)', E['title'] + ' ' + os.path.basename(f))[0])
        K = k_of(E['title'])
        probes = {nm: Probe(E[f'probe_{nm}']) for _, nm in DETS}
        row = {'run': run, 'K': K, 'title': E['title'], 'n_trig': E['n_trig'],
               'vetoed_frac': (E.get('bad_times') or {}).get('vetoed_frac'),
               'random': {nm: E['chambers'][nm]['random'] for _, nm in DETS}, 'R_mm': a.r, 'margin_mm': a.margin}
        for _, nm in DETS:
            row[nm] = measure(probes[nm], a.r, a.margin)
        rows.append(row)
        if not a.no_validation:
            validation_fig(run, E['title'], probes, a.r, a.margin, os.path.join(a.out_dir, f'validation_run{run}.png'))
        lab = f'K{K:+d}' if K is not None else f'run{run}'
        for _, nm in DETS:
            bylab[lab][nm].append(probes[nm])
        bylab[lab]['runs'].append(run)
    rows.sort(key=lambda r: r['run'])

    for lab, pl in bylab.items():
        runs = '+'.join(str(r) for r in sorted(pl['runs']))
        maps_fig(f'{lab} (run {runs})', pl, a.r, a.margin, a.target, a.max_kernel,
                 os.path.join(a.out_dir, f'maps_{lab}.png'))

    # table
    lines = [f'operating efficiency = read-out-area tags (margin {a.margin:g} mm) with a probe hit within {a.r:g} mm; '
             'masked / silent pads count as losses',
             'live-pad efficiency  = the same without tags near a not-live pad ("kept" = fraction of read-out-area tags left)',
             'breakdown of the read-out-area tags: within R / fired elsewhere / no hit; "out" = tags predicted outside the area', '']
    hdr = f'{"run":>4} {"K":>4} {"trig":>5} ' + ' | '.join(f'{nm:^62s}' for _, nm in DETS)
    sub = ' ' * 15 + ' | '.join(f'{"out":>4} {"operating":>12} {"live-pad":>12} {"kept":>4} {"in R/else/none":>15} {"pads":>4}'
                                for _ in DETS)
    lines += [hdr, sub]
    for r in rows:
        cells = []
        for _, nm in DETS:
            m = r[nm]
            b = m['breakdown']
            cells.append(f'{m["frac_outside"]:4.2f} {m["operating"][0]:6.3f}+-{m["operating"][1]:.3f} '
                         f'{m["live"][0]:6.3f}+-{m["live"][1]:.3f} {m["live_kept"]:4.2f} '
                         f'{b["within_R"]:4.2f}/{b["fired_elsewhere"]:4.2f}/{b["no_hit"]:4.2f} {m["live_pad_frac"]:4.2f}')
        lines.append(f'{r["run"]:>4} {r["K"]!s:>4} {r["n_trig"]:>5} ' + ' | '.join(cells))
    txt = '\n'.join(lines)
    print(txt)
    open(os.path.join(a.out_dir, 'scan.txt'), 'w').write(txt + '\n')
    json.dump(rows, open(os.path.join(a.out_dir, 'scan.json'), 'w'), indent=1, default=float)

    # scan figure: operating + live-pad efficiency against K, breakdown, live pads
    rk = [r for r in rows if r['K'] is not None]
    fig, axs = plt.subplots(1, 4, figsize=(18, 4.4), constrained_layout=True)
    for j, (dd, nm) in enumerate(DETS):
        ax = axs[j]
        style(ax)
        K = np.array([r['K'] for r in rk], float)
        op = np.array([r[nm]['operating'] for r in rk])
        lv = np.array([r[nm]['live'] for r in rk])
        an = np.array([r[nm]['any'] for r in rk])
        ax.errorbar(K, op[:, 0], yerr=op[:, 1], fmt='o', color=COL[dd], ms=6, lw=1.2,
                    label=f'operating (hit within {a.r:g} mm, masked pads = losses)')
        ax.errorbar(K + 1.0, lv[:, 0], yerr=lv[:, 1], fmt='o', color=COL[dd], mfc='white', ms=6, lw=1.2,
                    label='live pads only (diagnostic)')
        ax.plot(K - 1.0, an[:, 0], '_', color=MUTED, ms=10, mew=1.5, label='hit anywhere in the chamber')
        ax.set_ylim(0, 1)
        ax.invert_xaxis()
        ax.set_xlabel('threshold: safe + K [DAC]', fontsize=8)
        ax.set_ylabel(f'efficiency (read-out area, {a.margin:g} mm margin)', fontsize=8)
        ax.set_title(nm, fontsize=9, color=INK)
        if j == 0:
            ax.legend(fontsize=7, frameon=False, loc='lower left')
    ax = axs[3]
    style(ax)
    for dd, nm in DETS:
        K = np.array([r['K'] for r in rk], float)
        ax.plot(K, [r[nm]['live_pad_frac'] for r in rk], 'o', color=COL[dd], ms=5, label=f'{nm}: live pads / pads read')
        ax.plot(K, [r[nm]['live_kept'] for r in rk], 'D', color=COL[dd], mfc='white', ms=5,
                label=f'{nm}: tags kept by the live-pad cut')
    ax.set_ylim(0, 1.02)
    ax.invert_xaxis()
    ax.set_xlabel('threshold: safe + K [DAC]', fontsize=8)
    ax.set_ylabel('fraction', fontsize=8)
    ax.set_title('masking', fontsize=9, color=INK)
    ax.legend(fontsize=6.5, frameon=False, loc='lower left')
    fig.suptitle('Threshold scan, scintillator reference ((top OR) AND bottom), one point per run', fontsize=11)
    fig.savefig(os.path.join(a.out_dir, 'scan.png'), dpi=110)
    plt.close(fig)

    # per-VMM figure
    fig, axs = plt.subplots(1, 3, figsize=(17, 4.6), constrained_layout=True)
    for j, (dd, nm) in enumerate(DETS):
        ax = axs[j]
        style(ax)
        vmms = sorted({v for r in rk for v in r[nm]['per_vmm']})
        cm = plt.get_cmap('tab10')
        for i, v in enumerate(vmms):
            pts = [(r['K'], r[nm]['per_vmm'][v]) for r in rk if v in r[nm]['per_vmm'] and r[nm]['per_vmm'][v]['n'] >= 20]
            if not pts:
                continue
            K = np.array([p[0] for p in pts], float) + (i - len(vmms) / 2) * 0.5
            ax.errorbar(K, [p[1]['eff'] for p in pts], yerr=[p[1]['err'] for p in pts], fmt='o-', color=cm(i % 10),
                        ms=4, lw=0.9, label=f'VMM {v}')
        ax.set_ylim(0, 1)
        ax.invert_xaxis()
        ax.set_xlabel('threshold: safe + K [DAC]', fontsize=8)
        ax.set_ylabel(f'operating efficiency (hit within {a.r:g} mm)', fontsize=8)
        ax.set_title(f'{nm}: per VMM (tags assigned to the chip of the nearest pad; >= 20 tags)', fontsize=9, color=INK)
        ax.legend(fontsize=7, frameon=False, loc='lower left', ncol=3)
    fig.savefig(os.path.join(a.out_dir, 'per_vmm.png'), dpi=110)
    plt.close(fig)
    print('wrote', a.out_dir)


if __name__ == '__main__':
    sys.exit(main())
