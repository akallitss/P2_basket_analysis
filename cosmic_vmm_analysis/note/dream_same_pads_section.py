#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Section "The same pads with DREAM" of the note "Where do the missed muons go?" (22 Nov 2026).

Compares, pad by pad in the board frame, the VMM lab stand (REF configuration, per predicted pad) with the DREAM
cosmic bench (M3 reference tracks) on the read-out area of the stand (connectors c3-c5), and for det3 with the SPS beam
(DREAM eff_nominal_1 and VMM run_46 on P2_OUT) on the pads that both the beam spot and the stand cover.

Inputs
  --vmm      silent_tests.py output with the per-pad found share (pads rows [vmm, ch, cx, cy, n, silent, found])
  --dream    directory with dream_pads_<run_key>.csv from dream_pad_export.py
  --beam     sps_beam_analysis/dream_vmm_comparison/data/compare_dream_vmm_P2_OUT.csv (+ run_conditions_P2_OUT.json)
  --map      Detector_Mapping/P2_BASKET/P2_BASKET_mapping.csv
Output: an HTML fragment (--out), inserted by muon_funnel_page.py (or --note: in place, before "Open points").
"""
import argparse
import html
import json
import os
import re

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

# lab stand: station -> chamber; DREAM bench run used for the maps (and the extra run quoted in the table)
LAB = [('P2_IN', 'det3', 'det3_initial1', None), ('P2_MID', 'det2', 'det2_long1', None),
       ('P2_OUT', 'det1', 'det1_long5', 'det1_initial1')]
LAB_HV = {'P2_IN': '420 / 720', 'P2_MID': '430 / 730', 'P2_OUT': '430 / 730'}      # run_metadata_corrections.json
BENCH = {'det3_initial1': ('16 Jul', '420 / 820', 'before its 17 Jul gain collapse'),
         'det2_long1': ('9 Jul', '430 / 600', '10 h long run'),
         'det1_long5': ('19 Jul', '415 / 615', '12 h long run, after the repair'),
         'det1_initial1': ('19 Jul', '430 / 740', '1 h, after the repair')}
CONNS = (3, 4, 5)
PAD_W, PAD_H = 11.4, 11.0
MIN_N = 15            # pads with fewer VMM muons (or DREAM rays) are drawn empty
CELL = 36.0           # mm, cells of the difference map (~3 x 3 pads)
STYLE = """<style>
  :root { --dr:#2a78d6; --vm:#eb6834; --neg:#2a78d6; }
  @media (prefers-color-scheme: dark) { :root { --dr:#3987e5; --vm:#d95926; --neg:#3987e5; } }
  .dr-d { fill: var(--dr); } .vm-d { fill: var(--vm); }
  .dr-o { fill: none; stroke: var(--dr); stroke-width: 1.8; } .vm-o { fill: none; stroke: var(--vm); stroke-width: 1.8; }
  .pair { stroke: var(--ink2); stroke-width: 1; opacity: .6; }
  table.sp td { white-space: nowrap; } table.sp td.wr { white-space: normal; }
</style>"""


def wavg(v, w):
    w = np.asarray(w, float)
    return float((np.asarray(v, float) * w).sum() / w.sum()) if w.sum() > 0 else float('nan')


def pc(x, d=1):
    return '–' if x is None or not np.isfinite(x) else f'{100 * x:.{d}f} %'


def load(a):
    m = pd.read_csv(a.map)
    tree = cKDTree(m[['pad_cx', 'pad_cy']].to_numpy())
    J = json.load(open(a.vmm))
    V, D = {}, {}
    for st, det, key, key2 in LAB:
        P = pd.DataFrame(J['REF'][st]['pads'], columns=['vmm', 'ch', 'cx', 'cy', 'n', 'sil', 'found'])
        dd, ii = tree.query(P[['cx', 'cy']].to_numpy())
        if dd.max() > 0.5:
            raise SystemExit(f'{st}: a VMM pad centre is {dd.max():.2f} mm from every Gerber pad')
        P['cid'] = m.channel_id.to_numpy()[ii]
        P['conn'] = P.cid // 128 + 1
        P['anyh'] = 1 - P.sil
        V[st] = P
        for k in (key, key2):
            if k:
                D[k] = pd.read_csv(os.path.join(a.dream, f'dream_pads_{k}.csv'))
    ro = m[(m.channel_id // 128 + 1).isin(CONNS)][['channel_id', 'pad_cx', 'pad_cy']]
    B = pd.read_csv(a.beam)
    B = B[B.use]
    cond = json.load(open(os.path.join(os.path.dirname(a.beam), 'run_conditions_P2_OUT.json')))
    return m, ro, V, D, B, cond


# ---------------------------------------------------------------- maps
def bounds(xy, pad=10):
    return xy.pad_cx.min() - pad, xy.pad_cx.max() + pad, xy.pad_cy.min() - pad, xy.pad_cy.max() + pad


def svg_map(ro, val, bnd, title, sub, w=290, diff=False, cells=None, label='missed'):
    """val: {channel_id: (value, n, tooltip)}; value = missed share (0..0.30+) or, with diff, VMM - DREAM missed share."""
    x0, x1, y0, y1 = bnd
    sc = w / (x1 - x0)
    h = (y1 - y0) * sc
    X = lambda x: (x - x0) * sc
    Y = lambda y: h - (y - y0) * sc
    o = [f'<svg viewBox="0 0 {w} {h + 40:.0f}" class="dia" role="img" aria-label="{html.escape(title)}">']
    if cells is None:
        for cid, x, y in ro.itertuples(index=False):
            if not (x0 <= x <= x1 and y0 <= y <= y1):
                continue
            q = val.get(cid)
            if q:
                op = min(abs(q[0]) / 0.30, 1.0)
                col = 'var(--neg)' if q[0] < 0 else 'var(--k-sil)'
                st = f'style="fill:{col};fill-opacity:{max(op, 0.04):.2f};stroke:var(--pad-line);stroke-width:.5"'
                tt = f'<title>{html.escape(q[2])}</title>'
            else:
                st, tt = 'class="pad"', f'<title>pad {cid}: fewer than {MIN_N} muons</title>'
            o.append(f'<rect x="{X(x) - PAD_W / 2 * sc:.1f}" y="{Y(y) - PAD_H / 2 * sc:.1f}" width="{PAD_W * sc:.1f}" '
                     f'height="{PAD_H * sc:.1f}" {st}>{tt}</rect>')
    else:
        for cid, x, y in ro.itertuples(index=False):
            o.append(f'<rect x="{X(x) - PAD_W / 2 * sc:.1f}" y="{Y(y) - PAD_H / 2 * sc:.1f}" width="{PAD_W * sc:.1f}" '
                     f'height="{PAD_H * sc:.1f}" class="pad"/>')
        for (cx, cy), (v, n, tip) in cells.items():
            op = min(abs(v) / 0.30, 1.0)
            col = 'var(--neg)' if v < 0 else 'var(--k-sil)'
            o.append(f'<rect x="{X(cx - CELL / 2) + .6:.1f}" y="{Y(cy + CELL / 2) + .6:.1f}" width="{CELL * sc - 1.2:.1f}" '
                     f'height="{CELL * sc - 1.2:.1f}" style="fill:{col};fill-opacity:{max(op, 0.05):.2f};'
                     f'stroke:var(--ink2);stroke-width:.4"><title>{html.escape(tip)}</title></rect>')
    o.append(f'<text x="{w / 2}" y="{h + 16:.0f}" class="lbl" text-anchor="middle">{html.escape(title)}</text>'
             f'<text x="{w / 2}" y="{h + 32:.0f}" class="lbl small" text-anchor="middle">{html.escape(sub)}</text></svg>')
    return ''.join(o)


def svg_scale(w=900, H=34, diff=False):
    o = [f'<svg viewBox="0 0 {w} {H}" class="dia" role="img" aria-label="colour scale">']
    vals = [-0.10, -0.05, 0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30] if diff else [0, .05, .10, .15, .20, .25, .30]
    xs = 300
    for i, v in enumerate(vals):
        col = 'var(--neg)' if v < 0 else 'var(--k-sil)'
        lab = f'{100 * v:+.0f}' if diff else f'{100 * v:.0f}'
        o.append(f'<rect x="{xs + i * 52}" y="4" width="50" height="14" style="fill:{col};'
                 f'fill-opacity:{max(min(abs(v) / 0.30, 1), 0.04):.2f};stroke:var(--pad-line)"/>'
                 f'<text x="{xs + 25 + i * 52}" y="31" class="lbl small" text-anchor="middle">'
                 f'{lab}{"+" if i == len(vals) - 1 else ""}</text>')
    txt = ('VMM minus DREAM missed share (points)' if diff else 'missed share = 1 − any-hit efficiency')
    o.append(f'<text x="{xs - 10}" y="16" class="lbl small" text-anchor="end">{txt}</text>')
    if not diff:
        o.append(f'<rect x="{xs + 7 * 52 + 20}" y="4" width="16" height="14" class="pad"/>'
                 f'<text x="{xs + 7 * 52 + 42}" y="16" class="lbl small">&lt; {MIN_N} muons</text>')
    o.append('</svg>')
    return ''.join(o)


def pad_vals_vmm(P):
    return {int(r.cid): (1 - r.anyh, r.n, f'pad {r.cid} (VMM {r.vmm} ch {r.ch}): any hit {100 * r.anyh:.0f} %, found '
                         f'{100 * r.found:.0f} % of {r.n:.0f} muons') for r in P.itertuples() if r.n >= MIN_N}


def pad_vals_dream(Dk):
    return {int(r.channel_id): (1 - r.any, r.n, f'pad {r.channel_id} (connector {r.connector_N}): any hit '
                                f'{100 * r.any:.0f} %, within {100 * r.within:.0f} % of {r.n} rays')
            for r in Dk.itertuples() if r.n >= MIN_N and r.connector_N in CONNS}


def diff_cells(P, Dk, ro):
    """cells of CELL mm: VMM-tag-weighted any-hit of both read-outs on the pads of the cell that both measured."""
    k = P.merge(Dk, left_on='cid', right_on='channel_id', suffixes=('', '_d'))
    gx = np.floor(k.cx / CELL).astype(int)
    gy = np.floor(k.cy / CELL).astype(int)
    out = {}
    for (ix, iy), g in k.groupby([gx, gy]):
        if g.n.sum() < 60:
            continue
        v_ = 1 - wavg(g.anyh, g.n)
        d_ = 1 - wavg(g['any'], g.n)
        out[((ix + .5) * CELL, (iy + .5) * CELL)] = (
            v_ - d_, g.n.sum(), f'{len(g)} pads, {g.n.sum():.0f} VMM muons / {g.n_d.sum()} DREAM rays: any hit VMM '
            f'{100 * (1 - v_):.1f} %, DREAM {100 * (1 - d_):.1f} % (difference {100 * (d_ - v_):+.1f} points)')
    return out


# ---------------------------------------------------------------- per chip
def chip_table(V, D):
    rows = []
    for st, det, key, _ in LAB:
        k = V[st].merge(D[key], left_on='cid', right_on='channel_id', suffixes=('', '_d'))
        for vmm, g in k.groupby('vmm'):
            rows.append(dict(st=st, det=det, vmm=int(vmm), conn=int(g.conn.mode()[0]), n=g.n.sum(), nd=int(g.n_d.sum()),
                             v_found=wavg(g.found, g.n), v_any=wavg(g.anyh, g.n),
                             d_within=wavg(g.within, g.n), d_any=wavg(g['any'], g.n)))
    return pd.DataFrame(rows)


def svg_chips(C, w=900, H=300):
    x0, x1, yt, yb = 60, w - 15, 20, H - 60
    lo, hi = 0.65, 1.0
    Y = lambda v: yb - (yb - yt) * (min(max(v, lo), hi) - lo) / (hi - lo)
    n = len(C)
    gap = 18
    step = (x1 - x0 - 2 * gap) / n
    o = [f'<svg viewBox="0 0 {w} {H}" class="dia" role="img" aria-label="efficiency per VMM chip, DREAM and VMM">']
    for t in (0.7, 0.8, 0.9, 1.0):
        o.append(f'<line x1="{x0}" y1="{Y(t):.1f}" x2="{x1}" y2="{Y(t):.1f}" class="grid"/>'
                 f'<text x="{x0 - 6}" y="{Y(t) + 4:.1f}" class="lbl small" text-anchor="end">{100 * t:.0f} %</text>')
    for i, r in enumerate(C.itertuples()):
        grp = ['P2_IN', 'P2_MID', 'P2_OUT'].index(r.st)
        x = x0 + gap * grp + (i + .5) * step
        o.append(f'<line x1="{x:.1f}" y1="{Y(r.d_any):.1f}" x2="{x:.1f}" y2="{Y(r.v_any):.1f}" class="pair"/>')
        tip = (f'{r.st} ({r.det}) VMM {r.vmm}, connector {r.conn}: VMM any hit {100 * r.v_any:.1f} % / found '
               f'{100 * r.v_found:.1f} % ({r.n:.0f} muons); DREAM on the same pads any {100 * r.d_any:.1f} % / within '
               f'{100 * r.d_within:.1f} % ({r.nd} rays)')
        for v, cls, rad in ((r.d_any, 'dr-d', 5), (r.d_within, 'dr-o', 4), (r.v_any, 'vm-d', 5), (r.v_found, 'vm-o', 4)):
            o.append(f'<circle cx="{x:.1f}" cy="{Y(v):.1f}" r="{rad}" class="{cls}"><title>{html.escape(tip)}</title></circle>')
        o.append(f'<text x="{x:.1f}" y="{yb + 16}" class="lbl small" text-anchor="middle">{r.vmm}</text>'
                 f'<text x="{x:.1f}" y="{yb + 29}" class="lbl small" text-anchor="middle">c{r.conn}</text>')
    for grp, (st, det, _, _) in enumerate(LAB):
        idx = np.nonzero((C.st == st).to_numpy())[0]
        xa = x0 + gap * grp + idx.min() * step
        xb = x0 + gap * grp + (idx.max() + 1) * step
        o.append(f'<line x1="{xa + 4:.1f}" y1="{yb + 36}" x2="{xb - 4:.1f}" y2="{yb + 36}" class="axis"/>'
                 f'<text x="{(xa + xb) / 2:.1f}" y="{yb + 52}" class="lbl" text-anchor="middle">{st} = {det}</text>')
    o.append(f'<text x="{x0 - 44}" y="{(yt + yb) / 2:.0f}" class="lbl small" text-anchor="middle" '
             f'transform="rotate(-90 {x0 - 44} {(yt + yb) / 2:.0f})">efficiency</text></svg>')
    return ''.join(o)


# ---------------------------------------------------------------- beam
def beam_overlap(V, D, B):
    P = V['P2_IN']
    c = B.merge(P, left_on='pad_id', right_on='cid', suffixes=('_b', '')).merge(
        D['det3_initial1'], left_on='pad_id', right_on='channel_id', how='left', suffixes=('', '_d'))
    return c


def section(a):
    m, ro, V, D, B, cond = load(a)
    out = [STYLE, '<h2>9. The same pads with DREAM (cosmic bench and SPS beam)</h2>']
    out.append(
        '<p>The three stand chambers were all measured with DREAM before, on the M3 cosmic bench at Saclay (reference '
        'tracks from the M3 telescope), and det1 and det3 also at the SPS beam. The stand reads connectors <b>c3, c4, '
        'c5</b> of each chamber. These are the same connectors as <code>connector_N</code> 3-5 of the bench pipeline '
        '(both are 1-based, <code>channel_id = (N − 1)·128 + strip</code>; Gerber <code>connector_number</code> 2-4), '
        'so the comparison can be made <b>pad by pad</b> in the board frame. The stand has P2_IN = det3, P2_MID = det2 '
        'and P2_OUT = det1.</p>')
    out.append(
        '<p>Two definitions are compared. <b>Any hit</b> means the chamber fired anywhere: DREAM <code>has_any</code>, '
        'and for the VMM 1 − SILENT (no live hit in the chamber within ±100 ns of the tag). This one does not depend on '
        'how good the pointing is. <b>With position</b> means DREAM <code>within</code> (a cluster within 40 mm of the M3 '
        'track) and VMM FOUND (a hit within 26 / 14 / 26 mm and the time window). DREAM numbers on the VMM pads are '
        'weighted by the number of VMM muons on each pad, so both measurements describe the same area.</p>')

    # ---- table per connector
    T = ['<div style="overflow-x:auto"><table class="sp"><tr><th>stand position (chamber)</th><th>read-out</th><th>run (HV mesh / '
         'drift, V)</th><th>c3</th><th>c4</th><th>c5</th><th><b>c3-c5</b></th><th>c3-c5 any hit</th></tr>']
    summ = {}
    for st, det, key, key2 in LAB:
        P = V[st]
        cells = []
        for cn in CONNS:
            g = P[P.conn == cn]
            cells.append(wavg(g.found, g.n))
        summ[st] = dict(v_found=wavg(P.found, P.n), v_any=wavg(P.anyh, P.n))
        T.append(f'<tr><td class="wr" rowspan="{3 if key2 else 2}">{st} ({det})</td><td>VMM lab, REF</td><td class="wr">19-22 Nov '
                 f'({LAB_HV[st]})</td>' + ''.join(f'<td>{pc(x)}</td>' for x in cells) +
                 f'<td><b>{pc(summ[st]["v_found"])}</b></td><td>{pc(summ[st]["v_any"])}</td></tr>')
        for k in (key, key2):
            if not k:
                continue
            Dk = D[k]
            Dk3 = Dk[Dk.connector_N.isin(CONNS)]
            cl = [wavg(Dk[Dk.connector_N == cn].within, Dk[Dk.connector_N == cn].n) for cn in CONNS]
            day, hv, what = BENCH[k]
            s_w, s_a = wavg(Dk3.within, Dk3.n), wavg(Dk3['any'], Dk3.n)
            if k == key:
                kk = P.merge(Dk, left_on='cid', right_on='channel_id', suffixes=('', '_d'))
                summ[st].update(d_within=wavg(kk.within, kk.n), d_any=wavg(kk['any'], kk.n), d_rays=int(Dk3.n.sum()))
            T.append(f'<tr><td>DREAM bench</td><td class="wr"><code>{k}</code>, {day} ({hv}; {what})</td>' +
                     ''.join(f'<td>{pc(x)}</td>' for x in cl) + f'<td><b>{pc(s_w)}</b></td><td>{pc(s_a)}</td></tr>')
    T.append('</table></div>')
    out.append(''.join(T))
    out.append(
        '<p class="faint">VMM: FOUND per connector of the predicted pad (REF, all REF runs dumped by 22 Nov: '
        f'{V["P2_IN"].n.sum():.0f} / {V["P2_MID"].n.sum():.0f} / {V["P2_OUT"].n.sum():.0f} muons for IN / MID / OUT). '
        'These come from the per-tag dumps (scintillator-confirmed tags, random tag pairs subtracted, as in 5d), '
        'which are not exactly the funnel sample of section 3, so they differ from it by up to ~1 point. '
        'DREAM: <code>within</code> per connector of the pad nearest to the M3 track, spark-vetoed, inside the active '
        'area; rays in each connector are counted as the bench saw them (not reweighted). The bench HV is not the stand '
        'HV: det2 had its drift at 600 V (730 V on the stand) and det3 at 820 V (720 V). For det1 the 430 / 740 V run is '
        'the closer match, and it agrees with the long run at 415 / 615 V to within a point.</p>')

    # ---- weighted on the same pads
    S = ['<table><tr><th>on the pads the VMM tested (VMM-muon weighted)</th><th>DREAM any hit</th><th>VMM any hit</th>'
         '<th>gap</th><th>DREAM with position</th><th>VMM with position</th><th>gap</th></tr>']
    for st, det, key, _ in LAB:
        s = summ[st]
        S.append(f'<tr><td>{st} ({det}) vs <code>{key}</code></td><td>{pc(s["d_any"])}</td><td>{pc(s["v_any"])}</td>'
                 f'<td><b>{100 * (s["d_any"] - s["v_any"]):.1f} points</b></td><td>{pc(s["d_within"])}</td>'
                 f'<td>{pc(s["v_found"])}</td><td><b>{100 * (s["d_within"] - s["v_found"]):.1f} points</b></td></tr>')
    S.append('</table>')
    out.append(''.join(S))

    # ---- maps
    bnd = bounds(ro)
    out.append('<h3>9a. Maps of the missed share on the read-out area</h3>')
    out.append(
        f'<p>Each square is one pad of c3-c5, coloured by the share of muons that left <b>no hit at all</b> (1 − any '
        f'hit), on the same scale for both read-outs. Pads with fewer than {MIN_N} muons are left empty. The third '
        f'column groups pads into {CELL:.0f} mm cells (at least 60 VMM muons each) and shows VMM minus DREAM: red = the VMM '
        'misses more. Hover a pad or cell for its numbers. A VMM pad has 15-40 muons, so ±6-8 % per pad is statistics '
        'alone. det3 has only ~13 bench rays per pad (a 1 h run, the day before its gain collapse).</p>')
    for st, det, key, _ in LAB:
        P, Dk = V[st], D[key]
        day, hv, _w = BENCH[key]
        s = summ[st]
        cells = diff_cells(P, Dk, ro)
        out.append('<div class="row3">' +
                   svg_map(ro, pad_vals_dream(Dk), bnd, f'{det}: DREAM bench', f'{key}, {hv} V · any hit {pc(s["d_any"])}') +
                   svg_map(ro, pad_vals_vmm(P), bnd, f'{st} = {det}: VMM lab (REF)',
                           f'{LAB_HV[st]} V · any hit {pc(s["v_any"])}') +
                   svg_map(ro, {}, bnd, 'VMM − DREAM, missed share', f'gap {100 * (s["d_any"] - s["v_any"]):.1f} points '
                           'on the same pads', cells=cells) + '</div>')
    out.append(svg_scale())
    out.append(svg_scale(diff=True))

    # ---- per chip
    C = chip_table(V, D)
    out.append('<h3>9b. Per VMM chip: the same strip of pads, two read-outs</h3>')
    out.append(f'<figure>{svg_chips(C)}<figcaption>Each column is one VMM chip of the stand: its 64 pads form a strip '
               'along the band. Blue = DREAM bench on those pads, orange = VMM lab; filled = any hit, open = with '
               'position. DREAM is weighted with the VMM muons of each pad. Hover for the numbers.</figcaption></figure>')
    worst = C.assign(g=C.d_any - C.v_any).sort_values('g', ascending=False).iloc[0]
    cin = C[C.st == 'P2_IN']
    rest_in = cin[cin.vmm != worst.vmm]
    cmid, cout = C[C.st == 'P2_MID'], C[C.st == 'P2_OUT']
    out.append(
        '<ol>'
        f'<li><b>P2_MID (det2): the gap is the same everywhere.</b> On every chip the VMM is '
        f'{100 * (cmid.d_any - cmid.v_any).min():.0f}-{100 * (cmid.d_any - cmid.v_any).max():.0f} points (any hit) below '
        'DREAM on the same pads. DREAM sees no weak region there (92-95 % any hit on every chip), and the per-pad VMM and '
        'DREAM values are not correlated (r ≈ 0 over the pads with ≥ 25 entries in both). So the stand\'s chamber-wide '
        'silent share of P2_MID (section 5d) is <b>not a weak area of det2 that DREAM also sees</b>. Either the VMM sits '
        'further into the det2 signal than the DREAM threshold, or det2 has less gas gain on the stand than on the bench. '
        'The bench drift was 600 V, not 730 V, so the gain may differ.</li>'
        f'<li><b>P2_OUT (det1): a smaller, also flat gap.</b> {100 * (cout.d_any - cout.v_any).min():.0f}-'
        f'{100 * (cout.d_any - cout.v_any).max():.0f} points per chip. DREAM reaches 97-99 % any hit on these pads, so '
        'det1 is the best chamber with both read-outs.</li>'
        f'<li><b>P2_IN (det3): one chip.</b> VMM {worst.vmm} (connector c{worst.conn}) finds {pc(worst.v_found)} '
        f'({pc(worst.v_any)} any hit) where DREAM has {pc(worst.d_within)} ({pc(worst.d_any)}) on the same pads. The '
        f'other P2_IN chips are at {100 * rest_in.v_found.min():.0f}-{pc(rest_in.v_found.max(), 0)} found, level with or '
        'above DREAM. The c5 deficit is there in all configurations of the weekend, not only in REF. DREAM found nothing '
        'wrong with det3 there, and the beam (9c) read the same pads with other VMM chips at ~89 %. So this is a property '
        f'of chip {worst.vmm} or its channels on the stand (threshold, trim, a cable), not of the chamber. It is the '
        '"bad chip" of P2_IN in section 5b.</li>'
        '</ol>')

    # ---- beam
    c = beam_overlap(V, D, B)
    out.append('<h3>9c. det3 at the SPS beam, on the same pads</h3>')
    vb, db = cond['vmm'], cond['dream']
    hv = cond['hv']
    w = c.n
    lab_chips = sorted(c.vmm.unique())
    beam_chips = sorted(c.vmm_b.unique())
    out.append(
        f'<p>At the beam det3 sat at P2_OUT. The beam spot covered {len(B)} pads of its connectors c5 and c6, and '
        f'{len(c)} of them (all on c5) are also read on the stand, where det3 is P2_IN. On those pads there are four '
        f'measurements of the same chamber. The beam used DREAM run <code>{db["run"]}</code> ({db["start"]}) and VMM '
        f'<code>{vb["run"]}</code> ({vb["start"]}, gain {vb["gain_mv_per_fc"]} mV/fC, sdt {vb["sdt"]}), both at mesh '
        f'{hv["p2_mesh_v"]} / drift {hv["p2_drift_v"]} V with {cond["gas"]}. The beam probe radius was 15 mm, against '
        'tracks good to 0.3 mm. The beam VMM read these pads with chips '
        f'{", ".join(map(str, beam_chips))}; the stand reads them with chips {", ".join(map(str, lab_chips))}.</p>')
    rows = [('DREAM, SPS beam', wavg(c.eff_d, w), f'{c.n_track_d.sum():,} tracks'),
            ('VMM, SPS beam', wavg(c.eff_v, w), f'{c.n_track_v.sum():,} tracks'),
            ('DREAM, cosmic bench (det3_initial1)', wavg(c.within, w.where(c.within.notna(), 0)),
             f'{int(c.n_d.sum())} rays; any hit {pc(wavg(c["any"], w.where(c["any"].notna(), 0)))}'),
            ('VMM, lab stand (REF)', wavg(c.found, w), f'{w.sum():.0f} muons; any hit {pc(wavg(c.anyh, w))}')]
    out.append('<table><tr><th>det3, the ' + str(len(c)) + ' common pads (weighted by the stand muons)</th><th>efficiency '
               '(with position)</th><th>sample</th></tr>' +
               ''.join(f'<tr><td>{a_}</td><td><b>{pc(b_)}</b></td><td>{c_}</td></tr>' for a_, b_, c_ in rows) + '</table>')
    # mini maps on the beam-spot region
    reg = ro[ro.channel_id.isin(c.pad_id)]
    bb = bounds(reg, pad=14)
    roz = ro[(ro.pad_cx.between(bb[0], bb[1])) & (ro.pad_cy.between(bb[2], bb[3]))]
    mk = lambda col, n, txt: {int(r.pad_id): (1 - getattr(r, col), getattr(r, n),
                                              f'pad {r.pad_id}: {txt} {100 * getattr(r, col):.0f} % of {getattr(r, n):.0f}')
                              for r in c.itertuples() if np.isfinite(getattr(r, col)) and getattr(r, n) >= 1}
    out.append('<div class="row4">' +
               svg_map(roz, mk('eff_d', 'n_track_d', 'DREAM beam'), bb, 'DREAM, beam', f'{pc(rows[0][1])}', w=220) +
               svg_map(roz, mk('eff_v', 'n_track_v', 'VMM beam'), bb, 'VMM, beam', f'{pc(rows[1][1])}', w=220) +
               svg_map(roz, mk('within', 'n_d', 'DREAM bench within'), bb, 'DREAM, bench', f'{pc(rows[2][1])}', w=220) +
               svg_map(roz, mk('found', 'n', 'VMM stand found'), bb, 'VMM, stand', f'{pc(rows[3][1])}', w=220) + '</div>')
    out.append('<p class="faint">Missed share (1 − efficiency with position), same colour scale as 9a. Here every pad '
               'with at least one entry is coloured: the stand has 9-43 muons per pad and the bench 9-25 rays, so '
               'single pads are not significant there; the beam has thousands.</p>')
    r_vb = np.corrcoef(c.found, c.eff_v)[0, 1]
    out.append(
        '<ol>'
        f'<li><b>DREAM is {100 * (rows[0][1] - rows[2][1]):.0f} points higher at the beam</b> ({pc(rows[0][1])}) than on '
        f'the bench ({pc(rows[2][1])}). The conditions were different: mesh 450 V instead of 420 V, another gas, '
        'and a 0.3 mm track instead of M3 with a 40 mm match radius. The bench run was also taken the day before '
        'det3\'s gain collapse, and it has only ~400 rays on these pads.</li>'
        f'<li><b>The VMM gap is larger on the stand than at the beam</b>: {100 * (rows[0][1] - rows[1][1]):.1f} points '
        f'at the beam, {100 * (rows[2][1] - rows[3][1]):.1f} points on the stand. At the beam, the gap on this chamber '
        'came from the VMM threshold cutting into the Landau of the weak pads (the beam study, "Why the VMM loses the '
        f'weak pads"). On the stand, chip {worst.vmm} adds a loss of its own (9b).</li>'
        f'<li><b>Pad by pad</b> the stand cannot follow the beam\'s gain map: r = {r_vb:+.2f} between the stand and beam '
        'VMM per-pad efficiencies, with ~25 stand muons per pad. Testing whether the weak beam pads are also weak on '
        'the stand needs ~10x more stand statistics on these pads.</li>'
        '</ol>')
    out.append(
        '<div class="callout"><b>What the DREAM comparison adds.</b> On the same pads DREAM sees all three chambers at '
        f'{100 * summ["P2_IN"]["d_any"]:.0f}-{pc(summ["P2_OUT"]["d_any"], 0)} any hit, with no weak region in the stand\'s '
        'read-out area. The VMM stand is below DREAM by '
        f'{100 * (summ["P2_MID"]["d_any"] - summ["P2_MID"]["v_any"]):.0f} points any hit '
        f'({100 * (summ["P2_MID"]["d_within"] - summ["P2_MID"]["v_found"]):.0f} with position) on det2 (P2_MID), '
        'uniformly; by '
        f'{100 * (summ["P2_OUT"]["d_any"] - summ["P2_OUT"]["v_any"]):.0f} '
        f'({100 * (summ["P2_OUT"]["d_within"] - summ["P2_OUT"]["v_found"]):.0f}) points on det1 (P2_OUT), uniformly; and '
        f'on det3 (P2_IN) level with DREAM except on chip {worst.vmm}. So the P2_MID silent share is not a weak area of '
        'det2, and the P2_IN loss is one chip. Before a mesh scan of P2_MID, two cheap checks: '
        f'(1) the threshold and trims of chip {worst.vmm} against its P2_IN neighbours; '
        '(2) the det2 bench HV (drift 600 V) against the stand HV (730 V).</div>')
    out.append(
        '<p class="faint">Code: <code>cosmic_vmm_analysis/dream_pad_export.py</code> (DREAM per pad, from the stage-06 '
        'ray lists), <code>silent_tests.py</code> with the per-pad found share (scratch copy run on lxplus 22 Nov, '
        'groups REF/ALL), <code>note/dream_same_pads_section.py</code> (this section). Beam per-pad table: '
        '<code>sps_beam_analysis/dream_vmm_comparison/data/compare_dream_vmm_P2_OUT.csv</code>.</p>')
    return '\n'.join(out), summ, C, c


UPDATE = ('<div class="callout ok"><b>Update 22 Nov: the same pads with DREAM.</b> Section 9 compares the stand, pad by '
          'pad, with the DREAM cosmic bench, and for det3 with the SPS beam. DREAM sees no weak region in the stand\'s '
          'read-out area (c3-c5). The P2_MID (det2) gap to DREAM is uniform over the chamber, P2_OUT (det1) has a '
          'smaller uniform gap, and the P2_IN (det3) loss is one VMM chip.</div>')


def insert_html(doc, frag):
    """put the section before "Open points" (renumbering 9, 10 -> 10, 11) and the update callout at the top; idempotent."""
    doc = re.sub(r'\n?<!-- same-pads:(begin|top) -->.*?<!-- same-pads:(end|topend) -->\n?', '\n', doc, flags=re.S)
    doc = doc.replace('<h2>10. References</h2>', '<h2>11. References</h2>').replace(
        '<h2>9. Open points</h2>', '<h2>10. Open points</h2>')
    i = doc.index('<h2>10. Open points</h2>')
    doc = doc[:i] + '<!-- same-pads:begin -->\n' + frag + '\n<!-- same-pads:end -->\n' + doc[i:]
    j = doc.index('<p class="faint">Data: lab cosmic stand')
    return doc[:j] + '<!-- same-pads:top -->\n' + UPDATE + '\n<!-- same-pads:topend -->\n' + doc[j:]


def insert(note, frag):
    doc = insert_html(open(note).read(), frag)
    open(note, 'w').write(doc)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--vmm', required=True)
    ap.add_argument('--dream', required=True)
    ap.add_argument('--beam', required=True)
    ap.add_argument('--map', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--note', default=None)
    a = ap.parse_args()
    frag, summ, C, c = section(a)
    open(a.out, 'w').write(frag)
    for k, v in summ.items():
        print(k, {kk: round(vv, 4) for kk, vv in v.items()})
    print(C.round(3).to_string())
    if a.note:
        insert(a.note, frag)
        print('inserted into', a.note)
