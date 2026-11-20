#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build the note body "Where do the missed muons go?" from muon_funnel_summary.py output (20 Nov 2026).

Self-contained HTML (inline SVG, no scripts from outside). Replaces the <div class="note-body"> of NOTE_HTML.
Usage: muon_funnel_page.py SUMMARY.json NOTE_HTML [GAIN_SCAN_JSON [POINTING_JSON [OLD_SUMMARY_JSON]]]
  POINTING_JSON: pointing_resolution.py output (the per-probe search radius and edge margin SUMMARY was made with);
  OLD_SUMMARY_JSON: the same funnel with the fixed 15 mm radius and no margin, for the before / after table.
"""
import html
import json
import math
import re
import os
import sys
from collections import Counter

S = json.load(open(sys.argv[1]))
G = S['groups']
REF, ALL = G['REF'], G['ALL']
PADS = S['pads']
EX = S['examples']
GS = json.load(open(sys.argv[3])) if len(sys.argv) > 3 else None
PT = json.load(open(sys.argv[4])) if len(sys.argv) > 4 else None
OLD = json.load(open(sys.argv[5]))['groups'] if len(sys.argv) > 5 else None
RAD = {p: (PT['probes'][p]['probe_r'] if PT else 15.0) for p in ('P2_IN', 'P2_MID', 'P2_OUT')}
MARG = {p: (PT['probes'][p]['margin_mm'] if PT else 0.0) for p in ('P2_IN', 'P2_MID', 'P2_OUT')}
PILLAR = 0.0378  # bulk pillars, mask V1 (Saclay bulk: det1-4, i.e. all three stand chambers): 0.8 mm on a 4 mm grid
ST = json.load(open(sys.argv[6])) if len(sys.argv) > 6 else None   # silent_tests.py output
FIGDIR = os.path.basename(sys.argv[2]).replace('.html', '_files/')
PROBES = ['P2_MID', 'P2_OUT', 'P2_IN']
TAGS = {'P2_MID': ('P2_IN', 'P2_OUT'), 'P2_OUT': ('P2_IN', 'P2_MID'), 'P2_IN': ('P2_MID', 'P2_OUT')}
Z = {'P2_IN': 320, 'P2_MID': 630, 'P2_OUT': 940}

# colours (CSS variables defined in the page style, light + dark)
CLS = [  # key, short label, css var, step
    ('no_tag_both', 'both tag chambers silent', '--k-geo', 1),
    ('no_tag_A', 'tag A silent', '--k-geo2', 1),
    ('no_tag_B', 'tag B silent', '--k-geo3', 1),
    ('tag_rejected', 'tag group rejected', '--k-rej', 2),
    ('outside_area', 'prediction outside the read-out area', '--k-fid', 3),
    ('dead_pad', 'prediction on a dead / masked pad', '--k-fid2', 3),
    ('probe_silent', 'probe SILENT', '--k-sil', 4),
    ('probe_far', 'probe fired FAR from the prediction', '--k-far', 4),
    ('probe_near_late', 'probe fired near, but OUTSIDE the tag time window', '--k-late', 4),
    ('found', 'FOUND', '--k-ok', 5),
]
LAB = {k: l for k, l, _, _ in CLS}
SHORT = {'found': 'found', 'probe_silent': 'silent (no hit at all)', 'probe_far': 'far (outside the search radius)',
         'probe_near_late': 'near, but late'}


def legend4(x0, y, step=170):
    o = []
    for i, k in enumerate(['found', 'probe_silent', 'probe_far', 'probe_near_late']):
        x = x0 + i * step
        o.append(f'<rect x="{x}" y="{y - 10}" width="12" height="12" style="fill:var({COL[k]})"/>'
                 f'<text x="{x + 17}" y="{y + 1}" class="lbl small">{esc(SHORT[k])}</text>')
    return ''.join(o)
COL = {k: c for k, _, c, _ in CLS}


def esc(s):
    return html.escape(str(s))


def pct(x, n, d=1):
    return f'{100 * x / n:.{d}f} %' if n else '-'


def counts(g, p):
    on = g['probes'][p]['on']
    return {k: on.get(k, 0) for k, *_ in CLS}


# ------------------------------------------------------------------------------------------- diagrams (schematic)
def svg_stand():
    """side view of the stand: paddles, chambers with their read-out band, three example tracks."""
    W, H = 760, 470
    zx = lambda z: 40 + (z + 140) * (H - 70) / 1420          # z (mm, down) -> y px
    xx = lambda x: 380 + x * 0.62                               # x (mm) -> px
    o = [f'<svg viewBox="0 0 {W} {H}" class="dia" role="img" aria-label="side view of the cosmic stand">']
    # paddles
    for z0, z1, lab in ((-100, -50, 'top paddle 1'), (-50, 0, 'top paddle 2'), (1250, 1300, 'bottom paddle')):
        o.append(f'<rect x="{xx(-300)}" y="{zx(z0)}" width="{600 * 0.62}" height="{zx(z1) - zx(z0) - 1}" '
                 f'class="pad-scint"/><text x="{xx(310)}" y="{(zx(z0) + zx(z1)) / 2 + 4}" class="lbl">{lab}</text>')
    # chambers: full chamber + read-out band (the part that is read out on this stand)
    for name, z in Z.items():
        y = zx(z)
        o.append(f'<rect x="{xx(-267)}" y="{y - 5}" width="{534 * 0.62}" height="10" class="ch-full"/>')
        o.append(f'<rect x="{xx(-150)}" y="{y - 5}" width="{300 * 0.62}" height="10" class="ch-ro"/>')
        o.append(f'<text x="{xx(-275)}" y="{y + 4}" class="lbl" text-anchor="end">{name} (z = {z} mm)</text>')
    # tracks: (x at z=-100, x at z=1300, css class, label)
    tracks = [(-40, 60, 'tr-ok', 'A: crosses all three read-out areas'),
              (-60, 290, 'tr-miss', 'B: crosses IN and MID, misses the OUT read-out'),
              (240, 240, 'tr-miss2', 'C: paddles fire, crosses no read-out area')]
    for i, (x0, x1, cl, lab) in enumerate(tracks):
        o.append(f'<line x1="{xx(x0)}" y1="{zx(-100)}" x2="{xx(x1)}" y2="{zx(1300)}" class="{cl}"/>')
        o.append(f'<text x="14" y="{zx(90) + 18 * i}" class="lbl {cl}-t">— {esc(lab)}</text>')
    o.append(f'<text x="{xx(0)}" y="{zx(320) - 12}" class="lbl small" text-anchor="middle">orange = read-out area '
             f'(3 of 6 connectors), grey = not read out</text>')
    o.append('</svg>')
    return ''.join(o)


def svg_timeline():
    """time around one scintillator coincidence."""
    W, H = 760, 200
    tx = lambda t: 60 + (t + 100) * 560 / 600     # -100..500 ns -> px
    o = [f'<svg viewBox="0 0 {W} {H}" class="dia" role="img" aria-label="timing windows">']
    o.append(f'<line x1="{tx(-100)}" y1="150" x2="{tx(500)}" y2="150" class="axis"/>')
    for t in range(-100, 501, 100):
        o.append(f'<line x1="{tx(t)}" y1="150" x2="{tx(t)}" y2="156" class="axis"/>'
                 f'<text x="{tx(t)}" y="172" class="lbl small" text-anchor="middle">{t}</text>')
    o.append(f'<text x="{tx(500)}" y="192" class="lbl small" text-anchor="end">ns after the scintillator</text>')
    o.append(f'<rect x="{tx(140)}" y="40" width="{tx(310) - tx(140)}" height="110" class="win-mu"/>')
    o.append(f'<text x="{(tx(140) + tx(310)) / 2}" y="34" class="lbl" text-anchor="middle">muon window 140-310 ns</text>')
    o.append(f'<line x1="{tx(0)}" y1="40" x2="{tx(0)}" y2="150" class="scint"/>'
             f'<text x="{tx(0) + 5}" y="56" class="lbl">scintillator</text>')
    # chamber hits 165-265 (5-95 %)
    for i, (t, ch) in enumerate(((198, 'IN'), (212, 'MID'), (222, 'OUT'))):
        y = 80 + 20 * i
        o.append(f'<circle cx="{tx(t)}" cy="{y}" r="5" class="hit-{ch.lower()}"/>'
                 f'<text x="{tx(t) + 9}" y="{y + 4}" class="lbl small">{ch} hit</text>')
    # tag time window: +-100 ns around the tag chamber's earliest hit
    o.append(f'<rect x="{tx(198 - 100)}" y="122" width="{tx(298) - tx(98)}" height="14" class="win-tag"/>')
    o.append(f'<text x="{tx(98) - 4}" y="133" class="lbl small" text-anchor="end">tag-and-probe window: ±100 ns '
             f'around the tag hit</text>')
    o.append(f'<text x="{tx(440)}" y="80" class="lbl small" text-anchor="middle">random-hit check:</text>'
             f'<text x="{tx(440)}" y="96" class="lbl small" text-anchor="middle">same window 4 µs later</text>')
    o.append('</svg>')
    return ''.join(o)


def pad_grid(cells, hits, w=170, n=9, title='', extra=None, circle=None, cls_hit='hit-mid'):
    """small schematic pad grid n x n; hits = list of (i, j, css); circle = (i, j, r_cells)."""
    s = w / n
    o = [f'<svg viewBox="0 0 {w} {w + 38}" class="dia mini" role="img" aria-label="{esc(title)}">']
    for i in range(n):
        for j in range(n):
            o.append(f'<rect x="{i * s + 1}" y="{j * s + 1}" width="{s - 2}" height="{s - 2}" class="pad"/>')
    for i, j, c in hits:
        o.append(f'<rect x="{i * s + 1}" y="{j * s + 1}" width="{s - 2}" height="{s - 2}" class="{c}"/>')
    if circle:
        i, j, r = circle
        o.append(f'<circle cx="{(i + .5) * s}" cy="{(j + .5) * s}" r="{r * s}" class="circ"/>')
    if extra:
        o.append(extra(s))
    t1, _, t2 = title.partition(': ')
    o.append(f'<text x="{w / 2}" y="{w + 16}" class="lbl" text-anchor="middle" font-weight="700">{esc(t1)}</text>'
             f'<text x="{w / 2}" y="{w + 31}" class="lbl small" text-anchor="middle">{esc(t2)}</text></svg>')
    return ''.join(o)


def cluster_kinds():
    a = pad_grid(9, [(4, 4, 'hit-mid'), (5, 4, 'hit-mid')], title='CLEAN: 1-6 pads within 15 mm',
                 circle=(4.5, 4, 1.3))
    b = pad_grid(9, [(2, 3, 'hit-mid'), (3, 3, 'hit-mid'), (7, 7, 'hit-x')],
                 title='REJECTED: 2nd piece > 15 mm away', circle=(2.5, 3, 1.3))
    c = pad_grid(9, [(i, j, 'hit-mid') for i in range(3, 6) for j in range(3, 6)],
                 title='REJECTED: more than 6 pads', circle=(4, 4, 1.3))
    d = pad_grid(9, [(4, 4, 'hit-mid'), (5, 4, 'hit-dead')], title='REJECTED: hit on a dead pad',
                 circle=(4.5, 4, 1.3))
    return f'<div class="row4">{a}{b}{c}{d}</div>'


def probe_kinds():
    def pred(s):
        return ''
    a = pad_grid(9, [(4, 4, 'hit-mid')], title='FOUND: a pad within the radius', circle=(4, 4, 1.25))
    b = pad_grid(9, [], title='SILENT: no hit in the chamber', circle=(4, 4, 1.25))
    c = pad_grid(9, [(8, 1, 'hit-mid')], title='FAR: hit, outside the radius', circle=(4, 4, 1.25))
    d = pad_grid(9, [(4, 4, 'hit-late')], title='LATE: near, outside tag window', circle=(4, 4, 1.25))
    return f'<div class="row4">{a}{b}{c}{d}</div>'


def svg_prediction():
    W, H = 760, 300
    zy = lambda z: 30 + (z - 320) * 220 / 620
    o = [f'<svg viewBox="0 0 {W} {H}" class="dia" role="img" aria-label="prediction from two tags">']
    for name, z in Z.items():
        o.append(f'<line x1="140" y1="{zy(z)}" x2="620" y2="{zy(z)}" class="ch-line"/>'
                 f'<text x="130" y="{zy(z) + 4}" class="lbl" text-anchor="end">{name}</text>')
    # MID probe: tags IN at x=300, OUT at x=460
    xa, xb = 300, 460
    xm = xa + (xb - xa) * (630 - 320) / (940 - 320)
    o.append(f'<line x1="{xa}" y1="{zy(320)}" x2="{xb}" y2="{zy(940)}" class="tr-pred"/>')
    o.append(f'<circle cx="{xa}" cy="{zy(320)}" r="7" class="hit-in"/><circle cx="{xb}" cy="{zy(940)}" r="7" class="hit-out"/>')
    o.append(f'<circle cx="{xm}" cy="{zy(630)}" r="16" class="circ"/><line x1="{xm - 5}" y1="{zy(630)}" x2="{xm + 5}" '
             f'y2="{zy(630)}" class="axis"/><line x1="{xm}" y1="{zy(630) - 5}" x2="{xm}" y2="{zy(630) + 5}" class="axis"/>')
    o.append(f'<text x="{xa + 14}" y="{zy(320) - 6}" class="lbl">tag A (IN): clean cluster, its centre</text>')
    o.append(f'<text x="{xb + 14}" y="{zy(940) + 4}" class="lbl">tag B (OUT): clean cluster</text>')
    o.append(f'<text x="{xm + 24}" y="{zy(630) + 4}" class="lbl">prediction in the probe (MID) ± {RAD["P2_MID"]:.0f} mm</text>')
    o.append(f'<text x="380" y="{H - 8}" class="lbl small" text-anchor="middle">straight line through the two tag '
             f'centres (+ the run\'s median offset); for OUT or IN as probe the line is extrapolated, not interpolated</text>')
    o.append('</svg>')
    return ''.join(o)


# ------------------------------------------------------------------------------------------- data figures
def svg_funnel(g, p, w=900):
    """steps: coincidences -> both tags fired -> both clean -> expected on a live pad -> found."""
    c = counts(g, p)
    n0 = g['n_scint']
    n1 = n0 - c['no_tag_both'] - c['no_tag_A'] - c['no_tag_B']
    n2 = n1 - c['tag_rejected']
    n3 = n2 - c['outside_area'] - c['dead_pad']
    n4 = c['found']
    A, B = TAGS[p]
    steps = [
        (n0, 'scintillator coincidences (the reference: a muon went through the paddles)', [
            ('no_tag_both', f'neither {A[3:]} nor {B[3:]} fired'), ('no_tag_A', f'{A[3:]} did not fire'),
            ('no_tag_B', f'{B[3:]} did not fire')]),
        (n1, f'both tag chambers ({A[3:]} and {B[3:]}) fired', [('tag_rejected', 'a tag hit group is not a clean cluster')]),
        (n2, 'two clean tags = a tagged muon with a predicted point', [
            ('outside_area', f'prediction outside {p[3:]}\'s read-out area' +
                             (f' or < {MARG[p]:.0f} mm from its edge' if MARG[p] else '')), ('dead_pad', 'prediction on a dead pad')]),
        (n3, f'muon expected on a live pad of {p[3:]} = the efficiency denominator', [
            ('probe_silent', f'{p[3:]} SILENT'), ('probe_far', f'{p[3:]} fired far away'),
            ('probe_near_late', 'near, outside the tag window')]),
        (n4, f'FOUND: {p[3:]} fired within {RAD[p]:.0f} mm, in the window', []),
    ]
    row, x0, bw = 64, 250, w - 270
    H = row * len(steps) + 10
    o = [f'<svg viewBox="0 0 {w} {H}" class="dia funnel" role="img" aria-label="funnel {p}">']
    for i, (n, lab, losses) in enumerate(steps):
        y = 8 + i * row
        wid = bw * n / n0
        o.append(f'<rect x="{x0}" y="{y}" width="{wid:.1f}" height="24" rx="3" class="fbar{" ok" if i == 4 else ""}"/>')
        o.append(f'<text x="{x0 - 8}" y="{y + 17}" class="num" text-anchor="end">{n:,}</text>')
        o.append(f'<text x="{x0 - 8}" y="{y + 33}" class="lbl small" text-anchor="end">{pct(n, n0)} of the '
                 f'coincidences</text>')
        if wid > 7.2 * len(lab) + 12:
            o.append(f'<text x="{x0 + 6}" y="{y + 17}" class="lbl inbar">{esc(lab)}</text>')
        else:
            o.append(f'<text x="{x0 + wid + 8:.1f}" y="{y + 17}" class="lbl" font-weight="600">{esc(lab)}</text>')
        xl = x0 + bw * (n - sum(c[q] for q, _ in losses)) / n0   # losses drawn at the right end of this bar
        for k, t in losses:
            ww = bw * c[k] / n0
            if ww > 0:
                o.append(f'<rect x="{xl:.1f}" y="{y + 28}" width="{ww:.1f}" height="10" style="fill:var({COL[k]})"/>')
            xl += ww
        txt = ' · '.join(f'<tspan style="fill:var({COL[k]})">■</tspan> {esc(t)}: {c[k]:,} ({pct(c[k], n)})'
                         for k, t in losses)
        if losses:
            o.append(f'<text x="{x0 + 6}" y="{y + 52}" class="lbl small">↓ lost: {txt}</text>')
    o.append('</svg>')
    return ''.join(o), dict(n0=n0, n1=n1, n2=n2, n3=n3, n4=n4, **c)


def svg_split(rows, w=900):
    """100 % bars: of the expected muons, found / silent / far / late. rows = [(label, counts)]"""
    keys = ['found', 'probe_silent', 'probe_far', 'probe_near_late']
    rh, x0 = 34, 190
    H = rh * len(rows) + 40
    o = [f'<svg viewBox="0 0 {w} {H}" class="dia" role="img" aria-label="split of expected muons">']
    for i, (lab, c) in enumerate(rows):
        n = sum(c[k] for k in keys)
        y = 6 + i * rh
        o.append(f'<text x="{x0 - 8}" y="{y + 17}" class="lbl" text-anchor="end">{esc(lab)}</text>')
        x = x0
        for k in keys:
            ww = (w - x0 - 20) * c[k] / n
            o.append(f'<rect x="{x:.1f}" y="{y}" width="{max(ww - 2, 0):.1f}" height="24" style="fill:var({COL[k]})"/>')
            if ww > 38:
                o.append(f'<text x="{x + ww / 2:.1f}" y="{y + 17}" class="lbl inbar2" text-anchor="middle">'
                         f'{100 * c[k] / n:.1f} %</text>')
            x += ww
    y = 6 + len(rows) * rh + 14
    x = x0
    for k in keys:
        o.append(f'<rect x="{x}" y="{y - 10}" width="12" height="12" style="fill:var({COL[k]})"/>'
                 f'<text x="{x + 17}" y="{y + 1}" class="lbl small">{esc(SHORT[k])}</text>')
        x += 170
    o.append('</svg>')
    return ''.join(o)


def svg_chips(p, g):
    ch = g['probes'][p]['chip']
    keys = ['found', 'probe_silent', 'probe_far', 'probe_near_late']
    chips = sorted(ch, key=int)
    w, H, x0 = 900, 250, 60
    bw = (w - x0 - 20) / len(chips)
    o = [f'<svg viewBox="0 0 {w} {H}" class="dia" role="img" aria-label="per chip {p}">']
    for f in (0, .25, .5, .75, 1):
        y = 10 + (1 - f) * 170
        o.append(f'<line x1="{x0}" y1="{y}" x2="{w - 20}" y2="{y}" class="grid"/><text x="{x0 - 6}" y="{y + 4}" '
                 f'class="lbl small" text-anchor="end">{int(f * 100)} %</text>')
    for i, c in enumerate(chips):
        d = ch[c]
        n = sum(d.get(k, 0) for k in keys)
        y = 180
        for k in keys:
            hh = 170 * d.get(k, 0) / n if n else 0
            o.append(f'<rect x="{x0 + i * bw + 4:.1f}" y="{y - hh:.1f}" width="{bw - 8:.1f}" height="{max(hh - 1, 0):.1f}" '
                     f'style="fill:var({COL[k]})"><title>chip {c}: {k} {d.get(k, 0)} of {n}</title></rect>')
            y -= hh
        o.append(f'<text x="{x0 + (i + .5) * bw:.1f}" y="196" class="lbl small" text-anchor="middle">chip {c}</text>'
                 f'<text x="{x0 + (i + .5) * bw:.1f}" y="211" class="lbl small" text-anchor="middle">n={n}</text>')
        if n:
            o.append(f'<text x="{x0 + (i + .5) * bw:.1f}" y="{180 - 170 * d.get("found", 0) / n + 14:.1f}" '
                     f'class="lbl inbar2" text-anchor="middle">{100 * d.get("found", 0) / n:.0f} %</text>')
    o.append(legend4(x0, 240))
    o.append('</svg>')
    return ''.join(o)


def svg_hist(edges, vals, xlabel, w=430, H=210, cls='bar', vals2=None, cls2='bar2', log=False, mark=None, curve=None):
    x0, y0 = 50, 170
    f = (lambda v: math.log10(v + 1)) if log else (lambda v: v)
    m = max(f(max(vals)), 1e-9)
    bw = (w - x0 - 10) / len(vals)
    o = [f'<svg viewBox="0 0 {w} {H}" class="dia" role="img" aria-label="{esc(xlabel)}">']
    for i, v in enumerate(vals):
        hh = 150 * f(v) / m
        o.append(f'<rect x="{x0 + i * bw + 1:.1f}" y="{y0 - hh:.1f}" width="{max(bw - 2, 1):.1f}" height="{hh:.1f}" '
                 f'class="{cls}"><title>{edges[i]} to {edges[i + 1]}: {v}</title></rect>')
        if vals2:
            h2 = 150 * f(vals2[i]) / m
            o.append(f'<rect x="{x0 + i * bw + 1:.1f}" y="{y0 - h2:.1f}" width="{max(bw - 2, 1):.1f}" height="{h2:.1f}" '
                     f'class="{cls2}"/>')
    for i in range(0, len(edges), max(1, len(edges) // 6)):
        o.append(f'<text x="{x0 + i * bw:.1f}" y="{y0 + 15}" class="lbl small" text-anchor="middle">{edges[i]}</text>')
    if curve:      # expected counts per bin (e.g. a fit), drawn as a line through the bin centres
        pts = ' '.join(f'{x0 + (i + 0.5) * bw:.1f},{y0 - 150 * max(f(c), 0) / m:.1f}' for i, c in enumerate(curve))
        o.append(f'<polyline points="{pts}" fill="none" class="fitl"/>')
    marks = [] if mark is None else (mark if isinstance(mark, list) else [(mark, f'{mark} mm', 'markl')])
    for j, (mv, ml, mc) in enumerate(marks):
        xm = x0 + (mv - edges[0]) / (edges[1] - edges[0]) * bw
        o.append(f'<line x1="{xm:.1f}" y1="{y0 - 155}" x2="{xm:.1f}" y2="{y0}" class="{mc}"/>'
                 f'<text x="{xm + 5:.1f}" y="{y0 - 142 + 13 * j}" class="lbl small">{esc(ml)}</text>')
    ticks = [1, 10, 100, 1000, 10000] if log else [max(vals)]
    for t in ticks:
        if f(t) <= m * 1.001:
            y = y0 - 150 * f(t) / m
            o.append(f'<text x="{x0 - 6}" y="{y + 4:.1f}" class="lbl small" text-anchor="end">{t}</text>'
                     f'<line x1="{x0 - 3}" y1="{y:.1f}" x2="{x0}" y2="{y:.1f}" class="axis"/>')
    o.append(f'<line x1="{x0}" y1="{y0}" x2="{w - 10}" y2="{y0}" class="axis"/>'
             f'<text x="{(x0 + w) / 2}" y="{H - 6}" class="lbl small" text-anchor="middle">{esc(xlabel)}</text></svg>')
    return ''.join(o)


def svg_chip_channels(name, first=24, n=6, w=900):
    """one chip's pads (light) and channels first..first+n-1 numbered, joined in channel order."""
    x0, x1, y0, y1 = BOUNDS[name]
    sc = w / (x1 - x0)
    h = (y1 - y0) * sc
    X = lambda x: (x - x0) * sc
    Y = lambda y: h - (y - y0) * sc
    from collections import Counter as _C
    chip = _C(q[2] for q in PADS[name]).most_common(1)[0][0]
    o = [f'<svg viewBox="0 0 {w} {h + 30:.0f}" class="dia" role="img" aria-label="chip map">']
    for x, y, v, c in PADS[name]:
        o.append(f'<rect x="{X(x) - 6 * sc:.1f}" y="{Y(y) - 5.8 * sc:.1f}" width="{11.4 * sc:.1f}" height="{11 * sc:.1f}" '
                 f'class="{"chip-pad" if v == chip else "pad"}"><title>VMM {v} ch {c}</title></rect>')
    sel = sorted((c, X(x), Y(y)) for x, y, v, c in PADS[name] if v == chip and first <= c < first + n)
    for (c0, xa, ya), (c1, xb, yb) in zip(sel, sel[1:]):
        o.append(f'<line x1="{xa:.1f}" y1="{ya:.1f}" x2="{xb:.1f}" y2="{yb:.1f}" class="chl"/>')
    for c, x, y in sel:
        o.append(f'<rect x="{x - 6 * sc:.1f}" y="{y - 5.8 * sc:.1f}" width="{11.4 * sc:.1f}" height="{11 * sc:.1f}" '
                 f'class="hit-mid"/><text x="{x:.1f}" y="{y - 9:.1f}" class="lbl" text-anchor="middle" '
                 f'font-weight="700">{c}</text>')
    o.append(f'<text x="8" y="{h + 22:.0f}" class="lbl small">{name}: pads of VMM {chip} (light orange), channels '
             f'{first}-{first + n - 1} (numbered), joined in channel order</text></svg>')
    return ''.join(o)


# ------------------------------------------------------------------------------------------- event displays
BOUNDS = {}
for name, pads in PADS.items():
    xs = [q[0] for q in pads]; ys = [q[1] for q in pads]
    BOUNDS[name] = (min(xs) - 10, max(xs) + 10, min(ys) - 10, max(ys) + 10)


def svg_chamber(name, hits, role, pred=None, w=250):
    x0, x1, y0, y1 = BOUNDS[name]
    sc = w / (x1 - x0)
    h = (y1 - y0) * sc
    X = lambda x: (x - x0) * sc
    Y = lambda y: h - (y - y0) * sc
    o = [f'<svg viewBox="0 0 {w} {h + 44:.0f}" class="dia evd" role="img" aria-label="{name}">']
    for x, y, v, c in PADS[name]:
        o.append(f'<rect x="{X(x) - 6 * sc:.1f}" y="{Y(y) - 5.8 * sc:.1f}" width="{11.4 * sc:.1f}" height="{11 * sc:.1f}" '
                 f'class="pad"/>')
    nin = 0
    for x, y, dt, adc, v, c, ok in hits:
        if x is None:
            continue
        inwin = 140 <= dt <= 310
        nin += inwin
        cl = ('hit-' + name[3:].lower()) if inwin and ok else ('hit-dead' if inwin else 'hit-oow')
        o.append(f'<rect x="{X(x) - 6 * sc:.1f}" y="{Y(y) - 5.8 * sc:.1f}" width="{11.4 * sc:.1f}" height="{11 * sc:.1f}" '
                 f'class="{cl}"><title>VMM {v} ch {c}: {dt} ns after the scintillator, ADC {adc}'
                 f'{"" if ok else ", dead/masked pad"}</title></rect>')
    if pred:
        o.append(f'<circle cx="{X(pred[0]):.1f}" cy="{Y(pred[1]):.1f}" r="{RAD['P2_MID'] * sc:.1f}" class="circ"/>'
                 f'<circle cx="{X(pred[0]):.1f}" cy="{Y(pred[1]):.1f}" r="1.6" class="predc"/>')
    o.append(f'<text x="{w / 2}" y="{h + 16:.0f}" class="lbl" text-anchor="middle">{name}: {esc(role)}</text>'
             f'<text x="{w / 2}" y="{h + 32:.0f}" class="lbl small" text-anchor="middle">{nin} hit(s) in the muon window'
             f'</text></svg>')
    return ''.join(o)


def event_display(p, ev):
    A, B = TAGS[p]
    roles = {A: 'tag A', B: 'tag B', p: 'probe'}
    order = ['P2_IN', 'P2_MID', 'P2_OUT']
    panes = ''.join(svg_chamber(n, ev['hits'][n], roles[n], ev['pred'] if n == p else None) for n in order)
    lines = []
    for n in order:
        hs = [h for h in ev['hits'][n] if -200 <= h[2] <= 600]
        if hs:
            lines.append(f'{n}: ' + ', '.join(f'{h[2]} ns{"" if h[6] else " (dead pad)"}' for h in hs[:12]) +
                         (' …' if len(hs) > 12 else ''))
        else:
            lines.append(f'{n}: no hit from -200 to +600 ns')
    return (f'<figure class="evfig"><div class="row3">{panes}</div><figcaption>{esc(ev["run"])} · hit times after the '
            f'scintillator: ' + ' · '.join(esc(x) for x in lines) + '</figcaption></figure>')


def d_pred(ev, p):
    """distance of the nearest in-window live probe hit to the prediction."""
    if not ev['pred']:
        return math.inf
    d = [math.hypot(h[0] - ev['pred'][0], h[1] - ev['pred'][1]) for h in ev['hits'][p]
         if h[0] is not None and h[6] and 140 <= h[2] <= 310]
    return min(d) if d else math.inf


def pick(p, cls, sub=None):
    lst = EX.get(p, {}).get(cls, [])
    if sub == 'near':
        lst = sorted([e for e in lst if d_pred(e, p) < 30], key=lambda e: d_pred(e, p)) or lst
    elif sub == 'far':
        lst = sorted([e for e in lst if d_pred(e, p) >= 60], key=lambda e: -d_pred(e, p)) or lst
    if not lst:
        return '<p class="faint">(no example stored)</p>'
    ev = lst[0]
    extra = f' · nearest {p} hit {d_pred(ev, p):.0f} mm from the prediction' if cls == 'probe_far' else ''
    return event_display(p, ev).replace('</figcaption>', esc(extra) + '</figcaption>')


PBINS = ['inside >60', 'inside 30-60', 'inside 15-30', 'inside 0-15', 'misses 0-15', 'misses 15-30',
         'misses 30-60', 'misses 60-120', 'misses >120']
PKEYS = ['found', 'probe_silent', 'probe_far', 'probe_near_late']


def svg_paddle_diagram():
    W, H = 900, 380
    zy = lambda z: 30 + (z + 120) * (H - 60) / 1420
    xx = lambda x: 450 + x * 0.9
    o = [f'<svg viewBox="0 0 {W} {H}" class="dia" role="img" aria-label="tag line extrapolated to the paddles">']
    for z0, z1, lab, x0, x1 in ((-100, 0, 'top paddles', -100, 100), (1225, 1275, 'bottom paddle', -60, 140)):
        o.append(f'<rect x="{xx(x0)}" y="{zy(z0)}" width="{xx(x1) - xx(x0)}" height="{zy(z1) - zy(z0)}" class="pad-scint"/>'
                 f'<text x="{xx(x1) + 8}" y="{(zy(z0) + zy(z1)) / 2 + 4}" class="lbl">{lab}</text>')
    for name, z in Z.items():
        o.append(f'<line x1="{xx(-330)}" y1="{zy(z)}" x2="{xx(330)}" y2="{zy(z)}" class="ch-line"/>'
                 f'<text x="{xx(-340)}" y="{zy(z) + 4}" class="lbl" text-anchor="end">{name}</text>')
    # good: IN at 0, OUT at 40 -> line inside both paddles
    def line(xa, xb, cls):
        f = lambda z: xa + (xb - xa) * (z - 320) / (940 - 320)
        return (f'<line x1="{xx(f(-100))}" y1="{zy(-100)}" x2="{xx(f(1275))}" y2="{zy(1275)}" class="{cls}"/>'
                f'<circle cx="{xx(xa)}" cy="{zy(320)}" r="6" class="hit-in"/><circle cx="{xx(xb)}" cy="{zy(940)}" r="6" '
                f'class="hit-out"/>')
    o.append(line(0, 40, 'tr-ok'))
    o.append(line(-60, 230, 'tr-miss'))
    o.append(f'<text x="{xx(325)}" y="{zy(1205)}" class="lbl tr-miss-t" text-anchor="end">misses the bottom paddle</text>')
    o.append(f'<text x="{xx(12)}" y="{zy(200)}" class="lbl tr-ok-t">goes back through both</text>')
    o.append(f'<text x="{W / 2}" y="{H - 6}" class="lbl small" text-anchor="middle">the scintillators fired, so the '
             f'muon crossed both paddles; a tag line that misses them points to a wrong tag (or to the paddle-position '
             f'uncertainty)</text></svg>')
    return ''.join(o)


def svg_paddle_bins(g, p, w=900):
    P = g['probes'][p]['paddle']
    rh, x0, xr = 30, 150, 150
    H = rh * len(PBINS) + 46
    o = [f'<svg viewBox="0 0 {w} {H}" class="dia" role="img" aria-label="classes vs paddle miss {p}">']
    for i, b in enumerate(PBINS):
        c = P.get(b, {})
        n = sum(c.get(k, 0) for k in PKEYS)
        y = 6 + i * rh
        o.append(f'<text x="{x0 - 8}" y="{y + 16}" class="lbl" text-anchor="end">{esc(b)} mm</text>')
        if not n:
            continue
        x = x0
        for k in PKEYS:
            ww = (w - x0 - xr) * c.get(k, 0) / n
            o.append(f'<rect x="{x:.1f}" y="{y}" width="{max(ww - 1.5, 0):.1f}" height="22" style="fill:var({COL[k]})">'
                     f'<title>{b}: {k} {c.get(k, 0)} of {n}</title></rect>')
            if ww > 42:
                o.append(f'<text x="{x + ww / 2:.1f}" y="{y + 16}" class="lbl inbar2" text-anchor="middle">'
                         f'{100 * c.get(k, 0) / n:.0f} %</text>')
            x += ww
        o.append(f'<text x="{w - xr + 10}" y="{y + 16}" class="lbl small">n = {n:,} (+{c.get("outside_area", 0):,} '
                 f'outside area)</text>')
    o.append(legend4(x0, H - 12, 175))
    o.append('</svg>')
    return ''.join(o)


def paddle_table(groups):
    cuts = [('all tagged muons', None), ('line inside both paddles', 4), ('misses by at most 15 mm', 5),
            ('misses by at most 30 mm', 6), ('misses by at most 60 mm', 7)]
    rows = ['<table><tr><th>selection</th>' + ''.join(f'<th>{esc(lab)}: {p[3:]}</th>' for lab, _ in groups
                                                        for p in PROBES) + '</tr>']
    for lab, k in cuts:
        cells = []
        for _, g in groups:
            for p in PROBES:
                P = g['probes'][p]['paddle']
                bins = PBINS if k is None else PBINS[:k]
                tot = Counter()
                for b in bins:
                    tot.update(P.get(b, {}))
                n = sum(tot.get(q, 0) for q in PKEYS)
                cells.append(f'<td>{pct(tot.get("found", 0), n)}<br><span class="faint">{tot.get("probe_silent", 0) / n * 100 if n else 0:.1f} % silent'
                             f', n={n:,}</span></td>')
        rows.append(f'<tr><td>{esc(lab)}</td>' + ''.join(cells) + '</tr>')
    rows.append('</table>')
    return '<div style="overflow-x:auto">' + ''.join(rows) + '</div>'


def svg_xtalk(g, w=440, H=250):
    """P(neighbour channel c+-k fires within +-100 ns | muon pad with ADC in bin), on minus off, per chamber."""
    out = []
    for d in ('P2_IN', 'P2_MID', 'P2_OUT'):
        X_ = g['xtalk'].get(d)
        if not X_:
            continue
        x0, y0, ww, hh = 50, 200, w - 70, 175
        vals = []
        for k in range(3):
            row = []
            for b in range(11):
                n = X_['n'][b]
                row.append((X_['on'][k][b] - X_['off'][k][b]) / n if n >= 30 else None)
            vals.append(row)
        top = max([v for r in vals for v in r if v is not None] + [0.05])
        top = math.ceil(top * 20) / 20
        o = [f'<svg viewBox="0 0 {w} {H}" class="dia" role="img" aria-label="crosstalk {d}">']
        for f in (0, .25, .5, .75, 1):
            y = y0 - hh * f
            o.append(f'<line x1="{x0}" y1="{y:.1f}" x2="{x0 + ww}" y2="{y:.1f}" class="grid"/><text x="{x0 - 5}" '
                     f'y="{y + 4:.1f}" class="lbl small" text-anchor="end">{100 * top * f:.0f} %</text>')
        X = lambda b: x0 + (b + .5) * ww / 11
        for k, cls in zip(range(3), ('xk1', 'xk2', 'xk3')):
            pts = [(X(b), y0 - hh * v / top) for b, v in enumerate(vals[k]) if v is not None]
            if pts:
                o.append(f'<polyline points="{" ".join(f"{x:.1f},{y:.1f}" for x, y in pts)}" class="{cls}"/>')
                o.extend(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="3.5" class="{cls}d"/>' for x, y in pts)
        for b in range(11):
            lab = 'sat' if b == 10 else f'{b * 100}'
            o.append(f'<text x="{X(b):.1f}" y="{y0 + 14}" class="lbl small" text-anchor="middle">{lab}</text>')
        o.append(f'<text x="{x0 + ww / 2}" y="{y0 + 32}" class="lbl small" text-anchor="middle">ADC of the muon pad '
                 f'(bins of 100; sat = 1000+)</text><text x="{x0}" y="14" class="lbl" font-weight="700">{d}</text>')
        o.append('</svg>')
        out.append(''.join(o))
    leg = ('<p class="small-legend"><span class="sw" style="background:var(--xk1)"></span>channel ± 1 &nbsp; '
           '<span class="sw" style="background:var(--xk2)"></span>± 2 &nbsp; '
           '<span class="sw" style="background:var(--xk3)"></span>± 3 (chance rate from +4 µs subtracted)</p>')
    return '<div class="row3">' + ''.join(out) + '</div>' + leg


def xtalk_totals(g):
    out = {}
    for d, X_ in g['xtalk'].items():
        n = sum(X_['n'])
        out[d] = [(sum(X_['on'][k]) - sum(X_['off'][k])) / n if n else 0 for k in range(3)]
        out[d + ':sat'] = [(X_['on'][k][10] - X_['off'][k][10]) / X_['n'][10] if X_['n'][10] else 0 for k in range(3)]
        out[d + ':b0'] = [(X_['on'][k][0] - X_['off'][k][0]) / X_['n'][0] if X_['n'][0] else 0 for k in range(3)]
        out[d + ':low'] = [(sum(X_['on'][k][:3]) - sum(X_['off'][k][:3])) / max(sum(X_['n'][:3]), 1) for k in range(3)]
    return out


GAIN_COL = {3.0: '--g3', 4.5: '--g45', 6.0: '--g6', 9.0: '--g9'}


def svg_gain_scan(GS, d, y='silent', w=440, H=300, ymax=0.30):
    rows = [r for r in GS['configs'] if r['det'] == d and r['thr_fC'] is not None]
    x0, y0, hh = 52, H - 55, H - 80
    ww = w - 70
    X = lambda t: x0 + (t - 2.5) / (9.0 - 2.5) * ww
    Y = lambda v: y0 - v / ymax * hh
    o = [f'<svg viewBox="0 0 {w} {H}" class="dia" role="img" aria-label="{y} vs threshold {d}">']
    for f in range(0, int(ymax * 100) + 1, 5):
        o.append(f'<line x1="{x0}" y1="{Y(f / 100):.1f}" x2="{x0 + ww}" y2="{Y(f / 100):.1f}" class="grid"/>'
                 f'<text x="{x0 - 6}" y="{Y(f / 100) + 4:.1f}" class="lbl small" text-anchor="end">{f} %</text>')
    for t in range(3, 10):
        o.append(f'<text x="{X(t):.1f}" y="{y0 + 16}" class="lbl small" text-anchor="middle">{t}</text>')
    o.append(f'<line x1="{x0}" y1="{y0}" x2="{x0 + ww}" y2="{y0}" class="axis"/>'
             f'<text x="{x0 + ww / 2}" y="{y0 + 34}" class="lbl small" text-anchor="middle">threshold above baseline '
             f'/ gain (fC, nominal)</text><text x="{x0}" y="14" class="lbl" font-weight="700">{d}: {y.upper()}</text>')
    # reachable band: below ~4.4 fC
    o.append(f'<rect x="{X(2.5):.1f}" y="{y0 - hh}" width="{X(4.45) - X(2.5):.1f}" height="{hh}" class="reach"/>'
             f'<text x="{X(2.6):.1f}" y="{y0 - hh + 14}" class="lbl small">reachable (R = 2-50 Hz)</text>')
    fit = GS['fits'].get(f'{d}:silent:all') if y == 'silent' else GS['fits'].get(f'{d}:found:no_bscan')
    if fit:
        a, b = fit['slope_per_fC'], fit['at_0fC']
        lo, hi = fit['range_fC']
        o.append(f'<line x1="{X(lo):.1f}" y1="{Y(a * lo + b):.1f}" x2="{X(hi):.1f}" y2="{Y(a * hi + b):.1f}" '
                 f'class="fitl"/>')
    for r in rows:
        v, n = r[y], r['n']
        if v is None or not n:
            continue
        e = math.sqrt(v * (1 - v) / n)
        col = GAIN_COL.get(r['gain'], '--ink2')
        bs = r['family'] == 'bscan'
        ref = r['family'] == 'REF'
        o.append(f'<line x1="{X(r["thr_fC"]):.1f}" y1="{Y(v - e):.1f}" x2="{X(r["thr_fC"]):.1f}" y2="{Y(v + e):.1f}" '
                 f'style="stroke:var({col})" stroke-width="1.5"/>')
        shape = (f'<rect x="{X(r["thr_fC"]) - 5:.1f}" y="{Y(v) - 5:.1f}" width="10" height="10"' if ref else
                 f'<circle cx="{X(r["thr_fC"]):.1f}" cy="{Y(v):.1f}" r="5"')
        style = (f'style="fill:var(--bg-elev);stroke:var({col});stroke-width:2"' if bs else
                 f'style="fill:var({col});stroke:var(--bg-elev);stroke-width:1.5"')
        o.append(f'{shape} {style}><title>{r["label"]}: {r["thr_fC"]:.2f} fC, {y} {100 * v:.1f} ± {100 * e:.1f} % '
                 f'(n = {n})</title>{"</rect>" if ref else "</circle>"}')
    o.append('</svg>')
    return ''.join(o)


def gain_legend():
    items = ''.join(f'<span class="sw" style="background:var({c})"></span>gain {g:g} mV/fC &nbsp; '
                    for g, c in GAIN_COL.items())
    return (f'<p class="small-legend">{items} · square = REF · open circle = bscan (gain 4.5, threshold scan) · '
            'dashed = straight line through all points · error bars = statistics</p>')


def gain_section(GS):
    out = []
    A = out.append
    F = GS['fits']
    rows = sorted([r for r in GS['configs'] if r['det'] == 'P2_MID'], key=lambda r: r['thr_fC'])
    A('<h2>7. Does the electronics gain change the efficiency?</h2>')
    A('<p>The weekend grid sets the VMM gain to 3, 4.5, 6 and 9 mV/fC. At each gain every chip\'s global threshold '
      'and every channel\'s trim were set (<code>bscan_target.py</code>, from the HV-on noise scans of 18 Nov) so that '
      '<b>each channel sits where its own measured noise rate is R = 2, 10 or 50 Hz</b>: optimal for equal noise per '
      'channel, not chosen for efficiency. Channels whose trim would have to go beyond its range are clipped and stay '
      'higher: 0-24 at gain 3, 28-167 at gain 6, 360-708 at gain 9 (of 1152). "bscan" is the same at gain 4.5, REF '
      'the reference equalisation. For each '
      'setting the threshold is converted to charge: (threshold - baseline) in mV divided by the nominal gain. '
      'This is a nominal charge (the transfer function is not measured yet). The quantity plotted is the one that '
      'does not depend on positions: <b>SILENT</b>, the share of muons (two clean tags, prediction on a live pad) for '
      'which the probe records no hit at all.</p>')
    A(f'<figure>{svg_gain_scan(GS, "P2_MID", w=900, H=380)}</figure>')
    A('<div class="row2">' + svg_gain_scan(GS, 'P2_IN') + svg_gain_scan(GS, 'P2_OUT') + '</div>')
    A(gain_legend())
    A('<div style="overflow-x:auto"><table><tr><th>P2_MID, setting</th><th>runs</th><th>threshold (fC)</th>'
      '<th>muons</th><th>silent</th><th>found</th></tr>' +
      ''.join(f'<tr><td>{esc(r["label"])}</td><td>{r["runs"]}</td><td>{r["thr_fC"]:.1f}</td><td>{r["n"]:,}</td>'
              f'<td>{100 * r["silent"]:.1f} %</td><td>{100 * r["found"]:.1f} %'
              f'{" (positions suspect)" if r["family"] == "bscan" and r["far"] > 0.15 else ""}</td></tr>'
              for r in rows) + '</table></div>')
    fm, fl = F['P2_MID:silent:all'], F['P2_MID:silent:below5fC']
    fo = F['P2_OUT:silent:below5fC']
    A('<p><b>What it says.</b></p><ol>'
      '<li><b>The electronics gain itself does not matter, the threshold in charge does.</b> Gains 3, 4.5 and 6 '
      'put the threshold at 3.2-4.4 fC (the noise-rate rule gives about the same charge at every gain) and give the '
      'same silent share, 12.4-14.5 % for P2_MID. Gain 9 sits at 6.7-8 fC and loses clearly more (22-25 %). Two '
      'things push its thresholds up: the same channel needs ~1.5x more charge for the same noise rate, and on '
      '360-708 of the 1152 channels the trims could not reach their target (clipped). Section 8 shows that the '
      'clipped channels are not more often silent, so clipping is not what makes gain 9 worse.</li>'
      f'<li><b>Over the whole range</b> the silent share grows by {100 * fm["slope_per_fC"]:.1f} ± '
      f'{100 * fm["slope_err"]:.1f} points per fC (P2_MID; P2_OUT and P2_IN similar).</li>'
      f'<li><b>Inside the reachable range</b> (3.2-4.4 fC, the shaded band) the slope is only {100 * fl["slope_per_fC"]:.1f} ± '
      f'{100 * fl["slope_err"]:.1f} points per fC for P2_MID ({100 * fo["slope_per_fC"]:.1f} ± {100 * fo["slope_err"]:.1f} '
      'for P2_OUT). At gains 3 and 4.5, going from R = 2 Hz to 50 Hz (25x more noise) lowers the threshold by about '
      '0.5 fC and wins at most 1.4 points. The lowest settings we have, 3.2-3.5 fC, still leave 12.5-13 % of the '
      'muons silent in P2_MID.</li>'
      '<li><b>So with these chambers at this HV the VMM settings reach about 87 % "any hit" (P2_MID) and about 82 % '
      'with our position and time cuts; the electronics settings cannot close the gap to Dream</b> (93-94 % on the same '
      'pads). Lowering the threshold further is not possible without a large noise rate.</li>'
      '<li><b>How much detector gain is needed is not determined by this scan.</b> If the silent share depends only '
      'on charge / threshold, a slope of ~3 points per fC at a 4 fC threshold means each +25 % of detector gain '
      '(equivalent to -0.8 fC) wins ~2-3 points, and removing the ~7 points to the Dream level would need the '
      'equivalent of a 1.5 fC threshold, i.e. ~2.5x more detector gain. With the flatter slope seen at low threshold '
      'it would need more, and the straight-line intercepts (3-9 % silent at 0 fC, depending on the fit) leave room '
      'for a part of the 12.5 % that is not threshold at all. The direct measurement is a mesh-HV scan with REF '
      'electronics: silent vs mesh HV.</li></ol>')
    A('<p class="faint">Caveats: the threshold in fC uses the nominal gains and the DAC-vs-baseline distance, not a '
      'measured transfer function (the muon pulse height confirms the nominal gains, section 8). The '
      'bscan settings show up to 45 % FAR (probe fired, but far from the prediction) although random hits are below '
      '1 %: the tag positions are wrong there for a reason still to find; their silent share is normal and is used.</p>')
    return '\n'.join(out)


# ------------------------------------------------------------------------------------------- the page
def frac4(on):
    n = sum(on.get(k, 0) for k in PKEYS)
    return n, {k: on.get(k, 0) / n if n else float('nan') for k in PKEYS}


def pointing_section(exp):
    """3b: the per-probe search radius from the pointing resolution, and the before / after table."""
    o = ['<h3>3b. The search radius follows the pointing resolution (updated 20 Nov)</h3>']
    o.append('<p>The predicted point is the straight line through the two tag clusters. For P2_MID it is an '
             '<b>interpolation</b> (tags IN and OUT, half-way); for P2_IN and P2_OUT it is an <b>extrapolation</b> by one '
             'full gap beyond the nearer tag. With z = 320 / 630 / 940 mm the prediction is p = w<sub>A</sub> a + '
             'w<sub>B</sub> b and its error, per axis, σ<sub>point</sub><sup>2</sup> = w<sub>A</sub><sup>2</sup> '
             'σ<sub>A</sub><sup>2</sup> + w<sub>B</sub><sup>2</sup> σ<sub>B</sub><sup>2</sup>; the probe hit adds its '
             'own σ<sub>P</sub>. With the same resolution s in every chamber:</p>')
    o.append('<table><tr><th>probe</th><th>tags</th><th>prediction</th><th>σ<sub>point</sub></th>'
             '<th>σ<sub>res</sub> (hit − prediction)</th><th>measured σ<sub>res</sub> (REF)</th>'
             '<th>15 mm keeps</th><th>radius now (99 %)</th><th>edge margin (2 σ<sub>point</sub>)</th></tr>')
    for p in ('P2_IN', 'P2_MID', 'P2_OUT'):
        md = PT['model'][p]
        wa, wb = md['w']
        ms = PT['measured'].get('REF', {}).get(p, {})
        sg = ms.get('sigma_res')
        keep = f'{100 * (1 - math.exp(-15 ** 2 / (2 * sg ** 2))):.0f} %' if sg else '-'
        fs = math.hypot(wa, wb)
        o.append(f'<tr><td>{p}</td><td>{TAGS[p][0][3:]} + {TAGS[p][1][3:]}</td><td>{wa:+.1f} a {wb:+.1f} b</td>'
                 f'<td>{fs:.2f} s</td><td>{math.hypot(fs, 1):.2f} s</td>'
                 f'<td><b>{sg:.1f} ± {ms.get("err", 0):.1f} mm</b></td><td>{keep}</td>'
                 f'<td><b>{RAD[p]:.0f} mm</b></td><td>{MARG[p]:.0f} mm</td></tr>')
    o.append('</table>')
    mid = PT['measured'].get('REF', {}).get('P2_MID', {}).get('sigma_res', float('nan'))
    o.append('<p>So with one fixed 15 mm radius the extrapolated chambers lose about a fifth of their real hits to the '
             'tail of the pointing, while P2_MID loses almost none: that is the "FAR" excess of IN and OUT, a '
             '<b>calculation</b> effect, not the electronics. The widths are <b>measured</b> on this stand: the radial '
             'distance of the nearest probe hit (section 5a) is fitted with a 2-D Gaussian plus a flat part for unrelated '
             f'hits (<code>pointing_resolution.py</code>). P2_MID gives σ<sub>res</sub> = {mid:.1f} mm, i.e. '
             f's ≈ {mid / math.sqrt(1.5):.1f} mm per chamber - the 12 mm pads alone give 11.8/√12 = 3.4 mm, and the SPS '
             'beam measured a pad-limited σ<sub>core</sub> ≈ 3.65 mm. Because the stand is equally spaced, the model '
             'predicts σ<sub>res</sub>(IN, OUT) = 2 σ<sub>res</sub>(MID) for <i>any</i> mix of chamber resolutions; '
             'the data give 1.8. Multiple scattering (~1 mrad over 310 mm) is negligible.</p>')
    o.append('<p><b>New rule:</b> a probe hit is accepted within R<sub>P</sub> = 3.03 σ<sub>res,P</sub> (99 % of a '
             '2-D Gaussian), and the predicted point must lie at least 2 σ<sub>point,P</sub> inside the chamber '
             'outline (otherwise the muon may really have crossed outside the read-out area, and would count as '
             'silent). The margin muons are counted with "outside the area" in the funnel.</p>')
    if OLD:
        o.append('<table><tr><th>probe (REF)</th><th>radius / margin</th><th>tested</th><th>found</th><th>silent</th>'
                 '<th>far</th><th>late</th><th>removed by the margin (would have been: found / silent / far)</th></tr>')
        for p in ('P2_IN', 'P2_MID', 'P2_OUT'):
            for lab, G_, r_, m_ in (('before', OLD, 15.0, 0.0), ('now', G, RAD[p], MARG[p])):
                P = G_['REF']['probes'][p]
                n, f = frac4(P['on'])
                mc = P.get('margin_cut', {}) if lab == 'now' else {}
                cut = (f'{mc.get("found", 0):,} / {mc.get("probe_silent", 0):,} / {mc.get("probe_far", 0):,}'
                       if lab == 'now' else '')
                o.append(f'<tr><td>{p if lab == "before" else ""}</td><td>{lab}: {r_:.0f} mm / {m_:.0f} mm</td>'
                         f'<td>{n:,}</td><td><b>{100 * f["found"]:.1f} %</b></td>'
                         f'<td>{100 * f["probe_silent"]:.1f} %</td><td>{100 * f["probe_far"]:.1f} %</td>'
                         f'<td>{100 * f["probe_near_late"]:.1f} %</td><td>{cut}</td></tr>')
        o.append('</table>')
        ex = []
        for gk, lab in (('g3.0 gscan_hvon10_hvm', 'gain 3, R = 10 Hz'), ('g9.0 gscan_hvon10_hvm', 'gain 9, R = 10 Hz')):
            if gk in OLD and gk in G:
                ex.append(lab + ': ' + ', '.join(
                    f'{p[3:]} {100 * frac4(OLD[gk]["probes"][p]["on"])[1]["found"]:.1f} → '
                    f'{100 * frac4(G[gk]["probes"][p]["on"])[1]["found"]:.1f} %' for p in ('P2_IN', 'P2_OUT', 'P2_MID')))
        o.append('<p class="faint">Chance hits (the same test 4 µs later) stay at ~0.01 % with the larger radius. The '
                 'improvement is the same in every gain-scan setting (' + '; '.join(ex) + ').</p>')
    return ''.join(o)


def pillar_section(exp):
    """3c: the bulk pillars set the ceiling below 100 %."""
    o = ['<h3>3c. The ceiling is not 100 %: the bulk pillars</h3>']
    o.append('<p>The mesh of these chambers stands on insulating <b>pillars</b>; a muon whose charge lands on a pillar '
             'is not amplified. All three chambers of the stand were bulked at Saclay (det1-4) with insulation mask '
             '<b>V1</b> (<code>Detector_Drawings/Version_Apr26/Insulation_masks/V1/P2_BASKET-Mask_M2_V1.gbr</code>): '
             '12,477 pillars of ⌀ 0.8 mm on a 4 mm grid plus 5 big ⌀ 6.15 mm ones, <b>3.78 % of the pad area</b> - the '
             'same mask and number the cosmic-bench (DREAM) analysis uses (<code>p2_qa_config.MASK_GBR_PATH</code>). '
             '(The CERN-bulked det5 has a different mask, 0.5 mm on 2 mm = 4.87 %, and is not on this stand.) On det1 '
             'that bench measured a "no hit" share of 3.17 % against its Gerber pillar area of 3.43 %, so its efficiency '
             'is quoted against a ceiling of ~96.6 %, not 100 %. The SPS beam (DREAM, telescope at 0.2-0.3 mm) showed '
             'that the cost grows with the threshold - a pillar also removes part of the charge around it: 2.4-3.9 '
             'points per chamber at the high-gain point, 9.5 for P2_IN at the nominal one - so at the VMM threshold '
             'the area is a floor, not the whole pillar cost.</p>')
    o.append('<p>Here the pads (12 mm) and the pointing (3-8 mm) cannot resolve a 3 mm lattice, so the pillar loss '
             f'cannot be cut out muon by muon; it is a floor under the silent share (the five big ones can be, see 5d). Taking '
             'the mask\'s'
             f' {100 * PILLAR:.2f} % as the ceiling:</p>')
    o.append('<table><tr><th>probe (REF)</th><th>found</th><th>found / (1 − pillars)</th><th>silent</th>'
             '<th>silent beyond the pillars</th></tr>')
    for p in ('P2_IN', 'P2_MID', 'P2_OUT'):
        e = exp[p]
        fo, si = e['found'] / e['n3'], e['probe_silent'] / e['n3']
        o.append(f'<tr><td>{p}</td><td>{100 * fo:.1f} %</td><td><b>{100 * fo / (1 - PILLAR):.1f} %</b></td>'
                 f'<td>{100 * si:.1f} %</td><td>{100 * (si - PILLAR):.1f} points</td></tr>')
    o.append('</table>')
    o.append('<p>P2_IN and P2_OUT are then at ~89-91 % of what the pillars allow; P2_MID is ~7 points further down, '
             'all of it in the "no hit at all" class and spread evenly over its chips and pads (sections 5b, 5d). '
             'That is the part left to explain.</p>')
    return ''.join(o)


def silent_tests_section():
    """5d: angle, pad map, tag pulse height of the silent muons (silent_tests.py on the per-tag dumps)."""
    R = ST['REF']
    o = ['<h3>5d. Why silent? Angle, place, tag pulse height</h3>']
    o.append('<p>From the per-tag dumps of the same runs (scintillator-confirmed tags, random tag pairs subtracted, the '
             'same edge margin; silent = no live probe hit within ±100 ns of the tag + latency). First: <b>a silent '
             'muon really leaves nothing</b> - for ' + ', '.join(
                 f'{p[3:]} {100 * R[p]["dead1us"][0]:.1f} of the {100 * R[p]["silent"][0]:.1f} %'
                 for p in ('P2_IN', 'P2_MID', 'P2_OUT')) +
             ' the whole chamber has no hit on any chip within ±1 µs.</p>')
    hd = ''.join(f'<th>{b["lo"]}-{b["hi"]}</th>' for b in R['P2_MID']['vs_slope'][:5])
    o.append('<table><tr><th>silent vs tan θ (REF)</th>' + hd + '</tr>')
    for p in ('P2_IN', 'P2_MID', 'P2_OUT'):
        o.append(f'<tr><td>{p}</td>' + ''.join(
            f'<td>{100 * b["f"]:.1f} ± {100 * b["err"]:.1f} %<br><span class="faint">n {b["n"]:.0f}</span></td>'
            for b in R[p]['vs_slope'][:5]) + '</tr>')
    o.append('</table>')
    o.append('<p><b>1. Angle: flat.</b> A hole in the amplification (pillar, dead area) swallows mostly vertical '
             'tracks: with the 4 mm drift gap an inclined track spreads its charge over 4 tan θ mm of mesh, more than '
             'a 0.8 mm pillar from tan θ ≈ 0.2. The silent share drops by only ~2 points from the vertical bin to '
             'tan θ = 0.3 (about what the longer path, 1/cos θ, gives to a signal near threshold), so the silent muons '
             'are <b>not</b> a geometric hole.</p>')
    med = {}
    for p in ('P2_IN', 'P2_MID', 'P2_OUT'):
        fr = sorted(r[5] for r in R[p]['pads'])
        med[p] = fr[len(fr) // 2] if fr else float('nan')
    o.append('<p><b>2. Place: MID is high everywhere.</b> Per predicted pad (REF, pads with ≥ 15 tags) the '
             'median pad is silent for ' + ', '.join(f'{p[3:]} {100 * med[p]:.1f} %' for p in ('P2_IN', 'P2_MID',
                                                                                              'P2_OUT')) +
             ': IN and OUT sit near the pillar floor on a typical pad and reach their mean through a few bad pads; '
             'P2_MID\'s typical pad is already at 12 %. The pad-to-pad spread is only modestly above statistical (χ²/ndf ' +
             ', '.join(f'{p[3:]} {R[p]["pad_chi2_ndf"][0]:.0f}/{R[p]["pad_chi2_ndf"][1]}' for p in ('P2_IN', 'P2_MID',
                                                                                                    'P2_OUT')) +
             '). Two local effects exist but are small: predictions within 12 mm of a dead / masked pad are silent '
             f'{100 * R["P2_MID"]["vs_dead_dist"][0]["f"]:.0f} % (MID), ~0.3 points of the total; and the big ⌀ 6.15 mm '
             f'pillars are visible in P2_MID (within 8 mm: {100 * R["P2_MID"]["vs_big_pillar"][0]["f"]:.0f} % silent, '
             f'n = {R["P2_MID"]["vs_big_pillar"][0]["n"]:.0f}; all configurations: '
             f'{100 * ST["ALL"]["P2_MID"]["vs_big_pillar"][0]["f"]:.0f} %), worth &lt; 0.2 points. The worst pads '
             'repeat across chambers on the same channel numbers (e.g. channel 62 at (418, 139) mm on MID and OUT, '
             'channel 42 at (330, 115) mm): to follow up with more statistics.</p>')
    o.append('<p><b>3. Tag pulse height: no fake tags.</b> If some tags were noise in time with the scintillator, '
             'low-ADC tags would be followed by more silent probes. At REF the silent share is flat in the smaller '
             'tag\'s largest ADC (MID: ' + ', '.join(f'{b["lo"]}+: {100 * b["f"]:.0f} %' for b in
                                                     R['P2_MID']['vs_tag_adc'][:5]) +
             '). Over all configurations it even rises with tag ADC, because the high-gain settings have both larger '
             'pulses and more silent muons.</p>')
    o.append('<p class="faint">Conclusion of 5d: the P2_MID silent muons cross live, normal pads at normal angles '
             'with genuine tags, and the chamber records nothing at all within a microsecond. What is left is a signal '
             'below threshold over the whole chamber (lower gas gain, or a higher effective threshold, of P2_MID) - '
             'the mesh-HV scan, or DREAM on the same chamber, decides. Code: <code>silent_tests.py</code>.</p>')
    return ''.join(o)


def cuts_section():
    """3d: every cut of the current estimate, in order."""
    rows = [
        ('reference', 'scintillator coincidence (top AND bottom paddles, VMM 0 ch 40)', 'a muon crossed the stand; '
         'pulser configurations (_tpA / _tpB) excluded by name; first / last 10 µs of a file dropped'),
        ('chamber "fired"', 'any hit 140-310 ns after the scintillator', 'hits grouped in time with a 100 ns gap'),
        ('tag', 'both other chambers have a CLEAN group in the window: 1-6 hits, all on live pads, all within 15 mm of '
         'their mean, not inside a stack burst (&gt; 36 hits in one group of the whole stack)',
         'position = mean pad centre, time = earliest hit; otherwise "tag rejected" (not counted)'),
        ('prediction', 'straight line through the two tag positions (z 320 / 630 / 940 mm) + the run\'s median '
         'probe offset', ''),
        ('fiducial', 'nearest probe pad within 9 mm of the prediction and live (≥ 1 hit in the run, not on a dead '
         'chip, not vetoed); <b>new:</b> prediction at least 2 σ<sub>point</sub> inside the outline: ' +
         ', '.join(f'{p[3:]} {MARG[p]:.0f} mm' for p in ('P2_IN', 'P2_MID', 'P2_OUT')),
         'otherwise "outside the area" / "dead pad" (not counted)'),
        ('silent', 'no hit at all in the probe chamber 140-310 ns after the scintillator', 'counted as a miss'),
        ('far', '<b>new:</b> nearest live probe hit farther than R<sub>P</sub> = 3.03 σ<sub>res</sub>: ' +
         ', '.join(f'{p[3:]} {RAD[p]:.0f} mm' for p in ('P2_IN', 'P2_MID', 'P2_OUT')) + ' (was 15 mm for all)',
         'counted as a miss'),
        ('late', 'a live hit within R<sub>P</sub> in 140-310 ns, but not within ±100 ns of tag time + latency',
         'counted as a miss'),
        ('found', 'a live hit within R<sub>P</sub> and within ±100 ns of tag time + latency', ''),
        ('chance hits', 'every probe class re-evaluated with the window 4 µs later (same tags)', '~0.01 % found by '
         'chance: negligible'),
        ('efficiency', 'found / (found + silent + far + late)', f'ceiling: 1 − pillars ≈ {100 * (1 - PILLAR):.1f} % '
         '(section 3c)'),
    ]
    o = ['<h3>3d. The cuts of the current estimate, in order</h3>',
         '<table><tr><th>step</th><th>cut</th><th>note</th></tr>']
    o += [f'<tr><td><b>{a}</b></td><td>{b}</td><td>{c}</td></tr>' for a, b, c in rows]
    o.append('</table><p class="faint">Code: <code>muon_funnel.py --pointing pointing.json</code> (cuts and parameters '
             'from the <code>efficiency_tp.py</code> json of each chunk; <code>efficiency_tp.py --pointing</code> applies '
             'the same radius and margin to the self-tagged estimate). The efficiency definition itself is still on '
             'hold.</p>')
    return ''.join(o)


def body():
    fm, nm = svg_funnel(REF, 'P2_MID')
    fo, no = svg_funnel(REF, 'P2_OUT')
    fi, ni = svg_funnel(REF, 'P2_IN')
    exp = {p: n for p, n in (('P2_MID', nm), ('P2_OUT', no), ('P2_IN', ni))}
    rej = REF['rejects']
    rj = {d: (rej[d]['states'].get('rejected', 0), sum(rej[d]['states'].values())) for d in rej}

    def eff(n):
        return n['found'] / n['n3']
    off = {p: REF['probes'][p]['off'] for p in PROBES}
    out = []
    A = out.append
    A('<div class="callout"><b>In one paragraph.</b> We start from every <b>scintillator coincidence</b> (a muon went '
      'through the paddles) and follow it step by step through what the tag-and-probe needs: the two other chambers '
      'must fire (<i>tags</i>), their hits must look like one clean muon cluster, the line through them must point at a '
      'live pad of the chamber under test (<i>probe</i>), and the probe must fire near that point. Every coincidence '
      'ends in exactly one box. The steps before the probe decide <b>which muons we test</b>; only the last step '
      f'decides the <b>efficiency</b>. For P2_MID in the reference configuration: of {nm["n3"]:,} muons expected on a '
      f'live pad, {pct(nm["found"], nm["n3"])} are found, <b>{pct(nm["probe_silent"], nm["n3"])} leave no hit at all '
      f'in the chamber</b>, {pct(nm["probe_far"], nm["n3"])} fire far from the prediction, '
      f'{pct(nm["probe_near_late"], nm["n3"])} fire near but outside the time window.</div>')
    if PT:
        no_, ni_ = exp['P2_OUT'], exp['P2_IN']
        A('<div class="callout warn"><b>Update 20 Nov: two calculation effects corrected.</b> (1) The search radius now '
          'follows the <b>pointing resolution</b> of each probe (section 3b): the extrapolated P2_IN and P2_OUT '
          f'predictions are ~2x less precise than the interpolated P2_MID one, so they get {RAD["P2_IN"]:.0f} mm instead '
          f'of 15 mm, and a {MARG["P2_IN"]:.0f} mm edge margin. P2_IN and P2_OUT go to {pct(ni_["found"], ni_["n3"])} '
          f'and {pct(no_["found"], no_["n3"])} found (from 72-74 %); P2_MID does not move. (2) The ceiling is not '
          f'100 %: the bulk pillars (mask V1) cover {100 * PILLAR:.2f} % of the pad area (section 3c). The cuts in use are listed in '
          'section 3d. All numbers on this page use the new radius, except the PNG figures of section 8 (made with the '
          'fixed 15 mm; P2_MID, which they show, changes by &lt; 0.3 points).</div>')
    A('<p class="faint">Data: lab cosmic stand, weekend scan from 19 Nov 2026, runs '
      f'{min(ALL["runs"])}-{max(ALL["runs"])} without the pulser runs (configurations ending in _tpA / _tpB). '
      f'Reference configuration REF (gain 4.5 mV/fC, trimmed thresholds Km18; chip config file '
      f'<code>..._trimmed_Km18_v1_hvm.txt</code> in the logbook): runs {", ".join(map(str, sorted(set(REF["runs"]))))} '
      f'(268 = the 20-min end-to-end test), {REF["n_scint"]:,} scintillator coincidences; all configurations: '
      f'{ALL["n_scint"]:,}. Code: <code>cosmic_vmm_analysis/muon_funnel.py</code> (same cuts as '
      '<code>efficiency_tp.py</code>).</p>')

    A('<h2>1. The set-up and the words</h2>')
    A('<h3>1a. The stand</h3>')
    A(f'<figure>{svg_stand()}<figcaption>Side view, not to scale. Three P2 chambers between two stacked top paddles and '
      'one bottom paddle. A <b>scintillator coincidence</b> = top AND bottom paddles fired: a muon crossed the stand. '
      'Only 3 of the 6 read-out connectors of each chamber are connected here, so a chamber can only see a muon that '
      'crosses its <b>read-out area</b>. The 4 Nov trigger simulation found that only ~44 % of the coincidences cross '
      'all three read-out areas (track A); the others are like B or C. These are not inefficiency: nothing could '
      'have been read out.</figcaption></figure>')
    A('<h3>1b. Time</h3>')
    A(f'<figure>{svg_timeline()}<figcaption>All three chambers fire 165-265 ns after the scintillator (5-95 % of the muons). '
      'Everything on this page uses the <b>muon window 140-310 ns</b> after the scintillator to decide whether a chamber '
      '"fired". The tag-and-probe code itself (<code>efficiency_tp.py</code>) does not use the scintillator: it looks '
      'for the probe hit within <b>±100 ns of the tag hit</b>. To count chance hits (noise), every test is repeated '
      'with the window moved 4 µs later, where there is no muon.</figcaption></figure>')
    A('<h3>1c. Tag, clean cluster, rejected group</h3>')
    A('<p>A chamber\'s hits are grouped in time: hits closer than 100 ns to the previous one belong to the same '
      '<b>group</b>. A group is accepted as a <b>tag</b> only if it looks like one muon: 1 to 6 pads, all within 15 mm '
      'of their centre, all on working pads, and not part of a stack-wide burst. Otherwise the whole group is '
      '<b>rejected</b> and the muon is not used at all - it is neither a success nor a failure.</p>')
    A(cluster_kinds())
    A('<h3>1d. Prediction and probe</h3>')
    A(f'<figure>{svg_prediction()}</figure>')
    A(probe_kinds())

    A('<h2>2. The funnel: where every coincidence ends</h2>')
    A('<p>Read top to bottom. Each bar is what is left; the coloured strip under it is what was removed at that step, '
      'with the reason. The first three steps only decide <b>which muons enter the test</b>. The last step is the '
      'efficiency.</p>')
    A('<h3>P2_MID as probe (tags IN and OUT) - the clean case</h3>')
    A(f'<figure>{fm}</figure>')
    A('<details><summary>P2_OUT as probe (tags IN and MID)</summary>' + f'<figure>{fo}</figure></details>')
    A('<details><summary>P2_IN as probe (tags MID and OUT)</summary>' + f'<figure>{fi}</figure></details>')
    A('<p>P2_MID sits between its two tags, so its predicted point is an interpolation; for IN and OUT it is an '
      'extrapolation and more predictions fall outside the read-out area or near its edge.</p>')

    A('<h2>3. The efficiency step, alone</h2>')
    A('<p>Only the muons that reached the last step: expected on a live pad of the probe. This is where the '
      'inefficiency is.</p>')
    A(svg_split([(f'{p} (REF)', exp[p]) for p in PROBES] +
                [(f'{p} (all configs)', counts(ALL, p)) for p in PROBES]))
    A('<table><tr><th>probe</th><th>expected on a live pad</th><th>found</th><th>silent</th><th>far</th><th>late</th>'
      '<th>chance hits (window +4 µs)</th></tr>')
    for p in PROBES:
        n = exp[p]
        o_ = off[p]
        n_off = sum(o_.values())
        rnd = (o_.get('found', 0) + o_.get('probe_far', 0) + o_.get('probe_near_late', 0))
        A(f'<tr><td>{p}</td><td>{n["n3"]:,}</td><td><b>{pct(n["found"], n["n3"])}</b></td>'
          f'<td>{pct(n["probe_silent"], n["n3"])}</td><td>{pct(n["probe_far"], n["n3"])}</td>'
          f'<td>{pct(n["probe_near_late"], n["n3"])}</td><td>{pct(rnd, n_off)} of the tests have any probe hit</td></tr>')
    A('</table>')
    if PT:
        A(pointing_section(exp))
    A(pillar_section(exp))
    A(cuts_section())

    A('<h2>4. Each class, with a real event</h2>')
    A('<p>Each display is the three chambers seen from above; only the read-out area is drawn (the odd shape is the '
      'P2 pad fan, connectors 3-5). Light squares = pads. <b>Coloured</b> = a hit in the muon window (blue IN, orange '
      'MID, green OUT); <b>dark grey</b> = a hit on a dead/masked pad; <b>outlined</b> = a hit outside the muon window. '
      f'The dashed circle = the {RAD["P2_MID"]:.0f} mm search radius around the prediction in the probe (here P2_MID). Hover a pad for '
      'its chip, channel, time and ADC. All events are from the reference configuration.</p>')
    sections = [
        ('found', None, 'FOUND - the normal case', 'Both tags are one clean cluster each, the line through them points '
         'into P2_MID, and P2_MID has a hit inside the circle.'),
        ('probe_silent', None, 'SILENT probe - the chamber recorded nothing',
         'Two clean tags, the prediction is on a live pad, and P2_MID has <b>no hit at all</b> - no channel, anywhere in '
         'the chamber - in the muon window. No choice of radius, window or clustering can recover these. Read-out '
         'losses, dead time and sparks (19 Nov tests) explain at most 0.55 points of them together. The rest: the muon\'s '
         'signal stayed below the VMM threshold on every pad it crossed (what the beam autopsy found for the same '
         'chambers).'),
        ('probe_far', 'near', 'FAR (a): just outside the circle',
         'P2_MID fired 15-30 mm from the prediction. This is the muon, found by the chamber, but the prediction was off '
         'by a bit more than the radius (tag position resolution, multiple scattering, alignment). It is a '
         '<b>calculation</b> loss: a larger radius would count it.'),
        ('probe_far', 'far', 'FAR (b): a hit somewhere else entirely',
         'P2_MID fired, but far from the prediction. Most likely the muon pad was silent (as in SILENT) and an unrelated '
         'hit - noise or crosstalk - happened in the same 170 ns. This is an <b>electronics</b> loss that hides as a '
         'position mismatch.'),
        ('probe_near_late', None, 'LATE - near, but outside the tag-and-probe window',
         'P2_MID fired inside the circle and 140-310 ns after the scintillator, but more than 100 ns from the tag hit '
         'that the tag-and-probe uses as the clock. Rare (0.6 %).'),
        ('tag_rejected', None, 'Tag group REJECTED - the muon leaves the sample',
         'P2_IN (tag A) has the muon plus single pads far away, all within 100 ns: the group is "too wide", so this muon '
         'is not used. Section 5c shows what these extra pads are.'),
        ('outside_area', None, 'Prediction OUTSIDE the read-out area',
         'Both tags are clean, but the line points more than 9 mm beyond the nearest pad of the probe - here just past '
         'the edge of the read-out area. Correctly removed.'),
        ('no_tag_A', None, 'A tag chamber did not fire',
         'P2_OUT and P2_MID fired, P2_IN did not: the muon missed P2_IN\'s read-out area (geometry, like track B in 1a), '
         'or P2_IN was inefficient. This pair of tags cannot test P2_MID for this muon.'),
    ]
    for cls, sub, title, txt in sections:
        A(f'<h3><span class="sw" style="background:var({COL[cls]})"></span>{esc(title)}</h3><p>{txt}</p>')
        A(pick('P2_MID', cls, sub))

    A('<h2>5. Looking closer</h2>')
    A('<h3>5a. FAR: how far is the nearest probe hit?</h3>')
    edges = list(range(0, 155, 5))
    HG = OLD['REF'] if OLD else REF        # the distance is measured without a radius cut; OLD has no edge margin
    for p in PROBES:
        df, dfd = HG['probes'][p]['d_far'], HG['probes'][p]['d_found']
        vals = [df.get(str(e), 0) + dfd.get(str(e), 0) for e in edges[:-1]]
        beyond = sum(v for k, v in df.items() if k != 'inf' and int(k) >= 150)
        curve, mk, cap = None, 15, ''
        if PT and 'REF' in PT.get('measured', {}) and p in PT['measured']['REF']:
            ms = PT['measured']['REF'][p]
            sg, ns = ms['sigma_res'], ms['n_sig']
            curve = [ns * (math.exp(-lo ** 2 / (2 * sg ** 2)) - math.exp(-(lo + 5) ** 2 / (2 * sg ** 2)))
                     for lo in edges[:-1]]
            mk = [(15, 'old 15 mm', 'markl'), (RAD[p], f'new {RAD[p]:.0f} mm', 'markn')]
            cap = (f' Dashed line: the fitted pointing spread, σ<sub>res</sub> = {sg:.1f} mm (2-D Gaussian; the flat '
                   f'part is unrelated hits). Red dashed: the old 15 mm radius, which keeps '
                   f'{100 * (1 - math.exp(-15 ** 2 / (2 * sg ** 2))):.0f} % of the real muon hits; green: the new '
                   f'{RAD[p]:.0f} mm (99 %).')
        if p == 'P2_MID':
            near_mid = sum(vals[3:6])
        A(f'<figure>{svg_hist(edges, vals, f"{p}: distance of the nearest probe hit to the prediction (mm), REF - log scale", w=900, log=True, mark=mk, curve=curve)}'
          f'<figcaption>{p} as probe ({"interpolated" if p == "P2_MID" else "extrapolated"}). Muons expected on a live '
          f'pad whose probe fired in the muon window: {sum(vals[3:6]):,} at 15-30 mm, {sum(vals[6:]) + beyond:,} '
          f'further out.{cap}</figcaption></figure>')
    vals = [0, 0, 0, near_mid, 0, 0]
    A('<h3>5b. Per chip: is the loss everywhere or in a few chips?</h3>')
    def chip_rng(p, k='probe_silent'):
        fr = [v.get(k, 0) / sum(v.get(q, 0) for q in PKEYS) for v in REF['probes'][p]['chip'].values()
              if sum(v.get(q, 0) for q in PKEYS) >= 50]
        return f'{100 * min(fr):.0f}-{100 * max(fr):.0f} %' if fr else '-'
    A('<p>The chip is the one that reads the predicted pad. If the loss were a few broken chips it would show here. '
      f'P2_MID: every chip is silent for {chip_rng("P2_MID")} of its muons - a chamber-wide loss. P2_IN '
      f'({chip_rng("P2_IN")}) and P2_OUT ({chip_rng("P2_OUT")}) are not uniform: their best chips are silent for only '
      '3-5 %, about what the bulk pillars alone cost (section 3c).</p>')
    for p in PROBES:
        A(f'<figure>{svg_chips(p, REF)}<figcaption>{p} as probe (REF).</figcaption></figure>')
    if ST:
        A(silent_tests_section())
    A('<h3>5c. What are the extra pads that get tag groups rejected?</h3>')
    rows = []
    for d in ('P2_IN', 'P2_MID', 'P2_OUT'):
        E = rej[d]['extra']
        r, n = rj[d]
        act = n - rej[d]['states'].get('none', 0)
        one = E['size'].get('1', 0)
        near_ch = E['dch'][0] + E['dch'][1]
        rows.append(f'<tr><td>{d}</td><td>{pct(r, act)}</td><td>{E["n"]:,}</td><td>{pct(one, E["n"], 0)}</td>'
                    f'<td>{pct(E["same_chip"], E["n"], 0)}</td><td>{pct(near_ch, E["same_chip"], 0)}</td></tr>')
    A('<table><tr><th>chamber</th><th>muon windows rejected</th><th>extra pieces</th><th>single pad</th>'
      '<th>on the same chip as the muon</th><th>... and 1-3 channels away from a muon channel</th></tr>' +
      ''.join(rows) + '</table>')
    E = rej['P2_MID']['extra']
    A('<div class="row2">' + svg_hist(E['dch_edges'], E['dch'], 'P2_MID, same chip: channel number difference to the '
                                      'muon (bins of 2)') +
      svg_hist(E['sep_edges'], E['sep'], 'P2_MID: distance on the chamber to the muon (mm)', vals2=E['sep_same']) +
      '</div>')
    A('<div class="row2">' + svg_hist(E['dt_edges'], E['dt'], 'P2_MID: time of the extra pad minus the muon (ns)') +
      svg_hist(E['adc_edges'], E['adc'], 'P2_MID: ADC of the extra pad') + '</div>')
    A('<p>The extra pads are <b>one pad</b>, <b>small</b> (ADC 50-100, just above threshold), <b>at the same time</b> '
      'as the muon, and when on the same chip <b>1-3 channel numbers</b> away from a channel the muon fired - yet '
      '20-400 mm away on the chamber. By chance a random pad would be on the muon\'s chip ~17 % of the time (1 chip '
      'of 6) and 1-3 channels away ~10 % of the time; we see ~60 % and ~65 %. That is the signature of <b>electronic '
      'crosstalk between neighbouring VMM channels</b> (on the chip, the hybrid or the cable). It always looks like a '
      'far pad because, on the P2 board, consecutive channels of a chip are never neighbouring pads: they are '
      '<b>149 mm apart (median)</b>, every one of them more than 15 mm:</p>')
    A(f'<figure>{svg_chip_channels("P2_MID")}<figcaption>One chip of P2_MID: the pads it reads (a strip along the band) '
      'and six consecutive channels numbered. Crosstalk from channel 26 into 25 or 27 lights a pad far along the '
      'strip, so the tag group becomes "too wide" and the muon is dropped.</figcaption></figure>')
    A('<p>So for these muons the right cluster is the large one, and the extra pad can be recognised without using the '
      'probe: single pad, low ADC, in time, 1-3 channels from a muon channel on the same chip. Keeping them would '
      'add ~15 % more tagged muons. It would not change the efficiency much: the lost muons fire the probe about as '
      'often as the kept ones (+1-3 points).</p>')

    xt = xtalk_totals(REF)
    A('<h3>5d. How often does crosstalk fire in muon data?</h3>')
    A('<p>For every pad of a muon cluster (in any chamber, in the muon window): does the channel next to it on the same '
      'chip (± 1, ± 2, ± 3) fire within ± 100 ns? The same channels 4 µs later give the chance rate (it is zero to '
      '1 in 10 000). Plotted against the ADC of the muon pad:</p>')
    A(svg_xtalk(REF))
    A('<p>For a muon pad with a medium signal (ADC 200-700) a neighbour fires in about 1 % of the cases. When the '
      f'muon pad <b>saturates</b> (ADC 1000+) the ± 1 neighbour fires in {100 * xt["P2_IN:sat"][0]:.0f} % (IN), '
      f'{100 * xt["P2_MID:sat"][0]:.0f} % (MID), {100 * xt["P2_OUT:sat"][0]:.0f} % (OUT) of the cases, and the '
      f'<b>± 2 neighbour even more</b>: {100 * xt["P2_IN:sat"][1]:.0f} %, {100 * xt["P2_MID:sat"][1]:.0f} %, '
      f'{100 * xt["P2_OUT:sat"][1]:.0f} %; ± 3 stays at the per-cent level. The size agrees with the internal-pulser '
      'study (16 Nov: the neighbour\'s baseline moves by ~3.5 % of the pulse): a saturated muon pad carries a charge '
      'well above the ADC range, and a few per cent of it is above the threshold. So this is the crosstalk seen with '
      'the pulser, now with real signals; here it costs tags, not probe efficiency. Two things to understand: why '
      '± 2 couples more than ± 1 (the channel order on the hybrid / connector?), and the lowest-ADC bin (muon pad ADC '
      f'&lt; 100: ± 1 fires in {100 * xt["P2_OUT:b0"][0]:.0f}-{100 * xt["P2_IN:b0"][0]:.0f} %), where the "muon pad" '
      'is itself small and may not be the muon.</p>')

    A('<h3>5e. Do the tagged muons go back to the scintillators?</h3>')
    A(f'<figure>{svg_paddle_diagram()}</figure>')
    A('<p>The scintillators fired, so the muon went through the top and the bottom paddles. We extend the line '
      'through the two tag clusters up to the top paddles (z = -100 ... 0 mm) and down to the bottom paddle '
      '(z = 1225 ... 1275 mm) and measure how far it is inside or outside them. The paddle positions are '
      '<b>not measured</b>: they come from the 4 Nov fit of the October tag maps (statistical precision ~1-2 mm, '
      'systematic unknown; the stand may have been touched since). The extrapolation also magnifies the tag position '
      'errors (the bottom paddle is 335 mm below OUT). So a miss of a few cm is not yet proof of a wrong tag; '
      'the test is whether the classes change with the miss.</p>')
    A(f'<figure>{svg_paddle_bins(REF, "P2_MID")}<figcaption>P2_MID as probe (REF): the muons expected on a live pad, '
      'binned by how far the tag line is inside the paddles (top) or misses them (bottom).</figcaption></figure>')
    A('<details><summary>P2_OUT and P2_IN as probe</summary>' +
      f'<figure>{svg_paddle_bins(REF, "P2_OUT")}<figcaption>P2_OUT (REF)</figcaption></figure>' +
      f'<figure>{svg_paddle_bins(REF, "P2_IN")}<figcaption>P2_IN (REF)</figcaption></figure></details>')
    A('<p><b>Found share</b> (and the silent share below it) when only the tag lines that go back to the paddles are '
      'kept:</p>')
    A(paddle_table([('REF', REF), ('all', ALL)]))
    PM = REF['probes']['P2_MID']['paddle']
    def share(bins, key):
        t = Counter()
        for b in bins:
            t.update(PM.get(b, {}))
        n = sum(t.get(q, 0) for q in PKEYS)
        return 100 * t.get(key, 0) / n if n else 0, n
    fi, ni = share(PBINS[:4], 'found')
    fo, no_ = share(PBINS[4:], 'found')
    si, _ = share(PBINS[:4], 'probe_silent')
    so, _ = share(PBINS[4:], 'probe_silent')
    A(f'<p>For P2_MID: {no_:,} of {ni + no_:,} tagged muons ({100 * no_ / (ni + no_):.0f} %) have a line that misses a '
      f'paddle. They are found less often ({fo:.0f} % against {fi:.0f} %), mostly because more of them are FAR, '
      f'which is what a wrong tag would give. But the SILENT share barely changes ({so:.0f} % against {si:.0f} %), '
      'and keeping only the lines that go back through both paddles moves the efficiency by about one point. '
      '<b>Discarding these muons is a reasonable cleaning step, but it does not change the picture: the chamber '
      'stays silent for ~12-13 % of the muons that certainly crossed it.</b></p>')

    A('<h2>6. So where does the inefficiency arise?</h2>')
    m = exp['P2_MID']
    n3 = m['n3']

    def near_of(p):      # far hits 15-30 mm from the prediction, with the radius in use
        d = REF['probes'][p]['d_far']
        return sum(d.get(str(b), 0) for b in (15, 20, 25) if b >= RAD[p] - 5)
    near = near_of('P2_MID')
    A('<ol>'
      f'<li><b>Before the test (not inefficiency).</b> Of {m["n0"]:,} coincidences only {pct(m["n3"], m["n0"])} become a '
      'test of P2_MID. Most are lost because a tag chamber did not fire - mostly geometry: the muon did not cross its '
      'read-out area. Another ~15 % of the tagged muons are dropped for crosstalk pads in a tag. None of these count '
      'as a miss.</li>'
      f'<li><b>In the test, P2_MID: {pct(n3 - m["found"], n3)} of the expected muons are not found:</b><ul>'
      f'<li>{pct(m["probe_silent"], n3)} the chamber records nothing at all;</li>'
      f'<li>{pct(m["probe_far"] - near, n3)} only an unrelated hit far away (&gt; 30 mm) - the same loss, hidden;</li>'
      f'<li>{pct(near, n3)} the muon fired 15-30 mm away - beyond 3 σ of the pointing, so a non-Gaussian tail '
      '(wrong tag pad, crosstalk piece): a calculation loss;</li>'
      f'<li>{pct(m["probe_near_late"], n3)} late with respect to the tag hit - a calculation loss (clock).</li></ul></li>'
      f'<li><b>Of the {pct(m["probe_silent"], n3)} silent, ~{100 * PILLAR:.0f} points are expected from the bulk '
      'pillars</b> (section 3c; DREAM at the SPS lost 2.4-3.9 points to them at high gain, more at low gain). That '
      f'leaves ~{100 * m["probe_silent"] / n3 - 100 * PILLAR:.0f} points of P2_MID '
      'with no hit at all, not explained by the counting (read-out, dead time, sparks &lt; 0.6 points), the time '
      'window (&lt; 1.5 points, below) or the pointing.</li>'
      '<li><b>P2_IN and P2_OUT</b>, now with the radius scaled to their pointing: ' +
      '; '.join(f'{p} silent {pct(exp[p]["probe_silent"], exp[p]["n3"])}, far {pct(exp[p]["probe_far"], exp[p]["n3"])}'
                for p in ('P2_IN', 'P2_OUT')) +
      '. Their silent share is within a few points of the pillar floor; P2_MID\'s is not.</li></ol>')
    A('<p class="faint">Time window check (<code>time_window_check.py</code>, REF runs 268-269): for tagged muons the '
      'first probe hit lies in 140-310 ns for 85.6 % (MID), 71.8 % (IN), 70.7 % (OUT; no fiducial cut here) against '
      '86.9 / 73.4 / 71.7 % in -200 to +1500 ns: widening the window recovers at most 1.5 points in any chamber, and '
      'P2_MID\'s hit times (1-99 %: 150-450 ns) are as tight as the others\'.</p>')
    A('<div class="callout warn"><b>Corrections to what I said before.</b><ul>'
      '<li>I said a probe window built from the scintillator time would recover 1-4 %. With the random hits removed it '
      f'is {pct(m["probe_near_late"], n3)} (LATE). The 1-4 % came from widening the window on all tags without the '
      'random-hit subtraction.</li>'
      '<li>I proposed ranking clusters by time, paddle geometry and angle. The extra pieces turn out to be crosstalk '
      'pads, so a simpler rule (drop single, low-ADC, in-time pads 1-3 channels from a muon channel) should do; to be '
      'tested.</li></ul></div>')
    if GS:
        A(gain_section(GS))
    A('<h2>8. Checking that nothing is mixed up: four tests</h2>')
    A('<p>Before concluding that more electronics gain cannot help, four checks (code in '
      '<code>cosmic_vmm_analysis/</code>, products in <code>$W/analysis/gain_scan/out/</code>, see References). The '
      'figures are the ones of the weekly-meeting deck of 20 Nov.</p>')
    F = FIGDIR
    A(f'<div class="row2"><figure><img src="{F}c_gain_real.png" alt="gain real"><figcaption><b>1. The gain is real.</b> '
      'The muon pulse height (largest pad ADC minus ~40 ADC baseline, median) grows like the nominal gain: 1.8-2.2x at '
      'gain 6, 2.8-3.2x at gain 9, in all three chambers. The settings were applied and the conversion to fC is '
      'right.</figcaption></figure>'
      f'<figure><img src="{F}c_noise_fC.png" alt="noise in charge"><figcaption><b>2. Effect 1: more gain, more noise in '
      'charge.</b> For the same channel and the same noise rate, the threshold in charge is 1.19x (gain 6) and 1.47x '
      '(gain 9) the one at gain 3 (bscan_target reports of 18 Nov). The noise grows faster than the signal, so the '
      'threshold relative to the muon signal gets worse with gain.</figcaption></figure></div>')
    A(f'<div class="row2"><figure><img src="{F}c_trim_spread.png" alt="trim spread"><figcaption><b>3a. Effect 2: the trim '
      'range.</b> Inside one chip the needed thresholds spread over ~10 mV (gain 3), ~20 mV (gain 6) and 38-65 mV '
      '(gain 9); a trim corrects at most ~34 mV, so at gain 9 360-708 channels stay higher (clipped).</figcaption>'
      f'</figure><figure><img src="{F}c_clipped.png" alt="clipped"><figcaption><b>3b. But clipping does not cost '
      'muons.</b> Muons predicted on clipped channels are not more often silent than on well-trimmed ones (gain 9: '
      '19-23 % vs 20-24 %). Channels without a measured noise crossing (given the chip median) are worse at gain 6 '
      'R = 2 Hz and gain 9 R = 10 Hz: a calibration gap.</figcaption></figure></div>')
    A(f'<div class="row2"><figure><img src="{F}c_spectrum.png" alt="spectrum"><figcaption><b>4a. Is the silent share a '
      'threshold tail?</b> Largest-pad ADC of the found muons (full resolution, ADC non-linearity corrected), relative '
      'to the median. The threshold edge (dotted) is at 12-19 % of the typical signal; extending the spectrum below '
      'it gives 2-6 % of the muons, against 13-14 % silent. Charge sharing is not the cause: 93-97 % of the found '
      'muons fire one pad.</figcaption></figure>'
      f'<figure><img src="{F}c_late.png" alt="late"><figcaption><b>4b. At gain 9 some "silent" muons are late.</b> A hit '
      'near the prediction 100 ns - 1 µs after the window: ~3 % of the silent muons at REF, 14 % (P2_MID) and 25 % '
      '(P2_IN) at gain 9; the found muons are on time at every gain.</figcaption></figure></div>')
    A(f'<figure><img src="{F}c_silent_parts.png" alt="parts" style="max-width:720px"><figcaption><b>Putting it '
      'together (P2_MID).</b> Orange: the threshold tail predicted from the spectrum (error bar: the two '
      'extrapolations); navy: hits outside the time window; red: no hit at all, not explained by a smooth tail. '
      'Gains 3 to 6: 13-14 % silent, of which 3-5 points look like a threshold tail and 8-10 points do not. Gain 9: '
      '24 %, of which ~4 points are late hits.</figcaption></figure>')
    A('<div class="callout"><b>Conclusion.</b> Nothing is mixed up: the gain settings are real and were applied. '
      'Raising the VMM gain does not help, because the noise in charge grows faster than the signal (effect 1); the '
      'short trim range at gain 9 (effect 2) exists but does not cost muons. Within the noise we can accept, P2_MID '
      'reaches ~87-88 % "any hit". Only a part (3-5 points) of the 13 % silent muons behaves like signal just below '
      'threshold, which more detector gain would recover; the larger part (8-10 points) leaves no signal above '
      'threshold at all, and is the next thing to understand (where, which chips, does DREAM lose the same '
      'muons?).</div>')
    A('<h2>9. Open points</h2><ol>'
      '<li><b>Efficiency definition</b>: on hold until the points below are understood.</li>'
      '<li><b>P2_MID\'s extra silent share</b> (~8 points beyond the pillars, even over its chips, not the time window, '
      'not the pointing): compare with DREAM on the same chamber and pads; a chamber swap (MID ↔ OUT position) would '
      'tell a chamber property from a position-in-the-stack effect.</li>'
      '<li><b>Pillar term for the VMM</b>: the SPS fold showed the pillar cost grows with the threshold (2.4 → 9.5 '
      'points); the VMM threshold sits higher than DREAM\'s, so the 3.8 % of mask V1 may be a lower bound here.</li>'
      '<li><b>Detector gain</b>: a mesh-HV scan of P2_MID with REF electronics would measure directly how fast the '
      'silent share drops with signal (your call on the HV).</li>'
      '<li><b>bscan settings</b>: why the tag positions go wrong there (45 % FAR at R = 10 Hz).</li>'
      '<li><b>Paddle positions</b>: the paddle check uses positions fitted to the October tag maps. Has the stand '
      '(paddles or chambers) moved since? A refit on the weekend data, or a measurement with a ruler, would turn '
      'the paddle check from "indicative" into a cut.</li>'
      '<li><b>Crosstalk</b>: measured here per muon pad vs ADC; next is a filter that drops single, in-time '
      'neighbour-channel pads before the cluster cuts, and to count how many tags it recovers.</li>'
      '<li><b>Crosstalk pattern</b>: why channel ± 2 fires more often than ± 1 next to a saturated pad, and why '
      'the smallest "muon pads" (ADC &lt; 100) also show neighbours: to look at.</li>'
      '<li><b>The no-signal muons</b> (8-10 points of P2_MID): where on the chamber, which chips, at which times; and '
      'does DREAM lose the same muons on the same pads?</li>'
      '<li><b>Gain-9 late hits</b>: what delays them (large or saturated pulses?).</li>'
      '<li><b>Calibration gap</b>: channels without a measured noise-rate crossing (given the chip median).</li></ol>')
    A('<h2>10. References</h2>')
    A('<h3>Data</h3><ul>'
      '<li>Raw captures (EOS backup of the DAQ box): <code>/eos/experiment/ntof/data/x17/p2_sps_july/vmm_cosmic_run_oct26/runs/run_NNN</code></li>'
      '<li>Logbook (configuration of every run): <code>/eos/experiment/ntof/data/x17/p2_sps_july/vmm_cosmic_run_oct26/logbook/LOGBOOK.md</code> '
      '(DAQ box: <code>/local/p2/p2data/cosmics_3det_202610/logbook/</code>)</li>'
      '<li>Weekend analysis area <code>$W = /eos/experiment/ntof/data/x17/p2_sps_july/vmm_cosmic_run_oct26/analysis_lxplus/weekend_20261119</code>: '
      'per-tag dumps <code>$W/tags/</code>, pseudo-run map <code>$W/pseudo/chunks.json</code>, efficiency pipeline '
      '<code>$W/results/current/</code></li>'
      '<li>Muon funnel: <code>$W/analysis/muon_funnel/</code> (run*.json, summary.json, paddles_fit_20261104.json)</li>'
      '<li>Gain scan: <code>$W/analysis/gain_scan/out/</code> (gain_scan_eff.*, muon_adc_gain.csv, threshold_spectrum.json, '
      'silent_vs_clipped.txt, silent_late_hits.txt, silent_vs_pad_position.txt), threshold reports '
      '<code>$W/analysis/gain_scan/reports/</code> (copies of <code>~/trim_work/out/gain_20261118/hvon_g*_R*_report.json</code> and '
      '<code>~/trim_work/out/bscan_20261116/hvon*_v1_report.json</code> on the DAQ box), deck figures <code>$W/analysis/gain_scan/deck_figs/</code></li>'
      '<li>Read-out / dead time / sparks: <code>$W/analysis/deadtime_sparks/results/</code>, <code>$W/analysis/pulser_tick225/</code>, '
      'figures <code>$W/analysis/wrapup_20261119/out/</code></li>'
      '<li>Paddle positions: 4 Nov trigger-acceptance fit, <code>data/LAB_Measurements/VMM-cosmics/vmm_threshold_scan/analysis_wk2/trigger_acceptance_mc/fit.json</code></li>'
      '<li>DREAM references: cosmic_bench_analysis stage 06 <code>ray_hit_miss_list*.csv</code> (M3 bench), '
      '<code>sps_beam_analysis/dream_vmm_comparison/</code> (SPS beam)</li></ul>')
    A('<h3>Code</h3><p>Repository P2_basket_analysis, branch <code>decoder-tick-225</code>, directory <code>cosmic_vmm_analysis/</code>; '
      'all run on lxplus with <code>$W/pipeline/env.sh</code> and <code>PYTHONPATH=$W/pipeline/cva</code>.</p>'
      '<div style="overflow-x:auto"><table><tr><th>script</th><th>what</th></tr>'
      '<tr><td><code>muon_funnel.py</code> (+ <code>_summary.py</code>, <code>_job.sh</code>, <code>.sub</code>)</td><td>every scintillator coincidence through the tag-and-probe; paddle check; crosstalk rate</td></tr>'
      '<tr><td><code>efficiency_tp.py</code>, <code>channel_efficiency.py</code></td><td>the tag-and-probe and the per-tag dumps</td></tr>'
      '<tr><td><code>gain_scan_eff.py</code></td><td>funnel classes per gain-scan setting vs threshold in charge</td></tr>'
      '<tr><td><code>muon_adc_gain.py</code></td><td>real gain from the muon pulse height</td></tr>'
      '<tr><td><code>trim_reports.py</code></td><td>what the per-gain threshold calibration achieved</td></tr>'
      '<tr><td><code>threshold_spectrum.py</code></td><td>muon spectrum at the threshold edge, charge sharing</td></tr>'
      '<tr><td><code>silent_vs_clipped.py</code>, <code>silent_late_hits.py</code>, <code>silent_vs_pad_position.py</code></td><td>the silent-muon tests</td></tr>'
      '<tr><td><code>tag_rejects.py</code>, <code>wrapup_losses.py</code></td><td>rejected tags; read-out / dead time / sparks figures</td></tr>'
      '<tr><td><code>note/muon_funnel_page.py</code>, <code>note/gain_deck_figs.py</code>, <code>note/deck_update_20261120.py</code></td><td>this page, the deck figures, the deck update</td></tr>'
      '</table></div>')
    A('<h3>Documents</h3><ul>'
      '<li>Weekly-meeting deck: <code>LaCie …/1.General_Info/presentations/WeeklyMeetings/P2_Detector_Performance_Gain_2026-11-20.pptx</code> '
      '(and .pdf; build inputs in <code>.build_P2_Detector_Performance/update_20261120/</code>)</li>'
      '<li>Handoffs: <code>cosmic_vmm_analysis/HANDOFF_analysis_20261119d.md</code>, <code>HANDOFF_analysis_20261120.md</code>; '
      'wrap-up <code>WRAPUP_20261119.md</code></li>'
      '<li>Weekend run plan: <code>DAQ_Control_VMM_Beam/docs/HANDOFF_WEEKEND_RUNS_20261119.md</code>; threshold calibration: '
      '<code>vmm_monitor_tool/fastscan/bscan_target.py</code></li>'
      '<li>Related notes on this site: "Cosmic stand: 3D trigger geometry" (4 Nov), "VMM thresholds explained" (16 Nov), '
      '"Does a lower VMM threshold make the P2 chambers more efficient?" (17 Nov), "VMM efficiency autopsy" (SPS, 16 Aug)</li>'
      '</ul>')
    return '\n'.join(out)


STYLE = """
  :root { --k-geo:#b9b7b0; --k-geo2:#cfcdc6; --k-geo3:#a3a19a; --k-rej:#6250d6; --k-fid:#e87ba4; --k-fid2:#c25b87;
          --k-sil:#e34948; --k-far:#eda100; --k-late:#2a78d6; --k-ok:#1baf7a;
          --c-in:#2a78d6; --c-mid:#eb6834; --c-out:#1baf7a; --pad:#ecebe6; --pad-line:#d6d4cc; --ink:#1d1f24;
          --ink2:#5f6570; --dead:#55534e; --xk1:#e34948; --xk2:#eda100; --xk3:#6250d6; --g3:#2a78d6; --g45:#eb6834; --g6:#1baf7a; --g9:#6250d6; }
  @media (prefers-color-scheme: dark) {
    :root { --k-geo:#6f6e68; --k-geo2:#5b5a55; --k-geo3:#84837c; --k-rej:#9085e9; --k-fid:#d55181; --k-fid2:#e87ba4;
            --k-sil:#e66767; --k-far:#c98500; --k-late:#3987e5; --k-ok:#199e70;
            --c-in:#3987e5; --c-mid:#d95926; --c-out:#199e70; --pad:#2a2d33; --pad-line:#3a3e45; --ink:#e8e8e4;
            --ink2:#a5a9b1; --dead:#8c8a84; --xk1:#e66767; --xk2:#c98500; --xk3:#9085e9; --g3:#3987e5; --g45:#d95926; --g6:#199e70; --g9:#9085e9; }
  }
  main.wrap, .site-header .wrap, .site-footer .wrap { max-width: 1180px; }
  .note-body { max-width: none; }
  .note-body p, .note-body ul, .note-body .callout, .note-body details > p { max-width: 90ch; }
  svg.dia { width: 100%; height: auto; display: block; font-family: inherit; }
  svg.mini { max-width: 230px; }
  .row4 { display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 1rem; margin: 1rem 0; }
  .row3 { display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: .8rem; }
  .row2 { display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 1rem; margin: 1rem 0; }
  @media (max-width: 760px) { .row4, .row2 { grid-template-columns: repeat(2, minmax(0, 1fr)); }
                              .row3 { grid-template-columns: 1fr; } }
  .lbl { fill: var(--ink); font-size: 13px; } .small { font-size: 11px; fill: var(--ink2); }
  .num { fill: var(--ink); font-size: 15px; font-weight: 700; }
  .inbar { fill: #fff; font-size: 12.5px; font-weight: 600; paint-order: stroke; stroke: rgba(0,0,0,.35); stroke-width: 2px; }
  .inbar2 { fill: #fff; font-size: 11.5px; font-weight: 700; }
  .axis { stroke: var(--ink2); stroke-width: 1.2; } .grid { stroke: var(--pad-line); stroke-width: 1; }
  .pad { fill: var(--pad); stroke: var(--pad-line); stroke-width: .5; }
  .pad-scint { fill: #c9b8e8; stroke: #8f7cc0; } .ch-full { fill: var(--pad); stroke: var(--ink2); stroke-width: .6; }
  .ch-ro { fill: var(--c-mid); opacity: .85; } .ch-line { stroke: var(--ink2); stroke-width: 6; opacity: .35; }
  .tr-ok { stroke: var(--k-ok); stroke-width: 2.5; } .tr-ok-t { fill: var(--k-ok); }
  .tr-miss { stroke: var(--k-far); stroke-width: 2.5; stroke-dasharray: 7 4; } .tr-miss-t { fill: var(--k-far); }
  .tr-miss2 { stroke: var(--k-geo3); stroke-width: 2.5; stroke-dasharray: 3 4; } .tr-miss2-t { fill: var(--ink2); }
  .tr-pred { stroke: var(--ink2); stroke-width: 1.6; stroke-dasharray: 6 4; }
  .win-mu { fill: var(--k-ok); opacity: .14; } .win-tag { fill: var(--k-late); opacity: .3; }
  .scint { stroke: #8f7cc0; stroke-width: 3; }
  .hit-in { fill: var(--c-in); } .hit-mid { fill: var(--c-mid); } .hit-out { fill: var(--c-out); }
  .hit-x { fill: var(--k-rej); } .hit-late { fill: var(--k-late); }
  .hit-dead { fill: var(--dead); } .hit-oow { fill: none; stroke: var(--ink2); stroke-width: 1; }
  .circ { fill: none; stroke: var(--ink); stroke-width: 1.4; stroke-dasharray: 4 3; } .predc { fill: var(--ink); }
  .fbar { fill: var(--ink2); opacity: .55; } .fbar.ok { fill: var(--k-ok); opacity: 1; }
  .bar { fill: var(--c-in); } .bar2 { fill: var(--c-mid); }
  .xk1 { fill: none; stroke: var(--xk1); stroke-width: 2; } .xk1d { fill: var(--xk1); }
  .xk2 { fill: none; stroke: var(--xk2); stroke-width: 2; } .xk2d { fill: var(--xk2); }
  .xk3 { fill: none; stroke: var(--xk3); stroke-width: 2; } .xk3d { fill: var(--xk3); }
  .small-legend { font-size: .85rem; color: var(--fg-muted); }
  .reach { fill: var(--k-ok); opacity: .08; } .fitl { stroke: var(--ink2); stroke-width: 1.5; stroke-dasharray: 6 4; }
  .markl { stroke: var(--k-sil); stroke-width: 1.5; stroke-dasharray: 4 3; }
  .markn { stroke: var(--k-ok); stroke-width: 2; }
  .chl { stroke: var(--ink); stroke-width: 1.6; }
  .chip-pad { fill: var(--c-mid); opacity: .28; }
  .sw { display: inline-block; width: .8em; height: .8em; border-radius: 2px; margin-right: .45em; vertical-align: -.05em; }
  .evfig figcaption { font-size: .82rem; }
  .note-body figure img { width: 100%; height: auto; display: block; }
"""

page = body()
doc = open(sys.argv[2]).read()
doc = re.sub(r'<style>.*?</style>', '<style>' + STYLE.replace('\\', '\\\\') + '</style>', doc, count=1, flags=re.S)
i = doc.index('<div class="note-body">')
j = doc.index('</article>')
k = doc.rindex('</div>', i, j)
doc = doc[:i] + '<div class="note-body">\n' + page + '\n' + doc[k:]
open(sys.argv[2], 'w').write(doc)
print('written', sys.argv[2], len(doc))
