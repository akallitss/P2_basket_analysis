#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build the note body "Where do the missed muons go?" from muon_funnel_summary.py output (20 Nov 2026).

Self-contained HTML (inline SVG, no scripts from outside). Replaces the <div class="note-body"> of NOTE_HTML.
Usage: muon_funnel_page.py SUMMARY.json NOTE_HTML
"""
import html
import json
import math
import re
import sys

S = json.load(open(sys.argv[1]))
G = S['groups']
REF, ALL = G['REF'], G['ALL']
PADS = S['pads']
EX = S['examples']
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
SHORT = {'found': 'found', 'probe_silent': 'silent (no hit at all)', 'probe_far': 'far (> 15 mm)',
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
    a = pad_grid(9, [(4, 4, 'hit-mid')], title='FOUND: a pad within 15 mm', circle=(4, 4, 1.25))
    b = pad_grid(9, [], title='SILENT: no hit in the chamber', circle=(4, 4, 1.25))
    c = pad_grid(9, [(8, 1, 'hit-mid')], title='FAR: hit, but > 15 mm away', circle=(4, 4, 1.25))
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
    o.append(f'<text x="{xm + 24}" y="{zy(630) + 4}" class="lbl">prediction in the probe (MID) ± 15 mm</text>')
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
            ('outside_area', f'prediction outside {p[3:]}\'s read-out area'), ('dead_pad', 'prediction on a dead pad')]),
        (n3, f'muon expected on a live pad of {p[3:]} = the efficiency denominator', [
            ('probe_silent', f'{p[3:]} SILENT'), ('probe_far', f'{p[3:]} fired far away'),
            ('probe_near_late', 'near, outside the tag window')]),
        (n4, f'FOUND: {p[3:]} fired within 15 mm, in the window', []),
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


def svg_hist(edges, vals, xlabel, w=430, H=210, cls='bar', vals2=None, cls2='bar2', log=False, mark=None):
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
    if mark is not None:
        xm = x0 + (mark - edges[0]) / (edges[1] - edges[0]) * bw
        o.append(f'<line x1="{xm:.1f}" y1="{y0 - 155}" x2="{xm:.1f}" y2="{y0}" class="markl"/>'
                 f'<text x="{xm + 5:.1f}" y="{y0 - 142}" class="lbl small">{mark} mm</text>')
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
        o.append(f'<circle cx="{X(pred[0]):.1f}" cy="{Y(pred[1]):.1f}" r="{15 * sc:.1f}" class="circ"/>'
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


# ------------------------------------------------------------------------------------------- the page
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
    A('<p class="faint">Data: lab cosmic stand, weekend scan 19-20 Nov 2026, all runs 268-286 except the pulser runs 270 '
      f'and 283. Reference configuration (gain 4.5, trimmed thresholds Km18): runs {sorted(set(REF["runs"]))}, '
      f'{REF["n_scint"]:,} scintillator coincidences; all configurations: {ALL["n_scint"]:,}. Code: '
      '<code>cosmic_vmm_analysis/muon_funnel.py</code> (same cuts as <code>efficiency_tp.py</code>).</p>')

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

    A('<h2>4. Each class, with a real event</h2>')
    A('<p>Each display is the three chambers seen from above; only the read-out area is drawn (the odd shape is the '
      'P2 pad fan, connectors 3-5). Light squares = pads. <b>Coloured</b> = a hit in the muon window (blue IN, orange '
      'MID, green OUT); <b>dark grey</b> = a hit on a dead/masked pad; <b>outlined</b> = a hit outside the muon window. '
      'The dashed circle = the 15 mm search radius around the prediction in the probe (here P2_MID). Hover a pad for '
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
    df = REF['probes']['P2_MID']['d_far']
    dfd = REF['probes']['P2_MID']['d_found']
    edges = list(range(0, 155, 5))
    vals = [df.get(str(e), 0) + dfd.get(str(e), 0) for e in edges[:-1]]
    beyond = sum(v for k, v in df.items() if k != 'inf' and int(k) >= 150)
    A(f'<figure>{svg_hist(edges, vals, "P2_MID: distance of the nearest probe hit to the prediction (mm), REF - log scale", w=900, log=True, mark=15)}'
      f'<figcaption>Muons expected on a live pad whose probe fired in the muon window. Left of the 15 mm line = found '
      f'(or late). Right of it = FAR: {sum(vals[3:6]):,} muons at 15-30 mm (the muon, slightly mis-predicted), '
      f'{sum(vals[6:]) + beyond:,} further out (flat: unrelated hits).</figcaption></figure>')
    A('<h3>5b. Per chip: is the loss everywhere or in a few chips?</h3>')
    A('<p>The chip is the one that reads the predicted pad. If the loss were a few broken chips it would show here; '
      'for P2_MID it does not - every chip loses 16-22 %.</p>')
    for p in PROBES:
        A(f'<figure>{svg_chips(p, REF)}<figcaption>{p} as probe (REF).</figcaption></figure>')
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

    A('<h2>6. So where does the inefficiency arise?</h2>')
    m = exp['P2_MID']
    n3 = m['n3']
    near = sum(vals[3:6])
    A('<ol>'
      f'<li><b>Before the test (not inefficiency).</b> Of {m["n0"]:,} coincidences only {pct(m["n3"], m["n0"])} become a '
      'test of P2_MID. Most are lost because a tag chamber did not fire - mostly geometry: the muon did not cross its '
      'read-out area. Another ~15 % of the tagged muons are dropped for crosstalk pads in a tag. None of these count '
      'as a miss.</li>'
      f'<li><b>In the test, P2_MID: {pct(n3 - m["found"], n3)} of the expected muons are not found:</b><ul>'
      f'<li>{pct(m["probe_silent"], n3)} the chamber records nothing at all - electronics: signal below threshold '
      '(read-out, dead time and sparks are < 0.6 points of it);</li>'
      f'<li>{pct(m["probe_far"] - near, n3)} only an unrelated hit far away - the same electronics loss, hidden;</li>'
      f'<li>{pct(near, n3)} the muon fired 15-30 mm away - a calculation loss (radius);</li>'
      f'<li>{pct(m["probe_near_late"], n3)} late with respect to the tag hit - a calculation loss (clock).</li></ul></li>'
      f'<li><b>Electronics ≈ {pct(m["probe_silent"] + m["probe_far"] - near, n3)}, calculation ≈ '
      f'{pct(near + m["probe_near_late"], n3)}.</b> The silent part is the same 12-13 % "no hit at all" that the '
      'reconstruction-free check (scintillator + IN + OUT: 0.875) gave.</li></ol>')
    A('<div class="callout warn"><b>Corrections to what I said before.</b><ul>'
      '<li>I said a probe window built from the scintillator time would recover 1-4 %. With the random hits removed it '
      f'is {pct(m["probe_near_late"], n3)} (LATE). The 1-4 % came from widening the window on all tags without the '
      'random-hit subtraction.</li>'
      '<li>I proposed ranking clusters by time, paddle geometry and angle. The extra pieces turn out to be crosstalk '
      'pads, so a simpler rule (drop single, low-ADC, in-time pads 1-3 channels from a muon channel) should do; to be '
      'tested.</li></ul></div>')
    A('<h2>7. Questions for you</h2><ol>'
      '<li>Is the neighbour-channel crosstalk known for this set-up (chip, hybrid, cables)? Channel differences of 1-3 '
      'suggest it happens before the channel-to-pad routing.</li>'
      '<li>For the write-up, shall we use this chain as the efficiency definition: scintillator-confirmed muons, two '
      'clean tags, prediction on a live pad, probe hit within 15 mm in the muon window?</li>'
      f'<li>The reference set includes run 284 (same configuration name as 268/269). Is that right, or was something '
      'else changed for 284?</li></ol>')
    return '\n'.join(out)


STYLE = """
  :root { --k-geo:#b9b7b0; --k-geo2:#cfcdc6; --k-geo3:#a3a19a; --k-rej:#6250d6; --k-fid:#e87ba4; --k-fid2:#c25b87;
          --k-sil:#e34948; --k-far:#eda100; --k-late:#2a78d6; --k-ok:#1baf7a;
          --c-in:#2a78d6; --c-mid:#eb6834; --c-out:#1baf7a; --pad:#ecebe6; --pad-line:#d6d4cc; --ink:#1d1f24;
          --ink2:#5f6570; --dead:#55534e; }
  @media (prefers-color-scheme: dark) {
    :root { --k-geo:#6f6e68; --k-geo2:#5b5a55; --k-geo3:#84837c; --k-rej:#9085e9; --k-fid:#d55181; --k-fid2:#e87ba4;
            --k-sil:#e66767; --k-far:#c98500; --k-late:#3987e5; --k-ok:#199e70;
            --c-in:#3987e5; --c-mid:#d95926; --c-out:#199e70; --pad:#2a2d33; --pad-line:#3a3e45; --ink:#e8e8e4;
            --ink2:#a5a9b1; --dead:#8c8a84; }
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
  .markl { stroke: var(--k-sil); stroke-width: 1.5; stroke-dasharray: 4 3; }
  .chl { stroke: var(--ink); stroke-width: 1.6; }
  .chip-pad { fill: var(--c-mid); opacity: .28; }
  .sw { display: inline-block; width: .8em; height: .8em; border-radius: 2px; margin-right: .45em; vertical-align: -.05em; }
  .evfig figcaption { font-size: .82rem; }
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
