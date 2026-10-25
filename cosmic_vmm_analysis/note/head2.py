#!/usr/bin/env python3
"""Builds the unlisted note on the VMM cosmic threshold scan: real numbers in the charts,
labelled cartoons for the ideas. Output: the note file given as argv[1]."""
import math
import sys

NOTE = sys.argv[1]

# ----------------------------------------------------------------- real data
D = [50, 40, 34, 30, 27, 24, 22, 20, 18, 16]            # DAC above the assumed baseline
CHIPS = [2, 3, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
ROWS = {  # kHz per chip, first 40 MB of each point; runs 4..13 = d 50..16
    50: [.04, .04, .04, .04, .04, .01, .03, .01, .02, .02, .06, .13, .06, .09, .08, .10],
    40: [.04, .05, .05, .04, 10.92, .02, .03, .02, .02, .02, .06, .15, .08, .09, .09, .11],
    34: [.05, .06, .09, .04, 489.17, .02, .03, .02, .06, .03, .23, 3.94, .08, .11, .08, .21],
    30: [.03, .26, .61, .03, 332.36, .02, .02, .01, 3.44, .02, 4.65, 110.67, .11, .11, .10, 3.27],
    27: [.06, 3.89, 7.28, .05, 60.40, .19, .03, .03, 50.88, .12, 31.78, 532.86, .14, .46, .15, 25.75],
    24: [.03, 37.57, 62.16, .02, 4.22, 2.21, .03, .25, 395.25, 1.52, 154.42, 1284.53, .32, 4.22, .90, 132.51],
    22: [.03, 129.80, 181.12, .02, .51, 10.06, .13, 1.80, 1197.76, 7.85, 329.49, 2120.12, 1.82, 18.97, 4.49, 323.09],
    20: [.02, 357.91, 416.00, .03, .04, 34.59, .86, 10.48, 2960.24, 35.05, 696.14, 2777.16, 9.40, 67.43, 20.96, 692.32],
    18: [.08, 814.98, .33, .02, .05, 101.88, 4.85, 41.63, 1768.60, 131.30, 441.60, 3651.86, 43.98, 189.75, 79.88, 1248.86],
    16: [.33, 1536.91, 2.89, .02, .37, 246.99, 24.40, 148.50, 3142.34, 352.57, 1091.87, 5153.04, 160.98, 450.85, 270.51, 2196.68],
}
MB_MIN = {50: 5, 40: 9, 34: 179, 30: 185, 27: 267, 24: 752, 22: 1530, 20: 2867, 18: 2560, 16: 4300}
RUN = {50: 4, 40: 5, 34: 6, 30: 7, 27: 8, 24: 9, 22: 10, 20: 11, 18: 12, 16: 13}
BASE = {2: 208, 3: 182, 4: 199, 5: 212, 7: 200, 8: 197, 9: 191, 10: 210, 11: 204, 12: 210,
        13: 188, 14: 201, 15: 180, 16: 199, 17: 187, 18: 182, 19: 204}
CHAMBER = lambda v: 'det3 · IN' if v <= 7 else ('det2 · MID' if v <= 13 else 'det1 · OUT')


def esc(s):
    return s.replace('&', '&amp;').replace('<', '&lt;')


def fmt(v):
    return f'{v:.2f}' if v < 1 else (f'{v:.1f}' if v < 100 else f'{v:.0f}')


def binof(v):
    return 0 if v < .1 else 1 if v < 1 else 2 if v < 10 else 3 if v < 100 else 4 if v < 1000 else 5


# ----------------------------------------------------------------- figures
def fig_heatmap(rows=None, title='Hit rate per chip (kHz) at each threshold step', runs=None, foot6='* VMM 6: sdt fixed at 230 at every step (no baseline estimate), so d does not apply to it.'):
    rows = rows or ROWS
    runs = runs or RUN
    cw, ch, left, top = 51, 25, 126, 94
    chips = sorted(CHIPS + [4, 9])
    W = left + cw * len(D) + 14
    H = top + ch * len(chips) + 84
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Heat map of hit rate per chip against threshold">']
    o.append(f'<text class="t" x="{left - 100}" y="18" font-size="13" font-weight="600">{title}</text>')
    o.append(f'<text class="m" x="{left - 100}" y="36" font-size="11">columns: threshold step d above the assumed baseline, quiet → noisy</text>')
    o.append(f'<text class="m" x="{left - 100}" y="50" font-size="11">rate from the first 40 MB of each 2-minute point</text>')
    for j, d in enumerate(D):
        x = left + j * cw + cw / 2
        o.append(f'<text class="m" x="{x}" y="{top - 22}" font-size="11" text-anchor="middle">d={d}</text>')
        o.append(f'<text class="f" x="{x}" y="{top - 9}" font-size="10" text-anchor="middle">run_{runs[d]}</text>')
    for i, v in enumerate(chips):
        y = top + i * ch
        lab = f'VMM {v}'
        o.append(f'<text class="t" x="{left - 100}" y="{y + 17}" font-size="12">{lab}{"*" if v == 6 else ""}</text>')
        o.append(f'<text class="f" x="{left - 6}" y="{y + 17}" font-size="9" text-anchor="end">{CHAMBER(v)}</text>')
        if v in (4, 9):
            o.append(f'<rect x="{left}" y="{y + 1}" width="{cw * len(D) - 2}" height="{ch - 2}" rx="3" fill="url(#hatch)" stroke="var(--border)"/>')
            o.append(f'<text class="m" x="{left + cw * len(D) / 2}" y="{y + 17}" font-size="11" text-anchor="middle">fully masked: hopping burst (section 4)</text>')
            continue
        k = CHIPS.index(v)
        first = None
        for j, d in enumerate(D):
            val = rows[d][k]
            b = binof(val)
            x = left + j * cw
            if first is None and val >= 1:
                first = j
            stroke = ' stroke="var(--alert)" stroke-width="2.5"' if first == j and val >= 1 else ''
            o.append(f'<g><title>VMM {v}, d={d} (run_{runs[d]}): {val:g} kHz</title>'
                     f'<rect x="{x + 1}" y="{y + 1}" width="{cw - 2}" height="{ch - 2}" rx="3" class="b{b}"{stroke}/>'
                     f'<text class="x{b}" x="{x + cw / 2}" y="{y + 17}" font-size="11" text-anchor="middle">{fmt(val)}</text></g>')
    ly = top + ch * len(chips) + 20
    labels = ['&lt; 0.1', '0.1 – 1', '1 – 10', '10 – 100', '100 – 1000', '&gt; 1000']
    o.append(f'<text class="m" x="{left - 100}" y="{ly + 12}" font-size="11">kHz</text>')
    for b, lab in enumerate(labels):
        x = left + b * 90
        o.append(f'<rect x="{x}" y="{ly}" width="20" height="14" rx="3" class="b{b}"/>'
                 f'<text class="m" x="{x + 26}" y="{ly + 11}" font-size="11">{lab}</text>')
    o.append(f'<rect x="{left}" y="{ly + 24}" width="20" height="14" rx="3" fill="none" stroke="var(--alert)" stroke-width="2.5"/>'
             f'<text class="m" x="{left + 26}" y="{ly + 35}" font-size="11">first step above 1 kHz</text>')
    o.append(f'<text class="f" x="{left - 100}" y="{ly + 58}" font-size="10">{foot6}</text>')
    o.append('<defs><pattern id="hatch" width="7" height="7" patternUnits="userSpaceOnUse" patternTransform="rotate(45)">'
             '<line x1="0" y1="0" x2="0" y2="7" stroke="var(--fg-faint)" stroke-width="2"/></pattern></defs>')
    o.append('</svg>')
    return '\n'.join(o)


def fig_volume():
    W, H, l, r, t, b = 620, 300, 62, 20, 36, 52
    pw, ph = W - l - r, H - t - b
    lo, hi = math.log10(1), math.log10(10000)
    def Y(v): return t + ph * (1 - (math.log10(v) - lo) / (hi - lo))
    bw = pw / len(D)
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Data volume per minute against threshold step">']
    o.append(f'<text class="t" x="{l}" y="18" font-size="13" font-weight="600">Data written per minute, all chips together</text>')
    for e in range(0, 5):
        y = Y(10 ** e)
        o.append(f'<line x1="{l}" x2="{W - r}" y1="{y}" y2="{y}" stroke="var(--border)"/>'
                 f'<text class="m" x="{l - 8}" y="{y + 4}" font-size="11" text-anchor="end">{10 ** e:,} MB</text>')
    flood = Y(635)
    o.append(f'<line x1="{l}" x2="{W - r}" y1="{flood}" y2="{flood}" stroke="var(--alert)" stroke-dasharray="5 4"/>'
             f'<text class="m" x="{l + 6}" y="{flood - 5}" font-size="11">first test: hot channels on, 635 MB/min</text>')
    for j, d in enumerate(D):
        x = l + j * bw + 6
        v = MB_MIN[d]
        o.append(f'<g><title>d={d} (run_{RUN[d]}): {v} MB/min</title>'
                 f'<rect x="{x}" y="{Y(v)}" width="{bw - 12}" height="{t + ph - Y(v)}" rx="3" fill="var(--s1)"/></g>'
                 f'<text class="t" x="{x + (bw - 12) / 2}" y="{Y(v) - 5}" font-size="11" text-anchor="middle">{v:,}</text>'
                 f'<text class="m" x="{x + (bw - 12) / 2}" y="{t + ph + 16}" font-size="11" text-anchor="middle">{d}</text>')
    o.append(f'<text class="m" x="{l + pw / 2}" y="{H - 10}" font-size="12" text-anchor="middle">threshold step d (DAC above the assumed baseline): quiet → noisy</text>')
    o.append('</svg>')
    return '\n'.join(o)


def fig_onset():
    onset = {}
    for k, v in enumerate(CHIPS):
        if v == 6:
            continue
        for d in D:
            if ROWS[d][k] >= 1:
                onset[v] = d
                break
        else:
            onset[v] = None
    order = sorted(onset, key=lambda v: (-(onset[v] or 0), v))
    W, l, r, t = 620, 100, 100, 50
    rh = 22
    H = t + rh * len(order) + 52
    pw = W - l - r
    def X(d): return l + pw * (50 - d) / (50 - 14)
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="First threshold step where each chip exceeds 1 kHz">']
    o.append(f'<text class="t" x="{l - 80}" y="18" font-size="13" font-weight="600">Where does each chip start to fire (&gt; 1 kHz)?</text>')
    o.append(f'<text class="m" x="{l - 80}" y="34" font-size="11">noisier chips are further left: they start firing while the threshold is still high</text>')
    for d in (50, 40, 30, 20):
        o.append(f'<line x1="{X(d)}" x2="{X(d)}" y1="{t - 6}" y2="{t + rh * len(order)}" stroke="var(--border)"/>'
                 f'<text class="m" x="{X(d)}" y="{t + rh * len(order) + 16}" font-size="11" text-anchor="middle">d={d}</text>')
    for i, v in enumerate(order):
        y = t + i * rh + 11
        o.append(f'<text class="t" x="{l - 12}" y="{y + 4}" font-size="12" text-anchor="end">VMM {v}</text>')
        d = onset[v]
        if d is None:
            o.append(f'<g><title>VMM {v}: never above 1 kHz down to d=16</title><line x1="{X(50)}" x2="{X(16)}" y1="{y}" y2="{y}" stroke="var(--fg-faint)" stroke-dasharray="3 4"/>'
                     f'<text class="m" x="{X(16) + 8}" y="{y + 4}" font-size="11">quiet to d=16</text></g>')
        else:
            hollow = v in (7,)
            fill = 'var(--surface)' if hollow else 'var(--s1)'
            o.append(f'<g><title>VMM {v}: first above 1 kHz at d={d}{" (burst, not smooth noise)" if hollow else ""}</title>'
                     f'<line x1="{X(50)}" x2="{X(d)}" y1="{y}" y2="{y}" stroke="var(--fg-faint)" stroke-width="2"/>'
                     f'<circle cx="{X(d)}" cy="{y}" r="6" fill="{fill}" stroke="var(--s1)" stroke-width="2.5"/></g>')
    ly = t + rh * len(order) + 36
    o.append(f'<circle cx="{l}" cy="{ly}" r="6" fill="var(--s1)"/><text class="m" x="{l + 12}" y="{ly + 4}" font-size="11">first step above 1 kHz</text>'
             f'<circle cx="{l + 170}" cy="{ly}" r="6" fill="var(--surface)" stroke="var(--s1)" stroke-width="2.5"/><text class="m" x="{l + 182}" y="{ly + 4}" font-size="11">VMM 7: a burst, not smooth noise</text>')
    o.append('</svg>')
    return '\n'.join(o)


def fig_baselines():
    chips = sorted(BASE)
    W, H, l, r, t, b = 620, 280, 62, 120, 44, 46
    pw, ph = W - l - r, H - t - b
    lo, hi = 170, 220
    def Y(v): return t + ph * (1 - (v - lo) / (hi - lo))
    bw = pw / len(chips)
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Estimated baseline per chip with a common threshold of 200">']
    o.append(f'<text class="t" x="{l}" y="18" font-size="13" font-weight="600">The chips do not sit at the same baseline (July estimates, sdt DAC)</text>')
    o.append(f'<text class="m" x="{l}" y="33" font-size="11">VMM 6 had no estimate; VMM 4, 7, 9 are extrapolated</text>')
    for v in (170, 180, 190, 200, 210, 220):
        o.append(f'<line x1="{l}" x2="{W - r}" y1="{Y(v)}" y2="{Y(v)}" stroke="var(--border)"/>'
                 f'<text class="m" x="{l - 8}" y="{Y(v) + 4}" font-size="11" text-anchor="end">{v}</text>')
    o.append(f'<rect x="{l}" y="{t}" width="{pw}" height="{Y(200) - t}" fill="var(--alert)" opacity="0.10"/>')
    o.append(f'<line x1="{l}" x2="{W - r}" y1="{Y(200)}" y2="{Y(200)}" stroke="var(--alert)" stroke-width="2"/>'
             f'<text class="t" x="{W - r + 8}" y="{Y(200) + 4}" font-size="11" font-weight="600">common sdt = 200</text>')
    for j, v in enumerate(chips):
        x = l + j * bw + bw / 2
        extra = v in (4, 7, 9)
        bad = BASE[v] > 200
        o.append(f'<g><title>VMM {v}: baseline ≈ {BASE[v]}{" (extrapolated)" if extra else ""}</title>'
                 f'<line x1="{x}" x2="{x}" y1="{Y(lo)}" y2="{Y(BASE[v])}" stroke="var(--fg-faint)" stroke-width="2"/>'
                 f'<circle cx="{x}" cy="{Y(BASE[v])}" r="6" fill="{"var(--surface)" if extra else "var(--s1)"}" stroke="var(--s1)" stroke-width="2.5"/></g>'
                 f'<text class="m" x="{x}" y="{H - b + 16}" font-size="10" text-anchor="middle">{v}</text>')
    o.append(f'<text class="m" x="{l + pw / 2}" y="{H - 8}" font-size="12" text-anchor="middle">VMM number</text>')
    o.append('</svg>')
    return '\n'.join(o)


def fig_hop():
    # real: VMM 9 and VMM 4 channel shares, test run_1 (before masking) vs run_3 (after masking)
    W, H = 640, 250
    panels = [
        ('VMM 9, run_1', [(41, 32.2), (63, 19.5), (33, 17.1), (44, 15.1), (42, 13.6)], 30, 'before masking'),
        ('VMM 9, run_3', [(20, 82.3), (21, 17.2)], 30, 'after masking those five'),
        ('VMM 4, run_1', [(49, 99.2)], 30, 'before masking'),
        ('VMM 4, run_3', [(38, 72.6), (36, 27.0)], 30, 'after masking channel 49'),
    ]
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Hot channels before and after masking">']
    o.append(f'<text class="t" x="10" y="18" font-size="13" font-weight="600">Mask the hot channels and the burst moves to others</text>')
    o.append(f'<text class="m" x="10" y="34" font-size="11">bars: share of the chip\'s hits (%); number under a bar = channel; only channels above 10% shown</text>')
    pw = (W - 20) / 4
    for i, (name, bars, _, sub) in enumerate(panels):
        x0 = 10 + i * pw
        o.append(f'<text class="t" x="{x0 + pw / 2}" y="60" font-size="12" text-anchor="middle" font-weight="600">{name}</text>'
                 f'<text class="f" x="{x0 + pw / 2}" y="74" font-size="10" text-anchor="middle">{sub}</text>')
        base, top = 220, 92
        o.append(f'<line x1="{x0 + 8}" x2="{x0 + pw - 8}" y1="{base}" y2="{base}" stroke="var(--border)"/>')
        n = len(bars)
        bw = min(28, (pw - 30) / max(n, 1) - 6)
        total = n * (bw + 6) - 6
        sx = x0 + (pw - total) / 2
        for k, (ch, pct) in enumerate(bars):
            h = (base - top) * pct / 100
            bx = sx + k * (bw + 6)
            o.append(f'<g><title>channel {ch}: {pct}% of this chip\'s hits</title>'
                     f'<rect x="{bx}" y="{base - h}" width="{bw}" height="{h}" rx="3" fill="var(--s1)"/></g>'
                     f'<text class="t" x="{bx + bw / 2}" y="{base - h - 4}" font-size="10" text-anchor="middle">{pct:g}</text>'
                     f'<text class="m" x="{bx + bw / 2}" y="{base + 14}" font-size="10" text-anchor="middle">{ch}</text>')
    o.append('</svg>')
    return '\n'.join(o)


def fig_odd():
    chips = [5, 7]
    W, H, l, r, t, b = 620, 260, 62, 20, 52, 46
    pw, ph = W - l - r, H - t - b
    lo, hi = math.log10(0.01), math.log10(1000)
    def Y(v): return t + ph * (1 - (math.log10(max(v, 0.01)) - lo) / (hi - lo))
    def X(j): return l + pw * j / (len(D) - 1)
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Rate against threshold for VMM 5 and 7">']
    o.append(f'<text class="t" x="{l}" y="18" font-size="13" font-weight="600">Two chips whose rate falls again at low threshold</text>')
    o.append(f'<text class="m" x="{l}" y="33" font-size="11">a lower threshold should never give fewer hits: unexplained</text>')
    for e in (-2, -1, 0, 1, 2, 3):
        y = Y(10 ** e)
        o.append(f'<line x1="{l}" x2="{W - r}" y1="{y}" y2="{y}" stroke="var(--border)"/>'
                 f'<text class="m" x="{l - 8}" y="{y + 4}" font-size="11" text-anchor="end">{10 ** e:g}</text>')
    o.append(f'<text class="m" x="14" y="{t + ph / 2}" font-size="11" transform="rotate(-90 14 {t + ph / 2})" text-anchor="middle">kHz</text>')
    cols = {5: 'var(--s1)', 7: 'var(--s2)'}
    for v in chips:
        k = CHIPS.index(v)
        pts = [(X(j), Y(ROWS[d][k])) for j, d in enumerate(D)]
        o.append(f'<polyline fill="none" stroke="{cols[v]}" stroke-width="2.5" points="' + ' '.join(f'{x:.1f},{y:.1f}' for x, y in pts) + '"/>')
        for j, d in enumerate(D):
            o.append(f'<g><title>VMM {v}, d={d}: {ROWS[d][k]:g} kHz</title><circle cx="{pts[j][0]:.1f}" cy="{pts[j][1]:.1f}" r="4.5" fill="{cols[v]}" stroke="var(--surface)" stroke-width="2"/></g>')
        px, py = pts[4 if v == 7 else 7]
        dy = 18 if v == 7 else -10
        o.append(f'<text class="t" x="{px + 8}" y="{py + dy}" font-size="12" font-weight="600">VMM {v}</text>')
    for j, d in enumerate(D):
        o.append(f'<text class="m" x="{X(j)}" y="{t + ph + 16}" font-size="11" text-anchor="middle">{d}</text>')
    o.append(f'<text class="m" x="{l + pw / 2}" y="{H - 8}" font-size="12" text-anchor="middle">threshold step d (quiet → noisy)</text>')
    o.append('</svg>')
    return '\n'.join(o)


# ------------------------------------------------------------ cartoons
def cartoon_stack():
    W, H = 640, 430
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Schematic of the three-chamber cosmic stand">']
    cham = [(70, 'det3', 'P2_IN', 'H1 H2 H3', 'VMM 2–7'), (170, 'det2', 'P2_MID', 'H4 H5 H6', 'VMM 8–13'), (270, 'det1', 'P2_OUT', 'H7 H8 H9', 'VMM 14–19')]
    for x0, x1 in ((120, 135), (165, 155), (200, 190)):
        o.append(f'<line x1="{x0}" y1="44" x2="{x1}" y2="305" stroke="var(--alert)" stroke-width="2"/>')
    o.append('<text class="m" x="100" y="34" font-size="11">cosmic muons, from the sky</text>')
    for y, det, st, hyb, vm in cham:
        o.append(f'<text class="t" x="6" y="{y + 12}" font-size="12" font-weight="600">{det}</text>'
                 f'<text class="f" x="6" y="{y + 25}" font-size="10">{st}</text>')
        o.append(f'<rect x="60" y="{y}" width="240" height="26" rx="6" fill="var(--s1)" opacity="0.22"/>'
                 f'<rect x="60" y="{y}" width="240" height="26" rx="6" fill="none" stroke="var(--s1)" stroke-width="2"/>'
                 f'<text class="m" x="226" y="{y + 17}" font-size="11">410 / 710 V</text>')
        o.append(f'<rect x="336" y="{y - 2}" width="106" height="30" rx="5" fill="var(--bg-elev)" stroke="var(--border)"/>'
                 f'<text class="t" x="389" y="{y + 12}" font-size="11" text-anchor="middle">{hyb}</text>'
                 f'<text class="f" x="389" y="{y + 24}" font-size="10" text-anchor="middle">{vm}</text>')
        o.append(f'<path d="M300 {y + 13} H336" stroke="var(--fg-faint)" stroke-width="2"/>')
        o.append(f'<path d="M442 {y + 13} H470" stroke="var(--fg-faint)" stroke-width="2"/>')
    o.append('<rect x="470" y="60" width="152" height="250" rx="8" fill="var(--bg-elev)" stroke="var(--border)"/>'
             '<text class="t" x="546" y="90" font-size="13" text-anchor="middle" font-weight="600">ALINX / FEC</text>'
             '<text class="m" x="546" y="108" font-size="11" text-anchor="middle">collects the 9 hybrids</text>'
             '<text class="m" x="546" y="124" font-size="11" text-anchor="middle">sends UDP packets</text>')
    o.append('<path d="M546 310 V350" stroke="var(--fg-faint)" stroke-width="2"/>'
             '<rect x="470" y="350" width="152" height="56" rx="8" fill="var(--bg-elev)" stroke="var(--border)"/>'
             '<text class="t" x="546" y="373" font-size="12" text-anchor="middle" font-weight="600">DAQ PC (vmm_daplxa)</text>'
             '<text class="m" x="546" y="390" font-size="11" text-anchor="middle">dumpcap → pcapng files</text>')
    o.append('<text class="m" x="6" y="338" font-size="11">Samtec 2 m cables on every connector (the beam test used Hitachi)</text>'
             '<text class="m" x="6" y="356" font-size="11">no trigger scintillator: the chips decide for themselves when to record</text>'
             '<text class="m" x="6" y="374" font-size="11">same z spacing as at the SPS beam test: 320 / 630 / 940 mm</text>'
             '<text class="m" x="6" y="392" font-size="11">voltages: mesh 410 V, drift 710 V, written by hand in the logbook</text>')
    o.append('</svg>')
    return '\n'.join(o)

def cartoon_tagprobe():
    W, H = 640, 250
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Tag and probe idea">']
    o.append('<text class="t" x="10" y="18" font-size="13" font-weight="600">Tag and probe: how efficiency is measured with cosmics (next step, not done yet)</text>')
    for k, (y, nm, st) in enumerate(((60, 'IN', 'tag'), (120, 'MID', 'probe'), (180, 'OUT', 'tag'))):
        col = 'var(--s1)' if st == 'tag' else 'var(--s2)'
        o.append(f'<rect x="40" y="{y}" width="260" height="20" rx="5" fill="{col}" opacity="0.25" stroke="{col}" stroke-width="2"/>'
                 f'<text class="t" x="50" y="{y + 14}" font-size="12">{nm}</text><text class="m" x="260" y="{y + 14}" font-size="11">{st}</text>')
    o.append('<line x1="150" y1="40" x2="190" y2="215" stroke="var(--alert)" stroke-width="2"/>')
    o.append('<circle cx="156" cy="70" r="6" fill="var(--s1)"/><circle cx="183" cy="190" r="6" fill="var(--s1)"/>')
    o.append('<circle cx="168" cy="130" r="9" fill="none" stroke="var(--s2)" stroke-width="2.5" stroke-dasharray="3 3"/>')
    o.append('<text class="t" x="320" y="70" font-size="12">1. both outer chambers fire,</text>'
             '<text class="m" x="320" y="86" font-size="11">consistent in time and position:</text>'
             '<text class="m" x="320" y="101" font-size="11">a real muon (the tag)</text>'
             '<text class="t" x="320" y="130" font-size="12">2. extend the line to the middle chamber:</text>'
             '<text class="m" x="320" y="146" font-size="11">did it fire near that spot? (the probe)</text>'
             '<text class="t" x="320" y="175" font-size="12">3. repeat at every threshold step:</text>'
             '<text class="m" x="320" y="191" font-size="11">efficiency = probes that fired / probes</text>'
             '<text class="m" x="320" y="206" font-size="11">swap the probe chamber to cover all three</text>')
    o.append('</svg>')
    return '\n'.join(o)


def cartoon_trim():
    W, H = 640, 250
    import random
    random.seed(4)
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Why a per-channel calibration helps">']
    o.append('<text class="t" x="10" y="18" font-size="13" font-weight="600">Next: calibrate each channel\'s own baseline (cartoon, not measured)</text>')
    for p, (x0, title, spread) in enumerate(((20, 'before: channels differ', 14), (340, 'after: trimmed to a common level', 3))):
        o.append(f'<text class="t" x="{x0}" y="46" font-size="12" font-weight="600">{title}</text>')
        base = 200
        vals = [random.gauss(0, spread) for _ in range(32)]
        for k, v in enumerate(vals):
            x = x0 + 4 + k * 9
            h = 70 + v
            o.append(f'<rect x="{x}" y="{base - h}" width="6" height="{h}" rx="2" fill="var(--s1)"/>')
        ly = base - 70 - max(vals) - 5
        o.append(f'<line x1="{x0}" x2="{x0 + 300}" y1="{ly}" y2="{ly}" stroke="var(--alert)" stroke-width="2"/>')
        o.append(f'<text class="m" x="{x0 + 296}" y="{ly - 5}" font-size="11" text-anchor="end">global threshold</text>')
        o.append(f'<text class="m" x="{x0}" y="{base + 18}" font-size="11">one bar = the baseline of one channel</text>')
    o.append('<text class="m" x="20" y="232" font-size="11">The global threshold has to clear the highest channel. Trimming brings the channels together, so it can sit lower.</text>')
    o.append('</svg>')
    return '\n'.join(o)

def widget_threshold():
    # interactive cartoon: noise + signal, move the threshold
    return '''
<div class="viz-card">
  <svg id="thr-svg" viewBox="0 0 640 260" role="img" aria-label="Interactive threshold cartoon"></svg>
  <label class="slider">Threshold: <input id="thr" type="range" min="2" max="60" value="30"> <span id="thr-v"></span></label>
  <div class="grid-two">
    <div><span class="muted">Muon signals seen (efficiency)</span><div class="big" id="eff"></div></div>
    <div><span class="muted">Noise hits per second, relative</span><div class="big" id="noi"></div></div>
  </div>
  <p class="muted small">Cartoon with invented shapes: the point is the trade-off, not the numbers.</p>
</div>
<script>
(function(){
  var svg=document.getElementById('thr-svg'),sl=document.getElementById('thr');
  var W=640,H=260,l=20,r=20,t=20,b=40,pw=W-l-r,ph=H-t-b,N=80;
  function noise(x){return Math.exp(-0.5*Math.pow((x-8)/4.2,2));}
  function sig(x){var u=(x-30)/9; return x<4?0:0.55*Math.exp(-0.5*(u+Math.exp(-u)));}
  function X(x){return l+pw*x/N;} function Y(y){return t+ph*(1-y/1.05);}
  function path(f,scale){var p='M'+X(0)+','+Y(0);for(var x=0;x<=N;x+=1){p+=' L'+X(x)+','+Y(f(x)*scale);}return p+' L'+X(N)+','+Y(0)+' Z';}
  function tailSum(f,from){var s=0,a=0;for(var x=0;x<=N;x+=0.25){var v=f(x);a+=v;if(x>=from)s+=v;}return s/a;}
  function draw(){
    var th=+sl.value; document.getElementById('thr-v').textContent=th+' (arbitrary units)';
    var eff=tailSum(sig,th), noi=tailSum(noise,th);
    var h='';
    h+='<path d="'+path(noise,1)+'" fill="var(--s3)" opacity="0.35" stroke="var(--s3)" stroke-width="2"/>';
    h+='<path d="'+path(sig,1.5)+'" fill="var(--s1)" opacity="0.30" stroke="var(--s1)" stroke-width="2"/>';
    h+='<rect x="'+X(th)+'" y="'+t+'" width="'+(X(N)-X(th))+'" height="'+ph+'" fill="var(--alert)" opacity="0.08"/>';
    h+='<line x1="'+X(th)+'" x2="'+X(th)+'" y1="'+t+'" y2="'+(t+ph)+'" stroke="var(--alert)" stroke-width="3"/>';
    h+='<line x1="'+l+'" x2="'+(W-r)+'" y1="'+(t+ph)+'" y2="'+(t+ph)+'" stroke="var(--border)"/>';
    h+='<text class="t" x="'+X(8)+'" y="'+(t+ph+16)+'" font-size="12" text-anchor="middle">electronic noise</text>';
    h+='<text class="t" x="'+X(40)+'" y="'+(t+ph+16)+'" font-size="12" text-anchor="middle">muon signal</text>';
    h+='<text class="m" x="'+(W/2)+'" y="'+(H-6)+'" font-size="12" text-anchor="middle">pulse height seen by the chip →   (a channel fires when the pulse is right of the red line)</text>';
    svg.innerHTML=h;
    document.getElementById('eff').textContent=(100*eff).toFixed(0)+' %';
    var n=noi; document.getElementById('noi').textContent=n<0.0005?'≈ 0':(n<0.01?(1000*n).toFixed(1)+' ‰':(100*n).toFixed(1)+' %');
  }
  sl.addEventListener('input',draw); draw();
})();
</script>'''


def cartoon_flow():
    steps = [
        ('Smoke test', 'The July chip configuration at high voltage, 2 minutes: all 18 chips show up, but 1.5 GB of data.'),
        ('Find the flood', 'Almost all hits come from three chips, from single channels: VMM 4, 9 (and later 7).'),
        ('Try masking', 'Switched-off channels go silent, but the burst hops to other channels.'),
        ('Mask whole chips', 'VMM 4 and 9 fully masked: the data rate falls from 635 to 5 MB/min.'),
        ('Scan the threshold', '10 steps, d = 50 → 16, 2 minutes each; 15 chips give a curve (VMM 7 only part way).'),
    ]
    items = ''.join(f'<li><span class="stepn">{i + 1}</span><div><strong>{t}</strong><br><span class="muted">{d}</span></div></li>'
                    for i, (t, d) in enumerate(steps))
    return f'<ol class="steps">{items}</ol>'

