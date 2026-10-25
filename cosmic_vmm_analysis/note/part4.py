
# ---------------------------------------------------------------- 25 Oct: are these muons? (section 10c)
STR = json.load(open(LOCAL + 'analysis_muons/straight_track_rate_runs76-108.json'))['runs']
STR_SET = [('safe + 16', [92]), ('safe + 8', [93, 96, 100, 104, 108]), ('safe', list(range(76, 92)) + [97, 101, 105]),
           ('safe &minus; 4', [94, 98, 102, 106]), ('safe &minus; 6', [99, 103, 107]), ('safe &minus; 8', [95])]


def _comb(runs, key, ekey):
    rr = [STR[f'run_{r}'] for r in runs]
    h = sum(x['span_s'] for x in rr) / 3600
    return (sum(x[key] for x in rr) / h, math.sqrt(sum(x[ekey] ** 2 for x in rr)) / h, h,
            sum(x['groups'] for x in rr) / h)


STR_ROWS = [(nm, runs, _comb(runs, 'signal_chip', 'signal_chip_err'), _comb(runs, 'signal', 'signal_err')) for nm, runs in STR_SET]
# chamber-level tag-and-probe (tag_and_probe.py, det2 probe), from the lxplus logs
TP_ROWS = [('safe (16 repeats, runs 76&ndash;91)', '0.52&ndash;0.58'), ('safe + 8 (runs 93, 96)', '0.50, 0.46'),
           ('safe + 16 (run 92)', '0.49'), ('safe &minus; 4 (run 94)', '0.36'), ('safe &minus; 8 (run 95, a quarter of the files)', '0.04')]
# wiring test (orientation_check.py): fraction of single-hit 3-chamber events with det2 within 30 mm of the midpoint
WIRE = [('SPS map, held-out events', 13.6, 'var(--s1)'), ('best other wiring, events it was fitted on', 15.4, 'var(--c2)'),
        ('best other wiring, held-out events', 12.1, 'var(--s2)'), ('shuffled events (chance)', 6.9, 'var(--fg-faint)')]
# fixed noise pattern in the low-threshold "tracks": share carried by the 10 most frequent pad combinations
PATTERN = [('safe &minus; 8 (run 95)', 2443, 0.59), ('safe &minus; 6 (run 99)', 409, 0.44), ('safe &minus; 6 (run 103)', 381, 0.45),
           ('safe &minus; 6 (run 107)', 467, 0.46)]


def cartoon_midpoint():
    W, H = 640, 300
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="The straight-line test: the middle hit must sit half way between the outer hits">']
    o.append('<text class="t" x="12" y="18" font-size="13" font-weight="600">The straight-line test (cartoon)</text>')
    o.append('<text class="m" x="12" y="34" font-size="11">The chambers lie flat, 31 cm apart. A muon is a straight line, so its hit in the middle chamber is half way between the other two,</text>')
    o.append('<text class="m" x="12" y="48" font-size="11">whatever its angle. Noise hits that happen at the same moment land anywhere.</text>')
    for k, (x0, ttl) in enumerate(((40, 'a muon'), (350, 'noise at the same moment'))):
        o.append(f'<text class="t" x="{x0 + 115}" y="76" font-size="12" font-weight="600" text-anchor="middle">{ttl}</text>')
        ys = (100, 175, 250)
        for y, nm in zip(ys, ('det3 (P2_IN)', 'det2 (P2_MID)', 'det1 (P2_OUT)')):
            o.append(f'<rect x="{x0}" y="{y - 5}" width="230" height="10" rx="3" fill="var(--c1)"/>')
            if k == 0:
                o.append(f'<text class="f" x="{x0 + 236}" y="{y + 4}" font-size="9">{nm}</text>')
        if k == 0:
            xa, xb = x0 + 50, x0 + 170
            o.append(f'<line x1="{xa - 18}" y1="84" x2="{xb + 18}" y2="266" stroke="var(--s1)" stroke-width="2.2"/>')
            for x, y in ((xa, 100), ((xa + xb) / 2, 175), (xb, 250)):
                o.append(f'<circle cx="{x:.0f}" cy="{y}" r="6" fill="var(--s1)"/>')
            o.append(f'<text class="m" x="{(xa + xb) / 2 + 12:.0f}" y="168" font-size="10">half way: straight</text>')
        else:
            for x, y in ((x0 + 40, 100), (x0 + 190, 175), (x0 + 120, 250)):
                o.append(f'<circle cx="{x}" cy="{y}" r="6" fill="var(--s2)"/>')
            o.append(f'<line x1="{x0 + 40}" y1="100" x2="{x0 + 120}" y2="250" stroke="var(--s2)" stroke-width="1.5" stroke-dasharray="4 4"/>')
            o.append(f'<text class="m" x="{x0 + 70}" y="168" font-size="10" text-anchor="end">expected here</text>')
            o.append(f'<text class="m" x="{x0 + 200}" y="194" font-size="10" text-anchor="middle">but it is here</text>')
    o.append('</svg>')
    return '\n'.join(o)


def fig_wiring():
    W, H = 640, 200
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Fraction of events that are straight lines for different wirings">']
    o.append('<text class="t" x="12" y="18" font-size="13" font-weight="600">Is the channel-to-pad map right? Straight events with each wiring</text>')
    o.append('<text class="m" x="12" y="34" font-size="11">% of events (one hit per chamber) with the middle hit within 30 mm of the midpoint; 40 wirings tried per connector</text>')
    l, r_ = 300, 610
    for i, (nm, v, col) in enumerate(WIRE):
        y = 52 + i * 34
        w = (r_ - l) * v / 18
        o.append(f'<text class="t" x="{l - 8}" y="{y + 15}" font-size="11" text-anchor="end">{nm}</text>'
                 f'<g><title>{nm}: {v} %</title><rect x="{l}" y="{y}" width="{w:.1f}" height="22" rx="4" fill="{col}"/></g>'
                 f'<text class="t" x="{l + w + 6:.1f}" y="{y + 15}" font-size="11">{v} %</text>')
    o.append('</svg>')
    return '\n'.join(o)


def fig_straight_rate():
    W, H = 640, 380
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Straight tracks per hour above chance at each threshold setting">']
    o.append('<text class="t" x="12" y="18" font-size="13" font-weight="600">Straight three-chamber tracks per hour, above chance, at each setting</text>')
    o.append('<text class="m" x="12" y="34" font-size="11">Runs 76&ndash;108 (the long run included). Bars: the strict chance level (each middle hit kept on its chip, channel shuffled);</text>')
    o.append('<text class="m" x="12" y="48" font-size="11">open bars: the loose chance level (middle hits taken from another event). Lines: &plusmn;1 standard error.</text>')
    l, r_, t, bt = 60, 620, 70, 320
    hi = 160

    def Y(v):
        return bt - (bt - t) * max(min(v, hi), -10) / hi
    for a in (0, 40, 80, 120, 160):
        o.append(f'<line x1="{l}" x2="{r_}" y1="{Y(a):.1f}" y2="{Y(a):.1f}" stroke="var(--border)"/>'
                 f'<text class="f" x="{l - 6}" y="{Y(a) + 3:.1f}" font-size="9" text-anchor="end">{a}</text>')
    bw = (r_ - l) / len(STR_ROWS)
    for i, (nm, runs, (s, e, hrs, gph), (s2, e2, _, _)) in enumerate(STR_ROWS):
        x = l + i * bw + 10
        w = (bw - 20) / 2
        clipped = s > hi
        o.append(f'<g><title>{nm.replace("&minus;", "-")}: {s:.0f} +- {e:.0f} per hour (strict), {s2:.0f} +- {e2:.0f} (loose); {hrs:.1f} h, {len(runs)} runs</title>'
                 f'<rect x="{x:.1f}" y="{Y(max(s, 0)):.1f}" width="{w:.1f}" height="{abs(Y(max(s, 0)) - Y(0)):.1f}" rx="3" fill="{"var(--s2)" if clipped or s > 40 else "var(--s1)"}"/>'
                 f'<rect x="{x + w + 3:.1f}" y="{Y(max(s2, 0)):.1f}" width="{w - 3:.1f}" height="{abs(Y(max(s2, 0)) - Y(0)):.1f}" rx="3" fill="none" stroke="var(--fg-faint)" stroke-width="1.5"/></g>')
        cx = x + w / 2
        if not clipped:
            o.append(f'<line x1="{cx:.1f}" x2="{cx:.1f}" y1="{Y(s - e):.1f}" y2="{Y(s + e):.1f}" stroke="var(--fg)" stroke-width="1.3"/>')
        lab = f'{s:.0f} &plusmn; {e:.0f}' if not clipped else f'&#8593; {s:.0f}'
        o.append(f'<text class="t" x="{cx:.1f}" y="{(Y(s + e) if not clipped else t) - 6:.1f}" font-size="10" text-anchor="middle">{lab}</text>')
        o.append(f'<text class="m" x="{x + (bw - 20) / 2:.1f}" y="{bt + 14}" font-size="10" text-anchor="middle">{nm}</text>'
                 f'<text class="f" x="{x + (bw - 20) / 2:.1f}" y="{bt + 27}" font-size="9" text-anchor="middle">{hrs:.1f} h, {len(runs)} run{"s" if len(runs) > 1 else ""}</text>')
    o.append(f'<text class="f" x="{(l + r_) / 2}" y="{bt + 46}" font-size="9" text-anchor="middle">the orange bars are a fixed pattern of hot channels, not muons (next figure)</text>')
    o.append('</svg>')
    return '\n'.join(o)


def fig_pattern():
    W, H = 640, 200
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Share of the low-threshold straight events carried by ten pad combinations">']
    o.append('<text class="t" x="12" y="18" font-size="13" font-weight="600">The low-threshold &#8220;tracks&#8221; repeat the same few pads</text>')
    o.append('<text class="m" x="12" y="34" font-size="11">Share of the straight events carried by the 10 most frequent (det3 pad, det1 pad) combinations. Real muons would spread</text>')
    o.append('<text class="m" x="12" y="48" font-size="11">over hundreds of combinations; one combination alone (VMM 5 ch 46, VMM 11 ch 10, VMM 19 ch 20) is 22 % of run 95.</text>')
    l, r_ = 220, 600
    for i, (nm, n, f) in enumerate(PATTERN):
        y = 64 + i * 32
        o.append(f'<text class="t" x="{l - 8}" y="{y + 15}" font-size="11" text-anchor="end">{nm}</text>'
                 f'<rect x="{l}" y="{y}" width="{r_ - l}" height="22" rx="4" fill="var(--c0)"/>'
                 f'<g><title>{nm.replace("&minus;", "-")}: {n} straight events, {100 * f:.0f} % in 10 combinations</title>'
                 f'<rect x="{l}" y="{y}" width="{(r_ - l) * f:.1f}" height="22" rx="4" fill="var(--s2)"/></g>'
                 f'<text class="x3" x="{l + 8}" y="{y + 15}" font-size="11">{100 * f:.0f} % of {n}</text>')
    o.append('</svg>')
    return '\n'.join(o)


def tp_table():
    rows = ''.join(f'<tr><th>{a}</th><td>{b}</td></tr>' for a, b in TP_ROWS)
    return f'<table class="bltab"><tr><th>setting</th><th>det2 fired when det3 and det1 fired</th></tr>{rows}</table>'


def str_table():
    rows = ''.join(f'<tr><th>{nm}</th><td>{len(runs)}</td><td>{h:.1f}</td><td>{g:,.0f}</td><td>{s:.0f} &plusmn; {e:.0f}</td><td>{s2:.0f} &plusmn; {e2:.0f}</td></tr>'
                   for nm, runs, (s, e, h, g), (s2, e2, _, _) in STR_ROWS)
    return ('<table class="bltab"><tr><th>setting</th><th>runs</th><th>hours</th><th>3-chamber groups per hour</th>'
            f'<th>straight above chance, strict</th><th>loose</th></tr>{rows}</table>')
