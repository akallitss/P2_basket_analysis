
# ---------------------------------------------------------------- 26-30 Oct: the muons were there (section 15)
S1 = json.load(open(LOCAL + 'analysis_runs76_111/s1_summary_per_run.json'))
SC130 = json.load(open(LOCAL + 'analysis_run130/scint130.json'))
E124 = json.load(open(LOCAL + 'analysis_run124/eff124_dt.json'))
DETS = ((3, 'det3 P2_IN', 'var(--s1)'), (2, 'det2 P2_MID', 'var(--s2)'), (1, 'det1 P2_OUT', 'var(--s3)'))
# straight tracks per hour above the strict chance level, runs 76-91 pooled (8.9 h), against the connector shift
SHIFT = [(-3, 11.1, 5.7), (-2, 6.8, 5.3), (-1, 13.6, 5.1), (0, 10.8, 6.5), (1, 510.7, 9.7), (2, 8.1, 5.1), (3, 7.9, 5.2), (4, 10.8, 5.4)]
# straight tracks per hour (strict chance) in the new set-up, final map
NEW_STR = [('run_121', 'VMM 4, 7, 9 off, discharge vetoed', 604, 43), ('run_124', 'all chips', 813, 49)]
# efficiency of the normal and the high-baseline channels of the burst chips, runs 76-91 summed (chamber-only tags)
BLOCK_EFF = [(4, 0.31, 34, 0.20, 19), (7, 0.18, 22, 0.10, 29), (9, 0.37, 17, 0.12, 31)]
# run_130 trigger board (hybrid H0, VMM 0)
TRIG = [('ch 40', '0.28 Hz', 'coincidence of the three paddles', '99.8&nbsp;% of its pulses come with all three paddle signals, 15&ndash;30&nbsp;ns after them'),
        ('ch 41', '13 Hz', 'top paddle', '74&nbsp;% of its pulses with ch 45'),
        ('ch 45', '23 Hz', 'top paddle', '42&nbsp;% with ch 41'),
        ('ch 44', '15 Hz', 'bottom paddle', '2&ndash;4&nbsp;% with the top paddles; det1, the bottom chamber, answers it most')]
DREAM = {1: (0.81, 'DREAM, cosmic bench, 19 Jul, 430 V'), 2: (0.66, 'DREAM, cosmic bench, 9 Jul, 430 V (65&ndash;67&nbsp;%)')}


# weekend scan adwk1 (31 Oct - 1 Nov): 30-minute slices of the chamber response to the scintillator coincidence
WK_TL = json.load(open(LOCAL + 'analysis_wk1/timeline.json'))['slices'] + \
    [s for s in json.load(open(LOCAL + 'analysis_wk1/timeline_131_134.json'))['slices'] if s['run'] != 'run_130end']
WK_K = {130: 0, 131: 0, 132: 0, 133: 64, 134: 64, 135: 64, 140: 24, 143: 0, 146: -8, 149: -12, 154: 0, 157: 40, 160: 12,
        163: -4, 166: -10}
# (top OR) AND bottom coincidence against the hardware triple, run_130
TRIG_OR = [('hardware triple (ch 40)', 1002, 237), ('top 41 AND bottom 44', 1170, 263), ('top 45 AND bottom 44', 2107, 394),
           ('(top 41 OR top 45) AND bottom 44', 3277, 659)]


def fig_weekend_timeline():
    from datetime import datetime as _dt
    W, H = 640, 420
    t0 = _dt(2026, 10, 30, 18, 0)
    span = 40.0

    def hrs(s):
        return (_dt.fromisoformat(s['start']) - t0).total_seconds() / 3600 + s['files'] / 120
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Chamber response and muon pulse height against time during the weekend scan">']
    o.append('<text class="t" x="12" y="18" font-size="13" font-weight="600">The weekend scan: the chambers faded, whatever the threshold</text>')
    o.append('<text class="m" x="12" y="34" font-size="11">Each dot: 30 minutes (2 to 4 minutes for the short runs 131&ndash;134). Shaded bands: the 3-hour points, labelled with</text>')
    o.append('<text class="m" x="12" y="48" font-size="11">the threshold offset K. The scintillators triggered at the same rate all along (about 1000 per hour).</text>')
    l, r_ = 60, 620
    panels = ((70, 220, 'response to the scintillators', 'resp', 0.6), (250, 380, 'median ADC of the muon hits', 'adc_median', 360))

    def X(h):
        return l + (r_ - l) * h / span
    runs = {}
    for s in WK_TL:
        n = int(s['run'].split('_')[1])
        h = hrs(s)
        runs.setdefault(n, [h, h])
        runs[n][0] = min(runs[n][0], h - s['files'] / 120)
        runs[n][1] = max(runs[n][1], h + s['files'] / 120)
    for t, bt, lab, key, hi in panels:
        for k, (a, b) in sorted(runs.items()):
            if b - a > 1:
                o.append(f'<rect x="{X(a):.1f}" y="{t}" width="{X(b) - X(a):.1f}" height="{bt - t}" fill="var(--c0)" opacity="0.7"/>')
                if key == 'resp':
                    o.append(f'<text class="f" x="{(X(a) + X(b)) / 2:.1f}" y="{t + 10}" font-size="8.5" text-anchor="middle">K {WK_K[k]:+d}</text>'
                             if WK_K[k] else f'<text class="f" x="{(X(a) + X(b)) / 2:.1f}" y="{t + 10}" font-size="8.5" text-anchor="middle">K 0</text>')

        def Y(v):
            return bt - (bt - t) * max(min(v, hi), 0) / hi
        for v in ((0, 0.2, 0.4, 0.6) if key == 'resp' else (0, 100, 200, 300)):
            o.append(f'<line x1="{l}" x2="{r_}" y1="{Y(v):.1f}" y2="{Y(v):.1f}" stroke="var(--border)"/>'
                     f'<text class="f" x="{l - 6}" y="{Y(v) + 3:.1f}" font-size="9" text-anchor="end">{v:g}</text>')
        o.append(f'<text class="m" x="{l - 44}" y="{(t + bt) / 2:.1f}" font-size="10" transform="rotate(-90 {l - 44} {(t + bt) / 2:.1f})" text-anchor="middle">{lab}</text>')
        for d, nm, col in DETS:
            for s in WK_TL:
                v = s[nm][key]
                if v is None or (key == 'adc_median' and (s[nm]['n_hit'] < 15)):
                    continue
                o.append(f'<g><title>{s["run"]} {s["start"][5:16]}: {nm} {v:.3g}</title>'
                         f'<circle cx="{X(hrs(s)):.1f}" cy="{Y(v):.1f}" r="2.8" fill="{col}"/></g>')
    for h, lab in ((0, 'Fri 18:00'), (6, '24:00'), (12, 'Sat 06:00'), (18, '12:00'), (24, '18:00'), (30, '24:00'), (36, 'Sun 06:00')):
        o.append(f'<text class="f" x="{X(h):.1f}" y="{panels[1][1] + 14}" font-size="9" text-anchor="middle">{lab}</text>')
    for i, (d, nm, col) in enumerate(DETS):
        o.append(f'<circle cx="{l + 10 + i * 130}" cy="{H - 10}" r="4" fill="{col}"/><text class="m" x="{l + 18 + i * 130}" y="{H - 6}" font-size="10">{nm}</text>')
    o.append(f'<text class="f" x="{r_}" y="{H - 6}" font-size="9" text-anchor="end">ADC shown where a slice has at least 15 muon hits</text>')
    o.append('</svg>')
    return '\n'.join(o)


def trig_or_table():
    rows = ''.join(f'<tr><th>{a}</th><td>{b}</td><td>{c}</td><td>{c / TRIG_OR[0][2]:.1f}&times;</td></tr>' for a, b, c in TRIG_OR)
    return ('<table class="bltab"><tr><th>trigger (run_130, built offline from the recorded paddles)</th><th>per hour</th>'
            f'<th>with a hit in all three chambers, per hour</th><th>gain</th></tr>{rows}</table>')


def _sc(d, key):
    nm = {3: 'det3 P2_IN', 2: 'det2 P2_MID', 1: 'det1 P2_OUT'}[d]
    r = SC130[f'probe_{nm}'][key]
    return r['eff'], r['err']


def _tp124(d, key):
    P = E124[f'probe_det{d}']
    r = P['any_position'] if key == 'any' else P['vs_probe_r'][key]
    return r['eff'], r.get('err') or 0.0


def fig_shift_scan():
    W, H = 640, 250
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Straight tracks per hour against the connector assumed for the hybrids">']
    o.append('<text class="t" x="12" y="18" font-size="13" font-weight="600">Which connectors were the old runs really read on?</text>')
    o.append('<text class="m" x="12" y="34" font-size="11">Runs 76&ndash;91 (8.9 hours): straight tracks per hour above the strict chance level, if the three hybrids of each chamber</text>')
    o.append('<text class="m" x="12" y="48" font-size="11">sit on the connectors below. The run files said c4&ndash;c6. Lines: &plusmn;1 standard error.</text>')
    l, r_, t, bt = 60, 620, 70, 200
    hi = 550

    def Y(v):
        return bt - (bt - t) * max(min(v, hi), 0) / hi
    for a in (0, 100, 200, 300, 400, 500):
        o.append(f'<line x1="{l}" x2="{r_}" y1="{Y(a):.1f}" y2="{Y(a):.1f}" stroke="var(--border)"/>'
                 f'<text class="f" x="{l - 6}" y="{Y(a) + 3:.1f}" font-size="9" text-anchor="end">{a}</text>')
    bw = (r_ - l) / len(SHIFT)
    for i, (s, v, e) in enumerate(SHIFT):
        x = l + i * bw + 12
        w = bw - 24
        col = 'var(--s3)' if s == 1 else ('var(--s2)' if s == 0 else 'var(--c2)')
        lab = f'c{4 + s}&ndash;c{6 + s}'
        o.append(f'<g><title>{lab.replace("&ndash;", "-")}: {v:.0f} +- {e:.0f} per hour</title>'
                 f'<rect x="{x:.1f}" y="{Y(v):.1f}" width="{w:.1f}" height="{max(Y(0) - Y(v), 1):.1f}" rx="3" fill="{col}"/></g>'
                 f'<text class="t" x="{x + w / 2:.1f}" y="{Y(v) - 5:.1f}" font-size="10" text-anchor="middle">{v:.0f}</text>'
                 f'<text class="m" x="{x + w / 2:.1f}" y="{bt + 14}" font-size="10" text-anchor="middle">{lab}</text>')
        if s in (0, 1):
            o.append(f'<text class="f" x="{x + w / 2:.1f}" y="{bt + 27}" font-size="9" text-anchor="middle">{"run files" if s == 0 else "the data"}</text>')
    o.append(f'<text class="f" x="{(l + r_) / 2}" y="{bt + 44}" font-size="9" text-anchor="middle">connectors of the three hybrids (same shift on all three chambers)</text>')
    o.append('</svg>')
    return '\n'.join(o)


def fig_eff_per_run():
    W, H = 640, 300
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Efficiency per run against the threshold setting, three chambers">']
    o.append('<text class="t" x="12" y="18" font-size="13" font-weight="600">The threshold scan of 24&ndash;25 October, with the right map: one dot per run</text>')
    o.append('<text class="m" x="12" y="34" font-size="11">Efficiency from the chambers alone (tag-and-probe, probe hit within 15 mm), runs 76&ndash;111. These numbers are lower bounds</text>')
    o.append('<text class="m" x="12" y="48" font-size="11">(section 15c); compare the settings with each other, not with DREAM. Left of the dashed line the noise floods the tags.</text>')
    ks = [16, 8, 0, -4, -6, -8]
    pw = 190
    for j, (d, nm, col) in enumerate(DETS):
        l = 40 + j * (pw + 15)
        r_ = l + pw - 10
        t, bt = 70, 240

        def X(k):
            return r_ - (r_ - l) * (k + 9) / 26

        def Y(v):
            return bt - (bt - t) * v / 0.6
        o.append(f'<text class="t" x="{(l + r_) / 2:.1f}" y="66" font-size="11" text-anchor="middle">{nm}</text>')
        for a in (0, 0.2, 0.4, 0.6):
            o.append(f'<line x1="{l}" x2="{r_}" y1="{Y(a):.1f}" y2="{Y(a):.1f}" stroke="var(--border)"/>'
                     + (f'<text class="f" x="{l - 4}" y="{Y(a) + 3:.1f}" font-size="9" text-anchor="end">{a:.1f}</text>' if j == 0 else ''))
        o.append(f'<line x1="{X(-2):.1f}" x2="{X(-2):.1f}" y1="{t}" y2="{bt}" stroke="var(--fg-faint)" stroke-dasharray="3 3"/>')
        for k in ks:
            rr = [r for r in S1 if r['K'] == k and r.get(f'eff{d}') is not None]
            for i, r in enumerate(rr):
                x = X(k) + (i - (len(rr) - 1) / 2) * min(1.6, 14 / max(len(rr), 1))
                e = r[f'eff{d}']
                o.append(f'<g><title>run {r["run"]}, {r["setting"]}: {e:.3f} +- {r[f"err{d}"]:.3f}</title>'
                         f'<circle cx="{x:.1f}" cy="{Y(e):.1f}" r="2.6" fill="{col}" opacity="{0.9 if k >= 0 else 0.45}"/></g>')
            o.append(f'<text class="m" x="{X(k):.1f}" y="{bt + 13}" font-size="9" text-anchor="middle">{k:+d}' + '</text>' if k else
                     f'<text class="m" x="{X(k):.1f}" y="{bt + 13}" font-size="9" text-anchor="middle">safe</text>')
        o.append(f'<text class="f" x="{(l + r_) / 2:.1f}" y="{bt + 27}" font-size="9" text-anchor="middle">threshold: safe + K DAC</text>')
    o.append('</svg>')
    return '\n'.join(o)


def fig_scint():
    W, H = 640, 300
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Efficiency from the chambers alone and with the scintillator reference">']
    o.append('<text class="t" x="12" y="18" font-size="13" font-weight="600">The scintillator reference: how efficient are the chambers really?</text>')
    o.append('<text class="m" x="12" y="34" font-size="11">Same thresholds and voltages. Grey: chambers alone (run_124). Colour: the scintillator coincidence fired and the other two</text>')
    o.append('<text class="m" x="12" y="48" font-size="11">chambers have one clean cluster each (run_130, about 300 events per chamber). Solid: probe hit within 50 mm; light: anywhere.</text>')
    l, r_, t, bt = 50, 620, 66, 250

    def Y(v):
        return bt - (bt - t) * v
    for a in (0, 0.2, 0.4, 0.6, 0.8, 1.0):
        o.append(f'<line x1="{l}" x2="{r_}" y1="{Y(a):.1f}" y2="{Y(a):.1f}" stroke="var(--border)"/>'
                 f'<text class="f" x="{l - 6}" y="{Y(a) + 3:.1f}" font-size="9" text-anchor="end">{a:.1f}</text>')
    gw = (r_ - l) / 3
    for j, (d, nm, col) in enumerate(DETS):
        x0 = l + j * gw + 18
        bars = [('chambers alone', _tp124(d, '50'), _tp124(d, 'any'), 'var(--fg-faint)'),
                ('scintillator', _sc(d, 'within_50mm'), _sc(d, 'any_position'), col)]
        for i, (lab, (e50, s50), (ea, sa), c) in enumerate(bars):
            x = x0 + i * 72
            o.append(f'<g><title>{nm}, {lab}: within 50 mm {e50:.2f} +- {s50:.2f}, anywhere {ea:.2f} +- {sa:.2f}</title>'
                     f'<rect x="{x:.1f}" y="{Y(ea):.1f}" width="60" height="{Y(0) - Y(ea):.1f}" rx="3" fill="{c}" opacity="0.35"/>'
                     f'<rect x="{x:.1f}" y="{Y(e50):.1f}" width="60" height="{Y(0) - Y(e50):.1f}" rx="3" fill="{c}"/></g>'
                     f'<text class="t" x="{x + 30:.1f}" y="{Y(ea) - 5:.1f}" font-size="10" text-anchor="middle">{e50:.2f} / {ea:.2f}</text>'
                     f'<text class="m" x="{x + 30:.1f}" y="{bt + 13}" font-size="9" text-anchor="middle">{lab}</text>')
        if d in DREAM:
            v, txt = DREAM[d]
            o.append(f'<g><title>{txt}</title><line x1="{x0 - 6:.1f}" x2="{x0 + 138:.1f}" y1="{Y(v):.1f}" y2="{Y(v):.1f}" stroke="var(--fg)" stroke-width="1.6" stroke-dasharray="5 3"/></g>'
                     f'<text class="m" x="{x0 + 140:.1f}" y="{Y(v) + 3:.1f}" font-size="9">DREAM</text>')
        o.append(f'<text class="t" x="{x0 + 66:.1f}" y="{bt + 30}" font-size="11" text-anchor="middle">{nm}</text>')
    o.append('</svg>')
    return '\n'.join(o)


def fig_block_eff():
    W, H = 640, 190
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Efficiency of the normal and the high-baseline channels of VMM 4, 7 and 9">']
    o.append('<text class="t" x="12" y="18" font-size="13" font-weight="600">The high-baseline channels of VMM 4, 7 and 9 lose muons</text>')
    o.append('<text class="m" x="12" y="34" font-size="11">Runs 76&ndash;91 summed (chambers alone, so lower bounds; compare the two bars of a chip with each other).</text>')
    l, r_ = 150, 600
    for i, (v, en, nn, eb, nb) in enumerate(BLOCK_EFF):
        y = 50 + i * 44
        o.append(f'<text class="t" x="{l - 10}" y="{y + 22}" font-size="12" text-anchor="end">VMM {v}</text>')
        for k, (e, n, lab, col) in enumerate(((en, nn, 'normal channels', 'var(--s1)'), (eb, nb, 'high-baseline block', 'var(--s2)'))):
            yy = y + k * 18
            o.append(f'<g><title>VMM {v} {lab}: {e:.2f} ({n} channels)</title><rect x="{l}" y="{yy}" width="{(r_ - l) * e / 0.7:.1f}" height="15" rx="3" fill="{col}"/></g>'
                     f'<text class="m" x="{l + (r_ - l) * e / 0.7 + 6:.1f}" y="{yy + 12}" font-size="10">{e:.2f} &middot; {lab} ({n} ch)</text>')
    o.append('</svg>')
    return '\n'.join(o)


def trig_table():
    rows = ''.join(f'<tr><th>VMM 0 {a}</th><td style="text-align:left">{b}</td><td style="text-align:left">{c}</td><td style="text-align:left">{d}</td></tr>'
                   for a, b, c, d in TRIG)
    return ('<table class="bltab"><tr><th>channel</th><th>rate</th><th>what it is</th><th>evidence</th></tr>' + rows + '</table>')


def scint_table():
    rows = ''
    for d, nm, _ in DETS:
        sc = SC130['chambers'][nm]
        rows += (f'<tr><th>{nm}</th><td>{sc["any_hit"]:.2f}</td><td>{_sc(d, "any_position")[0]:.2f} &plusmn; {_sc(d, "any_position")[1]:.2f}</td>'
                 f'<td>{_sc(d, "within_50mm")[0]:.2f}</td><td>{_sc(d, "within_15mm")[0]:.2f}</td>'
                 f'<td>{_tp124(d, "50")[0]:.2f} / {_tp124(d, "any")[0]:.2f}</td></tr>')
    return ('<div style="overflow-x:auto"><table class="bltab"><tr><th></th><th>coincidence<br>only</th><th>+ clean tags:<br>anywhere</th><th>within<br>50 mm</th>'
            f'<th>within<br>15 mm</th><th>chambers alone (run_124):<br>50 mm / anywhere</th></tr>{rows}</table></div>')
