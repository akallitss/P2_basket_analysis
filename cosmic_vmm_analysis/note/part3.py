
# ------------------------------------------------------------------ neighbour read-out (sng 1) data and figures
SNG_DIR = LOCAL + 'analysis_sng1_hvoff/'
BL, BL_FIRST = {}, {}
for _r in csv.DictReader(open(SNG_DIR + 'neighbor_baseline.csv')):
    if _r['adc_nb_median'] != '' and int(_r['n_nb']) >= 50:
        BL[(int(_r['vmm']), int(_r['ch']))] = (int(_r['n_nb']), float(_r['adc_nb_median']), float(_r['adc_nb_std']))
        BL_FIRST[(int(_r['vmm']), int(_r['ch']))] = int(_r['first_file'])
BL_JSON = json.load(open(SNG_DIR + 'neighbor_baseline.json'))
SNG_CHIPS = [v for v in range(2, 20) if v not in (4, 9)]
SNG_STEPS = ['probe', 'o0', '−3', '−6', '−9']
BL_TOTAL = 64 * len(SNG_CHIPS)
BL_COV = sum(1 for k in BL if k[0] in SNG_CHIPS)
import statistics as _st
_chip_med = {v: _st.median([x[1] for k, x in BL.items() if k[0] == v]) for v in SNG_CHIPS if any(k[0] == v for k in BL)}
CHIP_MED_MIN, CHIP_MED_MAX = min(_chip_med.values()), max(_chip_med.values())
try:
    PEAK_ROWS = [r for r in csv.DictReader(open(SNG_DIR + 'baseline_vs_peak.csv')) if float(r['baseline_adc']) <= 100]
except OSError:
    PEAK_ROWS = []


def _spearman(x, y):
    def rk(v):
        o = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(o):
            j = i
            while j + 1 < len(o) and v[o[j + 1]] == v[o[i]]:
                j += 1
            for k in range(i, j + 1):
                r[o[k]] = (i + j) / 2
            i = j + 1
        return r
    a, b = rk(x), rk(y)
    ma, mb = _st.mean(a), _st.mean(b)
    sa = sum((p - ma) ** 2 for p in a)
    sb = sum((q - mb) ** 2 for q in b)
    return sum((p - ma) * (q - mb) for p, q in zip(a, b)) / (sa * sb) ** 0.5 if sa and sb else float('nan')


PEAK_N = len(PEAK_ROWS)
PEAK_R = _spearman([float(r['dx_adc_minus_chip_median']) for r in PEAK_ROWS],
                   [float(r['dy_peak_d_minus_chip_median']) for r in PEAK_ROWS]) if PEAK_N > 3 else float('nan')


def fig_baseline_heat():
    cw, ch_, left, top = 9.2, 19, 92, 62
    W = int(left + 64 * cw + 14)
    H = top + ch_ * len(SNG_CHIPS) + 70
    def cls(a):
        return 1 if a < 40 else 2 if a < 48 else 3 if a < 56 else 4 if a < 64 else 5
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Measured ADC baseline of every channel">']
    o.append(f'<text class="t" x="14" y="18" font-size="13" font-weight="600">Measured baseline (ADC) of every channel, from neighbour read-out</text>')
    o.append(f'<text class="m" x="14" y="34" font-size="11">one row per chip, one cell per channel 0-63; the colour is the median ADC of its neighbour read-outs</text>')
    o.append('<defs><pattern id="hatch3" width="6" height="6" patternUnits="userSpaceOnUse" patternTransform="rotate(45)"><line x1="0" y1="0" x2="0" y2="6" stroke="var(--fg-faint)" stroke-width="1.6"/></pattern></defs>')
    for c in range(0, 64, 8):
        o.append(f'<text class="m" x="{left + c * cw + cw / 2}" y="{top - 8}" font-size="10" text-anchor="middle">{c}</text>')
    o.append(f'<text class="m" x="14" y="{top - 8}" font-size="10">channel →</text>')
    for i, v in enumerate(SNG_CHIPS):
        y = top + i * ch_
        o.append(f'<text class="t" x="{left - 8}" y="{y + 14}" font-size="12" text-anchor="end">VMM {v}</text>')
        for c in range(64):
            x = left + c * cw
            if (v, c) in BL:
                n, med, sd = BL[(v, c)]
                o.append(f'<g><title>VMM {v} ch {c}: baseline ADC {med:g} (spread {sd:.1f}, {n} read-outs)</title>'
                         f'<rect x="{x:.1f}" y="{y + 1}" width="{cw - 1:.1f}" height="{ch_ - 2}" rx="1.5" class="b{cls(med)}"/></g>')
            else:
                o.append(f'<rect x="{x:.1f}" y="{y + 1}" width="{cw - 1:.1f}" height="{ch_ - 2}" rx="1.5" fill="url(#hatch3)" opacity="0.55"><title>VMM {v} ch {c}: not covered</title></rect>')
    ly = top + ch_ * len(SNG_CHIPS) + 22
    labs = ['&lt; 40', '40 – 47', '48 – 55', '56 – 63', '≥ 64']
    o.append(f'<text class="m" x="14" y="{ly + 11}" font-size="11">ADC</text>')
    for b, lab in enumerate(labs):
        x = left + b * 86
        o.append(f'<rect x="{x}" y="{ly}" width="18" height="13" rx="2" class="b{b + 1}"/><text class="m" x="{x + 24}" y="{ly + 11}" font-size="11">{lab}</text>')
    o.append(f'<rect x="{left + 5 * 86 + 6}" y="{ly}" width="18" height="13" rx="2" fill="url(#hatch3)" opacity="0.55"/><text class="m" x="{left + 5 * 86 + 30}" y="{ly + 11}" font-size="11">not covered yet</text>')
    o.append(f'<text class="f" x="14" y="{ly + 38}" font-size="10">VMM 4 and 9 are fully masked and have no data. The value used so far for every channel was 38.</text>')
    o.append('</svg>')
    return '\n'.join(o)


def fig_baseline_range():
    W, l, r_, t, rh = 620, 96, 80, 64, 22
    H = t + rh * len(SNG_CHIPS) + 44
    lo, hi = 20, 90
    pw = W - l - r_
    def X(a): return l + pw * (min(max(a, lo), hi) - lo) / (hi - lo)
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Range of measured baselines per chip">']
    o.append('<text class="t" x="14" y="18" font-size="13" font-weight="600">Baseline ADC of the covered channels, chip by chip</text>')
    o.append('<text class="m" x="14" y="34" font-size="11">bar: lowest to highest channel; dot: median of the chip; the dashed line is the 38 used so far for every channel</text>')
    for a in range(20, 91, 10):
        o.append(f'<line x1="{X(a)}" x2="{X(a)}" y1="{t - 6}" y2="{t + rh * len(SNG_CHIPS)}" stroke="var(--border)"/>'
                 f'<text class="m" x="{X(a)}" y="{t + rh * len(SNG_CHIPS) + 15}" font-size="11" text-anchor="middle">{a}</text>')
    o.append(f'<line x1="{X(38)}" x2="{X(38)}" y1="{t - 6}" y2="{t + rh * len(SNG_CHIPS)}" stroke="var(--alert)" stroke-width="2" stroke-dasharray="6 4"/>'
             f'<text class="t" x="{X(38)}" y="{t - 10}" font-size="11" text-anchor="middle" font-weight="600">38</text>')
    for i, v in enumerate(SNG_CHIPS):
        y = t + i * rh + 11
        vals = sorted(x[1] for k, x in BL.items() if k[0] == v)
        o.append(f'<text class="t" x="{l - 10}" y="{y + 4}" font-size="12" text-anchor="end">VMM {v}</text>')
        if not vals:
            o.append(f'<text class="m" x="{l + 6}" y="{y + 4}" font-size="11">no covered channel</text>')
            continue
        mn, mx, md = vals[0], vals[-1], _st.median(vals)
        o.append(f'<g><title>VMM {v}: {len(vals)} channels covered, ADC {mn:g} to {mx:g}, median {md:g}</title>'
                 f'<line x1="{X(mn)}" x2="{X(mx)}" y1="{y}" y2="{y}" stroke="var(--fg-faint)" stroke-width="5" stroke-linecap="round"/>'
                 f'<circle cx="{X(md)}" cy="{y}" r="5.5" fill="var(--s1)" stroke="var(--surface)" stroke-width="2"/></g>')
        o.append(f'<text class="f" x="{W - r_ + 8}" y="{y + 4}" font-size="10">{len(vals)}/64' + (f' (max {mx:g})' if mx > hi else '') + '</text>')
    o.append('</svg>')
    return '\n'.join(o)


def fig_coverage_table():
    cov = BL_JSON['covered_after_each_file']
    cw, ch_, left, top = 84, 22, 100, 62
    W = left + cw * len(cov) + 20
    H = top + ch_ * len(SNG_CHIPS) + 52
    def b(n): return 0 if n < 8 else 1 if n < 16 else 2 if n < 32 else 3 if n < 48 else 4
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Channels covered after each capture">']
    o.append('<text class="t" x="14" y="18" font-size="13" font-weight="600">How many of the 64 channels have a baseline, as captures are added</text>')
    o.append('<text class="m" x="14" y="34" font-size="11">counts are cumulative; a channel counts once it has 50 neighbour read-outs</text>')
    names = ['probe', 'step 0', 'step −3', 'step −6', 'step −9']
    for j in range(len(cov)):
        o.append(f'<text class="m" x="{left + j * cw + cw / 2}" y="{top - 10}" font-size="11" text-anchor="middle">{names[j] if j < len(names) else "+" + str(j)}</text>')
    for i, v in enumerate(SNG_CHIPS):
        y = top + i * ch_
        o.append(f'<text class="t" x="{left - 10}" y="{y + 15}" font-size="12" text-anchor="end">VMM {v}</text>')
        for j, c in enumerate(cov):
            n = c.get(str(v), 0)
            x = left + j * cw
            o.append(f'<g><title>VMM {v}, after {names[j] if j < len(names) else j}: {n}/64 channels covered</title>'
                     f'<rect x="{x + 1}" y="{y + 1}" width="{cw - 2}" height="{ch_ - 2}" rx="3" class="b{b(n) + 1}"/>'
                     f'<text class="x{b(n) + 1}" x="{x + cw / 2}" y="{y + 15}" font-size="11" text-anchor="middle">{n}</text></g>')
    ly = top + ch_ * len(SNG_CHIPS) + 20
    for k, lab in enumerate(['&lt; 8', '8 – 15', '16 – 31', '32 – 47', '≥ 48']):
        o.append(f'<rect x="{left + k * 84}" y="{ly}" width="18" height="13" rx="2" class="b{k + 1}"/><text class="m" x="{left + k * 84 + 24}" y="{ly + 11}" font-size="11">{lab}</text>')
    o.append('</svg>')
    return '\n'.join(o)


def fig_baseline_vs_peak():
    if PEAK_N < 4:
        return '<p class="muted">(not enough channels with both a baseline and a clear peak)</p>'
    W, H, l, r_, t, b = 560, 330, 84, 20, 54, 50
    xs = [float(r['dx_adc_minus_chip_median']) for r in PEAK_ROWS]
    ys = [float(r['dy_peak_d_minus_chip_median']) for r in PEAK_ROWS]
    xm = max(8, max(abs(x) for x in xs) + 1)
    ym = max(6, max(abs(y) for y in ys) + 1)
    pw, ph = W - l - r_, H - t - b
    def X(x): return l + pw * (x + xm) / (2 * xm)
    def Y(y): return t + ph * (1 - (y + ym) / (2 * ym))
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Baseline against the position of the noise peak">']
    o.append('<text class="t" x="14" y="18" font-size="13" font-weight="600">Does a higher baseline mean a noise peak at a higher threshold step?</text>')
    o.append(f'<text class="m" x="14" y="34" font-size="11">one dot per channel ({PEAK_N} channels); both axes relative to the chip\'s own median</text>')
    o.append(f'<line x1="{X(0)}" x2="{X(0)}" y1="{t}" y2="{t + ph}" stroke="var(--border)"/><line x1="{l}" x2="{l + pw}" y1="{Y(0)}" y2="{Y(0)}" stroke="var(--border)"/>')
    xt, yt = int(xm // 4 * 4), int(ym // 3 * 3)
    for x in range(-xt, xt + 1, 4):
        o.append(f'<text class="m" x="{X(x)}" y="{t + ph + 16}" font-size="11" text-anchor="middle">{x:+d}</text>')
    for y in range(-yt, yt + 1, 3):
        o.append(f'<text class="m" x="{l - 8}" y="{Y(y) + 4}" font-size="11" text-anchor="end">{y:+d}</text>')
    o.append(f'<text class="m" x="{l + pw / 2}" y="{H - 8}" font-size="12" text-anchor="middle">channel baseline ADC, minus the chip median (measured by neighbour read-out)</text>')
    o.append(f'<text class="m" x="14" y="{t + ph / 2}" font-size="12" text-anchor="middle" transform="rotate(-90 14 {t + ph / 2})">noise-peak step d, minus chip median</text>')
    for r in PEAK_ROWS:
        x, y = float(r['dx_adc_minus_chip_median']), float(r['dy_peak_d_minus_chip_median'])
        o.append(f'<g><title>VMM {r["vmm"]} ch {r["ch"]}: baseline ADC {r["baseline_adc"]}, noise peak at d={r["peak_d"]}</title>'
                 f'<circle cx="{X(x):.1f}" cy="{Y(y):.1f}" r="4.5" fill="var(--s1)" fill-opacity="0.65" stroke="var(--surface)" stroke-width="1.2"/></g>')
    mx_, my_ = _st.mean(xs), _st.mean(ys)
    sxx = sum((x - mx_) ** 2 for x in xs)
    if sxx:
        k_ = sum((x - mx_) * (y - my_) for x, y in zip(xs, ys)) / sxx
        a_ = my_ - k_ * mx_
        o.append(f'<line x1="{X(-xm)}" y1="{Y(a_ - k_ * xm)}" x2="{X(xm)}" y2="{Y(a_ + k_ * xm)}" stroke="var(--alert)" stroke-width="2.2" stroke-dasharray="6 4"/>')
        o.append(f'<text class="t" x="{l + pw - 6}" y="{t + 16}" font-size="12" text-anchor="end" font-weight="600">rank correlation {PEAK_R:+.2f}</text>')
    o.append('</svg>')
    return '\n'.join(o)


def cartoon_neighbour():
    W, H = 640, 300
    n = 11
    random.seed(21)
    base = [46 + random.gauss(0, 4) for _ in range(n)]
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="How neighbour read-out gives a baseline">']
    o.append('<text class="t" x="10" y="18" font-size="13" font-weight="600">Neighbour read-out: a firing channel makes its neighbours report their baseline (cartoon)</text>')
    x0, bw, gap, ybase = 40, 34, 14, 210
    fire = 5
    for k in range(n):
        x = x0 + k * (bw + gap)
        h = base[k] * 1.6
        col = 'var(--s2)' if k == fire else ('var(--s1)' if abs(k - fire) == 1 else 'var(--fg-faint)')
        extra = 0
        if k == fire:
            extra = 70
            o.append(f'<rect x="{x}" y="{ybase - h - extra}" width="{bw}" height="{h + extra}" rx="4" fill="{col}"/>')
            o.append(f'<text class="t" x="{x + bw / 2}" y="{ybase - h - extra - 8}" font-size="11" text-anchor="middle" font-weight="600">fires: pulse above threshold</text>')
        else:
            o.append(f'<rect x="{x}" y="{ybase - h}" width="{bw}" height="{h}" rx="4" fill="{col}" opacity="{1 if abs(k - fire) == 1 else 0.5}"/>')
        o.append(f'<text class="m" x="{x + bw / 2}" y="{ybase + 15}" font-size="10" text-anchor="middle">ch {20 + k}</text>')
    thr = ybase - 46 * 1.6 - 36
    o.append(f'<line x1="{x0 - 10}" x2="{x0 + n * (bw + gap)}" y1="{thr}" y2="{thr}" stroke="var(--alert)" stroke-width="2.2" stroke-dasharray="6 4"/>'
             f'<text class="m" x="{x0 - 10}" y="{thr - 6}" font-size="11">threshold</text>')
    for d in (-1, 1):
        x = x0 + (fire + d) * (bw + gap) + bw / 2
        o.append(f'<path d="M{x} {ybase - base[fire + d] * 1.6 - 12} v-24" stroke="var(--s1)" stroke-width="2.4" marker-end="url(#arr)"/>'
                 f'<text class="t" x="{x}" y="{ybase - base[fire + d] * 1.6 - 44}" font-size="10" text-anchor="middle">read out too</text>')
    o.append('<defs><marker id="arr" markerWidth="8" markerHeight="8" refX="4" refY="4" orient="auto"><path d="M0 0 L8 4 L0 8 z" fill="var(--s1)"/></marker></defs>')
    o.append('<text class="t" x="10" y="246" font-size="12">Each hit comes with a flag: 1 = this channel itself went over threshold; 0 = it was only read because a neighbour fired.</text>')
    o.append('<text class="m" x="10" y="264" font-size="11">The ADC of a flag-0 hit is the channel at rest, with no pulse and no threshold in the way: a direct measurement of its baseline.</text>')
    o.append('<text class="m" x="10" y="280" font-size="11">Switched on by one register setting, sng = 1, on every chip. Cost: about three times more data per firing channel.</text>')
    o.append('</svg>')
    return '\n'.join(o)


def cartoon_ladder():
    W, H = 640, 360
    random.seed(8)
    n = 16
    base = [random.gauss(0, 1) for _ in range(n)]
    steps = [1.6, 0.4, -0.9]
    covered = set()
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="A ladder of thresholds covers all channels">']
    o.append('<text class="t" x="10" y="18" font-size="13" font-weight="600">Why a ladder of thresholds: each step makes different channels fire (cartoon)</text>')
    x0, bw, gap = 214, 18, 8
    for p, th in enumerate(steps):
        y0 = 40 + p * 98
        yb = y0 + 74
        o.append(f'<text class="t" x="10" y="{y0 + 30}" font-size="12" font-weight="600">threshold step {p + 1}</text>')
        fires = {k for k in range(n) if abs(base[k] - th) < 0.55}
        nbrs = {k + d for k in fires for d in (-1, 1) if 0 <= k + d < n}
        for k in range(n):
            x = x0 + k * (bw + gap)
            h = 34 + base[k] * 9
            col = 'var(--s2)' if k in fires else ('var(--s1)' if k in nbrs else 'var(--fg-faint)')
            o.append(f'<rect x="{x}" y="{yb - h:.1f}" width="{bw}" height="{h:.1f}" rx="3" fill="{col}" opacity="{1 if (k in fires or k in nbrs) else 0.45}"/>')
        ty = yb - (34 + th * 9)
        o.append(f'<line x1="{x0 - 8}" x2="{x0 + n * (bw + gap)}" y1="{ty:.1f}" y2="{ty:.1f}" stroke="var(--alert)" stroke-width="2" stroke-dasharray="6 4"/>')
        covered |= fires | nbrs
        o.append(f'<text class="m" x="10" y="{y0 + 46}" font-size="11">fire: {len(fires)}; neighbours: {len(nbrs - fires)}</text>')
        o.append(f'<text class="t" x="10" y="{y0 + 62}" font-size="11" font-weight="600">covered so far: {len(covered)}/{n}</text>')
    o.append('<text class="m" x="10" y="%d" font-size="11">Bars: baselines of 16 channels of one chip (they differ). Red dashed: the chip\'s threshold.</text>' % (H - 38))
    o.append('<text class="m" x="10" y="%d" font-size="11">Orange: fires. Blue: read out as a neighbour. One step leaves most channels unmeasured; a few steps cover (nearly) all.</text>' % (H - 22))
    o.append('</svg>')
    return '\n'.join(o)


BL_MEDS = [x[1] for k, x in BL.items() if k[0] in SNG_CHIPS]
BL_MEDIAN = _st.median(BL_MEDS)
N_LE38 = sum(1 for m in BL_MEDS if m <= 38)
_nb_frac = [v['nb_fraction'] for k, v in BL_JSON['per_vmm'].items() if int(k) in SNG_CHIPS]
NB_LO, NB_HI = 100 * min(_nb_frac), 100 * max(_nb_frac)
COV_STEPS = [sum(v for k, v in c.items() if int(k) in SNG_CHIPS) for c in BL_JSON['covered_after_each_file']]
PEAK_ROWS_ALL = {}
for _mh in (30, 10, 5):
    try:
        _rows = [r for r in csv.DictReader(open(SNG_DIR + f'baseline_vs_peak_mh{_mh}.csv')) if float(r['baseline_adc']) <= 100]
        PEAK_ROWS_ALL[_mh] = (len(_rows), _spearman([float(r['dx_adc_minus_chip_median']) for r in _rows],
                                                    [float(r['dy_peak_d_minus_chip_median']) for r in _rows]) if len(_rows) > 3 else float('nan'))
    except OSError:
        pass


def cartoon_languages():
    W, H = 640, 190
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Three units for the same level">']
    o.append('<text class="t" x="10" y="18" font-size="13" font-weight="600">Three ways to say the same level: what we set, what is real, what the data shows</text>')
    boxes = [(10, 'threshold setting', 'sdt, in "DAC" steps', 'what we write to the chip', 'known (we choose it)', 'var(--s1)'),
             (230, 'millivolts', 'the real voltage level', 'via the chip monitor', 'not measured yet', 'var(--s2)'),
             (450, 'ADC counts', 'what appears in the data', 'e.g. from neighbour read-out', 'measured now', 'var(--s1)')]
    for x, ttl, l1, l2, st_, col in boxes:
        done = 'measured now' in st_ or 'known' in st_
        o.append(f'<rect x="{x}" y="36" width="180" height="118" rx="10" fill="{col}" opacity="0.14"/><rect x="{x}" y="36" width="180" height="118" rx="10" fill="none" stroke="{col}" stroke-width="2.4" stroke-dasharray="{"0" if done else "7 5"}"/>')
        o.append(f'<text class="t" x="{x + 90}" y="62" font-size="13" text-anchor="middle" font-weight="600">{ttl}</text>'
                 f'<text class="m" x="{x + 90}" y="82" font-size="11" text-anchor="middle">{l1}</text>'
                 f'<text class="m" x="{x + 90}" y="98" font-size="11" text-anchor="middle">{l2}</text>'
                 f'<text class="t" x="{x + 90}" y="130" font-size="11" text-anchor="middle" font-weight="600">{st_}</text>')
    for x in (192, 412):
        o.append(f'<path d="M{x} 95 h26" stroke="var(--fg-faint)" stroke-width="3" marker-end="url(#arr2)"/>')
    o.append('<defs><marker id="arr2" markerWidth="8" markerHeight="8" refX="4" refY="4" orient="auto"><path d="M0 0 L8 4 L0 8 z" fill="var(--fg-faint)"/></marker></defs>')
    o.append('<text class="m" x="10" y="176" font-size="11">Both ends are in hand; the middle link is the next measurement. Solid frame: have it. Dashed: still to do.</text>')
    o.append('</svg>')
    return '\n'.join(o)



def fig_baseline_panels():
    cols, pw, ph = 4, 156, 150
    rows_ = (len(SNG_CHIPS) + cols - 1) // cols
    W = cols * pw + 24
    H = 62 + rows_ * ph + 30
    lo, hi = 30, 72
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Baseline of every channel, one panel per chip">']
    o.append('<text class="t" x="12" y="18" font-size="13" font-weight="600">Baseline of every channel, one panel per chip</text>')
    o.append('<text class="m" x="12" y="32" font-size="11">x: channel 0 to 63; y: ADC baseline (neighbour read-out), shown from 30 to 72.</text>')
    o.append('<text class="m" x="12" y="45" font-size="11">Dashed orange: the 38 used so far. Grey line: the chip median. Ticks at the bottom: channels not covered.</text>')
    for i, v in enumerate(SNG_CHIPS):
        r_, c_ = divmod(i, cols)
        x0, y0 = 12 + c_ * pw, 54 + r_ * ph
        l, t, w, h = x0 + 24, y0 + 16, pw - 34, ph - 40
        def X(ch): return l + w * ch / 63
        def Y(a): return t + h * (1 - (min(max(a, lo), hi) - lo) / (hi - lo))
        pts = [(ch, BL[(v, ch)]) for ch in range(64) if (v, ch) in BL]
        o.append(f'<text class="t" x="{x0 + 24}" y="{y0 + 11}" font-size="12" font-weight="600">VMM {v}</text>'
                 f'<text class="f" x="{x0 + pw - 10}" y="{y0 + 11}" font-size="10" text-anchor="end">{len(pts)}/64</text>')
        for a in (40, 50, 60, 70):
            o.append(f'<line x1="{l}" x2="{l + w}" y1="{Y(a):.1f}" y2="{Y(a):.1f}" stroke="var(--border)"/>'
                     f'<text class="f" x="{l - 4}" y="{Y(a) + 3:.1f}" font-size="9" text-anchor="end">{a}</text>')
        o.append(f'<line x1="{l}" x2="{l + w}" y1="{Y(38):.1f}" y2="{Y(38):.1f}" stroke="var(--alert)" stroke-width="1.4" stroke-dasharray="4 3"/>')
        if pts:
            md = _st.median([x[1][1] for x in pts])
            o.append(f'<line x1="{l}" x2="{l + w}" y1="{Y(md):.1f}" y2="{Y(md):.1f}" stroke="var(--fg-faint)" stroke-width="1.4"/>')
        for ch in range(64):
            if (v, ch) not in BL:
                o.append(f'<line x1="{X(ch):.1f}" x2="{X(ch):.1f}" y1="{t + h + 2}" y2="{t + h + 6}" stroke="var(--fg-faint)" stroke-width="1"/>')
        for ch, (n, med, sd) in pts:
            clip = med > hi or med < lo
            o.append(f'<g><title>VMM {v} ch {ch}: baseline ADC {med:g} (spread {sd:.1f}, {n} read-outs)</title>'
                     f'<circle cx="{X(ch):.1f}" cy="{Y(med):.1f}" r="{2.4 if clip else 1.9}" fill="{"var(--s2)" if clip else "var(--s1)"}"/></g>')
        for tx in (0, 32, 63):
            o.append(f'<text class="f" x="{X(tx):.1f}" y="{t + h + 18}" font-size="9" text-anchor="middle">{tx}</text>')
    o.append('</svg>')
    return '\n'.join(o)


def baseline_tables():
    out = []
    for v in SNG_CHIPS:
        n_cov = sum(1 for ch in range(64) if (v, ch) in BL)
        rows = []
        for r_ in range(8):
            cells = []
            for c_ in range(8):
                ch = r_ * 8 + c_
                if (v, ch) in BL:
                    n, med, sd = BL[(v, ch)]
                    cells.append(f'<td title="ch {ch}: {n} read-outs, spread {sd:.1f}">{med:g}</td>')
                else:
                    cells.append('<td class="na" title="not covered">&ndash;</td>')
            rows.append(f'<tr><th>ch {r_ * 8}&ndash;{r_ * 8 + 7}</th>' + ''.join(cells) + '</tr>')
        out.append(f'<details><summary>VMM {v} &mdash; {n_cov} of 64 channels measured</summary>'
                   f'<table class="bltab">{"".join(rows)}</table></details>')
    return '\n'.join(out)


# ------------------------------------------------------------------ blocks A + C: VMM 4, 7, 9 unmasked (HV off)
BURST = json.load(open(LOCAL + 'analysis_n479_hvoff/burst_anatomy_A_C.json'))
BURST_A, BURST_C = {}, {}
for _k, _v in BURST.items():
    _m = _re.search(r'_([ns])479hvoff_d(\d+)$', _k)
    if not _m:
        continue
    _tgt = BURST_A if _m.group(1) == 'n' else BURST_C
    for _c, _r in _v.items():
        _top = list(_r['channel_share_pct'].items())[0]
        _tgt.setdefault(int(_m.group(2)), {})[int(_c)] = (_r['hits'] / max(_r['span_s'], 1e-9) / 1e3, int(_top[0]), _top[1],
                                                         [int(x) for x in list(_r['channel_share_pct'])[:3]])
BURST_D = sorted(set(BURST_A) | set(BURST_C), reverse=True)


def fig_burst_ladder():
    W, H = 640, 356
    t, b = 104, 52
    ph = H - t - b
    lo, hi = math.log10(0.01), math.log10(6000)
    pw = (W - 20) / 3
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Rate of VMM 4, 7 and 9 against threshold step">']
    o.append('<text class="t" x="10" y="18" font-size="13" font-weight="600">VMM 4, 7 and 9 unmasked, voltage off: rate against threshold step</text>')
    o.append('<text class="m" x="10" y="33" font-size="11">kHz per chip, log scale; first 300 MB of each 2-minute capture.</text>')
    o.append('<text class="m" x="10" y="47" font-size="11">Dots on the bottom edge: no hits, or fewer than 0.01 kHz.</text>')
    o.append(f'<circle cx="16" cy="66" r="5" fill="var(--s1)"/><text class="m" x="26" y="70" font-size="11">normal read-out (sng 0)</text>'
             f'<circle cx="190" cy="66" r="6.5" fill="none" stroke="var(--s2)" stroke-width="2.6"/><text class="m" x="203" y="70" font-size="11">neighbour read-out (sng 1)</text>')
    o.append('<text class="f" x="10" y="86" font-size="10">x: threshold step d (left: very high threshold, right: low). Hover a point for the hot channel.</text>')
    n = len(BURST_D)
    for i, v in enumerate((4, 7, 9)):
        x0 = 10 + i * pw
        l, r_ = x0 + 34, x0 + pw - 6
        def X(j): return l + (r_ - l) * (j + 0.5) / n
        def Y(kh): return t + ph * (1 - (math.log10(max(kh, 0.01)) - lo) / (hi - lo))
        o.append(f'<text class="t" x="{x0 + pw / 2}" y="{t - 6}" font-size="12" font-weight="600" text-anchor="middle">VMM {v}</text>')
        for e in (-2, 0, 2):
            y = Y(10 ** e)
            o.append(f'<line x1="{l}" x2="{r_}" y1="{y:.1f}" y2="{y:.1f}" stroke="var(--border)"/>'
                     f'<text class="f" x="{l - 4}" y="{y + 3:.1f}" font-size="9" text-anchor="end">{10 ** e:g}</text>')
        for j, d in enumerate(BURST_D):
            if j % 2 == 0:
                o.append(f'<text class="f" x="{X(j):.1f}" y="{t + ph + 14}" font-size="9" text-anchor="middle">{d}</text>')
        for series, dx, ring in ((BURST_A, -2.5, False), (BURST_C, 2.5, True)):
            for j, d in enumerate(BURST_D):
                rec = series.get(d, {}).get(v)
                if d not in series:
                    continue
                kh = rec[0] if rec else 0.0
                tip = (f'VMM {v}, d={d}, {"sng 1" if ring else "sng 0"}: {kh:.2f} kHz' + (f', hottest channel {rec[1]} ({rec[2]:.0f} %)' if rec else ', no hits'))
                cx, cy = X(j) + dx, Y(kh)
                o.append(f'<g><title>{tip}</title>' + (f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="4.4" fill="none" stroke="var(--s2)" stroke-width="2.2"/>' if ring else
                         f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="3.4" fill="var(--s1)"/>') + '</g>')
    o.append(f'<text class="m" x="{W / 2}" y="{H - 8}" font-size="11" text-anchor="middle">threshold step d, from very high (left) to low (right)</text>')
    o.append('</svg>')
    return '\n'.join(o)


def burst_table():
    rows = []
    for d in BURST_D:
        cells = []
        for v in (4, 7, 9):
            for series in (BURST_A, BURST_C):
                rec = series.get(d, {}).get(v) if d in series else None
                if d not in series:
                    cells.append('<td class="na">&ndash;</td>')
                elif rec is None:
                    cells.append('<td>no hits</td>')
                else:
                    kh = rec[0]
                    txt = '&lt; 0.01' if kh < 0.01 else (f'{kh:.2f}' if kh < 10 else (f'{kh:.0f}' if kh < 10000 else f'{kh / 1000:.1f}k'))
                    cells.append(f'<td title="top channels: {", ".join(str(c) for c in rec[3])}">{txt} kHz &middot; ch {rec[1]} ({rec[2]:.0f}%)</td>')
        rows.append(f'<tr><th>d = {d}</th>' + ''.join(cells) + '</tr>')
    head = '<tr><th></th>' + ''.join(f'<th colspan="2">VMM {v}</th>' for v in (4, 7, 9)) + '</tr><tr><th></th>' + '<th>sng 0</th><th>sng 1</th>' * 3 + '</tr>'
    return (f'<details><summary>Rate and hottest channel for every threshold step</summary><div style="overflow-x:auto">'
            f'<table class="bltab">{head}{"".join(rows)}</table></div></details>')


def cartoon_two_kinds():
    W, H = 640, 308
    random.seed(5)
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Two kinds of noise on the burst chips">']
    o.append('<text class="t" x="10" y="18" font-size="13" font-weight="600">Two kinds of noise on VMM 4, 7 and 9 (cartoon, not measured)</text>')
    for p, (x0, ttl, sub) in enumerate(((14, 'A: smooth onset (all chips have it)', 'the rate rises steadily as the threshold goes down'),
                                          (332, 'B: bursts (these three chips only)', 'the rate jumps up and down with each new configuration'))):
        o.append(f'<rect x="{x0}" y="32" width="296" height="206" rx="8" fill="var(--bg-elev)" stroke="var(--border)"/>')
        o.append(f'<text class="t" x="{x0 + 148}" y="52" font-size="12" font-weight="600" text-anchor="middle">{ttl}</text>'
                 f'<text class="m" x="{x0 + 148}" y="67" font-size="10" text-anchor="middle">{sub}</text>')
        l, r_, tp, bt = x0 + 26, x0 + 284, 84, 208
        o.append(f'<line x1="{l}" x2="{r_}" y1="{bt}" y2="{bt}" stroke="var(--border)"/>')
        if p == 0:
            pts = [(l + (r_ - l) * k / 11, bt - (bt - tp) * (math.exp(0.38 * k) - 1) / (math.exp(0.38 * 11) - 1) * 0.95) for k in range(12)]
            o.append('<polyline fill="none" stroke="var(--s1)" stroke-width="2.6" points="' + ' '.join(f'{x:.1f},{y:.1f}' for x, y in pts) + '"/>')
            for x, y in pts:
                o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="3.6" fill="var(--s1)"/>')
        else:
            heights = [0.0, 0.1, 0.05, 0.9, 0.0, 0.0, 0.55, 0.0, 0.12, 0.8, 0.0, 0.3]
            cols = ['var(--s2)', 'var(--s2)', 'var(--s3)', 'var(--s2)', 'var(--s2)', 'var(--s2)', 'var(--s3)', 'var(--s2)', 'var(--s2)', 'var(--s1)', 'var(--s2)', 'var(--s1)']
            for k, h in enumerate(heights):
                x = l + (r_ - l) * k / 11
                o.append(f'<rect x="{x - 8:.1f}" y="{bt - (bt - tp) * max(h, 0.015):.1f}" width="16" height="{(bt - tp) * max(h, 0.015):.1f}" rx="3" fill="{cols[k]}" opacity="0.9"/>')
        o.append(f'<text class="f" x="{l}" y="{bt + 14}" font-size="10">high threshold</text><text class="f" x="{r_}" y="{bt + 14}" font-size="10" text-anchor="end">low threshold</text>')
    o.append('<text class="m" x="332" y="252" font-size="11">colour = which channel is hot (changes each time)</text>')
    o.append('<text class="m" x="14" y="268" font-size="11">A is ordinary noise crossing the threshold. B appears even at very high thresholds (far above any baseline),</text>')
    o.append('<text class="m" x="14" y="282" font-size="11">so it is not noise crossing a threshold: it looks like a state a channel gets stuck in when configured.</text>')
    o.append('<text class="m" x="14" y="296" font-size="11">Whether the same configuration gives the same state again has still to be tested.</text>')
    o.append('</svg>')
    return '\n'.join(o)



# ------------------------------------------------------------------ 24 Oct: HV ON health check, sng 1 with cosmics, VMM 4/7/9 ladder with HV on
HEALTH = json.load(open(LOCAL + 'analysis_health/run_57_adc_spectra.json'))


def _load_bl(path, minnb=50):
    d = {}
    for r in csv.DictReader(open(path)):
        if r['adc_nb_median'] != '' and int(r['n_nb']) >= minnb:
            d[(int(r['vmm']), int(r['ch']))] = float(r['adc_nb_median'])
    return d


BL_ON = _load_bl(LOCAL + 'analysis_hvon/sng/neighbor_baseline.csv')
BL_OFFALL = _load_bl(LOCAL + 'analysis_sng1_hvoff/neighbor_baseline.csv')
for _k, _v in _load_bl(LOCAL + 'analysis_n479_hvoff/sng/neighbor_baseline.csv').items():
    BL_OFFALL.setdefault(_k, _v)
ALL_CHIPS = list(range(2, 20))
NORMAL_CHIPS = [v for v in ALL_CHIPS if v not in (4, 7, 9)]
ON_COV, ON_TOTAL = len(BL_ON), 64 * len(ALL_CHIPS)
OFF_COV = len(BL_OFFALL)
ON_MEDIAN = _st.median(BL_ON.values())
ON_LE38 = sum(1 for x in BL_ON.values() if x <= 38)
_core = [k for k in BL_ON if k in BL_OFFALL and k[0] in NORMAL_CHIPS]
_dd = [BL_ON[k] - BL_OFFALL[k] for k in _core]
CORE_N, CORE_MED = len(_core), _st.median(_dd)
CORE_MAD = _st.median([abs(x - CORE_MED) for x in _dd])
CORE_W3 = 100 * sum(abs(x) <= 3 for x in _dd) / len(_dd)
CORE_W5 = 100 * sum(abs(x) <= 5 for x in _dd) / len(_dd)
FULL64 = sum(1 for v in ALL_CHIPS if sum(1 for ch in range(64) if (v, ch) in BL_ON) == 64)
OTHER_HIGH = sorted((k, x) for k, x in BL_ON.items() if k[0] in NORMAL_CHIPS and x > 80)
_HOT_UNION = {4: [34, 35, 36, 37, 38, 39, 40, 41, 43, 45, 49], 7: [9, 10, 11, 22, 23, 24, 25, 26, 30, 35, 47, 48],
              9: [15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 33, 36, 41, 42, 44, 45, 63]}
_MASKED = {4: [34, 35, 36, 37, 38, 39], 7: [22, 23, 24, 25], 9: [15, 16, 17, 18, 19, 20, 21, 24, 25, 26, 27, 33, 44]}
BLOCK = {}
for _v in (4, 7, 9):
    _cov = {ch: BL_ON[(_v, ch)] for ch in range(64) if (_v, ch) in BL_ON and ch not in _MASKED[_v]}
    _hi = sorted(ch for ch, x in _cov.items() if x > 80)
    _start = None
    _chs = sorted(_cov)
    for _i in range(len(_chs) - 4):
        if all(_cov[c] > 80 for c in _chs[_i:_i + 5]):       # first run of 5 measured channels that are all high
            _start = _chs[_i]
            break
    BLOCK[_v] = {'n_cov': len(_cov), 'n_high': len(_hi), 'start': _start, 'vmax': max(_cov.values()),
                 'median_high': _st.median([_cov[c] for c in _hi]) if _hi else None}


def _table_ladder(path, tag):
    out = {}
    for k, v in json.load(open(path)).items():
        m = _re.search(r'_%s_d(\d+)$' % tag, k)
        if not m:
            continue
        for c, r in v.items():
            top = list(r['channel_share_pct'].items())[0]
            out.setdefault(int(m.group(1)), {})[int(c)] = (r['hits'] / max(r['span_s'], 1e-9) / 1e3, int(top[0]), top[1])
    return out


LAD_ON = _table_ladder(LOCAL + 'analysis_hvon/burst_anatomy_hvon.json', 'n479hvon')
LAD_OFF = _table_ladder(LOCAL + 'analysis_n479_hvoff/burst_anatomy_A_C.json', 'n479hvoff')
LAD_D = sorted(set(LAD_ON) & set(LAD_OFF), reverse=True)


def _cmp(a, b):
    ra = a[0] if a else 0.0
    rb = b[0] if b else 0.0
    if ra < 1 and rb < 1:
        return 'quiet'
    if ra >= 1 and rb >= 1 and 0.5 <= ra / rb <= 2 and a[1] == b[1]:
        return 'same'
    return 'differ'


LAD_CMP = {c: {'same': 0, 'quiet': 0, 'differ': 0} for c in (4, 7, 9)}
for _d in LAD_D:
    for _c in (4, 7, 9):
        LAD_CMP[_c][_cmp(LAD_ON[_d].get(_c), LAD_OFF[_d].get(_c))] += 1


def fig_health():
    W, rh = 640, 19
    l, r_, t = 130, 70, 64
    rows = [(n, v) for n, rng in (('det3 (IN)', range(2, 8)), ('det2 (MID)', range(8, 14)), ('det1 (OUT)', range(14, 20))) for v in rng]
    H = t + rh * len(rows) + 46 + 2 * 8
    pw = W - l - r_
    def X(a): return l + pw * min(a, 1023) / 1023
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="ADC spectrum of the hits of every chip, voltage on">']
    o.append('<text class="t" x="12" y="18" font-size="13" font-weight="600">Health check with the voltage on: pulse-height spectrum of every chip</text>')
    o.append('<text class="m" x="12" y="34" font-size="11">bar: 10th to 90th percentile of the ADC of the hits above threshold; dot: median; 2 minutes at the safe setting.</text>')
    o.append('<text class="m" x="12" y="48" font-size="11">Orange tick: the spectrum reaches the ADC limit (1023). Right: hits per second of the chip.</text>')
    for a in (0, 200, 400, 600, 800, 1000):
        o.append(f'<line x1="{X(a)}" x2="{X(a)}" y1="{t - 4}" y2="{t + rh * len(rows) + 6 * 0}" stroke="var(--border)"/>'
                 f'<text class="f" x="{X(a)}" y="{t + rh * len(rows) + 14}" font-size="10" text-anchor="middle">{a}</text>')
    prev = None
    for i, (name, v) in enumerate(rows):
        y = t + i * rh + rh / 2
        if name != prev:
            o.append(f'<text class="t" x="12" y="{y + 4}" font-size="11" font-weight="600">{name}</text>')
            prev = name
        c = HEALTH['chips'].get(str(v))
        o.append(f'<text class="m" x="{l - 8}" y="{y + 4}" font-size="11" text-anchor="end">VMM {v}</text>')
        if not c or not c.get('hits'):
            o.append(f'<text class="f" x="{l + 6}" y="{y + 4}" font-size="10">no hits</text>')
            continue
        o.append(f'<g><title>VMM {v}: {c["rate_hz"]:g} hits/s, median {c["median"]}, 10-90 % range {c["p10"]} to {c["p90"]}</title>'
                 f'<rect x="{X(c["p10"]):.1f}" y="{y - 5}" width="{max(X(c["p90"]) - X(c["p10"]), 2):.1f}" height="10" rx="3" fill="var(--s1)" opacity="0.35"/>'
                 f'<circle cx="{X(c["median"]):.1f}" cy="{y}" r="4" fill="var(--s1)"/></g>')
        if c['p90'] >= 1023:
            o.append(f'<line x1="{X(1023):.1f}" x2="{X(1023):.1f}" y1="{y - 7}" y2="{y + 7}" stroke="var(--s2)" stroke-width="3"/>')
        o.append(f'<text class="m" x="{W - 8}" y="{y + 4}" font-size="11" text-anchor="end">{c["rate_hz"]:.0f} Hz</text>')
    o.append(f'<text class="m" x="{l + pw / 2}" y="{H - 6}" font-size="11" text-anchor="middle">ADC counts</text>')
    o.append('</svg>')
    return '\n'.join(o)


def fig_cov_compare():
    W, l, r_, t, rh = 620, 70, 20, 54, 21
    H = t + rh * len(ALL_CHIPS) + 30
    pw = W - l - r_
    def X(n): return l + pw * n / 64
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Channels with a measured baseline per chip, voltage off and on">']
    o.append('<text class="t" x="12" y="18" font-size="13" font-weight="600">Channels with a measured baseline (of 64), voltage off and voltage on</text>')
    o.append('<text class="m" x="12" y="34" font-size="11">neighbour read-out; a channel counts when it has at least 50 read-outs</text>')
    o.append(f'<rect x="{l}" y="40" width="18" height="9" rx="2" fill="var(--fg-faint)"/><text class="m" x="{l + 24}" y="48" font-size="11">voltage off (all earlier runs)</text>'
             f'<rect x="{l + 200}" y="40" width="18" height="9" rx="2" fill="var(--s1)"/><text class="m" x="{l + 224}" y="48" font-size="11">voltage on, with cosmics</text>')
    for n in (0, 16, 32, 48, 64):
        o.append(f'<line x1="{X(n)}" x2="{X(n)}" y1="{t - 2}" y2="{t + rh * len(ALL_CHIPS)}" stroke="var(--border)"/>'
                 f'<text class="f" x="{X(n)}" y="{t + rh * len(ALL_CHIPS) + 14}" font-size="10" text-anchor="middle">{n}</text>')
    for i, v in enumerate(ALL_CHIPS):
        y = t + i * rh
        a = sum(1 for ch in range(64) if (v, ch) in BL_OFFALL)
        b = sum(1 for ch in range(64) if (v, ch) in BL_ON)
        o.append(f'<text class="t" x="{l - 8}" y="{y + 14}" font-size="12" text-anchor="end">VMM {v}</text>')
        o.append(f'<g><title>VMM {v}: {a} channels with the voltage off</title><rect x="{l}" y="{y + 2}" width="{X(a) - l:.1f}" height="7" rx="2" fill="var(--fg-faint)"/></g>'
                 f'<g><title>VMM {v}: {b} channels with the voltage on</title><rect x="{l}" y="{y + 10}" width="{X(b) - l:.1f}" height="7" rx="2" fill="var(--s1)"/></g>')
        o.append(f'<text class="m" x="{X(max(a, b)) + 5:.1f}" y="{y + 14}" font-size="10">{b}</text>')
    o.append('</svg>')
    return '\n'.join(o)


def fig_on_off_scatter():
    pts = [(BL_OFFALL[k], BL_ON[k]) for k in _core if 30 <= BL_OFFALL[k] <= 80 and 30 <= BL_ON[k] <= 80]
    W, H, l, r_, t, b = 520, 400, 64, 20, 44, 50
    pw, ph = W - l - r_, H - t - b
    lo, hi = 30, 80
    def X(a): return l + pw * (a - lo) / (hi - lo)
    def Y(a): return t + ph * (1 - (a - lo) / (hi - lo))
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Baseline with the voltage on against the voltage off">']
    o.append('<text class="t" x="12" y="18" font-size="13" font-weight="600">Does the voltage change a channel&#8217;s baseline? The 15 normal chips</text>')
    o.append(f'<text class="m" x="12" y="34" font-size="11">one dot per channel ({len(pts)} of the {CORE_N} channels fall inside the plotted range 30 to 80); dashed line: &#8220;no change&#8221;</text>')
    for a in range(30, 81, 10):
        o.append(f'<line x1="{X(a)}" x2="{X(a)}" y1="{t}" y2="{t + ph}" stroke="var(--border)"/><line x1="{l}" x2="{l + pw}" y1="{Y(a)}" y2="{Y(a)}" stroke="var(--border)"/>'
                 f'<text class="m" x="{X(a)}" y="{t + ph + 15}" font-size="11" text-anchor="middle">{a}</text><text class="m" x="{l - 8}" y="{Y(a) + 4}" font-size="11" text-anchor="end">{a}</text>')
    o.append(f'<line x1="{X(lo)}" y1="{Y(lo)}" x2="{X(hi)}" y2="{Y(hi)}" stroke="var(--alert)" stroke-width="2" stroke-dasharray="6 4"/>')
    for x, y in pts:
        o.append(f'<circle cx="{X(x):.1f}" cy="{Y(y):.1f}" r="3.6" fill="var(--s1)" fill-opacity="0.35"/>')
    o.append(f'<text class="m" x="{l + pw / 2}" y="{H - 8}" font-size="12" text-anchor="middle">baseline, voltage off (ADC counts)</text>')
    o.append(f'<text class="m" x="16" y="{t + ph / 2}" font-size="12" text-anchor="middle" transform="rotate(-90 16 {t + ph / 2})">baseline, voltage on (ADC counts)</text>')
    o.append('</svg>')
    return '\n'.join(o)


def fig_block_panels():
    chips = (4, 7, 9, 12)
    pw, ph = 156, 190
    W, H = 4 * pw + 24, 70 + ph + 40
    lo, hi = math.log10(20), math.log10(1100)
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Baseline of every channel of VMM 4, 7, 9 and a normal chip">']
    o.append('<text class="t" x="12" y="18" font-size="13" font-weight="600">Baseline of every channel, voltage on: the three burst chips and a normal chip (VMM 12)</text>')
    o.append('<text class="m" x="12" y="34" font-size="11">x: channel 0 to 63; y: ADC baseline on a log scale. Orange: above 80. Grey ticks at the bottom: masked or not measured.</text>')
    o.append('<text class="m" x="12" y="48" font-size="11">Orange triangles on top: channels that were hot (bursting) in at least 3 runs.</text>')
    for i, v in enumerate(chips):
        x0, y0 = 12 + i * pw, 58
        l, t, w, h = x0 + 30, y0 + 22, pw - 40, ph - 30
        def X(ch): return l + w * ch / 63
        def Y(a): return t + h * (1 - (math.log10(min(max(a, 20), 1100)) - lo) / (hi - lo))
        o.append(f'<text class="t" x="{x0 + 30}" y="{y0 + 12}" font-size="12" font-weight="600">VMM {v}</text>')
        for a in (30, 50, 100, 300, 1000):
            o.append(f'<line x1="{l}" x2="{l + w}" y1="{Y(a):.1f}" y2="{Y(a):.1f}" stroke="var(--border)"/><text class="f" x="{l - 4}" y="{Y(a) + 3:.1f}" font-size="9" text-anchor="end">{a}</text>')
        o.append(f'<line x1="{l}" x2="{l + w}" y1="{Y(80):.1f}" y2="{Y(80):.1f}" stroke="var(--alert)" stroke-width="1" stroke-dasharray="3 3"/>')
        for ch in range(64):
            if (v, ch) in BL_ON and BL_ON[(v, ch)] > 0:
                a = BL_ON[(v, ch)]
                o.append(f'<g><title>VMM {v} ch {ch}: baseline ADC {a:g}</title><circle cx="{X(ch):.1f}" cy="{Y(a):.1f}" r="2.2" fill="{"var(--s2)" if a > 80 else "var(--s1)"}"/></g>')
            else:
                o.append(f'<line x1="{X(ch):.1f}" x2="{X(ch):.1f}" y1="{t + h + 2}" y2="{t + h + 6}" stroke="var(--fg-faint)"/>')
        for ch in _MASKED.get(v, []):
            o.append(f'<path d="M{X(ch) - 3:.1f} {t - 8} h6 l-3 6 z" fill="var(--s2)"/>')
        for tx in (0, 32, 63):
            o.append(f'<text class="f" x="{X(tx):.1f}" y="{t + h + 18}" font-size="9" text-anchor="middle">{tx}</text>')
    o.append('</svg>')
    return '\n'.join(o)


def fig_ladder_onoff():
    W, H = 640, 360
    t, b = 104, 52
    ph = H - t - b
    lo, hi = math.log10(0.01), math.log10(6000)
    pw = (W - 20) / 3
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="VMM 4, 7 and 9 unmasked: rate against threshold with the voltage on and off">']
    o.append('<text class="t" x="10" y="18" font-size="13" font-weight="600">VMM 4, 7 and 9 unmasked: voltage on (dot) against voltage off (ring)</text>')
    o.append('<text class="m" x="10" y="33" font-size="11">kHz per chip, log scale; same configurations about an hour apart. Dots on the bottom edge: no hits, or fewer than 0.01 kHz.</text>')
    o.append('<text class="m" x="10" y="47" font-size="11">Where a dot sits inside its ring, the rate is the same with and without the voltage.</text>')
    o.append(f'<circle cx="16" cy="66" r="4.6" fill="var(--s1)"/><text class="m" x="26" y="70" font-size="11">voltage on (HV on, normal read-out)</text>'
             f'<circle cx="260" cy="66" r="6.6" fill="none" stroke="var(--s2)" stroke-width="2.6"/><text class="m" x="273" y="70" font-size="11">voltage off</text>')
    o.append('<text class="f" x="10" y="86" font-size="10">x: threshold step d (left: very high threshold, right: low). Hover for the hot channel.</text>')
    n = len(LAD_D)
    for i, v in enumerate((4, 7, 9)):
        x0 = 10 + i * pw
        l, r_ = x0 + 34, x0 + pw - 6
        def X(j): return l + (r_ - l) * (j + 0.5) / n
        def Y(kh): return t + ph * (1 - (math.log10(max(kh, 0.01)) - lo) / (hi - lo))
        o.append(f'<text class="t" x="{x0 + pw / 2}" y="{t - 6}" font-size="12" font-weight="600" text-anchor="middle">VMM {v}</text>')
        for e in (-2, 0, 2):
            y = Y(10 ** e)
            o.append(f'<line x1="{l}" x2="{r_}" y1="{y:.1f}" y2="{y:.1f}" stroke="var(--border)"/><text class="f" x="{l - 4}" y="{y + 3:.1f}" font-size="9" text-anchor="end">{10 ** e:g}</text>')
        for j, d in enumerate(LAD_D):
            if j % 2 == 0:
                o.append(f'<text class="f" x="{X(j):.1f}" y="{t + ph + 14}" font-size="9" text-anchor="middle">{d}</text>')
        for tbl, ring in ((LAD_OFF, True), (LAD_ON, False)):
            for j, d in enumerate(LAD_D):
                rec = tbl[d].get(v)
                kh = rec[0] if rec else 0.0
                tip = f'VMM {v}, d={d}, voltage {"off" if ring else "on"}: {kh:.3g} kHz' + (f', hottest channel {rec[1]} ({rec[2]:.0f} %)' if rec else ', no hits')
                cx, cy = X(j), Y(kh)
                o.append(f'<g><title>{tip}</title>' + (f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="5.6" fill="none" stroke="var(--s2)" stroke-width="2.2"/>' if ring else
                         f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="3.2" fill="var(--s1)"/>') + '</g>')
    o.append(f'<text class="m" x="{W / 2}" y="{H - 8}" font-size="11" text-anchor="middle">threshold step d, from very high (left) to low (right)</text>')
    o.append('</svg>')
    return '\n'.join(o)


def cartoon_block():
    W, H = 640, 330
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="A block of channels with a high baseline and the bursts at its edge">']
    o.append('<text class="t" x="10" y="18" font-size="13" font-weight="600">What the three burst chips look like (cartoon of the data above)</text>')
    x0, w = 24, 592
    cw = w / 64
    start = 34
    base = 214
    for ch in range(64):
        high = ch >= start
        h = 40 + (ch - start) * 2.2 if high else 36 + (ch % 5)
        h = min(h, 108)
        col = 'var(--s2)' if high else 'var(--s1)'
        o.append(f'<rect x="{x0 + ch * cw + 0.8:.1f}" y="{base - h:.1f}" width="{cw - 1.6:.1f}" height="{h:.1f}" rx="1.5" fill="{col}" opacity="{0.9 if high else 0.75}"/>')
    for ch in (28, 30, 32, 34):
        o.append(f'<path d="M{x0 + ch * cw + cw / 2 - 5:.1f} 150 h10 l-5 9 z" fill="var(--alert)"/>')
    o.append(f'<line x1="{x0}" x2="{x0 + w}" y1="{base}" y2="{base}" stroke="var(--border)"/>')
    o.append(f'<text class="m" x="{x0}" y="{base + 16}" font-size="11">channel 0</text><text class="m" x="{x0 + w}" y="{base + 16}" font-size="11" text-anchor="end">channel 63</text>')
    o.append(f'<text class="t" x="{x0 + 2 * cw}" y="130" font-size="12" font-weight="600">normal channels</text>'
             f'<text class="m" x="{x0 + 2 * cw}" y="145" font-size="11">baseline about 50</text>')
    o.append(f'<text class="t" x="{x0 + 54 * cw}" y="70" font-size="12" font-weight="600" text-anchor="middle">block of channels with a high baseline</text>'
             f'<text class="m" x="{x0 + 54 * cw}" y="85" font-size="11" text-anchor="middle">100 to 400 and more</text>')
    o.append(f'<text class="t" x="{x0 + 31 * cw}" y="138" font-size="11" text-anchor="middle">the bursts: at the edge of the block</text>')
    for k, t_ in enumerate(('VMM 4, 7 and 9 each have such a block in the upper part of the chip: about a third of the channels of',
                            'VMM 4, more than half of VMM 7 and two thirds of VMM 9. The channels that burst are the ones next to it.',
                            'Where we could compare, the block is there with the voltage off too, so it is not the chamber.',
                            'Why the block exists, and why the edge bursts, is not known: it looks like a property of the chip,',
                            'the hybrid or the connection of those channels.')):
        o.append(f'<text class="m" x="10" y="{250 + 15 * k}" font-size="11">{t_}</text>')
    o.append('</svg>')
    return '\n'.join(o)



# ------------------------------------------------------------------ 24 Oct: global thresholds per chip, first look at muon events
THR = list(csv.DictReader(open(LOCAL + 'analysis_thresholds/chip_thresholds_gain4.5.csv')))
EVM = json.load(open(LOCAL + 'analysis_muons/evmult_run76.json'))
COI = json.load(open(LOCAL + 'analysis_muons/coinc_run76.json'))


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def fig_thresholds():
    rows = [r for r in THR]
    W, l, r_, t, rh = 640, 90, 150, 62, 22
    H = t + rh * len(rows) + 50
    lo, hi = 196, 256
    pw = W - l - r_
    def X(a): return l + pw * (a - lo) / (hi - lo)
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Noise onset and safe threshold per chip, in sdt">']
    o.append('<text class="t" x="12" y="18" font-size="13" font-weight="600">Global threshold of each chip (sdt): where its noise starts, and what we run at</text>')
    o.append('<text class="m" x="12" y="34" font-size="11">sdt is the number written to the chip; lower means a lower threshold. Gain 4.5 mV/fC, peaking time 200 ns.</text>')
    o.append(f'<circle cx="16" cy="50" r="4.4" fill="var(--s1)"/><text class="m" x="26" y="54" font-size="11">noise passes 1 kHz, HV on</text>'
             f'<circle cx="190" cy="50" r="6.2" fill="none" stroke="var(--s2)" stroke-width="2.4"/><text class="m" x="202" y="54" font-size="11">same, HV off</text>'
             f'<rect x="300" y="45" width="9" height="9" fill="var(--s3)"/><text class="m" x="314" y="54" font-size="11">setting used in the muon runs (onset + 8)</text>')
    for a in range(200, 256, 10):
        o.append(f'<line x1="{X(a)}" x2="{X(a)}" y1="{t - 4}" y2="{t + rh * len(rows)}" stroke="var(--border)"/>'
                 f'<text class="f" x="{X(a)}" y="{t + rh * len(rows) + 14}" font-size="10" text-anchor="middle">{a}</text>')
    for i, r in enumerate(rows):
        v = int(r['vmm'])
        y = t + i * rh + rh / 2
        o.append(f'<text class="t" x="{l - 10}" y="{y + 4}" font-size="12" text-anchor="end">VMM {v}</text>')
        on, off, safe = _num(r['onset_sdt_HVon']), _num(r['onset_sdt_HVoff']), _num(r['safe_sdt_used'])
        if v == 6:
            o.append(f'<g><title>VMM 6: held at sdt 230 in every scan, no onset measured</title><rect x="{X(230) - 4.5:.1f}" y="{y - 4.5}" width="9" height="9" fill="var(--s3)"/></g>'
                     f'<text class="f" x="{X(230) + 12:.1f}" y="{y + 4}" font-size="10">held at 230: not measured</text>')
            continue
        pts = [p for p in (on, off, safe) if p]
        o.append(f'<line x1="{X(min(pts)):.1f}" x2="{X(max(pts)):.1f}" y1="{y}" y2="{y}" stroke="var(--fg-faint)" stroke-width="2"/>')
        if off:
            o.append(f'<g><title>VMM {v}: HV-off noise passes 1 kHz at sdt {off:g}</title><circle cx="{X(off):.1f}" cy="{y}" r="6.2" fill="none" stroke="var(--s2)" stroke-width="2.4"/></g>')
        if on:
            o.append(f'<g><title>VMM {v}: HV-on noise passes 1 kHz at sdt {on:g}</title><circle cx="{X(on):.1f}" cy="{y}" r="4.2" fill="var(--s1)"/></g>')
        if safe:
            o.append(f'<g><title>VMM {v}: setting used in the muon runs, sdt {safe:g}</title><rect x="{X(safe) - 4.5:.1f}" y="{y - 4.5}" width="9" height="9" fill="var(--s3)"/></g>')
        tag = ''
        if v in (4, 9):
            tag = 'HV-off onset only'
        elif v in (2, 7):
            tag = 'no onset found: a guess'
        elif v == 5:
            tag = 'on and off differ by 11'
        if tag:
            o.append(f'<text class="f" x="{W - r_ + 8}" y="{y + 4}" font-size="10">{tag}</text>')
    o.append('</svg>')
    return '\n'.join(o)


def thr_table():
    head = '<tr><th>VMM</th><th>assumed baseline</th><th>onset HV on (d / sdt)</th><th>onset HV off (d / sdt)</th><th>used in the muon runs (d / sdt)</th></tr>'
    rows = []
    for r in THR:
        def cell(d, s):
            return '&ndash;' if not d else f'{d} / {s}'
        rows.append(f'<tr><th>VMM {r["vmm"]}</th><td>{r["assumed_baseline_sdt"]}</td><td>{cell(r["onset_d_HVon"], r["onset_sdt_HVon"])}</td>'
                    f'<td>{cell(r["onset_d_HVoff"], r["onset_sdt_HVoff"])}</td><td>{cell(r["safe_d_used"], r["safe_sdt_used"])}</td></tr>')
    return ('<details><summary>The numbers for every chip (also as a CSV in the data folder)</summary><div style="overflow-x:auto">'
            f'<table class="bltab">{head}{"".join(rows)}</table></div></details>')


def fig_event_sizes():
    keys = list(EVM['events_by_size'].keys())
    labels = {'1-1': '1', '2-2': '2', '3-4': '3-4', '5-9': '5-9', '10-19': '10-19', '20-49': '20-49', '50-99': '50-99', '100-999': '100-999', '1000-999999999': '1000+'}
    vals = [EVM['events_by_size'][k] for k in keys]
    W, H, l, r_, t, b = 640, 330, 70, 20, 60, 58
    pw, ph = W - l - r_, H - t - b
    lo, hi = 0, 6.3
    def Y(v): return t + ph * (1 - (math.log10(max(v, 1)) - lo) / (hi - lo))
    bw = pw / len(keys)
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="How many hits belong together in time">']
    o.append('<text class="t" x="12" y="18" font-size="13" font-weight="600">Hits that belong together in time (within 100 ns), reference run, 30 minutes</text>')
    o.append('<text class="m" x="12" y="34" font-size="11">number of such groups (&#8220;events&#8221;) of each size, log scale. A muon through three chambers gives a small group.</text>')
    o.append('<text class="m" x="12" y="48" font-size="11">Orange: groups of more than 36 hits, which carry a quarter of all hits.</text>')
    for e in range(0, 7):
        y = Y(10 ** e)
        o.append(f'<line x1="{l}" x2="{W - r_}" y1="{y:.1f}" y2="{y:.1f}" stroke="var(--border)"/><text class="f" x="{l - 6}" y="{y + 3:.1f}" font-size="10" text-anchor="end">{10 ** e:,}</text>')
    for i, (k, v) in enumerate(zip(keys, vals)):
        x = l + i * bw + 6
        big = k in ('50-99', '100-999', '1000-999999999')
        o.append(f'<g><title>{v:,} events of {labels[k]} hits</title><rect x="{x:.1f}" y="{Y(v):.1f}" width="{bw - 12:.1f}" height="{t + ph - Y(v):.1f}" rx="3" fill="{"var(--s2)" if big else "var(--s1)"}"/></g>'
                 f'<text class="t" x="{x + (bw - 12) / 2:.1f}" y="{Y(v) - 5:.1f}" font-size="10" text-anchor="middle">{v:,}</text>'
                 f'<text class="m" x="{x + (bw - 12) / 2:.1f}" y="{t + ph + 15}" font-size="10" text-anchor="middle">{labels[k]}</text>')
    o.append(f'<text class="m" x="{l + pw / 2}" y="{H - 8}" font-size="11" text-anchor="middle">hits in the group</text>')
    o.append('</svg>')
    return '\n'.join(o)


def cartoon_three_kinds():
    W, H = 640, 300
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Three kinds of event in the three chambers">']
    o.append('<text class="t" x="10" y="18" font-size="13" font-weight="600">Three kinds of event in the stack (cartoon)</text>')
    kinds = [(14, 'a muon', 'one or a few hits in each chamber,', 'all at the same moment, in line', 'good: this is what we want', 'var(--s1)'),
             (222, 'noise', 'a single hit, nothing in the other', 'chambers at that moment', 'about half of all hits', 'var(--fg-faint)'),
             (430, 'a burst', 'hundreds of hits over many chips', 'and all chambers at once', 'a quarter of all hits: must be cut', 'var(--s2)')]
    for i, (x0, ttl, l1, l2, l3, col) in enumerate(kinds):
        o.append(f'<rect x="{x0}" y="30" width="196" height="224" rx="8" fill="var(--bg-elev)" stroke="var(--border)"/>')
        o.append(f'<text class="t" x="{x0 + 98}" y="50" font-size="12" font-weight="600" text-anchor="middle">{ttl}</text>')
        for k, y in enumerate((88, 138, 188)):
            o.append(f'<rect x="{x0 + 18}" y="{y}" width="160" height="14" rx="4" fill="var(--s1)" opacity="0.18" stroke="var(--s1)" stroke-width="1"/>')
        if i == 0:
            o.append(f'<line x1="{x0 + 70}" y1="66" x2="{x0 + 126}" y2="212" stroke="var(--alert)" stroke-width="2"/>')
            for (px, y) in ((86, 95), (101, 145), (116, 195)):
                o.append(f'<circle cx="{x0 + px}" cy="{y}" r="5" fill="{col}"/>')
        elif i == 1:
            o.append(f'<circle cx="{x0 + 70}" cy="145" r="5" fill="{col}"/>')
        else:
            random.seed(9)
            for y in (95, 145, 195):
                for _ in range(26):
                    o.append(f'<circle cx="{x0 + random.uniform(24, 172):.1f}" cy="{y + random.uniform(-3, 3):.1f}" r="2.4" fill="{col}" opacity="0.85"/>')
        o.append(f'<text class="m" x="{x0 + 98}" y="228" font-size="10" text-anchor="middle">{l1}</text><text class="m" x="{x0 + 98}" y="240" font-size="10" text-anchor="middle">{l2}</text>')
        o.append(f'<text class="t" x="{x0 + 98}" y="276" font-size="11" text-anchor="middle" font-weight="600">{l3}</text>')
    o.append('<text class="f" x="10" y="296" font-size="10">The three bars are the three chambers (det3, det2, det1); dots are hits.</text>')
    o.append('</svg>')
    return '\n'.join(o)



# ------------------------------------------------------------------ weekend runs (reference run, threshold scan)
_PR = {}
for _r in range(76, 96):
    _d = json.load(open(LOCAL + f'analysis_muons/per_run/ev_run_{_r}.json'))
    _PR[_r] = {'mu': _d['muon_like'] / _d['span_s'] * 60, 'hps': _d['hits'] / _d['span_s'], 'burst': _d['hits_in_large_events_pct'],
               'big': _d['large_events'] / _d['span_s']}
REF = [_PR[r] for r in range(76, 92)]
REF_MU = _st.mean(x['mu'] for x in REF)
REF_SD = _st.pstdev([x['mu'] for x in REF])
SCAN = [('safe + 16', _PR[92]), ('safe + 8', _PR[93]), ('safe', {'mu': REF_MU, 'burst': _st.mean(x['burst'] for x in REF), 'hps': _st.mean(x['hps'] for x in REF)}),
        ('safe &minus; 4', _PR[94]), ('safe &minus; 8', _PR[95])]


def fig_weekend():
    W, H = 640, 440
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Stability of the reference run and the muon-like yield against threshold">']
    o.append('<text class="t" x="12" y="18" font-size="13" font-weight="600">The weekend runs: stability over 8 hours, and the yield against threshold</text>')
    o.append('<text class="m" x="12" y="34" font-size="11">&#8220;muon-like&#8221; = a group of hits in all three chambers within 100 ns, at most 12 hits per chamber (a proxy, not an efficiency)</text>')
    # left-top: muon-like per minute per reference run
    for k, (ttl, key, lo, hi, tick, fmt) in enumerate((('muon-like events per minute, each 30-minute run', 'mu', 0, 30, (0, 10, 20, 30), '{:g}'),
                                                       ('all hits per second, each 30-minute run', 'hps', 0, 1600, (0, 400, 800, 1200, 1600), '{:g}'))):
        l, r_, t, bt = 56, 330, 78 + k * 175, 78 + k * 175 + 120
        o.append(f'<text class="t" x="{l}" y="{t - 10}" font-size="11" font-weight="600">{ttl}</text>')
        def X(i): return l + (r_ - l) * (i + 0.5) / 16
        def Y(v): return bt - (bt - t) * (v - lo) / (hi - lo)
        for a in tick:
            o.append(f'<line x1="{l}" x2="{r_}" y1="{Y(a):.1f}" y2="{Y(a):.1f}" stroke="var(--border)"/><text class="f" x="{l - 5}" y="{Y(a) + 3:.1f}" font-size="9" text-anchor="end">{fmt.format(a)}</text>')
        for i, x in enumerate(REF):
            o.append(f'<g><title>run_{76 + i}: {x[key]:.1f}</title><circle cx="{X(i):.1f}" cy="{Y(x[key]):.1f}" r="3.8" fill="var(--s1)"/></g>')
        for i in (0, 5, 10, 15):
            o.append(f'<text class="f" x="{X(i):.1f}" y="{bt + 13}" font-size="9" text-anchor="middle">{76 + i}</text>')
        o.append(f'<text class="f" x="{(l + r_) / 2}" y="{bt + 26}" font-size="9" text-anchor="middle">run number (8 hours in total)</text>')
    # right: yield against setting
    l, r_, t, bt = 380, 628, 78, 78 + 290
    lo, hi = math.log10(5), math.log10(1000)
    def Y(v): return bt - (bt - t) * (math.log10(max(v, 5)) - lo) / (hi - lo)
    o.append(f'<text class="t" x="{l - 30}" y="{t - 10}" font-size="11" font-weight="600">muon-like per minute against threshold (log)</text>')
    for a in (5, 10, 20, 50, 100, 200, 500, 1000):
        o.append(f'<line x1="{l}" x2="{r_}" y1="{Y(a):.1f}" y2="{Y(a):.1f}" stroke="var(--border)"/><text class="f" x="{l - 5}" y="{Y(a) + 3:.1f}" font-size="9" text-anchor="end">{a}</text>')
    bw = (r_ - l) / len(SCAN)
    for i, (nm, rec) in enumerate(SCAN):
        x = l + i * bw + 8
        big = rec['mu'] > 100
        o.append(f'<g><title>{nm.replace("&minus;", "-")}: {rec["mu"]:.1f} per minute, {rec["burst"]:.0f} % of the hits in large bursts</title>'
                 f'<rect x="{x:.1f}" y="{Y(rec["mu"]):.1f}" width="{bw - 16:.1f}" height="{bt - Y(rec["mu"]):.1f}" rx="3" fill="{"var(--s2)" if big else "var(--s1)"}"/></g>'
                 f'<text class="t" x="{x + (bw - 16) / 2:.1f}" y="{Y(rec["mu"]) - 5:.1f}" font-size="10" text-anchor="middle">{rec["mu"]:.0f}</text>'
                 f'<text class="m" x="{x + (bw - 16) / 2:.1f}" y="{bt + 13}" font-size="9" text-anchor="middle">{nm}</text>')
    o.append(f'<text class="f" x="{(l + r_) / 2}" y="{bt + 28}" font-size="9" text-anchor="middle">setting, relative to the &#8220;safe&#8221; threshold (lower = lower threshold)</text>')
    o.append('</svg>')
    return '\n'.join(o)

