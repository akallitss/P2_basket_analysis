
# ================================================================== NEW PART (v2 of the note)
import csv
import json
import random
import re as _re

LOCAL = '/local/home/ak271430/Documents/PostDocSaclay/data/Lab_Measurements/VMM-cosmics/vmm_threshold_scan/'


def _jload(sub):
    return json.load(open(LOCAL + sub + '/noise_scan_summary.json'))


def _cload(sub):
    return list(csv.DictReader(open(LOCAL + sub + '/noise_scan_summary.csv')))


J_ON, J_OFF = _jload('analysis_noise_scan'), _jload('analysis_noise_scan_hvoff')
C_ON, C_OFF = _cload('analysis_noise_scan'), _cload('analysis_noise_scan_hvoff')


def _by_d(J, tag):
    out, runs = {}, {}
    for run, v in J.items():
        m = _re.search(r'_%s_d(\d+)$' % tag, v['chip_config'])
        if m:
            out[int(m.group(1))] = v['khz_per_vmm']
            runs[int(m.group(1))] = int(run.split('_')[1])
    return out, runs


ON_D, ON_RUN = _by_d(J_ON, 'm2')
OFF_D, OFF_RUN = _by_d(J_OFF, 'm2hvoff')
OFF0_D, OFF0_RUN = _by_d(J_OFF, 'm0hvoff')
assert sorted(ON_D) == sorted(D) and sorted(OFF_D) == sorted(D), (sorted(ON_D), sorted(OFF_D))
ROWS_ON = {d: [float(ON_D[d].get(str(v), 0.0)) for v in CHIPS] for d in D}
ROWS_OFF = {d: [float(OFF_D[d].get(str(v), 0.0)) for v in CHIPS] for d in D}
# cross-check against the numbers typed in by hand for the first version of the note
for d in D:
    for k, v in enumerate(CHIPS):
        if abs(ROWS_ON[d][k] - ROWS[d][k]) > 0.011 + 0.01 * ROWS[d][k]:
            print('WARNING hand-typed vs file differ', d, v, ROWS[d][k], ROWS_ON[d][k])
MB_OFF = {50: 5, 40: 11, 34: 206, 30: 136, 27: 153, 24: 291, 22: 379, 20: 781, 18: 1778, 16: 3537}  # from config_scan.log


def onset_of(rows, k):
    for d in D:
        if rows[d][k] >= 1:
            return d
    return None


def fig_onset_pair():
    chips = [v for v in CHIPS if v != 6]
    on = {v: onset_of(ROWS_ON, CHIPS.index(v)) for v in chips}
    off = {v: onset_of(ROWS_OFF, CHIPS.index(v)) for v in chips}
    # VMM 7: both are a single-channel burst; keep it but mark it
    order = sorted(chips, key=lambda v: (-(on[v] or 0), -(off[v] or 0), v))
    W, l, r, t, rh = 620, 100, 120, 56, 22
    H = t + rh * len(order) + 64
    pw = W - l - r
    def X(d): return l + pw * (50 - d) / (50 - 14)
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="First step above 1 kHz, HV on and HV off">']
    o.append('<text class="t" x="20" y="18" font-size="13" font-weight="600">Where does each chip start to fire (&gt; 1 kHz)? Voltage on and voltage off</text>')
    o.append('<text class="m" x="20" y="34" font-size="11">same chip settings in both scans; further left = the chip starts firing earlier</text>')
    for d in (50, 40, 30, 20):
        o.append(f'<line x1="{X(d)}" x2="{X(d)}" y1="{t - 8}" y2="{t + rh * len(order)}" stroke="var(--border)"/>'
                 f'<text class="m" x="{X(d)}" y="{t + rh * len(order) + 16}" font-size="11" text-anchor="middle">d={d}</text>')
    for i, v in enumerate(order):
        y = t + i * rh + 11
        o.append(f'<text class="t" x="{l - 12}" y="{y + 4}" font-size="12" text-anchor="end">VMM {v}</text>')
        a, b = on[v], off[v]
        xa = X(a) if a else X(14)
        xb = X(b) if b else X(14)
        o.append(f'<line x1="{xa}" x2="{xb}" y1="{y}" y2="{y}" stroke="var(--fg-faint)" stroke-width="2"/>')
        if a:
            o.append(f'<g><title>VMM {v}, HV on: first above 1 kHz at d={a}</title><circle cx="{xa}" cy="{y}" r="6" fill="var(--s1)"/></g>')
        else:
            o.append(f'<text class="m" x="{X(14) + 12}" y="{y + 4}" font-size="10">on: never</text>')
        if b:
            o.append(f'<g><title>VMM {v}, HV off: first above 1 kHz at d={b}</title><circle cx="{xb}" cy="{y}" r="9.5" fill="none" stroke="var(--s2)" stroke-width="3"/></g>')
        else:
            o.append(f'<text class="m" x="{X(14) + (62 if not a else 12)}" y="{y + 4}" font-size="10">off: never</text>')
    ly = t + rh * len(order) + 38
    o.append(f'<circle cx="{l}" cy="{ly}" r="6" fill="var(--s1)"/><text class="m" x="{l + 12}" y="{ly + 4}" font-size="11">HV on</text>'
             f'<circle cx="{l + 100}" cy="{ly}" r="9.5" fill="none" stroke="var(--s2)" stroke-width="3"/><text class="m" x="{l + 116}" y="{ly + 4}" font-size="11">HV off (ring)</text>')
    o.append('</svg>')
    return '\n'.join(o)


def fig_volume_pair():
    W, H, l, r, t, b = 640, 320, 62, 20, 50, 52
    pw, ph = W - l - r, H - t - b
    lo, hi = 0, math.log10(10000)
    def Y(v): return t + ph * (1 - (math.log10(max(v, 1)) - lo) / (hi - lo))
    gw = pw / len(D)
    def lab(v): return f'{v / 1000:.1f}k' if v >= 1000 else str(v)
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Data per minute, HV on and off">']
    o.append('<text class="t" x="20" y="18" font-size="13" font-weight="600">Data written per minute, all chips together</text>')
    o.append(f'<rect x="{l}" y="28" width="14" height="10" rx="2" fill="var(--s1)"/><text class="m" x="{l + 20}" y="37" font-size="11">HV on</text>'
             f'<rect x="{l + 80}" y="28" width="14" height="10" rx="2" fill="var(--s2)"/><text class="m" x="{l + 100}" y="37" font-size="11">HV off</text>')
    for e in range(0, 5):
        y = Y(10 ** e)
        o.append(f'<line x1="{l}" x2="{W - r}" y1="{y}" y2="{y}" stroke="var(--border)"/>'
                 f'<text class="m" x="{l - 8}" y="{y + 4}" font-size="11" text-anchor="end">{10 ** e:,} MB</text>')
    for j, d in enumerate(D):
        x0 = l + j * gw + 4
        bw = (gw - 8) / 2
        for k, (v, col, nm) in enumerate(((MB_MIN[d], 'var(--s1)', 'HV on'), (MB_OFF[d], 'var(--s2)', 'HV off'))):
            x = x0 + k * bw
            o.append(f'<g><title>d={d}, {nm}: {v} MB/min</title><rect x="{x}" y="{Y(v)}" width="{bw - 1}" height="{t + ph - Y(v)}" rx="2" fill="{col}"/></g>'
                     f'<text class="m" x="{x + bw / 2}" y="{Y(v) - 4}" font-size="9" text-anchor="middle">{lab(v)}</text>')
        o.append(f'<text class="m" x="{x0 + bw}" y="{t + ph + 16}" font-size="11" text-anchor="middle">{d}</text>')
    o.append(f'<text class="m" x="{l + pw / 2}" y="{H - 10}" font-size="12" text-anchor="middle">threshold step d: quiet → noisy</text>')
    o.append('</svg>')
    return '\n'.join(o)


def fig_bursts_off():
    steps = [50, 40, 34, 30]
    data = {v: [float(J_OFF['run_%d' % OFF0_RUN[d]]['khz_per_vmm'].get(str(v), 0.0)) for d in steps] for v in (4, 7, 9)}
    W, H = 640, 270
    t, b = 62, 52
    ph = H - t - b
    lo, hi = math.log10(0.01), math.log10(3000)
    def Y(v): return t + ph * (1 - (math.log10(max(v, 0.01)) - lo) / (hi - lo))
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Unmasked chips with the HV off">']
    o.append('<text class="t" x="10" y="18" font-size="13" font-weight="600">Nothing masked, HV off: the three burst chips still burst</text>')
    o.append('<text class="m" x="10" y="34" font-size="11">kHz per chip on a log scale; x: threshold step d</text>')
    pw = (W - 20) / 3
    for i, v in enumerate((4, 7, 9)):
        x0 = 10 + i * pw
        o.append(f'<text class="t" x="{x0 + pw / 2}" y="54" font-size="12" font-weight="600" text-anchor="middle">VMM {v}</text>')
        for e in (-2, 0, 2):
            y = Y(10 ** e)
            o.append(f'<line x1="{x0 + 34}" x2="{x0 + pw - 6}" y1="{y}" y2="{y}" stroke="var(--border)"/>'
                     f'<text class="m" x="{x0 + 30}" y="{y + 4}" font-size="10" text-anchor="end">{10 ** e:g}</text>')
        bw = (pw - 60) / 4 - 6
        for j, d in enumerate(steps):
            val = data[v][j]
            x = x0 + 40 + j * (bw + 6)
            o.append(f'<g><title>VMM {v}, d={d}: {val:g} kHz</title><rect x="{x}" y="{Y(val)}" width="{bw}" height="{max(0, t + ph - Y(val))}" rx="3" fill="var(--s2)"/></g>'
                     f'<text class="t" x="{x + bw / 2}" y="{Y(val) - 4}" font-size="10" text-anchor="middle">{val:.0f}</text>'
                     f'<text class="m" x="{x + bw / 2}" y="{t + ph + 15}" font-size="10" text-anchor="middle">{d}</text>')
    o.append('</svg>')
    return '\n'.join(o)


def fig_class():
    rows = list(csv.DictReader(open(LOCAL + 'analysis_noise_scan_hvoff/hv_on_off_channels.csv')))
    cnt = {}
    for r in rows:
        cnt.setdefault(int(r['vmm']), {'HV-dependent': 0, 'partly': 0, 'HV-independent': 0})[r['class']] += 1
    chips = sorted(cnt)
    W, l, r_, t, rh = 620, 80, 30, 58, 22
    H = t + rh * len(chips) + 34
    mx = max(sum(c.values()) for c in cnt.values())
    pw = W - l - r_
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Which noisy channels need the voltage">']
    o.append('<text class="t" x="20" y="18" font-size="13" font-weight="600">Noisy channels: do they need the voltage?</text>')
    o.append(f'<rect x="{l}" y="28" width="14" height="10" rx="2" fill="var(--s1)"/><text class="m" x="{l + 20}" y="37" font-size="11">gone with HV off</text>'
             f'<rect x="{l + 140}" y="28" width="14" height="10" rx="2" fill="var(--fg-faint)"/><text class="m" x="{l + 160}" y="37" font-size="11">weaker</text>'
             f'<rect x="{l + 230}" y="28" width="14" height="10" rx="2" fill="var(--s2)"/><text class="m" x="{l + 250}" y="37" font-size="11">still there with HV off</text>')
    o.append('<text class="f" x="20" y="52" font-size="10">each unit = one channel at one threshold step that fired with HV on</text>')
    for i, v in enumerate(chips):
        y = t + i * rh
        o.append(f'<text class="t" x="{l - 10}" y="{y + 14}" font-size="12" text-anchor="end">VMM {v}</text>')
        x = l
        for key, col in (('HV-dependent', 'var(--s1)'), ('partly', 'var(--fg-faint)'), ('HV-independent', 'var(--s2)')):
            n = cnt[v][key]
            w = pw * n / mx
            if n:
                o.append(f'<g><title>VMM {v}: {n} channel-steps {key}</title><rect x="{x}" y="{y + 2}" width="{max(w - 1, 1)}" height="{rh - 5}" rx="2" fill="{col}"/></g>')
                if w > 26:
                    o.append(f'<text x="{x + w / 2}" y="{y + 15}" font-size="10" text-anchor="middle" fill="#fff">{n}</text>')
            x += w
    o.append('</svg>')
    return '\n'.join(o)


def fig_channel_curves():
    def series(C, tag, vmm, ch):
        out = {}
        for r in C:
            m = _re.search(r'_%s_d(\d+)$' % tag, r['chip_config'])
            if m and int(r['vmm']) == vmm and int(r['ch']) == ch:
                out[int(m.group(1))] = float(r['rate_hz'])
        return out
    S = [('VMM 7 ch 23, HV on', series(C_ON, 'm2', 7, 23), 'var(--s1)', False),
         ('VMM 7 ch 23, HV off', series(C_OFF, 'm2hvoff', 7, 23), 'var(--s1)', True),
         ('VMM 5 ch 46, HV on', series(C_ON, 'm2', 5, 46), 'var(--s2)', False)]
    W, H, l, r, t, b = 640, 346, 70, 20, 80, 50
    pw, ph = W - l - r, H - t - b
    lo, hi = 0, 6
    def X(d): return l + pw * (50 - d) / (50 - 16)
    def Y(v): return t + ph * (1 - (math.log10(max(v, 1)) - lo) / (hi - lo))
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Hit rate of single channels against threshold">']
    o.append('<text class="t" x="20" y="18" font-size="13" font-weight="600">One channel, many thresholds: the rate rises, peaks, and falls again</text>')
    o.append('<text class="m" x="20" y="34" font-size="11">real data; hits per second of one channel (log scale; points at 1 Hz mean 1 Hz or fewer)</text>')
    for k, (nm, _, col, dashed) in enumerate(S):
        x = l + k * 170
        o.append(f'<circle cx="{x + 6}" cy="52" r="5" fill="{"var(--surface)" if dashed else col}" stroke="{col}" stroke-width="2.5"/>'
                 f'<text class="m" x="{x + 16}" y="56" font-size="11">{nm}</text>')
    for e in range(0, 7):
        y = Y(10 ** e)
        o.append(f'<line x1="{l}" x2="{W - r}" y1="{y}" y2="{y}" stroke="var(--border)"/>'
                 f'<text class="m" x="{l - 8}" y="{y + 4}" font-size="11" text-anchor="end">{10 ** e:,}</text>')
    for d in D:
        o.append(f'<text class="m" x="{X(d)}" y="{t + ph + 16}" font-size="11" text-anchor="middle">{d}</text>')
    o.append(f'<text class="m" x="{l + pw / 2}" y="{H - 8}" font-size="12" text-anchor="middle">threshold step d (quiet → noisy)</text>')
    for nm, sr, col, dashed in S:
        pts = [(X(d), Y(sr[d])) for d in D if d in sr and sr[d] > 0]
        if len(pts) > 1:
            o.append(f'<polyline fill="none" stroke="{col}" stroke-width="2.5"' + (' stroke-dasharray="6 4"' if dashed else '') +
                     ' points="' + ' '.join(f'{x:.1f},{y:.1f}' for x, y in pts) + '"/>')
        for d in D:
            if d in sr and sr[d] > 0:
                o.append(f'<g><title>{nm}, d={d}: {sr[d]:.0f} Hz</title><circle cx="{X(d):.1f}" cy="{Y(sr[d]):.1f}" r="4.5" fill="{"var(--surface)" if dashed else col}" stroke="{col}" stroke-width="2.5"/></g>')
    o.append(f'<text class="t" x="{X(34)}" y="{Y(489137) - 10}" font-size="11" text-anchor="middle">peak</text>')
    o.append(f'<text class="t" x="{X(20)}" y="{Y(415938) - 10}" font-size="11" text-anchor="middle">peak</text>')
    o.append('</svg>')
    return '\n'.join(o)


def cartoon_chain():
    W, H = 640, 330
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Where noise can come from, voltage on and off">']
    def row(y, title, on):
        col = 'var(--s1)'
        o.append(f'<text class="t" x="10" y="{y}" font-size="13" font-weight="600">{title}</text>')
        o.append(f'<rect x="10" y="{y + 14}" width="170" height="92" rx="8" fill="{col if on else "var(--border)"}" opacity="{0.20 if on else 0.55}" stroke="{col if on else "var(--fg-faint)"}" stroke-width="2"/>'
                 f'<text class="t" x="95" y="{y + 38}" font-size="12" text-anchor="middle" font-weight="600">chamber (gas + pads)</text>')
        if on:
            o.append(f'<text class="m" x="95" y="{y + 58}" font-size="11" text-anchor="middle">voltage on: the gas multiplies</text>'
                     f'<text class="m" x="95" y="{y + 73}" font-size="11" text-anchor="middle">charge, so pads can give</text>'
                     f'<text class="m" x="95" y="{y + 88}" font-size="11" text-anchor="middle">signals, sparks, small noise</text>')
        else:
            o.append(f'<text class="m" x="95" y="{y + 58}" font-size="11" text-anchor="middle">voltage off: no gas gain,</text>'
                     f'<text class="m" x="95" y="{y + 73}" font-size="11" text-anchor="middle">the pads give no signals</text>'
                     f'<text class="m" x="95" y="{y + 88}" font-size="11" text-anchor="middle">at all</text>')
        o.append(f'<path d="M180 {y + 60} H222" stroke="var(--fg-faint)" stroke-width="3"/>'
                 f'<text class="f" x="201" y="{y + 52}" font-size="10" text-anchor="middle">cable</text>')
        o.append(f'<rect x="222" y="{y + 14}" width="190" height="92" rx="8" fill="var(--s2)" opacity="0.20" stroke="var(--s2)" stroke-width="2"/>'
                 f'<text class="t" x="317" y="{y + 38}" font-size="12" text-anchor="middle" font-weight="600">hybrid + VMM chip</text>'
                 f'<text class="m" x="317" y="{y + 58}" font-size="11" text-anchor="middle">64 channels, each with a</text>'
                 f'<text class="m" x="317" y="{y + 73}" font-size="11" text-anchor="middle">threshold; its own electronic</text>'
                 f'<text class="m" x="317" y="{y + 88}" font-size="11" text-anchor="middle">noise is always there</text>')
        o.append(f'<path d="M412 {y + 60} H440" stroke="var(--fg-faint)" stroke-width="3"/>')
        o.append(f'<rect x="440" y="{y + 30}" width="190" height="60" rx="8" fill="var(--bg-elev)" stroke="var(--border)"/>'
                 f'<text class="t" x="535" y="{y + 56}" font-size="12" text-anchor="middle" font-weight="600">hits we count</text>'
                 f'<text class="m" x="535" y="{y + 74}" font-size="11" text-anchor="middle">{"detector + electronics" if on else "electronics only"}</text>')
    row(20, 'Scan 1: voltage ON', True)
    row(170, 'Scan 2: voltage OFF (same chip settings)', False)
    o.append('<text class="m" x="10" y="318" font-size="11">Compare the two: what disappears needs the detector; what stays comes from the chip, the hybrid or the cable.</text>')
    o.append('</svg>')
    return '\n'.join(o)


def cartoon_scurve():
    random.seed(11)
    W, H = 640, 372
    n = 90
    base = [0.0]
    for _ in range(n - 1):
        base.append(0.55 * base[-1] + random.gauss(0, 0.9))
    sd = (sum(x * x for x in base) / n) ** 0.5
    base = [x / sd for x in base]
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Why a channel fires most when the threshold is at its baseline">']
    o.append('<text class="t" x="10" y="18" font-size="13" font-weight="600">Why a channel fires most when the threshold sits at its own baseline (cartoon)</text>')
    pw = (W - 40) / 3
    cases = [(3.0, 'A: threshold well above', 'noise almost never reaches it', 'quiet'),
             (0.3, 'B: threshold at the baseline', 'noise crosses it all the time', 'many hits'),
             (-3.0, 'C: threshold well below', 'always above it: no new crossings', 'quiet again')]
    for i, (th, ttl, sub, res) in enumerate(cases):
        x0 = 14 + i * (pw + 6)
        yc, sc = 130, 26
        o.append(f'<rect x="{x0}" y="36" width="{pw}" height="188" rx="8" fill="var(--bg-elev)" stroke="var(--border)"/>')
        o.append(f'<text class="t" x="{x0 + pw / 2}" y="56" font-size="12" font-weight="600" text-anchor="middle">{ttl}</text>')
        o.append(f'<text class="m" x="{x0 + pw / 2}" y="71" font-size="10" text-anchor="middle">{sub}</text>')
        pts = ' '.join(f'{x0 + 8 + k * (pw - 16) / (n - 1):.1f},{yc - sc * base[k] * 0.55:.1f}' for k in range(n))
        o.append(f'<polyline fill="none" stroke="var(--s1)" stroke-width="1.6" points="{pts}"/>')
        ty = yc - sc * th * 0.55
        o.append(f'<line x1="{x0 + 6}" x2="{x0 + pw - 6}" y1="{ty:.1f}" y2="{ty:.1f}" stroke="var(--alert)" stroke-width="2.5" stroke-dasharray="6 4"/>')
        cr = 0
        for k in range(1, n):
            if base[k - 1] < th <= base[k]:
                cr += 1
                xk = x0 + 8 + k * (pw - 16) / (n - 1)
                o.append(f'<circle cx="{xk:.1f}" cy="{ty:.1f}" r="3.5" fill="var(--alert)"/>')
        o.append(f'<text class="t" x="{x0 + pw / 2}" y="214" font-size="12" text-anchor="middle" font-weight="600">{res}</text>')
    o.append('<text class="m" x="14" y="242" font-size="11">blue: the channel\'s baseline wobbling with electronic noise. Red dashed: the threshold. Red dots: a hit is recorded</text>')
    o.append('<text class="m" x="14" y="256" font-size="11">each time the blue line crosses the threshold going up.</text>')
    # the bell
    bx0, bw, by0, bh = 120, 400, 292, 50
    o.append(f'<line x1="{bx0}" x2="{bx0 + bw}" y1="{by0 + bh}" y2="{by0 + bh}" stroke="var(--border)"/>')
    pts = []
    for k in range(0, 101):
        t_ = 4 - 8 * k / 100
        rate = math.exp(-0.5 * (t_ - 0.0) ** 2)
        pts.append((bx0 + bw * k / 100, by0 + bh * (1 - rate)))
    o.append('<polyline fill="none" stroke="var(--s1)" stroke-width="2.5" points="' + ' '.join(f'{x:.1f},{y:.1f}' for x, y in pts) + '"/>')
    for (tt, lab) in ((3.0, 'A'), (0.3, 'B'), (-3.0, 'C')):
        x = bx0 + bw * (4 - tt) / 8
        rate = math.exp(-0.5 * tt ** 2)
        o.append(f'<circle cx="{x:.1f}" cy="{by0 + bh * (1 - rate):.1f}" r="5" fill="var(--alert)"/><text class="t" x="{x:.1f}" y="{by0 + bh * (1 - rate) - 8:.1f}" font-size="11" text-anchor="middle">{lab}</text>')
    o.append(f'<text class="m" x="{bx0 - 10}" y="{by0 + 26}" font-size="11" text-anchor="end">hits per second</text>')
    o.append(f'<text class="m" x="{bx0 + bw / 2}" y="{by0 + bh + 14}" font-size="11" text-anchor="middle">threshold: high ← → low. The peak is the channel\'s own baseline.</text>')
    o.append('</svg>')
    return '\n'.join(o)


def cartoon_global_vs_channel():
    random.seed(3)
    W, H = 640, 236
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="One threshold per chip versus one per channel">']
    o.append('<text class="t" x="10" y="18" font-size="13" font-weight="600">One threshold for the whole chip, or a small correction for each channel? (cartoon)</text>')
    for p, (x0, title, spread, lowline) in enumerate(((14, 'Now: one threshold per chip', 11, False), (330, 'Next: one correction per channel', 2.5, True))):
        o.append(f'<text class="t" x="{x0}" y="46" font-size="12" font-weight="600">{title}</text>')
        base_y = 160
        vals = [random.gauss(0, spread) for _ in range(32)]
        if p == 0:
            vals[17] = 24  # the "VMM 7 ch 23" kind of channel
        for k, v in enumerate(vals):
            x = x0 + 4 + k * 9.3
            h = 70 + v
            hot = p == 0 and k == 17
            o.append(f'<rect x="{x:.1f}" y="{base_y - h:.1f}" width="6" height="{h:.1f}" rx="2" fill="{"var(--s2)" if hot else "var(--s1)"}"/>')
        ly = base_y - 70 - max(vals) - 5
        o.append(f'<line x1="{x0}" x2="{x0 + 300}" y1="{ly:.1f}" y2="{ly:.1f}" stroke="var(--alert)" stroke-width="2"/>')
        o.append(f'<text class="m" x="{x0 + 298}" y="{ly - 5:.1f}" font-size="11" text-anchor="end">the chip\'s threshold</text>')
        o.append(f'<text class="m" x="{x0}" y="{base_y + 16}" font-size="11">one bar = the baseline of one channel</text>')
    o.append('<text class="m" x="14" y="206" font-size="11">Left: the single tall channel (orange, like VMM 7 ch 23) forces the threshold up for everybody. Right: after the per-channel</text>')
    o.append('<text class="m" x="14" y="221" font-size="11">correction the channels line up, so the chip\'s threshold can be lower and the muons that were lost are counted again.</text>')
    o.append('</svg>')
    return '\n'.join(o)


def stepcard(n, title, did, saw, told, nxt, status='done'):
    pill = {'done': '<span class="pill ok">done</span>', 'next': '<span class="pill warn">next, not done yet</span>',
            'later': '<span class="pill">later</span>'}[status]
    lab = (('What we did', 'What we saw', 'What it tells us', 'So the next step is') if status == 'done' else
           ('Where we are', 'What is needed', 'What it will tell us', 'Then'))
    return (f'<div class="stepcard {status}"><div class="sc-head"><span class="stepn">{n}</span><strong>{title}</strong> {pill}</div>'
            f'<dl><dt>{lab[0]}</dt><dd>{did}</dd><dt>{lab[1]}</dt><dd>{saw}</dd>'
            f'<dt>{lab[2]}</dt><dd>{told}</dd><dt>{lab[3]}</dt><dd>{nxt}</dd></dl></div>')


