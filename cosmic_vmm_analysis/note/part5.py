
import re
# ---------------------------------------------------------------- 25 Oct: all runs at a glance (section 13)
# (first run, last run, HV, short name for the timeline, when, measurement, config labels, length)
RUN_GROUPS = [
    (1, 3, 'on', 'first points', 'Thu 22 Oct 18:08', 'First points: the flood of hits from VMM 4 and 9 (held at sdt 400), then masks on single channels', '<code>d50</code>, <code>d40</code>, <code>m1_d50</code>', '2 min'),
    (4, 13, 'on', 'noise scan', 'Thu 22 Oct 18:19', '<strong>Noise scan, voltage on</strong>: all chips together, d = 50, 40, 34, 30, 27, 24, 22, 20, 18, 16', '<code>m2_d50</code> &hellip; <code>m2_d16</code>', '2 min'),
    (14, 23, 'off', 'noise scan', 'Fri 23 Oct 09:54', '<strong>Noise scan, voltage off</strong>: the same ten points, to separate the electronics from the chambers', '<code>m2hvoff_d50</code> &hellip; <code>m2hvoff_d16</code>', '2 min'),
    (24, 27, 'off', 'no masks', 'Fri 23 Oct 10:28', 'Noise scan, voltage off, <strong>nothing masked</strong> (to see VMM 4, 7, 9), d = 50, 40, 34, 30', '<code>m0hvoff_d50</code> &hellip; <code>m0hvoff_d30</code>', '2 min'),
    (28, 32, 'off', 'baselines', 'Fri 23 Oct 13:31', '<strong>Neighbour read-out baselines, voltage off</strong>: a probe point (28), then 0, &minus;3, &minus;6, &minus;9 DAC from each chip&rsquo;s onset (29&ndash;32)', '<code>sng1hvoff_probe</code>, <code>_o0</code>, <code>_m3</code>, <code>_m6</code>, <code>_m9</code>', '2 min'),
    (33, 34, 'none', '', 'Fri 23 Oct 16:21', 'Numbers lost in the run-counter incident; two stray captures are kept apart in <code>orphan_2026-10-23_*</code>', '&ndash;', '&ndash;'),
    (35, 47, 'off', 'block A', 'Fri 23 Oct 16:28', '<strong>Block A</strong>: VMM 4, 7, 9 unmasked, normal read-out, d = 200, 100, 60, 50, 40, 34, 30, 27, 24, 22, 20, 18, 16', '<code>n479hvoff_d200</code> &hellip; <code>_d16</code>', '2 min'),
    (48, 56, 'off', 'block C', 'Fri 23 Oct 17:12', '<strong>Block C</strong>: the same with neighbour read-out, d = 100, 60, 50, 40, 34, 27, 22, 18, 14', '<code>s479hvoff_d100</code> &hellip; <code>_d14</code>', '2 min'),
    (57, 57, 'on', '', 'Fri 23 Oct 18:31', '<strong>Health check with the voltage on</strong>, at the safe setting', '<code>safehvon_v1</code>', '10 min'),
    (58, 70, 'on', 'VMM 4/7/9', 'Fri 23 Oct 18:44', '<strong>VMM 4, 7, 9 noise scan, voltage on</strong> (the same 13 points as block A)', '<code>n479hvon_d200</code> &hellip; <code>_d16</code>', '2 min'),
    (71, 75, 'on', 'baselines', 'Fri 23 Oct 19:27', '<strong>Neighbour read-out with cosmics, voltage on</strong>: safe, &minus;4, &minus;8, &minus;12, &minus;16 (baselines for 97&nbsp;% of the channels)', '<code>sngcosmic_o0</code>, <code>_o4</code>, <code>_o8</code>, <code>_o12</code>, <code>_o16</code>', '10 min'),
    (76, 91, 'on', 'reference run', 'Sat 24 Oct 18:32 &ndash; Sun 02:55', '<strong>Reference run</strong> at the safe setting', '<code>safehvon_v1</code>', '30 min (16 runs)'),
    (92, 95, 'on', 'scan', 'Sun 25 Oct 02:55 &ndash; 06:00', '<strong>Threshold scan with tracks</strong>: safe + 16 (92), safe + 8 (93), safe &minus; 4 (94), safe &minus; 8 (95)', '<code>thrhvon_p16</code>, <code>_p8</code>, <code>_m4</code>, <code>_m8</code>', '60 min'),
    (96, 111, 'on', 'interleaved long run', 'Sun 25 Oct 13:22 &ndash; 21:34', '<strong>Interleaved long run</strong>, stopped on request (run_111 after 21 minutes); settings below', '<code>thrhvon_p8</code>, <code>safehvon_v1</code>, <code>thrhvon_m4</code>, <code>thrhvon_m6</code> in rotation', '30 min'),
    (112, 118, 'on', 'new HV', 'Mon 26 Oct 12:26', '<strong>New voltages and read-out connectors</strong> (section 15): VMM 4, 7, 9 noise ladder', '<code>n479hvon_d50</code> &hellip; <code>_d16</code>', '2 min'),
    (119, 125, 'on', 'muons', 'Mon 26 &ndash; Thu 29 Oct', '<strong>Safe setting</strong>: 119 and 123 nothing masked; 120&ndash;122 VMM 4, 7, 9 parked (121 analysed); 124 and 125 with 12 channels of VMM 4 and 9 masked (124 analysed)', '<code>safe_nomask_oldHV</code>, <code>safe_hold479_*</code>, <code>safe_m9c16to21_m4c36to41</code>', '2&ndash;30 min'),
    (126, 130, 'on', 'scint.', 'Thu 29 &ndash; Fri 30 Oct', '<strong>Scintillators on hybrid H0</strong>: 128 first wiring (no coincidence, TTL), 130 the three-paddle coincidence (60 min)', '<code>safe_m9c16to21_m4c36to41_h0open</code>', '2&ndash;60 min'),
]


def fig_run_timeline():
    W, H = 640, 250
    o = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="All runs from 1 to 130 grouped by measurement, coloured by the chamber voltage">']
    o.append('<text class="t" x="12" y="18" font-size="13" font-weight="600">All 130 runs at a glance</text>')
    o.append('<text class="m" x="12" y="34" font-size="11">Each block is one measurement, placed by run number; the colour is the chamber voltage. Hover a block for details.</text>')
    l, r_ = 20, 620

    def X(n):
        return l + (r_ - l) * (n - 0.5) / 130
    rows = {0: 70, 1: 120}
    for k, (a, b, hv, short, when, what, cfg, ln) in enumerate(RUN_GROUPS):
        col = {'on': 'var(--s1)', 'off': 'var(--c2)', 'none': 'var(--border)'}[hv]
        x0, x1 = X(a - 0.45), X(b + 0.45)
        y = 82
        txt = re.sub('<[^>]+>', '', what).replace('&minus;', '-').replace('&rsquo;', "'").replace('&nbsp;', ' ').replace('&ndash;', '-').replace('&hellip;', '...')
        o.append(f'<g><title>runs {a}-{b}: {txt}</title><rect x="{x0:.1f}" y="{y}" width="{max(x1 - x0, 2):.1f}" height="34" rx="3" fill="{col}"/></g>')
        if short:
            ly = 72 if k % 2 == 0 else 136
            anchor = 'middle'
            cx = (x0 + x1) / 2
            o.append(f'<text class="m" x="{cx:.1f}" y="{ly}" font-size="9.5" text-anchor="{anchor}">{short}</text>')
    for n in (1, 14, 28, 35, 57, 76, 96, 112, 130):
        o.append(f'<text class="f" x="{X(n):.1f}" y="152" font-size="9" text-anchor="middle">{n}</text>' if n not in (14, 28) else '')
    o.append(f'<text class="f" x="{(l + r_) / 2}" y="166" font-size="9" text-anchor="middle">run number</text>')
    # day bands
    for a, b, lab in ((1, 13, 'Thu 22'), (14, 75, 'Fri 23'), (76, 91, 'Sat 24 &ndash; Sun'), (92, 111, 'Sun 25'), (112, 130, 'Mon 26 &ndash; Fri 30')):
        o.append(f'<line x1="{X(a - 0.45):.1f}" x2="{X(b + 0.45):.1f}" y1="184" y2="184" stroke="var(--fg-faint)" stroke-width="1.5"/>'
                 f'<text class="m" x="{(X(a) + X(b)) / 2:.1f}" y="198" font-size="10" text-anchor="middle">{lab}</text>')
    o.append(f'<text class="f" x="{X(76):.1f}" y="212" font-size="9" text-anchor="middle">(nothing recorded Fri 20:24 &ndash; Sat 18:32)</text>')
    for i, (lab, col) in enumerate((('voltage on', 'var(--s1)'), ('voltage off', 'var(--c2)'), ('number not used', 'var(--border)'))):
        o.append(f'<rect x="{20 + i * 150}" y="228" width="14" height="14" rx="3" fill="{col}"/><text class="m" x="{40 + i * 150}" y="239" font-size="11">{lab}</text>')
    o.append('</svg>')
    return '\n'.join(o)


def run_groups_table():
    rows = ''.join(f'<tr><td>{a}{"&ndash;" + str(b) if b != a else ""}</td><td>{when}</td><td>{ {"on": "on", "off": "off", "none": "&ndash;"}[hv] }</td>'
                   f'<td>{what}</td><td>{cfg}</td><td>{ln}</td></tr>'
                   for a, b, hv, short, when, what, cfg, ln in RUN_GROUPS)
    return ('<div style="overflow-x:auto"><table class="runs"><tr><th>runs</th><th>when</th><th>voltage</th><th>measurement</th>'
            f'<th>config label</th><th>length each</th></tr>{rows}</table></div>')
