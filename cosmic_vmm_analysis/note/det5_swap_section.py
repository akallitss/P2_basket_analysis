#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Section "det5 at P2_MID and the read-out swap" of the note "Where do the missed muons go?" (25 Nov 2026).

What changed on the stand on 23-24 Nov (det5 replaces det2 at P2_MID, Kapton on the det5 c4 card, P2_MID <-> P2_OUT
read-out swap, new DAQ PC) and what the funnel (pointing radius + edge margin + paddle cut, REF configuration) gives
before and after. The decisive table is per connector: the same pads of a chamber read by other chips.

Inputs (muon_funnel_summary.py outputs, same cuts):
  --weekend   weekend 19-23 Nov, det2 at P2_MID (runs 268-388)
  --night     night 23-24 Nov, det5 at P2_MID, default read-out (runs 389-406)
  --swap      24-25 Nov, det5 at P2_MID, MID <-> OUT read-out swapped (runs 407, 410-424; 408/409 excluded)
  --timing    optional JSON from in_timing_check.py (P2_IN hit-time check old vs new DAQ PC)
Output: an HTML fragment (--out), inserted by muon_funnel_page.py (SWAP_HTML) or in place with --note.
"""
import argparse
import html
import json
import math
import re

K = ('found', 'probe_silent', 'probe_far', 'probe_near_late')
PROBES = ('P2_IN', 'P2_MID', 'P2_OUT')
# read-out per period: station -> {connector: (bottom VMM, top VMM)}
DEFAULT = {'P2_IN': {3: (2, 3), 4: (4, 5), 5: (6, 7)}, 'P2_MID': {3: (8, 9), 4: (10, 11), 5: (12, 13)},
           'P2_OUT': {3: (14, 15), 4: (16, 17), 5: (18, 19)}}
SWAP = {'P2_IN': DEFAULT['P2_IN'], 'P2_MID': DEFAULT['P2_OUT'], 'P2_OUT': DEFAULT['P2_MID']}
GAINS = (('REF', 'REF (gain 4.5 mV/fC, trimmed)'), ('g3.0 gscan_hvon50_hvm', 'D50 (gain 3, 50 Hz noise level)'),
         ('g6.0 gscan_hvon50_hvm', 'B50 (gain 6, 50 Hz noise level)'))
STYLE = """<style>
  :root { --b4:#9aa3ad; --af:#2a78d6; --ok:#2b9a5a; }
  @media (prefers-color-scheme: dark) { :root { --b4:#7d8691; --af:#3987e5; --ok:#3cb371; } }
  .swfig svg { width: 100%; height: auto; max-width: 760px; display: block; }
  .swfig text { font-size: 11px; fill: currentColor; }
  .swtab td, .swtab th { text-align: center; }
  .swtab td.l, .swtab th.l { text-align: left; }
</style>"""


def probe(summary, group, p):
    g = json.load(open(summary))['groups'].get(group)
    if not g:
        return None
    o = g['probes'][p]['on']
    n = sum(o.get(k, 0) for k in K)
    chips = {}
    for c, v in g['probes'][p]['chip'].items():
        m = sum(v.get(k, 0) for k in K)
        if m:
            chips[int(c)] = (v.get('probe_silent', 0) / m, m)
    runs = g.get('runs')
    return dict(n=n, fr={k: o.get(k, 0) / n for k in K}, chips=chips,
                runs=len(runs) if isinstance(runs, list) else runs)


def err(p, n):
    return math.sqrt(p * (1 - p) / n) if n else float('nan')


def pct(x, d=1):
    return f'{100 * x:.{d}f}'


def per_connector(r, layout, station):
    """[(connector, half, vmm, silent, n)] for one station."""
    out = []
    for c, (b, t) in sorted(layout[station].items()):
        for half, v in (('bot', b), ('top', t)):
            s, m = r['chips'].get(v, (float('nan'), 0))
            out.append((c, half, v, s, m))
    return out


def svg_positions(before, after, title_b, title_a):
    """grouped bars: silent share per connector half, two periods."""
    W, H, L, B, T = 760, 230, 46, 40, 24
    rows = list(zip(before, after))
    n = len(rows)
    gw = (W - L - 10) / n
    ymax = 0.30
    y = lambda v: T + (H - T - B) * (1 - min(v, ymax) / ymax)
    s = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="silent share per connector half, {html.escape(title_b)} vs '
         f'{html.escape(title_a)}">']
    for t in (0, 0.1, 0.2, 0.3):
        s.append(f'<line x1="{L}" x2="{W - 10}" y1="{y(t):.1f}" y2="{y(t):.1f}" stroke="currentColor" '
                 f'stroke-opacity=".15"/><text x="{L - 6}" y="{y(t) + 4:.1f}" text-anchor="end">{int(t * 100)} %</text>')
    for i, ((c, half, vb, sb, nb), (_, _, va, sa, na)) in enumerate(rows):
        x0 = L + i * gw + gw * 0.15
        bw = gw * 0.33
        for j, (v, sv, col, lab) in enumerate(((vb, sb, 'var(--b4)', title_b), (va, sa, 'var(--af)', title_a))):
            if not math.isfinite(sv):
                continue
            x = x0 + j * (bw + 2)
            s.append(f'<rect x="{x:.1f}" y="{y(sv):.1f}" width="{bw:.1f}" height="{y(0) - y(sv):.1f}" fill="{col}">'
                     f'<title>c{c} {half}: VMM {v}, {pct(sv)} % silent ({lab}, n = {nb if j == 0 else na})</title></rect>'
                     f'<text x="{x + bw / 2:.1f}" y="{y(sv) - 3:.1f}" text-anchor="middle">{pct(sv, 0)}</text>'
                     f'<text x="{x + bw / 2:.1f}" y="{H - B + 14:.1f}" text-anchor="middle" opacity=".7">{v}</text>')
        s.append(f'<text x="{x0 + bw + 1:.1f}" y="{H - B + 29:.1f}" text-anchor="middle">c{c} {half}</text>')
    s.append(f'<rect x="{L + 4}" y="4" width="10" height="10" fill="var(--b4)"/><text x="{L + 18}" y="13">'
             f'{html.escape(title_b)}</text><rect x="{L + 290}" y="4" width="10" height="10" fill="var(--af)"/>'
             f'<text x="{L + 304}" y="13">{html.escape(title_a)}</text></svg>')
    return ''.join(s)


def section(a):
    W = {g: {p: probe(a.weekend, g, p) for p in PROBES} for g, _ in GAINS}
    N = {g: {p: probe(a.night, g, p) for p in PROBES} for g, _ in GAINS}
    S = {g: {p: probe(a.swap, g, p) for p in PROBES} for g, _ in GAINS}
    tim = json.load(open(a.timing)) if a.timing else None
    h = [STYLE, '<h2>10. det5 at P2_MID and the read-out swap (23-25 Nov)</h2>']
    h.append('<p>Section 9 left two candidates for the P2_MID silent muons: the chamber (det2) or the place - the MID '
             'position in the stack and its read-out (hybrids H4-H6, VMM 8-13). On 23-24 Nov both were changed in turn.</p>')
    h.append('<h3>10a. What changed on the stand</h3><div style="overflow-x:auto"><table class="swtab"><tr><th class="l">'
             'when</th><th class="l">change</th><th>runs</th></tr>'
             '<tr><td class="l">23 Nov, between run_388 and run_389</td><td class="l"><b>det5 replaces det2 at P2_MID</b>. '
             'det5 is bulked at CERN (CERN insulation mask: 0.5 mm pillars on a 2 mm grid, 4.87 % dead area, vs 3.78 % '
             'for the Saclay mask V1 of det1-det4). Same connectors c3-c5, same hybrids H4-H6 (VMM 8-13), same HV '
             'channels, 430 / 730 V like det2.</td><td>389-406</td></tr>'
             '<tr><td class="l">23 Nov</td><td class="l"><b>Kapton</b> on the det5 c4 adaptation card (a few edge '
             'contacts after pin 35 taped, so that det5 holds its HV).</td><td></td></tr>'
             '<tr><td class="l">24 Nov ~12:00 (setup entry 12:14)</td><td class="l"><b>P2_MID &harr; P2_OUT read-out '
             'swap</b>: the adaptation cards, with their cables, were moved between the two chambers connector for '
             'connector. Now det5 (MID) is read by H7-H9 (VMM 14-19) and det1 (OUT) by H4-H6 (VMM 8-13). HV cabling '
             'unchanged; the Kapton was redone on the card now on det5 c4 (H8).</td><td>407, 410-424</td></tr>'
             '<tr><td class="l">24 Nov 12:54</td><td class="l">run_407 stopped by hand (24 of 49 min kept) for the '
             '<b>DAQ PC change</b> (dedippce185 &rarr; dedippcy245, same hybrids, FEC and chip configurations); run_408 '
             '(HV going off) and the HV-off test run_409 are left out.</td><td></td></tr>'
             '<tr><td class="l">24 Nov 19:13 &rarr; 25 Nov 07:38</td><td class="l">night on the new PC: 11 REF, 2 D50, '
             '1 B50, 1 pulser run, 49 min each, no HV trip.</td><td>410-424</td></tr></table></div>')
    h.append('<h3>10b. The Kapton: which det5 pads are blind</h3><p>From the muon signals of the night REF runs '
             '391-406 (in-time hits per channel, det5 vs det2 on the same channels, each chip normalised to its median '
             'channel): only <b>two channels</b> lost their muons, VMM 11 ch 27 (no hit at all) and VMM 11 ch 23 (~85 % '
             'lost); the rest of the taped range is live. That matches the photo: the tape covers about 2&frac12; '
             'contacts. Three other det5 channels are silent all night and are not under the tape (VMM 9 ch 0, VMM 9 '
             'ch 35, VMM 11 ch 52): dead pads or contacts of det5 itself. The four channels without any hit are in the '
             'per-run dead list, so muons predicted on them are not tested; VMM 11 ch 23 still fires now and then, is '
             'tested, and adds a little to the det5 silent share (one pad of ~380 read).</p>')
    # chamber table
    rows = []
    for lab, D, chamber in (('19-23 Nov (det2 at MID)', W, {'P2_IN': 'det3', 'P2_MID': 'det2', 'P2_OUT': 'det1'}),
                            ('23-24 Nov (det5 at MID)', N, {'P2_IN': 'det3', 'P2_MID': 'det5', 'P2_OUT': 'det1'}),
                            ('24-25 Nov (det5 at MID, read-out swapped)', S,
                             {'P2_IN': 'det3', 'P2_MID': 'det5', 'P2_OUT': 'det1'})):
        cells = []
        for p in PROBES:
            r = D['REF'][p]
            if not r:
                cells.append('<td>-</td>')
                continue
            s = r['fr']['probe_silent']
            cells.append(f'<td>{chamber[p]}: {pct(r["fr"]["found"])} / <b>{pct(s)} &plusmn; {pct(err(s, r["n"]))}</b>'
                         f'<br><span class="faint">n = {r["n"]:,}</span></td>')
        nr = D['REF']['P2_MID']['runs'] if D['REF']['P2_MID'] else '-'
        rows.append(f'<tr><td class="l">{lab}<br><span class="faint">{nr} REF runs</span></td>{"".join(cells)}</tr>')
    h.append('<h3>10c. Chamber by chamber (REF): found / silent %</h3><div style="overflow-x:auto"><table class="swtab">'
             '<tr><th class="l">period</th><th>P2_IN</th><th>P2_MID</th><th>P2_OUT</th></tr>' + ''.join(rows) +
             '</table></div>')
    sm, so = S['REF']['P2_MID'], S['REF']['P2_OUT']
    nm, no = N['REF']['P2_MID'], N['REF']['P2_OUT']
    h.append(f'<p>Swapping the read-out changes nothing on either chamber: det5 at MID {pct(nm["fr"]["probe_silent"])} '
             f'&rarr; {pct(sm["fr"]["probe_silent"])} % silent, det1 at OUT {pct(no["fr"]["probe_silent"])} &rarr; '
             f'{pct(so["fr"]["probe_silent"])} %. In the funnel det5 at MID is about 4 points worse than det2 was '
             f'(the counting estimator had them equal), and still ~9 points above det1.</p>')
    # per connector
    bM, aM = per_connector(nm, DEFAULT, 'P2_MID'), per_connector(sm, SWAP, 'P2_MID')
    bO, aO = per_connector(no, DEFAULT, 'P2_OUT'), per_connector(so, SWAP, 'P2_OUT')
    h.append('<h3>10d. The decisive view: the same pads, read by other chips</h3><p>Each bar is one half-connector '
             '(64 strips) of a chamber; the number under a bar is the VMM that read it. Grey = before the swap, blue = '
             'after. If the electronics caused the silent muons, the bad bars would move with their VMM to the other '
             'chamber.</p>')
    h.append('<figure class="swfig">' + svg_positions(bM, aM, 'det5 before (VMM 8-13)', 'det5 after (VMM 14-19)') +
             '<figcaption>det5 (P2_MID): silent share per half-connector, REF.</figcaption></figure>')
    h.append('<figure class="swfig">' + svg_positions(bO, aO, 'det1 before (VMM 14-19)', 'det1 after (VMM 8-13)') +
             '<figcaption>det1 (P2_OUT): silent share per half-connector, REF.</figcaption></figure>')
    trs = []
    for v in range(8, 14):
        b = nm['chips'].get(v, (float('nan'), 0))
        af = so['chips'].get(v, (float('nan'), 0))
        trs.append(f'<tr><td>VMM {v}</td><td>{pct(b[0], 0)} % <span class="faint">({b[1]})</span></td>'
                   f'<td>{pct(af[0], 0)} % <span class="faint">({af[1]})</span></td></tr>')
    h.append('<p>The six chips that read P2_MID all along (VMM 8-13), before and after:</p><div style="overflow-x:auto">'
             '<table class="swtab"><tr><th>chip</th><th>on det5 at MID (before)</th><th>on det1 at OUT (after)</th></tr>'
             + ''.join(trs) + '</table></div>')
    # gains
    grs = []
    for g, lab in GAINS:
        cells = []
        for D in (N, S):
            for p in ('P2_MID', 'P2_OUT'):
                r = D[g][p]
                cells.append('<td>-</td>' if not r else
                             f'<td>{pct(r["fr"]["probe_silent"])} <span class="faint">({r["n"]:,})</span></td>')
        grs.append(f'<tr><td class="l">{lab}</td>{"".join(cells)}</tr>')
    h.append('<h3>10e. The other electronics settings (silent %)</h3><div style="overflow-x:auto"><table class="swtab">'
             '<tr><th class="l" rowspan="2">configuration</th><th colspan="2">before swap</th><th colspan="2">after swap'
             '</th></tr><tr><th>det5 MID</th><th>det1 OUT</th><th>det5 MID</th><th>det1 OUT</th></tr>' + ''.join(grs) +
             '</table></div><p class="faint">D50 has 2 runs per period, B50 2 before and 1 after the swap: their errors '
             'are 2-4 times those of REF (n in brackets).</p>')
    # new PC
    sI, nI = S['REF']['P2_IN'], N['REF']['P2_IN']
    t = (f'<p>P2_IN (det3, read-out not touched) moved after the PC change: late hits {pct(nI["fr"]["probe_near_late"])} '
         f'&rarr; {pct(sI["fr"]["probe_near_late"])} %, far {pct(nI["fr"]["probe_far"])} &rarr; '
         f'{pct(sI["fr"]["probe_far"])} %, VMM 2 silent {pct(nI["chips"].get(2, (0, 0))[0], 0)} &rarr; '
         f'{pct(sI["chips"].get(2, (0, 0))[0], 0)} %.</p>')
    if tim:
        t += '<p>' + tim.get('html', '') + '</p>'
    h.append('<h3>10f. The new DAQ PC</h3>' + t)
    h.append('<div class="callout"><b>What the swap shows.</b> The P2_MID silent excess does <b>not</b> follow the '
             'electronics: the chips that read P2_MID (VMM 8-13) are as good as the OUT chips once they read det1, and '
             'det5 keeps its silent share, and its bad region (c3 and c4), with other chips. Hybrids H4-H6 are cleared. '
             'What is left is the chamber or the MID position: both chambers that sat at MID show it (det2 13 %, det5 '
             '17 % vs ~8 % at OUT). One thing only MID has: <b>its HV comes from the iseg NHR</b>, IN and OUT from the '
             'CAEN N1471. Next test: MID on a CAEN channel and OUT on the NHR (no mechanics, 2-3 REF runs), or a '
             'physical MID &harr; OUT chamber swap.</div>')
    return '\n'.join(h)


UPDATE = ('<div class="callout ok"><b>Update 25 Nov: det5 at P2_MID and a read-out swap.</b> Section 10: det5 (CERN '
          'bulk) replaced det2 at P2_MID on 23 Nov, then the P2_MID and P2_OUT read-out (adaptation cards, cables, '
          'hybrids, VMMs) were swapped on 24 Nov. The P2_MID silent excess stays with det5 and its pads, not with the '
          'chips: the read-out is cleared. Remaining suspects: the chamber or the MID position (its HV supply).</div>')
OPEN2 = ('<li><b>P2_MID\'s extra silent share</b>: not the read-out (section 10: the MID &harr; OUT read-out swap moved '
         'nothing). Chamber or position: swap the HV supply channels of MID and OUT (MID is the only chamber on the '
         'iseg NHR), or the chambers themselves.</li>')


def insert_html(doc, frag):
    """section before "Open points" (renumbering), update callout at the top, open point 2 rewritten; idempotent."""
    doc = re.sub(r'\n?<!-- swap:(begin|top) -->.*?<!-- swap:(end|topend) -->\n?', '\n', doc, flags=re.S)
    if '<!-- swap:open2 -->' not in doc:
        doc = doc.replace('<h2>11. References</h2>', '<h2>12. References</h2>').replace(
            '<h2>10. Open points</h2>', '<h2>11. Open points</h2>')
        doc = re.sub(r"<li><b>P2_MID's extra silent share</b>.*?</li>", '<!-- swap:open2 -->' + OPEN2, doc,
                     count=1, flags=re.S)
    i = doc.index('<h2>11. Open points</h2>')
    doc = doc[:i] + '<!-- swap:begin -->\n' + frag + '\n<!-- swap:end -->\n' + doc[i:]
    j = doc.index('<p class="faint">Data: lab cosmic stand')
    return doc[:j] + '<!-- swap:top -->\n' + UPDATE + '\n<!-- swap:topend -->\n' + doc[j:]


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--weekend', required=True)
    ap.add_argument('--night', required=True)
    ap.add_argument('--swap', required=True)
    ap.add_argument('--timing', default=None)
    ap.add_argument('--out', required=True)
    ap.add_argument('--note', default=None)
    a = ap.parse_args()
    frag = section(a)
    open(a.out, 'w').write(frag)
    print('fragment', a.out, len(frag))
    if a.note:
        doc = insert_html(open(a.note).read(), frag)
        open(a.note, 'w').write(doc)
        print('inserted into', a.note)
