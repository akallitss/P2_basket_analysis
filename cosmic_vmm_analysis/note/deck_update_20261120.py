#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Update of the weekly-meeting deck P2_Detector_Performance_Gain_2026-11-18.pptx (20 Nov 2026): a new section
"Where the VMM Inefficiency Comes From" after the summary (muon funnel, non-causes, the electronics gain scan, the two
effects of more VMM gain, the threshold-tail test, late hits, what the silent muons are made of), with three
step-by-step build sequences in the deck's own style, and new backup slides at the end.

The drawing helpers are taken verbatim from the deck's build.py (same template, fonts, colours, boxes, tables).

Usage: python deck_update_20261120.py BUILD_PY OLD_PPTX FIGDIR FUNNEL_SUMMARY GAIN_SCAN_JSON OUT_PPTX
  (python with python-pptx 1.0.x, e.g. the deck build venv)
"""
import ast
import copy
import json
import sys

from lxml import etree
from PIL import Image  # noqa: F401  (used by the helper pic())
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_CONNECTOR, MSO_SHAPE
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.oxml.ns import qn
from pptx.util import Inches, Pt

BUILD, OLD, FIG, FUNNEL, GSJ, OUT = sys.argv[1:7]
FIG = FIG.rstrip('/') + '/'

# ---------------------------------------------------------------- the deck's helpers, verbatim from build.py
src = open(BUILD).read()
tree = ast.parse(src)
WANT = {'_bullet', '_nobullet', 'text', 'pic', 'box', 'arrow', 'caption', 'table', 'stat', '_cell_lines',
        'plain_table', 'outline', 'hline', 'vline', 'label', 'numbered'}
CONST = {'RED', 'DARK', 'NAVY', 'TEAL', 'GREY', 'LIGHT', 'PALERED', 'PALENAVY', 'GREEN', 'WHITE', 'FONT', 'LINEGREY'}
code = []
for node in tree.body:
    if isinstance(node, ast.FunctionDef) and node.name in WANT:
        code.append(ast.get_source_segment(src, node))
    elif isinstance(node, ast.Assign) and any(getattr(t, 'id', None) in CONST for t in node.targets):
        code.append(ast.get_source_segment(src, node))
exec('\n\n'.join(code), globals())
AMBER = RGBColor(0xE8, 0xA3, 0x3D)
NOISEFILL = RGBColor(0xD9, 0xD9, 0xD9)
PALEGREEN = RGBColor(0xE6, 0xF2, 0xE8)

prs = Presentation(OLD)
slides = list(prs.slides)
LAYOUT = [l for l in prs.slide_layouts if l.name == 'Title and Content'][0]
summary = next(s for s in slides if s.shapes.title is not None and s.shapes.title.text.startswith('Summary and Next'))
ref = summary
foot_xml = [copy.deepcopy(sh._element) for sh in ref.shapes
            if sh.is_placeholder and sh.placeholder_format.idx in (11, 12)]
div_src = next(s for s in slides if any(sh.has_text_frame and sh.text_frame.text.strip() == 'Detector Performance'
                                        for sh in s.shapes))
NEW_MAIN, NEW_BACKUP = [], []


def new_slide(title, notes=None, backup=False):
    s = prs.slides.add_slide(LAYOUT)
    for ph in list(s.placeholders):
        if ph.placeholder_format.idx == 1:
            ph._element.getparent().remove(ph._element)
    s.shapes.title.text = title
    for x in foot_xml:
        s.shapes._spTree.append(copy.deepcopy(x))
    if notes:
        s.notes_slide.notes_text_frame.text = notes
    (NEW_BACKUP if backup else NEW_MAIN).append(s)
    return s


def divider(label_text, notes=None):
    s = prs.slides.add_slide(div_src.slide_layout)
    for ph in list(s.placeholders):
        ph._element.getparent().remove(ph._element)
    for sh in div_src.shapes:
        s.shapes._spTree.append(copy.deepcopy(sh._element))
    for sh in s.shapes:
        if sh.has_text_frame and sh.text_frame.text.strip() == 'Detector Performance':
            r = sh.text_frame.paragraphs[0].runs
            r[0].text = label_text
            for extra in r[1:]:
                extra._r.getparent().remove(extra._r)
            sh.width = Inches(11.0)
    if notes:
        s.notes_slide.notes_text_frame.text = notes
    NEW_MAIN.append(s)
    return s


def steps(slide, k, n, x=10.4):
    tb = slide.shapes.add_textbox(Inches(x), Inches(1.0), Inches(2.5), Inches(0.3))
    p = tb.text_frame.paragraphs[0]; p.alignment = PP_ALIGN.RIGHT
    for j in range(1, n + 1):
        r = p.add_run(); r.text = f' {j} '
        r.font.name = FONT; r.font.size = Pt(13)
        r.font.bold = (j == k); r.font.color.rgb = RED if j == k else RGBColor(0xBF, 0xBF, 0xBF)


def rect(slide, x, y, w, h, fill, line=None, name=None):
    r = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(x), Inches(y), Inches(max(w, 0.001)), Inches(max(h, 0.001)))
    r.fill.solid(); r.fill.fore_color.rgb = fill
    if line is None:
        r.line.fill.background()
    else:
        r.line.color.rgb = line; r.line.width = Pt(1.0)
    r.shadow.inherit = False
    if name:
        r.name = name
    return r


# ---------------------------------------------------------------- numbers
FS = json.load(open(FUNNEL))['groups']
GE = json.load(open(GSJ))
NUM = json.load(open(FIG + 'deck_numbers.json'))
REF = FS['REF']
on = REF['probes']['P2_MID']['on']
n0 = REF['n_scint']
n1 = n0 - on['no_tag_both'] - on['no_tag_A'] - on['no_tag_B']
n2 = n1 - on['tag_rejected']
n3 = n2 - on.get('outside_area', 0) - on.get('dead_pad', 0)
nf = on['found']
pc = lambda a, b: 100.0 * a / b
sil = pc(on['probe_silent'], n3)
far = pc(on['probe_far'], n3)
late = pc(on.get('probe_near_late', 0), n3)
fnd = pc(nf, n3)
PARTS = NUM['silent_parts']
runs_ref = ', '.join(map(str, sorted(set(REF['runs']))))
X = REF['xtalk']
xsat = {d: [(X[d]['on'][k][10] - X[d]['off'][k][10]) / X[d]['n'][10] * 100 for k in range(3)] for d in X}
conf = {(r['det'], r['label']): r for r in GE['configs']}

# ================================================================ the new section
divider('Update 20 Nov: Where the VMM Inefficiency Comes From',
        'New since the 18 Nov meeting: the weekend cosmic runs (19-20 Nov, runs 268-307 without the pulser runs), '
        'every scintillator coincidence followed through the efficiency measurement, and the VMM gain x noise-rate '
        'scan. Full write-up with event displays: CERN site note "Where do the missed muons go?" (unlisted).')

# ---------------------------------------------------------------- the question
s = new_slide('The Question: Electronics or the Way We Count?',
              'VMM-read P2 chambers reach 0.85-0.90 on the cosmic stand; DREAM reaches 0.93-0.98 on the same chambers '
              '(M3 bench, connector-matched). Before asking the detector for more gain we checked that the gap is not an '
              'artefact of the efficiency calculation, and whether the VMM settings (gain, thresholds) can close it.')
stat(s, 0.6, 1.5, 3.8, '~80 %', 'P2_MID with VMM: muons found\n(cosmic stand, REF settings)')
stat(s, 4.7, 1.5, 3.8, '93–94 %', 'DREAM on the same pads\n(M3 bench, any hit)', color=NAVY)
text(s, 0.6, 3.4, 12.2, 3.2, [
    dict(t='Two questions for this update', head=True, size=16),
    dict(t='Is the missing 10–15 % real, or lost in how we count (read-out, dead time, tracking, cuts)?'),
    dict(t='If it is real: can the VMM settings recover it, by raising the electronics gain or lowering the '
           'thresholds, or do we need more gain from the detector?'),
    dict(t='Method', head=True, size=16),
    dict(t='Start from every scintillator coincidence (a muon certainly crossed the stand) and follow it through '
           'every step of the measurement; then compare the 13 settings of the weekend gain × noise-rate scan'),
], size=15, space=6)

# ---------------------------------------------------------------- animation A1: the funnel
FUN = [(n0, 'scintillator coincidences', 'a muon crossed the paddles: our reference'),
       (n1, 'both tag chambers fired', 'P2_IN and P2_OUT both have a hit 140–310 ns later'),
       (n2, 'two clean tags', 'each tag is one compact cluster (≤ 6 pads within 15 mm)'),
       (n3, 'expected on a live pad of P2_MID', 'the line through the tags points at a working pad'),
       (nf, 'found by P2_MID', 'a pad within 15 mm fired, in time')]
EXPL = [
    [dict(t='1. Start from the scintillators', head=True, size=16),
     dict(t=f'{n0:,} coincidences in the REF runs ({runs_ref})'),
     dict(t='P2_IN and P2_OUT are the "tags": they tell us where the muon went'),
     dict(t=f'{pc(n0 - n1, n0):.0f} % are lost here, mostly geometry: the paddles are larger than the read-out '
            'area (3 of 6 connectors), so many muons miss a tag chamber', c=GREY)],
    [dict(t='2. Clean tags only', head=True, size=16),
     dict(t=f'{pc(on["tag_rejected"], n1):.0f} % of the muons with both tags firing are dropped: one tag has an extra '
            'pad far away in the same 100 ns'),
     dict(t='These extra pads are crosstalk into the neighbouring VMM channel (see backup)'),
     dict(t='Dropped muons are not counted as missed: they leave the sample', c=GREY)],
    [dict(t='3. A predicted point', head=True, size=16),
     dict(t='Straight line through the two tag clusters → a point in P2_MID'),
     dict(t=f'{on.get("outside_area", 0) + on.get("dead_pad", 0)} muons point outside the read-out area or at a '
            'dead pad and are removed'),
     dict(t=f'{n3:,} muons are left: the denominator of the efficiency', b=True)],
    [dict(t='4. The probe', head=True, size=16),
     dict(t=[(f'{fnd:.1f} % found', {'b': True, 'c': GREEN})]),
     dict(t=[(f'{sil:.1f} % silent: ', {'b': True, 'c': RED}), ('P2_MID records no hit at all, anywhere', {})]),
     dict(t=f'{far:.1f} % fire, but more than 15 mm away; {late:.1f} % fire near but late'),
     dict(t='Only this step is the efficiency; the silent muons are the main loss', b=True)],
]
for k in range(1, 5):
    s = new_slide('Following Every Muon, Step by Step',
                  'Reference configuration REF (gain 4.5 mV/fC, trimmed thresholds Km18), P2_MID as probe. Every '
                  'coincidence ends in exactly one class. Code: cosmic_vmm_analysis/muon_funnel.py (same cuts as the '
                  'tag-and-probe, efficiency_tp.py).')
    steps(s, k, 4)
    x0, w0, y = 0.6, 6.4, 1.55
    for j, (n, lab, sub) in enumerate(FUN[:k + 1]):
        w = w0 * n / n0
        fill = GREEN if j == 4 else (NAVY if j == k else RGBColor(0xA9, 0xAE, 0xC9))
        rect(s, x0, y, w, 0.5, fill)
        wl = w0 * n3 / n0 if j == 4 else w      # the found row carries the loss strips after it
        if w > 3.6:
            label(s, x0 + 0.08, y + 0.08, w - 0.1, f'{n:,}  {lab}', size=13, color=WHITE, bold=True)
        else:
            label(s, x0 + wl + 0.1, y + 0.08, 5.6, f'{n:,}  {lab}', size=13, color=DARK, bold=True)
        label(s, x0, y + 0.53, 6.6, sub, size=10, color=GREY)
        y += 0.98
    if k == 4:
        xs = x0
        for nm, val, col in (('silent', on['probe_silent'], RED), ('far', on['probe_far'], AMBER),
                             ('late', on.get('probe_near_late', 0), NAVY)):
            ww = w0 * val / n0
            rect(s, x0 + w0 * nf / n0 + (xs - x0), 1.55 + 4 * 0.98, ww, 0.5, col)
            xs += ww
    text(s, 7.6, 1.45, 5.3, 5.3, EXPL[k - 1], size=14, space=6)

# ---------------------------------------------------------------- not the read-out, dead time, sparks
s = new_slide('Not the Read-out, Not Dead Time, Not Sparks',
              'Read-out: test pulses on one channel per chip during cosmics (run 270, gain 4.5): 99.82-99.87 % '
              'received on all six chips. Dead time: silent share vs time since the previous hit on the predicted '
              'chip, excess only below 10 us = 83 of 28 480 muons. Sparks: muons 1-100 ms after a spark of >= 3 '
              'saturated hits. Details: WRAPUP_20261119.md and the CERN-site note.')
pic(s, FIG + 'c_budget.png', 0.5, 1.3, w=7.2, name='Fig budget')
text(s, 8.0, 1.3, 4.9, 5.4, [
    dict(t='Three ways to lose a muon in the electronics', head=True, size=16),
    dict(t='Read-out drops data: test pulses arrive 99.8 % of the time → 0.14 points'),
    dict(t='Dead time after a hit: only below 10 µs → 0.29 points'),
    dict(t='Sparks blind the chamber: 1–100 ms after a spark → 0.12 points'),
    dict(t='Together 0.55 points, against 13.8 % silent muons in P2_MID', b=True),
    dict(t='→ the silent muons are not lost data: the chamber simply records nothing', b=True, bullet=False, c=RED),
], size=14, space=6)

# ---------------------------------------------------------------- the efficiency step alone
s = new_slide('The Efficiency Step Alone',
              'Muons expected on a live pad of the probe (REF). Found = a live pad within 15 mm of the prediction '
              'fired within +-100 ns of the tag hit. Silent = no hit of the probe chamber at all 140-310 ns after the '
              'scintillator. Far = hits, none within 15 mm. Late = within 15 mm but outside the tag window. Random '
              'hits (same test 4 us later): 0.1-0.2 %.')
pic(s, FIG + 'c_split.png', 0.5, 1.3, w=7.2, name='Fig split')
text(s, 8.0, 1.3, 4.9, 5.4, [
    dict(t='Four outcomes', head=True, size=16),
    dict(t=[('Silent ', {'b': True, 'c': RED}), ('12–14 % in all three chambers: the main loss, and it does not '
                                                  'depend on any position or cut', {})]),
    dict(t=[('Far ', {'b': True}), ('6 % for P2_MID, 14–16 % for P2_IN/P2_OUT: their prediction is an '
                                    'extrapolation, so less precise', {})]),
    dict(t=[('Late ', {'b': True}), ('< 1 %', {})]),
    dict(t='From now on we use the silent share: it is the cleanest measure of what the electronics loses', b=True),
], size=14, space=6)

# ---------------------------------------------------------------- animation A2: raising the VMM gain
SIG3 = [2.4, 3.1, 2.6, 3.6, 2.9, 3.3, 2.5, 3.9]        # mV noise sigma at gain 3, illustrative
OFF9 = [-6, 4, -2, 14, 0, 9, -8, 20]                     # extra spread at gain 9 (mV), illustrative
AX0, AX1, AY0, AY1, VMAX = 1.3, 7.4, 6.15, 1.6, 300.0


def ya(v):
    return AY0 - min(v, VMAX) / VMAX * (AY0 - AY1)


def xc(i):
    return AX0 + 0.5 + i * 0.75


def gain_frame(s, k):
    vline(s, AX0, AY1 - 0.05, AY0, DARK, 1.0)
    hline(s, AX0, AX1, AY0, DARK, 1.0)
    for v in (0, 50, 100, 150, 200, 250, 300):
        hline(s, AX0 - 0.06, AX0, ya(v), DARK, 0.75)
        label(s, AX0 - 0.6, ya(v) - 0.12, 0.5, str(v), size=10, color=GREY, align=PP_ALIGN.RIGHT)
    lb = text(s, -1.05, 3.75, 2.6, 0.3, [dict(t='mV above the baseline', bullet=False, c=GREY, size=11)],
              align=PP_ALIGN.CENTER)
    lb.rotation = 270
    for i in range(8):
        label(s, xc(i) - 0.2, AY0 + 0.03, 0.4, str(i + 1), size=10, color=GREY, align=PP_ALIGN.CENTER)
    label(s, AX0, AY0 + 0.28, AX1 - AX0, 'channels of one chip (8 shown, illustrative values)', size=11, color=GREY,
          align=PP_ALIGN.CENTER)
    steps(s, k, 4)


def channels_g(s, scale, thr_scale, extra=None, clip=None):
    """noise band (2 sigma), the threshold each channel needs (red tick), optionally clipped at 'clip'."""
    for i, sg in enumerate(SIG3):
        x = xc(i)
        top = 2 * sg * scale
        rect(s, x - 0.2, ya(top), 0.4, ya(0) - ya(top), NOISEFILL)
        need = 3.3 * sg * thr_scale + (extra[i] if extra else 0)
        got = max(need, clip) if clip is not None else need
        hline(s, x - 0.27, x + 0.27, ya(got), RED, 2.25)
        if clip is not None and got > need + 0.5:
            hline(s, x - 0.2, x + 0.2, ya(need), RED, 1.0, dash=4)
            vline(s, x + 0.31, ya(got), ya(need), RED, 1.0, arrow_both=True)


def signal_line(s, v, txt):
    hline(s, AX0, AX1, ya(v), GREEN, 2.0, dash=7)
    label(s, AX0 + 0.08, ya(v) - 0.3, 6.0, txt, size=11, color=GREEN, bold=True)


def legend_gain(s):
    y = AY1 - 0.35
    rect(s, AX0 + 0.1, y + 0.06, 0.25, 0.16, NOISEFILL)
    label(s, AX0 + 0.4, y, 1.6, 'noise (± 2σ)', size=10, color=GREY)
    hline(s, AX0 + 1.9, AX0 + 2.2, y + 0.14, RED, 2.25)
    label(s, AX0 + 2.25, y, 3.2, 'threshold for the same noise rate', size=10, color=GREY)


GAIN_EXPL = [
    [dict(t='1. Gain 3 mV/fC', head=True, size=16),
     dict(t='Each channel sits where its own noise rate is R (e.g. 10 Hz): ~3.5 fC, about 10 mV'),
     dict(t='A typical muon gives ~8× more on its pad'),
     dict(t='The threshold is at ~12 % of the typical muon signal', b=True)],
    [dict(t='2. Gain 9: what we hoped', head=True, size=16),
     dict(t='Three times the gain: the muon pulse is 3× taller'),
     dict(t='If the noise only scaled with the gain, the thresholds would also be 3× higher'),
     dict(t='Same ratio threshold / signal → no gain at all', b=True),
     dict(t='More electronics gain only helps if the noise grows less than the signal', c=GREY)],
    [dict(t='3. Gain 9: what we measured', head=True, size=16),
     dict(t='The muon pulse grows 3× (measured, the gain is real)'),
     dict(t='But for the same noise rate the same channel needs 1.5× more charge: the noise grows faster'),
     dict(t='Threshold / signal goes from 12 % to ~18 %: worse', b=True, c=RED)],
    [dict(t='4. And the trims run out', head=True, size=16),
     dict(t='One threshold per chip + a 5-bit trim per channel: a trim can lower a channel by ~34 mV at most'),
     dict(t='At gain 9 the channels of a chip need thresholds spread over 40–65 mV'),
     dict(t='Channels beyond the trim range stay higher ("clipped"): 360–708 of 1152'),
     dict(t='Measured: clipped channels are not more silent; effect 3 is the one that costs', c=GREY)],
]
for k in range(1, 5):
    s = new_slide('Raising the VMM Gain, Step by Step',
                  'Illustrative channels. Measured inputs: real gain from the muon pulse height (gain 6 / gain 3 = '
                  '1.8-2.2, gain 9 / gain 3 = 2.8-3.2); same-channel noise-rate target in charge: gain 6 / gain 3 = '
                  f'{NUM["noise_ratio_g6"]:.2f}, gain 9 / gain 3 = {NUM["noise_ratio_g9"]:.2f} (median, R = 10 Hz); '
                  'per-chip spread of the targets: median 10 mV (gain 3), 20 mV (gain 6), 39 mV (gain 9), trim range '
                  '31 x 1.1 mV. Sources: bscan_target reports of 18 Nov, muon_adc_gain.py, trim_reports.py.')
    gain_frame(s, k); legend_gain(s)
    if k == 1:
        channels_g(s, 1.0, 1.0)
        signal_line(s, 85, 'typical muon signal at gain 3')
    elif k == 2:
        channels_g(s, 3.0, 3.0)
        signal_line(s, 255, 'typical muon signal at gain 9 (3×)')
    elif k == 3:
        channels_g(s, 4.5, 4.5)
        signal_line(s, 255, 'typical muon signal at gain 9 (3×)')
    else:
        need = [3.3 * sg * 4.5 + o for sg, o in zip(SIG3, OFF9)]
        top = max(need)
        clip = top - 34
        rect(s, AX0 + 0.1, ya(top), AX1 - AX0 - 0.2, ya(clip) - ya(top), PALEGREEN)
        label(s, AX1 - 3.3, ya(top) - 0.32, 3.2, 'trim range of this chip (~34 mV)', size=10, color=GREEN,
              bold=True, align=PP_ALIGN.RIGHT)
        channels_g(s, 4.5, 4.5, extra=OFF9, clip=clip)
        signal_line(s, 255, 'typical muon signal at gain 9 (3×)')
    text(s, 7.75, 1.45, 5.15, 5.3, GAIN_EXPL[k - 1], size=14, space=6)

# ---------------------------------------------------------------- measured: the gain is real
s = new_slide('Measured: The Gain Is Real',
              'Median ADC of the muon pad (largest pad within 15 mm and +-150 ns of the prediction, scintillator-'
              'confirmed tags) minus a 40 ADC baseline, divided by the same at gain 3, for every gain-scan '
              'configuration. Code: muon_adc_gain.py.')
pic(s, FIG + 'c_gain_real.png', 0.5, 1.3, w=6.6, name='Fig gain real')
text(s, 7.5, 1.3, 5.4, 5.4, [
    dict(t='Before blaming the electronics: is the setting what we think?', head=True, size=16),
    dict(t='The muon pulse height grows exactly like the nominal gain: 2× at gain 6, 3× at gain 9'),
    dict(t='Same in all three chambers'),
    dict(t='→ the settings were applied, and converting thresholds to charge with the nominal gain is right',
         b=True, bullet=False),
], size=14, space=6)

# ---------------------------------------------------------------- effect 1
s = new_slide('Effect 1: More Gain, More Noise in Charge',
              'bscan_target.py reports of 18 Nov (HV-on noise scans): for every channel the threshold where its noise '
              'rate equals R, converted to charge with the nominal gain. Only channels with a measured crossing at '
              'both gains (765 for gain 6 / 3, 702 for gain 9 / 3 at R = 10 Hz).')
pic(s, FIG + 'c_noise_fC.png', 0.5, 1.3, w=6.8, name='Fig noise fC')
text(s, 7.6, 1.3, 5.3, 5.4, [
    dict(t='The same channel, the same noise rate', head=True, size=16),
    dict(t='At ratio 1 the threshold in charge would not change: more gain would leave more room below the muons'),
    dict(t=f'Measured: {NUM["noise_ratio_g6"]:.2f}× more charge at gain 6, {NUM["noise_ratio_g9"]:.2f}× at gain 9'),
    dict(t='The noise grows faster than the signal: raising the VMM gain does not open room below the muon signal',
         b=True),
], size=14, space=6)

# ---------------------------------------------------------------- effect 2
s = new_slide('Effect 2: The Trim Range Runs Out at Gain 9',
              'Left: per chip, max - min of the channels\' noise-rate targets (R = 10 Hz, channels with a measured '
              'crossing). Right: silent share of the muons whose predicted pad is a clipped channel, a channel '
              'without measured crossing (given the chip median) or a normal one. Code: trim_reports.py, '
              'silent_vs_clipped.py.')
pic(s, FIG + 'c_trim_spread.png', 0.4, 1.3, w=6.1, name='Fig trim spread')
pic(s, FIG + 'c_clipped.png', 6.7, 1.3, w=6.2, name='Fig clipped')
text(s, 0.5, 5.35, 12.3, 1.4, [
    dict(t='At gain 9 the needed thresholds of one chip spread over more than a trim can correct → 360–708 channels '
           'stay higher than their target'),
    dict(t='But their muons are not more often silent than on the well-trimmed channels: clipping is not what makes '
           'gain 9 worse', b=True),
    dict(t='Channels whose noise curve never crossed R (given the chip median) are worse at gain 6 R = 2 Hz and gain 9 '
           'R = 10 Hz: a calibration gap to close', c=GREY),
], size=14, space=5)

# ---------------------------------------------------------------- together: silent vs threshold
s = new_slide('Together: Silent Muons vs Threshold in Charge',
              'P2_MID, every weekend setting: gain 3 / 4.5 / 6 / 9 mV/fC x R = 2 / 10 / 50 Hz, the gain-4.5 threshold '
              'scan (open circles) and REF (square). Threshold in charge = (threshold - baseline) / nominal gain, '
              'median over the chamber. Fits: +2.9 +- 0.2 points per fC over all settings, +1.3 +- 0.8 inside the '
              'reachable band. P2_IN and P2_OUT in the backup.')
pic(s, FIG + 'c_silent_thr_mid.png', 0.4, 1.3, w=7.6, name='Fig silent vs thr')
text(s, 8.2, 1.3, 4.7, 5.4, [
    dict(t='What the electronics can reach', head=True, size=16),
    dict(t='Gains 3, 4.5 and 6 all land at 3.2–4.4 fC (green band) and lose the same: 12.4–14.5 %'),
    dict(t='25× more noise (R = 2 → 50 Hz) lowers the threshold by ~0.5 fC and wins ≤ 1.4 points'),
    dict(t='Gain 9 sits at 6.7–8 fC and loses 22–25 %'),
    dict(t='→ no setting brings P2_MID above ~87 % "any hit"', b=True, bullet=False, c=RED),
], size=14, space=6)

# ---------------------------------------------------------------- the grid
s = new_slide('Efficiency over the Gain × Noise-Rate Grid',
              'P2_MID: share of the muons expected on a live pad for which the chamber records any hit (= 1 - silent), '
              'and the threshold in charge of that setting. Gain 4.5 column = the threshold scan (bscan); REF (gain '
              '4.5, Km18) gives 86.2 % at 3.9 fC.')
pic(s, FIG + 'c_grid.png', 0.5, 1.3, w=7.4, name='Fig grid')
text(s, 8.2, 1.3, 4.7, 5.4, [
    dict(t='How to read it', head=True, size=16),
    dict(t='Columns: electronics gain; rows: how much noise we accept per channel'),
    dict(t='Gains 3 to 6: 83–88 % everywhere, best at 50 Hz'),
    dict(t='Gain 9: 75–78 %, whatever the noise rate'),
    dict(t='→ the electronics knobs move the efficiency by a few points at most; the best corner is ~88 %',
         b=True, bullet=False),
], size=14, space=6)

# ---------------------------------------------------------------- is it a threshold tail?
pm = PARTS['REF (gain 4.5)']
s = new_slide('Is the 13 % Silent a Threshold Tail?',
              'Largest pad ADC of the found muons, full ADC resolution with the period-16 ADC non-linearity corrected '
              'from the data, divided by the median muon signal (both minus a 40 ADC baseline). Threshold edge = 1st '
              'percentile. Below-edge extrapolation: density flat down to zero signal (upper) or falling linearly to '
              'zero (lower). Code: threshold_spectrum.py. Charge sharing: 93-97 % of found muons fire a single pad; '
              'two-pad clusters are not closer to the edge.')
pic(s, FIG + 'c_spectrum.png', 0.4, 1.3, w=7.2, name='Fig spectrum')
text(s, 7.9, 1.3, 5.0, 5.4, [
    dict(t='If the silent muons were a tail below threshold…', head=True, size=16),
    dict(t='…the spectrum of the found muons would run into the threshold edge (dotted) and continue below it'),
    dict(t='The edge sits at 12–19 % of the typical muon signal (P2_MID)'),
    dict(t=f'Extending the spectrum below the edge gives {pm["tail_lo"]:.0f}–{pm["tail_hi"]:.0f} % of the muons '
           f'(REF), against {pm["silent"]:.1f} % silent'),
    dict(t='→ only part of the silent muons is a threshold tail; most leave no signal above threshold at all',
         b=True, bullet=False),
    dict(t='Charge sharing between pads is not the cause (93–97 % single-pad muons)', c=GREY),
], size=14, space=6)

# ---------------------------------------------------------------- late hits at gain 9
s = new_slide('Gain 9: Some "Silent" Muons Are Late',
              'Silent muons (no probe hit within +-100 ns of the tag time): share with a live probe pad within 15 mm of '
              'the prediction in later (or earlier) time bins, minus the random rate from -2..-1 us. Code: '
              'silent_late_hits.py.')
pic(s, FIG + 'c_late.png', 0.4, 1.3, w=7.2, name='Fig late')
text(s, 7.9, 1.3, 5.0, 5.4, [
    dict(t='Some of them do fire, later', head=True, size=16),
    dict(t='REF: ~3 % of the silent muons have a hit near the prediction 100 ns – 1 µs later'),
    dict(t='Gain 9: 14 % (P2_MID), 25 % (P2_IN)'),
    dict(t='The found muons are on time at every gain: it is a late tail, only at high gain'),
    dict(t='→ about 4 of the ~10 extra silent points of gain 9', b=True, bullet=False),
], size=14, space=6)

# ---------------------------------------------------------------- animation A3: what the silent muons are made of
GROUPS = list(PARTS)
BX0, BW, BY0, BYT, PMAX = 1.55, 1.15, 6.1, 1.7, 26.0


def yb(v):
    return BY0 - v / PMAX * (BY0 - BYT)


def parts_frame(s, k):
    vline(s, BX0 - 0.2, BYT - 0.05, BY0, DARK, 1.0)
    hline(s, BX0 - 0.2, BX0 + 4 * (BW + 0.35), BY0, DARK, 1.0)
    for v in (0, 5, 10, 15, 20, 25):
        hline(s, BX0 - 0.26, BX0 - 0.2, yb(v), DARK, 0.75)
        label(s, BX0 - 0.8, yb(v) - 0.12, 0.5, f'{v}', size=10, color=GREY, align=PP_ALIGN.RIGHT)
    lb = text(s, -0.95, 3.75, 2.5, 0.3, [dict(t='P2_MID silent muons (%)', bullet=False, c=GREY, size=11)],
              align=PP_ALIGN.CENTER)
    lb.rotation = 270
    steps(s, k, 3)


PART_EXPL = [
    [dict(t='1. How many are silent', head=True, size=16),
     dict(t='Gains 3, 4.5 and 6: 13–14 %; gain 9: 24 %'),
     dict(t='What are these muons?', b=True)],
    [dict(t='2. A threshold tail', head=True, size=16),
     dict(t=[('Orange: ', {'b': True, 'c': AMBER}), ('the part the muon spectrum predicts below the threshold edge '
                                                     '(range from two extrapolations): 3–5 points', {})]),
     dict(t='More detector gain would recover this part', c=GREY)],
    [dict(t='3. Late, and the rest', head=True, size=16),
     dict(t=[('Navy: ', {'b': True, 'c': NAVY}), ('a hit near the prediction, but outside the time window: ~0.5 '
                                                  'point, 4 points at gain 9', {})]),
     dict(t=[('Red: ', {'b': True, 'c': RED}), ('no hit at all and not explained by a smooth tail: 8–10 points '
                                                '(16 at gain 9)', {})]),
     dict(t='This is the part to understand next: muons with no signal above threshold on any pad', b=True)],
]
for k in range(1, 4):
    s = new_slide('What the Silent Muons Are Made Of',
                  'P2_MID. Silent share from the funnel (n-weighted over the settings of each group: gain 3 = R 2/10/50, '
                  'gain 6 = R 10/50, gain 9 = R 2/10/50). Tail = middle of the two below-edge extrapolations '
                  '(error bar = their range). Late = share of silent muons with a near hit -1 us..-100 ns or '
                  '100 ns..1 us, random-subtracted. The silent definitions differ slightly (funnel: scintillator '
                  'window; late: tag window).')
    parts_frame(s, k)
    for i, g in enumerate(GROUPS):
        p = PARTS[g]
        x = BX0 + i * (BW + 0.35)
        label(s, x - 0.2, BY0 + 0.05, BW + 0.4, g, size=11, color=DARK, align=PP_ALIGN.CENTER)
        if k == 1:
            rect(s, x, yb(p['silent']), BW, yb(0) - yb(p['silent']), RGBColor(0xC8, 0xC8, 0xC8))
        else:
            rect(s, x, yb(p['tail']), BW, yb(0) - yb(p['tail']), AMBER)
            if k == 2:
                rect(s, x, yb(p['silent']), BW, yb(p['tail']) - yb(p['silent']), RGBColor(0xC8, 0xC8, 0xC8))
            else:
                rect(s, x, yb(p['tail'] + p['late']), BW, yb(p['tail']) - yb(p['tail'] + p['late']), NAVY)
                rect(s, x, yb(p['silent']), BW, yb(p['tail'] + p['late']) - yb(p['silent']), RED)
            vline(s, x + BW / 2, yb(p['tail_hi']), yb(p['tail_lo']), DARK, 1.25, arrow_both=False)
        label(s, x - 0.2, yb(p['silent']) - 0.32, BW + 0.4, f'{p["silent"]:.1f} %', size=12, color=DARK, bold=True,
              align=PP_ALIGN.CENTER)
    text(s, 7.75, 1.45, 5.15, 5.3, PART_EXPL[k - 1], size=14, space=6)

# ---------------------------------------------------------------- side checks
s = new_slide('Side Checks: Tags, Crosstalk, Paddles',
              'Crosstalk: for every pad of a muon cluster, does the channel c+-1/2/3 of the same chip fire within '
              '+-100 ns (chance rate from +4 us: ~0). Paddles: line through the two tag clusters extrapolated to the '
              'top (z -100..0) and bottom (1225..1275 mm) paddles, positions from the 4 Nov fit of the October tag '
              'maps (not measured). See backup.')
rows = [['Check', 'Result', 'Consequence'],
        ['Rejected tag groups', 'extra single pad, small ADC, in time, 1–3 channels from a muon channel '
                                '(60–65 %, chance 10–17 %)', 'crosstalk: costs tags, not efficiency'],
        ['Crosstalk vs signal', f'saturated muon pad: neighbour ± 1 fires {min(v[0] for v in xsat.values()):.0f}–'
                                f'{max(v[0] for v in xsat.values()):.0f} %, ± 2 {min(v[1] for v in xsat.values()):.0f}–'
                                f'{max(v[1] for v in xsat.values()):.0f} %; ~1 % for medium signals',
         'agrees with the pulser (~3.5 % of the pulse)'],
        ['Tag line back to the paddles', '13 % of the lines miss a paddle; P2_MID 79.6 % → 80.6 % found if removed',
         'cleaning step, silent share unchanged'],
        ['Chance hits', 'same test 4 µs later: 0.1–0.2 %', 'negligible']]
table(s, 0.5, 1.4, 12.3, [2.6, 6.1, 3.6], rows, size=12, row_h=0.62, name='Side checks')
text(s, 0.5, 5.3, 12.3, 1.2, [
    dict(t='None of these changes the silent share: the 13 % is a property of the chamber + VMM, not of the counting',
         b=True),
], size=14)

# ---------------------------------------------------------------- summary of the update
s = new_slide('Summary of the Update and Next Steps',
              'All numbers: REF = gain 4.5 Km18, runs ' + runs_ref + '; gain scan runs 268-307 (no pulser runs). '
              'CERN-site note with event displays, the funnel per chamber and a reference section with every data and '
              'code path.')
text(s, 0.5, 1.3, 6.1, 5.4, [
    dict(t='What we learned', head=True, size=16),
    dict(t='Counting is not the problem: read-out, dead time and sparks cost 0.55 points; P2_MID records no hit '
           'at all for 13–14 % of the muons that certainly crossed it'),
    dict(t='The VMM gain is real, but the noise in charge grows faster than the signal: more electronics gain '
           'does not help (gain 9 is worse: 24 % silent)'),
    dict(t='Within the noise we can accept, the thresholds move the efficiency by ≤ 1.5 points: best ~87–88 % '
           'any hit for P2_MID'),
    dict(t='Only 3–5 points of the silent muons look like a threshold tail; 8–10 points leave no signal above '
           'threshold at all; at gain 9 part of them come late', b=True),
], size=14, space=6)
text(s, 6.85, 1.3, 6.0, 5.4, [
    dict(t='Next steps', head=True, size=16),
    dict(t=[('Mesh-HV scan of P2_MID with REF electronics: ', {'b': True}),
            ('measures directly how fast the silent share falls with detector gain', {})]),
    dict(t=[('The no-signal muons: ', {'b': True}),
            ('where are they (map, chip, time), and does DREAM lose the same muons?', {})]),
    dict(t=[('Gain 9 late hits: ', {'b': True}), ('what delays them (large / saturated pulses?)', {})]),
    dict(t=[('Calibration: ', {'b': True}),
            ('channels without a measured noise crossing; crosstalk filter for the tags', {})]),
    dict(t='Not now: the VMM transfer function (test-pulse amplitude scan)', c=GREY),
    dict(t='For discussion: is ~88 % the VMM ceiling for these chambers at this HV?', b=True, bullet=False),
], size=14, space=6)

# ================================================================ new backups (at the end)
s = new_slide('Backup – Silent Muons vs Threshold, All Chambers',
              'Same as the main slide for P2_IN, P2_MID and P2_OUT.', backup=True)
pic(s, FIG + 'c_silent_thr_all.png', 0.9, 1.3, w=11.4)
caption(s, 0.5, 6.25, 12.3, 'Fits (silent vs threshold in charge): P2_IN +2.8, P2_MID +2.9, P2_OUT +2.4 points per fC '
        'over all settings; inside the reachable band +1.3 (P2_MID, P2_OUT), P2_IN noisy.', size=12)

s = new_slide('Backup – The Weekend Gain Scan, Numbers (P2_MID)',
              'From gain_scan_eff.py over the muon funnel; found includes the position and time cuts, silent does not.',
              backup=True)
rows = [['Setting', 'Runs', 'Threshold (fC)', 'Muons', 'Silent', 'Found']]
for r in sorted([r for r in GE['configs'] if r['det'] == 'P2_MID'], key=lambda r: r['thr_fC']):
    rows.append([r['label'], str(r['runs']), f'{r["thr_fC"]:.1f}', f'{r["n"]:,}', f'{100 * r["silent"]:.1f} %',
                 f'{100 * r["found"]:.1f} %' + (' (positions suspect)' if r['family'] == 'bscan' and r['far'] > 0.15
                                                else '')])
plain_table(s, 0.6, 1.3, [3.6, 0.8, 1.6, 1.2, 1.2, 2.9], rows, size=11, row_h=0.36)

s = new_slide('Backup – Charge Sharing and Pad Position',
              'threshold_spectrum.py (cluster size of the found muons), silent_vs_pad_position.py.', backup=True)
text(s, 0.5, 1.3, 12.3, 5.4, [
    dict(t='Found muons: 93–97 % fire one pad, 3–7 % two pads; the largest pad of two-pad clusters is as large as '
           'for one-pad ones, and they sit less often near the threshold edge (3–11 % vs 10–13 %)'),
    dict(t='Silent share vs where the prediction falls inside its pad: flat at REF and gain 3 (± 2–3 points); slight rise '
           'towards the edge at gain 9 (P2_MID 21 → 26 %)'),
    dict(t='Limitation: most tags are single pads, so the predictions sit on a grid of pad centres', c=GREY),
    dict(t='→ charge shared between pads does not explain the silent muons', b=True),
], size=15, space=8)

s = new_slide('Backup – Where Every Muon Ends (All Three Probes)',
              'Same funnel with P2_OUT (tags IN + MID) and P2_IN (tags MID + OUT) as probe: their predictions are '
              'extrapolations, so more fall outside the read-out area and more are FAR.', backup=True)
rows = [['Probe', 'Coincidences', 'Both tags', 'Clean tags', 'On a live pad', 'Found', 'Silent', 'Far']]
for d in ('P2_IN', 'P2_MID', 'P2_OUT'):
    o = REF['probes'][d]['on']
    a1 = n0 - o['no_tag_both'] - o['no_tag_A'] - o['no_tag_B']
    a2 = a1 - o['tag_rejected']
    a3 = a2 - o.get('outside_area', 0) - o.get('dead_pad', 0)
    rows.append([d, f'{n0:,}', f'{a1:,}', f'{a2:,}', f'{a3:,}', f'{pc(o["found"], a3):.1f} %',
                 f'{pc(o["probe_silent"], a3):.1f} %', f'{pc(o["probe_far"], a3):.1f} %'])
table(s, 0.6, 1.5, 12.1, [1.4, 1.6, 1.5, 1.5, 1.7, 1.4, 1.4, 1.6], rows, size=12, row_h=0.45)

# ================================================================ title + order
for sh in slides[0].shapes:
    if sh.has_text_frame and sh.text_frame.text.strip() == 'Wednesday, 18 November 2026':
        r = [r for p in sh.text_frame.paragraphs for r in p.runs if r.text.strip()]
        r[-1].text = 'Wednesday, 18 November 2026 – updated 20 November 2026'
lst = prs.slides._sldIdLst
ids = list(lst)
by_slide = {id(s_._element): sid for s_, sid in zip(prs.slides, ids)}
main_ids = [by_slide[id(s_._element)] for s_ in NEW_MAIN]
anchor = by_slide[id(summary._element)]
for sid in main_ids:
    lst.remove(sid)
pos = list(lst).index(anchor) + 1
for k, sid in enumerate(main_ids):
    lst.insert(pos + k, sid)
prs.save(OUT)
print('saved', OUT, len(prs.slides), 'slides;', len(NEW_MAIN), 'new main,', len(NEW_BACKUP), 'new backup')
