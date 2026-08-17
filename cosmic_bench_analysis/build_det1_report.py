#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""build_det1_report.py -- DET1_EFFICIENCY.md as a PDF, with the figures that
support each argument.

    python3 build_det1_report.py [-o OUT.pdf]

This is deliberately a CURATED figure set, not a dump of the ~245 products in
the det1 tree: one figure per claim the document makes, each with a caption
saying what it is evidence for. The full product tree stays the reference.

Text is read from Analysis/DET1_EFFICIENCY.md so the PDF cannot drift from the
markdown. A missing figure is reported and skipped, never silently dropped.
"""

import argparse
import datetime
import os
import re
import textwrap

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.backends.backend_pdf import PdfPages

import p2_qa_config as qa

INK = '#1f3a5f'
GREY = '#5b6b7f'
A4 = (8.27, 11.69)

D5 = 'p2_det1_long_run_efficiency_7-19-26/long_run_det1_415_615'   # det1_long5
MS = 'p2_det1_long_run_mesh_scan_7-19-26/mesh_scan'                # det1_meshscan1
DS = 'p2_det1_drift_scan_7-19-26/drift_scan'                       # det1_driftscan2
H2 = 'p2_det1_mesh_hv_scan_7-2-26/hv_scan'                         # det1_hvscan (P1)
L1 = 'p2_det1_long_run_6-30-26/efficiency_long_run_p2_det1'        # det1_long  (P1)
L3 = 'p2_det1_long_run_7-7-26/mesh_440V_drift_600V'                # det1_long3 (P1)
L4 = 'p2_det1_det2_long_run_mesh_scan_7-9-26/long_run_mesh_430V_drift_600V'

V2 = '_without_connectors_1_2_10_spark_vetoed.png'   # z=702 runs (conn 2 dropped)
V1 = '_without_connectors_1_10_spark_vetoed.png'     # P1 runs
V0 = '_without_connector_10_spark_vetoed.png'        # det1_long (only conn 10)

# (section, path, caption)  -- caption states the ARGUMENT the figure supports
FIGURES = [
    ('1. The measurement — det1 at the P2 plane (det1_long5, 24 h @ 415/615)',
     f'{D5}/06_efficiency/efficiency_map_sliding{V2}',
     'Local efficiency is flat at ~1.0 across the whole wedge. The dark spots '
     'are the five big pillars and a few dead pads; the third panel is M3 ray '
     'density, which is what produces the occupancy gradient seen in raw '
     'hitmaps — it divides out here. Evidence for §1 and §3.1.'),
    ('1. The measurement — det1 at the P2 plane (det1_long5, 24 h @ 415/615)',
     f'{D5}/06_efficiency/efficiency_breakdown{V2}',
     'The loss budget of §3: reco 93.4 % | fired-not-reco 3.4 % | no-hit 3.2 %.'),
    ('1. The measurement — det1 at the P2 plane (det1_long5, 24 h @ 415/615)',
     f'{D5}/06_efficiency/efficiency_vs_time{V2}',
     'Stability over the 24 h run — the quoted figure is not a duty-cycle '
     'average of an intermittent detector (contrast §4).'),
    ('1. The measurement — det1 at the P2 plane (det1_long5, 24 h @ 415/615)',
     f'{D5}/06_efficiency/radial_residual{V2}',
     'Residual distribution, core σ = 4.7 mm. This sets the match radius: '
     'R = 40 mm is ~8σ, comfortably past the ~3σ needed for containment.'),

    ('2. The cuts sit on plateaus (not tuned to a number)',
     f'{D5}/12_validation/validation_knob_plateau{V2}',
     'The central methodological plot. LEFT: eff vs match radius R tracks '
     'ideal Gaussian containment and is flat from ~20 mm, so R = 40 mm is on '
     'the plateau. RIGHT: eff vs active-area radius is flat over 5–9 mm and '
     'then falls; the shaded band ends at the pad half-diagonal (8.2 mm) and '
     'the grey curve (rays in the denominator) saturates at exactly the same '
     'place. Two independent anchors for active_r = 8.4 mm.'),
    ('2. The cuts sit on plateaus (not tuned to a number)',
     f'{D5}/03_m3_alignment/best_angle_corr{V2}',
     'P2↔M3 correlation r = 0.946 / 0.949 at θ = 89.0°. Before the hot-pad '
     'mask this was r = 0.651 / 0.574 at θ = 75.5° — six oscillating pads '
     'carrying 69 % of the hits had pinned the centroid and rotated the '
     'geometry fit by 13.5°.'),
    ('2. The cuts sit on plateaus (not tuned to a number)',
     f'{D5}/12_validation/validation_mapping_vs_dead{V2}',
     'Low-efficiency pads fire at 0 % of the good-pad rate, i.e. they are '
     'genuinely dead rather than mis-mapped. After removing connector 2 there '
     'are 0 pads with local eff < 15 %.'),
    ('2. The cuts sit on plateaus (not tuned to a number)',
     f'{D5}/02_map_validation/pad_occupancy{V2}',
     'Per-connector occupancy is uniform (60–90 k) once the hot pads are '
     'masked and connector 2 is dropped. Before, connector 9 alone held '
     '970 k hits.'),

    ('3. Where the missing 6.6 % goes',
     f'{D5}/05_detector_deep_qa/surface_hitmap{V2}',
     'Pad-level hitmap. The masked hot pads and the five big pillars are '
     'marked; the grey band along the bottom edge is connector 2 (dead). '
     'Supports §3.1.'),
    ('3. Where the missing 6.6 % goes',
     f'{D5}/06_efficiency/map_has_any{V2}',
     'Detection map (any pad fired, no position requirement). The residual '
     'holes are the pillars — this is the 96.8 % ceiling of §3.1.'),
    ('3. Where the missing 6.6 % goes',
     f'{D5}/06_efficiency/nonreco_ray_positions{V2}',
     'Where reconstruction fails. Not concentrated in any region, consistent '
     'with §3.2: the centroid is pulled off by a second cluster rather than '
     'the detector failing somewhere specific.'),

    ('4. Mesh (gain) scan — 7-19, P2 plane',
     f'{MS}/11_hv_scan_efficiency/efficiency_vs_hv{V2}',
     'Clean turn-on with a plateau from ~400 V at 93.7–94.6 % reco / '
     '97–98 % any-pad. The 415 V operating point of the long run sits on the '
     'plateau, so the quoted efficiency is not a point on a slope.'),
    ('4. Mesh (gain) scan — 7-19, P2 plane',
     f'{MS}/11_hv_scan_efficiency/amplitude_vs_hv{V2}',
     'Pad amplitude rising with mesh voltage — the gain curve behind the '
     'efficiency turn-on.'),
    ('4. Mesh (gain) scan — 7-19, P2 plane',
     f'{MS}/11_hv_scan_efficiency/resolution_vs_hv{V2}',
     'Position resolution across the scan; flat on the plateau, so the '
     'operating point is not trading resolution for efficiency.'),

    ('5. Drift scan — 7-19, P2 plane, mesh fixed at 415 V',
     f'{DS}/16_drift_scan_efficiency/efficiency_vs_drift{V2}',
     'Flat at ~94 % from 565 V. The 415 V point collapsing to 18 % is a null '
     'test, not a detector effect: drift = mesh = 415 V means zero field '
     'across the drift gap, so primaries never reach the mesh.'),
    ('5. Drift scan — 7-19, P2 plane, mesh fixed at 415 V',
     f'{DS}/16_drift_scan_efficiency/drift_velocity_vs_drift{V2}',
     'CAUTION — the "apparent" curve is NOT an independent measurement of the '
     'drift velocity. It is the measured p90−p10 peak-time spread times an '
     'assumed geometry factor of 0.8 (the p90−p10 of a uniform distribution, '
     'i.e. assuming the detected arrival times fill the gap uniformly). With '
     'the gap known at 3.0 mm and the Magboltz table verified to 0.3 % against '
     'an independent run, the data requires a factor of 0.63 ± 0.03 instead — '
     'the arrival-time distribution is more concentrated than uniform, because '
     'the ≥300 ADC signal-band cut removes the diffused late tail. The offset '
     'is an estimator artefact, not a gas effect (§6). Re-fitting the factor '
     'would make this agree with Magboltz by construction and must not be '
     'presented as a measurement.'),
    ('5. Drift scan — 7-19, P2 plane, mesh fixed at 415 V',
     f'{DS}/16_drift_scan_efficiency/amplitude_vs_drift{V2}',
     'Amplitude vs drift field — confirms the plateau is not an amplitude '
     'effect.'),
    ('5. Drift scan — 7-19, P2 plane, mesh fixed at 415 V',
     f'{DS}/16_drift_scan_efficiency/time_resolution_vs_drift{V2}',
     'Time resolution across the drift scan. Uses decoded_root, which was '
     'verified unchanged by the reprocessing (§5), so this measurement did '
     'not need repeating.'),

    ('6. Surface timing map (stage 19) — instrumental, not drift',
     f'{D5}/19_timing_surface/timing_surface{V2}',
     'LEFT: per-pad median peak time. In a parallel-plate geometry the drift '
     'time depends on the DEPTH z, not on (x,y), so the physics prediction for '
     'this map is flat — a Garfield/HEED drift simulation would predict no '
     'structure and cannot interpret one. The coherent diagonal banding is '
     'therefore instrumental (cable / front-end delays, following the channel '
     'ordering). It is REAL, not noise: splitting the events into even/odd '
     'halves gives two independent estimates per pad that correlate at '
     'r = 0.88, with a genuine pad-to-pad spread of 7.5 ns against only 2.8 ns '
     'of statistical error. RIGHT: per-pad time resolution, median 20.8 ns — '
     'the same sigma = (p84.1-p15.9)/2 definition stage 16 uses for its ~21 ns '
     'event-level figure, so the two are directly comparable. The same '
     'banding appears here, i.e. the affected channels are both shifted and '
     'slightly worse resolved.'),
    ('6. Surface timing map (stage 19) — instrumental, not drift',
     f'{D5}/19_timing_surface/timing_per_connector{V2}',
     'The same offsets grouped by connector: an ~11 ns range between '
     'connectors (5 at +3.3 ns, 7 at -7.6 ns), but with a within-connector '
     'spread of 5-9 ns, so the connector explains only part of it — the rest '
     'is finer, at the pad-row level. Total per-pad sigma = 7.7 ns over 885 '
     'pads, comparable to the 15-22 ns per-station timing resolution from the '
     'beam test, so a per-pad timing calibration from this map is worth '
     'having.'),

    ('6. Surface timing map (stage 19) — instrumental, not drift',
     f'{D5}/19_timing_surface/timing_map_sliding{V2}',
     'Continuous version of the same measurement: a sliding kernel (r = 6 mm) '
     'in the M3 reference frame, so it overlays the stage-10 efficiency map '
     'directly. LEFT: the offset structure as coherent lobes rather than pad '
     'speckle — a blue region near X=-25, Y=100 at about -8 ns against red '
     'through the lower-right at about +8 ns. MIDDLE: resolution, 18-22 ns '
     'across the plane. RIGHT: events per kernel, blank where there are none; '
     'the noisier sigma patches at the edges track the low-occupancy rim, so '
     'they are statistics rather than a real feature. Kernel is 6 mm and not '
     'smaller on purpose: each event carries the time of ONE pad, so nothing '
     'finer than the 11.8 mm pad pitch is resolved and the M3 track adds '
     '~5 mm of its own.'),

    ('7. The P1 / P2 question (open — §4)',
     f'{H2}/11_hv_scan_efficiency/efficiency_vs_hv{V1}',
     'THE key figure for §4. On 7-2, at the P1 plane, det1 reaches 96.1 % '
     'any-pad — matching its best P2-plane performance. P1 is therefore not '
     'inherently a low-efficiency position, and the P1 behaviour is '
     'intermittency, not a monotonic decline.'),
    ('7. The P1 / P2 question (open — §4)',
     f'{L1}/06_efficiency/efficiency_map_sliding{V0}',
     'det1_long (6-30, P1): 54.0 % reco / 64.5 % any. Compare the uniform '
     'P2-plane map in §1 — the loss here is spread over the detector, not '
     'localised.'),
    ('7. The P1 / P2 question (open — §4)',
     f'{L1}/06_efficiency/efficiency_vs_time{V0}',
     'Efficiency vs time within the 6-30 run — the within-run behaviour '
     'behind the intermittency argument.'),
    ('7. The P1 / P2 question (open — §4)',
     f'{L3}/06_efficiency/efficiency_map_sliding{V1}',
     'det1_long3 (7-7, P1): 16.1 % reco / 25.5 % any.'),
    ('7. The P1 / P2 question (open — §4)',
     f'{L4}/06_efficiency/efficiency_map_sliding{V1}',
     'det1_long4 (7-9, P1): 3.6 % reco / 10.8 % any — the same sub_run in '
     'which det2 read 77 %, which is the internal control showing the bench '
     'and the M3 telescope were fine while det1 was not.'),
]


def draw_line(fig, y, txt, size=9, weight='normal', color=INK, mono=False, x=0.07):
    fig.text(x, y, txt, fontsize=size, weight=weight, color=color,
             family='monospace' if mono else 'sans-serif', va='top')


def clean(t):
    t = re.sub(r'\[([^\]]+)\]\([^)]+\)', r'\1', t)
    return t.replace('**', '').replace('`', '').replace('*', '')


def render_markdown(pdf, md_path):
    lines = open(md_path).read().splitlines()
    fig = plt.figure(figsize=A4); y = 0.945
    in_code = False
    table = []

    def flush(fig, y):
        if not table:
            return fig, y
        rows = [[c.strip() for c in r.strip().strip('|').split('|')]
                for r in table if not re.match(r'^\|[\s:|-]+\|$', r)]
        if rows:
            nc = max(len(r) for r in rows)
            w = [max(len(r[i]) if i < len(r) else 0 for r in rows) for i in range(nc)]
            for ri, r in enumerate(rows):
                if y < 0.06:
                    pdf.savefig(fig); plt.close(fig)
                    fig = plt.figure(figsize=A4); y = 0.945
                cells = [(r[i] if i < len(r) else '').ljust(w[i]) for i in range(nc)]
                draw_line(fig, y, '  '.join(cells)[:120], size=6.8, mono=True,
                          weight='bold' if ri == 0 else 'normal',
                          color=INK if ri == 0 else GREY)
                y -= 0.0135
        table.clear()
        return fig, y - 0.008

    for raw in lines:
        s = raw.rstrip()
        if s.startswith('```'):
            in_code = not in_code
            continue
        if s.lstrip().startswith('|'):
            table.append(s)
            continue
        fig, y = flush(fig, y)
        if y < 0.07:
            pdf.savefig(fig); plt.close(fig)
            fig = plt.figure(figsize=A4); y = 0.945
        if not s.strip():
            y -= 0.010
            continue
        if in_code:
            draw_line(fig, y, s[:110], size=7.5, mono=True, color=GREY, x=0.09)
            y -= 0.014
            continue
        m = re.match(r'^(#{1,4})\s+(.*)', s)
        if m:
            lvl = len(m.group(1))
            if lvl <= 2 and y < 0.88:
                y -= 0.012
            draw_line(fig, y, clean(m.group(2))[:95],
                      size={1: 16, 2: 12.5, 3: 10.5, 4: 9.5}[lvl], weight='bold')
            y -= 0.020 + 0.006 * (4 - lvl)
            continue
        b = re.match(r'^(\s*)[-*]\s+(.*)', s)
        body = clean(b.group(2)) if b else clean(s)
        for i, chunk in enumerate(textwrap.wrap(body, 104) or ['']):
            if y < 0.06:
                pdf.savefig(fig); plt.close(fig)
                fig = plt.figure(figsize=A4); y = 0.945
            draw_line(fig, y, ('  • ' if (b and i == 0) else ('    ' if b else ''))
                      + chunk, size=8.4, color=GREY if b else INK)
            y -= 0.0155
        y -= 0.004
    fig, y = flush(fig, y)
    pdf.savefig(fig); plt.close(fig)


def figure_pages(pdf, base):
    missing = []
    last_section = None
    for section, rel, caption in FIGURES:
        path = os.path.join(base, rel)
        if not os.path.isfile(path):
            missing.append(rel)
            continue
        fig = plt.figure(figsize=A4)
        if section != last_section:
            fig.text(0.07, 0.975, section, fontsize=10.5, weight='bold',
                     color=INK, va='top')
            last_section = section
        else:
            fig.text(0.07, 0.975, section, fontsize=8, color=GREY, va='top')
        img = mpimg.imread(path)
        # Size the frame to the image's own aspect so wide plots are drawn as
        # large as the page allows instead of floating in a fixed tall box.
        ih, iw = img.shape[0], img.shape[1]
        w_fig = 0.90
        h_fig = w_fig * (ih / iw) * (A4[0] / A4[1])
        top = 0.945
        h_fig = min(h_fig, 0.62)
        ax = fig.add_axes([0.05, top - h_fig, w_fig, h_fig])
        ax.imshow(img, interpolation='antialiased')
        ax.axis('off')
        cap_top = top - h_fig - 0.030
        fig.text(0.07, cap_top, os.path.basename(rel).replace('.png', ''),
                 fontsize=6.5, color=GREY, va='top', family='monospace')
        y = cap_top - 0.022
        for chunk in textwrap.wrap(caption, 96):
            fig.text(0.07, y, chunk, fontsize=8.6, color=INK, va='top')
            y -= 0.017
        fig.text(0.07, 0.04, rel, fontsize=5.5, color=GREY, family='monospace')
        pdf.savefig(fig); plt.close(fig)
    return missing


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-o', '--out', default=None)
    ap.add_argument('--md', default=None)
    args = ap.parse_args()

    root = qa.ANALYSIS_ROOT
    base = os.path.join(root, 'det1')
    md = args.md or os.path.join(root, 'DET1_EFFICIENCY.md')
    out = args.out or os.path.join(base, 'DET1_EFFICIENCY_REPORT.pdf')

    with PdfPages(out) as pdf:
        fig = plt.figure(figsize=A4)
        fig.text(0.5, 0.66, 'det1 (P2_1)', ha='center', fontsize=26,
                 weight='bold', color=INK)
        fig.text(0.5, 0.60, 'Efficiency measurements and loss budget',
                 ha='center', fontsize=13, color=INK)
        fig.text(0.5, 0.555, 'long runs · mesh HV scans · drift scans',
                 ha='center', fontsize=10, color=GREY)
        fig.text(0.5, 0.51, 'August-2026 reprocessed data', ha='center',
                 fontsize=9.5, color=GREY)
        fig.text(0.5, 0.46, datetime.datetime.now().strftime('generated %Y-%m-%d %H:%M'),
                 ha='center', fontsize=9, color=GREY)
        pdf.savefig(fig); plt.close(fig)

        if os.path.isfile(md):
            render_markdown(pdf, md)
        else:
            print('!! markdown missing:', md)

        missing = figure_pages(pdf, base)

    n = len(FIGURES) - len(missing)
    print(f'wrote {out} ({os.path.getsize(out)/1e6:.1f} MB), {n}/{len(FIGURES)} figures')
    for m in missing:
        print('  !! MISSING figure:', m)


if __name__ == '__main__':
    main()
