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

# Standard run per detector (agreed 2026-08-12): one long run each, chosen for
# being healthy and having the best statistics. The others stay in the product
# tree as supporting material.
CROSS = {
    'charging': dict(
        md='CHARGING_UP.md',
        title='Charging-up and run-time stability — all four detectors',
        figs=[('Per-detector timelines', 'det1', 'p2_det1_long_run_efficiency_7-19-26/'
               'long_run_det1_415_615/20_charging_up/charging_timeline'
               '_without_connectors_1_2_10_spark_vetoed.png',
               'det1: +3.4 points over 18 h, plateaued. Run average within 0.4 '
               'points of the late window.'),
              ('Per-detector timelines', 'det2', 'p2_det1_det2_long_run_mesh_scan_7-9-26/'
               'long_run_mesh_430V_drift_600V/20_charging_up/charging_timeline'
               '_without_connectors_1_8_9_10_spark_vetoed.png',
               'det2: no charging-up at all, flat within ~1 point across 10 h.'),
              ('Per-detector timelines', 'det3', 'p2_det3_mesh_scan_det4_initial_7-16-26/'
               'initial_run_det3_420_820_det4_430_830/20_charging_up/charging_timeline'
               '_without_connectors_1_8_9_10_spark_vetoed.png',
               'det3: efficiency FALLS 90.6 -> 81.0 % inside 2 h. Not charging-up '
               '— the collapse had already begun. No window is det3s efficiency.'),
              ('Per-detector timelines', 'det4', 'p2_det4_long_run_7-20-26/'
               'long_run_det4_410_610/20_charging_up/charging_timeline'
               '_without_connectors_1_2_7_10_spark_vetoed.png',
               'det4: +17.6 points, rate up a factor 9, flat after ~4 h. The run '
               'average understates the detector by 4.0 points.'),
              ('Combined', None, 'combined/combined_efficiency.png',
               'All four detectors, run average vs charged-up late window. det3 '
               'hatched: it was degrading, so neither bar is its efficiency.'),
              ('Combined', None, 'combined/combined_mesh_scan.png',
               'Mesh (gain) scans. det3 excluded — its scan ran while the '
               'detector was collapsing, so it records the collapse, not gain.'),
              ('Combined', None, 'combined/combined_drift_scan.png',
               'Drift scans. LEFT: efficiency; zero field is a null test and must '
               'collapse — det4 at the P1 plane stays FLAT, i.e. no drift field '
               'was established. RIGHT: time resolution; hollow points carry too '
               'few events to be a resolution.')]),
}

STANDARD = {
    'det1': dict(
        long='p2_det1_long_run_efficiency_7-19-26/long_run_det1_415_615',
        mesh='p2_det1_long_run_mesh_scan_7-19-26/mesh_scan',
        drift='p2_det1_drift_scan_7-19-26/drift_scan',
        sfx='_without_connectors_1_2_10_spark_vetoed.png',
        title='det1 (P2_1) — 24 h @ 415/615, P2 plane',
        # det1 is the only detector with runs at BOTH bench planes, so it keeps
        # a P1/P2 section (document section 4). Suffixes differ per run.
        extra=[('7. The P1 / P2 question (open — see section 4)',
                'p2_det1_mesh_hv_scan_7-2-26/hv_scan/11_hv_scan_efficiency/'
                'efficiency_vs_hv_without_connectors_1_10_spark_vetoed.png',
                'On 7-2, at the P1 plane, det1 reached 96.1 % any-pad — as good '
                'as its best P2-plane performance. So P1 is not inherently a '
                'low-efficiency position and the P1 behaviour is intermittency, '
                'not a monotonic decline.'),
               ('7. The P1 / P2 question (open — see section 4)',
                'p2_det1_long_run_6-30-26/efficiency_long_run_p2_det1/'
                '06_efficiency/efficiency_map_sliding'
                '_without_connector_10_spark_vetoed.png',
                'det1_long (6-30, P1): 54.0 % reco / 64.5 % any. The loss is '
                'spread over the detector rather than localised.'),
               ('7. The P1 / P2 question (open — see section 4)',
                'p2_det1_long_run_7-7-26/mesh_440V_drift_600V/06_efficiency/'
                'efficiency_map_sliding_without_connectors_1_10_spark_vetoed.png',
                'det1_long3 (7-7, P1): 16.1 % reco / 25.5 % any.'),
               ('7. The P1 / P2 question (open — see section 4)',
                'p2_det1_det2_long_run_mesh_scan_7-9-26/'
                'long_run_mesh_430V_drift_600V/06_efficiency/'
                'efficiency_map_sliding_without_connectors_1_10_spark_vetoed.png',
                'det1_long4 (7-9, P1): 3.6 % reco / 10.8 % any — the same '
                'sub_run in which det2 read 90.3 %, the internal control '
                'showing the bench and M3 telescope were fine.')]),
    'det2': dict(
        long='p2_det1_det2_long_run_mesh_scan_7-9-26/long_run_mesh_430V_drift_600V',
        mesh='p2_det1_det2_long_run_mesh_scan_7-9-26/hv_scan',
        drift=None,                      # det2 has no drift scan
        sfx='_without_connectors_1_8_9_10_spark_vetoed.png',
        title='det2 (P2_2) — 10 h @ 430/600, P2 plane'),
    'det3': dict(
        long='p2_det3_mesh_scan_det4_initial_7-16-26/initial_run_det3_420_820_det4_430_830',
        mesh='p2_det3_mesh_scan_det4_initial_7-16-26/mesh_scan',
        drift='p2_det3_det4_drift_scan_7-16-26/drift_scan',
        sfx='_without_connectors_1_8_9_10_spark_vetoed.png',
        title='det3 (P2_3) — 2 h @ 420/820, P2 plane'),
    'det4': dict(
        long='p2_det4_long_run_7-20-26/long_run_det4_410_610',
        mesh=None,                       # det4 mesh scan was never taken
        drift='p2_det3_det4_drift_scan_7-16-26/drift_scan',
        sfx='_without_connectors_1_2_7_10_spark_vetoed.png',
        # The 7-16 drift scan is det4 at the P1 plane with a different cabling
        # (dead connectors 1,10), so its products carry a different suffix from
        # the 7-20 long run. Suffixes follow the RUN, not the detector.
        drift_sfx='_without_connectors_1_10_spark_vetoed.png',
        title='det4 (P2_4) — 10.4 h @ 410/610, P2 plane'),
}


def figure_spec(det):
    """(section, relative path, caption) for one detector's standard run."""
    c = STANDARD[det]
    L, M, D, V = c['long'], c['mesh'], c['drift'], c['sfx']
    VM = c.get('mesh_sfx', V)     # scan products can carry another suffix
    VD = c.get('drift_sfx', V)
    F = []
    F += [('1. The measurement — standard long run',
           f'{L}/06_efficiency/efficiency_map_sliding{V}',
           'Local efficiency across the plane. Dark spots are the five big '
           'mesh-support pillars and any dead pads; the third panel is M3 ray '
           'density, which produces the gradient seen in raw hitmaps and '
           'divides out here.'),
          ('1. The measurement — standard long run',
           f'{L}/06_efficiency/efficiency_breakdown{V}',
           'The loss budget: reco / fired-not-reco / no-hit.'),
          ('1. The measurement — standard long run',
           f'{L}/06_efficiency/efficiency_vs_time{V}',
           'Stability over the run — shows whether the quoted number is a '
           'steady value or a duty-cycle average of an intermittent detector.'),
          ('1. The measurement — standard long run',
           f'{L}/06_efficiency/radial_residual{V}',
           'Residual distribution; its core sigma sets the match radius R.'),
          ('2. Cuts on plateaus',
           f'{L}/12_validation/validation_knob_plateau{V}',
           'LEFT: efficiency vs match radius R against ideal Gaussian '
           'containment. RIGHT: efficiency vs active-area radius, with the pad '
           'half-diagonal marked and the denominator-occupancy curve; both '
           'anchor active_r = 8.4 mm.'),
          ('2. Cuts on plateaus',
           f'{L}/03_m3_alignment/best_angle_corr{V}',
           'P2 to M3 position correlation at the fitted rotation.'),
          ('2. Cuts on plateaus',
           f'{L}/12_validation/validation_mapping_vs_dead{V}',
           'Separates genuinely dead pads from mis-mapped ones: dead pads must '
           'also stop firing, not just reconstruct elsewhere.'),
          ('2. Cuts on plateaus',
           f'{L}/02_map_validation/pad_occupancy{V}',
           'Per-connector occupancy uniformity after hot-pad masking and dead '
           'connector removal.'),
          ('3. Where the inefficiency is',
           f'{L}/05_detector_deep_qa/surface_hitmap{V}',
           'Pad-level hitmap with masked hot pads and the big pillars marked.'),
          ('3. Where the inefficiency is',
           f'{L}/06_efficiency/map_has_any{V}',
           'Detection map (any pad fired, no position requirement) — the '
           'ceiling that position reconstruction works against.'),
          ('3. Where the inefficiency is',
           f'{L}/06_efficiency/nonreco_ray_positions{V}',
           'Where reconstruction fails; a diffuse pattern points at centroid '
           'confusion rather than a localised detector defect.')]
    if M:
        F += [('4. Mesh (gain) scan',
               f'{M}/11_hv_scan_efficiency/efficiency_vs_hv{VM}',
               'Efficiency vs mesh voltage. Long/initial/final runs sharing a '
               'working point are excluded from the curve.'),
              ('4. Mesh (gain) scan',
               f'{M}/11_hv_scan_efficiency/amplitude_vs_hv{VM}',
               'Pad amplitude vs mesh voltage. Read amp_median, not amp_mean: '
               'where efficiency is low the surviving signal-band hits are '
               'noise tails and the mean is misleading.'),
              ('4. Mesh (gain) scan',
               f'{M}/11_hv_scan_efficiency/resolution_vs_hv{VM}',
               'Position resolution across the scan.')]
    if D:
        F += [('5. Drift scan',
               f'{D}/16_drift_scan_efficiency/efficiency_vs_drift{VD}',
               'Efficiency vs drift field. The zero-field point (drift = mesh) '
               'is a null test: efficiency must collapse there. A FLAT curve '
               'means no drift field is being established.'),
              ('5. Drift scan',
               f'{D}/16_drift_scan_efficiency/time_resolution_vs_drift{VD}',
               'Time resolution vs drift field; should fall and plateau.'),
              ('5. Drift scan',
               f'{D}/16_drift_scan_efficiency/drift_velocity_vs_drift{VD}',
               'CAUTION: the "apparent" curve is not an independent velocity '
               'measurement — it is the p90-p10 peak-time spread times an '
               'assumed 0.8 geometry factor, which the data puts at 0.63.'),
              ('5. Drift scan',
               f'{D}/16_drift_scan_efficiency/amplitude_vs_drift{VD}',
               'Amplitude vs drift field.')]
    F += [('6. Surface timing map',
           f'{L}/19_timing_surface/timing_surface{V}',
           'LEFT: per-pad median peak time. Drift time depends on depth z, not '
           'on (x,y), so the physics prediction is flat and any structure is '
           'instrumental (cable / front-end delays). The split-half test in the '
           'title separates genuine pad-to-pad spread from statistics. RIGHT: '
           'per-pad resolution, sigma = (p84.1-p15.9)/2, the same definition '
           'stage 16 uses.'),
          ('6. Surface timing map',
           f'{L}/19_timing_surface/timing_per_connector{V}',
           'The same offsets grouped by connector, plus their distribution.'),
          ('6. Surface timing map',
           f'{L}/19_timing_surface/timing_map_sliding{V}',
           'Continuous version in the M3 reference frame, overlaying the '
           'stage-10 efficiency map. Kernel auto-widens on sparse runs and is '
           'stated on the plot; nothing finer than the 11.8 mm pad pitch is '
           'resolved.')]
    F += [('7. Charging-up and run-time stability',
           f'{L}/20_charging_up/charging_timeline{V}',
           'Efficiency, rate, arrival time and amplitude vs time. Early and '
           'late windows shaded. A monotonic rise that flattens is charging-up '
           '(quote the late window); a FALL means the detector was degrading '
           'and no single number is its efficiency; an ON/OFF rate with a '
           'matching arrival-time shift would instead indicate an intermittent '
           'drift field. See CHARGING_UP.md.')]
    F += c.get('extra', [])
    return F


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


def figure_pages(pdf, base, FIGURES):
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
    ap.add_argument('--det', default='det1',
                    choices=sorted(STANDARD) + sorted(CROSS))
    ap.add_argument('-o', '--out', default=None)
    ap.add_argument('--md', default=None)
    args = ap.parse_args()

    det = args.det
    root = qa.ANALYSIS_ROOT
    if det in CROSS:
        c = CROSS[det]
        base = root
        md = args.md or os.path.join(root, c['md'])
        out = args.out or os.path.join(root, 'combined',
                                       c['md'].replace('.md', '_REPORT.pdf'))
        os.makedirs(os.path.dirname(out), exist_ok=True)
        FIGURES = [(sec, (os.path.join(d, r) if d else r), cap)
                   for sec, d, r, cap in c['figs']]
        TITLE = c['title']
    else:
        base = os.path.join(root, det)
        md = args.md or os.path.join(root, f'{det.upper()}_EFFICIENCY.md')
        out = args.out or os.path.join(base, f'{det.upper()}_EFFICIENCY_REPORT.pdf')
        FIGURES = figure_spec(det)
        TITLE = STANDARD[det]['title']

    with PdfPages(out) as pdf:
        fig = plt.figure(figsize=A4)
        fig.text(0.5, 0.66, TITLE.split(' — ')[0], ha='center',
                 fontsize=22, weight='bold', color=INK)
        fig.text(0.5, 0.625, TITLE.split(' — ')[-1],
                 ha='center', fontsize=11, color=GREY)
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

        missing = figure_pages(pdf, base, FIGURES)

    n = len(FIGURES) - len(missing)
    print(f'wrote {out} ({os.path.getsize(out)/1e6:.1f} MB), {n}/{len(FIGURES)} figures')
    for m in missing:
        print('  !! MISSING figure:', m)


if __name__ == '__main__':
    main()
