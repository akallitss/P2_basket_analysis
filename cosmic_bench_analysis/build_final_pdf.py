#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_final_pdf.py

Compile the full P2 (BASKET) det1 cosmic-bench QA into a single multi-page PDF
that gathers ALL the investigations and observations from the long-run analysis
(stages 02..09) into one styled document. Adapted in spirit from
nTof_x17/mx_june_cosmic_qa/build_final_pdf.py, but P2 is a metallic pad
Micromegas read out on pads with a MESH HV (not resistive strips), so this is a
per-stage narrative report rather than one dense per-detector page.

Layout (A4 portrait):
  1. Cover  : operating point + headline stat cards + key observations
  2. Map    : pad occupancy/amplitude, dead-pad map, ordering validation (02)
  3. Align  : M3 reference + P2<->M3 rotation/z-theta alignment (04 + 03)
  4. Effic. : efficiency maps, breakdown, radial residual (06)
  5. Sliding: sliding-window efficiency map + non-reco rays -> diagonal dead band
  6. DeepQA : surface hitmap, pad firing, multiplicity (05)
  7. HV spk : mesh spark timeline / rate / DAQ cross-check (07)
  8. PadSpk : per-pad micro-spark flagging, review-only (08)
  9. Pedest.: pedestal mean/rms maps + RMS-vs-rate (09)

Everything is read from the det1_long OUT_BASE under .../Analysis/det1/... . Each
image is looked up suffix-tolerantly (spark-vetoed variant preferred).

Usage:
  python3 build_final_pdf.py [run_key] [--out=PATH]     (default key det1_long)
"""
import glob
import os
import re
import sys
import datetime

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec

import p2_qa_config as qa

# Suffix search order: prefer the spark-vetoed dead-connector product, then the
# un-vetoed dead-connector product, then a bare name. This is a module default
# for a single-dead-connector run; build() overwrites it from the actual run
# config so a run dropping e.g. connectors 1 & 10 resolves the right filenames.
_SUFFIXES = ('_without_connector_10_spark_vetoed', '_without_connector_10', '')


def suffixes_for(cfg):
    """Filename-suffix search order for a run: vetoed dead-connector product
    first, then the un-vetoed one, then a bare name. Derived from the run's own
    dead-connector set so it matches whatever cfg.product_suffix() wrote."""
    ordered = (cfg.product_suffix(veto_sparks=True), cfg.dead_suffix, '')
    return tuple(dict.fromkeys(ordered))   # de-dupe, preserve order

INK = '#1f3a5f'


def find_img(base, stem):
    """First existing PNG among the suffix variants of a stage-relative stem.
    The stem may contain a glob wildcard (e.g. 'map_within_*mm' — the match
    radius R is a per-run config knob, so the filename varies)."""
    for suf in _SUFFIXES:
        pat = os.path.join(base, stem + suf + '.png')
        if '*' in pat or '?' in pat:
            hits = sorted(glob.glob(pat))
            if hits:
                return hits[0]
        elif os.path.isfile(pat):
            return pat
    return None


def read_txt(base, *cands):
    for c in cands:
        p = os.path.join(base, c)
        if os.path.isfile(p):
            return open(p).read()
    return ''


def grab(txt, pat, default='—'):
    m = re.search(pat, txt)
    return m.group(1) if m else default


def headline(base, suffixes=_SUFFIXES):
    """Parse the stage summary text files into the cover stat cards. Summary
    filenames are resolved against `suffixes` (vetoed / dead-connector / bare)
    so the cover works for any dead-connector set, not just connector 10."""
    def cands(stem):
        return [stem + suf + '.txt' for suf in suffixes]
    eff = read_txt(base, *cands('06_efficiency/efficiency_summary'))
    deep = read_txt(base, *cands('05_detector_deep_qa/deep_qa_summary'))
    spk = read_txt(base, *cands('07_hv_spark_qa/spark_qa_summary'))
    ped = read_txt(base, *cands('09_pedestal_qa/pedestal_qa_summary'))
    pad = read_txt(base, *cands('08_pad_spark_qa/pad_spark_qa_summary'))
    return dict(
        within=grab(eff, r'within \d+(?:\.\d+)? mm\s*:\s*([\d.]+)%'),
        has_any=grab(eff, r'has_any \(fired\)\s*:\s*([\d.]+)%'),
        med_r=grab(eff, r'median \|r\| residual\s*:\s*([\d.]+)'),
        rays=grab(eff, r'clean M3 rays \(total\)\s*:\s*([\d,]+)'),
        rays_act=grab(eff, r'rays in active area\s*:\s*([\d,]+)'),
        rot=grab(eff, r'transform rotation/scale\s*:\s*([\d.]+) deg'),
        scale=grab(eff, r'/\s*([\d.]+) \(reflection'),
        no_hit=grab(eff, r'no-hit ([\d.]+)%'),
        fnr=grab(eff, r'fired-not-reco ([\d.]+)%'),
        pads=grab(deep, r'distinct pads fired:\s*(\d+)'),
        hpe=grab(deep, r'\(([\d.]+) hits/event\)'),
        sparks=grab(spk, r'sparks detected\s*:\s*(\d+)'),
        deadt=grab(spk, r'\(([\d.]+)% of run\)'),
        ped_rms=grab(ped, r'median ([\d.]+),'),
        spark_pads=grab(pad, r'SPARK pads[^:]*:\s*(\d+)'),
    )


# --------------------------------------------------------------------------- #
# generic page renderers
# --------------------------------------------------------------------------- #
def place(ax, img_path, label):
    ax.axis('off')
    if img_path and os.path.isfile(img_path):
        ax.imshow(mpimg.imread(img_path), interpolation='antialiased', resample=True)
    else:
        ax.add_patch(plt.Rectangle((0, 0), 1, 1, fc='whitesmoke', ec='grey', ls='--'))
        ax.text(0.5, 0.5, f'(missing)\n{label}', ha='center', va='center',
                fontsize=8, color='grey', transform=ax.transAxes)
    ax.set_title(label, fontsize=8.5, pad=3)


def notes_box(fig, text, y=0.045):
    fig.text(0.5, y, text, ha='center', va='center', fontsize=8.2, color='#222222',
             wrap=True,
             bbox=dict(boxstyle='round,pad=0.6', fc='#eef3f8', ec='#9db8d2', lw=0.8))


def grid_page(pdf, title, subtitle, tiles, notes=None, ncols=2):
    """A titled page: ncols grid of (img_path, caption) tiles + optional notes."""
    nrows = (len(tiles) + ncols - 1) // ncols
    fig = plt.figure(figsize=(8.27, 11.69))
    fig.text(0.06, 0.972, title, fontsize=19, fontweight='bold', color=INK, va='top')
    if subtitle:
        fig.text(0.06, 0.945, subtitle, fontsize=9.5, color='#555555', va='top')
    bottom = 0.17 if notes else 0.05
    gs = GridSpec(nrows, ncols, figure=fig, hspace=0.22, wspace=0.08,
                  left=0.05, right=0.95, top=0.915, bottom=bottom)
    for i, (img, cap) in enumerate(tiles):
        r, c = i // ncols, i % ncols
        if i == len(tiles) - 1 and c == 0 and len(tiles) % ncols:
            place(fig.add_subplot(gs[r, :]), img, cap)  # lone last tile -> full width
        else:
            place(fig.add_subplot(gs[r, c]), img, cap)
    if notes:
        notes_box(fig, notes)
    fig.text(0.96, 0.008, datetime.date.today().isoformat(), ha='right',
             fontsize=6, color='grey')
    pdf.savefig(fig, dpi=200)
    plt.close(fig)


def cover_page(pdf, cfg, h, rv_mesh, rv_drift):
    fig = plt.figure(figsize=(8.27, 11.69))
    fig.text(0.06, 0.965, 'P2 BASKET — Detector 1', fontsize=30, fontweight='bold',
             color=INK, va='top')
    fig.text(0.06, 0.918, 'Cosmic-bench characterisation (long run)', fontsize=14,
             color='#555555', va='top')
    fig.text(0.06, 0.893, f'{cfg.RUN} / {cfg.SUB_RUN}',
             fontsize=9.5, color='#333333', va='top', family='monospace')
    fig.text(0.06, 0.873,
             f'mesh {rv_mesh} V · drift {rv_drift} V     ·     '
             f'connector 10 dropped ({h["pads"]} active pads)',
             fontsize=9.5, color='#333333', va='top', family='monospace')

    # headline stat cards -------------------------------------------------- #
    hax = fig.add_axes([0.06, 0.70, 0.88, 0.16]); hax.axis('off')
    hax.set_xlim(0, 1); hax.set_ylim(0, 1)
    try:
        ev = float(h['within']); ecol = '#1a7f37' if ev >= 50 else '#b26a00' if ev >= 25 else '#b3261e'
    except ValueError:
        ecol = 'black'

    def card(x, value, unit, label, color='black'):
        sep = '' if unit == '%' else ' '
        val = f'{value}{sep}{unit}' if unit else f'{value}'
        hax.text(x, 0.62, val, fontsize=19, fontweight='bold', color=color, va='center')
        hax.text(x, 0.12, label, fontsize=8.0, color='dimgrey', va='center')

    card(0.00, h['within'], '%', 'Efficiency (reco ≤20 mm)', ecol)
    card(0.21, h['has_any'], '%', 'Any pad fired')
    card(0.41, h['med_r'], 'mm', 'Median |r| residual')
    card(0.63, h['deadt'], '%', 'HV-spark deadtime')
    card(0.83, h['ped_rms'], 'ADC', 'Pedestal RMS (median)')

    # observations --------------------------------------------------------- #
    obs = (
        'KEY OBSERVATIONS (last session)\n'
        f'• Connector 10 is disconnected on P2 det1 → dropped everywhere, leaving '
        f'{h["pads"]} pads (9×128). Removing it fixed the M3 alignment: '
        f'r_x 0.92 / r_y 0.87, clean {h["rot"]}° rotation, scale {h["scale"]}, z*≈247 mm.\n'
        f'• Integrated efficiency in the fixed active area: reco ≤20 mm = {h["within"]}%, '
        f'any pad fired = {h["has_any"]}%. Loss splits into fired-not-reco {h["fnr"]}% + '
        f'silent no-hit {h["no_hit"]}%. Median |r| residual {h["med_r"]} mm (pad-pitch limited, '
        f'σ ~ 10–11 mm).\n'
        f'• A real DIAGONAL near-zero-efficiency band cuts across the fan '
        '(reference X ~ 120–160 mm) plus a few dead spots — a dead/low-gain region or '
        'mapping seam, not yet traced back to pads/connectors.\n'
        f'• HV sparks on the mesh (ch 1:0): {h["sparks"]} events, {h["deadt"]}% deadtime, '
        'vetoed. They are NOT time-correlated with DAQ high-multiplicity bursts — with an '
        'external cosmic trigger the DAQ rarely records the discharge; the veto is a small clean '
        'correction.\n'
        f'• {h["spark_pads"]} pads carry a per-pad micro-spark signature (saturation + '
        'high-amplitude outliers), clustered in connector 6/7. Pedestal QA confirms these are '
        'REAL discharges, not noisy channels (RMS uniform, r(RMS,rate)=0.06). Flagged for review, '
        'NOT masked.'
    )
    fig.text(0.06, 0.635, obs, fontsize=9.4, color='#1a1a1a', va='top', ha='left',
             wrap=True, linespacing=1.5,
             bbox=dict(boxstyle='round,pad=0.9', fc='#f4f7fb', ec='#9db8d2', lw=1.0))

    fig.text(0.06, 0.045,
             f'clean M3 rays: {h["rays"]}   (in active area: {h["rays_act"]})     '
             f'{h["hpe"]} hits/event     data: {qa.DATA_ROOT}',
             fontsize=7.5, color='#555555', family='monospace')
    fig.text(0.96, 0.008, datetime.date.today().isoformat(), ha='right',
             fontsize=6, color='grey')
    pdf.savefig(fig, dpi=200)
    plt.close(fig)


# --------------------------------------------------------------------------- #
def build(cfg, out):
    global _SUFFIXES
    base = cfg.OUT_BASE
    if not os.path.isdir(base):
        print(f'No Analysis tree at {base}')
        return
    # Resolve product filenames against this run's dead-connector set.
    _SUFFIXES = suffixes_for(cfg)
    h = headline(base, _SUFFIXES)

    def T(stem, cap):
        return (find_img(base, stem), cap)

    with PdfPages(out) as pdf:
        cover_page(pdf, cfg, h, 420, 600)

        grid_page(pdf, 'Channel map & pad-response validation',
                  'Stage 02 — Gerber pad map, occupancy/amplitude uniformity, within-connector ordering.',
                  [T('02_map_validation/pad_occupancy', 'Pad occupancy'),
                   T('02_map_validation/pad_amplitude', 'Pad amplitude (charge)'),
                   T('02_map_validation/dead_pad_map', 'Dead / connected pad map'),
                   T('02_map_validation/strategy_compare', 'Within-connector ordering: reverse vs alternatives')],
                  notes='REVERSE within-connector ordering was validated against M3 (2-pin MEC8 gap; '
                        'reverse ≠ any simple pin orientation). Connector 10 dropped → 1152 pads.')

        grid_page(pdf, 'M3 reference & P2↔M3 alignment',
                  'Stages 04 + 03 — cosmic tracker reference and the rigid pad→M3 transform.',
                  [T('04_m3_reference_qa/m3_beam_profile_detz', 'M3 track profile at detector z'),
                   T('03_m3_alignment/rotation_scan', 'Rotation scan (correlation vs angle)'),
                   T('03_m3_alignment/best_angle_corr', 'Best-angle P2 vs M3 correlation'),
                   T('03_m3_alignment/z_theta_scan', 'z–θ alignment scan')],
                  notes=f'Best rigid transform: rotation {h["rot"]}°, isotropic scale {h["scale"]}, '
                        'z*≈247 mm, corr(x)=0.92 / corr(y)=0.88. Charge-weighted pad centroids.')

        grid_page(pdf, 'Efficiency',
                  'Stage 06 — track-matched efficiency in a fixed footprint-defined active area.',
                  [T('06_efficiency/map_within_*mm', 'Per-pad efficiency map (reco within R)'),
                   T('06_efficiency/map_has_any', 'Any-pad-fired map'),
                   T('06_efficiency/efficiency_breakdown', 'Efficiency breakdown (where do muons go?)'),
                   T('06_efficiency/radial_residual', 'Radial residual P2−M3'),
                   T('06_efficiency/efficiency_vs_time', 'Efficiency vs time (30-min bins)')],
                  notes=f'Reco ≤20 mm = {h["within"]}%, any-pad = {h["has_any"]}%; loss = '
                        f'fired-not-reco {h["fnr"]}% + silent {h["no_hit"]}%. Median |r| {h["med_r"]} mm. '
                        'The vs-time panel exposes gain dropouts: when the response is intermittent, '
                        'the integrated numbers are duty-cycle averages, not working-point efficiencies.')

        grid_page(pdf, 'Sliding-window efficiency map',
                  'Stage 06 (sliding) — smooth footprint-masked efficiency; the diagonal dead band.',
                  [T('06_efficiency/efficiency_map_sliding', 'Sliding-window efficiency (within / has_any / kernel)'),
                   T('06_efficiency/nonreco_ray_positions', 'Non-reconstructed ray positions')],
                  notes='FINDING: a real diagonal near-zero-efficiency band crosses the fan '
                        '(reference X ~ 120–160 mm) plus a few dead spots — dead/low-gain region '
                        'or mapping seam, an open thread to map back to pads/connectors.',
                  ncols=1)

        grid_page(pdf, 'Detector deep QA',
                  'Stage 05 — occupancy, per-pad firing, event multiplicity and its time stability.',
                  [T('05_detector_deep_qa/surface_hitmap', 'Surface hit map'),
                   T('05_detector_deep_qa/pad_firing_fraction', 'Per-pad firing fraction'),
                   T('05_detector_deep_qa/event_multiplicity', 'Event multiplicity'),
                   T('05_detector_deep_qa/multiplicity_vs_time', 'Multiplicity vs time'),
                   T('05_detector_deep_qa/centroid_hitmap', 'Event centroid hit map')],
                  notes=f'All {h["pads"]}/1152 pads fire; {h["hpe"]} hits/event (median 1). Hottest pads '
                        'cluster in connectors 6/7; the bright corner pad ch_id 27 (conn 1) is a separate '
                        'rate outlier with only a borderline discharge signature.')

        grid_page(pdf, 'HV spark veto',
                  'Stage 07 — mesh (ch 1:0) discharge detection from hv_monitor.csv and the event veto.',
                  [T('07_hv_spark_qa/hv_current_vmon_timeline', 'Mesh current / voltage timeline'),
                   T('07_hv_spark_qa/spark_rate_evolution', 'Spark-rate evolution'),
                   T('07_hv_spark_qa/spark_amplitude_dist', 'Spark amplitude distribution'),
                   T('07_hv_spark_qa/spark_daq_crosscheck', 'HV spark vs DAQ multiplicity cross-check')],
                  notes=f'{h["sparks"]} sparks (4.1/h), {h["deadt"]}% deadtime vetoed (±2/+10 s guard). '
                        'One 270-s discharge dominates the charge. HV sparks are NOT correlated with DAQ '
                        'high-multiplicity bursts (7/294 ≈ deadtime) → veto is a small clean correction.')

        grid_page(pdf, 'Per-pad micro-spark flag  (review-only)',
                  'Stage 08 — single-pad discharges invisible to the mesh-current HV veto.',
                  [T('08_pad_spark_qa/pad_spark_map', 'Flagged spark-pad map'),
                   T('08_pad_spark_qa/flagged_per_connector', 'Flagged pads per connector'),
                   T('08_pad_spark_qa/rate_vs_amplitude', 'Rate vs amplitude (flagging plane)'),
                   T('08_pad_spark_qa/firing_before_after', 'Firing before/after (spark signature)')],
                  notes=f'{h["spark_pads"]} pads flagged by discharge signature (saturation-rate AND '
                        'high-amp>1000 ADC outliers), clustered in connector 6 (12) + c7 (5) in two tight '
                        'groups; 0 noise pads. Identified but NOT masked (user decision) — review only.')

        grid_page(pdf, 'Pedestal QA',
                  'Stage 09 — per-channel pedestal mean/RMS from the hits_root files.',
                  [T('09_pedestal_qa/pedestal_mean_map', 'Pedestal mean map'),
                   T('09_pedestal_qa/pedestal_rms_map', 'Pedestal RMS map'),
                   T('09_pedestal_qa/pedestal_rms_dist', 'Pedestal RMS distribution'),
                   T('09_pedestal_qa/pedestal_rms_vs_rate', 'Pedestal RMS vs firing rate')],
                  notes=f'RMS median {h["ped_rms"]} ADC (uniform), thresholds healthy, r(RMS,rate)=0.06 '
                        '→ pedestal noise does NOT drive the hit rate: the hot pads are real discharges, '
                        'not noisy channels. Corner pad ch_id 27 has a normal pedestal.')

    print(f'Wrote final QA PDF -> {out}')


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('-')]
    key = args[0] if args else 'det1_long'
    cfg = qa.get_config(key)
    # Name the report after the run so multiple runs of the same detector (e.g.
    # two long runs) don't overwrite each other's PDF.
    default_out = os.path.join(qa.ANALYSIS_ROOT, cfg.DET_TAG,
                               f'{cfg.RUN}_final_qa.pdf')
    out = next((a.split('=', 1)[1] for a in sys.argv if a.startswith('--out=')), default_out)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    print(cfg)
    build(cfg, out)


if __name__ == '__main__':
    main()
