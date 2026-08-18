#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gather_figures.py -- collect every figure the MPGD26 talk uses into one place.

The figures live in five different trees (bench stage outputs, the beam report,
the gas study, mpgd2026/make_talk_figs.py, and the lxplus-rendered event
displays).  The deck on the CERN site holds a *copy* of each, with no record of
where it came from -- so a regenerated stage output silently diverges from the
slide.  This script is that record: one manifest, source of truth per figure,
re-runnable.

Layout written under figures/:

    1_intro/  2_act1_nov2025/  3_act2_bench/  4_act3_beam/
    5_optional/  6_backup/  9_spare/          <- gathered but not on a slide

Names are prefixed with the draft-deck slide number (s05_, s11_, ...) so the
directory reads in talk order.  PDFs are copied alongside when the source tree
has one -- prefer them for the slides, the PNGs are for the web deck.

Sources are stale-checked: if a source file is newer than the copy here, the
copy is refreshed and the change reported.

Usage:  python3 gather_figures.py [--check] [--figures DIR]
          --check   report what would change, copy nothing
"""

import argparse
import os
import shutil
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SACLAY = os.path.abspath(os.path.join(HERE, '..', '..'))

# The lxplus-rendered 3D event displays exist only as the copies embedded in
# the deck -- they are produced on lxplus (nTof_x17/mpgd26) from the merged
# ROOT file, not by anything in this repo.  Flagged in the inventory.
DECK = os.path.join(
    SACLAY, 'cern-site/notes/2026-08-17-mpgd26-talk-draft-deck_files')

# (act_dir, deck slide no., tag, figure name, slide title, source path)
MANIFEST = [
    ('1_intro', 5, 'CORE', 'pad_sector_layout.png', 'The detector: a metallic, non-resistive bulk Micromegas with a wedge pad anode',
     'SACLAY/data/SPS_Beam_Test/VMM-alinx-data/results/Full_detector_model_class_from_gerber.png'),
    ('2_act1_nov2025', 8, 'CORE', 'snr_matrix.png', 'Nov 2025: the VMM3a shaping optimum',
     'SACLAY/P2_basket_analysis/mpgd2026/figs/snr_matrix.png'),
    ('3_act2_bench', 11, 'CORE', '01_det1_efficiency_map.png', 'Bench efficiency, and where the missing few percent actually is',
     'SACLAY/data/Cosmic_Bench/Analysis/conference/01_det1_efficiency_map.png'),
    ('3_act2_bench', 11, 'CORE', '11_pillar_accounting.png', 'Bench efficiency, and where the missing few percent actually is',
     'SACLAY/data/Cosmic_Bench/Analysis/conference/11_pillar_accounting.png'),
    ('3_act2_bench', 12, 'CORE', '05_all_efficiency.png', 'Four detectors, one method — and the honesty that method forces',
     'SACLAY/data/Cosmic_Bench/Analysis/conference/05_all_efficiency.png'),
    ('3_act2_bench', 13, 'OPTIONAL', '02_det1_mesh_scan.png', 'Bench HV scans: the working point is on a plateau',
     'SACLAY/data/Cosmic_Bench/Analysis/conference/02_det1_mesh_scan.png'),
    ('3_act2_bench', 13, 'OPTIONAL', '03_det1_drift_scan.png', 'Bench HV scans: the working point is on a plateau',
     'SACLAY/data/Cosmic_Bench/Analysis/conference/03_det1_drift_scan.png'),
    ('3_act2_bench', 14, 'CORE', 'lifetime_master.png', 'What the bench caught before these detectors ever saw a beam',
     'SACLAY/P2_basket_analysis/reports/det_lifetime_autopsy_2026-07/figs/lifetime_master.png'),
    ('3_act2_bench', 15, 'OPTIONAL', '13_fem_planarity_vs_timing.png', 'Bench timing — and one hypothesis killed by simulation',
     'SACLAY/data/Cosmic_Bench/Analysis/conference/13_fem_planarity_vs_timing.png'),
    ('4_act3_beam', 17, 'CORE', 'urw_panels_P2_MID.png', 'The absolute measurement, in one figure',
     'SACLAY/P2_basket_analysis/reports/mpgd26_sps_beam_2026-08/figs/urw_panels_P2_MID.png'),
    ('4_act3_beam', 18, 'CORE', 'coincidence_side.png', 'One muon, five chambers',
     'SACLAY/cern-site/notes/2026-08-17-mpgd26-talk-draft-deck_files/coincidence_side.png'),
    ('4_act3_beam', 19, 'OPTIONAL', 'coincidence_hero.png', 'The same event, two other ways to look at it',
     'SACLAY/cern-site/notes/2026-08-17-mpgd26-talk-draft-deck_files/coincidence_hero.png'),
    ('4_act3_beam', 19, 'OPTIONAL', 'event_display.png', 'The same event, two other ways to look at it',
     'SACLAY/P2_basket_analysis/reports/mpgd26_sps_beam_2026-08/figs/event_display.png'),
    ('4_act3_beam', 20, 'CORE', 'padmap_working_point.png', '96–97 % absolute efficiency — and a full account of the rest',
     'SACLAY/P2_basket_analysis/reports/mpgd26_sps_beam_2026-08/figs/padmap_working_point.png'),
    ('4_act3_beam', 21, 'CORE', 'bench_beam_mesh.png', 'ε(HV): the bench predicted the beam',
     'SACLAY/P2_basket_analysis/mpgd2026/figs/bench_beam_mesh.png'),
    ('4_act3_beam', 22, 'CORE', 'timing_campaigns.png', 'Timing: the 20 ns goal is met on two of three stations',
     'SACLAY/P2_basket_analysis/mpgd2026/figs/timing_campaigns.png'),
    ('4_act3_beam', 22, 'CORE', 'timing_vs_drift_magboltz.png', 'Timing: the 20 ns goal is met on two of three stations',
     'SACLAY/P2_basket_analysis/reports/mpgd26_sps_beam_2026-08/figs/timing_vs_drift_magboltz.png'),
    ('4_act3_beam', 23, 'CORE', 'dream_vs_vmm.png', 'DREAM vs VMM3a on the same chambers',
     'SACLAY/P2_basket_analysis/mpgd2026/figs/dream_vs_vmm.png'),
    ('4_act3_beam', 24, 'CORE', 'vmm_threshold.png', 'That gap is the discriminator threshold — three independent handles',
     'SACLAY/P2_basket_analysis/mpgd2026/figs/vmm_threshold.png'),
    ('5_optional', 26, 'OPTIONAL', 'eff_vs_mesh_urw.png', 'ε(HV) on the beam, on its own axes',
     'SACLAY/P2_basket_analysis/reports/mpgd26_sps_beam_2026-08/figs/eff_vs_mesh_urw.png'),
    ('5_optional', 26, 'OPTIONAL', 'eff_vs_drift_urw.png', 'ε(HV) on the beam, on its own axes',
     'SACLAY/P2_basket_analysis/reports/mpgd26_sps_beam_2026-08/figs/eff_vs_drift_urw.png'),
    ('5_optional', 27, 'OPTIONAL', 'timing_ladder.png', 'How the timing number is actually made',
     'SACLAY/P2_basket_analysis/reports/mpgd26_sps_beam_2026-08/figs/timing_ladder.png'),
    ('5_optional', 27, 'OPTIONAL', 'timing_vs_mesh_bothruns.png', 'How the timing number is actually made',
     'SACLAY/P2_basket_analysis/reports/mpgd26_sps_beam_2026-08/figs/timing_vs_mesh_bothruns.png'),
    ('5_optional', 28, 'OPTIONAL', 'hv_campaign_timeline.png', 'The campaign in one picture',
     'SACLAY/P2_basket_analysis/reports/mpgd26_sps_beam_2026-08/figs/hv_campaign_timeline.png'),
    ('5_optional', 29, 'OPTIONAL', 'eff_vs_time_highstat_eff_1.png', 'Charging-up, caught in the act',
     'SACLAY/P2_basket_analysis/reports/mpgd26_sps_beam_2026-08/figs/eff_vs_time_highstat_eff_1.png'),
    ('5_optional', 30, 'OPTIONAL', 'vmm_turnon_gasAB.png', 'The gas change — predicted, then observed',
     'SACLAY/P2_basket_analysis/reports/mpgd26_sps_beam_2026-08/figs/vmm_turnon_gasAB.png'),
    ('5_optional', 31, 'OPTIONAL', 'eff_2d_drift_mesh_2d_2.png', 'The full (mesh, drift) efficiency surface',
     'SACLAY/P2_basket_analysis/reports/mpgd26_sps_beam_2026-08/figs/eff_2d_drift_mesh_2d_2.png'),
    ('5_optional', 31, 'OPTIONAL', 'timing_2d_drift_mesh_2d_2.png', 'The full (mesh, drift) efficiency surface',
     'SACLAY/P2_basket_analysis/reports/mpgd26_sps_beam_2026-08/figs/timing_2d_drift_mesh_2d_2.png'),
    ('5_optional', 32, 'OPTIONAL', 'eff_methods_comparison.png', 'Three methods, one answer',
     'SACLAY/P2_basket_analysis/reports/mpgd26_sps_beam_2026-08/figs/eff_methods_comparison.png'),
    ('5_optional', 33, 'OPTIONAL', 'spark_summary.png', 'Sparks and HV health across twelve days',
     'SACLAY/P2_basket_analysis/reports/mpgd26_sps_beam_2026-08/figs/spark_summary.png'),
    ('5_optional', 34, 'OPTIONAL', 'eff_p2in_swap.png', 'Two P2_IN chambers, one beam — including the CERN-built one',
     'SACLAY/P2_basket_analysis/reports/mpgd26_sps_beam_2026-08/figs/eff_p2in_swap.png'),
    ('5_optional', 35, 'OPTIONAL', 'vmm_config_scan.png', 'The VMM configuration scan — what to actually set',
     'SACLAY/P2_basket_analysis/reports/mpgd26_sps_beam_2026-08/figs/vmm_config_scan.png'),
    ('5_optional', 36, 'OPTIONAL', '07_all_mesh_scans.png', 'Bench detail: all four detectors, every scan',
     'SACLAY/data/Cosmic_Bench/Analysis/conference/07_all_mesh_scans.png'),
    ('5_optional', 36, 'OPTIONAL', '08_all_drift_scans.png', 'Bench detail: all four detectors, every scan',
     'SACLAY/data/Cosmic_Bench/Analysis/conference/08_all_drift_scans.png'),
    ('5_optional', 36, 'OPTIONAL', '06_coverage_table.png', 'Bench detail: all four detectors, every scan',
     'SACLAY/data/Cosmic_Bench/Analysis/conference/06_coverage_table.png'),
    ('6_backup', 40, 'BACKUP', '04_det1_timing_vs_drift.png', 'Q: "Why is your bench timing so much worse than your beam timing?"',
     'SACLAY/data/Cosmic_Bench/Analysis/conference/04_det1_timing_vs_drift.png'),
    ('6_backup', 41, 'BACKUP', '10_det1_loss_budget.png', 'Q: "How do you know the inefficiency is not the gas?"',
     'SACLAY/data/Cosmic_Bench/Analysis/conference/10_det1_loss_budget.png'),
    ('6_backup', 43, 'BACKUP', 'gas_drift_velocity.png', 'Reference: the gas study behind the mixture choice',
     'SACLAY/P2_basket_analysis/gas_studies/figs/drift_velocity_vs_dV.png'),
    ('6_backup', 43, 'BACKUP', 'gas_timing_floor.png', 'Reference: the gas study behind the mixture choice',
     'SACLAY/P2_basket_analysis/gas_studies/figs/timing_floor_vs_driftHV.png'),
    ('6_backup', 44, 'BACKUP', 'drift_ladder_trajectory.png', 'Reference: the det3 autopsy in one figure',
     'SACLAY/P2_basket_analysis/reports/det_lifetime_autopsy_2026-07/figs/drift_ladder_trajectory.png'),
    # --- new for this talk, not yet placed on a slide (2026-08-17) --------- #
    ('4_act3_beam', 0, 'NEW', 'bench_beam_drift.png',
     'WP-C, transport half: eps vs drift FIELD, bench and beam on one axis',
     'SACLAY/P2_basket_analysis/mpgd2026/figs/bench_beam_drift.png'),
    ('4_act3_beam', 0, 'NEW', 'bench_beam_maps.png',
     'wish-list 5a.4: bench efficiency map beside the beam map, same scale',
     'SACLAY/P2_basket_analysis/mpgd2026/figs/bench_beam_maps.png'),
    ('5_optional', 0, 'NEW', 'eff_2d_curves_drift_mesh_2d_2.png',
     'the (mesh, drift) efficiency surface as curves instead of boxes',
     'SACLAY/P2_basket_analysis/reports/mpgd26_sps_beam_2026-08/figs/'
     'eff_2d_curves_drift_mesh_2d_2.png'),
    ('5_optional', 0, 'NEW', 'timing_2d_curves_drift_mesh_2d_2.png',
     'the (mesh, drift) timing surface as curves instead of boxes',
     'SACLAY/P2_basket_analysis/reports/mpgd26_sps_beam_2026-08/figs/'
     'timing_2d_curves_drift_mesh_2d_2.png'),
    ('5_optional', 0, 'NEW', 'timing_2d_curves_vs_drift_drift_mesh_2d_2.png',
     'the timing surface transposed: sigma vs drift voltage, one curve per mesh setting',
     'SACLAY/P2_basket_analysis/reports/mpgd26_sps_beam_2026-08/figs/'
     'timing_2d_curves_vs_drift_drift_mesh_2d_2.png'),
    ('5_optional', 0, 'NEW', 'timing_drift_choice_table_drift_mesh_2d_2.png',
     'table figure: sigma at the 150 V gap, at the 250 V working point, and the scan best, per station',
     'SACLAY/P2_basket_analysis/reports/mpgd26_sps_beam_2026-08/figs/'
     'timing_drift_choice_table_drift_mesh_2d_2.png'),
    ('4_act3_beam', 0, 'NEW', 'p2_standalone_tracking.png',
     '5b: the P2-only track vs the uRWELL reference - localises to the pad, cannot measure angle',
     'SACLAY/P2_basket_analysis/reports/mpgd26_sps_beam_2026-08/figs/'
     'p2_standalone_tracking.png'),
    ('4_act3_beam', 0, 'NEW', 'p2_standalone_tracking_explained.png',
     '5b, told as geometry: the beam cone inside one pad column, the pad-step '
     'staircase, and the crossing point in the pad box',
     'SACLAY/P2_basket_analysis/reports/mpgd26_sps_beam_2026-08/figs/'
     'p2_standalone_tracking_explained.png'),
    ('4_act3_beam', 0, 'NEW', 'rate_performance.png',
     '5b.5: efficiency and fake-match rate vs beam load, within sub-runs',
     'SACLAY/P2_basket_analysis/reports/mpgd26_sps_beam_2026-08/figs/'
     'rate_performance.png'),
    ('4_act3_beam', 0, 'NEW', 'twotrack_side.png',
     '5b: the same muon, two independent tracks - uRWELL reference (red) and the P2-only fit (blue)',
     'SACLAY/nTof_x17/mpgd26/figures/twotrack_side_light.png'),
    ('4_act3_beam', 0, 'NEW', 'twotrack_hero.png',
     '5b: the two-track event display, hero angle',
     'SACLAY/nTof_x17/mpgd26/figures/twotrack_hero_light.png'),
    ('9_spare', 0, 'SPARE', '12_charging_up.png', '',
     'SACLAY/data/Cosmic_Bench/Analysis/conference/12_charging_up.png'),
    ('9_spare', 0, 'SPARE', 'coincidence_beam.png', '',
     'SACLAY/cern-site/notes/2026-08-17-mpgd26-talk-draft-deck_files/coincidence_beam.png'),
    ('9_spare', 0, 'SPARE', 'eff_drift_highstat.png', '',
     'SACLAY/P2_basket_analysis/reports/mpgd26_sps_beam_2026-08/figs/eff_drift_highstat.png'),
    ('9_spare', 0, 'SPARE', 'eff_vs_time_eff_drift_ab_1.png', '',
     'SACLAY/P2_basket_analysis/reports/mpgd26_sps_beam_2026-08/figs/eff_vs_time_eff_drift_ab_1.png'),
    ('9_spare', 0, 'SPARE', 'eff_vs_time_eff_nominal_1.png', '',
     'SACLAY/P2_basket_analysis/reports/mpgd26_sps_beam_2026-08/figs/eff_vs_time_eff_nominal_1.png'),
    ('9_spare', 0, 'SPARE', 'timing_vs_mesh.png', '',
     'SACLAY/P2_basket_analysis/reports/mpgd26_sps_beam_2026-08/figs/timing_vs_mesh.png'),
    ('9_spare', 0, 'SPARE', 'urw_summary_highstat.png', '',
     'SACLAY/P2_basket_analysis/reports/mpgd26_sps_beam_2026-08/figs/urw_summary_highstat.png'),
    ('9_spare', 0, 'SPARE', 'vmm_drift_gasAB.png', '',
     'SACLAY/P2_basket_analysis/reports/mpgd26_sps_beam_2026-08/figs/vmm_drift_gasAB.png'),
    ('9_spare', 0, 'SPARE', 'vmm_vs_dream_caveat.png', '',
     'SACLAY/P2_basket_analysis/reports/mpgd26_sps_beam_2026-08/figs/vmm_vs_dream_caveat.png'),
    ('9_spare', 0, 'SPARE', '09_all_timing_vs_drift.png', '',
     'SACLAY/data/Cosmic_Bench/Analysis/conference/09_all_timing_vs_drift.png'),
]


def resolve(src):
    """Manifest paths are SACLAY-relative; the deck copies are the fallback."""
    return os.path.join(SACLAY, src[len('SACLAY/'):])


def gather(figdir, check):
    copied, refreshed, missing, unchanged = [], [], [], []
    for act, slide, tag, name, title, src in MANIFEST:
        srcp = resolve(src)
        stem, ext = os.path.splitext(name)
        prefix = f's{slide:02d}_' if slide else ''
        dstdir = os.path.join(figdir, act)
        dst = os.path.join(dstdir, prefix + name)

        if not os.path.exists(srcp):
            missing.append((name, src))
            continue

        if not os.path.exists(dst):
            state, bucket = 'new', copied
        elif os.path.getmtime(srcp) > os.path.getmtime(dst):
            state, bucket = 'stale', refreshed
        else:
            state, bucket = 'ok', unchanged
        bucket.append(name)
        if state != 'ok' and not check:
            os.makedirs(dstdir, exist_ok=True)
            shutil.copy2(srcp, dst)

        # the vector version, when the producing stage wrote one
        pdf = os.path.splitext(srcp)[0] + '.pdf'
        if os.path.exists(pdf):
            dpdf = os.path.join(dstdir, prefix + stem + '.pdf')
            if not check and (not os.path.exists(dpdf)
                              or os.path.getmtime(pdf) > os.path.getmtime(dpdf)):
                os.makedirs(dstdir, exist_ok=True)
                shutil.copy2(pdf, dpdf)
    return copied, refreshed, missing, unchanged


def write_inventory(figdir):
    lines = [
        '# MPGD26 talk figures — inventory',
        '',
        f'Generated by `gather_figures.py` on '
        f'{time.strftime("%Y-%m-%d %H:%M")}.  Do not edit by hand; edit the',
        'manifest in the script and re-run.  Paths are relative to',
        '`~/Documents/PostDocSaclay/`.',
        '',
    ]
    acts = {}
    for row in MANIFEST:
        acts.setdefault(row[0], []).append(row)
    for act in sorted(acts):
        lines += [f'## {act}', '',
                  '| slide | tag | figure | source |', '|---|---|---|---|']
        for _, slide, tag, name, title, src in acts[act]:
            s = str(slide) if slide else '—'
            vec = '' if os.path.exists(
                os.path.splitext(resolve(src))[0] + '.pdf') else ' *(png only)*'
            lines.append(f'| {s} | {tag} | `{name}`{vec} | `{src[7:]}` |')
        lines.append('')
        for _, slide, tag, name, title, src in acts[act]:
            if title:
                lines.append(f'- **{name}** — slide {slide}: {title}')
        lines.append('')
    with open(os.path.join(os.path.dirname(figdir), 'INVENTORY.md'), 'w') as f:
        f.write('\n'.join(lines))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--figures', default=os.path.join(HERE, 'figures'))
    ap.add_argument('--check', action='store_true')
    a = ap.parse_args()

    new, stale, missing, ok = gather(a.figures, a.check)
    verb = 'would copy' if a.check else 'copied'
    print(f'{len(MANIFEST)} figures in the manifest -> {a.figures}')
    print(f'  {verb:>10}: {len(new)}   refreshed (source newer): {len(stale)}'
          f'   up to date: {len(ok)}')
    for n in stale:
        print(f'    stale -> {n}')
    if missing:
        print(f'  MISSING {len(missing)}:')
        for n, s in missing:
            print(f'    {n}  <-  {s}')
    if not a.check:
        write_inventory(a.figures)
        print('  wrote INVENTORY.md')


if __name__ == '__main__':
    main()
