#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_beam_pdf.py

Collate the SPS beam-test analysis stages into one styled PDF per run and
station. Adapted from cosmic_bench_analysis/build_hv_scan_pdf.py: it opens with
a stat-card cover (from 20's beam_mpv_vs_hv CSV) and then embeds the curated
summary figures produced by the pipeline, so a single document per station
carries the whole story:

  cover ........ MPV / HV range / points stat cards + run notes
  telescope .... 24 event-sync QA, 21 alignment residuals + offset summary
  spectra ...... 20 spectra overlay, Landau MPV vs HV, width, rate
  efficiency ... 22 tag-and-probe efficiency vs HV (+ per-pad map)
  beam profile . 23 beam spot, x/y core profiles, spill timing

Only figures that exist are added, so the PDF builds from whatever stages have
been run. The stat cover is drawn from the 20 CSV when present, else a plain
title page is used.

Usage:
  python3 build_beam_pdf.py [run_key] [--det P2_OUT] [--out PATH]
"""

import os
import glob
import argparse
import datetime

import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import sps_config as sc

INK = '#1f3a5f'


def _first(*globs):
    """First existing file matching any of the glob patterns (in order)."""
    for g in globs:
        hits = sorted(glob.glob(g))
        if hits:
            return hits[0]
    return None


def _all(*globs):
    out = []
    for g in globs:
        out.extend(sorted(glob.glob(g)))
    return out


def cover_page(pdf, cfg, det, mpv_csv):
    fig = plt.figure(figsize=(8.27, 11.69))
    fig.text(0.06, 0.955, 'P2 BASKET — SPS beam test', fontsize=26,
             fontweight='bold', color=INK, va='top')
    fig.text(0.06, 0.912, f'Station {det.name}  ({det.det_tag})', fontsize=15,
             color='#555555', va='top')
    fig.text(0.06, 0.888, f'{cfg.RUN}', fontsize=9.5, color='#333333',
             va='top', family='monospace')

    df = None
    if mpv_csv and os.path.isfile(mpv_csv):
        df = pd.read_csv(mpv_csv)

    hax = fig.add_axes([0.06, 0.78, 0.88, 0.08]); hax.axis('off')
    hax.set_xlim(0, 1); hax.set_ylim(0, 1)

    def card(x, value, unit, label):
        sep = '' if unit in ('%', '') else ' '
        hax.text(x, 0.66, f'{value}{sep}{unit}', fontsize=18, fontweight='bold',
                 color=INK, va='center')
        hax.text(x, 0.10, label, fontsize=8, color='dimgrey', va='center')

    if df is not None and df['mpv_adc'].notna().any():
        ok = df[df['mpv_adc'].notna()]
        ipk = ok['mpv_adc'].idxmax()
        card(0.00, f'{ok.loc[ipk, "mpv_adc"]:.0f}', 'ADC', 'peak Landau MPV')
        if ok['hv'].notna().any():
            card(0.28, f'{ok.loc[ipk, "hv"]:.0f}', 'V', 'at mesh HV')
            card(0.50, f'{int(df["hv"].min())}–{int(df["hv"].max())}', 'V',
                 'mesh HV range')
        card(0.78, f'{len(df)}', '', 'scan points')

    lines = [
        'CONTENTS',
        ' - Event-sync QA across the telescope FEUs (stage 24)',
        ' - Mutual plane alignment, residuals + offsets (stage 21)',
        ' - Beam cluster spectra + Landau (moyal) MPV vs HV (stage 20)',
        ' - Tag-and-probe efficiency vs HV (stage 22)',
        ' - Beam spot, profiles and spill timing (stage 23)',
        '',
        'NOTES',
        f' - {cfg.NOTE}' if cfg.NOTE else '',
        ' - Efficiency (stage 22) is relative to the tag selection and the',
        '   overlap acceptance only, not an absolute detector efficiency.',
        ' - "MPV" is a scipy.stats.moyal (Landau approximation) location fit.',
    ]
    fig.text(0.06, 0.70, '\n'.join(l for l in lines if l is not None),
             fontsize=10, color='#1a1a1a', va='top', ha='left', linespacing=1.5,
             family='monospace')
    fig.text(0.96, 0.02, datetime.date.today().isoformat(), ha='right',
             fontsize=7, color='grey')
    pdf.savefig(fig, dpi=150)
    plt.close(fig)


def image_page(pdf, png, caption):
    if not png or not os.path.isfile(png):
        return
    img = plt.imread(png)
    h, w = img.shape[:2]
    fig = plt.figure(figsize=(8.27, 11.69))
    ax = fig.add_axes([0.04, 0.06, 0.92, 0.86])
    ax.imshow(img); ax.axis('off')
    fig.text(0.06, 0.965, caption, fontsize=12, color=INK, fontweight='bold',
             va='top')
    fig.text(0.96, 0.02, os.path.basename(png), ha='right', fontsize=6,
             color='grey', family='monospace')
    pdf.savefig(fig, dpi=150)
    plt.close(fig)


def build_for_det(cfg, det, out):
    A = cfg.ANALYSIS_ROOT
    run = cfg.RUN
    subs = cfg.find_subruns()
    ref_sub = subs[0] if subs else ''
    prod_sub = 'scan' if len(subs) > 1 else ref_sub
    dtag = det.det_tag
    ttag = sc.TELESCOPE_TAG

    mpv_csv = _first(f'{A}/{dtag}/{run}/{prod_sub}/20_beam_spectra/'
                     'beam_mpv_vs_hv*.csv')

    with PdfPages(out) as pdf:
        cover_page(pdf, cfg, det, mpv_csv)

        # telescope-wide
        image_page(pdf, _first(f'{A}/{ttag}/{run}/{ref_sub}/24_event_sync_qa/'
                               'event_sync_*.png'),
                   'Event-sync QA (telescope)')
        for p in _all(f'{A}/{ttag}/{run}/{ref_sub}/21_telescope_align/'
                      'residuals_*.png'):
            image_page(pdf, p, 'Telescope alignment residuals')
        image_page(pdf, _first(f'{A}/{ttag}/{run}/{ref_sub}/21_telescope_align/'
                               'alignment*.png'),
                   'Telescope plane offsets')

        # spectra (stage 20)
        for name, cap in [('beam_spectra_overlay', 'Beam spectra overlay'),
                          ('beam_mpv_vs_hv', 'Landau MPV vs mesh HV'),
                          ('beam_width_vs_hv', 'Landau width vs mesh HV'),
                          ('beam_rate_vs_hv', 'Event rate vs mesh HV')]:
            image_page(pdf, _first(f'{A}/{dtag}/{run}/{prod_sub}/'
                                   f'20_beam_spectra/{name}*.png'),
                       cap)

        # efficiency (stage 22)
        image_page(pdf, _first(f'{A}/{dtag}/{run}/{prod_sub}/'
                               '22_tag_probe_efficiency/'
                               'tag_probe_efficiency*.png'),
                   'Tag-and-probe efficiency')
        for p in _all(f'{A}/{dtag}/{run}/{prod_sub}/22_tag_probe_efficiency/'
                      'eff_map_*.png')[:2]:
            image_page(pdf, p, 'Tag-and-probe per-pad efficiency map')

        # beam profile (stage 23) for the reference sub_run
        for name, cap in [('beam_spot', 'Beam spot'),
                          ('profiles', 'Beam profile projections'),
                          ('timing', 'Rate vs time / pile-up')]:
            image_page(pdf, _first(f'{A}/{dtag}/{run}/{ref_sub}/'
                                   f'23_beam_profile/{name}_*.png'),
                       cap)

    print(f'Wrote {out}')


def main():
    ap = argparse.ArgumentParser(description='Collate SPS beam stages into a '
                                             'per-station PDF.')
    ap.add_argument('run_key', nargs='?', default=sc.DEFAULT_RUN)
    ap.add_argument('--det', default=None,
                    help='station name/det_tag (default: all stations).')
    ap.add_argument('--out', default=None,
                    help='output PDF path (single-station only).')
    args = ap.parse_args()

    cfg = sc.get_config(args.run_key)
    print(cfg)
    dets = cfg.detectors()
    if args.det:
        dets = [d for d in dets if args.det in (d.name, d.det_tag)]
    for det in dets:
        out = args.out or os.path.join(cfg.ANALYSIS_ROOT, det.det_tag,
                                       f'beam_report_{det.det_tag}_{cfg.RUN}.pdf')
        os.makedirs(os.path.dirname(out), exist_ok=True)
        build_for_det(cfg, det, out)


if __name__ == '__main__':
    main()
