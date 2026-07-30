#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
31_detector_summary_tex.py

One LaTeX overview per detector, across the whole campaign.

The other stages answer "how did this sub_run / this scan behave?". This one
answers "how did this detector behave, everywhere we ran it?" -- so it reads
only the scan-level products the pipeline already persisted and never touches
raw data. It is therefore cheap (seconds) and can be re-run after every
merge_and_pull.

Sections, per detector:

  1 overview ............ runs covered, points, HV span, which stages exist
  2 hit maps ............ 20's per-point occupancy maps
  3 efficiency maps ..... 22's per-pad efficiency maps + eff vs HV
  4 mesh scans .......... every run whose mesh voltage was stepped
  5 drift scans ......... every run whose drift voltage was stepped
  6 timing resolution ... 29's sigma, two ways:
      6a the correction ladder  raw -> fine-timestamp -> time-walk, per
         algorithm, from waveform_timing_summary.json
      6b standalone sigma vs the value implied by station-pair correlation
         (sigma_pair/sqrt(2)-style single-station numbers), against drift
         and against mesh

Scan type is read off the data, not the run name: a run is a mesh scan if its
mesh voltage takes more than one value across sub_runs, a drift scan if the
drift voltage does, and a 2D scan if both.

Usage:
  python3 31_detector_summary_tex.py                 # all detectors
  python3 31_detector_summary_tex.py --det P2_MID    # one
  python3 31_detector_summary_tex.py --pdf           # also run pdflatex
  python3 31_detector_summary_tex.py --out-dir DIR
"""

import os
import re
import csv
import glob
import json
import shutil
import argparse
import datetime
import subprocess
from collections import defaultdict

ANALYSIS = os.environ.get(
    'SPS_ANALYSIS_ROOT',
    '/local/home/banco/P2_data/TB_July2026_H4/analysis')

DETECTORS = ['P2_IN', 'P2_MID', 'P2_OUT']

# Products are written with a suffix recording whether spark-contaminated
# points were vetoed. Prefer the vetoed copy; fall back to the plain one.
SUFFIXES = ['_spark_vetoed', '']

MAX_MAPS = 6        # figures per map gallery, evenly spread across the scan
FIG_W = r'0.32\textwidth'


# --------------------------------------------------------------------------
# small helpers
# --------------------------------------------------------------------------

def tex_escape(s):
    """Escape the characters that actually occur in run/sub_run names."""
    return str(s).replace('\\', r'\textbackslash{}').replace('_', r'\_') \
                 .replace('%', r'\%').replace('&', r'\&').replace('#', r'\#')


def read_csv(path):
    """CSV -> list of dicts, or [] if absent/empty. No pandas dependency."""
    if not path or not os.path.isfile(path):
        return []
    with open(path, newline='') as fh:
        try:
            return [r for r in csv.DictReader(fh)]
        except csv.Error:
            return []


def num(row, key):
    """float(row[key]) or None -- these CSVs carry empty cells for failed fits."""
    v = (row or {}).get(key, '')
    if v is None or v == '':
        return None
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def first_existing(dirpath, stem, exts=('csv',)):
    """<stem><suffix>.<ext> for the first suffix that exists."""
    for suf in SUFFIXES:
        for ext in exts:
            p = os.path.join(dirpath, f'{stem}{suf}.{ext}')
            if os.path.isfile(p):
                return p
    return None


def spread(items, n):
    """Up to n items, evenly spread across the list (keeps both endpoints)."""
    if len(items) <= n:
        return items
    step = (len(items) - 1) / (n - 1)
    return [items[round(i * step)] for i in range(n)]


def fmt(v, nd=2, dash='--'):
    return dash if v is None else f'{v:.{nd}f}'


# --------------------------------------------------------------------------
# gathering
# --------------------------------------------------------------------------

def scan_dir(det, run, stage):
    return os.path.join(ANALYSIS, det, run, 'scan', stage)


def runs_for(det):
    """Runs that have any scan-level product for this detector."""
    out = []
    for p in sorted(glob.glob(os.path.join(ANALYSIS, det, '*', 'scan'))):
        run = p.split(os.sep)[-2]
        if any(os.path.isdir(os.path.join(p, s)) for s in
               ('20_beam_spectra', '22_tag_probe_efficiency',
                '26_hv_spark_qa', '28_timing_qa', '30_raw_stream_efficiency')):
            out.append(run)
    return out


def timing_tables(run):
    """{'drift': rows, 'mesh': rows} from 29's scan-level CSVs (telescope tree)."""
    base = os.path.join(ANALYSIS, 'telescope', run, 'scan', '29_waveform_timing')
    return {axis: read_csv(os.path.join(base, f'timing_vs_{axis}.csv'))
            for axis in ('drift', 'mesh')}


def hv_span(rows, det, field):
    """(min, max, n_distinct) of mesh_v_<det> / drift_v_<det> across rows."""
    vals = sorted({v for v in (num(r, f'{field}_v_{det}') for r in rows)
                   if v is not None})
    return (vals[0], vals[-1], len(vals)) if vals else (None, None, 0)


def all_timing_rows(run):
    """
    Both stage-29 scan tables merged, deduplicated by sub_run.

    A 2D scan writes one table per axis and each holds only that axis's slice,
    so classifying off either alone reports the other axis as constant --
    which is how a 2D run gets mislabelled as a plain drift scan.
    """
    tabs = timing_tables(run)
    seen, rows = set(), []
    for r in (tabs['drift'] or []) + (tabs['mesh'] or []):
        k = r.get('sub_run')
        if k in seen:
            continue
        seen.add(k)
        rows.append(r)
    return rows


def fallback_points(det, run):
    """
    Point count and mesh span for runs with no scan-level stage-29 table
    (single-point runs, and runs where 29 produced only per-sub_run output).
    Stage 28's CSV is the widest-covering per-detector scan row we have.
    """
    for stage, stem in (('28_timing_qa', f'timing_qa_{det}'),
                        ('26_hv_spark_qa', f'spark_qa_{det}')):
        rows = read_csv(os.path.join(scan_dir(det, run, stage), stem + '.csv'))
        if rows:
            vals = sorted({v for v in (num(r, 'mesh_v') for r in rows)
                           if v is not None})
            span = (vals[0], vals[-1], len(vals)) if vals else (None, None, 0)
            return len(rows), span
    return 0, (None, None, 0)


def gap_span(rows, det):
    """(min, max, n_distinct) of the drift-minus-mesh gap, rounded to 1 V."""
    vals = sorted({round(d - m) for d, m in (
        (num(r, f'drift_v_{det}'), num(r, f'mesh_v_{det}')) for r in rows)
        if d is not None and m is not None})
    return (vals[0], vals[-1], len(vals)) if vals else (None, None, 0)


def classify(run, det):
    """
    ('mesh'|'drift'|'2D'|'fixed', mesh_span, drift_span, gap_span).

    Keyed on the drift-minus-mesh gap, not on the drift electrode alone.
    Stepping the mesh at a fixed drift *field* means moving the drift
    electrode in lockstep, so a pure mesh scan has both voltages varying and
    would otherwise be misread as 2D -- meshscan_fine_1 holds the gap at
    250 V across all 12 points, low_mesh_scan_1 at 300 V. A genuine 2D scan
    varies the gap as well.
    """
    rows = all_timing_rows(run)
    if not rows:
        return None, (None, None, 0), (None, None, 0), (None, None, 0)
    m = hv_span(rows, det, 'mesh')
    d = hv_span(rows, det, 'drift')
    g = gap_span(rows, det)
    if m[2] > 1 and g[2] > 1:
        kind = '2D'
    elif m[2] > 1:
        kind = 'mesh'
    elif d[2] > 1:
        kind = 'drift'
    else:
        kind = 'fixed'
    return kind, m, d, g


def efficiency_rows(det, run):
    """22's scan-level rows for this detector as the probe."""
    d = scan_dir(det, run, '22_tag_probe_efficiency')
    rows = read_csv(first_existing(d, 'tag_probe_efficiency'))
    return [r for r in rows if r.get('probe') in (det, '')]


def correction_ladder(run, det):
    """
    Per-sub_run raw -> ftst -> walk sigma for this detector, from every
    waveform_timing_summary.json of the run. Returns rows sorted by sub_run.
    """
    out = []
    pat = os.path.join(ANALYSIS, 'telescope', run, '*',
                       '29_waveform_timing', 'waveform_timing_summary.json')
    for p in sorted(glob.glob(pat)):
        try:
            with open(p) as fh:
                js = json.load(fh)
        except (OSError, ValueError):
            continue
        for st in js.get('stations', []):
            if st.get('detector') != det:
                continue
            best = st.get('best_algorithm')
            algo = (st.get('algorithms') or {}).get(best, {})
            out.append({
                'sub_run': js.get('sub_run', os.path.basename(p)),
                'algo': best,
                'n_hits': st.get('n_hits'),
                'raw': algo.get('sigma_raw_ns'),
                'ftst': algo.get('sigma_ftst_ns'),
                'walk': algo.get('sigma_walk_ns'),
                'slope': algo.get('ftst_slope_ns'),
                'algos': st.get('algorithms') or {},
            })
    return out


# --------------------------------------------------------------------------
# LaTeX emission
# --------------------------------------------------------------------------

PREAMBLE = r"""\documentclass[11pt,a4paper]{article}
\usepackage[margin=20mm]{geometry}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{longtable}
\usepackage{caption}
\usepackage[table]{xcolor}
\usepackage{hyperref}
\hypersetup{colorlinks=true,linkcolor=blue!50!black,urlcolor=blue!50!black}
\captionsetup{font=small,labelfont=bf}
\setlength{\parindent}{0pt}
\setlength{\parskip}{4pt}
\graphicspath{{%s/}}
\title{%s --- performance overview\\[2pt]\large SPS July 2026 test beam (H4)}
\author{Generated by \texttt{31\_detector\_summary\_tex.py}}
\date{%s}
\begin{document}
\maketitle
\tableofcontents
\clearpage
"""


def figure_grid(out, paths, captions, width=FIG_W, cap=None):
    """A row-wrapped grid of figures with a shared caption."""
    if not paths:
        return
    out.append(r'\begin{figure}[htbp]\centering')
    for p, c in zip(paths, captions):
        rel = os.path.relpath(p, ANALYSIS)
        out.append(r'\begin{minipage}[t]{%s}\centering' % width)
        out.append(r'\includegraphics[width=\linewidth]{%s}' % rel)
        out.append(r'\\{\scriptsize %s}' % tex_escape(c))
        out.append(r'\end{minipage}\hfill')
    if cap:
        out.append(r'\caption{%s}' % cap)
    out.append(r'\end{figure}')


def single_figure(out, path, cap, width=r'0.62\textwidth'):
    if not path or not os.path.isfile(path):
        return False
    rel = os.path.relpath(path, ANALYSIS)
    out.append(r'\begin{figure}[htbp]\centering')
    out.append(r'\includegraphics[width=%s]{%s}' % (width, rel))
    out.append(r'\caption{%s}' % cap)
    out.append(r'\end{figure}')
    return True


def section_overview(out, det, runs, meta):
    out.append(r'\section{Overview}')
    out.append('This document collates every scan-level product the pipeline '
               'holds for \\texttt{%s}. Each figure is the one the stage '
               'itself produced; nothing is recomputed here.' % tex_escape(det))
    out.append('Gap is the drift-minus-mesh voltage, i.e.\\ what sets the '
               'drift field. A run that steps the mesh at constant gap is a '
               'mesh scan even though both electrodes move; only a run that '
               'varies the gap as well is 2D.')

    def span(s):
        if s[2] > 1:
            return '%g--%g' % (s[0], s[1])
        return fmt(s[0], 0) if s[2] else '--'

    out.append(r'\begin{longtable}{llllrrl}')
    out.append(r'\toprule')
    out.append(r'Run & Type & Mesh [V] & Drift [V] & Gap [V] & Points '
               r'& Stages \\')
    out.append(r'\midrule\endhead')
    for run in runs:
        m = meta[run]
        out.append(r'%s & %s & %s & %s & %s & %d & %s \\' % (
            tex_escape(run), m['kind'] or 'n/a', span(m['mesh']),
            span(m['drift']), span(m['gap']), m['n_points'],
            ', '.join(m['stages']) or '--'))
    out.append(r'\bottomrule')
    out.append(r'\end{longtable}')


def section_hit_maps(out, det, runs):
    out.append(r'\section{Hit maps}')
    out.append('Per-pad occupancy from stage 20, one map per scan point. Where '
               'a scan has many points only a spread of them is shown.')
    any_shown = False
    for run in runs:
        d = os.path.join(scan_dir(det, run, '20_beam_spectra'), 'hit_maps')
        maps = sorted(glob.glob(os.path.join(d, 'hit_map_*.png')))
        if not maps:
            continue
        any_shown = True
        out.append(r'\subsection{%s}' % tex_escape(run))
        sel = spread(maps, MAX_MAPS)
        caps = [re.sub(r'^hit_map_|_spark_vetoed|\.png$', '',
                       os.path.basename(p)) for p in sel]
        figure_grid(out, sel, caps,
                    cap='%s --- occupancy, %d of %d scan points.'
                        % (tex_escape(run), len(sel), len(maps)))
    if not any_shown:
        out.append(r'\emph{No hit maps for this detector.}')


def section_eff_maps(out, det, runs):
    out.append(r'\section{Efficiency maps and efficiency vs.\ HV}')
    out.append('Stage 22, tag-and-probe. \\texttt{eff\\_corr} is corrected for '
               'DAQ overlap using the recorded-event count; measured DAQ loss '
               'was 0.0\\,\\% throughout this campaign, so it equals '
               '\\texttt{eff}.')
    any_shown = False
    for run in runs:
        d = scan_dir(det, run, '22_tag_probe_efficiency')
        maps = sorted(glob.glob(os.path.join(d, f'eff_map_{det}_*.png')))
        curve = first_existing(d, 'tag_probe_efficiency', exts=('png',))
        rows = efficiency_rows(det, run)
        if not (maps or curve or rows):
            continue
        any_shown = True
        out.append(r'\subsection{%s}' % tex_escape(run))
        if curve:
            single_figure(out, curve,
                          '%s --- efficiency vs.\\ scan variable.'
                          % tex_escape(run))
        if rows:
            best = max(rows, key=lambda r: num(r, 'eff_corr') or -1)
            out.append(r'\begin{center}\begin{tabular}{lrrr}')
            out.append(r'\toprule')
            out.append(r'Point & HV [V] & $\varepsilon$ & $\varepsilon_{corr}$ \\')
            out.append(r'\midrule')
            for r in rows:
                mark = r'\bfseries ' if r is best else ''
                out.append(r'%s%s & %s & %s & %s \\' % (
                    mark, tex_escape(r.get('sub_run', '')),
                    fmt(num(r, 'hv'), 0), fmt(num(r, 'eff'), 4),
                    fmt(num(r, 'eff_corr'), 4)))
            out.append(r'\bottomrule\end{tabular}\end{center}')
            out.append('Best point: \\textbf{%s} at %s\\,V, '
                       '$\\varepsilon_{corr} = %s$.'
                       % (tex_escape(best.get('sub_run', '')),
                          fmt(num(best, 'hv'), 0),
                          fmt(num(best, 'eff_corr'), 4)))
        if maps:
            sel = spread(maps, MAX_MAPS)
            caps = [re.sub(r'^eff_map_%s_|_spark_vetoed|\.png$' % det, '',
                           os.path.basename(p)) for p in sel]
            figure_grid(out, sel, caps,
                        cap='%s --- per-pad efficiency, %d of %d points.'
                            % (tex_escape(run), len(sel), len(maps)))
    if not any_shown:
        out.append(r'\emph{No tag-and-probe products for this detector.}')


def _scan_block(out, det, run):
    """The four curves a scan produces, whichever of them exist."""
    figs = [
        (first_existing(scan_dir(det, run, '20_beam_spectra'),
                        'beam_mpv_vs_hv', exts=('png',)),
         'Landau MPV vs.\\ scan variable (stage 20)'),
        (first_existing(scan_dir(det, run, '20_beam_spectra'),
                        'beam_spectra_overlay', exts=('png',)),
         'Charge spectra overlay (stage 20)'),
        (first_existing(scan_dir(det, run, '22_tag_probe_efficiency'),
                        'tag_probe_efficiency', exts=('png',)),
         'Tag-and-probe efficiency (stage 22)'),
        (os.path.join(scan_dir(det, run, '30_raw_stream_efficiency'),
                      f'raw_eff_vs_hv_{det}.png'),
         'Raw-stream efficiency vs.\\ HV (stage 30)'),
        (os.path.join(scan_dir(det, run, '26_hv_spark_qa'),
                      f'spark_vs_hv_{det}.png'),
         'Spark rate and current (stage 26)'),
        (os.path.join(scan_dir(det, run, '28_timing_qa'),
                      f'timing_vs_hv_{det}.png'),
         'Hit-time mean and spread (stage 28)'),
    ]
    present = [(p, c) for p, c in figs if p and os.path.isfile(p)]
    if not present:
        out.append(r'\emph{No scan curves for this run.}')
        return
    for p, c in present:
        single_figure(out, p, '%s --- %s.' % (tex_escape(run), c),
                      width=r'0.58\textwidth')


def section_scans(out, det, runs, meta, want):
    title = {'mesh': 'Mesh scans', 'drift': 'Drift scans'}[want]
    out.append(r'\section{%s}' % title)
    sel = [r for r in runs if meta[r]['kind'] in (want, '2D')]
    if not sel:
        out.append(r'\emph{No %s for this detector.}' % title.lower())
        return
    out.append('Runs in which the %s voltage was stepped (%d). '
               'Runs marked 2D stepped both and appear in this section and '
               'the other.' % (want, len(sel)))
    for run in sel:
        m = meta[run]
        span = m[want]
        out.append(r'\subsection{%s}' % tex_escape(run))
        out.append('Type: %s. %s %g--%g\\,V over %d distinct settings.'
                   % (m['kind'], want.capitalize(), span[0], span[1], span[2]))
        _scan_block(out, det, run)


def section_timing(out, det, runs, meta):
    out.append(r'\section{Timing resolution}')
    out.append(
        'Stage 29 measures the single-station time resolution two independent '
        'ways, and both are reported here because they answer different '
        'questions.')
    out.append(
        r'\textbf{(a) The correction ladder.} For one station in isolation, '
        r'$\sigma_{raw}$ is the spread of the hit time; $\sigma_{ftst}$ adds '
        r'the FEU fine-timestamp correction; $\sigma_{walk}$ additionally '
        r'corrects amplitude-dependent time walk. The gap between them is how '
        r'much of the raw spread was instrumental rather than physical.')
    out.append(
        r'\textbf{(b) Station-pair correlation.} Taking the $\Delta t$ between '
        r'two stations cancels the common trigger jitter, so the width of that '
        r'distribution converts to a per-station number that does \emph{not} '
        r'inherit the trigger term. Where the two disagree, the standalone '
        r'value is the pessimistic one.')

    # ---- 6a: the ladder --------------------------------------------------
    out.append(r'\subsection{Correction ladder, per scan point}')
    shown = False
    for run in runs:
        ladder = correction_ladder(run, det)
        if not ladder:
            continue
        shown = True
        out.append(r'\subsubsection{%s}' % tex_escape(run))
        out.append(r'\begin{center}\begin{longtable}{lrrrrr}')
        out.append(r'\toprule')
        out.append(r'Sub-run & Algo & $\sigma_{raw}$ & $\sigma_{ftst}$ & '
                   r'$\sigma_{walk}$ & gain \\')
        out.append(r' & & [ns] & [ns] & [ns] & [\%] \\')
        out.append(r'\midrule\endhead')
        best = None
        for r in ladder:
            gain = (None if not (r['raw'] and r['walk'])
                    else 100.0 * (1.0 - r['walk'] / r['raw']))
            if r['walk'] is not None and (best is None or r['walk'] < best['walk']):
                best = r
            out.append(r'%s & %s & %s & %s & %s & %s \\' % (
                tex_escape(r['sub_run']), tex_escape(r['algo'] or '--'),
                fmt(r['raw']), fmt(r['ftst']), fmt(r['walk']), fmt(gain, 1)))
        out.append(r'\bottomrule\end{longtable}\end{center}')
        if best:
            out.append('Best in this run: \\textbf{%s\\,ns} at %s '
                       '(%s, $\\sigma_{raw} = %s$\\,ns).'
                       % (fmt(best['walk']), tex_escape(best['sub_run']),
                          tex_escape(best['algo'] or '--'), fmt(best['raw'])))
    if not shown:
        out.append(r'\emph{No stage-29 summaries for this detector.}')

    # ---- 6b: standalone vs pair, against drift and mesh -------------------
    out.append(r'\subsection{Standalone vs.\ pair-correlated, '
               r'against drift and mesh}')
    for axis, label in (('drift', 'drift'), ('mesh', 'mesh')):
        sel = [r for r in runs if timing_tables(r)[axis]]
        if not sel:
            continue
        out.append(r'\subsubsection{Versus %s voltage}' % label)
        for run in sel:
            rows = timing_tables(run)[axis]
            fig = os.path.join(ANALYSIS, 'telescope', run, 'scan',
                               '29_waveform_timing', f'timing_vs_{axis}.png')
            out.append(r'\paragraph{%s}' % tex_escape(run))
            single_figure(out, fig,
                          '%s --- resolution vs.\\ %s voltage, all stations.'
                          % (tex_escape(run), label),
                          width=r'0.66\textwidth')
            pairs = [k[len('pair_'):-len('_single')] for k in rows[0]
                     if k.startswith('pair_') and k.endswith('_single')]
            mine = [p for p in pairs if det in p]
            out.append(r'\begin{center}\begin{longtable}{lrrrl}')
            out.append(r'\toprule')
            out.append(r'Sub-run & %s [V] & $\sigma_{standalone}$ & '
                       r'$\sigma_{pair}$ & pair \\' % label.capitalize())
            out.append(r' & & [ns] & [ns] & \\')
            out.append(r'\midrule\endhead')
            for r in rows:
                solo = num(r, f'{det}_sigma')
                cand = [(num(r, f'pair_{p}_single'), p) for p in mine]
                cand = [(v, p) for v, p in cand if v is not None]
                pv, pn = min(cand) if cand else (None, '')
                out.append(r'%s & %s & %s & %s & %s \\' % (
                    tex_escape(r.get('sub_run', '')),
                    fmt(num(r, f'{axis}_v_{det}'), 0),
                    fmt(solo), fmt(pv),
                    tex_escape(pn.replace('_', '-')) if pn else '--'))
            out.append(r'\bottomrule\end{longtable}\end{center}')


# --------------------------------------------------------------------------

def build(det, out_dir):
    runs = runs_for(det)
    meta = {}
    for run in runs:
        kind, m, d, g = classify(run, det)
        stages = [s.split('_')[0] for s in
                  ('20_beam_spectra', '22_tag_probe_efficiency',
                   '26_hv_spark_qa', '28_timing_qa', '30_raw_stream_efficiency')
                  if os.path.isdir(scan_dir(det, run, s))]
        if os.path.isdir(os.path.join(ANALYSIS, 'telescope', run, 'scan',
                                      '29_waveform_timing')):
            stages.append('29')
        n_points = len(all_timing_rows(run))
        if not n_points:
            # No scan-level stage-29 table: fall back to a per-detector stage
            # so single-point runs still report their point count and HV.
            n_points, m_fb = fallback_points(det, run)
            if m[2] == 0:
                m = m_fb
            if kind is None:
                kind = 'mesh' if m[2] > 1 else ('fixed' if n_points else None)
        meta[run] = {
            'kind': kind, 'mesh': m, 'drift': d, 'gap': g,
            'n_points': n_points,
            'stages': sorted(stages, key=int),
        }

    out = [PREAMBLE % (ANALYSIS, tex_escape(det),
                       datetime.date.today().isoformat())]
    section_overview(out, det, runs, meta)
    out.append(r'\clearpage')
    section_hit_maps(out, det, runs)
    out.append(r'\clearpage')
    section_eff_maps(out, det, runs)
    out.append(r'\clearpage')
    section_scans(out, det, runs, meta, 'mesh')
    out.append(r'\clearpage')
    section_scans(out, det, runs, meta, 'drift')
    out.append(r'\clearpage')
    section_timing(out, det, runs, meta)
    out.append(r'\end{document}')

    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, f'{det}_summary.tex')
    with open(path, 'w') as fh:
        fh.write('\n'.join(out) + '\n')
    return path, len(runs)


def compile_pdf(tex_path):
    if not shutil.which('pdflatex'):
        return None, 'pdflatex not installed'
    d = os.path.dirname(tex_path)
    for _ in range(2):          # twice, so the ToC resolves
        p = subprocess.run(
            ['pdflatex', '-interaction=nonstopmode', '-halt-on-error',
             os.path.basename(tex_path)],
            cwd=d, capture_output=True, text=True)
    pdf = tex_path[:-4] + '.pdf'
    if os.path.isfile(pdf):
        return pdf, None
    tail = '\n'.join(p.stdout.strip().splitlines()[-15:])
    return None, tail


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--det', action='append', choices=DETECTORS,
                    help='detector (repeatable); default all')
    ap.add_argument('--out-dir', default=os.path.join(ANALYSIS, 'summaries'))
    ap.add_argument('--pdf', action='store_true', help='also run pdflatex')
    args = ap.parse_args()

    for det in (args.det or DETECTORS):
        path, n = build(det, args.out_dir)
        print(f'{det}: {n} run(s) -> {path}')
        if args.pdf:
            pdf, err = compile_pdf(path)
            print(f'  {"-> " + pdf if pdf else "pdflatex failed: " + err}')


if __name__ == '__main__':
    main()
