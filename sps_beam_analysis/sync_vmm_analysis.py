#!/usr/bin/env python3
"""
Organise the VMM analysis of the July 2026 SPS test on the LaCie, next to (and
separate from) the DREAM tree that sync_dream_analysis.py builds.

The VMM results are not per station at the source: one capture reads all 32
VMM chips, so every EOS store (`analysis/vmm/<run>/<sub_run>/<capture>/`
counts.npz + meta.json + scalars.json) and every QA page covers P2_IN, P2_MID
and P2_OUT at once. Those stay per run. What is per detector are the derived
per-station tables, split here by chamber.

Local layout (default root: the LaCie TB_July26_H4_vmm/analysis/):
  captures/<run>/<sub_run>/<capture>/   EOS VMM stores (all stations)
  qa_plots/<run>/<sub_run>/             QA pages (all stations; render_from_store.py)
  tables/                               the full per-sub_run / per-capture tables
  det1/ det3/ det5/                     per-chamber rows of the per-sub_run and
                                        per-capture tables (efficiency, timing),
                                        + station.txt
  combined_dream_vmm/                   DREAM vs VMM on the same beam data:
      dream_vmm_comparison/             figures, data, reports, README
      vmm_dream_matching/               tables, notes, arrays
  INDEX.md

All VMM runs (run_2 .. run_71, 29 Jul - 3 Aug) are after the det5 install, so
P2_IN = det5, P2_MID = det1, P2_OUT = det3 throughout (chamber_history.py).
Gas: Ar/CO2/iso 93/5/2, Ar/CF4/iso 88/10/2 from the 3 Aug runs (see the tables).

    python3 sps_beam_analysis/sync_vmm_analysis.py [--dest DIR] [--skip-eos]
"""
import argparse
import datetime as dt
import os
import shutil
import subprocess
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import chamber_history as ch  # noqa: E402

HOST = 'lxplus'
EOS_VMM = '/eos/experiment/ntof/data/x17/p2_sps_july/analysis/vmm'
DEST = ('/media/ak271430/LaCie/Extras/Physics/Post-Doc-Saclay/data/'
        'SPS_beam_test/TB_July26_H4_vmm/analysis')
TABLES = os.path.expanduser('~/Documents/PostDocSaclay/data/SPS_Beam_Test/'
                            'mpgd26_workspace/report_data')
TABLE_FILES = ('vmm_subruns.csv', 'vmm_captures.csv')      # mpgd26_figs/aggregate.py
COMPARISON = {
    'dream_vmm_comparison': ['figures', 'data', 'README.md', 'report.html',
                             'report_dream_vmm.html', 'report_pillar.html',
                             'report_slide.html', 'report_track.html'],
    'vmm_dream_matching': ['README.md', 'COMPARISON_CLAIMS.md', 'EFFICIENCY_AUTOPSY.md',
                           'EFFICIENCY_TABLE.md', 'MATCHING_TABLE.md', 'TIMING.md',
                           'HANDOFF_VMM_PAD_ADC.md', 'efficiency_conditions.csv',
                           'efficiency_table.csv', 'matching_table.csv',
                           'vmm_timing_by_subrun.csv', 'vmm_trigger_timing_by_hv.csv',
                           'vmm_trigger_timing.csv', 'thresholds_run_46.json',
                           'model_P2_OUT.json', 'timing_peaks_gasB.npz',
                           'timing_peaks_nomA.npz', 'timing_peaks_outA.npz',
                           'timing_peaks_run33_gap150.npz',
                           'vmm_adc_P2_MID_run46.npz', 'vmm_adc_P2_OUT_run46.npz'],
}


def chamber_of(station, run):
    """All VMM runs are 29 Jul - 3 Aug: P2_IN = det5 (in from 28 Jul)."""
    return 'det5' if station == 'P2_IN' else ch.chamber(station, run)


def copy_item(src, dst):
    if os.path.isdir(src):
        shutil.copytree(src, dst, dirs_exist_ok=True)
    elif os.path.isfile(src):
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copy2(src, dst)
    else:
        print(f'  (missing: {src})')


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[1])
    ap.add_argument('--dest', default=DEST)
    ap.add_argument('--skip-eos', action='store_true', help='do not pull the captures')
    args = ap.parse_args()
    d = args.dest

    if not args.skip_eos:
        print('captures <- EOS analysis/vmm/')
        os.makedirs(os.path.join(d, 'captures'), exist_ok=True)
        subprocess.run(['rsync', '-a', '--exclude=*.sys.v#*', f'{HOST}:{EOS_VMM}/',
                        os.path.join(d, 'captures') + '/'], check=True)

    old = os.path.join(d, 'vmm_qa_plots')            # the first layout: move, once
    if os.path.isdir(old) and not os.path.exists(os.path.join(d, 'qa_plots')):
        os.rename(old, os.path.join(d, 'qa_plots'))
        print('vmm_qa_plots -> qa_plots')

    split = {}
    os.makedirs(os.path.join(d, 'tables'), exist_ok=True)
    for name in TABLE_FILES:
        shutil.copy2(os.path.join(TABLES, name), os.path.join(d, 'tables', name))
        t = pd.read_csv(os.path.join(TABLES, name))
        # capture-level rows ('mask', 'stations') belong to no station: tables/ only
        t = t[t['det'].isin(['P2_IN', 'P2_MID', 'P2_OUT'])].copy()
        t['chamber'] = [chamber_of(s, r) for s, r in zip(t['det'], t['run'])]
        for chb, g in t.groupby('chamber'):
            out = os.path.join(d, chb)
            os.makedirs(out, exist_ok=True)
            g.to_csv(os.path.join(out, name), index=False)
            split.setdefault(chb, {})[name] = (len(g), sorted(g['det'].unique()))
    for chb, info in split.items():
        st = sorted({s for _, ss in info.values() for s in ss})
        with open(os.path.join(d, chb, 'station.txt'), 'w') as f:
            f.write(f'station: {", ".join(st)}\nchamber: {chb}  (chamber_history.py; all '
                    f'VMM runs are 29 Jul - 3 Aug)\nsource:  {TABLES}/'
                    f'{{{",".join(TABLE_FILES)}}} (mpgd26_figs/aggregate.py), rows of '
                    f'this chamber\n')
        print(f'{chb}: ' + ', '.join(f'{n} {k} rows' for k, (n, _) in info.items()))

    for pkg, items in COMPARISON.items():
        for it in items:
            copy_item(os.path.join(HERE, pkg, it),
                      os.path.join(d, 'combined_dream_vmm', pkg, it))
        print(f'combined_dream_vmm/{pkg} <- sps_beam_analysis/{pkg}')

    write_index(d, split)


def write_index(d, split):
    runs = sorted(os.listdir(os.path.join(d, 'captures'))) if os.path.isdir(
        os.path.join(d, 'captures')) else []
    lines = ['# July 2026 SPS test -- VMM analysis', '',
             f'Built by `sps_beam_analysis/sync_vmm_analysis.py`, last run '
             f'{dt.datetime.now():%Y-%m-%d %H:%M}. The DREAM analysis is in '
             '`../../TB_July2026_H4/analysis/` (per detector).', '',
             'The VMM stores and QA pages cover all three stations at once (one capture '
             'reads all 32 VMM chips), so they stay per run; the per-station tables are '
             'split per chamber.', '',
             '| folder | content |', '|---|---|',
             f'| `captures/` | EOS `analysis/vmm/`: <run>/<sub_run>/<capture>/ counts.npz, '
             f'meta.json, scalars.json, all stations; {len(runs)} runs |',
             '| `qa_plots/` | QA pages per <run>/<sub_run>, all stations '
             '(`render_from_store.py`) |',
             '| `tables/` | the full `vmm_subruns.csv` / `vmm_captures.csv` '
             '(mpgd26_figs/aggregate.py), incl. capture-level rows |']
    for chb in sorted(split):
        lines.append(f'| `{chb}/` | rows of ' + ', '.join(
            f'`{k}` ({n})' for k, (n, _) in split[chb].items())
            + f' for this chamber; station.txt |')
    lines += ['| `combined_dream_vmm/` | DREAM vs VMM on the same beam data: '
              '`dream_vmm_comparison/` (figures, data, HTML reports) and '
              '`vmm_dream_matching/` (tables, notes, arrays); copies of the repo folders '
              'of the same name, where the code writes |', '',
              'Station -> chamber (all VMM runs, 29 Jul - 3 Aug): P2_IN = det5, '
              'P2_MID = det1, P2_OUT = det3.', '']
    with open(os.path.join(d, 'INDEX.md'), 'w') as f:
        f.write('\n'.join(lines))
    print(f'wrote {os.path.join(d, "INDEX.md")}')


if __name__ == '__main__':
    main()
