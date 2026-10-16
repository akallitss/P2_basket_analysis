#!/usr/bin/env python3
"""
Mirror the DREAM analysis of the July 2026 SPS test from EOS into a local
tree organised like the cosmic bench: per DETECTOR / run / sub_run / stage.

EOS (`p2_sps_july/analysis/`, the canonical copy, written by condor) is
organised per telescope STATION. Which chamber sat in which station changed
during the test (P2_IN: det2 -> det4 -> det2 -> det5; P2_MID = det1 and
P2_OUT = det3 throughout), so every station folder is renamed run by run with
chamber_history.chamber(). Each run folder gets a station.txt saying where the
chamber was and on what the assignment rests.

Local layout (default root: the LaCie TB_July2026_H4/analysis/):
  det1/ det3/ det2/ det4/ det5/  <run>/<sub_run or scan>/<NN_stage>/   station products
  P2_IN_empty/<run>/             P2_IN products of runs with no chamber powered
  telescope/<run>/...            telescope-wide (event sync, alignment)
  combined/<name>/               several chambers at once: low_field_timing,
                                 urw_timebins, pillars, selftrack (+ the LaCie-only
                                 urw_referenced_efficiency and urw_alignment)
  online_qa/<run>/<sub_run>/<detector>/   the DAQ's online QA plots
  logbook/                       campaign logbook notes
  INDEX.md                       generated: what is where, station <-> chamber
VMM results (EOS analysis/vmm/) are NOT mirrored here: they belong in the VMM
tree (TB_July26_H4_vmm/analysis/).

    python3 sps_beam_analysis/sync_dream_analysis.py [--dest DIR] [--dry-run]
            [--only det1,telescope,...]

Needs an open lxplus connection (ssh alias `lxplus`). rsync only moves what
changed, so re-running it brings the tree up to date.
"""
import argparse
import datetime as dt
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import chamber_history as ch  # noqa: E402

HOST = 'lxplus'
EOS = '/eos/experiment/ntof/data/x17/p2_sps_july/analysis'
DEST = ('/media/ak271430/LaCie/Extras/Physics/Post-Doc-Saclay/data/'
        'SPS_beam_test/TB_July2026_H4/analysis')
STATIONS = ('P2_IN', 'P2_MID', 'P2_OUT')
COMBINED = ('low_field_timing', 'urw_timebins', 'pillars', 'selftrack')
SKIP = ('vmm',)                       # -> VMM tree
LOCAL_ONLY = ('urw_referenced_efficiency', 'urw_alignment')   # in combined/, not on EOS
# Runs before the beam: chamber_history's "MID = det1, OUT = det3 throughout"
# does not hold (run_1, 18 Jul, was the det2/det3 Fe55 test while det1 was on
# the cosmic bench), so they keep their station names under pre_beam/<run>/.
PRE_BEAM = ('run_1',)


def ssh_ls(path):
    out = subprocess.run(['ssh', HOST, f'ls -1p {path}'], capture_output=True,
                         text=True, check=True).stdout.split()
    return [x for x in out if '.sys.v#' not in x]


def rsync(src, dst, dry):
    cmd = ['rsync', '-a', '--exclude=*.sys.v#*', f'{HOST}:{src}', dst]
    if dry:
        print('   ', ' '.join(cmd))
        return
    os.makedirs(dst, exist_ok=True)
    subprocess.run(cmd, check=True)


def station_note(station, run, det):
    lines = [f'station: {station}',
             f'chamber: {det}  (sps_beam_analysis/chamber_history.py, {dt.date.today()})',
             f'source:  EOS {EOS}/{station}/{run}/']
    if run in PRE_BEAM:
        lines[1] = 'chamber: not established (pre-beam run; see PRE_BEAM in sync_dream_analysis.py)'
    if station == 'P2_IN':
        if det == ch.EMPTY:
            lines.append('note:    no chamber powered at P2_IN during this run')
        if run == 'mesh_drift_scan_up_1':
            lines.append('note:    UNCONFIRMED -- det4 by the clock; flips to det2 if the '
                         'swap happened on Sunday afternoon (see chamber_history.py)')
    return '\n'.join(lines) + '\n'


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[1])
    ap.add_argument('--dest', default=DEST)
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--only', default='',
                    help='comma-separated top-level targets (det1, telescope, combined, ...)')
    args = ap.parse_args()
    only = set(filter(None, args.only.split(',')))

    def wanted(target):
        return not only or target in only

    top = ssh_ls(EOS)
    placed = {}                       # target -> [(run, station)]
    for entry in top:
        name = entry.rstrip('/')
        if not entry.endswith('/'):
            if name.endswith('.md') and wanted('logbook'):
                rsync(f'{EOS}/{name}', os.path.join(args.dest, 'logbook') + '/', args.dry_run)
            continue
        if name in SKIP:
            continue
        if name in STATIONS:
            for r in ssh_ls(f'{EOS}/{name}'):
                run = r.rstrip('/')
                det = ch.chamber(name, run)
                target = ('P2_IN_empty' if det == ch.EMPTY
                          else det if det != '?' else f'{name}_unknown')
                if run in PRE_BEAM:           # station -> chamber not established
                    det, target = '?', 'pre_beam'
                if not wanted(target):
                    continue
                dst = (os.path.join(args.dest, target, run, name) if run in PRE_BEAM
                       else os.path.join(args.dest, target, run))
                print(f'{name}/{run} -> {target}/{run}')
                rsync(f'{EOS}/{name}/{run}/', dst + '/', args.dry_run)
                if not args.dry_run:
                    with open(os.path.join(dst, 'station.txt'), 'w') as f:
                        f.write(station_note(name, run, det))
                placed.setdefault(target, []).append((run, name))
        elif name == 'telescope':
            if wanted('telescope'):
                print('telescope -> telescope')
                rsync(f'{EOS}/telescope/', os.path.join(args.dest, 'telescope') + '/', args.dry_run)
        elif name in COMBINED:
            if wanted('combined'):
                print(f'{name} -> combined/{name}')
                rsync(f'{EOS}/{name}/', os.path.join(args.dest, 'combined', name) + '/', args.dry_run)
        else:                         # run-level folders: the DAQ's online QA
            if wanted('online_qa'):
                print(f'{name} -> online_qa/{name}')
                rsync(f'{EOS}/{name}/', os.path.join(args.dest, 'online_qa', name) + '/', args.dry_run)

    if args.dry_run:
        return
    write_index(args.dest, placed)


def write_index(dest, placed):
    runs = sorted({r for v in placed.values() for r, _ in v})
    lines = ['# July 2026 SPS test -- DREAM analysis, per detector', '',
             f'Mirrored from EOS `{EOS}/` by '
             '`sps_beam_analysis/sync_dream_analysis.py`, last run '
             f'{dt.datetime.now():%Y-%m-%d %H:%M}. EOS is the canonical copy (condor '
             'writes there, per station); this tree renames each station to the chamber '
             'that sat in it, run by run (`chamber_history.py`). Each run folder has a '
             '`station.txt`.', '',
             '## Station -> chamber, per run', '',
             '| run | P2_IN | P2_MID | P2_OUT |', '|---|---|---|---|']
    for r in runs:
        if r in PRE_BEAM:
            lines.append(f'| {r} | pre-beam: not established, see `pre_beam/{r}/` | | |')
            continue
        lines.append(f'| {r} | ' + ' | '.join(str(ch.chamber(s, r)) for s in STATIONS) + ' |')
    lines += ['', 'P2_IN `EMPTY` = no chamber powered (products under `P2_IN_empty/`). '
              '`mesh_drift_scan_up_1` P2_IN = det4 is unconfirmed (see its station.txt).',
              '', '## Folders', '']
    for d in sorted(os.listdir(dest)):
        p = os.path.join(dest, d)
        if not os.path.isdir(p):
            continue
        sub = sorted(os.listdir(p))
        what = {'telescope': 'telescope-wide products (event sync, alignment)',
                'combined': 'analyses using several chambers at once',
                'online_qa': 'the DAQ\'s online QA plots, per run / sub_run / detector',
                'logbook': 'campaign logbook notes',
                'P2_IN_empty': 'P2_IN products of runs with no chamber powered',
                'pre_beam': 'pre-beam runs, kept per station (chamber not established)',
                '_archive_pre_condor_2026-07': 'the July station-based copy, before the '
                'condor reprocessing (superseded by EOS; kept, not maintained)'}.get(
                    d, 'chamber products: <run>/<sub_run or scan>/<NN_stage>/')
        lines.append(f'- `{d}/` -- {what}. ' + (f'{len(sub)} entries.' if sub else ''))
    lines += ['', 'Not here: VMM results (`TB_July26_H4_vmm/analysis/`). The LaCie-only '
              f'{", ".join(LOCAL_ONLY)} are in `combined/` (not on EOS).', '']
    with open(os.path.join(dest, 'INDEX.md'), 'w') as f:
        f.write('\n'.join(lines))
    print(f'wrote {os.path.join(dest, "INDEX.md")}')


if __name__ == '__main__':
    main()
