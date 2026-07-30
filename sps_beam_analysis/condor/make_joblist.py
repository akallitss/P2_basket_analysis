#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_joblist.py -- enumerate what is on EOS and write condor's joblist.txt.

    ./make_joblist.py [--group hits|wave] [--run RUN[,RUN...]] [-o joblist.txt]

One line per job: "<run> <sub_run> <group>", which analysis.sub consumes via
`queue run,sub_run,group from joblist.txt`.

A sub_run is only listed if the input the group actually needs is present on
EOS -- combined_hits_root/*.root for `hits`, decoded_root/*.root for `wave`.
That matters because the campaign's EOS mirror is uneven: several runs were
pushed across three different destinations while quotas were being sorted out,
and a sub_run whose hits never made it across would otherwise become a job
that starts, transfers nothing and fails.

The whole EOS tree is listed in ONE recursive xrdfs call rather than one per
sub_run: at ~170 sub_runs the per-call Kerberos handshake dominates everything
else (measured ~1 s each, i.e. minutes, vs a few seconds for the sweep).
"""

import argparse
import os
import subprocess
import sys
from collections import defaultdict

# The campaign is spread over three EOS locations: it moved twice in three days
# while quotas were being sorted out, so the early runs exist ONLY on salsachip.
# Keys here must match the `case "$SRC"` block in job.sh.
SOURCES = {
    'ntof': ('root://eospublic.cern.ch',
             '/eos/experiment/ntof/data/x17/p2_sps_july'),
    'salsachip': ('root://eosproject.cern.ch',
                  '/eos/project/s/salsachip/Data/T2_tests/P2_SPS_Dream_Data'),
    'user': ('root://eosuser.cern.ch',
             '/eos/user/a/akallits/P2_SPS_backup_temp'),
}

# What each stage group needs to exist on EOS before a job is worth queueing.
# 'rec' reads decoded_root's eventId branch remotely rather than staging it.
NEED = {'hits': 'combined_hits_root',
        'raweff': 'decoded_root',
        'wave': 'decoded_root',
        'rec': 'decoded_root'}

# A sub_run taken with no beam still produces .root files -- just empty ones,
# a ROOT header and nothing else. They are indistinguishable from real input by
# existence alone, so queue on total input SIZE instead. The two populations are
# cleanly separated on this campaign: the largest empty sub_run is 9.5 kB, the
# smallest real one 888 kB (drift_mesh_2d_1/dm_01_01_m430_d490, a genuine point
# from the run the P2_MID trip cut short). 100 kB sits ~10x above the former and
# ~9x below the latter. Do NOT raise this to 1 MB: that swallows the small but
# real drift_mesh_2d_1 points, which have already been analysed successfully.
#
# The DAQ also drops a NEEDS_RETAKE marker beside no-beam points, but that flag
# is not a usable skip signal on its own -- eff_nominal_1/eff_nominal_12 carries
# it and still holds 227 MB of perfectly good data.
MIN_INPUT_BYTES = 100 * 1024


def xrdfs_ls_recursive(url, base):
    """Every path under `base` with its size, as (path, size), one xrdfs call.

    `-l` as well as `-R`, because "the file exists" is not the same question as
    "the file has data in it". A sub_run whose only .root file is a zero-byte
    husk passes an existence check and then fails the job. Seen on
    p2_mesh_drift_2d_1/g450_m440, where the beam died part way through and left
    an empty file behind while its neighbours are ~500 MB.
    """
    # -l and -R must be SEPARATE flags: 'ls -lR' silently ignores the -l and
    # returns bare paths, which parses as zero usable entries.
    cmd = ['xrdfs', url, 'ls', '-l', '-R', base]
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    if p.returncode != 0:
        sys.exit(f'xrdfs ls -l -R failed ({p.returncode}) on {url}{base}:\n'
                 f'{p.stderr.strip()}')
    out = []
    for ln in p.stdout.splitlines():
        f = ln.split()
        # perms owner group size date time path
        if len(f) < 5:
            continue
        try:
            size = int(f[3])
        except ValueError:
            continue
        out.append((f[-1], size))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--group', default='hits', choices=sorted(NEED),
                    help='stage group to queue (default: hits)')
    ap.add_argument('--run', default=None,
                    help='comma-separated run names (default: every run found)')
    ap.add_argument('-o', '--out', default='joblist.txt')
    ap.add_argument('--print-runs', action='store_true',
                    help='just summarise what is on EOS and exit')
    ap.add_argument('--source', default='ntof',
                    help='comma-separated EOS sources to sweep: '
                         + ', '.join(sorted(SOURCES)) + ", or 'all'. A run "
                         'present on more than one source is taken from the '
                         'FIRST one listed that actually has its data.')
    args = ap.parse_args()

    srcs = (sorted(SOURCES) if args.source == 'all'
            else [s.strip() for s in args.source.split(',') if s.strip()])
    for s in srcs:
        if s not in SOURCES:
            sys.exit(f'unknown source {s!r}; known: {", ".join(sorted(SOURCES))}')

    want = NEED[args.group]
    # runs[run][sub_run] = source key that will supply it
    runs = defaultdict(dict)
    all_subs = defaultdict(set)
    claimed = {}                      # run -> source that won it
    empty = []                        # (run, sub_run, bytes) below MIN_INPUT_BYTES

    for src in srcs:
        url, base = SOURCES[src]
        runs_base = f'{base}/runs'
        print(f'listing {url}{runs_base} ...', flush=True)
        paths = xrdfs_ls_recursive(url, runs_base)
        print(f'  {len(paths)} paths')
        prefix = runs_base.rstrip('/') + '/'
        found = defaultdict(set)
        nbytes = defaultdict(int)         # (run, sub_run) -> input bytes
        for p, size in paths:
            if not p.startswith(prefix):
                continue
            parts = p[len(prefix):].split('/')
            # <run>/<sub_run>/<subdir>/<file>.root
            # Skip EOS atomic-version placeholders (.sys.v#.<name>.root):
            # they match *.root but are not copyable, so counting them would
            # queue jobs for sub_runs that have no real data -- and would
            # inflate the totals reported here.
            if (p.endswith('.root') and len(parts) == 4 and parts[2] == want
                    and not os.path.basename(p).startswith('.')):
                nbytes[(parts[0], parts[1])] += size
            # every sub_run dir that exists, so gaps can be reported rather
            # than silently dropped
            if len(parts) >= 2 and '.' not in parts[1]:
                all_subs[parts[0]].add(parts[1])
        # Size is judged on the sub_run TOTAL, not per file: a no-beam point can
        # hold several .root files that are each a bare header.
        for (run, sub), n in nbytes.items():
            if n >= MIN_INPUT_BYTES:
                found[run].add(sub)
            else:
                empty.append((run, sub, n))
        for run, subs in found.items():
            # First source listed that actually holds a run's data wins it.
            # This is what makes `--source ntof,salsachip` do the right thing:
            # take everything from nTOF, and fall back to salsachip only for
            # the early runs that never made it across.
            if run in claimed:
                continue
            claimed[run] = src
            for s in subs:
                runs[run][s] = src

    selected = sorted(runs) if not args.run else \
        [r.strip() for r in args.run.split(',') if r.strip()]

    lines, n_skip = [], 0
    print(f'\n{"run":<26} {"src":<10} {"jobs":>5} {"skipped":>8}   '
          f'(group={args.group}, needs {want})')
    print('-' * 72)
    for run in selected:
        subs = sorted(runs.get(run, {}))
        src = claimed.get(run, '-')
        missing = sorted(all_subs.get(run, set()) - set(subs))
        n_skip += len(missing)
        print(f'{run:<26} {src:<10} {len(subs):>5} {len(missing):>8}'
              + (f'   -> {", ".join(missing)}' if missing else ''))
        # prod_sub: where stages 20/22 file their per-point products. A worker
        # only ever sees ONE sub_run, so it cannot tell a 49-point scan from a
        # single-point run -- but that decides whether the eff maps belong
        # together under 'scan/' (where a serial pass would have put them, and
        # where the GIF builder and the GUI look) or under the sub_run's own
        # directory. Decided here, where the whole run is in view.
        prod_sub = 'scan' if len(subs) > 1 else (subs[0] if subs else '-')
        if args.group == 'raweff':
            # RUN-level: stage 30 sweeps every sub_run itself and writes one
            # curve per run, so one job covers the lot. The sub_run column is
            # still filled (job.sh logs and names by it) but goes unused.
            if subs:
                lines.append(f'{run} {subs[0]} {args.group} {prod_sub} '
                             f'{runs[run][subs[0]]}')
            continue
        for s in subs:
            lines.append(f'{run} {s} {args.group} {prod_sub} {runs[run][s]}')
    print('-' * 72)
    print(f'{"TOTAL":<26} {"":<10} {len(lines):>5} {n_skip:>8}')

    if args.print_runs:
        return
    with open(args.out, 'w') as fh:
        fh.write('\n'.join(lines) + ('\n' if lines else ''))
    print(f'\nwrote {args.out} ({len(lines)} jobs)')
    if n_skip:
        print(f'NOTE: {n_skip} sub_run(s) on EOS have no {want} and were '
              f'skipped — check whether their backup is incomplete.')
    if empty:
        print(f'NOTE: {len(empty)} sub_run(s) have a {want} holding less than '
              f'{MIN_INPUT_BYTES // 1024} kB — taken with no beam, not queued:')
        # dedupe: a run present on two EOS sources is listed once per source
        for run, sub, n in sorted(set(empty)):
            print(f'        {run}/{sub}  ({n} B)')


if __name__ == '__main__':
    main()
