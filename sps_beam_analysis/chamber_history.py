#!/usr/bin/env python3
"""Which physical chamber sat in which telescope station, and when.

Source of truth: the operator's account of the July 2026 H4 run (restated
2026-08-19), cross-checked run by run against the HV monitor -- the mesh
setpoint of the P2_IN channel tells you unambiguously whether ANY chamber was
powered there, and the wall-clock start of every run is in
report_data/hv_campaign_30s.csv.  Do NOT use instead:

  * `run_config.json` `description` fields -- they contradict each other
    between runs and were never updated on a swap;
  * `p2_qa_config.py:624` -- claims the 7-18-26 Fe55 run measured
    det2 = P2_OUT and det3 = P2_MID, which the logbook contradicts.

P2_MID = det1 and P2_OUT = det3 for the WHOLE test; neither was ever touched.
Every swap happened at P2_IN, and only there:

    Wed 22 - Thu 23 Jul    P2_IN = det2    dead: the tag-probe efficiency at
                                           mesh 430-450 V is 0.03-0.05, i.e.
                                           no signal (drift HV contact)
    Fri 24 Jul (midday)    P2_IN = det4    det4 has a leaky drift frame; runs
                                           from 12:28 on Fri 24
    Sun 26 Jul (morning)   P2_IN off       low_mesh_scan_1: P2_IN mesh at 0 V
                                           for the whole run, FEU 3 out
    Sun 26 Jul ~23:00      P2_IN = det2    det2 repaired and swapped back in,
                                           late Sunday evening
    Mon 27 Jul (midday)    P2_IN pulled    from 15:42 the P2_IN mesh is at 0 V
                                           in every run -- the telescope ran on
                                           MID + OUT alone for the afternoon
    Mon 27 Jul (night)     det5 installed  the CERN-built chamber goes in, but
                                           is NOT on HV: connection issues
    Tue 28 Jul ~10:00      P2_IN = det5    first HV, first data; in place for
                                           the rest of the campaign

`det5` is this module's name for the CERN-built chamber, so that every station
can be quoted with a detector number.  Its drift frame needs modification for
the old HV connection.

UNCONFIRMED, flagged 2026-08-19: `mesh_drift_scan_up_1` (Sun 26 Jul,
15:54-20:50).  The P2_IN mesh ramps 60 -> 450 V, so a chamber WAS powered --
but the run ends more than two hours before the ~23:00 det2 swap, so by the
timeline it must still be det4.  Its efficiency curve sits about 5 points
above BOTH the known det4 curve (25 Jul) and the known det2 curve (27 Jul) at
equal mesh voltage, so the curve shape does not settle it either way; the
drift gap differs between the runs.  Assigned det4 here on the strength of the
clock.  If the swap actually happened Sunday afternoon, this one run -- and
only this one -- flips to det2.
"""

# station -> chamber, for the stations that never changed
STATIC = {'P2_MID': 'det1', 'P2_OUT': 'det3'}

# sentinel: a chamber is installed or not, but nothing is powered / read out
EMPTY = 'EMPTY'

# run -> chamber sitting at P2_IN while that run was taken.  Wall-clock start
# from the HV monitor is quoted for every entry, because that is what the
# assignment is actually based on.
P2IN_BY_RUN = {
    'run_1':                   'det2',   # Sat 18 Jul 20:45 (pre-beam, Fe55)
    'beam_commissioning_1':    'det2',   # Thu 23 Jul 17:06  -- dead, eff 0.03
    'latency_scan_1':          'det2',   # Thu 23 Jul 17:19
    'beam_nominal_meshscan_1': 'det2',   # Thu 23 Jul 18:17  -- dead, eff 0.05
    'drift_scan_1':            'det4',   # Fri 24 Jul 12:28  -- det4 goes in
    'drift_scan_2':            'det4',   # Fri 24 Jul 14:37
    'meshscan_fine_1':         'det4',   # Fri 24 Jul 18:38
    'p2in_check_1':            'det4',   # Fri 24 Jul 20:49  -- the P2_IN debug
    'drift_scan_final':        'det4',   # Fri 24 Jul 23:43
    'drift_mesh_scan_1':       'det4',   # Sat 25 Jul 01:25  -- mesh 38-430 V
    'env_test_1':              'det4',   # Sat 25 Jul 10:54  -- NOT 22 Jul
    'highstat_eff_1':          'det4',   # Sat 25 Jul 11:37
    'drift_mesh_2d_1':         'det4',   # Sat 25 Jul 19:28
    'low_mesh_scan_1':          EMPTY,   # Sun 26 Jul 11:47  -- mesh 0 V, FEU3 out
    'mesh_drift_scan_up_1':    'det4',   # Sun 26 Jul 15:54  -- see UNCONFIRMED
    'drift_mesh_2d_2':         'det2',   # Mon 27 Jul 01:49  -- after the swap
    'eff_nominal_1':           'det2',   # Mon 27 Jul 10:40  -- ends 13:48
    'eff_drift_ab_1':           EMPTY,   # Mon 27 Jul 15:42  -- mesh 0 V
    'drift_transparency_1':     EMPTY,   # Mon 27 Jul 20:05  -- mesh 0 V
    'p2in_hvrange_1':          'det5',   # Tue 28 Jul 11:03  -- CERN, first HV
    'p2in_hvrange_2':          'det5',   # Tue 28 Jul 11:43
    'p2_mesh_drift_eff_1':     'det5',   # Tue 28 Jul 15:38
    'p2_mesh_drift_2d_1':      'det5',   # Wed 29 Jul 00:07
}

# when each chamber sat at P2_IN, for legends that span the swap
P2IN_WINDOW = {'det2': '22–23 and 26–27 Jul', 'det4': '24–26 Jul',
               'det5': 'from 28 Jul'}

NAME = {'det1': 'det1', 'det2': 'det2', 'det3': 'det3', 'det4': 'det4',
        'det5': 'det5 = CERN-built', EMPTY: 'no chamber powered', None: '?'}


# what P2_IN saw over the whole campaign, for figures that span many runs and
# cannot name a single chamber
P2IN_SEQUENCE = 'det2 → det4 → det2 → det5'


def chamber(station, run):
    """Physical chamber in `station` during `run` ('?' when unrecorded)."""
    if station in STATIC:
        return STATIC[station]
    if station != 'P2_IN':
        return '?'
    if run in P2IN_BY_RUN:
        return P2IN_BY_RUN[run]
    # the numbered runs are the August campaign: the CERN-built chamber went
    # in on 28 Jul and stayed for the rest of the test, so every run_NN from
    # run_21 (30 Jul) on is det5.  run_1 is 18 Jul, before the beam.
    import re as _re
    m = _re.fullmatch(r'run_(\d+)', str(run))
    if m:
        return 'det2' if int(m.group(1)) == 1 else 'det5'
    return '?'


def label_any(station):
    """Label for a figure spanning many runs: names the whole P2_IN sequence."""
    if station in STATIC:
        return f'{station} ({STATIC[station]})'
    if station == 'P2_IN':
        return f'P2_IN ({P2IN_SEQUENCE})'
    return station


def label(station, run, window=False):
    """Legend label: 'P2_MID (det1)', or 'P2_IN (det4, 24-26 Jul)'.

    Station first, chamber in brackets -- the station is the fixed thing on the
    beam line and the chamber is what changed.  `window` adds the dates, which
    is only worth the width when one figure shows P2_IN more than once.
    """
    c = chamber(station, run)
    name = NAME.get(c, c)
    if c == EMPTY:
        return f'{station} ({name})'
    if window and c in P2IN_WINDOW:
        return f'{station} ({name}, {P2IN_WINDOW[c]})'
    return f'{station} ({name})'


def bench_label(det, station=None):
    """The mirror image, for bench series: 'det1 on the bench (= P2_MID)'."""
    if station is None:
        station = {v: k for k, v in STATIC.items()}.get(det)
    return (f'{det} on the bench (= {station})' if station
            else f'{det} on the bench')


def powered(run):
    """Was any chamber powered at P2_IN during `run`?"""
    return chamber('P2_IN', run) not in (EMPTY, '?', None)


if __name__ == '__main__':
    print(f'{"run":<26} {"P2_IN":<22} P2_MID / P2_OUT')
    for r in P2IN_BY_RUN:
        print(f'{r:<26} {NAME.get(chamber("P2_IN", r)):<22} det1 / det3')
