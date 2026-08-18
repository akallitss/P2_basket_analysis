#!/usr/bin/env python3
"""Which physical chamber sat in which telescope station, and when.

Source of truth: the hardware logbook (entries of 2026-07-26 and 2026-07-28),
which is the ONLY record that is self-consistent.  Do not use instead:

  * `run_config.json` `description` fields -- they contradict each other
    between runs (`run_1` says P2_OUT = "det2, bulked 25-6-26"; every later run
    says P2_IN = "det2 in bulking order"), and were never updated on a swap;
  * `p2_qa_config.py:624` -- claims the 7-18-26 Fe55 run measured
    det2 = P2_OUT and det3 = P2_MID, which the logbook contradicts.

All four Saclay chambers travelled to H4.  P2_MID and P2_OUT were never touched
for the whole test; every swap happened at P2_IN, and only there.

    Wed 22 - Thu 23 Jul   P2_IN = det2      ~20x low signal on 23 Jul: the
                                            drift HV contact problem
    Fri 24 Jul (morning)  P2_IN = det4      det4 turns out to have a leaky
                                            drift frame
    Sun 26 Jul            P2_IN = det2      back in, drift HV contact believed
                                            fixed
    Mon 27 Jul ~15:00     P2_IN pulled
    Tue 28 Jul (morning)  P2_IN = CERN      the CERN-built chamber, in place
                                            for the rest of the campaign; its
                                            drift frame needs modification for
                                            the old HV connection

Run dates come from `analysis/2026-07-27_run_inventory.md`.
"""

# station -> chamber, for the stations that never changed
STATIC = {'P2_MID': 'det1', 'P2_OUT': 'det3'}

# run -> chamber sitting at P2_IN while that run was taken
P2IN_BY_RUN = {
    'run_1': 'det2',
    'env_test_1': 'det2',                  # ~22 Jul
    'beam_commissioning_1': 'det2',        # 23 Jul
    'latency_scan_1': 'det2',              # 23 Jul
    'beam_nominal_meshscan_1': 'det2/det4',  # 23-24 Jul, spans the swap
    'drift_scan_1': 'det4',                # 24 Jul
    'drift_scan_2': 'det4',                # 24 Jul
    'meshscan_fine_1': 'det4',             # 24 Jul
    'p2in_check_1': 'det4',                # 24 Jul, the P2_IN debug run
    'drift_scan_final': 'det4',            # 24 Jul 23:45
    'drift_mesh_scan_1': 'det4',           # 25 Jul am
    'highstat_eff_1': 'det4',              # 25 Jul pm
    'drift_mesh_2d_1': 'det4',             # 25 Jul eve
    'low_mesh_scan_1': None,               # 26 Jul am: P2_IN at 0 V, FEU 3 out
    'mesh_drift_scan_up_1': 'det2',        # 26 Jul pm
    'drift_mesh_2d_2': 'det2',             # 27 Jul night
    'eff_nominal_1': 'det2',               # 27 Jul 10:40
    'p2in_hvrange_1': 'CERN',              # 28 Jul am
    'p2in_hvrange_2': 'CERN',              # 28 Jul am
    'p2_mesh_drift_eff_1': 'CERN',         # 28 Jul pm
    'p2_mesh_drift_2d_1': 'CERN',
    'eff_drift_ab_1': 'CERN',
    'drift_transparency_1': 'CERN',
}

LABEL = {
    'det1': 'det1', 'det2': 'det2', 'det3': 'det3', 'det4': 'det4',
    'det2/det4': 'det2 → det4', 'CERN': 'CERN-built chamber', None: 'not read out',
}


def chamber(station, run):
    """Physical chamber in `station` during `run` ('?' when unrecorded)."""
    if station in STATIC:
        return STATIC[station]
    if station != 'P2_IN':
        return '?'
    return P2IN_BY_RUN.get(run, '?')


def label(station, run):
    """Human label for a legend entry, e.g. 'P2_IN = det4'."""
    c = chamber(station, run)
    return f'{station} = {LABEL.get(c, c)}'


if __name__ == '__main__':
    for r in P2IN_BY_RUN:
        print(f'{r:<26} P2_IN = {LABEL.get(chamber("P2_IN", r))}')
