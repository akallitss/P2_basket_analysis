#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
vmm_trigger_timing.py -- gather the trigger-referenced VMM timing and
efficiency of the whole campaign into one table.

Nothing here re-analyses anything: `vmm_efficiency.py` (P2_basket_online_analysis)
already ran the measurement on every capture during the beam test, and each
capture's `scalars.json` carries the result under `efficiency`:

    mu_ns, sigma_ns   Gaussian fit of the station-minus-trigger BCID phase
                      difference, hits paired only within one SRS marker
    contrast          peak height over the flat accidental background
    raw / accidental / corrected efficiency, and the binomial interval

That is the measurement made at the beam: efficiency relative to the TRIGGER
acceptance, accidental-subtracted from an equal-width sideband, with the
coincidence built on BCID phase inside a single marker interval because the
firmware leaves `offset` stuck and absolute time is ambiguous by multiples of
92.16 us.  See vmm_efficiency.py's header for why that grouping is necessary.

This script joins those per-capture results to each run's own run_config.json
(gas, mesh, drift, chip config) so the campaign can be sliced by working point
and by gas.

Products:  vmm_trigger_timing.csv       one row per capture per station
           vmm_trigger_timing_by_hv.csv median sigma/mu per (station, gas, HV)

Usage:  python3 vmm_trigger_timing.py [--analysis DIR] [--meta DIR]
"""

import argparse
import glob
import json
import os

import numpy as np
import pandas as pd

ANALYSIS = ('/local/home/ak271430/Documents/PostDocSaclay/data/SPS_Beam_Test/'
            'mpgd26_workspace/products/analysis/vmm')
META = ('/local/home/ak271430/Documents/PostDocSaclay/data/SPS_Beam_Test/'
        'mpgd26_workspace/eos_inventory/vmm_meta/runs')
STATIONS = ('P2_IN', 'P2_MID', 'P2_OUT')


def conditions(meta):
    """(run, sub_run) -> gas, start, chip config and per-station mesh/drift."""
    out = {}
    for run in sorted(os.listdir(meta)):
        p = os.path.join(meta, run, 'run_config.json')
        if not os.path.exists(p):
            continue
        c = json.load(open(p))
        ch = {d['name']: d.get('hv_channels', {}) for d in c.get('detectors', [])}
        for sub in c.get('sub_runs', []):
            hvs = sub.get('hvs', {})

            def hv(station, kind):
                pair = ch.get(station, {}).get(kind)
                if not pair:
                    return None
                return hvs.get(str(pair[0]), {}).get(str(pair[1]))

            row = dict(gas=c.get('gas'), start=c.get('start_time'),
                       chip=c.get('chip_config'))
            for st in STATIONS:
                row[f'mesh_{st}'] = hv(st, 'mesh')
                row[f'drift_{st}'] = hv(st, 'drift')
            out[(run, sub['sub_run_name'])] = row
    return out


def gather(analysis, cond):
    rows = []
    for sc in glob.glob(os.path.join(analysis, '*', '*', '*', 'scalars.json')):
        parts = sc.split(os.sep)
        run, sub, cap = parts[-4], parts[-3], parts[-2]
        try:
            d = json.load(open(sc))
        except (ValueError, OSError):
            continue
        eff = d.get('efficiency') or {}
        cnd = cond.get((run, sub), {})
        for st in STATIONS:
            v = eff.get(st)
            if not isinstance(v, dict):
                continue
            rows.append(dict(
                run=run, sub=sub, capture=cap, station=st,
                gas=cnd.get('gas'), start=cnd.get('start'),
                chip=cnd.get('chip'), mesh=cnd.get(f'mesh_{st}'),
                drift=cnd.get(f'drift_{st}'),
                mu_ns=v.get('mu_ns'), sigma_ns=v.get('sigma_ns'),
                contrast=v.get('contrast'), n_triggers=v.get('n_triggers'),
                eff=v.get('efficiency'), eff_corr=v.get('efficiency_corrected'),
                accidental=v.get('accidental_efficiency')))
    return pd.DataFrame(rows)


def quality(t):
    """Drop captures whose fit did not converge.

    A failed fit is unmistakable: mu runs off to the edge of the +-1 us search
    window and the peak sits at the accidental floor.  Keeping them would drag
    every median.
    """
    return t[(t.contrast > 5) & np.isfinite(t.contrast)
             & (t.mu_ns.abs() < 600) & (t.sigma_ns > 2) & (t.sigma_ns < 200)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--analysis', default=ANALYSIS)
    ap.add_argument('--meta', default=META)
    a = ap.parse_args()

    t = gather(a.analysis, conditions(a.meta))
    print(f'{len(t)} station-captures over {t.run.nunique()} runs')
    q = quality(t)
    print(f'{len(q)} survive the fit-quality cut '
          f'({100 * len(q) / max(len(t), 1):.0f} %)')
    q.to_csv('vmm_trigger_timing.csv', index=False)

    by = (q.groupby(['station', 'gas', 'mesh', 'drift'])
            .agg(n_captures=('sigma_ns', 'size'),
                 sigma_ns=('sigma_ns', 'median'),
                 mu_ns=('mu_ns', 'median'),
                 contrast=('contrast', 'median'),
                 eff=('eff', 'median'))
            .round(2).reset_index())
    by.to_csv('vmm_trigger_timing_by_hv.csv', index=False)
    print('\n=== sigma at mesh 450, by drift and gas ===')
    n = by[by.mesh == 450]
    for st in STATIONS:
        s = n[n.station == st]
        if not len(s):
            continue
        print(f'\n{st}')
        print(s.pivot_table(index='drift', columns='gas',
                            values='sigma_ns').round(1).to_string())


if __name__ == '__main__':
    main()
