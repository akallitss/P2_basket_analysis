#!/usr/bin/env python3
"""V1/V3/V5 -- read the run conditions off the VMM DAQ's own run_config.json.

The DREAM and VMM DAQs number runs independently and drift apart after run_58,
and 136 gas labels on EOS were wrong once already, so nothing here is inferred
from a run number: gas, HV, chip config and start time all come from the config
the VMM DAQ itself wrote at the start of the run.

Usage: python3 verify_configs.py [run ...]
"""
import json
import os
import sys

META = ('/local/home/ak271430/Documents/PostDocSaclay/data/SPS_Beam_Test/'
        'mpgd26_workspace/eos_inventory/vmm_meta/runs')

# HV card 8: channel pairs (drift, mesh) per station, read off the detector
# blocks of the config rather than assumed.
def hv_channels(cfg):
    out = {}
    for det in cfg.get('detectors', []):
        ch = det.get('hv_channels', {})
        if 'mesh' in ch and 'drift' in ch:
            out[det['name']] = {'mesh': ch['mesh'], 'drift': ch['drift']}
    return out


def hv_values(sub, chmap):
    hvs = sub.get('hvs', {})
    out = {}
    for name, ch in chmap.items():
        def get(pair):
            card, chan = str(pair[0]), str(pair[1])
            return hvs.get(card, {}).get(chan)
        out[name] = (get(ch['mesh']), get(ch['drift']))
    return out


def main(runs):
    rows = []
    for r in runs:
        p = os.path.join(META, r, 'run_config.json')
        if not os.path.exists(p):
            print(f'{r}: no config on disk')
            continue
        c = json.load(open(p))
        chmap = hv_channels(c)
        for sub in c.get('sub_runs', []):
            rows.append(dict(
                run=r, sub=sub['sub_run_name'],
                start=c.get('start_time'), gas=c.get('gas'),
                chip=c.get('chip_config'),
                plan=c.get('run_plan'),
                hv=hv_values(sub, chmap)))

    w = max(len(f"{x['run']}/{x['sub']}") for x in rows) if rows else 20
    print(f"{'run/sub_run':<{w}}  {'start':<19}  {'gas':<22}  "
          f"{'mesh/drift IN,MID,OUT':<28}  chip_config")
    print('-' * (w + 100))
    for x in rows:
        hv = ' '.join(f"{(x['hv'].get(s) or (None, None))[0]}/"
                      f"{(x['hv'].get(s) or (None, None))[1]}"
                      for s in ('P2_IN', 'P2_MID', 'P2_OUT'))
        print(f"{x['run']+'/'+x['sub']:<{w}}  {str(x['start']):<19}  "
              f"{str(x['gas']):<22}  {hv:<28}  {x['chip']}")
    return rows


if __name__ == '__main__':
    args = sys.argv[1:] or ['run_38', 'run_39', 'run_40', 'run_41', 'run_44',
                            'run_45', 'run_46', 'run_47', 'run_48', 'run_54']
    main(args)
