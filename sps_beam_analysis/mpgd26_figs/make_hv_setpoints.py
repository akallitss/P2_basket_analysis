#!/usr/bin/env python3
"""Per (run, sub_run, detector) HV setpoints, from the DAQ run configs.

`run_config.json` records the ladder the operator programmed, one `hvs` block
per sub_run keyed by CAEN card then channel.  The stage products only carry the
sub_run *name*, so every figure that wants a voltage axis (the 2D scans above
all) joins against this table.

Card 8 carries the whole telescope:
    8:0/8:1 = P2_IN drift/mesh      8:6 = uRW_front drift
    8:2/8:3 = P2_MID drift/mesh     8:7 = uRW_back drift
    8:4/8:5 = P2_OUT drift/mesh     12:0, 12:1 = uRWELL resistive
"""

import glob
import json
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import S  # noqa: E402

# detector -> (drift channel, mesh-or-resistive channel) on the named card
CHANNELS = {
    'P2_IN': ('8', '0', '1'),
    'P2_MID': ('8', '2', '3'),
    'P2_OUT': ('8', '4', '5'),
    'uRW_front': ('8', '6', None),
    'uRW_back': ('8', '7', None),
}
RESIST = {'uRW_front': ('12', '0'), 'uRW_back': ('12', '1')}


def main():
    rows = []
    for f in sorted(glob.glob(f'{S}/eos_inventory/runs/*/run_config.json')):
        run = os.path.basename(os.path.dirname(f))
        try:
            cfg = json.load(open(f))
        except (OSError, ValueError):
            continue
        for sr in cfg.get('sub_runs', []):
            hvs = sr.get('hvs') or {}
            for det, (card, dch, mch) in CHANNELS.items():
                block = hvs.get(card) or {}
                drift = block.get(dch)
                mesh = block.get(mch) if mch else None
                if mesh is None and det in RESIST:
                    rcard, rch = RESIST[det]
                    mesh = (hvs.get(rcard) or {}).get(rch)
                if drift is None and mesh is None:
                    continue
                rows.append(dict(run=run, sub_run=sr.get('sub_run_name'),
                                 det=det, drift=drift, mesh_or_resist=mesh))
    if not rows:
        sys.exit(f'no run configs under {S}/eos_inventory/runs -- '
                 'run fetch_from_eos.sh first')
    out = f'{S}/eos_inventory/hv_setpoints.csv'
    df = pd.DataFrame(rows)
    df.to_csv(out, index=False)
    print(f'hv_setpoints.csv     {len(df)} rows, runs={df.run.nunique()}')


if __name__ == '__main__':
    main()
