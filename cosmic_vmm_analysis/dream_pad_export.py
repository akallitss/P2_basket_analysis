#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Per-pad DREAM cosmic-bench efficiency in the board (Gerber pad) frame, for comparison with the VMM stand (22 Nov 2026).

Reads the stage-06 ray list of a cosmic_bench_analysis run key (spark-vetoed, default chi2 cut), keeps the rays inside
the active area, and gives every ray the pad its M3 prediction is nearest to (pad_idx = row of the mapped, de-duplicated
channel table that stage 06 built; rebuilt here with the same arguments). Output, one row per pad:
  channel_id, connector_N (1-based, = the VMM run_config's c<N>), pad_cx, pad_cy (board frame, same as the VMM Mapper),
  n (rays), within (hit within the run's MATCH_R of the prediction), any (any P2 pad fired), run_key, dead_connectors.

Runs locally (needs the Analysis tree of cosmic_bench_analysis/p2_qa_config.py and the LaCie run_config):
  dream_pad_export.py OUT_DIR det1_initial1 det1_long5 det2_long1 det3_initial1
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'cosmic_bench_analysis'))
import p2_qa_config as qa      # noqa: E402
import p2_mapping as pmap      # noqa: E402


def export(key, out):
    cfg = qa.get_config(key)
    ct = pmap.build_channel_table(cfg.run_config_path, cfg.MAP_CSV_PATH, det_type=cfg.DET_TYPE,
                                  det_name=cfg.DET_NAME, strategy='reverse', drop_connectors=cfg.DEAD_CONNECTORS,
                                  strategy_overrides=cfg.STRATEGY_OVERRIDES)
    padc = ct[ct['mapped']].drop_duplicates('channel_id').reset_index(drop=True)
    od = cfg.out_dir('06_efficiency')
    sfx = cfg.product_suffix(True) + qa.chi2_tag(qa.M3_CHI2_CUT)
    fp = pd.read_csv(f'{od}/pad_footprint{sfx}.csv')
    if len(fp) != len(padc):
        raise SystemExit(f'{key}: pad table ({len(padc)}) does not match the stage-06 footprint ({len(fp)})')
    d = pd.read_csv(f'{od}/ray_hit_miss_list{sfx}.csv')
    d = d[d.in_active]
    g = d.groupby('pad_idx').agg(n=('within', 'size'), within=('within', 'mean'), any=('has_any', 'mean'))
    r = pd.concat([padc.loc[g.index, ['channel_id', 'connector_N', 'pad_cx', 'pad_cy']].reset_index(drop=True),
                   g.reset_index(drop=True)], axis=1)
    r['run_key'] = key
    r['dead_connectors'] = ' '.join(map(str, cfg.DEAD_CONNECTORS))
    r.to_csv(f'{out}/dream_pads_{key}.csv', index=False)
    print(f'{key}: {len(r)} pads, {int(r.n.sum())} rays, within {100 * (r.within * r.n).sum() / r.n.sum():.1f} %')


if __name__ == '__main__':
    os.makedirs(sys.argv[1], exist_ok=True)
    for k in sys.argv[2:]:
        export(k, sys.argv[1])
