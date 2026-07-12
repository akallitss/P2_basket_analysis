#!/usr/bin/env python3
"""In-window vs out-of-window pad hitmap for det1_long3 (7-7-26 run)."""
import os, sys, glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib.collections import PolyCollection

sys.path.insert(0, '/local/home/ak271430/Documents/PostDocSaclay/P2_basket_analysis/cosmic_bench_analysis')
import p2_qa_config as qa
import p2_mapping as pmap
import uproot

FIGS = os.path.join(qa.DATA_ROOT, 'Analysis', 'reports',
                    'p2_efficiency_intermittency', 'figs')

cfg = qa.get_config('det1_long3')
ct = pmap.build_channel_table(cfg.run_config_path, cfg.MAP_CSV_PATH,
                              det_type=cfg.DET_TYPE, det_name=cfg.DET_NAME,
                              strategy='reverse',
                              drop_connectors=cfg.DEAD_CONNECTORS)
feu_set = set(ct.attrs['feus'])
parts = []
for fp in sorted(glob.glob(os.path.join(cfg.combined_hits_dir, '*.root'))):
    a = uproot.open(f'{fp}:hits').arrays(
        ['eventId', 'trigger_timestamp_ns', 'channel', 'amplitude', 'feu'],
        library='pd')
    parts.append(a[a['feu'].isin(feu_set)])
df = pd.concat(parts, ignore_index=True)
df = pmap.attach_pads_to_hits(df, ct)
df = df[df['mapped'].fillna(False).astype(bool)].copy()
npads = df.groupby('eventId')['channel'].transform('size')
df = df[npads < cfg.BURST_NPADS]
df['t_h'] = df['trigger_timestamp_ns'] / 1e9 / 3600.0

lo, hi = 4.0, 6.25
sel_in = (df['t_h'] >= lo) & (df['t_h'] <= hi)
pads, verts = pmap.pad_tiles(ct)

fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.4))
for ax, sub, ttl in [
        (axes[0], df[sel_in], f'inside window ({hi-lo:.2f} h, '
                              f'{int(sel_in.sum()):,} hits)'),
        (axes[1], df[~sel_in], f'outside window ({10.0-(hi-lo):.2f} h, '
                               f'{int((~sel_in).sum()):,} hits)')]:
    counts = sub.groupby('channel_id').size().reindex(
        pads['channel_id']).fillna(0).to_numpy()
    fired = counts > 0
    dead = PolyCollection(verts[~fired], facecolors='0.92',
                          edgecolors='0.7', linewidths=0.3)
    ax.add_collection(dead)
    pc = PolyCollection(verts[fired], array=counts[fired], cmap='inferno',
                        norm=matplotlib.colors.LogNorm(vmin=1),
                        edgecolors='face', linewidths=0.2)
    ax.add_collection(pc)
    fig.colorbar(pc, ax=ax, label='pad hits (log)')
    ax.set_xlabel('pad_cx [mm]'); ax.set_ylabel('pad_cy [mm]')
    ax.autoscale_view(); ax.set_aspect('equal')
    ax.set_title(ttl, fontsize=10)
fig.suptitle('P2_1  7-7-26 — pad hitmap inside vs outside the efficiency '
             'window (grey = never fired)', y=1.0)
fig.tight_layout()
fig.savefig(os.path.join(FIGS, 'hitmap_window_det1_long3.pdf'),
            bbox_inches='tight')
fig.savefig(os.path.join(FIGS, 'hitmap_window_det1_long3.png'), dpi=140,
            bbox_inches='tight')
n_in, n_out = int(sel_in.sum()), int((~sel_in).sum())
p_in = df[sel_in]['channel_id'].nunique()
p_out = df[~sel_in]['channel_id'].nunique()
print(f'pads fired inside window: {p_in}/{len(pads)}, outside: {p_out}/{len(pads)}')
