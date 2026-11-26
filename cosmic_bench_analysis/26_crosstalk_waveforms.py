"""26_crosstalk_waveforms.py -- electronics crosstalk from waveform shapes.

Usage: python3 26_crosstalk_waveforms.py <run_key> [--max-chunk N]

First used on det2_long_repos (2026-11-26): uniform same-shape coupling of
~0.85% to ch+-1/+-2 (0.16% at +-3, 0 at +-17..19), plus 8 SHORTED groups at
~100% (now cfg.SHORTED_CHANNELS). Needs decoded_root + the 06 ray list.

Source = channel whose pad lies under a clean M3 ray (<= 5 mm from its centre)
and that carries a clear pulse. For each source pulse we look at the SAME
event's waveforms on electronics neighbours ch +-1..+-3 (same FEU) whose pads
are > 30 mm from the ray (so no real charge can reach them), and on far control
channels (ch +-17..+-19, same 64-ch block -> same chip, same CNS block).

Per neighbour: peak (max) and trough (min) relative to the source peak, peak
time offset, correlation with the source shape and with its time derivative
(capacitive coupling ~ dV/dt -> bipolar), fraction above 5 sigma noise.
Also: average normalised shapes (png)."""
import sys
import os
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use('Agg')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib.pyplot as plt
import argparse
import p2_qa_config as qa
import p2_mapping as pmap
import p2_waveforms as pw

ap = argparse.ArgumentParser()
ap.add_argument('run_key')
ap.add_argument('--max-chunk', type=int, default=5, help='decoded chunks 0..N (default 5)')
args = ap.parse_args()
cfg = qa.get_config(args.run_key)
T = cfg.out_dir('26_crosstalk') + '/'
# NO dropped channels here: the study must see the shorted/dead ones too
ct = pmap.build_channel_table(cfg.run_config_path, cfg.MAP_CSV_PATH, det_type=cfg.DET_TYPE,
                              det_name=cfg.DET_NAME, drop_connectors=cfg.DEAD_CONNECTORS,
                              strategy_overrides=cfg.STRATEGY_OVERRIDES)
padc = ct[ct['mapped']].drop_duplicates('channel_id').reset_index(drop=True)
import p2_align as pa
sfx = cfg.product_suffix(True)
rays = pd.read_csv(cfg.out_dir('06_efficiency') + f'/ray_hit_miss_list{sfx}.csv')
# pad -> M3 transform: refit from stage 06's live-pad footprint (row-aligned
# with the live channel table), then apply it to ALL pads incl. dropped ones
ct_live = pmap.build_channel_table(cfg.run_config_path, cfg.MAP_CSV_PATH, det_type=cfg.DET_TYPE,
                                   det_name=cfg.DET_NAME, drop_connectors=cfg.DEAD_CONNECTORS,
                                   drop_channels=cfg.DEAD_CHANNELS,
                                   strategy_overrides=cfg.STRATEGY_OVERRIDES)
pl = ct_live[ct_live['mapped']].drop_duplicates('channel_id')
fp = pd.read_csv(cfg.out_dir('06_efficiency') + f'/pad_footprint{sfx}.csv')
coef, *_ = np.linalg.lstsq(np.c_[pl.pad_cx, pl.pad_cy, np.ones(len(pl))], fp[['x', 'y']].to_numpy(), rcond=None)
q = np.c_[padc.pad_cx, padc.pad_cy, np.ones(len(padc))] @ coef
padc['mx'], padc['my'] = q[:, 0], q[:, 1]
pos = {(int(f), int(c)): (x, y) for f, c, x, y in zip(padc.feu, padc.channel, padc.mx, padc.my)}
tree = cKDTree(padc[['mx', 'my']].to_numpy())
d, i = tree.query(rays[['x', 'y']].to_numpy())
rays['sf'], rays['sc'] = padc.feu.to_numpy()[i], padc.channel.to_numpy()[i]
rays = rays[(d <= 5.0)].set_index('eventId')
print(f'rays within 5 mm of a pad centre: {len(rays)}')

OFFS = [-3, -2, -1, 1, 2, 3]
CTRL = [-19, -18, -17, 17, 18, 19]
files = pw.list_decoded_files(cfg.decoded_root_dir)
files = files[files.chunk <= args.max_chunk]
rows = []
shapes = {k: [] for k in ['src'] + [f'n{o:+d}' for o in OFFS] + ['ctrl']}
for _, fr in files.iterrows():
    feu = int(fr.feu)
    ped = pw.pedestal_from_data(fr.path)
    rms = pw.noise_rms_from_data(fr.path, ped)
    for ev, ts, ft, wf in pw.iter_decoded_events(fr.path, step=2000):
        sel = np.isin(ev, rays.index.to_numpy())
        if not sel.any():
            continue
        w = pw.subtract_common_noise(wf[sel] - ped[None, :, None])
        for k, e in enumerate(ev[sel]):
            r = rays.loc[e]
            if isinstance(r, pd.DataFrame):
                r = r.iloc[0]
            if int(r.sf) != feu:
                continue
            c = int(r.sc)
            src = w[k, c]
            a_src = src.max()
            if a_src < 10 * max(rms[c], 1.0) or a_src > 3800:   # clear, unsaturated
                continue
            t_src = int(np.argmax(src))
            dsrc = np.gradient(src)
            shapes['src'].append(src / a_src)
            for o, kind in [(o, 'nb') for o in OFFS] + [(o, 'ctrl') for o in CTRL]:
                cn = c + o
                if cn < 0 or cn >= 512 or cn // 64 != c // 64:
                    continue
                p = pos.get((feu, cn))
                if p is None or np.hypot(p[0] - r.x, p[1] - r.y) < 30:
                    continue
                nb = w[k, cn]
                corr_s = np.corrcoef(nb, src)[0, 1] if nb.std() > 0 else np.nan
                corr_d = np.corrcoef(nb, dsrc)[0, 1] if nb.std() > 0 else np.nan
                rows.append(dict(feu=feu, src=c, off=o, kind=kind, a_src=a_src,
                                 pk=nb.max() / a_src, tr=nb.min() / a_src,
                                 at_src_peak=nb[t_src] / a_src,
                                 dt=int(np.argmax(nb)) - t_src,
                                 corr_shape=corr_s, corr_deriv=corr_d,
                                 above5=nb.max() > 5 * max(rms[cn], 1.0), sig=nb.max() / max(rms[cn], 1.0)))
                if kind == 'nb':
                    shapes[f'n{o:+d}'].append(nb / a_src)
                else:
                    shapes['ctrl'].append(nb / a_src)
    print(f'  done {fr.path.split("/")[-1]}: rows {len(rows)}')

t = pd.DataFrame(rows)
t.to_csv(T + 'xtalk_rows.csv', index=False)
print(f'\nsource pulses: {len(shapes["src"])}, neighbour samples {len(t)}')
g = t.groupby(['kind', 'off']).agg(n=('pk', 'size'),
                                  frac_above5sig=('above5', 'mean'),
                                  amp_at_src_peak_med=('at_src_peak', 'median'),
                                  peak_ratio_med=('pk', 'median'),
                                  trough_ratio_med=('tr', 'median'),
                                  dt_med=('dt', 'median'),
                                  corr_shape_med=('corr_shape', 'median'),
                                  corr_deriv_med=('corr_deriv', 'median'))
print(g.round(3).to_string())
# coupling coefficient: slope of neighbour value at the source peak vs source amplitude
print('\nlinear coupling k (nb value at source-peak sample = k * A_src), robust:')
for (kind, o), s in t.groupby(['kind', 'off']):
    y = s.at_src_peak * s.a_src; x = s.a_src
    k = np.median(y / x)
    print(f'  {kind:4s} {o:+d}: k = {100*k:+.2f}%   (IQR {100*np.percentile(y/x,25):+.2f} .. {100*np.percentile(y/x,75):+.2f}%)')
# per FEU / odd-even of the source channel for +-1
t['src_parity'] = np.where(t.src % 2 == 0, 'even', 'odd')
print('\nch+-1 by FEU and source parity (frac >5 sigma, median k at source peak):')
print(t[(t.kind == 'nb') & (t.off.abs() == 1)].groupby(['feu', 'src_parity', 'off'])
      .agg(n=('pk', 'size'), frac5=('above5', 'mean'), k=('at_src_peak', 'median')).round(3).to_string())

fig, ax = plt.subplots(1, 2, figsize=(13, 4.5))
xs = np.arange(32) * 60
for k_, lst in shapes.items():
    if not lst:
        continue
    m = np.median(np.array(lst), axis=0)
    if k_ == 'src':
        ax[0].plot(xs, m, 'k-', lw=2, label=f'source (n={len(lst)})')
    else:
        ax[1].plot(xs, m, label=f'{k_} (n={len(lst)})', lw=2 if k_ in ('n-1', 'n+1') else 1)
ax[0].set_title('source pulse (median, normalised to its peak)')
ax[1].set_title('same-event neighbour channels / source peak (median)')
for a in ax:
    a.set_xlabel('time [ns]'); a.grid(alpha=.3)
ax[1].axhline(0, color='grey', lw=.5); ax[1].legend(fontsize=8)
fig.tight_layout(); fig.savefig(T + 'xtalk_shapes.png', dpi=120)
print('saved', T + 'xtalk_shapes.png')
