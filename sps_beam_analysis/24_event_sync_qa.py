#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
24_event_sync_qa.py

Event-alignment QA across the FEUs of every telescope station, per sub_run.
Every station's FEU receives the SAME external (TCM / scintillator) trigger
stream, so the combined-hits `eventId` and `trigger_timestamp_ns` must line up
FEU-to-FEU. This stage is the guard that every correlation-based stage (21
alignment, 22 tag-and-probe) rides on: if two FEUs slip by N triggers, cluster
correlations silently pair the wrong events.

Method (streamed, one pass per sub_run)
---------------------------------------
For each station FEU, reduce the combined hits to a per-event record
(eventId -> n_hits, first trigger_timestamp_ns), memory-bounded via p2_io.
Then, for the reference FEU vs every other FEU:
  * event counts and the shared-eventId overlap fraction;
  * the per-common-event timestamp delta (should be ~0 -- same trigger);
  * a hit-count-per-event cross-correlation as a function of an INTEGER
    event-index shift, over the common eventId range. A clean, aligned stream
    peaks sharply at shift 0; a peak at shift N flags an off-by-N slip.

The eventId is the DAQ's shared trigger counter, so alignment is checked on it
directly; the cross-correlation over shifts is the independent, ID-agnostic
cross-check (it would catch a per-FEU counter reset that eventId alone hides).

Products (<Analysis>/telescope/<run>/<sub_run>/24_event_sync_qa/):
  event_sync_<sub_run>.png    counts, timestamp-delta hist, x-corr vs shift
  event_sync_<sub_run>.json   machine-readable verdict + numbers

Usage:
  python3 24_event_sync_qa.py [run_key] [--sub-run NAME] [--max-shift 5]
"""

import os
import json
import argparse

import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import sps_config as sc
import p2_io as p2io


def per_feu_events(hits_dir, feu):
    """Streamed per-event reduction for ONE FEU: DataFrame indexed by eventId
    with n_hits and the first trigger_timestamp_ns seen for that event."""
    parts = []
    for df in p2io.iter_hits(hits_dir, ['eventId', 'feu', 'trigger_timestamp_ns'],
                             feus=[feu], progress=False):
        if not len(df):
            continue
        g = df.groupby('eventId')
        parts.append(pd.DataFrame({'n_hits': g.size(),
                                   'ts': g['trigger_timestamp_ns'].first()}))
    if not parts:
        return pd.DataFrame(columns=['n_hits', 'ts'])
    out = pd.concat(parts)
    out = out.groupby(level=0).agg(n_hits=('n_hits', 'sum'), ts=('ts', 'min'))
    return out.sort_index()


def xcorr_vs_shift(a, b, max_shift):
    """Pearson correlation of two per-event hit-count series (aligned on a
    common eventId index) as a function of integer eventId shift of b."""
    common = a.index.intersection(b.index)
    if len(common) < 10:
        return np.arange(-max_shift, max_shift + 1), \
               np.full(2 * max_shift + 1, np.nan)
    lo, hi = int(common.min()), int(common.max())
    full = np.arange(lo, hi + 1)
    va = a['n_hits'].reindex(full, fill_value=0).to_numpy(float)
    vb = b['n_hits'].reindex(full, fill_value=0).to_numpy(float)
    shifts = np.arange(-max_shift, max_shift + 1)
    corr = np.full(len(shifts), np.nan)
    for i, s in enumerate(shifts):
        if s >= 0:
            x, y = va[s:], vb[:len(vb) - s] if s else vb
        else:
            x, y = va[:len(va) + s], vb[-s:]
        if len(x) > 10 and x.std() > 0 and y.std() > 0:
            corr[i] = float(np.corrcoef(x, y)[0, 1])
    return shifts, corr


def analyse_subrun(cfg, sub_run, max_shift):
    hits_dir = cfg.combined_hits_dir(sub_run)
    # (det_tag, feu) -> per-event frame
    feu_of = {}
    order = []
    for det in cfg.detectors():
        for feu in det.feus:
            print(f'    reducing {det.name} FEU {feu} ...', flush=True)
            feu_of[(det.name, feu)] = per_feu_events(hits_dir, feu)
            order.append((det.name, feu))
    if not order:
        return None

    counts = {f'{n} (FEU{f})': int(len(feu_of[(n, f)])) for n, f in order}
    # reference = the FEU with the most events (usually the illuminated plane)
    ref_key = max(order, key=lambda k: len(feu_of[k]))
    ref = feu_of[ref_key]

    pairs = []
    for key in order:
        if key == ref_key:
            continue
        b = feu_of[key]
        common = ref.index.intersection(b.index)
        dt = (ref.loc[common, 'ts'].astype(np.int64) -
              b.loc[common, 'ts'].astype(np.int64)) if len(common) else \
            pd.Series(dtype=np.int64)
        shifts, corr = xcorr_vs_shift(ref, b, max_shift)
        best_shift = int(shifts[np.nanargmax(corr)]) if np.isfinite(corr).any() \
            else None
        pairs.append({
            'name': f'{key[0]} FEU{key[1]}',
            'n_events': int(len(b)),
            'n_common': int(len(common)),
            'overlap_frac': float(len(common) / max(1, len(b))),
            'dt_ns_median': float(dt.median()) if len(dt) else None,
            'dt_ns_p95_abs': float(dt.abs().quantile(0.95)) if len(dt) else None,
            'best_shift': best_shift,
            'shifts': shifts.tolist(),
            'corr': [None if not np.isfinite(c) else float(c) for c in corr],
            'dt_ns': dt.to_numpy() if len(dt) else np.array([]),
        })

    aligned = all((p['best_shift'] == 0) for p in pairs if p['best_shift']
                  is not None) and \
        all((p['dt_ns_median'] is None or abs(p['dt_ns_median']) < 1e6)
            for p in pairs)
    return {'sub_run': sub_run, 'ref': f'{ref_key[0]} FEU{ref_key[1]}',
            'counts': counts, 'pairs': pairs, 'aligned': bool(aligned)}


def plot(cfg, res, out_png):
    pairs = res['pairs']
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))

    # 1) event counts per FEU
    ax = axes[0]
    labels = list(res['counts'].keys())
    ax.bar(range(len(labels)), [res['counts'][k] for k in labels],
           color='steelblue')
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=25, ha='right', fontsize=8)
    ax.set_ylabel('events with >=1 hit')
    ax.set_title('event counts per FEU')
    ax.grid(True, axis='y', alpha=0.3)

    # 2) timestamp-delta histogram (common events, ref - other)
    ax = axes[1]
    any_dt = False
    for p in pairs:
        dt = np.asarray(p['dt_ns'], dtype=float)
        if len(dt) > 5:
            any_dt = True
            ax.hist(dt / 1e3, bins=60, histtype='step', lw=1.4,
                    label=f'{p["name"]} (n={len(dt)})')
    ax.set_xlabel(r'$t_{ref} - t_{other}$ per common event [$\mu$s]')
    ax.set_ylabel('common events')
    ax.set_title('trigger-timestamp delta (should be ~0)')
    ax.grid(True, alpha=0.3)
    if any_dt:
        ax.legend(fontsize=8)

    # 3) cross-correlation vs event-index shift
    ax = axes[2]
    for p in pairs:
        s = np.asarray(p['shifts'], float)
        c = np.asarray([np.nan if v is None else v for v in p['corr']], float)
        ax.plot(s, c, 'o-', lw=1.4, ms=4, label=p['name'])
    ax.axvline(0, color='grey', ls='--', lw=1)
    ax.set_xlabel('event-index shift of other FEU')
    ax.set_ylabel('hit-count correlation')
    ax.set_title('x-corr vs shift (peak at 0 = aligned)')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    verdict = 'ALIGNED' if res['aligned'] else 'CHECK — possible slip'
    fig.suptitle(f'Event-sync QA — {cfg.RUN} / {res["sub_run"]}   '
                 f'[ref {res["ref"]}]   ->  {verdict}', fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out_png, dpi=150, bbox_inches='tight')
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description='Per-sub_run FEU event-sync QA.')
    ap.add_argument('run_key', nargs='?', default=sc.DEFAULT_RUN)
    ap.add_argument('--sub-run', default=None,
                    help='sub_run name (default: all discovered sub_runs).')
    ap.add_argument('--max-shift', type=int, default=5,
                    help='max event-index shift scanned for the x-corr check.')
    args = ap.parse_args()

    cfg = sc.get_config(args.run_key)
    print(cfg)
    subruns = [args.sub_run] if args.sub_run else cfg.find_subruns()
    if not subruns:
        print('No sub_runs with combined hits found.')
        return

    for sub_run in subruns:
        print(f'  sub_run {sub_run}:')
        res = analyse_subrun(cfg, sub_run, args.max_shift)
        if res is None:
            print('    no FEUs / no data, skipping')
            continue
        out_dir = cfg.out_dir(sc.TELESCOPE_TAG, sub_run, '24_event_sync_qa')
        plot(cfg, res, os.path.join(out_dir, f'event_sync_{sub_run}.png'))
        # strip the bulky raw dt arrays before dumping json
        js = {k: v for k, v in res.items() if k != 'pairs'}
        js['pairs'] = [{k: v for k, v in p.items() if k != 'dt_ns'}
                       for p in res['pairs']]
        with open(os.path.join(out_dir, f'event_sync_{sub_run}.json'), 'w') as f:
            json.dump(js, f, indent=2)
        print(f'    ref {res["ref"]}, counts {res["counts"]}')
        for p in res['pairs']:
            print(f'    {p["name"]}: {p["n_common"]} common '
                  f'({p["overlap_frac"]:.1%} of its events), '
                  f'median dt = '
                  f'{p["dt_ns_median"]/1e3 if p["dt_ns_median"] is not None else float("nan"):.1f} us, '
                  f'best shift = {p["best_shift"]}')
        print(f'    verdict: {"ALIGNED" if res["aligned"] else "CHECK"}')
        print(f'    -> {out_dir}')


if __name__ == '__main__':
    main()
