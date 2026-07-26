#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
record_mapping_alignment.py -- freeze the uRWELL mapping and alignment to disk.

Everything the P2 measurement rests on that is NOT in the raw data: the
channel -> strip wiring of the two uRWELLs, the front-to-back alignment, and
the uRWELL -> P2 pad-frame transform of each station.  The authoritative source
stays the code (`urw_lib.VIEW_MODE_DEFAULT` / `AXIS_FLIP_DEFAULT`, and the fits
inside `urw_p2_efficiency.py`), but the code is only useful to someone running
it.  This writes the same information out in three forms that outlive it:

  mapping_urwell.csv          one row per FEU channel: view, local position,
                              pitch, zone.  512 rows per detector.  This is the
                              mapping with no code in the way - anyone can join
                              on (detector, channel) and get a position.
  mapping_alignment.json      the wiring table, the front->back alignment, the
                              per-station frame transforms, and the beam optics,
                              with the measurement each number came from.
  MAPPING_AND_ALIGNMENT.md    the same thing for a human, with the evidence.

Run after urw_p2_efficiency.py, pointing at its output:

  unset PYTHONPATH
  $PY record_mapping_alignment.py \
      --eff-json out_eff/urw_p2_efficiency_highstat_eff_1.json \
      --out /local/home/banco/P2_data/TB_July2026_H4/analysis/urw_referenced_efficiency
"""
import os
import sys
import json
import argparse
from datetime import datetime

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import urw_lib as U  # noqa: E402

RUN_JSON_DEFAULT = ('/local/home/banco/P2_data/TB_July2026_H4/runs/'
                    'highstat_eff_1/run_config.json')
SUB_DEFAULT = 'beam_commissioning_00'


def channel_table(run_json, sub_run):
    """One row per (detector, FEU channel) with the position it maps to."""
    rows = []
    for name in U.URW_DETS:
        geo = U.UrwGeometry(name, run_json, sub_run_name=sub_run)
        for ch in range(U.N_CHAN):
            if not geo.mapped[ch]:
                continue
            rows.append(dict(
                detector=name, det_type=geo.det_type, feu=geo.feu_num,
                channel=ch, feu_connector=ch // 64 + 1, connector_channel=ch % 64,
                view=geo.view[ch], position_mm=round(float(geo.pos[ch]), 4),
                pitch_mm=float(geo.pitch[ch]),
                interpitch_mm=float(geo.interpitch[ch]),
                zone=int(geo.zone[ch]),
                view_mode=geo.view_mode[geo.view[ch]],
                axis_flipped=bool(geo.axis_flip.get(geo.view[ch], False)),
                z_mm=float(geo.center[2])))
    return pd.DataFrame(rows)


def summarise_frames(eff):
    """Per-station frame transform, averaged over the sub_runs, with spread."""
    out = {}
    for st in sorted({s['station'] for s in eff}):
        rows = [s for s in eff if s['station'] == st]
        f = [r['frame'] for r in rows]
        A = np.array([x['affine_A'] for x in f])
        t = np.array([x['affine_t'] for x in f])
        out[st] = dict(
            z_mm=rows[0]['z_mm'], n_sub_runs=len(rows),
            affine_A_mean=np.round(A.mean(0), 6).tolist(),
            affine_A_std=np.round(A.std(0), 6).tolist(),
            affine_t_mean_mm=np.round(t.mean(0), 4).tolist(),
            affine_t_std_mm=np.round(t.std(0), 4).tolist(),
            rotation_deg_mean=float(np.mean([x['affine_rotation_deg'] for x in f])),
            rotation_deg_std=float(np.std([x['affine_rotation_deg'] for x in f])),
            det_mean=float(np.mean([x['det'] for x in f])),
            stretch_xx_mean=float(np.mean([x['stretch_xx'] for x in f])),
            stretch_yy_mean=float(np.mean([x['stretch_yy'] for x in f])),
            stretch_xy_mean=float(np.mean([x['stretch_xy'] for x in f])),
            rigid_mean={k: float(np.mean([x['rigid'][k] for x in f]))
                        for k in ('dx', 'dy', 'theta_deg')},
            rigid_rmse_mm_mean=float(np.mean([x['rigid_rmse_mm'] for x in f])),
            residual_rms_x_mm=float(np.mean([r['residual']['x']['rms_mm'] for r in rows])),
            residual_rms_y_mm=float(np.mean([r['residual']['y']['rms_mm'] for r in rows])))
    return out


def summarise_front_back(eff):
    fb = [s['urwell']['front_back'] for s in eff]
    out = {}
    for ax in ('x', 'y'):
        out[ax] = {k: float(np.mean([f[ax][k] for f in fb]))
                   for k in ('slope', 'offset_mm', 'sigma_iqr_mm', 'rms_mm')}
        out[ax]['slope_std'] = float(np.std([f[ax]['slope'] for f in fb]))
    return out


def divergence(frames):
    """Fit stretch_ii(z) = 1 + z / L_i across the stations."""
    z = np.array([v['z_mm'] for v in frames.values()])
    out = {}
    for key, ax in (('stretch_xx_mean', 'x'), ('stretch_yy_mean', 'y')):
        v = np.array([frames[k][key] for k in frames])
        a, b = np.polyfit(z, v, 1)
        out[ax] = dict(d_scale_dz_per_mm=float(a), intercept=float(b),
                       virtual_source_m=float(1.0 / a / 1000.0))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run-json', default=RUN_JSON_DEFAULT)
    ap.add_argument('--sub-run', default=SUB_DEFAULT)
    ap.add_argument('--eff-json', default='out_eff/urw_p2_efficiency_highstat_eff_1.json')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    tab = channel_table(args.run_json, args.sub_run)
    csv_path = os.path.join(args.out, 'mapping_urwell.csv')
    tab.to_csv(csv_path, index=False)
    print(f'wrote {csv_path}  ({len(tab)} mapped channels)')

    with open(args.eff_json) as fh:
        eff = json.load(fh)
    frames = summarise_frames(eff)
    fb = summarise_front_back(eff)
    div = divergence(frames)

    rec = {
        'written': datetime.now().isoformat(timespec='seconds'),
        'source_run': eff[0]['run'],
        'n_sub_runs': len({s['sub_run'] for s in eff}),
        'urwell_wiring': {
            'view_mode': U.VIEW_MODE_DEFAULT,
            'axis_flip': U.AXIS_FLIP_DEFAULT,
            'meaning': {
                'AB': 'map connector 0 -> lower FEU connector, channel order as-is',
                'BA': 'the two connectors of the view interchanged',
                'AB_rev': 'as AB, channel order reversed inside each connector',
                'BA_rev': 'as BA, channel order reversed inside each connector'},
            'evidence': {
                'method': 'width of (back - front) per candidate wiring, with the '
                          'front as reference; explore4_back_map.py on '
                          'highstat_eff_1/beam_commissioning_00',
                'back_x_mm': {'AB': 45.6, 'BA': 5.42, 'AB_rev': 0.80, 'BA_rev': 47.0},
                'back_y_mm': {'AB': 35.2, 'BA': 6.40, 'AB_rev': 0.88, 'BA_rev': 38.8},
                'front_note': "the front's own data cannot separate BA/AB from the "
                              'mirror partner AB_rev/BA_rev; the choice recorded here '
                              'is the one that maps the front into the P2 pad frame by '
                              'a PROPER rotation (det > 0), which a reflection could '
                              'not do',
                'supersedes': 'CONNECTOR_SWAP_DEFAULT (removed 2026-07-26), which had '
                              "back x and back y as 'BA' and left the back pointing "
                              'at 4.4 mm'}},
        'urwell_front_to_back': dict(
            model='back = slope * front + offset, per axis, in local mm',
            dz_mm=eff[0]['urwell']['dz_mm'], **fb),
        'urwell_to_p2_frame': dict(
            model='(X, Y)_P2 = A @ (x, y)_uRWELL_track_at_z + t, both in mm',
            note='A is a free 2x2 fitted per station per sub_run; the values here '
                 'are the mean and spread over the sub_runs of the source run. '
                 'The rigid fit (dx, dy, theta) is reported alongside as the check '
                 'that A stayed orthogonal.',
            stations=frames),
        'beam_optics': dict(
            model='stretch_ii(z) = 1 + z / L_i, fitted across the three P2 z',
            note='the departure of A from orthogonality; a divergent beam, not a '
                 'detector distortion - see explore6_divergence.py',
            **div),
    }
    json_path = os.path.join(args.out, 'mapping_alignment.json')
    with open(json_path, 'w') as fh:
        json.dump(rec, fh, indent=2)
    print(f'wrote {json_path}')

    md_path = os.path.join(args.out, 'MAPPING_AND_ALIGNMENT.md')
    with open(md_path, 'w') as fh:
        fh.write(render_md(rec, tab, csv_path, json_path))
    print(f'wrote {md_path}')


def render_md(rec, tab, csv_path, json_path):
    L = []
    w = L.append
    w('# uRWELL mapping and alignment — the record\n')
    w(f'Written {rec["written"]} from `{rec["source_run"]}` '
      f'({rec["n_sub_runs"]} sub-runs).\n')
    w('Everything the uRWELL-referenced P2 measurement rests on that is not in the')
    w('raw data. The authoritative source is the code — `urw_lib.VIEW_MODE_DEFAULT`,')
    w('`urw_lib.AXIS_FLIP_DEFAULT`, and the fits in `urw_p2_efficiency.py` — and this')
    w('file is generated from it by `record_mapping_alignment.py`. Regenerate rather')
    w('than edit.\n')
    w('| file | what it is |')
    w('|---|---|')
    w(f'| `{os.path.basename(csv_path)}` | one row per FEU channel: view, local '
      'position, pitch, zone. Join on (detector, channel); no code needed. |')
    w(f'| `{os.path.basename(json_path)}` | the same numbers, machine readable. |')
    w('| this file | the same numbers, with the evidence. |\n')

    w('## 1. Channel → strip wiring\n')
    w('A uRWELL view is read out on **two** 64-channel Dream connectors, so the')
    w('wiring has two binary degrees of freedom and four possible answers:\n')
    w('| mode | meaning |')
    w('|---|---|')
    for k, v in rec['urwell_wiring']['meaning'].items():
        w(f'| `{k}` | {v} |')
    w('\nMeasured:\n')
    w('| detector | x view | y view | axis flip |')
    w('|---|---|---|---|')
    for det, m in rec['urwell_wiring']['view_mode'].items():
        fl = rec['urwell_wiring']['axis_flip'].get(det, {})
        flt = ', '.join(k for k, v in fl.items() if v) or 'none'
        w(f'| `{det}` | `{m["x"]}` | `{m["y"]}` | {flt} |')
    ev = rec['urwell_wiring']['evidence']
    w(f'\n**How.** {ev["method"]}. Width of `back − front` per candidate:\n')
    w('| view | `AB` | `BA` | `AB_rev` | `BA_rev` |')
    w('|---|---|---|---|---|')
    for view, key in (('back x', 'back_x_mm'), ('back y', 'back_y_mm')):
        d = ev[key]
        best = min(d, key=d.get)
        w(f'| {view} | ' + ' | '.join(
            (f'**{d[k]:.2f} mm**' if k == best else f'{d[k]:.2f} mm')
            for k in ('AB', 'BA', 'AB_rev', 'BA_rev')) + ' |')
    w(f'\n**Mirror ambiguity.** {ev["front_note"]}.\n')
    w(f'**Supersedes.** {ev["supersedes"]}.\n')
    w('The axis flip is a pure labelling choice: with the correct wiring the back')
    w('reads anti-parallel to the front on both views, so its two axes are mirrored')
    w('in software to keep the front→back slope at +1. No strip → position')
    w('assignment is touched.\n')

    fb = rec['urwell_front_to_back']
    w('## 2. Front → back alignment\n')
    w(f'`{fb["model"]}`, lever arm dz = {fb["dz_mm"]:.0f} mm.\n')
    w('| axis | slope | offset [mm] | core sigma [mm] |')
    w('|---|---|---|---|')
    for ax in ('x', 'y'):
        d = fb[ax]
        w(f'| {ax} | {d["slope"]:+.5f} ± {d["slope_std"]:.5f} | '
          f'{d["offset_mm"]:+.3f} | {d["sigma_iqr_mm"]:.2f} |')
    w('\nThe y slope is **not** a scale error — see §4. The core sigma is the two')
    w('planes\' resolutions in quadrature plus the beam\'s angular spread, so each')
    w('plane is better than ~0.6 mm.\n')

    w('## 3. uRWELL → P2 pad frame\n')
    w(f'`{rec["urwell_to_p2_frame"]["model"]}`\n')
    w(rec['urwell_to_p2_frame']['note'] + '\n')
    w('| station | z [mm] | rotation [deg] | det(A) | stretch xx / yy | shear | '
      'rigid rmse [mm] | residual rms x / y [mm] |')
    w('|---|---|---|---|---|---|---|---|')
    for st, v in rec['urwell_to_p2_frame']['stations'].items():
        w(f'| {st} | {v["z_mm"]:.0f} | {v["rotation_deg_mean"]:+.3f} ± '
          f'{v["rotation_deg_std"]:.3f} | {v["det_mean"]:+.4f} | '
          f'{v["stretch_xx_mean"]:.4f} / {v["stretch_yy_mean"]:.4f} | '
          f'{v["stretch_xy_mean"]:+.4f} | {v["rigid_rmse_mm_mean"]:.2f} | '
          f'{v["residual_rms_x_mm"]:.2f} / {v["residual_rms_y_mm"]:.2f} |')
    w('\nA proper rotation of about −60° with no reflection and no shear at every')
    w('station. Because a reflection is geometrically impossible here — the two')
    w('detectors are seen from the same side along the same beam — getting an')
    w('orthogonal `A` with det = +1 is a real test of both strip maps and the pad')
    w('map, not a fit artefact.\n')
    w('The full per-sub-run matrices are in the `frame` block of each entry of the')
    w('run\'s `urw_p2_efficiency_<run>.json`.\n')

    w('## 4. Beam optics\n')
    w(f'`{rec["beam_optics"]["model"]}` — {rec["beam_optics"]["note"]}.\n')
    w('| axis | d(scale)/dz [1/mm] | virtual source [m upstream] |')
    w('|---|---|---|')
    for ax in ('x', 'y'):
        d = rec['beam_optics'][ax]
        w(f'| {ax} | {d["d_scale_dz_per_mm"]:+.3e} | {d["virtual_source_m"]:+.1f} |')
    w('\nThe beam diverges in y and is essentially parallel in x. Extrapolating the')
    w('y term to dz = 1370 mm reproduces the front→back slope of §2, which is why')
    w('that slope is optics rather than a bad pitch. This is also why the applied')
    w('uRWELL → P2 transform is the affine and not the rigid one: otherwise a known')
    w('±1 mm optical term leaks into the residuals.\n')

    w('## 5. The mapping table\n')
    w(f'`{os.path.basename(csv_path)}` has {len(tab)} rows. Columns: `detector`,')
    w('`det_type`, `feu`, `channel` (0–511 global), `feu_connector`,')
    w('`connector_channel`, `view` (the coordinate the strip MEASURES), '
      '`position_mm`')
    w('(local, 0 at the low edge), `pitch_mm`, `interpitch_mm`, `zone`, `view_mode`,')
    w('`axis_flipped`, `z_mm`.\n')
    w('Zone summary:\n')
    z = (tab.groupby(['detector', 'view', 'zone', 'pitch_mm'], as_index=False)
         .agg(n_strips=('channel', 'size'), pos_min=('position_mm', 'min'),
              pos_max=('position_mm', 'max')))
    w('| detector | view | zone | pitch [mm] | strips | position range [mm] |')
    w('|---|---|---|---|---|---|')
    for _, r in z.iterrows():
        w(f'| {r["detector"]} | {r["view"]} | {r["zone"]} | {r["pitch_mm"]:.1f} | '
          f'{r["n_strips"]} | {r["pos_min"]:.2f} – {r["pos_max"]:.2f} |')
    w('\n> Reminder: in the raw map files `axis` is the direction the strip **runs**,')
    w('> so `axis=y` measures x. The `view` column here is already the measured')
    w('> coordinate, so no further inversion is needed.\n')
    return '\n'.join(L) + '\n'


if __name__ == '__main__':
    main()
