#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p2_mapping.py

Spatial mapping for the P2 (BASKET) pad detector: turn a DAQ hit address
``(feu, channel)`` into a physical pad position on the readout PCB.

Two independent pieces are joined here:

  1. The Gerber-derived pad map  (Detector_Mapping/P2_BASKET/P2_BASKET_mapping.csv)
     One row per pad, keyed by ``channel_id`` (0..1279), with
       - pad centre        : pad_cx, pad_cy            [mm]
       - polar coords      : radius, phi, delta_phi
       - logical address   : connector_number(=sector) 0..9, strip 1..128,
                             mec8_connector 0/1 (strip 1-64 vs 65-128), mec8_pin
     with the definition  channel_id = sector*128 + (strip-1).

  2. The DAQ electronics wiring  (run_config.json -> detector 'dream_feus')
     The P2 detector is read out with DREAM FEUs. Each of the 10 physical
     connectors is split into a 'bot' and a 'top' 64-channel half, each wired to
     one DREAM connector on some FEU, e.g.  c_1_bot -> (feu 3, dream_conn 1).
     A DREAM FEU channel decodes as
         dream_conn = channel // 64 + 1   (1-based)
         within     = channel %  64       (0..63 inside the connector)

Chain built by this module
---------------------------
    (feu, channel)
      -> dream_conn, within                      (DREAM decode)
      -> connector name 'c_<N>_<bot|top>'        (run_config dream_feus, inverted)
      -> physical connector N (1..10), half      (parse name)
      -> sector = N-1 , mec8_connector = 0/1     (bot=0/strips 1-64, top=1/strips 65-128)
      -> strip                                   (within-half order, see STRATEGY)
      -> channel_id = sector*128 + (strip-1)
      -> pad_cx, pad_cy, radius, phi             (map lookup)

The one genuinely uncertain link -- the order of ``within`` (0..63) inside a
64-strip half -- is what the mapping-validation stage exists to test, so it is
made an explicit, swappable STRATEGY:

  'linear'  : strip = half_base + within         (within 0->strip 1, 63->strip 64)
  'reverse' : strip = half_base + (63 - within)   ** VALIDATED DEFAULT **
  'pairswap': linear with adjacent pairs swapped (0<->1, 2<->3, ...), the
              pattern seen in the MEC8/VMM mec8_pin column; provided so the
              alternative can be rendered and compared.

'reverse' was validated against the M3 telescope on run p2_det1_long_run_6-30-26
(03_m3_alignment.py): with 'reverse' the per-event P2 pad centroid correlates
with the projected M3 track at r_x=0.93 / r_y=0.88 (clean 89 deg rotation,
scale ~1.0), whereas 'linear'/'pairswap' give r_y~0.15. The within-half wiring
is a fixed property of the DREAM/K59V readout, so 'reverse' should hold for any
P2 detector read out the same way. Flat cosmic occupancy alone cannot
distinguish the orderings (they all illuminate the fan uniformly) -- only the
track correlation can, which is why this stays swappable.
"""

import os
import json
import re

import numpy as np
import pandas as pd

CH_PER_CONNECTOR = 64      # one DREAM connector = 64 channels
STRIPS_PER_CONN = 128      # one physical connector = 128 strips (2 halves)
N_CONNECTORS = 10          # 10 physical connectors -> 1280 pads
N_PADS = N_CONNECTORS * STRIPS_PER_CONN


# --------------------------------------------------------------------------- #
# Pad map
# --------------------------------------------------------------------------- #
def load_pad_map(map_csv):
    """Load the Gerber pad map, indexed by channel_id.

    Returns a DataFrame indexed by channel_id (0..1279) with at least
    pad_cx, pad_cy, radius, phi, connector_number(=sector), strip,
    mec8_connector, mec8_pin, pad_area.
    """
    m = pd.read_csv(map_csv)
    if 'channel_id' not in m.columns:
        raise ValueError(f'{map_csv} has no channel_id column')
    m = m.set_index('channel_id').sort_index()
    # sanity: channel_id should be the full 0..1279 range
    if len(m) != N_PADS:
        print(f'[p2_mapping] WARNING: map has {len(m)} rows, expected {N_PADS}')
    return m


# --------------------------------------------------------------------------- #
# DAQ wiring from run_config.json
# --------------------------------------------------------------------------- #
_CONN_RE = re.compile(r'c_?(\d+)_(bot|top|tob)', re.IGNORECASE)


def parse_dream_wiring(run_config_path, det_type='P2', det_name=None):
    """Invert run_config 'dream_feus' into a (feu, dream_conn) -> (N, half) map.

    Returns
    -------
    wiring : dict {(feu:int, dream_conn:int): (connector_N:int, half:str)}
             half in {'bot','top'}.  connector_N is 1-based (1..10).
    feus   : sorted list of FEUs used by the detector.
    name   : the detector name from the config.
    """
    with open(run_config_path) as fh:
        cfg = json.load(fh)
    detectors = cfg.get('detectors', [])
    chosen = None
    for d in detectors:
        if det_name is not None and d.get('name') == det_name:
            chosen = d
            break
        if det_name is None and d.get('det_type') == det_type:
            chosen = d
            break
    if chosen is None:
        raise ValueError(f'No detector det_type={det_type!r}/name={det_name!r} '
                         f'in {run_config_path}')

    wiring = {}
    for conn_name, (feu, dream_conn) in chosen['dream_feus'].items():
        m = _CONN_RE.search(conn_name)
        if not m:
            print(f'[p2_mapping] WARNING: cannot parse connector name {conn_name!r}, skipping')
            continue
        n = int(m.group(1))
        half = m.group(2).lower()
        half = 'top' if half in ('top', 'tob') else 'bot'   # tolerate 'tob' typo
        wiring[(int(feu), int(dream_conn))] = (n, half)

    feus = sorted({feu for feu, _dc in wiring})
    return wiring, feus, chosen.get('name', det_type)


# --------------------------------------------------------------------------- #
# within-half ordering strategies
# --------------------------------------------------------------------------- #
def _within_to_strip(within, half, strategy='linear'):
    """Map within-connector DREAM index (0..63) + half to a strip (1..128)."""
    within = np.asarray(within)
    base = np.where(np.asarray(half) == 'top', 65, 1) if not isinstance(half, str) \
        else (65 if half == 'top' else 1)
    if strategy == 'linear':
        off = within
    elif strategy == 'reverse':
        off = 63 - within
    elif strategy == 'pairswap':
        off = within ^ 1                      # 0<->1, 2<->3, ...
    else:
        raise ValueError(f'unknown strategy {strategy!r}')
    return base + off


# --------------------------------------------------------------------------- #
# Full resolver: build a per-(feu,channel) pad lookup table
# --------------------------------------------------------------------------- #
def build_channel_table(run_config_path, map_csv, det_type='P2', det_name=None,
                        strategy='reverse', drop_connectors=()):
    """Return a DataFrame with one row per instrumented (feu, channel).

    Columns: feu, channel, dream_conn, within, connector_N, half, sector,
             strip, channel_id, pad_cx, pad_cy, radius, phi, mec8_connector,
             pad_w, pad_h, pad_angle, delta_phi (pad tile geometry, when the
             map CSV provides them), mapped (bool: channel_id found in the map).

    drop_connectors : iterable of physical connector numbers (1..10) that are
        disconnected/dead. Their channels are omitted entirely, so their hits do
        not map to pads and their pads never appear in the analysis.
    """
    wiring, feus, name = parse_dream_wiring(run_config_path, det_type, det_name)
    pad = load_pad_map(map_csv)
    drop = set(int(c) for c in drop_connectors)

    rows = []
    for (feu, dream_conn), (n, half) in wiring.items():
        if n in drop:
            continue
        sector = n - 1
        for within in range(CH_PER_CONNECTOR):
            channel = (dream_conn - 1) * CH_PER_CONNECTOR + within
            strip = int(_within_to_strip(within, half, strategy))
            channel_id = sector * STRIPS_PER_CONN + (strip - 1)
            rows.append((feu, channel, dream_conn, within, n, half,
                         sector, strip, channel_id))

    tab = pd.DataFrame(rows, columns=[
        'feu', 'channel', 'dream_conn', 'within', 'connector_N', 'half',
        'sector', 'strip', 'channel_id'])

    cols = ['pad_cx', 'pad_cy', 'radius', 'phi', 'mec8_connector',
            'pad_w', 'pad_h', 'pad_angle', 'delta_phi']
    have = [c for c in cols if c in pad.columns]
    tab = tab.merge(pad[have], left_on='channel_id', right_index=True,
                    how='left')
    tab['mapped'] = tab['pad_cx'].notna() if 'pad_cx' in tab else False
    tab.attrs['det_name'] = name
    tab.attrs['feus'] = feus
    tab.attrs['strategy'] = strategy
    tab.attrs['drop_connectors'] = sorted(drop)
    return tab


def pad_tiles(channel_table):
    """Real pad outlines for map rendering: each mapped pad as a rotated
    rectangle (pad_w x pad_h at pad_angle about its centre), drawn to scale.

    Returns (pads, verts): `pads` is the table deduplicated to one row per
    channel_id (mapped only), `verts` an (n_pads, 4, 2) corner array aligned
    row-for-row with `pads`. Requires the pad_w/pad_h/pad_angle columns
    (present when the map CSV carries the Gerber tile geometry).
    """
    t = channel_table[channel_table['mapped']].drop_duplicates('channel_id')
    a = np.radians(t['pad_angle'].to_numpy())
    w, h = t['pad_w'].to_numpy(), t['pad_h'].to_numpy()
    cx, cy = t['pad_cx'].to_numpy(), t['pad_cy'].to_numpy()
    ca, sa = np.cos(a), np.sin(a)
    corners = [np.stack([cx + dx * w * ca - dy * h * sa,
                         cy + dx * w * sa + dy * h * ca], axis=1)
               for dx, dy in ((-.5, -.5), (.5, -.5), (.5, .5), (-.5, .5))]
    return t, np.stack(corners, axis=1)


# --------------------------------------------------------------------------- #
# Insulation-mask pillars (mesh-support pillars = local dead spots)
# --------------------------------------------------------------------------- #
def load_pillars(gbr_path):
    """Pillar positions from the insulation-mask Gerber (KiCad, mm, same board
    frame as the pad map — verified by overlaying the mask fan on the pads).

    Pillars are drawn as full circles built from G02/G03 arc pairs: the arc
    centre is start + (I,J) and its radius |IJ|; the pen is wide enough that
    the stroke fills the disk, so the physical pillar radius is
    |IJ| + aperture/2. Two populations on the M2_V2 mask:
      small : arc r 0.20 mm, 0.4 mm pen -> 0.8 mm pillars (~11.7k, ~4 mm grid)
      big   : arc r 1.54 mm, 3.075 mm pen -> 6.15 mm pillars (exactly 5)
    Everything else in the file (fan outlines, text, fiducials) is not an arc
    at these radii and is ignored. Returns a DataFrame(x, y, r, big) or an
    empty one if the file is missing.
    """
    cols = ['x', 'y', 'r', 'big']
    if not gbr_path or not os.path.isfile(gbr_path):
        return pd.DataFrame(columns=cols)

    ap_def = re.compile(r'%ADD(\d+)C,([\d.]+)')          # circular apertures
    ap_sel = re.compile(r'^D(\d+)\*')
    move = re.compile(r'X(-?\d+)Y(-?\d+)D02\*')
    arc = re.compile(r'X(-?\d+)Y(-?\d+)I(-?\d+)J(-?\d+)D01\*')

    apertures, ap_r = {}, 0.0
    cur = None
    circles = {}                       # (cx,cy rounded) -> (cx, cy, r_outer)
    with open(gbr_path) as fh:
        for line in fh:
            m = ap_def.search(line)
            if m:
                apertures[int(m.group(1))] = float(m.group(2))
                continue
            m = ap_sel.match(line)
            if m:
                ap_r = apertures.get(int(m.group(1)), 0.0) / 2
                continue
            m = move.search(line)
            if m:
                cur = (int(m.group(1)) / 1e6, int(m.group(2)) / 1e6)
                continue
            m = arc.search(line)
            if m and cur is not None:
                i, j = int(m.group(3)) / 1e6, int(m.group(4)) / 1e6
                r_arc = float(np.hypot(i, j))
                # pillar circles only: small (~0.2) and big (~1.54) arc radii
                if 0.15 < r_arc < 0.30 or 1.0 < r_arc < 3.0:
                    cx, cy = cur[0] + i, cur[1] + j
                    key = (round(cx, 2), round(cy, 2))
                    circles.setdefault(key, (cx, cy, r_arc + ap_r))
                cur = (int(m.group(1)) / 1e6, int(m.group(2)) / 1e6)

    if not circles:
        return pd.DataFrame(columns=cols)
    df = pd.DataFrame(list(circles.values()), columns=['x', 'y', 'r'])
    df['big'] = df['r'] > 1.0
    return df


def draw_pillars(ax, pillars, transform=None, small=True, label=True):
    """Overlay the mask pillars on an axis. `transform` (optional) is a
    p2_align-style transform with .apply(x, y) and scale .s — used to carry the
    pad-frame pillar positions into the M3/reference frame of a plot. Big
    pillars are drawn at true size (red) plus a dashed locator ring; small
    pillars as faint dots."""
    from matplotlib.patches import Circle

    if pillars is None or not len(pillars):
        return
    x, y = pillars['x'].to_numpy(float), pillars['y'].to_numpy(float)
    r = pillars['r'].to_numpy(float)
    if transform is not None:
        x, y = transform.apply(x, y)
        r = r * transform.s
    big = pillars['big'].to_numpy(bool)
    if small and (~big).any():
        ax.scatter(x[~big], y[~big], s=0.4, c='deepskyblue', alpha=0.25,
                   linewidths=0, zorder=4,
                   label='small pillars' if label else None)
    first = True
    for xi, yi, ri in zip(x[big], y[big], r[big]):
        ax.add_patch(Circle((xi, yi), ri, facecolor='red', edgecolor='red',
                            alpha=0.9, zorder=5))
        ax.add_patch(Circle((xi, yi), 8.0 * (transform.s if transform else 1.0),
                            facecolor='none', edgecolor='red', ls='--', lw=1.0,
                            alpha=0.9, zorder=5,
                            label=('big pillars (5)' if label and first else None)))
        first = False


def has_tile_geometry(channel_table):
    """True when the channel table carries the pad tile geometry columns."""
    return {'pad_w', 'pad_h', 'pad_angle'} <= set(channel_table.columns)


def daq_lookup(channel_table):
    """Return a dict {(feu, channel): row-as-Series-dict} for fast hit joining,
    plus a convenience merge helper."""
    ct = channel_table.set_index(['feu', 'channel'])
    return ct


def attach_pads_to_hits(hits, channel_table):
    """Left-join a hits DataFrame (needs 'feu','channel') to pad positions.

    Returns hits with added columns from the channel table (channel_id, sector,
    strip, pad_cx, pad_cy, radius, phi, connector_N, half, mapped).
    """
    keep = ['feu', 'channel', 'channel_id', 'sector', 'strip', 'connector_N',
            'half', 'pad_cx', 'pad_cy', 'radius', 'phi', 'mapped']
    keep = [c for c in keep if c in channel_table.columns]
    return hits.merge(channel_table[keep], on=['feu', 'channel'], how='left')


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description='Inspect the P2 DAQ->pad mapping.')
    ap.add_argument('--run-config', required=True)
    ap.add_argument('--map-csv', required=True)
    ap.add_argument('--strategy', default='reverse',
                    choices=['linear', 'reverse', 'pairswap'])
    args = ap.parse_args()

    tab = build_channel_table(args.run_config, args.map_csv, strategy=args.strategy)
    print(f'detector {tab.attrs["det_name"]}  FEUs {tab.attrs["feus"]}  '
          f'strategy={tab.attrs["strategy"]}')
    print(f'{len(tab)} instrumented channels, '
          f'{int(tab["mapped"].sum())} resolved to a pad, '
          f'{int((~tab["mapped"]).sum())} unmapped')
    print(f'unique channel_id used: {tab["channel_id"].nunique()} '
          f'(collisions: {len(tab) - tab["channel_id"].nunique()})')
    print(tab.head(8).to_string(index=False))
