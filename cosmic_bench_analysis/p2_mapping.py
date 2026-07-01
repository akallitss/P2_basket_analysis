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
                        strategy='reverse'):
    """Return a DataFrame with one row per instrumented (feu, channel).

    Columns: feu, channel, dream_conn, within, connector_N, half, sector,
             strip, channel_id, pad_cx, pad_cy, radius, phi, mec8_connector,
             mapped (bool: channel_id found in the pad map).
    """
    wiring, feus, name = parse_dream_wiring(run_config_path, det_type, det_name)
    pad = load_pad_map(map_csv)

    rows = []
    for (feu, dream_conn), (n, half) in wiring.items():
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

    cols = ['pad_cx', 'pad_cy', 'radius', 'phi', 'mec8_connector']
    have = [c for c in cols if c in pad.columns]
    tab = tab.merge(pad[have], left_on='channel_id', right_index=True,
                    how='left')
    tab['mapped'] = tab['pad_cx'].notna() if 'pad_cx' in tab else False
    tab.attrs['det_name'] = name
    tab.attrs['feus'] = feus
    tab.attrs['strategy'] = strategy
    return tab


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
