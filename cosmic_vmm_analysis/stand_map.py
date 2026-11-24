#!/usr/bin/env python3
"""
Which VMM chips read which chamber of the lab cosmic stand, for the analyses that need chip -> chamber
(coincidence_check.CH / CHAMBER and through it efficiency_tp, channel_efficiency.DETV, probe_diagnostics.DETV).
Keys: 3 = P2_IN (z 320 mm), 2 = P2_MID (z 630), 1 = P2_OUT (z 940).

Pad positions are NOT affected by the layouts below: each chip keeps its connector number and the chambers have the
same design, so orientation_check.Mapper (pad x, y per chip) stays as it is. Only chamber membership (z, label) moves.

Layout from the environment variable VMM_STAND (set it in the pipeline area's env.sh, one layout per area):
  (unset) / 'default'          runs 1-406: chips 2-7 IN, 8-13 MID, 14-19 OUT
  'swap_mid_out_20261124'      runs >= 407 (setup entry 2026-11-24 12:14): adaptation cards + cables moved between
                               P2_MID (det5) and P2_OUT (det1): chips 14-19 read MID, 8-13 read OUT
@author: Alexandra Kallitsopoulou
"""
import os

LAYOUTS = {
    'default': {3: range(2, 8), 2: range(8, 14), 1: range(14, 20)},
    'swap_mid_out_20261124': {3: range(2, 8), 2: range(14, 20), 1: range(8, 14)},
}
LAYOUT = os.environ.get('VMM_STAND', '') or 'default'
if LAYOUT not in LAYOUTS:
    raise SystemExit(f'stand_map: unknown VMM_STAND={LAYOUT!r}; known: {sorted(LAYOUTS)}')
CHIPS = LAYOUTS[LAYOUT]
