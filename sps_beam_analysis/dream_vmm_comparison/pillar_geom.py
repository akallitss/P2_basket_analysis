#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""pillar_geom.py -- the P2 bulk-pillar geometry, straight from the mask that
was ordered.

`P2_Mask2.gbr` (CERN bulk soldermask, 2026-05-13) is the file the bulk was made
from, and it is the definition of the amplification-gap dead area.  Parsing it
is four lines of gerber and gives an exact answer, so nothing here is a design
recollection:

    small   D0.500 mm   41 366 of them, on an EXACT 2.000 mm square lattice
                        at even-integer mm in the board frame
    medium  D4.600 mm   40, around the rim and on the fan
    big     D6.150 mm   5

The board frame is the pad-map frame -- `P2_BASKET_mapping.csv` spans
x[70.2, 589.1], y[10.9, 504.7] against the small-pillar array's
x[70, 588], y[10, 504] -- so a track position in a station's pad frame can be
compared with a pillar position with no transform in between.  That is what
makes the lattice a *prediction* rather than a fit: the modulation is at
k = (pi, 0) and (0, pi) rad/mm with phase zero, no free parameters.

Coverage: pi*(0.25)^2 / 2.000^2 = 4.91 % of the amplification area.

    python3 pillar_geom.py --build     # re-parse the gerber into data/
"""
import os
import re

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
NPZ = os.path.join(DATA, "p2_pillars_gerber.npz")
GERBER = os.path.expanduser(
    "~/CLionProjects/p2_geant/design/gerbers/bulk_masks_CERN/P2_Mask2.gbr")

PITCH = 2.000            # mm, exact
DIA = 0.500              # mm, exact
COVERAGE = np.pi * (DIA / 2) ** 2 / PITCH ** 2      # 0.04909


# --------------------------------------------------------------------------- #
def parse_gerber(path=GERBER):
    """Flash positions by aperture.  FSLAX26Y26 + MOIN, so coordinates are
    inches at 1e-6, and the aperture definitions carry the true diameters."""
    adef = re.compile(r"%ADD(\d+)([CROP]),([\d.X]+)\*%")
    asel = re.compile(r"^D(\d+)\*")
    flash = re.compile(r"^(?:X(-?\d+))?(?:Y(-?\d+))?D03\*")
    scale = 25.4 / 1e6
    ap, cur, out, x, y = {}, None, {}, 0, 0
    for line in open(path):
        line = line.strip()
        m = adef.match(line)
        if m:
            ap[int(m.group(1))] = float(m.group(3).split("X")[0]) * 25.4
            continue
        m = asel.match(line)
        if m:
            cur = int(m.group(1))
            continue
        m = flash.match(line)
        if m:
            if m.group(1) is not None:
                x = int(m.group(1))
            if m.group(2) is not None:
                y = int(m.group(2))
            out.setdefault(cur, []).append((x * scale, y * scale))
    return {ap[k]: np.array(v) for k, v in out.items() if k in ap}


def build(path=GERBER, out=NPZ):
    by_dia = parse_gerber(path)
    small = np.array([])
    blocks = {}
    for d, pos in by_dia.items():
        tag = "small" if d < 1.0 else ("big" if d > 5.5 else "medium")
        blocks[f"{tag}_xy"] = pos
        blocks[f"{tag}_dia"] = np.array([d])
        if tag == "small":
            small = pos
    # the lattice claim, checked rather than asserted
    r = np.abs(small - PITCH * np.round(small / PITCH)).max()
    blocks["lattice_residual_mm"] = np.array([r])
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.savez_compressed(out, **blocks)
    print(f"{out}: " + ", ".join(
        f"{t} n={len(blocks[t + '_xy'])} D={blocks[t + '_dia'][0]:.3f} mm"
        for t in ("small", "medium", "big")))
    print(f"  lattice residual {r:.2e} mm; coverage {COVERAGE:.4%}")
    return blocks


def load(path=NPZ):
    if not os.path.isfile(path):
        raise SystemExit(f"{path} missing -- run `python3 pillar_geom.py --build`")
    z = np.load(path)
    return {k: z[k] for k in z.files}


# --------------------------------------------------------------------------- #
def big_pillar_mask(x, y, geom, margin=0.0):
    """True where a track is NOT shadowed by a medium or big pillar.

    Only these are masked: they are millimetres across, so a track landing on
    one is a genuine hole in the acceptance and shows up as a dead patch of a
    single pad.  The 0.5 mm ones are 4.9 % of the area everywhere and cannot be
    masked away -- they are the subject of the second half of this analysis,
    not a defect to be removed.
    """
    keep = np.ones(len(x), bool)
    for tag in ("medium", "big"):
        p, d = geom[f"{tag}_xy"], float(geom[f"{tag}_dia"][0])
        r = d / 2 + margin
        for px, py in p:
            if abs(px - np.median(x)) > 200 or abs(py - np.median(y)) > 200:
                continue
            keep &= (x - px) ** 2 + (y - py) ** 2 > r * r
    return keep


def nearest_node(x, y, pitch=PITCH):
    """Offset from the nearest lattice node, in [-pitch/2, pitch/2)."""
    return (x + pitch / 2) % pitch - pitch / 2, (y + pitch / 2) % pitch - pitch / 2


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", action="store_true")
    ap.add_argument("--gerber", default=GERBER)
    a = ap.parse_args()
    if a.build:
        build(a.gerber)
    else:
        g = load()
        for t in ("small", "medium", "big"):
            print(f"{t}: n={len(g[t + '_xy'])} D={g[t + '_dia'][0]:.3f} mm")
        print(f"lattice residual {g['lattice_residual_mm'][0]:.2e} mm")
