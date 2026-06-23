#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FanPadDetectorTest.py

Visual geometry tests for FanPadDetector.

Loads the *full* detector (all 10 FPC connectors, 1280 pads) from the
Gerber-derived CSV and runs two display tests:

  1. Sector layout  — pads coloured by FPC connector number (0–9)
  2. Strip ordering — pads coloured by strip number within each connector (0–127)

Run from the P2_basket_analysis directory:
    python3 FanPadDetectorTest.py

Note: hit heatmaps, efficiency maps, and orientation comparison are handled
by vmm_detector_efficiency.py with real run data.

@author: ak271430
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as mcm

from FanPadDetector import FanPadDetector

# ── Configuration ──────────────────────────────────────────────

MAP_CSV = os.path.join(
    os.path.dirname(__file__),
    "Detector_Mapping", "P2_BASKET", "P2_BASKET_mapping.csv",
)
HALF_WIDTH = 5.0   # mm — angular half-width of each pad
N_CONNECTORS = 10  # full detector


# ── Helpers ────────────────────────────────────────────────────

def _set_ax_limits(ax, pads, margin=20):
    all_x = [p.outer_x for p in pads] + [p.via_x for p in pads]
    all_y = [p.outer_y for p in pads] + [p.via_y for p in pads]
    ax.set_xlim(min(all_x) - margin, max(all_x) + margin)
    ax.set_ylim(min(all_y) - margin, max(all_y) + margin)
    ax.set_aspect("equal")
    ax.set_xlabel("X (mm)")
    ax.set_ylabel("Y (mm)")


# ── 1. Sector layout: pads coloured by FPC connector number ────

def plot_sector_layout(det):
    """Draw all pads coloured by FPC connector number (0–9, tab10)."""
    fig, ax = plt.subplots(figsize=(14, 10))
    cmap = plt.colormaps["tab10"].resampled(N_CONNECTORS)

    for pad in det.pads:
        verts = pad.vertices()
        poly  = plt.Polygon(verts,
                            facecolor=cmap(pad.vmm_id),
                            edgecolor="black", linewidth=0.2, alpha=0.85)
        ax.add_patch(poly)

    _set_ax_limits(ax, det.pads)

    handles = [mpatches.Patch(color=cmap(i), label=f"FPC {i}")
               for i in range(N_CONNECTORS)]
    ax.legend(handles=handles, fontsize=8, loc="lower right",
              ncol=2, title="Connector")
    ax.set_title(f"Full detector — sector layout\n"
                 f"({len(det.pads)} pads, {N_CONNECTORS} FPC connectors)",
                 fontsize=12, fontweight="bold")
    return fig, ax


# ── 2. Strip ordering: pads coloured by strip number (0–127) ───

def plot_strip_order(det):
    """Draw all pads coloured by strip index within the connector (0–127)."""
    fig, ax = plt.subplots(figsize=(14, 10))
    cmap = plt.cm.plasma
    norm = plt.Normalize(vmin=0, vmax=127)

    for pad in det.pads:
        verts = pad.vertices()
        poly  = plt.Polygon(verts,
                            facecolor=cmap(norm(pad.ch)),
                            edgecolor="black", linewidth=0.2, alpha=0.9)
        ax.add_patch(poly)

    _set_ax_limits(ax, det.pads)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ax=ax, shrink=0.75, label="Strip index within connector (0 = strip 1)")
    ax.set_title(f"Full detector — strip ordering within each connector\n"
                 f"({len(det.pads)} pads)",
                 fontsize=12, fontweight="bold")
    return fig, ax


# ── Main ───────────────────────────────────────────────────────

def main():
    det = FanPadDetector.from_mapping_csv_full(MAP_CSV, half_width_mm=HALF_WIDTH)
    assert len(det.pads) > 0, "No pads loaded — check CSV path"

    connectors_found = sorted(set(p.vmm_id for p in det.pads))
    print(f"Connectors in CSV: {connectors_found}")
    print(f"Strips per connector: {len(det.pads) // len(connectors_found)}")

    fig1, _ = plot_sector_layout(det)
    fig2, _ = plot_strip_order(det)

    plt.tight_layout()
    print("All tests passed")
    plt.show()


if __name__ == "__main__":
    main()