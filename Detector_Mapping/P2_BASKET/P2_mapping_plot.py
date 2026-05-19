#!/usr/bin/env python3
"""
Visualize and validate the P2 Basket channel mapping CSV.

Three panels:
  1. Full routing overview  — connector pads → strip endpoints, colored by sector
  2. Strip electrode layout — strip (x,y) colored by strip number within sector
  3. Per-sector crossing check — lines within one sector should be non-crossing

Usage:
  python3 P2_mapping_plot.py P2_BASKET_mapping.csv
  python3 P2_mapping_plot.py P2_BASKET_mapping.csv --sector 0   # zoom one sector
  python3 P2_mapping_plot.py P2_BASKET_mapping.csv --out map.pdf
"""

import argparse
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.collections import LineCollection

SECTOR_COLORS = plt.cm.tab10.colors   # 10 distinct colors


def load_csv(path):
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f):
            rows.append({
                'pad_number':       int(r['pad_number']),
                'x':                float(r['x']),
                'y':                float(r['y']),
                'connector_number': int(r['connector_number']),
                'connector_pin':    int(r['connector_pin']),
                'pin_name':         r['pin_name'],
                'sector':           int(r['sector']),
                'strip':            int(r['strip']),
                'channel_id':       int(r['channel_id']),
                'pin_x':            float(r['pin_x']),
                'pin_y':            float(r['pin_y']),
                'pad_angle':        float(r['pad_angle']),
            })
    return rows


def make_segments(rows):
    """Build (N,2,2) array of line segments: connector pin → strip endpoint."""
    segs = np.array([[[r['pin_x'], r['pin_y']], [r['x'], r['y']]] for r in rows])
    return segs


def plot_overview(ax, rows, alpha_lines=0.3, lw=0.4):
    """Panel 1: full routing overview."""
    segs = make_segments(rows)
    colors = [SECTOR_COLORS[r['sector']] for r in rows]
    lc = LineCollection(segs, colors=colors, linewidths=lw, alpha=alpha_lines)
    ax.add_collection(lc)

    # Strip endpoints
    sx = [r['x']     for r in rows]
    sy = [r['y']     for r in rows]
    sc = [SECTOR_COLORS[r['sector']] for r in rows]
    ax.scatter(sx, sy, c=sc, s=4, zorder=3, linewidths=0)

    # Connector pads (one dot per unique position)
    px = [r['pin_x'] for r in rows]
    py = [r['pin_y'] for r in rows]
    ax.scatter(px, py, c='k', s=6, marker='s', zorder=4, linewidths=0)

    # Sector legend
    handles = [plt.Line2D([0], [0], color=SECTOR_COLORS[s], lw=2, label=f'Sector {s}')
               for s in range(10)]
    ax.legend(handles=handles, fontsize=6, loc='upper right', ncol=2)
    ax.set_aspect('equal')
    ax.set_title('Full routing: connector pins (■) → strip endpoints', fontsize=9)
    ax.set_xlabel('x (mm)'); ax.set_ylabel('y (mm)')


def plot_strip_layout(ax, rows):
    """Panel 2: strip electrode positions colored by strip number."""
    sx = np.array([r['x']     for r in rows])
    sy = np.array([r['y']     for r in rows])
    sv = np.array([r['strip'] for r in rows])   # 1–128

    sc = ax.scatter(sx, sy, c=sv, cmap='plasma', s=5, linewidths=0)
    plt.colorbar(sc, ax=ax, label='Strip number within sector')
    ax.set_aspect('equal')
    ax.set_title('Strip electrode layout (color = strip number)', fontsize=9)
    ax.set_xlabel('x (mm)'); ax.set_ylabel('y (mm)')


def plot_sector_detail(ax, rows, sector):
    """Panel 3: one-sector zoom — non-crossing lines = correct mapping."""
    sub = [r for r in rows if r['sector'] == sector]
    if not sub:
        ax.set_visible(False)
        return

    segs = make_segments(sub)
    strips = np.array([r['strip'] for r in sub])
    norm = plt.Normalize(strips.min(), strips.max())
    colors = cm.viridis(norm(strips))

    lc = LineCollection(segs, colors=colors, linewidths=0.8, alpha=0.8)
    ax.add_collection(lc)

    # Strip endpoints labeled with strip number
    for r in sub:
        ax.plot(r['x'], r['y'], '.', color=cm.viridis(norm(r['strip'])), ms=4)

    ax.autoscale()
    ax.set_aspect('equal')
    ax.set_title(f'Sector {sector} detail (color = strip number; lines must not cross)',
                 fontsize=9)
    ax.set_xlabel('x (mm)'); ax.set_ylabel('y (mm)')

    sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
    plt.colorbar(sm, ax=ax, label='Strip number')


def plot_pin_to_strip_distance(ax, rows):
    """Bonus panel: histogram of connector-pin → strip distances."""
    dists = [np.hypot(r['x'] - r['pin_x'], r['y'] - r['pin_y']) for r in rows]
    ax.hist(dists, bins=40, color='steelblue', edgecolor='white', linewidth=0.3)
    ax.set_xlabel('Pin → strip distance (mm)')
    ax.set_ylabel('Channels')
    ax.set_title(f'Routing length distribution  '
                 f'(min {min(dists):.0f}  median {np.median(dists):.0f}  '
                 f'max {max(dists):.0f} mm)', fontsize=9)


def main():
    parser = argparse.ArgumentParser(description='Visualize P2 Basket channel mapping')
    parser.add_argument('csv', help='P2_BASKET_mapping.csv')
    parser.add_argument('--sector', type=int, default=0,
                        help='Sector to zoom in panel 3 (default 0)')
    parser.add_argument('--out', default=None,
                        help='Save figure to file (PDF/PNG); default: show interactively')
    args = parser.parse_args()

    rows = load_csv(args.csv)
    print(f"Loaded {len(rows)} channels from {args.csv}")

    fig, axes = plt.subplots(2, 2, figsize=(16, 14))
    fig.suptitle('P2 Basket – channel mapping validation', fontsize=12, fontweight='bold')

    plot_overview(axes[0, 0], rows)
    plot_strip_layout(axes[0, 1], rows)
    plot_sector_detail(axes[1, 0], rows, args.sector)
    plot_pin_to_strip_distance(axes[1, 1], rows)

    plt.tight_layout()
    if args.out:
        fig.savefig(args.out, dpi=150, bbox_inches='tight')
        print(f"Saved → {args.out}")
    else:
        plt.show()


if __name__ == '__main__':
    main()