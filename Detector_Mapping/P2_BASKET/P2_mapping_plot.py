#!/usr/bin/env python3
"""
Visualize and validate the P2 Basket channel mapping CSV.

Panels (always shown):
  1. Full routing overview  — connector pads → strip endpoints, colored by sector
  2. Strip electrode layout — strip (x,y) colored by strip number within sector
  3. Per-sector crossing check — lines within one sector should be non-crossing
  4. Routing length histogram

Extra panels (shown when mec8_connector / mec8_pin columns are present):
  5. Strip layout colored by MEC8 connector (JM_1 vs JM_2) — should split cleanly in space
  6. MEC8 pin vs strip number — should be two monotonic lines, one per connector

Usage:
  python3 P2_mapping_plot.py P2_BASKET_mapping.csv
  python3 P2_mapping_plot.py P2_BASKET_mapping.csv --sector 3
  python3 P2_mapping_plot.py P2_BASKET_mapping.csv --out map.pdf
"""

import argparse
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.collections import LineCollection

SECTOR_COLORS = plt.cm.tab10.colors
MEC8_COLORS   = ['#2196F3', '#FF5722']   # blue = JM_1, orange = JM_2


def load_csv(path):
    rows = []
    has_mec8 = False
    with open(path) as f:
        reader = csv.DictReader(f)
        has_mec8 = 'mec8_connector' in reader.fieldnames
        for r in reader:
            row = {
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
            }
            if has_mec8:
                row['mec8_connector'] = int(r['mec8_connector'])
                row['mec8_pin']       = int(r['mec8_pin'])
            rows.append(row)
    return rows, has_mec8


def make_segments(rows):
    return np.array([[[r['pin_x'], r['pin_y']], [r['x'], r['y']]] for r in rows])


# ── Panel 1: full routing overview ─────────────────────────────────────────────

def plot_overview(ax, rows, alpha_lines=0.3, lw=0.4):
    segs   = make_segments(rows)
    colors = [SECTOR_COLORS[r['sector']] for r in rows]
    ax.add_collection(LineCollection(segs, colors=colors, linewidths=lw, alpha=alpha_lines))
    ax.scatter([r['x']     for r in rows], [r['y']     for r in rows],
               c=[SECTOR_COLORS[r['sector']] for r in rows], s=4, zorder=3, linewidths=0)
    ax.scatter([r['pin_x'] for r in rows], [r['pin_y'] for r in rows],
               c='k', s=6, marker='s', zorder=4, linewidths=0)
    handles = [plt.Line2D([0], [0], color=SECTOR_COLORS[s], lw=2, label=f'Sector {s}')
               for s in range(10)]
    ax.legend(handles=handles, fontsize=6, loc='upper right', ncol=2)
    ax.set_aspect('equal')
    ax.set_title('Full routing: connector pins (■) → strip endpoints', fontsize=9)
    ax.set_xlabel('x (mm)'); ax.set_ylabel('y (mm)')


# ── Panel 2: strip layout colored by strip number ──────────────────────────────

def plot_strip_layout(ax, rows):
    sx = np.array([r['x']     for r in rows])
    sy = np.array([r['y']     for r in rows])
    sv = np.array([r['strip'] for r in rows])
    sc = ax.scatter(sx, sy, c=sv, cmap='plasma', s=5, linewidths=0)
    plt.colorbar(sc, ax=ax, label='Strip number within sector')
    ax.set_aspect('equal')
    ax.set_title('Strip electrode layout (color = strip number)', fontsize=9)
    ax.set_xlabel('x (mm)'); ax.set_ylabel('y (mm)')


# ── Panel 3: per-sector crossing check ────────────────────────────────────────

def plot_sector_detail(ax, rows, sector):
    sub = [r for r in rows if r['sector'] == sector]
    if not sub:
        ax.set_visible(False)
        return
    strips = np.array([r['strip'] for r in sub])
    norm   = plt.Normalize(strips.min(), strips.max())
    colors = cm.viridis(norm(strips))
    ax.add_collection(LineCollection(make_segments(sub), colors=colors,
                                     linewidths=0.8, alpha=0.8))
    ax.scatter([r['x'] for r in sub], [r['y'] for r in sub],
               c=cm.viridis(norm(strips)), s=8, zorder=3, linewidths=0)
    ax.autoscale(); ax.set_aspect('equal')
    ax.set_title(f'Sector {sector} — routing detail (non-crossing = correct)', fontsize=9)
    ax.set_xlabel('x (mm)'); ax.set_ylabel('y (mm)')
    plt.colorbar(plt.cm.ScalarMappable(cmap='viridis', norm=norm),
                 ax=ax, label='Strip number')


# ── Panel 4: routing length histogram ─────────────────────────────────────────

def plot_distance_hist(ax, rows):
    dists = [np.hypot(r['x'] - r['pin_x'], r['y'] - r['pin_y']) for r in rows]
    ax.hist(dists, bins=40, color='steelblue', edgecolor='white', linewidth=0.3)
    ax.set_xlabel('Pin → strip distance (mm)')
    ax.set_ylabel('Channels')
    ax.set_title(f'Routing length  (min {min(dists):.0f}  '
                 f'median {np.median(dists):.0f}  max {max(dists):.0f} mm)', fontsize=9)


# ── Panel 5: strip layout colored by MEC8 connector ───────────────────────────

def plot_mec8_spatial(ax, rows):
    """
    Strip positions colored by MEC8 connector (JM_1 vs JM_2).
    Validation: each sector should show a clean spatial split, not random mixing.
    """
    for mec_idx, label in [(0, 'JM_1 (connector 0)'), (1, 'JM_2 (connector 1)')]:
        sub = [r for r in rows if r['mec8_connector'] == mec_idx]
        ax.scatter([r['x'] for r in sub], [r['y'] for r in sub],
                   c=MEC8_COLORS[mec_idx], s=5, linewidths=0, label=label, alpha=0.8)
    ax.legend(fontsize=8, loc='upper right')
    ax.set_aspect('equal')
    ax.set_title('Strip layout by MEC8 connector  '
                 '(each sector should split cleanly in two)', fontsize=9)
    ax.set_xlabel('x (mm)'); ax.set_ylabel('y (mm)')


# ── Panel 6: MEC8 pin vs strip number ─────────────────────────────────────────

def plot_mec8_pin_order(ax, rows):
    """
    MEC8 pin number vs strip number, one point per unique strip (1–128).
    Each connector should show a monotonic line — disorder means wrong routing.
    """
    # Collapse to unique strip numbers (same across all sectors)
    strip_info = {}
    for r in rows:
        s = r['strip']
        if s not in strip_info:
            strip_info[s] = (r['mec8_connector'], r['mec8_pin'])

    for mec_idx, label in [(0, 'JM_1 (connector 0)'), (1, 'JM_2 (connector 1)')]:
        pts = sorted([(s, pin) for s, (c, pin) in strip_info.items() if c == mec_idx])
        strips_x, pins_y = zip(*pts)
        ax.plot(strips_x, pins_y, 'o', ms=4, color=MEC8_COLORS[mec_idx],
                label=label, alpha=0.8)

    ax.set_xlabel('Strip number (1–128)')
    ax.set_ylabel('MEC8 pin number')
    ax.set_title('MEC8 pin vs strip number  (should be monotonic per connector)',
                 fontsize=9)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description='Visualize P2 Basket channel mapping')
    parser.add_argument('csv', help='P2_BASKET_mapping.csv')
    parser.add_argument('--sector', type=int, default=0,
                        help='Sector to zoom in panel 3 (default 0)')
    parser.add_argument('--out', default=None,
                        help='Save to file (PDF/PNG); default: show interactively')
    args = parser.parse_args()

    rows, has_mec8 = load_csv(args.csv)
    print(f"Loaded {len(rows)} channels  |  MEC8 columns: {'yes' if has_mec8 else 'no'}")

    if has_mec8:
        fig, axes = plt.subplots(3, 2, figsize=(16, 20))
        fig.suptitle('P2 Basket – channel mapping validation (with K59V MEC8 mapping)',
                     fontsize=12, fontweight='bold')
        plot_overview(axes[0, 0], rows)
        plot_strip_layout(axes[0, 1], rows)
        plot_sector_detail(axes[1, 0], rows, args.sector)
        plot_distance_hist(axes[1, 1], rows)
        plot_mec8_spatial(axes[2, 0], rows)
        plot_mec8_pin_order(axes[2, 1], rows)
    else:
        fig, axes = plt.subplots(2, 2, figsize=(16, 14))
        fig.suptitle('P2 Basket – channel mapping validation',
                     fontsize=12, fontweight='bold')
        plot_overview(axes[0, 0], rows)
        plot_strip_layout(axes[0, 1], rows)
        plot_sector_detail(axes[1, 0], rows, args.sector)
        plot_distance_hist(axes[1, 1], rows)

    plt.tight_layout()
    if args.out:
        fig.savefig(args.out, dpi=150, bbox_inches='tight')
        print(f"Saved → {args.out}")
    else:
        plt.show()


if __name__ == '__main__':
    main()