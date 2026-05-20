#!/usr/bin/env python3
"""
Compare P2 BASKET channel mappings:
  - Ours  : P2_BASKET_mapping.csv  (Gerber net-name approach)
  - Theirs: connector_<n>_v2.txt   (external algorithm, per-connector txt files)

Panels:
  1. Both point clouds overlaid — colored by connector
  2. Per-connector centroid offset vectors (quiver)
  3. Nearest-neighbour residuals after global translation — map view
  4. Residual histogram + CDF

Usage:
  python3 P2_mapping_compare.py P2_BASKET_mapping.csv <mappingP2_v2_dir>
  python3 P2_mapping_compare.py P2_BASKET_mapping.csv <mappingP2_v2_dir> --out compare.png
"""

import argparse
import csv
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from matplotlib.collections import LineCollection
from scipy.spatial import KDTree

SECTOR_COLORS = plt.cm.tab10.colors


# ── Loaders ────────────────────────────────────────────────────────────────────

def load_ours(path):
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f):
            if not r['x']:
                continue
            rows.append({
                'connector': int(r['connector_number']),
                'strip':     int(r['strip']),
                'x':         float(r['x']),
                'y':         float(r['y']),
                'pin_x':     float(r['pin_x']),
                'pin_y':     float(r['pin_y']),
            })
    return rows


def load_theirs(base_dir):
    rows = []
    for c in range(10):
        p = f'{base_dir}/connector_{c}_v2.txt'
        with open(p, encoding='utf-8-sig') as f:
            lines = f.readlines()
        hdr = [h.strip() for h in lines[0].split('\t')]
        for line in lines[1:]:
            parts = [v.strip() for v in line.split('\t')]
            row = dict(zip(hdr, parts))
            if not row.get('PinName', ''):
                continue
            m = re.search(r'\d+', row['PinName'])
            if not m:
                continue
            rows.append({
                'connector': int(row['Connector']),
                'strip':     int(m.group()),
                'x':         float(row['X']),
                'y':         float(row['Y']),
            })
    return rows


# ── Geometry helpers ───────────────────────────────────────────────────────────

def per_connector_centroids(rows):
    """Return arrays (connectors, cx_ours, cy_ours) for connectors 0-9."""
    pts = {}
    for r in rows:
        pts.setdefault(r['connector'], []).append([r['x'], r['y']])
    cs = sorted(pts)
    cx = np.array([np.mean([p[0] for p in pts[c]]) for c in cs])
    cy = np.array([np.mean([p[1] for p in pts[c]]) for c in cs])
    return np.array(cs), cx, cy


def find_global_translation(ours_pts, theirs_pts):
    """Simple centroid-based translation (rotation is negligible, ~0.01 deg)."""
    return ours_pts.mean(axis=0) - theirs_pts.mean(axis=0)


# ── Panels ────────────────────────────────────────────────────────────────────

def panel_overlay(ax, ours, theirs, alpha=0.6):
    """Both point clouds overlaid, colored by connector."""
    for c in range(10):
        o = [(r['x'], r['y']) for r in ours  if r['connector'] == c]
        t = [(r['x'], r['y']) for r in theirs if r['connector'] == c]
        col = SECTOR_COLORS[c]
        if o:
            ox, oy = zip(*o)
            ax.scatter(ox, oy, c=[col], s=5, linewidths=0, alpha=alpha,
                       label=f'C{c} (ours)' if c == 0 else f'C{c}')
        if t:
            tx, ty = zip(*t)
            ax.scatter(tx, ty, c=[col], s=5, marker='x', linewidths=0.6,
                       alpha=alpha)

    from matplotlib.lines import Line2D
    legend_els = [
        Line2D([0], [0], marker='o', color='gray', ms=5, ls='', label='Ours (●)'),
        Line2D([0], [0], marker='x', color='gray', ms=5, ls='', lw=0.6,
               label='Theirs (×)'),
    ]
    ax.legend(handles=legend_els, fontsize=7, loc='upper right')
    ax.set_aspect('equal')
    ax.set_title('Both mappings overlaid  (● ours, × theirs, color = connector)',
                 fontsize=9)
    ax.set_xlabel('x (mm)'); ax.set_ylabel('y (mm)')


def panel_quiver(ax, ours, theirs):
    """Per-connector centroid offset vectors (ours − theirs)."""
    _, ocx, ocy = per_connector_centroids(ours)
    _, tcx, tcy = per_connector_centroids(theirs)
    dx = ocx - tcx
    dy = ocy - tcy
    mag = np.hypot(dx, dy)

    norm = Normalize(mag.min(), mag.max())
    colors = cm.plasma(norm(mag))

    ax.quiver(tcx, tcy, dx, dy, color=colors, angles='xy',
              scale_units='xy', scale=1, width=0.004)
    ax.scatter(tcx, tcy, c='k', s=20, zorder=5)

    for c, (tx, ty, ddx, ddy, m) in enumerate(zip(tcx, tcy, dx, dy, mag)):
        ax.annotate(f'C{c}\n{m:.1f}mm', (tx + ddx/2, ty + ddy/2),
                    fontsize=6, ha='center', va='center')

    sm = plt.cm.ScalarMappable(cmap='plasma', norm=norm)
    plt.colorbar(sm, ax=ax, label='Offset magnitude (mm)', shrink=0.7)
    ax.set_aspect('equal')
    ax.autoscale()
    ax.set_title('Per-connector centroid offset  (ours − theirs)', fontsize=9)
    ax.set_xlabel('x (mm)'); ax.set_ylabel('y (mm)')


def panel_residual_map(ax, ours, theirs, shift):
    """After global translation: nearest-neighbour residual map."""
    t_pts = np.array([[r['x'], r['y']] for r in theirs]) + shift
    o_pts = np.array([[r['x'], r['y']] for r in ours])

    tree = KDTree(t_pts)
    dists, idxs = tree.query(o_pts)

    sc = ax.scatter([r['x'] for r in ours], [r['y'] for r in ours],
                    c=dists, cmap='hot_r', s=8, linewidths=0,
                    vmin=0, vmax=np.percentile(dists, 95))
    plt.colorbar(sc, ax=ax, label='NN residual after shift (mm)', shrink=0.7)
    ax.set_aspect('equal')
    ax.set_title(f'NN residual map after global translation\n'
                 f'(shift Δx={shift[0]:.1f}, Δy={shift[1]:.1f} mm)', fontsize=9)
    ax.set_xlabel('x (mm)'); ax.set_ylabel('y (mm)')


def panel_histogram(ax, ours, theirs, shift):
    """Residual histogram + cumulative fraction."""
    t_pts = np.array([[r['x'], r['y']] for r in theirs]) + shift
    o_pts = np.array([[r['x'], r['y']] for r in ours])

    tree = KDTree(t_pts)
    dists, _ = tree.query(o_pts)

    ax2 = ax.twinx()
    n, bins, _ = ax.hist(dists, bins=40, color='steelblue',
                         edgecolor='white', linewidth=0.3, label='Count')
    cdf = np.cumsum(n) / len(dists)
    ax2.plot(bins[1:], cdf, color='tomato', lw=1.5, label='CDF')
    ax2.set_ylim(0, 1.05)
    ax2.set_ylabel('Cumulative fraction', color='tomato')

    for thr in (0.5, 1.0, 2.0):
        frac = (dists < thr).mean()
        ax.axvline(thr, ls='--', lw=0.8, color='gray')
        ax.text(thr + 0.05, ax.get_ylim()[1] * 0.85,
                f'{frac*100:.0f}%\n<{thr}mm', fontsize=6, color='gray')

    ax.set_xlabel('NN residual (mm)')
    ax.set_ylabel('Channels')
    ax.set_title(f'Residual distribution  '
                 f'(median={np.median(dists):.2f} mm, '
                 f'max={dists.max():.2f} mm)', fontsize=9)


def panel_per_connector_offset(ax, ours, theirs):
    """Offset magnitude and direction angle per connector."""
    _, ocx, ocy = per_connector_centroids(ours)
    _, tcx, tcy = per_connector_centroids(theirs)
    dx = ocx - tcx; dy = ocy - tcy
    mag   = np.hypot(dx, dy)
    angle = np.degrees(np.arctan2(dy, dx))

    connectors = np.arange(10)
    ax.bar(connectors - 0.2, mag,   width=0.4, label='Offset magnitude (mm)', color='steelblue')
    ax2 = ax.twinx()
    ax2.plot(connectors, angle, 'o-', color='tomato', ms=5, lw=1.2,
             label='Offset direction (°)')
    ax2.set_ylabel('Offset direction (°)', color='tomato')
    ax.set_xticks(connectors); ax.set_xticklabels([f'C{c}' for c in connectors])
    ax.set_xlabel('Connector'); ax.set_ylabel('Magnitude (mm)')
    ax.set_title('Per-connector offset: constant magnitude, rotating direction\n'
                 '→ systematic radial offset between the two algorithms', fontsize=9)
    ax.legend(loc='upper left', fontsize=7)
    ax2.legend(loc='upper right', fontsize=7)


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description='Compare two P2 BASKET channel mappings')
    parser.add_argument('csv',     help='P2_BASKET_mapping.csv (our mapping)')
    parser.add_argument('ref_dir', help='Directory with connector_<n>_v2.txt files')
    parser.add_argument('--out',   default=None,
                        help='Save figure (PNG/PDF); default: show interactively')
    args = parser.parse_args()

    ours   = load_ours(args.csv)
    theirs = load_theirs(args.ref_dir)
    print(f"Loaded: ours={len(ours)}, theirs={len(theirs)}")

    o_pts = np.array([[r['x'], r['y']] for r in ours])
    t_pts = np.array([[r['x'], r['y']] for r in theirs])
    shift = find_global_translation(o_pts, t_pts)
    print(f"Global centroid shift (ours − theirs): Δx={shift[0]:.3f} Δy={shift[1]:.3f} mm")

    tree_shifted = KDTree(t_pts + shift)
    dists, _ = tree_shifted.query(o_pts)
    print(f"After shift — median NN dist: {np.median(dists):.3f} mm  max: {dists.max():.3f} mm")

    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('P2 BASKET mapping comparison  (Gerber net-name vs. external algorithm)',
                 fontsize=12, fontweight='bold')

    panel_overlay(axes[0, 0], ours, theirs)
    panel_quiver(axes[0, 1], ours, theirs)
    panel_per_connector_offset(axes[0, 2], ours, theirs)
    panel_residual_map(axes[1, 0], ours, theirs, shift)
    panel_histogram(axes[1, 1], ours, theirs, shift)

    # Summary text panel
    ax = axes[1, 2]
    ax.axis('off')
    _, ocx, ocy = per_connector_centroids(ours)
    _, tcx, tcy = per_connector_centroids(theirs)
    dx = ocx - tcx; dy = ocy - tcy
    mags = np.hypot(dx, dy)

    tree_raw = KDTree(t_pts)
    d_raw, _ = tree_raw.query(o_pts)
    tree_sh  = KDTree(t_pts + shift)
    d_sh, _  = tree_sh.query(o_pts)

    lines = [
        'Summary',
        '─' * 36,
        f'Channels:  ours {len(ours)} / theirs {len(theirs)}',
        '',
        'Raw (no alignment):',
        f'  median NN = {np.median(d_raw):.2f} mm',
        f'  max    NN = {d_raw.max():.2f} mm',
        '',
        f'After global centroid shift',
        f'  Δx={shift[0]:.2f} mm  Δy={shift[1]:.2f} mm:',
        f'  median NN = {np.median(d_sh):.2f} mm',
        f'  max    NN = {d_sh.max():.2f} mm',
        '',
        'Per-connector offset magnitude:',
        f'  mean={mags.mean():.2f} mm  σ={mags.std():.2f} mm',
        f'  range [{mags.min():.2f}, {mags.max():.2f}] mm',
        '',
        'Interpretation:',
        '  Constant offset magnitude with rotating',
        '  direction → systematic radial offset.',
        '  Both algorithms locate the same strips',
        '  but at different positions along their',
        '  radial extent (~11 mm apart).',
    ]
    ax.text(0.05, 0.95, '\n'.join(lines), transform=ax.transAxes,
            fontsize=8.5, va='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='#f0f4ff', alpha=0.8))

    plt.tight_layout()
    if args.out:
        fig.savefig(args.out, dpi=150, bbox_inches='tight')
        print(f"Saved → {args.out}")
    else:
        plt.show()


if __name__ == '__main__':
    main()