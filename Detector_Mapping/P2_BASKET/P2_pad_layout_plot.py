#!/usr/bin/env python3
"""
P2 Basket fan prototype — readout pad layout with dimensions.

Draws the true pad outlines (rotated rectangles reconstructed from pad_cx/pad_cy,
pad_w/pad_h and pad_angle), filled with one colour per connector.  No routing
lines — this is the detector-geometry view, not the mapping-validation view.

Annotated dimensions:
  - overall active-area bounding box (width x height)
  - fan apex, opening angle, inner/outer radius from the apex
  - pad size, radial pitch, total active pad area
  - total pad count per prototype (10 connectors x 128 pads = 1280)

Usage:
  python3 P2_pad_layout_plot.py P2_BASKET_mapping.csv
  python3 P2_pad_layout_plot.py P2_BASKET_mapping.csv --out layout.pdf
  python3 P2_pad_layout_plot.py P2_BASKET_mapping.csv --no-dims
"""

import argparse
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.patches import Arc

from P2_mapping import apex_from_strips

CONNECTOR_COLORS = plt.cm.tab10.colors
DIM_COLOR = '#444444'


def load_csv(path):
    with open(path) as f:
        rows = [r for r in csv.DictReader(f) if r['pad_cx'] not in ('', None)]
    return rows


def pad_polygon(row):
    """Pad outline as a 4x2 array: pad_w (tangential) x pad_h (radial),
    rotated by pad_angle about the pad centroid."""
    ang = np.deg2rad(float(row['pad_angle']))
    w, h = float(row['pad_w']), float(row['pad_h'])
    c, s = np.cos(ang), np.sin(ang)
    local = np.array([[-w / 2, -h / 2], [w / 2, -h / 2],
                      [w / 2, h / 2], [-w / 2, h / 2]])
    rot = np.array([[c, -s], [s, c]])
    return local @ rot.T + np.array([float(row['pad_cx']), float(row['pad_cy'])])


def geometry_summary(rows, polys, apex):
    """Envelope and pad statistics, all derived from the real pad outlines."""
    pts = np.vstack(polys)
    r = np.hypot(pts[:, 0] - apex[0], pts[:, 1] - apex[1])
    phi = np.degrees(np.arctan2(pts[:, 1] - apex[1], pts[:, 0] - apex[0]))

    area = np.array([float(x['pad_area']) for x in rows])
    pw = np.array([float(x['pad_w']) for x in rows])
    ph = np.array([float(x['pad_h']) for x in rows])

    # radial pitch: spacing of pad rows inside one sector, in that sector's frame.
    # Cluster the radial coordinate first — pads in the same row differ by a few
    # tenths of a mm, so a plain unique() would report that jitter as the pitch.
    sec0 = [x for x in rows if int(x['sector']) == 0]
    a = np.deg2rad(float(sec0[0]['pad_angle']))
    c, s = np.cos(-a), np.sin(-a)
    cx = np.array([float(x['pad_cx']) for x in sec0])
    cy = np.array([float(x['pad_cy']) for x in sec0])
    v = np.sort(s * cx + c * cy)
    tol = 0.5 * ph.mean()
    rows_v, cur = [], [v[0]]
    for val in v[1:]:
        if val - cur[-1] > tol:
            rows_v.append(np.mean(cur))
            cur = []
        cur.append(val)
    rows_v.append(np.mean(cur))
    radial_pitch = float(np.median(np.diff(rows_v)))

    return {
        'n_pads': len(rows),
        'n_conn': len({int(x['connector_number']) for x in rows}),
        'n_per_conn': len(rows) // len({int(x['connector_number']) for x in rows}),
        'xmin': pts[:, 0].min(), 'xmax': pts[:, 0].max(),
        'ymin': pts[:, 1].min(), 'ymax': pts[:, 1].max(),
        'width': pts[:, 0].max() - pts[:, 0].min(),
        'height': pts[:, 1].max() - pts[:, 1].min(),
        'r_min': r.min(), 'r_max': r.max(),
        'phi_min': phi.min(), 'phi_max': phi.max(),
        'opening': phi.max() - phi.min(),
        'pad_w': pw.mean(), 'pad_h': ph.mean(),
        'pad_area': area.mean(), 'active_cm2': area.sum() / 100.0,
        'radial_pitch': radial_pitch,
    }


def draw_pads(ax, rows, polys):
    conn = np.array([int(r['connector_number']) for r in rows])
    ax.add_collection(PolyCollection(
        polys, facecolors=[CONNECTOR_COLORS[c % 10] for c in conn],
        edgecolors='white', linewidths=0.15, zorder=2))

    # label each sector just outside its outer pad row
    for c in sorted(set(conn)):
        sub = [p for p, cc in zip(polys, conn) if cc == c]
        cen = np.array([p.mean(axis=0) for p in sub])
        far = cen[np.argmax(np.hypot(cen[:, 0], cen[:, 1]))]
        ax.annotate(f'C{c}', far, xytext=(14, 6), textcoords='offset points',
                    fontsize=8, fontweight='bold', color=CONNECTOR_COLORS[c % 10],
                    zorder=6)


def draw_dimensions(ax, g, apex):
    x0, x1, y0, y1 = g['xmin'], g['xmax'], g['ymin'], g['ymax']
    pad = 0.055 * max(g['width'], g['height'])

    # active-area bounding box
    ax.add_patch(plt.Rectangle((x0, y0), g['width'], g['height'], fill=False,
                               ec=DIM_COLOR, ls='--', lw=0.7, alpha=0.6, zorder=1))

    arrow = dict(arrowstyle='<->', color=DIM_COLOR, lw=1.0,
                 shrinkA=0, shrinkB=0)
    tbox = dict(boxstyle='round,pad=0.25', fc='white', ec='none', alpha=0.9)

    # width dimension, below
    yd = y0 - pad
    ax.annotate('', (x0, yd), xytext=(x1, yd), arrowprops=arrow, zorder=5)
    for xv in (x0, x1):
        ax.plot([xv, xv], [y0, yd], color=DIM_COLOR, lw=0.5, ls=':', zorder=1)
    ax.text((x0 + x1) / 2, yd + 5, f"{g['width']:.1f} mm", ha='center', va='bottom',
            fontsize=9, color=DIM_COLOR, bbox=tbox, zorder=5)

    # height dimension, left
    xd = x0 - pad
    ax.annotate('', (xd, y0), xytext=(xd, y1), arrowprops=arrow, zorder=5)
    for yv in (y0, y1):
        ax.plot([x0, xd], [yv, yv], color=DIM_COLOR, lw=0.5, ls=':', zorder=1)
    ax.text(xd - 6, (y0 + y1) / 2, f"{g['height']:.1f} mm", ha='right', va='center',
            rotation=90, fontsize=9, color=DIM_COLOR, bbox=tbox, zorder=5)

    # fan apex and the two edge rays
    ax.plot(*apex, marker='x', ms=9, mew=1.8, color='k', zorder=6)
    ax.annotate('fan apex\n(%.1f, %.1f)' % (apex[0], apex[1]), apex,
                xytext=(6, 12), textcoords='offset points', ha='left',
                fontsize=8, color='k', zorder=6)
    for a in (g['phi_min'], g['phi_max']):
        t = np.deg2rad(a)
        ax.plot([apex[0], apex[0] + g['r_max'] * np.cos(t)],
                [apex[1], apex[1] + g['r_max'] * np.sin(t)],
                color=DIM_COLOR, lw=0.7, ls='-.', alpha=0.55, zorder=1)

    # opening angle arc
    rad = 0.30 * g['r_max']
    ax.add_patch(Arc(apex, 2 * rad, 2 * rad, theta1=g['phi_min'],
                     theta2=g['phi_max'], color=DIM_COLOR, lw=1.0, zorder=5))
    tm = np.deg2rad((g['phi_min'] + g['phi_max']) / 2)
    ax.text(apex[0] + rad * np.cos(tm) * 1.06, apex[1] + rad * np.sin(tm) * 1.06,
            f"{g['opening']:.1f}°", fontsize=9, color=DIM_COLOR,
            ha='center', va='center', bbox=tbox, zorder=6)

    # inner / outer radius arcs + radial dimension along the bisector
    for rr in (g['r_min'], g['r_max']):
        ax.add_patch(Arc(apex, 2 * rr, 2 * rr, theta1=g['phi_min'] - 2,
                         theta2=g['phi_max'] + 2, color=DIM_COLOR,
                         lw=0.7, ls=':', alpha=0.8, zorder=1))
    p_in = apex + g['r_min'] * np.array([np.cos(tm), np.sin(tm)])
    p_out = apex + g['r_max'] * np.array([np.cos(tm), np.sin(tm)])
    ax.annotate('', p_out, xytext=p_in, arrowprops=arrow, zorder=5)
    mid = (p_in + p_out) / 2
    ax.text(mid[0], mid[1], f"R {g['r_min']:.0f} → {g['r_max']:.0f} mm",
            fontsize=9, color=DIM_COLOR, ha='center', va='center',
            rotation=np.degrees(tm), rotation_mode='anchor', bbox=tbox, zorder=6)


def info_text(g):
    """Two-column summary block; kept short so it fits the bottom margin band."""
    left = [
        f"Pad size     {g['pad_w']:.2f} × {g['pad_h']:.2f} mm",
        f"Pad area     {g['pad_area']:.1f} mm² (mean)",
        f"Radial pitch {g['radial_pitch']:.2f} mm",
        '',
    ]
    right = [
        f"Envelope     {g['width']:.1f} × {g['height']:.1f} mm",
        f"Radial span  {g['r_min']:.1f} – {g['r_max']:.1f} mm",
        f"Opening      {g['opening']:.1f}°",
        f"Active area  {g['active_cm2']:.0f} cm²",
    ]
    head = (f"TOTAL PADS:  {g['n_pads']}"
            f"   ({g['n_conn']} connectors × {g['n_per_conn']} pads)")
    body = [f"{l:<32}{r}" for l, r in zip(left, right)]
    return '\n'.join([head, '─' * 62] + body)


def main():
    ap = argparse.ArgumentParser(description='P2 Basket fan pad-layout figure')
    ap.add_argument('csv', help='P2_BASKET_mapping.csv')
    ap.add_argument('--out', default=None, help='Save to PNG/PDF (default: show)')
    ap.add_argument('--no-dims', action='store_true',
                    help='Draw the pads only, without dimension annotations')
    args = ap.parse_args()

    rows = load_csv(args.csv)
    polys = [pad_polygon(r) for r in rows]
    apex = apex_from_strips(rows)
    g = geometry_summary(rows, polys, apex)

    print(f"Loaded {g['n_pads']} pads  |  {g['n_conn']} connectors "
          f"× {g['n_per_conn']} pads")
    print(f"Apex ({apex[0]:.2f}, {apex[1]:.2f}) mm  |  envelope "
          f"{g['width']:.1f} × {g['height']:.1f} mm  |  opening {g['opening']:.2f} deg")
    print(f"Radial span {g['r_min']:.1f}–{g['r_max']:.1f} mm  |  "
          f"active pad area {g['active_cm2']:.0f} cm²")

    fig, ax = plt.subplots(figsize=(12, 11))
    draw_pads(ax, rows, polys)
    if not args.no_dims:
        draw_dimensions(ax, g, apex)

    handles = [plt.Line2D([0], [0], marker='s', ls='', ms=7,
                          color=CONNECTOR_COLORS[c % 10], label=f'Connector {c}')
               for c in range(g['n_conn'])]
    ax.legend(handles=handles, fontsize=8, loc='upper right',
              ncol=2, framealpha=0.9, title='Readout connector',
              title_fontsize=8.5)

    # summary sits in the bottom margin band, clear of the dimension lines
    ax.text(0.015, 0.018, info_text(g), transform=ax.transAxes, fontsize=9,
            va='bottom', ha='left', family='monospace',
            bbox=dict(boxstyle='round', fc='#f0f4ff', ec='#b8c4e0', alpha=0.95))

    ax.set_aspect('equal')
    ax.autoscale_view()
    ax.set_xlim(min(g['xmin'], apex[0]) - 0.07 * g['width'],
                g['xmax'] + 0.10 * g['width'])
    ax.set_ylim(min(g['ymin'], apex[1]) - 0.28 * g['height'],
                g['ymax'] + 0.05 * g['height'])
    ax.set_xlabel('x (mm)')
    ax.set_ylabel('y (mm)')
    ax.set_title(f"P2 BASKET fan prototype — readout pad layout "
                 f"({g['n_pads']} pads)", fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.15, lw=0.5)

    fig.tight_layout()
    if args.out:
        fig.savefig(args.out, dpi=200, bbox_inches='tight')
        print(f"Saved → {args.out}")
    else:
        plt.show()


if __name__ == '__main__':
    main()
