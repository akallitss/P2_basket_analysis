"""
triangle_geometry_check.py
Interactive geometry check for the equilateral triangle FEM simulation.

Matches the area of the original P2 basket trapezoid shape (~169 256 mm²).
Draws the triangle boundary, wall BC zone, pillar grid, and centroid.

Pure Python — no FreeCAD needed.

Run:
    py FreeCAD/triangle_geometry_check.py

Output:
    FreeCAD/results/triangle_geometry_check.html
"""

import os
import math

# ---------------------------------------------------------------------------
# Triangle geometry — same area as the original P2 basket shape
# ---------------------------------------------------------------------------
AREA_TARGET_MM2 = 169_256          # mm²  measured from original mesh

SIDE_MM   = math.sqrt(4 * AREA_TARGET_MM2 / math.sqrt(3))   # 625.2 mm
HEIGHT_MM = SIDE_MM * math.sqrt(3) / 2                       # 541.4 mm

# Vertices: base at y=0, apex pointing up, left corner at origin
V = [
    (0.0,            0.0       ),   # bottom-left
    (SIDE_MM,        0.0       ),   # bottom-right
    (SIDE_MM / 2.0,  HEIGHT_MM ),   # apex
]

CENTROID = (
    sum(v[0] for v in V) / 3,
    sum(v[1] for v in V) / 3,
)

# ---------------------------------------------------------------------------
# BC parameters  (same as original simulation)
# ---------------------------------------------------------------------------
WALL_MM         = 1.0    # mm  border width → wall BC
PILLAR_PITCH_MM = 3.0    # mm  pillar centre-to-centre
PILLAR_SNAP_MM  = 1.5    # mm  snap radius for pillar BC

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(SCRIPT_DIR, "results")
OUTPUT_HTML = os.path.join(RESULTS_DIR, "triangle_geometry_check.html")
PLOTLY_CDN  = "https://cdn.plot.ly/plotly-2.35.2.min.js"


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def _cross_z(ax, ay, bx, by, cx, cy):
    """Z-component of (b-a) × (c-a)."""
    return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)


def point_in_triangle(px, py, vertices):
    """True if (px, py) is strictly inside or on the triangle."""
    ax, ay = vertices[0]
    bx, by = vertices[1]
    cx, cy = vertices[2]
    d1 = _cross_z(ax, ay, bx, by, px, py)
    d2 = _cross_z(bx, by, cx, cy, px, py)
    d3 = _cross_z(cx, cy, ax, ay, px, py)
    has_neg = (d1 < 0) or (d2 < 0) or (d3 < 0)
    has_pos = (d1 > 0) or (d2 > 0) or (d3 > 0)
    return not (has_neg and has_pos)


def dist_point_to_line(px, py, x1, y1, x2, y2):
    """Perpendicular distance from point to the infinite line through (x1,y1)-(x2,y2)."""
    dx, dy = x2 - x1, y2 - y1
    length = math.hypot(dx, dy)
    if length < 1e-12:
        return math.hypot(px - x1, py - y1)
    return abs(dx * (y1 - py) - dy * (x1 - px)) / length


def min_dist_to_boundary(px, py):
    """Minimum perpendicular distance from point to any edge of the outer triangle."""
    edges = [(V[0], V[1]), (V[1], V[2]), (V[2], V[0])]
    return min(dist_point_to_line(px, py, a[0], a[1], b[0], b[1]) for a, b in edges)


def in_active_area(px, py):
    """True if point is inside the triangle AND more than WALL_MM from every edge."""
    return (point_in_triangle(px, py, V) and
            min_dist_to_boundary(px, py) > WALL_MM)


def inset_triangle(dist):
    """
    Return vertices of the triangle inset by `dist` mm perpendicular to each edge.
    For an equilateral triangle the inset is a scaled copy about the centroid.
      inradius r = SIDE / (2*sqrt(3))
      scale_factor = (r - dist) / r
    """
    r = SIDE_MM / (2.0 * math.sqrt(3))
    if dist >= r:
        return None          # degenerate — inset larger than triangle
    scale = (r - dist) / r
    cx, cy = CENTROID
    return [(cx + (vx - cx) * scale, cy + (vy - cy) * scale) for vx, vy in V]


# ---------------------------------------------------------------------------
# Pillar centre generation
# ---------------------------------------------------------------------------

def generate_pillar_centers():
    """
    Rectangular grid of pillar centres at PILLAR_PITCH_MM spacing,
    keeping only those inside the active area (> WALL_MM from every edge).
    """
    xmin = min(v[0] for v in V)
    xmax = max(v[0] for v in V)
    ymin = min(v[1] for v in V)
    ymax = max(v[1] for v in V)

    pillars = []
    pcx = xmin + WALL_MM + PILLAR_PITCH_MM / 2.0
    while pcx <= xmax - WALL_MM:
        pcy = ymin + WALL_MM + PILLAR_PITCH_MM / 2.0
        while pcy <= ymax - WALL_MM:
            if in_active_area(pcx, pcy):
                pillars.append((pcx, pcy))
            pcy += PILLAR_PITCH_MM
        pcx += PILLAR_PITCH_MM
    return pillars


# ---------------------------------------------------------------------------
# Wall zone sample points (for shading)
# ---------------------------------------------------------------------------

def wall_zone_points(n=4000):
    """
    Random-ish sample of points that are inside the triangle but in the
    wall zone (within WALL_MM of any edge), for shading purposes.
    Uses a dense grid rather than random to avoid import dependencies.
    """
    xmin = min(v[0] for v in V)
    xmax = max(v[0] for v in V)
    ymin = min(v[1] for v in V)
    ymax = max(v[1] for v in V)

    step = 4.0   # mm grid spacing for wall zone markers
    xs, ys = [], []
    x = xmin
    while x <= xmax:
        y = ymin
        while y <= ymax:
            if point_in_triangle(x, y, V):
                d = min_dist_to_boundary(x, y)
                if d <= WALL_MM:
                    xs.append(x)
                    ys.append(y)
            y += step
        x += step
    return xs, ys


# ---------------------------------------------------------------------------
# JSON helper
# ---------------------------------------------------------------------------

def _jv(v):
    if v is None:
        return "null"
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float)):
        if math.isnan(v) or math.isinf(v):
            return "null"
        return repr(v)
    if isinstance(v, str):
        return '"' + v.replace('\\', '\\\\').replace('"', '\\"').replace('\n', '\\n') + '"'
    if isinstance(v, (list, tuple)):
        return "[" + ",".join(_jv(i) for i in v) + "]"
    if isinstance(v, dict):
        return "{" + ",".join('"' + k + '":' + _jv(vv) for k, vv in v.items()) + "}"
    return '"' + str(v) + '"'


# ---------------------------------------------------------------------------
# Build HTML
# ---------------------------------------------------------------------------

def build_html(pillars, wall_xs, wall_ys, inner_v):

    # --- Outer triangle boundary ---
    outer_xs = [v[0] for v in V] + [V[0][0]]
    outer_ys = [v[1] for v in V] + [V[0][1]]

    trace_outer = _jv({
        "type": "scatter", "mode": "lines",
        "name": "Triangle boundary",
        "x": outer_xs, "y": outer_ys,
        "line": {"color": "#2c3e50", "width": 2},
        "showlegend": True,
    })

    # --- Inner (active area) boundary ---
    inner_xs = [v[0] for v in inner_v] + [inner_v[0][0]]
    inner_ys = [v[1] for v in inner_v] + [inner_v[0][1]]

    trace_inner = _jv({
        "type": "scatter", "mode": "lines",
        "name": "Active area boundary (1 mm inset)",
        "x": inner_xs, "y": inner_ys,
        "line": {"color": "#e74c3c", "width": 1.5, "dash": "dash"},
        "showlegend": True,
    })

    # --- Wall zone shading ---
    trace_wall = _jv({
        "type": "scatter", "mode": "markers",
        "name": "Wall BC zone (1 mm border)",
        "x": wall_xs, "y": wall_ys,
        "marker": {"color": "#e74c3c", "size": 2, "opacity": 0.35},
        "showlegend": True,
    })

    # --- Pillar centres ---
    px = [p[0] for p in pillars]
    py = [p[1] for p in pillars]
    hover = ["Pillar ({:.1f}, {:.1f})".format(p[0], p[1]) for p in pillars]

    trace_pillars = _jv({
        "type": "scatter", "mode": "markers",
        "name": "Pillar centres ({} total)".format(len(pillars)),
        "x": px, "y": py,
        "marker": {"color": "#2980b9", "size": 3, "opacity": 0.7},
        "hovertext": hover,
        "hovertemplate": "%{hovertext}<extra></extra>",
        "showlegend": True,
    })

    # --- Centroid marker ---
    trace_centroid = _jv({
        "type": "scatter", "mode": "markers+text",
        "name": "Centroid",
        "x": [CENTROID[0]], "y": [CENTROID[1]],
        "marker": {"color": "#e67e22", "size": 12, "symbol": "cross"},
        "text": ["Centroid<br>({:.0f}, {:.0f})".format(CENTROID[0], CENTROID[1])],
        "textposition": "top center",
        "showlegend": True,
    })

    all_traces = "[{},{},{},{},{}]".format(
        trace_outer, trace_inner, trace_wall, trace_pillars, trace_centroid)

    # --- Stats table rows ---
    inradius = SIDE_MM / (2 * math.sqrt(3))
    perimeter_active = 3 * (SIDE_MM - 2 * WALL_MM * math.sqrt(3))
    stats = [
        ("Outer side length", "{:.1f} mm".format(SIDE_MM)),
        ("Outer height", "{:.1f} mm".format(HEIGHT_MM)),
        ("Target area", "{:,.0f} mm² = {:.0f} cm²".format(AREA_TARGET_MM2, AREA_TARGET_MM2 / 100)),
        ("Inradius", "{:.1f} mm".format(inradius)),
        ("Wall border width", "{:.1f} mm".format(WALL_MM)),
        ("Pillar pitch", "{:.1f} mm".format(PILLAR_PITCH_MM)),
        ("Pillar snap radius", "{:.1f} mm".format(PILLAR_SNAP_MM)),
        ("Number of pillar centres", "{:,}".format(len(pillars))),
        ("Wall zone area (approx)", "{:,.0f} mm²".format(
            AREA_TARGET_MM2 - math.sqrt(3) / 4 * (SIDE_MM - 2 * WALL_MM * math.sqrt(3)) ** 2)),
    ]
    table_rows = "\n".join(
        "<tr><td>{}</td><td><strong>{}</strong></td></tr>".format(k, v)
        for k, v in stats
    )

    return """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Triangle Geometry Check — P2 Basket FEM</title>
<script src="{cdn}"></script>
<style>
  body   {{ font-family: Arial, sans-serif; margin: 24px; background: #f7f7f7; }}
  h1     {{ color: #2c3e50; margin-bottom: 4px; }}
  .sub   {{ color: #666; font-size: 0.92em; margin-top: 0; }}
  h2     {{ color: #34495e; border-bottom: 1px solid #ccc;
            padding-bottom: 4px; margin-top: 32px; }}
  .box   {{ background: white; border-radius: 6px;
            box-shadow: 0 1px 4px #ccc; padding: 14px; margin-bottom: 24px; }}
  table  {{ border-collapse: collapse; }}
  td     {{ border: 1px solid #ddd; padding: 6px 14px; }}
  tr:nth-child(even) {{ background: #f4f4f4; }}
  .note  {{ font-size: 0.83em; color: #888; margin-top: 8px; }}
  .ok    {{ color: #27ae60; font-weight: bold; }}
</style>
</head>
<body>

<h1>Equilateral Triangle — Geometry Check</h1>
<p class="sub">
  Side = {side:.1f} mm &nbsp;|&nbsp;
  Height = {height:.1f} mm &nbsp;|&nbsp;
  Area = {area:,.0f} mm² &nbsp;|&nbsp;
  Pillar pitch = {pitch:.0f} mm &nbsp;|&nbsp;
  Wall border = {wall:.0f} mm
</p>

<h2>Boundary condition layout</h2>
<div id="geo" class="box" style="height:640px"></div>
<p class="note">
  <span style="color:#e74c3c">&#9632;</span> Red zone = wall BC (fully fixed, 1 mm border) &nbsp;|&nbsp;
  <span style="color:#2980b9">&#9632;</span> Blue dots = Dynamask pillar centres (3 mm pitch)
  &nbsp;|&nbsp; <span style="color:#e67e22">&#10010;</span> Orange cross = centroid
</p>

<h2>Parameters</h2>
<div class="box">
<table>{table_rows}</table>
<p class="note">
  The pillar centres shown are those that fall inside the active area
  (&gt; {wall:.0f} mm from every edge). In the FEM each pillar centre will snap
  nodes within {snap:.1f} mm radius to a fixed BC.
</p>
</div>

<h2>What to check</h2>
<div class="box">
<ul>
  <li>Pillar grid covers the full active area uniformly — no gaps or clusters.</li>
  <li>Wall zone is continuous around all three edges.</li>
  <li>Centroid is at the geometric centre of the triangle (it will be the point of
      maximum out-of-plane displacement in the FEM result if the mesh is uniform).</li>
  <li>Number of pillar centres is reasonable (~{n_pillars} expected for this area
      at 3 mm pitch).</li>
</ul>
<p class="ok">If all of the above look correct, proceed to run freecad_triangle_fem_setup.py.</p>
</div>

<script>
Plotly.newPlot("geo", {traces}, {{
  xaxis: {{
    title: "x (mm)",
    scaleanchor: "y",
    scaleratio: 1,
    zeroline: false,
  }},
  yaxis: {{
    title: "y (mm)",
    zeroline: false,
  }},
  legend: {{ x: 0.01, y: 0.99, bgcolor: "rgba(255,255,255,0.85)" }},
  margin: {{ t: 20, b: 60, l: 70, r: 20 }},
  paper_bgcolor: "white",
  plot_bgcolor: "#fafafa",
}}, {{responsive: true}});
</script>

</body>
</html>
""".format(
        cdn=PLOTLY_CDN,
        side=SIDE_MM,
        height=HEIGHT_MM,
        area=AREA_TARGET_MM2,
        pitch=PILLAR_PITCH_MM,
        wall=WALL_MM,
        snap=PILLAR_SNAP_MM,
        table_rows=table_rows,
        n_pillars=len(pillars),
        traces=all_traces,
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)

    print("Equilateral triangle geometry check")
    print("  Side:     {:.2f} mm".format(SIDE_MM))
    print("  Height:   {:.2f} mm".format(HEIGHT_MM))
    print("  Area:     {:,.0f} mm²  (target {:,.0f} mm²)".format(
        math.sqrt(3) / 4 * SIDE_MM ** 2, AREA_TARGET_MM2))
    print("  Centroid: ({:.1f}, {:.1f})".format(*CENTROID))

    inner_v = inset_triangle(WALL_MM)
    if inner_v is None:
        print("ERROR: wall inset larger than triangle inradius.")
        return

    print("\nGenerating pillar grid ...")
    pillars = generate_pillar_centers()
    print("  Pillar centres inside active area: {}".format(len(pillars)))

    print("\nGenerating wall zone shading ...")
    wall_xs, wall_ys = wall_zone_points()
    print("  Wall zone sample points: {}".format(len(wall_xs)))

    html = build_html(pillars, wall_xs, wall_ys, inner_v)

    with open(OUTPUT_HTML, "w", encoding="utf-8") as f:
        f.write(html)

    print("\nSaved: {}".format(OUTPUT_HTML))
    print("Open in any browser to verify the geometry before running the FEM.")


main()