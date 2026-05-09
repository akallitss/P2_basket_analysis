
"""
Geometry verification script for p2_basket_detector_mechanical_tension_v4.FCStd
Reads the Gmsh UNV mesh embedded in the FCStd ZIP, extracts node coordinates,
and writes a self-contained HTML with Plotly (CDN) showing:
  - Top-down XY view: scatter coloured by depth
  - Filled material overlay: convex-hull polygons per layer, stacked
  - Per-layer top-down panels (same XY scale)
  - Side view with annotated layer boundaries

Run via:
  "C:\Program Files\FreeCAD 0.20\bin\FreeCADCmd.exe" check_geometry.py
Or directly in any Python interpreter (standard library only).
"""

import zipfile
import json
import math
import os

FCSTD   = r"C:\Temp\fc_v4.zip"
OUT_HTML = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "geometry_check.html")

SAMPLE_EVERY = 5   # every Nth node for scatter (102k → ~20k)
HULL_SAMPLE  = 3   # every Nth node for hull computation

PILLAR_PITCH_MM = 3.0    # centre-to-centre spacing
PILLAR_DIA_MM   = 0.4    # 400 µm Dynamask pillars
WALL_MM         = 1.0    # peripheral mesh-frame wall thickness

# ---------------------------------------------------------------------------
# Layer definitions
# z boundaries from UNV: -0.4, -0.2, 0.0, +0.128, +0.628  (all mm)
# Model v4 vs. target for new sweep model
# ---------------------------------------------------------------------------
LAYER_BANDS = [
    # (name, z_lo, z_hi, fill_color, rgba_fill, note_model, note_target)
    ("PCB substrate",
     -0.45, -0.01, "#8B4513", "rgba(139,69,19,0.35)",
     "400 µm (2×200)", "200 µm  ⚠ fix in sweep"),
    ("Amplification gap",
      0.01,  0.15, "#1E90FF", "rgba(30,144,255,0.40)",
     "128 µm", "150 µm  ⚠ fix in sweep"),
    ("Mesh layer",
      0.15,  0.65, "#2E8B57", "rgba(46,139,87,0.45)",
     "500 µm (effective)", "500 µm  ✓"),
]

# ---------------------------------------------------------------------------
# UNV node parser
# ---------------------------------------------------------------------------
def parse_unv_nodes(fcstd_path):
    nodes = {}
    with zipfile.ZipFile(fcstd_path) as z:
        with z.open("FemMesh.unv") as f:
            in_section = False
            want_coord = False
            current_id = None
            for raw in f:
                line = raw.decode("utf-8", errors="replace").rstrip()
                s = line.strip()
                if s == "-1":
                    in_section = False
                    want_coord = False
                    continue
                if s == "2411":
                    in_section = True
                    want_coord = False
                    continue
                if not in_section:
                    continue
                if not want_coord:
                    parts = s.split()
                    if parts and parts[0].isdigit():
                        current_id = int(parts[0])
                        want_coord = True
                    continue
                parts = line.replace("D", "E").replace("d", "e").split()
                if len(parts) == 3:
                    try:
                        x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
                        nodes[current_id] = (x, y, z)
                    except ValueError:
                        pass
                want_coord = False
    return nodes

# ---------------------------------------------------------------------------
# Convex hull  (Andrew's monotone chain, pure Python, O(n log n))
# ---------------------------------------------------------------------------
def convex_hull_2d(xy_pairs):
    pts = sorted(set(xy_pairs))
    if len(pts) < 3:
        return pts

    def cross(O, A, B):
        return (A[0]-O[0])*(B[1]-O[1]) - (A[1]-O[1])*(B[0]-O[0])

    lower = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper = []
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]

# ---------------------------------------------------------------------------
# Trace builders
# ---------------------------------------------------------------------------
def build_scatter_all(nodes, sample):
    """Scatter of all nodes, colour = Z depth."""
    DEPTH_COLORS = {
        "PCB substrate":      "#CD853F",
        "Amplification gap":  "#87CEEB",
        "Mesh layer":         "#2E8B57",
    }
    def band_of(z):
        for name, z_lo, z_hi, *_ in LAYER_BANDS:
            if z_lo <= z <= z_hi:
                return name
        return "other"

    buckets = {}
    for i, (nid, (x, y, z)) in enumerate(nodes.items()):
        if i % sample != 0:
            continue
        lbl = band_of(z)
        if lbl not in buckets:
            buckets[lbl] = {"x": [], "y": [], "text": []}
        buckets[lbl]["x"].append(round(x, 2))
        buckets[lbl]["y"].append(round(y, 2))
        buckets[lbl]["text"].append("z={:.3f}".format(z))

    traces = []
    for name, *_ in LAYER_BANDS:
        if name not in buckets:
            continue
        d = buckets[name]
        color = DEPTH_COLORS.get(name, "#999")
        traces.append({
            "type": "scatter", "mode": "markers", "name": name,
            "x": d["x"], "y": d["y"], "text": d["text"],
            "hoverinfo": "text",
            "marker": {"size": 2, "color": color, "opacity": 0.55},
        })
    return traces


def build_hull_traces(nodes, sample, show_legend=True):
    """One filled convex-hull polygon per layer band."""
    traces = []
    for name, z_lo, z_hi, line_color, fill_color, *_ in LAYER_BANDS:
        xy = []
        for i, (nid, (x, y, z)) in enumerate(nodes.items()):
            if i % sample != 0:
                continue
            if z_lo <= z <= z_hi:
                xy.append((round(x, 1), round(y, 1)))
        hull = convex_hull_2d(xy)
        if len(hull) < 3:
            continue
        hx = [p[0] for p in hull] + [hull[0][0]]
        hy = [p[1] for p in hull] + [hull[0][1]]
        traces.append({
            "type": "scatter", "mode": "lines",
            "name": name,
            "x": hx, "y": hy,
            "fill": "toself", "fillcolor": fill_color,
            "line": {"color": line_color, "width": 1.5},
            "showlegend": show_legend,
            "hoverinfo": "name",
        })
    return traces


def build_scatter_band(nodes, z_lo, z_hi, color, sample):
    xs, ys, ts = [], [], []
    for i, (nid, (x, y, z)) in enumerate(nodes.items()):
        if i % sample != 0:
            continue
        if z_lo <= z <= z_hi:
            xs.append(round(x, 2))
            ys.append(round(y, 2))
            ts.append("z={:.3f}".format(z))
    return [{"type": "scatter", "mode": "markers", "showlegend": False,
             "x": xs, "y": ys, "text": ts, "hoverinfo": "text",
             "marker": {"size": 2, "color": color, "opacity": 0.5}}]


def build_side_traces(nodes, sample):
    COLORS = {"PCB substrate": "#CD853F",
              "Amplification gap": "#87CEEB",
              "Mesh layer": "#2E8B57"}
    def band_of(z):
        for name, z_lo, z_hi, *_ in LAYER_BANDS:
            if z_lo <= z <= z_hi:
                return name
        return "other"
    buckets = {}
    for i, (nid, (x, y, z)) in enumerate(nodes.items()):
        if i % sample != 0:
            continue
        lbl = band_of(z)
        if lbl not in buckets:
            buckets[lbl] = {"x": [], "y": []}
        buckets[lbl]["x"].append(round(x, 2))
        buckets[lbl]["y"].append(round(z, 3))
    traces = []
    for name, *_ in LAYER_BANDS:
        if name not in buckets:
            continue
        d = buckets[name]
        traces.append({
            "type": "scatter", "mode": "markers", "showlegend": False,
            "name": name, "x": d["x"], "y": d["y"], "hoverinfo": "x+y",
            "marker": {"size": 2, "color": COLORS.get(name, "#999"), "opacity": 0.4},
        })
    return traces


def unique_z(nodes):
    zs = set(round(v[2], 3) for v in nodes.values())
    return sorted(zs)


def get_layer_hull(nodes, z_lo, z_hi, sample=HULL_SAMPLE):
    """Return convex hull vertex list [(x,y),...] for nodes in a z band."""
    xy = []
    for i, (nid, (x, y, z)) in enumerate(nodes.items()):
        if i % sample != 0:
            continue
        if z_lo <= z <= z_hi:
            xy.append((round(x, 1), round(y, 1)))
    return convex_hull_2d(xy)


def point_in_polygon(px, py, poly):
    """Ray-casting point-in-polygon test."""
    inside = False
    n = len(poly)
    j = n - 1
    for i in range(n):
        xi, yi = poly[i]
        xj, yj = poly[j]
        if ((yi > py) != (yj > py)) and (px < (xj - xi) * (py - yi) / (yj - yi) + xi):
            inside = not inside
        j = i
    return inside


def shrink_convex_polygon(hull, amount):
    """
    Shrink a convex polygon inward by `amount` mm.
    Each vertex is moved toward the centroid by `amount` along the centroid ray.
    Accurate enough for large polygons with small offsets.
    """
    cx = sum(p[0] for p in hull) / len(hull)
    cy = sum(p[1] for p in hull) / len(hull)
    shrunk = []
    for x, y in hull:
        dx, dy = x - cx, y - cy
        dist = math.sqrt(dx * dx + dy * dy)
        if dist > amount:
            scale = (dist - amount) / dist
            shrunk.append((round(cx + dx * scale, 3), round(cy + dy * scale, 3)))
    return shrunk

# ---------------------------------------------------------------------------
# Boundary-condition helpers
# ---------------------------------------------------------------------------

def gen_pillar_centers(xmin, xmax, ymin, ymax, polygon=None):
    """
    Return (x, y) pillar-centre coordinates on the 3 mm grid.
    If polygon is given (list of (x,y) vertices), only keep points inside it.
    The grid origin is fixed at (xmin+WALL, ymin+WALL) + half-pitch so it
    stays consistent between the full view and the detail zoom.
    """
    centers = []
    x = xmin + WALL_MM + PILLAR_PITCH_MM / 2
    while x <= xmax - WALL_MM:
        y = ymin + WALL_MM + PILLAR_PITCH_MM / 2
        while y <= ymax - WALL_MM:
            if polygon is None or point_in_polygon(x, y, polygon):
                centers.append((round(x, 2), round(y, 2)))
            y += PILLAR_PITCH_MM
        x += PILLAR_PITCH_MM
    return centers


def build_bc_traces(xmin, xmax, ymin, ymax, hull_outer=None, hull_inner=None):
    """
    Full-mesh BC overview: peripheral wall polygon + pillar grid.
    hull_outer: convex hull of the mesh layer (drawn as wall BC boundary).
    hull_inner: hull shrunk by WALL_MM (defines valid pillar area).
    """
    centers = gen_pillar_centers(xmin, xmax, ymin, ymax, polygon=hull_inner)

    pillar_trace = {
        "type": "scatter", "mode": "markers",
        "name": "Pillars — Uz=0  ({} pts, {} mm pitch, {} µm dia)".format(
            len(centers), int(PILLAR_PITCH_MM), int(PILLAR_DIA_MM * 1000)),
        "x": [c[0] for c in centers],
        "y": [c[1] for c in centers],
        "marker": {"size": 3, "color": "#FF8C00", "opacity": 0.55},
        "hovertemplate": "x=%{x:.1f} mm, y=%{y:.1f} mm<extra>Pillar BC (Uz=0)</extra>",
    }

    if hull_outer is not None:
        wx = [p[0] for p in hull_outer] + [hull_outer[0][0]]
        wy = [p[1] for p in hull_outer] + [hull_outer[0][1]]
    else:
        r = WALL_MM
        wx = [xmin + r, xmax - r, xmax - r, xmin + r, xmin + r]
        wy = [ymin + r, ymin + r, ymax - r, ymax - r, ymin + r]

    wall_trace = {
        "type": "scatter", "mode": "lines",
        "name": "Mesh wall — Ux=Uy=Uz=0  (1 mm frame)",
        "x": wx, "y": wy,
        "fill": "toself", "fillcolor": "rgba(220,20,60,0.07)",
        "line": {"color": "#DC143C", "width": 2.5},
        "hoverinfo": "name",
    }
    return [pillar_trace, wall_trace]


def build_detail_data(cx, cy, xmin_global, ymin_global, half=10.5, polygon=None):
    """
    Returns (traces, shapes) for a zoomed 21×21 mm window centred on (cx, cy).
    Uses the same grid origin as gen_pillar_centers so positions match.
    polygon: if given, only pillars inside it are drawn (matches full-view filter).
    Circles drawn at actual 400 µm scale as Plotly shapes.
    """
    r = PILLAR_DIA_MM / 2
    pitch = PILLAR_PITCH_MM
    x0 = xmin_global + WALL_MM + pitch / 2
    y0 = ymin_global + WALL_MM + pitch / 2
    nx0 = math.ceil((cx - half - x0) / pitch)
    ny0 = math.ceil((cy - half - y0) / pitch)

    px, py, shapes = [], [], []
    nx = nx0
    while True:
        x = round(x0 + nx * pitch, 4)
        if x > cx + half:
            break
        ny = ny0
        while True:
            y = round(y0 + ny * pitch, 4)
            if y > cy + half:
                break
            in_window = (cx - half) <= x <= (cx + half) and (cy - half) <= y <= (cy + half)
            in_poly = polygon is None or point_in_polygon(x, y, polygon)
            if in_window and in_poly:
                px.append(x)
                py.append(y)
                shapes.append({
                    "type": "circle", "xref": "x", "yref": "y",
                    "x0": round(x - r, 4), "y0": round(y - r, 4),
                    "x1": round(x + r, 4), "y1": round(y + r, 4),
                    "fillcolor": "rgba(255,140,0,0.35)",
                    "line": {"color": "#CC6600", "width": 1.5},
                })
            ny += 1
        nx += 1

    traces = [{
        "type": "scatter", "mode": "markers",
        "name": "Pillar centres",
        "x": px, "y": py,
        "marker": {"size": 5, "color": "#CC6600", "symbol": "cross"},
        "hovertemplate": "x=%{x:.2f}, y=%{y:.2f} mm<extra>Pillar centre</extra>",
    }]
    return traces, shapes


# ---------------------------------------------------------------------------
# HTML
# ---------------------------------------------------------------------------
HTML = """\
<!DOCTYPE html><html>
<head>
<meta charset="utf-8">
<title>P2 Basket – Geometry Check</title>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
<style>
  body{{font-family:sans-serif;margin:20px;background:#f5f5f5}}
  h2{{color:#222;margin:6px 0}} h3{{color:#444;margin:14px 0 6px}}
  .row{{display:flex;gap:16px;flex-wrap:wrap;align-items:flex-start}}
  .panel{{background:#fff;border-radius:8px;padding:12px;
          box-shadow:0 2px 6px rgba(0,0,0,.12)}}
  table{{border-collapse:collapse;font-size:12px}}
  td,th{{border:1px solid #ccc;padding:3px 8px}}
  th{{background:#eee}} .warn{{color:#c00;font-weight:bold}}
  .ok{{color:#060}} hr{{border:none;border-top:1px solid #ddd;margin:16px 0}}
</style>
</head>
<body>
<h2>P2 Basket Bulk Micromegas – Geometry Verification</h2>

<div class="row"><div class="panel">
  <b>Layer stack — v4 model vs. sweep target</b>
  <table>
    <tr><th>Layer</th><th>Z (mm)</th><th>Model (µm)</th>
        <th>Target (µm)</th><th>Note</th></tr>
    {layer_rows}
  </table>
  <br>Nodes: {node_count} &nbsp;|&nbsp; Plotted 1/{sample}: {plotted}
  &nbsp;|&nbsp; XY extent: x={xmin:.0f}…{xmax:.0f} mm, y={ymin:.0f}…{ymax:.0f} mm
</div></div><hr>

<!-- ── Row 1: scatter top-down + side view ── -->
<h3>Scatter view — nodes coloured by layer</h3>
<div class="row">
  <div class="panel">
    <b>Top-down (XY)</b>
    <div id="p_xy" style="width:620px;height:560px"></div>
  </div>
  <div class="panel">
    <b>Side view (X vs Z)</b>
    <div id="p_xz" style="width:360px;height:560px"></div>
  </div>
</div><hr>

<!-- ── Row 2: filled material overlay + per-layer panels ── -->
<h3>Filled material view — convex hull per layer, same scale</h3>
<div class="row">

  <div class="panel">
    <b>All layers overlaid</b><br>
    <span style="font-size:11px;color:#666">
      Green = mesh &nbsp;|&nbsp; Blue = gap &nbsp;|&nbsp; Brown = PCB</span>
    <div id="p_hull" style="width:420px;height:520px"></div>
  </div>

  <div class="panel">
    <b>▲ Mesh layer</b>
    <div class="layer-label" style="font-size:11px;color:#2E8B57">
      z = +0.128 → +0.628 mm (500 µm ✓)</div>
    <div id="p_l2" style="width:300px;height:300px"></div>
  </div>

  <div class="panel">
    <b>── Amplification gap</b>
    <div style="font-size:11px;color:#1E90FF">
      z = 0 → +0.128 mm (model 128 µm → target 150 µm)</div>
    <div id="p_l1" style="width:300px;height:300px"></div>
  </div>

  <div class="panel">
    <b>▼ PCB substrate</b>
    <div style="font-size:11px;color:#8B4513">
      z = −0.4 → 0 mm (model 400 µm → target 200 µm)</div>
    <div id="p_l0" style="width:300px;height:300px"></div>
  </div>

</div>

<hr>
<h3>Boundary Conditions — mesh wall &amp; Dynamask pillar supports</h3>
<div class="row">
  <div class="panel">
    <b>Full mesh overview</b><br>
    <span style="font-size:11px;color:#666">
      <span style="color:#DC143C;font-weight:bold">Red border</span>
      = peripheral wall BC (Ux=Uy=Uz=0, 1 mm frame)&nbsp;|&nbsp;
      <span style="color:#FF8C00;font-weight:bold">Orange dots</span>
      = pillar centres (Uz=0 only), {n_pillars} pillars, {pillar_pitch} mm pitch
    </span>
    <div id="p_bc_full" style="width:620px;height:560px"></div>
  </div>
  <div class="panel">
    <b>Detail — 21×21 mm around mesh centroid (to scale)</b><br>
    <span style="font-size:11px;color:#666">
      Circles = Dynamask pillars at actual 400 µm diameter &amp; 3 mm pitch
    </span>
    <div id="p_bc_det" style="width:420px;height:420px"></div>
  </div>
</div>

<script>
var scatter_xy = {scatter_xy};
var scatter_xz = {scatter_xz};
var hull_all   = {hull_all};
var hull_l0    = {hull_l0};
var hull_l1    = {hull_l1};
var hull_l2    = {hull_l2};

var XRNG = [{xmin:.0f}-15, {xmax:.0f}+15];
var YRNG = [{ymin:.0f}-15, {ymax:.0f}+15];
var BASE = {{
  xaxis:{{title:'X (mm)',scaleanchor:'y',scaleratio:1,range:XRNG}},
  yaxis:{{title:'Y (mm)',range:YRNG}},
  margin:{{t:10,r:10,b:40,l:50}},
}};

Plotly.newPlot('p_xy', scatter_xy, Object.assign({{}}, BASE, {{
  legend:{{orientation:'h', y:-0.15, font:{{size:10}}}},
}}));

Plotly.newPlot('p_hull', hull_all, Object.assign({{}}, BASE, {{
  legend:{{orientation:'h', y:-0.12, font:{{size:10}}}},
}}));

Plotly.newPlot('p_l0', hull_l0, Object.assign({{}}, BASE, {{margin:{{t:5,r:5,b:30,l:40}}}}));
Plotly.newPlot('p_l1', hull_l1, Object.assign({{}}, BASE, {{margin:{{t:5,r:5,b:30,l:40}}}}));
Plotly.newPlot('p_l2', hull_l2, Object.assign({{}}, BASE, {{margin:{{t:5,r:5,b:30,l:40}}}}));

Plotly.newPlot('p_xz', scatter_xz, {{
  xaxis:{{title:'X (mm)'}},
  yaxis:{{title:'Z (mm)',range:[-0.55,0.80]}},
  margin:{{t:10,r:130,b:40,l:50}},
  shapes:[
    {{type:'line',x0:0,x1:1,xref:'paper',y0:-0.4, y1:-0.4, line:{{color:'#8B4513',dash:'dot',width:1.5}}}},
    {{type:'line',x0:0,x1:1,xref:'paper',y0:-0.2, y1:-0.2, line:{{color:'#CD853F',dash:'dot',width:1.5}}}},
    {{type:'line',x0:0,x1:1,xref:'paper',y0: 0.0, y1: 0.0, line:{{color:'#2d8a2d',dash:'dot',width:1.5}}}},
    {{type:'line',x0:0,x1:1,xref:'paper',y0: 0.128,y1:0.128,line:{{color:'#1E90FF',dash:'dot',width:1.5}}}},
    {{type:'line',x0:0,x1:1,xref:'paper',y0: 0.628,y1:0.628,line:{{color:'#2E8B57',dash:'dot',width:1.5}}}},
  ],
  annotations:[
    {{x:1.02,y:-0.4,  xref:'paper',yref:'y',text:'PCB base −0.4', showarrow:false,xanchor:'left',font:{{size:9}}}},
    {{x:1.02,y:-0.2,  xref:'paper',yref:'y',text:'PCB mid −0.2',  showarrow:false,xanchor:'left',font:{{size:9}}}},
    {{x:1.02,y: 0.0,  xref:'paper',yref:'y',text:'gap btm 0',     showarrow:false,xanchor:'left',font:{{size:9}}}},
    {{x:1.02,y: 0.128,xref:'paper',yref:'y',text:'gap top +0.128',showarrow:false,xanchor:'left',font:{{size:9}}}},
    {{x:1.02,y: 0.628,xref:'paper',yref:'y',text:'mesh top +0.628',showarrow:false,xanchor:'left',font:{{size:9}}}},
  ],
}});

// ── Boundary conditions ──────────────────────────────────────────
var bc_traces  = {bc_traces};
var det_traces = {det_traces};
var det_shapes = {det_shapes};

Plotly.newPlot('p_bc_full', bc_traces, Object.assign({{}}, BASE, {{
  title:{{text:'Full mesh — {n_pillars} pillars + peripheral wall', font:{{size:13}}}},
  legend:{{orientation:'h', y:-0.15, font:{{size:10}}}},
  margin:{{t:40,r:10,b:60,l:50}},
}}));

Plotly.newPlot('p_bc_det', det_traces, {{
  xaxis:{{title:'X (mm)', range:[{det_cx_lo:.2f}, {det_cx_hi:.2f}],
          scaleanchor:'y', scaleratio:1}},
  yaxis:{{title:'Y (mm)', range:[{det_cy_lo:.2f}, {det_cy_hi:.2f}]}},
  shapes: det_shapes,
  annotations:[
    {{x:{det_cx:.2f}, y:{det_cy_lo:.2f}+1.0, xref:'x', yref:'y',
      text:'← 3 mm →', showarrow:false, font:{{size:10, color:'#555'}}}},
    {{x:{det_cx_lo:.2f}+0.5, y:{det_cy:.2f}, xref:'x', yref:'y',
      text:'400 µm Ø', showarrow:true, ax:30, ay:0,
      font:{{size:9, color:'#CC6600'}}}},
  ],
  legend:{{orientation:'h', y:-0.15, font:{{size:10}}}},
  margin:{{t:10,r:10,b:60,l:55}},
  paper_bgcolor:'white', plot_bgcolor:'white',
}});
</script>
</body></html>
"""

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def _main():
    log_path = OUT_HTML.replace(".html", ".log")
    log = open(log_path, "w", encoding="utf-8")

    def info(msg):
        log.write(msg + "\n"); log.flush()
        try:
            import sys; print(msg); sys.stdout.flush()
        except Exception:
            pass

    try:
        info("Reading nodes ...")
        nodes = parse_unv_nodes(FCSTD)
        info("  {} nodes".format(len(nodes)))
        info("  Z values: {}".format(unique_z(nodes)))

        xs = [v[0] for v in nodes.values()]
        ys = [v[1] for v in nodes.values()]
        xmin, xmax = min(xs), max(xs)
        ymin, ymax = min(ys), max(ys)

        info("Building scatter traces ...")
        scatter_xy = build_scatter_all(nodes, SAMPLE_EVERY)
        scatter_xz = build_side_traces(nodes, SAMPLE_EVERY)
        plotted = sum(len(t["x"]) for t in scatter_xy)

        info("Computing convex hulls ...")
        hull_all = build_hull_traces(nodes, HULL_SAMPLE, show_legend=True)
        hull_l0  = build_hull_traces(nodes, HULL_SAMPLE, show_legend=False)[:1]
        hull_l1  = build_hull_traces(nodes, HULL_SAMPLE, show_legend=False)[1:2]
        hull_l2  = build_hull_traces(nodes, HULL_SAMPLE, show_legend=False)[2:3]
        info("  Hulls: {} traces".format(len(hull_all)))

        layer_rows = "\n".join(
            "<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>"
            .format(name, "{:.3f} → {:.3f}".format(z_lo, z_hi), m, t, note)
            for (name, z_lo, z_hi, _, __, m, t, note)
            in [b + (b[6],) for b in
                [(n,a,b,c,d,e,f) for n,a,b,c,d,e,f in LAYER_BANDS]]
        )

        # simpler row build
        rows = []
        for name, z_lo, z_hi, _, __, model_t, target_t in LAYER_BANDS:
            rows.append(
                "<tr><td>{}</td><td>{:.3f} → {:.3f}</td>"
                "<td>{}</td><td>{}</td></tr>".format(
                    name, z_lo, z_hi, model_t, target_t)
            )
        layer_rows = "\n".join(rows)

        info("Extracting mesh layer hull ...")
        hull_outer = get_layer_hull(nodes, 0.10, 0.70, sample=HULL_SAMPLE)
        hull_inner = shrink_convex_polygon(hull_outer, WALL_MM)
        info("  Hull vertices: {} outer, {} inner".format(
            len(hull_outer), len(hull_inner)))

        info("Building BC constraint traces ...")
        bc_traces = build_bc_traces(xmin, xmax, ymin, ymax,
                                    hull_outer=hull_outer, hull_inner=hull_inner)
        n_pillars = len(bc_traces[0]["x"])

        det_cx = sum(p[0] for p in hull_outer) / len(hull_outer)
        det_cy = sum(p[1] for p in hull_outer) / len(hull_outer)
        det_half = 10.5
        det_traces, det_shapes = build_detail_data(
            det_cx, det_cy, xmin, ymin, det_half, polygon=hull_inner)
        info("  {} pillar BCs, {} detail circles".format(
            n_pillars, len(det_shapes)))

        html = HTML.format(
            scatter_xy=json.dumps(scatter_xy),
            scatter_xz=json.dumps(scatter_xz),
            hull_all=json.dumps(hull_all),
            hull_l0=json.dumps(hull_l0),
            hull_l1=json.dumps(hull_l1),
            hull_l2=json.dumps(hull_l2),
            bc_traces=json.dumps(bc_traces),
            det_traces=json.dumps(det_traces),
            det_shapes=json.dumps(det_shapes),
            layer_rows=layer_rows,
            node_count=len(nodes),
            sample=SAMPLE_EVERY,
            plotted=plotted,
            xmin=xmin, xmax=xmax,
            ymin=ymin, ymax=ymax,
            n_pillars=n_pillars,
            pillar_pitch=int(PILLAR_PITCH_MM),
            det_cx=det_cx,  det_cy=det_cy,
            det_cx_lo=det_cx - det_half, det_cx_hi=det_cx + det_half,
            det_cy_lo=det_cy - det_half, det_cy_hi=det_cy + det_half,
        )

        with open(OUT_HTML, "w", encoding="utf-8") as f:
            f.write(html)

        info("SUCCESS -> " + OUT_HTML)

    except Exception as exc:
        import traceback
        info("ERROR: " + str(exc))
        info(traceback.format_exc())

    log.close()


_main()
