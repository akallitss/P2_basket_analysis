"""
visualize_triangle_sweep.py  –  Interactive 3D Uz surface for triangle FEM results

Usage (standard Python, no FreeCAD needed):
    python FreeCAD/visualize_triangle_sweep.py

Reads:
    FreeCAD/results/triangle_simulation/triangle_sweep_T????Nm.csv

Writes:
    FreeCAD/results/triangle_simulation/triangle_uz_surface.html
"""

import os
import glob
import csv
import re
import numpy as np
import matplotlib.tri as mtri

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(SCRIPT_DIR, "results", "triangle_simulation")
OUTPUT_HTML = os.path.join(RESULTS_DIR, "triangle_uz_surface.html")

GRID_N = 200   # grid resolution per axis


def load_csv(path):
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            rows.append({k: float(v) for k, v in row.items()})
    return rows


def midplane_data(rows):
    """Average Uz over both z-layers to get one value per (x, y) position."""
    by_xy = {}
    for r in rows:
        key = (round(r["x_mm"], 3), round(r["y_mm"], 3))
        by_xy.setdefault(key, []).append(r["Uz_mm"])
    xs = np.array([k[0] for k in by_xy])
    ys = np.array([k[1] for k in by_xy])
    uzs = np.array([sum(v) / len(v) for v in by_xy.values()]) * 1000.0  # → µm
    return xs, ys, uzs


def interpolate_to_grid(xs, ys, zs, n=GRID_N):
    """
    Delaunay-triangulate the scatter data (matplotlib.tri), then linearly
    interpolate onto a regular grid.  Points outside the convex hull get NaN.
    """
    triang = mtri.Triangulation(xs, ys)
    interp = mtri.LinearTriInterpolator(triang, zs)

    xi = np.linspace(xs.min(), xs.max(), n)
    yi = np.linspace(ys.min(), ys.max(), n)
    Xi, Yi = np.meshgrid(xi, yi)
    Zi = interp(Xi, Yi)          # masked array — NaN outside hull
    Zi = np.ma.filled(Zi, np.nan)
    return xi.tolist(), yi.tolist(), Zi.tolist()


def load_all_steps():
    pattern = os.path.join(RESULTS_DIR, "triangle_sweep_T*Nm.csv")
    paths = sorted(glob.glob(pattern))
    steps = []
    for p in paths:
        m = re.search(r"T(\d+)Nm", p)
        T_Nm = int(m.group(1)) if m else 0
        rows = load_csv(p)
        xs, ys, uzs_um = midplane_data(rows)
        print("  T={} N/m: {} unique xy nodes, max|Uz|={:.2f} µm".format(
            T_Nm, len(xs), float(np.max(np.abs(uzs_um)))))
        x_ax, y_ax, grid = interpolate_to_grid(xs, ys, uzs_um)
        steps.append((T_Nm, x_ax, y_ax, grid, float(np.max(np.abs(uzs_um)))))
    return steps


# ---------------------------------------------------------------------------
# JSON helpers (avoid json module dependency on large arrays)
# ---------------------------------------------------------------------------

def js_flat(lst):
    return "[" + ",".join("{:.4f}".format(v) for v in lst) + "]"


def js_2d(grid):
    def fmt(v):
        return "null" if (v != v) else "{:.4f}".format(v)   # NaN check

    rows = ["[" + ",".join(fmt(v) for v in row) + "]" for row in grid]
    return "[" + ",".join(rows) + "]"


# ---------------------------------------------------------------------------
# HTML builder
# ---------------------------------------------------------------------------

def build_html(steps):
    traces_js = []
    buttons = []

    for i, (T_Nm, x_ax, y_ax, grid, max_uz) in enumerate(steps):
        label = "{} N/m ({} N/cm)".format(T_Nm, T_Nm // 100)
        visible = "true" if i == 0 else "false"

        trace = """{{
  type: 'surface',
  x: {x},
  y: {y},
  z: {z},
  colorscale: 'RdBu',
  reversescale: true,
  colorbar: {{title: {{text: 'Uz (µm)', side: 'right'}}}},
  visible: {vis},
  name: '{lbl}',
  hovertemplate: 'X: %{{x:.1f}} mm<br>Y: %{{y:.1f}} mm<br>Uz: %{{z:.3f}} µm<extra></extra>'
}}""".format(
            x=js_flat(x_ax),
            y=js_flat(y_ax),
            z=js_2d(grid),
            vis=visible,
            lbl=label,
        )
        traces_js.append(trace)

        vis_arr = ["false"] * len(steps)
        vis_arr[i] = "true"
        vis_str = "[" + ",".join(vis_arr) + "]"
        buttons.append("""{{
  label: '{lbl}',
  method: 'update',
  args: [
    {{visible: {vis}}},
    {{title: 'Out-of-plane displacement Uz  ·  {lbl}  ·  max |Uz| = {mx:.2f} µm'}}
  ]
}}""".format(lbl=label, vis=vis_str, mx=max_uz))

    traces_str = "[" + ",\n".join(traces_js) + "]"
    buttons_str = "[" + ",\n".join(buttons) + "]"

    first_T, _, _, _, first_max = steps[0]
    first_label = "{} N/m ({} N/cm)".format(first_T, first_T // 100)

    return """<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>Triangle FEM — Uz surface</title>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
<style>
  body  {{ margin:0; background:#111; color:#eee; font-family:sans-serif; }}
  #hdr  {{ padding:14px 22px 2px; }}
  #hdr h2 {{ margin:0 0 3px; font-size:1.15rem; color:#ccc; }}
  #hdr p  {{ margin:0; font-size:0.82rem; color:#777; }}
  #plot {{ width:100vw; height:calc(100vh - 72px); }}
</style>
</head>
<body>
<div id="hdr">
  <h2>Triangle FEM — Out-of-plane displacement U<sub>z</sub></h2>
  <p>Equilateral triangle · area = 169 256 mm² · E = 1 GPa · t = 0.5 mm · pillar pitch = 3 mm · wall BC = 1 mm</p>
</div>
<div id="plot"></div>
<script>
var traces = {traces};

var layout = {{
  title: {{
    text: 'Out-of-plane displacement Uz  ·  {first_lbl}  ·  max |Uz| = {first_max:.2f} µm',
    font: {{ color: '#ccc', size: 13 }}
  }},
  paper_bgcolor: '#111',
  font: {{ color: '#ccc' }},
  scene: {{
    xaxis: {{ title: 'X (mm)', color: '#999', gridcolor: '#333', zerolinecolor: '#555' }},
    yaxis: {{ title: 'Y (mm)', color: '#999', gridcolor: '#333', zerolinecolor: '#555' }},
    zaxis: {{ title: 'Uz (µm)', color: '#999', gridcolor: '#333', zerolinecolor: '#555' }},
    bgcolor: '#0a0a0a',
    camera: {{ eye: {{ x: 0.7, y: -1.5, z: 1.1 }} }},
    aspectmode: 'manual',
    aspectratio: {{ x: 1.2, y: 1.0, z: 0.4 }}
  }},
  updatemenus: [{{
    type: 'buttons',
    direction: 'right',
    x: 0.01, y: 1.06,
    xanchor: 'left', yanchor: 'top',
    bgcolor: '#1e1e1e',
    bordercolor: '#444',
    font: {{ color: '#ddd', size: 12 }},
    buttons: {buttons}
  }}],
  margin: {{ l: 0, r: 0, t: 70, b: 0 }}
}};

Plotly.newPlot('plot', traces, layout, {{responsive: true}});
</script>
</body>
</html>
""".format(
        traces=traces_str,
        buttons=buttons_str,
        first_lbl=first_label,
        first_max=first_max,
    )


def main():
    print("Loading triangle FEM results ...")
    steps = load_all_steps()
    if not steps:
        print("No result CSVs found in", RESULTS_DIR)
        return

    print("Building HTML ...")
    html = build_html(steps)
    with open(OUTPUT_HTML, "w", encoding="utf-8") as f:
        f.write(html)
    print("Written:", OUTPUT_HTML)


if __name__ == "__main__":
    main()