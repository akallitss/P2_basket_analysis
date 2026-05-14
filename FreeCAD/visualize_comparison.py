"""
visualize_comparison.py  -  Side-by-side Uz surface: Triangle vs P2 Basket (Gmsh)

Usage (standard Python, no FreeCAD needed):
    python FreeCAD/visualize_comparison.py

Reads:
    FreeCAD/results/triangle_simulation/triangle_sweep_T????Nm.csv
    FreeCAD/results/basket_gmsh_simulation/basket_gmsh_T????Nm.csv

Writes:
    FreeCAD/results/comparison_basket_vs_triangle.html
"""

import os
import glob
import csv
import re
import numpy as np
import matplotlib.tri as mtri

SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
TRIANGLE_DIR = os.path.join(SCRIPT_DIR, "results", "triangle_simulation")
BASKET_DIR   = os.path.join(SCRIPT_DIR, "results", "basket_gmsh_simulation")
OUTPUT_HTML  = os.path.join(SCRIPT_DIR, "results", "comparison_basket_vs_triangle.html")

GRID_N = 150


# ---------------------------------------------------------------------------
# Data loading & interpolation
# ---------------------------------------------------------------------------

def load_csv(path):
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            rows.append({k: float(v) for k, v in row.items()})
    return rows


def midplane_data(rows):
    """Average Uz over z-layers; return scatter (x, y, Uz_um)."""
    by_xy = {}
    for r in rows:
        key = (round(r["x_mm"], 3), round(r["y_mm"], 3))
        by_xy.setdefault(key, []).append(r["Uz_mm"])
    xs  = np.array([k[0] for k in by_xy])
    ys  = np.array([k[1] for k in by_xy])
    uzs = np.array([sum(v) / len(v) for v in by_xy.values()]) * 1000.0  # -> um
    return xs, ys, uzs


def interpolate_to_grid(xs, ys, zs, n=GRID_N):
    triang = mtri.Triangulation(xs, ys)
    interp = mtri.LinearTriInterpolator(triang, zs)
    xi = np.linspace(xs.min(), xs.max(), n)
    yi = np.linspace(ys.min(), ys.max(), n)
    Xi, Yi = np.meshgrid(xi, yi)
    Zi = np.ma.filled(interp(Xi, Yi), np.nan)
    return xi.tolist(), yi.tolist(), Zi.tolist()


def load_dataset(directory, csv_prefix, label):
    pattern = os.path.join(directory, csv_prefix + "T*Nm.csv")
    paths   = sorted(glob.glob(pattern))
    steps   = {}
    for p in paths:
        m = re.search(r"T(\d+)Nm", p)
        T_Nm = int(m.group(1)) if m else 0
        rows = load_csv(p)
        xs, ys, uzs_um = midplane_data(rows)
        max_uz = float(np.max(np.abs(uzs_um)))
        print("  [{:>8}] T={:4d} N/m  {:5d} xy nodes  max|Uz|={:.2f} um".format(
            label, T_Nm, len(xs), max_uz))
        x_ax, y_ax, grid = interpolate_to_grid(xs, ys, uzs_um)
        steps[T_Nm] = (x_ax, y_ax, grid, max_uz)
    return steps


# ---------------------------------------------------------------------------
# JSON helpers
# ---------------------------------------------------------------------------

def js_flat(lst):
    return "[" + ",".join("{:.4f}".format(v) for v in lst) + "]"


def js_2d(grid):
    def fmt(v):
        return "null" if (v != v) else "{:.4f}".format(v)
    rows = ["[" + ",".join(fmt(v) for v in row) + "]" for row in grid]
    return "[" + ",".join(rows) + "]"


# ---------------------------------------------------------------------------
# HTML builder
# ---------------------------------------------------------------------------

def build_html(tri_steps, bsk_steps):
    tensions = sorted(set(tri_steps) & set(bsk_steps))
    if not tensions:
        raise RuntimeError("No matching tensions between the two datasets")

    traces_js = []
    buttons   = []
    n_t = len(tensions)

    for idx, T_Nm in enumerate(tensions):
        T_x, T_y, T_grid, T_max = tri_steps[T_Nm]
        B_x, B_y, B_grid, B_max = bsk_steps[T_Nm]

        # shared colorscale range across both geometries for this tension
        T_vals   = [v for row in T_grid for v in row if v == v]
        B_vals   = [v for row in B_grid for v in row if v == v]
        all_vals = T_vals + B_vals
        cmin_v   = min(all_vals)
        cmax_v   = max(all_vals)

        label   = "{} N/m ({} N/cm)".format(T_Nm, T_Nm // 100)
        visible = "true" if idx == 0 else "false"

        # trace 2*idx   : triangle (left  — scene)
        traces_js.append("""{{
  type: 'surface',
  scene: 'scene',
  x: {tx}, y: {ty}, z: {tz},
  colorscale: 'RdBu', reversescale: true,
  cmin: {cmin:.4f}, cmax: {cmax:.4f},
  showscale: false,
  visible: {vis},
  name: 'Triangle',
  hovertemplate: 'X: %{{x:.1f}} mm<br>Y: %{{y:.1f}} mm<br>Uz: %{{z:.3f}} um<extra>Triangle</extra>'
}}""".format(tx=js_flat(T_x), ty=js_flat(T_y), tz=js_2d(T_grid),
             cmin=cmin_v, cmax=cmax_v, vis=visible))

        # trace 2*idx+1 : basket  (right — scene2)
        traces_js.append("""{{
  type: 'surface',
  scene: 'scene2',
  x: {bx}, y: {by}, z: {bz},
  colorscale: 'RdBu', reversescale: true,
  cmin: {cmin:.4f}, cmax: {cmax:.4f},
  colorbar: {{title: {{text: 'Uz (um)', side: 'right'}}, x: 1.01, len: 0.65}},
  visible: {vis},
  name: 'Basket (Gmsh)',
  hovertemplate: 'X: %{{x:.1f}} mm<br>Y: %{{y:.1f}} mm<br>Uz: %{{z:.3f}} um<extra>Basket (Gmsh)</extra>'
}}""".format(bx=js_flat(B_x), by=js_flat(B_y), bz=js_2d(B_grid),
             cmin=cmin_v, cmax=cmax_v, vis=visible))

        # button: show only this tension's pair of traces
        vis_arr          = ["false"] * (2 * n_t)
        vis_arr[2 * idx]     = "true"
        vis_arr[2 * idx + 1] = "true"
        vis_str = "[" + ",".join(vis_arr) + "]"
        buttons.append("""{{
  label: '{lbl}',
  method: 'update',
  args: [
    {{visible: {vis}}},
    {{title: {{text: 'Uz comparison  |  {lbl}  |  Triangle max|Uz|={tmx:.2f} um  |  Basket max|Uz|={bmx:.2f} um'}}}}
  ]
}}""".format(lbl=label, vis=vis_str, tmx=T_max, bmx=B_max))

    traces_str  = "[\n" + ",\n".join(traces_js) + "\n]"
    buttons_str = "[\n" + ",\n".join(buttons) + "\n]"

    T0      = tensions[0]
    T0_lbl  = "{} N/m ({} N/cm)".format(T0, T0 // 100)
    T0_tmx  = tri_steps[T0][3]
    T0_bmx  = bsk_steps[T0][3]

    return """<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>Basket vs Triangle - Uz comparison</title>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
<style>
  body {{ margin:0; background:#111; color:#eee; font-family:sans-serif; }}
  #hdr {{ padding:12px 22px 2px; }}
  #hdr h2 {{ margin:0 0 2px; font-size:1.1rem; color:#ccc; }}
  #hdr p  {{ margin:0; font-size:0.80rem; color:#666; }}
  #plot   {{ width:100vw; height:calc(100vh - 68px); }}
</style>
</head>
<body>
<div id="hdr">
  <h2>Out-of-plane displacement Uz  |  Triangle (left) vs P2 Basket (right)  |  Both Gmsh C3D10</h2>
  <p>Same mesh settings: Gmsh C3D10, Lchar 2-5 mm, E=1 GPa, t=0.5 mm, pillar pitch=3 mm, wall BC=1 mm.
     Shared colorscale range per tension step.</p>
</div>
<div id="plot"></div>
<script>
var traces = {traces};

var scene_base = {{
  xaxis: {{ title: 'X (mm)', color: '#999', gridcolor: '#333', zerolinecolor: '#555' }},
  yaxis: {{ title: 'Y (mm)', color: '#999', gridcolor: '#333', zerolinecolor: '#555' }},
  zaxis: {{ title: 'Uz (um)', color: '#999', gridcolor: '#333', zerolinecolor: '#555' }},
  bgcolor: '#0a0a0a',
  camera: {{ eye: {{ x: 0.7, y: -1.5, z: 1.0 }} }},
  aspectmode: 'data'
}};

var layout = {{
  title: {{
    text: 'Uz comparison  |  {T0_lbl}  |  Triangle max|Uz|={T0_tmx:.2f} um  |  Basket max|Uz|={T0_bmx:.2f} um',
    font: {{ color: '#ccc', size: 12 }}
  }},
  paper_bgcolor: '#111',
  font: {{ color: '#ccc' }},

  scene:  Object.assign({{}}, scene_base, {{ domain: {{ x: [0.00, 0.47], y: [0.0, 1.0] }} }}),
  scene2: Object.assign({{}}, scene_base, {{ domain: {{ x: [0.53, 1.00], y: [0.0, 1.0] }} }}),

  annotations: [
    {{
      text: 'TRIANGLE', xref: 'paper', yref: 'paper',
      x: 0.235, y: 1.03, showarrow: false,
      font: {{ color: '#9bc', size: 13 }}
    }},
    {{
      text: 'P2 BASKET (Gmsh re-mesh)', xref: 'paper', yref: 'paper',
      x: 0.765, y: 1.03, showarrow: false,
      font: {{ color: '#c9b', size: 13 }}
    }}
  ],

  updatemenus: [{{
    type: 'buttons',
    direction: 'right',
    x: 0.01, y: 1.10,
    xanchor: 'left', yanchor: 'top',
    bgcolor: '#1e1e1e', bordercolor: '#444',
    font: {{ color: '#ddd', size: 12 }},
    buttons: {buttons}
  }}],

  margin: {{ l: 0, r: 70, t: 90, b: 0 }}
}};

Plotly.newPlot('plot', traces, layout, {{responsive: true}});
</script>
</body>
</html>
""".format(
        traces=traces_str,
        buttons=buttons_str,
        T0_lbl=T0_lbl,
        T0_tmx=T0_tmx,
        T0_bmx=T0_bmx,
    )


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    print("Loading triangle results ...")
    tri_steps = load_dataset(TRIANGLE_DIR, "triangle_sweep_", "triangle")

    print("Loading basket Gmsh results ...")
    bsk_steps = load_dataset(BASKET_DIR, "basket_gmsh_", "basket")

    if not tri_steps:
        print("No triangle CSVs found in", TRIANGLE_DIR)
        return
    if not bsk_steps:
        print("No basket Gmsh CSVs found in", BASKET_DIR)
        return

    print("Building HTML ...")
    html = build_html(tri_steps, bsk_steps)
    os.makedirs(os.path.dirname(OUTPUT_HTML), exist_ok=True)
    with open(OUTPUT_HTML, "w", encoding="utf-8") as f:
        f.write(html)
    print("Written:", OUTPUT_HTML)


if __name__ == "__main__":
    main()