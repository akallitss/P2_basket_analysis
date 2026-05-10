"""
buckling_geometry_check.py
Generates results/buckling_geometry_check.html

Interactive BC-zone map to verify the boundary conditions that will be used
in the buckling eigenvalue sweep (freecad_buckling_sweep.py).

No FreeCAD needed — reads node positions from an existing sweep CSV.
Run:
    py -3 FreeCAD/buckling_geometry_check.py
"""

import os
import csv
import glob
import math
import io
import base64

# ---------------------------------------------------------------------------
SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR  = os.path.join(SCRIPT_DIR, "results")
OUTPUT_HTML  = os.path.join(RESULTS_DIR, "buckling_geometry_check.html")

WALL_MM         = 1.0
PILLAR_PITCH_MM = 3.0
PILLAR_SNAP_MM  = 1.5
Z_TOL           = 0.05   # surface node threshold

# Sweep parameters for the upcoming simulation
TENSIONS_Nm     = [300, 500, 700, 1000]          # N/m
K_GLUE_LABELS   = ["Rigid (∞)", "1×10⁷", "1×10⁵", "1×10³", "Free (0)"]
K_GLUE_VALUES   = [1e10,        1e7,      1e5,      1e3,     0.0]   # N/m per node
N_EIGENMODES    = 5

PLOTLY_CDN = "https://cdn.plot.ly/plotly-2.35.2.min.js"

# ---------------------------------------------------------------------------

def load_surface_nodes():
    # Search both RESULTS_DIR and any immediate subdirectories
    paths = sorted(glob.glob(os.path.join(RESULTS_DIR, "sweep_T????Nm.csv")))
    if not paths:
        paths = sorted(glob.glob(os.path.join(RESULTS_DIR, "*", "sweep_T????Nm.csv")))
    if not paths:
        raise FileNotFoundError("No sweep CSVs found in " + RESULTS_DIR)
    xs, ys = [], []
    with open(paths[0], newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if abs(float(row["z_mm"])) < Z_TOL:
                xs.append(float(row["x_mm"]))
                ys.append(float(row["y_mm"]))
    print(f"  Surface nodes loaded: {len(xs)}  (from {os.path.basename(paths[0])})")
    return xs, ys


def make_pillar_grid(xmin, xmax, ymin, ymax):
    px_list, py_list = [], []
    px = xmin + WALL_MM + PILLAR_PITCH_MM / 2
    while px <= xmax - WALL_MM:
        py = ymin + WALL_MM + PILLAR_PITCH_MM / 2
        while py <= ymax - WALL_MM:
            px_list.append(px)
            py_list.append(py)
            py += PILLAR_PITCH_MM
        px += PILLAR_PITCH_MM
    return px_list, py_list


def classify(xs, ys, pillar_xs, pillar_ys, xmin, xmax, ymin, ymax):
    zones = []
    for x, y in zip(xs, ys):
        if (x < xmin + WALL_MM or x > xmax - WALL_MM or
                y < ymin + WALL_MM or y > ymax - WALL_MM):
            zones.append("wall")
        elif any(math.hypot(x - px, y - py) <= PILLAR_SNAP_MM
                 for px, py in zip(pillar_xs, pillar_ys)):
            zones.append("pillar")
        else:
            zones.append("free")
    return zones


def _jv(v):
    if v is None:               return "null"
    if isinstance(v, bool):     return "true" if v else "false"
    if isinstance(v, (int, float)):
        if math.isnan(v) or math.isinf(v): return "null"
        return repr(v)
    if isinstance(v, str):
        return '"' + v.replace("\\","\\\\").replace('"','\\"').replace("\n","\\n") + '"'
    if isinstance(v, (list, tuple)):
        return "[" + ",".join(_jv(i) for i in v) + "]"
    if isinstance(v, dict):
        return "{" + ",".join(f"{_jv(k)}:{_jv(val)}" for k, val in v.items()) + "}"
    return '"' + str(v) + '"'


def make_bc_schematic():
    """Matplotlib schematic: fixed wall BC vs. spring wall BC."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.lines import Line2D

    fig, axes = plt.subplots(1, 2, figsize=(7, 3))

    for ax, title, is_spring in zip(axes,
                                    ["Current study\n(static, fixed wall)",
                                     "Buckling sweep\n(spring wall — glue)"],
                                    [False, True]):
        ax.set_xlim(0, 6); ax.set_ylim(0, 5); ax.axis("off")
        ax.set_title(title, fontsize=10, pad=6)

        # Frame (hatched box on left)
        frame = mpatches.FancyBboxPatch((0.2, 0.5), 0.6, 4.0,
                                        boxstyle="square,pad=0",
                                        fc="#7f8c8d", ec="black", lw=1.5)
        ax.add_patch(frame)
        ax.text(0.5, 4.8, "Frame\n(SS)", ha="center", va="bottom", fontsize=8)

        # Mesh line
        ax.plot([0.8, 5.5], [2.5, 2.5], lw=3, color="#2980b9",
                solid_capstyle="round")
        ax.text(3.0, 2.8, "Mesh", ha="center", fontsize=9, color="#2980b9")

        # Pillar symbols
        for px in [2.0, 3.5, 5.0]:
            ax.plot([px, px], [0.5, 2.5], lw=2, color="#e74c3c")
            ax.plot(px, 0.5, "v", ms=8, color="#e74c3c")
        ax.text(3.5, 0.1, "Pillars  (Uz=0)", ha="center", fontsize=8,
                color="#e74c3c")

        if is_spring:
            # Spring symbol between frame and mesh end
            import numpy as np
            t = np.linspace(0, 1, 40)
            sx = 0.8 + t * 0.0          # constant x
            spring_x = [0.4] * 6        # left fixed
            # Draw a zigzag spring
            n_coils = 5
            seg_x = [0.4]
            seg_y = [2.5]
            for i in range(n_coils * 2):
                seg_x.append(0.4 + (0.2 if i % 2 == 0 else -0.2))
                seg_y.append(2.5 + (i + 1) * 4.0 / (n_coils * 2 + 1) - 2.0)
            seg_x.append(0.8); seg_y.append(2.5)
            ax.plot(seg_x, seg_y, lw=1.5, color="#f39c12")
            ax.text(0.5, 2.5, "k\nglue", ha="center", va="center",
                    fontsize=7, color="#f39c12",
                    bbox=dict(fc="white", ec="none", pad=1))
            # Fixed BC triangle at left
            ax.plot([0.2, 0.2], [1.5, 3.5], lw=0, marker="|",
                    ms=10, color="#2c3e50")
            ax.text(0.05, 2.5, "fixed\n(Uz)", ha="center", va="center",
                    fontsize=7, rotation=90, color="#2c3e50")
        else:
            # Fixed BC (triangle hatching)
            ax.plot([0.8, 0.8], [1.5, 3.5], lw=3, color="#2c3e50")
            for yh in [1.5, 2.0, 2.5, 3.0, 3.5]:
                ax.plot([0.8, 0.4], [yh, yh - 0.3], lw=1, color="#2c3e50")
            ax.text(0.4, 2.5, "fixed\nUx,Uy,Uz", ha="center", va="center",
                    fontsize=7, color="#2c3e50")

    fig.suptitle("Boundary condition comparison: static study vs. buckling sweep",
                 fontsize=11, y=1.02)
    fig.tight_layout()

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=140, bbox_inches="tight",
                facecolor="white")
    buf.seek(0)
    plt.close(fig)
    return base64.b64encode(buf.read()).decode("ascii")


def build_traces(xs, ys, zones, pillar_xs, pillar_ys,
                 xmin, xmax, ymin, ymax):
    xf = [x for x, z in zip(xs, zones) if z == "free"]
    yf = [y for y, z in zip(ys, zones) if z == "free"]
    xw = [x for x, z in zip(xs, zones) if z == "wall"]
    yw = [y for y, z in zip(ys, zones) if z == "wall"]
    xp = [x for x, z in zip(xs, zones) if z == "pillar"]
    yp = [y for y, z in zip(ys, zones) if z == "pillar"]

    n_wall   = len(xw)
    n_pillar = len(xp)
    n_free   = len(xf)

    traces = [
        {"type":"scatter","mode":"markers","x":xf,"y":yf,
         "name":f"Free interior  ({n_free:,} nodes)",
         "marker":{"color":"#bdc3c7","size":1.5,"opacity":0.6}},
        {"type":"scatter","mode":"markers","x":xp,"y":yp,
         "name":f"Pillar BC  Uz=0  ({n_pillar:,} nodes)",
         "marker":{"color":"#e74c3c","size":2.5,"opacity":0.8}},
        {"type":"scatter","mode":"markers","x":xw,"y":yw,
         "name":f"Wall spring BC  Ux,Uy=spring  Uz=0  ({n_wall:,} nodes)",
         "marker":{"color":"#2c3e50","size":2.5}},
        {"type":"scatter","mode":"markers",
         "x":pillar_xs,"y":pillar_ys,
         "name":f"Pillar centres  ({len(pillar_xs):,} pillars, {PILLAR_PITCH_MM} mm pitch)",
         "marker":{"symbol":"cross","color":"#c0392b","size":6,"line":{"width":1}}},
    ]
    return traces, n_wall, n_pillar, n_free


def build_html(xs, ys, zones, pillar_xs, pillar_ys, bounds, b64_schematic):
    xmin, xmax, ymin, ymax = bounds
    traces, n_wall, n_pillar, n_free = build_traces(
        xs, ys, zones, pillar_xs, pillar_ys, xmin, xmax, ymin, ymax)

    traces_json = "[" + ",".join(_jv(t) for t in traces) + "]"

    layout = {
        "title": "Buckling sweep — boundary condition zone map",
        "xaxis": {"title": "x (mm)", "scaleanchor": "y", "scaleratio": 1,
                  "showgrid": True, "gridcolor": "#eee"},
        "yaxis": {"title": "y (mm)", "showgrid": True, "gridcolor": "#eee"},
        "legend": {"x": 0.01, "y": 0.99, "bgcolor": "rgba(255,255,255,0.85)",
                   "bordercolor": "#ccc", "borderwidth": 1},
        "paper_bgcolor": "white", "plot_bgcolor": "white",
        "margin": {"t": 60, "b": 60, "l": 70, "r": 20},
    }

    # Sweep parameter table rows
    tension_rows = "".join(
        f"<tr><td>{t} N/m</td><td>{t/100:.0f} N/cm</td>"
        f"<td>{t/(1e9*0.5e-3*50e-6):.1f} K</td></tr>"
        for t in TENSIONS_Nm)
    kglue_rows = "".join(
        f"<tr><td>{lab}</td><td>{kv:.0e} N/m/node</td></tr>"
        for lab, kv in zip(K_GLUE_LABELS, K_GLUE_VALUES))

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Buckling geometry check</title>
<script src="{PLOTLY_CDN}"></script>
<style>
body{{font-family:Arial,sans-serif;margin:24px auto;max-width:1100px;
     background:#f4f6fa;color:#1a1a2e;}}
h1{{color:#1a252f;border-bottom:3px solid #e74c3c;padding-bottom:6px;}}
h2{{color:#2c3e50;margin-top:36px;border-left:5px solid #e74c3c;padding-left:10px;}}
.row{{display:flex;gap:24px;flex-wrap:wrap;}}
.card{{background:white;border-radius:8px;box-shadow:0 2px 8px #b0bec5;
       padding:16px;flex:1;min-width:280px;}}
table{{border-collapse:collapse;width:100%;margin-top:8px;}}
th,td{{border:1px solid #ccc;padding:7px 12px;font-size:0.9em;text-align:center;}}
th{{background:#2c3e50;color:white;}}
tr:nth-child(even){{background:#f0f4f8;}}
.plotbox{{background:white;border-radius:8px;box-shadow:0 2px 8px #b0bec5;
          margin:24px 0;padding:12px;}}
.note{{background:#fff8e1;border-left:5px solid #f39c12;
       padding:10px 16px;border-radius:4px;margin:12px 0;font-size:0.92em;}}
.ok {{background:#e8f5e9;border-left:5px solid #27ae60;
      padding:10px 16px;border-radius:4px;margin:12px 0;font-size:0.92em;}}
img{{max-width:100%;border-radius:6px;}}
</style>
</head>
<body>
<h1>Buckling eigenvalue sweep — geometry & BC check</h1>
<p>Verify that the boundary condition zones are correctly identified before
running <code>freecad_buckling_sweep.py</code>.
All node positions are read from the existing tension-sweep CSVs — no FreeCAD needed.</p>

<h2>1 &nbsp; BC zone map (interactive)</h2>
<div class="plotbox" id="main_plot" style="height:560px"></div>

<div class="row">
  <div class="card">
    <strong>Zone statistics</strong>
    <table>
      <tr><th>Zone</th><th>Colour</th><th>BC type</th><th>Node count</th></tr>
      <tr><td>Wall spring</td><td style="background:#2c3e50;color:white">■</td>
          <td>Ux,Uy: spring k_glue<br>Uz: fixed 0</td><td>{n_wall:,}</td></tr>
      <tr><td>Pillar snap</td><td style="background:#e74c3c;color:white">■</td>
          <td>Uz: fixed 0 (Ux,Uy free)</td><td>{n_pillar:,}</td></tr>
      <tr><td>Free interior</td><td style="background:#bdc3c7">■</td>
          <td>unconstrained</td><td>{n_free:,}</td></tr>
      <tr><td><strong>Pillar centres</strong></td><td>✕ red</td>
          <td>{PILLAR_PITCH_MM} mm pitch, snap r={PILLAR_SNAP_MM} mm</td>
          <td>{len(pillar_xs):,}</td></tr>
    </table>
    <p style="font-size:0.85em;color:#555;margin-top:8px">
    Active area: {xmax-xmin:.0f} × {ymax-ymin:.0f} mm &nbsp;|&nbsp;
    Wall strip: {WALL_MM} mm</p>
  </div>

  <div class="card">
    <strong>Physics</strong>
    <div class="ok">
    SS frame + SS mesh → same α → <strong>no differential thermal expansion</strong>
    during baking.  The wave mechanism is glue softening: the adhesive loses
    stiffness at temperature, partially releasing the pre-tension near the border.
    </div>
    <div class="note">
    The wall BC spring stiffness <em>k</em><sub>glue</sub> represents the in-plane
    stiffness of the softened adhesive per wall node.  Rigid (k→∞) reproduces the
    current fixed-wall model.  Free (k=0) models fully released glue.
    </div>
  </div>
</div>

<h2>2 &nbsp; BC schematic</h2>
<div style="text-align:center;margin:16px 0">
  <img src="data:image/png;base64,{b64_schematic}"
       alt="BC schematic" style="max-width:680px">
</div>

<h2>3 &nbsp; Planned sweep parameters</h2>
<div class="row">
  <div class="card">
    <strong>Pre-tension steps</strong>
    <table>
      <tr><th>T (N/m)</th><th>T (N/cm)</th><th>ΔT cooling (K)</th></tr>
      {tension_rows}
    </table>
  </div>
  <div class="card">
    <strong>Glue stiffness steps</strong>
    <table>
      <tr><th>Label</th><th>k_glue</th></tr>
      {kglue_rows}
    </table>
    <p style="font-size:0.85em;color:#555;margin-top:8px">
    Total combinations: {len(TENSIONS_Nm) * len(K_GLUE_VALUES)} runs &nbsp;|&nbsp;
    Eigenmodes per run: {N_EIGENMODES}</p>
  </div>
  <div class="card">
    <strong>Pillar BCs</strong>
    <table>
      <tr><th>Parameter</th><th>Value</th></tr>
      <tr><td>Pitch</td><td>{PILLAR_PITCH_MM} mm × {PILLAR_PITCH_MM} mm</td></tr>
      <tr><td>Snap radius</td><td>{PILLAR_SNAP_MM} mm</td></tr>
      <tr><td>DOF constrained</td><td>Uz = 0 only</td></tr>
      <tr><td>In-plane</td><td>free (Ux, Uy)</td></tr>
    </table>
  </div>
</div>

<script>
Plotly.newPlot("main_plot", {traces_json}, {_jv(layout)}, {{responsive:true}});
</script>
</body>
</html>"""
    return html, n_wall, n_pillar, n_free


def main():
    import sys
    print("Loading surface nodes ...")
    xs, ys = load_surface_nodes()

    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    bounds = (xmin, xmax, ymin, ymax)

    print("Computing BC zones ...")
    pillar_xs, pillar_ys = make_pillar_grid(xmin, xmax, ymin, ymax)
    zones = classify(xs, ys, pillar_xs, pillar_ys, xmin, xmax, ymin, ymax)

    print("Generating BC schematic ...")
    b64_schematic = make_bc_schematic()

    print("Building HTML ...")
    html, n_wall, n_pillar, n_free = build_html(
        xs, ys, zones, pillar_xs, pillar_ys, bounds, b64_schematic)

    os.makedirs(RESULTS_DIR, exist_ok=True)
    with open(OUTPUT_HTML, "w", encoding="utf-8") as f:
        f.write(html)

    print(f"""
Done.  {OUTPUT_HTML}

Zone summary:
  Wall spring BC : {n_wall:>6,} nodes  (Ux,Uy spring  +  Uz fixed)
  Pillar BC      : {n_pillar:>6,} nodes  (Uz fixed only)
  Free interior  : {n_free:>6,} nodes
  Pillar centres : {len(pillar_xs):>6,}

Sweep: {len(TENSIONS_Nm)} tensions × {len(K_GLUE_VALUES)} k_glue values = {len(TENSIONS_Nm)*len(K_GLUE_VALUES)} runs
If the zones look correct, run:
    py -3 FreeCAD/freecad_buckling_sweep.py
""")


if __name__ == "__main__":
    main()