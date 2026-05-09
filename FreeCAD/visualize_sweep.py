"""
visualize_sweep.py  –  Interactive HTML report for P2 basket FEM tension sweep

Usage (standard Python, no FreeCAD needed):
    python visualize_sweep.py

Reads:
    FreeCAD/results/sweep_T????Nm.csv   (one per tension step)
    FreeCAD/results/sweep_summary.csv

Writes:
    FreeCAD/results/sweep_report.html
"""

import os
import glob
import csv
import math

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR  = os.path.join(SCRIPT_DIR, "results")
OUTPUT_HTML  = os.path.join(RESULTS_DIR, "sweep_report.html")

THRESHOLD_UM = 50.0          # planarity requirement (µm)
THRESHOLD_MM = THRESHOLD_UM / 1000.0

# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_step_csv(path):
    """Return dict of lists: node_id, x_mm, y_mm, z_mm, Ux, Uy, Uz, vonMises."""
    data = {k: [] for k in ["node_id", "x_mm", "y_mm", "z_mm",
                              "Ux_mm", "Uy_mm", "Uz_mm", "vonMises_MPa"]}
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            for k in data:
                data[k].append(float(row[k]) if k != "node_id" else int(row[k]))
    return data


def load_summary_csv(path):
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append({k: float(v) for k, v in row.items()})
    return rows


def discover_steps(results_dir):
    """Return sorted list of (T_Nm, filepath) for all sweep CSVs found."""
    pattern = os.path.join(results_dir, "sweep_T????Nm.csv")
    paths = sorted(glob.glob(pattern))
    steps = []
    for p in paths:
        basename = os.path.basename(p)          # sweep_T0300Nm.csv
        try:
            t_nm = int(basename[7:11])
            steps.append((t_nm, p))
        except ValueError:
            continue
    return steps

# ---------------------------------------------------------------------------
# Grid interpolation (scatter → regular grid via nearest-node binning)
# ---------------------------------------------------------------------------

def scatter_to_grid(x_list, y_list, z_list, nx=60, ny=60):
    """
    Bin scattered (x, y, z) values onto a regular grid.
    Returns (x_grid_1d, y_grid_1d, z_grid_2d[ny][nx]).
    Empty bins get None (shown as gaps in Plotly surface).
    """
    x_arr = x_list
    y_arr = y_list
    z_arr = z_list

    x_min, x_max = min(x_arr), max(x_arr)
    y_min, y_max = min(y_arr), max(y_arr)

    dx = (x_max - x_min) / nx
    dy = (y_max - y_min) / ny
    if dx == 0:
        dx = 1
    if dy == 0:
        dy = 1

    grid_sum   = [[0.0] * nx for _ in range(ny)]
    grid_count = [[0]   * nx for _ in range(ny)]

    for xv, yv, zv in zip(x_arr, y_arr, z_arr):
        ix = min(int((xv - x_min) / dx), nx - 1)
        iy = min(int((yv - y_min) / dy), ny - 1)
        grid_sum[iy][ix]   += zv
        grid_count[iy][ix] += 1

    z_grid = []
    for iy in range(ny):
        row = []
        for ix in range(nx):
            if grid_count[iy][ix] > 0:
                row.append(grid_sum[iy][ix] / grid_count[iy][ix])
            else:
                row.append(None)
        z_grid.append(row)

    x_1d = [x_min + (ix + 0.5) * dx for ix in range(nx)]
    y_1d = [y_min + (iy + 0.5) * dy for iy in range(ny)]

    return x_1d, y_1d, z_grid

# ---------------------------------------------------------------------------
# JSON helpers (no external libs)
# ---------------------------------------------------------------------------

def _jv(v):
    """Encode a Python value to a JSON-safe string fragment."""
    if v is None:
        return "null"
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float)):
        if math.isnan(v) or math.isinf(v):
            return "null"
        return repr(v)
    if isinstance(v, str):
        return '"' + v.replace("\\", "\\\\").replace('"', '\\"').replace("\n", "\\n") + '"'
    if isinstance(v, (list, tuple)):
        return "[" + ",".join(_jv(i) for i in v) + "]"
    if isinstance(v, dict):
        return "{" + ",".join('"' + k + '":' + _jv(vv) for k, vv in v.items()) + "}"
    return '"' + str(v) + '"'

# ---------------------------------------------------------------------------
# Plot builders  (return Plotly trace dicts as Python objects → serialise later)
# ---------------------------------------------------------------------------

def build_surface_traces(steps_data):
    """One Surface trace per tension step; visibility toggled by slider."""
    traces = []
    for t_nm, data in steps_data:
        t_ncm = t_nm / 100.0
        uz_um = [v * 1000 for v in data["Uz_mm"]]   # mm → µm
        x1d, y1d, z2d = scatter_to_grid(data["x_mm"], data["y_mm"], uz_um)
        trace = {
            "type": "surface",
            "x": x1d,
            "y": y1d,
            "z": z2d,
            "colorscale": "RdBu",
            "colorbar": {"title": "Uz (µm)", "x": 1.02},
            "name": f"{t_ncm:.0f} N/cm",
            "visible": False,
            "showscale": True,
            "hovertemplate": "x: %{x:.1f} mm<br>y: %{y:.1f} mm<br>Uz: %{z:.2f} µm<extra></extra>",
        }
        traces.append(trace)
    if traces:
        traces[0]["visible"] = True
    return traces


def build_vonmises_traces(steps_data):
    """One Surface trace per step for von Mises stress."""
    traces = []
    for t_nm, data in steps_data:
        t_ncm = t_nm / 100.0
        x1d, y1d, z2d = scatter_to_grid(data["x_mm"], data["y_mm"], data["vonMises_MPa"])
        trace = {
            "type": "surface",
            "x": x1d,
            "y": y1d,
            "z": z2d,
            "colorscale": "Hot",
            "colorbar": {"title": "σ_vm (MPa)", "x": 1.02},
            "name": f"{t_ncm:.0f} N/cm",
            "visible": False,
            "showscale": True,
            "hovertemplate": "x: %{x:.1f} mm<br>y: %{y:.1f} mm<br>σ: %{z:.3f} MPa<extra></extra>",
        }
        traces.append(trace)
    if traces:
        traces[0]["visible"] = True
    return traces


def build_heatmap_trace(summary_rows, steps_data):
    """Heatmap: Y axis = tension, X axis = x-position, value = max |Uz| per column."""
    if not steps_data:
        return None

    # Use first step to get x-grid
    sample_x = steps_data[0][1]["x_mm"]
    x_min, x_max = min(sample_x), max(sample_x)
    nx = 80
    dx = (x_max - x_min) / nx
    x_centers = [x_min + (i + 0.5) * dx for i in range(nx)]

    y_labels = []
    z_rows   = []

    for t_nm, data in steps_data:
        y_labels.append(t_nm / 100.0)
        uz_um = [abs(v) * 1000 for v in data["Uz_mm"]]

        col_max   = [0.0] * nx
        for xv, uv in zip(data["x_mm"], uz_um):
            ix = min(int((xv - x_min) / dx), nx - 1)
            if uv > col_max[ix]:
                col_max[ix] = uv
        z_rows.append(col_max)

    return {
        "type": "heatmap",
        "x": x_centers,
        "y": y_labels,
        "z": z_rows,
        "colorscale": "Viridis",
        "colorbar": {"title": "max|Uz| (µm)"},
        "hovertemplate": "x: %{x:.1f} mm<br>Tension: %{y:.0f} N/cm<br>max|Uz|: %{z:.2f} µm<extra></extra>",
    }


def build_summary_line_traces(summary_rows):
    """max|Uz| and max vonMises vs tension line charts."""
    t_ncm = [r["T_Ncm"]          for r in summary_rows]
    max_uz = [r["max_abs_Uz_mm"] * 1000 for r in summary_rows]  # → µm
    max_vm = [r["max_vonMises_MPa"]     for r in summary_rows]

    uz_trace = {
        "type": "scatter",
        "mode": "lines+markers",
        "x": t_ncm,
        "y": max_uz,
        "name": "max |Uz| (µm)",
        "line": {"color": "#1f77b4", "width": 2},
        "marker": {"size": 7},
        "hovertemplate": "T: %{x:.0f} N/cm<br>max|Uz|: %{y:.2f} µm<extra></extra>",
    }
    vm_trace = {
        "type": "scatter",
        "mode": "lines+markers",
        "x": t_ncm,
        "y": max_vm,
        "name": "max σ_vm (MPa)",
        "yaxis": "y2",
        "line": {"color": "#d62728", "width": 2, "dash": "dash"},
        "marker": {"size": 7},
        "hovertemplate": "T: %{x:.0f} N/cm<br>max σ_vm: %{y:.3f} MPa<extra></extra>",
    }
    threshold_line = {
        "type": "scatter",
        "mode": "lines",
        "x": [min(t_ncm), max(t_ncm)],
        "y": [THRESHOLD_UM, THRESHOLD_UM],
        "name": f"{THRESHOLD_UM:.0f} µm threshold",
        "line": {"color": "orange", "width": 2, "dash": "dot"},
        "hoverinfo": "skip",
    }
    return uz_trace, vm_trace, threshold_line

# ---------------------------------------------------------------------------
# Slider builder
# ---------------------------------------------------------------------------

def build_slider(n_steps, labels):
    steps = []
    for i, label in enumerate(labels):
        vis = [False] * n_steps
        vis[i] = True
        steps.append({
            "label": label,
            "method": "update",
            "args": [{"visible": vis}],
        })
    return {
        "active": 0,
        "currentvalue": {"prefix": "Tension: ", "suffix": " N/cm"},
        "pad": {"t": 50},
        "steps": steps,
    }

# ---------------------------------------------------------------------------
# HTML assembly
# ---------------------------------------------------------------------------

PLOTLY_CDN = "https://cdn.plot.ly/plotly-2.35.2.min.js"

def render_summary_table(summary_rows):
    rows_html = ""
    for r in summary_rows:
        uz_um = r["max_abs_Uz_mm"] * 1000
        ok    = uz_um <= THRESHOLD_UM
        style = "" if ok else ' style="background:#ffe0e0"'
        flag  = "" if ok else " ⚠"
        rows_html += (
            f"<tr{style}>"
            f"<td>{r['T_Ncm']:.0f}</td>"
            f"<td>{r['T_Nm']:.0f}</td>"
            f"<td>{r['delta_T_K']:.1f}</td>"
            f"<td>{uz_um:.2f}{flag}</td>"
            f"<td>{r['max_vonMises_MPa']:.4f}</td>"
            "</tr>\n"
        )
    return f"""
<table>
<thead>
<tr>
  <th>Tension (N/cm)</th><th>Tension (N/m)</th><th>ΔT (K)</th>
  <th>max |Uz| (µm)</th><th>max σ_vm (MPa)</th>
</tr>
</thead>
<tbody>
{rows_html}
</tbody>
</table>
<p style="font-size:0.85em;color:#888">⚠ exceeds {THRESHOLD_UM:.0f} µm planarity threshold</p>
"""


def build_html(steps_data, summary_rows):
    n = len(steps_data)
    labels = [f"{t / 100:.0f}" for t, _ in steps_data]

    surf_traces = build_surface_traces(steps_data)
    vm_traces   = build_vonmises_traces(steps_data)
    hm_trace    = build_heatmap_trace(summary_rows, steps_data)
    uz_tr, vm_tr, thr_tr = build_summary_line_traces(summary_rows)

    slider_surf = build_slider(n, labels)
    slider_vm   = build_slider(n, labels)

    summary_table = render_summary_table(summary_rows)

    # Serialise traces
    surf_json  = "[" + ",".join(_jv(t) for t in surf_traces)  + "]"
    vmsurf_json= "[" + ",".join(_jv(t) for t in vm_traces)    + "]"
    hm_json    = "[" + _jv(hm_trace) + "]" if hm_trace else "[]"
    line_json  = "[" + ",".join(_jv(t) for t in [uz_tr, vm_tr, thr_tr]) + "]"
    sl_surf    = _jv(slider_surf)
    sl_vm      = _jv(slider_vm)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>P2 Basket FEM Tension Sweep</title>
<script src="{PLOTLY_CDN}"></script>
<style>
  body {{ font-family: Arial, sans-serif; margin: 20px; background: #f9f9f9; }}
  h1   {{ color: #2c3e50; }}
  h2   {{ color: #34495e; margin-top: 40px; }}
  .plotbox {{ background: white; border-radius: 6px; box-shadow: 0 1px 4px #ccc;
              margin-bottom: 30px; padding: 10px; }}
  table {{ border-collapse: collapse; width: 100%; max-width: 700px; }}
  th, td {{ border: 1px solid #ddd; padding: 8px 12px; text-align: center; }}
  th {{ background: #2c3e50; color: white; }}
  tr:nth-child(even) {{ background: #f2f2f2; }}
</style>
</head>
<body>
<h1>P2 Basket Bulk Micromegas — FEM Tension Sweep</h1>
<p>Tension range: 3–10 N/cm &nbsp;|&nbsp; Model: v4 (CalculiX via FreeCAD FEM)
   &nbsp;|&nbsp; Planarity threshold: {THRESHOLD_UM:.0f} µm</p>

<h2>Out-of-plane displacement U<sub>z</sub> (µm)</h2>
<div id="surf_uz" class="plotbox" style="height:520px"></div>

<h2>Von Mises stress σ<sub>vm</sub> (MPa)</h2>
<div id="surf_vm" class="plotbox" style="height:520px"></div>

<h2>max |U<sub>z</sub>| heatmap — tension vs. x-position</h2>
<div id="heatmap" class="plotbox" style="height:400px"></div>

<h2>Sweep summary</h2>
<div id="lineplot" class="plotbox" style="height:400px"></div>
{summary_table}

<script>
// --- Surface: Uz ---
Plotly.newPlot("surf_uz", {surf_json}, {{
  scene: {{
    xaxis: {{title: "x (mm)"}},
    yaxis: {{title: "y (mm)"}},
    zaxis: {{title: "Uz (µm)"}},
    camera: {{eye: {{x:1.4, y:-1.4, z:0.8}}}}
  }},
  sliders: [{sl_surf}],
  margin: {{t:40, b:10}},
  paper_bgcolor: "white"
}}, {{responsive: true}});

// --- Surface: von Mises ---
Plotly.newPlot("surf_vm", {vmsurf_json}, {{
  scene: {{
    xaxis: {{title: "x (mm)"}},
    yaxis: {{title: "y (mm)"}},
    zaxis: {{title: "σ_vm (MPa)"}},
    camera: {{eye: {{x:1.4, y:-1.4, z:0.8}}}}
  }},
  sliders: [{sl_vm}],
  margin: {{t:40, b:10}},
  paper_bgcolor: "white"
}}, {{responsive: true}});

// --- Heatmap ---
Plotly.newPlot("heatmap", {hm_json}, {{
  xaxis: {{title: "x-position (mm)"}},
  yaxis: {{title: "Tension (N/cm)", ticksuffix: " N/cm"}},
  margin: {{t:20, b:60, l:80, r:20}},
  paper_bgcolor: "white",
  plot_bgcolor: "white"
}}, {{responsive: true}});

// --- Line summary ---
Plotly.newPlot("lineplot", {line_json}, {{
  xaxis: {{title: "Tension (N/cm)"}},
  yaxis: {{title: "max |Uz| (µm)", color: "#1f77b4"}},
  yaxis2: {{title: "max σ_vm (MPa)", overlaying: "y", side: "right", color: "#d62728"}},
  legend: {{x: 0.05, y: 0.95}},
  margin: {{t:20, b:60, l:70, r:70}},
  paper_bgcolor: "white",
  plot_bgcolor: "white"
}}, {{responsive: true}});
</script>
</body>
</html>
"""
    return html

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    steps = discover_steps(RESULTS_DIR)
    if not steps:
        print("No sweep CSVs found in", RESULTS_DIR)
        print("Run freecad_fem_setup.py first.")
        return

    print(f"Found {len(steps)} tension step(s): {[t for t, _ in steps]}")

    steps_data = []
    for t_nm, path in steps:
        print(f"  Loading T={t_nm} N/m ...")
        steps_data.append((t_nm, load_step_csv(path)))

    summary_path = os.path.join(RESULTS_DIR, "sweep_summary.csv")
    if os.path.exists(summary_path):
        summary_rows = load_summary_csv(summary_path)
        print(f"  Summary: {len(summary_rows)} row(s)")
    else:
        # Build summary on the fly from the per-step data
        print("  sweep_summary.csv not found — computing from step CSVs")
        summary_rows = []
        for t_nm, data in steps_data:
            uzs = data["Uz_mm"]
            summary_rows.append({
                "T_Nm":            float(t_nm),
                "T_Ncm":           t_nm / 100.0,
                "delta_T_K":       t_nm / 25.0,
                "max_Uz_mm":       max(uzs),
                "min_Uz_mm":       min(uzs),
                "max_abs_Uz_mm":   max(abs(v) for v in uzs),
                "max_vonMises_MPa": max(data["vonMises_MPa"]),
            })

    html = build_html(steps_data, summary_rows)
    os.makedirs(RESULTS_DIR, exist_ok=True)
    with open(OUTPUT_HTML, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"\nReport saved: {OUTPUT_HTML}")
    print("Open it in any browser.")


main()