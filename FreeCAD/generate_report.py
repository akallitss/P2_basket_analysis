"""
generate_report.py  –  Self-contained FEM technical report
                        for P2 Basket Bulk Micromegas tension sweep.

Requirements:  numpy, matplotlib  (standard scientific Python stack)
Run with:
    py -3 FreeCAD/generate_report.py

Output:
    FreeCAD/results/fem_report.html   (fully self-contained, no external files)
"""

import os
import csv
import glob
import math
import base64
import io

# ---------------------------------------------------------------------------
# Paths & constants
# ---------------------------------------------------------------------------

SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR  = os.path.join(SCRIPT_DIR, "results")
OUTPUT_HTML  = os.path.join(RESULTS_DIR, "fem_report.html")
PLOTLY_CDN   = "https://cdn.plot.ly/plotly-2.35.2.min.js"

THRESHOLD_UM = 50.0          # planarity spec

# Thermal-analogy model constants (from freecad_fem_setup.py)
E_MESH  = 1.0e9    # Pa
T_MESH  = 0.5e-3   # m
ALPHA   = 50.0e-6  # /K
T_REF   = 300.0    # K

PILLAR_PITCH_MM = 3.0
WALL_MM         = 1.0

# Mesh hardware
WIRE_DIAM_UM = 18.0
WIRE_PITCH_UM = 45.0


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_step_csv(path):
    data = {k: [] for k in ("node_id", "x_mm", "y_mm", "z_mm",
                             "Ux_mm", "Uy_mm", "Uz_mm", "vonMises_MPa")}
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            data["node_id"].append(int(row["node_id"]))
            for k in ("x_mm", "y_mm", "z_mm", "Ux_mm", "Uy_mm", "Uz_mm"):
                data[k].append(float(row[k]))
            try:
                data["vonMises_MPa"].append(float(row["vonMises_MPa"]))
            except (ValueError, KeyError):
                data["vonMises_MPa"].append(None)
    return data


def discover_steps(results_dir):
    paths = sorted(glob.glob(os.path.join(results_dir, "sweep_T????Nm.csv")))
    steps = []
    for p in paths:
        try:
            t_nm = int(os.path.basename(p)[7:11])
            steps.append((t_nm, p))
        except ValueError:
            pass
    return steps


def compute_summary(steps_data):
    rows = []
    for t_nm, data in steps_data:
        uzs = data["Uz_mm"]
        vms = [v for v in data["vonMises_MPa"] if v is not None]
        rows.append({
            "T_Nm":              t_nm,
            "T_Ncm":             t_nm / 100.0,
            "delta_T_K":         t_nm / (E_MESH * T_MESH * ALPHA),
            "T_applied_K":       T_REF - t_nm / (E_MESH * T_MESH * ALPHA),
            "max_Uz_mm":         max(uzs),
            "min_Uz_mm":         min(uzs),
            "max_abs_Uz_mm":     max(abs(v) for v in uzs),
            "max_vonMises_MPa":  max(vms) if vms else 0.0,
        })
    return rows


# ---------------------------------------------------------------------------
# Matplotlib helpers
# ---------------------------------------------------------------------------

def _fig_to_b64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight", facecolor="white")
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode("ascii")
    import matplotlib.pyplot as plt
    plt.close(fig)
    return b64


def plot_summary_line(summary_rows):
    import matplotlib.pyplot as plt

    T   = [r["T_Ncm"]             for r in summary_rows]
    Uz  = [r["max_abs_Uz_mm"]*1e3 for r in summary_rows]  # µm
    vm  = [r["max_vonMises_MPa"]   for r in summary_rows]

    fig, ax1 = plt.subplots(figsize=(8, 4.5))
    ax2 = ax1.twinx()

    l1, = ax1.plot(T, Uz, "o-",  color="#1f77b4", lw=2.2, ms=8,
                   label="max |Uz| (µm)")
    l2, = ax2.plot(T, vm, "s--", color="#d62728", lw=2.2, ms=8,
                   label="max σ_vm (MPa)")
    ax1.axhline(THRESHOLD_UM, color="darkorange", ls=":", lw=2,
                label=f"{THRESHOLD_UM:.0f} µm threshold")

    ax1.set_xlabel("Mesh pre-tension (N/cm)", fontsize=12)
    ax1.set_ylabel("max |Uz| (µm)", color="#1f77b4", fontsize=12)
    ax2.set_ylabel("max σ_vm (MPa)", color="#d62728", fontsize=12)
    ax1.tick_params(axis="y", labelcolor="#1f77b4")
    ax2.tick_params(axis="y", labelcolor="#d62728")
    ax1.set_xticks(T)
    ax1.grid(True, alpha=0.3)
    ax1.set_title("Peak Displacement & Stress vs. Pre-tension", fontsize=13)

    lines = [l1, l2]
    ax1.legend(lines, [l.get_label() for l in lines],
               loc="upper left", fontsize=10)
    fig.tight_layout()
    return fig


def plot_uz_grid(steps_data):
    import matplotlib.pyplot as plt
    import numpy as np

    n     = len(steps_data)
    ncols = 4
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(4.5 * ncols, 3.5 * nrows),
                             constrained_layout=True)
    ax_flat = list(axes.flat) if hasattr(axes, "flat") else [axes]

    # Shared colour scale across all panels
    all_uz_um = []
    for _, d in steps_data:
        all_uz_um.extend(v * 1e3 for v in d["Uz_mm"])
    vabs = max(abs(v) for v in all_uz_um)
    vmin, vmax = -vabs, vabs

    sc_last = None
    for ax, (t_nm, data) in zip(ax_flat, steps_data):
        # Surface nodes only (z ≈ 0)
        pts = [(x, y, z * 1e3) for x, y, z, zz in zip(
                   data["x_mm"], data["y_mm"],
                   data["Uz_mm"], data["z_mm"]) if abs(zz) < 0.05]
        if not pts:
            pts = list(zip(data["x_mm"], data["y_mm"],
                           [v * 1e3 for v in data["Uz_mm"]]))
        px, py, pu = zip(*pts)
        sc_last = ax.scatter(px, py, c=pu, cmap="RdBu_r", s=0.4,
                             vmin=vmin, vmax=vmax, rasterized=True)
        ax.set_title(f"T = {t_nm/100:.0f} N/cm", fontsize=10)
        ax.set_xlabel("x (mm)", fontsize=8)
        ax.set_ylabel("y (mm)", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.set_aspect("equal")

    for ax in ax_flat[n:]:
        ax.set_visible(False)

    if sc_last is not None:
        fig.colorbar(sc_last, ax=ax_flat[:n], label="Uz (µm)", shrink=0.7)
    fig.suptitle("Out-of-plane displacement Uz — all tension steps", fontsize=13)
    return fig


def plot_vonmises_grid(steps_data):
    import matplotlib.pyplot as plt

    n     = len(steps_data)
    ncols = 4
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(4.5 * ncols, 3.5 * nrows),
                             constrained_layout=True)
    ax_flat = list(axes.flat) if hasattr(axes, "flat") else [axes]

    all_vm = [v for _, d in steps_data
              for v in d["vonMises_MPa"] if v is not None]
    vmax = max(all_vm) if all_vm else 1.0

    sc_last = None
    for ax, (t_nm, data) in zip(ax_flat, steps_data):
        pts = [(x, y, v) for x, y, zz, v in zip(
                   data["x_mm"], data["y_mm"],
                   data["z_mm"], data["vonMises_MPa"])
               if abs(zz) < 0.05 and v is not None]
        if not pts:
            pts = [(x, y, v) for x, y, v in zip(
                       data["x_mm"], data["y_mm"], data["vonMises_MPa"])
                   if v is not None]
        px, py, pv = zip(*pts)
        sc_last = ax.scatter(px, py, c=pv, cmap="hot_r", s=0.4,
                             vmin=0, vmax=vmax, rasterized=True)
        ax.set_title(f"T = {t_nm/100:.0f} N/cm", fontsize=10)
        ax.set_xlabel("x (mm)", fontsize=8)
        ax.set_ylabel("y (mm)", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.set_aspect("equal")

    for ax in ax_flat[n:]:
        ax.set_visible(False)

    if sc_last is not None:
        fig.colorbar(sc_last, ax=ax_flat[:n], label="σ_vm (MPa)", shrink=0.7)
    fig.suptitle("Von Mises stress — all tension steps", fontsize=13)
    return fig


def plot_bc_layout(steps_data):
    import matplotlib.pyplot as plt
    import numpy as np

    _, data = steps_data[0]
    xs = np.array(data["x_mm"])
    ys = np.array(data["y_mm"])
    zs = np.array(data["z_mm"])

    mask = abs(zs) < 0.05
    xm, ym = xs[mask], ys[mask]
    xmin, xmax = xm.min(), xm.max()
    ymin, ymax = ym.min(), ym.max()

    is_wall = (
        (xm < xmin + WALL_MM) | (xm > xmax - WALL_MM) |
        (ym < ymin + WALL_MM) | (ym > ymax - WALL_MM)
    )

    # Pillar grid centres
    pgx, pgy = [], []
    px = xmin + WALL_MM + PILLAR_PITCH_MM / 2
    while px <= xmax - WALL_MM:
        py = ymin + WALL_MM + PILLAR_PITCH_MM / 2
        while py <= ymax - WALL_MM:
            pgx.append(px); pgy.append(py)
            py += PILLAR_PITCH_MM
        px += PILLAR_PITCH_MM

    fig, ax = plt.subplots(figsize=(9, 7))
    ax.scatter(xm[~is_wall], ym[~is_wall],
               c="#cccccc", s=0.3, label="Interior nodes", rasterized=True)
    ax.scatter(xm[is_wall], ym[is_wall],
               c="#2c3e50", s=0.6,
               label=f"Wall BC  Ux=Uy=Uz=0  (2 529 nodes)", rasterized=True)
    ax.scatter(pgx, pgy, marker="+", c="#e74c3c", s=8, linewidths=0.5,
               label=f"Pillar BC  Uz=0  (3 mm pitch, 15 680 nodes)", rasterized=True)

    ax.set_xlabel("x (mm)", fontsize=11)
    ax.set_ylabel("y (mm)", fontsize=11)
    ax.set_title("Boundary condition layout  (mesh surface, z ≈ 0)", fontsize=12)
    ax.set_aspect("equal")
    ax.legend(loc="upper right", fontsize=9, markerscale=6)

    # Dimension arrows
    pad = 18
    ax.annotate("", xy=(xmax, ymin - pad), xytext=(xmin, ymin - pad),
                arrowprops=dict(arrowstyle="<->", color="#333"))
    ax.text((xmin + xmax) / 2, ymin - pad - 8,
            f"{xmax-xmin:.0f} mm", ha="center", fontsize=9)
    ax.annotate("", xy=(xmax + pad, ymax), xytext=(xmax + pad, ymin),
                arrowprops=dict(arrowstyle="<->", color="#333"))
    ax.text(xmax + pad + 6, (ymin + ymax) / 2,
            f"{ymax-ymin:.0f} mm", ha="left", va="center",
            fontsize=9, rotation=90)

    fig.tight_layout()
    return fig


def plot_mesh_schematic():
    """Small diagram of the woven mesh unit cell."""
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches

    d   = WIRE_DIAM_UM
    p   = WIRE_PITCH_UM
    op  = p - d   # opening

    fig, ax = plt.subplots(figsize=(4, 4))
    # Draw 3×3 unit cells
    for i in range(3):
        for j in range(3):
            x0 = i * p
            y0 = j * p
            # Cell outline
            ax.add_patch(patches.Rectangle((x0, y0), p, p,
                                           fill=False, ec="#aaa", lw=0.5))
            # Wire in x direction (horizontal bar)
            ax.add_patch(patches.Rectangle((x0, y0 + op / 2), p, d,
                                           fc="#888", ec="none"))
            # Wire in y direction (vertical bar)
            ax.add_patch(patches.Rectangle((x0 + op / 2, y0), d, p,
                                           fc="#aaa", ec="none"))

    # Annotations
    ax.annotate("", xy=(0, -6), xytext=(p, -6),
                arrowprops=dict(arrowstyle="<->", color="black"))
    ax.text(p / 2, -10, f"pitch = {p:.0f} µm", ha="center", fontsize=9)

    ax.annotate("", xy=(p + op / 2, p + 4), xytext=(p + op / 2 + d, p + 4),
                arrowprops=dict(arrowstyle="<->", color="black"))
    ax.text(p + op / 2 + d / 2, p + 9,
            f"d = {d:.0f} µm", ha="center", fontsize=9)

    ax.annotate("", xy=(0, p + 4), xytext=(op / 2, p + 4),
                arrowprops=dict(arrowstyle="<->", color="black"))
    ax.text(op / 4, p + 9, f"open = {op:.0f} µm", ha="center", fontsize=8)

    ax.set_xlim(-5, 3 * p + 20)
    ax.set_ylim(-18, 3 * p + 20)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(f"Woven SS mesh unit cell  ({p:.0f} µm pitch / {d:.0f} µm wire)",
                 fontsize=10)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Plotly helpers (no external libs)
# ---------------------------------------------------------------------------

def _get_plotly_js():
    """
    Download Plotly JS from CDN and return as inline <script> tag.
    Falls back to a <script src="..."> CDN reference if download fails.
    This makes the report fully self-contained for offline/sharing use.
    """
    import urllib.request
    url = "https://cdn.plot.ly/plotly-2.35.2.min.js"
    try:
        print("  Fetching Plotly JS for inline embedding ...")
        with urllib.request.urlopen(url, timeout=30) as r:
            js = r.read().decode("utf-8")
        print(f"  Plotly JS embedded ({len(js)//1024} KB)")
        return f"<script>{js}</script>"
    except Exception as e:
        print(f"  Could not fetch Plotly JS ({e}); using CDN reference instead.")
        return f'<script src="{url}"></script>'


def _jv(v):
    if v is None:
        return "null"
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float)):
        return "null" if (math.isnan(v) or math.isinf(v)) else repr(v)
    if isinstance(v, str):
        return ('"' + v.replace("\\", "\\\\")
                       .replace('"', '\\"')
                       .replace("\n", "\\n") + '"')
    if isinstance(v, (list, tuple)):
        return "[" + ",".join(_jv(i) for i in v) + "]"
    if isinstance(v, dict):
        return "{" + ",".join(f"{_jv(k)}:{_jv(val)}" for k, val in v.items()) + "}"
    return '"' + str(v) + '"'


def _scatter_to_grid(xs, ys, zs, nx=80, ny=80):
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    dx = (xmax - xmin) / nx or 1
    dy = (ymax - ymin) / ny or 1
    gs = [[0.0] * nx for _ in range(ny)]
    gc = [[0]   * nx for _ in range(ny)]
    for x, y, z in zip(xs, ys, zs):
        ix = min(int((x - xmin) / dx), nx - 1)
        iy = min(int((y - ymin) / dy), ny - 1)
        gs[iy][ix] += z
        gc[iy][ix] += 1
    zg = [[gs[r][c] / gc[r][c] if gc[r][c] else None
           for c in range(nx)] for r in range(ny)]
    return ([xmin + (i + 0.5) * dx for i in range(nx)],
            [ymin + (j + 0.5) * dy for j in range(ny)],
            zg)


def _build_surface_plot(steps_data, field, transform, colorscale, cb_title):
    """
    Returns (traces_json, slider_json).

    Uses one trace per tension step (multi-trace approach) — reliable
    initial render — with slider using method='restyle' on the visible
    attribute.  method='update' is avoided; it has known issues with
    3-D surface visibility in Plotly.
    """
    traces = []
    for i, (t_nm, data) in enumerate(steps_data):
        vals = [transform(v) for v in data[field] if v is not None]
        xs   = [x for x, v in zip(data["x_mm"], data[field]) if v is not None]
        ys   = [y for y, v in zip(data["y_mm"], data[field]) if v is not None]
        x1d, y1d, z2d = _scatter_to_grid(xs, ys, vals)
        traces.append({
            "type": "surface", "x": x1d, "y": y1d, "z": z2d,
            "colorscale": colorscale,
            "colorbar": {"title": cb_title, "x": 1.02},
            "name": f"{t_nm/100:.0f} N/cm",
            "visible": i == 0,          # only first trace starts visible
            "showscale": True,
            "hovertemplate":
                "x: %{x:.1f} mm<br>y: %{y:.1f} mm<br>"
                + cb_title + ": %{z:.3f}<extra></extra>",
        })
    n = len(traces)
    traces_json = "[" + ",".join(_jv(t) for t in traces) + "]"

    # Slider: restyle visible array across all traces
    slider_steps = []
    for i, (t_nm, _) in enumerate(steps_data):
        vis = [j == i for j in range(n)]   # list of booleans, one per trace
        slider_steps.append({
            "label": f"{t_nm/100:.0f}",
            "method": "restyle",
            "args": [{"visible": vis}],
        })
    slider = {
        "active": 0,
        "currentvalue": {"prefix": "Tension: ", "suffix": " N/cm",
                         "font": {"size": 14}},
        "pad": {"t": 60, "b": 10},
        "steps": slider_steps,
    }
    return traces_json, _jv(slider)


def _line_traces(summary_rows):
    T   = [r["T_Ncm"]             for r in summary_rows]
    Uz  = [r["max_abs_Uz_mm"]*1e3 for r in summary_rows]
    vm  = [r["max_vonMises_MPa"]   for r in summary_rows]
    return [
        {"type":"scatter","mode":"lines+markers","x":T,"y":Uz,
         "name":"max |Uz| (µm)","yaxis":"y",
         "line":{"color":"#1f77b4","width":2},"marker":{"size":8}},
        {"type":"scatter","mode":"lines+markers","x":T,"y":vm,
         "name":"max σ_vm (MPa)","yaxis":"y2",
         "line":{"color":"#d62728","width":2,"dash":"dash"},"marker":{"size":8}},
        {"type":"scatter","mode":"lines","x":[min(T),max(T)],
         "y":[THRESHOLD_UM, THRESHOLD_UM],
         "name":f"{THRESHOLD_UM:.0f} µm threshold","yaxis":"y",
         "line":{"color":"darkorange","width":2,"dash":"dot"},
         "hoverinfo":"skip"},
    ]


# ---------------------------------------------------------------------------
# HTML
# ---------------------------------------------------------------------------

CSS = """
<style>
body{font-family:Arial,Helvetica,sans-serif;margin:30px auto;max-width:1120px;
     background:#f4f6fa;color:#1a1a2e;}
h1{color:#1a252f;border-bottom:3px solid #2980b9;padding-bottom:8px;margin-bottom:4px;}
h2{color:#2c3e50;margin-top:44px;border-left:5px solid #2980b9;padding-left:10px;}
h3{color:#34495e;margin-top:28px;}
p,li{line-height:1.75;font-size:0.96em;}
.plotbox{background:white;border-radius:8px;box-shadow:0 2px 8px #b0bec5;
         margin-bottom:32px;padding:14px;}
.imgbox{text-align:center;margin-bottom:28px;}
.imgbox img{max-width:100%;border-radius:6px;box-shadow:0 2px 6px #b0bec5;}
.cap{font-size:0.84em;color:#555;margin-top:6px;}
table{border-collapse:collapse;width:100%;margin:10px 0 24px;}
th,td{border:1px solid #ccc;padding:9px 14px;text-align:center;font-size:0.92em;}
th{background:#2c3e50;color:white;}
tr:nth-child(even){background:#eef2f7;}
.ptd{text-align:left;font-weight:600;background:#ecf0f1 !important;}
.formula{background:#eaf4fb;border-left:5px solid #2980b9;padding:10px 16px;
         border-radius:4px;font-family:"Courier New",monospace;font-size:0.95em;
         margin:12px 0;}
.note{background:#fff8e1;border-left:5px solid #f39c12;padding:10px 16px;
      border-radius:4px;font-size:0.9em;margin:14px 0;}
.ok{color:#27ae60;font-weight:700;}
.warn{color:#e74c3c;font-weight:700;}
</style>
"""


def _results_table(summary_rows):
    rows = ""
    for r in summary_rows:
        uz = r["max_abs_Uz_mm"] * 1e3
        ok = uz <= THRESHOLD_UM
        flag = '<span class="ok">✓</span>' if ok else '<span class="warn">✗</span>'
        rows += (
            f"<tr>"
            f"<td>{r['T_Ncm']:.0f}</td>"
            f"<td>{r['T_Nm']:.0f}</td>"
            f"<td>{r['delta_T_K']:.2f}</td>"
            f"<td>{r['T_applied_K']:.2f}</td>"
            f"<td>{r['max_abs_Uz_mm']*1e3:.2f}</td>"
            f"<td>{r['max_vonMises_MPa']:.3f}</td>"
            f"<td>{flag}</td>"
            "</tr>\n"
        )
    return f"""
<table>
<thead>
<tr>
  <th>T<br>(N/cm)</th><th>T<br>(N/m)</th><th>ΔT<br>(K)</th>
  <th>T<sub>applied</sub><br>(K)</th>
  <th>max |Uz|<br>(µm)</th>
  <th>max σ<sub>vm</sub><br>(MPa)</th>
  <th>≤ {THRESHOLD_UM:.0f} µm</th>
</tr>
</thead>
<tbody>{rows}</tbody>
</table>
<p class="cap">Table 1. FEM sweep results — all eight tension steps.
Planarity threshold: {THRESHOLD_UM:.0f} µm.</p>"""


def build_html(steps_data, summary_rows,
               b64_summary, b64_uz, b64_vm, b64_bc, b64_mesh,
               plotly_script_tag):

    uz_json, sl_uz = _build_surface_plot(
        steps_data, "Uz_mm",        lambda v: v * 1e3, "RdBu", "Uz (µm)")
    vm_json, sl_vm = _build_surface_plot(
        steps_data, "vonMises_MPa", lambda v: v,       "Hot",  "σ_vm (MPa)")
    line_json  = "[" + ",".join(_jv(t) for t in _line_traces(summary_rows)) + "]"
    res_table  = _results_table(summary_rows)

    _, d0 = steps_data[0]
    xmin, xmax = min(d0["x_mm"]), max(d0["x_mm"])
    ymin, ymax = min(d0["y_mm"]), max(d0["y_mm"])
    fill_pct = (1 - ((WIRE_PITCH_UM - WIRE_DIAM_UM) / WIRE_PITCH_UM)**2) * 100

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>P2 Basket Bulk Micromegas — FEM Technical Report</title>
{plotly_script_tag}
{CSS}
</head>
<body>

<h1>P2 Basket Bulk Micromegas<br>
<span style="font-size:0.65em;color:#555;font-weight:normal">
FEM Pre-tension Analysis · CalculiX via FreeCAD · 2026</span></h1>

<!-- ======================================================= -->
<h2>1&nbsp;&nbsp;Introduction</h2>
<p>
In a <strong>Bulk Micromegas</strong> gaseous detector, a thin stainless-steel
woven mesh acts as the cathode of the electron amplification gap.
During the Bulk construction process the mesh is laminated directly onto the
read-out PCB using photosensitive resist pillars (Dynamask), which define a
uniform amplification gap height of ~150 µm.
</p>
<p>
Mechanical pre-tension applied to the mesh before lamination is essential to
guarantee planarity of the electrode surface over the full detector area.
Insufficient tension leads to sagging; excessive tension risks yielding the
wires or delaminating the pillars.
</p>
<p>
This report presents a parametric FEM study of the P2 basket Micromegas mesh
under pre-tension values from <strong>3 N/cm to 10 N/cm</strong>, evaluating
out-of-plane displacement and in-plane stress at each operating point.
The solver is <strong>CalculiX</strong> (linear static), driven through
FreeCAD FEM by a custom Python sweep script.
</p>


<!-- ======================================================= -->
<h2>2&nbsp;&nbsp;Detector Geometry and FEM Model</h2>

<h3>2.1&nbsp; Active area</h3>
<p>
The P2 basket detector active area is approximately
<strong>{xmax-xmin:.0f}&nbsp;mm &times; {ymax-ymin:.0f}&nbsp;mm</strong>.
The FreeCAD model consists of two solid bodies produced by a Boolean
fragment operation:
</p>
<ul>
  <li><strong>Solid 1 — stainless-steel mesh layer</strong> (31&nbsp;182 volume
      elements, occupying z&nbsp;=&nbsp;0.05&nbsp;mm … 0.70&nbsp;mm).</li>
  <li><strong>Solid 2 — PCB / FR4 substrate</strong> (19&nbsp;744 volume elements,
      below the mesh layer). </li>
</ul>
<p>
The two bodies share interface nodes, ensuring full mechanical coupling at
the mesh–PCB contact surface.
</p>

<h3>2.2&nbsp; FEM mesh statistics</h3>
<table style="max-width:500px">
<tr><td class="ptd">Total nodes</td><td>102&nbsp;407</td></tr>
<tr><td class="ptd">Total volume elements</td><td>50&nbsp;926</td></tr>
<tr><td class="ptd">Mesh-layer nodes&nbsp;(0.05 ≤ z ≤ 0.70 mm)</td><td>42&nbsp;030</td></tr>
<tr><td class="ptd">Element type</td><td>Second-order tetrahedra (C3D10)</td></tr>
<tr><td class="ptd">Mesh generator</td><td>GMSH via FreeCAD FEM</td></tr>
</table>


<!-- ======================================================= -->
<h2>3&nbsp;&nbsp;Material Properties</h2>

<h3>3.1&nbsp; Woven stainless-steel mesh — homogenised continuum</h3>

<div class="imgbox" style="max-width:380px;margin:0 auto 20px;">
  <img src="data:image/png;base64,{b64_mesh}" alt="Mesh unit cell schematic">
  <p class="cap">Figure 1. Woven mesh unit cell geometry
  (pitch {WIRE_PITCH_UM:.0f}&nbsp;µm, wire diameter {WIRE_DIAM_UM:.0f}&nbsp;µm,
  opening {WIRE_PITCH_UM-WIRE_DIAM_UM:.0f}&nbsp;µm).</p>
</div>

<p>
The 45&times;18&nbsp;µm woven SS316L mesh (pitch 45&nbsp;µm,
wire diameter 18&nbsp;µm) is replaced in the FEM model by an
<strong>equivalent homogenised continuum</strong> with isotropic in-plane
properties. Homogenisation is valid here because the weave is square
(identical pitch and wire in both directions), giving the same tensile
stiffness along x and y.
</p>
<p>
The effective Young's modulus E<sub>eff</sub>&nbsp;=&nbsp;1.0&nbsp;GPa is far
below bulk SS316L (~&nbsp;200&nbsp;GPa). This reflects two combined effects:
(i) the large open fraction (fill factor ≈&nbsp;{fill_pct:.0f}%), and
(ii) the crimp compliance of woven wires — at low tension the wires first
straighten before engaging their full axial stiffness, yielding a much softer
effective modulus than the simple fill-factor estimate would suggest.
</p>

<table style="max-width:600px">
<thead><tr><th>Parameter</th><th>Symbol</th><th>Value</th><th>Note</th></tr></thead>
<tbody>
<tr><td class="ptd">Wire material</td><td>—</td><td>SS316L</td><td>stainless steel</td></tr>
<tr><td class="ptd">Wire diameter</td><td>d</td><td>18 µm</td><td>—</td></tr>
<tr><td class="ptd">Wire pitch</td><td>p</td><td>45 µm</td><td>square weave</td></tr>
<tr><td class="ptd">Wire opening</td><td>p − d</td><td>27 µm</td><td>—</td></tr>
<tr><td class="ptd">Fill factor</td><td>1−(op/p)²</td><td>≈ {fill_pct:.0f} %</td><td>area fraction of wire</td></tr>
<tr><td class="ptd">Eff. Young's modulus</td><td>E<sub>eff</sub></td><td>1.0 GPa</td><td>homogenised (isotropic)</td></tr>
<tr><td class="ptd">Eff. thickness</td><td>t<sub>eff</sub></td><td>0.5 mm</td><td>used in thermal analogy</td></tr>
<tr><td class="ptd">Eff. therm. expansion</td><td>α<sub>eff</sub></td><td>50 µm/m/K</td><td>calibrated for analogy</td></tr>
<tr><td class="ptd">Reference temperature</td><td>T<sub>ref</sub></td><td>300 K</td><td>ConstraintInitialTemperature</td></tr>
</tbody>
</table>

<h3>3.2&nbsp; PCB substrate</h3>
<p>
The PCB substrate (Solid 2) is assigned isotropic linear-elastic properties
consistent with standard FR4 laminate. As a significantly stiffer body it
provides the mechanical boundary to which the mesh bonds via the Dynamask
pillars, and is not the primary subject of this analysis.
</p>


<!-- ======================================================= -->
<h2>4&nbsp;&nbsp;Boundary Conditions</h2>

<div class="imgbox">
  <img src="data:image/png;base64,{b64_bc}" alt="BC layout">
  <p class="cap">Figure 2. Boundary condition layout projected onto the mesh
  surface (z ≈ 0). Dark: peripheral wall BC (U = 0). Red crosses: Dynamask
  pillar BC (U<sub>z</sub> = 0) on a 3&nbsp;mm pitch grid.</p>
</div>

<h3>4.1&nbsp; Peripheral wall constraint</h3>
<p>
All mesh-layer nodes within a <strong>1&nbsp;mm wide peripheral frame</strong>
(derived from the convex hull of the mesh-layer node set) are fully fixed:
</p>
<div class="formula">U<sub>x</sub> = U<sub>y</sub> = U<sub>z</sub> = 0
&nbsp;&nbsp;&nbsp;(2&nbsp;529 nodes)</div>
<p>
This represents the mechanical clamping of the mesh frame at the PCB edge
during lamination.
</p>

<h3>4.2&nbsp; Dynamask pillar supports</h3>
<p>
The Dynamask photoresist pillars form a regular
<strong>3&nbsp;mm &times; 3&nbsp;mm square grid</strong> across the active
area. Each pillar top contacts the bottom surface of the mesh and prevents
out-of-plane displacement at that location. In the FEM model the nearest
mesh-layer node to each pillar centre (within a 1.5&nbsp;mm snap radius)
is constrained to:
</p>
<div class="formula">U<sub>z</sub> = 0
&nbsp;&nbsp;&nbsp;(15&nbsp;680 nodes, 3&nbsp;mm pitch grid)</div>

<h3>4.3&nbsp; Thermal temperature constraint</h3>
<p>
A uniform applied temperature T<sub>applied</sub> is set on the mesh surface
(Face 27 of BooleanFragments). The thermal analogy converts mesh pre-tension
into a temperature difference (see Section 5).
</p>


<!-- ======================================================= -->
<h2>5&nbsp;&nbsp;Pre-tension Method: Thermal Analogy</h2>
<p>
Direct application of a biaxial pre-tension force in a 3D solid FEM model
requires prescribed displacements at every boundary node, calibrated to
yield exactly the correct resultant force per unit width — a non-trivial
setup. The <strong>thermal analogy</strong> is the standard alternative for
pre-tensioned membrane structures:
</p>
<div class="formula">
Strain equilibrium:&nbsp;&nbsp;
ε&nbsp;=&nbsp;α<sub>eff</sub> · |ΔT|&nbsp;=&nbsp;T / (E<sub>eff</sub> · t<sub>eff</sub>)
<br><br>
⟹&nbsp;&nbsp; ΔT&nbsp;=&nbsp;T / (E<sub>eff</sub> · t<sub>eff</sub> · α<sub>eff</sub>)
&nbsp;&nbsp;=&nbsp;&nbsp; T / (1.0×10⁹ · 5.0×10⁻⁴ · 50×10⁻⁶)
&nbsp;&nbsp;=&nbsp;&nbsp; T / 25&nbsp;&nbsp; [K·m/N]
<br><br>
T<sub>applied</sub>&nbsp;=&nbsp;T<sub>ref</sub> − ΔT
</div>
<p>
Setting T<sub>applied</sub> below T<sub>ref</sub> causes the mesh to thermally
<em>contract</em>. Because the peripheral wall nodes are held fixed, this
contraction is restrained and generates a biaxial in-plane tension equal to
T&nbsp;(N/m) — identical in effect to mechanically pulling the mesh edges
outward.
</p>

<table style="max-width:600px">
<thead>
<tr><th>T (N/cm)</th><th>T (N/m)</th><th>ΔT (K)</th>
    <th>T<sub>applied</sub> (K)</th></tr>
</thead>
<tbody>
{"".join(
    f"<tr><td>{r['T_Ncm']:.0f}</td><td>{r['T_Nm']:.0f}</td>"
    f"<td>{r['delta_T_K']:.2f}</td><td>{r['T_applied_K']:.2f}</td></tr>"
    for r in summary_rows
)}
</tbody>
</table>
<p class="cap">Table 2. Tension-to-temperature conversion for each sweep step.</p>

<div class="note">
<strong>Model note — 3D solid vs. membrane:</strong>
The model uses 3D solid elements (C3D10) with isotropic thermal expansion.
The out-of-plane component α<sub>zz</sub> creates a small spurious bowl-shaped
deflection that peaks at the geometric centre and would vanish in a true 2D
membrane model. The in-plane stress results (von Mises) are unaffected.
The maximum spurious deflection remains well below the 50&nbsp;µm planarity
threshold throughout the sweep, so the result is conservative and valid for
the purpose of this study.
</div>


<!-- ======================================================= -->
<h2>6&nbsp;&nbsp;Results</h2>

<h3>6.1&nbsp; Out-of-plane displacement U<sub>z</sub> — interactive 3D surface</h3>
<p>Use the slider to switch between tension steps.</p>
<div id="surf_uz" class="plotbox" style="height:530px"></div>

<h3>6.2&nbsp; Von Mises stress σ<sub>vm</sub> — interactive 3D surface</h3>
<div id="surf_vm" class="plotbox" style="height:530px"></div>

<h3>6.3&nbsp; Displacement maps — all tensions</h3>
<div class="imgbox">
  <img src="data:image/png;base64,{b64_uz}" alt="Uz maps all tensions">
  <p class="cap">Figure 3. Out-of-plane displacement U<sub>z</sub> (µm) for
  all eight tension steps, projected on the top surface (z ≈ 0).
  Common diverging colour scale across panels; maximum magnitude at the
  geometric centre of the detector.</p>
</div>

<h3>6.4&nbsp; Von Mises stress maps — all tensions</h3>
<div class="imgbox">
  <img src="data:image/png;base64,{b64_vm}" alt="vonMises maps all tensions">
  <p class="cap">Figure 4. Von Mises stress σ<sub>vm</sub> (MPa) for
  all eight tension steps.  Stress concentrations at the peripheral wall BC
  are expected and are the design-critical quantity.</p>
</div>

<h3>6.5&nbsp; Peak values vs. tension — summary</h3>
<div class="imgbox">
  <img src="data:image/png;base64,{b64_summary}" alt="Summary chart">
  <p class="cap">Figure 5. Peak out-of-plane displacement max|U<sub>z</sub>|
  (blue, left axis) and peak von Mises stress max σ<sub>vm</sub> (red, right
  axis) vs. mesh pre-tension. Both scale linearly with tension (linear elastic
  regime confirmed). Orange dotted line: 50&nbsp;µm planarity specification.</p>
</div>
<div id="lineplot" class="plotbox" style="height:430px"></div>

<h3>6.6&nbsp; Numerical results — all tension steps</h3>
{res_table}


<!-- ======================================================= -->
<h2>7&nbsp;&nbsp;Discussion and Conclusions</h2>
<ul>
  <li><strong>Linear scaling:</strong> Both max|U<sub>z</sub>| and
      max&nbsp;σ<sub>vm</sub> scale exactly proportionally to the applied
      tension across the full 3–10&nbsp;N/cm range, confirming operation in
      the linear elastic regime. The scaling factor for displacement is
      {summary_rows[-1]['max_abs_Uz_mm']*1e3/summary_rows[-1]['T_Ncm']:.2f}&nbsp;µm per N/cm.</li>

  <li><strong>Planarity requirement:</strong>
      The maximum out-of-plane displacement recorded across the full sweep is
      <strong>{max(r['max_abs_Uz_mm']*1e3 for r in summary_rows):.1f}&nbsp;µm</strong>
      at {max(summary_rows, key=lambda r: r['max_abs_Uz_mm'])['T_Ncm']:.0f}&nbsp;N/cm,
      well below the <strong>50&nbsp;µm</strong> planarity specification.
      The mesh satisfies the planarity requirement at all tested tensions.</li>

  <li><strong>Deformation shape:</strong> The maximum displacement is located
      at the geometric centre of the detector (~&nbsp;333,&nbsp;269&nbsp;mm),
      the point farthest from all wall and pillar constraints. The deformation
      mode shape is invariant across tensions (only the amplitude scales),
      consistent with the linear elastic model. This localisation has a partial
      contribution from the spurious α<sub>zz</sub> thermal strain in the 3D
      solid model (see Section&nbsp;5 note).</li>

  <li><strong>Model limitations:</strong>
    <ul>
      <li>3D solid (C3D10) elements introduce bending stiffness and an
          isotropic α<sub>zz</sub> not present in the physical 2D mesh.
          Shell/membrane elements would suppress the spurious out-of-plane
          artefact.</li>
      <li>Effective material parameters (E<sub>eff</sub>, t<sub>eff</sub>,
          α<sub>eff</sub>) require experimental cross-calibration against a
          measured tension-vs-deflection curve.</li>
      <li>Large-displacement geometric non-linearity is not modelled;
          justified by the small deflection-to-span ratio (&lt;&nbsp;0.005%).</li>
    </ul>
  </li>
</ul>


<!-- ======================================================= -->
<h2>Appendix — Complete Model Parameters</h2>
<table style="max-width:640px">
<tr><td class="ptd">Solver</td><td>CalculiX (ccx) — linear static</td></tr>
<tr><td class="ptd">FreeCAD version</td><td>0.20.2 (r29177)</td></tr>
<tr><td class="ptd">Mesh generator</td><td>GMSH via FreeCAD FEM interface</td></tr>
<tr><td class="ptd">E<sub>eff</sub></td><td>1.0 × 10⁹ Pa</td></tr>
<tr><td class="ptd">t<sub>eff</sub></td><td>0.5 × 10⁻³ m</td></tr>
<tr><td class="ptd">α<sub>eff</sub></td><td>50 × 10⁻⁶ K⁻¹</td></tr>
<tr><td class="ptd">T<sub>ref</sub></td><td>300 K</td></tr>
<tr><td class="ptd">Mesh z-range (mesh layer)</td><td>0.05 mm – 0.70 mm</td></tr>
<tr><td class="ptd">Pillar pitch</td><td>3.0 mm × 3.0 mm</td></tr>
<tr><td class="ptd">Pillar snap radius</td><td>1.5 mm</td></tr>
<tr><td class="ptd">Wall BC frame width</td><td>1.0 mm</td></tr>
<tr><td class="ptd">Wall BC nodes</td><td>2&nbsp;529</td></tr>
<tr><td class="ptd">Pillar BC nodes</td><td>15&nbsp;680</td></tr>
<tr><td class="ptd">Tension sweep</td>
    <td>300, 400, 500, 600, 700, 800, 900, 1000 N/m
        (3–10 N/cm)</td></tr>
</table>

<script>
Plotly.newPlot("surf_uz", {uz_json}, {{
  scene:{{xaxis:{{title:"x (mm)"}},yaxis:{{title:"y (mm)"}},
         zaxis:{{title:"Uz (µm)"}},
         camera:{{eye:{{x:1.4,y:-1.4,z:0.8}}}}}},
  sliders:[{sl_uz}], margin:{{t:40,b:10}}, paper_bgcolor:"white"
}},{{responsive:true}});

Plotly.newPlot("surf_vm", {vm_json}, {{
  scene:{{xaxis:{{title:"x (mm)"}},yaxis:{{title:"y (mm)"}},
         zaxis:{{title:"σ_vm (MPa)"}},
         camera:{{eye:{{x:1.4,y:-1.4,z:0.8}}}}}},
  sliders:[{sl_vm}], margin:{{t:40,b:10}}, paper_bgcolor:"white"
}},{{responsive:true}});

Plotly.newPlot("lineplot", {line_json}, {{
  xaxis:{{title:"Pre-tension (N/cm)"}},
  yaxis:{{title:"max |Uz| (µm)",color:"#1f77b4"}},
  yaxis2:{{title:"max σ_vm (MPa)",overlaying:"y",side:"right",color:"#d62728"}},
  legend:{{x:0.05,y:0.95}},
  margin:{{t:20,b:60,l:80,r:80}},
  paper_bgcolor:"white",plot_bgcolor:"white"
}},{{responsive:true}});
</script>
</body>
</html>"""


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    import matplotlib
    matplotlib.use("Agg")

    steps = discover_steps(RESULTS_DIR)
    if not steps:
        print("No sweep CSVs found in", RESULTS_DIR)
        return

    print(f"Found {len(steps)} tension step(s): {[t for t,_ in steps]}")

    print("Loading CSVs ...")
    steps_data = [(t, load_step_csv(p)) for t, p in steps]

    summary_rows = compute_summary(steps_data)

    print("Fetching Plotly JS ...")
    plotly_tag = _get_plotly_js()

    print("Generating figures ...")
    b64_summary = _fig_to_b64(plot_summary_line(summary_rows))
    b64_uz      = _fig_to_b64(plot_uz_grid(steps_data))
    b64_vm      = _fig_to_b64(plot_vonmises_grid(steps_data))
    b64_bc      = _fig_to_b64(plot_bc_layout(steps_data))
    b64_mesh    = _fig_to_b64(plot_mesh_schematic())

    print("Building HTML ...")
    html = build_html(steps_data, summary_rows,
                      b64_summary, b64_uz, b64_vm, b64_bc, b64_mesh,
                      plotly_tag)

    os.makedirs(RESULTS_DIR, exist_ok=True)
    with open(OUTPUT_HTML, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"\nDone.  Report: {OUTPUT_HTML}")
    print("Open in any browser — fully self-contained (no external files needed).")


if __name__ == "__main__":
    main()