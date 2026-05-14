"""
visualize_baking_buckle.py  -  Interactive HTML report for baking buckling results
P2 Basket Bulk Micromegas mesh

Usage (standard Python, no FreeCAD needed):
    python visualize_baking_buckle.py

Reads:
    FreeCAD/results/baking_buckle_summary.csv
    FreeCAD/results/baking_buckle_T????Nm_mode??.csv

Writes:
    FreeCAD/results/baking_buckle_report.html
"""

import os
import glob
import csv
import math

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(SCRIPT_DIR, "results")
OUTPUT_HTML = os.path.join(RESULTS_DIR, "baking_buckle_report.html")
PLOTLY_CDN  = "https://cdn.plot.ly/plotly-2.35.2.min.js"

T_BAKE_C    = 146.0
DT_BAKE_REAL = 121.0


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def discover_mode_csvs(results_dir):
    pattern = os.path.join(results_dir, "baking_buckle_T????Nm_mode??.csv")
    results = []
    for p in sorted(glob.glob(pattern)):
        b = os.path.basename(p)
        try:
            t_nm = int(b[15:19])
            mode = int(b[26:28])
            results.append((t_nm, mode, p))
        except (ValueError, IndexError):
            continue
    return results


def load_mode_csv(path):
    data = {k: [] for k in ["node_id", "x_mm", "y_mm", "z_mm",
                              "Ux_norm", "Uy_norm", "Uz_norm",
                              "eigenvalue", "T_buckle_C"]}
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            data["node_id"].append(int(row["node_id"]))
            for k in ("x_mm", "y_mm", "z_mm", "Ux_norm", "Uy_norm", "Uz_norm",
                      "eigenvalue", "T_buckle_C"):
                data[k].append(float(row[k]))
    return data


def load_summary_csv(path):
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            ev = float(row["eigenvalue"])
            rows.append({
                "T_Nm":       float(row["T_Nm"]),
                "T_Ncm":      float(row["T_Ncm"]),
                "mode":       int(row["mode"]),
                "eigenvalue": ev,
                "T_buckle_C": float(row["T_buckle_C"]),
                "survived":   row["survived"].strip().lower() == "yes",
                "physical":   ev > 0,
            })
    return rows


def build_summary_from_modes(mode_records):
    rows = []
    for (t_nm, mode), data in mode_records.items():
        ev = data["eigenvalue"][0] if data["eigenvalue"] else float("nan")
        t_buckle = 25.0 + ev * DT_BAKE_REAL if not math.isnan(ev) else float("nan")
        rows.append({
            "T_Nm":       float(t_nm),
            "T_Ncm":      t_nm / 100.0,
            "mode":       mode,
            "eigenvalue": ev,
            "T_buckle_C": t_buckle,
            "survived":   ev >= 1.0,
            "physical":   ev > 0,
        })
    rows.sort(key=lambda r: (r["T_Nm"], r["mode"]))
    return rows


# ---------------------------------------------------------------------------
# Grid interpolation
# ---------------------------------------------------------------------------

def scatter_to_grid(x_list, y_list, z_list, nx=80, ny=80):
    x_min, x_max = min(x_list), max(x_list)
    y_min, y_max = min(y_list), max(y_list)
    dx = (x_max - x_min) / nx or 1.0
    dy = (y_max - y_min) / ny or 1.0

    grid_sum   = [[0.0] * nx for _ in range(ny)]
    grid_count = [[0]   * nx for _ in range(ny)]
    for xv, yv, zv in zip(x_list, y_list, z_list):
        ix = min(int((xv - x_min) / dx), nx - 1)
        iy = min(int((yv - y_min) / dy), ny - 1)
        grid_sum[iy][ix]   += zv
        grid_count[iy][ix] += 1

    z_grid = []
    for iy in range(ny):
        row = []
        for ix in range(nx):
            row.append(grid_sum[iy][ix] / grid_count[iy][ix]
                       if grid_count[iy][ix] > 0 else None)
        z_grid.append(row)

    x_1d = [x_min + (ix + 0.5) * dx for ix in range(nx)]
    y_1d = [y_min + (iy + 0.5) * dy for iy in range(ny)]
    return x_1d, y_1d, z_grid


# ---------------------------------------------------------------------------
# JSON helpers
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
        return '"' + v.replace("\\", "\\\\").replace('"', '\\"').replace("\n", "\\n") + '"'
    if isinstance(v, (list, tuple)):
        return "[" + ",".join(_jv(i) for i in v) + "]"
    if isinstance(v, dict):
        return "{" + ",".join('"' + k + '":' + _jv(vv) for k, vv in v.items()) + "}"
    return '"' + str(v) + '"'


# ---------------------------------------------------------------------------
# Chart builders
# ---------------------------------------------------------------------------

MODE_COLORS = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00"]
SURVIVED_COLOR = "#2ca02c"
BUCKLED_COLOR  = "#d62728"


def build_eigenvalue_bar_traces(summary_rows):
    tensions = sorted(set(r["T_Nm"] for r in summary_rows))
    modes    = sorted(set(r["mode"] for r in summary_rows))

    traces = []
    for mi, m in enumerate(modes):
        x, y, colors, hover = [], [], [], []
        for t in tensions:
            match = next((r for r in summary_rows if r["T_Nm"] == t and r["mode"] == m), None)
            if match:
                ev = match["eigenvalue"]
                is_phys = match.get("physical", ev > 0)
                x.append("{:.0f} N/m".format(int(t)))
                y.append(ev if is_phys else None)
                colors.append(SURVIVED_COLOR if match["survived"] else BUCKLED_COLOR)
                note = "(reverse-direction, not a baking mode)" if not is_phys else \
                       ("SURVIVED" if match["survived"] else "BUCKLES before full bake")
                hover.append(
                    "T={:.0f} N/m  Mode {}<br>λ={:.4f}<br>{}".format(t, m, ev, note)
                )
        traces.append({
            "type": "bar",
            "name": "Mode {}".format(m),
            "x": x,
            "y": y,
            "marker": {"color": MODE_COLORS[mi % len(MODE_COLORS)]},
            "hovertext": hover,
            "hovertemplate": "%{hovertext}<extra></extra>",
        })
    return traces


def build_tbuckle_scatter_traces(summary_rows):
    tensions = sorted(set(r["T_Nm"] for r in summary_rows))
    modes    = sorted(set(r["mode"] for r in summary_rows))

    traces = []
    for mi, m in enumerate(modes):
        x, y, hover = [], [], []
        for t in tensions:
            match = next((r for r in summary_rows if r["T_Nm"] == t and r["mode"] == m), None)
            if match and match.get("physical", match["eigenvalue"] > 0):
                x.append("{:.0f} N/m".format(int(t)))
                y.append(match["T_buckle_C"])
                hover.append("T={:.0f} N/m  Mode {}<br>T_buckle≈{:.1f} °C".format(
                    t, m, match["T_buckle_C"]))
        if x:
            traces.append({
                "type": "scatter",
                "name": "Mode {} T_buckle".format(m),
                "x": x,
                "y": y,
                "mode": "lines+markers",
                "yaxis": "y2",
                "line": {"color": MODE_COLORS[mi % len(MODE_COLORS)], "dash": "dot", "width": 1},
                "marker": {"size": 5},
                "hovertext": hover,
                "hovertemplate": "%{hovertext}<extra></extra>",
                "showlegend": False,
            })
    return traces


def build_mode_shape_traces(mode_records):
    traces = []
    labels = []
    tensions = sorted(set(t for t, _ in mode_records))
    for t_nm in tensions:
        modes = sorted(m for (tt, m) in mode_records if tt == t_nm)
        for mode in modes:
            data = mode_records[(t_nm, mode)]
            ev = data["eigenvalue"][0] if data["eigenvalue"] else float("nan")
            x1d, y1d, z2d = scatter_to_grid(data["x_mm"], data["y_mm"], data["Uz_norm"])
            is_phys = ev > 0
            survived = ev >= 1.0
            if not is_phys:
                status = "(reverse-direction, not a baking mode)"
            elif survived:
                status = "✓ SURVIVED"
            else:
                status = "✗ BUCKLES at ~{:.0f} °C".format(25.0 + ev * DT_BAKE_REAL)
            label = "T={} N/m  Mode {}  λ={:.4f}  {}".format(t_nm, mode, ev, status)
            traces.append({
                "type": "heatmap",
                "x": x1d,
                "y": y1d,
                "z": z2d,
                "colorscale": "RdBu",
                "zmid": 0,
                "colorbar": {"title": "Uz (norm)"},
                "visible": False,
                "hovertemplate": "x: %{x:.1f} mm<br>y: %{y:.1f} mm<br>Uz_norm: %{z:.4f}<extra></extra>",
                "name": label,
            })
            labels.append(label)

    if traces:
        traces[0]["visible"] = True
    return traces, labels


def build_combined_slider(n_traces, labels):
    steps = []
    for i, label in enumerate(labels):
        vis = [False] * n_traces
        vis[i] = True
        steps.append({
            "label": "M{}".format((i % 5) + 1),
            "method": "update",
            "args": [{"visible": vis}, {"title": label}],
        })
    return steps


# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------

def _first_positive_ev(summary_rows, t_nm):
    candidates = sorted(
        (r for r in summary_rows if r["T_Nm"] == t_nm and r.get("physical", r["eigenvalue"] > 0)),
        key=lambda r: r["mode"]
    )
    return candidates[0] if candidates else None


def render_summary_table(summary_rows):
    tensions = sorted(set(r["T_Nm"] for r in summary_rows))
    modes    = sorted(set(r["mode"] for r in summary_rows))

    header = "<tr><th>T (N/m)</th><th>T (N/cm)</th>"
    for m in modes:
        header += "<th>λ<sub>{}</sub></th>".format(m)
    header += "<th>λ<sub>crit</sub></th><th>T_buckle (°C)</th><th>Status</th></tr>"

    body = ""
    for t in tensions:
        ncm = t / 100.0
        row_evs = []
        for m in modes:
            match = next((r for r in summary_rows if r["T_Nm"] == t and r["mode"] == m), None)
            if match:
                ev = match["eigenvalue"]
                is_phys = match.get("physical", ev > 0)
                if is_phys:
                    row_evs.append("{:.4f}".format(ev))
                else:
                    row_evs.append(
                        '<span style="color:#aaa;font-style:italic" '
                        'title="Reverse-direction mode — not a baking buckling mode">'
                        '{:.3g}</span>'.format(ev)
                    )
            else:
                row_evs.append("—")

        crit = _first_positive_ev(summary_rows, t)
        if crit:
            survived   = crit["survived"]
            lambda_crit = "{:.4f}".format(crit["eigenvalue"])
            t_buckle    = "{:.1f}".format(crit["T_buckle_C"])
        else:
            survived    = False
            lambda_crit = "—"
            t_buckle    = "—"

        status_style = "" if survived else ' style="background:#ffe0e0"'
        status_text  = "SURVIVED" if survived else "BUCKLES ✗"
        status_color = 'style="color:#2ca02c;font-weight:bold"' if survived \
                       else 'style="color:#d62728;font-weight:bold"'

        body += "<tr{}>".format(status_style)
        body += "<td>{:.0f}</td><td>{:.1f}</td>".format(t, ncm)
        for ev_str in row_evs:
            body += "<td>{}</td>".format(ev_str)
        body += "<td><strong>{}</strong></td>".format(lambda_crit)
        body += "<td>{}</td>".format(t_buckle)
        body += '<td {}>{}</td>'.format(status_color, status_text)
        body += "</tr>\n"

    n_neg = sum(1 for r in summary_rows if not r.get("physical", r["eigenvalue"] > 0))
    neg_note = ""
    if n_neg:
        neg_note = (
            "<p style='font-size:0.82em;color:#888'>"
            "<em>Grey-italic values</em>: negative eigenvalues = reverse-direction modes. "
            "These indicate that this mode is stabilised by the reference load "
            "(baking compression) — buckling only occurs if the load is reversed. "
            "They are not baking failure modes. "
            "λ<sub>crit</sub> is the first positive eigenvalue per tension."
            "</p>"
        )

    return """
<table>
<thead>{}</thead>
<tbody>{}</tbody>
</table>
<p style="font-size:0.85em;color:#888">
  λ &lt; 1 → buckles before end of bake at ≈ (25 + λ×121) °C &nbsp;|&nbsp;
  λ = 1 → buckles exactly at 146 °C &nbsp;|&nbsp;
  λ &gt; 1 → survives full bake
</p>
{}
""".format(header, body, neg_note)


# ---------------------------------------------------------------------------
# HTML assembly
# ---------------------------------------------------------------------------

def build_html(summary_rows, mode_records):
    ev_bar_traces    = build_eigenvalue_bar_traces(summary_rows)
    tbuckle_traces   = build_tbuckle_scatter_traces(summary_rows)
    mode_traces, mode_labels = build_mode_shape_traces(mode_records)

    tensions = sorted(set(r["T_Nm"] for r in summary_rows))
    n_per_tension = len(set(r["mode"] for r in summary_rows))

    # Build tension-group slider steps (groups of n_per_tension traces)
    n_mode_traces = len(mode_traces)
    slider_steps  = build_combined_slider(n_mode_traces, mode_labels)

    # Build tension-level steps for the outer slider (jump by n_per_tension)
    tension_steps = []
    for ti, t in enumerate(tensions):
        vis = [False] * n_mode_traces
        base = ti * n_per_tension
        if base < n_mode_traces:
            vis[base] = True
        tension_steps.append({
            "label": "{:.0f} N/m".format(t),
            "method": "update",
            "args": [{"visible": vis}],
        })

    all_bar_traces = ev_bar_traces + tbuckle_traces
    bar_json       = "[" + ",".join(_jv(t) for t in all_bar_traces) + "]"
    mode_json      = "[" + ",".join(_jv(t) for t in mode_traces) + "]"
    tension_sl_j   = _jv(tension_steps)
    mode_sl_j      = _jv(slider_steps)
    summary_table  = render_summary_table(summary_rows)

    tensions_str = ", ".join("{:.0f}".format(t) for t in tensions)
    n_completed  = len(tensions)

    html = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>P2 Basket — Baking Buckling Analysis</title>
<script src="{cdn}"></script>
<style>
  body    {{ font-family: Arial, sans-serif; margin: 24px; background: #f7f7f7; }}
  h1      {{ color: #2c3e50; margin-bottom: 4px; }}
  .sub    {{ color: #666; margin-top: 0; font-size: 0.95em; }}
  h2      {{ color: #34495e; margin-top: 40px; border-bottom: 1px solid #ccc;
             padding-bottom: 4px; }}
  .box    {{ background: white; border-radius: 6px;
             box-shadow: 0 1px 4px #ccc; margin-bottom: 28px; padding: 12px; }}
  table   {{ border-collapse: collapse; width: 100%; max-width: 860px; }}
  th, td  {{ border: 1px solid #ddd; padding: 7px 12px; text-align: center; }}
  th      {{ background: #2c3e50; color: white; }}
  tr:nth-child(even) {{ background: #f4f4f4; }}
  .note   {{ font-size: 0.82em; color: #888; margin-top: 6px; }}
</style>
</head>
<body>

<h1>P2 Basket Bulk Micromegas — Baking Buckling Analysis</h1>
<p class="sub">
  Tensions: {tensions_str} N/m &nbsp;|&nbsp;
  Completed: {n_completed}/4 tensions &nbsp;|&nbsp;
  Bake: 25 °C → 146 °C &nbsp;|&nbsp;
  Δα = α<sub>SS</sub> − α<sub>FR4</sub> = 2×10<sup>-6</sup> /K &nbsp;|&nbsp;
  λ=1 ↔ buckles exactly at 146 °C
</p>

<h2>Eigenvalues vs. pre-tension</h2>
<div id="ev_bar" class="box" style="height:460px"></div>

<h2>Mode shape browser</h2>
<div id="mode_shape" class="box" style="height:520px"></div>

<h2>Summary table</h2>
<div class="box">
{summary_table}
</div>

<script>
// ---- Eigenvalue bar chart + T_buckle overlay ----
Plotly.newPlot("ev_bar", {bar_json}, {{
  barmode: "group",
  xaxis: {{ title: "Pre-tension" }},
  yaxis: {{ title: "Eigenvalue λ", zeroline: false }},
  yaxis2: {{
    title: "T_buckle (°C)",
    overlaying: "y",
    side: "right",
    showgrid: false,
  }},
  shapes: [{{
    type: "line",
    xref: "paper", x0: 0, x1: 1,
    yref: "y",     y0: 1, y1: 1,
    line: {{ color: "black", width: 2, dash: "dash" }},
  }}],
  annotations: [{{
    xref: "paper", x: 1.0,
    yref: "y",     y: 1.0,
    text: "λ = 1<br>(146 °C)",
    showarrow: false,
    xanchor: "right",
    font: {{ size: 11 }},
  }}],
  legend: {{ x: 0.01, y: 0.99, bgcolor: "rgba(255,255,255,0.8)" }},
  margin: {{ t: 20, b: 60, l: 60, r: 80 }},
  paper_bgcolor: "white",
  plot_bgcolor: "white",
}}, {{responsive: true}});

// ---- Mode shape heatmap with two-level slider ----
Plotly.newPlot("mode_shape", {mode_json}, {{
  xaxis: {{ title: "x (mm)" }},
  yaxis: {{ title: "y (mm)", scaleanchor: "x" }},
  sliders: [
    {{
      active: 0,
      currentvalue: {{ prefix: "Tension: " }},
      pad: {{ t: 10 }},
      steps: {tension_sl_j},
      y: 1.06,
    }},
    {{
      active: 0,
      currentvalue: {{ prefix: "Mode: " }},
      pad: {{ t: 55 }},
      steps: {mode_sl_j},
      y: 1.0,
    }}
  ],
  margin: {{ t: 110, b: 60, l: 60, r: 20 }},
  paper_bgcolor: "white",
  plot_bgcolor: "white",
}}, {{responsive: true}});
</script>

</body>
</html>
""".format(
        cdn=PLOTLY_CDN,
        tensions_str=tensions_str,
        n_completed=n_completed,
        bar_json=bar_json,
        mode_json=mode_json,
        tension_sl_j=tension_sl_j,
        mode_sl_j=mode_sl_j,
        summary_table=summary_table,
    )
    return html


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    mode_entries = discover_mode_csvs(RESULTS_DIR)
    if not mode_entries:
        print("No baking_buckle_T????Nm_mode??.csv files found in", RESULTS_DIR)
        print("Run freecad_baking_buckling.py first.")
        return

    tensions_found = sorted(set(t for t, _, _ in mode_entries))
    print("Found {} mode CSV(s) across tensions: {}".format(
        len(mode_entries), tensions_found))

    mode_records = {}
    for t_nm, mode, path in mode_entries:
        print("  Loading T={} N/m  Mode {} ...".format(t_nm, mode))
        mode_records[(t_nm, mode)] = load_mode_csv(path)

    summary_path = os.path.join(RESULTS_DIR, "baking_buckle_summary.csv")
    if os.path.isfile(summary_path):
        summary_rows = load_summary_csv(summary_path)
        print("  Summary CSV: {} row(s)".format(len(summary_rows)))
        # Fill in any tensions present in mode CSVs but absent from summary CSV
        csv_tensions = set(r["T_Nm"] for r in summary_rows)
        missing = {t for t, _ in mode_records if t not in csv_tensions}
        if missing:
            print("  Adding {} tension(s) missing from summary CSV: {}".format(
                len(missing), sorted(missing)))
            extra = build_summary_from_modes(
                {k: v for k, v in mode_records.items() if k[0] in missing})
            summary_rows = sorted(summary_rows + extra, key=lambda r: (r["T_Nm"], r["mode"]))
    else:
        print("  baking_buckle_summary.csv not found — building from mode CSVs")
        summary_rows = build_summary_from_modes(mode_records)

    html = build_html(summary_rows, mode_records)
    os.makedirs(RESULTS_DIR, exist_ok=True)
    with open(OUTPUT_HTML, "w", encoding="utf-8") as f:
        f.write(html)
    print("\nReport saved: {}".format(OUTPUT_HTML))
    print("Open in any browser.")


main()