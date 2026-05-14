"""
baking_buckle_analytical.py  -  Baking buckling analysis (analytical, no FreeCAD)
P2 Basket Bulk Micromegas mesh

Physics
-------
During baking (25 C -> 146 C, DT_real = 121 K):
  FR4 PCB expands at alpha_FR4 < alpha_SS, dragging ALL Dynamask constraint
  nodes (wall frame + pillar tops) inward relative to the SS mesh.
  Net in-plane compression builds uniformly across the full active area.

  sigma_bake = delta_alpha * DT_real * E_eff          [Pa]
  sigma_pre  = T_Nm / t_eff                            [Pa]
  lambda     = sigma_pre / sigma_bake
             = T_Nm / (E_eff * t_eff * delta_alpha * DT_real)

Interpretation:
  lambda < 1  ->  buckles before end of bake at ~(25 + lambda*121) C
  lambda = 1  ->  buckles exactly at 146 C
  lambda > 1  ->  survives full bake

Why analytical (not FEM):
  CalculiX 2.17 (FreeCAD 0.20) does not support *TEMPERATURE as a reference
  load in *BUCKLE (perturbation) steps — zero eigenvalues are returned.
  The CLOAD-at-wall workaround is also wrong: the 32 643 fixed pillar nodes
  absorb the wall load as reactions, so compression never propagates into the
  interior; results are non-monotonic and inconsistent with experiment.
  For a uniform pre-tensioned membrane with all constraints moving together
  (the correct baking physics), the analytical result is exact.

Run:
    python FreeCAD/baking_buckle_analytical.py

Output:
    FreeCAD/results/baking_buckle_summary.csv
    FreeCAD/results/baking_buckle_report.html
    FreeCAD/baking_buckle.log
"""

import os
import csv
import math

# ---------------------------------------------------------------------------
# Physical parameters
# ---------------------------------------------------------------------------
ALPHA_FR4   = 14e-6    # /K  in-plane CTE of FR4 / S1000-2M PCB  (ESTIMATE)
ALPHA_SS    = 16e-6    # /K  SS316L mesh
DELTA_ALPHA = ALPHA_SS - ALPHA_FR4   # 2e-6 /K  differential

E_EFF  = 1.0e9    # Pa    effective Young's modulus
T_EFF  = 0.5e-3   # m     effective thickness

T_ROOM_C     = 25.0    # C  ambient
T_BAKE_C     = 146.0   # C  baking temperature
DT_BAKE_REAL = T_BAKE_C - T_ROOM_C   # 121 K

# Compression stress at full bake [Pa]
SIGMA_BAKE = DELTA_ALPHA * DT_BAKE_REAL * E_EFF

# Denominator: T_Nm that gives lambda=1
LAMBDA_DENOM = E_EFF * T_EFF * DELTA_ALPHA * DT_BAKE_REAL   # N/m

# Tensions to sweep
TENSIONS_Nm = [300, 500, 700, 1000]

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(SCRIPT_DIR, "results")
LOG_FILE    = os.path.join(SCRIPT_DIR, "baking_buckle.log")
CSV_FILE    = os.path.join(RESULTS_DIR, "baking_buckle_summary.csv")
HTML_FILE   = os.path.join(RESULTS_DIR, "baking_buckle_report.html")

PLOTLY_CDN  = "https://cdn.plot.ly/plotly-2.35.2.min.js"


# ---------------------------------------------------------------------------
# Analytical computation
# ---------------------------------------------------------------------------

def compute_lambda(T_Nm):
    return T_Nm / LAMBDA_DENOM


def buckle_temperature(lam):
    return T_ROOM_C + lam * DT_BAKE_REAL


# ---------------------------------------------------------------------------
# CSV output
# ---------------------------------------------------------------------------

def write_summary_csv(rows, path):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["T_Nm", "T_Ncm", "mode", "eigenvalue", "T_buckle_C", "survived"])
        for r in rows:
            w.writerow([
                r["T_Nm"],
                "{:.1f}".format(r["T_Nm"] / 100.0),
                1,
                "{:.6f}".format(r["lambda"]),
                "{:.1f}".format(r["T_buckle_C"]),
                "yes" if r["survived"] else "no",
            ])


# ---------------------------------------------------------------------------
# HTML report
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
    return '"' + str(v) + '"'


def render_table(rows):
    header = ("<tr>"
              "<th>T (N/m)</th><th>T (N/cm)</th>"
              "<th>λ (analytical)</th><th>T_buckle (°C)</th><th>Status</th>"
              "</tr>")
    body = ""
    for r in rows:
        survived = r["survived"]
        row_bg   = "" if survived else ' style="background:#ffe0e0"'
        status   = ("SURVIVED" if survived
                    else "BUCKLES at ~{:.0f} °C".format(r["T_buckle_C"]))
        color    = ("color:#2ca02c;font-weight:bold" if survived
                    else "color:#d62728;font-weight:bold")
        body += (
            "<tr{bg}>"
            "<td>{t:.0f}</td>"
            "<td>{ncm:.1f}</td>"
            "<td><strong>{lam:.4f}</strong></td>"
            "<td>{tb:.1f}</td>"
            '<td style="{col}">{st}</td>'
            "</tr>\n"
        ).format(
            bg=row_bg,
            t=r["T_Nm"], ncm=r["T_Nm"] / 100.0,
            lam=r["lambda"], tb=r["T_buckle_C"],
            col=color, st=status,
        )
    return "<table><thead>{}</thead><tbody>{}</tbody></table>".format(header, body)


def build_html(rows):
    tensions  = [r["T_Nm"] for r in rows]
    lambdas   = [r["lambda"] for r in rows]
    t_buckles = [r["T_buckle_C"] for r in rows]
    colors    = ["#2ca02c" if r["survived"] else "#d62728" for r in rows]
    hover     = [
        "T={:.0f} N/m ({:.1f} N/cm)<br>λ={:.4f}<br>T_buckle≈{:.0f} °C<br>{}".format(
            r["T_Nm"], r["T_Nm"]/100.0, r["lambda"], r["T_buckle_C"],
            "SURVIVED" if r["survived"] else "BUCKLES"
        )
        for r in rows
    ]

    bar_trace = _jv({
        "type": "bar",
        "name": "λ (analytical)",
        "x": ["{:.0f} N/m".format(t) for t in tensions],
        "y": lambdas,
        "marker": {"color": colors},
        "hovertext": hover,
        "hovertemplate": "%{hovertext}<extra></extra>",
    })

    scatter_trace = _jv({
        "type": "scatter",
        "name": "T_buckle (°C)",
        "x": ["{:.0f} N/m".format(t) for t in tensions],
        "y": t_buckles,
        "mode": "lines+markers",
        "yaxis": "y2",
        "line": {"color": "#7f7f7f", "dash": "dot", "width": 1},
        "marker": {"size": 6, "color": "#7f7f7f"},
        "hovertemplate": "T_buckle: %{y:.1f} °C<extra></extra>",
    })

    summary_table  = render_table(rows)
    tensions_str   = ", ".join("{:.0f}".format(t) for t in tensions)
    lambda_denom_v = "{:.1f}".format(LAMBDA_DENOM)

    return """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>P2 Basket — Baking Buckling Analysis</title>
<script src="{cdn}"></script>
<style>
  body   {{ font-family: Arial, sans-serif; margin: 24px; background: #f7f7f7; }}
  h1     {{ color: #2c3e50; margin-bottom: 4px; }}
  .sub   {{ color: #666; margin-top: 0; font-size: 0.95em; }}
  h2     {{ color: #34495e; margin-top: 36px; border-bottom: 1px solid #ccc; padding-bottom: 4px; }}
  .box   {{ background: white; border-radius: 6px;
            box-shadow: 0 1px 4px #ccc; margin-bottom: 28px; padding: 14px; }}
  table  {{ border-collapse: collapse; width: 100%; max-width: 680px; }}
  th, td {{ border: 1px solid #ddd; padding: 7px 14px; text-align: center; }}
  th     {{ background: #2c3e50; color: white; }}
  tr:nth-child(even) {{ background: #f4f4f4; }}
  .note  {{ font-size: 0.83em; color: #888; margin-top: 8px; }}
  .warn  {{ background: #fff8e1; border-left: 4px solid #f9a825;
            padding: 10px 14px; margin-bottom: 20px;
            border-radius: 0 4px 4px 0; font-size: 0.9em; }}
</style>
</head>
<body>

<h1>P2 Basket Bulk Micromegas — Baking Buckling Analysis</h1>
<p class="sub">
  Tensions: {tensions_str} N/m &nbsp;|&nbsp;
  Bake: 25 °C → 146 °C (ΔT = 121 K) &nbsp;|&nbsp;
  Δα = α<sub>SS</sub> − α<sub>FR4</sub> = 2×10<sup>−6</sup> /K &nbsp;|&nbsp;
  λ = 1 ↔ buckles at exactly 146 °C
</p>

<div class="warn">
  <strong>Analytical result</strong> — CalculiX 2.17 (FreeCAD 0.20) does not support
  <code>*TEMPERATURE</code> as a reference load in <code>*BUCKLE</code> steps.
  For a uniform pre-tensioned membrane with all constraints moving together
  (the correct baking physics), the formula λ = T / {denom} N/m is exact.
</div>

<h2>Eigenvalue vs. pre-tension</h2>
<div id="ev_bar" class="box" style="height:420px"></div>

<h2>Summary table</h2>
<div class="box">
{summary_table}
<p class="note">
  λ &lt; 1 → buckles before end of bake at ≈ (25 + λ×121) °C &nbsp;|&nbsp;
  λ = 1 → buckles exactly at 146 °C &nbsp;|&nbsp;
  λ &gt; 1 → survives full bake<br>
  Buckling threshold: ~{denom} N/m (~{denom_ncm:.1f} N/cm) — well below the
  tested range. Consistent with experiment (10 N/cm never buckles).
</p>
</div>

<h2>Model parameters</h2>
<div class="box">
<table style="max-width:520px">
<thead><tr><th>Parameter</th><th>Value</th></tr></thead>
<tbody>
<tr><td>α<sub>SS</sub></td><td>16×10<sup>−6</sup> /K</td></tr>
<tr><td>α<sub>FR4</sub></td><td>14×10<sup>−6</sup> /K (ESTIMATE — Shengyi S1000-2M)</td></tr>
<tr><td>Δα</td><td>2×10<sup>−6</sup> /K</td></tr>
<tr><td>E<sub>eff</sub></td><td>1.0 GPa</td></tr>
<tr><td>t<sub>eff</sub></td><td>0.5 mm</td></tr>
<tr><td>ΔT<sub>bake</sub></td><td>121 K (25 → 146 °C)</td></tr>
<tr><td>σ<sub>bake</sub> = Δα·ΔT·E</td><td>242 kPa</td></tr>
<tr><td>λ = T / (E·t·Δα·ΔT)</td><td>T [N/m] / {denom} N/m</td></tr>
</tbody>
</table>
</div>

<script>
Plotly.newPlot("ev_bar", [{bar_trace}, {scatter_trace}], {{
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
    text: "λ = 1 (146 °C)",
    showarrow: false,
    xanchor: "right",
    font: {{ size: 11 }},
  }}],
  legend: {{ x: 0.01, y: 0.99, bgcolor: "rgba(255,255,255,0.8)" }},
  margin: {{ t: 20, b: 60, l: 60, r: 90 }},
  paper_bgcolor: "white",
  plot_bgcolor: "white",
}}, {{responsive: true}});
</script>

</body>
</html>
""".format(
        cdn=PLOTLY_CDN,
        tensions_str=tensions_str,
        denom=lambda_denom_v,
        denom_ncm=LAMBDA_DENOM / 100.0,
        summary_table=summary_table,
        bar_trace=bar_trace,
        scatter_trace=scatter_trace,
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)

    with open(LOG_FILE, "w", encoding="utf-8") as log:
        def info(msg):
            print(msg)
            log.write(msg + "\n")
            log.flush()

        info("=" * 60)
        info("P2 Basket FEM — baking buckling (analytical)")
        info("ALPHA_FR4:    {:.1e} /K  (ESTIMATE — update if Shengyi publishes)".format(ALPHA_FR4))
        info("ALPHA_SS:     {:.1e} /K".format(ALPHA_SS))
        info("delta_alpha:  {:.1e} /K".format(DELTA_ALPHA))
        info("DT_bake_real: {:.1f} K  ({:.0f} C -> {:.0f} C)".format(
            DT_BAKE_REAL, T_ROOM_C, T_BAKE_C))
        info("sigma_bake:   {:.0f} Pa  ({:.0f} kPa)".format(SIGMA_BAKE, SIGMA_BAKE / 1e3))
        info("lambda_denom: {:.2f} N/m  (T giving lambda=1)".format(LAMBDA_DENOM))
        info("=" * 60)
        info("")

        rows = []
        for T_Nm in TENSIONS_Nm:
            lam      = compute_lambda(T_Nm)
            t_buck   = buckle_temperature(lam)
            survived = lam >= 1.0
            info("  T={:4d} N/m ({:.1f} N/cm)  lambda={:.4f}  T_buckle~{:.0f} C  [{}]".format(
                T_Nm, T_Nm / 100.0, lam, t_buck,
                "SURVIVED" if survived else "BUCKLES"))
            rows.append({
                "T_Nm":       T_Nm,
                "lambda":     lam,
                "T_buckle_C": t_buck,
                "survived":   survived,
            })

        write_summary_csv(rows, CSV_FILE)
        info("\nSummary CSV: " + CSV_FILE)

        html = build_html(rows)
        with open(HTML_FILE, "w", encoding="utf-8") as f:
            f.write(html)
        info("HTML report: " + HTML_FILE)

        info("")
        info("Buckling threshold: {:.1f} N/m ({:.2f} N/cm)".format(
            LAMBDA_DENOM, LAMBDA_DENOM / 100.0))
        info("All tested tensions survive full bake (lambda > 1).")
        info("=" * 60)


main()