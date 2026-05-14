# FreeCAD — P2 Basket Bulk Micromegas FEM Simulations

Finite-element analysis of the P2 basket Bulk Micromegas PCB detector mesh under
mechanical pre-tension, plus buckling analysis for the baking cycle.

All simulations use **CalculiX** as the solver (bundled with FreeCAD 0.20).
Meshing uses either the FreeCAD GUI built-in mesher or **Gmsh** via the FreeCAD API.
Post-processing scripts are plain Python and require no FreeCAD installation.

---

## Prerequisites

| Requirement | Version | Notes |
|-------------|---------|-------|
| FreeCAD | 0.20 | Required only for simulation scripts |
| Python | 3.9+ | Bundled with FreeCAD; also needed standalone for post-processing |
| numpy | any | Post-processing only |
| matplotlib | any | Post-processing only (`matplotlib.tri` for interpolation) |
| python-docx | any | `generate_docx.py` only |

Install post-processing dependencies:
```
pip install numpy matplotlib python-docx
```

---

## Physical model

Pre-tension is applied via a **thermal analogy**: cooling the mesh material below its
stress-free reference temperature creates the target biaxial tensile stress.

```
ΔT = T [N/m] / (E × t × α)  =  T / 25
```

| Parameter | Value |
|-----------|-------|
| Young's modulus E | 1 GPa |
| Thickness t | 0.5 mm |
| Thermal expansion α | 50 µm / (m · K) |
| Stress-free temperature T₀ | 300 K |
| Pillar pitch | 3 mm |
| Wall BC strip | 1 mm from boundary |

Pillar boundary conditions (Uz = 0) are applied by snapping each 3 mm-pitch grid
centre to the nearest mesh node within 1.5 mm.

---

## Scripts

### Simulation — require FreeCADCmd.exe

Run as:
```
"C:\Program Files\FreeCAD 0.20\bin\FreeCADCmd.exe" FreeCAD\<script>.py
```

> **Note:** Scripts must call their entry point unconditionally at module level.
> FreeCADCmd does not set `__name__ == "__main__"`.

| Script | Purpose | Results directory |
|--------|---------|-------------------|
| `freecad_fem_setup.py` | Original basket tension sweep using FreeCAD GUI adaptive mesh. T = 300–1000 N/m in 100 N/m steps. | `results/simulation_with_pillars/` |
| `freecad_basket_gmsh_fem.py` | **Authoritative** basket sweep re-meshed with Gmsh (uniform C3D10, 2–5 mm). T = 300, 700, 1000 N/m. | `results/basket_gmsh_simulation/` |
| `freecad_triangle_fem_setup.py` | Equilateral triangle geometry (~169 000 mm²) with Gmsh mesh. Same workflow as basket Gmsh. T = 300, 700, 1000 N/m. | `results/triangle_simulation/` |
| `freecad_buckling_sweep.py` | Buckling eigenvalue sweep: 4 tensions × 5 glue-stiffness values (k_glue). Uses spring BCs at wall to model adhesive softening. | `results/baking_simulation/` |
| `freecad_baking_buckling.py` | Baking buckling sweep (CLOAD approach). **Known limitation** — CLOAD at wall nodes is absorbed by the first pillar row; global compression field is not correctly modelled. Use `baking_buckle_analytical.py` for correct results. | `results/baking_simulation/` |

### Analytical — plain Python, no FreeCAD

| Script | Purpose |
|--------|---------|
| `baking_buckle_analytical.py` | Correct baking buckling analysis. Computes the critical eigenvalue λ = T / 121 analytically. Buckling threshold ≈ 121 N/m — well below the operational range (300–1000 N/m), so the mesh never buckles during baking. |

### Visualisation & reports — plain Python

| Script | Output | Notes |
|--------|--------|-------|
| `visualize_sweep.py` | `results/…/sweep_report.html` | Interactive Plotly charts for the original basket tension sweep |
| `visualize_triangle_sweep.py` | `results/triangle_simulation/triangle_uz_surface.html` | Interactive 3D Uz surface for triangle results; tension switcher buttons |
| `visualize_comparison.py` | `results/comparison_basket_vs_triangle.html` | Side-by-side 3D Uz comparison: basket Gmsh vs triangle |
| `visualize_baking_buckle.py` | `results/baking_simulation/baking_buckle_report.html` | Buckling mode shapes and eigenvalue summary |
| `generate_report.py` | `results/fem_report.html` | Full technical HTML report for the original simulation |
| `generate_docx.py` | `results/fem_report.docx` | Word document version of the above |
| `generate_fem_comparison_report.py` | `results/fem_comparison_report.html` | **Comprehensive comparison report** for all three models — displacement maps, root-cause analysis, validation |

### Geometry checks & diagnostics — plain Python

| Script | Purpose |
|--------|---------|
| `check_geometry.py` | Visualises node layers and material zones from the original FCStd mesh |
| `buckling_geometry_check.py` | Interactive BC-zone map for the buckling sweep; verifies wall and pillar BC placement |
| `triangle_geometry_check.py` | Draws triangle boundary, wall BC zone, and pillar grid |
| `parse_frd_results.py` | Standalone CalculiX `.frd` reader; writes the same CSV format as the simulation scripts |
| `recover_from_frd.py` | Recovers sweep results from `.frd` files left in `%TEMP%` when a run is interrupted |
| `watch_baking_buckle.py` | Live log monitor — run in a second terminal while a sweep is in progress |

---

## Results directories

All output files (CSV, HTML, log, DOCX) are **gitignored** — they are regenerated
by running the scripts above.

| Directory | Content |
|-----------|---------|
| `results/simulation_with_pillars/` | Original adaptive-mesh basket sweep (T = 300–1000 N/m) |
| `results/simulation_without_pillars/` | Baseline sweep with no pillar BCs (free membrane) |
| `results/basket_gmsh_simulation/` | Gmsh re-mesh basket sweep — **use these as the reference** |
| `results/triangle_simulation/` | Equilateral triangle Gmsh sweep |
| `results/baking_simulation/` | Buckling eigenvalue results (glue sweep + baking sweep) |

---

## Key findings

### Why the original simulation shows a peak at (~317, 289) mm

The FreeCAD GUI adaptive mesher creates a coarse interior mesh (5–12 mm elements
in the centre vs 1–2 mm near the boundary). The pillar BC snapping algorithm
requires a mesh node within 1.5 mm of each 3 mm-pitch grid centre. In the coarse
interior only ~62 % of pillar positions find a node — the rest are silently skipped.

The resulting large unsupported spans (up to ~9 mm) in the central region produce
a spurious deformation peak at the basket's geometric interior. **This is a mesh
artefact, not a physical effect.**

### Authoritative reference values (Basket Gmsh, 100 % pillar coverage)

| T (N/m) | T (N/cm) | max \|Uz\| |
|---------|---------|-----------|
| 300 | 3 | 0.52 µm |
| 700 | 7 | 1.22 µm |
| 1000 | 10 | 1.74 µm |

The original adaptive-mesh result at T = 1000 N/m is **31.96 µm** — approximately
**18× too high** due to the missing pillar BCs.

### Baking buckling

Analytical critical eigenvalue: **λ = T / 121**

At all operational tensions (300–1000 N/m), λ > 1 — the mesh **does not buckle**
during a standard baking cycle (25 °C → 146 °C). This is consistent with experiment
(10 N/cm has never been observed to buckle).

---

## Original basket geometry

Irregular polygon defined in `C:\Temp\fc_v4_open.FCStd` (local copy):
- Apex: (315, 546) mm
- Base: x ≈ 110–630 mm along y ≈ 0
- Solid body label: `fr4`