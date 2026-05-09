"""
freecad_fem_setup.py  –  Tension sweep for P2 basket Bulk Micromegas PCB

Runs inside FreeCAD's headless interpreter:
    "C:\Program Files\FreeCAD 0.20\bin\FreeCADCmd.exe" freecad_fem_setup.py

What it does
------------
For each tension value T in TENSION_SWEEP_Nm (N/m):
  1. Converts T → equivalent thermal ΔT using the mesh material parameters
  2. Updates ConstraintTemperature in the loaded FreeCAD document
  3. Runs the CalculiX solver
  4. Reads node coordinates + DisplacementVectors + vonMises from the result
  5. Saves  results/sweep_T<XXXX>Nm.csv
After the sweep saves results/sweep_summary.csv

The geometry discrepancies vs. target (PCB 400 µm→200 µm, gap 128 µm→150 µm)
are noted below; fix the Pad lengths before the final production run.

Tension-to-temperature conversion
----------------------------------
The mesh material simulates pre-tension thermally:
  ε = α_mesh × |ΔT|  =  T / (E_mesh × t_mesh)
  → ΔT = T / (E_mesh × t_mesh × α_mesh)

Model v4 parameters (from Document.xml):
  E_mesh   = 1.0e9  Pa   (1 GPa)
  t_mesh   = 0.5e-3 m    (0.5 mm effective thickness)
  α_mesh   = 50e-6  /K   (50 µm/m/K)
  T_ref    = 300 K        (ConstraintInitialTemperature)
  T_applied = T_ref - ΔT  (applied to mesh surface, always < T_ref → contraction)

  → ΔT_per_Nm = 1 / (1e9 × 0.5e-3 × 50e-6) = 1/25  K·m/N
"""

import os
import sys
import shutil
import csv
import traceback

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

FCSTD_TEMPLATE = (
    r"C:\Users\alexa\OneDrive - "
    r"Αριστοτέλειο"
    r" Πανεπιστήμιο"
    r" Θεσσαλονίκης"
    r"\Documents\FreeCAD"
    r"\p2_basket_detector_mechanical_tension_v4.FCStd"
)
# Fallback to the local copy we already made
FCSTD_LOCAL_COPY = r"C:\Temp\fc_v4_open.FCStd"

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(SCRIPT_DIR, "results")
WORK_DIR    = r"C:\Temp\sweep_work"
LOG_FILE    = os.path.join(SCRIPT_DIR, "sweep_run.log")

# Tension sweep: 3 N/cm → 10 N/cm in 1 N/cm steps  (1 N/cm = 100 N/m)
TENSION_SWEEP_Nm = [300, 400, 500, 600, 700, 800, 900, 1000]

# Mesh material constants (from Document.xml)
E_MESH   = 1.0e9    # Pa
T_MESH   = 0.5e-3   # m  (effective thickness)
ALPHA    = 50.0e-6  # /K
T_REF    = 300.0    # K  (ConstraintInitialTemperature)

# FreeCAD object names (from Document.xml object list)
OBJ_ANALYSIS  = "Analysis"
OBJ_SOLVER    = "SolverCcxTools"
OBJ_TEMP_BC   = "ConstraintTemperature"
OBJ_RESULT    = "CCX_Results"
OBJ_RES_MESH  = "ResultMesh"

# ---------------------------------------------------------------------------
# Boundary-condition injection parameters
# ---------------------------------------------------------------------------
MESH_Z_LO_MM    = 0.05   # z lower bound for mesh-layer nodes (mm)
MESH_Z_HI_MM    = 0.70   # z upper bound
PILLAR_PITCH_MM = 3.0    # Dynamask pillar pitch, centre-to-centre (mm)
WALL_MM         = 1.0    # peripheral mesh-frame wall thickness (mm)
PILLAR_TOL_MM   = 1.5    # max snap distance: pillar centre → nearest node


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def tension_to_delta_T(T_Nm):
    """Convert mesh tension (N/m) to required thermal ΔT (K)."""
    return T_Nm / (E_MESH * T_MESH * ALPHA)


def tension_to_applied_T(T_Nm):
    """Temperature to set on ConstraintTemperature (K)."""
    return T_REF - tension_to_delta_T(T_Nm)


def open_log():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    log = open(LOG_FILE, "w", encoding="utf-8", buffering=1)
    return log


def info(log, msg):
    log.write(msg + "\n")
    try:
        print(msg)
        sys.stdout.flush()
    except Exception:
        pass


# ---------------------------------------------------------------------------
# FreeCAD solver runner
# ---------------------------------------------------------------------------

def run_one_step(doc, T_Nm, log, wall_nodes=None, pillar_nodes=None):
    """
    Update ConstraintTemperature, inject BCs, run CalculiX, return result arrays.

    Parameters
    ----------
    wall_nodes   : set[int] or None  — peripheral frame nodes (Ux=Uy=Uz=0)
    pillar_nodes : set[int] or None  — pillar-snapped nodes (Uz=0)

    Returns dict with keys:
        node_ids   : list[int]
        x, y, z    : list[float]  node coordinates (mm)
        Ux, Uy, Uz : list[float]  displacements (mm)
        von_mises  : list[float]  von Mises stress (MPa)
    """
    import FreeCAD
    from femtools import ccxtools

    # 1. Update thermal BC (tension analogy)
    T_applied = tension_to_applied_T(T_Nm)
    delta_T   = tension_to_delta_T(T_Nm)
    info(log, "  T={} N/m  →  ΔT={:.2f} K  →  T_applied={:.2f} K".format(
        T_Nm, delta_T, T_applied))

    bc = doc.getObject(OBJ_TEMP_BC)
    bc.Temperature = T_applied
    doc.recompute()

    # 2. Prepare CalculiX solver
    analysis = doc.getObject(OBJ_ANALYSIS)
    solver   = doc.getObject(OBJ_SOLVER)
    fea = ccxtools.FemToolsCcx(analysis, solver)
    fea.update_objects()
    fea.setup_working_dir()
    fea.setup_ccx()

    # 3. Inject wall + pillar displacement BCs into the .inp file
    if wall_nodes or pillar_nodes:
        try:
            inject_bcs_into_inp(fea, wall_nodes or set(), pillar_nodes or set(), log)
        except Exception as exc:
            info(log, "  WARNING: BC injection failed — {}".format(exc))

    info(log, "  Running CalculiX ...")
    message = fea.run()
    # fea.run() returns True on success in FreeCAD 0.20, or a non-empty string on error
    if message and not isinstance(message, bool):
        raise RuntimeError("CalculiX returned error:\n" + str(message))

    info(log, "  Loading results ...")
    fea.load_results()

    # 4. Extract results
    result = doc.getObject(OBJ_RESULT)
    res_mesh = doc.getObject(OBJ_RES_MESH)

    if result is None:
        raise RuntimeError("CCX_Results object not found after solve")

    node_numbers = list(result.NodeNumbers)   # list[int] - node ID order
    disp_vecs    = result.DisplacementVectors  # list[FreeCAD.Vector]
    von_mises    = list(result.vonMises)       # list[float], Pa

    if len(node_numbers) == 0:
        raise RuntimeError("NodeNumbers is empty - results not loaded")

    # Build node coordinate lookup from result mesh
    mesh_nodes = res_mesh.FemMesh.Nodes  # dict: node_id -> FreeCAD.Vector (mm)

    nids, xs, ys, zs = [], [], [], []
    uxs, uys, uzs    = [], [], []
    vms               = []

    for i, nid in enumerate(node_numbers):
        coord = mesh_nodes.get(nid)
        if coord is None:
            continue
        disp = disp_vecs[i]
        nids.append(nid)
        xs.append(round(float(coord.x), 4))
        ys.append(round(float(coord.y), 4))
        zs.append(round(float(coord.z), 4))
        uxs.append(round(float(disp.x), 6))
        uys.append(round(float(disp.y), 6))
        uzs.append(round(float(disp.z), 6))
        # CalculiX outputs stress in MPa when geometry is in mm
        vms.append(round(float(von_mises[i]), 4))

    info(log, "  Extracted {} nodes".format(len(nids)))

    return {
        "node_ids": nids,
        "x": xs, "y": ys, "z": zs,
        "Ux": uxs, "Uy": uys, "Uz": uzs,
        "von_mises_MPa": vms,
    }


# ---------------------------------------------------------------------------
# CSV writers
# ---------------------------------------------------------------------------

def save_step_csv(T_Nm, data, results_dir):
    fname = os.path.join(results_dir, "sweep_T{:04d}Nm.csv".format(T_Nm))
    with open(fname, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["node_id", "x_mm", "y_mm", "z_mm",
                    "Ux_mm", "Uy_mm", "Uz_mm", "vonMises_MPa"])
        for i in range(len(data["node_ids"])):
            w.writerow([
                data["node_ids"][i],
                data["x"][i], data["y"][i], data["z"][i],
                data["Ux"][i], data["Uy"][i], data["Uz"][i],
                data["von_mises_MPa"][i],
            ])
    return fname


def save_summary_csv(summary_rows, results_dir):
    fname = os.path.join(results_dir, "sweep_summary.csv")
    with open(fname, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["T_Nm", "T_Ncm", "delta_T_K",
                    "max_Uz_mm", "min_Uz_mm",
                    "max_abs_Uz_mm", "max_vonMises_MPa"])
        for row in summary_rows:
            w.writerow(row)
    return fname


# ---------------------------------------------------------------------------
# BC injection helpers  (pure Python, no external deps)
# ---------------------------------------------------------------------------

def _hull_2d(xy_pairs):
    """Andrew's monotone-chain convex hull."""
    pts = sorted(set((round(x, 2), round(y, 2)) for x, y in xy_pairs))
    if len(pts) < 3:
        return pts

    def _cross(O, A, B):
        return (A[0] - O[0]) * (B[1] - O[1]) - (A[1] - O[1]) * (B[0] - O[0])

    lower = []
    for p in pts:
        while len(lower) >= 2 and _cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper = []
    for p in reversed(pts):
        while len(upper) >= 2 and _cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]


def _pip(px, py, poly):
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


def _shrink_hull(hull, amount):
    """Move each vertex toward the centroid by `amount` mm."""
    import math
    cx = sum(p[0] for p in hull) / len(hull)
    cy = sum(p[1] for p in hull) / len(hull)
    out = []
    for x, y in hull:
        dx, dy = x - cx, y - cy
        dist = math.sqrt(dx * dx + dy * dy)
        if dist > amount:
            s = (dist - amount) / dist
            out.append((cx + dx * s, cy + dy * s))
    return out


def build_bc_node_sets(fem_mesh, log):
    """
    Derive wall-BC and pillar-BC node sets from the FreeCAD FemMesh object.

    Returns
    -------
    wall_nodes   : set[int]  nodes on the peripheral 1-mm frame (Ux=Uy=Uz=0)
    pillar_nodes : set[int]  one node per pillar centre, snapped to nearest
                             interior mesh node (Uz=0)
    """
    import math

    # --- collect mesh-layer nodes -------------------------------------------
    all_nodes = fem_mesh.Nodes          # dict: node_id -> FreeCAD.Vector (mm)
    mesh_nodes = {}
    for nid, v in all_nodes.items():
        if MESH_Z_LO_MM <= v.z <= MESH_Z_HI_MM:
            mesh_nodes[nid] = (v.x, v.y, v.z)

    if not mesh_nodes:
        raise RuntimeError(
            "No nodes found in mesh-layer z=[{}, {}] mm".format(
                MESH_Z_LO_MM, MESH_Z_HI_MM))
    info(log, "  Mesh-layer nodes: {}".format(len(mesh_nodes)))

    # --- convex hull of mesh layer ------------------------------------------
    xy_pairs = [(x, y) for x, y, z in mesh_nodes.values()]
    hull_outer = _hull_2d(xy_pairs)
    hull_inner = _shrink_hull(hull_outer, WALL_MM)
    info(log, "  Hull vertices: {} outer → {} inner (after {}-mm shrink)".format(
        len(hull_outer), len(hull_inner), WALL_MM))

    xs = [p[0] for p in xy_pairs]
    ys = [p[1] for p in xy_pairs]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)

    # --- wall nodes: mesh-layer nodes outside the inner hull ----------------
    wall_nodes = set()
    for nid, (x, y, z) in mesh_nodes.items():
        if not _pip(x, y, hull_inner):
            wall_nodes.add(nid)
    info(log, "  Wall BC nodes: {}".format(len(wall_nodes)))

    # --- spatial hash of interior nodes for fast nearest-node lookup --------
    cell = PILLAR_PITCH_MM
    grid = {}
    for nid, (x, y, z) in mesh_nodes.items():
        if nid in wall_nodes:
            continue
        key = (int(x / cell), int(y / cell))
        grid.setdefault(key, []).append((nid, x, y))

    # --- pillar nodes: nearest interior node to each pillar centre ----------
    pillar_nodes = set()
    tol2 = PILLAR_TOL_MM ** 2
    px = xmin + WALL_MM + PILLAR_PITCH_MM / 2
    while px <= xmax - WALL_MM:
        py = ymin + WALL_MM + PILLAR_PITCH_MM / 2
        while py <= ymax - WALL_MM:
            if _pip(px, py, hull_inner):
                cx_cell = int(px / cell)
                cy_cell = int(py / cell)
                best_nid, best_d2 = None, tol2
                for di in (-1, 0, 1):
                    for dj in (-1, 0, 1):
                        for nid, nx, ny in grid.get((cx_cell + di, cy_cell + dj), []):
                            d2 = (nx - px) ** 2 + (ny - py) ** 2
                            if d2 < best_d2:
                                best_d2, best_nid = d2, nid
                if best_nid is not None:
                    pillar_nodes.add(best_nid)
            py += PILLAR_PITCH_MM
        px += PILLAR_PITCH_MM

    info(log, "  Pillar BC nodes: {}".format(len(pillar_nodes)))
    return wall_nodes, pillar_nodes


def _find_inp_file(fea):
    """Return path to the CalculiX .inp file produced by setup_ccx()."""
    # FreeCAD 0.20 stores it in fea.inp_file_name
    candidate = getattr(fea, "inp_file_name", None)
    if candidate and os.path.isfile(candidate):
        return candidate
    # Fallback: search working dir for any .inp
    wdir = getattr(fea, "working_dir", None)
    if wdir and os.path.isdir(wdir):
        for fname in os.listdir(wdir):
            if fname.endswith(".inp"):
                return os.path.join(wdir, fname)
    raise RuntimeError("Cannot locate CalculiX .inp file")


def inject_bcs_into_inp(fea, wall_nodes, pillar_nodes, log):
    """
    Write *BOUNDARY cards for wall and pillar constraints into the
    CalculiX .inp file, immediately before the first *STEP card.
    """
    inp_file = _find_inp_file(fea)

    with open(inp_file, "r", encoding="utf-8", errors="replace") as f:
        content = f.read()

    lines = [
        "**",
        "** === Injected BCs: peripheral wall + Dynamask pillar supports ===",
    ]
    if wall_nodes:
        lines += [
            "** Peripheral frame wall — Ux=Uy=Uz=0  ({} nodes)".format(len(wall_nodes)),
            "*BOUNDARY",
        ]
        for nid in sorted(wall_nodes):
            lines.append("{}, 1, 3, 0.0".format(nid))

    if pillar_nodes:
        lines += [
            "** Dynamask pillars — Uz=0  ({} nodes, {} mm pitch)".format(
                len(pillar_nodes), int(PILLAR_PITCH_MM)),
            "*BOUNDARY",
        ]
        for nid in sorted(pillar_nodes):
            lines.append("{}, 3, 3, 0.0".format(nid))

    lines.append("** === end injected BCs ===")
    bc_block = "\n".join(lines) + "\n"

    idx = content.find("*STEP")
    if idx == -1:
        raise RuntimeError("*STEP marker not found in .inp: " + inp_file)
    line_start = content.rfind("\n", 0, idx) + 1
    new_content = content[:line_start] + bc_block + content[line_start:]

    with open(inp_file, "w", encoding="utf-8") as f:
        f.write(new_content)

    info(log, "  BCs injected → {}  ({} wall + {} pillar nodes)".format(
        inp_file, len(wall_nodes), len(pillar_nodes)))


# ---------------------------------------------------------------------------
# Main sweep
# ---------------------------------------------------------------------------

def run_sweep():
    log = open_log()
    info(log, "=" * 60)
    info(log, "P2 Basket FEM tension sweep")
    info(log, "Tensions: {} N/m".format(TENSION_SWEEP_Nm))
    info(log, "=" * 60)

    # Import FreeCAD here (only available inside FreeCADCmd)
    try:
        import FreeCAD
    except ImportError:
        info(log, "ERROR: FreeCAD module not available. "
             "Run this script via FreeCADCmd.exe")
        log.close()
        return

    # Open the template document
    fcstd_path = FCSTD_LOCAL_COPY
    if not os.path.exists(fcstd_path):
        info(log, "Local copy not found, trying OneDrive path ...")
        fcstd_path = FCSTD_TEMPLATE
    if not os.path.exists(fcstd_path):
        info(log, "ERROR: cannot find FCStd file at:\n  {}\n  {}".format(
            FCSTD_LOCAL_COPY, FCSTD_TEMPLATE))
        log.close()
        return

    info(log, "Opening: " + fcstd_path)
    doc = FreeCAD.openDocument(fcstd_path)

    # --- pre-compute BC node sets once (geometry is fixed across all steps) --
    wall_nids = pillar_nids = None
    try:
        info(log, "\nLocating FEM input mesh ...")
        mesh_obj = None
        for obj in doc.Objects:
            if hasattr(obj, "FemMesh") and obj.Name != OBJ_RES_MESH:
                mesh_obj = obj
                info(log, "  Found mesh object: '{}'  ({} nodes)".format(
                    obj.Name, len(obj.FemMesh.Nodes)))
                break
        if mesh_obj is None:
            info(log, "  WARNING: no FEM input mesh found — "
                 "wall/pillar BCs will not be injected")
        else:
            info(log, "Computing wall + pillar BC node sets ...")
            wall_nids, pillar_nids = build_bc_node_sets(mesh_obj.FemMesh, log)
    except Exception as exc:
        info(log, "WARNING: BC node-set computation failed: " + str(exc))
        info(log, traceback.format_exc())

    summary_rows = []
    failed = []

    for T_Nm in TENSION_SWEEP_Nm:
        info(log, "\n--- T = {} N/m  ({:.1f} N/cm) ---".format(
            T_Nm, T_Nm / 100.0))

        try:
            data = run_one_step(doc, T_Nm, log,
                                wall_nodes=wall_nids, pillar_nodes=pillar_nids)

            csv_path = save_step_csv(T_Nm, data, RESULTS_DIR)
            info(log, "  Saved: " + csv_path)

            uzs = data["Uz"]
            max_uz = max(uzs)
            min_uz = min(uzs)
            max_abs_uz = max(abs(v) for v in uzs)
            max_vm = max(data["von_mises_MPa"])

            info(log, "  max|Uz| = {:.4f} mm  ({:.1f} µm)".format(
                max_abs_uz, max_abs_uz * 1000))
            info(log, "  max vonMises = {:.4f} MPa".format(max_vm))

            summary_rows.append([
                T_Nm,
                round(T_Nm / 100.0, 1),
                round(tension_to_delta_T(T_Nm), 2),
                round(max_uz, 6),
                round(min_uz, 6),
                round(max_abs_uz, 6),
                round(max_vm, 4),
            ])

        except Exception as exc:
            info(log, "  FAILED: " + str(exc))
            info(log, traceback.format_exc())
            failed.append(T_Nm)

    # Save summary
    if summary_rows:
        scsv = save_summary_csv(summary_rows, RESULTS_DIR)
        info(log, "\nSummary saved: " + scsv)

    info(log, "\n" + "=" * 60)
    if failed:
        info(log, "FAILED steps: {}".format(failed))
    else:
        info(log, "All {} steps completed successfully.".format(
            len(TENSION_SWEEP_Nm)))
    info(log, "Results directory: " + RESULTS_DIR)
    info(log, "=" * 60)

    log.close()


run_sweep()