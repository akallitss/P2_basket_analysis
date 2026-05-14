"""
freecad_basket_gmsh_fem.py  —  P2 basket re-meshed with Gmsh, direct CalculiX

Re-meshes the ORIGINAL basket solid body (from fc_v4_open.FCStd) using Gmsh
at the same settings as freecad_triangle_fem_setup.py.  Runs CalculiX directly
without using FreeCAD's FEM framework.  Results go to a separate directory so
existing simulation_with_pillars / simulation_without_pillars results are
never touched.

Run inside FreeCAD's headless interpreter:
    "C:/Program Files/FreeCAD 0.20/bin/FreeCADCmd.exe"
        "C:/Users/alexa/PycharmProjects/P2_basket_analysis/FreeCAD/freecad_basket_gmsh_fem.py"

Purpose: compare the original P2 basket geometry against the equilateral
triangle using an identical Gmsh mesh (same element type, same characteristic
length, same direct-CalculiX workflow) so the two geometries can be fairly
compared without mesh-quality artefacts from the FreeCAD GUI mesher.
"""

import os
import sys
import csv
import math
import traceback
import subprocess
import tempfile

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR  = os.path.join(SCRIPT_DIR, "results", "basket_gmsh_simulation")
LOG_FILE     = os.path.join(RESULTS_DIR, "basket_gmsh_sweep.log")

# Original basket FCStd  (geometry only — we create a new Gmsh mesh from it)
BASKET_FCSTD_PATHS = [
    r"C:\Temp\fc_v4_open.FCStd",
]

# New FCStd for the Gmsh-remeshed basket
GMSH_FCSTD_PATH = r"C:\Temp\basket_gmsh_fem.FCStd"

CCX_CANDIDATES = [
    r"C:\Program Files\FreeCAD 0.20\bin\ccx.exe",
    r"C:\Program Files (x86)\FreeCAD 0.20\bin\ccx.exe",
    r"C:\Program Files\FreeCAD\bin\ccx.exe",
]

# ---------------------------------------------------------------------------
# Sweep parameters
# ---------------------------------------------------------------------------

TENSION_SWEEP_Nm = [300, 700, 1000]   # same three tensions as triangle script

# ---------------------------------------------------------------------------
# Mesh material constants  (identical to freecad_fem_setup.py / triangle script)
# ---------------------------------------------------------------------------

E_MESH   = 1.0e9    # Pa
T_MESH   = 0.5e-3   # m
ALPHA    = 50.0e-6  # /K
T_REF    = 300.0    # K  (stress-free reference temperature)
POISSON  = 0.3

# Pillar / wall BC geometry (same as all other scripts)
WALL_MM         = 1.0
PILLAR_PITCH_MM = 3.0
PILLAR_TOL_MM   = 1.5

# Gmsh characteristic length  (match triangle script exactly)
GMSH_LMIN = 2.0
GMSH_LMAX = 5.0

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def info(log, msg):
    print(msg)
    log.write(msg + "\n")
    log.flush()


def find_ccx_exe(log):
    for p in CCX_CANDIDATES:
        if os.path.isfile(p):
            info(log, "CalculiX: " + p)
            return p
    raise RuntimeError("ccx.exe not found — checked: " + str(CCX_CANDIDATES))


def tension_to_applied_T(T_Nm):
    """Return the CalculiX nodal temperature that produces pre-tension T_Nm."""
    eps = T_Nm / (E_MESH * T_MESH)
    delta_T = eps / ALPHA
    return T_REF - delta_T


# ---------------------------------------------------------------------------
# Hull-based geometry helpers  (ported from freecad_fem_setup.py)
# ---------------------------------------------------------------------------

def _hull_2d(xy_pairs):
    """Andrew's monotone-chain convex hull."""
    pts = sorted(set((round(x, 2), round(y, 2)) for x, y in xy_pairs))
    if len(pts) < 3:
        return pts
    def cross(O, A, B):
        return (A[0]-O[0])*(B[1]-O[1]) - (A[1]-O[1])*(B[0]-O[0])
    lo, hi = [], []
    for p in pts:
        while len(lo) >= 2 and cross(lo[-2], lo[-1], p) <= 0: lo.pop()
        lo.append(p)
    for p in reversed(pts):
        while len(hi) >= 2 and cross(hi[-2], hi[-1], p) <= 0: hi.pop()
        hi.append(p)
    return lo[:-1] + hi[:-1]


def _shrink_hull(hull, amount):
    """Move each vertex toward the centroid by `amount` mm."""
    cx = sum(p[0] for p in hull) / len(hull)
    cy = sum(p[1] for p in hull) / len(hull)
    out = []
    for x, y in hull:
        dx, dy = x - cx, y - cy
        dist = math.sqrt(dx*dx + dy*dy)
        if dist < 1e-9:
            out.append((x, y))
        else:
            f = max(0.0, dist - amount) / dist
            out.append((cx + dx*f, cy + dy*f))
    return out


def _pip(px, py, poly):
    """Ray-casting point-in-polygon test."""
    inside = False
    n = len(poly)
    j = n - 1
    for i in range(n):
        xi, yi = poly[i]
        xj, yj = poly[j]
        if ((yi > py) != (yj > py) and
                px < (xj - xi) * (py - yi) / (yj - yi) + xi):
            inside = not inside
        j = i
    return inside


# ---------------------------------------------------------------------------
# FreeCAD document creation with Gmsh mesh
# ---------------------------------------------------------------------------

def _find_basket_fcstd():
    for p in BASKET_FCSTD_PATHS:
        if os.path.isfile(p):
            return p
    return None


def create_basket_gmsh_fcstd(save_path, log):
    """
    Open the original basket FCStd, extract the solid body, re-mesh it with
    Gmsh (C3D10, CharacteristicLengthMax=5 mm), save to save_path.
    """
    import FreeCAD
    import Part
    import ObjectsFem
    from femmesh import gmshtools

    src = _find_basket_fcstd()
    if src is None:
        raise RuntimeError(
            "Cannot find original basket FCStd at:\n  " +
            "\n  ".join(BASKET_FCSTD_PATHS))
    info(log, "  Opening original basket FCStd: " + src)
    orig_doc = FreeCAD.openDocument(src)

    # Find the first object that has a non-empty solid shape
    solid_shape = None
    for obj in orig_doc.Objects:
        if hasattr(obj, "Shape") and len(obj.Shape.Solids) > 0:
            solid_shape = obj.Shape.copy()
            info(log, "  Found solid body: '{}' ({} solids)".format(
                obj.Label, len(obj.Shape.Solids)))
            break
    FreeCAD.closeDocument(orig_doc.Name)

    if solid_shape is None:
        raise RuntimeError("No solid body found in the basket FCStd")

    # New document — solid body + Gmsh mesh only
    doc = FreeCAD.newDocument("BasketGmsh")
    part = doc.addObject("Part::Feature", "BasketSolid")
    part.Shape = solid_shape
    doc.recompute()

    mesh_obj = ObjectsFem.makeMeshGmsh(doc, "FEMMeshGmsh")
    mesh_obj.Part = part
    mesh_obj.CharacteristicLengthMin = GMSH_LMIN
    mesh_obj.CharacteristicLengthMax = GMSH_LMAX
    mesh_obj.ElementOrder = "2nd"   # C3D10
    doc.recompute()

    info(log, "  Running Gmsh (char. length {}-{} mm, C3D10) ...".format(
        GMSH_LMIN, GMSH_LMAX))
    mr = gmshtools.GmshTools(mesh_obj)
    mr.create_mesh()
    doc.recompute()

    n_nodes = len(mesh_obj.FemMesh.Nodes)
    n_vols  = len(mesh_obj.FemMesh.Volumes)
    if n_nodes == 0:
        raise RuntimeError("Gmsh produced no nodes")
    info(log, "  Mesh: {} nodes, {} volume elements".format(n_nodes, n_vols))

    doc.saveAs(save_path)
    info(log, "  Saved: " + save_path)
    return doc


# ---------------------------------------------------------------------------
# Mesh data extraction
# ---------------------------------------------------------------------------

def _tet_signed_vol6(nodes, n0, n1, n2, n3):
    x0, y0, z0 = nodes[n0]
    x1, y1, z1 = nodes[n1]
    x2, y2, z2 = nodes[n2]
    x3, y3, z3 = nodes[n3]
    ax, ay, az = x1-x0, y1-y0, z1-z0
    bx, by, bz = x2-x0, y2-y0, z2-z0
    cx, cy, cz = x3-x0, y3-y0, z3-z0
    return (ax*(by*cz - bz*cy) - ay*(bx*cz - bz*cx) + az*(bx*cy - by*cx))


def extract_mesh_data(mesh_obj, log):
    fem_mesh = mesh_obj.FemMesh
    nodes = {}
    for nid, v in fem_mesh.Nodes.items():
        nodes[nid] = (v.x, v.y, v.z)

    raw = {}
    for eid in fem_mesh.Volumes:
        en = fem_mesh.getElementNodes(eid)
        raw[eid] = (en, _tet_signed_vol6(nodes, en[0], en[1], en[2], en[3]))

    pos = sum(1 for _, v in raw.values() if v > 1e-10)
    neg = sum(1 for _, v in raw.values() if v < -1e-10)
    info(log, "  Element vol sign: {} positive, {} negative — majority sign = {}".format(
        pos, neg, 1 if pos >= neg else -1))
    expected_sign = 1 if pos >= neg else -1

    elements = {}
    for eid, (en, vol6) in raw.items():
        if vol6 * expected_sign >= 1e-10:
            elements[eid] = en

    n_vert = len(next(iter(elements.values())))
    etype_map = {4: "C3D4", 10: "C3D10", 8: "C3D8", 20: "C3D20R"}
    info(log, "  Mesh nodes: {}  Volume elements: {}".format(
        len(nodes), len(elements)))
    info(log, "  Element type: {} ({}-node)".format(
        etype_map.get(n_vert, "unknown"), n_vert))
    return nodes, elements


# ---------------------------------------------------------------------------
# BC node-set builder  (hull-based — works for any convex polygon)
# ---------------------------------------------------------------------------

def build_bc_node_sets(mesh_nodes, log):
    """
    Wall nodes   — outside the hull shrunk by WALL_MM     → Ux=Uy=Uz=0
    Pillar nodes — nearest interior node per pillar centre → Uz=0
    """
    # Project all nodes to xy (use all z layers)
    xy_pairs = [(x, y) for x, y, z in mesh_nodes.values()]
    hull_outer = _hull_2d(xy_pairs)
    hull_inner = _shrink_hull(hull_outer, WALL_MM)
    info(log, "  Hull: {} outer vertices → {} inner vertices".format(
        len(hull_outer), len(hull_inner)))

    xs = [p[0] for p in xy_pairs]
    ys = [p[1] for p in xy_pairs]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)

    # Wall nodes: any node NOT inside the inner hull
    wall_nodes = {nid for nid, (x, y, _) in mesh_nodes.items()
                  if not _pip(x, y, hull_inner)}
    info(log, "  Wall BC nodes: {}".format(len(wall_nodes)))

    # Spatial hash of interior nodes
    cell = PILLAR_PITCH_MM
    grid = {}
    for nid, (x, y, _) in mesh_nodes.items():
        if nid not in wall_nodes:
            grid.setdefault((int(x / cell), int(y / cell)), []).append((nid, x, y))

    pillar_nodes = set()
    tol2 = PILLAR_TOL_MM ** 2
    pcx = xmin + WALL_MM + PILLAR_PITCH_MM / 2.0
    while pcx <= xmax - WALL_MM:
        pcy = ymin + WALL_MM + PILLAR_PITCH_MM / 2.0
        while pcy <= ymax - WALL_MM:
            if _pip(pcx, pcy, hull_inner):
                cx_cell = int(pcx / cell)
                cy_cell = int(pcy / cell)
                best_nid, best_d2 = None, tol2
                for di in (-1, 0, 1):
                    for dj in (-1, 0, 1):
                        for nid, nx, ny in grid.get(
                                (cx_cell + di, cy_cell + dj), []):
                            d2 = (nx - pcx)**2 + (ny - pcy)**2
                            if d2 < best_d2:
                                best_d2, best_nid = d2, nid
                if best_nid is not None:
                    pillar_nodes.add(best_nid)
            pcy += PILLAR_PITCH_MM
        pcx += PILLAR_PITCH_MM

    info(log, "  Pillar BC nodes: {}".format(len(pillar_nodes)))
    return wall_nodes, pillar_nodes


# ---------------------------------------------------------------------------
# CalculiX .inp writer
# ---------------------------------------------------------------------------

def write_inp(work_dir, mesh_nodes, mesh_elements, T_Nm, wall_nodes, pillar_nodes):
    n_vert = len(next(iter(mesh_elements.values())))
    etype  = {4: "C3D4", 10: "C3D10", 8: "C3D8", 20: "C3D20R"}.get(n_vert, "C3D4")
    T_applied = tension_to_applied_T(T_Nm)

    inp_path = os.path.join(work_dir, "basket.inp")
    with open(inp_path, "w") as f:
        f.write("** freecad_basket_gmsh_fem.py  T={} N/m\n\n".format(T_Nm))

        f.write("*NODE, NSET=NALL\n")
        for nid, (x, y, z) in sorted(mesh_nodes.items()):
            f.write("{}, {:.6f}, {:.6f}, {:.6f}\n".format(nid, x, y, z))

        # FreeCAD/Gmsh → CalculiX winding fix (same as triangle script)
        f.write("*ELEMENT, TYPE={}, ELSET=EALL\n".format(etype))
        for eid, n in sorted(mesh_elements.items()):
            if len(n) == 4:
                out = (n[0], n[2], n[1], n[3])
            elif len(n) == 10:
                out = (n[0], n[2], n[1], n[3], n[6], n[5], n[4], n[7], n[9], n[8])
            else:
                out = n
            f.write("{}, {}\n".format(eid, ", ".join(str(x) for x in out)))

        f.write("*MATERIAL, NAME=BASKETMESH\n")
        f.write("*ELASTIC\n")
        f.write("{:.1f}, {:.2f}\n".format(E_MESH / 1e6, POISSON))
        f.write("*EXPANSION, ZERO={:.1f}\n".format(T_REF))
        f.write("{:.6E}\n".format(ALPHA))
        f.write("*SOLID SECTION, ELSET=EALL, MATERIAL=BASKETMESH\n\n")

        f.write("** Wall BC (1 mm from boundary): Ux=Uy=Uz=0\n")
        f.write("*BOUNDARY\n")
        for nid in sorted(wall_nodes):
            f.write("{}, 1, 3, 0.0\n".format(nid))

        if pillar_nodes:
            f.write("** Pillar BC ({} mm pitch): Uz=0\n".format(int(PILLAR_PITCH_MM)))
            f.write("*BOUNDARY\n")
            for nid in sorted(pillar_nodes):
                f.write("{}, 3, 3, 0.0\n".format(nid))

        f.write("*INITIAL CONDITIONS, TYPE=TEMPERATURE\n")
        f.write("NALL, {:.1f}\n\n".format(T_REF))
        f.write("*STEP\n")
        f.write("*STATIC\n")
        f.write("*TEMPERATURE\n")
        f.write("NALL, {:.6f}\n".format(T_applied))
        f.write("*NODE FILE\n")
        f.write("U\n")
        f.write("*EL FILE\n")
        f.write("S\n")
        f.write("*END STEP\n")

    return inp_path


# ---------------------------------------------------------------------------
# FRD parser  (identical to triangle script)
# ---------------------------------------------------------------------------

def _parse_frd(path):
    nodes, disp, stress = {}, {}, {}
    state = None

    def _fields(line):
        vals, pos = [], 13
        raw = line.rstrip()
        while pos + 12 <= len(raw):
            try:
                vals.append(float(raw[pos:pos + 12]))
            except ValueError:
                break
            pos += 12
        if not vals:
            for tok in raw[13:].strip().split():
                try:
                    vals.append(float(tok))
                except ValueError:
                    pass
        return vals

    with open(path, "r", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip("\n\r")
            if not line:
                continue
            key6 = line[:6]
            if key6 == "    2C":
                state = "nodes"; continue
            if key6 in ("    3C", "    2P"):
                state = "elements"; continue
            if key6.startswith("  100C"):
                state = "pending"; continue
            if line[:3] == " -4":
                if state == "pending":
                    c = line[3:].upper().strip()
                    if any(k in c for k in ("DISP", "D1 ", "D2 ", "D3 ", " U ")):
                        state = "disp"
                    elif any(k in c for k in ("STRESS", "SXX", "SYY", "SZZ")):
                        state = "stress"
                    else:
                        state = "skip"
                continue
            if line[:3] == " -3" or line.startswith("  -3"):
                state = None; continue
            if line.strip() == "9999":
                break
            if line[:3] == " -1" and state not in (
                    None, "skip", "pending", "elements"):
                try:
                    nid = int(line[3:13])
                except ValueError:
                    continue
                vals = _fields(line)
                if state == "nodes" and len(vals) >= 3:
                    nodes[nid] = (vals[0], vals[1], vals[2])
                elif state == "disp" and len(vals) >= 3:
                    disp[nid] = (vals[0], vals[1], vals[2])
                elif state == "stress" and len(vals) >= 6:
                    stress[nid] = tuple(vals[:6])

    return nodes, disp, stress


# ---------------------------------------------------------------------------
# Run one tension step
# ---------------------------------------------------------------------------

def run_one_step(mesh_nodes, mesh_elements, T_Nm,
                 wall_nodes, pillar_nodes, ccx_exe, log):
    T_applied = tension_to_applied_T(T_Nm)
    delta_T   = T_REF - T_applied
    info(log, "  T={} N/m  →  ΔT={:.2f} K  →  T_applied={:.2f} K".format(
        T_Nm, delta_T, T_applied))

    work_dir = tempfile.mkdtemp(prefix="bask_T{}_".format(T_Nm))
    inp_path = write_inp(work_dir, mesh_nodes, mesh_elements,
                         T_Nm, wall_nodes, pillar_nodes)
    inp_stem = os.path.join(work_dir, "basket")
    info(log, "  Running CalculiX: " + inp_path)

    result = subprocess.run(
        [ccx_exe, "-i", inp_stem],
        cwd=work_dir,
        capture_output=True,
        text=True,
        timeout=3600,
    )
    if result.returncode != 0:
        info(log, "  CCX stdout (tail): " + result.stdout[-600:])
        info(log, "  CCX stderr: " + result.stderr[-200:])
        raise RuntimeError("CalculiX exited with code {}".format(result.returncode))

    frd_path = inp_stem + ".frd"
    if not os.path.isfile(frd_path):
        raise RuntimeError("FRD not found: " + frd_path)

    info(log, "  Parsing FRD: " + frd_path)
    nodes_frd, disp_frd, stress_frd = _parse_frd(frd_path)
    info(log, "  FRD: {} nodes, {} disp, {} stress".format(
        len(nodes_frd), len(disp_frd), len(stress_frd)))

    if not disp_frd:
        raise RuntimeError("No displacement data in FRD: " + frd_path)

    nids, xs, ys, zs, uxs, uys, uzs, vms = [], [], [], [], [], [], [], []
    for nid in sorted(set(nodes_frd) & set(disp_frd)):
        x, y, z     = nodes_frd[nid]
        ux, uy, uz  = disp_frd[nid]
        nids.append(nid)
        xs.append(round(x, 4));  ys.append(round(y, 4));  zs.append(round(z, 4))
        uxs.append(round(ux, 6)); uys.append(round(uy, 6)); uzs.append(round(uz, 6))
        if nid in stress_frd:
            s = stress_frd[nid]
            vm2 = (s[0]**2 + s[1]**2 + s[2]**2
                   - s[0]*s[1] - s[1]*s[2] - s[2]*s[0]
                   + 3.0*(s[3]**2 + s[4]**2 + s[5]**2))
            vms.append(round(math.sqrt(max(vm2, 0.0)), 4))
        else:
            vms.append(0.0)

    info(log, "  Extracted {} result nodes".format(len(nids)))
    return {"node_ids": nids, "x": xs, "y": ys, "z": zs,
            "Ux": uxs, "Uy": uys, "Uz": uzs, "von_mises_MPa": vms}


# ---------------------------------------------------------------------------
# CSV writers
# ---------------------------------------------------------------------------

def _save_step_csv(T_Nm, data):
    os.makedirs(RESULTS_DIR, exist_ok=True)
    fname = os.path.join(RESULTS_DIR, "basket_gmsh_T{:04d}Nm.csv".format(T_Nm))
    with open(fname, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["node_id", "x_mm", "y_mm", "z_mm",
                    "Ux_mm", "Uy_mm", "Uz_mm", "vonMises_MPa"])
        for i in range(len(data["node_ids"])):
            w.writerow([data["node_ids"][i],
                        data["x"][i], data["y"][i], data["z"][i],
                        data["Ux"][i], data["Uy"][i], data["Uz"][i],
                        data["von_mises_MPa"][i]])
    return fname


def _save_summary_csv(rows):
    os.makedirs(RESULTS_DIR, exist_ok=True)
    fname = os.path.join(RESULTS_DIR, "basket_gmsh_summary.csv")
    with open(fname, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["T_Nm", "T_Ncm", "delta_T_K",
                    "max_Uz_mm", "min_Uz_mm", "max_abs_Uz_mm",
                    "max_vonMises_MPa"])
        for row in rows:
            w.writerow(row)
    return fname


# ---------------------------------------------------------------------------
# Main sweep
# ---------------------------------------------------------------------------

def run_sweep():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    log = open(LOG_FILE, "w", encoding="utf-8", buffering=1)

    info(log, "=" * 60)
    info(log, "P2 basket — Gmsh re-mesh tension sweep")
    info(log, "Tensions : {} N/m  ({} N/cm)".format(
        TENSION_SWEEP_Nm, [t // 100 for t in TENSION_SWEEP_Nm]))
    info(log, "Gmsh     : CharacteristicLength {}-{} mm, C3D10".format(
        GMSH_LMIN, GMSH_LMAX))
    info(log, "=" * 60)

    try:
        import FreeCAD
        import ObjectsFem
    except ImportError:
        info(log, "ERROR: FreeCAD module not available. "
             "Run via FreeCADCmd.exe")
        log.close()
        return

    try:
        ccx_exe = find_ccx_exe(log)
    except RuntimeError as exc:
        info(log, "ERROR: " + str(exc))
        log.close()
        return

    # --- create or open the Gmsh-meshed basket FCStd ---
    doc = None
    if os.path.exists(GMSH_FCSTD_PATH):
        info(log, "\nOpening existing Gmsh basket FCStd: " + GMSH_FCSTD_PATH)
        doc = FreeCAD.openDocument(GMSH_FCSTD_PATH)
        n = len(next(iter([o for o in doc.Objects if hasattr(o, "FemMesh")]),
                     type("", (), {"FemMesh": type("", (), {"Nodes": {}})()})()
                     ).FemMesh.Nodes)
        if n == 0:
            info(log, "  Empty mesh — regenerating ...")
            FreeCAD.closeDocument(doc.Name)
            os.remove(GMSH_FCSTD_PATH)
            doc = None

    if doc is None:
        info(log, "\nCreating Gmsh basket mesh ...")
        try:
            doc = create_basket_gmsh_fcstd(GMSH_FCSTD_PATH, log)
        except Exception as exc:
            info(log, "ERROR: " + str(exc))
            info(log, traceback.format_exc())
            log.close()
            return

    # --- find the FemMesh object ---
    mesh_obj = None
    for obj in doc.Objects:
        if hasattr(obj, "FemMesh") and len(obj.FemMesh.Nodes) > 0:
            mesh_obj = obj
            break
    if mesh_obj is None:
        info(log, "ERROR: No FemMesh object found in document")
        log.close()
        return
    info(log, "\nMesh object: '{}' ({} nodes)".format(
        mesh_obj.Label, len(mesh_obj.FemMesh.Nodes)))

    # --- extract mesh data ---
    info(log, "\nExtracting mesh data ...")
    try:
        mesh_nodes, mesh_elements = extract_mesh_data(mesh_obj, log)
    except Exception as exc:
        info(log, "ERROR extracting mesh: " + str(exc))
        info(log, traceback.format_exc())
        log.close()
        return

    # --- build BC node sets once (geometry fixed across all tensions) ---
    info(log, "\nBuilding basket BC node sets ...")
    try:
        wall_nodes, pillar_nodes = build_bc_node_sets(mesh_nodes, log)
    except Exception as exc:
        info(log, "ERROR building BCs: " + str(exc))
        info(log, traceback.format_exc())
        log.close()
        return

    # --- tension sweep ---
    summary_rows = []
    failed = []

    for T_Nm in TENSION_SWEEP_Nm:
        info(log, "\n--- T = {} N/m  ({} N/cm) ---".format(T_Nm, T_Nm // 100))
        try:
            data = run_one_step(mesh_nodes, mesh_elements, T_Nm,
                                wall_nodes, pillar_nodes, ccx_exe, log)
            fname = _save_step_csv(T_Nm, data)
            info(log, "  Saved: " + fname)

            uzs = data["Uz"]
            max_uz     = max(uzs)
            min_uz     = min(uzs)
            max_abs_uz = max(abs(v) for v in uzs)
            max_vm     = max(data["von_mises_MPa"])
            info(log, "  max|Uz| = {:.4f} mm  ({:.1f} µm)".format(
                max_abs_uz, max_abs_uz * 1000))
            info(log, "  max vonMises = {:.4f} MPa".format(max_vm))
            summary_rows.append([
                T_Nm, round(T_Nm / 100.0, 1),
                round(T_REF - tension_to_applied_T(T_Nm), 2),
                round(max_uz, 6), round(min_uz, 6),
                round(max_abs_uz, 6), round(max_vm, 4),
            ])
        except Exception as exc:
            info(log, "  FAILED: " + str(exc))
            info(log, traceback.format_exc())
            failed.append(T_Nm)

    if summary_rows:
        fname = _save_summary_csv(summary_rows)
        info(log, "\nSummary saved: " + fname)

    info(log, "\n" + "=" * 60)
    if failed:
        info(log, "FAILED steps: {}".format(failed))
    else:
        info(log, "All {} steps completed.".format(len(TENSION_SWEEP_Nm)))
    info(log, "Results: " + RESULTS_DIR)
    info(log, "=" * 60)
    log.close()


run_sweep()