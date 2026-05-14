"""
freecad_triangle_fem_setup.py  —  Tension sweep for equilateral triangle FEM

Equilateral triangle with the same area as the original P2 basket trapezoid
(~169 256 mm²).  Runs three tensions: 3, 7, 10 N/cm (300, 700, 1000 N/m).

Run inside FreeCAD's headless interpreter:
    "C:\Program Files\FreeCAD 0.20\bin\FreeCADCmd.exe" FreeCAD/freecad_triangle_fem_setup.py

What it does
------------
1. Creates the triangle FreeCAD document + Gmsh mesh (C:\Temp\triangle_fem.FCStd)
   if it does not already exist.
2. Extracts mesh nodes and elements directly from the FemMesh object.
3. For each tension T in TENSION_SWEEP_Nm:
   a. Builds wall BC (1 mm from every triangle edge, Ux=Uy=Uz=0) and pillar BC
      (3 mm pitch inside the triangle active area, Uz=0) node sets using
      triangle-aware perpendicular-distance geometry.
   b. Writes a complete CalculiX .inp file from scratch (thermal pre-tension
      analogy, same constants as freecad_fem_setup.py).
   c. Runs ccx.exe directly as a subprocess.
   d. Parses displacement + von Mises from the .frd file.
   e. Saves results/triangle_simulation/triangle_sweep_T<XXXX>Nm.csv
4. Saves results/triangle_simulation/triangle_sweep_summary.csv

Approach
--------
FreeCAD is used only for geometry creation and Gmsh meshing.  The FEM analysis
objects (Analysis, Solver, etc.) are not needed — CalculiX is invoked directly.
This sidesteps the ObjectsFem API differences across FreeCAD builds.

Thermal analogy (identical constants to freecad_fem_setup.py)
--------------------------------------------------------------
  ε = α_mesh × |ΔT|  =  T / (E_mesh × t_mesh)
  → ΔT = T / (E_mesh × t_mesh × α_mesh)
  E_mesh = 1 GPa,  t_mesh = 0.5 mm,  α_mesh = 50 µm/m/K,  T_ref = 300 K
  → ΔT_per_Nm = 1/25 K·m/N
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
RESULTS_DIR  = os.path.join(SCRIPT_DIR, "results", "triangle_simulation")
FCSTD_PATH   = r"C:\Temp\triangle_fem.FCStd"
LOG_FILE     = os.path.join(SCRIPT_DIR, "triangle_sweep.log")

# Candidates for the CalculiX executable bundled with FreeCAD 0.20
CCX_CANDIDATES = [
    r"C:\Program Files\FreeCAD 0.20\bin\ccx.exe",
    r"C:\Program Files\FreeCAD 0.20\bin\ccx",
    r"C:\Program Files\FreeCAD\bin\ccx.exe",
    r"C:\Program Files\FreeCAD\bin\ccx",
]

# ---------------------------------------------------------------------------
# Tension sweep  (only 3, 7, 10 N/cm)
# ---------------------------------------------------------------------------
TENSION_SWEEP_Nm = [300, 700, 1000]   # N/m

# ---------------------------------------------------------------------------
# Triangle geometry — same area as original P2 basket trapezoid
# ---------------------------------------------------------------------------
AREA_TARGET_MM2  = 169_256
SIDE_MM          = math.sqrt(4 * AREA_TARGET_MM2 / math.sqrt(3))   # 625.2 mm
HEIGHT_MM        = SIDE_MM * math.sqrt(3) / 2                       # 541.4 mm
THICKNESS_MM     = 0.5   # effective plate thickness mm

# Vertices: base at y=0, apex pointing up, bottom-left corner at origin
TRI_VERTS = [
    (0.0,           0.0       ),
    (SIDE_MM,       0.0       ),
    (SIDE_MM / 2.0, HEIGHT_MM ),
]

# ---------------------------------------------------------------------------
# Mesh material constants  (must match freecad_fem_setup.py)
# ---------------------------------------------------------------------------
E_MESH  = 1.0e9    # Pa
T_MESH  = 0.5e-3   # m
ALPHA   = 50.0e-6  # /K
T_REF   = 300.0    # K  (stress-free reference temperature)
POISSON = 0.3

# ---------------------------------------------------------------------------
# BC parameters
# ---------------------------------------------------------------------------
WALL_MM          = 1.0   # mm  border → Ux=Uy=Uz=0
PILLAR_PITCH_MM  = 3.0   # mm  pillar pitch
PILLAR_TOL_MM    = 1.5   # mm  max snap distance pillar → nearest node

# Mesh node z-filter — wide range so the whole 0.5 mm solid is included
MESH_Z_LO_MM     = -0.1
MESH_Z_HI_MM     =  1.0


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
def info(log, msg):
    if log:
        log.write(msg + "\n")
        log.flush()
    try:
        print(msg)
        sys.stdout.flush()
    except Exception:
        pass


# ---------------------------------------------------------------------------
# Triangle geometry helpers  (mirrors triangle_geometry_check.py)
# ---------------------------------------------------------------------------

def _cross_z(ax, ay, bx, by, cx, cy):
    return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)


def point_in_triangle(px, py):
    (ax, ay), (bx, by), (cx, cy) = TRI_VERTS
    d1 = _cross_z(ax, ay, bx, by, px, py)
    d2 = _cross_z(bx, by, cx, cy, px, py)
    d3 = _cross_z(cx, cy, ax, ay, px, py)
    has_neg = (d1 < 0) or (d2 < 0) or (d3 < 0)
    has_pos = (d1 > 0) or (d2 > 0) or (d3 > 0)
    return not (has_neg and has_pos)


def _dist_to_edge(px, py, x1, y1, x2, y2):
    dx, dy = x2 - x1, y2 - y1
    length = math.hypot(dx, dy)
    if length < 1e-12:
        return math.hypot(px - x1, py - y1)
    return abs(dx * (y1 - py) - dy * (x1 - px)) / length


def min_dist_to_boundary(px, py):
    edges = [
        (TRI_VERTS[0], TRI_VERTS[1]),
        (TRI_VERTS[1], TRI_VERTS[2]),
        (TRI_VERTS[2], TRI_VERTS[0]),
    ]
    return min(_dist_to_edge(px, py, a[0], a[1], b[0], b[1]) for a, b in edges)


# ---------------------------------------------------------------------------
# Tension → temperature conversion
# ---------------------------------------------------------------------------

def tension_to_applied_T(T_Nm):
    return T_REF - T_Nm / (E_MESH * T_MESH * ALPHA)


# ---------------------------------------------------------------------------
# CalculiX executable discovery
# ---------------------------------------------------------------------------

def find_ccx_exe(log):
    # Try known install locations first
    for p in CCX_CANDIDATES:
        if os.path.isfile(p):
            info(log, "CalculiX: " + p)
            return p
    # Try FreeCAD's stored preference
    try:
        import FreeCAD
        pref = FreeCAD.ParamGet(
            "User parameter:BaseApp/Preferences/Mod/Fem/Ccx")
        stored = pref.GetString("ccxBinaryPath", "")
        if stored and os.path.isfile(stored):
            info(log, "CalculiX (from FreeCAD prefs): " + stored)
            return stored
    except Exception:
        pass
    raise RuntimeError(
        "CalculiX executable not found.  Expected at:\n  " +
        "\n  ".join(CCX_CANDIDATES))


# ---------------------------------------------------------------------------
# FreeCAD document creation  (solid + mesh only; no FEM analysis objects)
# ---------------------------------------------------------------------------

def create_triangle_fcstd(save_path, log):
    """
    Build a FreeCAD document with triangle solid and Gmsh mesh.
    No FEM analysis objects are created — CalculiX is run directly.
    """
    import FreeCAD
    import Part
    import ObjectsFem
    from femmesh import gmshtools

    info(log, "Creating new triangle FreeCAD document ...")
    doc = FreeCAD.newDocument("TriangleFEM")

    # Triangle solid: extruded face, z = 0 → THICKNESS_MM
    v = [FreeCAD.Vector(x, y, 0.0) for x, y in TRI_VERTS]
    edges = [
        Part.LineSegment(v[0], v[1]).toShape(),
        Part.LineSegment(v[1], v[2]).toShape(),
        Part.LineSegment(v[2], v[0]).toShape(),
    ]
    wire  = Part.Wire(edges)
    face  = Part.Face(wire)
    solid = face.extrude(FreeCAD.Vector(0, 0, THICKNESS_MM))

    solid_obj = doc.addObject("Part::Feature", "TriangleSolid")
    solid_obj.Shape = solid
    doc.recompute()
    info(log, "  Triangle solid: side={:.1f} mm, height={:.1f} mm, "
              "thickness={:.1f} mm".format(SIDE_MM, HEIGHT_MM, THICKNESS_MM))

    # Gmsh mesh — first-order (C3D4) to avoid mid-edge Jacobian issues
    mesh_obj = ObjectsFem.makeMeshGmsh(doc, "FEMMeshGmsh")
    mesh_obj.Part = solid_obj
    mesh_obj.CharacteristicLengthMin = 2.0
    mesh_obj.CharacteristicLengthMax = 5.0
    mesh_obj.ElementOrder = "1st"
    doc.recompute()

    info(log, "  Running Gmsh (char. length 2–5 mm) ...")
    gmsh_tool = gmshtools.GmshTools(mesh_obj)
    gmsh_msg  = gmsh_tool.create_mesh()
    if gmsh_msg:
        info(log, "  Gmsh note: {}".format(str(gmsh_msg).strip()[:300]))
    doc.recompute()

    n_nodes = len(mesh_obj.FemMesh.Nodes)
    if n_nodes == 0:
        raise RuntimeError("Gmsh produced an empty mesh")
    info(log, "  Mesh nodes: {}".format(n_nodes))

    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    doc.saveAs(save_path)
    info(log, "  Saved: " + save_path)
    return doc


# ---------------------------------------------------------------------------
# Mesh data extraction
# ---------------------------------------------------------------------------

def _tet_signed_vol6(nodes, n0, n1, n2, n3):
    """6× signed volume of a tet.  Returns 0 for degenerate/inverted element."""
    try:
        x0, y0, z0 = nodes[n0]
        x1, y1, z1 = nodes[n1]
        x2, y2, z2 = nodes[n2]
        x3, y3, z3 = nodes[n3]
    except KeyError:
        return 0.0
    ax, ay, az = x1-x0, y1-y0, z1-z0
    bx, by, bz = x2-x0, y2-y0, z2-z0
    cx, cy, cz = x3-x0, y3-y0, z3-z0
    return (ax*(by*cz - bz*cy) +
            ay*(bz*cx - bx*cz) +
            az*(bx*cy - by*cx))


def extract_mesh_data(mesh_obj, log):
    """
    Extract nodes {nid: (x,y,z)} and volume elements {eid: (n1,...)}
    from a FreeCAD FemMesh object.  Elements with near-zero or inverted
    Jacobians (ill-shaped tets from Gmsh) are filtered out so CalculiX
    does not abort.
    """
    fem_mesh = mesh_obj.FemMesh

    nodes = {}
    for nid, v in fem_mesh.Nodes.items():
        nodes[nid] = (v.x, v.y, v.z)

    # Compute signed volumes for all elements first, then filter by majority sign.
    # FreeCAD/Gmsh and CalculiX may use opposite winding conventions; the
    # majority sign is the correct one, and elements with the opposite sign
    # (plus near-zero degenerate ones) are the ill-shaped outliers to remove.
    raw = {}
    for eid in fem_mesh.Volumes:
        enodes = fem_mesh.getElementNodes(eid)
        raw[eid] = (enodes,
                    _tet_signed_vol6(nodes, enodes[0], enodes[1], enodes[2], enodes[3]))

    pos = sum(1 for _, v in raw.values() if v > 1e-10)
    neg = sum(1 for _, v in raw.values() if v < -1e-10)
    expected_sign = 1 if pos >= neg else -1
    info(log, "  Element vol sign: {} positive, {} negative — "
              "majority sign = {}".format(pos, neg, expected_sign))

    elements = {}
    skipped = 0
    for eid, (enodes, vol6) in raw.items():
        if vol6 * expected_sign < 1e-10:   # wrong sign or degenerate
            skipped += 1
            continue
        elements[eid] = enodes

    if skipped:
        info(log, "  Filtered {} ill-shaped element(s)".format(skipped))

    info(log, "  Mesh nodes: {}  Volume elements: {}".format(
        len(nodes), len(elements)))

    if not elements:
        raise RuntimeError(
            "No volume elements found in mesh — "
            "check Gmsh ElementDimension (expected 3)")

    # Report element type from node count of first element
    n_vert = len(next(iter(elements.values())))
    etype_map = {4: "C3D4", 10: "C3D10", 8: "C3D8", 20: "C3D20R"}
    info(log, "  Element type: {} ({}-node)".format(
        etype_map.get(n_vert, "unknown"), n_vert))

    return nodes, elements


# ---------------------------------------------------------------------------
# BC node-set builder  (triangle-aware geometry)
# ---------------------------------------------------------------------------

def build_bc_node_sets(mesh_nodes, log):
    """
    Wall nodes   — within WALL_MM of any triangle edge    → Ux=Uy=Uz=0
    Pillar nodes — nearest interior node per grid centre   → Uz=0

    Uses perpendicular-distance-to-edge geometry (not a bounding box).
    """
    # Restrict to the mesh-layer z range
    surf = {nid: (x, y, z) for nid, (x, y, z) in mesh_nodes.items()
            if MESH_Z_LO_MM <= z <= MESH_Z_HI_MM}
    if not surf:
        info(log, "  WARNING: z filter empty — using all nodes")
        surf = dict(mesh_nodes)
    info(log, "  Nodes after z-filter: {}".format(len(surf)))

    # Wall nodes: within WALL_MM of any edge
    wall_nodes = {nid for nid, (x, y, _) in surf.items()
                  if min_dist_to_boundary(x, y) <= WALL_MM}
    info(log, "  Wall BC nodes: {}".format(len(wall_nodes)))

    xs = [c[0] for c in surf.values()]
    ys = [c[1] for c in surf.values()]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)

    # Spatial hash of interior nodes
    cell = PILLAR_PITCH_MM
    grid = {}
    for nid, (x, y, _) in surf.items():
        if nid not in wall_nodes:
            key = (int(x / cell), int(y / cell))
            grid.setdefault(key, []).append((nid, x, y))

    pillar_nodes = set()
    tol2 = PILLAR_TOL_MM ** 2
    pcx = xmin + WALL_MM + PILLAR_PITCH_MM / 2.0
    while pcx <= xmax - WALL_MM:
        pcy = ymin + WALL_MM + PILLAR_PITCH_MM / 2.0
        while pcy <= ymax - WALL_MM:
            if (point_in_triangle(pcx, pcy) and
                    min_dist_to_boundary(pcx, pcy) > WALL_MM):
                cx_cell = int(pcx / cell)
                cy_cell = int(pcy / cell)
                best_nid, best_d2 = None, tol2
                for di in (-1, 0, 1):
                    for dj in (-1, 0, 1):
                        for nid, nx, ny in grid.get(
                                (cx_cell + di, cy_cell + dj), []):
                            d2 = (nx - pcx) ** 2 + (ny - pcy) ** 2
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

def write_inp(work_dir, mesh_nodes, mesh_elements,
              T_Nm, wall_nodes, pillar_nodes):
    """Write a complete CalculiX .inp for one tension step."""
    T_applied = tension_to_applied_T(T_Nm)

    # Determine CalculiX element type from vertex count
    n_vert = len(next(iter(mesh_elements.values())))
    etype  = {4: "C3D4", 10: "C3D10", 8: "C3D8", 20: "C3D20R"}.get(
        n_vert, "C3D4")

    inp_path = os.path.join(work_dir, "triangle.inp")
    with open(inp_path, "w", encoding="utf-8") as f:
        f.write("** Triangle FEM — T={} N/m ({:.0f} N/cm)\n".format(
            T_Nm, T_Nm / 100.0))
        f.write("** freecad_triangle_fem_setup.py\n\n")

        # Nodes
        f.write("*NODE, NSET=NALL\n")
        for nid, (x, y, z) in sorted(mesh_nodes.items()):
            f.write("{}, {:.6f}, {:.6f}, {:.6f}\n".format(nid, x, y, z))

        # Elements
        # FreeCAD/Gmsh returns tet nodes in left-handed order; CalculiX C3D4
        # expects right-handed.  Swapping nodes at positions 1 and 2 reverses
        # the winding (same fix that importCcxInpMesh.py applies internally).
        f.write("*ELEMENT, TYPE={}, ELSET=EALL\n".format(etype))
        for eid, enodes in sorted(mesh_elements.items()):
            if len(enodes) == 4:
                out = (enodes[0], enodes[2], enodes[1], enodes[3])
            else:
                out = enodes
            f.write("{}, {}\n".format(eid, ", ".join(str(n) for n in out)))

        # Material (E in MPa for CalculiX when geometry is in mm)
        f.write("*MATERIAL, NAME=TRIANGLEMESH\n")
        f.write("*ELASTIC\n")
        f.write("{:.1f}, {:.2f}\n".format(E_MESH / 1e6, POISSON))
        f.write("*EXPANSION, ZERO={:.1f}\n".format(T_REF))
        f.write("{:.6E}\n".format(ALPHA))

        f.write("*SOLID SECTION, ELSET=EALL, MATERIAL=TRIANGLEMESH\n\n")

        # Wall BCs
        f.write("** Wall BC (1 mm from each edge): Ux=Uy=Uz=0\n")
        f.write("*BOUNDARY\n")
        for nid in sorted(wall_nodes):
            f.write("{}, 1, 3, 0.0\n".format(nid))

        # Pillar BCs
        if pillar_nodes:
            f.write("** Pillar BC ({} mm pitch): Uz=0\n".format(
                int(PILLAR_PITCH_MM)))
            f.write("*BOUNDARY\n")
            for nid in sorted(pillar_nodes):
                f.write("{}, 3, 3, 0.0\n".format(nid))

        # Initial temperature (stress-free reference)
        f.write("*INITIAL CONDITIONS, TYPE=TEMPERATURE\n")
        f.write("NALL, {:.1f}\n\n".format(T_REF))

        # Static step with thermal load
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
# FRD parser  (identical to freecad_fem_setup.py)
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
    """Write .inp, run ccx, parse .frd, return result dict."""
    T_applied = tension_to_applied_T(T_Nm)
    delta_T   = T_REF - T_applied
    info(log, "  T={} N/m  →  ΔT={:.2f} K  →  T_applied={:.2f} K".format(
        T_Nm, delta_T, T_applied))

    work_dir = tempfile.mkdtemp(prefix="tri_T{}_".format(T_Nm))
    inp_path = write_inp(work_dir, mesh_nodes, mesh_elements,
                         T_Nm, wall_nodes, pillar_nodes)

    inp_stem = os.path.join(work_dir, "triangle")
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
        raise RuntimeError(
            "CalculiX exited with code {}".format(result.returncode))

    frd_path = inp_stem + ".frd"
    if not os.path.isfile(frd_path):
        raise RuntimeError("FRD not found: " + frd_path)

    info(log, "  Parsing FRD: " + frd_path)
    nodes_frd, disp_frd, stress_frd = _parse_frd(frd_path)
    info(log, "  FRD: {} nodes, {} disp, {} stress".format(
        len(nodes_frd), len(disp_frd), len(stress_frd)))

    if not disp_frd:
        raise RuntimeError("No displacement data in FRD: " + frd_path)

    nids, xs, ys, zs = [], [], [], []
    uxs, uys, uzs, vms = [], [], [], []

    for nid in sorted(set(nodes_frd) & set(disp_frd)):
        x, y, z    = nodes_frd[nid]
        ux, uy, uz = disp_frd[nid]
        nids.append(nid)
        xs.append(round(x,  4));  ys.append(round(y,  4));  zs.append(round(z,  4))
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
    fname = os.path.join(RESULTS_DIR,
                         "triangle_sweep_T{:04d}Nm.csv".format(T_Nm))
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
    fname = os.path.join(RESULTS_DIR, "triangle_sweep_summary.csv")
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
    info(log, "Triangle FEM tension sweep")
    info(log, "Geometry : equilateral triangle, side={:.1f} mm, "
              "area={:,.0f} mm²".format(SIDE_MM, AREA_TARGET_MM2))
    info(log, "Tensions : {} N/m  ({} N/cm)".format(
        TENSION_SWEEP_Nm, [t // 100 for t in TENSION_SWEEP_Nm]))
    info(log, "=" * 60)

    try:
        import FreeCAD
    except ImportError:
        info(log, "ERROR: FreeCAD module not available. "
             "Run this script via FreeCADCmd.exe")
        log.close()
        return

    # Find CalculiX executable
    try:
        ccx_exe = find_ccx_exe(log)
    except RuntimeError as exc:
        info(log, "ERROR: " + str(exc))
        log.close()
        return

    # Open or create the triangle document.
    # If the saved FCStd has a C3D10 mesh (quadratic — ill-shaped Jacobian risk),
    # delete it and regenerate with C3D4 (first-order, no mid-edge nodes).
    def _mesh_order(doc_obj):
        for obj in doc_obj.Objects:
            if hasattr(obj, "FemMesh"):
                vols = list(obj.FemMesh.Volumes)
                if vols:
                    return len(obj.FemMesh.getElementNodes(vols[0]))
        return 0

    if os.path.exists(FCSTD_PATH):
        info(log, "\nOpening: " + FCSTD_PATH)
        doc = FreeCAD.openDocument(FCSTD_PATH)
        order = _mesh_order(doc)
        if order > 4:
            info(log, "  Mesh is {}-node (C3D10) — regenerating as C3D4 ...".format(order))
            FreeCAD.closeDocument(doc.Name)
            try:
                os.remove(FCSTD_PATH)
            except OSError as exc:
                info(log, "ERROR: cannot delete old FCStd: " + str(exc))
                log.close()
                return
            doc = None
        else:
            info(log, "  Mesh is C3D4 — OK")
    else:
        doc = None

    if doc is None:
        info(log, "\nCreating triangle document from scratch ...")
        try:
            doc = create_triangle_fcstd(FCSTD_PATH, log)
        except Exception as exc:
            info(log, "ERROR: could not create FCStd: " + str(exc))
            info(log, traceback.format_exc())
            log.close()
            return

    # Locate the FEM mesh object
    mesh_obj = None
    for obj in doc.Objects:
        if hasattr(obj, "FemMesh") and "Result" not in obj.Name:
            mesh_obj = obj
            info(log, "\nMesh object: '{}' ({} nodes)".format(
                obj.Name, len(obj.FemMesh.Nodes)))
            break
    if mesh_obj is None:
        info(log, "ERROR: no FEM mesh found in document")
        log.close()
        return

    # Extract mesh data once — geometry is constant across all tension steps
    info(log, "\nExtracting mesh data ...")
    try:
        mesh_nodes, mesh_elements = extract_mesh_data(mesh_obj, log)
    except Exception as exc:
        info(log, "ERROR: " + str(exc))
        info(log, traceback.format_exc())
        log.close()
        return

    # Build BC node sets once
    info(log, "\nBuilding triangle BC node sets ...")
    try:
        wall_nids, pillar_nids = build_bc_node_sets(mesh_nodes, log)
    except Exception as exc:
        info(log, "ERROR: BC node-set computation failed: " + str(exc))
        info(log, traceback.format_exc())
        log.close()
        return

    summary_rows = []
    failed = []

    for T_Nm in TENSION_SWEEP_Nm:
        info(log, "\n--- T = {} N/m  ({:.0f} N/cm) ---".format(
            T_Nm, T_Nm / 100.0))
        try:
            data     = run_one_step(mesh_nodes, mesh_elements, T_Nm,
                                    wall_nids, pillar_nids, ccx_exe, log)
            csv_path = _save_step_csv(T_Nm, data)
            info(log, "  Saved: " + csv_path)

            uzs     = data["Uz"]
            max_uz  = max(uzs)
            min_uz  = min(uzs)
            max_abs = max(abs(v) for v in uzs)
            max_vm  = max(data["von_mises_MPa"])
            delta_T = T_REF - tension_to_applied_T(T_Nm)

            info(log, "  max|Uz| = {:.4f} mm  ({:.1f} µm)".format(
                max_abs, max_abs * 1000))
            info(log, "  max vonMises = {:.4f} MPa".format(max_vm))

            summary_rows.append([
                T_Nm,
                round(T_Nm / 100.0, 1),
                round(delta_T, 2),
                round(max_uz, 6),
                round(min_uz, 6),
                round(max_abs, 6),
                round(max_vm, 4),
            ])

        except Exception as exc:
            info(log, "  FAILED: " + str(exc))
            info(log, traceback.format_exc())
            failed.append(T_Nm)

    if summary_rows:
        spath = _save_summary_csv(summary_rows)
        info(log, "\nSummary saved: " + spath)

    info(log, "\n" + "=" * 60)
    if failed:
        info(log, "FAILED steps: {}".format(failed))
    else:
        info(log, "All {} steps completed.".format(len(TENSION_SWEEP_Nm)))
    info(log, "Results: " + RESULTS_DIR)
    info(log, "=" * 60)
    log.close()


run_sweep()
