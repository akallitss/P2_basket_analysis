"""
freecad_baking_buckling.py  -  Baking buckling eigenvalue sweep
P2 Basket Bulk Micromegas mesh

Physics
-------
During baking (25 C -> 146 C, DT_real = 121 K):
  SS mesh is bonded to Dynamask photoresist at both the wall frame and the
  pillar tops.  Dynamask is chemically anchored to the FR4 PCB (same UV step
  that creates the pillars).  As FR4 expands more slowly than SS
  (ALPHA_FR4 < ALPHA_SS), the constraints pull inward relative to the mesh:
  net in-plane compression builds up, pre-tension decreases, and the mesh
  buckles if compression exceeds the initial tension.

Equivalent model baking temperature (thermal analogy, alpha_eff = 50e-6/K):
  DT_bake_model = DT_bake_real * (delta_alpha / alpha_eff)
                = 121 * (2e-6 / 50e-6) = 4.84 K   (heating reference load)

Step 1  *STATIC      : establish pre-tension via thermal cooling
Step 2  *BUCKLE      : reference load = DT_bake_model K heating, all BCs fixed

Eigenvalue interpretation:
  lambda < 1  ->  buckles at ~(25 + lambda*121) C  (before end of bake)
  lambda = 1  ->  buckles exactly at 146 C
  lambda > 1  ->  survives full bake

Run:
    "C:\\Program Files\\FreeCAD 0.20\\bin\\FreeCADCmd.exe" FreeCAD/freecad_baking_buckling.py

Output:
    FreeCAD/results/baking_buckle_summary.csv
    FreeCAD/results/baking_buckle_T????Nm_mode??.csv
    FreeCAD/baking_buckle.log
"""

import os
import sys
import csv
import math
import traceback

# ---------------------------------------------------------------------------
# CTE parameters
# Update ALPHA_FR4 if Shengyi publishes the S1000-2M in-plane CTE.
# ---------------------------------------------------------------------------
ALPHA_FR4   = 14e-6   # /K  in-plane CTE of FR4 / S1000-2M PCB  (ESTIMATE)
ALPHA_SS    = 16e-6   # /K  SS316L mesh
DELTA_ALPHA = ALPHA_SS - ALPHA_FR4   # 2e-6 /K  driving differential

# ---------------------------------------------------------------------------
# Thermal analogy constants  (must match freecad_fem_setup.py)
# ---------------------------------------------------------------------------
E_EFF     = 1.0e9    # Pa    effective Young's modulus
T_EFF     = 0.5e-3   # m     effective thickness
ALPHA_EFF = 50.0e-6  # /K    effective CTE (artificial, pre-tension analogy only)
T_REF     = 300.0    # K     initial / reference temperature

# ---------------------------------------------------------------------------
# Baking parameters
# ---------------------------------------------------------------------------
T_ROOM_C      = 25.0    # C  ambient
T_BAKE_C      = 146.0   # C  baking temperature
DT_BAKE_REAL  = T_BAKE_C - T_ROOM_C   # 121 K  real baking delta-T

# Model equivalent: same thermal strain as delta_alpha * DT_BAKE_REAL
# but expressed with alpha_eff so CalculiX material needs no change
DT_BAKE_MODEL = DT_BAKE_REAL * (DELTA_ALPHA / ALPHA_EFF)  # 4.84 K

assert DELTA_ALPHA > 0, "ALPHA_SS must be > ALPHA_FR4 for baking to compress the mesh"
assert DT_BAKE_MODEL > 0

# ---------------------------------------------------------------------------
# Sweep
# ---------------------------------------------------------------------------
TENSIONS_Nm  = [300, 500, 700, 1000]   # N/m  pre-tension values
N_EIGENMODES = 5                        # eigenmodes extracted per run

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(SCRIPT_DIR, "results")
FCSTD_PATH  = r"C:\Temp\fc_v4_open.FCStd"
LOG_FILE    = os.path.join(SCRIPT_DIR, "baking_buckle.log")

# ---------------------------------------------------------------------------
# FreeCAD object names
# ---------------------------------------------------------------------------
OBJ_MESH   = "FEMMeshGmsh"
OBJ_FEM    = "Analysis"
OBJ_SOLVER = "SolverCcxTools"

# ---------------------------------------------------------------------------
# BC geometry
# ---------------------------------------------------------------------------
MESH_Z_LO_MM    = 0.05   # mm  lower z bound for mesh-layer nodes
MESH_Z_HI_MM    = 0.70   # mm  upper z bound
WALL_MM         = 1.0    # mm  peripheral wall border thickness
PILLAR_PITCH_MM = 3.0    # mm  Dynamask pillar centre-to-centre pitch
PILLAR_SNAP_MM  = 1.5    # mm  snap radius: nodes within this -> pillar BC


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def info(log, msg):
    print(msg)
    if log:
        log.write(msg + "\n")
        log.flush()


def tension_to_applied_T(T_Nm):
    delta_T = T_Nm / (E_EFF * T_EFF * ALPHA_EFF)
    return T_REF - delta_T


# ---------------------------------------------------------------------------
# BC node-set builder
# ---------------------------------------------------------------------------

def build_bc_node_sets(fem_mesh_obj, log):
    """
    Returns (wall_nodes, pillar_nodes) — both fully fixed (Ux=Uy=Uz=0).
    fem_mesh_obj : FreeCAD document object with a .FemMesh attribute.
    """
    all_nodes = fem_mesh_obj.FemMesh.Nodes   # dict: nid -> FreeCAD.Vector (mm)

    surf = {nid: (v.x, v.y, v.z) for nid, v in all_nodes.items()
            if MESH_Z_LO_MM <= v.z <= MESH_Z_HI_MM}
    if not surf:
        info(log, "  WARNING: no nodes in z=[{}, {}] mm — using all nodes".format(
            MESH_Z_LO_MM, MESH_Z_HI_MM))
        surf = {nid: (v.x, v.y, v.z) for nid, v in all_nodes.items()}
    info(log, "  Surface nodes: {}".format(len(surf)))

    xs = [c[0] for c in surf.values()]
    ys = [c[1] for c in surf.values()]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)

    # Wall nodes: 1 mm border, fully fixed
    wall_nodes = {nid for nid, (x, y, _) in surf.items()
                  if (x < xmin + WALL_MM or x > xmax - WALL_MM or
                      y < ymin + WALL_MM or y > ymax - WALL_MM)}
    info(log, "  Wall BC nodes:   {}".format(len(wall_nodes)))

    # Pillar nodes: all nodes within PILLAR_SNAP_MM of any pillar centre
    pillar_nodes = set()
    pcx = xmin + WALL_MM + PILLAR_PITCH_MM / 2.0
    while pcx <= xmax - WALL_MM:
        pcy = ymin + WALL_MM + PILLAR_PITCH_MM / 2.0
        while pcy <= ymax - WALL_MM:
            for nid, (x, y, _) in surf.items():
                if math.hypot(x - pcx, y - pcy) <= PILLAR_SNAP_MM:
                    pillar_nodes.add(nid)
            pcy += PILLAR_PITCH_MM
        pcx += PILLAR_PITCH_MM

    pillar_nodes -= wall_nodes
    info(log, "  Pillar BC nodes: {}".format(len(pillar_nodes)))
    return wall_nodes, pillar_nodes, surf


# ---------------------------------------------------------------------------
# CalculiX .inp injection
# ---------------------------------------------------------------------------

def _find_inp_file(fea):
    for attr in ("inp_file_name", "ccx_inp_file"):
        val = getattr(fea, attr, None)
        if val and isinstance(val, str) and val.endswith(".inp") and os.path.isfile(val):
            return val
    wdir = getattr(fea, "working_dir", None)
    if wdir and os.path.isdir(wdir):
        for root, _, files in os.walk(wdir):
            for fname in files:
                if fname.endswith(".inp"):
                    return os.path.join(root, fname)
    raise RuntimeError("Cannot locate .inp file — working_dir={!r}".format(
        getattr(fea, "working_dir", None)))


def inject_baking_buckle_inp(fea, wall_nodes, pillar_nodes, T_applied, surf, log):
    """
    Replace FreeCAD's *STEP sections with:
      Step 1 *STATIC  : pre-tension cooling to T_applied
      Step 2 *BUCKLE  : reference load = inward *CLOAD at wall nodes, calibrated
                        so that lambda=1 corresponds to the full bake (121 K).
    Wall + pillar BCs (both Ux=Uy=Uz=0) are injected before Step 1.

    *TEMPERATURE is NOT allowed in a CalculiX *BUCKLE step.  The equivalent
    in-plane compressive force is used instead:
      q = delta_alpha * DT_BAKE_REAL * E_EFF * T_EFF  [N/m per unit length]
    Distributed equally over all wall nodes in the inward radial direction.
    """
    inp_path = _find_inp_file(fea)

    with open(inp_path, "r", encoding="utf-8", errors="replace") as f:
        content = f.read()

    content_lower = content.lower()
    idx = content_lower.find("\n*step")
    if idx == -1:
        idx = content_lower.find("*step")
    preamble = content[:idx].rstrip() if idx != -1 else content.rstrip()

    wall_list   = sorted(wall_nodes)
    pillar_list = sorted(pillar_nodes)

    # Calibrate inward CLOAD at wall nodes so lambda=1 = full bake compression.
    # q = line force [N/m] equivalent to full bake differential expansion
    q = DELTA_ALPHA * DT_BAKE_REAL * E_EFF * T_EFF   # N/m
    wall_x = [surf[n][0] for n in wall_nodes if n in surf]
    wall_y = [surf[n][1] for n in wall_nodes if n in surf]
    xmin_w, xmax_w = min(wall_x), max(wall_x)
    ymin_w, ymax_w = min(wall_y), max(wall_y)
    cx = (xmin_w + xmax_w) / 2.0
    cy = (ymin_w + ymax_w) / 2.0
    perim_m = 2.0 * ((xmax_w - xmin_w) + (ymax_w - ymin_w)) * 1e-3   # mm -> m
    F_per_node = q * perim_m / max(len(wall_nodes), 1)   # N per wall node

    lines = [preamble, ""]

    # BCs before first *STEP — persist through all steps
    lines += [
        "**",
        "** === Injected BCs: wall frame + Dynamask pillar supports ===",
        "** Wall frame — Ux=Uy=Uz=0  ({} nodes)".format(len(wall_list)),
        "*BOUNDARY",
    ]
    for nid in wall_list:
        lines.append("{}, 1, 3, 0.0".format(nid))
    lines += [
        "** Dynamask pillars — Ux=Uy=Uz=0  ({} nodes, {} mm pitch)".format(
            len(pillar_list), int(PILLAR_PITCH_MM)),
        "*BOUNDARY",
    ]
    for nid in pillar_list:
        lines.append("{}, 1, 3, 0.0".format(nid))
    lines += ["** === end injected BCs ===", ""]

    # Step 1: static pre-tension
    lines += [
        "*STEP",
        "*STATIC",
        "*TEMPERATURE",
        "NALL, {:.6f}".format(T_applied),
        "*NODE FILE",
        "U",
        "*END STEP",
        "",
    ]

    # Step 2: perturbation buckling — inward CLOAD at wall nodes
    # lambda=1 -> full bake compression (Δα*121K*E*t per unit length)
    lines += [
        "*STEP, PERTURBATION",
        "*BUCKLE",
        "{}".format(N_EIGENMODES),
        "*CLOAD",
    ]
    for nid in wall_list:
        if nid not in surf:
            continue
        x, y = surf[nid][0], surf[nid][1]
        dx, dy = cx - x, cy - y
        dist = math.sqrt(dx * dx + dy * dy)
        if dist < 1e-9:
            continue
        fx = F_per_node * dx / dist
        fy = F_per_node * dy / dist
        lines.append("{}, 1, {:.6E}".format(nid, fx))
        lines.append("{}, 2, {:.6E}".format(nid, fy))
    lines += [
        "*NODE FILE",
        "U",
        "*END STEP",
        "",
    ]

    with open(inp_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    info(log, "  .inp patched: {} wall + {} pillar nodes  T_applied={:.2f} K  "
         "q={:.2f} N/m  F_per_node={:.4f} N  (lambda=1 = full bake)".format(
             len(wall_list), len(pillar_list), T_applied, q, F_per_node))


def patch_write_for_baking_buckle(fea, wall_nodes, pillar_nodes, T_applied, surf, log):
    """
    Monkey-patch fea.write_inp_file so the injection survives the internal
    write_inp_file() call inside fea.run().
    """
    orig = fea.write_inp_file

    def _write_and_inject():
        result = orig()
        try:
            inject_baking_buckle_inp(fea, wall_nodes, pillar_nodes, T_applied, surf, log)
        except Exception as exc:
            info(log, "  WARNING: .inp injection failed — {}".format(exc))
            info(log, traceback.format_exc())
        return result

    fea.write_inp_file = _write_and_inject


# ---------------------------------------------------------------------------
# FRD / DAT parsers
# ---------------------------------------------------------------------------

def _parse_buckle_frd(path):
    """
    Parse all displacement result blocks from a CalculiX .frd file.
    Returns: (nodes dict, list of mode dicts)
      nodes : {nid: (x, y, z)} mm
      modes : [{'mode': int, 'disp': {nid: (ux, uy, uz)}}]
    Includes Step 1 (static) block as well as Step 2 eigenvectors; the
    caller should trim static blocks before pairing with eigenvalues.
    """
    nodes = {}
    modes = []
    state = None
    current_mode = None
    current_disp = {}

    def _fields(line, start=13):
        vals, pos = [], start
        raw = line.rstrip()
        while pos + 12 <= len(raw):
            try:
                vals.append(float(raw[pos:pos + 12]))
            except ValueError:
                break
            pos += 12
        if not vals:
            for tok in raw[start:].strip().split():
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
                if current_mode is not None and current_disp:
                    modes.append({"mode": current_mode, "disp": current_disp})
                current_mode = None
                current_disp = {}
                state = "pending"; continue
            if line[:3] == " -4":
                if state == "pending":
                    tok = line[3:].upper().strip()
                    if any(k in tok for k in ("DISP", "D1 ", "D2 ", "D3 ", " U ")):
                        current_mode = len(modes) + 1
                        state = "disp"
                    else:
                        state = "skip"
                continue
            if line[:3] == " -3" or line.startswith("  -3"):
                state = None; continue
            if line.strip() == "9999":
                break
            if line[:3] == " -1" and state not in (None, "skip", "pending", "elements"):
                try:
                    nid = int(line[3:13])
                except ValueError:
                    continue
                vals = _fields(line)
                if state == "nodes" and len(vals) >= 3:
                    nodes[nid] = (vals[0], vals[1], vals[2])
                elif state == "disp" and len(vals) >= 3:
                    current_disp[nid] = (vals[0], vals[1], vals[2])

    if current_mode is not None and current_disp:
        modes.append({"mode": current_mode, "disp": current_disp})

    return nodes, modes


def _parse_buckle_dat(dat_path):
    """
    Parse buckling eigenvalues from CalculiX .dat file.
    Returns list of (mode_number, eigenvalue) tuples.

    CalculiX 2.17 writes the header with spaced letters:
        B U C K L I N G   F A C T O R S
    so we strip spaces before matching.
    """
    eigenvalues = []
    if not os.path.isfile(dat_path):
        return eigenvalues
    in_buckle = False
    header_lines = []
    with open(dat_path, "r", errors="replace") as f:
        for line in f:
            compact = line.upper().replace(" ", "")
            header_lines.append(line.rstrip())
            if len(header_lines) > 60:
                header_lines.pop(0)
            if "BUCKLINGFACTOR" in compact or ("BUCKL" in compact and "FACTOR" in compact):
                in_buckle = True
                continue
            if in_buckle:
                stripped = line.strip()
                if not stripped:
                    continue
                parts = stripped.split()
                if len(parts) >= 2:
                    try:
                        mode = int(parts[0])
                        val  = float(parts[1])
                        if math.isnan(val) or math.isinf(val):
                            print("  DAT WARNING: mode {} eigenvalue is {} — skipped".format(mode, val))
                            continue
                        eigenvalues.append((mode, val))
                    except ValueError:
                        continue  # skip sub-headers (MODE NO / FACTOR)
    if not eigenvalues:
        print("  DAT DIAGNOSTIC — last {} lines of {}:".format(
            len(header_lines), dat_path))
        for l in header_lines:
            print("    |" + l)
    return eigenvalues


# ---------------------------------------------------------------------------
# Run one tension value
# ---------------------------------------------------------------------------

def run_one_tension(doc, T_Nm, wall_nodes, pillar_nodes, surf, log):
    """
    Run Step 1 (static pretension) + Step 2 (perturbation buckle).
    Returns (nodes_frd, modes, eigenvalues).
    modes[i]['eigenvalue'] is set from the .dat file.
    """
    import femtools.ccxtools as ccx

    T_applied = tension_to_applied_T(T_Nm)
    delta_T   = T_REF - T_applied
    info(log, "  T={} N/m  dT={:.2f} K  T_applied={:.2f} K".format(
        T_Nm, delta_T, T_applied))

    fea = ccx.FemToolsCcx(analysis=doc.getObject(OBJ_FEM),
                          solver=doc.getObject(OBJ_SOLVER))
    fea.update_objects()
    fea.setup_working_dir()
    fea.setup_ccx()

    # Monkey-patch so injection survives the write_inp_file() call inside run()
    patch_write_for_baking_buckle(fea, wall_nodes, pillar_nodes, T_applied, surf, log)

    # Prevent FreeCAD from accumulating CCX_Results objects across runs
    fea.load_results = lambda: None

    info(log, "  Running CalculiX ...")
    msg = fea.run()
    if msg and not isinstance(msg, bool):
        raise RuntimeError("CalculiX error:\n" + str(msg))

    # Locate output files
    wdir = getattr(fea, "working_dir", None)
    frd_path = dat_path = None
    if wdir and os.path.isdir(wdir):
        for root, _, files in os.walk(wdir):
            for fname in files:
                fp = os.path.join(root, fname)
                if fname.endswith(".frd"):
                    frd_path = fp
                if fname.endswith(".dat"):
                    dat_path = fp

    if not frd_path:
        raise RuntimeError("Cannot find .frd in: " + str(wdir))

    info(log, "  Parsing FRD: {}".format(frd_path))
    nodes_frd, modes = _parse_buckle_frd(frd_path)
    info(log, "  FRD: {} nodes, {} displacement block(s)".format(
        len(nodes_frd), len(modes)))

    eigenvalues = _parse_buckle_dat(dat_path) if dat_path else []
    info(log, "  Eigenvalues: {}".format(eigenvalues))

    if not eigenvalues:
        raise RuntimeError("No eigenvalues found in .dat — check CalculiX output")
    if all(math.isnan(v) or math.isinf(v) for _, v in eigenvalues):
        raise RuntimeError("All eigenvalues are NaN/inf — CalculiX solver likely failed")

    # The .frd has one static block (Step 1) followed by N eigenvector blocks.
    # Discard any leading static blocks so indices align with eigenvalues.
    if len(modes) > len(eigenvalues):
        modes = modes[len(modes) - len(eigenvalues):]

    for i, (mode_n, ev) in enumerate(eigenvalues):
        if i < len(modes):
            modes[i]["mode"]       = mode_n
            modes[i]["eigenvalue"] = ev

    return nodes_frd, modes, eigenvalues


# ---------------------------------------------------------------------------
# CSV writers
# ---------------------------------------------------------------------------

def save_mode_csv(nodes, mode_data, T_Nm, mode_idx, results_dir):
    fname = "baking_buckle_T{:04d}Nm_mode{:02d}.csv".format(T_Nm, mode_idx + 1)
    path  = os.path.join(results_dir, fname)
    disp  = mode_data.get("disp", {})
    ev    = mode_data.get("eigenvalue", float("nan"))

    t_buckle = (25.0 + ev * DT_BAKE_REAL) if not math.isnan(ev) else float("nan")

    # Normalise mode shape to max|Uz| = 1 for comparability across tensions
    uzs   = [d[2] for d in disp.values()]
    scale = max(abs(v) for v in uzs) if uzs else 1.0
    if scale == 0.0:
        scale = 1.0

    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["node_id", "x_mm", "y_mm", "z_mm",
                    "Ux_norm", "Uy_norm", "Uz_norm",
                    "eigenvalue", "T_buckle_C"])
        for nid, (x, y, z) in nodes.items():
            if nid in disp:
                ux, uy, uz = disp[nid]
                w.writerow([
                    nid,
                    "{:.4f}".format(x), "{:.4f}".format(y), "{:.4f}".format(z),
                    "{:.6f}".format(ux / scale),
                    "{:.6f}".format(uy / scale),
                    "{:.6f}".format(uz / scale),
                    "{:.6f}".format(ev),
                    "{:.1f}".format(t_buckle),
                ])
    return path


def save_summary_csv(summary_rows, results_dir):
    path = os.path.join(results_dir, "baking_buckle_summary.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["T_Nm", "T_Ncm", "mode", "eigenvalue", "T_buckle_C", "survived"])
        for row in summary_rows:
            ev = row["eigenvalue"]
            if math.isnan(ev) or math.isinf(ev):
                continue
            t_buckle = 25.0 + ev * DT_BAKE_REAL
            w.writerow([
                row["T_Nm"],
                "{:.1f}".format(row["T_Nm"] / 100.0),
                row["mode"],
                "{:.6f}".format(ev),
                "{:.1f}".format(t_buckle),
                "yes" if ev >= 1.0 else "no",
            ])
    return path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    import FreeCAD  # available only inside FreeCADCmd

    os.makedirs(RESULTS_DIR, exist_ok=True)

    with open(LOG_FILE, "w", encoding="utf-8", buffering=1) as log:
        info(log, "=" * 60)
        info(log, "P2 Basket FEM — baking buckling sweep")
        info(log, "Tensions:      {} N/m".format(TENSIONS_Nm))
        info(log, "Eigenmodes:    {}".format(N_EIGENMODES))
        info(log, "ALPHA_FR4:     {:.1e} /K  (ESTIMATE — update if Shengyi publishes)".format(ALPHA_FR4))
        info(log, "ALPHA_SS:      {:.1e} /K".format(ALPHA_SS))
        info(log, "delta_alpha:   {:.1e} /K".format(DELTA_ALPHA))
        info(log, "DT_bake_real:  {:.1f} K  ({:.0f} C -> {:.0f} C)".format(
            DT_BAKE_REAL, T_ROOM_C, T_BAKE_C))
        info(log, "DT_bake_model: {:.4f} K  (= {:.1f} * {:.1e} / {:.1e})".format(
            DT_BAKE_MODEL, DT_BAKE_REAL, DELTA_ALPHA, ALPHA_EFF))
        info(log, "Interpretation: lambda=1 -> buckles at {:.0f} C".format(T_BAKE_C))
        info(log, "=" * 60)

        if not os.path.isfile(FCSTD_PATH):
            info(log, "ERROR: FCStd not found: " + FCSTD_PATH)
            return

        info(log, "\nOpening: " + FCSTD_PATH)
        doc = FreeCAD.openDocument(FCSTD_PATH)

        info(log, "\nBuilding BC node sets ...")
        fem_mesh_obj = doc.getObject(OBJ_MESH)
        if fem_mesh_obj is None:
            info(log, "ERROR: mesh object '{}' not found".format(OBJ_MESH))
            return
        info(log, "  Total mesh nodes: {}".format(len(fem_mesh_obj.FemMesh.Nodes)))
        wall_nodes, pillar_nodes, surf = build_bc_node_sets(fem_mesh_obj, log)

        summary_rows = []
        failed = []

        for T_Nm in TENSIONS_Nm:
            info(log, "\n--- T = {} N/m  ({:.1f} N/cm) ---".format(
                T_Nm, T_Nm / 100.0))
            try:
                nodes_frd, modes, eigenvalues = run_one_tension(
                    doc, T_Nm, wall_nodes, pillar_nodes, surf, log)

                for mi, md in enumerate(modes):
                    ev = md.get("eigenvalue", float("nan"))
                    if math.isnan(ev):
                        continue
                    t_buckle = 25.0 + ev * DT_BAKE_REAL
                    survived = ev >= 1.0
                    info(log, "    Mode {:2d}: lambda={:.4f}  T_buckle~{:.0f} C  [{}]".format(
                        md.get("mode", mi + 1), ev, t_buckle,
                        "SURVIVED" if survived else "BUCKLES"))
                    csv_path = save_mode_csv(nodes_frd, md, T_Nm, mi, RESULTS_DIR)
                    info(log, "    Saved: " + csv_path)
                    summary_rows.append({
                        "T_Nm":       T_Nm,
                        "mode":       md.get("mode", mi + 1),
                        "eigenvalue": ev,
                    })

            except Exception as exc:
                info(log, "  FAILED: " + str(exc))
                info(log, traceback.format_exc())
                failed.append(T_Nm)

        if summary_rows:
            scsv = save_summary_csv(summary_rows, RESULTS_DIR)
            info(log, "\nSummary: " + scsv)

        info(log, "\n" + "=" * 60)
        if failed:
            info(log, "FAILED tensions: {}".format(failed))
        else:
            info(log, "All {} tensions completed.".format(len(TENSIONS_Nm)))
        info(log, "Results: " + RESULTS_DIR)
        info(log, "Log:     " + LOG_FILE)
        info(log, "=" * 60)


main()