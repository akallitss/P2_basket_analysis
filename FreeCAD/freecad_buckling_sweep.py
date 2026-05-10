"""
freecad_buckling_sweep.py
Buckling eigenvalue sweep for the P2 basket Bulk Micromegas mesh.

Physics:
  SS frame + SS mesh → same α → no differential thermal expansion during baking.
  Waves are caused by glue softening → pre-tension loss near the border.
  The wall BC spring stiffness k_glue represents the softened adhesive.

Method:
  For each (T_Nm, k_glue) combination:
    1. Static step: apply pre-tension via thermal analogy (same as tension sweep).
    2. Perturbation BUCKLE step: find the first N eigenmodes and eigenvalues.
       Reference load = 1 K additional cooling (= 1 unit of further tension).
       Eigenvalue λ > 0: mesh is stable; would need λ × 25 N/m more tension loss
                         before buckling (with α_eff = 50e-6/K).
       Eigenvalue λ ≤ 1: mesh buckles before full design tension is reached.

Run from repo root:
    "C:\\Program Files\\FreeCAD 0.20\\bin\\FreeCADCmd.exe" FreeCAD/freecad_buckling_sweep.py
Or let the launcher handle it (same as freecad_fem_setup.py).
"""

import os
import sys
import csv
import math
import time
import shutil
import tempfile
import datetime

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(SCRIPT_DIR, "results")
FCSTD_PATH  = r"C:\Temp\fc_v4_open.FCStd"

# ---------------------------------------------------------------------------
# Sweep parameters
# ---------------------------------------------------------------------------
TENSIONS_Nm   = [300, 500, 700, 1000]      # N/m — pre-tension values
K_GLUE_VALUES = [1e10, 1e7, 1e5, 1e3, 0.0] # N/m/node (1e10 ≈ rigid wall)
N_EIGENMODES  = 5                           # eigenmodes to extract per run

# ---------------------------------------------------------------------------
# Thermal analogy constants (must match freecad_fem_setup.py)
# ---------------------------------------------------------------------------
E_EFF   = 1.0e9    # Pa
T_EFF   = 0.5e-3   # m
ALPHA   = 50.0e-6  # /K  (effective, inflated — pre-tension analogy only)
T_REF   = 300.0    # K
# Reference load for BUCKLE step: 1 K of additional cooling
# → equivalent tension increment = E_EFF × T_EFF × ALPHA × 1 K
DELTA_T_REF_K = 1.0   # K per unit eigenvalue

# Physical constants for BC geometry
WALL_MM        = 1.0
PILLAR_PITCH_MM = 3.0
PILLAR_SNAP_MM  = 1.5
MESH_Z_LO_MM    = 0.05   # mm  lower z bound for mesh-layer nodes
MESH_Z_HI_MM    = 0.70   # mm  upper z bound

# FreeCAD object names
OBJ_MESH       = "FEMMeshGmsh"
OBJ_FEM        = "Analysis"
OBJ_SOLVER     = "SolverCcxTools"
OBJ_MATERIAL   = "MechanicalMaterial"
OBJ_CONSTRAINT = "ConstraintTemperature"

# ---------------------------------------------------------------------------
# Logging helper
# ---------------------------------------------------------------------------
def info(log, msg):
    print(msg)
    if log:
        log.write(msg + "\n")
        log.flush()


# ---------------------------------------------------------------------------
# BC node-set builder (identical to freecad_fem_setup.py)
# ---------------------------------------------------------------------------
def build_bc_node_sets(fem_mesh, log):
    nodes = fem_mesh.FemMesh.Nodes   # dict: nid → FreeCAD.Vector
    info(log, f"  Mesh-layer nodes: {len(nodes)}")

    # Collect mesh-layer nodes by z range
    surf = {nid: (v.x, v.y, v.z) for nid, v in nodes.items()
            if MESH_Z_LO_MM <= v.z <= MESH_Z_HI_MM}
    if not surf:
        info(log, f"  WARNING: no nodes in z=[{MESH_Z_LO_MM}, {MESH_Z_HI_MM}] mm — using all nodes")
        surf = {nid: (v.x, v.y, v.z) for nid, v in nodes.items()}

    xs = [c[0] for c in surf.values()]
    ys = [c[1] for c in surf.values()]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)

    # Wall nodes (1 mm border)
    wall_nodes = {nid for nid, (x, y, _) in surf.items()
                  if (x < xmin + WALL_MM or x > xmax - WALL_MM or
                      y < ymin + WALL_MM or y > ymax - WALL_MM)}
    info(log, f"  Wall BC nodes: {len(wall_nodes)}")

    # Pillar grid centres (inset by wall + half pitch)
    pcx = xmin + WALL_MM + PILLAR_PITCH_MM / 2
    pillar_nodes = set()
    while pcx <= xmax - WALL_MM:
        pcy = ymin + WALL_MM + PILLAR_PITCH_MM / 2
        while pcy <= ymax - WALL_MM:
            for nid, (x, y, _) in surf.items():
                if math.hypot(x - pcx, y - pcy) <= PILLAR_SNAP_MM:
                    pillar_nodes.add(nid)
            pcy += PILLAR_PITCH_MM
        pcx += PILLAR_PITCH_MM

    pillar_nodes -= wall_nodes
    info(log, f"  Pillar BC nodes: {len(pillar_nodes)}")
    return wall_nodes, pillar_nodes, surf


# ---------------------------------------------------------------------------
# CalculiX .inp patcher — injects spring/BC/BUCKLE into the generated file
# ---------------------------------------------------------------------------
def _inject_buckling_inp(inp_path, wall_nodes, pillar_nodes,
                         T_applied, k_glue, n_modes, surf, log):
    """
    Reads the FreeCAD-generated .inp, removes the *STEP/*STATIC section,
    and replaces it with:
      - SPRING1 elements for wall nodes (Ux, Uy) with stiffness k_glue
      - *BOUNDARY: Uz=0 at wall nodes + Uz=0 at pillar nodes
      - *STEP/*STATIC — pre-tension via thermal load
      - *STEP,PERTURBATION/*BUCKLE — eigenvalue extraction
    """
    with open(inp_path, "r") as f:
        raw = f.read()

    # ---- Strip everything from the first *STEP onward ----
    idx = raw.lower().find("\n*step")
    if idx == -1:
        idx = raw.lower().find("*step")
    preamble = raw[:idx] if idx != -1 else raw

    wall_list   = sorted(wall_nodes)
    pillar_list = sorted(pillar_nodes)

    # Find max existing element ID so SPRING1 IDs are compact (not 10_000_000+).
    # CalculiX allocates arrays sized by max element ID, so high IDs → OOM crash.
    max_eid = 0
    in_elem = False
    for line in preamble.splitlines():
        s = line.strip()
        if s.upper().startswith("*ELEMENT") and not s.startswith("**"):
            in_elem = True; continue
        if s.startswith("*") and not s.startswith("**"):
            in_elem = False; continue
        if in_elem and s:
            try:
                eid = int(s.split(",")[0])
                if eid > max_eid:
                    max_eid = eid
            except ValueError:
                pass
    spring_base = max_eid + 1

    lines = [preamble.rstrip()]

    # ---- SPRING1 elements for Ux and Uy at wall nodes ----
    if k_glue > 0:
        # Ux springs  (IDs: spring_base … spring_base + len(wall)-1)
        lines.append("\n*ELEMENT, TYPE=SPRING1, ELSET=ESPRING_UX")
        for i, nid in enumerate(wall_list):
            lines.append(f"{spring_base + i}, {nid}")
        lines.append("*SPRING, ELSET=ESPRING_UX")
        lines.append("1,")
        lines.append(f"{k_glue:.6E}")
        # Uy springs  (IDs continue immediately after Ux block)
        uy_base = spring_base + len(wall_list)
        lines.append("*ELEMENT, TYPE=SPRING1, ELSET=ESPRING_UY")
        for i, nid in enumerate(wall_list):
            lines.append(f"{uy_base + i}, {nid}")
        lines.append("*SPRING, ELSET=ESPRING_UY")
        lines.append("2,")
        lines.append(f"{k_glue:.6E}")

    # ---- BOUNDARY: Uz=0 at wall nodes ----
    lines.append("*BOUNDARY")
    for nid in wall_list:
        if k_glue <= 0:
            # Free glue: only fix Uz, leave Ux/Uy unconstrained
            lines.append(f"{nid}, 3, 3, 0.0")
        else:
            # Spring BCs handle Ux/Uy; only fix Uz here
            lines.append(f"{nid}, 3, 3, 0.0")

    # ---- BOUNDARY: Uz=0 at pillar nodes ----
    for nid in pillar_list:
        lines.append(f"{nid}, 3, 3, 0.0")

    # ---- STEP 1: Static pre-tension via thermal cooling ----
    lines.append("*STEP")
    lines.append("*STATIC")
    lines.append("*TEMPERATURE")
    lines.append(f"NALL, {T_applied:.6f}")
    lines.append("*NODE FILE")
    lines.append("U")
    lines.append("*END STEP")

    # ---- STEP 2: Perturbation BUCKLE ----
    # *TEMPERATURE is not allowed in *BUCKLE (CalculiX restriction).
    # Use inward *CLOAD at wall nodes instead (1 N per node in the inward
    # radial direction).  Sign of lambda is the stability indicator:
    #   lambda > 0: stable under this compression direction
    #   lambda < 0: already past critical (buckled at equilibrium)
    wall_x = [surf[n][0] for n in wall_nodes if n in surf]
    wall_y = [surf[n][1] for n in wall_nodes if n in surf]
    cx = (min(wall_x) + max(wall_x)) / 2.0
    cy = (min(wall_y) + max(wall_y)) / 2.0
    lines.append("*STEP, PERTURBATION")
    lines.append("*BUCKLE")
    lines.append(f"{n_modes}")
    lines.append("*CLOAD")
    for nid in sorted(wall_nodes):
        if nid not in surf:
            continue
        x, y = surf[nid][0], surf[nid][1]
        dx, dy = cx - x, cy - y
        dist = math.sqrt(dx * dx + dy * dy)
        if dist < 1e-9:
            continue
        lines.append(f"{nid}, 1, {dx/dist:.6f}")
        lines.append(f"{nid}, 2, {dy/dist:.6f}")
    lines.append("*NODE FILE")
    lines.append("U")
    lines.append("*END STEP")

    with open(inp_path, "w") as f:
        f.write("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# FRD parser — mode shapes (eigenvectors from BUCKLE step)
# ---------------------------------------------------------------------------
def _parse_buckle_frd(path):
    """
    Parse displacement eigenvectors from a CalculiX buckling .frd file.
    Returns: list of dicts, one per mode:
      {'mode': int, 'nodes': {nid: (x,y,z)}, 'disp': {nid: (ux,uy,uz)}}
    Also reads node coordinates from the first 2C block.
    """
    nodes = {}        # nid → (x, y, z) mm
    modes = []        # list of {'mode': int, 'disp': {nid: (ux,uy,uz)}}
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
                # Start of a result block — commit previous mode if any
                if current_mode is not None and current_disp:
                    modes.append({"mode": current_mode, "disp": current_disp})
                current_mode = None
                current_disp = {}
                state = "pending"; continue
            if line[:3] == " -4":
                if state == "pending":
                    content = line[3:].upper().strip()
                    if any(k in content for k in ("DISP", "D1 ", "D2 ", "D3 ", " U ")):
                        # Try to extract mode number from the 100C line header
                        # (it's encoded in the step/increment counter)
                        current_mode = (len(modes) + 1) if current_mode is None else current_mode
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
    Parse buckling eigenvalues from the CalculiX .dat file.
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
            if "BUCKLINGFACTOR" in compact or "BUCKL" in compact and "FACTOR" in compact:
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
                        eigenvalues.append((mode, val))
                    except ValueError:
                        continue  # skip sub-headers (MODE NO / FACTOR)
    if not eigenvalues:
        # Diagnostic: dump last 60 lines of .dat so we can see the real header
        print("  DAT DIAGNOSTIC — last {} lines of {}:".format(
            len(header_lines), dat_path))
        for l in header_lines:
            print("    |" + l)
    return eigenvalues


# ---------------------------------------------------------------------------
# Run one (T, k_glue) combination
# ---------------------------------------------------------------------------
def run_one_case(doc, T_Nm, k_glue, wall_nodes, pillar_nodes, surf, log):
    import FreeCAD
    import femtools.ccxtools as ccx

    delta_T = T_Nm / (E_EFF * T_EFF * ALPHA)
    T_applied = T_REF - delta_T
    info(log, f"  T={T_Nm} N/m  k_glue={k_glue:.0E}  ΔT={delta_T:.2f} K  "
              f"T_applied={T_applied:.2f} K")

    # Set up FEA object
    fea = ccx.FemToolsCcx(analysis=doc.getObject(OBJ_FEM),
                           solver=doc.getObject(OBJ_SOLVER))
    fea.update_objects()
    fea.setup_working_dir()
    fea.setup_ccx()

    # Monkey-patch write_inp_file so the *BUCKLE injection survives the
    # internal write_inp_file() call that fea.run() makes before launching ccx.
    orig_write = fea.write_inp_file
    def _write_and_inject():
        result = orig_write()
        inp_path = getattr(fea, "inp_file_name", None)
        if not inp_path or not os.path.isfile(inp_path):
            wdir = getattr(fea, "working_dir", None)
            if wdir and os.path.isdir(wdir):
                for root, _, files in os.walk(wdir):
                    for fname in files:
                        if fname.endswith(".inp"):
                            inp_path = os.path.join(root, fname)
                            break
        if inp_path:
            info(log, f"  Patching .inp for buckling: {inp_path}")
            _inject_buckling_inp(inp_path, wall_nodes, pillar_nodes,
                                 T_applied, k_glue, N_EIGENMODES, surf, log)
        else:
            info(log, "  WARNING: cannot find .inp to patch")
        return result
    fea.write_inp_file = _write_and_inject

    # Prevent load_results from accumulating CCX_Results objects
    fea.load_results = lambda: None

    info(log, "  Running CalculiX (BUCKLE) ...")
    msg = fea.run()
    if msg and not isinstance(msg, bool):
        raise RuntimeError("CalculiX error:\n" + str(msg))

    # Locate .frd and .dat files
    wdir = getattr(fea, "working_dir", None)
    frd_path = dat_path = None
    if wdir and os.path.isdir(wdir):
        for root, _, files in os.walk(wdir):
            for fname in files:
                fp = os.path.join(root, fname)
                if fname.endswith(".frd"): frd_path = fp
                if fname.endswith(".dat"): dat_path = fp

    if not frd_path:
        raise RuntimeError("Cannot find .frd in working dir: " + str(wdir))

    info(log, f"  Parsing FRD: {frd_path}")
    nodes_frd, modes = _parse_buckle_frd(frd_path)
    info(log, f"  Found {len(modes)} mode(s), {len(nodes_frd)} nodes")

    info(log, f"  Parsing DAT: {dat_path}")
    eigenvalues = _parse_buckle_dat(dat_path)
    info(log, f"  Eigenvalues: {eigenvalues}")

    # Re-number modes from eigenvalues if available
    for i, (mode_n, ev) in enumerate(eigenvalues):
        if i < len(modes):
            modes[i]["mode"]       = mode_n
            modes[i]["eigenvalue"] = ev
            # Physical interpretation:
            # ev × DELTA_T_REF_K cooling → ev × E_EFF*T_EFF*ALPHA N/m more tension
            modes[i]["tension_margin_Nm"] = ev * E_EFF * T_EFF * ALPHA * DELTA_T_REF_K

    return nodes_frd, modes, eigenvalues


# ---------------------------------------------------------------------------
# Save results
# ---------------------------------------------------------------------------
def save_mode_csv(nodes, mode_data, T_Nm, k_glue, mode_idx, results_dir):
    k_tag  = f"{k_glue:.0E}".replace("+", "p").replace("-", "m")
    fname  = f"buckle_T{T_Nm:04d}Nm_k{k_tag}_mode{mode_idx+1:02d}.csv"
    path   = os.path.join(results_dir, fname)
    disp   = mode_data.get("disp", {})
    ev     = mode_data.get("eigenvalue", float("nan"))
    tm     = mode_data.get("tension_margin_Nm", float("nan"))

    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["node_id", "x_mm", "y_mm", "z_mm",
                    "Ux_norm", "Uy_norm", "Uz_norm",
                    "eigenvalue", "tension_margin_Nm"])
        # Normalise mode shape to max|Uz|=1 for comparability
        uzs  = [d[2] for d in disp.values()]
        scale = max(abs(v) for v in uzs) if uzs else 1.0
        if scale == 0: scale = 1.0
        for nid, (x, y, z) in nodes.items():
            if nid in disp:
                ux, uy, uz = disp[nid]
                w.writerow([nid, f"{x:.4f}", f"{y:.4f}", f"{z:.4f}",
                             f"{ux/scale:.6f}", f"{uy/scale:.6f}", f"{uz/scale:.6f}",
                             f"{ev:.6f}", f"{tm:.4f}"])
    return path


def save_summary_csv(summary_rows, results_dir):
    path = os.path.join(results_dir, "buckle_summary.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["T_Nm", "T_Ncm", "k_glue", "mode", "eigenvalue",
                    "tension_margin_Nm", "stable"])
        for row in summary_rows:
            w.writerow([row["T_Nm"], f"{row['T_Nm']/100:.1f}",
                        f"{row['k_glue']:.2E}",
                        row["mode"], f"{row['eigenvalue']:.6f}",
                        f"{row['tension_margin_Nm']:.4f}",
                        "yes" if row["eigenvalue"] > 0 else "no"])
    return path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    import FreeCAD
    import FreeCADGui   # noqa — needed to init module even headless

    os.makedirs(RESULTS_DIR, exist_ok=True)
    log_path = os.path.join(SCRIPT_DIR, "buckling_sweep.log")

    with open(log_path, "w", encoding="utf-8") as log:
        banner = ("=" * 60 + "\n" +
                  "P2 Basket FEM buckling sweep\n" +
                  f"Tensions: {TENSIONS_Nm} N/m\n" +
                  f"k_glue:   {K_GLUE_VALUES}\n" +
                  f"Modes:    {N_EIGENMODES}\n" +
                  "=" * 60)
        info(log, banner)

        info(log, f"\nOpening: {FCSTD_PATH}")
        doc = FreeCAD.openDocument(FCSTD_PATH)

        # Build BC node sets once (geometry doesn't change)
        info(log, "Computing wall + pillar BC node sets ...")
        fem_mesh = doc.getObject(OBJ_MESH)
        if fem_mesh is None:
            raise RuntimeError(f"Mesh object '{OBJ_MESH}' not found")
        wall_nodes, pillar_nodes, surf = build_bc_node_sets(fem_mesh, log)

        summary_rows = []
        total = len(TENSIONS_Nm) * len(K_GLUE_VALUES)
        run_n = 0

        for T_Nm in TENSIONS_Nm:
            for k_glue in K_GLUE_VALUES:
                run_n += 1
                info(log, f"\n--- Run {run_n}/{total}: "
                          f"T={T_Nm} N/m  k_glue={k_glue:.0E} ---")
                try:
                    nodes_frd, modes, eigenvalues = run_one_case(
                        doc, T_Nm, k_glue, wall_nodes, pillar_nodes, surf, log)

                    # Save mode shapes
                    for mi, md in enumerate(modes):
                        save_mode_csv(nodes_frd, md, T_Nm, k_glue, mi, RESULTS_DIR)
                        ev = md.get("eigenvalue", float("nan"))
                        tm = md.get("tension_margin_Nm", float("nan"))
                        summary_rows.append({
                            "T_Nm": T_Nm, "k_glue": k_glue,
                            "mode": md.get("mode", mi + 1),
                            "eigenvalue": ev,
                            "tension_margin_Nm": tm,
                        })
                        info(log, f"    Mode {md.get('mode', mi+1)}: "
                                  f"λ={ev:.4f}  "
                                  f"margin={tm:.1f} N/m  "
                                  f"({'STABLE' if ev > 0 else 'BUCKLED'})")
                except Exception as e:
                    info(log, f"  ERROR: {e}")
                    import traceback
                    info(log, traceback.format_exc())

        summary_path = save_summary_csv(summary_rows, RESULTS_DIR)
        info(log, "\n" + "=" * 60)
        info(log, f"Buckling sweep complete.  Summary: {summary_path}")
        info(log, f"Results directory: {RESULTS_DIR}")
        info(log, "=" * 60)

    print(f"\nLog: {log_path}")
    print("Next step: run visualize_buckling.py to plot stability map and mode shapes.")


main()