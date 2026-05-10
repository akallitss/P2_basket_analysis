"""
recover_from_frd.py  –  Extract sweep results from existing CalculiX FRD files

When freecad_fem_setup.py's result-reading step fails, CalculiX still writes
a valid .frd file to a temp directory.  This script:

  1. Scans all fcfem_* directories in %TEMP%
  2. Reads the .inp file to find the temperature BC  → identifies tension
  3. Runs the FRD parser on each directory
  4. Writes the per-step CSVs and sweep_summary.csv to FreeCAD/results/

Temperature → Tension mapping (from model parameters):
  T_applied = T_ref - ΔT = 300 - T_Nm/25
  → T_Nm = (300 - T_applied) * 25
"""

import os
import re
import sys
import glob
import csv
import math

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR  = os.path.join(SCRIPT_DIR, "results")
TEMP_DIR     = os.environ.get("TEMP", r"C:\Users\alexa\AppData\Local\Temp")

T_REF        = 300.0   # K (ConstraintInitialTemperature)
DT_PER_NM    = 1.0 / 25.0   # K·m/N

# Expected tensions in sweep (N/m) — used to snap detected values
EXPECTED_Nm  = [300, 400, 500, 600, 700, 800, 900, 1000]

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def find_fcfem_dirs(temp_dir):
    pattern = os.path.join(temp_dir, "fcfem_*")
    dirs = sorted(glob.glob(pattern))
    return [d for d in dirs if os.path.isdir(d)]


def read_tension_from_inp(fcfem_dir):
    """
    Read the .inp file in the directory, find the temperature BC,
    and return the corresponding tension in N/m.

    CalculiX temperature DOF is 11. FreeCAD writes lines like:
        <nodeset>, 11, 11, 288.0
    or (older format):
        NT11, <nodeset>, 288.0
    """
    inp_files = glob.glob(os.path.join(fcfem_dir, "*.inp"))
    if not inp_files:
        return None

    inp_path = inp_files[0]
    try:
        with open(inp_path, "r", errors="replace") as f:
            content = f.read()
    except Exception:
        return None

    # Pattern 1: "..., 11, 11, 288.0"  (BOUNDARY with DOF 11 = temperature)
    m = re.search(r",\s*11\s*,\s*11\s*,\s*([\d.eE+\-]+)", content)
    if m:
        t_applied = float(m.group(1))
        t_nm = (T_REF - t_applied) / DT_PER_NM
        return t_applied, t_nm

    # Pattern 2: "*TEMPERATURE" block
    m = re.search(r"\*TEMPERATURE[^\n]*\n[^,\n]+,\s*([\d.eE+\-]+)", content)
    if m:
        t_applied = float(m.group(1))
        t_nm = (T_REF - t_applied) / DT_PER_NM
        return t_applied, t_nm

    return None


def snap_to_expected(t_nm_raw, tolerance=10.0):
    """Round detected tension to nearest expected value (within tolerance N/m)."""
    best = min(EXPECTED_Nm, key=lambda x: abs(x - t_nm_raw))
    if abs(best - t_nm_raw) <= tolerance:
        return best
    return None


# ---------------------------------------------------------------------------
# FRD parser (inline, same logic as parse_frd_results.py)
# ---------------------------------------------------------------------------

def _parse_fixed_values(line, start=13):
    vals = []
    pos = start
    raw = line.rstrip("\n\r")
    while pos + 12 <= len(raw):
        field = raw[pos:pos + 12]
        try:
            vals.append(float(field))
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


def parse_frd(path):
    nodes  = {}
    disp   = {}
    stress = {}
    # 'nodes' | 'elements' | 'pending' | 'disp' | 'stress' | 'skip' | None
    state  = None

    with open(path, "r", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip("\n\r")
            if not line:
                continue
            key6 = line[:6]

            if key6 == "    2C":
                state = "nodes"
                continue
            if key6 in ("    3C", "    2P"):
                state = "elements"
                continue
            if key6.startswith("  100C"):
                # Result type NOT in this line — determined by the -4 descriptor below
                state = "pending"
                continue
            if line[:3] == " -4":
                if state == "pending":
                    content = line[3:].upper().strip()
                    if any(k in content for k in ("DISP", "D1 ", "D2 ", "D3 ", " U ")):
                        state = "disp"
                    elif any(k in content for k in ("STRESS", "SXX", "SYY", "SZZ")):
                        state = "stress"
                    else:
                        state = "skip"
                continue
            if line[:3] == " -3" or line.startswith("  -3"):
                state = None
                continue
            if line.strip() == "9999":
                break
            if line[:3] == " -1" and state not in (None, "skip", "pending", "elements"):
                try:
                    nid = int(line[3:13])
                except ValueError:
                    continue
                vals = _parse_fixed_values(line, start=13)
                if state == "nodes" and len(vals) >= 3:
                    nodes[nid] = (vals[0], vals[1], vals[2])
                elif state == "disp" and len(vals) >= 3:
                    disp[nid] = (vals[0], vals[1], vals[2])
                elif state == "stress" and len(vals) >= 6:
                    stress[nid] = tuple(vals[:6])

    return nodes, disp, stress


def von_mises(sxx, syy, szz, sxy, sxz, syz):
    vm2 = (sxx**2 + syy**2 + szz**2
           - sxx*syy - syy*szz - szz*sxx
           + 3.0*(sxy**2 + sxz**2 + syz**2))
    return math.sqrt(max(vm2, 0.0))


def write_step_csv(t_nm, nodes, disp, stress, results_dir):
    fname = os.path.join(results_dir, "sweep_T{:04d}Nm.csv".format(t_nm))
    common = sorted(set(nodes) & set(disp))
    if not common:
        return None, 0
    with open(fname, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["node_id", "x_mm", "y_mm", "z_mm",
                    "Ux_mm", "Uy_mm", "Uz_mm", "vonMises_MPa"])
        for nid in common:
            x, y, z = nodes[nid]
            ux, uy, uz = disp[nid]
            vm = von_mises(*stress[nid]) if nid in stress else ""
            w.writerow([nid,
                        round(x, 4), round(y, 4), round(z, 4),
                        round(ux, 6), round(uy, 6), round(uz, 6),
                        round(vm, 4) if vm != "" else ""])
    return fname, len(common)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)

    print("Scanning: {}".format(TEMP_DIR))
    fcfem_dirs = find_fcfem_dirs(TEMP_DIR)
    print("Found {} fcfem_* directories".format(len(fcfem_dirs)))

    # Map: t_nm -> (frd_path, t_applied)
    found = {}

    for d in fcfem_dirs:
        frd_files = glob.glob(os.path.join(d, "*.frd"))
        if not frd_files:
            continue
        frd_path = frd_files[0]

        result = read_tension_from_inp(d)
        if result is None:
            print("  [skip] {}: no .inp or no temperature BC found".format(
                os.path.basename(d)))
            continue

        t_applied, t_nm_raw = result
        t_nm = snap_to_expected(t_nm_raw)
        if t_nm is None:
            print("  [skip] {}: T_applied={:.1f} K → {:.0f} N/m not in sweep".format(
                os.path.basename(d), t_applied, t_nm_raw))
            continue

        frd_size_mb = os.path.getsize(frd_path) / 1e6
        print("  {} → T={} N/m ({:.0f} N/cm)  frd={:.1f} MB".format(
            os.path.basename(d), t_nm, t_nm / 100.0, frd_size_mb))

        # Keep newest file if duplicate tension
        if t_nm not in found or os.path.getmtime(frd_path) > os.path.getmtime(found[t_nm][0]):
            found[t_nm] = (frd_path, t_applied)

    if not found:
        print("\nNo matching FRD files found. Check TEMP_DIR and ensure the sweep ran.")
        sys.exit(1)

    print("\nProcessing {} unique tension steps ...".format(len(found)))
    summary_rows = []

    for t_nm in sorted(found):
        frd_path, t_applied = found[t_nm]
        print("\n  T={} N/m ({:.0f} N/cm) → {}".format(
            t_nm, t_nm / 100.0, frd_path))
        print("    Parsing FRD ...")
        nodes, disp, stress = parse_frd(frd_path)
        print("    nodes={} disp={} stress={}".format(
            len(nodes), len(disp), len(stress)))

        if len(disp) == 0:
            print("    WARNING: no displacement data – skipping")
            continue

        csv_path, n_rows = write_step_csv(t_nm, nodes, disp, stress, RESULTS_DIR)
        if csv_path is None:
            print("    WARNING: no common nodes – skipping")
            continue

        uzs = [disp[nid][2] for nid in nodes if nid in disp]
        max_abs_uz = max(abs(v) for v in uzs)
        max_uz = max(uzs)
        min_uz = min(uzs)
        vms = [von_mises(*stress[nid]) for nid in nodes if nid in stress]
        max_vm = max(vms) if vms else 0.0

        print("    max|Uz| = {:.4f} mm  ({:.1f} µm)".format(
            max_abs_uz, max_abs_uz * 1000))
        print("    max vonMises = {:.4f} MPa".format(max_vm))
        print("    Wrote {} rows → {}".format(n_rows, csv_path))

        summary_rows.append([
            t_nm,
            round(t_nm / 100.0, 1),
            round((T_REF - t_applied), 2),
            round(max_uz, 6),
            round(min_uz, 6),
            round(max_abs_uz, 6),
            round(max_vm, 4),
        ])

    if summary_rows:
        scsv = os.path.join(RESULTS_DIR, "sweep_summary.csv")
        with open(scsv, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["T_Nm", "T_Ncm", "delta_T_K",
                        "max_Uz_mm", "min_Uz_mm",
                        "max_abs_Uz_mm", "max_vonMises_MPa"])
            for row in summary_rows:
                w.writerow(row)
        print("\nSummary: {} → {} steps".format(scsv, len(summary_rows)))

    print("\nDone.  Results in: {}".format(RESULTS_DIR))


main()