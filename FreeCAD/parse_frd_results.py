"""
parse_frd_results.py  –  Standalone CalculiX .frd reader

No FreeCAD required.  Reads raw CalculiX binary result files and writes
the same CSV format that freecad_fem_setup.py produces, so visualize_sweep.py
can consume either output transparently.

Usage
-----
# Single file → CSV next to the .frd
python parse_frd_results.py path/to/ccx.frd

# Directory scan → one CSV per .frd found, summary printed to stdout
python parse_frd_results.py path/to/working_dir/

# Explicit output path
python parse_frd_results.py ccx.frd --out results/sweep_T0300Nm.csv

# Tell the script what tension this run was (written into summary line)
python parse_frd_results.py ccx.frd --tension 300

FRD format notes (CalculiX 2.x)
---------------------------------
All data lines start with ' -1' (3 chars), followed by:
  node_id  : 10-char right-justified integer  (cols 4-13)
  values   : one or more 12-char E12.5 fields (cols 14+)

Block headers:
  '    2C'   → node coordinate block
  '  100C'   → result block  (DISP = displacements, STRESS = stress tensor)
  '  -3 '    → end of block
  '  -4 '    → component descriptor inside a result block
  '9999'     → end of file

Von Mises from stress tensor (Voigt notation: Sxx Syy Szz Sxy Sxz Syz):
  σ_vm = sqrt( Sxx²+Syy²+Szz² - Sxx·Syy - Syy·Szz - Szz·Sxx
               + 3·(Sxy²+Sxz²+Syz²) )
"""

import os
import sys
import csv
import math
import glob
import argparse


# ---------------------------------------------------------------------------
# FRD parser
# ---------------------------------------------------------------------------

def _parse_fixed_values(line, start=13):
    """
    Extract 12-char E12.5 float fields starting at column `start` (0-indexed).
    Falls back to whitespace-split if fixed-width parsing yields nothing.
    """
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
        # fallback: split the tail
        tail = raw[start:].strip()
        for tok in tail.split():
            try:
                vals.append(float(tok))
            except ValueError:
                pass
    return vals


def parse_frd(path):
    """
    Parse a CalculiX .frd file.

    Returns
    -------
    nodes  : dict  node_id -> (x_mm, y_mm, z_mm)
    disp   : dict  node_id -> (Ux_mm, Uy_mm, Uz_mm)
    stress : dict  node_id -> (Sxx, Syy, Szz, Sxy, Sxz, Syz)  [MPa when geometry in mm]
    meta   : dict  with keys 'name', 'n_nodes', 'n_disp', 'n_stress'
    """
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

            # ---- block openers ------------------------------------------------
            if key6 == "    2C":
                state = "nodes"
                continue

            # Element connectivity block — skip data lines
            if key6 in ("    3C", "    2P"):
                state = "elements"
                continue

            if key6.startswith("  100C"):
                # Result type is NOT in this line — determined by the -4 descriptor below
                state = "pending"
                continue

            # ---- component descriptor: reveals the result type ----------------
            if line[:3] == " -4":
                if state == "pending":
                    content = line[3:].upper().strip()
                    if any(k in content for k in ("DISP", "D1 ", "D2 ", "D3 ", " U ")):
                        state = "disp"
                    elif any(k in content for k in ("STRESS", "SXX", "SYY", "SZZ")):
                        state = "stress"
                    elif "MISES" in content:
                        state = "skip"   # single-component; not needed with full tensor
                    else:
                        state = "skip"
                # subsequent -4 lines for same block: state already set
                continue

            # ---- end of block -------------------------------------------------
            if line[:3] == " -3" or line.startswith("  -3"):
                state = None
                continue

            # end-of-file marker
            if line.strip() == "9999":
                break

            # ---- data lines ---------------------------------------------------
            if line[:3] == " -1" and state not in (None, "skip", "pending", "elements"):
                try:
                    nid = int(line[3:13])
                except ValueError:
                    continue

                vals = _parse_fixed_values(line, start=13)

                if state == "nodes" and len(vals) >= 3:
                    nodes[nid] = (vals[0], vals[1], vals[2])

                elif state == "disp" and len(vals) >= 3:
                    # DISP block has D1 D2 D3 [ALL]
                    disp[nid] = (vals[0], vals[1], vals[2])

                elif state == "stress" and len(vals) >= 6:
                    # Sxx Syy Szz Sxy Sxz Syz
                    stress[nid] = tuple(vals[:6])

    meta = {
        "name":     os.path.basename(path),
        "n_nodes":  len(nodes),
        "n_disp":   len(disp),
        "n_stress": len(stress),
    }
    return nodes, disp, stress, meta


# ---------------------------------------------------------------------------
# Von Mises
# ---------------------------------------------------------------------------

def von_mises(sxx, syy, szz, sxy, sxz, syz):
    vm2 = (sxx**2 + syy**2 + szz**2
           - sxx*syy - syy*szz - szz*sxx
           + 3.0*(sxy**2 + sxz**2 + syz**2))
    return math.sqrt(max(vm2, 0.0))


# ---------------------------------------------------------------------------
# CSV writer (same columns as freecad_fem_setup.py)
# ---------------------------------------------------------------------------

def write_csv(csv_path, nodes, disp, stress):
    """
    Write node_id, x_mm, y_mm, z_mm, Ux_mm, Uy_mm, Uz_mm, vonMises_MPa.
    Only writes nodes that have both displacement and coordinate data.
    Stress is optional – if missing, vonMises_MPa is left blank.
    """
    common_ids = sorted(set(nodes) & set(disp))
    if not common_ids:
        raise RuntimeError(
            "No nodes have both coordinate and displacement data.\n"
            "  nodes={}, disp={}, stress={}".format(len(nodes), len(disp), len(stress))
        )

    os.makedirs(os.path.dirname(os.path.abspath(csv_path)), exist_ok=True)

    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["node_id", "x_mm", "y_mm", "z_mm",
                    "Ux_mm", "Uy_mm", "Uz_mm", "vonMises_MPa"])
        for nid in common_ids:
            x, y, z = nodes[nid]
            ux, uy, uz = disp[nid]
            if nid in stress:
                vm = von_mises(*stress[nid])
            else:
                vm = ""
            w.writerow([
                nid,
                round(x,  4), round(y,  4), round(z,  4),
                round(ux, 6), round(uy, 6), round(uz, 6),
                round(vm, 4) if vm != "" else "",
            ])

    return len(common_ids)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def _default_csv_path(frd_path, tension_nm=None):
    base_dir = os.path.dirname(os.path.abspath(frd_path))
    if tension_nm is not None:
        return os.path.join(base_dir, "sweep_T{:04d}Nm.csv".format(int(tension_nm)))
    stem = os.path.splitext(os.path.basename(frd_path))[0]
    return os.path.join(base_dir, stem + "_parsed.csv")


def process_one(frd_path, csv_path=None, tension_nm=None, verbose=True):
    if csv_path is None:
        csv_path = _default_csv_path(frd_path, tension_nm)

    if verbose:
        print("Parsing: {}".format(frd_path))

    nodes, disp, stress, meta = parse_frd(frd_path)

    if verbose:
        print("  nodes={n_nodes}  disp={n_disp}  stress={n_stress}".format(**meta))

    if meta["n_disp"] == 0:
        print("  WARNING: no displacement data found – check .frd contains DISP block")
        return None

    n_written = write_csv(csv_path, nodes, disp, stress)

    if verbose:
        uzs = [disp[nid][2] for nid in nodes if nid in disp]
        if uzs:
            max_abs = max(abs(v) for v in uzs) * 1000  # mm → µm
            print("  max|Uz| = {:.2f} µm".format(max_abs))
        print("  Wrote {} rows → {}".format(n_written, csv_path))

    return csv_path


def main():
    ap = argparse.ArgumentParser(
        description="Parse CalculiX .frd file(s) → CSV (no FreeCAD needed)"
    )
    ap.add_argument("path",
                    help=".frd file or directory containing .frd files")
    ap.add_argument("--out", metavar="CSV",
                    help="output CSV path (single-file mode only)")
    ap.add_argument("--tension", metavar="N_m", type=float,
                    help="mesh tension in N/m – used to name the output CSV automatically")
    args = ap.parse_args()

    target = args.path

    if os.path.isdir(target):
        frd_files = sorted(glob.glob(os.path.join(target, "*.frd")))
        if not frd_files:
            print("No .frd files found in", target)
            sys.exit(1)
        print("Found {} .frd file(s)".format(len(frd_files)))
        for frd in frd_files:
            process_one(frd, tension_nm=args.tension)
    elif os.path.isfile(target):
        process_one(target, csv_path=args.out, tension_nm=args.tension)
    else:
        print("ERROR: path not found:", target)
        sys.exit(1)


main()