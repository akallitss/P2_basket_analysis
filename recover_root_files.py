#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ── ROOT bootstrap ────────────────────────────────────────────────────────────
# This block runs before anything else. If ROOT is not importable it finds
# thisroot.sh automatically and re-launches the script with the correct env.
import sys, os, subprocess, shutil

def _find_thisroot():
    """Return path to thisroot.sh, or None if not found."""
    # 1) root-config is already on PATH → derive prefix from it
    rc = shutil.which("root-config")
    if rc:
        try:
            prefix = subprocess.check_output([rc, "--prefix"],
                                             stderr=subprocess.DEVNULL).decode().strip()
            candidate = os.path.join(prefix, "bin", "thisroot.sh")
            if os.path.isfile(candidate):
                return candidate
        except Exception:
            pass

    # 2) Common hard-coded install locations (add yours here if needed)
    common = [
        "/usr/local/root/bin/thisroot.sh",
        "/opt/root/bin/thisroot.sh",
        "/usr/share/root/bin/thisroot.sh",
        os.path.expanduser("~/root/bin/thisroot.sh"),
        os.path.expanduser("~/Software/root/bin/thisroot.sh"),
    ]
    for p in common:
        if os.path.isfile(p):
            return p

    # 3) Search under /usr, /opt, /local, ~/Software for any thisroot.sh
    search_roots = ["/usr", "/opt", "/local",
                    os.path.expanduser("~/Software"),
                    os.path.expanduser("~")]
    for base in search_roots:
        for dirpath, _, filenames in os.walk(base):
            if "thisroot.sh" in filenames:
                return os.path.join(dirpath, "thisroot.sh")

    return None

def _relaunch_with_root(thisroot_sh):
    """
    Source thisroot.sh in a subshell, capture the resulting environment,
    then re-exec this script inside that environment.
    """
    print(f"[ROOT bootstrap] Sourcing: {thisroot_sh}")
    # Capture the env after sourcing thisroot.sh
    cmd = f'source "{thisroot_sh}" > /dev/null 2>&1 && env -0'
    try:
        raw = subprocess.check_output(cmd, shell=True, executable="/bin/bash")
    except subprocess.CalledProcessError as e:
        print(f"[ROOT bootstrap] Failed to source thisroot.sh: {e}")
        sys.exit(1)

    env = {}
    for item in raw.split(b"\x00"):
        if b"=" in item:
            k, _, v = item.partition(b"=")
            env[k.decode(errors="replace")] = v.decode(errors="replace")

    # Re-exec with the enriched environment (same interpreter, same args)
    os.execve(sys.executable, [sys.executable] + sys.argv, env)

# Only bootstrap if ROOT is not already importable
try:
    import ROOT  # noqa: F401 — just checking
except ModuleNotFoundError:
    _thisroot = _find_thisroot()
    if _thisroot is None:
        print(
            "[ERROR] Cannot find thisroot.sh anywhere.\n"
            "  Please source it manually:\n"
            "    source $(root-config --prefix)/bin/thisroot.sh\n"
            "  Or add its directory to the search list in _find_thisroot()."
        )
        sys.exit(1)
    _relaunch_with_root(_thisroot)
    # os.execve above replaces the process — code below only runs
    # in the re-launched interpreter where ROOT is available.

# ── end ROOT bootstrap ────────────────────────────────────────────────────────
"""
recover_root_files.py
─────────────────────
Recovers ROOT files that were not properly closed (e.g. DAQ crash),
making them readable by uproot.

Problem:
    When a ROOT file is not closed cleanly, its file-header directory
    record is never written. uproot relies on that record and returns
    empty keys []. PyROOT/TBrowser can still read such files because
    ROOT's C++ layer runs an automatic recovery scan — this script
    replicates that and writes out a clean copy.

Usage (CLI):
    # Interactive mode — lists available runs and lets you pick
    python recover_root_files.py \
        --base   /drf/projets/clas12/cern_202511_p2_alinx \
        --output /drf/projets/clas12/cern_202511_p2_alinx_recovered

    # Single specific run (non-interactive)
    python recover_root_files.py \
        --input  /drf/projets/clas12/cern_202511_p2_alinx/run_101 \
        --output /drf/projets/clas12/cern_202511_p2_alinx_recovered/run_101

    # Single file
    python recover_root_files.py \
        --input  /path/to/bad.root \
        --output /path/to/fixed.root

    # All runs at once (non-interactive)
    python recover_root_files.py \
        --input  /drf/projets/clas12/cern_202511_p2_alinx \
        --output /drf/projets/clas12/cern_202511_p2_alinx_recovered \
        --recursive

    # Dry run (just report, don't write anything)
    python recover_root_files.py \
        --input  /drf/projets/clas12/cern_202511_p2_alinx/run_101 \
        --output /tmp/test/run_101 \
        --dry-run

    # Custom file glob pattern
    python recover_root_files.py \
        --input  /path/to/dir \
        --output /path/out \
        --pattern "*.root"

Requirements:
    PyROOT must be available. Activate it with:
        source $(root-config --prefix)/bin/thisroot.sh
    then run this script with the same Python interpreter.

Author: generated for P2-basket analysis
"""

import argparse
import os
import glob
import sys
import time


# ── ROOT import ───────────────────────────────────────────────────────────────

def _import_root():
    try:
        import ROOT
        ROOT.gROOT.SetBatch(True)
        ROOT.gErrorIgnoreLevel = ROOT.kWarning   # silence per-key recovery spam
        return ROOT
    except ModuleNotFoundError:
        print(
            "[ERROR] Cannot import ROOT.\n"
            "  Make sure PyROOT is available:\n"
            "    source $(root-config --prefix)/bin/thisroot.sh\n"
            "  then re-run this script with the same Python interpreter."
        )
        sys.exit(1)


# ── uproot check ─────────────────────────────────────────────────────────────

def _needs_recovery(fpath: str) -> bool:
    """
    Returns True if uproot cannot read the file (keys list is empty).
    Falls back to True if uproot is not installed.
    """
    try:
        import uproot
        with uproot.open(fpath) as f:
            return len(f.keys()) == 0
    except Exception:
        return True


# ── single-file recovery ──────────────────────────────────────────────────────

def recover_file(src: str, dst: str, ROOT, *,
                 dry_run: bool = False,
                 force: bool = False,
                 check_first: bool = True) -> str:
    """
    Recover a single ROOT file.

    Returns one of:
        "ok"               – file was recovered successfully
        "skipped_exists"   – dst already exists and --force not set
        "skipped_readable" – uproot can already read src, no action needed
        "dry_run"          – would have recovered but dry_run=True
        "failed"           – ROOT CloneTree() call failed
        "zombie"           – ROOT could not open src at all
    """
    if not os.path.isfile(src):
        return "failed"

    if os.path.exists(dst) and not force:
        return "skipped_exists"

    if check_first and not _needs_recovery(src):
        return "skipped_readable"

    if dry_run:
        return "dry_run"

    os.makedirs(os.path.dirname(os.path.abspath(dst)), exist_ok=True)

    f_in = ROOT.TFile.Open(src)
    if not f_in or f_in.IsZombie():
        return "zombie"

    # CloneTree rewrites a properly closed file that uproot can read.
    # TFile::Cp just copies the raw bytes and does NOT fix the broken header.
    try:
        tree = f_in.Get("hits")
        if not tree:
            f_in.Close()
            return "failed"
        f_out    = ROOT.TFile(dst, "RECREATE")
        new_tree = tree.CloneTree(-1, "fast")
        new_tree.Write()
        f_out.Close()
        success  = True
    except Exception as e:
        print(f"    [CloneTree error] {e}")
        success = False
    finally:
        f_in.Close()

    return "ok" if success else "failed"


# ── verification ──────────────────────────────────────────────────────────────

def verify_file(dst: str) -> bool:
    """Return True if uproot can open the file and find at least one key."""
    try:
        import uproot
        with uproot.open(dst) as f:
            return len(f.keys()) > 0
    except Exception:
        return False


# ── pair collection ───────────────────────────────────────────────────────────

def _collect_pairs(input_path: str, output_path: str,
                   pattern: str, recursive: bool):
    """Build a list of (src, dst) path pairs."""
    pairs = []

    if os.path.isfile(input_path):
        pairs.append((input_path, output_path))

    elif os.path.isdir(input_path):
        if recursive:
            subdirs = sorted(
                d for d in glob.glob(os.path.join(input_path, "run_*"))
                if os.path.isdir(d)
            )
            if not subdirs:
                subdirs = [input_path]
        else:
            subdirs = [input_path]

        for subdir in subdirs:
            rel = os.path.relpath(subdir, input_path)
            out_subdir = os.path.join(output_path, rel) if rel != "." else output_path
            for fpath in sorted(glob.glob(os.path.join(subdir, pattern))):
                fname = os.path.basename(fpath)
                pairs.append((fpath, os.path.join(out_subdir, fname)))
    else:
        print(f"[ERROR] Input path not found: {input_path}")
        sys.exit(1)

    return pairs


# ── run-level summary ─────────────────────────────────────────────────────────

def _run_info(run_dir: str, pattern: str) -> dict:
    """Gather basic info about a run directory for display."""
    files = sorted(glob.glob(os.path.join(run_dir, pattern)))
    n_total = len(files)
    n_needs = sum(1 for f in files if _needs_recovery(f))
    return {"path": run_dir, "n_total": n_total, "n_needs": n_needs}


# ── interactive run picker ────────────────────────────────────────────────────

def interactive_mode(base_dir: str, output_base: str, pattern: str,
                     args) -> None:
    """
    List all run_* directories under base_dir, show how many files need
    recovery in each, and let the user pick which runs to process.
    """
    run_dirs = sorted(
        d for d in glob.glob(os.path.join(base_dir, "run_*"))
        if os.path.isdir(d)
    )

    if not run_dirs:
        print(f"[ERROR] No run_* directories found under {base_dir}")
        sys.exit(1)

    print(f"\n{'─'*65}")
    print(f"  Base directory : {base_dir}")
    print(f"  Output base    : {output_base}")
    print(f"  File pattern   : {pattern}")
    print(f"{'─'*65}")
    print(f"  Scanning {len(run_dirs)} run directories for files needing recovery...")
    print(f"  (This does a quick uproot check on each file — may take a moment)\n")

    infos = []
    for i, run_dir in enumerate(run_dirs, 1):
        run_name = os.path.basename(run_dir)
        info = _run_info(run_dir, pattern)
        infos.append(info)
        status = f"{info['n_needs']}/{info['n_total']} need recovery" if info["n_total"] > 0 else "no matching files"
        marker = "  !" if info["n_needs"] > 0 else "   "
        print(f"{marker} [{i:3d}] {run_name:<20s}  {status}")

    print(f"\n{'─'*65}")
    print("  Enter run numbers to recover (comma-separated), e.g.: 1,3,5")
    print("  Or enter 'all' to recover all runs with files needing recovery.")
    print("  Or enter 'q' to quit.")
    print(f"{'─'*65}")

    while True:
        try:
            choice = input("\n  Your choice: ").strip().lower()
        except (KeyboardInterrupt, EOFError):
            print("\nAborted.")
            sys.exit(0)

        if choice == "q":
            print("Exiting.")
            sys.exit(0)

        if choice == "all":
            selected = [info for info in infos if info["n_needs"] > 0]
            if not selected:
                print("  No runs need recovery. Nothing to do.")
                sys.exit(0)
            break

        try:
            indices = [int(x.strip()) for x in choice.split(",")]
            selected = []
            valid = True
            for idx in indices:
                if idx < 1 or idx > len(infos):
                    print(f"  [ERROR] Index {idx} out of range (1–{len(infos)}). Try again.")
                    valid = False
                    break
                selected.append(infos[idx - 1])
            if valid:
                break
        except ValueError:
            print("  [ERROR] Invalid input. Enter numbers like: 1,3,5  or  all  or  q")

    print(f"\n  Selected {len(selected)} run(s):")
    for info in selected:
        run_name = os.path.basename(info["path"])
        print(f"    {run_name}  ({info['n_needs']}/{info['n_total']} files need recovery)")

    confirm = input("\n  Proceed? [y/N]: ").strip().lower()
    if confirm != "y":
        print("Aborted.")
        sys.exit(0)

    ROOT = _import_root()

    total_counts = {"ok": 0, "skipped_exists": 0, "skipped_readable": 0,
                    "dry_run": 0, "failed": 0, "zombie": 0}
    t0 = time.time()

    for info in selected:
        run_name = os.path.basename(info["path"])
        out_dir  = os.path.join(output_base, run_name)
        pairs    = _collect_pairs(info["path"], out_dir, pattern, recursive=False)

        print(f"\n  ── {run_name} ({len(pairs)} files) ──")
        counts = _process_pairs(pairs, ROOT, args, base_ref=info["path"])
        for k, v in counts.items():
            total_counts[k] += v

    _print_summary(total_counts, time.time() - t0)


# ── batch processing ──────────────────────────────────────────────────────────

def _process_pairs(pairs, ROOT, args, base_ref: str = "") -> dict:
    """Process a list of (src, dst) pairs and return status counts."""
    counts = {"ok": 0, "skipped_exists": 0, "skipped_readable": 0,
              "dry_run": 0, "failed": 0, "zombie": 0}

    for i, (src, dst) in enumerate(pairs, 1):
        fname   = os.path.basename(src)
        rel_dir = os.path.relpath(os.path.dirname(src), base_ref) if base_ref else ""
        prefix  = f"{rel_dir}/" if rel_dir and rel_dir != "." else ""

        status = recover_file(
            src, dst, ROOT,
            dry_run    = args.dry_run,
            force      = args.force,
            check_first= not args.no_check,
        )
        counts[status] += 1

        icon  = {"ok": "✓", "skipped_exists": "–", "skipped_readable": "✓",
                 "dry_run": "?", "failed": "✗", "zombie": "✗"}.get(status, "?")
        label = {
            "ok":               "recovered",
            "skipped_exists":   "skip (already exists)",
            "skipped_readable": "skip (uproot can read)",
            "dry_run":          "dry-run",
            "failed":           "FAILED (CloneTree error)",
            "zombie":           "FAILED (zombie file)",
        }.get(status, status)

        print(f"    [{i:4d}/{len(pairs)}] {icon}  {prefix}{fname}  →  {label}")

        if status == "ok" and args.verify:
            ok = verify_file(dst)
            if ok:
                print(f"              ✓ uproot verified")
            else:
                print(f"              ✗ uproot STILL cannot read recovered file!")
                counts["ok"]     -= 1
                counts["failed"] += 1

    return counts


def _print_summary(counts: dict, elapsed: float) -> None:
    print(f"\n{'─'*60}")
    print(f"  Done in {elapsed:.1f}s")
    print(f"  Recovered  : {counts['ok']}")
    print(f"  Skipped    : {counts['skipped_exists'] + counts['skipped_readable']}")
    print(f"  Dry-run    : {counts['dry_run']}")
    print(f"  Failed     : {counts['failed'] + counts['zombie']}")
    if counts["failed"] + counts["zombie"] > 0:
        print("  [WARNING] Some files could not be recovered — see above.")


# ── non-interactive batch run ─────────────────────────────────────────────────

def batch_mode(args) -> None:
    ROOT = _import_root()
    pairs = _collect_pairs(args.input, args.output, args.pattern, args.recursive)

    if not pairs:
        print("[WARNING] No files matched. Check --input and --pattern.")
        return

    prefix = "[DRY RUN] " if args.dry_run else ""
    print(f"{prefix}Found {len(pairs)} file(s) to process.\n")

    t0     = time.time()
    counts = _process_pairs(pairs, ROOT, args, base_ref=args.input)
    _print_summary(counts, time.time() - t0)


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Recover unclosed ROOT files so uproot can read them.\n\n"
            "Run with --base to enter interactive mode and pick individual runs.\n"
            "Run with --input for non-interactive (single file, single run, or all runs)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # ── mode selection ──
    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument(
        "--base", "-b",
        metavar="BASE_DIR",
        help=(
            "Interactive mode: base directory containing run_* subdirs. "
            "Lists all runs and lets you pick which to recover."
        ),
    )
    mode_group.add_argument(
        "--input", "-i",
        metavar="PATH",
        help=(
            "Non-interactive mode: single .root file, a run_* directory, "
            "or (with --recursive) a base directory of run_* dirs."
        ),
    )

    parser.add_argument(
        "--output", "-o", required=True,
        help=(
            "Output path. For interactive/recursive mode: base output directory "
            "(run_* structure is mirrored). For single run/file: direct output path."
        ),
    )
    parser.add_argument(
        "--pattern", "-p", default="enp*.root",
        help="Glob pattern for ROOT files inside directories (default: enp*.root).",
    )
    parser.add_argument(
        "--recursive", "-r", action="store_true",
        help="(--input mode only) Recurse into all run_* subdirectories at once.",
    )
    parser.add_argument(
        "--force", "-f", action="store_true",
        help="Overwrite existing output files.",
    )
    parser.add_argument(
        "--dry-run", "-n", action="store_true",
        help="Report what would be done without writing anything.",
    )
    parser.add_argument(
        "--verify", "-v", action="store_true",
        help="After recovery, verify each output file is readable by uproot.",
    )
    parser.add_argument(
        "--no-check", action="store_true",
        help="Skip uproot pre-check and recover all files unconditionally.",
    )

    args = parser.parse_args()

    if args.base:
        interactive_mode(args.base, args.output, args.pattern, args)
    else:
        batch_mode(args)


if __name__ == "__main__":
    # ──────────────────────────────────────────────────────────────────────
    # CONFIGURE HERE — edit these before hitting Run in PyCharm
    # ──────────────────────────────────────────────────────────────────────

    # Mode: "interactive", "single_run", "single_file", or "all_runs"
    RUN_MODE = "single_run"

    # Base directory containing all run_* subdirectories
    BASE_DIR = "/drf/projets/clas12/cern_202511_p2_alinx"

    # Where recovered files will be written (run_* structure is mirrored)
    OUTPUT_BASE = "/drf/projets/clas12/cern_202511_p2_alinx_recovered"

    # --- Only used in "single_run" mode ---
    # Change this to whichever run you want to recover
    RUN_NAME = "run_101"

    # --- Only used in "single_file" mode ---
    SINGLE_FILE_IN  = f"{BASE_DIR}/{RUN_NAME}/enp4s0f1_20251115-215717_0_20251115220717.root"
    SINGLE_FILE_OUT = f"{OUTPUT_BASE}/{RUN_NAME}/enp4s0f1_20251115-215717_0_20251115220717.root"

    # Options (True/False)
    VERIFY   = True    # check each recovered file with uproot afterwards
    FORCE    = False   # overwrite already-recovered files
    DRY_RUN  = False   # set True to preview without writing anything
    NO_CHECK = False   # set True to skip uproot pre-check
    PATTERN  = "enp*.root"  # file glob pattern inside run directories

    # ──────────────────────────────────────────────────────────────────────
    # Nothing to edit below this line
    # ──────────────────────────────────────────────────────────────────────

    class _Args:
        recursive = False

    args          = _Args()
    args.pattern  = PATTERN
    args.force    = FORCE
    args.dry_run  = DRY_RUN
    args.verify   = VERIFY
    args.no_check = NO_CHECK

    if RUN_MODE == "interactive":
        interactive_mode(BASE_DIR, OUTPUT_BASE, PATTERN, args)

    elif RUN_MODE == "single_run":
        args.input     = os.path.join(BASE_DIR, RUN_NAME)
        args.output    = os.path.join(OUTPUT_BASE, RUN_NAME)
        args.recursive = False
        batch_mode(args)

    elif RUN_MODE == "single_file":
        args.input     = SINGLE_FILE_IN
        args.output    = SINGLE_FILE_OUT
        args.recursive = False
        batch_mode(args)

    elif RUN_MODE == "all_runs":
        args.input     = BASE_DIR
        args.output    = OUTPUT_BASE
        args.recursive = True
        batch_mode(args)

    else:
        print(f"[ERROR] Unknown RUN_MODE: '{RUN_MODE}'. "
              f"Choose: interactive, single_run, single_file, all_runs")

print('bonzo')


if __name__ == '__main__':
    main()
