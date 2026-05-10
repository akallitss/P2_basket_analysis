"""
watch_baking_buckle.py  -  Live monitor for a running baking buckling sweep.

Run in a second terminal (plain Python 3, no FreeCAD needed) while
FreeCADCmd executes freecad_baking_buckling.py in another window:

    python FreeCAD/watch_baking_buckle.py

Polls baking_buckle.log every 2 s, printing new lines as they arrive.
Prints a summary table each time a tension case finishes.
Exits when the run completes (all tensions done or hard failure).
"""

import os
import sys
import time
import re

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
LOG_FILE    = os.path.join(SCRIPT_DIR, "baking_buckle.log")
RESULTS_DIR = os.path.join(SCRIPT_DIR, "results")
POLL_S      = 2.0   # seconds between polls

# Regex to pull eigenvalue lines out of the log
RE_MODE = re.compile(
    r"Mode\s+(\d+):\s+lambda=([-\d.]+)\s+T_buckle~([\d.]+)\s+C\s+\[(\w+)\]"
)
RE_TENSION = re.compile(r"--- T = (\d+) N/m")
RE_DONE    = re.compile(r"(All \d+ tensions completed|FAILED tensions:)")


def clear_line():
    sys.stdout.write("\r\033[K")


def print_summary(results):
    """results: list of dicts {T_Nm, mode, lam, T_buckle, status}"""
    if not results:
        return
    print()
    print("  {:<10} {:>6} {:>10} {:>12} {:>10}".format(
        "T (N/m)", "Mode", "lambda", "T_buckle (C)", "Status"))
    print("  " + "-" * 52)
    for r in results:
        print("  {:<10} {:>6} {:>10.4f} {:>12.0f} {:>10}".format(
            r["T_Nm"], r["mode"], r["lam"], r["T_buckle"], r["status"]))
    print()


def main():
    print("Watching: {}".format(LOG_FILE))
    print("Press Ctrl-C to stop.\n")

    # Wait up to 60 s for the log to appear
    waited = 0
    while not os.path.isfile(LOG_FILE):
        sys.stdout.write("\rWaiting for run to start ... {}s".format(waited))
        sys.stdout.flush()
        time.sleep(1.0)
        waited += 1
        if waited % 60 == 0:
            print("\r  (still waiting — launch FreeCADCmd when ready)       ")

    if waited:
        print()

    results = []
    current_T = None
    pos = 0

    try:
        while True:
            with open(LOG_FILE, "r", encoding="utf-8", errors="replace") as f:
                f.seek(pos)
                new_text = f.read()
                pos = f.tell()

            if new_text:
                for line in new_text.splitlines():
                    print(line)

                    m = RE_TENSION.search(line)
                    if m:
                        current_T = int(m.group(1))

                    m = RE_MODE.search(line)
                    if m and current_T is not None:
                        results.append({
                            "T_Nm":    current_T,
                            "mode":    int(m.group(1)),
                            "lam":     float(m.group(2)),
                            "T_buckle": float(m.group(3)),
                            "status":  m.group(4),
                        })

                    if RE_DONE.search(line):
                        # Final summary table
                        print()
                        print("=" * 56)
                        print("  SWEEP COMPLETE — eigenvalue summary")
                        print_summary(results)
                        # Check for summary CSV
                        scsv = os.path.join(RESULTS_DIR, "baking_buckle_summary.csv")
                        if os.path.isfile(scsv):
                            print("  Summary CSV: " + scsv)
                        print("=" * 56)
                        return

            time.sleep(POLL_S)

    except KeyboardInterrupt:
        print("\n\nInterrupted.")
        if results:
            print("Partial results so far:")
            print_summary(results)


main()