#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gases.py -- single source of truth for the P2 gas studies.

Each entry maps a short key (used for file names) to:
  label     human-readable name for plots/tables
  magboltz  list of (Garfield/Magboltz gas name, percentage) pairs
  role      "reference" (the gas we are calibrated on) or "candidate"
            (a gas we want to map onto the reference); free text is fine.

To add a NEW gas at the beam, just append an entry here and rerun
    ./run_lxplus.sh <newkey>
then
    python3 analyze.py
Nothing else needs to change -- the C++ scanner reads the composition from
the command line and the analysis reads whatever CSVs exist in results/.

Garfield/Magboltz gas names are lower-case tags, e.g.:
  ar (argon), co2, ch4 (methane), ic4h10 (isobutane), cf4, c2h6, n2, xe, ...
See the Magboltz gas list in the Garfield++ docs for the full set.
"""

# --- detector geometry (fixed by the P2 Micromegas) ------------------------
D_AMP_CM = 0.015    # amplification gap  = 150 um
D_DRIFT_CM = 0.300  # drift/conversion gap = 3 mm

# --- gas registry ----------------------------------------------------------
GASES = {
    "ar_iso_95_5": {
        "label": "Ar/iC4H10 95/5",
        "magboltz": [("ar", 95.0), ("ic4h10", 5.0)],
        "role": "reference",   # Saclay lab calibration gas
    },
    "ar_co2_iso_95_3_2": {
        "label": "Ar/CO2/iC4H10 95/3/2",
        "magboltz": [("ar", 95.0), ("co2", 3.0), ("ic4h10", 2.0)],
        "role": "candidate",   # NSW "NSW gas", SPS beam test
    },
    # --- append new beam gases below, e.g. ---------------------------------
    # "ar_co2_93_7": {
    #     "label": "Ar/CO2 93/7",
    #     "magboltz": [("ar", 93.0), ("co2", 7.0)],
    #     "role": "candidate",
    # },
}


def magboltz_args(key):
    """Return the argv tail for p2_gas_scan, e.g. ['ar', '95', 'ic4h10', '5']."""
    comp = GASES[key]["magboltz"]
    out = []
    for name, frac in comp:
        out += [name, ("%g" % frac)]
    return out


if __name__ == "__main__":
    import sys
    # `python3 gases.py --emit-args [key ...]` prints one line per gas:
    #   <key> <gasA> <fracA> [<gasB> <fracB> ...]
    # consumed by run_lxplus.sh (single source of truth for compositions).
    if len(sys.argv) >= 2 and sys.argv[1] == "--emit-args":
        keys = sys.argv[2:] or list(GASES.keys())
        for k in keys:
            if k not in GASES:
                sys.stderr.write("unknown gas key: %s\n" % k)
                sys.exit(2)
            print(k + " " + " ".join(magboltz_args(k)))
    else:
        for k, v in GASES.items():
            print("%-22s %-24s %s" % (k, v["label"], v["role"]))
