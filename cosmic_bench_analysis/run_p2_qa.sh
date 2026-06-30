#!/bin/bash
# run_p2_qa.sh
#
# Convenience wrapper for p2_channel_qa.py — the raw P2 channel QA on a cosmic
# bench run (pulse height, multiplicity, saturation, rate vs time, per-FEU
# occupancy/pulse-height, and the reference-coincidence "cosmic-tagged" view).
#
# Usage:
#   ./run_p2_qa.sh                        # uses the default run below
#   ./run_p2_qa.sh /path/to/combined_hits_root   # any other run's hits dir
#   ./run_p2_qa.sh /path/to/combined_hits_root  --coinc-window-ns 80
#
# Extra args after the hits dir are passed straight through to p2_channel_qa.py
# (e.g. --ref-det-type m3, --det-name, --out-dir, --coinc-window-ns).
#
# The script auto-finds run_config.json two dirs above the hits dir, and reads
# the P2 / reference FEU assignment from it, so nothing here is hard-coded per
# detector.

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# Default run: the long P2 det1 sanity-check overnight run (6-27-26).
DEFAULT_HITS_DIR="/local/home/ak271430/Documents/PostDocSaclay/data/Cosmic_Bench/mx17_det3_p2_det1_overnight_6-27-26/long_run_p2_det1_sanity_check/combined_hits_root"

# First positional arg overrides the hits dir; remaining args pass through.
if [[ $# -ge 1 && "$1" != -* ]]; then
    HITS_DIR="$1"
    shift
else
    HITS_DIR="$DEFAULT_HITS_DIR"
fi

# Use the project venv if present (uproot, pandas, matplotlib live there).
if [[ -f "${SCRIPT_DIR}/../.venv/bin/activate" ]]; then
    # shellcheck disable=SC1091
    source "${SCRIPT_DIR}/../.venv/bin/activate"
fi

echo "=== P2 channel QA ==="
echo "hits dir : ${HITS_DIR}"
echo "extra    : $*"
echo ""

python3 "${SCRIPT_DIR}/p2_channel_qa.py" --hits-dir "${HITS_DIR}" "$@"
