#!/bin/bash
# run_m3_pipeline.sh
#
# Runs the M3 tracking DataReader pipeline on both the non-ZS and ZS runs,
# producing analyse trees (cluster data) that can be compared by
# m3_signal_diagnostic.py.
#
# Pipeline:
#   non-ZS FDF  → read → signal tree
#   signal tree → ped  → Ped.dat + RMSPed.dat
#   non-ZS FDF + Ped.dat → analyse → analyse_non_zs.root
#   ZS FDF     + Ped.dat → analyse → analyse_zs.root

set -e

TRACKING_DIR="$(cd "$(dirname "$0")/../../cosmic_bench_m3_tracking" && pwd)"
DATA_BASE="$(cd "$(dirname "$0")/../../data/Cosmic_Bench/mx17_non_and_ZS_Test_6-17-26" && pwd)"

NON_ZS_BASENAME="${DATA_BASE}/run_non_zs/raw_daq_data/MX17_run_non_zs_datrun_260617_09H21_"
ZS_BASENAME="${DATA_BASE}/run_zs/raw_daq_data/MX17_run_zs_datrun_260617_09H32_"

cd "$TRACKING_DIR"
mkdir -p root_files

source /local/home/ak271430/Software/root/bin/thisroot.sh

# ---------------------------------------------------------------------------
# Build config for non-ZS (read + ped + analyse)
# Modifies Tree and signal_file paths; Ped/RMSPed stay at default locations.
# ---------------------------------------------------------------------------
# Strip // comments — handles quotes inside comments correctly (Boost read_json rejects them)
python3 - << 'PYEOF'
import sys
text = open('config_ref_2022.json').read()
in_str, result, i = False, [], 0
while i < len(text):
    c = text[i]
    if c == '"' and (i == 0 or text[i-1] != '\\'):
        in_str = not in_str
    if not in_str and c == '/' and i+1 < len(text) and text[i+1] == '/':
        while i < len(text) and text[i] != '\n':
            i += 1
        continue
    result.append(c)
    i += 1
open('/tmp/config_stripped.json', 'w').write(''.join(result))
PYEOF

sed 's|root_files/test_analyse.root|root_files/analyse_non_zs.root|g;
     s|root_files/test_signal.root|root_files/signal_non_zs.root|g' \
    /tmp/config_stripped.json > /tmp/config_base_non_zs.json

{
    cat /tmp/config_base_non_zs.json
    printf '    "data_file_first" : 0,\n'
    printf '    "data_file_last"  : 0,\n'
    printf '    "data_file_basename" : "%s"\n' "${NON_ZS_BASENAME}"
    printf '}\n'
} > config_diag_non_zs.json

# ---------------------------------------------------------------------------
# Build config for ZS (analyse only — reuses Ped.dat from non-ZS)
# ---------------------------------------------------------------------------
sed 's|root_files/test_analyse.root|root_files/analyse_zs.root|g;
     s|root_files/test_signal.root|root_files/signal_zs.root|g' \
    /tmp/config_stripped.json > /tmp/config_base_zs.json

{
    cat /tmp/config_base_zs.json
    printf '    "data_file_first" : 0,\n'
    printf '    "data_file_last"  : 0,\n'
    printf '    "data_file_basename" : "%s"\n' "${ZS_BASENAME}"
    printf '}\n'
} > config_diag_zs.json

# ---------------------------------------------------------------------------
# Run pipeline
# ---------------------------------------------------------------------------
echo ""
echo "=== [1/4] Read non-ZS FDF → raw signal tree ==="
./DataReader config_diag_non_zs.json read

echo ""
echo "=== [2/4] Compute pedestals from non-ZS signal tree ==="
./DataReader config_diag_non_zs.json ped

echo ""
echo "=== [3/4] Build non-ZS analyse tree (cluster data) ==="
./DataReader config_diag_non_zs.json analyse

echo ""
echo "=== [4/4] Build ZS analyse tree (using non-ZS pedestals) ==="
./DataReader config_diag_zs.json analyse

echo ""
echo "=== Done ==="
ls -lh root_files/analyse_non_zs.root root_files/analyse_zs.root \
        root_files/test_Ped.dat root_files/test_RMSPed.dat 2>/dev/null
