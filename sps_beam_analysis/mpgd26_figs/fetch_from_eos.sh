#!/bin/bash
# Populate the MPGD26 figure workspace from the nTOF EOS copy of the campaign.
#
# The stage products on EOS are the newest ones (the LaCie copy predates the
# condor re-processing and, e.g., still carries the old stage-29 column names),
# so EOS is the source of truth here.  Only small files travel: CSV/JSON stage
# products, run configs and the VMM scalars.  The HV monitor traces are 136 MB
# and only feed the campaign-timeline figure, so they are opt-in.
#
#   ./fetch_from_eos.sh            products + configs (~50 MB)
#   ./fetch_from_eos.sh --hv       also the HV monitor CSVs (~136 MB)
#
set -euo pipefail

W=${MPGD26_WORKSPACE:-$HOME/Documents/PostDocSaclay/data/SPS_Beam_Test/mpgd26_workspace}
LX=${LXPLUS_HOST:-lxplus}
E=/eos/experiment/ntof/data/x17/p2_sps_july
WANT_HV=0
[ "${1:-}" = "--hv" ] && WANT_HV=1

say() { echo "[fetch] $*"; }

# products/analysis is a real directory here (it was a LaCie symlink while the
# workspace was being set up -- replace it, keeping nothing stale).
if [ -L "$W/products/analysis" ]; then rm "$W/products/analysis"; fi
mkdir -p "$W/products/analysis" "$W/eos_inventory" \
         "$W/products/urw_timebins_x" "$W/products/hv"

say 'stage products (CSV/JSON) from analysis/'
ssh "$LX" "cd $E && tar -czf - \
      \$(find analysis -type f \\( -name '*.csv' -o -name '*.json' \\) \
         -not -path '*.sys.v#*' -printf '%p\n')" \
  | tar -xzf - -C "$W/products" --strip-components=1 --one-top-level=analysis

say 'uRWELL time-binned efficiency products'
ssh "$LX" "cd $E/analysis && tar -czf - urw_timebins 2>/dev/null" \
  | tar -xzf - -C "$W/products/urw_timebins_x"

say 'VMM per-capture scalars'
ssh "$LX" "cd $E && tar -czf - \
      \$(find vmm/runs -name 'scalars.json' -printf '%p\n' 2>/dev/null)" \
  | tar -xzf - -C "$W/products/analysis" --transform 's|^vmm/runs|vmm|'

say 'run configs (DREAM DAQ and VMM DAQ number their runs independently)'
ssh "$LX" "cd $E && tar -czf - \
      \$(find runs -maxdepth 2 -name 'run_config.json' -printf '%p\n') \
      \$(find vmm/runs -maxdepth 2 -name 'run_config.json' -printf '%p\n' \
         2>/dev/null)" | tar -xzf - -C "$W/eos_inventory"
mkdir -p "$W/eos_inventory/vmm_meta"
[ -d "$W/eos_inventory/vmm/runs" ] && \
  rm -rf "$W/eos_inventory/vmm_meta/runs" && \
  mv "$W/eos_inventory/vmm/runs" "$W/eos_inventory/vmm_meta/runs"

if [ "$WANT_HV" = 1 ]; then
  say 'HV monitor traces (136 MB)'
  ssh "$LX" "cd $E && tar -czf - \
        \$(find runs -name 'hv_monitor.csv' -printf '%p\n')" \
    | tar -xzf - -C "$W/products/hv" --strip-components=1
fi

say 'HV setpoints table'
python3 "$(dirname "$0")/make_hv_setpoints.py"
say "done -> $W"
