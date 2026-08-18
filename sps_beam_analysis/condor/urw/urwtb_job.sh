#!/bin/bash
# urwtb_job.sh -- one condor job = uRWELL-referenced P2 efficiency with time
# bins for one (run, sub_run). Stages combined_hits + hv_monitor +
# recorded_events + run_config + config/detectors from nTOF EOS, runs the
# patched urw_p2_efficiency.py --time-bins, pushes products to
#   analysis/urw_timebins/<run>/<sub_run>/
set -uo pipefail
RUN=${1:?}; SUB=${2:?}; TBINS=${3:-12}
EOS_URL=root://eospublic.cern.ch
EBASE=/eos/experiment/ntof/data/x17/p2_sps_july
SCRATCH=${_CONDOR_SCRATCH_DIR:-$PWD}
say() { echo "[$(date -u +%H:%M:%S)] $*"; }
fail() { echo "FATAL: $*" >&2; exit 1; }
UPLOAD_OK=0
pack_products() {
  if [ "$UPLOAD_OK" = "1" ]; then tar -czf "$SCRATCH/products.tgz" -T /dev/null 2>/dev/null; return; fi
  tar -czf "$SCRATCH/products.tgz" -C "$SCRATCH/out" . 2>/dev/null || tar -czf "$SCRATCH/products.tgz" -T /dev/null 2>/dev/null
}
mkdir -p "$SCRATCH/out"; pack_products; trap pack_products EXIT
say "=== urwtb  $RUN/$SUB  time_bins=$TBINS  host=$(hostname -f)"
if [ -n "${_CONDOR_CREDS:-}" ] && [ -d "${_CONDOR_CREDS}" ]; then
  cc=$(ls "$_CONDOR_CREDS"/*.cc 2>/dev/null | head -1); [ -n "$cc" ] && export KRB5CCNAME="FILE:$cc"
fi
LCG_VIEW=${LCG_VIEW:-/cvmfs/sft.cern.ch/lcg/views/LCG_110}
PLAT=""; for p in x86_64-el9-gcc13-opt x86_64-el9-gcc14-opt; do [ -f "$LCG_VIEW/$p/setup.sh" ] && { PLAT=$p; break; }; done
[ -n "$PLAT" ] || fail "no LCG platform"
set +u; source "$LCG_VIEW/$PLAT/setup.sh" || fail "LCG setup"; set -u
export MPLBACKEND=Agg MPLCONFIGDIR="$SCRATCH/.mpl" HOME="$SCRATCH"; mkdir -p "$MPLCONFIGDIR"
tar -xzf urwcode.tgz -C "$SCRATCH" || fail "unpack"
export SACLAY_MM_DIR="$SCRATCH/saclay_micromegas_EIC"

D=$SCRATCH/data/$RUN/$SUB; mkdir -p "$D"
xrdcp -f -s "$EOS_URL/$EBASE/runs/$RUN/run_config.json" "$SCRATCH/data/$RUN/run_config.json" || fail "run_config"
# The DAQ wrote det_type "P2" for the two uRWELL reference planes in
# p2in_hvrange_2 -- the only run in the campaign with this. Correct the STAGED
# copy so tracking can build the uRWELL geometry; EOS is left untouched.
python3 - "$SCRATCH/data/$RUN/run_config.json" <<'PYFIX'
import json, sys
p = sys.argv[1]
d = json.load(open(p))
want = {"EIC_uRWELL_front": "urw_inter", "EIC_uRWELL_back": "urw_strip"}
n = 0
for det in d.get("detectors", []):
    t = want.get(det.get("name"))
    if t and det.get("det_type") != t:
        det["det_type"] = t
        n += 1
if n:
    json.dump(d, open(p, "w"), indent=1)
    print(f"[fix] corrected {n} mis-typed uRWELL det_type entries")
PYFIX

for f in hv_monitor.csv recorded_events.npz run_time.txt; do
  xrdcp -f -s "$EOS_URL/$EBASE/runs/$RUN/$SUB/$f" "$D/$f" 2>/dev/null || say "note: no $f"
done
mkdir -p "$D/combined_hits_root"
for f in $(xrdfs $EOS_URL ls "$EBASE/runs/$RUN/$SUB/combined_hits_root" 2>/dev/null | grep '\.root$' | grep -v '.sys.v#'); do
  xrdcp -f -s "$EOS_URL/$f" "$D/combined_hits_root/" || fail "stage $f"
done
[ -n "$(ls -A "$D/combined_hits_root")" ] || fail "no combined hits"
# hits_root for the uRWELL FEU (feu 01) if present -- faster + safer than fallback
mkdir -p "$D/hits_root"
for f in $(xrdfs $EOS_URL ls "$EBASE/runs/$RUN/$SUB/hits_root" 2>/dev/null | grep '_01_hits\.root$' | grep -v '.sys.v#'); do
  xrdcp -f -s "$EOS_URL/$f" "$D/hits_root/" || say "note: hits_root stage failed"
done
mkdir -p "$SCRATCH/config/detectors"
for f in $(xrdfs $EOS_URL ls "$EBASE/config/detectors" 2>/dev/null | grep -v '.sys.v#'); do
  xrdcp -f -s "$EOS_URL/$f" "$SCRATCH/config/detectors/" || true
done
export URW_DET_DIR="$SCRATCH/config/detectors/"
df -h "$SCRATCH" | tail -1

cd "$SCRATCH/sps_beam_analysis/urw_reference"
python3 urw_p2_efficiency.py --run "$RUN" --data-root "$SCRATCH/data" \
    --sub-run "$SUB" --time-bins "$TBINS" --out "$SCRATCH/out" || fail "stage failed (exit 2)" 
OUT_BASE=$EBASE/analysis/urw_timebins/$RUN/$SUB
NUP=0
for f in "$SCRATCH"/out/*; do
  [ -f "$f" ] || continue
  xrdcp -f -s "$f" "$EOS_URL/$OUT_BASE/$(basename "$f")" && NUP=$((NUP+1)) || say "upload failed: $f"
done
say "uploaded $NUP products"
[ "$NUP" -gt 0 ] && UPLOAD_OK=1 || exit 1
say "=== done"
