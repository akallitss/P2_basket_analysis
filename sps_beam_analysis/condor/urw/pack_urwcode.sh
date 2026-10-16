#!/usr/bin/env bash
# pack_urwcode.sh -- build urwcode.tgz for the uRWELL condor pipeline, with
# its version stamp. Until 2026-10-15 this tarball was packed by hand; the
# 2026-08-16 one on lxplus is identical to tag mpgd26-urw + saclay_micromegas_EIC
# 2eaf34e.
#
#   pack_urwcode.sh [saclay_micromegas_EIC dir]
#
# Then copy urwcode.tgz to lxplus ~/p2_condor_urw/ before condor_submit.
# urwtb_job.sh uploads CODE_VERSION.txt next to the products.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$HERE/../../.." && pwd)                       # .../P2_basket_analysis
SMM=$(cd "${1:-$REPO/../saclay_micromegas_EIC}" && pwd)
P2_PATHS=(sps_beam_analysis cosmic_bench_analysis
          Detector_Mapping/P2_BASKET/P2_BASKET_mapping.csv)
# the package's source, as in the 2026-08-16 tarball: top-level files plus
# Detector_Classes/ and quick_scripts/; no notebooks, poc/, tests/ or binaries
mapfile -t SMM_FILES < <(cd "$SMM" && git ls-files -- '*.py' '*.cxx' '*.h' \
                         Makefile README.md requirements.txt \
                         | grep -E '^([^/]+|Detector_Classes/[^/]+|quick_scripts/[^/]+)$')

STAMP="$HERE/stamp"
rm -rf "$STAMP"; mkdir -p "$STAMP"
{ "$HERE/../code_version.sh" "$REPO" "$STAMP/CODE_DIFF.patch" "${P2_PATHS[@]}"
  echo
  "$HERE/../code_version.sh" "$SMM" "$STAMP/CODE_DIFF_saclay_micromegas.patch" \
      "${SMM_FILES[@]}"
} > "$STAMP/CODE_VERSION.txt"
cat "$STAMP/CODE_VERSION.txt"

SMM_NAME=$(basename "$SMM")
tar --exclude=__pycache__ --exclude='*.pyc' --exclude=logbook \
    --exclude=patches --exclude=condor --exclude=mpgd26_figs \
    --exclude=vmm_dream_matching --exclude=dream_vmm_comparison \
    -czf "$HERE/urwcode.tgz" \
    -C "$REPO" "${P2_PATHS[@]}" \
    -C "$(dirname "$SMM")" "${SMM_FILES[@]/#/$SMM_NAME/}" \
    -C "$STAMP" $(ls "$STAMP")
echo "urwcode.tgz  $(du -h "$HERE/urwcode.tgz" | cut -f1)"
