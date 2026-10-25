#!/bin/bash
# Rebuilds (in place: it replaces the style and the note-body of an existing file) the unlisted CERN-site note "Measuring the VMM noise threshold with cosmic muons" from the parts in this
# folder (concatenated in this order: they share one namespace) and the analysis outputs under
# ~/Documents/PostDocSaclay/data/Lab_Measurements/VMM-cosmics/vmm_threshold_scan/ (the LOCAL path in head2.py).
# Then: cd ~/Documents/PostDocSaclay/cern-site && python3 scripts/build_hub.py && git commit && scripts/deploy.sh
# BODY in part2b.py is an f-string: braces in its HTML must be doubled, but the CSS string above it is not (its }} are real).
set -e
cd "$(dirname "$0")"
NOTE=${1:-$HOME/Documents/PostDocSaclay/cern-site/notes/unlisted/2026-10-22-measuring-the-vmm-noise-threshold-with-cosmic-muons-what-we-did-and-what-we-saw.html}
cat head2.py part2a.py part3.py part4.py part5.py part2b.py > /tmp/build_vmm_cosmic_note.py
python3 /tmp/build_vmm_cosmic_note.py "$NOTE"
