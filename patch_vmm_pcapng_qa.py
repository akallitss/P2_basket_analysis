#!/usr/bin/env python3
"""Generate the patched vmm_pcapng_qa.py: fast decoder + fast histograms.

Two independent changes, both reversible:

  1. Parsing goes through vmm_decode.decode() instead of the in-file scapy /
     parse_block loop. The legacy loop is kept and is used automatically if the
     decoder is missing or raises, and on demand via --legacy-parser.
  2. Every ax.hist() call gets histtype=HIST_TYPE ('step'). matplotlib's
     default 'bar' builds one Rectangle per bin -- ~20k patches per 20-VMM
     figure -- which is the single largest cost in the script.

Neither change alters a plotted value: the decoder is validated bit-identical
to the legacy parser, and histtype only changes how the same bin contents are
drawn.

Usage:
    python patch_vmm_pcapng_qa.py                  # write patched copy + diff
    python patch_vmm_pcapng_qa.py --hist stepfilled
"""

import argparse
import difflib
import os
import re
import sys

SRC = "/local/p2/DAQ_Control_VMM_Beam/vmm_qa/vmm_pcapng_qa.py"

# --- 1. new CLI flag -------------------------------------------------------
ANCHOR_ARG = '_args = _ap.parse_args()'
NEW_ARG = '''_ap.add_argument("--legacy-parser", action="store_true",
                 help="Use the original in-file scapy/parse_block loop instead of "
                      "vmm_decode. Slower; kept as a fallback and for A/B checks.")
_args = _ap.parse_args()'''

# --- 2. histtype constant --------------------------------------------------
ANCHOR_CONST = "TDC_RANGE       = 255   # TDC full-scale bin count (matches vmm-sdat SRSTime::tdc_range)"
NEW_CONST = ANCHOR_CONST + '''

# matplotlib's default histtype='bar' emits one Rectangle patch per bin, so a
# 20-VMM x 1024-bin figure builds ~20k artists and takes ~19 s. 'step' draws
# the identical bin contents as a single line in ~2.5 s.
#
# 'stepfilled' keeps the filled look and is still ~6x faster than 'bar', but
# measured end-to-end it lands at 44.3-44.8 s against the 44.4 s dumpcap
# rotation, i.e. no margin -- so 'step' is the default. Bin contents are the
# same for all three; only the rendering differs.
HIST_TYPE = "@HIST@"'''

# --- 3. the decoder splice -------------------------------------------------
SPLICE = '''# ---- Parse ---------------------------------------------------------------
# Fast path: vmm_decode is a vectorised decoder (no scapy, no per-word Python
# loop) validated to reproduce this script's hits DataFrame bit for bit --
# every column, every dtype -- across SRS and TRG on the beam campaign data.
# It is ~5-12x faster and uses less memory.
#
# The legacy loop below is kept deliberately: if vmm_decode is missing or
# raises for any reason, QA degrades to its previous speed instead of failing.
# --legacy-parser forces it, which is also how to A/B the two.
hits = None
if not _args.legacy_parser:
    try:
        import vmm_decode as _vd
        _t_parse = time.time()
        print(f"Reading: {pcap_file}  (vmm_decode)")
        hits, _dmeta = _vd.decode(
            pcap_file,
            data_format=data_format,
            src_ips=src_ips,
            max_packets=_args.max_packets,
            calibration=((_cal_to, _cal_ts, _cal_ao, _cal_as)
                         if _has_calibration else None),
        )
        pkt_count = _dmeta["n_packets"]
        vm3_count = _dmeta["n_vm3_packets"]
        print(f"Done: {pkt_count} packets | {vm3_count} VM3 packets | "
              f"{len(hits):,} total hits  ({time.time() - _t_parse:.1f}s)")
    except Exception as _e:
        print(f"WARNING: vmm_decode failed ({type(_e).__name__}: {_e}) -- "
              f"falling back to the legacy parser.")
        hits = None

if hits is None:
'''


def build(hist_type="step"):
    src = open(SRC).read()
    lines = src.split("\n")

    # locate the legacy parse region: marker dict .. end of DataFrame build
    i0 = next(i for i, l in enumerate(lines) if l.startswith("markers = {}"))
    i1 = next(i for i, l in enumerate(lines) if l.startswith("del (fec_buf,"))
    i2 = next(i for i, l in enumerate(lines) if i > i1 and l.startswith("mem_mb"))

    legacy = lines[i0:i2]
    # indent the legacy block so it sits under `if hits is None:`
    legacy = [("    " + l if l.strip() else l) for l in legacy]

    patched = lines[:i0] + SPLICE.split("\n") + legacy + lines[i2:]
    out = "\n".join(patched)

    # CLI flag
    assert out.count(ANCHOR_ARG) == 1, "argparse anchor not unique"
    out = out.replace(ANCHOR_ARG, NEW_ARG)

    # histtype constant
    assert out.count(ANCHOR_CONST) == 1, "constant anchor not unique"
    out = out.replace(ANCHOR_CONST, NEW_CONST.replace("@HIST@", hist_type))

    out, n_hist = _patch_hist_calls(out)
    return out, n_hist


def _patch_hist_calls(text):
    """Add histtype=HIST_TYPE to every .hist( call, single- or multi-line.

    Walks paren depth from the `.hist(` to the line that closes it, so calls
    whose arguments span several lines get the kwarg on the right line.
    """
    lines = text.split("\n")
    n = 0
    i = 0
    while i < len(lines):
        if re.search(r"\.hist\(", lines[i]):
            depth = lines[i].count("(") - lines[i].count(")")
            j = i
            while depth > 0 and j + 1 < len(lines):
                j += 1
                depth += lines[j].count("(") - lines[j].count(")")
            block = "\n".join(lines[i:j + 1])
            if "histtype" not in block:
                close = lines[j].rstrip()
                if close.endswith(")"):
                    lines[j] = close[:-1] + ", histtype=HIST_TYPE)"
                    n += 1
            i = j
        i += 1
    return "\n".join(lines), n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hist", default="step",
                    choices=["stepfilled", "step", "bar"])
    ap.add_argument("--out", default=None, help="write patched file here")
    args = ap.parse_args()

    out_path = args.out or os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "vmm_pcapng_qa.patched.py")
    patched, n_hist = build(args.hist)
    open(out_path, "w").write(patched)

    src = open(SRC).read()
    diff = list(difflib.unified_diff(
        src.split("\n"), patched.split("\n"),
        fromfile="vmm_pcapng_qa.py (live)", tofile="vmm_pcapng_qa.py (patched)",
        lineterm="", n=2))
    diff_path = out_path.replace(".py", ".diff")
    open(diff_path, "w").write("\n".join(diff))

    n_all = len(re.findall(r"\.hist\(", patched))
    n_ht = patched.count("histtype=HIST_TYPE")
    print(f"patched   : {out_path}")
    print(f"diff      : {diff_path}  ({len(diff)} lines)")
    print(f"hist calls: {n_all} total, {n_ht} now carry histtype=HIST_TYPE")
    if n_ht != n_all:
        print(f"WARNING: {n_all - n_ht} .hist( call(s) were not patched")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
