# FEC scripts from colleagues — kept for reference, not used

The first pcapng readers for the VMM3a/FEC chain, written by colleagues in
Nov 2025 and carried in `P2_basket_online_analysis` until it was retired (tag
`legacy-final`). **Nothing uses them any more**: decoding is
[p2-decode](https://github.com/akallitss/p2-decode), the channel → pad map is
[p2-mapping](https://github.com/akallitss/p2-mapping), live monitoring is
[p2-monitor](https://github.com/akallitss/p2-monitor). They parse with scapy
one packet at a time and are not maintained.

| File | What it was |
|---|---|
| `FEC_Datadescrambler_Default.py` | follows a pcapng in real time, live ADC and channel histograms; assumes 4 VMMs (ids 0–3) |
| `FEC_Datadescrambler.py` | the same, and also starts `dumpcap` on a hardcoded interface (`enp2s0`) |
| `FEC_Analysis.py` | offline analysis with strip remapping from an Excel pinout and a 2D hit map (D. Baudin; Windows paths, not adapted to hybrid-only data) |

Copied unchanged from `P2_basket_online_analysis` at `legacy-final` (0f3043e).
