# MPGD26 conference figures

One place holding a copy of every figure the talk uses (4 Sep 2026, 15 min + 5 Q),
gathered from the five trees that actually produce them.

```
python3 gather_figures.py            # copy / refresh, rewrite INVENTORY.md
python3 gather_figures.py --check    # report what is stale, copy nothing
```

`figures/` is regenerated output and is git-ignored (`*.png`, `*.pdf`); the
manifest inside `gather_figures.py` and `INVENTORY.md` are the tracked record.
Destination names carry the draft-deck slide number, so each act directory reads
in talk order. `9_spare/` holds figures that were gathered for the deck but do
not currently sit on a slide.

## Where the figures come from

| tree | figures | produced by |
|---|---|---|
| `data/Cosmic_Bench/Analysis/conference/` | 12 bench figures (`01_`–`13_`) | `cosmic_bench_analysis/22_conference_figures.py` (PNG + PDF) |
| `P2_basket_analysis/reports/mpgd26_sps_beam_2026-08/figs/` | 26 beam figures | `sps_beam_analysis/mpgd26_figs/fig_*.py` (restored 2026-08-17) |
| `P2_basket_analysis/mpgd2026/figs/` | the 8 talk-only figures | `mpgd2026/make_talk_figs.py` |
| `P2_basket_analysis/gas_studies/figs/` | 2 gas-study figures | `gas_studies/` |
| `cern-site/.../deck_files/` | 3 × `coincidence_*` | rendered on lxplus (`nTof_x17/mpgd26`), copy only |

**The beam figures now have a producer again.** They were made in a scratch tree
that no longer existed, so nothing in the repo could re-make them. The scripts
are back under `sps_beam_analysis/mpgd26_figs/` (`fetch_from_eos.sh` →
`aggregate.py` → `fig_*.py`), reading a workspace built from the EOS copy — which
is *newer* than the LaCie one, e.g. the old LaCie stage-29 products still carry
the pre-rename timing columns. Every one of them writes PDF beside PNG now.

Still PNG-only, 12 of 55: the three lxplus-rendered `coincidence_*`;
`urw_panels_P2_MID` / `urw_summary_highstat` (stage outputs of
`urw_p2_efficiency.py` — a re-run would be needed); the two gas-study and two
lifetime-autopsy figures (their stages are in-repo and could emit PDF);
`pad_sector_layout`; and `event_display` / `timing_vs_mesh_bothruns`, the only
two with **no** recoverable producer.

## Roadmap §5a scorecard — the bench ↔ beam / DREAM ↔ VMM comparisons

| item | figure | state |
|---|---|---|
| 5a.1 DREAM vs VMM, same detector & HV | `s23_dream_vs_vmm`, `s24_vmm_threshold` | **done** — same uRWELL tracks, same fiducial, same probe radius, 78 sub-runs. 0.96/0.97/0.95 DREAM vs 0.85 VMM at P2_OUT, and the gap is shown to be the discriminator threshold by three independent handles. |
| 5a.2 beam ε(HV) overlaid on the bench | `s21_bench_beam_mesh`, **`bench_beam_drift`** | **done, both halves.** The drift overlay is new: plotted against drift *field*, four chambers over two gases fall on one curve and agree to 2–4 points at the plateau. |
| 5a.3 VMM optimum vs non-optimal config | `s35_vmm_config_scan`, `s08_snr_matrix` | **done differently, and better** — a full July (gain, peaking) scan on the new detectors, plus the Nov-2025 SNR grid, plus the resolution of the apparent contradiction (SNR-optimal ≠ efficiency-optimal). |
| 5a.4 efficiency map, beam vs bench | **`bench_beam_maps`** | **done** — bench sliding map beside the beam per-pad map on one colour scale, medians 0.950 vs 0.969. |
| 5a.5 beam timing vs the bench | `s22_timing_campaigns`, `s22_timing_vs_drift_magboltz` | **done** — 15.5 / 18.4 / 22.4 ns per station against the bench and the Garfield floor, with the gas/electronics decomposition. |

## Chamber identity across the two campaigns

From the **2026-07-28 logbook**, which is the record to trust here — the DAQ
`run_config.json` descriptions disagree with each other between runs and
`p2_qa_config.py:624` disagrees with both:

| beam station | chamber | note |
|---|---|---|
| **P2_MID** | **det1** | unchanged for the whole beam test — *the same chamber the cosmic bench used as its worked example* |
| **P2_OUT** | **det3** | unchanged for the whole beam test |
| P2_IN | a first chamber 23–27 Jul, then the **CERN-built** one **from 28 Jul** | the CERN chamber's drift frame needs modification for the old HV connection |
| — | det2, det4 | never reached the beam: leaky drift frames |

So `bench_beam_mesh` and `bench_beam_drift` are **same-chamber** comparisons, not
cross-chamber ones: det1 ↔ P2_MID in both, det3 ↔ P2_OUT in the drift figure.
Only the gas (Ar/iC₄H₁₀ 95/5 vs Ar/CO₂/iC₄H₁₀ 93/5/2) and the track reference
(M3 cosmics vs uRWELL muons) differ — which is why the drift curves land within
2–4 points of each other.

`fe55_bench_beam_P2_OUT` is **parked**: it was built on the Fe55 registry's claim
that the 18 Jul scan measured "P2_OUT" as det2, which the logbook contradicts
(det2 never left Saclay). The analysis is sound once the label is; the function
is kept in `make_talk_figs.py` and must be called explicitly.

Nothing in §5a needed new *analysis*, but it did need new *processing*: the
mesh turn-on now runs 330–450 V because `low_mesh_scan_1`, `p2in_hvrange_1` and
`p2in_hvrange_2` were pushed through the uRWELL-referenced stage on condor
(30 sub_runs, 2026-08-17).
