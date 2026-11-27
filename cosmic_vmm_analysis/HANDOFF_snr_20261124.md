# Handoff: section 8 of the missed-muons note as signal-to-noise (24-27 Nov 2026)

Session "VMM_CosmicTB Analysis". This picks up from `HANDOFF_analysis_20261122.md`, whose state is otherwise
unchanged. It covers only the S/N work on section 8 and the note repairs made at the same time.

## 1. Question and answer

The question was that effects 1 (noise in charge) and 2 (trim range / clipped channels) of section 8, "Checking that
nothing is mixed up", were not clear from the deck figures. The request was to express them as a signal-to-noise
ratio.

- **S** = median muon pulse height (largest pad of the found muons, ADC - 40 baseline) / 1.07 ADC/mV. The scale is
  the vmmsc Offline_ADC calibration of Nov 2025; it is gain-independent because the VMM ADC reads voltage.
- **T** = threshold above baseline for a noise rate R per channel, median over the chamber. Two versions:
  - `T_need` = per-channel targets in the bscan_target reports, i.e. what the noise asks for (effect 1 only);
  - `T_set` = `dist_mV` of the efficiency pipeline (`gain_scan_eff.csv` `thr_mV`), i.e. what was applied after the
    trim (effect 1 + clipping).
- Fixed R → fixed T/σ for Gaussian noise (R ∝ exp(-T²/2σ²)), so S/T = S/N × a constant.

P2_MID, R = 10 Hz:

| | g3 | g6 | g9 |
|---|---|---|---|
| S (mV) | 100 | 195 | 287 |
| S / S(g3) | 1 | 1.95 | 2.87 |
| T_need (mV) | 10.4 | 24.5 | 46.1 |
| T_need / T_need(g3) | 1 | 2.35 | 4.43 |
| S/T_need | 9.6 | 8.0 | 6.2 |
| S/T_set | 9.0 | 7.9 | 4.8 |
| clipped channels | 0/257 | 19/353 | 168/289 |
| silent | 13.6 % | 14.2 % | 24.5 % |

S/T_set for all chambers (g3 / g6 / g9): IN 10.7 / 9.7 / 5.9, MID 9.0 / 7.9 / 4.8, OUT 6.3 / 6.4 / 4.4.

- **Effect 1:** the signal follows the nominal gain, but the noise threshold grows FASTER than the gain. S/N on MID
  goes 1 : 0.83 : 0.65 (IN 1 : 0.86 : 0.68). In fC the signal is a flat 32-33 fC and the threshold goes 3.5 → 4.1 →
  5.1 fC (per channel 1.19× / 1.47×).
- **Effect 2:** at g9, clipping costs a further 24-25 % of S/T on MID and IN (7 % on OUT); at g3 and g6 it costs 1-8 %.
  It does not cost muons, because the clipped channels are mostly the quiet ones of their chip (median target 0.8× the
  chip median).
- **S/N does not explain the 13-14 % silent share.** At g3-g6 OUT runs at S/T ≈ 6.3 and MID at 8-9, yet both lose
  12-14 %. The silent share rises only at g9, where S/T < 5.
- **Not understood: why the noise grows faster than the gain.** The preamplifier noise is before the gain stage, so
  S/N should be flat. Candidates: non-Gaussian tails (then 10 Hz is not a fixed number of σ), a baseline shift with
  gain, or a change of noise bandwidth. The direct test is the pedestal width in mV at each gain, or the slope of
  per-channel baseline-scan curves (rate vs threshold) at g3/g6/g9.

The user's own confusion was "signal should grow with gain too". It does; the first table in the note was in fC,
which hid that. The note now uses mV with "relative to gain 3" rows.

## 2. What changed

### Code (all merged into `main` and pushed, HEAD 619c24c)
- `note/snr_figs.py` (new): `snr_figs.py GS_OUT_DIR REPORT_DIR OUT_DIR` → `c_snr.png` (S/T vs gain, need vs set, 3
  chambers) and `snr_numbers.csv` (all g × R × chamber). Needs pandas; locally use `~/miniconda3/bin/python`.
  Inputs = `$W/analysis/gain_scan/out/*` and `$W/analysis/gain_scan/reports/hvon_g{3,6,9}.0_R{2,10,50}_report.json`.
  Last run copies are in the job scratch only; re-fetch with scp from lxplus.
- `note/muon_funnel_page.py`:
  - `SNR` dict + `svg_snr_pulses()`, `svg_snr_rate()`, `svg_snr_ratio()` (inline SVG, CSS vars --g3/--g6/--g9,
    --k-sil, --k-ok; light and dark both checked);
  - in section 8, a block between `<!-- snr:begin -->` and `<!-- snr:end -->`: steps 1-4, the c_snr.png figure, the
    mV table, three bullets, and a numbers footnote.
  - The pulse and noise shapes in step 1 are drawn by hand and the rate-curve shape in step 2 is schematic (Gaussian,
    100 kHz prefactor). Only S and T are measured; the captions say so.

### CERN note (cern-site, all deployed and pushed)
`notes/unlisted/2026-11-19-where-do-the-missed-muons-go-the-vmm-cosmic-stand-step-by-step.html`
- `3238500`: first S/N block + `c_snr.png` in `_files/`
- `30248db`: **restored section 9 "The same pads with DREAM"** and its top "Update 22 Nov" box. They were dropped by
  the 22 Nov rebuild (aefb192). The fragment was copied verbatim from `c53558c`, so its VMM numbers predate the 92-run
  / paddle-cut rebuild.
- `bda2348`: S/N block replaced by the pedagogical version (snr markers).
- Since then, another session added section 10 "det5 at P2_MID and the read-out swap" (`28f1453`). Current order: …
  8 (with S/N), 9 same pads with DREAM, 10 det5 at MID, 11 Open points, 12 References.

## 3. TRAP: do not full-rebuild the note from main

`muon_funnel_page.py` on main generates only sections 1-8 + Open points + References.
- Section 9 is added by `dream_same_pads_section.py` (`insert_html`, env `SAME_PADS_HTML`) and lives only on branch
  `note-dream-vmm-same-pads` (commit 4bf96d6, worktree `.claude/worktrees/note-dream-vmm-same-pads`).
- Section 10 (det5) comes from the det5 session (see memory `project_det5_mid_night_20261123.md` /
  `HANDOFF_det5_night_20261124.md`, branch note-dream-vmm-same-pads / det5-mid-chamber).

A rebuild with the recipe in HANDOFF_analysis_20261122.md would silently drop sections 9 and 10 again. Until those
branches are merged into main and the generator inserts both, edit the live HTML by splicing:
- the S/N block can be regenerated in place with the snr markers (the pattern was: exec the `SNR`..`svg_snr_ratio`
  definitions + the `A('<!-- snr:begin -->')`..`A('<!-- snr:end -->')` lines of the generator and replace between
  the markers);
- section 9 sits between `<!-- same-pads:begin/end -->` and `<!-- same-pads:top/topend -->`.

Then `cd ~/Documents/PostDocSaclay/cern-site && python3 scripts/build_hub.py && git commit && scripts/deploy.sh &&
git push`.

## 4. Open / next
1. **Noise faster than gain:** measure the pedestal (noise) width in mV at g3/g6/g9, or fit log-rate vs T² of the
   bscan baseline-scan curves per channel; it tests the Gaussian assumption behind S/T ≈ S/N. This needs the scan
   points on the DAQ box (`~/trim_work/out/gain_20261118/`, bscan_analyze curves).
2. Merge `note-dream-vmm-same-pads` (+ det5 section code) into main so a full rebuild keeps sections 9-10.
3. The section 9 VMM numbers predate the 92-run / paddle-cut funnel; recompute with `dream_same_pads_section.py` on
   `analysis/muon_funnel_physical`.
4. The user's last message was cut off: "I have moved det2 from the P2_MID position to the …". The new position is
   unknown. Ask, then check chamber↔position in the logbook before any per-chamber numbers on new runs (the det5 swap
   of 23-25 Nov put det5 at MID; this is a further change).
5. Small cosmetic: in step 1 of the S/N part, the gain-9 "noise ×3 only: 31 mV" label sits slightly over the noise
   trace.
