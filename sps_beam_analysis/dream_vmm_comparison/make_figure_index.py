#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""make_figure_index.py -- an index page for the figure directory, so the PNGs
can be browsed and picked up one at a time from the web.

The notes inline their figures as data: URIs, which is right for offline
reading but useless when you want the actual file to drop into a talk. This
publishes the originals.

    python3 make_figure_index.py            # writes figures/index.html
    python3 make_figure_index.py --deploy   # ... and rsyncs to EOS

Deploy target is /eos/user/d/dneff/www/p2_sps/, served at
https://dylan-neff.web.cern.ch/p2_sps/. That directory is hand-published: it is
NOT in dylan-cern-site's PAYLOAD, so the site deploy neither writes nor removes
it.
"""
import argparse
import os
import subprocess

from PIL import Image

import make_report as R

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
REMOTE = "lxplus:/eos/user/d/dneff/www/p2_sps/"

# name -> caption. Order and grouping are the order they are meant to be used.
GROUPS = [
    ("The MPGD2026 slides", "16:9, 2132&times;1200 px, sized for a projector "
     "(type ~2&times; a printed figure's; maps carry a scale bar, not tick "
     "labels). Colour separates the two jobs: DREAM orange / VMM blue is "
     "categorical and rides the dots, bars and lines, while the map ramp is "
     "sequential and changes with the quantity &mdash; purple for efficiency, "
     "amber for gain. Built by <code>figures_deck.py</code>.", [
         ("deck_1_deficit.png",
          "Same pads, same beam, same working point &mdash; DREAM 95.6 %, "
          "VMM 85.3 %. Two per-pad efficiency maps on one purple ramp, and "
          "the pad-by-pad pairing that shows the deficit is a place, not a "
          "scale factor."),
         ("deck_2_gainmap.png",
          "The gain rolls off &times;3.9 across the beam spot and both "
          "readouts measure the same roll-off &mdash; per-pad r = +0.94, same "
          "gradient direction to 5&deg;. The dashed ellipse is the corner "
          "marked on slide 1."),
         ("deck_3_threshold.png",
          "The mechanism: 53 pads in six gain bands, spectra on a log axis so "
          "a gain factor is a sideways shift, one discriminator level drawn "
          "through them, and what each band actually records &mdash; DREAM "
          "98 &rarr; 92 %, the VMM 93 &rarr; 63 %. The blue line is that one "
          "level fitted to all 53 pads at once; each green tick is the same "
          "model inverted for a single pad, so the audit of the claim is on "
          "the slide itself. No chrome &mdash; the talk supplies the title."),
         ("setup_conditions.png",
          "The backing slide for \u201csame chamber, same working "
          "point\u201d: both runs, their HV, gas and trigger, the VMM's gain "
          "/ peaking / threshold registers &mdash; and DREAM's OWN threshold "
          "written into each pad. 26&ndash;29 ADC, flat, because it is 5&sigma; "
          "of each channel's own pedestal; the VMM's is 162 ADC and scatters "
          "26 %."),
         ("setup_groups.png",
          "The backing slide for \u201ca row is a gain band\u201d: which "
          "pads are in each band, and the sort they came from. Sorted on the "
          "pad's median DREAM pulse height &mdash; not efficiency, not "
          "position, not the VMM's own median &mdash; then cut into bands of "
          "equal pad count."),
     ]),
    ("P2 tracking against the uRWELL reference",
     "The three P2 stations judged against external reference tracks: what one "
     "station knows, what the three know together, and why P2 cannot check "
     "its own tracking. Built by <code>figures_track.py</code> from "
     "<code>p2_selftrack.py</code>.", [
         ("track_1_pointing.png",
          "One station against the reference: a 12 mm box, the same on all "
          "three, and a two-pad centroid does not improve it."),
         ("track_2_selftrack.png",
          "The P2-only track against the reference &mdash; offset at the exit "
          "plane, angle, and how the error grows with the lever arm. The "
          "satellites in the angle sit at &plusmn;1 pad over the 620 mm "
          "between the outer stations."),
         ("track_3_illusion.png",
          "The headline: P2_MID checked against the other two P2 stations "
          "reads more than 20&times; better than the same events checked against the "
          "reference, because 70 % of the time all three report the identical "
          "pad."),
         ("track_4_maps.png",
          "Efficiency against the reference track, binned in each station's "
          "own pad frame."),
         ("track_5_purity.png",
          "How much of the efficiency is a real hit: the track-to-cluster "
          "distance per unit area, and the same as efficiency against probe "
          "radius."),
         ("track_6_inpad.png",
          "Every track folded onto the face of the pad it pointed at &mdash; "
          "efficiency and pulse height across a single pad, and the drop at "
          "the edge where the avalanche is split between two."),
     ]),
    ("Dead areas and the bulk pillar lattice",
     "What the reference telescope can see of the amplification gap itself: "
     "one &oslash;6.15 mm support pillar imaged directly, and the "
     "&oslash;0.75 mm bulk pillars measured as a lattice. Built by "
     "<code>figures_pillar.py</code> from <code>p2_pillars.py</code>.", [
         ("pill_1_bigpillar.png",
          "A support pillar, imaged at 0.2 mm from reference tracks alone. "
          "The dashed circle is the gerber, drawn not fitted; efficiency "
          "inside is 1.6&ndash;2.3 % against 91&ndash;98 % outside."),
         ("pill_2_kplane.png",
          "The whole k-plane of the efficiency map: six spots 60&deg; apart at "
          "|k| = 2.42 rad/mm, i.e. a triangular lattice at d = 3.00 mm, at "
          "15&ndash;19&sigma;. No pitch is assumed anywhere in this figure."),
         ("pill_3_avgpillar.png",
          "The average bulk pillar &mdash; every one in the beam spot folded "
          "onto one cell, because no single one has enough tracks on it."),
         ("pill_4_resolution.png",
          "The reference's pointing resolution at the P2 planes, measured two "
          "independent ways, and the contrast that leaves at each scale."),
         ("pill_5_cost.png",
          "What the lattice costs: the charge footprint is the same on all "
          "three chambers, the efficiency footprint follows the gain margin."),
         ("pill_6_mask.png",
          "The dead areas that can be masked, and what masking them does to "
          "the stack &mdash; in this beam spot, almost nothing."),
     ]),
    ("Backup slides", "Kept off the main line so each slide carries one idea. "
     "Built by <code>figures_slide.py</code>.", [
         ("slide_proof.png",
          "That the threshold is quantitatively the whole story: predicted vs "
          "measured efficiency over all 53 pads, and the spread waterfall "
          "0.423 &rarr; 0.355 &rarr; 0.260."),
         ("slide_fix.png",
          "What it costs to fix: efficiency against signal-over-threshold, "
          "all pads and the weakest pad, with &times;3.4 needed to reach "
          "DREAM's 95 %."),
         ("slide_ridge.png",
          "The ridgeline on its own, without the slide chrome &mdash; the "
          "source of deck slide 3."),
         ("slide_ridge_perpad.png",
          "The same ridgeline with each pad's own fitted threshold ticked on "
          "it in green against the one global line in blue, and both colours "
          "named on the figure &mdash; the desk-sized version of "
          "<code>deck_3_threshold.png</code>, eight bands instead of six and "
          "the numbers in a footer."),
         ("slide_full.png",
          "All three ideas on one 16:9 slide, for a single-slide version of "
          "the story."),
         ("slide_ridge_check.png",
          "The ridgeline with each pad's own independently fitted threshold "
          "ticked against the single global line, and those fits plotted "
          "against pad gain."),
     ]),
    ("Does one threshold really fit every pad?",
     "The test the joint fit cannot do on itself: invert the model separately "
     "for each of the 53 pads and compare the answers. Built by "
     "<code>figures_perpad.py</code>.", [
         ("perpad_ridge.png",
          "The ridgeline with the band averaging removed &mdash; 53 rows, one "
          "pad each, every one carrying its own fitted threshold next to the "
          "global line."),
         ("perpad_check.png",
          "What the spread in those 53 answers is worth: 77 % of it lives "
          "inside a gain band, it is 46 ADC against 4 from counting noise, "
          "and at 10 ADC per efficiency point it is the 5.1-point residual "
          "of the global fit in other units."),
     ]),
    ("VMM vs DREAM, the comparison", "Built by <code>figures_dv.py</code>.", [
        ("dv_headline.png", "The headline comparison of the two readouts."),
        ("dv_maps.png",
         "The gradient as two maps and as one profile along its own axis."),
        ("dv_spectra.png",
         "Pulse-height spectra, VMM and DREAM, on the pads they share."),
        ("dv_efficiency.png",
         "Per-pad efficiency against per-pad pulse height, both readouts."),
    ]),
    ("The VMM-only study", "Built by <code>figures.py</code> and "
     "<code>figures_p2out.py</code> from <code>run_46</code>.", [
         ("p2out_headline.png", "P2_OUT: efficiency against per-pad ADC."),
         ("p2out_spectra_map.png",
          "P2_OUT per-pad spectra laid out on the pad map."),
         ("eff_vs_adc.png",
          "Efficiency against per-pad median ADC, all three stations."),
         ("eff_vs_adc_centered.png",
          "The same, with each chip's mean removed &mdash; the within-chip "
          "relation, free of the between-chip threshold trend."),
         ("chamber_maps.png", "Per-pad maps for the three P2 stations."),
         ("spectra_low_high.png",
          "Pooled spectra of the weakest and strongest pads."),
         ("tracked_vs_untracked.png",
          "The untracked <code>adc_vs_ch</code> proxy against the tracked "
          "per-pad median: r = 0.975 on P2_OUT, no decoding needed."),
         ("dnl.png",
          "The VMM ADC's differential non-linearity: a period of 16 codes, "
          "exactly two bins of the <code>adc&gt;&gt;3</code> histogram."),
     ]),
]

CSS_EXTRA = """
main{max-width:1100px}
.figrow{display:flex;gap:1.1rem;align-items:flex-start;padding:1.1rem 0;
border-bottom:1px solid var(--line)}
.figrow:last-child{border-bottom:0}
.figrow a.thumb{flex:0 0 340px;display:block}
.figrow img{width:100%;height:auto;border:1px solid var(--line);
border-radius:6px;display:block}
.figrow .meta{flex:1;min-width:0}
.figrow h3{margin:0 0 .3rem;font-size:1rem;font-family:ui-monospace,
SFMono-Regular,Menlo,monospace}
.figrow p{margin:0 0 .5rem;color:var(--fg2);font-size:.92rem}
.figrow .dl{font-size:.85rem;color:var(--mut)}
.figrow .dl a{margin-right:.9rem}
@media(max-width:700px){.figrow{display:block}.figrow a.thumb{width:100%}}
"""


def rows(items):
    out = []
    for name, cap in items:
        p = os.path.join(FIG, name)
        if not os.path.isfile(p):
            print(f"  missing, skipped: {name}")
            continue
        w, h = Image.open(p).size
        kb = os.path.getsize(p) / 1024
        out.append(f"""<div class="figrow">
<a class="thumb" href="{name}"><img src="{name}" alt="{name}" loading="lazy">
</a>
<div class="meta"><h3>{name}</h3><p>{cap}</p>
<div class="dl"><a href="{name}">open full size</a>
<a href="{name}" download>download</a>
{w}&times;{h} px &middot; {kb:.0f} kB</div></div></div>""")
    return "\n".join(out)


def build():
    body = []
    for title, blurb, items in GROUPS:
        body.append(f"<h2>{title}</h2>\n<p class='note'>{blurb}</p>\n"
                    f"{rows(items)}")
    h = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<meta name="robots" content="noindex, nofollow">
<title>P2_OUT figures &mdash; VMM3a vs DREAM</title>
<meta name="description" content="The original PNGs behind the P2_OUT VMM/DREAM
comparison: the MPGD2026 slides, the backup slides and the analysis figures.">
<style>{R.CSS}{CSS_EXTRA}</style></head><body>
<main>
<h1>P2_OUT figures &mdash; VMM3a vs DREAM</h1>
<p class="sub">The originals, one file each. P2 basket, SPS H4, July 2026. VMM3a
<code>run_46 / cfg_gain4.5_peaktime200</code> against DREAM
<code>eff_nominal_1</code>, matched at mesh 450&nbsp;V / drift 750&nbsp;V,
53 pads under the beam.</p>
<p class="sub">The written-up versions:
<a href="../notes/p2out-vmm-three-slides.html">three slides for MPGD2026</a>,
<a href="../notes/p2out-vmm-threshold-weak-pads.html">why the VMM loses the
weak pads</a>,
<a href="../notes/p2out-vmm-vs-dream-pad-gain.html">VMM vs DREAM pad gain</a>.
Source: <code>nTof_x17/p2_sps/</code>, branch
<code>p2-sps-vmm-vs-dream</code>.</p>
{"".join(body)}
</main></body></html>
"""
    out = os.path.join(FIG, "index.html")
    open(out, "w").write(h)
    print(f"wrote {out}  ({os.path.getsize(out) / 1024:.0f} kB)")
    return out


def deploy():
    """rsync the PNGs and the index. No --delete: this directory is published by
    hand and nothing else writes it, but the habit is worth keeping."""
    files = [os.path.join(FIG, f) for f in sorted(os.listdir(FIG))
             if f.endswith((".png", ".html"))]
    subprocess.run(["ssh", "lxplus", "mkdir -p /eos/user/d/dneff/www/p2_sps"],
                   check=True)
    subprocess.run(["rsync", "-vz", "--no-perms", "--no-owner", "--no-group",
                    "--omit-dir-times", *files, REMOTE], check=True)
    print("\nhttps://dylan-neff.web.cern.ch/p2_sps/")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--deploy", action="store_true")
    a = ap.parse_args()
    build()
    if a.deploy:
        deploy()
