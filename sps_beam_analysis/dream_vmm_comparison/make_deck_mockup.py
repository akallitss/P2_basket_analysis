#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""make_deck_mockup.py -- the three-slide MPGD2026 sequence, mocked up in
running order with what to say over each one.

The slides themselves come from figures_deck.py; this only frames them, so the
numbers in the speaker notes are read from the same JSON the figures use and a
re-run keeps prose and pictures in step.
"""
import json
import os

import pandas as pd

import figures_deck as D
import figures_slide as FS
import make_report as R

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")

EXTRA = """
main{max-width:1240px}
.slide{margin:1.6rem 0 0}
.slide img{border:1px solid var(--line);border-radius:8px}
.say{background:var(--card);border-left:3px solid var(--warn);
padding:.85rem 1.1rem;border-radius:6px;margin:.9rem 0 0}
.say b{color:var(--warn)}
.order{background:var(--card);border:1px solid var(--line);border-radius:8px;
padding:.9rem 1.1rem;margin:1.2rem 0;font-size:.94rem;line-height:1.9}
.order .now{color:var(--acc);font-weight:700}
.order .prev{color:var(--mut)}
h2 small{font-weight:400;color:var(--fg2);font-size:.8em}
"""


def main():
    n = json.load(open(os.path.join(DATA, "threshold_model_P2_OUT.json")))
    g = pd.read_csv(os.path.join(DATA, "compare_dream_vmm_P2_OUT.csv"))
    g = g[g["use"]]
    fix = {r["factor"]: r for r in n["fix"]}
    f2, f15 = fix[2.0], fix[1.5]
    span = g["amp_med_d"].max() / g["amp_med_d"].min()
    ev = g["eff_v"].sort_values().to_numpy()
    worst_live = ev[1]
    # the band numbers quoted below are read off the slide's own bands, so they
    # cannot drift from the figure when NGROUP changes
    _g, _H, _bw, _S, _n = FS.load()
    bands = FS.groups(_g, _H, _bw, ngroup=D.NGROUP)
    b_lo, b_hi = bands[0], bands[-1]

    h = f"""<title>Three slides: the VMM efficiency deficit</title>
<meta name="description" content="A three-slide MPGD2026 sequence for P2_OUT:
the VMM is less efficient than DREAM, the deficit is a place on the chamber,
and one discriminator level inside the Landau explains it.">
<style>{R.CSS}{EXTRA}</style>
<main>
<h1>Three slides: why the VMM is less efficient than DREAM</h1>
<p class="sub">A mock-up for MPGD2026. Three 16:9 slides, almost no body text
&mdash; each one is a picture with a headline and a single grey line of
provenance. They are built by <code>figures_deck.py</code> from the same tables
as <a href="p2out-vmm-threshold-weak-pads.html">the full note</a>, so the
numbers move together.</p>

<p class="sub">Sized for a projector, not a desk: type is ~2&times; a printed
figure's, the pad maps drop their tick labels for a scale bar, and slide 3 runs
{D.NGROUP} gain bands instead of eight so the rows are tall enough to read from
the back. Colour does two separate jobs and they are kept apart &mdash;
<b>which readout</b> is categorical (DREAM orange, VMM blue, on the dots, bars
and lines), while <b>what is mapped</b> is a sequential ramp that changes with
the quantity: <b style="color:#7250a4">purple for efficiency</b> on slide 1,
<b style="color:#ac5f1a">amber for gain</b> on slide 2. Both maps on a slide
share one ramp and one colorbar, which is what makes them comparable. Each ramp
passes the ordinal checks against the slide background, including the one that
matters in a hall: the pale end clears 2:1 contrast, so it does not wash out to
paper.</p>

<div class="order">
<span class="prev">&hellip;&nbsp; DREAM on the SPS beam &mdash; efficiency,
timing, space</span><br>
<span class="prev">&hellip;&nbsp; VMM3a &mdash; timing, space: it matches</span>
<br>
<span class="now">&rarr;&nbsp; 1. &nbsp;but the efficiency does not</span><br>
<span class="now">&rarr;&nbsp; 2. &nbsp;the deficit is a place on the chamber,
and both readouts see that place</span><br>
<span class="now">&rarr;&nbsp; 3. &nbsp;it is one threshold sitting inside the
Landau</span><br>
<span class="prev">&hellip;&nbsp; what it would take to fix</span>
</div>

<h2>Slide 1 &nbsp;<small>&mdash; the observation</small></h2>
<div class="slide">
<img src="figures/deck_1_deficit.png" alt="Two per-pad efficiency maps of
P2_OUT, DREAM and VMM3a on one colour scale, and the pad-by-pad pairing of the
two efficiencies.">
</div>
<div class="say">
<b>Say:</b> same chamber, same beam, same mesh and drift voltage, the same
uRWELL tracks and the same selection &mdash; only the electronics differ.
DREAM records {n['eff_dream_all'] * 100:.1f}&nbsp;% of the tracks that point at
a pad; the VMM records {n['eff_obs_all'] * 100:.1f}&nbsp;%. Then point at the
right-hand panel: the deficit is <b>not</b> a uniform scale factor. DREAM holds
every live pad above 92&nbsp;%, while the VMM runs from 96&nbsp;% down to
{worst_live * 100:.0f}&nbsp;% &mdash; and the pads it loses are all in one
corner.
</div>

<h2>Slide 2 &nbsp;<small>&mdash; whose fault is it?</small></h2>
<div class="slide">
<img src="figures/deck_2_gainmap.png" alt="Per-pad gain maps of P2_OUT for
DREAM and VMM3a on one colour scale, with the fitted gradient direction, and a
pad-by-pad scatter of the two relative gains.">
</div>
<div class="say">
<b>Say:</b> the same corner, now in pulse height. The gas gain rolls off by a
factor {span:.1f} across the beam spot, one smooth gradient, and
<b>both readouts measure the same map</b> &mdash; pad for pad, r&nbsp;=&nbsp;+0.94,
same direction to 5&deg;. So the variation is the chamber's, not the
electronics'. Leave them with the one loose end that sets up slide&nbsp;3: the
VMM's copy of that map is systematically <i>flatter</i>.
</div>

<h2>Slide 3 &nbsp;<small>&mdash; the mechanism, and the close</small></h2>
<div class="slide">
<img src="figures/deck_3_threshold.png" alt="Ridgeline of the per-pad DREAM
pulse-height spectra banded by gain on a logarithmic axis, with the VMM
discriminator level drawn through them and the per-band recorded fraction
alongside.">
</div>
<div class="say">
<b>Say:</b> every row is the same Landau; the axis is logarithmic, so a gain
factor is a pure sideways shift. The blue line is the VMM's discriminator
&mdash; one level, all six chips at <code>sdt&nbsp;=&nbsp;224</code>. Walk down
the rows: the solid orange area is what falls below it. On the strong pads that
is the bottom of the tail; on the weak pads the peak itself has slid onto the
line. The bars on the right are the consequence &mdash; DREAM holds
{b_hi['eff_d'] * 100:.0f}&nbsp;&rarr;&nbsp;{b_lo['eff_d'] * 100:.0f}&nbsp;%
across the {D.NGROUP} gain bands, the VMM falls
{b_hi['eff_v'] * 100:.0f}&nbsp;&rarr;&nbsp;{b_lo['eff_v'] * 100:.0f}&nbsp;%.
<br><br>
One fitted level, {n['T']:.0f} DREAM ADC, cut into DREAM's own per-pad spectra
reproduces the VMM's efficiency pad by pad: r&nbsp;=&nbsp;{n['r_pred']:+.2f},
{n['eff_pred_all'] * 100:.1f}&nbsp;% predicted against
{n['eff_obs_all'] * 100:.1f}&nbsp;% measured. And it closes slide&nbsp;2's loose
end: cutting the spectra at that level is also what flattens the VMM's copy of
the gain map.
</div>

<h2>What to say next &mdash; the options, cheapest first</h2>
<ul>
<li><b>Raise the VMM front-end gain.</b> A register change. It halves the
threshold in charge and buys back most of the deficit &mdash;
&times;2&nbsp;signal-over-threshold takes the fleet to
{f2['eff_all'] * 100:.1f}&nbsp;% and the weakest pad from
{fix[1.0]['eff_min'] * 100:.0f}&nbsp;% to
{f2['eff_min'] * 100:.0f}&nbsp;%. Noise scales with it and the ADC saturates
earlier.</li>
<li><b>Lower <code>sdt</code>.</b> Already at the campaign minimum on these
chips; how much further it can go is set by the noise, which this measurement
does not constrain. Needs a threshold scan (224 / 208 / 192) with beam.</li>
<li><b>Raise the mesh.</b> Works &mdash; &times;1.5 alone is worth
{f15['eff_all'] * 100:.1f}&nbsp;% &mdash; but it is multiplicative, so the
gradient survives and the strong pads move too.</li>
<li><b>Use the per-channel trim DACs against the measured gain map.</b> The
right shape of fix, but the trim range is a few mV against a factor
{span:.1f} in gain, and buying uniformity that way costs noise.</li>
<li><b>Fix the chamber.</b> {n['rms_dream_raw']:.2f} relative rms, about
13&nbsp;% per 10&nbsp;mm, on the amplification side. This is the only option
that helps <i>both</i> readouts.</li>
</ul>

<h2>Deliberately not on these three</h2>
<p>Kept as backup slides so the main line stays one idea per slide:</p>
<ul>
<li><b>The proof panel</b> (<code>figures/slide_proof.png</code>) &mdash; the
predicted-vs-measured scatter over all {n['npad']} pads, and the spread
waterfall {n['rms_dream_raw']:.2f} &rarr; {n['rms_dream_cut']:.2f} &rarr;
{n['rms_vmm']:.2f} that shows the threshold plus a
+{n['offset_adc']:.0f}-count additive ADC offset accounts for the compression
with nothing left over.</li>
<li><b>The fix curve</b> (<code>figures/slide_fix.png</code>) &mdash;
efficiency against signal-over-threshold, all pads and the weakest pad, with
{n['factor_for_95']:.1f}&times; needed to reach DREAM's 95&nbsp;%.</li>
<li><b>Caveats</b> &mdash; the two runs are five days apart; the fitted level is
a charge threshold and would absorb anything else that bites at the low end;
the additive offset is measured, not explained; this is P2_OUT only; and
nothing here measures the noise occupancy a lower threshold would bring.
They are written out in
<a href="p2out-vmm-threshold-weak-pads.html">the full note</a>.</li>
</ul>

<h2>Getting the files</h2>
<p>The figures on this page are inlined, which is no use when you want the PNG
itself. The originals are published one file each, with captions and pixel
sizes, at <a href="../p2_sps/"><b>dylan-neff.web.cern.ch/p2_sps/</b></a> &mdash;
the three slides, the backup slides and every analysis figure behind them.</p>

<p class="note">Slides: <code>p2_sps/figures/deck_{{1_deficit,2_gainmap,
3_threshold}}.png</code>, 13.33&times;7.5&nbsp;in at 160&nbsp;dpi
(2132&times;1200&nbsp;px), built by <code>p2_sps/figures_deck.py</code> on
branch <code>p2-sps-vmm-vs-dream</code>.
Sources: VMM3a <code>run_46 / cfg_gain4.5_peaktime200</code>, DREAM
<code>eff_nominal_1</code>, both at mesh 450&nbsp;V / drift 750&nbsp;V.</p>
"""
    out = os.path.join(HERE, "deck_mockup.html")
    open(out, "w").write(h)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
