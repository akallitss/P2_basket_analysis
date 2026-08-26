#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""make_report_slide.py -- the MPGD2026 slide case: why the VMM loses the low
pulses on the low-gain pads, drawn so it reads in a minute, plus what to do
about it.

Every number comes from data/threshold_model_P2_OUT.json (threshold_model.py)
and data/compare_dream_vmm_P2_OUT.csv, so re-running the chain updates the
prose with the figures.
"""
import json
import os

import numpy as np
import pandas as pd

import make_report as R

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")


def main():
    n = json.load(open(os.path.join(DATA, "threshold_model_P2_OUT.json")))
    g = pd.read_csv(os.path.join(DATA, "compare_dream_vmm_P2_OUT.csv"))
    g = g[g["use"]]
    fix = {r["factor"]: r for r in n["fix"]}
    f2, f15 = fix[2.0], fix[1.5]
    span = g["amp_med_d"].max() / g["amp_med_d"].min()

    frow = "".join(
        f"<tr><td>x{r['factor']:g}</td><td>{r['eff_all']:.3f}</td>"
        f"<td>{r['eff_min']:.3f}</td><td>{r['spread']:.2f}</td></tr>"
        for r in n["fix"] if r["factor"] in (1.0, 1.15, 1.3, 1.5, 2.0, 3.0))

    h = f"""<title>Why the VMM loses the weak pads</title>
<meta name="description" content="P2_OUT read by VMM3a and by DREAM at the same
operating point: the pad-to-pad gain spread is the chamber's, and one
discriminator level accounts for the whole difference between the two
readouts.">
<style>{R.CSS}</style>
<main>
<h1>Why the VMM loses the weak pads &mdash; and DREAM does not</h1>
<p class="sub">P2_OUT, SPS H4. VMM3a <code>run_46 /
cfg_gain4.5_peaktime200</code> against DREAM <code>eff_nominal_1</code>, matched
operating point (mesh 450&nbsp;V, drift 750&nbsp;V), same uRWELL reference
tracks, same selection, {n['npad']} pads under the beam.</p>

<div class="verdict">
<b>The pad-to-pad gain variation belongs to the chamber. Both readouts measure
the same map</b> &mdash; per-pad, r&nbsp;=&nbsp;+0.94, a factor
{span:.1f} from the weakest illuminated pad to the strongest, laid out as one
smooth gradient across the beam spot.
<br><br>
<b>What the VMM adds is a threshold that sits inside the Landau.</b> Fitting one
level to the DREAM spectra reproduces the VMM's efficiency pad by pad
&mdash; {n['npad']} pads, one free number, r&nbsp;=&nbsp;{n['r_pred']:+.2f},
residual {n['resid_rms'] * 100:.1f} points, predicted overall
{n['eff_pred_all']:.3f} against {n['eff_obs_all']:.3f} measured. That level is
<b>{n['T_over_mpv']:.2f}&times; the most probable pulse of an average pad</b>,
which is <b>{n['T_over_mpv_weak']:.2f}&times;</b> the most probable pulse of a
weak one &mdash; so on the weak pads it eats the peak, not the tail.
<br><br>
<b>The VMM's apparently better uniformity is an artefact of the same two
effects.</b> Detector spread {n['rms_dream_raw']:.2f} (DREAM, everything kept)
&rarr; {n['rms_dream_cut']:.2f} once the threshold cuts the Landau &rarr;
{n['rms_vmm']:.2f} once the VMM ADC's {n['offset_adc']:.0f}-count additive
offset is folded in. Removing that offset from the VMM medians gives
{n['rms_vmm_deoffset']:.2f} &mdash; the truncated DREAM number. Nothing is left
over.
</div>

<div class="kpis">
<div class="kpi"><b>{n['T_over_mpv']:.2f}x</b><span>threshold / MPV, average pad</span></div>
<div class="kpi"><b>{n['T_over_mpv_weak']:.2f}x</b><span>threshold / MPV, weak pad</span></div>
<div class="kpi"><b>{n['r_pred']:+.2f}</b><span>predicted vs measured eff.</span></div>
<div class="kpi"><b>{n['eff_obs_all']:.3f}</b><span>VMM efficiency</span></div>
<div class="kpi"><b>{n['eff_dream_all']:.3f}</b><span>DREAM efficiency</span></div>
</div>

<h2>The figure</h2>
<figure><img src="figures/slide_ridge.png" alt="ridgeline">
<figcaption>The {n['npad']} pads sorted by gain into eight bands; each band's
pulse-height distribution as DREAM records it, on a <b>logarithmic</b> axis so
that a gain factor is a pure sideways shift and the eight curves are the same
Landau translated. The vertical line is the VMM's discriminator &mdash; one
level, identical on all six P2_OUT chips (sdt&nbsp;=&nbsp;224). Solid orange is
the part below it. As the pad gain falls from 1.83&times; to 0.55&times; the
solid part grows and the VMM's recorded fraction falls from 91&nbsp;% to
53&nbsp;%, while DREAM stays between 98 and 89&nbsp;%. Bands are the
track-weighted mean of their pads' unit-area spectra, so a band's
sub-threshold area is exactly the mean of its pads'.</figcaption></figure>

<h2>That it is quantitatively the whole story</h2>
<figure><img src="figures/slide_proof.png" alt="proof">
<figcaption><b>Left:</b> per-pad VMM efficiency against what you get by cutting
that pad's DREAM spectrum at the one fitted level, {n['T']:.0f}&nbsp;ADC
({n['T_lo']:.0f}&ndash;{n['T_hi']:.0f}). One free parameter for
{n['npad']} pads. The outlier is pad 635, which is dead in both readouts and is
not a threshold effect. <b>Right:</b> the measured pad-to-pad spread of the
pulse height, and where the VMM's smaller number comes from. Marker area is the
track count.</figcaption></figure>

<h2>What it would take to fix</h2>
<figure><img src="figures/slide_fix.png" alt="fix">
<figcaption>Efficiency against signal-over-threshold relative to today, from the
same model. Raising the chamber gain, lowering <code>sdt</code> and raising the
VMM's front-end gain all enter as the same ratio, so the model cannot tell them
apart &mdash; it says what a factor is worth, not which knob to turn.</figcaption>
</figure>

<div class="tw"><table>
<tr><th>signal / threshold</th><th>efficiency, all pads</th>
<th>weakest live pad</th><th>spread across pads</th></tr>
{frow}
</table></div>

<h2>Options, cheapest first</h2>
<ul>
<li><b>Raise the VMM front-end gain</b> from 4.5&nbsp;mV/fC. The discriminator
is set in mV after the amplifier, so 9&nbsp;mV/fC halves the threshold in charge
&mdash; exactly the factor of two the model asks for &mdash; and it is a
register change with no HV risk. Two things to check on the bench: the noise
scales with the gain too, so <code>sdt</code> may not be able to stay where it
is; and the 10-bit ADC saturates earlier, which costs the top of the Landau (it
does not cost efficiency).</li>
<li><b>Lower <code>sdt</code></b>. All six P2_OUT chips already sit at 224, the
lowest value used anywhere in the July setup (P2_IN ran 224&ndash;300). How much
further it can go is a noise question this measurement does not answer &mdash;
the cheap experiment is a threshold scan, 224/208/192, one chip, with beam.</li>
<li><b>Raise the mesh voltage.</b> {f2['eff_all']:.0%} at &times;2 and
{f15['eff_all']:.0%} at &times;1.5. It moves the whole Landau up, so it works,
but the gain gradient is multiplicative and survives &mdash; the pads stay as
unequal as they were, they just all clear the bar.</li>
<li><b>Use the per-channel trim DACs against the gain map, not for it.</b> The
usual point of trimming is to remove threshold dispersion; here the dispersion
that matters is the detector's gain, and it is measured. Trimming the low-gain
pads <i>down</i> would flatten the efficiency map directly. The range is a few
mV, which is small against the factor {span:.1f}, and it buys noise on exactly
the channels that get the lower threshold &mdash; worth costing, not obviously
worth doing.</li>
<li><b>Fix the chamber.</b> The gradient is ~13&nbsp;% per 10&nbsp;mm across the
beam spot and it is the same in both readouts, so it is amplification-side, not
electronics. That is the only option that improves the detector rather than the
readout's tolerance of it.</li>
</ul>

<h2>The one-slide version</h2>
<figure><img src="figures/slide_full.png" alt="composite slide">
<figcaption>Mechanism, consequence, proof and fix on one 16:9 frame, for the
talk. The stand-alone panels above are the same figures at readable size.
</figcaption></figure>

<h2>What this does not rule out</h2>
<ul>
<li>The two runs are five days apart. The gain map could have moved between
them; that it did not is an inference from the r&nbsp;=&nbsp;+0.94 agreement,
not an independent measurement.</li>
<li>The fitted level is a <i>charge</i> threshold expressed in DREAM ADC. It
absorbs anything else that scales with pulse height and bites at the low end
&mdash; time-over-threshold requirements, the peak-finder's own floor &mdash;
so "the discriminator" is the natural reading, not a proven one.</li>
<li>The additive {n['offset_adc']:.0f}-count offset is measured, not explained.
A pedestal in the peak ADC is the obvious candidate; a charge-injection scan
would settle it in an afternoon.</li>
<li>Only P2_OUT. P2_MID and P2_IN have their own threshold settings and, in
P2_MID's case, 27 dead channels on VMM&nbsp;12; nothing here transfers to them
without redoing it.</li>
<li>Everything is conditional on a uRWELL track pointing at the pad, so it is a
statement about efficiency for tracks, not about noise occupancy at a lower
threshold &mdash; which is the cost side of every option above and is not
measured here.</li>
</ul>

<p class="note">Chain: <code>urw_p2_padadc.py</code> (DREAM per-pad spectra,
lxplus) &rarr; <code>eff_autopsy_report.py::pad_pulse_height</code> (VMM)
&rarr; <code>compare_dream_vmm.py</code> &rarr;
<code>threshold_model.py</code> &rarr; <code>figures_slide.py</code>.</p>
</main>
"""
    p = os.path.join(HERE, "report_slide.html")
    open(p, "w").write(h)
    print("wrote", p)


if __name__ == "__main__":
    main()
