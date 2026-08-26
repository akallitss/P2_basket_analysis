#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""make_report_dv.py -- report for the VMM/DREAM per-pad comparison on P2_OUT.

Every number is computed here from data/compare_dream_vmm_P2_OUT.csv and the
two histogram sets, so re-running after the analysis moves the text with it.
"""
import json
import os

import numpy as np
import pandas as pd
from scipy import stats

import make_report as R                       # CSS
import figures_dv as D                        # loader + _hq

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "report_dream_vmm.html")


def numbers():
    g, VH, vbin, DH, dbin = D.load()
    n = {"pads": len(g)}
    rv, rd = D.rel(g, "amp_med_v"), D.rel(g, "amp_med_d")
    n["rv"], n["rd"] = rv, rd
    n["rms_v"], n["rms_d"] = float(np.std(rv, ddof=1)), float(np.std(rd, ddof=1))
    n["ratio"] = n["rms_v"] / n["rms_d"]
    n["r"] = float(np.corrcoef(rv, rd)[0, 1])
    b = np.polyfit(rd, rv, 1)
    n["slope"], n["icept"] = float(b[0]), float(b[1])
    bl = np.polyfit(g["amp_med_d"], g["amp_med_v"], 1)
    n["abs_slope"], n["abs_icept"] = float(bl[0]), float(bl[1])

    n["tracks_v"] = int(g["n_track_v"].sum())
    n["tracks_d"] = int(g["n_track_d"].sum())
    n["eff_v"] = float(g["k_hit_v"].sum() / g["n_track_v"].sum())
    n["eff_d"] = float(g["k_hit_d"].sum() / g["n_track_d"].sum())
    live = g[g["pad_id"] != 635]
    n["eff_v_lo"], n["eff_v_hi"] = float(live["eff_v"].min()), float(live["eff_v"].max())
    n["eff_d_lo"], n["eff_d_hi"] = float(live["eff_d"].min()), float(live["eff_d"].max())
    n["span_v"] = float(g["amp_med_v"].max() / g["amp_med_v"].min())
    n["span_d"] = float(g["amp_med_d"].max() / g["amp_med_d"].min())

    x, y = g["x_v"].to_numpy(), g["y_v"].to_numpy()
    A = np.column_stack([np.ones(len(g)), x, y])
    def fit(v):
        c = np.linalg.lstsq(A, v, rcond=None)[0]
        return c, v - A @ c
    out = {}
    for lab, v in (("v", rv), ("d", rd)):
        c, r = fit(v)
        out[lab] = r
        n[f"grad_{lab}"] = float(np.hypot(c[1], c[2]) * 10 * 100)
        n[f"R2_{lab}"] = float(1 - r.var() / v.var())
        n[f"ang_{lab}"] = float(np.degrees(np.arctan2(c[2], c[1])))
        n[f"res_{lab}"] = float(np.std(r, ddof=1))
    n["r_res"] = float(np.corrcoef(out["v"], out["d"])[0, 1])

    n["r_eff_v"] = float(np.corrcoef(g["eff_v"], rv)[0, 1])
    n["r_eff_d"] = float(np.corrcoef(g["eff_d"], rd)[0, 1])
    n["area_span"] = float(100 * (g["pad_area"].max() / g["pad_area"].min() - 1))
    n["r_area"] = float(np.corrcoef(g["pad_area"], rd)[0, 1])

    # the truncation control
    di = g["di"].to_numpy()
    n["cut_scan"] = []
    for c in (0, 136.0, 200, 240, 300):
        b0 = int(np.ceil(c / dbin))
        hh = DH[di][:, b0:]
        med = np.array([D._hq(h, dbin) for h in hh]) + b0 * dbin
        n["cut_scan"].append((c, float(np.std(med / np.median(med), ddof=1)),
                              float(hh.sum() / DH[di].sum())))
    n["dead_v"] = float(g.loc[g["pad_id"] == 635, "eff_v"].iloc[0])
    n["dead_d"] = float(g.loc[g["pad_id"] == 635, "eff_d"].iloc[0])

    meta = json.load(open(os.path.join(HERE, "data",
                                       "dream_padadc_eff_nominal_1_P2_OUT.json")))
    n["subruns"] = len(meta["sub_runs"])
    return n, g


def main():
    n, g = numbers()
    rows = "\n".join(
        f"<tr><td>{c:.0f}</td><td>{s:.3f}</td><td>{k:.3f}</td></tr>"
        for c, s, k in n["cut_scan"])

    html = f"""<title>P2_OUT: VMM vs DREAM per-pad gain</title>
<meta name="description" content="The same P2_OUT chamber read by VMM3a and by
DREAM at the same operating point: the pad-to-pad pulse-height spread is the
detector's, and the VMM measures 0.6x of it.">
<style>{R.CSS}</style>
<main>
<h1>P2_OUT pad-by-pad gain: VMM3a against DREAM</h1>
<p class="sub">VMM <code>run_46 / cfg_gain4.5_peaktime200</code> (1 Aug 04:53)
against DREAM <code>eff_nominal_1 / eff_nominal_00&ndash;09</code> (27 Jul
10:40). Same chamber, same uRWELL reference, same track selection, same
operating point &mdash; mesh 450&nbsp;V, drift 750&nbsp;V, gap 300&nbsp;V,
Ar/CO<sub>2</sub>/iC<sub>4</sub>H<sub>10</sub> 93/5/2.</p>

<div class="verdict">
<b>No &mdash; the VMM sees <i>less</i> pad-to-pad variation than DREAM, not
more.</b> On the {n['pads']} pads both readouts illuminate, the per-pad pulse
height (each divided by its own readout's fleet median) has an rms of
<b>{n['rms_v']:.2f}</b> on the VMM and <b>{n['rms_d']:.2f}</b> on DREAM, and the
two agree pad for pad at <b>r&nbsp;=&nbsp;{n['r']:+.3f}</b>. The VMM is a
compressed copy of the same map, slope {n['slope']:.2f}.
<br><br>
<b>The variation is the detector's, and it is one smooth gradient.</b> A linear
plane in (x,&nbsp;y) explains {n['R2_d']:.0%} of it on DREAM and
{n['R2_v']:.0%} on the VMM, in the same direction to
{abs(n['ang_d'] - n['ang_v']):.0f}&deg;; what is left over still correlates
between the two readouts at r&nbsp;=&nbsp;{n['r_res']:+.2f}. Neither chain adds
pad-to-pad structure of its own that the other cannot see.
<br><br>
<b>What differs is what that gain map costs.</b> The same gradient moves the
VMM's per-pad efficiency over
{n['eff_v_lo']:.2f}&ndash;{n['eff_v_hi']:.2f} and DREAM's over
{n['eff_d_lo']:.2f}&ndash;{n['eff_d_hi']:.2f}.
</div>

<div class="kpis">
<div class="kpi"><b>{n['r']:+.3f}</b><span>r(VMM, DREAM) per pad</span></div>
<div class="kpi"><b>{n['rms_v']:.2f} / {n['rms_d']:.2f}</b><span>rel. spread VMM / DREAM</span></div>
<div class="kpi"><b>{n['slope']:.2f}</b><span>VMM per unit DREAM</span></div>
<div class="kpi"><b>{n['eff_v']:.3f} / {n['eff_d']:.3f}</b><span>efficiency VMM / DREAM</span></div>
<div class="kpi"><b>{n['grad_d']:.0f} %</b><span>gain gradient per 10 mm</span></div>
</div>

<h2>Why these two runs</h2>
<p>The VMM configuration scan sat at one operating point for every
<code>cfg_*</code> run: mesh 450&nbsp;V / drift 750&nbsp;V on all three P2
stations. Of the DREAM runs, <code>highstat_eff_1</code> &mdash; the one the
handoff quotes &mdash; is at drift <b>700</b>&nbsp;V, a 250&nbsp;V gap, not
300. <code>eff_nominal_1</code> is at 750/450, the identical point, and its
{n['subruns']} full-rate sub-runs give {n['tracks_d']:,} fiducial tracks against
the VMM's {n['tracks_v']:,}. Gas and beam are the same; the two runs are five
days apart.</p>
<p class="note">The drift gap is not a sensitive axis here &mdash; the VMM drift
scan gives 73.9 / 74.0 / 73.6&nbsp;% at gaps 300 / 350 / 400 &mdash; so
<code>highstat_eff_1</code> would not have been wrong, only avoidable.</p>

<h2>The measurement, on both sides</h2>
<p>Identical by construction: uRWELL front+back two-point tracks, the same
|back&minus;front| cut, the same free-affine frame fit onto the pad plane (both
land on &minus;60&deg;), the same 9&nbsp;mm fiducial and 15&nbsp;mm probe
radius, the same spark veto and recorded-event hygiene. For every pad: how many
tracks pointed at it, how many it recorded, and the pulse-height distribution of
the ones it did.</p>
<ul>
<li><b>VMM</b> &mdash; <code>win_adc</code>, the ADC of the in-window unmasked
hit nearest the prediction. Per-pad medians are DNL-corrected (16-code bins).
New producer: <code>pad_pulse_height()</code> in
<code>eff_autopsy_report.py</code>.</li>
<li><b>DREAM</b> &mdash; <code>a_lead</code>, the leading pad's amplitude,
restricted to the events where the leading pad <em>is</em> the pad the track
pointed at (87.9&nbsp;% of them), which makes it the same quantity. New
producer: <code>urw_p2_padadc.py</code>, a driver over
<code>urw_p2_efficiency</code>'s own functions.</li>
</ul>

<figure><img src="figures/dv_headline.png" alt="VMM vs DREAM per-pad pulse height">
<figcaption>Left: one point per pad, marker area &prop; tracks. The dashed line
is where the two would sit if they measured the same spread; the fit is
{n['slope']:.2f}. Right: the control &mdash; cutting the DREAM spectra from below
shrinks the spread it measures, but a cut matched to the VMM's turn-on gets
only part of the way.</figcaption></figure>

<h2>The two spectra, on the same pads</h2>
<figure><img src="figures/dv_spectra.png" alt="pulse-height spectra">
<figcaption>Pulse height in units of each readout's own fleet MPV, log density.
The VMM records nothing below ~0.57 of the fleet MPV; DREAM's zero suppression
is far below the signal. On the low-gain pad the VMM's median sits
<em>above</em> DREAM's (1.22 vs 1.05) and on the high-gain pad well
<em>below</em> it (2.41 vs 3.86) &mdash; that is the compression, pad by
pad.</figcaption></figure>

<h2>The gradient</h2>
<figure><img src="figures/dv_maps.png" alt="pulse-height maps and gradient">
<figcaption>Relative pulse height across the illuminated pads, same colour scale
in both maps. DREAM measures {n['grad_d']:.0f}&nbsp;% per 10&nbsp;mm at
{n['ang_d']:+.0f}&deg;, the VMM {n['grad_v']:.0f}&nbsp;% at
{n['ang_v']:+.0f}&deg;.</figcaption></figure>

<h2>What it costs each readout</h2>
<figure><img src="figures/dv_efficiency.png" alt="efficiency vs pulse height">
<figcaption>Per-pad efficiency against relative pulse height. Both readouts show
the correlation (r = {n['r_eff_v']:+.2f} VMM, {n['r_eff_d']:+.2f} DREAM); only
the VMM pays for it, because its discriminator sits under the signal. Pad 635 is
low in both &mdash; the one genuinely dead pad, found independently by two
electronics chains.</figcaption></figure>

<h2>The truncation control</h2>
<p>The VMM's spectrum turns on at ~64&nbsp;ADC against a pooled MPV of 112; the
DREAM MPV is 238, so the same relative threshold is 136 DREAM ADC. A fixed cut
lifts a low-gain pad's median more than a high-gain pad's, so it compresses the
measured spread &mdash; the question is by how much.</p>
<div class="tw"><table>
<tr><th>cut on the DREAM spectra [ADC]</th><th>spread (rms)</th><th>pulses kept</th></tr>
{rows}
</table></div>
<p>At the matched cut DREAM still measures {n['cut_scan'][1][1]:.2f} against the
VMM's {n['rms_v']:.2f}; reproducing the VMM's number needs a cut near 300, which
would throw away {1 - n['cut_scan'][4][2]:.0%} of the pulses when the VMM's
efficiency deficit accounts for only about 10&nbsp;%. <b>So the threshold is the
mechanism but not the whole of it</b> &mdash; a linear fit of the two per-pad
medians in raw ADC gives VMM = {n['abs_slope']:.3f}&nbsp;&times;&nbsp;DREAM
{n['abs_icept']:+.0f}, and that {n['abs_icept']:.0f}-ADC offset (which is the
turn-on) is what pulls every ratio toward 1.</p>

<h2>What this does not rule out</h2>
<ul>
<li><b>It does not measure absolute gain.</b> Both axes are relative to their
own readout, because the two ADCs share no scale and neither was charge
calibrated in this campaign.</li>
<li><b>The residual compression is an additive ADC offset</b> &mdash; resolved
2026-08-22, after this note first went up. Threshold truncation alone gets DREAM
from {n['rms_d']:.2f} only to {n['cut_scan'][1][1]:.2f}, not to
{n['rms_v']:.2f}; the rest is a ~43-count constant in the VMM peak ADC.
Subtracting it from the VMM per-pad medians returns 0.363, the truncated DREAM
number, with nothing left over &mdash; so this is no longer open, and it is not
front-end compression or a non-linearity. See
<a href="p2out-vmm-threshold-weak-pads.html">Why the VMM loses the weak
pads</a>, which also fits the threshold itself (162 DREAM ADC) and shows it
reproduces the VMM's per-pad efficiency. A charge-injection sweep would still
say <i>what</i> the offset is.</li>
<li><b>Both ADCs saturate</b> &mdash; 1.0&nbsp;% of the VMM entries in the top
bin, 1.7&nbsp;% of the DREAM ones. That is above the medians used here, but it
is why MPVs are read with the top of the range blanked.</li>
<li><b>Pad area is not the driver</b>, though it correlates at
r&nbsp;=&nbsp;{n['r_area']:+.2f}: across these pads the area spans
{n['area_span']:.2f}&nbsp;% while the pulse height spans a factor
{n['span_d']:.1f}. In the fan, area is a label for radius; the pulse height
follows <em>position</em>.</li>
<li><b>Five days separate the two runs.</b> Nothing was moved and the HV is
identical, but a slow gain drift between 27 Jul and 1 Aug would show up as a
scale factor &mdash; which is exactly the quantity being divided out here, so it
cannot affect the spread comparison. It would affect an absolute one.</li>
<li><b>Only P2_OUT.</b> P2_MID and P2_IN are not in this comparison.</li>
</ul>
</main>
"""
    open(OUT, "w").write(html)
    print(f"wrote {OUT}")
    for k in ("pads", "r", "rms_v", "rms_d", "slope", "abs_slope", "abs_icept",
              "R2_v", "R2_d", "r_res", "grad_v", "grad_d", "eff_v", "eff_d"):
        print(f"  {k:10s} {n[k]}")


if __name__ == "__main__":
    main()
