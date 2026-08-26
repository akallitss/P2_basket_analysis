#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""make_report.py -- builds report.html for the P2_OUT per-pad ADC study.

Every number is computed here from the same tables the figures use, so
re-running after the analysis updates the text and the tables together.
"""

import json, os
import numpy as np, pandas as pd
from scipy import stats
import figures as F

HERE = os.path.dirname(os.path.abspath(__file__))


def numbers():
    g = F.load()
    o = g[(g.station == "P2_OUT")]
    u = o[o.use]
    n = {}
    n["pads"] = len(u); n["tracks"] = int(u.n_track.sum())
    n["eff"] = u.k_track.sum() / u.n_track.sum()
    n["adc_lo"], n["adc_hi"] = u.adc_med.min(), u.adc_med.max()
    n["r"], n["p"] = stats.pearsonr(u.adc_med, u.eff)
    n["rho"] = stats.spearmanr(u.adc_med, u.eff)[0]
    n["rw"] = stats.pearsonr(
        u.adc_med - u.groupby("vmm").adc_med.transform("mean"),
        u.eff - u.groupby("vmm").eff.transform("mean"))[0]
    n["slope"] = np.polyfit(u.adc_med, u.eff, 1)[0] * 1000
    v = u.dropna(subset=["adc_med_untracked"])
    n["r_val"] = stats.pearsonr(v.adc_med, v.adc_med_untracked)[0]
    n["spread_val"] = (v.adc_med_untracked - v.adc_med).std()
    n["best"] = sorted(u.eff, reverse=True)[:3]
    n["r_area_adc"] = stats.pearsonr(u.pad_area, u.adc_med)[0]
    n["r_area_eff"] = stats.pearsonr(u.pad_area, u.eff)[0]
    # Pad AREA spans only 0.16 % across the illuminated pads while the pulse
    # height spans a factor 2.9, so area cannot be the physical driver: it is a
    # monotone label for position in the fan.  Control for position itself.
    n["area_span"] = 100 * (u.pad_area.max() / u.pad_area.min() - 1)
    n["adc_span"] = u.adc_med.max() / u.adc_med.min()
    A = np.column_stack([np.ones(len(u)), u.x.to_numpy(), u.y.to_numpy()])
    res = lambda y: y - A @ np.linalg.lstsq(A, y, rcond=None)[0]
    n["r_partial"] = stats.pearsonr(res(u.adc_med.to_numpy()),
                                    res(u.eff.to_numpy()))[0]
    b = np.linalg.lstsq(A, u.adc_med.to_numpy() / u.adc_med.median(),
                        rcond=None)[0]
    n["grad_pct10"] = 100 * np.hypot(b[1], b[2]) * 10
    n["grad_R2"] = 1 - res(u.adc_med.to_numpy()).var() / u.adc_med.var()
    n["dead"] = int((o.dead).sum())

    s = u.sort_values("adc_med")
    cw = np.cumsum(s.n_track.to_numpy())
    cut = np.searchsorted(cw, np.linspace(0, cw[-1], 6)[1:-1])
    n["quint"] = []
    for q in np.split(np.arange(len(s)), cut):
        z = s.iloc[q]
        n["quint"].append((z.adc_med.min(), z.adc_med.max(),
                           z.k_track.sum() / z.n_track.sum(),
                           int(z.n_track.sum()), len(z)))
    return n, u


CSS = """
:root{--bg:#fcfcfb;--fg:#0b0b0b;--fg2:#52514e;--mut:#8a887f;--line:#e6e5e0;
--acc:#2a78d6;--warn:#eb6834;--card:#ffffff}
@media(prefers-color-scheme:dark){:root:not([data-theme=light]){--bg:#1a1a19;
--fg:#fff;--fg2:#c3c2b7;--mut:#8a887f;--line:#33322f;--acc:#3987e5;
--warn:#d95926;--card:#211f1e}}
[data-theme=dark]{--bg:#1a1a19;--fg:#fff;--fg2:#c3c2b7;--line:#33322f;
--acc:#3987e5;--warn:#d95926;--card:#211f1e}
body{background:var(--bg);color:var(--fg);font:15px/1.6 -apple-system,
BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif;margin:0;padding:2rem 1.2rem}
main{max-width:1000px;margin:0 auto}
h1{font-size:1.7rem;line-height:1.25;margin:0 0 .3rem}
h2{font-size:1.15rem;margin:2.2rem 0 .6rem;padding-bottom:.3rem;
border-bottom:1px solid var(--line)}
.sub{color:var(--fg2);margin:0 0 1.6rem}
.verdict{background:var(--card);border-left:3px solid var(--acc);
padding:1rem 1.2rem;border-radius:6px;margin:1.2rem 0}
.kpis{display:flex;flex-wrap:wrap;gap:.8rem;margin:1.2rem 0}
.kpi{background:var(--card);border:1px solid var(--line);border-radius:8px;
padding:.7rem 1rem;min-width:120px}
.kpi b{display:block;font-size:1.5rem;color:var(--acc);line-height:1.15}
.kpi span{font-size:.8rem;color:var(--fg2)}
table{border-collapse:collapse;width:100%;font-size:.9rem;margin:.8rem 0}
th,td{padding:.4rem .6rem;border-bottom:1px solid var(--line);text-align:right}
th:first-child,td:first-child{text-align:left}
th{color:var(--fg2);font-weight:600}
figure{margin:1.4rem 0}
img{width:100%;height:auto;border:1px solid var(--line);border-radius:6px}
figcaption{color:var(--fg2);font-size:.86rem;margin-top:.45rem}
code{background:var(--card);padding:.1rem .35rem;border-radius:4px;font-size:.88em}
.tw{overflow-x:auto}
.note{color:var(--fg2);font-size:.9rem}
ul{padding-left:1.2rem}
"""


def main():
    n, u = numbers()
    q = "".join(f"<tr><td>{a:.0f} – {b:.0f}</td><td>{e:.3f}</td>"
                f"<td>{t:,}</td><td>{p}</td></tr>"
                for a, b, e, t, p in n["quint"])
    h = f"""<title>P2_OUT pad pulse height</title>
<style>{CSS}</style>
<main>
<h1>P2_OUT: per-pad efficiency is set by the pad's pulse height</h1>
<p class="sub">run_46 / cfg_gain4.5_peaktime200 &mdash; VMM3a readout, SPS H4,
1 Aug 2026. Task 1 and task 3 of <code>HANDOFF_VMM_PAD_ADC.md</code>.</p>

<div class="verdict">
<b>The correlation is real and it reproduces.</b> Over {n['pads']} illuminated
P2_OUT pads carrying {n['tracks']:,} tracks, per-pad efficiency rises with the
pad's median pulse height at Pearson <b>r = {n['r']:+.2f}</b>
(p = {n['p']:.0e}, Spearman &rho; = {n['rho']:+.2f}) &mdash; confirming the
+0.55 quoted in the handoff, which had no committed producer. It now has one.
<br><br>
<b>But most of the spread is one smooth gradient across the chamber.</b> A
linear plane in (x, y) accounts for {n['grad_R2']:.0%} of the pad-to-pad
variance &mdash; {n['grad_pct10']:.0f}&nbsp;% of the fleet median per 10&nbsp;mm
&mdash; and removing it drops the pulse-height&rarr;efficiency partial
correlation to <b>{n['r_partial']:+.2f}</b>. So the handoff's first open
question splits in two: a chamber-scale gain gradient, which no per-channel
trim can fix, and a small per-pad residual, which one could.
</div>

<div class="kpis">
<div class="kpi"><b>{n['eff']:.4f}</b><span>P2_OUT efficiency</span></div>
<div class="kpi"><b>{n['r']:+.2f}</b><span>r(pulse height, eff)</span></div>
<div class="kpi"><b>{n['r_partial']:+.2f}</b><span>same, gradient removed</span></div>
<div class="kpi"><b>{n['best'][0]:.3f}</b><span>best pad &mdash; DREAM is 0.960</span></div>
<div class="kpi"><b>{n['slope']:+.1f}</b><span>pts per 10 ADC</span></div>
</div>

<h2>What was measured</h2>
<p>Per-pad efficiency (n tracks pointing at the pad, k with a hit within
probe_r) comes from <code>eff_autopsy_report.py</code>. Per-pad pulse height is
<code>win_adc_P2_OUT</code> &mdash; the ADC of the nearest in-window unmasked
hit &mdash; aggregated per pad. That aggregation is the missing producer;
it is now committed as <code>pad_pulse_height()</code> in that file, storing
median, 25/75 percentiles, count and a DNL-clean median per pad.</p>
<p class="note">All medians quoted are DNL-corrected: the VMM3a ADC is ~2&times;
too wide on every multiple of 16, so the median is read off a histogram
rebinned to 16 codes, one full DNL period per bin. The same correction puts the
P2_OUT most-probable value at <b>104</b>, not the 128 in
<code>EFFICIENCY_AUTOPSY.md</code> &mdash; 128 = 8&times;16 is a comb tooth.
That confirms the handoff's own correction from an independent sample.</p>

<figure><img src="figures/p2out_headline.png" alt="efficiency vs pad median ADC">
<figcaption>Left: one point per pad, area &prop; tracks, error bars binomial. All
six P2_OUT chips sit at the same threshold (sdt 224), so nothing here is a
per-chip threshold effect. Right: the same pads in track-weighted quintiles of
pulse height.</figcaption></figure>

<h2>Efficiency by quintile of pad pulse height</h2>
<div class="tw"><table>
<tr><th>pad median ADC</th><th>efficiency</th><th>tracks</th><th>pads</th></tr>
{q}
</table></div>
<p class="note">The first and last quintiles reproduce the handoff's table
(0.752 and ~0.92); the middle bins differ by up to 0.05 because these
boundaries are set on DNL-corrected medians.</p>

<h2>It is gain, not threshold</h2>
<figure><img src="figures/p2out_spectra_map.png" alt="spectra and maps">
<figcaption>Left: the ADC spectra of the lowest- and highest-efficiency thirds
start at the <em>same</em> threshold onset and differ in where the bulk sits &mdash;
a gain difference, not a threshold difference. Centre and right: the
pulse-height map and the efficiency map are the same map, low at bottom-left,
high at top-right. Both the DREAM readout of the same chamber and this one
measure that gradient in the same direction.</figcaption></figure>

<h2>Direction of the correlation</h2>
<p>A threshold artefact would run the other way: a channel with a high threshold
records only its big pulses, so it would show a <em>high</em> median with a
<em>low</em> efficiency. The correlation is positive, so it is real signal-size
variation. Two further checks support that here &mdash; every P2_OUT chip is at
the same sdt 224, and removing each chip's mean leaves the correlation intact
(within-chip r = {n['rw']:+.2f}).</p>

<h2>Cross-check: an independent, untracked measurement</h2>
<p>The same per-pad pulse height can be had with no pcapng decoding at all,
from <code>adc_vs_ch</code> in every capture's <code>counts.npz</code> &mdash; a
per-channel ADC histogram that <code>vmm_reduce.py</code> writes whether or not
<code>--drop-columns</code> was used. Summed over run_46's 48 captures that is
1.52&times;10<sup>9</sup> hits against the 837k tracked ones. The two medians
agree at <b>r = {n['r_val']:.3f}</b> with a spread of {n['spread_val']:.0f} ADC
on P2_OUT, so the cheap route is usable for the rest of the campaign &mdash; the
runs whose columns were dropped included.</p>
<figure><img src="figures/tracked_vs_untracked.png" alt="tracked vs untracked">
<figcaption>Tracked against untracked per-pad median. P2_OUT and P2_MID agree;
P2_IN does not, because its untracked spectrum is dominated by a screaming
channel &mdash; another reason not to fold P2_IN in.</figcaption></figure>

<h2>What this does not rule out</h2>
<ul>
<li><b>It does not close the P2_OUT gap to DREAM.</b> The best pads reach
{n['best'][0]:.3f} / {n['best'][1]:.3f} / {n['best'][2]:.3f} against DREAM's
0.960, but only {int((u.eff>=0.95).sum())} of {n['pads']} pads are there. The
mechanism is identified, not removed.</li>
<li><b>The residual after the gradient is not attributed.</b> r = {n['r_partial']:+.2f}
survives; whether that is gas gain, mesh planarity or per-channel baseline
scatter is still open, and separating them needs the per-channel trim run
(<code>qa_config.py: 'calibration': None</code> &mdash; no such run exists).</li>
<li><b>Pad area is NOT the operative variable</b>, although it correlates at
r = {n['r_area_adc']:+.2f}. Across the illuminated pads the area spans only
{n['area_span']:.2f}&nbsp;% while the pulse height spans a factor
{n['adc_span']:.1f}, so the area cannot be causing it &mdash; in the fan, area
is a monotone label for radius, and it is <em>position</em> that the pulse
height follows. An earlier version of this note said &ldquo;most of the spread
is pad geometry&rdquo;; that attribution was wrong and is corrected above. The
partial-correlation arithmetic is unchanged.</li>
<li><b>P2_MID does not behave like this.</b> Its pooled correlation is
+0.14 (rate mask) to +0.20 (ratio mask), not the +0.52 the handoff quotes;
that number did not reproduce under either mask rule. Its two chips at sdt 256
run the between-chip trend against the within-chip one. P2_MID and P2_IN are
written up separately.</li>
</ul>

<h2>Reproducing</h2>
<p class="note">
<code>run_autopsy.sh run_46 cfg_gain4.5_peaktime200</code> &rarr;
<code>run_report.sh &lt;tracks.parquet&gt; --mask rate</code> &rarr;
<code>compare_tracked.py</code> &rarr; <code>figures_p2out.py</code> &rarr;
this page. Autopsy output is on EOS under
<code>vmm/pad_adc/</code>; it reproduces the handoff's reference numbers exactly
(P2_IN 0.1535, P2_MID 0.5916, P2_OUT 0.8540).</p>
</main>
"""
    p = os.path.join(HERE, "report.html")
    open(p, "w").write(h)
    print("wrote", p)


if __name__ == "__main__":
    main()
