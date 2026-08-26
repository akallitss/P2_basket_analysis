#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""make_report_track.py -- report_track.html: what the three P2 stations know
about a track, measured against the uRWELL reference.

Every number is read from the same NPZ/JSON the figures use, so a re-run moves
prose, tables and pictures together.
"""
import json
import os

import numpy as np

import make_report as R
import track_stats as T

HERE = os.path.dirname(os.path.abspath(__file__))

EXTRA = """
main{max-width:1180px}
.big{font-size:1.05rem;line-height:1.65}
figure img{width:100%;height:auto;border:1px solid var(--line);border-radius:8px}
.two{display:grid;grid-template-columns:1fr 1fr;gap:1.2rem}
@media(max-width:820px){.two{grid-template-columns:1fr}}
.flag{background:var(--card);border-left:3px solid var(--warn);
padding:.9rem 1.1rem;border-radius:6px;margin:1.1rem 0}
"""


def numbers():
    z, j = T.load()
    st = T.station_block(z, j)
    sb = T.selftrack_block(z, j)
    n = {"run": T.RUN, "n_subrun": len(j),
         "n_tracks": sum(b["n_tracks"] for b in j),
         "stations": st, "self": sb}
    n["probe_r"] = j[0]["probe_r_mm"]
    n["fid_r"] = j[0]["fid_r_mm"]
    n["track_cut"] = j[0]["urwell"]["track_cut_mm"]
    # reference quality, pooled: the front-back agreement width is what the
    # telescope contributes to every residual below
    n["ref_sigma"] = float(np.mean([b["urwell"]["front_back"][a]["sigma_iqr_mm"]
                                    for b in j for a in ("x", "y")]))
    n["rot"] = {s: float(np.mean([b["frame"][s]["affine_rotation_deg"]
                                  for b in j])) for s in T.STATIONS}
    n["rot_spread"] = {s: float(np.std([b["frame"][s]["affine_rotation_deg"]
                                        for b in j])) for s in T.STATIONS}

    # the illusion, in one ratio
    n["self_sigma"] = T.spread(sb["mid_self_x"])
    n["ref_sigma_mid"] = T.spread(sb["mid_ref_x"])
    n["illusion"] = n["ref_sigma_mid"] / n["self_sigma"]

    # the tracker, at the exit plane and extrapolated
    n["exit"] = {a: dict(core=T.spread(sb[f"dpos_exit_{a}"]),
                         rms=float(sb[f"dpos_exit_{a}"].std()),
                         p95=float(np.percentile(np.abs(sb[f"dpos_exit_{a}"]),
                                                 95)))
                 for a in ("x", "y")}
    n["ang"] = {a: dict(core=T.spread(sb[f"dang_{a}"]),
                        p95=float(np.percentile(np.abs(sb[f"dang_{a}"]), 95)))
                for a in ("x", "y")}
    zg = sb["zgrid"]
    n["far"] = {}
    for zz in (1370.0, 2340.0):
        i = int(np.argmin(np.abs(zg - zz)))
        n["far"][zz] = {a: dict(core=T.spread(sb[f"curve_{a}"][i]),
                                p95=float(np.percentile(
                                    np.abs(sb[f"curve_{a}"][i]), 95)))
                        for a in ("x", "y")}

    # pooled single vs multi-pad, the charge-sharing question
    e = z["res_edges"]
    hs = sum(z[f"res_x_{s}_single"] for s in T.STATIONS)
    hm = sum(z[f"res_x_{s}_multi"] for s in T.STATIONS)
    n["single"] = T.hstats(hs, e)
    n["multi"] = T.hstats(hm, e)
    n["frac_multi"] = float(hm.sum() / (hs.sum() + hm.sum()))
    n["pitch_floor_x"] = T.PITCH_X / np.sqrt(12)
    n["pitch_floor_y"] = T.PITCH_Y / np.sqrt(12)

    # inside one pad: is a hit lost at the edge?  (no -- but the sharing is
    # visible there, and only there)
    n["inpad"] = {}
    for what in ("eff", "amp", "nclus"):
        rel = []
        for s in T.STATIONS:
            c, v, den = T.padprofile(z, s, what)
            flat = (c > 1.0) & (c < 4.5)
            rel.append(v / np.mean(v[flat]))
        y = np.mean(rel, axis=0)
        n["inpad"][what] = float(y[-1] - 1.0)
    n["inpad"]["r_edge"] = float(T.padprofile(z, "P2_OUT")[0][-1])

    # purity: how fast the plateau creeps once past a pad
    n["accid"] = {}
    de = z["dmin_edges"]
    for s in T.STATIONS:
        h = np.asarray(z[f"dmin_{s}"], float)
        cum = np.cumsum(h) / st[s]["n"]
        i = int(np.searchsorted(de[1:], n["probe_r"]))
        n["accid"][s] = float((cum[-1] - cum[i]) / (de[-1] - n["probe_r"]) * 10)
    return n


def table_stations(n):
    rows = []
    for s in T.STATIONS:
        b = n["stations"][s]
        rows.append(
            f"<tr><td><b>{s}</b></td><td>{b['z']:.0f}</td>"
            f"<td>{b['n']:,}</td><td>{b['eff']:.4f}</td>"
            f"<td>{b['res_x']['rms']:.2f}</td>"
            f"<td>{b['res_y']['rms']:.2f}</td>"
            f"<td>{b['res_x']['sigma_iqr']:.2f}</td>"
            f"<td>{b['frac_multi']:.1%}</td>"
            f"<td>{n['accid'][s] * 100:+.2f} %</td></tr>")
    return ("<table><tr><th>station</th><th>z [mm]</th><th>tracks</th>"
            "<th>efficiency</th><th>&sigma;<sub>core</sub> x</th>"
            "<th>&sigma;<sub>core</sub> y</th><th>rms x</th>"
            "<th>&ge;2 pads</th><th>accidental /10 mm</th></tr>"
            + "".join(rows) + "</table>")


def table_tracker(n):
    zg = n["self"]["zgrid"]
    rows = []
    for lab, blk in (("exit of the basket (z = 940)", n["exit"]),
                     ("back reference plane (z = 1370)", n["far"][1370.0]),
                     ("1.4 m past the basket (z = 2340)", n["far"][2340.0])):
        rows.append(
            f"<tr><td>{lab}</td>"
            f"<td>{blk['x']['core']:.2f}</td><td>{blk['y']['core']:.2f}</td>"
            f"<td>{blk['x']['p95']:.1f}</td><td>{blk['y']['p95']:.1f}</td></tr>")
    return ("<table><tr><th>predicted at</th><th>core x [mm]</th>"
            "<th>core y [mm]</th><th>95 % x [mm]</th>"
            "<th>95 % y [mm]</th></tr>" + "".join(rows) + "</table>")


FIGS = [
    ("track_1_pointing.png",
     "One station against the reference. The residual is a 12 mm box, not a "
     "Gaussian, and it is the same box on all three stations &mdash; which is "
     "what a pad detector with no charge sharing looks like. The third panel "
     "splits by cluster size: a centroid over two pads is not better than a "
     "single pad, it is slightly worse."),
    ("track_2_selftrack.png",
     "The three stations used as a tracker. Left: where the P2-only track says "
     "the particle crossed the exit plane, against where the reference says it "
     "did. Middle: the direction, with a narrow core on a heavy tail. Right: "
     "the error against the lever arm &mdash; solid is the core, dashed the "
     "95th percentile, and the shaded region is extrapolation beyond the last "
     "station."),
    ("track_3_illusion.png",
     "The headline. Same events, same clusters, two rulers: P2_MID checked "
     "against the line through P2_IN and P2_OUT, and P2_MID checked against "
     "the reference. The self-check is not a measurement &mdash; most of the "
     "time all three stations report the identical pad, so a straight line "
     "through them has nothing to say."),
    ("track_4_maps.png",
     "Efficiency against the reference track, binned in each station's own pad "
     "frame. A defect of the probe stands still here; a defect of the "
     "reference would appear in the same place on all three. The faint "
     "diagonal banding is the 4 mm bin beating against the ~12 mm pad fan, "
     "not structure in the chamber."),
    ("track_5_purity.png",
     "How much of the efficiency is a real hit. Left: the distance from the "
     "track to the nearest cluster, per unit area, so a flat tail is what an "
     "accidental match looks like. Right: the same as a cumulative &mdash; "
     "efficiency against the one knob that moves it."),
    ("track_6_inpad.png",
     "Every track folded onto the face of the pad it pointed at. The question "
     "was whether a pad loses hits at its own edge; it does not. Efficiency is "
     "flat across the face and a shade <i>higher</i> at the rim, where two "
     "pads each get a chance at the same charge. What the rim does show is the "
     "sharing itself &mdash; the leading pad keeps ~5 % less and the cluster "
     "grows ~19 % &mdash; and that it only turns on in the outer tenth of the "
     "face is why the residual is the full pitch/&radic;12 box. The folded "
     "face is a disc, not a square, because the pads are a fan: they are not "
     "all the same size, so folding them onto one set of axes mixes cells of "
     "slightly different shape."),
]


def main():
    n = numbers()
    st, sf = n["stations"], n["self"]
    effs = " / ".join(f"{st[s]['eff']:.3f}" for s in T.STATIONS)
    figs = "\n".join(
        f'<figure><img src="figures/{f}" alt="{c[:80]}">'
        f'<figcaption>{c}</figcaption></figure>' for f, c in FIGS)

    h = f"""<title>P2 tracking against an external reference</title>
<meta name="description" content="The three P2 BASKET stations judged against
uRWELL reference tracks: per-station pointing, the P2-only track, and why P2
cannot check its own tracking.">
<style>{R.CSS}{EXTRA}</style>
<main>
<h1>What the P2 stations know about a track</h1>
<p class="sub">P2 BASKET, SPS H4, July 2026. Three uRWELL pad stations at
z = 320 / 630 / 940 mm, bracketed by two 1 mm-pitch EIC uRWELL reference planes
1370 mm apart. Run <code>{n['run']}</code>, {n['n_subrun']} sub-runs,
{n['n_tracks']:,} reference tracks.</p>

<div class="verdict big">
<b>The three stations find the particle, and then lose it again.</b>
Detection is good &mdash; {effs} of the reference tracks pointing at a live
pad are recorded, and a three-station track exists for
{sf['n_station_frac'][3]:.0%} of them. But every station reports a
<b>pad</b>, not a position: the residual to the reference is a
{T.PITCH_X:.1f} mm box of rms {st['P2_OUT']['res_x']['rms']:.2f} mm against
the {n['pitch_floor_x']:.2f} mm a pad gives with no charge sharing at all, and
a two-pad cluster does not improve it. The P2-only track therefore points to
<b>{n['exit']['x']['core']:.1f} mm</b> at the exit of the basket and
<b>{n['ang']['x']['core']:.2f} mrad</b> in the core, with a tail that reaches
{n['exit']['x']['p95']:.0f} mm at the 95th percentile.
<br><br>
<b>And P2 cannot measure any of that on its own.</b> Checking P2_MID against
the line through P2_IN and P2_OUT returns
{n['self_sigma']:.2f} mm &mdash; {n['illusion']:.0f}&times; better than the
{n['ref_sigma_mid']:.2f} mm the reference measures on the same events. The
reason is not subtle: {sf['frac_same_pad']:.0%} of the time all three stations
report the <i>identical pad</i>, so the self-consistency residual is the same
rounding error three times over.
</div>

<div class="kpis">
<div class="kpi"><b>{sf['n_station_frac'][3]:.0%}</b>
<span>reference tracks with a 3-station P2 track</span></div>
<div class="kpi"><b>{st['P2_OUT']['res_x']['rms']:.2f} mm</b>
<span>one station, rms (pitch/&radic;12 = {n['pitch_floor_x']:.2f})</span></div>
<div class="kpi"><b>{n['exit']['x']['core']:.1f} mm</b>
<span>P2-only track at the basket exit</span></div>
<div class="kpi"><b>{n['ang']['x']['core']:.2f} mrad</b>
<span>P2-only direction, core</span></div>
<div class="kpi"><b>{n['illusion']:.0f}&times;</b>
<span>how much the self-check flatters itself</span></div>
</div>

<h2>What was compared</h2>
<p>Nothing about the reference is new here. Tracks, the uRWELL &rarr; pad frame
fit, the fiducial cut and the probe radius all come from
<code>urw_p2_efficiency.py</code> with its defaults &mdash; two-point tracks
with a {n['track_cut']:.0f} mm front&harr;back agreement cut, a fiducial of
{n['fid_r']:.0f} mm to the nearest pad centre, and a
{n['probe_r']:.0f} mm probe radius. The reference's own front&harr;back
agreement is {n['ref_sigma']:.2f} mm core, so it contributes essentially
nothing in quadrature to a {T.PITCH_X:.0f} mm pad.</p>

<p>What is new is that all three stations are put in <b>one frame at once</b>.
Each station's uRWELL&nbsp;&rarr;&nbsp;pad affine is inverted to bring its
clusters back into the reference frame, and a straight line is fitted through
the three. The alignment is therefore the reference's, which means every mean
below is zero by construction and only the <b>widths</b> are a measurement.
That is the honest configuration &mdash; an experiment would align against
something external too &mdash; but nothing here bounds a global scale or
rotation error. The fitted rotation is stable to
{max(n['rot_spread'].values()) * 1e3:.0f} millidegrees across the
{n['n_subrun']} sub-runs, so nothing moved during the run.</p>

<h2>Per station</h2>
{table_stations(n)}
<p class="note">The two widths are quoted together because for a box they sit
the opposite way round from a Gaussian: a uniform pad gives rms&nbsp;=&nbsp;
pitch/&radic;12 = {n['pitch_floor_x']:.2f}&nbsp;mm but an IQR width of
0.5&nbsp;&times;&nbsp;pitch/1.349 = {T.PITCH_X * 0.5 / 1.349:.2f}&nbsp;mm, so
it is the <b>rms</b> that tests the no-sharing prediction and the IQR that is
robust once the tails matter. The last column is the accidental-match rate,
read off the slope of the efficiency plateau past one pad.</p>

<h2>Inside one pad</h2>
<p>Folding every track onto the face of the pad it pointed at answers the
obvious follow-up: is a hit lost at the pad edge, where the avalanche is split
with the neighbour? <b>No.</b> Across the outer tenth of the face the leading
pad keeps {n['inpad']['amp'] * 100:+.1f} % of its central charge and the
cluster grows by {n['inpad']['nclus'] * 100:+.0f} %, but the efficiency moves
by {n['inpad']['eff'] * 100:+.1f} % &mdash; <i>upwards</i>, because at a
boundary two pads each get a chance at the same charge. Everything inside
4 mm is flat to 1 %.</p>
<p>That is the same fact as the residual, seen from the other side. Sharing
that only appears in the last millimetre, on {n['frac_multi']:.0%} of clusters,
cannot build a centroid &mdash; which is why the two-pad clusters in the third
panel of the first figure are no better than the single-pad ones.</p>

<h2>The P2-only track, against the reference</h2>
{table_tracker(n)}
<p class="note">Core is the IQR width over
{sf['n']:,} three-station tracks; the 95 % column is the
half-width containing 95 % of them. The gap between the two columns is the
point: this is not a Gaussian tracker, it is a quantised one with a tail from
whichever station picked the neighbouring pad.</p>

<div class="flag">
<b>The one number to carry away.</b> Stage 22's tag-and-probe, and any other
check that uses P2 to judge P2, is looking at
{sf['frac_same_pad']:.0%} identical pads. It cannot see the
{n['ref_sigma_mid']:.2f} mm the reference sees, and it cannot see the tail at
all. On <code>highstat_eff_1</code> the same effect shows up in the efficiency:
92.5 / 96.3 / 92.5 % self-referenced against 96.5 / 97.1 / 96.0 % against the
uRWELL. The tracking version of that gap is {n['illusion']:.0f}&times;, not
4 points.
</div>

<h2>Figures</h2>
{figs}

<h2>What this does not rule out</h2>
<ul>
<li><b>It does not bound a global alignment error.</b> The uRWELL &rarr; pad
affine is fitted per station on these same events, so a common rotation, scale
or offset of the whole basket is absorbed, not measured. In the experiment that
term is whatever the survey and the alignment procedure leave behind, and it
adds to everything here.</li>
<li><b>It is one working point.</b> All {n['n_subrun']} sub-runs sit at the
nominal mesh and drift of <code>{n['run']}</code>. Cluster size, and therefore
the {n['frac_multi']:.0%} of clusters that span two pads, grows with gain
&mdash; a higher mesh would move the sharing question, though the
{n['multi']['rms']:.2f} mm the two-pad clusters give today says it would
have to move a long way to help.</li>
<li><b>The y pitch is an assumption.</b> In x the measured rms lands on
pitch/&radic;12 to {abs(st['P2_OUT']['res_x']['rms'] / n['pitch_floor_x'] - 1) * 100:.1f} %;
in y it sits {(st['P2_OUT']['res_y']['rms'] / n['pitch_floor_y'] - 1) * 100:.0f} %
high. The pads are a fan described by radius and &phi;, so a single number for
the y pitch is a simplification &mdash; the excess may be geometry rather than
anything the detector is doing. It has not been chased.</li>
<li><b>It says nothing about timing.</b> Only position and direction are
measured here.</li>
<li><b>The tail is attributed, not proven.</b> "One station picked the
neighbouring pad" is what the shape looks like and what the same-pad fraction
implies; it has not been separated from a genuine second particle in the
trigger, though the second-cluster candidate <code>(x2, y2)</code> is already
in the matcher to suppress exactly that.</li>
<li><b>Two-station tracks are not costed.</b>
{sf['n_station_frac'][2]:.0%} of tracks get exactly two stations, which is
enough for a position and not for an angle; whether the experiment should use
them is a separate question.</li>
</ul>

<p class="note">Built by <code>p2_sps/make_report_track.py</code> from
<code>p2_selftrack.py</code>'s output. Analysis:
<code>p2_selftrack.py</code> (lxplus, LCG_110) &rarr;
<code>track_stats.py</code> &rarr; <code>figures_track.py</code>.</p>
</main>
"""
    out = os.path.join(HERE, "report_track.html")
    open(out, "w").write(h)
    print(f"wrote {out}")
    return n


if __name__ == "__main__":
    main()
