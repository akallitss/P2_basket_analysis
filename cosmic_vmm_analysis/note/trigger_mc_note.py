#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Write the CERN-site note "Cosmic stand in 3D: what the trigger geometry alone predicts" from the products of
trigger_acceptance_mc.py (geometry_3d.json, fit.json, PNGs). Self-contained HTML (no CDN): the 3D view is a small
canvas renderer written here (orbit / zoom / pan, preset views, layer toggles).

Usage: trigger_mc_note.py MC_DIR NOTE_HTML   (figures are copied to NOTE_HTML without .html + '_files/')
"""
import json
import os
import shutil
import sys

mc_dir, note = sys.argv[1], sys.argv[2]
G = json.load(open(os.path.join(mc_dir, 'geometry_3d.json')))
R = json.load(open(os.path.join(mc_dir, 'fit.json')))
slug = os.path.basename(note)[:-5]
fdir = os.path.join(os.path.dirname(note), slug + '_files')
os.makedirs(fdir, exist_ok=True)
FIGS = ['kernel_K+0.png', 'projections_all.png', 'density_geometry_x_tag_eff.png', 'density_geometry.png',
        'projections_K+0.png']
for f in FIGS:
    shutil.copy(os.path.join(mc_dir, f), os.path.join(fdir, f))

NM = ['det3 P2_IN', 'det2 P2_MID', 'det1 P2_OUT']
geo, eff = R['geometry'], R['geometry_x_tag_eff']
off = geo['offset_bottom_minus_top_mm']
offe = eff['offset_bottom_minus_top_mm']
import math
ang = R['band_angle_deg']
ua = (math.cos(math.radians(ang)), math.sin(math.radians(ang)))
along = off[0] * ua[0] + off[1] * ua[1]
along_e = offe[0] * ua[0] + offe[1] * ua[1]
across = -off[0] * ua[1] + off[1] * ua[0]
tilt = math.degrees(math.atan(abs(along) / (R['zfaces_mm'][2] + 25 - (R['zfaces_mm'][1] + R['zfaces_mm'][0]) / 2)))
pp = R['per_probe']
K0 = R['K0']
f = lambda v, n=2: f'{v:.{n}f}'


def tr(cells, th=False):
    t = 'th' if th else 'td'
    return '<tr>' + ''.join(f'<{t}>{c}</{t}>' for c in cells) + '</tr>'


rows_abs = '\n'.join(tr([nm, f(pp[nm]['readout_crossing_per_trigger_mc']), f(pp[nm]['clean_cluster_per_trigger_data']),
                         f(pp[nm]['implied_clean_cluster_eff']), f(pp[nm]['tags_per_trigger_mc']),
                         f(pp[nm]['tags_per_trigger_data']), f(pp[nm]['implied_tag_pair_eff_product']),
                         f"{K0['n_readout_tags_eff1_at_measured_triggers'][nm]:.0f}"]) for nm in NM)
gfit = geo['fit']
rows_geo = '\n'.join([
    tr(['top paddle pair', f"({gfit['top_x']:.0f}, {gfit['top_y']:.0f})", f"{gfit['top_angle']:.1f}&deg;",
        f"{R['zfaces_mm'][1]:.0f} &hellip; {R['zfaces_mm'][0]:.0f}", 'x-y fitted; z from the logbook (top at z&nbsp;=&nbsp;0)']),
    tr(['det3 P2_IN', '&mdash;', '&mdash;', f"{R['z_chambers_mm'][NM[0]]:.0f}", 'analysis code (SPS spacing, &asymp;&nbsp;31&nbsp;cm)']),
    tr(['det2 P2_MID', '&mdash;', '&mdash;', f"{R['z_chambers_mm'][NM[1]]:.0f}", 'idem']),
    tr(['det1 P2_OUT', '&mdash;', '&mdash;', f"{R['z_chambers_mm'][NM[2]]:.0f}", 'idem']),
    tr(['bottom paddle', f"({gfit['bottom_x']:.0f}, {gfit['bottom_y']:.0f})", f"{gfit['bottom_angle']:.1f}&deg;",
        f"{R['zfaces_mm'][2]:.0f} &hellip; {R['zfaces_mm'][3]:.0f}", 'x-y fitted; z = det1 + 31&nbsp;cm (you, 4 Nov)']),
])

G3 = json.dumps(G, separators=(',', ':'))

HTML = r'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<!-- ===== hub metadata (read by scripts/build_hub.py) ===================== -->
<title>Cosmic stand in 3D: what the trigger geometry alone predicts</title>
<meta name="description" content="Interactive 3D model of the lab cosmic stand (2+1 scintillator paddles, three P2 chambers read on c3-c5) and a toy Monte Carlo: where the scintillator-confirmed tags fall if every chamber were 100 % efficient, compared with the wk2 efficiency-map statistics.">
<meta name="category" content="Detector R&amp;D">
<meta name="keywords" content="VMM,cosmics,P2,lab,scintillator,acceptance,Monte Carlo,geometry,3D">
<meta name="date" content="2026-11-04">
<meta name="updated" content="2026-11-04">
<meta name="visibility" content="unlisted">
<!-- ======================================================================= -->
<link rel="stylesheet" href="../../assets/site.css">
<style>
  .note-body { max-width: none; }
  .note-body > p, .note-body > ul, .note-body > ol, .note-body > h2, .note-body > h3, .note-body > .callout { max-width: 72ch; }
  .viewer { position: relative; border: 1px solid var(--border); border-radius: var(--radius); background: var(--bg-elev);
            overflow: hidden; touch-action: none; }
  .viewer canvas { display: block; width: 100%; height: 560px; cursor: grab; }
  .viewer canvas:active { cursor: grabbing; }
  .vbar { display: flex; flex-wrap: wrap; gap: 6px; margin: 10px 0 4px; align-items: center; }
  .vbar button, .vbar label { font: inherit; font-size: .85rem; border: 1px solid var(--border); background: var(--bg-elev);
            color: var(--fg); border-radius: 6px; padding: 4px 10px; cursor: pointer; }
  .vbar button:hover { border-color: var(--accent); }
  .vbar label { display: inline-flex; gap: 5px; align-items: center; }
  .vhelp { position: absolute; left: 10px; bottom: 8px; font-size: .78rem; color: var(--fg-muted); pointer-events: none; }
  .legend3d { position: absolute; right: 10px; top: 8px; font-size: .78rem; color: var(--fg-muted); line-height: 1.5;
            background: color-mix(in srgb, var(--bg-elev) 85%, transparent); padding: 4px 8px; border-radius: 6px; }
  .sw { display: inline-block; width: 10px; height: 10px; border-radius: 2px; margin-right: 5px; vertical-align: -1px; }
  figure img { max-width: 100%; height: auto; border: 1px solid var(--border); border-radius: 6px; }
  td, th { font-variant-numeric: tabular-nums; }
  @media (max-width: 640px) { .viewer canvas { height: 420px; } }
</style>
</head>
<body>
<header class="site-header"><div class="wrap">
  <a class="brand" href="../../">Alexandra Kallitsopoulou</a>
  <nav class="nav"><a href="../../">Home</a><a href="../../notes/">Notes</a></nav>
</div></header>

<main class="wrap">
  <article>
    <header class="note-header">
      <h1>Cosmic stand in 3D: what the trigger geometry alone predicts</h1>
      <p class="muted">Lab cosmic stand at Saclay, scintillator-referenced wk2 scan (runs 181&ndash;208). Question: is the
<em>statistics</em> pattern of the efficiency maps (kernel radius in <code>maps_K+0.png</code>) a detector effect or the
geometry of the trigger with respect to the cosmic-muon distribution? Written 4 Nov 2026.</p>
      <div class="meta-row">
        <span class="pill">Detector R&amp;D</span>
        <time datetime="2026-11-04">2026-11-04</time>
        <span class="faint">VMM,cosmics,P2,lab,scintillator,acceptance,Monte Carlo,3D</span>
      </div>
    </header>

    <div class="note-body">

<div class="callout ok">
<p><strong>Short answer.</strong> For det3 and det2 the statistics pattern is geometry. For det1 it is not, and the cause is in the other two chambers, not in det1.</p>
<ul>
<li><strong>The probe's efficiency never enters.</strong> A tag for probe P needs the trigger and a clean cluster in each of the two <em>other</em> chambers. So the bottom row of the maps (kernel radius = tag density) cannot show the probe's own efficiency.</li>
<li><strong>det3 and det2: geometry.</strong> With every chamber 100&nbsp;% efficient, the simulation reproduces the measured pattern: densest in the middle of the read-out band, thinning towards its ends and long edges, empty in the corners (section 4a). A muon predicted near an end of the band in one chamber usually misses the band in one of the other two. The paddles are bigger than the band, so the trigger only shifts this a little. The fit puts the bottom paddle about __ALONG__&nbsp;mm lower-left of the top pair along the band.</li>
<li><strong>det1: not geometry.</strong> Geometry predicts det1 looks like the other two, with the dense region upper-right. The data have it lower-left, with the upper-right end almost empty. Since det1's own response doesn't enter, this comes from how det3 and det2 tag or position the track.
  <ul>
  <li>The tag-chamber efficiencies don't explain it. Every VMM runs the whole length of the band, so no single weak chip can make a gradient along the band.</li>
  <li>The likeliest cause is det3's measured position, which enters det1's prediction with weight &minus;1 (det1&nbsp;=&nbsp;2&middot;det2&nbsp;&minus;&nbsp;det3). This is still to check (section 5).</li>
  </ul></li>
</ul>
</div>

<h2>1. The stand in 3D</h2>
<p>Drag to rotate, scroll or pinch to zoom, right-drag (or shift-drag) to pan. Pad frame of the P2 BASKET Gerber map
(x, y in mm, fan apex at the origin side); z points <em>down</em>, from the top paddles to the bottom paddle. Filled pads
= the 384 read-out pads per chamber (Gerber connectors 3&ndash;5); faint = the rest of the 1280-pad fan, not read out.
Tracks are a random sample of simulated triggered muons that cross at least two read-out areas.</p>
<div class="vbar">
  <button data-view="iso">3D</button><button data-view="top">plan view (as the maps)</button>
  <button data-view="along">side, looking along the band</button><button data-view="across">side, looking across the band</button>
  <label><input type="checkbox" id="t-tracks" checked> muons</label>
  <label><input type="checkbox" id="t-fan" checked> whole fan</label>
  <label><input type="checkbox" id="t-paddles" checked> paddles</label>
  <label><input type="checkbox" id="t-spin"> spin</label>
</div>
<div class="viewer" id="viewer">
  <canvas id="c3d" aria-label="3D view of the cosmic stand"></canvas>
  <div class="legend3d">
    <div><span class="sw" style="background:var(--p-paddle)"></span>scintillator paddles (200&times;600&times;50 mm)</div>
    <div><span class="sw" style="background:var(--p-det3)"></span>det3 P2_IN&nbsp; <span class="sw" style="background:var(--p-det2)"></span>det2 P2_MID&nbsp; <span class="sw" style="background:var(--p-det1)"></span>det1 P2_OUT</div>
    <div><span class="sw" style="background:var(--p-trk3)"></span>muon through all three read-out areas&nbsp; <span class="sw" style="background:var(--p-trk2)"></span>through two</div>
  </div>
  <div class="vhelp">drag: rotate &middot; wheel/pinch: zoom &middot; right/shift-drag: pan &middot; double-click: reset</div>
</div>

<h2>2. Geometry used</h2>
<div style="overflow-x:auto"><table>
<thead>__GEO_HEAD__</thead>
<tbody>
__ROWS_GEO__
</tbody></table></div>
<ul>
<li><strong>Paddles:</strong> 200&times;600&times;50&nbsp;mm (your numbers). The two top paddles lie exactly on top of each other, so the "top 0:41 OR top 0:45" trigger has the footprint of one paddle. A paddle counts as crossed when the track is inside its rectangle on its top or bottom face.</li>
<li><strong>Chambers:</strong> same footprint, no lateral offset (the straight-track residuals gave det2 within about 9&nbsp;mm and rotations below 0.4&deg;). The read-out area is the 384 pads per chamber of the scint_efficiency json, which are Gerber connectors 3&ndash;5. The band centre is at (__BCX__, __BCY__)&nbsp;mm and its long axis is at __BANG__&deg; to x.</li>
<li><strong>Not recorded anywhere:</strong> where the paddles sit in x&ndash;y. You said they are roughly aligned with the read-out connectors. I started there and fitted both paddle centres and angles to the measured tag positions (section 3).</li>
</ul>

<h2>3. The simulation</h2>
<ul>
<li><strong>Muons:</strong> I(&theta;)&nbsp;&prop;&nbsp;cos&sup2;&theta;, with 1&nbsp;/cm&sup2;/min through a horizontal plane (the usual sea-level number). __NTHROWN__ muons are thrown uniformly on z&nbsp;=&nbsp;0 over 1.5&times;1.4&nbsp;m and propagated as straight lines.</li>
<li><strong>Tag for probe P</strong> (as <code>scint_efficiency.py</code>): the trigger fires, and the muon crosses a read-out pad in each of the two other chambers. Each chamber is 100&nbsp;% efficient. The predicted point is the line through the two measured points, smeared by 3.4&nbsp;mm (pad pitch / &radic;12) to stand in for the cluster centroid.</li>
<li><strong>Fit of the paddle positions:</strong> a binned Poisson shape likelihood of the tag points on each probe's read-out footprint, all three probes together, with the MC normalised to the data per probe. It starts from a grid over the top and bottom displacement along the band, then Nelder&ndash;Mead.</li>
<li><strong>Two models:</strong>
  <ul>
  <li>(a) geometry only;</li>
  <li>(b) geometry &times; the measured local efficiency of the two <em>tag</em> chambers, taken from their own probe maps (hit within 30&nbsp;mm, 200 nearest tags). This shows how much of the tag pattern their inefficiency can explain.</li>
  </ul></li>
</ul>
<p><strong>Fit results:</strong></p>
<ul>
<li>Bottom paddle relative to the top pair: (__OFFX__, __OFFY__)&nbsp;mm with geometry only, and (__OFFEX__, __OFFEY__)&nbsp;mm with model (b). Along the band that is __ALONG__ and __ALONGE__&nbsp;mm toward the lower-left; across the band, __ACROSS__&nbsp;mm.</li>
<li>Selected muons are therefore tilted by about __TILT__&deg; along the band on average.</li>
<li>The statistical errors from the likelihood curvature are 1&ndash;2&nbsp;mm and are meaningless here. The real uncertainty is systematic. Fitting each probe's along-band distribution alone gives a bottom-paddle offset of &minus;150 to &minus;200&nbsp;mm for det2, around &minus;50&nbsp;mm for det1, and a weak preference near &minus;250&nbsp;mm for det3. The det1 residual pulls the combined fit.</li>
<li>Paddles centred on the band with no offset is clearly worse: &Delta;(&minus;ln&nbsp;L)&nbsp;=&nbsp;__DNLL__.</li>
</ul>
<p class="muted">Code: <code>cosmic_vmm_analysis/trigger_acceptance_mc.py</code> (P2_basket_analysis, branch
<code>worktree-trigger-acceptance-mc</code>); products in
<code>data/LAB_Measurements/VMM-cosmics/vmm_threshold_scan/analysis_wk2/trigger_acceptance_mc/</code>.</p>

<h2>4. Results</h2>
<h3>4a. The statistics map, data against geometry</h3>
<figure><img src="__F__/kernel_K+0.png" alt="kernel radius maps: data and MC">
<figcaption>
<p>Kernel radius of <code>scint_eff_study.py</code> (the radius holding 150 tags, the bottom row of <code>maps_K+0.png</code>) for K&nbsp;=&nbsp;0 (runs 187+196). Rows:</p>
<ol>
<li>data;</li>
<li>MC geometry only, scaled to the same number of read-out-area tags;</li>
<li>MC with the tag-chamber efficiencies;</li>
<li>the number of tags 100&nbsp;% efficient chambers would give in the same 5.6&nbsp;h at the measured trigger rate: __NFULL3__ / __NFULL2__ / __NFULL1__ for det3 / det2 / det1, against 1890 / 2333 / 1428 measured.</li>
</ol>
<p>For det3 and det2, rows 2 and 3 reproduce the data: same dense region, same empty corners, similar radii. That pattern is geometry. For det1 they do not: geometry puts the densest tags upper-right (around x&nbsp;&asymp;&nbsp;450&nbsp;mm), while the data have them lower-left (x&nbsp;&asymp;&nbsp;250&nbsp;mm), with the upper-right end nearly empty. Row 4 is what your plot would look like with perfect chambers: the same shape, at roughly half the kernel radius.</p>
</figcaption></figure>

<h3>4b. Projections on the band axes</h3>
<figure><img src="__F__/projections_all.png" alt="tag projections along and across the band">
<figcaption>
<p>All wk2 runs pooled. Curves:</p>
<ul>
<li>grey dashed: paddles centred on the band;</li>
<li>black: fitted offset;</li>
<li>red: fitted offset with the tag-chamber efficiencies.</li>
</ul>
<p>Across the band, the width of the band sets the distributions and they match well. Along the band, det3 matches and det2 is close. det1 has too many tags at &minus;250 to &minus;100&nbsp;mm and too few at 0 to +200&nbsp;mm. No paddle position fixes that without spoiling det3 and det2.</p>
</figcaption></figure>

<h3>4c. Where data and geometry differ</h3>
<figure><img src="__F__/density_geometry_x_tag_eff.png" alt="2D tag density, data, MC and ratio">
<figcaption>Tag density in 20&nbsp;mm cells: data, MC (geometry &times; tag-chamber efficiency) and their ratio. det3 is
flat within the statistics. det2 has a mild gradient. det1 has a clear excess at its lower-left end and a deficit at
its upper-right end. Geometry-only version: <a href="__F__/density_geometry.png">density_geometry.png</a>.</figcaption></figure>

<h3>4d. Absolute numbers (side result)</h3>
<p><strong>Trigger rate:</strong> the paddle telescope alone would trigger at __TRIGMC__&nbsp;Hz with 1&nbsp;/cm&sup2;/min. The measured rate is __TRIGDATA__&nbsp;Hz after the spark/trip veto, a ratio of __TRIGRATIO__.</p>
<p>That is plausible: the flux in the building is below the open-sky value, and the paddles and discriminators are not fully efficient.</p>
<p><strong>Coverage:</strong> because the paddles are larger than the read-out band, only __ALL3MC__ of the triggers cross all three read-out areas, even with perfect chambers. In the data __ALL3DATA__ of the triggers have a clean cluster in all three chambers.</p>
<div style="overflow-x:auto"><table>
<thead>__ABS_HEAD__</thead>
<tbody>
__ROWS_ABS__
</tbody></table></div>
<p class="muted">The "clean-cluster efficiency" is the fraction of triggered muons that cross a chamber's read-out area and leave a clean cluster there (1&ndash;6 pads within 15&nbsp;mm, the tag condition). It is lower than the operating efficiency (det3 0.73, det2 0.67, det1 0.80 at K&nbsp;=&nbsp;0) because it also needs the cluster to be clean, and because the trigger includes some non-muon or non-crossing coincidences. That is why the measured tag count is about a third of the geometric expectation.</p>

<h2>5. What I think the det1 residual is, and how to check</h2>
<p>Each VMM reads a strip that runs the <em>whole length</em> of the band; VMM 2 to 7 sit side by side across it. So no single weak or noisy chip, such as det2 VMM&nbsp;9 or det3 VMM&nbsp;4/7, can make a gradient <em>along</em> the band. That is also why model (b) changes almost nothing.</p>
<p>The pattern instead follows how much each chamber's measured position enters the prediction:</p>
<ul>
<li>det3's tags come from det2 and det1 and look perfect;</li>
<li>det2's prediction takes det3 with weight &frac12; and shows a mild residual;</li>
<li>det1's prediction takes det3 with weight &minus;1 (det1&nbsp;=&nbsp;2&middot;det2&nbsp;&minus;&nbsp;det3) and shows the strongest one.</li>
</ul>
<p>The common factor is det3's measured position. Two things could cause it:</p>
<ul>
<li>a local position error in det3 along the band, for example in the connectors where the P2_IN order corrections were applied (c3 top and c5 top set to "lin");</li>
<li>a z of det3 or det2 that differs from the assumed 320 / 630&nbsp;mm. This is less likely: it stretches the extrapolation rather than shifting it.</li>
</ul>
<p>If it is a position error, the det1 <em>efficiency</em> map (top row of the maps) is affected too. A mis-predicted tag counts as a det1 miss, so det1's efficiency in its upper-right is the least trustworthy number on the page.</p>
<p>Checks, cheapest first:</p>
<ol>
<li>Measure the real chamber heights on the stand.</li>
<li>Plot the det1 probe residual (prediction minus nearest hit) against the along-band position. A det3 position error shows up as a residual shift that changes along the band.</li>
<li>Repeat the half-sample order scan for P2_IN c3 and c5 with the scintillator-confirmed tracks.</li>
<li>Note down where the paddles actually sit on the stand. That would turn the fitted offset into a test.</li>
</ol>

<h2>6. Caveats</h2>
<ul>
<li>The paddle positions are fitted, not measured. Their along-band offset is the least constrained number (&minus;50 to &minus;200&nbsp;mm depending on the probe).</li>
<li>No muon energy spectrum or multiple scattering: muons are straight lines. Showers and accidental coincidences are not simulated.</li>
<li>A "clean cluster" is a single 100&nbsp;% efficient hit. Cluster size, noise clusters and masked pads are not modelled. They only change the absolute numbers and the tag-chamber factors, not the geometric shape.</li>
<li>The data tails far outside the read-out area (a few % of the tags, mis-predicted) are not in the fit.</li>
</ul>
    </div>
  </article>
</main>

<footer class="site-footer"><div class="wrap">
  <span>Alexandra Kallitsopoulou</span>
  <span><a href="../../notes/">&larr; all notes</a></span>
</div></footer>

<script>
(function () {
const G = __G3__;
const cv = document.getElementById('c3d');
const ctx = cv.getContext('2d');
const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
// palette (theme-aware via CSS variables, defaults here)
const root = document.documentElement.style;
const dark = matchMedia('(prefers-color-scheme: dark)');
function setPalette() {
  const d = dark.matches;
  root.setProperty('--p-paddle', d ? '#d9a54a' : '#c98a1b');
  root.setProperty('--p-det3', d ? '#6fa8ff' : '#2463c9');
  root.setProperty('--p-det2', d ? '#ff8a5c' : '#d9531e');
  root.setProperty('--p-det1', d ? '#4fd1a5' : '#17966b');
  root.setProperty('--p-trk3', d ? '#f0f0f0' : '#202226');
  root.setProperty('--p-trk2', d ? '#8a8f98' : '#9aa0aa');
}
setPalette();
const NM = ['det3 P2_IN', 'det2 P2_MID', 'det1 P2_OUT'];
const DCOL = ['--p-det3', '--p-det2', '--p-det1'];
const ZC = NM.map(n => G.chambers[n].z);
// scene centre = read-out band centre, mid height
let cx = 0, cy = 0, nn = 0;
for (const n of NM) { const p = G.readout_pads[n]; for (let i = 0; i < p.cx.length; i++) { cx += p.cx[i]; cy += p.cy[i]; nn++; } }
cx /= nn; cy /= nn;
const zMid = (G.paddles.top.z[0] + G.paddles.bottom.z[1]) / 2;
const bandAng = __BANGRAD__;

function padQuad(x, y, w, h, a, z) {
  const r = a * Math.PI / 180, c = Math.cos(r), s = Math.sin(r);
  const pts = [[-w / 2, -h / 2], [w / 2, -h / 2], [w / 2, h / 2], [-w / 2, h / 2]];
  return pts.map(([u, v]) => [x + u * c - v * s, y + u * s + v * c, z]);
}
// geometry primitives in world coords (x, y, z down)
const prims = [];   // {pts:[[x,y,z]..], fill, alpha, stroke, layer}
NM.forEach((n, k) => {
  const p = G.readout_pads[n];
  for (let i = 0; i < p.cx.length; i++)
    prims.push({pts: padQuad(p.cx[i], p.cy[i], p.w[i] * 0.92, p.h[i] * 0.92, p.a[i], ZC[k]), fill: DCOL[k], alpha: 0.78, layer: 'ro'});
  const q = G.all_pads;
  for (let i = 0; i < q.cx.length; i++)
    prims.push({pts: padQuad(q.cx[i], q.cy[i], q.w[i] * 0.92, q.h[i] * 0.92, q.a[i], ZC[k] + 0.5), fill: DCOL[k], alpha: 0.12, layer: 'fan'});
});
function boxFaces(corners, z0, z1) {
  const top = corners.map(c => [c[0], c[1], z0]), bot = corners.map(c => [c[0], c[1], z1]);
  const f = [top, bot];
  for (let i = 0; i < 4; i++) { const j = (i + 1) % 4; f.push([top[i], top[j], bot[j], bot[i]]); }
  return f;
}
const pz = G.paddles;
const topZ = [pz.top.z[0], pz.top.z[0] + 50, pz.top.z[1]];
for (const [z0, z1] of [[topZ[0], topZ[1]], [topZ[1], topZ[2]]])
  for (const f of boxFaces(pz.top.corners, z0, z1)) prims.push({pts: f, fill: '--p-paddle', alpha: 0.22, stroke: true, layer: 'pad'});
for (const f of boxFaces(pz.bottom.corners, pz.bottom.z[0], pz.bottom.z[1])) prims.push({pts: f, fill: '--p-paddle', alpha: 0.22, stroke: true, layer: 'pad'});
// tracks: split into segments for depth sorting
const zA = topZ[0] - 60, zB = pz.bottom.z[1] + 60, NSEG = 14;
const TR = G.tracks.filter(t => t.hits.every(h => h)).slice(0, 120).concat(G.tracks.filter(t => !t.hits.every(h => h)).slice(0, 40));
for (const t of TR) {
  const all3 = t.hits.every(h => h);
  for (let s = 0; s < NSEG; s++) {
    const za = zA + (zB - zA) * s / NSEG, zb = zA + (zB - zA) * (s + 1) / NSEG;
    prims.push({line: true, pts: [[t.x0 + za * t.tx, t.y0 + za * t.ty, za], [t.x0 + zb * t.tx, t.y0 + zb * t.ty, zb]],
                fill: all3 ? "--p-trk3" : "--p-trk2", alpha: all3 ? 0.35 : 0.3, layer: "trk"});
  }
}

// camera
const V0 = {yaw: -0.6, pitch: 0.42, dist: 4000, px: 0, py: 0};
let V = Object.assign({}, V0);
const VIEWS = {
  iso: {yaw: -0.6, pitch: 0.42},
  top: {yaw: 0, pitch: -Math.PI / 2 + 0.0001},      // plan view with x right, y up (as the 2D maps)
  along: {yaw: Math.PI / 2 - bandAng, pitch: 0.0},   // band axis = viewing direction
  across: {yaw: -bandAng, pitch: 0.0},               // band axis = screen horizontal
};
function project(p, W, H, f) {
  // world (x, y, z-down) -> view: up axis = -z
  let X = p[0] - cx, Y = p[1] - cy, Z = -(p[2] - zMid);
  const cyw = Math.cos(V.yaw), syw = Math.sin(V.yaw);
  let x1 = X * cyw - Y * syw, y1 = X * syw + Y * cyw;          // rotate about vertical
  const cp = Math.cos(V.pitch), sp = Math.sin(V.pitch);
  // camera looks horizontally along +y1 when pitch = 0; pitch tilts it down
  const yc = y1 * cp + Z * sp, zc = -y1 * sp + Z * cp;          // depth = yc, up = zc
  const d = V.dist + yc;
  const s = f / Math.max(d, 50);
  return [W / 2 + (x1 * s) + V.px, H / 2 - (zc * s) + V.py, d];
}
function draw() {
  const dpr = window.devicePixelRatio || 1;
  const W = cv.clientWidth, H = cv.clientHeight;
  if (cv.width !== Math.round(W * dpr)) { cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr); }
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.fillStyle = css('--bg-elev') || '#fff';
  ctx.fillRect(0, 0, W, H);
  const f = Math.min(W, H) * 1.9;
  const show = {ro: true, fan: document.getElementById('t-fan').checked, pad: document.getElementById('t-paddles').checked,
                trk: document.getElementById('t-tracks').checked};
  const list = [];
  for (const p of prims) {
    if (!show[p.layer]) continue;
    const q = p.pts.map(a => project(a, W, H, f));
    let d = 0; for (const a of q) d += a[2]; d /= q.length;
    list.push([d, p, q]);
  }
  list.sort((a, b) => b[0] - a[0]);
  const cache = {};
  const col = n => cache[n] || (cache[n] = css(n));
  for (const [, p, q] of list) {
    ctx.globalAlpha = p.alpha;
    if (p.line) {
      ctx.strokeStyle = col(p.fill); ctx.lineWidth = 1.1;
      ctx.beginPath(); ctx.moveTo(q[0][0], q[0][1]); ctx.lineTo(q[1][0], q[1][1]); ctx.stroke();
      continue;
    }
    ctx.fillStyle = col(p.fill);
    ctx.beginPath(); ctx.moveTo(q[0][0], q[0][1]);
    for (let i = 1; i < q.length; i++) ctx.lineTo(q[i][0], q[i][1]);
    ctx.closePath(); ctx.fill();
    if (p.stroke) { ctx.globalAlpha = 0.7; ctx.strokeStyle = col(p.fill); ctx.lineWidth = 1; ctx.stroke(); }
  }
  ctx.globalAlpha = 1;
  // labels
  ctx.font = '12px ' + (css('--font') || 'sans-serif');
  ctx.fillStyle = css('--fg-muted') || '#555';
  const lab = (t, p) => { const q = project(p, W, H, f); ctx.fillText(t, q[0] + 6, q[1]); };
  const bx = Math.max(...G.all_pads.cx) + 30, by = Math.min(...G.all_pads.cy);
  NM.forEach((n, k) => lab(n + '  (z = ' + ZC[k] + ' mm)', [bx, by, ZC[k]]));
  const tc = pz.top.corners, bc = pz.bottom.corners;
  lab('top paddles (2, stacked)', [tc[1][0], tc[1][1], topZ[0]]);
  lab('bottom paddle', [bc[1][0], bc[1][1], pz.bottom.z[1]]);
  // axis triad at the fan apex
  const o = [0, 0, ZC[0]];
  const ax = [[[150, 0, 0], 'x'], [[0, 150, 0], 'y'], [[0, 0, 150], 'z (down)']];
  const qo = project(o, W, H, f);
  ctx.strokeStyle = css('--fg-faint') || '#999'; ctx.lineWidth = 1.2;
  for (const [d, t] of ax) {
    const q = project([o[0] + d[0], o[1] + d[1], o[2] + d[2]], W, H, f);
    ctx.beginPath(); ctx.moveTo(qo[0], qo[1]); ctx.lineTo(q[0], q[1]); ctx.stroke();
    ctx.fillText(t, q[0] + 3, q[1] - 3);
  }
}
let raf = 0;
const redraw = () => { if (!raf) raf = requestAnimationFrame(() => { raf = 0; draw(); }); };
// interaction
let drag = null;
cv.addEventListener('contextmenu', e => e.preventDefault());
cv.addEventListener('pointerdown', e => {
  cv.setPointerCapture(e.pointerId);
  drag = {x: e.clientX, y: e.clientY, pan: e.button === 2 || e.shiftKey, id: e.pointerId};
});
const pointers = new Map();
cv.addEventListener('pointerdown', e => pointers.set(e.pointerId, [e.clientX, e.clientY]));
cv.addEventListener('pointerup', e => { pointers.delete(e.pointerId); drag = null; pinch = null; });
cv.addEventListener('pointercancel', e => { pointers.delete(e.pointerId); drag = null; pinch = null; });
let pinch = null;
cv.addEventListener('pointermove', e => {
  if (pointers.has(e.pointerId)) pointers.set(e.pointerId, [e.clientX, e.clientY]);
  if (pointers.size === 2) {
    const [a, b] = [...pointers.values()];
    const d = Math.hypot(a[0] - b[0], a[1] - b[1]);
    if (pinch) V.dist = Math.min(9000, Math.max(800, V.dist * pinch / d));
    pinch = d; redraw(); return;
  }
  if (!drag) return;
  const dx = e.clientX - drag.x, dy = e.clientY - drag.y;
  drag.x = e.clientX; drag.y = e.clientY;
  if (drag.pan) { V.px += dx; V.py += dy; }
  else { V.yaw -= dx * 0.008; V.pitch = Math.min(Math.PI / 2 - 0.0001, Math.max(-Math.PI / 2 + 0.0001, V.pitch + dy * 0.008)); }
  redraw();
});
cv.addEventListener('wheel', e => { e.preventDefault(); V.dist = Math.min(9000, Math.max(800, V.dist * Math.exp(e.deltaY * 0.001))); redraw(); },
                    {passive: false});
cv.addEventListener('dblclick', () => { V = Object.assign({}, V0); redraw(); });
document.querySelectorAll('[data-view]').forEach(b => b.addEventListener('click', () => {
  Object.assign(V, VIEWS[b.dataset.view], {px: 0, py: 0}); redraw();
}));
document.querySelectorAll('.vbar input').forEach(i => i.addEventListener('change', redraw));
let spinT = 0;
function spin(t) {
  if (document.getElementById('t-spin').checked) { if (spinT) V.yaw += (t - spinT) * 0.00025; spinT = t; draw(); }
  else spinT = 0;
  requestAnimationFrame(spin);
}
requestAnimationFrame(spin);
dark.addEventListener('change', () => { setPalette(); redraw(); });
new ResizeObserver(redraw).observe(cv);
draw();
})();
</script>
</body>
</html>
'''
sub = {
    '__G3__': G3, '__F__': slug + '_files',
    '__GEO_HEAD__': tr(['element', 'centre (x, y) [mm]', 'long axis to x', 'z [mm] (down)', 'source'], th=True),
    '__ROWS_GEO__': rows_geo,
    '__ABS_HEAD__': tr(['probe', 'triggers crossing its read-out (MC)', 'clean cluster per trigger (data)',
                        '&rArr; clean-cluster efficiency', 'tags per trigger, MC eff = 1', 'tags per trigger, data',
                        '&rArr; tag-pair efficiency product', 'K0 read-out-area tags at eff = 1'], th=True),
    '__ROWS_ABS__': rows_abs,
    '__BCX__': f"{R['band_centre_mm'][0]:.0f}", '__BCY__': f"{R['band_centre_mm'][1]:.0f}",
    '__BANG__': f'{ang:.1f}', '__BANGRAD__': f'{math.radians(ang):.5f}',
    '__NTHROWN__': f"{R['params']['n']:.0e}".replace('e+0', '&times;10<sup>') + '</sup>',
    '__OFFX__': f'{off[0]:+.0f}', '__OFFY__': f'{off[1]:+.0f}', '__OFFEX__': f'{offe[0]:+.0f}', '__OFFEY__': f'{offe[1]:+.0f}',
    '__ALONG__': f'{abs(along):.0f}', '__ALONGE__': f'{abs(along_e):.0f}', '__ACROSS__': f'{across:+.0f}',
    '__TILT__': f'{tilt:.0f}', '__DNLL__': f"{geo['nll_aligned'] - geo['nll']:.0f}",
    '__NFULL3__': f"{K0['n_readout_tags_eff1_at_measured_triggers'][NM[0]]:.0f}",
    '__NFULL2__': f"{K0['n_readout_tags_eff1_at_measured_triggers'][NM[1]]:.0f}",
    '__NFULL1__': f"{K0['n_readout_tags_eff1_at_measured_triggers'][NM[2]]:.0f}",
    '__TRIGMC__': f"{R['trigger_rate_hz_eff1']:.2f}", '__TRIGDATA__': f"{R['data']['trigger_rate_hz']:.2f}",
    '__TRIGRATIO__': f"{R['data']['trigger_rate_hz'] / R['trigger_rate_hz_eff1']:.2f}",
    '__ALL3MC__': f"{R['all_three_per_trigger_mc']:.2f}", '__ALL3DATA__': f"{R['all_three_per_trigger_data']:.2f}",
}
for k, v in sub.items():
    HTML = HTML.replace(k, v)
assert '__' not in HTML.replace('__init__', ''), [l for l in HTML.split('\n') if '__' in l][:3]
open(note, 'w').write(HTML)
print('wrote', note, f'{len(HTML) / 1024:.0f} kB')
