#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Write the CERN-site note "Cosmic stand in 3D: the trigger geometry, and why the tag statistics is the paddles'
efficiency" from the products of trigger_acceptance_mc.py (geometry_3d.json, fit.json, PNGs) and of paddle_study.py
(paddle_geometry.json, paddle_efficiency.json, closure.json, trigger_composition.json, PNGs). Self-contained HTML
(no CDN): the 3D view is a small canvas renderer written here (orbit / zoom / pan, preset views, layer toggles); the
paddles are drawn at their measured rectangles, with the efficient region filled.

Usage: trigger_mc_note.py MC_DIR PADDLE_STUDY_DIR NOTE_HTML   (figures are copied to NOTE_HTML without .html + '_files/')
"""
import json
import os
import shutil
import sys

mc_dir, ps_dir, note = sys.argv[1], sys.argv[2], sys.argv[3]
G = json.load(open(os.path.join(mc_dir, 'geometry_3d.json')))
R = json.load(open(os.path.join(mc_dir, 'fit.json')))
PG = json.load(open(os.path.join(ps_dir, 'paddle_geometry.json')))
PE = json.load(open(os.path.join(ps_dir, 'paddle_efficiency.json')))
CL = json.load(open(os.path.join(ps_dir, 'closure.json')))
TC = json.load(open(os.path.join(ps_dir, 'trigger_composition.json')))
slug = os.path.basename(note)[:-5]
fdir = os.path.join(os.path.dirname(note), slug + '_files')
os.makedirs(fdir, exist_ok=True)
FIGS = ['kernel_K+0.png', 'projections_all.png', 'density_geometry_x_tag_eff.png', 'density_geometry.png',
        'projections_K+0.png']
for f in FIGS:
    shutil.copy(os.path.join(mc_dir, f), os.path.join(fdir, f))
PS_FIGS = ['paddle_efficiency.png', 'paddle_maps.png', 'closure.png', 'paddle_images.png']
for f in PS_FIGS:
    shutil.copy(os.path.join(ps_dir, f), os.path.join(fdir, f))

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

# the 3D paddles: measured rectangles (extent where the efficiency > 0.1) and the efficient region (> half plateau)
bc = PG['band_centre_xy']
def band_to_xy(al, ac):
    return [bc[0] + al * ua[0] - ac * ua[1], bc[1] + al * ua[1] + ac * ua[0]]
PADZ = {'p0_45': [-100.0, -50.0], 'p0_41': [-50.0, 0.0], 'p0_44': [PG['z_bottom'] - 25, PG['z_bottom'] + 25]}
PNAME = {'p0_45': 'top paddle 45', 'p0_41': 'top paddle 41', 'p0_44': 'bottom paddle 44'}
G['paddles'] = []
for k in ('p0_45', 'p0_41', 'p0_44'):
    g = PG['paddles'][k]
    e0, e1 = g['efficient_along_mm']
    c0, c1 = g['efficient_across_mm']
    G['paddles'].append({'name': PNAME[k], 'corners': [[round(v, 1) for v in c] for c in g['corners_xy_mm']],
                         'eff_corners': [[round(v, 1) for v in band_to_xy(a, c)] for a, c in ((e0, c0), (e1, c0), (e1, c1), (e0, c1))],
                         'z': PADZ[k]})
G3 = json.dumps(G, separators=(',', ':'))
pg = PG['paddles']
pe = PE['paddles']
def padrow(k, zlab, src):
    g = pg[k]
    return tr([PNAME[k], f"({g['centre_xy_mm'][0]:.0f}, {g['centre_xy_mm'][1]:.0f})", f"{ang:.1f}&deg;",
               f"{g['length_mm']:.0f} &times; {g['width_mm']:.0f}", zlab,
               f"{pe[k]['eff_lower_left_half']:.2f} / {pe[k]['eff_upper_right_half']:.2f}", src])
rows_geo2 = '\n'.join([
    padrow('p0_45', '&minus;100 &hellip; 0 (pair; which of the two is on top is not recorded)', 'extent and efficiency from chamber tracks (section 5)'),
    padrow('p0_41', 'idem', 'idem'),
    tr(['det3 P2_IN', '(' + ', '.join(f'{v:.0f}' for v in R['band_centre_mm']) + ')', f'{ang:.1f}&deg;', 'wedge 470 long, 30&ndash;180 wide',
        f"{R['z_chambers_mm'][NM[0]]:.0f}", '&mdash;', 'read-out pads; z from the analysis code (SPS spacing, &asymp;&nbsp;31&nbsp;cm)']),
    tr(['det2 P2_MID', 'idem', 'idem', 'idem', f"{R['z_chambers_mm'][NM[1]]:.0f}", '&mdash;', 'idem']),
    tr(['det1 P2_OUT', 'idem', 'idem', 'idem', f"{R['z_chambers_mm'][NM[2]]:.0f}", '&mdash;', 'idem']),
    padrow('p0_44', f"{PG['z_bottom'] - 25:.0f} &hellip; {PG['z_bottom'] + 25:.0f} (assumed; its measured length suggests &asymp;&nbsp;1150, section 6)",
           'extent and efficiency from chamber tracks (section 5)'),
])
cl = CL['per_probe']
rows_cl = '\n'.join(tr([nm, f"{cl[nm]['data_mean_along']:+.0f}", f"{cl[nm]['mc_measured_mean_along']:+.0f}",
                         f"{cl[nm]['chi2_along_measured']:.0f} / {cl[nm]['bins_along_measured']}",
                         f"{cl[nm]['mc_uniform_mean_along']:+.0f}", f"{cl[nm]['chi2_along_uniform']:.0f} / {cl[nm]['bins_along_uniform']}"]) for nm in NM)
tc = TC['per_category']
rows_tc = '\n'.join(tr([k, f"{100 * TC['fraction'][k]:.0f}&nbsp;%", f"{tc[k]['slope_along_mean']:+.2f}",
                         f"{tc[k]['tag_median_along_det3']:+.0f}", f"{tc[k]['tag_median_along_det2']:+.0f}", f"{tc[k]['tag_median_along_det1']:+.0f}"])
                     for k in ('45 only', 'both', '41 only'))
rows_pe = '\n'.join(tr([PNAME[k], f"{PE['eff_on_straight_tracks'][k]:.2f}", f"{pe[k]['eff_lower_left_half']:.2f}", f"{pe[k]['eff_upper_right_half']:.2f}",
                         f"{pg[k]['efficient_along_mm'][0]:+.0f} &hellip; {pg[k]['efficient_along_mm'][1]:+.0f}",
                         f"{pe[k]['delay_ns_per_m_along']:+.1f} (vs {pe[k]['delay_reference']})"]) for k in ('p0_45', 'p0_41', 'p0_44'))

HTML = r'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<!-- ===== hub metadata (read by scripts/build_hub.py) ===================== -->
<title>Cosmic stand in 3D: the trigger geometry, and why the tag statistics is the paddles' efficiency</title>
<meta name="description" content="Interactive 3D model of the lab cosmic stand (2+1 scintillator paddles, three P2 chambers read on c3-c5), a toy Monte Carlo of where the scintillator-confirmed tags fall with perfect chambers, and the answer: the paddles are efficient only near their PMT end (44: 0.90 / 0.26 on its two halves), measured with the chambers as the reference; with the measured paddle maps the MC reproduces the tag distributions with no free parameter.">
<meta name="category" content="Detector R&amp;D">
<meta name="keywords" content="VMM,cosmics,P2,lab,scintillator,paddle efficiency,acceptance,Monte Carlo,geometry,3D">
<meta name="date" content="2026-11-04">
<meta name="updated" content="2026-11-05">
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
      <h1>Cosmic stand in 3D: the trigger geometry, and why the tag statistics is the paddles' efficiency</h1>
      <p class="muted">Lab cosmic stand at Saclay, scintillator-referenced wk2 scan (runs 181&ndash;208). Question: is the
<em>statistics</em> pattern of the efficiency maps (kernel radius in <code>maps_K+0.png</code>) a detector effect or the
geometry of the trigger with respect to the cosmic-muon distribution? Written 4 Nov 2026; resolved and
rewritten 5 Nov 2026 (sections 2, 5, 6 and the 3D paddles are from the resolution; sections 3&ndash;4 are the
geometry-only study as first done).</p>
      <div class="meta-row">
        <span class="pill">Detector R&amp;D</span>
        <time datetime="2026-11-04">2026-11-04</time>
        <span class="faint">VMM,cosmics,P2,lab,scintillator,paddle efficiency,Monte Carlo,3D</span>
      </div>
    </header>

    <div class="note-body">

<div class="callout ok">
<p><strong>Short answer.</strong> The statistics pattern of the maps is the trigger, but not its geometry: it is the
<em>efficiency</em> of the scintillator paddles, which is high only near each paddle's PMT end.</p>
<ul>
<li><strong>The probe's efficiency never enters.</strong> A tag for probe P needs the trigger and a clean cluster in each of the two <em>other</em> chambers. So the bottom row of the maps (kernel radius = tag density) cannot show the probe's own efficiency.</li>
<li><strong>Geometry alone gets det3 and det2 roughly right and det1 wrong</strong> (sections 3&ndash;4): with perfect chambers and uniform paddles the dense region of det1 should be upper-right; the data have it lower-left, with a sharp step halfway along the band.</li>
<li><strong>The cause is the bottom paddle 44</strong>, which is in every trigger: measured with the chambers as the reference (no paddle in the selection), it is __E44LL__ efficient on its lower-left half and __E44UR__ on its upper-right half (section 5). Det1, 31&nbsp;cm above it, sees that step unsmeared. Top paddle 41 is the mirror image (__E41LL__ / __E41UR__), so the tracks that fire it are tilted; top paddle 45 is fine (__E45LL__ / __E45UR__). The pulse timing puts the PMT of 44 at the lower-left end and those of 41 and 45 at the upper-right end: light attenuation against the discriminator threshold.</li>
<li><strong>Closure:</strong> the same Monte Carlo with the measured paddle maps instead of uniform paddles, and no free parameter, reproduces the three tag distributions (mean along-band positions __CLD3__ / __CLD2__ / __CLD1__&nbsp;mm in the data against __CLM3__ / __CLM2__ / __CLM1__&nbsp;mm; uniform paddles gave __CLU3__ / __CLU2__ / __CLU1__).</li>
<li><strong>The chambers are fine.</strong> Neither their positions, nor the det3 map, nor the z spacing is involved. The chamber efficiency numbers are unaffected; only the statistics per unit area is 3&ndash;4 times lower on the upper-right half of every chamber, worst for det1.</li>
</ul>
</div>

<h2>1. The stand in 3D</h2>
<p>Drag to rotate, scroll or pinch to zoom, right-drag (or shift-drag) to pan. Pad frame of the P2 BASKET Gerber map
(x, y in mm, fan apex at the origin side); z points <em>down</em>, from the top paddles to the bottom paddle. Filled pads
= the 384 read-out pads per chamber (Gerber connectors 3&ndash;5); faint = the rest of the 1280-pad fan, not read out.
Tracks are a random sample of simulated triggered muons that cross at least two read-out areas (uniform-paddle
trigger of section 3). The paddle rectangles are the measured ones (section 5): each is drawn at its extent, and the
part where its efficiency is above half of its plateau is filled. Paddle 41 is drawn above 45; which of the two is
physically on top is not recorded.</p>
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
    <div><span class="sw" style="background:var(--p-paddle)"></span>paddles: outline = measured extent, filled = region where the efficiency is above half its plateau</div>
    <div><span class="sw" style="background:var(--p-det3)"></span>det3 P2_IN&nbsp; <span class="sw" style="background:var(--p-det2)"></span>det2 P2_MID&nbsp; <span class="sw" style="background:var(--p-det1)"></span>det1 P2_OUT</div>
    <div><span class="sw" style="background:var(--p-trk3)"></span>muon through all three read-out areas&nbsp; <span class="sw" style="background:var(--p-trk2)"></span>through two</div>
  </div>
  <div class="vhelp">drag: rotate &middot; wheel/pinch: zoom &middot; right/shift-drag: pan &middot; double-click: reset</div>
</div>

<h2>2. Geometry</h2>
<div style="overflow-x:auto"><table>
<thead>__GEO_HEAD__</thead>
<tbody>
__ROWS_GEO__
</tbody></table></div>
<ul>
<li><strong>Paddles:</strong> 200&times;600&times;50&nbsp;mm (your numbers), long axis along the read-out band, the two top paddles on top of each other; the table gives their rectangles as measured from the chamber tracks (extent = where the efficiency exceeds 0.1; about 10&nbsp;mm resolution, so 630&ndash;640&nbsp;mm for a 600&nbsp;mm paddle). The measured bottom-paddle length of __L44__&nbsp;mm is longer than the others because its plane is extrapolated furthest: it would be 600&nbsp;mm with the paddle at z&nbsp;&asymp;&nbsp;1150 instead of the assumed 1250 (section 6).</li>
<li><strong>Chambers:</strong> same footprint, no lateral offset (the straight-track residuals give det2 within 2&nbsp;mm of the det3&ndash;det1 line in every run). The read-out area is the 384 pads per chamber of the scint_efficiency json, Gerber connectors 3&ndash;5: a wedge 470&nbsp;mm long, 30&nbsp;mm wide at its lower-left end and 180&nbsp;mm at its upper-right end. The band centre is at (__BCX__, __BCY__)&nbsp;mm and its long axis at __BANG__&deg; to x; "along" below is this axis, positive toward the wide upper-right end.</li>
<li><strong>Efficiency halves:</strong> mean efficiency of the paddle on the lower-left and on the upper-right half of its length (section 5).</li>
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
<p><strong>Fit results</strong> (kept for the record; superseded by section 5: a paddle <em>position</em> moves an edge, it cannot halve the density over half of a paddle):</p>
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

<h2>5. The cause: the paddles are efficient only near their PMT end</h2>
<p>The capture files record the paddle discriminator pulses on hybrid 0 (channels 41 and 45 = the two top paddles,
44 = the bottom one, 40 = the hardware triple coincidence; the ADC is a constant logic level, so the pulse height
itself is not available). That allows two measurements that need nothing but the data in hand.</p>

<h3>5a. The triggered sample, split by which top paddle fired</h3>
<p><code>scint_dump.py</code> writes every trigger with its paddle pulses and the clean clusters. The two top paddles
are far from equivalent: 45 fires in 93&nbsp;% of the triggers, 41 in only 52&nbsp;%, and the hardware triple
(channel 40) is exactly 41&nbsp;AND&nbsp;45&nbsp;AND&nbsp;44. Fitting a line through the clean clusters of each
trigger:</p>
<div style="overflow-x:auto"><table>
<thead>__TC_HEAD__</thead>
<tbody>
__ROWS_TC__
</tbody></table></div>
<p>Tracks that fire 45 only are vertical on average; tracks that fire 41 are tilted by 0.15&ndash;0.22 along the band.
They enter at the upper-right end, where 41 works, and leave through 44 at the lower-left end, where 44 works. This
tilted sub-population is what the paddle-position fit of section 3 was trying to absorb.
<a href="__F__/paddle_images.png">paddle_images.png</a>: the triggered tracks extrapolated to the paddle planes, and
the fraction of them with a 41 / 45 pulse.</p>

<h3>5b. Paddle efficiency with the chambers as the reference</h3>
<p><code>paddle_efficiency.py</code> finds three-chamber straight tracks in the raw hit stream with <em>no paddle in
the selection</em> (a det2 hit with det3 and det1 hits within &plusmn;150&nbsp;ns, one clean cluster each, det2
within 25&nbsp;mm of the line through det3 and det1), extrapolates each to the paddle planes and asks whether the
paddle channel has a pulse at the chamber time minus the 242&nbsp;ns latency, within &plusmn;120&nbsp;ns. The
chance level, from a window 5&nbsp;&micro;s earlier, is 0.000. __NSTRAIGHT__ such tracks in the ten wk2 runs.</p>
<figure><img src="__F__/paddle_efficiency.png" alt="paddle efficiency maps and profiles">
<figcaption>Top: efficiency of 45, 41 and 44 on their planes, in band coordinates. Middle: the profile along the band
(|across|&nbsp;&lt;&nbsp;60&nbsp;mm). Bottom: the median pulse time against position (44 and 41 relative to 45; 45
relative to the chamber time). 44 is 0.9 efficient on its lower-left half and about 0.25 on its upper-right half,
with the step between &minus;150 and &minus;30&nbsp;mm; 41 is the mirror image; 45 is 0.9 everywhere, mildly lower at
its lower-left end.</figcaption></figure>
<div style="overflow-x:auto"><table>
<thead>__PE_HEAD__</thead>
<tbody>
__ROWS_PE__
</tbody></table></div>
<p>The pulse timing says where the PMTs are. The 44 pulse arrives later toward the upper right and the 45 pulse
earlier, __D44__&nbsp;ns/m relative to each other, i.e. about 5.5&nbsp;ns/m each, the light propagation speed in a
scintillator slab; 41 and 45 are in step. So the PMT of 44 is at the <strong>lower-left</strong> end and the PMTs of
41 and 45 at the <strong>upper-right</strong> end, and each paddle is efficient near its own PMT and loses efficiency
at the far end: light attenuation against the discriminator threshold, with 45 having the best light yield or the
lowest threshold. The flat low plateau far from the PMT (0.15&ndash;0.25) is the fraction of the Landau tail still
above threshold.</p>
<figure><img src="__F__/paddle_maps.png" alt="adaptive-kernel paddle efficiency maps">
<figcaption>The same, as adaptive-kernel maps in the pad x&ndash;y frame, exactly as the chamber maps of
<code>scint_eff_study.py</code> (each point = mean of the 150 nearest tracks; kernel radius below). Solid outline =
extent (efficiency above 0.1), dashed = the region above half of the plateau. The read-out wedge is drawn for
reference. With __NSTRAIGHT__ tracks the kernel radius is 15&ndash;25&nbsp;mm over the paddles, as on the chamber
maps.</figcaption></figure>

<h3>5c. Closure</h3>
<figure><img src="__F__/closure.png" alt="tag distributions: data against the MC with the measured paddle maps">
<figcaption>Tag positions per probe, all wk2 runs: data (points), the geometric MC of section 3 with the measured
2-D paddle efficiency maps (black; every chamber 100&nbsp;% efficient, <strong>no free parameter</strong>, no
position fit), and the same MC with uniform 600&times;200 paddles centred on the band (grey dashed). The two-level
det1 distribution, the det2 shape and the across-band widths are reproduced.</figcaption></figure>
<div style="overflow-x:auto"><table>
<thead>__CL_HEAD__</thead>
<tbody>
__ROWS_CL__
</tbody></table></div>
<p>What remains (det3 below the MC at &minus;100 to 0&nbsp;mm and above it at +50 to +120; the narrower across-band
peak of det3 and det1; the tails beyond the wedge ends) is the chambers' own efficiency non-uniformity and the
mis-predicted tags, neither of which this MC models.</p>

<h3>5d. What to do</h3>
<ol>
<li><strong>Fix the trigger:</strong> lower the discriminator thresholds (or raise the HV) of 44 and 41 until their profile matches 45's; or at least turn 44 end for end so that its good half lies under the half where 41 is good. <code>paddle_efficiency.py</code> measures the result without any trigger, before and after.</li>
<li>Until then, prefer 45&nbsp;AND&nbsp;44 to the hardware triple: on chamber-only straight tracks the triple fires in __E40__, the offline (41&nbsp;OR&nbsp;45)&nbsp;AND&nbsp;44 in __ETRIG__, 45 alone in __E45__.</li>
<li>The kernel-radius maps and the "empty upper-right of det1" are the paddles, not the chambers; the chamber efficiency numbers stand.</li>
</ol>
<p class="muted">Code: <code>cosmic_vmm_analysis/scint_dump.py</code>, <code>paddle_efficiency.py</code>,
<code>paddle_study.py</code> (P2_basket_analysis, branch <code>worktree-trigger-acceptance-mc</code>); write-up
<code>FINDINGS_paddle_efficiency.md</code>; products in
<code>analysis_wk2/paddle_study/</code>, per-trigger dumps in <code>analysis_wk2/dumps/</code>.</p>

<h2>6. Caveats</h2>
<ul>
<li>The paddle z positions are not measured: the top pair at z&nbsp;=&nbsp;&minus;100&hellip;0 from the logbook, the bottom paddle at det1&nbsp;+&nbsp;31&nbsp;cm from memory. The measured length of 44 (__L44__&nbsp;mm against 630&ndash;640 for the top paddles at the same resolution) suggests the bottom paddle is nearer, z&nbsp;&asymp;&nbsp;1150, i.e. about 21&nbsp;cm below det1. A wrong z stretches its map along the band but does not change the efficiency step.</li>
<li>The chamber z spacing (310&nbsp;mm, SPS value) is also unmeasured; it scales all extrapolations alike.</li>
<li>No muon energy spectrum or multiple scattering: muons are straight lines. Showers and accidental coincidences are not simulated.</li>
<li>A "clean cluster" is a single 100&nbsp;% efficient hit in the MC. Cluster size, noise clusters and masked pads are not modelled; the chambers' own efficiency non-uniformity is what remains in the closure.</li>
<li>The paddle-position fit of section 3 is kept for the record only; its &minus;73&nbsp;mm bottom-paddle offset is an artefact of fitting an efficiency step with a position.</li>
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
const zMid = (Math.min(...G.paddles.map(p => p.z[0])) + Math.max(...G.paddles.map(p => p.z[1]))) / 2;
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
const pz = G.paddles;   // measured: {name, corners (extent), eff_corners (above half plateau), z: [z0, z1]}
for (const p of pz) {
  for (const f of boxFaces(p.corners, p.z[0], p.z[1])) prims.push({pts: f, fill: '--p-paddle', alpha: 0.12, stroke: true, layer: 'pad'});
  for (const f of boxFaces(p.eff_corners, p.z[0] + 2, p.z[1] - 2)) prims.push({pts: f, fill: '--p-paddle', alpha: 0.45, stroke: false, layer: 'pad'});
}
const topZ0 = Math.min(...pz.map(p => p.z[0])), botZ1 = Math.max(...pz.map(p => p.z[1]));
// tracks: split into segments for depth sorting
const zA = topZ0 - 60, zB = botZ1 + 60, NSEG = 14;
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
  for (const p of pz) lab(p.name, [p.corners[1][0], p.corners[1][1], p.z[0]]);
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
    '__GEO_HEAD__': tr(['element', 'centre (x, y) [mm]', 'long axis to x', 'size [mm]', 'z [mm] (down)', 'efficiency, lower-left / upper-right half', 'source'], th=True),
    '__ROWS_GEO__': rows_geo2,
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
    '__E44LL__': f"{pe['p0_44']['eff_lower_left_half']:.2f}", '__E44UR__': f"{pe['p0_44']['eff_upper_right_half']:.2f}",
    '__E41LL__': f"{pe['p0_41']['eff_lower_left_half']:.2f}", '__E41UR__': f"{pe['p0_41']['eff_upper_right_half']:.2f}",
    '__E45LL__': f"{pe['p0_45']['eff_lower_left_half']:.2f}", '__E45UR__': f"{pe['p0_45']['eff_upper_right_half']:.2f}",
    '__CLD3__': f"{cl[NM[0]]['data_mean_along']:+.0f}", '__CLD2__': f"{cl[NM[1]]['data_mean_along']:+.0f}", '__CLD1__': f"{cl[NM[2]]['data_mean_along']:+.0f}",
    '__CLM3__': f"{cl[NM[0]]['mc_measured_mean_along']:+.0f}", '__CLM2__': f"{cl[NM[1]]['mc_measured_mean_along']:+.0f}", '__CLM1__': f"{cl[NM[2]]['mc_measured_mean_along']:+.0f}",
    '__CLU3__': f"{cl[NM[0]]['mc_uniform_mean_along']:+.0f}", '__CLU2__': f"{cl[NM[1]]['mc_uniform_mean_along']:+.0f}", '__CLU1__': f"{cl[NM[2]]['mc_uniform_mean_along']:+.0f}",
    '__L44__': f"{pg['p0_44']['length_mm']:.0f}", '__NSTRAIGHT__': f"{PE['n_straight']}",
    '__D44__': f"{pe['p0_44']['delay_ns_per_m_along']:+.0f}",
    '__E40__': f"{PE['eff_on_straight_tracks']['p0_40']:.2f}", '__ETRIG__': f"{PE['eff_trigger_41or45_and_44']:.2f}", '__E45__': f"{PE['eff_on_straight_tracks']['p0_45']:.2f}",
    '__TC_HEAD__': tr(['top paddle(s) fired', 'share of triggers', 'mean slope along the band', 'tag median along, det3', 'det2', 'det1'], th=True),
    '__ROWS_TC__': rows_tc,
    '__PE_HEAD__': tr(['paddle', 'efficiency on straight tracks (whole wedge)', 'lower-left half', 'upper-right half', 'region above half plateau, along [mm]', 'pulse delay along the band [ns/m]'], th=True),
    '__ROWS_PE__': rows_pe,
    '__CL_HEAD__': tr(['probe', 'data: mean along [mm]', 'MC, measured paddle maps', '&chi;&sup2; / bins (along)', 'MC, uniform paddles', '&chi;&sup2; / bins'], th=True),
    '__ROWS_CL__': rows_cl,
}
for k, v in sub.items():
    HTML = HTML.replace(k, v)
assert '__' not in HTML.replace('__init__', ''), [l for l in HTML.split('\n') if '__' in l][:3]
open(note, 'w').write(HTML)
print('wrote', note, f'{len(HTML) / 1024:.0f} kB')
