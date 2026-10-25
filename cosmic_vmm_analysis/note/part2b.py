CSS = '''
<style>
  :root{--s1:#2a78d6;--s2:#eb6834;--s3:#1baf7a;--alert:#eb6834;--surface:#fcfcfb;
    --c0:#eef0f3;--c1:#cde2fb;--c2:#86b6ef;--c3:#3987e5;--c4:#1c5cab;--c5:#0d366b;
    --t0:#52514e;--t1:#1d1f24;--t2:#0b0b0b;--t3:#fff;--t4:#fff;--t5:#fff;}
  @media (prefers-color-scheme: dark){:root:not([data-theme="light"]){--s1:#3987e5;--s2:#d95926;--s3:#199e70;--alert:#e8743f;--surface:#1c1f25;
    --c0:#262a31;--c1:#184f95;--c2:#256abf;--c3:#3987e5;--c4:#86b6ef;--c5:#cde2fb;
    --t0:#a5a9b1;--t1:#e8e8e4;--t2:#e8e8e4;--t3:#0b0b0b;--t4:#0b0b0b;--t5:#0b0b0b;}}
  :root[data-theme="dark"]{--s1:#3987e5;--s2:#d95926;--s3:#199e70;--alert:#e8743f;--surface:#1c1f25;
    --c0:#262a31;--c1:#184f95;--c2:#256abf;--c3:#3987e5;--c4:#86b6ef;--c5:#cde2fb;
    --t0:#a5a9b1;--t1:#e8e8e4;--t2:#e8e8e4;--t3:#0b0b0b;--t4:#0b0b0b;--t5:#0b0b0b;}
  .viz-card{background:var(--bg-elev);border:1px solid var(--border);border-radius:var(--radius);padding:14px 14px 10px;margin:18px 0;box-shadow:var(--shadow)}
  .viz-card svg{width:100%;height:auto;display:block}
  svg text.t{fill:var(--fg)} svg text.m{fill:var(--fg-muted)} svg text.f{fill:var(--fg-faint)}
  svg .b0{fill:var(--c0)} svg .b1{fill:var(--c1)} svg .b2{fill:var(--c2)} svg .b3{fill:var(--c3)} svg .b4{fill:var(--c4)} svg .b5{fill:var(--c5)}
  svg text.x0{fill:var(--t0)} svg text.x1{fill:var(--t1)} svg text.x2{fill:var(--t2)} svg text.x3{fill:var(--t3)} svg text.x4{fill:var(--t4)} svg text.x5{fill:var(--t5)}
  .caption{font-size:13px;color:var(--fg-muted);margin:6px 2px 4px}
  .slider{display:block;margin:6px 4px 10px;font-size:14px} .slider input{width:60%;vertical-align:middle}
  .grid-two{display:grid;grid-template-columns:1fr 1fr;gap:12px;margin:6px 4px}
  .big{font-size:30px;font-weight:650} .small{font-size:12.5px}
  .tldr{border-left:4px solid var(--s1);background:var(--accent-soft);padding:12px 16px;border-radius:0 8px 8px 0;margin:14px 0}
  .tldr li{margin:3px 0}
  .stepn{flex:none;width:26px;height:26px;border-radius:50%;background:var(--s1);color:#fff;font-weight:650;font-size:14px;display:inline-flex;align-items:center;justify-content:center}
  .stepcard{border:1px solid var(--border);border-left:5px solid var(--s1);border-radius:8px;background:var(--bg-elev);padding:10px 14px 4px;margin:14px 0}
  .stepcard.next{border-left-color:var(--s2)} .stepcard.later{border-left-color:var(--fg-faint)}
  .sc-head{display:flex;gap:10px;align-items:center;margin-bottom:6px;flex-wrap:wrap}
  .stepcard dl{margin:0;display:grid;grid-template-columns:150px 1fr;gap:4px 12px}
  .stepcard dt{color:var(--fg-muted);font-size:13.5px} .stepcard dd{margin:0 0 6px;font-size:15px}
  @media (max-width:640px){.stepcard dl{grid-template-columns:1fr}.stepcard dt{margin-top:4px}}
  .two{display:grid;grid-template-columns:1fr 1fr;gap:14px} @media (max-width:640px){.two{grid-template-columns:1fr}}
  .two>div{border:1px solid var(--border);border-radius:8px;padding:8px 14px;background:var(--bg-elev)}
  .two h3{margin:4px 0 6px}
  details{margin:6px 0;border:1px solid var(--border);border-radius:8px;background:var(--bg-elev);padding:4px 12px}
  details summary{cursor:pointer;font-weight:600;padding:4px 0}
  table.bltab{border-collapse:collapse;margin:6px 0 10px;font-variant-numeric:tabular-nums;font-size:13px}
  table.bltab th{font-weight:500;color:var(--fg-muted);text-align:left;padding:2px 10px 2px 0;white-space:nowrap}
  table.bltab td{padding:2px 10px;text-align:right;border-top:1px solid var(--border)} table.bltab td.na{color:var(--fg-faint)}
  table.runs td,table.runs th{padding:4px 8px;font-size:13.5px}
  dl.gloss{display:grid;grid-template-columns:170px 1fr;gap:4px 14px;margin:8px 0 14px} dl.gloss dt{font-weight:600;font-size:14px} dl.gloss dd{margin:0;font-size:14.5px}
  @media (max-width:640px){dl.gloss{grid-template-columns:1fr} dl.gloss dt{margin-top:6px}}
</style>'''


def card(svg, caption):
    return f'<div class="viz-card">{svg}<p class="caption">{caption}</p></div>'


BODY = f'''
<div class="tldr"><strong>In one minute.</strong>
<ul>
<li><strong>The question:</strong> how low can we set the threshold of each VMM chip before electronic noise takes over? The beam test could only <em>guess</em> the answer. Here, with cosmic muons, we can <em>measure</em> it.</li>
<li><strong>What we did, in order:</strong> (1) a first test, which showed a flood of data from a few chips; (2) a noise scan with the chamber voltage on; (3) the same scan with the voltage <em>off</em>, to separate the detector from the electronics; (4) a closer look at single channels; (5) a direct measurement of the baseline of each channel, using a chip setting where a firing channel makes its neighbours report; (6) the same scans with the three chips we had masked (VMM 4, 7, 9) switched back on; (7) the first measurements with the high voltage on: a health check, and the neighbour read-out with cosmic muons, which gave a baseline for 97&nbsp;% of the channels.</li>
<li><strong>What we found:</strong> the chips do <em>not</em> start to be noisy at the same threshold, and this spread stays when the voltage is off, so it comes from the electronics, not the chambers. Three chips (VMM 4, 7, 9) have noisy channels that also stay with the voltage off; when we switch them back on they show two kinds of noise: an ordinary smooth onset like the other chips, plus bursts that jump up and down with each new configuration, even at very high thresholds. And each channel has a noise peak of its own, so the channels inside one chip differ too. The baselines we could measure directly (about half of the channels so far) are mostly 45 to 65 ADC counts, not the 38 used so far for every channel, and they differ from chip to chip. With the voltage on we now have a baseline for 97&nbsp;% of the channels, the voltage does not change them, and VMM 4, 7 and 9 each have a large block of channels with baselines of 100 to several hundred, next to the channels that burst.</li>
<li><strong>What this means for the next step:</strong> one threshold per chip is not enough. We need the baseline of <em>every channel</em> (half is done), in threshold units (millivolts, not yet done), and then to correct each channel. Then we repeat the scan and, only then, measure efficiency. The measurement with muon tracks (reference run and threshold scan, high voltage on) is running now, with a global threshold for each chip read from the noise scans; the reference run (8 hours) and the threshold scan are done and show a steady muon rate of about 21 per minute for the stack, but also bursts that carry a quarter of the hits and have to be cut. An interleaved 24-hour run at four thresholds is taking the statistics for the efficiency.</li>
<li><strong>Update, Sunday 25 October evening (step 10c):</strong> the tag-and-probe analysis, with the beam-test pad map, shows that the &ldquo;muon-like&rdquo; events are <em>not</em> muon tracks: their hits do not line up in a straight line through the three chambers. Counting the events that do form a straight line, minus what chance gives, leaves about 16 per hour at the safe setting, where muons should give hundreds to thousands. Below the safe setting the extra &ldquo;tracks&rdquo; are a fixed pattern of a few hot channels. The pad map itself is fine (no other wiring does better, no chamber is flipped). So the stack sees almost no cosmic muons above the noise at these settings, and the best threshold cannot be chosen from efficiency yet. The long run was stopped at 21:34 (runs 96&ndash;111). Next: check the high voltage and the gas against the beam test, and a source test.</li>
</ul></div>

<h2>0. The idea in one picture</h2>
<p>A VMM channel records a hit when the pulse it sees goes above a <strong>threshold</strong>. If the threshold is high, small muon pulses are missed (low efficiency). If it is low, electronic noise also goes above it (many false hits). Move the slider to feel the trade-off:</p>
{widget_threshold()}
<p>In the SPS beam data the VMM efficiency was about 85&nbsp;% against 95.6&nbsp;% for the DREAM readout of the same chambers, and the analysis said the VMM threshold was probably too high. That was a <em>prediction</em>: the beam gave no time to vary the threshold. In the lab we can, because cosmic muons arrive at a steady, gentle rate and we can wait.</p>

<h2>1. The set-up</h2>
{card(cartoon_stack(), 'Three chambers stacked like at the beam test. Nine hybrids with two VMM chips each read them out (VMM 2 to 19). The chips decide for themselves when to record a hit; there is no trigger detector.')}
<p>Two things to keep in mind: the voltages (410&nbsp;V mesh, 710&nbsp;V drift) are set by hand and are not recorded by the DAQ, so we note them in a logbook; and the gas (Ar/isobutane 95/5) is <strong>not yet confirmed</strong>.</p>

<h2>2. A word about the threshold knob</h2>
<p>The threshold of a chip is a number (<code>sdt</code>, in &ldquo;DAC units&rdquo;, just a counting scale). Each chip rests at a slightly different level, its <strong>baseline</strong>: a threshold below the baseline makes the chip fire all the time. So one common number for all chips is a trap.</p>
{card(fig_baselines(), 'Estimated baselines from the July beam data. With one common threshold of 200, the chips above the band would fire all the time.')}
<p>That is why the scan is written as <strong>d = threshold &minus; baseline</strong> for each chip: d = 50 is a comfortable margin, d = 16 is close to the noise. These baselines are only <em>estimates</em> from July. Remember this: it turns out to matter.</p>

<h2>3. Step by step</h2>
<p>Each step below says what we did, what we saw, what it told us, and why the next step follows. Nothing is skipped.</p>

{stepcard(1, 'A first test with the July settings',
  'We ran the July chip settings for 2 minutes with the chambers at voltage.',
  '1.5&nbsp;GB of data in 2 minutes: about 1.6 million hits per second, almost all from 3 chips (VMM 4, 9 and later 7).',
  'Something fires non-stop on those chips. A scan with this flood would fill the disk and hide the other chips.',
  'find <em>which channels</em> fire, and try to switch them off (mask them).')}

{stepcard(2, 'Which channels? Try masking',
  'We counted hits per channel, switched off the hottest ones, and tried again.',
  'A single channel gave nearly everything (VMM 4 ch 49: 99&nbsp;% of that chip). Once masked, the burst jumped to other channels (it &ldquo;walks&rdquo;). Masking whole chips VMM 4 and 9 brought the data from 635 to 5&nbsp;MB per minute.',
  'Masking channel by channel does not work for these chips, so we leave VMM 4 and 9 out of the scan (they give no curve). The others are quiet enough.',
  'scan the threshold on all the other chips.')}
{card(fig_hop(), 'Real counts. Switching off the hot channels did not help: the burst moved to channels 20 and 21 on VMM 9, and 36 and 38 on VMM 4.')}

{stepcard(3, 'Noise scan with the chamber voltage ON',
  'Ten threshold steps, d = 50 down to 16, 2 minutes each. The chip settings are the same except the threshold.',
  'Almost nothing at d = 50. As d falls, chips switch on one after the other, and then climb very fast. VMM 15 starts at d = 34, VMM 10 only at d = 18.',
  'The chips do not start to be noisy at the same distance above the baselines we assumed (a spread of 16 steps), so our baselines are not right. But is the spread a property of the electronics or of the chambers?',
  'repeat the same scan with the voltage <em>off</em>.')}
{card(fig_heatmap(rows=ROWS_ON, runs=ON_RUN, title='Voltage ON: hit rate per chip (kHz) at each threshold step'), 'Real hit rates. Each row is a chip, each column a threshold step (quiet on the left). The orange outline is the first step above 1&nbsp;kHz. Hover a cell for the number.')}

<h2>4. Why switch the voltage off?</h2>
<p>With the voltage on, a hit can come from the electronics <em>or</em> from the detector (a muon, a discharge, small pad noise). With the voltage off, the gas gives no signal, so every hit must come from the electronics. Doing the same scan twice separates the two.</p>
{card(cartoon_chain(), 'The idea of the two scans. Same chips, same settings, same cables; only the chamber voltage changes.')}

{stepcard(4, 'The same scan with the voltage OFF',
  'You switched all six high-voltage channels off by hand. We ran the same ten steps, and then four more steps with nothing masked at all, to see the burst chips.',
  'At d = 50 every chip is at zero. Chips start to fire in nearly the same order as before, usually 0 to 6 steps later (VMM 5 is the biggest exception). VMM 4, 7 and 9 still burst. Some channels are still noisy.',
  'The spread between chips and the bursts are properties of the electronics (chip, hybrid or cable). The small background at d = 50 and most noisy channels need the detector.',
  'understand a strange feature we saw: some chips&rsquo; rate <em>falls</em> when we lower the threshold further.')}
{card(fig_heatmap(rows=ROWS_OFF, runs=OFF_RUN, title='Voltage OFF: hit rate per chip (kHz) at each threshold step', foot6='* VMM 6: silent at every step with the voltage off (no hits at all), so it shows 0.00.'), 'Real hit rates with the voltage off, same layout as before. VMM 4 and 9 are masked as in the first scan.')}
{card(fig_onset_pair(), 'The first step above 1&nbsp;kHz for each chip, voltage on (blue dot) and off (orange ring). Where they sit together, the ring surrounds the dot. The order of the chips is almost the same, and the spread of about 16 steps is still there with the voltage off.')}
{card(fig_volume_pair(), 'Data per minute, voltage on and off. At moderate thresholds the voltage adds a lot of hits; at the lowest thresholds the electronics dominate and the two get closer.')}
{card(fig_bursts_off(), 'The three burst chips with nothing masked and the voltage off. They burst all the same, so the cause is not the pads: it is the chip, the hybrid or the cable.')}
{card(fig_class(), 'Channel by channel: for every channel that fired with the voltage on, did it still fire with the voltage off? Most disappear (blue). VMM 12 and VMM 15 keep many channels (orange): the whole chip is noisy by itself.')}

<h2>5. The strange curves, explained</h2>
<p>In the voltage-on scan, some chips got <em>quieter</em> at the lowest thresholds (for example VMM 5: 416&nbsp;kHz at d = 20, then 0.3&nbsp;kHz at d = 18). Lowering a threshold should never give fewer hits&hellip; unless we look at <em>one channel</em> at a time.</p>
{card(cartoon_scurve(), 'A channel records a hit each time its (noisy) baseline crosses the threshold going up. Far above, the noise rarely reaches it. Far below, the channel is always above the threshold, so there are no new crossings. The most crossings happen when the threshold sits at the channel&rsquo;s own baseline.')}
{card(fig_channel_curves(), 'Real data: single channels follow that bell shape. VMM 7 channel 23 peaks at d = 34 and VMM 5 channel 46 at d = 20, with or without voltage. The place of the peak tells where the threshold meets that channel&rsquo;s level.')}

{stepcard(5, 'Look at single channels',
  'We followed individual channels across the threshold steps.',
  'Each channel has a bell-shaped curve with a peak. Different channels of the same chip peak at different steps (for VMM 7: channel 23 peaks at step 34, while channel 10 is still rising at step 16).',
  'Inside one chip the channels peak at different steps, so the level that the threshold has to clear differs from channel to channel. One high channel can make the whole chip look noisy, and a single threshold per chip cannot suit all channels. Our picture is that the peak sits where the threshold meets the channel&rsquo;s own level; the next section tests this.',
  'measure the baseline of every channel directly (next section).')}

<h2>6. Measuring the baselines directly</h2>
<p>Up to here the baselines of the channels were only <em>guesses</em>: estimates from the July beam data, and the positions of the noise peaks. Can we measure them? Yes, with a setting of the chip called <code>sng = 1</code> (&ldquo;neighbour read-out&rdquo;). When one channel fires, its two neighbours are read out too, even though they are below threshold. A neighbour that is read out without any pulse reports the level it rests at: its baseline.</p>
{card(cartoon_neighbour(), 'A firing channel (orange) makes its two neighbours (blue) report their ADC value. The hit word carries a flag that tells the two cases apart, so the data contain both the hits and the baselines.')}
<p>The catch: with noise as the only trigger, only the neighbours of channels that fire are covered. And a channel fires most at its own bell-curve peak (section 5), so each threshold step makes <em>different</em> channels fire. A short <strong>ladder of thresholds</strong> therefore adds up to many channels:</p>
{card(cartoon_ladder(), 'Each step lets other channels fire, and their neighbours get measured. The cartoon uses 16 channels and three steps.')}

{stepcard(6, 'Neighbour read-out with the voltage OFF',
  'We switched on <code>sng = 1</code> on every chip and took 2 minutes at five settings: one first test, and then four thresholds per chip, 0, 3, 6 and 9 steps lower than that chip&rsquo;s noise onset from step 4. VMM 4 and 9 stay masked. The voltage stayed off, so the baselines are measured without any detector effect.',
  f'Between {NB_LO:.0f} and {NB_HI:.0f}&nbsp;% of the hits on every chip are neighbour read-outs, as expected. After adding the five captures, <strong>{BL_COV} of the {BL_TOTAL} channels</strong> on the 16 working chips had a baseline (at least 50 read-outs); the count after each capture was {", ".join(str(c) for c in COV_STEPS)}, and the two lowest thresholds added the most. The whole of each 2-minute capture was used (599&nbsp;s in total). The median baseline is {BL_MEDIAN:g}&nbsp;ADC; chip medians range from {CHIP_MED_MIN:g} to {CHIP_MED_MAX:g}; only {N_LE38} of the {BL_COV} channels are at or below 38.',
  'The single value used so far for every channel (38) is wrong for almost all of them, and the baselines differ between chips and, by a few ADC counts, inside each chip. About {100 * (BL_TOTAL - BL_COV) / BL_TOTAL:.0f}&nbsp;% of the channels still have no measurement (VMM 6 and VMM 7 the most). All the data taken have now been used, so more coverage needs new data: more steps, longer runs, or cosmic muons with the voltage on (they hit all channels evenly).',
  'put the baselines in threshold units (millivolts), which needs one more measurement (step 12). First, though, the three chips we had masked must be put back into the study (next section).')}
{card(fig_baseline_heat(), 'Real data: every cell is one channel. Hatched cells have no measurement yet. Hover a cell for the value and the number of read-outs. VMM 6 and VMM 7 are mostly hatched: they have the fewest measurements.')}
{card(fig_baseline_range(), 'The same, chip by chip: the bar runs from the lowest to the highest channel, the dot is the chip median. Every chip sits well to the right of the old value of 38. Three channels are far above the plot and not understood yet: VMM 16 channel 33 (ADC 704), VMM 6 channel 14 (186) and VMM 14 channel 34 (125).')}
{card(fig_baseline_panels(), 'The same data, but one panel per chip, so every channel can be read. Each blue dot is one channel; the ticks at the bottom mark channels without a measurement yet. Orange dots lie outside the plot range (below 30 or above 72) and are drawn at its edge. Hover a dot for the channel, its value and the number of read-outs.')}
<h3>The number for every channel</h3>
<p>The exact values are in the tables below, one per chip, eight channels per row (open a chip to see it; hover a number for its read-out count). The same data are in a file you can download: <a href="2026-10-22-measuring-the-vmm-noise-threshold-with-cosmic-muons-what-we-did-and-what-we-saw_files/neighbor_baseline_all_channels.csv">neighbor_baseline_all_channels.csv</a> (columns: <code>vmm</code>, <code>ch</code>, <code>n_over</code> = hits where the channel itself fired, <code>n_nb</code> = neighbour read-outs, <code>adc_nb_median / mean / std</code> = the baseline and its spread, <code>adc_over_median</code> = median ADC of its own hits, <code>first_file</code> = the first capture in which it reached 50 read-outs, &minus;1 if never). A channel counts as measured when <code>n_nb</code> is at least 50; rows with fewer are kept but should not be used.</p>
{baseline_tables()}
{card(fig_coverage_table(), 'How the coverage grew as captures were added (cumulative): the first test, then the four steps. The last steps, at the lowest thresholds, still added many channels, so the ladder works.')}

<h3>Does the baseline predict where a channel&rsquo;s noise peaks?</h3>
<p>If a channel fires most when the threshold sits at its baseline, then a channel with a higher measured baseline should peak at a higher threshold step. We can now test this. The voltage-off scan does not reach its peaks (the rates still rise at the last step), so we used the peak positions from the voltage-on scan and the baselines measured with the voltage off (the baseline is an electronics property).</p>
{card(fig_baseline_vs_peak(), f'Each dot is a channel with both a measured baseline and a clear noise peak ({PEAK_N} channels here, peaks with at least 10 hits). A line through the dots is almost flat. The rank correlation is {", ".join(f"{PEAK_ROWS_ALL[m][1]:+.2f} ({PEAK_ROWS_ALL[m][0]} channels, at least {m} hits)" for m in (30, 10, 5) if m in PEAK_ROWS_ALL)}; it gets smaller as more channels are included, and the ordinary correlation is small (+0.04 to +0.21). Channels above 100 ADC (the three outliers) are left out of the plot and of these numbers.')}
<p><strong>So the answer is no, at least not clearly.</strong> An earlier look with fewer channels showed a weak positive trend; with the whole data it largely disappears. Three reasons are possible, and none is tested yet: (1) the ADC baseline is not the same thing as the level the threshold acts on, because each channel&rsquo;s threshold comparator has its own small offset; (2) the threshold steps are coarse (3 to 10 units) and many peaks have few hits; (3) the peaks come from the voltage-on scan and the baselines from the voltage-off scan. Reason (1) is exactly why the chip monitor measurement of step 12 matters: it measures the pedestal <em>and</em> the threshold of each channel, so it gives the missing link.</p>
{card(cartoon_languages(), 'The chain we are building. We set the threshold in DAC steps and see the baseline in ADC counts. To correct channel by channel we also need the level in millivolts, which the chip can report through its monitor.')}


<h2>7. Putting VMM 4, 7 and 9 back in</h2>
<p>For the scans above we had masked VMM 4 and VMM 9 completely (and one channel of VMM 7), because their bursts filled the disk. That left them out exactly where the comparison with the other chips matters. So we repeated the measurement with the voltage off and these three chips unmasked and following the threshold, from very high thresholds down to low ones, once with the normal read-out and once with neighbour read-out. All other chips sat at a quiet setting (d = 50). The three chips were read at 22 thresholds in total; each point is 2 minutes.</p>
{card(fig_burst_ladder(), 'Real data, voltage off, nothing masked. For each chip the rate does <em>not</em> go up smoothly as the threshold goes down: it jumps by orders of magnitude between neighbouring settings. At the low-threshold end, though, a steady rise appears (VMM 4 and VMM 9).')}
{burst_table()}
{card(cartoon_two_kinds(), 'The two behaviours that the figure above contains. The steady rise at low threshold is the same as on the other chips and can be placed on the same scale. The bursts are something else.')}

{stepcard(7, 'The three burst chips, unmasked',
  'We repeated the voltage-off scan with VMM 4, 7 and 9 switched back on, at 13 thresholds with the normal read-out (d = 200 down to 16) and 9 with neighbour read-out (d = 100 down to 14), the other chips at d = 50. We also read the temperature of all chips (31 to 38 &deg;C, nothing unusual).',
  'Two things. (a) A smooth noise onset like the other chips: VMM 4 passes 1&nbsp;kHz at d = 22 (hottest channels 37 and 35) and VMM 9 at d = 27 (channel 16); VMM 7 shows none down to d = 16. (b) Bursts of up to millions of hits per second that do not follow the threshold: VMM 9 gives 2.7&nbsp;MHz at d = 200, far above any baseline, and VMM 4 gives 0 at d = 60 but 1.8&nbsp;MHz at d = 40. The hot channel changes with each configuration (VMM 4: channel 45, 39, 38/36, then 37/35). Within one run it does not change, and its pulse height is almost constant.',
  'The chips can be placed on the same threshold scale as the rest (the smooth part), and the bursts are not ordinary noise: they look like a state a channel gets stuck in when the chip is configured. We have not yet tested whether applying the same configuration twice gives the same state.',
  'repeat a few configurations several times (a short run, about 15 minutes, with the voltage off) to see whether the burst state is repeatable.')}

<h2>8. First measurements with the voltage on</h2>
<p>The high voltage was switched on after the runs above (we declare it in the run files; the voltages themselves are not read by the DAQ and are still to be written in the logbook). Three things were done with it on: a health check of the chambers, the neighbour read-out measurement again but now with cosmic muons, and the noise ladder of VMM 4, 7 and 9 again.</p>

<h3>8a. Is every chamber alive?</h3>
<p>We took 10 minutes at the &ldquo;safe&rdquo; setting of each chip (8 steps above the noise onset we measured with the voltage on) and looked at the pulse heights of the hits.</p>
{card(fig_health(), 'Real data, 2 minutes of the 10-minute run. Every chip has hits and a broad spectrum. VMM 2, 3, 4, 6 and 7 of chamber det3 reach the ADC limit; the other two chambers stay between about 500 and 850.')}
<p>No chamber looks like the July failure of det3 (a collapsed gain). The det3 chips reaching the ADC limit is unexplained: it could be a higher gain than the other chambers or discharges, and the spectrum shape with more statistics will tell. And 6 to 100 hits per second and chip is well above the cosmic rate we expect (about 4 per second), so at this setting most hits are noise or background, and muon tracks will be picked out by requiring a coincidence of the three chambers.</p>

<h3>8b. A baseline for (nearly) every channel</h3>
<p>With the voltage off only channels that fire by themselves were covered, about 40&nbsp;% of the channels were still missing. With the voltage on, cosmic muons and chamber activity make every pad fire sometimes, so every channel becomes the neighbour of a firing one. We ran the neighbour read-out at five thresholds (10 minutes each, 50 minutes in total).</p>
{card(fig_cov_compare(), f'Channels with a measured baseline: {ON_COV} of {ON_TOTAL} ({100 * ON_COV / ON_TOTAL:.0f}&nbsp;%) with the voltage on, against {OFF_COV} with the voltage off. {FULL64} of the 18 chips are complete. The missing ones are mostly the masked channels of the three burst chips.')}
<p>Does the voltage change the baseline? We compared, for the 15 normal chips, the {CORE_N} channels measured both ways:</p>
{card(fig_on_off_scatter(), f'The median difference is {CORE_MED:+.1f} ADC count (typical spread {CORE_MAD:.0f}), {CORE_W3:.0f}&nbsp;% of the channels agree within 3 counts and {CORE_W5:.0f}&nbsp;% within 5. So the voltage does not move the baselines of these chips, and the voltage-off measurement was fine for them.')}
<p>Over all {ON_COV} channels the median baseline is {ON_MEDIAN:.0f} ADC, and only {ON_LE38} channels ({100 * ON_LE38 / ON_COV:.0f}&nbsp;%) are at or below the 38 we had assumed for all of them.</p>

<h3>8c. VMM 4, 7 and 9 have a block of channels with a high baseline</h3>
<p>With complete coverage the three burst chips look different from every other chip. Among the 15 normal chips just one channel out of about {len(NORMAL_CHIPS) * 64} has a baseline above 80 ADC (VMM 16 channel 33, which reads 704). On VMM 4, 7 and 9 a large block of channels does:</p>
{card(fig_block_panels(), f'Baseline of every channel, voltage on. VMM 12 (right) is a normal chip. On VMM 4, VMM 7 and VMM 9 a block of channels in the upper part has baselines of 100 to several hundred ADC: {BLOCK[4]["n_high"]} of {BLOCK[4]["n_cov"]} measured, unmasked channels on VMM 4 (from channel {BLOCK[4]["start"]}), {BLOCK[7]["n_high"]} of {BLOCK[7]["n_cov"]} on VMM 7 (from {BLOCK[7]["start"]}, up to {BLOCK[7]["vmax"]:.0f}), {BLOCK[9]["n_high"]} of {BLOCK[9]["n_cov"]} on VMM 9 (from {BLOCK[9]["start"]}).')}
{card(cartoon_block(), 'The picture in one drawing. The bursts we have followed since the first day sit at the edge of the block: the channels that were hot in at least three runs are 34 to 39 on VMM 4 (the block starts at 40), 22 to 25 on VMM 7 (block from 26) and 15 to 27 on VMM 9, plus 33 and 44 (block from 23).')}
<p>Where a channel of the block was also measured with the voltage off, it was high then too: 6 of the 7 channels we could compare (VMM 4 channels 40 and 41, VMM 7 channel 26, VMM 9 channels 23, 28 and 31; the seventh read 0, probably because it was masked then). So the block does not come from the chamber. This is a small sample, and we don&rsquo;t know why the block exists: it looks like a property of the chip, the hybrid or the connection of those channels. It is a good thing to hand to the people who built these hybrids.</p>

<h3>8d. The noise ladder of the three chips, with the voltage on</h3>
{card(fig_ladder_onoff(), f'The same configurations as before, voltage on (dot) and off (ring), taken about an hour apart on the same evening. At most thresholds the two agree, including the hot channel. VMM 4 matches in {LAD_CMP[4]["same"]} of 13 thresholds (and both are quiet in {LAD_CMP[4]["quiet"]}), VMM 7 in {LAD_CMP[7]["same"]} (quiet in {LAD_CMP[7]["quiet"]}), VMM 9 in {LAD_CMP[9]["same"]} (quiet in {LAD_CMP[9]["quiet"]}).')}
<p>Two things stand out. First, a burst at a high threshold is the same in both runs: VMM 4 at d = 40 gives 1740&nbsp;kHz with the voltage on and 1770 with it off, both on channel 36; VMM 7 at d = 34 gives 546 and 547&nbsp;kHz, both on channel 23. So the burst state follows the threshold setting and is not a matter of chance at each configuration, at least at these settings; the planned Monday repeat test will check this directly. Second, near the noise onset (d from 34 down to 27) the chips are noisier with the voltage on (VMM 4 at d = 30: 85&nbsp;kHz against almost nothing with it off). That is the chamber adding hits near the onset, as on every other chip.</p>

{stepcard(8, 'First measurements with the voltage on',
  'The voltage was switched on by hand, and we took a 10-minute health check at the safe setting, the noise ladder of VMM 4, 7 and 9 again (13 points), and the neighbour read-out with cosmic muons at five thresholds. The runs are <code>run_57</code> to <code>run_75</code>.',
  f'No chamber is dead. Baselines now exist for {ON_COV} of {ON_TOTAL} channels, and the voltage does not change them on the 15 normal chips. VMM 4, 7 and 9 each have a block of channels with baselines of 100 to several hundred ADC, next to the channels that burst. The bursts at high thresholds are the same with the voltage on and off.',
  'The baseline map is essentially complete. The three burst chips are different in a clear and measurable way, and their problem looks like a property of the hardware, not of the chambers or of chance.',
  'set a global threshold for each chip from the noise scans (step 9), then take the efficiency measurement with muon tracks (step 11), keeping the burst chips in it, with the 23 hot channels masked.')}

<h2>9. A global threshold for each chip</h2>
<p>Each chip has one threshold number, <code>sdt</code>, and the noise scans told us at which value each chip&rsquo;s noise starts. That number does not depend on our baseline estimates (they only fix the d scale), so we can read the global threshold of each chip straight from the scan. For the muon runs we use <strong>8 steps above the noise onset measured with the voltage on</strong>: low enough to catch small pulses, with a margin against the noise.</p>
{card(fig_thresholds(), 'Real data. For most chips the voltage-on and voltage-off onsets are close (0 to 6 steps apart); VMM 5 is the exception, 11 steps apart. For VMM 2 and VMM 7 no onset was found in the range scanned, so their setting is a guess; for VMM 4 and 9 only the voltage-off onset can be read, because the bursts hide the voltage-on one; VMM 6 was held at 230 in every scan.')}
{thr_table()}

{stepcard(9, 'A global threshold for each chip',
  'We read, for every chip, the threshold at which its noise passes 1&nbsp;kHz, with the voltage on and off (gain 4.5&nbsp;mV/fC, peaking time 200&nbsp;ns), and chose onset + 8 steps for the muon runs.',
  'Thirteen chips have a clean onset with the voltage on (VMM 3, 5, 8, 10 to 19); the two sets of onsets agree within a few steps. Five chips are special (VMM 2, 4, 6, 7, 9, see the figure).',
  'The global threshold of 13 of the 18 chips is determined from data. It is only a floor: a noise onset says where the noise starts, not where the efficiency is best, and one noisy channel sets it for the whole chip. Only the muon data can say whether the margin of 8 steps is right, and a per-channel correction would let a chip go lower.',
  'compare efficiency and noise at several settings with muon tracks (step 11), and keep the five special chips under watch.')}

<h2>10. First look at the muon events</h2>
<p>The reference run is already taking data with the thresholds above. Can we find the muons in it? We took the first 30 minutes (the voltage on, the safe setting, the 23 hot channels of VMM 4, 7 and 9 masked) and grouped the hits that happen at the same moment.</p>
{card(cartoon_three_kinds(), 'The three kinds of event we expect to see in the three chambers.')}
{card(fig_event_sizes(), f'Real data. {EVM["hits"]:,} hits in 30 minutes, {EVM["events"]:,} groups. More than half of the groups are a single, isolated hit; {EVM["large_events"]:,} groups have more than 36 hits and carry {EVM["hits_in_large_events_pct"]:.0f}&nbsp;% of all hits.')}
<p><strong>Yes, the muons are there:</strong> {EVM["three_chambers"]:,} groups have hits in all three chambers, and {EVM["muon_like"]} of them are small (at most 12 hits per chamber, typically exactly one per chamber), which is what a muon looks like. That is about {EVM["muon_like"] / EVM["span_s"] * 60:.0f} per minute for the whole stack. The time difference between the chambers in these events is centred on zero (the first hit of each chamber, median 0 to &minus;7 ns), so there is no big constant timing offset between the hybrids, and the spread is about 40&nbsp;ns, because the fine timing is not calibrated yet.</p>
<div class="callout warn"><strong>Correction (25 October, step 10c).</strong> &ldquo;Yes, the muons are there&rdquo; was too quick. The later analysis put each hit at its pad position and asked whether the three hits form a straight line: almost none of these events do, beyond what chance gives. Most of the small three-chamber events are noise that fires the chambers at the same moment. The text above is kept as it was written.</div>
<p><strong>The complications we have to deal with</strong> (this is what makes the correlation between the chambers harder than it looks):</p>
<ul>
<li><strong>Bursts.</strong> A quarter of all hits sit in a few thousand groups with dozens to hundreds of hits across many chips and chambers at the same moment. If we just ask for hits in all three chambers within 100&nbsp;ns, these bursts give millions of false pairs, against a few thousand real ones. They must be cut by size before anything else. We do not know what they are (discharges in the chambers, or pick-up in the electronics); the repeat test with the voltage off would tell us.</li>
<li><strong>Noise.</strong> Half of all hits are isolated singles. They cost nothing in a three-chamber coincidence, but they fake a probe hit: the efficiency must be corrected for the chance that a noise hit lands next to the muon&rsquo;s track.</li>
<li><strong>Few tracks.</strong> About 19 per minute for the whole stack, fewer than the 40 per chip we hoped for (not yet corrected for the masked channels, the pad size or the angles). An hour gives about 1100 tracks for the stack, so an efficiency per chip is limited to a few per cent per hour. We need hours of data per setting.</li>
<li><strong>Timing.</strong> The timing is not calibrated (about 40&nbsp;ns spread), so the coincidence window must be about &plusmn;100&nbsp;ns. This is fine at these rates.</li>
<li><strong>Geometry.</strong> We still have to map the pads of each chip to positions, check the orientation of each chamber, and use the beam-test distances. Cosmic muons arrive at wide angles, so the line between two chambers must allow a wide search box in the third.</li>
<li><strong>The three burst chips.</strong> Their masked channels leave holes in the efficiency map, and their bursts must be rejected by the size cut as well.</li>
<li><strong>det3.</strong> Its pulse heights reach the ADC limit (step 8). That is a question of gain or discharges, and it matters for how we define a hit.</li>
</ul>

<h3>10b. The weekend runs: what the muon data look like over a day</h3>
<p>Since the first look the stack took 8 hours at the &ldquo;safe&rdquo; setting (16 runs of 30 minutes, all complete) and then four one-hour points at other thresholds. The same grouping of hits was applied to every run.</p>
{card(fig_weekend(), f'Real data. Left: over the 8 hours the muon-like rate stays at {REF_MU:.1f} &plusmn; {REF_SD:.1f} per minute for the stack, with no drift; the total hit rate, however, rises by about half (914 to 1430 per second) and the large bursts go from 0.7 to 1.2 per second, for a reason we do not know. Right: the number of muon-like events per minute at five settings.')}
<p>Reading the right-hand panel: from 16 steps above the safe setting down to 4 steps below it, the yield rises from 13 to 22 per minute (steps of +22&nbsp;%, +27&nbsp;%, +6&nbsp;%), and the share of the hits that sit in large bursts falls from 27&nbsp;% to 13&nbsp;% because more ordinary hits join. At 8 steps below the safe setting the number jumps to 469 per minute: the chips are now so noisy that chance coincidences between the three chambers look like muons. So the useful range ends somewhere between 4 and 8 steps below the safe setting.</p>
<p><strong>This does not tell us the best threshold.</strong> The yield counts real muons and chance coincidences together. The flattening between 0 and &minus;4 hints that the efficiency is levelling off, but only the tag-and-probe analysis, with the noise correction, can show it. What the runs do show is that the muon rate is steady, the time offsets between the chambers are stable (&minus;11 to +1 ns, spread 35 to 48 ns in every run), and that the interesting settings are the ones from 8 steps above the safe setting down to about 6 steps below it.</p>

{stepcard(10, 'A first look at the muon events',
  'We took the first 30 minutes of the reference run and grouped the hits that occur at the same time (within 100&nbsp;ns).',
  f'{EVM["muon_like"]} muon-like events in 30 minutes (about 19 per minute), a time spread of about 40&nbsp;ns between chambers with no big offset, and a quarter of all hits in large bursts.',
  'The muon data can be analysed, and no hardware change is needed. But the analysis needs a size cut for the bursts, a coincidence window of about &plusmn;100&nbsp;ns, a noise correction and a pad-to-position map before it gives an efficiency.',
  'write the tag-and-probe analysis with these cuts, and run it on the reference run and on the threshold scan (step 11).')}

<div class="callout warn"><strong>Correction (25 October).</strong> The &ldquo;muon-like&rdquo; yield above is mostly noise that hits the three chambers at the same moment (step 10c): only about 16 per hour of these events are straight tracks beyond chance. The steady rate over 8 hours and the stable timing are real, but they describe this shared noise, not muons.</div>

<h3>10c. Are these really muons? The tag-and-probe analysis</h3>
<p><strong>The idea.</strong> The two outer chambers &ldquo;tag&rdquo; a muon; the middle chamber is the &ldquo;probe&rdquo;: did it see the muon too? The fraction of yes answers is the efficiency. This is the same method as at the beam test, with the outer chambers playing the role of the reference detectors. It only works if the tags really are muons.</p>
<p><strong>First try, without positions.</strong> We just asked: when det3 and det1 both fire, does det2 fire at the same moment?</p>
{tp_table()}
<p>This looks like an efficiency of about 55&nbsp;% at the safe setting that falls at lower thresholds. It is not one. At lower thresholds the outer chambers fire together on noise much more often, so the &ldquo;tags&rdquo; are mostly fake and the fraction drops; and a det2 hit anywhere in the chamber counts, not only where the muon went.</p>
<p><strong>Second try, with positions.</strong> We put every hit on its pad with the beam-test map and asked whether the det2 hit is where the line through the outer hits says it should be.</p>
{card(cartoon_midpoint(), 'Lying flat, 31 cm apart: a straight line crosses the middle chamber half way between its crossings of the outer two, whatever its angle. This test needs no other detector.')}
<p>Only 3.6&nbsp;% of the tags had a det2 hit at the right place, although 63&nbsp;% had a det2 hit somewhere at the right time. Either the map is wrong for this set-up (new Samtec cables, chambers turned flat), or the events are not muons. We checked both.</p>
<p><strong>Is the map wrong?</strong> We followed the method used on the cosmic bench: try every reasonable way the 64 channels of each chip could be wired to the pads (40 per connector: in order, reversed, pairs swapped, the two halves swapped, odd and even pins split), and keep the one that makes the most events straight. To be fair, we fitted on half of the events and tested on the other half.</p>
{card(fig_wiring(), 'Real data, runs 76&ndash;91. The beam-test map does best on events it has never seen. A different wiring can look better on the events it was tuned on, but it does worse on new events: it learned the noise. We also scanned rotations and mirror images between the chambers: no chamber is flipped.')}
<p>One detail is useful for the future: at the beam test two connector halves of P2_IN needed a correction because of the Hitachi cables. With the Samtec cables they follow the plain order.</p>
<p><strong>So are the events muons?</strong> We counted, for every run, the groups that contain a straight triplet (one hit per chamber, middle hit within 30&nbsp;mm of the midpoint), and subtracted what chance gives. Chance is measured by moving each middle-chamber hit to another channel <em>of the same chip</em>: this keeps the noise that is shared between chips, and only real straight lines survive.</p>
{card(fig_straight_rate(), 'Real data, runs 76&ndash;108, all settings of the threshold scan and of the long run. At the safe setting about 16 per hour, from 10.6 hours. For comparison, a rough estimate: the stack covers about 0.05 m&sup2; per chamber, so muons crossing all three should be of order one per second, that is a few thousand per hour times the efficiency of each chamber to the third power (about 400 per hour at 50&nbsp;% per chamber).')}
{str_table()}
<p>The rise at 6 and 8 steps below the safe setting looks like muons appearing, but it is not:</p>
{card(fig_pattern(), 'Real data. The straight events at low threshold come again and again from the same few channels in the three chambers firing together: hot channels with correlated noise, not particles.')}
<p>We also looked at the big bursts (more than 36 hits): about 190 hits in under a microsecond, spread evenly over every chip of a chamber, mostly det3 and det2. They are whole-chamber events, not tracks. And the pulse heights of the hits in three-chamber events are no larger than those of single noise hits: no sign of a muon signal there either.</p>
<p><strong>What this means.</strong> With these settings the stack sees almost no cosmic muons above the noise: the chambers fire together, but on noise. The best threshold cannot be chosen from the efficiency until a muon signal is there. At the beam test the same chambers and chips gave 67&ndash;85&nbsp;%, so something is different here. The first suspects are the high voltage and the gas (the voltage is set by hand and not logged; the gas is not confirmed), then the new cables. A radioactive source (Fe-55) on one chamber would show directly whether signals come through, and would check the map at the same time.</p>
{stepcard('10c', 'Tag-and-probe: are these muons?',
  'Tag-and-probe without and with pad positions; the beam-test map against 40 alternative wirings per connector and rotations / mirrors between chambers; straight tracks above chance for every run, runs 76&ndash;108.',
  'The map is the best one available and no chamber is flipped. Straight tracks above chance: about 16 per hour at the safe setting, nothing more at the other settings except a fixed pattern of hot channels at the lowest ones.',
  'The three-chamber coincidences are shared noise, not muons. The efficiency, and so the best threshold, cannot be measured until the chambers show a muon signal.',
  'check the high voltage and the gas against the beam test, and do a source test (step 12b).')}

<h2>11. What we know, and what we do not</h2>
<div class="two">
<div><h3>We know</h3><ul>
<li>With no voltage, the electronics are silent at d = 50.</li>
<li>The chip-to-chip spread of the noise onset (about 16 steps) is an electronics effect.</li>
<li>VMM 4, 7 and 9 have noisy channels that need no voltage. Unmasked, they also show an ordinary smooth noise onset (VMM 4 at d = 22, VMM 9 at d = 27) plus bursts whose size and channel change with each configuration.</li>
<li>VMM 12 and VMM 15 are noisy over many channels without voltage.</li>
<li>A single channel&rsquo;s noise rate has a bell-shaped curve with a peak.</li>
<li>With neighbour read-out we can measure a channel&rsquo;s baseline directly. The value used so far (38) is wrong for nearly all channels; the chips differ (medians {CHIP_MED_MIN:g} to {CHIP_MED_MAX:g}).</li>
<li>The measured ADC baselines do <em>not</em> predict where a channel&rsquo;s noise peaks (correlation between 0 and 0.4).</li>
<li>With the voltage on, 97&nbsp;% of the channels have a baseline, and the voltage does not change them on the 15 normal chips.</li>
<li>VMM 4, 7 and 9 each have a block of channels with a high baseline (also with the voltage off); the bursts are at its edge and, at high thresholds, are the same with the voltage on and off.</li>
<li>No chamber is dead with the voltage on.</li>
<li>A global threshold for 13 of the 18 chips can be read from the noise scans (onset + 8 steps); the five others are special.</li>
<li>Muon-like events are in the data (about 21 per minute for the stack, steady over 8 hours), with no large timing offset between the chambers, and the usable range of thresholds ends between 4 and 8 steps below the safe setting, where noise coincidences take over. <em>Corrected on 25 October: these are mostly shared noise, not muons (next item).</em></li>
<li>Only about 16 per hour of the three-chamber events are straight tracks beyond chance, at every setting from 16 steps above to 4 below the safe one; the extra &ldquo;tracks&rdquo; at lower thresholds are a fixed pattern of hot channels.</li>
<li>The beam-test channel-to-pad map is the best one for this set-up; no chamber is flipped or mirrored; the two P2_IN corrections of the Hitachi cables are not needed with the Samtec cables.</li>
</ul></div>
<div><h3>We do not know yet</h3><ul>
<li>The <strong>efficiency</strong>: this was the noise side and a first look at the muon events.</li>
<li>Which threshold gives the best efficiency (the yield rises down to 4 steps below the safe setting, but it mixes muons and chance coincidences), and what the large bursts (a quarter of all hits) are.</li>
<li><strong>Why the chambers show almost no muon signal here</strong>, when they gave 67&ndash;85&nbsp;% at the beam test: high voltage, gas, cables?</li>
<li>Why the total hit rate rose by half over the 8 hours of the reference run while the muon rate did not.</li>
<li>Why the three burst chips have that block of channels, and whether it is the chip, the hybrid or the connection.</li>
<li>The link between ADC counts and millivolts, channel by channel.</li>
<li>Whether the ADC baseline and the level the threshold acts on are the same thing, and why three channels read very high values (VMM 16 ch 33: 704, VMM 6 ch 14: 186, VMM 14 ch 34: 125).</li>
<li>Whether applying the very same configuration twice on the same day gives the same burst (the Monday repeat test), and what the det3 pulse heights that reach the ADC limit mean.</li>
<li>Whether the lab noise looks like the beam-hall noise (about 100 times higher in the beam hall).</li>
<li>What the voltage actually is: it is set by hand and not logged.</li>
</ul></div>
</div>

<h2>12. What comes next</h2>
{card(cartoon_global_vs_channel(), 'Today each chip has one threshold. The next step gives each channel a small correction so that all channels start at the same level (cartoon, not measured).')}

{stepcard(11, 'Reference run, threshold scan and the interleaved long run',
  'Reference run (8 hours at the safe setting), threshold scan (four points of one hour), and the interleaved long run from Sunday 13:22 to 21:34: 16 runs of 30 minutes (run_96 to run_111) in rotation at 8 steps above the safe setting, the safe setting, 4 and 6 steps below it. Stopped on request after the analysis of step 10c.',
  'About 33 hours of data at six settings in total. The tag-and-probe analysis (step 10c) found no muon signal at any of them, so more hours at the same settings would not help.',
  'The data are good for noise and stability studies (the total rate drift, the bursts, the hot channels), but not yet for efficiency.',
  'find out why there is no muon signal (step 12b) before taking more muon data.')}

{stepcard('12b', 'Why no muon signal? High voltage, gas, source test',
  'The high voltage was switched on by hand on 23 October and is still on; its values are not in the logbook. The gas mixture is not confirmed in the run configuration.',
  'The voltage of each chamber and the gas, compared with the beam test; then an Fe-55 source on one chamber for a few minutes, with the same chip settings.',
  'Whether signals come through at all: a clear excess of hits under the source, at the right pads, would show the chambers work and confirm the map; no excess would point to the gain (voltage or gas).',
  'once there is a muon signal, repeat the threshold scan and the tag-and-probe analysis (step 14).', status='next')}

{stepcard(12, 'Put the baselines in millivolts',
  'Not started. The new slow-control library can already read the chip&rsquo;s monitor for each channel, but the command-line tool does not offer it.',
  'A way to call the library: either ask the developer to add a command, or write a small program of our own that calls it (without changing his software). We would test it on one channel of one hybrid first, then run all channels of all hybrids.',
  'The level of every channel in millivolts. Together with the ADC baselines of step 6 this gives the full chain, threshold setting to millivolts to ADC, channel by channel.',
  'compute the corrections (step 13).', status='next')}

{stepcard(13, 'Compute and apply the corrections, then scan again',
  'Not started.',
  'The measurements of steps 6 to 12.',
  'If the corrections work, the noise onset of all chips lines up and the usable threshold window opens for every chip.',
  'repeat the noise scan with the corrections, voltage on and off, then measure the efficiency (step 14).', status='later')}

{stepcard(14, 'Measure the efficiency',
  'Not started.',
  'A night of cosmic muons at a safe threshold, then the tag-and-probe method: the two outer chambers say a muon went through, and the middle one is asked whether it saw it. The analysis is written and was tried on 33 hours of data (step 10c); it needs a muon signal first (step 12b).',
  'The efficiency against the threshold, chip by chip: the real answer to the first question.',
  'compare with the beam test and with the DREAM readout.', status='later')}

<h2>13. All runs at a glance, and the words we use</h2>
<p>This is the map of the whole campaign so far: which run numbers belong to which measurement, and what the names mean. The same overview now opens the logbook.</p>
{card(fig_run_timeline(), 'From the logbook. 111 run numbers in four days; numbers 33 and 34 were lost in the run-counter incident of Friday afternoon.')}

<h3>The words</h3>
<dl class="gloss">
<dt>run_N</dt><dd>One start-to-stop of the DAQ, with one chip configuration. The data are in <code>runs/run_N/cfg_&lt;config&gt;/raw_daq_data/</code>.</dd>
<dt>config label</dt><dd>The chip configuration files are called <code>p2b-config-cern-ext_gain4.5_peaktime200_opt_cosmic_&lt;label&gt;.txt</code>. Only the label changes between runs, so the label is what we quote: e.g. <code>m2hvoff_d34</code>. Every run used gain 4.5&nbsp;mV/fC and 200&nbsp;ns peaking time.</dd>
<dt>sdt</dt><dd>The threshold register of a chip, in DAC counts. One value per chip.</dd>
<dt>d</dt><dd>The threshold step of the noise scans: d = sdt &minus; the chip&rsquo;s assumed baseline (the July estimate per chip). <code>_d34</code> means every chip at its baseline + 34. Lower d = lower threshold = more noise. VMM 6 has no baseline estimate and stays at sdt 230.</dd>
<dt>onset</dt><dd>For each chip, the d at which its hit rate first passes 1&nbsp;kHz in a noise scan.</dd>
<dt>safe</dt><dd>The working threshold of each chip: its onset with the voltage on, plus 8 DAC (VMM 6 at sdt 230). Config <code>safehvon_v1</code>.</dd>
<dt>safe + k, safe &minus; k</dt><dd>Every chip k DAC above or below its safe value. Configs <code>thrhvon_p16</code>, <code>_p8</code> (+16, +8) and <code>_m4</code>, <code>_m6</code>, <code>_m8</code> (&minus;4, &minus;6, &minus;8).</dd>
<dt>masks</dt><dd>A masked channel is switched off (<code>sm 1</code>). <code>m1</code>: single hot channels (did not work, the burst hops); <code>m2</code>: VMM 4 and 9 fully masked plus VMM 7 ch 24; <code>m0</code>: nothing masked; from run 57 on (all <code>safehvon</code> and <code>thrhvon</code> configs): 23 channels of VMM 4, 7 and 9.</dd>
<dt>hvoff, hvon</dt><dd>The chamber voltage when the config was made. The voltage is set by hand and not read by the DAQ.</dd>
<dt>n479, s479</dt><dd>VMM 4, 7 and 9 unmasked and following the threshold d, all other chips parked at d = 50. <code>n</code> = normal read-out, <code>s</code> = neighbour read-out.</dd>
<dt>neighbour read-out (sng 1)</dt><dd>A chip setting that makes the neighbours of a firing channel report their pulse height too; a neighbour&rsquo;s value is its baseline. Labels <code>sng1hvoff_*</code> (0/&minus;3/&minus;6/&minus;9 DAC from each chip&rsquo;s onset, voltage off), <code>s479*</code> and <code>sngcosmic_oK</code> (safe &minus; K, voltage on).</dd>
<dt>burst</dt><dd>VMM 4, 7 and 9 fire a band of channels at up to MHz rates; the band moves with each new configuration. Not the same as the &ldquo;large events&rdquo; of the muon runs (more than 36 hits within 100&nbsp;ns, about 1 per second).</dd>
<dt>reference run, threshold scan, long run</dt><dd>Long data at the safe setting; one-hour points at other settings; 30-minute runs cycling through four settings so that a slow drift affects them all alike.</dd>
<dt>det3, det2, det1</dt><dd>The three chambers: det3 = P2_IN (VMM 2&ndash;7), det2 = P2_MID (VMM 8&ndash;13), det1 = P2_OUT (VMM 14&ndash;19), 31&nbsp;cm apart, lying flat.</dd>
</dl>

<h3>The groups</h3>
{run_groups_table()}
<p><strong>The long run by setting:</strong> safe + 8: runs 96, 100, 104, 108; safe: 97, 101, 105, 109; safe &minus; 4: 98, 102, 106, 110; safe &minus; 6: 99, 103, 107, 111.<br>
<strong>All muon data (voltage on) by setting:</strong> safe: 76&ndash;91, 97, 101, 105, 109 (20 runs, about 11 hours); safe + 16: 92; safe + 8: 93, 96, 100, 104, 108; safe &minus; 4: 94, 98, 102, 106, 110; safe &minus; 6: 99, 103, 107, 111; safe &minus; 8: 95.</p>
<div class="callout"><strong>Reading the times in the logbook.</strong> Summer time ended during run 92 (03:00 became 02:00 on 25 October): run 92 is logged at 02:55 and ran its full hour, run 93 is logged at 02:56 but started one hour later. A note typed in the logbook is filed under the run that was current at that moment, so the note announcing a new block sits under the last run of the block before.</div>

<h2>14. Where the data and code are</h2>
<table class="runs">
<tr><th>Item</th><th>Where</th></tr>
<tr><td>Voltage-on scan: <code>run_4</code> &hellip; <code>run_13</code>; first tests: <code>run_1</code> &hellip; <code>run_3</code></td><td><code>/local/p2/p2data/cosmics_3det_202610/runs/</code> on vmm_daplxa, copied to EOS <code>&hellip;/p2_sps_july/vmm_cosmic_run_oct26/</code></td></tr>
<tr><td>Voltage-off scan: <code>run_14</code> &hellip; <code>run_23</code> (same configs), <code>run_24</code> &hellip; <code>run_27</code> (nothing masked)</td><td>same place; configs end in <code>m2hvoff_d&lt;d&gt;</code> and <code>m0hvoff_d&lt;d&gt;</code></td></tr>
<tr><td>Neighbour read-out (sng 1), voltage off: <code>run_28</code> (first test), <code>run_29</code> &hellip; <code>run_32</code> (thresholds 0, &minus;3, &minus;6, &minus;9 from each chip&rsquo;s onset)</td><td>same place; configs <code>&hellip;_cosmic_sng1hvoff_probe / _o0 / _m3 / _m6 / _m9</code>. These runs are large: 2 to 9&nbsp;GB each, about 25&nbsp;GB in total. Analysed on lxplus from the EOS backup, so that the DAQ machine was not loaded</td></tr>
<tr><td>VMM 4, 7, 9 unmasked, voltage off: <code>run_35</code> &hellip; <code>run_47</code> (normal read-out), <code>run_48</code> &hellip; <code>run_56</code> (neighbour read-out)</td><td>same place; configs <code>&hellip;_cosmic_n479hvoff_d&lt;d&gt;</code> and <code>&hellip;_s479hvoff_d&lt;d&gt;</code>. Two stray captures from a software problem at the start are kept apart in <code>orphan_2026-10-23_*</code> (not part of any run)</td></tr>
<tr><td>Voltage ON, 24 Oct: <code>run_57</code> (health check), <code>run_58</code> &hellip; <code>run_70</code> (noise ladder of VMM 4, 7, 9), <code>run_71</code> &hellip; <code>run_75</code> (neighbour read-out with cosmics, 4 to 44&nbsp;GB each)</td><td>same place; configs <code>&hellip;_safehvon_v1</code>, <code>&hellip;_n479hvon_d&lt;d&gt;</code>, <code>&hellip;_sngcosmic_o&lt;k&gt;</code>. A sixth neighbour point was refused by the disk guard. Nothing was recorded from Friday 20:24 to Saturday 18:33: the reference run had not been chained behind these scans, and the connection to the machine was lost (then the machine&rsquo;s account service stopped accepting the login)</td></tr>
<tr><td>First look at the muon events: the first 30 minutes of the reference run (<code>run_76</code>)</td><td><code>analysis_muons/</code> in the data folder; per-chip thresholds in <code>analysis_thresholds/</code> (chip table and a per-channel baseline table, CSV)</td></tr>
<tr><td>Reference run and threshold scan (done)</td><td><code>run_76</code> &hellip; <code>run_91</code>: 16 &times; 30&nbsp;min at <code>&hellip;_safehvon_v1</code>; <code>run_92</code> &hellip; <code>run_95</code>: <code>&hellip;_thrhvon_p16 / _p8 / _m4 / _m8</code> (60&nbsp;min each); per-run results in <code>analysis_muons/per_run/</code></td></tr>
<tr><td>Interleaved long run (Sunday 13:22 to 21:34, stopped on request)</td><td><code>run_96</code> &hellip; <code>run_111</code>: 30&nbsp;min each in rotation <code>&hellip;_thrhvon_p8</code>, <code>&hellip;_safehvon_v1</code>, <code>&hellip;_thrhvon_m4</code>, <code>&hellip;_thrhvon_m6</code> (<code>run_111</code> stopped after 21 minutes)</td></tr>
<tr><td>Tag-and-probe, map and straight-track analysis (step 10c)</td><td>Code in <code>P2_basket_analysis/cosmic_vmm_analysis/</code> (linked in the DAQ repo as <code>cosmic_vmm_analysis</code>): <code>tag_and_probe.py</code>, <code>tag_and_probe_pads.py</code>, <code>muon_signal_check.py</code>, <code>dump_candidates.py</code>, <code>orientation_check.py</code>, <code>straight_track_rate.py</code>. Results in <code>analysis_muons/</code> (<code>straight_track_rate_runs76-108.json</code>, <code>orientation_runs76-91.json</code>, <code>signal_check/</code>, candidate events in <code>cand/</code>)</td></tr>
<tr><td>Tables used for the figures</td><td><code>analysis_noise_scan/</code> and <code>analysis_noise_scan_hvoff/</code> (chip tables, per-channel CSV, voltage on/off comparison), <code>analysis_n479_hvoff/</code> and <code>analysis_hvon/</code> (rate and hot channels of VMM 4, 7, 9 with the voltage off and on; baselines with the voltage on), <code>analysis_health/</code> (pulse-height spectra of the health check), <code>analysis_sng1_hvoff/</code> (baseline of every covered channel, coverage per capture, baseline against noise peak)</td></tr>
<tr><td>Chip settings</td><td>4.5&nbsp;mV/fC, 200&nbsp;ns peaking time; only the threshold lines and the masks change between files</td></tr>
<tr><td>Code and feedback notes</td><td>branch <code>lab-cosmics-oct26</code>, <code>docs/LAB_COSMICS_OBSERVATIONS.md</code> (items 1&ndash;64); the analysis scripts moved on 25 October from the DAQ repo&rsquo;s <code>quick_scripts/</code> to <code>P2_basket_analysis/cosmic_vmm_analysis/</code></td></tr>
</table>
'''

html = open(NOTE).read()
a = html.index('<style>')
b = html.index('</style>') + len('</style>')
html = html[:a] + CSS + html[b:]
a = html.index('<div class="note-body">')
b = html.index('</div>\n  </article>') + len('</div>')
html = html[:a] + '<div class="note-body">' + BODY + '</div>' + html[b:]
html = _re.sub(r'(<meta name="updated" content=")[^"]*(">)', r'\g<1>2026-10-25\2', html)
html = _re.sub(r'(<meta name="description" content=")[^"]*(">)',
               r'\g<1>A picture book of the lab threshold scans on the three-chamber cosmic stand: what we did step by step, with the chamber voltage on and off, and what we found per chip.\2', html)
open(NOTE, 'w').write(html)
print('written', NOTE, len(html))
