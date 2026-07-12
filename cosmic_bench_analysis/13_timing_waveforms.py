#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
13_timing_waveforms.py

Timing resolution of the P2 pad detector from raw decoded waveforms, comparing
time-of-arrival (TOA) algorithms against the scintillator trigger.

Reference frame: the DAQ window opens on the 2-PMT scintillator coincidence
(~5 ns resolution per PMT), so the extracted waveform time *within the window*
(plus the DREAM fine-timestamp correction, see p2_waveforms) is the detector
time relative to the trigger. Two complementary benchmarks:

  vs-trigger : sigma of the per-hit TOA distribution across events. Contains
               detector resolution (+) trigger jitter (+) residual clock terms.
  pad-pad    : sigma of (t_pad1 - t_pad2)/sqrt(2) for 2-pad events of the same
               muon. The trigger cancels event-by-event -> intrinsic per-pad
               resolution (assuming equal, independent pads).

For every algorithm both benchmarks are computed before and after an
amplitude-slewing (time-walk) correction. Multi-pad events additionally test
cluster-TOA combination strategies (max-amp pad / earliest / amplitude-weighted
/ inverse-variance-weighted) against the trigger.

The ftst sign convention is validated empirically: mean TOA vs ftst for the
three variants (+, -, none) — the correct one is flat (stage prints the slopes
and uses the flattest).

Products (<Analysis>/<det>/<run>/<sub_run>/13_timing/):
  waveform_gallery<sfx>.png       example waveforms + algorithm markers
  ftst_validation<sfx>.png        mean TOA vs ftst per sign convention
  time_walk_<algo><sfx>.png       TOA vs amplitude + slewing correction
  toa_distributions<sfx>.png      per-algorithm TOA distributions (walk-corr.)
  resolution_comparison<sfx>.png  sigma summary: vs-trigger and pad-pad
  combo_comparison<sfx>.png       multi-pad cluster-TOA combination benchmark
  timing_benchmarks<sfx>.csv      all numbers
  timing_summary<sfx>.txt

Usage:
  python3 13_timing_waveforms.py [run_key] [--max-files N] [--max-fits N]
          [--amp-min 60] [--no-veto-sparks]
"""

import os
import json
import argparse
import matplotlib
matplotlib.use('Agg')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import p2_qa_config as qa
import p2_mapping as pmap
import p2_sparks as ps
import p2_waveforms as pw


def sig_label(s, n):
    return f'{s:.1f} ns (n={n:,})' if np.isfinite(s) else 'n/a'


def main():
    ap = argparse.ArgumentParser(description='P2 waveform timing benchmarks.')
    ap.add_argument('run_key', nargs='?', default=qa.DEFAULT_RUN)
    ap.add_argument('--max-files', type=int, default=None,
                    help='cap on decoded files (chunk x FEU) to process.')
    ap.add_argument('--max-fits', type=int, default=40000,
                    help='cap on waveforms for the fit-based algorithms '
                         '(sigmoid); cheap algorithms always see all hits.')
    ap.add_argument('--amp-min', type=float, default=60.0,
                    help='min baseline-subtracted amplitude [ADC] for the '
                         'benchmark sample (kills noise-edge hits).')
    ap.add_argument('--veto-sparks', action=argparse.BooleanOptionalAction,
                    default=True)
    ap.add_argument('--sub-run', default=None,
                    help='override the config sub_run (e.g. one HV-scan point '
                         'scan_mesh_410V_drift_580V); outputs go to that '
                         'sub_run\'s Analysis dir.')
    args = ap.parse_args()

    cfg = qa.get_config(args.run_key)
    if args.sub_run:
        cfg.SUB_RUN = args.sub_run
        cfg.OUT_BASE = os.path.join(cfg.DATA_ROOT, 'Analysis', cfg.DET_TAG,
                                    cfg.RUN, args.sub_run)
    print(cfg)
    out_dir = cfg.out_dir('13_timing')
    sfx = cfg.product_suffix(args.veto_sparks)

    # sampling period from the DAQ config (60 ns on the cosmic bench)
    with open(cfg.run_config_path) as f:
        rc = json.load(f)
    tps = float(rc.get('dream_daq_info', {}).get('sample_period', 60.0))
    print(f'  time per sample: {tps:g} ns, ftst unit: {pw.TIME_PER_FTST:g} ns')

    # channel table -> this detector's (feu, channel) -> pad
    ct = pmap.build_channel_table(cfg.run_config_path, cfg.MAP_CSV_PATH,
                                  det_type=cfg.DET_TYPE, det_name=cfg.DET_NAME,
                                  strategy='reverse',
                                  drop_connectors=cfg.DEAD_CONNECTORS)
    feus = ct.attrs['feus']
    files = pw.list_decoded_files(cfg.decoded_root_dir, feus=feus)
    if args.max_files:
        files = files.head(args.max_files)
    if not len(files):
        raise SystemExit(f'no decoded files for FEUs {feus} in '
                         f'{cfg.decoded_root_dir} — fetch decoded_root from rays.')
    print(f'  decoded files: {len(files)} '
          f'(chunks {sorted(files["chunk"].unique())}, FEUs {sorted(files["feu"].unique())})')

    # ---------------- hit extraction ----------------
    parts, wparts = [], []
    for feu, gf in files.groupby('feu'):
        first = gf.iloc[0]['path']
        ped = pw.pedestal_from_data(first)
        rms = pw.noise_rms_from_data(first, ped)
        print(f'  FEU {feu}: pedestal median {np.median(ped):.0f} ADC, '
              f'noise median {np.median(rms):.2f} ADC')
        for _, row in gf.iterrows():
            df, wv = pw.extract_hits(row['path'], ped, rms)
            df['feu'] = feu
            df['chunk'] = row['chunk']
            parts.append(df)
            wparts.append(wv)
            print(f'    {os.path.basename(row["path"])}: {len(df):,} hits')
    hits = pd.concat(parts, ignore_index=True)
    waves = np.vstack(wparts)

    # pad mapping: keep only channels wired to this detector's live connectors
    hits = hits.reset_index(drop=True)
    hits['_row'] = np.arange(len(hits))
    hits = hits.merge(ct[['feu', 'channel', 'channel_id', 'connector_N',
                          'pad_cx', 'pad_cy', 'mapped']],
                      on=['feu', 'channel'], how='inner')
    hits = hits[hits['mapped']].reset_index(drop=True)
    waves = waves[hits['_row'].to_numpy()]
    print(f'  detector hits after pad mapping: {len(hits):,} '
          f'({hits["eventId"].nunique():,} events)')

    # HV-spark / burst veto (same convention as the other stages)
    if args.veto_sparks:
        sv = ps.SparkVeto.from_cfg(cfg)
        bad = sv.vetoed_ids_from_hits(cfg.combined_hits_dir, feus)
        keep = ~hits['eventId'].isin(bad)
        hits, waves = hits[keep].reset_index(drop=True), waves[keep.to_numpy()]
        print(f'  spark veto: {len(sv.sparks)} sparks + {sv.last_burst_events} '
              f'bursts -> {len(hits):,} hits kept')

    # benchmark sample: real pulses, not saturated, decent amplitude
    good = (~hits['saturated']) & (hits['amp_max'] >= args.amp_min)
    npad_ev = hits.groupby('eventId')['channel_id'].transform('nunique')
    good &= npad_ev <= 6            # burst remnants out of the benchmark
    hb, wb = hits[good].reset_index(drop=True), waves[good.to_numpy()]
    print(f'  benchmark sample: {len(hb):,} hits '
          f'(amp >= {args.amp_min:g}, not saturated, <=6 pads/event)')

    # ---------------- TOA extraction (all algorithms) ----------------
    algos = dict(pw.ALGORITHMS)
    t_samp = {}
    for name, fn in algos.items():
        if name == 'sigmoid' and args.max_fits and len(wb) > args.max_fits:
            sub = np.random.RandomState(42).choice(
                len(wb), args.max_fits, replace=False)
            t = np.full(len(wb), np.nan)
            t[sub] = pw.t_sigmoid(wb[sub])
            t_samp[name] = t
        else:
            t_samp[name] = fn(wb, hb)
        nok = np.isfinite(t_samp[name]).sum()
        print(f'  {name:9s}: {nok:,} TOAs ({100*nok/max(len(hb),1):.0f}%)')

    # template algorithm: build the average pulse shape, then correlate
    t_up, template = pw.build_template(wb, hb['amp_max'].to_numpy(),
                                       sat=hb['saturated'].to_numpy())
    t_samp['template'] = pw.t_template(wb, t_up, template)
    nok = np.isfinite(t_samp['template']).sum()
    print(f'  {"template":9s}: {nok:,} TOAs ({100*nok/max(len(hb),1):.0f}%)')

    # ---------------- ftst sign validation (on frac30) ----------------
    ft = hb['ftst'].to_numpy(float)
    base = t_samp['frac30'] * tps
    variants = {'+ftst (reference)': base + ft * pw.TIME_PER_FTST,
                '-ftst': base - ft * pw.TIME_PER_FTST,
                'no correction': base}
    slopes = {}
    fig, ax = plt.subplots(figsize=(8, 5))
    for lab, t in variants.items():
        mu = [np.nanmedian(t[ft == v]) for v in range(8)]
        mu = np.array(mu) - np.nanmean(mu)
        ok = np.isfinite(mu)
        slope = np.polyfit(np.arange(8)[ok], mu[ok], 1)[0] if ok.sum() > 2 else np.nan
        slopes[lab] = slope
        ax.plot(range(8), mu, 'o-', label=f'{lab}  (slope {slope:+.2f} ns/ftst)')
    ax.set_xlabel('ftst [10 ns units]'); ax.set_ylabel('median TOA - mean [ns]')
    ax.set_title(f'{cfg.DET_NAME} ftst sign validation (frac30)\n'
                 f'{cfg.RUN}/{cfg.SUB_RUN}')
    ax.grid(True, alpha=0.3); ax.legend()
    fig.tight_layout()
    fig.savefig(f'{out_dir}/ftst_validation{sfx}.png', dpi=150,
                bbox_inches='tight')
    plt.close(fig)
    best = min(slopes, key=lambda k: abs(slopes[k]))
    ft_sign = {'+ftst (reference)': 1.0, '-ftst': -1.0, 'no correction': 0.0}[best]
    print(f'  ftst validation: flattest = "{best}" '
          f'(slopes: ' + ', '.join(f'{k} {v:+.2f}' for k, v in slopes.items()) + ')')

    # TOA in ns with the validated ftst convention
    t_ns = {a: t * tps + ft_sign * ft * pw.TIME_PER_FTST
            for a, t in t_samp.items()}

    # persist the per-hit TOA table for downstream kinematics analysis
    # (15_timing_kinematics: sigma vs amplitude / cosmic angle / chunk)
    toa_df = hb[['eventId', 'chunk', 'feu', 'channel_id', 'connector_N',
                 'pad_cx', 'pad_cy', 'amp_max', 'ftst', 'saturated']].copy()
    for a, t in t_ns.items():
        toa_df[f't_{a}'] = t
    toa_df.to_csv(f'{out_dir}/toa_hits{sfx}.csv', index=False)

    # ---------------- benchmarks per algorithm ----------------
    amp = hb['amp_max'].to_numpy()
    ev = hb['eventId'].to_numpy()
    rows = []
    t_corr = {}
    for a, t in t_ns.items():
        s_raw, q_raw = pw.robust_sigma(t), pw.q68_half_width(t)
        tc, walk_amp, walk_t = pw.slewing_correction(t, amp)
        t_corr[a] = tc
        s_cor, q_cor = pw.robust_sigma(tc), pw.q68_half_width(tc)

        # pad-pad: two highest-amp pads in >=2-pad events
        dfa = pd.DataFrame({'eventId': ev, 't': tc, 'amp': amp})
        dfa = dfa[np.isfinite(dfa['t'])]
        dfa = dfa.sort_values(['eventId', 'amp'], ascending=[True, False])
        g = dfa.groupby('eventId')
        two = dfa[dfa['eventId'].isin(g.size()[g.size() >= 2].index)]
        top2 = two.groupby('eventId').head(2)
        dt = top2.groupby('eventId')['t'].agg(lambda x: x.iloc[0] - x.iloc[1])
        s_pp = pw.robust_sigma(dt) / np.sqrt(2)
        rows.append(dict(algo=a, n=int(np.isfinite(t).sum()),
                         sigma_raw_ns=s_raw, q68_raw_ns=q_raw,
                         sigma_walkcorr_ns=s_cor, q68_walkcorr_ns=q_cor,
                         n_pairs=int(len(dt)), sigma_padpad_ns=s_pp))

        # time-walk plot
        fig, axs = plt.subplots(1, 2, figsize=(12, 4.5))
        ok = np.isfinite(t)
        axs[0].hexbin(amp[ok], t[ok] - np.nanmedian(t), gridsize=60,
                      xscale='log', cmap='viridis', mincnt=1)
        if len(walk_amp):
            axs[0].plot(walk_amp, walk_t - np.nanmedian(t), 'r.-', lw=1.5,
                        label='slewing profile')
            axs[0].legend()
        axs[0].set_xlabel('amplitude [ADC]')
        axs[0].set_ylabel('TOA - median [ns]')
        axs[0].set_title(f'{a}: time walk')
        okc = np.isfinite(tc)
        axs[1].hexbin(amp[okc], tc[okc], gridsize=60, xscale='log',
                      cmap='viridis', mincnt=1)
        axs[1].set_xlabel('amplitude [ADC]')
        axs[1].set_ylabel('TOA (walk-corrected) [ns]')
        axs[1].set_title(f'{a}: after slewing correction '
                         f'(sigma {s_cor:.1f} ns)')
        fig.suptitle(f'{cfg.DET_NAME} {a} — {cfg.RUN}/{cfg.SUB_RUN}')
        fig.tight_layout()
        fig.savefig(f'{out_dir}/time_walk_{a}{sfx}.png', dpi=130,
                    bbox_inches='tight')
        plt.close(fig)

    bench = pd.DataFrame(rows).sort_values('sigma_walkcorr_ns')
    bench.to_csv(f'{out_dir}/timing_benchmarks{sfx}.csv', index=False)

    # ---------------- multi-pad combination strategies ----------------
    best_algo = bench.iloc[0]['algo']
    tb = t_corr[best_algo]
    dfa = pd.DataFrame({'eventId': ev, 't': tb, 'amp': amp})
    dfa = dfa[np.isfinite(dfa['t'])]
    multi = dfa[dfa['eventId'].isin(
        dfa.groupby('eventId').size().loc[lambda s: s >= 2].index)]

    # per-amplitude-bin sigma for inverse-variance weights
    la = np.log10(multi['amp'])
    edges = np.nanpercentile(la, np.linspace(0, 100, 9))
    edges = np.unique(edges)
    sig_bin = {}
    ib_all = np.clip(np.digitize(la, edges) - 1, 0, len(edges) - 2)
    for b in range(len(edges) - 1):
        s = pw.robust_sigma(multi['t'].to_numpy()[ib_all == b])
        sig_bin[b] = s if np.isfinite(s) and s > 0 else 1e3

    def combos(g):
        t, a = g['t'].to_numpy(), g['amp'].to_numpy()
        ib = np.clip(np.digitize(np.log10(a), edges) - 1, 0, len(edges) - 2)
        w = 1.0 / np.array([sig_bin[b] ** 2 for b in ib])
        return pd.Series({
            'maxamp': t[a.argmax()],
            'earliest': t.min(),
            'amp_weighted': float((t * a).sum() / a.sum()),
            'invvar_weighted': float((t * w).sum() / w.sum()),
            'n_pads': len(t)})

    cmb = multi.groupby('eventId').apply(combos)
    combo_rows = [dict(combo=c, sigma_ns=pw.robust_sigma(cmb[c]),
                       q68_ns=pw.q68_half_width(cmb[c]), n=len(cmb))
                  for c in ['maxamp', 'earliest', 'amp_weighted',
                            'invvar_weighted']]
    # single-pad reference on the same footing
    single = dfa[~dfa['eventId'].isin(multi['eventId'])]
    combo_rows.append(dict(combo='single-pad events',
                           sigma_ns=pw.robust_sigma(single['t']),
                           q68_ns=pw.q68_half_width(single['t']),
                           n=len(single)))
    combo = pd.DataFrame(combo_rows)

    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.barh(combo['combo'], combo['sigma_ns'], color='mediumpurple')
    for i, r in combo.iterrows():
        ax.text(r['sigma_ns'] + 0.1, i, f'{r["sigma_ns"]:.1f} ns', va='center')
    ax.set_xlabel('sigma vs trigger [ns]')
    ax.set_title(f'{cfg.DET_NAME} multi-pad cluster TOA combinations '
                 f'({best_algo}, walk-corrected)\n{cfg.RUN}/{cfg.SUB_RUN}  '
                 f'({len(cmb):,} multi-pad events)')
    ax.grid(True, alpha=0.3, axis='x')
    fig.tight_layout()
    fig.savefig(f'{out_dir}/combo_comparison{sfx}.png', dpi=150,
                bbox_inches='tight')
    plt.close(fig)

    # ---------------- summary figures ----------------
    # TOA distributions (walk-corrected)
    fig, ax = plt.subplots(figsize=(9.5, 5.5))
    for a in bench['algo']:
        t = t_corr[a]
        t = t[np.isfinite(t)]
        med = np.median(t)
        ax.hist(t - med, bins=np.arange(-60, 60.1, 2), histtype='step', lw=1.5,
                label=f'{a} (sigma {bench.set_index("algo").loc[a, "sigma_walkcorr_ns"]:.1f} ns)')
    ax.set_yscale('log')
    ax.set_xlabel('TOA - median [ns]'); ax.set_ylabel('hits')
    ax.set_title(f'{cfg.DET_NAME} walk-corrected TOA vs trigger — '
                 f'{cfg.RUN}/{cfg.SUB_RUN}')
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(f'{out_dir}/toa_distributions{sfx}.png', dpi=150,
                bbox_inches='tight')
    plt.close(fig)

    # resolution comparison bars
    fig, ax = plt.subplots(figsize=(10, 5.5))
    x = np.arange(len(bench))
    ax.bar(x - 0.3, bench['sigma_raw_ns'], 0.3, label='vs trigger (raw)',
           color='lightsteelblue')
    ax.bar(x, bench['sigma_walkcorr_ns'], 0.3,
           label='vs trigger (walk-corrected)', color='steelblue')
    ax.bar(x + 0.3, bench['sigma_padpad_ns'], 0.3,
           label='pad-pad / sqrt(2)  (trigger-free)', color='seagreen')
    ax.set_xticks(x); ax.set_xticklabels(bench['algo'], rotation=30)
    ax.set_ylabel('sigma [ns]')
    ax.set_title(f'{cfg.DET_NAME} timing resolution by algorithm — '
                 f'{cfg.RUN}/{cfg.SUB_RUN}\n'
                 'trigger = 2-PMT scintillator coincidence (5 ns each)')
    ax.legend(); ax.grid(True, alpha=0.3, axis='y')
    fig.tight_layout()
    fig.savefig(f'{out_dir}/resolution_comparison{sfx}.png', dpi=150,
                bbox_inches='tight')
    plt.close(fig)

    # waveform gallery with algorithm markers
    rng = np.random.RandomState(7)
    show = rng.choice(np.nonzero(np.isfinite(t_samp['frac30']))[0],
                      min(12, len(wb)), replace=False)
    fig, axes = plt.subplots(3, 4, figsize=(16, 9), sharex=True)
    tt = np.arange(pw.N_SAMP) * tps
    for ax, k in zip(axes.ravel(), show):
        ax.plot(tt, wb[k], 'k.-', ms=3, lw=1)
        for a, c in [('frac30', 'tab:red'), ('dcfd', 'tab:blue'),
                     ('template', 'tab:green'), ('parabola', 'tab:orange')]:
            v = t_samp[a][k]
            if np.isfinite(v):
                ax.axvline(v * tps, color=c, lw=1, alpha=0.8, label=a)
        ax.set_title(f'amp {hb["amp_max"].iloc[k]:.0f} ADC, '
                     f'ftst {hb["ftst"].iloc[k]}', fontsize=8)
        ax.grid(True, alpha=0.3)
    axes[0, 0].legend(fontsize=7)
    for ax in axes[-1]:
        ax.set_xlabel('time in window [ns]')
    fig.suptitle(f'{cfg.DET_NAME} waveform gallery — {cfg.RUN}/{cfg.SUB_RUN}')
    fig.tight_layout()
    fig.savefig(f'{out_dir}/waveform_gallery{sfx}.png', dpi=130,
                bbox_inches='tight')
    plt.close(fig)

    # ---------------- text summary ----------------
    lines = [
        f'P2 waveform timing — {cfg.DET_TAG} {cfg.RUN}/{cfg.SUB_RUN}',
        f'  decoded files      : {len(files)} (FEUs {sorted(files["feu"].unique())}, '
        f'chunks {sorted(files["chunk"].unique())})',
        f'  benchmark hits     : {len(hb):,} (amp >= {args.amp_min:g} ADC, '
        f'not saturated, <=6 pads/event)',
        f'  time per sample    : {tps:g} ns   ftst convention: {best}',
        f'  trigger reference  : 2-PMT scintillator coincidence (5 ns each)',
        '',
        '  algorithm ranking (sigma vs trigger, walk-corrected | pad-pad/sqrt2):']
    for _, r in bench.iterrows():
        lines.append(f'    {r["algo"]:9s}: {r["sigma_walkcorr_ns"]:6.1f} ns '
                     f'(raw {r["sigma_raw_ns"]:6.1f})  |  '
                     f'{r["sigma_padpad_ns"]:6.1f} ns  '
                     f'({r["n_pairs"]:,} pairs)')
    lines += ['', '  multi-pad cluster combinations '
              f'(best algo = {best_algo}):']
    for _, r in combo.iterrows():
        lines.append(f'    {r["combo"]:18s}: {r["sigma_ns"]:6.1f} ns '
                     f'(q68 {r["q68_ns"]:5.1f})  n={r["n"]:,}')
    txt = '\n'.join(lines)
    print('\n' + txt)
    with open(f'{out_dir}/timing_summary{sfx}.txt', 'w') as f:
        f.write(txt + '\n')
    print(f'\nWritten to: {out_dir}')


if __name__ == '__main__':
    main()
