#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Why do the P2 chambers read out by the VMM miss muons? Diagnostics on the scintillator-confirmed tags of
efficiency_tp.py --dump-diag dumps (channel_efficiency.py dump --diag), together with the HV-monitor bad times.

Inputs: TAGS_DIR (tags_run<pid>.npz + eff_run<pid>.json, pid = run x 100 + chunk), PSEUDO_DIR/chunks.json (pid -> run,
configuration), VETO_DIR/veto_run_<run>.json (spark_veto.py). Every quantity is computed on ON-time minus r x OFF-time
tags (r = r_scale of the chunk), i.e. for muons, with the scintillator window --scint (as channel_efficiency measure).

  1. sparks   are the HV-monitor flags real? Per flag class (hv_spark / hv_small / hv_vdev / lv_glitch, electrode):
              chamber hits/s and saturated hits/s (ADC >= 1000) around the flag against random times, from the
              per-second rates of the dump. Data-side discharges: seconds whose saturated-hit count is far above the
              run's median; how many fall inside the HV veto, and how many HV flags have one nearby.
  2. eff      a simple nearest-hit efficiency (hit = live probe pad within --probe-r mm and +-win ns: efficiency_tp's
              definition) with and without the HV veto, per configuration and chamber. The official counting
              estimator with the veto is channel_efficiency.py measure --bad-times.
  3. adc      pulse height of the probe hits that make the efficiency: per muon the largest ADC of the hits within
              --probe-r mm in the window, and the predicted pad's own hit ADC; quantiles, saturated fraction,
              spectra per configuration / chamber, medians per chip.
  4. timing   time of the probe hits within --probe-r mm, from tag + latency, over +-2 us: spectrum, fraction outside the
              window used for the efficiency, efficiency vs window half-width (random-hit level from 1.5-2 us).
  5. misses   every muon without a hit, by the first matching cause:
                hv_bad_time      inside an HV/LV veto interval
                discharge        a saturated hit (ADC >= 1000) in the probe chamber within +-100 us
                out_of_time      a hit within --probe-r mm, but outside +-win (within +-2 us)
                own_dead_time    the predicted pad's channel fired < --dead-us us before
                next_to_masked   the prediction is within --probe-r mm of a dead/masked pad
                chamber_silent   no live hit anywhere in the probe chamber within +-win
                near_miss        hits in the window, nearest 15-40 mm away (tracking tails, map, charge spread)
                far_only         hits in the window, all > 40 mm away (unrelated activity)
              per configuration and chamber, and per chip.

Usage: probe_diagnostics.py TAGS_DIR PSEUDO_DIR VETO_DIR OUT_DIR [--exclude-runs 233] [--scint 140:310]
Outputs OUT_DIR/*.csv, *.png, summary.txt (lxplus, LCG view; light: reads the dumps only).
"""
import argparse
import glob
import json
import os
import re
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
NAME = {3: 'P2_IN', 2: 'P2_MID', 1: 'P2_OUT'}
DETV = {3: range(2, 8), 2: range(8, 14), 1: range(14, 20)}
ELEC_DET = {'IN': 'P2_IN', 'MID': 'P2_MID', 'OUT': 'P2_OUT'}
CAUSES = ['hv_bad_time', 'discharge', 'out_of_time', 'own_dead_time', 'next_to_masked', 'chamber_silent',
          'near_miss', 'far_only']
WINDOWS = [25, 50, 75, 100, 150, 200, 300, 500, 1000, 1500]


def label(cfg):
    try:
        from common import label as L
        return L(cfg)
    except ImportError:
        return cfg


def real_run(pid):
    return pid if pid < 1000 else (pid // 1000 if pid >= 100000 else pid // 100)


def load_veto_json(vdir, run):
    f = f'{vdir}/veto_run_{run}.json'
    return json.load(open(f)) if os.path.exists(f) else None


def in_iv(t, iv):
    if not len(iv):
        return np.zeros(len(t), bool)
    k = np.searchsorted(iv[:, 0], t, side='right') - 1
    return (k >= 0) & (t <= iv[np.maximum(k, 0), 1])


def wmean(x, w):
    return float(np.sum(x * w) / np.sum(w)) if np.sum(w) else np.nan


def wquant(x, w, qs):
    """quantiles of x with (possibly negative: on - r off) weights; negative weights clipped after cumulation."""
    if not len(x) or np.sum(w) <= 0:
        return [np.nan] * len(qs)
    o = np.argsort(x)
    c = np.cumsum(w[o])
    c = np.maximum.accumulate(np.clip(c, 0, None)) / c[-1]
    return [float(x[o][min(np.searchsorted(c, q), len(x) - 1)]) for q in qs]


# ------------------------------------------------------------------------------------------------------------ loading
def load_all(a):
    chunks = {int(k): v for k, v in json.load(open(f'{a.pseudo}/chunks.json')).items()}
    excl = {int(x) for x in a.exclude_runs.split(',') if x}
    from channel_efficiency import pad_geometry, run_info
    X, Y, _ = pad_geometry()
    lo, hi = map(float, a.scint.split(':'))
    T, NR, RATES = [], [], []
    vetoes = {}
    for fz in sorted(glob.glob(f'{a.tags}/tags_run*.npz')):
        pid = int(re.search(r'run(\d+)\.npz', fz)[1])
        if pid not in chunks:
            continue
        run = chunks[pid]['run']
        if run in excl:
            continue
        Z = np.load(fz)
        if 'p3_wall' not in Z.files:
            print(f'{fz}: not a --dump-diag dump, skipped')
            continue
        cfg = label(chunks[pid]['cfg'])
        if run not in vetoes:
            vj = load_veto_json(a.veto, run)
            vetoes[run] = (np.array([[x, y] for x, y, _ in vj['intervals']], float).reshape(-1, 2) if vj else
                           np.zeros((0, 2)), vj)
        iv = vetoes[run][0]
        info = run_info(a.tags, pid)
        if 'rates_t0' in Z.files:
            RATES.append(dict(run=run, pid=pid, t0=float(Z['rates_t0']),
                              **{f'{k}_{d}': Z[f'rates_{k}_p{d}'] for k in ('hits', 'sat') for d in (3, 2, 1)}))
        for p in (3, 2, 1):
            g = lambda k: Z[f'p{p}_{k}']
            dts = g('dt_scint')
            sc = (dts >= lo) & (dts <= hi)
            n = len(dts)
            idx = np.nonzero(sc)[0]
            if not len(idx):
                continue
            r = float(g('r_scale'))
            w = np.where(g('on')[idx], 1.0, -r)
            vmm, ch = g('vmm')[idx], g('ch')[idx]
            # dead/masked pads of the probe chamber within --probe-r of the prediction
            dead = [(v, c) for v, c in info['dead'] if v in DETV[p] and np.isfinite(X[v, c])]
            px, py = g('px')[idx], g('py')[idx]
            if dead:
                dxy = np.array([(X[v, c], Y[v, c]) for v, c in dead])
                from scipy.spatial import cKDTree
                ddead = cKDTree(dxy).query(np.c_[px, py])[0]
            else:
                ddead = np.full(len(idx), np.inf)
            wall = g('wall')[idx]
            T.append(pd.DataFrame(dict(
                run=run, pid=pid, cfg=cfg, det=NAME[p], w=w, on=g('on')[idx], vmm=vmm, ch=ch,
                key=vmm * 64 + ch, wall=wall, bad=in_iv(wall, iv), d_win=g('d_win')[idx], n_any=g('n_any')[idx],
                n100us=g('n100us')[idx], sat100us=g('sat100us')[idx], prev_chan=g('prev_chan')[idx],
                prev_chip=g('prev_chip')[idx], d_dead=ddead, tid=idx + p * 10_000_000 + pid * 100_000_000)))
            # near hits of the selected tags, re-indexed to the global tag id
            fi = g('fired_nr_i')
            keep = sc[fi] if len(fi) else np.zeros(0, bool)
            if keep.any():
                NR.append(pd.DataFrame(dict(
                    tid=fi[keep] + p * 10_000_000 + pid * 100_000_000, dt=g('fired_nr_dt')[keep],
                    adc=g('fired_nr_adc')[keep], d=g('fired_nr_d')[keep], hkey=g('fired_nr_key')[keep])))
            win = float(g('win'))
    tags = pd.concat(T, ignore_index=True)
    near = pd.concat(NR, ignore_index=True) if NR else pd.DataFrame(columns=['tid', 'dt', 'adc', 'd', 'hkey'])
    near = near.merge(tags[['tid', 'cfg', 'det', 'w', 'key', 'bad', 'vmm', 'run']], on='tid')
    return tags, near, RATES, vetoes, win


# --------------------------------------------------------------------------------------------------------- 1. sparks
def spark_check(RATES, vetoes, out, lag_max=60):
    """per flag class: mean saturated / all hits per second per chamber around the flag, vs random times."""
    rows, prof, data_sparks = [], {}, []
    by_run = {}
    for R in RATES:
        by_run.setdefault(R['run'], []).append(R)
    rng = np.random.default_rng(1)
    for run, Rs in sorted(by_run.items()):
        t0 = min(R['t0'] for R in Rs)
        t1 = max(R['t0'] + len(R['hits_3']) for R in Rs)
        nsec = int(t1 - t0) + 1
        S = {f'{k}_{d}': np.full(nsec, np.nan) for k in ('hits', 'sat') for d in (3, 2, 1)}
        for R in Rs:
            o = int(R['t0'] - t0)
            for k in S:
                v = R[k].astype(float)
                seg = S[k][o:o + len(v)]
                S[k][o:o + len(v)] = np.where(np.isnan(seg), v, seg + v)   # chunks share a second at their edge
        vj = vetoes.get(run, (None, None))[1]
        if vj is None:
            continue
        iv = vetoes[run][0]
        flags = [(t, w, el) for t, w, el, _ in vj['flags']]
        cls = {}
        for t, w, el in flags:
            part = el.split('_')[1] if el.startswith('P2_') else el.split()[0]
            elec = el.split('_')[2] if el.startswith('P2_') else ''
            cls.setdefault(f'{w} {part} {elec}'.strip(), []).append(t)
        cls['random'] = list(rng.uniform(t0 + lag_max, t1 - lag_max, 200))
        lag = np.arange(-lag_max, lag_max + 1)
        for c, ts in cls.items():
            for t in ts:
                k = int(np.floor(t - t0))
                if k + lag[0] < 0 or k + lag[-1] >= nsec:
                    continue
                for d in (3, 2, 1):
                    for q in ('hits', 'sat'):
                        prof.setdefault((c, q, d), []).append(S[f'{q}_{d}'][k + lag])
        # data-side discharges: seconds with saturated hits far above the run median (per chamber)
        for d in (3, 2, 1):
            s = S[f'sat_{d}']
            ok = np.isfinite(s)
            med = np.median(s[ok])
            mad = np.median(np.abs(s[ok] - med)) * 1.4826 + 1
            cut = med + 10 * mad
            ks = np.nonzero(ok & (s > cut))[0]
            for k in ks:
                t = t0 + k + 0.5
                near_flag = min((abs(t - tf) for tf, _, _ in flags), default=np.inf)
                data_sparks.append(dict(run=run, det=NAME[d], t=t, sat=int(s[k]), median=float(med), cut=float(cut),
                                        in_hv_veto=bool(in_iv(np.array([t]), iv)[0]), dt_nearest_flag=near_flag))
    lag = np.arange(-lag_max, lag_max + 1)
    near = np.abs(lag) <= 10
    for (c, q, d), L in prof.items():
        if q != 'sat':
            continue
        A = np.nanmean(np.array(L), 0)
        H = np.nanmean(np.array(prof[(c, 'hits', d)]), 0)
        Rr = np.nanmean(np.array(prof[('random', 'sat', d)]), 0)
        Rh = np.nanmean(np.array(prof[('random', 'hits', d)]), 0)
        rows.append(dict(cls=c, det=NAME[d], n=len(L), sat_peak=float(np.nanmax(A[near])),
                         sat_random=float(np.nanmean(Rr)), sat_ratio=float(np.nanmax(A[near]) / max(np.nanmean(Rr), 1e-9)),
                         hits_peak=float(np.nanmax(H[near])), hits_random=float(np.nanmean(Rh)),
                         hits_ratio=float(np.nanmax(H[near]) / max(np.nanmean(Rh), 1e-9))))
    F = pd.DataFrame(rows).sort_values(['cls', 'det'])
    F.to_csv(f'{out}/spark_flag_check.csv', index=False)
    DS = pd.DataFrame(data_sparks)
    DS.to_csv(f'{out}/data_discharges.csv', index=False)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    classes = sorted({c for c, _, _ in prof})
    fig, axes = plt.subplots(2, 3, figsize=(16, 8), sharex=True)
    for j, d in enumerate((3, 2, 1)):
        for i, q in enumerate(('hits', 'sat')):
            ax = axes[i, j]
            for c in classes:
                L = prof.get((c, q, d))
                if not L or (len(L) < 3 and c != 'random'):
                    continue
                ax.plot(lag, np.nanmean(np.array(L), 0), lw=2 if c == 'random' else 1,
                        color='k' if c == 'random' else None, label=f'{c} ({len(L)})')
            ax.set_yscale('log')
            ax.axvline(0, color='0.6', lw=0.8)
            ax.grid(alpha=0.3)
            ax.set_title(f'{NAME[d]}: {"all hits" if q == "hits" else "saturated hits (ADC >= 1000)"} per s')
    axes[0, 0].legend(fontsize=6, ncol=2)
    axes[1, 1].set_xlabel('time from the HV/LV monitor flag [s]')
    fig.suptitle('chamber activity around the HV / LV monitor flags (all runs), against random times')
    fig.tight_layout()
    fig.savefig(f'{out}/spark_flag_check.png', dpi=100)
    plt.close(fig)
    return F, DS


# ------------------------------------------------------------------------------------------------------------ 2. eff
def eff_table(tags, probe_r):
    rows = []
    for (cfg, det), g in tags.groupby(['cfg', 'det'], sort=False):
        for veto in (False, True):
            m = ~g.bad if veto else np.ones(len(g), bool)
            gg = g[m]
            n = gg.w.sum()
            k = (gg.w * (gg.d_win < probe_r)).sum()
            n_on = int(gg.on.sum())
            rows.append(dict(cfg=cfg, det=det, veto=veto, muons=round(n, 1), eff=k / n if n > 0 else np.nan,
                             err=np.sqrt(max(k / n * (1 - k / n), 0) / n) if n > 0 else np.nan, on_tags=n_on,
                             vetoed_on=int((g.on & g.bad).sum()) if veto else 0))
    E = pd.DataFrame(rows)
    return E


# ------------------------------------------------------------------------------------------------------------ 3. adc
def adc_tables(near, tags, probe_r, win, out):
    inw = near[(near.d < probe_r) & (near.dt.abs() <= win) & ~near.bad]
    mx = inw.groupby('tid').agg(adc_max=('adc', 'max')).reset_index()
    own = inw[inw.hkey == inw.key].groupby('tid').agg(adc_own=('adc', 'max')).reset_index()
    t = tags[~tags.bad][['tid', 'cfg', 'det', 'w', 'vmm']].merge(mx, on='tid').merge(own, on='tid', how='left')
    rows = []
    bins = np.arange(0, 1040, 16)
    spectra = {}
    for (cfg, det), g in t.groupby(['cfg', 'det'], sort=False):
        q = wquant(g.adc_max.to_numpy(float), g.w.to_numpy(), [0.05, 0.1, 0.25, 0.5, 0.75, 0.9])
        go = g[g.adc_own.notna()]
        qo = wquant(go.adc_own.to_numpy(float), go.w.to_numpy(), [0.1, 0.5, 0.9])
        rows.append(dict(cfg=cfg, det=det, muons_with_hit=round(g.w.sum(), 1),
                         q05=q[0], q10=q[1], q25=q[2], median=q[3], q75=q[4], q90=q[5],
                         sat_frac=wmean((g.adc_max >= 1000).to_numpy(float), g.w.to_numpy()),
                         below100=wmean((g.adc_max < 100).to_numpy(float), g.w.to_numpy()),
                         own_q10=qo[0], own_median=qo[1], own_q90=qo[2]))
        spectra[(cfg, det)] = np.histogram(g.adc_max, bins, weights=g.w)[0]
    A = pd.DataFrame(rows)
    chip = t.groupby(['cfg', 'det', 'vmm']).apply(
        lambda g: pd.Series(dict(muons=round(g.w.sum(), 1), median=wquant(g.adc_max.to_numpy(float), g.w.to_numpy(),
                                                                           [0.5])[0],
                                 sat_frac=wmean((g.adc_max >= 1000).to_numpy(float), g.w.to_numpy())))).reset_index()
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(17, 5), sharey=True)
    for ax, det in zip(axes, ('P2_IN', 'P2_MID', 'P2_OUT')):
        for (cfg, d), h in spectra.items():
            if d != det:
                continue
            ax.step(bins[:-1], h / max(h.sum(), 1), where='post', lw=1, label=cfg)
        ax.set_title(f'{det}: largest ADC of the probe hits within {probe_r:.0f} mm (muons, HV veto applied)')
        ax.set_xlabel('ADC')
        ax.grid(alpha=0.3)
    axes[0].set_ylabel('fraction of muons with a hit per 16 ADC')
    axes[0].legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(f'{out}/adc_spectra.png', dpi=100)
    plt.close(fig)
    return A, chip


# --------------------------------------------------------------------------------------------------------- 4. timing
def timing_tables(near, tags, probe_r, win, out):
    nr = near[(near.d < probe_r) & ~near.bad]
    tg = tags[~tags.bad]
    rows, hists = [], {}
    edges = np.arange(-2000, 2001, 25)
    for (cfg, det), g in tg.groupby(['cfg', 'det'], sort=False):
        h = nr[(nr.cfg == cfg) & (nr.det == det)]
        nmu = g.w.sum()
        # random-hit density from 1.5-2 us (per ns per muon), to subtract from the window counts
        side = h[h.dt.abs() >= 1500]
        dens = side.groupby('tid').size().reindex(g.tid).fillna(0).to_numpy() @ g.w.to_numpy() / nmu / 1000.0
        hists[(cfg, det)] = np.histogram(h.dt, edges, weights=h.w)[0]
        rec = dict(cfg=cfg, det=det, muons=round(nmu, 1), random_hits_per_us=dens * 1000)
        best = h.assign(a=h.dt.abs()).sort_values('a').drop_duplicates('tid')
        bw = best.set_index('tid').a.reindex(g.tid).to_numpy()
        for W in WINDOWS:
            k = ((bw <= W) * g.w.to_numpy()).sum()
            acc = min(1.0, 2 * W * dens)
            rec[f'eff_w{W}'] = ((k / nmu) - acc) / (1 - acc) if nmu > 0 else np.nan
        inside = (h.dt.abs() <= win)
        out_ = (h.dt.abs() > win) & (h.dt.abs() < 1500)
        rec['hits_in_win'] = float((h.w * inside).sum())
        rec['hits_out_win_lt1500'] = float((h.w * out_).sum() - 2 * (1500 - win) * dens * nmu)
        q = wquant(h.dt[inside].to_numpy(float), h.w[inside].to_numpy(), [0.05, 0.5, 0.95])
        rec.update(dt_q05=q[0], dt_median=q[1], dt_q95=q[2])
        rows.append(rec)
    Tm = pd.DataFrame(rows)
    per_run = nr[nr.dt.abs() <= win].groupby(['run', 'det']).apply(
        lambda g: wquant(g.dt.to_numpy(float), g.w.to_numpy(), [0.5])[0]).rename('dt_median').reset_index()
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(17, 5), sharey=True)
    c = 0.5 * (edges[1:] + edges[:-1])
    for ax, det in zip(axes, ('P2_IN', 'P2_MID', 'P2_OUT')):
        for (cfg, d), h in hists.items():
            if d == det:
                ax.step(c, h / max(h.sum(), 1), where='mid', lw=1, label=cfg)
        ax.axvspan(-win, win, color='0.85', zorder=0)
        ax.set_yscale('log')
        ax.set_xlabel('probe hit time - (tag + latency) [ns]')
        ax.set_title(f'{det}: hits within {probe_r:.0f} mm of the prediction (grey: window used)')
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(f'{out}/timing.png', dpi=100)
    plt.close(fig)
    return Tm, per_run


# --------------------------------------------------------------------------------------------------------- 5. misses
def classify(tags, near, probe_r, win, dead_us):
    nr = near[near.d < probe_r]
    has_in = set(nr.tid[nr.dt.abs() <= win])
    has_out = set(nr.tid[nr.dt.abs() > win])
    t = tags.copy()
    hit = t.tid.isin(has_in) | (t.d_win < probe_r)
    cause = np.full(len(t), '', object)
    rules = [('hv_bad_time', t.bad.to_numpy()),
             ('discharge', (t.sat100us > 0).to_numpy()),
             ('out_of_time', t.tid.isin(has_out).to_numpy()),
             ('own_dead_time', (t.prev_chan < dead_us * 1e3).to_numpy()),
             ('next_to_masked', (t.d_dead < probe_r).to_numpy()),
             ('chamber_silent', (t.n_any == 0).to_numpy()),
             ('near_miss', (t.d_win <= 40).to_numpy()),
             ('far_only', np.ones(len(t), bool))]
    left = ~hit.to_numpy()
    for name, m in rules:
        sel = left & m
        cause[sel] = name
        left &= ~m
    cause[hit.to_numpy()] = 'hit'
    t['cause'] = cause

    def table(keys):
        rows = []
        for key, g in t.groupby(keys, sort=False):
            n = g.w.sum()
            rec = dict(zip(keys, key if isinstance(key, tuple) else (key,)))
            rec['muons'] = round(n, 1)
            rec['eff'] = (g.w * (g.cause == 'hit')).sum() / n if n > 0 else np.nan
            for c in CAUSES:
                rec[c] = (g.w * (g.cause == c)).sum() / n if n > 0 else np.nan
            rows.append(rec)
        return pd.DataFrame(rows)
    return table(['cfg', 'det']), table(['det', 'vmm']), table(['cfg', 'det', 'vmm'])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('tags')
    ap.add_argument('pseudo')
    ap.add_argument('veto')
    ap.add_argument('out')
    ap.add_argument('--exclude-runs', default='233')
    ap.add_argument('--scint', default='140:310')
    ap.add_argument('--probe-r', type=float, default=15)
    ap.add_argument('--dead-us', type=float, default=10)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    tags, near, RATES, vetoes, win = load_all(a)
    pd.set_option('display.width', 250)
    pd.set_option('display.max_columns', 40)
    S = open(f'{a.out}/summary.txt', 'w')

    def say(*x):
        print(*x)
        print(*x, file=S)
    say(f'{tags.run.nunique()} runs, {tags.pid.nunique()} chunks, {int(tags.on.sum())} on-time scintillator tags; '
        f'window +-{win:.0f} ns, probe radius {a.probe_r:.0f} mm')
    F, DS = spark_check(RATES, vetoes, a.out)
    say('\n== 1. HV/LV flags vs chamber activity (saturated hits/s within +-10 s of the flag, peak / random)')
    say(F.to_string(index=False, float_format=lambda x: f'{x:.3g}'))
    if len(DS):
        say(f'\ndata-side discharges (seconds with saturated hits > median + 10 MAD): {len(DS)}; inside the HV veto '
            f'{DS.in_hv_veto.mean() * 100:.0f} %; an HV flag within 10 s: {(DS.dt_nearest_flag <= 10).mean() * 100:.0f} %')
        say(DS.groupby('det').agg(n=('t', 'size'), in_veto=('in_hv_veto', 'mean'),
                                  flag_10s=('dt_nearest_flag', lambda x: (x <= 10).mean())).to_string())
    E = eff_table(tags, a.probe_r)
    E.to_csv(f'{a.out}/eff_veto.csv', index=False)
    say('\n== 2. nearest-hit efficiency (muons = on - r off), without / with the HV veto')
    P = E.pivot_table(index=['cfg', 'det'], columns='veto', values='eff', sort=False)
    P.columns = ['all', 'veto']
    P['delta'] = P['veto'] - P['all']
    V = E[E.veto].set_index(['cfg', 'det'])
    P['vetoed_on_frac'] = V.vetoed_on / (V.on_tags + V.vetoed_on)
    P['err'] = V.err
    say(P.to_string(float_format=lambda x: f'{x:.3f}'))
    A, chipA = adc_tables(near, tags, a.probe_r, win, a.out)
    A.to_csv(f'{a.out}/adc.csv', index=False)
    chipA.to_csv(f'{a.out}/adc_per_chip.csv', index=False)
    say('\n== 3. pulse height of the hit muons (largest ADC within the probe radius and window; own = predicted pad)')
    say(A.to_string(index=False, float_format=lambda x: f'{x:.3g}'))
    Tm, per_run = timing_tables(near, tags, a.probe_r, win, a.out)
    Tm.to_csv(f'{a.out}/timing.csv', index=False)
    per_run.to_csv(f'{a.out}/timing_per_run.csv', index=False)
    say('\n== 4. timing: efficiency vs window half-width (ns), random hits subtracted; hits outside the window')
    say(Tm.to_string(index=False, float_format=lambda x: f'{x:.3g}'))
    C, Cchip, Ccfgchip = classify(tags, near, a.probe_r, win, a.dead_us)
    C.to_csv(f'{a.out}/misses.csv', index=False)
    Cchip.to_csv(f'{a.out}/misses_per_chip.csv', index=False)
    Ccfgchip.to_csv(f'{a.out}/misses_per_cfg_chip.csv', index=False)
    say('\n== 5. muons without a hit, by first matching cause (fractions of all muons; eff + causes = 1)')
    say(C.to_string(index=False, float_format=lambda x: f'{x:.3f}'))
    say('\nper chip (all configurations):')
    say(Cchip.sort_values(['det', 'vmm']).to_string(index=False, float_format=lambda x: f'{x:.3f}'))
    S.close()


if __name__ == '__main__':
    main()
