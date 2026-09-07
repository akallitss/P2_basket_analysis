#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""threshold_model.py -- does ONE number (the VMM discriminator) turn DREAM's
per-pad spectra into the VMM's per-pad efficiency map?

The comparison in compare_dream_vmm.py established that the two readouts see the
same pad-to-pad gain structure (r = +0.94) but that the VMM records a *smaller*
spread and a much *wider* efficiency spread.  This module tests the obvious
mechanism end to end and, if it holds, costs the fixes.

Model.  A pad's true pulse-height distribution is what DREAM records: its own
threshold sits far below the Landau, and it keeps 96 % of the tracks that point
at it.  The VMM keeps only the part of that same distribution above one common
discriminator level T, identical for all six P2_OUT chips (sdt = 224 on every
one of them).  So

    eff_VMM(pad) = eff_DREAM(pad) x  F_pad(T),      F_pad = fraction above T

with T the single free parameter, fitted in DREAM ADC by track-weighted least
squares against the measured per-pad VMM efficiency.

The same T predicts what truncation does to the measured pulse height, which is
the second half of the story: cutting the DREAM spectra at T reproduces the
VMM's compressed per-pad spread once the VMM ADC's additive offset is removed.
"""
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")


def load():
    g = pd.read_csv(os.path.join(DATA, "compare_dream_vmm_P2_OUT.csv"))
    g = g[g["use"]].reset_index(drop=True)
    z = np.load(os.path.join(DATA,
                             "dream_padadc_eff_nominal_1_P2_OUT_hist.npz"))
    # hist_own = the pulses whose LEADING pad is the pad the track pointed at,
    # i.e. the same quantity the VMM side reports.
    H = z["hist_own"][g["di"].to_numpy()].astype(float)
    return g, H, float(z["amp_bin"])


class Spectra:
    """The 53 per-pad DREAM spectra, with interpolating cuts."""

    def __init__(self, H, binw):
        self.H, self.binw = H, binw
        self.tot = H.sum(1)
        self.cum = np.cumsum(H, 1)

    def frac_above(self, T):
        b = T / self.binw
        i = int(np.floor(b))
        below = self.cum[:, i - 1] if i > 0 else 0.0
        return 1.0 - (below + self.H[:, i] * (b - i)) / self.tot

    def median_above(self, T):
        """Median of what survives the cut -- the truncated estimator the VMM
        is actually reporting."""
        b0 = int(np.ceil(T / self.binw))
        out = []
        for h in self.H[:, b0:]:
            c = np.cumsum(h) / h.sum()
            i = int(np.searchsorted(c, 0.5, side="left"))
            lo = c[i - 1] if i > 0 else 0.0
            out.append((i + (0.5 - lo) / max(c[i] - lo, 1e-12)) * self.binw)
        return np.array(out) + b0 * self.binw

    def mpv(self, lo=0.0, hi=None):
        """Pooled Landau peak.  Blanked above `hi` because the top bin of both
        readouts is an overflow pile-up, and a raw argmax lands on it."""
        h = self.H.sum(0).copy()
        h[:int(lo / self.binw)] = 0
        if hi is not None:
            h[int(hi / self.binw):] = 0
        k = np.convolve(h, np.ones(5) / 5, "same")
        return float((np.argmax(k) + 0.5) * self.binw)


def relrms(v):
    return float(np.std(np.asarray(v) / np.median(v), ddof=1))


def _pearson_p(r, n):
    """Two-sided p for a Pearson r on n points.  scipy is imported lazily --
    the rest of this module is numpy/pandas only and is run on lxplus."""
    from scipy import stats
    dof = max(n - 2, 1)
    t = r * np.sqrt(dof / max(1.0 - r * r, 1e-12))
    return float(2.0 * stats.t.sf(abs(t), dof))


NBAND = 8       # the gain bands slide_ridge.png averages its rows over


def gain_bands(g, nband=NBAND):
    """Band index per pad, pads sorted by DREAM gain.  Lives here rather than
    in the figure module because the per-pad summary below has to decompose
    against exactly the split the eight-row figure averages over."""
    order = np.argsort(g["amp_med_d"].to_numpy())
    b = np.empty(len(g), int)
    for k, idx in enumerate(np.array_split(order, nband)):
        b[idx] = k
    return b


def fit_threshold(g, S):
    effv, effd = g["eff_v"].to_numpy(), g["eff_d"].to_numpy()
    w = g["n_track_v"].to_numpy().astype(float)
    Ts = np.arange(4.0, 400.0, 1.0)
    chi = np.array([np.sum(w * (effd * S.frac_above(T) - effv) ** 2) / w.sum()
                    for T in Ts])
    T = float(Ts[np.argmin(chi)])
    # 1-sigma-ish band: where the weighted rms residual grows by 10 %
    ok = Ts[chi <= chi.min() * 1.21]
    return T, float(ok.min()), float(ok.max()), Ts, np.sqrt(chi)


def fit_threshold_per_pad(g, S, Ts=None):
    """The same model, inverted PAD BY PAD instead of fit jointly: solve
    eff_v(pad) = eff_d(pad) x frac_above(T_pad) for T_pad using only that
    pad's own spectrum and its own two measured efficiencies.

    This is the test of the global-threshold claim, not a restatement of it:
    fit_threshold() can always find SOME T that minimises a joint residual --
    that's just least squares. If the mechanism is really one common
    discriminator, the 53 independent T_pad should cluster tightly around the
    global T; if something else is going on (per-channel noise, timewalk,
    gain-dependent effects the model doesn't have) they'll scatter or trend
    with gain instead.

    Each T_pad carries an approximate statistical uncertainty, so the
    comparison to the global fit can be judged against counting noise rather
    than by eye alone.  Two independent sources go in, both propagated through
    the local slope of frac_above at T_pad:
      * the finite track counts behind eff_v and eff_d (binomial), and
      * the finite number of pulses in the pad's OWN DREAM histogram, which
        makes frac_above(T) itself a binomial estimate.  The weakest pads have
        ~1e3 pulses, so this is not negligible against the first term.
    The slope is taken over three bins rather than one: a single 8-ADC bin of
    one pad's spectrum is noisy enough to distort sigma_T on its own.

    Returns
    -------
    Tpad, sigT : per-pad threshold and its ~1 sigma uncertainty (nan where
        undefined -- i.e. a boundary-clipped pad, see below)
    lo_clip, hi_clip : boundary flags. lo_clip = pad's eff_v/eff_d ratio is
        higher than frac_above(Ts[0]) -- no threshold in the scanned range is
        low enough (usually eff_v ~ eff_d, a noise fluctuation, not a real
        low threshold). hi_clip = the opposite, ratio below frac_above(Ts[-1])
        -- would need an implausibly high threshold.
    """
    if Ts is None:
        Ts = np.arange(4.0, 400.0, 1.0)
    effv, effd = g["eff_v"].to_numpy(), g["eff_d"].to_numpy()
    nv = g["n_track_v"].to_numpy().astype(float)
    nd = g["n_track_d"].to_numpy().astype(float)
    target = np.clip(effv / np.clip(effd, 1e-9, None), 0.0, 1.0)

    # binomial errors on the two efficiencies, propagated onto the ratio
    sv = np.sqrt(np.clip(effv * (1 - effv), 0, None) / np.clip(nv, 1, None))
    sd = np.sqrt(np.clip(effd * (1 - effd), 0, None) / np.clip(nd, 1, None))
    starget = target * np.sqrt((sv / np.clip(effv, 1e-6, None)) ** 2
                               + (sd / np.clip(effd, 1e-6, None)) ** 2)

    F = np.array([S.frac_above(T) for T in Ts])   # (nT, npad), decreasing in T
    npad = effv.size
    Tpad = np.empty(npad)
    sigT = np.full(npad, np.nan)
    lo_clip = np.zeros(npad, bool)
    hi_clip = np.zeros(npad, bool)
    for j in range(npad):
        f = F[:, j]
        if target[j] >= f[0]:
            Tpad[j] = Ts[0]; lo_clip[j] = True
            continue
        if target[j] <= f[-1]:
            Tpad[j] = Ts[-1]; hi_clip[j] = True
            continue
        Tj = np.interp(target[j], f[::-1], Ts[::-1])
        Tpad[j] = Tj
        # frac_above at the solution is a binomial estimate over the pad's own
        # pulses; that error adds to the two efficiencies' in quadrature
        fj = float(np.interp(Tj, Ts, f))
        s_hist = np.sqrt(max(fj * (1.0 - fj), 0.0) / max(S.tot[j], 1.0))
        s = np.hypot(starget[j], s_hist)
        # local slope of frac_above at T_j is minus the (normalised) pad
        # spectrum density there -- convert the ratio's error into a T error
        i = int(np.clip(Tj / S.binw, 0, S.H.shape[1] - 1))
        dens = S.H[j, max(i - 1, 0):i + 2].mean() / S.binw / max(S.tot[j], 1.0)
        sigT[j] = s / dens if dens > 1e-9 else np.nan
    return Tpad, sigT, lo_clip, hi_clip


def main():
    g, H, binw = load()
    S = Spectra(H, binw)
    effv, effd = g["eff_v"].to_numpy(), g["eff_d"].to_numpy()
    w = g["n_track_v"].to_numpy().astype(float)

    T, Tlo, Thi, Ts, rmsc = fit_threshold(g, S)
    pred = effd * S.frac_above(T)
    n = dict(T=T, T_lo=Tlo, T_hi=Thi,
             mpv_dream=S.mpv(hi=3660.0),
             eff_pred_all=float((pred * w).sum() / w.sum()),
             eff_obs_all=float((effv * w).sum() / w.sum()),
             eff_dream_all=float((effd * g["n_track_d"]).sum()
                                 / g["n_track_d"].sum()),
             r_pred=float(np.corrcoef(pred, effv)[0, 1]),
             resid_rms=float(np.sqrt(np.average((pred - effv) ** 2,
                                                weights=w))),
             npad=int(len(g)))
    n["T_over_mpv"] = T / n["mpv_dream"]
    # ...and what that same level is in units of the WEAKEST pad's own peak
    mpv_pad = g["amp_med_d"].to_numpy() / np.median(g["amp_med_d"]) \
        * n["mpv_dream"]
    n["T_over_mpv_weak"] = float(T / np.percentile(mpv_pad, 10))

    # --- the pulse-height half: truncation + an additive ADC offset --------- #
    ampv = g["amp_med_v"].to_numpy()
    mt = S.median_above(T)
    b = np.polyfit(mt, ampv, 1)
    n["rms_dream_raw"] = relrms(g["amp_med_d"])
    n["rms_dream_cut"] = relrms(mt)
    n["rms_vmm"] = relrms(ampv)
    n["offset_adc"] = float(b[1])
    n["gain_adc"] = float(b[0])
    n["rms_vmm_deoffset"] = relrms(ampv - b[1])
    n["r_cut"] = float(np.corrcoef(mt, ampv)[0, 1])

    # --- the same model inverted pad by pad, as a test of the global fit ----- #
    # A joint least squares always returns SOME T; whether one common
    # discriminator is really the mechanism is a question about the 53
    # independent answers, so summarise them here and let the report quote it.
    Tpad, sigT, lo_c, hi_c = fit_threshold_per_pad(g, S)
    ok = ~(lo_c | hi_c) & np.isfinite(sigT)
    b = gain_bands(g)
    tp, bb = Tpad[ok], b[ok]
    betw = sum((bb == q).sum() * (tp[bb == q].mean() - tp.mean()) ** 2
               for q in range(NBAND) if (bb == q).any())
    with_ = sum(((tp[bb == q] - tp[bb == q].mean()) ** 2).sum()
                for q in range(NBAND) if (bb == q).any())
    # how many ADC of threshold one efficiency point is worth, per pad: the
    # exchange rate that turns the T spread below into the residual above
    i0 = np.clip((Tpad / binw).astype(int), 1, H.shape[1] - 2)
    dens = np.array([H[j, i - 1:i + 2].mean() / binw / max(S.tot[j], 1.0)
                     for j, i in enumerate(i0)])
    n["perpad"] = dict(
        n_ok=int(ok.sum()), n_clip=int((~ok).sum()),
        T_med=float(np.median(Tpad[ok])),
        T_iqr_lo=float(np.percentile(Tpad[ok], 25)),
        T_iqr_hi=float(np.percentile(Tpad[ok], 75)),
        T_rms=float(np.std(Tpad[ok], ddof=1)),
        T_relrms=float(np.std(Tpad[ok], ddof=1) / np.median(Tpad[ok])),
        sig_med=float(np.median(sigT[ok])),
        chi2_dof=float(np.sum(((Tpad[ok] - T) / sigT[ok]) ** 2)
                       / max(int(ok.sum()) - 1, 1)),
        frac_within_band=float(with_ / (with_ + betw)),
        adc_per_point=float(np.median(0.01 / (effd[ok] * dens[ok]))),
        # a dispersion that TRACKED gain would be a different mechanism (and
        # would matter for the fixes below); this one does not, at p ~ 0.1
        r_gain=float(np.corrcoef(np.log(g["amp_med_d"].to_numpy()[ok]),
                                 Tpad[ok])[0, 1]))
    n["perpad"]["p_gain"] = float(_pearson_p(n["perpad"]["r_gain"],
                                             int(ok.sum())))

    # --- costing the fixes -------------------------------------------------- #
    # Lowering the discriminator and raising the chamber gain are the SAME move
    # in this model: both are the ratio T/gain.  Report it as a gain factor.
    rows = []
    for f in [1.0, 1.15, 1.3, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0]:
        p = effd * S.frac_above(T / f)
        rows.append(dict(factor=f,
                         eff_all=float((p * w).sum() / w.sum()),
                         eff_min=float(np.sort(p)[1]),   # skip the dead pad
                         eff_p10=float(np.percentile(p, 10)),
                         spread=float(p.max() - np.sort(p)[1])))
    n["fix"] = rows
    n["factor_for_95"] = float(np.interp(
        0.95, [r["eff_all"] for r in rows], [r["factor"] for r in rows]))

    out = os.path.join(DATA, "threshold_model_P2_OUT.json")
    json.dump(n, open(out, "w"), indent=1)
    print(json.dumps(n, indent=1))
    print("wrote", out)


if __name__ == "__main__":
    main()
