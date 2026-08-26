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
