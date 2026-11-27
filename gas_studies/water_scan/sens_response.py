#!/usr/bin/env python3
"""sens_response.py -- systematics of the pulse-response MC (att = 0).
For each variant: c(T) = M(T)/T at a few T, and the full drift time implied by the chamber
data y = median - t0 at dV = 25, 50, 75, 100, 200, 300 V (dt0 = 0).
Output results/response_systematics.csv"""
import os
import numpy as np
import pandas as pd
from scipy.optimize import brentq
import pulse_response as pr

HERE = os.path.dirname(os.path.abspath(__file__))
Y = {25: 365.31, 50: 173.06, 75: 101.89, 100: 74.55, 200: 58.04, 300: 63.23}   # chamber data, ns
Ts = np.array([0, 60, 114, 170, 250, 350, 500, 700, 900])


def run(name, ser_tref=125.0, **kw):
    if ser_tref != pr.T_REF:
        old = pr.T_REF
        pr.T_REF = ser_tref
        f = pr.single_electron_response
        pr.single_electron_response = lambda method="param", **k: f(t_ref=ser_tref, method=method)
        mc = pr.PulseMC(n_events=12000, **kw)
        pr.single_electron_response = f
        pr.T_REF = old
    else:
        mc = pr.PulseMC(n_events=12000, **kw)
    t0 = mc.reference_t0()

    def M(T):
        s = mc._simulate(T)
        sel = (s["amp"] >= mc.cut) & ~s["trunc"]
        return np.median(s["tom"][sel]) - t0
    Mt = np.array([M(T) for T in Ts])
    row = dict(variant=name, ser_peak_ns=mc.ser_peak)
    for T, m in zip(Ts[1:], Mt[1:]):
        row[f"c_T{T}"] = m / T
    for dv, y in Y.items():
        row[f"T_impl_dV{dv}"] = brentq(lambda T: np.interp(T, Ts, Mt) - y, 1, 899) if y < Mt[-1] else np.nan
    print(row)
    return row


if __name__ == "__main__":
    rows = [run("baseline"),
            run("SER from T_ref=100", ser_tref=100.0), run("SER from T_ref=150", ser_tref=150.0),
            run("clusters 2.0/mm", clusters_per_mm=2.0), run("clusters 3.6/mm", clusters_per_mm=3.6),
            run("Polya theta 0.5", polya_theta=0.5), run("Polya theta 3", polya_theta=3.0),
            run("diffusion a=0.03", diff_a=0.03), run("diffusion a=0.12", diff_a=0.12),
            run("amp cut 150", amp_cut=150.0), run("amp cut 600", amp_cut=600.0),
            run("exact deconvolution SER", ser="exact")]
    pd.DataFrame(rows).to_csv(os.path.join(HERE, "results", "response_systematics.csv"), index=False)
