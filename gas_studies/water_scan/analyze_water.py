#!/usr/bin/env python3
"""analyze_water.py -- which gas reproduces the det2 drift-field scan timing?

Data (44_deform_scan_regions, det2 drift scan 24-25 Nov 2026, mesh 430 V):
  y(dV)   = median time_of_max (leading pad, A >= 300 ADC) - zero-drift t0     [ns]
  s68(dV) = sigma68 of time_of_max                                              [ns]
  amp(dV) = median leading-pad amplitude (all), ioa(dV) = median integral/amplitude
for (a) the whole chamber and (b) the reg24 regions whose high-field gap estimate is
within 5 % of 4.0 mm (and (c) the same for reg8, as a check).

Model for a gas g, gap d:  E = dV/d,  T = d / v_d^g(E),  att = eta^g(E) d,
  y_model = s * M(T, att) + dt0,   s68_model = sqrt(S(T, att)^2 + sigma0^2)
with M, S from the cluster/pulse MC (pulse_response.py -> results/pulse_response_table.csv).
s = 1 is "c fixed" (response from the MC); s free is "c free" (empirical factor on top);
dt0 has a Gaussian prior of 5 ns (zero-drift t0 uses small pulses: walk).
Contaminant fractions are fitted continuously by interpolating ln v_d (and eta) linearly
between the Magboltz grid points of a family.
"""
import glob
import os
import sys

import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
MAG = os.path.join(RES, "magboltz")
SCAN = ("/local/home/ak271430/Documents/PostDocSaclay/data/Cosmic_Bench/Analysis/det2/"
        "p2_det2_alignment_mesh_drift_scan_long_run_430_11-24-26/44_deform_scan_regions/")
DT0_PRIOR = 5.0
GAP_TOL = 0.05

# ----------------------------------------------------------------------------- gases
FAMILIES = {
    "h2o": [("dry_95_5", 0.0), ("h2o_0p05", 0.05), ("h2o_0p1", 0.1), ("h2o_0p2", 0.2),
            ("h2o_0p3", 0.3), ("h2o_0p5", 0.5), ("h2o_0p75", 0.75), ("h2o_1p0", 1.0),
            ("h2o_1p5", 1.5), ("h2o_2p0", 2.0), ("h2o_3p0", 3.0)],
    "air": [("dry_95_5", 0.0), ("air_0p5", 0.5), ("air_1p0", 1.0), ("air_2p0", 2.0),
            ("air_5p0", 5.0)],
    "iso": [("dry_95_5", 5.0), ("iso_90_10", 10.0), ("iso_85_15", 15.0), ("iso_80_20", 20.0)],
}
UNIT = {"h2o": "% H2O", "air": "% air", "iso": "% iC4H10"}


def load_gases():
    gas = {}
    for f in sorted(glob.glob(os.path.join(MAG, "*.csv"))):
        k = os.path.basename(f)[:-4]
        if os.path.getsize(f) == 0:
            continue
        d = pd.read_csv(f)
        if len(d) < 5:
            continue
        gas[k] = d.sort_values("E_Vcm").reset_index(drop=True)
    return gas


def write_combined(gas):
    rows = []
    for k, d in gas.items():
        rows.append(d.assign(gas=k))
    allg = pd.concat(rows)[["gas"] + list(rows[0].columns[:-1])]
    allg.to_csv(os.path.join(RES, "magboltz_all_gases.csv"), index=False)
    return allg


class GasModel:
    """v_d(E), eta(E) for a (possibly interpolated) gas; log-log interpolation in E."""

    def __init__(self, tab):
        self.E = tab["E_Vcm"].to_numpy(float)
        self.lv = np.log(tab["vz_cm_ns"].to_numpy(float))
        self.eta = np.clip(tab["eta_cm"].to_numpy(float), 0, None)
        self.dl = tab["dl_sqrtcm"].to_numpy(float)

    @classmethod
    def mix(cls, fam_tabs, f):
        """fam_tabs: list of (frac, table); linear in f of ln v and eta."""
        fr = np.array([a for a, _ in fam_tabs])
        i = np.clip(np.searchsorted(fr, f) - 1, 0, len(fr) - 2)
        w = (f - fr[i]) / (fr[i + 1] - fr[i])
        a, b = cls(fam_tabs[i][1]), cls(fam_tabs[i + 1][1])
        m = cls.__new__(cls)
        m.E = a.E
        m.lv = (1 - w) * a.lv + w * np.interp(a.E, b.E, b.lv)
        m.eta = (1 - w) * a.eta + w * np.interp(a.E, b.E, b.eta)
        m.dl = (1 - w) * a.dl + w * np.interp(a.E, b.E, b.dl)
        return m

    def v(self, E):
        return np.exp(np.interp(np.log(E), np.log(self.E), self.lv))

    def etaf(self, E):
        return np.interp(np.log(E), np.log(self.E), self.eta)


# ----------------------------------------------------------------------------- response
class Response:
    def __init__(self, path=os.path.join(RES, "pulse_response_table.csv")):
        t = pd.read_csv(path)
        self.T = np.unique(t.T_ns)
        self.A = np.unique(t.att)
        g = lambda c: t.pivot(index="att", columns="T_ns", values=c).loc[self.A, self.T].to_numpy().T
        kw = dict(bounds_error=False, fill_value=None)
        self.M = RegularGridInterpolator((self.T, self.A), g("M_ns"), **kw)
        self.S = RegularGridInterpolator((self.T, self.A), g("S68_ns"), **kw)
        self.amp = RegularGridInterpolator((self.T, self.A), g("amp_med"), **kw)
        self.ioa = RegularGridInterpolator((self.T, self.A), g("ioa_med"), **kw)
        self.sel = RegularGridInterpolator((self.T, self.A), g("frac_sel"), **kw)

    def __call__(self, what, T, att):
        T = np.clip(T, self.T[0], self.T[-1])
        att = np.clip(att, self.A[0], self.A[-1])
        return getattr(self, what)(np.column_stack([T, att]))


# ----------------------------------------------------------------------------- data
def load_data():
    r = pd.read_csv(SCAN + "region_stats.csv")
    z = pd.read_csv(SCAN + "zero_drift_t0.csv")
    gap = pd.read_csv(SCAN + "gap_from_lambda_highfield.csv")
    dr = r[(r.kind == "drift") & (r.dV > 10)].copy()
    z0 = z[z.amin == 0]
    t0ch, t0ch_e = z0.ch0.iloc[0], z0.ch0_e.iloc[0]
    out = {}
    ch = dr[dr.scheme == "all"].sort_values("dV")
    out["chamber"] = pd.DataFrame(dict(
        dV=ch.dV.values, y=ch.tom_med_300.values - t0ch, ye=ch.tom_med_e_300.values,
        s68=ch.tom_s68_300.values, s68e=ch.tom_s68_e_300.values,
        amp=ch.amp_med.values, ampe=ch.amp_med_e.values, ioa=ch.ioa_med_300.values,
        ioae=ch.ioa_med_e_300.values, n=ch.n_300.values, nall=ch.n_all.values,
        ept=ch.ev_per_trig.values))
    out["chamber"].attrs.update(t0=t0ch, t0e=t0ch_e, label="chamber (all pads)")
    for sch in ("reg24", "reg8"):
        gs = gap[gap.scheme == sch]
        good = gs[np.abs(gs.d_mm / 4.0 - 1) < GAP_TOL].region.astype(int).tolist()
        sub = dr[(dr.scheme == sch) & dr.region.isin(good)]
        t0r = z0[z0.scheme == sch].set_index("region").t0
        rows = []
        for dv, g in sub.groupby("dV"):
            y = g.tom_med_300.values - t0r.loc[g.region].values
            rows.append(dict(dV=dv, y=y.mean(), ye=np.sqrt((g.tom_med_e_300**2).sum()) / len(g),
                             s68=g.tom_s68_300.mean(), s68e=np.sqrt((g.tom_s68_e_300**2).sum()) / len(g),
                             amp=g.amp_med.mean(), ampe=np.sqrt((g.amp_med_e**2).sum()) / len(g),
                             ioa=g.ioa_med_300.mean(), ioae=np.sqrt((g.ioa_med_e_300**2).sum()) / len(g),
                             n=g.n_300.sum(), nall=g.n_all.sum(), ept=g.ev_per_trig.sum()))
        d = pd.DataFrame(rows).sort_values("dV").reset_index(drop=True)
        t0e = np.sqrt((z0[z0.scheme == sch].set_index("region").t0_e.loc[good] ** 2).sum()) / len(good)
        d.attrs.update(t0=None, t0e=t0e, label=f"{sch} regions with |d/4mm-1|<5 %: {good}")
        out[sch + "_good"] = d
    return out


# ----------------------------------------------------------------------------- fit
def model_y(gm, resp, dV, d_mm, s, dt0, att_fn=None):
    """att_fn(dV) overrides the Magboltz attachment (empirical eta*d from the amplitudes)."""
    E = dV / (d_mm / 10.0)
    T = (d_mm / 10.0) / gm.v(E)
    att = gm.etaf(E) * d_mm / 10.0 if att_fn is None else att_fn(dV)
    return s * resp("M", T, att) + dt0, T, att


def fit_one(gm, resp, data, d_free=False, s_free=False, dvmin=20.0, sys_ns=0.0, d0=4.0, att_fn=None):
    m = data.dV.values >= dvmin
    dV, y = data.dV.values[m], data.y.values[m]
    ye = np.sqrt(data.ye.values[m] ** 2 + sys_ns**2)

    def unpack(p):
        i = 0
        d = d0
        if d_free:
            d = p[i]; i += 1
        s = 1.0
        if s_free:
            s = p[i]; i += 1
        return d, s, p[i]

    def chi2(p):
        d, s, dt0 = unpack(p)
        if d <= 1.0 or s <= 0.2:
            return 1e12
        ym, _, _ = model_y(gm, resp, dV, d, s, dt0, att_fn)
        return np.sum(((y - ym) / ye) ** 2) + (dt0 / DT0_PRIOR) ** 2

    p0 = ([d0] if d_free else []) + ([1.0] if s_free else []) + [0.0]
    best = None
    starts = [p0]
    if d_free:
        starts += [[dd] + p0[1:] for dd in (6.0,)]
    for st in starts:
        r = minimize(chi2, st, method="Nelder-Mead", options=dict(xatol=1e-4, fatol=1e-4, maxiter=4000))
        if best is None or r.fun < best.fun:
            best = r
    d, s, dt0 = unpack(best.x)
    npar = len(p0) - 1           # dt0 is constrained by its prior: counts as a measurement
    ndf = int(m.sum()) - npar
    return dict(chi2=best.fun, ndf=ndf, d=d, s=s, dt0=dt0, x=best.x, chi2f=chi2)


def profile_fraction(fam, gas, resp, data, fgrid, **kw):
    tabs = [(f, gas[k]) for k, f in FAMILIES[fam] if k in gas]
    tabs.sort(key=lambda a: a[0])
    fr = np.array([a for a, _ in tabs])
    fgrid = fgrid[(fgrid >= fr[0]) & (fgrid <= fr[-1])]
    prof = []
    for f in fgrid:
        gm = GasModel.mix(tabs, f)
        prof.append((f, fit_one(gm, resp, data, **kw)))
    c2 = np.array([p[1]["chi2"] for p in prof])
    ib = int(np.argmin(c2))
    fb = fgrid[ib]
    # Delta chi2 = 1 interval, and scaled by chi2/ndf (PDG) when > 1
    ndf = prof[ib][1]["ndf"] - 1
    sc = max(1.0, c2[ib] / max(ndf, 1))
    lo = fgrid[:ib + 1][c2[:ib + 1] <= c2[ib] + sc]
    hi = fgrid[ib:][c2[ib:] <= c2[ib] + sc]
    return dict(f=fb, f_lo=lo.min(), f_hi=hi.max(), chi2=c2[ib], ndf=ndf, fit=prof[ib][1],
                prof=(fgrid, c2), gm=GasModel.mix(tabs, fb), scale=sc)


if __name__ == "__main__":
    print("run make_results.py")
