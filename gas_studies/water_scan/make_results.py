#!/usr/bin/env python3
"""make_results.py -- fits, tables and figures of the det2 gas (water) study.
Run after the Magboltz CSVs are in results/magboltz/ and pulse_response.py has written
results/pulse_response_table.csv.   Output: results/*.csv, results/*.png, results/summary.txt
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze_water import (FAMILIES, RES, UNIT, GasModel, Response, fit_one, load_data,
                           load_gases, model_y, profile_fraction, write_combined)

COL = {"dry": "#257de0", "h2o": "#f47230", "air": "#00b894", "iso": "#f2a600",
       "combo": "#ee7fac", "data": "#222222"}
LOG = []


def log(s=""):
    print(s)
    LOG.append(s)


def variants():
    return [("d=4 fixed, c from MC (s=1)", dict(d_free=False, s_free=False)),
            ("d=4 fixed, s free", dict(d_free=False, s_free=True)),
            ("d free, s=1", dict(d_free=True, s_free=False)),
            ("d free, s free", dict(d_free=True, s_free=True))]


def main():
    gas = load_gases()
    allg = write_combined(gas)
    log(f"Magboltz gases loaded: {len(gas)}: {', '.join(gas)}")
    resp = Response()
    data = load_data()
    for k, d in data.items():
        log(f"dataset {k}: {d.attrs['label']}  ({len(d)} drift points)")

    # ------------------------------------------------------------------ per-gas fits
    rows = []
    for dk, d in data.items():
        for gk, tab in gas.items():
            gm = GasModel(tab)
            for vname, kw in variants():
                for dvmin in (20.0, 200.0):
                    f = fit_one(gm, resp, d, dvmin=dvmin, **kw)
                    rows.append(dict(dataset=dk, gas=gk, variant=vname, dVmin=dvmin,
                                     chi2=f["chi2"], ndf=f["ndf"], chi2ndf=f["chi2"] / f["ndf"],
                                     d_mm=f["d"], s=f["s"], dt0=f["dt0"]))
    pg = pd.DataFrame(rows)
    pg.to_csv(os.path.join(RES, "fits_per_gas.csv"), index=False)

    # ------------------------------------------------------------------ family profiles
    fgrids = {"h2o": np.round(np.arange(0, 3.0001, 0.02), 3),
              "air": np.round(np.arange(0, 5.0001, 0.05), 3),
              "iso": np.round(np.arange(5, 20.0001, 0.25), 3)}
    fam_rows, best = [], {}
    fams = [f for f in FAMILIES if sum(k in gas for k, _ in FAMILIES[f]) >= 2]
    for dk, d in data.items():
        for fam in fams:
            for vname, kw in variants():
                for sys_ns in ((0.0, 2.0) if dk == "chamber" else (0.0,)):
                    p = profile_fraction(fam, gas, resp, d, fgrids[fam], sys_ns=sys_ns, **kw)
                    fr = p["fit"]
                    fam_rows.append(dict(dataset=dk, family=fam, variant=vname, sys_ns=sys_ns,
                                         f_best=p["f"], f_lo=p["f_lo"], f_hi=p["f_hi"],
                                         chi2=p["chi2"], ndf=p["ndf"], chi2ndf=p["chi2"] / p["ndf"],
                                         d_mm=fr["d"], s=fr["s"], dt0=fr["dt0"]))
                    best[(dk, fam, vname, sys_ns)] = p
    fp = pd.DataFrame(fam_rows)
    fp.to_csv(os.path.join(RES, "fits_family_profiles.csv"), index=False)

    log("\n=== dry gas, chamber ===")
    sub = pg[(pg.gas == "dry_95_5")]
    for _, r in sub.iterrows():
        log(f"  {r.dataset:12s} {r.variant:28s} dV>={r.dVmin:4.0f}: chi2 {r.chi2:8.1f}/{r.ndf:2d}"
            f"  d {r.d_mm:5.2f}  s {r.s:5.3f}  dt0 {r.dt0:+5.1f}")
    log("\n=== family profiles (all dV >= 25 V; fraction error: Delta chi2 = max(1, chi2/ndf)) ===")
    for _, r in fp.iterrows():
        log(f"  {r.dataset:12s} {r.family:4s} {r.variant:28s} sys{r.sys_ns:3.0f}: "
            f"f = {r.f_best:6.3f} [{r.f_lo:6.3f},{r.f_hi:6.3f}] {UNIT[r.family]:9s} chi2 {r.chi2:8.1f}/{r.ndf:2d}"
            f"  d {r.d_mm:5.2f}  s {r.s:5.3f}  dt0 {r.dt0:+5.1f}")

    # combos (discrete)
    combo = pg[pg.gas.str.startswith("h2oair") & (pg.dVmin == 20.0)]
    if len(combo):
        log("\n=== water + air combinations (discrete Magboltz points) ===")
        for _, r in combo.sort_values(["dataset", "variant", "chi2"]).iterrows():
            log(f"  {r.dataset:12s} {r.gas:18s} {r.variant:28s} chi2 {r.chi2:8.1f}/{r.ndf:2d} d {r.d_mm:5.2f} s {r.s:5.3f}")

    # ------------------------------------------------------------------ figures
    ch = data["chamber"]
    key = ("chamber", "h2o", "d=4 fixed, c from MC (s=1)", 0.0)
    pw = best[key]
    fig_data_vs_models(data, gas, resp, best, pg)
    fig_profiles(best)
    fig_vd(gas, resp, data, best)
    amp_check(data, gas, resp, best)
    attachment_required(data, resp)
    water_with_empirical_attachment(data, gas, resp, best)
    loss_mechanism(data, resp)
    conversion(gas, best)
    with open(os.path.join(RES, "summary.txt"), "w") as f:
        f.write("\n".join(LOG) + "\n")


def implied_T(resp, y, att=0.0):
    """invert M(T, att) for the full drift time."""
    out = []
    for yy in np.atleast_1d(y):
        g = lambda T: resp("M", np.array([T]), np.array([att]))[0] - yy
        try:
            out.append(brentq(g, 1.0, 1999.0))
        except ValueError:
            out.append(np.nan)
    return np.array(out)


def fig_data_vs_models(data, gas, resp, best, pg):
    fig, axs = plt.subplots(2, 3, figsize=(16, 8.5), sharex="col",
                            gridspec_kw=dict(height_ratios=[3, 1.3]))
    for j, (dk, vname) in enumerate([("chamber", "d=4 fixed, c from MC (s=1)"),
                                     ("chamber", "d free, s free"),
                                     ("reg24_good", "d=4 fixed, c from MC (s=1)")]):
        d = data[dk]
        ax, axr = axs[0, j], axs[1, j]
        E = d.dV / 0.4
        ax.errorbar(E, d.y, d.ye, fmt="o", ms=4, color=COL["data"], label="data: median t_max - t0", zorder=5)
        dvf = np.linspace(20, 530, 300)
        kw = dict(d_free="d free" in vname, s_free="s free" in vname)
        curves = [("dry Ar/iC4H10 95/5", "dry", GasModel(gas["dry_95_5"]),
                   fit_one(GasModel(gas["dry_95_5"]), resp, d, **kw))]
        for fam in ("h2o", "air", "iso"):
            if (dk, fam, vname, 0.0) not in best:
                continue
            p = best[(dk, fam, vname, 0.0)]
            curves.append((f"{fam}: {p['f']:.2f} {UNIT[fam]}", fam, p["gm"], p["fit"]))
        for lab, c, gm, f in curves:
            ym, _, _ = model_y(gm, resp, dvf, f["d"], f["s"], f["dt0"])
            ax.plot(dvf / 0.4, ym, color=COL[c], lw=2,
                    label=f"{lab}  chi2 {f['chi2']:.0f}/{f['ndf']}" + (f", d={f['d']:.2f}" if kw["d_free"] else "")
                    + (f", s={f['s']:.2f}" if kw["s_free"] else ""))
            ymd, _, _ = model_y(gm, resp, d.dV.values, f["d"], f["s"], f["dt0"])
            axr.plot(E, (d.y - ymd) / d.ye, "o-", ms=3, lw=1, color=COL[c])
        ax.set_yscale("log")
        ax.set_ylabel("median time_of_max - t0 [ns]")
        ax.set_title(f"{d.attrs['label'][:40]}\n{vname}", fontsize=9)
        ax.legend(fontsize=7, frameon=False)
        ax.grid(alpha=0.25)
        axr.axhline(0, color="#999999", lw=1)
        axr.set_ylabel("pull")
        axr.set_ylim(-15, 15)
        axr.set_xlabel("nominal E = dV / 4 mm [V/cm]")
        axr.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "fit_data_vs_models.png"), dpi=130)
    plt.close(fig)


def fig_profiles(best):
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.2))
    for ax, fam in zip(axs, ("h2o", "air", "iso")):
        if ("chamber", fam, "d free, s free", 0.0) not in best:
            continue
        for (dk, ls) in (("chamber", "-"), ("reg24_good", "--"), ("reg8_good", ":")):
            for vname, c in (("d=4 fixed, c from MC (s=1)", "#257de0"), ("d=4 fixed, s free", "#f47230"),
                             ("d free, s=1", "#00b894"), ("d free, s free", "#ee7fac")):
                fg, c2 = best[(dk, fam, vname, 0.0)]["prof"]
                ax.plot(fg, c2 / best[(dk, fam, vname, 0.0)]["ndf"], ls, color=c, lw=1.5,
                        label=f"{dk}, {vname}" if fam == "h2o" else None)
        ax.set_yscale("log")
        ax.set_xlabel(UNIT[fam])
        ax.set_ylabel("chi2 / ndf")
        ax.grid(alpha=0.25)
    axs[0].legend(fontsize=6, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "chi2_profiles.png"), dpi=130)
    plt.close(fig)


def fig_vd(gas, resp, data, best):
    fig, ax = plt.subplots(figsize=(8, 5.5))
    E = np.geomspace(20, 1500, 300)
    order = ["dry_95_5", "h2o_0p1", "h2o_0p3", "h2o_0p5", "h2o_0p75", "h2o_1p0", "h2o_2p0",
             "air_1p0", "air_5p0", "iso_80_20"]
    sty = {"dry": ("-", COL["dry"]), "h2o": ("-", COL["h2o"]), "air": ("--", COL["air"]),
           "iso": (":", COL["iso"]), "h2oair": ("-.", COL["combo"])}
    for k in list(order) + [k for k in gas if k.startswith("h2oair")]:
        if k not in gas:
            continue
        fam = k.split("_")[0]
        ls, c = sty.get(fam, sty["dry"])
        g = GasModel(gas[k])
        lw = 2.5 if k == "dry_95_5" else 1.2
        a = 1.0 if fam in ("dry", "air", "iso", "h2oair") else 0.35 + 0.65 * float(k.split("_")[1].replace("p", ".")) / 2.0
        ax.plot(E, g.v(E) * 1e3, ls, color=c, lw=lw, alpha=min(a, 1.0), label=k)
    for dk, mk in (("chamber", "o"), ("reg24_good", "s")):
        d = data[dk]
        f = best[(dk, "h2o", "d=4 fixed, c from MC (s=1)", 0.0)]["fit"]
        T = implied_T(resp, d.y.values - f["dt0"])
        Te = np.abs(implied_T(resp, d.y.values - f["dt0"] + d.ye.values) - T)
        ax.errorbar(d.dV / 0.4, 0.4 / T * 1e3, 0.4 / T**2 * Te * 1e3, fmt=mk, ms=5,
                    mfc="none" if mk == "s" else COL["data"], color=COL["data"],
                    label=f"data ({dk}), T from MC response, d = 4 mm")
    ax.set_xscale("log")
    ax.set_xlabel("E [V/cm]")
    ax.set_ylabel("v_d [cm/us]")
    ax.set_title("Magboltz v_d (293 K, 760 Torr) vs data-implied v_d = 4 mm / T", fontsize=10)
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=7, ncol=2, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "vd_models_vs_data.png"), dpi=130)
    plt.close(fig)
    rows = []
    for dk, d in data.items():
        f = best[(dk, "h2o", "d=4 fixed, c from MC (s=1)", 0.0)]["fit"]
        T = implied_T(resp, d.y.values - f["dt0"])
        rows.append(pd.DataFrame(dict(dataset=dk, dV=d.dV, E_nom=d.dV / 0.4, y=d.y, ye=d.ye,
                                      T_implied_ns=T, v_implied_cm_us=0.4 / T * 1e3,
                                      T_dry_ns=0.4 / GasModel(gas["dry_95_5"]).v(d.dV / 0.4))))
    pd.concat(rows).to_csv(os.path.join(RES, "data_implied_drift_time.csv"), index=False)


def amp_check(data, gas, resp, best):
    """amplitude, integral/amplitude and sigma68 predicted by the best timing models."""
    d = data["chamber"]
    i750 = np.argmin(np.abs(d.dV - 300))
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.6))
    E = d.dV / 0.4
    axs[0].errorbar(E, d.amp / d.amp[i750], d.ampe / d.amp[i750], fmt="o", color=COL["data"], ms=4, label="data")
    axs[1].errorbar(E, d.ioa, d.ioae, fmt="o", color=COL["data"], ms=4, label="data")
    axs[2].errorbar(E, d.s68, d.s68e, fmt="o", color=COL["data"], ms=4, label="data")
    log("\n=== amplitude / shape check (chamber, d = 4 mm, s = 1, timing-fit models) ===")
    log("   model                     amp/amp(750): 62  125  187  250 V/cm | ioa: 62 125 187 | s68 chi2 (sigma0 free)")
    log(f"   data                               " + " ".join(f"{a:5.2f}" for a in (d.amp / d.amp[i750])[:4])
        + " | " + " ".join(f"{a:4.2f}" for a in d.ioa[:3]))
    cands = [("dry", GasModel(gas["dry_95_5"]), "dry")]
    for fam in ("h2o", "air", "iso"):
        if ("chamber", fam, "d=4 fixed, c from MC (s=1)", 0.0) not in best:
            continue
        p = best[("chamber", fam, "d=4 fixed, c from MC (s=1)", 0.0)]
        cands.append((f"{fam} {p['f']:.2f}", p["gm"], fam))
    for k in gas:
        if k.startswith("h2oair"):
            cands.append((k, GasModel(gas[k]), "combo"))
    rows = []
    for lab, gm, c in cands:
        T = 0.4 / gm.v(E.values)
        att = gm.etaf(E.values) * 0.4
        a = resp("amp", T, att)
        a750 = resp("amp", 0.4 / gm.v(np.array([750.0])), gm.etaf(np.array([750.0])) * 0.4)[0]
        ioa = resp("ioa", T, att)
        S = resp("S", T, att)
        s0s = np.linspace(0, 20, 201)
        c2 = [np.sum(((d.s68 - np.sqrt(S**2 + s0**2)) / d.s68e) ** 2) for s0 in s0s]
        s0 = s0s[int(np.argmin(c2))]
        log(f"   {lab:24s} " + " ".join(f"{x:5.2f}" for x in (a / a750)[:4]) + " | "
            + " ".join(f"{x:4.2f}" for x in ioa[:3]) + f" | {min(c2):7.0f}/{len(d) - 1} sigma0 {s0:.1f}"
            + f" | att(62 V/cm) = {att[0]:.3f}")
        rows.append(pd.DataFrame(dict(model=lab, dV=d.dV, E=E, T=T, att=att, amp_rel=a / a750,
                                      ioa=ioa, s68=np.sqrt(S**2 + s0**2))))
        style = "-" if c != "combo" else "-."
        axs[0].plot(E, a / a750, style, color=COL[c], label=lab)
        axs[1].plot(E, ioa, style, color=COL[c], label=lab)
        axs[2].plot(E, np.sqrt(S**2 + s0**2), style, color=COL[c], label=f"{lab} (sigma0 {s0:.0f})")
    pd.concat(rows).to_csv(os.path.join(RES, "shape_amplitude_check.csv"), index=False)
    for ax, yl in zip(axs, ("median amplitude / value at 750 V/cm", "median integral / amplitude [samples]",
                            "sigma68(time_of_max) [ns]")):
        ax.set_xscale("log")
        ax.set_xlabel("nominal E [V/cm]")
        ax.set_ylabel(yl)
        ax.grid(alpha=0.25, which="both")
        ax.legend(fontsize=7, frameon=False)
    axs[2].set_yscale("log")
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "amplitude_shape_check.png"), dpi=130)
    plt.close(fig)


def attachment_required(data, resp):
    """Model-light: T from the timing (MC response, dt0 = 0), then the attenuation att = eta d
    that the MC needs to reproduce the median amplitude relative to 750 V/cm; check with ioa."""
    log("\n=== attachment required by the amplitudes (T from timing, d = 4 mm, amp normalised at 750 V/cm) ===")
    rows = []
    for dk, d in data.items():
        T = implied_T(resp, d.y.values)
        i750 = int(np.argmin(np.abs(d.dV - 300)))
        a750 = resp("amp", np.array([T[i750]]), np.array([0.0]))[0]
        for i in range(len(d)):
            r = d.amp[i] / d.amp[i750]
            re = d.ampe[i] / d.amp[i750]
            f = lambda a, rr=r: resp("amp", np.array([T[i]]), np.array([a]))[0] / a750 - rr
            def solve(rr):
                try:
                    return brentq(lambda a: f(a, rr), 0.0, 3.0)
                except ValueError:
                    return 0.0 if rr > 0.5 else np.nan
            att = solve(r)
            att_lo, att_hi = solve(r + re), solve(r - re)
            ioa_m = resp("ioa", np.array([T[i]]), np.array([att if np.isfinite(att) else 0.0]))[0]
            ioa_0 = resp("ioa", np.array([T[i]]), np.array([0.0]))[0]
            rows.append(dict(dataset=dk, dV=d.dV[i], E_nom=d.dV[i] / 0.4, T_ns=T[i], amp_rel=r, amp_rel_e=re,
                             att=att, att_lo=att_lo, att_hi=att_hi, eta_cm=att / 0.4,
                             ioa_data=d.ioa[i], ioa_mc_att0=ioa_0, ioa_mc_att=ioa_m))
    a = pd.DataFrame(rows)
    a.to_csv(os.path.join(RES, "attachment_required.csv"), index=False)
    for _, r in a[(a.dataset == "chamber") & (a.dV < 260)].iterrows():
        log(f"   E {r.E_nom:6.1f}  T {r.T_ns:5.0f} ns  amp/amp750 {r.amp_rel:.3f}  eta*d {r.att:.2f} [{r.att_lo:.2f},{r.att_hi:.2f}]"
            f"  eta {r.eta_cm:.2f}/cm  ioa data {r.ioa_data:.2f}  MC(no att) {r.ioa_mc_att0:.2f}  MC(att) {r.ioa_mc_att:.2f}")


def water_with_empirical_attachment(data, gas, resp, best):
    """water fraction refitted with the attenuation eta*d(dV) required by the amplitudes
    (attachment_required.csv, same data set) inserted in the timing response."""
    a = pd.read_csv(os.path.join(RES, "attachment_required.csv"))
    fg = np.round(np.arange(0, 3.0001, 0.02), 3)
    log("\n=== water fraction with the empirical attenuation from the amplitudes in the response ===")
    rows = []
    for dk, d in data.items():
        ad = a[a.dataset == dk].sort_values("dV")
        att = np.nan_to_num(ad.att.to_numpy(), nan=0.0)
        fn = lambda dV, x=ad.dV.to_numpy(), y=att: np.interp(dV, x, y)
        for vname, kw in [("d=4 fixed, c from MC (s=1)", dict(d_free=False, s_free=False)),
                          ("d=4 fixed, s free", dict(d_free=False, s_free=True)),
                          ("d free, s=1", dict(d_free=True, s_free=False))]:
            for sys_ns in (0.0, 2.0):
                p = profile_fraction("h2o", gas, resp, d, fg, sys_ns=sys_ns, att_fn=fn, **kw)
                fr = p["fit"]
                best[(dk, "h2o+att", vname, sys_ns)] = p
                rows.append(dict(dataset=dk, variant=vname, sys_ns=sys_ns, f_best=p["f"], f_lo=p["f_lo"],
                                 f_hi=p["f_hi"], chi2=p["chi2"], ndf=p["ndf"], d_mm=fr["d"], s=fr["s"], dt0=fr["dt0"]))
                log(f"  {dk:12s} {vname:28s} sys{sys_ns:3.0f}: f = {p['f']:5.2f} [{p['f_lo']:5.2f},{p['f_hi']:5.2f}] % H2O"
                    f"  chi2 {p['chi2']:7.1f}/{p['ndf']}  d {fr['d']:5.2f}  s {fr['s']:5.3f}  dt0 {fr['dt0']:+5.1f}")
    pd.DataFrame(rows).to_csv(os.path.join(RES, "fits_water_empirical_attachment.csv"), index=False)
    # scale k of the attachment (shape = what the amplitudes require), water fraction free
    log("\n=== timing-only limit on depth-dependent attachment: att = k * att_required(dV), water free (d=4, s=1, sys 2 ns) ===")
    rows = []
    for dk in ("chamber", "reg24_good"):
        d = data[dk]
        ad = a[a.dataset == dk].sort_values("dV")
        att = np.nan_to_num(ad.att.to_numpy(), nan=0.0)
        c0 = None
        for k in (0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0):
            fn = lambda dV, x=ad.dV.to_numpy(), y=att, kk=k: kk * np.interp(dV, x, y)
            p = profile_fraction("h2o", gas, resp, d, fg, sys_ns=2.0, att_fn=fn)
            c0 = p["chi2"] if c0 is None else c0
            rows.append(dict(dataset=dk, k=k, att62=k * att[0], f_best=p["f"], chi2=p["chi2"], ndf=p["ndf"],
                             dchi2=p["chi2"] - c0, dt0=p["fit"]["dt0"]))
            log(f"  {dk:12s} k {k:4.2f} (eta*d at 62 V/cm {k * att[0]:4.2f}): f {p['f']:4.2f} %  chi2 {p['chi2']:7.1f}/{p['ndf']}"
                f"  Delta {p['chi2'] - c0:+7.1f}  dt0 {p['fit']['dt0']:+5.1f}")
    pd.DataFrame(rows).to_csv(os.path.join(RES, "attachment_scale_scan.csv"), index=False)


def loss_mechanism(data, resp):
    """Can the amplitude loss be depth-dependent attachment? Compare, at the data-implied T,
    attachment exp(-att z/d) (att from the amplitudes) with a depth-independent loss g that gives
    the same median amplitude: timing y, amplitude, ioa, selected fraction (direct MC, 20k events)."""
    import pulse_response as pr
    mc = pr.PulseMC(n_events=20000)
    t0 = mc.reference_t0()
    a = pd.read_csv(os.path.join(RES, "attachment_required.csv"))
    a = a[a.dataset == "chamber"].sort_values("dV")
    d = data["chamber"]
    a750 = np.median(mc._simulate(float(a.T_ns.iloc[int(np.argmin(np.abs(a.dV - 300)))]))["amp"])
    rows = []
    log("\n=== loss mechanism (chamber): depth-dependent attachment vs depth-independent loss, MC at the data-implied T ===")
    for i in range(4):
        r = a.iloc[i]
        base = mc._simulate(r.T_ns)
        g = r.amp_rel / (np.median(base["amp"]) / a750)
        for lab, s, gg in (("no loss", base, 1.0), ("attachment", mc._simulate(r.T_ns, att=r.att), 1.0),
                           ("uniform loss", base, g)):
            amp = s["amp"] * gg
            sel = (amp >= 300) & ~s["trunc"]
            row = dict(E_nom=r.E_nom, T_ns=r.T_ns, model=lab, att=r.att if lab == "attachment" else 0.0,
                       gain_factor=gg, y_mc=np.median(s["tom"][sel]) - t0, y_data=d.y[i], y_data_e=d.ye[i],
                       amp_rel_mc=np.median(amp) / a750, amp_rel_data=r.amp_rel,
                       ioa_mc=np.median(s["ioa"][sel]), ioa_data=d.ioa[i], sel_mc=sel.mean(), sel_data=d.n[i] / d.nall[i])
            rows.append(row)
            log(f"   E {r.E_nom:6.1f} T {r.T_ns:4.0f} {lab:13s}: y {row['y_mc']:6.1f} (data {row['y_data']:6.1f})"
                f"  amp {row['amp_rel_mc']:.3f} ({r.amp_rel:.3f})  ioa {row['ioa_mc']:.2f} ({d.ioa[i]:.2f})"
                f"  sel {row['sel_mc']:.2f} ({row['sel_data']:.2f})")
    pd.DataFrame(rows).to_csv(os.path.join(RES, "loss_mechanism_check.csv"), index=False)


def conversion(gas, best):
    p = best[("chamber", "h2o", "d=4 fixed, c from MC (s=1)", 0.0)]
    fs = [v["f"] for k, v in best.items() if k[1] == "h2o"]   # water-only fits, all data sets, variants, sys 0 and 2
    f_lo, f_hi = min(fs), max(fs)
    tabs = [(f, gas[k]) for k, f in FAMILIES["h2o"] if k in gas]
    log(f"\n   water fraction: reference {p['f']:.2f} % (chamber, d = 4, s = 1); spread over data sets and"
        f" variants (water only, d / s free, sys 0 / 2 ns): {f_lo:.2f} ... {f_hi:.2f} %")
    models = [("dry_95_5", GasModel(gas["dry_95_5"])), (f"h2o_{p['f']:.2f}pct_bestfit", p["gm"]),
              (f"h2o_lo_{f_lo:.2f}pct", GasModel.mix(tabs, f_lo)),
              (f"h2o_hi_{f_hi:.2f}pct", GasModel.mix(tabs, f_hi))]
    p = dict(p, f_lo=f_lo, f_hi=f_hi)
    # v_d tables
    E = np.array([20, 30, 40, 50, 62.5, 75, 87.5, 100, 125, 150, 187.5, 225, 250, 300, 350, 400, 450,
                  500, 562.5, 625, 700, 750, 875, 1000, 1150, 1312.5, 1500.])
    vt = pd.DataFrame({"E_Vcm": E})
    for k, gm in models:
        vt[f"vd_cm_us_{k}"] = gm.v(E) * 1e3
    vt.to_csv(os.path.join(RES, "vd_table_dry_and_best.csv"), index=False)
    # T(d) at dV = 200, 300, 400
    dd = np.round(np.arange(2.5, 9.0001, 0.25), 3)
    rows = []
    for k, gm in models:
        for dV in (200., 300., 400.):
            for d in dd:
                E1 = dV / (d / 10.)
                rows.append(dict(gas=k, dV=dV, d_mm=d, E_Vcm=E1, vd_cm_us=gm.v(np.array([E1]))[0] * 1e3,
                                 T_ns=(d / 10.) / gm.v(np.array([E1]))[0]))
    tt = pd.DataFrame(rows)
    tt.to_csv(os.path.join(RES, "T_vs_gap_dry_and_best.csv"), index=False)
    # inverse
    log("\n=== drift time T -> gap d (T = d / v_d(dV/d)) ===")
    log(f"   best water fraction (chamber, d = 4 fixed, s = 1): {p['f']:.2f} % [{p['f_lo']:.2f}, {p['f_hi']:.2f}]")
    inv = []
    for dV in (200., 300., 400.):
        for T in (107., 125., 170., 190.):
            line = dict(dV=dV, T_ns=T)
            for k, gm in models:
                g = lambda d: (d / 10.) / gm.v(np.array([dV / (d / 10.)]))[0] - T
                try:
                    line[f"d_mm_{k}"] = brentq(g, 0.5, 20.0)
                except ValueError:
                    line[f"d_mm_{k}"] = np.nan
            band = []
            for f in np.arange(p["f_lo"], p["f_hi"] + 1e-6, 0.02):
                gm = GasModel.mix(tabs, f)
                try:
                    band.append(brentq(lambda d: (d / 10.) / gm.v(np.array([dV / (d / 10.)]))[0] - T, 0.5, 20.0))
                except ValueError:
                    pass
            line["d_mm_band_min"], line["d_mm_band_max"] = min(band), max(band)
            inv.append(line)
    inv = pd.DataFrame(inv)
    kb = models[1][0]
    inv["delta_mm_best_minus_dry"] = inv[f"d_mm_{kb}"] - inv["d_mm_dry_95_5"]
    inv["delta_pct"] = 100 * inv.delta_mm_best_minus_dry / inv["d_mm_dry_95_5"]
    inv.to_csv(os.path.join(RES, "gap_from_T_dry_vs_best.csv"), index=False)
    for _, r in inv.iterrows():
        log(f"   dV {r.dV:3.0f} V  T {r.T_ns:5.0f} ns : d_dry {r['d_mm_dry_95_5']:5.3f}  d_best {r[f'd_mm_{kb}']:5.3f}"
            f"  [band {r.d_mm_band_min:5.3f}, {r.d_mm_band_max:5.3f}]  delta {r.delta_mm_best_minus_dry:+6.3f} mm ({r.delta_pct:+5.1f} %)")
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.5))
    for k, gm in models[:2]:
        c = COL["dry"] if k.startswith("dry") else COL["h2o"]
        for dV, ls in ((200., ":"), (300., "-"), (400., "--")):
            s = tt[(tt.gas == k) & (tt.dV == dV)]
            axs[0].plot(s.d_mm, s.T_ns, ls, color=c, label=f"{k}, dV = {dV:.0f} V")
        Ef = np.geomspace(20, 1500, 300)
        axs[1].plot(Ef, gm.v(Ef) * 1e3, color=c, label=k)
    axs[0].set_xlabel("gap d [mm]"); axs[0].set_ylabel("full drift time T [ns]")
    axs[1].set_xscale("log"); axs[1].set_xlabel("E [V/cm]"); axs[1].set_ylabel("v_d [cm/us]")
    for ax in axs:
        ax.grid(alpha=0.25); ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "T_to_gap_conversion.png"), dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    main()
