#!/usr/bin/env python3
"""
QA of the P2 photon-background simulation histograms (BkgHistograms.root).

Content of the file (no axis titles / units are stored, inferred here):
  hEnergy_<i>           photon spectrum reaching MM layer i, 200 x 5 keV, 0-1 MeV
                        (energy unit MeV: Pb K/L lines + 511 keV land where expected)
  hEnergy_response_<i>  same photons weighted by the layer response (detection prob.)
  h_rate_MM / hr_rate_MM  per-layer integrals of the two above (bin i+1 = layer i)
  h_map_xy_<i>          TH2Poly, 6 x 1280 BASKET pads (full disk), response-weighted
                        rate; scale is 1000 x hr_rate_MM
  htmpbkg (3 cycles)    leftover temporary: layer-0 spectrum in 2 keV bins

Writes figures + a qa_summary.txt next to this script (figs/).

    python3 background_sim/bkg_qa.py [--root PATH] [--out DIR]
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot
from matplotlib.collections import PolyCollection
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
DEFAULT_ROOT = ("/local/home/ak271430/Documents/PostDocSaclay/data/Simulations/"
                "background_mattieu_data/BkgHistograms.root")
BOARD_CSV = os.path.join(REPO, "Detector_Mapping/P2_BASKET/P2_BASKET_mapping.csv")

LAYERS = (0, 1, 2)
COLORS = ("#2a78d6", "#eb6834", "#1baf7a")
# reference lines used to identify the energy unit
LINES_KEV = {"Pb L": 10.55, "Pb K$\\alpha$": 74.0,
             "Pb K$\\beta$": 85.5, "e$^+$e$^-$": 511.0}

plt.rcParams.update({"font.size": 12, "axes.grid": True, "grid.alpha": 0.25})


def th1(f, name):
    h = f[name]
    v_all = h.values(flow=True)
    e_all = h.errors(flow=True)
    return dict(edges=h.axis().edges(), v=v_all[1:-1], e=e_all[1:-1],
                under=v_all[0], over=v_all[-1], entries=h.member("fEntries"))


def th2poly(f, name):
    p = f[name]
    bins = p.member("fBins")
    polys = [np.column_stack([b.member("fPoly").member("fX"),
                              b.member("fPoly").member("fY")]) for b in bins]
    c = np.array([b.member("fContent") for b in bins])
    # fArea is stored as 0 (ROOT fills it lazily) -> shoelace area of each polygon
    a = np.array([0.5 * abs(np.dot(q[:, 0], np.roll(q[:, 1], 1)) - np.dot(q[:, 1], np.roll(q[:, 0], 1)))
                  for q in polys])
    sumw2 = np.asarray(p.member("fSumw2"))[9:]          # first 9 = TH2Poly overflow regions
    over = np.asarray(p.member("fOverflow"))
    return dict(polys=polys, c=c, area=a, sumw2=sumw2, over=over,
                entries=p.member("fEntries"), title=p.member("fTitle"))


def neff(w, w2):
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(w2 > 0, w * w / w2, 0.0)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=DEFAULT_ROOT)
    ap.add_argument("--out", default=os.path.join(HERE, "figs"))
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    f = uproot.open(args.root)
    lines = [f"file: {args.root}", f"keys: {', '.join(f.keys())}", ""]

    E = {i: th1(f, f"hEnergy_{i}") for i in LAYERS}
    R = {i: th1(f, f"hEnergy_response_{i}") for i in LAYERS}
    M = {i: th2poly(f, f"h_map_xy_{i}") for i in LAYERS}
    rate, rrate = th1(f, "h_rate_MM"), th1(f, "hr_rate_MM")

    # ---------------- consistency checks ----------------
    lines.append("== consistency ==")
    for i in LAYERS:
        e, r, m = E[i], R[i], M[i]
        w_photon = (e["v"].sum() + e["over"]) / e["entries"]
        n_resp = neff(r["v"].sum() + r["over"], (r["e"] ** 2).sum())
        map_tot = m["c"].sum() + m["over"].sum()
        lines += [
            f"MM{i}: MC photons {e['entries']:.0f} (Neff of energy hist "
            f"{neff(e['v'].sum(), (e['e']**2).sum()):.0f} -> uniform weight "
            f"{w_photon:.5f}/photon), overflow >1 MeV {e['over']/w_photon:.0f} photons",
            f"      rate  sum(hEnergy)={e['v'].sum():.1f}  h_rate_MM={rate['v'][i]:.1f}"
            f" +- {rate['e'][i]:.1f}  -> {'OK' if np.isclose(e['v'].sum(), rate['v'][i], rtol=1e-6) else 'MISMATCH'}",
            f"      resp  sum(hEnergy_response)={r['v'].sum():.3f}  hr_rate_MM={rrate['v'][i]:.3f}"
            f" +- {rrate['e'][i]:.3f}  -> {'OK' if np.isclose(r['v'].sum(), rrate['v'][i], rtol=1e-6) else 'MISMATCH'}",
            f"      mean detection prob. hr/h = {rrate['v'][i]/rate['v'][i]*100:.3f} %,"
            f"  response Neff {n_resp:.0f} of {r['entries']:.0f}",
            f"      map: {m['entries']:.0f} entries, sum(pads)+overflow = {map_tot:.1f}"
            f" = {map_tot/(r['v'].sum()+r['over']):.2f} x response total"
            f" (outside pads: {m['over'].sum():.1f} = {m['over'].sum()/map_tot*100:.2f} %)",
        ]
    t = [th1(f, f"htmpbkg;{c}") for c in (1, 2, 3)]
    same = all(np.array_equal(t[0]["v"], x["v"]) for x in t[1:])
    lines += ["", f"htmpbkg: 3 cycles identical={same}; total "
              f"{t[0]['v'].sum()+t[0]['over']:.1f} vs hEnergy_0 total "
              f"{E[0]['v'].sum()+E[0]['over']:.1f} (= layer-0 spectrum, 2 keV bins, 0-400 keV)"]

    # ---------------- geometry check vs our pad table ----------------
    board = pd.read_csv(BOARD_CSV)
    cen = np.array([p.mean(axis=0) for p in M[0]["polys"]])
    r_c, phi_c = np.hypot(*cen.T), np.degrees(np.arctan2(cen[:, 1], cen[:, 0])) % 360
    sector = (phi_c // 60).astype(int)
    d, idx = cKDTree(board[["pad_cx", "pad_cy"]].values).query(cen[:1280])
    lines += ["", "== geometry ==",
              f"{len(cen)} pads = {len(cen)//1280} x 1280; r {r_c.min():.1f}-{r_c.max():.1f} mm",
              f"sector 0 vs P2_BASKET_mapping.csv: {len(set(idx))}/1280 unique matches,"
              f" max centroid offset {d.max():.2f} mm",
              f"pad area {M[0]['area'].min():.1f}-{M[0]['area'].max():.1f} mm2"]

    # ---------------- energy spectra ----------------
    ed = E[0]["edges"] * 1e3                       # keV
    fig, axs = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    for i, col in zip(LAYERS, COLORS):
        for ax, H in ((axs[0], E[i]), (axs[1], R[i])):
            ax.stairs(H["v"], ed, color=col, lw=1.6, label=f"MM{i}")
            ctr = 0.5 * (ed[1:] + ed[:-1]); ok = H["v"] > 0
            ax.errorbar(ctr[ok], H["v"][ok], H["e"][ok], fmt="none", ecolor=col, alpha=0.5, lw=1)
    for ax in axs:
        ax.set_yscale("log")
        for lab, x in LINES_KEV.items():
            ax.axvline(x, color="0.4", ls=":", lw=1)
    for lab, x in LINES_KEV.items():
        axs[0].text(x * 1.03, axs[0].get_ylim()[1] * 0.4, lab, fontsize=9, color="0.3", rotation=90, va="top")
    axs[0].set_ylabel("photon rate / 5 keV\n(file units)")
    axs[1].set_ylabel("detected rate / 5 keV\n(response-weighted)")
    axs[1].set_xlabel("photon energy [keV]  (file axis in MeV)")
    axs[0].legend(); axs[1].set_xlim(0, 800)
    fig.tight_layout(); fig.savefig(os.path.join(args.out, "energy_spectra.png"), dpi=150); plt.close(fig)

    # zoom on the low-energy part in log x, fine htmpbkg binning for layer 0
    fig, ax = plt.subplots(figsize=(10, 5))
    tb = t[2]
    ax.stairs(tb["v"], tb["edges"] * 1e3, color=COLORS[0], lw=1.4, label="MM0, 2 keV bins (htmpbkg)")
    for lab, x in LINES_KEV.items():
        if x < 400:
            ax.axvline(x, color="0.4", ls=":", lw=1); ax.text(x * 1.02, tb["v"].max() * 0.6, lab, fontsize=9, color="0.3", rotation=90)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(3, 400)
    ax.set_xlabel("photon energy [keV]"); ax.set_ylabel("photon rate / 2 keV (file units)")
    ax.legend(); fig.tight_layout(); fig.savefig(os.path.join(args.out, "energy_spectrum_MM0_fine.png"), dpi=150); plt.close(fig)

    # detection probability vs energy (response / incident)
    fig, ax = plt.subplots(figsize=(10, 5))
    lines += ["", "== detection probability vs energy (response/incident, all layers) =="]
    ctr = 0.5 * (ed[1:] + ed[:-1])
    for i, col in zip(LAYERS, COLORS):
        ok = E[i]["v"] > 0
        p = R[i]["v"][ok] / E[i]["v"][ok]
        ax.plot(ctr[ok], p * 100, "o", ms=4, color=col, label=f"MM{i}")
    ax.set_yscale("log"); ax.set_xscale("log")
    ax.set_xlabel("photon energy [keV]"); ax.set_ylabel("response / incident  [%]")
    ax.legend(); fig.tight_layout(); fig.savefig(os.path.join(args.out, "detection_probability.png"), dpi=150); plt.close(fig)
    for lo, hi in ((5, 20), (20, 70), (70, 90), (90, 200), (200, 1000)):
        s = (ctr > lo) & (ctr < hi)
        num = sum(R[i]["v"][s].sum() for i in LAYERS); den = sum(E[i]["v"][s].sum() for i in LAYERS)
        lines.append(f"  {lo:4d}-{hi:4d} keV: {num/den*100 if den else float('nan'):.3f} %  "
                     f"(share of incident {den/sum(E[i]['v'].sum() for i in LAYERS)*100:.1f} %,"
                     f" of detected {num/sum(R[i]['v'].sum() for i in LAYERS)*100:.1f} %)")

    # ---------------- rate summary ----------------
    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    x = np.arange(3)
    for ax, H, lab in ((axs[0], rate, "h_rate_MM (incident)"), (axs[1], rrate, "hr_rate_MM (detected)")):
        ax.bar(x, H["v"], yerr=H["e"], color=COLORS, capsize=4)
        ax.set_xticks(x, [f"MM{i}" for i in LAYERS]); ax.set_title(lab, fontsize=12); ax.set_ylabel("rate (file units)")
    fig.tight_layout(); fig.savefig(os.path.join(args.out, "rates_per_layer.png"), dpi=150); plt.close(fig)

    # ---------------- pad maps ----------------
    lines += ["", "== pad maps =="]
    fig, axs = plt.subplots(2, 3, figsize=(16, 10))
    for i in LAYERS:
        m = M[i]
        nf = neff(m["c"], m["sumw2"])
        dens = m["c"] / m["area"] * 100          # per cm2
        vmax = np.percentile(dens[dens > 0], 99)
        for ax, val, lab, cm, vm in ((axs[0, i], dens, "rate density / cm$^2$ (map units)", "viridis", vmax),
                                     (axs[1, i], nf, "MC photons per pad (Neff)", "magma", max(3, nf.max()))):
            pc = PolyCollection(m["polys"], array=np.ma.masked_equal(val, 0), cmap=cm,
                                edgecolors="none", clim=(0, vm))
            pc.cmap.set_bad("0.92")
            ax.add_collection(pc); ax.set_aspect("equal"); ax.autoscale_view(); ax.grid(False)
            ax.set_title(f"MM{i}", fontsize=12); fig.colorbar(pc, ax=ax, shrink=0.75, label=lab)
        lines.append(f"MM{i}: non-empty pads {(m['c']>0).sum()}/{len(m['c'])};"
                     f" Neff/pad mean {nf.mean():.2f}, max {nf.max():.1f};"
                     f" pads with Neff>=10: {(nf>=10).sum()}; per-pad stat err at mean Neff ~"
                     f" {100/np.sqrt(max(nf.mean(),1e-9)):.0f} %")
    for ax in axs.flat:
        ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]")
    fig.tight_layout(); fig.savefig(os.path.join(args.out, "pad_maps.png"), dpi=130); plt.close(fig)

    # radial and azimuthal profiles (what the statistics can actually support)
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.5))
    rb = np.linspace(120, 590, 11)
    lines += ["", "== radial / sector profiles (6 sectors summed) =="]
    for i, col in zip(LAYERS, COLORS):
        m = M[i]
        k = np.digitize(r_c, rb) - 1
        w = np.bincount(k, m["c"], len(rb))[:-1]; w2 = np.bincount(k, m["sumw2"], len(rb))[:-1]
        a = np.bincount(k, m["area"], len(rb))[:-1] / 100
        rc = 0.5 * (rb[1:] + rb[:-1])
        axs[0].errorbar(rc, w / a, np.sqrt(w2) / a, fmt="o-", color=col, ms=5, capsize=3, label=f"MM{i}")
        ws = np.bincount(sector, m["c"], 6); ws2 = np.bincount(sector, m["sumw2"], 6)
        axs[1].errorbar(np.arange(6) + (i - 1) * 0.1, ws, np.sqrt(ws2), fmt="o", color=col, ms=6, capsize=3, label=f"MM{i}")
        chi2 = (((ws - ws.mean()) ** 2) / ws2).sum()
        lines.append(f"MM{i}: inner/outer density ratio {w[0]/a[0]/(w[-1]/a[-1]):.2f};"
                     f" sector sums {np.round(ws, 0).tolist()}; chi2/ndf flat = {chi2:.1f}/5")
    axs[0].set_xlabel("r [mm]"); axs[0].set_ylabel("rate density / cm$^2$ (map units)"); axs[0].legend()
    axs[1].set_xlabel("60° sector"); axs[1].set_ylabel("rate per sector (map units)"); axs[1].legend()
    fig.tight_layout(); fig.savefig(os.path.join(args.out, "radial_sector_profiles.png"), dpi=150); plt.close(fig)

    txt = "\n".join(lines)
    with open(os.path.join(args.out, "qa_summary.txt"), "w") as fo:
        fo.write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
