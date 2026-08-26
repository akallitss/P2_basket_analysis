#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""pillar_numbers.py -- run the analysis once and cache every number.

The lattice scan is a few seconds per station per observable, and the figures
and the report must not be allowed to compute it separately and drift apart.
So it happens here, and both of them read the cache.

    python3 pillar_numbers.py --run eff_nominal_1
"""
import argparse
import json
import os

import numpy as np

import pillar_geom as G
import pillar_stats as P

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")


def jsonable(o):
    if isinstance(o, dict):
        return {k: jsonable(v) for k, v in o.items() if k not in ("power",)}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return jsonable(o.tolist())
    if isinstance(o, complex):
        return [o.real, o.imag]
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    return o


def build(run):
    z, j = P.load(run)
    geom = G.load()
    out = {"run": run, "n_subrun": len(j),
           "n_track": sum(b["n_tracks"] for b in j),
           "stations": {}}
    arrays = {}

    c = P.sample(z)
    masks = {s: P.build_mask(z, s) for s in P.STATIONS}
    keep = P.stack_masks(c, z, geom, masks)

    kg = np.arange(-4.4, 4.401, 0.04)
    for s in P.STATIONS:
        f = P.fine(z, s)
        blk = {"z": P.Z[s]}
        # One lattice per station, fitted on the amplitude map.  Amplitude is
        # the stronger and steadier channel of the two (17-18 sigma on all
        # three against 15-19 for efficiency, and it is linear in the charge
        # the pillar removes rather than folded through a threshold), so it
        # gives the cleanest phase -- and there is only one pillar lattice, so
        # letting the two channels fit their own would just be two estimates of
        # the same thing, with the weaker one occasionally landing a cell away.
        lat = P.hex_fit(f, "amp")
        blk["lattice"] = dict(
            d=lat["scan"]["d"], theta_deg=lat["scan"]["theta_deg"],
            cell_area=lat["cell_area"], x0=lat["x0"],
            field_min_ratio=lat["field_min_ratio"], fitted_on="amp")
        A, _ = P.modulation(f, kg, kg, "eff")
        arrays[f"kplane_{s}"] = np.abs(A).astype(np.float32)
        for what in ("eff", "amp"):
            dv, av = P.d_profile(f, lat["scan"]["theta_deg"], what)
            arrays["dscan_d"] = dv
            arrays[f"dscan_{what}_{s}"] = av
        for what in ("eff", "amp"):
            null, _, _ = P.null_scale(f, what, nk=10)
            sA = float(np.sqrt(null / 2))
            g1, a1 = P.hex_at(f, lat["scan"]["d"], lat["scan"]["theta_deg"],
                              what, order=1)
            g2, a2 = P.hex_at(f, lat["scan"]["d"], lat["scan"]["theta_deg"],
                              what, order=2)
            A1, A2 = float(np.abs(a1).mean()), float(np.abs(a2).mean())
            G1 = float(np.hypot(*g1.T).mean())
            G2 = float(np.hypot(*g2.T).mean())
            sol = P.solve_lattice(A1, G1, A2, G2, lat["cell_area"])
            dfc = P.lattice_deficit(f, lat, what)
            blk[what] = dict(
                d=lat["scan"]["d"], theta_deg=lat["scan"]["theta_deg"],
                cell_area=lat["cell_area"],
                field_min_ratio=lat["field_min_ratio"],
                A1=A1, A2=A2, abs_A=np.abs(a1), abs_A2=np.abs(a2),
                G1=G1, G2=G2, sigma_A=sA, sig1=A1 / sA, sig2=A2 / sA,
                x0=lat["x0"], solve=sol, deficit=dfc)
            # the fold bin must not be finer than the map bin (0.20 mm) or
            # the cell picture aliases into stripes
            e, val, N, _ = P.fold_hex(f, lat, what, nb=18, reach=1.8)
            r, v, n = P.fold_profile(f, lat, what, rmax=1.73, nb=30)
            arrays[f"fold_{what}_{s}"] = val
            arrays[f"foldn_{what}_{s}"] = N
            arrays[f"fold_edges_{what}_{s}"] = e
            arrays[f"prof_r_{what}_{s}"] = r
            arrays[f"prof_v_{what}_{s}"] = v
            arrays[f"prof_n_{what}_{s}"] = n
            if what == "eff":
                arrays["kplane_axis"] = kg
                arrays[f"kpower_{s}"] = lat["scan"]["power"]
                arrays[f"kd_{s}"] = lat["scan"]["d_grid"]
                arrays[f"kt_{s}"] = lat["scan"]["theta_grid"]

        # cross-check: this station's lattice read on the other two maps
        blk["cross"] = {}
        for s2 in P.STATIONS:
            if s2 == s:
                continue
            _, A = P.hex_at(P.fine(z, s2), blk["lattice"]["d"],
                            blk["lattice"]["theta_deg"], "amp")
            blk["cross"][s2] = float(np.abs(A).mean())

        # the independent resolution: the edge of the pad-residual box
        if "edge_edges" in z.files:
            blk["edge"] = {}
            for ax in ("rw", "rh"):
                fit = P.edge_fit(z, s, ax)
                if fit:
                    blk["edge"][ax] = {k: fit[k] for k in
                                       ("width", "sigma", "chi2_ndf", "n")}
                    arrays[f"edge_c_{ax}_{s}"] = fit["centres"]
                    arrays[f"edge_h_{ax}_{s}"] = fit["hist"]
                    arrays[f"edge_f_{ax}_{s}"] = fit["curve"]

        # dead areas
        blk["mask"] = dict(
            n_dead=len(masks[s]["dead"]), n_weak=len(masks[s]["weak"]),
            weak_blocks=masks[s]["weak_blocks"],
            median_block_eff=masks[s]["median_block_eff"],
            eff_all=P.station_eff(z, s), eff_masked=P.station_eff(
                z, s, masks[s]["drop"]),
            keep_frac=float(keep[s].mean()))
        out["stations"][s] = blk

    # the one big pillar inside the window, imaged
    out["big_pillar"] = {}
    for s in P.STATIONS:
        f = P.fine(z, s)
        for pxy in geom["big_xy"]:
            if not (f["xc"][0] + 8 < pxy[0] < f["xc"][-1] - 8
                    and f["yc"][0] + 8 < pxy[1] < f["yc"][-1] - 8):
                continue
            X, Y = np.meshgrid(f["xc"], f["yc"], indexing="ij")
            r = np.hypot(X - pxy[0], Y - pxy[1])
            e = np.arange(0.0, 10.01, 0.2)
            i = np.digitize(r.ravel(), e) - 1
            ok = (i >= 0) & (i < len(e) - 1)
            n = np.bincount(i[ok], f["n"].ravel()[ok], len(e) - 1)
            k = np.bincount(i[ok], f["k"].ravel()[ok], len(e) - 1)
            arrays[f"big_r_{s}"] = 0.5 * (e[1:] + e[:-1])
            arrays[f"big_n_{s}"] = n
            arrays[f"big_k_{s}"] = k
            sx = (f["xc"] > pxy[0] - 12) & (f["xc"] < pxy[0] + 12)
            sy = (f["yc"] > pxy[1] - 12) & (f["yc"] < pxy[1] + 12)
            arrays[f"bigmap_n_{s}"] = f["n"][np.ix_(sx, sy)]
            arrays[f"bigmap_k_{s}"] = f["k"][np.ix_(sx, sy)]
            arrays[f"bigmap_x_{s}"] = f["xc"][sx] - pxy[0]
            arrays[f"bigmap_y_{s}"] = f["yc"][sy] - pxy[1]
            core = n > 0
            inner = (arrays[f"big_r_{s}"] < 2.0) & core
            outer = (arrays[f"big_r_{s}"] > 7.5) & core
            out["big_pillar"][s] = dict(
                xy=pxy.tolist(), dia=float(geom["big_dia"][0]),
                eff_core=float(k[inner].sum() / max(n[inner].sum(), 1)),
                eff_far=float(k[outer].sum() / max(n[outer].sum(), 1)),
                n_core=float(n[inner].sum()))
            break

    # the stack, with and without the dead areas
    all_true = np.ones(len(c["xf"]), bool)
    out["stack"] = {"all": P.stack_block(c, all_true),
                    "masked": P.stack_block(c, keep["all"]),
                    "keep_frac": float(keep["all"].mean())}

    with open(os.path.join(DATA, f"pillar_numbers_{run}.json"), "w") as fh:
        json.dump(jsonable(out), fh, indent=1)
    np.savez_compressed(os.path.join(DATA, f"pillar_arrays_{run}.npz"),
                        **arrays)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default=P.RUN)
    a = ap.parse_args()
    o = build(a.run)
    for s, b in o["stations"].items():
        e, m = b["eff"], b["amp"]
        print(f"\n{s}  (z = {b['z']:.0f} mm)")
        print(f"  lattice   d = {e['d']:.4f} mm, theta = {e['theta_deg']:+.2f} deg, "
              f"field min ratio {e['field_min_ratio']:.4f}")
        print(f"  eff       A1 = {e['A1']:.2%} ({e['sig1']:.0f} sigma)  "
              f"A2 = {e['A2']:.2%} ({e['sig2']:.0f} sigma)  "
              f"-> D_eff {2*e['solve']['a']:.3f} mm, sigma {e['solve']['sigma']:.3f} mm")
        print(f"  amp       A1 = {m['A1']:.2%} ({m['sig1']:.0f} sigma)  "
              f"-> D_eff {2*m['solve']['a']:.3f} mm, sigma {m['solve']['sigma']:.3f} mm")
        print(f"  deficit   eff {e['deficit']['deficit']:.4f} "
              f"(plateau {e['deficit']['plateau']:.4f} vs {e['deficit']['mean']:.4f}); "
              f"amp {m['deficit']['rel_deficit']:.2%}")
        if b.get("edge"):
            for ax, v in b["edge"].items():
                print(f"  box {ax}   width {v['width']:.3f} mm, "
                      f"sigma {v['sigma']:.3f} mm (chi2/ndf {v['chi2_ndf']:.2f})")
        bp = o["big_pillar"].get(s)
        if bp:
            print(f"  big pillar D{bp['dia']:.2f} at {np.round(bp['xy'],2).tolist()}: "
                  f"core eff {bp['eff_core']:.4f} vs {bp['eff_far']:.4f} far")
    st = o["stack"]
    print(f"\nstack: keep {st['keep_frac']:.4f} of tracks after masking; "
          f"3-of-3 {st['all']['frac_3of3']:.4f} -> {st['masked']['frac_3of3']:.4f}")


if __name__ == "__main__":
    main()
