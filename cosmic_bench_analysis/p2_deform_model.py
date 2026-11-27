#!/usr/bin/env python3
"""p2_deform_model.py -- forward model of a P2 muon pulse, for the deformation study.

    pulse(t) = sum over primary electrons  G_e * SER(t - z_e / v_d)
    SER      = [ f_e * delta(t) + (1 - f_e) * box(0, T_ion) / T_ion ]  (x)  shaper(t; t_peak)

  * primary clusters: HEED (garfield_inputs/heed_clusters.csv.gz, Ar/iC4H10 95/5):
    cluster density along the track and the cluster-size distribution are taken
    from the file and re-thrown for any gap thickness and angle (Poisson number of
    clusters, uniform in z), so the gap is a free parameter;
  * drift: v_d(E) and longitudinal diffusion from the Magboltz table
    (garfield_inputs/gas_table.csv + the 50 V/cm point of gas_studies), E = dV / d;
  * amplification: every electron is multiplied at the mesh; the induced current is a
    fast electron part (fraction f_e ~ 1/ln G) plus the ions drifting back across
    the amplification gap at constant velocity -> a box of length
    T_ion = g^2 / (mu_ion * V_mesh); single-electron gain fluctuation = Polya(theta=1);
  * shaper: AFTER/DREAM-like  h(t) ~ exp(-3 t/tau) (t/tau)^3 sin(t/tau), tau from the
    peaking time;
  * sampling: 32 samples of 60 ns with a random trigger phase, white noise.

The point of the model is NOT absolute prediction; it gives the RESPONSE of each
measured observable to a change of the drift gap and of the ion transit time
(d observable / d T_drift, d observable / d T_ion), and shows which observables
are sensitive to only one of them.
"""
import os

import numpy as np
import pandas as pd

import p2_deform as pdf

_HERE = os.path.dirname(os.path.abspath(__file__))
TS = 60.0
NS = 32


# --------------------------------------------------------------------------- #
# gas
# --------------------------------------------------------------------------- #
def gas_table():
    g = pd.read_csv(os.path.join(_HERE, 'garfield_inputs', 'gas_table.csv'))[['E_Vcm', 'vz', 'dl', 'dt']]
    extra = os.path.join(_HERE, '..', 'gas_studies', 'results', 'ar_iso_95_5.csv')
    if os.path.isfile(extra):
        e = pd.read_csv(extra)
        e = e[(e.region == 'drift') & (e.E_Vcm < g.E_Vcm.min())]
        g = pd.concat([pd.DataFrame({'E_Vcm': e.E_Vcm, 'vz': e.vz_cm_ns, 'dl': e.dl, 'dt': e.dt}), g])
    g = g.sort_values('E_Vcm').reset_index(drop=True)
    return g


_GAS = None


def vdrift(E):
    """Drift velocity [mm/ns] at field E [V/cm]; linear to zero below the table."""
    global _GAS
    if _GAS is None:
        _GAS = gas_table()
    E = np.asarray(E, float)
    x = np.concatenate([[0.0], _GAS.E_Vcm.to_numpy()])
    y = np.concatenate([[0.0], _GAS.vz.to_numpy() * 10.0])      # cm/ns -> mm/ns
    return np.interp(E, x, y)


def dlong(E):
    """Longitudinal diffusion coefficient [sqrt(cm)] at E."""
    global _GAS
    if _GAS is None:
        _GAS = gas_table()
    return np.interp(E, _GAS.E_Vcm.to_numpy(), _GAS.dl.to_numpy())


def t_drift(d_mm, dV, V0=0.0):
    """Maximum drift time [ns] across a gap d [mm] at (dV - V0) volts."""
    E = np.clip((dV - V0), 0, None) / (np.asarray(d_mm, float) / 10.0)
    v = vdrift(E)
    return np.where(v > 0, np.asarray(d_mm, float) / np.where(v > 0, v, 1), np.inf)


# --------------------------------------------------------------------------- #
# primary ionisation
# --------------------------------------------------------------------------- #
class Heed:
    def __init__(self):
        h = pd.read_csv(os.path.join(_HERE, 'garfield_inputs', 'heed_clusters.csv.gz'))
        h0 = h[h.angle_deg == 0]
        self.gap_cm = 0.3
        self.lam = len(h0) / h0.track.nunique() / self.gap_cm      # clusters per cm of track
        self.ne = h0['ne'].to_numpy()

    def throw(self, n, d_mm, theta_deg, rng, ne_max=None):
        """n tracks through a gap d: returns (track index, z [mm], n_e) per cluster."""
        L = d_mm / 10.0 / np.cos(np.radians(theta_deg))
        k = rng.poisson(self.lam * L, n)
        it = np.repeat(np.arange(n), k)
        z = rng.uniform(0, d_mm, len(it))
        ne = self.ne[rng.integers(0, len(self.ne), len(it))]
        if ne_max:
            ne = np.minimum(ne, ne_max)
        return it, z, ne


# --------------------------------------------------------------------------- #
# single-electron response
# --------------------------------------------------------------------------- #
def shaper(t, t_peak):
    """AFTER/DREAM-like shaper, unit peak, peaking at t_peak [ns]."""
    # exp(-3x) x^3 sin(x) peaks at x = 1.0196 (root of tan x = x / (x - 1) ... numerically)
    x_pk = 1.01957
    x = np.clip(t, 0, None) / t_peak * x_pk
    h = np.exp(-3 * x) * x ** 3 * np.sin(x)
    return h / (np.exp(-3 * x_pk) * x_pk ** 3 * np.sin(x_pk))


def ser(T_ion, t_peak, f_e=0.1, dt=1.0, t_max=2400.0):
    """Single-electron response on a dt grid: (t, s). Area-normalised input charge 1."""
    t = np.arange(0, t_max, dt)
    h = shaper(t, t_peak)
    n_ion = max(int(round(T_ion / dt)), 1)
    cur = np.zeros_like(t)
    cur[:n_ion] = (1 - f_e) / n_ion
    cur[0] += f_e
    s = np.convolve(cur, h)[:len(t)]
    return t, s


# --------------------------------------------------------------------------- #
# pulses
# --------------------------------------------------------------------------- #
def simulate(n, d_mm=4.0, dV=300.0, T_ion=180.0, t_peak=123.0, theta_deg=0.0, V0=0.0,
             f_e=0.1, noise=0.0, adc_per_e=40.0, t_first=230.0, polya_theta=1.0,
             heed=None, rng=None, z_range=None, phase=True, v_mm_ns=None, dl_sqrtcm=None):
    """Simulate n muon pulses. Returns (WF [n, 32] ADC, ftst_ns [n], truth DataFrame).

    z_range=(z0, z1) [fractions of the gap]: keep only the charge created in that
    slice of the gap (for the two pads of a boundary-crossing inclined track).
    t_first: time of an electron created AT the mesh, relative to sample 0 [ns].
    """
    rng = rng or np.random.default_rng(1)
    heed = heed or Heed()
    E = max(dV - V0, 0.0) / (d_mm / 10.0)
    v = float(vdrift(E)) if v_mm_ns is None else float(v_mm_ns)      # v_mm_ns / dl_sqrtcm: another gas than the dry table
    it, z, ne = heed.throw(n, d_mm, theta_deg, rng)
    if z_range is not None:
        k = (z >= z_range[0] * d_mm) & (z < z_range[1] * d_mm)
        it, z, ne = it[k], z[k], ne[k]
    # one row per electron
    ie = np.repeat(it, ne)
    ze = np.repeat(z, ne)
    sig_l = (float(dlong(E)) if dl_sqrtcm is None else float(dl_sqrtcm)) * np.sqrt(ze / 10.0) * 10.0          # mm
    te = (ze + rng.normal(0, 1, len(ze)) * sig_l) / v
    ge = rng.gamma(polya_theta + 1.0, 1.0 / (polya_theta + 1.0), len(ze))
    ts, s = ser(T_ion, t_peak, f_e)
    ph = rng.uniform(0, TS, n) if phase else np.zeros(n)        # trigger phase within a sample
    tgrid = np.arange(NS)[None, :] * TS                          # sample times
    WF = np.zeros((n, NS))
    # sum SER(t_sample + phase - t_first - te) over electrons, chunked
    order = np.argsort(ie, kind='stable'); ie, te, ge = ie[order], te[order], ge[order]
    CH = 200000
    for a in range(0, len(ie), CH):
        sl = slice(a, a + CH)
        tt = tgrid + ph[ie[sl], None] - t_first - te[sl, None]
        val = np.interp(tt, ts, s, left=0.0, right=0.0) * ge[sl, None]
        np.add.at(WF, ie[sl], val)
    WF *= adc_per_e
    if noise > 0:
        WF += rng.normal(0, noise, WF.shape)
    q = np.bincount(ie, weights=ge, minlength=n)
    zq = np.bincount(ie, weights=ge * te, minlength=n) / np.where(q > 0, q, 1)
    truth = pd.DataFrame({'q': q, 't_q': zq, 'phase': ph})
    # ftst-like correction: the time of sample i relative to the trigger is i*60 + phase
    return WF, ph, truth, dict(v=v, E=E, T_d=d_mm / v if v > 0 else np.inf)


def observables(WF, ph, amp_min=300.0, amp_max=3500.0):
    """Medians / robust widths of the stage-47 features for a simulated sample."""
    F = pdf.features(WF)
    for c in ('t_pk', 't_l10', 't_l20', 't_l30', 't_l50', 't_t50', 't_t30', 't_cen'):
        F[c] = F[c] + ph
    k = (F.amp >= amp_min) & (F.amp <= amp_max)
    F = F[k]
    out = {'n': len(F), 'eff': float(k.mean()), 'amp_med': F.amp.median(),
           'amp_mpv': mpv(F.amp.to_numpy())}
    for c in ('t_pk', 't_l10', 't_l20', 't_l30', 't_l50', 't_t50', 't_cen'):
        out[f'med_{c}'] = F[c].median()
        out[f's68_{c}'] = pdf.s68(F[c])
    for c in ('rise', 'r2050', 'r1050', 'fwhm', 'fall', 'wid'):
        out[f'med_{c}'] = F[c].median()
        out[f's68_{c}'] = pdf.s68(F[c])
    return out, F


def mpv(a, lo=None, hi=None, nb=60):
    """Robust most-probable value: mode of a smoothed log-amplitude histogram."""
    a = np.asarray(a, float); a = a[np.isfinite(a) & (a > 0)]
    if len(a) < 30:
        return np.nan
    la = np.log(a)
    lo = np.percentile(la, 1) if lo is None else lo
    hi = np.percentile(la, 99) if hi is None else hi
    from scipy.stats import gaussian_kde
    k = gaussian_kde(la, bw_method=0.25)
    x = np.linspace(lo, hi, 400)
    return float(np.exp(x[np.argmax(k(x))]))


def dtrans(E):
    """Transverse diffusion coefficient [sqrt(cm)] at E."""
    global _GAS
    if _GAS is None:
        _GAS = gas_table()
    return np.interp(E, _GAS.E_Vcm.to_numpy(), _GAS.dt.to_numpy())


def _sample(ie, te, ge, n, ph, T_ion, t_peak, f_e, t_first):
    """Sampled pulses (n, 32) from electrons (event index, arrival time, gain)."""
    ts, s = ser(T_ion, t_peak, f_e)
    tgrid = np.arange(NS)[None, :] * TS
    WF = np.zeros((n, NS))
    CH = 200000
    for a in range(0, len(ie), CH):
        sl = slice(a, a + CH)
        tt = tgrid + ph[ie[sl], None] - t_first - te[sl, None]
        np.add.at(WF, ie[sl], np.interp(tt, ts, s, left=0.0, right=0.0) * ge[sl, None])
    return WF


def simulate_pair(n, slope, d_mm=4.0, dV=300.0, T_ion=180.0, t_peak=123.0, f_e=0.1, adc_per_e=40.0,
                  t_first=230.0, x_half=6.0, v_scale=1.0, noise=0.0, heed=None, rng=None, dt_scale=1.0):
    """Inclined tracks across the boundary (x = 0) between pad A (x < 0) and pad B (x > 0).

    The track is at x_m + slope * z at height z above the mesh (x_m uniform in
    +-x_half mm), so for slope > 0 pad B collects the charge created higher in the
    gap. Electrons drift straight down with transverse diffusion. v_scale scales the
    drift velocity at fixed gap (slow-drift hypothesis). Returns (WF_A, WF_B, phase).
    """
    rng = rng or np.random.default_rng(1)
    heed = heed or Heed()
    E = dV / (d_mm / 10.0)
    v = float(vdrift(E)) * v_scale
    theta = np.degrees(np.arctan(abs(slope)))
    it, z, ne = heed.throw(n, d_mm, theta, rng)
    xm = rng.uniform(-x_half, x_half, n)
    ie = np.repeat(it, ne); ze = np.repeat(z, ne)
    xe = xm[ie] + slope * ze + rng.normal(0, 1, len(ze)) * dt_scale * float(dtrans(E)) * np.sqrt(ze / 10.0) * 10.0
    te = (ze + rng.normal(0, 1, len(ze)) * float(dlong(E)) * np.sqrt(ze / 10.0) * 10.0) / v
    ge = rng.gamma(2.0, 0.5, len(ze))
    ph = rng.uniform(0, TS, n)
    inB = xe > 0
    WA = _sample(ie[~inB], te[~inB], ge[~inB], n, ph, T_ion, t_peak, f_e, t_first) * adc_per_e
    WB = _sample(ie[inB], te[inB], ge[inB], n, ph, T_ion, t_peak, f_e, t_first) * adc_per_e
    if noise > 0:
        WA += rng.normal(0, noise, WA.shape); WB += rng.normal(0, noise, WB.shape)
    return WA, WB, ph, dict(v=v, E=E, T_d=d_mm / v)
