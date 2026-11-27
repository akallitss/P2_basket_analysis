#!/usr/bin/env python3
"""pulse_response.py -- drift-time response of the det2 leading-pad timing observables.

1. Single-electron response (SER) h(t) by exact deconvolution of the chamber mean pulse
   m = box(T_ref) (x) h, T_ref = 125 ns (stage 47 argument):
       h(t) - h(t - T_ref) = T_ref m'(t)   =>   h(t) = T_ref * sum_k m'(t - k T_ref).
2. Cluster-level MC of one leading pad: clusters Poisson along the gap (uniform z),
   arrival time t = T z/d (+ longitudinal diffusion sigma = a T sqrt(z/d)), cluster size
   from an Ar-like distribution, Polya gain per electron, optional attachment (survival
   exp(-eta*(d-z)) ... z measured from the mesh, so the path is z), waveform sampled every
   60 ns with a random phase, time_of_max = parabola through the 3 samples around the max
   (as WaveformAnalyzer). Observables per T:
       M(T)   = median(tom | A >= cut) - median(tom | T = 0)       (data: median - t0)
       S(T)   = sigma68(tom | A >= cut)                           (data: sigma68 (-) sigma0)
       A(T)   = median amplitude (all), I/A(T) = median integral/amplitude (A >= cut)
   The amplitude scale is calibrated so that median A at T_cal (dry gas, 750 V/cm) = A_cal ADC.
"""
import numpy as np

MEAN_PULSE = ("/local/home/ak271430/Documents/PostDocSaclay/data/Cosmic_Bench/Analysis/det2/"
              "p2_det2_long_run_430_repositioned_11-25-26/long_run_det2_430_730/"
              "47_deform_pulse_maps/mean_pulse_ref.csv")
T_REF = 125.0
DT = 2.5            # ns, fine grid
SAMPLE = 60.0       # ns, ADC sampling period
WINDOW_NS = 32 * SAMPLE   # readout window (p2_waveforms.N_SAMP = 32)
T0_ABS = 393.1      # ns, zero-drift median time_of_max (chamber, 44_deform_scan_regions)
TRUNC_ADC = 40.0    # ADC, "still above threshold at the last sample" (assumed ~5 sigma noise)


def _crrc(t, n, tau):
    x = np.clip(t, 0, None) / tau
    return np.where(t > 0, x**n * np.exp(-x), 0.0)


def ser_param(t, t0, n, tau, u, tau2):
    """CR-RC^n main lobe minus a slow undershoot lobe, normalised to peak 1 later."""
    return _crrc(t - t0, n, tau) / _crrc(np.array([n * tau]), n, tau)[0] \
        - u * _crrc(t - t0, n, tau2) / _crrc(np.array([n * tau2]), n, tau2)[0]


def single_electron_response(path=MEAN_PULSE, t_ref=T_REF, method="param"):
    """Return (t rel. to SER start, h, t_meanpulse_fine, m_fine, info).
    method 'exact': h = T_ref * sum_k m'(t - k T_ref) (rings in the tail);
    method 'param': CR-RC^n with undershoot, fitted through box(T_ref) (x) h to m."""
    from scipy.interpolate import CubicSpline
    from scipy.optimize import least_squares
    d = np.genfromtxt(path, delimiter=",", names=True)
    t, m = d["t_ns"], d["m_ref"]
    tf = np.arange(t[0], t[-1] + 1e-9, DT)
    mf = CubicSpline(t, m)(tf)
    k = int(round(t_ref / DT))
    if method == "exact":
        dm = np.gradient(mf, DT)
        h = np.zeros_like(mf)
        for i in range(len(h)):
            h[i] = t_ref * dm[i] + (h[i - k] if i >= k else 0.0)
        h = h / h.max()
        t0 = tf[np.argmax(h > 0.01)]
        return tf - t0, h, tf, mf, {}
    box = np.ones(k) / k

    def model(p):
        a = p[0]
        hh = ser_param(tf, *p[1:])
        return a * np.convolve(hh, box)[:len(tf)]

    p0 = [1.0, 180.0, 3.0, 65.0, 0.02, 200.0]
    r = least_squares(lambda p: model(p) - mf, p0,
                      bounds=([0.1, 0, 0.5, 10, 0, 50], [10, 400, 12, 300, 0.5, 2000]))
    a, t0, n, tau, u, tau2 = r.x
    th = np.arange(-50, 2500, DT)
    h = ser_param(th, 0.0, n, tau, u, tau2)
    h = h / h.max()
    info = dict(t0=t0, n=n, tau=tau, u=u, tau2=tau2, rms=float(np.sqrt(np.mean(r.fun**2))),
                maxerr=float(np.abs(r.fun).max()), peak_ns=float(th[np.argmax(h)]))
    return th, h, tf, mf, info


class PulseMC:
    def __init__(self, n_events=20000, clusters_per_mm=2.8, diff_a=0.06, polya_theta=1.0,
                 amp_cut=300.0, A_cal=1296.0, T_cal=114.0, seed=1, ser="param"):
        self.th, self.h, self.tm, self.m, self.ser_info = single_electron_response(method=ser)
        self.n = n_events
        self.cpm = clusters_per_mm
        self.a = diff_a
        self.theta = polya_theta
        self.cut = amp_cut
        self.rng = np.random.default_rng(seed)
        # Ar-like cluster-size distribution (Fischle et al. order of magnitude), tail 1/n^2
        n = np.arange(1, 61)
        p = np.where(n == 1, 0.656, np.where(n == 2, 0.15, np.where(n == 3, 0.064,
                     np.where(n == 4, 0.035, np.where(n == 5, 0.0225, 0.0225 * 25 / n**2)))))
        self.cs_n, self.cs_p = n, p / p.sum()
        self.ser_peak = float(self.th[np.argmax(self.h)])
        self.scale = 1.0
        self.scale = A_cal / np.median(self._simulate(T_cal)["amp"])

    def _hval(self, t):
        return np.interp(t, self.th, self.h, left=0.0, right=0.0)

    def _simulate(self, T, d_mm=4.0, att=0.0, n=None, fixed_rng=True):
        """att = eta * d (attenuation over the full gap); survival of an electron
        created at z (from the mesh) = exp(-att z/d)."""
        rng = np.random.default_rng(12345) if fixed_rng else self.rng
        n = n or self.n
        lam = self.cpm * d_mm
        ncl = np.maximum(rng.poisson(lam, n), 1)
        tot = ncl.sum()
        ev = np.repeat(np.arange(n), ncl)
        u = rng.random(tot)                       # z/d, z from the mesh
        size = rng.choice(self.cs_n, tot, p=self.cs_p)
        if att > 0:                               # binomial survival over the path z
            ps = np.exp(-att * u)
            size = rng.binomial(size, ps)
        q = rng.gamma(self.theta * np.maximum(size, 1), 1.0 / self.theta) * (size > 0)
        t = T * u + self.a * T * np.sqrt(u) * rng.standard_normal(tot)
        # sample grid: phase in [0, 60), samples from -120 ns to T + 1200 ns
        ns = int((T + 1300) // SAMPLE) + 3
        phase = rng.random(n) * SAMPLE
        ts = phase[:, None] - 120.0 + SAMPLE * np.arange(ns)[None, :]
        w = np.zeros((n, ns))
        # accumulate cluster by cluster in blocks to bound memory
        order = np.argsort(ev, kind="stable")
        ev, t, q = ev[order], t[order], q[order]
        blk = 200000
        for s in range(0, tot, blk):
            e = ev[s:s + blk]
            np.add.at(w, e, q[s:s + blk, None] * self._hval(ts[e] - t[s:s + blk, None]))
        w *= self.scale
        im = np.clip(np.argmax(w, axis=1), 1, ns - 2)
        r = np.arange(n)
        y0, y1, y2 = w[r, im - 1], w[r, im], w[r, im + 1]
        den = y0 - 2 * y1 + y2
        off = np.where(den < 0, 0.5 * (y0 - y2) / np.where(den < 0, den, -1), 0.0)
        off = np.clip(off, -0.5, 0.5)
        tom = ts[r, im] + off * SAMPLE
        amp = y1 - 0.25 * (y0 - y2) * off
        integ = np.where(w > 0, w, 0).sum(axis=1)
        # readout window: 32 x 60 ns from the trigger; a T = 0 pulse peaks at T0_ABS ns.
        # right-truncated = the waveform is still above TRUNC_ADC at the last sample.
        t_last = WINDOW_NS - SAMPLE - (T0_ABS - self.ser_peak)        # in the t(arrival) frame
        wl = np.zeros(n)
        blk2 = 200000
        for s0 in range(0, tot, blk2):
            np.add.at(wl, ev[s0:s0 + blk2], q[s0:s0 + blk2] * self._hval(t_last - t[s0:s0 + blk2]))
        wl *= self.scale
        return {"tom": tom, "amp": amp, "ioa": integ / np.maximum(amp, 1e-9),
                "trunc": wl > TRUNC_ADC}

    def reference_t0(self):
        """median tom for T = 0 (all charge at the mesh), A >= 0 (zero-drift data use all)."""
        s = self._simulate(0.0)
        return np.median(s["tom"])

    def table(self, Ts, atts=(0.0,), d_mm=4.0):
        t0 = self.reference_t0()
        rows = []
        for att in atts:
            for T in Ts:
                s = self._simulate(T, d_mm, att)
                sel = (s["amp"] >= self.cut) & ~s["trunc"]
                tm = s["tom"][sel]
                p16, p50, p84 = np.percentile(tm, [16, 50, 84])
                rows.append((T, att, p50 - t0, 0.5 * (p84 - p16), np.median(s["amp"]),
                             np.median(s["ioa"][sel]), sel.mean(), s["trunc"].mean()))
        return np.array(rows)


if __name__ == "__main__":
    import sys
    out = sys.argv[1] if len(sys.argv) > 1 else "pulse_response_table.csv"
    mc = PulseMC()
    Ts = np.concatenate([np.arange(0, 300, 20), np.arange(300, 1000, 50), np.arange(1000, 2001, 100)])
    atts = [float(a) for a in sys.argv[2].split(",")] if len(sys.argv) > 2 else \
        [0.0, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0]
    tab = mc.table(Ts, atts)
    np.savetxt(out, tab, delimiter=",",
               header="T_ns,att,M_ns,S68_ns,amp_med,ioa_med,frac_sel,frac_trunc",
               comments="", fmt="%.4f")
    print("SER", mc.ser_info)
    for r in tab[::7]:
        print(" ".join(f"{x:9.3f}" for x in r))
