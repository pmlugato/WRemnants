#!/usr/bin/env python3
"""Does coarse mass binning hide the peak S-shape -- and does it remove the bias?

The 2018 tensors all used 10 mass bins of 20 MeV. This scans bin width on BOTH
samples with the same model, and reports three things per width:

  chi2/ndf        -- whether the mis-modelling is still visible
  fitted shift    -- whether it is still there. The S-shape is antisymmetric
                     about the peak, which is what a MASS SCALE looks like, so
                     if the shift survives coarsening then the bias survives
                     and only the evidence for it is gone.
  N_sig           -- whether the yield moves

Model: MC signal template, free Gaussian smear, free shift, exp*poly2, Poisson.
"""
import numpy as np
import hist
from scipy.optimize import minimize
from wremnants.postprocessing.rabbit_cvhnano_helpers import open_results

H = "/ceph/submit/data/user/p/pmlugato/mz/histograms/"
P = print


def vv(h):
    a = np.asarray(h.view())
    return a["value"].astype(float) if a.dtype.names else a.astype(float)


def rebin(y, n):
    k = (len(y) // n) * n
    return y[:k].reshape(-1, n).sum(1)


def fit(data, sig, centers, nseed=10):
    lo, hi = centers[0], centers[-1]
    t = 2 * (centers - lo) / (hi - lo) - 1.0
    B = []
    for k in (1, 2):
        v = t ** float(k); v -= v.mean()
        for u in B:
            v -= (v @ u) / (u @ u) * u
        B.append(v / np.sqrt(np.mean(v * v)))
    B = np.array(B)
    dd = centers[:, None] - centers[None, :]

    def build(p):
        s = max(abs(p[1]), 1e-6)
        K = np.exp(-0.5 * (dd / s) ** 2); K /= K.sum(axis=0, keepdims=True)
        sh = K @ sig
        dm = p[2] * 1e-3
        if abs(dm) > 1e-9:
            sh = np.interp(centers - dm, centers, sh, left=0.0, right=0.0)
        tot = sh.sum()
        m = np.exp(p[0]) * sh / (tot if tot > 0 else 1.0)
        return m + np.exp(p[3] + B.T @ p[4:6])

    def nll(p):
        if p[0] > np.log(max(data.sum(), 1)):
            return 1e12
        mu = build(p)
        if not np.all(np.isfinite(mu)) or np.any(mu <= 0):
            return 1e12
        return float(np.sum(mu - data * np.log(mu)))

    rng = np.random.RandomState(5)
    bf, bx = 1e30, None
    for k in range(nseed):
        p0 = np.array([np.log([0.2, 0.4, 0.7][k % 3] * data.sum()),
                       [0.003, 0.005, 0.007][k % 3], [0.0, -1.0, 1.0][k % 3],
                       np.log(max(data.mean(), 1)), *(rng.randn(2) * 0.05)])
        x = p0
        for _ in range(2):
            r = minimize(nll, x, method="L-BFGS-B",
                         options=dict(maxiter=4000, maxfun=40000, ftol=1e-15))
            r = minimize(nll, r.x, method="Nelder-Mead",
                         options=dict(maxiter=12000, maxfev=12000,
                                      xatol=1e-11, fatol=1e-9))
            x = r.x
        if r.fun < bf:
            bf, bx = r.fun, x
    mu = build(bx)
    c2 = float(np.sum((data - mu) ** 2 / np.maximum(mu, 1)))
    return c2, len(data) - 6, np.exp(bx[0]), abs(bx[1]) * 1000, bx[2]


def scan(label, data0, sig0, edges0, factors):
    P(f"MARKER")
    P(f"MARKER {label}")
    P(f"MARKER {'bin width':>10s}{'bins':>6s}{'chi2/ndf':>10s}{'N_sig':>12s}"
      f"{'smear MeV':>11s}{'shift MeV':>11s}")
    w0 = (edges0[-1] - edges0[0]) / (len(edges0) - 1)
    for f in factors:
        d = rebin(data0, f); s = rebin(sig0, f)
        e = edges0[::f][: len(d) + 1]
        c = 0.5 * (e[:-1] + e[1:])
        if len(d) < 10:
            continue
        s = s / s.sum()
        c2, ndf, ns, sm, dm = fit(d, s, c)
        P(f"MARKER {1000*w0*f:>9.0f} {len(d):>6d}{c2/ndf:>10.2f}{ns:>12,.0f}"
          f"{sm:>11.2f}{dm:>11.3f}")


with open_results(H + "btojpsik_sel2018full.hdf5") as r:
    hd = r["Charmonium_2018"]["output"]["nominal_bkmm_jpsimc_mass"].get()
    hm = r["BuToJpsiK_2018"]["output"]["nominal_bkmm_jpsimc_mass"].get()
    e18 = np.asarray(hd.axes[0].edges)
    d18, s18 = vv(hd), vv(hm)
scan("2018 (BMM nano), mass axis 5.20-5.40", d18, s18, e18, [1, 2, 5, 10, 20])

LO, HI = 5.05, 5.45
with open_results(H + "btojpsik_cvhnano_fullstats.hdf5") as r:
    sl = {"jointCvh_mass": slice(hist.loc(LO), hist.loc(HI))}
    hd = r["Charmonium_2016PostVFP_cvhNano"]["output"]["nominal"].get()[sl]
    hm = r["InclusiveBToJpsiX_2016PostVFP"]["output"]["nominal"].get()[sl]
    eo = np.asarray(hd.axes["jointCvh_mass"].edges)
    do = vv(hd)[:, 7]; so = vv(hm)[:, 0]
scan("2016postVFP (CVH nano), mass axis 5.05-5.45", do, so, eo, [1, 2, 4, 8, 10])
