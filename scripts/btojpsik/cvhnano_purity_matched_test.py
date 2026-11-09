#!/usr/bin/env python3
"""Is purity the reason 2018's chi2 is good and ours is not?

The only fair test is same WINDOW, same MODEL, and then a purity ladder on our
side. 2018's chi2/ndf 1.44 was measured on 5.20-5.40 with
(MC signal template) x (Gaussian smear) x (free mass shift) + exp*poly2. So the
same model is applied to our data on the same window, inclusively and then
restricted to progressively purer kaon-pT bands.
"""
import numpy as np
import hist
from scipy.optimize import minimize
from wremnants.postprocessing.rabbit_cvhnano_helpers import open_results

H = "/ceph/submit/data/user/p/pmlugato/mz/histograms/"
LO, HI = 5.20, 5.40          # the 2018 window
P = print


def vv(h):
    a = np.asarray(h.view())
    return a["value"].astype(float) if a.dtype.names else a.astype(float)


def fit(data, sig, centers, nseed=16):
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
        sh = np.interp(centers - p[2] * 1e-3, centers, K @ sig, left=0.0, right=0.0)
        tot = sh.sum()
        return (np.exp(p[0]) * sh / (tot if tot > 0 else 1.0)
                + np.exp(p[3] + B.T @ p[4:6]))

    def nll(p):
        if p[0] > np.log(max(data.sum(), 1)):
            return 1e12
        mu = build(p)
        if not np.all(np.isfinite(mu)) or np.any(mu <= 0):
            return 1e12
        return float(np.sum(mu - data * np.log(mu)))

    rng = np.random.RandomState(11)
    bf, bx = 1e30, None
    for k in range(nseed):
        p0 = np.array([np.log([0.2, 0.45, 0.7][k % 3] * data.sum()),
                       [0.002, 0.004, 0.006, 0.009][k % 4],
                       [0.0, -1.0, -2.0][k % 3],
                       np.log(max(data.mean(), 1)), *(rng.randn(2) * 0.05)])
        x = p0
        for _ in range(2):
            r = minimize(nll, x, method="L-BFGS-B",
                         options=dict(maxiter=4000, maxfun=40000, ftol=1e-15))
            r = minimize(nll, r.x, method="Nelder-Mead",
                         options=dict(maxiter=15000, maxfev=15000,
                                      xatol=1e-11, fatol=1e-9))
            x = r.x
        if r.fun < bf:
            bf, bx = r.fun, x
    mu = build(bx)
    c2 = float(np.sum((data - mu) ** 2 / np.maximum(mu, 1)))
    return c2, len(data) - 6, np.exp(bx[0]), 1000 * abs(bx[1]), bx[2]


P(f"MARKER Same window 5.20-5.40, same model: signal (x)Gauss, free shift, exp*poly2")
P(f"MARKER {'sample':<34s}{'N data':>11s}{'purity':>8s}{'chi2/ndf':>10s}"
  f"{'smear':>8s}{'shift':>8s}")

# ---- 2018 reference, rebinned to 5 MeV so the bin width matches ours
with open_results(H + "btojpsik_sel2018full.hdf5") as r:
    hd = r["Charmonium_2018"]["output"]["nominal_bkmm_jpsimc_mass"].get()
    hm = r["BuToJpsiK_2018"]["output"]["nominal_bkmm_jpsimc_mass"].get()
    e = np.asarray(hd.axes[0].edges)
    d18, s18 = vv(hd), vv(hm)
f = 1
c = 0.5 * (e[:-1] + e[1:])
m = (c > LO) & (c < HI)
d, s = d18[m], s18[m] / s18[m].sum()
c2, ndf, ns, sm, dm = fit(d, s, c[m])
P(f"MARKER {'2018, full selection (2 MeV bins)':<34s}{d.sum():>11,.0f}"
  f"{ns/d.sum():>8.3f}{c2/ndf:>10.2f}{sm:>8.2f}{dm:>8.2f}")

# ---- ours, inclusively and then on purity-ordered kaon-pT bands
with open_results(H + "btojpsik_cvhnano_fullstats.hdf5") as r:
    sl = {"jointCvh_mass": slice(hist.loc(LO), hist.loc(HI))}
    hd = r["Charmonium_2016PostVFP_cvhNano"]["output"]["fitbins"].get()[sl]
    hm = r["InclusiveBToJpsiX_2016PostVFP"]["output"]["fitbins"].get()[sl]
    pc = np.asarray(hd.axes["fit_bachelor_pt"].centers)
    cm = np.asarray(hd.axes["jointCvh_mass"].centers)
    D = vv(hd)[..., 7]
    M = vv(hm)

for lab, lo_pt in (("all kaon pT", 1.0), ("kaon pT > 1.8", 1.8),
                   ("kaon pT > 3.2", 3.2), ("kaon pT > 5.0", 5.0)):
    k = pc >= lo_pt
    d = D[k].sum(axis=(0, 1, 2))
    s = M[k][..., 0].sum(axis=(0, 1, 2))
    if s.sum() < 100:
        continue
    c2, ndf, ns, sm, dm = fit(d, s / s.sum(), cm)
    P(f"MARKER {'2016 CVH, ' + lab:<34s}{d.sum():>11,.0f}"
      f"{ns/d.sum():>8.3f}{c2/ndf:>10.2f}{sm:>8.2f}{dm:>8.2f}")
