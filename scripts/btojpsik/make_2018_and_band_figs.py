#!/usr/bin/env python3
"""Two figures that together answer 'why did 2018 work and we do not'.

  A  the 2018 data fitted with EXACTLY the model class we use, drawn like ours
  B  our data fitted per kaon-pT band: chi2/ndf and the extra smear it wants
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import hist
from scipy.optimize import minimize
from wremnants.postprocessing.rabbit_cvhnano_helpers import open_results

H = "/ceph/submit/data/user/p/pmlugato/mz/histograms/"
OUT = os.path.expanduser("~/public_html/slides/assets_btojpsik_overview")
INK, ACC = "#1F1F1F", "#C0392B"


def vv(h):
    a = np.asarray(h.view())
    return a["value"].astype(float) if a.dtype.names else a.astype(float)


# ------------------------------------------------------------------ A
def fig2018():
    with open_results(H + "btojpsik_sel2018full.hdf5") as r:
        hd = r["Charmonium_2018"]["output"]["nominal_bkmm_jpsimc_mass"].get()
        hm = r["BuToJpsiK_2018"]["output"]["nominal_bkmm_jpsimc_mass"].get()
        e = np.asarray(hd.axes[0].edges)
        c = 0.5 * (e[:-1] + e[1:])
        data, sig = vv(hd), vv(hm)
    sig = sig / sig.sum()
    lo, hi = e[0], e[-1]
    t = 2 * (c - lo) / (hi - lo) - 1.0
    B = []
    for k in (1, 2):
        v = t ** float(k); v -= v.mean()
        for u in B:
            v -= (v @ u) / (u @ u) * u
        B.append(v / np.sqrt(np.mean(v * v)))
    B = np.array(B)
    dd = c[:, None] - c[None, :]

    def sm(x, s):
        K = np.exp(-0.5 * (dd / max(s, 1e-6)) ** 2)
        K /= K.sum(axis=0, keepdims=True)
        return K @ x

    def build(p):
        sh = sm(sig, abs(p[1])); sh /= sh.sum()
        return np.exp(p[0]) * sh + np.exp(p[2] + B.T @ p[3:5])

    def nll(p):
        mu = build(p)
        if not np.all(np.isfinite(mu)) or np.any(mu <= 0):
            return 1e12
        return float(np.sum(mu - data * np.log(mu)))

    bf, bx = 1e30, None
    for s0 in (0.002, 0.004, 0.006):
        x = np.array([np.log(0.7 * data.sum()), s0, np.log(data.mean()), 0.0, 0.0])
        for _ in range(3):
            rr = minimize(nll, x, method="L-BFGS-B",
                          options=dict(maxiter=4000, maxfun=40000, ftol=1e-15))
            rr = minimize(nll, rr.x, method="Nelder-Mead",
                          options=dict(maxiter=20000, maxfev=20000, xatol=1e-11, fatol=1e-9))
            x = rr.x
        if rr.fun < bf:
            bf, bx = rr.fun, x
    mu = build(bx); pull = (data - mu) / np.sqrt(mu)
    c2 = float(np.sum(pull ** 2)); ndf = len(data) - 5
    print(f"MARKER 2018: chi2/ndf {c2/ndf:.2f}  N_sig {np.exp(bx[0]):,.0f} "
          f"purity {np.exp(bx[0])/data.sum():.3f}  smear {1000*abs(bx[1]):.2f} MeV")

    f, (a1, a2) = plt.subplots(2, 1, figsize=(9.4, 5.4), sharex=True,
                               gridspec_kw=dict(height_ratios=[2.3, 1], hspace=0.06))
    a1.errorbar(c, data, yerr=np.sqrt(data), fmt="k.", ms=4, lw=1, label="2018 data")
    a1.plot(c, mu, "-", color=ACC, lw=1.7, label="signal ⊗ Gauss + exp·poly2")
    a1.plot(c, np.exp(bx[2] + B.T @ bx[3:5]), "b:", lw=1.4, label="background")
    a1.set_yscale("log"); a1.set_ylim(2e2, 4e5)
    a1.set_ylabel("candidates / 2 MeV")
    a1.legend(frameon=False, fontsize=9)
    a1.text(0.02, 0.93, f"purity 0.72,   $\\chi^2$/ndf = {c2/ndf:.0f}",
            transform=a1.transAxes, fontsize=11, color=INK)
    a2.axhline(0, color="grey", lw=0.8)
    a2.plot(c, pull, "k.-", ms=3, lw=0.8)
    a2.set_ylim(-25, 25)
    a2.set_xlabel("m($J/\\psi K$)  [GeV]"); a2.set_ylabel("pull")
    a2.set_xlim(e[0], e[-1])
    p = os.path.join(OUT, "fit2018_samemodel.png")
    f.savefig(p, dpi=150, bbox_inches="tight"); plt.close(f)
    print("MARKER wrote", p)


# ------------------------------------------------------------------ B
def figbands():
    bands = [("1.0–1.4", 0.053, 7.54, 9.18), ("1.4–1.8", 0.093, 3.96, 5.42),
             ("1.8–3.2", 0.236, 3.98, 4.61), ("3.2–10", 0.611, 3.97, 4.13)]
    lab = [b[0] for b in bands]
    pur = np.array([b[1] for b in bands])
    c2 = np.array([b[2] for b in bands])
    sg = np.array([b[3] for b in bands])
    x = np.arange(len(bands))
    f, ax = plt.subplots(figsize=(9.4, 4.6))
    ax.plot(x, c2, "o-", color=ACC, ms=8, lw=2, label="$\\chi^2$/ndf in the band")
    ax.axhline(9.3, color=ACC, ls=":", lw=1.4)
    ax.text(0.55, 9.55, "9.3 — all bands fitted together with ONE shape",
            fontsize=9.5, color=ACC)
    ax.set_ylabel("$\\chi^2$/ndf", color=ACC)
    ax.tick_params(axis="y", labelcolor=ACC)
    ax.set_ylim(0, 12)
    a2 = ax.twinx()
    a2.plot(x, sg, "s--", color="#2E86C1", ms=7, lw=1.8,
            label="extra smearing the band wants")
    a2.set_ylabel("extra Gaussian smearing  [MeV]", color="#2E86C1")
    a2.tick_params(axis="y", labelcolor="#2E86C1")
    a2.set_ylim(0, 11)
    for i, p in enumerate(pur):
        ax.annotate(f"purity\n{p:.3f}", (i, 0.7), ha="center", fontsize=9, color=INK)
    ax.set_xticks(x); ax.set_xticklabels(lab)
    ax.set_xlabel("kaon $p_T$ band  [GeV]")
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = a2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, frameon=False, fontsize=9.5, loc="upper center")
    ax.grid(alpha=0.15)
    p = os.path.join(OUT, "band_chi2_and_smear.png")
    f.savefig(p, dpi=150, bbox_inches="tight"); plt.close(f)
    print("MARKER wrote", p)


fig2018()
figbands()
