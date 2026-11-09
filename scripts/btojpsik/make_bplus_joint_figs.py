#!/usr/bin/env python3
"""What adding B+ does to the joint Z + J/psi Asimov fit.

Two panels, because the marginal A/e/M errors are NOT the figure of merit: the
reference already constrains them to ~0.1 of the prior, and B+ standalone is far
weaker. B+'s claim is the lever arm -- k = 1/pT reaching 1.0 where J/psi stops at
0.33 -- so what has to move is the A-e degeneracy.

    python scripts/btojpsik/make_bplus_joint_figs.py COMBINED.hdf5 REFERENCE.hdf5
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rabbit.io_tools as io

OUT = os.path.expanduser("~/public_html/slides/assets_btojpsik_overview")
INK, ACC, BLU = "#1F1F1F", "#C0392B", "#2E86C1"


def read(path, n_eta):
    fr = io.get_fitresult(path)
    hp = fr["parms"].get()
    names = [x.decode() if isinstance(x, bytes) else str(x) for x in hp.axes[0]]
    sig = {n: float(np.sqrt(v)) for n, v in zip(names, hp.variances())}
    hc = fr["cov"].get()
    cn = [x.decode() if isinstance(x, bytes) else str(x) for x in hc.axes[0]]
    idx = {n: i for i, n in enumerate(cn)}
    C = hc.values()
    rho = np.full(n_eta, np.nan)
    for i in range(n_eta):
        ia, ie = idx.get(f"A{i}"), idx.get(f"e{i}")
        if ia is None or ie is None:
            continue
        den = np.sqrt(C[ia, ia] * C[ie, ie])
        if den > 0:
            rho[i] = C[ia, ie] / den
    return sig, rho


def main():
    p = argparse.ArgumentParser()
    p.add_argument("combined")
    p.add_argument("reference")
    p.add_argument("--nEta", type=int, default=24)
    a = p.parse_args()
    sn, rn = read(a.combined, a.nEta)
    sr, rr = read(a.reference, a.nEta)
    x = np.arange(a.nEta)
    os.makedirs(OUT, exist_ok=True)

    # --- rho(A,e)
    f, ax = plt.subplots(figsize=(9.8, 4.5))
    ax.plot(x, rr, "o-", color=INK, ms=6, lw=1.8, label="Z + J/ψ")
    ax.plot(x, rn, "s-", color=ACC, ms=6, lw=2.0, label="Z + J/ψ + B$^+$")
    ax.set_xlabel("nuisance η bin")
    ax.set_ylabel(r"$\rho(A_i,\, e_i)$")
    ax.legend(frameon=False, fontsize=10)
    ax.grid(alpha=0.15)
    q = os.path.join(OUT, "bplus_rho_Ae.png")
    f.savefig(q, dpi=150, bbox_inches="tight"); plt.close(f)
    print("MARKER wrote", q)

    # --- sigma ratio per family
    f, ax = plt.subplots(figsize=(9.8, 4.5))
    for param, col, mk in (("A", INK, "o"), ("e", ACC, "s"), ("M", BLU, "^")):
        r = np.array([sn.get(f"{param}{i}", np.nan) / sr[f"{param}{i}"]
                      if f"{param}{i}" in sr and sr[f"{param}{i}"] else np.nan
                      for i in range(a.nEta)])
        ax.plot(x, r, mk + "-", color=col, ms=6, lw=1.7, label=f"${param}$")
    ax.axhline(1.0, color="grey", lw=0.9, ls=":")
    ax.set_xlabel("nuisance η bin")
    ax.set_ylabel("σ with B$^+$  /  σ without")
    ax.legend(frameon=False, fontsize=10, ncol=3)
    ax.grid(alpha=0.15)
    q = os.path.join(OUT, "bplus_sigma_ratio.png")
    f.savefig(q, dpi=150, bbox_inches="tight"); plt.close(f)
    print("MARKER wrote", q)
    print(f"MARKER rho(A,e) mean: reference {np.nanmean(rr):+.4f}  "
          f"with B+ {np.nanmean(rn):+.4f}  change {np.nanmean(rn-rr):+.4f}")
    solo_panel(a.nEta)
    return 0




def solo_panel(n_eta=24):
    """Why the gain is small: B+ alone against Z + J/psi, same basis, log scale."""
    P = "/ceph/submit/data/user/p/pmlugato/mz/fitresults/"
    def sig(p):
        fr = io.get_fitresult(P + p)
        h = fr["parms"].get()
        nm = [x.decode() if isinstance(x, bytes) else str(x) for x in h.axes[0]]
        return dict(zip(nm, np.sqrt(h.variances())))
    solo = sig("bplus_solo_joint24_asimov_fit.hdf5")
    floor = sig("bplus_solo_joint24_nobbb_fit.hdf5")
    ref = sig("joint_mz_jpsi_repro_fit.hdf5")
    comb = sig("joint_mz_jpsi_bplus_btojpsik_asimov_fit.hdf5")
    x = np.arange(n_eta)
    g = lambda d, p: np.array([d[f"{p}{i}"] for i in range(n_eta)])

    f, axs = plt.subplots(1, 3, figsize=(11.4, 3.9), sharey=True)
    for ax, param in zip(axs, ("A", "e", "M")):
        ax.plot(x, g(solo, param), "s-", color=ACC, ms=4, lw=1.6,
                label="B$^+$ alone")
        ax.plot(x, g(floor, param), "s--", color=ACC, ms=3, lw=1.2, alpha=0.6,
                label="B$^+$ alone, no MC-stat")
        ax.plot(x, g(ref, param), "o-", color=INK, ms=4, lw=1.6, label="Z + J/ψ")
        ax.plot(x, g(comb, param), "^-", color=BLU, ms=4, lw=1.6, label="all three")
        ax.axhline(1.0, color="grey", lw=0.8, ls=":")
        ax.set_yscale("log")
        ax.set_title(f"${param}$", fontsize=12)
        ax.set_xlabel("nuisance η bin")
        ax.grid(alpha=0.15)
    axs[0].set_ylabel("σ  [prior units]")
    axs[0].legend(frameon=False, fontsize=8, loc="lower left")
    q = os.path.join(OUT, "bplus_solo_vs_joint.png")
    f.savefig(q, dpi=150, bbox_inches="tight"); plt.close(f)
    print("MARKER wrote", q)


if __name__ == "__main__":
    sys.exit(main())
