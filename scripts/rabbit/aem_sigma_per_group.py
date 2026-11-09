#!/usr/bin/env python3
"""Per-eta posterior sigma of A/e/M (prior units) for two or more fits, with a
ratio panel to the first. Also prints the summary table used in the decks:
sigma(mZ), mean sigma(A/e/M), rho(A_i,e_i). Run inside the container.

    python scripts/rabbit/aem_sigma_per_group.py --fits 'Z+Jpsi=ref.hdf5' 'Z+Jpsi+B=fit.hdf5' -o out.png [--table out.md]
"""
import argparse, sys
import numpy as np
sys.path.insert(0, "rabbit"); sys.path.insert(0, "wums")
from rabbit import io_tools
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--fits", nargs="+", required=True, help="LABEL=fit.hdf5; the first is the ratio reference")
ap.add_argument("--nEta", type=int, default=24)
ap.add_argument("-o", "--output", required=True)
ap.add_argument("--table", default=None)
ap.add_argument("--title", default="")
ap.add_argument("--poi", default="massShiftZ100MeV", help="POI for the sigma/rho columns; skipped if absent from a fit")
ap.add_argument("--poiScale", type=float, default=100.0)
ap.add_argument("--poiLabel", default="mZ")
a = ap.parse_args()
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]

def load(path):
    r = io_tools.get_fitresult(path)
    names = np.array([str(x) for x in r["parms"].get().axes["parms"]]); cov = np.asarray(r["cov"].get().values())
    sig = np.sqrt(np.diag(cov)); rho = cov / np.outer(sig, sig)
    ix = lambda n: int(np.flatnonzero(names == n)[0])
    idx = {p: [ix(f"{p}{i}") for i in range(a.nEta)] for p in "AeM"}
    has = a.poi in names
    iz = ix(a.poi) if has else None
    return dict(sig={p: sig[idx[p]] for p in "AeM"}, mz=a.poiScale * sig[iz] if has else np.nan,
                rhoAe=np.array([rho[i, j] for i, j in zip(idx["A"], idx["e"])]),
                rho_mz={p: (rho[iz, idx[p]].mean() if has else np.nan) for p in "AeM"})

fits = [(s.split("=", 1)[0], load(s.split("=", 1)[1])) for s in a.fits]
fig, axs = plt.subplots(2, 3, figsize=(13, 6.2), sharex=True, gridspec_kw={"height_ratios": [2.2, 1], "hspace": 0.08})
x = np.arange(a.nEta)
for k, p in enumerate("AeM"):
    ref = fits[0][1]["sig"][p]
    for (lab, f), c in zip(fits, C):
        axs[0, k].plot(x, f["sig"][p], "o-", color=c, ms=4, lw=1.5, label=lab)
        axs[1, k].plot(x, f["sig"][p] / ref, "o-", color=c, ms=4, lw=1.5)
    axs[0, k].axhline(1.0, color="#52514e", ls=":", lw=0.8); axs[0, k].set_title(f"sigma({p}) per eta group, prior units", fontsize=11)
    axs[1, k].set_xlabel("eta group"); axs[1, k].set_ylabel(f"ratio to {fits[0][0]}"); axs[0, k].set_ylim(0, 1.05)
    for ax in axs[:, k]: ax.grid(color="#e6e6e3")
axs[0, 0].legend(frameon=False, fontsize=9.5)
if a.title: fig.suptitle(a.title, fontsize=11.5)
fig.savefig(a.output, dpi=130, bbox_inches="tight"); print("wrote", a.output)

PL = a.poiLabel
L = [f"| fit | σ({PL}) [MeV] | mean σ(A_i) | mean σ(e_i) | mean σ(M_i) | mean ρ(A_i,e_i) | ρ({PL},A_i) | ρ({PL},e_i) |", "|---|--:|--:|--:|--:|--:|--:|--:|"]
for lab, f in fits:
    L.append(f"| {lab} | {f['mz']:.4f} | {f['sig']['A'].mean():.3f} | {f['sig']['e'].mean():.3f} | {f['sig']['M'].mean():.3f} | {f['rhoAe'].mean():+.3f} | {f['rho_mz']['A']:+.3f} | {f['rho_mz']['e']:+.3f} |")
L.append(""); L.append("σ of scale parameters in units of the prior width; ρ(A_i,e_i) averaged over the 24 eta groups.")
print("\n".join(L))
if a.table: open(a.table, "w").write("\n".join(L) + "\n")
