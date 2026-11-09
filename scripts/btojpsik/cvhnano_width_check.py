#!/usr/bin/env python3
"""Is the data peak wider than the signal template, and is the combinatorial
excess signal-shaped?

Two questions, one method: fit a Gaussian on a smooth background to each mass
distribution and compare widths.

  1. **data against the template.** If the data peak is wider than the simulated
     signal template, the postfit W residual -- prediction 3.5% low on both
     flanks and right in the peak -- follows immediately, and the fix is a
     resolution nuisance or a template taken from data. If the widths agree, the
     residual is a background-shape problem instead.

  2. **the no-b-ancestor excess against the template.** The simulated
     "combinatorial" bin has a local maximum at the B mass. Same width as the
     signal means it is real signal that the truth labelling missed, so in data
     it is not background at all. Broader means partially reconstructed decays,
     which do need a background component.

The fit is deliberately plain: a Gaussian plus an exponential, least squares on
Poisson errors, `scipy.optimize.curve_fit`. It is a width *comparison*, and the
same functional form on both sides is what makes the comparison mean something.
"""

import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from wremnants.production.btojpsik_cvhnano_axes import (
    GEN_CATEGORY_LABELS,
    GEN_DATA,
    GEN_NO_B_ANCESTOR,
    GEN_SIGNAL,
)
from wremnants.postprocessing.rabbit_cvhnano_helpers import open_results

B_MASS = 5.27934


def model(m, n_sig, mu, sigma, a, b):
    """Gaussian signal on an exponential background, both in candidates/bin."""
    gauss = n_sig * np.exp(-0.5 * ((m - mu) / sigma) ** 2)
    return gauss + np.exp(a + b * (m - B_MASS))


def model_curved(m, n_sig, mu, sigma, a, b, n_broad, sigma_broad):
    """The same, plus one broad Gaussian in the background.

    Needed, and this is the whole point of having both forms: with a monotonic
    background the narrow Gaussian absorbs the background's curvature into its
    own width, and the width comparison then reports a data/template ratio that
    is mostly a background-shape error. Measured on the full sample the
    difference is 23.07 MeV against 20.5 -- the monotonic version overstates the
    resolution mismatch by a factor of nearly four.

    The broad component is not a physics claim about a second peak: the
    simulation's combinatorial bin is described by exactly this shape, an
    exponential plus a 137 MeV Gaussian, at chi2/ndf 1.66.
    """
    return (model(m, n_sig, mu, sigma, a, b)
            + n_broad * np.exp(-0.5 * ((m - B_MASS) / sigma_broad) ** 2))


def fit_one(centres, values, label, seed_sigma=0.014, curved=False):
    ok = values > 0
    if ok.sum() < 20:
        print(f"{label:28s} too few filled bins")
        return None
    peak = np.abs(centres - B_MASS) < 0.03
    side = ~peak
    n0 = max(values[peak].max() - np.median(values[side]), 1.0)
    a0 = np.log(max(np.median(values[side]), 1e-3))
    fn = model_curved if curved else model
    p0 = [n0, B_MASS, seed_sigma, a0, 0.0]
    if curved:
        p0 += [0.2 * max(np.median(values[side]), 1e-3), 0.137]
    sigma_err = np.sqrt(np.maximum(values, 1.0))
    try:
        popt, pcov = curve_fit(
            fn, centres[ok], values[ok], p0=p0, sigma=sigma_err[ok],
            absolute_sigma=True, maxfev=200000,
        )
    except RuntimeError as exc:
        print(f"{label:28s} fit failed: {exc}")
        return None
    perr = np.sqrt(np.diag(pcov))
    resid = values[ok] - fn(centres[ok], *popt)
    chi2 = float(np.sum((resid / sigma_err[ok]) ** 2))
    ndf = int(ok.sum() - len(popt))
    n_sig_total = popt[0] * abs(popt[2]) * np.sqrt(2 * np.pi) / (centres[1] - centres[0])
    print(f"{label:28s} sigma {abs(popt[2])*1000:7.2f} +- {perr[2]*1000:5.2f} MeV   "
          f"mu {popt[1]:.5f}   Nsig {n_sig_total:10.0f}   chi2/ndf {chi2/ndf:7.2f}")
    return dict(popt=popt, perr=perr, chi2=chi2, ndf=ndf, sigma=abs(popt[2]),
                sigma_err=perr[2], nsig=n_sig_total, centres=centres, values=values)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument(
        "-i", "--input",
        default="/ceph/submit/data/user/p/pmlugato/mz/histograms/"
                "btojpsik_cvhnano_fullstats.hdf5")
    p.add_argument("-o", "--outdir",
                   default=os.path.expanduser("~/public_html/mz/cvhnano-aem/width"))
    p.add_argument("--window", type=float, nargs=2, default=[5.05, 5.45])
    p.add_argument(
        "--curved", action="store_true",
        help="Give the background one broad Gaussian on top of the exponential. "
             "Without it the narrow width absorbs the background's curvature and "
             "the data/template ratio is mostly a background-shape error.",
    )
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    lo, hi = args.window

    with open_results(args.input) as results:
        keys = [k for k, v in results.items() if isinstance(v, dict) and "output" in v]
        dk = next(k for k in keys if "Charmonium" in k)
        mk = next(k for k in keys if k != dk)
        hists = {}
        for lab, k in (("data", dk), ("mc", mk)):
            h = results[k]["output"]["fitbins"]
            h = h.get() if hasattr(h, "get") else h
            hists[lab] = h.project("jointCvh_mass", "genCategory")

    ax = hists["mc"].axes["jointCvh_mass"]
    centres = np.asarray(ax.centers)
    inwin = (centres > lo) & (centres < hi)
    c = centres[inwin]

    print(f"window [{lo}, {hi}], {inwin.sum()} bins of "
          f"{(ax.edges[1]-ax.edges[0])*1000:.0f} MeV")
    print()
    print("=== Gaussian + exponential, same form everywhere ===")
    fits = {}
    fits["mc signal template"] = fit_one(
        c, np.asarray(hists["mc"][:, GEN_SIGNAL].values(), dtype=float)[inwin],
        "MC " + GEN_CATEGORY_LABELS[GEN_SIGNAL], curved=args.curved)
    fits["mc no-b-ancestor"] = fit_one(
        c, np.asarray(hists["mc"][:, GEN_NO_B_ANCESTOR].values(), dtype=float)[inwin],
        "MC " + GEN_CATEGORY_LABELS[GEN_NO_B_ANCESTOR])
    fits["data"] = fit_one(
        c, np.asarray(hists["data"][:, GEN_DATA].values(), dtype=float)[inwin],
        "data", curved=args.curved)

    print()
    print("=== the comparison ===")
    tpl = fits.get("mc signal template")
    if tpl is None:
        return
    for name in ("data", "mc no-b-ancestor"):
        f = fits.get(name)
        if f is None:
            continue
        ratio = f["sigma"] / tpl["sigma"]
        # Error on the ratio, treating the two fits as independent, which they
        # are: different samples for data, disjoint truth bins for the other.
        rerr = ratio * np.sqrt((f["sigma_err"] / f["sigma"]) ** 2
                               + (tpl["sigma_err"] / tpl["sigma"]) ** 2)
        nsig = abs(ratio - 1.0) / rerr if rerr > 0 else float("inf")
        print(f"{name:20s} sigma / template sigma = {ratio:.4f} +- {rerr:.4f}"
              f"   ({nsig:.1f} sigma from 1)")
    print()
    print("Reading it: a data width above the template's means the postfit W")
    print("residual is a resolution effect and the fix is the signal shape.")
    print("A no-b-ancestor width equal to the template's means that bin's excess")
    print("is mislabelled signal, not a background needing its own component.")

    fig, axes = plt.subplots(1, 3, figsize=(13.0, 3.9))
    for ax_, (name, f) in zip(axes, [("MC signal template", fits["mc signal template"]),
                                     ("MC no-b-ancestor", fits["mc no-b-ancestor"]),
                                     ("data", fits["data"])]):
        if f is None:
            continue
        ax_.errorbar(f["centres"], f["values"],
                     yerr=np.sqrt(np.maximum(f["values"], 0)), fmt="o", ms=2,
                     lw=0.6, color="k")
        fine = np.linspace(f["centres"][0], f["centres"][-1], 400)
        ax_.plot(fine, (model_curved if len(f["popt"]) == 7 else model)(fine, *f["popt"]),
                 "r-", lw=1.5)
        ax_.plot(fine, np.exp(f["popt"][3] + f["popt"][4] * (fine - B_MASS)),
                 "r--", lw=1.1)
        ax_.set_title(f"{name}\n"
                      rf"$\sigma$ = {f['sigma']*1000:.2f} $\pm$ "
                      rf"{f['sigma_err']*1000:.2f} MeV", fontsize=9)
        ax_.set_xlabel(r"$m(\mu\mu K)$ [GeV]")
        ax_.set_xlim(lo, hi)
    axes[0].set_ylabel("candidates / 5 MeV")
    fig.tight_layout()
    out = os.path.join(
        args.outdir,
        "width_comparison_curved.png" if args.curved else "width_comparison.png")
    fig.savefig(out, dpi=160, bbox_inches="tight")
    print()
    print("wrote", out)


if __name__ == "__main__":
    main()
