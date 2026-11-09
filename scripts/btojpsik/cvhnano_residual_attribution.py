#!/usr/bin/env python3
"""Which component owns the postfit residual: the signal shape or the background?

The W-shaped residual could in principle come from either -- a signal template
that is too narrow, or a background whose curvature is wrong -- and the two call
for different fixes. Two tests, both on the postfit per-process histograms the
fit already wrote.

  1. **Where the residual sits relative to the signal.** If the residual tracks
     the signal's own derivative structure it is the signal shape; if it tracks
     the background it is the background.

  2. **A widened signal template, tested directly.** Convolve the postfit signal
     with a Gaussian of the width the peak comparison says is missing,
     sqrt(sigma_data^2 - sigma_template^2), rescale to the same integral, and ask
     how much of the residual it removes. This is not a fit, so it cannot
     flatter itself by absorbing the residual into a free parameter: the
     smearing width and the normalisation are both fixed from outside.
"""

import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import rabbit.io_tools as iot

B_MASS = 5.27934


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument(
        "fit",
        default="/ceph/submit/data/user/p/pmlugato/mz/fitresults/"
                "cvhnano_full_exp_data_eta12_normgrouped_aem_full_fit.hdf5",
        nargs="?",
    )
    p.add_argument("-o", "--outdir",
                   default=os.path.expanduser("~/public_html/mz/cvhnano-aem/attribution"))
    # From cvhnano_width_check.py on the full sample. Passed in rather than
    # refitted here so that this script cannot quietly choose a width that
    # flatters the conclusion.
    p.add_argument("--sigmaData", type=float, default=23.07e-3)
    p.add_argument("--sigmaTemplate", type=float, default=19.60e-3)
    return p.parse_args()


def smear(values, centres, width):
    """Gaussian convolution on a uniform grid, normalisation preserved."""
    step = centres[1] - centres[0]
    if width <= 0:
        return values.copy()
    half = int(np.ceil(4.0 * width / step))
    off = np.arange(-half, half + 1) * step
    kern = np.exp(-0.5 * (off / width) ** 2)
    kern /= kern.sum()
    # 'same' truncates at the window edges, which loses a little of the tail;
    # the integral is restored explicitly so the comparison is at equal yield.
    out = np.convolve(values, kern, mode="same")
    total = values.sum()
    return out * (total / out.sum()) if out.sum() > 0 else out


def best_norms(d, comps, var):
    """Least-squares signal and background normalisations for fixed shapes.

    The question is about SHAPE, so the normalisations must not be allowed to
    carry the answer. Profiling them out analytically at every trial width is
    what makes the comparison a shape comparison: a 1% overall offset -- which
    the fit file's own postfit histogram has relative to the plot-ready one,
    because the two treat the bin-by-bin statistical nuisances differently --
    would otherwise dominate the chi-square and decide the scan.
    """
    A = np.stack(comps, axis=1) / np.sqrt(var)[:, None]
    b = d / np.sqrt(var)
    coef, *_ = np.linalg.lstsq(A, b, rcond=None)
    pred = np.stack(comps, axis=1) @ coef
    return coef, pred, float(np.sum((d - pred) ** 2 / var))


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    fitres = iot.get_fitresult(args.fit)

    key = "Project btojpsik_cvhnano jointCvh_mass"
    base = fitres["mappings"][key]["channels"]["btojpsik_cvhnano"]

    def g(name):
        o = base[name]
        return o.get() if hasattr(o, "get") else o

    nobs = g("hist_nobs")
    post = g("hist_postfit")
    procs = list(post.axes["processes"])
    centres = np.asarray(nobs.axes["jointCvh_mass"].centers)
    d = np.asarray(nobs.values(), dtype=float)
    per = {p: np.asarray(post[{"processes": p}].values(), dtype=float) for p in procs}

    sig_name = next(p for p in procs if "Bkg" not in p)
    bkg_name = next(p for p in procs if "Bkg" in p)
    sig, bkg = per[sig_name], per[bkg_name]
    var = np.maximum(d, 1.0)

    print(f"processes: {procs}")
    print(f"signal = {sig_name}   background = {bkg_name}")
    print(f"postfit signal {sig.sum():12.0f}   background {bkg.sum():12.0f}   "
          f"data {d.sum():12.0f}")
    print()
    print("The per-process histograms live only in the fit file, whose postfit")
    print("total is 6 141 387 against the plot-ready file's 6 204 940 -- a 1.0%")
    print("difference in how the bin-by-bin statistical nuisances are folded in.")
    print("Both normalisations are therefore profiled out at every trial width,")
    print("so that offset cannot decide the answer.")
    print()

    print("=== the shape scan: extra Gaussian width on the SIGNAL only ===")
    extra = np.sqrt(max(args.sigmaData**2 - args.sigmaTemplate**2, 0.0))
    print(f"the peak comparison says the missing width is "
          f"sqrt({args.sigmaData*1e3:.2f}^2 - {args.sigmaTemplate*1e3:.2f}^2) = "
          f"{extra*1e3:.2f} MeV")
    print(f"{'smearing [MeV]':>16s} {'chi2':>12s} {'chi2/bin':>10s} "
          f"{'sig scale':>10s} {'bkg scale':>10s} {'max|r-1|':>10s}")
    best = None
    for frac in (0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0):
        w = frac * extra
        comps = [smear(sig, centres, w), bkg]
        coef, pred, c2 = best_norms(d, comps, var)
        mx = float(np.max(np.abs(d / np.where(pred > 0, pred, 1.0) - 1.0)))
        print(f"{w*1e3:16.2f} {c2:12.1f} {c2/len(d):10.1f} {coef[0]:10.4f} "
              f"{coef[1]:10.4f} {mx:10.4f}")
        if best is None or c2 < best[2]:
            best = (w, coef, c2, pred)
    c2_0 = best_norms(d, [sig, bkg], var)[2]
    print()
    print(f"unsmeared, normalisations profiled: chi2 {c2_0:.0f} "
          f"({c2_0/len(d):.1f}/bin)")
    print(f"best: {best[0]*1e3:.2f} MeV -> chi2 {best[2]:.0f} "
          f"({best[2]/len(d):.1f}/bin), a factor {c2_0/max(best[2],1e-9):.2f}")
    print()

    # The same scan, but smearing the BACKGROUND's curvature instead: a control.
    # If widening the signal is the fix, adding curvature to the background
    # should do less; if not, the residual is not a signal-shape effect.
    print("=== control: one extra free shape in the background instead ===")
    # A quadratic in mass, added to the two-component fit. Three normalisations
    # profiled, nothing about the signal touched.
    x = (centres - B_MASS) / 0.2
    for label, extra_comps in (
        ("+ linear", [x * bkg]),
        ("+ linear + quadratic", [x * bkg, x**2 * bkg]),
        ("+ broad Gaussian (137 MeV)", [np.exp(-0.5 * (x * 0.2 / 0.137) ** 2)]),
    ):
        coef, pred, c2 = best_norms(d, [sig, bkg] + extra_comps, var)
        mx = float(np.max(np.abs(d / np.where(pred > 0, pred, 1.0) - 1.0)))
        print(f"{label:30s} chi2 {c2:12.1f} ({c2/len(d):8.1f}/bin)  "
              f"max|r-1| {mx:.4f}")
    print()
    print("And both together, which is the model actually being proposed:")
    for frac in (0.0, 0.5, 1.0):
        w = frac * extra
        comps = [smear(sig, centres, w), bkg, x * bkg, x**2 * bkg]
        coef, pred, c2 = best_norms(d, comps, var)
        mx = float(np.max(np.abs(d / np.where(pred > 0, pred, 1.0) - 1.0)))
        print(f"  signal +{w*1e3:5.2f} MeV, background + linear + quadratic: "
              f"chi2 {c2:10.1f} ({c2/len(d):7.1f}/bin)  max|r-1| {mx:.4f}")

    tot = best_norms(d, [sig, bkg], var)[1]

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.3))
    ax = axes[0]
    ax.errorbar(centres, d, yerr=np.sqrt(np.maximum(d, 0)), fmt="o", ms=2.4,
                color="k", lw=0.7, label="data")
    ax.step(centres, tot, where="mid", color="tab:red", lw=1.6,
            label="postfit shapes, norms profiled")
    ax.step(centres, best[3], where="mid", color="tab:blue", lw=1.6,
            label=f"signal smeared by {best[0]*1e3:.1f} MeV")
    ax.step(centres, bkg, where="mid", color="0.55", lw=1.1, ls="--",
            label="postfit background")
    ax.set_xlabel(r"$m(\mu\mu K)$ [GeV]")
    ax.set_ylabel("candidates / 5 MeV")
    ax.legend(fontsize=8)

    ax = axes[1]
    ax.axhline(1.0, color="0.5", lw=0.9)
    c2_0 = float(np.sum((d - tot) ** 2 / np.maximum(d, 1.0)))
    ax.errorbar(centres, d / np.where(tot > 0, tot, 1.0),
                yerr=np.sqrt(np.maximum(d, 0)) / np.where(tot > 0, tot, 1.0),
                fmt="o", ms=2.4, color="tab:red", lw=0.7,
                label=f"postfit shapes, χ²/bin {c2_0/len(d):.1f}")
    ax.errorbar(centres, d / np.where(best[3] > 0, best[3], 1.0),
                yerr=np.sqrt(np.maximum(d, 0)) / np.where(best[3] > 0, best[3], 1.0),
                fmt="s", ms=2.4, color="tab:blue", lw=0.7,
                label=f"signal smeared, χ²/bin {best[2]/len(d):.1f}")
    ax.set_xlabel(r"$m(\mu\mu K)$ [GeV]")
    ax.set_ylabel("data / prediction")
    ax.set_ylim(0.94, 1.08)
    ax.legend(fontsize=8)

    fig.tight_layout()
    out = os.path.join(args.outdir, "residual_attribution.png")
    fig.savefig(out, dpi=160, bbox_inches="tight")
    print()
    print("wrote", out)


if __name__ == "__main__":
    main()
