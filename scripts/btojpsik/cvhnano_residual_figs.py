#!/usr/bin/env python3
"""Per-cell prefit and postfit residuals on one canvas.

rabbit writes one PNG per cell, which is the right output for a readout and the
wrong one for reading a *pattern*. The pattern is the whole point here: a model
that is too narrow leaves the same W-shaped residual in every cell, and a model
that is merely mis-normalised leaves a flat offset. Sixteen separate plots
cannot be compared by eye; one canvas can.

Also prints, per cell, the chi-square the plots do not carry: on a data fit run
with --noChi2 the value rabbit prints on the postfit plot is uninitialised
garbage (-1.1e15 was observed), so the number has to be recomputed here.
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
    p.add_argument("plotready", nargs="+")
    p.add_argument("--labels", nargs="*", default=None)
    p.add_argument("-o", "--outdir",
                   default=os.path.expanduser("~/public_html/mz/cvhnano-aem/residuals"))
    p.add_argument("--channel", default="btojpsik_cvhnano")
    p.add_argument("--tag", default="cells")
    p.add_argument(
        "--onlyPt", type=int, default=None,
        help="Keep one bachelor-pT bin only. A 4x4 grid of cells is too tall to "
             "read on a slide: at 16:9 it displays 6.9 in wide and the panel "
             "titles land at under 5 pt. One pT bin gives a 4x2 grid, which is "
             "the same information at a legible size.",
    )
    return p.parse_args()


def get(fitres, key):
    """Walk a nested key path and materialise the leaf.

    `H5PickleProxy.get()` is the read; a plain dict also has `.get`, so the
    proxy is identified by its type rather than by duck-typing on the name.
    """
    o = fitres
    for part in key:
        o = o[part]
    if isinstance(o, dict):
        return o
    return o.get() if hasattr(o, "get") else o


def cell_mapping(fitres):
    """The mapping whose histogram carries the cell axes, and its axis names."""
    best, best_n = None, -1
    for name in fitres["mappings"].keys():
        h = get(fitres, ("mappings", name, "channels"))
        for ch in h.keys():
            hh = get(fitres, ("mappings", name, "channels", ch, "hist_postfit_inclusive"))
            n = len(hh.axes) - 1
            if n > best_n:
                best, best_n = (name, ch), n
    return best


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    labels = args.labels or [os.path.basename(f).replace("cvhnano_", "").replace("_plotready.hdf5", "")
                             for f in args.plotready]

    for path, label in zip(args.plotready, labels):
        fitres = iot.get_fitresult(path)
        name, ch = cell_mapping(fitres)
        base = ("mappings", name, "channels", ch)
        # hist_nobs, NOT hist_data_obs. On an Asimov fit both are present and
        # they are different things: hist_nobs is the pseudodata that was
        # actually fitted, hist_data_obs is the real observation, which the
        # Asimov fit never saw. Using data_obs silently compares an Asimov
        # postfit against the data and reports the *prefit* disagreement as if
        # it were postfit -- prefit and postfit chi-square came out equal to
        # four digits, which is what gave it away.
        data = get(fitres, base + ("hist_nobs",))
        pre = get(fitres, base + ("hist_prefit_inclusive",))
        post = get(fitres, base + ("hist_postfit_inclusive",))

        cell_axes = [a for a in data.axes if a.name != "jointCvh_mass"]
        mass = data.axes["jointCvh_mass"]
        centres = np.asarray(mass.centers)

        if args.onlyPt is not None:
            sel = {"fit_bachelor_pt": args.onlyPt}
            data, pre, post = data[sel], pre[sel], post[sel]
            cell_axes = [a for a in cell_axes if a.name != "fit_bachelor_pt"]

        shape = tuple(a.size for a in cell_axes)
        ncell = int(np.prod(shape))
        dv = data.values().reshape(ncell, mass.size)
        # Poisson errors on the observation. hist_data_obs carries variances on
        # an Asimov dataset but not necessarily on the observation, so sqrt(N)
        # is used throughout for one consistent definition.
        de = np.sqrt(np.maximum(dv, 0.0))
        prv = pre.values().reshape(ncell, mass.size)
        pov = post.values().reshape(ncell, mass.size)

        ncol = min(4, ncell)
        nrow = int(np.ceil(ncell / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(3.05 * ncol, 2.35 * nrow),
                                 sharex=True, sharey=True, squeeze=False)

        print(f"=== {label}   {ncell} cells, {mass.size} mass bins")
        print(f"{'cell':22s} {'chi2/ndf pre':>13s} {'chi2/ndf post':>14s} "
              f"{'max |resid|':>12s}")
        tot_pre = tot_post = 0.0
        for icell in range(ncell):
            idx = np.unravel_index(icell, shape)
            ax = axes[icell // ncol][icell % ncol]
            ok = pov[icell] > 0
            ratio = np.divide(dv[icell], pov[icell], out=np.ones(mass.size), where=ok)
            rerr = np.divide(de[icell], pov[icell], out=np.zeros(mass.size), where=ok)
            rpre = np.divide(dv[icell], prv[icell], out=np.ones(mass.size),
                             where=prv[icell] > 0)
            ax.axhline(1.0, color="0.5", lw=0.9)
            ax.axvspan(B_MASS - 0.03, B_MASS + 0.03, color="0.9", zorder=0)
            ax.step(centres, rpre, where="mid", color="tab:red", lw=1.0,
                    label="prefit" if icell == 0 else None)
            ax.errorbar(centres, ratio, yerr=rerr, fmt="o", ms=1.8, lw=0.7,
                        color="k", label="postfit" if icell == 0 else None)
            cellname = ", ".join(
                f"{a.name.replace('fit_bachelor_', '')}{i}"
                for a, i in zip(cell_axes, idx)
            )
            if args.onlyPt is not None:
                cellname = f"pt{args.onlyPt}, " + cellname
            ax.set_title(cellname, fontsize=8)
            ax.set_ylim(0.90, 1.12)
            ax.set_xlim(centres[0], centres[-1])
            ax.tick_params(labelsize=7)

            var = np.maximum(dv[icell], 1.0)
            c_pre = float(np.sum((dv[icell] - prv[icell]) ** 2 / var))
            c_post = float(np.sum((dv[icell] - pov[icell]) ** 2 / var))
            tot_pre += c_pre
            tot_post += c_post
            print(f"{cellname:22s} {c_pre/mass.size:13.1f} {c_post/mass.size:14.2f} "
                  f"{np.max(np.abs(ratio - 1.0)):12.4f}")

        ndf = ncell * mass.size
        print(f"{'TOTAL':22s} {tot_pre/ndf:13.1f} {tot_post/ndf:14.2f}   "
              f"(ndf {ndf}, no parameter count subtracted)")

        # The same residual on the cell-summed mass projection. It is the more
        # demanding test, not a weaker one: summing the cells cancels nothing
        # because the mismatch has the same sign in every cell, while the
        # statistical error falls as sqrt(16). A shape error that is coherent
        # across cells shows up here and can hide per cell.
        dsum, prsum, posum = dv.sum(0), prv.sum(0), pov.sum(0)
        vsum = np.maximum(dsum, 1.0)
        c_pre = float(np.sum((dsum - prsum) ** 2 / vsum))
        c_post = float(np.sum((dsum - posum) ** 2 / vsum))
        print(f"{'cell-summed':22s} {c_pre/mass.size:13.1f} "
              f"{c_post/mass.size:14.2f}   (ndf {mass.size})")
        print(f"{'':22s} max |data/pred - 1| postfit, summed: "
              f"{np.max(np.abs(dsum/np.where(posum>0,posum,1) - 1.0)):.4f}")

        for j in range(ncell, nrow * ncol):
            axes[j // ncol][j % ncol].set_visible(False)
        for j in range(ncol):
            axes[nrow - 1][j].set_xlabel(r"$m(\mu\mu K)$ [GeV]", fontsize=8)
        for i in range(nrow):
            axes[i][0].set_ylabel("data / pred.", fontsize=8)
        axes[0][0].legend(fontsize=7, loc="upper left")
        fig.tight_layout()
        out = os.path.join(args.outdir, f"residual_{args.tag}_{label}.png")
        fig.savefig(out, dpi=160, bbox_inches="tight")
        print("wrote", out)
        print()


if __name__ == "__main__":
    main()
