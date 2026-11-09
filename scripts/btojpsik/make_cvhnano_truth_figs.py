"""The two figures the central plotting chain cannot make.

The stacks come from `makeDataMCStackPlot.py --noData` and the data prefit from
`btojpsik_cvhnano_prefit.py`. This adds the two comparisons that decide what the
prefit residual means.

**fig_bkgmodel** -- the same data, the same signal template, three background
choices: the true simulated non-signal shape, an exponential, and a polynomial.
This is the measurement that separates "the background model is wrong" from
"the signal template is wrong", and it answers the first.

**fig_width** -- the signal template against data with the *fitted true
background* subtracted. Not a linear subtraction: the simulated background is a
broad maximum sitting on the B mass, so a straight line under a 200 MeV window
leaves a positive residual that inflates the apparent data width by a factor
1.7 and is pure artefact. The first version of this figure made exactly that
mistake. The control that catches it is in `widths.json`: apply the identical
treatment to the simulated TOTAL, whose width must then agree with the data's.

    python scripts/btojpsik/make_cvhnano_truth_figs.py <histmaker.hdf5>
"""

import argparse
import json
import os

import matplotlib as mpl
import numpy as np

mpl.use("Agg")
import matplotlib.pyplot as plt

from wremnants.postprocessing import rabbit_cvhnano_helpers as H
from wremnants.production.btojpsik_cvhnano_axes import (
    GEN_DATA,
    GEN_MC_CATEGORIES,
    GEN_SIGNAL,
    GEN_TEMPLATE_CATEGORIES,
)

DATA = "Charmonium_2016PostVFP_cvhNano"
MC = "InclusiveBToJpsiX_2016PostVFP"
B_MASS, JPSI_MASS = 5.27934, 3.09690
FIT = (5.15, 5.42)
WIDTH_WINDOW = (5.19, 5.37)
C_DATA, C_MC, C_TRUE, C_EXP = "black", "#2a78d6", "#27ae60", "#c0392b"
SURFACE, INK_MUTED = "#fcfcfb", "#6b6b6b"

# (key, histogram, label, peak centre, width-fit window)
#
# The dimuon entries are the point, not extras. If the B width mismatch were
# generic tracking resolution, the unconstrained DIMUON mass -- the same two muon
# tracks, no bachelor -- would be mismatched by a similar factor.
ARMS = [
    ("jointCvh", "nominal_BuJpsiK_jointCvhMass", "joint N-body CVH", B_MASS, WIDTH_WINDOW),
    ("kvfCvh", "nominal_BuJpsiK_kvfCvhMass", "KVF on CVH tracks", B_MASS, WIDTH_WINDOW),
    ("kvfRaw", "nominal_BuJpsiK_kvfRawMass", "KVF on raw tracks", B_MASS, WIDTH_WINDOW),
    ("raw", "nominal_BuJpsiK_mass", "raw, pre-fit", B_MASS, WIDTH_WINDOW),
    ("dimuonCvh", "nominal_BuJpsiK_kvfCvhDimuonMass", "dimuon, CVH tracks",
     JPSI_MASS, (3.01, 3.19)),
    ("dimuonRaw", "nominal_BuJpsiK_kvfRawDimuonMass", "dimuon, raw tracks",
     JPSI_MASS, (3.01, 3.19)),
]


def gauss_width(c, y, centre, quad=False):
    """Gaussian width on a linear (optionally quadratic) background.

    The width is scanned and everything else solves exactly at each grid point,
    so this is deterministic and cannot fail to converge. 0.125 MeV steps,
    against a 5 MeV bin width.
    """
    best = None
    for sigma in np.linspace(0.004, 0.060, 449):
        g = np.exp(-0.5 * ((c - centre) / sigma) ** 2)
        cols = [g, np.ones_like(c), c - centre]
        if quad:
            cols.append((c - centre) ** 2)
        A = np.vstack(cols).T
        p, *_ = np.linalg.lstsq(A, y, rcond=None)
        res = float(np.sum((A @ p - y) ** 2))
        if best is None or res < best[0]:
            best = (res, sigma, p)
    return best[1], best[2]


def lsq(components, target, sigma):
    A = np.vstack([c / sigma for c in components]).T
    p, *_ = np.linalg.lstsq(A, target / sigma, rcond=None)
    model = sum(a * c for a, c in zip(p, components))
    return p, model, float(np.sum(((model - target) / sigma) ** 2))


def get(o, name):
    h = o[name]
    return h.get() if hasattr(h, "get") else h


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("infile")
    p.add_argument(
        "-o",
        "--outdir",
        default="/home/submit/pmlugato/public_html/slides/assets_cvhnano_truth",
    )
    args = p.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    out = {"widths": {}, "background_models": {}}

    with H.open_results(args.infile) as r:
        do, mo = r[DATA]["output"], r[MC]["output"]

        # ---- widths, data against the simulated TOTAL under one treatment ----
        for key, hname, label, centre, window in ARMS:
            hd, hm = get(do, hname), get(mo, hname)
            c = hd.axes[0].centers
            d = hd.values()[:, GEN_DATA]
            sig = hm.values()[:, GEN_SIGNAL]
            total = hm.values()[:, GEN_MC_CATEGORIES].sum(axis=1)
            w = (c > window[0]) & (c < window[1])
            row = {"label": label}
            for tag, quad in (("lin", False), ("linquad", True)):
                sd, _ = gauss_width(c[w], d[w], centre, quad)
                st, _ = gauss_width(c[w], total[w], centre, quad)
                ss, _ = gauss_width(c[w], sig[w], centre, quad)
                row[tag] = dict(
                    data_mev=sd * 1000.0,
                    mc_total_mev=st * 1000.0,
                    mc_signal_only_mev=ss * 1000.0,
                    ratio_to_total=sd / st if st > 0 else float("nan"),
                )
            out["widths"][key] = row
            r0 = row["lin"]
            print(f"{label:20s} data {r0['data_mev']:6.2f}  MC total "
                  f"{r0['mc_total_mev']:6.2f}  ratio {r0['ratio_to_total']:5.3f}   "
                  f"(MC signal alone {r0['mc_signal_only_mev']:6.2f})")

        # ---- the background-model comparison, on the nominal arm -------------
        hd, hm = get(do, "nominal_BuJpsiK_jointCvhMass"), get(mo, "nominal_BuJpsiK_jointCvhMass")
        c, edges = hd.axes[0].centers, hd.axes[0].edges
        m = (c >= FIT[0]) & (c <= FIT[1])
        d, dv = hd.values()[:, GEN_DATA], hd.variances()[:, GEN_DATA]
        real_b = hm.values()[:, GEN_TEMPLATE_CATEGORIES].sum(axis=1)
        nonsig = hm.values()[:, [k for k in GEN_MC_CATEGORIES if k != GEN_SIGNAL]].sum(axis=1)
        sig_only = hm.values()[:, GEN_SIGNAL]

        S = real_b[m] / real_b[m].sum()
        Btrue = nonsig[m] / nonsig[m].sum()
        e = np.sqrt(np.maximum(d[m], 1.0))
        x = (c[m] - c[m][0]) / (c[m][-1] - c[m][0])

        p_true, mod_true, chi_true = lsq([S, Btrue], d[m], e)
        best = None
        for b in np.linspace(-30.0, 30.0, 1201):
            comps = [S, np.exp(b * x)]
            pp, mm, cc = lsq(comps, d[m], e)
            if np.any(pp < 0):
                continue
            if best is None or cc < best[0]:
                best = (cc, pp, mm, b)
        chi_exp, p_exp, mod_exp, slope = best
        p_poly, mod_poly, chi_poly = lsq(
            [S, np.ones_like(x), x, x ** 2], d[m], e)

        ndf_true, ndf_exp, ndf_poly = m.sum() - 2, m.sum() - 3, m.sum() - 4
        out["background_models"] = {
            "true_shape": dict(chi2ndf=chi_true / ndf_true, signal=float(p_true[0] * S.sum())),
            "exponential": dict(chi2ndf=chi_exp / ndf_exp, signal=float(p_exp[0] * S.sum()),
                                slope=float(slope)),
            "quadratic": dict(chi2ndf=chi_poly / ndf_poly, signal=float(p_poly[0] * S.sum())),
            "bins": int(m.sum()), "data": float(d[m].sum()),
        }
        for k, v in out["background_models"].items():
            if isinstance(v, dict):
                print(f"  {k:14s} chi2/ndf {v['chi2ndf']:8.2f}   signal {v['signal']:9.0f}")

        # figure: the two models against the data
        ce = c[m]
        fig, (ax, rax) = plt.subplots(
            2, 1, figsize=(7.6, 5.4), sharex=True,
            gridspec_kw={"height_ratios": [3, 1.15], "hspace": 0.07})
        for a in (ax, rax):
            a.set_facecolor(SURFACE)
            a.grid(alpha=0.22)
            a.set_axisbelow(True)
        fig.patch.set_facecolor(SURFACE)
        ax.errorbar(ce, d[m], yerr=np.sqrt(dv[m]), fmt="o", ms=3.2, color=C_DATA,
                    lw=0, elinewidth=1.0, label="data")
        ax.plot(ce, mod_true, color=C_TRUE, lw=2.2,
                label=f"template + true simulated shape, $\\chi^2$/ndf {chi_true/ndf_true:.2f}")
        ax.plot(ce, mod_exp, color=C_EXP, lw=2.0, ls="--",
                label=f"template + exponential, $\\chi^2$/ndf {chi_exp/ndf_exp:.1f}")
        ax.set_ylabel("candidates / bin")
        ax.set_xlim(ce[0], ce[-1])
        ax.legend(fontsize=8.5, frameon=False)
        rax.axhline(1.0, color=INK_MUTED, lw=1.0)
        rax.plot(ce, d[m] / mod_true, "o", ms=3.0, color=C_TRUE)
        rax.plot(ce, d[m] / mod_exp, "s", ms=3.0, color=C_EXP)
        rax.set_ylabel("data / model", fontsize=9)
        rax.set_xlabel("$m(\\mu\\mu K)$, joint CVH [GeV]")
        rax.set_ylim(0.7, 1.35)
        fig.tight_layout()
        fig.savefig(f"{args.outdir}/fig_bkgmodel.png", bbox_inches="tight", dpi=140)
        plt.close(fig)

        # figure: the template against data with the FITTED TRUE background removed
        w = (c > WIDTH_WINDOW[0]) & (c < WIDTH_WINDOW[1])
        bkg_full = np.zeros_like(d)
        bkg_full[m] = p_true[1] * Btrue
        dsub = np.maximum(d - bkg_full, 0.0)
        sd, _ = gauss_width(c[w], dsub[w], B_MASS)
        ss, _ = gauss_width(c[w], sig_only[w], B_MASS)
        out["width_after_true_subtraction"] = dict(
            data_mev=sd * 1000.0, signal_mev=ss * 1000.0, ratio=sd / ss)
        print(f"  after true-shape subtraction: data {sd*1000:.2f} MeV, "
              f"template {ss*1000:.2f} MeV, ratio {sd/ss:.3f}")

        ee = edges[:-1][w]
        ee = np.append(ee, ee[-1] + (ee[-1] - ee[-2]))
        dn = dsub[w] / dsub[w].sum()
        dne = np.sqrt(dv[w]) / dsub[w].sum()
        sn = sig_only[w] / sig_only[w].sum()
        fig, ax = plt.subplots(figsize=(7.4, 4.4))
        fig.patch.set_facecolor(SURFACE)
        ax.set_facecolor(SURFACE)
        ax.stairs(sn, ee, color=C_MC, lw=2.0, baseline=None,
                  label=f"simulated B$^+\\to$J/$\\psi$K$^+$,  $\\sigma$ = {ss*1000:.1f} MeV")
        ax.errorbar(0.5 * (ee[:-1] + ee[1:]), dn, yerr=dne, fmt="o", ms=3.6,
                    color=C_DATA, lw=0, elinewidth=1.0,
                    label=f"data, true background removed,  $\\sigma$ = {sd*1000:.1f} MeV")
        ax.set_xlabel("$m(\\mu\\mu K)$, joint CVH [GeV]")
        ax.set_ylabel("fraction of the peak / bin")
        ax.set_xlim(ee[0], ee[-1])
        ax.set_ylim(0.0, None)
        ax.grid(alpha=0.22)
        ax.set_axisbelow(True)
        ax.legend(fontsize=9, frameon=False)
        ax.set_title(f"ratio {sd/ss:.2f}", fontsize=9, color=INK_MUTED, loc="right")
        fig.tight_layout()
        fig.savefig(f"{args.outdir}/fig_width.png", bbox_inches="tight", dpi=140)
        plt.close(fig)

    json.dump(out, open(f"{args.outdir}/widths.json", "w"), indent=1)
    print(f"wrote fig_bkgmodel.png, fig_width.png and widths.json to {args.outdir}")


if __name__ == "__main__":
    main()
