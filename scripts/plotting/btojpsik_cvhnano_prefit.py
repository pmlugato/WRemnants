"""The data-side prefit for the CVH NanoAOD B+ -> J/psi K channel.

This is the second of two plots per mass variable, and the division of labour is
deliberate.

The simulation is an inclusive b -> J/psi X production. Every component it can
offer is a **real b-hadron decay**; it contains none of the prompt-charmonium or
fake-muon combinatorial that data also has, and the combinatorial it does
contain is 87% of the peak window. So the simulated stack is not a model of the
data, and normalising it to the data yield would inflate the signal by roughly
eight to cover background the sample does not have.

  * The **simulation-only stack** is therefore produced by the central
    `makeDataMCStackPlot.py --noData`, with no normalisation to data. Nothing
    here duplicates it.
  * This script draws the **data**, with the signal template plus an analytic
    background fitted to it, so the signal normalisation is *measured on the
    data* -- which is what a prefit plot is supposed to show, and what the
    central tool cannot do.

Only mass variables get one, because a signal peak and a smooth background are
separable nowhere else.

Groups, labels and colours come from the central machinery: `datagroups_cvhnano`
slices the `genCategory` axis and `styles` names the components.

    python scripts/plotting/btojpsik_cvhnano_prefit.py <histmaker.hdf5> \
        -o /home/user/public_html/mz --outfolder cvhnano_prefit
"""

import argparse
import os

import matplotlib as mpl
import numpy as np

mpl.use("Agg")
import matplotlib.pyplot as plt

from wremnants.postprocessing.datagroups.datagroups import Datagroups
from wremnants.postprocessing.datagroups.datagroups_cvhnano import TEMPLATE_GROUP_NAMES
from wums import logging

logger = logging.child_logger(__name__)

INK_MUTED, SURFACE = "#6b6b6b", "#fcfcfb"
C_DATA, C_BKG, C_TOT = "black", "#6b6b6b", "#c0392b"

# variable -> (histogram name, axis name, axis label, is a mass)
VARIABLES = [
    (
        "nominal_BuJpsiK_jointCvhMass",
        "jointCvh_mass",
        r"$m(\mu\mu K)$, joint CVH [GeV]",
        True,
    ),
    (
        "nominal_BuJpsiK_kvfCvhMass",
        "kvfCvh_mass",
        r"$m(\mu\mu K)$, KVF on CVH tracks [GeV]",
        True,
    ),
    (
        "nominal_BuJpsiK_kvfRawMass",
        "kvfRaw_mass",
        r"$m(\mu\mu K)$, KVF on raw tracks [GeV]",
        True,
    ),
    ("nominal_BuJpsiK_mass", "cand_mass", r"$m(\mu\mu K)$, pre-fit [GeV]", True),
    (
        "nominal_BuJpsiK_kvfCvhDimuonMass",
        "kvfCvh_dimuon_mass",
        r"$m(\mu\mu)$, unconstrained, CVH tracks [GeV]",
        False,
    ),
    (
        "nominal_BuJpsiK_jointCvhPt",
        "jointCvh_pt",
        r"$p_{T}(B)$, joint CVH [GeV]",
        False,
    ),
    ("nominal_BuJpsiK_kaonPt", "bachelor_pt", r"$p_{T}(K)$ [GeV]", False),
    ("nominal_BuJpsiK_kaonEta", "bachelor_eta", r"$\eta(K)$", False),
    (
        "nominal_BuJpsiK_jointCvhJpsiPt",
        "jointCvh_jpsi_pt",
        r"$p_{T}(J/\psi)$, joint CVH [GeV]",
        False,
    ),
    ("nominal_cand_mu0Pt", "mu0_pt", r"$p_{T}(\mu_{0})$ [GeV]", False),
    ("nominal_cand_bachelorD0", "bachelor_d0", r"$d_{0}(K)$ w.r.t. PV [cm]", False),
    (
        "nominal_cand_bachelorD0Sig",
        "bachelor_d0Sig",
        r"$d_{0}(K)/\sigma_{d_{0}}$",
        False,
    ),
    (
        "nominal_cand_bachelorNValidHits",
        "bachelor_nValidHits",
        r"$N_{\mathrm{hits}}(K)$",
        False,
    ),
    (
        "nominal_BuJpsiK_jointCvhSl3d",
        "jointCvh_sl3d",
        r"$L_{3D}/\sigma$ ($B$), joint CVH",
        False,
    ),
    (
        "nominal_BuJpsiK_jointCvhAlpha3dPV",
        "jointCvh_alpha3dPV",
        r"$\alpha_{3D}$ ($B$, PV) [rad]",
        False,
    ),
    ("nominal_nPV", "nPV", r"$N_{PV}$", False),
    ("nominal_Pileup_nTrueInt", "nTrueInt", r"$N_{\mathrm{true\ int}}$", False),
]

# Fit window for panel B. Wider than the eventual fit window on purpose: the
# background has to be seen turning over on both sides of the peak, or an
# exponential describing it cannot be judged.
FIT_WINDOW = (5.15, 5.42)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("infile", help="Histmaker output")
    p.add_argument("-o", "--outpath", default=os.path.expanduser("~/public_html/mz"))
    p.add_argument("--outfolder", default="cvhnano_prefit")
    p.add_argument("--rebin", type=int, default=1)
    p.add_argument(
        "--vars", nargs="*", default=None, help="Subset of histogram names to plot"
    )
    return p.parse_args()


def lstsq_amplitudes(components, target, sigma):
    """Weighted linear least squares for the amplitudes of fixed shapes.

    Closed form on purpose. `scipy.optimize.least_squares` in this build
    intermittently raises "Memory error in scipy.linalg.svd" on systems of this
    size, and there is no reason to run an optimiser on a problem that has an
    exact solution.
    """
    A = np.vstack([c / sigma for c in components]).T
    return np.linalg.lstsq(A, target / sigma, rcond=None)[0]


def fit_signal_plus_exp(templates, data, sigma, x):
    """Signal template(s) plus an exponential background.

    Conditionally linear: for a fixed slope every amplitude solves exactly, so
    the slope is scanned on a grid and then refined. Deterministic, and it cannot
    fail to converge -- which matters here, because this is a diagnostic plot and
    a plot that sometimes has no curve on it is worse than no plot.
    """
    best = None
    for grid in (np.linspace(-25.0, 25.0, 501),):
        for b in grid:
            comps = list(templates) + [np.exp(b * x)]
            p = lstsq_amplitudes(comps, data, sigma)
            if np.any(p < 0):
                continue
            model = sum(a * c for a, c in zip(p, comps))
            chi2 = float(np.sum(((model - data) / sigma) ** 2))
            if best is None or chi2 < best[0]:
                best = (chi2, b, p)
    if best is None:
        return None
    _, b0, _ = best
    for b in np.linspace(b0 - 0.1, b0 + 0.1, 41):
        comps = list(templates) + [np.exp(b * x)]
        p = lstsq_amplitudes(comps, data, sigma)
        if np.any(p < 0):
            continue
        model = sum(a * c for a, c in zip(p, comps))
        chi2 = float(np.sum(((model - data) / sigma) ** 2))
        if chi2 < best[0]:
            best = (chi2, b, p)
    return best


def project(h, axis, rebin=1):
    hp = h.project(axis)
    if rebin > 1:
        hp = hp[{axis: hist_rebin(rebin)}]
    v = hp.values()
    var = hp.variances()
    return hp.axes[0].edges, v, (v if var is None else var)


def hist_rebin(n):
    import hist

    return hist.rebin(n)


def panel_b(groups_hists, axis, label, outfile, rebin):
    """Data, with the signal template plus an exponential fitted to it."""
    if "Data" not in groups_hists or "BuJpsiK" not in groups_hists:
        return None
    edges, d, dvar = project(groups_hists["Data"], axis, rebin)
    _, s, _ = project(groups_hists["BuJpsiK"], axis, rebin)
    if d.sum() <= 0 or s.sum() <= 0:
        return None

    centres = 0.5 * (edges[:-1] + edges[1:])
    m = (centres >= FIT_WINDOW[0]) & (centres <= FIT_WINDOW[1])
    if m.sum() < 8:
        return None

    # The other real-B categories are NOT free parameters. Their yields relative
    # to the signal are fixed by the simulation, so they enter as one shape
    # scaled with the signal amplitude. Letting a 0.5% component float would
    # give the fit a second peak-like degree of freedom it has not earned.
    # TEMPLATE_GROUP_NAMES, not "everything that is not Data" -- the datagroups
    # also define a NoBAncestor diagnostic group, and folding the combinatorial
    # into the signal template would be silent and catastrophic.
    other = np.zeros_like(s)
    for name in TEMPLATE_GROUP_NAMES:
        if name == "BuJpsiK" or name not in groups_hists:
            continue
        other += project(groups_hists[name], axis, rebin)[1]
    ratio_other = other.sum() / s.sum() if s.sum() > 0 else 0.0
    norm = s.sum() + other.sum()
    sig_shape = (s + other) / norm
    # Fraction of the fitted template that is the signal decay itself, so the
    # quoted yield is the signal and not signal-plus-everything-else.
    denom = s[m].sum() + other[m].sum()
    if denom <= 0:
        logger.warning("no template entries inside the fit window for %s", axis)
        return None
    signal_fraction = float(s[m].sum() / denom)

    x = (centres[m] - centres[m][0]) / (centres[m][-1] - centres[m][0])
    sigma = np.sqrt(np.maximum(d[m], 1.0))
    fit = fit_signal_plus_exp([sig_shape[m]], d[m], sigma, x)
    if fit is None:
        logger.warning("prefit did not find a positive solution for %s", axis)
        return None
    chi2, slope, pars = fit
    a_sig, a_bkg = float(pars[0]), float(pars[1])
    bkg = a_bkg * np.exp(slope * x)
    sig = a_sig * sig_shape[m]
    ndf = int(m.sum()) - 3

    fig, (ax, rax) = plt.subplots(
        2,
        1,
        figsize=(7.0, 5.6),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.07},
    )
    for a in (ax, rax):
        a.set_facecolor(SURFACE)
        a.grid(alpha=0.22)
        a.set_axisbelow(True)
    fig.patch.set_facecolor(SURFACE)

    ax.errorbar(
        centres,
        d,
        yerr=np.sqrt(dvar),
        fmt="o",
        ms=3.2,
        color=C_DATA,
        lw=0,
        elinewidth=1.0,
        label="data",
    )
    ax.plot(
        centres[m],
        bkg,
        color=C_BKG,
        lw=1.6,
        ls="--",
        label="fitted background (exponential)",
    )
    ax.plot(centres[m], bkg + sig, color=C_TOT, lw=2.0, label="templates + background")
    ax.fill_between(
        centres[m],
        bkg,
        bkg + sig,
        color=C_TOT,
        alpha=0.20,
        lw=0,
        label="simulated b-hadron templates",
    )
    ax.set_ylabel("candidates / bin")
    ax.set_xlim(edges[0], edges[-1])
    ax.set_ylim(0.0, None)
    ax.legend(fontsize=8.5, frameon=False)
    n_template = a_sig * sig_shape[m].sum()
    ax.set_title(
        f"fitted B$^+\\to$J/$\\psi$K$^+$ {n_template * signal_fraction:.0f} "
        f"in window, $\\chi^2$/ndf {chi2 / max(ndf, 1):.2f}",
        fontsize=9,
        color=INK_MUTED,
        loc="right",
    )

    model = bkg + sig
    rax.axhline(1.0, color=C_TOT, lw=1.2)
    rax.errorbar(
        centres[m],
        d[m] / model,
        yerr=np.sqrt(dvar[m]) / model,
        fmt="o",
        ms=3.2,
        color=C_DATA,
        lw=0,
        elinewidth=1.0,
    )
    rax.set_ylabel("data / model", fontsize=9)
    rax.set_xlabel(label)
    lim = max(0.12, float(np.nanmax(np.abs(d[m] / model - 1.0))) * 1.25)
    rax.set_ylim(1 - lim, 1 + lim)

    fig.tight_layout()
    fig.savefig(outfile, bbox_inches="tight", dpi=140)
    plt.close(fig)
    return dict(
        template_yield=float(n_template),
        signal=float(n_template * signal_fraction),
        signal_fraction_of_template=signal_fraction,
        background=float(bkg.sum()),
        data=float(d[m].sum()),
        chi2ndf=float(chi2 / max(ndf, 1)),
        slope=float(slope),
        other_over_signal=float(ratio_other),
    )


def main():
    args = parse_args()
    outdir = os.path.join(args.outpath, args.outfolder)
    os.makedirs(outdir, exist_ok=True)

    groups = Datagroups(args.infile)
    logger.info("mode %s, groups %s", groups.mode, list(groups.groups))

    wanted = [v for v in VARIABLES if v[3]]
    if args.vars:
        wanted = [v for v in VARIABLES if v[0] in args.vars]

    summary = {}
    for hname, axis, label, is_mass in wanted:
        try:
            groups.setNominalName(hname)
            groups.loadHistsForDatagroups(hname, syst="", label="nominal")
        except Exception as exc:
            logger.warning("skipping %s: %s", hname, exc)
            continue
        hists = {}
        for name, group in groups.groups.items():
            h = group.hists.get("nominal")
            if h is None:
                continue
            if axis not in [a.name for a in h.axes]:
                continue
            hists[name] = h
        if not hists:
            logger.warning("no group carries axis %s for %s", axis, hname)
            continue

        base = hname.replace("nominal_", "")
        yields = {
            name: float(project(h, axis, args.rebin)[1].sum())
            for name, h in hists.items()
        }
        fitres = panel_b(
            hists, axis, label, os.path.join(outdir, f"{base}_prefit.png"), args.rebin
        )
        summary[base] = dict(yields=yields, prefit=fitres)
        logger.info("%s: %s", base, fitres if fitres else "yields only")

    import json

    with open(os.path.join(outdir, "prefit_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    print(f"wrote {len(summary)} variables to {outdir}")


if __name__ == "__main__":
    main()
