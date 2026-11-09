"""Read out the per-eta resolution study (openspec add-cvhnano-per-eta-resolution).

Covers only what rabbit's own plotting scripts cannot make. Prefit/postfit,
pulls, impacts and correlation matrices come from `rabbit_plot_hists.py`,
`rabbit_plot_pulls_and_impacts.py` and `rabbit_plot_cov.py`, run by the fit
driver. This script produces:

  * the fitted resolution parameter per eta bin, with its uncertainty, for every
    arm that carries `sigRes_eta{i}` nuisances -- on Asimov the central values
    are zero by construction and the UNCERTAINTY is the result: it says how well
    each eta bin's width is determined by that bin's own data;
  * the same over the (eta, pT) grid for the arm that has `sigRes_eta{i}_pt{j}`;
  * a cross-arm table: convergence, NLL, mean sigma(A/e/M), rho(A, e), the
    resolution uncertainties, and the largest |correlation| between any
    resolution nuisance and any scale nuisance, with the pair named;
  * the free-versus-constrained A/e/M comparison per eta group.

The resolution nuisance is log-interpolated, (smeared/nominal)^theta, with the
smeared template a convolution by the reference width sigma_ref. To first order
in sigma_ref^2 / w^2 a theta of t corresponds to an added width of
sigma_ref * sqrt(t) in quadrature (t < 0 means narrower than the simulation);
the implied width is quoted on that basis and marked as first-order.

Convergence is checked from edmval, two-sided, and a non-converged arm is
labelled as such in every table and plot.

Example:
    python scripts/btojpsik/cvhnano_resolution_report.py -o ~/public_html/mz/per_eta_resolution \\
        primary=.../cvhnano_..._reseta_..._fit.hdf5 constrained=.../..._res0005c_reseta_..._fit.hdf5 \\
        --free primary --constrained constrained
"""

import argparse
import os
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rabbit.io_tools as io

AEM = re.compile(r"^(?:scale)?([AeM])(\d+)$")
RES = re.compile(r"^sigRes(?:_eta(\d+))?(?:_pt(\d+))?(?:_q(\d+))?$")
REFW = re.compile(r"_res(\d{4,5})c?_")  # _res0005_ -> 0.005 GeV
EDM_MAX = 1e-5

# Categorical palette in fixed slot order (dataviz skill reference instance).
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
INK = "#1a1a19"
TITLES = True
PREFIX = ""
MUTED = "#6b6a63"


def read(path):
    fr = io.get_fitresult(path)
    h = fr["parms"].get()
    names = [str(x) for x in h.axes[0]]
    vals = np.asarray(h.values(), dtype=float)
    sig = np.sqrt(np.asarray(h.variances(), dtype=float))
    cov = np.asarray(fr["cov"].get().values()) if "cov" in fr else None
    edm = float(np.asarray(fr["edmval"]).reshape(-1)[0]) if "edmval" in fr else None
    nll = (
        float(np.asarray(fr["nllvalreduced"]).reshape(-1)[0])
        if "nllvalreduced" in fr
        else None
    )
    ndf = int(np.asarray(fr["ndfsat"]).reshape(-1)[0]) if "ndfsat" in fr else None
    m = REFW.search(os.path.basename(path))
    refw = int(m.group(1)) / 10 ** (len(m.group(1)) - 1) / 1.0 if m else None
    # "_res0005_" is written by tr -d '.' from 0.005 -> digits "0005" mean 0.005.
    if m:
        digits = m.group(1)
        refw = float("0." + digits[1:]) if digits.startswith("0") else float(digits)
    return dict(
        path=path, names=names, vals=vals, sig=sig, cov=cov, edm=edm, nll=nll,
        ndf=ndf, refw=refw,
        converged=(edm is not None and 0.0 <= edm <= EDM_MAX),
        index={n: i for i, n in enumerate(names)},
    )


def aem(d):
    out = {"A": {}, "e": {}, "M": {}}
    for n, s, v in zip(d["names"], d["sig"], d["vals"]):
        m = AEM.match(n)
        if m:
            out[m.group(1)][int(m.group(2))] = (v, s)
    return out


def rho_Ae(d):
    if d["cov"] is None:
        return np.nan
    ia = {int(AEM.match(n).group(2)): i for i, n in enumerate(d["names"])
          if AEM.match(n) and AEM.match(n).group(1) == "A"}
    ie = {int(AEM.match(n).group(2)): i for i, n in enumerate(d["names"])
          if AEM.match(n) and AEM.match(n).group(1) == "e"}
    c = d["cov"]
    vals = [c[ia[g], ie[g]] / np.sqrt(c[ia[g], ia[g]] * c[ie[g], ie[g]])
            for g in sorted(set(ia) & set(ie)) if c[ia[g], ia[g]] > 0 and c[ie[g], ie[g]] > 0]
    return float(np.mean(vals)) if vals else np.nan


def rho_Ae_dataonly(d):
    """Mean A-e correlation per eta group with the unit prior stripped.

    A and e are constrained nuisances with a unit Gaussian prior, so the
    posterior precision is P = I + F with F the data information. The posterior
    correlation read off the covariance is diluted toward zero whenever F is
    small compared to the prior. F = V^-1 - I on the 2x2 (A, e) marginal gives
    the data's own correlation, which is what should be compared across
    configurations with different events per nuisance group.
    """
    if d["cov"] is None:
        return np.nan
    ia = {int(AEM.match(n).group(2)): i for i, n in enumerate(d["names"])
          if AEM.match(n) and AEM.match(n).group(1) == "A"}
    ie = {int(AEM.match(n).group(2)): i for i, n in enumerate(d["names"])
          if AEM.match(n) and AEM.match(n).group(1) == "e"}
    c = d["cov"]
    vals = []
    for g in sorted(set(ia) & set(ie)):
        idx = [ia[g], ie[g]]
        V = c[np.ix_(idx, idx)]
        F = np.linalg.inv(V) - np.eye(2)
        if np.all(np.linalg.eigvals(F) > 0):
            Cd = np.linalg.inv(F)
            vals.append(Cd[0, 1] / np.sqrt(Cd[0, 0] * Cd[1, 1]))
    return float(np.mean(vals)) if vals else np.nan


def res_params(d):
    """{(eta, pt): (value, sigma, name)}; eta/pt None when the nuisance is global."""
    out = {}
    for n, v, s in zip(d["names"], d["vals"], d["sig"]):
        m = RES.match(n)
        if m:
            eta = int(m.group(1)) if m.group(1) is not None else None
            pt = int(m.group(2)) if m.group(2) is not None else None
            out[(eta, pt)] = (v, s, n)
    return out


def max_res_scale_corr(d):
    if d["cov"] is None:
        return (np.nan, "", "")
    c = d["cov"]
    ri = [i for i, n in enumerate(d["names"]) if RES.match(n)]
    si = [i for i, n in enumerate(d["names"]) if AEM.match(n)]
    if not ri or not si:
        return (np.nan, "", "")
    best = (0.0, "", "")
    for i in ri:
        for j in si:
            den = np.sqrt(c[i, i] * c[j, j])
            if den <= 0:
                continue
            r = c[i, j] / den
            if abs(r) > abs(best[0]):
                best = (float(r), d["names"][i], d["names"][j])
    return best


def eta_centres(nbins, lo=-2.4, hi=2.4):
    e = np.linspace(lo, hi, nbins + 1)
    return 0.5 * (e[:-1] + e[1:])


def implied_width_mev(theta, refw):
    if refw is None:
        return np.nan
    return np.sign(theta) * np.sqrt(abs(theta)) * refw * 1e3


def label_of(lab, d):
    return lab if d["converged"] else f"{lab} (NOT CONVERGED)"


def plot_eta(arms, labels, out, nbins):
    """Theta per eta bin with errors (top) and its uncertainty per bin (bottom)."""
    sel = [(l, arms[l]) for l in labels
           if any(k[0] is not None and k[1] is None for k in res_params(arms[l]))]
    if not sel:
        return None
    x = eta_centres(nbins)
    fig, (ax, ax2) = plt.subplots(
        2, 1, figsize=(12.0, 6.0), sharex=True,
        gridspec_kw=dict(height_ratios=[2.2, 1.0], hspace=0.06),
    )
    for k, (lab, d) in enumerate(sel):
        rp = res_params(d)
        th = np.array([rp.get((i, None), (np.nan, np.nan, ""))[0] for i in range(nbins)])
        er = np.array([rp.get((i, None), (np.nan, np.nan, ""))[1] for i in range(nbins)])
        off = (k - (len(sel) - 1) / 2) * 0.03
        col = PALETTE[k % len(PALETTE)]
        ax.errorbar(x + off, th, yerr=er, fmt="o", ms=4, lw=1.5, capsize=2,
                    color=col, label=label_of(lab, d))
        ax2.plot(x + off, er, "o-", ms=4, lw=1.5, color=col)
    ax.axhline(0.0, color=MUTED, lw=0.8, ls="--")
    ax.axhline(1.0, color=MUTED, lw=0.8, ls=":")
    refw = next((d["refw"] for _, d in sel if d["refw"]), None)
    ref_txt = f"; theta = 1 is a {refw*1e3:.0f} MeV smear in quadrature (first order)" if refw else ""
    ax.set_ylabel("fitted resolution parameter theta")
    if TITLES:
        ax.set_title(f"Per-eta resolution nuisance, fitted value and uncertainty{ref_txt}",
                     loc="left", fontsize=11, color=INK)
    ax.legend(frameon=False, fontsize=9, ncol=2)
    ax2.set_ylabel("sigma(theta)")
    ax2.set_xlabel("kaon eta (fit bin centre)")
    ax2.set_ylim(bottom=0)
    for a in (ax, ax2):
        a.grid(True, color="#e6e5df", lw=0.6)
        a.spines[["top", "right"]].set_visible(False)
    fn = os.path.join(out, PREFIX + "resolution_vs_eta.png")
    fig.savefig(fn, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return fn


def plot_etapt(arms, labels, out, nbins):
    sel = [(l, arms[l]) for l in labels
           if any(k[0] is not None and k[1] is not None for k in res_params(arms[l]))]
    if not sel:
        return None
    lab, d = sel[0]
    rp = res_params(d)
    npt = 1 + max(k[1] for k in rp if k[1] is not None)
    x = eta_centres(nbins)
    fig, (ax, ax2) = plt.subplots(
        2, 1, figsize=(12.0, 6.0), sharex=True,
        gridspec_kw=dict(height_ratios=[2.2, 1.0], hspace=0.06),
    )
    for j in range(npt):
        th = np.array([rp.get((i, j), (np.nan, np.nan, ""))[0] for i in range(nbins)])
        er = np.array([rp.get((i, j), (np.nan, np.nan, ""))[1] for i in range(nbins)])
        off = (j - (npt - 1) / 2) * 0.03
        col = PALETTE[j % len(PALETTE)]
        ax.errorbar(x + off, th, yerr=er, fmt="o", ms=4, lw=1.5, capsize=2, color=col,
                    label=f"kaon pT bin {j}")
        ax2.plot(x + off, er, "o-", ms=4, lw=1.5, color=col)
    ax.axhline(0.0, color=MUTED, lw=0.8, ls="--")
    if TITLES:
        ax.set_title(f"(eta, pT) resolution nuisances, arm '{label_of(lab, d)}'",
                     loc="left", fontsize=11, color=INK)
    ax.set_ylabel("fitted theta")
    ax.legend(frameon=False, fontsize=9, ncol=3)
    ax2.set_ylabel("sigma(theta)")
    ax2.set_xlabel("kaon eta (fit bin centre)")
    ax2.set_ylim(bottom=0)
    for a in (ax, ax2):
        a.grid(True, color="#e6e5df", lw=0.6)
        a.spines[["top", "right"]].set_visible(False)
    fn = os.path.join(out, PREFIX + "resolution_vs_eta_pt.png")
    fig.savefig(fn, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return fn


def plot_aem_per_group(arms, labels, out):
    """sigma(A), sigma(e), sigma(M) per eta group, one line per arm; ratio to the first."""
    fig, axes = plt.subplots(2, 3, figsize=(13.0, 6.0), sharex=True,
                             gridspec_kw=dict(height_ratios=[2.0, 1.0], hspace=0.08, wspace=0.25))
    base = aem(arms[labels[0]])
    for c, par in enumerate(("A", "e", "M")):
        ax, axr = axes[0, c], axes[1, c]
        g0 = sorted(base[par])
        b = np.array([base[par][g][1] for g in g0])
        for k, lab in enumerate(labels):
            d = arms[lab]
            gg = aem(d)[par]
            y = np.array([gg.get(g, (np.nan, np.nan))[1] for g in g0])
            col = PALETTE[k % len(PALETTE)]
            ax.plot(g0, y, "o-", ms=3.5, lw=1.4, color=col, label=label_of(lab, d))
            axr.plot(g0, y / b, "o-", ms=3.5, lw=1.4, color=col)
        ax.axhline(1.0, color=MUTED, lw=0.8, ls=":")
        ax.set_title(f"sigma({par}) per eta group, prior units", loc="left", fontsize=10, color=INK)
        axr.axhline(1.0, color=MUTED, lw=0.8, ls="--")
        axr.set_ylabel(f"ratio to {labels[0]}")
        axr.set_xlabel("eta group")
        for a in (ax, axr):
            a.grid(True, color="#e6e5df", lw=0.6)
            a.spines[["top", "right"]].set_visible(False)
    axes[0, 0].legend(frameon=False, fontsize=8)
    fn = os.path.join(out, PREFIX + "aem_sigma_per_group.png")
    fig.savefig(fn, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return fn


def plot_aem_pulls(arms, labels, out):
    """A, e, M fitted values with errors per eta group, for arms whose values
    are not identically zero (i.e. observed fits); Asimov arms are skipped."""
    sel = [(l, arms[l]) for l in labels
           if any(abs(v) > 1e-9 for n, v in zip(arms[l]["names"], arms[l]["vals"]) if AEM.match(n))]
    if not sel:
        return None
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.6), sharey=False)
    for c, par in enumerate(("A", "e", "M")):
        ax = axes[c]
        for k, (lab, d) in enumerate(sel):
            g = aem(d)[par]
            gs = sorted(g)
            off = (k - (len(sel) - 1) / 2) * 0.15
            ax.errorbar([x + off for x in gs], [g[x][0] for x in gs], yerr=[g[x][1] for x in gs],
                        fmt="o", ms=4, lw=1.5, capsize=2, color=PALETTE[k % len(PALETTE)],
                        label=label_of(lab, d))
        ax.axhline(0.0, color=MUTED, lw=0.8, ls="--")
        ax.axhspan(-1, 1, color="#e6e5df", alpha=0.5, lw=0)
        ax.set_title(f"scale{par}: fitted value, prior units", loc="left", fontsize=10, color=INK)
        ax.set_xlabel("eta group")
        ax.grid(True, color="#e6e5df", lw=0.6)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].legend(frameon=False, fontsize=8)
    fn = os.path.join(out, PREFIX + "aem_pulls_per_group.png")
    fig.savefig(fn, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return fn


def table(arms, labels):
    rows = []
    hdr = ["arm", "conv", "edmval", "NLL(red)", "ndf", "mean sA", "mean se", "mean sM",
           "rho(A,e)", "rho(A,e) data-only", "nRes", "mean s(res)", "max s(res)", "max|rho(res,scale)|", "pair"]
    for lab in labels:
        d = arms[lab]
        a = aem(d)
        rp = res_params(d)
        rs = np.array([v[1] for v in rp.values()]) if rp else np.array([np.nan])
        r, n1, n2 = max_res_scale_corr(d)
        rows.append([
            lab, "yes" if d["converged"] else "NO",
            f"{d['edm']:.2e}" if d["edm"] is not None else "n/a",
            f"{d['nll']:.2f}" if d["nll"] is not None else "n/a",
            str(d["ndf"]) if d["ndf"] is not None else "n/a",
            *(f"{np.mean([v[1] for v in a[p].values()]):.4f}" if a[p] else "n/a" for p in "AeM"),
            f"{rho_Ae(d):.4f}", f"{rho_Ae_dataonly(d):.3f}", str(len(rp)),
            f"{np.nanmean(rs):.4f}" if rp else "-", f"{np.nanmax(rs):.4f}" if rp else "-",
            f"{r:.3f}" if rp else "-", f"{n1} / {n2}" if rp else "-",
        ])
    w = [max(len(h), *(len(r[i]) for r in rows)) for i, h in enumerate(hdr)]
    lines = ["| " + " | ".join(h.ljust(w[i]) for i, h in enumerate(hdr)) + " |",
             "|" + "|".join("-" * (x + 2) for x in w) + "|"]
    for r in rows:
        lines.append("| " + " | ".join(r[i].ljust(w[i]) for i in range(len(hdr))) + " |")
    return "\n".join(lines)


def free_vs_constrained(arms, free, con):
    df, dc = arms[free], arms[con]
    af, ac = aem(df), aem(dc)
    lines = [f"Free ('{free}') against constrained ('{con}'): per eta group, "
             "(theta_free - theta_con)/sigma_free and sigma_con/sigma_free",
             "| group | A shift | A sig ratio | e shift | e sig ratio | M shift | M sig ratio |",
             "|---|---|---|---|---|---|---|"]
    groups = sorted(set(af["A"]) & set(ac["A"]))
    worst = 0.0
    for g in groups:
        cells = []
        for p in "AeM":
            vf, sf = af[p][g]
            vc, sc = ac[p][g]
            sh = (vf - vc) / sf if sf > 0 else np.nan
            worst = max(worst, abs(sh) if np.isfinite(sh) else 0.0)
            cells += [f"{sh:+.3f}", f"{sc / sf:.4f}" if sf > 0 else "n/a"]
        lines.append(f"| {g} | " + " | ".join(cells) + " |")
    lines.append(f"\nLargest |shift| in units of the free-arm sigma: {worst:.3f}")
    return "\n".join(lines)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("inputs", nargs="+", help="label=fitresult.hdf5 ...")
    p.add_argument("-o", "--outfolder", required=True)
    p.add_argument("--etaBins", type=int, default=24, help="fit eta bins (48/etaRebin)")
    p.add_argument("--free", default=None, help="label of the free-resolution arm")
    p.add_argument("--constrained", default=None, help="label of the constrained arm")
    p.add_argument("--title", default="per_eta_resolution")
    p.add_argument("--plotArms", nargs="*", default=None,
                   help="labels to draw in the figures (default: all, at most 6 -- "
                   "the categorical palette has six hues and is never cycled)")
    p.add_argument("--noTitles", action="store_true",
                   help="omit in-figure titles (the slide title carries them)")
    p.add_argument("--figPrefix", default="", help="prefix for figure filenames")
    a = p.parse_args()
    os.makedirs(a.outfolder, exist_ok=True)

    labels, arms = [], {}
    for spec in a.inputs:
        lab, path = spec.split("=", 1)
        try:
            arms[lab] = read(path)
            labels.append(lab)
        except Exception as exc:
            print(f"{lab}: could not read {path} ({type(exc).__name__}: {exc})")
    if not labels:
        return 1

    md = [f"# {a.title}\n", table(arms, labels), ""]
    for lab in labels:
        d = arms[lab]
        rp = res_params(d)
        if rp and any(k[0] is not None for k in rp):
            md.append(f"\n## Resolution parameters, arm '{lab}'"
                      + ("" if d["converged"] else " (NOT CONVERGED)"))
            md.append("| nuisance | theta | sigma | implied added width [MeV] (first order) |")
            md.append("|---|---|---|---|")
            for key in sorted(rp, key=lambda k: (k[1] if k[1] is not None else -1, k[0])):
                v, s, n = rp[key]
                md.append(f"| {n} | {v:+.4f} | {s:.4f} | {implied_width_mev(v, d['refw']):+.2f} |")
    if a.free and a.constrained and a.free in arms and a.constrained in arms:
        md.append("\n## Free versus constrained\n")
        md.append(free_vs_constrained(arms, a.free, a.constrained))

    plot_labels = [l for l in (a.plotArms or labels) if l in arms]
    if len(plot_labels) > len(PALETTE):
        raise SystemExit(
            f"{len(plot_labels)} arms to draw but the palette has {len(PALETTE)} "
            "hues; pass --plotArms with at most six labels (hues are not cycled)"
        )
    global TITLES, PREFIX
    TITLES = not a.noTitles
    PREFIX = a.figPrefix
    figs = [f for f in (plot_eta(arms, plot_labels, a.outfolder, a.etaBins),
                        plot_etapt(arms, plot_labels, a.outfolder, a.etaBins),
                        plot_aem_per_group(arms, plot_labels, a.outfolder),
                        plot_aem_pulls(arms, plot_labels, a.outfolder)) if f]
    md.append("\n## Figures\n" + "\n".join(f"- {os.path.basename(f)}" for f in figs))
    text = "\n".join(md)
    with open(os.path.join(a.outfolder, "report.md"), "w") as f:
        f.write(text + "\n")
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
