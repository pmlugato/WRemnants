#!/usr/bin/env python3
"""Figures for the A/e/M unit-fix deck (notes/aem_unit_mismatch.md, section 9).

Every number is read from its source: the histograms' recorded arguments (via
the joint-fit guard), the closure outputs in public_html, and the fit results.
Run inside the container:

    python scripts/btojpsik/make_aem_unit_fix_figs.py
"""

import importlib.util
import json
import os
import re
import sys

import numpy as np

sys.path.insert(0, "rabbit")
sys.path.insert(0, "wums")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from rabbit import io_tools

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "..", "slides", "aem_unit_fix", "assets")
CEPH = "/ceph/submit/data/user/p/pmlugato/mz"
H, F = f"{CEPH}/histograms", f"{CEPH}/fitresults"
PUB = os.path.expanduser("~/public_html/mz/combined_fits_and_mB_prefit_scale_uncs_matched")

INK, INK2, GRID = "#1F1F1F", "#595959", "#e6e6e3"
REF, OLD, NEW = "#8a8985", "#eb6834", "#2a78d6"
plt.rcParams.update({"font.size": 12, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.titlesize": 13})


def style(ax):
    ax.grid(color=GRID)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    fig.savefig(os.path.join(OUT, name), dpi=170, bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


def load_guard():
    spec = importlib.util.spec_from_file_location(
        "cjt", os.path.join(HERE, "..", "rabbit", "check_joint_tensor.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


# ---------------------------------------------------------------- 1. units
def fig_units():
    g = load_guard()
    rows = {
        "J/ψ full (mB fits)": ("jpsi", f"{H}/dimuon_resonances_calinput_jpsi_2016PostVFP_etaBins_24_joint_full.hdf5", REF, "-"),
        "Z, J/ψ repro (mZ fits)": ("z", f"{H}/mz_dilepton_scetlib_dyturbo_CT18Z_N3p0LL_N2LO_Corr_joint_repro.hdf5", REF, "--"),
        "B⁺ as built (massfit widths)": ("bplus", f"{H}/btojpsik_cvhnano_kpiA_kpi.hdf5", OLD, "-"),
    }
    eta = np.linspace(-2.4, 2.4, 49)
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.4))
    for label, (kind, path, color, ls) in rows.items():
        c = g.Checker()
        w = g.implied_widths(c, label, kind, g.read_meta_args(path))
        for ax, i in zip(axs, (0, 1)):
            ax.stairs(w[:, i], eta, color=color, lw=2.2, ls=ls, label=label)
    for ax, p, lim in zip(axs, ("A (fractional scale)", "e [GeV]"), ((5e-5, 1.5e-2), (5e-4, 1.5e-1))):
        ax.set_yscale("log")
        ax.set_ylim(*lim)
        ax.set_xlabel("η (48 calibration bins)")
        ax.set_title(f"one unit of {p}")
        ax.set_xlim(-2.4, 2.4)
        style(ax)
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, frameon=False, fontsize=11, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.07))
    save(fig, "units_per_eta.png")


# ------------------------------------------------- 2. template response per unit
def parse_shift(path):
    out, tag = {}, None
    for line in open(path):
        m = re.match(r"=== (.*)", line)
        if m:
            tag = m.group(1).strip()
            out[tag] = {}
            continue
        m = re.match(r"\s+(\w): pT bin (\d): induced mean shift per unit = ([-+.\d]+)", line)
        if m and tag:
            out[tag].setdefault(m.group(1), [0, 0, 0])[int(m.group(2))] = float(m.group(3))
    return out


def fig_slopes():
    s = parse_shift(f"{PUB}/02_closure/template_induced_slope_theta0p05.txt")
    old = s["old unit, mB chain (realmc)"]
    new = s["matched 1/1/1 stretched A,e x5 (mB chain)"]
    x = np.arange(3)
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.4))
    for ax, p, title in zip(axs, "Ae", ("A", "e")):
        vo, vn = np.abs(old[p]), np.abs(new[p])
        ax.bar(x - 0.19, vo, 0.36, color=OLD, edgecolor="white", lw=2, label="as built (old unit)")
        ax.bar(x + 0.19, vn, 0.36, color=NEW, edgecolor="white", lw=2, label="matched unit")
        for xi, a, b in zip(x, vo, vn):
            ax.text(xi - 0.19, a * 1.12, f"{a:.2f}", ha="center", fontsize=10.5, color=INK)
            ax.text(xi + 0.19, b * 1.12, f"{b:.1f}", ha="center", fontsize=10.5, color=INK)
        ax.set_yscale("log")
        ax.set_ylim(0.1, 80)
        ax.set_xticks(x, ["kaon pT bin 1\n(softest)", "bin 2", "bin 3"])
        ax.set_title(f"|B⁺ mass shift| per unit of {title} [MeV]")
        style(ax)
    axs[0].legend(frameon=False, fontsize=10.5, loc="upper left")
    save(fig, "template_slope_per_unit.png")


# ------------------------------------------------------ 3. reweight linearity
def fig_linearity():
    j = json.load(open(f"{PUB}/02_closure/closure_summary.json"))["linearity"]
    bands = ["1–2", "2–2.8", "2.8–5.7", "> 5.7"]
    x = np.arange(4)
    fig, ax = plt.subplots(figsize=(8, 4.2))
    for p, color, mk in (("A", INK2, "s"), ("e", NEW, "o")):
        v = [j[f"{p}_{i}"]["weighted_mean"] for i in range(4)]
        lab = f"A  (0.99 in every band)" if p == "A" else "e"
        ax.plot(x, v, marker=mk, ms=9, lw=2, color=color, label=lab)
        if p == "e":
            for xi, vi in zip(x[:3], v[:3]):
                ax.text(xi + 0.08, vi - 0.035, f"{vi:.2f}", fontsize=11, color=INK)
    ax.axhline(1, color=INK2, lw=1, ls="--")
    ax.set_xticks(x, bands)
    ax.set_xlabel("kaon pT [GeV]")
    ax.set_ylabel("built at 5/5/1, relative to linear")
    ax.set_ylim(0.55, 1.08)
    ax.legend(frameon=False, loc="lower right")
    style(ax)
    save(fig, "reweight_linearity.png")


# ----------------------------------------------------------- fit summaries
def summarise(path, poi, scale, n_eta=24):
    r = io_tools.get_fitresult(path)
    names = np.array([str(x) for x in r["parms"].get().axes["parms"]])
    cov = np.asarray(r["cov"].get().values())
    idx = {p: np.array([int(np.flatnonzero(names == f"{p}{i}")[0]) for i in range(n_eta)]) for p in "AeM"}
    sig = np.sqrt(np.diag(cov))
    out = dict(sigA=sig[idx["A"]].mean(), sige=sig[idx["e"]].mean())
    hit = np.flatnonzero(names == poi)
    if hit.size:
        ip = int(hit[0])
        fix = np.r_[idx["A"], idx["e"], idx["M"]]
        c = cov[ip, fix]
        cond = np.sqrt(cov[ip, ip] - c @ np.linalg.solve(cov[np.ix_(fix, fix)], c))
        out.update(sigma=scale * sig[ip], floor=scale * cond,
                   impact=scale * np.sqrt(max(sig[ip] ** 2 - cond ** 2, 0)))
    return out


def hbars(ax, labels, vals, colors, fmt, hatches=None, floors=None):
    y = np.arange(len(vals))[::-1]
    for i, (yi, v, c) in enumerate(zip(y, vals, colors)):
        ax.barh(yi, v, 0.62, color=c, edgecolor="white", lw=2,
                hatch=(hatches[i] if hatches else None))
        ax.text(v, yi, " " + fmt.format(v), va="center", fontsize=11, color=INK)
        if floors is not None and floors[i] is not None:
            ax.plot([floors[i]] * 2, [yi - 0.36, yi + 0.36], color=INK, lw=2.4)
    ax.set_yticks(y, labels)
    ax.set_xlim(0, max(vals) * 1.32)
    style(ax)


def fig_mB():
    f = {k: summarise(p, "massShiftB5MeV", 5) for k, p in {
        "old": f"{F}/realmc_A_data_fit.hdf5",
        "new": f"{F}/joint_mz_jpsi_bplus_matched_mB_noZ_data_jpsimu_jpsifull_fit.hdf5",
        "jpsi": f"{F}/realmc_D_jpsionly_data_fit.hdf5"}.items()}
    lab = ["old unit", "matched"]
    fig, axs = plt.subplots(1, 2, figsize=(12, 3.4))
    hbars(axs[0], lab, [f["old"]["sigma"], f["new"]["sigma"]], [OLD, NEW], "{:.4f}",
          floors=[f["old"]["floor"], f["new"]["floor"]])
    axs[0].set_title("σ(mB) [MeV]   (black tick: A, e, M fixed)")
    hbars(axs[1], lab, [f["old"]["impact"], f["new"]["impact"]], [OLD, NEW], "{:.4f}")
    axs[1].set_title("impact of A+e+M on mB [MeV]")
    axs[1].set_yticklabels([])
    save(fig, "mB_before_after.png")

    fig, axs = plt.subplots(1, 2, figsize=(12, 3.6))
    lab = ["J/ψ alone", "J/ψ + B⁺, old unit", "J/ψ + B⁺, matched"]
    cols = [REF, OLD, NEW]
    hbars(axs[0], lab, [f[k]["sigA"] for k in ("jpsi", "old", "new")], cols, "{:.3f}")
    axs[0].set_title("mean σ(A_i)  [units of 5×10⁻³]")
    hbars(axs[1], lab, [f[k]["sige"] for k in ("jpsi", "old", "new")], cols, "{:.3f}")
    axs[1].set_title("mean σ(e_i)  [units of 5×10⁻² GeV]")
    axs[1].set_yticklabels([])
    save(fig, "aem_sigma_mB_chain.png")


def fig_mZ():
    f = {k: summarise(p, "massShiftZ100MeV", 100) for k, p in {
        "ref": f"{F}/joint_mz_jpsi_bplus_res5_260914_noB_fit.hdf5",
        "old": f"{F}/joint_mz_jpsi_bplus_res5_260914_fit.hdf5",
        "new": f"{F}/joint_mz_jpsi_bplus_matched_mz_fit.hdf5"}.items()}
    lab = ["Z + J/ψ", "Z + J/ψ + B⁺, old unit", "Z + J/ψ + B⁺, matched"]
    cols = [REF, OLD, NEW]
    fig, axs = plt.subplots(1, 2, figsize=(12, 3.6))
    hbars(axs[0], lab, [f[k]["sigma"] for k in ("ref", "old", "new")], cols, "{:.3f}",
          floors=[f[k]["floor"] for k in ("ref", "old", "new")])
    axs[0].set_title("σ(mZ) [MeV]   (black tick: A, e, M fixed)")
    hbars(axs[1], lab, [f[k]["impact"] for k in ("ref", "old", "new")], cols, "{:.3f}")
    axs[1].set_title("impact of A+e+M on mZ [MeV]")
    axs[1].set_yticklabels([])
    save(fig, "mZ_before_after.png")


def fig_postfit():
    fits = [("J/ψ alone", REF, f"{F}/realmc_D_jpsionly_data_fit.hdf5"),
            ("J/ψ + B⁺, old unit", OLD, f"{F}/realmc_A_data_fit.hdf5"),
            ("J/ψ + B⁺, matched", NEW, f"{F}/joint_mz_jpsi_bplus_matched_mB_noZ_data_jpsimu_jpsifull_fit.hdf5")]
    eta = np.linspace(-2.4, 2.4, 25)
    etac = 0.5 * (eta[1:] + eta[:-1])
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.6), sharey=True)
    for k, (label, color, path) in enumerate(fits):
        h = io_tools.get_fitresult(path)["parms"].get()
        names = np.array([str(x) for x in h.axes["parms"]])
        v, s = h.values(), np.sqrt(h.variances())
        for ax, p in zip(axs, "Ae"):
            ii = [int(np.flatnonzero(names == f"{p}{i}")[0]) for i in range(24)]
            ax.errorbar(etac + (k - 1) * 0.045, v[ii], s[ii], fmt="o", ms=5, lw=1.5, color=color, label=label)
    for ax, t in zip(axs, ("postfit A_i [units of 5×10⁻³]", "postfit e_i [units of 5×10⁻² GeV]")):
        ax.axhline(0, color=INK2, lw=0.8)
        ax.set_title(t)
        ax.set_xlabel("η group")
        style(ax)
    axs[0].legend(frameon=False, fontsize=10.5, loc="upper center")
    save(fig, "postfit_Ae_data.png")


# ------------------------------------------------- 7. the real B+ templates
def tensor_cell(fname, syst, ipt=0, ieta=10, theta=1.0):
    import h5py, hdf5plugin  # noqa: F401
    h = h5py.File(f"{CEPH}/tensors/{fname}", "r")
    systs = [x.decode() for x in h["hsysts"][...]]
    procs = [x.decode() for x in h["hprocs"][...]]
    isig = procs.index("BuJpsiK")
    norm = np.asarray(h["hnorm"][...]).reshape(-1, len(procs))
    logk = np.asarray(h["hlogk"][...]).reshape(norm.shape[0], len(procs), -1)
    js = systs.index(syst)
    asym = logk.shape[2] == 2 * len(systs)
    shape = (3, 24, 2, 20)
    nom = norm[:, isig].reshape(shape)[ipt, ieta].sum(axis=0)
    if asym:   # rabbit stores (avg, halfdiff) pairs for asymmetric tensors
        lk_up = (logk[:, isig, js] + logk[:, isig, js + len(systs)]).reshape(shape)
        lk_dn = (logk[:, isig, js] - logk[:, isig, js + len(systs)]).reshape(shape)
    else:
        lk_up = logk[:, isig, js].reshape(shape); lk_dn = lk_up
    n0 = norm[:, isig].reshape(shape)[ipt, ieta]
    up = (n0 * np.exp(theta * lk_up[ipt, ieta])).sum(axis=0)
    dn = (n0 * np.exp(-theta * lk_dn[ipt, ieta])).sum(axis=0)
    return nom, up, dn


def fig_templates():
    edges = np.linspace(5.05, 5.45, 21)
    TH = 0.2
    fig, axs = plt.subplots(2, 2, figsize=(12, 5.2), sharex=True, sharey="row",
                            gridspec_kw={"height_ratios": [2.3, 1]})
    for col, (fname, title) in enumerate((
            ("btojpsik_cvhnano_tensor_joint24_mB_realmc_noZ_jpsifull.hdf5", "old unit"),
            ("btojpsik_cvhnano_tensor_joint24_matched_mB_noZ_jpsifull.hdf5", "matched unit"))):
        nom, up, dn = tensor_cell(fname, "e10", theta=TH)
        ax, rx = axs[0, col], axs[1, col]
        ax.stairs(nom, edges, color=INK, lw=2, label="nominal (θ = 0)")
        ax.stairs(up, edges, color=NEW, lw=1.8, label=f"θ(e10) = +{TH:g}")
        ax.stairs(dn, edges, color=OLD, lw=1.8, ls="--", label=f"θ(e10) = −{TH:g}")
        ax.set_title(title)
        rx.stairs(up / nom, edges, color=NEW, lw=1.8)
        rx.stairs(dn / nom, edges, color=OLD, lw=1.8, ls="--")
        rx.axhline(1, color=INK2, lw=0.8)
        rx.set_ylim(0.2, 2.8)
        rx.set_ylabel("ratio to nominal")
        rx.set_xlabel("m(μμK) [GeV]")
        for a in (ax, rx):
            style(a)
        ax.set_xlim(5.15, 5.40)
    axs[0, 0].set_ylabel("expected events")
    axs[0, 0].legend(frameon=False, fontsize=10)
    save(fig, "templates_e10_old_vs_matched.png")


# ------------------------------------------ 8. information adds: parabolas
def fig_parabolas():
    sj = summarise(f"{F}/realmc_D_jpsionly_asimov_fit.hdf5", "x", 1)["sige"]
    so = summarise(f"{F}/realmc_C_asimov_fit.hdf5", "x", 1)["sige"]
    sm = summarise(f"{F}/joint_mz_jpsi_bplus_matched_mB_noZ_jpsimu_jpsifull_fit.hdf5", "x", 1)["sige"]
    th = np.linspace(-0.2, 0.2, 401)
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.3), sharey=True)
    for ax, sc, title in ((axs[0], so, "old unit"), (axs[1], sm, "matched unit")):
        ib = max(1 / sc ** 2 - 1 / sj ** 2, 1e-12)
        sb = 1 / np.sqrt(ib)
        ax.plot(th, (th / sj) ** 2, color=REF, lw=2, label=f"J/ψ alone   σ = {sj:.3f}")
        ax.plot(th, (th / sb) ** 2, color=OLD if title == "old unit" else NEW, lw=2, ls="--",
                label=f"B⁺ term     σ = {sb:.3f}" if sb < 10 else f"B⁺ term     σ = {sb:.1f}")
        ax.plot(th, (th / sc) ** 2, color=INK, lw=2.4, label=f"combined   σ = {sc:.3f}")
        ax.axhline(1, color=INK2, lw=0.8, ls=":")
        ax.text(0.19, 1.15, "−2ΔlnL = 1  →  ±1σ", ha="right", fontsize=10, color=INK2)
        ax.set_ylim(0, 6)
        ax.set_xlim(-0.2, 0.2)
        ax.set_title(title)
        ax.set_xlabel("θ(e_i) [units of 5×10⁻² GeV], mean over η")
        ax.legend(frameon=False, fontsize=10, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=1)
        style(ax)
    axs[0].set_ylabel("−2 Δ ln L")
    save(fig, "information_parabolas_e.png")
    print(f"  parabolas: sigma J/psi {sj:.4f}, old {so:.4f}, matched {sm:.4f}")


# ------------------------------------- 9. degeneracy toy: (mB, theta) ellipses
def fig_degeneracy():
    slopes = np.array([30.3, 17.0, 7.2])      # MeV per unit of e, small-theta slopes (02_closure)
    sig_c = 0.069 * np.sqrt(3)               # per-cell peak precision -> 0.069 MeV combined
    sig_j = 0.068                            # J/psi-alone constraint on theta (prior units)

    def cov(k, with_jpsi=True):
        Fm = np.zeros((2, 2))
        for s in k * slopes:
            g = np.array([1.0, s])
            Fm += np.outer(g, g) / sig_c ** 2
        if with_jpsi:
            Fm[1, 1] += 1 / sig_j ** 2
        return np.linalg.inv(Fm)

    def ellipse(C, n=200):
        w, v = np.linalg.eigh(C)
        t = np.linspace(0, 2 * np.pi, n)
        return (v @ (np.sqrt(w)[:, None] * np.vstack([np.cos(t), np.sin(t)]))).T

    fig, axs = plt.subplots(1, 2, figsize=(12, 4.6), sharey=True)
    for ax, k, title in ((axs[0], 1 / 41, "old unit (slopes ÷ 41)"), (axs[1], 1.0, "matched unit")):
        Cb = cov(k, with_jpsi=False)
        Cc = cov(k)
        ax.axhspan(-sig_j, sig_j, color=REF, alpha=0.25, lw=0, label="J/ψ: ±1σ on θ")
        if np.all(np.isfinite(Cb)) and np.sqrt(Cb[1, 1]) < 10:
            eb = ellipse(Cb)
            ax.plot(eb[:, 0], eb[:, 1], color=OLD if k < 1 else NEW, lw=5 if k == 1 else 2,
                    ls="-" if k == 1 else "--", alpha=0.75, label="B⁺ alone, 1σ")
        ec = ellipse(Cc)
        ax.plot(ec[:, 0], ec[:, 1], color=INK, lw=2.4, label="combined, 1σ")
        smb, sth = np.sqrt(Cc[0, 0]), np.sqrt(Cc[1, 1])
        ax.annotate("", xy=(smb, -0.135), xytext=(-smb, -0.135),
                    arrowprops=dict(arrowstyle="<->", color=INK, lw=1.4))
        ax.text(0, -0.128, f"σ(mB) = {smb:.3f} MeV", ha="center", fontsize=10.5, color=INK)
        ax.set_title(f"{title}:  σ(θ) = {sth:.3f}")
        ax.set_xlim(-0.6, 0.6)
        ax.set_ylim(-0.15, 0.15)
        ax.set_xlabel("δmB [MeV]")
        style(ax)
        ax.legend(frameon=False, fontsize=10, loc="upper left")
        print(f"  toy k={k:.4f}: sigma(mB) {smb:.4f} MeV, sigma(theta) {sth:.4f}, "
              f"B+ alone sigma(theta) {np.sqrt(Cb[1, 1]):.4f}")
    axs[0].set_ylabel("θ (one global e) [prior units]")
    save(fig, "degeneracy_toy.png")


# ------------------------------------------------ 10. rendered display equations
# Each block: lines of mathtext, centred, drawn at the slide's display width so
# the deck places them at scale ~1 and every block shows at the same size.
EQUATIONS = {
    "eq_kaon_scale": [
        r"$m_B^2 = (E_\psi + E_K)^2 - |\vec{p}_\psi + \vec{p}_K|^2,\qquad \vec{p}_K \to (1+\delta)\,\vec{p}_K$",
        r"$\dfrac{\partial m_B}{\partial \delta} = \dfrac{(E_\psi + E_K)\,p_K^2/E_K \;-\; p_K^2 \;-\; \vec{p}_\psi\cdot\vec{p}_K}{m_B}$",
        r"B rest frame:  $\vec{p}_\psi\cdot\vec{p}_K = -p_K^2\,,\qquad E_\psi + E_K = m_B$",
        r"$\Rightarrow\;\; \dfrac{\partial m_B}{\partial \delta} = \dfrac{p_K^2}{E_K} = \dfrac{(1683\ \mathrm{MeV})^2}{1754\ \mathrm{MeV}} = 1615\ \mathrm{MeV}$",
        r"$\delta = 10^{-3}\;\;\Rightarrow\;\; \Delta m_B = 1.6\ \mathrm{MeV}$",
    ],
    "eq_nuisance": [
        r"$\nu_b(\theta) = \nu_b^0\,\exp(\theta\,\ln\kappa_b),\qquad \kappa_b = \mathrm{up}_b\,/\,\nu_b^0$",
        r"$L(\theta) = \prod_b\, \mathrm{Pois}\left(n_b \,|\, \nu_b(\theta)\right)\;\times\;\mathcal{N}(\theta\,|\,0,1)$",
        r"$A:\; k \to (1+A)\,k \qquad e:\; \delta k \propto e\,k^2 \qquad M:\; \delta k = q\,M \qquad (k = 1/p_T)$",
    ],
    "eq_stretch": [
        r"$\ln\kappa_5 = 5\,\ln\kappa_1 \qquad (A,\,e;\;\; M\ \mathrm{unchanged})$",
        r"$\nu(\theta) = \nu^0\, e^{\theta\,\cdot\, 5\ln\kappa_1} = \nu^0\, e^{(5\theta)\,\ln\kappa_1}$",
        r"the 1/1/1 nuisance, with $\theta \to 5\,\theta$:  one unit becomes $5\times10^{-3}$ / $5\times10^{-2}$ GeV",
    ],
    "eq_curvature": [
        r"$-2\,\Delta\ln L(\theta) \;\approx\; \dfrac{(\theta - \hat{\theta})^2}{\sigma^2},\qquad \dfrac{1}{\sigma^2} = \dfrac{\partial^2(-\ln L)}{\partial\theta^2}$",
        r"$\dfrac{1}{\sigma^2} = \dfrac{1}{\sigma^2_{J/\psi}} + \dfrac{1}{\sigma^2_{B^+}} + 1\;\;(\mathrm{prior})$",
        r"$\dfrac{1}{\sigma^2_{B^+}} \;\propto\; \sum_b \dfrac{1}{\nu_b}\left(\dfrac{\partial \nu_b}{\partial\theta}\right)^2 \qquad \mathrm{slope}\ \div 41 \;\Rightarrow\; \mathrm{information}\ \div 1700$",
    ],
    "eq_profiling": [
        r"$C = H^{-1},\qquad H_{ij} = \dfrac{\partial^2(-\ln L)}{\partial\theta_i\,\partial\theta_j},\qquad \rho_{ij} = \dfrac{C_{ij}}{\sigma_i\,\sigma_j}$",
        r"$\mathrm{impact}(G) = \sqrt{\sigma^2(\mathrm{POI}) - \sigma^2(\mathrm{POI}\,|\,G\ \mathrm{fixed})}$",
        r"$\sigma^2(\mathrm{POI}) \;\approx\; \sigma^2_{\mathrm{rest}} + g^{T} C_\theta\, g,\qquad g = \partial\,\mathrm{POI}/\partial\theta$",
    ],
    "eq_two_params": [
        r"$p_c = m_B + s_c\,\theta \qquad (\mathrm{cell}\ c = p_T \times \eta \times q)$",
        r"$\bar{p} = m_B + \bar{s}\,\theta\;\;(\mathrm{degenerate\ with}\ m_B) \qquad p_c - p_{c'} = (s_c - s_{c'})\,\theta\;\;(\mathrm{measures}\ \theta)$",
        r"$\sigma^2(m_B) \;\approx\; \sigma^2_{\mathrm{stat}} + \bar{s}^{T} C_\theta\,\bar{s}$",
    ],
    "eq_weighted_average": [
        r"$\hat{\theta} = \dfrac{\theta_J/\sigma_J^2 + \theta_B/\sigma_B^2}{1/\sigma_J^2 + 1/\sigma_B^2 + 1}$",
        r"$\Delta\chi^2 = \dfrac{(\theta_J - \theta_B)^2}{\sigma_J^2 + \sigma_B^2}$",
    ],
    "eq_gof": [
        r"$2\Delta\mathrm{NLL}_{\mathrm{sat}} = -2\,\ln\dfrac{L(\hat{\theta})}{L_{\mathrm{sat}}},\qquad \mathrm{good\ model:}\;\; \approx n_{\mathrm{dof}} \pm \sqrt{2\,n_{\mathrm{dof}}}$",
        r"$\mathrm{same\ data,\ same}\ n_{\mathrm{dof}}:\;\; \Delta(2\Delta\mathrm{NLL}) = 64\,266 - 61\,987 = +2\,280$",
    ],
}


def fig_equations():
    import matplotlib as mpl
    with mpl.rc_context({"mathtext.fontset": "cm", "font.family": "serif",
                         "font.serif": ["DejaVu Serif"]}):
        for name, lines in EQUATIONS.items():
            tall = [any(k in ln for k in ("dfrac", "sqrt", "sum", "prod")) for ln in lines]
            heights = [1.05 if t else 0.62 for t in tall]
            H = sum(heights) + 0.1
            fig = plt.figure(figsize=(11.0, H))
            y = H - 0.05
            for ln, h in zip(lines, heights):
                y -= h / 2
                fig.text(0.5, y / H, ln, ha="center", va="center", fontsize=21, color=INK)
                y -= h / 2
            fig.savefig(os.path.join(OUT, f"{name}.png"), dpi=200, transparent=True)
            plt.close(fig)
            print("wrote", f"{name}.png")


if __name__ == "__main__" and len(sys.argv) > 1:
    for name in sys.argv[1:]:
        globals()[name]()
elif __name__ == "__main__":
    fig_units()
    fig_slopes()
    fig_linearity()
    fig_mB()
    fig_mZ()
    fig_postfit()
    fig_templates()
    fig_parabolas()
    fig_degeneracy()
    fig_equations()
