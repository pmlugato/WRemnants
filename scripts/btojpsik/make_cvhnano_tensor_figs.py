"""Figures for the CVH NanoAOD B+ -> J/psi K tensor+fit deck.

Every number is measured on the full data (22 327 files) and the full simulation
(2 520 files); none is projected.

Simulation uncertainties come from the `fitbins` Weight() storage, so sigma is
sqrt(sum w^2). Where a simulated curve is drawn, its band is drawn with it -- on
the signal template, on the background template, and on their sum in quadrature.

For the *signal* template that reduces to sqrt(N), and the figures say so rather
than implying more: the forced-decay weight is normalised to the B+ precisely so
the signal category keeps weight 1, and these histograms were made without
pileup weights. Measured: sum = 77 413 exactly, effective entries = 77 413. The
variance is still read rather than assumed, because it does not reduce for the
non-signal categories -- 402 814 weighted events there carry the statistical
power of 278 650 -- and because a pileup weight would break the equality for the
signal too.
"""

import json

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import least_squares

from wremnants.postprocessing import rabbit_cvhnano_helpers as H
from wremnants.production.btojpsik_cvhnano_axes import (
    GEN_DATA,
    GEN_NON_SIGNAL_CATEGORIES,
    GEN_SIGNAL,
)

OUT = "/home/submit/pmlugato/public_html/slides/assets_cvhnano_tensor"
HIST = "/ceph/submit/data/user/p/pmlugato/mz/histograms/btojpsik_cvhnano_full.hdf5"

# The quantile edges the tensor writer derived on the full simulation; hard-coded
# here so a figure and a tensor cannot silently disagree about a cell.
PT_EDGES = [0.0, 2.3, 5.2, 10.0]
KW = dict(
    mass_low=5.20,
    mass_high=5.35,
    mass_bins=10,
    signal_sigma=0.0138,
    pt_edges=PT_EDGES,
    eta_rebin=1,
)

plt.rcParams.update(
    {
        "font.size": 12,
        "figure.dpi": 140,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "axes.axisbelow": True,
    }
)
C_DATA, C_SIG, C_BKG, C_TOT = "black", "#c0392b", "#2980b9", "#27ae60"

facts = {}


# --------------------------------------------------------------------------
# uncertainty helpers
# --------------------------------------------------------------------------
def vals_errs(h, keep_axis):
    """Sum a hist down to one axis, returning (values, sigma).

    sigma is sqrt of the summed Weight() variance, which is the only correct
    error for a weighted histogram.
    """
    v, var = h.values(), h.variances()
    if var is None:
        raise ValueError("histogram carries no variances; the fit input must be Weight()")
    axes = tuple(i for i in range(v.ndim) if i != keep_axis)
    return v.sum(axis=axes), np.sqrt(var.sum(axis=axes))


def frac_errs(n, var):
    """Uncertainty on n_i / sum(n), propagated through the shared denominator.

    Ignoring the denominator's correlation with the numerator overstates the
    error on a bin holding most of the sample -- which is exactly the soft-pT
    bin this is used on.
    """
    N = n.sum()
    f = n / N
    return np.sqrt(((1.0 - f) / N) ** 2 * var + (f / N) ** 2 * (var.sum() - var))


def eff_entries(n, var):
    """Effective unweighted entries, (sum w)^2 / sum w^2."""
    return float(n.sum() ** 2 / var.sum()) if var.sum() > 0 else 0.0


def hstep(ax, edges, y, color, label=None, lw=1.8, ls="-"):
    """Draw a binned quantity as a histogram.

    These are histograms, so they are drawn as histograms. A line through the bin
    centres asserts an interpolation between them that the content does not have,
    and on a ten-bin mass window it makes a step function look like a smooth peak.
    """
    return ax.stairs(y, edges, color=color, lw=lw, ls=ls, label=label, baseline=None)


def band(ax, edges, y, err, color, alpha=0.28, label=None):
    """Uncertainty band on a histogram, stepped to match it."""
    return ax.stairs(y + err, edges, baseline=y - err, fill=True, color=color,
                     alpha=alpha, lw=0, label=label)


# --------------------------------------------------------------------------
# load
# --------------------------------------------------------------------------
d_raw, m_raw, keys = H.load_fitbins(HIST, "Charmonium", "InclusiveBToJpsiX")
facts["datasets"] = list(keys)

# --------------------------------------------------------------------------
# 1. the pT overflow fold
# --------------------------------------------------------------------------
m_fold = H.strip_flow(H.fold_overflow(m_raw, H.PT_AXIS))
m_nofold = H.strip_flow(m_raw)
sig_fold = H.select_category(m_fold, GEN_SIGNAL)
sig_nofold = H.select_category(m_nofold, GEN_SIGNAL)
ipt = list(sig_fold.axes.name).index(H.PT_AXIS)
pt_ax = sig_fold.axes[H.PT_AXIS]

spec_fold, spec_fold_e = vals_errs(sig_fold, ipt)
spec_nofold, spec_nofold_e = vals_errs(sig_nofold, ipt)
facts["signal_fold"] = float(spec_fold.sum())
facts["signal_nofold"] = float(spec_nofold.sum())
facts["overflow_frac"] = float((spec_fold.sum() - spec_nofold.sum()) / spec_fold.sum())
facts["signal_fold_err"] = float(np.sqrt((spec_fold_e**2).sum()))
facts["signal_eff_entries"] = eff_entries(
    spec_fold, sig_fold.variances().sum(axis=tuple(i for i in range(sig_fold.values().ndim) if i != ipt))
)

fig, ax = plt.subplots(figsize=(7.4, 4.1))
ctr = pt_ax.centers
pt_edges_stored = pt_ax.edges
hstep(ax, pt_edges_stored, spec_nofold, C_BKG, lw=1.6, label="stored bins only (what every earlier count used)")
band(ax, pt_edges_stored, spec_nofold, spec_nofold_e, C_BKG)
hstep(ax, pt_edges_stored, spec_fold, C_SIG, lw=1.6, ls="--", label="with $p_T>10$ GeV folded into the top bin")
band(ax, pt_edges_stored, spec_fold, spec_fold_e, C_SIG)
ax.set_yscale("log")
ax.set_xlabel("bachelor $p_T$ [GeV]")
ax.set_ylabel("simulated signal events")
ax.annotate(
    f"folded overflow = {spec_fold[-1]-spec_nofold[-1]:,.0f} events\n"
    f"{100*facts['overflow_frac']:.1f}% of all simulated signal",
    xy=(9.9, spec_fold[-1]),
    xytext=(3.2, spec_fold[-1] * 0.014),
    arrowprops=dict(arrowstyle="->", color=C_SIG),
    color=C_SIG,
    fontsize=9,
)
ax.legend(fontsize=8.5, loc="upper right")
fig.tight_layout()
fig.savefig(f"{OUT}/fig_overflow.png", bbox_inches="tight")
plt.close(fig)

# --------------------------------------------------------------------------
# reduce to the fit binning
# --------------------------------------------------------------------------
d_red = H.reduce_fitbins(d_raw, **KW)
m_red = H.reduce_fitbins(m_raw, **KW)
D3 = H.select_category(d_red, GEN_DATA)
S3 = H.select_category(m_red, GEN_SIGNAL)
O3 = H.sum_categories(m_red, GEN_NON_SIGNAL_CATEGORIES)
im = list(S3.axes.name).index(H.MASS_AXIS)
mass_c = S3.axes[H.MASS_AXIS].centers
mass_e = S3.axes[H.MASS_AXIS].edges

# --------------------------------------------------------------------------
# 2. pT spectrum mismatch, with errors on both
# --------------------------------------------------------------------------
d_pt, _ = vals_errs(D3, 0)
s_pt, _ = vals_errs(S3, 0)
d_pt_var = D3.variances().sum(axis=(1, 2, 3))
s_pt_var = S3.variances().sum(axis=(1, 2, 3))
fd, fs = d_pt / d_pt.sum(), s_pt / s_pt.sum()
fd_e, fs_e = frac_errs(d_pt, d_pt_var), frac_errs(s_pt, s_pt_var)
facts["pt_frac_data"] = fd.tolist()
facts["pt_frac_signal"] = fs.tolist()
facts["pt_frac_signal_err"] = fs_e.tolist()
facts["pt_eff_entries"] = [eff_entries(np.array([s_pt[i]]), np.array([s_pt_var[i]])) for i in range(3)]

fig, ax = plt.subplots(figsize=(7.4, 4.1))
x = np.arange(3)
w = 0.38
lbl = [f"{PT_EDGES[i]:.1f}–{PT_EDGES[i+1]:.1f}" for i in range(3)]
ax.bar(x - w / 2, fd, w, yerr=fd_e, color=C_DATA, label="data", capsize=3, error_kw=dict(ecolor="#777777", lw=1.2))
ax.bar(x + w / 2, fs, w, yerr=fs_e, color=C_SIG, label="simulated signal", capsize=3, error_kw=dict(ecolor="#555555", lw=1.2))
for i in range(3):
    ax.text(i - w / 2, fd[i] + 0.02, f"{fd[i]:.3f}", ha="center", fontsize=9)
    ax.text(i + w / 2, fs[i] + 0.02, f"{fs[i]:.3f}", ha="center", fontsize=9, color=C_SIG)
ax.set_xticks(x)
ax.set_xticklabels(lbl)
ax.set_xlabel("bachelor $p_T$ [GeV]")
ax.set_ylabel("fraction of events in window")
ax.set_ylim(0, 1.0)
ax.legend(fontsize=9)
fig.tight_layout()
fig.savefig(f"{OUT}/fig_ptmismatch.png", bbox_inches="tight")
plt.close(fig)

# --------------------------------------------------------------------------
# inclusive shapes
# --------------------------------------------------------------------------
ax_other = tuple(i for i in range(D3.values().ndim) if i != im)
D = D3.values().sum(axis=ax_other)
Dv = D3.variances().sum(axis=ax_other)
S = S3.values().sum(axis=ax_other)
Sv = S3.variances().sum(axis=ax_other)
O = O3.values().sum(axis=ax_other)
Ov = O3.variances().sum(axis=ax_other)
Se, Oe = np.sqrt(Sv), np.sqrt(Ov)

facts.update(
    data_window=float(D.sum()),
    mc_signal_window=float(S.sum()),
    mc_other_window=float(O.sum()),
    mc_signal_window_err=float(np.sqrt(Sv.sum())),
    mc_signal_eff_entries=eff_entries(S, Sv),
    mc_other_eff_entries=eff_entries(O, Ov),
)

# --------------------------------------------------------------------------
# 3. the background shape, with its own errors
# --------------------------------------------------------------------------
x01 = (mass_c - mass_c[0]) / (mass_c[-1] - mass_c[0])
On = O / O.sum()
On_e = frac_errs(O, Ov)


def fit_curve(model, p0):
    r = least_squares(lambda p: (model(p, x01) - On) / On_e, p0)
    return model(r.x, x01), r


f_expo = lambda p, x: np.exp(p[0] + p[1] * x)
f_line = lambda p, x: p[0] * (1 - x) + p[1] * x
expo, r_expo = fit_curve(f_expo, [np.log(On.mean()), -0.5])
line, r_line = fit_curve(f_line, [On[0], On[-1]])
# Drawn on a fine grid: these two are genuinely smooth functions of the mass, and
# stepping them would hide that it is the *simulated histogram* that has structure.
mass_fine = np.linspace(mass_e[0], mass_e[-1], 400)
x01_fine = (mass_fine - mass_c[0]) / (mass_c[-1] - mass_c[0])
facts["bkg_shape"] = On.tolist()
facts["bkg_shape_err"] = On_e.tolist()
facts["bkg_shape_max_over_min"] = float(On.max() / On.min())
facts["bkg_argmax_mass"] = float(mass_c[np.argmax(On)])

fig, ax = plt.subplots(figsize=(7.6, 4.4))
hstep(ax, mass_e, On, C_BKG, lw=2.0, label="simulated non-signal (the actual background shape)")
band(ax, mass_e, On, On_e, C_BKG, alpha=0.45)
ax.plot(mass_fine, f_expo(r_expo.x, x01_fine), color=C_SIG, lw=1.8, label="best-fit exponential  $e^{a+bx}$")
ax.plot(mass_fine, f_line(r_line.x, x01_fine), color="#8e44ad", lw=1.8, ls="--", label="best-fit straight line (order-1 Bernstein)")
ax.axvline(H.B_MASS, color="grey", lw=1, ls=":")
lo, hi = ax.get_ylim()
ax.set_ylim(lo, hi + 0.30 * (hi - lo))  # headroom for the legend, which otherwise sits on a data point
lo, hi = ax.get_ylim()
ax.text(H.B_MASS + 0.0012, lo + 0.55 * (hi - lo), "signal peak", fontsize=8.5, color="grey", rotation=90, va="center")
ax.set_xlabel("$m(J/\\psi K)$ [GeV]")
ax.set_ylabel("fraction of background per bin")
ax.legend(fontsize=8.5, loc="upper left")
fig.tight_layout()
fig.savefig(f"{OUT}/fig_bkgshape.png", bbox_inches="tight")
plt.close(fig)

# --------------------------------------------------------------------------
# 4. best achievable fit of each background family, two panels
# --------------------------------------------------------------------------
r1 = least_squares(lambda p: (p[0] * S + np.exp(p[1] + p[2] * x01) * np.mean(D) - D) / np.sqrt(D), [15.0, 0.0, -0.5])
m1 = r1.x[0] * S + np.exp(r1.x[1] + r1.x[2] * x01) * np.mean(D)
s1, s1e = r1.x[0] * S, r1.x[0] * Se
m1e = s1e  # analytic background carries no simulated statistical error

r2 = least_squares(lambda p: (p[0] * S + (p[1] * (1 - x01) + p[2] * x01) * np.mean(D) - D) / np.sqrt(D), [15.0, 1.0, 1.0])
m2 = r2.x[0] * S + (r2.x[1] * (1 - x01) + r2.x[2] * x01) * np.mean(D)

r3 = least_squares(lambda p: (p[0] * S + p[1] * O - D) / np.sqrt(D), [15.0, 20.0])
m3 = r3.x[0] * S + r3.x[1] * O
s3v, s3e = r3.x[0] * S, r3.x[0] * Se
b3, b3e = r3.x[1] * O, r3.x[1] * Oe
m3e = np.sqrt(s3e**2 + b3e**2)


def chi2ndf(model, npar):
    return float(np.sum(((model - D) / np.sqrt(D)) ** 2)) / (len(D) - npar)


facts["inclusive_fits"] = {
    "signal + exponential": dict(chi2_ndf=chi2ndf(m1, 3), signal_yield=float(s1.sum()), frac=float(s1.sum() / D.sum())),
    "signal + straight line": dict(chi2_ndf=chi2ndf(m2, 3), signal_yield=float(r2.x[0] * S.sum()), frac=float(r2.x[0] * S.sum() / D.sum())),
    "signal + simulated background": dict(chi2_ndf=chi2ndf(m3, 2), signal_yield=float(s3v.sum()), frac=float(s3v.sum() / D.sum())),
    "mc_bkg_fitted_scale": float(r3.x[1]),
}

panels = [
    ("signal + exponential  (analytic)", m1, m1e, s1, s1e, None, None, chi2ndf(m1, 3)),
    ("signal + simulated background", m3, m3e, s3v, s3e, b3, b3e, chi2ndf(m3, 2)),
]
fig, axs = plt.subplots(2, 2, figsize=(11.0, 5.4), sharex=True, gridspec_kw={"height_ratios": [2.4, 1]})
for j, (t, tot, tote, sg, sge, bk, bke, c2) in enumerate(panels):
    A, B = axs[0, j], axs[1, j]
    hstep(A, mass_e, tot / 1e6, C_TOT, lw=2.2, label="total")
    band(A, mass_e, tot / 1e6, tote / 1e6, C_TOT, alpha=0.35)
    hstep(A, mass_e, sg / 1e6, C_SIG, lw=1.7, ls="--", label="signal")
    band(A, mass_e, sg / 1e6, sge / 1e6, C_SIG)
    bkv = (tot - sg) if bk is None else bk
    hstep(A, mass_e, bkv / 1e6, C_BKG, lw=1.7, ls=":", label="background")
    if bke is not None:
        band(A, mass_e, bkv / 1e6, bke / 1e6, C_BKG)
    A.errorbar(mass_c, D / 1e6, yerr=np.sqrt(Dv) / 1e6, fmt="o", color=C_DATA, ms=4, label="data", zorder=5)
    A.set_title(f"{t}\n$\\chi^2$/ndf = {c2:,.0f}   signal = {sg.sum()/1e6:.2f}M", fontsize=11)
    A.set_ylim(0, 1.85)
    if j == 0:
        A.set_ylabel("events / 15 MeV  [$10^6$]")
        A.legend(fontsize=9, loc="upper left")
    B.axhline(1, color=C_TOT, lw=1.5)
    band(B, mass_e, np.ones_like(tot), tote / tot, C_TOT, alpha=0.35)
    B.errorbar(mass_c, D / tot, yerr=np.sqrt(Dv) / tot, fmt="o", color=C_DATA, ms=4)
    B.set_ylim(0.86, 1.14)
    B.set_xlabel("$m(J/\\psi K)$ [GeV]")
    if j == 0:
        B.set_ylabel("data / model")
fig.tight_layout()
fig.savefig(f"{OUT}/fig_models2.png", bbox_inches="tight")
plt.close(fig)

# --------------------------------------------------------------------------
# 5. prefit at two edge cells, with the template's own error band
# --------------------------------------------------------------------------
sc = H.signal_seed_scale(D3, S3, 0.03)
sc, cap = H.cap_seed_scale(sc, D3, S3)
facts["seed_used"] = float(sc)
facts["seed_cap"] = float(cap)
res = np.clip(D3.values().sum(axis=im, keepdims=True) - sc * S3.values().sum(axis=im, keepdims=True), 0, None)
bkg = np.broadcast_to(
    np.maximum(res, 0.05 * D3.values().sum(axis=im, keepdims=True)) / len(mass_c), D3.values().shape
)

cells = [
    (0, 1, 0, "softest $p_T$ (0.0–2.3 GeV),  $\\eta \\approx -2.35$,  $q=-1$"),
    (2, 1, 0, "hardest $p_T$ (5.2–10 GeV),  $\\eta \\approx -2.35$,  $q=-1$"),
]
facts["edge_cells"] = []
fig, axs = plt.subplots(2, 2, figsize=(11.0, 5.4), sharex=True, gridspec_kw={"height_ratios": [2.4, 1]})
for j, (a_, b_, c_, t) in enumerate(cells):
    dd = D3.values()[a_, b_, c_, :]
    dde = np.sqrt(D3.variances()[a_, b_, c_, :])
    sraw = S3.values()[a_, b_, c_, :]
    sraw_e = np.sqrt(S3.variances()[a_, b_, c_, :])
    ss, sse = sc * sraw, sc * sraw_e
    bb = bkg[a_, b_, c_, :]
    tt, tte = ss + bb, sse  # only the template contributes a simulated error
    facts["edge_cells"].append(
        dict(
            label=t,
            data=float(dd.sum()),
            mc_signal=float(sraw.sum()),
            mc_signal_err=float(np.sqrt((sraw_e**2).sum())),
            mc_eff_entries=eff_entries(sraw, S3.variances()[a_, b_, c_, :]),
            max_rel_err=float(np.nanmax(np.where(sraw > 0, sraw_e / sraw, np.nan))),
        )
    )
    A, B = axs[0, j], axs[1, j]
    hstep(A, mass_e, tt, C_TOT, lw=2.2, label="prefit total")
    band(A, mass_e, tt, tte, C_TOT, alpha=0.35)
    hstep(A, mass_e, ss, C_SIG, lw=1.7, ls="--", label="signal, seeded")
    band(A, mass_e, ss, sse, C_SIG, label="signal, simulated stat.")
    hstep(A, mass_e, bb, C_BKG, lw=1.7, ls=":", label="background, seeded")
    A.errorbar(mass_c, dd, yerr=dde, fmt="o", color=C_DATA, ms=4, label="data", zorder=5)
    A.set_title(
        f"{t}\n{dd.sum():,.0f} data,  {sraw.sum():.0f} simulated signal "
        f"($N_{{\\rm eff}}$ = {facts['edge_cells'][-1]['mc_eff_entries']:.0f})",
        fontsize=10.5,
    )
    if j == 0:
        A.set_ylabel("events / 15 MeV")
        A.legend(fontsize=8.5)
    B.axhline(1, color=C_TOT, lw=1.5)
    band(B, mass_e, np.ones_like(tt), tte / tt, C_TOT, alpha=0.35)
    B.errorbar(mass_c, dd / tt, yerr=dde / tt, fmt="o", color=C_DATA, ms=4)
    B.set_xlabel("$m(J/\\psi K)$ [GeV]")
    if j == 0:
        B.set_ylabel("data / prefit")
fig.tight_layout()
fig.savefig(f"{OUT}/fig_edge2.png", bbox_inches="tight")
plt.close(fig)

# --------------------------------------------------------------------------
# 6. occupancy, in effective entries as well as raw
# --------------------------------------------------------------------------
sig_cell = S3.values().sum(axis=im).ravel()
sig_cell_var = S3.variances().sum(axis=im).ravel()
with np.errstate(divide="ignore", invalid="ignore"):
    neff_cell = np.where(sig_cell_var > 0, sig_cell**2 / sig_cell_var, 0.0)
pk = H.peak_mask(S3, 0.03)
pk_cell = np.compress(pk, S3.values(), axis=im).sum(axis=im).ravel()
tot_cell = (S3.values() + O3.values()).sum(axis=im).ravel()
facts["occ"] = {
    k: dict(
        sum=float(v.sum()),
        min=float(v.min()),
        median=float(np.median(v)),
        max=float(v.max()),
        above=int((v >= 100).sum()),
        n=int(v.size),
    )
    for k, v in [
        ("total", tot_cell),
        ("signal", sig_cell),
        ("signal_peak", pk_cell),
        ("signal_neff", neff_cell),
    ]
}

fig, ax = plt.subplots(figsize=(7.4, 4.2))
bins = np.linspace(0, max(sig_cell.max(), pk_cell.max()) * 1.02, 34)
# The N_eff overlay is deliberately not drawn: it lies exactly on the weighted
# histogram, because the signal category is unweighted by construction. Stating
# that is worth more than a duplicate curve.
ax.hist(
    sig_cell, bins=bins, color=C_SIG, alpha=0.5,
    label=f"signal category (median {np.median(sig_cell):.0f}, and $N_{{\\rm eff}}$ is identical)",
)
ax.hist(pk_cell, bins=bins, color=C_BKG, alpha=0.5, label=f"signal in mass peak (median {np.median(pk_cell):.0f})")
ax.axvline(100, color="black", lw=1.6, ls="--")
ax.text(104, ax.get_ylim()[1] * 0.72, "100-event floor", fontsize=9)
ax.set_xlabel("simulated events in a $(p_T,\\eta,q)$ cell")
ax.set_ylabel("cells")
ax.legend(fontsize=8.5)
fig.tight_layout()
fig.savefig(f"{OUT}/fig_occupancy.png", bbox_inches="tight")
plt.close(fig)


# --------------------------------------------------------------------------
# 7. the commands, as monospace figures
#
# The slides skill has no code-block builder and a 200-character command line
# overflows an 18pt text box, so these are rendered as images. They are the only
# slides in the deck that are a picture of text, and that is deliberate: line
# breaks in a shell command have to survive exactly.
# --------------------------------------------------------------------------
def command_figure(path, lines, width=11.0, size=9.0):
    """Render a monospace block. A line starting with '#' is drawn grey."""
    fig_h = 0.20 * len(lines) + 0.24
    fig, ax = plt.subplots(figsize=(width, fig_h))
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, len(lines))
    for i, ln in enumerate(lines):
        stripped = ln.lstrip()
        ax.text(
            0.0,
            len(lines) - i - 0.5,
            # A shell variable is not a maths delimiter. Left unescaped,
            # matplotlib pairs the dollars and silently renders `$HIST -o $MZ`
            # as italic mathtext with the flag eaten.
            ln.replace("$", r"\$"),
            family="monospace",
            fontsize=size,
            color="#6b6b6b" if stripped.startswith("#") else "#1f1f1f",
            va="center",
            ha="left",
        )
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


command_figure(
    f"{OUT}/fig_cmd_hist.png",
    [
        "# Everything runs in the wmassdevrolling container, after `source setup.sh`.",
        "# MZ=/ceph/submit/data/user/p/pmlugato/mz",
        "",
        "# 1. histmaker, both datasets in ONE pass -> ONE file",
        "#    data 22 327 files / 18.46M candidates, simulation 2 520 / 887 100",
        "python scripts/histmakers/btojpsik_cvhnano.py --era 2016PostVFP \\",
        "    --dataPath /ceph/submit/data/group/cms/store \\",
        "    --filterProcs Charmonium_2016PostVFP_cvhNano InclusiveBToJpsiX_2016PostVFP \\",
        "    --allaxes --fitBinsHist --cutflow \\",
        "    --saveCutflow $MZ/cutflows/cvhnano_onepass \\",
        "    --noScaleToData --theoryCorr none \\",
        "    -o $MZ/histograms/ -p onepass",
        "",
        "# --filterProcs is variadic, so there is no merge step and no merge script.",
        "# The first campaign ran the two passes separately, purely so the simulation",
        "# could be iterated on without repeating the long data pass, and that is the",
        "# only reason a merge ever existed.",
        "",
        "# --arm defaults to jointCvh, so the fit mass is jointCvh_mass throughout.",
        "# --fitBinsHist books the 5D `fitbins` histogram the tensor is written from.",
        "",
        "# 2. the figures in this deck",
        "python scripts/btojpsik/make_cvhnano_tensor_figs.py",
    ],
)

command_figure(
    f"{OUT}/fig_cmd_scan.png",
    [
        "# The selection scan: one 300-file pass per working point, on MATCHED files,",
        "# so a shape difference cannot come from which files happened to be read.",
        "SEL='--era 2016PostVFP --dataPath /ceph/submit/data/group/cms/store'",
        "SEL=\"$SEL --filterProcs Charmonium_2016PostVFP_cvhNano InclusiveBToJpsiX_2016PostVFP\"",
        "SEL=\"$SEL --allaxes --fitBinsHist --noScaleToData --theoryCorr none --maxFiles 300\"",
        "H=scripts/histmakers/btojpsik_cvhnano.py",
        "",
        "python $H $SEL -p ref300                                    # the reference",
        "python $H $SEL -p combolo --minBachelorAbsD0 0.01 \\",
        "                          --bachelorMinValidHits 10         # RECOMMENDED",
        "python $H $SEL -p d0lo    --minBachelorAbsD0 0.01",
        "python $H $SEL -p d0hi    --minBachelorAbsD0 0.02",
        "python $H $SEL -p hits13  --bachelorMinValidHits 13",
        "python $H $SEL -p niter5  --maxJointCvhNiter 5              # fit-quality family",
        "python $H $SEL -p vtx03   --minVtxProb 0.3",
        "python $H $SEL -p chisqn  --maxChisqNdof 1.2",
        "",
        "# the 2018-equivalent configuration (the two non-transferable cuts excluded)",
        "python $H $SEL -p sel2018_300 \\",
        "    --muonMaxEta 1.4 --muonMinPt 4 --bachelorMaxEta 1.4 --bachelorMaxPt 8 \\",
        "    --minDimuonPt 7 --minDimuonSl3d 4 --minDimuonVtxProb 0.1",
        "",
        "# the comparison table, and the scan figure",
        "python scripts/btojpsik/compare_cvhnano_selections.py \\",
        "    --label ref     $MZ/histograms/selscan/..._ref300.hdf5 \\",
        "    --label combolo $MZ/histograms/selscan/..._combolo.hdf5",
        "python scripts/btojpsik/make_cvhnano_selection_figs.py",
    ],
)

command_figure(
    f"{OUT}/fig_cmd_fit.png",
    [
        "HIST=$MZ/histograms/btojpsik_cvhnano_full.hdf5",
        "",
        "# 5. one line per configuration; the driver does tensor -> fit -> plot-ready",
        "scripts/rabbit/btojpsik_cvhnano_fit.sh $HIST full                 # 288 cells, closes",
        "scripts/rabbit/btojpsik_cvhnano_fit.sh $HIST full --observed      # 288 cells, does not",
        "scripts/rabbit/btojpsik_cvhnano_fit.sh $HIST inclusive --observed # 2 cells,   does not",
        "BKGMODEL=bernstein ASIMOV=0 ...  $HIST full                      # order-1 straight line",
        "ASIMOV=0 ...  $HIST full --pseudoData mcTotal                    # simulation closure",
        "",
        "# what one of those expands to, written verbatim to",
        "#   $MZ/logs/cvhnano_full_exp_data_commands.sh",
        "python scripts/rabbit/btojpsik_cvhnano_tensor.py -i $HIST -o $MZ/tensors \\",
        "    -p full_exp_data --ptBins 3 --etaRebin 1",
        "",
        "python rabbit/bin/rabbit_fit.py $MZ/tensors/btojpsik_cvhnano_tensor_full_exp_data.hdf5 \\",
        "    -v 4 --diagnostics -o $MZ/fitresults --outname cvhnano_full_exp_data_fit.hdf5 \\",
        "    --saveHists --saveHistsPerProcess --computeHistErrors \\",
        "    --computeHistErrorsPerProcess --computeHistCov \\",
        "    --paramModel AxisNormModel btojpsik_cvhnano BuJpsiK \\",
        "                 fit_bachelor_pt,fit_bachelor_eta,fit_bachelor_charge \\",
        "    --paramModel AxisExpModel  btojpsik_cvhnano BuJpsiKBkg jointCvh_mass \\",
        "                 fit_bachelor_pt,fit_bachelor_eta,fit_bachelor_charge fit_bachelor_eta \\",
        "    -m Project btojpsik_cvhnano jointCvh_mass \\",
        "    -m Project btojpsik_cvhnano fit_bachelor_pt fit_bachelor_eta \\",
        "                                fit_bachelor_charge jointCvh_mass \\",
        "    -t 0                    #  -t 0 fits data;  -t -1 is Asimov, and the default",
    ],
)

json.dump(facts, open(f"{OUT}/facts.json", "w"), indent=1)
print(
    json.dumps(
        {k: v for k, v in facts.items() if k not in ("bkg_shape", "bkg_shape_err", "pt_frac_data", "pt_frac_signal")},
        indent=1,
    )
)
