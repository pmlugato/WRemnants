"""Figures for the CVH NanoAOD selection study.

Two things the earlier deck got wrong and these figures correct:

  1. The structure under the peak is not a partially reconstructed satellite. Over
     the whole stored mass range it is one broad maximum sitting on the B mass, and
     the fit window only ever showed its top.
  2. Purity and background flatness pull against each other, but only for cuts
     conditional on the B fit. Cuts on the bachelor track alone improve both.

Inputs are the selection-scan histmaker outputs. Run
`scripts/btojpsik/compare_cvhnano_selections.py` for the numbers these plot.
"""

import argparse
import json
import os

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from wremnants.postprocessing import rabbit_cvhnano_helpers as H
from wremnants.production.btojpsik_cvhnano_axes import (
    GEN_DATA,
    GEN_NON_SIGNAL_CATEGORIES,
    GEN_SIGNAL,
)

OUT = "/home/submit/pmlugato/public_html/slides/assets_cvhnano_tensor"
SCAN = "/ceph/submit/data/user/p/pmlugato/mz/histograms/selscan"
FULL = "/ceph/submit/data/user/p/pmlugato/mz/histograms/btojpsik_cvhnano_full.hdf5"
WINDOW = (5.20, 5.35)

plt.rcParams.update({"font.size": 12, "figure.dpi": 140, "axes.grid": True,
                     "grid.alpha": 0.25, "axes.axisbelow": True})
C_DATA, C_SIG, C_BKG, C_TOT = "black", "#c0392b", "#2980b9", "#27ae60"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--scan", default=SCAN, help="directory of selection-scan outputs")
args = parser.parse_args()
facts = {}


def get(o, n):
    h = o[n]
    return h.get() if hasattr(h, "get") else h


# --------------------------------------------------------------------------
# 1. the background over the whole stored range
# --------------------------------------------------------------------------
with H.open_results(FULL) as r:
    mo = r["InclusiveBToJpsiX_2016PostVFP"]["output"]
    hm = get(mo, "nominal")
    ax = hm.axes[0]
    ce, ed = ax.centers, ax.edges
    S = hm[:, GEN_SIGNAL].values()
    O1 = sum(hm[:, c].values() for c in GEN_NON_SIGNAL_CATEGORIES[:-1])
    O2 = hm[:, GEN_NON_SIGNAL_CATEGORIES[-1]].values()
    Ov = sum(hm[:, c].variances() for c in GEN_NON_SIGNAL_CATEGORIES)

O = O1 + O2
facts["fullrange_argmax"] = float(ce[np.argmax(O)])
facts["fullrange_maxmin"] = float(O.max() / O.min())
facts["combinatorial_frac"] = float(O2.sum() / O.sum())

fig, ax_ = plt.subplots(figsize=(7.8, 4.3))
ax_.stairs(O, ed, color=C_BKG, lw=2.0, label="simulated non-signal", baseline=None)
ax_.stairs(O + np.sqrt(Ov), ed, baseline=O - np.sqrt(Ov), fill=True, color=C_BKG,
           alpha=0.35, lw=0)
ax_.stairs(S, ed, color=C_SIG, lw=1.8, ls="--", label="simulated signal", baseline=None)
ax_.axvspan(*WINDOW, color="grey", alpha=0.16, lw=0)
ax_.text(np.mean(WINDOW), O.max() * 0.72, "the fit window\nonly saw this part",
         ha="center", va="top", fontsize=9, color="#444444")
ax_.axvline(H.B_MASS, color="grey", lw=1, ls=":")
ax_.annotate(f"background maximum at {facts['fullrange_argmax']:.4f}\n"
             f"— on the B mass, and {facts['fullrange_maxmin']:.0f}x the low edge",
             xy=(facts["fullrange_argmax"], O.max()),
             xytext=(5.345, O.max() * 0.80), fontsize=9, color=C_BKG,
             arrowprops=dict(arrowstyle="->", color=C_BKG))
ax_.set_xlabel("$m(J/\\psi K)$ [GeV]")
ax_.set_ylabel("simulated events / 5 MeV")
ax_.set_xlim(ed[0], ed[-1])
ax_.legend(fontsize=9, loc="upper left")
fig.tight_layout()
fig.savefig(f"{OUT}/fig_fullrange.png", bbox_inches="tight")
plt.close(fig)

# --------------------------------------------------------------------------
# 2. the selection scan: yield bias against signal efficiency
# --------------------------------------------------------------------------
# Recomputed here rather than parsed from the comparison script's stdout, so the
# figure cannot drift from the table.

def lstsq_amplitudes(components, target, weights):
    """Weighted linear least squares for the amplitudes of fixed shapes.

    Both fits here are linear in their amplitudes, so they are solved rather than
    iterated. `scipy.optimize.least_squares` in this build intermittently raises
    "Memory error in scipy.linalg.svd" on a 10x2 system, which is not a memory
    problem and not something to retry around -- there is no reason to run an
    optimiser on a problem with a closed form.
    """
    A = np.vstack([c / weights for c in components]).T
    b = target / weights
    return np.linalg.lstsq(A, b, rcond=None)[0]


def fit_signal_plus_exp(s, dd, w, x):
    """Signal amplitude with an exponential background.

    Conditionally linear: for a fixed slope the amplitudes solve exactly, so the
    slope is scanned on a grid and then refined. Deterministic, and it cannot fail
    to converge.
    """
    best = None
    for coarse in (np.linspace(-8.0, 8.0, 161),):
        for b in coarse:
            p = lstsq_amplitudes([s, np.exp(b * x)], dd, w)
            if p[0] < 0 or p[1] < 0:
                continue
            chi2 = float(np.sum(((p[0] * s + p[1] * np.exp(b * x) - dd) / w) ** 2))
            if best is None or chi2 < best[0]:
                best = (chi2, b, p)
    if best is None:  # every slope needed a negative amplitude
        p = lstsq_amplitudes([s, np.ones_like(x)], dd, w)
        return float(p[0])
    _, b0, _ = best
    for b in np.linspace(b0 - 0.1, b0 + 0.1, 41):
        p = lstsq_amplitudes([s, np.exp(b * x)], dd, w)
        if p[0] < 0 or p[1] < 0:
            continue
        chi2 = float(np.sum(((p[0] * s + p[1] * np.exp(b * x) - dd) / w) ** 2))
        if chi2 < best[0]:
            best = (chi2, b, p)
    return float(best[2][0])


# (label, filename tag, family)
POINTS = [
    ("no cuts", "ref300", "reference"),
    ("2018 selection", "sel2018_300", "kinematic"),
    ("nIter <= 5", "niter5", "fit quality"),
    ("B vtx prob > 0.3", "vtx03", "fit quality"),
    ("chisq/ndof < 1.2", "chisqn", "fit quality"),
    ("|d0| > 0.01", "d0lo", "track quality"),
    ("|d0| > 0.02", "d0hi", "track quality"),
    ("|d0| > 0.05", "d0hi2", "track quality"),
    ("hits > 13", "hits13", "track quality"),
    ("|d0|>0.01, hits>10", "combolo", "track quality"),
    ("|d0|>0.02, hits>13", "combo", "track quality"),
]
COLOUR = {"reference": "black", "kinematic": "#7f8c8d",
          "fit quality": "#c0392b", "track quality": "#2980b9"}

kw = dict(mass_low=WINDOW[0], mass_high=WINDOW[1], mass_bins=10, signal_sigma=0.0138,
          pt_edges=[0.0, 10.0], eta_rebin=48,
          collapse_axes=(H.PT_AXIS, H.ETA_AXIS, H.CHARGE_AXIS))

scan = []
for label, tag, family in POINTS:
    path = os.path.join(args.scan, f"btojpsik_cvhnano_maxFiles_300_{tag}.hdf5")
    if not os.path.exists(path):
        print(f"  missing, skipping: {tag}")
        continue
    d_raw, m_raw, _ = H.load_fitbins(path, "Charmonium", "InclusiveBToJpsiX")
    m, d = H.reduce_fitbins(m_raw, **kw), H.reduce_fitbins(d_raw, **kw)
    Sh = H.select_category(m, GEN_SIGNAL)
    Oh = H.sum_categories(m, GEN_NON_SIGNAL_CATEGORIES)
    Dh = H.select_category(d, GEN_DATA)
    im = list(Sh.axes.name).index(H.MASS_AXIS)
    oth = tuple(i for i in range(Sh.values().ndim) if i != im)
    s, o, dd = (h.values().sum(axis=oth) for h in (Sh, Oh, Dh))
    if s.sum() <= 0 or o.sum() <= 0 or dd.sum() <= 0:
        # A seed of inf reaches scipy's SVD as a memory error rather than as a
        # convergence failure, so an empty selection has to be caught here.
        print(f"  empty after selection, skipping: {tag} "
              f"(signal {s.sum():.0f}, bkg {o.sum():.0f}, data {dd.sum():.0f})")
        continue
    x = np.linspace(0.0, 1.0, len(dd))
    w = np.sqrt(np.maximum(dd, 1.0))
    a_true = float(lstsq_amplitudes([s, o], dd, w)[0])
    a_exp = fit_signal_plus_exp(s, dd, w, x)
    n = o / o.sum()
    scan.append(dict(label=label, family=family, signal=float(s.sum()),
                     bias=a_exp / a_true if a_true > 0 else np.nan,
                     sb=float(s.sum() / o.sum()),
                     maxmin=float(n.max() / n.min())))

ref = next(p for p in scan if p["family"] == "reference")
for p in scan:
    p["eff"] = p["signal"] / ref["signal"]
facts["scan"] = scan

# One panel, not two. The max/min panel told the same story with worse label
# collisions, and its numbers are in the table; this way the decisive plot is
# displayed at nearly twice the size on a slide.
OFFSET = {
    "no cuts": (-52, 10),
    "2018 selection": (10, -4),
    "nIter <= 5": (-58, 10),
    "B vtx prob > 0.3": (8, 6),
    "chisq/ndof < 1.2": (8, -14),
    "|d0| > 0.01": (8, 6),
    "hits > 13": (8, 6),
    "|d0| > 0.02": (-70, -4),
    "|d0|>0.01, hits>10": (10, -4),
    "|d0| > 0.05": (10, 4),
    "|d0|>0.02, hits>13": (10, 6),
}
fig, a = plt.subplots(figsize=(8.6, 4.7))
seen = set()
for p in scan:
    lbl = p["family"] if p["family"] not in seen else None
    seen.add(p["family"])
    a.plot(p["eff"], p["bias"], "o", ms=10, color=COLOUR[p["family"]], label=lbl, zorder=4)
    a.annotate(p["label"], (p["eff"], p["bias"]), textcoords="offset points",
               xytext=OFFSET.get(p["label"], (8, 6)), fontsize=8.5, color="#333333")
a.axhline(1.0, color=C_TOT, lw=1.6, ls="--")
a.text(0.315, 1.04, "no bias", fontsize=9, color=C_TOT)
a.set_xlabel("signal efficiency")
a.set_ylabel("yield bias: exponential / true shape")
a.set_ylim(0.92, 2.95)
a.set_xlim(0.30, 1.10)
a.legend(fontsize=9.5, loc="center left")
fig.tight_layout()
fig.savefig(f"{OUT}/fig_selscan.png", bbox_inches="tight")
plt.close(fig)

json.dump(facts, open(f"{OUT}/facts_selection.json", "w"), indent=1)
print(json.dumps({k: v for k, v in facts.items() if k != "scan"}, indent=1))
print(f"{'selection':<22}{'family':<15}{'sig eff':>9}{'S/B':>8}{'max/min':>9}{'bias':>8}")
for p in scan:
    print(f"{p['label']:<22}{p['family']:<15}{p['eff']:>9.3f}{p['sb']:>8.3f}"
          f"{p['maxmin']:>9.3f}{p['bias']:>7.2f}x")
