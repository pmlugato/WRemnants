"""Compare what a selection does to the background shape the fit has to model.

The figure of merit is not purity. It is whether the surviving background inside
the fit window can be described by the models rabbit offers, both of which are
monotonic: the maximum-to-minimum ratio, the number of turning points, and the
chi-square of an exponential and a straight line fitted to it. Signal efficiency
is reported as the cost, to be minimised among selections that pass.

Usage:
  compare_cvhnano_selections.py --label ref FILE --label sel2018 FILE [...]
"""

import argparse

import numpy as np
from scipy.optimize import least_squares

from wremnants.postprocessing import rabbit_cvhnano_helpers as H
from wremnants.production.btojpsik_cvhnano_axes import (
    GEN_DATA,
    GEN_NON_SIGNAL_CATEGORIES,
    GEN_SIGNAL,
)

WINDOW = (5.20, 5.35)
NBINS = 10


class LabelledFile(argparse.Action):
    def __call__(self, parser, ns, values, option_string=None):
        ns.inputs = getattr(ns, "inputs", None) or []
        ns.inputs.append(tuple(values))


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--label", nargs=2, metavar=("NAME", "FILE"), action=LabelledFile,
                    help="repeatable: a name and a histmaker output")
args = parser.parse_args()


def turning_points(y, smooth=3):
    """Sign changes of the smoothed derivative -- structure, not statistical flips."""
    if len(y) < smooth + 2:
        return 0
    k = np.ones(smooth) / smooth
    d = np.sign(np.diff(np.convolve(y, k, mode="valid")))
    d = d[d != 0]
    return int(np.sum(d[1:] * d[:-1] < 0)) if len(d) > 1 else 0


MODELS = (
    ("exp", lambda p, t: np.exp(p[0] + p[1] * t), lambda n: [np.log(n.mean()), -0.5]),
    ("lin", lambda p, t: p[0] * (1 - t) + p[1] * t, lambda n: [n[0], n[-1]]),
)


def fit_monotonic(y, err):
    """Fit the two monotonic families; return chi2/ndf and a scale-free misfit.

    chi2 alone cannot be compared between selections, because a tighter selection
    has fewer events, larger relative errors and therefore a smaller chi2 for the
    *same* shape mismatch -- which would read as an improvement that is really just
    a loss of statistics. `rms` is the root-mean-square fractional residual of the
    best fit, which depends on the shape and not on how many events measured it.
    """
    x = np.linspace(0.0, 1.0, len(y))
    n = y / y.sum()
    e = np.where(err > 0, err / y.sum(), n.max() * 1e-3)
    out = {}
    for name, model, seed in MODELS:
        r = least_squares(lambda p: (model(p, x) - n) / e, seed(n))
        fit = model(r.x, x)
        out[f"c2_{name}"] = float(np.sum(((fit - n) / e) ** 2)) / (len(n) - len(r.x))
        out[f"rms_{name}"] = float(np.sqrt(np.mean(((n - fit) / fit) ** 2)))
    return out


rows = []
for name, path in args.inputs:
    d_raw, m_raw, _ = H.load_fitbins(path, "Charmonium", "InclusiveBToJpsiX")
    # Inclusive in the cell axes: this compares mass shapes, so pT and eta are
    # summed over and cannot contribute a statistics argument.
    kw = dict(mass_low=WINDOW[0], mass_high=WINDOW[1], mass_bins=NBINS,
              signal_sigma=0.0138, pt_edges=[0.0, 10.0], eta_rebin=48,
              collapse_axes=(H.PT_AXIS, H.ETA_AXIS, H.CHARGE_AXIS))
    m = H.reduce_fitbins(m_raw, **kw)
    d = H.reduce_fitbins(d_raw, **kw)

    S = H.select_category(m, GEN_SIGNAL)
    O = H.sum_categories(m, GEN_NON_SIGNAL_CATEGORIES)
    D = H.select_category(d, GEN_DATA)
    im = list(S.axes.name).index(H.MASS_AXIS)
    other = tuple(i for i in range(S.values().ndim) if i != im)

    s = S.values().sum(axis=other)
    o = O.values().sum(axis=other)
    ov = O.variances().sum(axis=other)
    dd = D.values().sum(axis=other)
    centres = S.axes[H.MASS_AXIS].centers

    # The figure of merit that actually decides the measurement: how far the
    # fitted signal yield moves when the data is fitted with a monotonic
    # background instead of the true simulated shape. It combines the two things
    # that matter -- how curved the background is, and how much of the window it
    # occupies -- so a large flat background and a small curved one both score
    # well, which is the correct behaviour and is not true of S/B or of max/min.
    xx = np.linspace(0.0, 1.0, len(dd))
    w = np.sqrt(np.maximum(dd, 1.0))
    r_true = least_squares(
        lambda p: (p[0] * s + p[1] * o - dd) / w, [dd.sum() / max(s.sum(), 1), 1.0])
    r_exp = least_squares(
        lambda p: (p[0] * s + np.exp(p[1] + p[2] * xx) * dd.mean() - dd) / w,
        [dd.sum() / max(s.sum(), 1), 0.0, -0.5])
    y_true = float(r_true.x[0] * s.sum())
    y_exp = float(r_exp.x[0] * s.sum())
    bias = y_exp / y_true if y_true > 0 else np.nan

    n = o / o.sum()
    nerr = np.sqrt(ov) / o.sum()
    row = dict(
        yield_true=y_true, yield_exp=y_exp, bias=bias,
        name=name, signal=s.sum(), bkg=o.sum(), data=dd.sum(),
        maxmin=float(n.max() / n.min()),
        # Uncertainty on max/min, so a change can be told from noise.
        maxmin_err=float((n.max() / n.min()) * np.hypot(
            nerr[np.argmax(n)] / n.max(), nerr[np.argmin(n)] / n.min())),
        argmax=float(centres[np.argmax(n)]),
        turns=turning_points(o), sb=float(s.sum() / o.sum()),
    )
    row.update(fit_monotonic(o, np.sqrt(ov)))
    rows.append(row)

ref = rows[0]
print(f"Fit window [{WINDOW[0]}, {WINDOW[1]}], {NBINS} bins. "
      f"Reference is '{ref['name']}'.\n")
hdr = (f"{'selection':<13}{'sigeff':>8}{'S/B':>8}{'purity':>8}"
       f"{'max/min':>15}{'turns':>6}{'rms':>8}{'yield bias':>12}")
print(hdr)
print("-" * len(hdr))
for r in rows:
    mm = f"{r['maxmin']:.3f}+-{r['maxmin_err']:.3f}"
    purity = r["signal"] / (r["signal"] + r["bkg"])
    print(f"{r['name']:<13}{r['signal']/ref['signal']:>8.3f}{r['sb']:>8.4f}{purity:>8.3f}"
          f"{mm:>15}{r['turns']:>6d}{r['rms_exp']:>8.4f}{r['bias']:>11.2f}x")

print("\n'yield bias' is the ratio of the signal yield fitted with an exponential")
print("background to the one fitted with the true simulated shape. It is the number")
print("that decides the measurement: it folds together how curved the background is")
print("AND how much of the window it occupies, so a large flat background and a small")
print("curved one both score well. S/B and max/min each capture only one half.\n")
print("Verdict:")
for r in rows[1:]:
    dmm = ref["maxmin"] - r["maxmin"]
    sig = dmm / np.hypot(ref["maxmin_err"], r["maxmin_err"])
    print(f"  {r['name']}: signal eff {r['signal']/ref['signal']:.3f}, "
          f"S/B {ref['sb']:.3f} -> {r['sb']:.3f}, "
          f"max/min {ref['maxmin']:.3f} -> {r['maxmin']:.3f} ({sig:+.1f} sigma), "
          f"yield bias {ref['bias']:.2f}x -> {r['bias']:.2f}x")
    better_bias = abs(r["bias"] - 1) < abs(ref["bias"] - 1) - 0.05
    worse_bias = abs(r["bias"] - 1) > abs(ref["bias"] - 1) + 0.05
    if r["turns"] < ref["turns"]:
        print("    ACCEPT: turning point removed, so a monotonic model can describe it.")
    elif better_bias:
        print("    IMPROVES the bias, though the background still has a maximum. "
              "Worth keeping\n            if the signal efficiency is affordable.")
    elif worse_bias:
        print("    REJECT: makes the yield bias WORSE. These cuts select B-like "
              "candidates, and\n            a background of B-like candidates peaks at the "
              "B mass too -- so purity\n            and background flatness pull against "
              "each other here.")
    else:
        print("    REJECT: bias unchanged. Signal efficiency was paid for nothing, since "
              "a large\n            flat background is fittable and this one is no flatter.")
