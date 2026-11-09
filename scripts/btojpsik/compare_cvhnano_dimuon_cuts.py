"""What the 2018-style dimuon, muon and B-vertex cuts do -- each alone, all together.

The question this exists to answer: the CVH NanoAOD has never been compared to
2018 on a level footing. 2018 imposed a tight dimuon and muon selection *and* a
1 GeV kaon pT floor at production, and only the kaon floor has been tested here.
If the dimuon cuts flatten the combinatorial on their own, then blaming the soft
kaons is premature -- and soft kaons are worth keeping, because they carry the
best B mass resolution and the lever that breaks the A/e correlation.

Per configuration this reports:

  * signal efficiency, and the efficiency **against bachelor pT**, because a cut
    that removes soft kaons preferentially is not acceptable however clean it
    looks inclusively;
  * peak-window purity;
  * the flatness of the surviving combinatorial -- its max/min inside the fit
    window and how many turning points it has;
  * what the background model then costs: chi2/ndf and fitted signal for the true
    simulated shape, an exponential and a quadratic.

    python scripts/btojpsik/compare_cvhnano_dimuon_cuts.py --tag baseline FILE \
        --tag "dimuon pT > 7" FILE ...
"""

import argparse
import json

import hist
import numpy as np

from wremnants.postprocessing import rabbit_cvhnano_helpers as H
from wremnants.production.btojpsik_cvhnano_axes import (
    GEN_DATA,
    GEN_MC_CATEGORIES,
    GEN_SIGNAL,
    GEN_TEMPLATE_CATEGORIES,
)

B_MASS = 5.27934
PEAK = (5.24, 5.32)
FIT = (5.15, 5.42)
PT_EDGES = [0.0, 0.5, 1.0, 1.5, 2.5, 4.0, 20.0]


class Tag(argparse.Action):
    def __call__(self, parser, ns, values, option_string=None):
        ns.inputs = getattr(ns, "inputs", None) or []
        ns.inputs.append(tuple(values))


def get(o, n):
    h = o[n]
    return h.get() if hasattr(h, "get") else h


def lsq(components, target, sigma):
    A = np.vstack([c / sigma for c in components]).T
    p, *_ = np.linalg.lstsq(A, target / sigma, rcond=None)
    model = sum(a * c for a, c in zip(p, components))
    return p, model, float(np.sum(((model - target) / sigma) ** 2))


def turning_points(y, smooth=3):
    if len(y) < smooth + 2:
        return 0
    k = np.ones(smooth) / smooth
    d = np.sign(np.diff(np.convolve(y, k, mode="valid")))
    d = d[d != 0]
    return int(np.sum(d[1:] * d[:-1] < 0)) if len(d) > 1 else 0


def measure(path):
    out = {}
    with H.open_results(path) as r:
        mo = r["InclusiveBToJpsiX_2016PostVFP"]["output"]
        do = r["Charmonium_2016PostVFP_cvhNano"]["output"]
        hm = get(mo, "nominal_BuJpsiK_jointCvhMass")
        hd = get(do, "nominal_BuJpsiK_jointCvhMass")
        c = hm.axes[0].centers
        peak = (c > PEAK[0]) & (c < PEAK[1])
        v = hm.values()
        sig = v[:, GEN_SIGNAL]
        tot = v[:, GEN_MC_CATEGORIES].sum(axis=1)
        nons = v[:, [k for k in GEN_MC_CATEGORIES if k != GEN_SIGNAL]].sum(axis=1)
        out["signal_all"] = float(sig.sum())
        out["signal_peak"] = float(sig[peak].sum())
        out["purity_peak"] = float(sig[peak].sum() / tot[peak].sum()) if tot[peak].sum() else 0.0
        realb = v[:, GEN_TEMPLATE_CATEGORIES].sum(axis=1)

        # flatness of the surviving combinatorial inside the fit window
        m = (c >= FIT[0]) & (c <= FIT[1])
        nb = nons[m]
        out["bkg_maxmin"] = float(nb.max() / max(nb.min(), 1e-9))
        out["bkg_turns"] = turning_points(nb)
        out["bkg_argmax"] = float(c[m][int(np.argmax(nb))])

        # background-model comparison on the data
        d = hd.values()[:, GEN_DATA][m]
        e = np.sqrt(np.maximum(d, 1.0))
        S = realb[m] / realb[m].sum()
        x = (c[m] - c[m][0]) / (c[m][-1] - c[m][0])
        p_t, _, chi_t = lsq([S, nb / nb.sum()], d, e)
        best = None
        for b in np.linspace(-30, 30, 601):
            comps = [S, np.exp(b * x)]
            pp, _, cc = lsq(comps, d, e)
            if np.any(pp < 0):
                continue
            if best is None or cc < best[0]:
                best = (cc, pp)
        p_q, _, chi_q = lsq([S, np.ones_like(x), x, x ** 2], d, e)
        n = int(m.sum())
        out["chi2_true"] = chi_t / (n - 2)
        out["chi2_exp"] = best[0] / (n - 3) if best else float("nan")
        out["chi2_quad"] = chi_q / (n - 4)
        out["yield_true"] = float(p_t[0] * S.sum())
        out["yield_exp"] = float(best[1][0] * S.sum()) if best else float("nan")
        out["bias_exp"] = out["yield_exp"] / out["yield_true"] if out["yield_true"] else np.nan
        out["data_window"] = float(d.sum())

        # signal against bachelor pT, for the bias check
        hf = get(mo, "fitbins")
        hs = hf[{"genCategory": hist.loc(GEN_SIGNAL)}].project(
            "fit_bachelor_pt", "jointCvh_mass")
        pt = hs.axes["fit_bachelor_pt"]
        mm = hs.axes["jointCvh_mass"].centers
        pk = (mm > PEAK[0]) & (mm < PEAK[1])
        per = []
        for lo, hi in zip(PT_EDGES[:-1], PT_EDGES[1:]):
            sel = (pt.centers >= lo) & (pt.centers < hi)
            per.append(float(hs.values()[sel][:, pk].sum()))
        out["signal_vs_pt"] = per
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tag", nargs=2, metavar=("NAME", "FILE"), action=Tag)
    p.add_argument("-o", "--output", default=None)
    args = p.parse_args()

    rows = []
    for name, path in args.inputs:
        try:
            m = measure(path)
        except Exception as exc:
            print(f"  skipping {name}: {exc}")
            continue
        m["name"] = name
        rows.append(m)
    base = rows[0]
    for m in rows:
        m["eff"] = m["signal_peak"] / base["signal_peak"]
        m["eff_vs_pt"] = [
            a / b if b > 0 else float("nan")
            for a, b in zip(m["signal_vs_pt"], base["signal_vs_pt"])
        ]

    hdr = (f"{'configuration':22s}{'sig eff':>9s}{'purity':>8s}{'bkg max/min':>12s}"
           f"{'turns':>6s}{'chi2 exp':>10s}{'chi2 quad':>10s}{'bias':>7s}")
    print(f"\nPeak window {PEAK}, fit window {FIT}. Reference is '{base['name']}'.\n")
    print(hdr)
    print("-" * len(hdr))
    for m in rows:
        print(f"{m['name']:22s}{m['eff']:9.3f}{m['purity_peak']:8.3f}"
              f"{m['bkg_maxmin']:12.2f}{m['bkg_turns']:6d}{m['chi2_exp']:10.1f}"
              f"{m['chi2_quad']:10.1f}{m['bias_exp']:6.2f}x")

    print(f"\nSignal efficiency against bachelor pT -- a flat row is an unbiased cut")
    print(f"{'configuration':22s}" + "".join(
        f"{f'{lo:g}-{hi:g}':>9s}" for lo, hi in zip(PT_EDGES[:-1], PT_EDGES[1:])))
    print("-" * (22 + 9 * (len(PT_EDGES) - 1)))
    for m in rows[1:]:
        print(f"{m['name']:22s}" + "".join(f"{e:9.3f}" for e in m["eff_vs_pt"]))

    if args.output:
        json.dump(rows, open(args.output, "w"), indent=1)
        print(f"\nwrote {args.output}")


if __name__ == "__main__":
    main()
