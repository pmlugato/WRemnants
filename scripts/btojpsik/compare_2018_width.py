"""Compare the 2018 BMM-NanoAOD peak with the CVH NanoAOD peak, on one method.

The question this answers: in 2018 the data peak visibly matched the exclusive
BuToJpsiK simulation before any fitting, and in the CVH NanoAOD it appears not
to. Two explanations compete, and they are distinguishable.

  (a) **Purity.** 2018 was ~85% signal in the window, so its data peak was
      nearly background-free and its width could be read straight off. The CVH
      NanoAOD is ~21% signal, so the same naive read is contaminated by a broad
      combinatorial sitting under the peak, and looks wider for that reason
      alone.
  (b) **The two discriminants.** `softMuonMva` and the BMM BDT have no
      counterpart in our AlCaReco. If they are what produced the agreement,
      dropping them from the 2018 selection must break it.

So the 2018 histmaker is run twice, with and without those two cuts, and every
number below is computed the same way on every input: a Gaussian on a linear
term, width scanned, over one common window around the B mass.

The comparison that matters is **data against the simulated TOTAL**, never
against the signal component alone -- comparing a spectrum that contains
background to a template that does not is what made the first version of this
measurement report a spurious factor 1.71.

    python scripts/btojpsik/compare_2018_width.py \
        --label "2018, full selection" FILE --label "2018, no BDT/softMVA" FILE \
        --label "CVH NanoAOD" FILE
"""

import argparse
import json

import numpy as np

from wremnants.postprocessing import rabbit_cvhnano_helpers as H

B_MASS = 5.27934
# Half-width of the common fit window. 74 MeV, because the 2018 mass axis stops
# at 5.2 and a window symmetric about the B mass cannot be wider than that.
HALF = 0.074

# (hist, axis, data key, simulated keys, signal-only key) per input flavour.
LAYOUTS = {
    "cvhnano": dict(
        hist="nominal_BuJpsiK_jointCvhMass",
        data="Charmonium_2016PostVFP_cvhNano",
        mc=["InclusiveBToJpsiX_2016PostVFP"],
        gen_axis=True,
    ),
    "bmm2018": dict(
        hist="nominal_bkmm_jpsimc_mass",
        data="Charmonium_2018",
        mc=["BuToJpsiK_2018"],
        gen_axis=False,
    ),
}


class Labelled(argparse.Action):
    def __call__(self, parser, ns, values, option_string=None):
        ns.inputs = getattr(ns, "inputs", None) or []
        ns.inputs.append(tuple(values))


def gauss_width(c, y, centre=B_MASS):
    best = None
    for sigma in np.linspace(0.004, 0.060, 449):
        g = np.exp(-0.5 * ((c - centre) / sigma) ** 2)
        A = np.vstack([g, np.ones_like(c), c - centre]).T
        p, *_ = np.linalg.lstsq(A, y, rcond=None)
        res = float(np.sum((A @ p - y) ** 2))
        if best is None or res < best[0]:
            best = (res, sigma, p)
    return best[1], best[2]


def get(o, name):
    h = o[name]
    return h.get() if hasattr(h, "get") else h


def spectra(path, layout):
    """Data and simulated-total mass spectra, plus the signal component."""
    L = LAYOUTS[layout]
    with H.open_results(path) as r:
        hd = get(r[L["data"]]["output"], L["hist"])
        c = hd.axes[0].centers
        if L["gen_axis"]:
            from wremnants.production.btojpsik_cvhnano_axes import (
                GEN_DATA,
                GEN_MC_CATEGORIES,
                GEN_SIGNAL,
            )

            d = hd.values()[:, GEN_DATA]
            hm = get(r[L["mc"][0]]["output"], L["hist"])
            total = hm.values()[:, GEN_MC_CATEGORIES].sum(axis=1)
            sig = hm.values()[:, GEN_SIGNAL]
        else:
            d = hd.values()
            total = np.zeros_like(d)
            for key in L["mc"]:
                total = total + get(r[key]["output"], L["hist"]).values()
            # The 2018 simulation IS exclusive signal, so the total is the signal.
            sig = total.copy()
    return c, d, total, sig


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--label", nargs=2, metavar=("NAME", "FILE"), action=Labelled,
                   help="repeatable")
    p.add_argument("--layout", nargs=2, metavar=("NAME", "LAYOUT"), action="append",
                   default=[], help="repeatable: name -> cvhnano | bmm2018")
    p.add_argument("-o", "--output", default=None)
    args = p.parse_args()
    layouts = dict(args.layout)

    rows = []
    for name, path in args.inputs:
        layout = layouts.get(name, "bmm2018" if "2018" in name else "cvhnano")
        c, d, total, sig = spectra(path, layout)
        w = np.abs(c - B_MASS) < HALF
        if w.sum() < 8 or d[w].sum() <= 0:
            print(f"  skipping {name}: {w.sum()} bins, {d[w].sum():.0f} data")
            continue
        sd, _ = gauss_width(c[w], d[w])
        st, _ = gauss_width(c[w], total[w])
        ss, _ = gauss_width(c[w], sig[w])
        # Simulated purity is only meaningful where the simulation models the
        # background. The 2018 sample is EXCLUSIVE BuToJpsiK, so sig/total is 1.0
        # by construction and reporting it would look like a 100%-pure dataset.
        # Its real in-window purity, ~85%, comes from the normalised stack plot
        # (data 355 869, BuToJpsiK 301 554) and cannot be recovered from an
        # unnormalised histmaker output.
        tight = np.abs(c - B_MASS) < 0.04
        purity = (
            float(sig[tight].sum() / total[tight].sum())
            if layout == "cvhnano" and total[tight].sum()
            else np.nan
        )
        rows.append(dict(
            name=name, layout=layout,
            data=float(d[w].sum()), mc_total=float(total[w].sum()),
            data_sigma_mev=sd * 1000.0, mc_total_sigma_mev=st * 1000.0,
            mc_signal_sigma_mev=ss * 1000.0,
            ratio_to_total=sd / st if st > 0 else np.nan,
            ratio_to_signal=sd / ss if ss > 0 else np.nan,
            simulated_purity=purity,
        ))

    hdr = (f"{'input':30s}{'data':>11s}{'data σ':>9s}{'MCtot σ':>9s}"
           f"{'σ ratio':>9s}{'sim purity':>12s}")
    print(f"\nCommon window |m - {B_MASS}| < {HALF*1000:.0f} MeV, one fit method.\n")
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        pur = ("n/a" if np.isnan(r["simulated_purity"])
               else f"{r['simulated_purity']:.3f}")
        print(f"{r['name']:30s}{r['data']:11.0f}{r['data_sigma_mev']:9.2f}"
              f"{r['mc_total_sigma_mev']:9.2f}{r['ratio_to_total']:9.3f}{pur:>12s}")
    print("\n'σ ratio' is data / simulated total. It is the only comparison that is")
    print("fair when purity differs between samples: both sides then contain the")
    print("same background, so a shared inflation cancels. For the 2018 inputs the")
    print("simulation is EXCLUSIVE signal with no background at all, so their ratio")
    print("is an UPPER bound on a fair comparison -- the data side carries background")
    print("the simulated side does not, which can only inflate it.")
    if args.output:
        json.dump(rows, open(args.output, "w"), indent=1)
        print(f"\nwrote {args.output}")


if __name__ == "__main__":
    main()
