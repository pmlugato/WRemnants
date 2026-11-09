#!/usr/bin/env python3
"""What does adding B+ to the joint Z + J/psi Asimov fit buy?

Compares the combined fit against the Z + J/psi reference, parameter by
parameter. Everything here is Asimov, so it is an UPPER BOUND on the real gain:
a wrong model fits itself perfectly, and B+'s data fit does not yet converge.
"""

import argparse
import re
import sys

import numpy as np

import rabbit.io_tools as io

AEM = re.compile(r"^([AeM])(\d+)$")


def parms(path):
    fitres = io.get_fitresult(path)
    h = fitres["parms"].get()
    names = [n.decode() if isinstance(n, bytes) else str(n) for n in h.axes[0]]
    return {
        n: (float(v), float(np.sqrt(e)))
        for n, v, e in zip(names, h.values(), h.variances())
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("fit")
    p.add_argument("--reference", required=True)
    p.add_argument("--nEta", type=int, default=24)
    a = p.parse_args()

    new, ref = parms(a.fit), parms(a.reference)
    print(f"combined : {len(new)} parameters   {a.fit.split('/')[-1]}")
    print(f"reference: {len(ref)} parameters   {a.reference.split('/')[-1]}")

    print("\n--- A/e/M, sigma in prior units. sigma<1 means the data informs it.")
    print(
        f"{'param':<7s}{'eta':>4s}{'ref':>9s}{'with B+':>9s}{'ratio':>8s}"
        f"{'informed by':>22s}"
    )
    gains = {"A": [], "e": [], "M": []}
    for param in ("A", "e", "M"):
        for i in range(a.nEta):
            n = f"{param}{i}"
            if n not in new or n not in ref:
                continue
            r, w = ref[n][1], new[n][1]
            gains[param].append(w / r if r else np.nan)
            tag = ""
            if w < 1.0 and r >= 1.0:
                tag = "B+ (newly)"
            elif w < 1.0:
                tag = "both"
            elif r < 1.0:
                tag = "reference only"
            else:
                tag = "neither"
            if i < 4 or w / r < 0.97:
                print(f"{param:<7s}{i:>4d}{r:>9.3f}{w:>9.3f}{w/r:>8.3f}{tag:>22s}")
    print("\n--- summary per family (mean sigma ratio, combined / reference)")
    for param, g in gains.items():
        g = np.asarray([x for x in g if np.isfinite(x)])
        if not len(g):
            continue
        nimp = sum(
            1
            for i in range(a.nEta)
            if f"{param}{i}" in new and new[f"{param}{i}"][1] < 1.0
        )
        print(
            f"  {param}: mean ratio {g.mean():.3f}   best {g.min():.3f}   "
            f"eta bins with sigma<1: {nimp}/{a.nEta}"
        )

    # --- the figure of merit is NOT the marginal A/e/M errors. The reference
    # already constrains them to ~0.1 of the prior; B+ is far weaker in absolute
    # terms. Its claim is the LEVER ARM -- k = 1/pT reaching 1.0 where J/psi
    # stops at 0.33 -- so what has to move is the A-e degeneracy and, through it,
    # the calibration's impact on mZ.
    def corr_and_mz(path, tag):
        fr = io.get_fitresult(path)
        hc = fr["cov"].get()
        nm = [x.decode() if isinstance(x, bytes) else str(x) for x in hc.axes[0]]
        idx = {n: i for i, n in enumerate(nm)}
        C = hc.values()
        out = {}
        rho = []
        for i in range(a.nEta):
            ia, ie = idx.get(f"A{i}"), idx.get(f"e{i}")
            if ia is None or ie is None:
                continue
            den = np.sqrt(C[ia, ia] * C[ie, ie])
            if den > 0:
                rho.append(C[ia, ie] / den)
        out["rho_Ae"] = np.asarray(rho)
        gi = fr.get("impacts_grouped")
        out["impacts"] = {}
        if gi is not None:
            hg = gi.get()
            # The 'parms' axis is NOT ordered with the POI first -- row 0 is a
            # J/psi normalisation. Select the POI BY NAME; taking row 0 reports
            # the impacts on whatever happened to be listed first, which is how
            # the first version of this script produced a table showing every
            # systematic falling 40x.
            pax = [
                x.decode() if isinstance(x, bytes) else str(x) for x in hg.axes["parms"]
            ]
            gax = [
                x.decode() if isinstance(x, bytes) else str(x)
                for x in hg.axes["impacts"]
            ]
            poi = next(
                (n for n in ("massShiftZ100MeV", "zmass", "massShiftZ") if n in pax),
                None,
            )
            if poi is None:
                print(
                    f"  WARN: no mZ POI on the impacts 'parms' axis; "
                    f"saw e.g. {pax[:3]}"
                )
            else:
                row = np.asarray(hg.values())[pax.index(poi)]
                out["impacts"] = dict(zip(gax, row))
                out["poi"] = poi
        print(
            f"\n[{tag}] rho(A_i, e_i): mean {out['rho_Ae'].mean():+.3f}  "
            f"min {out['rho_Ae'].min():+.3f}  max {out['rho_Ae'].max():+.3f}  "
            f"(n={len(out['rho_Ae'])})"
        )
        return out

    cr = corr_and_mz(a.reference, "reference  Z + J/psi")
    cn = corr_and_mz(a.fit, "combined   Z + J/psi + B+")
    if len(cr["rho_Ae"]) == len(cn["rho_Ae"]) and len(cr["rho_Ae"]):
        d = cn["rho_Ae"] - cr["rho_Ae"]
        print(
            f"  change in rho(A,e): mean {d.mean():+.4f}  "
            f"largest reduction {d.min():+.4f}"
        )
    if cr.get("impacts") and cn.get("impacts"):
        print(
            "\n--- grouped impacts on "
            + str(cr.get("poi", "the POI"))
            + "  [x100 for MeV]"
        )
        keys = sorted(set(cr["impacts"]) | set(cn["impacts"]))
        print(f"  {'group':<28s}{'reference':>11s}{'with B+':>11s}{'ratio':>8s}")
        for k in keys:
            r, w = cr["impacts"].get(k, np.nan), cn["impacts"].get(k, np.nan)
            if not (np.isfinite(r) and np.isfinite(w)):
                continue
            if abs(r) < 1e-6 and abs(w) < 1e-6:
                continue
            print(
                f"  {k:<28s}{r:>11.4f}{w:>11.4f}"
                + (f"{w/r:>8.3f}" if r else f"{'--':>8s}")
            )

    for poi in ("massShiftZ100MeV", "zmass", "massShiftZ"):
        if poi in new and poi in ref:
            print(
                f"\n--- {poi}: reference {ref[poi][1]:.5f}  "
                f"with B+ {new[poi][1]:.5f}  ratio {new[poi][1]/ref[poi][1]:.4f}"
            )
            break
    else:
        cands = [n for n in new if "mass" in n.lower() or "zmass" in n.lower()]
        print(f"\n(no standard mZ POI name matched; candidates: {cands[:6]})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
