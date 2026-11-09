"""Compare standalone B+ Asimov fits across response-model configurations.

Reports sigma(A), sigma(e), sigma(M) -- averaged over the eta groups and per
group -- plus the A-e correlation, for each configuration side by side. Every
sigma is in PRIOR units, so sigma < 1 means the channel informs the parameter
and a SMALLER sigma means more constraint.

The direction of that last point is the trap this whole study exists around: a
model asked to extrapolate outside its training range produces a LARGER
template variation, which reads as a SMALLER sigma. An out-of-domain model
therefore looks like improved sensitivity in exactly this table. That is why
section 4's acceptance gate runs first and why a configuration that fails it is
not quoted here.

Convergence is checked, not assumed: a fit that writes output has not
necessarily converged.

Example:
    python scripts/btojpsik/kpi_solo_report.py \\
        muon_rec=.../cvhnano_full_..._kpi_muon_rec_fit.hdf5 \\
        kpi=.../cvhnano_full_..._kpi_kpi_fit.hdf5
"""

import argparse
import re
import sys

import numpy as np
import rabbit.io_tools as io

AEM = re.compile(r"^(?:scale)?([AeM])(\d+)$")
EDM_MAX = 1e-5


def read(path):
    fitres = io.get_fitresult(path)
    h = fitres["parms"].get()
    names = [n.decode() if isinstance(n, bytes) else str(n) for n in h.axes[0]]
    sig = {n: float(np.sqrt(v)) for n, v in zip(names, h.variances())}
    cov = None
    for key in ("cov", "covariance_matrix"):
        if key in fitres:
            cov = np.asarray(fitres[key].get().values())
            break
    edm = None
    for key in ("edmval", "edm"):
        if key in fitres:
            try:
                edm = float(np.asarray(fitres[key]).reshape(-1)[0])
            except Exception:
                edm = None
            break
    return names, sig, cov, edm


def group(sig):
    out = {"A": {}, "e": {}, "M": {}}
    for n, s in sig.items():
        m = AEM.match(n)
        if m:
            out[m.group(1)][int(m.group(2))] = s
    return out


def rho_Ae(names, cov):
    """Mean correlation between A_i and e_i across eta groups."""
    if cov is None:
        return np.nan
    idx = {}
    for i, n in enumerate(names):
        m = AEM.match(n)
        if m:
            idx.setdefault(m.group(1), {})[int(m.group(2))] = i
    if "A" not in idx or "e" not in idx:
        return np.nan
    vals = []
    for g in sorted(set(idx["A"]) & set(idx["e"])):
        ia, ie = idx["A"][g], idx["e"][g]
        d = np.sqrt(cov[ia, ia] * cov[ie, ie])
        if d > 0:
            vals.append(cov[ia, ie] / d)
    return float(np.mean(vals)) if vals else np.nan


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("inputs", nargs="+", help="label=fitresult.hdf5 ...")
    a = p.parse_args()

    labels, data = [], {}
    for spec in a.inputs:
        lab, path = spec.split("=", 1)
        labels.append(lab)
        try:
            data[lab] = read(path)
        except Exception as exc:
            print(f"{lab}: could not read ({type(exc).__name__}: {exc})")
            labels.remove(lab)

    if not labels:
        return 1

    w = 14
    print(f"{'':<12}" + "".join(f"{l:>{w}}" for l in labels))
    conv = {}
    for l in labels:
        edm = data[l][3]
        # Two-sided: -0.0 <= 1e-5 is True, and a failed fit can report -0.0.
        conv[l] = (edm is not None) and (0.0 <= edm <= EDM_MAX)
    print(f"{'edmval':<12}" + "".join(
        f"{('%.2e' % data[l][3]) if data[l][3] is not None else 'n/a':>{w}}"
        for l in labels))
    print(f"{'converged':<12}" + "".join(
        f"{('yes' if conv[l] else 'NO'):>{w}}" for l in labels))

    print()
    for par in ("A", "e", "M"):
        row = []
        for l in labels:
            g = group(data[l][1])[par]
            row.append(np.mean(list(g.values())) if g else np.nan)
        print(f"{'mean sigma(' + par + ')':<12}" +
              "".join(f"{v:>{w}.4f}" for v in row))
    print(f"{'rho(A,e)':<12}" +
          "".join(f"{rho_Ae(data[l][0], data[l][2]):>{w}.4f}" for l in labels))

    print(f"\n{'n nuisances':<12}" +
          "".join(f"{len(data[l][0]):>{w}}" for l in labels))

    base = labels[0]
    print(f"\nRatio to '{base}' (below 1 = tighter)")
    for par in ("A", "e", "M"):
        b = np.mean(list(group(data[base][1])[par].values()) or [np.nan])
        print(f"{'sigma(' + par + ')':<12}" + "".join(
            f"{np.mean(list(group(data[l][1])[par].values()) or [np.nan]) / b:>{w}.4f}"
            for l in labels))
    if not all(conv.values()):
        print("\nAT LEAST ONE FIT DID NOT CONVERGE -- its numbers are not results.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
