#!/usr/bin/env python3
"""Does extending the bachelor down to 0.1 GeV buy A/e separation?

Three Asimov fits, identical except for the bachelor pT floor and the response
network's conditioning clamp:

  pT > 1.0, clamped at 2 GeV   -- the configuration everything else used
  pT > 0.1, clamped at 2 GeV   -- extended range, network conditioning frozen
  pT > 0.1, clamp OFF          -- extended range, network extrapolated

Below 2 GeV the reported sensitivity is a property of the clamp rather than of
the detector, so the SPREAD between the last two is the answer's uncertainty. A
conclusion that does not survive both is not a conclusion.
"""
import sys

import numpy as np
import rabbit.io_tools as io

FITS = "/ceph/submit/data/user/p/pmlugato/mz/fitresults/"
CFGS = [
    ("pT > 1.0, clamped", "cvhnano_full_exppoly_asimov_eta12_normgrouped_aem_n2_pt1clamp_fit.hdf5"),
    ("pT > 0.1, clamped", "cvhnano_full_exppoly_asimov_eta12_normgrouped_aem_n2_softclamp_fit.hdf5"),
    ("pT > 0.1, clamp OFF", "cvhnano_full_exppoly_asimov_eta12_normgrouped_aem_n2_softnoclamp_fit.hdf5"),
]


def load(path):
    fr = io.get_fitresult(FITS + path)
    h = fr["parms"].get()
    nm = [x.decode() if isinstance(x, bytes) else str(x) for x in h.axes[0]]
    sig = {n: float(np.sqrt(v)) for n, v in zip(nm, h.variances())}
    e = fr["edmval"]
    try:
        e = float(np.asarray(e))
    except Exception:
        e = float(np.asarray(e.get()))
    hc = fr["cov"].get()
    cn = [x.decode() if isinstance(x, bytes) else str(x) for x in hc.axes[0]]
    idx = {n: i for i, n in enumerate(cn)}
    C = hc.values()
    rho = []
    for i in range(4):
        ia, ie = idx.get(f"scaleA{i}"), idx.get(f"scalee{i}")
        if ia is None or ie is None:
            continue
        d = np.sqrt(C[ia, ia] * C[ie, ie])
        if d > 0:
            rho.append(C[ia, ie] / d)
    return sig, e, np.asarray(rho)


rows = []
for lab, f in CFGS:
    try:
        rows.append((lab,) + load(f))
    except Exception as exc:
        print(f"MARKER {lab:<22s}  MISSING ({type(exc).__name__})")

print(f"MARKER {'configuration':<22s}{'edmval':>11s}{'sigma A':>9s}{'sigma e':>9s}"
      f"{'sigma M':>9s}{'rho(A,e)':>10s}")
base = None
for lab, sig, edm, rho in rows:
    g = lambda p: np.mean([sig[f"scale{p}{i}"] for i in range(4)
                           if f"scale{p}{i}" in sig])
    a, e_, m = g("A"), g("e"), g("M")
    ok = "" if (0 <= edm <= 1e-5) else "  NOT CONVERGED"
    print(f"MARKER {lab:<22s}{edm:>11.2e}{a:>9.3f}{e_:>9.3f}{m:>9.3f}"
          f"{np.mean(rho):>10.3f}{ok}")
    if base is None:
        base = (a, e_, m, np.mean(rho))
print("MARKER")
print(f"MARKER relative to pT > 1.0 clamped:")
for lab, sig, edm, rho in rows[1:]:
    g = lambda p: np.mean([sig[f"scale{p}{i}"] for i in range(4)
                           if f"scale{p}{i}" in sig])
    a, e_, m = g("A"), g("e"), g("M")
    print(f"MARKER   {lab:<22s} sigma A x{a/base[0]:.3f}   sigma e x{e_/base[1]:.3f}"
          f"   sigma M x{m/base[2]:.3f}   rho {np.mean(rho)-base[3]:+.3f}")
