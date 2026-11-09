#!/usr/bin/env python3
"""Do the peaking-background shapes vary across fit cells?

This decides whether more inclusive simulation is needed. If a component's mass
shape is the same in every kaon-pT / eta slice, one cell-shared shape with a
floating normalisation is sufficient and the sample we have already supports it.
If the shape moves, each slice needs its own, and the statistics required go up
by the number of slices.

Per-bin template NOISE is not the criterion and is worth stating once: the ratio
of a template's per-bin noise to the data's own is sqrt(share x 18.7), which is
independent of binning -- splitting into cells divides template and data alike.
What breaks at fine binning is SPARSITY: a template bin with zero entries asserts
zero rate, and no bin-by-bin nuisance repairs that.
"""
import numpy as np
from wremnants.postprocessing.rabbit_cvhnano_helpers import open_results

H = "/ceph/submit/data/user/p/pmlugato/mz/histograms/"
LO, HI = 5.05, 5.45
COMPONENTS = {1: "B+ -> J/psi K X", 2: "B+ -> J/psi pi X", 3: "B0", 6: "combinatorial"}
P = print

with open_results(H + "btojpsik_cvhnano_fullstats.hdf5") as r:
    h = r["InclusiveBToJpsiX_2016PostVFP"]["output"]["fitbins"].get()
    pe = np.asarray(h.axes["fit_bachelor_pt"].edges)
    ee = np.asarray(h.axes["fit_bachelor_eta"].edges)
    me = np.asarray(h.axes["jointCvh_mass"].edges)
    a = np.asarray(h.view())
    y = a["value"].astype(float) if a.dtype.names else a.astype(float)

pc, ec, mc = (0.5 * (x[:-1] + x[1:]) for x in (pe, ee, me))
keep = (mc > LO) & (mc < HI)
y = y[..., keep, :]
mc = mc[keep]


def chi2_shape(u, v):
    """Poisson-ish chi2 between two raw histograms compared as SHAPES."""
    nu, nv = u.sum(), v.sum()
    if nu < 20 or nv < 20:
        return np.nan, 0
    k = nv / nu
    m = (u + v) > 0
    num = (v[m] - k * u[m]) ** 2
    den = v[m] + k * k * u[m]
    return float(np.sum(num[den > 0] / den[den > 0])), int(m.sum()) - 1


def split(axis_c, y, axis, thresh):
    lo = np.take(y, np.where(axis_c < thresh)[0], axis=axis).sum(axis=axis)
    hi = np.take(y, np.where(axis_c >= thresh)[0], axis=axis).sum(axis=axis)
    return lo, hi


P(f"MARKER Shape stability of each component, window {LO}-{HI}, 80 bins")
P(f"MARKER {'component':<20s}{'split':<22s}{'N low':>9s}{'N high':>9s}"
  f"{'chi2':>9s}{'ndf':>5s}{'chi2/ndf':>10s}{'<m> shift [MeV]':>17s}")
for ic, name in COMPONENTS.items():
    comp = y[..., ic]                                  # pt, eta, charge, mass
    for label, arr, ax, thr in (
        ("kaon pT  1.6 GeV", comp, 0, 1.6),
        ("kaon pT  3.0 GeV", comp, 0, 3.0),
        ("kaon |eta| 1.0", comp, 1, None),
    ):
        if thr is None:
            idx_lo = np.where(np.abs(ec) < 1.0)[0]
            idx_hi = np.where(np.abs(ec) >= 1.0)[0]
            lo = np.take(arr, idx_lo, axis=1)
            hi = np.take(arr, idx_hi, axis=1)
            lo = lo.sum(axis=(0, 1, 2))
            hi = hi.sum(axis=(0, 1, 2))
        else:
            lo, hi = split(pc, arr, 0, thr)
            lo = lo.sum(axis=(0, 1))
            hi = hi.sum(axis=(0, 1))
        c2, ndf = chi2_shape(lo, hi)
        if np.isnan(c2):
            P(f"MARKER {name:<20s}{label:<22s}{lo.sum():>9,.0f}{hi.sum():>9,.0f}"
              f"{'--':>9s}{'--':>5s}{'too few':>10s}")
            continue
        ml = float((lo * mc).sum() / lo.sum())
        mh = float((hi * mc).sum() / hi.sum())
        P(f"MARKER {name:<20s}{label:<22s}{lo.sum():>9,.0f}{hi.sum():>9,.0f}"
          f"{c2:>9.1f}{ndf:>5d}{c2/ndf:>10.2f}{1000*(mh-ml):>17.1f}")

P("MARKER")
P("MARKER Events per mass bin at each granularity (need >~15/bin for a usable shape)")
P(f"MARKER {'component':<20s}{'inclusive':>11s}{'3 pT bands':>12s}{'3pT x 2q':>10s}{'288 cells':>11s}")
for ic, name in COMPONENTS.items():
    n = y[..., ic].sum()
    P(f"MARKER {name:<20s}{n/80:>11,.0f}{n/3/80:>12,.0f}{n/6/80:>10,.0f}{n/288/80:>11,.2f}")
