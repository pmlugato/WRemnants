#!/usr/bin/env python3
"""How much signal sits below the ONNX response network's 2 GeV gen-pT floor,
and whether the A/e/M variations behave there.

The selection floors the bachelor at reco pT 1 GeV; the response model was
trained on a sample floored at gen pT 2 GeV. `--genPtReweightSaturation 2.0`
in mode `condition` clamps only the network's conditioning input, leaving the
response residual and the physical shift at the true kinematics -- so the
question is not whether the network is extrapolated (it is not, it is clamped)
but whether clamping distorts the variation where it bites.
"""
import numpy as np
from wremnants.postprocessing.rabbit_cvhnano_helpers import open_results

H = "/ceph/submit/data/user/p/pmlugato/mz/histograms/"
GEN_SIGNAL = 0
FLOOR = 2.0
P = print


def vals(h):
    a = np.asarray(h.view())
    return (a["value"].astype(float) if a.dtype.names else a.astype(float))


# ---------------------------------------------------------------- gen pT
with open_results(H + "btojpsik_cvhnano_maxFiles_m1_gendiag_full.hdf5") as r:
    o = r["InclusiveBToJpsiX_2016PostVFP"]["output"]
    hg = o["genTruthKaonReco"].get()          # gen_truth_pt, gen_truth_eta, genCategory
    ge = np.asarray(hg.axes["gen_truth_pt"].edges)
    g = vals(hg)[..., GEN_SIGNAL].sum(axis=1)
    hall = o["genTruthKaon"].get()            # gen_truth_pt, eta, inAcc  (denominator)
    gall = vals(hall).sum(axis=1)

gc = 0.5 * (ge[:-1] + ge[1:])
below = gc < FLOOR
P(f"MARKER gen pT axis: {len(gc)} bins over [{ge[0]:.2f}, {ge[-1]:.2f}]")
P(f"MARKER selected signal, gen pT < {FLOOR}: {g[below].sum():,.0f} of {g.sum():,.0f}"
  f"  = {100*g[below].sum()/g.sum():.1f}%")
P(f"MARKER   ... produced (pre-selection) below the floor: "
  f"{100*gall[...,1][below].sum()/gall[...,1].sum():.1f}% of in-acceptance gen signal")
P("MARKER gen pT profile of SELECTED signal:")
for i in range(len(gc)):
    if g[i] == 0:
        continue
    flag = "  <-- below floor" if gc[i] < FLOOR else ""
    P(f"MARKER   [{ge[i]:5.2f},{ge[i+1]:5.2f})  {g[i]:>9,.0f}  {100*g[i]/g.sum():5.2f}%{flag}")

# ------------------------------------------------- variation size vs reco pT
with open_results(H + "btojpsik_cvhnano_fullstats.hdf5") as r:
    o = r["InclusiveBToJpsiX_2016PostVFP"]["output"]
    hv = o["nominal_muonScaleSyst_responseWeights"].get()
    hn = o["fitbins"].get()
    pe = np.asarray(hv.axes["fit_bachelor_pt"].edges)
    nom = vals(hn)[..., GEN_SIGNAL]                 # pt, eta, charge, mass
    var = vals(hv)[..., GEN_SIGNAL, :, :]           # pt, eta, charge, mass, unc, updn
    unc_n = len(hv.axes["unc"])

pc = 0.5 * (pe[:-1] + pe[1:])
mc = np.asarray(hn.axes["jointCvh_mass"].centers)


def mean_mass(w):
    s = w.sum()
    return float((w * mc).sum() / s) if s > 0 else np.nan


P("MARKER")
P("MARKER Peak SHIFT induced by each A/e/M variation, per reco-pT band [MeV]")
P(f"MARKER {'pT band':<12s}{'signal':>10s}{'share':>8s}" +
  "".join(f"{f'unc{u}':>9s}" for u in range(min(unc_n, 6))))
bands = [(1.0, 1.5), (1.5, 2.0), (2.0, 3.0), (3.0, 5.0), (5.0, 10.0)]
rows = []
for lo, hi in bands:
    m = (pc >= lo) & (pc < hi)
    if not m.any():
        continue
    n0 = nom[m].sum(axis=(0, 1, 2))
    if n0.sum() == 0:
        continue
    m0 = mean_mass(n0)
    cells = []
    for u in range(min(unc_n, 6)):
        up = var[m][..., u, 1].sum(axis=(0, 1, 2))
        dn = var[m][..., u, 0].sum(axis=(0, 1, 2))
        cells.append(0.5 * (mean_mass(up) - mean_mass(dn)) * 1e3)
    rows.append((lo, hi, n0.sum(), cells))
tot = sum(r[2] for r in rows)
for lo, hi, n, cells in rows:
    P(f"MARKER {f'{lo:.1f}-{hi:.1f}':<12s}{n:>10,.0f}{100*n/tot:>7.1f}%" +
      "".join(f"{c:>9.3f}" for c in cells))

P("MARKER")
P("MARKER Variation SIZE as a fraction of yield (max |up/nom - 1| over mass bins)")
P(f"MARKER {'pT band':<12s}" + "".join(f"{f'unc{u}':>9s}" for u in range(min(unc_n, 6))))
for lo, hi in bands:
    m = (pc >= lo) & (pc < hi)
    n0 = nom[m].sum(axis=(0, 1, 2))
    if n0.sum() == 0:
        continue
    cells = []
    for u in range(min(unc_n, 6)):
        up = var[m][..., u, 1].sum(axis=(0, 1, 2))
        ok = n0 > 0.05 * n0.max()
        cells.append(float(np.max(np.abs(up[ok] / n0[ok] - 1.0))) if ok.any() else np.nan)
    P(f"MARKER {f'{lo:.1f}-{hi:.1f}':<12s}" + "".join(f"{100*c:>8.2f}%" for c in cells))

# ------------------------------------------------ finiteness / throw guard
bad = ~np.isfinite(var)
P("MARKER")
P(f"MARKER non-finite variation entries: {int(bad.sum()):,} of {var.size:,}")
neg = (var < 0).sum()
P(f"MARKER negative variation entries:   {int(neg):,}")


def floor_vs_nofloor():
    """Clamped conditioning (floor 2.0) against the network extrapolated below it.

    `full_nofloor` is the same full-statistics run with
    --genPtReweightSaturation -1, so the ONLY difference is whether the
    network's conditioning input log pt_gen is clamped at 2 GeV.
    """
    out = {}
    for tag, fn in (("floor 2.0", "btojpsik_cvhnano_fullstats.hdf5"),
                    ("no floor", "btojpsik_cvhnano_full_nofloor.hdf5")):
        with open_results(H + fn) as r:
            o = r["InclusiveBToJpsiX_2016PostVFP"]["output"]
            hv = o["nominal_muonScaleSyst_responseWeights"].get()
            names = [a.name for a in hv.axes]
            v = vals(hv)
            gi = names.index("genCategory")
            v = np.take(v, GEN_SIGNAL, axis=gi)
            pe = np.asarray(hv.axes["fit_bachelor_pt"].edges)
            mi = [a.name for a in hv.axes].index("jointCvh_mass")
            out[tag] = (pe, v, names)
    (pe, vf, names), (_, vn, _) = out["floor 2.0"], out["no floor"]
    pc = 0.5 * (pe[:-1] + pe[1:])
    # axes after removing genCategory: pt, eta, charge, mass, unc, downUpVar
    mc_ = np.asarray(np.linspace(5.0, 5.5, vf.shape[3] + 1))
    mc_ = 0.5 * (mc_[:-1] + mc_[1:])

    def shift(w):
        s = w.sum()
        return float((w * mc_).sum() / s) if s > 0 else np.nan

    P("MARKER")
    P("MARKER Clamped conditioning vs the network extrapolated below 2 GeV")
    P("MARKER (half up-down peak shift in MeV; difference is the clamping bias)")
    P(f"MARKER {'pT band':<12s}{'unc':>5s}{'floor 2.0':>12s}{'no floor':>11s}{'diff':>9s}")
    for lo, hi in ((1.0, 1.5), (1.5, 2.0), (2.0, 3.0), (3.0, 5.0)):
        m = (pc >= lo) & (pc < hi)
        for u in (0, 1, 4):
            a = 0.5 * (shift(vf[m][..., u, 1].sum(axis=(0, 1, 2)))
                       - shift(vf[m][..., u, 0].sum(axis=(0, 1, 2)))) * 1e3
            b = 0.5 * (shift(vn[m][..., u, 1].sum(axis=(0, 1, 2)))
                       - shift(vn[m][..., u, 0].sum(axis=(0, 1, 2)))) * 1e3
            P(f"MARKER {f'{lo:.1f}-{hi:.1f}':<12s}{u:>5d}{a:>12.4f}{b:>11.4f}{b-a:>9.4f}")


floor_vs_nofloor()
