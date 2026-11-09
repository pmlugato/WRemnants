#!/usr/bin/env python3
"""Per-leg energy-loss mass rule: regression test for the merge that unified
our ``particle_masses`` path with havyn's ``reweight_mass`` path.

Three properties, each of which a plausible resolution of that merge breaks:

1. An empty mass list and an all-zero mass list give **bit-identical**
   variations. This is what says the massless legs kept the behaviour they had
   before the mass-aware overload existed -- the 6-argument overload has no
   ``1/cosh(eta)`` on the e term and the 7-argument one does, so routing a
   zero-mass leg to the wrong overload changes its answer by up to 5.6x.
2. ``[m_K, 0, 0]`` differs from all-zero only through leg 0.
3. ``[m_K, 0, 0]`` differs from ``[m_K, m_K, m_K]``. This is the broadcast trap:
   a one-element ``[m_K]`` means "all legs" under the unified rule, so an
   implicit list would silently give both muons a 494 MeV energy-loss term.
   In the pt convention this difference is *small in size* -- the muon legs here
   carry 21 and 32 GeV of total momentum, where E/p differs from 1 by 3e-4 -- so
   the assertion is that it is non-zero, not that it is large. Its smallness is
   itself a check: the mass enters only through E/p, so it must vanish for a hard
   track.

Run inside the container, with setup.sh sourced:

    python scripts/tests/test_perleg_mass_rule.py
"""

import hist
import numpy as np
import ROOT

import narf
import wremnants.production.muon_calibration  # noqa: F401  (declares the headers)
from wremnants.utilities import binning

M_K = 0.493677
NETA = 4
NPARM = 3


def _make_helper(masses):
    """A JpsiCorrectionsUncHelperSplines over NETA eta bins with one variation
    per (eta, parameter) and a known uncertainty in each."""
    nvars = NETA * NPARM
    axis_eta = hist.axis.Regular(NETA, -2.4, 2.4, name="scale_eta")
    axis_parm = hist.axis.Integer(
        0, NPARM, underflow=False, overflow=False, name="scale_params"
    )
    axis_unc = hist.axis.Integer(0, nvars, underflow=False, overflow=False, name="unc")

    # A = 1e-4, e = 1e-2 GeV, M = 1e-5: the same order as the real correction
    # file, so the test exercises a realistic shift rather than a huge one.
    uncs = np.array([1e-4, 1e-2, 1e-5])
    var_mat = np.zeros((NETA, NPARM, nvars))
    for ieta in range(NETA):
        for iparm in range(NPARM):
            var_mat[ieta, iparm, ieta * NPARM + iparm] = uncs[iparm]

    h = hist.Hist(axis_eta, axis_parm, axis_unc)
    h[...] = var_mat
    h_cpp = narf.hist_to_pyroot_boost(h, tensor_rank=2)

    vec = ROOT.std.vector["double"]()
    for m in masses:
        vec.push_back(float(m))
    helper = ROOT.wrem.JpsiCorrectionsUncHelperSplines[type(h_cpp).__cpp_name__](
        ROOT.std.move(h_cpp), vec
    )
    helper.tensor_axes = (h.axes["unc"], binning.down_up_axis)
    return helper


# One candidate: a 1.5 GeV bachelor at eta = 0.8 and two muons of 6 and 4 GeV.
# The muon etas are deliberately large, where the 1/cosh(eta) convention
# difference is biggest and property 1 has the most to catch.
REC_PT = [1.5, 6.0, 4.0]
REC_ETA = [0.8, -1.9, 2.1]
REC_Q = [1, 1, -1]


def _evaluate(masses):
    helper = _make_helper(masses)
    rec_pt = ROOT.VecOps.RVec["float"]()
    rec_eta = ROOT.VecOps.RVec["float"]()
    rec_q = ROOT.VecOps.RVec["int"]()
    for pt, eta, q in zip(REC_PT, REC_ETA, REC_Q):
        rec_pt.push_back(pt)
        rec_eta.push_back(eta)
        rec_q.push_back(q)
    # gen == reco, and a unit response derivative per leg, so the returned
    # tensor is 1 + delta_qop and isolates the mass handling under test.
    resp = ROOT.VecOps.RVec["std::pair<double,double>"]()
    for _ in REC_PT:
        resp.push_back(ROOT.std.pair["double", "double"](1.0, 0.0))
    out = helper(rec_pt, rec_eta, rec_q, rec_pt, rec_eta, rec_q, resp, 1.0)
    # Eigen::TensorFixedSize does not convert to a numpy array through the
    # buffer protocol, so read it element by element.
    nvars = NETA * NPARM
    return np.array(
        [[out(ivar, idownup) for idownup in range(2)] for ivar in range(nvars)]
    )


def main():
    empty = _evaluate([])
    allzero = _evaluate([0.0, 0.0, 0.0])
    kaon_only = _evaluate([M_K, 0.0, 0.0])
    broadcast = _evaluate([M_K, M_K, M_K])

    failures = []

    # 1. massless is massless, exactly.
    if not np.array_equal(empty, allzero):
        worst = np.max(np.abs(empty - allzero))
        failures.append(
            f"empty mass list != all-zero mass list (max |diff| {worst:.3e}); "
            "a zero-mass leg is being routed to the mass-aware overload"
        )

    # 2. the bachelor's mass has to do something.
    if np.array_equal(kaon_only, allzero):
        failures.append(
            "[m_K, 0, 0] is identical to all-zero; the mass is not reaching leg 0"
        )

    # 3. and it must not reach the muon legs.
    if np.array_equal(kaon_only, broadcast):
        failures.append(
            "[m_K, 0, 0] is identical to [m_K, m_K, m_K]; the per-leg list is "
            "being broadcast and the muons have acquired a kaon energy-loss term"
        )

    scale = np.max(np.abs(allzero - 1.0))
    print(f"nominal variation size, massless : {scale:.6e}")
    print(
        "  kaon mass on leg 0 is (E/p - 1) on the e term; bachelor p = 2.006 GeV "
        "here, so E/p = 1.0299 and the expected size is 3.0%"
    )
    print(
        f"max |kaon_only - allzero|        : {np.max(np.abs(kaon_only - allzero)):.6e}"
    )
    print(
        f"max |broadcast - kaon_only|      : {np.max(np.abs(broadcast - kaon_only)):.6e}"
    )

    if failures:
        for f in failures:
            print(f"FAIL: {f}")
        raise SystemExit(1)
    print("PASS: per-leg mass rule behaves as documented")


if __name__ == "__main__":
    main()
