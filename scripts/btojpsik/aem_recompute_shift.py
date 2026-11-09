#!/usr/bin/env python3
"""Analytic recompute of the A/e/M induced mass shift (worklist 4.2).

Applies the A/e/M shift directly to the reconstructed four-vectors of truth
B+ -> J/psi(mumu) K+ candidates and rebuilds m(mu mu K), instead of reweighting
the template. Four configurations, whose differences are the point:

  legs = all | kaon    which legs the shift is applied to
  jpsi = free | fixed  'fixed' rescales the dimuon momenta after the shift so
                       m(mu mu) is back at the J/psi mass, which is what the
                       joint CVH refit imposes. A common muon scale error should
                       then leave m(mu mu K) almost unchanged; the reweighted
                       templates cannot know this, which is what we are testing.

The shift itself is muon_calibration.hpp's calculateShiftedPtExact, transcribed:
A and M act on 1/pt, e is a shift of the total energy.

    pt' = sqrt((sqrt((pt cosh eta)^2 + m^2) + e)^2 - m^2) / cosh eta
    k'  = (1 + A)/pt' + q M ,   pt_final = 1/k'

Reports d m(mu mu K) / d(parameter) in MeV per unit of A, e and M, per kaon pT
band, from the half up-down difference of the mean at +/- the probe shift.
Run in the suephist environment (uproot).
"""
import argparse, glob, os
import numpy as np
import uproot

MU, MK, MJPSI = 0.1056583745, 0.493677, 3.0969
C = "BuJpsiK"

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--inputDir", default="/ceph/submit/data/group/cms/store/mc/inclusive_btojpsix_2016postvfp_v3_NanoAOD_v2")
ap.add_argument("--maxFiles", type=int, default=150)
ap.add_argument("--probeA", type=float, default=1e-4)
ap.add_argument("--probeE", type=float, default=1e-3)
ap.add_argument("--probeM", type=float, default=1e-5)
ap.add_argument("--ptEdges", type=float, nargs="+", default=[1.0, 2.0, 2.8, 5.7, 10.0])
ap.add_argument("-o", "--outdir", required=True)
a = ap.parse_args(); os.makedirs(a.outdir, exist_ok=True)


def shifted_pt(pt, eta, q, dA, dE, dM, mass):
    ch = np.cosh(eta); p = pt * ch
    E = np.sqrt(p * p + mass * mass) + dE
    p2 = np.clip(E * E - mass * mass, 1e-12, None)
    pt1 = np.sqrt(p2) / ch
    k = (1.0 + dA) / pt1 + q * dM
    return 1.0 / np.clip(k, 1e-12, None)


def p4(pt, eta, phi, mass):
    return (pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta),
            np.sqrt((pt * np.cosh(eta)) ** 2 + mass * mass))


def mass_of(parts):
    px = sum(p[0] for p in parts); py = sum(p[1] for p in parts)
    pz = sum(p[2] for p in parts); E = sum(p[3] for p in parts)
    return np.sqrt(np.clip(E * E - px * px - py * py - pz * pz, 0, None))


def build(d, dA, dE, dM, legs, jpsi):
    """m(mu mu K) after shifting `legs` ('all' or 'kaon') by (dA, dE, dM)."""
    sh = dict(all=("mu0", "mu1", "k"), kaon=("k",))[legs]
    out = {}
    for tag, mass in (("mu0", MU), ("mu1", MU), ("k", MK)):
        pt = d[f"{tag}_pt"]
        if tag in sh:
            pt = shifted_pt(pt, d[f"{tag}_eta"], d[f"{tag}_q"], dA, dE, dM, mass)
        out[tag] = p4(pt, d[f"{tag}_eta"], d[f"{tag}_phi"], mass)
    if jpsi == "fixed":
        # rescale both muon momenta by a common factor so m(mu mu) = m_Jpsi,
        # the constraint the joint CVH refit imposes
        mmm = mass_of([out["mu0"], out["mu1"]])
        s = np.ones_like(mmm)
        ok = mmm > 2 * MU
        # m^2 ~ 2(E1E2 - p1.p2); solve for the scale on |p| numerically-free:
        # for a common scale s on both momenta, m^2(s) is monotonic, so bisect
        lo, hi = np.full_like(s, 0.5), np.full_like(s, 2.0)
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            m0 = [mid * out["mu0"][i] for i in range(3)]
            m1 = [mid * out["mu1"][i] for i in range(3)]
            e0 = np.sqrt(sum(x * x for x in m0) + MU ** 2)
            e1 = np.sqrt(sum(x * x for x in m1) + MU ** 2)
            mm = np.sqrt(np.clip((e0 + e1) ** 2 - sum((m0[i] + m1[i]) ** 2 for i in range(3)), 0, None))
            lo = np.where(mm < MJPSI, mid, lo); hi = np.where(mm < MJPSI, hi, mid)
        s = np.where(ok, 0.5 * (lo + hi), 1.0)
        for tag in ("mu0", "mu1"):
            v = out[tag]; sp = [s * v[i] for i in range(3)]
            out[tag] = (sp[0], sp[1], sp[2], np.sqrt(sum(x * x for x in sp) + MU ** 2))
    return mass_of([out["mu0"], out["mu1"], out["k"]])


cols = [f"{C}_mu0TrackIdx", f"{C}_mu1TrackIdx", f"{C}_kaonTrackIdx",
        f"{C}_leg0GenIdx", f"{C}_leg1GenIdx", f"{C}_leg2GenIdx", f"{C}_nLegsGenMatched",
        "Track_pt", "Track_eta", "Track_phi", "Track_charge"]
acc = {k: [] for k in ("mu0_pt", "mu0_eta", "mu0_phi", "mu0_q", "mu1_pt", "mu1_eta", "mu1_phi", "mu1_q",
                       "k_pt", "k_eta", "k_phi", "k_q")}
files = sorted(glob.glob(os.path.join(a.inputDir, "*.root")))[: a.maxFiles]
nev = 0
for f in files:
    try:
        t = uproot.open(f)["Events"]
        arr = t.arrays(cols, library="ak")
    except Exception:
        continue
    import awkward as ak
    # nLegsGenMatched is per CANDIDATE, so the mask selects candidates, not
    # events; the track indices it yields then index the full per-event Track
    # collection. Masking the Track arrays with a candidate-level mask is the
    # obvious wrong turn and awkward will not broadcast it.
    cand = (arr[f"{C}_nLegsGenMatched"] == 3) & (arr[f"{C}_mu0TrackIdx"] >= 0) \
         & (arr[f"{C}_mu1TrackIdx"] >= 0) & (arr[f"{C}_kaonTrackIdx"] >= 0)
    for tag, idxcol in (("mu0", f"{C}_mu0TrackIdx"), ("mu1", f"{C}_mu1TrackIdx"), ("k", f"{C}_kaonTrackIdx")):
        idx = arr[idxcol][cand]
        for v, src in (("pt", "Track_pt"), ("eta", "Track_eta"), ("phi", "Track_phi"), ("q", "Track_charge")):
            acc[f"{tag}_{v}"].append(ak.to_numpy(ak.flatten(arr[src][idx])).astype(float))
    nev += int(ak.sum(ak.num(arr[f"{C}_mu0TrackIdx"][cand])))
d = {k: np.concatenate(v) for k, v in acc.items()}
print(f"{len(files)} files, {nev} truth-matched candidates")

edges = np.array(a.ptEdges)
lines = ["# Analytic recompute of the induced mass shift (worklist 4.2)", "",
         f"{nev} truth-matched B+ candidates from {len(files)} files. "
         f"Probe shifts: A = {a.probeA:g}, e = {a.probeE:g} GeV, M = {a.probeM:g}. "
         "Values are d m(mumuK)/d(parameter) in MeV per unit parameter, from the half "
         "up-down difference of the mean.", "",
         "| kaon pT [GeV] | N | config | dm/dA [MeV] | dm/de [MeV/GeV] | dm/dM [MeV] |",
         "|---|--:|---|--:|--:|--:|"]
res = {}
for lo, hi in zip(edges[:-1], edges[1:]):
    m = (d["k_pt"] >= lo) & (d["k_pt"] < hi)
    sub = {k: v[m] for k, v in d.items()}
    for legs in ("all", "kaon"):
        for jpsi in ("free", "fixed"):
            row = []
            for p, probe in (("A", a.probeA), ("e", a.probeE), ("M", a.probeM)):
                dA, dE, dM = (probe if p == "A" else 0), (probe if p == "e" else 0), (probe if p == "M" else 0)
                up = build(sub, dA, dE, dM, legs, jpsi).mean()
                dn = build(sub, -dA, -dE, -dM, legs, jpsi).mean()
                row.append(1e3 * 0.5 * (up - dn) / probe)
            res[(lo, hi, legs, jpsi)] = row
            lines.append(f"| {lo:g}-{hi:g} | {m.sum()} | legs={legs}, J/psi {jpsi} | "
                         f"{row[0]:.1f} | {row[1]:.1f} | {row[2]:.1f} |")
lines += ["", "Reading: comparing 'J/psi free' with 'J/psi fixed' at legs=all shows how much of the "
          "muon response the mass constraint removes. Comparing legs=all with legs=kaon at "
          "J/psi fixed shows what the muon legs still add once the constraint is respected. "
          "A sign change of dm/de between pT bands is the signature the 2 GeV clamp would leave "
          "in the REWEIGHTED templates; the recompute here has no clamp, so a sign change in this "
          "table would be kinematic and a disagreement with the reweight would be the clamp."]
open(os.path.join(a.outdir, "recompute_summary.md"), "w").write("\n".join(lines) + "\n")
print("\n".join(lines))

import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
COL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
ctr = np.sqrt(edges[:-1] * edges[1:])
fig, axs = plt.subplots(1, 3, figsize=(14, 4.2))
for k, (pi, lab) in enumerate([(0, "dm/dA [MeV per unit A]"), (1, "dm/de [MeV per GeV of e]"), (2, "dm/dM [MeV per unit M]")]):
    for j, (legs, jpsi) in enumerate([("all", "free"), ("all", "fixed"), ("kaon", "free"), ("kaon", "fixed")]):
        y = [res[(lo, hi, legs, jpsi)][pi] for lo, hi in zip(edges[:-1], edges[1:])]
        axs[k].plot(ctr, y, "o-", color=COL[j], label=f"legs={legs}, J/psi {jpsi}", ms=6, lw=1.6)
    axs[k].axhline(0, color="#52514e", lw=0.8); axs[k].axvline(2.0, color="#52514e", ls=":", lw=1.0)
    axs[k].set_xscale("log"); axs[k].set_xlabel("kaon pT [GeV]"); axs[k].set_title(lab, fontsize=10)
    axs[k].grid(color="#e6e6e3")
axs[0].legend(frameon=False, fontsize=8.5)
fig.suptitle("Recomputed induced mass shift (dotted line: the 2 GeV ONNX clamp)", fontsize=11)
fig.tight_layout(); fig.savefig(os.path.join(a.outdir, "recompute_shift.png"), dpi=130)
print("wrote", os.path.join(a.outdir, "recompute_shift.png"))
