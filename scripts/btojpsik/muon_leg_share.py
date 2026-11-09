"""How much of the A/e/M mass variation do the MUON legs carry?

Done directly on the joint-CVH FITTED leg momenta, which is the point. Every
previous statement about the muon share came through the reweight machinery,
and that is parameterised in the pre-fit TRACK pT/eta -- `BuJpsiK_kaonPt` is
`Track_pt[kaonTrackIdx]` to float precision, while `BuJpsiK_jointCvhLeg2Pt`
differs from it by a median +0.61% with a 5.1% RMS. So the reweight answer and
a direct kinematic answer are not the same measurement.

For each candidate the three fitted legs are shifted by the per-eta A/e/M
uncertainty and m(mu mu K) is recomputed, three ways: bachelor only, the two
muons only, all three. The shift convention is the one the C++ uses,

    kUnc = (A - e*k)*k + q*M ,   k = 1/pt   ->   dpt/pt = -kUnc/k

and the uncertainties are the per-eta sigmas from the correction file.

TWO variants, and the difference between them is the physics:

  unconstrained  -- recompute the invariant mass from the shifted legs.
  J/psi-constrained -- after shifting, rescale the two muon three-momenta by a
      common factor so m(mu mu) returns to the PDG value, then recompute. The
      joint fit constrains m(mu mu), so to first order this is what it does. A
      COMMON muon scale shift is therefore absorbed entirely, and only the part
      that differs between the two muons survives. That is exactly the
      suppression the reweight cannot express, because a reweight is not a
      recomputation.

Example:
    python scripts/btojpsik/muon_leg_share.py -o out/ --maxFiles 300
"""

import argparse
import glob
import os
import sys

import awkward as ak
import matplotlib as mpl
import numpy as np
import uproot

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

MC = ("/ceph/submit/data/group/cms/store/mc/"
      "inclusive_btojpsix_2016postvfp_v3_NanoAOD_v2")
CALIB = ("/work/submit/pmlugato/REALmz/mz/real/WRemnants/wremnants-data/data/"
         "calibration/calibrationJDATA_MCstat_inclusive_binsel.root")

M_MU = 0.1056583745
M_K = 0.493677
M_JPSI = 3.0969

SURFACE = "#faf9f7"
INK = "#1F1F1F"
INK_MUTED = "#595959"
COLORS = {"bachelor": "#2a78d6", "muons": "#eb6834", "all": "#4a3aa7"}

BRANCHES = [
    "BuJpsiK_jointCvhLeg0Pt", "BuJpsiK_jointCvhLeg0Eta", "BuJpsiK_jointCvhLeg0Phi",
    "BuJpsiK_jointCvhLeg1Pt", "BuJpsiK_jointCvhLeg1Eta", "BuJpsiK_jointCvhLeg1Phi",
    "BuJpsiK_jointCvhLeg2Pt", "BuJpsiK_jointCvhLeg2Eta", "BuJpsiK_jointCvhLeg2Phi",
    "BuJpsiK_mu0TrackIdx", "BuJpsiK_mu1TrackIdx", "BuJpsiK_kaonTrackIdx",
    "Track_charge", "Track_eta", "Track_pt",
    "BuJpsiK_jointCvhOk", "BuJpsiK_jointCvhMass", "BuJpsiK_jointCvhMassErr",
    "BuJpsiK_jointCvhEdmRef", "BuJpsiK_jointCvhChisq", "BuJpsiK_jointCvhNdof",
    "BuJpsiK_nLegsRefit", "BuJpsiK_kaonPt", "BuJpsiK_kaonEta",
    "BuJpsiK_genCategory" if False else "BuJpsiK_leg2GenPt",
]


def load_sigmas():
    f = uproot.open(CALIB)
    edges = f["A"].axis().edges()
    return edges, {k: f[k].errors() for k in ("A", "e", "M")}


def eta_bin(eta, edges):
    return np.clip(np.digitize(eta, edges) - 1, 0, len(edges) - 2)


def read(nfiles):
    files = sorted(glob.glob(os.path.join(MC, "*.root")))[:nfiles]
    parts = []
    for fn in files:
        try:
            t = uproot.open(fn)["Events"]
            parts.append(t.arrays(BRANCHES))
        except Exception:
            continue
    return ak.concatenate(parts) if parts else None


def flatten_selected(a):
    """Apply the histmaker's active selection and return flat per-leg arrays."""
    # Cross-link BEFORE masking so the Track indices still line up.
    q = [a["Track_charge"][a[f"BuJpsiK_{n}TrackIdx"]] for n in ("mu0", "mu1", "kaon")]
    trk_eta = [a["Track_eta"][a[f"BuJpsiK_{n}TrackIdx"]] for n in ("mu0", "mu1", "kaon")]

    sent = -98.0
    ok = (
        (a["BuJpsiK_jointCvhOk"] == 1)
        & (a["BuJpsiK_jointCvhLeg0Pt"] > sent)
        & (a["BuJpsiK_jointCvhLeg1Pt"] > sent)
        & (a["BuJpsiK_jointCvhLeg2Pt"] > sent)
        & (a["BuJpsiK_jointCvhMass"] > sent)
        & (a["BuJpsiK_jointCvhEdmRef"] >= 0) & (a["BuJpsiK_jointCvhEdmRef"] < 1e-2)
        & (a["BuJpsiK_jointCvhChisq"] > 0) & (a["BuJpsiK_jointCvhChisq"] < 1e6)
        & (a["BuJpsiK_jointCvhNdof"] > 0) & (a["BuJpsiK_jointCvhNdof"] < 200)
        & (a["BuJpsiK_jointCvhMassErr"] > 0) & (a["BuJpsiK_jointCvhMassErr"] < 1.0)
        & (q[0] * q[1] < 0)
        & (a["BuJpsiK_kaonPt"] > 1.0) & (abs(a["BuJpsiK_kaonEta"]) < 2.4)
        & (a["BuJpsiK_jointCvhMass"] > 5.0) & (a["BuJpsiK_jointCvhMass"] < 5.5)
        & (a["BuJpsiK_nLegsRefit"] >= 3)
    )
    # muon pT > 3 on the TRACK pT, as the histmaker does
    mupt = [a["Track_pt"][a[f"BuJpsiK_{n}TrackIdx"]] for n in ("mu0", "mu1")]
    ok = ok & (mupt[0] > 3.0) & (mupt[1] > 3.0)
    ok = ok & (abs(trk_eta[0]) < 2.4) & (abs(trk_eta[1]) < 2.4)

    f = lambda x: ak.flatten(x[ok]).to_numpy()
    legs = []
    for i in range(3):
        legs.append(dict(
            pt=f(a[f"BuJpsiK_jointCvhLeg{i}Pt"]),
            eta=f(a[f"BuJpsiK_jointCvhLeg{i}Eta"]),
            phi=f(a[f"BuJpsiK_jointCvhLeg{i}Phi"]),
            q=f(q[i]).astype(float),
            trk_eta=f(trk_eta[i]),
            mass=M_MU if i < 2 else M_K,
        ))
    return legs, f(a["BuJpsiK_jointCvhMass"])


def p4(pt, eta, phi, m):
    px, py, pz = pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta)
    return px, py, pz, np.sqrt(pt * pt * np.cosh(eta) ** 2 + m * m)


def dpt_over_pt(leg, param, sig, edges):
    """dpt/pt for a one-sigma shift of `param` in this leg's eta bin.

    kUnc = (A - e*k)*k + q*M with k = 1/pt, and dpt/pt = -kUnc/k.
    """
    s = sig[param][eta_bin(leg["trk_eta"], edges)]
    k = 1.0 / leg["pt"]
    if param == "A":
        kunc = s * k
    elif param == "e":
        kunc = -s * k * k
    else:
        kunc = leg["q"] * s
    return -kunc / k


def masses(legs, scales):
    px = py = pz = e = 0.0
    for leg, sc in zip(legs, scales):
        a, b, c, d = p4(leg["pt"] * sc, leg["eta"], leg["phi"], leg["mass"])
        px, py, pz, e = px + a, py + b, pz + c, e + d
    return np.sqrt(np.maximum(e * e - px * px - py * py - pz * pz, 0.0))


def constrain_jpsi(legs, scales):
    """Rescale both muon three-momenta by a common lambda so m(mu mu) = PDG.

    Solved by bisection on lambda: m(mu mu) is monotonic in a common momentum
    scale, so a bracketed bisection is exact to machine precision and needs no
    derivative. 40 iterations takes the bracket below 1e-12.
    """
    def mmumu(lam):
        px = py = pz = e = 0.0
        for i in (0, 1):
            a, b, c, d = p4(legs[i]["pt"] * scales[i] * lam, legs[i]["eta"],
                            legs[i]["phi"], legs[i]["mass"])
            px, py, pz, e = px + a, py + b, pz + c, e + d
        return np.sqrt(np.maximum(e * e - px * px - py * py - pz * pz, 0.0))

    lo = np.full_like(legs[0]["pt"], 0.5)
    hi = np.full_like(legs[0]["pt"], 2.0)
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        too_big = mmumu(mid) > M_JPSI
        hi = np.where(too_big, mid, hi)
        lo = np.where(too_big, lo, mid)
    lam = 0.5 * (lo + hi)
    return [scales[0] * lam, scales[1] * lam, scales[2]]


def shift(legs, param, sig, edges, which, constrained):
    """Mean-mass shift in MeV for a one-sigma variation of `param`."""
    d = [dpt_over_pt(l, param, sig, edges) for l in legs]
    use = {"bachelor": [0, 0, 1], "muons": [1, 1, 0], "all": [1, 1, 1]}[which]
    out = []
    for s in (+1.0, -1.0):
        sc = [1.0 + s * d[i] * use[i] for i in range(3)]
        if constrained:
            sc = constrain_jpsi(legs, sc)
        out.append(masses(legs, sc))
    nom = masses(legs, constrain_jpsi(legs, [1.0, 1.0, 1.0]) if constrained
                 else [1.0, 1.0, 1.0])
    return 0.5 * ((out[0] - nom) - (out[1] - nom)) * 1e3, nom


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("-o", "--outfolder", required=True)
    p.add_argument("--maxFiles", type=int, default=300)
    args = p.parse_args()

    os.makedirs(args.outfolder, exist_ok=True)
    edges, sig = load_sigmas()
    a = read(args.maxFiles)
    legs, mfit = flatten_selected(a)
    n = len(mfit)
    print(f"{n} selected candidates from {args.maxFiles} files\n")

    res = {}
    for param in ("A", "e", "M"):
        for constrained in (False, True):
            for which in ("bachelor", "muons", "all"):
                s, _ = shift(legs, param, sig, edges, which, constrained)
                res[(param, constrained, which)] = s

    print("Mean-mass shift per one-sigma variation [MeV], and the MUON SHARE")
    print("share = |muons| / (|muons| + |bachelor|), so 50% means equal.\n")
    hdr = f"{'':<22}{'bachelor':>11}{'muons':>11}{'all':>11}{'muon share':>13}"
    for constrained in (False, True):
        tag = "J/psi CONSTRAINED" if constrained else "unconstrained"
        print(f"-- {tag}")
        print(hdr)
        for param in ("A", "e", "M"):
            b = res[(param, constrained, "bachelor")].mean()
            m = res[(param, constrained, "muons")].mean()
            t = res[(param, constrained, "all")].mean()
            share = abs(m) / (abs(m) + abs(b)) if (abs(m) + abs(b)) else np.nan
            print(f"{param:<22}{b:>11.4f}{m:>11.4f}{t:>11.4f}{share:>12.1%}")
        print()

    # ---- figure: the distribution of the induced shift, per leg set
    for constrained in (False, True):
        tag = "constrained" if constrained else "unconstrained"
        fig, axes = plt.subplots(1, 3, figsize=(12.4, 3.9))
        for ax, param in zip(axes, ("A", "e", "M")):
            allv = np.concatenate([res[(param, constrained, w)]
                                   for w in ("bachelor", "muons", "all")])
            lim = np.percentile(np.abs(allv), 99.5)
            bins = np.linspace(-lim, lim, 80)
            for w in ("bachelor", "muons", "all"):
                v = res[(param, constrained, w)]
                ax.hist(v, bins=bins, histtype="step", linewidth=1.6,
                        color=COLORS[w],
                        label=f"{w}:  mean {v.mean():+.4f} MeV")
            ax.axvline(0, color=INK_MUTED, linewidth=0.9)
            ax.set_xlabel(f"induced mass shift, 1σ {param} [MeV]")
            ax.set_yscale("log")
            ax.legend(fontsize=7.5, frameon=False)
            ax.grid(alpha=0.5)
            ax.set_axisbelow(True)
        axes[0].set_ylabel("candidates")
        fig.suptitle(
            f"A/e/M variation on the joint-CVH FITTED legs — {tag}"
            + ("  (m(μμ) pinned to the PDG value, as the fit does)"
               if constrained else "  (invariant mass recomputed, no constraint)"),
            fontsize=10.5, color=INK, x=0.005, ha="left")
        out = os.path.join(args.outfolder, f"muon_leg_share_{tag}.png")
        fig.savefig(out, dpi=160, bbox_inches="tight", facecolor=SURFACE)
        plt.close(fig)
        print("wrote", out)


if __name__ == "__main__":
    main()
