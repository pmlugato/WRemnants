#!/usr/bin/env python3
"""Study: apply the REAL B->J/psi K selection to the AlCaReco NanoAOD.

Standalone study -- does NOT touch the histmaker or the production chain. It
mirrors the two authoritative selections, both applied to the FITTED candidate
(JpsiXKinematicFitProducer output), not the raw four-vector sum:

  * Bmm5 nano production (../Bmm5, DileptonPlusX_cff.py) -- the candidate-level
    baseline applied at nano-build time: muon/kaon pT>1, |eta|<2.4, B(->Kll)
    mass in (4,6), two-track DOCA<0.1, kinematic vertex fit -> vertex prob.
  * btojpsik histmaker analysis path (get_bkmm_selections) -- reads the fitted
    bkmm_jpsimc_* branches; the key background rejector is
    bkmm vtx_prob (>0.1 nomc / 0.3 default) plus the fitted mass window
    5.3 +/- 0.1, muon |eta|<1.4 & pT>4, dimuon pT>7, kaon |eta|<1.4 & pT<8.

Handles on our NanoAOD: fitted mass/vtxprob = cvhFitMass / cvhFitVtxProb; muon
legs via the mu0/mu1TrackIdx cross-links into the Track table; kaon via the
candidate columns.

Cuts in the real selection that are NOT in the NanoAOD yet (listed, not applied,
so the yields are upper bounds): muon softMVA, dimuon vtx-prob / alphaBS / sl3d,
the bmm BDT, and the two-track DOCA. The point here is to check that the fitted
mass + fit-quality + kinematics sculpt the distribution correctly; these will be
loosened or dropped anyway.

Usage:  study_alcareco_nano_selection.py <nano.root> <outdir> [vtxprob]
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import uproot

MU = 0.1056583745
ORANGE, BLUE, INK, MUTED, RULE = "#E87722", "#5A7BB8", "#111111", "#6B6B6B", "#E5E5E5"


def p4(pt, eta, phi, m):
    px, py, pz = pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta)
    e = np.sqrt(px * px + py * py + pz * pz + m * m)
    return np.stack([px, py, pz, e], axis=-1)


def dimuon_pt(mu0, mu1):
    s = mu0 + mu1
    return np.hypot(s[..., 0], s[..., 1])


def style(ax):
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(RULE)
    ax.tick_params(colors=MUTED, labelsize=9)
    ax.grid(True, color=RULE, lw=0.8)
    ax.set_axisbelow(True)


def main():
    if len(sys.argv) not in (3, 4):
        sys.exit(__doc__)
    infile, outdir = sys.argv[1], sys.argv[2]
    vtxprob_cut = float(sys.argv[3]) if len(sys.argv) == 4 else 0.1
    os.makedirs(outdir, exist_ok=True)
    t = uproot.open(infile)["Events"]

    br = ["BuJpsiK_kaonPt", "BuJpsiK_kaonEta", "BuJpsiK_kaonPhi",
          "BuJpsiK_mu0TrackIdx", "BuJpsiK_mu1TrackIdx",
          "BuJpsiK_cvhFitMass", "BuJpsiK_cvhFitOk", "BuJpsiK_cvhFitVtxProb",
          "BuJpsiK_dimuonVtxProb", "BuJpsiK_dimuonAlphaBS", "BuJpsiK_dimuonSxy",
          "Track_pt", "Track_eta", "Track_phi", "Track_charge"]
    have_mm = "BuJpsiK_dimuonVtxProb" in t.keys()
    if not have_mm:
        br = [b for b in br if not b.startswith("BuJpsiK_dimuon")]
    a = t.arrays(br, library="np")

    keys = ["fitM", "vtxP", "fitOk", "mu0Pt", "mu1Pt", "mu0Eta", "mu1Eta",
            "mmPt", "kPt", "kEta", "chMul", "mmVtxP", "mmAlpha", "mmSxy"]
    acc = {k: [] for k in keys}
    for ev in range(len(a["BuJpsiK_kaonPt"])):
        tp, te, tf, tc = (a[f"Track_{x}"][ev]
                          for x in ("pt", "eta", "phi", "charge"))
        for i in range(len(a["BuJpsiK_kaonPt"][ev])):
            i0 = a["BuJpsiK_mu0TrackIdx"][ev][i]
            i1 = a["BuJpsiK_mu1TrackIdx"][ev][i]
            if not (0 <= i0 < len(tp) and 0 <= i1 < len(tp)):
                continue
            mu0 = p4(tp[i0], te[i0], tf[i0], MU)
            mu1 = p4(tp[i1], te[i1], tf[i1], MU)
            acc["fitM"].append(a["BuJpsiK_cvhFitMass"][ev][i])
            acc["vtxP"].append(a["BuJpsiK_cvhFitVtxProb"][ev][i])
            acc["fitOk"].append(a["BuJpsiK_cvhFitOk"][ev][i])
            acc["mu0Pt"].append(tp[i0])
            acc["mu1Pt"].append(tp[i1])
            acc["mu0Eta"].append(abs(te[i0]))
            acc["mu1Eta"].append(abs(te[i1]))
            acc["mmPt"].append(dimuon_pt(mu0, mu1))
            acc["kPt"].append(a["BuJpsiK_kaonPt"][ev][i])
            acc["kEta"].append(abs(a["BuJpsiK_kaonEta"][ev][i]))
            acc["chMul"].append(tc[i0] * tc[i1])
            acc["mmVtxP"].append(a["BuJpsiK_dimuonVtxProb"][ev][i] if have_mm else 1.0)
            acc["mmAlpha"].append(a["BuJpsiK_dimuonAlphaBS"][ev][i] if have_mm else 0.0)
            acc["mmSxy"].append(a["BuJpsiK_dimuonSxy"][ev][i] if have_mm else 99.0)
    d = {k: np.array(v) for k, v in acc.items()}

    n0 = len(d["fitM"])
    # The applicable real selection, on the FITTED candidate. Values from
    # get_bkmm_selections (analysis) + DileptonPlusX_cff (production kaon pT>1).
    cuts = [
        ("all fitted candidates", d["fitOk"] == 1),
        ("opposite-sign dimuon", d["chMul"] < 0),
        ("muon |eta| < 1.4", (d["mu0Eta"] < 1.4) & (d["mu1Eta"] < 1.4)),
        ("muon pT > 4", (d["mu0Pt"] > 4) & (d["mu1Pt"] > 4)),
        ("dimuon pT > 7", d["mmPt"] > 7),
        ("dimuon vtx prob > 0.1", d["mmVtxP"] > 0.1),
        ("dimuon alphaBS < 0.4", d["mmAlpha"] < 0.4),
        ("dimuon Sxy > 4 (~sl3d)", d["mmSxy"] > 4),
        ("kaon pT in (1,8)", (d["kPt"] > 1) & (d["kPt"] < 8)),
        ("kaon |eta| < 1.4", d["kEta"] < 1.4),
        (f"fitted vtx prob > {vtxprob_cut}", d["vtxP"] > vtxprob_cut),
        ("fitted mass in (5.2,5.4)", (d["fitM"] > 5.2) & (d["fitM"] < 5.4)),
    ]
    mask = np.ones(n0, bool)
    print(f"{'cut':34s} {'pass':>8s} {'cum.eff':>8s}")
    stages = {}
    for name, c in cuts:
        mask = mask & c
        print(f"{name:34s} {mask.sum():8d} {mask.sum()/max(n0,1):8.3f}")
        stages[name] = mask.copy()
    mm_note = ("dimuon Sxy is the 2D Lxy-significance proxy for the analysis's "
               "3D sl3d (true 3D needs the PV). " if have_mm else
               "dimuon vtx-prob / alphaBS / Sxy columns ABSENT -- not applied. ")
    print("NOTE: " + mm_note + "Still absent from the NanoAOD (not applied): "
          "muon softMVA>0.45, bmm BDT>0.1, two-track DOCA<0.1.")

    kin = stages["kaon |eta| < 1.4"]              # kinematics only
    full = stages["fitted mass in (5.2,5.4)"]     # full selection incl. mass band
    # mass band excluded from the plot so a peak, if any, is visible
    full_nomass = stages[f"fitted vtx prob > {vtxprob_cut}"]

    fig, ax = plt.subplots(figsize=(6.0, 3.8))
    style(ax)
    b = np.linspace(5.0, 5.6, 60)
    ax.hist(d["fitM"][d["fitOk"] == 1], bins=b, histtype="stepfilled",
            color=RULE, label=f"all fitted ({int((d['fitOk']==1).sum())})")
    klabel = "+ kin & dimuon quality" if have_mm else "+ kinematics"
    ax.hist(d["fitM"][kin], bins=b, histtype="step", color=BLUE, lw=2,
            label=f"{klabel} ({int(kin.sum())})")
    ax.hist(d["fitM"][full_nomass], bins=b, histtype="step", color=ORANGE, lw=2,
            label=f"+ B vtxProb>{vtxprob_cut} ({int(full_nomass.sum())})")
    ax.axvline(5.27934, color=INK, ls="--", lw=1.2)
    ax.axvspan(5.2, 5.4, color=ORANGE, alpha=0.06)
    ax.set_xlabel(r"fitted $m(\mu\mu K)$ [GeV]")
    ax.set_ylabel("candidates")
    ax.set_title("Real B->J/psi K selection on the fitted NanoAOD candidate",
                 color=INK, fontsize=11, loc="left", pad=8)
    ax.legend(frameon=False, fontsize=9)
    fig.savefig(os.path.join(outdir, "selection_fitB.png"), dpi=200,
                bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"final yield {int(full.sum())}/{n0}; "
          f"fitted-mass median after kin+vtxProb: "
          f"{np.median(d['fitM'][full_nomass]):.4f}" if full_nomass.sum() else "")
    print("wrote selection_fitB.png")


if __name__ == "__main__":
    main()
