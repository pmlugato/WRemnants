#!/usr/bin/env python3
"""Study: apply the btojpsik AlCaReco-path selections to the AlCaReco NanoAOD.

This is a STANDALONE study, not part of the production chain and it does not
touch the histmaker. It replays `get_bkmm_alcareco_selections` (preset A and
preset B, `btojpsik.py`) directly on the NanoAOD this stream writes, to see the
mass shape after selection and to inform a *looser* in-chain pre-filter that
could shrink the output and speed up histmaker iteration.

The selection uses the candidate columns plus the daughter->Track cross-links
(mu0/mu1/kaonTrackIdx), so it also exercises those links on real analysis cuts.

Cuts replayed (raw, class-1/2 only -- same as the histmaker's AlCaReco path):
  preset A : opposite-sign dimuon, m(mumu) in (2.95,3.25), m(mumuK) in (5.0,5.5)
  preset B : A + muon pT>4, kaon pT>1.5, |kaon eta|<1.8, mumu pT>3, B pT>5
  NOTE: the kaon-muon DOCA<0.03 cut (preset B) is NOT applied -- DOCA is not yet
  in the NanoAOD. Flagged in the cutflow so the B numbers are an upper bound.

Usage:  study_alcareco_nano_selection.py <nano.root> <outdir>
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import uproot

MU = 0.1056583745
K = 0.493677
ORANGE, BLUE, INK, MUTED, RULE = "#E87722", "#5A7BB8", "#111111", "#6B6B6B", "#E5E5E5"


def p4(pt, eta, phi, m):
    px, py, pz = pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta)
    e = np.sqrt(px * px + py * py + pz * pz + m * m)
    return np.stack([px, py, pz, e], axis=-1)


def inv(*vs):
    s = sum(vs)
    m2 = s[..., 3] ** 2 - (s[..., 0] ** 2 + s[..., 1] ** 2 + s[..., 2] ** 2)
    return np.sqrt(np.clip(m2, 0, None)), np.hypot(s[..., 0], s[..., 1])


def style(ax):
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(RULE)
    ax.tick_params(colors=MUTED, labelsize=9)
    ax.grid(True, color=RULE, lw=0.8)
    ax.set_axisbelow(True)


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    infile, outdir = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)
    t = uproot.open(infile)["Events"]

    br = ["BuJpsiK_mass", "BuJpsiK_kaonPt", "BuJpsiK_kaonEta", "BuJpsiK_kaonPhi",
          "BuJpsiK_mu0TrackIdx", "BuJpsiK_mu1TrackIdx", "BuJpsiK_cvhFitMass",
          "BuJpsiK_cvhFitOk", "Track_pt", "Track_eta", "Track_phi", "Track_charge"]
    a = t.arrays(br, library="np")

    # flatten per candidate, resolving the muon legs through the Track table
    cols = ["rawB", "mmM", "mmPt", "bPt", "muPtMin", "kPt", "kEta", "chMul",
            "fitM", "fitOk"]
    acc = {c: [] for c in cols}
    for ev in range(len(a["BuJpsiK_mass"])):
        tp, te, tf, tc = (a[f"Track_{x}"][ev]
                          for x in ("pt", "eta", "phi", "charge"))
        for i in range(len(a["BuJpsiK_mass"][ev])):
            i0 = a["BuJpsiK_mu0TrackIdx"][ev][i]
            i1 = a["BuJpsiK_mu1TrackIdx"][ev][i]
            if not (0 <= i0 < len(tp) and 0 <= i1 < len(tp)):
                continue
            mu0 = p4(tp[i0], te[i0], tf[i0], MU)
            mu1 = p4(tp[i1], te[i1], tf[i1], MU)
            kpt_i = a["BuJpsiK_kaonPt"][ev][i]
            keta_i = a["BuJpsiK_kaonEta"][ev][i]
            kp = p4(kpt_i, keta_i, a["BuJpsiK_kaonPhi"][ev][i], K)
            mm, mmpt = inv(mu0, mu1)
            _, bpt = inv(mu0, mu1, kp)
            acc["rawB"].append(a["BuJpsiK_mass"][ev][i])
            acc["mmM"].append(mm)
            acc["mmPt"].append(mmpt)
            acc["bPt"].append(bpt)
            acc["muPtMin"].append(min(tp[i0], tp[i1]))
            acc["kPt"].append(kpt_i)
            acc["kEta"].append(abs(keta_i))
            acc["chMul"].append(tc[i0] * tc[i1])
            acc["fitM"].append(a["BuJpsiK_cvhFitMass"][ev][i])
            acc["fitOk"].append(a["BuJpsiK_cvhFitOk"][ev][i])
    rawB, mmM, mmPt, bPt, muPtMin, kPt, kEta, chMul, fitM, fitOk = (
        np.array(acc[c]) for c in cols)

    n0 = len(rawB)
    cuts = [
        ("all candidates", np.ones(n0, bool)),
        ("opposite-sign dimuon", chMul < 0),
        ("m(mumu) in (2.95,3.25)", (mmM > 2.95) & (mmM < 3.25)),
        ("m(mumuK) in (5.0,5.5)  [preset A]", (rawB > 5.0) & (rawB < 5.5)),
        ("muon pT > 4", muPtMin > 4.0),
        ("kaon pT > 1.5", kPt > 1.5),
        ("|kaon eta| < 1.8", kEta < 1.8),
        ("mumu pT > 3", mmPt > 3.0),
        ("B pT > 5  [preset B, no DOCA]", bPt > 5.0),
    ]
    mask = np.ones(n0, bool)
    print(f"{'cut':40s} {'pass':>8s} {'cum.eff':>8s}")
    presetA = presetB = None
    for name, c in cuts:
        mask = mask & c
        print(f"{name:40s} {mask.sum():8d} {mask.sum()/n0:8.3f}")
        if "preset A" in name:
            presetA = mask.copy()
    presetB = mask.copy()
    print("NOTE: kaon-mu DOCA<0.03 (preset B) NOT applied -- absent from NanoAOD; "
          "preset-B yield is an upper bound.")

    # shapes: raw B mass, no selection vs preset A vs preset B
    fig, ax = plt.subplots(figsize=(5.8, 3.7))
    style(ax)
    b = np.linspace(5.0, 5.5, 60)
    ax.hist(rawB, bins=b, histtype="stepfilled", color=RULE, label=f"all ({n0})")
    ax.hist(rawB[presetA], bins=b, histtype="step", color=BLUE, lw=2,
            label=f"preset A ({presetA.sum()})")
    ax.hist(rawB[presetB], bins=b, histtype="step", color=ORANGE, lw=2,
            label=f"preset B* ({presetB.sum()})")
    ax.axvline(5.27934, color=INK, ls="--", lw=1.2)
    ax.set_xlabel(r"raw $m(\mu\mu K)$ [GeV]"); ax.set_ylabel("candidates")
    ax.set_title("AlCaReco-path selection replayed on the NanoAOD",
                 color=INK, fontsize=11, loc="left", pad=8)
    ax.legend(frameon=False, fontsize=9)
    fig.savefig(os.path.join(outdir, "selection_rawB.png"), dpi=200,
                bbox_inches="tight", facecolor="white")
    plt.close(fig)

    # fitted mass after preset B
    good = presetB & (fitOk == 1)
    if good.sum() > 5:
        fig, ax = plt.subplots(figsize=(5.8, 3.7))
        style(ax)
        ax.hist(fitM[good], bins=np.linspace(5.0, 5.6, 60), histtype="stepfilled",
                color=ORANGE, alpha=0.7)
        ax.axvline(5.27934, color=INK, ls="--", lw=1.2)
        ax.set_xlabel(r"fitted $m(\mu\mu K)$ [GeV]"); ax.set_ylabel("candidates")
        ax.set_title(f"Fitted B mass after preset B* ({good.sum()} cand.)",
                     color=INK, fontsize=11, loc="left", pad=8)
        ax.annotate(f"median {np.median(fitM[good]):.4f}", xy=(0.97, 0.94),
                    xycoords="axes fraction", ha="right", va="top",
                    color=MUTED, fontsize=9)
        fig.savefig(os.path.join(outdir, "selection_fitB.png"), dpi=200,
                    bbox_inches="tight", facecolor="white")
        plt.close(fig)
    print(f"wrote selection_rawB.png"
          + (" + selection_fitB.png" if good.sum() > 5 else ""))


if __name__ == "__main__":
    main()
