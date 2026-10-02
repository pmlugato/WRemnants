#!/usr/bin/env python3
"""Gen-level truth per J/psi+X channel: purity, overlap, and K/pi sharing.

Tasks 1.7, 1.8 and the input to 5.2a, on the merged multi-channel MC nano.

THREE TRAPS THIS SCRIPT EXISTS TO AVOID.

1. `genBPdgId` does NOT identify the channel. It is the nearest last-copy
   b-hadron to the candidate, nothing more. Measured 2026-10-01: B0Kstar
   candidates match |genBPdgId| = 521 (a B+) MORE often than 511 (a B0).
   So "true signal" needs the right mother AND the right leg species, which
   is what `--require-legs` enforces here.

2. The excited states. |pdgId| 523 and 513 are B*+ and B*0, which decay to
   B pi / B gamma, so a candidate matched to one is still a B decay but not
   the one we want. Counted separately rather than folded in.

3. Per-leg gen columns used to stop at leg2, so the SECOND hadron of a
   four-body channel had no gen information at all. The K/pi momentum
   sharing below is impossible without leg3; if `leg3GenPdgId` is missing
   from the input this script says so rather than silently reporting the
   three-leg answer.

  jpsix_gen_truth.py -i mc_merged.root
"""

import argparse
import collections

import numpy as np
import uproot

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", default="/tmp/mc_all.root")
ap.add_argument("--dr", type=float, default=0.1, help="genBDR cut")
ap.add_argument(
    "--require-legs",
    action="store_true",
    default=True,
    help="demand the leg species match the channel, not just the mother",
)
args = ap.parse_args()

CH = {
    0: ("BPlus", 521, (321,)),
    1: ("B0Kstar", 511, (321, 211)),
    2: ("BsPhi", 531, (321, 321)),
}

f = uproot.open(args.input)["Events"]
have_leg3 = "JpsiX_leg3GenPdgId" in f.keys()
print(f"leg3 gen columns present: {have_leg3}")
if not have_leg3:
    print(
        "  -> the second hadron of every four-body channel has NO gen info.\n"
        "     Rerun the nano with the leg-3 gen columns wired before trusting\n"
        "     any K/pi sharing number below."
    )

cols = [
    "JpsiX_channelId",
    "JpsiX_genBPdgId",
    "JpsiX_genBDR",
    "JpsiX_kvfCvhMass",
    "JpsiX_kaonTrackIdx",
    "JpsiX_had1TrackIdx",
    "JpsiX_mu0TrackIdx",
    "JpsiX_mu1TrackIdx",
]
for L in (2, 3):
    for q in ("GenPdgId", "GenPt", "GenMotherPdgId"):
        c = f"JpsiX_leg{L}{q}"
        if c in f.keys():
            cols.append(c)
a = f.arrays(cols, library="np")


def cat(name):
    return np.concatenate([np.asarray(x) for x in a[name]]) if name in a else None


cid, pdg, dr = cat("JpsiX_channelId"), cat("JpsiX_genBPdgId"), cat("JpsiX_genBDR")

print(
    f"\n{'channel':9s} {'cands':>6} {'|genB|=sig':>11} {'B* (523/513)':>13} "
    f"{'other b':>9} {'no match':>9}"
)
print("-" * 66)
for c, (nm, sig, _) in CH.items():
    s = cid == c
    if s.sum() == 0:
        continue
    close = s & (dr >= 0) & (dr < args.dr)
    ap_ = np.abs(pdg)
    n_sig = int(np.sum(close & (ap_ == sig)))
    n_exc = int(np.sum(close & np.isin(ap_, [523, 513, 533])))
    n_oth = int(close.sum()) - n_sig - n_exc
    print(
        f"{nm:9s} {int(s.sum()):6d} {n_sig:11d} {n_exc:13d} {n_oth:9d} "
        f"{int(s.sum()) - int(close.sum()):9d}"
    )

# ---- K/pi momentum sharing, the f_L input (task 5.2a) --------------------
if have_leg3:
    p2, p3 = cat("JpsiX_leg2GenPt"), cat("JpsiX_leg3GenPt")
    g2, g3 = cat("JpsiX_leg2GenPdgId"), cat("JpsiX_leg3GenPdgId")
    print("\n=== gen-level hadron momentum sharing (the f_L systematic input) ===")
    print("  Restricted to GEN-MATCHED SIGNAL. Without that the sample is")
    print("  combinatorial and the hadron legs are random pions, which would")
    print("  make this measure the background's momentum sharing, not the")
    print("  resonance's. Reported both ways so the difference is visible.")
    ap_ = np.abs(pdg)
    for c, (nm, sig, want) in CH.items():
        if c == 0:
            continue
        base = (cid == c) & (p2 > 0) & (p3 > 0)
        # signal: right mother, and the two hadron legs are the right species
        legs_ok = np.isin(np.abs(g2), want) & np.isin(np.abs(g3), want)
        if len(want) == 2 and want[0] != want[1]:  # K pi: one of each
            legs_ok &= np.abs(g2) != np.abs(g3)
        sel_sig = base & (ap_ == sig) & (dr >= 0) & (dr < args.dr) & legs_ok
        for tag, s in (
            ("ALL (combinatorial-dominated)", base),
            ("gen-matched SIGNAL", sel_sig),
        ):
            if s.sum() == 0:
                print(f"  {nm:9s} {tag:30s} n=0")
                continue
            hard = np.maximum(p2[s], p3[s])
            soft = np.minimum(p2[s], p3[s])
            r = soft / hard
            print(
                f"  {nm:9s} {tag:30s} n={int(s.sum()):4d}  "
                f"soft/hard median {np.median(r):.3f}  "
                f"[16,84]% {np.percentile(r, 16):.3f}-{np.percentile(r, 84):.3f}"
            )
        continue
    for c, (nm, sig, _) in CH.items():
        if c == 0:
            continue
        s = (cid == c) & (p2 > 0) & (p3 > 0)
        if s.sum() == 0:
            continue
        print(
            f"  {nm:9s} leg species, ALL: "
            f"leg2 {collections.Counter(np.abs(g2[s]).tolist()).most_common(3)}"
        )
        print(
            f"            leg3 {collections.Counter(np.abs(g3[s]).tolist()).most_common(3)}"
        )
