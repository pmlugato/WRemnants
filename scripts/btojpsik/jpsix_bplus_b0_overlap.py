#!/usr/bin/env python3
"""B+ and B0 -> J/psi K*0 are not independent samples: measure the overlap.

Task 1.8 / design.md D11. Drop the pion from a B0 -> J/psi K*0 and what
remains is m(J/psi K) BELOW the B+ mass -- which is exactly the
partially-reconstructed-B0 edge template the B+ fit already needed. So in a
joint fit the B0 signal and the B+ background template would be built from
the same tracks, and treating the channels as independent double-counts them.

The merged multi-channel nano makes this a two-row query: a B+ candidate and
a B0 candidate in the same event, sharing the same two muon tracks AND a
hadron track, are built from the same objects.

WHAT THIS DOES AND DOES NOT SHOW. On data it establishes that the candidates
share TRACKS, and that the sharing ones sit systematically lower in
m(J/psi K). It does NOT establish that the underlying decay is a true B0 --
that needs generator matching, which is task 1.7 on MC. The mass shift is
strong circumstantial evidence, not proof.

  jpsix_bplus_b0_overlap.py -i merged_nano.root
"""

import argparse
import os

import matplotlib as mpl
import numpy as np
import uproot

mpl.use("Agg")
import matplotlib.pyplot as plt

M_BU = 5.27941
M_B0 = 5.27972

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", default="/tmp/legs_all.root")
ap.add_argument(
    "-o",
    "--output",
    default="/home/submit/pmlugato/public_html/jpsix_second_channel/jpsix_bp_b0_overlap.png",
)
ap.add_argument("--bplus-id", type=int, default=0)
ap.add_argument("--b0-id", type=int, default=1)
args = ap.parse_args()

f = uproot.open(args.input)["Events"]
br = [
    "JpsiX_channelId",
    "JpsiX_kaonTrackIdx",
    "JpsiX_had1TrackIdx",
    "JpsiX_mu0TrackIdx",
    "JpsiX_mu1TrackIdx",
    "JpsiX_kvfCvhMass",
]
a = f.arrays(br, library="np")

shared, unshared = [], []
n_both = n_all_bp = 0
for ev in range(len(a["JpsiX_channelId"])):
    cid = np.asarray(a["JpsiX_channelId"][ev])
    if cid.size == 0:
        continue
    k0 = np.asarray(a["JpsiX_kaonTrackIdx"][ev])
    k1 = np.asarray(a["JpsiX_had1TrackIdx"][ev])
    m0 = np.asarray(a["JpsiX_mu0TrackIdx"][ev])
    m1 = np.asarray(a["JpsiX_mu1TrackIdx"][ev])
    mm = np.asarray(a["JpsiX_kvfCvhMass"][ev])
    ib = np.where(cid == args.bplus_id)[0]
    ik = np.where(cid == args.b0_id)[0]
    n_all_bp += ib.size
    if ib.size == 0 or ik.size == 0:
        continue
    n_both += 1
    for p in ib:
        mu_p = {int(m0[p]), int(m1[p])}
        same = False
        for q in ik:
            if mu_p == {int(m0[q]), int(m1[q])} and int(k0[p]) in {
                int(k0[q]),
                int(k1[q]),
            }:
                same = True
                break
        (shared if same else unshared).append(float(mm[p]))

shared, unshared = np.array(shared), np.array(unshared)
print(f"events with both channels : {n_both}")
print(f"B+ candidates, total       : {n_all_bp}")
print(
    f"  sharing dimuon + a hadron track with a B0 candidate : {shared.size} "
    f"({100 * shared.size / max(n_all_bp, 1):.1f}% of all B+)"
)
for lab, v in (("shared", shared), ("not shared", unshared)):
    if v.size:
        print(f"  m(J/psi K) {lab:11s}: n={v.size:4d} median {np.median(v):.4f}")
if shared.size and unshared.size:
    print(
        f"  median shift            : "
        f"{1e3 * (np.median(shared) - np.median(unshared)):+.0f} MeV"
    )

fig, ax = plt.subplots(figsize=(7.6, 4.4))
bins = np.linspace(4.95, 5.6, 45)
ax.hist(
    unshared,
    bins=bins,
    histtype="step",
    lw=1.9,
    color="C0",
    label=f"no shared track (n={unshared.size})",
)
ax.hist(
    shared,
    bins=bins,
    histtype="step",
    lw=1.9,
    color="C3",
    label=f"shares tracks with a B$^0$ cand (n={shared.size})",
)
ax.axvline(M_BU, color="k", ls=":", lw=1.6)
ax.text(M_BU + 0.008, ax.get_ylim()[1] * 0.78, f"m(B$^+$)\n{M_BU:.3f}", fontsize=9)
ax.set_xlabel(r"$m(\mathrm{J}/\psi\,K)$  [GeV]")
ax.set_ylabel("B$^+$ candidates / 14 MeV")
ax.legend(fontsize=9, loc="upper left")
fig.tight_layout()
os.makedirs(os.path.dirname(args.output), exist_ok=True)
fig.savefig(args.output, dpi=140, bbox_inches="tight")
print(f"[figure] {args.output}")
