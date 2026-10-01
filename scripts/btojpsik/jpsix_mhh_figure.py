#!/usr/bin/env python3
"""Two-hadron invariant mass per channel, from the merged J/psi+X nano.

Task 1.2 validation. Until 2026-10-01 the four-body nano recorded the
COMPOSITE resonance rather than its two tracks -- `kaonPdgId` read 313 for
K*0 and 333 for phi -- and the second hadron had no columns at all, because
the driver's per-leg loop stopped at 3. So `m(hh)` could not be formed, and
with it neither the K-pi swap veto nor the phi mass anchor.

This script forms `m(hh)` from the newly exposed `kaonTrackIdx` and
`had1TrackIdx` cross-links into the Track table, under each mass assignment,
and checks it against the published resonance masses. If the phi appears at
1.0195 GeV and the K*0 near 0.896 GeV, the four-body chain works end to end.

  jpsix_mhh_figure.py -i merged_nano.root
"""

import argparse
import os

import matplotlib as mpl
import numpy as np
import uproot

mpl.use("Agg")
import matplotlib.pyplot as plt

M_K, M_PI = 0.493677, 0.13957039
M_PHI, M_KSTAR = 1.019461, 0.89555
CHAN = {
    0: "B$^+\\to$J/$\\psi$K$^+$",
    1: "B$^0\\to$J/$\\psi$K$^{*0}$",
    2: "B$_s\\to$J/$\\psi\\phi$",
}

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", default="/tmp/legs_all.root")
ap.add_argument(
    "-o",
    "--output",
    default="/home/submit/pmlugato/public_html/jpsix_second_channel/jpsix_mhh.png",
)
args = ap.parse_args()

f = uproot.open(args.input)["Events"]
br = [
    "JpsiX_channelId",
    "JpsiX_kaonTrackIdx",
    "JpsiX_had1TrackIdx",
    "Track_pt",
    "Track_eta",
    "Track_phi",
]
a = f.arrays(br, library="np")


def mass_pair(pt1, eta1, phi1, m1, pt2, eta2, phi2, m2):
    def c(pt, eta, phi, m):
        px, py, pz = pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta)
        return np.sqrt(px * px + py * py + pz * pz + m * m), px, py, pz

    E1, x1, y1, z1 = c(pt1, eta1, phi1, m1)
    E2, x2, y2, z2 = c(pt2, eta2, phi2, m2)
    s = (E1 + E2) ** 2 - ((x1 + x2) ** 2 + (y1 + y2) ** 2 + (z1 + z2) ** 2)
    return np.sqrt(max(s, 0.0))


acc = {1: {"KK": [], "Kpi": [], "piK": []}, 2: {"KK": [], "Kpi": [], "piK": []}}
for ev in range(len(a["JpsiX_channelId"])):
    cid = np.asarray(a["JpsiX_channelId"][ev])
    if cid.size == 0:
        continue
    i0 = np.asarray(a["JpsiX_kaonTrackIdx"][ev])
    i1 = np.asarray(a["JpsiX_had1TrackIdx"][ev])
    tp = np.asarray(a["Track_pt"][ev])
    te = np.asarray(a["Track_eta"][ev])
    tf = np.asarray(a["Track_phi"][ev])
    for j in range(cid.size):
        c = int(cid[j])
        if c not in acc or i0[j] < 0 or i1[j] < 0:
            continue
        A, B = int(i0[j]), int(i1[j])
        if A >= len(tp) or B >= len(tp):
            continue
        acc[c]["KK"].append(
            mass_pair(tp[A], te[A], tf[A], M_K, tp[B], te[B], tf[B], M_K)
        )
        acc[c]["Kpi"].append(
            mass_pair(tp[A], te[A], tf[A], M_K, tp[B], te[B], tf[B], M_PI)
        )
        acc[c]["piK"].append(
            mass_pair(tp[A], te[A], tf[A], M_PI, tp[B], te[B], tf[B], M_K)
        )

fig, ax = plt.subplots(1, 2, figsize=(11.0, 4.3))

a0 = ax[0]
d = acc[1]
a0.hist(
    d["Kpi"],
    bins=40,
    range=(0.6, 1.3),
    histtype="step",
    lw=1.8,
    color="C0",
    label="K$^+\\pi^-$ assignment",
)
a0.hist(
    d["piK"],
    bins=40,
    range=(0.6, 1.3),
    histtype="step",
    lw=1.8,
    color="C1",
    ls="--",
    label="$\\pi^+$K$^-$ assignment",
)
a0.axvline(M_KSTAR, color="k", ls=":", lw=1.6)
a0.text(
    M_KSTAR + 0.012, a0.get_ylim()[1] * 0.80, f"K*(892)$^0$\n{M_KSTAR:.4f}", fontsize=9
)
a0.set_xlabel("$m(hh)$ [GeV]")
a0.set_ylabel(f"candidates / 17.5 MeV   (n={len(d['Kpi'])})")
a0.set_title(CHAN[1] + ", both assignments")
a0.legend(fontsize=9)

a1 = ax[1]
d = acc[2]
a1.hist(
    d["KK"],
    bins=30,
    range=(0.98, 1.10),
    histtype="step",
    lw=2.0,
    color="C2",
    label="K$^+$K$^-$ assignment",
)
a1.axvline(M_PHI, color="k", ls=":", lw=1.6)
a1.text(
    M_PHI + 0.003, a1.get_ylim()[1] * 0.72, f"$\\phi$(1020)\n{M_PHI:.4f}", fontsize=9
)
a1.set_xlabel("$m(hh)$ [GeV]")
a1.set_ylabel(f"candidates / 4 MeV   (n={len(d['KK'])})")
a1.set_title(CHAN[2] + ", no assignment ambiguity")
a1.legend(fontsize=9)

fig.tight_layout()
os.makedirs(os.path.dirname(args.output), exist_ok=True)
fig.savefig(args.output, dpi=140, bbox_inches="tight")
print(f"[figure] {args.output}")
for c in (1, 2):
    for k, v in acc[c].items():
        if v:
            print(f"  ch{c} m({k}): n={len(v)} median {np.median(v):.4f}")
