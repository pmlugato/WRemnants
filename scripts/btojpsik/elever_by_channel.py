#!/usr/bin/env python3
"""The e-term lever arm, per J/psi+X channel, at generator level.

Task 0.7b of `add-second-jpsix-channel-end-to-end`. D1(a) claims the channel
choice is driven by how far E/p departs from 1, because the material term
enters the curvature bias as

    dk = e_i(eta) * (E/p) * k^2 ,      E/p = sqrt(1 + m^2/p^2)

and muons sit at E/p - 1 < 2e-4 over the whole Z/J/psi range. That claim has
been asserted with a table of E/p at fixed momenta; what it needs is the same
number folded over each channel's ACTUAL hadron momentum spectrum.

Three channels, same B production spectrum, same generator:

    B+  -> J/psi K+                 one hard kaon (the current situation)
    B0  -> J/psi K*0(-> K+ pi-)     a kaon and a near-null pion probe
    Bs  -> J/psi phi(-> K+ K-)      two kaons at ~half the momentum

The B+ arm is validated against the measured truth-signal kaon pT spectrum in
the produced histmaker output, so the generator is not trusted blind.

Usage:
  elever_by_channel.py -n 1000000
"""

import argparse
import os

import matplotlib as mpl
import numpy as np

mpl.use("Agg")
import matplotlib.pyplot as plt

M_B0, M_BU, M_BS = 5.27972, 5.27941, 5.36692
M_JPSI = 3.096900
M_KSTAR, W_KSTAR = 0.89555, 0.0473
M_PHI, W_PHI = 1.019461, 0.004249
M_K, M_PI = 0.493677, 0.13957039

HIST = "/ceph/submit/data/user/p/pmlugato/mz/histograms/btojpsik_cvhnano_softk_noclamp.hdf5"
MC_DS = "InclusiveBToJpsiX_2016PostVFP"
GEN_SIGNAL = 0

ap = argparse.ArgumentParser()
ap.add_argument("-n", "--nevents", type=int, default=1_000_000)
ap.add_argument("--seed", type=int, default=20260930)
ap.add_argument("--hist", default=HIST)
ap.add_argument(
    "-o", "--outdir", default="/home/submit/pmlugato/public_html/jpsix_second_channel"
)
args = ap.parse_args()
os.makedirs(args.outdir, exist_ok=True)
rng = np.random.default_rng(args.seed)


def p4(E, px, py, pz):
    return np.stack([E, px, py, pz], axis=-1)


def boost(v, beta):
    b2 = np.sum(beta**2, axis=-1)
    gamma = 1.0 / np.sqrt(1.0 - b2)
    bp = np.sum(beta * v[..., 1:], axis=-1)
    g2 = np.where(b2 > 1e-12, (gamma - 1.0) / np.where(b2 > 1e-12, b2, 1.0), 0.0)
    E = gamma * (v[..., 0] + bp)
    p = v[..., 1:] + (g2 * bp)[..., None] * beta + (gamma * v[..., 0])[..., None] * beta
    return p4(E, p[..., 0], p[..., 1], p[..., 2])


def two_body_p(M, m1, m2):
    return np.sqrt(
        np.clip((M**2 - (m1 + m2) ** 2) * (M**2 - (m1 - m2) ** 2), 0, None)
    ) / (2 * M)


def rel_bw(n, m0, g, lo, hi):
    out = np.empty(0)
    while out.size < n:
        x = rng.uniform(lo, hi, size=2 * (n - out.size) + 1000)
        f = x**2 * g**2 / ((x**2 - m0**2) ** 2 + m0**2 * g**2)
        out = np.concatenate([out, x[rng.random(x.size) < f / f.max()]])
    return out[:n]


def isotropic(n):
    c = rng.uniform(-1, 1, n)
    s = np.sqrt(1 - c**2)
    ph = rng.uniform(-np.pi, np.pi, n)
    return np.stack([s * np.cos(ph), s * np.sin(ph), c], axis=-1)


def make_B(n, mB):
    """Power-law pt, trigger-like threshold (see kstar_swap_toy.py: the data
    and parametric spectra give the same conclusions)."""
    pt = np.empty(0)
    while pt.size < n:
        x = rng.uniform(5.0, 100.0, size=3 * (n - pt.size) + 1000)
        f = x / (1.0 + (x / 5.5) ** 2) ** 3.2
        pt = np.concatenate([pt, x[rng.random(x.size) < f / f.max()]])
    pt = pt[:n]
    eta = rng.uniform(-2.2, 2.2, n)
    phi = rng.uniform(-np.pi, np.pi, n)
    p = pt * np.cosh(eta)
    return p4(
        np.sqrt(p**2 + mB**2), pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta)
    )


def decay_to_X(B, mB, mX, n):
    """B -> J/psi X in the B rest frame, isotropic; return X in the lab."""
    q = np.broadcast_to(np.asarray(two_body_p(mB, M_JPSI, mX), float), (n,))
    mX = np.broadcast_to(np.asarray(mX, float), (n,))
    nh = isotropic(n)
    X_rest = p4(np.sqrt(q**2 + mX**2), *(q[:, None] * nh).T)
    return boost(X_rest, B[:, 1:] / B[:, 0:1])


def split_X(X, mX, m1, m2, n):
    """X -> h1 h2 isotropically in the X frame; return both in the lab."""
    qq = two_body_p(mX, m1, m2)
    d = isotropic(n)
    h1 = p4(np.sqrt(qq**2 + m1**2), *(qq[:, None] * d).T)
    h2 = p4(np.sqrt(qq**2 + m2**2), *(-qq[:, None] * d).T)
    b = X[:, 1:] / X[:, 0:1]
    return boost(h1, b), boost(h2, b)


def pmag(v):
    return np.linalg.norm(v[..., 1:], axis=-1)


def eoverp(v, m):
    return np.sqrt(1.0 + m**2 / pmag(v) ** 2)


N = args.nevents
print(f"[elever] {N:,} decays per channel")

# --- B+ -> J/psi K+ : the kaon IS the X
Bu = make_B(N, M_BU)
Ku = decay_to_X(Bu, M_BU, M_K, N)

# --- B0 -> J/psi K*0(-> K pi)
B0 = make_B(N, M_B0)
mKst = rel_bw(N, M_KSTAR, W_KSTAR, M_K + M_PI + 1e-3, 1.35)
Kst = decay_to_X(B0, M_B0, mKst, N)
K0, PI0 = split_X(Kst, mKst, M_K, M_PI, N)

# --- Bs -> J/psi phi(-> K K)
Bs = make_B(N, M_BS)
mPhi = rel_bw(N, M_PHI, W_PHI, 2 * M_K + 1e-5, 1.10)
Phi = decay_to_X(Bs, M_BS, mPhi, N)
Ks1, Ks2 = split_X(Phi, mPhi, M_K, M_K, N)

legs = [
    (r"B$^+\to$J/$\psi$K$^+$  (K)", Ku, M_K, "C3", "-"),
    (r"B$^0\to$J/$\psi$K$^{*0}$  (K)", K0, M_K, "C0", "-"),
    (r"B$^0\to$J/$\psi$K$^{*0}$  ($\pi$)", PI0, M_PI, "C0", ":"),
    (r"B$_s\to$J/$\psi\phi$  (K, both)", np.concatenate([Ks1, Ks2]), M_K, "C2", "-"),
]

print(
    f"\n{'leg':38s} {'<p> GeV':>9} {'median p':>9} "
    f"{'<E/p-1>':>10} {'median':>9} {'frac p<2GeV':>12}"
)
print("-" * 92)
summary = {}
for name, v, m, _, _ in legs:
    p = pmag(v)
    ep = np.sqrt(1.0 + m**2 / p**2) - 1.0
    summary[name] = (np.mean(ep), np.mean(p))
    print(
        f"{name:38s} {np.mean(p):9.2f} {np.median(p):9.2f} "
        f"{np.mean(ep):10.5f} {np.median(ep):9.5f} {np.mean(p < 2.0):12.3f}"
    )

base = summary[legs[0][0]][0]
print(f"\ne-lever relative to the B+ kaon  (<E/p-1> ratio):")
for name, _, _, _, _ in legs:
    print(f"  {name:38s} {summary[name][0]/base:6.2f}x")
# per-EVENT lever: Bs has two kaons
print(f"\nper-EVENT e-lever (sum over hadron legs), relative to B+:")
tot_b0 = summary[legs[1][0]][0] + summary[legs[2][0]][0]
tot_bs = 2 * summary[legs[3][0]][0]
print(f"  B0 (K + pi) : {tot_b0/base:6.2f}x")
print(f"  Bs (K + K)  : {tot_bs/base:6.2f}x")

# ---------------------------------------------- validate B+ against measured
meas = None
try:
    import h5py

    import narf.ioutils as io

    with h5py.File(args.hist, "r") as f:
        h = io.pickle_load_h5py(f[MC_DS])["output"]["fitbins"].get()
    sel = {"genCategory": GEN_SIGNAL}
    proj = h[sel].project("fit_bachelor_pt")
    meas = (proj.axes[0].edges, proj.values())
    print(
        f"\n[validate] measured truth-signal kaon pT loaded "
        f"({meas[1].sum():.0f} entries, "
        f"edges {meas[0][0]:.2f}-{meas[0][-1]:.2f})"
    )
except Exception as e:
    print(f"\n[validate] could not load measured spectrum: {e}")

# ------------------------------------------------------------------ figures
fig, ax = plt.subplots(1, 3, figsize=(15.0, 4.6))

a = ax[0]
bins = np.linspace(0, 12, 120)
for name, v, m, c, ls in legs:
    a.hist(
        pmag(v),
        bins=bins,
        histtype="step",
        lw=1.8,
        color=c,
        ls=ls,
        density=True,
        label=name,
    )
a.set_xlabel("hadron momentum $p$ [GeV]")
a.set_ylabel("normalised")
a.set_title("Hadron momentum per channel")
a.legend(fontsize=8)

a = ax[1]
bins = np.logspace(-4, 0.3, 120)
for name, v, m, c, ls in legs:
    p = pmag(v)
    a.hist(
        np.sqrt(1 + m**2 / p**2) - 1,
        bins=bins,
        histtype="step",
        lw=1.8,
        color=c,
        ls=ls,
        density=True,
        label=name,
    )
a.axvline(2e-4, color="k", ls="--", lw=1.2)
a.text(
    2.3e-4, a.get_ylim()[1] * 0.5, "muon ceiling\n($E/p-1<2\\times10^{-4}$)", fontsize=8
)
a.set_xscale("log")
a.set_xlabel("$E/p - 1$")
a.set_ylabel("normalised")
a.set_title("The e-term lever: $E/p-1$ per leg")
a.legend(fontsize=8, loc="upper left")

a = ax[2]
if meas is not None:
    e, v = meas
    ctr = 0.5 * (e[1:] + e[:-1])
    norm = np.sum(v * np.diff(e))
    a.step(
        ctr,
        v / norm,
        where="mid",
        color="k",
        lw=1.8,
        label="measured truth signal (MC)",
    )
pt_u = np.linalg.norm(Ku[:, 1:3], axis=-1)
lo, hi = (meas[0][0], meas[0][-1]) if meas is not None else (0.0, 10.0)
sel_t = (pt_u > lo) & (pt_u < hi)
a.hist(
    pt_u[sel_t],
    bins=np.linspace(lo, hi, 60),
    histtype="step",
    lw=1.8,
    color="C3",
    density=True,
    label="toy B$^+$ kaon",
)
if meas is not None:
    e, v = meas
    c = 0.5 * (e[1:] + e[:-1])
    cdf = np.cumsum(v) / v.sum()
    med_meas = np.interp(0.5, cdf, c)
    med_toy = np.median(pt_u[sel_t])
    print(
        f"[validate] median kaon pT: measured {med_meas:.2f} GeV, "
        f"toy {med_toy:.2f} GeV  (ratio {med_toy/med_meas:.2f})"
    )
    a.text(
        0.40,
        0.55,
        f"median $p_T$\nmeasured {med_meas:.2f} GeV\n" f"toy {med_toy:.2f} GeV",
        transform=a.transAxes,
        fontsize=9,
    )
a.set_xlabel(r"kaon $p_T$ [GeV]")
a.set_ylabel("normalised")
a.set_title("Toy validated against measured B$^+$ signal")
a.legend(fontsize=8)

fig.tight_layout()
out = os.path.join(args.outdir, "elever_by_channel.png")
fig.savefig(out, dpi=140)
print(f"\n[figure] {out}")
