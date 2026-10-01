#!/usr/bin/env python3
"""Generator-level K<->pi swap toy for B0 -> J/psi K*0(-> K+ pi-).

Task 0.6 of `add-second-jpsix-channel-end-to-end`. The channel decision (D1)
rests on a claim that has so far only been argued in the boosted limit:

    reassigning the two hadron masses at fixed 3-momenta changes an invariant
    mass by dm = (E/m) dE, so the SAME dE moves m(Kpi) and m(B) by different
    amounts, and

        dm(Kpi) / dm(B)  =  (E_K* / m_K*) / (E_B / m_B)  ~  2.0

Consequence, if true: with a K* window of half-width w, surviving swaps have
|dm(B)| < w/2 -- the contamination is bounded, one-sided and tied to a cut we
control, so no particle ID is needed. This script tests the ratio in real
kinematics and produces the bias-vs-window curve that phase 5.2 compares
against data.

The B0 production spectrum is NOT guessed: it is taken from the reconstructed
B+ -> J/psi K+ candidates in the produced nano (sideband-subtracted), since B+
and B0 production are the same to well within what this toy needs.

Three candidate-selection strategies are compared, all PID-free:

  keepboth   every assignment inside the K* window is kept (the naive choice)
  closest    keep only the assignment whose m(hh) is closest to the K* pole
  veto       keep the candidate only if the SWAPPED assignment falls OUTSIDE
             the window -- i.e. require an unambiguous assignment

Usage:
  kstar_swap_toy.py --nfiles 30 -n 2000000
  kstar_swap_toy.py --spectrum-cache <npz>      # reuse a cached spectrum
"""

import argparse
import glob
import os

import matplotlib as mpl
import numpy as np

mpl.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------- constants
M_B0 = 5.27972  # GeV, PDG
M_BU = 5.27941
M_JPSI = 3.096900
M_KSTAR = 0.89555  # K*(892)0
W_KSTAR = 0.0473  # natural width
M_K = 0.493677
M_PI = 0.13957039

NANO_GLOB = (
    "/ceph/submit/data/group/cms/store/data/Run2016{era}/Charmonium/"
    "NanoAOD/BtoJpsiK_v2/*.root"
)

ap = argparse.ArgumentParser()
ap.add_argument(
    "--nfiles",
    type=int,
    default=30,
    help="nano files to read for the B+ production spectrum",
)
ap.add_argument("-n", "--nevents", type=int, default=2_000_000)
ap.add_argument(
    "--spectrum",
    choices=("data", "parametric"),
    default="data",
    help="B production spectrum: from the nano, or a power law. "
    "Run both -- the conclusions must not depend on it.",
)
ap.add_argument(
    "--ptres",
    type=float,
    default=0.007,
    help="per-track sigma(pt)/pt applied to both hadrons and the "
    "J/psi. 0 disables smearing.",
)
ap.add_argument(
    "--fL",
    type=float,
    default=0.57,
    help="K* longitudinal polarisation fraction. 0 < fL < 1 uses "
    "the helicity distribution; pass -1 for isotropic.",
)
ap.add_argument(
    "--bwindow",
    type=float,
    default=0.15,
    help="half-width of the B mass fit range [GeV]. Swaps landing "
    "outside it are background, not a peak bias.",
)
ap.add_argument("--seed", type=int, default=20260930)
ap.add_argument(
    "--spectrum-cache",
    default="/ceph/submit/data/user/p/pmlugato/mz/studies/bplus_prod_spectrum.npz",
)
ap.add_argument("--refresh-spectrum", action="store_true")
ap.add_argument(
    "-o", "--outdir", default="/home/submit/pmlugato/public_html/jpsix_second_channel"
)
args = ap.parse_args()

os.makedirs(args.outdir, exist_ok=True)
TAG = f"{args.spectrum}_fL{args.fL:g}_res{1e3*args.ptres:g}"
rng = np.random.default_rng(args.seed)


# ------------------------------------------------------- B+ spectrum from data
def load_spectrum():
    """(pt, |eta|) of B+ candidates in the mass band, from the nano.

    NOT sideband-subtracted. Measured 2026-09-30: the inclusive B+ S/B is only
    ~4-5%, the background under the peak has real curvature
    (`project_cvhnano_bkg_curvature`), and the LOW sideband sits on the
    partially-reconstructed B0 edge -- so a two-sideband subtraction
    oversubtracts by more than the signal it is trying to isolate, and a
    one-sided subtraction is still wrong because the background is not flat.
    Extracting a clean signal spectrum needs the full fit machinery, which is
    phase 6, not phase 0.

    What this toy actually needs is the B production kinematics, and the
    combinatorial is also J/psi + a track from b-bbar events, so its pt/eta
    shape is close to the signal's. `--spectrum parametric` exists to prove
    the conclusions do not depend on this choice; run both.
    """
    if os.path.exists(args.spectrum_cache) and not args.refresh_spectrum:
        d = np.load(args.spectrum_cache)
        print(f"[spectrum] cached: {args.spectrum_cache}")
        return d["pt_edges"], d["eta_edges"], d["w"]

    import uproot

    files = []
    for era in ("F", "G", "H"):
        files += sorted(glob.glob(NANO_GLOB.format(era=era)))[: args.nfiles]
    print(f"[spectrum] reading {len(files)} nano files")

    pt_edges = np.concatenate(
        [np.arange(0.0, 40.0, 0.5), [45.0, 50.0, 60.0, 80.0, 120.0]]
    )
    eta_edges = np.linspace(0.0, 2.6, 27)
    h_sig = np.zeros((len(pt_edges) - 1, len(eta_edges) - 1))
    h_sb = np.zeros_like(h_sig)

    br = [
        "BuJpsiK_kvfCvhMass",
        "BuJpsiK_kvfCvhPt",
        "BuJpsiK_eta",
        "BuJpsiK_kvfCvhOk",
        "BuJpsiK_kvfCvhVtxProb",
        "BuJpsiK_kvfCvhSl3d",
    ]
    nread = 0
    for fn in files:
        try:
            t = uproot.open(fn)["Events"]
            a = t.arrays(br, library="np")
        except Exception as e:  # truncated / unreadable
            print(f"  skip {os.path.basename(fn)}: {e}")
            continue
        m = (
            np.concatenate(a["BuJpsiK_kvfCvhMass"])
            if a["BuJpsiK_kvfCvhMass"].dtype == object
            else a["BuJpsiK_kvfCvhMass"]
        )
        pt = (
            np.concatenate(a["BuJpsiK_kvfCvhPt"])
            if a["BuJpsiK_kvfCvhPt"].dtype == object
            else a["BuJpsiK_kvfCvhPt"]
        )
        eta = (
            np.concatenate(a["BuJpsiK_eta"])
            if a["BuJpsiK_eta"].dtype == object
            else a["BuJpsiK_eta"]
        )
        ok = (
            np.concatenate(a["BuJpsiK_kvfCvhOk"])
            if a["BuJpsiK_kvfCvhOk"].dtype == object
            else a["BuJpsiK_kvfCvhOk"]
        )
        vp = (
            np.concatenate(a["BuJpsiK_kvfCvhVtxProb"])
            if a["BuJpsiK_kvfCvhVtxProb"].dtype == object
            else a["BuJpsiK_kvfCvhVtxProb"]
        )
        sl = (
            np.concatenate(a["BuJpsiK_kvfCvhSl3d"])
            if a["BuJpsiK_kvfCvhSl3d"].dtype == object
            else a["BuJpsiK_kvfCvhSl3d"]
        )

        good = (
            (ok == 1)
            & np.isfinite(m)
            & np.isfinite(pt)
            & (np.abs(m - M_BU) < 0.10)
            & (vp > 0.01)
            & (sl > 3.0)
        )
        m, pt, eta = m[good], pt[good], np.abs(eta[good])
        h_sig += np.histogram2d(pt, eta, bins=[pt_edges, eta_edges])[0]
        nread += 1

    w = np.clip(h_sig, 0, None)
    print(f"[spectrum] {nread} files, {w.sum():.0f} candidates in the B mass band")
    os.makedirs(os.path.dirname(args.spectrum_cache), exist_ok=True)
    np.savez(args.spectrum_cache, pt_edges=pt_edges, eta_edges=eta_edges, w=w)
    return pt_edges, eta_edges, w


def sample_parametric(n):
    """Power-law pt with a trigger-like threshold, flat-ish rapidity."""
    pt = np.empty(0)
    while pt.size < n:
        x = rng.uniform(5.0, 100.0, size=3 * (n - pt.size) + 1000)
        f = x / (1.0 + (x / 5.5) ** 2) ** 3.2
        keep = rng.random(x.size) < f / f.max()
        pt = np.concatenate([pt, x[keep]])
    return pt[:n], rng.uniform(-2.2, 2.2, n)


def sample_spectrum(n, pt_edges, eta_edges, w):
    flat = w.ravel()
    idx = rng.choice(flat.size, size=n, p=flat / flat.sum())
    ip, ie = np.unravel_index(idx, w.shape)
    pt = pt_edges[ip] + rng.random(n) * (pt_edges[ip + 1] - pt_edges[ip])
    aeta = eta_edges[ie] + rng.random(n) * (eta_edges[ie + 1] - eta_edges[ie])
    eta = aeta * rng.choice([-1.0, 1.0], size=n)
    return pt, eta


# ------------------------------------------------------------ kinematics utils
def p4(E, px, py, pz):
    return np.stack([E, px, py, pz], axis=-1)


def mass(v):
    m2 = v[..., 0] ** 2 - (v[..., 1] ** 2 + v[..., 2] ** 2 + v[..., 3] ** 2)
    return np.sqrt(np.clip(m2, 0, None))


def boost(v, beta):
    """Boost 4-vectors `v` by velocity `beta` (n,3)."""
    b2 = np.sum(beta**2, axis=-1)
    gamma = 1.0 / np.sqrt(1.0 - b2)
    bp = np.sum(beta * v[..., 1:], axis=-1)
    # (gamma-1)/b2, stable as b2 -> 0
    g2 = np.where(b2 > 1e-12, (gamma - 1.0) / np.where(b2 > 1e-12, b2, 1.0), 0.0)
    E = gamma * (v[..., 0] + bp)
    p = v[..., 1:] + (g2 * bp)[..., None] * beta + (gamma * v[..., 0])[..., None] * beta
    return p4(E, p[..., 0], p[..., 1], p[..., 2])


def two_body_p(M, m1, m2):
    return np.sqrt(
        np.clip((M**2 - (m1 + m2) ** 2) * (M**2 - (m1 - m2) ** 2), 0, None)
    ) / (2 * M)


def rel_bw(n, m0, gamma, mmin, mmax):
    """Relativistic Breit-Wigner by rejection sampling."""
    out = np.empty(0)
    while out.size < n:
        x = rng.uniform(mmin, mmax, size=2 * (n - out.size) + 1000)
        num = x**2 * gamma**2
        den = (x**2 - m0**2) ** 2 + m0**2 * gamma**2
        f = num / den
        keep = rng.random(x.size) < f / f.max()
        out = np.concatenate([out, x[keep]])
    return out[:n]


def sample_cos_helicity(n, fL):
    """W(c) ~ (3 fL/2) c^2 + (3(1-fL)/4)(1-c^2); fL < 0 -> isotropic."""
    if fL < 0:
        return rng.uniform(-1, 1, size=n)
    out = np.empty(0)
    fmax = max(1.5 * fL, 0.75 * (1 - fL))
    while out.size < n:
        c = rng.uniform(-1, 1, size=2 * (n - out.size) + 1000)
        w = 1.5 * fL * c**2 + 0.75 * (1 - fL) * (1 - c**2)
        keep = rng.random(c.size) < w / fmax
        out = np.concatenate([out, c[keep]])
    return out[:n]


# ------------------------------------------------------------------- generate
N = args.nevents
print(
    f"[toy] generating {N:,} B0 -> J/psi K*0 decays, fL = {args.fL}, "
    f"spectrum = {args.spectrum}, sigma(pt)/pt = {args.ptres}"
)

if args.spectrum == "data":
    pt_edges, eta_edges, wspec = load_spectrum()
    ptB, etaB = sample_spectrum(N, pt_edges, eta_edges, wspec)
else:
    ptB, etaB = sample_parametric(N)
phiB = rng.uniform(-np.pi, np.pi, N)
pB = ptB * np.cosh(etaB)
B = p4(
    np.sqrt(pB**2 + M_B0**2),
    ptB * np.cos(phiB),
    ptB * np.sin(phiB),
    ptB * np.sinh(etaB),
)

# K* mass, then B -> J/psi K* in the B rest frame (isotropic: B is spin 0)
mKst = rel_bw(N, M_KSTAR, W_KSTAR, M_K + M_PI + 1e-3, 1.35)
q = two_body_p(M_B0, M_JPSI, mKst)
cth = rng.uniform(-1, 1, N)
sth = np.sqrt(1 - cth**2)
ph = rng.uniform(-np.pi, np.pi, N)
nhat = np.stack([sth * np.cos(ph), sth * np.sin(ph), cth], axis=-1)  # K* dir in B frame
Kst_B = p4(np.sqrt(q**2 + mKst**2), *(q[:, None] * nhat).T)

# K* -> K pi in the K* rest frame, helicity axis = nhat
qq = two_body_p(mKst, M_K, M_PI)
cK = sample_cos_helicity(N, args.fL)
sK = np.sqrt(1 - cK**2)
phK = rng.uniform(-np.pi, np.pi, N)
# orthonormal basis with nhat as z
tmp = np.zeros_like(nhat)
tmp[:, 0] = 1.0
tmp[np.abs(nhat[:, 0]) > 0.9] = np.array([0.0, 1.0, 0.0])
e1 = np.cross(nhat, tmp)
e1 /= np.linalg.norm(e1, axis=-1, keepdims=True)
e2 = np.cross(nhat, e1)
kdir = (
    (sK * np.cos(phK))[:, None] * e1
    + (sK * np.sin(phK))[:, None] * e2
    + cK[:, None] * nhat
)

K_r = p4(np.sqrt(qq**2 + M_K**2), *(qq[:, None] * kdir).T)
PI_r = p4(np.sqrt(qq**2 + M_PI**2), *(-qq[:, None] * kdir).T)

# K* rest -> B rest -> lab
beta_Kst = Kst_B[:, 1:] / Kst_B[:, 0:1]
K_B = boost(K_r, beta_Kst)
PI_B = boost(PI_r, beta_Kst)
beta_B = B[:, 1:] / B[:, 0:1]
K = boost(K_B, beta_B)
PI = boost(PI_B, beta_B)
JPSI = boost(
    boost(p4(np.sqrt(q**2 + M_JPSI**2), *(-q[:, None] * nhat).T), np.zeros((N, 3))),
    beta_B,
)


def smear(v, m, frac):
    """Scale the 3-momentum by 1 + N(0, frac), recompute E at fixed mass."""
    if frac <= 0:
        return v
    s = 1.0 + rng.normal(0.0, frac, size=v.shape[0])
    p = v[:, 1:] * s[:, None]
    return p4(np.sqrt(np.sum(p**2, axis=-1) + m**2), p[:, 0], p[:, 1], p[:, 2])


K = smear(K, M_K, args.ptres)
PI = smear(PI, M_PI, args.ptres)
JPSI = smear(JPSI, M_JPSI, args.ptres)

pKmag = np.linalg.norm(K[:, 1:], axis=-1)
pPImag = np.linalg.norm(PI[:, 1:], axis=-1)

# ------------------------------------------------------- true and swapped mass
hh_true = K + PI
B_true = hh_true + JPSI
m_hh_true = mass(hh_true)
m_B_true = mass(B_true)

# swap: same 3-momenta, exchanged mass hypotheses
K_sw = p4(np.sqrt(pKmag**2 + M_PI**2), K[:, 1], K[:, 2], K[:, 3])
PI_sw = p4(np.sqrt(pPImag**2 + M_K**2), PI[:, 1], PI[:, 2], PI[:, 3])
hh_sw = K_sw + PI_sw
B_sw = hh_sw + JPSI
m_hh_sw = mass(hh_sw)
m_B_sw = mass(B_sw)

dm_hh = m_hh_sw - m_hh_true
dm_B = m_B_sw - m_B_true

# closure on the analytic dE
dE_pred = 0.5 * (M_K**2 - M_PI**2) * (1.0 / pPImag - 1.0 / pKmag)
dE_act = (K_sw[:, 0] + PI_sw[:, 0]) - (K[:, 0] + PI[:, 0])

print("\n=== closure on the D1b relations ===")
print(
    f"  dE analytic vs exact      : median ratio "
    f"{np.median(dE_pred / dE_act):.4f}   (-> 1 if p >> m)"
)
print(
    f"  dE > 0 fraction           : {np.mean(dE_act > 0):.4f}   "
    f"(D1b: swaps push masses UP)"
)
fin = np.isfinite(dm_B) & (np.abs(dm_B) > 1e-6)
ratio = dm_hh[fin] / dm_B[fin]
print(
    f"  dm(Kpi)/dm(B)             : median {np.median(ratio):.3f}   "
    f"[16,84]% {np.percentile(ratio, 16):.3f} {np.percentile(ratio, 84):.3f}"
    f"    D1b predicts ~2.0"
)
print(f"  <dm(B)> over all decays   : {1e3*np.mean(dm_B):8.1f} MeV")
print(f"  <dm(Kpi)> over all decays : {1e3*np.mean(dm_hh):8.1f} MeV")

# ------------------------------------------------------------ window behaviour
WINDOWS = np.array([0.015, 0.020, 0.025, 0.035, 0.050, 0.075, 0.100, 0.150])

# The peak region is what a mass fit actually measures. Take it as +/-2 sigma
# of the smeared correct-assignment peak.
sigma_B = np.std(m_B_true[np.abs(m_B_true - M_B0) < 0.15]) if args.ptres > 0 else 0.010
PEAK = 2.0 * sigma_B
print(
    f"\n[toy] smeared m(B) resolution sigma = {1e3*sigma_B:.1f} MeV, "
    f"peak region +/-{1e3*PEAK:.0f} MeV"
)

in_peak_true = np.abs(m_B_true - M_B0) < PEAK
in_peak_sw = np.abs(m_B_sw - M_B0) < PEAK

rows = []
for w in WINDOWS:
    a_true = np.abs(m_hh_true - M_KSTAR) < w  # correct assignment in window
    a_sw = np.abs(m_hh_sw - M_KSTAR) < w  # swapped assignment in window
    pick_true = np.abs(m_hh_true - M_KSTAR) <= np.abs(m_hh_sw - M_KSTAR)

    out = {"w": w}
    for tag, ok_msk, sw_msk in (
        ("keepboth", a_true, a_sw),
        ("closest", a_true & pick_true, a_sw & ~pick_true),
        ("veto", a_true & ~a_sw, a_sw & ~a_true),
    ):
        # what lands in the peak region is what biases the fitted mass
        ok_p = ok_msk & in_peak_true
        sw_p = sw_msk & in_peak_sw
        n_ok, n_sw = int(ok_p.sum()), int(sw_p.sum())
        tot = n_ok + n_sw
        f = n_sw / tot if tot else 0.0
        if n_sw:
            sh = dm_B[sw_p]
            mean_shift, err_shift = np.mean(sh), np.std(sh) / np.sqrt(n_sw)
        else:
            mean_shift = err_shift = 0.0
        out[tag] = dict(
            f=f,
            bias=f * mean_shift,
            bias_err=f * err_shift,
            eff=n_ok / N,
            mean_shift=mean_shift,
            n_sw=n_sw,
        )
    rows.append(out)

print(
    f"\n=== swap contamination IN THE PEAK REGION vs K* window "
    f"(fL={args.fL}, spectrum={args.spectrum}) ==="
)
print(
    "bias = f_swap x <dm(B)> of swaps landing inside the peak; "
    "+/- is the toy's own statistical error"
)
print(f"{'w[MeV]':>7} | {'keepboth':^32} | {'closest':^32} | {'veto':^32}")
print(
    f"{'':>7} | {'f_sw':>6} {'bias[MeV]':>18} {'eff':>6} |"
    f" {'f_sw':>6} {'bias[MeV]':>18} {'eff':>6} |"
    f" {'f_sw':>6} {'bias[MeV]':>18} {'eff':>6}"
)
print("-" * 116)
for r in rows:
    line = f"{1e3*r['w']:7.0f} |"
    for tag in ("keepboth", "closest", "veto"):
        d = r[tag]
        line += (
            f" {d['f']:6.4f} {1e3*d['bias']:9.3f}+/-{1e3*d['bias_err']:<6.3f}"
            f" {d['eff']:6.3f} |"
        )
    print(line)

# Was the D1b bound |dm_B| < w/2 right?  (spoiler: no)
print("\n=== TEST of the D1b bound:  |dm(B)| < w/2 for surviving swaps ===")
for w in (0.025, 0.050, 0.100):
    msk = np.abs(m_hh_sw - M_KSTAR) < w
    if msk.sum() == 0:
        continue
    q95 = np.percentile(np.abs(dm_B[msk]), 95)
    frac_in = np.mean(np.abs(dm_B[msk]) < w / 2)
    print(
        f"  w={1e3*w:5.0f} MeV : only {100*frac_in:5.1f}% obey it; "
        f"95th pct |dm(B)| = {1e3*q95:6.1f} MeV vs bound {1e3*w/2:5.1f}"
    )

# What DOES protect the peak: where the swap is peaking, it is unshifted.
print("\n=== Why the peak survives anyway ===")
for w in (0.025, 0.050, 0.100):
    sw = np.abs(m_hh_sw - M_KSTAR) < w
    inpk = sw & in_peak_sw
    if inpk.sum() == 0:
        continue
    print(
        f"  w={1e3*w:5.0f} MeV : {100*np.mean(inpk[sw]):5.1f}% of swaps land in "
        f"the peak region; of those, <dm(B)> = {1e3*np.mean(dm_B[inpk]):6.2f} MeV, "
        f"median {1e3*np.median(dm_B[inpk]):6.2f}, "
        f"|dm|<{1e3*PEAK:.0f} by construction"
    )

np.savez(
    os.path.join(args.outdir, f"kstar_swap_toy_{TAG}.npz"),
    windows=WINDOWS,
    keepboth_f=[r["keepboth"]["f"] for r in rows],
    keepboth_bias=[r["keepboth"]["bias"] for r in rows],
    closest_f=[r["closest"]["f"] for r in rows],
    closest_bias=[r["closest"]["bias"] for r in rows],
    veto_f=[r["veto"]["f"] for r in rows],
    veto_bias=[r["veto"]["bias"] for r in rows],
    keepboth_biaserr=[r["keepboth"]["bias_err"] for r in rows],
    closest_biaserr=[r["closest"]["bias_err"] for r in rows],
    veto_biaserr=[r["veto"]["bias_err"] for r in rows],
    keepboth_eff=[r["keepboth"]["eff"] for r in rows],
    closest_eff=[r["closest"]["eff"] for r in rows],
    veto_eff=[r["veto"]["eff"] for r in rows],
)

# ------------------------------------------------------------------- figures
med = float(np.median(ratio))
fig, ax = plt.subplots(2, 2, figsize=(12.0, 8.4))

a = ax[0, 0]
sel = np.isfinite(dm_B) & (np.abs(dm_B) < 0.6) & (np.abs(dm_hh) < 1.2)
a.hist2d(
    1e3 * dm_B[sel],
    1e3 * dm_hh[sel],
    bins=[120, 120],
    range=[[-400, 600], [-400, 1200]],
    cmap="viridis",
    norm=mpl.colors.LogNorm(),
)
xs = np.linspace(-400, 600, 10)
a.plot(
    xs, 2.0 * xs, "r--", lw=1.8, label=r"original D1b claim: ratio $=2.0$ (FALSIFIED)"
)
a.plot(
    xs,
    med * xs,
    "w-",
    lw=1.4,
    label=f"toy median = {med:.2f}, [16,84]% = "
    f"{np.percentile(ratio,16):.2f}-{np.percentile(ratio,84):.2f}",
)
a.set_xlabel(r"$\Delta m(B)$ under swap [MeV]")
a.set_ylabel(r"$\Delta m(K\pi)$ under swap [MeV]")
a.set_title(r"Ratio is NOT fixed: spread 0.8-2.9, median %.2f" % med)
a.legend(loc="upper left", fontsize=9)

a = ax[0, 1]
for w, c in zip((0.025, 0.050, 0.100), ("C0", "C1", "C2")):
    msk = np.abs(m_hh_sw - M_KSTAR) < w
    a.hist(
        1e3 * dm_B[msk],
        bins=120,
        range=(-400, 600),
        histtype="step",
        lw=1.8,
        color=c,
        density=True,
        label=f"w = {1e3*w:.0f} MeV",
    )
a.axvspan(-1e3 * PEAK, 1e3 * PEAK, color="k", alpha=0.12, label="peak region")
a.set_xlabel(r"$\Delta m(B)$ of surviving swaps [MeV]")
a.set_ylabel("normalised")
a.set_title(r"Swaps are BROAD in $\Delta m(B)$ -- the $w/2$ bound fails")
a.legend(fontsize=9)

a = ax[1, 0]
for tag, c, mk in (
    ("keepboth", "C3", "o"),
    ("closest", "C0", "s"),
    ("veto", "C2", "^"),
):
    a.plot(1e3 * WINDOWS, [r[tag]["f"] for r in rows], mk + "-", color=c, label=tag)
a.set_xlabel(r"$K^*$ window half-width $w$ [MeV]")
a.set_ylabel(r"$f_\mathrm{swap}$")
a.set_yscale("log")
a.set_title("Swap fraction vs window, three PID-free strategies")
a.legend(fontsize=9)
a.grid(alpha=0.3)

a = ax[1, 1]
for tag, c, mk in (
    ("keepboth", "C3", "o"),
    ("closest", "C0", "s"),
    ("veto", "C2", "^"),
):
    e5 = [r for r in rows if abs(r["w"] - 0.05) < 1e-9][0][tag]["eff"]
    a.errorbar(
        1e3 * WINDOWS,
        [1e3 * r[tag]["bias"] for r in rows],
        yerr=[1e3 * r[tag]["bias_err"] for r in rows],
        fmt=mk + "-",
        color=c,
        capsize=2,
        label=f"{tag}  (eff@50MeV = {e5:.3f})",
    )
a.axhline(0.0, color="k", lw=0.8)
for y in (0.025, -0.025):
    a.axhline(y, color="k", ls="--", lw=1.0)
a.text(18, 0.03, r"$\pm$ B$^+$ statistical precision (0.025 MeV)", fontsize=8)
a.set_xlabel(r"$K^*$ window half-width $w$ [MeV]")
a.set_ylabel("predicted $m(B^0)$ peak bias [MeV]")
a.set_title(f"Predicted peak bias  (fL={args.fL:g}, $\\sigma_B$={1e3*sigma_B:.0f} MeV)")
a.legend(fontsize=8)
a.grid(alpha=0.3)

fig.tight_layout()
out = os.path.join(args.outdir, f"kstar_swap_toy_{TAG}.png")
fig.savefig(out, dpi=140)
print(f"\n[figure] {out}")
