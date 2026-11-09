#!/usr/bin/env python3
"""Reco kaon pT spectrum at full statistics, with the 1 GeV fit cut marked.

Reads the soft-kaon histogram (`--bachelorMinPt 0.1`, every other selection
nominal) and shows, per 0.1 GeV bin of reco kaon pT:

  top     data candidates in the fit mass window -- the combinatorial the soft
          region would cost;
  bottom  truth signal in the simulation, on its own scale.

The two are shown separately on purpose. The measured data/MC signal ratio
(18.24) is an integral above 1 GeV; the simulated kaon pT spectrum is harder
than the data, so applying that single number per bin overshoots above ~5 GeV
and would make a per-bin S/B meaningless. It is used here only for the
integrals quoted in the annotations.

Counts are candidate counts in 5.05 < m(mumuK) < 5.45 GeV.
"""
import argparse

import h5py
import matplotlib as mpl
import numpy as np

mpl.use("Agg")
import matplotlib.pyplot as plt
import narf.ioutils as io

HIST = "/ceph/submit/data/user/p/pmlugato/mz/histograms/btojpsik_cvhnano_softk_noclamp.hdf5"
MC_DS = "InclusiveBToJpsiX_2016PostVFP"
DATA_DS = "Charmonium_2016PostVFP_cvhNano"
GEN_SIGNAL, GEN_DATA = 0, 7

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", default=HIST)
ap.add_argument("-o", "--output", default="/home/submit/pmlugato/public_html/kaon_pt_spectrum_fullstats.png")
ap.add_argument("--signalScale", type=float, default=18.2433,
                help="data/MC signal ratio, measured as an integral above 1 GeV")
ap.add_argument("--ptMax", type=float, default=10.0)
ap.add_argument("--logy", action="store_true")
args = ap.parse_args()

with h5py.File(args.input, "r") as f:
    hmc = io.pickle_load_h5py(f[MC_DS])["output"]["fitbins"].get()
    hda = io.pickle_load_h5py(f[DATA_DS])["output"]["fitbins"].get()

pt_edges = hmc.axes["fit_bachelor_pt"].edges
m_edges = hmc.axes["jointCvh_mass"].edges
win = slice(int(np.searchsorted(m_edges, 5.05)), int(np.searchsorted(m_edges, 5.45)))


def project(h, cat):
    """Sum eta, charge and the mass window; the last pT entry is the >10 GeV overflow."""
    v = h.values(flow=True)[..., cat]
    e = h.variances(flow=True)[..., cat] if h.variances() is not None else v
    return v[:, :, :, win].sum(axis=(1, 2, 3)), e[:, :, :, win].sum(axis=(1, 2, 3))


sig, sig_var = project(hmc, GEN_SIGNAL)
dat, _ = project(hda, GEN_DATA)

ctr = 0.5 * (pt_edges[:-1] + pt_edges[1:])
n_pt = len(ctr)
below = ctr < 1.0
n_below = sig[:n_pt][below].sum()
n_above = sig[:n_pt][~below].sum() + sig[n_pt]  # keep the >10 GeV overflow in the fraction
frac = n_below / (n_below + n_above)
R = args.signalScale

fig, (top, bot) = plt.subplots(
    2, 1, figsize=(9, 7.4), sharex=True, gridspec_kw=dict(height_ratios=[1, 1], hspace=0.08)
)

top.stairs(dat[:n_pt], pt_edges, color="0.35", lw=1.8, fill=True, alpha=0.30)
top.stairs(dat[:n_pt], pt_edges, color="0.25", lw=1.8)
top.set_ylabel("data candidates / 0.1 GeV")
top.set_title("Reco kaon $p_T$ in $5.05 < m(\\mu\\mu K) < 5.45$ GeV — full statistics, 2016 postVFP CVH NanoAOD",
              fontsize=10.5, pad=12)
top.text(0.42, 0.93,
         f"data, all candidates\n"
         f"  $p_T < 1$ GeV : {dat[:n_pt][below].sum():12,.0f}\n"
         f"  $p_T > 1$ GeV : {dat[:n_pt][~below].sum() + dat[n_pt]:12,.0f}",
         transform=top.transAxes, va="top", fontsize=9.5, family="monospace")

bot.errorbar(ctr, sig[:n_pt], yerr=np.sqrt(sig_var[:n_pt]), fmt="o", ms=3.4, lw=1.1, color="#2f6fba")
bot.set_ylabel("simulated $B^+\\to J/\\psi K^+$ / 0.1 GeV")
bot.set_xlabel("reco kaon $p_T$ [GeV]")
bot.text(0.42, 0.40 if args.logy else 0.93,
         f"truth signal (simulation)\n"
         f"  $p_T < 1$ GeV : {n_below:10,.0f}   ({frac:.1%})\n"
         f"  $p_T > 1$ GeV : {n_above:10,.0f}\n"
         f"data-equivalent below 1 GeV\n"
         f"  $\\times${R:.2f} : {n_below * R:12,.0f}",
         transform=bot.transAxes, va="top", fontsize=9.5, family="monospace")

for ax in (top, bot):
    ax.axvspan(0, 1.0, color="#d9534f", alpha=0.09, lw=0)
    ax.axvline(1.0, color="#d9534f", lw=1.6, ls="--")
    ax.grid(alpha=0.25, lw=0.6)
    ax.set_axisbelow(True)
    ax.set_xlim(0, args.ptMax)
    if args.logy:
        ax.set_yscale("log")
    else:
        ax.set_ylim(0, None)
bot.annotate("fit selection: $p_T > 1$ GeV", xy=(1.08, 0.04), xycoords=("data", "axes fraction"),
             color="#d9534f", fontsize=9.5)

fig.savefig(args.output, dpi=140, bbox_inches="tight")
print("wrote", args.output)
print(f"signal below 1 GeV: MC {n_below:.0f}  data-equivalent {n_below * R:,.0f}  ({frac:.2%})")
print(f"signal above 1 GeV: MC {n_above:.0f}  data-equivalent {n_above * R:,.0f}")
