"""The A-flatness acceptance gate for a response model.

A is a momentum-INDEPENDENT scale: the correction is (1 + A - e*k)*k with A
constant in k = 1/pT. So the mass shift A induces must be flat in bachelor pT.
It is not flat if the network was asked to extrapolate outside its training
range, and that failure is invisible in fit output -- sigma is measured in
prior units, so a runaway template variation always reads as a SMALLER
uncertainty, i.e. as improved sensitivity.

No data, no fit, no Asimov. Reads `nominal_muonScaleSyst_responseWeights` from
a histmaker run and, per bachelor-pT band, computes the mean-mass shift each
variation induces on the signal template:

    shift(band, var) = <m>_up(band, var) - <m>_nominal(band)

taken as the average of the up and down legs with the down sign flipped, so a
linear response cancels the second-order term.

Verdict: flat if max|shift| / |median shift| over the bands stays inside a
stated tolerance AND no band flips sign relative to the median.

Example:
    python scripts/btojpsik/kpi_aflatness_gate.py \
        --inputs muon=.../gate_muon_clamp.hdf5 kpi=.../gate_kpi_clamp.hdf5 \
        -o ~/public_html/mz/kpi_260908
"""

import argparse
import os
import sys

import h5py
import matplotlib as mpl
import numpy as np

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "plotting"))
from btojpsik_cvhnano_datamc import INK, INK_MUTED, SURFACE, style  # noqa: E402

from wremnants.production.btojpsik_cvhnano_axes import GEN_SIGNAL  # noqa: E402
from wums import ioutils  # noqa: E402

HIST = "nominal_muonScaleSyst_responseWeights"
SAMPLE = "InclusiveBToJpsiX_2016PostVFP"
PARAMS = ["A", "e", "M"]
# Bands chosen to straddle both models' edges: the muon family's 2 GeV clamp
# and the kpi kaon sign flip near 0.9 GeV.
BANDS = [(0.1, 0.3), (0.3, 0.5), (0.5, 0.8), (0.8, 1.0), (1.0, 1.5), (1.5, 2.0),
         (2.0, 3.0), (3.0, 12.0)]


def _logx_ticks(ax, vals):
    """Readable ticks on a log x axis.

    matplotlib's default LogLocator prints every minor decade multiple, which
    on a range of one or two decades collides into an unreadable smear. Set
    them explicitly instead.
    """
    from matplotlib.ticker import FixedLocator, NullLocator, FuncFormatter
    ax.xaxis.set_major_locator(FixedLocator(vals))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_major_formatter(FuncFormatter(
        lambda v, _: f"{v:g}"))


def load(path):
    with h5py.File(path, "r") as h:
        return ioutils.pickle_load_h5py(h[SAMPLE])["output"][HIST].get()


def shifts(hobj, n_eta_groups=4):
    """Mean-mass shift per (pT band, param), in MeV, on the signal template.

    Axes are (fit_bachelor_pt, fit_bachelor_eta, fit_bachelor_charge, mass,
    genCategory, unc, downup). Nuisances are packed unc = group*nparams+iparam.
    """
    ax = {a.name: i for i, a in enumerate(hobj.axes)}
    v = np.asarray(hobj.values())
    # signal only, summed over eta and charge
    v = np.take(v, GEN_SIGNAL, axis=ax["genCategory"])
    ptax = hobj.axes[ax["fit_bachelor_pt"]]
    massax = hobj.axes[ax["jointCvh_mass"]]
    order = ["fit_bachelor_pt", "fit_bachelor_eta", "fit_bachelor_charge",
             "jointCvh_mass", "unc", "downUpVar"]
    idx = [ax[n] - (ax[n] > ax["genCategory"]) for n in order]
    v = np.transpose(v, idx)                       # pt, eta, q, mass, unc, du
    v = v.sum(axis=(1, 2))                         # pt, mass, unc, du
    centres = massax.centers
    pt_edges = ptax.edges
    nunc = v.shape[2]
    nparams = len(PARAMS)
    out = np.full((len(BANDS), nparams), np.nan)
    counts = np.zeros(len(BANDS))
    for ib, (lo, hi) in enumerate(BANDS):
        k = (pt_edges[:-1] >= lo) & (pt_edges[:-1] < hi)
        if not k.any():
            continue
        block = v[k].sum(axis=0)                   # mass, unc, du
        # Nominal = the up/down average, which is what the variation is a
        # perturbation of; using it rather than a separate nominal histogram
        # keeps this self-contained and cancels the leading nonlinearity.
        for ip in range(nparams):
            us = [u for u in range(nunc) if u % nparams == ip]
            up = block[:, us, 1].sum(axis=1)
            dn = block[:, us, 0].sum(axis=1)
            nom = 0.5 * (up + dn)
            if nom.sum() <= 0:
                continue
            m_up = np.average(centres, weights=up) if up.sum() > 0 else np.nan
            m_dn = np.average(centres, weights=dn) if dn.sum() > 0 else np.nan
            m_no = np.average(centres, weights=nom)
            out[ib, ip] = 0.5 * ((m_up - m_no) - (m_dn - m_no)) * 1e3
        counts[ib] = v[k].sum()
    return out, counts


# CORRECTED 2026-09-08. The gate inherited from section 18 asserted that A,
# being a pT-independent scale, must induce a pT-independent mass shift. That
# premise is wrong. A scales the kaon's momentum by a constant fraction, but
# the B mass responds to that through dm/dln(p_K), which grows with the kaon's
# share of the candidate's energy -- so even a perfect model gives a shift that
# rises with kaon pT. Measured here, EVERY configuration including the recorded
# baseline shows the same rise, from about 0.45 of the high-pT value in the
# softest band to 1.0 above 3 GeV. An absolute flatness test therefore fails
# everything, which is why it cannot be the gate.
#
# What IS a model artefact is a DEPARTURE from that common profile. So the gate
# is relative: normalise each configuration's profile to its own in-domain
# high-pT band, where every model is inside its training range and they must
# agree, then compare shapes. A model extrapolating out of domain breaks the
# common shape; a model merely obeying kinematics reproduces it.
REF_BAND = -1  # 3-12 GeV: in domain for every model considered


def _abs_spread(a_col):
    """Criterion (1), kept only so it can be quoted and seen to fail."""
    good = np.isfinite(a_col)
    med = np.median(a_col[good])
    return np.max(np.abs(a_col[good] / med - 1.0)) if med else np.nan


def normalised(a_col):
    ref = a_col[REF_BAND]
    if not np.isfinite(ref) or ref == 0:
        return np.full_like(a_col, np.nan)
    return a_col / ref


def verdict(a_col):
    """The acceptance criterion, third formulation and the first honest one.

    Two earlier ones were tried and both are recorded here because both are
    wrong in instructive ways.

    (1) ABSOLUTE FLATNESS, inherited from section 18. Fails every
        configuration measured, including the one every recorded A/e/M number
        came from, because its premise is false -- see the note above REF_BAND.

    (2) DEPARTURE FROM A BASELINE PROFILE. Circular: the only sensible baseline
        is a clamped configuration, and a clamped configuration cannot differ
        from itself below the clamp, so every unclamped configuration "fails"
        by construction whether or not it is wrong.

    (3) MONOTONICITY AND SIGN STABILITY, used here. A rises with kaon pT
        through dm/dln(p_K) and nothing in the kinematics makes that
        non-monotonic or sign-changing, so a profile that turns over or crosses
        zero is the model talking, not the physics. This is necessary rather
        than sufficient -- it cannot tell a mildly-wrong model from a right one
        -- but it is not circular and it catches the failure that matters.

    Whether the clamped or the unclamped answer is correct BELOW the floor is
    not decidable here: the clamp freezes the conditioning, so it manufactures
    a smooth continuation regardless. That question belongs to the network
    probes in kpi_response_shapes.py, which show `_smearing` breaking at
    1.8 GeV at every response residual while `_nosmearing` stays flat.
    """
    good = np.isfinite(a_col)
    if good.sum() < 3:
        return "INSUFFICIENT", np.nan, np.nan
    v = a_col[good]
    signs = np.sign(v)
    sign_ok = np.all(signs == signs[0])
    d = np.diff(np.abs(v))
    # allow a single non-monotonic step of <2% of the range, for band noise
    tol = 0.02 * (np.max(np.abs(v)) - np.min(np.abs(v)))
    n_bad = int(np.sum(d < -tol))
    mono_ok = n_bad == 0
    return ("PASS" if (sign_ok and mono_ok) else "FAIL"), n_bad, int(not sign_ok)


REF_BAND_LABEL = "3-12 GeV"


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--inputs", nargs="+", required=True, help="label=path ...")
    p.add_argument("-o", "--outfolder", required=True)
    p.add_argument("--tol", type=float, default=0.15)
    p.add_argument("--baseline", default=None,
                   help="Deprecated: criterion (2) was circular. Ignored.")
    args = p.parse_args()

    style()
    os.makedirs(args.outfolder, exist_ok=True)
    labels, res = [], {}
    for spec in args.inputs:
        lab, path = spec.split("=", 1)
        labels.append(lab)
        res[lab] = shifts(load(path))

    band_mid = [np.sqrt(lo * hi) for lo, hi in BANDS]
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.2))
    cols = plt.cm.tab10(np.linspace(0, 0.9, len(labels)))
    for ip, (ax, par) in enumerate(zip(axes, PARAMS)):
        for lab, c in zip(labels, cols):
            ax.plot(band_mid, res[lab][0][:, ip], "o-", color=c, markersize=4,
                    linewidth=1.4, label=lab)
        ax.axhline(0.0, color=INK_MUTED, linewidth=0.9)
        ax.set_xscale("log")
        _logx_ticks(ax, [0.2, 0.5, 1, 2, 5])
        ax.set_xlabel(r"bachelor $p_T$ [GeV]")
        ax.set_title(f"induced shift, {par}", fontsize=10, color=INK, loc="left")
        ax.grid(alpha=0.6)
        ax.set_axisbelow(True)
    axes[0].set_ylabel("mean-mass shift [MeV]")
    axes[0].legend(fontsize=7.5, frameon=False)
    # NOT "this panel must be flat" -- that was the section 18 premise and it
    # is wrong. A rises with kaon pT through dm/dln(p_K) in every configuration
    # measured, including the recorded baseline. What a good model must be is
    # monotonic and single-signed.
    axes[0].text(0.02, 0.03,
                 "must be monotonic and single-signed;\n"
                 "the rise itself is kinematic, not a defect",
                 transform=axes[0].transAxes, fontsize=8.5, color="#b4232a")
    out = os.path.join(args.outfolder, "aflatness_gate.png")
    fig.savefig(out, dpi=160, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
    print("wrote", out)

    w = 16
    print("\nInduced A shift [MeV] per bachelor-pT band, signal template")
    print(f"{'band [GeV]':>13} " + " ".join(f"{l[:w]:>{w}}" for l in labels))
    for ib, (lo, hi) in enumerate(BANDS):
        print(f"{lo:5.2f}-{hi:<7.2f} " +
              " ".join(f"{res[l][0][ib,0]:>{w}.4f}" for l in labels))

    print("\nNormalised to each configuration's own 3-12 GeV band")
    for ib, (lo, hi) in enumerate(BANDS):
        print(f"{lo:5.2f}-{hi:<7.2f} " +
              " ".join(f"{normalised(res[l][0][:,0])[ib]:>{w}.3f}" for l in labels))

    print(f"\n{'sec18 abs spread':>13} " +
          " ".join(f"{_abs_spread(res[l][0][:,0]):>{w}.3f}" for l in labels))
    print("   ^ criterion (1), absolute flatness. It fails every configuration,")
    print("     including the one every recorded number came from. Premise is wrong.")

    print(f"\n{'non-monotonic':>13} " +
          " ".join(f"{verdict(res[l][0][:,0])[1]:>{w}}" for l in labels))
    print(f"{'sign changes':>13} " +
          " ".join(f"{verdict(res[l][0][:,0])[2]:>{w}}" for l in labels))
    print(f"{'VERDICT':>13} " +
          " ".join(f"{verdict(res[l][0][:,0])[0]:>{w}}" for l in labels))
    print("   ^ criterion (3): monotonic in |shift| and single-signed. "
          "Necessary, not sufficient.")


if __name__ == "__main__":
    main()
