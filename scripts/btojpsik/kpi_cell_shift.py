"""Section 2.4: the induced B-mass shift per fit cell, per response model.

The network probes in kpi_response_shapes.py say what the model does; this says
what the FIT sees. Same quantity as the acceptance gate, but binned the way the
tensor is -- bachelor pT x eta -- rather than in wide pT bands, because a
variation that is well behaved inclusively can still be disordered cell by cell,
and the cells are what carry the nuisances.

Example:
    python scripts/btojpsik/kpi_cell_shift.py -o out/ \\
        --inputs muon_rec=.../kpiA_muon_rec.hdf5 kpi=.../kpiA_kpi.hdf5
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


def cell_shifts(path, eta_rebin=6):
    """Mean-mass shift in MeV on a (pT, eta) grid, per parameter."""
    with h5py.File(path, "r") as h:
        obj = ioutils.pickle_load_h5py(h[SAMPLE])["output"][HIST].get()
    ax = {a.name: i for i, a in enumerate(obj.axes)}
    v = np.take(np.asarray(obj.values()), GEN_SIGNAL, axis=ax["genCategory"])
    order = ["fit_bachelor_pt", "fit_bachelor_eta", "fit_bachelor_charge",
             "jointCvh_mass", "unc", "downUpVar"]
    idx = [ax[n] - (ax[n] > ax["genCategory"]) for n in order]
    v = np.transpose(v, idx).sum(axis=2)          # pt, eta, mass, unc, du
    centres = obj.axes[ax["jointCvh_mass"]].centers
    npt, neta = v.shape[0], v.shape[1]
    ne = neta // eta_rebin
    v = v[:, : ne * eta_rebin].reshape(npt, ne, eta_rebin, *v.shape[2:]).sum(axis=2)
    nunc = v.shape[-2]
    out = np.full((npt, ne, len(PARAMS)), np.nan)
    for ip in range(len(PARAMS)):
        us = [u for u in range(nunc) if u % len(PARAMS) == ip]
        up = v[..., us, 1].sum(axis=-1)
        dn = v[..., us, 0].sum(axis=-1)
        nom = 0.5 * (up + dn)
        for i in range(npt):
            for j in range(ne):
                if nom[i, j].sum() <= 0:
                    continue
                mu = np.average(centres, weights=up[i, j])
                md = np.average(centres, weights=dn[i, j])
                mn = np.average(centres, weights=nom[i, j])
                out[i, j, ip] = 0.5 * ((mu - mn) - (md - mn)) * 1e3
    pt_edges = obj.axes[ax["fit_bachelor_pt"]].edges
    return out, pt_edges


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--inputs", nargs="+", required=True)
    p.add_argument("-o", "--outfolder", required=True)
    p.add_argument("--etaRebin", type=int, default=6)
    a = p.parse_args()

    style()
    os.makedirs(a.outfolder, exist_ok=True)
    labels, res = [], {}
    for spec in a.inputs:
        lab, path = spec.split("=", 1)
        labels.append(lab)
        res[lab] = cell_shifts(path, a.etaRebin)

    ncols = len(labels)
    for ip, par in enumerate(PARAMS):
        vals = np.concatenate([res[l][0][..., ip].ravel() for l in labels])
        vals = vals[np.isfinite(vals)]
        if not vals.size:
            continue
        lim = np.nanpercentile(np.abs(vals), 99)
        fig, axes = plt.subplots(1, ncols, figsize=(3.4 * ncols + 1.2, 3.6),
                                 sharey=True)
        axes = np.atleast_1d(axes)
        for ax, lab in zip(axes, labels):
            m, pt_edges = res[lab]
            im = ax.pcolormesh(pt_edges, np.arange(m.shape[1] + 1),
                               m[..., ip].T, cmap="RdBu_r",
                               vmin=-lim, vmax=lim, shading="flat")
            ax.set_title(lab, fontsize=9.5, color=INK, loc="left")
            ax.set_xlabel(r"bachelor $p_T$ [GeV]")
        axes[0].set_ylabel(f"eta group (rebin {a.etaRebin})")
        fig.colorbar(im, ax=axes.tolist(), label=f"{par} shift [MeV]",
                     fraction=0.02, pad=0.01)
        out = os.path.join(a.outfolder, f"cell_shift_{par}.png")
        fig.savefig(out, dpi=160, bbox_inches="tight", facecolor=SURFACE)
        plt.close(fig)
        print("wrote", out)

    print("\nCell-level summary: mean |shift| and the fraction of cells whose "
          "sign\ndisagrees with the majority -- a disordered variation is one "
          "a fit cannot use.")
    print(f"{'':<14}" + "".join(f"{p_:>26}" for p_ in PARAMS))
    print(f"{'':<14}" + "".join(f"{'mean|s|':>13}{'sign flip':>13}"
                                for _ in PARAMS))
    for lab in labels:
        m = res[lab][0]
        row = ""
        for ip in range(len(PARAMS)):
            x = m[..., ip][np.isfinite(m[..., ip])]
            if not x.size:
                row += f"{'--':>13}{'--':>13}"
                continue
            maj = np.sign(np.median(x))
            row += f"{np.mean(np.abs(x)):>13.4f}{np.mean(np.sign(x) != maj):>13.1%}"
        print(f"{lab:<14}" + row)


if __name__ == "__main__":
    main()
