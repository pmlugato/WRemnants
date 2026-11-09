#!/usr/bin/env python3
"""Postfit B+ mass in every (pT, eta, charge) cell from a rabbit fit result that
saved the Project mapping over fit_bachelor_pt fit_bachelor_eta fit_bachelor_charge
jointCvh_mass. One PNG per cell with a stacked prediction, data, ratio panel and a
per-cell chi2 = sum (d - p)^2 / (sigma_d^2 + sigma_p^2) over the mass bins, plus a
summary of chi2 per cell and an index page. Run inside the container.
"""
import argparse, os, sys
import numpy as np
sys.path.insert(0, "rabbit"); sys.path.insert(0, "wums")
from rabbit import io_tools
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import chi2 as chi2dist

ap = argparse.ArgumentParser()
ap.add_argument("fit"); ap.add_argument("-o", "--outdir", required=True)
ap.add_argument("--channel", default="btojpsik_cvhnano")
ap.add_argument("--procs", nargs="+", default=["BuJpsiKBkg", "BuJpsiKMCPeak", "BuJpsiK"], help="stack order, bottom first")
ap.add_argument("--labels", nargs="+", default=["combinatorial (exp x quad)", "peaking MC template", "B+ signal"])
ap.add_argument("--colors", nargs="+", default=["#eda100", "#4a3aa7", "#1baf7a"])
ap.add_argument("--dataName", default="Data")
ap.add_argument("--title", default="")
a = ap.parse_args(); os.makedirs(a.outdir, exist_ok=True)
INK2, GRID = "#52514e", "#e6e6e3"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True, "legend.frameon": False})

r = io_tools.get_fitresult(a.fit)
key = [k for k in r["mappings"] if k.startswith("Project") and "fit_bachelor_pt" in k and "fit_bachelor_eta" in k and "fit_bachelor_charge" in k][0]
# The postfit prediction covariance (saved by rabbit_fit --computeHistCov) makes the
# per-cell chi2 a proper generalised chi2 over the cell's 20 mass bins, res^T
# (C_pred + C_data)^-1 res, i.e. rabbit's own non-profiled definition
# (fitter._residuals / _chi2) restricted to one cell. Without it the script falls
# back to the diagonal, which is 6% lower in total. NOTE rabbit's *postfit* mapping
# chi2 uses the PROFILED residual covariance instead (it accounts for the prediction
# having been fitted to this same data), which is ~1.37x larger and is not saved
# per bin, so it cannot be decomposed per cell.
cov_full = None
if "hist_postfit_inclusive_cov" in r["mappings"][key]:
    cov_full = np.asarray(r["mappings"][key]["hist_postfit_inclusive_cov"].get().values())
ch = r["mappings"][key]["channels"][a.channel]
hd = ch["hist_data_obs"].get(); hi = ch["hist_postfit_inclusive"].get(); hp = ch["hist_postfit"].get()
procs_axis = list(hp.axes["processes"])
pt_e = hi.axes["fit_bachelor_pt"].edges; eta_e = hi.axes["fit_bachelor_eta"].edges; q_e = hi.axes["fit_bachelor_charge"].edges; m_e = hi.axes["jointCvh_mass"].edges
mc = 0.5 * (m_e[1:] + m_e[:-1]); w = np.diff(m_e)
npt, neta, nq, nm = len(pt_e) - 1, len(eta_e) - 1, len(q_e) - 1, len(m_e) - 1
flat_idx = np.arange(npt * neta * nq * nm).reshape(npt, neta, nq, nm)
dv_flat = np.asarray(hd.variances()).ravel()
rows = []
for ipt in range(len(pt_e) - 1):
    for ieta in range(len(eta_e) - 1):
        for iq in range(len(q_e) - 1):
            sel = {"fit_bachelor_pt": ipt, "fit_bachelor_eta": ieta, "fit_bachelor_charge": iq}
            d = hd[sel].values(); dv = hd[sel].variances(); p = hi[sel].values(); pv = hi[sel].variances()
            if pv is None: pv = np.zeros_like(p)
            stacks = [hp[dict(sel, processes=procs_axis.index(pr))].values() for pr in a.procs]
            nb = len(d)
            if cov_full is not None:
                ii = flat_idx[ipt, ieta, iq]
                C = cov_full[np.ix_(ii, ii)] + np.diag(dv_flat[ii])
                chi2 = float((d - p) @ np.linalg.solve(C, d - p))
            else:
                chi2 = float(np.sum((d - p) ** 2 / np.where(dv + pv > 0, dv + pv, 1.0)))
            pval = float(chi2dist.sf(chi2, nb))
            q = "+" if q_e[iq] >= 0 else "-"
            label = f"kaon pT [{pt_e[ipt]:g}, {pt_e[ipt+1]:g}] GeV\nkaon eta [{eta_e[ieta]:+.1f}, {eta_e[ieta+1]:+.1f}]\ncharge {q}"
            rows.append((ipt, ieta, iq, chi2, nb, pval, float(d.sum()), float(stacks[-1].sum())))
            fig, (ax, axr) = plt.subplots(2, 1, figsize=(6.4, 5.6), sharex=True, gridspec_kw={"height_ratios": [3, 1.1], "hspace": 0.05})
            bottom = np.zeros_like(p)
            for s, lab, c in zip(stacks, a.labels, a.colors):
                ax.bar(mc, s, width=w, bottom=bottom, color=c, label=lab, lw=0); bottom = bottom + s
            ax.errorbar(mc, d, yerr=np.sqrt(dv), fmt="o", color="black", ms=3.5, lw=1.2, label=a.dataName, zorder=5)
            ax.set_ylim(0, 1.45 * max(d.max(), p.max())); ax.set_ylabel("events / 20 MeV")
            ax.legend(loc="upper right", fontsize=8.5)
            ax.text(0.02, 0.97, f"postfit  {a.title}\n{label}\nχ² / bins = {chi2:.1f} / {nb}  (p = {pval:.2f})", transform=ax.transAxes, va="top", ha="left", fontsize=8.5, color=INK2)
            ratio = np.where(p > 0, d / p, np.nan); rerr = np.where(p > 0, np.sqrt(dv) / p, np.nan); band = np.where(p > 0, np.sqrt(pv) / p, 0)
            axr.fill_between(mc, 1 - band, 1 + band, step="mid", color="#d9d9d6", label="pred. unc.")
            axr.errorbar(mc, ratio, yerr=rerr, fmt="o", color="black", ms=3.5, lw=1.2); axr.axhline(1, color=INK2, lw=0.8)
            axr.set_ylim(0.9, 1.1); axr.set_ylabel(f"{a.dataName} / pred."); axr.set_xlabel("m(μμK), joint CVH refit [GeV]"); axr.legend(loc="upper left", fontsize=8)
            fig.savefig(os.path.join(a.outdir, f"cell_pt{ipt}_eta{ieta}_q{iq}.png"), dpi=120, bbox_inches="tight"); plt.close(fig)

# summary: chi2 per cell as a table + heat map, and the index
rows.sort(); c2 = np.array([x[3] for x in rows]); nb = rows[0][4]
kind = "generalised, 20x20 covariance block (C_pred + C_data)" if cov_full is not None else "diagonal (C_pred + C_data), no off-diagonals"
L = [f"# Per-cell chi2 (data vs postfit, {nb} mass bins each, {kind}; no ndf correction for shared parameters)", "",
     f"cells {len(rows)}, chi2 per cell: mean {c2.mean():.1f}, median {np.median(c2):.1f}, min {c2.min():.1f}, max {c2.max():.1f}; total {c2.sum():.0f} over {len(rows)*nb} bins", "",
     "| pT bin | eta group | charge | chi2 | p | data | fitted signal |", "|--:|--:|:-:|--:|--:|--:|--:|"]
for ipt, ieta, iq, chi2, n, pval, nd, ns in rows: L.append(f"| {ipt} | {ieta} | {'+' if iq else '-'} | {chi2:.1f} | {pval:.2f} | {nd:.0f} | {ns:.0f} |")
open(os.path.join(a.outdir, "chi2_summary.md"), "w").write("\n".join(L) + "\n")
fig, axs = plt.subplots(1, 2, figsize=(12, 4), sharey=True)
for iq, ax in enumerate(axs):
    M = np.array([[next(x[3] for x in rows if x[:3] == (ipt, ieta, iq)) for ieta in range(24)] for ipt in range(3)])
    im = ax.imshow(M, aspect="auto", cmap="viridis", vmin=0, vmax=max(60, c2.max())); ax.set_title(f"χ² per cell ({nb} bins), charge {'+' if iq else '-'}", fontsize=10)
    ax.set_xlabel("eta group"); ax.set_yticks(range(3)); ax.set_yticklabels([f"pT {pt_e[i]:g}-{pt_e[i+1]:g}" for i in range(3)])
    ax.set_xticks(range(0, 24, 2))
fig.colorbar(im, ax=axs, fraction=0.02); fig.savefig(os.path.join(a.outdir, "chi2_map.png"), dpi=130, bbox_inches="tight")
H = ["<html><head><meta charset='utf-8'><title>postfit per cell</title><style>body{font-family:sans-serif;max-width:1600px;margin:auto} .grid{display:grid;grid-template-columns:repeat(4,1fr);gap:6px} img{max-width:100%;border:1px solid #ddd} pre{background:#f6f6f4;padding:10px;white-space:pre-wrap;font-size:12px}</style></head><body>",
     f"<h1>Postfit B⁺ mass per kinematic cell {a.title}</h1><figure><img src='chi2_map.png'></figure><pre>" + "\n".join(L[:3]) + "</pre>"]
for iq in (1, 0):
    for ipt in range(3):
        H.append(f"<h2>charge {'+' if iq else '-'}, kaon pT [{pt_e[ipt]:g}, {pt_e[ipt+1]:g}] GeV</h2><div class='grid'>")
        H += [f"<figure><img src='cell_pt{ipt}_eta{ieta}_q{iq}.png'></figure>" for ieta in range(24)]
        H.append("</div>")
open(os.path.join(a.outdir, "index.html"), "w").write("\n".join(H) + "</body></html>")
print("\n".join(L[:3])); print("wrote", len(rows), "cell plots to", a.outdir)
