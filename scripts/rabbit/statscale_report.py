#!/usr/bin/env python3
"""Tabulate what the B+ statistics scaling does to the joint mZ fits.

For each fit result: sigma(mZ); the conditional chain sigma(mZ | e fixed),
(| A fixed), (| A,e fixed), (| A,e,M fixed) and the implied marginal group
impacts; mean posterior sigma of A_i / e_i / M_i in prior units; mean
rho(A_i, e_i), rho(mZ, A_i), rho(mZ, e_i). Writes a markdown table and, with
--plot, one figure of sigma(mZ) and the scale-group impacts against the data
scale. Run inside the container (needs rabbit + wums).

    python scripts/rabbit/statscale_report.py --fits TAG=path ... -o report.md
"""
import argparse, os, sys
import numpy as np
sys.path.insert(0, "rabbit"); sys.path.insert(0, "wums")
from rabbit import io_tools

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--fits", nargs="+", required=True, help="LABEL=fitresult.hdf5, in display order")
ap.add_argument("--nEta", type=int, default=24)
ap.add_argument("--poi", default="massShiftZ100MeV", help="parameter of interest (unconstrained NOI)")
ap.add_argument("--poiScale", type=float, default=100.0, help="MeV per unit of the POI")
ap.add_argument("--poiLabel", default="mZ")
ap.add_argument("-o", "--output", required=True, help="markdown file to write")
ap.add_argument("--plot", default=None, help="optional png with sigma(mZ) and scale impacts")
args = ap.parse_args()


def analyse(path, n_eta, poi=None, scale=100.0):
    poi = poi or args.poi
    r = io_tools.get_fitresult(path)
    names = np.array([str(x) for x in r["parms"].get().axes["parms"]])
    cov = np.asarray(r["cov"].get().values())
    iz = int(np.flatnonzero(names == poi)[0])
    res_idx = np.array([i for i, n in enumerate(names) if n.startswith("sigRes_")])
    idx = {p: np.array([int(np.flatnonzero(names == f"{p}{i}")[0]) for i in range(n_eta)]) for p in "AeM"}
    sig = np.sqrt(np.diag(cov)); rho = cov / np.outer(sig, sig)

    def cond(fix):
        c = cov[iz, fix]; V = cov[np.ix_(fix, fix)]
        return scale * np.sqrt(cov[iz, iz] - c @ np.linalg.solve(V, c))

    s = scale * sig[iz]
    out = dict(sigma_mZ=s, edm=float(r["edmval"]) if "edmval" in r else np.nan)
    chain = {"e": idx["e"], "A": idx["A"], "A+e": np.r_[idx["A"], idx["e"]], "A+e+M": np.r_[idx["A"], idx["e"], idx["M"]]}
    if len(res_idx): chain["res"] = res_idx; chain["A+e+M+res"] = np.r_[idx["A"], idx["e"], idx["M"], res_idx]
    for k, fix in chain.items():
        sc = cond(fix); out[f"cond_{k}"] = sc; out[f"imp_{k}"] = np.sqrt(max(s * s - sc * sc, 0.0))
    for p in "AeM":
        out[f"sig_{p}"] = float(sig[idx[p]].mean()); out[f"rho_mZ_{p}"] = float(rho[iz, idx[p]].mean())
    out["rho_Ae"] = float(np.mean([rho[a, e] for a, e in zip(idx["A"], idx["e"])]))
    # bin-by-bin group impact if stored
    try:
        gi = r["impacts_grouped"].get()
        gl = [str(x) for x in gi.axes["impacts"]]; pi = [str(x) for x in gi.axes["parms"]].index(poi)
        for g in ("binByBinStat", "stat", "theory", "resolution"):
            if g in gl: out[f"gimp_{g}"] = scale * float(gi.values()[pi, gl.index(g)])
    except Exception:
        pass
    return out


rows = []
for spec in args.fits:
    label, path = spec.split("=", 1)
    if not os.path.exists(path):
        print(f"missing {path}, skipped"); continue
    rows.append((label, analyse(path, args.nEta, args.poi, args.poiScale)))

P = args.poiLabel
L = [f"| fit | σ({P}) [MeV] | σ({P} \\| e fixed) | σ({P} \\| A fixed) | σ({P} \\| A,e fixed) | σ({P} \\| A,e,M fixed) | σ({P} \\| res fixed) | σ({P} \\| A,e,M,res fixed) | impact e | impact A | impact A+e+M | impact res | Simulation stat. | Data stat. | Theory |",
     "|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|"]
for label, o in rows:
    L.append(f"| {label} | {o['sigma_mZ']:.4f} | {o['cond_e']:.4f} | {o['cond_A']:.4f} | {o['cond_A+e']:.4f} | {o['cond_A+e+M']:.4f} | {o.get('cond_res', np.nan):.4f} | {o.get('cond_A+e+M+res', np.nan):.4f} | "
             f"{o['imp_e']:.4f} | {o['imp_A']:.4f} | {o['imp_A+e+M']:.4f} | {o.get('imp_res', np.nan):.4f} | {o.get('gimp_binByBinStat', np.nan):.4f} | {o.get('gimp_stat', np.nan):.4f} | {o.get('gimp_theory', np.nan):.4f} |")
L += ["", f"| fit | mean σ(A_i) | mean σ(e_i) | mean σ(M_i) | mean ρ(A_i,e_i) | mean ρ({P},A_i) | mean ρ({P},e_i) | mean ρ({P},M_i) | edmval |", "|---|--:|--:|--:|--:|--:|--:|--:|--:|"]
for label, o in rows:
    L.append(f"| {label} | {o['sig_A']:.3f} | {o['sig_e']:.3f} | {o['sig_M']:.3f} | {o['rho_Ae']:+.3f} | {o['rho_mZ_A']:+.3f} | {o['rho_mZ_e']:+.3f} | {o['rho_mZ_M']:+.3f} | {o['edm']:.1e} |")
L += ["", "Posterior σ of the scale parameters are in units of their prior width. Impacts are rabbit's marginal group impacts, "
      "recomputed here as sqrt(σ²(mZ) − σ²(mZ | group fixed)); the conditional columns fix the group at its best fit with everything else floating."]
open(args.output, "w").write("\n".join(L) + "\n"); print("\n".join(L))

if args.plot and rows:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#4a3aa7"]
    # x = statistics point; one line per channel combination (label containing "Z+B" is the no-J/psi arm)
    order = ["res5", "data1", "data10", "data50"]; xlab = ["MC as produced\ndata x1", "MC 2x data\ndata x1", "MC 2x data\ndata x10", "MC 2x data\ndata x50"]
    def xpos(label):
        for i, k in reversed(list(enumerate(order))):   # longest keys first: data50, data10, then data1
            if k in label: return i
        return None
    arms = {"Z + J/ψ + B⁺": [(xpos(l), o) for l, o in rows if "Z+B" not in l and "only" not in l and xpos(l) is not None],
            "Z + B⁺": [(xpos(l), o) for l, o in rows if "Z+B" in l and xpos(l) is not None]}
    ref = [o for l, o in rows if "only" in l]
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.3))
    panels = [("sigma_mZ", "σ(mZ) [MeV]"), ("imp_A+e+M", "impact of A+e+M on mZ [MeV]"), ("sig_e", "mean σ(e_i) [prior units]")]
    for ax, (key, title) in zip(axs, panels):
        for (arm, pts), c, mk in zip(arms.items(), C, "os"):
            pts = sorted(pts, key=lambda t: t[0]); ax.plot([x for x, _ in pts], [o[key] for _, o in pts], marker=mk, color=c, lw=1.8, ms=7, label=arm)
        if ref: ax.axhline(ref[0][key], color="#52514e", ls="--", lw=1.0, label="Z + J/ψ (no B⁺)")
        if key == "sigma_mZ": ax.axhline(rows[0][1]["cond_A+e+M"], color="#52514e", ls=":", lw=1.0, label="A, e, M fixed (floor)")
        ax.set_xticks(range(4)); ax.set_xticklabels(xlab, fontsize=8.5); ax.set_title(title, fontsize=11); ax.grid(color="#e6e6e3"); ax.set_ylim(bottom=0)
    axs[0].legend(frameon=False, fontsize=9)
    fig.suptitle("B+ statistics scaling: what changes when the B+ simulation and data grow", fontsize=11.5)
    fig.tight_layout(); fig.savefig(args.plot, dpi=130); print("wrote", args.plot)
