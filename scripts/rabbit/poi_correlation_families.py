#!/usr/bin/env python3
"""Correlation of one POI with every parameter of the B+ nuisance families, per
family, as bar charts ordered by eta (and pT/charge where the family has them).
Families: resolution widths, signal norms, peaking-template norms, background
lnAmpl (per cell), background poly1 / poly2 (per eta), A/e/M. Run in the container.
"""
import argparse, re, sys
import numpy as np
sys.path.insert(0, "rabbit"); sys.path.insert(0, "wums")
from rabbit import io_tools
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("fit"); ap.add_argument("--poi", default="massShiftB5MeV"); ap.add_argument("--poiLabel", default="mB")
ap.add_argument("-o", "--output", required=True); ap.add_argument("--table", default=None)
a = ap.parse_args()
r = io_tools.get_fitresult(a.fit); names = np.array([str(x) for x in r["parms"].get().axes["parms"]])
cov = np.asarray(r["cov"].get().values()); sig = np.sqrt(np.diag(cov)); rho = cov / np.outer(sig, sig)
iz = int(np.flatnonzero(names == a.poi)[0])
fam = {
    "resolution widths sigRes_eta{i}": r"^sigRes_eta(\d+)$",
    "signal norm per eta": r"^norm_BuJpsiK_fit_bachelor_eta(\d+)$",
    "peaking-template norm per eta": r"^norm_BuJpsiKMCPeak_fit_bachelor_eta(\d+)$",
    "background poly1 per eta": r"^poly1_BuJpsiKBkg_fit_bachelor_eta(\d+)$",
    "background poly2 per eta": r"^poly2_BuJpsiKBkg_fit_bachelor_eta(\d+)$",
    "A per eta": r"^A(\d+)$", "e per eta": r"^e(\d+)$", "M per eta": r"^M(\d+)$",
}
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#4a3aa7", "#2a78d6", "#eb6834"]
fig, axs = plt.subplots(3, 3, figsize=(15, 10)); axs = axs.ravel()
lines = [f"| family | n | mean ρ({a.poiLabel},x) | max |ρ| | at |", "|---|--:|--:|--:|---|"]
for k, (lab, pat) in enumerate(fam.items()):
    sel = [(int(re.match(pat, n).group(1)), i) for i, n in enumerate(names) if re.match(pat, n)]
    sel.sort(); x = [s[0] for s in sel]; y = [rho[iz, s[1]] for s in sel]
    ax = axs[k]; ax.bar(x, y, color=C[k], width=0.8); ax.axhline(0, color="#52514e", lw=0.8)
    ax.set_title(f"ρ({a.poiLabel}, {lab})", fontsize=10); ax.set_xlabel("eta group"); ax.grid(color="#e6e6e3", axis="y")
    if sel:
        j = int(np.argmax(np.abs(y))); lines.append(f"| {lab} | {len(sel)} | {np.mean(y):+.4f} | {abs(y[j]):.4f} | {names[sel[j][1]]} |")
# background amplitudes: 144 cells, plotted against cell index grouped by (pt, charge)
ax = axs[8]; pat = r"^lnAmpl_BuJpsiKBkg_fit_bachelor_pt(\d+)_fit_bachelor_eta(\d+)_fit_bachelor_charge(\d+)$"
sel = [(int(m.group(1)), int(m.group(3)), int(m.group(2)), i) for i, n in enumerate(names) for m in [re.match(pat, n)] if m]
sel.sort(); y = [rho[iz, s[3]] for s in sel]
cols = {(0, 0): "#2a78d6", (0, 1): "#7fb0e8", (1, 0): "#eb6834", (1, 1): "#f3a98a", (2, 0): "#1baf7a", (2, 1): "#8fd7bc"}
ax.bar(range(len(sel)), y, color=[cols[(s[0], s[1])] for s in sel], width=0.9); ax.axhline(0, color="#52514e", lw=0.8)
ax.set_title(f"ρ({a.poiLabel}, background lnAmpl per cell)  colour: pT bin (blue/orange/green), light = charge +", fontsize=9); ax.set_xlabel("cell (pt, charge, eta)"); ax.grid(color="#e6e6e3", axis="y")
if sel:
    j = int(np.argmax(np.abs(y))); lines.append(f"| background lnAmpl per cell | {len(sel)} | {np.mean(y):+.4f} | {abs(y[j]):.4f} | {names[sel[j][3]]} |")
ymax = max(abs(v) for ax_ in axs for v in (ax_.get_ylim())) if False else None
for ax_ in axs: lo, hi = ax_.get_ylim(); m = max(abs(lo), abs(hi), 0.02); ax_.set_ylim(-m, m)
fig.suptitle(f"Correlation of {a.poiLabel} with the B+ nuisance families  ({a.fit.split('/')[-1]})", fontsize=11)
fig.tight_layout(); fig.savefig(a.output, dpi=130); print("wrote", a.output); print("\n".join(lines))
if a.table: open(a.table, "w").write("\n".join(lines) + "\n")
