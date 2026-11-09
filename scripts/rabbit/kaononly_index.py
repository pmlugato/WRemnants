#!/usr/bin/env python3
"""Index pages for the kaon-only vs all-legs study directory."""
import os, glob, html
P = os.path.expanduser("~/public_html/mz/btojpsik_kaononly_260921")
HEAD = ("<html><head><meta charset='utf-8'><title>{t}</title><style>body{{font-family:sans-serif;max-width:1500px;"
        "margin:auto}} img{{max-width:100%;border:1px solid #ddd}} figcaption{{color:#52514e;font-size:14px}} "
        "pre{{background:#f6f6f4;padding:12px;white-space:pre-wrap}}</style></head><body>")
SUB = {
 "01_bplus_alone": ("Standalone B+ Asimov", "The B+ channel fitted on its own, so the muon dilution shows up undiluted by the J/psi. Compare <b>alllegs</b> (A/e/M applied to bachelor and both muons, the configuration behind every earlier fit) with <b>kaononly</b> (bachelor only). If the dilution is real, rho(A,e) falls and sigma(e) improves."),
 "02_mz_joint": ("Z + J/psi + B+ Asimov", "The combined mZ fit. <b>alllegs</b> is res5_260914 at data x1 and the statscale runs at x10 and x50; <b>kaononly</b> is the same chain on the kaon-only histogram. Watch sigma(mZ) and the per-eta sigma(e) gain over Z + J/psi."),
 "03_mB_joint": ("J/psi + B+ with the B+ mass floating", "The fit we intend to turn into a measurement of mB. <b>alllegs</b> is realmc_C (Asimov) and realmc_A (blinded data); <b>kaononly</b> is the same on the kaon-only histogram. Watch sigma(mB), its scale term, and rho(mB, A/e)."),
 "04_recompute": ("Analytic recompute vs the reweighted templates", "The A/e/M shift applied directly to the reconstructed four-vectors and the mass rebuilt, with and without the J/psi mass constraint. This is what the reweighted templates should reproduce and the reference for deciding which leg configuration is right."),
}
for d, (title, blurb) in SUB.items():
    p = os.path.join(P, d)
    if not os.path.isdir(p): continue
    L = [HEAD.format(t=title), f"<h1>{title}</h1><p>{blurb}</p><p><a href='../'>back to the study index</a></p>"]
    for md in sorted(glob.glob(f"{p}/*.md")):
        L.append(f"<h2>{os.path.basename(md)}</h2><pre>{html.escape(open(md).read())}</pre>")
    for png in sorted(glob.glob(f"{p}/*.png")):
        L.append(f"<figure><img src='{os.path.basename(png)}'><figcaption>{os.path.basename(png)}</figcaption></figure>")
    open(f"{p}/index.html", "w").write("\n".join(L) + "</body></html>")
L = [HEAD.format(t="B+ scale variations: kaon only vs all legs"),
     "<h1>B+ A/e/M variations: kaon only vs all three legs (2026-09-21)</h1>",
     "<p>Every earlier B+ fit applied the A/e/M variations to the bachelor AND both muons "
     "(<code>scaleVarLegs = all</code>). The joint CVH refit constrains the dimuon mass to the J/psi, so a common "
     "muon scale error should largely cancel in m(&mu;&mu;K); but the variations are event reweights, not "
     "recomputations, so they never see the constraint. This study measures what that costs. Plan and running log: "
     "<code>notes/jpsi_bplus_mB_fit.md</code>, option 0.</p>",
     "<p>Inputs: histogram <code>btojpsik_cvhnano_kpiB_kaononly.hdf5</code> (identical to "
     "<code>kpiA_kpi</code> except <code>--scaleVarLegs bachelor</code>; closure verified, every nominal histogram "
     "bit-identical, only the variation content differs). Fit tags begin with <code>kaononly_</code> or "
     "<code>alllegs_</code>.</p>", "<ul>"]
for d, (title, blurb) in SUB.items():
    if os.path.isdir(os.path.join(P, d)):
        n = len(glob.glob(f"{P}/{d}/*.png"))
        L.append(f"<li><a href='{d}/'>{d}</a> &mdash; {title} ({n} figures)</li>")
L.append("</ul></body></html>")
open(f"{P}/index.html", "w").write("\n".join(L))
print("indexes written")
