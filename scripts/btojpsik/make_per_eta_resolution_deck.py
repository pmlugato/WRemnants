#!/usr/bin/env python3
"""Deck for openspec/changes/add-cvhnano-per-eta-resolution.

Context first, then one section per deliverable. Fit plots are rabbit's own
(prefit/postfit, pulls and impacts, correlations); the bespoke figures are the
resolution-against-eta plots from cvhnano_resolution_report.py.

    python scripts/btojpsik/make_per_eta_resolution_deck.py
    python ~/.claude/skills/slides/check.py ~/public_html/slides/<date>_btojpsik_per_eta_resolution.pptx
"""

import os
import sys

sys.path.insert(0, os.path.expanduser("~/.claude/skills/slides"))
from deck import Deck  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ASSETS = os.path.join(HERE, "..", "..", "slides", "per_eta_resolution", "assets")


def has(png):
    return os.path.exists(os.path.join(ASSETS, png))


def build():
    d = Deck(assets=ASSETS)

    d.title_slide(
        "B⁺ → J/ψ K⁺: the lineshape enters the fit",
        "Per-eta resolution nuisances and the peaking template, at 3 pT × 24 eta × 2 charge",
        "Pietro Lugato  ·  14 September 2026  ·  mZ / kaon calibration")

    # ---------------------------------------------------------------- context
    d.bullets("Context: what the B⁺ channel is for", [
        "One hadron track at 1–10 GeV, reconstructed by the CVH joint refit, in a regime no dimuon channel reaches",
        "Its job: constrain the kaon momentum scale A/e/M, 24 eta groups, shared by name with Z and J/ψ",
        "Joint Asimov (2026-09-01): B⁺ moves σ(mZ) by 0.3 %; it is a cross-check, and that needs the data fit",
        "No B⁺ data fit at 144 cells has converged: the model was wrong at the 1 % level while data is at 0.36 %",
    ])
    d.bullets("Context: what the 2026-09-08 fit had, and lacked", [
        "Had: signal MC template, exp·poly2 background per eta, 72 A/e/M nuisances, 10 × 15 MeV bins, 144 cells",
        "Lacked a resolution degree of freedom: worth χ² 6151 → 3312 on the 1-cell data spectrum",
        "Lacked the peaking-b shapes: the B⁰ edge in the low shoulder carried 52 % of the reference χ²",
        "Asimov signal 5–20 % below the data-fitted yield, because the seed cap fired on every tensor",
        ("A single global smear was known to be wrong: the deficit runs 9 MeV (soft kaons) to 4 MeV (hard)", 1),
    ])
    d.bullets("What this study adds", [
        "24 free resolution nuisances, one per kaon eta bin, shared over pT and charge (the physics choice)",
        "The merged peaking-b template (B⁺→J/ψ K X, B⁺→J/ψ π X, B⁰, Bs, other b) with a free norm per eta",
        "Asimov signal pinned to the data-fitted 1 327 932 events: scale 18.24, cap bypassed",
        "20 × 20 MeV bins over 5.05–5.45 GeV; kpi response model, un-smeared muon model as cross-check",
        "Eight Asimov arms, then the observation; every fit gated on edmval, never on the output file",
    ])
    d.table("Why eta, and what each arm asks", [
        ["arm", "resolution", "template", "question"],
        ["primary", "per eta, free", "peak", "the result"],
        ["constrained", "per eta, unit prior", "peak", "does the prior move A/e/M"],
        ["eta × pT", "72, free", "peak", "what eta-only leaves out"],
        ["global / none", "1 / 0", "peak", "cost of the resolution freedom"],
        ["no template", "per eta", "none", "cost of the peaking template"],
        ["split", "per eta", "comb + peak", "the known comb degeneracy"],
        ["3 MeV ref.", "per eta", "peak", "reference-width dependence"],
    ], lead="Multiple scattering follows the material, a function of eta; the pT dependence is kinematic.",
        widths=[0.18, 0.24, 0.18, 0.40])

    # ---------------------------------------------------------------- asimov
    d.table("Asimov: lineshape freedom barely costs the scale", [
        ["arm", "mean σ(A)", "mean σ(e)", "mean σ(M)", "ρ(A,e)"],
        ["primary (per-eta res + peak)", "0.9544", "0.6777", "0.9942", "0.227"],
        ["no resolution", "0.9542", "0.6764", "0.9941", "0.228"],
        ["eta × pT resolution", "0.9548", "0.6800", "0.9942", "0.226"],
        ["no template", "0.9520", "0.6725", "0.9936", "0.233"],
        ["split (comb + peak)", "0.9608", "0.8141", "0.9948", "0.152"],
        ["3 MeV reference", "0.9544", "0.6779", "0.9942", "0.227"],
        ["bare (2026-09-08 model)", "0.9518", "0.6714", "0.9936", "0.234"],
    ], lead="kpi histogram, 144 cells, 20 × 20 MeV, pinned scale; 8 of 8 arms converged. Prior units.",
        widths=[0.40, 0.15, 0.15, 0.15, 0.15],
        note="24 free widths cost 0.2 % on σ(e); the template 0.8 %; the combinatorial template 20 %.")
    d.figure("Asimov: σ(A/e/M) per eta group, arm by arm", "asimov_aem_per_group.png",
             lead="Every arm but split lies on the primary to < 1 %; the ratio panels are the whole story.")
    d.figure("Asimov: how well each eta bin determines its own width", "asimov_res_vs_eta.png",
             lead="σ(θ) is 0.15 at |η| < 0.3 and 0.4–0.5 at |η| > 1.9; a unit prior only tightens the edges.")
    d.bullets_figure("Asimov: eta × pT widths are unconstrained at the edges", [
        "72 nuisances: σ(θ) 0.62 on average, 1.35 in the soft-pT edge cells",
        "Costs 0.3 % on σ(e) against eta-only; largest res–scale |ρ| 0.12 → 0.21",
    ], "asimov_res_vs_eta_pt.png")
    d.figure("Asimov: resolution and scale do not talk", "asimov_corr_scale_vs_res.png",
             lead="Mid-eta A/e/M against all 24 widths: largest |ρ| anywhere is 0.12.")

    return d


def add_data(d):
    """Data-side slides."""
    d.bullets("Data: the fits stalled, and not because of the model", [
        "Every arm: loss 986k → 204k in 350 iterations, then flat for 67 000; Hessian not positive-definite",
        "Same stall without the resolution nuisances, without the template, with a sideband-exponential seed",
        "Cause: rabbit blinds every normalisation on an observed fit by exp(N(0, 5)) unless --unblind is passed",
        "Unblinded prefit loss 46 318 against 986 182 blinded; every B⁺ data fit on record ran blinded",
        ("unblinded prefit Hessian: 3 negative eigenvalues of 336, each one background amplitude in a sparse edge cell", 1),
    ])
    d.bullets("Data: unblinded, the primary fit converges", [
        "edmval 1.2e-12 after 20 iterations; prefit loss 8 048 → 1 782",
        "2ΔNLL to the saturated model 3 564 over 2 616 degrees of freedom: 1.36 per d.o.f. (1-cell models: 53.6, 10.1)",
        "First converged B⁺ data fit at 3 pT × 24 eta × 2 charge on record",
        "Signal normalisations 0.71–1.11, eta-shaped; peaking-template normalisations 0.69–1.31",
        "σ(e) on data 0.679 against 0.678 on Asimov; ρ(A,e) 0.225",
        ("un-smeared muon model instead of kpi: same θ to 0.04, e pull −0.94 vs −0.90, A +0.80 vs +0.73", 1),
    ])
    d.figure("Data: resolution per eta bin is 1.3 to 2.5, edges higher", "data_res_vs_eta_arms.png",
             lead="Mean θ 1.66 at a 5 MeV reference ≈ 6.4 MeV added width; the 1-cell study found 6.26 MeV.")
    d.figure("Data: A pulls positive, e pulls negative, M at zero", "data_aem_pulls.png",
             lead="Mean A +0.7, mean e −0.9, M 0.0 ± 0.1 prior units: a coherent data/simulation scale offset.")
    d.table("Data: what each ingredient is worth, and what it moves", [
        ["arm", "NLL (reduced)", "2ΔNLL vs primary", "mean e pull", "max |Δe|/σ"],
        ["primary: per-eta res + peak", "1782", "—", "−0.90", "—"],
        ["constrained resolution", "1813", "+63", "−0.92", "0.12"],
        ["eta × pT resolution", "1746", "−71 (48 params)", "−0.89", "0.29"],
        ["no resolution", "2223", "+882", "−1.03", "0.87"],
        ["no template", "24 978", "+46 400", "−0.12", "2.15"],
        ["bare (2026-09-08 model)", "25 512", "+47 500", "−0.23", "2.34"],
    ], lead="All arms converged unblinded; 2 880 bins. Shifts in units of the primary's per-group σ.",
        widths=[0.36, 0.16, 0.20, 0.14, 0.14],
        note="The peaking template carries the likelihood and moves e by 0.8 prior units; the resolution nuisances move it by at most 0.9σ and A by 0.4σ.")
    d.figure("Data: eta × pT widths add little", "data_res_vs_eta_pt.png",
             lead="Soft-pT bin θ ≈ 2, hard bins ≈ 1 (the 9 → 4 MeV trend); 2ΔNLL gain 71 for 48 params.")
    d.figure_row("Data, all 144 cells summed: prefit and postfit",
                 ["data_prefit_inclusive.png", "data_postfit_inclusive.png"],
                 lead="Left: flat seed, narrow template. Right: within 1 % in every bin.")
    d.figure_row("Soft kaons, |η| = 2.3, q = −: prefit and postfit",
                 ["data_prefit_cell_soft_edge.png", "data_postfit_cell_soft_edge.png"],
                 lead="Purity 5 %; flat seed, narrow template. Cell χ²/20 = 0.51 (data variance only).")
    d.figure_row("Hard kaons, |η| = 2.3, q = +: prefit and postfit",
                 ["data_prefit_cell_hard_edge.png", "data_postfit_cell_hard_edge.png"],
                 lead="Purity 60 %; the steps are single simulated events × 18. Cell χ²/20 = 0.06.")
    d.figure_row("Soft kaons, |η| < 0.1, q = −: prefit and postfit",
                 ["data_prefit_cell_soft_central.png", "data_postfit_cell_soft_central.png"],
                 lead="Highest statistics: width deficit visible prefit, absorbed postfit. Cell χ²/20 = 1.31.")
    d.figure_row("Hard kaons, |η| < 0.1, q = −: prefit and postfit",
                 ["data_prefit_cell_hard_central.png", "data_postfit_cell_hard_central.png"],
                 lead="Purity 60 %; the peaking template carries the shoulder at the seed. Cell χ²/20 = 0.29.")
    d.figure("What the fit sees: the e variation at η ≈ 0", "var_scalee12.png",
             lead="±1σ of e moves the mean by 0.34 MeV at soft pT and 0.08 at hard pT: the k lever arm.")
    d.figure("The A variation is flat in pT, as a scale must be", "var_scaleA12.png",
             lead="Same cells: Δ⟨m⟩ −0.07 MeV soft, −0.09 hard, opposite in sign to e. Hatched: < 3 MC events.")
    d.figure("The e variation at |η| = 2.3 is ten times larger", "var_scalee0.png",
             lead="+3.2 MeV soft, +0.8 MeV hard for the same prior width: more material, more energy loss.")
    d.figure("The per-eta width variation, η ≈ 0", "var_sigRes_eta12.png",
             lead="A pure shape change: the peak bin loses, the neighbours gain; the mean moves ≤ 0.4 MeV.")
    d.figure("Data: every cell's own χ²/ndf, postfit", "data_cell_chi2.png",
             lead="Data variance only; median 0.25; 5 of 144 cells above 2, all soft kaons at |η| < 0.5.")
    d.figure_row("Data: A/e/M correlations, barrel and endcap",
                 ["data_corr_barrel.png", "data_corr_endcap.png"],
                 lead="Postfit: ρ(A,e) 0.25 barrel, 0.23 endcap; M decoupled from both at the 10⁻³ level.")
    d.figure_row("Data: leading impacts on the signal yield",
                 ["data_impacts_barrel.png", "data_impacts_endcap.png"],
                 lead="Top 12 on the η ≈ 0 (left) and η ≈ −2.3 (right) yields; e, curvature, width and A lead.")
    return d


def add_backup(d):
    d.section("Backup")
    d.bullets("Two kernel flaws found and fixed on the way", [
        "Gaussian sampled at 20 MeV bin centres with a 5 MeV width is exp(−8) off-diagonal: a delta function",
        ("first per-eta tensor moved the peak bin by 0.02 %", 1),
        "Exact bin-transfer matrix on 20 MeV bins assumes a flat density inside a bin holding a 14 MeV core",
        ("over-smeared by ~50 %; σ(θ) scaled as σ_ref, not σ_ref²; peak shift 15.6 % → 10.5 % once fixed", 1),
        "Now: transfer matrix on the histmaker's 5 MeV bins, rebinned by 4; rebinned un-smeared = nominal asserted",
    ])
    d.bullets("Pinned Asimov scale, and what the cap had been doing", [
        "signalScale = 1 327 932 / S_MC(5.05–5.45) = 1 327 932 / 72 790 = 18.2433 on both kpiA histograms",
        "Window-independent: a global multiplier, so the 10-bin fallback would take the same number",
        "Default seed on this histogram: 19.35 capped to 17.81 (−2.4 %); the 2026-09-08 tensor sat at 14.64 (−19.7 %)",
        "Asimov signal integral 1 327 929.8 against 1 327 932: ratio 0.9999983",
    ])
    d.bullets("The --mcBkg options, precisely", [
        "none: signal + analytic exp·poly2 only (2026-09-08 configuration)",
        "peak (new): categories 1–5, B⁺→J/ψ K X, B⁺→J/ψ π X, B⁰, Bs, other b — 15 393 raw MC events",
        "comb: category 6, no common b ancestor — 192 686 raw; 3.3× noisier than the data it would describe",
        "split: comb + peak as two templates; all: categories 1–6 as one",
        "'Peaking' is a label, not literal: peak fraction 0.04 vs 0.83 for signal; it carries the B⁰ edge",
    ])
    d.table("Why ρ(A,e) reads 0.25 here and 0.6 on 2026-08-25", [
        ["configuration", "σ(A)", "σ(e)", "ρ posterior", "ρ data-only"],
        ["Z + J/ψ joint Asimov", "0.12", "0.11", "0.85", "0.86"],
        ["B⁺, 4 η groups × 6 pT, 2026-08-25", "0.80", "0.34", "0.59", "0.80"],
        ["B⁺, 24 η groups × 3 pT, Asimov (this)", "0.95", "0.68", "0.23", "0.81"],
        ["B⁺, 24 η groups × 3 pT, data (this)", "0.95", "0.68", "0.23", "0.82"],
    ], lead="Posterior precision P = I + F; the diagonal unit prior wins where F is small.",
        widths=[0.44, 0.12, 0.12, 0.16, 0.16],
        note="Strip the prior, F = V⁻¹ − I on the (A, e) block, and every configuration gives the same physics: ρ ≈ 0.8–0.9, the k = 1/pT lever arm.")
    d.bullets("The same prior dilution, one step further", [
        "A and e both shift the peak; within an η group only the spread of k across the pT cells separates them, so the data's own ρ is ~0.85",
        "24 groups instead of 4 put six times fewer events behind each nuisance: F drops by ~6, the posterior slides toward the prior",
        "Scaling the 2026-08-25 information by 1/6–1/8 predicts σ(A) 0.94–0.95, σ(e) 0.60–0.65, ρ 0.24–0.28: observed 0.955, 0.679, 0.225",
        "So 0.10 (200 files) → 0.60 (full, 4 groups) → 0.25 (full, 24 groups) is one curve, F/(I+F), not three physics results",
        ("Compare configurations on the prior-stripped ρ; the report now prints it.", 1),
    ])
    d.figure_row("The same variations as rabbit draws them",
                 ["rabbit_var_pt0_eta12.png", "rabbit_var_pt2_eta12.png"],
                 lead="rabbit_plot_inputdata: e and width on the signal, soft (left) and hard (right), η ≈ 0.")
    d.bullets("Reference width: 5 MeV, checked at 3", [
        "A/e/M and ρ(A,e) identical to 0.05 % between the two references",
        "σ(θ) scales as 0.540 / 0.341 = 1.58, near σ_ref not σ_ref²: on single-event templates the smear is first order",
        "So θ reads as a width in units of the reference; 'implied added width' in the report is first-order only",
        "θ range check [−1, 2] is a data-side test; on Asimov every θ is 0 by construction",
    ])
    return d


if __name__ == "__main__":
    deck = build()
    deck = add_data(deck)
    deck = add_backup(deck)
    print(deck.save("btojpsik_per_eta_resolution.pptx"))
