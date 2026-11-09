#!/usr/bin/env python3
"""Deck: the B+ A/e/M unit bug, the fix, what changes, and the statistics behind it.

Source: notes/aem_unit_mismatch.md, section 9. Build:

    python scripts/btojpsik/make_aem_unit_fix_figs.py           # container: own figures
    bash   scripts/btojpsik/make_aem_unit_fix_rabbit_plots.sh   # container: rabbit plots
    python scripts/btojpsik/make_aem_unit_fix_deck.py           # writes public_html/slides/
"""

import os
import sys

sys.path.insert(0, os.path.expanduser("~/.claude/skills/slides"))
from deck import Deck  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ASSETS = os.path.join(HERE, "..", "..", "slides", "aem_unit_fix", "assets")
R = "rabbit/"
IMP_MB = R + "traditional_impacts_grouped_massShiftB5MeV_mB_data.png"
IMP_MB_ASIMOV = R + "traditional_impacts_grouped_massShiftB5MeV_mB_asimov.png"
IMP_MZ = R + "traditional_impacts_grouped_massShiftZ100MeV_mZ_asimov.png"
PULLS_E = R + "pulls_and_traditional_impacts_massShiftB5MeV_e_pulls_data.png"
PULLS_E_BARREL = R + "pulls_and_traditional_impacts_massShiftB5MeV_e_pulls_barrel_data.png"
PULLS_A = R + "pulls_and_traditional_impacts_massShiftB5MeV_A_pulls_data.png"
CORR_OLD = R + "hist_corr_postfit_parms_x_parms_y_params_corr4x4_old_eta12.png"
CORR_NEW = R + "hist_corr_postfit_parms_x_parms_y_params_corr4x4_matched_eta12.png"


def build():
    d = Deck(assets=ASSETS)
    d.title_slide(
        "The B⁺ momentum-scale units bug, explained",
        "Why A/e/M did not move mB, what the fix changes, and the statistics behind it",
        "Pietro Lugato  ·  8 October 2026  ·  mZ / B⁺ → J/ψ K⁺")

    d.bullets("Outline", [
        "Physics: why mB has to depend on the kaon momentum scale",
        "Primer: a nuisance parameter is defined by its template",
        "The bug, what the B⁺ templates looked like, and the fix",
        "Primer: uncertainties as curvature; profiling, correlations, impacts",
        "Results on Asimov: mB and mZ, before and after",
        "Primer: Asimov, data, blinding, weighted averages, goodness of fit",
        "Data: the scale values move, and the B⁺ and J/ψ disagree",
    ])

    # ------------------------------------------------------------ motivation
    d.bullets("Last time: the scale looked irrelevant to mB", [
        "J/ψ + B⁺ data fit, mB floating: σ(mB) = 0.0715 MeV, of which A+e+M only 0.0014 MeV",
        "Correlation of mB with every A_i, e_i: |ρ| ≈ 0.003",
        "Adding the B⁺ to Z + J/ψ improved σ(mZ) by only 0.2 %",
        "An absolute mass insensitive to the momentum scale is not physical",
        "→ something in the model was switching the kaon's lever arm off",
    ])

    d.bullets_figure("First principles: the kaon momentum sets mB", [
        "Muons drop out: the CVH refit fixes m(μμ) = m(J/ψ) (checked to 4 digits)",
        "Full recompute, 236k truth-matched B⁺: ∂mB/∂A = −1200 to −1570 MeV",
    ], "eq_kaon_scale.png")

    # ------------------------------------------------------- primer: nuisances
    d.bullets_figure("Primer: a nuisance is defined by its template", [
        "θ = 1 means “the up template”; the Gaussian prior says θ is of order 1",
        "The physical size of θ = 1 is not stored: it is whatever the histmaker shifted",
        "merge_inputs correlates by NAME: A17 must be the same shift in Z, J/ψ and B⁺",
    ], "eq_nuisance.png")

    # ----------------------------------------------------------------- the bug
    d.bullets_figure("The bug: “A17” meant a different shift per channel", [
        "Z, J/ψ histmakers: flat fitted-scale widths (--fitMuonScaleAndResolution)",
        "B⁺ histmaker: flag parsed but never passed on → massfit widths, 20–60× smaller",
        "η median: one B⁺ unit of A / e / M was 34 / 41 / 20× too small (7 / 8 / 20× vs repro Z)",
    ], "units_per_eta.png")

    d.figure("What the B⁺ template looked like, before and after", "templates_e10_old_vs_matched.png",
             lead="Signal, kaon pT < 2.8 GeV, −0.4 < η < −0.2, at θ(e10) = ±0.2, a typical postfit value")

    d.figure("The response per unit was 30× too small", "template_slope_per_unit.png",
             lead="B⁺ mass shift per unit of the shared nuisance, at small θ: 0.2–1 MeV before, 5–30 MeV after")

    # ----------------------------------------------------------------- the fix
    d.bullets("The fix: one flag, a guard, and a linear build", [
        "btojpsik_cvhnano.py now passes fit_muon_scale=--fitMuonScaleAndResolution",
        "Guard: check_joint_tensor.py computes one unit of A/e/M per η for every channel",
        ("joint fit driver refuses to fit on a mismatch; every past B⁺ merge would now be refused", 1),
        "Closure on 95.5M events: every nominal histogram identical (to 4×10⁻¹⁵)",
        "Variations scale as the width ratio: A to 1 %, M to 1 % (median); e has a catch (next)",
    ])

    d.bullets_figure("Primer: a reweight cannot move events it does not have", [
        "Template = MC reweighted by w = P(reco | truth shifted) / P(reco | truth)",
        "Fine while the shift ≪ resolution; at 5/5/1 one unit of e is a 5 % shift on a 1 GeV kaon",
        "Built directly at 5/5/1, e responds at 0.67 of linear at 1–2 GeV (73 % of the signal)",
    ], "reweight_linearity.png")

    d.bullets_figure("A stretch in ln κ is exactly a change of unit", [
        "Build the B⁺ at 1/1/1, where the reweight is linear; stretch A and e in ln κ",
        "Measured small-θ slopes: stretched / 1/1/1 = 4.997–5.000 in every kaon-pT bin",
        "A linear stretch ν⁰ + 5(up − ν⁰) gave 57 negative bins and was refused",
    ], "eq_stretch.png")

    # ------------------------------------------------- primer: information
    d.bullets_figure("Primer: an uncertainty is a curvature", [
        "σ says how fast the fit gets worse as θ moves away from its best value",
        "Independent channels add information, i.e. curvature, not σ",
        "Information goes as the template slope squared",
    ], "eq_curvature.png")

    d.figure("Information adds: J/ψ, B⁺ and combined, for e", "information_parabolas_e.png",
             lead="Asimov, mean over η; B⁺ term = 1/σ²(combined) − 1/σ²(J/ψ alone), one-parameter view")

    d.bullets_figure("Primer: correlations, profiling and impacts", [
        "Profiled σ(POI): every other parameter floats; correlations with the POI inflate it",
        "Conditional σ(POI | G fixed): the group G frozen at its best fit",
        "Impact of G: what letting G float costs the POI",
    ], "eq_profiling.png")

    d.bullets_figure("Two parameters, one peak: what each B⁺ cell measures", [
        "Old unit: s_c ≈ 0, so every cell measures mB alone and says nothing about θ",
        "Matched: cell differences (pT: A vs e; charge: M; η) pin θ, so σ(mZ) falls",
        "The cell average is indistinguishable from mB, so σ(mB) rises",
    ], "eq_two_params.png")

    d.figure("Toy: why the B⁺ tightens θ while σ(mB) grows", "degeneracy_toy.png",
             lead="3 kaon-pT cells, measured e slopes, one global e, J/ψ ±0.068; shape, not values")

    # ---------------------------------------------------------- Asimov results
    d.figure("mB: the scale is now a real uncertainty", IMP_MB_ASIMOV,
             lead="J/ψ + B⁺ Asimov (rabbit): σ(mB) 0.072 → 0.089 MeV; A+e+M 0.0014 → 0.051 MeV")

    d.figure("…and the B⁺ now calibrates the scale itself", "aem_sigma_mB_chain.png",
             lead="J/ψ + B⁺ DATA: the B⁺ halves the J/ψ's σ(A_i), σ(e_i); data and Asimov agree to 5 %")

    d.figure_row("mB now correlates with the scale it depends on", [CORR_OLD, CORR_NEW],
                 lead="DATA, η group 12 (rabbit), old (left) vs matched (right): ρ(mB, e12) −0.005 → −0.23")

    d.figure("mZ: the B⁺ now buys 9.5 %", IMP_MZ,
             lead="Z + J/ψ + B⁺ Asimov, MC as produced (rabbit): σ(mZ) 2.20 → 2.00 MeV, A 1.70 → 1.42")

    # --------------------------------------------------- primer: data fits
    d.bullets("Primer: Asimov, data and blinding", [
        "Asimov: the “data” are the expected yields at θ = 0, so θ̂ = 0 exactly",
        "Asimov σ is the expected uncertainty: every Asimov change so far is in σ, never in θ̂",
        "Data: θ̂ moves to wherever the observed counts prefer",
        "Blinding: rabbit adds a hidden offset to mB only; all other parameters are released",
        ("verified in the fit log: massShiftB5MeV is the only blinded parameter", 1),
    ])

    d.bullets_figure("Primer: a postfit value is a weighted average", [
        "Old unit: B⁺ weight ≈ 3 % of the J/ψ's, so θ̂ ≈ θ_J (moved by at most 0.24 units)",
        "Matched: B⁺ weight 2.5–3.3× the J/ψ's, so θ̂ follows the B⁺",
        "One-parameter guide only: the real fit has 72 correlated scale parameters",
    ], "eq_weighted_average.png")

    d.figure("On data the scale values move to the B⁺'s preference", "postfit_Ae_data.png",
             lead="Barrel A, e move +0.3 to +0.55 units (7–12 σ of the old fit); M about −0.4 everywhere")

    d.figure("Seen in rabbit: e_i pulls tighten, and now carry mB impact", PULLS_E_BARREL,
             lead="DATA, barrel e6–e17 (rabbit): bars = impact on mB; markers = matched θ ± σ, boxes = old")

    d.bullets_figure("Primer: goodness of fit, and the tension it shows", [
        "J/ψ alone: 51,934 / 5,758 = 9.0 per dof; its model does not describe its data",
        "Joint J/ψ + B⁺ (ndof 8,375): 61,987 old unit → 64,266 matched, for the same data",
        "At 9.0 per dof the J/ψ's σ(θ) understates it: which side is right is open",
    ], "eq_gof.png")

    d.bullets("Summary and next steps", [
        "The scale was never negligible for mB: the B⁺ was told its kaons barely react to it",
        "Fixed and guarded; B⁺ built at 1/1/1, log-stretched × 5 for the mB chain",
        "mB (data, blinded): σ 0.0715 → 0.0839 MeV, scale ×34; old σ(mB) was an underestimate",
        "mZ (Asimov): σ 2.208 → 1.998 MeV (−9.5 %); the B⁺ halves σ(A), σ(e) in J/ψ + B⁺",
        "Next: J/ψ vs B⁺ compatibility test, then a common-mode scale nuisance",
    ])

    # ----------------------------------------------------------------- backup
    d.section("Backup")

    d.table("The guard: every past B⁺ merge would now be refused", [
        ["merge tested", "expected", "result"],
        ["J/ψ full + B⁺ old (all mB fits)", "FAIL", "FAIL, 34 / 41 / 20×"],
        ["Z + J/ψ repro + B⁺ old (mZ fits)", "FAIL", "FAIL, 6.9 / 8.1 / 20×"],
        ["Z repro + J/ψ repro", "PASS", "PASS"],
        ["Z repro + J/ψ full", "FAIL", "FAIL, 0.2× in A, e"],
        ["J/ψ full + B⁺ 1/1/1, unstretched", "FAIL", "FAIL"],
        ["J/ψ full + B⁺ 1/1/1, stretched × 5", "PASS", "PASS"],
    ], lead="check_joint_tensor.py, stage 3 and stage-4 gate of joint_with_btojpsik.sh",
       widths=[0.5, 0.18, 0.32],
       note="In this container sys.exit(1) after opening a ROOT file returns 0; the guard exits via os._exit.")

    d.table("σ(mB) decomposition, before and after [MeV]", [
        ["fit", "σ(mB)", "sim. stat", "data stat", "A+e+M"],
        ["old unit, Asimov", "0.0724", "0.0692", "0.0212", "0.0014"],
        ["matched, Asimov", "0.0887", "0.0692", "0.0212", "0.0513"],
        ["old unit, DATA", "0.0715", "0.0685", "0.0204", "0.0014"],
        ["matched, DATA", "0.0839", "0.0665", "0.0200", "0.0472"],
    ], lead="J/ψ (full tensor, μ model) + B⁺, MC as produced, mB blinded on data",
       widths=[0.32, 0.17, 0.17, 0.17, 0.17],
       note="The note's projection (0.1–0.3 MeV) assumed the B⁺ would not tighten A, e; it halves them.")

    d.figure("mB impacts on DATA (rabbit), matched vs old unit", IMP_MB,
             lead="J/ψ + B⁺ DATA, mB blinded: σ(mB) 0.0715 → 0.0839 MeV; A+e+M 0.0014 → 0.047 MeV")

    d.figure("All 24 e_i pulls and impacts on DATA (rabbit)", PULLS_E,
             lead="As the barrel slide, all η groups: boxes = old unit, markers = matched; bars = impact on mB")

    d.figure("A_i pulls and impacts on DATA (rabbit)", PULLS_A,
             lead="Same as the e_i plot: boxes = old unit, markers = matched; bars = impact on mB")

    d.figure("σ(mZ) and the scale impact, as bars", "mZ_before_after.png",
             lead="Z + J/ψ + B⁺ Asimov; black tick = σ(mZ) with A, e, M fixed (1.40 MeV in all three)")
    return d


if __name__ == "__main__":
    print(build().save("aem_unit_fix.pptx"))
