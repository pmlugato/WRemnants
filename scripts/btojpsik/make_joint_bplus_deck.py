#!/usr/bin/env python3
"""Deck for worklist item 5: Z + J/psi against Z + B+ and Z + J/psi + B+, on
the shared 24-eta A/e/M basis, with the 2026-09-14 B+ model.

    python scripts/btojpsik/make_joint_bplus_deck.py
"""

import os
import sys

sys.path.insert(0, os.path.expanduser("~/.claude/skills/slides"))
from deck import Deck  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ASSETS = os.path.join(HERE, "..", "..", "slides", "joint_mz_jpsi_bplus", "assets")


def build():
    d = Deck(assets=ASSETS)
    d.title_slide(
        "mZ + J/ψ, mZ + B⁺, mZ + J/ψ + B⁺",
        "Combined Asimov fits on one A/e/M basis, with the converged B⁺ lineshape model",
        "Pietro Lugato  ·  14 September 2026  ·  mZ / kaon calibration")

    d.bullets("Context: one basis, three channels", [
        "72 scale nuisances A/e/M in 24 eta bins, shared by NAME across channels through merge_inputs",
        "Z: mll templates with the 144-parameter scale variations; J/ψ: prompt and displaced dimuon mass",
        "B⁺: kaon momentum scale through m(μμK), 3 pT × 24 eta × 2 charge cells, one hadron track at 1–10 GeV",
        "August result on the old B⁺ template: B⁺ improves σ(mZ) by 0.3 %, ρ(A,e) 0.854 → 0.852",
    ])
    d.bullets("Context: what changed in the B⁺ arm since August", [
        "Peaking-b template (B⁺→J/ψ K X, B⁺→J/ψ π X, B⁰, Bs) with a free norm per eta",
        "24 free per-eta resolution nuisances at a 5 MeV reference, smeared at the stored 5 MeV granularity",
        "Signal scale pinned to the data-fitted 1.33 M events; 20 × 20 MeV bins; kpi response model",
        "This model converges on data (2ΔNLL/ndf 1.36) once rabbit's blinding is switched off; Asimov is unblinded anyway",
        ("Z and J/ψ tensors are the validated August Asimov tensors, refitted here with identical rabbit options", 1),
    ])
    d.bullets("What to look at", [
        "σ(A), σ(e), σ(M) per eta bin and their means, in prior units, across the three configurations",
        "ρ(A_i, e_i): posterior, and with the unit prior stripped (F = V⁻¹ − I), which is the physics number",
        "σ(mZ) and the impacts on it, grouped and by nuisance; prefit and postfit of every channel",
        "All plots from rabbit's own tools; the cross-configuration figures from the study's report script",
    ])
    return d


def add_results(d):
    d.table("Three converged Asimov fits on the shared 24-eta basis", [
        ["configuration", "σ(mZ) [MeV]", "mean σ(A)", "mean σ(e)", "mean σ(M)", "ρ(A,e)", "ρ(A,e) data-only"],
        ["Z + J/ψ (nominal)", "2.208", "0.120", "0.106", "0.069", "0.854", "0.858"],
        ["Z + B⁺", "5.079", "0.213", "0.653", "0.110", "0.767", "0.852"],
        ["Z + J/ψ + B⁺", "2.204", "0.120", "0.105", "0.069", "0.853", "0.857"],
    ], lead="edmval ~1e-27 in all three; 72 of 72 basis parameters correlate. Prior units.",
        widths=[0.26, 0.13, 0.12, 0.12, 0.12, 0.11, 0.14],
        note="B⁺ on top of Z + J/ψ: σ(mZ) −0.2 %, σ(e) −0.5 %, σ(A) −0.3 %, σ(M) 0. Z + B⁺ alone leaves e at 0.65: the J/ψ carries e.")
    d.figure("σ(A), σ(e), σ(M) per eta group in the three configurations", "aem_sigma_per_group.png",
             lead="Without the J/ψ, e is 2.5–10× less constrained and A 1.2–2.5×.")
    d.figure("What B⁺ adds on top of Z + J/ψ, per eta group", "gain_aem_sigma_per_group.png",
             lead="Ratio panels: e −3 % at |η| = 2.3, −0.2 % centrally; A −1.3 % → −0.1 %; M untouched.")
    d.bullets("The correlations between the scale parameters", [
        "ρ(A_i, e_i) posterior: 0.854 (Z + J/ψ), 0.767 (Z + B⁺), 0.853 (Z + J/ψ + B⁺); range 0.58–0.91 across eta",
        "Prior-stripped, F = V⁻¹ − I, all three give 0.85–0.86: the degeneracy is the physics, the posterior is the prior's dilution",
        "A–M and e–M stay below 0.16 everywhere; adjacent eta pairs: ≈ 0 for A and e, +0.1 to +0.25 for M (backup)",
        "mZ correlates with A at +0.21 and with e at +0.18 in both J/ψ configurations, with M at 0",
        ("Z + B⁺: ρ(mZ, A) +0.31 and ρ(mZ, e) 0.00 — with e unconstrained, the Z's own scale information goes into A", 1),
    ])
    d.figure_row("A/e/M correlations at η ≈ 0, with and without B⁺",
                 ["corr_barrel_noB.png", "corr_barrel_full.png"],
                 lead="Z + J/ψ left, Z + J/ψ + B⁺ right, group 12: identical to the third digit; ρ(A,e) 0.91.")
    d.figure_row("A/e/M correlations at |η| = 2.3: Z + J/ψ vs Z + B⁺",
                 ["corr_endcap_noB.png", "corr_endcap_noJpsi.png"],
                 lead="Group 0: ρ(A,e) 0.76 with the J/ψ (left), 0.58 with B⁺ alone (right), where σ(e) = 0.65.")
    d.figure_row("Grouped impacts on mZ, with and without B⁺",
                 ["impacts_grouped_noB.png", "impacts_grouped_full.png"],
                 lead="Z + J/ψ left, Z + J/ψ + B⁺ right: A 1.70 MeV, e 1.57 MeV in both; nothing else moves.")
    d.figure("Grouped impacts on mZ, Z + B⁺ alone", "impacts_grouped_noJpsi.png",
             lead="σ(mZ) 5.1 MeV: without the J/ψ the scale groups are the whole budget.")
    d.figure_row("Leading nuisances on mZ, with and without B⁺",
                 ["impacts_mZ_noB.png", "impacts_mZ_full.png"],
                 lead="Top 25: the central-eta A's at 0.5 MeV each lead in both; the ordering does not change.")
    d.figure_row("Z channel, Z + J/ψ + B⁺: prefit and postfit", ["prefit_Z_full.png", "postfit_Z_full.png"],
                 lead="Asimov: data equals prediction; the postfit band is the residual scale uncertainty.")
    d.figure_row("J/ψ prompt and displaced, postfit, Z + J/ψ + B⁺", ["postfit_jpsi_prompt_full.png", "postfit_jpsi_displaced_full.png"])
    d.figure_row("B⁺ channel postfit: Z + B⁺ (left), Z + J/ψ + B⁺ (right)",
                 ["postfit_B_noJpsi.png", "postfit_B_full.png"],
                 lead="Summed over 144 cells; peaking template and per-eta widths in both.")
    d.figure_row("Constraints against eta, with and without B⁺",
                 ["constraints_noB.png", "constraints_full.png"])
    d.bullets("Findings", [
        "B⁺ with the converged lineshape changes the combined result by 0.2 % on σ(mZ), 0.5 % on mean σ(e), 0.3 % on σ(A): the August conclusion holds",
        "Where it helps is the endcap e (−3 %), where the e variation per prior sigma is ten times larger than centrally",
        "Z + B⁺ without the J/ψ is not a substitute: σ(mZ) 5.1 MeV against 2.2, σ(e) 0.65 against 0.11",
        "The scale correlations are unchanged by B⁺ and, prior-stripped, sit at 0.85–0.86 in every configuration",
        ("Every fit unblinded by construction (Asimov); plots by rabbit with --titlePos 0, per-channel legends upper left", 1),
    ])
    return d


def add_backup(d):
    d.section("Backup")
    d.figure("The full 72-parameter scale block, Z + J/ψ + B⁺", "corr_scale_all_full.png",
             lead="The off-diagonal stripe is the A_i–e_i correlation; A–M and e–M are blank.")
    d.figure_row("Adjacent-pair and mZ correlations vs eta, Z + J/ψ + B⁺",
                 ["pair_corr_full.png", "mz_corr_full.png"])
    d.bullets("What was done to produce the plots", [
        "Driver: joint_with_btojpsik.sh gained --noB, --refTag, --plotDir; the B⁺ tensor stage carries the converged model",
        "Plot stage is a separate script (joint_bplus_plots.sh): one rabbit_plot_hists call per channel with --filterProcs",
        ("a single call legends every process of the merged tensor on every panel — eight entries over the Z peak", 1),
        "--titlePos 0, legend upper left, one column; impacts without printed numbers (they overlap the bars)",
        "Cross-configuration σ and ρ figures from cvhnano_resolution_report.py, which strips the prior for ρ(A,e)",
        "Merged fits: 13 470 / 13 734 parameters, --hessianChunkSize 256, about 35 minutes each on one core",
    ])
    return d


if __name__ == "__main__":
    deck = build()
    deck = add_results(deck)
    deck = add_backup(deck)
    print(deck.save("joint_mz_jpsi_bplus.pptx"))
