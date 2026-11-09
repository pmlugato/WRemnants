"""Deck for the kaon/pion ONNX response model study.

openspec/changes/add-kpi-onnx-reweights section 7. Plot-led throughout: the
study's results are shapes and comparisons, and the numeric scans belong in
backup.

    python scripts/btojpsik/make_kpi_deck.py
"""

import os
import sys

sys.path.insert(0, "/home/submit/pmlugato/.claude/skills/slides")
from deck import Deck  # noqa: E402

ASSETS = os.path.join(os.path.dirname(__file__), "..", "..", "slides",
                      "kpi_assets")


def build():
    d = Deck(assets=ASSETS)

    d.title_slide(
        "The kaon/pion ONNX response model",
        "What it changes, and what it says about the model we were using",
        "Pietro Lugato  ·  8 September 2026  ·  mZ / B → J/ψ K",
    )

    d.bullets("What this study set out to do", [
        "Havyn trained a response-reweight network carrying charged hadrons: "
        "wremnants-data d5b8414",
        "Every B⁺ conclusion so far used a network with no hadron class, "
        "trained an order of magnitude harder than our kaons",
        "The question: does the new model change any of those conclusions",
        ("It does — but the largest change came from something else the "
         "study turned up on the way", 1),
    ])

    d.bullets("Three structural facts, before any fit", [
        "The class alphabet changed completely: {1, 15, 443} became "
        "{211, 321, other}, and our adaptor turns 321 into 1 silently",
        "Training ranges are disjoint: mean gen pᴛ 28.2 GeV for the muon "
        "family, 3.04 GeV for kpi",
        "Valid windows overlap only around 5–15 GeV, and there the two "
        "disagree by a factor 2–5",
    ])

    d.figure(
        "Response slope: the two models cover different momenta",
        "slope_vs_pt.png",
        lead="kpi reverses sign below ~0.9 and above ~25 GeV: not for Z.",
    )

    d.section("What we found about the model we were already using")

    d.bullets("The B⁺ channel was loading the wrong ONNX variant", [
        "The histmaker's own comment says the channel smears nothing and needs "
        "the un-smeared model",
        "It then computed smearing = not args.noSmearing, which defaults True",
        "Every production before 2026-09-08 ran the SMEARED-trained model on "
        "an unsmeared sample",
        ("Fixed. --useSmearedTrainedModel reproduces the old numbers for "
         "comparison", 1),
    ])

    d.figure_row(
        "The 2 GeV cliff belongs to that wrong variant",
        ["cliff_vs_yraw.png"],
    )

    d.bullets_figure(
        "And y_raw = 0 is the worst place to probe",
        ["Steep odd function of the residual, amplitude ±250, zero crossing at y_raw = 0",
         "A pᴛ trend there is what is left after a large offset cancels",
         "_smearing breaks at every residual; _nosmearing is flat at every one"],
        "slope_vs_yraw.png",
    )

    d.section("The acceptance gate")

    d.bullets("The gate we inherited has a false premise", [
        "Section 18: “A is pᴛ-independent, so its induced shift must "
        "be flat in pᴛ”",
        "That does not follow: A scales the kaon momentum by a constant "
        "fraction, but the mass responds through dm/dln(p_K)",
        "Measured: the shift rises from ~0.45 to 1.0 across pᴛ in EVERY "
        "configuration, including the recorded one",
        "An absolute-flatness test fails everything, so it cannot be the gate",
    ])

    d.figure(
        "One configuration is not like the others",
        "aflatness_gate.png",
        lead="Only _smearing unclamped turns over and changes sign.",
    )

    d.bullets("Three criteria tried; the third is the one that holds", [
        "Absolute flatness — fails everything, premise is wrong",
        "Departure from a baseline profile — circular: the baseline must "
        "be clamped, so unclamped configs fail by construction",
        "Monotonicity and sign stability — not circular, and it catches "
        "the failure that matters",
        ("kpi passes with no conditioning clamp at all", 1),
    ])

    d.figure(
        "Per fit cell, the variation is well ordered in A and e",
        "cell_shift_A.png",
        lead="Sign disagreement with the majority: 0.8–2.6% in A, ~49% in M.",
    )

    d.section("Does it change any conclusion?")

    d.figure(
        "Both fits say the same thing: kpi costs a little",
        "kpi_fit_summary.png",
        lead="All eight fits converged. The recorded joint result reproduces "
             "exactly at 0.32%.",
    )

    d.bullets("Standalone B⁺ Asimov, four configurations", [
        "kpi gives σ(e) 5.1% LARGER — less constraining, not more",
        "Predicted before the fits: kpi's per-cell shifts are 15–16% smaller, "
        "and a smaller variation reads as a larger σ",
        "It does break the A–e degeneracy a little better, ρ 0.229 → 0.200",
        ("The smearing bug moves σ by under 2%: real, but it does not "
         "invalidate the recorded numbers", 1),
    ])

    d.bullets("Config D answers the muon-leg question outright", [
        "Swapping the entire response model on both muon legs changes "
        "σ(e) by 0.05% and ρ by 0.001",
        "In the joint fit, D and B agree to four digits",
        "So whether kpi's class 0 is a trained muon class or a residual "
        "bucket does not matter",
        ("Cleaner than --scaleVarLegs bachelor: that removes the legs, "
         "this re-describes them and still finds nothing", 1),
    ])

    d.table(
        "What survives and what is overturned",
        [["conclusion", "status"],
         ["B⁺ improves σ(mZ) by ~0.3%", "survives, 0.25–0.35%"],
         ["B⁺ does not break the A–e degeneracy", "survives, 0.854 → 0.852"],
         ["72/72 basis parameters correlate", "survives, all configs"],
         ["2 GeV cliff is the training floor", "OVERTURNED — _smearing only"],
         ["443-vs-1 swing is a factor 31", "OVERTURNED — 1.4 corrected"],
         ["A's induced shift must be flat in pᴛ", "OVERTURNED — rises kinematically"],
         ["soft unclamped kaons are an artefact", "stands, wrong stated reason"]],
        widths=[0.55, 0.45],
    )

    d.bullets("Where that leaves us", [
        "kpi is not a precision win for B⁺, and it was never going to be: "
        "B⁺ is a cross-check",
        "It is the right model on physics grounds — real kaon class, trained "
        "where our kaons live, no clamp needed",
        "The real lever remains the DATA fit, and the selection: the kaon ID "
        "is measured at χ²/ndf 41.2 → 25.4 and is in no fit",
        "Open: kpi on the J/ψ legs (config C), and Havyn's remaining questions",
    ])

    d.section("Backup")

    d.figure("Slope against eta", "slope_vs_eta.png")
    d.figure("log W against the variation size", "logw_shape.png")
    d.figure("Per-cell shift, e", "cell_shift_e.png")
    d.figure("Per-cell shift, M — noise in every configuration",
             "cell_shift_M.png")

    d.table(
        "Response slope, eta = 0.6, y_raw = 0",
        [["pᴛ [GeV]", "_smear 443", "_nosm 443", "kpi 321", "kpi 13"],
         ["1.0", "509", "68.4", "3.1", "11.5"],
         ["1.5", "387", "46.4", "13.9", "16.2"],
         ["1.95", "18", "33.8", "16.5", "16.2"],
         ["5.0", "6.1", "5.8", "9.9", "11.4"],
         ["25", "1.7", "-0.6", "0.5", "0.6"],
         ["45", "4.2", "-4.7", "-20.1", "-19.5"]],
        note="The factor-21 break at 1.5 → 1.95 is _smearing only.",
    )

    d.table(
        "Plumbing regression vs the pre-change HEAD",
        [["comparison", "cells differing", "max rel", "verdict"],
         ["old vs --useSmearedTrainedModel", "10 / 184.3M", "2.2e-16", "same"],
         ["old vs new default", "58,560", "8.2e-2", "different"]],
        note="Bit-exactness is unattainable with -j > 1.",
    )

    return d


if __name__ == "__main__":
    print(build().save("260908_kpi_onnx_reweights.pptx"))
