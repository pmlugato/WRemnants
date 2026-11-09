#!/usr/bin/env python3
"""Deck: A/e/M momentum-scale variations in the CVH NanoAOD B+ -> J/psi K channel.

Covers how the ONNX response reweight works and what it is entitled to, the
merge of havyn's implementation into our fork, the e-convention question and its
resolution, and the first measurements: variation size, what constrains what, and
the expected data yield.
"""

import sys

sys.path.insert(0, "/home/submit/pmlugato/.claude/skills/slides")

from deck import Deck  # noqa: E402

ASSETS = "/home/submit/pmlugato/public_html/slides/assets_aem_variations"

d = Deck(assets=ASSETS)

d.title_slide(
    "The kaon scale variations now exist",
    "ONNX response reweight, mass-aware, fitted on the full 2016 sample",
    "Pietro Lugato  ·  25 August 2026  ·  mZ / B → J/ψ K",
)

d.bullets(
    "Read this first: no data fit has converged",
    [
        "The reference data fit stopped at iteration 7199 with \"Maximum number of "
        "iterations has been exceeded\" and edmval 72 668, against a 1e-5 scale",
        "It wrote results anyway: rabbit writes whenever the Hessian passes Cholesky, "
        "and a Hessian away from the minimum factorises fine while being no covariance",
        "So every data-side number here is provisional — the postfit χ², the "
        "correlations, the pulls, the impacts",
        "Unaffected: all Asimov results (edmval ~1e-28), the background family and "
        "window scans, the 17% width deficit, and everything gen-level",
        "A convergence gate is now in the fit script, and toys identify the cause: "
        "model mis-specification, not the parameterisation",
    ],
)

d.bullets(
    "Outline",
    [
        "How the reweight works, and what it is entitled to",
        "The merge, and the kaon mass it silently deleted",
        "The 1/cosh η question, resolved from our own notes",
        "What the variations look like, and what constrains what",
        "The fits: sensitivity, the degeneracy that is not there, and the data",
    ],
)

d.bullets(
    "Where this started: a fit with nothing to fit",
    [
        "The CVH-NanoAOD chain ran end to end and carried not one calibration nuisance",
        "The channel exists to constrain the kaon A, e and M — and had no A, e or M",
        "Gap G2 closed (per-leg gen matching), and havyn had already built the "
        "machinery — for D0 → K π, a hadronic two-body decay with per-leg masses",
        "So: merge his work rather than rebuild it, then measure",
        ("Everything below is the full sample: 24 847 files, 95.5M events, "
         "1.16M signal candidates in data.", 1),
    ],
)

# --- how the thing works -------------------------------------------------

d.figure(
    "How the reweight works",
    "flow_schematic.png",
    lead="A trained density ratio, not a linearised derivative.",
)

d.bullets(
    "What the network is, and what it is entitled to",
    [
        "It learned p(y | c): the density of the response residual "
        "y = (r_κ, δλ, δφ) given a conditioning vector c",
        "Asking for a variation means asking for a ratio: how much more likely is "
        "this response if the A/e/M shift u had been applied. That is exp(log r).",
        "It never recomputes the reconstructed mass. It reweights events.",
        "c carries log pt_gen, charge, λ, sin φ, cos φ and an integer class tag "
        "— 1 prompt W/Z, 15 τ-decay, 443 J/ψ leg. There is no hadron class.",
        ("Two consequences, both measured later: the class tag matters enormously, "
         "and the mass constraint does not suppress the muon legs.", 1),
    ],
)

d.figure(
    "The conditioning floor, and what is well behaved",
    "floor_and_linearity.png",
    lead="Training floors gen pt at 2 GeV, so a 1–2 GeV kaon is out of domain.",
)

d.figure(
    "Why kaon pT > 1 GeV, twice over",
    "massaware_and_margin.png",
    lead="Soft kaons carry the e sensitivity, and are where the exact shift breaks.",
)

# --- the merge and the mass ----------------------------------------------

d.bullets(
    "The merge compiled and deleted the kaon mass",
    [
        "70 havyn-only commits against 235 of ours, and only 12 files touched by "
        "both — of which two mattered",
        "Both sides had added a 7-argument mass-aware calculateQopUnc with an "
        "identical signature, so the header would not compile at all",
        "Resolving it left particle_masses a dead parameter: btojpsik.py still passed "
        "it, the constructor defaulted it to {} — massless, no error, no log line",
        "The 2018 histmaker had lost its kaon energy-loss term and nothing said so",
        ("Fixed, unified to one rule on both backends, and pinned by a test.", 1),
    ],
)

# --- the convention ------------------------------------------------------

d.bullets(
    "The 1/cosh η question, and what it actually is",
    [
        "Two implementations of the mass-aware e term differed by 1/cosh η — 1.0 at "
        "η 0, 2.15 at 1.4, 5.56 at 2.4. Not a rounding difference",
        "e is an energy loss, and Δp = (E/p) ΔE. Mapping that to a curvature shift "
        "brings a cos λ = 1/cosh η — so the physics does contain the factor",
        "But e is fitted per η bin, so that factor was already absorbed into "
        "e_i(η) at fit time. Reapplying it double counts",
        "The file's own central-value corrector settles it: it applies "
        "(1 + A − e·k)·k with no cosh η. So e is a shift of pT, not of energy",
        ("Only E/p is genuinely mass-dependent and has to be reinstated: 1.115 at "
         "pT 1 GeV, 1.029 at 2 GeV, 1.002 at 8 GeV.", 1),
    ],
)

d.figure(
    "So the pT convention is the right one, and it matters",
    "e_convention.png",
    lead="Now a switch, --eConvention, defaulting to pT. Everything below uses it.",
)

# --- results -------------------------------------------------------------

d.figure(
    "The variation looks like a mass scale should",
    "template_e_variation.png",
    lead="e in η group 2: an antisymmetric ±0.6% tilt across the peak.",
)

d.figure(
    "It constrains e, weakly A, and barely sees M",
    "sensitivity_and_muons.png",
    lead="S = 1 means constrained to one σ of its input uncertainty.",
)

d.bullets(
    "The muon legs are not a small correction",
    [
        "The expectation was: the joint CVH fit constrains the dimuon to the J/ψ "
        "mass, so varying the muons cannot move the B mass",
        "That is about the right quantity and still fails — the reweight varies the "
        "muon q/p and asks how likely that response is. The weight exists anyway.",
        "Measured: muons carry 38–46% of A, 60–68% of M, and 11–18% of e",
        "So on A and M this channel largely re-measures what Z and J/ψ already do "
        "much better; on the kaon e term the bachelor carries 82–89%",
        ("Read together with the table above: e is the measurement.", 1),
    ],
)

d.figure(
    "The network is conditioned on a label we get wrong",
    "class_tag.png",
    lead="443 (J/ψ leg) against 1 (prompt W/Z): S(e) moves by 280–1050%.",
)

# --- the fits ------------------------------------------------------------

d.bullets(
    "The variations are now nuisances in the fit",
    [
        "btojpsik_cvhnano_tensor.py gained --scaleVariations: one constrained "
        "nuisance per (A/e/M, η group), 12 at the default 4 groups",
        "The variation histogram goes through the *same* reduce_fitbins call as the "
        "nominal — overflow fold, mass window, pT rebin — so the two cannot drift",
        "Largest bin shift in the tensor: 3.5–6.2% for A, 10.0–18.8% for e, "
        "1.1–1.5% for M, with template-integral shifts of 0.02–0.4%",
        "So these are shape nuisances, which matters because the per-cell signal "
        "normalisation is free and would otherwise eat them",
        ("Asimov fit closes: χ²/ndf = 0, all 12 pulls at zero, no boundary.", 1),
        "Data: 1 918 760 candidates in the window, 1 162 080 signal, purity 0.61",
    ],
)

d.figure(
    "The simulated sample is now what limits us",
    "fit_constraints.png",
    lead="Full sample, 16 cells. The simulated template is what limits us now.",
)

d.figure(
    "The A–e degeneracy is there — I read the wrong number",
    "ae_correlation.png",
    lead="Data-only ρ is stable at 0.87; the posterior moved 0.10 → 0.60.",
)

d.bullets(
    "Why the postfit correlation hid it",
    [
        "These are constrained nuisances with a unit Gaussian prior, so the "
        "posterior precision is P = I + F with F the information from data",
        "When F is small the posterior is prior-dominated and its correlations are "
        "suppressed — the prior is diagonal and it wins",
        "Strip it: F = V⁻¹ − I, then the data-only covariance is F⁻¹. Its "
        "correlation is 0.86–0.88, against the 2018 value of about 0.9",
        "The 2D contours are Hessian-exact, so the likelihood is Gaussian and this "
        "algebra is exact rather than approximate",
        ("A–M and e–M stay below 0.16, so it is specifically the A–e pair.", 1),
    ],
)

d.figure(
    "And you can see it in the variations",
    "variation_A_eta.png",
    lead="A-down tilts one way, e-down the other, so +A and +e partly cancel.",
)

d.figure(
    "The e variation, all four η groups",
    "variation_e_eta.png",
    lead="±5–9% antisymmetric tilt, crossing zero at the B mass. 10 MeV bins.",
)

d.figure(
    "The e variation by pT and charge",
    "variation_e_cells.png",
    lead="Central η group. Soft cells carry ~2× the amplitude; K⁺ and K⁻ agree.",
)

d.figure(
    "M is a few tenths of a percent, everywhere",
    "variation_M_eta.png",
    lead="Same axes as the e figure. This is why M never gets constrained.",
)

d.bullets(
    "What the fit costs, and what it does not",
    [
        "Removing the simulated template's own statistics raises S(A) by 3.4×, "
        "S(e) by 2.3× and S(M) by 3.8× — it is the dominant term now",
        "62 448 simulated signal candidates against 1 162 080 in the data. More "
        "data buys almost nothing; more simulation buys a factor 2–4",
        "Not the background lineshape: exponential against Bernstein gives S(e) "
        "identical to three figures on the same Asimov",
        "Finer cells are affordable: 4 pT × 8 η × 2 charge gives S(e) 2.52–2.91 "
        "against 2.34–2.69, with the biggest gain on A",
        ("At 200 files the template cost 20–25%. I flagged it as worth watching; "
         "it is now the answer.", 1),
    ],
)

d.table(
    "The measurement, on the full sample",
    [
        ["parameter", "posterior S", "data only, A–e free", "200-file guess"],
        ["A", "0.58 – 0.63", "1.78 – 1.96 × input", "0.62 – 0.72"],
        ["e", "2.34 – 2.69", "0.58 – 0.67 × input", "0.24 – 0.27"],
        ["M", "0.25 – 0.28", "3.63 – 3.96 × input", "1.28 – 1.38"],
    ],
    lead="Measured, not projected. The middle column is the honest one.",
    note="This channel measures the kaon e term to about 60% of its input "
    "uncertainty, and does not measure A or M. The 200-file projection was "
    "optimistic by 2.7× across the board, because it scaled by the data file count "
    "(112×) while the simulation only grew 12.6×.",
)

d.figure(
    "Only e is a measurement",
    "fit_projection.png",
    lead="Below the line is a measurement. Only e is.",
)

d.table(
    "The observation still does not fit, and it is not stats",
    [
        ["configuration", "cells", "free norms", "outcome"],
        ["wide, exponential", "16", "16", "Hessian not pos.-def"],
        ["wide, Bernstein", "16", "16", "Hessian not pos.-def"],
        ["narrow, exponential", "16", "16", "Hessian not pos.-def"],
        ["η collapsed, wide", "4", "4", "edm 3e−13, but nonsense"],
        ["wide, norms grouped on η", "16", "4", "Hessian fine, edm 7e+04"],
    ],
    lead="At 200 files I blamed sparse cells. 159× more data did not fix it.",
    note="The last row is decisive: same 16 cells, only the signal normalisations "
    "cut from 16 to 4, and the Hessian becomes well behaved. The degenerate "
    "direction is per-cell signal normalisation against per-cell background "
    "amplitude.",
)

d.figure(
    "The model is percent-level wrong, and 6.2M notice",
    "postfit_data_full.png",
    lead="χ²/bin = 203 cell-summed, 34 per cell. The scale nuisances pull 6–8σ.",
)

d.bullets(
    "Open questions for the meeting",
    [
        "The e convention: I now read the correction file as pT-level from its own "
        "central-value corrector. Worth confirming with David and havyn, since havyn "
        "calls his the 'new cosh η convention' and may have refitted",
        "Is a kaon-trained response model worth doing? The class tag swings S(e) by "
        "up to a factor 11, and the floor by a factor 3 where it acts",
        "Is calinput-on-pr692 the branch to build on, or will it be rebased?",
        "Should this be presented as a kaon-e measurement rather than an A/e/M one?",
        "Which lineshape fix: a wider signal template, or a non-monotonic "
        "background? The no-b-ancestor bin decides it and the measurement is one "
        "run away",
    ],
)

d.bullets(
    "Next, in order",
    [
        "The lineshape, both halves of it: one degree of curvature in the "
        "background (worth χ² factor 2.0) and a resolution nuisance on the signal "
        "(1.16). Bernstein with the norms grouped on η, plus a smearing parameter",
        "More simulation, or a smoothed signal template. The template's own "
        "statistics cost a factor 2–4 and 19× more data will not help",
        "Then re-read the pulls, take the finer 64-cell binning, and quantile-bin "
        "in mass to keep the sparse cells populated",
        ("Only then is the 48-η-bin basis the joint mZ + J/ψ fit wants worth "
         "asking about.", 1),
    ],
)

d.bullets(
    "Three new readouts, and they agree with each other",
    [
        "The nuisance correlation matrix, from rabbit's own tool, sliced three "
        "ways so it can be read and quoted",
        "The postfit residual in all sixteen cells on one canvas, instead of "
        "sixteen files nobody compares",
        "The simulated non-signal mass shape, decomposed by truth category",
        ("All three say the same thing: the model's shape is wrong, and it is "
         "wrong in one specific way.", 1),
    ],
)

d.figure_row(
    "A–e is +0.6 on Asimov and +0.63 on data, in the same place",
    ["corr_asimov_eta0.png", "corr_data_eta0.png"],
    lead="η group 0. A–M and e–M are below 0.005; M is measured by something else.",
)

d.figure(
    "Cross-η elements are zero: four separate measurements",
    "corr_data_Ae_all.png",
    lead="A(0–3) then e(0–3), data. A_i correlates with e_i and with nothing else.",
)

d.bullets(
    "Why the posterior correlation moved from 0.10 to 0.60",
    [
        "P = I + F: the posterior mixes the unit prior with the data information",
        "At 200 files F was small, so the posterior looked like the prior — diagonal",
        "At full statistics F dominates and the data's own +0.86 shows through",
        "Stripping the prior gives +0.86 to +0.88 at both sample sizes, which is "
        "the 2018 number",
        ("The 0.10 I first reported was the prior, not the physics.", 1),
    ],
)

d.figure_row(
    "One cell, prefit and postfit, as rabbit draws it",
    ["cell_prefit_data.png", "cell_postfit_data.png"],
    lead="pT bin 0, η bin 0, charge 0. The seed is crude; the fit fixes the apex "
         "and not the flanks.",
)

d.figure(
    "The postfit residual is one W shape, repeated in every cell",
    "residual_lowpt_data.png",
    lead="+3.5% on both flanks, 1.000 in the peak. Identical across η and charge.",
)

d.table(
    "The residual, in numbers",
    [["projection", "prefit χ²/ndf", "postfit χ²/ndf", "max |ratio − 1|"],
     ["16 cells × 80 mass bins", "122.1", "34.05", "1.15"],
     ["cell-summed mass, 80 bins", "853.5", "203.1", "0.226"],
     ["4 cells, eta collapsed, summed", "853.5", "82.7", "0.160"],
     ["Asimov, any projection", "0.00", "0.00", "0.000"]],
    lead="Converged norm-grouped data fit, 1.9M candidates in the peak.",
    widths=[0.40, 0.20, 0.20, 0.20],
    note="Summing cells is the harder test, not the weaker one. The Asimov row is "
         "what makes this a data/model statement, not a convergence one.",
)

d.figure(
    "The simulated background is a bump, not an exponential",
    "bkg_composition.png",
    lead="χ²/ndf = 28.2 against its own best-fit exponential; +16.4% in the peak.",
)

d.table(
    "Who carries it: the no-b-ancestor bin, and nothing else",
    [["category", "in window", "fraction", "in peak", "peak/window"],
     ["BuJpsiK (signal)", "72 790", "0.259", "65 296", "0.897"],
     ["NoBAncestor", "192 686", "0.686", "46 539", "0.242"],
     ["Bd", "7 375", "0.026", "163", "0.022"],
     ["BuJpsiKX", "3 662", "0.013", "8", "0.002"],
     ["BuJpsiPiX", "3 573", "0.013", "534", "0.150"],
     ["Bs", "693", "0.003", "101", "0.145"]],
    lead="Simulation, [5.05, 5.45], forced BRs undone.",
    widths=[0.30, 0.18, 0.16, 0.18, 0.18],
    note="The exponential's residual in the peak is +7 654 candidates, all of it "
         "in NoBAncestor: a 12% peaking component the fit cannot represent.",
)

d.figure(
    "Which component owns the residual? Profile the norms out",
    "residual_attribution.png",
    lead="Postfit shapes, both normalisations profiled at every trial, so only "
         "shape can win.",
)

d.table(
    "The background's curvature is the larger error",
    [["change to the postfit model", "χ² (80 bins)", "factor"],
     ["nothing", "43 047", "1.00"],
     ["signal convolved with 7.3 MeV", "37 050", "1.16"],
     ["background × (1 + linear + quadratic)", "21 311", "2.02"],
     ["background + a 137 MeV Gaussian", "25 730", "1.67"],
     ["both: curvature + 6.1 MeV on the signal", "20 150", "2.14"]],
    lead="Visible without fitting: the postfit background is flat, the data's arcs.",
    widths=[0.52, 0.26, 0.22],
    note="The background gets three free shape parameters here and the signal one, "
         "so this is the direction, not a likelihood ratio.",
)

d.bullets(
    "A correction to what I said earlier today",
    [
        "I first read the whole residual as the signal being 18% too narrow",
        "It is not: the signal width buys a factor 1.16 and the background "
        "curvature a factor 2.02",
        "The width mismatch is real and separately measured — it is just not the "
        "dominant term",
        ("And convolving by the 12.17 MeV the first fit suggested made things "
         "worse, which is what gave the error away.", 1),
    ],
)

d.figure(
    "The signal core is 17% too narrow, and that survives",
    "width_comparison_curved.png",
    lead="Narrow Gaussian + exponential + one broad Gaussian, same form on all three.",
)

d.table(
    "The width ratio is robust to the background form",
    [["background model", "MC template σ", "data σ", "ratio"],
     ["exponential only", "19.60 ± 0.06", "23.07 ± 0.03", "1.1767 ± 0.0041"],
     ["+ one broad Gaussian", "14.81 ± 0.11", "17.29 ± 0.04", "1.1679 ± 0.0093"]],
    lead="Both fits absorb the curvature symmetrically.",
    widths=[0.30, 0.24, 0.22, 0.24],
    note="Missing convolution width √(17.29² − 14.81²) = 8.9 MeV, matching the "
         "6–9 MeV the profiled scan prefers. The uncurved 12.17 MeV was too big.",
)

d.bullets(
    "The scale nuisances stand in for a missing parameter",
    [
        "A, e and M are scale parameters: they shift the peak, they do not widen it",
        "Resolution variations are deferred, so the fit owns no resolution handle",
        "Handed a 17% width mismatch and a curved background, it spends the only "
        "shape handles it has",
        ("That is the 6–8σ pulls. Nothing downstream is a measurement until both "
         "are fixed.", 1),
    ],
)

d.bullets(
    "The good news, in the same fit",
    [
        "The fitted means differ by −1.14 MeV, data below simulation: 2.2×10⁻⁴",
        "That is a momentum-scale difference, and it is exactly what A/e/M exist "
        "to describe",
        "So the signal the channel is looking for is there, under two shape "
        "errors sitting on top of it",
        ("Fix the lineshape and the −1.14 MeV becomes the measurement.", 1),
    ],
)

d.bullets(
    "What to build, and what not to",
    [
        "Do not tighten selections: the background level is already right to 2%, "
        "and no selection fixes a shape",
        "Do not add peaking channels: B → J/ψπ and Bs together are 635 candidates "
        "in the peak against 65 296 signal",
        "The no-b-ancestor structure is 137 MeV wide — nine times the signal core "
        "— so it is the combinatorial's own arc, not a peak",
        "Do build both: one degree of curvature in the background, and a "
        "resolution nuisance on the signal",
        ("Even both together leave χ²/bin = 252 and a 13.8% maximum residual. The "
         "5.10–5.20 shoulder is still unmodelled.", 1),
    ],
)

d.bullets(
    "The composition is right; only the shape is wrong",
    [
        "Simulated S/B in the peak 1.379, data's sideband-extrapolated S/B 1.354",
        "So the simulation predicts the data's signal fraction to 2%",
        "And gets the shape wrong by 3.5%, coherently, in every cell",
        "Tightening selections lowers the background level, which is the thing "
        "that is already right",
        ("The exponential is rejected by the data on the sidebands alone, peak "
         "excluded: χ²/ndf = 87 705/62 = 1415.", 1),
    ],
)

d.bullets(
    "Is the signal template biased by how we define signal?",
    [
        "Worth asking: the template is the set of candidates whose legs the "
        "producer's matcher accepted, and the matcher cuts on |Δpt|/pt < 0.1",
        "That is a cut on resolution, and resolution is what we measure",
        "Measured: the bachelor residual has rms 0.0180, so the cut is at 5.5σ",
        "0.45% sits in the last fifth before the cut, 0 beyond it, and the rms is "
        "flat in pT — 5.2 to 5.7σ in every one of eighteen bins",
        ("Answered and closed. The match is not sculpting the template.", 1),
    ],
)

d.figure(
    "Gen truth with no reco requirement, against what we select",
    "gen_vs_reco_spectrum.png",
    lead="714 190 gen decays, 572 082 in acceptance, 62 422 called signal.",
)

d.table(
    "The informative η cells are the least populated",
    [["η group", "|η|", "cosh η at outer edge", "efficiency"],
     ["1, 2", "0.0 – 1.2", "1.81", "0.1417, 0.1406"],
     ["0, 3", "1.2 – 2.4", "5.56", "0.0884, 0.0905"]],
    lead="Efficiency with respect to acceptance, from the gen-level denominator.",
    widths=[0.20, 0.24, 0.30, 0.26],
    note="The forward groups carry the cosh η lever arm on e and are reconstructed "
         "at 0.63× the central rate. An argument for keeping them, not collapsing η.",
)

d.figure(
    "The no-b-ancestor bin, and the cut it did not fail",
    "gen_truth_diag.png",
    lead="89.7% is three legs matched with no ancestor, which is correct here.",
)

d.table(
    "Our gen matching against Bmm5, line by line",
    [["Bmm5", "ours", "verdict"],
     ["|Δpt|/pt < 0.1 and dR < 0.02", "same two cuts, same values", "identical"],
     ["is_acceptable: skip |pdg|<10, 2212, 21", "same three", "identical"],
     ["get_mother, find_common_ancestor, depth 10", "same, plus an n=0 guard", "ours safer"],
     ["packedGenParticles (status 1 by construction)", "genParticles with status≠1 skip", "equivalent"],
     ["no charge requirement", "requires charge agreement", "ours tighter"],
     ["first passing candidate", "smallest dR among passing", "ours tighter"],
     ["prefers the PAT embedded gen ref", "no such ref exists", "the one real difference"]],
    lead="The producer says 'ported from Bmm5' with line citations; this checks it.",
    widths=[0.38, 0.36, 0.26],
    note="Bmm5's loose mode (|Δpt|/pt<0.2, dR<0.1) is used only for a Ks composite, "
         "never for a track leg — so on legs it is exactly as tight as we are.",
)

d.section("Backup")

d.table(
    "The merge: small overlap, two files that mattered",
    [
        ["", "count", "note"],
        ["havyn-only commits", "70", "calinput-on-pr692"],
        ["ours-only commits", "235", ""],
        ["files touched by both", "12", "of which 2 mattered"],
        ["submodules", "3", "narf, wremnants-data theirs; rabbit ours"],
        ["shared utilities", "3", "ours 2–7 lines each"],
        ["histmakers we do not use", "3", "recorded before resolving"],
    ],
    lead="A real merge on a branch, not a cherry-pick.",
    note="Branch onnx-weight-catchup, merge commit 844df06. havyn is still developing "
    "this branch, so later updates stay ordinary merges. rabbit kept as ours: theirs "
    "is older and not reachable from the declared submodule URL.",
)

d.table(
    "The regression test that would have caught it",
    [
        ["comparison", "result", "what it catches"],
        ["massless variation size", "3.32e−3", "reference"],
        ["[] vs [0,0,0]", "bit-identical", "zero-mass leg, wrong overload"],
        ["[m_K,0,0] vs [0,0,0]", "9.9e−5  =  3.0%", "mass not reaching leg 0"],
        ["[m_K,0,0] vs [m_K,m_K,m_K]", "6.7e−8", "the broadcast trap"],
    ],
    lead="One candidate: bachelor 1.5 GeV at η 0.8, muons 6 and 4 GeV.",
    note="The mass enters only as (E/p − 1) on the e term: the bachelor has "
    "p = 2.006 GeV, so E/p = 1.0299 and 3.0% is what should come out. The broadcast "
    "trap is tiny because the muons carry 21 and 32 GeV.",
)

d.bullets(
    "How much to trust the 1–2 GeV slice",
    [
        "The floor moves only log pt_gen in c. The response residual and the "
        "physical shift keep the true kinematics — so nothing unphysical is fed in",
        "But it is control, not validation: with the floor off, S(e) in the 1–2 GeV "
        "slice goes 0.35 → 1.05 and single weights move by 138%",
        "Above 2 GeV the floor is inactive and the two agree to 8e−4, which confirms "
        "the switch does nothing where it should do nothing",
        "Scaling the e uncertainty ×0.1 reproduces S to 0.1%, so the reweight is a "
        "smooth local expansion — not chaotic extrapolation",
        ("Linearity does not validate the gradient's size. A wrong-but-smooth "
         "response would pass this test too.", 1),
    ],
)

d.table(
    "The pre-fit estimate, which the fit revised down",
    [
        ["η group", "A [σ]", "e [σ]", "M [σ]"],
        ["0  (η < −1.2)", "5.2", "13.1", "1.1"],
        ["1  (−1.2 to 0)", "5.5", "13.1", "1.3"],
        ["2  (0 to 1.2)", "5.9", "13.5", "1.4"],
        ["3  (η > 1.2)", "3.7", "8.4", "0.8"],
    ],
    lead="607 000 signal candidates expected: 40 of 22 327 data files, S/B 0.51.",
    note="From the histogram shape distance alone, before any fit: 304 candidates "
    "scaled as √N with a √(S/(S+B)) dilution. The fit gives 7σ on e rather than 13, "
    "the difference being the free per-cell normalisation and the template's own "
    "statistics. Not the A–e degeneracy, which turned out to be negligible.",
)

d.bullets(
    "The gen-match assertion, checked not assumed",
    [
        "Exclusive signal: 366 of 366 candidates fully three-leg matched, 12 files",
        "Other real-b decays: 100%. Combinatorial: 90.0%",
        "Structural, not lucky — the chain walk needs three matched legs to reach a "
        "common ancestor, so a partially matched signal candidate cannot exist",
        "A counted audit histogram, not a throw: a throw inside a Define kills the "
        "event loop, same reason as the shift guard",
    ],
)

d.bullets(
    "The throw guard earns its keep",
    [
        "calculateShiftedPtExact raises std::domain_error once the shifted energy "
        "falls below m_K — inside a Define that ends the job, hours in",
        "Replaced by a non-throwing variant that falls back to the massless "
        "linearisation for that (leg, variation) and counts it",
        "At pT > 1 GeV: zero fallbacks. At pT > 0.1 with e = 50 MeV: 20 on 12 files",
        "A local divergence from havyn's header, marked in the source, worth "
        "offering upstream",
    ],
)

d.table(
    "Variation size as a mass shift",
    [
        ["η group", "A [keV]", "e [keV]", "M [keV]"],
        ["0", "−69", "70", "3"],
        ["1", "−57", "103", "9"],
        ["2", "−33", "55", "−3"],
        ["3", "9", "1", "−13"],
    ],
    lead="Half the down-to-up span of the mean B mass.",
    note="Measured before the convention fix, so the forward groups are understated "
    "by up to cosh η. The down/up asymmetry of the exact shift is under 1.3 keV. The "
    "mean mass throws away the shape, which is why S is the figure of merit.",
)

d.bullets(
    "Why the splines cross-check is deferred, not skipped",
    [
        "The 2018 splines route does not use the muon tflite — it uses a "
        "kaon-trained one, kaon_response_studyanchorplusresidualpchipcells_v1.tflite",
        "So it is the only response model we have that saw kaons, and comparing "
        "against it is more than a numerical consistency check",
        "But it predates this selection and this NanoAOD, so a comparison now would "
        "confound three changes at once",
        "Plan: regenerate it, then cross-check. The histmaker rejects the splines "
        "backend with that explanation rather than running something misleading.",
    ],
)

d.table(
    "Configuration, as it now stands",
    [
        ["setting", "value", "why"],
        ["ONNX model", "nosmearing", "this channel applies no smearing"],
        ["bachelor pT", "> 1 GeV", "background, and the throw margin"],
        ["e convention", "pt", "matches the central-value corrector"],
        ["η groups", "4 of 48", "144 nuisances on 10⁴ MC candidates"],
        ["covariance", "off", "coarsening needs it off"],
        ["masses", "[m_K, 0, 0]", "explicit; one element is broadcast"],
        ["cond_pt_gen_min", "2.0 GeV", "matches the training floor"],
        ["muon_source", "443", "no hadron class exists"],
        ["mass window", "[5.05, 5.45]", "25% more S(e) than ±75 MeV"],
        ["mass bins", "5 MeV", "σ = 13.8 MeV, so 6.5 bins per FWHM"],
    ],
    note="48 η bins remains the end state: the joint mZ + J/ψ fit correlates by "
    "nuisance name on that basis, so anything coarser cannot be combined.",
)

d.bullets(
    "What is left",
    [
        "Full-sample pass: the match audit and the fallback counter travel with every "
        "run, so those numbers arrive with the first production job",
        "Template-tail truncation by the truth definition, re-measured here — a badly "
        "measured kaon fails the dR match and leaves the signal category (3.6%)",
        "Resolution variations, which need the smearing helper this channel does not "
        "run, and would flip the ONNX model to the smeared one",
        "Grep emanca's 21 rabbit-only commits for anything the ONNX path or "
        "d0_tensor.py needs, and cherry-pick rather than adopt by pointer",
    ],
)


# ---------------------------------------------------------------- background
d.bullets(
    "The background family, decided by measurement not by taste",
    [
        "Both component normalisations profiled out by lstsq at every trial, so only "
        "shape can win — the same guard as the attribution scan",
        "exp(a + b·x), what the model has now: χ² 43 047 over 80 bins",
        "exp(a + b·x + c·x²), one parameter more: 21 781, a factor 1.98",
        "Bernstein order 2, two parameters more: 21 311, a factor 2.02",
        "One parameter buys 98% of what two buy — and exp() of a real number needs "
        "no positivity constraint, so no boundary and no c² → 0 degeneracy",
    ],
)

d.table(
    "Why the basis is centred: free conditioning, same optimum",
    [
        ["basis", "corr(linear, quadratic)", "χ² at best fit"],
        ["x, x² on [0, 1]", "+0.968", "21 780.7"],
        ["u, u² on [−½, +½]", "−0.000", "21 780.7"],
    ],
    lead="Same functions, so the fit cannot tell them apart. The Hessian can.",
    widths=[0.34, 0.33, 0.33],
    note="AxisBernsteinModel is FIRST order, c₀(1−x) + c₁x — a linear background "
    "with no curvature. BKGMODEL=bernstein was never the curved alternative; "
    "AxisExpQuadModel is.",
)

d.table(
    "Moving the window is not the fix",
    [
        ["window", "bins", "exp, χ²/bin", "expquad, χ²/bin"],
        ["5.050–5.450", "80", "551.9", "282.9"],
        ["5.150–5.450", "60", "422.2", "353.2"],
        ["5.200–5.450", "50", "439.0", "419.0"],
    ],
    lead="Cutting at 5.15 makes χ² per bin worse, once the family is right.",
    widths=[0.28, 0.14, 0.29, 0.29],
    note="Cutting removes sideband that constrains the background and keeps the "
    "peak-region disagreement, which is where 46% of the χ² actually lives.",
)

d.table(
    "Where the disagreement actually is",
    [
        ["mass range", "share of χ²", "mean data/pred"],
        ["5.050–5.100", "28.2%", "0.855"],
        ["5.100–5.250", "7.8%", "0.99–1.07"],
        ["5.250–5.300", "45.9%", "1.119"],
        ["5.300–5.400", "0.5%", "1.00"],
        ["5.400–5.450", "17.6%", "0.883"],
    ],
    lead="Both edges over-predicted, the peak under: an arch, not a slope.",
    widths=[0.34, 0.33, 0.33],
)

d.bullets(
    "After the curvature, what is left is the signal's business",
    [
        "Best fit of the curved background: residual within ±4% across the whole window",
        "The largest remaining terms are +4.2% at 5.250–5.300 and −3.8% at 5.300–5.350",
        "That is an antisymmetric pair straddling m(B⁺) = 5.279 — a position and width "
        "residual, not a background one",
        "So the ordering background → resolution is not just cheaper, it is separable: "
        "fitting resolution first would have absorbed background curvature into it",
    ],
)

# ------------------------------------------------------------------- Havyn
d.bullets(
    "Havyn's D⁰ → Kπ deck settles two of our open questions",
    [
        "Energy loss for a massive track: he writes δk = −εk²/β, no cosh η",
        "That is our calculateQopUncPtConvention exactly, with 1/β = E/p — the pT "
        "convention, not the cosh η one, so openspec 2b.2 resolves in favour of "
        "what we already run",
        "ONNX reweights were not trained below 2 GeV; he floors the conditioning pT "
        "at 3 GeV for margin, we floor at 2.0, the boundary itself",
        "His A–e correlation is 0.81 (4-channel) and 0.90 (3-channel) — ours at 0.60 "
        "and 0.63 is the less correlated of the two",
    ],
)

d.table(
    "The pT floor is not a corner case for us",
    [
        ["sample", "median pT", "< 2 GeV", "< 3 GeV"],
        ["gen decays, all", "2.25 GeV", "44.4%", "62.8%"],
        ["gen decays in acceptance", "2.75 GeV", "31.0%", "53.8%"],
        ["selected, signal category", "3.25 GeV", "22.7%", "43.2%"],
    ],
    lead="Our 2.0 GeV floor freezes 22.7% of the signal; his 3.0 would freeze 43.2%.",
    widths=[0.40, 0.20, 0.20, 0.20],
    note="Neither choice is obviously right — inside the domain but near the edge is "
    "where a network is least reliable. The honest treatment is to float the floor "
    "between 2 and 3 GeV and take the A/e/M difference as a systematic.",
)

d.table(
    "What a third low-mass channel buys, from his Asimov",
    [
        ["impact on mZ [MeV]", "with D⁰", "without D⁰"],
        ["total", "1.93", "2.61"],
        ["magnetic field (A)", "1.28", "2.17"],
        ["energy loss (e)", "1.04", "2.01"],
        ["data stat", "1.03", "1.08"],
        ["misalignment (M)", "0.04", "0.06"],
    ],
    lead="One more low-momentum channel nearly halves the A and e impacts.",
    widths=[0.40, 0.30, 0.30],
    note="His caveat, quoted: status update, not final — low MC statistics and "
    "hacked ONNX reweights.",
)

# ------------------------------------------------------------- gen matching
d.bullets(
    "The PAT embedded gen ref: neither available nor better",
    [
        "It exists only on pat::Muon and pat::Electron — Bmm5's hadron constructor "
        "leaves gen_particle_ null, so the bachelor is ALWAYS matched by the cone "
        "fallback, in Bmm5 as in ours",
        "For muons it is MCMatcher at ΔR < 0.5 and |Δpt|/pt < 0.5 — 25× and 5× looser "
        "than Bmm5's own dr_match, matched in collection order",
        "Its own config comment: \"cuts are NOT tuned … old values from TQAF, "
        "january 2008\"",
        "We could not use it regardless: the AlCaReco carries "
        "reco::RecoChargedCandidate, which has no gen reference at all",
    ],
)

d.table(
    "Its two real advantages, measured on our sample",
    [
        ["what PAT adds", "would it change our answer?"],
        ["resolveAmbiguities (one-to-one)", "0 of 219 995 three-leg candidates"],
        ["  leg0=leg1, leg0=leg2, leg1=leg2", "0, 0, 0"],
        ["mcPdgId on leg0 (muon)", "27 = 0.012% (16 K, 11 π)"],
        ["mcPdgId on leg1 (muon)", "21 = 0.010% (14 K, 7 π)"],
        ["mcPdgId on the bachelor", "applied at the category stage instead"],
    ],
    lead="Both are things we lack. Neither costs us anything.",
    widths=[0.48, 0.52],
    note="The species requirement on the bachelor is not missing but moved: the "
    "category demands leg2's gen pdgId be a kaon, so a pion bachelor is LABELLED "
    "GEN_BU_PIX rather than silently dropped to unmatched. That is the better place "
    "for it.",
)


d.table(
    "Why the data arm fails: a ladder in mis-specification",
    [
        ["what the \"data\" is", "how wrong the model is", "edmval"],
        ["Asimov", "not at all, and at the truth", "3.6e−29"],
        ["toy (model + Poisson)", "not at all, plus fluctuation", "1e−17"],
        ["mcTotal pseudodata", "arched combinatorial vs monotonic bkg", "986"],
        ["real data", "that, plus the 17% width deficit", "72 668"],
    ],
    lead="Toys converge, so the parameterisation is sound.",
    widths=[0.30, 0.45, 0.25],
    note="Threshold is 1e−5. A toy is generated from the model itself, so the model "
    "is right by construction; mcTotal is NOT that, since it carries the simulation's "
    "arched combinatorial and so mis-specifies the analytic background too.",
)

d.bullets(
    "One reading I had wrong: the scipy message is a red herring",
    [
        "\"A bad approximation caused failure to predict improvement\" appears on all "
        "three CONVERGED toys, at edmval 1e−17",
        "trust-krylov emits it whenever the trust region shrinks at its stopping "
        "point, converged or not",
        "edmval is the only quantity that says whether a fit converged — which is what "
        "the new gate tests, so the gate is right even though my reading was not",
        "And edmval is not comparable ACROSS models: it is a distance scaled by the "
        "Hessian, so a more flexible model can stop further out at equal fit quality",
    ],
)

d.table(
    "The curvature parameter is free",
    [
        ["quantity", "exp", "expquad"],
        ["scaleA0 uncertainty", "0.856", "0.856"],
        ["scalee0 uncertainty", "0.386", "0.386"],
        ["scaleM0 uncertainty", "0.968", "0.968"],
        ["A–e correlation, η0", "+0.5911", "+0.5882"],
        ["A–e correlation, η1", "+0.6156", "+0.6153"],
        ["max |corr| curv with any scale", "—", "0.0588"],
    ],
    lead="Converged Asimov, both models. Mean ratio 1.000, worst 1.002.",
    widths=[0.46, 0.27, 0.27],
    note="The EXISTING slope is more entangled with e than the new curvature is: "
    "slope–scalee0 = −0.1131 against curv–scalee0 = −0.0588. And A–e = 0.59–0.62 "
    "now comes from a converged fit, against Havyn's 0.81 and 0.90.",
)

d.table(
    "Also ruled out, one flag at a time",
    [
        ["change, one at a time", "outcome"],
        ["--minimizerMaxiter 50 000", "terminates itself, Cholesky fails"],
        ["--minimizerMethod trust-exact", "edmval 410 524, worse"],
        ["--noBinByBinStat", "Cholesky fails"],
        ["trust-exact + noBinByBinStat", "edmval −1.5e15, negative"],
        ["L-BFGS-B", "no edmval, fails"],
        ["expquad instead of exp", "Asimov 2.2e-29, data still fails"],
    ],
    lead="A negative edmval is an indefinite Hessian: a saddle, not a hard-to-reach "
    "minimum.",
    widths=[0.45, 0.55],
    note="Also ruled out: a cell wanting zero background, which would send its lnAmpl "
    "to −∞ with vanishing gradient and would have explained the Asimov/data split. "
    "The tensor log says 0 of 16 cells hit the positivity floor.",
)

d.bullets(
    "So the next step IS the background model after all",
    [
        "Toys converge at 1e−17, so nothing is wrong with the parameterisation or the "
        "minimiser — the fit fails because the model does not describe the data",
        "Curvature first: it is free in uncertainty, converges on Asimov, and removes "
        "the larger of the two known shape errors (χ² factor 1.98 against 1.16)",
        "Then the resolution nuisance, which is the remaining data/MC difference and "
        "is absent from both toys and mcTotal",
        "The MC combinatorial template is still worth building — as the nominal shape "
        "under a free tilt, and as a stricter closure test than mcTotal gives",
        "Judge each step by edmval, not by how the postfit plot looks",
    ],
)


# ------------------------------------------------- stage A: the analytic model
d.bullets(
    "Stage A: a general-order background, order measured",
    [
        "AxisExpPolyModel: exp(lnAmpl + Σₖ cₖ Pₖ), any order, orthogonal basis",
        "Order 2 buys χ² factor 2.00; order 3 adds 0.03, order 4 adds 0.04, order 5 "
        "adds nothing — so order 2 is the model and 3/4 are cross-checks",
        "Free cross-check: order 2 Asimov gives edmval 2.2339e−29, identical to the "
        "expquad Asimov to five figures — independently the same fit",
        "All three orders converge on Asimov. ALL THREE FAIL ON DATA, at nit 263, "
        "311 and 473 with an indefinite Hessian",
        "More background freedom does not rescue the data fit — as the toy ladder "
        "predicted, the remaining error is the resolution",
    ],
)

d.table(
    "The basis must be orthogonal to the constant",
    [
        ["basis on the mass grid", "max overlap", "consequence"],
        ["x, x² on [0,1]", "0.968", "hopeless conditioning"],
        ["u, u² centred", "⟨u²⟩ = 1/12", "constant explains 55.6% of u²"],
        ["nominal Legendre, centre-anchored", "1.3e−2", "still overlaps lnAmpl"],
        ["Gram-Schmidt on the real grid", "**5e−17**", "orthogonal by construction"],
    ],
    lead="lnAmpl is a free constant per cell, so overlap with it stalls the fit.",
    widths=[0.40, 0.22, 0.38],
    note="Anchor on the axis EDGES, not the first and last bin centre: centres at "
    "exactly ±1 break the discrete orthogonality. The existing exp and expquad "
    "models sit at 1e−6 from float32 basis constants; this one is at machine "
    "precision.",
)

# -------------------------------------------------- stage B: the MC templates
d.table(
    "The \"peaking\" backgrounds do not peak",
    [
        ["process", "window", "peak ±30 MeV", "peak fraction"],
        ["DATA", "6 206 572", "1 918 758", "0.309"],
        ["signal", "1 408 371", "1 162 085", "**0.825**"],
        ["MC combinatorial", "3 728 169", "677 791", "0.182"],
        ["MC \"peaking\" b", "297 826", "11 651", "**0.039**"],
    ],
    lead="They lose a pion, so they sit in the low sideband.",
    widths=[0.34, 0.22, 0.22, 0.22],
    note="Peak fraction 0.039 against the signal's 0.825, so including them cannot "
    "bias A/e/M through the peak — which was the reason to worry about them. It also "
    "explains why the old \"635 candidates in the peak\" was so small.",
)

d.bullets(
    "Three variants; the analytic component always stays",
    [
        "The simulation has no prompt and no non-b fake background, so it cannot "
        "describe the data alone — the analytic component is not optional",
        "all: every non-signal category in one template. comb: combinatorial only. "
        "split: combinatorial and peaking-b independently normalised",
        "Templates seeded at the simulation's OWN prediction, so a fitted norm of 1 "
        "means the simulation was right — the fit measures composition, not a scale",
        "Warning already visible at the seed: all and split push 5 of 16 cells onto "
        "the analytic positivity floor, comb pushes none",
        "That is the two-normalisation tension, and the fitted norms are the test",
    ],
)

# ---------------------------------------------- stage C: resolution, from basics
d.bullets(
    "From basics: what a resolution nuisance is for",
    [
        "What we see is the true lineshape convolved with the detector response: "
        "m_obs = m_true ⊗ R(σ)",
        "Measured: data peak σ = 17.29 ± 0.04 MeV, template σ = 14.81 ± 0.11, "
        "ratio 1.168 ± 0.009 — the simulated response is 17% too narrow",
        "A, e and M SHIFT the peak; the normalisations SCALE it. Nothing in the "
        "model WIDENS it",
        "A shift is an odd distortion of the lineshape and a width error is an even "
        "one, so the width error has no parameter to land on",
        "It lands on the scale parameters instead — which is what the multi-σ pulls "
        "and the indefinite Hessian are",
    ],
)

d.bullets(
    "Stage C: how the nuisance is built",
    [
        "One alternative template: the signal convolved with a Gaussian of "
        "δ = √(17.29² − 14.81²) = 8.9 MeV, so the fitted θ reads as multiples of the "
        "measured deficit",
        "Integral preserved per cell to 1.000000000, so it is a pure SHAPE variation "
        "and cannot overlap the free signal normalisation",
        "Kernel renormalised bin by bin at the window edges — zero-padding would "
        "quietly leak yield out of the window",
        "rabbit interpolates multiplicatively, (up/nominal)^θ: θ=0 nominal, θ=1 the "
        "smeared template, θ<0 narrower by log-extrapolation",
        "Free first (no prior width is defensible for a measured 17% error), then "
        "constrained for the uncertainty budget",
    ],
)

d.bullets(
    "One construction I got wrong first: the additive mirror",
    [
        "I first built the down variation as 2×nominal − up, the usual symmetric "
        "linear construction",
        "It went NEGATIVE in 274 of 1280 bins: smearing by 1.78 bins moves more yield "
        "into the tails than twice the nominal covers there",
        "A single template with rabbit's mirror=True is positive for every θ, because "
        "the interpolation is multiplicative",
        "And log-symmetric is the right shape anyway — the width discrepancy is a "
        "RATIO (1.168), not a difference",
    ],
)


# ------------------------------------------------------ 2018 vs now: answered
d.table(
    "Why 2018 worked: purity, and it was already on disk",
    [
        ["configuration", "purity", "data σ", "MCtot σ", "MCsig σ", "d/tot", "d/sig"],
        ["2018, full selection", "n/a", "15.25", "14.87", "14.87", "1.025", "1.025"],
        ["2018, no BDT / softMVA", "n/a", "15.25", "14.87", "14.87", "1.025", "1.025"],
        ["CVH baseline", "0.214", "23.62", "22.25", "15.12", "1.062", "**1.562**"],
        ["CVH, kaon pT > 1", "0.579", "17.50", "16.12", "15.75", "1.085", "**1.111**"],
    ],
    lead="In 2018 the simulated TOTAL and the simulated SIGNAL have the same width.",
    widths=[0.28, 0.11, 0.11, 0.13, 0.13, 0.12, 0.12],
    note="At that purity they are the same object, so signal-plus-monotonic is "
    "nearly exact. At purity 0.214 the total is 47% wider than the signal, so the "
    "model must describe the background shape — which a monotonic form cannot.",
)

d.bullets(
    "And it revises the resolution claim I have been quoting",
    [
        "compare_2018_width.py warns in its own docstring: comparing a spectrum "
        "with background to a template without one gave a spurious factor 1.71",
        "cvhnano_width_check.py does exactly that — it fits the signal category "
        "alone against data containing everything. The 17% is that comparison",
        "Against the simulated TOTAL: 1.025 (2018), 1.062 (CVH), 1.085 (kaon pT>1) "
        "— six to nine percent, not seventeen",
        "d/sig tracks purity, 1.562 at 0.214 and 1.111 at 0.579 — the signature of "
        "contamination, not of resolution",
        "The BDT/softMVA hypothesis is dead: 2018 with and without gives identical "
        "widths, only the yield changes",
    ],
)

d.bullets(
    "So 8.9 MeV is a starting scale, not a measurement",
    [
        "The deficit is bracketed between ~6% and ~17%, and the spread IS the "
        "background-shape problem one level down",
        "A background flexible enough to absorb the combinatorial says composition "
        "contributes ±0.5 MeV; a linear one over ±74 MeV says the total is 47% wider",
        "Which is why the nuisance is floated freely rather than pinned to 8.9 MeV",
        "The fits agree: the resolution nuisance moved the data NLL by 0.2% "
        "(316 816 → 316 154), an MC background template moved it by a factor of 6",
        "Reorders the work — the background shape is the whole story so far, and "
        "resolution is a correction on top of it",
    ],
)


d.table(
    "What each stage costs, on converged Asimov fits",
    [
        ["configuration", "scaleA0", "scalee0", "scaleM0", "mean", "worst"],
        ["order 2, analytic only", "0.856", "0.386", "0.968", "1.000", "—"],
        ["+ MC combinatorial", "0.877", "0.548", "0.975", "1.153", "**1.419**"],
        ["+ MC all non-signal", "0.882", "0.551", "0.977", "1.157", "1.425"],
        ["+ MC comb & peak split", "0.882", "0.551", "0.977", "1.158", "1.426"],
        ["+ resolution (free or constr.)", "0.856", "0.386", "0.968", "1.0001", "1.0003"],
    ],
    lead="The resolution nuisance is free; the MC template costs e 42%.",
    widths=[0.34, 0.13, 0.13, 0.13, 0.13, 0.14],
    note="Opposite to what the NLL alone suggests. And the MC template is "
    "degenerate with the analytic curvature: the order-2 coefficient's "
    "uncertainty inflates 23–27× (0.0016 → 0.0364, 0.0440) when a template is added.",
)

d.bullets(
    "Asimov and toys converge everywhere; data never does",
    [
        "Toys from the analytic model: edmval 1e−17. Toys from the MC-template "
        "model: 1e−13 to 1e−16. Both with zero Cholesky failures",
        "Removing the A/e/M nuisances entirely (SCALEVAR=0) changes nothing — the "
        "data fit fails identically, so it is NOT the scale parameters",
        "Nor the minimiser, the iteration budget, the bin-by-bin nuisances, the "
        "polynomial order, the MC template, or the resolution nuisance",
        "The parameterisation is sound everywhere it can be tested; only real data "
        "fails, and it fails the same way in every configuration",
    ],
)

d.bullets(
    "A process failure: the submodule ate the model",
    [
        "rabbit is a git submodule; a reset to havyn/upstream-main-havyn removed "
        "AxisExpPolyModel and its registry entry mid-session",
        "Six fits had already used it. The next run failed with \"Class "
        "AxisExpPolyModel not found\" — and three queued fits reported \"done\" in "
        "two minutes having done nothing",
        "The class now lives OUTSIDE the submodule with an idempotent installer; "
        "run it after any submodule update",
        "A run that produces output has not necessarily done the work — the same "
        "lesson as edmval, one level up",
    ],
)


d.table(
    "The MC template is statistics-limited, and badly",
    [
        ["", "total", "per bin", "stat per bin"],
        ["DATA", "6 206 572", "4 849", "**1.44%**"],
        ["signal template", "1 408 371", "1 100", "**13.26%**"],
        ["MC background template", "4 025 995", "3 145", "**8.09%**"],
        ["analytic background", "822 252", "642", "0 by construction"],
    ],
    lead="Each simulated background event stands for ~19 data events.",
    widths=[0.34, 0.22, 0.20, 0.24],
    note="The MC background template is 5.5× noisier than the data it describes, "
    "and the SIGNAL template is noisier still. 1280 bin-by-bin nuisances then "
    "absorb 8% per bin on the background and 13% on the signal — which is the "
    "mechanism behind the 42% inflation of e. Only more simulation fixes this.",
)

d.table(
    "Higher purity does not rescue the data fit either",
    [
        ["window", "nit", "outcome"],
        ["5.23–5.33 (±50 MeV)", "8000", "stalls, last step exactly 0"],
        ["5.20–5.36 (±80 MeV)", "8000", "stalls, last step exactly 0"],
        ["5.15–5.41 (±130 MeV)", "149", "Cholesky"],
        ["5.05–5.45 (±200 MeV)", "263", "Cholesky"],
    ],
    lead="The failure mode changes; the failure does not.",
    widths=[0.36, 0.16, 0.48],
    note="So purity explains why 2018's MODEL was adequate — it does not explain "
    "the non-convergence. Those are two different questions and only the first is "
    "now answered.",
)


d.table(
    "MC statistics is 81.5% of the A/e/M uncertainty",
    [
        ["parameter", "with bbb", "no bbb", "MC stat", "MC share"],
        ["scaleA0–3", "0.856", "0.452", "0.727", "**85–88%**"],
        ["scalee0–3", "0.386", "0.202", "0.329", "**85–90%**"],
        ["scaleM0–3", "0.968", "0.718", "0.649", "67–74%"],
    ],
    lead="Two converged Asimov fits differing only by --noBinByBinStat.",
    widths=[0.24, 0.19, 0.19, 0.19, 0.19],
    note="The analytic background carries zero variance by construction, so this "
    "is the SIGNAL template alone. It supersedes the per-bin rule of thumb (85× for "
    "noise parity), which overstates the cost because the fit aggregates 1280 bins.",
)

d.table(
    "6× captures 79% of the available gain",
    [
        ["more signal MC", "MC share", "total unc", "reduction"],
        ["1×", "81.5%", "1.000", "—"],
        ["2×", "70.5%", "0.817", "18%"],
        ["4×", "57.5%", "0.708", "29%"],
        ["**6×**", "**49.8%**", "**0.668**", "**33%**"],
        ["25×", "27.1%", "0.602", "40%"],
        ["100×", "13.9%", "0.585", "42%"],
    ],
    lead="The data-statistical floor is 57.9% of today, so 42% is the whole prize.",
    widths=[0.28, 0.24, 0.24, 0.24],
    note="Everything beyond ~25× is chasing the last two percent.",
)

d.bullets(
    "And 6× is cheap, if the production is exclusive",
    [
        "Today: 6.77M nano events → 714k contain a gen B⁺→J/ψK → 72.8k selected "
        "signal, so 1.08% of nano events yield one",
        "An exclusive B⁺→J/ψK production puts the decay in every event: the same "
        "chain then gives 10.19%, a measured 9.48× gain per event",
        "So 6× the template statistics needs 4.29M exclusive events — 0.63× the "
        "full simulation already run — against 40.6M inclusive-style",
        "BUT the background template needs the INCLUSIVE sample — 70.4% of selected "
        "candidates are NoBAncestor, which an exclusive sample does not contain",
        "So the background strategy decides the production, and the analytic route "
        "is 9.5× cheaper",
    ],
)

d.bullets(
    "What more MC does not fix",
    [
        "A template whose SHAPE is wrong: the 6–17% width discrepancy is a modelling "
        "error, and more events reproduce it more precisely, not less",
        "The forward-η efficiency of 0.63× central — that is acceptance, and it is "
        "the cells with the longest cosh η lever arm on e",
        "A–e = 0.59–0.62 — Havyn sees 0.81–0.90 on the same kind of fit, so that is "
        "physics, not a defect",
        "The non-convergence: excluded so far are the minimiser, iterations, "
        "bin-by-bin nuisances, polynomial order, all MC variants, the resolution "
        "nuisance, the scale nuisances entirely, and window widths ±50 to ±200 MeV",
    ],
)


d.table(
    "The failure is model error against data precision",
    [
        ["dataScale", "nit", "edmval", "Hessian"],
        ["1.0 (full)", "263", "none", "indefinite"],
        ["0.3", "128", "none", "indefinite"],
        ["0.03", "201", "none", "indefinite"],
        ["**0.003**", "1536", "**44.4**", "**positive-definite**"],
    ],
    lead="Same data shape, weaker Poisson weights — only the significance changes.",
    widths=[0.26, 0.16, 0.24, 0.34],
    note="At 0.3% of the statistics the Hessian turns positive-definite and edmval "
    "is a factor 1637 better than the full-statistics 72 668. So more MC does NOT "
    "fix the non-convergence — only a better model shape does. And the failure is "
    "a symptom of having 6.2M events, which is why 2018 never saw it.",
)


# ---------------------------------------------------------------- 2026-08-28
d.section("Update, 28 August")

d.table(
    "The event accounting, all three numbers that get confused",
    [
        ["", "data", "inclusive MC"],
        ["Pythia events generated", "—", "4.03×10⁸"],
        ["nano events", "87,341,047", "6,770,705"],
        ["…containing a gen B⁺→J/ψK⁺", "—", "714,190"],
        ["selected candidates", "6,246,614", "282,753"],
        ["**signal**", "**1,363,434**", "**72,800**"],
        ["purity", "0.218", "0.257"],
    ],
    lead="Full 2016 sample, 16.89 fb⁻¹, fit window 5.05–5.45.",
    widths=[0.40, 0.30, 0.30],
    note="6.77M is the nano INPUT, not a yield — earlier slides also quoted 62k "
    "selected signal from a looser accounting. The gen-signal count is 72,916 over "
    "the full 5.0–5.5 and 72,800 in the fit window. And the ratio that matters is "
    "18.7 on signal, not the 22.1 you get from all candidates.",
)

d.bullets(
    "The signal MC is worth 0.90 fb⁻¹ against 16.89 fb⁻¹",
    [
        "N_sig(data) = 1,363,434 is measured, not assumed, and stable to ±4% "
        "across every background variant",
        "Cross-check: 500 μb × 16.89 fb⁻¹ × 2f_u × BRs = 4.19×10⁸ produced, so "
        "ε_total = 0.33% — right for a displaced dimuon trigger",
        "Three routes, one answer: 18.7× fewer signal events, MC stat 80.8% of "
        "the A/e/M uncertainty, template 2.0× noisier per bin",
        "6× costs 4.29M exclusive full-sim events — less than the 6.77M already run",
    ],
)

d.bullets_figure(
    "The model is wrong at 1%; the data has 0.36% per bin",
    [
        "One cell, 6.2M events, 80 bins: χ²/ndf = 53.6 — no A/e/M, no 288 cells",
        "S-shape across the peak: +8.9σ at 5.2775, then −11.8 and −13.7 above",
        "Wave in the low shoulder: +4% at 5.08–5.13, then −5% at 5.15–5.20",
    ],
    "inclusive_datafit.png",
)

d.table(
    "Two fixes, factor 5.6 — and they beat pure flexibility 2.2×",
    [
        ["model", "npar", "χ²/ndf", "N_sig", "σ [MeV]"],
        ["exp·poly2, no resolution", "4", "79.9", "1.314M", "—"],
        ["+ Gaussian resolution", "5", "53.6", "1.363M", "6.26"],
        ["+ asymmetric instead", "6", "44.2", "1.364M", "5.66"],
        ["Gaussian + B⁰/J/ψπ/J/ψKX", "8", "26.4", "1.330M", "5.57"],
        ["**asymmetric + those**", "**9**", "**10.1**", "1.328M", "4.77"],
        ["exp·poly6 + Gaussian res", "9", "22.4", "1.280M", "4.55"],
    ],
    lead="Asymmetric = a Gaussian whose low-side width is 1.27× the high-side one.",
    widths=[0.36, 0.12, 0.16, 0.18, 0.18],
    note="The last row is the control: same 9 parameters, 2.2× worse. Resolution is "
    "worth Δχ² 2839 and the peaking templates a further 2585 — each dwarfs any "
    "background order (2→3 buys 623, 3→4 buys 3).",
)

d.figure(
    "The low shoulder is a B⁰ edge, not combinatorial",
    "inclusive_datafit_physics.png",
    lead="χ²/ndf 10.1. B⁰ and J/ψKX die at 5.15 — an edge carrying 52% of the χ².",
)

d.bullets(
    "Two corrections to what this deck said earlier",
    [
        "\"Don't add peaking channels\" is WITHDRAWN. The peak-fraction measurement "
        "(0.039 vs signal's 0.825) was right; the inference was not",
        "Whether a component peaks at the B mass is the wrong question — whether it "
        "puts structure where a smooth background cannot follow is the right one",
        "RESVAR must come down from 0.0089 to ~0.005: the 8.9 MeV is a "
        "data-vs-MC-signal width comparison, the one the 2018 script warns against",
        "Fitting the smear with the background free gives 6.26 MeV Gaussian, or "
        "4.77 MeV with r = 1.27 once the peaking templates stop stealing the tail",
    ],
)

d.table(
    "Template statistics, per bin, against the data's own 0.36%",
    [
        ["component", "MC evts", "→ data", "MC/bin", "of total", "vs data"],
        ["signal", "72,800", "1.363M", "899", "0.73%", "**2.0×**"],
        ["combinatorial (NoB)", "194,034", "3.634M", "2,395", "1.19%", "**3.3×**"],
        ["B⁰", "7,687", "144k", "95", "0.24%", "0.7×"],
        ["B⁺→J/ψKX", "3,829", "72k", "47", "0.17%", "0.5×"],
        ["B⁺→J/ψπX", "3,600", "67k", "44", "0.16%", "0.4×"],
        ["Bs", "714", "13k", "9", "0.07%", "0.2×"],
    ],
    lead="Fractional template error 1/√N per bin, weighted by that component's "
    "share of the data.",
    widths=[0.28, 0.15, 0.14, 0.13, 0.15, 0.15],
    note="Per CELL the simulation holds 27 B⁰, 13 J/ψKX and 12 J/ψπX events over 81 "
    "mass bins — unusable. The peaking templates are a second reason for more MC, "
    "and an exclusive B⁺ production does not serve it.",
)


d.table(
    "Matched to 2018's acceptance, the CVH fit is 5% BETTER",
    [
        ["acceptance", "2016 σ", "2018 σ", "ratio"],
        ["CVH selection, |η| < 2.4, pT > 1", "17.20 MeV", "16.30 MeV", "1.055"],
        ["**2018 acceptance, |η| < 1.4, 1 < pT < 8**", "**15.50 MeV**", "16.30 MeV", "**0.951**"],
    ],
    lead="Signal only: 2016 is the gen-matched slice, 2018 the exclusive BuToJpsiK.",
    widths=[0.46, 0.19, 0.19, 0.16],
    note="Like-for-like — the joint maker applies its own J/ψ constraint, verified "
    "in runCvhBplusJpsiK.py, so both are J/ψ-constrained masses. The apparent "
    "degradation is entirely the forward region 2018 never had: its kaon η axis "
    "stops at 1.4 and it capped the kaon at 8 GeV.",
)

d.figure(
    "Narrower core, smaller tails, inside the same acceptance",
    "cmp1618_mass_2018acc.png",
    lead="Ratio 1.3–1.4 in both tails and 0.93 at the peak — 2018 has more tail.",
)

d.figure(
    "2018 starts at 1.2 GeV and is 1.3–1.4× softer up to 2.2",
    "cmp1618_kaonpt_cvh.png",
    lead="Any low-pT comparison between the two analyses has to carry this.",
)

d.bullets(
    "Yes, we use 1–2 GeV kaons — 22.7% of the signal",
    [
        "11.4% in 1.0–1.5 and 11.2% in 1.5–2.0, against a response network whose "
        "training sample floors gen pT at 2.0 GeV",
        "Nothing is pathological: zero non-finite and zero negative in 23,040,000 "
        "variation entries, and the shifts are smooth through the boundary",
        "The physics is right too — A is flat at −0.052 to −0.058 MeV at every pT, "
        "and e grows 0.071 → 0.352 MeV as pT falls, as an energy-loss term must",
        "But smooth and finite is not the same as correct, and the next slide is why",
    ],
)

d.figure(
    "Below the floor, the clamp flips the sign of e",
    "pt_floor_variations.png",
    lead="Same run, only the clamp differs. Above 2 GeV they agree to 0.0000 MeV.",
)

d.table(
    "The clamp is worth up to 132% of the variation",
    [
        ["pT band", "param", "clamped", "extrapolated", "difference"],
        ["1.0–1.5", "**e (unc1)**", "**+0.352**", "**−0.111**", "**−0.463**"],
        ["1.0–1.5", "e (unc4)", "+0.236", "+0.062", "−0.174"],
        ["1.0–1.5", "A (unc0)", "−0.052", "+0.016", "+0.068"],
        ["1.5–2.0", "e (unc1)", "+0.272", "+0.158", "−0.114"],
        ["≥ 2.0", "all", "identical", "identical", "0.0000"],
    ],
    lead="Half up–down peak shift, MeV. unc1/unc4 are e, from the packing order.",
    widths=[0.17, 0.21, 0.20, 0.22, 0.20],
    note="Neither is validated: clamping asserts a flat residual below 2 GeV, "
    "extrapolating runs the network outside its domain. Take the difference as a "
    "shape nuisance on e — both outputs already exist. Cutting at 2 GeV instead "
    "loses 22.7% of the signal AND halves e's lever arm.",
)

d.table(
    "Which peaking shapes vary across cells",
    [
        ["component", "χ²/ndf, split at kaon pT 1.6", "⟨m⟩ shift", "treatment"],
        ["B⁺→J/ψKX", "1.23", "−0.4 MeV", "one shared shape"],
        ["B⁰", "3.38", "−34 MeV", "per pT band"],
        ["B⁺→J/ψπX", "**9.53**", "**+115 MeV**", "per pT band"],
        ["combinatorial", "2.67 (η split 5.12)", "−3.7 MeV", "analytic, slope per η"],
    ],
    lead="Shape χ² between kaon-pT slices, 80 bins, normalisations divided out.",
    widths=[0.24, 0.30, 0.20, 0.26],
    note="The 115 MeV is the pion-as-kaon mass misassignment, ∝ (m_K²−m_π²)/2p — "
    "~75 MeV at p = 1.5 GeV against ~22 MeV at 5. It is physics, so that component "
    "can eventually be built from the signal template with a computed shift.",
)

d.table(
    "The production ask, split by what each sample is for",
    [
        ["purpose", "route", "nano events", "vs sim already run"],
        ["signal template, 6×", "exclusive B⁺→J/ψK⁺", "4.29M", "0.63×"],
        ["peaking shapes, 4×", "4× more inclusive", "27M", "4.0× (1.6×10⁹ Pythia)"],
        ["**peaking shapes, 4×**", "**exclusive J/ψπ⁺ + J/ψK*⁰**", "**~0.4M**", "**0.06×**"],
    ],
    lead="Two different problems; the exclusive route is ~68× cheaper for the second.",
    widths=[0.26, 0.30, 0.20, 0.24],
    note="Caveat: the inclusive B⁰ category is every B⁰→J/ψX, so an exclusive "
    "sample covers only the dominant mode of a mixture we cannot decompose "
    "(gap G1). Exclusive samples for the SHAPES, the existing inclusive one for "
    "the mixture composition — a normalisation question it already answers.",
)

print(d.save("260825_cvhnano_aem_variations.pptx"))
