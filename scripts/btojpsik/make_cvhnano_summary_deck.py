#!/usr/bin/env python3
"""Consolidated summary of the first look at the CVH NanoAOD B+ -> J/psi K+ data.

Condenses four working decks into one meeting deck:
  260812_btojpsik_cvhnano_first_look     -- the histmaker, the peak, the gap record
  260813_btojpsik_cvhnano_diagnostics    -- disagreements chased to causes, statistics
  260817_cvhnano_tensor_first_fits       -- the tensor, the fits, the background finding
  260818_cvhnano_truth_and_prefits       -- exclusive truth, composition, the kaon ID

Every number here is carried over from those decks; nothing is recomputed. Figures
are extracted from the source decks into assets_cvhnano_summary/.
"""
import sys

sys.path.insert(0, "/home/submit/pmlugato/.claude/skills/slides")
from deck import Deck  # noqa: E402

d = Deck(assets="/home/submit/pmlugato/public_html/slides/assets_cvhnano_summary")

# ----------------------------------------------------------------- 1
d.title_slide(
    "First look at the CVH NanoAOD B⁺ → J/ψ K⁺",
    "The chain works end to end. Purity is the problem, and the kaon track is the handle.",
    "Pietro Lugato  ·  19 August 2026  ·  mZ / B → J/ψ K",
)

# ----------------------------------------------------------------- 2
d.bullets(
    "Where this stands",
    [
        "Our own AlCaReco → CVH refit → NanoAOD → histmaker → tensor → fit runs end to end",
        "The B⁺ peak is there, at the right mass, and the three fit arms agree",
        "Exclusive B⁺ → J/ψ K⁺ truth is now definable, so every simulated component is identified",
        "The fit closes on Asimov and fails on data: the analytic background shape is wrong",
        "2018's advantage was purity, not resolution — and two of its cuts we cannot reproduce",
    ],
)

# ----------------------------------------------------------------- 3
d.bullets_figure(
    "The chain works: 18.5M data candidates, one clean peak",
    [
        "Full samples: 22 327 data files, 18.46M selected; 2 520 simulated, 887 100 selected",
        "Peak at 5.2798 GeV in data, 5.2800 in simulation — 0.2 MeV apart, both within 0.5 MeV of PDG",
        "Event loop 404 s on data, 152 s on simulation; plots from the central stack plotter",
        ("Figure is the earlier 400-file pass; its lumped 'combinatorial' is superseded overleaf", 1),
    ],
    "peak_arms.png",
)

# ----------------------------------------------------------------- 4
d.table(
    "Re-running the NanoAOD closed eight gaps",
    [
        ["gap", "what arrived", "evidence"],
        ["G1", "Gen_genPartIdxMother", "0.135% at −1, i.e. record roots only"],
        ["G2", "per-leg gen match", "kaon reco/gen pT 1.00001, rms 0.018"],
        ["G4", "Pileup_nTrueInt", "mean 21.82 sim vs 24.31 data"],
        ["G6", "track errors, hit detail", "finite and positive, 100% of tracks"],
        ["G11", "unconstrained m(μμ)", "median 3.09385 vs PDG 3.09690"],
        ["G13", "PV-referenced d0", "φ modulation 0.187 → 0.002 cm"],
    ],
    lead="No AlCaReco reproduction: everything came from files already on disk.",
    widths=[0.10, 0.38, 0.52],
    note="G7 closes as NOT A GAP: candidate muons are already 100% global-and-tracker, and "
    "softMuonMva cannot be evaluated at all — innerTrack() is unavailable for 100% of persisted muons.",
)

# ----------------------------------------------------------------- 5
d.table(
    "There is no peaking background",
    [
        ["category", "candidates", "% of peak window"],
        ["B⁺ → J/ψ K⁺", "72 974", "21.43%"],
        ["B⁺ → J/ψ K⁺X", "3 899", "0.00%"],
        ["B⁺ → J/ψ π⁺X", "6 149", "0.40%"],
        ["B⁰, B⁰ₛ, other b hadron", "11 508", "0.39%"],
        ["no b ancestor — combinatorial", "769 696", "77.77%"],
    ],
    lead="Joint-arm mass, gated, window 5.24 < m < 5.32. Full simulation.",
    widths=[0.46, 0.27, 0.27],
    note="Every real-b background together is 0.79% of the window, 3.7% of the signal, and each "
    "peaks 100–200 MeV BELOW the B mass. So no peaking template is needed — signal plus ONE background component is the right model. What that component must look like is the next slide.",
)

# ----------------------------------------------------------------- 6
d.figure(
    "The exponential misses by 30%, and the ratio says how",
    "prefit_jointCvh_mass.png",
    lead="Full statistics. The background peaks 3 MeV from the signal, so the residual double-humps.",
)

# ----------------------------------------------------------------- 7
d.table(
    "The background model is what breaks the fit",
    [
        ["background", "χ²/ndf", "fitted signal"],
        ["true simulated non-signal shape", "5.62", "132 391"],
        ["exponential", "745.79", "371 918"],
        ["quadratic", "85.47", "253 439"],
        ["cubic", "85.22", "—"],
    ],
    lead="54 bins over 5.15 < m < 5.42, 1 452 352 data candidates.",
    widths=[0.50, 0.25, 0.25],
    note="The exponential inflates the yield 2.8×. A cubic is no better than the quadratic, so this "
    "is not one more polynomial order away: rabbit's models are monotonic and the background has a "
    "maximum. The 288-cell fit closes on Asimov and fails only on data.",
)

# ----------------------------------------------------------------- 8
d.figure_row(
    "Mass, without and with kaon pT > 1 GeV",
    ["stack_residual_combofirst_truthv2.png", "stack_residual_combofirst_kpt1.png"],
    lead="Combinatorial at the bottom, log scale. Only the bachelor pT floor changes, 0.15 → 1 GeV.",
)

# ----------------------------------------------------------------- 9
d.table(
    "That cut buys almost everything — and 2018 had it",
    [
        ["", "baseline", "kaon pT > 1 GeV"],
        ["signal, % of peak window", "21.43%", "57.87%"],
        ["signal candidates kept", "65 992", "63 893"],
        ["exponential χ²/ndf", "745.8", "41.2"],
        ["quadratic χ²/ndf", "85.5", "8.5"],
        ["exponential yield bias", "2.8×", "1.11×"],
    ],
    lead="96.8% of the signal kept for 2.7× the purity.",
    widths=[0.46, 0.27, 0.27],
    note="Applied at PRODUCTION in 2018; our nano keeps kaons to 0.15 GeV. So every selection "
    "comparison before this one was blind to the single biggest difference between the samples.",
)

# ----------------------------------------------------------------- 10
d.bullets(
    "Why 2018 was immune, and the two cuts we cannot reproduce",
    [
        "Same window, same measure: the 2018 sample is 80% signal, ours is 12% — a factor seven",
        "At 20% background a wrong background shape barely moves the yield; at 88% it doubles it",
        "The BMM BDT and softMuonMVA remove 64% of what survives everything else in 2018",
        "But they do not change the peak shape — core 10.07 → 10.06 MeV. They buy purity, only",
        "softMuonMVA has no counterpart here: innerTrack() is absent for 100% of persisted muons",
    ],
)

# ----------------------------------------------------------------- 11
d.table(
    "The 2018 dimuon block does nothing here",
    [
        ["configuration", "sig eff", "purity", "χ² exp", "bias"],
        ["baseline", "1.000", "0.214", "749", "2.82×"],
        ["dimuon pT > 7", "1.000", "0.214", "750", "2.82×"],
        ["dimuon vtxProb > 0.1", "0.977", "0.215", "711", "2.82×"],
        ["dimuon αBS < 0.4", "0.347", "0.296", "178", "2.08×"],
        ["B vtxProb > 0.1", "0.855", "0.262", "1239", "2.95×"],
        ["muon pT > 4", "0.999", "0.214", "753", "2.83×"],
        ["all five together", "0.290", "0.376", "274", "2.15×"],
    ],
    lead="χ² exp is the exponential fit to data; bias is its yield against the true shape.",
    widths=[0.36, 0.16, 0.16, 0.16, 0.16],
    note="Three are no-ops. B vtxProb makes the background WORSE — it sculpts toward the B mass. "
    "αBS is the only one that helps and it costs 65% of the signal. The combinatorial contains a "
    "REAL J/ψ from a real b hadron, so no dimuon variable can separate it.",
)

# ----------------------------------------------------------------- 12
d.figure(
    "So the handle has to be the bachelor track",
    "fig_id_d0sig.png",
    lead="Combinatorial spikes at zero d0 significance; signal is broad. Every pT slice.",
)

# ----------------------------------------------------------------- 13
d.table(
    "The kaon ID, measured — and it beats the pT floor",
    [
        ["", "baseline", "kaon pT > 1", "|d0/σ| > 1.5", "2018 dimuon"],
        ["peak purity", "0.214", "0.579", "0.677", "0.376"],
        ["signal efficiency", "1.000", "0.968†", "0.858", "0.290"],
        ["exponential χ²/ndf", "745.8", "41.2", "25.4", "273.6"],
        ["quadratic χ²/ndf", "85.5", "8.5", "6.8", "58.1"],
        ["yield bias", "2.8×", "1.11×", "1.073×", "2.15×"],
    ],
    lead="Displacement, not particle identity: the background is prompt and the signal is not.",
    widths=[0.32, 0.17, 0.17, 0.17, 0.17],
    note="†of candidates above 1 GeV only. The kaon ID does not sculpt the signal — template width "
    "15.37 → 15.50 MeV. dE/dx was tested and dropped: rejection 0.108 against 0.370 for d0 "
    "significance, and the estimator is missing for 26.9% of signal tracks below 0.4 GeV.",
)

# ----------------------------------------------------------------- 14
d.table(
    "But the soft kaons are the best-measured signal",
    [
        ["kaon pT [GeV]", "signal n", "core width", "% of signal", "% of bkg"],
        ["0.15 – 0.50", "1 182", "7.64 MeV", "0.51", "33.08"],
        ["0.50 – 1.00", "5 183", "7.84 MeV", "2.58", "36.15"],
        ["1.00 – 1.50", "6 495", "8.03 MeV", "4.02", "15.04"],
        ["2.00 – 3.00", "11 648", "10.41 MeV", "16.16", "4.78"],
        ["5.00 – 10.0", "16 411", "12.44 MeV", "48.42", "0.46"],
    ],
    lead="σ(p)/p is roughly flat, so a soft kaon carries a smaller ABSOLUTE error into the B mass.",
    widths=[0.24, 0.19, 0.21, 0.18, 0.18],
    note="A factor 1.6 across the range, so a pT cut discards the best-measured signal — the "
    "template got WIDER, 10.93 → 11.55 MeV. Below 1 GeV sits 7.1% of the signal and 84.3% of the "
    "background. A fixed |d0/σ| threshold keeps only 38.5% below 0.4 GeV against 88.4% above 1.5.",
)

# ----------------------------------------------------------------- 15
d.bullets(
    "Open questions for the meeting",
    [
        "Keep the soft kaons? 7% of the signal, 84% of the background, and the best resolution",
        ("Answering it needs ∂m/∂A and ∂m/∂e per candidate — not yet computed", 1),
        "If we keep them: what functional form has a maximum? rabbit's backgrounds are monotonic",
        "Is a data/simulation resolution scale, floated in the fit, the right way to settle the width?",
        "Do we chase 2:1 at all? On signal we are at 15:1, so it needs ~30× more simulation",
        ("At the measured efficiencies that is ~10¹⁰ generated events, against 4×10⁸ so far", 1),
        "Is a pT-biased kaon ID acceptable, or must the threshold hold efficiency constant?",
    ],
)

# ----------------------------------------------------------------- backup
d.section("Backup")

d.table(
    "Every width number in this deck, and why none is final",
    [
        ["input", "data σ", "sim total σ", "ratio"],
        ["2018, full selection", "15.25 MeV", "14.87 MeV", "1.025"],
        ["2018, no BDT / softMVA", "15.25 MeV", "14.87 MeV", "1.025"],
        ["CVH baseline", "23.62 MeV", "22.25 MeV", "1.062"],
        ["CVH, kaon pT > 1", "17.50 MeV", "16.12 MeV", "1.085"],
        ["joint CVH, full range", "26.50 MeV", "24.62 MeV", "1.076"],
        ["dimuon only, CVH tracks", "24.37 MeV", "23.12 MeV", "1.054"],
    ],
    lead="Common window, one method. All line-subtracted, so absolute values are inflated.",
    widths=[0.40, 0.20, 0.20, 0.20],
    note="Four treatments of the same data gave 1.71, 1.076, 1.101 and 0.957, differing only in how "
    "the background under the peak was modelled. I retracted the 1.71. The gap is real and NOT yet "
    "reliably measured — settle it inside the fit, not with subtractions.",
)

d.table(
    "The gen match: dropping the pT window costs nothing",
    [
        ["", "dR<0.02 AND pT window", "dR<0.02 only"],
        ["three legs matched", "78.5%", "79.7%"],
        ["max |resp − 1|", "0.0981", "0.5283"],
        ["matched kaons beyond 10%", "0.0000", "0.0139"],
        ["response rms", "0.0250", "0.0327"],
        ["exclusive signal candidates", "383", "390"],
        ["wrong-match rate", "0.4%", "0.4%"],
    ],
    lead="Same four files through both configurations, 13 279 gated candidates each.",
    widths=[0.42, 0.29, 0.29],
    note="The pT window was truncating the kaon response at exactly 10% by construction, which is "
    "the quantity the A/e/M nuisances need. Removing it costs no purity and gains 1.2 points of "
    "matching efficiency. Production setting: dR < 0.02, no pT window.",
)

d.table(
    "How short the sample is depends which target you pick",
    [
        ["target", "its data anchor", "MC signal needed", "short by"],
        ["phase-1, 1:1", "~5×10⁵ published", "~5×10⁵", "5.7×"],
        ["the note's 2:1", "~5×10⁵ published", "10⁶", "11.4×"],
        ["2:1 at OUR selection", "1.17×10⁶ measured", "2.34×10⁶", "30×"],
    ],
    lead="We have 87 776 selected simulated signal events, 77 413 of them in the fit window.",
    widths=[0.28, 0.26, 0.26, 0.20],
    note="The ~5×10⁵ anchor is a PUBLISHED fitted yield after a tight selection; ours is far looser, "
    "so rows 1–2 compare our MC to somebody else's data selection. The note flagged this. And 2:1 "
    "was never this production's target — phase-1 asked for 1:1 and 2:1 was scoped as phase-2.",
)

d.table(
    "The generation hit its number; the efficiency model did not",
    [
        ["factor", "planned", "actual", "ratio"],
        ["events generated (phase-1)", "4.0 × 10⁸", "4.03 × 10⁸", "1.01"],
        ["generator filter efficiency", "0.0208", "0.016802", "0.81"],
        ["reco + selection efficiency", "absent from the model", "0.184", "0.18"],
    ],
    lead="Every one of the 2 520 files was processed, so nothing is left on the table.",
    widths=[0.37, 0.24, 0.21, 0.18],
    note="The plan's model was generated = target / (inclusive fraction × filter eff), which assumes "
    "EVERY gen-filtered signal event becomes a selected candidate. Measured it is 0.184, and the "
    "trigger alone is 68% of that loss. A future plan must carry a reco+selection term.",
)

d.bullets(
    "And 2:1 is not the requirement anyway",
    [
        "2:1 is a rule of thumb; nothing derives it from a target on the scale uncertainty",
        "What binds a template fit is the template's statistical noise PER CELL, and that is measured",
        ("Signal template σ/value: 0.36% over the whole window, 9.4–14.0% in the edge cells of 288", 1),
        "So size the next production against per-cell noise at the binning we intend to fit",
        "Occupancy already supports 3–4 kaon pT bins, not the ~9 projected before the full run",
    ],
)

d.figure_row(
    "The impact parameter was referenced to the origin",
    ["d0_before.png", "d0_after.png"],
    lead="Gap G13, fixed at source: shape χ² per bin 2340 → 3.92, a factor 597.",
)

d.figure_row(
    "By eye: the 2018 dimuon block against the kaon ID",
    ["stack_residual_s18all.png", "stack_residual_kaonid.png"],
    lead="Left: all five 2018 dimuon cuts. Right: |d0/σ| > 1.5. Combinatorial at the bottom.",
)

d.figure(
    "The background is 500 MeV wide and has a maximum",
    "bkg_fullrange.png",
    lead="Over the whole stored range: one broad maximum, sitting on the B mass.",
)

d.figure(
    "The prefit an exponential can almost carry",
    "prefit_kpt1.png",
    lead="Kaon pT > 1 GeV: 143 371 signal, χ²/ndf 41.02 against 745.63 before.",
)

d.figure(
    "Swap only the background, and the fit closes",
    "bkg_models.png",
    lead="Same data, same template. Green is the true simulated shape; red is an exponential.",
)

d.bullets(
    "What is implemented, and what is not",
    [
        "Histmaker, datagroups, tensor writer, fit driver, selection scans, figure scripts",
        "Two processes, no systematics — every free parameter belongs to a rabbit param model",
        "Kaon A/e/M nuisances are now unblocked: per-leg gen kinematics arrived with G2",
        "Not done: full-statistics pass with a chosen working point, then the differential fit",
        "Not done: absolute normalisation — the cross section still needs a GEN-only job",
    ],
)

print(d.save("260819_cvhnano_first_look_summary.pptx"))
