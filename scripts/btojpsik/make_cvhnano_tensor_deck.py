"""Deck: the CVH NanoAOD B+ -> J/psi K tensor and fit, first results.

Figures come from the scratchpad script that produced
public_html/slides/assets_cvhnano_tensor/. Every number here is measured on the
full data (22 327 files) and the full simulation (2 520 files); none is projected.
"""

import sys

sys.path.insert(0, "/home/submit/pmlugato/.claude/skills/slides")
from deck import Deck

ASSETS = "/home/submit/pmlugato/public_html/slides/assets_cvhnano_tensor"

d = Deck(assets=ASSETS)

d.title_slide(
    "A tensor and a fit for B⁺ → J/ψ K⁺",
    "The machinery works. The analytic background does not — and a track cut fixes most of it.",
    "Pietro Lugato  ·  17 August 2026  ·  mZ / CVH NanoAOD",
)

d.bullets(
    "What this was built to answer",
    [
        "How many events there are, in data and in the simulated signal",
        "What statistical power the simulation has, per cell and per mass bin",
        "What binning the sample can support",
        "What is missing, named with the measurement that shows it",
        ("A converged fit was never on that list — and did not arrive", 1),
    ],
)

d.bullets(
    "Where this goes",
    [
        "Inputs, the signal template, and two histmaker fixes",
        "Binning, occupancy, statistical power, seeding — what the sample supports",
        "The fits: Asimov closes, data does not, and the background model is why",
        "What the background actually is, and why 2018 did not have this problem",
        "The selection study: one cut family fixes it, the other makes it worse",
    ],
)

d.table(
    "The inputs, both samples in full for the first time",
    [
        ["", "data", "simulation"],
        ["files", "22 327", "2 520"],
        ["event loop", "404 s", "152 s"],
        ["selected events", "18.46M", "887 100"],
        ["in the [5.20, 5.35] window", "9 841 674", "480 227"],
        ["of which signal category", "—", "77 413"],
    ],
    lead="First full-data pass; Phase 2 used 4 000 files and extrapolated.",
    widths=[0.44, 0.28, 0.28],
    note="data / simulation = 20.8, against the 18.9 projected. The extrapolation was 12% low. "
    "Simulation reproduces the 887 100 baseline exactly.",
)

d.bullets(
    "The signal template, precisely",
    [
        "`jointCvh_mass`: the joint 3-body CVH refit, which is the `--arm` default",
        "genCategory bin 0 of the inclusive sample: candidate matched to a B⁺ at dR < 0.01",
        "Excited states folded first — B*⁺ wins the dR match in a quarter of true signal",
        "Not exclusive B⁺ → J/ψ K⁺: Gen has no mother index (G1), so the chain is unwalkable",
        ("Shape held fixed; only the normalisation floats, one parameter per cell", 1),
    ],
)

d.bullets_figure(
    "Two histmaker fixes, and one recovered 12.8% of the signal",
    [
        "It was filling the 5-bin candidate-charge diagnostic axis, not a kaon charge axis",
        "It carried no variances, so the template's statistical cost was unmeasurable",
        "Folding the pT overflow adds 11 275 signal events: earlier counts were 13% low",
    ],
    "fig_overflow.png",
)

d.table(
    "The mass binning, and a guard that caught its own default",
    [
        ["window", "MeV/bin", "bins per FWHM", "verdict"],
        ["[5.00, 5.50]", "50.0", "0.65", "rejected"],
        ["[5.15, 5.40]", "25.0", "1.30", "warned"],
        ["[5.20, 5.40]", "20.0", "1.62", "warned — the 2018 binning"],
        ["[5.20, 5.35]", "15.0", "2.17", "default"],
    ],
    lead="Signal core 13.8 MeV, so the peak is 32.5 MeV wide at half maximum.",
    widths=[0.26, 0.18, 0.24, 0.32],
    note="The proposal said 'bin width below the core width' and then chose 20 MeV bins, so the "
    "guard rejected its own default on first use. The threshold is the peak width, not σ.",
)

d.figure(
    "288 cells, and 90% of them clear the 100-event floor",
    "fig_occupancy.png",
    lead="3 pT quantile bins × 48 η × 2 charge, on the full simulation.",
)

d.table(
    "The simulated statistical power, where it is spent",
    [
        ["simulated component", "events", "effective", "σ / value"],
        ["signal template, whole window", "77 413", "77 413", "0.36%"],
        ["non-signal, whole window", "402 814", "278 650", "0.19%"],
        ["signal, softest-pT edge cell", "112", "112", "9.4%"],
        ["signal, hardest-pT edge cell", "51", "51", "14.0%"],
    ],
    lead="Effective entries are (Σw)²/Σw², so they measure power rather than count.",
    widths=[0.40, 0.20, 0.20, 0.20],
    note="The signal category is unweighted by construction — the forced-decay weight is "
    "normalised to the B⁺ — so its σ is √N. The non-signal categories are not, and lose 31% "
    "of their power to the weighting. Every band in this deck reads Σw², not N.",
)

d.table(
    "Seeding: two corrections, both forced by running it",
    [
        ["signal scale from", "value", "cells with no room for background"],
        ["peak-region data/MC ratio", "29.16 *", "24 of 36"],
        ["… sideband-subtracted", "2.59 *", "22 of 36"],
        ["… capped at the physical bound", "1.11 *", "0 of 36"],
        ["full samples: estimate, then cap", "29.16 → 13.94", "1 of 288"],
    ],
    lead="A signal template cannot exceed the data, so the cap is a bound, not a knob.",
    widths=[0.42, 0.20, 0.38],
    note="* diagnosed on a 40-file pair, where the failure was first visible. The full-sample "
    "sideband estimate implies a peak purity of 0.368 and is still cut in half by the bound.",
)

d.figure(
    "Why one global scale cannot start every cell near unity",
    "fig_ptmismatch.png",
    lead="Data is 89% soft-kaon background; the template is a third in each bin.",
)

d.figure(
    "Prefit at two edge cells of the 288",
    "fig_edge2.png",
    lead="Bands are the template's own √(Σw²): 51 events on the right, ±58% in one bin.",
)

d.table(
    "The fits: Asimov closes, data does not",
    [
        ["configuration", "dataset", "edmval", "χ²/ndf"],
        ["288 cells", "Asimov", "2.0 × 10⁻²⁷", "0"],
        ["288 cells", "data", "status 2", "1.6 × 10⁵"],
        ["2 cells (charge)", "Asimov", "converged", "0"],
        ["2 cells (charge)", "data", "1462", "2.0 × 10⁴"],
    ],
    lead="Same model, same templates, same 288 cells — only the observation changes.",
    widths=[0.30, 0.20, 0.24, 0.26],
    note="No Cholesky failure and no boundary parameter, so the degeneracy documented at 448 "
    "cells is not what limits this channel. Failing on data only rules out conditioning.",
)

d.figure(
    "It is not a satellite peak. It is 500 MeV wide.",
    "fig_fullrange.png",
    lead="Over the whole stored range: one broad maximum, sitting on the B mass.",
)

d.figure(
    "Inside the window, that top looks like a bump",
    "fig_bkgshape.png",
    lead="Rises 70%, peaks 3 MeV from the signal, falls. Both models are monotonic.",
)

d.figure(
    "Two orders of magnitude, and a factor 2.2 on the yield",
    "fig_models2.png",
    lead="Best achievable fit of each family; the simulated band is 0.4%, invisible here.",
)

d.table(
    "What each background model does to the answer",
    [
        ["background", "χ²/ndf", "signal yield", "fraction of data"],
        ["exponential", "6 250", "2.55M", "0.259"],
        ["straight line", "6 228", "2.54M", "0.258"],
        ["simulated non-signal", "61", "1.17M", "0.119"],
    ],
    lead="The two analytic families are indistinguishable from each other and both wrong.",
    widths=[0.34, 0.18, 0.24, 0.24],
    note="The simulated background's own scale comes out 21.5 against an overall data/simulation "
    "ratio of 20.5 — a 5% consistency check that its shape describes the data.",
)

d.bullets(
    "Where that leaves the yield, and the production",
    [
        "Signal in the window is ~1.17M events, from the two-template estimate, not from a fit",
        "The note assumed ~5 × 10⁵, so 2:1 costs roughly 2.3× more simulation than planned",
        "An analytic background would have reported 2.55M and been wrong by a factor 2.2",
        "More simulation would not make the analytic fit work — the shape is wrong, not the errors",
    ],
)

d.bullets(
    "Why it worked in 2018: purity, not shape",
    [
        "Same measure, same window, on the 2018 histograms: the sample is 80% signal",
        "The CVH NanoAOD is 12% signal in the same window — a factor of seven",
        "At 20% background a wrong background shape barely moves the yield; at 88% it doubles it",
        "So the 2018 background may well have curved too. It did not matter there.",
        ("The target is therefore purity, and no window scan can deliver it", 1),
    ],
)

d.table(
    "The 2018 selection is not the difference",
    [
        ["2018 cut", "effect here", "verdict"],
        ["dimuon pT > 7", "signal 1.000, gain 1.000", "trigger already does it"],
        ["dimuon sl3d > 4", "signal 0.948, gain 1.003", "AlCaReco already does it"],
        ["dimuon vtx prob > 0.1", "signal 0.978, gain 1.001", "no-op, median 0.69"],
        ["muon, kaon |η| < 1.4", "signal 0.79, 0.72", "bites — pure acceptance"],
        ["all of it together", "signal 0.505, bias 2.40→2.47×", "rejected"],
    ],
    lead="Ten of twelve cuts have a counterpart here. Most were already imposed upstream.",
    widths=[0.30, 0.36, 0.34],
    note="Two thresholds do not transfer: our dimuon αBS has signal median 0.614 rad, so "
    "2018's < 0.4 cuts two thirds of the signal — the dimuon does not point back to the "
    "beamspot, because the missing kaon carries real momentum.",
)

d.figure(
    "One cut family fixes it; the other makes it worse",
    "fig_selscan.png",
    lead="Bias = yield with an exponential ÷ yield with the true shape. 300 files, matched.",
)

d.bullets(
    "Why the two families behave oppositely",
    [
        "Fit-quality cuts — χ²/ndof, nIter, vertex probability — ask how well the B fit worked",
        "They raise S/B, and pile the background up at the B mass: bias 2.40 → 2.6–2.7×",
        "Track-quality cuts — |d0|, hit count — judge the kaon alone, with no B hypothesis",
        "They raise S/B and flatten the background: |d0| > 0.01 with hits > 10 gives 1.17×",
        ("So purity and flatness pull against each other, but only for the first family", 1),
    ],
)

d.table(
    "The working points, and what each costs",
    [
        ["selection", "signal eff", "S/B", "max/min", "yield bias"],
        ["no cuts", "1.000", "0.195", "1.749", "2.40×"],
        ["|d0| > 0.01", "0.811", "0.248", "1.462", "1.61×"],
        ["|d0| > 0.01, hits > 10", "0.720", "0.400", "1.317", "1.17×"],
        ["|d0| > 0.02, hits > 13", "0.516", "0.657", "1.221", "1.06×"],
        ["|d0| > 0.05", "0.358", "0.687", "1.226", "1.08×"],
    ],
    lead="Recommendation: |d0| > 0.01 with hits > 10 — 72% of the signal, bias 2.40 → 1.17×.",
    widths=[0.32, 0.17, 0.15, 0.17, 0.19],
    note="Signal efficiency is the cost, not a constraint, but it is not free either: 72% "
    "raises the simulated sample the 2:1 target needs by 1.4×. Chosen over the 1.06× point "
    "because that one costs half the signal to remove a bias already near one.",
)

d.table(
    "The quick refit: 5× better, and now it does not converge",
    [
        ["inclusive fit, data", "no cuts", "|d0| > 0.01, hits > 10"],
        ["saturated 2ΔNLL", "360.1", "72.7"],
        ["edmval", "4.3 × 10⁻¹⁵", "3.26"],
        ["minimiser", "converged", "did not"],
        ["seeded signal scale", "2.85", "1.68"],
    ],
    lead="Two charge cells, no pT or eta binning — so statistics explains nothing here.",
    widths=[0.34, 0.30, 0.36],
    note="Goodness of fit improves 5×, which is the point. The minimiser stops closing because "
    "a smaller background fraction leaves the exponential's slope weakly constrained — a "
    "conditioning problem, not the shape problem. rabbit's 'linear chi2' returns 1.7e16 here.",
)

d.bullets(
    "Where this leaves us",
    [
        "The fit model does not need changing: signal template plus analytic background",
        "It needs a purer sample, and one track cut gets the bias from 2.40× to 1.17×",
        "softMuonMVA and the BMM BDT are not implicated — the available cuts explain it",
        "Next: full-statistics pass with |d0| > 0.01, hits > 10, then the differential fit",
        ("Open: the pT bias of a d0 cut, which is exactly what a soft kaon suffers", 1),
    ],
)

d.section("Backup")

d.figure(
    "The commands: histograms",
    "fig_cmd_hist.png",
    lead="One pass, one file. --filterProcs is variadic, so there is no merge step.",
)

d.figure(
    "The commands: the selection scan",
    "fig_cmd_scan.png",
    lead="Matched 300-file passes, so a shape difference cannot come from file choice.",
)

d.figure(
    "The commands: tensor and fit",
    "fig_cmd_fit.png",
    lead="The driver writes every expanded command to a .sh beside its log.",
)

d.bullets(
    "Two rabbit interface traps, copied from the old chain",
    [
        "rabbit takes -o/--outpath; --output no longer exists and the chain died parsing it",
        "-t is inverted from combinetf: -t -1 is Asimov and the default, -t 0 fits data",
        ("The old run_full_chain.sh passes -t 0, so every fit would have been mislabelled", 1),
        "With -t -1 the Asimov fit closes at χ²/ndf 0 — the check that it is now right",
    ],
)

d.table(
    "Per-cell occupancy at 288 cells, full simulation",
    [
        ["quantity", "sum", "min", "median", "≥ 100"],
        ["all candidates", "480 227", "31", "404", "268/288"],
        ["signal category", "77 413", "31", "307", "260/288"],
        ["signal in mass peak", "63 668", "17", "242", "235/288"],
    ],
    lead="Phase 2 predicted 92.0% and 80.2% at 3 pT bins; measured 90.3% and 81.6%.",
    widths=[0.32, 0.18, 0.14, 0.18, 0.18],
    note="Occupancy is not the binding constraint: every cell fits on Asimov, none on data.",
)

d.bullets(
    "Three bugs the selection study turned up",
    [
        "`--minVtxProb` was declared and never wired in — silently off in every run so far",
        "`|d0| > x` on a −98 sentinel *passes*: the one cut taken on an absolute value",
        "matplotlib parsed `$HIST -o $MZ` as mathtext and ate the flag on a command slide",
        ("scipy's least_squares raises a bogus \"memory error\" on a 10×2 SVD here", 1),
        ("Both fits are conditionally linear, so they are now solved, not iterated", 1),
    ],
)

d.bullets(
    "What is implemented, and what is not",
    [
        "Tensor writer, helpers, fit driver, selection comparison and figure scripts",
        "Two processes, no systematics — every free parameter belongs to a rabbit param model",
        "Kaon A/e/M nuisances still blocked by G2: no per-leg gen kinematics in this NanoAOD",
        "Not done: full pass with the chosen cuts, the pT-bias check, the differential fit",
    ],
)

print(d.save("cvhnano_tensor_first_fits.pptx"))
