#!/usr/bin/env python3
"""Overview deck: where B+ -> J/psi K stands, in the context of a combined
J/psi + Z + D0 + B+ momentum-scale fit.

Sources: the six CVH NanoAOD decks (12/13/17/18/19 August), the A/e/M deck
(25 August, now 115 slides), the joint mZ + J/psi proposal (27 August), and
openspec/changes/add-cvhnano-kaon-aem-variations.
"""
import sys

sys.path.insert(0, "/home/submit/pmlugato/.claude/skills/slides")

from deck import Deck  # noqa: E402

ASSETS = "/home/submit/pmlugato/public_html/slides/assets_btojpsik_overview"
d = Deck(assets=ASSETS)

d.title_slide(
    "B⁺ → J/ψ K: where we stand",
    "One channel of a combined J/ψ + Z + D⁰ + B⁺ momentum-scale fit — "
    "what works, what does not, and what is genuinely open",
    "Pietro Lugato  ·  28 August 2026  ·  mZ / muon calibration",
)

d.bullets(
    "Outline",
    [
        "Why B⁺ is in the combined fit at all — the lever arm no other channel has",
        "The chain, end to end: AlCaReco → CVH refit → NanoAOD → histmaker → tensor → fit",
        "The fit problem, which is the whole story: no data fit has ever converged",
        "What is genuinely open — resolution, the pT floor, the A–e degeneracy",
        "What limits us, what to produce, and where this joins the other channels",
    ],
)

# ---------------------------------------------------------------- why B+
d.section("Why this channel")

d.figure(
    "Only the soft bachelor separates A from e",
    "lever_arm.png",
    lead="δk/k = A − e·k + qM/k. Z lives at k ≈ 0.02, the B⁺ kaon reaches k ≈ 1.",
)

d.bullets(
    "0.85% of the Hessian trace, 76% of the information",
    [
        "A trace share is not a measure of value — the channel is small and not redundant",
        "It reaches curvatures no dimuon channel does, so it breaks the A–e degeneracy "
        "that a Z-only or J/ψ-only fit cannot",
        "It carries a mass constraint the K⁰s channel has none of, and its kaon is a "
        "hadron, so it also probes the energy-loss term directly",
        "The cost of being wrong here is a bias that propagates into mZ with the "
        "calibration — which is exactly the coupling the combined fit exists to expose",
    ],
)

d.bullets_figure(
    "What the channel constrains: e, weakly A, barely M",
    [
        "Shape sensitivity per nuisance, one representative cell of 304 signal candidates",
        "e reaches S ≈ 0.5, A ≈ 0.2, M ≈ 0.04 — only e is close to a measurement",
        "The muon legs carry 13–68% of S — the J/ψ constraint does not suppress them",
    ],
    "sensitivity_and_muons.png",
)

# ---------------------------------------------------------------- chain
d.section("The chain")

d.bullets(
    "Everything from AlCaReco to a fitted tensor runs",
    [
        "Our own TkAlJpsiX AlCaReco → joint 3-body CVH refit → NanoAOD → histmaker → "
        "tensor → rabbit, on the full 2016postVFP sample",
        "87.3M data nano events (16.89 fb⁻¹) and 6.77M simulated; 6.25M and 283k "
        "selected candidates in the fit window",
        "Peak at 5.2798 GeV in data against 5.2800 in simulation, both within 0.5 MeV "
        "of PDG, and all three fit arms agree",
        "Exclusive B⁺ → J/ψ K⁺ truth is definable per leg, so every simulated "
        "component is identified rather than lumped",
        "The A/e/M kaon nuisances exist, as an ONNX response reweight, and are in the fit",
    ],
)

d.figure_row(
    "The peak is clean and the three fit arms agree",
    ["peak_arms.png", "stack_jointCvh_mass.png"],
)

d.bullets(
    "The numbers that frame everything else",
    [
        "**1,363,434** signal events in data (fitted), against **72,800** in simulation "
        "— a factor **18.7**",
        "Equivalently: the signal simulation is worth **0.90 fb⁻¹** against **16.89 fb⁻¹** "
        "of data",
        "Purity 0.218 in data — against 2018's ~80%, which is the single number that "
        "explains why 2018's fit worked and ours does not",
        "Cross-check on the yield: σ(bb̄) × L × BRs gives 4.19×10⁸ produced, so "
        "ε_total = 0.33% — right for a displaced dimuon trigger on a soft b spectrum",
    ],
)

# ---------------------------------------------------------------- the fit
d.section("The fit problem")

d.bullets(
    "Read this first: no data fit has ever converged",
    [
        "The reference fit stopped at iteration 7199 with edmval **72,668**, against a "
        "1×10⁻⁵ threshold — and wrote a perfectly readable 442 kB results file",
        "So every data-side number from before that check — the postfit χ²/ndf 34, "
        "A–e = 0.63 on data, the 6–8σ pulls, the impacts — is not yet a measurement",
        "Asimov converges. Toys converge, at 10⁻¹³ to 10⁻¹⁷. Only the observation fails",
        "Ruled out one at a time: the minimiser, the iteration budget, bin-by-bin "
        "nuisances, polynomial order 2–5, all three MC-template variants, the "
        "resolution nuisance, the scale nuisances entirely, windows ±50 to ±200 MeV",
    ],
)

d.bullets(
    "It is model error measured against data precision",
    [
        "Weakening the data by a factor 333 (`dataScale` 0.003, same shape) turns the "
        "Hessian positive-definite and edmval finite at 44.4 — 1637× better",
        "So the failure is not conditioning and not statistics: the model is wrong by a "
        "fixed fraction, and 6.2M events are enough to notice",
        "**More simulation does not fix this. Only a better model shape does.**",
        "I had also said this explains 2018. Measured since, it does not — next slide",
    ],
)

d.figure(
    "2018's residual was scale. Ours is shape.",
    "shift_and_binwidth.png",
    lead="One free mass shift takes 2018 from χ²/ndf 66 to 1.44. It takes us 54 to 44.",
)

d.bullets(
    "Two things were confounded: the shift, and the window",
    [
        "2018's entire residual was a **−1.44 MeV offset** between data and template — "
        "a momentum-scale difference, which is the measurement, not a defect",
        "Δχ² = 6145 for that one parameter, on 1.9M events. **2018's lineshape model is "
        "correct**; ours has a −0.81 MeV offset and removing it leaves 43.9",
        "Coarse binning is not the explanation: the 2018 tensors did use 20 MeV bins, "
        "but rebinning data and model together makes χ²/ndf **worse**, not invisible",
        "But those two χ² were measured on **different windows** (±100 vs ±200 MeV), so "
        "they were never comparable — corrected on the next slide",
    ],
)

d.table(
    "Matched window and model: purity does most of it",
    [
        ["sample", "N data", "purity", "χ²/ndf", "shift"],
        ["2018, full selection", "2.03M", "0.717", "**1.44**", "−1.44"],
        ["2016 CVH, all kaon pT", "3.78M", "0.288", "5.33", "−1.03"],
        ["2016 CVH, pT > 1.8", "1.59M", "0.547", "3.57", "−1.21"],
        ["2016 CVH, **pT > 3.2**", "0.74M", "**0.755**", "**2.89**", "−1.41"],
        ["2016 CVH, pT > 5.0", "0.37M", "0.819", "3.16", "−1.46"],
    ],
    lead="Same window and model: signal ⊗ Gauss, free shift, exp·poly2.",
    widths=[0.32, 0.15, 0.15, 0.18, 0.20],
    note="Purity is N_sig/N_total in the DATA inside the window, so it depends on the "
    "window as much as the selection. At matched purity a factor 2 remains and stops "
    "improving — that is the genuine shape work. The −1.4 MeV offset is common to both.",
)

d.bullets(
    "So we can have 2018's purity or B⁺'s lever arm, not both",
    [
        "Purity differs for three reasons: the window (0.288 vs 0.218 on identical "
        "data), two discriminants we cannot compute, and a harder kaon pT cut",
        "The BMM BDT has no counterpart here, and `softMuonMva` needs `innerTrack()` — "
        "absent for 100% of persisted muons in our AlCaReco. Together, 64% of survivors",
        "We CAN reach purity 0.755 with kaon pT > 3.2 — but that keeps 19.5% of the data "
        "and removes everything above k ≈ 0.31, which is why B⁺ is in the fit at all",
        "A bachelor displacement cut is the only handle that might buy purity without "
        "spending the lever arm",
        "Also: at matched purity the fitted smear is 3.34 MeV against 2018's 2.88 — so "
        "the resolution mismatch is smaller than the 6.26 MeV I measured on the wide window",
    ],
)

d.figure(
    "The model is wrong at 1%; the data has 0.36% per bin",
    "inclusive_datafit.png",
    lead="One cell, 6.2M events, 80 bins, seconds per fit — χ²/ndf = 53.6.",
)

d.bullets(
    "The residual names its own two causes",
    [
        "An **S-shape across the peak** — +8.9σ at 5.2775, then −11.8 and −13.7 just "
        "above. The data falls off faster on the high side than a symmetric smear can",
        "A **wave in the low shoulder** — +4% at 5.08–5.13, then −5% at 5.15–5.20, "
        "sitting exactly where the simulated B⁰ falls from 17% to 1% of the MC",
        "That second one is a partially-reconstructed kinematic edge, and no polynomial "
        "follows an edge",
        "χ² shares: shoulder 52.3%, peak 32.7%, upper 11.0%",
    ],
)

d.figure(
    "Two fixes, factor 5.6: asymmetric resolution and a B⁰ edge",
    "inclusive_datafit_physics.png",
    lead="χ²/ndf 53.6 → 10.1. B⁰ and J/ψKX die at 5.15; J/ψπ⁺ collapses to zero.",
)

d.figure(
    "And they beat pure flexibility at equal parameter count",
    "model_ladder.png",
    lead="The grey line is the control: nine parameters of polynomial reach only 22.4.",
)

d.bullets(
    "What that does and does not settle",
    [
        "Resolution is worth Δχ² 2839 for two parameters; the peaking templates a "
        "further 2585 — comparable, and each dwarfs any background order",
        "Both were things I had argued against. \"Don't add peaking channels\" is "
        "withdrawn: they do not peak at the B mass, they put an edge in the shoulder",
        "**χ²/ndf 10.1 is not 1.** The remaining pulls redistribute rather than vanish, "
        "so the residual is now spread rather than localised",
        "Neither fix is in the 288-cell fit yet — that is the next thing to run, and the "
        "first change in the campaign with a measured factor 5.6 behind it",
    ],
)

d.figure(
    "Half of what is left is fitting 288 cells with one shape",
    "band_chi2_and_smear.png",
    lead="Split by kaon pT, χ²/ndf falls 9.3 → ~4, and the smear runs 9.2 → 4.1 MeV.",
)

# ---------------------------------------------------------------- open
d.section("Genuinely open")

d.figure(
    "The response network's floor flips the sign of e",
    "pt_floor_variations.png",
    lead="22.7% of the signal sits below the 2 GeV gen-pT training floor.",
)

d.bullets(
    "Why that one matters more than it looks",
    [
        "Above the floor the two treatments agree to 0.0000 MeV; at 1.0–1.5 GeV they "
        "differ by 0.463. Clamping asserts a flat residual, extrapolating leaves the "
        "training domain — neither is validated",
        "The spread lands on **e** — the one parameter whose entire physics is low-pT "
        "energy loss — and it is 132% of the variation in the softest band",
        "Immediate: take the difference as a shape nuisance. Both histmaker outputs "
        "already exist, so it costs a tensor option and no new running",
        "Do NOT cut at 2 GeV instead: that loses 22.7% of the signal AND halves e's "
        "lever arm, since the shift range collapses from [0.071, 0.352] to [0.071, 0.191]",
        "The real fix is retraining the response model down to 1 GeV",
    ],
)

d.figure(
    "Matched to 2018's acceptance, the CVH fit is 5% better",
    "cmp1618_mass_2018acc.png",
    lead="15.50 MeV against 16.30, in the acceptance 2018 actually had.",
)

d.bullets(
    "Resolution, after four wrong answers",
    [
        "Four treatments of the same data gave width ratios 1.71, 1.076, 1.101 and "
        "0.957, differing only in how the background under the peak was modelled",
        "The 8.9 MeV wired into `RESVAR` came from √(17.29²−14.81²) — a "
        "data-vs-MC-**signal** comparison, the one the 2018 script warns against",
        "Fitting the smear with the background free gives **6.26 MeV** as a Gaussian, "
        "or **4.77 MeV with a low-side width ratio 1.27**",
        "Settle it inside the fit with a floated, asymmetric resolution — never by "
        "subtracting a background and measuring what is left",
    ],
)

d.figure(
    "The A–e degeneracy is real and it is where it should be",
    "scan2d_A_e_full.png",
    lead="Asimov: ρ(A,e) ≈ +0.6, cross-η elements zero — four separate measurements.",
)

d.bullets(
    "The soft-kaon question, still not decided",
    [
        "Below 1 GeV sits 7.1% of the signal and 84.3% of the background — and the "
        "**best-measured** signal, because σ(p)/p is flat so a soft kaon carries a "
        "smaller absolute error into the B mass",
        "The pT > 1 cut buys 2.7× purity for 3.2% of the signal, and 2018 had it at "
        "production. But it discards exactly the curvature that separates A from e",
        "A displacement cut on the bachelor (|d0/σ| > 1.5) reaches similar purity at "
        "full efficiency — the background is prompt, the signal is not",
        "Not yet answered: is the A/e leverage of that 7% worth the background it brings? "
        "That needs ∂m/∂A and ∂m/∂e per candidate, which we can now compute",
    ],
)

# ---------------------------------------------------------------- limits
d.section("What limits us")

d.figure(
    "MC statistics is 80.8% of the A/e/M uncertainty",
    "mc_scaling.png",
    lead="Measured on a converged Asimov pair differing only by --noBinByBinStat.",
)

d.bullets(
    "The production ask, split by what each sample is for",
    [
        "**Signal**: exclusive B⁺ → J/ψ K⁺, **4.3M** full-sim events for 6× — which is "
        "0.63× the simulation already run, and buys a 32% uncertainty reduction",
        "An exclusive sample yields a selected signal candidate for 10.19% of its nano "
        "events against 1.08% for the inclusive one — a measured 9.48× gain",
        "**Peaking shapes**: B⁰ and J/ψπX need per-kaon-pT shapes (χ²/ndf 3.4 and 9.5 "
        "between pT halves); J/ψKX is stable and needs one shared shape",
        "Those cost 4× more inclusive (27M events, 1.6×10⁹ Pythia) — or ~0.4M events as "
        "two small exclusive samples, **68× cheaper**",
        "Combinatorial stays analytic: its MC template is 3.3× noisier per bin than the "
        "data and the simulation has no prompt or non-b fake component at all",
    ],
)

# ---------------------------------------------------------------- combined
d.section("Toward the combined fit")

d.figure(
    "The joint mZ + J/ψ machinery works and reproduces havyn",
    "joint_s19.png",
    lead="The mZ / calibration correlation, a smooth arch peaking at 0.24.",
)

d.bullets(
    "How the four channels join, and what is already in place",
    [
        "One tensor, one basis of 144 nuisances — 48 η bins × (A, e, M) — correlated "
        "**by name and nothing else** in `merge_inputs`; ρ(A11,e11) reproduces havyn's "
        "0.90 to 0.91",
        "Z + J/ψ: done and validated on Asimov, edmval 9.5×10⁻²⁸, contract check passing; "
        "the two-step chain's ×2.1 inflation and closure systematics become unnecessary",
        "D⁰ → Kπ: havyn's arm, same basis, and it answers two of our questions — it hits "
        "the same pT floor and shows what a third low-mass channel buys",
        "B⁺ → J/ψ K: the tensor and nuisances exist and the Asimov fit converges. "
        "**The data fit is the only thing standing between us and merging it in**",
        "Two silent traps already found: `--scale_e` defaults to 5.0 on the J/ψ side and "
        "1.0 on the Z side, and `--calibrationStatScaling` is 2.1 on the Z side only",
    ],
)

d.bullets(
    "The merge was correlating nothing, and now it does",
    [
        "B⁺ wrote `scaleA0…scaleA3`; the joint tensor writes `A0…A23`; "
        "`merge_inputs` correlates **by name and nothing else**",
        "So merging B⁺ as it stood would have correlated **zero** parameters — and "
        "converged, and looked entirely normal",
        "Fixed with `--scaleParamNaming joint`, guarded to refuse a wrong basis size "
        "or a partial basis, and tested by being made to fail both ways",
        "Now **72 of 72** correlate. The decisive check: adding a 72-systematic "
        "channel changed the merged systematic count by **exactly zero**",
        "All four channels converged — edmval 1.7×10⁻²⁷, 9.5×10⁻²⁸, 2.8×10⁻²⁹",
    ],
)

d.bullets(
    "And the answer is that B⁺ does not measurably improve mZ",
    [
        "σ(mZ) 0.022082 → 0.022011, a **0.3%** improvement. The muon-calibration "
        "impact goes 1.72 → 1.71 MeV",
        "σ(e) improves **1.1%** on average and 3.9% in the most forward η bin; σ(A) "
        "0.4%; σ(M) not at all",
        "ρ(A,e) moves +0.854 → +0.852. The degeneracy is **not** broken",
        "Every non-calibration group is unchanged to 1% — the sanity check that B⁺ "
        "acts where it should and nowhere else",
    ],
)

d.figure(
    "Why: B⁺ alone is a factor 5.5 weaker on e",
    "bplus_solo_vs_joint.png",
    lead="Same basis. σ(A) 0.933, σ(e) 0.587, σ(M) 0.990 against 0.120 / 0.106 / 0.069.",
)

d.bullets(
    "The lever arm was the right idea about the wrong quantity",
    [
        "Reaching k = 1/pT ≈ 1 does break the A–e degeneracy in principle — but "
        "breaking it still needs precision ALONG that direction, and B⁺ has 1/5.5",
        "In information terms B⁺ carries **3% of the total on e**, which is exactly "
        "the ~1% it moves σ(e) by",
        "The combined result IS the quadrature sum: predicted 0.1187/0.1031/0.0686 "
        "against measured 0.1195/0.1048/0.0690. Nothing unexpected is happening",
        "**Ceiling, measured:** with `--noBinByBinStat` B⁺'s own σ(e) floor is 0.334, "
        "so even with unlimited MC the combined gain is **4.7% on e, ~1% on mZ**",
        "B⁺ is limited by its own 1.36M signal events, not by the template",
    ],
)

d.bullets(
    "So B⁺ is a cross-check, not a precision contribution",
    [
        "An independent handle on a **hadron** track, in a momentum regime no dimuon "
        "channel reaches, through a completely different reconstruction",
        "A 0.3% improvement in mZ is not the point. Whether B⁺ **agrees** is",
        "And that test needs the data fit — on Asimov every channel agrees by "
        "construction",
        "Which puts the lineshape work back at the top: not a prerequisite for a "
        "better mZ uncertainty, a prerequisite for B⁺ meaning anything at all",
        "More MC is still worth requesting — for the standalone measurement and the "
        "peaking shapes. It is not justified by what it adds to mZ",
    ],
)

d.bullets(
    "Next, in order",
    [
        "1. Put the asymmetric resolution and the cell-shared peaking templates into "
        "the 288-cell fit, with resolution varying by kaon pT, and re-test on data",
        "2. Add the floor clamp-vs-extrapolate difference as a shape nuisance on e "
        "— both histmaker outputs already exist, so no new running",
        "3. Re-run this combined Asimov after 1 and 2; today's numbers use the "
        "current template and will move",
        "4. Request exclusive B⁺ → J/ψ K⁺ (4.3M, priority) and the two small peaking "
        "samples — for the standalone measurement, not for mZ",
        "5. Raise the dropped `pixel_multiplicity_stat_*` with the joint change: the "
        "merge drops 192 of them, identically with and without B⁺",
    ],
)

d.bullets(
    "With an asymmetric resolution, the shift buys nothing",
    [
        "asym res + 3 templates: χ²/ndf 10.10. Adding a free shift as well: **10.09**",
        "Δχ² = 10 for that parameter, and the fitted shift runs to +2.27 MeV without "
        "changing the fit — so it is unconstrained",
        "**The asymmetric resolution and the mass scale describe the same thing here.** "
        "Once the low-side tail is modelled there is no offset left to absorb",
        "Consequence: the resolution nuisance and A/e/M will compete for the same "
        "information, and a resolution nuisance that moves the scale is worse than none",
    ],
)

d.bullets(
    "Why D⁰ moves mZ and B⁺ does not",
    [
        "It is not the lever arm — D⁰→Kπ daughters and our bachelor both run pT 1–10, "
        "so both reach k ≈ 1. The difference is statistical power at that lever arm",
        "BR(D⁰→Kπ) = 3.95×10⁻² against BR(B⁺→J/ψK)×BR(J/ψ→μμ) = 6.1×10⁻⁵ — a factor "
        "**649**, before the ~20× larger charm cross section",
        "D⁰ constrains **two hadron legs**; B⁺ has one hadron and two muons, and the "
        "muons carry 13–68% of its sensitivity — information J/ψ already has",
        "To match D⁰'s −48% on σ(e), B⁺ would need 83× more information. Its data floor "
        "is 27× short — 37M signal events against 1.36M. Not reachable",
        "havyn: total mZ 2.61 → 1.93, e 2.01 → 1.04, ρ(A,e) 0.90 → 0.81. His baseline "
        "is ~28% weaker than ours, so the exact factor is uncertain; the direction is not",
    ],
)

d.figure(
    "Going below 1 GeV buys 17% on σ(e)",
    "lever_arm_clamp_test.png",
    lead="Three Asimov fits, identical but for the bachelor floor and the clamp.",
)

d.bullets(
    "The lever arm is real; there are few events there",
    [
        "Dropping the floor 1.0 → 0.1 GeV buys **σ(e) ×0.832** and ρ(A,e) 0.586 → 0.540, "
        "for **10.2%** more signal — a high return per event",
        "But not the 25× the linear-in-k extrapolation suggested: 89% of candidates "
        "below 1 GeV are background, and only 0.2% of signal is below 0.3 GeV",
        "With the clamp OFF the fit reports σ(e) ×0.215 — **a factor 4.6 better, and it "
        "is an artefact**",
        "The proof needs no fit: **A is a pT-independent scale, so its induced shift must "
        "be flat**. Clamped it is, −0.038 to −0.057 over two decades. Unclamped it swings "
        "to −0.380, and the e variation changes sign twice",
        "A large template variation always reads as small σ in prior units — so σ(e) = "
        "0.074 measures the extrapolation, not the detector",
    ],
)

d.bullets(
    "Which gives a free acceptance test for a retrained model",
    [
        "**A-flatness in pT gates any retrained response model** — no data, no fit, no "
        "Asimov required. It passes above the floor in both treatments and fails below "
        "it whenever the network extrapolates",
        "And it is the cheapest argument for the retrain: the defensible prize is 17% on "
        "σ(e), the unclamped fit hints at much more, and only a model valid below 2 GeV "
        "can tell the two apart",
        "Zero non-finite and zero negative entries in both — this is not a numerical "
        "failure, it is the network confidently answering a question it was never asked",
    ],
)

d.section("Backup")

d.table(
    "The event accounting, end to end",
    [
        ["", "data", "inclusive MC"],
        ["Pythia events generated", "—", "4.03×10⁸"],
        ["nano events", "87,341,047", "6,770,705"],
        ["…containing a gen B⁺→J/ψK⁺", "—", "714,190"],
        ["selected candidates", "6,246,614", "282,753"],
        ["**signal**", "**1,363,434**", "**72,800**"],
        ["purity", "0.218", "0.257"],
    ],
    lead="Full 2016postVFP, 16.894 fb⁻¹, fit window 5.05–5.45.",
    widths=[0.40, 0.30, 0.30],
)

d.table(
    "The model ladder, in numbers",
    [
        ["model", "npar", "χ²/ndf", "N_sig", "σ [MeV]"],
        ["exp·poly2, no resolution", "4", "79.9", "1.314M", "—"],
        ["+ Gaussian resolution", "5", "53.6", "1.363M", "6.26"],
        ["+ asymmetric instead", "6", "44.2", "1.364M", "5.66"],
        ["Gaussian + B⁰/J/ψπ/J/ψKX", "8", "26.4", "1.330M", "5.57"],
        ["**asymmetric + those**", "**9**", "**10.1**", "1.328M", "4.77"],
        ["exp·poly6 + Gaussian res", "9", "22.4", "1.280M", "4.55"],
    ],
    widths=[0.34, 0.12, 0.16, 0.19, 0.19],
    note="N_sig is stable to ±4% across every row, so the yield is not a "
    "background-model artefact.",
)

d.table(
    "The data-statistics ladder that located the failure",
    [
        ["dataScale", "nit", "edmval", "Hessian"],
        ["1.0 (full)", "263", "none", "indefinite"],
        ["0.3", "128", "none", "indefinite"],
        ["0.03", "201", "none", "indefinite"],
        ["**0.003**", "1536", "**44.4**", "**positive-definite**"],
    ],
    lead="Same data shape, weaker Poisson weights — only the significance changes.",
    widths=[0.26, 0.16, 0.24, 0.34],
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
    ],
    widths=[0.28, 0.15, 0.14, 0.13, 0.15, 0.15],
    note="The ratio is √(share × 18.7), independent of binning. What breaks at fine "
    "binning is sparsity, not noise: a zero bin asserts zero rate.",
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
    widths=[0.24, 0.30, 0.20, 0.26],
    note="The 115 MeV is the pion-as-kaon mass misassignment, ∝ (m_K²−m_π²)/2p.",
)

d.figure(
    "The A/e/M variation, all four η groups",
    "variation_e_eta.png",
    lead="The e variation looks like a mass scale should, and grows toward low pT.",
)

d.figure(
    "The composition is right; the shape was wrong",
    "bkg_composition.png",
)

d.figure(
    "2018's kaon starts at 1.2 GeV and runs softer",
    "cmp1618_kaonpt_cvh.png",
)

d.figure(
    "Grouped impacts on the joint mZ + J/ψ fit",
    "joint_s18.png",
    lead="Leaving rename_groups empty is what recovers the per-parameter rows.",
)

d.bullets(
    "Things found the hard way, worth not repeating",
    [
        "A fit that writes output has not converged — grep `edmval`, and two-sided, "
        "because −0.0 passes a one-sided test",
        "`genBPdgId` is the nearest LAST-COPY b hadron, so B*⁺ beats its own B⁺ one time "
        "in four; requiring |pdgId| = 521 alone threw away 25% of signal",
        "The nano's impact parameter was referenced to the origin, not the PV — shape "
        "χ² per bin 2340 → 3.92 once fixed",
        "The rabbit submodule gets reset and eats our param models; the source of truth "
        "lives in `scripts/rabbit/patches/` with an idempotent installer",
        "The mass window is a positional extra, not an env var — one pair of fits ran on "
        "a 10-bin axis before that was noticed",
    ],
)

print(d.save("260828_btojpsik_overview.pptx"))
