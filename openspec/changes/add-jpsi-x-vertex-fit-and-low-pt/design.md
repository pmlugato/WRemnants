# Design: low-pT + vertex-fit preset C for the non-V0 J/ψ+X channels

## Context

The sibling change `add-jpsi-x-selection-comparison` deliberately ruled out Kalman vertex fitting at the AlCaReco stage on cost grounds (David's 2026-05 design review, slides 18–28). That decision was made with the assumption that the kinematic + bachelor-pT-floor cuts in preset B were sufficient to keep cands/event in budget on the four non-V0 channels.

Phase 2 of the sibling change confirmed the budget holds: preset B sits at O(0.05–0.5) cands/event vs preset A's O(10–125) cands/event, and the bytes/event ratio is 12×. But it holds *because of* the bachelor-pT floor at 1.5 GeV — exactly the cut that defeats the alignment use case the stream was designed for. The tracker-alignment fit needs the low-pT track regime (~0.1 GeV) where reconstruction biases are largest; cutting at 1.5 GeV throws away the diagnostic value of having J/ψ+X in the first place.

This change reopens the cost question with a quantitative target: if the only thing standing between us and a 0.1 GeV bachelor floor is the Kalman fit, we can pay the Kalman cost on the four non-V0 channels — they're the channels where the central V0Producer's Kalman fit is *not* doing the job for us. The V0-mode channels (B0→Ks, Λb, ψ(2S)) consume the central V0Producer (via the `ALCARECOTkAlV0Candidates` clone already running at `tkPtCut=0.1`); we do not add Kalman work there.

The original cost argument cited "every event in the parent dataset" as the reason. Two things have changed:
1. Preset B is now combinatorially sparse — most events have 0 or 1 non-V0 candidates surviving kinematic + DOCA cuts, so the Kalman fit runs on a *very small* per-event pool. Cost scales with surviving candidates, not events.
2. The V0Producer in central RECO is itself a Kalman fit, and the AlCaReco stream consumes ~10⁵ V0Producer fits per Run2016H file with no observed cost issue. Adding an analogous fit on a handful-per-event non-V0 candidates is a tiny incremental cost.

The proposal is therefore: lower bachelor pT in preset B *and* add preset C as a Kalman-fitted alternative; let Phase 2 measure the actual budget impact under both and decide per-channel.

## Architectural decisions

### Single multi-body B-level Kalman fit (not two 2-body fits)

The 2026-06-11 chat discussion proposed two separate 2-body fits for B0→K*0 and Bs→φ (μμ for J/ψ + sub-resonance daughter pair for K*0/φ), with quality cuts on each. On reviewing the existing `JpsiXCandidateProducer.cc`, we found the C++ already implements a single B-level Kalman fit (`applyBLevelVertexFit()` at line 299) that takes all leaf tracks (3 in track mode, 4 in VCC mode) and finds a common vertex. The proposal switches to using the existing implementation rather than rewriting for two-2-body fits, on three grounds:

- **Zero C++ regression risk.** The producer's compiled, tested behaviour under presets A/B is unchanged because the Kalman branch is gated on parameter existence in the cms.PSet, not on a value. Under preset C we add the parameters; under A/B we don't.
- **Standard B-physics pattern.** Single multi-body Kalman fit on all leaf tracks at a common parent vertex is the textbook B-physics approach. The αBS and Lxy/σ cuts are naturally interpreted against the parent (B) vertex, which is what the alignment use case wants to constrain.
- **Stronger topological discriminant.** The two-2-body approach loses the cross-fit constraint that all 4 tracks share a common B parent — exactly the constraint that suppresses combinatorial backgrounds where the K*0 and J/ψ are independent. A χ² on the 4-track fit penalizes mis-paired combinations far more aggressively than separate per-pair χ²'s would.

The remaining concern from the original discussion — that a failed multi-body fit drops the candidate entirely — is mitigated by the producer's behaviour: invalid fits return false from `applyBLevelVertexFit()`, which is the correct "reject this candidate" semantics for an alignment producer (we want clean candidates, not best-effort retries).

### Bc topology under preset C — same as B+

Bc → J/ψ π is a 3-body decay; with pion mass in place of kaon mass, the producer's `track` mode (used for both B+ and Bc) handles it identically. Per user clarification 2026-06-11, Bc gets the same Kalman parameters as B+ under preset C — no separate code path, no deferral.

### Bc cut values — same as B+ as starting point, asymmetry in Lxy/σ acceptance flagged

The physical asymmetry is real but small enough to absorb. Bc cτ ≈ 137 μm vs B+ cτ ≈ 491 μm; at typical sample pT, mean lab-frame Bc Lxy is ~25% of B+ Lxy (boost factor difference contributes a further ~10%, lifetime ratio dominates). With tracker σ_xy ≈ 30 μm, the same `Lxy/σ > 3` cut translates to:

- B+: Lxy > 90 μm acceptance is P(exp, mean 930 μm) ≈ 91%.
- Bc: same cut threshold, but P(exp, mean 220 μm) ≈ 66%.

A 25 percentage-point asymmetry. Similarly, `cos α > 0.99` requires angle < 0.14 rad, and Bc's shorter Lxy means the per-event displacement-direction resolution is wider (σ_dir ≈ σ_xy / Lxy ≈ 30/220 = 0.15 rad for Bc, vs 30/930 = 0.03 rad for B+). cos α > 0.99 sits at the ~1σ resolution mark for Bc and the ~3–5σ mark for B+.

Three reasons to keep the cuts identical for the first Phase-2 run anyway:

- **Alignment-grade ≠ signal-efficiency-optimal.** For the alignment fit, we care about *quality* of surviving candidates, not raw signal yield. Tighter cuts on Bc bias toward better-measured vertices, which the alignment fit prefers.
- **Insufficient Bc MC to tune separately.** BMuonFilter MC is signal-only for B+ → J/ψ K; Bc events appear only through generic combinatorial. Without dedicated Bc MC, we can't measure the signal-efficiency / fake-rate tradeoff for Bc-specific cuts.
- **Iteration is cheap.** If Phase-2 data shows Bc yields are acceptance-limited (e.g. < 50 candidates per 1000 events on Run2016H), the obvious knobs are `minBLxyOverSigma` (3.0 → 2.0 recovers ~10 pp Bc acceptance) and `maxMotherAlphaBS` (acos(0.99) → acos(0.98) ≈ 0.20 rad recovers ~5 pp). A Bc-specific preset C' is one ~5-line cff edit away.

The cff implements the symmetric starting point; the iteration trigger lives in `tasks.md` § 7.

Other Bc decay modes (J/ψ 3π, J/ψ μν, J/ψ τν) remain out of scope (the cff only exposes Bc → J/ψ π).

### V0-mode channels stay untouched

The central V0Producer (via the local clone at `ALCARECOTkAlV0Candidates`) already runs a Kalman fit on Ks/Λ with `tkPtCut=0.1`. Routing B0→Ks, Λb, ψ(2S) through that clone gets us alignment-grade output for free. No AlCaReco-stage Kalman work is added on V0-mode channels under any preset.

The sibling spec scenario "V0-mode channels are preset-invariant" therefore continues to hold across A, B, and C. This is preserved in the spec deltas.

### C++ parameter no-op defaults

The new `minBVtxProb`, `minCosAlphaBS`, `minBLxyOverSigma` parameters on `JpsiXCandidateProducer` default to `-1.0` (interpreted as "no Kalman fit; do not ESConsume `TransientTrackBuilder` / `MagneticField`"). This keeps presets A and B C++-behaviour-identical to the sibling-change baseline. The Kalman fit branch only activates when at least one of the three parameters is set to a non-default value. The cff under preset C is the only place that sets them.

The reason this matters: if preset C is found not worthwhile after Phase 2, the C++ parameters can stay (no-op) without affecting production. The wiring is reversible by config alone — no C++ revert needed if the alignment community decides the cost is too high.

### Where the Phase-1 Kalman fit lives

`scripts/histmakers/btojpsik.py` runs on NanoAOD via narf RDataFrame. NanoAOD does *not* carry track-level information sufficient to re-do a Kalman fit downstream; what NanoAOD has are the `BToJpsiK` derived quantities (`bkmm_jpsimc_*`) which already pass through BPark / BMM-Tools upstream selections (track pre-filter, dimuon mass-window, implicit vertex-quality requirement). Two consequences make NanoAOD unsuitable for Phase-1 gate evaluation:

1. **Baseline bias.** The "no-cut" baseline against which Phase-1 signal efficiency and fake-rate are measured is, on Nano, already a low-cut baseline. Apparent signal efficiency is inflated (denominator pre-cleaned); apparent fake rate is deflated (worst combinatorics already gone). Cut values tuned against these biased numbers under-predict the work preset C's Kalman fit needs to do on real ALCAReco.
2. **Cut-surface mismatch.** Even if we re-used `bkmm_jpsimc_bvtxProb`, `bkmm_jpsimc_bAlphaBS`, `bkmm_jpsimc_bLxyOverSigma` as Kalman proxies, those branches come from a different `KalmanVertexFitter` configuration (different track collection, different beamspot version) than what `JpsiXCandidateProducer` will run in cmsRun. Cut values measured on the proxy may shift between Phase 1 and Phase 2.

Decision: Phase 1 is run on MiniAOD via the same FWLite script as Phase 0 (`scripts/btojpsik/mini_mass_demo.py`), extended with gen-matching against `prunedGenParticles` and a multi-preset event-loop. The vertex quality at preset C is computed from a **closed-form 3-track common-vertex estimate** — analytic least-squares of three helix points-of-closest-approach to a candidate vertex, with χ² from track impact-parameter sigmas. This is **not bit-identical** to the iterative `KalmanVertexFitter` in `JpsiXCandidateProducer`, but it is a much closer proxy than NanoAOD's BMM-Tools branches: it operates on the same MiniAOD tracks, the same `offlineBeamSpot`, and produces the same three observables (vtxProb-equivalent from χ²/ndf, cos(α_xy), Lxy/σ). A small offset on cut thresholds between Phase 1 and Phase 2 is expected; the Phase-2 ALCAReco pass on real data is the definitive measurement and one iteration cycle on Phase-2 results is budgeted in `tasks.md` § 7 if needed.

Rationale for choosing the closed-form path over a full `KalmanVertexFitter` invocation: pure-Python FWLite cannot easily obtain `MagneticField` from `EventSetup` in CMSSW 10.6; getting the real Kalman would require writing a small C++ EDAnalyzer plugin and running via cmsRun, which is significant additional scope for a study whose decisive measurement is Phase 2 anyway. The closed-form vertex is exact in the limit of small helix curvature over the track-vertex distance (i.e. for low-mass parents and central-tracker decays, both true here) and degrades gracefully — the worst case is a small migration in cut values, which we measure and correct against Phase 2.

Cost of this choice: gen-matching in FWLite is ~30 lines (mother-daughter walk + ΔR-match); the closed-form vertex is ~50 lines; event-loop speed is single-threaded; 50k events × 3 presets × vertex computation completes in O(5 min) on one core — adequate for gate-tuning. The histmaker is not touched.

### Phase 0 sample choice

The user's data already has `BuToJpsiK_BMuonFilter` UL18 MINIAODSIM staged via xrootd. A B mass-spectrum demo on this sample shows:
- the falling-exp tail (no cuts);
- the kinematic-cut shape (preset B);
- the Kalman-fit shape (preset C — `KalmanVertexFitter` accessible via FWLite + `DataFormats.RecoVertex`).

We get the visual that motivates the whole stream and a side-by-side validation of preset C's Kalman cuts in MC — without writing C++ or running cmsRun. The MC sample is signal-only (no actual B+ → other-channel decay), so it cannot stand in for a multi-channel demo. The other channels' MC sourcing is deferred to a follow-up if/when needed for the public-facing deck.

## Tradeoffs

- **Kalman fit cost on data**: estimated at ~1 ms per candidate × ≤ 5 cands/event × 10⁵ events/file = O(10 minutes) per file. Acceptable but real. Measured precisely in Phase 2b via FastTimerService.
- **Lower bachelor pT under preset B without Kalman**: the trade. If updated preset B passes all four non-V0 gates on its own (Phase 2a), preset C is unnecessary and we save the Kalman cost. The decision is data-driven, not pre-committed.
- **C++ wiring complexity**: the conditional ESConsume + Kalman branch adds ~80 lines to `JpsiXCandidateProducer.cc`. The cost is paid once; the parameter no-op defaults keep the wiring invisible under presets A and B.
- **Single-threaded FWLite for Phase 1**: re-using the Phase-0 standalone script for gate evaluation gives up narf RDataFrame's multi-thread speedup. Cost is O(10 min) per cut-value sweep on one MiniAOD file — acceptable for a one-off tuning study, would matter if we ran Phase 1 over many MC datasets. We don't (B+ MiniAOD only).
- **Worklog deck length**: extending the sibling's worklog rather than starting new keeps the iteration history linear and readable. The deck will grow another ~10 frames; budget acceptable.

## Alternatives considered

1. **Make preset B's bachelor pT a runtime parameter and skip preset C**. Rejected: lowering the floor without topological cleanup breaks Phase-2 gates on B0→K*0 (combinatorial blowup). A floor parameter without a vertex fit doesn't solve the problem.
2. **Fork preset B into B_loose (pT > 0.5) and B_tight (pT > 1.5)**. Rejected: gives one knob too many, doesn't address the alignment use case (still no 0.1 GeV access).
3. **Add the Kalman fit at the Stage-2 CVH stage instead**. Rejected: Stage-2 takes cloned tracks as input and depends on AlCaReco-stage output being of usable quality. If AlCaReco outputs O(100 cands/event) on B0→K*0 because we lowered pT without cleanup, Stage 2 inherits the explosion; the savings have to happen at AlCaReco.
4. **One full 4-body Kalman fit on B0→K*0 / Bs→φ instead of two 2-body fits**. Rejected for the locality / CPU / robustness reasons in "Architectural decisions" above.

## Out of scope

- V0-mode channel modifications (any).
- Bc-under-preset-C topology decision (deferred to follow-up).
- Stage-2 CVH refit work.
- Frontend / fit-stage code changes (histmaker tensor / rabbit fit) beyond the `--alcarecoPreset C` extension.
- Public-facing deck restructuring (only one new headline frame).
