# Change: Lower bachelor pT and add a Kalman-vertex-fit preset for the non-V0 J/ψ+X channels

## Why

`add-jpsi-x-selection-comparison` (Phase-1/2 complete, decision locked at preset B for the four non-V0 channels) left two problems unresolved.

**Problem 1 — preset B's bachelor pT is too high for an alignment producer.** The locked preset B has `minBachelorPt=1.5 GeV` on B+/Bc and `minDaughterPt=1.0 GeV` on K*0/φ. The whole point of the J/ψ+X stream is to feed low-pT tracks into the tracker-alignment fit. The central V0Producer clone we route Ks/Λ through (`ALCARECOTkAlV0Candidates`) already runs at `tkPtCut=0.1 GeV` — the same floor used by the sibling `TkAlKsToPiPi` and `TkAlLambdaToProtonPi` streams. The non-V0 channels currently sit a factor of 15 above that floor. This is the asymmetry the alignment use case asks us to fix.

**Problem 2 — the only reason preset B survives at 1.5 GeV is that it has no topological cleanup.** Preset A (kinematic-loose) blows up to 13–125 cands/event on data; preset B trims that to O(1) using bachelor-pT and bachelor-to-muon DOCA cuts alone. Lowering the bachelor pT to 0.1 GeV reopens the combinatorial floodgates: the suppressed-by-pT background returns, and without a topological handle the AlCaReco bytes/event budget breaks. Phase 2 results (`add-jpsi-x-selection-comparison/results.md`) show preset A at 12× preset B's output size driven mostly by B0→K*0 daughter combinatorics — exactly the regime we re-enter when we drop the bachelor pT.

The natural answer is a third preset, **C**, that lowers the bachelor pT and recovers the cands/event budget through a Kalman vertex fit. The 2026-06 cost decision ("no AlCaReco-stage Kalman vertex fitting") in the sibling proposal is **partially reversed** under preset C only — and only on the four non-V0 channels — for the explicit, quantitative reason that lowering bachelor pT without a topological handle is structurally unworkable on data. Preset A and the updated preset B remain Kalman-free; preset C is the alignment-grade alternative if the budget-and-yield gates clear.

## What changes

Three things, in order of dependency.

### 1. Lower the bachelor / daughter pT floor in preset B itself

Preset B is updated **in place** (not forked) so the production default after the sibling change archives lands at the alignment-aligned pT floor. Specifically:

| Channel | Parameter | Old preset B | New preset B |
|---|---|---:|---:|
| B+, Bc | `minBachelorPt` | 1.5 | **0.1** |
| B0→K*0 | `Kstar.minDaughterPt` | 1.0 | **0.1** |
| Bs→φ | `Phi.minDaughterPt` | 1.0 | **0.1** |
| (all) | `maxBachelorEta`, `maxBachelorMuTrackDOCA`, `minJpsiPt`, `minMotherPt` | unchanged | unchanged |

No other parameter moves. The η, DOCA, J/ψ pT, and mother pT cuts remain at their Phase-1-locked values. The V0-mode channels remain preset-invariant (and continue to consume the existing `ALCARECOTkAlV0Candidates` clone at `tkPtCut=0.1`).

This delivers Problem 1's fix immediately, *before* preset C is evaluated. If the updated preset B passes the Phase-2 gates on its own (cands/event ≤ 5, tight-window fraction ≥ 30%, tight-window yield ≥ 100 / 1000 events), preset C becomes optional: bachelor pT is now where alignment wants it, and the topological cleanup may be unnecessary. Phase 2 of this change is the empirical test.

### 2. Add preset C: updated preset B + Kalman vertex-fit topology on non-V0 channels

Preset C inherits the updated preset B's kinematic + geometric cuts unchanged, and adds a single B-level Kalman vertex fit on each of the four non-V0 producers:

| Channel | Vertex-fit topology | Cuts at fit stage |
|---|---|---|
| B+ → J/ψ K+ | **3-body** Kalman fit on (μ⁺, μ⁻, K⁺) at a common B+ vertex | `vtxProb > 0.01`, `α_BS < acos(0.99) ≈ 0.142 rad` (3D angle from beamspot to fitted B vertex), `Lxy/σ > 3` (from first offlinePrimaryVertices entry) |
| Bc → J/ψ π | **3-body** Kalman fit on (μ⁺, μ⁻, π) at a common Bc vertex (same topology as B+, with pion mass in place of kaon mass per user clarification 2026-06-11) | Same as B+ |
| B0 → J/ψ K*0(→K⁺π⁻) | **4-body** Kalman fit on all leaf tracks (μ⁺, μ⁻, K, π) at a common B0 vertex | Same as B+ |
| Bs → J/ψ φ(→K⁺K⁻) | **4-body** Kalman fit on all leaf tracks (μ⁺, μ⁻, K⁺, K⁻) at a common Bs vertex | Same as B+ |

**Implementation note — single fit per channel, not two 2-body fits**: an earlier discussion (chat 2026-06-11) proposed two separate 2-body fits for B0→K*0 and Bs→φ (μμ for J/ψ + (K,π)/(K,K) for sub-resonance). The existing C++ in `JpsiXCandidateProducer.cc` already implements a single multi-body Kalman fit on all leaf tracks (the `applyBLevelVertexFit()` method, line 299). This change uses the existing implementation rather than rewriting it for two 2-body fits, on three grounds: (i) the C++ is compiled, tested, and unchanged in this proposal — zero regression risk; (ii) the single multi-body fit is the standard B-physics pattern and gives a B vertex that the alphaBS and Lxy/σ cuts are naturally interpreted against; (iii) the two-2-body approach loses the cross-fit constraint that the 4 tracks share a common parent, which is the strongest topological discriminant against combinatorial background. The K*0 / φ sub-resonance `TwoBodyDecayCandidateProducer.applyVertexFit` flag stays `False` under all three presets (the B-level fit already constrains the sub-resonance tracks).

**Bc topology (resolved 2026-06-11 in chat)**: Bc → J/ψ π is treated identically to B+ → J/ψ K, with pion mass in place of kaon mass. No additional code path needed — the C++ track-mode handles both 3-body cases. Bc is **no longer deferred** under preset C.

**Bc Kalman cut values — identical to B+ for first Phase-2 pass, follow-up if signal-acceptance-limited.** Physical differences between Bc and B+ (Bc cτ ≈ 137 μm vs B+ cτ ≈ 491 μm; Bc lab-frame Lxy ~25% of B+ at sample pT) make the same `Lxy/σ > 3` cut asymmetric: P(Lxy > 90 μm | Bc) ≈ 66% vs 91% for B+ (25 pp signal-acceptance loss). Similarly, at Bc's shorter Lxy the per-event displacement-direction angular resolution is wider, so cos α > 0.99 pinches Bc harder. We **keep the cuts identical** for the first Phase-2 run on three grounds: (i) the alignment producer optimizes for track quality of surviving candidates, not signal efficiency; (ii) we have very little Bc MC, so tuning Bc cuts separately at this stage is empirically poorly-grounded; (iii) if Phase-2 data shows Bc yields are acceptance-limited (e.g. < 50 / 1000 events on Run2016H), the obvious follow-up knobs are `minBLxyOverSigma` (3.0 → 2.0 recovers ~10 pp acceptance) and `maxMotherAlphaBS` (acos(0.99) → acos(0.98) ≈ 0.20 rad recovers ~5 pp). Recorded in `tasks.md` § 7 as a Phase-2 iteration trigger.

**No C++ changes in this proposal**: the optional Kalman parameters (`minBVtxProb`, `maxMotherAlphaBS`, `minBLxyOverSigma`) already exist on `JpsiXCandidateProducer` with no-op defaults. Conditional ESConsume of `TransientTrackBuilder` is already gated on `applyBVtxFit_` (true iff any of the three is set). The C++ producer was authored with this preset-C use case in mind (the source's "Phase 3" docstring section, lines 47-52, anticipates this activation).

Preset C does **not** touch V0-mode channels (B0→Ks, Λb, ψ(2S)). They continue to ride the central V0Producer (clone at `tkPtCut=0.1`) — no additional Kalman work at AlCaReco stage on V0 modes.

### 3. Phase-0 MC tutorial figures on B+ MiniAOD (standalone FWLite)

Before any AlCaReco run, produce the falling-exponential → mass-peak demonstration plots that motivate the whole stream. NanoAOD is unsuitable as the no-cut baseline (BPark Nano pre-filters on dimuon and track quality); MiniAOD has all `slimmedMuons` + `packedPFCandidates` + `lostTracks` unselected. The standalone script reads the existing `BuToJpsiK_BMuonFilter` UL18 MINIAODSIM dataset (`mc_minilist.txt`, 437 files; one file suffices for a tutorial figure) over xrootd.

Stages plotted (one canvas, four overlays):
1. **No cuts** — every (μ⁺μ⁻, charged-hadron) triplet from the event; expected shape: falling combinatorial under the B+ mass window.
2. **+ minimal kinematic** — muons soft-ID + pT > 2 GeV, track pT > 0.1 GeV, |η| < 2.4, dimuon mass within ±200 MeV of J/ψ.
3. **+ updated preset B kinematic** — bachelor pT > 0.1, |η_K| < 1.8, DOCA-to-muon < 0.03 cm (computed by hand from straight-line track-track DCA; no Kalman fit here).
4. **+ preset C** — preset B + a Kalman vertex fit on (μμK), cuts on vtxProb / cos(α) / Lxy/σ as listed above.

Output goes to (a) the worklog deck `slides/jpsix-selection-comparison-worklog.tex` as a chronological record, (b) `~/public_html/mz/alcareco/btojpsik-mass-demo/` for shareable browsing. Other channels (B0→K*0, Bs→φ, Bc) are out of scope for Phase 0 — extending the demo will need separate MC samples and is recorded as a follow-up.

The same script is extended in Phase 1 (Studies section) with gen-matching + a `gates.json` emitter to evaluate the Phase-1 gates on MiniAOD directly — no NanoAOD histmaker involvement. Running both Phase 0 and Phase 1 from one FWLite tool ensures the cut surface in Phase 1 matches what preset C will do at the AlCaReco stage in Phase 2.

## Decision criterion

Two phases, mirroring the sibling change.

**Phase 1 (NanoAOD MC, B+ only)** — re-evaluate the existing histmaker gates under both the updated preset B and preset C:
- signal efficiency on gen-matched B+ ≥ 70% relative to current preset B (i.e. lowering the pT floor must not throw away signal — it should *gain* it);
- per-event fake rate under (updated B vs current B): the bachelor-pT relaxation will *increase* fakes; quantify by how much. Preset C must restore fake rate to within 1.5× current preset B's value (the Kalman fit's job);
- mass-shape width and mean preserved across all three (current B, updated B, C).

**Phase 2 (Run2016H Charmonium RAW, all non-V0 channels)** — the original three gates from the sibling change, applied to updated B and to C:
- cands/event ≤ 5 in the channel's mass window;
- tight-window fraction ≥ 30%;
- tight-window absolute yield ≥ 100 / 1000 events.

Three outcomes are possible after Phase 2:
- **(i) updated B passes all gates on all four non-V0 channels** → preset C is unnecessary; updated B becomes the production preset; preset C is recorded as "evaluated, not needed" in `results.md` and the C++ vertex-fit hooks are retained (no-op by default) for future use.
- **(ii) updated B fails on some channels; preset C passes on those** → per-channel default in the cff: each channel picks the cheapest passing preset.
- **(iii) Both fail on a channel** → flag for Stage-2 cleanup, mirror the sibling change's handling of currently-failing channels.

The cands/event-vs-CPU tradeoff is reported alongside the gate evaluation. Preset C's KalmanVertexFitter cost (estimated O(1 ms) per surviving candidate) is measured via FastTimerService.

## Studies

The phases are sequenced strictly because each gates the next.

### Phase 0 — MC tutorial figures (standalone script, no cmsRun)

`scripts/btojpsik/mini_mass_demo.py` — new, FWLite-only, no narf/RDataFrame. Reads one MiniAOD file from `mc_minilist.txt` via xrootd, runs over `--max-events 5000` (configurable), emits one PNG with the four-stage overlay described above plus a `summary.json` with per-stage candidate counts. The Kalman fit at stage (4) uses FWLite-accessible `TransientTrackBuilder` and `KalmanVertexFitter` via `DataFormats.RecoVertex` — the same modules `JpsiXCandidateProducer` will call in preset C, so what stage (4) shows in MC is exactly what preset C will do on data.

### Phase 1 — MiniAOD gate evaluation via the same FWLite script

Phase 1 reuses the same FWLite script introduced in Phase 0 (`scripts/btojpsik/mini_mass_demo.py`) rather than extending the NanoAOD histmaker. NanoAOD's `bkmm_jpsimc_*` branches pass through BPark / BMM-Tools upstream selections — a track pre-filter, a dimuon mass-window, and (for some variants) an implicit vertex-quality requirement. The resulting "no-cut" baseline on Nano is a low-cut baseline, inflating apparent signal efficiency and deflating apparent fake rate. Cut values tuned against those biased numbers under-predict the work preset C's Kalman fit needs to do on real ALCAReco, and we wouldn't catch the mismatch until Phase 2 — exactly the regime where iteration is expensive.

Phase 1 therefore extends `mini_mass_demo.py` with: gen-matching against `prunedGenParticles` (B+ → J/ψ(μμ) K+ at PDG IDs 521 / 443 / 321; ΔR < 0.03 with charge alignment); a multi-preset loop that runs A / updated B / C in a single pass over the event loop; and a `gates.json` emitter for the three Phase-1 gates. The vertex quality at preset C is computed from a **closed-form 3-track common-vertex estimate** (analytic least-squares of three helix points-of-closest-approach, χ² from impact-parameter sigmas). This is not bit-identical to the iterative `KalmanVertexFitter` in `JpsiXCandidateProducer`, but operates on the same MiniAOD tracks and the same `offlineBeamSpot` — a much closer proxy than NanoAOD BMM-Tools branches. A small offset in cut thresholds between Phase 1 and Phase 2 is expected; the Phase-2 ALCAReco pass is the definitive measurement and one iteration cycle on Phase-2 results is budgeted (`tasks.md` § 7). The rationale for not using the full Kalman in Phase 1 — pure-Python FWLite cannot easily access `MagneticField` from `EventSetup` in CMSSW 10.6 — is documented in `design.md`.

Iteration on Kalman cut values is cheap: the gen-matched candidate pool is cached in memory across cut-value sweeps, so a full vtxProb × cos(α_xy) × Lxy/σ grid scan is one pass over the file. Each iteration produces a worklog frame.

The histmaker's existing `--alcarecoPreset {A,B}` from the sibling change is **not** extended to C in this proposal — the preset-comparison evaluation lives in the FWLite tool where the cut surface matches Phase 2 exactly. The histmaker continues to serve its analysis-level role unchanged.

### Phase 2 — ALCAReco confirmation on Run2016H Charmonium RAW

Same pattern as the sibling change's Phase 2: cmsDriver three configs (`TKALJPSIX_SELECTION_PRESET={A,B,C}`), cmsRun each on the same 1000-event sample, run `jpsix_preset_compare.py` (extended to accept three preset inputs instead of two), and emit:
- 2a — per-channel mass quality and gate evaluation;
- 2b — bytes/event + CPU-per-event from FastTimerService;
- 2c — dedup compression (already implemented from sibling change).

Preset A is kept in the comparison only as a reference baseline (structurally unusable, established by the sibling change); the operationally meaningful comparison is updated B vs C.

### Slide deck and public-html record

Same worklog discipline as the sibling change: every iteration, anomaly, and decision produces a frame in `slides/jpsix-selection-comparison-worklog.tex` (extending the existing deck, not a new one). The Phase-0 MC figures are copied to `~/public_html/mz/alcareco/btojpsik-mass-demo/` immediately on production; the Phase-1 + Phase-2 numbers are summarized into one headline frame in the public-facing `slides/jpsix-alcareco-producer.tex` after the gate decision is final.

## Resolved questions

The four open design questions before drafting:

1. **MC tier for the mass-spectrum demo** — MiniAOD, not NanoAOD. (Resolved 2026-06-11 in chat; rationale: BPark Nano is pre-filtered, so the no-cut falling-exp baseline is unavailable.)
2. **Sample for Phase 0** — `BuToJpsiK_BMuonFilter` UL18 MINIAODSIM, list at `mc_minilist.txt`. B+ only at Phase 0; other channels deferred until matching MC is located.
3. **Bachelor pT floor under updated preset B** — 0.1 GeV (matches existing `ALCARECOTkAlV0Candidates` clone's `tkPtCut`).
4. **Vertex-fit topology on B0→K*0 and Bs→φ** — two 2-body fits (μμ + sub-resonance), no subsequent 4-body refit. (Resolved 2026-06-11 in chat.)
5. **Phase-1 evaluation tier** — MiniAOD via the standalone Phase-0 script, not NanoAOD via the histmaker. (Resolved 2026-06-11 in chat; rationale: NanoAOD BPark pre-filtering biases the "no-cut" baseline, making gate evaluation untrustworthy; running Phase 1 on MiniAOD via the same FWLite path as Phase 0 + Phase 2 cmsRun eliminates the cut-surface mismatch entirely.)

## Remaining open items (non-blocking)

- Bc topology under preset C (3-body for Bc → J/ψ π; multi-body for Bc → J/ψ 3π / μν; user-requested follow-up).
- Sister-channel MC for B0→K*0, Bs→φ, Bc to extend the Phase-0 mass-demo gallery.
- CPU-budget sign-off for preset C from the alignment review (Kalman fit per surviving candidate is the cost we're adding).
- Whether the FWLite Kalman fit in `mini_mass_demo.py` should use `KalmanVertexFitter` (geometric vertex only) or `KinematicConstrainedVertexFitter` (mass-constrained on the J/ψ); preset C in cmsRun will use the same choice, so they should match. Default: unconstrained `KalmanVertexFitter` for both, on cost grounds.

## Impact

- **Affected specs**: `alcareco-jpsi-x` (MODIFIES the Selection-Preset Switching and No-Kalman-Vertex-Fitting requirements; ADDS new Kalman-vertex-fit-under-preset-C and Phase-0-MC-tutorial requirements).
- **Affected code**:
  - `Alignment/CommonAlignmentProducer/python/ALCARECOTkAlJpsiX_cff.py` — update preset B `minBachelorPt` / `minDaughterPt` to 0.1 GeV; add `'C'` entry to `_NON_V0_PRESETS`; extend env-var validation to accept `'A','B','C'`; under preset C only, inject `minBVtxProb`, `maxMotherAlphaBS`, `minBLxyOverSigma` onto the four non-V0 `JpsiXCandidateProducer` instances.
  - `Alignment/CommonAlignmentProducer/plugins/JpsiXCandidateProducer.cc` — **no changes** in this proposal. The Kalman vertex-fit branch already exists (method `applyBLevelVertexFit()`, line 299; activated via existence-check on `minBVtxProb` / `maxMotherAlphaBS` / `minBLxyOverSigma` in the cms.PSet). No new ESConsume needed (`TransientTrackBuilder` already conditionally consumed via `applyBVtxFit_`).
  - `CMSSW_10_6_17_patch1/jpsix_preset_compare.py` — extend to three-preset mode.
  - `scripts/btojpsik/mini_mass_demo.py` — new standalone FWLite script serving **both** Phase 0 (four-stage visual demo) and Phase 1 (gen-matched gate evaluation with multi-preset Kalman). No changes to `scripts/histmakers/btojpsik.py` or `wremnants/production/btojpsik_selections.py` — the histmaker is not extended in this change.
  - `scripts/btojpsik/figs/` + `~/public_html/mz/alcareco/btojpsik-mass-demo/` — figure output dirs.
  - `slides/jpsix-selection-comparison-worklog.tex` — extended worklog frames.
- **Dependency on sibling change**: `add-jpsi-x-selection-comparison` introduces the preset-switching mechanism, `JpsiXCandidateProducer` parameter set, `_NON_V0_PRESETS` dict, and the Phase-2 gate framework. This change builds on top of all of them. If the sibling has not archived before this change is applied, the deltas here MUST be reviewed for compatibility with the sibling's final form.
- **Out of scope**: any change to V0-mode channels; any Bc preset-C topology; any Stage-2 (CVH refit) work; any frontend / fit-stage code (histmaker / tensor / rabbit fit) — the histmaker is not touched in this change.
