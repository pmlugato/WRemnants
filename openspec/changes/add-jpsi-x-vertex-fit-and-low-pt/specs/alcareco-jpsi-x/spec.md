## MODIFIED Requirements

### Requirement: Selection-Preset Switching (non-V0 channels only)

`ALCARECOTkAlJpsiX_cff.py` SHALL expose a single selection-preset switch (`TKALJPSIX_SELECTION_PRESET` environment variable, with the values `'A'`, `'B'`, or `'C'`; default `'B'`) that determines the cut parameters of the four non-V0-mode `JpsiXCandidateProducer` instances (B+, B0→K*0, Bs→φ, Bc) plus the `applyVertexFit` / `minVtxProb` / `minDaughterPt` values on the K*0 and φ `TwoBodyDecayCandidateProducer` instances. The three V0-mode `JpsiXCandidateProducer` instances (B0→Ks, Λb, ψ(2S)) SHALL have their cuts fixed at "V0 quality + minimal kinematic" (J/ψ pT, αBS shared with the rest of the stream; no DCA; no Kalman vertex-fit cuts) regardless of the preset value — the V0Producer's Kalman vertex fit on the Ks/Λ sub-resonance (via the `ALCARECOTkAlV0Candidates` local clone running at `tkPtCut=0.1`) already provides clean candidates at ~1 / event without any further AlCaReco-stage selection. Mass windows for every channel SHALL NOT depend on the preset (mass windows define the channels). All non-mass cuts on non-V0 channels (J/ψ pT, mother pT, bachelor pT, bachelor η, bachelor-to-J/ψ-vertex DCA, αBS-at-J/ψ, plus the preset-C-only Kalman vertex-fit cuts) SHALL be derived from one of three preset profiles defined in a single `_NON_V0_PRESETS` dict in the cff. A single CMSSW build SHALL be sufficient to run all three presets without recompilation. The env var is consumed at cff-import time; since `ALCARECOTkAlJpsiX_cff.py` is transitively loaded by `Configuration.StandardSequences.AlCaRecoStreams_cff` inside the generated cmsRun config, the env var MUST be set both at cmsDriver invocation AND at cmsRun invocation, per the bug documented in the sibling change (2026-06-11).

#### Scenario: Preset B is the default for non-V0 channels
- **WHEN** `cmsRun` is invoked without `TKALJPSIX_SELECTION_PRESET` set in the environment
- **THEN** the cuts on every non-V0 per-channel producer are the **updated** Phase-1+2 values (bachelor pT 0.1 GeV; daughter pT 0.1 GeV on K*0 and φ; other kinematic / DOCA cuts unchanged from the sibling change), and the K*0 and φ `TwoBodyDecayCandidateProducer` instances have `applyVertexFit=False`

#### Scenario: Preset A disables non-mass cuts on non-V0 channels
- **WHEN** `TKALJPSIX_SELECTION_PRESET=A` is set
- **THEN** every kinematic cut (`minJpsiPt`, `minMotherPt`, `maxBachelorEta`, `maxBachelorMuTrackDOCA`) on the four non-V0 `JpsiXCandidateProducer` instances is configured to its no-op value (0 or +∞ as appropriate). The bachelor-pT cut on B+ and Bc producers remains at the HLT-floor value (0.5 / 0.3 GeV) by design. The K*0 and φ producers have `applyVertexFit=False`

#### Scenario: Preset C inherits updated preset B kinematics and adds Kalman vertex-fit cuts on all four non-V0 channels
- **WHEN** `TKALJPSIX_SELECTION_PRESET=C` is set
- **THEN** every kinematic + geometric cut value on the four non-V0 producers is bit-identical to its updated preset B value, AND on all four non-V0 producers (`ALCARECOTkAlJpsiX{BPlus,Bc,B0Kstar,BsPhi}Candidates`) the three Kalman parameters are set with values `minBVtxProb=0.01`, `maxMotherAlphaBS=acos(0.99)≈0.1415 rad`, `minBLxyOverSigma=3.0`. Bc receives the same parameters as B+ (Bc → J/ψ π is a 3-body track-mode decay analogous to B+ → J/ψ K). The K*0 and φ `TwoBodyDecayCandidateProducer` instances retain `applyVertexFit=False` (the B-level multi-body fit in `JpsiXCandidateProducer` already constrains the sub-resonance tracks)

#### Scenario: V0-mode channels are preset-invariant across all three presets
- **WHEN** the preset is changed between A, B, and C
- **THEN** the cut parameters of `ALCARECOTkAlJpsiXB0KsCandidates`, `ALCARECOTkAlJpsiXLambdabCandidates`, and `ALCARECOTkAlJpsiXPsi2SCandidates` remain bit-identical, and the candidates emitted by these three producers on the same input are bit-identical across all three preset settings. The `ALCARECOTkAlV0Candidates` clone parameters (`tkPtCut=0.1` etc.) are likewise preset-invariant

#### Scenario: Env var must be set on cmsRun (not just cmsDriver)
- **WHEN** the env var is set during cmsDriver config generation but not during the subsequent cmsRun invocation
- **THEN** the cmsRun job consumes the default preset because the cff is re-imported by the generated config's `process.load('Configuration.StandardSequences.AlCaRecoStreams_cff')` at cmsRun startup. Operators MUST set `TKALJPSIX_SELECTION_PRESET=<preset>` on both invocations of cmsDriver and cmsRun for a given preset run

---

### Requirement: No Kalman Vertex Fitting at the AlCaReco Stage

The AlCaReco-stage producers in `ALCARECOTkAlJpsiX_cff.py` SHALL NOT invoke any Kalman vertex fit **under presets A or B**. Specifically, under presets A and B, `applyVertexFit` on the K*0 and φ `TwoBodyDecayCandidateProducer` instances SHALL be `False`, and `JpsiXCandidateProducer` SHALL run with `minBVtxProb`, `minCosAlphaBS`, `minBLxyOverSigma` at their no-op default values (so the producer does NOT ESConsume `TransientTrackBuilder` or `MagneticField` and the Kalman fit branch is not entered). **Under preset C only**, Kalman vertex fitting IS enabled on the four non-V0 channels per the separate "Kalman Vertex-Fit Topology Under Preset C" requirement below. The V0Producer's Kalman fit on Ks/Λ via `ALCARECOTkAlV0Candidates` is exempt under all presets — it runs centrally in standard RECO (the AlCaReco stream consumes its output, not adds to it) and the local clone at `tkPtCut=0.1` is preset-invariant.

#### Scenario: K*0 sub-resonance vertex fit disabled under preset A or B
- **WHEN** an event passes the J/ψ + K*0 + B0 candidate chain under preset A or B
- **THEN** `ALCARECOTkAlJpsiXKstarCandidates` emits all opposite-sign track pairs whose K±π∓ mass falls in the K*0 mass window with NO vertex-probability cut applied, and the stored candidate vertex is the midpoint of the two track reference points (not a Kalman-fitted vertex)

#### Scenario: No B-level Kalman fit invoked under preset A or B
- **WHEN** any of the seven `JpsiXCandidateProducer` instances runs under preset A or B
- **THEN** it does NOT ESConsume `TransientTrackBuilder`, does NOT call `KalmanVertexFitter`, and the cuts on each candidate are exhaustively the kinematic + geometric cuts in the preset profile (mass windows, pT, η, αBS, DCA)

## ADDED Requirements

### Requirement: B-Level Kalman Vertex Fit Under Preset C

Under preset C, each of the four non-V0 `JpsiXCandidateProducer` instances (B+, Bc, B0→K*0, Bs→φ) SHALL invoke a single B-level Kalman vertex fit on all of its candidate's leaf tracks at a common parent vertex. The fit is performed by the producer's existing `applyBLevelVertexFit()` method (`JpsiXCandidateProducer.cc` line 299), which is activated when ANY of the three parameters `minBVtxProb`, `maxMotherAlphaBS`, `minBLxyOverSigma` exists in the cms.PSet. No C++ changes are required by this proposal — the Kalman branch is wired and tested; under presets A and B those parameters are simply absent from the cms.PSet and the branch is inactive.

The K*0 and φ `TwoBodyDecayCandidateProducer.applyVertexFit` flag SHALL remain `False` under all three presets — the B-level multi-body fit in `JpsiXCandidateProducer` already constrains the sub-resonance tracks; a separate sub-resonance fit is redundant and is not introduced in this change.

**Topology per channel under preset C** (single fit per candidate, not two separate fits):

| Channel | Producer mode | Tracks in fit | Common vertex |
|---|---|---|---|
| B+ → J/ψ K+ | track | (μ⁺, μ⁻, K⁺) | B+ |
| Bc → J/ψ π | track | (μ⁺, μ⁻, π) | Bc |
| B0 → J/ψ K*0(→K⁺π⁻) | vcc | (μ⁺, μ⁻, K, π) | B0 |
| Bs → J/ψ φ(→K⁺K⁻) | vcc | (μ⁺, μ⁻, K⁺, K⁻) | Bs |

**Cuts** (same values for all four channels):

- `minBVtxProb = 0.01` — Kalman fit probability from χ²/ndf
- `maxMotherAlphaBS = acos(0.99) ≈ 0.1415 rad` — max 3D angle between mother p⃗ and the (beamspot → fitted-vertex) direction; equivalent to `cos α ≥ 0.99`
- `minBLxyOverSigma = 3.0` — transverse flight significance from the first `offlinePrimaryVertices` entry

**ESConsume**: `JpsiXCandidateProducer` ESConsumes `TransientTrackBuilder` via `TransientTrackRecord` only when at least one of the three Kalman parameters exists in its cms.PSet. Under presets A/B those parameters are absent and the producer does not consume the record.

#### Scenario: B+ 3-body Kalman fit applied under preset C
- **WHEN** under preset C a (J/ψ, K⁺) combination passes all kinematic + geometric cuts
- **THEN** `JpsiXCandidateProducer` collects two muon track refs and the bachelor track ref, builds three `reco::TransientTrack` instances via the consumed `TransientTrackBuilder`, calls `KalmanVertexFitter::vertex(tracks)`, and rejects the candidate if the fit is invalid OR `vtxProb < 0.01` OR α_BS (3D angle, beamspot to fitted vertex) > `acos(0.99)` OR `Lxy/σ < 3.0` (using the first `offlinePrimaryVertices` entry)

#### Scenario: Bc 3-body Kalman fit applied under preset C
- **WHEN** under preset C a (J/ψ, π⁺) combination passes all kinematic + geometric cuts
- **THEN** the same fit + cut sequence as B+ runs, with the pion track in place of the kaon track; Bc is not deferred under preset C

#### Scenario: B0→K*0 4-body Kalman fit applied under preset C
- **WHEN** under preset C a (J/ψ, K*0) combination passes all kinematic + geometric cuts
- **THEN** `JpsiXCandidateProducer` (VCC mode) collects the J/ψ's two muon track refs AND the K*0's two daughter (K, π) track refs, builds four `reco::TransientTrack` instances, calls `KalmanVertexFitter::vertex(tracks)`, and applies the same three Kalman cuts as B+ to reject. No separate K*0 sub-resonance Kalman fit is performed (the K*0 producer keeps `applyVertexFit=False`)

#### Scenario: Bs→φ 4-body Kalman fit applied under preset C
- **WHEN** under preset C a (J/ψ, φ) combination passes all kinematic + geometric cuts
- **THEN** the same fit + cut sequence as B0→K*0 runs, with (K⁺, K⁻) in place of (K, π)

#### Scenario: K*0 / φ sub-resonance producers retain `applyVertexFit=False` under preset C
- **WHEN** the preset is C
- **THEN** `ALCARECOTkAlJpsiXKstarCandidates.applyVertexFit` and `ALCARECOTkAlJpsiXPhiCandidates.applyVertexFit` both remain `False`. The B-level fit in `JpsiXCandidateProducer` already constrains the sub-resonance daughter tracks at the B vertex; no additional sub-resonance Kalman fit is added

#### Scenario: Presets A and B do not activate the Kalman branch
- **WHEN** the preset is A or B
- **THEN** none of `minBVtxProb`, `maxMotherAlphaBS`, `minBLxyOverSigma` is set as a parameter on any `JpsiXCandidateProducer` instance in the cms.PSet (the C++ existence-check returns false), the producer's `applyBVtxFit_` flag is false, `TransientTrackBuilder` is not consumed, and the candidate selection on the same input is bit-identical to the sibling-change baseline (up to the lowered bachelor / daughter pT floor in preset B)

---

### Requirement: V0Producer Private Clone with Relaxed Track-pT Cut (Preset-Invariant)

The J/ψ+X stream SHALL consume V0 candidates (Ks, Λ) exclusively via the local `ALCARECOTkAlV0Candidates` clone of `generalV0Candidates`, which SHALL run with `tkPtCut=0.1` GeV (relaxed from the central default of 0.35 GeV) and otherwise inherit all other parameters from the central producer. This clone is **shared with the sibling streams** `TkAlKsToPiPi` and `TkAlLambdaToProtonPi` via standard CMSSW framework deduplication on `cms.InputTag` — within one cmsRun job, only one instance of `ALCARECOTkAlV0Candidates` is constructed regardless of how many AlCaReco streams reference it. The V0 clone parameters SHALL be preset-invariant — they are not driven by `_NON_V0_PRESETS` and do not change between presets A, B, and C.

The three V0-mode `JpsiXCandidateProducer` instances (B0→Ks, Λb, ψ(2S)) consume `ALCARECOTkAlV0Candidates:Kshort` and `ALCARECOTkAlV0Candidates:Lambda` as their X candidate input, exactly as the sibling change specifies. No additional AlCaReco-stage Kalman vertex fitting SHALL be performed on V0-mode channels under any preset.

#### Scenario: V0 clone is the sole V0 input for the stream
- **WHEN** any of the three V0-mode `JpsiXCandidateProducer` instances runs
- **THEN** its `xSrc` InputTag is `ALCARECOTkAlV0Candidates:Kshort` (B0→Ks, ψ(2S)) or `ALCARECOTkAlV0Candidates:Lambda` (Λb), NOT `generalV0Candidates:*`

#### Scenario: V0 clone parameters are preset-invariant
- **WHEN** the preset is changed across A, B, C
- **THEN** the `tkPtCut`, `tkChi2Cut`, `tkIPSigXYCut`, `tkIPSigZCut`, `vtxChi2Cut` (and any other parameter) of `ALCARECOTkAlV0Candidates` remain bit-identical, and the candidates it emits on a given input are bit-identical across all three preset settings

#### Scenario: Framework dedup across sibling streams
- **WHEN** a cmsRun job runs both `seqALCARECOTkAlJpsiX` and `seqALCARECOTkAlKsToPiPi` in the same process
- **THEN** only one `ALCARECOTkAlV0Candidates` cms.EDProducer is constructed; downstream modules in both streams consume from the same producer instance

---

### Requirement: Phase-0 Standalone MC Mass-Spectrum Demo (MiniAOD)

A standalone FWLite script `scripts/btojpsik/mini_mass_demo.py` SHALL be provided that reads one `BuToJpsiK_BMuonFilter` UL18 MINIAODSIM file (via xrootd or a local path), reconstructs B+ → J/ψ(μμ) K+ candidates from `slimmedMuons` + (`packedPFCandidates` ∪ `lostTracks`), and emits a single overlaid B-mass-distribution PNG with four selection stages plus a `summary.json` reporting per-stage candidate counts. The four stages SHALL be: (1) no cuts, (2) minimal kinematic, (3) updated preset B kinematic + DOCA, (4) updated preset B + vertex-quality cuts approximating preset C. The vertex computation at stage (4) SHALL be a closed-form 3-track common-vertex estimate (analytic least-squares of three helix points-of-closest-approach, with χ² from track impact-parameter sigmas) operating on MiniAOD tracks and the event's `offlineBeamSpot`; this is a Phase-1 approximation to the iterative `KalmanVertexFitter` used by preset C in cmsRun, with one iteration cycle on cut thresholds budgeted at Phase 2 (`tasks.md` § 7) to align Phase 1 and Phase 2. The stage-(4) cut values SHALL be `vtxChi2Prob > 0.01`, `cos(α_xy) > 0.99`, `Lxy/σ > 3`, matching the preset C cff values but interpreted on the closed-form vertex. The output figure SHALL be written to both `scripts/btojpsik/figs/btojpsik-mass-demo/B_mass_stages.png` AND `~/public_html/mz/alcareco/btojpsik-mass-demo/B_mass_stages.png` (the second copy with a one-line `index.html` listing the file and date). The script SHALL NOT depend on cmsRun, narf, or RDataFrame — it is pure FWLite + matplotlib.

#### Scenario: Stage 1 ("no cuts") shows the falling-exp combinatorial tail
- **WHEN** the script runs over ≥ 1000 events and stage (1) is plotted
- **THEN** the B mass distribution in [4.5, 6.0] GeV exhibits a smoothly falling shape across the full range with no statistically significant peak at the PDG B+ mass (5.279 GeV); the stage-1 candidate count is at least an order of magnitude higher than the stage-4 candidate count

#### Scenario: Stage 4 (preset C analogue) shows a clean B+ mass peak
- **WHEN** the script runs over ≥ 1000 events and stage (4) is plotted
- **THEN** the B mass distribution shows a statistically significant peak at the PDG B+ mass, the per-event candidate count after stage (4) is ≤ 5, and the tight-window (PDG ± 75 MeV) fraction is ≥ 30%

#### Scenario: Public-facing copy is updated on each run
- **WHEN** the script is invoked with `--public-html-dir ~/public_html/mz/alcareco/btojpsik-mass-demo/`
- **THEN** the output PNG is copied to that directory and an `index.html` containing one line (filename + ISO-8601 date) is written or appended
