## MODIFIED Requirements

### Requirement: Selection-Preset Switching (non-V0 channels only)

`ALCARECOTkAlJpsiX_cff.py` SHALL expose a single selection-preset switch (`TKALJPSIX_SELECTION_PRESET` environment variable, with the values `'A'`, `'B'`, or `'C'`; default `'B'`) that determines the cut parameters of the four non-V0-mode `JpsiXCandidateProducer` instances (B+, B0→K*0, Bs→φ, Bc) plus the `applyVertexFit` / `minVtxProb` / `minDaughterPt` values on the K*0 and φ `TwoBodyDecayCandidateProducer` instances. The three V0-mode `JpsiXCandidateProducer` instances (B0→Ks, Λb, ψ(2S)) SHALL have their cuts fixed at "V0 quality + minimal kinematic" (J/ψ pT, αBS shared with the rest of the stream; no DCA; no Kalman vertex-fit cuts; no J/ψ mass constraint) regardless of the preset value — the V0Producer's Kalman vertex fit on the Ks/Λ sub-resonance (via the `ALCARECOTkAlV0Candidates` local clone running at `tkPtCut=0.1`) already provides clean candidates at ~1 / event without any further AlCaReco-stage selection. Mass windows for every channel SHALL NOT depend on the preset (mass windows define the channels). All non-mass cuts on non-V0 channels (J/ψ pT, mother pT, bachelor pT, bachelor η, bachelor-to-J/ψ-vertex DCA, αBS-at-J/ψ, plus the preset-C-only Kalman vertex-fit cuts) SHALL be derived from one of three preset profiles defined in a single `_NON_V0_PRESETS` dict in the cff. The `applyJpsiMassConstraint` flag on every `JpsiXCandidateProducer` instance SHALL be derived from the same preset switch: `False` under presets A and B (raw track-sum dimuon p4 used for the mother mass), `True` under preset C only (the constraint is structurally required by the multi-track Kalman fit's `TwoTrackMassKinematicConstraint` and used as the dimuon-p4 fallback when the multi-track fit fails). A single CMSSW build SHALL be sufficient to run all three presets without recompilation. The env var is consumed at cff-import time; since `ALCARECOTkAlJpsiX_cff.py` is transitively loaded by `Configuration.StandardSequences.AlCaRecoStreams_cff` inside the generated cmsRun config, the env var MUST be set both at cmsDriver invocation AND at cmsRun invocation, per the bug documented in the sibling change (2026-06-11).

#### Scenario: Preset B is the default for non-V0 channels
- **WHEN** `cmsRun` is invoked without `TKALJPSIX_SELECTION_PRESET` set in the environment
- **THEN** the cuts on every non-V0 per-channel producer are the **updated** Phase-1+2 values (bachelor pT 0.1 GeV; daughter pT 0.1 GeV on K*0 and φ; other kinematic / DOCA cuts unchanged from the sibling change), the K*0 and φ `TwoBodyDecayCandidateProducer` instances have `applyVertexFit=False`, AND every `JpsiXCandidateProducer` instance has `applyJpsiMassConstraint=False` (raw track-sum dimuon p4 used in the mother mass computation)

#### Scenario: Preset A disables non-mass cuts on non-V0 channels
- **WHEN** `TKALJPSIX_SELECTION_PRESET=A` is set
- **THEN** every kinematic cut (`minJpsiPt`, `minMotherPt`, `maxBachelorEta`, `maxBachelorMuTrackDOCA`) on the four non-V0 `JpsiXCandidateProducer` instances is configured to its no-op value (0 or +∞ as appropriate), the K*0 and φ producers have `applyVertexFit=False`, AND every `JpsiXCandidateProducer` instance has `applyJpsiMassConstraint=False`. The bachelor-pT cut on B+ and Bc producers remains at the HLT-floor value (0.5 / 0.3 GeV) by design

#### Scenario: Preset C inherits updated preset B kinematics, adds Kalman vertex-fit cuts, and enables the J/ψ mass constraint
- **WHEN** `TKALJPSIX_SELECTION_PRESET=C` is set
- **THEN** every kinematic + geometric cut value on the four non-V0 producers is bit-identical to its updated preset B value, on all four non-V0 producers (`ALCARECOTkAlJpsiX{BPlus,Bc,B0Kstar,BsPhi}Candidates`) the three Kalman parameters are set with values `minBVtxProb=0.01`, `maxMotherAlphaBS=acos(0.99)≈0.1415 rad`, `minBLxyOverSigma=3.0`, AND every `JpsiXCandidateProducer` instance has `applyJpsiMassConstraint=True`. Bc receives the same parameters as B+. The K*0 and φ `TwoBodyDecayCandidateProducer` instances retain `applyVertexFit=False`

#### Scenario: V0-mode channels are preset-invariant across all three presets
- **WHEN** the preset is changed between A, B, and C
- **THEN** the cut parameters of `ALCARECOTkAlJpsiXB0KsCandidates`, `ALCARECOTkAlJpsiXLambdabCandidates`, and `ALCARECOTkAlJpsiXPsi2SCandidates` remain bit-identical (including `applyJpsiMassConstraint=False` under all three presets), and the candidates emitted by these three producers on the same input are bit-identical across all three preset settings. The `ALCARECOTkAlV0Candidates` clone parameters (`tkPtCut=0.1` etc.) are likewise preset-invariant

#### Scenario: Env var must be set on cmsRun (not just cmsDriver)
- **WHEN** the env var is set during cmsDriver config generation but not during the subsequent cmsRun invocation
- **THEN** the cmsRun job consumes the default preset because the cff is re-imported by the generated config's `process.load('Configuration.StandardSequences.AlCaRecoStreams_cff')` at cmsRun startup. Operators MUST set `TKALJPSIX_SELECTION_PRESET=<preset>` on both invocations of cmsDriver and cmsRun for a given preset run

---

### Requirement: Merged-Collection Track Filter pT Threshold

The final `AlignmentTrackSelectorWithIndexMap` instance (`ALCARECOTkAlJpsiX`) in `ALCARECOTkAlJpsiX_cff.py` SHALL apply a track pT threshold of `0.1 GeV` (`ptMin = 0.1`). This value SHALL match the soft-pT floor used by every upstream per-candidate threshold in the stream (`minBachelorPt = 0.1` on `JpsiXCandidateProducer`, `minDaughterPt = 0.1` on the `TwoBodyDecayCandidateProducer` instances, `tkPtCut = 0.1` on `ALCARECOTkAlV0Candidates`). Tracks below 0.1 GeV SHALL NOT appear in the saved `ALCARECOTkAlJpsiX` collection.

#### Scenario: V0 daughter pion survives all upstream cuts and is kept
- **WHEN** a Ks candidate's daughter pion has reconstructed `pT = 0.15 GeV` and passes every upstream filter (V0Producer kinematics, `JpsiXCandidateProducer` αBS, etc.)
- **THEN** that track survives the final `AlignmentTrackSelectorWithIndexMap` filter and is written to the output collection

#### Scenario: Very soft track is rejected
- **WHEN** a track has reconstructed `pT = 0.05 GeV`
- **THEN** the track is filtered out by `ALCARECOTkAlJpsiX` and does NOT appear in the output collection

---

### Requirement: J/ψ Mass Treatment — Window-Only Under Preset B

Under preset B, the J/ψ → μμ candidates produced by `TwoBodyDecayCandidateProducer` (`ALCARECOTkAlJpsiXJpsiCandidates`) SHALL be selected solely by the dimuon mass window (`minMass = 2.95 GeV`, `maxMass = 3.25 GeV`), opposite-sign requirement, and the upstream `TkAlGoodIdMuonSelector` muon ID. No `MassKinematicConstraint` SHALL be applied to the dimuon four-momentum used in any downstream mother mass calculation. Specifically, `JpsiXCandidateProducer::constrainJpsi4Momentum` SHALL NOT be invoked under preset B; the mother mass m(J/ψ + X) used for the mother mass-window cut and for the stored mother candidate SHALL be computed as the raw sum of the J/ψ candidate's track-based p4() and the bachelor / sub-resonance p4().

Under preset C, the constraint IS applied, both as the dimuon-p4 input to the multi-track Kalman fit (`TwoTrackMassKinematicConstraint`) and as the fallback dimuon p4 when the multi-track fit fails. The mass window remains identical between presets B and C.

#### Scenario: Preset B uses raw track-sum dimuon p4 in B+ mass computation
- **WHEN** preset B is active and a B+ candidate is built from a J/ψ candidate (dimuon mass 3.07 GeV from track fit) and a bachelor K+ track
- **THEN** the producer computes m(B+) = (lvJpsi + lvK).M() where lvJpsi is the unmodified `jpsi.p4()` from `ALCARECOTkAlJpsiXJpsiCandidates`. The `n_jpsi_constraint_attempted_` counter logged at job end SHALL be 0 under preset B

#### Scenario: Preset C invokes the constraint
- **WHEN** preset C is active
- **THEN** `n_jpsi_constraint_attempted_` SHALL equal the number of J/ψ candidates passed to `JpsiXCandidateProducer.produceVccMode` and `produceTrackMode`, and the dimuon p4 used downstream SHALL be the constrained 4-momentum (when the constraint converges) or the raw track sum (fallback, counted in `n_jpsi_constraint_fallback_`)

---

### Requirement: V0 Mass Treatment — Window-Only End-to-End

The V0-mode channels (B0→Ks, Λb, ψ(2S)) SHALL use mass-window selection on the V0 sub-resonance (Ks, Λ), NOT a `MassKinematicConstraint`. This requirement holds end-to-end and under every preset:

- The shared V0Producer clone `ALCARECOTkAlV0Candidates` (`ALCARECOTkAlV0Candidates_cff.py`) SHALL invoke `V0Producer` with its standard `kShortMassCut` / `lambdaMassCut` window cuts on each V0 candidate's reconstructed mass; no `MassKinematicConstraint` or `TwoTrackMassKinematicConstraint` SHALL be inserted into the V0 fit. The clone differs from the central `generalV0Candidates` only in `tkPtCut=0.1`; all mass-handling parameters SHALL be inherited unchanged from the standard module.
- The three V0-mode `JpsiXCandidateProducer` instances SHALL combine each J/ψ candidate with each V0 candidate by summing their reco::Candidate p4() values: `lvM = lvJpsi + xCand.p4()`. They SHALL NOT apply any kinematic constraint to the V0 sub-resonance's 4-momentum.
- The downstream mother mass window (`minMotherMass` / `maxMotherMass` on each producer) is the mother-level selection; it SHALL NOT use a kinematic-constrained dimuon p4 under preset B (per the J/ψ requirement above) and SHALL NOT constrain the V0 p4 under any preset.

#### Scenario: Ks daughter pion masses are not constrained by the V0Producer clone
- **WHEN** a Ks candidate is built from two oppositely-charged tracks
- **THEN** the V0Producer clone applies the pion mass hypothesis to each daughter but does NOT pin the reconstructed Ks mass to `kShortMass`; the reconstructed Ks mass is kept as a free output observable, only required to fall within `kShortMassCut`

#### Scenario: B0 → J/ψ Ks mother mass computation uses raw V0 p4
- **WHEN** preset B is active and a B0 → J/ψ Ks candidate is built
- **THEN** m(B0) = (lvJpsi + lvKs).M() where lvJpsi is the raw track-sum dimuon p4 (per the J/ψ requirement) and lvKs is the raw V0 p4 from `ALCARECOTkAlV0Candidates`. Neither input is mass-constrained

---

### Requirement: Helix Propagation of Daughter Tracks to a Common Candidate Point

Before computing the mother 4-vector sum and the mother invariant mass, `JpsiXCandidateProducer` SHALL propagate each leaf track (the two J/ψ muons; the bachelor track in track mode; the two sub-resonance daughters in vcc mode) from that track's native point of closest approach to the beam line **to a common candidate-level point** along the track's helix, and rebuild that daughter's 4-momentum from the propagated `(p_x, p_y, p_z)` together with the daughter's species mass hypothesis.

The propagation SHALL use the following existing CMSSW classes — no custom helix-propagation code SHALL be introduced:

- **`ClosestApproachInRPhi`** (`TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h`) for every 2-body sub-resonance pair: the two J/ψ muons (all channels), and the two sub-resonance daughters in vcc mode (K*0 → Kπ, φ → KK, Ks → ππ, Λ → pπ). The producer SHALL build two `FreeTrajectoryState`s from the daughter `TrackRef`s, construct a `ClosestApproachInRPhi`, call `calculate(...)`, and on `status() == true` read both propagated daughter momenta from `trajectoryParameters()` and the pair's vertex estimate from `crossingPoint()`.
- **`AnalyticalImpactPointExtrapolator`** (`TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h`) for the bachelor track in track mode (B+ K, Bc π). The producer SHALL build a `FreeTrajectoryState` from the bachelor `TrackRef`, construct the extrapolator with the `MagneticField` available via the `TransientTrackBuilder::field()` ESData, call `extrapolate(fts, target_point)` with `target_point` set to the dimuon `crossingPoint()` obtained from the J/ψ `ClosestApproachInRPhi`, and on `result.isValid() == true` read the propagated momentum from `result.globalMomentum()`.

The propagation SHALL NOT run a Kalman vertex fit, SHALL NOT iterate, and SHALL NOT update the track or vertex covariance.

This requirement holds under **every preset (A, B, and C)**. Under preset C the multi-track Kalman fit subsequently overwrites the propagated mother 4-vector with its own constrained-fit result; the propagation is then used only for the initial mass-window cut and for the fallback dimuon p4 when the Kalman fit fails. Under presets A and B the propagated 4-vector IS the stored mother 4-vector.

This requirement is a kinematic refinement, not a selection. It SHALL NOT change the candidate count emitted by any `JpsiXCandidateProducer` instance on a fixed input dataset: the mass-window cut tolerates the ~1–3 MeV shift, the kinematic and DOCA cuts use raw track-level helix parameters unchanged by the propagation, and the geometric `vtxM` is unchanged by the propagation step (it is the *input* to the propagator).

#### Scenario: B+ → J/ψ K+ mother 4-vector built from ClosestApproachInRPhi + AnalyticalImpactPointExtrapolator
- **WHEN** a B+ candidate is built in `produceTrackMode` from a J/ψ candidate and a bachelor K+ track
- **THEN** the producer constructs `lvJpsi_prop` as the sum of the two muon 4-vectors evaluated at the propagated `(p_x, p_y, p_z)` returned by `ClosestApproachInRPhi::trajectoryParameters()` applied to the two muon helices, propagates the bachelor K+ to the dimuon `crossingPoint()` via `AnalyticalImpactPointExtrapolator::extrapolate(...)` to obtain `lvBach_prop` from the propagated state's `globalMomentum()`, and stores `lvM = lvJpsi_prop + lvBach_prop` as the candidate's 4-vector (subject to overwrite by the preset-C Kalman fit). The mass-window cut at `[5.0, 5.5] GeV` is applied to `lvM.M()` so computed

#### Scenario: V0-mode sub-resonance 4-vector also built from ClosestApproachInRPhi
- **WHEN** a Bs → J/ψ φ candidate (or any vcc-mode candidate) is built
- **THEN** the producer applies `ClosestApproachInRPhi` once on the dimuon helices and once on the sub-resonance daughter helices (φ → K⁺K⁻ in this case), reads the propagated 4-vectors for all four daughters from each call's `trajectoryParameters()`, and stores `lvM` as the sum of the four propagated 4-vectors. `AnalyticalImpactPointExtrapolator` is NOT invoked in vcc mode

#### Scenario: Propagation does not change the candidate count
- **WHEN** the producer is run twice on the same input — once with the propagation enabled, once with it disabled (call sites reverted to the pre-change raw track-momentum construction)
- **THEN** the number of candidates emitted in each of the 7 output collections is bit-identical between the two runs. Only the per-candidate `lvM` value (and therefore the candidate's reported `mass()`, `pt()`, `eta()`, `phi()`) differs

#### Scenario: ESConsumes the TransientTrackBuilder under every preset
- **WHEN** the producer is constructed under any preset
- **THEN** it consumes the `TransientTrackBuilder` EventSetup record unconditionally (the propagator requires it). The previous gating on preset-C / J/ψ-constraint paths is removed

---

### Requirement: Per-Particle Track Reference Layout (Stage-2 Hand-Off Contract)

Every leaf daughter of every saved mother candidate (post-remap, in the `ALCARECOTkAlJpsiX{BPlus,B0Kstar,B0Ks,BsPhi,Lambdab,Psi2S,Bc}Resonances` collections) SHALL be a `reco::RecoChargedCandidate` carrying BOTH:

- A species identifier sufficient for Stage-2 to select the correct mass hypothesis. **Two regimes coexist by construction in CMSSW**, both acceptable, both verified by the `_diag_daughter_layout.py` diagnostic on the 100k Condor output:
  - **PDG regime** — `|daughter->pdgId()|` carries the expected species number:
    - Muon daughters (the two leaves of every J/ψ `daughter(0)`): `|pdgId()| == 13`.
    - B+ bachelor (track mode): `|pdgId()| == 321` (kaon).
    - Bc bachelor (track mode): `|pdgId()| == 211` (pion).
    - B0 → J/ψ K*0 sub-resonance daughters: PDG-id pair `{321, 211}` (Kπ).
    - Bs → J/ψ φ sub-resonance daughters: both `|pdgId()| == 321` (KK).
    - These four channels are built by `JpsiXCandidateProducer.cc` (track mode) and `TwoBodyDecayCandidateProducer` (K*0, φ), which set the daughter PDG explicitly.
  - **Mass regime** — `daughter->mass()` carries the species mass within 1 MeV of the PDG value:
    - B0 → J/ψ Ks sub-resonance daughters: both daughters' `mass()` ≈ 0.139570 GeV (pion).
    - ψ(2S) → J/ψ Ks sub-resonance daughters: both daughters' `mass()` ≈ 0.139570 GeV (pion).
    - Λb → J/ψ Λ sub-resonance daughters: one daughter's `mass()` ≈ 0.938272 GeV (proton), the other's ≈ 0.139570 GeV (pion).
    - These three V0-mode channels are built by the central `V0Producer`, which constructs `RecoChargedCandidate` daughters with `pdgId() == 0` but with the correct species `mass()` set via `setMass()`.

  Stage-2 SHALL prefer the PDG regime when `|daughter->pdgId()| != 0` and fall back to the mass regime when `|daughter->pdgId()| == 0`. Both regimes resolve to the same set of species masses; the choice between them is dictated by which upstream producer built the candidate, not by the Stage-2 consumer.

- A resolvable `reco::TrackRef` to a hit-bearing track in the saved `ALCARECOTkAlJpsiX` track collection (i.e. `daughter->bestTrack() != nullptr` AND the dereferenced track's `extra()` is non-null AND `recHitsSize() > 0`).

The mother candidate's first daughter (`daughter(0)`) SHALL always be a `reco::VertexCompositeCandidate` representing the J/ψ with two `RecoChargedCandidate` muon daughters in the **PDG regime** (J/ψ is built by `TwoBodyDecayCandidateProducer`). The second daughter (`daughter(1)`) SHALL be either a `reco::RecoChargedCandidate` (track-mode: B+, Bc — PDG regime) or a `reco::VertexCompositeCandidate` (vcc-mode: K*0, φ — PDG regime; Ks, Λ — mass regime) whose own daughters obey the contract.

This layout is the hand-off contract to Stage-2 CVH. The downstream splitter (e.g. `JpsiKCandidateSplitter` from `add-jpsi-x-stage2-bplus-cvh`) and any future N-track CVH refit MUST be able to read each daughter's species (via PDG or mass, per the regime in effect) to select the correct mass hypothesis (muon, kaon, pion, proton) and dereference the `TrackRef` to feed the leaf track's rec hits into `Geant4ePropagator`. Any change to producer code that breaks this layout (replacing a `RecoChargedCandidate` with a generic `LeafCandidate`, dropping the explicit `setTrack(...)` call, dropping BOTH the PDG tag and the mass tag, or returning a `MasterClone`-only reference without a resolvable `TrackRef`) constitutes a regression.

Smoke verification on the 50k preset-B Condor output (file
`jpsix_alcareco_presetB_0029A508-5D93-E611-9D38-02163E013482.root`,
4457 events, 17760 mother candidates inspected; output preserved in
`condor/jpsix_alcareco/_layout_check_smoke.txt`):

| Channel | Regime | n_candidates | issues | species histogram |
|---|---|---:|---:|---|
| BPlus | pdg | 4947 | 0 | {321: 4947} |
| Bc | pdg | 8161 | 0 | {211: 8161} |
| B0Kstar | pdg | 4013 | 0 | {321: 4013, 211: 4013} |
| BsPhi | pdg | 377 | 0 | {321: 754} |
| B0Ks | mass | 135 | 0 | {211: 270} |
| Psi2S | mass | 101 | 0 | {211: 202} |
| Lambdab | mass | 26 | 0 | {2212: 26, 211: 26} |

#### Scenario: B+ → J/ψ K+ saved candidate carries per-particle (PDG, track) pairs in the PDG regime
- **WHEN** a saved B+ candidate in `ALCARECOTkAlJpsiXBPlusResonances` is inspected
- **THEN** `numberOfDaughters() == 2`; `daughter(0)` is a `VertexCompositeCandidate` whose two leaf daughters are `RecoChargedCandidate`s with `|pdgId()| == 13` and resolvable `bestTrack()`s; `daughter(1)` is a `RecoChargedCandidate` with `|pdgId()| == 321` and a resolvable `bestTrack()`. Each of those three tracks resolves into the saved `ALCARECOTkAlJpsiX` collection with `recHitsSize() > 0`

#### Scenario: B0 → J/ψ K*0 saved candidate carries per-particle (PDG, track) pairs in the PDG regime
- **WHEN** a saved B0→K*0 candidate in `ALCARECOTkAlJpsiXB0KstarResonances` is inspected
- **THEN** `daughter(1)` is a `VertexCompositeCandidate` with two `RecoChargedCandidate` daughters whose PDG-id pair is exactly `{|321|, |211|}` (one kaon + one pion), each with a resolvable `bestTrack()` and rec hits

#### Scenario: Λb → J/ψ Λ saved candidate carries per-particle (mass, track) pairs in the mass regime
- **WHEN** a saved Λb candidate in `ALCARECOTkAlJpsiXLambdabResonances` is inspected
- **THEN** `daughter(1)` is a `VertexCompositeCandidate` from `V0Producer` with two `RecoChargedCandidate` daughters whose `pdgId() == 0` (V0Producer does not set per-daughter PDG) but whose `mass()` values are exactly `{0.938272 ± 0.001, 0.139570 ± 0.001} GeV` (one proton + one pion), each with a resolvable `bestTrack()` and rec hits

#### Scenario: Stage-2 splitter can pair (track, mass) unambiguously across both regimes
- **WHEN** a downstream splitter walks the leaf daughters of any saved mother candidate
- **THEN** for each leaf it can read either `pdgId()` (PDG regime; non-V0 channels) or `mass()` (mass regime; V0 channels) to obtain the species, AND read `bestTrack()` to obtain the rec hits, with no further bookkeeping or external mapping required. The splitter SHALL apply the rule "use `pdgId` when non-zero, otherwise use `mass`" uniformly across all 7 channels
