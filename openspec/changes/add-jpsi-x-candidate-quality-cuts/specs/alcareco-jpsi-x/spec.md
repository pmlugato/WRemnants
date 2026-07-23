## MODIFIED Requirements

### Requirement: Stage-1 Shared Candidate Construction
The AlCaReco sequence SHALL build all Stage-1 candidates before any per-channel B-meson construction. No HLT filter SHALL be applied in the default sequence; an `ALCARECOTkAlJpsiXHLT` module SHALL be defined in the cff but excluded from `seqALCARECOTkAlJpsiX`. J/ψ → μ+μ- candidates SHALL be built by `TwoBodyDecayCandidateProducer` with muon filter, mass window **2.95–3.25 GeV** (±5σ at 30 MeV CMS dimuon resolution; replaces 2.7–3.4 GeV), and daughter mass 0.105 GeV (muon). K*0(892) → K±π∓ candidates SHALL be built by `TwoBodyDecayCandidateProducer` with `firstDaughterMass = 0.494 GeV` (K), `secondDaughterMass = 0.140 GeV` (π), `tryBothChargeAssignments = True`, mass window **0.80–0.99 GeV** (±95 MeV ≈ ±2.0 Γ around PDG K*0 mass 895.55 MeV; replaces 0.75–1.05 GeV), and `maxDaughterEta = 2.4` (Phase 2); all passing candidates from both charge assignments SHALL be emitted. φ(1020) → K+K- candidates SHALL be built by `TwoBodyDecayCandidateProducer` (symmetric, daughter mass 0.494 GeV) with mass window **0.990–1.040 GeV** (lower edge at K+K- threshold 0.987 GeV; upper edge tightened from 1.060 GeV; replaces 0.99–1.06 GeV) and `maxDaughterEta = 2.4` (Phase 2). Ks and Λ candidates SHALL come from `ALCARECOTkAlV0Candidates` (the local V0Producer clone with `tkPtCut = 0.1`), shared with existing `TkAlKsToPiPi` and `TkAlLambdaToProtonPi` streams via CMS framework deduplication. No separate bachelor track pre-filter module SHALL be created; pT thresholds for bachelor tracks SHALL be applied internally by `JpsiXCandidateProducer` (see Stage-2 requirement). Any `AlignmentTrackSelectorModule`-based pre-filter rewrites TrackExtraRefs in the cloned output, which would break the `CompositeDaughterTrackProducer` → `VertexCompositeCandidateRemapper` lookup chain.

#### Scenario: J/ψ candidate built from good muon tracks
- **WHEN** an event contains two generalTracks that are the innerTrack of a muon in `ALCARECOTkAlJpsiXGoodMuons` and have invariant mass in [2.95, 3.25] GeV under muon mass hypothesis
- **THEN** `ALCARECOTkAlJpsiXJpsiCandidates` contains at least one nested `VertexCompositeCandidate` with two `RecoChargedCandidate` daughters whose `TrackRef`s point into `generalTracks`

#### Scenario: J/ψ candidate outside tightened window rejected
- **WHEN** two muon tracks have invariant mass in [2.7, 2.95) GeV or (3.25, 3.4] GeV
- **THEN** the combination does NOT appear in `ALCARECOTkAlJpsiXJpsiCandidates`

#### Scenario: K*0 both-assignment combinatorics
- **WHEN** a track pair has charges (+,-) and the K+π- assignment gives invariant mass in [0.80, 0.99] GeV but the π+K- assignment does not
- **THEN** `ALCARECOTkAlJpsiXKstarCandidates` contains exactly one candidate for that pair

#### Scenario: K*0 both assignments accepted
- **WHEN** a track pair has charges (+,-) and both K+π- and π+K- assignments give invariant mass in [0.80, 0.99] GeV
- **THEN** `ALCARECOTkAlJpsiXKstarCandidates` contains two candidates for that pair (one per assignment)

#### Scenario: φ above-threshold pair accepted at lower edge
- **WHEN** a K+K- pair has invariant mass in [0.990, 1.040] GeV
- **THEN** the pair appears in `ALCARECOTkAlJpsiXPhiCandidates` (lower edge at K+K- threshold; no change from original)

#### Scenario: φ over-threshold pair rejected at upper edge
- **WHEN** a K+K- pair has invariant mass in (1.040, 1.060] GeV (above the tightened upper edge)
- **THEN** the pair does NOT appear in `ALCARECOTkAlJpsiXPhiCandidates`

#### Scenario: K*0 forward-eta daughter rejected (Phase 2)
- **WHEN** either daughter track of a K*0 candidate has |η| > 2.4
- **THEN** the candidate does NOT appear in `ALCARECOTkAlJpsiXKstarCandidates`

---

### Requirement: Stage-2 Per-Channel B-meson and Quarkonium Candidate Construction
`JpsiXCandidateProducer` SHALL operate in two modes. In **track mode** (B+, Bc), it iterates over `generalTracks` directly, skipping tracks with pT below a configurable `minBachelorPt` (**0.5 GeV for B+, unchanged; 0.3 GeV for Bc**, lowered from 0.5 GeV), and combines each passing track with each J/ψ candidate under a configurable mass hypothesis and mother mass window. The stored bachelor `TrackRef` points into `generalTracks` (NOT an intermediate clone), preserving the TrackExtraRef chain required by downstream modules. In **VCC mode** (all other channels), it combines each J/ψ candidate with each candidate from a second VCC (X) and applies a mother mass window. In both modes, the following optional parameters (all with backward-compatible defaults) SHALL be supported:

- `minJpsiPt` (default 0 GeV): cut on J/ψ candidate pT before the combination loop; configured **3.0 GeV** for all 7 channels (Phase 2; largely redundant with minMotherPt > 5 GeV; provides early loop exit)
- `minMotherPt` (default 0 GeV): cut on reconstructed mother pT after 4-momentum sum; configured **5.0 GeV** for B+, B0→K*0, Bs, Bc (Phase 2)
- `maxBachelorEta` (default +∞, track mode only): |η| cut on bachelor track; configured **2.4** for B+ and Bc (Phase 2)
- `maxJpsiAlphaBS` (default +∞): pointing angle of J/ψ to beamspot using VCC vertex midpoint; configured **1.0** for all channels (Phase 2); requires `offlineBeamSpot` InputTag consumed when parameter is present
- `maxBachelorIPToJpsiVertex` (default +∞, track mode only): 3D DCA of bachelor helix to J/ψ vertex midpoint; configured **10.0 mm** for B+ and Bc as initial value to be tuned from data (Phase 2)

In both modes the output SHALL be a `VertexCompositeCandidateCollection` where each candidate has daughter(0) = J/ψ (nested VCC with original daughter structure preserved) and daughter(1) = X. The J/ψ four-momentum used in the mother mass calculation SHALL be the stored candidate p4() (raw track-based, consistent with `TwoBodyDecayCandidateProducer` convention). The seven required channel instances and their parameters are defined in the proposal channel catalogue.

Note: `TwoBodyDecayCandidateProducer` does not perform a kinematic vertex fit and does not store vertex quality in the VCC; J/ψ vertex probability cuts are deferred to Phase 3.

#### Scenario: B+ candidate within mass window
- **WHEN** a J/ψ candidate with pT ≥ 3.0 GeV and a bachelor track with pT ≥ 0.5 GeV and |η| < 2.4 and kaon mass hypothesis (0.494 GeV) combine to give invariant mass in [5.0, 5.5] GeV and the mother pT ≥ 5.0 GeV
- **THEN** `ALCARECOTkAlJpsiXBPlusCandidates` contains a nested candidate with PDG ID 521, daughter(0) = J/ψ sub-candidate, daughter(1) = bachelor `RecoChargedCandidate`

#### Scenario: Soft Bc bachelor track rejected
- **WHEN** a bachelor track has pT < 0.3 GeV (Bc channel)
- **THEN** the track is not combined with any J/ψ candidate in `ALCARECOTkAlJpsiXBcCandidates`

#### Scenario: Very low-pT J/ψ candidate skipped (Phase 2)
- **WHEN** a J/ψ candidate has pT < 3.0 GeV
- **THEN** the J/ψ is not used in any B-meson combination loop across all channels

#### Scenario: Low-pT mother candidate rejected
- **WHEN** a J/ψ + bachelor track combination falls in the B+ mass window but the reconstructed mother pT < 5.0 GeV
- **THEN** the combination does NOT appear in `ALCARECOTkAlJpsiXBPlusCandidates`

#### Scenario: B+ candidate outside mass window rejected
- **WHEN** a J/ψ + kaon track combination has invariant mass outside [5.0, 5.5] GeV
- **THEN** the combination does NOT appear in `ALCARECOTkAlJpsiXBPlusCandidates`

#### Scenario: J/ψ pointing angle cut applied (Phase 2)
- **WHEN** a J/ψ candidate has alphaBS > 1.0 (angle between momentum vector and beamspot flight direction)
- **THEN** the J/ψ is rejected before any B-meson combination in all channels

#### Scenario: Bs → J/ψ φ candidate built from φ sub-candidate
- **WHEN** a J/ψ candidate (pT ≥ 5.0 GeV) and a φ candidate (K+K- mass in [1.000, 1.040] GeV) combine within [5.2, 5.6] GeV and mother pT ≥ 5.0 GeV
- **THEN** `ALCARECOTkAlJpsiXBsPhiCandidates` contains a 4-body nested candidate: daughter(0)=J/ψ(μ+,μ-), daughter(1)=φ(K+,K-)

#### Scenario: ψ(2S) and B0 channels share Ks input without conflict
- **WHEN** a Ks candidate from `ALCARECOTkAlV0Candidates` is used by both `ALCARECOTkAlJpsiXPsi2SCandidates` and `ALCARECOTkAlJpsiXB0KsCandidates`
- **THEN** both VCC collections are independently produced without double-counting and can independently pass or fail their respective mass windows

#### Scenario: minMotherPt default preserves V0-mode channel behaviour
- **WHEN** `minMotherPt` is not set (default 0 GeV) for B0→Ks, Λb, ψ(2S) channels
- **THEN** all candidates passing the mass window are stored, regardless of mother pT
