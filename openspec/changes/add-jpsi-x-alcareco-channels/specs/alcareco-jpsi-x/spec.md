## ADDED Requirements

### Requirement: Stage-1 Shared Candidate Construction
The AlCaReco sequence SHALL build all Stage-1 candidates before any per-channel B-meson construction. No HLT filter SHALL be applied in the default sequence; an `ALCARECOTkAlJpsiXHLT` module SHALL be defined in the cff but excluded from `seqALCARECOTkAlJpsiX`. J/ψ → μ+μ- candidates SHALL be built by `TwoBodyDecayCandidateProducer` with muon filter, mass window 2.7–3.4 GeV, and daughter mass 0.105 GeV (muon). K*0(892) → K±π∓ candidates SHALL be built by `TwoBodyDecayCandidateProducer` with `firstDaughterMass = 0.494 GeV` (K), `secondDaughterMass = 0.140 GeV` (π), `tryBothChargeAssignments = True`, and mass window 0.75–1.05 GeV; all passing candidates from both charge assignments SHALL be emitted. φ(1020) → K+K- candidates SHALL be built by `TwoBodyDecayCandidateProducer` (symmetric, daughter mass 0.494 GeV) with mass window 0.99–1.06 GeV. Ks and Λ candidates SHALL come from `ALCARECOTkAlV0Candidates` (the local V0Producer clone with `tkPtCut = 0.1`), shared with existing `TkAlKsToPiPi` and `TkAlLambdaToProtonPi` streams via CMS framework deduplication. No separate bachelor track pre-filter module SHALL be created; the pT > 0.5 GeV cut for bachelor tracks in the B+ and Bc channels SHALL be applied internally by `JpsiXCandidateProducer` (see Stage-2 requirement). Any `AlignmentTrackSelectorModule`-based pre-filter rewrites TrackExtraRefs in the cloned output, which would break the `CompositeDaughterTrackProducer` → `VertexCompositeCandidateRemapper` lookup chain.

#### Scenario: J/ψ candidate built from good muon tracks
- **WHEN** an event contains two generalTracks that are the innerTrack of a muon in `ALCARECOTkAlJpsiXGoodMuons` and have invariant mass in [2.7, 3.4] GeV under muon mass hypothesis
- **THEN** `ALCARECOTkAlJpsiXJpsiCandidates` contains at least one nested `VertexCompositeCandidate` with two `RecoChargedCandidate` daughters whose `TrackRef`s point into `generalTracks`

#### Scenario: K*0 both-assignment combinatorics
- **WHEN** a track pair has charges (+,-) and the K+π- assignment gives invariant mass in [0.75, 1.05] GeV but the π+K- assignment does not
- **THEN** `ALCARECOTkAlJpsiXKstarCandidates` contains exactly one candidate for that pair

#### Scenario: K*0 both assignments accepted
- **WHEN** a track pair has charges (+,-) and both K+π- and π+K- assignments give invariant mass in [0.75, 1.05] GeV
- **THEN** `ALCARECOTkAlJpsiXKstarCandidates` contains two candidates for that pair (one per assignment)

---

### Requirement: Stage-2 Per-Channel B-meson and Quarkonium Candidate Construction
`JpsiXCandidateProducer` SHALL operate in two modes. In **track mode** (B+, Bc), it iterates over `generalTracks` directly, skipping tracks with pT below a configurable `minBachelorPt` (0.5 GeV), and combines each passing track with each J/ψ candidate under a configurable mass hypothesis and mother mass window. The stored bachelor `TrackRef` points into `generalTracks` (NOT an intermediate clone), preserving the TrackExtraRef chain required by downstream modules. In **VCC mode** (all other channels), it combines each J/ψ candidate with each candidate from a second VCC (X) and applies a mother mass window. In both modes the output SHALL be a `VertexCompositeCandidateCollection` where each candidate has daughter(0) = J/ψ (nested VCC with original daughter structure preserved) and daughter(1) = X. The J/ψ four-momentum used in the mother mass calculation SHALL be the stored candidate p4() (raw track-based, consistent with `TwoBodyDecayCandidateProducer` convention). The mass windows in the channel catalogue are wide enough to accommodate J/ψ mass resolution spread. The seven required channel instances and their parameters are defined in the proposal channel catalogue.

#### Scenario: B+ candidate within mass window
- **WHEN** a J/ψ candidate and a bachelor track with kaon mass hypothesis (0.494 GeV) combine to give invariant mass in [5.0, 5.5] GeV
- **THEN** `ALCARECOTkAlJpsiXBPlusCandidates` contains a nested candidate with PDG ID 521, daughter(0) = J/ψ sub-candidate, daughter(1) = bachelor `RecoChargedCandidate`

#### Scenario: B+ candidate outside mass window rejected
- **WHEN** a J/ψ + kaon track combination has invariant mass outside [5.0, 5.5] GeV
- **THEN** the combination does NOT appear in `ALCARECOTkAlJpsiXBPlusCandidates`

#### Scenario: Bs → J/ψ φ candidate built from φ sub-candidate
- **WHEN** a J/ψ candidate and a φ candidate combine within [5.2, 5.6] GeV
- **THEN** `ALCARECOTkAlJpsiXBsPhiCandidates` contains a 4-body nested candidate: daughter(0)=J/ψ(μ+,μ-), daughter(1)=φ(K+,K-)

#### Scenario: ψ(2S) and B0 channels share Ks input without conflict
- **WHEN** a Ks candidate from `ALCARECOTkAlV0Candidates` is used by both `ALCARECOTkAlJpsiXPsi2SCandidates` and `ALCARECOTkAlJpsiXB0KsCandidates`
- **THEN** both VCC collections are independently produced without double-counting and can independently pass or fail their respective mass windows

---

### Requirement: Recursive Track Extraction and Deduplication
`CompositeDaughterTrackProducer` SHALL accept a configurable list of input `VertexCompositeCandidateCollection` labels (`srcs`) and produce one `TrackCollection` containing the unique set of leaf `RecoChargedCandidate` tracks reachable by recursing through all candidate trees in all input collections. Deduplication SHALL be by `TrackExtraRef.key()` (which maps uniquely to the original generalTracks index). Each output Track SHALL be a shallow copy preserving its `TrackExtraRef` so that downstream `AlignmentTrackSelectorWithIndexMapModule` and `VertexCompositeCandidateRemapper` can correctly navigate the chain `selected track → intermediate track → generalTracks key → candidate daughter`.

#### Scenario: B+ candidate extracts 3 unique tracks
- **WHEN** `CompositeDaughterTrackProducer` processes a B+ candidate tree with daughter(0)=J/ψ(μ+,μ-) and daughter(1)=K+
- **THEN** 3 tracks appear in the output collection for that candidate (μ+, μ-, K+)

#### Scenario: B0 → J/ψ K*0 extracts 4 unique tracks
- **WHEN** a B0 candidate tree has daughter(0)=J/ψ(μ+,μ-) and daughter(1)=K*0(K+,π-)
- **THEN** 4 tracks appear in the output for that candidate (μ+, μ-, K+, π-)

#### Scenario: Shared muon track deduplicated across two channels
- **WHEN** a μ+ track is a leaf in both a B+ candidate and a ψ(2S) candidate in the same event
- **THEN** that μ+ track appears exactly once in the merged output `TrackCollection`

---

### Requirement: Single Combined AlignmentTrackSelector and dE/dx
One `AlignmentTrackSelectorWithIndexMapModule` instance (`ALCARECOTkAlJpsiX`) SHALL run on the merged track collection from `CompositeDaughterTrackProducer`, applying pT > 0.4 GeV, |η| < 3.5, nHit ≥ 0, `GlobalSelector.applyGlobalMuonFilter = False`, `filter = True`. Three `DeDxValueMapProjector` instances SHALL project dE/dx values (strip Harmonic2, pixel Harmonic2, joint alcaDedxJointEstimator) onto the cloned track collection, following the pattern of `TkAlKsToPiPi`. Because `DeDxValueMapProjector` operates on the single merged track collection, dE/dx values are produced for every leaf track regardless of which channel contributed it; physically meaningful PID information is expected only for kaon-bearing channels (B+, B0→K*0, Bs) and the proton-bearing channel (Λb). The `originalIndex` ValueMap emitted by `AlignmentTrackSelectorWithIndexMapModule` SHALL be consumed by all seven `VertexCompositeCandidateRemapper` instances.

#### Scenario: Empty event rejected
- **WHEN** no B-meson candidates from any of the 7 channels survive their respective mass windows
- **THEN** `CompositeDaughterTrackProducer` produces an empty collection, `AlignmentTrackSelectorWithIndexMapModule` returns `false` (filter=True), and the event does not appear in the `TkAlJpsiX` ALCARECO output

#### Scenario: Single-channel event saved with all leaf tracks
- **WHEN** exactly one B+ → J/ψ K+ candidate is found and passes the mass window
- **THEN** `ALCARECOTkAlJpsiX` contains exactly 3 cloned tracks with full hit and cluster information

---

### Requirement: Per-Channel Candidate Remapping
Seven `VertexCompositeCandidateRemapper` instances SHALL re-key all leaf `RecoChargedCandidate` daughter `TrackRef`s in each channel's B-meson candidate collection onto the cloned `ALCARECOTkAlJpsiX` track collection. Each remapper SHALL use `selectedTracks = ALCARECOTkAlJpsiX`, `intermediateTracks = ALCARECOTkAlJpsiXAllTracks`, and `originalIndexMap = (ALCARECOTkAlJpsiX, originalIndex)`. Candidates whose leaf daughters did not all survive `AlignmentTrackSelectorWithIndexMapModule` SHALL be silently dropped.

#### Scenario: Remapped B+ candidate navigates to cloned track
- **WHEN** a remapped B+ candidate is accessed from `ALCARECOTkAlJpsiXBPlusResonances`
- **THEN** `candidate.daughter(1).get<reco::RecoChargedCandidate>().track()` returns a `TrackRef` into `ALCARECOTkAlJpsiX` (not `generalTracks`)

#### Scenario: Remapped ψ(2S) candidate with nested Ks navigates correctly
- **WHEN** a remapped ψ(2S) candidate has daughter(0)=J/ψ(μ+,μ-) and daughter(1)=Ks(π+,π-)
- **THEN** all four leaf `TrackRef`s point into `ALCARECOTkAlJpsiX`

---

### Requirement: ALCARECO Output Content and Stream Registration
The `TkAlJpsiX` stream SHALL write cloned tracks (+ extras + hits + clusters), all three dE/dx ValueMaps, all seven remapped candidate collections, and standard event-level products (trigger results, DCS status, primary vertices), following the `_noDrop` + `drop *` output command convention. The stream SHALL be registered in `AlCaRecoStreams_cff.py`, `autoAlca.py` (under `Charmonium` and `MuOniaParked`), and `AlCaRecoOutput_cff.py`.

#### Scenario: ALCARECO output navigable from candidates to tracks
- **WHEN** a `TkAlJpsiX` ALCARECO event is read back
- **THEN** each remapped B-meson candidate's leaf TrackRefs are valid within the ALCARECO file (no cross-file dangling refs), and the dE/dx ValueMap is accessible for each cloned track

#### Scenario: Stream activated via autoAlca
- **WHEN** a cmsRun job targets a `Charmonium` dataset with `ALCA:TkAlJpsiX`
- **THEN** the framework imports `ALCARECOTkAlJpsiX_cff.py`, runs `pathALCARECOTkAlJpsiX`, and writes output via `ALCARECOStreamTkAlJpsiX`
