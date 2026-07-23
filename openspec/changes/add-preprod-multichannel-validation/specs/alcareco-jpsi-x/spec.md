## ADDED Requirements

### Requirement: Single Canonical Preset B Configuration

The `ALCARECOTkAlJpsiX` cff SHALL carry exactly one set of numeric
cut values (preset B) and SHALL NOT expose a runtime preset switch.
Reads of `TKALJPSIX_SELECTION_PRESET`, the `_NON_V0_PRESETS` dict,
the `_PRESETC_KALMAN` dict, and any `preset == 'C'` late-mutation
stanza SHALL be removed. Numeric literals SHALL appear inline in
each per-channel producer clone.

#### Scenario: cff carries no preset switch

- **WHEN** the cff is loaded with no `TKALJPSIX_SELECTION_PRESET`
  environment variable set
- **THEN** no `_TKALJPSIX_SELECTION_PRESET` name exists at module
  scope and every per-channel producer holds numeric literals that
  reproduce preset B exactly

#### Scenario: env var is ignored if set

- **WHEN** `TKALJPSIX_SELECTION_PRESET=C` is set in the environment
- **THEN** the cff ignores it and the output matches any other run

### Requirement: Candidate-Level Event Gate

The `pathALCARECOTkAlJpsiX` path SHALL NOT gate on the tight-muon
collection. `ALCARECOTkAlJpsiXGoodMuons.filter` SHALL be `False`,
`ALCARECOTkAlJpsiXLooseMuons.filter` SHALL be `False`. The event
gate SHALL be a `CandViewCountFilter` requiring at least one
candidate in the merged 8-channel candidate collection
`ALCARECOTkAlJpsiXAnyCandidate`, placed after all per-channel
candidate producers and before `CompositeDaughterTrackProducer`.

#### Scenario: Tracker-only dimuon event is accepted

- **WHEN** an event has two `reco::Muon`s, both `isTrackerMuon =
  True`, `isGlobalMuon = False`, pairing to a dimuon mass in
  `[2.9, 3.3]` GeV
- **THEN** `ALCARECOTkAlJpsiXGoodMuons` is empty,
  `ALCARECOTkAlJpsiXLooseMuons` contains 2 muons,
  `ALCARECOTkAlJpsiXJpsiOnlyCandidates` contains 1 candidate,
  the gate passes, and the event is written

#### Scenario: Zero-candidate event is rejected

- **WHEN** every per-channel candidate producer emits zero
  candidates
- **THEN** `ALCARECOTkAlJpsiXAnyCandidate` is empty, the count
  filter fails, and the event is not written

### Requirement: Merged 8-Channel Candidate Collection

The cff SHALL construct `ALCARECOTkAlJpsiXAnyCandidate` by merging
the 8 per-channel VCC outputs (B+, B0->K*0, B0->Ks, Bs->phi,
Lambda_b, psi(2S), Bc, JpsiOnly). This collection exists to feed
the candidate-level event gate. It SHALL NOT be persisted in the
production output stream — the per-channel `Resonances` collections
already carry the same information with per-channel provenance.

#### Scenario: Merger produces correct size

- **WHEN** an event has 2 B+ candidates and 3 JpsiOnly candidates
- **THEN** `ALCARECOTkAlJpsiXAnyCandidate` has size 5

### Requirement: Rate Estimate Before Candidate Gate Adoption

The candidate-level gate SHALL be adopted only after a smoke run
measures the accept-rate ratio to the previous tight-gate smoke on
the same seed file. If the ratio is < 2x, the switch SHALL be
adopted. If >= 2x, adoption SHALL be blocked pending user review.

#### Scenario: Ratio under 2x — adopt

- **WHEN** the pre-adoption smoke shows accept rate 1.35x the
  previous smoke
- **THEN** the candidate gate is committed and the proposal
  proceeds to smoke tests

#### Scenario: Ratio at or above 2x — block

- **WHEN** the pre-adoption smoke shows accept rate 3.1x
- **THEN** the switch is not committed and the user is asked to
  review before landing
