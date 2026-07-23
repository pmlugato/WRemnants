## ADDED Requirements

### Requirement: TkAlJpsiMuMu Persists dE/dx Value Maps

`ALCARECOTkAlJpsiMuMu` SHALL persist three dE/dx value maps
projected onto its cloned track collection: Harmonic2 (strip),
PixelHarmonic2, and Joint (via `alcaDedxJointEstimator`), matching
the persist policy already followed by `ALCARECOTkAlJpsiX` and
`ALCARECOTkAlDstToD0Pi`. The maps SHALL be produced by
`DeDxValueMapProjector` clones keyed by the `originalIndex` map
that `AlignmentTrackSelectorWithIndexMap` already emits. Production-
permanent.

#### Scenario: dE/dx maps present in production output

- **WHEN** the TkAlJpsiMuMu recoskim runs
- **THEN** `edmDumpEventContent` on the output lists
  `ALCARECOTkAlJpsiMuMuDeDxHarmonic2`,
  `ALCARECOTkAlJpsiMuMuDeDxPixelHarmonic2`, and
  `ALCARECOTkAlJpsiMuMuDeDxAllHarmonic2`, each non-empty on
  events with at least one persisted track

### Requirement: TkAlJpsiMuMu Persists Its Tight Muon Selector Output

`ALCARECOTkAlJpsiMuMu` SHALL persist the `reco::MuonCollection`
produced by its tight muon selector (`TkAlGoodIdMuonSelector`,
already present in the cff as `ALCARECOTkAlJpsiMuMuGoodMuons`) with
a `keep *_ALCARECOTkAlJpsiMuMuGoodMuons_*_*` line in its output cff.
This is the tight-selector reference used by the loose-vs-tight
muon comparison plot. Production-permanent.

#### Scenario: Tight muons present in production output

- **WHEN** the TkAlJpsiMuMu recoskim runs
- **THEN** `edmDumpEventContent` lists
  `ALCARECOTkAlJpsiMuMuGoodMuons`, non-empty on every accepted
  event (accepted events by definition have >= 2 tight muons)

### Requirement: TkAlJpsiMuMu Track -> Muon Association

`ALCARECOTkAlJpsiMuMu` SHALL emit and persist a track-to-muon
`edm::Association` mapping the cloned track collection to the
persisted tight muon collection, using the same producer
(`AlignmentTrackToMuonAssociator`) that the J/psi+X stream already
uses. Production-permanent.

#### Scenario: Association resolves persisted tracks to persisted muons

- **WHEN** an event has 2 tight muons and the cloned track
  collection contains the 2 corresponding inner tracks
- **THEN** the association resolves each track to the correct
  muon by `innerTrack().key()`, and unmapped tracks (e.g. tracks
  from a J/psi candidate whose muons dropped out of the tight
  collection due to the muon selector's `filter = False`
  emission constraint) resolve to a null muon key
