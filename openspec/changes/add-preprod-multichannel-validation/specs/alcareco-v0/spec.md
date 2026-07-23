## ADDED Requirements

### Requirement: Standalone V0 AlCaReco Streams Run Under Our Build

`ALCARECOTkAlKShortTracks` and `ALCARECOTkAlLambdaTracks` SHALL run
under `CMSSW_10_6_17_patch1` sharing the local loosened
`ALCARECOTkAlV0Candidates` V0 producer clone (`vtxChi2Cut = 15`,
`vtxDecaySigXYCut = 10`, `cosThetaXYCut = 0.995`, `tkIPSigXYCut = 1`,
`tkPtCut = 0.1`). Each stream SHALL run its species-specific
`VertexCompositeCandidateSelector` (KShort / Lambda) into
`AlignmentTrackSelectorWithIndexMap`, produce dE/dx-projected value
maps against the cloned track collection, and persist the outputs
listed in the matching `_Output_cff.py`.

#### Scenario: Both V0 streams produce non-empty output

- **WHEN** both V0 recoskims run over 1000 events of a Run2016H
  SingleMuon RAW file
- **THEN** the Kshort stream persists its cloned tracks, its VCC,
  and its dE/dx maps; the Lambda stream persists the analogous
  products; and cmsRun exits 0

### Requirement: V0 Streams Reuse the Loosened Clone

The two standalone V0 streams SHALL consume `ALCARECOTkAlV0Candidates`
(the local loosened clone), NOT the shared-CMSSW
`generalV0Candidates`. The same underlying V0 collection is used by
the J/psi+X V0-mode channels.

#### Scenario: Underlying V0 collection shared

- **WHEN** J/psi+X (V0-mode channels) and the two standalone V0
  streams all run in the same job
- **THEN** `ALCARECOTkAlV0Candidates` is built exactly once per
  event and its object count is the same for every consumer

### Requirement: Persist Alignment With Other Streams

The V0 stream output cffs SHALL be audited and, where the physics
allows, aligned to the persist policy: cloned tracks with
`originalIndex`, dE/dx maps (Harmonic2 strip, PixelHarmonic2,
Joint), candidate VCC per species, standard L1 / DCS / vertex
bookkeeping. Missing categories SHALL be added and are production-
permanent.

#### Scenario: Audit + alignment recorded

- **WHEN** the audit task completes
- **THEN** the V0 output cffs list the same categories of persisted
  content (where the physics allows) as J/psi+X and D*, with any
  additions committed to the cff (not gated by any flag)
