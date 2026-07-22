# alcareco-to-nanoaod

## ADDED Requirements

### Requirement: AlCaReco input, NanoAOD output

The stream SHALL read a produced AlCaReco EDM file as input and write
a NanoAOD file containing a single flat `Events` tree, using
`PhysicsTools/NanoAOD` table producers and `NanoAODOutputModule`. The
stream SHALL NOT require MINIAOD input and SHALL NOT depend on any
event content absent from the AlCaReco.

#### Scenario: Run the stream on one AlCaReco file

- **WHEN** the `cmsRun` config runs on a single `ALCARECOTkAlJpsiX` file
- **THEN** it produces a NanoAOD file with an `Events` tree
- **AND** the tree is readable with bare ROOT / `uproot` with no CMSSW
  dependency.

#### Scenario: Cross-release read is validated before batch use

- **WHEN** the processing release differs from the AlCaReco production release
- **THEN** a one-file smoke SHALL confirm successful deserialization of
  `reco::Track`, `reco::Muon`, and `reco::VertexCompositeCandidate`
  before any batch production is launched.

### Requirement: CVH refit runs in application mode, not calibration

The stream SHALL run the persisted candidates' tracks through the CVH
refit with `fillGrads` disabled and `fillRunTree` disabled, so that no
gradient or Hessian sidecar payload is produced. The stream SHALL NOT
perform any gradient aggregation or correction solve.

#### Scenario: No calibration payload is emitted

- **WHEN** the stream runs over an AlCaReco file
- **THEN** no gradient/Hessian sidecar tree is written
- **AND** the job produces only the NanoAOD output.

#### Scenario: Refit configuration for the initial version

- **WHEN** the stream runs before real corrections are available
- **THEN** the refit SHALL run with `useIdealGeometry` disabled and
  with no correction file applied
- **AND** the corrected columns SHALL NOT be consumed by a calibration
  fit while in this configuration.

### Requirement: Single-job refit without split-then-recombine

The refit SHALL run all of a channel's per-subsystem makers within a
single `cmsRun` job, sharing one Geant4 master supplied as an
EventSetup product. The stream SHALL NOT require a separate job per
maker, and SHALL NOT require an offline file-based join to pair the
per-subsystem refit results.

#### Scenario: Both makers coexist in one process

- **WHEN** the job schedules the two-track and single-track makers together
- **THEN** both run to completion in one process without a Geant4
  master conflict
- **AND** their per-candidate results are paired in memory by the
  candidate index.

#### Scenario: Un-refit legs remain visible

- **WHEN** one subsystem of a candidate fails to refit
- **THEN** the candidate SHALL still be written with a status column
  identifying which legs were refit
- **AND** the candidate SHALL NOT be silently dropped.

### Requirement: Generic candidate decomposition without per-channel splitting

The refit SHALL derive each candidate's fit subsystems by traversing
the `VertexCompositeCandidate` daughter structure directly, treating a
composite daughter as a joint-fit subsystem and a leaf charged-candidate
daughter as a single track. The stream SHALL NOT require a per-channel
or per-topology splitter module to adapt candidates for the makers.

#### Scenario: Nested candidate is decomposed generically

- **WHEN** a candidate has a composite daughter and a leaf daughter
- **THEN** the composite daughter's leaf tracks are fit as a joint
  subsystem
- **AND** the leaf daughter is fit as a single track
- **AND** no channel-specific adapter module is involved.

#### Scenario: One code path serves every channel

- **WHEN** the refit runs across all persisted channels
- **THEN** the same decomposition logic handles both
  one-bachelor and two-bachelor topologies
- **AND** channels whose candidate is a single two-track resonance.

### Requirement: Cross-reference index columns

The stream SHALL emit integer index columns linking related objects
across tables, since a flat tree cannot carry EDM references: each
candidate's daughters SHALL index into the `Track` table, and each
track SHALL index into the `Muon` table and the `PV` table where the
corresponding association is persisted.

#### Scenario: Candidate daughter resolves to its track row

- **WHEN** a candidate table row is written
- **THEN** it carries an index column per daughter giving that
  daughter's row in the `Track` table
- **AND** reading that row yields the daughter's dE/dx, hit counts and
  corrected track quantities.

#### Scenario: Missing association is representable

- **WHEN** a track has no associated muon or primary vertex
- **THEN** its index column carries a sentinel value marking the
  absence rather than a valid row index.

### Requirement: Event-level trigger and detector-condition branches

The stream SHALL persist the L1 trigger decision bits from the
`L1GlobalTriggerReadoutRecord` and the `DcsStatus` detector-condition
quantities, including the magnet current, as event-level branches.

#### Scenario: Magnet current available per event

- **WHEN** the `Events` tree is written
- **THEN** it contains an event-level magnet-current branch derived
  from `DcsStatus`
- **AND** L1 decision columns derived from the
  `L1GlobalTriggerReadoutRecord`.

### Requirement: Candidate tables carry raw and corrected quantities

For each persisted candidate collection the stream SHALL emit a flat
table carrying both the raw AlCaReco quantities taken from the
`reco::VertexCompositeCandidate` and the CVH-corrected quantities
taken from the refit makers' per-candidate `ValueMap` products,
attached as external variables on the same table. Corrected columns
SHALL use names distinct from the raw columns.

#### Scenario: Raw and corrected columns coexist

- **WHEN** a candidate table is written for a channel
- **THEN** each row carries the raw candidate mass and kinematics
- **AND** the corresponding CVH-corrected mass and kinematics
- **AND** the corrected row count equals the raw candidate count.

### Requirement: Coverage of all persisted channels

The stream SHALL cover every persisted candidate channel across all
five AlCaReco streams, not a single channel.

#### Scenario: All channels produce tables

- **WHEN** the stream runs over each of the five AlCaReco streams
- **THEN** a candidate table is emitted for every persisted
  `VertexCompositeCandidate` collection in that stream.

### Requirement: Branch names match the histmaker raw contract

The stream SHALL emit raw per-candidate and per-daughter kinematic
branch names matching the contract already consumed by the `btojpsik`
histmaker's AlCaReco selection path, so that the existing histmaker
runs over this NanoAOD without a new input path.

#### Scenario: Histmaker runs unchanged

- **WHEN** the `btojpsik` histmaker is pointed at the produced NanoAOD
  with its AlCaReco selection path enabled
- **THEN** it resolves every branch it reads on that path
- **AND** requires no new code path to do so.

### Requirement: Track table with dE/dx and original-index columns

The stream SHALL emit a `Track` table from the AlCaReco
`vector<reco::Track>` collection, with the three dE/dx
`ValueMap<float>` estimators and the `originalIndex`
`ValueMap<unsigned int>` attached as columns on the same table.

#### Scenario: dE/dx columns present and keyed to tracks

- **WHEN** the `Track` table is written
- **THEN** each track row carries its `DeDxHarmonic2`,
  `DeDxAllHarmonic2`, `DeDxPixelHarmonic2`, and `originalIndex` values
  from the corresponding `ValueMap`s.

### Requirement: Primary-vertex table

The stream SHALL emit a `PV` table from `offlinePrimaryVertices`
carrying at least the vertex position, position uncertainties, and fit
quality (`chi2`, `ndof`).

#### Scenario: PV table row count matches the vertex collection

- **WHEN** the `PV` table is written
- **THEN** its row count per event equals the number of
  `offlinePrimaryVertices` in that event.

### Requirement: Muon table for the J/psi streams

The stream SHALL emit a `Muon` table for any AlCaReco stream that
persists a `reco::Muon` collection, exposing muon kinematics and the
identification flags used in the AlCaReco selection (`isGlobalMuon`,
`isTrackerMuon`, `numberOfMatches`).

#### Scenario: Muon table populated on a J/psi stream

- **WHEN** the config runs on a `TkAlJpsiX` or `TkAlJpsiMuMu` file
- **THEN** the `Muon` table row count per event equals the persisted
  `reco::Muon` collection size for that event.

### Requirement: HLT trigger decision columns

The stream SHALL expose the HLT trigger decisions from the AlCaReco
`TriggerResults` as event-level boolean columns.

#### Scenario: Trigger flags written

- **WHEN** the `Events` tree is written
- **THEN** it contains HLT decision columns derived from the AlCaReco
  HLT `TriggerResults`.
