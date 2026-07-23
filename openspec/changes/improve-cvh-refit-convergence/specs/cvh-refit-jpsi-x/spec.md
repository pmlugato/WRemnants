## ADDED Requirements

### Requirement: Configurable CVH joint-refit iteration cap and convergence threshold

`ResidualGlobalCorrectionMakerTwoTrackG4e` SHALL expose its Gauss–Newton iteration cap and EDM convergence threshold as ParameterSet inputs (`nIters`, default `10`; `edmConvergence`, default `1e-5`) instead of compile-time constants. The defaults SHALL bit-identically reproduce the published 81,948-event Run2016H baseline (`/ceph/.../cvh_outputs/run2016H_2a/joined_cvh_bplus_jpsik_plimit_0p05.root`) — any change to the defaults requires an explicit migration of the affected production configs.

#### Scenario: Default config matches published baseline
- **WHEN** the producer is instantiated with no `nIters` / `edmConvergence` overrides
- **THEN** the iteration loop uses `niters = 10` and the convergence test uses `edmval < 1e-5`, exactly as in the published baseline

#### Scenario: Operator overrides the iteration cap from the cfi
- **WHEN** the cfi sets `nIters = cms.uint32(100)` and `edmConvergence = cms.double(1e-2)`
- **THEN** the producer runs at most 100 Gauss–Newton iterations per icons phase and breaks out when `edmval < 1e-2`
- **AND** the values used are reported once per job at `LogInfo("ResidualGlobalCorrectionMakerTwoTrackG4e")`

#### Scenario: Bit-identical reproduction of baseline under explicit defaults
- **WHEN** the producer is instantiated with `nIters = 10` and `edmConvergence = 1e-5` explicitly set in the cfi
- **THEN** the per-event output (`Jpsi_*`, `Jpsicons_*`, `Mu*_*`, `chisqval`, `edmval`, `niter`) is bit-identical to the same producer with defaults

### Requirement: Selectable starting state for the joint refit linearisation

`ResidualGlobalCorrectionMakerTwoTrackG4e` SHALL accept a `useStartingState` ParameterSet input (string, default `"perigee"`) controlling the iteration-0 reference state used for the Gauss–Newton linearisation. The values `"perigee"` and `"midPropagated"` MUST be supported. `"perigee"` SHALL reproduce the current behaviour exactly (each muon's FTS taken at its own perigee). `"midPropagated"` SHALL propagate each muon's perigee FTS to the midpoint of the two perigee positions and use that as the iter-0 reference.

#### Scenario: Default perigee mode preserves baseline
- **WHEN** the producer is instantiated with no `useStartingState` override
- **THEN** the iter-0 reference state is taken at each muon's own perigee and the output is bit-identical to the published baseline

#### Scenario: Mid-propagated mode selects the alternative linearisation
- **WHEN** the cfi sets `useStartingState = cms.string("midPropagated")`
- **THEN** the iter-0 reference state for each muon is the result of propagating that muon's perigee FTS to the midpoint of the two perigee positions, and the GBL hessian, gradient, and per-iter updates are computed at that point

#### Scenario: Unknown starting-state value is rejected
- **WHEN** the cfi sets `useStartingState = cms.string("invalidValue")`
- **THEN** the producer throws an `edm::Exception` at construction time with a message listing the supported values

### Requirement: Per-iteration debug dump for convergence diagnostics

`ResidualGlobalCorrectionMakerTwoTrackG4e` SHALL accept a `debugPerIterDump` ParameterSet input (bool, default `false`). When `true`, the producer SHALL write per-iteration vector branches `chisqval_iter`, `edmval_iter`, `deltachisqval_iter`, `mu_qoverp_iter`, and `Jpsi_mass_iter` to the output TTree, recording the full sequence of values across both `icons=0` and `icons=1` phases for every event.

#### Scenario: Debug mode off by default
- **WHEN** the producer is instantiated with no `debugPerIterDump` override
- **THEN** the `*_iter` branches are NOT present in the output TTree

#### Scenario: Debug mode writes per-iter vectors
- **WHEN** the cfi sets `debugPerIterDump = cms.bool(True)` and an event with `niter = 7` is processed
- **THEN** the output TTree contains branches `chisqval_iter`, `edmval_iter`, `deltachisqval_iter`, `Jpsi_mass_iter` each with 7 entries for that event (one per iteration of the icons=1 phase), and `mu_qoverp_iter` with 14 entries (two muons × 7 iterations)

### Requirement: Per-muon final-hit count branches in joint refit output

`ResidualGlobalCorrectionMakerTwoTrackG4e` SHALL write per-muon counts of hits that survive the alignment-pass quality check as branches `Muplus_nvalidFinal`, `Muplus_nvalidpixelFinal`, `Muminus_nvalidFinal`, `Muminus_nvalidpixelFinal`. Each counter SHALL be incremented once per hit that satisfies the existing `morehitquality` condition (currently always true).

#### Scenario: Final-hit counters are symmetric with existing nvalid counters
- **WHEN** the producer runs on an event with `Muplus_nvalid = 16` and `Muplus_nvalidpixel = 2`
- **THEN** under the current `morehitquality = true` policy, `Muplus_nvalidFinal = 16` and `Muplus_nvalidpixelFinal = 2`
- **AND** if a future change introduces hit rejection in `morehitquality`, the `*Final` branches diverge from `nvalid` by exactly the count of rejected hits

### Requirement: Single-track CVH refit dead-branch cleanup

`ResidualGlobalCorrectionMakerG4e` (single-track CVH refit, used by the kaon leg) SHALL NOT ship branches whose value is always zero. The currently-declared `nValidHitsFinal` and `nValidPixelHitsFinal` branches in this producer SHALL either be incremented analogously to the two-track producer (symmetric implementation) or removed from the branch list entirely. Mixed-state declarations (declared, initialised, never incremented) are not allowed.

#### Scenario: Branch is meaningful or absent
- **WHEN** the joined output `joined_cvh_bplus_jpsik*.root` is read
- **THEN** `Kbach_nValidHitsFinal` either equals `Kbach_nValidHits` modulo per-hit rejection (option A) or is absent from the tree entirely (option B); the always-zero state observed in the published baseline SHALL NOT recur

### Requirement: Driver-level exposure of new CVH knobs

`CMSSW_15_0_19_patch2/src/Analysis/HitAnalyzer/test/runCvhBplusJpsiK.py` SHALL expose `nIters`, `edmConvergence`, `useStartingState`, and `debug` as VarParsing command-line arguments and forward them into the producer pset for the dimuon mode. The defaults SHALL match the producer defaults (`10`, `1e-5`, `"perigee"`, `False`).

#### Scenario: Defaults reproduce baseline
- **WHEN** `cmsRun runCvhBplusJpsiK.py mode=dimuon input=<file>` is run with no new knobs
- **THEN** the producer instantiation in the cmsRun config uses `nIters=10`, `edmConvergence=1e-5`, `useStartingState="perigee"`, `debugPerIterDump=False` and the output is bit-identical to the published baseline

#### Scenario: All four knobs forwarded
- **WHEN** `cmsRun runCvhBplusJpsiK.py mode=dimuon input=<file> nIters=100 edmConvergence=1e-2 useStartingState=midPropagated debug=True` is run
- **THEN** the producer pset contains the four overrides verbatim and the per-iter debug dump branches are written to the output TTree
