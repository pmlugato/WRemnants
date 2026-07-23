## ADDED Requirements

### Requirement: Composed B+ → J/ψ K+ Stage-2 CVH refit

Stage-2 for the B+ → J/ψ K+ channel SHALL refit the candidate's three leaf
tracks by composing two existing CVH kernels rather than a single
multi-track fit: the μ⁺μ⁻ pair SHALL be refit with
`ResidualGlobalCorrectionMakerTwoTrackG4e` under the muon mass hypothesis,
and the bachelor kaon SHALL be refit with the single-track
`ResidualGlobalCorrectionMakerG4e` configured with
`trackParticleName="kaon"`. No 3-track joint fitter SHALL be introduced.
Neither the B+ mass constraint nor the B+ vertex constraint SHALL be
applied inside Stage-2, and no J/ψ → μμ mass constraint SHALL be applied
inside the CVH joint fit (`doMassConstraint=False` on the two-track maker
cfi); the common-vertex coupling of the two muons inside
`ResidualGlobalCorrectionMakerTwoTrackG4e` is structural to the joint
state and is NOT a mass constraint. The kaon SHALL be refit independently
at its own PCA, with no link to the dimuon vertex inside Stage-2. The two
CVH fit kernels SHALL NOT be modified in their fitting logic. The bachelor
leg MUST be propagated with the kaon mass so that its energy-loss /
material Jacobian is correct for a kaon.

#### Scenario: Bachelor refit uses the kaon hypothesis
- **WHEN** the single-track maker runs on the bachelor track of a B+ → J/ψ K+ candidate
- **THEN** it is configured with `trackParticleName="kaon"`, the Geant4e propagation uses the kaon Geant4 particle name and the PDG kaon mass from `ParticleProperties`, and the emitted per-track Jacobian reflects the kaon energy-loss model (not the default muon)

#### Scenario: Bachelor hypothesis is taken from the collection name, not the daughter
- **WHEN** the Stage-2 chain runs on an `ALCARECOTkAlJpsiXBPlusResonances` collection
- **THEN** the splitter and the single-track maker hardcode the kaon hypothesis from the collection name and do NOT read `daughter(1).mass()` or `daughter(1).pdgId()`; the same code therefore works whether the upstream producer is storing the diagnostic mass/pdgId stamps or has dropped them (final postVFP run)

#### Scenario: No kinematic constraint applied in Stage-2
- **WHEN** the two legs are refit
- **THEN** the two-track CVH maker runs with `doMassConstraint=False` (no J/ψ mass constraint); no B+ mass constraint and no B+ vertex constraint couple the muon pair to the bachelor kaon; the kaon is refit independently at its own PCA; the B+ vertex/mass constraint is left for the downstream GBL fit (out of scope here)

#### Scenario: Fit kernels reused unmodified
- **WHEN** the Stage-2 B+ chain is built
- **THEN** the muon-pair and bachelor refits are performed by the existing `ResidualGlobalCorrectionMakerTwoTrackG4e` and `ResidualGlobalCorrectionMakerG4e` kernels with no change to their fitting code; the only new C++ is the candidate splitter and an additive per-candidate index

---

### Requirement: B+ candidate splitter

A new `JpsiKCandidateSplitter` EDProducer SHALL read the Stage-1 B+
`VertexCompositeCandidate` collection and, for each B+ candidate, emit (1) a
`reco::VertexCompositeCandidateCollection` containing the J/ψ (dimuon)
sub-candidate and (2) a collection of the bachelor kaon (the track or
`RecoChargedCandidate`) consumable by the single-track maker. Each emitted
dimuon and each emitted bachelor SHALL carry a B+ candidate index that
uniquely identifies which B+ candidate it came from within the event. The
splitter SHALL perform no vertex fit and no kinematic fit. The splitter
SHALL skip a B+ candidate (rather than abort the job) if its expected
daughter layout (a J/ψ VCC daughter plus a bachelor `RecoChargedCandidate`
daughter) is not present.

#### Scenario: One dimuon and one bachelor emitted per B+ candidate
- **WHEN** the splitter processes a B+ candidate whose daughters are a J/ψ VCC and a bachelor kaon `RecoChargedCandidate`
- **THEN** it emits exactly one dimuon sub-candidate and one bachelor element, both tagged with the same B+ candidate index for that event

#### Scenario: Malformed candidate is skipped, not fatal
- **WHEN** a B+ candidate does not expose the expected J/ψ-VCC + bachelor daughter layout
- **THEN** the splitter skips that candidate and continues processing the remaining candidates without raising an exception

---

### Requirement: Joinable per-candidate output for the GBL fit

The two refit makers SHALL each propagate the B+ candidate index into their
output so that the two output TTrees can be joined 1:1 offline on
`(run, lumiBlock, event, bCandIdx)`. An offline combination step SHALL
produce one flat ROOT TTree with one row per B+ candidate, carrying: the
refit parameters at PCA, full covariance, and per-track Jacobians for all
three tracks (two muons + bachelor kaon); the raw B+ kinematics; and the
m(μμK) invariant mass computed under the kaon hypothesis. The output format
SHALL be the existing flat-TTree CVH format (the GBL-fit input), NOT custom
NanoAOD. The per-candidate index plumbing SHALL be additive: existing
J/ψ → μμ, Υ → μμ, Z → μμ, Ks, Λ, and D0 maker configurations SHALL behave
identically when no candidate-index input is supplied.

#### Scenario: Two refit trees join 1:1
- **WHEN** the offline combination step runs on the two-track and single-track output trees from one Stage-2 job
- **THEN** each B+ candidate's two muons and bachelor kaon are matched into a single output row by `(run, lumiBlock, event, bCandIdx)`, and a candidate whose muon-pair refit or bachelor refit failed is reported (not silently merged with another candidate)

#### Scenario: Existing configs unaffected by the index plumbing
- **WHEN** the J/ψ → μμ production config (`runCvhJpsi.py`) is run after the per-candidate-index change is added
- **THEN** its output is unchanged from before the change (the index branch is absent or defaulted, and no other behaviour differs)

#### Scenario: Output carries all three tracks' Jacobians
- **WHEN** a B+ candidate is refit and combined
- **THEN** the output row contains the per-track Jacobians ∂(track params)/∂a for both muons and the bachelor kaon, the full per-track covariances, and m(μμK) under the kaon hypothesis

---

### Requirement: Cross-release Stage-2 driver

A cmsRun driver SHALL run the Stage-2 B+ chain (splitter → two-track maker →
single-track maker) in `CMSSW_15_0_19_patch2`, reading the finalized
preset-B 2016 `TkAlJpsiX` ALCARECO produced in `CMSSW_10_6_17_patch1`
(`preset_B_final/Run2016H/`) re-written at branch split-level 1 (the
`alcarecoSplitLevel` customise, commit `b8ebaa3d9f2`, MUST be applied to
the production cfg before Stage-2 can read the inputs; otherwise the
SiStripCluster v11→v14 I/O rule does not fire on split branches and
clusters silently return empty, breaking the CVH refit). The driver SHALL
use the `Run2_2016` era and a Run-2 data global tag, the scalar-potential
B-field model, and SHALL store per-track Jacobians (`fillJac=True`). It
SHALL be modeled on the existing `runCvhJpsi.py` and SHALL NOT require a
10_6_x release for Stage-2.

#### Scenario: Driver reads 2016 ALCARECO in the 15_0_19 release
- **WHEN** the driver is run on a finalized preset-B 2016 `TkAlJpsiX` ALCARECO file re-written at splitLevel=1
- **THEN** the B+ `VertexCompositeCandidate` collection deserializes, SiStripCluster bytes load non-empty through the I/O-rule chain, the splitter consumes the candidates, and the chain produces the two refit trees with Jacobians stored — without porting Stage-2 to a 10_6_x release

#### Scenario: SplitLevel=99 input is rejected with a diagnostic
- **WHEN** the driver is run on a `TkAlJpsiX` ALCARECO file whose `SiStripCluster` branch has split level 99 (no `setAlcaRecoSplitLevel` customise applied at production time)
- **THEN** the driver SHOULD detect the condition (e.g. via a startup branch-split-level probe or by counting empty hit collections in the first event) and emit a clear error pointing at the splitLevel customise, rather than silently producing refit trees with empty hits
