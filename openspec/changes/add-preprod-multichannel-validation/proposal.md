# Change: Pre-production validation of J/psi+X, TkAlJpsiMuMu, D*->D0pi and V0 AlCaRecos across Charmonium and SingleMuon

## Why

We are one build away from the full 2016 postVFP launch. Before we
commit compute, we run one end-to-end validation pass on every
AlCaReco stream that will ship in production, on both input primary
datasets, from the same build. The validation runs on the exact
production configuration — no runtime flag, no debug-only branches.

The launch spans two PDs:

- **Charmonium** — `ALCARECOTkAlJpsiX` (8 channels) and
  `ALCARECOTkAlJpsiMuMu` (reference).
- **SingleMuon** — `ALCARECOTkAlDstToD0Pi`, `ALCARECOTkAlKShortTracks`,
  `ALCARECOTkAlLambdaTracks`.

Two structural issues in the current cff must be resolved before any
smoke or 100-file batch runs (both are pre-validation, both change
production behaviour permanently).

## What Changes

### Pre-validation refactors (must land BEFORE any smoke run)

**PV1. Delete preset A/B/C logic from `ALCARECOTkAlJpsiX_cff`.**

Preset B is the only shipped config. The current cff still reads
`TKALJPSIX_SELECTION_PRESET` from the environment and threads a
per-channel dictionary lookup through ~30 sites, with a late-mutation
stanza that only fires on preset `C`. The multi-preset comparison
finished long ago; keeping the switch alive is dead weight and forces
a reader to trace an env-var to know what any parameter actually is.

- Delete `_TKALJPSIX_SELECTION_PRESET`, the env-var reads,
  `_NON_V0_PRESETS`, `_PRESETC_KALMAN`, and the `preset == 'C'`
  mutation stanza.
- Inline preset B numeric literals directly into each per-channel
  producer clone.
- Purge `TKALJPSIX_SELECTION_PRESET` from the Condor tarball
  builder, recoskim wrappers, and memory / spec text.
- No C++ change.

**PV2. Fix the loose-muon event gate.**

The current sequence gates the whole path on the tight-muon
collection: `ALCARECOTkAlJpsiXGoodMuons` is an `EDFilter` with
`filter = True`, so a zero-tight-muon event stops the path before
`ALCARECOTkAlJpsiXLooseMuons` and `ALCARECOTkAlJpsiXJpsiOnlyCandidates`
ever run. This structurally excludes the population the J/psi-only
channel was created to recover: J/psi events with both muons
tracker-only. Only mixed-pair events (one tight leg + one loose leg)
survive.

Fix — move the gate from muon-level to candidate-level:

- Set `ALCARECOTkAlJpsiXGoodMuons.filter = False`. It still runs
  and still emits the tight muon collection consumed by the 7
  bachelor channels; it just no longer gates the path.
- Keep `ALCARECOTkAlJpsiXLooseMuons.filter = False` (already).
- Merge the 8 per-channel VCC outputs into
  `ALCARECOTkAlJpsiXAnyCandidate` (stock
  `VertexCompositeCandidateCollectionMerger` if available, otherwise
  a small custom `Merger<>` plugin).
- Insert `ALCARECOTkAlJpsiXCandFilter`
  (`CandViewCountFilter`, `minNumber = 1`) after the merger and
  before `CompositeDaughterTrackProducer`.

Rate impact must be measured before the switch lands (task 0.6).
The change only widens acceptance; it does not modify per-event
content on events already accepted by the tight gate.

### Landed upstream (David's `444b20f`) — pre-conditions, not proposed

- `signedPdgId` convention: lepton pdgId sign is opposite to
  charge (mu- carries +13); hadron pdgId follows track charge
  (K+ carries +321). Downstream readers of `daughter->pdgId()`
  rely on this.
- `JpsiXCandidateProducer` VCC-mode X-candidate veto on tracks
  reused from the J/psi muons.
- `VertexCompositeCandidateRemapper` now carries vertex
  covariance and fit chi2/ndof into the remapped candidates —
  enables plotting vertex quality per channel from the persisted
  VCC.
- `TwoBodyDecayCandidateProducer` guards `VertexException` from
  Kalman fits.
- `AlignmentThreeBodyDecayTrackSelector` end-of-job summary via
  `MessageLogger`.
- ψ(2S) V0-mode stale-comment cleanup in the cff.

### Stream inventory that must run and be validated

- **`ALCARECOTkAlJpsiX`** (Charmonium, 8 channels). Post-PV1/PV2.
- **`ALCARECOTkAlJpsiMuMu`** (Charmonium). Standard CMSSW stream,
  tight muon selector. Used as the tight-muon reference for the
  loose-vs-tight muon comparison.
- **`ALCARECOTkAlDstToD0Pi`** (SingleMuon). New under our build
  from the cherry-pick. No HLT prefilter.
- **`ALCARECOTkAlKShortTracks`** and **`ALCARECOTkAlLambdaTracks`**
  (SingleMuon). Sharing the local loosened V0 clone
  (`ALCARECOTkAlV0Candidates`). No HLT prefilter.

### Stream persist-alignment (production-permanent)

Every stream shall persist, in production, the same categories of
output where the physics allows it. Today the four streams diverge:
J/psi+X keeps candidates + tracks + dE/dx + muon selector output +
track-muon association; TkAlJpsiMuMu keeps only candidates + tracks;
D* keeps tracks + dE/dx but no standalone candidate collection; V0
streams to be audited after the sparse-checkout add.

Alignment target for production output:

| Category                    | J/psi+X | JpsiMuMu | D*        | KS/Λ    |
| --------------------------- | ------- | -------- | --------- | ------- |
| Cloned tracks + originalIdx | ✓       | ✓        | ✓         | audit   |
| dE/dx (Harmonic2 x 3)       | ✓       | **add**  | ✓         | audit   |
| Candidate VCC               | ✓ x 8   | ✓        | **add**   | audit   |
| Muon selector output        | ✓ loose | **add**  | N/A       | N/A     |
| Track -> muon association   | ✓       | **add**  | N/A       | N/A     |
| Standard L1/HLT/DCS/PV      | ✓       | ✓        | ✓         | audit   |

"add" means a `keep` line + (for D*) a small candidate producer,
committed to the production `_cff.py` / `_Output_cff.py`. No runtime
flag, no diagnostic mode — the production stream always emits these
starting from this change.

Rationale for each addition:

- **JpsiMuMu dE/dx maps**: without them the tight sample cannot be
  compared to the J/psi+X sample on the dE/dx axis. Cheap to add
  (three `keep` lines + three `DeDxValueMapProjector` clones).
- **JpsiMuMu tight-muon collection**: needed for the tight-vs-loose
  muon comparison. Single `keep *_ALCARECOTkAlJpsiMuMuGoodMuons_*_*`.
- **JpsiMuMu track -> muon association**: needed to key persisted
  muons to persisted tracks for offline analyses that want that
  cross-reference (same pattern as J/psi+X).
- **D\* candidate VCC**: the current D* stream keeps only the
  cloned tracks that survived the internal `ThreeBodyDecaySelector`.
  There is no persisted D* candidate object to plot D* mass /
  Q-value / vertex chi2 from. Add a small candidate emitter that
  runs the same selection logic and writes a
  `VertexCompositeCandidateCollection` for D* (with the D0 as
  its intermediate accessible via `daughter()->daughter()`).

### Validation deliverables (read from production output only)

- **Loose-vs-tight muon comparison** (David's ask): per-observable
  overlay of TkAlJpsiMuMu tight muons and J/psi+X loose muons on
  the same events. Reads the tight-muon and loose-muon collections
  directly.
- **Per-channel kinematics** for all 11 candidate streams (8
  J/psi+X + D* + KS + Λ): mass, pT, eta, y, helicity angle,
  vertex chi2/ndof. Reads the persisted VCCs.
- **Cut-efficiency numbers**: parse the `LogInfo` counter lines
  that `JpsiXCandidateProducer` and `TwoBodyDecayCandidateProducer`
  already emit to cmsRun stderr on every job. No new persisted
  content required.
- **V0 loosening yield**: one-shot offline re-run of the
  CMSSW-default V0 producer on the same input files as a
  side-study, not part of the AlCaReco. Records the mean
  loosened / default yield ratio per species.
- **dE/dx distributions** per stream (from the maps kept by each
  stream after PV alignment).
- **Resource scatter** (runtime, peak RSS, output size / event)
  from the 100-file batch, split by PD.
- **Bit-invariance check**: preset B J/psi+X output on events
  accepted by both the old tight gate and the new candidate gate
  matches the last smoke's per-channel candidate counts and
  persisted track counts.

### Smoke tests and 100-file batch

- **Smoke**: 1 file per PD, ~1000 events. Charmonium job schedules
  J/psi+X + JpsiMuMu in one cmsRun; SingleMuon job schedules
  D* + KS + Λ in one cmsRun.
- **100-file batch**: two Condor clusters (one per PD) on the T3,
  same tarball, shared `run.sh` dispatching on recoskim filename.

### Full 2016 postVFP launch (planned, not implemented here)

Two Condor batches sharing one payload tarball:

- Charmonium: ~27k jobs on ~27k files.
- SingleMuon: N jobs, N measured by the DAS query in task 1.2.

## Impact

- Affected specs:
  - `alcareco-jpsi-x` — preset removal, candidate-gate event
    filter, and new persist alignments.
  - `alcareco-jpsi-mumu` (new) — persist alignments (dE/dx,
    tight muons, track-muon association).
  - `alcareco-dstar-d0pi` — new candidate VCC persist.
  - `alcareco-v0` — audit and align KShort / Lambda outputs.
  - `preprod-validation` — smoke, 100-file batch, plotting
    deliverables, projections.
- Affected code:
  - `CMSSW_10_6_17_patch1/src/Alignment/CommonAlignmentProducer/python/ALCARECOTkAlJpsiX_cff.py`
    — preset removal, muon gate demotion, candidate merger +
    filter.
  - `CMSSW_10_6_17_patch1/src/Alignment/CommonAlignmentProducer/python/ALCARECOTkAlJpsiMuMu_cff.py`
    (+ `_Output_cff.py`) — dE/dx projectors, muon collection
    persist, association producer.
  - `CMSSW_10_6_17_patch1/src/Alignment/CommonAlignmentProducer/python/ALCARECOTkAlDstToD0Pi_cff.py`
    (+ `_Output_cff.py`) — D* candidate emitter.
  - Sparse-add and audit
    `CMSSW_10_6_17_patch1/src/Alignment/CommonAlignmentProducer/python/ALCARECOTkAlKShortTracks_cff.py`
    and `ALCARECOTkAlLambdaTracks_cff.py`; align keeps.
  - `condor/multichannel_alcareco/` — new directory with
    tarball builder, per-recoskim wrapper, two `.sub` files, and
    a shared `run.sh`.
  - `scripts/btojpsik/plots_preprod_validation.py` — new
    plotting driver.
- No changes in `narf`, `rabbit`, `wremnants/`.

## Non-goals

- Stage-2 CVH refit for D* / V0. Deferred.
- Retuning the D* mass window, Q-value cut, or V0 clone cuts.
- Any behavioural change to `TkAlJpsiMuMu` beyond persisting
  additional collections for alignment with J/psi+X.
