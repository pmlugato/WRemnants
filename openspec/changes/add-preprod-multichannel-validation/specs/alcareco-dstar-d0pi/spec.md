## ADDED Requirements

### Requirement: TkAlDstToD0Pi AlCaReco Runs Under Our Build

`ALCARECOTkAlDstToD0Pi` SHALL run under `CMSSW_10_6_17_patch1` after
the cherry-pick. The Stage-1 sequence SHALL execute the DCS filter,
the `AlignmentTrackSelectorWithIndexMap` running the
`AlignmentThreeBodyDecayTrackSelector` (D* -> D0(K pi) pi_s with
K = 0.493677 GeV, pi = 0.139570 GeV, D0 window [1.75, 1.98] GeV,
D* window [1.89, 2.13] GeV, mass difference [0.140, 0.152] GeV,
`applyChargeFilter = True` with `charge = 1` and
`useUnsignedCharge = True`), the `alcaDedxJointEstimator`, and the
three dE/dx projectors. It SHALL persist the outputs listed in
`ALCARECOTkAlDstToD0Pi_Output_cff.py`.

#### Scenario: Non-empty output on SingleMuon RAW

- **WHEN** the D* recoskim runs over 1000 events of a Run2016H
  SingleMuon RAW file
- **THEN** the persisted collections (cloned tracks, three dE/dx
  maps, the new D* candidate VCC below) are all present and
  non-empty, and cmsRun exits 0

### Requirement: D* Candidate VCC Persisted

The stream SHALL persist a `VertexCompositeCandidateCollection` for
D* candidates. Each candidate SHALL expose the D* as the mother
with the D0 as an intermediate daughter (accessible via
`daughter()->daughter()`) so that D* mass, D0 mass, Q-value, and
vertex quality can be computed from the persisted candidate alone.
This is production-permanent — the stream always emits this
collection.

#### Scenario: Candidate content matches selector-surviving triplets

- **WHEN** the D* recoskim runs
- **THEN** the persisted D* candidate collection contains exactly
  the triplets whose tracks passed the internal three-body
  selector, and each candidate's daughter chain resolves to
  (D*, D0 = (K, pi), pi_s) with `signedPdgId` conventions
  (K carries +321 for a positive-charge track)

### Requirement: Cut Values Preserved As Merged

The D* cff SHALL preserve the numeric cut values as merged from
the cherry-pick: `trackQualities = ["highPurity"]`, `ptMin = 0.35`,
`etaMin/etaMax = -3.5 / 3.5`, `nHitMin = 0`, D* mass window
`[1.89, 2.13]` GeV, D0 mass window `[1.75, 1.98]` GeV, mass
difference `[0.140, 0.152]` GeV, `applyChargeFilter = True` with
`charge = 1` and `useUnsignedCharge = True`.

#### Scenario: Cut values unchanged

- **WHEN** the cff is loaded
- **THEN** each numeric cut parameter equals the merged value
