## ADDED Requirements

### Requirement: Pilot fixes reconciled into our tree

The MC production tree at `condor/mc_inclusive_btojpsix_2016postvfp/` SHALL
contain david_w's `pietro_fixes` (five items, per
`/work/submit/david_w/ZMass/BtoJpsiX_MCprod/pietro_fixes/PIETRO_FIXES_README.md`)
with the physics/tuning content bit-identical to `fixed_files/`, and SHALL
have `DEST_DIR` retargeted to `/ceph/submit/data/group/cms/store/mc/...`
before the full ~15 fb⁻¹ inclusive B→J/ψ+X MC submission is issued. David's
operator identity (proxy path, `+AccountingGroup`, transient site fences)
remains authoritative for the actual submission and is not overwritten by
this reconciliation.

#### Scenario: Physics fixes present, operator identity untouched
- **WHEN** the reconciliation described in tasks §1 is complete
- **THEN** `diff` against `pietro_fixes/fixed_files/`, excluding
  operator-identity lines and `DEST_DIR`, is empty
- **AND** operator-identity fields are unchanged on either side

#### Scenario: DEST_DIR retargeted to group MC store
- **WHEN** the reconciliation is complete
- **THEN** `DEST_DIR` in `run.sh`, `SUBMIT.sh`, `RECONCILE.sh`, `quota.sh`,
  and `HANDOFF_README.md` points at
  `/ceph/submit/data/group/cms/store/mc/inclusive_btojpsix_2016postvfp/`
  (or the equivalent agreed subpath), and the operator note for David
  is drafted for delivery after pilot approval

### Requirement: Pilot physics-health gate before full submission

The full ~15 fb⁻¹ MC production SHALL NOT be submitted until the pilot batch
(~108 jobs at
`/ceph/submit/data/user/d/david_w/mz/mc/inclusive_btojpsix_2016postvfp/`)
has been validated end-to-end through Stage-2 CVH refit + joint J/ψ+K
kinematic-vertex fit, and the user has approved the written
`mc_pilot_validation_report.md`.

#### Scenario: Stage-2 CVH ditrack J/ψ sanity met
- **WHEN** the pilot outputs are run through Stage-2 CVH with
  `applyIdealGeometry=False`
- **THEN** the ditrack J/ψ mass distribution shows a median near 3.08 GeV
  with the bulk of events inside the [3.0, 3.2] GeV core and zero
  propagation aborts (matching david_w's expectation for a correctly
  split=1 sample)

#### Scenario: B⁺ candidate formed via joint kinematic-vertex fit
- **WHEN** the CVH two-track J/ψ and CVH single-track K± are combined
  into a B⁺ candidate
- **THEN** an offline joint kinematic-vertex fit on the CVH-refit tracks
  produces a per-candidate χ²/ndof used as the primary selection knob
- **AND** the χ²/ndof cut visibly suppresses combinatoric background in
  the m(B⁺) distribution

#### Scenario: Gen-matched signal peak identified
- **WHEN** reco B⁺ candidates are split by gen matching (ΔR < 0.03 at both
  muons and the bachelor)
- **THEN** the m(B⁺) peak in the [5.25, 5.31] GeV window sits in the
  gen-matched sample, with the non-matched sample forming a flat
  combinatoric baseline

#### Scenario: MC-vs-data overlay produced
- **WHEN** the pilot MC and the data ALCARECOs at
  `/ceph/submit/data/user/p/pmlugato/mz/alcareco/full_2016postvfp/charmonium/Run2016*`
  have both been processed through the same scanner + joint-fit +
  selection stack
- **THEN** m(J/ψ) and m(B⁺) shape overlays are included in the report as
  a visual sanity check (no numerical agreement threshold applied)

#### Scenario: Go/no-go recommendation signed off
- **WHEN** the pilot validation report is drafted and the physics-health
  checks above are met
- **THEN** the user signs off on the report before the full ~15 fb⁻¹
  submission is issued to David
