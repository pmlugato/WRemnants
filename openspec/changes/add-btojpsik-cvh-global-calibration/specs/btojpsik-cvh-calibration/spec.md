## ADDED Requirements

### Requirement: B⁺→J/ψK⁺ CVH makers emit gradients and Hessians for calibration

The B⁺→J/ψK⁺ Stage-2 CVH makers SHALL be configurable to write per-candidate
gradient and Hessian sidecar output suitable for global-correction aggregation.
The dimuon maker (`globalCorJpsiK`, mass-constraint ON) and the bachelor-kaon
single-track maker (`globalCorJpsiKKaon`) SHALL each expose `fillGrads`,
`fillGradsFactored`, and `fillRunTree` knobs through `runCvhBplusJpsiK.py`.
The defaults SHALL remain OFF so the existing diagnostic offline-join behavior
is bit-identical unless calibration output is explicitly requested. When
calibration output is enabled, both makers SHALL use the grouped global
material model (`globalMaterialModel=True` + a shared `materialGroupsFile`,
with `skipHitlessSurfaces=True`) so material is fit as per-group scales
(parmtype 15).

#### Scenario: Calibration output enabled on the kaon side

- **WHEN** `runCvhBplusJpsiK.py` is run with `mode=kaon` and the grads flags enabled
- **THEN** the kaon sidecar tree contains `nParms`, `globalidxv[nParms]`,
  `gradv[nParms]`, and a factored Hessian (`hessfactorv`) or packed Hessian
  (`hesspackedv`), plus `run`/`lumi`/`event` and per-candidate quality branches

#### Scenario: Defaults preserve diagnostic-join behavior

- **WHEN** `runCvhBplusJpsiK.py` is run without the grads flags
- **THEN** `fillGrads` and `fillRunTree` remain `False` and the sidecar output
  is unchanged from the current diagnostic-join configuration

### Requirement: Parameter-catalog consistency across makers

Before any grad/Hess accumulation is trusted, the system SHALL verify that the
`runtree` global-parameter catalog is identical across the baseline
`TkAlJpsiMuMu` maker, the B⁺→J/ψK⁺ dimuon maker, and the B⁺→J/ψK⁺ kaon maker.
The check SHALL compare the number of parameters and the ordered
`(parmtype, subdet, layer, rawdetid)` tuples. A mismatch SHALL be a hard error
that stops the pipeline.

#### Scenario: Catalogs match

- **WHEN** the baseline and B⁺→J/ψK⁺ productions share the same CMSSW build,
  field-mode set, and `materialGroups` file
- **THEN** the catalog-consistency check passes and aggregation may proceed

#### Scenario: Catalogs mismatch

- **WHEN** any maker's `runtree` differs in length or parameter ordering
- **THEN** the check raises a hard error and no combined solve is produced

### Requirement: J/ψ-consistency sanity check

The system SHALL produce two dimuon-only (no kaon) global-correction fits — the
nominal mass-constrained J/ψ→μμ method on `TkAlJpsiMuMu_*.root` and the
J/ψ-from-B dimuon-only fit on the B⁺→J/ψK⁺ candidates — for Run2016H (subset
permitted), and SHALL confirm their b-field/material corrections agree within
uncertainties. This validates that the B dimuons behave like the standard J/ψ
sample before the within-B kaon comparison is trusted.

#### Scenario: Both J/ψ calibrations produced

- **WHEN** the `TkAlJpsiMuMu` and J/ψ-from-B dimuon-only CVH grads are each
  aggregated and solved with the shared prior + NMR-constraint configuration
- **THEN** a `correctionResults` file is written for each, containing the
  fitted parameters and their uncertainties

#### Scenario: J/ψ consistency confirmed

- **WHEN** the two dimuon-only fits are compared
- **THEN** their b-field and material corrections are tabulated side by side
  and their agreement within uncertainties is reported

### Requirement: B⁺→J/ψK⁺ kaon aggregation into the shared parameter space

The system SHALL aggregate the B⁺→J/ψK⁺ bachelor-kaon single-track grad/Hess
into the same `nparms` global index space used by the dimuon side, applying
kaon-appropriate quality filters (valid-hit and pixel-hit counts, χ²/ndof, EDM
convergence) rather than the dimuon-specific filters. The kaon side SHALL be
run with a low-momentum propagation floor (`plimit=0.05`) and the surviving
candidate fraction SHALL be reported.

#### Scenario: Kaon grads aggregated

- **WHEN** the kaon sidecar files are aggregated
- **THEN** a `combinedgrads` HDF5 is produced over the same parameter catalog
  as the dimuon side, and the number of accepted kaon candidates is logged

### Requirement: Within-B dimuon-only vs dimuon+kaon fits

The system SHALL produce two global-correction fits over the same B⁺→J/ψK⁺
events: (i) a dimuon-only fit from the jpsicons dimuon grad/Hess, and (ii) a
dimuon+kaon fit formed by summing the dimuon grad/Hess with the kaon grad/Hess
over the shared index space. Both fits SHALL use an identical prior and
NMR-constraint configuration so that only the kaon contribution differs between
them.

#### Scenario: Dimuon-only fit produced

- **WHEN** the B⁺→J/ψK⁺ dimuon grads are aggregated and solved
- **THEN** a dimuon-only `correctionResults` file is written

#### Scenario: Dimuon+kaon fit produced

- **WHEN** the dimuon and kaon aggregated files are summed and solved with the
  same configuration
- **THEN** the total gradient and Hessian are the element-wise sums over the
  shared index space and a dimuon+kaon `correctionResults` file is written

### Requirement: Report b-field and material corrections, uncertainties, and correlation

The system SHALL report, for the dimuon-only fit and the dimuon+kaon fit, the
b-field block (parmtype 14 = B-field scalar-potential modes) and the material
block (parmtype 15 = per-group material scale). The b-field report SHALL
separate **mode 0 (the absolute scale)** from **modes 1 and above (the field
shape)**, because a soft displaced hadron is expected to act on the shape rather
than the scale. For each block the report SHALL include the per-parameter
constraint σ, the fitted magnitude, and the b-field↔material correlation (the
parmtype-14 × parmtype-15 sub-covariance) as a diagnostic. The comparison SHALL
make the change from dimuon-only to dimuon+kaon visible.

#### Scenario: Covariance enabled

- **WHEN** the solve is run for reporting
- **THEN** the covariance (at least the diagonal and the b-field↔material block
  sub-covariance) is computed and non-zero uncertainties are returned

#### Scenario: Dimuon-only-vs-dimuon+kaon report

- **WHEN** both fits are available
- **THEN** the per-parameter constraint σ, magnitudes and the field↔material
  correlation are tabulated/plotted side by side for dimuon-only and
  dimuon+kaon, with the field scale (mode 0) reported separately from the field
  shape (modes 1+)

### Requirement: Alignment held fixed as the nominal configuration

The nominal fits SHALL float only the 92 global parameters (parmtype 14 + 15),
holding the alignment block (parmtypes 0-5) fixed, matching how the reference
fit is run. The alignment-marginalised variant SHALL also be produced, as the
honest figure for any absolute uncertainty statement.

#### Scenario: Nominal fit floats only the global parameters

- **WHEN** a nominal fit is solved
- **THEN** the alignment parameters are decoupled from the system and the
  number held fixed is logged

#### Scenario: Marginalised variant available

- **WHEN** the alignment-marginalised variant is solved on the same grads
- **THEN** its b-field/material uncertainties are reported alongside the
  alignment-fixed ones

### Requirement: Species closure check via the muon mass hypothesis

The system SHALL rerun the bachelor-kaon CVH refit on the same tracks under the
*muon* mass hypothesis (`kaonAsMuon=True`) and compare the resulting material
block against the kaon-hypothesis result. Because the track sample is identical
and only the assumed dE/dx changes, the difference isolates the species /
energy-loss effect.

#### Scenario: Muon-hypothesis variant produced and compared

- **WHEN** the bachelor tracks are refit with `trackParticleName='mu'` and
  aggregated and solved with the same configuration
- **THEN** the material group scalings from the kaon and muon hypotheses are
  reported side by side, and any shift is attributed to the energy-loss model

### Requirement: Documentation slides

Once the study is complete, the system SHALL produce a detailed slide deck
documenting the whole procedure: the three-stage pipeline, the field-mode /
grouped-material parameterization, the iterations and decisions taken, the
J/ψ-consistency sanity check, and the final results.

#### Scenario: Slides produced

- **WHEN** the dimuon-only, dimuon+kaon, and sanity-check results are available
- **THEN** a slide deck is produced that walks through the procedure,
  iterations, and results (including the b-field/material magnitude,
  uncertainty, and correlation comparison)
