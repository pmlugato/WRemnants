# nbody-cvh-maker

## ADDED Requirements

### Requirement: Joint N-track CVH fit in a new plugin

The maker SHALL fit all N charged tracks of a candidate jointly through
the Geant4e/GBL machinery with a common vertex, for N of at least 3, and
SHALL be implemented as a new plugin that leaves the existing two-track
and single-track makers unmodified.

#### Scenario: Three-track candidate is fit jointly

- **WHEN** the maker runs on a B+ -> J/psi K+ candidate
- **THEN** all three tracks enter one joint fit with a common vertex
- **AND** the existing two-track maker's behaviour is unchanged.

#### Scenario: State dimension scales with track count

- **WHEN** the maker is configured for N tracks
- **THEN** the fitted state has 5 parameters per track plus 5 per hit.

### Requirement: Mass constraints on track subsets

The maker SHALL accept mass constraints applied to a specified subset of
a candidate's tracks, so that an intermediate resonance can be
constrained without constraining the whole candidate.

#### Scenario: J/psi constrained inside a three-track B candidate

- **WHEN** a J/psi mass constraint is configured on the two muon tracks
  of a three-track candidate
- **THEN** the fit constrains only that pair's invariant mass
- **AND** the bachelor track is unconstrained.

### Requirement: Fitted mass and its global-parameter derivative

The maker SHALL emit, per candidate, the fitted mass and its
uncertainty, and SHALL emit the derivative of that mass with respect to
the global correction parameters together with the corresponding global
parameter indices when the global-fit payload is enabled.

#### Scenario: Mass jacobian available for calibration

- **WHEN** the global-fit payload is enabled
- **THEN** each candidate carries the fitted mass, its global-parameter
  derivative, and the matching global indices
- **AND** those are sufficient to use the mass as a calibration
  observable.

### Requirement: Two-track closure

The maker SHALL, when configured for two tracks, reproduce the existing
two-track maker's fitted results on the same input events within a
documented tolerance.

#### Scenario: N=2 reproduces the two-track maker

- **WHEN** the maker is configured for two tracks on events the
  two-track maker has processed
- **THEN** the fitted masses and per-leg kinematics agree within the
  documented tolerance.
