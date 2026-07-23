## ADDED Requirements

### Requirement: m32 Fit Observable
The histmaker SHALL accept a `--histToFit` flag that selects the fit observable used for
`nominal_HistToFit`. The default `mB` reproduces the current B-mass observable. The
`m32` option uses Q = m(B, vertex) − m(J/ψ, kinematic) as the fit observable, where
both masses are unconstrained (nomc-consistent).

#### Scenario: default mB is unchanged
- **WHEN** `btojpsik.py` is invoked without `--histToFit`
- **THEN** `nominal_HistToFit` uses `bkmm_{vertex_prefix}_mass` as its observable axis,
  identical to pre-change behaviour

#### Scenario: m32 observable selected
- **WHEN** `btojpsik.py` is invoked with `--histToFit m32`
- **THEN** `nominal_HistToFit` uses `bkmm_{vertex_prefix}_m32` as its observable axis,
  whose values are `bkmm_{vertex_prefix}_mass_scalar − mm_kin_mass[bkmm_mm_index[bkmm_best_idx]]`

#### Scenario: m32 axis range covers signal peak
- **WHEN** `--histToFit m32` is active
- **THEN** the fit observable axis spans [1.8, 2.6] GeV with 100 bins, covering the
  signal Q-value of approximately 2.182 GeV

#### Scenario: m32 with --nomc is the primary use case
- **WHEN** `--histToFit m32 --nomc` are both active
- **THEN** both the B mass (`bkmm_nomc_mass`) and the dimuon mass (`mm_kin_mass`) are
  unconstrained, making the Q-value self-consistent

#### Scenario: kaon scale variations unaffected
- **WHEN** `--histToFit m32` and `--includeKaonScaleVariations` are both active
- **THEN** `nominal_muonScaleSyst_responseWeights` is produced using the m32 fit axes
  (same pt/eta/charge axes; only the mass axis changes)
