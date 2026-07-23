## ADDED Requirements

### Requirement: Kaon energy in the e-term of A/e/M uncertainty weights
The btojpsik uncertainty computation SHALL use the kaon total energy E_K = sqrt(p_K² + m_K²)
via the factor E_K/p_K = sqrt(1 + m_K²/p_K²) in the `e` (energy-loss-in-material) term of
`calculateQopUnc` when producing A/e/M alternative weights for kaon scale variations. The A
and M terms SHALL remain unchanged. The kaon mass SHALL be m_K = 0.49368 GeV. The muon
uncertainty computation with particle_mass = 0 SHALL produce identical results to pre-change.

#### Scenario: e-term is energy-corrected for kaons
- **WHEN** `JpsiCorrectionsUncHelperSplines` is constructed with `particle_mass = 0.49368`
- **AND** `calculateQopUnc(recPt, recEta, recCharge, AUnc, eUnc, MUnc)` is called
- **THEN** the `kUnc` expression uses `eUnc * k * E_over_p` where `E_over_p = sqrt(1 + m_K²/p_total²)`
- **AND** the A and M contributions are unchanged

#### Scenario: Massless limit preserves muon behaviour
- **WHEN** `JpsiCorrectionsUncHelperSplines` is constructed with default `particle_mass = 0`
- **THEN** `E_over_p = 1.0` and `calculateQopUnc` output is identical to the pre-change result

#### Scenario: e-variation alternative weights differ for low-pT kaons
- **WHEN** kaon pT is small (e.g. pT ~ 1–2 GeV) where m_K/p_K is non-negligible
- **THEN** the e-term alternative weight deviates from the massless approximation by a
  factor ≈ E_K/p_K > 1
- **AND** A and M alternative weights for the same event are unchanged

#### Scenario: nomc combined-block correctness
- **WHEN** `--nomc --includeKaonScaleVariations` is used (kaon + μ1 + μ2 concatenated)
- **THEN** the kaon entry in the concatenated kinematic vectors uses the mass-corrected
  e-term uncertainty
- **AND** the muon entries continue to use the massless computation (particle_mass = 0)
