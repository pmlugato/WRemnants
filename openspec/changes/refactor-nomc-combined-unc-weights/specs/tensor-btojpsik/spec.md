## MODIFIED Requirements

### Requirement: nomc Tensor Mode
The tensor writer SHALL accept a `--nomc` flag that switches the default fit-observable
axis from `bkmm_jpsimc_mass` to `bkmm_nomc_mass`. When `--nomc` is active, A/e/M
nuisances are emitted from the single combined `nominal_muonScaleSyst_responseWeights`
histogram only; the separate per-muon histograms
(`nominal_muonScaleSyst_mu1_responseWeights`, `nominal_muonScaleSyst_mu2_responseWeights`)
are NOT loaded and no per-muon nuisances are written.

#### Scenario: massAxis default changes with --nomc
- **WHEN** `btojpsik_tensor.py` is invoked with `--nomc` and no explicit `--massAxis`
- **THEN** the tensor is built using `bkmm_nomc_mass` as the fit-observable axis

#### Scenario: explicit --massAxis overrides --nomc default
- **WHEN** `btojpsik_tensor.py` is invoked with `--nomc --massAxis bkmm_custom_mass`
- **THEN** the tensor is built using `bkmm_custom_mass` as the fit-observable axis

#### Scenario: combined A/e/M nuisances in nomc mode
- **WHEN** `--nomc` is active
- **THEN** the tensor contains exactly one set of A/e/M nuisances (`A_eta{i}`,
  `e_eta{i}`, `M_eta{i}`) sourced from `nominal_muonScaleSyst_responseWeights`,
  whose alt-weight tensor already encodes the kaon, mu1, and mu2 combined response

#### Scenario: no per-muon nuisances in nomc mode
- **WHEN** `--nomc` is active
- **THEN** no `_mu1_*` or `_mu2_*` nuisances are written to the tensor

#### Scenario: default mode unaffected
- **WHEN** `--nomc` is not passed
- **THEN** the tensor output is bit-identical to the pre-change behaviour: only kaon
  A/e/M nuisances are present and `bkmm_jpsimc_mass` is the fit-observable axis
