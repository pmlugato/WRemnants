## ADDED Requirements

### Requirement: nomc Tensor Mode
The tensor writer SHALL accept a `--nomc` flag that switches the default fit-observable
axis from `bkmm_jpsimc_mass` to `bkmm_nomc_mass` and adds A/e/M nuisances for both
J/ψ muons in addition to the kaon nuisances.

#### Scenario: massAxis default changes with --nomc
- **WHEN** `btojpsik_tensor.py` is invoked with `--nomc` and no explicit `--massAxis`
- **THEN** the tensor is built using `bkmm_nomc_mass` as the fit-observable axis

#### Scenario: explicit --massAxis overrides --nomc default
- **WHEN** `btojpsik_tensor.py` is invoked with `--nomc --massAxis bkmm_custom_mass`
- **THEN** the tensor is built using `bkmm_custom_mass` as the fit-observable axis

#### Scenario: kaon A/e/M nuisances always present
- **WHEN** `--nomc` is active
- **THEN** kaon nuisances `A_eta{i}`, `e_eta{i}`, `M_eta{i}` are still written to
  the tensor unchanged

#### Scenario: mu1 and mu2 A/e/M nuisances added in nomc mode
- **WHEN** `--nomc` is active
- **THEN** nuisances `A_mu1_eta{i}`, `e_mu1_eta{i}`, `M_mu1_eta{i}` and
  `A_mu2_eta{i}`, `e_mu2_eta{i}`, `M_mu2_eta{i}` are written for every eta bin `i`

#### Scenario: mu variation hists sourced from HDF5
- **WHEN** `--nomc` is active
- **THEN** `nominal_muonScaleSyst_mu1_responseWeights` and
  `nominal_muonScaleSyst_mu2_responseWeights` are loaded from the signal dataset in the
  input HDF5 file

#### Scenario: default mode unaffected
- **WHEN** `--nomc` is not passed
- **THEN** the tensor output is bit-identical to the pre-change behaviour: only kaon
  A/e/M nuisances are present and `bkmm_jpsimc_mass` is the fit-observable axis
