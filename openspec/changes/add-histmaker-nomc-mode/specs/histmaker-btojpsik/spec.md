## ADDED Requirements

### Requirement: nomc Vertex Mode
The histmaker SHALL accept a `--nomc` flag that switches all vertex-fit column
references from `bkmm_jpsimc_*` to `bkmm_nomc_*`, covering the fit observable,
kaon kinematic inputs to the J/ψ momentum correction, selection thresholds, and the
best-candidate picker. The default behaviour (no flag) MUST be unchanged.

#### Scenario: nomc fit observable
- **WHEN** `btojpsik.py` is invoked with `--nomc`
- **THEN** the `nominal_HistToFit` histogram uses `bkmm_nomc_mass` as its observable
  column instead of `bkmm_jpsimc_mass`

#### Scenario: nomc kaon kinematic inputs
- **WHEN** `--nomc` is active
- **THEN** the J/ψ momentum correction helper receives `bkmm_nomc_kaon1pt` and
  `bkmm_nomc_kaon1eta` as inputs instead of `bkmm_jpsimc_kaon1pt` /
  `bkmm_jpsimc_kaon1eta`

#### Scenario: nomc selection functions
- **WHEN** `--nomc` is active
- **THEN** `select_kaon_eta`, `select_kaon_pt`, `select_bkmm_vtx_prob`,
  `select_bkmm_mass_window`, and the best-candidate picker all evaluate
  `bkmm_nomc_*` columns instead of `bkmm_jpsimc_*`

#### Scenario: default mode unchanged
- **WHEN** `--nomc` is not passed
- **THEN** all column references remain `bkmm_jpsimc_*` and output histograms are
  identical to the pre-change behaviour

### Requirement: Muon Scale Variations in nomc Mode
The histmaker SHALL compute per-muon scale variation response weights when both
`--nomc` and `--includeKaonScaleVariations` are active and the dataset has gen
kinematics. Each muon's weight SHALL use `SplinesDifferentialWeightsHelper` loaded
from `calib_filepaths["tflite_file"]` (`muon_response.tflite`) with reco values
passed as gen for all six kinematic inputs (pt, eta, charge), consistent with the
kaon convention. Two sets of scale variation histograms SHALL be written with
distinct suffixes (`_mu1`, `_mu2`) to avoid naming collisions.

#### Scenario: per-muon response weights computed
- **WHEN** `--nomc` and `--includeKaonScaleVariations` are both active and the
  dataset has gen kinematics
- **THEN** `mm_mu1_stuff_response_weight` and `mm_mu2_stuff_response_weight` columns
  are defined using `muon_diff_weights_helper`, sourcing kinematics from
  `mm_kin_mu1pt`, `mm_kin_mu1eta`, `Muon_charge[mm_mu1_index]` (and mu2 equivalents)
  with reco values passed as gen

#### Scenario: muon scale variation histograms written
- **WHEN** `--nomc` and `--includeKaonScaleVariations` are both active
- **THEN** the output HDF5 contains `nominal_muonScaleSyst_mu1_responseWeights` and
  `nominal_muonScaleSyst_mu2_responseWeights` histograms in addition to the existing
  kaon scale variation histograms

#### Scenario: muon variations absent in standard mode
- **WHEN** `--nomc` is not set, even if `--includeKaonScaleVariations` is set
- **THEN** no muon-specific scale variation histograms are written
