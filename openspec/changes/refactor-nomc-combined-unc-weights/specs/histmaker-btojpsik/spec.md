## MODIFIED Requirements

### Requirement: Muon Scale Variations in nomc Mode
The histmaker SHALL compute per-muon scale variation response weights when both
`--nomc` and `--includeKaonScaleVariations` are active and the dataset has gen
kinematics. Each muon's weight SHALL use `SplinesDifferentialWeightsHelper` loaded
from `calib_filepaths["tflite_file"]` (`muon_response.tflite`) with reco values
passed as gen for all six kinematic inputs (pt, eta, charge), consistent with the
kaon convention. The kaon, mu1, and mu2 alt-weight tensors SHALL be multiplied
element-wise into a single combined tensor, and a single
`nominal_muonScaleSyst_responseWeights` histogram SHALL be written covering all three
particles. Per-particle histograms (`_kaon`, `_mu1`, `_mu2`) SHALL be written only
when `--validationHists` is active.

#### Scenario: combined weight tensor written
- **WHEN** `--nomc` and `--includeKaonScaleVariations` are both active
- **THEN** the output HDF5 contains exactly one `nominal_muonScaleSyst_responseWeights`
  histogram whose alt-weight tensor is the element-wise product of the kaon, mu1,
  and mu2 tensors

#### Scenario: per-particle histograms available as validation
- **WHEN** `--nomc`, `--includeKaonScaleVariations`, and `--validationHists` are all active
- **THEN** the output HDF5 additionally contains `nominal_muonScaleSyst_kaon_responseWeights`,
  `nominal_muonScaleSyst_mu1_responseWeights`, and `nominal_muonScaleSyst_mu2_responseWeights`
  alongside the combined histogram

#### Scenario: per-particle histograms absent without --validationHists
- **WHEN** `--nomc` and `--includeKaonScaleVariations` are both active but
  `--validationHists` is not set
- **THEN** no per-particle `_kaon`, `_mu1`, or `_mu2` scale variation histograms
  are written

#### Scenario: muon variations absent in standard mode
- **WHEN** `--nomc` is not set, even if `--includeKaonScaleVariations` is set
- **THEN** no combined or per-particle muon scale variation histograms are written;
  only the standard kaon `nominal_muonScaleSyst_responseWeights` is produced
