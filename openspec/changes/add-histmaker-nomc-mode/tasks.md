## 1. btojpsik_selections.py — parametrize vertex_prefix

- [x] 1.1 Add `vertex_prefix: str = "jpsimc"` parameter to `select_kaon_eta` and
      `select_kaon_pt`; replace the hardcoded `bkmm_jpsimc_kaon1eta` /
      `bkmm_jpsimc_kaon1pt` strings with `f"bkmm_{vertex_prefix}_kaon1eta"` /
      `f"bkmm_{vertex_prefix}_kaon1pt"`
- [x] 1.2 Add `vertex_prefix: str = "jpsimc"` to `select_bkmm_vtx_prob` and
      `select_bkmm_mass_window`; replace `bkmm_jpsimc_vtx_prob` /
      `bkmm_jpsimc_mass` with f-string equivalents
- [x] 1.3 Add `vertex_prefix: str = "jpsimc"` to
      `select_only_passing_bkmm_candidates`; replace the hardcoded
      `bkmm_jpsimc_vtx_prob` inside the C++ best-candidate picker f-string

## 2. btojpsik.py — --nomc flag and variable substitution

- [x] 2.1 Add `--nomc` / `--no-mc` boolean flag to the argparser with a short
      help string
- [x] 2.2 Derive `vertex_prefix = "nomc" if args.nomc else "jpsimc"` immediately
      after argument parsing; thread into all relevant call sites
- [x] 2.3 Update `get_bkmm_selections()` to accept and forward `vertex_prefix` in
      the four lambda closures that call selection functions
- [x] 2.4 Replace hardcoded `bkmm_jpsimc_mass_scalar`, `bkmm_jpsimc_kaon1pt`,
      `bkmm_jpsimc_kaon1eta` column references in `build_graph` and
      `build_fit_pt_quantile_hists` with f-string equivalents using `vertex_prefix`
- [x] 2.5 Update `configure_fit_histogram` to accept `vertex_prefix` and use it
      for the fit mass column name

## 3. muon_calibration.py — hist_suffix parameter

- [x] 3.1 Add `hist_suffix: str = ""` parameter to
      `add_jpsi_crctn_stats_unc_hists`; append it to all hardcoded histogram name
      strings inside the function (e.g. `"nominal_muonScaleSyst_responseWeights"`
      → `f"nominal_muonScaleSyst{hist_suffix}_responseWeights"`)

## 4. btojpsik.py — muon scale variations in --nomc mode

- [x] 4.1 Instantiate `muon_diff_weights_helper` at module load from
      `calib_filepaths["tflite_file"]`, guarded by the same
      `smearingWeightsSplines or validationHists` condition used for
      `diff_weights_helper`; set to `None` otherwise
- [x] 4.2 In `build_graph`, after candidate selection, when `args.nomc` and
      `has_gen_kinematics`: define length-1 RVec columns for each muon —
      `mm_mu1_stuff_recoPt` / `recoEta` / `recoCharge` from `mm_kin_mu1pt`,
      `mm_kin_mu1eta`, `Muon_charge[mm_mu1_index]`, and equivalents for mu2
- [x] 4.3 Compute per-muon response weights:
      `mm_mu1_stuff_response_weight` and `mm_mu2_stuff_response_weight` by calling
      `muon_diff_weights_helper` with reco values passed as gen for all six inputs
- [x] 4.4 When `args.nomc` and `args.includeKaonScaleVariations` and
      `has_gen_kinematics`: call `add_jpsi_crctn_stats_unc_hists` twice using
      `reco_sel_GF="mm_mu1_stuff"` with `hist_suffix="_mu1"` and
      `reco_sel_GF="mm_mu2_stuff"` with `hist_suffix="_mu2"`

## 5. Validation

- [x] 5.1 Syntax-check all modified files:
      `python -m py_compile scripts/histmakers/btojpsik.py
      wremnants/production/btojpsik_selections.py
      wremnants/production/muon_calibration.py`
- [ ] 5.2 Confirm default mode (no `--nomc`) is unaffected by running the
      existing histmaker invocation and verifying `nominal_HistToFit` uses
      `bkmm_jpsimc_mass`
- [ ] 5.3 Smoke-test `--nomc` by verifying `nominal_HistToFit` observable column
      is `bkmm_nomc_mass`
