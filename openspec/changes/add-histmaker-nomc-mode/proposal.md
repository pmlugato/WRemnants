# Change: Add --nomc mode to btojpsik histmaker with muon scale variations

## Why
The `bkmm_nomc_*` vertex fit omits the J/ψ MC mass constraint, so the reconstructed B
mass is directly sensitive to muon momentum scale. Running the histmaker in this mode
requires (a) switching all vertex-fit column references to `bkmm_nomc_*` and (b)
additionally propagating kaon-style scale variations to the two J/ψ muons, which is
unnecessary in the standard `bkmm_jpsimc_*` mode where the constraint absorbs muon
scale errors to first order.

## What Changes
- Add `--nomc` boolean CLI flag to `scripts/histmakers/btojpsik.py`; when set, a
  `vertex_prefix` variable switches from `"jpsimc"` to `"nomc"` and is threaded
  through all column-name references and selection-function calls.
- Parametrize the four selection functions in `btojpsik_selections.py` that currently
  hardcode `bkmm_jpsimc_*`: `select_kaon_eta`, `select_kaon_pt`,
  `select_bkmm_vtx_prob`, `select_bkmm_mass_window`. Also parametrize the
  best-candidate picker inside `select_only_passing_bkmm_candidates`, which hardcodes
  `bkmm_jpsimc_vtx_prob`.
- When `--nomc` and `--includeKaonScaleVariations` are both active, additionally
  compute scale variation response weights for the two J/ψ muons using
  `calib_filepaths["tflite_file"]` (`muon_response.tflite`). Muon kinematics are
  accessed from `mm_kin_mu1pt`/`mm_kin_mu1eta` and `mm_kin_mu2pt`/`mm_kin_mu2eta`;
  charges via `Muon_charge[mm_mu1_index]` and `Muon_charge[mm_mu2_index]`. Reco
  values are passed as gen (same convention already used for kaon eta/charge).
- Add a `hist_suffix` parameter to `add_jpsi_crctn_stats_unc_hists` (default `""`) so
  the two per-muon calls produce distinctly named output histograms (`_mu1`, `_mu2`)
  without colliding with the kaon histograms or each other.

## Backward Compatibility
All changes are strictly additive and opt-in. No existing call sites are broken:
- `--nomc` is a new flag; omitting it leaves behaviour identical to today.
- The `vertex_prefix` parameter added to selection functions defaults to `"jpsimc"`,
  so every existing caller continues to work without modification.
- The `hist_suffix` parameter added to `add_jpsi_crctn_stats_unc_hists` defaults to
  `""`, preserving the existing kaon histogram names exactly.
- The muon `diff_weights_helper` and per-muon scale variation code paths are gated
  behind `args.nomc` checks and are never reached in the default mode.

## Impact
- Affected specs: `histmaker-btojpsik` (new capability)
- Affected code:
  - `scripts/histmakers/btojpsik.py`
  - `wremnants/production/btojpsik_selections.py`
  - `wremnants/production/muon_calibration.py` (add `hist_suffix` param)
