## 1. rabbit_btojpsik_helpers.py — load_variation_hist helper

- [x] 1.1 Add `load_variation_hist(filename, dataset, hist_key)` that opens the HDF5
      file, resolves the dataset key via `_resolve_dataset_key`, and returns
      `results[dataset_key]["output"][hist_key].get()`

## 2. btojpsik_tensor.py — --nomc flag and massAxis default

- [x] 2.1 Add `--nomc` boolean flag to `parse_args`
- [x] 2.2 Change `--massAxis` default to `None`; resolve after parsing:
      `args.massAxis = args.massAxis or ("bkmm_nomc_mass" if args.nomc else "bkmm_jpsimc_mass")`

## 3. btojpsik_tensor.py — load mu variation hists

- [x] 3.1 After `load_histogram(...)`, when `args.nomc`: call `load_variation_hist`
      for `nominal_muonScaleSyst_mu1_responseWeights` and
      `nominal_muonScaleSyst_mu2_responseWeights`
- [x] 3.2 Apply the same `rebin_variation_unc_axis` + `rebin_histogram` rebinning to
      the mu1/mu2 variation hists as to the kaon `variation_hist`

## 4. btojpsik_tensor.py — factor A/e/M loop into helper

- [x] 4.1 Extract the existing A/e/M NOI loop (lines ~688–749) into a module-level
      helper `_add_aem_systematics(writer, variation_hist, signal_hist, args, name_suffix="")`
      that names each systematic `f"{label}{name_suffix}_eta{eta_idx}"`
- [x] 4.2 Replace the existing call site with `_add_aem_systematics(..., name_suffix="")`
- [x] 4.3 When `args.nomc`: call `_add_aem_systematics` with `variation_hist_mu1`
      and `name_suffix="_mu1"`, then with `variation_hist_mu2` and `name_suffix="_mu2"`

## 5. Validation

- [x] 5.1 Syntax-check modified files:
      `python -m py_compile scripts/rabbit/btojpsik_tensor.py
      wremnants/postprocessing/rabbit_btojpsik_helpers.py`
- [ ] 5.2 Confirm default mode (no `--nomc`) produces the same nuisance names
      (`A_eta0`, `e_eta0`, `M_eta0`, …) as before
- [ ] 5.3 Smoke-test `--nomc` by verifying nuisance names include
      `A_mu1_eta0`, `e_mu1_eta0`, `M_mu1_eta0` alongside the kaon set
