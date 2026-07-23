## 1. btojpsik_axes.py — m32 axis definitions

- [x] 1.1 Add `"bkmm_nomc_m32"` axis: `Regular(100, 1.8, 2.6, name="bkmm_nomc_m32",
      underflow=False, overflow=False)` near the existing `bkmm_nomc_mass` entry
- [x] 1.2 Add `"bkmm_jpsimc_m32"` axis with the same binning near `bkmm_jpsimc_mass`

## 2. btojpsik.py — --histToFit argument

- [x] 2.1 Add `--histToFit` / `--hist-to-fit` argument with
      `choices=["mB", "m32"]`, `default="mB"`, and a short help string

## 3. btojpsik.py — m32 column Define in build_graph

- [x] 3.1 After the `bkmm_{vertex_prefix}_mass_scalar` Define, when
      `args.histToFit == "m32"`: add
      ```python
      df = df.Define(
          f"bkmm_{vertex_prefix}_m32_scalar",
          f"bkmm_{vertex_prefix}_mass_scalar - static_cast<double>(mm_kin_mass)",
      )
      ```
      Note: `mm_kin_mass` is already a scalar after `select_only_passing_bkmm_candidates`
      (all `mm_*` columns are reindexed to scalars via `[mm_best_idx]`; all `bkmm_*`
      columns become length-1 RVecs, hence `[0]` for the B mass scalar).

## 4. btojpsik.py — configure_fit_histogram

- [x] 4.1 Add `hist_to_fit: str = "mB"` parameter to `configure_fit_histogram`
- [x] 4.2 When `hist_to_fit == "m32"`:
      set `fit_mass_col = f"bkmm_{vprefix}_m32_scalar"` and
      `fit_mass_axis = all_butojpsik_axes[f"bkmm_{vprefix}_m32"]`
- [x] 4.3 Update the call site in `build_graph` to pass `args.histToFit`

## 5. Validation

- [x] 5.1 Syntax-check: `python -m py_compile scripts/histmakers/btojpsik.py
      wremnants/production/btojpsik_axes.py`
- [ ] 5.2 Confirm default `--histToFit mB` produces `nominal_HistToFit` with
      `bkmm_{vertex_prefix}_mass` axis (unchanged)
- [ ] 5.3 Smoke-test `--histToFit m32 --nomc` by verifying `nominal_HistToFit`
      observable axis is `bkmm_nomc_m32` and values are in [1.8, 2.6]
