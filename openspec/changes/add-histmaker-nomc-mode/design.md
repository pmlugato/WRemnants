## Context

`btojpsik.py` hardcodes `bkmm_jpsimc_*` column names in ~6 places (fit observable,
kaon pT/eta inputs, `Define` for the mass scalar, the `configure_fit_histogram`
helper, and the call to `build_fit_pt_quantile_hists`). Four selection functions in
`btojpsik_selections.py` also hardcode those names. A fifth hardcoded reference lives
in the best-candidate picker C++ snippet inside
`select_only_passing_bkmm_candidates`.

Muon scale variations are currently only computed for the kaon via
`diff_weights_helper` (kaon_response.tflite). In nomc mode the J/ψ mass constraint is
absent, so the muon momenta enter the reconstructed B mass directly. The muon
response tflite is already registered in `calib_filepaths["tflite_file"]`.

## Goals / Non-Goals

- Goals:
  - `--nomc` produces the same histogram/tensor structure as the default mode, but
    sourced from `bkmm_nomc_*`.
  - Kaon-style scale variation histograms are produced for each J/ψ muon when both
    `--nomc` and `--includeKaonScaleVariations` are set.
  - No duplication of `build_graph` logic.
  - **Zero impact on existing behaviour**: running the histmaker without `--nomc`
    MUST produce bit-identical output to the pre-change code.
- Non-Goals:
  - Changing the selection thresholds or axis definitions.
  - Muon scale variations in the standard `bkmm_jpsimc_*` mode.
  - Extending tensor writer or fit scripts (those consume the same HDF5 schema).

## Backward-Compatibility Guarantee

Every interface change in this proposal uses a default value that reproduces the
current behaviour:

| Modified interface | New parameter | Default | Effect of default |
|-|-|-|-|
| `select_kaon_eta` / `select_kaon_pt` | `vertex_prefix` | `"jpsimc"` | unchanged column name |
| `select_bkmm_vtx_prob` / `select_bkmm_mass_window` | `vertex_prefix` | `"jpsimc"` | unchanged column name |
| `select_only_passing_bkmm_candidates` | `vertex_prefix` | `"jpsimc"` | unchanged C++ snippet |
| `add_jpsi_crctn_stats_unc_hists` | `hist_suffix` | `""` | unchanged histogram names |

The `muon_diff_weights_helper` object and all per-muon column `Define` calls are
placed inside `if args.nomc:` blocks and are never executed unless `--nomc` is
passed.

## Decisions

### 1. vertex_prefix string threading

Introduce `vertex_prefix = "nomc" if args.nomc else "jpsimc"` at the top of
`btojpsik.py`. Thread it into:
- `get_bkmm_selections()` lambdas that call the four parametrized selection functions
- `select_only_passing_bkmm_candidates(..., vertex_prefix=vertex_prefix)`
- The `Define("bkmm_{prefix}_mass_scalar", ...)` and alias lines in `build_graph`
- `configure_fit_histogram` (the mass column name)
- `build_fit_pt_quantile_hists` (the kaon pT/eta column names)

This avoids any code duplication and is easy to grep.

### 2. Parametrizing selection functions

Add `vertex_prefix: str = "jpsimc"` to `select_kaon_eta`, `select_kaon_pt`,
`select_bkmm_vtx_prob`, `select_bkmm_mass_window`, and
`select_only_passing_bkmm_candidates`. The C++ snippet in the best-candidate picker
that references `bkmm_jpsimc_vtx_prob` becomes an f-string substitution.

Alternative considered: duplicate the selections into a separate `_nomc` variant.
Rejected because it doubles maintenance surface with identical logic.

### 3. muon_diff_weights_helper instantiation

Create `muon_diff_weights_helper` at module load alongside `diff_weights_helper`,
guarded by the same condition (`smearingWeightsSplines` or `validationHists`). It
loads `calib_filepaths["tflite_file"]`. This matches the kaon pattern exactly.

### 4. Per-muon column extraction

After `select_only_passing_bkmm_candidates` the `mm_*` columns are reduced to
scalars. Wrap them as length-1 RVecs (e.g.
`ROOT::RVec<float>{static_cast<float>(mm_kin_mu1pt)}`) to match the interface
expected by `SplinesDifferentialWeightsHelper` and the uncertainty helper.

Charge is derived from `Muon_charge[mm_mu1_index]` (Muon array is not narrowed by
the candidate selection and remains a full vector).

Reco values are passed as gen for all six inputs (pt, eta, charge), consistent with
the kaon convention where gen eta and charge are already the reco values.

### 5. Naming collision resolution in add_jpsi_crctn_stats_unc_hists

Add `hist_suffix: str = ""` parameter. Histogram names `nominal_muonScaleSyst_*`
and `muonScaleSyst_*` become `nominal_muonScaleSyst{hist_suffix}_*`. Existing callers
pass no argument; the default preserves the current names. The two per-muon calls
use `hist_suffix="_mu1"` and `hist_suffix="_mu2"`.

## Risks / Trade-offs

- The C++ f-string substitution for the best-candidate picker is a minor readability
  concern, but the snippet is short and the pattern is already used elsewhere in the
  file.
- Two extra `add_jpsi_crctn_stats_unc_hists` calls per event loop add some overhead,
  but this is gated on `--nomc` and `--includeKaonScaleVariations` together.

## Open Questions

- None at this stage; clarifications confirmed: (a) reco values used as gen for
  muons, (b) J/ψ correction still applied on top of `bkmm_nomc_kaon1pt`.
