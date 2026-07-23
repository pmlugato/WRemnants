## Context

`btojpsik_tensor.py` currently calls `load_histogram(infile, datasetSignal)` which
returns `(signal_hist, data_hist, variation_hist)` where `variation_hist` is always
`nominal_muonScaleSyst_responseWeights` (kaon). In nomc mode the histmaker writes two
additional response tensors: `nominal_muonScaleSyst_mu1_responseWeights` and
`nominal_muonScaleSyst_mu2_responseWeights`. The tensor writer must load and process
these alongside the kaon tensor.

## Goals / Non-Goals

- Goals:
  - `--nomc` produces the same tensor structure as default mode (kaon A/e/M NOIs
    present) plus additional per-muon A/e/M NOIs named `{label}_mu1_eta{i}` /
    `{label}_mu2_eta{i}`.
  - `--nomc` changes the `--massAxis` default to `bkmm_nomc_mass`.
  - **Zero impact on existing behaviour**: running without `--nomc` MUST produce
    bit-identical output to the pre-change code.
- Non-Goals:
  - Changing the kaon A/e/M nuisance names.
  - Changing `load_histogram`'s signature or return value.
  - Adding muon diagnostic plots (curvature response, variation projection).

## Backward-Compatibility Guarantee

| Modified interface | New parameter | Default | Effect of default |
|-|-|-|-|
| `btojpsik_tensor.py` | `--nomc` | `False` | unchanged `--massAxis`; no mu variation hists loaded; no mu NOIs added |

All nomc code paths are gated on `args.nomc`.

## Decisions

### 1. `load_variation_hist` helper

Add `load_variation_hist(filename, dataset, hist_key)` to
`rabbit_btojpsik_helpers.py`. It re-opens the HDF5 file (like `load_histogram` does)
and returns a single variation hist by key. This keeps `load_histogram` unchanged and
avoids baking a variable-length return into an existing function.

Alternative considered: extend `load_histogram` to accept a list of extra hist keys.
Rejected because it changes the public return signature and forces callers to unpack a
variable-length tuple.

### 2. Per-muon A/e/M loop

After the existing kaon A/e/M loop, when `args.nomc`, run the same loop twice more
— once with `variation_hist_mu1` and suffix `_mu1`, once with `variation_hist_mu2`
and suffix `_mu2`. The loop body is factored into a helper
`_add_aem_systematics(writer, variation_hist, signal_hist, args, name_suffix="")` to
avoid copy-paste. The kaon call passes `name_suffix=""` (produces `A_eta{i}` etc.),
muon calls pass `name_suffix="_mu1"` / `"_mu2"`.

### 3. Rebinning

Apply the same `rebin_variation_unc_axis` + `rebin_histogram` steps to
`variation_hist_mu1` and `variation_hist_mu2` as currently applied to the kaon
`variation_hist`, guarded by `args.nomc`.

### 4. massAxis default in --nomc

`--nomc` sets `args.massAxis = args.massAxis or "bkmm_nomc_mass"` — but more
precisely, the argument default is changed to `None` and the tensor writer resolves it:
`"bkmm_nomc_mass" if args.nomc else "bkmm_jpsimc_mass"` when `args.massAxis` is None.
This keeps explicit `--massAxis` overrides working in both modes.

## Risks / Trade-offs

- Two extra passes through the A/e/M loop add a small amount of tensor-writing
  overhead, gated on `--nomc`.
- Re-opening the HDF5 file twice for mu hists is a minor inefficiency; acceptable
  given the file is small and this is a one-shot script.

## Open Questions

- None.
