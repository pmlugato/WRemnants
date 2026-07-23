# Change: Add --nomc mode to btojpsik tensor writer

## Why

The btojpsik tensor writer hardcodes `bkmm_jpsimc_mass` as the fit observable and
loads only the kaon response-weight tensor (`nominal_muonScaleSyst_responseWeights`).
When histograms are produced with `--nomc`, the fit observable is `bkmm_nomc_mass` and
the histmaker additionally writes per-muon response tensors
(`nominal_muonScaleSyst_mu1_responseWeights`, `_mu2_*`). The tensor writer must be
extended to consume those tensors and emit A/e/M nuisances for each J/ψ muon.

## What Changes

- Add `--nomc` flag to `btojpsik_tensor.py`; changes the `--massAxis` default from
  `bkmm_jpsimc_mass` to `bkmm_nomc_mass`
- When `--nomc`: load `nominal_muonScaleSyst_mu1_responseWeights` and
  `nominal_muonScaleSyst_mu2_responseWeights` from the HDF5 file
- Add a `load_variation_hist(filename, dataset, hist_key)` helper to
  `rabbit_btojpsik_helpers.py` for loading individual variation hists
- Run the A/e/M NOI loop for mu1 and mu2, naming systematics
  `A_mu1_eta{i}`, `e_mu1_eta{i}`, `M_mu1_eta{i}` (and `_mu2_*`)
- Kaon A/e/M nuisances (`A_eta{i}`, …) are always emitted regardless of `--nomc`
- Apply the same `etaBins` rebinning to mu1/mu2 variation hists as to the kaon hist

## Backward Compatibility

Every change uses a default that reproduces existing behaviour:

| Modified interface | New parameter | Default | Effect of default |
|-|-|-|-|
| `btojpsik_tensor.py` | `--nomc` | `False` | `--massAxis` stays `bkmm_jpsimc_mass`; no mu variation hists loaded |

The `load_variation_hist` helper is additive; `load_histogram` is unchanged.

## Impact

- Affected specs: tensor-btojpsik (new)
- Affected code: `scripts/rabbit/btojpsik_tensor.py`,
  `wremnants/postprocessing/rabbit_btojpsik_helpers.py`
