## Context
Three separate histogram outputs are currently produced for the nomc calibration
uncertainty, one per particle (kaon, mu1, mu2). This change merges them into one
combined histogram by feeding the existing `JpsiCorrectionsUncHelperSplines` helper
with concatenated per-particle kinematics so that its internal product loop covers
all three particles in a single call.

## Goals / Non-Goals
- Goals: one combined A/e/M NOI set per (eta bin, A/e/M label) for nomc mode;
  per-particle tensors available as validation output only.
- Non-Goals: changing anything outside the `--nomc --includeKaonScaleVariations`
  code path; changing tensor axes or NOI naming for the combined histogram;
  any interface changes to `muon_calibration.py`.

## Decisions

**Concatenated kinematics approach**

`JpsiCorrectionsUncHelperSplines::operator()` takes `RVec<float>` for recPts,
recEtas, etc. and iterates over `recPts.size()`, multiplying each particle's alt
weight into the running product. Feeding it 3-element RVecs — formed by
`ROOT::VecOps::Concatenate` of the already-defined per-particle columns — is
sufficient; the product falls out of the existing loop with no additional code.

The three response-weight columns are pre-computed from different TFLite models
(kaon: `diff_weights_helper`; muons: `muon_diff_weights_helper`) and are already
defined as 1-element `RVec<std::pair<double,double>>` columns before this logic runs.
Concatenating them gives a 3-element response-weight RVec accepted by the helper.

- Alternative considered: call `add_jpsi_crctn_stats_unc_hists` three times with
  `write_hist=False` and multiply the output tensors in the histmaker. Rejected
  because it introduces a new parameter, duplicates the loop, and obscures that the
  C++ helper already does the product.

**New combined prefix: `btojpsik_combined_sel`**

A new column prefix `btojpsik_combined_sel` is used for the concatenated columns.
The kaon-only path (`--nomc` absent) continues to call `add_jpsi_crctn_stats_unc_hists`
with `reco_sel_GF` unchanged.

**Validation per-particle histograms**

When `--validationHists` is active, the three per-particle HistoBoosts are written
by calling `add_jpsi_crctn_stats_unc_hists` separately for each particle with the
appropriate suffix. These calls are purely additive and do not affect the primary
combined histogram.

## Risks / Trade-offs
- Histfiles produced before this change contain `_mu1`/`_mu2` histograms; the updated
  tensor writer ignores them (extra HDF5 keys are harmless). New histfiles cannot be
  used with the old tensor writer in `--nomc` mode because `nominal_muonScaleSyst_responseWeights`
  was previously kaon-only.

## Open Questions
- None.
