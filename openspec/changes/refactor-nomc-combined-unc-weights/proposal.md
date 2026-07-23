# Change: Combine kaon + muon scale-variation weights into one tensor in nomc mode

## Why
When `--nomc` and `--includeKaonScaleVariations` are both active, the histmaker
currently produces three independent uncertainty histograms (kaon, mu1, mu2) that
feed three independent sets of A/e/M NOIs in the fit tensor. The physically correct
treatment mirrors `JpsiCorrectionsUncHelperSplines` in `muon_calibration.hpp`:
the combined observable's response to a single calibration shift is the element-wise
product of the individual particle responses. Three independent NOIs overcount the
degrees of freedom and may introduce degeneracies.

## What Changes
- When `--nomc` and `--includeKaonScaleVariations` are both active, the histmaker
  SHALL define a set of concatenated kinematic columns (kaon + mu1 + mu2) under a
  single combined prefix and call `add_jpsi_crctn_stats_unc_hists` **once**. Because
  `JpsiCorrectionsUncHelperSplines` already loops over the input RVec length and
  multiplies alt weights in-place, feeding it 3-element RVecs naturally produces the
  combined product tensor. This applies regardless of `--histToFit`.
- The result is a single `nominal_muonScaleSyst_responseWeights` histogram covering
  all three particles, identical in name and tensor axes to the existing kaon-only
  histogram.
- The separate per-call approach (calling `add_jpsi_crctn_stats_unc_hists` once per
  particle) is removed; `muon_calibration.py` requires no interface changes.
- The separate `_mu1` / `_mu2` histograms move to optional `--validationHists`
  output (separate per-particle calls gated on that flag).
- **BREAKING** (nomc mode only): `nominal_muonScaleSyst_mu1_responseWeights` and
  `nominal_muonScaleSyst_mu2_responseWeights` are no longer written to the primary
  HDF5 output. The tensor writer must not attempt to load them.
- `btojpsik_tensor.py`: remove the `variation_hist_mu1` / `variation_hist_mu2`
  loading block and the two extra `_add_aem_systematics` calls; the single combined
  histogram drives all A/e/M NOIs.

## Backward Compatibility
- Without `--nomc`, the kaon-only path is unchanged; `add_jpsi_crctn_stats_unc_hists`
  is called once with `reco_sel_GF` as before.
- The combined histogram has the same key name and tensor axes as the old kaon-only
  histogram, so the tensor writer needs only the removal of the mu1/mu2 loading block.
- Without `--nomc`, the tensor writer output is bit-identical to pre-change behaviour.

## Impact
- Affected specs: `histmaker-btojpsik`, `tensor-btojpsik`
- Affected code:
  - `scripts/histmakers/btojpsik.py` — concatenated kinematics + single combined call
  - `scripts/rabbit/btojpsik_tensor.py` — remove mu1/mu2 loading + NOI calls
  - `wremnants/production/muon_calibration.py` — no interface changes required
