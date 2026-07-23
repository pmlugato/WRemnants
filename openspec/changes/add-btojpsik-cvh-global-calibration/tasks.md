> **PAUSED 2026-07-22** — blocked on `add-alcareco-to-nanoaod`, which deletes
> `JpsiKCandidateSplitter` and collapses the two makers into one job. Refactor
> §4–§7 onto that path before resuming (one-job grads ⇒ same events by
> construction ⇒ §5 file-matching is deleted). See design.md STATUS.

## 1. Fork and set up the solve framework

- [x] 1.1 Fork `MuonMomentumScaleCalibration` (`aggregategrads.py`, `aggregategrads.cpp`, `globalfith5pynoreduction.py`, `data/`) from `/work/submit/david_w/ZMass/MuonMomentumScaleCalibration/` to `/work/submit/pmlugato/REALmz/mz/real/calibration/`
- [x] 1.2 Confirm the fork runs end-to-end on a small existing input (smoke: aggregate a few `globalcor_*.root`, solve, write `correctionResults`)
- [x] 1.3 Add thin repo-versioned wrappers under `scripts/btojpsik/calibration/` that invoke the fork with study-specific config paths

## 2. Enable calibration output on the B⁺→J/ψK⁺ makers

- [x] 2.1 Add `fillGrads`/`fillGradsFactored`/`fillRunTree` VarParsing knobs to `runCvhBplusJpsiK.py`, forwarded to both maker clones (default OFF)
- [x] 2.2 Verify defaults leave `ResidualGlobalCorrectionMakerTwoTrackJpsiKMuMuG4e_cfi.py` and `ResidualGlobalCorrectionMakerJpsiKSingleTrackKaonG4e_cfi.py` bit-identical for the diagnostic-join path
- [x] 2.3 Configure both makers with the grouped material model (`globalMaterialModel=True` + shared `materialGroupsFile` + `skipHitlessSurfaces=True`) so material is parmtype 15; field stays parmtype 14
- [x] 2.4 200-event smoke with grads ON (both `mode=dimuon` and `mode=kaon`); confirm `nParms`, `globalidxv`, `gradv`, `hessfactorv` branches present and the `runtree` contains parmtype 14 (~50) and parmtype 15 (~50), no parmtype 7

## 3. Parameter-catalog consistency check

- [x] 3.1 Write a check that loads the `runtree` from the J/ψ→μμ, B dimuon, and B kaon productions and asserts identical length + ordered `(parmtype, subdet, layer, rawdetid)` (all built with the same `materialGroupsFile` + field-mode set)
- [x] 3.2 Run the check across the three productions; treat any mismatch as a hard stop and document the shared build/config that guarantees a match

## 4. B⁺→J/ψK⁺ CVH grads production (Run2016H, subset OK)

- [ ] 4.1 Run `runCvhBplusJpsiK.py mode=dimuon` (grads ON, `fillRunTree=True` on the first job) on `TkAlJpsiX_*.root`
- [ ] 4.2 Run `runCvhBplusJpsiK.py mode=kaon` (grads ON, `plimit=0.05`) on the same files (separate jobs — single G4 master per process)
- [ ] 4.3 Report kaon survival fraction (plimit/maxlen losses, convergence)

## 5. Aggregation into the shared parameter space

- [ ] 5.1 Aggregate the B dimuon grads with dimuon quality filters → `combinedgrads_Bdimuon.hdf5`
- [x] 5.2 Add a kaon-aggregation variant of `aggregategrads.py` with kaon-appropriate filters (nvalid, pixel, χ²/ndof, edm; no `Jpsi_mass`) scattering into the same `nparms`
- [ ] 5.3 Produce `combinedgrads_kaon.hdf5`; log accepted candidate count

## 6. Primary fits — within-B dimuon-only vs dimuon+kaon

- [ ] 6.1 Solve the dimuon-only fit from `combinedgrads_Bdimuon.hdf5` with `CVH_FIX_ALIGNMENT=1` (nominal; only the 92 globals float) → `correctionResults_Bdimuon.root`
- [ ] 6.2 Solve the dimuon+kaon fit from `combinedgrads_Bdimuon.hdf5` + `combinedgrads_kaon.hdf5` (identical config to 6.1) → `correctionResults_Bdimuon_kaon.root`
- [ ] 6.3 Repeat both with alignment marginalised (`CVH_FIX_ALIGNMENT=0`) as the variant for absolute-uncertainty statements
- [ ] 6.4 Confirm covariance extraction (diagonal + parmtype-14 x 15 sub-covariance) in all four solves

## 7. Sanity check — J/ψ consistency

- [ ] 7.1 Run the standalone J/ψ→μμ CVH grads on `TkAlJpsiMuMu_*.root` (Run2016H, subset OK) with the same grouped-material config, aggregate → `combinedgrads_jpsimumu.hdf5`, solve → `correctionResults_jpsimumu.root`
- [ ] 7.2 Compare b-field/material corrections: `correctionResults_jpsimumu` vs `correctionResults_Bdimuon`; confirm agreement within uncertainties (flag any discrepancy)

## 7b. Species closure check (muon mass hypothesis)

- [ ] 7b.1 Rerun the bachelor-kaon grads with `kaonAsMuon=True` on the same input files (identical everything else)
- [ ] 7b.2 Aggregate → `combinedgrads_kaonAsMuon.hdf5`; solve dimuon+kaonAsMuon with the nominal config
- [ ] 7b.3 Compare the material group scalings (parmtype 15) between the kaon and muon hypotheses; attribute any shift to the energy-loss model

## 8. Reporting

- [ ] 8.1 Extract per-parameter constraint σ and magnitudes for dimuon-only and dimuon+kaon, splitting field mode 0 (scale) from modes 1+ (shape) and the material groups
- [ ] 8.2 Extract the b-field↔material correlation from the block sub-covariance for both fits (diagnostic; expected small)
- [ ] 8.3 Produce a side-by-side table/plot (dimuon-only vs dimuon+kaon) of magnitude, uncertainty, and correlation under `scripts/btojpsik/calibration/`; note which blocks are prior-dominated
- [ ] 8.4 Write a short summary: what the kaon actually constrains (shape? material?), how the species closure check came out, and what that implies for the channel

## 9. Slides

- [x] 9.1 Produce a detailed slide deck covering: motivation and the b-field↔material degeneracy; the three-stage pipeline; the field-mode (14) + grouped-material (15) parameterization; the within-B comparison design and why (vs full-inclusive, vs prompt-selection); the iterations/decisions taken along the way
- [ ] 9.2 Include results: kaon survival fraction, the J/ψ-consistency sanity check, and the dimuon-only vs dimuon+kaon comparison of b-field/material magnitude, uncertainty, and correlation (plots/tables from §8)

## 10. Validation

- [x] 10.1 `openspec validate add-btojpsik-cvh-global-calibration --strict --no-interactive` passes
- [ ] 10.2 Dimuon-only and dimuon+kaon fits use identical priors/NMR (only the kaon grads differ), confirmed from solve config
