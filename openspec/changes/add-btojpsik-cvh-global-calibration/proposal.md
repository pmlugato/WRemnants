# Change: Run the CVH global-correction calibration on B⁺→J/ψK⁺ and measure its effect on b-field / material corrections

## Why

The muon momentum-scale calibration is derived by re-fitting tracks with the
Geant4e CVH refit (`ResidualGlobalCorrectionMaker*`), accumulating each
candidate's gradient and Hessian in a shared global-parameter space, and
solving the resulting least-squares system for the correction parameters
(alignment, b-field, material energy loss, resolution, plus the global
scalar-potential field modes and per-group material scales). Historically this
was done on the standalone J/ψ→μμ (`TkAlJpsiMuMu`) sample.

The B⁺→J/ψK⁺ channel adds a **bachelor kaon single-track** probe: a soft,
displaced K± carrying a **kaon mass hypothesis** (`trackParticleName='kaon'`,
0.4937 GeV) that sets its dE/dx and Geant4 propagation. It samples the field
and material differently from the two muons, and — being a different species —
gives a handle on the **energy-loss model** that pions and muons do not.

The relevant channel roles are established (David Walter, 2026-07-21): J/ψ→μμ
sets the absolute momentum **scale** (stiff, high-pT muons, where energy loss
vanishes and the field scale is cleanly isolated); soft displaced hadrons set
the field **shape** and the **material**, because they curl through a wide
range of (r, φ, z) and respond strongly to the local field. The bachelor kaon
is soft (median pT 0.53 GeV) and displaced, so **shape and material — not
scale — are the plausible targets**.

This is an **exploratory study**: the question is what this channel actually
contributes on its own, not a claim that it resolves the calibration. K⁰s is
deliberately **excluded** so the B⁺→J/ψK⁺ contribution is seen in isolation,
and the fit is run **independently** rather than joined to the existing
combined fit. The AlCaReco exists (`TkAlJpsiX`), the Stage-2 CVH refit driver
(`runCvhBplusJpsiK.py`) is in place, and the solve framework is understood;
this change wires them together and measures the effect.

Not a production correction release.

## What Changes

- **Enable grad/Hess output** on the two B⁺→J/ψK⁺ CVH makers, which are
  currently configured for the diagnostic offline-join only
  (`fillGrads=False`, `fillRunTree=False`):
  `ResidualGlobalCorrectionMakerTwoTrackJpsiKMuMuG4e_cfi.py` (dimuon,
  mass-constraint ON) and `ResidualGlobalCorrectionMakerJpsiKSingleTrackKaonG4e_cfi.py`
  (bachelor kaon). Wire `fillGrads`/`fillRunTree`/`fillGradsFactored` knobs
  through `runCvhBplusJpsiK.py`.
- **Use the grouped global material model** (`globalMaterialModel=True` +
  `materialGroupsFile`, and its required `skipHitlessSurfaces=True`) on all
  three makers (B dimuon, B kaon, and the standalone J/ψ→μμ used for the sanity
  check), so material appears as ~50 grouped scales (parmtype 15) rather than
  the per-module energy-loss form (parmtype 7). The b-field is always the ~50
  global scalar-potential modes (parmtype 14); there is no per-detid b-field
  (parmtype 6) in this build. The three makers MUST share the identical
  `materialGroupsFile` and field-mode set so their parameter catalogs align.
- **Fork the solve framework** (`aggregategrads.py`, `aggregategrads.cpp`,
  `globalfith5pynoreduction.py` and support data) from David W's
  `MuonMomentumScaleCalibration` package into a user-owned fork under
  `/work/submit/pmlugato/REALmz/mz/real/calibration/`.
- **Run the B⁺→J/ψK⁺ CVH grads production** on `TkAlJpsiX_*.root` (Run2016H)
  in both `mode=dimuon` (jpsicons mass constraint) and `mode=kaon` (single
  track) as two separate `cmsRun` jobs per file.
- **Add a kaon single-track aggregation** variant with kaon-appropriate
  quality filters (the existing aggregation hard-codes dimuon branches such as
  `Jpsi_mass`, `Muplus_nvalid`), scattering into the **same** `nparms` global
  index space as the dimuon side.
- **Verify parameter-catalog consistency**: the dimuon and kaon makers MUST
  emit `globalidxv` into an identical `runtree` catalog so the summed grad/Hess
  are meaningful. Add a validation check.
- **Primary measurement — within-B degeneracy breaking.** Produce two fits over
  the *same* B⁺→J/ψK⁺ events: (i) **dimuon-only** (jpsicons) and (ii)
  **dimuon + kaon** (sum the two makers' grad/Hess). The kaon is a
  statistically independent track (different hits) on the same events, so the
  sum is clean and involves no double-counting. Comparing (i) vs (ii) isolates
  the kaon's effect with zero sample-composition confound and maximal
  sensitivity.
- **Sanity check — J/ψ consistency.** Separately run the nominal standalone
  J/ψ→μμ calibration on `TkAlJpsiMuMu_*.root` (Run2016H) and the J/ψ-from-B
  dimuon-only calibration, and confirm their b-field/material corrections
  agree. This validates that the B dimuons behave like the standard J/ψ sample,
  so the primary within-B comparison is meaningful.
- **Enable covariance output** in the solve (currently `docov=False`, errors
  returned as zeros) so uncertainties and off-diagonal correlations are
  available.
- **Hold alignment fixed** (`CVH_FIX_ALIGNMENT=1`) as the nominal configuration,
  so only the 92 global parameters float — matching how the reference fit is
  run. A measured diagnostic (design.md Decision 5) shows floating the 81,732
  alignment parameters inflates the b-field uncertainties ~2× through
  marginalisation and *suppresses* the field↔material correlation.
- **Report per-parameter constraint σ**, separating the two field regimes:
  **mode 0 = scale** (expected to stay with the J/ψ) from **modes 1+ = shape**
  (where a soft displaced hadron can plausibly contribute), plus the material
  groups (parmtype 15). Magnitude, uncertainty, and the field↔material
  correlation are reported for dimuon-only vs dimuon+kaon, plus the
  J/ψ-consistency sanity check.
- **Species closure check**: rerun the same bachelor tracks with the *muon*
  mass hypothesis (`kaonAsMuon=True`, already wired in the driver) and compare.
  This isolates the dE/dx/species effect from everything else — a direct test
  of the energy-loss model at fixed track sample.
- **Produce detailed slides** once the study is complete, walking through the
  whole procedure (three-stage pipeline, grouped-material/field-mode
  parameterization), the iterations and decisions taken along the way, the
  sanity check, and the final results.

## Impact

- Affected specs: new capability `btojpsik-cvh-calibration`.
- Affected code:
  - `CMSSW_15_0_19_patch2/src/Analysis/HitAnalyzer/python/ResidualGlobalCorrectionMakerTwoTrackJpsiKMuMuG4e_cfi.py`
  - `CMSSW_15_0_19_patch2/src/Analysis/HitAnalyzer/python/ResidualGlobalCorrectionMakerJpsiKSingleTrackKaonG4e_cfi.py`
  - `CMSSW_15_0_19_patch2/src/Analysis/HitAnalyzer/test/runCvhBplusJpsiK.py`
  - New fork: `/work/submit/pmlugato/REALmz/mz/real/calibration/` (aggregation, solve, reporting; not versioned in this repo)
  - New driver/reporting scripts under `scripts/btojpsik/calibration/` (thin wrappers + plots that live in this repo)
- Inputs: `/ceph/submit/data/user/p/pmlugato/mz/alcareco/full_2016postvfp/charmonium/Run2016H/` (`TkAlJpsiMuMu_*.root`, `TkAlJpsiX_*.root`).
- No production correction file is released; outputs are study artifacts (`correctionResults_*` and comparison plots/tables).
