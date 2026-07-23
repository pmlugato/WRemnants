## Context

The CVH producer investigation
(`slides/cvh_vs_alcareco_jpsi_mass_investigation.pdf`) surfaced a
σ~80 MeV dimuon-mass broadening under `useIdealGeometry=True`, and
observed that flipping to `useIdealGeometry=False` recovers σ~40
MeV. The AN2021_131_v8 §3.3–3.4 review reframes that observation:
the σ~80 MeV signature is what CVH looks like when Stage 1 of the
AN's documented two-stage pipeline runs *without* Stage 2 loaded.
Stage 2 is a linearised correction applied via
`corFiles = cms.vstring(<correctionResults_*.root>)` whose
alignment-like term absorbs the ~100 μm displacement between
design and aligned geometry. Files exist on disk under
`wremnants-data/data/calibration/` (`v718`, `v721`,
`lbl2017`, `lbl2018` variants).

Two development checkouts exist:

- `CMSSW_10_6_26/` (WMass fork) — used for the debug study only.
  Now redundant.
- `CMSSW_15_0_19_patch2/` — the target production checkout.

Two open loops from the previous proposal
(`improve-cvh-refit-convergence`) remain relevant *only* as
context — the ideal-geometry-vs-corrections diagnosis falsified
the algorithmic hypotheses that motivated them:

- `nIters` / `edmConvergence` are exposed as configurable
  ParameterSet inputs (task 1.1–1.3 of that proposal), but the
  matrix sweep is no longer necessary — under the AN-canonical
  config we expect cap-hit fraction ≈ 1% and median `chisqval` ~ 30.
  Defaults stay.
- `useStartingState="midPropagated"` was deferred (task 1.4) and
  should remain deferred: no evidence it's needed.

David W's fork (`davidwalter2/cmssw@62becafa038`) still has 7
commits ahead of the WMass tag; those may be useful (candidate-
driven driver, V0CandidateProducer removal) but adopting them
is postponed until we push our own tree and take a diff.

## Goals / Non-Goals

**Goals:**
- Make `CMSSW_15_0_19_patch2/` the single home for CVH producer
  development going forward.
- Adopt the AN-canonical CVH config in btojpsik drivers:
  `useIdealGeometry=True` + `corFiles=[correctionResults_v721_recjpsidata.root]`.
  Scale-consistent with the published mZ measurement.
- Close out the alcareco-producer validation with three remaining
  frames: (P1/P2/P3 comparison, (A) vs (B) direct overlay, four-flavor
  B+ mass overlay) and the AN-two-stage explainer.
- Inventory all launchers touching this producer so we know
  the blast radius of the driver-level wiring change and can
  flag any non-btojpsik launcher whose config looks suspect.
- Confirm no GT mismatch between the CVH job and the ALCAReco
  input at RECO time.

**Non-Goals:**
- Adopting David's fork or candidate-driven driver — separate
  follow-up.
- Optimising per-file join overhead — separate proposal.
- Any downstream Stage-3 integration — none exists yet.
- Full 100k-event production reprocessing — a
  production-management change-id, not this one.
- Investigating whether `useIdealGeometry=True` is ever
  physically correct for alignment training use cases —
  we leave that path available via explicit override.

## Decisions

### Decision 1: Adopt AN-canonical config in btojpsik drivers; cfi defaults untouched

- **What**: In btojpsik-scoped drivers, keep
  `useIdealGeometry=True` and populate
  `corFiles=[correctionResults_v721_recjpsidata.root]`. Wire
  both as VarParsing options so operators can override per-job
  without editing sources. The maker cfi defaults
  (`ResidualGlobalCorrectionMakerTwoTrackG4e_cfi.py`,
  `ResidualGlobalCorrectionMakerG4e_cfi.py`) stay unchanged:
  `useIdealGeometry=cms.bool(True)` and `corFiles=cms.vstring()`.
- **Why**: Three rationales.
  1. **Scale consistency with the published mZ measurement**.
     AN2021_131_v8 §3.4 states explicitly that the *ideal
     geometry version of the track fit with truth-assisted
     global corrections* is used as the nominal simulation
     reference. Adopting (A) keeps our btojpsik pipeline aligned
     with that reference so any downstream comparison against
     the mZ scale works without a config-level re-architecture.
  2. **Correction files exist already**.
     `correctionResults_v721_recjpsidata.root` is committed to
     `wremnants-data/data/calibration/`; no new derivation is
     needed. Rejecting (A) on plumbing cost isn't justified.
  3. **Cfi defaults preserved for non-btojpsik consumers**.
     Empty-`corFiles` cfi default is the appropriate default for
     an alignment-training consumer whose *purpose* is to
     produce Jacobians for a subsequent alignment fit, not to
     load pre-derived corrections. Modifying the cfi to require
     `corFiles` populated would silently break those consumers.
- **Alternatives considered**:
  - *(B) `useIdealGeometry=False + corFiles=[]`*. Rejected —
    physics-quality out-of-box but decoupled from the mZ scale
    reference. If we ever want to compare our btojpsik scale to
    the mZ measurement, (B) requires an extra alignment-scale
    reconciliation step.
  - *Blanket cfi flip*. Rejected — silently breaks non-btojpsik
    consumers.
  - *Do nothing (keep current empty-`corFiles`)*. Rejected — we
    demonstrated that config produces σ~80 MeV / 33% tail.
- **Risk**: The `correctionResults_v721_recjpsidata.root` file
  may not be the appropriate variant for our Run2016H sample
  (the AN's mZ analysis is UL 2016–2018; v721 may be a specific
  vintage). Mitigation: §3b (A)-vs-(B) direct overlay catches
  any disagreement empirically; if (A) and (B) don't agree on
  our sample, Open Question #4 is raised as a follow-up to
  identify a per-era file (`lbl2016...` if it exists, or a
  derivation request).

### Decision 2: Keep the previous proposal's iter/edmConvergence knobs

- **What**: Leave the `nIters` / `edmConvergence` / `debugPerIterDump`
  parameters in place (already partially implemented via
  `improve-cvh-refit-convergence`), but do not run the matrix
  sweep.
- **Why**: The knobs are cheap to keep, they cost nothing at
  runtime with defaults, and they leave the diagnostic path
  open for future investigations. Ripping them out would create
  merge churn without value.

### Decision 3: Per-leg knobs for both `useIdealGeometry` and `corFiles`

- **What**: Add `useIdealGeometryMuon`, `useIdealGeometryKaon`,
  `corFilesMuon`, `corFilesKaon` as VarParsing options in
  `runCvhBplusJpsiK.py`. Each defaults to the corresponding
  base value (`useIdealGeometry`, `corFiles`) if unset.
- **Why**: §3 sweeps three primary points (P1 AN-canonical, P2
  aligned baseline, P3 broken baseline) and a 2×2 muon × kaon
  sub-scan within P1. Independent per-leg knobs let this happen
  in a single driver on a single event sample.
- **Alternatives**: Two runs with the module cfi patched manually.
  Rejected — noisy diff history and error-prone.

### Decision 4: Publication surface is the existing
`cvh_producer_validation.tex` deck, not a new one

- **What**: Extend the current 13-frame deck with §3 and §4
  frames; do not spawn a fresh deck.
- **Why**: The deck already tracks the validation story
  end-to-end; adding two frames keeps context together and
  avoids the "which deck is current?" problem.

### Decision 5: Publish figures to `~/public_html/mz/cvh/` in
addition to the LaTeX-embedded PNG

- **What**: Every new figure lands both in
  `slides/figs/cvh_producer_validation/` (LaTeX include path)
  AND in `~/public_html/mz/cvh/`.
- **Why**: David and the rest of the group need a shareable
  URL without cloning the repo; the LaTeX include is for the
  deck; both come from the same PNG.

### Decision 6: GT verification is per-era, not per-lumi (first pass)

- **What**: Compare the driver's per-era `GlobalTag.globaltag`
  string against the GT observed in one ALCAReco file per era.
- **Why**: Per-lumi GT variation is rare and only matters if
  the alignment IOV rolled inside the era. First pass catches
  the coarse case; escalate only if a mismatch surfaces.
- **Escalation path**: If the per-era check flags a mismatch,
  read the `ProcessHistory` per-lumi and identify the exact
  alignment record IOVs used, then choose the driver GT that
  covers all of them.

## Risks / Trade-offs

- **`correctionResults_v721_recjpsidata.root` may not match our era**:
  the AN's mZ analysis is UL 2016-2018; the v721 tag suggests a
  specific derivation vintage. If it was tuned on a different
  era's alignment, our Run2016H sample could see a residual
  mismatch. Mitigation: §3b (A)-vs-(B) direct overlay catches
  this empirically; disagreement raises Open Question #4 as
  a follow-up before full-scale production.

- **Producer's own `corFiles` handling vs downstream Jacobian
  application**: it's unclear whether populating `corFiles`
  makes the producer apply the correction to output branches
  or whether the Jacobians are only written for downstream
  analysis-time consumption. Mitigation: §2.7 smoke run
  answers this empirically (if the mass distribution is
  physics-quality with just `corFiles=[v721]`, the producer
  applies; if not, we need a downstream step).

- **Silent config drift on non-btojpsik consumers**: A
  downstream consumer we don't grep for could still be running
  `useIdealGeometry=True + corFiles=[]` (the debug-broken
  config) unintentionally. Mitigation: §7 inventory + the new
  `LogInfo` at construction. Every job now announces both
  flags.

- **David's fork divergence**: We postpone adopting his
  candidate-driven driver. Risk: he keeps evolving, and our
  diff grows harder to reconcile. Mitigation: user's stated
  plan — push our tree first, diff soon after.

- **Preset-B collapse in §4**: Under our operational default
  (preset B), the "raw" and "ALCAReco output" mass flavors
  coincide (no kinematic constraint applied). Trade-off: keep
  both flavors on the plot for the record, with the caption
  noting they collapse.

- **`LogInfo` verbosity**: Adding one `LogInfo` line per
  construction is negligible, but if the maker is instantiated
  many times per job (per-lumi reloads etc.), we could get
  noise. Mitigation: log at construction is one-shot per job,
  not per event; safe.

## Migration Plan

1. Land §1 (consolidation notes) and §2 (default flip + LogInfo)
   first; verify smoke test still shows σ ≈ 40 MeV / 0 tail.
2. Land §3–§4 (validation plots) next; extend the deck.
3. §5–§7 (publication, GT check, inventory) as a single closing
   PR after the plots exist.
4. Rollback plan: `git revert` the cfi-default and driver-default
   flips restores `True`, no data-format change. §3–§4 plot
   scripts are additive.

## Open Questions

1. **RESOLVED (2026-07-02)**: adopt (A) AN-canonical. btojpsik
   drivers set `useIdealGeometry=True + corFiles=[v721]`; cfi
   defaults untouched.
2. **Per-lumi vs per-era GT check**: is per-era sufficient?
   Assumed yes; escalate on mismatch.
3. **Raw vs ALCAReco output mass under preset B**: keep both
   flavors on the four-flavor overlay for the record, or
   collapse them? Assumed keep both, note the collapse.
4. **Is `correctionResults_v721_recjpsidata.root` right for our
   Run2016H sample?** Answered empirically by §3b (A)-vs-(B)
   overlay. If disagreement, per-era file selection follow-up.
5. **Does the producer apply the `corFiles` correction to output
   branches, or only write Jacobians for downstream consumption?**
   Answered by the §2.7 smoke run (LogInfo dump + mass
   distribution). Determines whether §4 plots need an extra
   downstream Jacobian-application step.
