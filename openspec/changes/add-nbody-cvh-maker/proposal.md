# Change: N-body CVH maker (new plugin, 3-body first)

## Why

The B mass we can build today is not a fit. Stage-1 sets the candidate
`p4()` as a **raw four-vector sum** ("raw track-based", deliberately, so
the CVH refit sees unconstrained inputs), and the two-track maker's
`corMass` is the refit **dimuon** mass, not a B mass. So nothing in the
chain currently produces a fitted `m(mu mu K)`.

Two routes give one, and they are **not** substitutes:

1. **Refit tracks -> `KinematicConstrainedVertexFitter`** (the Bmm5
   pattern). Cheap, and lands in the sibling change
   `add-alcareco-to-nanoaod` via `emitRefitTracks`. Gives an
   analysis-quality fitted mass and vertex.
2. **A joint N-body CVH fit** (this change). Fits all N tracks together
   through Geant4e with a common vertex and subset mass constraints.

Only route 2 yields **`jacMass` = d(mass)/d(global correction params)**,
alongside `globalIdxs` and the factored Hessian. That derivative is what
makes a mass usable as a **calibration observable**. Route 1 is a second,
independent fit that knows nothing about the CVH global parameters, so it
can never feed a B-mass-based calibration.

If we want to calibrate on the B mass, route 2 is the only path.

## What Changes

**A new plugin file.** We do *not* generalize
`ResidualGlobalCorrectionMakerTwoTrackG4e`. Keeping it untouched avoids a
bit-identity regression burden on a fitter that is actively developed
upstream, and avoids merge friction. The new maker may share helpers, but
it is its own file and its own module.

### Feasibility (assessed against the existing two-track code)

- **The fit core already generalizes.** The linear algebra is dynamic
  Eigen sized off `nstateparms`
  (`MatrixXd::Zero(ncons, nstateparms)`) -- nothing in the GBL/Geant4e
  math is 2-specific.
- **The dimensioning is one literal.** Two-track uses
  `nstateparms = 10 + 5*nhits`, where the 10 is the two-track vertex PCA
  block (`idx 0..9`); single-track uses `5*(nhits+1)`. The N-body form is
  `5*ntracks + 5*nhits`.
- **The 2-ness is bookkeeping**: ~118 fixed-size sites
  (`std::array<...,2>`, `[2]`, `muPlus`/`muMinus` naming, per-leg outputs).
- **The one genuine design piece** is the constraint structure. Today
  `nicons = doMassConstraint_ ? 2 : 1` means "the two tracks". N-body
  needs constraints on a **subset**: for B+ -> J/psi K+, a J/psi mass
  constraint on 2 of the 3 tracks, plus optionally a mother-mass
  constraint. Constraints become a list of `(track subset, mass, width)`.

### Scope

3-body first, but written **N-general** so 4-body is configuration, not a
second refactor. Channel coverage:

| N | Channels |
|---|---|
| 3 | BuJpsiK (mumuK), Bc (mumupi), D*->D0pi |
| 4 | B0K*0, BsPhi, B0Ks, Lambda_b, psi(2S) |

3-body alone therefore covers 3 of the 12 persisted channels; the
J/psi + (2-body X) channels all need N=4.

### Outputs

Same contract as the two-track maker: per-candidate `ValueMap`s for the
fitted mass, its error, kinematics, per-leg refit quantities, `edmval`,
and -- when the global-fit payload is on -- `globalIdxs`, `jacRef`,
`jacMass`, factored Hessian.

## Impact

- New: `Analysis/HitAnalyzer/plugins/ResidualGlobalCorrectionMakerNTrackG4e.cc`
  plus per-channel cfis.
- Unchanged: the two-track and single-track makers.
- New capability spec `nbody-cvh-maker`.

## Validation

- **N=2 closure**: configured for two tracks, the new maker must
  reproduce `ResidualGlobalCorrectionMakerTwoTrackG4e` on the same
  events, within the documented tolerance. This is the regression test
  that makes the whole thing trustworthy.
- **3-body physics**: fitted `m(mu mu K)` peaks at the B+ mass with a
  resolution better than the raw four-vector sum.
- **Comparison against route 1**: joint 3-body CVH vs
  refit-tracks + `KinematicConstrainedVertexFitter` on identical events
  -- agreement of the central mass, and the resolution difference. Note
  the discriminator is not resolution but `jacMass`, which only route 2
  provides.

## Non-goals

- Generalizing or refactoring the existing two-track maker.
- 4-body validation (the code is N-general; validating it comes later).
- Wiring the result into a calibration fit.
