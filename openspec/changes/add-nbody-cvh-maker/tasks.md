# Tasks — N-body CVH maker

## 0. Design

- [ ] 0.1 Fix the constraint representation: list of
  `(track subset, mass, width)` replacing the boolean `doMassConstraint`.
- [ ] 0.2 Decide what is shared with the two-track maker (helpers/headers)
  vs duplicated, keeping the two-track file untouched.
- [ ] 0.3 Fix the per-candidate output contract (mirror the two-track
  ValueMaps, N-general per-leg quantities).

## 1. Plugin

- [ ] 1.1 New `ResidualGlobalCorrectionMakerNTrackG4e.cc`.
- [ ] 1.2 State layout `nstateparms = 5*ntracks + 5*nhits`.
- [ ] 1.3 N-track candidate decomposition (reuse the generic nested-VCC
  descent rule already landed in the two-track maker).
- [ ] 1.4 Subset mass constraints.
- [ ] 1.5 Per-candidate ValueMaps incl. `jacMass` / `globalIdxs` /
  factored Hessian.
- [ ] 1.6 Per-channel cfis for the 3-track channels.

## 2. Validation

- [ ] 2.1 **N=2 closure against the two-track maker** on the same events.
- [ ] 2.2 3-body fitted `m(mu mu K)`: peak position and resolution vs the
  raw four-vector sum.
- [ ] 2.3 Comparison against refit-tracks +
  `KinematicConstrainedVertexFitter` on identical events.
- [ ] 2.4 Resource cost vs the 2-track + 1-track pair.

## 3. Follow-up

- [ ] 3.1 4-body validation (code is N-general).
- [ ] 3.2 Feed `jacMass` into a B-mass calibration.
