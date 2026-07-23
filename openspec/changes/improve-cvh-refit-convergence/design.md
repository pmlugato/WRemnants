# Design: improve CVH refit convergence

## Context

The CVH joint two-track refit in `ResidualGlobalCorrectionMakerTwoTrackG4e.cc`
implements a sparse Generalised Broken Lines (GBL) χ² minimisation
over a Geant4e-propagated trajectory: each material plane contributes a
$\chi^2$ term from the scattering/E-loss model and one or two hit
constraints (planar / pixel). The iteration is Gauss–Newton; per-iter
update solves $\delta x = -H^{-1} g$ with $H = 2 F_s^\top V^{-1} F_s$
and $g = 2 F_s^\top V^{-1} r$, sparse-factored via Eigen.

The current configuration has:

- iteration cap `niters = 10` (line 1418, with two muons + joint vertex constraint + optional J/ψ mass constraint)
- convergence test `edmval < 1e-5` where `edmval = -delta chi²` from the last Gauss–Newton step (line 3114)
- linearisation point: each muon's perigee parameters from the input
  ALCAReco track, refit'd into the FreeTrajectoryState at perigee
- `morehitquality = true` unconditionally (line 2197) — no per-hit
  rejection at the alignment-pass stage

The investigation (`slides/cvh_vs_alcareco_jpsi_mass_investigation.pdf`)
showed that with these settings:

- 0 % of events meet `edmval < 1e-5`
- 20 % hit the `niters = 10` cap; 80 % leave via some other break path
- per-leg $\Delta(q/p)$ is $\sim 10\times 10^{-3}$ GeV$^{-1}$ MAD vs the input
- 33 % catastrophic-outlier tail with $|m - 3.097| > 0.2$ GeV

The shape of the data does not by itself decide whether this is a
threshold problem (the fit is fine, it just never declares itself
converged), a step-size problem (the GN steps are too large and the
fit oscillates / diverges for some events), or a linearisation problem
(the iter-0 reference is too far from the true minimum so the first
step overshoots). The matrix is designed to distinguish these three.

## Why three knobs, not one

Each candidate fix addresses a specific failure mode that the per-iter
trace will fingerprint differently:

| Knob | Targets | Expected fingerprint if it works | Expected fingerprint if not |
|---|---|---|---|
| `nIters` | threshold-only failure | `edmval` decays monotonically, just slowly; bumping the cap pushes more events under threshold | bumping the cap does not change σ_m or tail fraction; `edmval` plateaus or grows after some iter |
| `edmConvergence` | unit mismatch on the threshold | relaxing to $10^{-2}$ produces same per-leg q/p as the unconverged 10-iter solution; "converged" fraction jumps to ~100 % | relaxing the threshold *and* keeping σ_m unchanged means the events were essentially at the minimum already and the broadening lives elsewhere |
| `useStartingState` | iter-0 linearisation point too far from the joint vertex | first-iter Δχ² shrinks; the 33 % tail population stops appearing | starting state has no effect → the failure is in the GN step itself |

This gives a 3-way classification of the convergence failure that a
single-knob test cannot.

## Why not just relax `edmConvergence`

A loose threshold + an unconverged fit gives the same per-leg q/p as
the current 10-iter cap output. The threshold is meaningful only if
the iter-by-iter trace shows the fit actually approaching a minimum
within the iteration budget. The per-iter dump in §1 of the proposal
distinguishes "converged but threshold tight" from "not converging."

## Why a 24-point matrix (not exhaustive)

The matrix is small enough to run on one host overnight (24 × ~5 min
≈ 2 hours), large enough to span the three failure-mode axes
independently. The `useStartingState` axis has only 2 values (cheap to
add the second; expensive to add a third) because the proposal only
proposes one alternative (midpoint-propagated). If neither value works,
the per-iter trace will give the next state idea explicitly rather than
us blindly searching.

## Why the cross-driver sanity check (§5 of proposal)

The CVH convergence story has so far been read off `runCvhBplusJpsiK.py`'s
output exclusively. That driver wraps the dimuon-side maker behind a
custom splitter (`JpsiKCandidateSplitter`) and reads `TkAlJpsiX` ALCAReco.
The legacy `runCvhJpsi.py` invokes the *same* maker on the *legacy*
`TkAlJpsiMuMu` ALCAReco, using the producer's in-module track-pair loop
instead of an external candidate input.

For a J/ψ candidate that exists in both ALCAReco streams, the dimuon
refit's output is a function of: the two input track parameters, the
beam-spot, the magnetic-field map, the particle list, and the producer's
configurable knobs (mass constraint, plimit, iter cap, etc.). None of
these *should* depend on which test driver loads the maker. If we
observe a per-branch difference between the two drivers, one of those
configuration channels is silently disagreeing — and that hidden offset
would contaminate the §4 experiment matrix, since the matrix is run only
through the `runCvhBplusJpsiK.py` chain.

This is a one-off check, not a recurring test: once the two drivers
agree on a shared candidate population, the matrix in §4 can be read
without worrying that the choice of driver is part of the story.

## Why nominal CMSSW field everywhere in this study

The investigation that motivated this change ran with
`useScalarPot3D=False` (nominal CMSSW field). Re-running the §4 matrix
with the scalar-potential or Opera-3D field could mask or move the
convergence pathology; any change in `chisqval` distribution could then
be attributed to the field choice rather than the iter/edm knob. We
pin the field model to the nominal CMSSW value throughout the matrix
and the cross-driver check so the only varied axes are the three knobs
in scope.

## Risk: changing default convergence in production

The matrix runs are on debug-mode jobs with the per-iter dump on; the
production default is changed *only* after the §8 validation gate in
`tasks.md` is met. The defaults in the producer constructor reproduce
the published 81948-event sample bit-identically (task 1.7) so that
the existing slides remain a meaningful baseline.

## Out of scope and why

- **Material model / Geant4e step size**. The per-leg Δ(q/p) is ten
  times the Kalman uncertainty. If the issue were just MS modelling
  the offset would scale with the integrated material (and with η);
  the bottom-right panel of `slides/figs/cvh_jpsi_investigation/03_dpt_per_leg.png`
  shows σ(Δ(q/p)) is roughly flat in p_T, weakly η-dependent — more
  consistent with a fit-numerics issue than a material modelling one.
  Revisit only if the matrix shows the per-leg shift survives all three
  fixes.
- **B-field map**. Same argument: an off-magnitude B map would shift
  $\langle q/p\rangle$, not just inflate $\sigma(q/p)$; the data show no
  per-leg bias (median $\Delta p_T/p_T \approx 0$). Revisit only if the
  per-iter trace reveals an iter-0 q/p shift that doesn't depend on the
  starting state choice.
- **Single-track CVH refit on the kaon**. The kaon side is not implicated
  in the dimuon broadening. The proposal touches the kaon producer only
  to clean up the always-zero `nValidHitsFinal` branch (task 2).

## Open question on `edmval` units

Reading `ResidualGlobalCorrectionMakerTwoTrackG4e.cc:2469`:
```
const double deltachisq = (rfull.transpose() * VinvF * dxfree)(0, 0);
```
This is a scalar with the dimensions of the χ² scale, *not* normalised
by ndof. Standard GBL convergence threshold (`gbl::GblTrajectory::fit`)
is on $\Delta\chi^2 < 10^{-2}$ for a *single hit's chi²* — orders of
magnitude looser than `1e-5`. The 1e-5 threshold may simply be a
copy/paste from a different optimiser convention. The matrix test of
{1e-5, 1e-3, 1e-2} brackets this exactly.

## Matrix findings (2026-06-29, all 12 points complete)

Sweep on 200 events / one Run2016H ALCAReco file, with debug=True so
the per-iter trace is in the output. Mass-constrained branches
(`Jpsicons_mass`, `Jpsikin_mass`) are bit-identical across all points;
only the unconstrained `Jpsi_mass` reacts.

| tag                | sigma_wide MeV | sigma_core MeV | tail% | cap-hit % |
|--------------------|---------------:|---------------:|------:|----------:|
| n10_e1em5_perigee  |         79.0   |         27.7   | 28.2  |   18.7    |
| n10_e1em3_perigee  |         78.9   |         27.6   | 28.2  |   10.6    |
| n10_e1em2_perigee  |         78.9   |         27.6   | 28.2  |    8.0    |
| n20_e1em5_perigee  |         79.0   |         26.9   | 28.2  |    9.6    |
| n20_e1em3_perigee  |         78.9   |         26.9   | 28.2  |    7.8    |
| n20_e1em2_perigee  |         78.9   |         26.8   | 28.2  |    6.5    |
| n50_e1em5_perigee  |         79.0   |         27.4   | 28.2  |    9.1    |
| n50_e1em3_perigee  |         79.0   |         27.3   | 28.2  |    7.5    |
| n50_e1em2_perigee  |         79.0   |         27.3   | 28.2  |    6.0    |
| n100_e1em5_perigee |         79.1   |         27.7   | 28.2  |    9.1    |
| n100_e1em3_perigee |         79.0   |         27.7   | 28.2  |    7.0    |
| n100_e1em2_perigee |         79.0   |         27.7   | 28.2  |    6.0    |

**Both knobs falsified**:

1. **`edmConvergence`** at 10⁻⁵ → 10⁻² shifts σ_m by <0.2 MeV and leaves
   converged% at 0%. The typical `edmval` lands at $\sim 10^5$, so the
   threshold is unreachable at any of these values. It's a bookkeeping
   cut, not a step-size lever.

2. **`nIters`** at 10 → 20 cuts cap-hit events almost in half (18.7% →
   9.6%) — some previously-cap-hit events now "converge" between iter
   11 and 19 — but $\sigma_m$ on [2.95, 3.25] and the 28.2% tail
   fraction are unchanged. The events that no longer hit the cap are
   landing at the SAME bad solution they would have at iter 10; the
   extra iterations don't move them to a better minimum.

The per-iter trace explains why: bad events' dimuon mass *progressively
walks away from PDG* through the unconstrained icons=0 phase, sometimes
out to 6 GeV by iter 9. The Gauss–Newton step is too large for the
local curvature on these events; adding more iterations just lets the
divergence run further.

**Conclusion**: the convergence failure is **step-size instability** in
the GN update, not iteration count or threshold. The cure is a
step-size control mechanism (Marquardt damping, trust-region step, or
line search) — out of scope for this change; subject of the next
proposal.

The n50 and n100 rows are likely to confirm the trend without changing
the conclusion. They're left running for completeness but the
investigation can already act on the n10/n20 verdict.

## Preliminary per-iter findings (2026-06-29, partial matrix)

The first matrix point (`n10_e1em5_sperigee`, baseline-equivalent) on 200
events of one Run2016H ALCAReco file with `debugPerIterDump=True` shows
the per-iter `Jpsi_mass` trace doing two qualitatively different things:

- For events with final `Jpsi_mass` close to PDG (best 100, "good"), the
  per-iter mass stays in a tight band around 3.097 GeV through all
  iterations of the unconstrained phase.
- For events with final `Jpsi_mass` far from PDG (worst 100, "bad"), the
  per-iter mass **wanders progressively further** with each iteration —
  reaching values as large as 6 GeV by iter 8–9 — until the mass-
  constrained icons=1 phase yanks them back to PDG by construction.

This is direct evidence that the Gauss–Newton step is *too large* for the
local curvature on these events: the fit is in a regime where successive
linearisations drift further from the true minimum. The implications for
the matrix are:

- Bumping `nIters` from 10 to 100 will not reduce the tail; bad events
  would simply diverge further before icons=1 catches them.
- Relaxing `edmConvergence` is a bookkeeping change that doesn't touch
  the step-size dynamics — the tail fraction should be unchanged.
- The cure is a step-size control mechanism (Marquardt damping,
  trust-region step, or line search). That is the natural next OpenSpec
  change after this one closes.

## Failure mode if the matrix doesn't find a fix

The per-iter trace from §5 of tasks gives the next investigation
direction explicitly:

- If `chisqval` oscillates between iterations → the Gauss–Newton step
  is too large for the curvature; add a Marquardt damping factor or a
  line search.
- If `chisqval` decreases monotonically but plateaus far above
  `ndof` → the model is wrong (wrong B field, wrong material, wrong
  hit covariance scaling). Cross-check against `Jpsikin` q/p, which is
  free of the CVH side.
- If `chisqval` is fine but `edmval` doesn't go to zero → the EDM
  computation has a bug. Recheck the sparse identity that gives the
  formula at line 2467–2469.
- If the per-iter q/p trace shows the q/p starts at the input value at
  iter 0 then immediately leaves it at iter 1 → the iter-0 reference
  state is being corrupted before the first solve.

These are diagnoses, not fixes. Each becomes the next proposal.
