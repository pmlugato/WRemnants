# Change: Finalize preset-B selection and produce the 100k headline sample

## Why

Preset B is the operational AlCaReco default for the J/ψ + X stream and has
been validated end-to-end in `add-jpsi-x-condor-production`. Three small
clean-ups remain before we treat it as production-final:

1. The J/ψ mass constraint applied by `JpsiXCandidateProducer` runs under
   both preset B and preset C (`JpsiXCandidateProducer.cc:359–411`,
   call sites at `:644` vcc-mode and `:573-ish` track-mode). Under preset B
   there is no Kalman B-vertex fit; the constraint is the *only* kinematic
   refinement the producer applies, and it pins m(μμ) to PDG with σ = 1 keV
   — i.e. exactly. For the alignment use case we want preset B to deliver
   **raw, track-only kinematics** (mass window selection, no kinematic
   constraint) so downstream Stage-2 CVH refits see unconstrained inputs.
   The constraint stays on under preset C, where it is intrinsic to the
   multi-track Kalman fit via `TwoTrackMassKinematicConstraint`.

2. The merged-collection track filter
   (`AlignmentTrackSelectorWithIndexMap`, `ALCARECOTkAlJpsiX_cff.py:375`)
   has `ptMin = 0.4 GeV`. Every per-candidate bachelor/daughter pT cut
   inside the stream is already at 0.1 GeV (sibling change
   `add-jpsi-x-vertex-fit-and-low-pt`), and the V0 clone runs at
   `tkPtCut=0.1`. The global 0.4 GeV cap throws away V0 daughter pions
   that survived everything upstream. Lower it to **0.1 GeV** so the
   downstream alignment fit sees the full soft-pT population the producers
   were tuned to keep.

3. The user-facing characterization of preset B as "mass-window only on
   the dimuon" needs to be made true and *checked* against what the code
   does. After change (1) it will be true. The V0 (Ks, Λ) leg already
   uses **window selection only** (`V0Producer`'s standard
   `kShortMassCut` / `lambdaMassCut`), with no mass constraint applied
   anywhere downstream — confirm and document so the same statement holds
   verbatim for the V0-mode channels.

4. The producer currently builds each daughter 4-vector from
   `(track.px, track.py, track.pz)` at **the track's own reference
   point** (CMSSW default: point of closest approach to the beam
   line — `JpsiXCandidateProducer.cc:573` for the bachelor,
   `TwoBodyDecayCandidateProducer` for the J/ψ muons). These reference
   points are *different spatial points along each daughter's helix*,
   not the common decay vertex. In a B field the track momentum
   *direction* rotates along the helix; the magnitude is conserved.
   For displaced decays (B+ cτ ≈ 491 μm, similar for Bc/B0/Bs/Λb) the
   per-track angle at the PCA-to-beamline differs from the per-track
   angle at the true mother decay vertex by ~track-arclength / ρ_helix
   ≈ 0.1–0.5 mrad for our kinematics. The resulting opening-angle
   mismatch enters the invariant mass through
   `m² = m₁² + m₂² + 2(E₁E₂ − |p₁||p₂|cosθ₁₂)`. Effect on m(B+) under
   preset B kinematics: **O(1–2 MeV) per-candidate bias on average,
   scaling linearly with the mother cτ instance value** (~5× bigger on
   maximally displaced candidates). This is **small for resolution**
   (1–2 MeV ⊕ 30 MeV ≈ 30 MeV) but it is a real **systematic that
   correlates with B kinematics** (Lorentz boost → cτ → displacement →
   bias), which matters for the alignment use case where Stage-2 ingests
   raw preset-B masses as a check that the per-track momenta are
   internally consistent.

   Fix: before computing `lvM = lvJpsi + lvBach`, run a **rough helix
   propagation** of each leaf track to a common point (the candidate's
   decay vertex estimate — already in hand as `vtxM`, the midpoint of
   the J/ψ vertex and the bachelor PCA, lines 582–585) and rebuild
   each daughter 4-vector at that point. This is a geometric
   propagation only — no fit, no covariance update, no χ² — using
   the existing `TransientTrack::trajectoryStateClosestToPoint(...)`
   machinery already used by the preset-C Kalman path. The Kalman
   fit under preset C does the same thing AND minimizes a vertex
   χ² AND updates covariances; preset B gets only the first piece,
   which captures the bulk of the mass-bias removal at a tiny fraction
   of the CPU cost.

Once those three are in, run preset B over **100 RAW files × 1k events =
100 k events** on the Run2016H Charmonium RAW dataset via Condor, in
parallel jobs, to produce the headline sample for the final talk. Write
up motivation → implementation → run → validation → projections in a
single cohesive deck: `slides/jpsix-producer-final.tex`.

## What Changes

- **MODIFIED** Preset B no longer applies the J/ψ → μμ
  `MassKinematicConstraint`. New `applyJpsiMassConstraint` cfg parameter
  on `JpsiXCandidateProducer` (default `False`); preset B sets it
  `False`, preset C sets it `True`. The constraint stays on under preset
  C because it is required for the multi-track Kalman fit's
  `TwoTrackMassKinematicConstraint` and used as the fallback dimuon p4
  if the multi-track fit fails.
- **MODIFIED** The dimuon mass window stays at the current ±5σ value
  (`minMass=2.95, maxMass=3.25`) — the user-permitted tightening is
  unused for now; the window is already tight relative to the ~30 MeV
  detector resolution, and tightening further would bite into truth
  signal without a corresponding combinatorial-reduction case.
- **MODIFIED** `AlignmentTrackSelectorWithIndexMap` `ptMin`: 0.4 → 0.1
  GeV. Every other pT cut in the stream is already 0.1; this aligns the
  final filter with the rest.
- **ADDED** Confirmation requirement that V0-mode channels use mass
  *windows*, not constraints, end-to-end: V0Producer's standard window
  cuts upstream, no `MassKinematicConstraint` or
  `TwoTrackMassKinematicConstraint` applied to the V0 sub-resonance
  anywhere in `JpsiXCandidateProducer.produceVccMode`. Verified by
  reading the code; locked by spec so it cannot regress.
- **ADDED** `JpsiXCandidateProducer` SHALL, under all presets,
  propagate each leaf track to a candidate-level common point along
  its helix before computing daughter 4-vectors and the mother 4-vector
  sum. Implementation uses two **existing** CMSSW propagator classes:
  - **`ClosestApproachInRPhi`**
    (`TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h`)
    for every 2-body sub-resonance — the two J/ψ muons, and (under vcc
    mode) the two sub-resonance daughters (Kπ for K*0, KK for φ, ππ
    for Ks, pπ for Λ). One `calculate(FTS, FTS)` call returns
    `trajectoryParameters()` — a `pair<GlobalTrajectoryParameters,
    GlobalTrajectoryParameters>` evaluated at the two helices' PCA
    points — plus `crossingPoint()` (their average) as the
    sub-resonance vertex estimate. We read each daughter's propagated
    `(p_x, p_y, p_z)` directly from its `GlobalTrajectoryParameters`.
  - **`AnalyticalImpactPointExtrapolator`**
    (`TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h`)
    for the bachelor track in track mode (B+ K, Bc π). One
    `extrapolate(FTS, GlobalPoint)` call propagates the kaon/pion
    helix to the dimuon crossing point (the J/ψ vertex estimate from
    the dimuon `ClosestApproachInRPhi`); we read the propagated
    momentum from the returned `TrajectoryStateOnSurface`'s
    `globalMomentum()`.

  No new helix-propagation code is written; both classes are
  production-tested CMSSW helpers. No fit, no covariance update, no
  iteration. Removes the ~0.1–0.5 mrad per-track angle bias at
  displaced vertices. Cost: 1 `ClosestApproachInRPhi` call per
  sub-resonance + 1 `AnalyticalImpactPointExtrapolator` call per
  track-mode bachelor — total ≈ 1 call per leaf, ~µs each, negligible
  vs the existing ~1 ms / surviving B+ candidate cost. Efficiency
  impact: zero (mass window is wide enough that the ~1 MeV mass shift
  never moves a candidate across the window boundary; geometric and
  DOCA cuts use raw track-level parameters that are unchanged).
- **NEW PRODUCTION** Run preset B over 100 RAW files × 1k events on the
  Run2016H Charmonium dataset via existing Condor scaffolding in
  `condor/jpsix_alcareco/`. Validate locally on 1–2 files first; submit
  the full 100-job batch only after the smoke files compile and run
  cleanly.
- **NEW DECK** `slides/jpsix-producer-final.tex` — cohesive final-state
  deck covering motivation, implementation, the three changes above,
  validation on the 100k sample, projections to postVFP 2016 F+G+H.

## Impact

- **Affected specs**: capability `alcareco-jpsi-x` (existing; MODIFIED
  preset-B contract).
- **Affected code**:
  - `CMSSW_10_6_17_patch1/src/Alignment/CommonAlignmentProducer/plugins/JpsiXCandidateProducer.cc`
    — gate `constrainJpsi4Momentum` on a new `applyJpsiMassConstraint`
    config knob; use `ClosestApproachInRPhi` on every 2-body
    sub-resonance pair (dimuon; vcc-daughter pair) to obtain both the
    pair-vertex estimate and each daughter's propagated
    `GlobalTrajectoryParameters`; use `AnalyticalImpactPointExtrapolator`
    on the bachelor in track mode to propagate it to the dimuon
    crossing point. Call sites: replace the bachelor-4-vector
    construction at lines 572–573 and the corresponding J/ψ p4 use,
    plus the vcc-mode mother-sum logic at line 651. (No change to the
    preset-C multi-track Kalman path — it already does its own
    correct per-track propagation to the fitted vertex.)
  - `CMSSW_10_6_17_patch1/src/Alignment/CommonAlignmentProducer/python/ALCARECOTkAlJpsiX_cff.py`
    — pass `applyJpsiMassConstraint=False` for preset B, `True` for
    preset C, on all 7 `JpsiXCandidateProducer` instances; drop
    `AlignmentTrackSelectorWithIndexMap.ptMin` to 0.1.
  - `condor/jpsix_alcareco/` — new submit file and per-job arg list for
    the 100-file × 1k-event batch.
  - NEW `slides/jpsix-producer-final.tex`.
- **Reused unchanged**: `TwoBodyDecayCandidateProducer` (already supports
  no-fit / window-only mode under preset B; no code change here),
  `V0Producer` clone (`ALCARECOTkAlV0Candidates`), the cosα/Lxy preset-C
  cut machinery, the closest-in-z PV selector.
- **No changes** to Stage-2 CVH refit code, narf/rabbit fit code, or any
  WRemnants downstream tooling. Preset C behaviour is unchanged.
