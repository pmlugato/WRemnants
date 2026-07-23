# Design: Finalize preset B and produce the 100k sample

## Context

Preset B is the operational AlCaReco default for the J/ψ + X stream.
Three earlier changes
(`add-jpsi-x-selection-comparison`, `add-jpsi-x-vertex-fit-and-low-pt`,
`add-jpsi-x-condor-production`) settled the per-candidate cuts, the
Condor scaffolding, and the 5-file × 10k = 50k-event Run2016H
production. What remains is a small set of clean-ups that surfaced
during the final deck review, plus a larger statistically-significant
final run.

The three clean-ups are independent of each other but live in the same
two files (`JpsiXCandidateProducer.cc` + `ALCARECOTkAlJpsiX_cff.py`),
so we batch them in one proposal to minimize cmsRun configuration
churn.

## Goals / Non-Goals

**Goals**:
- Make preset B "raw track kinematics, mass-window selection only" —
  no implicit kinematic constraint between channels.
- Preserve preset C behaviour bit-for-bit (the J/ψ mass constraint
  stays on under C because it is required by the Kalman multi-track
  fit).
- Preserve V0-mode channel behaviour bit-for-bit (no constraint
  anywhere on the V0 path, today or after this change).
- Capture a 100 k-event headline sample for the talk, in parallel for
  speed, with the finalized preset B.

**Non-Goals**:
- No new fit kernels or selection variables. The preset B physics is
  done.
- No changes to Stage-2 CVH (`add-jpsi-x-stage2-bplus-cvh`) — that's a
  separate change.
- No tightening of the dimuon mass window beyond what's already in the
  cff. The user gave permission to tighten ("can be relatively tight")
  but it's unnecessary: ±5σ is already tight and tightening would only
  reduce truth efficiency without a paired combinatorial-reduction
  argument.

## Decisions

### Decision 1: How to gate the J/ψ mass constraint

A new boolean cfg parameter `applyJpsiMassConstraint` on
`JpsiXCandidateProducer`, default `False`.

- Default `False` so that any caller that doesn't know about the flag
  (re-uses of the producer in other contexts) gets the raw-kinematics
  behaviour the alignment use case wants. This matches the producer's
  documented convention "Optional quality parameters (all
  backward-compatible, default = no cut)" at line 32.
- The cff sets it `True` only under preset C, alongside the other
  Kalman-fit-related parameters.
- The constraint is structurally needed in the multi-track Kalman fit
  for preset C (`TwoTrackMassKinematicConstraint` line 459). We do NOT
  decouple the two — under preset C the dimuon constraint and the
  Kalman fit are one physics package.

**Alternatives considered**:
- Hard-code `False` and delete the constraint code. Rejected: preset C
  still needs it as a fallback for the multi-track-fit-failed path,
  and we want preset C bit-identical to its previous behaviour.
- Gate on the existing `applyBVtxFit_` flag. Rejected: that flag means
  "do the Kalman fit"; the dimuon mass constraint is conceptually a
  separate physics object, and the multi-track fit owner may want to
  toggle them independently in future studies.

### Decision 2: Where to lower the track pT floor

Single line change in `ALCARECOTkAlJpsiX_cff.py`:
`AlignmentTrackSelectorWithIndexMap.ptMin = 0.4` → `0.1`.

This is the *final* filter applied to the merged-and-deduplicated
track collection (`ALCARECOTkAlJpsiXAllTracks`) before output. Every
upstream per-candidate pT cut is already at 0.1 GeV (per
`add-jpsi-x-vertex-fit-and-low-pt`), and the V0 producer clone runs at
`tkPtCut=0.1`. The global 0.4 GeV cap throws away V0 daughter pions
that survived everything upstream, defeating the purpose of the lower
V0 cuts.

**Trade-off**: more saved tracks per event → larger output. Estimated
~50% more tracks (the V0-daughter soft tail) → ~150 MB extra per 100k
events. Manageable within the AlCaReco size budget.

### Decision 3a: Rough helix propagation to a common point before mass sum

#### Why a track is a helix and where the 5 vs 6 mismatch comes from

A charged particle in a uniform magnetic field follows a helix: circular
motion in the transverse plane (radius ρ = p_T / 0.3 B, in m·GeV·T units,
so for B = 3.8 T a 5 GeV track has ρ ≈ 4.4 m) plus uniform motion along
the field axis. The geometric helix has **5 free parameters**: e.g.
(curvature, dip, azimuth, transverse impact parameter, longitudinal
impact parameter), or equivalently the CMSSW perigee set
`(q/p, λ, φ₀, d₀, z₀)`. Those 5 parameters fix the helix as a geometric
locus.

But a real particle has **6 phase-space coordinates** at any given
instant: position `(x, y, z)` and momentum `(p_x, p_y, p_z)`. The
"missing" 6th parameter is **where along the helix the particle
currently is** — equivalently, the arclength `s` along the helix
trajectory measured from some chosen reference point. The 5 helix
parameters describe the *locus*; the 6th coordinate locates the
*particle* on that locus.

Standard CMSSW convention resolves the 6th coordinate by reporting the
track at the **point of closest approach (PCA) to the beam line** —
i.e. by choosing `s = s_PCA` such that the transverse distance from the
helix to the beamline is minimized. The track's stored
`(px, py, pz, vx, vy, vz)` are therefore the 4-momentum and position at
*that* point on the helix, not at any other point.

#### What changes as you propagate along the helix

In a B field, the momentum **magnitude** `|p|` is conserved along the
helix (in the absence of energy loss). The momentum **direction**
rotates: in the transverse plane, the direction rotates by angle
`Δφ = Δs / ρ` for arclength `Δs`. The longitudinal `p_z` is unchanged.
So as you move from one reference point to another along the same
helix, you get the **same `|p|`** but a different `(p_x, p_y)`
re-projected by the rotation angle.

For our kinematics:

| Channel-typical pT | ρ [m] in 3.8 T | Rotation over 500 μm arc |
|---|---|---|
| 1 GeV (soft V0 daughter) | 0.88 | **5.7 × 10⁻⁴ rad ≈ 0.6 mrad** |
| 2 GeV (typical K from B+) | 1.75 | 2.9 × 10⁻⁴ rad ≈ 0.3 mrad |
| 5 GeV (typical muon from J/ψ) | 4.39 | 1.1 × 10⁻⁴ rad ≈ 0.1 mrad |
| 10 GeV (hard muon) | 8.77 | 5.7 × 10⁻⁵ rad ≈ 0.06 mrad |

#### Why this introduces a mass bias for displaced decays

For a J/ψ → μ⁺μ⁻ (followed by, say, B+ → J/ψ K+) the muons and the kaon
were all **born at the B+ decay vertex**, which sits at displacement
`r_B` ≈ cτ_B+ · (p_B / m_B) ≈ 491 μm × γβ ≈ a few hundred μm transverse
from the PV in a typical event.

The invariant mass of the final state depends on the 4-vector sum:

```
   m² = m_μ² + m_μ² + m_K²
        + 2 (E_μ⁺ E_μ⁻ − p_μ⁺ · p_μ⁻)
        + 2 (E_μ⁺ E_K  − p_μ⁺ · p_K )
        + 2 (E_μ⁻ E_K  − p_μ⁻ · p_K )
```

Each cross term `p_i · p_j = |p_i| |p_j| cos θ_ij` depends on the
**angle between the two daughter momentum vectors** at the spacetime
point where the calculation is performed. The Lorentz-invariant
quantity `m` does not depend on *where* you evaluate it — but only if
**all daughter 4-vectors are evaluated at the same spacetime point on
each particle's worldline**.

When the producer reads `(px, py, pz)` from each track at *that
track's own PCA-to-beamline*, the three muons-and-kaon's momenta come
from **three different points** along three different helices. The
opening angle between muon-i's `(px, py, pz)` at its PCA-i and muon-j's
`(px, py, pz)` at its PCA-j is not the opening angle the two particles
had at their common production vertex (the B+ vertex).

The angular mismatch is ~`Δs / ρ` per track, where Δs is the helix
arclength between each track's PCA-to-beamline and the true B+ vertex
(typically ~ the transverse displacement of the vertex, i.e. ~500 μm
for an average B+). For the kinematics above, that's **~0.1–0.5 mrad
per track**.

The mass derivative `dm/dθ_ij` for J/ψ → μμ at typical lab opening angle
~30° and dimuon mass 3.1 GeV is `≈ |p_μ⁺| |p_μ⁻| sin θ / m_J/ψ ≈ 2–5
GeV/rad`. So a 0.2 mrad angular mismatch shifts m(μμ) by
**~0.5–1 MeV**. The shift on m(B+) (a 3-body sum where the kaon also
has its own mismatch) is **~1–3 MeV** per candidate on average,
scaling linearly with the actual cτ value of the B+ in that event.

Average mass-resolution impact: 1–3 MeV ⊕ 30 MeV ≈ 30 MeV. **Resolution
gain is invisible.** Average mass-bias impact:
- ~0 MeV on prompt mesons (no displacement).
- ~1–3 MeV on typical-lifetime B+.
- ~5–15 MeV on the long-lifetime tail (cτ > 2× average).

The bias correlates with the B+ momentum (boost → cτ → displacement), so
it is a **kinematics-correlated systematic**, not a uniform offset. For
the alignment use case (Stage-2 CVH consumes raw preset-B tracks), it's
worth removing.

#### Why "rough propagation" cures it

The fix is to propagate each track from its native PCA-to-beamline to a
**common candidate-level point** along its own helix, and rebuild each
daughter's 4-momentum at that point. Two existing CMSSW classes do
exactly this — no custom propagator code is needed.

**For every 2-body sub-resonance pair** (the two J/ψ muons in all
channels; the two sub-resonance daughters in vcc mode — Kπ for K*0,
KK for φ, ππ for Ks, pπ for Λ): use
**`ClosestApproachInRPhi`** from
`TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h`. One
`calculate(FreeTrajectoryState_a, FreeTrajectoryState_b)` call solves
the two-helix PCA problem in the transverse plane analytically and
returns:
- `trajectoryParameters()` — a `pair<GlobalTrajectoryParameters,
  GlobalTrajectoryParameters>` evaluated at the two helices' PCA
  points. Each `GlobalTrajectoryParameters` carries `.momentum()`,
  i.e. the propagated `(p_x, p_y, p_z)` at that PCA, with the helix
  geometry handled exactly (the constructor solves for the new helix
  origin and rotates the momentum direction accordingly — see
  `ClosestApproachInRPhi::newTrajectory(...)` declared at line 76 of
  the header).
- `crossingPoint()` — the arithmetic mean of the two PCA points, used
  as the sub-resonance vertex estimate.
- `status()` — false if the calculation failed (e.g. concentric
  circles); we fall back to the un-propagated daughter 4-vectors in
  that case and increment a `n_propagation_fallback_` counter.

**For the bachelor track in track mode** (the kaon in B+, the pion
in Bc): use **`AnalyticalImpactPointExtrapolator`** from
`TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h`.
`extrapolate(FreeTrajectoryState, GlobalPoint)` propagates the
bachelor's helix to the 3D PCA with respect to a given target point
— in our case the dimuon crossing point we just obtained from the
J/ψ `ClosestApproachInRPhi`. The returned `TrajectoryStateOnSurface`'s
`.globalMomentum()` gives the propagated `(p_x, p_y, p_z)`.

Sum the propagated daughter 4-vectors to get the mother 4-vector.

`ClosestApproachInRPhi` is the same helix-helix PCA primitive that
`KalmanVertexFitter`, `V0Producer`, and several other CMSSW vertex
fitters use as their **initial seed step** before iterating on the
χ² fit. We are using exactly the seed step, without the χ²
iteration. Preset C does the full Kalman; preset B does only the
seed step — but that's where most of the mass-bias removal sits, at
a tiny fraction of the cost (~µs per call vs ~1 ms per full fit).

#### Efficiency impact

**Zero**, by construction:

- The mass window cut (5.0–5.5 GeV for B+, similar for other channels)
  is ~500 MeV wide. A 1–3 MeV shift moves no candidate across the
  boundary except at the extreme tails of the cτ distribution, where
  the candidates are already in the wings of the mass window for other
  reasons.
- Kinematic cuts (`minJpsiPt`, `minMotherPt`, `maxBachelorEta`,
  `minBachelorPt`) use **track-level** quantities (track `pt()`,
  `eta()`) — read directly from the track's helix parameters at the
  track's own reference point. These are unchanged.
- Geometric cuts (`maxBachelorMuTrackDOCA`, the 3D DCA to J/ψ vertex)
  also use raw track-level helix parameters. Unchanged.

So the candidate population entering Stage-2 and the merged-track
collection size are identical before and after this change. The only
difference is the value of `lvM.M()` stored on the mother — shifted by
~1–3 MeV per candidate toward the true B+ mass — and the corresponding
mother `lvM.pt()` / `lvM.eta()` (used only for the mother-pT cut, which
is currently a no-op at 0 GeV under preset B).

#### Alternatives considered

- **Do nothing.** Accept the ~1–3 MeV systematic. Rejected because the
  fix is trivial in code (one CMSSW helper call per leaf) and removes
  a kinematics-correlated bias that would otherwise have to be
  re-derived at the analysis stage.
- **Run the full Kalman vertex fit under preset B.** Rejected because
  it costs ~1000× more CPU per candidate and updates covariances we
  don't need at the AlCaReco stage (Stage-2 CVH does its own
  covariance via Geant4e). Preset C remains the choice for full
  fitting.
- **Propagate the bachelor to the dimuon midpoint instead of to the
  dimuon crossing point.** The midpoint of the muons' track reference
  points (the current `vtxM` definition at line 582) is a worse
  estimate of the B vertex than the dimuon helix-helix crossing point
  from `ClosestApproachInRPhi` — the midpoint is just half-way between
  two PCAs-to-beamline, whereas the crossing point is the actual
  helix-helix PCA in the transverse plane. We use the
  helix-helix crossing point.
- **Custom helix-propagation helper around
  `TransientTrack::trajectoryStateClosestToPoint(GlobalPoint)`.**
  Rejected in favour of the two existing CMSSW classes named above —
  fewer lines of code, more obvious physics intent (the dimuon PCA
  is *named* `ClosestApproachInRPhi`), and the bachelor case maps
  exactly to `AnalyticalImpactPointExtrapolator`'s docstring
  ("extrapolate to impact point with respect to vtx, i.e. to the
  point of closest approach to vtx in 3D").

#### Implementation plan

Two new private helpers in `JpsiXCandidateProducer.cc`, both thin
wrappers around the CMSSW classes named in "Why rough propagation
cures it" above:

1. **`propagatePair(const reco::TrackRef& a, const reco::TrackRef& b,
   double mA, double mB, const TransientTrackBuilder& ttb,
   reco::Particle::LorentzVector& lvA_out,
   reco::Particle::LorentzVector& lvB_out,
   GlobalPoint& crossing_out) const`**.
   Builds two `FreeTrajectoryState`s from the two `TrackRef`s via the
   `TransientTrackBuilder`. Constructs a `ClosestApproachInRPhi`,
   calls `calculate(fts_a, fts_b)`. On `status() == false`, returns
   `false` (caller falls back to raw daughter p4s, increments a
   counter). On success, reads `trajectoryParameters()`, builds
   `lvA_out` and `lvB_out` as `(px, py, pz, sqrt(p² + m²))` from each
   `GlobalTrajectoryParameters::momentum()`, and writes
   `crossingPoint()` into `crossing_out`.

2. **`propagateSingleToPoint(const reco::TrackRef& t, double m,
   const GlobalPoint& target, const TransientTrackBuilder& ttb,
   const MagneticField& field,
   reco::Particle::LorentzVector& lv_out) const`**.
   Builds a `FreeTrajectoryState` from `t`. Constructs an
   `AnalyticalImpactPointExtrapolator(&field)`, calls
   `extrapolate(fts, target)`. On `!result.isValid()`, returns false
   (caller falls back, increments a counter). On success, builds
   `lv_out` from `result.globalMomentum()` and the species mass `m`.

Call sites:

- In **`produceTrackMode`** (lines 519–622, B+ and Bc): replace the
  current `lvJpsi = jpsi.p4()` block + bachelor-4-vector construction
  at lines 572–574 with:
  1. `propagatePair(muRefs[0], muRefs[1], mμ, mμ, ttb, lvMuPlus, lvMuMinus, vJpsi)` → propagated muons + dimuon crossing point.
  2. `lvJpsi_prop = lvMuPlus + lvMuMinus`.
  3. `propagateSingleToPoint(bachelorRef, mBach, vJpsi, ttb, field, lvBach_prop)`.
  4. `lvM = lvJpsi_prop + lvBach_prop`. Use this `lvM` for the
     mass-window cut, the mother-pT cut, and the stored mother
     candidate. The preset-C multi-track Kalman path at lines 587–608
     subsequently overwrites `lvM` and `vtxM` from its own fit.

- In **`produceVccMode`** (lines 624–700, K*0/φ/Ks/Λ/ψ(2S)):
  1. `propagatePair(muRefs[0], muRefs[1], mμ, mμ, ttb, lvMuPlus, lvMuMinus, vJpsi)`.
  2. `lvJpsi_prop = lvMuPlus + lvMuMinus`.
  3. `propagatePair(xRefs[0], xRefs[1], mD1, mD2, ttb, lvD1, lvD2, vX)`
     where `mD1`/`mD2` are the sub-resonance daughters' PDG-mass
     hypotheses read from each daughter's stored
     `RecoChargedCandidate::mass()` (set upstream by
     `TwoBodyDecayCandidateProducer` / `V0Producer`).
  4. `lvX_prop = lvD1 + lvD2`.
  5. `lvM = lvJpsi_prop + lvX_prop`.

- The Kalman path's own per-track propagation under preset C is
  unchanged (it uses `KinematicParticleFactoryFromTransientTrack`
  internally, which calls the same primitives — but in fit context).

ESConsumes changes:

- `TransientTrackBuilder` is already consumed under preset C and the
  J/ψ-constraint branch. After this change it is consumed
  unconditionally — one-line edit to the constructor.
- `IdealMagneticField` (or rather the `MagneticField` ES record) is
  needed for `AnalyticalImpactPointExtrapolator`. Add the
  `esConsumes<MagneticField, IdealMagneticFieldRecord>()` token (the
  TransientTrackBuilder ESData already exposes a `field()` method, so
  in practice we can reuse the field already retrieved through the
  builder — single-line addition, no extra token).

#### Expected outcome on the 100k Condor sample

- m(B+) peak position shifts by ~1–2 MeV toward the PDG value
  (5.27934 GeV) relative to the preset-B-without-propagation run.
- m(B+) RMS unchanged within statistics on the 100k sample (the
  ~1 MeV scale is far below the 30 MeV detector resolution; visible
  only if we slice in B+ cτ).
- Per-channel candidate counts unchanged (efficiency impact zero).

### Decision 3: V0 mass — window confirmation, no code change

The V0 path uses windows end-to-end already:
- `V0Producer` (cloned as `ALCARECOTkAlV0Candidates`,
  `ALCARECOTkAlV0Candidates_cff.py`) applies its standard
  `kShortMassCut` / `lambdaMassCut` mass-window cuts during Kalman
  vertex fitting. No `MassKinematicConstraint` is built into V0Producer.
- `JpsiXCandidateProducer.produceVccMode` (lines 624–700) computes
  m(J/ψ + V0) as the sum of the J/ψ and V0 4-momenta, with the J/ψ
  constraint applied only to the dimuon side (this change disables that
  under preset B too). The V0 4-momentum is never constrained.

So after this change, the V0-mode channels are mass-window-only on
*both* the V0 sub-resonance and the J/ψ — exactly the user's
specification. No code change is required; the requirement is locked
by spec so it cannot regress.

### Decision 4: Production scale and parallelism

100 RAW files × 1k events = 100k events read, ~22k events written
(preset-B filter ~22% acceptance from prior 50k run). Submit as
100 Condor jobs of 1 file × 1000 events each, in parallel.

Why 100 × 1k and not 10 × 10k?

- **Wall-clock**: 100 jobs run in parallel finish in roughly the same
  wall-clock time as 1 job. The 5-file × 10k Condor batch (the
  previous production) took ~50 min per job. Each per-job RAW2DIGI +
  L1Reco + RECO + ALCA cost is dominated by the per-event RECO cost
  (~1.2 s/event), so 1k events ≈ 20 min, 10k events ≈ 3 hours.
  Splitting to 1k/job → ~20 min for all 100k events vs ~3 hours for
  10 × 10k. Sixfold wall-clock win.
- **Job retries**: smaller jobs are cheaper to re-run if any fail.
- **No state**: ALCA output is per-job-independent; the merge step
  collects the 100 ROOT files into one analysis tree.

The downside is more files to manage (100 instead of 5), but the
existing `merge_and_report.py` already iterates per-file so this is
trivial.

### Decision 5: One deck, not a patch on the existing one

Make a **new** `slides/jpsix-producer-final.tex` rather than editing
`jpsix-producer-progress.tex`. The progress deck has chronological
weight (it traces the bug-fix arc); the final deck should be
narrative-first and stand alone. Both stay in the repo.

## Risks / Trade-offs

- **Constraint removal degrades visible peak in raw plots**: m(B+)
  resolution goes from ~15 MeV to ~30 MeV. The peak is still
  comfortably above combinatorial in the 100k sample. Mitigation:
  call this out in the deck as a deliberate decision (Stage-2 CVH
  does the kinematic refit downstream; AlCaReco's job is raw data,
  not pre-fit data).
- **0.1 GeV track floor inflates file size**: estimated ~150 MB
  extra per 100k events. Mitigation: budget tracked in task 6.4;
  not expected to breach the AlCaReco bytes/event soft ceiling.
- **Condor parallel submission failures**: 100 jobs > 5 jobs;
  more opportunities for transient grid/storage issues. Mitigation:
  submit the first 5 as smoke (task 5.3), inspect, then release the
  remaining 95.

## Migration Plan

1. Land the cfg-parameter + cff-wiring + ptMin change behind a single
   `scram b`.
2. Smoke test on 1 file locally. If anything fails, the per-event
   producer-summary log (`n_jpsi_constraint_attempted_=0` under
   preset B) is the immediate diagnostic.
3. Smoke 5 Condor jobs from the 100-job batch. If clean, release the
   rest.
4. Merge + plot + draft the deck.

Rollback: revert the two file edits in
`JpsiXCandidateProducer.cc` + `ALCARECOTkAlJpsiX_cff.py`. The
Condor outputs are in a separate `preset_B_final/` directory and do
not collide with the existing `preset_B/` directory from the prior
production.

## Open Questions

None blocking. The user has permission-granted the optional dimuon
mass-window tightening ("can be relatively tight") and we are
declining it for this change — if combinatorial in the 100k sample
turns out worse than expected, we can tighten in a follow-up rather
than couple it to this change.
