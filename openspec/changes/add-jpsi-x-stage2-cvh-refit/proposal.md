# Change: TkAlJpsiX Phase 3 — Kalman vertex fits in Stage-1 AlCaReco

## Why

Phase 1+2 validation (19,619 events, Run2016H Charmonium) shows flat mother-mass
distributions in all non-V0 channels (B+, B0→K*0, Bs→φ, Bc) regardless of
statistics. The V0 channels (B0→Ks, Λb, ψ(2S)) are clean at ~1 cand/event.

The difference is entirely vertex quality. The V0Producer uses a `KalmanVertexFitter`
on the daughter pair and applies vertex probability + flight significance. The
non-V0 channels currently use only kinematic cuts — no vertex fit of any kind.

This change adds Kalman vertex fits at Stage 1, following the V0Producer pattern:

1. **B-level Kalman vertex fit** inside `JpsiXCandidateProducer`: all leaf
   tracks of the B candidate fit to a common vertex, giving B vertex probability,
   B flight significance from PV, and a proper B-level αBS (replacing the
   current incorrect J/ψ-level αBS).

2. **Sub-resonance vertex fit** inside `TwoBodyDecayCandidateProducer`: optional
   `KalmanVertexFitter` on K*0 and φ daughter pairs, giving a vertex probability
   cut even for prompt resonances. No flight significance (cτ ~ fm), but the
   vertex p-value still discriminates real K*0/φ pairs (tracks from a common
   space point) from random combinatorial pairs.

## Why V0Producer does not extend to K*0 / φ

The V0Producer's primary discriminant at Stage 1 is flight significance Lxy/σ.
K*0(892) has Γ ≈ 50 MeV (cτ ~ 4 fm) and φ(1020) has Γ ≈ 4 MeV (cτ ~ 50 fm):
both are prompt. Lxy is identically zero — the flight significance cut would
reject everything.

The vertex probability from the Kalman fit (the χ²/ndof of the two tracks at
the fitted vertex point) does work for prompt decays. Adding it to
`TwoBodyDecayCandidateProducer` via an optional parameter achieves the same
discrimination as V0Producer provides for K*0/φ daughters, without the
inapplicable flight significance step.

## What changes

### 1. `JpsiXCandidateProducer.cc` — B-level Kalman vertex fit

Consume `TransientTrackBuilder` from EventSetup (same as V0Producer).

For each B candidate passing kinematic pre-selection, collect all leaf
`TrackRef`s recursively (the existing `collectLeafTrackKeys` already traverses
the tree; extend it to collect refs), build `TransientTrack`s, and run
`KalmanVertexFitter`. Apply cuts on the result:

| New parameter | Default | Applied to | Meaning |
|---|---|---|---|
| `minBVtxProb` | 0 | both modes | B Kalman vertex fit probability (χ²/ndof) |
| `minBLxyOverSigma` | 0 | both modes | B transverse flight significance from nearest offline PV |
| `maxMotherAlphaBS` | +∞ | both modes | αBS of B candidate using fitted B vertex (replaces `maxJpsiAlphaBS`) |

The B vertex position stored in the output `VertexCompositeCandidate` is
updated to the Kalman-fitted position (currently it is the midpoint of track
reference points).

`maxJpsiAlphaBS` remains in the C++ as a deprecated backward-compatible parameter
but is removed from the 7 channel configurations.

Consume `offlinePrimaryVertices` (already in the ALCARECO output) when
`minBLxyOverSigma` is set; consume `offlineBeamSpot` when `maxMotherAlphaBS`
is set (same pattern as the existing `applyAlphaBS_` flag).

### 2. `TwoBodyDecayCandidateProducer.cc` — sub-resonance vertex fit

Add an optional Kalman vertex fit on each candidate pair before emitting:

| New parameter | Default | Applied to | Meaning |
|---|---|---|---|
| `applyVertexFit` | False | both modes | Run KalmanVertexFitter on the two daughters |
| `minVtxProb` | 0 | when `applyVertexFit` | Minimum vertex fit probability |

When `applyVertexFit = True`, the candidate vertex position is updated from
the track-midpoint to the Kalman-fitted position. Candidates failing
`minVtxProb` are dropped before the mass window cut.

This is used for K*0 and φ at Stage 1. The J/ψ producer is explicitly left
without a vertex fit: global muon ID already makes it combinatorially clean,
and the B-level KVF in `JpsiXCandidateProducer` fits the same muon tracks,
making a J/ψ sub-vertex fit redundant.

### 3. `ALCARECOTkAlJpsiX_cff.py` — configuration updates

**J/ψ producer** (`ALCARECOTkAlJpsiXJpsiCandidates`):
No change. J/ψ → μ⁺μ⁻ with global muons is already combinatorially clean; the
B-level KVF in `JpsiXCandidateProducer` subsumes any J/ψ sub-vertex quality.
Adding a vertex fit here gives no gain and costs CPU.

**K*0 producer** (`ALCARECOTkAlJpsiXKstarCandidates`):
```python
applyVertexFit = cms.bool(True),
minVtxProb    = cms.double(0.01),
```

**φ producer** (`ALCARECOTkAlJpsiXPhiCandidates`):
```python
applyVertexFit = cms.bool(True),
minVtxProb    = cms.double(0.01),
```

**B+ and Bc** (track mode, replace `maxJpsiAlphaBS`):
```python
maxMotherAlphaBS      = cms.double(0.3),
minBVtxProb           = cms.double(0.01),
minBLxyOverSigma      = cms.double(2.0),
```

**B0→K*0, Bs→φ** (VCC mode, non-V0):
```python
maxMotherAlphaBS  = cms.double(0.3),
minBVtxProb       = cms.double(0.01),
minBLxyOverSigma  = cms.double(2.0),
```

**B0→Ks, Λb, ψ(2S)** (V0 mode, already clean):
```python
maxMotherAlphaBS  = cms.double(0.3),   ## replaces maxJpsiAlphaBS
## minBVtxProb and minBLxyOverSigma not applied (V0 quality sufficient)
```

## Implementation notes

- The B-level Kalman fit in track mode collects: 2 muon TrackRefs from J/ψ
  VCC daughters + 1 bachelor TrackRef = 3 tracks.
- In VCC mode the J/ψ daughters (2 muons) + sub-resonance daughters (2 tracks
  for K*0/φ/Ks/Λ) = 4 tracks. All 4 are passed to the `KalmanVertexFitter`
  directly; no hierarchical sub-vertex used (K*0/φ are prompt; Ks/Λ flight
  already accounted for by the V0Producer cut, so re-fitting all 4 to a single
  B vertex is the standard B-physics approach).
- `KalmanVertexFitter` is already a dependency of `Alignment/CommonAlignmentProducer`
  (used transitively through other producers in the same package). The
  `TransientTrackBuilder` is available in the full RECO EventSetup.
- CPU cost: `KalmanVertexFitter` on 3–4 tracks is O(10–100 μs) per candidate.
  After Phase 1+2 kinematic cuts the candidate multiplicity is:
  B+ ~5/event, K*0 (sub) ~many before vtx cut, B0→K*0 ~104/event (before cut),
  Bs ~10/event, Bc ~9/event. The vertex fit on K*0 sub-candidates is the
  dominant cost; expected to reduce K*0 B-candidates from 104 to O(1–5)/event,
  making the per-B-candidate fit cheap.

## Expected impact

| Channel | Cands/event now | Expected after Phase 3 |
|---------|----------------|------------------------|
| B+ | 4.8 | O(1) |
| B0→K*0 | 104 | O(1–5) |
| B0→Ks | 1.0 | unchanged (V0 already clean) |
| Bs→φ | 10.3 | O(1) |
| Λb | 1.0 | unchanged |
| ψ(2S) | 1.0 | unchanged |
| Bc | 9.0 | O(1) |

Peaks should become visible once combinatorics are reduced to O(1) cand/event
and sufficient statistics are accumulated (O(1000) signal events per channel).

## Affected files

- `Alignment/CommonAlignmentProducer/plugins/JpsiXCandidateProducer.cc`
  (new: `TransientTrackBuilder` ESConsumes, Kalman vertex fit, new parameters)
- `Alignment/CommonAlignmentProducer/plugins/TwoBodyDecayCandidateProducer.cc`
  (new: optional `KalmanVertexFitter`, `applyVertexFit`, `minVtxProb`)
- `Alignment/CommonAlignmentProducer/python/ALCARECOTkAlJpsiX_cff.py`
  (new parameters for all 7 channels + K*0 and φ sub-producers; J/ψ producer unchanged)

## Open questions

1. **K*0 vertex p-value threshold**: K*0 daughters come from the B decay vertex
   (displaced ~mm from PV), not the PV itself. The vertex p-value cut on K*0
   daughters tests consistency of the two tracks with any common point in space.
   Is 0.01 the right threshold, or should it be higher to be more discriminating?
   Tune from data.
2. **B Lxy/σ threshold**: Using the Kalman-fitted B vertex (not the midpoint)
   gives better resolution. At pT > 5 GeV, B flight distance is ~1–5 mm; PV
   transverse resolution is ~10–20 μm → expected significance O(50–500). A cut
   of > 2 is very loose; > 5–10 may be more appropriate. Tune from data once
   peaks are visible.
3. **B0→Ks and Λb vertex fit topology**: should the Ks/Λ vertex (from
   V0Producer) be used as a point constraint in the 4-track B Kalman fit, or
   should all 4 tracks be fit independently to the B vertex? The latter ignores
   the known Ks displacement; the former is more physically correct but requires
   a `VertexConstraint`. Start with the simpler 4-track free fit.
4. **ψ(2S)→J/ψKs**: the ψ(2S) is also a prompt resonance (cτ ~ fm). Apply
   `minBVtxProb` here too even though V0 quality is already applied on the Ks
   sub-decay?
