# Design: TkAlJpsiX Phase 3 — Kalman Vertex Fits in Stage-1 AlCaReco

## Post-Phase-3 Stage-1 Pipeline Architecture

```
RAW data (Run2016H Charmonium)
        │
        ▼
  RawToDigi + L1Reco + Reconstruction
  (generalTracks, muons, offlinePrimaryVertices, offlineBeamSpot)
        │
        ▼
┌─────────────────────────────────────────────────────────────────────┐
│  ALCARECOTkAlJpsiXJpsiCandidates  (TwoBodyDecayCandidateProducer)  │
│                                                                     │
│  input: generalTracks + globalMuons                                 │
│  cuts:  mass ∈ [2.95, 3.25] GeV, opposite charge                   │
│  NO vertex fit: global muon ID already suppresses combinatorics;    │
│  B-level KVF in JpsiXCandidateProducer subsumes J/ψ sub-vertex.    │
│  output: reco::VertexCompositeCandidateCollection (J/ψ)             │
└──────────────────────┬──────────────────────────────────────────────┘
                       │  J/ψ candidates
          ┌────────────┼──────────────────┐
          │            │                  │
          ▼            ▼                  ▼
┌─────────────────┐  ┌──────────────┐  ┌──────────────────────────────┐
│ ALCARECOTkAlJpsi│  │ ALCARECOTkAl │  │ ALCARECOTkAlJpsiXPhiCandida- │
│ XKstarCandidates│  │ JpsiXKsCands │  │ tes (TwoBodyDecayCandidate-   │
│ (TwoBodyDecay-  │  │ (V0Producer) │  │ Producer)                    │
│ CandidateProdu- │  │              │  │                              │
│ cer)            │  │              │  │                              │
│                 │  │              │  │                              │
│ input: general  │  │ input: gen-  │  │ input: generalTracks         │
│   Tracks        │  │   Tracks     │  │ cuts: mass ∈ [0.990,1.040]   │
│ cuts: mass ∈    │  │ cuts: mass ∈ │  │       maxDaughterEta < 2.4   │
│  [0.80,0.99]    │  │  [0.43,0.56] │  │ NEW:  KalmanVertexFitter     │
│  maxDaughterEta │  │  flight sig. │  │       → vtx prob ≥ 0.01      │
│  < 2.4          │  │  ≥ 15σ       │  │ output: VCC (φ)              │
│ tryBothCharge-  │  │  vtx prob,   │  └────────────┬─────────────────┘
│ Assignments     │  │  pointing    │               │ φ candidates
│ NEW: KalmanVtxF │  │ output: VCC  │               │
│  → vtx p ≥ 0.01 │  │  (Ks,Λ)     │               │
│ output: VCC(K*0)│  └──────┬───────┘               │
└────────┬────────┘         │ Ks,Λ candidates        │
         │ K*0 candidates   │                        │
         │                  │                        │
         ▼                  │                        │
┌─────────────────┐         │                        │
│ B0→K*0          │         │                        │
│ JpsiXCandidate  │         │                        │
│ Producer (VCC)  │         │                        │
│                 │         │                        │
│ input: J/ψ,K*0  │         │                        │
│ mass:[5.1,5.4]  │         │                        │
│ minJpsiPt ≥ 3   │         │                        │
│ minMotherPt ≥ 5 │         │                        │
│ NEW: KVF on all │         │                        │
│  4 leaf tracks  │         │                        │
│  (μ⁺,μ⁻,K,π)   │         │                        │
│  minBVtxProb    │         │                        │
│  ≥ 0.01         │         │                        │
│  minBLxy/σ ≥2.0 │         │                        │
│  maxMotherAlpha │         │                        │
│  BS ≤ 0.3       │         │                        │
│  (replaces      │         │                        │
│  maxJpsiAlphaBS)│         │                        │
│ output: VCC (B0)│         │                        │
└─────────────────┘         │                        │
                            ▼                        ▼
                   ┌─────────────────┐    ┌──────────────────┐
                   │ B0→Ks           │    │ Bs→φ             │
                   │ JpsiXCandidate  │    │ JpsiXCandidate   │
                   │ Producer (VCC)  │    │ Producer (VCC)   │
                   │                 │    │                  │
                   │ input: J/ψ,Ks   │    │ input: J/ψ,φ     │
                   │ mass:[5.1,5.4]  │    │ mass:[5.2,5.6]   │
                   │ minJpsiPt ≥ 3   │    │ minJpsiPt ≥ 3    │
                   │ NEW: KVF on     │    │ minMotherPt ≥ 5  │
                   │  4 leaf tracks  │    │ NEW: KVF on      │
                   │  (μ⁺,μ⁻,π⁺,π⁻) │    │  4 leaf tracks   │
                   │  maxMotherAlpha │    │  (μ⁺,μ⁻,K⁺,K⁻)  │
                   │  BS ≤ 0.3       │    │  minBVtxProb     │
                   │  (no minBVtx-   │    │  ≥ 0.01          │
                   │  Prob since V0  │    │  minBLxy/σ ≥ 2.0 │
                   │  already clean) │    │  maxMotherAlpha  │
                   │ output: VCC(B0) │    │  BS ≤ 0.3        │
                   └─────────────────┘    │ output: VCC(Bs)  │
                                          └──────────────────┘

    J/ψ + generalTracks (track mode)
          │
          ▼
┌─────────────────────────────────────────────────────────────┐
│ B+ / Bc  JpsiXCandidateProducer (track mode)                │
│                                                             │
│ input: J/ψ candidates + generalTracks                       │
│ cuts (B+): minBachelorPt ≥ 0.5, maxBachelorEta < 2.4       │
│            maxBachelorIPToJpsiVertex < 1.0 cm               │
│            mass ∈ [5.1, 5.4] GeV, minMotherPt ≥ 5          │
│ cuts (Bc): minBachelorPt ≥ 0.3, maxBachelorEta < 2.4       │
│            maxBachelorIPToJpsiVertex < 1.0 cm               │
│            mass ∈ [5.9, 6.6] GeV, minMotherPt ≥ 5          │
│ NEW: KVF on 3 leaf tracks (μ⁺, μ⁻, bachelor K/π)           │
│      minBVtxProb ≥ 0.01                                     │
│      minBLxy/σ ≥ 2.0                                        │
│      maxMotherAlphaBS ≤ 0.3  (replaces maxJpsiAlphaBS)      │
│ output: reco::VertexCompositeCandidateCollection             │
└─────────────────────────────────────────────────────────────┘

    J/ψ + Ks (track mode, psi(2S) -> J/psi Ks)
          │
          ▼
┌─────────────────────────────────────────────────────────────┐
│ psi(2S) / Lambda_b  JpsiXCandidateProducer (VCC mode)      │
│                                                             │
│ These channels already clean from V0 quality on Ks/Lambda.  │
│ Phase 3: add maxMotherAlphaBS ≤ 0.3 only.                   │
│ minBVtxProb and minBLxy/sigma: not applied (V0 sufficient). │
└─────────────────────────────────────────────────────────────┘

All 7 B-candidate collections
          │
          ▼
┌──────────────────────────────────────────┐
│  ALCARECOTkAlJpsiX                       │
│  (AlignmentTrackSelectorModule)          │
│                                          │
│  per-track quality: hits, chi2, pt, eta  │
│  reads from VertexCompositeCandidate     │
│  daughters across all 7 collections      │
│  deduplicates tracks (set of TrackRefs)  │
│  output: reco::TrackCollection           │
└──────────┬───────────────────────────────┘
           │
           ▼
  ┌─────────────────────────┐
  │  CompositeDaughter-     │
  │  TrackProducer          │
  │  (TrackRef remapper)    │
  └──────────┬──────────────┘
             │
             ▼
  ┌──────────────────────────────────────────────────────┐
  │  ALCARECOStreamTkAlJpsiX output (PoolOutputModule)   │
  │                                                      │
  │  saved objects:                                      │
  │  - reco::TrackCollection (deduplicated)              │
  │  - TrackExtras + RecHits                             │
  │  - SiPixelClusters + SiStripClusters                 │
  │  - reco::VertexCompositeCandidateCollection (×7)     │
  │  - offlinePrimaryVertices                            │
  │  - L1 trigger + TriggerResults + DcsStatus           │
  └──────────────────────────────────────────────────────┘
```

## What KalmanVertexFitter Does at Each Node

```
TwoBodyDecayCandidateProducer (J/psi, K*0, phi):
  TransientTrack(trPos) + TransientTrack(trNeg)
        │
        └──► KalmanVertexFitter::vertex({trPos, trNeg})
                  │
                  ├─ isValid()? → drop if not
                  ├─ normalisedChiSquared() → chi2/ndof
                  ├─ ChiSquaredProbability(chi2, ndof) ≥ minVtxProb → drop if not
                  └─ fitted vertex position → overwrites midpoint in VCC

JpsiXCandidateProducer (B-level, all channels):
  collect all leaf TrackRefs recursively from B candidate
  build TransientTracks from those refs
        │
        └──► KalmanVertexFitter::vertex({tr1, tr2, ..., trN})   N=3 or N=4
                  │
                  ├─ isValid()? → drop if not
                  ├─ ChiSquaredProbability(chi2, ndof) ≥ minBVtxProb → drop if not
                  ├─ Lxy(fitted B vertex, nearest offline PV) / sigma_Lxy
                  │    ≥ minBLxyOverSigma → drop if not
                  └─ alphaBS = angle(B p4, beamspot→fitted_B_vertex)
                       ≤ maxMotherAlphaBS → drop if not
```

## Key Design Decisions

### D1. Track source for KVF in JpsiXCandidateProducer

The leaf tracks collected by `collectLeafTrackKeys` currently traverse the
`VertexCompositeCandidate` tree stored in the event. In Phase 3 we need the
same tracks as `TransientTrack` objects for the KVF. Two options:

**Option A** (recommended): traverse the VCC tree to collect `reco::TrackRef`s,
then build `TransientTrack` for each via `TransientTrackBuilder::build(ref)`.
This is the same pattern as V0Producer and is guaranteed to use the same
magnetic field state as the original reconstruction.

**Option B**: use the `reco::Track` momentum/perigee stored in the daughter
`RecoChargedCandidate`. Cheaper but bypasses the builder and loses curvature
covariance from the track fitter — the KVF would give wrong uncertainties.

Decision: Option A. Extend `collectLeafTrackKeys` to return `TrackRef`s (not
just keys).

### D2. Which offline PV to use for Lxy

The `offlinePrimaryVertices` collection is sorted by sum-pT². Use `[0]` (the
hardest PV). Guard against empty collection with `isValid()` check; drop the
Lxy cut gracefully if PV is missing (or require at least one PV).

### D3. 2D vs. 3D Lxy significance

Use transverse (2D) Lxy: `sqrt((Bvtx.x-PV.x)^2 + (Bvtx.y-PV.y)^2)`.
This avoids PV z-resolution contaminating the significance for boosted B
mesons. Compute sigma from the KVF vertex error matrix (diagonal elements xx,
yy) propagated with PV uncertainty.

### D4. KVF non-convergence handling

`KalmanVertexFitter::vertex()` can return `!isValid()` or
`normalisedChiSquared() == 0`. Both cases: drop the candidate silently (count
in `LogDebug`). Do not fall back to midpoint — a failed fit means the tracks
are incompatible and the candidate is combinatorial noise.

### D5. B0→Ks and Lambda_b vertex topology

All 4 tracks (μ⁺, μ⁻, π⁺, π⁻ from Ks; or μ⁺, μ⁻, p, π from Λ) are fit to a
single common vertex. This ignores the known Ks/Λ displacement but is the
correct B-physics approximation: the Ks/Λ tracks originate at the B decay
vertex (the Ks/Λ flight distance is small compared to the B vertex resolution
in these events). The V0Producer cut already validates the Ks/Λ decay
geometry; the B-level KVF is an additional constraint.

Alternative (not implemented): use the Ks vertex position as a point constraint
(`VertexConstraint`) in the B fit. More physically correct but significantly
more complex and not standard in CMS B-physics workflows. Deferred to
future work.

### D6. minBVtxProb for V0-only channels (B0→Ks, Lambda_b, psi(2S))

These channels are already ~1 cand/event after V0 quality. Adding `minBVtxProb`
would be redundant and risks unnecessary signal loss from KVF non-convergence
on 4-track combinations where Ks/Λ tracks are displaced. Do not apply
`minBVtxProb` or `minBLxyOverSigma` to these channels in Phase 3.
`maxMotherAlphaBS` (the pointing angle correction from J/ψ-level to B-level)
is still applied.

### D7. Backward compatibility of TwoBodyDecayCandidateProducer

`applyVertexFit` defaults to `False`. All existing configs (TkAlJpsiMuMu,
TkAlUpsilonMuMu, TkAlZMuMu, and the TkAlJpsiX J/ψ producer) are unaffected.
Only the K*0 and φ sub-producers in TkAlJpsiX will set `applyVertexFit = True`.
The J/ψ sub-producer is explicitly left without vertex fit: global muon ID
already makes it clean, and the B-level KVF subsumes it.

## Open Design Questions

**Q1. K*0 vtx p-value threshold (0.01 or higher?)**
K*0 daughters originate at the displaced B decay vertex (~1 mm from PV),
so two tracks from the true K*0→Kπ decay will always form a good vertex.
The 0.01 threshold is very permissive. Data-driven tuning needed after Phase 3
is deployed. Candidate thresholds: 0.01 (loose), 0.05 (moderate), 0.10
(tight). Check signal efficiency vs. background rejection in B0→K*0 mass plot.

**Q2. B Lxy/σ threshold (2.0σ or higher?)**
At pT > 5 GeV, B flight distance is ~1–5 mm; PV transverse resolution ~10–20
μm. Expected significance is O(50–500). A cut at 2σ is extremely loose and
essentially passes everything. Consider 5–10σ once peaks are visible.
The initial 2σ value is intentionally conservative to avoid signal loss while
tuning from data.

**Q3. psi(2S)→J/psi Ks treatment**
The psi(2S) is prompt (cτ ~ fm). Should `minBVtxProb` be applied here (as for
B+/B0/Bs/Bc) even though the Ks sub-decay provides some V0 quality?
Currently excluded from B-level vertex p-value cut (same as B0→Ks, Lambda_b).
Revisit once statistics allow comparison of psi(2S) peak with and without cut.

**Q4. J/psi vertex fit and alphaBS in JpsiXCandidateProducer**
After Phase 3, the J/ψ VCC already has a Kalman-fitted vertex position
(from TwoBodyDecayCandidateProducer with `applyVertexFit=True`). The
`maxJpsiAlphaBS` cut in JpsiXCandidateProducer is replaced by
`maxMotherAlphaBS` at the B level. However, should the J/ψ-level alphaBS cut
be kept as an additional pre-filter (using the now-corrected J/ψ KVF vertex)?
For B→J/ψX, the J/ψ vertex is essentially the B decay vertex (J/ψ cτ ~ 0),
so the J/ψ alphaBS and B alphaBS are the same observable. Keep only the
B-level cut to avoid double-counting.

**Q5. TransientTrackBuilder in EDStream vs ESConsumes**
`edm::stream::EDProducer` can use `ESConsumes` via
`esConsumes<TransientTrackBuilder, TransientTrackRecord>()` in the constructor.
The TransientTrackBuilder record is available whenever the full event setup is
loaded (full RECO, not FASTSIM). Confirm this is available in the
`TkAlJpsiX` ALCARECO re-running context (it should be — we load full RECO).

**Q6. CPU budget**
After Phase 3 the expected hot path is:
- K*0 sub-candidates: O(many) per event before vtx cut → cut reduces to O(1–5)
- K*0 is the dominant CPU consumer; KVF on 2 tracks at ~10–50 μs each
- B0→K*0 B-level candidates: fall from 104 to O(1–5) after K*0 vtx cut
- B-level KVF: O(10–100 μs) × O(1–5) candidates = cheap after sub-cut
Profile after first Phase 3 test run to confirm K*0 vtx p-value cut is
doing the expected work.
