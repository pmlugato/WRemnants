# Design: J/ψ + X AlCaReco channels (`TkAlJpsiX`)

## Context

The existing `TkAlJpsiMuMu` AlCaReco saves only J/ψ → μμ daughter tracks, using:
- `TwoBodyDecayCandidateProducer` (symmetric daughter mass, muon filter) → `VertexCompositeCandidateCollection`
- `V0DaughterTrackProducer` (flat daughter extraction)
- `AlignmentTrackSelectorWithIndexMapModule` (track cloning + `originalIndex` ValueMap)
- `VertexCompositeCandidateRemapper` (updates daughter `TrackRef`s; already recurse-capable for nested trees)

The existing `TkAlKsToPiPi` / `TkAlLambdaToProtonPi` use `ALCARECOTkAlV0Candidates` (a local `V0Producer` clone with looser pT cut, shared via deduplication).

## Goals / Non-Goals

**Goals**:
- Single ALCARECO stream `TkAlJpsiX` covering all 7 B-meson/quarkonium channels
- Save all leaf tracks from every reconstructed candidate tree (de-duplicated across channels in the same event)
- Re-key all candidate daughter `TrackRef`s onto the one cloned track collection
- Reuse shared upstream modules via CMS framework deduplication
- Follow the established AlCaReco pattern: DCS filter → candidates → track extraction → AlignmentTrackSelector → candidate remapper → output

**Non-Goals**:
- Vertex fitting / kinematic refit at any stage
- Primary-vertex IP significance cuts at B-meson level
- DQM monitoring modules (can be added later)
- MC truth matching

## Architecture Decision Summary

| # | Decision | Choice | Rationale |
|---|---|---|---|
| 1 | Stream topology | Single `TkAlJpsiX` | Events firing multiple channels save tracks once; follows user intent |
| 2 | K*0 assignment | Both ± assignments | K+π- and K-π+ both tried; all mass-window-passing candidates kept |
| 3 | V0Producer flight sig. | Keep defaults | Ks/Λ from B decays are *more* displaced, not less; cut easily satisfied |
| 4 | HLT filter | Not in default sequence | Targets reprocessing/calibration; HLT module defined but excluded from seq |
| 5 | Bachelor track pre-filter | pT > 0.5 GeV | Reduces J/ψ × track combinatorics; applied before `JpsiXCandidateProducer` |
| 6 | Bc channel | Include | Same infrastructure, only mass window differs |
| 7 | dE/dx output | Yes, all tracks in merged collection | dE/dx projects onto one merged TrackCollection — cannot be channel-selective; physically meaningful for kaon (B+, B0→K*0, Bs) and proton (Λb) channels; follows `TkAlKsToPiPi` pattern |
| 8 | Test dataset | `SingleMuon` first | Charmonium dataset to be tested later |
| 9 | HLT file clone | Not needed | Resolved by Decision 4 |

---

## Detailed Design

### Module Naming Scheme

All new module labels use the prefix `ALCARECOTkAlJpsiX` to avoid clashes with the existing `ALCARECOTkAlJpsiMuMu*` modules.

```
Shared Stage-1 modules:
  ALCARECOTkAlJpsiXGoodMuons          (TkAlMuonSelectors clone)
  ALCARECOTkAlJpsiXJpsiCandidates     (TwoBodyDecayCandidateProducer, dimuon)
  ALCARECOTkAlV0Candidates            (shared with TkAlKsToPiPi/Lambda; dedup by framework)
  ALCARECOTkAlJpsiXKstarCandidates    (TwoBodyDecayCandidateProducer, asymmetric mode)
  ALCARECOTkAlJpsiXPhiCandidates      (TwoBodyDecayCandidateProducer, sym.)
  [no separate bachelor pre-filter — pT>0.5 applied inside JpsiXCandidateProducer]

Per-channel Stage-2 modules:
  ALCARECOTkAlJpsiXBPlusCandidates
  ALCARECOTkAlJpsiXB0KstarCandidates
  ALCARECOTkAlJpsiXB0KsCandidates
  ALCARECOTkAlJpsiXBsPhiCandidates
  ALCARECOTkAlJpsiXLambdabCandidates
  ALCARECOTkAlJpsiXPsi2SCandidates
  ALCARECOTkAlJpsiXBcCandidates

Stage-3 (merged track extraction):
  ALCARECOTkAlJpsiXAllTracks          (CompositeDaughterTrackProducer, all 7 VCC inputs)

Stage-4 (selection + dE/dx):
  ALCARECOTkAlJpsiX                   (AlignmentTrackSelectorWithIndexMapModule)
  ALCARECOTkAlJpsiXDeDxHarmonic2
  ALCARECOTkAlJpsiXDeDxPixelHarmonic2
  ALCARECOTkAlJpsiXDeDxAllHarmonic2

Stage-5 (candidate remapping):
  ALCARECOTkAlJpsiXBPlusResonances
  ALCARECOTkAlJpsiXB0KstarResonances
  ALCARECOTkAlJpsiXB0KsResonances
  ALCARECOTkAlJpsiXBsPhiResonances
  ALCARECOTkAlJpsiXLambdabResonances
  ALCARECOTkAlJpsiXPsi2SResonances
  ALCARECOTkAlJpsiXBcResonances
```

---

### Decision 1: Single Combined Stream — TrackRef Chain

The key correctness question for a combined stream is whether the `VertexCompositeCandidateRemapper` can correctly remap candidate daughter `TrackRef`s when the merged intermediate track collection comes from multiple VCC inputs.

**TrackRef chain (full trace)**:

```
TwoBodyDecayCandidateProducer / JpsiXCandidateProducer:
  daughter.track() = TrackRef into generalTracks  (key K = generalTracks index)

CompositeDaughterTrackProducer (merged output):
  shallow-copy of generalTracks Track objects
  each merged track preserves TrackExtra → extra().key() = K (same generalTracks index)

AlignmentTrackSelectorWithIndexMapModule:
  src = merged collection (intermH)
  selH[j] → intermH[origIdx] → intermH[origIdx].extra().key() = K
  originalIndexMap: selH index j → intermH index origIdx

VertexCompositeCandidateRemapper (per channel):
  origToSel[K] = j  (built from selH + intermH + originalIndexMap)
  daughter.track().key() = K → origToSel[K] = j → newTrackRef(selH, j)
```

This chain is identical to the existing flat-V0 pattern. The only difference is that `intermH` is now a merged multi-VCC track collection instead of a single-VCC one. Because the shallow copy preserves `TrackExtraRef`, and `TrackExtraRef.key()` maps uniquely back to `generalTracks`, the remapper works correctly regardless of how many VCC collections contributed to the merged collection.

---

### Decision 2: K*0 Asymmetric Reconstruction — Extend `TwoBodyDecayCandidateProducer`

Rather than a new plugin, `TwoBodyDecayCandidateProducer` gains three backward-compatible optional parameters:

| New parameter | Default | Meaning |
|---|---|---|
| `firstDaughterMass` | value of `daughterMass` | mass assigned to positive track |
| `secondDaughterMass` | value of `daughterMass` | mass assigned to negative track |
| `tryBothChargeAssignments` | `False` | if True, try both +(first)/-(second) and +(second)/-(first) per pair |

All existing configs that pass only `daughterMass` see no behavior change. K*0 configuration:

```python
ALCARECOTkAlJpsiXKstarCandidates = TwoBodyDecayCandidateProducer.clone(
    src                      = cms.InputTag('generalTracks'),
    muonSrc                  = cms.InputTag(''),         # no muon filter
    firstDaughterMass        = cms.double(0.493677),     # K+
    secondDaughterMass       = cms.double(0.139570),     # π-
    firstDaughterPdgId       = cms.int32(321),
    secondDaughterPdgId      = cms.int32(211),
    motherPdgId              = cms.int32(313),           # K*0
    minMass                  = cms.double(0.75),
    maxMass                  = cms.double(1.05),
    tryBothChargeAssignments = cms.bool(True),
    applyChargeFilter        = cms.bool(True),
    charge                   = cms.int32(0),
    useUnsignedCharge        = cms.bool(True),
    applyAcoplanarityFilter  = cms.bool(False),
    acoplanarDistance        = cms.double(1.0),
)
```

φ → K+K- is symmetric (`daughterMass = 0.493677` only, `tryBothChargeAssignments = False`) so uses the existing interface unchanged.

---

### Decision 3: V0Producer Flight Significance

The default `generalV0Candidates` applies `vtxSignificance3D > 15` for Kshort. `ALCARECOTkAlV0Candidates` only modifies `tkPtCut = 0.1`. No change is needed.

For Ks from B decays: the B meson has cτ ≈ 500 μm → decay vertex typically 1–5 mm from PV. The Ks is produced at this B vertex, so the Ks flight distance (measured from PV) includes both the B flight and the Ks proper flight (cτ_Ks ≈ 27 mm). Flight significance is very large — the 15σ cut is trivially satisfied.

---

### Decision 4: HLT Filter — Not in Default Sequence

No HLT filtering by default. `ALCARECOTkAlJpsiXHLT` is defined in the cff (reusing `eventSetupPathsKey = 'TkAlJpsiMuMu'`, `throw = False`) but is NOT included in `seqALCARECOTkAlJpsiX`. Users running online/express production can prepend it manually. This differs from the existing `TkAlJpsiMuMu` stream, which does include the HLT filter; the difference is intentional — `TkAlJpsiX` targets reprocessing and calibration use cases where all J/ψ candidates are wanted regardless of trigger.

---

### Decision 5: Bachelor Track Pre-filter — Internal to `JpsiXCandidateProducer`

**Why not a separate pre-filter module**: Any `AlignmentTrackSelectorModule`-based producer uses `helper::TrackCollectionStoreManager` internally, which **rewrites** `TrackExtraRef`s in the cloned output. If a bachelor `TrackRef` pointed into such a cloned collection, `CompositeDaughterTrackProducer` would shallow-copy that track and find `extra().key()` pointing into the pre-filter collection's extras, not `generalTracks`. The `VertexCompositeCandidateRemapper` lookup (`origToSel.find(origRef.key())`) would fail because the candidate's `TrackRef.key()` is still the `generalTracks` index.

**Correct approach**: `JpsiXCandidateProducer` in track mode takes `src = generalTracks` directly (same as `TwoBodyDecayCandidateProducer`) and applies a configurable `minBachelorPt` cut internally when iterating over tracks. The bachelor `TrackRef` stored in the daughter `RecoChargedCandidate` points into `generalTracks`. The TrackRef chain is then identical to the dimuon/V0 cases.

```python
# In JpsiXCandidateProducer track-mode configuration:
ALCARECOTkAlJpsiXBPlusCandidates = cms.EDProducer('JpsiXCandidateProducer',
    xMode         = cms.string('track'),
    jpsiSrc       = cms.InputTag('ALCARECOTkAlJpsiXJpsiCandidates'),
    trackSrc      = cms.InputTag('generalTracks'),   # NOT a pre-filtered clone
    minBachelorPt = cms.double(0.5),                 # pT cut applied internally
    bachelorMass  = cms.double(0.493677),            # kaon
    bachelorPdgId = cms.int32(321),
    motherPdgId   = cms.int32(521),
    minMotherMass = cms.double(5.0),
    maxMotherMass = cms.double(5.5),
)
```

For VCC-mode channels (B0, Bs, Λb, ψ(2S)), the combinatorics are already limited by the small number of K*0/φ/Ks/Λ candidates, so no pT cut is needed in the producer.

---

### Decision 7: dE/dx Output

`DeDxValueMapProjector` operates on the single cloned track collection produced by `AlignmentTrackSelectorWithIndexMapModule`. Because all channels share one merged track collection (`ALCARECOTkAlJpsiXAllTracks`), dE/dx cannot be projected channel-selectively — it is computed for every leaf track regardless of which channel contributed it.

Physical motivation by channel:
- **B+, B0→K*0, Bs**: kaon daughter; dE/dx distinguishes K from π at momenta < 1 GeV/c
- **Λb**: proton daughter; dE/dx distinguishes p from K/π at momenta < 2 GeV/c
- **B0→Ks, ψ(2S)→Ks, Bc**: pion daughters only; no useful PID from dE/dx, but values are saved as a consequence of the merged-collection design
- **All channels**: muon daughters; dE/dx MIP — not used for PID

Three ValueMaps are produced following the `TkAlKsToPiPi` pattern exactly:
```python
ALCARECOTkAlJpsiXDeDxHarmonic2 = DeDxValueMapProjector.clone(
    selectedTracks     = 'ALCARECOTkAlJpsiX',
    intermediateTracks = 'ALCARECOTkAlJpsiXAllTracks',
    sourceTracks       = 'generalTracks',
    sourceValueMap     = 'dedxHarmonic2',
    originalIndexMap   = ('ALCARECOTkAlJpsiX', 'originalIndex'),
)
# plus dedxPixelHarmonic2 and alcaDedxJointEstimator clones
```

---

## Data Flow Diagram

```
generalTracks ─────────────────────────────────────────────────────────────────────────┐
GoodMuons ──► TwoBodyDecayCandidateProducer ──────────────────► JpsiCandidates         │
               (dimuon, 2.7–3.4 GeV)                                 │                 │
V0Producer clone ───────────────────────────────────► Ks, Λ cands    │                 │
TwoBodyDecayCandidateProducer (K*0: asymmetric mode) ► K*0 cands      │                 │
TwoBodyDecayCandidateProducer (φ:   symmetric mode) ─► φ cands        │                 │
                                                                      │                 │
                              ┌───────────────────────────────────────┤                 │
                              │   JpsiXCandidateProducer (×7)         │                 │
                              │   B+: Jpsi+generalTracks(K,pT>0.5)◄──────────────────────┤
                              │   B0→K*0: Jpsi+Kstar                  │                 │
                              │   B0→Ks: Jpsi+Ks                      │                 │
                              │   Bs→φ: Jpsi+Phi                      │                 │
                              │   Λb→Λ: Jpsi+Lambda                   │                 │
                              │   ψ(2S)→Ks: Jpsi+Ks                   │                 │
                              │   Bc: Jpsi+generalTracks(π,pT>0.5)◄──────────────────────┤
                              └──────────────┬────────────────────────┘                 │
                                             │ 7 B-meson VCC collections                │
                                             ▼                                           │
                              CompositeDaughterTrackProducer                             │
                              (merge all 7, recurse, dedup by extra().key())            │
                                             │                                           │
                                             ▼                                           │
                              AlignmentTrackSelectorWithIndexMapModule ◄────────────────┘
                              (src = merged, clones tracks+extras+hits+clusters,
                               emits originalIndex ValueMap)
                                             │
                                  ┌──────────┴──────────┐
                                  ▼                     ▼
                          DeDxValueMapProjector    VertexCompositeCandidateRemapper (×7)
                          (×3: strip, pixel, all)  (each uses same selectedTracks +
                                                    originalIndexMap; different srcCandidates)
                                                         │
                                                         ▼
                                                   ALCARECO output:
                                                   - cloned track collection
                                                   - 7 remapped B-meson collections
                                                   - 3 dE/dx ValueMaps
                                                   - trigger + DCS + PV
```

---

## Risks / Trade-offs

- **Bachelor-track combinatorics (B+, Bc)**: Even with pT > 0.5 GeV pre-filter, a J/ψ × (all tracks in dense PU event) loop is O(n_μμ × n_tracks). The B+ mass window [5.0, 5.5] GeV is very effective — most combinations will be outside. In practice n_jpsi is small (~1–5 per event) and the mass window suppresses ~99.9% of combinations.
- **Module count**: The single sequence has ~20 module instances. This is verbose but follows the established CMS pattern and is self-documenting.
- **Shared V0Candidates**: `ALCARECOTkAlV0Candidates` is shared with `TkAlKsToPiPi` and `TkAlLambdaToProtonPi`. If those streams run in the same job, the V0Producer runs once (CMS deduplication). If only `TkAlJpsiX` runs, it also runs once. No issue.
