# Change: Add J/ψ + X AlCaReco channels for B-meson and quarkonium alignment

## Why

The existing `TkAlJpsiMuMu` AlCaReco stream saves only the two muon tracks from J/ψ → μμ candidates. For track-based alignment using B-meson and quarkonium decays to J/ψ, the additional bachelor tracks (kaon, pion, proton) and intermediate resonance daughters (K*0, φ, Ks, Λ) are equally important — they probe different detector regions at different momenta and constrain curvature biases complementary to the muon-only signal. Events satisfying multiple channels in the same event should save their tracks only once.

## What Changes

### New/Modified C++ Plugins

- **Extend `TwoBodyDecayCandidateProducer`** with three backward-compatible optional parameters: `firstDaughterMass` (default: value of `daughterMass`), `secondDaughterMass` (default: value of `daughterMass`), and `tryBothChargeAssignments` (default: `False`). When `tryBothChargeAssignments = True`, each opposite-sign track pair is tried under both (positive=first mass, negative=second mass) and the charge-conjugate assignment; all candidates inside the mass window are emitted. All existing configs that pass only `daughterMass` are unaffected.
- **`JpsiXCandidateProducer`** (new) — combines a J/ψ `VertexCompositeCandidateCollection` with either:
  - a track from `generalTracks` under a configurable `bachelorMass` hypothesis, with internal `minBachelorPt` cut (B+, Bc channels), or
  - a second `VertexCompositeCandidateCollection` representing an intermediate resonance (all other channels).
  Applies a configurable mother mass window. Output: nested `VertexCompositeCandidateCollection` where daughter(0) = J/ψ (nested VCC) and daughter(1) = X (RecoChargedCandidate or nested VCC). The bachelor `TrackRef` points directly into `generalTracks` — no intermediate cloned collection — preserving the `TrackExtraRef` chain required by the downstream remapper.
- **`CompositeDaughterTrackProducer`** (new) — takes one or more `VertexCompositeCandidateCollection` inputs (each may be nested), recursively extracts all leaf `RecoChargedCandidate` tracks, deduplicates by `TrackExtraRef.key()`, and produces one merged `TrackCollection`. Each output `Track` is a shallow copy preserving its `TrackExtraRef` pointing into `generalTracks`.

### New Python Configuration

- **`ALCARECOTkAlJpsiX_cff.py`** — single combined sequence:
  - Stage 0: DCS filter, good-muon selector. No HLT filter in the default sequence; `ALCARECOTkAlJpsiXHLT` is defined in the cff but excluded from `seqALCARECOTkAlJpsiX` and can be prepended for online/express use.
  - Stage 1 (shared candidates): J/ψ candidates, V0 candidates (Ks, Λ), K*0 candidates, φ candidates.
  - Stage 2 (per-channel B-meson candidates): one `JpsiXCandidateProducer` instance per channel (7 total).
  - Stage 3 (merged tracks): one `CompositeDaughterTrackProducer` consuming all 7 B-meson VCC collections → one merged, deduplicated `TrackCollection`.
  - Stage 4 (selection + dE/dx): one `AlignmentTrackSelectorWithIndexMapModule` on the merged collection + three `DeDxValueMapProjector` instances (strip Harmonic2, pixel Harmonic2, joint). dE/dx is projected onto all tracks in the merged collection; it is physically meaningful for kaon-bearing channels (B+, B0→K*0, Bs) and the proton-bearing channel (Λb). For pion and muon tracks the values are saved but provide no useful PID.
  - Stage 5 (remapping): one `VertexCompositeCandidateRemapper` per channel (7 total), all sharing the same `selectedTracks` and `originalIndexMap`.
- **`ALCARECOTkAlJpsiX_Output_cff.py`** — `OutALCARECOTkAlJpsiX`: keep cloned tracks, all 7 remapped candidate collections, 3 dE/dx ValueMaps, trigger results, DCS status, primary vertices.

### Modified CMSSW Configuration

- `Configuration/StandardSequences/python/AlCaRecoStreams_cff.py` — import new cff; add `pathALCARECOTkAlJpsiX`; add `ALCARECOStreamTkAlJpsiX` `FilteredStream`
- `Configuration/AlCa/python/autoAlca.py` — add `TkAlJpsiX` to `AlCaRecoMatrix` under `Charmonium` and `MuOniaParked`
- `Configuration/EventContent/python/AlCaRecoOutput_cff.py` — extend `ALCARECOEventContent` with `OutALCARECOTkAlJpsiX_noDrop.outputCommands`

## Channel Catalogue

| Module label suffix | Mother | X | X daughters | Mother window (GeV) | Disp. req. | dE/dx useful? |
|---|---|---|---|---|---|---|
| `BPlusToJpsiK` | B+/- | K+/- track | — | 5.0–5.5 | none | yes (kaon) |
| `B0ToJpsiKstar` | B0 | K*0(892) | K± π∓ | 5.0–5.5 | none | yes (kaon) |
| `B0ToJpsiKs` | B0 | Ks | π+π- | 5.0–5.5 | V0Producer Ks | no (pions) |
| `BsToJpsiPhi` | Bs | φ(1020) | K+K- | 5.2–5.6 | none | yes (kaon) |
| `LambdabToJpsiLambda` | Λb | Λ | p π | 5.3–6.0 | V0Producer Λ | yes (proton) |
| `Psi2SToJpsiKs` | ψ(2S) | Ks | π+π- | 3.5–3.9 | V0Producer Ks | no (pions) |
| `BcToJpsiPi` | Bc± | π+/- track | — | 5.9–6.6 | none | no (pion) |

Intermediate resonance windows:
- K*0(892): 0.75–1.05 GeV (tried under both K+π- and K-π+ assignment)
- φ(1020): 0.99–1.06 GeV (symmetric, existing `TwoBodyDecayCandidateProducer`)
- Ks, Λ: handled by `ALCARECOTkAlV0Candidates` (shared with `TkAlKsToPiPi`/`TkAlLambdaToProtonPi`)

Bachelor track pT cut: pT > 0.5 GeV applied **internally** by `JpsiXCandidateProducer` when iterating over `generalTracks` for B+ and Bc channels. No separate pre-filter module — any `AlignmentTrackSelector`-based pre-filter rewrites `TrackExtraRef`s, breaking the remapper's `extra().key()` lookup.

## Impact

- Affected specs: `alcareco-jpsi-x` (new)
- New/modified files in `Alignment/CommonAlignmentProducer/`:
  - `plugins/TwoBodyDecayCandidateProducer.cc` (modified: add `firstDaughterMass`, `secondDaughterMass`, `tryBothChargeAssignments`)
  - `plugins/JpsiXCandidateProducer.cc` (new)
  - `plugins/CompositeDaughterTrackProducer.cc` (new)
  - `python/ALCARECOTkAlJpsiX_cff.py` (new)
  - `python/ALCARECOTkAlJpsiX_Output_cff.py` (new)
- Modified CMSSW: `AlCaRecoStreams_cff.py`, `autoAlca.py`, `AlCaRecoOutput_cff.py`
- Existing `TkAlJpsiMuMu` stream: unchanged
