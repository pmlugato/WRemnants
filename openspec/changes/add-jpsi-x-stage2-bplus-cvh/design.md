# Design: Stage-2 CVH refit for B+ → J/ψ K+ (reduced/composed)

## Context

The CVH framework in `CMSSW_15_0_19_patch2/src/Analysis/HitAnalyzer` already
provides everything needed to refit individual tracks and dimuon pairs with
Geant4e per-particle propagation under a configurable mass hypothesis, and
to emit the per-track Jacobians the GBL global-correction fit consumes. What
it does not provide is a multi-track (≥3) joint fitter: the two-track maker
reads a `VertexCompositeCandidate` and takes exactly `daughter(0)` and
`daughter(1)` as `RecoChargedCandidate` leaves
(`ResidualGlobalCorrectionMakerTwoTrackG4e.cc:878–885`).

B+ → J/ψ K+ has three leaf tracks. The decision (taken with the user) is to
avoid building an N-track fitter and instead **compose** the two existing
kernels: refit the μ⁺μ⁻ pair with the two-track maker, refit the bachelor K
with the single-track maker under the kaon hypothesis, and join their
outputs offline. A 3-track fitter is explicitly out of scope for this
program. The B+ vertex/mass constraint is not applied in Stage-2 and is
left for the downstream GBL fit (option a). The dimuon mass constraint
*inside* the two-track CVH refit is **on** for this change — see
"Mass-constraint summary" below for the distinction.

## Why composition is sufficient for the GBL fit

The GBL fit accumulates, per track, the Jacobian ∂(track params)/∂a and the
track's contribution to a χ² that includes the measurement residuals and
any external constraints. The B+ mass and B+ vertex are *constraints that
relate the three tracks*; they are not needed to compute each track's
Jacobian. Each track's Jacobian depends only on that track's hits, its mass
hypothesis (through energy loss), and the field/alignment model — all
available from the per-track refit. So:

- the two-track maker gives correct μ Jacobians (it already does, in
  production for J/ψ → μμ) — we additionally enable the dimuon mass
  constraint inside the joint fit (`doMassConstraint=True`) so the refit
  pulls m(μμ) to the PDG J/ψ mass;
- the single-track maker gives a correct K Jacobian *provided it propagates
  with the kaon mass* — which `trackParticleName="kaon"` does;
- the B+ mass/vertex constraint is then a Lagrange term added in the GBL
  fit, linking the three already-computed Jacobians.

This matches David's "with mass constraint" output variant for the dimuon
(slide 20), composed with an independent kaon refit for the bachelor.

### Mass-constraint summary (explicit, since easy to confuse)

The two-track CVH maker is *not* two independent single-track refits — it
fits a **common-vertex 10-dim PCA state** for the two muons. Common-vertex
coupling and a *mass* constraint are distinct knobs:

- **Common-vertex coupling** of the two muons in the two-track maker:
  always on (it's the joint state, not a toggle). The two muons share a
  vertex in the refit; their Jacobians come out *coupled* through that
  shared state.
- **`doMassConstraint` on the two-track maker** (the dimuon mass
  constraint inside the joint fit): **ON** for this change. The refit
  pulls m(μμ) to the PDG J/ψ mass inside the joint vertex fit. The
  per-track Jacobians delivered to the downstream GBL fit are the
  mass-constrained ones.
- **Kaon refit**: fully independent at its own PCA — no link to the
  dimuon vertex in Stage-2.
- **B+ vertex/mass constraint coupling μμ to K**: not applied in Stage-2;
  not applied per-event at the GBL stage either. GBL solves for global
  module corrections **a** from per-track Jacobians; a B-mass/B-vertex
  Lagrange term linking the three tracks could be added downstream but
  is out of scope for the "see a sane B+ peak after refit" milestone.

## Module boundary: splitter + two makers + offline join

Two designs were considered.

**(A) One combined driver module** that internally invokes both refit
kernels per B+ candidate and writes a single combined tree. Cleaner output
(no offline join) but duplicates the makers' substantial setup
(EventSetup tokens, G4 master-thread plumbing, dedx lookups) inside a third
module, and couples the B+-specific logic to the kernel internals.

**(B) Splitter + two existing makers + offline join** (chosen). A thin
`JpsiKCandidateSplitter` reads the B+ VCC collection and emits the dimuon
sub-candidates and the bachelor tracks as separate collections, each tagged
with the B+ candidate index. The two existing makers run unchanged on those
collections. A Python/ROOT script joins the two output trees on
`(run, lumi, event, bCandIdx)`.

(B) is chosen because it reuses the makers as-is (no fit code touched), keeps
the new C++ to a thin, fit-free splitter, and matches the user's "offline
combine" framing. The cost is an offline join step and a per-candidate index
that must travel through both makers.

```
TkAlJpsiX ALCARECO (2016, cross-release)
        │  B+ VertexCompositeCandidateCollection
        ▼
┌──────────────────────────────────────────────┐
│ JpsiKCandidateSplitter  (NEW, fit-free)        │
│  for each B+ cand i:                           │
│   • emit dimuon VCC (J/ψ daughter)  → tag i   │
│   • emit bachelor K track/RCC       → tag i   │
└───────────┬───────────────────┬────────────────┘
   dimuon VCCs            bachelor K tracks
            │                   │
            ▼                   ▼
┌────────────────────┐  ┌────────────────────────┐
│ TwoTrackG4e        │  │ G4e (single-track)      │
│  muon hypothesis    │  │  trackParticleName=kaon │
│  massConstraint=ON  │  │                         │
│  → tree: Jpsi_/     │  │  → tree: Kbach_* + jac  │
│    Mu±_* + jacRef   │  │    + bCandIdx           │
│    + bCandIdx       │  │                         │
└─────────┬──────────┘  └───────────┬─────────────┘
          └──────── offline join on ─┘
              (run, lumi, event, bCandIdx)
                        │
                        ▼
        one flat tree per B+ candidate
        (3 tracks' CVH params + cov + Jacobians,
         m(μμK), raw B+ kinematics)  →  GBL fit
```

## Per-candidate index

A robust 1:1 offline join needs an explicit B+ candidate index, because
either leg's refit can fail (the maker drops the row) and row order is not
guaranteed to survive. Two implementations:

- **Output branch**: each maker writes a `bCandIdx` int branch, fed from the
  split collection. Smallest change to the join logic; needs an additive
  branch in both makers.
- **`edm::ValueMap<int>`**: the splitter associates an index with each
  emitted element; the maker reads and writes it. Slightly more wiring.

Either is additive and no-op for the existing J/ψ/Υ/Z/Ks/Λ/D0 configs (no
index input → branch defaults to −1 / not written). Choice deferred to apply.

## Cross-release input

Stage-1 ALCARECO is the produced finalized preset-B Condor output in
`CMSSW_10_6_17_patch1` (Run2-2016, `preset_B_final/Run2016H/`); Stage-2
runs in `CMSSW_15_0_19_patch2`. This is already an established pattern:
`runCvhJpsi.py` reads a 2016 `ALCARECOTkAlJpsiMuMu` with `Run2_2016` /
`auto:run2_data` in this exact release. The only new wrinkle is that the B+
VCC is a *nested* composite (a VCC whose daughter is itself a VCC). Standard
`reco::VertexCompositeCandidate` is release-portable, so the risk is low; a
deserialization smoke test on the produced ALCARECO is a Stage-0 task.

**SplitLevel=1 prerequisite** (cross-cutting; affects Stage-0 readiness).
The current `preset_B_final` files were written without
`setAlcaRecoSplitLevel` and therefore have `SiStripCluster` /
`TrackingRecHit` / `recoTracks` / `recoTrackExtras` branches at split-level
99. Reading those branches in CMSSW_15_0_19 hits ROOT cms-sw#19773: the
SiStripCluster v11→v14 schema rule is not applied to split branches, and
clusters silently return empty when dereferenced through the rec-hit
chain. The CVH refit cannot reconstruct the bachelor's hits in that
state. The fix is David's commit `b8ebaa3d9f2` (adds
`Alignment/CommonAlignmentProducer/python/alcarecoSplitLevel.py` with
`setAlcaRecoSplitLevel`, wired into the cmsDriver `--customise` list);
that commit is **not on the current `update-btojpsik-fitting` branch**.
Cherry-pick + re-produce (or re-stream the 96 files through a one-pass
`cmsRun` cfg that re-writes splitLevel=1) before Stage-2 can run.

## Validation strategy

1. Read one finalized preset-B 100k file
   (`/ceph/submit/data/user/p/pmlugato/mz/alcareco/jpsi-x/preset_B_final/Run2016H/`)
   in `CMSSW_15_0_19_patch2` **with the splitLevel=1 prerequisite resolved**
   (cherry-pick or re-stream — see `tasks.md` 0.0). Confirm the B+ VCC
   deserializes; confirm a sample bachelor track's SiStripCluster bytes
   load non-empty when dereferenced through the I/O-rule path (NOT just via
   FWLite `dataSize()`, which bypasses the rule chain); confirm the
   splitter emits the expected dimuon + bachelor counts (per-event cand
   rate from `finalize-jpsi-x-preset-b-production` headline numbers).
2. Run the chain; confirm both refit trees populate and join 1:1.
3. Sanity-check refit vs. raw: per-track refit params close to the input
   track params; covariances finite; `Kbach_jacRef` non-trivial.
4. m(μμK) under the kaon hypothesis peaks at the B+ mass on the joined tree.
5. Jacobian closure: use the framework's existing `runFDClosure` /
   `epsilonFDClosure` finite-difference check on the kaon leg to confirm the
   per-mode chain-rule columns are correct for a kaon (not just a muon).

## Risks

- **Cross-release VCC** — low; mitigated by the Stage-0 smoke test.
- **Bachelor maker reads tracks, not candidates** — the splitter output type
  must match what the single-track maker consumes; confirm at apply.
- **Index plumbing** — additive, but touches two makers; keep it behind an
  absent-input no-op so existing configs are untouched.
- **B+ vertex/mass constraint not applied in Stage-2** — by design (option a);
  flagged for the GBL-fit owner so the B-mass/vertex Lagrange term is applied
  downstream, not forgotten. (The dimuon J/ψ mass constraint inside the
  two-track refit *is* applied at Stage-2; it's a separate knob.)
