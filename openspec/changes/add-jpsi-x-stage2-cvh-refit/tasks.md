# Tasks: TkAlJpsiX Phase 3 — Kalman Vertex Fits

Prerequisite: proposal.md and design.md approved.

---

## 1. TwoBodyDecayCandidateProducer.cc — optional sub-resonance vertex fit

- [ ] **1.1** Add includes: `TrackingTools/TransientTrack/interface/TransientTrackBuilder.h`, `TrackingTools/Records/interface/TransientTrackRecord.h`, `RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h`, `RecoVertex/VertexPrimitives/interface/TransientVertex.h`, `TMath.h` (for `TMath::Prob`)
- [ ] **1.2** Change base class from `edm::stream::EDProducer<>` to `edm::stream::EDProducer<>` with ESConsumes support — add `esConsumes<TransientTrackBuilder, TransientTrackRecord>()` token in constructor
- [ ] **1.3** Add parameters: `applyVertexFit` (bool, default false), `minVtxProb` (double, default 0.0) — read in constructor with `existsAs` guard
- [ ] **1.4** Add member: `edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttbToken_`; `bool applyVertexFit_`; `double minVtxProb_`
- [ ] **1.5** In `produce()`: when `applyVertexFit_` is true, get `TransientTrackBuilder` from EventSetup; after mass window passes, build `TransientTrack` for each daughter and run `KalmanVertexFitter::vertex({ttPos, ttNeg})`; drop if `!vtx.isValid()` or `ChiSquaredProbability < minVtxProb_`; update VCC vertex position to `vtx.position()`
- [ ] **1.6** Verify no change to output type or TrackRef ownership — the TrackRefs in daughter `RecoChargedCandidate` objects must still point into the original track collection (vertex fit does not rebind them)
- [ ] **1.7** Add `LogDebug` counter for candidates dropped by vertex fit

## 2. JpsiXCandidateProducer.cc — B-level Kalman vertex fit

- [ ] **2.1** Add same includes as 1.1 plus `DataFormats/VertexReco/interface/Vertex.h`, `DataFormats/VertexReco/interface/VertexFwd.h`
- [ ] **2.2** Promote base class to carry ESConsumes; add `esConsumes<TransientTrackBuilder, TransientTrackRecord>()` token
- [ ] **2.3** Add parameters: `minBVtxProb` (double, default 0.0), `minBLxyOverSigma` (double, default 0.0), `maxMotherAlphaBS` (double, default +inf) — read in constructor
- [ ] **2.4** Keep `maxJpsiAlphaBS` as deprecated backward-compatible parameter (read it if present, warn in constructor via `edm::LogWarning`)
- [ ] **2.5** Add token for `offlinePrimaryVertices` (consume only when `minBLxyOverSigma > 0`)
- [ ] **2.6** Add member variables for all new parameters and tokens
- [ ] **2.7** Extend `collectLeafTrackKeys` (or add a parallel `collectLeafTrackRefs`) to return `std::vector<reco::TrackRef>` instead of (or in addition to) the key set
- [ ] **2.8** In `produceTrackMode()`: after kinematic pre-selection passes, collect 3 TrackRefs (2 muons from J/ψ daughters + 1 bachelor), build TransientTracks, run KVF; apply `minBVtxProb`, `minBLxyOverSigma`, `maxMotherAlphaBS` cuts; update B VCC vertex to KVF position if all cuts pass
- [ ] **2.9** In `produceVccMode()`: after kinematic pre-selection passes, collect 4 TrackRefs (2 muons from J/ψ + 2 daughters from X VCC), build TransientTracks, run KVF; apply same cuts; update VCC vertex position
- [ ] **2.10** For `computeAlphaBS` at B level: use the KVF-fitted B vertex position (not midpoint) and the B candidate 4-vector — update the helper signature or add a new overload that takes a `GlobalPoint`
- [ ] **2.11** Lxy significance computation: `Lxy = sqrt((Bx-PVx)^2 + (By-PVy)^2)`; sigma from KVF vertex error (xx, yy) + PV error (xx, yy); guard against divide-by-zero
- [ ] **2.12** Use `offlinePrimaryVertices[0]` as the reference PV; guard against empty collection

## 3. ALCARECOTkAlJpsiX_cff.py — configuration updates

- [ ] **3.1** `ALCARECOTkAlJpsiXJpsiCandidates`: no change (global muon ID already clean; B-level KVF subsumes it)
- [ ] **3.2** `ALCARECOTkAlJpsiXKstarCandidates`: add `applyVertexFit = cms.bool(True)`, `minVtxProb = cms.double(0.01)`
- [ ] **3.3** `ALCARECOTkAlJpsiXPhiCandidates`: add `applyVertexFit = cms.bool(True)`, `minVtxProb = cms.double(0.01)`
- [ ] **3.4** `ALCARECOTkAlJpsiXBPlusResonances`: remove `maxJpsiAlphaBS`, add `maxMotherAlphaBS = cms.double(0.3)`, `minBVtxProb = cms.double(0.01)`, `minBLxyOverSigma = cms.double(2.0)`
- [ ] **3.5** `ALCARECOTkAlJpsiXB0KstarResonances`: same as 3.4
- [ ] **3.6** `ALCARECOTkAlJpsiXBsPhiResonances`: same as 3.4
- [ ] **3.7** `ALCARECOTkAlJpsiXBcResonances`: same as 3.4
- [ ] **3.8** `ALCARECOTkAlJpsiXB0KsResonances`: remove `maxJpsiAlphaBS`, add `maxMotherAlphaBS = cms.double(0.3)` only (no minBVtxProb, no minBLxyOverSigma — V0 quality sufficient)
- [ ] **3.9** `ALCARECOTkAlJpsiXLambdabResonances`: same as 3.8
- [ ] **3.10** `ALCARECOTkAlJpsiXPsi2SResonances`: same as 3.8

## 4. Build and validation

- [ ] **4.1** `scram b -j8` — confirm clean build, no new warnings
- [ ] **4.2** Run on 500 events from Run2016H Charmonium; check `LogDebug` output for candidates dropped at each vertex fit stage
- [ ] **4.3** Check per-channel candidate multiplicity vs. Phase 1+2 baseline (expect B0→K*0 to drop from 104 to O(1–5))
- [ ] **4.4** Confirm B+, B0→K*0, Bs→φ, Bc mass distributions show peaks above combinatorial background
- [ ] **4.5** Check B0→Ks, Λb, ψ(2S) are unchanged (still ~1 cand/event)
- [ ] **4.6** Profile CPU per-event; confirm K*0 vtx p-value cut is the dominant reduction step

## 5. Slides update

- [ ] **5.1** Add Phase 3 slide to `jpsix-alcareco-producer.tex`: pipeline diagram, new parameters table, expected cand/event table
- [ ] **5.2** Replace alphaBS alert block with corrected description (B-level αBS using KVF vertex)
- [ ] **5.3** Add post-Phase-3 mass plots once validation data are available
