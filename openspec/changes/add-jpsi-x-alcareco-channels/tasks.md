## 1. New C++ Plugins

- [ ] 1.1 Extend `TwoBodyDecayCandidateProducer.cc` — add `firstDaughterMass` (default: `daughterMass`), `secondDaughterMass` (default: `daughterMass`), and `tryBothChargeAssignments` (default: `False`) parameters; when `tryBothChargeAssignments = True`, try both charge assignments per opposite-sign pair and emit all mass-window-passing candidates; all existing configs with only `daughterMass` must produce identical output; update `fillDescriptions` accordingly
- [ ] 1.2 Implement `JpsiXCandidateProducer.cc` (track-bachelor mode, controlled by `xMode = "track"` parameter) — takes J/ψ VCC + `generalTracks` directly; applies internal `minBachelorPt` cut (NOT a separate pre-filter module, which would rewrite TrackExtraRefs); stores `TrackRef(generalTracksH, idx)` for the bachelor daughter so downstream chain (CompositeDaughterTrackProducer → extra().key() → remapper) is intact; parameters: `bachelorMass`, `bachelorPdgId`, `minBachelorPt`, mother mass window
- [ ] 1.3 Implement `JpsiXCandidateProducer.cc` (VCC-bachelor mode, `xMode = "vcc"`) — takes J/ψ VCC + X VCC + mother mass window; emits nested VCC where daughter(0)=J/ψ, daughter(1)=X sub-candidate; J/ψ four-momentum uses the stored candidate p4() (raw track-based, NOT PDG-mass-constrained — consistent with TwoBodyDecayCandidateProducer convention; mass window is wide enough to accommodate J/ψ resolution spread)
- [ ] 1.4 Implement `CompositeDaughterTrackProducer.cc` — accepts `srcs = cms.VInputTag(...)` (multiple VCC collections); recurses into nested composites; deduplicates by `TrackExtraRef.key()`; outputs one merged `TrackCollection` of shallow-copied tracks preserving `TrackExtraRef` pointing into `generalTracks`
- [ ] 1.5 Register both new plugins and the modified plugin in `Alignment/CommonAlignmentProducer/plugins/BuildFile.xml` (two new: `JpsiXCandidateProducer`, `CompositeDaughterTrackProducer`; one modified: `TwoBodyDecayCandidateProducer`)

## 2. Python Configuration: `ALCARECOTkAlJpsiX_cff.py`

- [ ] 2.1 Stage 0 — DCS + muon selector (no HLT filter in default sequence):
  - `ALCARECOTkAlJpsiXHLT` (hltHighLevel, `eventSetupPathsKey = 'TkAlJpsiMuMu'`) — defined but NOT added to `seqALCARECOTkAlJpsiX`; can be prepended by users who want online/express-style HLT gating
  - `ALCARECOTkAlJpsiXDCSFilter` (dcsstatus, same detector set as `TkAlJpsiMuMu`)
  - `ALCARECOTkAlJpsiXGoodMuons` (TkAlGoodIdMuonSelector clone)
- [ ] 2.2 Stage 1 — shared candidates:
  - `ALCARECOTkAlJpsiXJpsiCandidates` (TwoBodyDecayCandidateProducer: muonSrc=GoodMuons, 2.7–3.4 GeV, daughterMass=0.105, pdgId=443)
  - Import `ALCARECOTkAlV0Candidates` from `ALCARECOTkAlV0Candidates_cff`
  - `ALCARECOTkAlJpsiXKstarCandidates` (extended `TwoBodyDecayCandidateProducer` with `firstDaughterMass=0.494` (K), `secondDaughterMass=0.140` (π), `tryBothChargeAssignments=True`, mass window 0.75–1.05 GeV, pdgId=313)
  - `ALCARECOTkAlJpsiXPhiCandidates` (unmodified `TwoBodyDecayCandidateProducer`: sym K 0.494, 0.99–1.06 GeV, pdgId=333)
  - No separate bachelor pre-filter module — pT > 0.5 GeV cut applied internally by `JpsiXCandidateProducer` (see task 1.2; any `AlignmentTrackSelectorModule`-based pre-filter rewrites TrackExtraRefs, breaking the remapper)
- [ ] 2.3 Stage 2 — 7 × JpsiXCandidateProducer instances per channel table in proposal
- [ ] 2.4 Stage 3 — `ALCARECOTkAlJpsiXAllTracks` (CompositeDaughterTrackProducer, srcs=all 7 B-meson VCC labels)
- [ ] 2.5 Stage 4 — `ALCARECOTkAlJpsiX` (AlignmentTrackSelectorWithIndexMapModule, src=AllTracks, pT>0.4, filter=True) + 3 dE/dx projectors (`ALCARECOTkAlJpsiXDeDxHarmonic2`, `ALCARECOTkAlJpsiXDeDxPixelHarmonic2`, `ALCARECOTkAlJpsiXDeDxAllHarmonic2`) following TkAlKsToPiPi pattern; dE/dx is projected onto the single merged track collection (all leaf tracks from all channels); physically meaningful for kaon-bearing channels (B+, B0→K*0, Bs) and proton-bearing channel (Λb); pion and muon tracks also get values but these carry no useful PID information
- [ ] 2.6 Stage 5 — 7 × VertexCompositeCandidateRemapper instances; each uses:
  - `srcCandidates` = its channel's B-meson VCC
  - `selectedTracks` = `ALCARECOTkAlJpsiX`
  - `intermediateTracks` = `ALCARECOTkAlJpsiXAllTracks`
  - `originalIndexMap` = `('ALCARECOTkAlJpsiX', 'originalIndex')`
- [ ] 2.7 Assemble `seqALCARECOTkAlJpsiX` cms.Sequence in correct dependency order (Stage 0 → 1 → 2 → 3 → 4 → 5)

## 3. Output Configuration: `ALCARECOTkAlJpsiX_Output_cff.py`

- [ ] 3.1 Write `OutALCARECOTkAlJpsiX_noDrop` with `SelectEvents = 'pathALCARECOTkAlJpsiX'` and `outputCommands`:
  - `keep *_ALCARECOTkAlJpsiX_*_*` (cloned tracks + extras + hits + clusters)
  - `keep *_ALCARECOTkAlJpsiXDeDxHarmonic2_*_*` (strip, kaon-bearing channels; also present for other tracks)
  - `keep *_ALCARECOTkAlJpsiXDeDxPixelHarmonic2_*_*` (pixel)
  - `keep *_ALCARECOTkAlJpsiXDeDxAllHarmonic2_*_*` (joint alcaDedxJointEstimator)
  - `keep *_ALCARECOTkAlJpsiXBPlusResonances_*_*`
  - (similarly for the other 6 channel resonance collections)
  - `keep L1AcceptBunchCrossings_*_*_*`
  - `keep L1GlobalTriggerReadoutRecord_gtDigis_*_*`
  - `keep *_TriggerResults_*_*`
  - `keep DcsStatuss_scalersRawToDigi_*_*`
  - `keep *_offlinePrimaryVertices_*_*`
- [ ] 3.2 Create `OutALCARECOTkAlJpsiX` as `copy.deepcopy(OutALCARECOTkAlJpsiX_noDrop)` with `"drop *"` prepended

## 4. CMSSW Integration

- [ ] 4.1 `AlCaRecoStreams_cff.py`: import `ALCARECOTkAlJpsiX_cff`, add `pathALCARECOTkAlJpsiX = cms.Path(seqALCARECOTkAlJpsiX)`, add `ALCARECOStreamTkAlJpsiX` FilteredStream
- [ ] 4.2 `autoAlca.py`: add `'TkAlJpsiX'` to `Charmonium` and `MuOniaParked` entries in `AlCaRecoMatrix` and `AlCaRecoMatrix2017`
- [ ] 4.3 `AlCaRecoOutput_cff.py`: extend `ALCARECOEventContent.outputCommands` with `OutALCARECOTkAlJpsiX_noDrop.outputCommands`

## 5. Testing

- [ ] 5.1 Syntax-validate all new Python files with `python -m py_compile` inside CMSSW environment (scram build or cmsenv)
- [ ] 5.2 Write test config extending `recoskim_Run2016G_SingleMuon_AllResonances.py` to add `pathALCARECOTkAlJpsiX` and a `ALCARECOStreamTkAlJpsiX` output module; run on existing SingleMuon RAW input
- [ ] 5.3 Verify that `CompositeDaughterTrackProducer` extracts the correct number of unique tracks per candidate type (2 for J/ψ alone, 3 for B+/Bc, 4 for B0/Bs/ψ(2S)/Λb) via event printout or a small EDAnalyzer
- [ ] 5.4 Verify that all 7 `VertexCompositeCandidateRemapper` outputs are non-null and that remapped daughter `TrackRef`s point into the cloned `ALCARECOTkAlJpsiX` collection (not into `generalTracks`)
- [ ] 5.5 Check event counts: `TkAlJpsiX` events should be a subset of events that have at least one J/ψ candidate; B+ and ψ(2S) channels should dominate on Charmonium data

## 6. Documentation

- [ ] 6.1 File-level comment in `ALCARECOTkAlJpsiX_cff.py` documenting the 5-stage architecture and TrackRef chain
- [ ] 6.2 Doxygen class doc in all three new C++ plugin files (match style of `TwoBodyDecayCandidateProducer.cc`)
