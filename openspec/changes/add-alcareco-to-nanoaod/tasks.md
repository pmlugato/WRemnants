# Tasks — AlCaReco -> CVH refit -> NanoAOD

## 0. Decisions

- [x] 0.1 Q1 cross-release read — **DONE**, 500-event smoke passed.
- [x] 0.2 Q2 schema — bespoke, named to the histmaker's raw contract;
  corrected columns additive.
- [x] 0.3 Q3 muon table — option A (`reco::Muon` typedef).
- [x] 0.4 Q4 scope — all five streams / all channels.
- [x] 0.5 Q5 — hits/clusters dropped.
- [x] 0.6 Refit mode — `useIdealGeometry=False`, no corrections
  applied (identity) for v1.
- [x] 0.7 Confirm the single-G4-master limitation is lifted —
  **CONFIRMED** (`e33c8d440ac`, shared master ESProducer on
  `CvhMasterRecord`; `CvhMasterThread.h` "shared across all CVH
  residual makers in a job").
- [x] 0.8 Q6 — **DECIDED**: emit all three cross-link index columns
  (daughter->Track, Track->Muon, Track->PV); cheap and they are what
  make the flat tables joinable.
- [x] 0.9 Q7 — **DECIDED**: persist L1 bits and `DcsStatus`
  (confirmed lightweight: 2 floats + a 25-bit mask).
- [x] 0.10 Q8 — **DECIDED**: delete the splitter; makers descend the
  nested VCC (VCC format is common across channels; splitters would
  need per-topology variants).

## 1. Collapse the two-job refit into one job

- [x] 1.1 `test/runCvhBplusJpsiK.py`: `mode=both` (now the default)
  schedules both makers in one `reconstruction_step`; either/or branch
  replaced by additive sequence building.
- [x] 1.2 Stale "cannot co-exist / single G4 master per process" text
  removed from the driver header, the `mode` help, the assert, and
  `runCvhBplusJpsiK.md`; replaced with the `e33c8d` ES-product
  explanation.
- [x] 1.3 **VALIDATED**: both makers ran in one process on one AlCaReco
  file; both end-of-job fit summaries printed, no G4 conflict/abort.
  Also fixed a latent driver bug -- `cvhMasterESProducer` was only
  created inside the `useScalarPot3D` branch although every maker
  unconditionally `esConsumes` the record; now created always.
- [x] 1.4 Confirm both makers emit corrected ValueMaps — **CONFIRMED**:
  two-track behind `produceValueMaps` (set it True); single-track
  emits `corPt/corEta/corPhi/corCharge/corDxy/corDz/edmval/
  nValidHits/nValidPixelHits/momCov` unconditionally
  (`ResidualGlobalCorrectionMakerG4e.cc:273-280`). No new C++ needed
  for emission.
- [ ] 1.5 In-job combiner producer: pair the subsystem ValueMaps per
  candidate (keyed directly to the candidate collection once 2.4
  lands), emit the corrected parent mass and a
  `joinFlags`-equivalent orphan column.
- [ ] 1.6 Retire `scripts/btojpsik/join_cvh_bplus_jpsik.py` from the
  production path.

## 2. Refit configuration

- [ ] 2.1 `useIdealGeometry = False`; no correction file (identity).
- [ ] 2.2 `fillGrads = False`, `fillRunTree = False`,
  `fillTrackTree = False` (nano replaces the sidecar tree).
- [x] 2.3 **DONE + VALIDATED**: two-track maker descends a nested VCC daughter
  (~10 lines at `...TwoTrackG4e.cc:1142`): when the
  `dynamic_cast<RecoChargedCandidate*>` fails, descend into the
  composite daughter. Added `leafTrack()` / `subsystemPair()` helpers +
  a `subsystemDaughter` config knob (-1 = auto/legacy, >=0 = descend
  that daughter, for channels with two composite daughters).
  Reading the nested `ALCARECOTkAlJpsiXBPlusResonances` directly is
  **bit-identical** to the legacy pre-split input: both
  `attempted=15 succeeded=15`, pixel `seen=110 sizeX1=29 demoted=29`.
- [ ] 2.4 Route leaf-daughter tracks to the single-track maker from
  the same candidate collection (no splitter input).
- [ ] 2.5 Key both makers' ValueMaps to the original candidate
  collection (retires `bCandIdx`).
- [ ] 2.6 Delete `JpsiKCandidateSplitter` + its `_cfi` and unit test
  once 2.3-2.5 are green; push the maker change upstream to David's
  branch rather than letting it diverge.

## 3. NanoAOD plugins

- [x] 3.1 **DONE + VALIDATED**: `SimpleMuonFlatTableProducer`
  (`SimpleFlatTableProducer<reco::Muon>`) added + registered. Smoke on
  `ALCARECOTkAlJpsiXLooseMuons`: nMuon=2, `Muon_pt=[11.80, 9.34]`
  (matches the two leading `Track_pt`), `isGlobal`/`isTracker` true,
  `nMatches=4`.
- [x] 3.2 **DONE + VALIDATED**: daughter-accessor gap closed by the
  concrete `SimpleVertexCompositeCandidateFlatTableProducer`
  (`SimpleFlatTableProducer<reco::VertexCompositeCandidate>`).
  `daughter(0).pt` / `daughter(0).pdgId` now resolve. On B+:
  `dau0pdg=443` (J/psi composite), `dau1pdg=+-321` (bachelor kaon,
  signed convention holding), `dau0pt` constant across the candidates
  of one J/psi while `dau1pt` varies -- the expected nesting.
- [ ] 3.3 Index-map helper producer(s) emitting `ValueMap<int>` for
  all three cross-links: daughter->Track, Track->Muon (invert the
  persisted `TrackToMuon` Association), Track->PV.
- [x] 3.4 `scram b` from inside the touched packages only -- both
  `Analysis/HitAnalyzer` and `PhysicsTools/NanoAOD` build green.

## 4. Tables and config

- [ ] 4.1 Per-channel candidate tables: raw VCC vars + corrected
  ValueMaps as `externalVariables`.
- [ ] 4.2 `Track` table + dE/dx + `originalIndex` externalVariables.
- [ ] 4.3 `PV` and `Muon` tables.
- [ ] 4.4 HLT trigger-flag table + L1 decision bits from
  `L1GlobalTriggerReadoutRecord`.
- [ ] 4.5 Event-level `DcsStatus` branches (`magnetCurrent`,
  `magnetTemperature`, `ready`).
- [ ] 4.6 Cross-link index columns on the candidate and Track tables.
- [ ] 4.7 Branch names matching the histmaker raw contract
  (`bkmm_kaon_*`, `mm_mu*`, `Muon_*`).
- [ ] 4.8 Top-level config + `NanoAODOutputModule` outputCommands.

## 5. Validation

- [ ] 5.1 One-file job on `TkAlJpsiX`: both makers in one process.
- [ ] 5.2 `uproot`: raw + corrected columns present, counts consistent.
- [ ] 5.3 In-job corrected `m(mu mu K)` == old offline-join output on
  the same file.
- [ ] 5.4 Repeat per channel and per stream (all five).
- [ ] 5.5 Resource check: peak RSS + runtime vs. the two-job baseline.
- [ ] 5.6 Run `btojpsik.py` over the produced nano via its existing
  AlCaReco path.

## 6. Slides — expand the NanoAOD-production deck

Source: the 2-slide background deck already built; extend in place via
the `mit-slides` skill so it supersedes it at the same path.

- [ ] 6.1 Slide: how ours differs — AlCaReco is not MINIAOD; content
  inventory; which tables can/cannot be filled; why `mz_dilepton` /
  `mw_*` cannot consume this nano.
- [ ] 6.2 Slide: what we needed to do — refit in the loop (application
  mode); two jobs -> one after the shared-G4-master fix; generic
  nested-VCC decomposition replacing the splitter; missing
  `reco::Muon` producer; cross-links as the flat-tree ref substitute.
- [ ] 6.3 Slide: what we did — one-job pipeline diagram, table
  inventory (raw + corrected), histmaker branch-name contract.
- [ ] 6.4 Slide: validation — cross-release smoke, raw-vs-corrected
  closure vs. the retired offline join, per-channel counts, resources.
- [ ] 6.5 Slide: status / next steps — identity corrections now, real
  2016 corrections later, batch production.
- [ ] 6.6 Rebuild the PDF and confirm it renders.

## 7. Follow-up (separate changes)

- [ ] 7.1 Wire real 2016 postVFP corrections once solved.
- [ ] 7.2 Batch production over the 38,708-file AlCaReco set.
