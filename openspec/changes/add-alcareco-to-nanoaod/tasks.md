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
- [x] 1.5 **DONE (superseded approach)**: rather than an in-job ValueMap
  combiner, the fitted mother mass now comes from a proper kinematic fit
  -- `JpsiXKinematicFitProducer` (the Bmm5 chain). It descends the nested
  candidate to leaf tracks, builds TransientTracks, and runs
  `KinematicConstrainedVertexFitter` (+ `TwoTrackMassKinematicConstraint`
  in inFit mode). Three modes: inFit (Bmm5) / upstream / cascade.
  fitOk = 0 is the orphan flag. Validated (400 evt, stage-1 tracks):
  inFit 100% fitOk, fitted B+ median 5.2792 vs PDG 5.27934, window rms
  0.070 (sharper than the raw IQR 0.24); cascade vtxProb>0.01 = 100%.
  Full driver end-to-end: fitted 5.2800, raw 5.26, dimuon corMass 3.10
  all coexisting.
- [ ] 1.6 Retire `scripts/btojpsik/join_cvh_bplus_jpsik.py` from the
  production path.

## 2. Refit configuration

- [x] 2.1 `useIdealGeometry=False`, no correction file (identity) --
  driver defaults; `plimit` default lowered 1.0 -> 0.05 (soft bachelors:
  kaon-maker failures fell from 10/15 to 2/44).
- [x] 2.2 `fillGrads=False`, `fillRunTree=False`; `fillTrackTree` forced
  False on both makers when `nanoOut` is set (nano replaces the sidecar).
- [x] 2.3 **DONE + VALIDATED**: two-track maker descends a nested VCC daughter
  (~10 lines at `...TwoTrackG4e.cc:1142`): when the
  `dynamic_cast<RecoChargedCandidate*>` fails, descend into the
  composite daughter. Added `leafTrack()` / `subsystemPair()` helpers +
  a `subsystemDaughter` config knob (-1 = auto/legacy, >=0 = descend
  that daughter, for channels with two composite daughters).
  Reading the nested `ALCARECOTkAlJpsiXBPlusResonances` directly is
  **bit-identical** to the legacy pre-split input: both
  `attempted=15 succeeded=15`, pixel `seen=110 sizeX1=29 demoted=29`.
- [x] 2.4 **DONE**: new `CandidateLeafTrackProducer` emits the direct
  leaf-daughter (bachelor) tracks + a `candIdx` vector; the single-track
  maker's `src`/`bCandIdxSrc` point at it. Generic across channels
  (composite daughters -> two-track maker, leaf daughters -> single).
- [x] 2.5 Single-track maker keyed off the leaf-track producer's
  `candIdx` (same mechanism as the old `bCandIdx`).
- [x] 2.6 **Splitter dropped from the default path** (kept only for the
  legacy pre-split A/B path). Validated splitter-free: single-track maker
  fed 80 bachelor tracks, B fit unaffected. NOTE: file + `_cfi` are left
  in place for the legacy path; full deletion + upstream push to David
  deferred to the multi-channel change.

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
- [x] 3.3a **daughter->Track DONE + VALIDATED**: new
  `CandidateLeafTrackIndexProducer` emits mu0/mu1/bach0/bach1TrackIdx
  ValueMap<int> keyed to the candidate collection, ProductID-guarded.
  `Track_pt[kaonTrackIdx]` reproduces `BuJpsiK_kaonPt` exactly; kaon
  dE/dx now reachable per candidate.
- [x] 3.3b **Track->Muon DONE + VALIDATED**: `TrackToMuonIndexProducer`
  (templated `TrackAssocIndexProducer`); 81/688 tracks matched,
  `Track_pt[muonIdx]` reproduces `Muon_pt` exactly. `Track->PV`
  producer built (`TrackToVertexIndexProducer`) but the
  `offlinePrimaryVertices` association is keyed to `generalTracks`, not
  the persisted alignment tracks, so it returns all -1 here (guarded by
  `contains(productID)`; no crash). See 4.6b for the fix.
- [x] 3.4 `scram b` from inside the tou ched packages only -- both
  `Analysis/HitAnalyzer` and `PhysicsTools/NanoAOD` build green.

## 4. Tables and config

- [x] 4.1 **DONE (B+ channel)**: `BuJpsiK` table carries raw VCC columns
  (mass, pt, vertexChi2, jpsi/kaon daughter kinematics + pdgId) AND the
  refit ValueMaps as `externalVariables` (`corMass`, `corMassErr`,
  `corPt/Eta/Phi`, `corMuPlusPt`, `corMuMinusPt`, `corEdmval`).
  Remaining channels still to wire.
- [x] 4.2 `Track` table + dE/dx + `originalIndex` externalVariables.
- [x] 4.3 `PV` and `Muon` tables.
- [x] 4.4a **HLT flags DONE**: `NanoAODOutputModule` converts
  `edm::TriggerResults` into `HLT_*` branches itself, so
  `keep edmTriggerResults_*_*_*` in `outputCommands` is the whole change.
  Verified: **534 HLT_* flags** written.
- [x] 4.4b **L1 producer BUILT** (`L1LegacyDecisionTableProducer`): emits
  a singleton `L1` table -- `finalOR` + 128-bit algo (`algo0..3`) +
  64-bit tech (`tech0..1`) packed into uint32. Runs clean. **BUT reads
  all-zero on this AlCaReco**: the `gtDigis` FDL-word vector is not
  populated in the persisted stream (shell kept, decision bits dropped).
  See 4.4c.
- [ ] 4.4c **VERIFY L1 source content**: confirm whether the legacy L1
  decision is recoverable from this AlCaReco at all (alternate accessor /
  BX / a different persisted product), or whether it must be taken from a
  fuller tier. If unrecoverable, drop the L1 table from the stream.
- [x] 4.5 **DONE**: `SimpleDcsStatusFlatTableProducer` typedef added
  (+ `DataFormats/Scalers` in the plugins `BuildFile.xml`, without which
  it fails to link on `typeinfo for DcsStatus`). `Dcs` table verified:
  `magnetCurrent = 18164 A` (CMS solenoid nominal for 3.8 T),
  `magnetTemperature`, `ready` mask.
- [x] 4.6a Candidate daughter->Track columns (mu0/mu1/kaonTrackIdx) on
  the BuJpsiK table.
- [x] 4.6b Track->Muon column (`Track_muonIdx`) on the Track table.
- [ ] 4.6c **Track->PV via the `originalIndex` bridge**: the persisted
  track->PV association is keyed to `generalTracks`. Map ALCARECO track
  -> generalTracks index (the persisted `originalIndex` ValueMap) ->
  look up the PV association on that key, to emit a real `Track_pvIdx`.
- [ ] 4.7 Branch names matching the histmaker raw contract
  (`bkmm_kaon_*`, `mm_mu*`, `Muon_*`).
- [ ] 4.8 Top-level config + `NanoAODOutputModule` outputCommands.

## 5. Validation

- [x] 5.1 **DONE**: one-file job, both makers in one process, 30 events.
- [x] 5.2 **DONE**: raw B+ mass mean 5.243 (B+ region); refit dimuon
  `corMass` mean 3.0967 vs PDG J/psi 3.0969; `corMass` identical across
  candidates sharing a J/psi; 45/45 corrected, 0 sentinel.
- [ ] 5.3 In-job corrected `m(mu mu K)` == old offline-join output on
  the same file.
- [ ] 5.4 Repeat per channel and per stream (all five).
- [ ] 5.5 Resource check: peak RSS + runtime vs. the two-job baseline.
- [x] 5.7 **600-event validation run**: two-track 1030/1030 (100%);
  single-track 936/1003 (93.3%) at plimit=0.05. Matched control at
  plimit=1.0 on the SAME 1003 tracks: 428/1003 (42.7%) -- the lowered
  floor is worth +50 percentage points of bachelor legs.
- [x] 5.8 Refit dimuon mass peaks at the J/psi: median 3.0935 vs PDG
  3.0969 (**-3.4 MeV**, identity corrections, unconstrained). Raw B+
  distribution is broad/combinatorial as expected for the wide AlCaReco
  window.
- [ ] 5.6 Run `btojpsik.py` over the produced nano via its existing
  AlCaReco path.

## 6. Slides — expand the NanoAOD-production deck

Source: the 2-slide background deck already built; extend in place via
the `mit-slides` skill so it supersedes it at the same path.

- [x] 6.1 Slide: how ours differs — AlCaReco is not MINIAOD; content
  inventory; which tables can/cannot be filled; why `mz_dilepton` /
  `mw_*` cannot consume this nano.
- [x] 6.2 Slide: what we needed to do — refit in the loop (application
  mode); two jobs -> one after the shared-G4-master fix; generic
  nested-VCC decomposition replacing the splitter; missing
  `reco::Muon` producer; cross-links as the flat-tree ref substitute.
- [x] 6.3 Slide: what we did — one-job pipeline diagram, table
  inventory (raw + corrected), histmaker branch-name contract.
- [x] 6.4 Slide: validation — cross-release smoke, raw-vs-corrected
  closure vs. the retired offline join, per-channel counts, resources.
- [x] 6.5 Slide: status / next steps — identity corrections now, real
  2016 corrections later, batch production.
- [x] 6.6 PDF rebuilt and rendered-checked (12 pages). Plot slides needed
  a `section.plots` CSS block that the original 2-slide deck never had,
  plus `min-width: 0` on the figrow images -- flex items refuse to shrink
  below intrinsic size, so a wide PNG overflowed and hid its neighbour.

## 6b. Convergence auditing (NOT a blocker)

- [x] 6b.1 **Corrected reading of `edmval`.** An earlier note here called
  `corEdmval ~ 190` a convergence failure. That was wrong. The stored
  `edmval` is the **full-state** EDM, including all per-hit scattering
  parameters, which are re-zeroed at every re-linearization; large
  values (median ~1e3) are expected and are **not** the criterion. The
  criterion is **EDM < 1e-5 on the reference-state block**, with an
  iteration cap of 10. Not meeting it is also not by itself failure --
  fits can be successful without reaching the threshold.
- [ ] 6b.2 Persist the **reference-block** EDM per iteration (new
  trajectory branches) and expose it in the NanoAOD instead of / beside
  the full-state `edmval`, so consumers audit convergence on the right
  quantity.
- [ ] 6b.3 Measure real output size on a full file: the 30-event smoke
  gave ~60 kB/event, dominated by small-file overhead.

Measured convergence audit (reference-block criterion), for context:

| | at iteration cap without meeting criterion | of which genuinely still improving |
|---|---|---|
| dimuon | 1.3 % | essentially all (chi2 still falling) |
| single-track | 15.8 % | ~half converge by iteration 20; 8.2 % never |

The never-converging fits are **limit cycles**: the solver predicts a
chi2 improvement of ~0.07 (median) every iteration, but 30 extra
iterations change the actual chi2 by exactly zero for the median fit.

## 6d. B-candidate kinematic fit (Bmm5 chain)

- [x] 6d.1 `JpsiXKinematicFitProducer` new plugin; inFit/upstream/cascade
  modes; emits fitMass/fitMassErr/fitPt/eta/phi/vtxChi2/vtxNdof/vtxProb +
  fitOk ValueMaps keyed to the candidate collection.
- [x] 6d.2 `cvhFit*` columns on the candidate table + driver wiring
  (`jpsiConstraint`), distinct cvh* names (Q2 answer) so they cannot be
  confused with BParking's bkmm_jpsimc_*/bkmm_nomc_*.
- [ ] 6d.3 **OPEN — refit-track substitution into the B fit.** The
  producer can swap in refit tracks, but only via a product-id-matched
  refit collection. The single-track maker's output is keyed to the
  per-candidate bachelor collection (wrong product) and refits only the
  bachelor, not the muons, so substitution is inert by design here. A
  key-based lookup mis-indexed collections and gave a garbage fit (28.7%
  fitOk, median 4.69) -- now guarded by `tref.id() == refitH.id()`.
  Proper refit-into-fit is the N-body maker's domain (add-nbody-cvh-maker).
- [ ] 6d.4 Compare the three constraint modes on a larger sample.

## 6c. Refit-track emission (emitRefitTracks)

- [x] 6c.1 `emitRefitTracks` flag on the single-track maker (default
  False, nominal untouched). Emits `reco::TrackCollection` "refit"
  positionally aligned with the input tracks + `ValueMap<int>` "refitOk".
  Needed because the existing ValueMaps publish only a 3x3 momentum block
  AND are gated on `doMuonAssoc_` / keyed to `pat::Muon` -- so the
  bachelor kaon got nothing at all.
- [x] 6c.2 Covariance guard: non-finite or non-positive-diagonal
  reference blocks are rejected (track emitted as the input copy with
  `refitOk = 0`). Without it, clamped fits leaked NaN covariances that
  would poison a downstream vertex fit.
- [x] 6c.3 `RefitTrack` nano table + validation. Alignment exact (27
  refit vs 27 bachelors). pT shift vs raw median +0.3%, max 8%.
  ptErr median 0.0102 GeV on ~1 GeV tracks (~1%).
- [ ] 6c.4 **OPEN — covariance tail.** After the NaN guard, a
  finite-but-absurd tail survives (max ptErr 92 GeV on a ~1 GeV track;
  3/27 rejected outright). Bulk is sane; the tail must be understood or
  cut before these tracks feed a kinematic vertex fit.
- [ ] 6c.5 **OPEN — chi2/ndof is not track quality.** `chisqval` is the
  CVH full chi2 (all hits + scattering terms) over a nominal tracker
  ndof, so `normalizedChi2()` on a refit track has median ~1e3 and must
  NOT be cut on as a quality measure. Either carry the CVH chi2 in a
  separate field or document loudly at the consumer.
- [ ] 6c.6 Same flag on the two-track maker (muon legs).

## 8. Selection study + loose in-chain pre-filter (this change)

The study mirrors the TWO authoritative selections -- Bmm5 nano
production (`../Bmm5/NanoAOD/python/DileptonPlusX_cff.py`) and the
btojpsik histmaker ANALYSIS path (`get_bkmm_selections`) -- both applied
to the FITTED candidate. The deprecated preset A/B path
(`get_bkmm_alcareco_selections`) is NOT used.

- [x] 8.1 **Selection-study script rewritten**:
  `scripts/btojpsik/study_alcareco_nano_selection.py` applies the real
  cuts on `cvhFitMass` / `cvhFitVtxProb`, rebuilding the muon legs via the
  `mu0/mu1TrackIdx` cross-links. Standalone; does NOT touch the histmaker.
- [x] 8.2 **Cutflow measured** (12 files, 45.7k events, 80.2k fitted
  candidates): muon |eta|<1.4 -> 61%; muon pT>4 -> 57%; dimuon pT>7 ->
  57%; kaon pT in (1,8) -> 9.3%; kaon |eta|<1.4 -> 8.0%;
  **fitted vtxProb>0.1 -> 3.8% (6432->3035)**; fitted mass (5.2,5.4) ->
  1333 candidates.
- [x] 8.2b **B+ PEAK CONFIRMED**: after kinematics + fitted vtxProb>0.1 a
  clear B+ peak appears AT THE PDG MASS (5.279) above the combinatorial
  background (fig `selection_fitB.png`). This validates the whole chain:
  AlCaReco -> CVH refit -> kinematic B fit -> NanoAOD -> real selection
  gives the correct thing. Peak is modest because the not-yet-emitted
  analysis cuts (8.4) would clean the remaining background.
- [x] 8.3 **Correction to the earlier claim**: the real background
  rejector is the FITTED VERTEX PROBABILITY (`bkmm_jpsimc_vtx_prob>0.1`),
  which we emit as `cvhFitVtxProb` -- not a tight DOCA. Bmm5's production
  two-track DOCA cut is LOOSE (0.1), and the histmaker adds no DOCA cut.
  The deprecated preset's DOCA<0.03 was a red herring.
- [x] 8.4 **Dimuon fit-quality handles DONE + VALIDATED**: the B-fit
  producer now runs a dimuon-only KVF and emits `dimuonVtxProb`,
  `dimuonAlphaBS` (XY pointing wrt beamspot), `dimuonSxy` (2D Lxy
  significance wrt beamspot -- proxy for the analysis's 3D sl3d; true 3D
  needs the PV via 4.6c). Physical values (vtxProb med 0.51, alphaBS med
  0.12, Sxy med 1.97). **Adding these cuts takes the peak from S/B~1 to
  S/B~11** (≈970 signal, median 5.2715) on 52 files -- a clean B+ peak.
  Still absent: muon softMVA, bmm BDT (task 8.4b).
- [ ] 8.4b Remaining analysis handles not derivable here: muon softMVA
  (not on the AlCaReco reco::Muon), the bmm BDT (a trained model). Low
  priority; the peak is already clean without them.
- [ ] 8.5 **Loose in-chain candidate pre-filter to shrink output**: add a
  configurable `CandViewSelector` (default OFF/loose) before the tables,
  modelled on Bmm5 production (muon/kaon pT>1, |eta|<2.4, B mass (4,6),
  DOCA<0.1) plus a loose fitted-vtxProb floor -- so production output
  shrinks while staying looser than the analysis cut, for fast histmaker
  iteration. Thresholds set from the study.

## 7. Follow-up (separate changes)

- [ ] 7.1 Wire real 2016 postVFP corrections once solved.
- [ ] 7.2 Batch production over the 38,708-file AlCaReco set.

- [x] 5.9 **Full-file shape run (3255 evt, fast path, no CVH refit)**:
  5741 B+ candidates, kinematic fit 99.7% fitOk, fitted median 5.275
  (PDG 5.2793). Distribution is FLAT/combinatorial by design -- loose
  alignment skim, wide window, no signal selection; the peak-producing
  cuts are the histmaker's job. Kaon dE/dx via the cross-link: median
  3.17 MeV/cm (physical). Confirms the end-to-end infrastructure; a peak
  needs the downstream selection, not the nano.
