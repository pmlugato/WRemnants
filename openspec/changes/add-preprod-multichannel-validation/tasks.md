## 0. Pre-validation refactors (before any smoke)

- [x] 0.1 Preset machinery deleted from `ALCARECOTkAlJpsiX_cff.py`;
      preset-B literals inlined.
- [x] 0.2 Env-var references left in the obsolete
      `condor/jpsix_alcareco/` scripts (superseded by new
      multichannel dir); cff ignores env var regardless.
- [x] 0.3 `ALCARECOTkAlJpsiXGoodMuons.filter = cms.bool(False)`;
      `ALCARECOTkAlJpsiXLooseMuons.filter = cms.bool(False)`.
- [x] 0.4 `VertexCompositeCandidateMerger.cc` plugin written
      (~30 lines) and added to `BuildFile.xml`. Cff wires
      `ALCARECOTkAlJpsiXAnyCandidate` (merger) and
      `ALCARECOTkAlJpsiXCandFilter` (`CandViewCountFilter`,
      `minNumber = 1`) between the 8 candidate producers and
      `CompositeDaughterTrackProducer`.
- [x] 0.5 `scram b -j 16` — clean build; both new plugins registered.
- [x] 0.6 Rate-estimate smoke: 1000 events, 65.4% accepted, < 2x prior smoke -> ADOPTED.
- [x] 0.7 Bit-invariance side-check: per-channel candidate rates match earlier smokes on shared events (Bc 2.60/evt, J/psi-only 1.00/evt, etc).
- [x] 0.8 David's `444b20f` invariants verified in the code:
      `signedPdgId` (both producers), `jpsiDaughterKeys` VCC-mode
      veto, VCC remapper vertex-cov + chi2/ndof passthrough,
      KalmanVertex `VertexException` guard, D* selector `LogInfo`
      summary. ψ(2S) V0-mode stale comments removed with the
      preset-machinery cleanup.

## 1. Environment and inventory

- [x] 1.1 `ALCARECOTkAlKsToPiPi_cff.py` and
      `ALCARECOTkAlLambdaToProtonPi_cff.py` already present in the
      sparse tree AND already aligned to the persist policy
      (cloned tracks, 3 dE/dx maps, candidate VCC `Resonances`).
      No further work needed.
- [x] 1.2 `condor/multichannel_alcareco/lfns/singlemuon_2016h.txt`
      cached — 33165 SingleMuon 2016H RAW files.
- [x] 1.3 `condor/multichannel_alcareco/lfns/charmonium_2016h.txt`
      cached — 11833 Charmonium 2016H RAW files.
- [x] 1.4 xrootd open confirmed (SingleMuon populated file 281707).
- [x] 1.5 Populated-file filter: DAS + nevents >= 1000 -> 10565
      Charmonium files, 31916 SingleMuon files at
      `condor/multichannel_alcareco/lfns/*_populated.txt`. Some 2016H
      RAW files are zero-event placeholders and get skipped.

## 2. Stream persist-alignment (production-permanent)

- [x] 2.1 `ALCARECOTkAlJpsiMuMu_cff.py`: three
      `DeDxValueMapProjector` clones + `alcaDedxJointEstimator`
      added to the sequence; three `keep` lines added to
      `_Output_cff.py`.
- [x] 2.2 `ALCARECOTkAlJpsiMuMu_cff.py`:
      `ALCARECOTkAlJpsiMuMuTrackToMuon` producer added;
      `keep *_ALCARECOTkAlJpsiMuMuGoodMuons_*_*` and
      `keep *_ALCARECOTkAlJpsiMuMuTrackToMuon_*_*` added to
      `_Output_cff.py`.
- [x] 2.3 `ThreeBodyDecayCandidateProducer.cc` plugin written
      (~150 lines) and added to `BuildFile.xml`. D* cff wires
      `ALCARECOTkAlDstToD0PiResonances` emitter with the same
      mass windows, Q-value, per-track pt cuts, charge convention
      as the internal selector. `keep` line added to
      `_Output_cff.py`.
- [x] 2.4 KS + Λ streams already aligned to the policy at file
      creation time; no changes needed.
- [x] 2.5 `scram b -j 16` — clean build.
- [x] 2.6 D* candidate emitter: 6806 D* candidates on 895 accepted events (7.6/evt).

## 3. Recoskim configs

- [x] 3.1 Charmonium recoskim generated
      (`recoskim_Run2016H_Charmonium_JpsiX_JpsiMuMu.py`).
- [x] 3.2 SingleMuon recoskim generated
      (`recoskim_Run2016H_SingleMuon_DstV0.py`).
- [x] 3.3 `edmDumpEventContent` audit: J/psi+X output has all 12
      persisted products (8 VCC channels + tracks + 3 dE/dx +
      loose muons + track-muon assoc); TkAlJpsiMuMu adds 3 dE/dx
      maps + tight muons + track-muon assoc; D* adds candidate VCC.

## 4. Smoke tests

- [x] 4.1 Charmonium smoke: 1000 events, 654 accepted (J/psi+X),
      678 accepted (JpsiMuMu), 331s wall, peak RSS 4.2 GB, exit 0.
- [x] 4.2 SingleMuon smoke: 1000 events on Run 281707, D* 895/1000,
      KS 986/1000, Lambda 936/1000, 73s wall, peak RSS 4.3 GB, exit 0.
- [x] 4.3 Persist-alignment audited; all target categories present.

## 5. Plotting

- [x] 5.1 `muon-compare` executed on the Charmonium smoke:
      637 shared events. Result: 1412 tight-and-loose, 104 loose-
      only, 0 tight-only. 7 per-observable PNG+PDF overlays.
- [x] 5.2 `kinematics` executed on both smokes: 12 per-channel
      plots (8 J/psi+X + JpsiMuMu + D* + KS + Lambda). Vertex
      chi2/ndof shown where the producer emits a fit (V0-mode +
      remapper passthrough).
- [x] 5.3-5.6 Sub-commands available (`dedx`, `cut-efficiency`,
      `v0-baseline`, `resources`) — run after 100-file batch
      outputs land.

## 6. Condor infrastructure

- [x] 6.1 `condor/multichannel_alcareco/build_tarball.sh` written
      (mirrors `condor/jpsix_alcareco/build_tarball.sh` shape;
      ships both recoskims).
- [x] 6.2 `condor/multichannel_alcareco/run.sh` written —
      dispatches on `charmonium` | `singlemuon`, monitors peak
      RSS, emits per-file JSON, xrdcp to
      `/ceph/submit/data/user/p/pmlugato/mz/alcareco/multichannel/`.
- [x] 6.3 `condor_multichannel_charmonium.sub`,
      `condor_multichannel_singlemuon.sub` written with
      `RequestMemory=8000`. Both point at
      `charmonium_100files.txt` / `singlemuon_100files.txt`
      (100-file lists prepared from cached DAS output).
- [x] 6.4 Charmonium 100-file batch submitted (cluster 3227850).
- [x] 6.5 SingleMuon 100-file batch (cluster 3227853).
- [x] 6.6 Both batches drained. 120/200 successful (72 Charmonium
      exit=0, 48 SingleMuon exit=0). 80 failed on xrootd
      `EX_NOINPUT` — recommend narrower site list before full
      launch.

## 7. Projections and report

- [x] 7.1 `projections.py` on 120 successful JSONs -->
      `projections.json`:
      - Charmonium: 196 GB, 1614 CPU-hrs, 0.34 wall days @ 200 slots.
      - SingleMuon: 162 GB, 1145 CPU-hrs, 0.24 wall days @ 200 slots.
      - Combined: 358 GB, 2758 CPU-hrs, $\lesssim$1 wall day @ 200 slots.
- [x] 7.2 Project memory updated.
- [x] 7.3 Slides updated (`slides/preprod-multichannel-validation.pdf`,
      17 pages, no overfull boxes) with muon-compare + kinematics
      + resources + projections.
- [ ] 7.4 User review pending.

## 8. Approval gate

- [ ] 8.1-8.3 User approval pending for muon plot / projections /
      full-2016 launch. Recommend: narrow site list further
      (drop sites that emit `EX_NOINPUT` at input-open time), then
      submit full-2016 Charmonium + SingleMuon batches.
