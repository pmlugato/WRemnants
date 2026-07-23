## 1. Phase 1 — Python configuration (no recompile)

- [x] 1.1 Tighten J/ψ mass window in `ALCARECOTkAlJpsiXJpsiCandidates`: `minMass = 2.95`, `maxMass = 3.25` GeV
- [x] 1.2 Tighten K*0 mass window in `ALCARECOTkAlJpsiXKstarCandidates`: `minMass = 0.80`, `maxMass = 0.99` GeV
- [x] 1.3 Tighten φ mass window upper edge in `ALCARECOTkAlJpsiXPhiCandidates`: `minMass = 0.990`, `maxMass = 1.040` GeV
- [x] 1.4 B+ bachelor pT: no change (keep `minBachelorPt = 0.5` GeV, matches HLT floor)
- [x] 1.5 Lower Bc bachelor pT threshold: `minBachelorPt = 0.3` GeV in `ALCARECOTkAlJpsiXBcCandidates`
- [ ] 1.6 Re-run test job and confirm visible B+ and Bs→φ peaks; document yields

## 2. Phase 2 — C++ extension to JpsiXCandidateProducer and TwoBodyDecayCandidateProducer

### JpsiXCandidateProducer.cc

- [ ] 2.1 Add `minJpsiPt_` member and optional `minJpsiPt` cfg parameter (default 0) to `JpsiXCandidateProducer.cc`
- [ ] 2.2 Apply `minJpsiPt` cut on J/ψ candidate before the bachelor/VCC loop; configure **3.0 GeV** for all channels
- [ ] 2.3 Add `minMotherPt_` member and optional `minMotherPt` cfg parameter (default 0)
- [ ] 2.4 Apply `minMotherPt` cut after mother 4-momentum is computed (both modes)
- [ ] 2.5 Add `maxBachelorEta_` member and optional `maxBachelorEta` cfg parameter (default +∞, track mode only)
- [ ] 2.6 Apply `maxBachelorEta` cut on bachelor track |η| before combining (track mode only)
- [ ] 2.7 Add `maxJpsiAlphaBS_` member and optional `maxJpsiAlphaBS` cfg parameter (default +∞); consume `offlineBeamSpot` InputTag when parameter is present
- [ ] 2.8 Compute J/ψ pointing angle to beamspot using J/ψ VCC vertex midpoint; apply `maxJpsiAlphaBS` cut
- [ ] 2.9 Add `maxBachelorIPToJpsiVertex_` member and optional `maxBachelorIPToJpsiVertex` cfg parameter (default +∞, track mode only)
- [ ] 2.10 Compute 3D DCA of bachelor helix to J/ψ vertex midpoint; apply cut (track mode only)
- [ ] 2.11 `scram b` and confirm no regressions in V0-mode channels

### TwoBodyDecayCandidateProducer.cc

- [ ] 2.12 Add `maxDaughterEta_` member and optional `maxDaughterEta` cfg parameter (default +∞)
- [ ] 2.13 Apply `maxDaughterEta` cut on both daughter tracks before mass-window combinatorics
- [ ] 2.14 `scram b` and confirm K*0 and φ channel yields are reduced as expected

### Python configuration

- [ ] 2.15 Update `ALCARECOTkAlJpsiX_cff.py`: set `minJpsiPt = 5.0` GeV for all 7 channels
- [ ] 2.16 Update `ALCARECOTkAlJpsiX_cff.py`: set `minMotherPt = 5.0` GeV for B+, B0→K*0, Bs, Bc channels
- [ ] 2.17 Update `ALCARECOTkAlJpsiX_cff.py`: set `maxBachelorEta = 2.4` for B+ and Bc channels
- [ ] 2.18 Update `ALCARECOTkAlJpsiX_cff.py`: set `maxJpsiAlphaBS = 1.0` for all 7 channels (requires `offlineBeamSpot` token)
- [ ] 2.19 Update `ALCARECOTkAlJpsiX_cff.py`: set `maxBachelorIPToJpsiVertex = 1.0` cm (= 10 mm) for B+ and Bc channels
- [ ] 2.20 Update `ALCARECOTkAlJpsiXKstarCandidates` and `ALCARECOTkAlJpsiXPhiCandidates`: set `maxDaughterEta = 2.4`
- [ ] 2.21 Re-run test job; document per-channel yields and compare to Phase 1 baseline

## 3. Spec and documentation updates

- [ ] 3.1 Update spec delta in this change (Stage-1 and Stage-2 requirements reflect new cut values) — DONE
- [ ] 3.2 Run `openspec validate add-jpsi-x-candidate-quality-cuts --strict --no-interactive`
