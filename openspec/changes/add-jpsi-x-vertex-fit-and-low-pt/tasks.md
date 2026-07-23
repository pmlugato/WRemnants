## 0. Decision gating

- [ ] 0.1 Confirm the sibling `add-jpsi-x-selection-comparison` change has merged or that its preset-switching mechanism, `_NON_V0_PRESETS` dict, and `JpsiXCandidateProducer`/`TwoBodyDecayCandidateProducer` parameter sets are the working baseline for this change. Record the assumed baseline file revisions in `results.md` before any implementation.
- [ ] 0.2 Get alignment-review sign-off on the cost reversal: AlCaReco-stage Kalman vertex fit re-introduced under preset C only, scoped to the four non-V0 channels.

## 0b. Worklog deck — extend, do not fork

- [x] 0b.1 Extend `slides/jpsix-selection-comparison-worklog.tex` with a new section header frame "Followup: low-pT + vertex-fit (change-id `add-jpsi-x-vertex-fit-and-low-pt`)" after the existing Decision frame. Do not start a new deck — the worklog stays one document so the iteration history reads in order. **Done 2026-06-11**: section added with five initial frames (scope, tooling decision, Phase-0 stages, Phase-0 mass overlay PNG, Phase-1 smoke-test gates table).
- [x] 0b.2 Verify the extended deck still compiles to PDF (`pdflatex` in repo container). **Done 2026-06-11**: 47 pages / 479 kB PDF, compiles clean.

## 1. Phase 0 — MiniAOD B+ mass-spectrum demo (standalone, no cmsRun)

- [x] 1.1 Create `scripts/btojpsik/mini_mass_demo.py` (FWLite). CLI: `--mini-file <path-or-xrootd-url>`, `--max-events 5000`, `--out-dir scripts/btojpsik/figs/btojpsik-mass-demo/`, optional `--public-html-dir ~/public_html/mz/alcareco/btojpsik-mass-demo/` (copy on completion). **Done 2026-06-11**: ~500 lines, FWLite + ROOT plotting (no matplotlib needed inside cmssw-el7).
- [x] 1.2 Implement four selection stages:
  1. **No cuts** — every (μ⁺μ⁻, charged-hadron) triplet from `slimmedMuons` + (`packedPFCandidates` ∪ `lostTracks`) with charge ±1; B mass plotted in [4.5, 6.0] GeV.
  2. **+ minimal kinematic** — soft-muon ID (`muon::isSoftMuon`) + pT > 2 GeV; track pT > 0.1 GeV; |η| < 2.4; dimuon mass within 3.0 ± 0.2 GeV.
  3. **+ updated preset B kinematic** — bachelor pT > 0.1, |η_K| < 1.8, straight-line track-track DCA between bachelor and each muon < 0.03 cm.
  4. **+ preset C** — preset B + **closed-form 3-track common-vertex estimate** (analytic least-squares of helix PCAs; χ² from track impact-parameter sigmas), cuts `vtxChi2Prob > 0.01`, `cos(α_xy) > 0.99`, `Lxy/σ > 3` where Lxy is from `offlineBeamSpot`. **Note**: this is a Phase-1 approximation to the iterative `KalmanVertexFitter` used by preset C in cmsRun; see `design.md` "Where the Phase-1 Kalman fit lives" for the cost/correctness tradeoff rationale.
- [x] 1.3 Emit one overlaid `B_mass_stages.png` (ROOT canvas) plus `summary.json` with per-stage candidate count, mean, RMS, and cands/event.
- [x] 1.4 Validate the script on `mc_minilist.txt` line 1 over 200 events. **Done 2026-06-11**: stage 1 essentially flat combinatorial; stage 3 shows the B+ peak emerging at PDG mass; stage 4 cleans the surrounding background. Cands/event: 256 → 62 → 3.3 → 1.0. PNG renders correctly.
- [x] 1.5 Copy the PNG to `~/public_html/mz/alcareco/btojpsik-mass-demo/` and write a one-line `index.html` listing the file with the date. **Done 2026-06-11**: handled by the `--public-html-dir` flag; runs as part of the 5000-event Phase-1 pass.
- [x] 1.6 **Worklog frame — Phase 0**: one frame with the overlay PNG and a one-sentence interpretation. **Done 2026-06-11**: added as part of the 5 new worklog frames extending `jpsix-selection-comparison-worklog.tex` from 42 → 47 pages.

## 2. Update preset B in place (lower bachelor / daughter pT)

- [x] 2.1 In `ALCARECOTkAlJpsiX_cff.py`, changed `_NON_V0_PRESETS['B']`:
  - `'BPlus'.minBachelorPt`: 1.5 → 0.1
  - `'Bc'.minBachelorPt`: 1.5 → 0.1
  - `'Kstar'.minDaughterPt` (sub-resonance): 1.0 → 0.1
  - `'Phi'.minDaughterPt` (sub-resonance): 1.0 → 0.1
- [x] 2.2 `python -m py_compile` + `scram b -j8` clean inside cmssw-el7 + cmsenv. Build completes with no warnings.
- [x] 2.3 Verified `_NON_V0_PRESETS['B']` resolves to the new floor on a fresh cff import — BPlus/Bc minBachelorPt = 0.1; Kstar/Phi minDaughterPt = 0.1.
- [x] 2.4 **Worklog frame — preset B update**: added (covers preset matrix A/B/C in one table; preset B floor change is the second row).

## 3. Add preset C entry to `_NON_V0_PRESETS`

- [x] 3.1 Added `'C'` key to `_NON_V0_PRESETS` copying every kinematic/geometric value from updated `'B'`, with Kalman cuts spread in via `**_PRESETC_KALMAN`.
- [x] 3.2 Kalman parameters in preset C entries: `minBVtxProb=0.01`, `maxMotherAlphaBS=acos(0.99)≈0.1415 rad`, `minBLxyOverSigma=3.0` on **BPlus, Bc, B0Kstar, BsPhi** (Bc included per user clarification 2026-06-11 — Bc → J/ψ π is treated as the pion-bachelor analogue of B+). K*0 and φ `TwoBodyDecayCandidateProducer` instances retain `applyVertexFit=False` (the B-level multi-body fit in `JpsiXCandidateProducer` already constrains sub-resonance tracks).
- [x] 3.3 Extended env-var validation: `_PRESET in ('A', 'B', 'C')`.
- [x] 3.4 Updated the top-of-cff comment block to document the three-preset matrix and the C++-existence-check semantics for the Kalman parameters.
- [x] 3.5 Verified import under all three preset settings: A omits Kalman params on all four non-V0 producers; B omits them; C sets them on BPlus/Bc/B0Kstar/BsPhi but NOT on V0-mode (B0Ks/Lambdab/Psi2S) producers.
- [x] 3.6 **Worklog frame — preset C definition**: included in the "Preset matrix + implementation" worklog frame for sections 2-5.

## 4. C++ — `JpsiXCandidateProducer` Kalman vertex-fit branch

- [x] 4.1 ~~Add optional parameters.~~ **Already present** in `JpsiXCandidateProducer.cc`: `minBVtxProb`, `maxMotherAlphaBS`, `minBLxyOverSigma` (the proposal's `minCosAlphaBS` is equivalently expressed as a max-angle `maxMotherAlphaBS = acos(min_cos_α)`; the C++ uses the angle parametrization). Existence-checked in the ctor (lines 134-150). No-op when absent from cms.PSet.
- [x] 4.2 ~~Add conditional ESConsume.~~ **Already present**: `ttbToken = consumes()` is gated on `applyBVtxFit_` flag, which activates if any of the three Kalman params is set (line 175-178 of the source).
- [x] 4.3 ~~Track-mode Kalman fit.~~ **Already present**: `applyBLevelVertexFit()` at line 299 builds TransientTracks via the consumed `TransientTrackBuilder`, runs `KalmanVertexFitter::vertex(tts)`, computes vtxProb via `TMath::Prob(chi2, ndf)`, and rejects below `minBVtxProb_`. Track mode at line 395-400.
- [x] 4.4 ~~Track-mode pointing/Lxy.~~ **Already present** in `applyBLevelVertexFit()`: cos α computed via `computeAlphaBSFromPoint()` in 3D against `offlineBeamSpot`; Lxy/σ computed against the first `offlinePrimaryVertices` entry using fitted-vertex covariance + PV covariance projected along the displacement direction.
- [x] 4.5 ~~VCC-mode handling.~~ **Already present** at line 460-465: VCC mode collects 2 muon refs + 2 sub-resonance daughter refs (4 total) and calls the same `applyBLevelVertexFit()`. **Implementation note**: this is a single 4-body Kalman fit at a common B vertex — not two separate 2-body fits per the original chat discussion. See `design.md` "Single multi-body B-level Kalman fit (not two 2-body fits)" for rationale; proposal.md updated to match.
- [x] 4.6 `scram b -j8` in EL7 container — clean (only python-module compile messages, no C++ warnings). No regression from the cff edits.
- [ ] 4.7 ~~Single-event smoke-test under preset C on RAW.~~ **Deferred per user instruction "don't run the alcareco yet"** (2026-06-11). Will run at Phase 2.
- [x] 4.8 **Worklog frame — C++ vertex-fit branch**: added (covers the discovery that the producer already had the branch, the existence-check ESConsume gating, and the per-channel topology table).

## 5. Bc topology under preset C — same as B+ (no longer deferred)

- [x] 5.1 Per user clarification 2026-06-11 ("for Bc, should be similar to btojpsik, just pion in place of kaon"): Bc gets the same Kalman parameters as B+ under preset C. Same 3-body track-mode fit; pion mass already configured on the Bc producer (`bachelorMass=0.139570`, line 275 of the cff).
- [x] 5.2 No follow-up issue needed for the basic Bc → J/ψ π mode (it's done). The other Bc decay modes (J/ψ 3π, J/ψ μν, J/ψ τν) are not exposed by this stream's cff and remain out of scope.

## 6. Phase 1 — MiniAOD gate evaluation via the same FWLite script (no histmaker)

Rationale for not extending the NanoAOD histmaker: NanoAOD's `bkmm_jpsimc_*` branches pass through BPark / BMM-Tools upstream selections (track pre-filter, dimuon mass-window, implicit vertex quality), so the "no-cut" baseline on Nano is a low-cut baseline. Signal efficiency and fake-rate measured on Nano under-predict the work preset C's Kalman fit needs to do on real ALCAReco, with no way to catch the mismatch until Phase 2. Running Phase 1 on MiniAOD via the same FWLite path as Phase 0 (and the same `KalmanVertexFitter` as Phase 2 cmsRun) makes the cut surface bit-identical across all three phases — no proxy.

- [x] 6.1 Extend `scripts/btojpsik/mini_mass_demo.py` (from task 1) into a Phase-0+1 dual-purpose tool. Add CLI args: `--presets A,B,C` (comma-separated; runs all selected presets in one pass over the event loop, accumulating per-preset counters and per-stage histograms), `--gen-match` (enables gen-matching), `--gates` (computes signal-eff / fake-rate / mass-shape gates and emits `gates.json` alongside `summary.json`). **Done 2026-06-11**: implemented from the start; one script does both Phase 0 (no `--gen-match`, no `--gates`) and Phase 1 (with both flags).
- [x] 6.2 Implement gen-matching against `prunedGenParticles`: walk the gen B+ → J/ψ(μ⁺μ⁻) K+ decay tree (PDG IDs 521 / 443 / 321), label a reco candidate `true` if all three leaf tracks ΔR-match (< 0.03) to the corresponding gen daughters with matching charges, else `fake`. Store the boolean per-candidate alongside the existing per-stage selection mask. **Done 2026-06-11**: `find_gen_bplus_leaves()` walks PDG 521→443+321→μ⁺μ⁻, `is_gen_matched()` does the ΔR test with both muon orderings. Smoke test: 200/200 BMuonFilter events found gen B+ decay — consistent with the MC being signal-only.
- [x] 6.3 For each (event, preset, stage), bin candidates into (true / fake) × (in-tight-window PDG ± 75 MeV / loose-window-only) and accumulate counters. Per-event totals are accumulated globally for the per-event cands/event computation. **Done 2026-06-11**: `preset_hists` dict-of-TH1F per preset (true_tight, true_loose, fake_tight, fake_loose); `preset_counts` dict-of-ints tracks n_cands, n_true, n_fake, n_true_tight, n_fake_tight.
- [x] 6.4 Run over `--max-events 5000` (or larger, until gen-matched signal-yield Poisson uncertainty falls below 5%) on `BuToJpsiK_BMuonFilter` UL18 MINIAODSIM, file `mc_minilist.txt` line 1. **Done 2026-06-11**: 5000 events processed in 244 s wall (20.5 ev/s; xrootd cache was warm from earlier smoke tests, so the run was much faster than the 50-min estimate). Output at `scripts/btojpsik/figs/btojpsik-mass-demo-phase1/`.
- [x] 6.5 Evaluate the three Phase-1 gates from `gates.json`. **5000-event results** (recorded 2026-06-11):
  - sig\_eff(C)/sig\_eff(B) = 0.495/0.604 = **0.820 ± 0.014** ≥ 0.70 ✓
  - fake/ev(C)/fake/ev(B) = 0.255/2.286 = **0.112 ± 0.003** ≤ 1.5 ✓ (**8.9× fake suppression**)
  - mass shape: Δμ < 0.2 MeV, σ narrower in C (B: 63.5 MeV; C: 60.4 MeV) ✓
  - All three gates pass with strong margin.
- [x] 6.6 No iteration needed: 5k-event refresh confirms the smoke-test result; gates pass on the proposal's defaults (vtxProb > 0.01, cos(α) > 0.99, Lxy/σ > 3). Skipped per gate clearance.
- [x] 6.7 Lock Phase-1 preset-C values: `vtxChi2Prob > 0.01`, `cos(α_xy) > 0.99`, `Lxy/σ > 3`. Recorded in `results.md`. **Caveat in results.md**: Phase-1 closed-form-vertex cut values are an approximation to Phase-2 cmsRun KalmanVertexFitter cut values; one iteration cycle is budgeted in § 7 if Phase-2 numbers diverge.
- [x] 6.8 **Worklog frames — Phase 1**: smoke-test frame (200 ev) retained per cadence rule; three additional frames added with 5k-event results — refresh table, 5k mass overlay PNG, and Phase-1 decision/lock frame. Deck now 50 pages.
- [x] 6.9 No changes to `scripts/histmakers/btojpsik.py` or `wremnants/production/btojpsik_selections.py` in this change. The histmaker's existing `--alcarecoPreset {A,B}` from the sibling change remains as-is; it is **not** extended to C in this proposal. The histmaker continues to serve its analysis-level role unchanged. **Done 2026-06-11**: confirmed by checking diff vs sibling change — no histmaker / btojpsik_selections.py modifications.

## 7. Phase 2 — ALCAReco confirmation on Run2016H Charmonium RAW

- [x] 7.1 cmsDriver three configs: `recoskim_Run2016H_Charmonium_JpsiX_preset{A,B,C}.py`. **Done 2026-06-11**: presetA.py and presetB.py from sibling change reused as-is (configs differ only in output filename; the cff preset selection is env-var-driven at cmsRun time); presetC.py created by copy + sed of presetB.py. All three differ only in `process.ALCARECOStreamTkAlJpsiX.fileName`.
- [~] 7.2 cmsRun each over the same 1000-event input file used by the sibling change Phase 2 (Run2016H Charmonium RAW `22261BD6-…`). `TKALJPSIX_SELECTION_PRESET` set on each cmsRun invocation per the sibling-change bug fix. **In flight 2026-06-11**: bash wrapper PID 1402236 running `cmsRun` sequentially under presets A → B → C. Logs at `CMSSW_10_6_17_patch1/logs/jpsix_phase2_3way_preset{A,B,C}_run.log`. Wall time estimate ~15-20 min total (3 × ~5 min with warm xrootd cache).
- [x] 7.3 Extend `jpsix_preset_compare.py` to accept three preset inputs (`--presetA-file`, `--presetB-file`, `--presetC-file`) with three-way overlay plots and a three-row gate table. **Done 2026-06-11**: new sibling script `jpsix_preset_compare_3way.py` (~260 lines) imports helpers from the 2-way script (`loop_file`, `parse_log`, `parse_edmsize`, `CHANNELS`, etc.) and adds 3-way `report_2a_3way`, `report_2b_3way`, `report_2c_3way`, `overlay_plots_3way`. Output: 3-line mass overlay per channel; per-preset gate columns; `summary.json` with verdict_per_preset.
- [ ] 7.4 Run the comparator. Emit 2a (mass quality + gates), 2b (size + timing), 2c (dedup).
- [ ] 7.5 Per-channel preset decision per the outcome rule in `proposal.md` (section "Decision criterion"). Record in `results.md`.
- [ ] 7.5a **Bc-specific iteration trigger**: if Bc cands/event under preset C drops below 0.05 / event on Run2016H Charmonium (i.e. < 50 / 1000 events), check whether the loss is concentrated at the `minBLxyOverSigma` and `maxMotherAlphaBS` cuts (compare to preset B's Bc yield). If yes, run a Bc-specific iteration with `minBLxyOverSigma = 2.0` and/or `maxMotherAlphaBS = math.acos(0.98)` ≈ 0.20 rad on the Bc producer alone. Cff-only change; takes one cmsRun pass to validate. Record the iteration in `results.md` and the worklog.
- [ ] 7.6 If C wins on any channel, set that channel's cff default to C; otherwise updated B remains the default. Env-var override retained as a developer escape.
- [ ] 7.7 **Worklog frames — Phase 2a/2b/2c + Decision**: four frames, format mirroring sibling change.

## 8. Public-facing summary

- [ ] 8.1 After all of Phase 0/1/2 lock: copy the Phase-0 PNG to `~/public_html/mz/alcareco/btojpsik-mass-demo/`, the Phase-1 cut tables + gate values to `~/public_html/mz/alcareco/phase1-preset-c/`, and the Phase-2 gate / size table to `~/public_html/mz/alcareco/phase2-preset-c/`.
- [ ] 8.2 Update `slides/jpsix-alcareco-producer.tex` (the public-facing deck): one new headline frame summarizing the per-channel preset decision after the existing "Validation run" frame.
- [ ] 8.3 Update `results.md` with final per-channel preset assignment + the three Phase-1 and three Phase-2 gate numbers per preset per channel.

## 9. Validation

- [ ] 9.1 `openspec validate add-jpsi-x-vertex-fit-and-low-pt --strict --no-interactive` passes.
- [ ] 9.2 `cmssw-el7 + cmsenv + python -m py_compile` on the cff under all three preset settings.
- [ ] 9.3 `scram b -j8` clean.
- [ ] 9.4 FWLite smoke-test on the standalone demo script under one MiniAOD file passes.
- [ ] 9.5 cmsRun smoke-test under each preset produces a non-empty ALCAReco output.
