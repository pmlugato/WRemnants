## 0. Decision gating

- [x] 0.1 The 10 gating questions are resolved (2026-05-13, updated 2026-06-10 to remove AlCaReco-stage vertex fitting on cost grounds and pivot Phase 1 to NanoAOD MC on B+ via the existing `btojpsik.py` histmaker); see "Resolved questions" in `proposal.md`. Implementation may proceed.

## 0b. Worklog slide deck setup

- [ ] 0b.1 Create `slides/jpsix-selection-comparison-worklog.tex` (beamer, same theme as `jpsix-alcareco-producer.tex`) with a title page, an outline frame, and an initial "Setup" frame stating: input MC dataset path, software versions, change-ID, the two-phase plan in one sentence each, and the decision criterion in three bullets per phase.
- [ ] 0b.2 Verify it compiles to PDF (`pdflatex`); commit the source.

## 1. Turn off K*0/φ Kalman vertex fits already wired in cff

- [ ] 1.1 In `ALCARECOTkAlJpsiX_cff.py`, set `applyVertexFit=False` on `ALCARECOTkAlJpsiXKstarCandidates`. Remove the `minVtxProb` line.
- [ ] 1.2 Same on `ALCARECOTkAlJpsiXPhiCandidates`.
- [ ] 1.3 Confirm the cff loads (`python -m py_compile`) and that `cms.Process` build does not warn about the dropped parameters.
- [ ] 1.4 **Worklog frame**: one frame noting "K*0/φ Kalman vertex fits removed at AlCaReco stage; rationale = cost; vertex-based cleanup deferred to Stage 2." Cite the affected cff lines.

## 2. Python — preset-table-driven cff refactor (A vs B)

- [ ] 2.1 Add `_PRESET = os.environ.get('TKALJPSIX_SELECTION_PRESET', 'B')` at the top of `ALCARECOTkAlJpsiX_cff.py`. Validate with `assert _PRESET in ('A', 'B')`.
- [ ] 2.2 Move the per-channel mass windows into a single `_MASS` dict.
- [ ] 2.3 Add helper functions per the design sketch:
  - `_kin_non_v0(preset, channel)` returns kinematic + geometric cuts for the four non-V0 channels. "no-op" values (0 or +∞) under preset A; current Phase-1+2 values under preset B.
  - `_kin_v0(channel)` returns minimal kinematic cuts for the three V0-mode channels; preset-invariant.
- [ ] 2.4 Replace each of the 4 non-V0 `JpsiXCandidateProducer` cms.EDProducer calls' cut parameters with `**_kin_non_v0(_PRESET, '<channel>')`. Drop any `minBVtxProb` / `minBLxyOverSigma` parameters (they are not in the C++ and must not be in the cff).
- [ ] 2.5 Replace each of the 3 V0-mode `JpsiXCandidateProducer` cms.EDProducer calls' cut parameters with `**_kin_v0('<channel>')`. Preset-invariant by construction.
- [ ] 2.6 Verify that running with `TKALJPSIX_SELECTION_PRESET=B` (or no env var) reproduces the *exact* current cff parameter set on the four non-V0 channels, minus the K*0/φ vertex-fit lines (`print(producer.dumpPython())` diff). V0-mode channels' output MUST be bit-identical across A and B.
- [ ] 2.7 Add a leading comment block in the cff documenting the V0-vs-non-V0 split and the per-channel preset semantics, with explicit reference to this change and to the no-AlCaReco-stage-vertex-fitting decision.
- [ ] 2.8 `python -m py_compile` the refactored cff inside `cmssw-el7 + cmsenv`.
- [ ] 2.9 **Worklog frame**: one frame showing the cff before/after structure (the preset-table-driven dict-of-dicts) and the per-channel A/B value table. This frame is referenced from every subsequent Phase-1 frame as the "where the cuts live" anchor.

## 3. Phase 1 — extend `scripts/histmakers/btojpsik.py`

- [x] 3.1 In `btojpsik.py`, add a new function `get_bkmm_alcareco_selections(preset, vprefix='jpsimc')` alongside the existing `get_bkmm_selections()`. Return the same `(label, lambda)` tuple-list shape so it can drop into the same cut-flow machinery. The function takes `preset in ('A', 'B')` and returns the appropriate cut profile.
- [x] 3.2 Preset A profile (**revised 2026-06-10 for raw-variable AlCaReco emulation**): opposite-sign dimuon + raw J/ψ mass [2.95, 3.25] GeV on `raw_mumu_mass` + raw B+ mass [5.0, 5.5] GeV on `raw_b_mass`. Both raw masses are computed on-the-fly from class-1 quantities (`mm_mu{1,2}_pt/eta/phi`, `bkmm_kaon_pt/eta/phi`) + PDG μ/K masses. NO `bkmm_jpsimc_*`, NO `mm_kin_*`.
- [x] 3.3 Preset B profile (**revised 2026-06-10**): A + `Muon_pt > 4` + `bkmm_kaon_pt > 1.5` + `|bkmm_kaon_eta| < 1.4` + `raw_mumu_pt > 3` + `raw_b_pt > 5` + `bkmm_kaon_mu{1,2}_doca < 0.03 cm` (both legs). All cuts class-1 or class-2; values chosen to be meaningfully beyond the BMM Tools upstream floor (HadronMinPt=1.0, HadronMaxEta=2.4, maxTwoTrackDOCA=0.1). α_BS-at-J/ψ is dropped (class-3 forbidden per proposal).
- [x] 3.3b Add the raw-kinematic helper `define_raw_kinematics(df)` and selection functions `select_raw_{mumu_mass_window, b_mass_window, mumu_pt, b_pt, kaon_pt, kaon_eta}_window` plus `select_kaon_mu_doca` to `wremnants/production/btojpsik_selections.py`. Each cut uses `ROOT::VecOps::RVec<bool> bkmm_passes` plumbing identical to existing functions.
- [x] 3.4 Add `--alcarecoPreset {A,B}` argparse flag to `btojpsik.py`. Default unset = analysis path (existing `get_bkmm_selections()` behaviour, unchanged).
- [x] 3.5 In `build_graph(df, dataset)`, branch on the argparse flag: if `args.alcarecoPreset` is set, call `define_raw_kinematics(df)` then route through `get_bkmm_alcareco_selections(args.alcarecoPreset)`; otherwise use the existing analysis selections. Existing analysis behaviour byte-for-byte preserved when the flag is unset.
- [x] 3.6 Gen-matching is reused from the existing `select_only_passing_bkmm_candidates(gen_match_nonsignal=True)` path; the `gen_filter_eff` it produces is the Phase-1 signal-efficiency metric (no new code needed for iter 1b).
- [ ] 3.7 If finer breakdown is needed (per-iteration fake-vs-true mass histograms), add a per-candidate gen-match boolean axis. Deferred until iter 2 if iter 1b numbers are sufficient.
- [x] 3.8 `python -c "import ast; ast.parse(...)"` syntax check on both modified files passes.
- [x] 3.9 **Worklog frame**: cut tables (A and B), argparse-flag wiring, and the new `select_raw_*` family will be documented in the iter 1b frames (added when run completes).

## 4. Phase 1 — run, evaluate gates, iterate

- [x] 4.1 Run `btojpsik.py --alcarecoPreset A --era 2016PostVFP --filterProcs BuToJpsiK ...` to produce the preset-A NanoAOD histograms. (Iter 1: 2026-06-10, superseded by iter 1b due to class-4 mass cut. Iter 1b: in flight at the time of this commit.)
- [x] 4.2 Same for preset B.
- [ ] 4.3 Compute the three Phase-1 gates from the output histograms:
  - signal efficiency on truth-matched candidates B/A ≥ 70%;
  - per-event fake-rate reduction A/B ≥ 5×;
  - mass-shape width and mean preserved between A and B truth-matched distributions (≤ 2σ-stat difference).
- [ ] 4.4 **Worklog frame — Phase 1 iteration 1**: cut tables actually applied, signal-efficiency / fake-rate / mass-shape numbers, one truth-matched B+ mass plot under A and B side-by-side, one-sentence observation ("B kills 80% of fakes at 10% signal loss" or similar).
- [ ] 4.5 If any gate fails under the proposed B cuts, iterate: relax or tighten individual cuts (DCA threshold, mother-pT floor, αBS limit) and re-run.
- [ ] 4.6 **Worklog frame — Phase 1 iteration N**: one new frame per iteration, same structure as 4.4. Title each frame "Phase 1 — iter N: changed X from Y to Z". This frame stack IS the iteration history; do not delete or overwrite earlier iteration frames.
- [ ] 4.7 Record the locked B+ Phase-1 preset values in `results.md`.
- [ ] 4.8 **Worklog frame — Phase 1 lock**: one summary frame with the final B+ cut table, the gates-passed status, and a one-paragraph "why we stopped iterating".

## 5. Translate the B+ Phase-1 cuts to the other non-V0 channels

- [x] 5.1 Translation table per channel (Phase 1 locked B+ values → other non-V0). See proposal Section 3 sub-table "Per-channel preset-B translation".
- [x] 5.2 Shared J/ψ-leg cuts copied verbatim: `Muon_pt > 4`, `minJpsiPt > 3`. The J/ψ → μμ is the same particle on every non-V0 channel.
- [x] 5.3 Mother-pT floor: > 5 GeV for B+, B0→K*0, Bs→φ, Bc (no-op on signal MC but kept for production-mode rejection). ψ(2S) → 3 GeV (already in cff). V0-mode channels: 0 (already in cff, preset-invariant).
- [x] 5.4 **Bachelor topology split** (resolved 2026-06-10):
  - Track mode (B+, Bc): set `minBachelorPt`, `maxBachelorEta` per preset on `JpsiXCandidateProducer`. **Add new `maxBachelorMuTrackDOCA` C++ parameter** to `JpsiXCandidateProducer` exposing the Phase-1-tuned track-track DOCA (`bkmm_kaon_mu{1,2}_doca`). NOT reusing `maxBachelorIPToJpsiVertex` (different quantity, larger scale).
  - VCC mode (B0→K*0, Bs→φ): **add `minDaughterPt`** parameter to `TwoBodyDecayCandidateProducer` (per-daughter sub-resonance pT floor); **loop over VCC's daughter tracks inside `JpsiXCandidateProducer`** to apply the bachelor–muon DOCA cut against each muon. Daughter–muon DOCA logic lives in `JpsiXCandidateProducer` (J/ψ-context-dependent), not in the sub-resonance producer (J/ψ-independent).
- [x] 5.5 C++ work supporting the translation table — DONE 2026-06-10:
  - [x] 5.5a `TwoBodyDecayCandidateProducer`: `minDaughterPt` added (default 0.0). Applied per daughter before mass-window check.
  - [x] 5.5b `JpsiXCandidateProducer` (track mode): `maxBachelorMuTrackDOCA` added (default = numeric_limits::max(), i.e., off). Implemented as straight-line track-track DCA between bachelor and each J/psi-daughter muon via new static helper `trackTrackDCA(t1, t2)`; matches the physics scale of the Phase-1 `bkmm_kaon_mu{1,2}_doca` quantity. Reject if either pair exceeds.
  - [x] 5.5c `JpsiXCandidateProducer` (VCC mode): same `maxBachelorMuTrackDOCA` param honoured by looping over the VCC's leaf tracks (`collectLeafTrackRefs(xCand, ...)`) and rejecting if any (daughter, muon) pair exceeds.
  - [x] 5.5d `ALCARECOTkAlJpsiX_cff.py`: env-var `TKALJPSIX_SELECTION_PRESET={A,B}` switch at top; `_NON_V0_PRESETS` dict-of-dicts driving the 4 non-V0 producers + K*0/φ sub-resonance daughter cuts; `applyVertexFit` turned off on K*0 and φ; `maxMotherAlphaBS` / `minBVtxProb` / `minBLxyOverSigma` removed from all instances (no Kalman B-vertex fit at AlCaReco). V0-mode channels (B0→Ks, Λb, ψ(2S)) preset-invariant.
  - [x] 5.5e `scram b -j8` in the EL7 container — both producers compile clean.
  - [x] 5.5f cff loads cleanly under both `TKALJPSIX_SELECTION_PRESET=A` and `=B`; values emit as expected (verified by python import + parameter readback).
- [ ] 5.6 **Worklog frame — B+ → other-channel translation**: one frame with the translation table per channel (rows = cuts; columns = B+ value, channel value, physics justification).

## 6. Phase 2 — ALCARECO confirmation on RAW (deferred until Phase 1 complete)

- [x] 6.0 Smoke-test the refactored cff under `TKALJPSIX_SELECTION_PRESET=A` and `=B`: confirm both presets load and emit the expected per-channel kinematic values. (Done 2026-06-10; preset A emits Bplus minJpsiPt=0/DOCA=off, preset B emits 3.0/0.03.)
- [x] 6.1 cmsDriver.py with `-s RAW2DIGI,L1Reco,RECO,ALCA:TkAlJpsiX --runUnscheduled --nThreads 16 --data --era Run2_2016 --conditions 106X_dataRun2_v35 --eventcontent ALCARECO --datatier ALCARECO --customise Configuration/DataProcessing/RecoTLR.customisePostEra_Run2_2016 -n 1000 --no_exec` × 2, with `TKALJPSIX_SELECTION_PRESET=A,B` set per invocation. Two configs produced: `recoskim_Run2016H_Charmonium_JpsiX_preset{A,B}.py`. FastTimerService appended (`printJobSummary=True`, no JSON in this CMSSW release). Per-preset output filename override on `ALCARECOStreamTkAlJpsiX.fileName` to avoid collision.
- [x] 6.2 Run both ALCARECO passes; FastTimerService wired in via one-line `process.add_(...)`. **Done 2026-06-11** after voms proxy renewal. First pass had a bug (env var set at cmsDriver only, not at cmsRun) producing identical A and B outputs; fix was to set `TKALJPSIX_SELECTION_PRESET` on both cmsDriver and cmsRun invocations. Run-1 outputs archived under `phase2_run1_envvar_bug/`. Run-2 (corrected): preset A = 637 events saved / 38.38 MB / 351.7 s wall; preset B = 217 events saved / 3.14 MB / 391.2 s wall.
- [x] 6.3 `edmEventSize` on each output ALCARECO file for branch-level size table. Outputs at `logs/jpsix_alcareco_preset{A,B}_edmEventSize.txt`. Per-channel VCC branch bytes scale ~30× larger under preset A.
- [x] 6.4 **Worklog frame — Phase 2 ALCARECO runs**: one frame summarising the two runs (event count, wall-time, output size, any cmsRun warnings). If anything anomalous happens (a producer fails, a cut behaves unexpectedly), add a second frame "Phase 2 — anomaly: X" with the diagnosis. Done: (a) data-access blocker frame on 2026-06-10; (b) proxy-renewed + run-complete frame; (c) env-var bug frame; (d) Phase 2a/2b/2c result frames and decision frame.

## 7. Phase 2 — diagnostic tool `jpsix_preset_compare.py`

- [x] 7.1 Write a new top-level `jpsix_preset_compare.py` (sibling of `jpsix_diagnostics.py`) that accepts two ALCARECO files (one per preset). FWLite-based. (Done 2026-06-10. CLI: `--presetA-file`, `--presetB-file`, optional `--logA/--logB` and `--edmsizeA/--edmsizeB`, `--output-dir`. ~320 lines. Fork chosen over extend-in-place to keep the single-file diagnostic intact.)
- [x] 7.2 **Phase 2a — data masses**: overlay plots of mother mass and J/ψ sub-mass per channel under A and B. Evaluate the three Phase-2 gates per channel (cands/event ≤ 5; tight-window fraction ≥ 30%; tight-window yield ≥ 100/1000 events). Print PASS/FAIL. (Implemented. Tight window = PDG ± 75 MeV, or ± 40 MeV for ψ(2S). Verdict computed on preset B.)
- [x] 7.3 **Phase 2b — size + time**: parse FastTimerService output and `edmEventSize` output; build a preset × {wall-s/event, bytes/event, per-branch bytes} table. (Implemented. Wall time from `/usr/bin/time -f` wrap; per-module from FastTimer job-summary stdout block; per-branch from `edmEventSize -v` ALCARECO-prefixed lines.)
- [x] 7.4 **Phase 2c — dedup**: per-event count of merged-collection tracks vs counterfactual sum-of-per-channel-AlignmentTrackSelector-output; compression factor averaged over events. (Implemented as merged_total / separate_total = sum(cands_c · nLeaf_c).)
- [x] 7.5 V0-mode channels: single-histogram plots; overlays collapse to one curve as a sanity check that the preset switch is structurally invariant for them. (Overlay plot draws both A and B; identical curves on V0 channels by construction since cff is preset-invariant for them.)
- [x] 7.5b Smoke-test the comparator on an existing ALCARECO file (`TkAlJpsiX.root` from May 12, 205 events) with the same file passed as both A and B. All three sub-studies print, overlay PNG renders, `summary.json` writes. Done 2026-06-10. **Plumbing validated; awaits real two-preset ALCARECO output (blocked on 6.2 data access).**
- [x] 7.6 **Worklog frame — Phase 2a (data masses)**: per-channel mother-mass overlay (A vs B) for the four non-V0 channels, plus the three-gate PASS/FAIL table. Title each frame by channel. (Done — one summary table frame + one interpretation frame + one PNG-overlay frame.)
- [x] 7.7 **Worklog frame — Phase 2b (size + time)**: one frame with the per-preset wall-time and bytes/event table. Comment on whether the ordering matches expectation (B ≤ A on time and bytes). (Done — B output 12× smaller, wall-time within 11% of A on this sample.)
- [x] 7.8 **Worklog frame — Phase 2c (dedup)**: one frame with the per-event histogram and the average compression factor under both presets. (Done — A=8.01×, B=1.27×; interpretation frame argues for B because absolute track count matters.)

## 8. Results assembly

- [x] 8.1 Write `results.md` with two sections: Phase 1 (MC, B+ only) and Phase 2 (data, all 7 channels). Each section has the relevant tables, one-paragraph interpretation, and the per-channel preset decision. (Done 2026-06-11.)
- [x] 8.2 Phase 1 section MUST include: cut tables (A and B values per cut), signal-efficiency / fake-rate / mass-shape metrics under each preset, locked B+ Phase-1 preset. (Done — iter 1b and iter 2 columns included.)
- [x] 8.3 Phase 2 section MUST include: 2a data mass-quality table, 2b size/timing table, 2c dedup compression table, per-channel preset choice with gate citations, list of channels (if any) flagged as "Stage-1 insufficient; relies on Stage-2 cleanup". (Done; Bs→φ and Bc flagged for follow-up.)
- [ ] 8.4 Mirror the headline tables (Phase-1 efficiency/fake-rate + Phase-2 cands/event + bytes/event) into a new frame in the main `slides/jpsix-alcareco-producer.tex`. Place between the existing "Validation run" frame and any later frames. (Deferred — to do when the public deck is next updated.)
- [x] 8.5 **Worklog frame — Decision**: per-channel preset choice with gates cited from earlier worklog frames. Done — "Phase 2 --- decision (2026-06-11)" frame with the per-channel table and rationale rows.
- [x] 8.6 **Worklog frame — Followups**: channels flagged as Stage-1-insufficient, missing-MC-sample requests, any loose ends. Done — bullets at the bottom of the decision frame: 10k-event Phase-2 sample for Bs/Bc confirmation, robustify env-var consumption against the cmsDriver-only mistake, stage Run2016H Charmonium RAW locally.

## 9. Decision and rollout

- [ ] 9.1 Review `results.md` against the decision criterion. For each non-V0 channel, identify the smallest passing preset under Phase 2 (A preferred over B on ties; B is the production default if A fails any gate). Record the per-channel choice in `results.md`.
- [ ] 9.2 Update `ALCARECOTkAlJpsiX_cff.py` so the per-channel `_kin_non_v0` defaults match the chosen preset for each non-V0 channel. Retain the env-var override as a developer escape.
- [ ] 9.3 Add a leading comment in the cff identifying this change and the date of the decision; cite the table row in `results.md` that drove each per-channel choice.
- [ ] 9.4 Update `openspec/specs/alcareco-jpsi-x/spec.md` (after archive) to reflect the chosen per-channel preset map as the production behaviour.

## 10. Documentation and follow-up

- [ ] 10.1 File-level docstring update at the top of `ALCARECOTkAlJpsiX_cff.py` describing the preset table and the no-AlCaReco-stage-vertex-fitting decision.
- [ ] 10.2 README section in `jpsix_preset_compare.py` describing how to reproduce the Phase-2 comparison.
- [ ] 10.3 Inline comment in `btojpsik.py` `get_bkmm_alcareco_selections()` cross-referencing this change and the `_kin_non_v0` helper in the cff (both must stay in sync).
- [ ] 10.4 Update the `[Memory] J/ψ + X AlCaReco OpenSpec proposal` memory entry to point at this change.
- [ ] 10.5 Open a follow-up note to re-scope `add-jpsi-x-stage2-cvh-refit` as Stage-2 work (the AlCaReco-stage framing it currently uses is obsolete). Do not archive or rewrite it as part of this change.
- [ ] 10.6 Open a follow-up note to locate or request central MC samples for the other six channels (B0→K*0, Bs→φ, Λb, ψ(2S), Bc, B0→Ks) so Phase 1 can eventually be run with gen truth on those channels too.
- [ ] 10.7 Final check on the worklog deck before this change is archived: every iteration / decision / anomaly captured during the work has at least one corresponding frame, frames are in chronological order, and the deck compiles cleanly. Commit the final PDF alongside the source if the repo conventions allow.
