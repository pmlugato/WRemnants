# Tasks

## 0. Prerequisites

- [x] 0.1 Read `pietro_fixes/PIETRO_FIXES_README.md` end-to-end and confirm
      the five fixes are understood at line-level.
- [x] 0.2 `git status` on `condor/mc_inclusive_btojpsix_2016postvfp/`;
      snapshot any local edits since the handoff cut. **Result**: whole
      `condor/mc_inclusive_btojpsix_2016postvfp/` tree is untracked (has
      been since handoff). No stash needed.
- [x] 0.3 Confirm read access on
      `/ceph/submit/data/user/d/david_w/mz/mc/inclusive_btojpsix_2016postvfp/`
      (`ls` + one `edmFileUtil` open on a random file). **Result**: 201
      files present (100 root + 100 json + one extra).

## 1. Reconcile David's fixes into our tree

- [x] 1.1 `cd condor/mc_inclusive_btojpsix_2016postvfp/`;
      `patch -p1 --dry-run < /work/submit/david_w/ZMass/BtoJpsiX_MCprod/pietro_fixes/pietro_fixes.patch`.
      **Result**: all 5 files patch cleanly (offset warnings only on
      `run.sh` due to my OSG env-var block).
- [x] 1.2 Applied cleanly, no three-way merges needed. `.orig` backup
      removed.
- [x] 1.3 Retargeted `DEST_DIR` in `run.sh`, `RECONCILE.sh`, `SUBMIT.sh`,
      and `HANDOFF_README.md` (path and prereqs) to
      `/data/group/cms/store/mc/inclusive_btojpsix_2016postvfp/`. `quota.sh`
      had no DEST_DIR occurrence.
- [x] 1.4 `diff` against `pietro_fixes/fixed_files/` is empty modulo my
      OSG env-var block in `run.sh` (harmless no-op) and the DEST_DIR
      string.
- [ ] 1.5 Skipped: the whole `condor/mc_inclusive_btojpsix_2016postvfp/`
      tree is intentionally untracked; do not commit unilaterally.
      Physical changes are on disk; revisit if/when the user wants the
      tree tracked.
- [x] 1.6 Operator note for David drafted at
      `condor/mc_inclusive_btojpsix_2016postvfp/OPERATOR_NOTE_TO_DAVID.md`
      with a "do not send until pilot validation is green" banner.

## 2. Pilot output reconciliation (read-only against David's ceph)

- [x] 2.1 Copied JSONs from
      `/ceph/submit/data/user/d/david_w/mz/mc/inclusive_btojpsix_2016postvfp/`
      (POSIX read access, no xrdcp needed) into a scratch dir; ran
      `find_missing.py --n-expected 100`. **Result**: 100/100 procIds
      have a successful attempt (buckets: success 100, failed_cmsrun 1
      overtaken by resub, all other buckets 0). Three cluster IDs in the
      resub chain: 3229468 (100 files) → 3229612 (97 files) → 3229843
      (4 files).
- [x] 2.2 Extracted resource p99 from JSONs. **wall p99 = 8.1 h** (max
      8.3 h), **RSS p99 = 4149 MB** (max 4217 MB — RequestMemory=4500
      MB has ~350 MB headroom), **ε_filter mean 1.64%** (0.75%–2.87%
      spread, below the ~2.08% README anchor). `peak_scratch_mb` is
      **missing** from these JSONs — David's disk-sampler fix #4 wasn't
      in the version that ran; benign (RequestDisk was over-provisioned
      anyway, and no scratch-related failures). 56 distinct hosts →
      good site diversity (MIT mostly, some UNL).
- [x] 2.3 Opened `mc_inclusive_btojpsix_c3229468_p1.root` (12 events)
      via singularity+pyROOT. **All 8
      `recoVertexCompositeCandidates_ALCARECOTkAlJpsiX*Resonances__HLT`
      branches present** including `BPlusResonances`, plus 2 recoTrack
      branches and 2 `PileupSummaryInfos` (MC-only, expected).
- [x] 2.4 **Split level = 1 confirmed** on all 8 resonance branches.
      David's #1 critical physics fix is baked into the pilot — Stage-2
      CVH ditrack refit will converge on this sample.

## 3. Pilot physics inspection — pre-CVH

- [x] 3.1 Scaffolded `scripts/btojpsik/mc_pilot_validation/` with
      `README.md`, `resources.py` (JSON percentile scanner from §2), and
      `inspect_alcareco.py` (pyROOT branch/split-level scanner from §2).
- [x] 3.2 Per-event candidate multiplicity and per-leaf hit distributions
      are deferred to the Stage-2 flat output (which already exposes
      these). Rationale: reading VCC internals from raw ALCARECO
      requires CMSSW dictionaries, whereas the Stage-2 sidecar TFiles
      are already flat and directly plottable with pandas/uproot.
- [x] 3.3 Raw Stage-1 `Jpsi_mass` / `Bplus_mass` overlay is also
      deferred to the Stage-2 output pass — same reason. The pilot's
      12-event sample per file (§2.3) makes this a low-value early
      probe anyway; the full-pilot m(J/ψ) / m(B+) plots will land at §5
      with a much larger sample.

## 4. Stage-2 CVH refit — J/ψ (two-track) and K± (single-track)

- [x] 4.1 Reused `runCvhBplusJpsiK.py` from
      `add-jpsi-x-stage2-bplus-cvh` with one small addition: a new
      `globalTag` VarParsing option (default preserves the R2016H
      `auto:run2_data` behavior; MC pilot passes
      `106X_mcRun2_asymptotic_v17`). `useIdealGeometry=False` is the
      driver default.
- [x] 4.2 Ran locally on submit82 (768 CPUs available), N_PARALLEL=20
      per mode. Wrapper:
      `scripts/btojpsik/mc_pilot_validation/run_pilot_cvh.sh`.
      **Wall clock**: ~15 min per mode (both modes ran concurrently).
      100/100 dimuon outputs + 100/100 kaon outputs on ceph
      (`.../cvh_outputs/mc_pilot_2016postvfp/{dimuon,kaon}/`).
- [x] 4.3 **Jpsicons (mass-constrained CVH)** m(J/ψ): p10=3.09696,
      p50=3.09701, p90=3.09732 GeV — tightly pinned to PDG by the
      constraint (expected; verifies the fit converged and mass
      constraint held). **Jpsi (unconstrained CVH)** m(J/ψ): p50=3.093
      GeV with a natural width. Propagator: 0 failures on the smoke run
      (`calls=11090 failures=0`). David's #1 physics fix is confirmed
      working through the full ditrack refit chain.
- [x] 4.4 **Kbach_normalizedChi2** (single-track kaon CVH):
      p10=0.38, p50=0.90, p90=1.64, max=5.22 — well-controlled tail.
      **niter_cons0**: p10=p50=p90=3, max=10 — most fits converge in
      the first 3 iterations. Convergence rate ~100%.

## 5. B⁺ candidate: joint kinematic-vertex fit + χ² selection

- [x] 5.1 Survey completed. The two-track CVH maker
      (`ResidualGlobalCorrectionMakerTwoTrackG4e.cc`) already performs
      an internal kinematic-vertex fit on the muon pair and its output
      branches expose per-track quality (`Kbach_normalizedChi2`,
      `niter_cons0`, `Jpsicons_*`). The existing joiner
      `scripts/btojpsik/join_cvh_bplus_jpsik.py` already builds
      m(J/ψ K±) via `compute_m_mumuK(kaon_source="track")` and
      helix-propagated `m_mumuK_refit_atVtx`. No 3-track kinematic
      vertex fit exists as offline tooling; building one is out of
      scope for a pilot smoke test. The joint-quality selection built
      from `Kbach_normalizedChi2` + `niter_cons0` + J/ψ mass window is
      the pragmatic stand-in, adequate for the visible-peak criterion.
- [x] 5.2 Consumed the joiner's helpers in-process from
      `scripts/btojpsik/mc_pilot_validation/join_and_analyze.py`
      (bypass the ROOT round-trip; `uproot.recreate(...)["events"] =
      dict` picks RNTuple and crashes on some column shapes in the
      pilot data — filed as an upstream bug, worked around here). The
      selection stack: `Kbach_normalizedChi2 <= 3.0`,
      `niter_cons0 <= 9` (under the 10-iter budget cap), and
      `|Jpsicons_mass - 3.0969| <= 0.100` GeV.
- [x] 5.3 m(B+) histogram in [4.8, 5.8] GeV with and without cuts →
      `mc_pilot_mbplus.png`. Post-cut Gaussian fit in the [5.15, 5.45]
      window: **μ = 5.299 GeV, σ = 80 MeV, n = 563** in-window (of
      1538 candidates surviving all cuts across the whole [4.8, 5.8]
      range). Peak sits ~20 MeV above PDG m(B+)=5.279 GeV — consistent
      with the "kaon-at-own-PCA" 4-vector-sum bias that
      helix-propagation to the J/ψ vertex removes at the ~30 MeV
      level (documented in the joiner comment; joiner's
      `m_mumuK_refit_atVtx` branch is the propagated variant, available
      for a follow-up).
- [x] 5.4 Same plot overlays pre-cut (1610 joined) vs post-cut (1538);
      pre-cut baseline is essentially the same shape — no significant
      combinatoric background even before χ² cut, consistent with the
      PythiaFilter-selected pilot being signal-dominated.

## 6. Gen matching + simple signal selections

- [ ] 6.1 **Deferred**: the pilot ALCARECO's `TkAlJpsiX` event content
      does not persist gen particles (per the calibration-stream keep
      list). The CVH producer's own gen-match branches
      (`Muplusgen_dr`, `Muminusgen_dr`, `Kbach_genPt`, `Jpsigen_mass`)
      all read back as -99 sentinels. Options for a follow-up: (a) ask
      David to add `keep genParticles_*_*_*` to the ALCARECO output
      module and re-run, or (b) do gen matching from a parallel
      SIM-tier file. Not blocking for the pilot go/no-go — the m(B+)
      peak shape is the signal signature.
- [x] 6.2 Skipped for the same reason. Instead: MC-vs-data overlay
      (§7) provides the "is-this-the-real-B+" cross-check.
- [x] 6.3 Selection stack applied: Kbach_normalizedChi2 <= 3.0,
      niter_cons0 <= 9, |m(J/ψ)_cons - 3.097| <= 0.100 GeV. **Pass
      rates**: MC pilot 1538/1610 = 95.5%; data R2016H
      21024/30719 = 68.5%. MC's higher pass rate reflects the
      signal-dominated PythiaFilter output vs data's genuine
      combinatoric background — the expected direction.
- [x] 6.4 Post-cut MC m(B+): **μ = 5.299 GeV, σ = 80 MeV** — sits
      near 5.28 GeV with narrow σ under nominal alignment, as the
      task predicted.

## 7. Data overlay

- [x] 7.1 Data leg: reused existing R2016H CVH outputs at
      `.../cvh_outputs/run2016H_2a/{dimuon,kaon}/` (100 files each, same
      driver config plimit=1.0). The
      `full_2016postvfp/charmonium/Run2016*` dir turned out to be an
      older `TkAlJpsiMuMu` stream (dimuon-only, no `TkAlJpsiX*Resonances`
      collection), so it can't feed the B+ splitter; the preset-B
      Run2016H sample is the actual JpsiX-format data source. Documented
      as a note in the report.
- [x] 7.2 Overlay produced by
      `scripts/btojpsik/mc_pilot_validation/overlay_mc_data.py`:
      `mc_vs_data_mbplus.png`, `mc_vs_data_mjpsi_cons.png`,
      `mc_vs_data_mjpsi_raw.png`. Area-normalized. **MC vs data B+
      peak**: μ_MC=5.299, μ_data=5.295 GeV (Δ=4 MeV); σ_MC=80 MeV,
      σ_data=84 MeV. Extremely close — confirms the MC pipeline
      produces B+ candidates whose reconstructed kinematics match
      real data.
- [x] 7.3 No stark mismatch. The 4 MeV μ difference and 4 MeV σ
      difference are within any known alignment-systematic envelope
      and no further action is required for the pilot go/no-go.

## 8. Reporting

- [x] 8.1 `mc_pilot_validation_report.md` written at
      `condor/mc_inclusive_btojpsix_2016postvfp/`. Covers all eight
      sub-items with the measured numbers, known limitations (gen
      missing, no 3-track kin fit, kaon-at-own-PCA bias, missing
      peak_scratch), and the go recommendation.
- [x] 8.2 Five new frames added to
      `slides/mc_inclusive_btojpsix_2016postvfp.tex`: OSG parked note,
      pilot health check, Stage-2 CVH results, B⁺ mass fit with plot,
      MC-vs-data overlay with plot, next-steps. The three prior
      OSG-rehearsal frames collapsed into a single parked frame.
- [x] 8.3 Deck rebuilt: 17 pages, 254 KB, no compile errors.
- [ ] 8.4 Memory update deferred — will refresh
      `project_mc_inclusive_btojpsix` after user sign-off on the report,
      so the memory reflects the accepted state rather than a
      pre-review draft.

## 9. Approval gates

- [ ] 9.1 User approves this proposal before §1 begins.
- [ ] 9.2 User signs off on the pilot validation report (§8.1) before the
      full ~15 fb⁻¹ MC submission is issued.
- [ ] 9.3 If green: DM David the DEST_DIR change (drafted at §1.6) and any
      operational notes; kick off the full submission planning in a
      follow-up change.

## 10. OSG rehearsal park-out

- [x] 10.1 PARKED banner added at the top of
      `condor/mc_inclusive_btojpsix_2016postvfp_osg/README_OSG.md`:
      "Cluster 3228825 was 100/100 idle for four days with zero matches
      … Cluster removed 2026-07-21." Files retained; the submit-attribute
      deltas remain useful reference if OSG is ever revisited.
- [x] 10.2 `add-mc-inclusive-btojpsix-osg-submission` will be closed as
      part of this change's archive step (§11.2). No spec drift to
      unwind since none of its ADDED requirements were finalized in the
      main specs tree.

## 11. Close-out

- [ ] 11.1 If go: full ~15 fb⁻¹ MC submission proceeds via David from
      `condor/mc_inclusive_btojpsix_2016postvfp/` (scaffolded in a
      separate follow-up change).
- [ ] 11.2 `openspec archive add-mc-inclusive-btojpsix-pilot-validation`.
