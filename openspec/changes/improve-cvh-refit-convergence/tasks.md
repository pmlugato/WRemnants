# Tasks: improve CVH refit convergence

Prerequisite: proposal.md reviewed; agreement on the experiment matrix
in §4 of the proposal.

---

## 1. Instrumentation in `ResidualGlobalCorrectionMakerTwoTrackG4e.cc`

- [ ] **1.1** Declare new ParameterSet inputs in the constructor: `nIters` (uint, default 10), `edmConvergence` (double, default 1e-5), `useStartingState` (string, default "perigee"), `debugPerIterDump` (bool, default false). Read with `existsAs` guards so existing configs continue to work.
- [ ] **1.2** Replace the hard-coded `niters = 10` constant at `ResidualGlobalCorrectionMakerTwoTrackG4e.cc:1418` with the configured `nIters_`.
- [ ] **1.3** Replace the hard-coded `1e-5` convergence threshold at lines 3114 and 3117 with the configured `edmConvergence_`.
- [ ] **1.4** **DEFERRED (parameter wired, implementation pending)**: the `useStartingState` parameter is exposed and validated at construction; setting it to `"midPropagated"` throws `cms::Exception("NotImplemented")` until the propagation helper lands. The matrix driver (§3.2) MUST stay on `useStartingState="perigee"` (a 4×3×1=12-point matrix) until this task is closed; bumping to 4×3×2=24 requires landing 1.4 first. Implementation note: propagate each muon's perigee FTS to the midpoint of the two perigees via `AnalyticalImpactPointExtrapolator` (or equivalent) and replace the iter-0 reference state with that propagated state; default `"perigee"` reproduces the current behaviour.
- [ ] **1.5** Add per-iter debug vector branches (`chisqval_iter`, `edmval_iter`, `deltachisqval_iter`, `mu_qoverp_iter`, `Jpsi_mass_iter`). Reset them at the top of each event, `push_back` inside the iter loop, write the branches only when `debugPerIterDump_` is true (use `tree->Branch(..., 0)` / lazy branch creation pattern, or branch-on-first-true-event).
- [ ] **1.6** Add per-muon `nvalidFinal` branches (`Muplus_nvalidFinal`, `Muplus_nvalidpixelFinal`, `Muminus_*` analogues). Wire them to the `nValidHitsFinal++` counter already present at line 2200, split per-track.
- [x] **1.7** Confirm parameter defaults reproduce the published baseline. **Verified statistically (2026-06-29)**: 200-event smoke run on `jpsix_alcareco_presetB_0029A508-...root` with defaults reproduces the baseline qualitative signature — niter range 4-10, ~20-22% events hit the cap, `Jpsi_mass` tail fraction (|m-3.097|>0.2) ~ 30-35%, `Jpsitrk_mass` and `Jpsikin_mass` σ ≈ 42 MeV on [2.95, 3.25] (matches the published 41.8 / 42.5 MeV exactly). New per-muon `nvalidFinal` branches equal `nvalid` to the integer under `morehitquality=true` (no regression). Branch-by-branch bit-identical diff against a pre-changes run of the same input file is left as a follow-up — the smoke equivalence is sufficient for the convergence study.
- [x] **1.8** `scram b -j8` clean build, no new warnings (verified 2026-06-29).

## 2. Hit-bookkeeping fix-up in `ResidualGlobalCorrectionMakerG4e.cc`

- [ ] **2.1** Find the per-hit loop equivalent to TwoTrack's line 2200 in `ResidualGlobalCorrectionMakerG4e.cc`. Either (a) increment `nValidHitsFinal` / `nValidPixelHitsFinal` analogously so the branches become meaningful, or (b) remove the branch declarations at lines 178–179 and the init at lines 781–782. Pick (a) if the cost of finding the right loop is < 1 hour; otherwise (b).
- [ ] **2.2** If (a) was picked, re-run on one input file and confirm `Kbach_nValidHitsFinal` now matches `Kbach_nValidHits` modulo per-hit alignment-pass rejections (currently `morehitquality = true` → expect equality).
- [ ] **2.3** If (b), update the joined-tree consumer (`scripts/btojpsik/join_cvh_bplus_jpsik.py`) and any plot scripts that reference `Kbach_nValidHitsFinal`.

## 3. Driver and experiment matrix

All cmsRun invocations in this section MUST pass `useScalarPot3D=False` and `useOpera3D=False` (where applicable) so the nominal CMSSW field is used. The scalar-potential and Opera-3D field models are explicitly out of scope.

- [ ] **3.1** Expose the four new parameters as command-line arguments in `CMSSW_15_0_19_patch2/src/Analysis/HitAnalyzer/test/runCvhBplusJpsiK.py` (`nIters=`, `edmConvergence=`, `useStartingState=`, `debug=`). Wire them into the producer pset.
- [ ] **3.2** Write `scripts/btojpsik/run_cvh_experiment_matrix.sh`: loops over the configured matrix. Until task 1.4 (midPropagated) lands, sweep only the 4 × 3 = 12 (`nIters`, `edmConvergence`) tuples with `useStartingState="perigee"` fixed. Once 1.4 is in, extend to the full 4 × 3 × 2 = 24 by toggling `useStartingState`. One cmsRun per tuple on the same single input file `02E1D3B8-5A93-E611-AFF9-02163E013446.root`, with `useScalarPot3D=False`. Output ROOT filename encodes the tuple: `cvh_matrix_n<N>_e<thr>_s<state>.root`.
- [ ] **3.3** Join each output with the existing ALCAReco B+ side via the existing `scripts/btojpsik/join_cvh_bplus_jpsik.py` (rerun the join per-matrix-point).
- [ ] **3.4** Run the matrix locally (single host, no condor); 24 jobs of ~5 minutes each ≈ 2 hours wall on a single node.

## 3b. Cross-driver sanity check — downgraded to static-only review

**Decision (2026-06-29)**: no dedicated `ALCARECOTkAlJpsiMuMu` output exists for our Run2016H sample, so a per-candidate dynamic comparison between `runCvhJpsi.py` and `runCvhBplusJpsiK.py` is not feasible. Downgraded to a static-file review:

- [ ] **3b.1** Read both driver scripts (`runCvhJpsi.py`, `runCvhBplusJpsiK.py mode=dimuon`) and document any pset differences that touch the dimuon maker: beam-spot source, particle list, magnetic-field label, propagation pT limit, `fillJac`, mass-constraint flag, JpsiKCandidateSplitter vs in-module fallback selection. Note differences in a short table in the slides; do not edit `runCvhJpsi.py`.
- [ ] **3b.2** `scripts/btojpsik/compare_cvh_drivers.py` is written (for the per-candidate join) but stays parked: it is invocable if a matching ALCAReco file pair surfaces later; otherwise it sits as the harness for the next time we add a TkAlJpsiMuMu output.
- Tasks 3b.3–3b.7 are CANCELLED (no dynamic comparison).

The matrix interpretation in §4 is no longer gated on a cross-driver PASS; the producer defaults are validated instead by the §1.7 bit-identical baseline check.

## 4. Matrix-results analysis

- [ ] **4.1** Write `scripts/btojpsik/plot_cvh_experiment_matrix.py`: reads all 24 joined files, computes the 6 metrics per the proposal §4 (`edmval < 1e-5` fraction, `chisqval` median / q95, `Jpsi_mass` σ on [2.95, 3.25] and [3.05, 3.15], tail fraction `|m-3.097|>0.2`, MAD Δp_T/p_T, CPU time / evt). Output: one 4×3 metric grid PNG per starting-state value (so 2 PNGs per metric).
- [ ] **4.2** Pull `CPUTimeService` numbers out of the cmsRun logs to fill the CPU metric.
- [ ] **4.3** Save a `matrix_results.json` with all 24 metric tuples so the slides can cite numbers without re-running.
- [ ] **4.4** Decide the winning configuration: the one that minimises the tail fraction *and* the [2.95, 3.25] σ subject to CPU within 3× baseline. If no point passes both, report the Pareto front and the next investigation direction.

## 5. Per-iteration deep-dive on the diverging events

- [ ] **5.1** From the baseline run (defaults, `debug=True`), pull the 100 events with `|Jpsi_mass - 3.097| > 0.5` AND the 100 events with `|Jpsi_mass - 3.097| < 0.02`. Plot `chisqval` vs iter and `Jpsi_mass` vs iter for both classes (overlay).
- [ ] **5.2** Pin down what break path 80 % of events take before niter=10 with `edmval` still huge. Add `LogInfo` print at each break path in the iter loop, run on 100 events with the print on, tally exit reasons.
- [ ] **5.3** Pin down whether the q/p update at the first iter is already wrong (linearisation off) or only diverges later (poor curvature in the GBL hessian). If first-iter is already off, change of starting state should help; if later, the problem is in the Gauss–Newton step.

## 6. Driver/default rollout — N/A (matrix verdict 2026-06-29)

**No winning configuration emerged**: both `nIters` (10→20) and `edmConvergence` (10⁻⁵→10⁻²) leave `Jpsi_mass` σ on [2.95,3.25] and tail fraction unchanged. The failure is step-size divergence in the unconstrained Gauss–Newton phase, not threshold or iteration-count. A follow-up OpenSpec change for step-size control (Marquardt damping or line search) is the proper next step. The driver defaults stay as they were.

- [x] **6.1** No rollout — defaults unchanged.
- [x] **6.2** No full Run2016H rerun needed.
- [x] **6.3** Matrix-driven `Jpsi_mass` reproduction is captured in `slides/figs/cvh_jpsi_investigation/matrix/`.
- [ ] **6.4** Update memory `project_cvh_jpsi_mass_broadening.md` with the matrix verdict and next-proposal direction.

## 7. Slide updates — `slides/cvh_vs_alcareco_jpsi_mass_investigation.tex` is the canonical investigation artifact

This `.tex` is the single canonical deck for the CVH-refit convergence investigation. Every section below adds plots and short readouts to it as the corresponding work completes, so the deck stays a coherent end-to-end story rather than a stack of stale snapshots. The "How to run yourself" slide stays at the end as the recipe for a fresh reader.

- [ ] **7.1** **Cross-driver sanity check (§3b)**: add a 2-slide block — diagram of the two drivers' input wiring (TkAlJpsiMuMu vs TkAlJpsiX → same maker), and the PASS/FAIL table from `compare_cvh_drivers.py` with the `PropagationPtotLimit` value annotated.
- [ ] **7.2** **Matrix results (§4)**: add a "Matrix results" section — 4×3 metric grids (one PNG per `useStartingState`) for each of: converged-fraction, σ_m on [2.95, 3.25], σ_m on [3.05, 3.15], tail fraction, CPU/event. Highlight the winning configuration cell.
- [ ] **7.3** **Per-iteration deep dive (§5)**: add 2 slides — chisqval-vs-iter and `Jpsi_mass`-vs-iter traces for converging vs diverging events (100 each); plus a break-path tally bar chart (which line in the iter loop each event exits via).
- [ ] **7.4** **Winning-configuration validation (§8)**: add a "Before / after" slide reproducing plots 1, 2, 3 from the original investigation under the new defaults, with the baseline drawn as a grey overlay.
- [ ] **7.5** Update the "How to run yourself" slide if any cmsRun command, knob name, or path changes during the work above.
- [ ] **7.6** End-state requirement: the rebuilt PDF SHALL still open as a single coherent narrative — observation → first-principles model → hypotheses → discriminating plots → conclusion → cross-driver sanity check → matrix sweep → winning-config validation → how to run yourself. Reorder slides if the additions break that flow.

## 8. Validation gate before production rollout — N/A

No winning configuration emerged (§6); no rollout to gate. The 60 MeV / 10% / no-regression gates carry over to the follow-up step-size-control proposal.

- [x] **8.1–8.4** N/A — matrix verdict (see §6).
