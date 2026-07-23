# Tasks: fix-jpsi-x-cosalpha-bug-and-jpsi-mass-constraint

## 1. Producer code — cosα + Lxy + PV matching (priority 1)

- [x] 1.1 Removed `computeAlphaBSFromPoint()`. Added `computeAlphaPV()`: 2D transverse cosα, matched PV as origin.
- [x] 1.2 Added `bestPVForCandidate()`: closest-in-z PV, fallback to `pvs.front()` beyond 1 cm dz with `edm::LogWarning`.
- [x] 1.3 In `applyBLevelVertexFit()`, matched PV resolved once and reused for Lxy AND cosα. Also extended the constructor's PV-token consume and `produce()`'s PV fetch to fire on `applyBLxyOverSigma_ || applyMotherAlphaBS_` (cosα now needs PV too).
- [x] 1.4 Removed the `lxy < 1e-10` silent-bypass. When `applyBLxyOverSigma_` is true and the fitted B vertex coincides with the matched PV to within machine precision, the candidate is rejected (no silent pass).
- [x] 1.5 Built with `scram b -j 8`. **Local smoke (preset C, run 283270, 10k events) shows preset C recovering from ~2.8 % of preset B (pre-fix) to ~28 % (task-1 fixed)** — bplus 139 → 1 421; bzero_kstar ~60 → 1 370; bs_phi ~6 → 146; bc ~236 → 1 755. V0 channels preset-invariant as expected.

## 2. Producer code — J/ψ mass constraint (priority 1)

- [x] 2.1 Added includes for KinematicParticleFitter, KinematicParticleVertexFitter, MassKinematicConstraint, KinematicConstrainedVertexFitter, TwoTrackMassKinematicConstraint, and KinematicParticleFactoryFromTransientTrack + the primitives headers (`RefCountedKinematicTree`, `RefCountedKinematicParticle`, `RefCountedKinematicVertex`).
- [x] 2.2 Added `constrainJpsi4Momentum()`: dimuon-only `KinematicParticleFitter` + `MassKinematicConstraint(m_J/ψ)`. Returns constrained J/ψ 4-momentum, or false on failure. Used under preset B.
- [x] 2.3 Added `constrainedBVertexFit()`: multi-track `KinematicConstrainedVertexFitter` with `TwoTrackMassKinematicConstraint(m_J/ψ)` applied to the muon pair. Returns 4-momentum, vertex, vertex error, χ², ndof. Used under preset C — replaces the previous `KalmanVertexFitter` invocation.
- [x] 2.4 Wired into `produceTrackMode()`. J/ψ constraint computed once per J/ψ (outside the bachelor loop) for efficiency. Multi-track fit runs inside the bachelor loop under preset C; on failure, fall back to J/ψ-constrained sum + midpoint vertex without dropping the candidate.
- [x] 2.5 Same wiring in `produceVccMode()`. VCC daughter masses read from each `xCand.daughter(i)->mass()` (correct per-channel K/π/p hypothesis).
- [x] 2.6 Removed `applyBLevelVertexFit()` and the `KalmanVertexFitter` include.
- [x] 2.7 Added `RecoVertex/KinematicFit` and `RecoVertex/KinematicFitPrimitives` to `BuildFile.xml`. Compile succeeded after two const-ref fixes (CMSSW's `KinematicParticleFactory::particle(...)` and `TwoTrackMassKinematicConstraint(...)` take non-const references by API; worked around with local non-const copies).
- [x] 2.8 Added `~JpsiXCandidateProducer()` destructor that emits an `edm::LogInfo` summary with `n_jpsi_constraint_attempted_/_fallback_` and `n_b_vertex_fit_attempted_/_fallback_` counters.

## 3. Local iteration loop

Cycle: edit code -> `scram b -j 8` -> `cmsRun config.py` on 1 file from EOS -> diagnose. Converged in 2 cycles (task-1 then task-2). Numbers below are post-both-fixes on the smoke file.

- [x] 3.1 Wrote `condor/jpsix_alcareco/_local_run.sh`. Uses an inner script in the scratch dir (passed to `cmssw-el7` Singularity via `bash inner.sh`) rather than nested `bash -c "..."` heredocs - those were fragile under the wrapper's own `bash -c` indirection.
- [x] 3.2 Preset B local: producer ran clean. **n_cands_bplus = 4904** (vs pre-fix 4946, -0.85 %); other channels within +-1 % of pre-fix. J/psi-constraint fallback rate well below 1 %.
- [x] 3.3 Re-extracted B+ signal. **Sig-window RMS 17.2 MeV** (was broad ~30 MeV), **mean 5.280 GeV** (PDG 5.279). Wide/tight ratio 1.56 (was 1.9). Wide signal eff 24.1 % (was 31.6 % - slightly lower as the constraint moves some events from the +-30 MeV wide window into the narrower core). Peak now clearly visible.
- [x] 3.4 Preset C local: producer ran clean. **n_cands_bplus = 1396** (vs pre-fix ~139 - a 10x recovery). C/B ratio across non-V0 channels: 22-38 %. V0 channels bit-identical to preset B.
- [x] 3.5 Preset C B+ signal: **sig-window RMS 17.5 MeV**, **mean 5.278 GeV**, wide signal eff 12.6 % (was 0.6 % pre-fix). Signal purity 45 % in wide window vs 32 % under preset B - preset C now has higher signal purity, exactly the alignment use case.
- [x] 3.6 Iteration converged. Both fixes validated. Preset C is ~4.5x slower per event than preset B (487 ms vs 108 ms), driven by the multi-track constrained vertex fit.

## 4. Validation report

- [x] 4.1 Wrote `results.md` with per-preset candidate counts, B+ signal sideband-subtracted, mass-shape PNG, and per-channel pre-vs-post table.
- [x] 4.2 `openspec validate` passes strict.

## 5. Condor scale-up (preset B only, 50k events)

User instruction (2026-06-15): after the local validation, scale up preset B to ~50 k events via Condor (5 files, identical to the deprecated Condor production layout). Do **NOT** scale up preset C — its 4.5x runtime makes scale-up disproportionate to its scientific value for this iteration.

- [x] 5.1 Submitted `condor_jpsix_alcareco_presetB_fix.sub` (5 jobs, cluster 3118303). Output route to `preset_B_fix/Run2016H/` (parallel to the deprecated `preset_B/` path) via the wrapper's new `DEST_SUFFIX="_fix"` 3rd arg.
- [x] 5.2 All 5 jobs exit 0 (workers UNL/Vanderbilt/Wisconsin, 37-78 min/job). Aggregated via `merge_and_report.py --preset-suffix _fix` (added the suffix CLI arg to support the parallel output layout). Per-channel sideband-subtracted eff via new `_diag_efficiencies_one.py`. Outputs at `condor/jpsix_alcareco/results_5files_fix.md`, plots dir `~/public_html/mz/alcareco/condor_5files_fix/`.
- [x] 5.3 Appended Phase-4 section to `results.md` (50k numbers + projections + alignment-goal alignment table + open questions). Worklog slides updated with two new frames covering scale-up numbers and projections.
