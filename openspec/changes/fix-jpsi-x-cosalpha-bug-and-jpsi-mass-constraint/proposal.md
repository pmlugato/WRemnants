# Change: Fix the cosα bug and add a J/ψ mass constraint in `JpsiXCandidateProducer`

## Why

Two producer-level defects were discovered in `add-jpsi-x-condor-production`'s 5-file Run2016H run + a one-file diagnostic re-run (Cdiag, cluster 3107800):

**Defect 1 — cosα cut is non-physical.** `JpsiXCandidateProducer::computeAlphaBSFromPoint()` uses the beamspot (LHC bunch center) as the origin for the displacement direction. The bunch z-spread is ~5–15 cm; the PV-vs-BS z-offset dominates `r = sqrt(dvx² + dvy² + dvz²)` (B+ flight is ~500 μm). The resulting cosα ≈ `sign(pz · dvz)` is essentially random vs the B's true pointing. Cutting at cosα > 0.95 rejects candidates pseudo-randomly with a slight anti-bias against signal (real B+ have a different pz distribution than combinatorial). Diagnostic data: with cosα + the two other Kalman cuts enabled, preset C kills 99.3 % of real B+ signal (signal eff 0.6 % vs preset B's 31.6 %); with the same cuts at no-op (Cdiag), Kalman fit convergence is 88–96 % per channel — the fit machinery is healthy, only the cuts on top of it are broken.

**Defect 2 — no J/ψ mass constraint.** The B candidate's 4-momentum is computed as `jpsi.p4() + bachelor.p4()` without constraining `m(μμ) = m_J/ψ`. The measured μμ mass resolution (~10–15 MeV) propagates directly into the B+ mass, broadening the peak from the ~15 MeV CMS detector floor to ~30+ MeV observed in our data. This is why neither preset B nor preset C produced a visible mass peak in the 5-file production: the signal IS there (923–1 447 sideband-subtracted B+ candidates in preset B, consistent with ~30 % efficiency × 4 580 expected), but smeared 2× too broad to stand above combinatorial.

Two smaller logic issues compound the problem:

- `applyBLevelVertexFit()` silently skips the `Lxy/σ` cut when `lxy < 1e-10` (line 327), creating a hole for prompt-vertex candidates.
- The PV used for Lxy and (post-fix) cosα is `pvs->front()` — the highest-sumpT² PV — which may not be the actual production vertex of the B candidate.

## What changes

The producer is patched in four places, all in `JpsiXCandidateProducer.cc`. No cut value changes; no cff changes; no preset re-design (modulo dropping preset A from active testing — see scope below). The cuts already specified in `add-jpsi-x-vertex-fit-and-low-pt` (iter-1: `vtxProb > 0.01`, `cosα > 0.95`, `Lxy/σ > 1.0`) are left in place and re-evaluated against the *corrected* underlying variables.

### 1. Replace `computeAlphaBSFromPoint(... bs)` → `computeAlphaPV(... pv)`

The cosα-from-beamspot helper is removed. A new helper takes the matched PV's position as the origin. The implementation uses the **2D transverse** angle by default — best practice in B physics analyses, robust to PV z-resolution, matches what published CMS B+ → J/ψ K+ measurements use. The call site in `applyBLevelVertexFit()` is updated accordingly.

### 2. Fix the `lxy < 1e-10` silent-pass

If the Kalman-fitted B vertex coincides with the PV to within 1 μm (lxy < 1e-10 cm), the candidate fails the Lxy/σ cut (`return false`), it doesn't pass silently. Prompt-like fits should not satisfy a "displaced" requirement.

### 3. Use the closest-in-z PV instead of `pvs->front()`

A small helper picks the PV with `|vz_B − vz_PV|` minimum among the offline PV collection (with a fallback to `pvs->front()` if no PV is within ~1 cm in z). This is the standard B-physics PV-matching convention.

### 4. Apply a J/ψ mass constraint in the B-candidate 4-momentum (preset B and C only)

The J/ψ mass constraint is *independent* of the B-level vertex fit. Two code paths, chosen by `applyBVtxFit_`:

- **Preset B path** (no B-vertex fit; `applyBVtxFit_ = false`): refit *only* the dimuon (μ⁺, μ⁻) with `KinematicParticleFitter` + `MassKinematicConstraint(m_J/ψ = 3.0969 GeV)`. Use the resulting constrained J/ψ 4-momentum in `lvM = jpsi_constrained.p4() + lvBach`. Vertex stays as the existing midpoint (no B-vertex fit performed). Effective cost: one fit on 2 tracks.
- **Preset C path** (B-vertex fit ON; `applyBVtxFit_ = true`): one `KinematicConstrainedVertexFitter` over all leaf tracks (μ⁺, μ⁻, bachelor) — or (μ⁺, μ⁻, K, π) etc. for VCC mode — with the same J/ψ mass constraint applied to the muon pair. Returns *both* the B vertex position and the constrained B 4-momentum in one call. Replaces the existing `KalmanVertexFitter` invocation.

On fit failure (expected rate ~1 %), fall back to the existing unconstrained `jpsi.p4() + bachelor.p4()` + midpoint vertex, without dropping the candidate. Fallback counters get logged in the producer's summary.

Expected: B+ peak σ ~30 MeV → ~15 MeV under both presets; this is essential to making the peak *visible* on real data.

### 5. Re-run validation locally (no Condor for the iteration cycle)

While iterating on the producer fixes, the validation runs locally on the submit machine (interactive `cmsRun` via `cmssw-el7` Singularity), 10 000 events from the smoke file (run 283270). Cycle time: edit code → `scram b -j 8` → `cmsRun` 10 k events (~25–45 min CPU) → diagnose → repeat. This is the fastest iteration loop while the code changes are still fluid.

Condor scale-up (5 files × 2 presets) is **out of scope** for this change — once the local iteration converges on a correct-looking producer, a sibling scale-up change runs the production. The existing Condor wrapper from `add-jpsi-x-condor-production` is reused without modification at that point.

### 6. Update OpenSpec results + slides at scale-up time

After the local iteration converges, the sibling scale-up change owns the new `results.md` and slide updates. This change only writes a brief `results.md` summarizing the local-validation efficiency numbers + before/after mass-shape comparison.

## Scope and explicit non-scope

**In scope.**

- The four producer-code fixes above (cosα, Lxy bypass, PV matching, J/ψ constraint).
- Local validation runs (1-file × {preset B, preset C}) sufficient to confirm the fixes work and produce sensible numbers.
- A short `results.md` documenting the local-validation outcome.
- The deprecation note already added to `add-jpsi-x-condor-production/results.md`.

**Out of scope.**

- **Preset A.** No active testing under preset A in this change. The J/ψ constraint applies under preset A in the *code* (it's a producer-level correctness fix), but no preset-A run is performed — preset A's defining purpose was the unselected-baseline shape, which is no longer a useful comparison point given the broader analysis trajectory.
- **Bachelor track quality cuts** (highPurity, pixel hits, normalized χ²). These would pre-filter the very tracks alignment wants to study — they bias the alignment sample toward already-clean tracks, defeating the use case. Explicitly *not* added in this change. Future-work note only.
- **Bachelor PID** (dE/dx-based K/π separation). Same reason as bachelor quality cuts — biases the alignment sample. Out of scope.
- **Iter-2 cut-value tuning.** Iter-1 values stay in place; if the fixed numbers warrant tuning, that's a sibling change after this one's data is in.
- **Condor scale-up.** Local iteration only for this change; scale-up is a sibling.
- **Multi-era expansion** (2016B-G, 2017, 2018).

## Impact

- Touched files: `JpsiXCandidateProducer.cc` (~120 added/changed lines), `BuildFile.xml` (one new dep: `RecoVertex/KinematicFit`), this change's docs.
- Producer recompile required (`scram b -j 8` in `CMSSW_10_6_17_patch1`).
- Cluster cost: **zero** in this change. Local CPU: ~25–45 min per validation cycle; expect 2–4 cycles to convergence.
- Output volume: a handful of local `.root` files (~30 MB each), kept under `condor/jpsix_alcareco/_local_validation/` so they're parallel to (not mixed with) the production layout.

## Open decisions

Decisions committed in the draft, but you may want to revisit:

1. **cosα flavor: 2D transverse vs 3D-from-PV.** Committed to 2D — simpler, matches published B+ → J/ψ K+ analyses, robust against PV z-resolution. 3D-from-PV uses more info but is more sensitive to PV z-residuals; the gain over 2D is small in practice.
2. **J/ψ constraint mechanism: preset B uses dimuon-only refit; preset C uses full `KinematicConstrainedVertexFitter`.** Committed as described above — keeps preset B identity (no B-level vertex fit) intact while still gaining the resolution improvement, and unifies preset C's vertex fit with the mass constraint.
3. **Iter-2 cut tuning: separate change after this.** Iter-1 cut values are unchanged here; new numbers from this change inform whether iter-2 is needed.
4. **Local validation only.** No Condor in this change; sibling scale-up after iteration converges.
