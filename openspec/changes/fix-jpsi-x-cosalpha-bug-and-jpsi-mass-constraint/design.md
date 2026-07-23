# Design: fix-jpsi-x-cosalpha-bug-and-jpsi-mass-constraint

## Goal recap

Fix three producer-level defects exposed by `add-jpsi-x-condor-production`'s 5-file run and Cdiag follow-up:

1. cosα cut origin is the beamspot (LHC bunch center), not the PV — randomizes the cut.
2. The Lxy/σ cut silently passes any candidate with `lxy < 1e-10 cm`.
3. The B-candidate mass is computed without a J/ψ mass constraint — peak σ ~30 MeV instead of ~15 MeV.

After the fixes, validate locally (1 file × 2 presets, no Condor) and compare per-channel efficiency + mass shape against the deprecated pre-fix run. If the iteration converges, file a sibling change for Condor scale-up.

## Core decisions

### 1. cosα: 2D transverse, not 3D-from-PV

The 2D transverse cosα uses `(px, py)` and `(dvx, dvy)` only — no z-component. Two reasons:

- **Robust against PV z-residuals.** Even with PV-matching by closest-z, the residual PV z-resolution (~10–50 μm) combined with B-vertex z-resolution (~50–100 μm) gives ~100 μm in z, while the B's pz contribution to its momentum is large. A small z-residual maps to a noticeable cosα-3D-from-PV degradation. The 2D version uses the transverse plane where PV-vs-B-vertex differences are dominated by *real* B flight (cτβγ projected on transverse).
- **Matches published CMS B+ → J/ψ K+ measurements.** Standard B-physics analysis machinery (e.g., `Bmm5`, BPH analyses) uses 2D pointing.

The 3D-from-PV variant is correct in principle but harder to tune. We commit to 2D unless data shows we need 3D.

### 2. PV matching: closest-in-z, fallback to `front()`

The PV collection in CMS is sorted by decreasing `sumPt²`. For B+ events with a "soft" parton-level interaction, the actual production PV may be the second or third entry, not the first. Closest-in-z matching is the standard choice — picks the PV whose z is closest to the B-vertex z, where the B-vertex z is what the vertex fit returned.

We fall back to `front()` if no PV is within 1 cm in z (rare; indicates the event has no nearby PV — should be discarded by upstream filters anyway, but the fallback prevents a hard crash).

### 3. J/ψ mass constraint — two code paths

The constraint is *independent* of the B-level vertex fit. The producer chooses the path by the same `applyBVtxFit_` flag that already gates the Kalman branch:

**Preset B path** (`applyBVtxFit_ = false`):

```
constrainJpsi4Momentum(muRefs, ttb)
  → builds 2 TransientTracks + 2 KinematicParticles for the muons
  → KinematicParticleFitter with MassKinematicConstraint(3.0969 GeV)
  → returns constrained J/ψ 4-momentum (or empty on failure)

lvM = jpsi_constrained.p4() + lvBach        // bachelor unchanged
vtxM = midpoint of (jpsi.vertex, bachelor.vertex)  // unchanged
```

Cost: one mass-constrained refit on 2 tracks. ~1 ms per candidate. **Preset B's selection identity (no B-level vertex fit) is preserved.**

**Preset C path** (`applyBVtxFit_ = true`):

```
constrainedBVertexFit(allRefs, muRefs, ttb)
  → builds TransientTracks for all leaf tracks (mu, mu, bach) or (mu, mu, K, pi)
  → builds KinematicParticles for each
  → KinematicConstrainedVertexFitter with MassKinematicConstraint(3.0969 GeV) on the (mu,mu) pair
  → returns:
      - constrained B 4-momentum
      - B vertex position + error
      - chi2 / ndof
      (or empty on failure)

lvM, vtxM, chi2/ndof come from the fit. Existing three cuts (vtxProb, cosα via the new PV-based 2D helper, Lxy/σ) are applied on the fit output.
```

Cost: one constrained-vertex fit on 3 or 4 tracks. ~5–10 ms per candidate.

The old `applyBLevelVertexFit()` method (using `KalmanVertexFitter` without mass constraint) is removed — the new `constrainedBVertexFit` subsumes both vertex fit AND mass constraint in one call.

**Fallback on failure** (either fit type): fall back to the existing unconstrained `jpsi.p4() + bachelor.p4()` + midpoint vertex. Log at debug. Don't drop the candidate. Maintain counters for the producer's summary log.

### 4. No bachelor track quality cuts

Explicitly NOT added. CMSSW `reco::TrackBase::highPurity`, `numberOfValidPixelHits`, and `normalizedChi2` are the standard B-physics bachelor quality knobs — but the alignment use case *requires* low-quality / low-pT / edge-case tracks for the tracker fit to gain information from. Pre-filtering on these knobs biases the alignment sample toward already-clean tracks, defeating the producer's purpose. If iter-2 ever shows we need bachelor-side discriminants for signal selection, they enter via a separate change with explicit alignment-impact justification.

### 5. Local iteration loop, no Condor in this change

Cycle: edit `JpsiXCandidateProducer.cc` → `scram b -j 8` → run `cmsRun` interactively on 10 000 events of the smoke file via the `cmssw-el7` Singularity wrapper → diagnose with `_diag_bplus_shape.py`. Expected per-cycle: ~30 min wall. 2–4 cycles to convergence.

The existing Condor wrapper, tarball, filelist, and `merge_and_report.py` from `add-jpsi-x-condor-production` are unchanged and ready to use once iteration converges — at which point a *separate, sibling* change drives the scale-up run.

This keeps the iteration loop as fast as possible while the code is fluid, and avoids the ~7-iteration Condor smoke shakedown we already paid for once.

## What this change does NOT do

- **No cut-value changes.** The iter-1 values (`vtxProb > 0.01`, `cosα > 0.95`, `Lxy/σ > 1.0`) from `add-jpsi-x-vertex-fit-and-low-pt` stay. The point of this change is to make those cuts operate on correct underlying variables.
- **No bachelor PID or quality cuts.** See §4 above.
- **No preset A testing.** Preset A's baseline-shape role is no longer load-bearing in the analysis trajectory.
- **No Condor scale-up.** Sibling change.
- **No multi-era expansion.**

## Risks and mitigations

| Risk | Mitigation |
|---|---|
| J/ψ constraint code has a bug — affects all 7 channels (signal AND background dropped or shifted) | Local 1-file validation before declaring done. Compare candidate counts to pre-fix smoke; shifts > ±20 % flag a bug. |
| 2D cosα doesn't recover preset C — maybe the issue isn't only cosα | Cdiag data already proves Kalman convergence is healthy (88–96 %); only the three cuts can be the cause. cosα via beamspot is by far the most likely culprit. If 2D cosα is still wrong, switch to 3D-from-PV — one-line follow-up. |
| `KinematicConstrainedVertexFitter` is slower than `KalmanVertexFitter` | ~5–10 ms / candidate extra. For ~4 500 candidates per event, that's ~25 s/event added — non-trivial but acceptable for a local 10 000-event run. Production at scale gets a per-event cost re-measurement. |
| The fixed cuts now reject TOO LITTLE — preset C ends up looking like preset B | Indicates iter-1 cut values were correct in spirit but operating on a random axis; with the fix, even cosα > 0.95 should be meaningful. If preset C retains > 50 % of preset B post-fix, iter-2 (sibling change) tightens the threshold. Healthy outcome either way. |
| `RecoVertex/KinematicFit` BuildFile dependency breaks the build | Standard CMSSW package, no external deps. Test locally before claiming done. Low risk. |
| Local CMSSW area gets dirty during iteration | All edits are to one source file + one BuildFile. `git diff CMSSW_10_6_17_patch1/src/Alignment/CommonAlignmentProducer/plugins/JpsiXCandidateProducer.cc` is the audit trail. If iteration goes badly we can revert from git. |

## Validation strategy

Two nested gates:

1. **Local preset-B validation** (task 3.2-3): does the producer build and run with the J/ψ constraint? Does the B+ peak narrow (wide/tight window ratio → 1.4)?
2. **Local preset-C validation** (task 3.4-5): does ε_C jump from 0.6 % to ≥ 10 %? Does the B+ mass peak emerge?

If both pass, file the Condor scale-up sibling. If either fails, diagnose and loop within this change.

## Why this is a separate proposal from `add-jpsi-x-vertex-fit-and-low-pt`

That sibling change owns the *cut design* (which variables, which thresholds, why). This change owns the *implementation correctness* of those cuts (do the variables actually compute what their names suggest). The two have orthogonal failure modes: a perfect cut design will fail if the producer mis-computes the variables, and a perfect producer can still be configured with bad cuts.

After this change lands and validates, `add-jpsi-x-vertex-fit-and-low-pt` is empirically validated for the first time and can be archived. Iter-2 cut tuning, if any, gets its own change.
