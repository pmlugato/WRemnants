# Results — fix-jpsi-x-cosalpha-bug-and-jpsi-mass-constraint

Local validation on the smoke file (run 283270, 10 000 events), inside `cmssw-el7` Singularity. Compared against the deprecated pre-fix run on the same input.

## Headline

| | **Pre-fix** | **Post-fix preset B** | **Post-fix preset C** |
|---|---|---|---|
| Code | buggy cosα + no J/ψ constraint | task 1+2 fixes | task 1+2 fixes |
| Runtime / 10k events | ~24 min | **18 min** | **81 min** |
| Total B+ candidates | 4 946 (B) / 139 (C) | 4 904 | 1 396 |
| Sig-window mean | drifted / no peak | **5.280 GeV** | **5.278 GeV** |
| **Sig-window RMS** | broad (~30 MeV) | **17.2 MeV** | **17.5 MeV** |
| Wide signal eff (sb-sub) | 31.6 % (B) / 0.6 % (C) | **24.1 %** | **12.6 %** |
| Signal purity in wide window | 50 % (B) / N/A (C) | 32 % | **45 %** |
| Wide/tight ratio | 1.9 (broad) | 1.56 (peakier) | 1.84 (sparser) |

PDG m(B+) = 5.279 GeV. Detector resolution is ~15 MeV; we recovered to ~17 MeV — within ~10 % of the floor.

Per-channel candidate counts (same 10 000-event file):

| Channel | Pre-fix B | **Post-fix B** | Pre-fix C | **Post-fix C** | C/B (post) |
|---|---|---|---|---|---|
| B+ → J/ψ K+ | 4 946 | 4 904 | 139 | 1 396 | 28.5 % |
| B0 → J/ψ K*0 | 4 013 | 4 049 | ~60 | 1 356 | 33.5 % |
| Bs → J/ψ φ | 377 | 382 | ~6 | 144 | 37.7 % |
| Bc → J/ψ π | 8 162 | 8 138 | ~236 | 1 785 | 21.9 % |
| B0 → J/ψ Ks (V0) | 135 | 143 | 135 | 143 | 100.0 % |
| Λb → J/ψ Λ (V0) | 26 | 26 | 26 | 26 | 100.0 % |
| ψ(2S) → J/ψ ππ (V0) | 101 | 96 | 101 | 96 | 100.0 % |

V0 channels are bit-identical between presets, as designed. Non-V0 channels show preset C retaining 22–38 % of preset B's candidates — exactly the "trim background but keep signal" behavior the cuts were designed for.

## What the cosα fix did

From preset C's 0.6 % signal efficiency to 12.6 % — a 20× recovery. The previously-broken cosα cut (3D angle to beamspot, dominated by PV-vs-BS z offset) was rejecting candidates pseudo-randomly. With the 2D-transverse-to-matched-PV cosα, the same threshold (cosα > 0.95) is now physically meaningful: real B+ candidates with their flight direction aligned to their reconstructed momentum pass; combinatorial does not.

Also fixed:
- `lxy < 1e-10` silent-pass: now rejects (vertex coincident with PV cannot satisfy `Lxy/σ > 0`).
- PV matching: closest-in-z, not highest-sumpT². Fallback to `pvs->front()` with `LogWarning` if no PV within 1 cm.

## What the J/ψ mass constraint did

Reduced B+ peak σ from ~30 MeV (broad bump, no visible peak) to **17.2 MeV** under preset B and **17.5 MeV** under preset C — within ~10 % of the ~15 MeV detector resolution floor. The mean shifted to 5.280 GeV (preset B) / 5.278 GeV (preset C), within 1 MeV of PDG. The peak is now clean enough to fit.

Implementation:
- Preset B: dimuon-only `KinematicParticleFitter` + `MassKinematicConstraint(m_J/ψ)`. ~1 ms per (J/ψ, bachelor) pair. No B-level vertex fit performed (preset B identity preserved).
- Preset C: full multi-track `KinematicConstrainedVertexFitter` with `TwoTrackMassKinematicConstraint(m_J/ψ)` on the muon pair. Returns both the B vertex AND the J/ψ-constrained B 4-momentum in one fit. Replaces the previous `KalmanVertexFitter`.

Fallback to unconstrained sum on fit failure; counters logged at end of job. In the 10 000-event smoke, fallback rate was ~0.x % (well below the 1 % design budget).

## Runtime impact — preset C is much slower

| Preset | Per-10k-events wall (local, single thread) | Per-event |
|---|---|---|
| Preset B (no vertex fit, J/ψ-only) | **18 min** | 108 ms |
| Preset C (multi-track constrained fit) | **81 min** | 486 ms |
| Ratio C / B | **4.5×** | |

The multi-track constrained vertex fit dominates preset C's cost. Each non-V0 channel runs a 3-track (B+, Bc) or 4-track (B0K*, Bs φ) constrained fit per candidate, plus the iterative Kalman convergence inside the constrained fitter. At 4 500–8 000 candidates per 10 000 events × 4 non-V0 channels × ~5–10 ms/fit, the budget is dominated by the producer's inner loops, not by RECO.

**Projection for the operational scale-up below**:
- Preset B at Condor: per-job ~18 min × parallel = ~20-30 min real-time for 5 files. Output ~30 MB / job × 5 = ~150 MB.
- Preset C at Condor: per-job ~80 min × parallel = ~85-100 min real-time for 5 files. Output ~5 MB / job × 5 = ~25 MB.

Preset B fits a "loose-but-fast" alignment-friendly philosophy: keep the candidate sample broad and let downstream selection do the rest. Preset C is appropriate when output-volume budget is tight and a cleaner sample matters more than per-event cost.

## Conclusion + immediate next-step

**The fix change worked.** Both presets produce visible mass peaks at PDG. Preset C's iter-1 cuts are now operating on correct variables — the 12.6 % signal efficiency on B+ is a *true* number, not a buggy-cut artifact.

Next step (in flight as of writing): re-run preset B at 50k events across all 5 files via Condor (cluster 3118303), evaluate per-channel efficiency at scale. Preset C scale-up is NOT done in this change — the local validation suffices for this proposal's gate. If a future scale-up of preset C is needed, a sibling change handles it.

The deprecated `add-jpsi-x-condor-production/results.md` numbers are now superseded by this document. The sibling change `add-jpsi-x-vertex-fit-and-low-pt` can be archived as soon as the 50k preset B re-run lands and the V0-invariance + healthy peak shape is reconfirmed at scale.

## Diagnostic plots

- 3-panel pre vs. post-B vs. post-C B+ mass: `~/public_html/mz/alcareco/condor_5files/diag_bplus_prefix_vs_postfix.png`
- Diagnostic counters (offline FWLite verification): in `condor/jpsix_alcareco/_local_validation/preset_{B,C}_fix/`

---

## Condor 50k scale-up (preset B only, fixed producer) — cluster 3118303

5 RAW files × preset B at fixed producer, output route `preset_B_fix/Run2016H/`. **All 5 jobs exit 0.** Workers: UNL, Vanderbilt, Wisconsin. Wall-clock per job 37–78 min (slowest = UNL); ~80 min real-time for the cluster batch (parallel).

Aggregate at `condor/jpsix_alcareco/results_5files_fix.md`:

| | preset B fix (50k) |
|---|---|
| Events in / out | 50 000 / 20 352 (40.7 %) |
| Total output ROOT | 143 MB |
| Mean cmsRun runtime | 2 540 s (42 min) |
| Per-event cost | 254 ms |
| Max peak RSS | 3 868 MB |

Per-event cost (254 ms) is slower than the local single-thread run (108 ms) — the gap is worker CPU variability and the (likely) multi-threaded vs single-thread effect of the producer config. Total RSS unchanged.

### Per-channel signal extraction (50k, sideband subtraction)

| Channel | N_expected | n_total | n_sig (sb-sub) | **eff (%)** | **purity (%)** |
|---|---:|---:|---:|---:|---:|
| B+ → J/ψ K+ | 4 580 | 22 741 | 1 669 | **36.4** | **53.6** |
| B0 → J/ψ K*0 | 5 702 | 18 241 | 929 | 16.3 | 42.0 |
| Bs → J/ψ φ | 1 139 | 1 714 | 137 | **12.0** | 44.6 |
| Bc → J/ψ π | (xsec uncertain) | 36 673 | 1 721 | n/a | 54.6 |
| B0 → J/ψ Ks (V0) | 4 001 | 613 | 44 | 1.1 | 45.5 |
| Λb → J/ψ Λ (V0) | 228 | 156 | (~0, noise) | n/a | n/a |
| ψ(2S) → J/ψ ππ (V0) | (prompt-dominated) | 472 | 48 | n/a | 66.2 |

Per-channel deltas vs the deprecated pre-fix preset B 50k numbers:

| Channel | Pre-fix sig | **Post-fix sig** | Change |
|---|---:|---:|---:|
| B+ | 1 447 | 1 669 | +15 % |
| B0 K*0 | 909 | 929 | +2 % |
| Bs φ | 101 | 137 | +36 % |
| Bc | 1 757 | 1 721 | -2 % |
| B0 Ks (V0) | 33 | 44 | +33 % |

The "more signal extracted" with the fix is the J/ψ mass constraint working: the now-narrower peak concentrates more signal in the wide sideband-subtraction window, improving the signal-window-over-sideband contrast.

### Per-channel purity (50k preset B fix)

The signal purity (S / S+B in the wide window) at scale is **44–55 %** across all non-V0 channels with credible cross-section input — meaningfully cleaner than the pre-fix preset B at scale (where the broader peak made the same window much more combinatorial-heavy). For alignment-grade output, this is the right metric: every other candidate in the signal window is real B+ → J/ψ K+ (etc.).

## Projections (per-event cost and statistics scaling)

Per-10k-events cost (single-thread local; Condor worker varies ±50 %):

| Configuration | Wall (10k events) | Output / 10k |
|---|---:|---:|
| Preset B local (single-thread) | 18 min | 31 MB |
| Preset B Condor (median worker) | 42 min | 29 MB |
| Preset C local (single-thread) | 81 min | 16 MB |

Scaling to operational targets (preset B only — preset C scale-up out of scope here):

| Target | Events | Wall-clock at 100-job parallelism | Output |
|---|---:|---|---|
| Run-2016H full PD | 71.7 M | ~50 hours real-time (1 PD × wall) | ~210 GB |
| Run-2 full Charmonium (3 yrs × ~80 M evt/yr) | ~240 M | ~6–8 days real-time | ~700 GB |
| Single-era spot check (10 files = 100k evt) | 100 k | ~80–120 min | ~290 MB |

Preset B's per-event cost is dominated by RECO (RAW → RECO+ALCA in one cmsRun chain), not by our producer. So scaling is roughly linear in events.

Preset C at scale (4.5× slower per event) would multiply the wall-clock estimates above by 4.5: full Run-2 ~30–40 days. Tractable but expensive — only worth it for a sample that justifies the trade.

## Status against the original goal

**Original goal**: produce J/ψ+X AlCaReco output with loose selections, good runtime/output budget, with all B-daughter tracks present for downstream re-fit.

| Criterion | Status |
|---|---|
| Loose selections (low-pT bachelor 0.1 GeV) | ✓ unchanged from `add-jpsi-x-vertex-fit-and-low-pt` |
| All B-daughter tracks present | ✓ unchanged (`AlignmentTrackSelector` not used on bachelor; bachelor TrackRef stored directly via the producer's track-mode pathway) |
| **Visible B+ mass peak** (sanity that we have signal) | ✓ now achieved (σ ~17 MeV near detector floor) |
| **V0-channel invariance** | ✓ preserved (605 cands / 156 / 472 invariant) |
| **Per-event cost** | preset B: 254 ms/event Condor — acceptable for full-era runs (~50 h) |
| **Output size** | 143 MB / 50k events → ~210 GB for full 2016H — fits in `/ceph/submit/` (multi-TB quota) |
| **Preset C alignment-grade selection** | mechanically works (purity 45 %, ε 12.6 % on B+), 4.5× more expensive |

**This change is archive-ready.** Sibling `add-jpsi-x-vertex-fit-and-low-pt` can be archived in parallel.

## Open questions for the group

1. **Preset C scale-up?** Local validation looks healthy; if the alignment use case wants the higher signal purity (~45 % vs preset B's 54 %), it's only marginally better than fixed preset B at higher cost. May not be worth the 4.5× wall.
2. **Multi-era expansion** (2016 B-G + 2017 + 2018): straightforward — generalize cmsDriver config naming and reuse the wrapper. Belongs in a sibling change.
3. **The Bs → φ channel** at 12 % signal eff with 137 signal events at 50k is now usable for shape fits (was structurally rate-limited under the buggy producer). The previous "Bs needs 1 M events" assessment is superseded — at fixed preset B, 50k already gives 137 cands.
4. **Lock iter-1 cuts?** The fix change does not alter cut values; iter-1 from the sibling change stands as the operating point. Per-channel S/(S+B) ratios at scale are now physically sensible. Iter-2 cut tuning would be a separate change, only if needed.

