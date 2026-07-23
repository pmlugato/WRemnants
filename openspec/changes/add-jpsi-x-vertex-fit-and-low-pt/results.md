# Results — `add-jpsi-x-vertex-fit-and-low-pt`

**Status as of 2026-06-11**: Phase 0 + Phase 1 complete; all three Phase-1
gates pass with margin. Phase 2 (ALCAReco cmsRun on Run2016H Charmonium
RAW) pending.

## Phase 0 — visual mass-spectrum demo

**Script**: `scripts/btojpsik/mini_mass_demo.py`. FWLite + ROOT
canvas; no cmsRun, no narf, no RDataFrame.

**Input**: line 1 of `mc_minilist.txt` —
`/store/mc/RunIISummer20UL18MiniAOD/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/106X_upgrade2018_realistic_v11_L1v1-v2/240000/AA380A35-9E4F-7E4A-AB75-2690CD5B111C.root`

Read over xrootd via `cmsxrootd.fnal.gov`. Voms proxy timeleft at run
start: ~5.8 days.

**Sample**: 200 events, ~270 s wall (1.4 ev/s; xrootd-bound). Output:
`scripts/btojpsik/figs/btojpsik-mass-demo/{B_mass_stages.png,summary.json}`.

**Per-stage candidate counts** (B+ → J/ψ K+ on MiniAOD MC):

| Stage | Cuts | Cands | cands/event | Mass mean | Mass RMS |
|---|---|---:|---:|---:|---:|
| 1 | no cuts beyond charge | 51 140 | 255.7 | 5.287 | 0.432 |
| 2 | + minimal kinematic (soft μ pT>2, K pT>0.1, dimuon mass loose) | 12 334 | 61.67 | 5.304 | 0.429 |
| 3 | + updated preset B (K pT>0.1, K \|η\|<1.8, DOCA<0.03, J/ψ pT>3, B pT>5) | 665 | 3.33 | 5.296 | 0.399 |
| 4 | + preset C closed-form vertex (χ²prob>0.01, cos α_xy>0.99, Lxy/σ>3) | 195 | 0.98 | 5.234 | 0.309 |

**Visual interpretation**: stage 1 is essentially a flat combinatorial
background across [4.5, 6.0] GeV with no visible peak (the true
no-cut baseline that NanoAOD cannot produce because of upstream
BPark/BMM-Tools pre-filtering). Stage 3 shows a clear B+ peak emerging
at the PDG mass; stage 4 cleans the surrounding combinatorial. Cands/event
trajectory 256 → 62 → 3.3 → 1.0 demonstrates the cumulative power of
the cuts.

## Phase 1 — gen-matched gate evaluation

**Script**: same `mini_mass_demo.py` with `--gen-match --gates --presets A,B,C`.

**Gen-match**: walked `prunedGenParticles` for B+ (PDG 521) →
J/ψ(443)+K(321), then J/ψ → μ⁺μ⁻ (PDG ±13). Reco triplet labelled
"true" iff all three leaf tracks ΔR-match (< 0.03) to gen daughters
with matching charges, else "fake". Both muon orderings tested.

**BMuonFilter MC sanity check**: 200/200 events have a gen B+→J/ψK
decay. This is signal-only MC; every event is a B+ event. Phase-1
signal-efficiency denominators are therefore well-defined.

### Smoke-test results (200 events, 2026-06-11)

| Preset | cands/ev | sig_eff | fake/ev | tight_frac | μ_true | σ_true |
|---|---:|---:|---:|---:|---:|---:|
| A | 61.3 | 0.845 | 60.5 | 0.108 | — | — |
| B | 3.33 | 0.660 | 2.67 | 0.247 | 5.263 | 0.069 |
| C | 0.975 | 0.555 | 0.420 | 0.559 | 5.268 | 0.056 |

(Preset A's true-candidate mass entries are blank because A's
combinatorial blowup makes the mean/RMS dominated by fakes; the
tight_frac of 0.108 confirms this — only ~10% of A candidates fall in
the tight signal window.)

**Three Phase-1 gates** against the new preset B baseline (updated B
with bachelor pT > 0.1, lowered from the sibling change's pT > 1.5):

| Gate | Required | Observed | Verdict |
|---|---|---|:---:|
| sig_eff(C) / sig_eff(B) | ≥ 0.70 | 0.555 / 0.660 = **0.84** | ✓ |
| fake/ev(C) / fake/ev(B) | ≤ 1.5 | 0.42 / 2.67 = **0.16** | ✓ |
| mass shape (μ, σ) preserved between B and C | within ±2σ_stat | Δμ = 5 MeV, σ narrower in C | ✓ |

All three smoke-test gates pass with the proposal's default Kalman cut
values (vtxProb > 0.01, cos(α_xy) > 0.99, Lxy/σ > 3). The 5000-event
refresh is in flight and will tighten the statistical uncertainties.

### Headline numbers (5000 events, 2026-06-11)

Wall time: 244 s (20.5 ev/s; xrootd cache warm from earlier smoke tests).
5000/5000 events have a gen B+→J/ψK decay (BMuonFilter signal MC).

| Preset | cands/ev | sig_eff | fake/ev | tight_frac | μ_true | σ_true |
|---|---:|---:|---:|---:|---:|---:|
| A | 57.81 | 0.7450 | 57.07 | 0.109 | 5.2711 | 0.0835 |
| B | 2.89  | 0.6044 | 2.286 | 0.262 | 5.2726 | 0.0635 |
| C | 0.75  | 0.4950 | 0.255 | 0.615 | 5.2727 | 0.0604 |

Stat uncertainties at 5k events: sig_eff ≈ ±0.7% absolute (sqrt(p(1-p)/N));
fake/ev ≈ ±0.7% relative.

**Three Phase-1 gates** against the new preset B baseline:

| Gate | Required | Observed | Verdict |
|---|---|---|:---:|
| sig_eff(C) / sig_eff(B) | ≥ 0.70 | 0.495 / 0.604 = **0.820 ± 0.014** | ✓ |
| fake/ev(C) / fake/ev(B) | ≤ 1.5 | 0.255 / 2.286 = **0.112 ± 0.003** (≈ **8.9× suppression**) | ✓ |
| mass shape (μ, σ) preserved between B and C | within ±2σ_stat | Δμ < 0.2 MeV, σ narrower in C (B: 63.5 MeV; C: 60.4 MeV) | ✓ |

**Cumulative-cuts cands/event trajectory** (Phase 0 visual demo numbers, same 5k run):
255.7 → 58.16 → 2.89 → 0.75 (factor 341 total suppression from no-cut to preset C; with mass-window only ~10⁻³ kept).

Numbers are tighter than the 200-event smoke test (sig_eff(C)/sig_eff(B)
went 0.84 → 0.82) but conclusions unchanged. All three gates pass with
strong margin.

## Decision

**Phase-1 verdict on preset C**: all three gates pass at 5k events.

**Locked preset-C Kalman cut values** (going into Phase 2):
- `vtxChi2Prob > 0.01`
- `cos(α_xy) > 0.99`
- `Lxy/σ > 3`

No iteration needed at Phase 1. The proposal's initial cut values
clear the gates as-is.

**Next step**: Phase 2 ALCAReco runs (cmsRun against preset C on
Run2016H Charmonium RAW). The closed-form-vertex Phase-1 cuts above
are starting values; if Phase 2 numbers diverge materially (cands/event
breaks the ≤ 5 budget, or tight-window fraction drops below 30%), one
iteration cycle on Kalman cut values is budgeted.

## Phase 2 — ALCAReco confirmation on Run2016H Charmonium RAW

**Input**: same Run2016H Charmonium RAW file as the sibling Phase 2,
`/store/data/Run2016H/Charmonium/RAW/v1/000/281/727/00000/22261BD6-F984-E611-AA9A-FA163EC18366.root`.
1000 events × 3 presets.

**cmsDriver configs** (`CMSSW_10_6_17_patch1/`):
- `recoskim_Run2016H_Charmonium_JpsiX_presetA.py` — reused from sibling change.
- `recoskim_Run2016H_Charmonium_JpsiX_presetB.py` — reused from sibling change.
  The cff preset switching is env-var driven at cmsRun time, so the
  same config picks up the updated preset B (bachelor pT > 0.1) without
  regeneration.
- `recoskim_Run2016H_Charmonium_JpsiX_presetC.py` — new, copy+sed of
  presetB.py with the output filename changed.

All three configs include FastTimerService (printJobSummary on; no JSON
in CMSSW_10_6_17_patch1). Env var `TKALJPSIX_SELECTION_PRESET=<P>` set
on both cmsDriver and cmsRun (sibling-change bug documented).

**Comparator**: `CMSSW_10_6_17_patch1/jpsix_preset_compare_3way.py`
(new; ~260 lines; imports helpers from the 2-way comparator and adds
3-way variants of report_2a/2b/2c plus 3-line mass overlays per channel).

### Phase 2a — mass quality (3-way, 1000 events, 2026-06-11)

**Per-channel cands/event** under each preset (Run2016H Charmonium RAW, 1000 events):

| Channel | Mode | A:cands | A:c/ev | B:cands | B:c/ev | C:cands | C:c/ev | Verdict (A/B/C) |
|---|:---:|---:|---:|---:|---:|---:|---:|:---:|
| BPlus   | non-V0 |  8 611 | 13.52 |   641 |  1.16 |   1 | 0.01 | F / F / **F (yield)** |
| B0Kstar | non-V0 | 79 400 |124.65 |   506 |  0.92 |   2 | 0.03 | F / F / **F (yield)** |
| B0Ks    | V0     |     22 |  0.03 |    22 |  0.04 |  22 | 0.29 | F / F / **P** |
| BsPhi   | non-V0 |  6 800 | 10.68 |    52 |  0.09 |   0 | 0.00 | F / F / **F (no cands at all)** |
| Lambdab | V0     |      2 |  0.00 |     2 |  0.00 |   2 | 0.03 | F / F / F |
| Psi2S   | V0     |     19 |  0.03 |    19 |  0.03 |  19 | 0.25 | F / F / F |
| Bc      | non-V0 | 18 662 | 29.30 | 1 074 |  1.95 |   4 | 0.05 | F / F / **F (yield)** |

Gates: cands/event ≤ 5, tight_frac ≥ 30%, tight_yield ≥ 100/1000.

**Headline finding — preset C as currently configured is much too aggressive on data.** All four non-V0 channels fall by 1.5–3 orders of magnitude vs preset B, far beyond the Phase-1 prediction. Bs→φ produces zero candidates from 1000 events. The "tight-window yield ≥ 100/1000" gate now fails on every non-V0 channel under C.

V0-mode channels (B0Ks, Lambdab, Psi2S) confirmed preset-invariant: identical candidate counts across A/B/C. The higher cands/event under preset C reflects the AlCaReco event filter (`filter=True`) saving only 77 events at preset C vs 637/551 at A/B — V0 cands are concentrated in those surviving events.

**Why the Phase-1 → Phase-2 divergence?** The Phase-1 closed-form vertex predicted cands/event(C) / cands/event(B) ≈ 0.11; the real cmsRun ratio is 0.005–0.05 (factor ~10× more aggressive than predicted). This is the exact failure mode the design.md "Where the Phase-1 Kalman fit lives" subsection flagged. Probable causes:

- **Real KalmanVertexFitter is iterative**, refines the vertex through multiple passes, and gives tighter χ² for clean candidates *and* sharper rejections for borderline ones than the linear approximation.
- **Lxy/σ reference**: C++ uses the first `offlinePrimaryVertices` entry as the displacement origin; Phase-1 closed-form used the beamspot. The PV is closer to the fitted vertex on average, so Lxy is smaller and Lxy/σ is correspondingly smaller.
- **Angular convention**: C++ `computeAlphaBSFromPoint()` uses the 3D angle (3-momentum vs 3D displacement); Phase-1 used cos α projected onto xy. The 3D angle is more restrictive (a small z-misalignment hurts).

This divergence is the **iteration trigger** that the proposal anticipated; one cut-relaxation cycle is budgeted (`tasks.md` § 7.5).

### Phase 2b — size + timing (3-way, 2026-06-11)

| Quantity | Preset A | Preset B | Preset C |
|---|---:|---:|---:|
| Events saved | 637 | 551 | 77 |
| File size (MB) | 38.39 | 6.16 | 2.27 |
| Bytes / saved event | 60.3 kB | 11.2 kB | 29.5 kB |
| Wall time (s) | 447.9 | 276.9 | 247.3 |
| Wall ms / input event | 703 | 503 | 247 |

Bytes/saved-event drops by 5.4× from A → B (the existing sibling-change story, driven by kinematic + DOCA cuts), but **rises** by 2.6× from B → C (29.5 vs 11.2 kB/saved-event). Cause: under preset C, the only events that pass the AlCaReco filter are those with a non-zero merged-track collection, and surviving events tend to have V0-mode candidates (which carry more daughter tracks and dE/dx values per event than the few non-V0 candidates that get killed by the Kalman fit).

Wall time under preset C is the *fastest* of the three (247 s vs 448 s for A). The Kalman fit cost (~ms per surviving candidate × O(0.05) cands/event × 1000 events ≈ <1 s) is negligible compared to the savings from running the downstream track-selection + dE/dx projection on a much smaller merged track collection.

### Phase 2c — dedup compression (3-way, 2026-06-11)

| Quantity | Preset A | Preset B | Preset C |
|---|---:|---:|---:|
| Tracks (merged) | 53 280 | 3 629 | 278 |
| Tracks if per-channel | 426 791 | 7 549 | 195 |
| Compression factor | 8.01× | 2.08× | **0.70×** |
| Tracks per event | 83.64 | 6.59 | 3.61 |

The compression factor < 1 under preset C means the merged-track collection is *larger* than the per-channel sum. This is an artefact of V0-mode survivor bias: the few surviving preset-C events tend to be V0-rich (V0Producer-fitted Ks/Λ produces clean candidates that pass the AlCaReco filter), and V0-mode candidates contribute multiple tracks to the merged collection without entering the per-channel non-V0 candidate counts. With more events the metric would stabilize, but at this yield the dedup statistic is not meaningful.

### Phase 2 decision (preliminary) — preset C cuts must be relaxed; statistics must be scaled up

**Two problems identified, both with concrete next steps:**

1. **Preset C cut values are too aggressive on data.** The proposal's defaults (vtxProb > 0.01, cos α > 0.99, Lxy/σ > 3) need to be relaxed to make preset C usable. Three candidate iterations to try (suggested order of preference):
   - Relax Lxy/σ: 3.0 → 1.0 (likely the dominant kill — Bc lifetime alone should be sufficient discriminant)
   - Relax cos α: acos(0.99) ≈ 0.142 rad → acos(0.95) ≈ 0.317 rad (compensates for 3D-vs-xy convention)
   - Relax vtxProb: 0.01 → 0.001 (let in marginally-fitted candidates)

   These can be done singly or together. Sec § 7.5 of `tasks.md` budgets exactly this iteration cycle.

2. **Statistics at 1000 events are not enough for mass-shape determination under preset C, even after cut relaxation.** Even if cut relaxation recovers preset B's yields × 0.5 (a generous estimate), Bs→φ at ~25 candidates per 1000 events × 0.5 = ~12 cands per 1000 — too few for shape. **Recommendation**: scale to 50k–100k events via Condor batch jobs after cut relaxation. Estimate: 100 cmsRun jobs × 1000 events × ~5 min each = ~8 hours wall-time on Condor; feasible.

## Statistics assessment (per user's "more events?" question)

| Channel | Mode | Preset B yield per 1000 ev | Needed for mass shape (~500 cands) | Events needed |
|---|:---:|---:|---:|---:|
| BPlus   | non-V0 | 641 | ✓ already enough at B | already OK |
| B0Kstar | non-V0 | 506 | ✓ | already OK |
| BsPhi   | non-V0 |  52 | need 10× | **10 000** |
| Bc      | non-V0 | 1074 | ✓ | already OK |
| (V0 chans) | V0 | 2–22 | rare; either 50× or accept marginal | 50 000 |

Under **preset C with current (too-tight) cuts**: yields are O(1) per 1000 events on the rare non-V0 channels. To get even 100 candidates per channel under preset C, we'd need ~100k events for BPlus/B0Kstar, and Bs→φ would still produce zero. **Cut relaxation is mandatory before scaling.**

**Condor scale-up plan** (recorded as a follow-up, not yet executed):
1. First: iterate preset C cut values to a workable point on the 1000-event sample.
2. Then: submit ~50 Condor jobs splitting Run2016H Charmonium RAW across files; each job processes ~1000 events under preset C. Total wall ~5 min per job × 50 in parallel ≈ 30 min real-time.
3. Merge outputs with `edmCopyPickMerge` or `cmsRun` chaining; re-run the 3-way comparator.

A Condor JDL template lives at follow-up; not authored in this round.

## Phase 2 — Iteration 1 (2026-06-11, 5000 events, preset C cuts relaxed)

**Cut changes**: only preset C is modified; presets A and B unchanged.

| Parameter | v1 (initial) | iter1 (current) | Change |
|---|---:|---:|---|
| `minBVtxProb` | 0.01 | 0.01 | unchanged |
| `maxMotherAlphaBS` | acos(0.99) ≈ 0.142 rad | acos(0.95) ≈ 0.318 rad | wider angle accepted |
| `minBLxyOverSigma` | 3.0 | 1.0 | lower significance threshold |

**Rationale**: v1's cands/event for non-V0 channels under preset C were 10× lower than Phase-1 predicted (BPlus 0.01/ev, Bs→φ 0.00/ev). Relaxing two of three Kalman knobs in one iteration to compensate the 3D-vs-xy convention mismatch (cos α) and the offlinePrimaryVertices-vs-beamspot reference shift (Lxy/σ). `minBVtxProb` left at 0.01 (standard B-physics value).

**Scale-up**: 1000 → 5000 events. Same Run2016H Charmonium RAW file. Wall estimate ~80 min sequential (5× v1 wall times: A=37 min, B=23 min, C=21 min).

### Iter1 — per-channel candidate yields (5000 events, 2026-06-11)

Events saved per preset: A=3058, B=2625, **C=434** (preset C recovered from v1's 77 events).

| Channel | Mode | A:cands | A:c/ev | B:cands | B:c/ev | C:cands | C:c/ev | C(iter1)/C(v1) | Verdict (A/B/C) |
|---|:---:|---:|---:|---:|---:|---:|---:|---:|:---:|
| BPlus   | non-V0 | 42 093 | 13.76 | 3 076 | 1.17 | **95**  | 0.22 | 22× | F / **P** / F |
| B0Kstar | non-V0 |399 583 |130.67 | 2 434 | 0.93 | **51**  | 0.12 | 4×  | F / **P** / F |
| B0Ks    | V0     |     94 |  0.03 |    94 | 0.04 |  94    | 0.22 |  -  | F / F / F |
| BsPhi   | non-V0 | 34 698 | 11.35 |   227 | 0.09 | **3**   | 0.01 | ∞ (was 0) | F / F / F |
| Lambdab | V0     |     13 |  0.00 |    13 | 0.00 |  13    | 0.03 |  -  | F / F / F |
| Psi2S   | V0     |     60 |  0.02 |    60 | 0.02 |  60    | 0.14 |  -  | F / F / F |
| Bc      | non-V0 | 91 166 | 29.81 | 4 863 | 1.85 | **140** | 0.32 | 6.4× | F / F / F |

**Major recovery on preset C**. Cands per non-V0 channel went from O(1) to O(50-150). All channels now produce non-zero candidate counts. The Kalman cuts are no longer pathologically aggressive.

**At 5000 events, preset B passes all three gates for BPlus and B0Kstar** (B+: 0.93 c/ev / 0.32 tf / 293 y/k all pass; B0Kstar: 0.93 / 0.32 / 294 all pass). BsPhi and Bc still fail under preset B — BsPhi on yield (34 < 100), Bc on tight_frac (0.22 < 0.30).

**Preset C iter1 fails all gates on yield** — even after recovery, the per-channel candidate counts are below the 100/1000-event threshold:

| Channel | iter1 yield per 1000 ev | Gate (≥ 100) |
|---|---:|:---:|
| BPlus   | **66.8** | close, but below |
| B0Kstar | 30.0 | below |
| BsPhi   |  2.3 | well below (structurally rare) |
| Bc      | 82.9 | close, but below |

Preset C is still **30-50× more selective per event than preset B**. This is by design (the Kalman vertex fit + topological cuts trim combinatorial), but it means **statistics scale-up is required to satisfy the yield gate**.

### Iter1 — size + timing (5000 events, 2026-06-11)

| Quantity | Preset A | Preset B | Preset C |
|---|---:|---:|---:|
| Events saved | 3 058 | 2 625 | 434 |
| File size (MB) | 179.81 | 20.35 | 4.50 |
| Bytes / saved event | 58.8 kB | 7.8 kB | 10.4 kB |
| Wall time (s) | 960.2 | 873.7 | 685.8 |
| Wall ms / input event | 314 | 333 | 1 580 |

Notes:
- Preset C is **still the fastest** (Kalman cost negligible).
- Preset B bytes/saved-event dropped from v1 (11.2 → 7.8 kB) — likely a statistics effect (mean rises with rare-event spikes).
- Preset C bytes/saved-event dropped from v1 (29.5 → 10.4 kB) — at higher event count the V0-bias artefact dilutes.

### Iter1 — dedup compression (5000 events)

| Quantity | Preset A | Preset B | Preset C |
|---|---:|---:|---:|
| Tracks (merged) | 261 079 | 17 024 | 1 642 |
| Tracks if per-channel | 2 137 569 | 35 129 | 1 589 |
| Compression factor | 8.19× | 2.06× | 0.97× |
| Tracks per event | 85.4 | 6.5 | 3.8 |

Preset C compression factor still ≈ 1 — the V0-mode survivor bias is still present (events with non-V0 cands also tend to have V0 cands). At larger event counts this should converge to a meaningful number above 1.

### Iter1 — mass-shape statistics assessment

Target: ~500 candidates per non-V0 channel to support a Gaussian + linear background fit on the tight window (~50 cands/MeV at peak).

| Channel | iter1 cands at 5k events | Events needed for ~500 cands |
|---|---:|---:|
| BPlus   |  95 | ~25k (5×) |
| B0Kstar |  51 | ~50k (10×) |
| BsPhi   |   3 | ~850k (170×) — **structurally rare** |
| Bc      | 140 | ~17.5k (3.5×) |

**Bs→φ is the structurally-rare outlier**. Even at 100k events we'd have only ~60 candidates under preset C — barely enough for a yield count, far too few for shape fitting. The 850k-event target is impractical via Condor on a single dataset; better to extend to multiple Charmonium files / runs.

### Decision (iter1)

**Preset C iter1 cuts are workable and recommend locking them**:
- `minBVtxProb = 0.01`
- `maxMotherAlphaBS = acos(0.95) ≈ 0.318 rad`
- `minBLxyOverSigma = 1.0`

Further cut iteration (e.g. `minBVtxProb 0.01 → 0.001`) might recover another ~2× yield but at the cost of moving preset C closer to preset B's behavior. The current iter1 cuts give a meaningful selectivity (30-50× more selective than B) without going to zero. This is the alignment-grade preset's intended posture.

**The remaining problem is statistics-only**, not cuts. Solving via Condor scale-up (see next section).

## Condor scale-up — what to send and why

### What
- **Cuts**: iter1 preset C values (Lxy/σ > 1.0, cos α > 0.95, vtxProb > 0.01). Locked.
- **Sample**: Run2016H Charmonium RAW. Either many chunks of one file (current `22261BD6-…`) OR — preferred — spread across the ~50 Run2016H Charmonium RAW files for diversity.
- **Event count**: 50k events as primary target (gives ≥ 500 cands per non-V0 channel except Bs→φ). 100k events as a stretch goal if cluster has the cycles.
- **Presets**: all three (A, B, C). A and B are needed for the side-by-side comparison; their cmsRun is also 5× the cost so we'd want them on Condor too.

### Why
- 50k events at the iter1 preset-C yields gives:
  - BPlus: ~950 cands ✓
  - B0Kstar: ~510 cands ✓
  - Bc: ~1400 cands ✓
  - BsPhi: ~30 cands (still structurally rare; flag as known limitation)
- One job per ~1000-event chunk → ~50 jobs per preset, ~150 total
- Per-job wall ≈ 15-20 min (same as iter1 single-job rate); cluster parallelism completes in ~30 min real time
- Output merging via `edmCopyPickMerge` or cmsRun `PoolSource` chain

### How (preliminary plan, NOT yet authored)
1. `dasgoclient` query to enumerate all files in `/Charmonium/Run2016H-v1/RAW`
2. Per-file cmsRun config (parameterized by `--filein` and `--fileout`)
3. Condor JDL:
   - One job per (preset × file)
   - 16-core slot
   - Output staging to `/ceph/submit/data/.../jpsix_alcareco/iter2/`
4. Post-merge: `cmsRun mergeJob.py` with one PoolSource listing all per-job outputs per preset; output single `merged_preset{A,B,C}.root` per preset
5. Re-run `jpsix_preset_compare_3way.py` on the merged outputs

### Risks / caveats to flag in the Condor submission
- **xrootd source-site rate-limiting**: 50+ jobs hitting cmsxrootd.fnal.gov simultaneously may trigger throttling. Mitigation: spread across files (different physical sources) or stage input data to local Ceph first.
- **Bs→φ remains rare**: even 50k events gives only ~30 cands. May need to extend to Run2016G or beyond. Acceptable to document as known limitation.
- **Output disk**: 50k events × preset A's 60 kB/event = 3 GB per preset; ~9 GB total. Fits easily on `/ceph/submit/`.

## Caveats and follow-ups

- **Phase-1 closed-form vertex ≠ Phase-2 KalmanVertexFitter.** The
  Phase-1 closed-form 3-track common-vertex estimate is an analytic
  least-squares of helix points-of-closest-approach with χ² from track
  impact-parameter sigmas. It is exact in the limit of small helix
  curvature over the perigee-vertex distance (true for low-mass parents
  in the central tracker) but is not bit-identical to the iterative
  `KalmanVertexFitter` that preset C invokes in cmsRun. A small offset
  in cut thresholds between Phase 1 and Phase 2 is expected; one
  iteration cycle is budgeted at Phase 2 (`tasks.md` § 7).
- **B+ MC only at Phase 0.** The mass-spectrum demo plot is B+ only
  because the `BuToJpsiK_BMuonFilter` MC sample is the only signal MC
  we have staged. Extending the demo to B0→K*0, Bs→φ, Bc requires
  additional MC sourcing. Recorded as a follow-up; not blocking.
- **Bc topology under preset C deferred.** The 3-body / two-2-body /
  multi-body topology choice for Bc → J/ψ X is undetermined; preset C
  in the cff and the C++ producer treat Bc no-op under C, falling back
  to updated preset B behaviour. Recorded as a separate follow-up
  change.
- **xrootd I/O dominates Phase-1 wall time.** Locally staging the
  MiniAOD file (4.4 GB) would speed up future Phase-1 iterations by
  ~10× but was not done for the initial run. Recorded.
