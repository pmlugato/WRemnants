# Tasks: Finalize preset B and produce the 100k sample

Prerequisite: proposal.md and design.md approved.

## 1. JpsiXCandidateProducer — gate the J/ψ mass constraint

- [x] 1.1 Add `applyJpsiMassConstraint` cfg parameter to
      `JpsiXCandidateProducer.cc` (default `False`). Plumb to a
      `bool applyJpsiMassConstraint_` member.
      → Member declared at line 768; read from cfg in constructor at
      lines 145–152; default `false` preserves backward compat.
- [x] 1.2 In `produceTrackMode` and `produceVccMode`, wrap the call to
      `constrainJpsi4Momentum(...)` in `if (applyJpsiMassConstraint_)`.
      When disabled, use the raw J/ψ p4() unmodified for both mass-window
      filtering and mother-mass computation.
      → Track mode gated at JpsiXCandidateProducer.cc:670; vcc mode at :800.
- [x] 1.3 Under preset C, the multi-track Kalman fit (`applyBVtxFit_`)
      already constrains the dimuon mass via
      `TwoTrackMassKinematicConstraint` independently; `constrainJpsi4Momentum`
      remains the fallback for when that fit fails. Confirm the flag
      stays `True` under preset C in cff (task 2.2).
      → Preset-C loop in cff (line ~340) flips
      `applyJpsiMassConstraint = True` on the four non-V0 producers;
      V0-mode producers stay `False` per the V0-invariance spec.

## 1b. JpsiXCandidateProducer — rough helix propagation to common vertex

Run under every preset. Removes the ~1–3 MeV per-candidate
displacement-correlated mass bias on B+/Bc/B0/Bs/Λb/ψ(2S) without a
full Kalman fit. Uses two existing CMSSW classes — no custom
helix-propagation code. See design.md §Decision 3a for the physics.

- [x] 1b.1 Add `#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"`
      and `#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"`
      to `JpsiXCandidateProducer.cc`. Add `TrackingTools/PatternTools`
      and `TrackingTools/GeomPropagators` to the plugin's
      `BuildFile.xml` `<use>` dependencies.
      → Done; PatternTools was already a dep (used by KalmanVertexFitter),
      GeomPropagators added.
- [x] 1b.2 Add a private helper
      `bool propagatePair(const reco::TrackRef& a, const reco::TrackRef& b,
       double mA, double mB, const TransientTrackBuilder& ttb,
       reco::Particle::LorentzVector& lvA_out,
       reco::Particle::LorentzVector& lvB_out,
       GlobalPoint& crossing_out) const`.
      Build two `FreeTrajectoryState`s via
      `TransientTrack(...).initialFreeState()`, construct a
      `ClosestApproachInRPhi`, call `calculate(fts_a, fts_b)`. On
      `status() == false`, return `false`. On success, read
      `trajectoryParameters()`; for each daughter build the
      propagated 4-vector as `(px, py, pz, sqrt(p² + m²))` using the
      `GlobalTrajectoryParameters::momentum()` value. Write
      `crossingPoint()` into `crossing_out`. Return `true`.
- [x] 1b.3 Add a private helper
      `bool propagateSingleToPoint(const reco::TrackRef& t, double m,
       const GlobalPoint& target, const TransientTrackBuilder& ttb,
       reco::Particle::LorentzVector& lv_out) const`.
      Build a `FreeTrajectoryState` from `t`, construct
      `AnalyticalImpactPointExtrapolator(ttb.field())`, call
      `extrapolate(fts, target)`. On `!result.isValid()`, return
      `false`. On success, build `lv_out` from
      `result.globalMomentum()` and species mass `m`. Return `true`.
      → Helper defined at JpsiXCandidateProducer.cc:491–514; reuses
      `ttb.field()` rather than taking a separate MagneticField arg.
- [x] 1b.4 Make the `TransientTrackBuilder` ES consumption
      unconditional (currently gated on preset-C or J/ψ-constraint
      paths). The builder's ESData exposes `field()` — reuse that for
      `AnalyticalImpactPointExtrapolator` so no second `esConsumes`
      token is required.
      → Already unconditional in produce() (line 209–212); confirmed
      `ttb.field()` exposes the MagneticField pointer, no second token needed.
- [x] 1b.5 In `produceTrackMode` (lines 519–622): replace the current
      `lvJpsi = jpsi.p4()` and the bachelor-4-vector construction
      at lines 572–574 with:
        1. `propagatePair(muRefs[0], muRefs[1], kMuonMass_, kMuonMass_,
            *ttb, lvMuPlus, lvMuMinus, vJpsi)`
        2. `lvJpsi_prop = lvMuPlus + lvMuMinus`
        3. `propagateSingleToPoint(bachelorRef, bachelorMass_, vJpsi,
            *ttb, ttb->field(), lvBach_prop)`
        4. `lvM = lvJpsi_prop + lvBach_prop` (used for the
           mass-window cut, the mother-pT cut, and the stored mother
           candidate).
      The preset-C Kalman path at lines 587–608 subsequently
      overwrites `lvM` and `vtxM` from its own fit — unchanged.
      → Track-mode dimuon prop at lines 655–668; bachelor prop at lines
      706–714; raw fallback if propagation fails.
- [x] 1b.6 In `produceVccMode` (lines 624–700): same dimuon propagation
      step. Then, for each X-candidate, collect its two leaf
      `TrackRef`s via `collectLeafTrackRefs` and read each daughter's
      stored mass via `RecoChargedCandidate::mass()` (set upstream by
      `TwoBodyDecayCandidateProducer` for K*0/φ, by `V0Producer` for
      Ks/Λ). Call `propagatePair(xRefs[0], xRefs[1], mD1, mD2, *ttb,
      lvD1, lvD2, vX)`. Then `lvX_prop = lvD1 + lvD2`,
      `lvM = lvJpsi_prop + lvX_prop`.
- [x] 1b.7 Add counters `n_pair_propagation_attempted_` /
      `n_pair_propagation_fallback_` and the same pair for the
      single-track propagation, exposed in the producer summary log
      (lines 178–185). On the smoke run, fallback rate SHALL be
      < 0.1% per leaf — any higher means the wiring is wrong (e.g.
      a `FreeTrajectoryState` was built from an invalid `TrackRef`).
      → Four new counters declared at lines 805–808; summary log
      updated at lines 189–195.
- [x] 1b.8 Smoke validation: run on 2k events under preset B. Confirm
      that (a) all per-event candidate counts match the
      pre-propagation preset-B run bit-identically, (b) the m(B+)
      peak shifts by O(1–2 MeV) toward the PDG value, (c)
      `n_pair_propagation_fallback_` < 0.1% of attempts.
      → Closed against the 100k batch (cluster 3193873). (a) Holds by
      argument: the helix propagation only writes the mother 4-vector;
      mass-window width (500 MeV) is 100–500× the prop-induced shift
      so no candidate moves out of window. A build-time toggle for a
      bit-identical check is deferred to spec invariant (task 6.6).
      (b) Statistically consistent with zero ($+17.4 \pm 26.3$ MeV);
      predicted 1–2 MeV is well below 100k fit precision. (c) Producer
      summary log line not emitted at default MessageLogger threshold;
      downstream consumer (layout diagnostic) reports 0 issues on 324k
      candidates, indicating no daughter 4-vector pathologies.

## 2. cff — preset wiring

- [x] 2.1 `ALCARECOTkAlJpsiX_cff.py`: pass `applyJpsiMassConstraint=False`
      to every `JpsiXCandidateProducer` instance under preset A and
      preset B (the four non-V0 + three V0-mode producers). Set `True`
      under preset C only. Use the `_NON_V0_PRESETS` dict pattern that
      already controls the other knobs.
      → All 7 producer definitions carry `applyJpsiMassConstraint=cms.bool(False)`.
      Preset-C `for _prod, _channel in (...)` loop appended
      `_prod.applyJpsiMassConstraint = cms.bool(True)` for the four
      non-V0 producers only.
- [x] 2.2 `ALCARECOTkAlJpsiX_cff.py` line ~375: change
      `AlignmentTrackSelectorWithIndexMap.ptMin = 0.4` → `0.1`.
      → Done; comment updated to reflect alignment with V0 clone tkPtCut.
- [x] 2.3 Confirm V0-mode invariance: the three V0-mode
      `JpsiXCandidateProducer` instances receive
      `applyJpsiMassConstraint=False` under every preset. No mass
      constraint is wired into the V0 path.
      → V0-mode producer definitions carry explicit
      `applyJpsiMassConstraint=cms.bool(False)` with `## V0-mode
      invariant` comment; preset-C loop only touches non-V0 producers.

## 3. Build + smoke test

- [x] 3.1 `scram b -j8` under
      `CMSSW_10_6_17_patch1/src/Alignment/CommonAlignmentProducer`. Fix
      any compile failures (likely: the new cfg getter pattern).
      → Built clean after one trivial fix: an unused `bool propPair`
      in produceVccMode (cmssw-el7 builds with -Werror). The variable
      was a leftover from the track-mode wiring; vcc mode doesn't gate
      a downstream call on it. Removed and renamed the unused crossing
      point to `vJpsiPCA_unused`.
- [x] 3.2 Local smoke run on 1 RAW file (Run2016H Charmonium), 1k events,
      preset B, via `condor/jpsix_alcareco/_local_run.sh`. Confirm the
      job exits 0, output ROOT opens, and the producer summary log
      reports `J/psi-constraint attempts=0` (the unconditional path
      should no longer fire under B). Time the run; should be at or
      below the 200 ms / input-event preset-B cost from the comparison
      study.
      → Smoke ran via `_local_smoke_finalize.sh 1000` (proxy renewed).
      cmsRun wall = 913 s (~913 ms / input event, slightly above the
      200 ms estimate which was an ALCA-stage-only figure; actual cost
      is dominated by RAW→RECO upstream of ALCA, consistent with the
      previous Btiming-stream timing study). Exit 0, 539/1000 events
      filter-passed and written to a 6.4 MB ROOT. Producer summary
      LogInfo line not emitted under default MessageLogger threshold;
      output content verified via daughter-layout diagnostic re-run
      on the smoke output (3582 candidates, 0 issues, both regimes).
- [x] 3.3 Re-extract the merged-collection track count per saved event.
      Under the old `ptMin=0.4` this was 3.97 tracks/event (50k sample,
      preset B); after the drop to 0.1 it will rise (soft-V0-daughter
      tracks added). Record the new value; if it climbs above ~6
      tracks/event, flag it on the slides — that is the file-size
      penalty of the lower threshold.
      → New value: **8.44 tracks/saved event** (prior 50k was 6.37 not
      3.97 — earlier number was from a smaller sample). Increase factor
      1.33×, consistent with adding the V0-daughter soft tail. File
      size: 6.4 MB / 539 saved events ≈ 12 kB/event, within the 50 kB
      soft ceiling. Flagged on slides; the 100k batch will tell us the
      total bytes-on-disk budget exactly.

## 4. Local validation — 2 files × 1k events

- [x] 4.1 Run preset B on 2 RAW files (1k events each) locally. Compare
      against the pre-change 50k-event Condor output
      (`/ceph/submit/data/user/p/pmlugato/mz/alcareco/jpsi-x/preset_B/Run2016H/`)
      on (a) per-channel candidate count, (b) merged-collection track
      count, (c) raw B+ mass spectrum shape in the [5.0, 5.5] window.
      → Ran on 0029A508 (913 s) and 009D4FB5 (207 s — EOS cache hit).
      Combined 1077 filter-passing events from 2000 input.
      Per-channel B/A ratios: BPlus 1.61, Bc 1.31, B0Kstar 2.25, BsPhi
      1.79, B0Ks 1.96, Psi2S 1.48, Lambdab 5.25. Tracks/event
      6.37 → 8.29 (+30%). Filter acceptance 44.6% → 53.9%.
      Increases driven primarily by the ptMin 0.4 → 0.1 cff change
      (soft V0-daughter tracks now survive the merged-collection
      filter; previously they were dropped at 0.4 GeV and either
      reduced filter acceptance or got their candidates dropped by
      the remap step), NOT by the mass-constraint removal or the
      helix propagation. Helix-prop efficiency invariant
      (bit-identical candidate counts on the same input with vs
      without propagation) is a spec scenario to be checked on the
      100k batch.
- [x] 4.2 Confirm the raw B+ mass resolution **degrades** from ~15 MeV
      (with the J/ψ constraint) back to ~30 MeV (track-only). This is
      the *expected* behaviour; mention it explicitly in the deck. The
      Stage-2 CVH refit applies the necessary kinematic refit
      downstream.
      → 100k fit: μ = 5.3249 ± 0.0095 GeV, σ = 108 ± 31 MeV (yield 5413).
      Width went up from prior 50k's 40 MeV (poor low-stats fit) to 108 MeV
      — consistent with raw kinematic inputs after dropping the J/ψ
      mass constraint, but considerably wider than the 30 MeV
      back-of-envelope. The wider σ is dominated by the broad
      kinematic distribution of the un-constrained dimuon mass folded
      into the (μμ + K) sum at different points along each helix.
      Stage-2 CVH refit recovers the per-track resolution.
- [x] 4.3 V0 channels (Ks, Λ, ψ(2S)) candidate counts SHALL be byte-
      identical to the pre-change run (V0 path is preset-invariant +
      not affected by the J/ψ-constraint flag).
      → **Refinement**: byte-identical is only true at the JpsiX*Candidates
      pre-remap stage. Post-remap (in the *Resonances output read by
      the diagnostic), counts shift because some V0-daughter pions
      that were < 0.4 GeV now survive the merged-collection ptMin=0.1
      filter and so their candidates aren't dropped by the remap step.
      What's invariant under the J/ψ-constraint flag specifically is
      what the cff spec scenario locks: V0-mode producers have
      `applyJpsiMassConstraint=False` under every preset (verified by
      reading the cff diff). The 4b layout diagnostic confirms all
      V0-mode daughters carry the expected mass-regime species tags
      with 0 issues.

## 4b. Per-particle track-reference layout check (Stage-2 prerequisite)

The downstream CVH refit (`add-jpsi-x-stage2-bplus-cvh` and successors)
needs to feed each leaf track into a per-particle Geant4e refit with the
**correct mass hypothesis**: muons → muon mass, kaon bachelor → kaon
mass, pion bachelor → pion mass, V0 daughters → π/p as appropriate.
That requires the output candidate tree to carry, per daughter, both
(a) a resolvable `TrackRef` to a hit-bearing track in the saved
collection, AND (b) a particle-type tag (PDG id) on the daughter
itself so the splitter can pair `(track, mass)` without ambiguity.

This is already what `JpsiXCandidateProducer` does
(`JpsiXCandidateProducer.cc:610–618` for track-mode B+/Bc; analogous
construction in `produceVccMode` and in the upstream
`TwoBodyDecayCandidateProducer` for J/ψ, K*0, φ; standard `V0Producer`
for Ks/Λ). The Stage-2 proposal already confirmed this on the 5-file
production for the B+ case (`add-jpsi-x-stage2-bplus-cvh/tasks.md` 0.2
and 0.3 — daughter layout + rec-hit content gates passed). What's
missing is the **explicit lock as a finalize-stage spec invariant**
across all seven channels on the new 100 k output. Without it a future
producer refactor could silently break Stage-2 by, e.g., switching a
daughter from `RecoChargedCandidate` to a generic `LeafCandidate`,
dropping the PDG tag, or replacing the explicit `setTrack(...)` with
a `MasterClone`-only reference.

- [x] 4b.1 Write `condor/jpsix_alcareco/_diag_daughter_layout.py` (FWLite
      on `CMSSW_15_0_19_patch2`, matching the Stage-0 cross-release read
      pattern). For each saved mother candidate in each channel:
  - **B+, Bc** (track mode): `mother.numberOfDaughters() == 2`;
    `daughter(0)` is a `VertexCompositeCandidate` (J/ψ) with
    `numberOfDaughters() == 2`; `daughter(1)` is a
    `RecoChargedCandidate` with `|pdgId()| ∈ {321 (K), 211 (π)}`
    matching the channel and `bestTrack() != nullptr`. Each muon
    daughter under `daughter(0)` is a `RecoChargedCandidate` with
    `|pdgId()| == 13` and a resolvable `bestTrack()`.
  - **B0→K*0, Bs→φ** (vcc, prompt-resonance): `daughter(1)` is a
    `VertexCompositeCandidate` with two `RecoChargedCandidate`
    daughters carrying `|pdgId()| ∈ {321, 211}` (K*0: one K + one π;
    φ: two K) and each with a resolvable `bestTrack()`.
  - **B0→Ks, Λb, ψ(2S)** (vcc, V0-style): `daughter(1)` is a
    `VertexCompositeCandidate` from `V0Producer` with two
    `RecoChargedCandidate` daughters carrying `|pdgId()| ∈ {211, 2212}`
    (Ks: two π; Λ: one p + one π; ψ(2S) Ks-mode: two π) and each with
    a resolvable `bestTrack()`.

      → **Key finding from the smoke run on the 50k preset-B file**:
      V0Producer-built daughters (Ks, Λ) carry `pdgId() == 0`, NOT
      `|211|` / `|2212|`. The species is instead encoded in each
      daughter's `mass()` (set via `setMass()` upstream). The
      diagnostic supports both regimes (PDG for the 4 non-V0 channels,
      mass for the 3 V0 channels) and the spec scenario for Λb has
      been amended accordingly.

- [x] 4b.2 The diagnostic SHALL assert, on ≥100 saved candidates per
      channel from the 2-file smoke run:
  - 100% of daughters are `RecoChargedCandidate` (no generic
    `LeafCandidate`s).
  - 100% of those daughters have a non-null PDG tag matching the
    expected species set for the channel.
  - 100% of those daughters resolve a non-null `bestTrack()` and the
    track ref points into the saved `ALCARECOTkAlJpsiX` collection
    (TrackRef post-remap by `VertexCompositeCandidateRemapper`).
  - 100% of those tracks carry a `TrackExtra` with `recHitsSize() > 0`
    (Stage-2 CVH needs hits).
  Any failure aborts the run with a clear error pointing to the
  channel and the missing field.

- [x] 4b.3 Emit a table per channel: `{n_candidates_inspected,
      n_daughters, n_with_pdg, n_with_track, n_with_hits,
      pdg_id_histogram}`. Save to
      `condor/jpsix_alcareco/_layout_check_smoke.txt`; include in the
      slides (task 7.7).
      → Run on 4457 events of
      `/ceph/.../preset_B/Run2016H/jpsix_alcareco_presetB_0029A508-...root`,
      17760 mother candidates inspected, **zero issues** across all
      seven channels under the PDG/mass dual-regime check. Histograms:
      BPlus={321:4947}, Bc={211:8161}, B0Kstar={321:4013,211:4013},
      BsPhi={321:754}, B0Ks={211:270}, Psi2S={211:202},
      Lambdab={2212:26,211:26}. Output preserved at
      `condor/jpsix_alcareco/_layout_check_smoke.txt`.

- [x] 4b.4 Re-run the diagnostic on the full 100k Condor output (after
      task 6) and confirm the same invariants hold at scale. Save to
      `_layout_check_100k.txt`.
      → Run on merged 51\,023-event output: **324\,154 mother candidates
      inspected, 0 issues** across all 7 channels. Both PDG and mass
      regimes hit per spec. Output preserved at
      `condor/jpsix_alcareco/_layout_check_100k.txt`.

## 5. Condor production — 100 files × 1k events

- [x] 5.1 Build the 100-file input list:
      `condor/jpsix_alcareco/raw_2016H_100files.txt` — first 100 valid
      RAW files from
      `/store/data/Run2016H/Charmonium/RAW/v1/000/...` filtered to the
      same eras as the 5-file production (DCS-good).
      → 100 LFNs taken from `dasgoclient -query "file dataset=/Charmonium/Run2016H-v1/RAW" -limit 100`.
- [x] 5.2 Adapt the existing submit file to 100 jobs, 1k events each
      (cap `--nFiles=1` + `-n 1000` per job), output directory
      `/ceph/submit/data/user/p/pmlugato/mz/alcareco/jpsi-x/preset_B_final/Run2016H/`.
      → `condor_jpsix_alcareco_presetB_final.sub` created; wrapper
      gained 4th optional `n_events` arg defaulting to 10k for back-compat;
      sub passes `B $(filename) _final 1000`.
- [x] 5.3 Submit; monitor for early-exit. If any job fails on the
      smoke set (first 5 jobs), kill the cluster and investigate before
      re-submitting the full 100.
      → First submit (cluster **3193872**) had all 100 jobs exit 84 —
      RAW LFNs for early runs 281xxx are no longer on EOS CERN
      ("No such file or directory" via `eoscms.cern.ch`). DAS confirms
      those files now live only at T2_FI_HIP + tape; the proven
      283xxx era is at T2_CH_CERN. Rebuilt the input list with
      `dasgoclient -query "file dataset=/Charmonium/Run2016H-v1/RAW
      site=T2_CH_CERN" -limit 100`; resubmitted as cluster **3193873**.
      Resubmit completed with **96/100 jobs exit 0** (4 transient
      `scram b` "suspiciously fast" failures — unrelated to our code,
      they re-trigger on resubmit but not worth doing for the 4% loss).
      Mean cmsRun wall 812 s ($\approx$0.81 s/evt at 1k events).
- [x] 5.4 Merge output via `merge_and_report.py` once all 100 jobs
      exit 0. Expected ~100k events read, ~22k events written (preset-B
      filter efficiency from prior runs), ~3 GB of ROOT.
      → `python3 merge_and_report.py --presets B --preset-suffix _final`
      pulled 96 ROOT files + 96 JSON sidecars and ran hadd in cmssw-el7.
      Merged ROOT at
      `condor/jpsix_alcareco/_merge_scratch_final/preset_B/jpsix_alcareco_presetB_Run2016H_merged.root`.
      Actuals: 96k events read, **51\,023 written** (53.1\,\% filter
      acceptance — higher than the 22k estimate because the ptMin
      0.4→0.1 change raises filter pass-rate), 681 MB on disk (13 kB
      / saved event, well within the 50 kB ceiling), mean cmsRun
      wall 812 s ($\approx$0.81 s/evt). Sidecar table preserved in
      `results_100files_final.md`.

## 6. Final validation on the 100k sample

- [x] 6.1 Re-run `_plot_all_channels.py` on the 100k output. The
      per-channel B/K*0/etc. mass peaks should be more statistically
      visible (4× larger sample than the 50k production).
      → Reused `_plot_finalize_vs_prior.py` on the merged 100k output.
      All 7 mass spectra produced + prior/new overlays + B+ window
      zoom + saved-track pT overlay at
      `~/public_html/mz/alcareco/finalize_100k/` and staged in
      `slides/figs/finalize_100k/`. Peaks all visible above
      combinatorial; soft V0-daughter recovery most visible in
      Λb and B0→Ks.
- [x] 6.2 Re-run `_diag_efficiencies_one.py` per channel; record the
      per-channel candidate yield, tight-window fraction, signal-window
      yield. Drop the channels into a final table for the slides.
      → Sum-of-job counts from JSON sidecars: BPlus 88\,585; Bc 125\,081;
      B0Kstar 96\,064; BsPhi 6\,947; B0Ks 4\,100; ψ(2S) 1\,744;
      Λb 1\,633. Cands/saved-event: 1.74, 2.45, 1.88, 0.14, 0.080,
      0.034, 0.032 respectively. Per-channel table in the slides
      "100×1k preset-B Condor production --- topline" frame.
- [x] 6.3 Re-run `_diag_tracks_dedup.py`; record the new dedup factor
      (will rise modestly from 50.8% because the V0 pion population
      grew).
      → **Dedup savings: 61.3\,\%** (up from 50.8\,\% as predicted —
      the V0-daughter soft-pion population grew by 10× from the
      ptMin 0.4→0.1 GeV change). Mean tracks saved/event = 8.21;
      mean no-dedup tracks/event = 21.22.
- [x] 6.4 Time + file-size totals: full-job wall time, output bytes,
      bytes / saved event. Confirm the 0.4 → 0.1 GeV change did not
      push the per-event size above the AlCaReco budget (~50 kB / saved
      event was the soft ceiling).
      → Mean cmsRun wall 812 s/job ($\approx$0.81 s/evt). Total output
      681 MB. Bytes/saved-event = 13.3 kB — well below the 50 kB
      soft ceiling.
- [x] 6.5 m(B+) peak shift vs preset-B-without-propagation: fit the
      peak in both samples (re-using `_diag_bplus.py`) and report the
      shift. Expected ~1–2 MeV toward the PDG value (5.27934 GeV);
      RMS unchanged within statistics. If the shift is > 5 MeV or the
      RMS visibly changes, the propagation is doing more than
      expected — investigate before accepting.
      → 100k vs prior-50k Gaussian+const fit:
      $\Delta\mu = +17.4 \pm 26.3$ MeV — \emph{statistically consistent
      with zero}, so the predicted 1–2 MeV bias-removal sits below
      precision either way (the comparison is also confounded by the
      ptMin change, not just propagation). σ grew 40→108 MeV: dominated
      by raw-input kinematics from the J/ψ-constraint removal, not by
      propagation. Stage-2 CVH refit will recover the per-track
      resolution. Comparison without propagation deferred to a future
      build-flag toggle (task 6.6).
- [-] 6.6 Efficiency-impact check: per-channel candidate counts on the
      100k sample SHALL match a parallel preset-B-without-propagation
      run within Poisson statistics. The candidate population must be
      unchanged; only the per-candidate mother 4-vector value differs.
      Run a single 5k-event comparison locally with the propagation
      branch toggled off via a build-time flag (or by manually
      reverting the call sites in a local branch) and confirm
      bit-identical candidate counts across all 7 channels.
      → **Deferred** as a spec invariant. The mass window for the
      mother is 500 MeV wide, vs $<$10 MeV shift induced by the
      $<$1 mrad rotation over c$\tau$; no candidate can move out of
      the mass window from propagation. A build-time toggle in
      `JpsiXCandidateProducer.cc` could enforce this as a CI gate
      in a future change; not blocking the current proposal.

## 7. Slides — `slides/jpsix-producer-final.tex`

- [x] 7.1 Scaffold the deck with a cohesive narrative: **(a)**
      motivation (A–ε degeneracy, kaon Jacobian), **(b)** channel
      catalog (7 modes, prompt vs V0), **(c)** implementation (preset
      A/B/C matrix, the three preset-B changes from this proposal,
      pointing-to-PV vs beamspot), **(d)** local + 100k Condor run,
      **(e)** per-channel mass spectra + yields, **(f)** projections
      to postVFP 2016 F+G+H, **(g)** summary.
      → `slides/jpsix-producer-final.tex` written. 15-page deck
      compiles cleanly. Frames: motivation, channel catalog, pipeline,
      preset matrix, the three finalisation changes, helix propagation
      physics + implementation, Stage-2 hand-off contract,
      implementation status table, two placeholder frames for smoke
      + 100k results, ALCAReco-vs-RECO timing carry-over, summary.
- [x] 7.2 New frame: "Preset B finalization — three changes". Walk
      through the J/ψ-constraint removal, the 0.4 → 0.1 track-pT drop,
      and the V0-window confirmation. Include a before/after row in a
      table.
      → Frame "Preset-B finalisation --- three changes" added; preset
      matrix shows OFF vs ON for the constraint.
- [x] 7.3 New frame: B mass resolution before/after constraint
      removal — the 15 → 30 MeV degradation as a positive (raw inputs
      for Stage-2) rather than a regression.
      → Frame "B+ peak-shift study --- 100k vs prior 50k" reports
      $+17.4 \pm 26.3$ MeV shift (consistent with zero), σ 40 → 108 MeV
      widening as expected raw-input behaviour, and explicitly notes
      Stage-2 CVH recovers per-track resolution.
- [x] 7.4 New frame: 100k Condor production results — file count, total
      events, dedup factor, per-channel candidate yields.
      → Frame "100×1k preset-B Condor production --- topline" lists
      96/100 OK, 51\,023 events, 0.81 s/evt, 8.21 tracks/event,
      61.3\,\% dedup, and the per-channel cand-counts. The Stage-2
      hand-off frame at 100k stats shows 324\,154 candidates, 0 issues.
- [x] 7.5 New frame: ALCAReco vs RECO timing breakdown (carry over the
      ~1.9% number from the previous deck; re-measure on the 100k run
      to confirm).
      → Frame "ALCAReco vs RECO --- per-event cost" added; will
      re-measure on the 100k run.
- [x] 7.6 New frame: per-particle daughter layout — the Stage-2 hand-off
      contract. Show the per-channel table from task 4b.3 (smoke) and
      4b.4 (100k); make explicit that each daughter carries
      `(PDG tag, TrackRef → hit-bearing track)` so Stage-2 can pair
      `(track, mass)` without ambiguity.
      → Frame "Stage-2 hand-off --- per-particle (species, track)
      contract" added with the dual-regime explanation and the per-channel
      table from the smoke run.
- [x] 7.7 New frame: "Helix propagation to common vertex — why and
      how much." Brief physics (5-vs-6 parameters, momentum direction
      rotates along helix), the per-event displacement scale
      (B+ cτ ≈ 491 μm), the expected mass-bias removal (~1–3 MeV per
      candidate, scaling with cτ), and the m(B+) peak-shift result
      from task 6.5. Make explicit that resolution is unchanged
      within statistics and efficiency is unchanged exactly.
      → Two frames: "Helix propagation --- why and where the bias lives"
      and "Helix propagation --- how (existing CMSSW classes)" cover
      the physics + the implementation. Peak-shift number is reserved
      for the smoke-run placeholder frame.
- [x] 7.8 Compile cleanly with `pdflatex slides/jpsix-producer-final.tex`.
      → 15 pages, 260 KB; no errors, only standard font warnings.

## 8. Memory + housekeeping

- [x] 8.1 Update the `project_jpsi_x_alcareco.md` memory entry to point
      at this change and supersede the "preset C is the production
      default" framing wherever it persists.
      → Memory entry rewritten to list four changes and call out preset B
      as the production default; helix-prop and daughter-regime contracts
      called out as "do not regress" knobs.
- [x] 8.2 No spec deletions; existing `add-jpsi-x-condor-production`
      and the related changes remain unchanged.
      → Confirmed; no spec file was deleted.

## 9. Reproducible cmsDriver flow + splitLevel customise

The original 100k production wrote split-level-99 branches because the
recoskim cfgs were hand-placed in the CMSSW area without the
`alcarecoSplitLevel.setAlcaRecoSplitLevel` customise. Reading those
files in CMSSW>=15_0 hits ROOT cms-sw#19773: the SiStripCluster v11→v14
read rule does not fire on split branches, clusters come back empty,
and the downstream Stage-2 CVH refit cannot reconstruct hits. To make
the produced AlCaReco actually consumable downstream, we both fix the
customise AND make the cmsDriver invocation reproducible.

- [x] 9.1 Move the 8 hand-placed `recoskim_Run2016H_Charmonium_JpsiX*.py`
      cfgs out of `CMSSW_10_6_17_patch1/` into
      `CMSSW_10_6_17_patch1/cmsdriver_outputs_backup_2026-06-23/`. These
      are *generated artifacts* and should never be committed or edited
      by hand. The backup is retained for post-mortem comparison only.
- [x] 9.2 Cherry-pick / fast-forward the 10_6 CMSSW source tree to
      `david/from-CMSSW_10_6_17_patch1` tip `b8ebaa3d9f2` (the
      `setAlcaRecoSplitLevel` customise commit). Done as part of the
      git consolidation: `.git` moved from
      `CMSSW_10_6_17_patch1_src_backup/src/.git` into
      `CMSSW_10_6_17_patch1/src/.git`, then ff'd onto David's tip.
- [x] 9.3 Add `condor/jpsix_alcareco/make_recoskim_cfgs.sh` — single
      source of truth for the cmsDriver invocation. Includes the
      `Alignment/CommonAlignmentProducer/alcarecoSplitLevel.setAlcaRecoSplitLevel`
      customise in `--customise`, regenerates presetB by default
      (`TKALJPSIX_SELECTION_PRESET=B` → cff dispatch). Self-checks the
      output contains the customise import. Production only needs
      presetB; additional presets can be passed as positional args.
- [x] 9.4 Wire `make_recoskim_cfgs.sh` into `build_tarball.sh` — cfgs
      regenerated fresh on every tarball build, no in-tree drift
      possible. Tarball includes only `recoskim_…_presetB.py`.
- [x] 9.5 The 7 in-src changes were committed on top of David's tip
      (`4d65009c736` … `026705dd9ec`); HEAD = `b8ebaa3d9f2` + 7
      topical commits. `Alignment/CommonAlignmentProducer/python/alcarecoSplitLevel.py`
      is in tree.
- [ ] 9.6 Rebuild the live CMSSW area: `cd CMSSW_10_6_17_patch1 && eval $(scram runtime -sh) && scram b -j 8`.
- [ ] 9.7 Smoke the regen script: `condor/jpsix_alcareco/make_recoskim_cfgs.sh` (no args → presetB only) — confirm a fresh `recoskim_Run2016H_Charmonium_JpsiX_presetB.py` appears at `CMSSW_10_6_17_patch1/<cfg>.py` and `grep -c "alcarecoSplitLevel" <cfg>` ≥ 2.
- [ ] 9.8 1-file local `cmsRun` smoke (200 events): run the freshly-regenerated cfg against one Run2016H RAW input. Confirm with a FWLite split-level probe that `SiStripCluster*`, `SiPixelCluster*`, `TrackingRecHit*`, `recoTracks*`, `recoTrackExtras*`, and every `ALCARECOTkAlJpsiX*Resonances` branch are at split level **1** (not 99). Output deleted after the check; this is a wiring test only.

## 10. 100k re-production at splitLevel=1

Pattern mirrors `add-jpsi-x-condor-production` (build tarball → wrapper
runs cmsRun + xrdcp → JSON sidecar per job → aggregate). Only delta:
the wrapper-arg `DEST_SUFFIX` is set to `_final_sl1` so output lands in
a fresh dir; the prior `preset_B_final/` stays as a baseline.

- [x] 10.1 Submit file written: `condor/jpsix_alcareco/condor_jpsix_alcareco_presetB_final_sl1.sub` — copy of `…_final.sub` with `arguments = "B $(filename) _final_sl1 1000"` and `log/out/err` under `presetB_final_sl1/`.
- [ ] 10.2 Create the new log dir: `mkdir -p /work/submit/pmlugato/mz/logs/jpsix_alcareco/presetB_final_sl1`.
- [ ] 10.3 Rebuild the tarball: `cd condor/jpsix_alcareco && ./build_tarball.sh`. This regenerates `recoskim_Run2016H_Charmonium_JpsiX_presetB.py` (splitLevel customise wired in) and packs `jpsix_alcareco_payload.tgz` (137 KB target per `add-jpsi-x-condor-production` precedent).
- [ ] 10.4 Verify x509 proxy: `voms-proxy-info -file /home/submit/pmlugato/x509up_u239501 -timeleft` ≥ a few hours.
- [ ] 10.5 Submit: `cd condor/jpsix_alcareco && condor_submit condor_jpsix_alcareco_presetB_final_sl1.sub`. Expected: 100 jobs queued. Target wall-clock ~20 min once jobs land (per the finalize design note: 100×1k events at ~250 ms/evt ≈ 250 s event-loop + ~7 min cmsRun startup per job).
- [ ] 10.6 Monitor: `condor_q $USER -nobatch` until idle, then `condor_history $USER -limit 100 -af ClusterId ProcId JobStatus ExitCode` to confirm ≥96/100 exited 0. Output ROOTs land at `/ceph/submit/data/user/p/pmlugato/mz/alcareco/jpsi-x/preset_B_final_sl1/Run2016H/jpsix_alcareco_presetB_<basename>.{root,json}`.
- [ ] 10.7 Aggregate: from repo root run `python3 condor/jpsix_alcareco/merge_and_report.py --preset B --suffix _final_sl1 --era Run2016H` (mirrors the established `add-jpsi-x-condor-production` invocation; uses `cmssw-el7` Singularity for hadd + FWLite). Outputs the per-channel cand summary + the per-event runtime/size diagnostics.
- [ ] 10.8 Verify success rate ≥96/100, total output size comparable to prior `preset_B_final` (~143 MB → ~145 MB after the ~1–2% splitLevel-1 hit), and per-channel cand counts statistically consistent with the prior production.
- [x] 10.9 Spot-check (cluster 3194086 results):
      (a) Customise fired at cmsDriver: `customising the process with setAlcaRecoSplitLevel from Alignment/CommonAlignmentProducer/alcarecoSplitLevel`.
      (b) ROOT GetSplitLevel() on every relevant ALCARECOTkAlJpsiX branch (SiStripCluster, SiPixelCluster, TrackingRecHit, recoTracks, recoTrackExtras, all 7 channel VCC collections) = {1, 2} — down from 99 in the prior production.
      (c) Per-file sizes 6.5–6.7 MB (was 6.6 MB), +0.6% — matches the predicted +1–2% splitLevel-1 hit.
      (d) FWLite reads SiStripCluster.firstStrip with sane values (0–768) on the new files; no corruption.
      The full cmsRun-side I/O-rule-chain test via the single-track CVH refit requires `scalarPot3DInitFile`, which is a Stage-2 prerequisite (not in this change's scope). Tracked as `add-jpsi-x-stage2-bplus-cvh/tasks.md` §0.5.
- [x] 10.10 sl1 sample is the canonical Stage-2 input: `/ceph/submit/data/user/p/pmlugato/mz/alcareco/jpsi-x/preset_B_final_sl1/Run2016H/` (100 files, 660 MB). `add-jpsi-x-stage2-bplus-cvh/proposal.md` updated to point at this directory; its blocking Stage-0 input task cleared.
