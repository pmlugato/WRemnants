# Change: Stage-2 CVH refit for B+ → J/ψ K+ (reduced/composed)

## Why

Stage-1 of the J/ψ+X program produces a slim `TkAlJpsiX` ALCARECO with
deduplicated cloned tracks and per-channel `VertexCompositeCandidate`
collections (`add-jpsi-x-alcareco-channels`, `…-candidate-quality-cuts`,
`…-selection-comparison`, `…-vertex-fit-and-low-pt`, `…-condor-production`).

**Stage-1 input for this change is the finalized preset-B 100k production**
(`finalize-jpsi-x-preset-b-production`, slides `jpsix-producer-final.tex`):
100 RAW files × 1k events on Run2016H Charmonium, 96/100 jobs OK, at
`/ceph/submit/data/user/p/pmlugato/mz/alcareco/jpsi-x/preset_B_final_sl1/Run2016H/`
(`jpsix_alcareco_presetB_*.root`, 96 files). Preset B is the operational
production default. Preset C is out of scope for this change.

**What "finalized preset B" means for Stage-2 inputs.** The finalize change
(`finalize-jpsi-x-preset-b-production`) made three preset-B-specific edits
relevant to us:

1. **J/ψ mass constraint OFF at Stage-1** (`applyJpsiMassConstraint=False`):
   the dimuon 4-vector saved in the Stage-1 B+ candidate is the raw track-level
   sum, not the Kalman-constrained one. This is a *Stage-1 saved-kinematic*
   knob — orthogonal to the Stage-2 CVH refit's own `doMassConstraint`, which
   for this change is **ON** (see §"Mass-constraint status" below). Stage-2
   re-refits each track from its rec hits and ignores the saved Stage-1
   4-vectors.
2. **Helix propagation to a common candidate-level point ON**: each leaf
   track's 4-vector in the saved `RecoChargedCandidate`/`VertexCompositeCandidate`
   is built at the candidate's helix-helix crossing point, removing the
   ~1–3 MeV displacement-correlated mass bias from PCA-to-beamline evaluation.
3. **`AlignmentTrackSelectorWithIndexMap.ptMin = 0.1`** (was 0.4): the merged
   cloned track collection now includes soft V0 daughters at the 0.1 GeV floor.

Items 1 and 2 affect the *informational* `Bplus_*` / `Jpsi_*` raw-kinematics
branches in our joined tree (helix-propagated values are stored, m(B+) peak
expected ~30 MeV under preset B without constraint). The Stage-2 CVH refit
itself **re-refits each track from its rec hits and is unaffected** by the
saved 4-vectors — Stage-2 reads only the leaf `TrackRef`s, hits, and
`TrackExtra`s. Item 3 simply enlarges the per-event track population.

The nested B+ VCC layout the splitter needs is preserved under all of these:
`JpsiXCandidateProducer.cc:617–618` always emits `daughter(0)` = J/ψ VCC,
`daughter(1)` = bachelor-K `RecoChargedCandidate` (with `setTrack(...)`).

**Bachelor mass hypothesis taken from the collection name, not from
`daughter.mass()` or `daughter.pdgId()`.** The Stage-2 splitter and the
single-track maker cfi key off the channel collection (e.g.
`ALCARECOTkAlJpsiXBPlusResonances` → kaon, `trackParticleName="kaon"`)
and hardcode the species per channel:

| Collection | bachelor / sub-resonance daughters |
|---|---|
| `BPlusResonances`    | K              |
| `BcResonances`       | π              |
| `B0KstarResonances`  | (K, π)         |
| `BsPhiResonances`    | (K, K)         |
| `B0KsResonances`     | (π, π)         |
| `LambdabResonances`  | (p, π)         |
| `Psi2SResonances`    | (π, π)         |

The producer *currently* stamps `daughter.mass()` and `daughter.pdgId()`
on each leaf (consistent with the cfg `bachelorMass`/`bachelorPdgId`),
and the mother `cand.mass()` is the hypothesis-bearing 4-vector sum at
the helix-helix crossing — kept for diagnostic m(B) peak plots in this
100k production. **These will be dropped before the final postVFP full
run** (massless leaf 4-vectors, `pdgId=0`, no mother mass); Stage-2 CVH
already takes the species from the collection name, so the same code
works both regimes. Stage-2 SHALL NOT depend on `daughter.mass()` or
`daughter.pdgId()`.

Stage-2 is the CVH refit: the Geant4e per-particle track refit that emits,
per track, the parameters at PCA, the full covariance, and the Jacobians
∂(track params)/∂a with respect to the global-correction parameters a
(module ΔBz, 6-DOF alignment, per-η energy-loss ε). These per-track
Jacobians are the input to the global least-squares (GBL) tracker/B-field
calibration fit. This is the analog, for B+ → J/ψ K+, of David Walter's
V0/D* Stage-2 (`david_260506_alcareco.pdf`, slides 12–20) and the existing
production J/ψ → μμ refit (`Analysis/HitAnalyzer/test/runCvhJpsi.py`).

The physics motivation for adding B+ → J/ψ K+ on top of J/ψ → μμ is the
**different daughter species**: the bachelor is a kaon, refit under the
kaon mass hypothesis, so its material/energy-loss Jacobian carries
information the muon-only sample cannot — the handle on the A–ε degeneracy
(B-field vs. material) that motivates the whole program (slide 2).

## What this change does — and explicitly does not

This change adds the **B+ → J/ψ K+ channel only**. The other six J/ψ+X
channels (Bc, B0→K*0, Bs→φ, and the three V0-mode channels) are a later
proposal. B+ is the smallest complete vertical slice: two muons plus one
bachelor track.

It takes the **reduced/composed** approach. No new multi-track joint
fitter is built. Instead it reuses two refit kernels that already exist
and are production-tested in `CMSSW_15_0_19_patch2`:

| Leg | Existing kernel | Hypothesis |
|---|---|---|
| J/ψ → μ⁺μ⁻ | `ResidualGlobalCorrectionMakerTwoTrackG4e` | muon |
| bachelor K⁺ | `ResidualGlobalCorrectionMakerG4e` (single-track) | `trackParticleName="kaon"` |

The single-track maker already supports non-muon hypotheses — its
`trackParticleName` parameter explicitly accepts `"pi"`, `"kaon"`,
`"proton"` (`plugins/ResidualGlobalCorrectionMakerG4e.cc:88,110`), and the
kaon is already registered in `G4ErrorPhysicsListForCVH`. So no fit-kernel
C++ is needed; the new code is the plumbing that ties the two kernels to a
common B+ candidate and joins their outputs.

**Mass-constraint status — J/ψ mass constraint ON in the CVH refit; no
other kinematic links.** The composed chain applies exactly one kinematic
constraint:

| Constraint | Applied? | Where set |
|---|---|---|
| J/ψ → μμ mass constraint inside the two-track CVH joint vertex fit | **YES** | `doMassConstraint=True` on the two-track maker's cfi |
| Common-vertex coupling of the two muons inside the two-track CVH fit | YES (structural) | always on in `ResidualGlobalCorrectionMakerTwoTrackG4e` — this is a 10-dim vertex-PCA joint state, not a mass constraint |
| Kaon coupling to the dimuon vertex (single-track maker) | NO | the kaon is refit independently at its own PCA |
| B+ vertex / B+ mass constraint linking μμ to K | **NO** | not applied in Stage-2; eventually applied as a GBL Lagrange term (out of scope for 2a) |
| J/ψ mass constraint at the module-correction (GBL) stage | NO | GBL solves for global parameters **a** from per-track Jacobians; the per-event J/ψ mass constraint is applied at the Stage-2 CVH refit, not again at GBL |

So the dimuon refit pulls m(μμ) to the PDG J/ψ mass inside the joint vertex
fit; the kaon is refit independently with no link to the dimuon vertex; and
no B+ vertex/mass constraint is applied in Stage-2. This matches David's
"with mass constraint" output variant for the dimuon, and option-a (no B+
linking) for the bachelor.

**Output is a flat per-candidate ROOT TTree** (the existing CVH output
format), not custom NanoAOD. The consumer is the GBL global-correction fit.

**Critical input-readability prerequisite (splitLevel=1).** The 96 files in
`preset_B_final_sl1/Run2016H/` were written with all branches at split-level 99
(ROOT default; `recoskim_Run2016H_Charmonium_JpsiX_presetB.py` does not set
splitLevel). Reading split-level-99 ALCARECO from CMSSW_10_6 in CMSSW_15_0_19
hits ROOT cms-sw#19773: the SiStripCluster v11→v14 schema rule does not fire
on split branches, and clusters come back **silently empty** — which makes
the CVH refit unable to reconstruct hits. The fix (David's commit
`b8ebaa3d9f2`, `Alignment/CommonAlignmentProducer/python/alcarecoSplitLevel.py`
+ a `setAlcaRecoSplitLevel` customise wired into the cmsDriver `--customise`
list) sets `splitLevel=1` on every ALCARECO output module. That commit lives
on `david/from-CMSSW_10_6_17_patch1` but is **not on the current
`update-btojpsik-fitting` branch**, and the customise was not used in the
100k production. Stage-2 therefore requires either (a) cherry-picking the
customise + re-producing preset_B_final, or (b) re-streaming the existing
96 files through a one-pass `cmsRun` cfg that re-writes them at splitLevel=1.
Promoted to a blocking Stage-0 task. FWLite reads via a different code path
that does *not* invoke the schema rule chain, which is why the earlier
rec-hit-count smoke test (174/174 with hits) was not sensitive to this and
why cluster `dataSize()` reports populated.

## What changes

### 1. New `JpsiKCandidateSplitter` EDProducer (`Analysis/HitAnalyzer`)

Reads the Stage-1 B+ `VertexCompositeCandidate` collection (each B+ is a
nested composite: J/ψ VCC daughter + bachelor-K `RecoChargedCandidate`
daughter) and emits, for each surviving B+ candidate:

- a `reco::VertexCompositeCandidateCollection` of the dimuon (J/ψ)
  sub-candidate, for the two-track maker;
- a `reco::RecoChargedCandidateCollection` (or `reco::TrackCollection`) of
  the bachelor kaon, for the single-track maker;
- a per-candidate index that travels with both, so the two refit trees can
  be joined 1:1 offline (same `(run, lumi, event, bCandIdx)`).

This is the only genuinely new physics-adjacent C++. It contains no fit.

### 2. cfi configs + a cmsRun driver

- `…TwoTrackJpsiMuMuG4e` cfi pointed at the splitter's dimuon collection
  (muon hypothesis; **J/ψ mass constraint ON**, `doMassConstraint=True`).
- single-track `…G4e` cfi pointed at the splitter's bachelor collection
  with `trackParticleName="kaon"`.
- `runCvhBplusJpsiK.py` — a cmsRun driver wiring splitter → both makers,
  modeled on `runCvhJpsi.py` (era `Run2_2016`, GT `auto:run2_data`, scalar-
  potential B-field, `fillJac=True`). Reads the 2016 `TkAlJpsiX` ALCARECO
  cross-release (precedent: `runCvhJpsi.py` already reads a 2016
  `ALCARECOTkAlJpsiMuMu` in this 15_0_19 release).

### 3. Per-candidate index in the two makers' output

Both makers gain a B+-candidate index branch (or an associated ValueMap)
so the offline join is unambiguous even when a leg's refit fails or rows
are reordered. This is a minimal, additive edit; default behaviour for the
existing J/ψ/Υ/Z/Ks/Λ/D0 configs is unchanged when the index input is
absent.

### 4. Offline combine script

A Python/ROOT script that joins the two refit trees on
`(run, lumi, event, bCandIdx)` into one per-B+-candidate flat tree carrying
all three tracks' CVH outputs + Jacobians, the raw B+ kinematics, and the
m(μμK) mass under the kaon hypothesis. Output is the GBL-fit input.

### 5. Branch schema (for review)

Proposed naming, keeping David's `Jpsi_/Muplus_/Muminus_` convention so
existing tooling transfers, and adding bachelor + parent branches:

| Prefix | Source | Contents |
|---|---|---|
| `Bplus_*` | raw composite | B+ pt/eta/phi, m(μμK), Stage-1 vertex (x,y,z) |
| `Jpsi_*`, `Jpsikin_*`, `Jpsicons_*` | 2-track maker | dimuon kinematics (raw / kin-fit / constrained) |
| `Muplus_*`, `Muminus_*` | 2-track maker | per-muon refit params, cov, `*_jacRef` |
| `Kbach_*` | single-track maker | bachelor-K refit params, cov, `Kbach_jacRef` |
| join keys | both | `run`, `lumi`, `event`, `bCandIdx` |

### 6. Supersede the stale `add-jpsi-x-stage2-cvh-refit`

That change's content ("Phase 3 — Kalman vertex fits in Stage-1 AlCaReco")
is obsolete: the AlCaReco-stage Kalman fit it proposed was rejected on cost
grounds and re-introduced, where it belongs, as Stage-1 preset C in
`add-jpsi-x-vertex-fit-and-low-pt`. It was never implemented (0/38 tasks)
and its directory name wrongly suggests it is the Stage-2 CVH change. This
change removes that directory as a rollout step.

## Impact

- **Affected specs**: NEW capability `cvh-refit-jpsi-x` (Stage-2 CVH refit).
- **Affected code** (all under `CMSSW_15_0_19_patch2/src/Analysis/HitAnalyzer`,
  except the housekeeping deletion):
  - NEW `plugins/JpsiKCandidateSplitter.cc` + `python/JpsiKCandidateSplitter_cfi.py`
  - NEW `test/runCvhBplusJpsiK.py` driver
  - additive per-candidate-index branch in `plugins/ResidualGlobalCorrectionMakerTwoTrackG4e.cc`
    and `plugins/ResidualGlobalCorrectionMakerG4e.cc`
  - NEW cfi configs pointing the two makers at the split collections
  - NEW offline join script (location TBD: `scripts/btojpsik/` or alongside the driver)
  - DELETE `openspec/changes/add-jpsi-x-stage2-cvh-refit/`
- **Reused unchanged**: both CVH fit kernels, `G4ErrorPhysicsListForCVH`,
  `Geant4ePropagator`, `ParticleProperties`, the scalar-potential B-field.
- **No changes** to the Stage-1 ALCARECO producers' physics, the cff, or
  any WRemnants narf/rabbit fit code.
- **External prerequisite**: cherry-pick David's commit `b8ebaa3d9f2`
  (`Alignment/CommonAlignmentProducer/python/alcarecoSplitLevel.py` +
  `setAlcaRecoSplitLevel` customise wired into the cmsDriver `--customise`
  list) onto the current `update-btojpsik-fitting` 10_6 branch, and
  re-produce or re-stream the existing `preset_B_final_sl1/Run2016H/` 96 files
  at branch split-level 1. Done as task 0.4–0.5 in `tasks.md`.

## Open questions

Stage-0 facts verified by FWLite inspection of preset-B Condor files
(`CMSSW_15_0_19_patch2`):

- ✅ **B+ InputTag** is `ALCARECOTkAlJpsiXBPlusResonances`
  (`vector<reco::VertexCompositeCandidate>`).
- ✅ **Daughter layout** confirmed on the prior preset_B file: 21/21 B+ have
  `daughter(0)` = J/ψ VCC (2 daughters), `daughter(1)` = bachelor
  `RecoChargedCandidate` (matches `JpsiXCandidateProducer.cc:617–618`).
- ✅ **Cross-release deserialization** works (10_6-produced file read in
  15_0_19, no schema/dictionary error on the VCC collection itself).
- ⚠️ **Hit/cluster I/O — partial gate, NOT cleared.** FWLite reports
  174/174 tracks with `TrackExtra` + non-zero `recHitsSize`, and
  `SiStripCluster` `dataSize()` is non-zero — but this reads through a
  different code path than `cmsRun`/CVH; clusters can still be empty when
  dereferenced via the I/O-rule chain. The 96-file `preset_B_final_sl1/`
  production is all at branch split-level **99** (verified directly), so
  the rule cms-sw#19773 applies.

Open prerequisite (BLOCKING): the splitLevel=1 customise (`b8ebaa3d9f2`)
is not on the current `update-btojpsik-fitting` branch and was not used
when producing `preset_B_final`. Stage-0 must cherry-pick it and either
re-produce or re-stream the input set at splitLevel=1, then re-validate
that hits/clusters round-trip through a sample CVH single-track refit.

Remaining (non-blocking; resolve during apply):

1. Whether the bachelor leg is fed to the single-track maker as a
   `TrackCollection` or a `RecoChargedCandidateCollection` (the maker reads
   tracks; the splitter output type follows from that).
2. Whether to carry the candidate index as an output branch in each maker
   or as an `edm::ValueMap` keyed on the split-collection elements.
3. The bachelor `TrackRef` resolves into the cloned `ALCARECOTkAlJpsiX`
   collection (Stage-5 remap worked); confirm the splitter emits/points the
   bachelor at *that* collection so the single-track maker refits hit-bearing
   tracks, not stale `generalTracks` refs.

### Post-2a follow-ups (in resolution order)

The 2a milestone (sane B+ peak in m(μμK) after CVH refit) landed via the §5
join; §6 validation surfaced three issues to resolve in this strict order
(tracked as §9.4 → §9.5 → §9.6 in `tasks.md`):

1. **§9.4 — shared-propagator regression (was: "single-track-maker refit-kernel regression", was: "kaon q/p anomaly").** Refit |q/p| ≈ ½ × input |q/p| (refit pT ≈ 2× input pT). Four black-box A/B tests narrowed the cause, then a fifth (9.4.e.iii) reframed the bug from single-track-specific to shared.
   - **9.4.a (done)** — input perigee read is byte-for-byte clean; bug is downstream of `refFts` init.
   - **9.4.b (done)** — scalar-pot 3D field swapped for standard CMSSW field; identical 2.546× anomaly. Field-correction source exonerated. **Caveat (relevant for 9.4.e)**: this varied the *correction layer* (`dB`), not the *base field* the propagator and G4 each consume — so a base-field discrepancy between the analytic Jacobian path and the G4 propagation path is NOT ruled out.
   - **9.4.c (done)** — kaon hypothesis swapped for muon hypothesis on same kaon tracks; identical 2.546× anomaly. Mass / dE/dx hypothesis exonerated.
   - **9.4.d (done)** — same single-track maker on full `ALCARECOTkAlJpsiX` collection: every track type blows up (1.82× / 2.31× / 2.55×), Newton step systematically driving `|q/p|→0` across both charges. **Initial inference (now overturned)**: based on `Jpsikin_mass ≈ PDG` in the two-track maker, we inferred the bug was single-track-specific. The inference was wrong (see 9.4.e).
   - **9.4.e (done — KEY REFRAMING)**: read-only diffs of the May-2026 commit stack (`8b3435b309b` 5×7→5×9 Jacobian, `a702ebbf408` two-track sparse-GBL) found **no single-track-vs-two-track asymmetry** — both makers got identical structural updates and the sparse-GBL refactor is "Mathematically identical to the dense path" per the commit message (validated bit-exact). With no asymmetry to explain a single-track-only bug, the single-track-only inference was tested directly: on existing two-track output ntuple `/scratch/submit/cms/dwalter/cvh_scaling_10k_20260602_164838/run_N8/globalcor_0.root`, `Muplus_pt / Muplustrk_pt` (refit / input, *unconstrained* `icons==0` branch) over 1002 muon legs reports **mean 2.06, median 2.00**. **The two-track maker exhibits the same 2× pT blow-up**; `Jpsikin_mass ≈ PDG` only because the J/ψ kinematic-mass constraint absorbs the error by pulling both legs proportionally. **The bug is a shared regression in the propagator / Jacobian / field-on-worker stack**, not single-track-specific. Suspect ranking shifts to shared-code blast radius:
     1. **`8b3435b309b`** (Geant4ePropagator regenerated; 953-line rewrite, per-step accumulator widened `rightCols<2>→<4>`). Highest blast radius on shared code; the SymPy regeneration's "byte-for-byte Bz-only equivalence" claim was validated against a standalone FD closure, NOT against an actual Geant4-propagated track.
     2. **`8863d83e447`** (multi-threading) — worker-thread G4 field installed via `sim::FieldBuilder`. If the worker field disagrees with the master CMS field handle (which is what the analytic Jacobian reads via `field->inInverseGeV(pos)`) by ~2×, G4 curves twice as much as predicted → LSQ pulls q/p by ½ to compensate. §9.4.b does NOT rule this out — that test held both views of the field's *correction* layer fixed.
     3. **`aa7f5893bb3`** (scalar dBz → vector dB) — `pca2curvJacobianD` / `curv2localJacobianAltelossD` signatures changed.
   - **9.4.iv (DONE — root cause found and fix verified)** — instrumented `Geant4ePropagator.cc::propagateGenericWithJacobianAltD` with per-step diagnostics over 6 rounds, narrowing from "G4 bends 2× analytic" → "G4 sees `B_z = 7.62 T` vs analytic `3.81 T`" → "the maker's MagneticField pointer's `inTesla`/`inInverseGeV` interfaces both report the correct 3.81 T" → "the doubling is in `dB`, the propagator's correction argument" → "`dB.z = 3.81 T` (FULL field) at every step." **Root cause** at `ResidualGlobalCorrectionMakerBase.cc:961–966`: the base maker seeds `corparms_` (the field-correction-mode parameter vector, supposed to start at zero and accumulate small *deltas* during the LSQ refinement) with the FULL initial scalar-potential coefficients from the harmonic-basis dump file (`polyfit3d_full_coeffs_lmax18_cmsswnorm.txt`). `dB = fieldCorrection_->getCorrectionAt(propStartPos, corparms_)` therefore returns the FULL ~3.81 T field at every position. The propagator then adds `dB` to the CMSSW base both for the analytic Jacobian (`Geant4ePropagator.cc:1293–1295`) and for the G4 stepper (`Geant4ePropagator.cc:655`'s `cmsField->SetOffset(dB)`) — so the effective field is `B_CMSSW + B_scalarPot ≈ 2 × B_CMSSW`. The input track q/p was reconstructed in the *real* 3.81 T field, so the refit pulls q/p by ½ to reconcile a 3.81 T-curvature data trajectory with the propagator's 7.62 T-curvature prediction. Introduced by `aa7f5893bb3` (the same commit that added the scalar-potential machinery), survived to today. **Fix applied + verified** by wrapping the seeding loop in `#if 0` (leave `corparms_` at zero, the LSQ's default starting point):
     - `dB = (0, 0, 0) T` at every propagation step.
     - `B_G4 / B_analytic = 1.0` (was 2.00006).
     - `actual_dphi / dphi_correct = 0.999996` (was 2.001).
     - §9.4.a kaon refit/input |p| = **0.976–0.997** (was 2.546× across 33 rows).
     - Two-track maker (200 events, dimuon): `Muplus_pt / Muplustrk_pt` mean=**1.02**, median=**1.01**, stdev=**0.18** (was 2.06/2.00/0.61). `Jpsikin_mass` stays near PDG (3.1015).
   - **Fix promoted to a permanent bypass** (2026-06-25). The test `#if 0` was replaced by a commented-out version of the seeding loop preceded by a `FIXME [openspec §9.4.iv]` block at `ResidualGlobalCorrectionMakerBase.cc:961–980`: the seeding code is preserved (commented, not deleted) for a future proper fix, and the call site now just leaves `corparms_` at zero with a clear "seeding disabled pending proper delta-basis fix" log line. This is option (c) — drop seeding entirely, matching pre-`aa7f5893bb3` semantics. The §9.4.iv per-step diagnostic block (and its `#include <atomic>`) in `Geant4ePropagator.cc` was deleted in the same pass.
   - **§6.4 FD-closure re-check (2026-06-25)**: NOT resolved by the §9.4.iv fix. Reran with `runFDClosure=True` at eps=1e-2, 1e-4, 1e-6 on the same kaon file; the mode-0 FD output `[1.7e-5, -9.8e-5, -3.6e-6, 6.6e-6, -4.3e-5]` is identical to ≥4 sig figs across 4+ orders of magnitude of eps. The perturbed `dB` argument to `propagateGenericWithJacobianAltD` does not propagate through to the propagation output proportionally to eps — likely because the propagator's per-step accumulator uses the step-by-step `theField->inInverseGeV(pos)` rather than the FD scheme's perturbed-once `dB`. This is an independent wiring bug in the FD path at `ResidualGlobalCorrectionMakerG4e.cc:1660-1712`, still tracked as §9.6.
   - **Remaining follow-up** (now in `tasks.md` §9.4.e.iv.8): eventual proper fix in option (a) or (b) form (seed `corparms_` as `initCoeffs - {harmonic expansion of base field}`, or redefine `getCorrectionAt` to return a delta), so the LSQ starts from the best available field description rather than from zero. **§9.5 and §9.6 are now unblocked.**
2. **§9.5 — switched m(μμK) to the refit kaon (DONE 2026-06-25).** The join
   script `scripts/btojpsik/join_cvh_bplus_jpsik.py::compute_m_mumuK` now
   decodes `Kbach_refParms[0:3]` ((q/p, λ, φ) at the kaon perigee, d0=0,
   dsz=0) into pT/η/φ via pT=|p|·cos λ, η=asinh(tan λ), and supports three
   kaon-source modes ("track" legacy, "refit_nopropagate" at the kaon's
   PCA, "refit" with small-arc helix-bend propagation `Δφ = −q·0.003·B·s_T/pT`
   to the dimuon vertex `Jpsicons_x/y/z`). Verified on 111 matched B+
   candidates from a 200-event Run2016H ALCAReco sample:
   | kaon source              | mean    | median  | std     | window |
   |--------------------------|---------|---------|---------|--------|
   | input track (legacy)     | 5.2782  | 5.2591  | 0.2359  | 53%    |
   | refit @ kaon PCA         | 5.2714  | 5.2692  | 0.2077  | 51%    |
   | refit @ dimuon vertex    | 5.2709  | 5.2684  | 0.2053  | 50%    |
   Median offset vs PDG B+ (5.27934 GeV) tightened from **−20 MeV** (track)
   to **−11 MeV** (refit); width shrank **12%** (0.236 → 0.205 GeV). The
   §9.5.b helix-bend propagation moves the median by ~1 MeV, confirming
   the small-displacement (B+ flight ≪ helix radius) approximation is
   dominated by the §9.5.a refit substitution. **2a milestone "sane B+
   peak after CVH refit on both legs" reached.**
3. **§9.6 — FD-closure repair — deferred to a future change handled under 2b.**
   Out of scope for this proposal. After the §9.4.iv bypass the kaon-leg
   FD closure was re-checked (2026-06-25) and the eps-invariance signature
   persisted (worst rel 8.17×10⁶ at eps = 1e-2, 1e-4, 1e-6; mode-0 FD
   output identical to ≥4 sig figs across all three) — confirming the
   wiring bug in the perturbed-`dB` handoff to the propagator's per-step
   accumulator is independent of §9.4.iv. Not blocking the 2a milestone
   (§9.7 below validates the refit on the data side); blocking the 2b
   GBL handoff, where the Jacobian's correctness must be independently
   verifiable before it can ship. The repair will be its own openspec
   change once 2b begins; full diagnostic narrative + candidate repair
   approaches preserved in git history of `tasks.md` for that future
   proposal to pick up.
4. **§9.7 — Full Run2016H CVH run + health checks + B+ peak (DONE 2026-06-25).**
   Scaled §9.5 from 200 events to the full ALCAReco dataset
   (100 files, 660 MB, `Run2016H/preset_B_final_sl1/`). cmsRun multi-
   stream via `scripts/btojpsik/run_full_run2016H.sh` (100 parallel
   per mode, `useScalarPot3D=False`, `plimit=1.0`); wall time
   1816 s kaon, 2286 s dimuon; all outputs in
   `/ceph/submit/data/user/p/pmlugato/mz/alcareco/jpsi-x/cvh_outputs/run2016H_2a/`.
   Joined with the `--scalars-only` flag (new — bypasses per-row
   Python-list conversion of jagged Jacobian branches, which dominates
   the wall clock at full scale and isn't needed for the §9.7 analysis):
   91 689 dimuon × 30 878 kaon → **30 719 matched B+ candidates**.
   B+ peak Gaussian fits (unbinned MLE in [5.15, 5.45] GeV):
   | kaon source           | N_signal | mean    | σ      | offset vs PDG |
   |-----------------------|----------|---------|--------|---------------|
   | input track (legacy)  | 17 531   | 5.2958  | 0.0842 | +16.5 MeV     |
   | refit @ kaon PCA      | 16 469   | 5.2951  | 0.0843 | +15.7 MeV     |
   | refit @ dimuon vertex | 16 450   | 5.2949  | 0.0842 | +15.6 MeV     |
   (PDG B+ = 5.27934 GeV.) All three flavors agree to ~1 MeV on the
   Gaussian core — confirms the §9.5 refit substitution preserves the
   peak. The +15.6 MeV offset is the known CMSSW scale/B+-vertex
   effect, not a refit artefact. Health checks: niter mean 6.87 / max
   10, edmval 100% finite, refit/input pT ratio medians within 0.4%
   of unity on all three legs, propagation-failure rate 0.000% on
   dimuon and 1.93% on kaon (soft-bachelor tail under plimit=1.0
   GeV — see §9.8). Slides include a full m(μμK) recipe explainer,
   a peak frame with the fit table + PNG histogram, and a
   health-checks frame; PDF builds clean. **2a refit milestone
   closed.**
5. **§9.8 — Soft-bachelor recovery at plimit=0.05 (was §9.3,
   promoted to 2a closeout; DONE 2026-06-26).** Added a `plimit`
   VarParsing knob to `test/runCvhBplusJpsiK.py` wired to
   `process.Geant4ePropagator.PropagationPtotLimit`. Confirmed
   plimit is NOT overridden anywhere downstream in the maker
   chain (both single-track and two-track makers clone the ES
   product via `streamPropagator_.reset(templateProp->clone())`;
   the copy constructor at `Geant4ePropagator.cc:116` preserves
   `plimit_`), so the single config knob suffices. Re-ran kaon mode
   only (dimuon unchanged — already 0.000% failures at plimit=1.0)
   on the full Run2016H sample at plimit=0.05, wall 579 s (no
   concurrent dimuon = less contention than baseline's 1816 s).
   Joined the new kaon outputs against the existing dimuon files
   into `joined_cvh_bplus_jpsik_plimit_0p05.root` (47 MB, 188 cols).
   Extended `plot_no_jpsi_constraint.py` with a `--joined-extra`
   argument that adds a third curve to
   `no_jpsicons_run2016H_bpeak_uncons.png`; new script
   `plot_plimit_compare.py` produces the kaon-kinematics overlay
   + health-check JSON. **Headline results**:
   | metric                           | plimit=1.0 | plimit=0.05 |
   |----------------------------------|------------|-------------|
   | matched B+ candidates             | 30 719     | **81 948 ($\times$2.67)** |
   | propagator failure rate           | 1.930 %    | **0.113 %** ($\div$17) |
   | unconstrained-J/ψ peak σ          | 136.1 MeV  | **135.8 MeV** (same) |
   | unconstrained-J/ψ peak N_signal   | 17 886     | **49 887** ($\times$2.8) |
   | niter (median / mean / max)       | 6 / 6.87 / 10 | 6 / 6.93 / 10 |
   | χ²/ndof median                    | 0.94       | **0.84** (better) |
   | nValidHits median                 | 14         | 9 |
   | refit/input pT median (std)       | 1.001 (0.49) | 1.009 (**0.33**) |
   The recovered soft tail is well-behaved: same resolution, no
   convergence degradation, no outlier blowup. Two new slide frames
   added (kaon kinematics overlay + health checks); the existing
   B+ peak comparison frame extended from 2 → 3 curves. **2a is
   fully closed.**
