# Change: Finalize CVH producer in CMSSW_15_0_19_patch2 — adopt AN-canonical `useIdealGeometry` + `corFiles` and complete alcareco-producer validation

## Why

The J/ψ→μμ dimuon-mass broadening investigation
(`slides/cvh_vs_alcareco_jpsi_mass_investigation.pdf`, memory
`project_cvh_jpsi_mass_broadening.md`) initially traced the σ~80 MeV,
33% catastrophic-tail failure of the CVH joint refit to
`useIdealGeometry=True`, and demonstrated that flipping to `False`
recovers σ~40 MeV, tail 0/386. That fix was correct *symptomatically*
but the AN2021_131_v8 §3.4 review reframes the diagnosis:

- The published mZ measurement runs the CVH refit with
  `useIdealGeometry=True` — the AN calls this out explicitly as
  the nominal simulation reference.
- Under that config, per-hit residuals encode the ~100 μm module
  displacement between design (ideal) geometry and aligned
  geometry (data) or the deliberate misalignment scenario (MC).
- Those residuals are absorbed by a **linearised second stage**:
  per-track Jacobians of the reference-point track state with
  respect to a set of global correction parameters (magnetic
  field, energy loss, alignment) are written by the CVH maker
  (`fillJac=True`), and applied at analysis time from
  `corFiles = cms.vstring(<correctionResults_*.root>)`
  (see AN §3.3, lines 420–433).
- **Our runs to date used `useIdealGeometry=True` with
  `corFiles = cms.vstring()` (empty)** — i.e. we ran the first
  stage of the two-stage AN pipeline without loading the second
  stage. That is not an algorithm bug; it is a workflow-mismatch
  by construction. The σ~80 MeV, 33% tail is exactly what the
  first stage looks like without the second.

Both configurations recover physics-quality tracks:

- **(A) AN-canonical**: `useIdealGeometry=True` +
  `corFiles = ['correctionResults_v721_recjpsidata.root']`. Same
  setup the published mZ measurement used. Scale-consistent with
  the mZ downstream framework. Correction files already exist on
  disk under `wremnants-data/data/calibration/`.
- **(B) Aligned-out-of-box**: `useIdealGeometry=False` +
  `corFiles=[]`. Simpler; what we tested during the debug study;
  physics-quality tracks directly. Decoupled from the mZ scale
  reference.

**Decision (2026-07-02, initial)**: adopt (A) as the btojpsik
default. Rationale: the published measurement uses (A); the
correction file (`correctionResults_v721_recjpsidata.root`) is
already committed to `wremnants-data/data/calibration/`; keeping
our btojpsik pipeline scale-consistent with the mZ analysis lets
the downstream comparison work later without a config-level
re-architecture.

**Decision revised (2026-07-02, mid-implementation)**: adopt (B)
as the btojpsik default; (A) becomes a follow-up. Rationale: the
smoke run against v721 crashed at
`ResidualGlobalCorrectionMakerBase.cc:1017`'s
`assert(nparms == parmset.size())`. Our 15_0_19_patch2 producer
builds a parmset of **95,392** entries while v721's `parmtree`
carries **108,332** — a structural mismatch between the WMass-era
producer that derived v721 and our current build. v721 does carry
an `idxmaptree` (108,332 × `idx/i`) that a properly-written
consumer could use to remap and load only the entries applicable
to the current parmset, but our maker's `corFiles` loader
(lines 1008–1034) does not read `idxmaptree` — it assumes 1:1
alignment. Options considered:

1. Adopt (B); park (A) — **selected**. The debug-study fix
   already validated (B) at σ ≈ 40 MeV, 0% tail on our sample.
   Fastest path to a physics-quality btojpsik pipeline; the
   AN-canonical target becomes a follow-up change-id.
2. Patch the maker's `corFiles` loader to use `idxmaptree` for
   remapping. Medium scope, medium risk — partial loads (missing
   modules → zero correction on those modules) may still not
   reproduce the mZ paper's scale.
3. Re-derive a `correctionResults_*.root` file from our current
   producer. Largest scope; a multi-day project.

The scope of this change is retargeted accordingly.

Two things follow from this:

1. **Development consolidates in `CMSSW_15_0_19_patch2/`.** The
   WMass 10_6_26 fork was only kept alive as a bit-identical reference
   for the debug study. Now that the pathology and its fix are
   understood, all producer edits, drivers, and new features SHALL
   land in `CMSSW_15_0_19_patch2/src/Analysis/HitAnalyzer/` and the
   10_6_26 checkout is treated as read-only reference material.

2. **The alcareco → CVH pipeline is not yet fully validated for
   production.** With the geometry fix in place we can now close out
   the producer-validation deck (`slides/cvh_producer_validation.tex`)
   with the two remaining plots the user asked for: the kaon-side
   ideal-geometry comparison (why did the kaon single-track refit
   look *less* broken than the joint dimuon refit?) and the four-flavor
   B+ mass overlay (raw input tracks, ALCAReco output candidate,
   CVH refit without mass constraint, CVH refit with mass constraint).
   These are validation of the alcareco producer itself — nothing
   downstream yet consumes the refit, so scope stays on producer
   correctness.

Adopting David's fork and candidate-driven driver is **postponed**:
the user will push our current work first, then a diff against
`davidwalter2/cmssw@62becafa038` will drive that follow-up
proposal. Per-file join overhead and stage-3 downstream integration
are explicitly deferred to separate change-ids.

## What Changes

### 1. Consolidate development in `CMSSW_15_0_19_patch2/`

- All new producer edits, drivers, and cfi changes SHALL be made in
  `CMSSW_15_0_19_patch2/src/Analysis/HitAnalyzer/`.
- `CMSSW_10_6_26/` remains checked in as a read-only reference for
  the debug/validation study; no new edits.
- Any producer-level edits (`ResidualGlobalCorrectionMakerBase.cc`,
  `ResidualGlobalCorrectionMakerTwoTrackG4e.cc`,
  `ResidualGlobalCorrectionMakerG4e.cc`) that landed in `10_6_26/`
  during the debug study SHALL be mirrored into `15_0_19_patch2/`.
- No new `README.md` is created for the `10_6_26` checkout; the
  scope note in `slides/cvh_producer_validation.tex` is the single
  place declaring the read-only status.

### 2. Adopt option (B) as btojpsik default: `useIdealGeometry=False` + `corFiles=[]`

btojpsik-scoped drivers SHALL default to
`useIdealGeometry=False` and `corFiles=[]`. This is
option (B) — aligned geometry loaded directly from the CMSSW
GlobalTag, with no linearised Stage-2 correction on top. The
debug-study fix already validated this config: σ ≈ 40 MeV on
[2.95, 3.25], tail 0/386 for the μμ mass on our smoke sample.

Both `useIdealGeometry` and `corFiles` SHALL be exposed as
VarParsing options so operators can override per-job to attempt
option (A) — however, (A) currently crashes at
`ResidualGlobalCorrectionMakerBase.cc:1017` with any WMass-era
`correctionResults_*.root`, and adopting (A) as the physics
default requires the follow-up work described in §9.
The maker cfi defaults themselves stay unchanged (both keep
`useIdealGeometry=True` and `corFiles=cms.vstring()`) — the driver
is the layer that populates `corFiles`, because the choice of
correction file is analysis-workflow-specific.

Concretely:

- `CMSSW_15_0_19_patch2/src/Analysis/HitAnalyzer/test/runCvhBplusJpsiK.py`:
  - Change the `opts.register('useIdealGeometry', True, ...)`
    default from `True` to `False`.
  - Add a `corFiles` VarParsing option (list-of-strings), default
    `[]`. Non-empty values pass through untouched; an operator
    trying option (A) supplies a matching correction file
    explicitly on the CLI.
  - Retain the `_resolve_default_corfile()` helper as a callable
    (not the default value) so a future change activating (A) can
    switch back to it without touching the CLI grammar.
  - Pass the resolved lists into both
    `globalCorJpsiK.clone(corFiles=cms.vstring(*_muon_corfiles), ...)`
    and `globalCorJpsiKKaon.clone(corFiles=cms.vstring(*_kaon_corfiles), ...)`.
  - Keep the `useIdealGeometry=cms.bool(_muon_ideal_geom)` explicit
    forward into `.clone(...)` on both makers (already added).
- For every other btojpsik-scoped launcher enumerated by §7,
  apply the same wiring: `corFiles` VarParsing → explicit forward
  into `.clone(...)`, `useIdealGeometry` explicit forward
  (default `True`).
- `ResidualGlobalCorrectionMakerTwoTrackG4e_cfi.py` and
  `ResidualGlobalCorrectionMakerG4e_cfi.py` cfi defaults remain
  unchanged. Both keep `useIdealGeometry=cms.bool(True)` and
  `corFiles=cms.vstring()`. Non-btojpsik consumers (alignment
  training, other physics channels) that load the cfi without
  overriding stay on the current behaviour.
- Add a one-line `LogInfo("ResidualGlobalCorrectionMakerBase")` at
  maker construction echoing the effective `useIdealGeometry_`
  value **and** the number of `corFiles` entries loaded (e.g.
  `"useIdealGeometry=true, corFiles.size()=1"`). Every future
  job log states its config so operators can audit without
  reading the driver.
- The empty-`corFiles` **with `useIdealGeometry=true`** failure
  mode is now self-diagnosing: any job that logs
  `useIdealGeometry=true, corFiles.size()=0` is running Stage 1
  of the AN two-stage pipeline without Stage 2 and is expected
  to produce σ~80 MeV / 33% tail. With the new (B) defaults this
  configuration is not the accident case (it is only reached by
  an explicit CLI override) so the alarm is louder than the
  silent bug it used to be, but it is still a valid diagnostic
  signal in a job log.

### 3. Kaon-side ideal-geometry × corFiles study

Re-run the B+ chain on the same file used for the dimuon smoke
(`preset_B_final_sl1/Run2016H/jpsix_alcareco_presetB_0029A508-...`)
across three configurations that together separate the effect of
geometry from the effect of the second-stage corrections:

| Point | muon geom | muon corFiles | kaon geom | kaon corFiles | expected |
|---|---|---|---|---|---|
| P1 (AN-canonical) | ideal | v721 | ideal | v721 | physics-quality |
| P2 (aligned baseline) | aligned | empty | aligned | empty | physics-quality |
| P3 (broken baseline) | ideal | empty | ideal | empty | σ~80 MeV / 33% tail |

Also runs the 2×2 muon-geom × kaon-geom scan **within the
AN-canonical column** (four sub-points inside P1 keeping the
correction file loaded) so we can answer:
**why did the kaon leg look less broken than the muon leg under
ideal geometry + empty `corFiles`?** Suspected: fewer hits per
track and larger intrinsic q/p resolution absorb the ~100 μm
alignment displacement without pushing the χ² tail catastrophic;
the effect scales as (alignment displacement / hit resolution) × √N_hits.
The numbers close the loop.

Metrics to report (per point):

- kaon Δp_T/p_T median and MAD (refit vs ALCAReco track);
- kaon Δ(q/p) MAD;
- `chisqval` median and cap-hit fraction, **and `edmval` median +
  fraction meeting `edmval < 1e-5`** (both convergence diagnostics
  reported together);
- B+ candidate mass σ on [5.15, 5.40] and tail fraction
  `|m − 5.279| > 0.20`.

### 3b. Cross-check: (A) vs (B) equivalence

On the same file, produce a **direct overlay** of the μμ mass
and B+ candidate mass under:

- (A) `useIdealGeometry=True` + `corFiles=v721`
- (B) `useIdealGeometry=False` + `corFiles=[]`

The two workflows should produce statistically-indistinguishable
physics-quality distributions (this is exactly the AN's claim: the
corrections make ideal-geometry residuals equivalent to
aligned-geometry residuals for downstream analysis). Any residual
disagreement between (A) and (B) is a real systematic worth
documenting — it would indicate the correction file's alignment-like
term is imperfect for our era / dataset, or that the aligned
geometry loaded by (B) differs from the geometry the v721
correction file was derived against.

Deliverable: one plot with (A) and (B) overlaid for μμ mass and one
for μμK mass, both with their Gaussian σ / mean / tail-fraction
reported.

### 4. Four-flavor B+ mass overlay

Build a validation plot of the B+ candidate mass with four
definitions, on the same event sample, running under the
AN-canonical config (§2: `useIdealGeometry=True` + v721 `corFiles`):

- **raw**: 4-vector sum of the three ALCAReco input tracks (2 μ + K),
  taking the muon/kaon masses from PDG;
- **ALCAReco output**: the persisted B+ candidate mass on
  `ALCARECOTkAlJpsiMuMuResonances` (the `VertexCompositeCandidate`
  from `JpsiXCandidateProducer`);
- **CVH refit, no mass constraint**: 4-vector sum of the CVH-refit
  muons (`Jpsi_*`, i.e. `Muplus_*` + `Muminus_*` after the
  unconstrained joint refit) + CVH-refit kaon (`Kbach_*`);
- **CVH refit, with J/ψ mass constraint**: 4-vector sum of the
  mass-constrained CVH muons (`Jpsicons_*`) + CVH-refit kaon.

Overlay on one axis with distinct colours, fit each with a
Gaussian on [5.15, 5.40], and report σ + mean shift per branch.
Distributions that agree well confirm the alcareco → CVH pipeline
is internally consistent; disagreement pinpoints the layer where
scale/resolution differs.

Note: the AN's linearised-correction application via
`corFiles + fillJac` writes the *Jacobians* to the tree at
CVH-refit time; the actual scalar correction to q/p is applied
downstream in the joined-tree consumer (or in the histmaker).
For this plot we take the CVH-refit branches AS WRITTEN by the
producer under the AN-canonical config. If a corresponding
downstream analysis step is expected to consume the Jacobians,
we note it in the slide caption; verifying that consumer is out
of scope for this change.

### 5. Publication surface

- Extend `slides/cvh_producer_validation.tex` with a **kaon
  ideal-geometry** frame (§3 plots + short table) and a **four-flavor
  B+ mass overlay** frame (§4).
- Dump all figures (PNG) into `~/public_html/mz/cvh/` so David and
  the rest of the group can look at them without cloning the repo.
- Update the pre-production checklist frames with the outcome of
  §6 (GT verification) and §5 (launcher inventory).

### 6. Verify GT consistency

**Clarifying**: the check we run is that the CMSSW GlobalTag used
by the CVH producer job (currently `106X_dataRun2_v37` in
`runCvhBplusJpsiK.py`, `132X_dataRun3_...` in the Run3 branch) matches
the GT that was used to produce the ALCAReco input file at RECO time.
The alignment records loaded by the CVH maker (`TrackerAlignmentRcd`,
`TrackerSurfaceDeformationRcd`) are pulled from that GT; if the two
differ, we could still be propagating against a slightly-wrong
"aligned" geometry even after §2. Concretely:

- read `edmProvDump <input>.root` on one ALCAReco file per era and
  extract the GT string used at RECO;
- diff against the GT set by the CVH driver for that era;
- if they differ, record which alignment IOVs are covered by both
  and pin the driver's GT to the RECO GT (or the closest superset).

Deliverable: a small table in `cvh_producer_validation.tex`
listing (era, ALCAReco GT, CVH driver GT, action taken).

### 7. Inventory all cmsRun launchers touching this producer

Grep for `cmsRun` scripts referencing
`ResidualGlobalCorrectionMakerTwoTrackG4e`,
`ResidualGlobalCorrectionMakerG4e`, `useIdealGeometry`, or
`corFiles` across both CMSSW checkouts and `scripts/`,
`condor/`, `old_condor_stuff/`, `cooper/` (including
`process.load()` paths that inherit the cfi default). For each
entry record (path, mode, effective `useIdealGeometry` value,
effective `corFiles` value, suspected purpose).

Classify each entry:

- **btojpsik-scoped** — part of the J/ψ+X / btojpsik measurement
  pipeline. Update to the AN-canonical config per §2
  (`useIdealGeometry=True`, `corFiles=v721`). Annotation:
  "btojpsik — updated to AN-canonical".
- **not btojpsik-scoped** — belongs to some other workflow
  (alignment training, other physics channel, unclear purpose).
  **Do not modify.** Record the effective config in the inventory
  table with the annotation "flag — do not change" and note the
  suspected workflow. If unclear, mark as "flag — needs owner
  review".

The inventory itself is a mandatory deliverable; the *modifications*
are limited to btojpsik-scoped entries.

### 8. Follow-up: enable AN-canonical (A) as production default

The AN2021_131_v8 §3.3–3.4 canonical config remains the eventual
target for scale consistency with the published mZ measurement.
A separate change-id `enable-an-canonical-corrections` SHALL be
opened as follow-up covering the required work:

- Patch `ResidualGlobalCorrectionMakerBase.cc:1008-1034` to read
  the `idxmaptree` sidecar and remap `parmtree` rows by
  global-index correspondence rather than by row position.
  Partial loads (v721 entries with no corresponding entry in the
  current parmset) SHALL log a per-file summary of skipped
  entries; loaded entries SHALL match the current build's
  `detidparms` map exactly.
- Under the remapped loader, re-run the AN-canonical smoke on
  Run2016H with `useIdealGeometry=True + corFiles=[v721]` and
  compare against the option-(B) result from this change (§3b
  overlay revived). If the two agree within statistical error,
  set the btojpsik driver's `useIdealGeometry` default back to
  `True` and populate `corFiles` with the resolved default.
- If they disagree, escalate: derive a new
  `correctionResults_*.root` from our current producer's parmset
  (the highest-scope option we deliberately deferred here).

Nothing in the present change forecloses (A). Both driver knobs
remain exposed as VarParsing options; the follow-up is a
producer-plugin patch, not a driver rewrite.

### 9. Note on the AN two-stage pipeline (previously "Stage-3 concern")

AN2021_131_v8 §3.3–3.4 documents the two-stage CVH workflow that
motivated the debug study's diagnosis reframe:

1. **Stage 1** — CVH refit with `useIdealGeometry=True`. The
   `fillJac=True` mode writes, for every track, the Jacobian of
   the reference-point track state with respect to a set of
   global correction parameters (magnetic field, energy loss,
   alignment). These are stored on the output tree.
2. **Stage 2** — linearised application of a
   `correctionResults_*.root` file at analysis time (or via the
   maker's own `corFiles` loading path): the stored Jacobian is
   contracted with the correction values to update the track
   state at the reference point. The alignment-like term of the
   correction absorbs the ~100 μm module displacement between
   design and aligned geometry; the other terms absorb residual
   magnetic-field / energy-loss biases.

The debug study observed the σ~80 MeV, 33% tail because we were
running **Stage 1 only** with `corFiles = cms.vstring()` (empty).
That is expected: the alignment residual isn't absorbed by anything
until Stage 2 runs. The AN's published mZ measurement uses (A)
end-to-end.

This proposal originally adopted (A); the pivot to (B) is
documented in the "Decision revised" block above and in §8
follow-up.

Consumers outside the btojpsik pipeline that currently rely on
the empty-`corFiles` cfi default (alignment training, other
channels) are unaffected by this change. §7 inventory flags
them but does not modify.

## Impact

- **Affected specs**: `cvh-refit-jpsi-x` (extended — AN-canonical
  config requirement + `corFiles` wiring + producer-validation
  requirements).
- **Affected code**:
  - `CMSSW_15_0_19_patch2/src/Analysis/HitAnalyzer/test/runCvhBplusJpsiK.py`
    (new `corFiles` VarParsing option, forwarded into both `.clone()` calls)
  - `CMSSW_15_0_19_patch2/src/Analysis/HitAnalyzer/plugins/ResidualGlobalCorrectionMakerBase.cc`
    (add `LogInfo` echoing `useIdealGeometry_` and `corFiles_.size()`)
  - `scripts/btojpsik/plot_bplus_mass_overlay.py` (new)
  - `scripts/btojpsik/plot_kaon_ideal_geom_compare.py` (new — now
    also carries the (A)-vs-(B) cross-check from §3b)
  - `scripts/btojpsik/inventory_cvh_launchers.sh` (new)
  - `slides/cvh_producer_validation.tex` (extended)
  - **Explicitly NOT edited**:
    `CMSSW_15_0_19_patch2/src/Analysis/HitAnalyzer/python/ResidualGlobalCorrectionMakerTwoTrackG4e_cfi.py`,
    `CMSSW_15_0_19_patch2/src/Analysis/HitAnalyzer/python/ResidualGlobalCorrectionMakerG4e_cfi.py`
    (cfi defaults preserved).
- **Affected inputs**:
  - `wremnants-data/data/calibration/correctionResults_v721_recjpsidata.root`
    becomes a load-time dependency of every btojpsik CVH job.
- **Affected outputs**:
  - New CVH outputs under `/ceph/.../cvh_outputs/kaon_ideal_geom_compare/`
    covering configurations P1/P2/P3 (§3) and the (A)-vs-(B)
    cross-check (§3b).
  - Figures under `~/public_html/mz/cvh/`.

## Out of scope

- Adopting David's fork (`davidwalter2/cmssw@62becafa038`) or the
  candidate-driven driver. Postponed pending user's push +
  diff-and-decide follow-up.
- Per-file join overhead in
  `scripts/btojpsik/join_cvh_bplus_jpsik.py`. Will be addressed in
  a separate proposal.
- Any Stage-3 downstream integration (btojpsik histmaker inputs,
  tensor updates) — the alcareco → CVH pipeline is being validated
  in isolation.
- WMass `CMSSW_10_6_26` code changes. That checkout is now
  read-only reference.
- Reproducing the σ 80 → 40 MeV smoke on the full 100k sample.
  Once the checklist is green we launch full-scale AlCaReco in
  a separate production-management change-id.

## Open questions

1. **RESOLVED (2026-07-02)**: adopt option (A) — AN-canonical
   config. btojpsik drivers set `useIdealGeometry=True` +
   `corFiles=[correctionResults_v721_recjpsidata.root]`. cfi
   defaults untouched.
2. For §6, the ALCAReco GT is technically recorded per-lumi in
   the ProcessHistory. Do we need the strict per-lumi GT, or is
   the driver's per-era GT close enough for the alignment records
   we care about? First pass: per-era check; escalate only if
   §6 uncovers a mismatch.
3. §4 four-flavor overlay — should the "raw" flavor use ALCAReco
   input-track four-vectors *before* any mass constraint (i.e. the
   pre-refit `RefTrack`s), or the ALCAReco output candidate
   daughters (which may already have the J/ψ mass constraint
   embedded from `JpsiXCandidateProducer` under preset C)? Under
   preset B — our operational default — the constraint is off, so
   this collapses to a single definition; call it out in the slide
   caption.
4. **NEW — is `correctionResults_v721_recjpsidata.root` accurate
   for our Run2016H sample?** The file was chosen provisionally
   (v721 is one of several vintages committed to `wremnants-data`
   alongside `v718`, `lbl2017`, `lbl2018`); v721 does NOT need to
   be the "ideal" file, only *accurate*. Confirmation happens in
   two places: (i) task 2.9 traces provenance in the git log /
   any README, (ii) §3b (A)-vs-(B) direct overlay is the empirical
   check — if they agree within statistical error on Run2016H,
   v721 is confirmed accurate for our sample. Disagreement raises
   this as a follow-up (picking `lbl2016...` if it exists, or
   requesting a derivation) but does NOT block this proposal.
5. **NEW — does downstream analysis apply the reference-point
   Jacobian correction?** The CVH producer writes the Jacobian
   with `fillJac=True`; the actual q/p update from
   `corFiles` values happens either inside the producer (via the
   `corFiles` load) or downstream in the analysis-time consumer.
   The maker's behaviour under `corFiles=[file]` needs a quick
   read: does the producer apply the correction to the output
   branches, or does it only write Jacobians for downstream
   consumption? If the latter, our validation plots need a
   corresponding downstream step for the (A) column to match
   (B). This is a code-reading task, not blocking — will be
   answered by §3's smoke run and the LogInfo dump.
