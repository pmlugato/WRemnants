## ADDED Requirements

### Requirement: CVH producer development home is `CMSSW_15_0_19_patch2`

New CVH producer code SHALL be developed in `CMSSW_15_0_19_patch2/src/Analysis/HitAnalyzer/`, covering plugin sources, cfi defaults, driver scripts, and cross-driver validation harnesses. The `CMSSW_10_6_26/src/Analysis/HitAnalyzer/` checkout SHALL be treated as read-only reference for the debug/validation study and SHALL NOT receive new edits.

#### Scenario: Producer edit lands in the target checkout
- **WHEN** a change is proposed to `ResidualGlobalCorrectionMakerTwoTrackG4e.cc` or `ResidualGlobalCorrectionMakerG4e.cc`
- **THEN** the edit SHALL target the file under `CMSSW_15_0_19_patch2/src/Analysis/HitAnalyzer/plugins/`
- **AND** the corresponding path under `CMSSW_10_6_26/` SHALL remain unchanged

#### Scenario: Developer looks up the current source
- **WHEN** a developer greps for the current CVH producer implementation
- **THEN** `CMSSW_15_0_19_patch2/src/Analysis/HitAnalyzer/README.md` (or the equivalent scope note in `slides/cvh_producer_validation.tex`) SHALL identify `CMSSW_15_0_19_patch2` as the development branch and `CMSSW_10_6_26` as reference-only

### Requirement: btojpsik drivers default to option (B) — `useIdealGeometry=False` + `corFiles=[]`

Every btojpsik-scoped CVH driver SHALL default to
`useIdealGeometry=False` and `corFiles=[]` (option (B), aligned
geometry loaded from the CMSSW GlobalTag with no Stage-2
correction). Both values SHALL remain exposed as VarParsing options
so operators can attempt option (A) per-job (see follow-up
`enable-an-canonical-corrections` for the plugin work required to
make (A) production-viable). The maker cfi files
(`ResidualGlobalCorrectionMakerTwoTrackG4e_cfi.py`,
`ResidualGlobalCorrectionMakerG4e_cfi.py`) SHALL retain their
existing `useIdealGeometry=cms.bool(True)` and `corFiles=cms.vstring()`
defaults so non-btojpsik consumers (alignment training, other
channels) are not silently modified.

#### Scenario: Default btojpsik driver invocation uses aligned geometry
- **WHEN** `cmsRun runCvhBplusJpsiK.py mode=dimuon input=<file>` is run with no `useIdealGeometry` or `corFiles` override
- **THEN** the producer pset for both `globalCorJpsiK` and `globalCorJpsiKKaon` SHALL contain `useIdealGeometry = cms.bool(False)` and `corFiles = cms.vstring()`
- **AND** the smoke-run outputs (σ_wide on [2.95, 3.25], tail fraction `|m − 3.097| > 0.20`) SHALL match the previously-recorded (B) numbers (σ ≈ 40 MeV, tail ≈ 0%) within statistical error
- **AND** the constructor-time `LogInfo` SHALL record `useIdealGeometry=false, corFiles.size()=0`

#### Scenario: cfi defaults are preserved
- **WHEN** a hypothetical driver `process.load('Analysis.HitAnalyzer.ResidualGlobalCorrectionMakerTwoTrackG4e_cfi')` and adds no `.clone(useIdealGeometry=..., corFiles=...)` override
- **THEN** the resulting maker instance SHALL see `useIdealGeometry = True` and `corFiles.size() = 0` (both inherited from the unchanged cfi defaults)
- **AND** the constructor-time `LogInfo` SHALL record `useIdealGeometry=true, corFiles.size()=0`

#### Scenario: Explicit override to attempt (A) surfaces the parmset-mismatch failure clearly
- **WHEN** a driver sets `useIdealGeometry=True` and `corFiles=['/path/to/correctionResults_v721_recjpsidata.root']` explicitly on the CLI
- **THEN** the producer SHALL crash at `ResidualGlobalCorrectionMakerBase.cc:1017`'s `assert(nparms == parmset.size())` with `nparms=108332, parmset.size()=95392` (or equivalent) — this failure IS the diagnostic pointer to the follow-up `enable-an-canonical-corrections` change-id
- **AND** the driver SHALL NOT paper over this failure with a partial load

#### Scenario: btojpsik driver forwards resolved values into every clone
- **WHEN** a btojpsik driver instantiates more than one maker (e.g. `globalCorJpsiK` and `globalCorJpsiKKaon`)
- **THEN** each `.clone(...)` call SHALL pass both `useIdealGeometry` and `corFiles` explicitly, so no maker in the driver ever falls through to the cfi defaults silently

### Requirement: Effective `useIdealGeometry` and `corFiles` values are logged at construction

`ResidualGlobalCorrectionMakerBase` SHALL emit exactly one `LogInfo("ResidualGlobalCorrectionMakerBase")` message at constructor time reporting both the effective `useIdealGeometry_` value and the number of `corFiles_` entries (e.g. `"useIdealGeometry=true, corFiles.size()=1"`).

#### Scenario: Log line appears in job output
- **WHEN** a cmsRun job instantiates any subclass of `ResidualGlobalCorrectionMakerBase`
- **THEN** the CMSSW log SHALL contain a line matching `useIdealGeometry=(true|false), corFiles.size()=[0-9]+` under the `ResidualGlobalCorrectionMakerBase` category

#### Scenario: One line per maker instance, not per event
- **WHEN** a cmsRun job processes N events with two maker instances (`globalCorJpsiK` and `globalCorJpsiKKaon`)
- **THEN** the `LogInfo` SHALL appear exactly twice per job (once per maker instance), not once per event

#### Scenario: Broken-config alarm is self-diagnosing
- **WHEN** a cmsRun job logs `useIdealGeometry=true, corFiles.size()=0`
- **THEN** operators reading the job log SHALL recognise this as Stage-1-only (broken) configuration, matching the σ~80 MeV / 33% tail signature documented in the debug study

### Requirement: Per-leg `useIdealGeometry` and `corFiles` knobs in the B+ chain driver

The B+ chain driver SHALL expose `useIdealGeometryMuon`, `useIdealGeometryKaon`, `corFilesMuon`, and `corFilesKaon` as VarParsing options in `CMSSW_15_0_19_patch2/src/Analysis/HitAnalyzer/test/runCvhBplusJpsiK.py`, in addition to the base `useIdealGeometry` and `corFiles` options. If any per-leg option is set, it SHALL take precedence for that leg; if unset, each leg SHALL inherit the corresponding base value.

#### Scenario: Independent per-leg control of geometry
- **WHEN** `cmsRun runCvhBplusJpsiK.py mode=dimuon input=<file> useIdealGeometryMuon=False useIdealGeometryKaon=True` is run
- **THEN** `globalCorJpsiK.useIdealGeometry` is `False` and `globalCorJpsiKKaon.useIdealGeometry` is `True`

#### Scenario: Independent per-leg control of `corFiles`
- **WHEN** the driver is invoked with `corFilesMuon=[foo.root]` and `corFilesKaon=[]`
- **THEN** `globalCorJpsiK.corFiles` is `cms.vstring('foo.root')` and `globalCorJpsiKKaon.corFiles` is `cms.vstring()`

#### Scenario: Base value inheritance
- **WHEN** only `useIdealGeometry=True` and `corFiles=[v721_path]` are set and no per-leg option is set
- **THEN** both maker instances see `useIdealGeometry = True` and `corFiles` populated with `[v721_path]`

### Requirement: Alcareco → CVH pipeline validation deliverables

Before full-scale AlCaReco production kicks off, `slides/cvh_producer_validation.tex` SHALL contain the following frames documenting alcareco → CVH producer correctness:

1. Prompt J/ψ→μμ dimuon mass comparison (raw / ALCAReco / CVH refit)
   under the previously-tested `useIdealGeometry=False` config;
2. B+ → J/ψK dimuon mass comparison under the previously-tested
   `useIdealGeometry=False` config;
3. AN two-stage pipeline explainer frame describing Stage 1
   (CVH refit with `useIdealGeometry=True` writing per-track
   Jacobians via `fillJac=True`) and Stage 2 (linearised correction
   applied via `corFiles = ['correctionResults_v721_recjpsidata.root']`),
   with a citation to AN2021_131_v8 §3.3–3.4;
4. P1 / P2 / P3 comparison frame: (P1) AN-canonical
   `useIdealGeometry=True + v721`, (P2) aligned baseline
   `useIdealGeometry=False + []`, (P3) broken baseline
   `useIdealGeometry=True + []`; μμ mass and B+ candidate mass
   overlays with fit table;
5. (A) vs (B) direct-overlay frame: (A) AN-canonical vs
   (B) aligned baseline on the same event sample, with per-fit
   σ / mean / tail-fraction comparison;
6. Kaon-side geometry sub-scan within P1 (muon geom × kaon geom
   with `corFiles` loaded on both legs) plus the explanation of
   why the kaon leg looked less broken than the muon leg under P3;
7. Four-flavor B+ candidate mass overlay (raw input tracks,
   ALCAReco output candidate, CVH refit without J/ψ mass
   constraint, CVH refit with J/ψ mass constraint) run under the
   AN-canonical config, with per-flavor Gaussian fit parameters;
8. Inventory table of all cmsRun launchers touching
   `ResidualGlobalCorrectionMakerTwoTrackG4e` /
   `ResidualGlobalCorrectionMakerG4e`, listing the effective
   `useIdealGeometry` and `corFiles` values per launcher;
9. GT verification table (era, ALCAReco RECO GT, CVH driver GT,
   action taken).

Frames 1–2 already exist. Frames 3–9 SHALL be added by this
change before the deck is considered production-ready.

#### Scenario: Deck rebuilt after change lands
- **WHEN** `slides/cvh_producer_validation.tex` is compiled
- **THEN** the resulting PDF SHALL contain frames titled (or clearly showing) each of the nine items above
- **AND** no LaTeX `Overfull hbox`/`vbox` warning SHALL cause visible content overflow past the printable area

#### Scenario: Figures published for group access
- **WHEN** a frame includes a figure from `slides/figs/cvh_producer_validation/`
- **THEN** the same PNG SHALL also be present under `~/public_html/mz/cvh/` (accessible via `https://submit.mit.edu/~pmlugato/mz/cvh/`)

### Requirement: GT consistency check between CVH driver and ALCAReco input

The CMSSW GlobalTag string set in the CVH driver SHALL be verified against the GT recorded in the corresponding ALCAReco input file's `ProcessHistory` (at RECO time), for every era in use (2016H, 2017, 2018, plus Run3 as applicable). Mismatches SHALL be resolved before the corresponding era's data enters full-scale AlCaReco production.

#### Scenario: Match — no action
- **WHEN** the CVH driver GT for era E equals the ALCAReco RECO GT for era E
- **THEN** the corresponding row in the GT verification table SHALL read "match — no action"

#### Scenario: Mismatch — driver GT pinned
- **WHEN** the CVH driver GT for era E differs from the ALCAReco RECO GT for era E
- **THEN** the driver GT SHALL be updated to the ALCAReco RECO GT (or a GT that is a strict superset with respect to the alignment records), and the GT verification table SHALL record both values plus the action taken

### Requirement: Launcher inventory captured before btojpsik driver wiring lands

An inventory SHALL be produced before any btojpsik-scoped driver's AN-canonical wiring lands, enumerating every cmsRun launcher and cfi consumer referencing `ResidualGlobalCorrectionMakerTwoTrackG4e`, `ResidualGlobalCorrectionMakerG4e`, `useIdealGeometry`, or `corFiles` under `CMSSW_15_0_19_patch2/`, `CMSSW_10_6_26/`, and top-level `scripts/`, `condor/`, `old_condor_stuff/`, and `cooper/`. Each row SHALL record the file path, the mode (dimuon / single-track / both), the effective `useIdealGeometry` value (explicit `True`, explicit `False`, or relies-on-default), the effective `corFiles` value (explicit list / empty / default), and a classification (btojpsik-scoped or not-btojpsik-scoped).

#### Scenario: Every launcher classified
- **WHEN** `scripts/btojpsik/inventory_cvh_launchers.sh` is run
- **THEN** its output SHALL enumerate every consumer of the CVH maker under the search roots
- **AND** each row SHALL identify the current `useIdealGeometry` and `corFiles` state (explicit or default)

#### Scenario: Non-btojpsik consumers flagged, not modified
- **WHEN** the inventory finds a launcher outside the btojpsik pipeline (alignment training, other physics channel, or unclear purpose)
- **THEN** the entry SHALL be flagged in the deck's inventory table with the annotation "flag — do not change" and its suspected workflow
- **AND** no edit SHALL be made to that launcher in this change

#### Scenario: btojpsik consumers updated to AN-canonical config
- **WHEN** the inventory finds a btojpsik-scoped launcher (path under `scripts/btojpsik/`, `condor/`, `old_condor_stuff/`, `cooper/`, or a CVH driver known to feed the btojpsik pipeline)
- **THEN** the entry SHALL be annotated "btojpsik — updated to AN-canonical" and the corresponding driver SHALL be updated to wire `useIdealGeometry=True` and `corFiles=[correctionResults_v721_recjpsidata.root]` (or their per-leg equivalents) with explicit `.clone(...)` forwarding
