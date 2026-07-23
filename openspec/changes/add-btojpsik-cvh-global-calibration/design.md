## Context

The CVH global-correction calibration has three stages, all understood from
the existing code:

1. **CVH refit (grads production).** A `ResidualGlobalCorrectionMaker*` EDProducer
   re-fits each track with Geant4e through the real field and material, and —
   when `fillGrads=True` — writes a per-candidate sidecar tree with:
   - `nParms`, `globalidxv[nParms]` — the global parameter indices this
     candidate touches,
   - `gradv[nParms]` — the gradient of the candidate χ² w.r.t. those parameters,
   - `hessfactorv` (row-major `nRank×nParms`, H = BᵀB) when
     `fillGradsFactored=True`, else `hesspackedv[nSym]` (packed upper triangle),
   - quality/diagnostic branches (`edmvalref`, `chisqval`, `ndof`, masses, …).
   The parameter catalog itself is written once as `runtree` (one row per
   global parameter: `parmtype`, `subdet`, `layer`, `eta`, `phi`, `rawdetid`).
   Ref: `ResidualGlobalCorrectionMakerBase.cc` branch setup (`nParms`,
   `globalidxv`, `gradv`, `hessfactorv`, `hesspackedv`, `run/lumi/event`).

2. **Aggregation** (`aggregategrads.py` + `aggregategrads.cpp`). TChain all
   `globalcor_*.root`, apply quality filters, and scatter-accumulate each
   candidate's `gradv`/`hess*` into a dense global `grad[nparms]` and
   `hess[nparms,nparms]` via `GradHelper` / `HessHelperFactored`, weighted by a
   per-candidate `massweightval` (currently `1.0`). Output: `combinedgrads.hdf5`.

3. **Solve** (`globalfith5pynoreduction.py`). Load one or more
   `combinedgrads.hdf5` (each with a scalar weight), sum them, read the
   `runtree` catalog, add per-parmtype Gaussian priors and NMR absolute-|B|
   constraint rows on the field-mode block, Cholesky-factor and solve
   `H x = −g`, optionally invert `H` for the covariance, and write
   `correctionResults_*.root` (`parmtree` with `x`, `err`; `idxmaptree`).

The key realization behind this proposal: **combining J/ψ and kaon is additive
accumulation in a shared parameter space, not a new approximation.** Each maker
already linearizes (Gauss–Newton) about the current parameters and emits its
grad/Hess in the same global index space. Independent tracks contribute
independent, additive terms to the total least-squares normal equations:
`g_total = g_Jψμμ + g_kaon`, `H_total = H_Jψμμ + H_kaon`. The solve stage
already sums a list of aggregated files (`fgradnames`/`weights`), so the
mechanism exists; the work is producing the kaon contribution correctly and
enabling the reporting.

Parmtype legend (from `globalfith5pynoreduction.py` prior widths):
0/1 local translations, 2–5 alignment rotations/translations,
6 = per-detid b-field (σ=0.038, **legacy — not emitted by this build**),
7 = per-module material energy loss (σ=ln2),
8/9 hit resolution, 10 material resolution, 11 ionization resolution,
**14 = B-field scalar-potential modes** (~50, global), **15 = per-group
material scale** (priors from `materialGroups50.txt`).

**Parameterization of this build (verified against the grads50 `runtree`,
95,082 params, and `ResidualGlobalCorrectionMakerBase.cc`):** the b-field is
fit *only* as the ~50 global scalar-potential modes (parmtype 14) — there is no
per-detid b-field (parmtype 6). Material is fit as *either* per-module eloss
(parmtype 7, `globalMaterialModel=False`, what grads50 used) *or* ~50 grouped
scales (parmtype 15, `globalMaterialModel=True` + `materialGroupsFile`), never
both — this is a production-time config choice. **This study uses grouped
material (parmtype 15)** on all makers.

## Goals / Non-Goals

**Goals**
- **Primary:** over the same B⁺→J/ψK⁺ events, produce a dimuon-only (jpsicons)
  fit and a dimuon+kaon fit, and quantify what the kaon adds — per-parameter
  constraint σ split into field **scale** (mode 0), field **shape** (modes 1+)
  and **material** (parmtype 15), plus magnitudes and the field↔material
  correlation as a diagnostic.
- **Species closure check:** rerun the same bachelor tracks under the *muon*
  mass hypothesis (`kaonAsMuon=True`) to isolate the dE/dx/species effect on the
  material block at fixed track sample.
- **Sanity check:** confirm the standalone J/ψ→μμ (`TkAlJpsiMuMu`) calibration
  and the J/ψ-from-B dimuon-only calibration give consistent b-field/material
  corrections, so the primary comparison rests on a sound J/ψ baseline.

**Non-Goals**
- **Not** an attempt to set the absolute momentum scale — that is the J/ψ's
  role and is not what this channel is being asked to do.
- **No K⁰s channel.** Deliberately excluded so the B⁺→J/ψK⁺ contribution is seen
  in isolation, and the fit is run independently rather than joined to the
  existing combined fit.
- Not a claim that this channel resolves the calibration; it is an exploratory
  measurement of how it behaves.
- No production correction-file release; no application of the derived
  corrections downstream (no rabbit/mZ propagation).
- No full-dataset (F/G/H) run — Run2016H (or a subset) only.
- No "full inclusive J/ψ→μμ + B kaon" combined fit — the primary comparison is
  within the B sample (same events), which isolates the kaon effect better.
- No prompt/non-prompt J/ψ selection — B candidates carry only DOCA cuts (not a
  displacement cut), so the split is not clean, and it is irrelevant to the
  mass-constrained momentum-scale corrections anyway.
- No change to the CVH refit algorithm, geometry config, or the field/material
  parameterization itself.

## Decisions

### Decision 1: The comparison is within the B sample (same events)

The primary measurement uses the **same B⁺→J/ψK⁺ events** on both sides:
dimuon-only (jpsicons) vs dimuon+kaon. Within the combined fit, the two muon
tracks and the kaon track are independent measurements (different hits), so
their grad/Hess add cleanly; the shared global correction parameters are
exactly what they co-constrain. The dimuons appear once in each fit — this is
not double-counting, it is the controlled variable held fixed while the kaon is
switched on.

- **Why within-B rather than "full inclusive J/ψ→μμ + B kaon":** if the kaon
  (~80k B⁺ candidates) is added on top of the full inclusive J/ψ→μμ sample
  (millions of candidates), the b-field/material blocks are already so tightly
  pinned that the kaon barely moves the covariance — the degeneracy-breaking
  signal is diluted. Within-B makes the kaon a comparable-statistics addition,
  so the correlation change is large and visible.
- **Alternatives considered:**
  - Full inclusive J/ψ→μμ + B kaon → diluted sensitivity; kept only as an
    optional realistic-setting cross-check, out of scope for this change.
  - Sum inclusive J/ψ→μμ dimuon + B dimuon → double-counts overlapping muons
    (the B dimuons are a subset of the inclusive sample); rejected.
  - Prompt/non-prompt selection to make samples disjoint → rejected (not a
    clean split under DOCA-only cuts; irrelevant to the momentum-scale physics).

### Decision 1b: J/ψ-consistency sanity check

To confirm the primary comparison rests on a sound J/ψ baseline, separately
solve the nominal standalone J/ψ→μμ calibration on `TkAlJpsiMuMu` and the
J/ψ-from-B dimuon-only calibration, and check that their b-field/material
corrections agree within uncertainties. Disagreement would mean the B dimuons
are not behaving like the standard J/ψ sample (e.g. a selection or reco
difference), which must be understood before trusting the within-B result.
Prompt vs non-prompt is not expected to matter here: the mass-constrained J/ψ
fit pulls m(μμ)→3.097 identically, and the per-track hit residuals that feed
the corrections are independent of the production vertex.

### Decision 1c: Grouped material model (parmtype 15)

Material is fit as ~50 grouped scales (parmtype 15, `globalMaterialModel=True`
+ `materialGroupsFile` + `skipHitlessSurfaces=True`) rather than the per-module
form (parmtype 7). This gives an interpretable per-group material magnitude and
a compact ~50×~50 b-field↔material (parmtype 14 × 15) correlation block, which
is exactly what the study reports. The b-field is always the ~50 scalar-
potential modes (parmtype 14); parmtype 6 does not exist in this build.

- **Consequence:** the standalone J/ψ→μμ sanity-check calibration must also be
  run with grouped material so its catalog matches (Decision 2) and the
  comparison is apples-to-apples. This is automatic — we generate that
  calibration ourselves; we do NOT compare to David's per-module grads50.
- **Alternative considered:** per-module material (parmtype 7) matches grads50
  but yields ~13k material params with no single interpretable "magnitude";
  rejected for this study. Neither form co-exists with the other in one fit, so
  reporting both would require a second production — out of scope.

### Decision 2: Parameter-catalog consistency is a hard precondition

Summed grad/Hess are only meaningful if `globalidxv` means the same thing for
every maker. The dimuon and kaon makers, and the baseline `TkAlJpsiMuMu`
maker, MUST be built in the same CMSSW area with the same alignment/field/
material parameterization (same `materialGroups50.txt`, same field-mode set).
- **Decision:** add an explicit validation step that asserts the `runtree`
  catalogs (`GetEntries()` and the ordered `(parmtype, subdet, layer,
  rawdetid)` tuples) are identical across the baseline, dimuon, and kaon
  productions before any aggregation is trusted. A mismatch is a hard stop.

### Decision 3: Covariance / correlation extraction

`docov=False` today returns zero errors. Full dense inversion of `H`
(`nparms`~10⁵) is expensive. For a study focused on the b-field/material
blocks:
- **Decision:** enable covariance but extract only what is needed — the
  diagonal errors for all reported parameters and the **block sub-covariance**
  coupling the b-field field-modes (parmtype 14) to the material groups
  (parmtype 15), from which the correlation matrix is formed. Prefer the
  targeted triangular-inverse path (`dtrtri` on the Cholesky factor) already
  sketched in the solve script over the full identity solve when memory is
  tight.

### Decision 4: Fork location and repo footprint

- **Decision:** fork `MuonMomentumScaleCalibration` to
  `/work/submit/pmlugato/REALmz/mz/real/calibration/`. The heavy solve code lives in the fork
  (not versioned here). Thin, repo-versioned wrappers and the reporting/plot
  scripts live under `scripts/btojpsik/calibration/` so the analysis entry
  points and comparison outputs are tracked in this repo.

## Risks / Trade-offs

- **Catalog drift between makers** → summed indices misalign silently.
  Mitigation: Decision 2 validation gate.
- **Kaon single-track convergence / soft-kaon losses** — recall `plimit`
  matters for soft bachelors; low-momentum kaons may fail propagation.
  Mitigation: run the kaon side with `plimit=0.05` (as in prior studies) and
  report the surviving fraction; the kaon quality filter must mirror the
  dimuon one in spirit (nvalid, pixel hits, χ²/ndof, edm).
- **Low B-dimuon statistics** → the dimuon-only fit is prior-dominated on the
  alignment blocks (~80k candidates, limited kinematics). Acceptable: the
  b-field/material blocks (global 14/15; per-detid 6/7 along traversed modules)
  are still constrained, and the reported correlation is set by the kinematic
  lever arm, not raw statistics. Report which blocks are prior-dominated.
- **NMR / prior consistency** — the dimuon-only and dimuon+kaon solves must
  apply the identical prior + NMR-constraint configuration, else the comparison
  conflates the kaon effect with a prior change. Mitigation: single shared
  solve config; vary only the input grad list.
- **`useScalarPot3D=False` would silently remove the field block.** The existing
  diagnostic launcher `scripts/btojpsik/run_full_run2016H.sh` passes
  `useScalarPot3D=False` (stock CMSSW field) — correct for the mass-shape
  studies it was written for, but it must NOT be copied into the grads
  production: the b-field parameters this study reports are the scalar-potential
  modes. Mitigation: the calibration launcher keeps the driver default
  (`True`) and says so in its header; the smoke test asserts parmtype 14 is
  present in the `runtree`.
- **`fillGradsFactored` is an untracked parameter.** Passing it as a tracked
  `cms.bool` lands in a different ParameterSet namespace and silently falls back
  to `False` — the packed Hessian instead of the factored one, with no error.
  Mitigation: the driver passes `cms.untracked.bool` with a comment.
- **Compute** — full Run2016H CVH grads is many `cmsRun` jobs (two per file for
  B: dimuon + kaon). Mitigation: subset allowed; condor batch reuse from prior
  btojpsik productions.

## Migration Plan

Additive study; no migration. The maker cfi changes default the new grad flags
OFF (bit-identical to current diagnostic-join behavior) and are turned on only
by the calibration driver invocation.

### Decision 5: Hold alignment fixed (measured, 2026-07-22)

99.9% of the parameters are alignment (81,732 of 81,824). A diagnostic ran the
*same* B-dimuon grads (4,881 candidates) through two solves — alignment
marginalised vs alignment fixed (`CVH_FIX_ALIGNMENT=1`):

| block | mean err (marg) | mean err (fixed) | ratio |
|---|---|---|---|
| b-field modes (14) | 1.975 | 1.304 | **0.48** (median 0.43) |
| material groups (15) | 4.303e-2 | 4.291e-2 | 1.00 |

| | max abs rho (field↔material) |
|---|---|
| alignment marginalised | 0.134 |
| alignment fixed | **0.315** |

Findings:
- **b-field uncertainties are ~2× inflated by alignment marginalisation** — over
  half of the b-field error is alignment leakage, which no amount of kaon data
  can remove. Measuring the kaon's effect through that is diluted.
- **Material is unaffected by alignment but is 100% prior-dominated** at this
  statistics: the fitted error equals the prior width exactly (median 0.0200).
  Dimuons carry essentially no energy-loss information — which is precisely the
  gap the kaon is meant to fill, so this is the expected "before" picture.
- **Alignment marginalisation *suppresses* the field↔material correlation**
  (0.13 vs 0.31) — it masks the very observable this study targets.

**Decision:** run with **alignment fixed** (`CVH_FIX_ALIGNMENT=1`) as the
nominal, so only the 92 global parameters float. This matches the reference fit,
which does not float alignment. No separate alignment stage is needed: with no
prior solution to import, alignment is held at the GT geometry (zero residual),
which is the same implicit choice the reference fit makes.

- Consequence: the reported field/material uncertainties are *conditional* on
  the GT alignment. That is the intended comparison here — the same condition
  applies to both sides of dimuon-only vs dimuon+kaon, so the delta is the
  kaon's. The marginalised numbers are also produced, as the honest figure for
  any absolute statement.
- If a future iteration wants alignment-marginalised absolute errors, a Stage-1
  alignment determination on the central ALCARECO
  (`/ceph/submit/data/group/cms/store/data/Run2016{F,G,H}/Charmonium/ALCARECO/`,
  ~46M events for Run2016H vs ~3.6M in our partial re-reco) is the natural
  source; CVH fits *residual* corrections relative to the GT geometry, so it
  transfers provided the GT and `useIdealGeometry=False` match. Out of scope for
  this change.

### Decision 6: Selection retuned from measured cut-flows (2026-07-22)

The inherited thresholds came from the prompt-J/ψ calibration:
- **Muon pT 5 → 3 GeV.** The 5 GeV cut alone kept 24% of B dimuons while every
  other cut kept ~100%, and dimuon χ²/ndof is healthy (median 1.10). Lowering to
  3.0 gives **2.28×** (survival 24.4% → 79.6%); below that returns diminish.
- **Kaon pT 0.5 → 0.2 GeV.** Gains +22%, and those tracks pass χ²<3 — 223 good
  kaons below 0.5 GeV, exactly the low-momentum ones carrying energy-loss
  information.
- **Kaon χ²/ndof stays at 3.** Initially assumed too tight for hadrons; the data
  refutes this — the kaon χ²/ndof distribution has median 938 and 75th percentile
  18,715 (short tracks, few dof, hadronic interactions, decay in flight). The cut
  removes genuinely broken fits; loosening 3 → 10 recovers only +3%.

All are env-overridable so the original values run as a systematic.

## STATUS: PAUSED pending `add-alcareco-to-nanoaod` (2026-07-22)

That change generalises how the CVH refit is run, and two of its decisions
directly rewrite parts of this one. Production is stopped; this change should be
**refactored on top of it** rather than finished as-is.

**What it changes for us**
- **Q8: `JpsiKCandidateSplitter` is deleted** — the makers descend the nested VCC
  themselves. Our `CVH_SRCCANDS` pin (`jpsiKCandidateSplitter:dimuon`) becomes
  obsolete, and every dimuon grads file produced against the splitter must be
  regenerated against the nested-candidate path.
- **One job, both makers** (shared G4 master, `e33c8d`). The dimuon and kaon
  grads then come from the *same job on the same events* **by construction** —
  which deletes the entire file-matching apparatus (`match_fileids.sh`,
  `CVH_GRADS_FILEIDS`, the matched-set fingerprinting in `run_study.sh`) and with
  it the class of bug that silently mismatched a 17-file dimuon aggregation with
  a 29-file kaon one.
- **Generalisation to all channels** is a *physics* win here, not just plumbing.
  Today this study has one soft bachelor (K⁺ from B⁺). The generalised path
  exposes bachelors from every persisted channel — K/π from K*⁰, K/K from φ,
  π/π from K⁰_S, p/π from Λ_b, π/π from ψ(2S). That is a much denser population
  of **soft displaced hadrons**, which is precisely the lever the reference fit
  identifies for field *shape* and *material*. A single 80k-candidate kaon sample
  was always going to be marginal; the multi-channel version is the version of
  this study that can actually move the blocks it targets.

**What survives the refactor unchanged** (hard-won; do not re-derive)
- Catalog gate (`check_parm_catalog.py`) — identical `runtree` across every
  production is still the precondition for summing.
- 50-mode scalar-potential dump pinned (NMR constraints are 4×50; the 360-mode
  default fails the solve's assert).
- Grouped material (parmtype 15) + field modes (parmtype 14) = the 92-parameter
  model; alignment **not** floated (Decisions 5, 6).
- Targeted covariance extraction — upstream had `docov=False` and an
  `assert(0)`; the ~100-column solve is what makes uncertainties/correlations
  available at all.
- Reporting split: field mode 0 (scale) vs modes 1+ (shape) vs material.
- Selection retuned from measured cut-flows (muon pT>3, kaon pT>0.2, χ²/ndof<3).
- Worker robustness: `fillGradsFactored` must be `cms.untracked.bool`; stage the
  workdir on cephfs (`${OUT_DIR}/.scratch`), never node `/tmp`; and **promote an
  output only if cmsRun returned 0 AND the ROOT file opens with a valid
  `tree`+`runtree`**.

**Production state at pause** (`/ceph/.../mz/cvh_calib/run2016H`, 672 GB)

| channel | good | corrupt | note |
|---|---|---|---|
| bkaon | 718 | 0 | ~379 MB/file |
| bkaon_asmuon | 760 | 5 | ~379 MB/file |
| bdimuon | 43 | 29 (40%) | ~800 MB/file |
| jpsimumu | (survey incomplete) | | |

Corruption scales with output size: large-output jobs were killed during
teardown, leaving fully sized but never-closed ROOT files that the pre-fix
worker promoted with an `OK`. The worker now rejects them. All of this data is
pre-refactor and will be regenerated.

## Context from the reference fit (David Walter, 2026-07-21)

`260721_magnet_field_david.pdf` fits the **same 92-parameter model** we do
(50 field modes parmtype 14 + 42 material groups parmtype 15; the group material
model *replaces* the ~13k per-module parmtype-7 parameters). Our smoke
reproduces exactly 50 + 42, confirming we are on that model.

Established channel roles, with the mechanism (slide 6):
`κ_cor/κ = 1 + A − εκ`, where A is the constant field scale and εκ is energy
loss, which grows at low pT.

| channel | role | why |
|---|---|---|
| J/ψ→μμ | absolute **scale** | stiff high-pT muons, εκ→0 isolates A; minimum-ionising, mass known to keV |
| soft displaced hadrons (his K⁰s) | field **shape** + **material** | displaced vertices, soft tracks curl through a wide (r,φ,z) range and respond strongly to the *local* field |
| cosmic μ | alignment / geometry | full-detector span; stiff single tracks, so ~nothing for B_z |
| NMR probes | boundary / absolute | genuinely independent; anchors r beyond the tracker |

**How this positions the bachelor kaon.** It is soft (median pT 0.53 GeV) and
displaced, so by the same mechanism the plausible targets are **field shape and
material, not scale** — scale stays with the J/ψ, which is not what this study
is after. The "single tracks add nothing to B_z" observation applies to
*cosmics*, which are **stiff**; it does not transfer to a soft curling track.

**Mass hypothesis, not mass constraint.** The kaon carries
`trackParticleName='kaon'` (0.4937 GeV), which sets its dE/dx and G4
propagation. There is no invariant-mass *constraint*, so the channel cannot set
the absolute scale — but the mass *hypothesis* is exactly what makes it a
species probe of the energy-loss model. At equal momentum a kaon and a pion have
different βγ and therefore different dE/dx.

**Expected size of the field↔material correlation.** Slide 9 reports the two
blocks decoupled (≈0) in his combined fit, and our own measurement agrees
(mean |rho| = 0.0008 marginalised / 0.0016 fixed). The correlation is therefore
reported as a diagnostic rather than as the headline; the more informative
observable is the **per-parameter constraint σ**, split into mode 0 (scale) vs
modes 1+ (shape) vs material — the analogue of his slide 5.

**Scope note.** K⁰s is deliberately excluded and the fit is run independently
(not joined to his combined fit), so that whatever B⁺→J/ψK⁺ contributes is seen
in isolation. Statistics are correspondingly modest: ~1.2M kaons for Run2016H
(~10M for F+G+H) after cuts, against 86M J/ψ and 177M K⁰s in the reference fit.
This is an exploratory measurement of a channel's behaviour, not a competitive
calibration.

## Open Questions

- Alignment is **not floated** in the reference fit, and this study follows that
  (Decision 5). With no separate alignment stage, alignment is held at the GT
  geometry (zero residual) — the same implicit choice. Worth confirming that the
  reference fit does the same rather than importing a cosmics-derived alignment.
- Whether a future iteration should fit B⁺→J/ψK⁺ **mass-constrained** (B mass on
  the J/ψ+K system, which the kaon cfi anticipates as a downstream GBL Lagrange
  term). That would make it scale-capable, but it is out of scope here.
