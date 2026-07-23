# Tasks: finalize CVH producer in `CMSSW_15_0_19_patch2`

Prerequisite: proposal.md reviewed and approved. Clarifying questions
in the proposal's "Open questions" section resolved before starting
§2 (cfi default flip scope) and §6 (per-era vs per-lumi GT check).

---

## 1. Consolidate development in `CMSSW_15_0_19_patch2/`

- [ ] **1.1** Add a scope note to
  `slides/cvh_producer_validation.tex` declaring
  `CMSSW_15_0_19_patch2` the development home and `CMSSW_10_6_26`
  read-only reference. No new README file is created.
- [ ] **1.2** Mirror any producer-level edits (`.cc` / `.h` under
  `src/Analysis/HitAnalyzer/plugins/`) that landed in `10_6_26/`
  during the debug study into `15_0_19_patch2/`:
  - `diff -r CMSSW_10_6_26/src/Analysis/HitAnalyzer/plugins CMSSW_15_0_19_patch2/src/Analysis/HitAnalyzer/plugins` to enumerate divergence;
  - port any 10_6_26-only edits forward; skip differences that
    are 15_0_19-only or reflect legitimate CMSSW-version drift.
  - Report the diff summary in the deck's scope frame.

## 2. Adopt option (B) as btojpsik default — `useIdealGeometry=False` + `corFiles=[]`

**Scope pivot (2026-07-02, mid-implementation)**: option (A)
was falsified during smoke-test — v721's `parmtree` has 108,332
entries while our producer's parmset has 95,392, and the
maker's `corFiles` loader (`ResidualGlobalCorrectionMakerBase.cc:1017`)
asserts strict one-to-one alignment. See §8 follow-up. §2 now
targets option (B) — physics-quality via aligned geometry, no
Stage-2 corrections.

**Scope guard**: cfi defaults are NOT changed. Non-btojpsik
launchers are §7-inventoried and flagged, never modified.

- [x] **2.1** In
  `CMSSW_15_0_19_patch2/src/Analysis/HitAnalyzer/test/runCvhBplusJpsiK.py`,
  add a `corFiles` VarParsing option (list-of-strings), default `[]`.
  The `_resolve_default_corfile()` helper is retained as a
  callable (not the default value) so a future §8 follow-up
  can activate (A) without CLI grammar changes.
- [x] **2.2** Change the `opts.register('useIdealGeometry', True, ...)`
  default from `True` to `False`. Forward the resolved
  `useIdealGeometry` and `corFiles` into both `.clone()` calls
  explicitly.
- [ ] **2.3** For every additional btojpsik-scoped launcher
  enumerated by §7 (`runCvhJpsi.py`, `runCvhJpsiXSmoke.py`), apply
  the same defaults: `useIdealGeometry=False`, `corFiles=[]`. Do
  NOT modify non-btojpsik launchers (David's
  `runCvhJpsiCandidateDriven*.py` already runs option (B); other
  channels' Ks/Lambda/D0 refits already run option (B); alignment-
  training configs like `RunGlobalCorRecJpsiAlca.py` stay at
  their explicit `True` setting).
- [x] **2.4** cfi defaults stay unchanged:
  `python/ResidualGlobalCorrectionMakerTwoTrackG4e_cfi.py` and
  `python/ResidualGlobalCorrectionMakerG4e_cfi.py` remain at
  `useIdealGeometry=cms.bool(True)`, `corFiles=cms.vstring()`.
  Verified by the §7 inventory tsv.
- [x] **2.5** In
  `CMSSW_15_0_19_patch2/src/Analysis/HitAnalyzer/plugins/ResidualGlobalCorrectionMakerBase.cc`,
  add a one-line `LogInfo("ResidualGlobalCorrectionMakerBase")`
  at construction time echoing both flags (LogInfo added at
  line 219; MessageLogger include added at line 100).
- [x] **2.6** `scram b -j8` clean build under `cmssw-el9`; only
  pre-existing Eigen `-Warray-bounds` warnings, no new warnings
  from the LogInfo addition.
- [ ] **2.7** 200-event B+ chain smoke under the new option-(B)
  default (no CLI overrides). Confirm the `LogInfo` reports
  `useIdealGeometry=false, corFiles.size()=0` and the physics
  metrics reproduce the previously-recorded σ ≈ 40 MeV, tail 0/386,
  median `chisqval` ≈ 30, `niter` cap-hit ≈ 1.3%. **Also record
  `edmval` median and `edmval < 1e-5` fraction** for the record.
- [x] **2.9** v721 provenance confirmed. `git log` shows a single
  2023-07-12 migration commit; `wremnants/production/muon_calibration.py:29`
  uses it as the mZ analysis's data-side correction file
  (`data_filename` in `make_muon_calibration_helpers`). The file
  is structurally sound (108,332-entry parmtree + idxmaptree).
  **But it does not fit our current producer's parmset (95,392
  entries)**; using it via the current `corFiles` loader crashes
  the maker at the `assert(nparms == parmset.size())`. Confirmed
  by the crash log in
  `/ceph/.../cvh_outputs/an_canonical_smoke/cvh_bplus_dimuon_an_canonical.log`.
  Punt to §8 follow-up.
- [ ] **2.10** Rerun the smoke with new (B) defaults; confirm
  physics-quality metrics as in 2.7.

## 3. Kaon-side ideal-geometry × corFiles study + (A) vs (B) cross-check

- [ ] **3.1** Prepare inputs: one Run2016H preset-B alcareco file
  (same one as the dimuon smoke). Confirm it has both dimuon and
  bachelor-kaon RefTracks populated (already validated in the smoke).
- [ ] **3.2** Expose per-leg knobs in the driver:
  `useIdealGeometryMuon`, `useIdealGeometryKaon`, `corFilesMuon`,
  `corFilesKaon` (each defaults to the base value if unset). This
  lets §3.3 and §3.4 vary muon and kaon legs independently on the
  same event sample.
- [ ] **3.3** Run the three primary points on the same file:
  - **P1 (AN-canonical)**: muon `ideal + v721`, kaon `ideal + v721`
  - **P2 (aligned baseline)**: muon `aligned + []`, kaon `aligned + []`
  - **P3 (broken baseline)**: muon `ideal + []`, kaon `ideal + []`

  `plimit=0.05` (confirmed by user, keeps kaon soft tail). Output
  tuple encoded in filename. Each point takes ≈5 min on 200 events.
- [ ] **3.4** Within P1 (AN-canonical), run the 2×2 muon-geom × kaon-geom
  sub-scan keeping the correction file loaded on both legs (i.e.
  toggle only geometry, not `corFiles`). Four sub-points; four
  outputs. This isolates residual asymmetry between the two legs
  once corrections are applied.
- [ ] **3.5** Join each output with the ALCAReco B+ side (§8
  says join overhead is out of scope; use existing join script
  once per point).
- [ ] **3.6** Compute per point:
  - kaon `Δp_T/p_T` median and MAD (refit vs ALCAReco);
  - kaon `Δ(q/p)` MAD;
  - `chisqval` median and cap-hit fraction;
  - `edmval` median, q95, and fraction meeting `edmval < 1e-5`
    (reported alongside `chisqval` — both are convergence diagnostics
    and disagreeing hints between them are informative);
  - μμ mass σ on [2.95, 3.25] and tail fraction `|m − 3.097| > 0.20`;
  - B+ candidate mass σ on [5.15, 5.40] and tail fraction
    `|m − 5.279| > 0.20`.
- [ ] **3.7** Write `scripts/btojpsik/plot_kaon_ideal_geom_compare.py`
  producing:
  - Primary panel: P1 / P2 / P3 side-by-side B+ candidate mass
    overlay + metric table;
  - Sub-panel: 2×2 muon-geom × kaon-geom within P1 (kaon Δ(q/p)
    scatter/hist per cell + B+ mass σ table);
  - Interpretation bullet: does the kaon leg under `ideal +
    empty corFiles` look better than the muon leg under the same
    config? If yes, quantify — the ratio
    `(alignment displacement / intrinsic σ(q/p))² × N_hits`
    should explain the observed asymmetry.

## 3b. (A) vs (B) direct overlay

- [ ] **3b.1** On the same file used in §3, produce a direct
  overlay of (A) `useIdealGeometry=True + corFiles=[v721]` and
  (B) `useIdealGeometry=False + corFiles=[]` for both μμ mass
  and B+ candidate mass.
- [ ] **3b.2** Fit each with a Gaussian; report σ, mean, tail
  fraction per branch. Also report per-branch `chisqval` median
  + `edmval` median + `edmval < 1e-5` fraction so the (A) vs (B)
  comparison captures convergence quality alongside mass shape.
- [ ] **3b.3** If (A) and (B) agree within statistical error on
  Run2016H, the AN's claim holds and we can consider (A)
  production-ready. If they disagree, document the disagreement
  in the deck and flag Open Question #4 (correct `correctionResults_*`
  file for our sample) as needing follow-up before full-scale
  production.
- [ ] **3b.4** Include the overlay in
  `scripts/btojpsik/plot_kaon_ideal_geom_compare.py` (same script
  as §3.7, additional figure).

## 4. Four-flavor B+ mass overlay

- [ ] **4.1** Write `scripts/btojpsik/plot_bplus_mass_overlay.py`:
  reads the joined-CVH ROOT with `useIdealGeometry=False`
  (both legs), computes the four mass definitions per candidate,
  overlays on one axis on [5.10, 5.45] with 100 bins.
- [ ] **4.2** Fit each of the four with a Gaussian on [5.20, 5.36],
  report μ, σ, and integral / total-tail fraction on the plot.
- [ ] **4.3** Verify (in a docstring, not in production code):
  under preset B the ALCAReco output candidate mass equals the
  raw 4-vector-sum mass to floating-point precision (no
  kinematic constraint applied by `JpsiXCandidateProducer` for
  preset B). Under preset C they differ by the J/ψ-constraint
  effect. If we're only running preset-B files this quarter,
  document the collapse of the "raw" and "ALCAReco output"
  flavors and either drop one or keep both for the record.
- [ ] **4.4** Sanity: the CVH-with-mass-constraint σ should equal
  the no-constraint σ within ~few MeV (the constraint only tightens
  the μμ subsystem, not the kaon), while the mean shifts by the
  refit systematic. Document.

## 5. Publication surface

- [ ] **5.1** Extend `slides/cvh_producer_validation.tex`:
  - new frame "AN two-stage pipeline explained" summarising the
    AN §3.3–3.4 architecture (Stage 1 = CVH with ideal geometry
    writing Jacobians; Stage 2 = linearised correction from
    `corFiles`) and pointing to `correctionResults_v721_recjpsidata.root`
    as the file we adopt;
  - new frame "P1 / P2 / P3 comparison" with §3 plot + metric
    table + the "kaon looks better under P3" explanation;
  - new frame "(A) vs (B) direct overlay" with §3b plot + fit
    numbers;
  - new frame "B+ mass, four flavors" with §4 overlay + fit table;
  - update pre-production checklist frames to reflect the outcome
    of §5 (launcher inventory) and §6 (GT verification).
- [ ] **5.2** Rebuild the PDF; verify no overfull hbox/vbox
  overflows past the printable area (rasterise the added frames
  with pdftoppm as we did during the previous polish pass).
- [ ] **5.3** Dump the two new figures (PNG) into
  `~/public_html/mz/cvh/` and confirm they render via the
  `submit.mit.edu/~pmlugato/` URL.

## 6. Verify GT consistency

- [ ] **6.1** Pick one ALCAReco file per era in use (2016H, 2017,
  2018 as applicable) and run
  `edmProvDump <file>.root | grep -A2 -i 'globalTag'` (or the
  equivalent). Record the GT string per era.
- [ ] **6.2** Compare against the GT hard-coded in the CVH
  driver `runCvhBplusJpsiK.py` for that era.
- [ ] **6.3** If any era mismatches, either:
  (a) pin the driver's GT to the ALCAReco RECO GT, or
  (b) verify that the alignment records referenced by both GTs
  point to identical IOVs. Prefer (a) if in doubt.
- [ ] **6.4** Table in `slides/cvh_producer_validation.tex`:
  (era, ALCAReco GT, CVH driver GT before, CVH driver GT after,
  action). If all match, one row per era saying "match — no
  action".

## 7. Inventory all cmsRun launchers

- [ ] **7.1** Write `scripts/btojpsik/inventory_cvh_launchers.sh`:
  greps `CMSSW_15_0_19_patch2/`, `CMSSW_10_6_26/`, top-level
  `scripts/`, `condor/`, `old_condor_stuff/`, and `cooper/` for
  references to `ResidualGlobalCorrectionMakerTwoTrackG4e`,
  `ResidualGlobalCorrectionMakerG4e`, `useIdealGeometry`, or
  `corFiles` (including `process.load()` paths that inherit
  the cfi default). Prints one line per file:
  path, mode (dimuon/singletrack/both), effective
  `useIdealGeometry` (explicit True / explicit False / default),
  effective `corFiles` (explicit list / empty / default), and a
  suspected purpose (btojpsik / alignment-training / other /
  unknown).
- [ ] **7.2** Classify each row into:
  - **btojpsik-scoped** — driver gets §2 wiring
    (`useIdealGeometry=True` + `corFiles=[v721]`);
  - **not btojpsik-scoped** — flag in the inventory table,
    do NOT modify. Annotate with the suspected workflow
    (alignment training vs other) and whether the launcher is
    running with `corFiles` populated or empty. An empty-`corFiles`
    non-btojpsik launcher is flagged as "may be broken like
    ours was" for the owner to review; a populated-`corFiles`
    non-btojpsik launcher is flagged as "AN-canonical for
    another workflow — leave alone".
- [ ] **7.3** Add the inventory as a table in
  `slides/cvh_producer_validation.tex` with columns:
  path, mode, `useIdealGeometry`, `corFiles`, classification,
  action taken. Any non-btojpsik launcher with an unusual config
  (either `useIdealGeometry=False` while running alignment
  training, or empty `corFiles` while running production) gets
  a "flag — do not change; review with owner" note.

## 8. Wrap-up

- [ ] **8.1** All figures dumped to `~/public_html/mz/cvh/`,
  linked from the slide-deck LaTeX log for reproducibility.
- [ ] **8.2** Update `memory/project_cvh_jpsi_mass_broadening.md`
  with the final validation numbers (kaon-side confirmation of
  the alignment-displacement hypothesis, four-flavor B+ mass
  agreement, GT check outcome).
- [ ] **8.3** Signal readiness for the David-fork diff pass
  (separate change-id).
