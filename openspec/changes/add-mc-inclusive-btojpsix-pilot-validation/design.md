## Context

David_w has produced ~108 pilot jobs of the MC inclusive B→J/ψ+X 2016 postVFP
sample using our handoff (`condor/mc_inclusive_btojpsix_2016postvfp/`) with
his own five fixes applied on top. His patch is at
`/work/submit/david_w/ZMass/BtoJpsiX_MCprod/pietro_fixes/pietro_fixes.patch`
(213 lines, 5 files); the critical physics fix (ALCARECO output-module
split=1) was in his tree at the time the pilot ran, so the pilot outputs are
the "good sample", not a diagnostic sample. The outputs live at
`/ceph/submit/data/user/d/david_w/mz/mc/inclusive_btojpsix_2016postvfp/` and
we have read access.

David is the operator for the eventual full ~15 fb⁻¹ MC submission, so his
tree — including proxy path, `+AccountingGroup`, and transient site fences —
stays authoritative. The only cross-tree operator change is the DEST_DIR: the
full production writes to the shared group MC store
`/ceph/submit/data/group/cms/store/mc/...`, not our user areas. We tell
David about the DEST_DIR change manually.

Before launching, we need to (a) apply David's fixes to our own tree so we
stop drifting from the physics-correct state, retarget DEST_DIR, and (b)
prove the pilot batch is physics-healthy end-to-end using existing Stage-2
CVH refit and B⁺ reconstruction machinery. The parked OSG rehearsal
(`add-mc-inclusive-btojpsix-osg-submission`) is superseded.

## Goals / Non-Goals

**Goals**
- Our tree at `condor/mc_inclusive_btojpsix_2016postvfp/` matches David's
  fixed tree for the five physics/tuning fixes; `DEST_DIR` points at the
  group MC store.
- Verified physics health of David's pilot output: reconstructed J/ψ mass
  peak, reconstructed B⁺ mass peak after joint J/ψ+K vertex fit with
  χ²/ndof cut, gen-matched vs non-matched split showing signal cleanly.
- Overlay against data ALCARECOs at
  `/ceph/submit/data/user/p/pmlugato/mz/alcareco/full_2016postvfp/charmonium/Run2016*`
  for visual sanity.
- Written go/no-go recommendation for the full submission, landed in the
  slide deck.

**Non-Goals**
- No fullchain cfg physics changes beyond David's split=1 output-module fix.
- No new AlCaReco channel work; we use `TkAlJpsiXBPlusResonances` as it
  exists.
- No mZ-fit-side (tensor / rabbit) changes.
- No new preset; we run on preset B, nominal alignment.
- No hard MC-vs-data agreement threshold — visual comparison only.
- No payload push-back to David unless the pilot surfaces a *new* bug beyond
  his five fixes.
- No alcareco→nano converter in this change; deferred to a follow-up.

## Decisions

- **D1: Apply David's fixes via `patch -p1 --dry-run` first, then real
  application against `fixed_files/`-diffed hand-merges for anything the
  patch can't cleanly land.** For operator-identity fields (proxy path,
  `+AccountingGroup`, transient site fences) we keep our own local values
  since we only dry-run on our tree; David's identity is authoritative for
  the real submission. `DEST_DIR` gets one intentional edit on top of the
  patch: retarget to `/ceph/submit/data/group/cms/store/mc/...`.
  Alternatives: (a) wholesale copy from `fixed_files/` — rejected, would
  overwrite our operator-identity fields with David's (irrelevant for our
  dry-runs but noisy); (b) skip the patch and rewrite from David's README
  description — rejected, error-prone for the split-level cfg-generation
  logic.

- **D2: No payload rebuild unless we find a new bug.** David's tree already
  contains the split=1 fix and produced the pilot correctly. Rebuilding just
  to have our own copy adds churn and forks the artifact identity. If pilot
  validation surfaces a bug requiring code changes, then we rebuild and
  hand-off a new tarball with a fresh sha256.
  Alternatives: (a) proactive rebuild "just in case" — rejected, no benefit
  and forks artifact identity; (b) never rebuild — rejected, closes the
  door on shipping a follow-up fix.

- **D3: Reuse `add-jpsi-x-stage2-bplus-cvh` Stage-2 machinery unchanged.**
  Two-track CVH refit for J/ψ + single-track CVH refit for the kaon is
  exactly what that change delivers on `TkAlJpsiXBPlusResonances` inputs. No
  new stage-2 code should be needed. If the pilot ALCARECO branch layout
  differs from what stage-2 expects, fix stage-2 code in place rather than
  fork.
  Alternatives: (a) write a fresh minimal refit — rejected, duplicates work;
  (b) skip stage-2 and use the raw saved ALCARECO kinematics — rejected,
  discards the reason we're producing the sample.

- **D4: `applyIdealGeometry=False` for the Stage-2 CVH refit on this
  validation.** We're validating the production sample under nominal
  alignment, which is what phase-1 will run. Ideal geometry is a separate
  study (systematic scan, not a production gate).
  Alternatives: run both — rejected as scope creep; run ideal only —
  rejected, doesn't match what phase-1 delivers.

- **D5: B⁺ candidate formation = joint kinematic-vertex fit of the CVH J/ψ
  track pair and the CVH bachelor K± track**, done offline on the Stage-2
  CVH outputs. The vertex-fit χ²/ndof is the primary selection knob; cut
  chosen at fit time from the χ² distribution (start with a loose
  χ²/ndof < 10 or similar). ALCARECO carries only track refs + hits — no
  fits happen at that stage — so the joint fit is a new step in the
  pilot-validation tool. Prefer wiring in whatever existing offline
  kinematic-vertex-fit machinery the JpsiX chain already provides for
  CVH-refit tracks (survey at §5.1); if nothing suitable exists, land a
  small tool under `scripts/btojpsik/mc_pilot_validation/`.
  Alternatives: pure 4-vector sum, no χ² cut — rejected, user explicitly
  called out the χ² cut as the key selection; reuse a pre-CVH chi² —
  rejected, no fits happen at ALCARECO stage so there is nothing to reuse.

- **D6: Gen-matching by ΔR < 0.03** at both muons and the bachelor kaon to
  a true B⁺→J/ψK⁺ decay chain in the gen record (post-PythiaFilter).
  Alternatives: match by reconstructed mass window — rejected, biases the
  efficiency measurement; ΔR looser — rejected, admits combinatoric
  contamination in the "matched" bucket.

- **D7: Data overlay is required, not stretch.** With
  `/ceph/submit/data/user/p/pmlugato/mz/alcareco/full_2016postvfp/charmonium/Run2016*`
  already produced, the incremental cost is one RDataFrame pass on both
  legs. No hard agreement threshold — the plot is the artifact. Once the
  alcareco→nano pipeline lands (out of scope here), this pass gets replaced
  by the standard histmakers.
  Alternatives: skip data — rejected, misses a cheap sanity signal; require
  formal agreement — rejected, over-scoped for a validation-of-production.

- **D8: Pilot validation lives in a dedicated
  `scripts/btojpsik/mc_pilot_validation/` directory, not the main
  `btojpsik.py` histmaker.** Keeps the main pipeline clean and this
  validation as a one-shot analysis. When the alcareco→nano pipeline lands
  we retire this dir in favor of the standard histmakers.
  Alternatives: bolt onto `btojpsik.py` — rejected, mixes calibration
  pipeline with production validation; put in `scripts/plotting/` —
  rejected, the scanner+fit+plotter are one unit and should live together.

- **D9: OSG rehearsal directory is parked with a README banner but not
  deleted.** The code is reusable if we later want OSG again; deleting
  forfeits the reconciliation work in
  `add-mc-inclusive-btojpsix-osg-submission`. Cluster 3228825 gets
  `condor_rm`'d after user confirmation.
  Alternatives: `git rm -r` — rejected, loses the OSG-specific submit-attribute
  reference for the future.

- **D10: Format for consuming CVH outputs — whatever is fastest to get to
  plots.** User signaled indifference since the alcareco→nano pipeline will
  standardize this later. Default choice: RDataFrame directly on the
  Stage-2 ROOT trees (fits the wider histmaker idiom and needs no new
  flattener).
  Alternatives: flat NanoAOD-style ntuple via a separate flattener —
  rejected for this one-shot, no need to build infrastructure that gets
  retired.

## Risks / Trade-offs

- **Risk**: The joint J/ψ+K kinematic-vertex fit tool may need to be built
  from scratch — ALCARECO has no fits, and the Stage-2 CVH refit itself
  operates on individual tracks (single or ditrack), not the full B⁺
  three-track topology.
  Mitigation: survey the JpsiX chain (`add-jpsi-x-vertex-fit-and-low-pt`
  and neighbors) at §5.1 to reuse; failing that, land a small offline tool
  under `scripts/btojpsik/mc_pilot_validation/`. Timeline impact: this
  becomes the largest single work item in the change if nothing is reusable.

- **Risk**: David's pilot ceph area is under `/data/user/d/david_w/`; if his
  quota tightens or he cleans it, our validation loses inputs.
  Mitigation: pin a small (~1–2 file) copy on our side for the plot-making
  step; keep pointing at his area for the reconcile counts.

- **Trade-off**: Running only the B⁺→J/ψK⁺ channel now (not the full
  inclusive X) gives us one focused signal but leaves the other J/ψ+X
  channels un-validated. Accepted — if B⁺→J/ψK⁺ is clean the others inherit
  confidence from sharing the same ALCARECO substrate.

## Migration Plan

1. Snapshot the current `mc_inclusive_btojpsix_2016postvfp/` tree
   (`git stash --include-untracked` or a tag) before applying David's patch.
2. `patch -p1 --dry-run` inside the handoff dir; iterate on conflicts.
3. Apply the patch; retarget `DEST_DIR` to the group MC store; verify each
   of the 5 fixes visually against `PIETRO_FIXES_README.md`.
4. Run the pilot validation on David's ceph outputs.
5. If green: DM David the DEST_DIR change and any additional operational
   notes; land the report in the slide deck; kick off the full ~15 fb⁻¹
   submission (in a follow-up change).
6. If red: open a follow-up change with the specific issue, rebuild the
   payload only if a code fix is warranted, and ship it back to David.

## Open Questions

1. Does the JpsiX chain (`add-jpsi-x-vertex-fit-and-low-pt` and its
   neighbors) already expose an offline joint kinematic-vertex-fit tool
   that consumes CVH-refit tracks and returns χ²/ndof? If yes, reuse; if
   no, land a small offline tool as part of this change.
2. Exact subpath under `/ceph/submit/data/group/cms/store/mc/` for this
   sample — the proposal assumes `inclusive_btojpsix_2016postvfp/`;
   confirm before shipping the DEST_DIR note to David.
