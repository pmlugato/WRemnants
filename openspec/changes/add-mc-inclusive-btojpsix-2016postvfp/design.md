## Context

The `btojpsik` calibration compares data ALCARECO to MC ALCARECO on
per-track observables (pT-response, dE/dx, vertex-fit residuals) in
the B → J/ψ + X final state. `launch-full-alcareco-2016postvfp`
produces the data leg. This proposal produces the MC leg for the same
era with the same reco.

Two forces shape the design:

1. **Data-MC release consistency.** Release-level RECO/ALCARECO
   differences (constants, module wiring, hit fitters) are
   calibration systematics that we cannot correct downstream.
   Ideal: build the MC in the same CMSSW as the data was reco'd
   (`CMSSW_10_6_17_patch1`). **Blocker discovered in local smoke
   test**: `10_6_17_patch1` does not ship
   `PythiaFilterMultiAncestor` (used by the GEN fragment). Fall
   back to `CMSSW_10_6_20_patch1` (the McM campaign release for
   this fragment; 3 Pythia/EvtGen patches ahead of the data-reco
   release, no tracking/RECO source diff).
2. **Content-level parity.** The raw-data ALCARECO writes a specific
   `outputCommands` block (7 resonance producers + 3 dE/dx maps +
   cloned track collection + auxiliary). MC ALCARECO must write the
   exact same branch set — otherwise per-branch histmaker code either
   crashes on MC or silently misses branches.

## Goals

- Produce 4×10⁸ generated events (phase-1; ~8.3×10⁶ accepted at
  ε≈2.08%) in `CMSSW_10_6_20_patch1`, GEN through ALCARECO in a
  single cmsRun per job. Sized to reach 1:1 MC:data on the priority
  channel B⁺→J/ψK⁺ (~4×10⁵ final-selection events at either side).
- Match raw-data ALCARECO event content byte-for-byte on the
  `ALCARECOTkAlJpsiX*` branches.
- Fix the McM fragment's missing anti-Λ_b entries so both particle
  and antiparticle Λ_b decays are forced into J/ψ Λ (Λ̄).
- Deliver reconcile hygiene equivalent to the raw-data launch
  (`find_missing.py`, per-job JSON, resume filter).

## Non-Goals

- Any dedicated per-channel MC production. Bs→J/ψφ, Bc→J/ψπ, dedicated
  B⁺→J/ψK⁺, and dedicated Λ_b→J/ψΛ are all documented in
  `mc_prod_preparation.md` but scoped to follow-up proposals so this
  proposal stays landable.
- Any change to the raw-data ALCARECO producers or output cff. This
  proposal imports them unchanged from the raw-data payload.
- pre-VFP or other-era MC. Postponed until postVFP MC-vs-data closure
  is demonstrated.

## Decisions

### D1. Same CMSSW release everywhere (`10_6_20_patch1`)

Alternatives considered:
- (a) `CMSSW_10_6_17_patch1` for everything (original plan). **Ruled
  out during implementation**: the data-reco release doesn't ship
  `PythiaFilterMultiAncestor` — the GEN fragment can't load.
- (b) Two-release setup: `10_6_20_patch1` for GEN-SIM-DIGI-HLT,
  then bridge to `10_6_17_patch1` for RECO-ALCARECO. Mirrors CMS
  central production practice.
- (c) `10_6_20_patch1` for everything (this proposal, adopted after
  (a) was ruled out).
- (d) `10_6_49_patch1` (latest McM tag) for GEN + `10_6_17_patch1`
  for RECO. Would pull in current Pythia/EvtGen bug fixes but at the
  cost of a bridge and unknown GEN-level diffs vs the frozen McM
  campaign tag.

Rationale for (c): keep it single-release, no bridge cfg to
maintain. Delta between 10_6_17 and 10_6_20 is 3 patches of
Pythia/EvtGen bug fixes; no changes to `RecoTracker/`, `Alignment/`,
`Configuration/AlCa/` (verified by inspecting the CMSSW release
notes / package diffs). MC-vs-data ALCARECO comparison remains
valid on the observables used by the calibration (track pT, dE/dx,
vertex residuals). Cost of (b): a second scram area, an
`edmCopyPickMerge`-style bridge between DIGI output and RECO input,
extra failure modes at the tier boundary. (c) is the simplest
correct choice given (a) is unavailable.

If a specific GEN-side bug in 10_6_20 is later found to bias any
observable in the calibration, we revisit either by escalating to
(d) or by patching just the offending module into the tarball.

### D2. Locally-checked-in patched fragment, not McM `curl`

Alternatives:
- (a) `sed`-patch the McM-served fragment inside `run.sh` at job
  start (the pattern used in `MCPRODTEST.sh`).
- (b) Locally-checked-in patched fragment, `cp` into place at
  build_tarball time (this proposal).
- (c) Register the patched fragment on McM under a new prepid.

Rationale for (b): the anti-Λ_b patch is a multi-block edit (Alias +
ChargeConj + Decay/CDecay + `list_forced_decays`) with textual
context that `sed` alone handles fragilely. Keeping the patched
fragment as a real file in the repo makes the diff reviewable in git
and lets the smoke-test parallel driver reuse the same file. McM
registration is overkill for a private production.

Downside: the fragment now has two sources of truth (McM version vs
our patched version). Task 1.4 mitigates by diff-auditing at each
build_tarball.

### D3. Fix anti-Λ_b instead of shipping a dedicated Λ_b fragment

The McM fragment defines `MyLambda_b0` (particle only) and includes it
in `list_forced_decays`. Anti-Λ_b (the antiparticle) is untouched: no
alias, no ChargeConj, no CDecay, not in `list_forced_decays`. Compare
Ξ_b and Ω_b, where both charge conjugates are aliased and forced.

Two options to make anti-Λ_b usable for the calibration:
- (a) Ship a dedicated Λ_b fragment forcing all Λ_b (particle+anti)
  into J/ψ Λ. Kills the ~98% non-Λ_b events at GEN → higher effective
  Λ_b yield per accepted event.
- (b) Fix anti-Λ_b in the inclusive fragment so it's forced alongside
  the particle. Effect: ~doubles the Λ_b contribution to filtered
  events (was 2.3% pilot, → ~4–5%). Marginally shortens all other
  channels by ~1–2%. (this proposal.)

Rationale for (b): the calibration wants a common inclusive sample
with all b-hadron channels present so per-channel truth-tags share the
same nuisance model. Dedicated Λ_b (option a) is a separate proposal
for a headroom-limited use case; not the base plan.

### D4. Single cmsRun full-chain per job (not multi-cmsRun with bridging)

Alternatives:
- (a) One cmsRun per job doing GEN→SIM→DIGI→HLT→RECO→ALCA (this
  proposal).
- (b) Two cmsRuns per job: GEN→DIGI, then bridge cfg RECO→ALCA.
- (c) Split into 3+ cmsRuns per stage.

Rationale for (a): simpler failure mode (one exit code), no
intermediate ROOT files on disk per job, no bridge cfg to maintain.
Cost: one long cmsRun per job (~7 h). Watchdog set to 10 h to catch
runaway processes without truncating legitimate long runs.

(b) becomes attractive only if the DIGI→RECO output would fit in
`$_CONDOR_SCRATCH_DIR` for extended debugging, which is not a current
need.

### D5. Seed strategy

Per-job initial seed = `base_seed + clusterId × 100000 + procId`,
injected by the wrapper. Recorded in the per-job JSON as `seed_used`.

Alternatives considered:
- (a) `Math.random()`-style non-reproducible seed. Rejected: makes
  resubmit-vs-original comparison for the determinism check
  (task 3.7) impossible.
- (b) Sequential seeds `procId`-only. Rejected: a resubmit of one
  cluster after another would collide seeds with the earlier
  cluster.

The chosen scheme guarantees (i) reproducibility — same
`(clusterId, procId)` always yields the same output; (ii) collision
avoidance across resubmits (a resubmit lives in a new cluster with a
different `clusterId`).

### D6b. Local-validate + hand off, do not run condor ourselves

The condor submission is executed by an external operator, not by us.
We produce a self-contained handoff package
(`mc_inclusive_handoff_<TS>.tgz`) that the operator unpacks and runs
via `./SUBMIT.sh <rehearsal|full>`. We locally exercise everything
that can be exercised without condor:

- Full-chain cmsRun on this workstation (multi-parallel, ~1000
  accepted events).
- Branch-content parity against raw-data ALCARECO.
- Seed determinism.
- Peak disk / RSS measurement per job.
- Anti-Λ_b patch effect on species fractions.

The `.sub` files we ship carry `RequestDisk` and `RequestMemory`
values derived from the local measurement (p99 + headroom), so the
operator's first-attempt condor sizing is grounded, not guessed.

What only the operator can validate:
- Real-world per-slot success rate.
- Site-level xrootd variance.
- Wall-clock in their condor pool.
- The reconcile/resubmit tooling under actual failure modes.

The rehearsal-then-full pattern shifts one step downstream: we sign
off before the operator submits the full 50k batch, based on their
100-job rehearsal reconcile report.

Alternative considered: we submit to our own local condor pool for
rehearsal, then hand off for full. Rejected: extra failure surface
(two condor pools, two site-list configurations) and we still can't
verify the operator's pool from ours. The local multi-parallel test
covers the same failure modes without the condor variable.

### D6. Reuse `jpsix_alcareco_payload.tgz` verbatim

The raw-data ALCARECO payload contains patched copies of
`Alignment/`, `Configuration/AlCa/`, `Configuration/EventContent/`,
`Configuration/StandardSequences/` that wire the
`ALCARECOTkAlJpsiX*Resonances` producers into the ALCA path. The MC
job needs the exact same wiring — otherwise the resonance branches
won't appear on MC. Reusing the tarball verbatim removes a class of
"MC and data have different ALCARECO content because someone forgot
to sync the payload" failures.

`build_tarball.sh` takes the raw-data payload as a build input, adds
the patched fragment under `Configuration/GenProduction/python/`, and
adds the full-chain cfg. No cff edits.

## Risks / Trade-offs

- **Risk**: `10_6_20_patch1` GEN-side bugs that were fixed in later
  10_6 patches bias the physics we care about (b-hadron production
  fractions, χ_cJ feed-down). **Mitigation**: task 6.2 re-runs
  `classify_jpsi.py` on 10⁶ accepted events; compare species
  fractions to PDG-driven expectations. If a >2σ deviation from PDG
  fractions turns up, escalate to option (a) of D1.
- **Risk**: HLT tag drift between MC and data. **Mitigation**:
  task 0.2 pins the HLT tag from the raw-data reco explicitly and
  fails the build if the choice can't be resolved.
- **Risk**: Premix library gets purged from cvmfs / DAS between
  build_tarball and full launch. **Mitigation**: task 0.3 caches the
  resolved LFN list under `manifests/premix_library.txt`. If a
  purge happens, resubmit with the DAS-cached list.
- **Risk**: 20,000 jobs saturates the operator queue. **Mitigation**:
  `submit_batch.sh` takes an `--limit` flag; operator batches
  N=5,000 at a time if their pool constraints require.
- **Trade-off**: 4×10⁸ hits 1:1 on B⁺→J/ψK⁺ but only ~0.5:1 for
  2:1 on the same channel; Bs→J/ψφ remains ~40× low; Bc→J/ψπ
  inaccessible. Accepted for this iteration; phase-2 scale-up
  and dedicated channel productions are follow-up proposals.

## Migration Plan

Not a migration — a new production. The existing
`btojpsix2016mcprod.sh` and `btojpsix2016mcprod_parallel.sh` are
retained as GEN-only smoke tests with the release bumped to
`10_6_20_patch1`. No user-visible interfaces move.

## Open Questions

- Exact postVFP-specific GT: `106X_mcRun2_asymptotic_v13` (McM UL16
  default) vs `106X_mcRun2_asymptotic_v17` (postVFP-specific). Resolve
  in task 0.1 by inspecting the raw-data reco tags.
- Whether the ALCARECO output should carry additional MC-only
  branches (`genParticles`, `PileupSummaryInfo`) for the calibration
  side. Not in the raw-data event content by construction; MC-only
  additions belong in a separate branch of the output cff (out of
  scope for this proposal — deferred to the histmaker consumer).
