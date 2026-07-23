# Change: Add inclusive B → J/ψ + X MC production for 2016 postVFP, matched to raw-data ALCARECO

## Why

`launch-full-alcareco-2016postvfp` will produce raw-data ALCARECO for
the 2016 postVFP era. To make MC-vs-data momentum-scale calibration
apples-to-apples, we need an inclusive B → J/ψ + X MC sample:

- generated in the **same CMSSW release** (`CMSSW_10_6_20_patch1`) as
  the raw-data ALCARECO pipeline;
- produced end-to-end through
  **GEN → SIM → DIGI-Premix-HLT → RAW2DIGI-RECO → ALCARECO**;
- with **the same `ALCARECOTkAlJpsiX*` event content** the raw-data
  ALCARECO writes (7 resonance producers + 3 dE/dx maps + cloned track
  collection);
- **2:1 MC:data** on the workhorse inclusive channels
  (B⁺→J/ψK⁺, B⁰→J/ψK*⁰, B⁰→J/ψKS, Λ_b→J/ψΛ, ψ(2S)→J/ψππ).

Two documented issues in the current setup must be fixed before scale:

1. **Release choice — corrected during implementation.** The
   original plan targeted `CMSSW_10_6_17_patch1` (the raw-data
   ALCARECO release) end-to-end. Blocker discovered in local smoke
   test: `10_6_17_patch1` does NOT ship `PythiaFilterMultiAncestor`
   (used by the GEN fragment) — the data-reco release strips all
   MC-side GEN filter plugins. Fall-back adopted: use
   `CMSSW_10_6_20_patch1` for the whole MC chain. This is the McM
   campaign release for `BPH-RunIISummer20UL16GEN-00017`. Delta vs
   data-reco = 3 patch releases (Pythia8/EvtGen bug fixes); no
   tracking/RECO/ALCARECO source changes between the two, so
   calibration comparison stays apples-to-apples on the observables
   this analysis uses.
2. The McM fragment (`BPH-RunIISummer20UL16GEN-00017`) is missing
   `Alias`/`ChargeConj`/`CDecay`/`list_forced_decays` entries for
   `anti-Lambda_b0`. Effect: ~half the Λ_b sample (the anti-baryons)
   fall through to Pythia's default decay tables and are not forced
   into `J/ψ Λ̄`. Compare Ξ_b/Ω_b: both charge conjugates are present.
   Fix in a locally patched fragment.

Dedicated per-channel productions (Bs→J/ψφ, Bc→J/ψπ, dedicated Λ_b)
are **out of scope**. They will be follow-up proposals.

### Operational model

**We do not run condor.** This proposal produces a self-contained
handoff package that an external operator runs on their side.
Concretely:

- **We do** locally on this machine: build the payload, run the full
  chain multi-parallel on ~1000 accepted events, validate ALCARECO
  branch content, measure per-job disk/RSS footprint, confirm
  seed determinism.
- **We produce** a `condor/mc_inclusive_btojpsix_2016postvfp/`
  directory containing everything the operator needs: payload
  tarball, `.sub` files with tuned RequestDisk/RequestMemory,
  `run.sh` worker, `submit_batch.sh` driver, reconcile tooling, and
  a one-line entry-point script (`SUBMIT.sh`).
- **The operator does**: unpack, run `SUBMIT.sh <rehearsal|full>`,
  drain, run reconcile, resubmit failures, report back.

We do NOT self-verify condor mechanics (queue behavior, per-slot
success rates, real-world wall-clock, xrdcp bandwidth at scale).
Those are validated by the operator on the first rehearsal batch.

## What Changes

### Scope (locked in with the user)

- **Physics**: inclusive B → J/ψ + X, 2016 postVFP MC.
- **Volume**: **4×10⁸ generated events (phase-1)** at ε ≈ 2.08%
  ⇒ ~8.3×10⁶ accepted events through GEN filter. Sized to reach
  **1:1 MC:data on the priority channel B⁺→J/ψK⁺** (~4×10⁵
  reconstructible at final-selection vs ~5×10⁵ data yield).
  Expected per-channel yields at final-selection stage after
  ALCARECO + downstream cuts: B⁺→J/ψK⁺ ~4×10⁵, B⁰→J/ψK*⁰ ~3×10⁵,
  Λ_b→J/ψΛ→pπ ~2×10⁵ (comfortably 2:1 with the anti-Λ_b fix),
  ψ(2S)→J/ψππ ~4×10⁵. Scaling to 10⁹ (**phase-2**, +6×10⁸ on top
  for 2:1 on B⁺→J/ψK⁺) is out-of-scope for this proposal —
  evaluated as a follow-up only if the phase-1 calibration
  demonstrates statistics-limited systematics. Bs→J/ψφ falls
  short at either scale (~40× low even at 10⁹) and Bc→J/ψπ is
  inaccessible from inclusive; both deferred to dedicated
  fragment proposals.
- **CMSSW release**: `CMSSW_10_6_20_patch1` for the entire chain
  (GEN through ALCARECO). Matches the raw-data ALCARECO release.
- **Beamspot / GT / era**: `Realistic25ns13TeV2016Collision`,
  **`106X_mcRun2_asymptotic_v17`** (MC counterpart of the
  `106X_dataRun2_v35` used by the raw-data ALCARECO run), `Run2_2016`.
  The McM prepid ships `_v13`; we bump to `_v17` to pick up the
  postVFP-specific conditions.
- **HLT simulation**: `--step ...,HLT:@relval2016,...`. Runs the 2016
  RelVal HLT menu on the digi'd events and writes
  `TriggerResults_HLT_*` into the RAW payload so downstream analysis
  can apply the same HLT bit selection to MC that it applies to data.
  This is a *simulation*, not a selection: no events are dropped by
  the HLT step. Skipping it would leave MC ALCARECO without a
  `TriggerResults` branch and violate the "content parity with
  data ALCARECO" requirement.
- **Pileup**: Premix, `RunIISummer20UL16` UL16 minbias premix
  library. Pin the exact secondary-input filelist at build_tarball
  time.
- **Event content**: same as the raw-data ALCARECO
  (`ALCARECOTkAlJpsiX_Output_cff.OutALCARECOTkAlJpsiX`) **plus** an
  MC-only `keep PileupSummaryInfos_*_*_*` line so downstream can
  PU-reweight against data. The MC-only addition is documented
  explicitly in the output cff, not hidden.
- **Output layout**:
  `/data/user/p/pmlugato/mz/mc/inclusive_btojpsix_2016postvfp/`,
  written via `root://submit50.mit.edu/` xrdcp.
- **Retry**: manual reconcile via `find_missing.py`, same shape as
  raw-data launch. No DAGman.

### New launch directory `condor/mc_inclusive_btojpsix_2016postvfp/`

Mirrors `condor/full_alcareco_2016postvfp/`, adapted to MC (no input
LFN list; per-job seed strategy replaces per-job input LFN).

```
condor/mc_inclusive_btojpsix_2016postvfp/
├── fragments/
│   └── BPH-RunIISummer20UL16GEN-00017-fragment.py    # patched (anti-Λ_b + annotation)
├── build_tarball.sh              # CMSSW area + patched fragment + full-chain cfg
├── make_fullchain_cfg.sh         # cmsDriver: GEN→SIM→DIGI→HLT→RECO→ALCARECO
├── run.sh                        # worker (seed handling, cmsRun, xrdcp out)
├── submit_batch.sh               # renders queue of N_JOBS jobs w/ distinct seeds
├── find_missing.py               # reconcile: expected {seed} vs successful outputs
├── build_manifest.py             # per-cluster manifest
├── merge_manifests.py            # roll-up
├── quota.sh                      # ceph headroom preflight
├── manifests/                    # per-cluster JSONs, cluster logs
├── submits/                      # per-batch .sub files
└── README.md                     # runbook
```

Rationale for a dedicated launch directory: MC job shape is materially
different from the data ALCARECO job (no input LFN; DIGI premix is the
disk-heavy step; single job runs the full chain). Reusing the raw-data
directory would blur two operational contracts.

### Patched fragment

Kept in `fragments/BPH-RunIISummer20UL16GEN-00017-fragment.py`. Diffs
against the McM-served version:

- **Add `Myanti-Lambda_b0` alias, ChargeConj, and forced-decay entry**:
  ```
  Alias      Myanti-Lambda_b0 anti-Lambda_b0
  ChargeConj MyLambda_b0      Myanti-Lambda_b0
  ...
  Decay MyLambda_b0
  0.44 MyJ/psi Lambda0 PHSP;
  0.56 Mypsi(2S) Lambda0 PHSP;
  Enddecay
  CDecay Myanti-Lambda_b0
  ...
  list_forced_decays = cms.vstring(
      ...
      'MyLambda_b0',
      'Myanti-Lambda_b0',   # NEW
      ...
  )
  ```
- **Fix stale annotation**: replace the misleading
  `'Jpsi->mumu (no  kin cuts on muons)'` string with the correct
  filter descriptor.

The `curl`-from-McM step in `btojpsix2016mcprod.sh` is replaced by a
`cp` from `fragments/`. The `request_fragment_check.py` step is
skipped for the patched fragment (it validates against McM, which
doesn't know about the patched version).

### Full-chain cmsDriver

Single-cmsRun end-to-end chain: `GEN,SIM,DIGI,L1,DIGI2RAW,HLT:@relval2016,RAW2DIGI,L1Reco,RECO,PAT,ALCA:TkAlJpsiX_step`.

The ALCARECO step must schedule
`Alignment.CommonAlignmentProducer.ALCARECOTkAlJpsiX_cff` and write the
`OutALCARECOTkAlJpsiX_Output` event content. Uses the same patched
`Alignment/`, `Configuration/AlCa/`, `Configuration/EventContent/`,
`Configuration/StandardSequences/` payload that the raw-data ALCARECO
run carries (`jpsix_alcareco_payload.tgz` — reuse verbatim). This
guarantees the 7 per-channel resonance producers, the 3 dE/dx maps,
and the cloned/deduplicated track collection are all present in the
MC ALCARECO output.

### Per-job structure

- **1 condor job = 1 cmsRun** invocation running the full chain for
  `N_gen_per_job` generated events, producing `~N_gen × 2.08%`
  accepted events.
- Target per-job wall-time: ≤ 8 h. At ~55 s/accepted × ~500 accepted
  + ~0.8 s × ~24k generated ≈ 3.4 h + 5.3 h ≈ 8.7 h — trim to
  N_gen_per_job = 20,000 (tuned off local validation) to land at
  ~7.5 h.
- **Phase-1 total: 4×10⁸ / 2×10⁴ = 20,000 jobs.**
  CPU budget ~240k CPU-h; ~3 wall-days on 3k operator slots.
  ALCARECO storage ~1.9 TB (200 kB/evt × 8.3×10⁶ accepted +
  MC-only PU branch overhead).
- **Per-job seed**: deterministic from `(clusterId, procId)` via
  `initialSeed = base_seed + clusterId × 100000 + procId` with a
  fixed `base_seed = 20260716` recorded in
  `manifests/base_seed.txt`. Reproducible; no seed collisions across
  a resubmit; a year-N rerun reproduces bit-identical outputs given
  the same `(base_seed, clusterId, procId)` triple.

### Worker script (`run.sh`)

Same skeleton as
`condor/full_alcareco_2016postvfp/run.sh`, modified for MC:

- **No stage-in**. No input LFN; no `xrdcp` on the input side.
- **Seed injection**: worker appends a
  `process.RandomNumberGeneratorService.generator.initialSeed = ...`
  block to the cfg based on `clusterId` and `procId` at job start.
- **cmsRun** runs the full chain to a single ALCARECO ROOT output.
- **`xrdcp` output** to
  `/data/user/p/pmlugato/mz/mc/inclusive_btojpsix_2016postvfp/`.
- **Wall-time watchdog**: flat **10-hour cap** (larger than the
  raw-data 2h cap because the MC chain is heavier per accepted event).
- **Structured exit tags**: `ok | cmsrun_crash | xrdcp_out_failed |
  wrapper_error | wrapper_timeout`. (No `input_stagein_failed`
  bucket — no input file.)
- **Per-file JSON**: adds `seed_used`, `n_generated`, `n_accepted`,
  `filter_efficiency`, plus the standard timing/RSS/host fields.

### Reconcile contract

- Expected work = the closed set of `(clusterId, procId)` pairs
  scheduled (recorded in the queue file).
- `find_missing.py` compares expected vs the manifest of successful
  outputs (ceph output + JSON exit 0), buckets misses into
  `unrun ∪ failed_cmsrun ∪ failed_xrdcp ∪ failed_wrapper`.
- Auto-safe resubmit = `unrun ∪ failed_xrdcp ∪ failed_wrapper`.
  `failed_cmsrun` is held for manual triage (same rule as raw-data
  launch).
- Resume filter: `submit_batch.sh` skips
  `(clusterId, procId)` pairs already recorded successful.

### Local pre-handoff validation (mandatory)

Run on this machine with `K` parallel single-threaded cmsRun processes
(same shape as `btojpsix2016mcprod_parallel.sh`), targeting **~1000
accepted events total** (~50k generated at ε≈2.08%). Validate:

- **(a) branch content**: `edmDumpEventContent` on one output file
  lists every branch in `ALCARECOTkAlJpsiX_Output_cff.py` (7
  resonance producers, 3 dE/dx maps, cloned track collection,
  `TriggerResults_HLT_*`, `offlinePrimaryVertices`, auxiliary).
- **(b) parity with data ALCARECO**: diff the `edmDumpEventContent`
  output against a raw-data ALCARECO file from
  `condor/full_alcareco_2016postvfp/` on the `ALCARECOTkAlJpsiX*`
  prefix — must be empty.
- **(c) per-job scratch footprint**: measure peak disk usage across
  the local runs. Set condor `RequestDisk` = p99 + 2× headroom.
- **(d) per-job RSS**: measure peak RSS. Set condor `RequestMemory`
  = p99 + 1 GB headroom.
- **(e) seed determinism**: run the same seed twice, byte-diff the
  ROOT output payload branches (metadata timestamps excepted) — must
  be identical.
- **(f) anti-Λ_b fix effect**: rerun `classify_jpsi.py` — Λ_b
  species fraction rises from ~2.3% to ~4–5%.

The `RequestDisk` / `RequestMemory` values measured here are the
values baked into the `.sub` files handed to the operator.

## Impact

- Affected specs:
  - `mc-production` (new capability) — MC launch mechanics: patched
    fragment provenance, CMSSW-release consistency, full-chain job
    contract, seed strategy, ALCARECO content match, reconcile
    contract.
- Affected code:
  - `condor/mc_inclusive_btojpsix_2016postvfp/` — new directory.
  - `btojpsix2016mcprod.sh` — retained as a standalone GEN-only
    smoke test; version bumped `10_6_20_patch1 → 10_6_17_patch1`.
    Sourced from the patched fragment, not McM.
  - `btojpsix2016mcprod_parallel.sh` — same release bump.
  - Reuses `jpsix_alcareco_payload.tgz` from
    `condor/full_alcareco_2016postvfp/` verbatim.
- No changes to `wremnants/`, `narf/`, `rabbit/`, or the analysis
  pipeline.

## Non-goals

- **Executing condor jobs from this workstation.** We produce the
  handoff package and validate locally; the operator runs it.
- Storing intermediate SIM / DIGI / RECO outputs (would double
  storage; user declined). One cmsRun per job, ALCARECO out only.
- Dedicated per-channel MC productions (Bs→J/ψφ, Bc→J/ψπ, dedicated
  Λ_b, dedicated B⁺→J/ψK⁺). Written up as follow-up proposals.
- Multiple MC eras. Only 2016 postVFP this iteration.
- Modifying the raw-data ALCARECO producers or event content.
- Downstream stage-2 CVH consumption of the MC ALCARECO outputs
  (handled by the CVH change family).
- Central McM re-registration of the patched fragment (private
  production only).

## Open questions (to resolve during apply)

1. **Premix library snapshot**. `dasgoclient` query at build_tarball
   time to pin the UL16 postVFP minbias premix DAS-resolved LFN list;
   cache under `manifests/premix_library.txt` for reproducibility.
   Escalates to the operator only if DAS returns empty.
2. **Per-job N_gen**. Provisional 20,000 gen/job (≈ 7.5 h wall
   projected). Tuned off the 100-job rehearsal (§3.6) before the
   50k-job full submit.

Resolved during scoping:
- Release strategy → same `10_6_17_patch1` for the full chain.
- MC GT → `106X_mcRun2_asymptotic_v17` (postVFP-specific MC
  counterpart of the raw-data `106X_dataRun2_v35`).
- HLT step → `HLT:@relval2016` (simulation, not selection; needed
  for `TriggerResults` branch parity with data ALCARECO).
- Pileup handling → add MC-only `PileupSummaryInfos_*_*_*` to the
  ALCARECO output; PU-reweighting happens downstream in the
  histmaker, not at MC production time.
- Filter kinematic thresholds → keep pT(J/ψ)>6, pT(μ)>3.8,
  |η(μ)|<2.52. Compatible with the ALCARECO TkAlJpsiX loose muon
  selector (which cuts μ pT > 4 GeV).
- Anti-Λ_b decay table → mirror the McM particle-side decay
  block via `CDecay Myanti-Lambda_b0`. Physical BRs; no Λ→pπ
  forcing (Λ decay follows its natural ~64% pπ⁻ / 36% nπ⁰ split;
  the downstream `ALCARECOTkAlJpsiXLambdabResonances` producer
  reconstructs the pπ⁻ subset).
- Volume phasing → **phase-1 at 4×10⁸ generated events (20,000
  jobs)**. Sized to hit 1:1 MC:data on the priority channel
  B⁺→J/ψK⁺; Λ_b, K*⁰, KS, ψ(2S) all comfortably ≥ 2:1 at this
  scale. Phase-2 scale-up (+6×10⁸ for 2:1 on B⁺→J/ψK⁺) is a
  separate follow-up decision after phase-1 outputs are on disk
  and the first-pass calibration exposes any statistics-limited
  systematics.
- `base_seed` → `20260716`, recorded in
  `manifests/base_seed.txt` for reproducibility.
- Λ_b coverage → in-scope via inclusive with the anti-Λ_b fragment
  fix, no dedicated Λ_b production this iteration.
- Bs / Bc / dedicated B⁺ / dedicated Λ_b → out-of-scope this
  iteration.
- Data-yield anchor → final-analysis fitted yields (BPH-published),
  not ALCARECO-loose numbers. MC-vs-data ratio is conservative;
  ALCARECO-stage ratio will be larger than the nominal.
