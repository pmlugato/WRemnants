# Spec Delta: mc-production (ADDED)

## ADDED Requirements

### Requirement: Single CMSSW Release For Whole MC Chain

The MC inclusive B → J/ψ + X production SHALL be built and executed
end-to-end in `CMSSW_10_6_20_patch1` (the McM campaign release for
`BPH-RunIISummer20UL16GEN-00017`, which ships
`PythiaFilterMultiAncestor`; the raw-data ALCARECO's
`CMSSW_10_6_17_patch1` does not). No release-bridging between GEN
and RECO is permitted. Delta vs raw-data reco = 3 patch releases of
Pythia8/EvtGen bug fixes; no tracking/RECO/ALCARECO source changes
between the two, so calibration comparison stays apples-to-apples.

#### Scenario: Payload uses the campaign release

- **GIVEN** the McM campaign for the fragment is
  `RunIISummer20UL16GEN`
- **WHEN** `build_tarball.sh` runs
- **THEN** `mc_inclusive_payload.tgz` contains
  `CMSSW_10_6_20_patch1/src/` and no other CMSSW release area
- **AND** the packaged full-chain cfg's `--conditions` and cmsRun
  environment resolve inside `10_6_20_patch1`

#### Scenario: Wrong release refused

- **GIVEN** an operator attempts to build the payload with
  `SCRAM_ARCH` pointing at a release different from
  `CMSSW_10_6_20_patch1`
- **WHEN** `build_tarball.sh` runs
- **THEN** the build exits non-zero with a message identifying the
  release mismatch

### Requirement: Patched Fragment Provenance

The fragment used for GEN SHALL be
`fragments/BPH-RunIISummer20UL16GEN-00017-fragment.py` checked into
the launch directory, NOT the McM-served fragment. The patched
fragment MUST include `Alias`, `ChargeConj`, `CDecay`, and
`list_forced_decays` entries for `Myanti-Lambda_b0`, symmetric with
the treatment of `Myanti-Xi_b+` and `Myanti-Omega_b+` already present
in the McM-served version.

#### Scenario: Anti-Λ_b entries present

- **GIVEN** the patched fragment
- **WHEN** grep is run for `anti-Lambda_b0`
- **THEN** at least one `Alias Myanti-Lambda_b0`, one
  `ChargeConj MyLambda_b0 Myanti-Lambda_b0`, one
  `CDecay Myanti-Lambda_b0`, and one entry
  `'Myanti-Lambda_b0'` in `list_forced_decays` are found

#### Scenario: Configuration annotation is not stale

- **GIVEN** the patched fragment
- **WHEN** grep is run for
  `'Jpsi->mumu (no  kin cuts on muons)'`
- **THEN** no match is found

#### Scenario: McM curl bypass

- **GIVEN** `build_tarball.sh` prepares the payload
- **WHEN** the fragment placement step runs
- **THEN** `cp` from `fragments/…` is used, NOT
  `curl` from `cms-pdmv-prod.web.cern.ch`
- **AND** `request_fragment_check.py` against McM is NOT invoked
  during the build

### Requirement: Full-Chain Single cmsRun Per Job

Each condor job SHALL run
`GEN,SIM,DIGI,L1,DIGI2RAW,HLT,RAW2DIGI,L1Reco,RECO,PAT,ALCA:TkAlJpsiX_step`
in a single cmsRun invocation and emit exactly one ALCARECO ROOT
file per job.

#### Scenario: One cmsRun, one output

- **GIVEN** a successfully completed job
- **WHEN** the worker finalizes
- **THEN** exactly one ALCARECO ROOT file is xrdcp'd to ceph
- **AND** the per-job JSON records `exit_code: 0` and `exit_reason:
  ok`

#### Scenario: cmsRun crash before ALCA step

- **GIVEN** cmsRun crashes mid-chain (e.g. during DIGI premix)
- **WHEN** the wrapper finalizes
- **THEN** no partial ROOT file is xrdcp'd to ceph
- **AND** the per-job JSON records
  `exit_reason: cmsrun_crash` with a stderr tail

### Requirement: ALCARECO Content Parity With Raw-Data ALCARECO

The MC ALCARECO output SHALL contain the same set of branch names as
the raw-data ALCARECO on the `ALCARECOTkAlJpsiX*` prefix, at
minimum: `ALCARECOTkAlJpsiX`, `ALCARECOTkAlJpsiXDeDxHarmonic2`,
`ALCARECOTkAlJpsiXDeDxPixelHarmonic2`,
`ALCARECOTkAlJpsiXDeDxAllHarmonic2`, and all seven
`ALCARECOTkAlJpsiX{BPlus,B0Kstar,B0Ks,BsPhi,Lambdab,Psi2S,Bc}Resonances`
collections. The MC ALCARECO output MAY additionally contain the
MC-only `PileupSummaryInfos_*_*_*` branch for downstream PU
reweighting; this is the only permitted MC-vs-data additive
difference.

#### Scenario: All branches present on MC

- **GIVEN** an MC ALCARECO output ROOT file
- **WHEN** `edmDumpEventContent` is run
- **THEN** every branch name in
  `Alignment/CommonAlignmentProducer/python/ALCARECOTkAlJpsiX_Output_cff.py`
  is listed exactly once
- **AND** the `PileupSummaryInfos_*_*_*` branch is also listed
- **AND** no `ALCARECOTkAlJpsiX*Resonances` collection has zero
  candidates on B⁺/B⁰ channels aggregated across the rehearsal
  batch

#### Scenario: Data ALCARECO byte-set match on shared prefix

- **GIVEN** a raw-data ALCARECO output from
  `condor/full_alcareco_2016postvfp/` and an MC ALCARECO output
  from this launch
- **WHEN** the two `edmDumpEventContent` outputs are diffed on the
  `ALCARECOTkAlJpsiX*` prefix
- **THEN** the diff is empty
- **AND** the only non-empty additive diff (MC minus data) is on
  the `PileupSummaryInfos_*_*_*` line

### Requirement: Deterministic Per-Job Seed

The worker SHALL inject an initial random seed derived from
`(clusterId, procId)` into every cfg it runs
(`initialSeed = base_seed + clusterId × 100000 + procId`), and MUST
record the seed used in the per-job JSON.

#### Scenario: Seed appears in JSON

- **GIVEN** a completed job with `clusterId=12345` and `procId=42`
- **WHEN** the wrapper finalizes
- **THEN** the per-job JSON `seed_used` field equals
  `base_seed + 12345 × 100000 + 42`

#### Scenario: Reproducibility on rerun

- **GIVEN** a job with a given `(clusterId, procId)` that produced
  ROOT output `O`
- **WHEN** the same `(clusterId, procId)` is re-executed
- **THEN** the new ROOT output has the same generated-event count
  and the same accepted-event count as `O`

#### Scenario: Resubmit collision safety

- **GIVEN** cluster A used `procId ∈ [0, 1000)` and cluster B is a
  resubmit of the same batch
- **WHEN** cluster B assigns its own `clusterId`
- **THEN** cluster B's `seed_used` values do not overlap cluster A's

### Requirement: Per-Job JSON Manifest With Structured Exit

Every job SHALL write a per-job JSON summary containing at minimum:
`cluster_id`, `proc_id`, `seed_used`, `n_generated`, `n_accepted`,
`filter_efficiency`, `output_file`, `output_size_mb`,
`runtime_s_cmsrun`, `peak_rss_mb`, `exit_code`, `exit_reason`,
`host`. The JSON MUST be xrdcp'd to the same ceph destination
directory as the ROOT output.

#### Scenario: Successful job JSON

- **GIVEN** a completed job with cmsRun exit 0 and successful xrdcp
- **WHEN** the wrapper finalizes
- **THEN** the JSON `exit_reason` field equals `ok`
- **AND** all required fields are populated with non-null values
  except `n_accepted` may be zero if the filter rejects everything

#### Scenario: xrdcp of output failed

- **GIVEN** cmsRun exited zero but the output xrdcp to ceph failed
  on all attempts
- **WHEN** the wrapper finalizes
- **THEN** the JSON `exit_reason` field equals `xrdcp_out_failed`
- **AND** the local ROOT file has been removed from scratch

### Requirement: Wall-Time Watchdog (10-Hour Cap)

The worker SHALL kill cmsRun with SIGTERM if cmsRun wall-time exceeds
a flat 10-hour cap and MUST record `exit_reason: wrapper_timeout`
in the per-job JSON. The cap is chosen to sit well above the
projected ~7 h per-job runtime and exists only to catch a hung
process.

#### Scenario: cmsRun spins past the cap

- **GIVEN** a job whose cmsRun has been running for 36000 s
- **WHEN** the wall-time watchdog fires
- **THEN** cmsRun receives SIGTERM
- **AND** the per-job JSON `exit_reason` equals `wrapper_timeout`

### Requirement: Manual Reconcile-and-Resubmit Contract

`find_missing.py` SHALL classify every scheduled
`(clusterId, procId)` pair into one of {`success`, `unrun`,
`failed_cmsrun`, `failed_xrdcp`, `failed_wrapper`}, and MUST emit
an auto-resubmit list containing ONLY
`unrun ∪ failed_xrdcp ∪ failed_wrapper`.

#### Scenario: Mixed-outcome cluster reconcile

- **GIVEN** a cluster's manifest contains 80 successes, 12
  failed_xrdcp, 3 failed_cmsrun, 5 unrun
- **WHEN** `find_missing.py` runs on the cluster
- **THEN** the resubmit list contains the 12 + 5 = 17 auto-safe
  entries
- **AND** the 3 `failed_cmsrun` entries are NOT included in the
  resubmit list
- **AND** the 3 `failed_cmsrun` entries are reported to stdout for
  human triage

#### Scenario: Resume filter prevents double compute

- **GIVEN** `manifests/success.txt` records 40 successful
  `(clusterId, procId)` pairs
- **WHEN** `submit_batch.sh` is invoked with the same batch size
- **THEN** the rendered queue file excludes those 40 pairs

### Requirement: Output Layout Under Ceph

MC ALCARECO output SHALL be written under
`/data/user/p/pmlugato/mz/mc/inclusive_btojpsix_2016postvfp/` via
`root://submit50.mit.edu/` xrdcp, and output ROOT filenames MUST
include `clusterId` and `procId` to prevent collision across jobs.

#### Scenario: Two workers finalize concurrently

- **GIVEN** two workers, `(cluster=1, proc=7)` and
  `(cluster=1, proc=8)`, both writing to the same ceph directory
- **WHEN** both workers finalize
- **THEN** the ceph directory contains two ROOT files whose names
  each embed the respective `(clusterId, procId)` pair (no
  collision)

### Requirement: Preflight Ceph Quota Check

`quota.sh` SHALL be executed before any batch of ≥ 500 jobs, and
the batch MUST NOT be submitted if the free space at
`/data/user/p/pmlugato/` is less than the projected batch output
volume plus 500 GB headroom.

#### Scenario: Insufficient headroom

- **GIVEN** free space is 200 GB and the projected batch output is
  1000 GB
- **WHEN** the operator runs the preflight check
- **THEN** `quota.sh` exits non-zero
- **AND** prints the required-vs-available numbers

### Requirement: Local Pre-Handoff Validation

A local multi-parallel run targeting ~1000 accepted events SHALL be
executed on this workstation before the handoff tarball is shipped,
and its measurements MUST be recorded in
`local_validation_report.md`. The report MUST contain: p99 per-job
scratch usage, p99 peak RSS, mean wall-clock per accepted event, the
byte-diff result against a raw-data ALCARECO file on the
`ALCARECOTkAlJpsiX*` branch prefix, the seed-determinism result, and
the Λ_b species fraction from `classify_jpsi.py` (target ~4–5% with
the anti-Λ_b patch, up from ~2.3% without).

#### Scenario: Local run completes and report exists

- **GIVEN** the local multi-parallel driver has finished
- **WHEN** the operator inspects the launch directory
- **THEN** `local_validation_report.md` exists and lists all six
  measurements
- **AND** `RequestDisk` and `RequestMemory` in the shipped `.sub`
  files match the report's chosen values

#### Scenario: Branch-content parity fails locally

- **GIVEN** the local ALCARECO output byte-diffs against the
  raw-data reference on the `ALCARECOTkAlJpsiX*` prefix
- **WHEN** the diff is non-empty
- **THEN** the handoff tarball MUST NOT be built
- **AND** the failing branches are logged in
  `local_validation_report.md` for triage

### Requirement: Self-Contained Operator Handoff Package

The change SHALL produce a single tarball
`mc_inclusive_handoff_<TS>.tgz` containing every artifact the
operator needs to run the production without further input from
this side. At minimum: the CMSSW payload tarball, all `.sub` files
with tuned `RequestDisk`/`RequestMemory`, `run.sh` worker,
`submit_batch.sh` driver, `find_missing.py` reconcile,
`build_manifest.py` and `merge_manifests.py` roll-ups, `SUBMIT.sh`
and `RECONCILE.sh` entry-point scripts, `quota.sh` preflight, and
a `HANDOFF_README.md` runbook.

#### Scenario: Operator can execute rehearsal from the tarball alone

- **GIVEN** the handoff tarball as the only input on the
  operator's submit host
- **WHEN** the operator unpacks the tarball and runs
  `./SUBMIT.sh rehearsal`
- **THEN** the rehearsal cluster is submitted with no missing
  file or command
- **AND** no reference to paths under our workstation appears in
  any operator-executed script

#### Scenario: Rehearsal reconcile is our sign-off gate

- **GIVEN** the operator has drained the 100-job rehearsal and run
  `./RECONCILE.sh`
- **WHEN** the operator sends the reconcile output back to us
- **THEN** we evaluate (i) success rate ≥ 95% after one resubmit
  cycle, (ii) sample-output branch content matches our local
  reference, (iii) no unexpected `failed_cmsrun` patterns
- **AND** only after a positive sign-off from us does the
  operator run `./SUBMIT.sh full`
