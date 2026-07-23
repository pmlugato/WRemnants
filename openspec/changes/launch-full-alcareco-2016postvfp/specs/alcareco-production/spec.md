# Spec Delta: alcareco-production (ADDED)

## ADDED Requirements

### Requirement: Era Boundary from Recorded Convention

The full 2016 postVFP launch SHALL define the era membership from a
recorded convention file
(`condor/full_alcareco_2016postvfp/lfns/postvfp_run_boundary.txt`),
NOT from a hard-coded value in a script.

#### Scenario: Boundary file exists and is used

- **GIVEN** the boundary file contains a run number `R` with a
  provenance comment
- **WHEN** `submit_batch.sh charmonium 2016F_postVFP` is invoked
- **THEN** the filelist consumed is
  `lfns/charmonium_2016F_postVFP_populated.txt`
- **AND** that filelist contains only files whose runs are ≥ `R`

#### Scenario: Boundary file missing

- **GIVEN** the boundary file is absent
- **WHEN** any launch step that needs the era boundary runs
- **THEN** the step exits non-zero with a message pointing at the
  expected file path

### Requirement: Filelist Provenance via DAS with Populated Filter

Every launch filelist SHALL be derived from a `dasgoclient` query
whose raw TSV output is cached under `lfns/raw/`, and MUST be
filtered to files with `nevents ≥ 1000` before being fed to the
submit driver.

#### Scenario: New era filelist is built

- **GIVEN** a target `(pd, era)` with no cached filelist
- **WHEN** the filelist-generation step runs
- **THEN** `lfns/raw/${pd}_${era}_<TS>.tsv` is written with the
  full DAS output including per-file event counts
- **AND** `lfns/${pd}_${era}_populated.txt` is written containing
  only LFNs with `nevents ≥ 1000`

#### Scenario: 2016H validation cache is reused

- **GIVEN** validation already produced 2016H populated filelists
- **WHEN** the launch prepares its 2016H filelist
- **THEN** the launch reuses (symlink or copy) the validation
  cache and does NOT re-issue the DAS query

### Requirement: Worker Stage-in Via xrdcp With AAA Fallback

The worker script SHALL stage the input LFN to
`$_CONDOR_SCRATCH_DIR/<basename>.root` via `xrdcp` before invoking
cmsRun, and cmsRun MUST read the input as `file:./<basename>.root`
— streaming reads from within cmsRun are prohibited. On xrdcp
failure the wrapper SHALL retry with a secondary redirector, then
a tertiary redirector. An end-of-job trap MUST `rm -f` the local
copy on both success and failure paths.

#### Scenario: Global redirector fails, US regional succeeds

- **GIVEN** an input LFN and `root://cms-xrd-global.cern.ch/` as
  the primary redirector
- **WHEN** the first xrdcp attempt exits non-zero
- **THEN** the wrapper retries xrdcp with
  `root://cmsxrootd.fnal.gov/`
- **AND** on success, cmsRun is invoked reading
  `file:./<basename>.root`
- **AND** the job completes successfully

#### Scenario: All three redirectors fail

- **GIVEN** an input LFN
- **WHEN** xrdcp exits non-zero on the primary, secondary, and
  tertiary redirectors in sequence
- **THEN** the wrapper writes `exit_reason: input_stagein_failed`
  into the per-file JSON and exits non-zero
- **AND** the job is classified as `failed_input_open` by
  `find_missing.py`

#### Scenario: Local copy is cleaned up on wrapper failure

- **GIVEN** stage-in succeeded and cmsRun subsequently crashed
- **WHEN** the wrapper exits
- **THEN** `$_CONDOR_SCRATCH_DIR/<basename>.root` no longer exists

### Requirement: Per-File JSON Summary With Structured Exit Reason

Every job SHALL write a per-file JSON summary containing at
minimum: `recoskim`, `input_lfn`, `output_files`,
`output_size_mb_total`, `runtime_s_cmsrun`, `peak_rss_mb`,
`exit_code`, `exit_reason`, `host`, `run_min`, `run_max`,
`lumi_count`, and this JSON MUST be xrdcp'd to the same ceph
destination directory as the ROOT outputs.

#### Scenario: Successful job JSON

- **GIVEN** a completed job with cmsRun exit 0 and successful
  xrdcp
- **WHEN** the wrapper finalizes
- **THEN** the JSON `exit_reason` field equals `ok`
- **AND** all required fields are populated

#### Scenario: cmsRun crash JSON

- **GIVEN** cmsRun exits with a non-zero code other than 66
- **WHEN** the wrapper finalizes
- **THEN** the JSON `exit_reason` field equals `cmsrun_crash`
- **AND** the full cmsRun stderr tail is included in the JSON
  (last 40 lines)

#### Scenario: Run/lumi extraction unavailable

- **GIVEN** `edmFileUtil -f -e` cannot open the input file at
  wrapper start
- **WHEN** the wrapper proceeds
- **THEN** `run_min`, `run_max`, `lumi_count` are set to `null`
- **AND** the wrapper does NOT fail

### Requirement: Wall-Time Watchdog (Flat 2-Hour Cap)

The worker SHALL kill cmsRun with SIGTERM if cmsRun wall-time
exceeds a flat 2-hour cap, applied uniformly to both PDs, and MUST
record `exit_reason: wrapper_timeout` in the per-file JSON. The
cap is chosen to sit well above the observed p99 of both PDs so
that legitimate large-file runs are not truncated; it exists only
to catch a truly hung process.

#### Scenario: cmsRun spins past the cap

- **GIVEN** a job whose cmsRun has been running for 7200 s
- **WHEN** the wall-time watchdog fires
- **THEN** cmsRun receives SIGTERM
- **AND** the per-file JSON `exit_reason` equals `wrapper_timeout`

#### Scenario: Legitimately long file completes under the cap

- **GIVEN** a Charmonium file that takes 1800 s (well above the
  mean but under the cap)
- **WHEN** cmsRun completes normally
- **THEN** the watchdog does NOT fire
- **AND** the per-file JSON `exit_reason` equals `ok`

### Requirement: Manual Reconcile-and-Resubmit Contract

`find_missing.py` SHALL classify every dispatched input into one
of {`success`, `unrun`, `failed_input_open`, `failed_cmsrun`,
`failed_xrdcp`}, and MUST emit an auto-resubmit list containing
ONLY `unrun`, `failed_input_open`, and `failed_xrdcp`.

#### Scenario: Mixed-outcome cluster reconcile

- **GIVEN** a cluster's manifest contains 80 successes, 12
  failed_input_open, 3 failed_cmsrun, 5 unrun
- **WHEN** `find_missing.py --pd charmonium --era 2016G` runs
- **THEN** the resubmit list contains the 12 + 5 = 17 auto-safe
  entries
- **AND** the 3 `failed_cmsrun` entries are NOT included in the
  resubmit list
- **AND** the 3 `failed_cmsrun` entries are reported to stdout for
  human triage

#### Scenario: Resume filter prevents double compute

- **GIVEN** `manifests/charmonium_2016G_success.txt` records 40
  successful LFNs
- **WHEN** `submit_batch.sh charmonium 2016G` is invoked with the
  full populated filelist
- **THEN** the queue file passed to `condor_submit` excludes those
  40 LFNs

### Requirement: Per-Cluster Manifest at Close-Out

`build_manifest.py` SHALL produce exactly one manifest JSON per
Condor cluster ID, and this JSON MUST contain per-input records
with the fields: `input_lfn`, `output_files_on_ceph`, `exit_code`,
`exit_reason`, `run_min`, `run_max`, `runtime_s_cmsrun`,
`peak_rss_mb`, `host`.

#### Scenario: Cluster drains and manifest is built

- **GIVEN** a Condor cluster with 200 jobs, all drained
- **WHEN** `build_manifest.py <clusterId>` runs
- **THEN** `manifests/<pd>_<era>_cluster_<clusterId>.json` contains
  exactly 200 records
- **AND** the per-input JSON fields on ceph and the condor log
  records are joined into each record

### Requirement: Site Deny List Interpolated at Submit Time

Submit files SHALL derive `DESIRED_Sites` by subtracting entries in
`site_deny.txt` from the base allowed site list at
`condor_submit` time, and MUST NOT hard-code the offender list into
`.sub` files.

#### Scenario: Site added to deny list

- **GIVEN** `site_deny.txt` gains a new entry `T2_XX_Example`
- **WHEN** the next `submit_batch.sh` invocation runs
- **THEN** the submitted cluster's `+DESIRED_Sites` classad does
  NOT contain `T2_XX_Example`

### Requirement: Output Layout Under Ceph

Per-PD per-era outputs SHALL be written under
`/data/user/p/pmlugato/mz/alcareco/full_2016postvfp/${pd}/Run2016${era}/`
via `root://submit50.mit.edu/` xrdcp, and output ROOT filenames
MUST include the input LFN basename to prevent collision.

#### Scenario: Two workers process files in the same PD+era

- **GIVEN** two workers, one processing file `A.root` and the
  other `B.root`, both writing to
  `/data/user/p/pmlugato/mz/alcareco/full_2016postvfp/charmonium/Run2016G/`
- **WHEN** both workers finalize
- **THEN** the ceph directory contains outputs named with `A` and
  `B` embedded (one output per input, no collision)

### Requirement: Preflight Ceph Quota Check

`quota.sh` SHALL be executed before any batch of ≥ 500 jobs, and
the batch MUST NOT be submitted if the free space at
`/data/user/p/pmlugato/` is less than 2× the projected batch
output volume.

#### Scenario: Insufficient headroom

- **GIVEN** free space is 100 GB and the projected batch output is
  200 GB
- **WHEN** the operator runs the preflight check
- **THEN** `quota.sh` exits non-zero
- **AND** prints the required-vs-available numbers
