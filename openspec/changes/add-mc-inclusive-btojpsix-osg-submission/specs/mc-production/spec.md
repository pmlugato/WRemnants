# Spec Delta: mc-production (ADDED — OSG rehearsal path)

Layers on top of the requirements introduced by the parent change
[`add-mc-inclusive-btojpsix-2016postvfp`](../../../add-mc-inclusive-btojpsix-2016postvfp/specs/mc-production/spec.md).
All parent requirements (patched fragment, seed determinism, ALCARECO
content parity, per-job JSON, 10-hour watchdog, reconcile buckets,
ceph output layout, self-contained handoff) continue to hold for jobs
routed to OSG; the additions below constrain the routing-specific
mechanics.

## ADDED Requirements

### Requirement: OSG Rehearsal Launch Directory Reuses Parent Artifacts

The OSG launch directory SHALL obtain its physics configuration by
symlink from the parent directory
`condor/mc_inclusive_btojpsix_2016postvfp/`. Reused artifacts MUST
include `run.sh`, `mc_inclusive_payload.tgz`,
`mc_inclusive_btojpsix_fullchain.py`, `find_missing.py`,
`build_manifest.py`, `merge_manifests.py`, `quota.sh`,
`manifests/base_seed.txt`, and `manifests/reco_tags.txt`. The OSG
directory MUST NOT copy or fork these files.

#### Scenario: Physics config is byte-identical to parent

- **GIVEN** the OSG launch directory has been scaffolded
- **WHEN** `md5sum` is run on `mc_inclusive_btojpsix_fullchain.py`,
  `run.sh`, `mc_inclusive_payload.tgz`, and
  `manifests/base_seed.txt` in both directories
- **THEN** the checksums match on every file (via symlink)

#### Scenario: Parent bugfix propagates automatically

- **GIVEN** a bugfix is applied to `run.sh` in the parent directory
- **WHEN** the OSG directory is inspected
- **THEN** `readlink -f run.sh` in the OSG directory resolves to
  the parent's file and the bugfix is present without a second edit

### Requirement: OSG-Specific Submit Attributes

The OSG submit template `submits/base_mc_inclusive_osg.sub` SHALL
declare, at minimum:

- `+ProjectName = "MIT_submit"`
- `+JobDurationCategory = "Medium"`
- `+SingularityImage = "/cvmfs/singularity.opensciencegrid.org/cmssw/cms:rhel7"`
- `+SingularityBindCVMFS = True`

and MUST NOT declare `+AccountingGroup` (which routes to the CMS
global pool rather than OSPool).

#### Scenario: OSG routing attributes present

- **GIVEN** `submits/base_mc_inclusive_osg.sub`
- **WHEN** `grep` is run for each of the four required attributes
- **THEN** each is present exactly once with the value above

#### Scenario: CMS-pool accounting attribute absent

- **GIVEN** `submits/base_mc_inclusive_osg.sub`
- **WHEN** `grep -c AccountingGroup` is run
- **THEN** the count is zero

### Requirement: OSG Requirements Clause

The `Requirements` clause SHALL contain, at minimum, the four ANDed
conditions:

- `OSGVO_OS_STRING == "RHEL 7"` (SLC7 for CMSSW_10_6_20_patch1)
- `HAS_SINGULARITY == TRUE`
- `HAS_CVMFS_singularity_opensciencegrid_org == True`
- `HAS_CVMFS_cms_cern_ch == True`

and MUST NOT contain the parent's `BOSCOCluster =!= …` exclusions
(which target MIT-T3 CEs on the CMS global pool, not relevant to
OSPool matching).

#### Scenario: All four OSG requirement clauses present

- **GIVEN** the OSG submit file
- **WHEN** the `Requirements` line is parsed
- **THEN** the four required conjuncts above are all present

#### Scenario: CMS-pool exclusion clauses absent

- **GIVEN** the OSG submit file
- **WHEN** `grep BOSCOCluster` is run
- **THEN** no line matches

### Requirement: Rehearsal-Only Scope

The OSG launch directory SHALL NOT support a phase-1 target.
`SUBMIT_OSG.sh` MUST accept `rehearsal` (100 jobs) as its sole
mode and MUST exit non-zero for any other argument, including
`phase1` or `full`. The rationale is that phase-1 routing (OSG
vs MIT T3 vs split) is a decision made after both rehearsal
reconciles have landed and is executed from the parent directory,
not from the OSG rehearsal directory.

#### Scenario: Rehearsal mode accepted

- **GIVEN** a valid CMS proxy on the submit host
- **WHEN** `./SUBMIT_OSG.sh rehearsal` is invoked
- **THEN** a 100-job cluster is submitted and its ID appended
  to `manifests/osg_rehearsal_clusters.txt`

#### Scenario: Phase-1 mode refused

- **GIVEN** a valid CMS proxy
- **WHEN** `./SUBMIT_OSG.sh phase1` is invoked
- **THEN** the script exits non-zero without submitting
- **AND** stderr explains that phase-1 launches from the parent
  directory after retune

### Requirement: OSG Output Isolation Under Ceph Subdir

MC ALCARECO output from the OSG rehearsal SHALL be written under
`/data/user/p/pmlugato/mz/mc/inclusive_btojpsix_2016postvfp/osg_rehearsal/`
(a subdirectory of the parent's ceph root), and MUST NOT be written
to the parent root directly. The subdirectory MUST be created via
`xrdfs submit50.mit.edu mkdir` before the first OSG job is queued
if it does not already exist.

#### Scenario: OSG output lands in the subdir

- **GIVEN** a successfully completed OSG rehearsal job
- **WHEN** the ceph destination directory is listed
- **THEN** the output ROOT and per-job JSON appear under
  `.../inclusive_btojpsix_2016postvfp/osg_rehearsal/` and not
  under `.../inclusive_btojpsix_2016postvfp/` directly

#### Scenario: OSG reconcile does not touch parent subdir

- **GIVEN** David's MIT-T3 rehearsal has produced outputs under
  the parent's ceph root and the OSG rehearsal has produced
  outputs under `osg_rehearsal/`
- **WHEN** `./RECONCILE_OSG.sh` is run
- **THEN** the reconcile emits `(clusterId, procId)` pairs only
  from OSG clusters recorded in
  `manifests/osg_rehearsal_clusters.txt`
- **AND** does not list, delete, or resubmit any of David's
  cluster/proc pairs

### Requirement: Placeholder Resource Requests Match Parent Rehearsal

The OSG rehearsal `.sub` SHALL ship with resource requests matching
the parent's shipped MIT `.sub` (`RequestMemory=8000`,
`RequestDisk=16 GiB`, `RequestCpus=4`, and N_GEN argument = 800)
so the two rehearsals produce directly comparable wall/RSS/scratch
distributions. Retune of these four numbers is out of scope for
this change and MUST NOT be applied inside the OSG launch
directory; retune happens after both rehearsals report and lands
in the parent directory (or a follow-up change).

#### Scenario: OSG placeholder equals parent

- **GIVEN** the parent `.sub` at
  `condor/mc_inclusive_btojpsix_2016postvfp/submits/base_mc_inclusive.sub`
  and the OSG `.sub` at
  `condor/mc_inclusive_btojpsix_2016postvfp_osg/submits/base_mc_inclusive_osg.sub`
- **WHEN** the two files' `RequestMemory`, `RequestDisk`,
  `RequestCpus`, and N_GEN argument values are extracted
- **THEN** all four values match

#### Scenario: In-rehearsal memory drop is a documented deviation

- **GIVEN** the OSG rehearsal is stalled in matching because of
  a narrow-slot RequestMemory penalty (Open Question §1 in the
  proposal)
- **WHEN** the operator drops `RequestMemory` to 4000 and
  resubmits the stalled portion
- **THEN** the deviation is recorded in `osg_rehearsal_report.md`
  with the observed matching latency before and after

### Requirement: Two-Cycle Rehearsal Success Rate Sign-Off

The OSG rehearsal SHALL be considered a passing routing test only
if the two-cycle (initial submission + one auto-safe resubmit)
success rate is ≥ 95%. A failing rehearsal blocks the OSG option
for phase-1 routing until the failure taxonomy is triaged.

#### Scenario: Rehearsal passes

- **GIVEN** 100 jobs submitted, 88 successful on first drain,
  9 auto-safe resubmitted and all succeeding, 3 `failed_cmsrun`
  held for triage
- **WHEN** the two-cycle success rate is computed
- **THEN** the rate is 97/100 = 97%, ≥ 95%
- **AND** `osg_rehearsal_report.md` records the pass and the
  three held cases are documented

#### Scenario: Rehearsal fails

- **GIVEN** 100 jobs submitted with a 60% two-cycle success rate
- **WHEN** the report is compiled
- **THEN** the report flags OSG as not-yet-ready for phase-1
  routing
- **AND** the parent's phase-1 decision defaults to David's MIT-T3
  route unless the OSG failure taxonomy identifies a fixable
  cause

### Requirement: OSG Rehearsal Report

An `osg_rehearsal_report.md` SHALL be produced under
`condor/mc_inclusive_btojpsix_2016postvfp_osg/` at the end of the
rehearsal. The report MUST document, at minimum: submission
cluster IDs, matching latency (idle-to-run p50 and p99),
two-cycle success rate, p99 wall/RSS/scratch across the 100
jobs, the reconcile failure taxonomy, and recommended `.sub`
deltas for phase-1 sizing based on the measurements.

#### Scenario: Report exists and is complete

- **GIVEN** the OSG rehearsal has drained and reconciled
- **WHEN** the launch directory is inspected
- **THEN** `osg_rehearsal_report.md` exists and lists all six
  required items above
- **AND** the report is the input to the phase-1 sizing decision
  (jointly with David's MIT-T3 reconcile output)
