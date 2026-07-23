# Change: Launch full 2016 postVFP AlCaReco production across Charmonium and SingleMuon PDs

## Why

`add-preprod-multichannel-validation` closed with all four AlCaReco
streams validated on 2016H at 100 files/PD. Projections came in at
≲1 wall day @ 200 slots for the combined 2016H workload
(358 GB / 2758 CPU-hrs). We are ready to launch the real production
over the full 2016 postVFP era. The validation directory
(`condor/multichannel_alcareco/`) is scoped to a smoke workload and
lacks the pieces we need for a launch that spans ~44k Charmonium and
~130k+ SingleMuon files with realistic failure handling.

This proposal stands up a dedicated launch directory
(`condor/full_alcareco_2016postvfp/`) and defines the operational
contract for filelist provenance, per-PD outputs, retry/reconcile,
and post-batch bookkeeping. No physics/producer change; no cff
change; no C++ change.

## What Changes

### Scope (locked in with the user)

- **Era**: 2016 postVFP = (2016F ∩ runs ≥ HIP-fix boundary) ∪ 2016G
  ∪ 2016H. The exact boundary run must be pulled from the WRemnants
  postVFP convention, not guessed (task 0.1).
- **PDs**: `/Charmonium/Run2016{F,G,H}-v1/RAW` and
  `/SingleMuon/Run2016{F,G,H}-v1/RAW`.
- **Streams per PD** (unchanged from validation):
  - Charmonium: `ALCARECOTkAlJpsiX` + `ALCARECOTkAlJpsiMuMu`.
  - SingleMuon: `ALCARECOTkAlDstToD0Pi` + `ALCARECOTkAlKShortTracks`
    + `ALCARECOTkAlLambdaTracks`.
- **Output**:
  `/data/user/p/pmlugato/mz/alcareco/full_2016postvfp/{charmonium,singlemuon}/Run2016{F,G,H}/`,
  written via `root://submit50.mit.edu/` xrdcp (same host as
  validation).
- **Read policy**: **stage-in, no streaming**. Worker `xrdcp`s the
  input to local scratch (`$_CONDOR_SCRATCH_DIR`), cmsRun reads
  `file:./<basename>.root`, wrapper deletes the local copy at
  end-of-job. AAA fallback (`root://cms-xrd-global.cern.ch/`,
  `root://cmsxrootd.fnal.gov/`, `root://xrootd-cms.infn.it/`) is
  applied at the xrdcp stage, not at cmsRun stage — cleaner error
  surface, no mid-file stream hangs.
- **Retry**: manual reconcile-and-resubmit via
  `find_missing.py`. No DAGman, no HTCondor `max_retries`.

### New launch directory `condor/full_alcareco_2016postvfp/`

Structured for a two-PD × three-era launch with resubmit hygiene.

```
condor/full_alcareco_2016postvfp/
├── build_tarball.sh                 # bundles CMSSW + both recoskim cfgs
├── make_recoskim_cfgs.sh            # cmsDriver wrappers (per-era ok)
├── run.sh                           # worker (AAA fallback + per-file JSON)
├── submit_batch.sh                  # wraps condor_submit over one PD+era
├── find_missing.py                  # reconcile submitted vs successful
├── build_manifest.py                # per-batch manifest (input, output, exit, run/lumi)
├── merge_manifests.py               # roll up per-batch manifests -> era-wide + PD-wide
├── lfns/
│   ├── charmonium_2016F_postVFP_populated.txt
│   ├── charmonium_2016G_populated.txt
│   ├── charmonium_2016H_populated.txt   (reuses validation cache)
│   ├── singlemuon_2016F_postVFP_populated.txt
│   ├── singlemuon_2016G_populated.txt
│   ├── singlemuon_2016H_populated.txt   (reuses validation cache)
│   └── raw/                              # unfiltered DAS output kept for provenance
├── submits/
│   ├── charmonium_Run2016F_postVFP.sub
│   ├── charmonium_Run2016G.sub
│   ├── charmonium_Run2016H.sub
│   ├── singlemuon_Run2016F_postVFP.sub
│   ├── singlemuon_Run2016G.sub
│   └── singlemuon_Run2016H.sub
├── manifests/                            # one JSON per submit cluster (auto-created)
├── site_deny.txt                         # sites that emitted EX_NOINPUT in validation
└── README.md                             # runbook (submit, watch, resubmit, close-out)
```

Rationale for splitting by era in the submit layer: cluster IDs stay
scoped, ceph output dir is naturally partitioned by era, and a
resubmit of "SingleMuon 2016G" doesn't need to touch F or H.

### Filelist provenance and populated-file filter

- Rebuild filelists from DAS with `dasgoclient` — the EOS namespace at
  CERN does not hold 2016 RAW replicas (Rucio-managed, distributed
  across Tier-2s). `ls /eos/cms/store/data/Run2016H/...` at CERN
  won't be authoritative and may miss populated files that only live
  at US Tier-2s.
- For 2016F we additionally apply a run-range filter with
  `dasgoclient --query "run,dataset=/Charmonium/Run2016F-v1/RAW"`
  and drop runs below the postVFP boundary before writing
  `charmonium_2016F_postVFP.txt`.
- Populated-file filter: DAS `file dataset=... | grep file.name,
  file.nevents | awk '$2 >= 1000'`. 2016H had ~11% zero-event
  placeholders; same filter applied uniformly across F/G/H.
- All DAS output is cached under `lfns/raw/` with a timestamp so we
  can re-derive the filtered list without re-querying.

### Worker script (`run.sh`) changes vs validation

Additive changes only; body of the recoskim invocation is unchanged.

- **Stage-in via xrdcp with AAA fallback**: before calling `cmsRun`,
  the wrapper xrdcp's the input LFN to
  `$_CONDOR_SCRATCH_DIR/<basename>.root` using
  `root://cms-xrd-global.cern.ch/` first; on non-zero exit it retries
  with `root://cmsxrootd.fnal.gov/`, then
  `root://xrootd-cms.infn.it/`. cmsRun then reads from
  `file:./<basename>.root` — no xrootd inside cmsRun. Worker
  `rm -f` the local copy in a trap at end-of-job so failures still
  clean up. Requires bumping `RequestDisk` to accommodate the
  largest RAW (~5 GB) plus the two output ROOT files.
- **Per-file JSON extended**: adds `run_min`, `run_max`, `lumi_count`
  (from `edmFileUtil -f -e <local-file>` after stage-in, best-effort;
  falls back to `null` on error). Enables run/lumi-level bookkeeping
  downstream.
- **Structured exit tags**: `exit_reason` field in JSON is one of
  `ok | input_stagein_failed | cmsrun_crash | xrdcp_out_failed |
  wrapper_error | wrapper_timeout`. Makes `find_missing.py`
  classification trivial.
- **Wall-time watchdog (generous)**: soft kill at a flat **2-hour**
  cap regardless of PD. Well above the observed p99 for both
  Charmonium (~1200 s) and SingleMuon (~300 s) with headroom for
  larger-than-average RAW files. Purpose is only to catch a truly
  hung process (network stall in xrdcp of the output, dead cmsRun
  spin), not to police per-file runtime variance. Records
  `exit_reason: wrapper_timeout` in the JSON.

### Submit driver (`submit_batch.sh`)

One entry point: `./submit_batch.sh <pd> <era> [--dry-run] [--limit N]`.
Reads the corresponding `lfns/${pd}_${era}_populated.txt`, applies
the resume filter (skip inputs already present in
`manifests/${pd}_${era}_success.txt`), writes a per-cluster
`queue_${pd}_${era}_${TIMESTAMP}.txt`, and calls `condor_submit` on
the matching `.sub`. Records the cluster ID into
`manifests/${pd}_${era}_clusters.log`.

### Failure handling: `find_missing.py` + resubmit

- `find_missing.py --pd charmonium --era 2016G` walks
  `manifests/*.json` for that (pd, era), joins against
  `lfns/charmonium_2016G_populated.txt`, prints:
  - `unrun`   — input never had a job dispatched
  - `failed_input_open` — dispatched, EX_NOINPUT on all three redirectors
  - `failed_cmsrun`     — cmsRun non-zero exit (needs investigation, not blind resubmit)
  - `failed_xrdcp`      — cmsRun ok, output never landed on ceph
  - `success`           — output on ceph and JSON exit 0
- Writes `resubmit_${pd}_${era}_${TIMESTAMP}.txt` containing only the
  `unrun ∪ failed_input_open ∪ failed_xrdcp` set. `failed_cmsrun` is
  held back for manual review (never auto-resubmitted).
- Resubmit is just `./submit_batch.sh charmonium 2016G
  --list resubmit_charmonium_2016G_....txt`; the resume filter takes
  care of the "don't rerun successful" invariant.

### Site-list refinement (input to submit files)

Validation showed 40% failure at input-open time from a subset of
Tier-2s. `site_deny.txt` starts with the validation offenders (task
0.5 pulls them from the per-cluster condor logs) and is intersected
with the current site list at submit time via a small `sed`.

### Bookkeeping and closeout

- `build_manifest.py` runs after each cluster drains: reads the
  per-file JSONs on ceph + the condor log, emits
  `manifests/${pd}_${era}_cluster_${clusterId}.json`.
- `merge_manifests.py` produces era-level and PD-level manifests
  after all clusters close.
- Success criterion: for each (pd, era), the union of successful
  outputs covers 100% of the populated filelist, no `failed_cmsrun`
  bucket remains unexamined.

### Ceph quota check (must land before first submit)

The full-era output volume estimated from validation:

- Charmonium: 44k files × 17 MB = ~750 GB
- SingleMuon: ~130k files × 5 MB = ~650 GB
- Combined: ~1.4 TB

Add a preflight step (`quota.sh`) that queries current usage and
free space at `/data/user/p/pmlugato/` before any large batch.

## Impact

- Affected specs:
  - `alcareco-production` (new) — launch mechanics: era boundary
    definition, filelist provenance, worker retry contract,
    per-file manifest, resubmit contract, output layout.
- Affected code:
  - `condor/full_alcareco_2016postvfp/` — new directory (see
    layout above); no changes under existing
    `condor/multichannel_alcareco/`.
  - No changes under `CMSSW_10_6_17_patch1/`.
  - No changes under `narf/`, `rabbit/`, `wremnants/`, `scripts/`.
- Reuses without modification:
  - Payload contents (same `Alignment/CommonAlignmentProducer` +
    `Configuration/*` tree that validation shipped).
  - Both recoskim cfgs, unchanged.

## Non-goals

- Physics/producer/cff changes. Locked as of
  `add-preprod-multichannel-validation`.
- Downstream stage-2 CVH consumption of these outputs. Handled in
  the CVH change family.
- Golden JSON / lumi mask application at recoskim time. Deferred to
  downstream stage; recoskim writes all events and downstream
  filters by validated lumi.
- Merging per-file outputs into per-run/per-era ROOT files
  (`hadd` step). Downstream consumers ingest per-file directly.
- Auto-retry inside HTCondor (max_retries) or DAGman orchestration.
  Manual `find_missing.py` chosen — simpler and transparent.
- Rucio tape recall automation. If any 2016F/G RAW files are
  TAPE-only in Rucio, they surface as `failed_input_open` and the
  user issues the recall separately.

## Open questions (to resolve during apply)

1. **PostVFP boundary run**. Task 0.1 will pull from `wremnants/`
   or CMS analysis convention. Confirmed common choice: run ≥
   278820 (HIP-fix boundary). Value gets recorded in
   `lfns/postvfp_run_boundary.txt` for reproducibility.
2. **RequestDisk sizing under stage-in**. RAW files can reach
   ~5 GB; outputs are ~20 MB combined. Bump `RequestDisk` from
   the validation value (8 GiB) to 16 GiB to cover the largest
   RAW + build tree + outputs with headroom. Confirmed as part of
   task 2.5.

Resolved during this iteration:
- Redirector strategy → **stage-in via xrdcp, no streaming** (user);
  AAA fallback lives at the xrdcp step.
- Watchdog → **flat 2-hour cap** regardless of PD (user); no per-PD
  scaling.
- JpsiMuMu redundancy → accepted (user); we re-run to pick up
  persist-alignment additions.
- Rehearsal size → 200 files/PD/era (user).
- Raw TSV storage → kept in-repo under `lfns/raw/` (user).
