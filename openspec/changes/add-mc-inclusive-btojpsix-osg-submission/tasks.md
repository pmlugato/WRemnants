# Tasks

## 0. Prerequisites

- [x] 0.1 Confirm valid CMS x509 proxy on submit50:
      `voms-proxy-info -e -valid 24:00 --file /home/submit/pmlugato/x509up_u239501`.
      If missing/expired: `voms-proxy-init -voms cms -valid 168:00`.
- [x] 0.2 Verify parent directory `condor/mc_inclusive_btojpsix_2016postvfp/`
      exists and contains `run.sh`, `mc_inclusive_payload.tgz`,
      `mc_inclusive_btojpsix_fullchain.py`, `find_missing.py`,
      `build_manifest.py`, `merge_manifests.py`, `quota.sh`, and
      `manifests/base_seed.txt` + `manifests/reco_tags.txt`.
- [x] 0.3 Read-only check on the parent's shipped `.sub`
      (`submits/base_mc_inclusive.sub`) — confirm the delta list
      in design.md D3–D6 still matches what's live.
- [x] 0.4 Ceph preflight for the OSG rehearsal subdir:
      `df -h /data/user/p/pmlugato/`. Confirm ≥ 5 GB free
      (100 jobs × ~50 MB ALCARECO output ≈ 5 GB rehearsal volume).

## 1. Scaffold the OSG launch directory

- [x] 1.1 Create
      `condor/mc_inclusive_btojpsix_2016postvfp_osg/` with
      subdirectories `submits/` and `manifests/`.
- [x] 1.2 Symlink from parent directory (relative paths so a
      moved parent doesn't break the OSG dir):
      `run.sh`, `mc_inclusive_payload.tgz`,
      `mc_inclusive_btojpsix_fullchain.py`, `find_missing.py`,
      `build_manifest.py`, `merge_manifests.py`, `quota.sh`.
- [x] 1.3 Symlink `manifests/base_seed.txt` and
      `manifests/reco_tags.txt` from the parent — these are
      shared invariants across both rehearsals.
- [x] 1.4 Verify all symlinks resolve
      (`ls -la` shows targets exist) and are readable/executable
      as their parent counterparts.

## 2. OSG submit template

- [x] 2.1 Copy `submits/base_mc_inclusive.sub` from the parent to
      `submits/base_mc_inclusive_osg.sub`. Apply deltas per
      design.md D2–D6:
      * ADD `+ProjectName = "MIT_submit"`.
      * ADD `+JobDurationCategory = "Medium"`.
      * ADD `+SingularityBindCVMFS = True`.
      * REPLACE the `BOSCOCluster =!= …` requirements clause with
        `Requirements = ( OSGVO_OS_STRING == "RHEL 7" && HAS_SINGULARITY == TRUE && HAS_CVMFS_singularity_opensciencegrid_org == True && HAS_CVMFS_cms_cern_ch == True )`.
      * REMOVE `+AccountingGroup = "analysis.plugato"`.
      * KEEP: `+SingularityImage`, `use_x509userproxy`,
        `x509userproxy`, `executable`, `should_transfer_files`,
        `when_to_transfer_output`, `transfer_input_files`,
        `transfer_output_files`, `RequestCpus`, `RequestMemory`,
        `RequestDisk`, `Universe`, `arguments`.
- [x] 2.2 Add a header comment documenting the delta from the
      parent `.sub` (design.md links) and the rehearsal-only
      scope.
- [x] 2.3 Static-check the `.sub` for syntax by running
      `condor_submit -dry-run /dev/stdout -queue 1
      submits/base_mc_inclusive_osg.sub`; confirm it prints a
      valid submit description with no error.

## 3. OSG entry-point and reconcile scripts

- [x] 3.1 `submit_batch_osg.sh`: thin wrapper. Reads N_JOBS from
      argv, ceph subdir `osg_rehearsal/`, submit template
      `submits/base_mc_inclusive_osg.sub`. Renders the queue
      block, submits, appends the returned cluster ID to
      `manifests/osg_rehearsal_clusters.txt`. Behavior matches
      the parent's `submit_batch.sh` except for the two paths.
- [x] 3.2 `SUBMIT_OSG.sh`: one-liner entry. Only accepts
      `rehearsal` as argument (100 jobs); any other argument
      exits non-zero with a message that phase-1 is not
      supported from this directory (parent's SUBMIT.sh is
      where phase-1 lives after retune).
- [x] 3.3 `RECONCILE_OSG.sh`: wraps `find_missing.py` with
      `--ceph-dir root://submit50.mit.edu//data/user/p/pmlugato/mz/mc/inclusive_btojpsix_2016postvfp/osg_rehearsal/`
      and `--clusters-file manifests/osg_rehearsal_clusters.txt`.
      Emits `manifests/osg_resubmit_<TS>.txt` scoped to only
      the OSG rehearsal's `(clusterId, procId)` pairs. Does NOT
      touch David's MIT-T3 clusters or the parent's ceph subdir.
- [x] 3.4 Set the executable bit on all three scripts and confirm
      shebangs are `#!/bin/bash`.

## 4. Ceph subdir plumbing

- [x] 4.1 Determine minimal `run.sh` diff needed to write output
      under `osg_rehearsal/` (options: (a) `run.sh` reads
      `OUTPUT_CEPH_SUBDIR` env var, (b) `submit_batch_osg.sh`
      passes the subdir as a 4th positional argument that
      `run.sh` already accepts). Pick whichever changes the
      fewest lines in the shared `run.sh`.
- [x] 4.2 Implement the chosen path. If `run.sh` needs a diff:
      apply it in the parent's `run.sh` (since the OSG dir
      symlinks to it) — the diff must be a no-op for the parent's
      MIT flow (default subdir = "" preserves the current path).
- [x] 4.3 `mkdir -p` the ceph rehearsal subdir via a one-shot
      `xrdfs submit50.mit.edu mkdir` inside `SUBMIT_OSG.sh`
      before the first job is queued.

## 5. Documentation

- [x] 5.1 `README_OSG.md`: brief runbook (us-facing, not
      operator-facing). Sections: purpose (routing rehearsal,
      not phase-1), the one command (`./SUBMIT_OSG.sh rehearsal`),
      the one reconcile command (`./RECONCILE_OSG.sh`), the
      symlink map from the parent, the four-tuple that will be
      retuned after the rehearsal.
- [x] 5.2 Cross-link this OSG directory from the parent's
      `HANDOFF_README.md`? No — the operator handoff is
      independent of what we do here; do NOT modify the parent's
      shipped README.

## 6. Dry-run and submission

- [x] 6.1 `./SUBMIT_OSG.sh rehearsal --dry-run` (add `--dry-run`
      flag to `SUBMIT_OSG.sh` in §3.2 if not already there):
      confirm the rendered queue has exactly 100 entries and the
      cluster placeholder is correct. Do not actually submit.
- [x] 6.2 First live submission: `./SUBMIT_OSG.sh rehearsal`.
      Record the returned cluster ID and expected drain time
      (Medium hard cap 20 h; expected 5–10 h steady state).
- [ ] 6.3 `condor_q -analyze <clusterId>` after 15 min: if idle
      count is > 90 (matching stalled), inspect the analyze
      output for narrow-slot warnings; if `RequestMemory` is
      the flagged blocker, hold the cluster, edit the `.sub` to
      `RequestMemory=4000`, resubmit, log the change in the
      reconcile write-up (Open Question §1 in the proposal).

## 7. Reconcile and write-up

- [ ] 7.1 Watch the cluster to completion with periodic
      `condor_q -all -run` snapshots into
      `manifests/osg_progress.log`.
- [ ] 7.2 After drain: `./RECONCILE_OSG.sh` and inspect the
      output. Auto-safe resubmit list must be empty or small
      (<10% of 100).
- [ ] 7.3 If auto-safe resubmits are present, resubmit once via
      `submit_batch_osg.sh --resubmit manifests/osg_resubmit_<TS>.txt`.
      Drain, reconcile again. Report the two-cycle success rate
      (target ≥ 95%).
- [ ] 7.4 Aggregate per-job JSONs into a single manifest via
      `python build_manifest.py --clusters manifests/osg_rehearsal_clusters.txt --out manifests/osg_manifest.json`.
- [ ] 7.5 Extract p99 wall, p99 RSS, p99 scratch, and mean
      s/event across the 100 jobs. Compute what N_GEN would
      look like at target_wall=8 h given these numbers.
- [ ] 7.6 `osg_rehearsal_report.md` under
      `condor/mc_inclusive_btojpsix_2016postvfp_osg/`:
      documents (a) submission cluster IDs, (b) matching latency
      (idle-to-run p50/p99), (c) two-cycle success rate, (d) p99
      wall/RSS/scratch, (e) failure taxonomy from the reconcile,
      (f) recommended `.sub` deltas for phase-1 based on the
      measurements, (g) any deviation from D6 placeholder values
      (e.g. if `RequestMemory` had to drop).
- [ ] 7.7 If phase-1 numbers from David's MIT-T3 rehearsal have
      already landed, add a comparison section: MIT vs OSG on
      the four-tuple. Otherwise defer that section until they do.

## 8. Approval gates

- [ ] 8.1 User approves this proposal before §1 begins.
- [ ] 8.2 User approves the rendered `.sub` (via §2.3 dry-run)
      before §6.2 live submission.
- [ ] 8.3 User approves the `osg_rehearsal_report.md` before
      any sizing retune is applied to the parent's shipped `.sub`.

## 9. Close-out

- [ ] 9.1 Consolidate `osg_rehearsal_report.md` findings into the
      project memory ([[project-mc-inclusive-btojpsix]]) so the
      phase-1 decision has both rehearsals' data in one place.
- [ ] 9.2 `openspec archive add-mc-inclusive-btojpsix-osg-submission`
      once the report is approved and (per §5.2) no parent
      artifacts have been mutated.
