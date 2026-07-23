# Design: OSG-Pool submission path for the MC inclusive B → J/ψ + X production

## Context

The parent change `add-mc-inclusive-btojpsix-2016postvfp` ships a
handoff to David (MIT-side operator) submitting to the CMS global
pool via MIT T3. This design adds a parallel routing option to the
OSG Open Science Pool (OSPool), which we run ourselves from
`submit50.mit.edu` under the MIT-provided `MIT_submit` project.

The OSG rehearsal is scoped narrowly: **100 jobs, no phase-1**.
Its purpose is (a) prove the routing works, (b) collect a second
independent per-job wall/RSS/scratch distribution to feed the
eventual phase-1 sizing retune. Everything downstream of routing
(physics config, seeds, reconcile, ceph output) is reused verbatim
from the parent.

References:
- MIT submit-users-guide "Batch computing" page (URL:
  `submit.mit.edu/submit-users-guide/running.html`),
  §"Jobs submission to OSG pool".
- OSG documentation on HTCondor submission
  (URL: `portal.osg-htc.org/documentation/htc_workloads/workload_planning/htcondor_job_submission/`).
- Parent proposal and its `mc-production` spec deltas.

## Decisions

### D1: Separate launch directory, aggressive reuse via symlink

`condor/mc_inclusive_btojpsix_2016postvfp_osg/` is a new
sibling directory to the parent's
`condor/mc_inclusive_btojpsix_2016postvfp/`.

**Rationale for the new directory** (vs a second `.sub` inside
the parent):

- The parent's handoff tarball has already shipped to David. Its
  layout is now a public artifact — we must not modify it in
  ways that would invalidate the tarball's referential integrity
  (base_seed, reco_tags, README).
- OSG-specific entry points (`SUBMIT_OSG.sh`, `RECONCILE_OSG.sh`,
  `submits/base_mc_inclusive_osg.sub`) are cleaner in their own
  root than smashed into a directory the operator is actively
  running from.

**Rationale for symlink-heavy reuse** (vs copy-and-drift):

- Any bugfix in the shared `run.sh` should propagate to both
  routing paths in one edit. Copies drift.
- CMSSW payload tarball (3.7 MB) does not benefit from copying;
  a symlink is functionally identical to submit hosts.
- `find_missing.py`, `build_manifest.py`, `merge_manifests.py`
  are stable and shared.

**Files symlinked from the parent**:
`run.sh`, `mc_inclusive_payload.tgz`,
`mc_inclusive_btojpsix_fullchain.py`, `find_missing.py`,
`build_manifest.py`, `merge_manifests.py`, `quota.sh`,
`manifests/base_seed.txt`, `manifests/reco_tags.txt`.

**Files unique to the OSG directory**:
`submits/base_mc_inclusive_osg.sub`, `submit_batch_osg.sh`,
`SUBMIT_OSG.sh`, `RECONCILE_OSG.sh`, `README_OSG.md`, and its
own `manifests/` output store (for JSON manifests of OSG-run
clusters).

### D2: Container = `cmssw/cms:rhel7` (validated locally)

Two candidates from the MIT submit guide:

- `/cvmfs/singularity.opensciencegrid.org/cmssw/cms:rhel7` —
  the CMS-blessed RHEL7 container. What our local validation
  used and what the parent `.sub` already sets. Ships CMSSW
  runtime deps; `cmsset_default.sh` is on the standard path.
- `/cvmfs/singularity.opensciencegrid.org/opensciencegrid/osgvo-el7:latest` —
  MIT-guide default for OSG. Lighter; we would need to source
  `cmsset_default.sh` explicitly and confirm the OSG image
  ships glibc/libc versions CMSSW_10_6_20 wants.

We pick `cmssw/cms:rhel7`: guaranteed to work with the same code
that ran locally, no divergence in the runtime environment
compared to the David rehearsal. `+SingularityBindCVMFS = True`
still needed for `/cvmfs/cms.cern.ch` access.

### D3: `+JobDurationCategory = "Medium"`

OSG requires this field. Choices:

- `Medium` — target < 10 h, hard cap 20 h.
- `Long` — target < 20 h, hard cap 40 h.

At N_GEN=800 and the local NoPU per-event time of 8.5 s (premix
expected 25–40 s), estimated wall per job is 5.5–9 h. `Medium`
fits with margin against the 20 h hard cap. `Long` would broaden
the slot pool and give more room for slow OSG sites, but at the
cost of tighter matching and longer time-in-queue.

We use `Medium` for the rehearsal and revisit for phase-1 after
the rehearsal reports actual p99 wall.

### D4: OSG Requirements clause

```
Requirements = ( OSGVO_OS_STRING == "RHEL 7"
                 && HAS_SINGULARITY == TRUE
                 && HAS_CVMFS_singularity_opensciencegrid_org == True
                 && HAS_CVMFS_cms_cern_ch == True )
```

Why each clause:

- `OSGVO_OS_STRING == "RHEL 7"` — CMSSW_10_6_20_patch1 is SLC7.
  Landing on RHEL 8/9 without CMSSW's SLC7-compat singularity
  would fail at cmsRun.
- `HAS_SINGULARITY == TRUE` — needed to run the container;
  guarantees the slot supports singularity/apptainer.
- `HAS_CVMFS_singularity_opensciencegrid_org == True` —
  the CVMFS repo that hosts our `+SingularityImage`.
- `HAS_CVMFS_cms_cern_ch == True` — required to run CMSSW
  (release area, conditions, generator data). Without this,
  cmsRun cannot resolve `/cvmfs/cms.cern.ch/…`.

Removed from the parent's Requirements:

- `BOSCOCluster =!= …` exclusions — those exclude MIT T3 CEs
  from the CMS global pool; irrelevant for OSPool.

### D5: Output isolation — `osg_rehearsal/` subdir

Both David's rehearsal and our OSG rehearsal write to
`/data/user/p/pmlugato/mz/mc/inclusive_btojpsix_2016postvfp/`.
We isolate OSG outputs to
`.../inclusive_btojpsix_2016postvfp/osg_rehearsal/` so:

- The two rehearsals can drain in parallel without file-name
  collision (they can't collide anyway — deterministic seeds
  encode `clusterId + procId` and clusterIds differ — but the
  subdir makes triage trivial).
- The reconcile scripts on each side see a clean, single-source
  directory.
- Deleting one rehearsal's outputs (e.g. cleanup after retune)
  is one `rm -rf` on the subdir.

The subdir is set at `run.sh` xrdcp time via the env var
`OUTPUT_CEPH_SUBDIR` that `submit_batch_osg.sh` exports (or by
appending it to the `xrdcp` destination in `run.sh` — whichever
requires the smaller `run.sh` diff; determined at apply time).

### D6: Resource requests — placeholder-match to David's rehearsal

The shipped MIT `.sub` (from parent) has:
`RequestMemory=8000`, `RequestDisk=16 GiB`, `RequestCpus=4`,
`N_GEN=800`.

The OSG rehearsal `.sub` ships with the **same** values. Rationale:

- If we tune them independently, we can't compare wall/RSS/scratch
  distributions between the two rehearsals meaningfully.
- Two variables changed (routing AND sizing) makes it harder to
  attribute any difference.
- MIT guide warns that high-memory requests "make it harder to
  find the slot" on OSG; if OSG matching stalls, we drop
  `RequestMemory` to 4000 for the rehearsal only and record the
  change in the reconcile write-up. Physics is unaffected.

The retune of these four numbers for phase-1 is deferred to a
follow-up (see §D7).

### D7: Post-rehearsal retune is out of scope for this change

Both rehearsals will each produce ~100 JSON manifests
(wall-clock, RSS, scratch, host, exit). The phase-1 retune
combines both distributions and re-issues:

- `RequestMemory` = p99(RSS) + 1 GB
- `RequestDisk` = p99(scratch) + 2 GB
- `N_GEN` = floor((target_wall × 3600) / mean(s_per_evt))
- `+JobDurationCategory` = Medium if p99(wall) < 10 h else Long

That retune is a downstream decision. It happens in either the
parent change's tasks (as a §8 addendum) or a fresh change,
depending on how big the delta is. Not modeled here.

### D8: Proxy — reuse existing `x509up_u239501`

`submit50.mit.edu` already has `/home/submit/pmlugato/x509up_u239501`
in use by the parent `.sub`. Same proxy path in the OSG `.sub`.
User is responsible for keeping the proxy valid
(`voms-proxy-init -voms cms -valid 168:00`) before submission.

### D9: Reconcile — reuse `find_missing.py` verbatim

The parent's `find_missing.py` operates on `(clusterId, procId)`
pairs and a ceph destination directory. Both are already
parameterizable. `RECONCILE_OSG.sh` is a thin wrapper that:

- points `--ceph-dir` at `.../osg_rehearsal/`,
- reads the OSG cluster IDs from `manifests/osg_rehearsal_clusters.txt`
  (written by `SUBMIT_OSG.sh` at submission time),
- emits the auto-safe resubmit list scoped to the OSG rehearsal
  (never touches David's clusters or the parent's ceph subdir).

No changes to `find_missing.py` itself.

## Alternatives considered

### A1: A second `.sub` inside the parent's directory (rejected)

Simpler on paper, but:

- Parent tarball has shipped; modifying its directory now means
  a new tarball, means David has an outdated copy.
- Blurs "operator-facing" and "we-run-it" contracts.

### A2: Fully independent OSG directory with copies (rejected)

Rejected because bugfixes in `run.sh` would need to be applied
twice. Symlinks give the same isolation-of-entry-points without
the drift.

### A3: Use `osgvo-el7` as container (rejected)

MIT-guide default and lighter. Rejected because we already have
`cmssw/cms:rhel7` validated locally; a container swap adds an
uncontrolled variable to the rehearsal.

### A4: Route rehearsal via CMS global pool (`+DESIRED_Sites = "T3_US_OSG"`) rather than the OSPool (rejected)

The MIT guide notes that the CMS global pool can include OSG-flavored
sites via `T3_US_OSG` etc. Rejected because:

- Different accounting (CMS global pool vs OSPool proper).
- Different matching semantics (`+AccountingGroup` vs
  `+ProjectName`).
- The user's stated goal is OSG *itself*, not "OSG-like sites in
  the CMS pool".

### A5: OSDF/stashcp for input (rejected)

Payload is 3.7 MB (well under OSG's 250 MB `transfer_input_files`
cap). No benefit from OSDF at this size; adds complexity.

## Risks and mitigations

| Risk | Mitigation |
|---|---|
| Matching stalls due to narrow-slot Requirements+Memory combo | Drop `RequestMemory` to 4000 for rehearsal (physics unchanged); noted in §Open Questions of proposal |
| Preemption on OSG kills a chunk of the rehearsal | Deterministic seed already makes resubmit safe; reconcile buckets `preempted → unrun` for auto-safe resubmit |
| Wall exceeds 10 h `Medium` target on slow sites | Category=Medium hard cap is 20 h; N_GEN=800 leaves margin; if p99 wall > 10 h in the reconcile, switch to `Long` for phase-1 |
| ceph auth from OSG worker to submit50 fails | xrootd path is identical to the MIT rehearsal (same `root://submit50.mit.edu/`); if this fails on OSG but works on MIT it means the OSG site's outbound network is restricted — noted, escalate to submit50 admins |
