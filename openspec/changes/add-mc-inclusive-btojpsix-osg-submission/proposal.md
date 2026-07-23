# Change: Add OSG-Pool submission path for the MC inclusive B → J/ψ + X production

## Why

The parent change [`add-mc-inclusive-btojpsix-2016postvfp`](../add-mc-inclusive-btojpsix-2016postvfp/proposal.md)
ships a handoff tarball that a MIT-side operator (David) runs on the
CMS global pool via MIT T3. That rehearsal has been shipped; we are
now waiting for the operator's reconcile output to firm up the
per-job sizing (`RequestMemory`, `RequestDisk`, `N_GEN`, wall).

While waiting, we can independently exercise the exact same
`run.sh` + payload + physics config against the **OSG Open Science
Pool** (OSPool) from `submit50.mit.edu` using the MIT-provided
`+ProjectName = "MIT_submit"` allocation. Value added:

1. **Independent routing validation.** OSG has a materially
   different slot ecology (different sites, container-enforced
   worker environment, preemption possible). If David's rehearsal
   ships back tight numbers, we still don't know whether OSG can
   match them; running our own 100-job rehearsal is the cheapest
   way to find out **before** we depend on OSG for the phase-1
   4×10⁸-event scale.
2. **A second data point for the sizing retune.** OSG worker
   nodes span a wider hardware range than MIT T3. Wall-clock,
   RSS, and scratch p99 measured on 100 OSPool jobs are a stronger
   input to the phase-1 `.sub` than 100 measurements from one MIT
   T3 CE.
3. **Zero blocking on David.** OSG rehearsal proceeds in parallel;
   we compare notes when both reconciles are in hand.

The scope is strict: **routing/rehearsal only**, no phase-1 launch
from this directory. Phase-1 (which pool, whether OSG or MIT T3
or split) is a follow-up decision after both rehearsals report.

## What Changes

### Scope (locked in with the user)

- **We run this.** Unlike the parent change, no external operator
  is in the loop. Submission is from `submit50.mit.edu`; reconcile
  is on our side.
- **Rehearsal only.** 100 jobs, N_GEN placeholder per job matches
  the shipped MIT `.sub` (currently 800). No `phase1` target in
  this directory.
- **Physics identical to parent.** Same patched fragment, same
  `mc_inclusive_payload.tgz`, same `mc_inclusive_btojpsix_fullchain.py`,
  same base seed, same seed injection math, same `run.sh`. Nothing
  physics-related may diverge — the point is to exercise the
  routing under the exact same job that David's operator is running.
- **Output isolated.** Rehearsal outputs land under
  `/data/user/p/pmlugato/mz/mc/inclusive_btojpsix_2016postvfp/osg_rehearsal/`
  (subdir under the parent's ceph root). Keeps David's and our
  outputs cleanly separable without a second ceph root.
- **Resource-request placeholder.** The OSG `.sub` ships with the
  same numbers as the MIT `.sub` for the rehearsal
  (`RequestMemory=8000`, `RequestDisk=16 GiB`, `RequestCpus=4`,
  N_GEN=800). We deliberately do NOT re-tune here; retune happens
  once (both for MIT and OSG) after David's reconcile lands. The
  100-job OSG rehearsal is a routing test, not a sizing test.

### New launch directory `condor/mc_inclusive_btojpsix_2016postvfp_osg/`

Mirrors the shipped MIT directory but with an OSG-specific submit
template and a `SUBMIT_OSG.sh` / `RECONCILE_OSG.sh` pair scoped
to `rehearsal`. Reused verbatim from the parent directory:

- `run.sh` — worker (unchanged)
- `mc_inclusive_payload.tgz` — CMSSW payload (unchanged)
- `mc_inclusive_btojpsix_fullchain.py` — full-chain cfg (unchanged)
- `find_missing.py`, `build_manifest.py`, `merge_manifests.py` — reconcile
- `quota.sh` — preflight
- `manifests/base_seed.txt` and `manifests/reco_tags.txt` — seed + reco tags

New/OSG-specific:

- `submits/base_mc_inclusive_osg.sub` — OSG submit template (deltas
  from `base_mc_inclusive.sub` listed in design.md).
- `submit_batch_osg.sh` — thin wrapper around `submit_batch.sh`
  that points at the OSG `.sub` and the OSG output subdir.
- `SUBMIT_OSG.sh` — one-liner: `./SUBMIT_OSG.sh rehearsal` (100 jobs).
  No `phase1` target.
- `RECONCILE_OSG.sh` — thin wrapper over `find_missing.py`
  scoped to the OSG output subdir.
- `README_OSG.md` — brief runbook for our own reference.

Reuse strategy is by symlink-or-copy from the parent directory,
so any bugfix in `run.sh` etc. propagates via a single edit
followed by re-linking (documented in tasks §2).

### Submit-file deltas (details in design.md)

Delta from parent `base_mc_inclusive.sub`:

- **ADD**: `+ProjectName = "MIT_submit"` (MIT's OSG allocation).
- **ADD**: `+JobDurationCategory = "Medium"` (mandatory on OSG;
  Medium = target < 10 h, hard cap 20 h — fits N_GEN=800 at
  25–40 s/evt).
- **ADD**: `+SingularityBindCVMFS = True` (needed for
  `/cvmfs/cms.cern.ch` inside the container).
- **REPLACE**: the CMS-pool `Requirements` clause
  (`BOSCOCluster =!= …`) with an OSG clause:
  `Requirements = ( OSGVO_OS_STRING == "RHEL 7" && HAS_SINGULARITY == TRUE && HAS_CVMFS_singularity_opensciencegrid_org == True && HAS_CVMFS_cms_cern_ch == True )`
- **REMOVE**: `+AccountingGroup = "analysis.plugato"` (this is
  the CMS global pool route; `+ProjectName` supersedes it for
  OSPool).
- **KEEP**: `+SingularityImage = "/cvmfs/singularity.opensciencegrid.org/cmssw/cms:rhel7"`
  (validated locally, resolves on OSG).
- **KEEP**: `use_x509userproxy`, `x509userproxy`, executable,
  transfer_input_files, `RequestCpus`, `RequestMemory`,
  `RequestDisk`, universe, arguments — all unchanged.
- The output ceph subdir is set in `submit_batch_osg.sh`
  (`osg_rehearsal/`), not in the `.sub`; `run.sh` reads it
  from an env var already.

### Sizing retune (post-rehearsal, out of this change)

Both rehearsals (David's MIT-T3 and our OSG) will produce per-job
wall/RSS/scratch distributions. The retune of `RequestMemory`,
`RequestDisk`, `N_GEN`, and (potentially) `+JobDurationCategory`
for phase-1 happens **after both reconciles land** and is not
executed inside this change. It will be documented and applied in
a follow-up change or as an addendum to the parent change's tasks.

## Impact

- Affected specs:
  - `mc-production` (extends parent) — adds OSG-routing
    requirements (submit-file attributes, requirements-clause
    contents, output-subdir isolation, we-run-not-operator
    contract, rehearsal-only scope).
- Affected code:
  - `condor/mc_inclusive_btojpsix_2016postvfp_osg/` — new
    directory, mostly symlinks/reuses parent artifacts.
  - `condor/mc_inclusive_btojpsix_2016postvfp/` — unchanged. The
    parent handoff tarball already shipped to David is unaffected.
- No changes to `wremnants/`, `narf/`, `rabbit/`, or the analysis
  pipeline.
- No changes to physics configuration, fragment, seed strategy,
  or ALCARECO content.

## Non-goals

- **No phase-1 launch from this directory.** Rehearsal only.
- **No physics-config divergence from the parent.** Any
  divergence would defeat the "second data point for the same
  job" purpose.
- **No re-tune of resource requests inside this change.**
  Placeholder values match the parent's shipped `.sub`; retune
  happens after both rehearsals report.
- **No OSDF / stashcp / pelican integration.** Payload is 3.7 MB
  (well under OSG's 250 MB `transfer_input_files` cap), so plain
  input transfer works. Output continues to xrdcp back to
  `submit50.mit.edu`, same path shape as the parent, just under
  a `osg_rehearsal/` subdir.
- **No production-grade OSG project.** We use `MIT_submit`; a
  CMS-specific OSG project would be a follow-up if phase-1 is
  decided to route via OSG.
- **No changes to the MIT-side reconcile tooling.** OSG reconcile
  reuses the parent's `find_missing.py` verbatim.

## Open questions (to resolve during apply)

1. **Match rate on OSPool with our requirements.** `OSGVO_OS_STRING == "RHEL 7"`
   + `HAS_CVMFS_cms_cern_ch == True` + `RequestMemory=8000`
   + `RequestCpus=4` may match a narrow slot pool and produce
   large idle times. If steady-state matching stalls, temporarily
   drop `RequestMemory` to 4000 or `RequestCpus` to 1 for the
   rehearsal only, and note the change in the reconcile write-up.
   No physics impact — same code, same seed math, same output.
2. **Preemption behavior.** OSPool jobs can be preempted mid-run.
   Our deterministic seed already makes resubmit safe (the parent
   spec's `Deterministic Per-Job Seed` requirement covers this).
   The rehearsal will tell us empirically what fraction gets
   preempted.
