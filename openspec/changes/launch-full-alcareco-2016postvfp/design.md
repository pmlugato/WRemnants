# Design: Launch full 2016 postVFP AlCaReco production

## Context

Validation (`add-preprod-multichannel-validation`) established that
the four AlCaReco streams run cleanly on 100 files per PD for
2016H. This change lifts that from a one-shot rehearsal to a
recurring launch pattern that spans three eras, handles ~180k
input files, and survives typical Grid failure modes without
duplicated compute.

## Decisions

### D1. Manual `find_missing.py` reconcile, not DAGman

DAGman gives per-node retry with different site preferences and
post-scripts that can check output presence. `find_missing.py` gives
none of that automatically but is transparent — the user reads a
plain resubmit list and re-issues a submit with the same driver.

Rationale for the manual path:

- Failure surface is bimodal. `EX_NOINPUT` is the dominant failure
  and is site-local, not input-local; a blind retry from HTCondor
  hits the same site pool. Our fix must reshape `DESIRED_Sites`,
  which DAG per-node retry cannot express.
- `cmsRun` non-zero exits demand human inspection — an auto-retry
  would either mask a real bug or produce indistinguishable output
  from a benign transient. Manual triage matches the physics-review
  expectation.
- Adds no new dependency (DAGman needs a scheduler-side setup at
  submit06 we haven't audited).

Cost: one extra command per resubmit cycle. Acceptable at the
expected cycle count (2–3 passes).

### D2. Stage-in via xrdcp, no streaming reads inside cmsRun

The wrapper xrdcp's the input LFN to `$_CONDOR_SCRATCH_DIR` before
cmsRun starts; cmsRun reads a local `file:./<basename>.root`. An
end-of-job trap `rm -f`s the local copy on success or failure.

Rationale for the switch away from cmsRun-side streaming:

- **Streaming failure modes are messy**. Xrootd stream stalls
  mid-file surface as `EX_NOINPUT`, `EX_IOERR`, or occasional
  cmsRun spin — no single error class covers them, and cmsRun
  logs make the exact cause hard to bucket.
- **Retry is atomic under stage-in**. `xrdcp` is one call; success
  or failure is unambiguous. AAA fallback (`cms-xrd-global.cern.ch`
  → `cmsxrootd.fnal.gov` → `xrootd-cms.infn.it`) is a simple
  `for` loop around `xrdcp` — no cfg patching, no cmsRun restart.
- **Local read is faster** for the ~2000-event RAW files we care
  about; cmsRun's TTreeCache is more effective on local storage.
- **Cleanup is explicit**. Scratch fills would be catastrophic on
  a shared worker; the trap removes the local copy even on wrapper
  exit paths.

Cost: `RequestDisk` grows (RAW ≤ ~5 GB + build tree + outputs).
Priced into the sub files.

Rejected alternative: pin to one regional redirector (loses global
availability if that redirector is degraded); rely on cmsRun's own
xrootd retry logic (opaque, hard to reason about).

### D3. Per-cluster manifest built at close-out, not per-job

Two alternatives:

- **Per-job**: worker writes manifest entry after xrdcp. Simple, but
  produces 100k+ tiny JSONs on ceph that a single reconcile has to
  fetch (slow xrdfs ls).
- **Per-cluster** (chosen): after `condor_wait`, `build_manifest.py`
  reads the local condor logs and the ceph directory in one pass,
  producing one JSON per cluster. Fast reconcile, single source of
  truth per submission.

The per-file JSONs from `run.sh` are still xrdcp'd (unchanged from
validation) — they carry the runtime and RSS numbers we need for
resource plots. But the manifest is a separate, cluster-level roll-up.

### D4. `site_deny.txt` seeded from validation, updated on drain

The failing-site set from cluster 3227852/3227853 is captured once
into `site_deny.txt`. Each subsequent cluster: after drain, we grep
the condor logs for `EX_NOINPUT` and add new offenders. Submit
files interpolate the file into `DESIRED_Sites` via `sed`. This is
the simplest closed-loop way to shrink the failure rate without
hand-editing sub files across submissions.

### D5. Filelist provenance from DAS, not EOS

The user asked whether `ls /eos/cms/store/data/Run2016H/{PD}/RAW/v1/`
would suffice. Answer: no.

- 2016 RAW is Rucio-managed. The CERN EOS namespace may not hold a
  replica of every file — replicas are distributed across US and EU
  Tier-2s. `ls` at CERN EOS undercounts.
- Zero-event RAW placeholder files exist in the namespace and are
  invisible to `ls` (indistinguishable from populated files by
  filename).
- DAS carries the authoritative dataset membership plus per-file
  event counts in a single query.

Query pattern locked to `dasgoclient --query "file dataset=... |
grep file.name, file.nevents"`, with `awk '$2 >= 1000'` filter for
populated files (same threshold as validation).

### D6. Era boundary in a data file, not hard-coded

The postVFP run boundary is a physics decision, not a launch
decision. Committing it as `lfns/postvfp_run_boundary.txt` (a single
integer + a comment describing where it came from) makes the
convention explicit and diffable; hard-coding into a shell script
buries it.

### D7. Per-era submit files, one payload tarball

Six `.sub` files (2 PDs × 3 eras) sharing one payload tarball.
Alternatives considered:

- **One .sub per PD, era as filter in run.sh** — hurts cluster-ID
  scoping (a resubmit of just 2016G becomes ambiguous).
- **One giant queue over all eras** — same problem plus 45k-line
  queue files are cumbersome to reason about.

Six .sub files, each with a scoped `queue filename from
${pd}_${era}.txt`, keeps clusters clean and lets us resubmit at era
granularity without disturbing the others.

### D8. Wall-time watchdog set at a flat 2-hour cap

Charmonium median was 491 s/file (p99 ~1200 s), SingleMuon 124 s
(p99 ~300 s), on the 100-file validation. RAW file sizes vary
non-trivially — using 3× the observed mean would leave too little
headroom for a legitimately larger file. A flat 2-hour cap sits
comfortably above the observed p99 for both PDs with enough slack
that a normal long-tail file does not trip the watchdog; the
watchdog's only role is to catch a truly hung process (xrdcp
network stall on output, dead cmsRun spin, misbehaving worker
node) rather than to police per-file runtime variance. Applied
uniformly across both PDs — no per-PD tuning. Writes
`exit_reason: wrapper_timeout` into the JSON so `find_missing.py`
can bucket it.

## Rejected alternatives

- **Copy RAW to ceph first, then read locally.** ~40 TB of RAW,
  hostile to storage quota. Xrootd read from Tier-2s is fine at
  our IO rates.
- **Rewrite as a Snakemake DAG.** Overkill for a two-stage
  submit-and-reconcile flow.
- **Move outputs to Tier-2 SE instead of MIT `/data/user/`.** The
  downstream stage-2 consumer runs at MIT; local storage is the
  right target. Ceph quota headroom is checked before launch (D1
  in tasks.md §0).

## Risks

- **Ceph fills mid-launch.** Mitigation: preflight quota check,
  emit a warning at 80%. Manual pause + free space if triggered.
- **Payload tarball drift between validation and launch.** Rebuild
  from scratch at launch time, not reused from
  `condor/multichannel_alcareco/`. Content is identical (both
  bundle the same CMSSW area) but the launch tarball lives in its
  own directory for provenance.
- **AAA fallback masks a real config bug.** If all three
  redirectors fail identically, `find_missing.py` reports it as
  `failed_input_open` and the human decides — not silently.
- **Rucio TAPE-only 2016F/G files.** Surfacing as
  `failed_input_open`. If the pattern shows up cluster-wide on a
  specific era, we pause and issue a `rucio` recall before
  resubmitting.

## Cross-references

- Physics/producer contract: locked upstream in
  `finalize-jpsi-x-preset-b-production`,
  `add-jpsi-x-muons-and-preprod-refinements`, and
  `add-preprod-multichannel-validation`. This change does not touch
  any of it.
- Validation directory: `condor/multichannel_alcareco/`. Retired
  after this change ships; kept in-repo for provenance.
