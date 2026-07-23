# Design: add-jpsi-x-condor-production

## Goal recap

Run preset B and preset C of the J/ψ+X ALCAReco producer on 2016H Charmonium RAW under HTCondor at the MIT submit cluster, at ~50 000 events per preset, producing the mass-shape plots and accurate per-event runtime/output diagnostics the Phase-2 interactive 5 000-event runs couldn't deliver.

## Core decisions

### 1. RAW input, not MiniAOD

`recoskim_Run2016H_Charmonium_JpsiX_preset{B,C}.py` are cmsDriver outputs whose `process.source.fileNames` points at a `/Charmonium/Run2016H-v1/RAW` file by construction. ALCAReco runs as part of the RAW→RECO+ALCA chain — it needs RECO-level objects from raw digi unpacking. The MiniAOD list in `old_condor_stuff/data_minilist.txt` (which fed the nano production) does **not** apply here. The Condor filelist `raw_2016H_list.txt` is therefore generated against the `/Charmonium/Run2016H-v1/RAW` dataset, sourced from CERN EOS at `/eos/cms/store/data/Run2016H/Charmonium/RAW/v1/000/283/<subdir>/00000/*.root` and accessed remotely via `root://eoscms.cern.ch/`. **Out of the 19 run-number subdirs under `283/`, only 4 actually contain data on EOS** — runs 283270 / 283283 / 283305 / 283682 — so the 5-file selection spans those 4 runs (the other 15 subdirs are empty leftovers from invalidated/aborted runs).

This is the **single most important difference** from the nano-production Condor setup the user referenced. Reuse the submit-file layout and the wrapper skeleton; do not reuse the filelist or assume MiniAOD-style input.

### 2. CMSSW payload — package tarball, not full release

The modified `Alignment/CommonAlignmentProducer/` package is what diverges from upstream `CMSSW_10_6_17_patch1`. Two options:

| Strategy | Tarball size | Build time on worker |
|---|---|---|
| Full release tarball | ~3–5 GB | 0 (no build) |
| Package-only + fresh `scram p` on worker + `scram b -j 8` | ~3–5 MB | ~3–5 min |

The package-only strategy wins on transfer cost by three orders of magnitude. The 3–5 min `scram b` on the worker is dwarfed by the cmsRun runtime (~30 min at `maxEvents=10000`). This is the same pattern the nano production used (`scram p CMSSW_10_6_26` on the worker, then untar `Bmm5.tgz` over the source tree).

### 3. `maxEvents=10000` per job, 5 files × 2 presets = 10 jobs

Two motivations.

First, **predictable accounting**: 5 files × 10 000 events = 50 000 events per preset exactly, the headline figure of this scale-up. With `maxEvents=-1` the per-preset total depends on the file sizes (100–300 k events/file each, so 500 k–1.5 M total per preset), which is 10–30× the stated target.

Second, **predictable wall-clock**: at the Phase-2 measured ~10 ms / event (averaged across the producer chain), 10 000 events ≈ 100 s of pure event-loop time. Adding ~3 min `scram b` and ~2 min CMSSW startup/conditions = ~7 min per job. This sits well below typical T2 walltime budgets and lets the cluster turn the entire 10-job batch around in ~10–15 min wall-clock once jobs start.

If the cluster is fast and the data is well-positioned, scaling to `maxEvents=-1` for the next-scale follow-up is a one-line wrapper change. If we want to go to 500 k events/preset, we keep `maxEvents=10000` and grow the filelist to 50 entries — same job size, more jobs.

### 4. Sequential build steps on worker, not pre-built tarball

The `scram b -j 8` happens on the worker, after untarring the package. This is intentional:

- Avoids ABI mismatches between the submit machine and the worker (different glibc, gcc, etc., possible on T2 sites — though `+SingularityImage = cms:rhel7` mitigates this).
- Makes the payload portable across sites (no compiled `.so` in the tarball — only `.py`, `.cc`, `.h`).
- Keeps the tarball reproducible (no build timestamps in the artefact).

Same pattern as the nano-production wrapper.

### 5. JSON sidecar per job for runtime / output diagnostics

The 5 000-event interactive runs gave one runtime number per preset, total. At Condor scale we need per-job numbers to characterize the variance and to detect outliers (slow workers, slow input pull, slow output xrdcp).

The sidecar is tiny and trivial to aggregate:

```json
{
  "preset": "B",
  "input_lfn": "/store/data/Run2016H/Charmonium/RAW/v1/...",
  "input_size_mb": 1234.5,
  "n_in": 10000,
  "n_out": 8421,
  "n_cands_bplus": 712,
  "n_cands_bzero_kstar": 358,
  "n_cands_bs_phi": 31,
  "n_cands_bc": 42,
  "n_cands_b0_ks": 89,
  "n_cands_lambdab": 12,
  "n_cands_psi2s": 28,
  "runtime_s_cmsrun": 612.4,
  "runtime_s_total_job": 728.1,
  "peak_rss_mb": 2410,
  "exit_code": 0,
  "host": "submit50.t2.something",
  "scram_b_seconds": 184
}
```

(`n_cands_*` are read out of the producer's `EventCounter` plugin or by post-loop `edmFileUtil`/PyROOT parse of the output ROOT — implementation detail in tasks.md item 6.1.)

### 6. Output path mirrors the nano layout

`/data/user/p/pmlugato/mz/alcareco/jpsi-x/preset_<X>/Run2016H/<basename>.{root,json}` — parallel to `/data/user/p/pmlugato/mz/<era>/...` used by the nano production. Same xrdcp redirector (`root://submit50.mit.edu/`), same proxy file, same accounting group, same site list. The user already operates this path layout for nano — minimizes new operational surface.

## What this proposal explicitly does NOT do

- **No cut iteration.** Iter-1 values from `add-jpsi-x-vertex-fit-and-low-pt` are frozen. If the 50 k run shows preset C is still too tight (or too loose), the response is a sibling cut-iteration change, not edits inside this one.
- **No producer code changes.** Bit-identical to what the sibling change ships.
- **No new cmsDriver invocation.** The two existing `recoskim_Run2016H_Charmonium_JpsiX_preset{B,C}.py` configs are the cmsRun inputs; only `process.source.fileNames` and `process.maxEvents.input` are patched per-job by the wrapper.
- **No multi-era expansion.** 2016H only for the first scale-up. 2016B–G and 2017/18 are explicitly out of scope.
- **No CRAB.** Local HTCondor at MIT submit only — the wrapper-script pattern is established and easier to debug at this scale than full CRAB submission.

## Risks and mitigations

| Risk | Mitigation |
|---|---|
| RAW input file not accessible from worker | Primary redirector `root://eoscms.cern.ch/`; on persistent failure fall back to `cmsxrootd.fnal.gov`. xrdcp retry loop wraps the primary; on final failure, JSON sidecar `exit_code != 0` so aggregator sees it. Worker x509 proxy is the same CMS VOMS proxy used by the nano production — already authenticates to CERN EOS as a CMS user. |
| Worker `scram b` failure (mismatched gcc, missing tool) | `+SingularityImage = cms:rhel7` (matches nano-production setup); if a worker is fundamentally incompatible the job exits non-zero and the aggregator flags it. |
| Env-var bug recurrence (`TKALJPSIX_SELECTION_PRESET` not propagating to `cmsRun`) | Wrapper exports it on the same shell line as `cmsRun` invocation; smoke-test job in task 4.1 validates this before the full submission. Also explicitly assertable: the JSON sidecar's first job dump of `process.ALCARECOTkAlJpsiXBPlusCandidates` cut values during a one-off interactive verify run. |
| Output `xrdcp` to `submit50.mit.edu` fails | Retry loop; on persistent failure, write to local Condor scratch and signal in JSON sidecar so we can pick it up manually. |
| Per-job runtime exceeds T2 walltime | Capped by `maxEvents=10000` — at 10 ms/evt that's < 2 min event-loop. If it doesn't fit in the worker's walltime budget the issue is CMSSW startup / conditions, not the producer; falls back to running locally on submit machines. |
| Preset C output size larger than expected at scale | Phase-2 estimate: ~3 MB per 100 events for preset C. 50 k events → ~1.5 GB per preset, well under quota. If actual is 10× the estimate (rare combinatorial blowup), the aggregator catches it and we revisit cuts. |

## Why this is a separate proposal from `add-jpsi-x-vertex-fit-and-low-pt`

The sibling change is about the *cut design* — what makes preset C physically meaningful. This change is about the *production infrastructure* — how to run preset B and C at scale on a cluster and report what they actually do at scale. The two specs have different ownership, different validation paths, and different archive cadences. Co-locating them would make the sibling change's specs validate against operational artefacts that have nothing to do with the cut decisions.

The single requirement added to `alcareco-jpsi-x` is specifically about the production-wrapper contract (env-var propagation, file-name override, output layout, JSON sidecar) — not about cuts or producer behavior. That requirement stands or falls on the wrapper script and the `.sub` files only.
