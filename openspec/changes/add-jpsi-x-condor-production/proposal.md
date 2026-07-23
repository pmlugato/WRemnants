# Change: Condor production for J/ψ+X ALCAReco — presets B and C at ~50k events/preset on 2016H

## Why

The Phase-2 work in `add-jpsi-x-vertex-fit-and-low-pt` (5 000-event interactive runs on 2016H Charmonium RAW) established that the iter-1 cut relaxation makes preset C operationally viable (O(50–150) cands per non-V0 channel, runtime ~50 s / 5 000 evts, output ~3× preset B). It did **not** produce the dimuon-plus-bachelor mass-peak demonstrations the group wants to see — 5 000 events is below the per-channel statistics needed for a clean shape, especially in the rarer modes (Bs→φ, Bc).

The natural next step is to scale up under Condor and run **both preset B and preset C in parallel**, for two complementary reasons:

1. **Mass shapes** — at ~50 k events per preset we expect O(500–1 500) cands in B+ and B0→K*0 across both presets, enough for a believable mass peak on top of the falling combinatorial; Bs→φ and Bc remain stats-limited but visible as a starting point.
2. **Accurate runtime / output diagnostics at scale** — the 5 000-event numbers extrapolated linearly, but interactive jobs are dominated by CMSSW startup, ESSetup, and conditions cache warm-up. A real Condor job profiles the per-event cost properly and exposes any I/O / memory hotspots that don't show up at 5 000 events.

Preset A (no-Kalman baseline) was fully characterized at 5 000 events and is **not** included in the scale-up: its per-event cost and yield are well-known, and burning cluster cycles re-confirming the baseline buys nothing.

## What changes

This change introduces a Condor production wrapper around the existing `recoskim_Run2016H_Charmonium_JpsiX_preset{B,C}.py` cmsRun configs. **No physics-level cuts change** — the producer code, the cff, the presets, and the env-var dispatch are all bit-identical to what `add-jpsi-x-vertex-fit-and-low-pt` ships. The only new artefacts are submission scripts, a per-job wrapper, an input filelist, an aggregator, and the runtime/output reporting at scale.

### 1. New `condor/jpsix_alcareco/` directory at repo root

Five files, modeled on `old_condor_stuff/` (which served the nano production):

- `condor_jpsix_alcareco_presetB.sub` — `.sub` template for preset B, `queue filename from miniaod_2016H_list.txt` (despite the name the list contains **RAW** files; see note below)
- `condor_jpsix_alcareco_presetC.sub` — identical except `arguments` carries `presetC` and the output subpath is `/preset_C/`
- `run_jpsix_alcareco.sh` — worker script: fresh `scram p CMSSW CMSSW_10_6_17_patch1`, untar modified `Alignment/CommonAlignmentProducer` over the package, `scram b -j 8`, export `TKALJPSIX_SELECTION_PRESET`, patch `fileNames` in the cmsRun config to the per-job input file (prefixed `root://eoscms.cern.ch/`), `cmsRun`, `xrdcp` the output and a small JSON sidecar
- `raw_2016H_list.txt` — newline-separated list of 5 (initially) RAW 2016H Charmonium files. Source: CERN EOS at `/eos/cms/store/data/Run2016H/Charmonium/RAW/v1/000/283/<subdir>/00000/*.root` (~2.0 TB total across 19 run-number subdirs). Probing the 19 subdirs via `xrdfs eoscms.cern.ch ls` showed that **only 4 actually contain data** (runs 283270, 283283, 283305, 283682) — the other 15 are empty containers from invalidated / aborted runs. Final 5-file selection spans those 4 runs (1 file each from 283270/283/305 plus 2 from 283682) — still gives 4 different runs of diversity for free. NOT MiniAOD — ALCAReco recoskim consumes RAW. Listing is done remotely from the submit machine via `xrdfs eoscms.cern.ch`, no lxplus dependency.
- `build_tarball.sh` — packages `CMSSW_10_6_17_patch1/src/Alignment/CommonAlignmentProducer/` plus the two cmsRun configs into `jpsix_alcareco_payload.tgz` for `transfer_input_files`

### 2. Per-job output and aggregation

Each job writes `jpsix_alcareco_preset<X>_<inputbasename>.root` (the AlCaReco output) and `jpsix_alcareco_preset<X>_<inputbasename>.json` (counters: n_in, n_out, runtime_s, peak_RSS_MB, exit_code). Both `xrdcp`'d to `root://submit50.mit.edu//data/user/p/pmlugato/mz/alcareco/jpsi-x/preset_<X>/Run2016H/`.

A new `condor/jpsix_alcareco/merge_and_report.py` does, after all jobs complete:

- `hadd` per-preset ROOTs into one preset-merged ROOT
- aggregate the JSON sidecars into a markdown table (mean/total/per-event runtime, output size, n_cands / channel)
- emit mass-overlay plots (extending `CMSSW_10_6_17_patch1/jpsix_preset_compare_3way.py` to many input files; B-only path)

### 3. New OpenSpec capability requirement

A single ADDED requirement under the existing `alcareco-jpsi-x` capability captures the production-wrapper contract: env var consumed on `cmsRun` line of the wrapper script (not just the `.sub` `environment` field), per-job RAW-input override of `process.source.fileNames`, output path layout, runtime/output JSON sidecar schema, and the no-physics-change invariant (Condor wrapper MUST NOT modify any cut value or producer parameter; deviations from the existing cff are out of scope and rejected).

## Scope and explicit non-scope

**In scope.** Condor `.sub` files for presets B and C; worker wrapper script; CMSSW-package tarball strategy; input filelist generation; per-job JSON sidecar; aggregation + reporting; updates to `slides/jpsix-producer-progress.tex` to add a "Phase-3 at scale" frame once results are in.

**Out of scope.**

- Any change to producer code or cff (those live in `add-jpsi-x-vertex-fit-and-low-pt` and earlier).
- Any further cut iteration (`add-jpsi-x-vertex-fit-and-low-pt` iter-1 values are frozen for this scale-up; if the 50 k run reveals a new problem, it goes in a sibling change).
- Preset A on Condor (deliberately excluded).
- 2016 eras other than H, and other years (deliberately excluded for the first scale-up; a follow-up change can add them once 2016H succeeds).
- Bigger scale (500 k, full era, full year) — explicit follow-up work after the 50 k run is validated.

## Impact

- New artefacts: 5 files in `condor/jpsix_alcareco/`, 1 spec delta, 1 new `slides/jpsix-producer-progress.tex` frame after results land.
- No edits to existing producer code, cff, or cmsRun configs.
- Cluster impact: 5 files × 2 presets × ~1 worker/job × ~30–60 min/job ≈ 5–10 worker-hours on T2/T3.
- Output volume: from Phase-2 sizing, preset B at ~1 MB/100-evt and preset C at ~3 MB/100-evt extrapolated → preset B ~500 MB total, preset C ~1.5 GB total at 50 k events/preset.

## Open question deferred to design

Whether `maxEvents` per job should be `-1` (full input file, ~100–300 k events/file → 30–60 min/job) or `10000` (consistent ~50 k events/preset across exactly 5 files). The proposal commits to `maxEvents=10000` per job as the starting point so the resulting per-preset total is exactly the figure quoted in this proposal (5 × 10 k = 50 k); the wrapper script makes this trivially overridable for the next-scale follow-up. See `design.md`.
