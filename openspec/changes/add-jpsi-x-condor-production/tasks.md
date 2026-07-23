# Tasks: add-jpsi-x-condor-production

## 1. Inputs and tarball

- [x] 1.1 Generated the RAW input list from CERN EOS via `xrdfs eoscms.cern.ch ls`. **Surprise finding when probing the 19 run-number subdirs under `/eos/cms/store/data/Run2016H/Charmonium/RAW/v1/000/283/`:** only 4 actually contain data — runs 283270, 283283, 283305, 283682 (the rest are empty leftovers from invalidated / aborted runs). The 2.0 TB the user mentioned is concentrated in those 4. Final 5-file selection: 1 file each from 283270 / 283283 / 283305 plus 2 from 283682 — still gives 4 different runs of diversity at no cost. Filelist at `condor/jpsix_alcareco/raw_2016H_list.txt` (bare LFNs, no `root://...` prefix — wrapper prepends `root://eoscms.cern.ch/` at cmsRun-config-patch time).
- [x] 1.2 Wrote `condor/jpsix_alcareco/build_tarball.sh` and ran it. Tarball includes `BuildFile.xml` + `python/` + `plugins/` + `interface/` + `src/` under the package, plus both cmsRun configs. Output: `condor/jpsix_alcareco/jpsix_alcareco_payload.tgz` — **137 KB, 239 entries** (well under the 5 MB target). Re-runnable.
- [x] 1.3 x509 proxy at `/home/submit/pmlugato/x509up_u239501` verified valid (~118 hours remaining at time of writing); matches what `condor_data_nano.sub` already uses.

## 2. Worker wrapper

- [x] 2.1 Wrote `condor/jpsix_alcareco/run_jpsix_alcareco.sh`. Argument: `$1` = preset letter (`B` or `C`), `$2` = RAW input LFN. Exec bit set.
- [x] 2.2 Worker steps (all implemented in `run_jpsix_alcareco.sh`):
  1. Source CMSSW, `scram p CMSSW CMSSW_10_6_17_patch1`, `cd CMSSW_10_6_17_patch1/src`, `cmsenv`.
  2. Untar `../../jpsix_alcareco_payload.tgz` (creates `Alignment/CommonAlignmentProducer/*` and `recoskim_*.py` two dirs up).
  3. `scram b -j 8`.
  4. Patch the cmsRun config: replace the `fileNames = cms.untracked.vstring('...')` line with the per-job input (`root://eoscms.cern.ch/$2`), and set `process.maxEvents.input = cms.untracked.int32(10000)`. Use `sed -i` for the file-name swap and append `process.maxEvents.input = cms.untracked.int32(10000)` at end of config (idempotent override).
  5. `export TKALJPSIX_SELECTION_PRESET=$1` — on the same line as `cmsRun`, per the env-var-set-on-cmsRun-only bug documented in `add-jpsi-x-vertex-fit-and-low-pt`.
  6. `time cmsRun recoskim_Run2016H_Charmonium_JpsiX_preset$1.py` → capture stderr and exit code.
  7. Emit JSON sidecar: parse n_in (`maxEvents` or stop time), n_out (`edmFileUtil -e` on the output), runtime (the `time` output `real`), peak RSS (`grep VmPeak /proc/<pid>/status` if available, else `ru_maxrss` from `/usr/bin/time -v`).
  8. `xrdcp output.root root://submit50.mit.edu//data/user/p/pmlugato/mz/alcareco/jpsi-x/preset_$1/Run2016H/jpsix_alcareco_preset$1_<basename>.root`. Same for the JSON sidecar.
- [x] 2.3 Worker robustness (implemented):
  - Input pull goes through CERN EOS (`root://eoscms.cern.ch/`) — wrap the cmsRun input access in a retry loop; on persistent failure, fall back to `cmsxrootd.fnal.gov` as a secondary redirector for cases where the file has been replicated.
  - On non-zero `cmsRun` exit, still emit the JSON sidecar with `exit_code != 0` so the aggregator sees the failure.
  - Clean up the scratch CMSSW work area at end (Condor scratch is auto-wiped but be explicit).

## 3. Submit files

- [x] 3.1 Wrote `condor/jpsix_alcareco/condor_jpsix_alcareco_presetB.sub`:
  - Modeled on `old_condor_stuff/condor_data_nano.sub` (sites list, BOSCO exclusions, accounting group, x509, log paths).
  - `executable = run_jpsix_alcareco.sh`
  - `transfer_input_files = jpsix_alcareco_payload.tgz`
  - `arguments = "B $(filename)"`
  - `RequestMemory = 4000`, `RequestDisk = 8 GB` (CMSSW scratch is heavy)
  - `output / error / log = /work/submit/pmlugato/mz/logs/jpsix_alcareco/presetB/$(ClusterId).$(ProcId).{out,err,log}`
  - `queue filename from raw_2016H_list.txt`
- [x] 3.2 Wrote `condor/jpsix_alcareco/condor_jpsix_alcareco_presetC.sub`: identical except preset letter in `arguments` and output subpath in log directives.
- [x] 3.3 Created log dirs: `/work/submit/pmlugato/mz/logs/jpsix_alcareco/preset{B,C}`.

**Also written ahead of schedule (originally task 6.1) for completeness of the wrapper-set:**
- [x] 6.1 (stub) `condor/jpsix_alcareco/merge_and_report.py` — I/O scaffolding, JSON aggregation, hadd, markdown summary all functional. `_emit_plots()` is a TODO marker until we have real outputs to debug the plotting against.

## 4. Smoke test — 1 file × 1 preset

- [x] 4.1 Smoke test passed on cluster 2998602 (preset B, single RAW file from run 283270 — the largest 3.8 GB file). Ran on `c0710a-s8.ufhpc` (T2_US_Florida). Producer healthy.
  - **Runtime**: 5133 s wall (~85 min), `cmsRun` 2403 s (~40 min) for 10 000 events ≈ 0.24 s/event. `scram b` 49 s. Peak RSS 4.0 GB. Output 31 MB.
  - **Event yield**: 4 457 events written out of 10 000 input (44.6% event-level pass).
  - **Candidate counts** (verified offline via FWLite — sidecar counter was broken in this run, fixed for next): bplus 4 947 / bzero_kstar 4 013 / bc 8 161 / bs_phi 377 / b0_ks 135 / lambdab 26 / psi2s 101. All 7 channels populated, yields consistent with Phase-2 5 000-event interactive numbers (~1.1–1.8 cands/event on B+/B0/Bc, sparser on Bs and V0 channels — expected from physics).
- [x] 4.2 Five iterations needed before the smoke passed. Bugs found and fixed (worker-side artefacts of going from interactive to Condor):
  - **2998489**: tarball untarred inside `src/` → `src/src/...` double-nest, `scram b` did nothing meaningful. Fix: untar at CMSSW root.
  - **2998490, 2998499**: still failing, but `cmsrun.err` was hidden inside worker scratch. Fix: wrapper now surfaces tail of `cmsrun.err` and `cmsrun.out` to Condor stderr on non-zero exit.
  - **2998499**: `AttributeError: 'Process' object has no attribute 'OutALCARECOTkAlJpsiX_noDrop'` — modified `Configuration/StandardSequences/AlCaRecoStreams_cff.py` and `Configuration/EventContent/AlCaRecoOutput_cff.py` were not in the tarball; the cmsDriver-generated config refers to the JpsiX output module which only exists when the central stream registry is overridden. Fix: add the modified python files to the tarball.
  - **2998500**: `ImportError: No module named Services_cff` — when shipping ONE python file in `src/Configuration/StandardSequences/python/`, the override directory becomes a regular Python package (has `__init__.py` after `scram b`) and **shadows the release-base entirely** — sibling files like `Services_cff.py` become invisible. Fix: ship the entire `python/` dir of every Configuration sub-package whose `__init__.py` we override (~1.4 MB total across StandardSequences, EventContent, AlCa).
  - **Wrapper counter bug found post-smoke** (not blocking; offline-verified): the per-channel candidate counter used PyROOT TTree::Draw, but `edm::Wrapper<vector<...>>::size()` isn't reachable that way. Fix: replaced with FWLite `Events`/`Handle` iteration. The `n_out` regex was also broken; fix: parse `edmFileUtil`'s `... N events ...` line with `grep -oE`.

## 5. Real submission — 5 files × 2 presets

**Filelist amendment before submission**: the original 5-file selection included `/283/305/.../04581E7A...root` at 1.1 MB. Probing the rest of `283/305/00000/` confirmed every file in that subdir is ~1.1 MB — the run was very short (probably a few thousand events total). Swapped that entry for a second file from `283/283`, so the filelist still represents 3 different runs (283270, 283283, 283682) but with all 5 files in the 3.4–3.8 GB range, matching the smoke-validated input.

- [x] 5.1 `condor_submit condor_jpsix_alcareco_presetB.sub` → 5 jobs queued as cluster **3031959.{0-4}**.
- [x] 5.2 `condor_submit condor_jpsix_alcareco_presetC.sub` → 5 jobs queued as cluster **3031961.{0-4}**.
- [x] 5.3 All 10 jobs landed in history with exit 0. Wall-clock per job 21–33 min (DESY/RAL workers, much faster than the Florida smoke worker). Total real-time elapsed ~33 min from job-start to last-completion.

## 6. Aggregation and reporting

- [x] 6.1 `condor/jpsix_alcareco/merge_and_report.py` (functional aggregator) + `_plot_mass_overlays.py` (FWLite-based mass-overlay plotter) shipped. Quirk: both `hadd` (CMSSW slc7 binary) and FWLite (cmsenv) need the `cmssw-el7` Singularity wrapper because the submit machine is EL9. The aggregator handles this transparently — must be run from the repo root so the bind mount includes the CMSSW area.
- [x] 6.2 Ran from repo root. Outputs:
  - `condor/jpsix_alcareco/results_5files.md` — per-preset diagnostics + candidate-count tables
  - `~/public_html/mz/alcareco/condor_5files/mass_overlay_<channel>.png` × 7 — per-channel B-vs-C mass overlays
  - `~/public_html/mz/alcareco/condor_5files/overlay_summary.json` — mean/RMS per channel per preset
  - `condor/jpsix_alcareco/_merge_scratch/preset_{B,C}/*.{root,json}` — per-job pulled outputs
  - `condor/jpsix_alcareco/_merge_scratch/preset_{B,C}/jpsix_alcareco_preset{B,C}_Run2016H_merged.root` — hadd'd per-preset merged ROOTs

## 7. Update OpenSpec change with results

- [x] 7.1 Wrote `results.md` with operational summary, per-channel cand counts, mass-shape diagnostics, compliance vs proposed numbers, and next-step recommendation.
- [x] 7.2 The spec requirements as written survived contact with the smoke test. The two wrapper-contract bug fixes (cmsrun.err capture, FWLite-based candidate counter) are wrapper-internal — they don't alter what the spec specifies, just how the wrapper achieves it. No spec deltas needed.
- [x] 7.3 `openspec validate add-jpsi-x-condor-production --strict --no-interactive` — valid.

## 8. Update slides

- [x] 8.1 Replaced the speculative "Next step --- Condor scale-up" frame in `slides/jpsix-producer-progress.tex` with three frames: (a) Phase 3 operational summary with actual numbers, (b) per-channel candidate count + mass-shape qualifier table, (c) B$^+$ mass-overlay PNG. Updated the Summary frame to mark Phase 3 done and rephrase open questions with informed recommendations.
- [x] 8.2 Recompiled: 16 pages, 253 KB. Two minor overflow warnings remain (one cosmetic table-edge 18.7 pt, one negligible 0.7 pt vbox) — visually clean.

## 9. Decision on next scale

- [ ] 9.1 **Awaiting user decision.** Recommendation: (a) lock the cuts (iter-1 from sibling change holds at scale — mass shapes clean, V0-invariant preserved, no surprises); (b) BsPhi the one channel below threshold, fixable by 25× scale-up (~70 min wall-clock with the wrapper as-is, just edit `raw_2016H_list.txt` to ~125 entries) rather than channel-specific cuts. Sibling change `add-jpsi-x-vertex-fit-and-low-pt` can be archived once the user approves; this change (`add-jpsi-x-condor-production`) can be archived once the BsPhi follow-up is decided.
