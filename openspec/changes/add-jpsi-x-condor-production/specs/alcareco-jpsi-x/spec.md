## ADDED Requirements

### Requirement: Condor Production Wrapper for J/ψ+X ALCAReco

A wrapper script + Condor submit-file pair SHALL drive the `recoskim_Run2016H_Charmonium_JpsiX_preset{B,C}.py` cmsRun configs at scale on the MIT submit HTCondor cluster, without altering any producer-code or cff parameter. The wrapper contract is exhaustively the file-name override, the `maxEvents` override, the per-job preset-letter dispatch via the `TKALJPSIX_SELECTION_PRESET` environment variable, the output path layout, and the JSON sidecar schema. No other behavior is delegated to the wrapper — in particular, the wrapper SHALL NOT modify any cut value, producer parameter, or any line of the cff. Wrapper modifications that change candidate-emission behavior on identical input MUST be rejected at code review.

#### Scenario: One submit file per preset, one preset letter per job
- **WHEN** the operator submits `condor_jpsix_alcareco_preset<X>.sub` for `<X> ∈ {B, C}`
- **THEN** every job spawned by that submission SHALL run with `TKALJPSIX_SELECTION_PRESET=<X>` exported on the same shell line as the `cmsRun` invocation in `run_jpsix_alcareco.sh`, regardless of any environment configured in the `.sub` file's `environment` directive — this honors the sibling-change-documented bug that the cff is re-imported at cmsRun startup and only the cmsRun-time env var actually selects the preset

#### Scenario: Per-job RAW-input override of process.source.fileNames
- **WHEN** the worker script `run_jpsix_alcareco.sh` runs with argument `$2` set to a RAW LFN (`/store/data/Run2016H/Charmonium/RAW/v1/...`)
- **THEN** the wrapper SHALL replace the single hardcoded `fileNames = cms.untracked.vstring('...')` in `recoskim_Run2016H_Charmonium_JpsiX_preset<X>.py` with `cms.untracked.vstring('root://eoscms.cern.ch/$2')` AND append `process.maxEvents.input = cms.untracked.int32(10000)` to the config (idempotent override of the cmsDriver-default `maxEvents`). The `raw_2016H_list.txt` filelist SHALL contain bare LFNs starting `/store/data/Run2016H/Charmonium/RAW/v1/...` with no redirector prefix — the wrapper supplies the `root://eoscms.cern.ch/` prefix at config-patch time

#### Scenario: Output and JSON-sidecar layout
- **WHEN** a job completes (success OR `cmsRun` non-zero exit)
- **THEN** the wrapper SHALL `xrdcp` the AlCaReco ROOT output (if produced) to `root://submit50.mit.edu//data/user/p/pmlugato/mz/alcareco/jpsi-x/preset_<X>/Run2016H/jpsix_alcareco_preset<X>_<basename>.root`, AND a JSON sidecar to the same directory with extension `.json`, containing at minimum the keys `preset`, `input_lfn`, `n_in`, `n_out`, `n_cands_<channel>` for each of the seven channels, `runtime_s_cmsrun`, `runtime_s_total_job`, `peak_rss_mb`, `exit_code`, `host`, `scram_b_seconds`. On `cmsRun` non-zero exit the JSON sidecar SHALL still be written with `exit_code != 0` so the aggregator detects the failure

#### Scenario: CMSSW-package tarball, not full release
- **WHEN** the worker pulls `transfer_input_files = jpsix_alcareco_payload.tgz`
- **THEN** the tarball SHALL contain only `Alignment/CommonAlignmentProducer/{python,plugins,interface}/*` plus the two `recoskim_Run2016H_Charmonium_JpsiX_preset{B,C}.py` configs (target tarball size < 5 MB). The worker SHALL run `scram p CMSSW CMSSW_10_6_17_patch1` from a clean scratch directory, untar the payload over `src/`, run `scram b -j 8`, then `cmsRun`. Pre-built `.so` files MUST NOT be shipped in the tarball — workers compile fresh

#### Scenario: No-cut-change invariant
- **WHEN** any reviewer compares the wrapper's emitted cmsRun config (after `sed` patching) against the unmodified `recoskim_Run2016H_Charmonium_JpsiX_preset<X>.py`
- **THEN** the diff SHALL consist exclusively of two lines: the `process.source.fileNames` replacement and the appended `process.maxEvents.input` override. No `process.ALCARECOTkAlJpsiX*Candidates.*` parameters, no `process.maxEvents` other than `input`, no `process.options`, and no producer-module-level config lines shall be touched by the wrapper

#### Scenario: Smoke test before full submission
- **WHEN** the wrapper or `.sub` is first deployed (or modified)
- **THEN** the operator SHALL run a single-job test (preset B, single RAW file from `raw_2016H_list.txt`) and verify (a) the job runs to completion, (b) the AlCaReco ROOT arrives at the expected `/data/user/p/pmlugato/mz/alcareco/...` path, (c) the JSON sidecar arrives with `exit_code == 0` and `n_cands_*` consistent with Phase-2 interactive numbers (within ~2× for the same preset and event count, accounting for input-file luminosity variance), BEFORE running the full 10-job submission across both presets

### Requirement: Aggregation and Reporting for Condor-Produced Outputs

A standalone Python module `condor/jpsix_alcareco/merge_and_report.py` SHALL aggregate per-job outputs into a per-preset summary (one merged ROOT and one markdown report per preset) and emit mass-overlay plots reusing the existing `CMSSW_10_6_17_patch1/jpsix_preset_compare_3way.py` plotting helpers. The aggregator SHALL NOT redo the producer's work — it consumes the wrapper's outputs only.

#### Scenario: Per-preset hadd
- **WHEN** all jobs for preset `<X>` complete
- **THEN** `merge_and_report.py` SHALL `hadd` every ROOT under `/data/user/p/pmlugato/mz/alcareco/jpsi-x/preset_<X>/Run2016H/jpsix_alcareco_preset<X>_*.root` into a single `jpsix_alcareco_preset<X>_Run2016H_merged.root`, preserving the EDM file structure

#### Scenario: Per-preset diagnostics table
- **WHEN** the aggregator runs against a per-preset set of JSON sidecars
- **THEN** it SHALL emit a markdown table with rows: total n_in, total n_out, total n_cands per channel, mean runtime_s_cmsrun per job, mean s/event, total output_size_mb, max peak_rss_mb, number of failed jobs (exit_code != 0). Output written to `condor/jpsix_alcareco/results_5files.md` (or analogous `results_<N>files.md` for next-scale runs)

#### Scenario: Mass-overlay plot reuse
- **WHEN** the aggregator emits per-channel mass-overlay plots
- **THEN** it SHALL import the existing plotting helpers from `CMSSW_10_6_17_patch1/jpsix_preset_compare_3way.py` (or call it as a subprocess with the merged ROOTs as input), NOT reimplement the histogramming logic in `merge_and_report.py`. Plots written to `~/public_html/mz/alcareco/condor_5files/` (or analogous path)
