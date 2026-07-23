# Tasks

## 0. Prerequisites (before scaffolding)

- [ ] 0.1 Record the locked reco tags into
      `condor/mc_inclusive_btojpsix_2016postvfp/manifests/reco_tags.txt`:
      GT `106X_mcRun2_asymptotic_v17`, beamspot
      `Realistic25ns13TeV2016Collision`, era `Run2_2016`, HLT step
      `HLT:@relval2016`.
- [ ] 0.2 Resolve the UL16 postVFP premix minbias library on DAS,
      cache the LFN list under
      `condor/mc_inclusive_btojpsix_2016postvfp/manifests/premix_library.txt`.
- [ ] 0.3 Ceph preflight: `df` + `du` at `/data/user/p/pmlugato/`.
      Confirm ≥ 3 TB free (phase-1 ALCARECO volume ≈ 1.9 TB at
      4×10⁸ generated, plus headroom for logs/manifests/resubmits).
- [ ] 0.4a Record `base_seed = 20260716` in
      `condor/mc_inclusive_btojpsix_2016postvfp/manifests/base_seed.txt`
      with a one-line provenance comment.
- [ ] 0.4 Sanity check: unpack `jpsix_alcareco_payload.tgz` from
      `condor/full_alcareco_2016postvfp/` and confirm the
      `Alignment/CommonAlignmentProducer/python/ALCARECOTkAlJpsiX_*`
      and `Configuration/AlCa/python/AlCaRecoStreams_cff.py` patches
      are present.

## 1. Patched fragment

- [ ] 1.1 Copy the McM fragment
      `Configuration/GenProduction/python/BPH-RunIISummer20UL16GEN-00017-fragment.py`
      into
      `condor/mc_inclusive_btojpsix_2016postvfp/fragments/BPH-RunIISummer20UL16GEN-00017-fragment.py`.
- [ ] 1.2 Add anti-Λ_b support:
      `Alias Myanti-Lambda_b0 anti-Lambda_b0`,
      `ChargeConj MyLambda_b0 Myanti-Lambda_b0`,
      `CDecay Myanti-Lambda_b0`, and add `Myanti-Lambda_b0` to
      `list_forced_decays`. The `Decay MyLambda_b0` block is
      NOT modified — anti-Λ_b inherits it via `CDecay`, and the
      physical BR ratios (J/ψ Λ vs ψ(2S) Λ) as shipped by McM are
      preserved. Λ→pπ vs Λ→nπ⁰ split follows the physical Lambda0
      natural BR (~64% pπ⁻); no forcing.
- [ ] 1.3 Fix the misleading `configurationMetadata.annotation` string
      (`'Jpsi->mumu (no  kin cuts on muons)'` → correct filter
      description).
- [ ] 1.4 Diff-audit: `diff` the patched fragment against the McM
      original; verify only §1.2 + §1.3 differ.
- [ ] 1.5 Smoke test: run the patched fragment through
      `btojpsix2016mcprod_parallel.sh` (bumped to `10_6_20_patch1`)
      with 40 jobs × 25 filtered ≈ 1000 accepted; re-run
      `classify_jpsi.py`. Confirm Λ_b fraction rises from ~2.3% to
      ~4–5% (particle + antiparticle), consistent with the fix.

## 2. Launch directory scaffolding

- [ ] 2.1 Create
      `condor/mc_inclusive_btojpsix_2016postvfp/` with the structure
      from the proposal.
- [ ] 2.2 `build_tarball.sh`: mirror the raw-data version. Bundles
      (a) `CMSSW_10_6_20_patch1/src/` from the raw-data payload,
      (b) the patched fragment placed at
      `CMSSW_10_6_20_patch1/src/Configuration/GenProduction/python/`,
      (c) the full-chain cfg produced by `make_fullchain_cfg.sh`.
      Emits `mc_inclusive_payload.tgz`.
- [ ] 2.3 `make_fullchain_cfg.sh`: single `cmsDriver.py` invocation
      producing GEN,SIM,DIGI,L1,DIGI2RAW,HLT,RAW2DIGI,L1Reco,RECO,PAT,ALCA
      chain. `--customise` includes (i) the standard
      `Utils.addMonitoring`, (ii) the `AlCaRecoStreams_cff` patch
      that wires the `ALCARECOTkAlJpsiX*Resonances` producers into
      the output, and (iii) a small `customise` snippet appending
      `keep PileupSummaryInfos_*_*_*` to the ALCARECO output
      commands (MC-only branch for downstream PU reweighting).
- [ ] 2.4 `run.sh`: seed injection (based on `clusterId × 100000 +
      procId`), cmsRun invocation, output xrdcp, per-job JSON with
      `seed_used`, `n_generated`, `n_accepted`, `filter_efficiency`,
      exit tags {`ok | cmsrun_crash | xrdcp_out_failed |
      wrapper_error | wrapper_timeout`}, flat 10 h watchdog.
- [ ] 2.5 `submits/mc_inclusive.sub`: single .sub template. Bump
      `RequestDisk` to 24 GiB (DIGI-Premix output is disk-heavy);
      `RequestMemory` 6 GB. Uses same submit50.mit.edu xrdcp
      redirector as raw-data.
- [ ] 2.6 `submit_batch.sh <N_JOBS> [--dry-run] [--list <procId list>]`:
      renders a queue of `N_JOBS` distinct `(clusterId, procId)`
      pairs. Skip pairs already recorded successful in
      `manifests/success.txt` (resume filter).
- [ ] 2.7 `find_missing.py`: expected work = union of scheduled
      `(clusterId, procId)` pairs across all clusters. Bucket:
      `success | unrun | failed_cmsrun | failed_xrdcp |
      failed_wrapper`. Emit `resubmit_<TS>.txt` containing only
      `unrun ∪ failed_xrdcp ∪ failed_wrapper` procId lists.
- [ ] 2.8 `build_manifest.py`: per-drained-cluster JSON manifest.
      Records `(clusterId, procId, seed_used, n_generated, n_accepted,
      output_size_mb, runtime_s, peak_rss_mb, exit_reason, host)`.
- [ ] 2.9 `merge_manifests.py`: cluster → run-level roll-up.
- [ ] 2.10 `quota.sh`: `df -h /data/user/p/pmlugato/` + `du -sh` of
       the MC output tree. Preflight before every batch.
- [ ] 2.11 One-line entry-point script `SUBMIT.sh` for the external
       operator: `./SUBMIT.sh rehearsal` (100 jobs) or
       `./SUBMIT.sh phase1` (20,000 jobs, 4×10⁸ generated). Runs
       quota preflight, calls `submit_batch.sh`, records the
       cluster ID. No `phase2` target — that's a separate future
       proposal.
- [ ] 2.12 `RECONCILE.sh`: wraps `find_missing.py` + emits the
       auto-resubmit list. Operator runs after each drain.
- [ ] 2.13 `HANDOFF_README.md`: operator-facing runbook. Assumes
       the operator has condor + xrdcp + cvmfs. Documents:
       (a) where to unpack the tarball, (b) `SUBMIT.sh rehearsal`
       first, (c) what a successful rehearsal reconcile looks
       like, (d) the go/no-go we (analysis side) will send them
       before full submit, (e) `SUBMIT.sh full` and drain, (f)
       reconcile cycle, (g) how to hand the outputs back
       (ceph path + manifest tarball).

## 3. Local validation (mandatory, runs on this machine)

- [ ] 3.1 Build the payload tarball (`./build_tarball.sh`).
- [ ] 3.2 Launch the local multi-parallel driver (adapt
      `btojpsix2016mcprod_parallel.sh`) with `K` single-threaded
      cmsRun processes running the full-chain cfg. Target ~1000
      accepted events total (~50k generated).
- [ ] 3.3 Validate ALCARECO event content on one output file:
      `edmDumpEventContent <file>` MUST list all branches in
      `ALCARECOTkAlJpsiX_Output_cff.py`.
- [ ] 3.4 Byte-diff the branch list against a raw-data ALCARECO
      output from `condor/full_alcareco_2016postvfp/`; confirm zero
      diff on the `ALCARECOTkAlJpsiX*` prefix.
- [ ] 3.5 Per-job scratch footprint: measure peak disk usage
      (`du -sm $CMSRUN_DIR` sampled during each run). Set
      `.sub` files' `RequestDisk` = p99 + 2× headroom.
- [ ] 3.6 Per-job RSS: measure peak RSS (from
      `FrameworkJobReport.xml` `PeakValueRss`). Set
      `.sub` files' `RequestMemory` = p99 + 1 GB headroom.
- [ ] 3.7 Wall-clock per accepted event: measure mean s/accepted.
      Compute N_gen_per_job that lands at ≤ 8 h wall (target ≤ 7 h).
      Bake into `submits/*.sub`.
- [ ] 3.8 Seed determinism: rerun the same
      `(clusterId=0, procId=k)` seed; byte-diff the ROOT payload
      branches (metadata timestamps excepted) — MUST be identical.
- [ ] 3.9 Anti-Λ_b sanity: rerun `classify_jpsi.py` on the merged
      local output; Λ_b species fraction rises from ~2.3% to
      ~4–5% (consistent with adding the anti-Λ_b decay table).
- [ ] 3.10 Package `local_validation_report.md` under the launch
       directory summarizing the six measurements above. Cite the
       exact `.sub` values chosen (RequestDisk, RequestMemory,
       N_gen_per_job, wall-time watchdog).

## 4. Handoff preparation

- [ ] 4.1 Tarball the launch directory into
      `mc_inclusive_handoff_<YYYY-MM-DD>.tgz` (includes the
      payload tarball, submit files, worker script, reconcile
      tooling, `SUBMIT.sh`, `RECONCILE.sh`, `HANDOFF_README.md`,
      and `local_validation_report.md`).
- [ ] 4.2 Confirm with the operator the target ceph destination
      is writable via `root://submit50.mit.edu/` from their
      submit host.
- [ ] 4.3 Send the tarball + a one-page cover note stating:
      (i) rehearsal-first, (ii) reconcile report format we
      expect back, (iii) go/no-go decision on our side before
      full submit.

## 5. Operator-side execution (informational; not run by us)

The following happens on the operator's side. No task in this
list is executed by us. We document expected outputs so that when
the operator's reconcile output comes back, we know how to read it.

- Rehearsal (`./SUBMIT.sh rehearsal`) → 100-job cluster; operator
  drains, runs `./RECONCILE.sh`, sends back the reconcile JSON.
- We review: (a) success rate ≥ 95% after one resubmit cycle,
  (b) sample output ROOT branch content matches our local
  reference, (c) no unexpected `failed_cmsrun` patterns.
- Go/no-go from us. If go: operator runs `./SUBMIT.sh full`.
- Operator sends the final reconcile + `LAUNCH_MANIFEST.json`.

## 6. Close-out

- [ ] 6.1 Verify final manifest: 100% coverage or a documented
      skip list.
- [ ] 6.2 Run `classify_jpsi.py` on a merged 10⁶-accepted-event
      subset of the operator-returned outputs. Record channel
      composition into `manifests/channel_composition.txt`.
- [ ] 6.3 Slides: add an "inclusive MC production complete"
      section.
- [ ] 6.4 Update project memory `[[project-jpsi-x-alcareco]]` with
      the final MC ALCARECO location and stats.
- [ ] 6.5 `openspec archive add-mc-inclusive-btojpsix-2016postvfp`.

## 7. Approval gates

- [ ] 7.1 User approves this proposal before §2 begins.
- [ ] 7.2 User approves local validation report (§3.10) before
      the handoff tarball ships (§4).
- [ ] 7.3 User approves rehearsal reconcile before we send
      go-signal to the operator for full submit.
- [ ] 7.4 User approves close-out before archival.
