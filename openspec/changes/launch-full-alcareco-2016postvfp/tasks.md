# Tasks

## 0. Prerequisites (before any submit)

- [ ] 0.1 Confirm postVFP boundary run (target ≥ 278820).
      Record as `condor/full_alcareco_2016postvfp/lfns/postvfp_run_boundary.txt`
      with a one-line provenance comment (path where it was read from
      in `wremnants/` or the CMS convention).
- [ ] 0.2 Ceph quota preflight: `df` + `du` at
      `/data/user/p/pmlugato/`. Confirm ≥ 2 TB free before launch
      (validation output already occupies ~30 GB; expected
      full-launch delta ≈ 1.4 TB).
- [ ] 0.3 Extract failing-site set from validation clusters
      3227852 / 3227853 condor logs into
      `condor/full_alcareco_2016postvfp/site_deny.txt`.
- [ ] 0.4 Rebuild populated-file filter on 2016H
      (`charmonium_2016H_populated.txt`,
      `singlemuon_2016H_populated.txt`) by symlinking or copying
      from `condor/multichannel_alcareco/lfns/`. No re-query.
- [ ] 0.5 Rehearsal test: submit 5 files each (Charmonium + SingleMuon)
      through the new `submit_batch.sh`, one file per era, to smoke
      the AAA-fallback + manifest + reconcile pipeline before
      committing to the rehearsal batch.

## 1. Filelist generation

- [ ] 1.1 DAS query for Charmonium 2016F, G, H:
      `dasgoclient --query "file dataset=/Charmonium/Run2016{F,G,H}-v1/RAW
      | grep file.name, file.nevents"` → `lfns/raw/charmonium_2016{F,G,H}_<TS>.tsv`.
- [ ] 1.2 DAS query for SingleMuon 2016F, G, H (same shape).
- [ ] 1.3 2016F run-range filter:
      `dasgoclient --query "run file=<LFN>"` per file, or the
      cheaper `dasgoclient --query "file dataset=... run=..."` for
      each run in the era; drop runs below the postVFP boundary.
      Emit `charmonium_2016F_postVFP.txt` and
      `singlemuon_2016F_postVFP.txt`.
- [ ] 1.4 Populated-file filter: `awk '$2 >= 1000'` on the raw TSVs,
      output `<pd>_<era>_populated.txt` (one LFN per line).
- [ ] 1.5 Sanity check: `wc -l` on all six populated lists; document
      counts in the runbook.
- [ ] 1.6 Commit `lfns/raw/*.tsv` (or add to `.gitignore` and stash
      elsewhere — the raw TSVs are ~10 MB total; keeping them in-repo
      is the reproducibility-friendly option).

## 2. Launch directory scaffolding

- [ ] 2.1 Create `condor/full_alcareco_2016postvfp/` with the
      structure from the proposal.
- [ ] 2.2 `build_tarball.sh`: mirror
      `condor/multichannel_alcareco/build_tarball.sh` shape, write
      `full_alcareco_payload.tgz`. Emit both recoskim cfgs freshly
      via `make_recoskim_cfgs.sh`.
- [ ] 2.3 `make_recoskim_cfgs.sh`: two `cmsDriver.py` invocations
      (one per PD). Era-independent — the recoskim cfg accepts any
      input file at cmsRun time via `sed`-patch inside `run.sh`.
- [ ] 2.4 `run.sh`: copy from validation, add stage-in xrdcp loop
      (`cms-xrd-global` → `cmsxrootd.fnal.gov` → `xrootd-cms.infn.it`)
      writing to `$_CONDOR_SCRATCH_DIR/<basename>.root`; patch
      recoskim cfg to read `file:./<basename>.root`; end-of-job
      trap to `rm -f` the local copy on any exit path;
      `edmFileUtil -f -e <local>` for run/lumi extraction
      (best-effort, non-fatal); flat 2-hour watchdog;
      `exit_reason` tag in JSON.
- [ ] 2.5 Six `.sub` files under `submits/`. Each interpolates
      `site_deny.txt` into a narrowed `DESIRED_Sites`. All share
      the payload tarball and `run.sh`. Bump `RequestDisk` to
      16 GiB (RAW ≤ ~5 GB + build tree + outputs + headroom).
- [ ] 2.6 `submit_batch.sh`: `<pd> <era> [--dry-run] [--limit N]
      [--list <filelist>]`. Skips inputs already recorded successful
      in `manifests/${pd}_${era}_success.txt`. Records cluster ID
      into `manifests/${pd}_${era}_clusters.log`.
- [ ] 2.7 `find_missing.py`: buckets unrun / failed_input_open /
      failed_cmsrun / failed_xrdcp / success; emits
      `resubmit_${pd}_${era}_<TS>.txt` covering only the auto-safe
      buckets.
- [ ] 2.8 `build_manifest.py`: run per drained cluster; produces
      `manifests/${pd}_${era}_cluster_${clusterId}.json`.
- [ ] 2.9 `merge_manifests.py`: era-level + PD-level roll-ups.
- [ ] 2.10 `README.md` runbook: submit → watch → reconcile →
      resubmit → close-out, with the exact commands.
- [ ] 2.11 `quota.sh`: `df -h /data/user/p/pmlugato/` + `du -sh`
      of the output tree. Preflight before every batch.

## 3. Rehearsal batch (200 files × 2 PDs × 3 eras optional per era)

- [ ] 3.1 Rebuild payload tarball (`./build_tarball.sh`).
- [ ] 3.2 Charmonium rehearsal: 200 files each across F-postVFP, G, H.
      Track failure rate per era.
- [ ] 3.3 SingleMuon rehearsal: 200 files each across F-postVFP, G, H.
- [ ] 3.4 Drain, reconcile, resubmit failures once. Confirm final
      success rate ≥ 95% per era.
- [ ] 3.5 Cross-check: mean runtime, output size, RSS per era match
      the 2016H validation numbers within a factor of 2. If not,
      flag before full launch.

## 4. Full launch

- [ ] 4.1 Submit Charmonium 2016F-postVFP (`submit_batch.sh
      charmonium 2016F_postVFP`).
- [ ] 4.2 Submit Charmonium 2016G.
- [ ] 4.3 Submit Charmonium 2016H (resume filter skips validation's
      120 already-successful outputs — no double compute).
- [ ] 4.4 Submit SingleMuon 2016F-postVFP.
- [ ] 4.5 Submit SingleMuon 2016G.
- [ ] 4.6 Submit SingleMuon 2016H.
- [ ] 4.7 Watch queue; drain per era; run `build_manifest.py` after
      each cluster.

## 5. Reconciliation cycle

- [ ] 5.1 For each (pd, era): `find_missing.py` → reconcile.
- [ ] 5.2 Resubmit `unrun ∪ failed_input_open ∪ failed_xrdcp` via
      `submit_batch.sh --list …`.
- [ ] 5.3 Manually triage `failed_cmsrun` bucket (grep stderr,
      identify pattern, decide keep/skip/patch).
- [ ] 5.4 Repeat until every (pd, era) shows 100% coverage or an
      accepted skip list is documented.

## 6. Close-out

- [ ] 6.1 `merge_manifests.py` — produce final era-level and
      PD-level manifests + a top-level `LAUNCH_MANIFEST.json`.
- [ ] 6.2 Final resource summary: total wall-days, total CPU-hrs,
      total output volume, failure rate breakdown.
- [ ] 6.3 Slides update: append a "full launch" section to
      `slides/preprod-multichannel-validation.tex` (or make a new
      `slides/full-launch-2016postvfp.tex` — decide during close-out).
- [ ] 6.4 Update project memory `[[project-jpsi-x-alcareco]]` with
      the final numbers and the fact that
      `condor/multichannel_alcareco/` is superseded.
- [ ] 6.5 `openspec archive launch-full-alcareco-2016postvfp`.

## 7. Approval gates

- [ ] 7.1 User approves this proposal before task §2 begins.
- [ ] 7.2 User approves the rehearsal batch results before task §4
      begins.
- [ ] 7.3 User approves close-out before archival.
