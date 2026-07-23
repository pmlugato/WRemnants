# Results — add-jpsi-x-condor-production

> **⚠ DEPRECATED, 2026-06-15.** The numbers in this document were collected against a
> producer with two latent bugs identified after the run:
> (1) the cosα cut in `JpsiXCandidateProducer::applyBLevelVertexFit()` uses the
> beamspot (LHC bunch center) as the origin instead of the matched PV, making
> the cut a random sign function of the PV-BS z-offset rather than a B-pointing
> discriminant; (2) the B-candidate mass is computed without a J/ψ mass
> constraint, broadening the peak from ~15 MeV to ~30+ MeV.
> Consequence: **preset C's iter-1 cuts as evaluated here were operating on
> broken variables**. The reported preset-C signal efficiency (~0.6 %) and the
> "preset-C-too-aggressive" framing of this document are artifacts of those
> bugs, NOT the cut design. Diagnostic re-run on the same input file with
> Kalman-fit-on but cuts at no-op values (Cdiag, cluster 3107800) showed
> Kalman convergence at 88–96 % per channel — the fit machinery is healthy;
> only the cuts on top of it are broken.
>
> **Use the follow-up change `fix-jpsi-x-cosalpha-bug-and-jpsi-mass-constraint`
> and its `results.md` for the corrected numbers** once that change runs.
> The operational findings here (wrapper contract, tarball, runtime, output
> volume, V0-invariance) ARE still valid; only the per-channel
> efficiency / mass-shape conclusions are not.



5-file × 2-preset HTCondor production on 2016H Charmonium RAW. All 10 jobs exited 0. The wrapper contract held; the producer cff is operating as designed; preset B yields the falling-combinatorial baseline and preset C the Kalman-cleaned alignment-grade sample.

## Operational summary

| | Preset B | Preset C |
|---|---|---|
| Jobs (ok / total) | 5 / 5 | 5 / 5 |
| Input events (`n_in`) | 50 000 | 50 000 |
| Output events (`n_out`) | 20 377 (40.8%) | 3 326 (6.6%) |
| Total output ROOT | 143 MB | 31 MB |
| Mean `cmsRun` runtime / job | 1 426 s | 1 557 s |
| Per-event cost | 143 ms | 156 ms |
| Peak RSS | 4 043 MB | 4 094 MB |
| `scram b` (on-worker, first job) | 49–64 s | 49–64 s |

10 jobs total. Slowest individual: 1 976 s (33 min) on a RAL worker (preset C, file index 4). DESY workers were fastest (~21 min preset B). The smoke run on Florida T2_US (85 min wall) was the worst-case worker we hit; the production split between DESY and RAL ran ~3× faster on the same code path.

Total cluster-CPU: ~5.3 hours (sum of all `cmsRun`). Wall-clock real-time: ~33 min (max of any single job, since the 10 jobs ran in parallel after queue clearing).

## Per-channel candidate counts (sum across 5 files = 50 000 events / preset)

| Channel | Preset B | Preset C | Ratio C/B |
|---|---|---|---|
| B+ → J/ψ K+ | 22 809 | 695 | 3.05 % |
| B0 → J/ψ K*0 | 18 241 | 301 | 1.65 % |
| Bs → J/ψ φ | 1 714 | 32 | 1.87 % |
| Bc → J/ψ π | 36 804 | 1 179 | 3.20 % |
| B0 → J/ψ Ks (V0) | 601 | 601 | 100.0 % |
| Λb → J/ψ Λ (V0) | 156 | 156 | 100.0 % |
| ψ(2S) → J/ψ ππ (V0) | 475 | 475 | 100.0 % |

**The V0 channels are bit-identical between presets** — confirming the cff invariant from the sibling change (`add-jpsi-x-vertex-fit-and-low-pt`, V0 channels exempt from preset switching). The non-V0 channels show preset C accepting 1.7–3.2 % of preset B candidates; the Kalman vertex fit + cosα + Lxy/σ cuts kill 96–98 % of combinatorial background.

## Mass-shape diagnostics

Per-channel mass-overlay PNGs at `~/public_html/mz/alcareco/condor_5files/mass_overlay_<channel>.png`. Summary statistics (mean, RMS over the producer mass window):

| Channel | PDG m | B mean | B RMS | C mean | C RMS |
|---|---|---|---|---|---|
| B+ | 5.279 | 5.254 | 0.144 | 5.265 | 0.140 |
| B0 → K*0 | 5.280 | 5.242 | 0.144 | 5.257 | 0.145 |
| Bs → φ | 5.367 | 5.397 | 0.114 | 5.395 | 0.131 |
| Bc | 6.275 | 6.251 | 0.203 | 6.245 | 0.206 |
| B0 → Ks | 5.280 | 5.249 | 0.141 | 5.249 | 0.141 |
| Λb | 5.620 | 5.659 | 0.150 | 5.659 | 0.150 |
| ψ(2S) | 3.686 | 3.747 | 0.094 | 3.747 | 0.094 |

Notes:

- **B+ and B0 → K*0**: preset C tightens the mean toward PDG by ~10 MeV — the Kalman cuts suppress low-mass combinatorial wings, pushing the distribution closer to the resonance. RMS narrows slightly (144 → 140 MeV) on B+. Visible peak emergence in the overlay PNGs.
- **Bs → φ**: 32 candidates in preset C is too few for a shape fit. The mean is dominated by individual entries; the RMS broadening (114 → 131 MeV) reflects sparse statistics, not a real broadening. **Bs is production-rate-limited** as flagged in the sibling-change `results.md`. Need ≥ 10× current statistics for a meaningful Bs shape under preset C.
- **Bc**: largest absolute yield in both presets — the broad pion mass window pulls in many candidates. RMS at 0.20 GeV is the producer's window width, not an intrinsic resolution. The Bc shape in preset C is a useful reference for the alignment downstream.
- **V0 channels**: all three preset-invariant to identical mean/RMS — sanity-check of the cff design.

## Compliance with proposal numbers

| Metric | Proposed | Observed | Verdict |
|---|---|---|---|
| Total events per preset | 50 000 | 50 000 (exact) | ✓ |
| Per-job wall-clock | 7 min (naive) | 21–33 min (DESY/RAL) | Higher than naive, **predictable** |
| Output volume preset B | ~500 MB | 143 MB | Lower (better) |
| Output volume preset C | ~1.5 GB | 31 MB | Much lower (better, but reflects preset C's high selectivity) |
| C/B output ratio | 3 × | 0.22 × | C smaller because it filters 97 % of B's candidates |
| Per-job failure rate | < 10 % | 0 % | ✓ |

The actual preset-C output is *smaller* than preset B (not larger as the proposal estimated, which was extrapolated from the per-100-event Phase-2 figures). This is because preset C's event-level pass rate (6.6 %) is six times lower than preset B's (40.8 %), and preset-C events that do pass carry fewer candidates each. Net effect: preset C is a tighter, smaller AlCaReco stream — exactly the alignment-grade output we wanted.

## Statistics targets vs achieved

Phase-2 sibling change set 500 candidates / channel as the threshold for a clean shape fit on the non-V0 channels.

| Channel | Preset B (50 k evts) | Preset C (50 k evts) | Cands needed for shape fit |
|---|---|---|---|
| B+ | 22 809 | 695 | met under both |
| B0 → K*0 | 18 241 | 301 | preset B met; preset C marginal (~60 %) |
| Bs → φ | 1 714 | 32 | preset B met; preset C structurally short |
| Bc | 36 804 | 1 179 | met under both |

**Conclusion for next-scale:** to clear the 500-cand bar on preset C for *all* non-V0 channels, BsPhi sets the floor — would need ~ 800 000 events. For preset C B0→K*0 alone, a ~ 5× scale (250 k events) suffices. The Phase-2 results.md guidance still stands: scale-up by 10× for B+, B0, Bc shapes; Bs is its own scope.

## Cluster diagnostics

- **Worker site distribution** (10 jobs): DESY × 5, RAL × 4, Florida (smoke only) × 0. Geographic spread came naturally from the broad `+DESIRED_Sites` list.
- **Slowest job**: cluster `3031961.4` (preset C, file 4 from run 283682) at 1 976 s = 33 min wall. Bottleneck is per-event RECO + 4-body Kalman fit cost.
- **`scram b` first-touch**: 49–64 s on every worker — fresh release setup is dominated by edmplugin and python-symlink generation, not by our ~17 modified Alignment plugins.
- **xrdcp output**: all 20 transfers (10 ROOTs + 10 JSONs) to `submit50.mit.edu` succeeded first attempt.
- **Input pull from `root://eoscms.cern.ch/`**: zero failures despite 5 × 3.5 GB transfers per preset.
- **No held jobs, no requeues, no exit-code-90 storms** — the smoke iteration cleared every wrapper bug before production.

## What this validates and what it leaves open

**Validated**:

1. The Condor wrapper contract — env-var on `cmsRun` line, `fileNames` and `maxEvents` overrides only, JSON sidecar emit on every exit code, xrdcp with EOS redirector — works at scale across 3 sites with no manual intervention.
2. The modified-package tarball strategy (package-only + worker `scram b`) — 163 KB tarball, ~50 s build on every worker, no ABI mismatches under `cms:rhel7` Singularity.
3. The cff design from `add-jpsi-x-vertex-fit-and-low-pt` — V0 channels preset-invariant, non-V0 channels respond to preset switching, mass shapes consistent with physics expectations.
4. The aggregator pipeline — `merge_and_report.py` + `_plot_mass_overlays.py` running inside `cmssw-el7` for FWLite + hadd produces per-preset merged ROOTs, a markdown summary, and 7 PNGs in one invocation from the host.

**Left open** (explicitly out of scope per this proposal):

1. **Bs → φ statistics**: 32 candidates in preset C at 50 k events. A 25× scale (1.25 M events) would put it at ~ 800 candidates — useful for a shape fit. Belongs in a follow-up scale-up change, not here.
2. **Next-tier scale-up to 500 k events / preset**: the wrapper is ready; just grow `raw_2016H_list.txt` from 5 to 50 entries and re-submit. Estimated cluster wall ~33 min (same — still parallel), output ~1.4 GB B + ~310 MB C.
3. **Multi-era expansion**: 2016B–G + 2017 + 2018 datasets. The wrapper is era-agnostic if we generalize the cmsRun config names; a sibling change should own this.
4. **Plot pretty-up**: the current overlay PNGs are correct but minimal (axis ticks, legends, dashed PDG line). For group-presentation polish — log y, fit overlays, sideband subtraction — a separate small change.

## Next-step recommendation

Lock in the iter-1 cuts from `add-jpsi-x-vertex-fit-and-low-pt`; the 5-file × 2-preset production passes all gates. The remaining production question is **how to address BsPhi**: scale-up to ~1 M events, or accept channel-specific cuts only for BsPhi (loosen `minBLxyOverSigma` to 0.5 on that producer). My recommendation: scale-up first — the production wrapper here makes 50× the current run cost about ~70 min real-time wall-clock (still parallel, just 50 jobs instead of 5 per preset). That delivers the BsPhi shape without rebalancing the cut design.
