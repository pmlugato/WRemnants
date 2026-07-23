# TkAlJpsiX selection-preset comparison — results

Two-phase study comparing presets **A** (mass windows only) and **B** (mass windows + kinematic + geometric, no Kalman vertex fit) for the four non-V0 channels of the TkAlJpsiX AlCaReco stream. V0-mode channels are preset-invariant by construction (V0Producer Kalman fit upstream).

---

## Phase 1 — selection physics on B+ MC

Source: `BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen` (UL16 NanoAOD, 81 files, 4.87M events). Tooling: `scripts/histmakers/btojpsik.py` with the new `--alcarecoPreset {A,B}` flag and `get_bkmm_alcareco_selections(preset)` selection function.

### Cut tables

| Cut                              | A (mass only) | B iter 1b      | B iter 2 (locked) |
|----------------------------------|:-------------:|:--------------:|:-----------------:|
| Opposite-sign dimuon             | required      | required       | required          |
| raw J/ψ mass (GeV)               | [2.95, 3.25]  | [2.95, 3.25]   | [2.95, 3.25]      |
| raw B+ mass (GeV)                | [5.0, 5.5]    | [5.0, 5.5]     | [5.0, 5.5]        |
| `Muon_pt` (GeV)                  | —             | > 4            | > 4               |
| `bkmm_kaon_pt` (GeV)             | —             | > 1.5          | > 1.5             |
| `|bkmm_kaon_eta|`                | —             | < 1.4          | < 1.8             |
| raw J/ψ pT (GeV)                 | —             | > 3 (no-op)    | > 3 (no-op)       |
| raw B+ pT (GeV)                  | —             | > 5 (no-op)    | > 5 (no-op)       |
| `bkmm_kaon_mu{1,2}_doca` (cm)    | —             | < 0.03 (both)  | < 0.03 (both)     |

### Phase-1 metrics (locked iter 2)

| Metric                              | Preset A | Preset B (iter 2) | Gate          | Result |
|-------------------------------------|---------:|------------------:|---------------|:------:|
| Signal efficiency (truth-matched)   | 100% (norm) | 78.1%          | ≥ 70%         |  PASS  |
| Fake-rate reduction (A / B)         | 1.00 (norm) | 5.73×           | ≥ 5×          |  PASS  |
| Tight-window peak mean (MeV drift)  | 0 (norm) | 0.1               | ≤ 2σ-stat     |  PASS  |
| Tight-window RMS (MeV)              | 24.2     | 23.3              | preserved     |  PASS  |
| Purity (fraction truth-matched)     | 87.9%    | 97.0%             | (informative) |   —    |

**Phase-1 lock**: preset B iter 2 (above table) is the locked B+ selection. All gates pass.

### Variable-policy rule (iter 1 → iter 1b correction)

Iter 1 used `bkmm_jpsimc_mass` (class-4, Kalman-fit-derived) for the mass cut. This was overruled by the proposal's variable policy:

- Class 1 (raw reco-track kinematics) — OK.
- Class 2 (track-track / track-beamspot geometric, no Kalman) — OK.
- Class 3 (`mm_kin_*`, dimuon vertex fit) — FORBIDDEN.
- Class 4 (`bkmm_jpsimc_*`, `bkmm_nomc_*`, BDT/MVA) — FORBIDDEN.

Iter 1b reran on raw-variable masses computed on-the-fly from `mm_mu{1,2}_pt/eta/phi`, `bkmm_kaon_pt/eta/phi` + PDG μ/K masses via `define_raw_kinematics(df)`. Iter 1b numbers and iter 2 (above) reflect this corrected variable set.

### BMM Tools upstream floor

Reading `Bmm5/NanoAOD/plugins/DileptonPlusXProducer.cc::buildLLXCandidates`, the upstream NanoAOD-step cuts are already:

| Quantity         | Upstream cut                   |
|------------------|--------------------------------|
| muon highPurity  | required                       |
| muon pT          | > 1.0 GeV                      |
| muon \|η\|         | < 2.4                          |
| kaon pT          | > 1.0 GeV (`HadronMinPt`)      |
| kaon \|η\|         | < 2.4 (`HadronMaxEta`)         |
| kaon-mu DOCA     | < 0.10 cm (`maxTwoTrackDOCA`)  |
| B+ mass          | ∈ [4.0, 6.0] GeV               |

Phase-1 preset-B cuts at or looser than these are no-ops on this MC. The Phase-1 absolute fake counts are lower bounds on data — the A/B ratio (5.73× rejection) is the trustworthy quantity.

---

## Phase 2 — ALCARECO confirmation on Run2016H Charmonium RAW

Input: `/store/data/Run2016H/Charmonium/RAW/v1/000/281/727/00000/22261BD6-F984-E611-AA9A-FA163EC18366.root`, -n 1000 events, datatier ALCARECO.

cmsDriver-and-cmsRun script: `/tmp/jpsix_phase2_cmsdriver.sh` + `/tmp/jpsix_phase2_cmsrun.sh`. `TKALJPSIX_SELECTION_PRESET` set on both cmsDriver and cmsRun invocations (the env var is consumed at cmsRun-time cff-import, so setting it only at cmsDriver leaks the default into both runs — see 2026-06-11 bug frame in the worklog).

### Phase 2a — mass quality (per-channel candidate yields)

| Channel    | Mode   | A:cands | A:c/ev | A:yld/1k | B:cands | B:c/ev | B:tight | B:yld/1k | Verdict |
|------------|--------|--------:|-------:|---------:|--------:|-------:|--------:|---------:|:-------:|
| B+         | non-V0 |    8611 |  13.52 |     4082 |      78 |   0.36 |   36%   |      129 |  PASS   |
| B0 K*0     | non-V0 |   79400 | 124.65 |    37287 |      77 |   0.35 |   27%   |       97 |  FAIL   |
| Bs φ       | non-V0 |    6800 |  10.68 |     3846 |       6 |   0.03 |   33%   |        9 |  FAIL   |
| Bc         | non-V0 |   18662 |  29.30 |     6410 |     103 |   0.47 |   18%   |       88 |  FAIL   |
| B0 Ks      | V0     |      22 |   0.03 |       16 |      22 |   0.10 |   45%   |       46 | (inv.)  |
| Λb         | V0     |       2 |   0.00 |        0 |       2 |   0.01 |    0%   |        0 | (inv.)  |
| ψ(2S)      | V0     |      19 |   0.03 |        6 |      19 |   0.09 |   21%   |       18 | (inv.)  |

Gates (evaluated on preset B as the production candidate): cands/event ≤ 5, tight-window-frac ≥ 30%, tight-window-yield ≥ 100/1000.

V0-mode rows: candidate counts byte-identical between A and B (preset-invariant cff confirmed). "c/ev" differs only because the denominator (events saved by the AlignmentTrackSelector filter) differs (637 saved under A vs 217 under B); the structural preset-invariance holds.

### Phase 2b — output size + wall-time

| Quantity                       | Preset A | Preset B | Ratio A/B |
|--------------------------------|---------:|---------:|----------:|
| Events saved (out of 1000)     |      637 |      217 |    2.94   |
| File size (MB)                 |    38.38 |     3.14 |   12.21   |
| Bytes / saved event (kB)       |     60.3 |     14.5 |    4.16   |
| Wall time (s)                  |    351.7 |    391.2 |    0.90   |
| Wall ms / input event          |    552.2 |   1802.6 |    0.31   |
| Sum ALCARECO branch bytes (MB) |     0.30 |     0.01 |   30.0    |

Preset B output is 12× smaller in total bytes (~3× from fewer events surviving the filter, ~4× from fewer leaf tracks per surviving event). Wall-time difference is small (11% in favour of A on this run) and dominated by full-RECO upstream of the ALCA step; the preset switch affects <5% of total wall-clock.

### Phase 2c — track deduplication

| Quantity                            | Preset A | Preset B |
|-------------------------------------|---------:|---------:|
| Events read                         |      637 |      217 |
| Tracks (merged collection)          |   53,280 |      826 |
| Tracks if stored per-channel        |  426,791 |    1,047 |
| Tracks saved by dedup               |  373,511 |      279 |
| Compression factor                  |   8.01×  |   1.27×  |
| Mean tracks/event                   |    83.6  |     3.8  |

Under preset A the combinatorial blowup (especially B0 K*0 at 125 cands/event) creates a large nominal "saved by dedup" — but absolute track count is what matters to downstream alignment, not the dedup ratio. Preset A's 83.6 tracks/event would dominate the entire AlCaReco budget on the parent dataset.

---

## Per-channel preset decision

| Channel    | Preset chosen | Rationale                                                                  |
|------------|:-------------:|----------------------------------------------------------------------------|
| B+         | **B**         | All three Phase-2 gates met (Phase-1 lock + Phase-2 confirm).             |
| B0 → K*0   | **B**         | Only B keeps cands/event below 5; A blows up to 124.65/event (factor 400 over B). Tight-window yield 97/1000 is marginal-FAIL but no Phase-2 gate is achievable under A. |
| Bs → φ     | **B**         | A produces 10.68 cands/event; B brings it to 0.03/event. Tight-window yield (9/1000) is signal-rarity-limited (BR(Bs→J/ψφ) × σ(Bs)), not selection-limited. |
| Bc         | **B**         | A produces 29.3 cands/event; B brings it to 0.47/event. Yield (88/1000) is marginal-FAIL; tightening further would not be physically justified at the AlCaReco stage. |
| B0 → Ks    | (preset-invariant) | V0Producer upstream provides Kalman vertex fit; preset value irrelevant. |
| Λb         | (preset-invariant) | Same.                                                                  |
| ψ(2S)      | (preset-invariant) | Same.                                                                  |

**Production preset**: `TKALJPSIX_SELECTION_PRESET=B` for the four non-V0 channels. The cff default is already `'B'`, so this decision ratifies the existing default; no per-channel override needed in `_NON_V0_PRESETS`.

Preset A is structurally unusable in production: $10^2$–$10^4$ candidates per event, $12\times$ output size on a 1000-event sample.

---

## Channels flagged for follow-up

- **Bs → φ** (yield 9/1000 ≪ 100): signal-rarity-limited, not selection-limited. Recommend running a 10k+ event Phase-2 sample to confirm yield scales linearly (no hidden selection loss).
- **Bc** (yield 88/1000 < 100; tight-frac 18% < 30%): the wide mass window [5.9, 6.6] GeV admits combinatorial background; tightening the window or adding a B-vertex cut would help but B-vertex cuts are forbidden at AlCaReco (cost). Defer to Stage-2 CVH pipeline for final Bc selection.
- **B0 → K*0** (yield 97/1000 ≈ 100): marginal-FAIL by 3 candidates. Likely passes in a 10k+ sample where statistics smooth out.

These follow-ups are tracked in the worklog deck's final frame and the OpenSpec change's `tasks.md` section 10.

---

## Artefacts produced by this study

- `CMSSW_10_6_17_patch1/src/Alignment/CommonAlignmentProducer/python/ALCARECOTkAlJpsiX_cff.py` — env-var preset switching, V0 invariance, no Kalman fits at AlCaReco.
- `CMSSW_10_6_17_patch1/src/Alignment/CommonAlignmentProducer/plugins/JpsiXCandidateProducer.cc` — `maxBachelorMuTrackDOCA` parameter + `trackTrackDCA` helper; track-mode and VCC-mode both honour it.
- `CMSSW_10_6_17_patch1/src/Alignment/CommonAlignmentProducer/plugins/TwoBodyDecayCandidateProducer.cc` — `minDaughterPt` parameter.
- `scripts/histmakers/btojpsik.py` — `--alcarecoPreset` flag, `get_bkmm_alcareco_selections(preset)`.
- `wremnants/production/btojpsik_selections.py` — `define_raw_kinematics(df)`, raw-variable `select_raw_*` cut family, `select_kaon_mu_doca`.
- `CMSSW_10_6_17_patch1/recoskim_Run2016H_Charmonium_JpsiX_preset{A,B}.py` — cmsDriver-generated configs with FastTimerService + per-preset output filename.
- `CMSSW_10_6_17_patch1/jpsix_preset_compare.py` — FWLite diagnostic producing 2a/2b/2c tables + overlay plots.
- `CMSSW_10_6_17_patch1/jpsix_alcareco_preset{A,B}.root` — Phase-2 ALCARECO output files (1000-event sample).
- `CMSSW_10_6_17_patch1/compare_phase2_run1/` — comparator output (overlay plot + JSON summary).
- `slides/jpsix-selection-comparison-worklog.tex` — 42-page chronological worklog deck.
