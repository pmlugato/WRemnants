# Change: Compare two AlCaReco-stage selection levels (A vs B) for TkAlJpsiX and pick the production preset

## Why

The TkAlJpsiX stream currently runs with Phase 1+2 kinematic + DCA cuts
(`add-jpsi-x-candidate-quality-cuts`). The David-Walter design review
(`david_260506_alcareco.pdf`, slides 18–28) frames the V0/D* AlCaReco
path as **AlCaReco production only** at Stage 1 — any kinematic vertex
fitting is deferred to Stage 2 (CVH refit, CMSSW_10_6_26). For TkAlJpsiX
we adopt the same constraint: **no Kalman vertex fitting at the AlCaReco
stage** (any new vertex fit is too expensive for a stream that runs
event-by-event on the full parent dataset).

Validation so far (19,619 events, Run2016H Charmonium; jpsix-alcareco-producer
slides 13–15) shows the V0-mode channels clean at ~1 cand/event while the
non-V0 channels remain flat at O(5–100) cand/event. The pending
`add-jpsi-x-stage2-cvh-refit` proposal hypothesised that a B-level Kalman
vertex fit would collapse those distributions; that work is now
explicitly **deferred to Stage 2** (out of scope here).

This change runs a quantitative comparison of two AlCaReco-stage
selection levels for the four non-V0 channels and uses the result to
pick the production preset. Channels split into two groups by which
selection question is asked:

### V0-mode channels (B0→Ks, Λb, ψ(2S)) — **decision already structurally fixed**

The Ks and Λ sub-resonances are reconstructed by the V0Producer, which
runs as part of standard RECO (we consume it; we do not add a vertex
fit). V0Producer applies a Kalman vertex fit to opposite-sign track
pairs, a flight significance > 15σ cut, a pointing-angle cut, and a
post-fit mass window (±70 MeV for Ks, ±50 MeV for Λ). On the X side of
the channel this gives an already-clean candidate at ~1 / event.

The `JpsiXCandidateProducer` instances for the three V0-mode channels
stay minimal — only the channel-defining B mass window plus J/ψ pT and
αBS-at-J/ψ cuts shared across the stream. No A/B comparison is run for
these channels; they are fixed at "V0 quality + minimal kinematic" and
serve as the within-run reference against which the non-V0 results are
interpreted.

### Non-V0 channels (B+, B0→K*0, Bs→φ, Bc) — **A vs B is the question**

These channels have **no vertex-quality handle at Stage 1** — they are
purely kinematic. The validation run (slide 15) shows them flat at
5–105 cand/event. Since vertex fitting is off the table at this stage,
the only available levers are kinematic + geometric (pT, η, αBS, DCA,
and the channel-defining mass windows). The question is *how much*
kinematic+geometric tightening buys, and how much candidate
multiplicity remains for Stage 2 to clean up.

| Preset | Definition (non-V0 channels only) |
|---|---|
| **A — raw mass only** (reference baseline) | Catalogue mass windows only: B mass, sub-resonance mass (K*0, φ), J/ψ mass. No pT, η, αBS, or DCA cuts. HLT-floor bachelor pT (0.5/0.3 GeV) retained because it is a trigger property, not a selection refinement. |
| **B — A + kinematic + geometric** | A plus every kinematic / geometric handle available without a Kalman fit: J/ψ pT > 3, daughter \|η\| < 2.4, mother pT > 5 (3 for ψ(2S)), αBS-at-J/ψ < 1.0, **bachelor-to-J/ψ-vertex DCA < 1.0 cm** (the explicit geometric substitute for B-vertex consistency). This is what `add-jpsi-x-candidate-quality-cuts` calls Phase 1+2. |

(A previously-defined "Preset C" would have added a B-level Kalman
vertex fit and turned on K*0/φ sub-resonance Kalman fits. That entire
direction is **dropped from the AlCaReco stage** per the cost
constraint. Any vertex-fit work belongs to Stage 2.)

**Decision rule**:
- V0-mode channels (B0→Ks, Λb, ψ(2S)) → fixed at "V0 quality + minimal kinematic" (no preset).
- Non-V0 channels (B+, B0→K*0, Bs→φ, Bc) → choose between A and B.
  Production candidate is B unless A meets the gates with comparable
  output cost (a fallback that is unlikely in practice; the gate
  evaluation will quantify it).

**Decision criterion** (two-phase; see Studies section for the split):

- **Phase 1 (MC, B+ only)** — gates against gen truth:
  - signal efficiency on gen-matched B+ candidates ≥ 70% relative to preset A;
  - per-event fake rate (non-gen-matched candidates surviving the preset) reduced by ≥ 5× from A to B;
  - mass-distribution width on gen-matched candidates is preserved (no migration outside the catalogue B mass window).
- **Phase 2 (data + RAW ALCARECO, all channels)** — original three operating-point gates:
  - cands/event ≤ 5 in the channel's mass window;
  - fraction of mass-window candidates inside the tight signal sub-window (B ± 50 MeV around PDG, ψ(2S) ± 30 MeV) ≥ 30%;
  - tight-window absolute yield ≥ 100 candidates per 1000 events.

Channels that fail the Phase-2 gates under B are flagged in `results.md`
as "insufficient at Stage 1 alone; relies on Stage-2 cleanup." This is
the explicit hand-off to David's Stage-2 pipeline. No further AlCaReco-
stage cuts are introduced to fix them.

## What changes

### 1. Disable K*0/φ sub-resonance Kalman vertex fits already wired in the cff

The current `ALCARECOTkAlJpsiX_cff.py` has `applyVertexFit=True,
minVtxProb=0.01` on both the K*0 and φ `TwoBodyDecayCandidateProducer`
instances. These SHALL be turned off (`applyVertexFit=False`, parameter
removed or set to `0`). The K*0 and φ sub-resonances are then
discriminated by their mass window alone at the AlCaReco stage.

### 2. Selection preset switching in `ALCARECOTkAlJpsiX_cff.py`

A single boolean-block or env-var-driven switch at the top of the cff:

```python
_TKALJPSIX_SELECTION_PRESET = os.environ.get('TKALJPSIX_SELECTION_PRESET', 'B')
assert _TKALJPSIX_SELECTION_PRESET in ('A', 'B')
```

The four non-V0 `JpsiXCandidateProducer` instances SHALL set their cut
parameters from a small `dict`-of-`dict`s keyed by preset, while the
three V0-mode instances SHALL be preset-invariant. Mass windows are the
same across both presets (they define the channels). The only
differences are pT/η/αBS/DCA (B vs A).

**No C++ changes are required for this change.** The
`JpsiXCandidateProducer` already supports all the kinematic + DCA
parameters needed for B (they are currently in use). The proposal
specifically does NOT add `minBVtxProb`, `minBLxyOverSigma`, or any
TransientTrackBuilder ESConsumes.

### 3. Studies — two-phase measurement

#### Execution environment (applies to all Phase-1 invocations)

Phase 1 commands SHALL be issued from `/work/submit/pmlugato/REALmz/mz/real/WRemnants` inside the project's Singularity container. Workflow:

1. From `/work/submit/pmlugato/REALmz/mz/real/WRemnants`, run `sing` — opens the container shell with the WRemnants environment.
2. Inside the container, invoke the histmaker with `python scripts/histmakers/btojpsik.py …`.

**Nominal (analysis-path) invocation** — the template Phase 1 will adapt:

```bash
python scripts/histmakers/btojpsik.py \
    --allaxes \
    --era 2018 \
    --dataPath '/scratch/submit/cms/zmass/' \
    --filterProcs BuToJpsiK 'Charmonium_2018A' 'Charmonium_2018B' 'Charmonium_2018C' 'Charmonium_2018D' \
    -p 'E_for_e_m32nomc-allaxes' \
    -o '/ceph/submit/data/user/p/pmlugato/mz/histograms/' \
    --muonScaleVariation smearingWeightsSplines \
    --muonCorrData 'lbl_massfit' \
    --muonCorrMC 'idealMC_lbltruth_massfit' \
    --includeKaonScaleVariations \
    --theoryCorr none \
    --fitPtQuantiles 8 \
    --nomc \
    --histToFit m32
```

**Phase 1 adaptations** (resolved 2026-06-10):

- `--era 2016PostVFP` — matches the `RunIISummer20UL16` MC vintage; corrections/calibrations align with the MC reconstruction era. (Drop `--fitPtQuantiles`: with AlCaReco selections, the quantile-axis builder receives too few candidates and fails; the fit-observable axis is irrelevant for Phase-1 metrics.)
- `--filterProcs BuToJpsiK` only — drop all `Charmonium_*` data procs.
- Drop `--nomc` — Phase-1 gen-truth metrics require the MC events to actually be processed.
- Add `--alcarecoPreset {A,B}` (the new flag this change introduces to `btojpsik.py`).
- `-p` tag is **iteration-stamped** so prior iterations are not overwritten: `jpsix-preset{A,B}-iter{N}` (N starting at 1). Each iteration's tag SHALL be cited explicitly in the worklog deck so the histogram tree is reachable from the frame.
- `-o '/ceph/submit/data/user/p/pmlugato/mz/histograms/'` retained — histograms land alongside the analysis-path output for side-by-side checks.

**Build-graph routing — disjoint routes**: `btojpsik.py` SHALL be extended so that `--alcarecoPreset {A,B}` routes the build-graph through `get_bkmm_alcareco_selections(preset)`; when the flag is unset the analysis path runs unchanged (bit-for-bit identical to the current behaviour). The two routes do NOT emit cut-flows in the same pass — one invocation produces one selection profile.

Each Phase-1 iteration SHALL log its exact command line in the worklog deck so iterations can be re-run verbatim.

#### AlCaReco-emulation variable policy (mandatory for Phase 1)

Phase 1 uses NanoAOD as a proxy for what the real AlCaReco-stage
`JpsiXCandidateProducer` would see. The NanoAOD branches fall into
four classes; only the first two are admissible for Phase-1 cuts.

| Class | Status | Examples |
|---|---|---|
| 1. Raw reco-track kinematics | **OK** | `Muon_pt/eta/phi/charge`, `mm_mu{1,2}_pt/eta/phi`, `bkmm_kaon_pt/eta/phi/charge` |
| 2. Track-to-track / track-to-beamspot geometric (no Kalman vertex fit) | **OK** | `bkmm_kaon_mu1_doca`, `bkmm_kaon_mu2_doca`, `bkmm_kaon_dxy_bs`, `bkmm_kaon_sdxy_bs` |
| 3. Dimuon Kalman vertex fit (J/ψ vertex) | **FORBIDDEN** | `mm_kin_mass`, `mm_kin_pt`, `mm_kin_alphaBS`, `mm_kin_vtx_prob`, `mm_kin_sl3d` |
| 4. B-vertex Kalman fit / J/ψ-mass-constrained refit | **FORBIDDEN** | `bkmm_jpsimc_*` (all), `bkmm_nomc_*` (all), `bkmm_bmm_bdt`, `bkmm_bmm_mva` |

Class 4 is forbidden because the AlCaReco stage does not perform a
B-vertex Kalman fit (cost constraint, Section 1). Any cut on a
class-4 variable would not be reproducible at the real AlCaReco
stage. Iteration 1 (2026-06-10) inadvertently cut on
`bkmm_jpsimc_mass` and `bkmm_jpsimc_kaon1eta` — its mass-shape and
fake-rate numbers are **superseded** and must be re-run with raw
variables. The result is preserved in the worklog as a documented
dead end (per the worklog cadence rule).

Class 3 is **also forbidden** (resolved 2026-06-10): the dimuon
Kalman vertex fit lives in `mm_kin_*` and was performed by the BMM
Tools NanoAOD producer at the mini-to-nano step. Although in CMSSW
the dimuon vertex fit by `JpsiMuMuCandidateProducer` is upstream of
AlCaReco, the strict policy for Phase 1 is to validate selections
using only quantities that don't depend on any vertex fit performed
in the mini-to-nano pipeline. Cuts that the real
`JpsiXCandidateProducer` would derive from a J/ψ vertex (α_BS@J/ψ,
J/ψ vertex probability, SL3D) are not validated in Phase 1. If
they are needed for the production preset they are introduced in
the C++ producer directly — Phase 1 will not have measured them on
this MC.

##### BMM Tools upstream pre-filter on `bkmm_*` (from source)

`Bmm5/NanoAOD/plugins/DileptonPlusXProducer.cc::buildLLXCandidates`
imposes the following pre-filter when writing `bkmm_*` branches.
Phase-1 sees its candidates only AFTER these cuts (no way around it
without re-running NanoAOD production):

| Quantity | Upstream BMM Tools pre-cut |
|---|---|
| Muon track quality | highPurity inner track |
| Muon pT | > 1.0 GeV (`MuonMinPt`) |
| Muon \|η\| | < 2.4 (`MuonMaxEta`) |
| Kaon track quality | highPurity, charge ≠ 0, `hasTrackDetails` |
| Kaon pdgId | \|pdgId\| == 211 (PF charged hadron) |
| Kaon pT | > 1.0 GeV (`HadronMinPt`) |
| Kaon \|η\| | < 2.4 (`HadronMaxEta`) |
| Kaon-muon overlap | rejected (PF cleaning) |
| Kaon-to-muon DOCA (both legs) | < 0.10 cm (`maxTwoTrackDOCA = 0.1`) |
| B+ candidate mass m(μμK) | ∈ [4.0, 6.0] GeV (`minBKllMass=4.0`, `maxBKllMass=6.0`) |

**Implication for preset B choices**: any Phase-1 cut at or looser
than the upstream floor is a no-op. The genuine handles available
beyond the BMM Tools floor are: tighter kaon pT (above 1.0 GeV),
tighter kaon \|η\| (below 2.4), tighter kaon-mu DOCA (below 0.10
cm), tighter B+ mass window (inside [4.0, 6.0]), tighter J/ψ mass
window, plus any cut on quantities BMM Tools doesn't pre-filter
(raw J/ψ pT, raw B+ pT, tighter muon pT than the HLT floor of 4
GeV). Iteration 1 saw most preset-B cuts as no-ops precisely because
they were at-or-below this floor.

**Bias acknowledgment**: the bkmm candidates we measure are a
subset of what real AlCaReco's `JpsiXCandidateProducer` would build
(AlCaReco iterates all `generalTracks` with its own quality cut, not
PF-cleaned charged hadrons with the BMM Tools pre-filter). The
A→B fake-rate **ratio** is robust to this — the same upstream
filter applies to both presets — but absolute within-event fake
counts in Phase 1 are a lower bound on what AlCaReco sees on data.
Absolute rates are measured in Phase 2 on RAW.

`IsoTrack_*` is not a usable alternative: empirically, that
collection has pT > 5 GeV + miniPFRelIso < 0.2 applied at NanoAOD
production, leaving zero kaon-like candidates in 50k B+ MC events.
Custom NanoAOD with the full track collection would be required to
remove the BMM Tools pre-filter bias — out of scope.

**Mandatory derived quantities** computed from raw 4-vectors via
RDataFrame `Define` (in `btojpsik_selections.py`):

- `raw_mumu_mass[i]` — m(μ⁺μ⁻) from `mm_mu{1,2}_pt/eta/phi` and the
  PDG muon mass 0.10566 GeV. Used for the J/ψ mass window.
- `raw_b_mass[i]` — m(μ⁺μ⁻K⁺) from the three raw 4-vectors and
  m_K = 0.49368 GeV. Used for the B-mass window.
- `raw_b_pt[i]` — pₜ of the raw 4-vector sum. Used for any "mother
  pT" cut.
- `raw_mumu_pt[i]` — pₜ of the raw dimuon sum. Used for the J/ψ pT
  cut (replacing `mm_kin_pt`, which is class-3 / forbidden).

A new `select_raw_*` family of selection functions in
`btojpsik_selections.py` cuts on these derived columns and on the
class-1/class-2 branches above. The existing `select_*` family
(which cuts on class-3 and class-4 branches) is unchanged and remains
the analysis path.

**Revised preset definitions (Phase 1, B+ channel)**:

Values chosen to be **meaningfully beyond the BMM Tools upstream
floor** (otherwise the cut is a no-op):

| Cut | A | B | BMM Tools floor |
|---|---|---|---|
| Opposite-sign dimuon | ✓ | ✓ | — |
| J/ψ mass on `raw_mumu_mass` (2.95–3.25 GeV) | ✓ | ✓ | none |
| B mass on `raw_b_mass` (5.0–5.5 GeV) | ✓ | ✓ | [4.0, 6.0] |
| `Muon_pt > 4` | — | ✓ | > 1.0 (HLT > 4) |
| `bkmm_kaon_pt > 1.5` | — | ✓ | > 1.0 |
| `|bkmm_kaon_eta| < 1.8` | — | ✓ | < 2.4 |
| `raw_mumu_pt > 3` (J/ψ pT) | — | ✓ | none |
| `raw_b_pt > 5` (mother pT) | — | ✓ | none |
| `bkmm_kaon_mu{1,2}_doca < 0.03` cm | — | ✓ | < 0.10 (both legs) |

`|Muon_eta|` is not in preset B (no-op vs the BMM `MuonMaxEta=2.4`
floor + HLT-trigger fiducial). Future iterations may add a stricter
`|Muon_eta| < 1.4` if useful.

**Preset B locked (2026-06-10) — iter 2 confirms**: iter 1b ran with
`|bkmm_kaon_eta| < 1.4` and yielded 68.6% signal efficiency
(1.4% below the 70% gate) and 6.77× fake-rate reduction (above the
5× gate). The dominant signal-killer was the kaon `|η|` cut (21%
loss). Iter 2 loosened to `|bkmm_kaon_eta| < 1.8` — locked in the
table above — and **measured 78.1% signal efficiency, 5.73× fake
reduction, 97.0% purity, 23.3 MeV mass-peak RMS**. All three Phase-1
gates met. Locked as the production preset B for Phase 2 ALCARECO
confirmation.

##### Per-channel preset-B translation (B+ → other non-V0)

The locked B+ values translate to the three remaining non-V0 channels
by topology and Q-value scaling. The shared J/ψ → μμ leg keeps its
B+ values verbatim; the bachelor side scales by available
phase-space; the mother-pT floor is a no-op on signal MC for all
heavy-B channels.

| Cut | B+ (locked) | B0→J/ψK*0(Kπ) | Bs→J/ψφ(KK) | Bc→J/ψπ | Physics motivation |
|---|---|---|---|---|---|
| Topology | track | VCC (K*0) | VCC (φ) | track | — |
| Sub-resonance mass | n/a | 0.80–0.99 (K*0 ±2Γ) | 0.99–1.04 (φ near-threshold) | n/a | channel-defining; preset-invariant |
| B mass window (GeV) | 5.0–5.5 | 5.0–5.5 | 5.2–5.6 | 5.9–6.6 | channel-defining; preset-invariant |
| `Muon_pt` | > 4 | > 4 | > 4 | > 4 | shared J/ψ leg; HLT-floor anyway |
| `minJpsiPt` | > 3 | > 3 | > 3 | > 3 | shared J/ψ leg |
| `minMotherPt` | > 5 | > 5 | > 5 | > 5 | floor for heavy-B selection; signal-MC no-op on all four |
| Bachelor / daughter pT | K: > 1.5 | K, π each > 1.0 | K1, K2 each > 1.0 | π: > 1.5 | one-prong gets full Q (B+: 1.69 GeV, Bc: 3.03 GeV); two-prong VCC daughters share resonance momentum, individual pT ~1/2 the one-prong scale |
| Bachelor / daughter \|η\| | < 1.8 | < 1.8 each | < 1.8 each | < 1.8 | tracker geometry; same for all charged daughters |
| Bachelor–muon DOCA (cm) | < 0.03 (both legs) | < 0.03 each daughter, each muon | < 0.03 each daughter, each muon | < 0.03 (both legs) | Phase-1-tuned for B+; same geometric scale on the other channels (B-daughter tracks come from the displaced B-vertex region, similar distance to dimuon tracks) |
| α_BS-at-J/ψ vertex | (n/a Phase-1) | (n/a Phase-1) | (n/a Phase-1) | (n/a Phase-1) | class-3 in Phase 1; at AlCaReco-C++ stage the cff already applies `maxMotherAlphaBS=0.3` from upstream RECO — keep current cff value across A and B; not preset-tuned |

**Topology notes**:
- **Track-mode** (B+, Bc): bachelor is a single `generalTracks` track.
  The `JpsiXCandidateProducer` already exposes `minBachelorPt`,
  `maxBachelorEta`; preset switching applies to these directly. A
  new parameter **`maxBachelorMuTrackDOCA`** SHALL be added to the
  `JpsiXCandidateProducer` to expose the Phase-1-tuned track-track
  DOCA (`bkmm_kaon_mu{1,2}_doca`) at the C++ stage — implemented by
  computing DOCA via `kalmanVertexFitter`-free track-track impact
  parameter (`TrajectoryStateClosestToPointBuilder` or
  `TwoTrackMinimumDistance`), rejecting candidates where either
  bachelor–muon pair fails. Reusing the existing
  `maxBachelorIPToJpsiVertex` (bachelor IP to the J/ψ vertex; a
  different geometric quantity, ~10× larger scale) was considered
  and rejected: the two quantities differ in physics, and Phase 1
  measured the track-track quantity specifically.
- **VCC-mode** (B0→K*0, Bs→φ): bachelor is a pre-built sub-resonance
  VCC produced by `TwoBodyDecayCandidateProducer`. Two C++ additions:
  - `TwoBodyDecayCandidateProducer` SHALL gain a **`minDaughterPt`**
    parameter (currently only `maxDaughterEta` exists). Floor cut
    applied per daughter at sub-resonance build time.
  - `JpsiXCandidateProducer` (VCC mode) SHALL **loop over the VCC's
    daughter tracks** when applying the bachelor–muon DOCA cut,
    rather than treating the VCC as a single composite. Each
    daughter track is checked against each of the two J/ψ muons
    using the same DOCA computation as track mode; candidate is
    rejected if any (daughter, muon) DOCA exceeds the threshold.
    Chosen over delegating to `TwoBodyDecayCandidateProducer`
    because the cut is intrinsically J/ψ-context-dependent (different
    dimuon → different muon tracks) while the sub-resonance build
    is J/ψ-independent — applying the cut in the JpsiX producer
    keeps DOCA logic centralized and avoids coupling
    `TwoBodyDecayCandidateProducer` to the dimuon collection.

**V0-mode channels** (B0→Ks, Λb, ψ(2S)) stay **preset-invariant**.
They ride on the central `V0Producer` (consumed, not added at
AlCaReco), which provides Kalman-vertex-fitted V0 candidates with
flight significance and pointing-angle cuts already applied. The
existing cff values for these instances (`minJpsiPt=3`,
`minMotherPt=5` for B0→Ks/Λb or `=3` for ψ(2S), `maxMotherAlphaBS=0.3`)
are kept across A and B.

**Forbidden cff parameters** (proposal: no Kalman vertex fitting at
AlCaReco): `minBVtxProb`, `minBLxyOverSigma` on every non-V0
`JpsiXCandidateProducer` instance must be disabled (either removed
or set to a sentinel that turns the cut off). `applyVertexFit=True`
on the K*0 and φ `TwoBodyDecayCandidateProducer` instances must be
turned to `False`. These edits are independent of the preset
switch (they apply unconditionally).

α_BS-at-J/ψ is **not present** in preset B because it requires a
J/ψ vertex (class 3, forbidden). If the real AlCaReco production
needs it, it is added in the C++ `JpsiXCandidateProducer` directly,
outside the Phase-1 validation scope.

These are the iteration-1b starting values; later iterations may
tighten further per the iteration plan in the worklog.

#### Phase 1 — Selection-physics validation on B+ MC (NanoAOD)

**Input**: the only currently-available NanoAOD MC,
`/scratch/submit/cms/zmass/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL16MiniAOD-106X_mcRun2_asymptotic_v13-v1+MINIAODSIM`.
Signal-only (every event is a generated B+ → J/ψ K+).

**Tool**: extend `scripts/histmakers/btojpsik.py`. Add a sibling selections
function `get_bkmm_alcareco_selections(preset)` (returning the same
`(label, lambda)` tuple list as the existing `get_bkmm_selections()`)
that encodes the A and B preset cut profiles. The existing analysis-
level selections SHALL remain unchanged. Wire a new argparse flag
`--alcarecoPreset {A,B}` (default unset = analysis path) into the
script; when set, the build-graph function routes through the AlCaReco
selections instead of the analysis ones. Two histmaker invocations per
preset produce the cut-flow + mass + kinematic histograms.

**Gen matching**: gen-match each surviving B+ candidate by ΔR matching
of the three leaf tracks (μ+, μ-, K+) to the truth decay products from
the generated B+. A candidate is "true" if all three leaf tracks gen-
match; "fake" otherwise.

**Metrics** (Phase-1 gates):
- **Signal efficiency**: surviving truth-matched candidates / truth-matched
  candidates that survived preset A. Required ≥ 70% under B.
- **Fake-rate reduction**: (non-truth-matched candidates per event under A)
  / (non-truth-matched candidates per event under B). Required ≥ 5×.
- **Mass-shape preservation**: width and mean of the truth-matched
  candidate mass distribution under B match those under A within
  statistical error.

**Caveat (documented in `results.md`)**: signal-only MC underestimates
data combinatorial background — the fake-rate metric measures within-B+
event random-track combinatorics only, which is a *subset* of what real
Charmonium data contains. Phase-1 results are sufficient to calibrate
cut tightness on B+ but cannot predict absolute cands/event on data;
Phase 2 closes that gap.

**Extrapolation to the other six channels**: B+ is the only channel
with available MC. Once Phase 1 settles the B+ preset, the same cut
philosophy SHALL be translated to the other non-V0 channels (B0→K*0,
Bs→φ, Bc) on physics grounds — same J/ψ pT, daughter η, mother pT
floors; channel-appropriate DCA, αBS, mass windows. Each translation
SHALL be logged in `results.md` with the per-channel justification.

#### Phase 2 — ALCARECO confirmation on RAW (deferred)

Once Phase 1 picks the B+ preset and the cuts are translated to the
other channels, **one** cmsRun ALCARECO pass per preset on Run2016H
Charmonium RAW
(`/store/data/Run2016H/Charmonium/RAW/v1/000/281/727/00000/22261BD6-F984-E611-AA9A-FA163EC18366.root`)
at `-n 1000` per preset. Covers all 7 channels (4 non-V0 under preset
switch; 3 V0-mode preset-invariant). Three sub-studies on the resulting
ALCARECO files:

- **Phase 2a — Mass distributions on data**: per-channel data mass
  distributions under A and B. Apply the original three Phase-2 gates
  (cands/event ≤ 5; tight-window fraction ≥ 30%; tight-window yield ≥
  100 / 1000 events). Channels that fail under B → flagged for Stage-2.
- **Phase 2b — Output size + compute time**: wall-clock per event
  (FastTimerService), output bytes per event, per-branch size from
  `edmEventSize`. Confirms the architectural cost of the chosen preset.
- **Phase 2c — Track deduplication**: per-event metrics — count of
  cloned tracks in merged collection; sum of unique-leaf-tracks per
  channel if each stream had its own AlignmentTrackSelector; compression
  factor. Validates the single-combined-stream architecture (independent
  of preset choice).

**Diagnostic tool for Phase 2**: a sibling `jpsix_preset_compare.py`
(top-level script, FWLite, reading two ALCARECO files) implements 2a +
2b + 2c. Smaller scope than the original plan because the
selection-physics question is already settled by Phase 1. Extends or
forks the existing `jpsix_diagnostics.py` (already in
`CMSSW_10_6_17_patch1/`).

##### Phase 2 execution environment

ALCARECO production needs the SL7 CMSSW container, NOT the same
Singularity image used for Phase 1. Workflow:

```bash
# 1. From any shell:
cmsset                                # source cmsset_default.sh
cmssw-el7                             # enters SL7 Singularity container
cmsset                                # re-source inside the container
cd CMSSW_10_6_17_patch1               # the working area (NOT _src_backup)
cmsenv                                # set CMSSW env vars
scram b -j8                           # build (picks up preset-table cff)

# 2. Generate the config (one per preset; the cff reads TKALJPSIX_SELECTION_PRESET)
TKALJPSIX_SELECTION_PRESET=A cmsDriver.py RECO \
  -s RAW2DIGI,L1Reco,RECO,ALCA:TkAlJpsiX \
  --runUnscheduled --nThreads 16 \
  --data --era Run2_2016 --scenario pp \
  --conditions 106X_dataRun2_v35 \
  --eventcontent ALCARECO --datatier ALCARECO \
  --customise Configuration/DataProcessing/RecoTLR.customisePostEra_Run2_2016 \
  --filein /store/data/Run2016H/Charmonium/RAW/v1/000/281/727/00000/22261BD6-F984-E611-AA9A-FA163EC18366.root \
  --python_filename recoskim_Run2016H_Charmonium_JpsiX_presetA.py \
  --no_exec -n 1000

# Same for preset B (env var = B, --python_filename presetB.py)

# 3. Run
TKALJPSIX_SELECTION_PRESET=A cmsRun recoskim_Run2016H_Charmonium_JpsiX_presetA.py
TKALJPSIX_SELECTION_PRESET=B cmsRun recoskim_Run2016H_Charmonium_JpsiX_presetB.py

# 4. Diagnostic + comparison
python jpsix_preset_compare.py \
  --presetA-file <ALCARECO_presetA.root> \
  --presetB-file <ALCARECO_presetB.root>
```

The env var `TKALJPSIX_SELECTION_PRESET` is read at the top of the
preset-table-driven `ALCARECOTkAlJpsiX_cff.py` and selects between A
and B preset values for the four non-V0 `JpsiXCandidateProducer`
instances. V0-mode instances ignore the env var (preset-invariant).

### 4. Results artifact

A `results.md` SHALL be added to this change directory holding the
final tables of (channel × preset) for: cands/event, mass-window-pass
fraction, events-with-candidate, output bytes/event, wall-time/event,
dedup compression factor. The same headline tables (cands/event and
bytes/event per preset, per channel) SHALL be mirrored into the main
public-facing slides `jpsix-alcareco-producer.tex` as a new frame
between "Validation run" and any successor frames, **only after** the
Phase-1 cuts are locked and Phase-2 numbers are in.

### 5. Documentation cadence — continuous slide-as-you-go

A second slide deck, `slides/jpsix-selection-comparison-worklog.tex`
(beamer, same theme as the existing deck), SHALL be created at the
start of this change and **grown frame-by-frame as work happens**, not
written at the end. The cadence rule:

- Every iteration that produces a histogram or a number → at least one
  new frame.
- Every "tried X, didn't work, switched to Y" → one frame with the
  before/after numbers.
- Every decision or conclusion → one frame stating it explicitly with
  the supporting evidence from earlier frames cited.

The worklog deck is the scratchpad / lab notebook; the main
`jpsix-alcareco-producer.tex` deck is the polished public summary.
Findings flow worklog → main only after they stabilize (typically at a
phase boundary). The worklog accumulates the full iteration history so
later readers (and future-you) can reconstruct *why* a cut value was
chosen, not just *what* it ended up being.

Structure of the worklog deck (frames added as work proceeds, not
pre-committed):

- **Setup** — input dataset, software versions, what's being asked.
- **Phase 1 — iteration N** — for each iteration: cuts tried, signal
  efficiency / fake-rate / mass-shape numbers, one mass-plot screenshot,
  one-sentence observation.
- **Phase 1 lock** — final B+ preset values + summary table.
- **B+ → other-channel translation** — per-channel cut translation
  table with physics justification per row.
- **Phase 2a / 2b / 2c** — data masses, size+time, dedup; one frame
  per sub-phase (more if anomalies surface).
- **Decision** — per-channel preset choice, gates cited.
- **Followups** — channels flagged for Stage-2 cleanup; MC-sample
  requests; loose ends.

Slide updates SHALL be committed alongside the code or measurement
change that produced them, in the same git commit when practical. A
commit that says "ran Phase 1 iteration 3" without a corresponding
worklog frame is incomplete by the rules of this change.

### 6. Decision and rollout

After Study 1+2+3 are complete and `results.md` is filled, the final
task is to **set the per-channel default preset** in the cff. Each
non-V0 channel SHALL be assigned A or B based on the three-gate
decision criterion. V0-mode channels remain at V0-quality (no
preset switching). The env-var override SHALL be retained as a
developer escape; the cff default literal SHALL match the decision
recorded in `results.md`.

## Impact

- Affected specs: `alcareco-jpsi-x` (modifies and adds requirements).
- Affected code:
  - `Alignment/CommonAlignmentProducer/python/ALCARECOTkAlJpsiX_cff.py`
    (turn off K*0/φ Kalman vertex fits; refactor to preset-table-driven
    cuts for the four non-V0 channels).
  - `scripts/histmakers/btojpsik.py` (add `get_bkmm_alcareco_selections(preset)`
    function alongside the existing `get_bkmm_selections()`; add
    `--alcarecoPreset {A,B}` argparse flag; route through build-graph
    when the flag is set). Existing analysis-level selections SHALL
    remain unchanged.
  - `CMSSW_10_6_17_patch1/jpsix_preset_compare.py` (new sibling of
    `jpsix_diagnostics.py`, for Phase 2 only).
  - New `slides/jpsix-selection-comparison-worklog.tex` (worklog deck,
    grown frame-by-frame as iterations happen).
  - New frames in `slides/jpsix-alcareco-producer.tex` (headline
    summary, populated only after Phase-2 numbers are in).
  - New `openspec/changes/add-jpsi-x-selection-comparison/results.md`.
- **No C++ changes.** All AlCaReco-stage vertex-fit work is out of scope.
- **No supersession of `add-jpsi-x-stage2-cvh-refit`.** That change
  remains pending and should be re-scoped as Stage-2 work (out of scope
  here) before any further action on it.

## Resolved questions (2026-05-13, updated 2026-06-10)

Original gating questions resolved with the user. Decisions:

1. **Decision rule** — per-channel, with structural defaults. V0-mode channels (B0→Ks, Λb, ψ(2S)) inherit V0Producer quality; no preset comparison. Non-V0 channels compare A vs B.
2. **K*0/φ sub-resonance vertex fit** — **off in both presets** (no AlCaReco-stage Kalman vertex fitting; cost rules it out). Sub-resonance discrimination = mass window only.
3. **MC samples for non-B+ channels** — proceed data-only; gen-vs-reco on B+ MC is a cross-check. Missing samples documented in `results.md` as a caveat.
4. **Event budget** — pilot 1000 events × 2 presets, scale to 5000 only if needed.
5. **Decision criterion** — three combined gates (cands/event ≤ 5; tight-window fraction ≥ 30%; tight-window yield ≥ 100 / 1000 events). Smallest passing preset wins; A preferred over B on ties (but B is the production default if any non-V0 channel fails A's gates).
6. **Plot publication quality** — analysis-level only.
7. **Dedup counterfactual** — per-stream `AlignmentTrackSelector` as the realistic alternative architecture; per-channel leaf-track union reported as upper-bound sanity check.
8. **Selection A bachelor pT** — keep the HLT-floor values (0.5 / 0.3 GeV). HLT floor is a trigger property, not a selection refinement.
9. **C-grid scan** — n/a (preset C dropped).
10. **Stage-1 / Stage-2 coupling** — independent for the preset decision. Stage 2 reconstructs candidates afresh from cloned tracks.

## Remaining open items (non-blocking)

- Locate or request central MC for B0→K*0, Bs→φ, Λb, ψ(2S), Bc, B0→Ks if downstream needs gen-matched signal yields.
- Re-scope `add-jpsi-x-stage2-cvh-refit` as Stage-2 work (the AlCaReco-stage framing it currently uses is obsolete; the actual vertex fitting belongs to David's Stage-2 CVH pipeline).
- If B leaves any non-V0 channel above the gate thresholds (specifically B0→K*0, where the current data shows 105 cand/event with Phase 1+2), document the residual combinatorics in `results.md` and flag for Stage-2 cleanup. Do not attempt to re-introduce vertex fitting at Stage 1 to fix it.
