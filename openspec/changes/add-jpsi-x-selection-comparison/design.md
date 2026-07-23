# Design: TkAlJpsiX AlCaReco-stage selection comparison (A vs B)

## Context

Two plausible AlCaReco-stage selection levels for the four non-V0
TkAlJpsiX channels:

- **A** — keep only the mass windows that define each channel.
- **B** — add kinematic + geometric cuts (J/ψ pT, daughter η, mother pT,
  αBS, bachelor-to-J/ψ-vertex DCA), i.e. the Phase 1+2 cuts in
  `add-jpsi-x-candidate-quality-cuts`.

Kalman vertex fitting at the AlCaReco stage (the original "preset C"
and the K*0/φ sub-resonance vertex fits currently wired in the cff)
**is dropped**. Cost rules it out for a stream that runs event-by-event
on the full parent dataset. Any vertex-based cleanup belongs to
David-Walter's Stage-2 CVH pipeline.

Validation to date (slides 13–15 of `jpsix-alcareco-producer.tex`)
shows B at 105 cand/event for B0→K*0 and 10 cand/event for Bs→φ — flat
distributions. Without vertex fitting we cannot collapse those numbers
further at Stage 1; the role of this study is to quantify exactly how
much candidate multiplicity B leaves and to confirm B's mass-shape
quality is enough for Stage 2 to operate on.

V0-mode channels (B0→Ks, Λb, ψ(2S)) ride on the V0Producer (which we
consume, not add). They are already ~1 cand/event with a clean mass
distribution. They are not preset-switched here.

## Goals / Non-Goals

**Goals**:
- One CMSSW build that can run both presets via a config switch.
- Three independent measurement campaigns (Studies 1, 2, 3) producing
  numeric tables comparing A and B.
- A documented decision and resulting default preset in the cff.

**Non-Goals**:
- Any new Kalman vertex fitting at the AlCaReco stage (cost).
- Stage-2 CVH refit changes (follow David's pipeline; out of scope).
- HLT pre-filter changes (already decided in
  `add-jpsi-x-alcareco-channels`: no HLT in default sequence).
- Tuning thresholds beyond a single point per preset.
- NanoAOD output format (independent of preset).

## Architectural decisions

| # | Decision | Choice | Rationale |
|---|---|---|---|
| 1 | Preset switching mechanism | Env-var-driven dict-of-dicts in cff | Single CMSSW build runs both; no recompile per preset |
| 2 | Vertex fitting at AlCaReco stage | **None** (drop preset C; turn off K*0/φ sub-resonance Kalman fits already wired in cff) | Cost: this stream runs event-by-event on the full parent dataset; any new Kalman fit is too expensive |
| 3 | Scope of A/B comparison | Non-V0 channels only (B+, B0→K*0, Bs→φ, Bc) | V0-mode channels inherit V0Producer quality (consumed, not added); no preset comparison needed |
| 4 | Selection A definition (non-V0) | Channel mass windows only; HLT-floor bachelor pT retained | Minimal-selection baseline; HLT floor is a trigger property, not a selection |
| 5 | Selection B definition (non-V0) | Mass + kinematic (J/ψ pT, daughter η, mother pT) + αBS + bachelor-to-J/ψ-vertex DCA | The current Phase 1+2 cuts in `add-jpsi-x-candidate-quality-cuts`; explicitly geometric — no Kalman fit |
| 6 | Mass windows | **Same across A and B** | Mass windows are part of the channel definition, not the selection level |
| 7 | Phase split | Phase 1 = NanoAOD MC on B+ (gen truth); Phase 2 = ALCARECO RAW (data, all channels, post-decision confirmation) | MC truth gives ground truth on signal-vs-fake without the RAW→RECO step; ALCARECO needed only to confirm size/timing/dedup |
| 8 | Phase 1 input format | NanoAOD (the format `btojpsik.py` already consumes) | Reuses existing narf/RDataFrame infrastructure; single MC sample available is in this format |
| 9 | Phase 1 tool | Extend `scripts/histmakers/btojpsik.py` with `get_bkmm_alcareco_selections(preset)` + `--alcarecoPreset` argparse flag; analysis-level selections untouched | Reuses validated cut-flow machinery; minimal new code; analysis path is preserved |
| 10 | Phase 1 scope | B+ channel only; other six channels' cuts translated by physics analogy and confirmed in Phase 2 | Only MC sample available is BuToJpsiK; other channels' MC is a future blocker for direct gen-truth validation |
| 11 | Phase 1 metrics | Signal efficiency on gen-matched candidates (≥ 70%); per-event fake-rate reduction A→B (≥ 5×); mass-shape preservation | Gen-truth-based; signal-only MC cannot give absolute data rates but rigorously answers the cut-tuning question |
| 12 | Phase 2 metrics | Original three gates: cands/event ≤ 5; tight-window fraction ≥ 30%; tight-window yield ≥ 100 / 1000 events | Data-level operating-point question; needs real Charmonium data to resolve |
| 13 | Phase 2 event budget | 1000 events × 2 presets on Run2016H Charmonium RAW | Confirmation pass only; selection physics already locked by Phase 1 |
| 14 | Phase 2 tool | Sibling `jpsix_preset_compare.py` (top-level FWLite script) | Reads two ALCARECO files; produces 2a/2b/2c tables; smaller than originally planned because Phase 1 has settled selection physics |
| 15 | Results persistence | `results.md` in this change dir + new frame in main slides deck | Durable spec artifact + polished public summary; separates Phase-1 and Phase-2 findings |
| 16 | **Slide-as-you-go cadence** | Sibling worklog deck `slides/jpsix-selection-comparison-worklog.tex`, grown frame-by-frame at every iteration / decision; main deck gets the polished summary only when findings stabilize | Continuous documentation lets reviewers reconstruct *why* a cut value was picked, not just *what*; commits without a corresponding worklog frame are incomplete by the rules of this change |
| 17 | Default after study | Per-channel preset assignment in cff | Cff is already per-channel; per-channel default is structurally free |
| 18 | Stage-1 / Stage-2 coupling | Independent for the preset decision | Stage 2 reconstructs candidates afresh from cloned tracks |

## TrackRef chain — unchanged

Neither preset touches the leaf `TrackRef` chain. The
`JpsiXCandidateProducer → CompositeDaughterTrackProducer →
AlignmentTrackSelectorWithIndexMapModule → VertexCompositeCandidateRemapper`
chain is the same one already validated in
`add-jpsi-x-alcareco-channels`. Removing the K*0/φ
`applyVertexFit=True` parameter likewise does not affect track
references — it only changes whether the candidate's vertex position
is the Kalman-fitted point or the track-midpoint, and whether
candidates failing the sub-vertex probability cut are dropped. With
the parameter off, all sub-vertex-pass-fail logic disappears and the
midpoint vertex is stored unconditionally.

## Preset-table-driven cff refactor sketch

```python
# In ALCARECOTkAlJpsiX_cff.py (top)
import os
_PRESET = os.environ.get('TKALJPSIX_SELECTION_PRESET', 'B')
assert _PRESET in ('A', 'B')

# Mass windows shared across both presets
_MASS = {
    'jpsi':    (2.95, 3.25),
    'kstar':   (0.80, 0.99),
    'phi':     (0.990, 1.040),
    'BPlus':   (5.0, 5.5),
    'B0Ks':    (5.0, 5.5),
    'B0Kstar': (5.0, 5.5),
    'BsPhi':   (5.2, 5.6),
    'Lambdab': (5.3, 6.0),
    'Psi2S':   (3.5, 3.9),
    'Bc':      (5.9, 6.6),
}

# Kinematic + geometric cuts for non-V0 channels
def _kin_non_v0(preset, channel):
    if preset == 'A':
        return dict(minJpsiPt=0., minMotherPt=0., maxBachelorEta=1e9,
                    maxBachelorIPToJpsiVertex=1e9, maxMotherAlphaBS=1e9)
    # preset == 'B'
    mother_pt = {'BPlus':5., 'B0Kstar':5., 'BsPhi':5., 'Bc':5.}[channel]
    return dict(minJpsiPt=3.0, minMotherPt=mother_pt, maxBachelorEta=2.4,
                maxBachelorIPToJpsiVertex=1.0, maxMotherAlphaBS=0.3)

# V0-mode channels stay minimal regardless of preset
def _kin_v0(channel):
    mother_pt = {'B0Ks':0., 'Lambdab':0., 'Psi2S':3.}[channel]
    return dict(minJpsiPt=3.0, minMotherPt=mother_pt,
                maxBachelorEta=1e9, maxBachelorIPToJpsiVertex=1e9,
                maxMotherAlphaBS=0.3)
```

K*0 and φ `TwoBodyDecayCandidateProducer` instances drop
`applyVertexFit` and `minVtxProb` entirely (or set
`applyVertexFit=False`). Sub-resonance discrimination is by mass
window alone.

## Risks / Trade-offs

- **Signal-only MC underestimates data combinatorial**. The Phase-1
  fake-rate metric measures within-B+-event random-track combinatorics
  only, not the full Charmonium background. Phase-1 rigorously answers
  the cut-tuning question (a cut that rejects within-event fakes will
  also reject the broader data fakes that include this subset) but
  cannot predict absolute cands/event on data. Phase 2 closes that gap.
- **B may leave B0→K*0 above the Phase-2 gate**. The validation run
  shows 105 cand/event with Phase 1+2. Removing the K*0 sub-resonance
  Kalman fit cleans that number *up* (more candidates pass), not down.
  If B0→K*0 fails the cands/event ≤ 5 gate even under B, document the
  residual and rely on Stage 2 to clean it. Do not re-introduce vertex
  fitting at Stage 1.
- **Cut translation B+ → other six channels is by physics analogy, not
  measurement**. The Phase-1 result drives B+ cut values directly;
  applying the same philosophy to other channels (same J/ψ pT, daughter
  η, mother pT, αBS floors; channel-appropriate DCA and mass windows)
  is an inference. Phase 2 confirms it operates correctly on data but
  cannot validate signal efficiency for those channels without their
  MC samples.
- **NanoAOD track-level info varies by channel**. The available
  BuToJpsiK NanoAOD has the leaf tracks needed for B+ reconstruction
  (the analysis already uses them). For other channels (K*0 → Kπ,
  φ → KK, etc.) the available NanoAOD tracks may not include both
  daughters at adequate pT. This is one reason the other channels
  defer to Phase-2 ALCARECO validation rather than Phase-1 MC.

## Migration plan (within this change)

1. Turn off K*0/φ Kalman vertex fits in the cff (`applyVertexFit=False`).
2. Refactor `ALCARECOTkAlJpsiX_cff.py` to the preset-table form (A/B).
3. **Phase 1** — extend `scripts/histmakers/btojpsik.py`: add
   `get_bkmm_alcareco_selections(preset)`, `--alcarecoPreset {A,B}`
   argparse flag, gen-match leaf tracks, produce histograms under
   each preset. Iterate cut tightness on B+ MC until Phase-1 gates
   (signal efficiency ≥ 70%, fake-rate reduction ≥ 5×, mass-shape
   preserved) are met.
4. Translate the locked B+ cuts to the other non-V0 channels by
   physics analogy; record each translation in `results.md`.
5. **Phase 2** — one cmsRun ALCARECO pass per preset on Run2016H RAW
   at `-n 1000`. Run Phase 2a (data masses), 2b (size+time), 2c (dedup)
   via the new `jpsix_preset_compare.py`.
6. Fill `results.md` with both phases' findings; mirror headline tables
   to the slides.
7. Pick the per-channel default preset; commit it.

## Open questions

Resolved (2026-05-13, updated 2026-06-10); see "Resolved questions" and
"Remaining open items" in `proposal.md` for the canonical record.
