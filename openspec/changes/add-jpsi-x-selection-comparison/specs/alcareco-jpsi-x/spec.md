## ADDED Requirements

### Requirement: Selection-Preset Switching (non-V0 channels only)

`ALCARECOTkAlJpsiX_cff.py` SHALL expose a single selection-preset switch (`TKALJPSIX_SELECTION_PRESET` environment variable, with the values `'A'` or `'B'`; default `'B'`) that determines the cut parameters of the four non-V0-mode `JpsiXCandidateProducer` instances (B+, B0→K*0, Bs→φ, Bc). The three V0-mode `JpsiXCandidateProducer` instances (B0→Ks, Λb, ψ(2S)) SHALL have their cuts fixed at "V0 quality + minimal kinematic" (J/ψ pT, αBS shared with the rest of the stream; no DCA; no Kalman vertex-fit cuts) regardless of the preset value — the V0Producer's Kalman vertex fit on the Ks/Λ sub-resonance already provides clean candidates at ~1 / event without any further AlCaReco-stage selection. Mass windows for every channel SHALL NOT depend on the preset (mass windows define the channels). All non-mass cuts on non-V0 channels (J/ψ pT, mother pT, bachelor pT, bachelor η, bachelor-to-J/ψ-vertex DCA, αBS-at-J/ψ) SHALL be derived from one of two preset profiles defined in helper functions in the cff. A single CMSSW build SHALL be sufficient to run both presets without recompilation. The env var is consumed at cff-import time; since `ALCARECOTkAlJpsiX_cff.py` is transitively loaded by `Configuration.StandardSequences.AlCaRecoStreams_cff` inside the generated cmsRun config, the env var MUST be set both at cmsDriver invocation (so the dumped config's process_load resolves correctly under the chosen preset) AND at cmsRun invocation (since cmsRun re-imports the cff at job startup). Setting it only at cmsDriver time silently leaks the default preset into the cmsRun pass — a failure mode observed on 2026-06-11 and fixed by setting the env var on both invocation lines.

#### Scenario: Preset B is the default for non-V0 channels
- **WHEN** `cmsRun` is invoked without `TKALJPSIX_SELECTION_PRESET` set in the environment
- **THEN** the cuts on every non-V0 per-channel producer are identical to the Phase-1+2 values documented in `add-jpsi-x-candidate-quality-cuts`, augmented with the bachelor-to-J/ψ-vertex DCA cut

#### Scenario: Env var must be set on cmsRun (not just cmsDriver)
- **WHEN** the env var is set during cmsDriver config generation but not during the subsequent cmsRun invocation
- **THEN** the cmsRun job consumes the default preset because the cff is re-imported by the generated config's `process.load('Configuration.StandardSequences.AlCaRecoStreams_cff')` at cmsRun startup. Operators MUST set `TKALJPSIX_SELECTION_PRESET=<preset>` on both invocations of cmsDriver and cmsRun for a given preset run

#### Scenario: V0-mode channels are preset-invariant
- **WHEN** the preset is changed between A and B
- **THEN** the cut parameters of `ALCARECOTkAlJpsiXB0KsCandidates`, `ALCARECOTkAlJpsiXLambdabCandidates`, and `ALCARECOTkAlJpsiXPsi2SCandidates` remain bit-identical, and the candidates emitted by these three producers on the same input are bit-identical across the two preset settings

#### Scenario: Preset A disables non-mass cuts on non-V0 channels
- **WHEN** `TKALJPSIX_SELECTION_PRESET=A` is set
- **THEN** every kinematic cut (`minJpsiPt`, `minMotherPt`, `maxBachelorEta`, `maxBachelorIPToJpsiVertex`, `maxMotherAlphaBS`) on the four non-V0 `JpsiXCandidateProducer` instances is configured to its no-op value (0 or +∞ as appropriate). The bachelor-pT cut on B+ and Bc producers remains at the HLT-floor value (0.5 / 0.3 GeV) by design.

---

### Requirement: No Kalman Vertex Fitting at the AlCaReco Stage

The AlCaReco-stage producers in `ALCARECOTkAlJpsiX_cff.py` SHALL NOT invoke any Kalman vertex fit. Specifically, `applyVertexFit` on the K*0 and φ `TwoBodyDecayCandidateProducer` instances (`ALCARECOTkAlJpsiXKstarCandidates`, `ALCARECOTkAlJpsiXPhiCandidates`) SHALL be `False`, and `JpsiXCandidateProducer` SHALL NOT be configured with `minBVtxProb` or `minBLxyOverSigma` parameters. This is a cost decision: the AlCaReco stream runs event-by-event on the full parent dataset; any added Kalman fit is too expensive at this stage. Vertex-based cleanup belongs to Stage-2 CVH processing. The V0Producer's Kalman fit on Ks/Λ is exempt because it runs centrally in standard RECO (the AlCaReco stream consumes its output, not adds to it).

#### Scenario: K*0 sub-resonance vertex fit disabled
- **WHEN** an event passes the J/ψ + K*0 + B0 candidate chain
- **THEN** `ALCARECOTkAlJpsiXKstarCandidates` emits all opposite-sign track pairs whose K±π∓ mass falls in the K*0 mass window with NO vertex-probability cut applied, and the stored candidate vertex is the midpoint of the two track reference points (not a Kalman-fitted vertex)

#### Scenario: No B-level Kalman fit invoked
- **WHEN** any of the seven `JpsiXCandidateProducer` instances runs
- **THEN** it does NOT ESConsume `TransientTrackBuilder`, does NOT call `KalmanVertexFitter`, and the cuts on each candidate are exhaustively the kinematic + geometric cuts in the preset profile (mass windows, pT, η, αBS, DCA)

---

### Requirement: Phase-1 Selection-Physics Validation via Extended Histmaker

`scripts/histmakers/btojpsik.py` SHALL be extended with a sibling
selections function `get_bkmm_alcareco_selections(preset)` (with
`preset` in `('A', 'B')`) returning a `(label, lambda)` tuple-list of
the same shape as the existing `get_bkmm_selections()`. The existing
analysis-level selections function SHALL remain unchanged. A new
`--alcarecoPreset {A,B}` argparse flag (default unset) SHALL be added
to the script; when set, the build-graph function routes through the
AlCaReco selections instead of the analysis ones, and the analysis
path SHALL be byte-for-byte preserved when the flag is unset. The
script's cut-flow SHALL classify each surviving B+ candidate as "true"
(all three leaf tracks gen-matched to the generated B+ decay products)
or "fake" (otherwise), and the gen-match flag SHALL be an axis on the
output histograms so signal efficiency and per-event fake rate can be
computed under each preset.

#### Scenario: Analysis path is preserved when the flag is unset
- **WHEN** `btojpsik.py` is invoked without `--alcarecoPreset`
- **THEN** the output histograms are byte-for-byte identical to the pre-change behaviour for the same input dataset

#### Scenario: AlCaReco preset routes through the new selections
- **WHEN** `btojpsik.py --alcarecoPreset B` is invoked
- **THEN** the cut-flow applies the cuts returned by `get_bkmm_alcareco_selections('B')`, the existing analysis cuts are bypassed, and the output histograms carry the gen-match flag

#### Scenario: Phase-1 gate evaluation produces a decision-ready summary
- **WHEN** the histograms from preset-A and preset-B invocations are post-processed
- **THEN** the signal efficiency (truth-matched B/A), per-event fake-rate reduction (A/B), and mass-shape preservation (width and mean within 2σ-stat) are computed per preset, with PASS/FAIL labels against the gates (≥ 70%, ≥ 5×, preserved)

---

### Requirement: Phase-2 ALCARECO Confirmation Tooling

The repository SHALL provide a Phase-2 diagnostic tool (`jpsix_preset_compare.py`, sibling of `jpsix_diagnostics.py`) that reads two ALCARECO files (one per preset, from cmsRun on Run2016H Charmonium RAW). The tool SHALL operate in FWLite and SHALL produce: overlaid mother-mass distributions per channel under A and B; per-channel summary tables of cands/event, events with ≥1 cand, and tight-window fraction; a wall-time-per-event and bytes-per-event table from FastTimerService and `edmEventSize` output; and a track-deduplication compression metric averaged over events. The tool SHALL print PASS/FAIL against the three Phase-2 gates (cands/event ≤ 5; tight-window fraction ≥ 30%; tight-window yield ≥ 100 / 1000 events) for each of the four non-V0 channels.

#### Scenario: Side-by-side mass overlay for both presets
- **WHEN** the tool is invoked with two ALCARECO ROOT files corresponding to presets A and B of the same input data
- **THEN** it emits one mother-mass plot per channel with two overlaid histograms (A, B) and one summary table per channel listing the operating-point metrics

#### Scenario: Size and timing tables reported per preset
- **WHEN** the tool consumes the FastTimerService and `edmEventSize` outputs alongside the two ALCARECO files
- **THEN** it produces a per-preset table of wall-time per event and ALCARECO bytes per event with per-branch breakdown

#### Scenario: Dedup-compression metric reported per event
- **WHEN** the tool processes an ALCARECO file
- **THEN** it produces a per-event histogram of (number of cloned tracks in merged collection) and (sum-over-channels of unique-leaf-tracks if each channel had its own AlignmentTrackSelector), and reports the average ratio (compression factor) over all processed events

#### Scenario: V0-mode channels are preset-invariant in the comparison
- **WHEN** the tool processes preset-A and preset-B files
- **THEN** for B0→Ks, Λb, ψ(2S) the overlaid histograms collapse to a single curve (numerically identical across presets within statistical floor)

---

### Requirement: Continuous Worklog Slide Deck

A worklog slide deck `slides/jpsix-selection-comparison-worklog.tex` (beamer, same theme as the existing `jpsix-alcareco-producer.tex`) SHALL be created at the start of work on this change and grown frame-by-frame as iterations, runs, and decisions happen. The deck SHALL accumulate the full iteration history; earlier frames SHALL NOT be deleted or rewritten when later iterations supersede them, so the deck constitutes a chronological lab notebook of the change. Every iteration that produces a histogram or a number SHALL be accompanied by at least one new worklog frame; every "tried X, didn't work, switched to Y" SHALL produce a frame with before/after numbers; every decision or conclusion SHALL produce a frame stating it explicitly with citations to the supporting earlier frames. The main public-facing deck `jpsix-alcareco-producer.tex` SHALL receive a new headline-summary frame only after Phase-2 numbers are in, not during iteration.

#### Scenario: Worklog frame accompanies each iteration
- **WHEN** a Phase-1 cut iteration is run on the BuToJpsiK MC and produces signal-efficiency / fake-rate / mass-shape numbers
- **THEN** at least one new frame is added to `jpsix-selection-comparison-worklog.tex` recording the cut values applied, the numbers produced, and a one-sentence observation, in the same git commit (or an immediately-adjacent one) as the histogram run

#### Scenario: Decision frame cites earlier worklog frames
- **WHEN** the per-channel preset decision is made
- **THEN** the worklog deck contains a "Decision" frame stating the per-channel preset assignment and citing the earlier iteration / Phase-2 frames whose numbers justify the choice (by frame number, hyperlink, or quoted summary)

#### Scenario: Worklog deck is complete at change close
- **WHEN** this change is archived (proposal moves to `openspec/changes/archive/`)
- **THEN** the worklog deck compiles cleanly to PDF, frames are in chronological order, and every Phase-1 iteration / Phase-2 sub-phase / decision / anomaly produced during the work has at least one corresponding frame

---

### Requirement: Comparison Results Artifact

This change directory SHALL contain a `results.md` file populated after the comparison campaign completes. The file SHALL be split into a Phase-1 section (MC, B+ only) and a Phase-2 section (data, all 7 channels). The Phase-1 section SHALL include the A and B cut tables, the signal-efficiency / fake-rate-reduction / mass-shape metrics under each preset, and the locked B+ preset values. The Phase-2 section SHALL include three tables — data mass quality (Phase 2a), output size + compute time (Phase 2b), dedup compression (Phase 2c) — each indexed by channel and preset where applicable, with one paragraph per sub-phase describing the interpretation. The per-channel preset assignment chosen for production SHALL be recorded in `results.md` with a one-sentence rationale anchored to the Phase-2 gates. Channels that fail the Phase-2 gates under both A and B SHALL be flagged as "insufficient at Stage 1 alone; relies on Stage-2 cleanup", with no further AlCaReco-stage cuts attempted.

#### Scenario: Results document per-channel preset choices
- **WHEN** the comparison campaign concludes and `results.md` is finalised
- **THEN** `results.md` names a preset (A or B) for each of the four non-V0 channels (B+, B0→K*0, Bs→φ, Bc), citing the gate evaluations from Phase 1 (for B+) and Phase 2a (for all four) as justification, and the corresponding per-channel defaults in `ALCARECOTkAlJpsiX_cff.py` match those choices

#### Scenario: V0-mode channels are not preset-switched in results
- **WHEN** `results.md` is finalised
- **THEN** the V0-mode channels (B0→Ks, Λb, ψ(2S)) are documented as fixed at "V0 quality + minimal kinematic" with no preset choice, and Phase 2a's V0-mode rows show identical numbers across A and B as a structural check

#### Scenario: Cut translation B+ → other non-V0 channels documented
- **WHEN** Phase 1 settles the B+ cuts and the same cuts are translated to B0→K*0, Bs→φ, Bc
- **THEN** `results.md` contains a one-row-per-cut translation table per channel, each row including the source B+ value, the channel-translated value, and a one-sentence physics justification for the translation
