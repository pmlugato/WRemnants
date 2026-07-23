## Context

Four AlCaReco streams will ship in production across two PDs
(Charmonium: J/psi+X, TkAlJpsiMuMu; SingleMuon: D*, KShort, Lambda).
The validation described here runs on the same configuration
production will run. There is no diagnostic flag, no debug-only
branch, no runtime switch. Whatever the validation reads is
production output.

## Goals / Non-Goals

Goals:

- One CMSSW build that runs all four streams.
- One production configuration per stream — no runtime toggles.
- Persist alignment across the four streams so the same categories
  of output (tracks + dE/dx + candidate VCC + selector outputs +
  associations) are available for downstream analysis on every
  stream where the physics allows it.
- Cross-dataset Condor tarball reuse: build once, submit twice.
- Bit-invariance of J/psi+X on events shared with the last smoke
  (i.e. both the tight-muon gate and the new candidate gate
  accept the event) — proves PV2 only widens acceptance.

Non-goals:

- Retuning any candidate cut.
- Stage-2 CVH refit for D* / V0.
- A diagnostic runtime mode.

## Decisions

### D1: Preset removal — inline literals, no compatibility shim

Delete the preset switch outright. Inline preset B numeric values
directly into each per-channel producer clone. Do not keep a
"preset name" comment above each block linking back to the
selection-comparison change proposal — the numeric values speak
for themselves.

Alternative — "leave the switch, always default to B":
Rejected. Env-var dependencies leak in silently and are impossible
to reason about from a code read alone.

### D2: Candidate-level event gate — merge then count-filter

Demote `ALCARECOTkAlJpsiXGoodMuons` to `filter = False`. Merge the
8 per-channel VCC outputs into `ALCARECOTkAlJpsiXAnyCandidate`
using a stock CMSSW merger. Terminate the sequence with
`ALCARECOTkAlJpsiXCandFilter` (`CandViewCountFilter`,
`minNumber = 1`) between the merger and
`CompositeDaughterTrackProducer`.

Alternative — "split into per-channel paths OR'd on the output
module's `SelectEvents`":
Rejected. Would require the RECO chain to be sequenced in every
path; the merge + count filter is a single sequence, no
duplication.

Alternative — "custom EDFilter reading N InputTags":
Deferred. Only build one if no stock merger accepts a
`VInputTag` for VCC collections.

Rate impact: a 1-file smoke measures the accept-rate ratio to
the last tight-gate smoke before adoption (task 0.6).

### D3: Persist alignment across streams — production-permanent, no flag

The four streams keep different categories today (see proposal
table). Adopt one policy across all four:

- Cloned tracks with `originalIndex` map.
- dE/dx (strip Harmonic2, pixel Harmonic2, joint).
- Candidate VCC — one per emitted channel.
- Muon selector output collection (streams with muons).
- Track -> muon association (streams with muons).
- Standard L1 / HLT / DCS / offlinePrimaryVertices.

Adjustments to the current cff files:

- **TkAlJpsiMuMu**: add three `DeDxValueMapProjector` clones and
  their `keep` lines; persist the `ALCARECOTkAlJpsiMuMuGoodMuons`
  collection (the tight muon selector output); add a
  `AlignmentTrackToMuonAssociator` producer.
- **D***: add a candidate emitter that writes a
  `VertexCompositeCandidateCollection` for D* (with the D0 as the
  intermediate). The selection logic is the same three-body decay
  criteria the current stream already applies internally to filter
  tracks; we surface the surviving candidates as first-class
  objects.
- **V0 streams**: audit after sparse-add. Align keeps to the
  policy where applicable.

### D4: D* candidate emitter — new small producer, not a keep-line

`AlignmentTrackSelectorWithIndexMap` running
`AlignmentThreeBodyDecayTrackSelector` filters tracks; it does not
emit candidates. To persist a D* candidate we either:

- (a) Add a `RecoChargedCandidate` / VCC emitter external to the
  selector, running the same three-body-decay logic and writing
  surviving triplets as candidates.
- (b) Modify `AlignmentThreeBodyDecayTrackSelector` to emit a VCC
  side-collection alongside its filter decision.

Choose (a). The selector code came from upstream; keeping it
untouched preserves its role as a track filter and localises the
new logic in a new plugin
(`ThreeBodyDecayCandidateProducer`, or reuse
`TwoBodyDecayCandidateProducer`'s existing three-body branch if
one exists) or in a python-only clone if a stock producer already
does the job.

Task 2.3 confirms which option is minimal.

### D5: Two Condor tarballs vs one

One tarball, two `.sub` files, one shared `run.sh` dispatching on
recoskim filename. Same reasoning as the earlier J/psi+X-only
batch: binaries are dataset-agnostic; only the recoskim `.py` and
input LFN differ.

### D6: Job packing

One cmsRun per input file per PD, with multiple AlCaReco paths in
`process.schedule` (Charmonium: J/psi+X + JpsiMuMu; SingleMuon:
D* + KShort + Lambda). RAW2DIGI + L1Reco + RECO are shared across
the streams in the same job, avoiding duplicate work.

RSS ceiling risk: bump Condor `RequestMemory` from 6 GB to 8 GB
for the multichannel jobs. If any smoke run overshoots, split
into two cmsRuns.

### D7: Plotting

`uproot` + `hist` + `matplotlib`, single driver at
`scripts/btojpsik/plots_preprod_validation.py`. No ROOT
dependency. Outputs `.png` + `.pdf` + `summary.json`.

Cut-efficiency numbers come from cmsRun stderr `LogInfo` counter
lines — the C++ producers already emit per-cut yield reductions
on every job. The plotting script parses the log tail; no
persisted collection is needed for this deliverable.

V0 loosening yield: measured once by a separate offline job that
re-runs the CMSSW-default V0 producer on the same input files
and records the ratio to the loosened yield. Not part of the
AlCaReco.

### D8: Input files

- Charmonium: reuse the last smoke's seed file for the pre-PV2
  rate estimate and the bit-invariance check; then 100 files from
  the same 2016H run range.
- SingleMuon: 1 file from Run2016H SingleMuon for smoke, then
  100 files covering the same run range as the Charmonium 100.

Cached via `dasgoclient` under `condor/multichannel_alcareco/lfns/`.

## Risks / Trade-offs

- **PV2 rate blow-up**: if the candidate gate accepts >= 2x the
  previous rate on the smoke seed file, the launch stops for user
  review before it lands (task 0.6).
- **RSS with two AlCaReco paths per cmsRun**: bump memory to 8 GB;
  fall back to two cmsRuns if smoke overshoots.
- **D\* candidate emitter is new code**: needs its own bit-
  invariance vs the pre-emitter cloned-track content (adding the
  emitter must not change which events / tracks the stream keeps,
  only add a candidate collection).
- **V0 sparse-checkout unknowns**: KShort / Lambda cff files may
  need extra plumbing (originalIndex map, dE/dx projection) to
  align with the persist policy. Task 1.1 confirms scope.

## Migration Plan

Not a data migration. Stage-2 CVH refit continues to consume the
existing J/psi+X output collections; new persisted collections in
the other streams are additive.

Rollback: revert the cff changes and delete
`condor/multichannel_alcareco/`. J/psi+X preset B production is
unaffected by cff-level rollback because PV1 is a no-op for
callers (env var already defaulted to B) and PV2 is behind a
smoke gate.

## Open Questions

1. Do `ALCARECOTkAlKShortTracks_cff.py` and
   `ALCARECOTkAlLambdaTracks_cff.py` need any plumbing beyond a
   sparse-checkout add (originalIndex map for dE/dx projection,
   etc.)? Task 1.1 verifies.
2. For the D* candidate emitter, is there already a stock
   `ThreeBodyDecayCandidateProducer` in CMSSW, or do we need to
   write one? Task 2.3 confirms.
3. Should the tight-muon collection be persisted only under
   TkAlJpsiMuMu, or also under J/psi+X (in addition to the loose
   one)? Current plan: only TkAlJpsiMuMu keeps tight; J/psi+X
   keeps loose. The muon comparison plot reads both streams,
   which run on the same events. If the alignment / calibration
   downstream ever wants tight muons keyed against J/psi+X's
   dedup tracks, that is a follow-up.
