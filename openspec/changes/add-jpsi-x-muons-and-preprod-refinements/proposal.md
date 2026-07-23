# Change: Save J/psi muons, add J/psi-only channel, and pre-production refinements

## Why

Preset B is validated end-to-end and the 100k Condor headline sample landed
under `finalize-jpsi-x-preset-b-production`. Before we commit resources to a
full-2016H (and later multi-era) production, six small items should be
locked in — none of them are large but they all affect what the persisted
AlCaReco stream contains, so they must land in **one** production build:

1. **Store the J/psi muons (not just the muon tracks).**
   Currently only tracks survive into the output stream
   (`ALCARECOTkAlJpsiX_Output_cff.py` keeps `*_ALCARECOTkAlJpsiX_*` for the
   deduplicated `AlignmentTrackSelectorWithIndexMap` track collection and
   the per-channel VCC collections whose daughters hold only `TrackRef`s).
   The `reco::Muon` objects that were used by `TkAlGoodIdMuonSelector` to
   filter the J/psi dimuon candidates are dropped after Stage 1. We need
   the muon-level information (`isGlobalMuon`, `isTrackerMuon`,
   `numberOfMatches`, `globalTrack.hitPattern.numberOfValidMuonHits`,
   station-level segment info, isolation R03) to study impacts of muon
   quality on the alignment / calibration fit and to enable the new
   J/psi-only channel below. To avoid duplicating storage, the J/psi
   candidate collections must **not** value-copy the muon information —
   they must retain a lightweight pointer (index/ref/ValueMap) into a
   single persisted muon collection kept side-by-side with the track
   collection.

2. **Add a J/psi-only channel with looser muon quality, run in
   production alongside the seven existing channels.**
   All seven current channels require a *bachelor / X*. We add an
   eighth channel — dimuon-only, no X — run over a **looser** muon
   selector (`isGlobalMuon | isTrackerMuon` instead of
   `isGlobalMuon & isTrackerMuon`, with the `globalTrack.*` sub-cuts
   dropped since a tracker-only muon has no `globalTrack`). Concretely:
   `(isGlobalMuon | isTrackerMuon) & abs(eta) < 2.5 &
   numberOfMatches > 1`. This channel ships **in production**, not just
   as a smoke test: it is the source of the J/psi sample we intend to
   feed to alignment / calibration Stage-2 workflows that want higher
   statistics than the bachelor channels can supply and want the muon-
   quality axis available at fit time. It also happens to give us a
   same-event tight-vs-loose J/psi comparison, so a follow-up proposal
   can decide whether to globalise the loose selector across all
   channels; but the channel exists on its own merits regardless of
   that follow-up outcome. Until that follow-up lands, the seven
   existing channels stay on the tight selector so their headline
   sample remains bit-invariant.

3. **Use daughter-track PDG-id instead of the collection name for the
   downstream mass hypothesis.**
   `JpsiXCandidateProducer` and `TwoBodyDecayCandidateProducer` already
   set the correct `pdgId` on every emitted `RecoChargedCandidate`
   daughter (`JpsiXCandidateProducer.cc:755-757`,
   `TwoBodyDecayCandidateProducer.cc:239-253`). `V0Producer` does the
   same for Ks / Lambda daughters. Downstream — specifically
   `JpsiKCandidateSplitter` — still disambiguates species **by the input
   collection name** (`JpsiKCandidateSplitter.cc:22-25`), which is
   fragile (name changes break Stage 2, and any per-track species logic
   at Stage 2 has no channel-local info once the tracks are packed into
   a flat collection). The daughter pdgId is already present; downstream
   must trust it. This also requires that every producer we ship writes
   a **non-zero** pdgId on every daughter it emits, so the invariant is
   safe to rely on.

4. **Raise the bachelor / K*0-daughter / phi-daughter η cut from 1.8 to
   2.5 under preset B (and C).**
   Preset A currently uses 2.4; we take everything to the tracker
   acceptance edge (|η| = 2.5). Capping the bachelor / K* pion / phi
   kaon at 1.8 throws away large-η bachelor-track alignment leverage.
   This affects `_NON_V0_PRESETS['BPlus'|'Bc']['maxBachelorEta']`
   (currently 1.8) and `_NON_V0_PRESETS['Kstar'|'Phi']['maxDaughterEta']`
   (currently 1.8) under preset B, and correspondingly under preset C
   which inherits B. Preset A also moves from 2.4 → 2.5 to keep the
   three presets consistent on the η boundary — the bachelor track can
   travel up to the tracker edge in every preset.

5. **Loosen the V0 selection on the local `ALCARECOTkAlV0Candidates`
   clone once, for all V0 channels (B0→Ks, Λb, ψ(2S)).**
   `V0Producer` defaults are tuned for narrow, high-purity Ks / Lambda
   for reconstruction, not for downstream Kalman-refit-only alignment
   candidates. Our local clone (`ALCARECOTkAlV0Candidates_cff.py`) only
   touches `tkPtCut = 0.1`; everything else — `vtxChi2Cut = 6.63`
   (CL ≈ 1 %), `vtxDecaySigXYCut = 15`, `cosThetaXYCut = 0.998`,
   `tkIPSigXYCut = 2`, `tkChi2Cut = 10` — sits at the default tight
   point. All three V0-mode channels (B0→Ks, Λb, ψ(2S)) share this one
   clone, and all three are yield-limited under the current tight
   defaults. We commit to a **single, defensible looser working point**
   applied unconditionally on the clone, not a study: `vtxChi2Cut =
   15`, `vtxDecaySigXYCut = 10`, `cosThetaXYCut = 0.995`,
   `tkIPSigXYCut = 1`, `tkChi2Cut = 10` (unchanged). These are moderate
   step-outs — each still tighter than "no cut" — and pick up the
   long-lived Λ / Ks tails whose displacement significance is
   comfortably above default without needing to see the mass fit.
   No optimisation study, no env-var switch; a proper working-point
   scan is deferred to a follow-up if the loosened default proves
   inadequate.

6. **Fix the ψ(2S) → J/ψ π+π- wiring: two-track dipion, not V0 Ks.**
   ψ(2S) → J/ψ π+π- (BR ≈ 34.7 %) is the correct hadronic decay of
   ψ(2S) to J/ψ + two pions. The two pions share the ψ(2S) vertex;
   they do not form a Ks. The current
   `ALCARECOTkAlJpsiXPsi2SCandidates` PSet has
   `xSrc = cms.InputTag('ALCARECOTkAlV0Candidates', 'Kshort')`, i.e.
   requires the π+π- to pass V0Producer's Ks-quality displaced-vertex
   cuts (Lxy significance ≳ 15, cosθ_XY ≳ 0.998, Ks-mass window). That
   throws away essentially all real signal — the 200-event smoke found
   only 7 ψ(2S) candidates and those are almost entirely combinatoric
   fakes.
   Fix: add a new `TwoBodyDecayCandidateProducer` instance
   `ALCARECOTkAlJpsiXPiPiCandidates` — analogous to the K*0 / φ
   producers — building a π+π- pair from `generalTracks` with only a
   physical dipion mass window (0.28, 0.65 GeV) and no vertex fit; then
   repoint the ψ(2S) producer at it. The V0-Ks path is removed
   entirely; there is no fallback. ψ(2S) itself may be prompt (~80 %
   of the ψ(2S) rate, from direct production at the PV) or non-prompt
   (~20 %, from B decays); the new wiring accepts both, matching preset
   B's philosophy of raw kinematics with displacement filtering deferred
   downstream. Side effect: ψ(2S) leaves the V0-mass regime and joins
   the PDG regime for downstream species dispatch (pions carry
   `pdgId = ±211`).

All six items are code-only and land as one CMSSW build.

## What Changes

- **New:** `ALCARECOTkAlJpsiXMuons` — a persisted `reco::MuonCollection`
  in the output stream, carrying every muon that was selected by the
  Stage-0 muon selectors (both the existing tight selector and the new
  loose one — see below), so J/psi daughter tracks can be resolved to
  their `reco::Muon` without duplication.
- **New:** an `edm::Association<reco::MuonCollection>` (or equivalent
  ValueMap) named `ALCARECOTkAlJpsiXTrackToMuon` mapping each track in
  the persisted `ALCARECOTkAlJpsiX` track collection to the
  corresponding `reco::Muon` in `ALCARECOTkAlJpsiXMuons` (null for
  non-muon tracks: bachelors, V0 daughters, K*0/phi daughters).
- **New:** a second, looser `MuonSelector` — `ALCARECOTkAlJpsiXLooseMuons`
  — with cut string `(isGlobalMuon | isTrackerMuon) & abs(eta) < 2.5 &
  numberOfMatches > 1`. The `globalTrack.*` sub-cuts are dropped
  because tracker-only muons carry no `globalTrack`. Consumed by the
  new production J/psi-only channel (see item 2); a follow-up may
  promote it to be the global selector once the impact on the seven
  bachelor channels is measured.
- **New:** a J/psi-only per-channel path,
  `ALCARECOTkAlJpsiXJpsiOnlyCandidates`, produced by a
  `TwoBodyDecayCandidateProducer` on `ALCARECOTkAlJpsiXLooseMuons` (no
  bachelor / X on top). Emitted as its own persisted collection
  `ALCARECOTkAlJpsiXJpsiOnlyResonances` on par with the seven existing
  channels. Trigger paths and sequence membership updated accordingly.
- **Modified:** `JpsiXCandidateProducer.cc` and
  `TwoBodyDecayCandidateProducer.cc` — no daughter is emitted with
  `pdgId == 0`; the daughter species is fully carried by
  `daughter->pdgId()` on the persisted `RecoChargedCandidate`
  (already true today; this proposal ratifies it as a spec invariant).
- **Modified:** `JpsiKCandidateSplitter.cc` (Stage-2 CVH) — species /
  mass hypothesis lookup uses the daughter `pdgId()` instead of the
  input collection name. Collection layout continues to define the
  daughter *order*; only the species tag changes source.
- **Modified:** presets A, B, and C — `maxBachelorEta` on B+ / Bc and
  `maxDaughterEta` on K*0 / phi all raised to **2.5** across every
  preset (from 1.8 under B/C and from 2.4 under A). V0-mode channels
  unaffected (V0 daughters go through the V0Producer, not through the
  preset's per-channel η cap).
- **Modified:** `ALCARECOTkAlV0Candidates_cff.py` — sets `vtxChi2Cut =
  15`, `vtxDecaySigXYCut = 10`, `cosThetaXYCut = 0.995`,
  `tkIPSigXYCut = 1` on the local clone (overriding
  `generalV0Candidates` defaults). Applied unconditionally to the two
  remaining V0-mode channels (B0→Ks, Λb) via the shared clone. ψ(2S)
  moves off the V0 path (see below) so is no longer affected.
- **New:** `ALCARECOTkAlJpsiXPiPiCandidates` — a
  `TwoBodyDecayCandidateProducer` instance building a two-track π+π-
  pair from `generalTracks` (dipion mass window (0.28, 0.65),
  `daughterPdgId = 211`, no vertex fit, `maxTrackTrackDOCA = 0.03` cm).
  Consumed by the ψ(2S) producer instead of the V0 Ks output.
- **Modified:** `ALCARECOTkAlJpsiXPsi2SCandidates.xSrc` — repointed
  from `ALCARECOTkAlV0Candidates:Kshort` to
  `ALCARECOTkAlJpsiXPiPiCandidates`.
- **Modified:** `TwoBodyDecayCandidateProducer.cc` gains an optional
  `maxTrackTrackDOCA` PSet field (default = infinity → no cut, so
  K*0 / φ / J/psi / J/psi-only remain bit-invariant). Uses the same
  static straight-line 3D DCA helper `JpsiXCandidateProducer` already
  applies for its `maxBachelorMuTrackDOCA` cut — no Kalman fit, no
  `TransientTrackBuilder` dependency.
- **Modified:** `ALCARECOTkAlJpsiX_Output_cff.py` — adds the new muon
  collection, the track-to-muon association, and the J/psi-only
  candidate collection to `outputCommands`.

## Impact

- **Affected specs:** `alcareco-jpsi-x`.
- **Affected code:**
  - `Alignment/CommonAlignmentProducer/python/ALCARECOTkAlJpsiX_cff.py`
    (Stage-0 selector clone, new J/psi-only path, preset eta knobs,
    muon association producer wiring)
  - `Alignment/CommonAlignmentProducer/python/ALCARECOTkAlJpsiX_Output_cff.py`
    (new keep statements)
  - `Alignment/CommonAlignmentProducer/python/ALCARECOTkAlV0Candidates_cff.py`
    (expose V0 quality knobs; landing values from the study)
  - `Alignment/CommonAlignmentProducer/python/TkAlMuonSelectors_cfi.py`
    (add loose selector)
  - `Alignment/CommonAlignmentProducer/plugins/JpsiXCandidateProducer.cc`
    and `TwoBodyDecayCandidateProducer.cc` (assert daughter pdgId is
    always non-zero; no functional change expected)
  - Optionally a small new plugin
    `AlignmentTrackToMuonAssociator` under
    `Alignment/CommonAlignmentProducer/plugins/` to emit the
    `edm::Association<reco::MuonCollection>` from a track collection
    (deduplicated `ALCARECOTkAlJpsiX`) + a muon collection.
  - Stage-2 CVH: `Analysis/HitAnalyzer/plugins/JpsiKCandidateSplitter.cc`
    switches species lookup to `daughter->pdgId()`.
- **Not affected:** the Stage-1 → Stage-2 wiring, the preset switch
  mechanism itself, the V0 output collection layout (Ks / Lambda names
  unchanged), Stage-1 producers' physics (except (3) & (4)), the
  existing seven channels.
- **Production impact:** one CMSSW rebuild before the full-data
  production launch; all changes are single-preset (preset B), no new
  preset needed.
