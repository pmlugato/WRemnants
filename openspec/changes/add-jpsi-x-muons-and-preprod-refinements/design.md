# Design: Save J/psi muons, add J/psi-only channel, pre-production refinements

## Context

This proposal collapses five items into one change so the resulting stream
layout (one output list, one CMSSW build) ships as a single production
version:

| # | Item | Kind |
|---|------|------|
| 1 | Persist J/psi muons + track→muon lookup | New capability |
| 2 | J/psi-only channel with loose muon selector — production channel, also a study vehicle | New capability |
| 3 | Use `daughter->pdgId()` for species downstream | Spec ratification + Stage-2 code change |
| 4 | Raise bachelor / X-daughter η → 2.5 (all presets) | Preset value change |
| 5 | Loosen V0 quality on the shared AlCaReco clone (one-shot, all V0 channels) | Preset value change |
| 6 | ψ(2S) → J/ψ π+π- via a two-track dipion (leave the V0 Ks path) | Physics bug fix |

Items 3, 4, and 5 are essentially trivial once the values are agreed.
Items 1 and 2 warrant longer design notes.

## Decision 1 — Muon-storage layout

**Chosen:** save a per-event `reco::MuonCollection`
(`ALCARECOTkAlJpsiXMuons`) *plus* a track→muon
`edm::Association<reco::MuonCollection>`
(`ALCARECOTkAlJpsiXTrackToMuon`), keyed on the deduplicated persisted
track collection `ALCARECOTkAlJpsiX`.

**Alternatives considered:**

- *Extend `RecoChargedCandidate` with `MuonRef`.* Would require a
  custom candidate type at Stage 1, and every downstream module that
  handles our candidates would need to know about it. Rejected: too
  invasive.
- *Duplicate the muon-level info onto each J/psi daughter
  `RecoChargedCandidate` (e.g. as user-int / user-float on a custom
  candidate).* Rejected explicitly by the user ("as to not duplicate
  its storage").
- *Persist the muons but rely on candidate → daughter → track →
  `edm::View<reco::Muon>::refByInnerTrack` at Stage 2.* Works if the
  original `reco::Muon` handle is intact, but the ALCARECOTkAlJpsiX
  track collection is a **cloned** collection produced by
  `AlignmentTrackSelectorWithIndexMap` — `muon.innerTrack()` points
  into the original `generalTracks`, not into our clone. We would
  still need an explicit association from the clone to the muon;
  might as well store it.

**Why an association rather than storing the muon collection
side-by-side and letting Stage 2 look up by inner-track pointer:**
the deduplicated track collection is a clone whose provenance is not
`generalTracks`, so the inner-track ref stored on `reco::Muon`
does not directly match the tracks that survive into the persisted
stream. A one-shot associator built at Stage 1 (before dedup / clone,
against `generalTracks`) plus the index map already emitted by
`AlignmentTrackSelectorWithIndexMap` gives Stage 2 a direct
`Muon = association->at(TrackRef into ALCARECOTkAlJpsiX)` lookup at
zero cost.

**Which muons go in?** Union of the tight and loose selector outputs.
- Tight (`ALCARECOTkAlJpsiXGoodMuons`, existing) — those seed the six
  bachelor-plus-J/psi channels + the K*0/φ dimuon inputs.
- Loose (`ALCARECOTkAlJpsiXLooseMuons`, new) — a superset, seeding the
  J/psi-only channel.
Because the loose selector is a superset of the tight one **at the
value level** (relaxed cuts), a single `keep` on the loose selector's
output would in principle suffice. To avoid coupling the persisted
collection name to selector semantics that may change later, we
instead persist a dedicated `ALCARECOTkAlJpsiXMuons` (the union,
built by a small merger EDProducer that concatenates unique
`reco::Muon` pointers from both selectors' outputs). Cost: negligible
compared to the tracks already in the stream (a few muons per event).

**Why muons are needed at all in this AlCaReco stream:** the J/psi-only
channel with a **looser** muon selector is only useful if Stage 2 can
recompute muon-quality-conditioned fits (`isGlobalMuon`,
`isTrackerMuon`, `nMatches`, `chi2`, etc.) from the persisted event.
Without the `reco::Muon` object those attributes are unreachable.

## Decision 2 — J/psi-only channel wiring, and the loose-muon config

**Framing:** the J/psi-only channel is a **production channel**, not
just a study vehicle. It ships as an eighth persisted per-channel
collection on every production event and is intended for downstream
Stage-2 workflows that need a higher-statistics dimuon sample and want
the muon-quality axis available. As a side benefit, having it on the
loose selector while the seven bachelor channels still run on the
tight selector gives us a same-event tight-vs-loose J/psi comparison
that a follow-up proposal can use to decide whether to globalise the
loose selector. Whether or not that follow-up ever lands, the
J/psi-only production channel stays. Until it lands, the seven
existing channels stay on the tight selector so their headline sample
remains bit-invariant.

**Chosen:** run a second, dedicated `TwoBodyDecayCandidateProducer`
(`ALCARECOTkAlJpsiXJpsiOnlyCandidates`) whose `muonSrc` points at the
new loose selector. It emits a `VertexCompositeCandidateCollection`
of J/psi (μμ) candidates, no bachelor / X requirement.
`ALCARECOTkAlJpsiXJpsiOnlyResonances` = that producer's output made
persistable (its daughter track refs will be re-keyed onto the
deduplicated track collection by the existing
`VertexCompositeCandidateRemapper` sequence).

**Loose-muon cut string.** `TkAlGoodIdMuonSelector` today reads
`isGlobalMuon & isTrackerMuon & numberOfMatches > 1 &
globalTrack.hitPattern.numberOfValidMuonHits > 0 & abs(eta) < 2.5 &
globalTrack.normalizedChi2 < 20`. Two of those sub-cuts —
`globalTrack.hitPattern...` and `globalTrack.normalizedChi2` — dereference
`globalTrack` and are only meaningful for **global** muons. Under the
proposed loose selector we accept tracker-only muons as well, and a
tracker-only muon has `globalTrack.isNull()` — the ExpressionEvaluator
either rejects (safe) or dereferences null. We handle this at the cut
level by dropping the two `globalTrack.*` sub-cuts entirely. Final
cut string:
```
(isGlobalMuon | isTrackerMuon) & abs(eta) < 2.5 & numberOfMatches > 1
```
`numberOfMatches` is a `reco::Muon` attribute defined regardless of
whether the muon is tracker- or global-flavoured, so this cut is
well-formed for the union population. The `numberOfMatches > 1`
threshold matches the tight selector so we're not doubly-loosening
(loose = OR of flavours only). If study results argue for it, a
follow-up can drop `numberOfMatches > 1` as well.

**Why not "just save all μμ within [2.95, 3.25]" as a separate
side product of `ALCARECOTkAlJpsiXJpsiCandidates`:** the tight-selector
J/psi input is what feeds the seven downstream channels. Keeping it
untouched and adding an *independent* J/psi-only producer means the
existing production is bit-invariant on all seven bachelor channels;
we cannot regress the headline sample.

**Sequence membership:** the J/psi-only producer is appended to
`seqALCARECOTkAlJpsiX`. `pathALCARECOTkAlJpsiX` remains the single
production path so operators do not need to add new event-selection
paths in the AlCaReco stream config.

**Trigger / DCS gating:** unchanged. Same Stage-0 gate for both the
existing seven channels and the new J/psi-only channel.

**Mass window:** same as the existing J/psi window (`minMass = 2.95`,
`maxMass = 3.25`) — the *muon quality* is what's loosened, not the
resonance window.

**Follow-up study (in `tasks.md`, feeds a separate future proposal):**
on the first full-data batch, compare tight-selector J/psi and
loose-selector J/psi mass shapes on the *same events*. Metrics: peak
μ shift, σ change, S/√B in the ±3σ window, and the fraction of J/psi's
in the loose sample that are NOT in the tight sample. Decision rule
for the follow-up proposal that would drop the tight selector: σ
change ≤ +20 %, μ shift ≤ 5 MeV, peak S/√B ≥ 90 % of the tight-only
sample's. This study does not gate the current proposal — the loose
selector and the J/psi-only channel ship in production regardless.

## Decision 3 — PDG-id-first species dispatch

The daughter pdgId is already set at Stage-1 emission:

- `TwoBodyDecayCandidateProducer.cc:239-253` sets
  `RecoChargedCandidate(charge, lv, vtx, pdgPos|pdgNeg)` where
  `pdgPos = ±firstDaughterPdgId_`, `pdgNeg = ±secondDaughterPdgId_`
  based on the swap flag. Both are non-zero for every configured
  producer instance.
- `JpsiXCandidateProducer.cc:755-757` (track mode) sets the bachelor
  pdgId from the `bachelorPdgId` PSet field (non-zero for B+ / Bc).
- `V0Producer` sets pdgId = ±211 (pion), ±2212 (proton) on the
  Ks / Λ daughter `RecoChargedCandidate`s.

The spec change is to **assert** — as a producer-level invariant — that
no daughter ever ships with `pdgId == 0`, and to have downstream
(`JpsiKCandidateSplitter.cc`) read `daughter->pdgId()` in place of the
current input-collection-name check.

**Backwards-compat:** the collection-name check in
`JpsiKCandidateSplitter.cc:22-25` is a comment + no runtime check; the
splitter currently discards the daughter species and re-derives it in
downstream cfi. After this change the splitter propagates the pdgId
through its output (e.g. by decoration on the emitted track collection
or by an int branch), so the Stage-2 kaon-vs-pion mass hypothesis is
set from a single source of truth.

**Explicit invariant:** the seven existing per-channel
`JpsiXCandidateProducer` instances have well-defined
`bachelor / firstDaughter / secondDaughterPdgId` values in the cff:

| Channel | daughter 0 (from J/psi) | daughter 1 | daughter species |
|---------|-------------------------|------------|------------------|
| B+ | J/psi(μμ) | K± (bachelor) | ±321 |
| Bc | J/psi(μμ) | π± (bachelor) | ±211 |
| B0→K*0 | J/psi(μμ) | K*0(Kπ) | ±321, ∓211 |
| Bs→φ | J/psi(μμ) | φ(KK) | ±321, ∓321 |
| B0→Ks | J/psi(μμ) | Ks(ππ) | ±211, ∓211 |
| Λb | J/psi(μμ) | Λ(pπ) | +2212 (or +211), ∓211 (or ∓2212) — Λ vs Λ̄ |
| ψ(2S)→J/ψ Ks | J/psi(μμ) | Ks(ππ) | ±211, ∓211 |
| **new: J/ψ-only** | μ+ | μ− | ±13 |

## Decision 4 — Eta cut → 2.5 (all presets)

Mechanical. All three presets get the same tracker-edge value:

- `_NON_V0_PRESETS['BPlus']['maxBachelorEta']`: 1.8 (B/C) or 2.4 (A) → 2.5
- `_NON_V0_PRESETS['Bc']['maxBachelorEta']`: 1.8 (B/C) or 2.4 (A) → 2.5
- `_NON_V0_PRESETS['Kstar']['maxDaughterEta']`: 1.8 (B/C) or 2.4 (A) → 2.5
- `_NON_V0_PRESETS['Phi']['maxDaughterEta']`: 1.8 (B/C) or 2.4 (A) → 2.5

The value matches the CMS tracker acceptance edge; a bachelor track
that reconstructs beyond |η| = 2.5 has effectively no hits and is
already gone. V0-mode channels have no per-channel η cap (only the
V0Producer's kinematics apply), so no change there.

## Decision 5 — V0 quality knobs, one-shot loosening, all channels

**Not a study — a single defensible working point.**  Λb, B0→Ks, and
ψ(2S) all share one V0Producer clone
(`ALCARECOTkAlV0Candidates`), and all three are yield-limited under
`generalV0Candidates` defaults. Rather than sweep, we pick one
moderately looser point and ship it. If it turns out to be
insufficient (or too permissive) a follow-up can revisit.

**Applied values on the local clone:**

| Knob | Default | New | Rationale |
|------|---------|-----|-----------|
| `tkChi2Cut` | 10 | 10 (unchanged) | already generous |
| `vtxChi2Cut` | 6.63 | **15** | CL ≈ 0.2 % → recovers tail vertices that fit poorly under tight-vertex hypothesis but still make sense as V0's |
| `vtxDecaySigXYCut` | 15 | **10** | Λ / Ks are long-lived (cτ ≈ 7.9 cm / 2.7 cm) — displacement significance ≫ 10 easily; the cut was tightening on the low-Lxy tail where signal still lives |
| `cosThetaXYCut` | 0.998 | **0.995** | Λb → J/ψ Λ inherits B decay length + Λ decay length; the Λ pointing back to *the primary vertex* is smeared by that B flight distance, so 0.998 preferentially rejects Λb-produced Λ's. 0.995 is looser but still purity-preserving on prompt Ks |
| `tkIPSigXYCut` | 2 | **1** | daughter IP significance; softer values are nearly free on V0 daughters that are already displaced by cτ ≫ 1 sigma |

Same clone drives all three V0 channels, so the changes hit B0→Ks,
Λb, and ψ(2S) uniformly — no per-channel dispatch.

**Landing:** the values above go straight into
`ALCARECOTkAlV0Candidates_cff.py` as the new defaults. No env-var
switch, no preset selection layer, no study script. Full-data
production runs with these on. A working-point scan is filed as a
follow-up proposal *after* the first full-data batch is landed and we
have real yield numbers to argue over.

## Decision 6 — ψ(2S) → J/psi π+π- via a two-track dipion (not V0 Ks)

**Bug being fixed.** `ALCARECOTkAlJpsiXPsi2SCandidates.xSrc` in the
current cff is `cms.InputTag('ALCARECOTkAlV0Candidates', 'Kshort')`.
The V0Producer Ks output requires a Kalman-fit *displaced* π+π-
vertex (default `vtxDecaySigXYCut = 15`, i.e. Lxy significance ≳ 15
before the vertex fit is accepted; Ks cτ ≈ 2.7 cm sets the natural
scale). ψ(2S) → J/psi π+π- (BR ≈ 34.7 %) does not have such a
displaced π+π- — the two pions share the ψ(2S) decay vertex which
sits either at the primary vertex (prompt ψ(2S), ~80 % of ψ(2S)
rate at CMS) or displaced by ~B cτ ≈ 500 μm (non-prompt ψ(2S) from
B decays, ~20 %). Neither case looks like a Ks. Result: the V0-Ks
wiring rejects essentially all real ψ(2S) → J/psi π+π- signal, and
what survives is combinatoric π+π- fakes.

**Fix.** Add a new `TwoBodyDecayCandidateProducer` instance
`ALCARECOTkAlJpsiXPiPiCandidates` — architecturally identical to
`ALCARECOTkAlJpsiXPhiCandidates` (symmetric two-track producer),
just with pion daughter mass / pdgId and no vertex fit — and repoint
the ψ(2S) producer at it. Do NOT keep a V0-Ks fallback path — the
previous wiring is removed entirely (the input tag is repointed,
not augmented).

**Prompt vs displaced ψ(2S).** Both are accepted. The producer chain
applies no displacement cut anywhere:
- The dipion producer sets `applyVertexFit = False`, so the π+π-
  pair is a raw two-track combination with only a mass window.
- The downstream `JpsiXCandidateProducer` under preset B applies no
  Kalman B-vertex fit and no α_BS / Lxy cut on the ψ(2S) vertex.
This matches the treatment of the other V0-adjacent channels
(B0→Ks, Λb) under preset B, and matches preset B's overall
philosophy of raw kinematics with cuts deferred to downstream
Stage-2 code. If a downstream user wants displaced-only ψ(2S)
(i.e. the from-B-decay fraction), they can filter on the ψ(2S)
vertex vs the primary vertex themselves. A future proposal could
add a `maxMotherAlphaBS` / `minBLxyOverSigma` cut on the ψ(2S)
producer under preset C only, mirroring what preset C already does
for the four non-V0 channels.

**Combinatoric-background filter: track-track DCA.** The 200-event
smoke with the naive dipion producer (no per-pair filter beyond the
mass window) gave 91 ψ(2S) candidates per accepted event — a
combinatoric background rate that would be problematic for downstream
Stage-2 CVH CPU. The natural filter is a per-pair *distance of
closest approach* cut: two pions from the same physical vertex have
DCAs at the beamspot + tracking-resolution scale (few tens of μm),
while incoherent pairs from random tracks have DCAs broadly
distributed up to cm-scale.

`JpsiXCandidateProducer` already applies exactly this style of cut on
the (bachelor, J/psi muon) track pair via the `maxBachelorMuTrackDOCA`
knob (value 300 μm under preset B) — it uses a static straight-line
3D DCA helper `trackTrackDCA(t1, t2)`, no Kalman fit, no
`TransientTrackBuilder` ES dependency. That helper is ported into
`TwoBodyDecayCandidateProducer` behind a new optional PSet field
`maxTrackTrackDOCA` (default = infinity → no cut). The dipion
producer sets `maxTrackTrackDOCA = 0.03` cm, matching the physical
scale of the bachelor-muon DOCA cut. Other instances of
`TwoBodyDecayCandidateProducer` (K*0, φ, J/psi, J/psi-only) leave
the knob at its default so they stay bit-invariant relative to any
build without the code change.

Alternatives considered but not applied in this proposal:
- **2-track Kalman vertex fit + probability cut** (`applyVertexFit
  = True, minVtxProb = 0.005`) — a stronger filter that also handles
  z-mismatched tracks better than a straight-line DCA, but it
  contradicts preset B's "no Kalman fits at AlCaReco" philosophy
  and adds a `TransientTrackBuilder` ES dependency to a producer
  that currently does not need it.
- **Per-pion `|dxy(BS)|` or `|dxy(BS)| / σ` cut** — a track-level
  displacement filter using data already on `reco::Track`. Useful
  for *prompt vs displaced* separation but not for combinatoric
  suppression (both prompt and displaced ψ(2S) have small pair DCA;
  their individual track IPs differ). Deferred as a follow-up if
  ψ(2S) needs to be split into prompt-only / displaced-only samples.

**Dipion mass window.** m(π+π-) from ψ(2S) → J/psi π+π- is bounded
above by `m_ψ(2S) − m_J/psi ≈ 0.589 GeV` and above the 2m_π threshold
of ≈ 0.279 GeV. The measured spectrum peaks near 0.5 GeV (chiral-
dynamics enhancement at large dipion mass). Window (0.28, 0.65) covers
the physical range with a small sideband margin.

**Daughter cuts.** Match K*0 / φ under the current cff: `minDaughterPt
= 0.1`, `maxDaughterEta = 2.5`, no vertex fit. Prompt π+π- come from
the primary vertex — no per-pair Kalman displacement to fit.

**Consequence for the PDG/mass regime split.** Historically the
daughter-species dispatch had two regimes: PDG (four non-V0 channels)
and mass-on-daughter (three V0 channels including ψ(2S)). After this
fix ψ(2S) daughters carry `pdgId = ±211` from
`TwoBodyDecayCandidateProducer`, so ψ(2S) joins the PDG regime.
Only two channels remain in the V0-mass regime: B0 → J/psi Ks and
Λb → J/psi Λ. Both still consume the V0Producer clone. The daughter-
layout diagnostic (`condor/jpsix_alcareco/_diag_daughter_layout.py`)
needs its per-channel regime tables updated.

**Ordering note.** `ALCARECOTkAlJpsiXPiPiCandidates` must be inserted
into `seqALCARECOTkAlJpsiX` before `ALCARECOTkAlJpsiXPsi2SCandidates`
so the CMS framework's dependency check is satisfied.

**Why now, in this proposal.** Anything shipped under
`add-jpsi-x-muons-and-preprod-refinements` is what runs in the next
production. Discovering this before Condor launch means we don't
waste a full-2016H batch on the broken wiring; putting the fix in
the same CMSSW build means we don't need a separate production
cycle.

## Cross-item interactions

- (1) + (2): the muon collection has to include the *loose*-selector
  muons or the J/psi-only candidates have no muon-level info at
  Stage 2. The dedicated `ALCARECOTkAlJpsiXMuons` collection is
  designed as the union of tight + loose selectors' outputs.
- (3) + (4): unrelated but shipped together.
- (5) + (6): (5) loosens V0 cuts on the shared clone which feeds
  B0→Ks and Λb; (6) removes ψ(2S) from that clone entirely. Net
  effect: B0→Ks + Λb get the V0 loosening; ψ(2S) gets the correct
  prompt π+π- source. No overlap or ordering conflict.
- (5)+(6) + everything: expected to raise Λb / B0→Ks yield (5) and
  ψ(2S) yield (6); the four non-V0-Ks/Λ channels (B+, Bc, K*0, φ)
  are untouched. Cross-check on the first full-data batch that
  per-event candidate multiplicity stays sane (a runaway π+π- pair
  or V0 multiplicity would blow up event size).

## Validation scope

Preset B is the production default (locked in by
`finalize-jpsi-x-preset-b-production`). All smoke runs and integration
tests in this proposal SHALL exercise **preset B only** — the
`TKALJPSIX_SELECTION_PRESET` env var is left unset (equivalent to
`= B`) throughout the validation tasks. Presets A and C are untouched
by this proposal in their non-numeric behaviour, and the η-cut update
lands on all three presets purely by editing the shared
`_NON_V0_PRESETS` block, so preset-A / preset-C runs are not needed to
gain confidence.

## Non-goals

- No new preset. Everything lands on top of preset B / C.
- No change to Stage-2 CVH physics (only the species-dispatch code).
- No change to the Stage-1 → Stage-2 hand-off wiring.
- No change to the tight muon selector's cuts (loose selector is a
  new sibling, not a replacement).
