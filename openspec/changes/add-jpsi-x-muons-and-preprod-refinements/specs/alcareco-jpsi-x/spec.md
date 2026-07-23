## ADDED Requirements

### Requirement: Persist Selected Muons in the AlCaReco Stream

The `ALCARECOTkAlJpsiX` stream SHALL persist, on every accepted event,
a `reco::MuonCollection` named `ALCARECOTkAlJpsiXMuons` containing every
`reco::Muon` that was accepted by any Stage-0 J/psi-X muon selector on
that event. Muons accepted by more than one selector SHALL appear at
most once in the collection (by `edm::Ptr` identity, not by track
identity). The collection SHALL be added to
`OutALCARECOTkAlJpsiX_noDrop.outputCommands` (`keep
*_ALCARECOTkAlJpsiXMuons_*_*`) and SHALL NOT depend on the
`TKALJPSIX_SELECTION_PRESET` value (the set of muons that pass the
per-channel selectors depends on the preset, but the persisted
collection SHALL be the union of what any selector emits — never a
subset).

#### Scenario: Event with tight and loose muons

- **WHEN** an event contains 3 `reco::Muon`s, 2 of which pass the
  existing tight selector (`isGlobalMuon & isTrackerMuon & ...`) and
  all 3 of which pass the new loose selector
  (`isGlobalMuon | isTrackerMuon & ...`)
- **THEN** `ALCARECOTkAlJpsiXMuons` on that event contains 3 muons
  (the union), each appearing exactly once, with `edm::Ptr` provenance
  pointing back to the original input `reco::MuonCollection`

#### Scenario: Event with no muons passing either selector

- **WHEN** an event has no muon passing any of the two selectors
- **THEN** the event is rejected by the Stage-0 filter cascade
  BEFORE any AlCaReco module runs, so `ALCARECOTkAlJpsiXMuons` is
  never emitted (the standard `filter=cms.bool(True)` behaviour is
  preserved)

---

### Requirement: Track-to-Muon Association in the AlCaReco Stream

The `ALCARECOTkAlJpsiX` stream SHALL persist, on every accepted event,
an `edm::Association<reco::MuonCollection>` named
`ALCARECOTkAlJpsiXTrackToMuon` keyed on the deduplicated persisted
track collection `ALCARECOTkAlJpsiX` and valued into the persisted
muon collection `ALCARECOTkAlJpsiXMuons` (see requirement above). For
every track `t` in `ALCARECOTkAlJpsiX`, the association SHALL return a
non-null `reco::MuonRef` if `t` is the `innerTrack()` of some muon
in `ALCARECOTkAlJpsiXMuons`, and a null ref otherwise. The association
SHALL be added to `OutALCARECOTkAlJpsiX_noDrop.outputCommands` (`keep
*_ALCARECOTkAlJpsiXTrackToMuon_*_*`).

#### Scenario: J/psi daughter track resolves to its muon

- **WHEN** a B+ candidate's J/psi daughter has a positive-charge muon
  daughter whose `RecoChargedCandidate::track()` returns a `TrackRef`
  into `ALCARECOTkAlJpsiX`
- **THEN** `ALCARECOTkAlJpsiXTrackToMuon->at(that TrackRef)` returns a
  non-null `MuonRef` into `ALCARECOTkAlJpsiXMuons` whose `innerTrack()`
  (dereferenced via the `edm::Ptr` provenance chain into the original
  `generalTracks`) corresponds to the same physical track

#### Scenario: Bachelor kaon track resolves to null

- **WHEN** a B+ candidate's bachelor kaon has a `TrackRef` into
  `ALCARECOTkAlJpsiX` and is not the inner track of any muon in
  `ALCARECOTkAlJpsiXMuons`
- **THEN** `ALCARECOTkAlJpsiXTrackToMuon->at(that TrackRef)` returns a
  null `MuonRef`

---

### Requirement: Loose Muon Selector for the Production J/psi-Only Channel

The AlCaReco stream SHALL expose a second muon selector
`TkAlLooseIdMuonSelector` (in
`Alignment/CommonAlignmentProducer/python/TkAlMuonSelectors_cfi.py`)
with the exact cut string
`(isGlobalMuon | isTrackerMuon) & abs(eta) < 2.5 &
numberOfMatches > 1`. The `globalTrack.*` sub-cuts from
`TkAlGoodIdMuonSelector` (namely
`globalTrack.hitPattern.numberOfValidMuonHits > 0` and
`globalTrack.normalizedChi2 < 20`) SHALL be omitted, because a
tracker-only muon (accepted by the OR clause) carries no `globalTrack`
and those sub-cuts would misbehave under `edm::View<reco::Muon>`
expression evaluation. The loose selector SHALL emit a superset
(by muon identity) of what the tight selector `TkAlGoodIdMuonSelector`
emits on every event. `ALCARECOTkAlJpsiX_cff.py` SHALL instantiate
this selector as `ALCARECOTkAlJpsiXLooseMuons`. It feeds the new
production J/psi-only channel (see the "J/psi-Only Per-Channel Path"
requirement) on every production event; the same-event tight-vs-loose
J/psi comparison this makes available is a side benefit that a
follow-up proposal MAY use to consider globalising the loose selector,
but the loose selector and the J/psi-only channel SHALL ship in
production regardless of any such follow-up.

#### Scenario: A tracker-only muon reaches the loose selector

- **WHEN** an event has a muon with `isGlobalMuon = false`,
  `isTrackerMuon = true`, `|η| < 2.5`, and `numberOfMatches = 2`
- **THEN** the muon passes `ALCARECOTkAlJpsiXLooseMuons` and does NOT
  pass `ALCARECOTkAlJpsiXGoodMuons`

#### Scenario: A tracker-only muon with numberOfMatches = 1 is rejected by both

- **WHEN** an event has a muon with `isGlobalMuon = false`,
  `isTrackerMuon = true`, `|η| < 2.5`, and `numberOfMatches = 1`
- **THEN** the muon passes neither `ALCARECOTkAlJpsiXLooseMuons` nor
  `ALCARECOTkAlJpsiXGoodMuons` (both require `numberOfMatches > 1`)

#### Scenario: The loose selector is a superset of the tight selector

- **WHEN** on any event `M_tight = ALCARECOTkAlJpsiXGoodMuons` and
  `M_loose = ALCARECOTkAlJpsiXLooseMuons`
- **THEN** `M_tight ⊆ M_loose` by `edm::Ptr` identity

---

### Requirement: J/psi-Only Production Channel

The AlCaReco stream SHALL emit a per-event
`VertexCompositeCandidateCollection`
`ALCARECOTkAlJpsiXJpsiOnlyResonances` on every accepted production
event, produced from a dedicated `TwoBodyDecayCandidateProducer`
(`ALCARECOTkAlJpsiXJpsiOnlyCandidates`) with
`muonSrc = ALCARECOTkAlJpsiXLooseMuons`, mass window
[2.95, 3.25] GeV, `daughterMass = 0.105`, `daughterPdgId = 13`,
`motherPdgId = 443`, and `applyVertexFit = False`. Its daughter track
refs SHALL be re-keyed onto the deduplicated `ALCARECOTkAlJpsiX`
collection via the existing `VertexCompositeCandidateRemapper`
mechanism (same as the other seven channels). This channel is an
eighth production channel on par with the existing seven; it is not a
smoke-only artifact and SHALL be included in every production build
that ships this proposal's changes. The producer SHALL be part of
`seqALCARECOTkAlJpsiX` and SHALL NOT depend on
`TKALJPSIX_SELECTION_PRESET`. The output SHALL be persisted:
`keep *_ALCARECOTkAlJpsiXJpsiOnlyResonances_*_*`.

#### Scenario: J/psi-only channel is populated when tight-selector J/psi is empty

- **WHEN** an event has two tracker-only muons forming a J/psi in the
  [2.95, 3.25] mass window, but zero muons pass the tight selector
- **THEN** `ALCARECOTkAlJpsiXJpsiOnlyResonances` on that event has at
  least one candidate, while `ALCARECOTkAlJpsiXJpsiCandidates` (the
  tight-selector J/psi feeding the seven bachelor channels) is empty

#### Scenario: The six bit-invariant existing channels

- **WHEN** the same event is processed with the J/psi-only channel
  wired in vs a pre-change reference build without it
- **THEN** the emitted candidate contents of the six existing
  channels that are NOT touched by this proposal's other items
  (`BPlus`, `B0Kstar`, `B0Ks`, `BsPhi`, `Lambdab`, `Bc`) are
  bit-identical between the two builds. `Psi2S` is EXCLUDED from
  the bit-invariance scenario because item 6 (Prompt Dipion
  Producer) intentionally changes its source, and B+/Bc/K*0/φ
  candidate yields may increase due to the η bump (item 4)

---

### Requirement: Loosened V0 Selection on the Local V0Producer Clone

`ALCARECOTkAlV0Candidates_cff.py` SHALL override the following
`generalV0Candidates` defaults on the `ALCARECOTkAlV0Candidates`
clone, in addition to the already-set `tkPtCut = 0.1`:
- `vtxChi2Cut = 15` (default 6.63)
- `vtxDecaySigXYCut = 10` (default 15)
- `cosThetaXYCut = 0.995` (default 0.998)
- `tkIPSigXYCut = 1` (default 2)
- `tkChi2Cut = 10` (unchanged)

These values SHALL be applied unconditionally — there is no env-var
switch, no preset variant, no per-channel override. Two V0-mode
channels (`ALCARECOTkAlJpsiXB0KsCandidates`,
`ALCARECOTkAlJpsiXLambdabCandidates`) consume the clone via its
`Kshort` / `Lambda` output labels, so the loosening reaches both
together. The former V0-mode ψ(2S) channel MOVES OFF the V0 path
under the "Prompt π+π- Producer" requirement below and is therefore
no longer affected by V0Producer cuts. A future proposal MAY
revisit the loosened values with a proper working-point scan; this
proposal explicitly does NOT run one.

#### Scenario: V0 yield rises on the two remaining V0-mode channels

- **WHEN** a fixed 100-file 2016H Charmonium input is processed with
  the loosened V0 clone vs a reference build with the same code
  minus the V0 overrides
- **THEN** the number of `ALCARECOTkAlJpsiXB0KsResonances` and
  `ALCARECOTkAlJpsiXLambdabResonances` candidates each rises
  relative to the reference build (yield gain is monitored, not
  gated — precise magnitude is a follow-up-proposal question)

#### Scenario: A Λb candidate whose Λ has vertex chi2/ndf between 6.63 and 15 now survives

- **WHEN** a Λb candidate has a Λ V0 fit with normalised vertex χ²
  in the range (6.63, 15] and every other cut passing at default
- **THEN** the Λ passes `ALCARECOTkAlV0Candidates`, and downstream
  `ALCARECOTkAlJpsiXLambdabCandidates` receives it as an input

---

### Requirement: Dipion Producer for the ψ(2S) Channel

The AlCaReco stream SHALL build ψ(2S) → J/psi π+π- candidates from
a two-track π+π- pair (a new `TwoBodyDecayCandidateProducer` instance
`ALCARECOTkAlJpsiXPiPiCandidates`) whose only geometric requirement is
that the pair fall in the physical dipion mass window — no Kalman
vertex fit, no displacement significance cut. The previous wiring
that fed the ψ(2S) producer from `ALCARECOTkAlV0Candidates:Kshort`
(a V0Producer output with Ks-quality displaced-vertex cuts) SHALL be
removed entirely; there is no fallback path to the V0 Ks input.

Physics rationale for accepting both prompt and non-prompt ψ(2S):
the two pions from ψ(2S) → J/psi π+π- (BR ≈ 34.7 %) always share
the ψ(2S) vertex regardless of how the ψ(2S) itself was produced
(direct ψ(2S) production at the primary vertex ≈ 80 % of the total
ψ(2S) rate at LHC energies; non-prompt ψ(2S) from B decays ≈ 20 %,
displaced by ~B cτ). Preset B does not apply Kalman-fit-based
displacement cuts anywhere in the stream, so both populations are
accepted; downstream analyses that want to isolate the non-prompt
(from-B) fraction can filter on the ψ(2S) vertex-to-PV displacement
themselves.

The dipion producer SHALL be configured with:
- `src = 'generalTracks'`, `muonSrc = ''` (no muon filter)
- `daughterMass = 0.139570` (pion), `daughterPdgId = 211`
- `motherPdgId = 100443` (ψ(2S) — used only as a candidate tag; the
  actual mother mass is reconstructed by the downstream
  `JpsiXCandidateProducer` from the J/psi + dipion sum)
- `minMass = 0.28`, `maxMass = 0.65` — physical dipion mass window
  (the ψ(2S) → J/psi π+π- dipion mass distribution peaks near
  0.5-0.58 GeV and cuts off at `m_ψ(2S) − m_J/psi ≈ 0.589 GeV`)
- `applyChargeFilter = True`, `charge = 0`, `useUnsignedCharge = True`
- `minDaughterPt = 0.1`, `maxDaughterEta = 2.5` (matches K*0 / φ)
- `applyVertexFit = False` (no Kalman fit at AlCaReco — the π+π-
  come from the ψ(2S) vertex which may itself be prompt or displaced;
  either way there is no per-pair *additional* displacement to fit)
- `maxTrackTrackDOCA = 0.03` cm (300 μm) — a static geometric
  distance-of-closest-approach cut between the two pion tracks,
  used to suppress the combinatoric background rate. Real π+π-
  from a common vertex (whether prompt at the PV or displaced at a
  B decay point) have DCAs at the beamspot + track-resolution
  scale (few tens of μm), well below the cut. Random-track
  combinatorics have DCAs broadly distributed up to cm-scale.
  Matches the physical scale of the `maxBachelorMuTrackDOCA` cut
  used on the four non-V0 preset-B `JpsiXCandidateProducer`
  instances. This filter does NOT bias prompt vs displaced ψ(2S)
  (both are single-vertex topologies with small pair DCA).

Additionally, `TwoBodyDecayCandidateProducer` SHALL accept a
`maxTrackTrackDOCA` PSet field (double, default = infinity / no
cut). When set, the producer SHALL compute a static straight-line
3D DCA between the two daughter tracks using the same
`trackTrackDCA(t1, t2)` primitive already used by
`JpsiXCandidateProducer` (i.e. no Kalman fit, no
`TransientTrackBuilder` dependency), and SHALL drop any pair whose
DCA exceeds the threshold. The `nDroppedByDOCA` counter SHALL be
incremented on each such drop and logged at end-of-job alongside
`nDroppedByVtx`. Other instances of `TwoBodyDecayCandidateProducer`
in the stream (K*0, φ, J/ψ, J/ψ-only) SHALL NOT set this field, so
their behaviour is bit-invariant under the code change.

`ALCARECOTkAlJpsiXPsi2SCandidates.xSrc` SHALL be
`cms.InputTag('ALCARECOTkAlJpsiXPiPiCandidates')` instead of
`cms.InputTag('ALCARECOTkAlV0Candidates', 'Kshort')`. The ψ(2S)
mother mass window (3.5, 3.9) SHALL be unchanged.
`ALCARECOTkAlJpsiXPiPiCandidates` SHALL be inserted into
`seqALCARECOTkAlJpsiX` before `ALCARECOTkAlJpsiXPsi2SCandidates` so
the producer graph is well-ordered. The pdgId invariant applies:
every ψ(2S) daughter π now carries `pdgId = ±211` (from
`TwoBodyDecayCandidateProducer`), so ψ(2S) moves from the historical
V0-mass regime into the PDG regime for downstream Stage-2 species
dispatch. Daughter-layout diagnostics that categorize channels by
regime SHALL be updated accordingly.

#### Scenario: A prompt-produced ψ(2S) is kept

- **WHEN** an event contains a J/psi + two pions where the (J/psi +
  π+π-) invariant mass lies in (3.5, 3.9), the π+π- invariant mass
  lies in (0.28, 0.65), and the ψ(2S) vertex is at the primary
  vertex (i.e. directly produced prompt ψ(2S), no B-parent flight)
- **THEN** the event yields at least one candidate in
  `ALCARECOTkAlJpsiXPsi2SResonances`. Under the pre-change V0-Ks
  wiring the same event would have been rejected because the
  primary-vertex π+π- fail the V0Producer Ks vertex-quality and
  displacement-significance cuts

#### Scenario: A non-prompt ψ(2S) from a B decay is also kept

- **WHEN** an event contains a J/psi + two pions where the (J/psi +
  π+π-) invariant mass lies in (3.5, 3.9), the π+π- invariant mass
  lies in (0.28, 0.65), and the ψ(2S) vertex is displaced from the
  PV by ~B cτ ≈ 500 μm (i.e. non-prompt ψ(2S) from a B → ψ(2S) X
  decay)
- **THEN** the event yields at least one candidate in
  `ALCARECOTkAlJpsiXPsi2SResonances`. The preset-B ψ(2S) producer
  applies no displacement cut, so both prompt and non-prompt ψ(2S)
  populations are represented in the output

#### Scenario: A displaced Ks + J/psi is now rejected by the ψ(2S) channel

- **WHEN** an event contains a J/psi + a displaced Ks (long Lxy;
  Ks-quality vertex fit passes) whose (J/psi + Ks) invariant mass
  happens to lie in the ψ(2S) window (3.5, 3.9)
- **THEN** the event MAY still yield a `B0Ks` candidate (via the
  `B0Ks` producer, which still consumes the V0 Ks) but SHALL NOT
  yield a `Psi2S` candidate, because the ψ(2S) producer no longer
  sources π+π- from the V0Ks output

#### Scenario: Random π+π- pair rejected by track-track DOCA

- **WHEN** two unrelated tracks in the same event happen to satisfy
  the (0.28, 0.65) dipion mass window but come from geometrically
  distinct origins (their straight-line 3D DCA exceeds 300 μm)
- **THEN** the pair is rejected by
  `ALCARECOTkAlJpsiXPiPiCandidates` before it can seed a ψ(2S)
  candidate, and the `nDroppedByDOCA` counter is incremented

#### Scenario: Real ψ(2S) dipion survives track-track DOCA

- **WHEN** a real ψ(2S) → J/psi π+π- decay in the event yields two
  pion tracks that share the ψ(2S) vertex (whether prompt at PV or
  displaced from a B decay); their 3D DCA is at the beamspot +
  tracking resolution scale (few tens of μm)
- **THEN** the pair passes the `maxTrackTrackDOCA = 0.03` cm cut
  and is kept in `ALCARECOTkAlJpsiXPiPiCandidates`

#### Scenario: K*0 / phi / J/psi / J/psi-only producers unaffected by the DOCA knob

- **WHEN** an event is processed with the code that supports the new
  `maxTrackTrackDOCA` PSet field vs a hypothetical baseline without
  the field
- **THEN** the candidate contents of `ALCARECOTkAlJpsiXKstarCandidates`,
  `ALCARECOTkAlJpsiXPhiCandidates`,
  `ALCARECOTkAlJpsiXJpsiCandidates`, and
  `ALCARECOTkAlJpsiXJpsiOnlyCandidates` are bit-identical between the
  two builds (these instances leave the knob at its default = no cut)

#### Scenario: ψ(2S) daughters carry non-zero pdgId

- **WHEN** a `Psi2S` candidate is emitted after this change
- **THEN** `daughter(1)` (the π+π- VCC) has `daughter(0)->pdgId()`
  and `daughter(1)->pdgId()` both equal to ±211 (from
  `TwoBodyDecayCandidateProducer.cc:239-240`), not 0 as under the
  V0-Ks wiring where the pions carried `pdgId==0` and mass≈0.1396

## MODIFIED Requirements

### Requirement: Bachelor and X-Daughter Pseudorapidity Cuts

The per-channel pseudorapidity caps SHALL be set to **2.5** under
every preset (A, B, C) — up from 1.8 under B/C and from 2.4 under A —
on all four non-V0 channels, specifically
`_NON_V0_PRESETS['BPlus']['maxBachelorEta']`,
`_NON_V0_PRESETS['Bc']['maxBachelorEta']`,
`_NON_V0_PRESETS['Kstar']['maxDaughterEta']`, and
`_NON_V0_PRESETS['Phi']['maxDaughterEta']`. This is the same value
across all three presets so the η boundary no longer varies with the
`TKALJPSIX_SELECTION_PRESET` switch. V0-mode channels (`B0Ks`,
`Lambdab`, `Psi2S`) have no per-channel `maxBachelorEta` /
`maxDaughterEta` in the cff (V0Producer kinematics dominate), so they
are unchanged. The user-facing comment block at
`ALCARECOTkAlJpsiX_cff.py:79-80` SHALL be updated to reflect the new
2.5 value.

#### Scenario: A B+ candidate with |η(K)| = 2.3 is kept under preset B

- **WHEN** an event contains a B+ candidate whose bachelor kaon has
  reconstructed `|η| = 2.3`, all other preset-B cuts pass, and preset
  B is active (`TKALJPSIX_SELECTION_PRESET` unset or `= B`)
- **THEN** the B+ candidate survives and is written to
  `ALCARECOTkAlJpsiXBPlusResonances`

#### Scenario: A B+ candidate with |η(K)| = 2.6 is rejected under every preset

- **WHEN** an event contains a B+ candidate whose bachelor kaon has
  reconstructed `|η| = 2.6` and any preset is active
- **THEN** the B+ candidate is rejected by the
  `maxBachelorEta = 2.5` cut

---

### Requirement: Daughter PDG-Id as the Species Source of Truth

The AlCaReco J/psi-X producers SHALL emit every `RecoChargedCandidate`
daughter with a **non-zero** `pdgId()` reflecting its assigned mass
hypothesis. This applies to `JpsiXCandidateProducer` (bachelor tracks
in track-mode; X sub-resonance daughters in vcc-mode) and to
`TwoBodyDecayCandidateProducer` (both symmetric and asymmetric modes). The producer ctors
SHALL reject (`edm::Exception`) any configuration where the relevant
PDG-id PSet field (`bachelorPdgId`, `firstDaughterPdgId`,
`secondDaughterPdgId`, or `daughterPdgId`) is 0. Downstream consumers
(specifically `JpsiKCandidateSplitter` and the Stage-2 CVH refit
config) SHALL use `daughter->pdgId()` to look up the mass hypothesis
for the invariant-mass / vertex-fit computation, instead of relying on
the input collection name. Daughter ORDERING within a candidate
(`daughter(0)` = J/psi; `daughter(1)` = bachelor / X) SHALL still be
determined by the collection layout — only the species tag SHALL move
to the pdgId.

#### Scenario: B+ bachelor track carries pdgId = ±321

- **WHEN** a B+ candidate is emitted by
  `ALCARECOTkAlJpsiXBPlusCandidates` with a `+` charge bachelor
- **THEN** `candidate.daughter(1)->pdgId() == -321` (positive-charge =
  anti-particle in the existing sign convention)

#### Scenario: K*0 candidate carries mixed pdgId on its daughters

- **WHEN** a K*0 candidate is emitted by
  `ALCARECOTkAlJpsiXKstarCandidates` with `swapped = false`
  (positive-charge = kaon, negative-charge = pion)
- **THEN** `daughter(0)->pdgId() = -321` and
  `daughter(1)->pdgId() = +211` (matching the existing sign convention
  in `TwoBodyDecayCandidateProducer.cc:239-240`)

#### Scenario: Downstream Stage-2 mass hypothesis derived from pdgId

- **WHEN** `JpsiKCandidateSplitter` receives a B+ candidate
- **THEN** the bachelor mass hypothesis passed to the Stage-2
  ResidualGlobalCorrectionMaker is derived from
  `daughter(1)->pdgId()` (`abs == 321 => 0.493677 GeV`, etc.) and NOT
  from a hard-coded `bachelorMass` in the splitter cfi

#### Scenario: Configuration with pdgId = 0 is rejected at ctor time

- **WHEN** an `ALCARECOTkAlJpsiXBPlusCandidates` PSet is constructed
  with `bachelorPdgId = cms.int32(0)`
- **THEN** the producer's constructor throws an `edm::Exception`
  during cmsRun configuration and the job fails fast
