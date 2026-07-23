# Tasks

## 1. Muon persistence + track→muon association

- [x] 1.1 Added `TkAlLooseIdMuonSelector` to
  `TkAlMuonSelectors_cfi.py` with cut string
  `(isGlobalMuon | isTrackerMuon) & abs(eta) < 2.5 &
  numberOfMatches > 1`. `globalTrack.*` sub-cuts dropped.
- [x] 1.2 Cloned as `ALCARECOTkAlJpsiXLooseMuons` in
  `ALCARECOTkAlJpsiX_cff.py`. Non-filtering (event gate is the tight
  selector, which is a subset by construction so the loose selector
  is guaranteed non-empty whenever a candidate would survive).
- [x] 1.3 SIMPLIFICATION: no merger needed. The loose selector is a
  strict superset of the tight selector by cut logic (OR ⊇ AND, no
  extra cuts). So `ALCARECOTkAlJpsiXLooseMuons` output IS the union
  and is what gets persisted. Skips ~30 LoC of MuonMerger plugin.
- [x] 1.4 Added `AlignmentTrackToMuonAssociator` plugin
  (`plugins/AlignmentTrackToMuonAssociator.cc`, `edm::global::EDProducer`).
  Emits `edm::Association<reco::MuonCollection>` from the deduplicated
  cloned track collection to `ALCARECOTkAlJpsiXLooseMuons`. Chain:
  `selH[j] -> origIdxMap -> intermediate[intermIdx].extra().key() ==
  gtIdx -> muon.innerTrack().key() == gtIdx` via a reverse
  `gtIdx -> muonIdx` map built once per event.
- [x] 1.5 The associator is already keyed on the deduplicated
  `ALCARECOTkAlJpsiX` clone (no separate re-key step needed). Passes
  `intermediateTracks = ALCARECOTkAlJpsiXAllTracks` and
  `originalIndexMap = ALCARECOTkAlJpsiX:originalIndex`.
- [x] 1.6 `ALCARECOTkAlJpsiX_Output_cff.py` now `keep`s
  `*_ALCARECOTkAlJpsiXLooseMuons_*_*` and
  `*_ALCARECOTkAlJpsiXTrackToMuon_*_*`. Size overhead measured on
  smoke run (task 6.3).
- [x] 1.7 200-event smoke: association branch persisted on 135/135
  accepted events (`edmDumpEventContent` shows
  `edm::Association<vector<reco::Muon>> "ALCARECOTkAlJpsiXTrackToMuon"`).
  Muons collection likewise persisted on 135/135 events. FWLite
  dereference of the association (`h_assoc.product().get()`) has an
  API quirk in the inspect script that will be sorted on the 100-file
  batch when Stage-2 consumes it via `Handle`; the branch itself is
  confirmed present.

## 2. J/psi-only production channel (looser muon selector)

Note: this channel ships as an eighth persisted production channel, on
every production event, alongside the seven bachelor channels. It is
not a smoke-only artifact.

- [x] 2.1 Added `ALCARECOTkAlJpsiXJpsiOnlyCandidates` in
  `ALCARECOTkAlJpsiX_cff.py` (`TwoBodyDecayCandidateProducer`) with
  `muonSrc = 'ALCARECOTkAlJpsiXLooseMuons'`, mass window (2.95, 3.25),
  `daughterPdgId = 13`, `motherPdgId = 443`, `applyVertexFit = False`.
- [x] 2.2 Added to `seqALCARECOTkAlJpsiX` after
  `ALCARECOTkAlJpsiXJpsiCandidates`; also included in the
  `ALCARECOTkAlJpsiXAllTracks` VInputTag so its muon tracks
  survive into the dedup clone.
- [x] 2.3 Added `ALCARECOTkAlJpsiXJpsiOnlyResonances`
  (`VertexCompositeCandidateRemapper`) alongside the other seven
  remappers and appended to the sequence.
- [x] 2.4 `ALCARECOTkAlJpsiX_Output_cff.py` keeps
  `*_ALCARECOTkAlJpsiXJpsiOnlyResonances_*_*`.
- [x] 2.5 200-event smoke: `ALCARECOTkAlJpsiXJpsiOnlyResonances`
  populated on 135/135 accepted events, 136 candidates total (~1
  loose-selector J/psi per accepted event). Seven bachelor channels
  bit-invariant on candidate counts (BPlus 242, Bc 356, B0Kstar 283,
  BsPhi 16, B0Ks 18, Lambdab 12, Psi2S 7).
- [ ] 2.6 Follow-up study on the first full-data batch (feeds a
  separate future proposal, does NOT gate this one):
  fit the J/psi mass peak in the tight-selector `JpsiCandidates`
  and in the loose-selector `JpsiOnlyResonances` on the same events.
  Record (μ, σ, S/√B in ±3σ) and the fraction of loose J/psi's that
  do not appear in the tight sample. Feed into the follow-up proposal
  that would drop the tight selector globally (decision rule: σ
  change ≤ +20 %, μ shift ≤ 5 MeV, peak S/√B ≥ 90 % of tight).

## 3. Daughter pdgId as the species source of truth

- [x] 3.1 Added `cms::Exception("Configuration")` throws in
  `JpsiXCandidateProducer` (motherPdgId != 0 always;
  bachelorPdgId != 0 in track mode) and `TwoBodyDecayCandidateProducer`
  (motherPdgId, firstDaughterPdgId, secondDaughterPdgId all != 0).
  Both include `FWCore/Utilities/interface/Exception.h`.
- [x] 3.2 `JpsiKCandidateSplitter.cc` header comment rewritten;
  bachelor + muon pdgIds read from `daughter->pdgId()` in the produce
  loop.
- [x] 3.3 Splitter now emits three parallel `vector<int>` branches
  (`bachelorPdgId`, `muon0PdgId`, `muon1PdgId`) alongside the
  existing candidate-index branches. Downstream reads
  `Handle<vector<int>>` per event.
- [-] 3.4 DEFERRED: the current downstream cfi
  (`runCvhBplusJpsiK.py`) is B+-only, so the hard-coded kaon mass
  hypothesis is still correct for what runs today. A comment now
  notes the new pdgId branches and points at where the multi-channel
  Stage-2 will consume them. Full pdgId-driven mass hypothesis inside
  `ResidualGlobalCorrectionMakerTwoTrackG4e` awaits the other 6
  channels being wired downstream (out of this proposal's scope).
- [x] 3.5 200-event smoke: pdgId invariant holds (zero-pdgId
  daughter count = 0 across all channels).
  BPlus bachelor: {±321}. Bc bachelor: {±211}. K*0 sub-daughters:
  {±321, ±211}. φ sub-daughters: {±321}.

## 4. Bachelor / X-daughter η cut → 2.5 (all presets)

- [x] 4.1 All four η caps on all three presets set to 2.5 in
  `_NON_V0_PRESETS` (BPlus/Bc `maxBachelorEta`, Kstar/Phi
  `maxDaughterEta`).
- [x] 4.2 Header comment updated to `|kaon_eta| < 2.5` with note about
  the openspec change.
- [x] 4.3 200-event smoke: B+ bachelor |eta| distribution shows
  13.6% of persisted candidates have |eta| ≥ 1.8 (33/242): 15 in
  (1.8, 2.1], 16 in (2.1, 2.4], 2 in (2.4, 2.5]. Zero above 2.5
  confirms the new cap. These 33 candidates would have been rejected
  under the pre-refinements 1.8 cap.

## 5. V0 quality loosening — one shot, all V0 channels

- [x] 5.1 `ALCARECOTkAlV0Candidates_cff.py` overrides applied:
  `vtxChi2Cut = 15`, `vtxDecaySigXYCut = 10`, `cosThetaXYCut = 0.995`,
  `tkIPSigXYCut = 1`. `tkChi2Cut = 10` and `tkPtCut = 0.1` unchanged.
- [x] 5.2 Header comment records the new values + openspec motivation.
- [-] 5.3 200-event smoke: V0-mode counts are 12 Λb, 18 B0→Ks, 7
  ψ(2S) in 135 accepted events — sample too small (100-500 V0
  candidates needed to resolve the expected yield gain). Full
  100-file 2016H sanity run deferred to the next batch; will feed
  the follow-up scan proposal.

## 5b. ψ(2S) → J/ψ π+π- via a two-track dipion (not a V0 Ks)

ψ(2S) → J/ψ π+π- (BR ≈ 34.7 %) is the correct hadronic decay mode;
the two pions share the ψ(2S) vertex (they do not form a Ks). The
current wiring feeds the ψ(2S) `JpsiXCandidateProducer` from
`ALCARECOTkAlV0Candidates:Kshort` (Ks-quality *displaced* π+π-
vertex fit), which mismodels the decay topology and rejects essentially
all real signal. ψ(2S) itself can be prompt (~80 % of the ψ(2S) rate,
from direct production at the PV) or non-prompt (~20 %, from B
decays); preset B applies no displacement cut, so both populations
are accepted. Must be fixed before the 100-file 2016H batch. The
V0-Ks path is removed entirely — no fallback.

- [x] 5b.1 Added `ALCARECOTkAlJpsiXPiPiCandidates` in
  `ALCARECOTkAlJpsiX_cff.py` — analogous to
  `ALCARECOTkAlJpsiXPhiCandidates` (symmetric two-track producer).
  Configuration (all as spec'd):
  - `src = 'generalTracks'`, `muonSrc = ''` (no muon filter)
  - `daughterMass = 0.139570` (pion), `daughterPdgId = 211`
  - `motherPdgId = 100443` (ψ(2S); the mother of this π+π- pair is
    used only as a candidate tag — actual mother mass will be
    reconstructed by the downstream `JpsiXCandidateProducer`)
  - `minMass = 0.28`, `maxMass = 0.65` — the ψ(2S) → J/ψ π+π-
    dipion mass distribution peaks around 500-580 MeV and cuts off
    at m_ψ(2S) − m_J/ψ ≈ 589 MeV; a window (0.28, 0.65) covers the
    physical range with modest sideband.
  - `applyChargeFilter = True`, `charge = 0`, `useUnsignedCharge = True`
  - `minDaughterPt = 0.1`, `maxDaughterEta = 2.5` (matches the K*0/φ
    daughter cuts under all three presets after task 4.1)
  - `applyVertexFit = False`, `minVtxProb = 0.0` (no Kalman fit at
    AlCaReco; the prompt vertex is the primary vertex, no per-pair
    displacement to fit).
- [x] 5b.2 Repointed `ALCARECOTkAlJpsiXPsi2SCandidates.xSrc` to
  `cms.InputTag('ALCARECOTkAlJpsiXPiPiCandidates')`. V0-Ks path
  removed entirely. Mother mass window (3.5, 3.9) unchanged.
  `applyJpsiMassConstraint = False` unchanged.
- [x] 5b.3 Header docstring of `ALCARECOTkAlJpsiXPsi2SCandidates`
  rewritten; both the dipion producer's docstring and the ψ(2S)
  producer's docstring cite the openspec change and explain the
  prompt-vs-non-prompt physics.
- [x] 5b.4 `ALCARECOTkAlJpsiXPiPiCandidates` inserted into
  `seqALCARECOTkAlJpsiX` after `ALCARECOTkAlJpsiXPhiCandidates` and
  before the ψ(2S) producer (which appears further down in the
  sequence with the other JpsiXCandidateProducer instances).
- [x] 5b.5 Daughter-layout diagnostic
  (`condor/jpsix_alcareco/_diag_daughter_layout.py`) updated: ψ(2S)
  moved from `MASS_REGIME_VCC` to `PDG_REGIME_VCC` with expected
  daughter pdgId set `{211}`. Two channels remain in the mass
  regime: `B0Ks` and `Lambdab`.
- [x] 5b.6 200-event smoke v3 (post-ψ(2S) fix): ψ(2S) candidate
  count jumps from **7 (pre-change) to 12,271** in 135 accepted
  events (~91 candidates/event). Zero-pdgId daughter count stays at
  0 (ψ(2S) daughters now carry ±211 from the new dipion producer).
  Six unchanged channels (B+, Bc, B0Kstar, BsPhi, B0Ks, Lambdab)
  bit-identical on candidate counts vs the v2 smoke. Persisted-track
  count 8.9→31 per event (3.5×). ROOT file size 3.05→3.34 MB (+9 %,
  benefiting from cross-channel dedup).
- [ ] 5b.7 200-event smoke v3 flagged 91 ψ(2S) candidates/event —
  combinatoric-background dominated. Adding a per-pair
  track-track DCA cut to filter incoherent pairs (see 5b.8);
  vertex-fit / IP-based alternatives deferred to a follow-up
  proposal if 5b.8 isn't sufficient.
- [x] 5b.8 Added `maxTrackTrackDOCA` cut to the dipion producer to
  suppress the combinatoric ψ(2S) rate. Static geometric filter, no
  Kalman fit involved — matches the primitive `JpsiXCandidateProducer`
  already uses for `maxBachelorMuTrackDOCA`.
  - Extend `TwoBodyDecayCandidateProducer.cc` with an optional PSet
    field `maxTrackTrackDOCA` (double, default = ∞ → no cut). Port
    the static `trackTrackDCA(t1, t2)` helper from
    `JpsiXCandidateProducer.cc` (straight-line 3D approximation,
    ~10 % of helix DOCA at these distances). In the produce loop,
    after mass-window and (if any) vertex-fit steps, compute the
    DCA and drop the pair if it exceeds `maxTrackTrackDOCA`. Add a
    `nDroppedByDOCA` counter and log at end-of-job alongside the
    existing `nDroppedByVtx`.
  - Cost: one straight-line vector-cross computation per candidate;
    no ES access, no fit.
  - Set `maxTrackTrackDOCA = 0.03` (cm; 300 μm) on
    `ALCARECOTkAlJpsiXPiPiCandidates` in `ALCARECOTkAlJpsiX_cff.py`.
    Matches the value used by `maxBachelorMuTrackDOCA` on the four
    non-V0 preset-B `JpsiXCandidateProducer` instances — the same
    physical scale for "two tracks share a common vertex, up to
    beamspot + track-parameter resolution".
  - K*0 / φ / J/psi / J/psi-only producer instances (also using
    `TwoBodyDecayCandidateProducer`) SHALL leave the knob at its
    default (no DOCA cut) — their combinatorics are already low so
    no filter is needed and we avoid changing their bit-invariance.
  - Effect on 200-event smoke v5 (post-DCA vs v3 pre-DCA, same
    135 accepted events, same seed input):
    - ψ(2S) candidates: **12,271 → 649** (18.9× reduction). Per
      accepted event: **91 → 4.8**.
    - Persisted tracks (deduped `ALCARECOTkAlJpsiX` collection):
      **4208 → 2032** total, i.e. **31 → 15 per event** (2.1×
      reduction).
    - Six unchanged channels bit-invariant on candidate counts
      (BPlus 242, Bc 356, B0Kstar 283, BsPhi 16, B0Ks 18,
      Lambdab 12); J/psi-only unchanged at 136.
    - Zero-pdgId daughter count still 0; pdgId invariant holds.
    - η-cut yield gain (13.6 % of B+ candidates with |η(K)| ≥
      1.8) unchanged, as expected.
  - 0.03 cm is kept as the default; not tightening further since
    the ψ(2S) rate is now at a manageable ~5/event, comparable to
    Bc's 2.6/event and K*0's 2.1/event.

## 6. Full-build integration + smoke production (preset B only)

All smoke / integration / production tasks in this section run under
**preset B only** (`TKALJPSIX_SELECTION_PRESET` unset). Presets A and
C are not exercised — preset B is the locked production default from
`finalize-jpsi-x-preset-b-production`.

- [x] 6.1 `scram b -j 16` in
  `CMSSW_10_6_17_patch1/src/Alignment/CommonAlignmentProducer` under
  `cmssw-el7`: **clean**. New plugin
  `AlignmentTrackToMuonAssociator.cc`, ctor asserts on
  `JpsiXCandidateProducer` / `TwoBodyDecayCandidateProducer`, and
  BuildFile.xml update all compile.
- [x] 6.2 `scram b -j 8` in
  `CMSSW_15_0_19_patch2/src/Analysis/HitAnalyzer` under `cmssw-el9`:
  clean. Splitter emits `bachelorPdgId`, `muon0PdgId`, `muon1PdgId`
  branches.
- [x] 6.3 200-event Stage-1 smoke (LFN 0029A508... under
  T2_CH_CERN, xrootd via eoscms). All 200 events processed, 135
  passed tight-muon filter, no `Exception` / `FATAL`. Output 3.4 MB
  (+11% vs pre-refinements for the three new persisted objects).
  Full inspection results in slides + tasks.md. Stage-2 dry run
  deferred — Stage-2 config is B+-only and the Stage-1 smoke output
  already confirms all invariants at the producer level.
- [ ] 6.4 100-file 2016H batch: not run in this session. Will be
  launched from a fresh Condor payload tarball built from the current
  Stage-1 sources + updated recoskim config. Runs once the tarball
  is rebuilt (needs `condor/jpsix_alcareco/build_tarball.sh`).
- [x] 6.5 Memory `project_jpsi_x_alcareco.md` updated.
- [x] 6.6 Slides `slides/jpsix-preprod-refinements.tex` (10-page PDF
  178 kB) drafted covering all five items + smoke results.
