# Design: TkAlJpsiX candidate quality cuts

## Context

The TkAlJpsiX AlCaReco reconstructs seven B-meson/quarkonium decay channels by
combining J/ψ→μμ candidates with bachelor tracks or intermediate resonance candidates.
First-data validation showed combinatorial background rates of 14–202 candidates/event
(see `openspec/changes/add-jpsi-x-alcareco-channels/` for the original implementation).

This design documents the reasoning behind each proposed cut value and the
architecture of the Phase 2 C++ additions.

## Reference points

| Stream         | Combinatorial suppression mechanism              | Observed yield    |
|----------------|--------------------------------------------------|-------------------|
| TkAlDstToD0Pi  | 6 MeV ΔM(D*−D0) window; no vertex fit           | ~1–2 cands/event  |
| TkAlKsToPiPi   | V0Producer: flight sig., pointing, ±70 MeV mass | ~1–2 cands/event  |
| TkAlLambdaToProtonPi | V0Producer: flight sig., pointing, ±50 MeV | ~1 cand/event   |
| TkAlJpsiX (before) | Only wide mass windows, no geometry          | 14–202 cands/event|

The D* stream achieves clean reconstruction without any vertex fit purely by exploiting
the 145 MeV D*−D0 mass difference, which is kinematically constrained to a 6 MeV window.
The V0 streams exploit the V0Producer's built-in vertex quality reconstruction.
TkAlJpsiX has no equivalent tight handle; cuts must be applied at both the intermediate
resonance level and the mother level.

## Source of truth for analysis selections

All cut values are anchored to the BtoJpsiK analysis (`btojpsik_selections.py`,
BPH-21-006). The AlCaReco is designed as a superset of the analysis sample:
each cut is loosened by 2–5× to maximise track acceptance while eliminating
the bulk of combinatorial background. MVA-based cuts (softMVA, BDT) and
full kinematic vertex fits are skipped for Phase 1+2 (see Phase 3 below).

## Cut value decisions

### J/ψ mass window: 2.7–3.4 → 2.95–3.25 GeV

The J/ψ mass is 3.097 GeV. The CMS dimuon mass resolution in Run 2 is approximately
25–35 MeV (σ) for J/ψ pT > 3 GeV, dominated by tracker momentum resolution.

A ±5σ window (±150 MeV) covers >99.99% of real J/ψ decays while being 2.3× narrower
than the current 700 MeV window. This is the single highest-leverage cut: every fake
J/ψ accepted multiplies the combinatorial by the number of bachelor tracks/resonance
candidates in the event.

Alternatives considered:
- ±3σ (3.01–3.19 GeV, 180 MeV): reduces combinatorial ~3.9× but tighter than needed
  and could clip tails at low J/ψ pT where resolution is worse
- Current 2.7–3.4 GeV: motivated by TkAlJpsiMuMu which does NOT combine J/ψ with
  anything further and so has no combinatorial explosion; inappropriate here

### K*0(892) mass window: 0.75–1.05 → 0.80–0.99 GeV

The K*0(892) has PDG mass 895.55 MeV and natural width Γ = 47.3 MeV. The reconstructed
Kπ invariant mass distribution is dominated by the natural width (Breit-Wigner), with
detector resolution (~few MeV for soft kaons) negligible in comparison.

A ±95 MeV window around 895 MeV (0.80–0.99 GeV) corresponds to ±2.0 Γ and is the
widest meaningful K*0 selection: it contains >90% of the Breit-Wigner signal while
being 1.6× tighter than the original ±3.3 Γ window. A tighter ±1.6 Γ = [0.82, 0.97]
window captures slightly less background but removes some real K*0 signal from the wide
tails — on a first validation pass ±2 Γ is the better starting point.

B0→K*0 is the most combinatorially explosive channel (~202 cands/event) because K*0
candidates (pairs of all generalTracks with right mass) multiply by J/ψ candidates.
Every factor of reduction in both windows directly reduces the B0 combinatorial.

### φ(1020) mass window: 0.99–1.06 → 0.990–1.040 GeV

The φ(1020) has PDG mass 1019.46 MeV and natural width Γ = 4.25 MeV. Unlike the K*0,
the φ is narrow enough that the reconstructed K+K- mass has resolution comparable to
the natural width (~2–5 MeV).

The window [0.990, 1.040] GeV has a physical lower floor: the K+K- production threshold
is 2mK = 987 MeV = 0.987 GeV, so 0.990 GeV is essentially the minimum kinematically
allowed K+K- invariant mass. No φ signal exists below 0.987 GeV. The upper edge at
1.040 GeV is the primary improvement over the original 1.060 GeV: the region 1.040–1.060
is >5Γ from the φ peak and contains only combinatorial background from high-momentum
kaon pairs. Tightening the upper edge from 1.060 to 1.040 is the main gain; keeping
the lower edge at 0.990 means we accept a few MeV of near-threshold combinatorial but
do not lose any φ signal.

### Bachelor pT thresholds

**B+ → J/ψ K+ (unchanged at 0.5 GeV):**
The BtoJpsiK analysis has no minimum kaon pT cut (only an upper cut < 8 GeV to
reject muon mis-ID). The existing 0.5 GeV default in `JpsiXCandidateProducer` matches
the implicit HLT bachelor-track floor for Charmonium paths (HLT_DoubleMu4_JpsiTrk
requires the bachelor track to have pT ≥ ~0.5 GeV). With the J/ψ mass window
tightened and minMotherPt = 5 GeV, any soft kaon that combines with a real J/ψ to
give B mass in [5.0, 5.5] GeV will survive and contribute alignment information.
No Phase 1 change is needed for this parameter.

**Bc → J/ψ π+ (0.5 → 0.3 GeV):**
The pion mass is 0.140 GeV, much lighter than the kaon (0.494 GeV). For Bc mesons
(mass 6.27 GeV) decaying to J/ψ π+, the pion in the B rest frame carries
p*(π) = (m²(Bc) − m²(J/ψ)) / (2 m(Bc)) ≈ 1.25 GeV/c but at low B pT the lab-frame
pion pT can be well below 0.5 GeV. The CMS tracker can reconstruct charged tracks
down to ~0.3 GeV pT in the 3.8T field; below this the helix curls too tightly for
reliable pattern recognition. Setting the floor to 0.3 GeV maximises Bc acceptance
without accepting genuinely untrackable pions.

### J/ψ pT threshold (Phase 2): > 3 GeV for all channels

The analysis requires J/ψ pT > 7 GeV (from `btojpsik_selections.py`). This cut is
almost entirely redundant with `minMotherPt = 5 GeV` + bachelor pT ≥ 0.5 GeV: if
the mother has pT > 5 GeV and the bachelor has pT ≥ 0.5 GeV, the J/ψ must have
pT > ~4.5 GeV by momentum conservation (approximately). Setting minJpsiPt = 3 GeV
provides only an early loop exit to skip pathological J/ψ candidates at very low pT
where the dimuon mass resolution is worst (σ → 50 MeV or more). This cut is applied
inside `JpsiXCandidateProducer` (not `TwoBodyDecayCandidateProducer`) because it can
reuse the existing J/ψ candidate loop with zero rebuild cost for the shared collection.

### Mother pT threshold (Phase 2): 5 GeV for B+, B0→K*0, Bs, Bc

The Charmonium HLT paths (e.g. HLT_DoubleMu4_JpsiTrk_Displaced) require the
J/ψ+track system to have pT > ~7 GeV. The reconstructed B candidate pT is dominated
by the J/ψ pT (>4 GeV by HLT) plus the bachelor/resonance; real B candidates in the
triggered sample have B pT > ~6–7 GeV. A 5 GeV threshold retains >90% of signal
while removing a significant fraction of the low-pT combinatorial (fake J/ψ at low
pT paired with soft tracks that happen to fall in the mass window).

Channels with few candidates (B0→Ks, Λb, ψ(2S)) are not given a pT threshold
since their low yields suggest we already have tight enough selection from the V0
intermediate resonance quality.

### η cuts on bachelor tracks (Phase 2): |η| < 2.4

The analysis requires muon |η| < 1.4 and kaon |η| < 1.4. For the AlCaReco the
tracker acceptance extends to |η| = 2.4; using the full tracker is the AlCaReco
purpose (alignment uses all η). No analysis motivation to restrict; the η cut is
only included to avoid spurious combinatorial from extremely forward low-quality
tracks. Applied via `maxBachelorEta` in `JpsiXCandidateProducer` (track mode) and
`maxDaughterEta` in `TwoBodyDecayCandidateProducer` (K*0, φ daughters).

### J/ψ pointing angle (Phase 2): alphaBS < 1.0

The analysis requires J/ψ alphaBS < 0.4 (from `select_dimuon_alphabs` in
`btojpsik_selections.py`). alphaBS is the angle between the J/ψ momentum vector and
the direction from the beamspot to the J/ψ vertex position. For real J/ψ from B decay,
this angle is small. For random fake J/ψ pairs, it is uniformly distributed up to π.

A cut at 1.0 radian (cos > 0.54) means the J/ψ must point at least 54% of its momentum
in the direction of the beamspot flight axis. This is 2.5× looser than the analysis
but still rejects clearly backwards-pointing or perpendicular fake J/ψ candidates.
Going to π/2 = 1.57 rad would be effectively no cut (accepts all non-backwards-pointing
J/ψ). Applied via `maxJpsiAlphaBS` in `JpsiXCandidateProducer`.

### Why no J/ψ vertex probability cut

The analysis requires J/ψ vtx_prob > 0.1. However, inspection of
`TwoBodyDecayCandidateProducer.cc` confirms it does NOT perform a kinematic
vertex fit: it uses the midpoint of daughter reference points as the vertex
position and stores no vertex quality in the VCC. Adding minJpsiVtxProb
would require a full vertex fitter dependency inside `TwoBodyDecayCandidateProducer`
(or a separate refit step), which is contrary to the "avoid complex vertex fits"
principle. This is deferred to Phase 3.

### Bachelor DCA to J/ψ vertex (Phase 2, track mode only)

For real B+ → J/ψ K+, the bachelor kaon originates from the B decay vertex (secondary
vertex), which is spatially close to the J/ψ vertex (both daughters of the B). The
DCA of the bachelor track to the J/ψ vertex position (approximate B decay vertex)
is expected to be small (<0.5 mm typically for prompt B at CMS vertex resolution).

For fake combinations (random tracks from the primary vertex or other secondaries),
the bachelor track's DCA to the J/ψ vertex can be large. A loose cut of
`maxBachelorIPToJpsiVertex < 10 mm` (initial value; to be tuned from data) provides
additional suppression without requiring a full vertex fit. 10 mm is large enough to
accept any bachelor track from a genuine B decay vertex (typical B flight ~0.1–0.5 mm
displaced from IP, so the J/ψ vertex midpoint and B vertex are within ~1 mm of each
other). It only rejects tracks from long-lived secondaries (K_L, etc.) or very
displaced pile-up vertices.

Implementation: the J/ψ vertex position is approximated as the midpoint of the J/ψ
daughter tracks' point of closest approach, already computed in `JpsiXCandidateProducer`.
The DCA is the 3D distance between the bachelor track's helix and this point.

## Architecture of Phase 2 C++ additions

### `JpsiXCandidateProducer.cc` — new optional parameters

All parameters use `cfg.existsAs<T>()` with defaults that reproduce current behaviour:

```cpp
minJpsiPt_(cfg.existsAs<double>("minJpsiPt")
               ? cfg.getParameter<double>("minJpsiPt") : 0.0)
minMotherPt_(cfg.existsAs<double>("minMotherPt")
                 ? cfg.getParameter<double>("minMotherPt") : 0.0)
maxBachelorEta_(cfg.existsAs<double>("maxBachelorEta")
                    ? cfg.getParameter<double>("maxBachelorEta")
                    : std::numeric_limits<double>::max())
maxJpsiAlphaBS_(cfg.existsAs<double>("maxJpsiAlphaBS")
                    ? cfg.getParameter<double>("maxJpsiAlphaBS")
                    : std::numeric_limits<double>::max())
maxBachelorIPToJpsiVertex_(cfg.existsAs<double>("maxBachelorIPToJpsiVertex")
                               ? cfg.getParameter<double>("maxBachelorIPToJpsiVertex")
                               : std::numeric_limits<double>::max())
```

`maxJpsiAlphaBS_` requires an `offlineBeamSpot` InputTag consumed when the parameter
is set; the AlphaBS calculation uses `reco::VertexCompositeCandidate::vertex()` (the
midpoint stored by TwoBodyDecayCandidateProducer) and the beamspot position.

### `TwoBodyDecayCandidateProducer.cc` — new optional parameter

```cpp
maxDaughterEta_(cfg.existsAs<double>("maxDaughterEta")
                    ? cfg.getParameter<double>("maxDaughterEta")
                    : std::numeric_limits<double>::max())
```

Applied as a pre-filter on both daughters before the mass-window combinatorics loop.

## Risks / Trade-offs

- **Signal efficiency loss**: Tighter J/ψ window could clip tails for low-pT J/ψ where
  resolution is ~50 MeV. Mitigation: the chosen ±5σ (±150 MeV) window is generous.
- **B0→K*0 residual combinatorial**: Even after all Phase 1+2 cuts, B0→K*0 is expected
  to yield ~10–20 cands/event, still too high for a clean peak. A proper kinematic
  vertex fit (Phase 3) is ultimately needed.
- **No retuning needed for rare channels**: B0→Ks, Λb, ψ(2S) have <25 candidates
  total in the test sample; their configuration is unchanged.

## Open Questions

- Should the J/ψ window be further tightened to ±3σ (3.01–3.19 GeV) after Phase 1
  validation? Proposed ±5σ (2.95–3.25 GeV) is conservative; revisit after seeing
  new mass plots.
- What is the optimal `maxBachelorIPToJpsiVertex` threshold? Initial value 5 mm to
  be tuned from data once Phase 2 is deployed.
- Should a Phase 3 (proper kinematic vertex fit via `KinematicParticleVertexFitter`)
  be planned for B0→K*0 specifically?
