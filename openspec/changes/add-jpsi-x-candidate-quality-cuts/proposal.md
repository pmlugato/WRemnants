# Change: Add kinematic quality cuts to TkAlJpsiX candidate reconstruction

## Why

First-data validation of TkAlJpsiX (678 events, Run2016H Charmonium) shows
combinatorial background dominating every channel with a wide-window intermediate
resonance: B+→J/ψK+ ~14/event, B0→J/ψK*0 ~202/event, Bs→J/ψφ ~18.5/event,
Bc→J/ψπ+ ~22.8/event. No signal peak is visible. V0-mode channels (B0→Ks, Λb,
ψ(2S)) are clean because the V0Producer already applies vertex quality.

Root causes:
1. J/ψ mass window 2.7–3.4 GeV (700 MeV) accepts many fake J/ψ from random OS
   track pairs; every fake multiplies the combinatorial downstream.
2. No pT, η, or geometry-based cuts applied at the intermediate or mother level.
   The downstream BtoJpsiK analysis applies a rich set of selections; the AlCaReco
   should mirror those, loosened by a well-defined factor.

This change derives the AlCaReco selections from `scripts/histmakers/btojpsik.py`
(anchored to BPH-21-006 selections), loosens each by ~2–5× to maximise track
acceptance, skips MVA-based or expensive vertex-fit requirements for Phase 1+2,
and migrates the resulting selection set to all 7 channels.

## Analysis selections (BtoJpsiK, source of truth)

From `btojpsik_selections.py` / `get_bkmm_selections()`:

| Variable | Analysis cut | Type |
|----------|-------------|------|
| dimuon charge | OS | kinematic |
| muon \|η\| | < 1.4 | kinematic |
| muon pT | > 4 GeV (both) | kinematic |
| muon softMVA | > 0.45 (both) | MVA — skip AlCaReco |
| J/ψ pT | > 7 GeV | kinematic |
| J/ψ alphaBS | < 0.4 | geometry (beamspot) |
| J/ψ vtx_prob | > 0.1 | vertex fit |
| J/ψ sl3d | > 4 | vertex fit + flight |
| B vtx_prob | > 0.3 (jpsimc) / > 0.1 (nomc) | vertex fit |
| B mass | \|m − 5.3\| < 0.1 GeV → [5.2, 5.4] | kinematic |
| kaon \|η\| | < 1.4 | kinematic |
| kaon pT upper | < 8 GeV | kinematic |
| BDT | > 0.10 | MVA — skip AlCaReco |

## Proposed AlCaReco selections (B+→J/ψK+) and mapping to other channels

### Design principles
- Start loose: relax each cut by ≥ 2× vs analysis; tighten in subsequent proposals
  once data yields are known
- Skip MVA (softMVA, BDT) and expensive B-vertex fits (Phase 1+2); flag as Phase 3
- Mirror analysis structure: per-particle pT+η, intermediate resonance mass,
  mother mass, mother pT; omit displacement/vertex fit for now
- All cuts go in Python config (Phase 1) or new C++ parameters (Phase 2)

### Cut table by channel

#### B+→J/ψK+ (primary reference)

| Variable | Analysis | AlCaReco Phase 1+2 | Rationale |
|----------|----------|---------------------|-----------|
| Muon \|η\| | < 1.4 | < 2.4 | Full tracker acceptance; alignment needs all η |
| Muon pT | > 4 GeV | > 3.5 GeV | Slightly below analysis; HLT floor already ~4 GeV so near-free |
| J/ψ mass | (from vtx fit) | **2.95–3.25 GeV** | ±5σ at 30 MeV resolution; 2.3× tighter than current |
| J/ψ pT | > 7 GeV | **> 3 GeV** | Phase 2; largely redundant with mother pT > 5 GeV; provides early loop exit |
| J/ψ alphaBS | < 0.4 | **< 1.0** | Phase 2; 2.5× looser; cos(1.0) = 0.54; still rejects backwards-pointing J/ψ |
| J/ψ vtx_prob | > 0.1 | skip Phase 1+2 | TwoBodyDecayCandidateProducer stores no vtx quality; defer to Phase 3 |
| J/ψ sl3d | > 4 | skip Phase 1+2 | Needs PV; defer to Phase 3 |
| B vtx_prob | > 0.3 | skip Phase 1+2 | Needs full B vertex fit; defer to Phase 3 |
| B mass | [5.2, 5.4] | **[5.0, 5.5]** | ±250 MeV, 2.5× looser than current analysis window |
| Kaon \|η\| | < 1.4 | **< 2.4** | Full tracker acceptance |
| Kaon pT min | (HLT floor ~0.5) | **> 0.5 GeV** | HLT bachelor-track implicit floor; no minimum in analysis |
| Kaon pT upper | < 8 GeV | skip | No motivation for upper cut in AlCaReco |
| Mother pT | implied by J/ψ pT | **> 5 GeV** | Primary combinatorial suppressor; below HLT trigger floor (~7 GeV) |

#### B0→J/ψK*0 (analogous, replace kaon→Kπ pair)

| Variable | AlCaReco | Derivation |
|----------|----------|-----------|
| J/ψ | same as B+ | shared Stage-1 module |
| K*0 mass | **[0.80, 0.99] GeV** | ±95 MeV = ±2.0 Γ; still tighter than original ±3.3Γ [0.75, 1.05] |
| K pT from K*0 | **> 0.5 GeV** | softer than B+ kaon; starts loosely |
| π pT from K*0 | **> 0.3 GeV** | pion is light; V0-like floor |
| B0 mass | **[5.0, 5.5] GeV** | same as B+ |
| B0 pT | **> 5 GeV** | same as B+ |
| B0 kaon \|η\| | **< 2.4** | same |
| B0 pion \|η\| | **< 2.4** | same |

Note: K and π η cuts applied inside `TwoBodyDecayCandidateProducer` for K*0
(new `maxDaughterEta` parameter) OR inside `JpsiXCandidateProducer` for VCC mode.

#### Bs→J/ψφ (analogous, replace kaon→K+K- pair)

| Variable | AlCaReco | Derivation |
|----------|----------|-----------|
| φ mass | **[0.990, 1.040] GeV** | lower edge at K+K- threshold (2mK = 0.987 GeV); upper tightened from 1.060 |
| K pT from φ | **> 0.5 GeV** | φ daughters; looser than B+ kaon |
| Bs mass | **[5.2, 5.6] GeV** | current window; already looser than B+ analysis |
| Bs pT | **> 5 GeV** | same as B+ |

#### B0→J/ψKs (V0-mode)

V0Producer already provides flight significance, pointing angle, ±70 MeV mass cut.
Only add:

| Variable | AlCaReco | Derivation |
|----------|----------|-----------|
| B0 mass | **[5.0, 5.5] GeV** | same as B+→K+ |
| B0 pT | **> 5 GeV** | same |

#### Λb→J/ψΛ (V0-mode)

V0Producer already provides flight significance, pointing angle, ±50 MeV mass cut.

| Variable | AlCaReco | Derivation |
|----------|----------|-----------|
| Λb mass | **[5.3, 6.0] GeV** | current; generous given kinematics |
| Λb pT | **> 5 GeV** | same |

#### ψ(2S)→J/ψKs (V0-mode)

| Variable | AlCaReco | Derivation |
|----------|----------|-----------|
| ψ(2S) mass | **[3.5, 3.9] GeV** | current; already reasonable |
| ψ(2S) pT | **> 3 GeV** | lower than B channels; ψ(2S) mass is 3.7 GeV |

#### Bc→J/ψπ+ (track mode)

| Variable | AlCaReco | Derivation |
|----------|----------|-----------|
| π pT | **> 0.3 GeV** | tracker reconstruction floor in 3.8T field; pion is lighter than kaon |
| π \|η\| | **< 2.4** | same |
| Bc mass | **[5.9, 6.6] GeV** | current; already reasonable |
| Bc pT | **> 5 GeV** | same |

## What Changes

### Phase 1 — Python configuration only (no recompile)

1. J/ψ mass window: 2.7–3.4 → **2.95–3.25 GeV** in `ALCARECOTkAlJpsiXJpsiCandidates`
2. K*0 mass window: 0.75–1.05 → **0.80–0.99 GeV** in `ALCARECOTkAlJpsiXKstarCandidates`
3. φ mass window: 0.99–1.06 → **0.990–1.040 GeV** in `ALCARECOTkAlJpsiXPhiCandidates` (upper edge only tightened)
4. B+ bachelor pT: unchanged at **0.5 GeV** in `ALCARECOTkAlJpsiXBPlusCandidates` (no change from current default)
5. Bc bachelor pT: 0.5 → **0.3 GeV** in `ALCARECOTkAlJpsiXBcCandidates`
6. B mass windows reviewed; current [5.0, 5.5] / [5.2, 5.6] etc. already looser than
   analysis; **no change needed** for B+, B0, Bs, Λb, ψ(2S), Bc

### Phase 2 — C++ extensions (new parameters, backward-compatible defaults)

All new parameters use `cfg.existsAs<T>()` with defaults that reproduce current behaviour.

**`JpsiXCandidateProducer.cc`:**
- `minJpsiPt` (default 0): J/ψ pT cut before combining → configure **3.0 GeV** for all channels (early loop exit; largely redundant with minMotherPt > 5 GeV)
- `minMotherPt` (default 0): mother pT cut → configure **5.0 GeV** for B+, B0→K*0, Bs, Bc; **3.0 GeV** for ψ(2S); **0** for Λb, B0→Ks
- `maxBachelorEta` (default +inf, track mode): |η| cut on bachelor → configure **2.4** for B+, Bc
- `maxJpsiAlphaBS` (default +inf): pointing angle cut on J/ψ; requires `beamspot` InputTag → configure **1.0** for all (2.5× looser than analysis 0.4; cos(1.0) = 0.54)
- `maxBachelorIPToJpsiVertex` (default +inf): 3D DCA of bachelor to J/ψ vertex midpoint → configure **10 mm** for B+, Bc (initial loose value; to be tuned from data)
- `minJpsiVtxProb`: **not implemented** — TwoBodyDecayCandidateProducer stores no vertex quality; defer to Phase 3

**`TwoBodyDecayCandidateProducer.cc`:**
- `maxDaughterEta` (default +inf): |η| cut on each daughter track → configure **2.4** for K*0 (K and π daughters) and φ (K daughters)

## Impact

- Affected specs: `alcareco-jpsi-x`
- Affected code:
  - `ALCARECOTkAlJpsiX_cff.py` (Phase 1 + Phase 2 configuration)
  - `JpsiXCandidateProducer.cc` (Phase 2)
  - `TwoBodyDecayCandidateProducer.cc` (Phase 2, maxDaughterEta)
- Not breaking: all new C++ parameters have backward-compatible defaults

## Expected improvement

| Channel | Current /event | After Phase 1 | After Phase 1+2 (estimate) |
|---------|---------------|---------------|---------------------------|
| B+ | ~14 | ~4–6 | ~1–2 |
| B0→K*0 | ~202 | ~50–70 | ~15–25 |
| Bs→φ | ~18.5 | ~5–7 | ~1–3 |
| Bc | ~23 | ~6–8 | ~2–4 |
| B0→Ks | ~0.03 | unchanged | unchanged |
| Λb | ~0.004 | unchanged | unchanged |
| ψ(2S) | ~0.03 | unchanged | unchanged |

Phase 1 alone is expected to make the B+ and Bs→φ peaks visible; mass distributions
should be non-flat from the J/ψ mass window alone. B0→K*0 will require Phase 2
(minMotherPt, alphaBS) for a clean distribution.
