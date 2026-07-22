# Change: AlCaReco -> CVH refit -> NanoAOD stream (one job, all channels)

## Why

The full 2016 postVFP AlCaReco production
(`/ceph/submit/data/user/p/pmlugato/mz/alcareco/full_2016postvfp/`)
gives small, alignment-oriented EDM files: cloned `reco::Track`s with
hits + dE/dx, `offlinePrimaryVertices`, per-channel
`reco::VertexCompositeCandidate` resonances, and (J/psi streams) a
`reco::Muon` collection with a track->muon association. Downstream
(the `btojpsik` histmaker, `narf`, `uproot`) wants flat NanoAOD.

This stream does two things in **one `cmsRun` job**:

1. Runs the persisted candidates' tracks through the **CVH refit**
   (`ResidualGlobalCorrectionMaker*`, Geant4e) in *application* mode —
   **no gradients, no runtree, no aggregate/solve**. This is not the
   calibration; it is just re-fitting the candidate tracks and getting
   corrected track/candidate quantities out.
2. Writes a proper **NanoAOD** carrying, per candidate, **both** the
   raw AlCaReco quantities (from the persisted VCC) **and** the
   CVH-corrected quantities (from the maker's EDM ValueMaps), via
   `PhysicsTools/NanoAOD` table producers and `NanoAODOutputModule`.

Two enabling facts make this mostly wiring, not new fit code:

- The two-track maker already emits per-candidate EDM ValueMaps
  (`corMass/corMassErr/corPt/corEta/corPhi` + per-leg refit
  kinematics) behind a `produceValueMaps` flag, and David already
  wired a `SimpleCandidateFlatTableProducer` "Dimuon" table
  (`CandVars + cvh* externalVariables`) on MINIAOD (`cd9950c`). We
  generalize that pattern to AlCaReco input and all channels.
- The **single-G4-master limitation is fixed** (`e33c8d`: shared
  Geant4 master as an ESProducer on `CvhMasterRecord`, consumed by all
  makers in a job). So the dimuon + bachelor refits **coexist in one
  job**. The old split-into-two-jobs + offline join
  (`join_cvh_bplus_jpsik.py`) is no longer required.

## What Changes

### One-job refit + nano pipeline (replaces the two-job + offline join)

```
PoolSource(AlCaReco)
  -> JpsiKCandidateSplitter            # trivial in-memory adapter, no fit
  -> globalCor<chan>      (two-track dimuon refit, produceValueMaps=True)
  -> globalCor<chan>Kaon  (single-track bachelor refit, produceValueMaps=True)
     [both share the CvhMasterRecord ES product]
  -> per-channel Simple*FlatTableProducer  (raw VCC vars + corrected ValueMaps as externalVariables)
  -> Track / PV / Muon tables (+ dE/dx, originalIndex externalVariables)
  -> NanoAODOutputModule
```

- **No offline recombine.** Both makers' per-candidate ValueMaps are
  present in the same event, keyed to the candidate via the splitter's
  `bCandIdx`. The corrected `m(mu mu K)` is formed in-job (small
  combiner producer or a derived table column), not by
  `join_cvh_bplus_jpsik.py`.
- **One job, not two.** `reconstruction_step` schedules the splitter
  and *both* makers. Remove the either/or `--maker` branch and the
  stale "the two makers cannot coexist" comments in
  `test/runCvhBplusJpsiK.py`.
- **Refit config for v1:** `useIdealGeometry = False`, **no
  correction file applied** (identity corrections). This still yields
  a distinct refit output and exercises the full
  AlCaReco -> refit -> NanoAOD infrastructure; the real 2016 postVFP
  corrections drop in later (they were only just produced).

### Tables

| Content | NanoAOD table | Source |
|---|---|---|
| candidate raw + corrected | one per channel | VCC vars (`mass`,`pt`,daughters) + maker ValueMaps (`corMass`,`corPt`,...) as `externalVariables` |
| tracks | `Track` | `SimpleTrackFlatTableProducer` + dE/dx + `originalIndex` externalVariables |
| primary vertices | `PV` | `SimpleVertexFlatTableProducer` |
| muons | `Muon` | `reco::Muon` typedef table (Q3 = option A) |
| HLT decisions | `HLT_*` | stock trigger-flag machinery |

Branch names follow the histmaker's raw contract (`bkmm_kaon_*`,
`mm_mu*`, `Muon_*`) so `btojpsik.py` runs over this nano via its
existing `get_bkmm_alcareco_selections` path (see design D0), with the
corrected columns added alongside.

### Generalization to all channels

All persisted channels get the refit->table path, not just B+:
`TkAlJpsiX` (BuJpsiK, B0Kstar, B0Ks, BsPhi, Lambdab, Bc, Psi2S,
JpsiOnly), `TkAlJpsiMuMu`, `TkAlDstToD0Pi`, `TkAlKsToPiPi`,
`TkAlLambdaToProtonPi`. Multi-bachelor channels (K*0->Kpi, phi->KK,
psi2S->J/psi pipi, ...) need the bachelor refit applied per bachelor
track; the splitter already tags species by daughter `pdgId`.

### Documentation deliverable — expand the NanoAOD-production deck

The two background slides already built
(`~/public_html/slides/260722_central_production.pdf`: the central
RAW->AOD->MINIAOD->NANOAOD chain with the parallel ALCA branch, and
how `SimpleFlatTableProducer` turns collections into a flat `Events`
tree) grow into a full narrative deck covering this change:

1. *(existing)* Central tier chain; where AlCaReco sits.
2. *(existing)* How NanoAOD is built — table producers, `variables` /
   `externalVariables`, the flat `Events` layout.
3. **How ours differs.** AlCaReco is not MINIAOD: no PF, jets, MET,
   e/gamma, no PAT. Inventory of what the AlCaReco *does* hold and,
   consequently, which tables can and cannot be filled. Why
   `mz_dilepton` / `mw_*` cannot run on this nano (they need the
   custom CVH muon production).
4. **What we needed to do.** The refit in the loop (application mode,
   not calibration); collapsing two jobs into one after the shared-G4
   -master fix; generic nested-VCC decomposition replacing the
   splitter; the missing `reco::Muon` table producer; cross-link index
   columns as the flat-tree substitute for EDM refs.
5. **What we did.** The one-job pipeline diagram, the table inventory
   (raw + corrected), and the branch-name contract that lets the
   existing `btojpsik` histmaker read this nano unchanged.
6. **Validation.** Cross-release read smoke; raw-vs-corrected closure
   against the retired offline join; per-channel candidate counts;
   resource scatter vs. the two-job baseline.
7. **Status and next steps.** Identity corrections now, real 2016
   postVFP corrections later; batch production over the full
   AlCaReco set.

Built with the `mit-slides` skill, same deck source, so it supersedes
the 2-slide version at the same path.

## Impact

- Affected code (`CMSSW_15_0_19_patch2/src`):
  - `Analysis/HitAnalyzer/test/runCvhBplusJpsiK.py` — collapse to one
    job (both makers), drop the offline-join branch, delete stale
    coexistence comments; add the nano table producers + output.
  - `Analysis/HitAnalyzer/plugins/ResidualGlobalCorrectionMakerTwoTrackG4e.cc`
    and `...MakerG4e.cc` (single-track) — `produceValueMaps` path used
    for all channels (already exists for two-track; confirm/extend for
    single-track).
  - small in-job combiner producer for `m(mu mu K)` from the two
    makers' ValueMaps (replaces `join_cvh_bplus_jpsik.py`).
  - `PhysicsTools/NanoAOD/plugins/SimpleFlatTableProducerPlugins.cc` —
    `reco::Muon` typedef (Q3=A); `reco::VertexCompositeCandidate`
    typedef if daughter accessors are needed (see design finding).
  - per-channel nano `_cff.py` + top-level config.
- Deprecated by this change: `scripts/btojpsik/join_cvh_bplus_jpsik.py`
  for the production path (kept only if still used for diagnostics).
- No changes in `narf`, `rabbit`, `wremnants/` (histmaker already has
  the AlCaReco path; new corrected columns are additive).

## Open Questions

- **Q1 [RESOLVED].** Cross-release read works — 500-event smoke in
  `CMSSW_15_0_19_patch2` on a `10_6` `TkAlJpsiX` file passed
  (Track+dE/dx+PV+candidate populated; `BuJpsiK_mass` 5.05-5.39 GeV).
- **Q2 [DECIDED].** Bespoke schema named to the histmaker's raw
  contract; corrected columns added alongside (design D0).
- **Q3 [DECIDED].** Muon table = option A (`reco::Muon` typedef).
  Confirmed the CVH muon nano branches `mz_dilepton`/`mw_*` need are
  a *separate* custom-production path, not derivable here; option A is
  the most capable and costs one line.
- **Q4 [DECIDED].** All channels from v1.
- **Q5 [DECIDED].** Hits/clusters dropped from the flat ntuple.
- **Q6 [DECIDED].** Emit **all three** cross-link index columns:
  daughter->Track (`<chan>_kaonTrackIdx`, `_mu1TrackIdx`,
  `_mu2TrackIdx`), Track->Muon (`Track_muonIdx`), Track->PV
  (`Track_pvIdx`). A flat tree has no pointers, so these integer
  columns are the only thing that makes the tables joinable — e.g.
  reading the dE/dx of a selected candidate's kaon requires
  `Track_dedxHarmonic2[<chan>_kaonTrackIdx[i]]`. Cost is one
  `int16`/`int` per row (tens of bytes/event), negligible against the
  float kinematics, and it is standard NanoAOD practice
  (`Muon_jetIdx`, `Muon_svIdx`, `Muon_genPartIdx`).
- **Q7 [DECIDED].** Persist L1 decision bits, and `DcsStatus` as
  event-level branches — confirmed lightweight: the object is two
  floats (`magnetCurrent`, `magnetTemperature`) plus a 25-bit `ready`
  mask, ~4-12 bytes/event.
- **Q8 [DECIDED].** **Delete `JpsiKCandidateSplitter`**; teach the
  makers to descend the nested VCC. The VCC format is the common
  interface across all channels, whereas splitters would need a
  per-topology variant (0/1/2 bachelors, nested vs. flat). The maker
  already does generic daughter traversal with a `dynamic_cast` that
  is itself the leaf-vs-composite discriminator, so this is ~10 lines
  in one block, not a refactor. See design D2a for the generic
  decomposition rule and the per-channel table.

## Non-goals

- The CVH calibration itself (gradients, aggregate, solve) — separate,
  existing pipeline.
- Applying real 2016 postVFP corrections in v1 (identity for now).
- Batch production at scale (follow-up once the per-file config is
  validated on a smoke).
- A unified N-track maker (no such maker exists; out of scope).
