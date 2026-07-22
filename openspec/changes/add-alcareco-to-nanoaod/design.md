# Design — AlCaReco -> CVH refit -> NanoAOD

## Context: what the AlCaReco holds

Verified with `edmDumpEventContent` on the 2016 postVFP files. Common
to every stream (`TkAlJpsiX`, `TkAlJpsiMuMu`, `TkAlDstToD0Pi`,
`TkAlKsToPiPi`, `TkAlLambdaToProtonPi`):

- `vector<reco::Track>` `ALCARECOTkAl<S>` (+ `TrackExtra`, hits, Si
  clusters) — the cloned alignment tracks.
- `ValueMap<float>` `...DeDx{,All,Pixel}Harmonic2`; `ValueMap<uint>`
  `ALCARECOTkAl<S>:originalIndex`.
- `vector<reco::Vertex>` `offlinePrimaryVertices` (+ track->PV
  `Association`).
- `vector<reco::VertexCompositeCandidate>` — the resonance(s); 8 for
  `TkAlJpsiX`.
- HLT + RECO `TriggerResults`, `L1GlobalTriggerReadoutRecord`,
  `DcsStatus`.
- J/psi streams: `vector<reco::Muon>` + track->muon `Association`.

## D0 — Schema reconciliation (why the histmaker runs over both)

`scripts/histmakers/btojpsik.py` already carries two input paths:

- **analysis path** (`get_bkmm_selections`) — reads the Kalman-fit
  BParking branches `bkmm_jpsimc_*`, `mm_kin_*`.
- **AlCaReco path** (`get_bkmm_alcareco_selections`) — reads **only
  raw track-level** branches: `bkmm_kaon_pt/eta/phi/charge`,
  `bkmm_kaon_mu{1,2}_doca`, `bkmm_mm_index`, `mm_mu{1,2}_pt/eta/phi`,
  `Muon_pt/eta/charge`, and recomputes masses in
  `define_raw_kinematics` (`btojpsik_selections.py:523`).

So the AlCaReco NanoAOD emits a **bespoke schema whose names match the
raw contract**. The same histmaker then runs over both inputs with no
new path. The CVH-corrected columns are **additive** — new names
alongside the raw ones, so nothing existing breaks. This is the
reconciliation surface; it is not a naming coincidence to be
maintained by hand but the deliberate contract of the AlCaReco path.

Note the AlCaReco candidates carry **no vertex fit** under preset B —
the smoke shows `vertexChi2 == 0`. The raw path is built for exactly
that (raw kinematics + DOCA, never a fit), so this is consistent.

## D1 — CVH refit in application mode

The refit is `ResidualGlobalCorrectionMaker*` (Geant4e). For this
stream it runs with:

- `fillGrads = False`, `fillRunTree = False` — **no calibration**.
  No gradient sidecars, no aggregate, no solve.
- `useIdealGeometry = False`.
- **No correction file applied** for v1 (identity). This still
  produces a distinct refit output and exercises the whole
  AlCaReco -> refit -> NanoAOD path. The real 2016 postVFP
  corrections drop in later by pointing at a
  `correctionResults_*.root` (the maker already reads `parmtree` from
  it in `ResidualGlobalCorrectionMakerBase.cc`).

## D2 — One job, no split-then-recombine

**The constraint that forced two jobs is gone.** `CvhMasterThread.h`:
the G4 master is "shared across all CVH residual makers in a job as an
EventSetup product on `CvhMasterRecord` (produced by
`CvhMasterESProducer`)", and each stream owns its own
`Geant4ePropagator` clone plus a `CvhWorker` for per-thread G4 setup
(`ResidualGlobalCorrectionMakerBase.h:626`). Landed in `e33c8d440ac`
("CVH: shared Geant4 master via ESProducer"). `runCvhBplusJpsiK.py:359`
already wires it.

The **stale** state to fix: `runCvhBplusJpsiK.py` lines 20 and 54 still
say "the two makers cannot co-exist ... single G4 master per process",
and the schedule is an either/or — run the job once per maker, write
two sidecar TFiles, then inner-join offline on
`(run, lumi, event, bCandIdx)` via
`scripts/btojpsik/join_cvh_bplus_jpsik.py`.

New shape — one `cmsRun`:

```
splitter -> globalCor<chan> (two-track, produceValueMaps=True)
         -> globalCor<chan>Kaon (single-track, produceValueMaps=True)
         -> in-job combiner (corrected m(mu mu K))
         -> tables -> NanoAODOutputModule
```

Both makers see the same event, so the corrected quantities are paired
per candidate by `bCandIdx` in memory. `join_cvh_bplus_jpsik.py` is no
longer on the production path. Its `joinFlags` orphan accounting
(dimuon-only / kaon-only) should survive as a **column**, not a
script — a leg that failed to refit must stay visible, not silently
drop.

### D2a — Delete the splitter; makers descend the nested VCC (Q8)

**Decision: retire `JpsiKCandidateSplitter`.** The VCC format is the
common interface across all channels; a splitter approach would need a
per-topology variant (0 / 1 / 2 bachelors, nested vs. flat), which does
not scale to the 12 persisted channels.

This is cheap because the maker **already does generic daughter
traversal** (`ResidualGlobalCorrectionMakerTwoTrackG4e.cc:1142`):

```cpp
if (cand.numberOfDaughters() < 2) continue;
const auto* d0 = dynamic_cast<const reco::RecoChargedCandidate*>(cand.daughter(0));
const auto* d1 = dynamic_cast<const reco::RecoChargedCandidate*>(cand.daughter(1));
if (!d0 || !d1 || d0->track().isNull() || d1->track().isNull()) continue;
trackPairs.push_back({{&*d0->track(), &*d1->track()}});
```

The `dynamic_cast` is *already* the leaf-vs-composite discriminator: a
nested J/psi VCC fails the cast and is silently skipped today. Adding a
recursion branch when the cast fails is ~10 lines in this one block.

**Generic decomposition rule** (falls out of the VCC nesting):

- composite daughter (VCC) -> a 2-track joint-fit subsystem -> two-track maker
- leaf daughter (RCC) -> a single track -> single-track maker

Covers every channel with the two existing makers:

| Channel | Decomposition |
|---|---|
| BuJpsiK, Bc | J/psi(2) + bachelor(1) |
| B0Kstar, BsPhi, B0Ks, Lambdab, Psi2S | J/psi(2) + X(2) |
| JpsiOnly, JpsiMuMu, KsToPiPi, LambdaToProtonPi | (2) |
| DstToD0Pi | D0(2) + pi(1) |

Consequences: no new maker, no per-channel code, and the `bCandIdx`
bookkeeping disappears — the makers key their ValueMaps directly to the
original candidate collection, so the nano table joins for free.

**Coordination caveat:** these are David's actively-developed files
(several recent commits). The change is small and localized but should
go upstream rather than be forked.

## D3 — Table construction

The two-track maker already emits per-candidate EDM ValueMaps behind
`produceValueMaps` (`cd9950c`): `corMass`, `corMassErr`, `corPt`,
`corEta`, `corPhi`, per-leg refit `pt/eta/phi`, `edmval`. These attach
to a candidate table as `externalVariables` — the *same* mechanism the
smoke proved for the dE/dx maps. David's MINIAOD "Dimuon" table
(`SimpleCandidateFlatTableProducer`, `CandVars + cvh*
externalVariables`) is the template.

The single-track maker (`ResidualGlobalCorrectionMakerG4e`) needs no
flag at all — it **unconditionally** produces `corPt`, `corEta`,
`corPhi`, `corCharge`, `corDxy`, `corDz`, `edmval`, `nValidHits`,
`nValidPixelHits`, plus `globalIdxs` / `jacRef` / `momCov` vector
maps (`plugins/ResidualGlobalCorrectionMakerG4e.cc:273-280`).

So **both** makers already emit the corrected quantities as EDM
ValueMaps. Wiring the refit into NanoAOD requires no new fit code and
no new emission code: set `produceValueMaps = True` on the two-track
maker, take the single-track maps as-is, and attach both as
`externalVariables`. The `momCov` map is the analogue of
`Muon_cvhMomCov` and can be carried as a flattened vector column.

## D4 — Muon table (Q3 = option A)

Add `typedef SimpleFlatTableProducer<reco::Muon>
SimpleMuonFlatTableProducer;` + `DEFINE_FWK_MODULE` beside the
existing typedefs and rebuild `PhysicsTools/NanoAOD`. Gives the full
`reco::Muon` string surface (`isGlobalMuon`, `isTrackerMuon`,
`numberOfMatches`).

Survey result justifying this: the `Muon_cvh*` / `Muon_standalone*` /
`Muon_nTrackerLayers` / `Muon_genPartIdx` / `Muon_pfRelIso04_all`
branches that `mz_dilepton` and `mw_*` are built around come from the
**WRemnants custom CVH production** (`wremnants/production/
muon_calibration.{py,hpp}`), not from AlCaReco. Those histmakers
cannot run on this nano regardless of muon-table choice. Option A is
chosen because it is the most capable `reco::Muon` option for one
line, not because it unlocks those histmakers.

## Findings from the smoke (all verified live)

- **Cross-release read works.** 500 events of a `10_6`-produced
  `TkAlJpsiX` file processed in `CMSSW_15_0_19_patch2`.
  `Track_pt`/`Track_dedxHarmonic2`/`Track_originalIndex` populated,
  `BuJpsiK_mass` 5.05-5.39 GeV, `PV` populated.
- **Daughter accessor gap.** `SimpleCandidateFlatTableProducer` (a
  `reco::Candidate` view) throws on `daughter(0).pt()` —
  `method "daughter" returned void`. Candidate-level vars are fine.
  The raw contract needs per-daughter kinematics, so resolve with a
  concrete `SimpleFlatTableProducer<reco::VertexCompositeCandidate>`
  typedef, or emit daughter columns from the splitter/combiner (which
  already walks daughters and emits `bachelorPdgId` / `muon0PdgId` /
  `muon1PdgId` vectors). The latter also produces
  `bkmm_kaon_mu{1,2}_doca` and `bkmm_mm_index`, so it is the likely
  path.

## Operational note — `plimit` on soft bachelors

The one-job validation run showed the single-track maker losing 10 of 15
bachelor tracks to `Geant4e fail[plimit]` at the driver default
`plimit = 1` GeV: B+ bachelor kaons are soft (p ~ 0.3-0.9 GeV), so the
propagation floor rejects most of them. The two-track (dimuon) maker was
unaffected (15/15, zero propagation failures).

This is pre-existing, not a consequence of the one-job collapse, and the
calibration production already runs `plimit=0.05` for the kaon channel.
The nano production must set a lowered `plimit` or it will silently
publish candidates whose bachelor leg was never refit — which is exactly
what the orphan/status column (D2) has to make visible.

## Reading `edmval` correctly

The `edmval` the makers store is the **full-state** EDM: it includes
every per-hit scattering parameter, and those are re-zeroed at each
re-linearization. Its median is ~1e3, so a value of order 1e2 is
unremarkable and says nothing about convergence.

The actual convergence criterion is **EDM < 1e-5 on the reference-state
block**, iteration cap 10. Failing to reach it is not equivalent to a
failed fit. Measured rates: dimuon fits hit the cap without meeting the
criterion 1.3 % of the time (essentially all still genuinely improving,
chi2 still falling); single-track 15.8 %, of which about half converge
by iteration 20 and 8.2 % never do. The never-converging class are limit
cycles -- the solver keeps predicting a ~0.07 chi2 gain while 30 further
iterations move the actual chi2 by exactly zero.

**Consequence for the NanoAOD:** exposing the full-state `edmval` alone
invites exactly the misreading above. The reference-block EDM is now
recordable per iteration, and that is the quantity the tables should
carry for convergence auditing.

## Risks

- **R1.** Multi-bachelor channels (K*0->Kpi, phi->KK, psi2S->J/psi
  pipi, Lambda->p pi) need the single-track refit applied per bachelor
  track. The splitter's one-bachelor assumption must be generalized.
- **R2.** Running both makers in one job doubles per-event Geant4e
  work in a single process; watch peak RSS and per-file runtime
  against the two-job baseline before batch launch.
- **R3.** Identity corrections in v1 mean the corrected columns are
  *not* physics-final. They must not be fed to a calibration fit until
  real corrections are wired.

## Validation plan

1. `scram b` from inside the touched packages only
   (`Analysis/HitAnalyzer`, `PhysicsTools/NanoAOD`) — never from
   `src/`, which scans thousands of packages.
2. One-file job on `TkAlJpsiX`: confirm both makers run in one
   process (no G4 master conflict), tables written.
3. `uproot` check: raw and corrected columns both present; corrected
   candidate count == raw candidate count; orphan/`joinFlags` column
   accounts for un-refit legs.
4. Cross-check the in-job corrected `m(mu mu K)` against the old
   offline-join output on the same input file — must agree.
5. Repeat per channel and per stream.
