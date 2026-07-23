# Change: Improve CVH joint-refit convergence for B$^+$→J/ψK+ stage-2

## Why

The Run2016H investigation (`cvh_vs_alcareco_jpsi_mass_investigation.md`,
slides `slides/cvh_vs_alcareco_jpsi_mass_investigation.pdf`) on 81,948
matched B+ candidates established that the CVH unconstrained joint refit
`Jpsi_*` produces a dimuon mass distribution whose Gaussian σ on
[2.95, 3.25] is 80 MeV vs ~42 MeV for the raw input tracks and the CMSSW
`KinematicParticleVertexFitter`. The same joined tree shows:

1. **No event meets the documented convergence criterion** —
   `edmval < 10^{-5}` is satisfied by **0 / 81948** candidates. Median
   `edmval` $= 3.6\times 10^5$, q99 $= 6.8\times 10^6$. 20 % of events
   hit the hard `niters = 10` cap; the rest leave via other break paths.
2. **`chisqval` is far from $\chi^2/\text{ndof} \sim 1$** — median 5780,
   q95 26687 at $\sim$500 fit parameters.
3. **Per-leg q/p shifts are an order of magnitude larger than the
   standard CMSSW Kalman per-track $\sigma(q/p)$**. MAD of $\Delta(q/p)$
   per leg is $\sim 10 \times 10^{-3}$ GeV$^{-1}$ vs $\sim 10^{-3}$
   GeV$^{-1}$ typical. The mass-formula propagation
   $\delta m / m \approx \tfrac12 \sum \delta p_i / p_i$ accounts
   quantitatively for the broadening within a factor 2.
4. **A 33 % catastrophic-outlier population** populates a tail running
   from $m(\mu\mu) = 1.3$ to 102 GeV. For those same events the raw
   input track sum and the CMSSW vertex fit remain pristine, so this is
   a refit-side pathology, not a sample property.

H1 (CVH algorithmic / convergence) is the dominant cause; H2 (material
propagation) and H3 (joint-vertex scrambling of combinatorial μμ) are
both ruled out by the agreement between `Jpsitrk_mass` and `Jpsikin_mass`
on the same matched events. Diagnosing the convergence failure and
fixing it is the highest-leverage next step.

## What changes

This change is a single coordinated investigation + iteration on the
CVH refit. It both adds the diagnostic instrumentation needed to
characterise the failure and applies the three most-likely candidate
fixes, so that the impact of each can be measured on the same event
sample. Nothing is committed to production until the experiment matrix
has been read.

### 1. Per-iteration diagnostic dump

Today only the final-iter values of `chisqval`, `edmval`, `edmvalref`,
and `niter` reach the tree, and only after the icons=1 (mass-constrained)
phase. The icons=0 (unconstrained) phase writes only `niter_cons0` and
`edmval_cons0`. To understand whether the fit is failing to make
progress, plateauing, or oscillating, we need the per-iteration sequence
for a small debug sample.

Add a debug mode (`debugPerIterDump = cms.bool(False)`, default off) to
`ResidualGlobalCorrectionMakerTwoTrackG4e.cc`. When on, writes one
TTree-friendly vector branch per quantity:

```
chisqval_iter        std::vector<double>   per icons=0 + icons=1 iterations
edmval_iter          std::vector<double>
deltachisqval_iter   std::vector<double>
mu_qoverp_iter       std::vector<double>   per-iter q/p of each muon (2 entries per iter)
Jpsi_mass_iter       std::vector<double>   per-iter dimuon mass
```

These are switched on only for the debug runs in the experiment matrix
below (kept off for production to avoid bloating the file).

### 2. Hit-bookkeeping fix-up

The current outputs ship two branches that are misleading: in the
single-track `ResidualGlobalCorrectionMakerG4e.cc` the
`nValidHitsFinal` and `nValidPixelHitsFinal` branches are declared and
initialised to 0 but **never incremented**. The joined tree therefore
reports `Kbach_nValidHitsFinal = 0` for all 81948 events, which suggests
(falsely) that all kaon hits are dropped. The actual kaon refit uses
`nValidHits` (median 9, q01 1) — those hits go in and stay in.

Either fill the `*Final` counters properly (preferred — they let us see
if the alignment-pass quality cut ever rejects hits) or remove them
from the branch list. Pick one; do not ship dead branches.

The dimuon producer increments `nValidHitsFinal` per hit but does not
write it to a branch. Add per-muon branches `Muplus_nvalidFinal`,
`Muminus_nvalidFinal` (and pixel variants) — symmetric with the existing
`nvalid` and `nvalidpixel` branches — so we can verify there too.

### 3. Three candidate fixes, applied as a configurable matrix

Each is a one-line constant change in
`ResidualGlobalCorrectionMakerTwoTrackG4e.cc`; they are kept as
configurable parameters under module psets so the experiment matrix can
sweep over them without rebuilding.

| Parameter | Default | Range to test | Reason |
|---|---|---|---|
| `nIters` | 10 | {10, 20, 50, 100} | The 20 % cap-hit fraction may mask events that need more iterations. Bumping from 10 → 100 directly tests this. |
| `edmConvergence` | 1e-5 | {1e-5, 1e-3, 1e-2} | Standard GBL convergence is $\Delta\chi^2/\text{ndof} < 10^{-2}$. The unitless `edmval` threshold may be too tight. |
| `useStartingState` | "perigee" | {"perigee", "midPropagated"} | Today the linearisation starts at each muon's perigee parameters; for high-IP B+ topologies that is far from the joint vertex. Try starting from the propagated state at the joint vertex midpoint. |

The default values remain those used in the published 81948-event
sample, so the same job script + same input file can reproduce the
baseline before any change is rolled forward.

### 4. Experiment matrix and metric set

All checks below SHALL use the nominal CMSSW magnetic field
(`useScalarPot3D=False`, `useOpera3D=False`); the scalar-potential and
Opera-3D field models are explicitly out of scope for this convergence
study so that any change in behaviour is attributable to the fit
algorithm and not to a swapped field model.

Run a 4 × 3 × 2 = 24-point sweep (`nIters` × `edmConvergence`
× `useStartingState`) on a single Run2016H ALCAReco file
(`02E1D3B8-5A93-E611-AFF9-02163E013446.root`, ~1500 candidates), with
the per-iter debug dump on, and report:

- fraction of events meeting `edmval < 10^{-5}`
- median and q95 of `chisqval`, `chisqval / nParms`
- σ of `Jpsi_mass` on [2.95, 3.25] and core [3.05, 3.15]
- fraction with `|Jpsi_mass - 3.0969| > 0.2` (the 33 % outlier tail)
- MAD of per-leg Δp_T/p_T and Δ(q/p) vs the raw input track
- CPU time per event

The winning configuration becomes the new default in the
`runCvhBplusJpsiK.py` driver. If no configuration kills both the
convergence failure and the 33 % outlier tail, the diagnostic dump
isolates which iteration the fit diverges in, which becomes the next
investigation step.

### 5. Cross-driver sanity check — static-only review

**Downgraded 2026-06-29**: a per-candidate dynamic comparison was
originally planned, requiring matched J/ψ candidates in both
`ALCARECOTkAlJpsiMuMu` and `ALCARECOTkAlJpsiX`. No dedicated
`TkAlJpsiMuMu` output exists for our Run2016H sample, so the dynamic
check is not feasible without producing one — out of scope for this
change. The downgraded review:

- Reads both driver scripts (`runCvhJpsi.py`,
  `runCvhBplusJpsiK.py mode=dimuon`) and documents any pset differences
  that touch the dimuon maker (beam-spot source, particle list,
  magnetic-field label, propagation pT limit, mass-constraint flag,
  candidate-input wiring).
- Notes those differences as a short table in the investigation slides.
- Does NOT modify `runCvhJpsi.py`.

The per-candidate join harness `scripts/btojpsik/compare_cvh_drivers.py`
is written but parked — invocable if a TkAlJpsiMuMu ALCAReco file pair
surfaces later. The §4 matrix interpretation is no longer gated on a
cross-driver PASS; baseline correctness is validated instead by the
§1.7 bit-identical baseline check (defaults reproduce the published
81,948-event Run2016H sample).

### 6. Documentation

Update `slides/cvh_vs_alcareco_jpsi_mass_investigation.pdf` with a new
"experiment matrix results" section once the matrix runs, and a
cross-driver-sanity-check call-out (pass/fail line per branch).

## Affected files

- `CMSSW_15_0_19_patch2/src/Analysis/HitAnalyzer/plugins/ResidualGlobalCorrectionMakerTwoTrackG4e.cc`
  (new parameters `nIters`, `edmConvergence`, `useStartingState`,
  `debugPerIterDump`; new per-iter branch fills under the debug flag;
  new per-muon `nvalidFinal` branches)
- `CMSSW_15_0_19_patch2/src/Analysis/HitAnalyzer/plugins/ResidualGlobalCorrectionMakerG4e.cc`
  (either fill or drop `nValidHitsFinal` / `nValidPixelHitsFinal`)
- `CMSSW_15_0_19_patch2/src/Analysis/HitAnalyzer/test/runCvhBplusJpsiK.py`
  (expose the three new parameters as command-line knobs;
  `debug=False` switches on the per-iter dump for the matrix runs)
- `scripts/btojpsik/run_cvh_experiment_matrix.sh` (new) — 24-point
  driver, one cmsRun per config, joined output filename encodes the
  parameter tuple
- `scripts/btojpsik/plot_cvh_experiment_matrix.py` (new) — reads the
  24 joined files, produces the 6-metric grid, writes JSON for the
  slides
- `slides/cvh_vs_alcareco_jpsi_mass_investigation.tex` (extended with
  a results section after the matrix runs)

## Out of scope

- Refactoring the GBL formulation or the Geant4e propagator. Those are
  too risky to touch before the matrix has identified which level of
  the fit is failing.
- Changes to the kaon (single-track) refit. The kaon side is not
  implicated in the dimuon broadening; if the matrix reveals a kaon-side
  fit pathology this becomes a separate proposal.
- Production reprocessing. Once a winning configuration is identified,
  rolling it forward to the 100k-event production is a separate
  change-id (it touches sample provenance, not the refit code).

## Open questions

1. Is `edmval`'s unit really $-\delta\chi^2$ (where the GBL convergence
   threshold is dimensionless O(1)), or is it $-\delta\chi^2/\text{ndof}$
   (where $10^{-5}$ is sensible)? Reading the code suggests the former
   (the comment at line 2469 reads `Sparse GBL deltachi2: dense g^T dx
   + 0.5 dx^T H dx with g=2Fs^TVinvr, ...`). If so the threshold is
   ~5 orders of magnitude too tight.
2. Why does the iter loop break for 80 % of events before niter=10 with
   `edmval` still ~$10^5$? There must be another break path that the
   investigation missed; the per-iter dump will pin it down.
3. Are the 33 % catastrophic-tail events all in one phase-space corner
   (e.g. forward muons, low-pT bachelor, large dimuon Lxy)? Or scattered?
   The matrix's per-event diagnostic output should answer this.
4. Does an analogous "1-track-at-a-time" CVH fit (the original
   `ResidualGlobalCorrectionMakerG4e.cc` applied to each muon
   independently) reproduce the same per-leg Δ(q/p) ~ 10⁻² shift? If
   yes, the joint-vertex aspect is irrelevant. If no, the joint
   formulation itself is the issue.
