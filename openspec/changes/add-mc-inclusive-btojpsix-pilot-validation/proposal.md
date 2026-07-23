# Change: Pilot validation of the 108-job MC B→J/ψ+X pilot before full production

## Why

Before we submit the full ~15 fb⁻¹-equivalent MC inclusive B→J/ψ+X production
to the CMS grid, we need to convince ourselves the pilot batch (~108 jobs
produced by david_w at
`/ceph/submit/data/user/d/david_w/mz/mc/inclusive_btojpsix_2016postvfp/`,
using his patched tree with the split=1 output-module fix already in place) is
physics-healthy end-to-end. David is the operator for the full submission, so
his tree and operator identity remain authoritative; we only need our own tree
patched enough to (a) reflect the split=1 physics fix so nothing regresses on
our side, and (b) redirect `DEST_DIR` to the shared group MC store
(`/ceph/submit/data/group/cms/store/mc/`) that the full production will write
to. The OSG rehearsal (cluster 3228825, four days idle) is being parked in
favor of routing the full production through the CMS grid via David.

## What Changes

- Reconcile the five fixes in
  `/work/submit/david_w/ZMass/BtoJpsiX_MCprod/pietro_fixes/` into
  `condor/mc_inclusive_btojpsix_2016postvfp/` (patch or hand-merge from
  `fixed_files/`). **David submits, so his operator identity — proxy path,
  `+AccountingGroup`, transient site fences — is authoritative**; we keep our
  local versions of those fields only for our own dry-runs and do not push
  them back. The one operator field we do change is `DEST_DIR`, which becomes
  the shared group MC store `/ceph/submit/data/group/cms/store/mc/...` (final
  subpath TBD in §1); we communicate this DEST_DIR change to David manually.
- **No payload push-back to David unless we find bugs.** David's tree already
  has the split=1 fix baked into his running payload. We only rebuild + ship
  him a new tarball if pilot validation surfaces a bug beyond the five fixes
  he already applied.
- Run the reconcile script against David's ~108 pilot outputs at his ceph URL
  (read-only): produce a per-job JSON census (success fraction,
  wall/RSS/scratch p99, OOM tally after the exit-137 fix).
- **Physics inspection of the pilot ALCARECOs — B⁺→J/ψK⁺ path only for now**:
  event counts, ALCARECO branch presence, `TkAlJpsiXBPlusResonances` candidate
  multiplicity per event, per-leaf hit counts.
- **Stage-2 CVH refit on the pilot** — reuse the machinery from
  `add-jpsi-x-stage2-bplus-cvh` (two-track refit for the J/ψ, single-track
  refit for the kaon). Run with `applyIdealGeometry=False` (nominal alignment).
  Sanity-check: ditrack J/ψ mass median lands at ~3.08 GeV; the pilot was
  produced with split=1 correct so this should just work.
- **B⁺ candidate formation from CVH tracks with a joint kinematic-vertex
  fit** — combine the CVH two-track J/ψ and CVH single-track K± into a B⁺
  candidate via a kinematic vertex fit (reuse whatever the preset-B / JpsiX
  chain already provides for this; if a new small tool is needed, land it
  under `scripts/btojpsik/mc_pilot_validation/`). The vertex-fit χ²/ndof is
  the primary selection knob to isolate real B⁺'s from combinatoric.
- **Simple signal selections**: J/ψ mass window, K± quality, B⁺ mass window,
  and **the joint-fit χ²/ndof cut**. Nothing beyond this for the pilot.
- **Gen matching**: identify true B⁺→J/ψK⁺ candidates via ΔR at the two muons
  + bachelor; split reco distributions into matched vs non-matched to
  visually separate signal from combinatoric.
- **Data overlay** against
  `/ceph/submit/data/user/p/pmlugato/mz/alcareco/full_2016postvfp/charmonium/Run2016*`:
  m(J/ψ) and m(B⁺) shape comparison after selections. No hard MC-vs-data
  agreement threshold — the point is a visual sanity check, not a
  measurement. (Once the alcareco→nano pipeline lands we'll switch to the
  standard histmakers; for now this is a one-shot RDataFrame pass.)
- Publish plots + report into the pilot validation slide deck.

## Impact

- **Affected specs**: `mc-production` (adds a pilot-validation acceptance gate
  before the full 15 fb⁻¹ MC submission).
- **Affected code / directories**:
  - `condor/mc_inclusive_btojpsix_2016postvfp/` — patched with david_w's five
    fixes; `DEST_DIR` retargeted to `/ceph/submit/data/group/cms/store/mc/...`;
    David's operator identity kept authoritative (we don't overwrite our own
    local fields for dry-run either — just don't propagate).
  - `condor/mc_inclusive_btojpsix_2016postvfp_osg/` — parked (README banner
    added; no code changes needed).
  - `slides/mc_inclusive_btojpsix_2016postvfp.tex` — add pilot-validation
    section (data-quality + physics-quality plots).
  - New analysis directory (proposed):
    `scripts/btojpsik/mc_pilot_validation/` — housing the pilot-scanner, the
    joint-fit wrapper (if needed), gen-matching helper, and MC-vs-data
    overlay plotters. Kept out of the main btojpsik histmaker path; this is a
    one-shot validation, not a permanent pipeline. When the alcareco→nano
    pipeline lands, we retire this dir in favor of the standard histmakers.
- **Approval gate**: user signs off on the pilot validation report before the
  full ~15 fb⁻¹ MC submission is issued.
- **Non-goals**:
  - No changes to the physics content of the fullchain cfg beyond the split=1
    output-module fix (already in David's tree).
  - No mZ tensor / rabbit-fit changes.
  - No OSG rehearsal rework — parked.
  - No new preset development; we run on preset-B alignment.
  - No hard MC-vs-data agreement threshold; visual comparison only.
  - No alcareco→nano converter in this change; deferred.
