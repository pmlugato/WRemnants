## ADDED Requirements

### Requirement: Multi-Stream Smoke Per Primary Dataset

The Charmonium smoke SHALL run `ALCARECOTkAlJpsiX` and
`ALCARECOTkAlJpsiMuMu` in a single cmsRun; the SingleMuon smoke
SHALL run `ALCARECOTkAlDstToD0Pi`, `ALCARECOTkAlKShortTracks`, and
`ALCARECOTkAlLambdaTracks` in a single cmsRun. Each smoke SHALL
process at least 1000 events, exit 0, and produce non-empty output
for every persisted product listed in each stream's
`_Output_cff.py`.

#### Scenario: Charmonium smoke succeeds

- **WHEN** the Charmonium smoke runs 1000 events of a Run2016H
  Charmonium RAW file
- **THEN** the output contains all persisted J/psi+X products, the
  TkAlJpsiMuMu cloned tracks, VCC, dE/dx maps (added under this
  change), tight muons (added under this change), and track-muon
  association (added under this change); cmsRun exits 0

#### Scenario: SingleMuon smoke succeeds

- **WHEN** the SingleMuon smoke runs 1000 events of a Run2016H
  SingleMuon RAW file
- **THEN** the output contains the D* cloned tracks, dE/dx maps,
  and candidate VCC (added under this change); both V0 streams'
  cloned tracks, VCCs, and dE/dx maps; cmsRun exits 0

### Requirement: Bit-Invariance of Preset B Under the New Gate

Preset B `ALCARECOTkAlJpsiX` output SHALL be bit-identical (per-
channel candidate counts, persisted track counts, persisted muon
counts) to the last smoke's output on the intersection of events
accepted by both the previous tight-muon gate and the new
candidate gate. A check script SHALL read both outputs and assert
equality on the shared subset.

#### Scenario: Bit-invariance holds on shared events

- **WHEN** the multichannel smoke and the previous smoke both
  accept the same event
- **THEN** every J/psi+X collection has the same object count on
  that event

#### Scenario: Bit-invariance fails

- **WHEN** any collection differs
- **THEN** the check script emits a per-collection diff table and
  the launch is blocked

### Requirement: 100-File Multichannel Batch and Projection

The project SHALL submit a 100-file Condor batch per PD and record
per-file `runtime_s_cmsrun`, `peak_rss_mb`, `output_size_mb`, per-
channel candidate counts, persisted track / muon counts, and exit
code. A projection JSON SHALL be emitted with wall days, output
GB, and CPU hours at N in {100, 200, 400} slots, per PD.

#### Scenario: Batches complete and projections recorded

- **WHEN** both 100-file batches exit cleanly for at least 95 / 100
  jobs each
- **THEN** `projections.json` is written with per-PD fields and a
  short markdown report

### Requirement: Loose-vs-Tight Muon Comparison Plot

The project SHALL produce, from production output alone, a plot
set comparing the persisted TkAlJpsiMuMu tight-muon collection to
the persisted J/psi+X loose-muon collection on the same events.
Per observable there SHALL be three overlaid histograms: tight-
only, loose-only, and the loose-and-tight intersection.

#### Scenario: Muon comparison plot produced

- **WHEN** `scripts/btojpsik/plots_preprod_validation.py
  --muon-compare` runs against the smoke outputs
- **THEN** PNG + PDF are produced for pT, eta, phi, charge,
  `isGlobalMuon`, `isTrackerMuon`, `numberOfMatches`, and dimuon
  mass, each with three overlaid histograms and a legend showing
  counts

### Requirement: Per-Channel Kinematic Validation Plots

The project SHALL produce, from the persisted candidate VCCs alone,
histograms of candidate mass, pT, eta, rapidity, helicity angle
(where applicable), and vertex chi2/ndof for every emitted
candidate collection across the four streams.

#### Scenario: All emitted candidate collections plotted

- **WHEN** `plots_preprod_validation.py --kinematics` runs
- **THEN** it emits per-channel plots for B+, B0->K*0, B0->Ks, Bs,
  Lambda_b, psi(2S), Bc, JpsiOnly, TkAlJpsiMuMu, D*, KShort,
  Lambda, reading `daughter()->p4()` and the remapped VCC vertex
  info

### Requirement: Cut-Efficiency Numbers From cmsRun Logs

Cut-efficiency numbers SHALL be sourced from cmsRun stderr
`LogInfo` counter lines that the C++ producers already emit; no
new persisted collection SHALL be introduced for this purpose.

#### Scenario: Counter parse produces yield table

- **WHEN** `plots_preprod_validation.py --cut-efficiency` runs on
  a job's stderr
- **THEN** the parser emits a per-channel yield-reduction table
  identifying the dominant killer cut per channel

### Requirement: V0 Loosening Yield As Offline Side-Study

The V0 loosening yield gain SHALL be measured by a one-shot
offline job that re-runs the CMSSW-default V0 producer on the
same input files as the smoke. It SHALL NOT be measured by adding
a baseline V0 collection to the production stream.

#### Scenario: Yield ratio recorded

- **WHEN** the offline side-study completes on the same 1-file
  smoke inputs
- **THEN** it emits mean loosened / default V0 yield ratios per
  species (KShort, Lambda), and the plotting script consumes them
