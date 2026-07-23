## ADDED Requirements

### Requirement: Axis-driven normalization POI model
The system SHALL provide an `AxisNormModel` class in
`rabbit/rabbit/poi_models/poi_model.py` alongside the existing built-in models.
For each process in proc_spec, the model assigns one independent unconstrained
normalization POI per bin-combination of the caller-specified axes within the
named channel. Processes never share POIs — each gets its own complete set.
All other channels and processes are left at scale factor 1.

The model MUST:
- Accept exactly 3 positional args via
  `parse_args(indata, channel, proc_spec, axes_csv, **kwargs)` and raise
  `ValueError` if a different number is supplied.
- Resolve `proc_spec` as: `"all"` → all processes in the channel; any other
  value → comma-separated list of named processes.
- Raise `ValueError` at construction time if the channel name is absent from
  `indata.channel_info`, if any requested axis name is absent from that
  channel's axes, or if a named process in `proc_spec` is not found in
  `indata.procs`.
- Set `npoi = n_procs * n_cell` where `n_cell` is the product of the sizes of
  the requested axes.
- Produce POI names of the form `norm_{procname}_{ax1name}{i}_{ax2name}{j}...`
  (bytes), ordered by process then by cell index.
- In `compute(poi, full=False)`, return an `rnorm` tensor of shape
  `[nbins, nproc]` where each target process column is scaled by its own
  per-cell POI values (broadcast across all non-selected axes), and all other
  process columns are 1.0.
- Follow the same `allowNegativePOI` / `set_poi_default` / `is_linear`
  contract as the existing built-in models.

#### Scenario: Single named process
- **WHEN** `AxisNormModel` is constructed with channel `btojpsik_stuff`,
  proc_spec `signal`, and axes `bkmm_kaon_pt,bkmm_kaon_charge`
  for a channel with axes `[mass(10), pt(8), eta(7), charge(2)]`
- **THEN** `model.npoi == 16` (1 process × 8 × 2), and only the `signal`
  column of `rnorm` is scaled

#### Scenario: Multiple processes get independent POIs
- **WHEN** proc_spec is `"signal,flatBkg"` and the channel has axes
  `[mass(10), pt(8), eta(7), charge(2)]` with axes `pt,charge` requested
- **THEN** `model.npoi == 32` (2 processes × 8 × 2), `signal` and `flatBkg`
  columns are each scaled by their own 16 independent POI values, and
  no POI is shared between the two processes

#### Scenario: proc_spec=all gives independent POIs per process
- **WHEN** `proc_spec` is `"all"` and the channel has 3 processes and
  `n_cell == 16`
- **THEN** `model.npoi == 48` (3 × 16) and each process column is scaled
  independently by its own 16 POIs

#### Scenario: compute() broadcasts across non-selected axes
- **WHEN** a channel has axes `[mass(10), pt(4), eta(7), charge(2)]` and the
  model is constructed over `pt,charge` only
- **THEN** `compute()` applies each process's (pt, charge) POI value uniformly
  across all 7 eta bins and all 10 mass bins for that (pt, charge) cell

#### Scenario: Wrong number of args raises ValueError
- **WHEN** `parse_args` is called with 2 or 4 positional args
- **THEN** a `ValueError` is raised

#### Scenario: Unknown channel raises ValueError
- **WHEN** the channel name is not present in `indata.channel_info`
- **THEN** a `ValueError` is raised at construction time naming the channel

#### Scenario: Unknown axis raises ValueError
- **WHEN** an axis name in `axes_csv` is not present in the channel's axes
- **THEN** a `ValueError` is raised at construction time naming the axis

#### Scenario: Unknown named process raises ValueError
- **WHEN** `proc_spec` names a process absent from `indata.procs`
- **THEN** a `ValueError` is raised at construction time naming the process

### Requirement: Registration in built-in loader
The system SHALL register `AxisNormModel` in
`rabbit/rabbit/poi_models/helpers.py` so that `--poiModel AxisNormModel ...`
resolves without a full dotted module path.

#### Scenario: Short-name resolution
- **WHEN** `load_model("AxisNormModel", indata, "btojpsik_stuff", "signal,flatBkg", "pt,charge")`
  is called
- **THEN** it returns an `AxisNormModel` instance without raising `ImportError`
  or `KeyError`

### Requirement: SignalNormModel removed
The system SHALL remove `SignalNormModel` from
`wremnants/postprocessing/btojpsik_poi_models.py` after callers migrate to
`AxisNormModel`.

#### Scenario: Removal does not break surviving callers
- **WHEN** `SignalNormModel` is removed
- **THEN** no remaining import or `--poiModel` reference to `SignalNormModel`
  exists in the codebase
