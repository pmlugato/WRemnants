## ADDED Requirements

### Requirement: Axis-driven exponential background POI model
The system SHALL provide an `AxisExpModel` class in
`rabbit/rabbit/poi_models/poi_model.py` that assigns two independent
unconstrained POIs `(A, B)` per (process, cell-bin-combination), where the
cell bins are defined by the caller-specified cell axes. In `compute()` the
model produces `rnorm = A·e^{−B·x_m}` for each mass bin `m` of the shape axis,
independently per process and per cell. All other channels and processes are
left at 1.0.

The model MUST:
- Accept exactly 4 positional args via
  `parse_args(indata, channel, proc_spec, shape_axis, cell_axes_csv, **kwargs)`
  and raise `ValueError` if a different number is supplied.
- Resolve `proc_spec` identically to `AxisNormModel`: `"all"` → all processes
  in the channel; otherwise a comma-separated list of named processes.
- Raise `ValueError` at construction time if the channel is absent from
  `indata.channel_info`, if `shape_axis` or any cell axis name is absent from
  that channel's axes, if `shape_axis` appears in `cell_axes_csv`, or if any
  named process is absent from `indata.procs`.
- Normalize shape-axis bin centers to `[0, 1]` and store them as a constant
  TF tensor for use in `compute()`.
- Set `npoi = n_procs * 2 * n_cell` where `n_cell` is the product of cell axis
  sizes. POI names follow `expA_{procname}_{cell_label}` then
  `expB_{procname}_{cell_label}` for each (process, cell) pair.
- In `compute(poi, full=False)`, return `rnorm` of shape `[nbins, nproc]`
  where for each target process the per-mass-bin scaling is
  `A_ijk · e^{−B_ijk · x_m}`, broadcast across all non-selected axes, and
  all other processes are 1.0.
- Set `is_linear = False` and default `allowNegativePOI = False` (both A and B
  constrained positive via the `x²` reparametrization).
- Follow the same `set_poi_default` / `xpoidefault` contract as other built-in
  models.

#### Scenario: Construction and POI count
- **WHEN** `AxisExpModel` is constructed with channel `btojpsik_stuff`,
  proc_spec `bkgExp`, shape_axis `bkmm_jpsimc_mass`, cell_axes
  `bkmm_kaon_pt,bkmm_kaon_eta,bkmm_kaon_charge` for a channel with axes
  `[mass(10), pt(8), eta(7), charge(2)]`
- **THEN** `model.npoi == 224` (1 process × 2 × 8 × 7 × 2) and POI names
  begin with `expA_bkgExp_...` and `expB_bkgExp_...`

#### Scenario: compute() produces mass-dependent rnorm
- **WHEN** all A POIs are set to 2.0 and all B POIs are set to 0.0
- **THEN** the `bkgExp` column of `rnorm` equals 2.0 in every bin (flat
  exponential with zero decay rate), and all other process columns are 1.0

#### Scenario: compute() applies exponential decay across shape axis
- **WHEN** A = 1.0 and B > 0.0 for a given cell
- **THEN** the `bkgExp` column of `rnorm` in that cell decreases monotonically
  from mass bin 0 to the last mass bin

#### Scenario: compute() broadcasts across non-cell non-shape axes
- **WHEN** the channel has an axis not in cell_axes and not the shape_axis
- **THEN** the same (A, B) POI values are applied uniformly across all bins
  of that axis

#### Scenario: shape_axis in cell_axes raises ValueError
- **WHEN** `shape_axis` is also listed in `cell_axes_csv`
- **THEN** a `ValueError` is raised at construction time

#### Scenario: Wrong number of args raises ValueError
- **WHEN** `parse_args` is called with 3 or 5 positional args
- **THEN** a `ValueError` is raised

#### Scenario: Unknown channel, axis, or process raises ValueError
- **WHEN** any of channel, shape_axis, a cell axis name, or a process name
  is not found in the tensor
- **THEN** a `ValueError` is raised at construction time naming the missing item

### Requirement: AxisExpModel registration in built-in loader
The system SHALL register `AxisExpModel` in
`rabbit/rabbit/poi_models/helpers.py` so that `--poiModel AxisExpModel ...`
resolves without a full dotted module path.

#### Scenario: Short-name resolution
- **WHEN** `load_model("AxisExpModel", indata, "ch", "bkgExp", "mass", "pt,eta")` is called
- **THEN** it returns an `AxisExpModel` instance without raising `ImportError`
  or `KeyError`

### Requirement: btojpsik tensor writer exp_poi background mode
The system SHALL add `exp_poi` as a `--bkgModel` choice in
`scripts/rabbit/btojpsik_tensor.py`. When selected:
- A `flatBkg` process is added with an all-ones histogram scaled by a per-cell
  flat yield estimate.
- A `bkgExp` process is added with an all-ones histogram scaled by a per-cell
  exponential yield estimate.
- No normalization NOIs are added for either background process.
- The signal process is added with `signal=True`.

#### Scenario: exp_poi mode writes all-ones templates
- **WHEN** `--bkgModel exp_poi --signalNormPOI` is passed
- **THEN** the tensor contains `flatBkg` and `bkgExp` processes whose templates
  are all-ones (times a per-cell scalar) and no `norm_flatBkg_*` or
  `norm_bkgExp_*` NOIs exist in the tensor
