## MODIFIED Requirements

### Requirement: Axis-driven exponential background POI model
The system SHALL provide an `AxisExpModel` class in
`rabbit/rabbit/poi_models/poi_model.py` that assigns per-(process, cell)
amplitude POIs and optionally grouped slope POIs, where cells are defined by
caller-specified cell axes. In `compute()` the model produces
`rnorm = exp(lnAmpl_ijk + slope_g · x_m)` for each mass bin `m`, where `g` is
the slope-group index for cell `(i,j,k)`. All other channels and processes are
left at 1.0.

The model MUST:
- Accept 4 or 5 positional args via
  `parse_args(indata, channel, proc_spec, shape_axis, cell_axes_csv[, slope_axes_csv], **kwargs)`
  and raise `ValueError` if fewer than 4 or more than 5 are supplied.
- When `slope_axes_csv` is omitted or equals `cell_axes_csv`, assign one
  independent slope per cell (fully backward-compatible with the 4-arg form).
- Validate that every axis name in `slope_axes_csv` is also present in
  `cell_axes_csv`; raise `ValueError` naming any axis that is not.
- Resolve `proc_spec` identically to `AxisNormModel`: `"all"` → all processes;
  otherwise a comma-separated list of named processes.
- Raise `ValueError` at construction time if the channel is absent from
  `indata.channel_info`, if `shape_axis` or any cell/slope axis name is absent
  from that channel's axes, if `shape_axis` appears in `cell_axes_csv`, or if
  any named process is absent from `indata.procs`.
- Normalize shape-axis bin centers to `[0, 1]` and store them as a constant TF
  tensor for use in `compute()`.
- Set `npoi = n_procs * (n_cell + n_slope_groups)` where `n_cell` is the product
  of cell axis sizes and `n_slope_groups` is the product of slope axis sizes.
  POI names follow `lnAmpl_{procname}_{cell_label}` (one per cell) then
  `slope_{procname}_{slope_label}` (one per slope group) for each process.
- In `compute(poi, full=False)` lay out POIs per process as
  `[lnAmpl_0..n_cell-1, slope_0..n_slope_groups-1]`. Reshape `lnAmpl` with all
  cell-axis sizes and 1 for the shape axis; reshape `slope` with slope-axis
  sizes and 1 for all non-slope cell axes and the shape axis. Return `rnorm` of
  shape `[nbins, nproc]` where the per-bin scaling is
  `exp(lnAmpl + slope * x_m)`, using TF broadcasting to tile slope across
  non-slope cell axes.
- Set `allowNegativePOI = True` (both POIs are unconstrained reals;
  `exp(lnAmpl + slope·x)` is always positive) and `is_linear = False`.
- Set `xpoidefault = tf.zeros([npoi])` so that the default prediction is a
  flat background with unit amplitude (`lnAmpl=0, slope=0` → `exp(0)=1`).

#### Scenario: Construction and POI count with per-cell slopes (4 args)
- **WHEN** `AxisExpModel` is constructed with channel `btojpsik_stuff`,
  proc_spec `bkgExp`, shape_axis `bkmm_jpsimc_mass`, cell_axes
  `bkmm_kaon_pt,bkmm_kaon_eta,bkmm_kaon_charge` for a channel with axes
  `[mass(10), pt(8), eta(7), charge(2)]`
- **THEN** `model.npoi == 224` (1 process × (112 lnAmpl + 112 slope), where
  112 = 8×7×2), POI names begin with `lnAmpl_bkgExp_...` then
  `slope_bkgExp_...`, and `model.n_slope_groups == 112`

#### Scenario: Construction and POI count with grouped slopes (5 args)
- **WHEN** `AxisExpModel` is constructed with the same channel and a 5th arg
  `slope_axes_csv = "bkmm_kaon_eta"` (eta has size 7)
- **THEN** `model.npoi == 119` (1 process × (112 lnAmpl + 7 slope)),
  `model.n_slope_groups == 7`, and slope POI names are
  `slope_bkgExp_bkmm_kaon_eta0` through `slope_bkgExp_bkmm_kaon_eta6`

#### Scenario: compute() produces flat unit rnorm at default POIs
- **WHEN** all POIs are 0.0 (xpoidefault)
- **THEN** the `bkgExp` column of `rnorm` equals 1.0 in every bin and all
  other process columns are 1.0

#### Scenario: compute() applies exponential shape across shape axis
- **WHEN** `lnAmpl = 0.0` and `slope > 0.0` for a given cell
- **THEN** the `bkgExp` column of `rnorm` in that cell increases monotonically
  from mass bin 0 to the last mass bin (rising exponential)

#### Scenario: grouped slope is shared across non-slope cell axes
- **WHEN** slope_axes = ["eta"] and `slope_bkgExp_eta3` is set to 0.5
- **THEN** `rnorm` for `bkgExp` at eta=3 is identical for all pt and charge
  values (slope tiled via broadcasting), while cells at other eta indices use
  their own slope value

#### Scenario: compute() broadcasts across non-cell non-shape axes
- **WHEN** the channel has an axis not in cell_axes and not the shape_axis
- **THEN** the same (lnAmpl, slope) values are applied uniformly across all
  bins of that axis

#### Scenario: slope axis not in cell axes raises ValueError
- **WHEN** `slope_axes_csv` contains an axis name absent from `cell_axes_csv`
- **THEN** a `ValueError` is raised at construction time naming the offending axis

#### Scenario: shape_axis in cell_axes raises ValueError
- **WHEN** `shape_axis` is also listed in `cell_axes_csv`
- **THEN** a `ValueError` is raised at construction time

#### Scenario: Wrong number of args raises ValueError
- **WHEN** `parse_args` is called with 3 or 6 positional args
- **THEN** a `ValueError` is raised

#### Scenario: Unknown channel, axis, or process raises ValueError
- **WHEN** any of channel, shape_axis, a cell axis name, or a process name
  is not found in the tensor
- **THEN** a `ValueError` is raised at construction time naming the missing item
