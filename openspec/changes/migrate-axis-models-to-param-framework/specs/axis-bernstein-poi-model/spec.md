## MODIFIED Requirements

### Requirement: Axis-driven first-order Bernstein background POI model
The system SHALL provide an `AxisBernsteinModel` class in
`rabbit/rabbit/param_models/param_model.py` that assigns two independent
non-negative POIs `(c₀, c₁)` per `(process, cell)`, where cells are defined by
caller-specified axes. In `compute()` the model produces

    rnorm(x_m) = c₀ · (1 − x_m) + c₁ · x_m

for each bin `m` of the shape axis, where `x_m ∈ [0, 1]` is the normalized
bin center. All other channels and processes are left at 1.0.

The model MUST:
- Accept exactly 4 positional args via
  `parse_args(indata, channel, proc_spec, shape_axis, cell_axes_csv, **kwargs)`
  and raise `ValueError` if a different number is supplied.
- Resolve `proc_spec` identically to `AxisNormModel` and `AxisExpModel`:
  `"all"` → all processes; otherwise a comma-separated list of named processes.
- Raise `ValueError` at construction time if the channel is absent from
  `indata.channel_info`, if `shape_axis` or any cell axis name is absent from
  that channel's axes, if `shape_axis` appears in `cell_axes_csv`, or if any
  named process is absent from `indata.procs`.
- Normalize shape-axis bin centers to `[0, 1]` and store them as a constant
  TF tensor for use in `compute()`.
- Set `npoi = n_procs * 2 * n_cell` where `n_cell` is the product of cell axis
  sizes, and `npou = 0`. Parameter names (`self.params`) follow
  `c0_{procname}_{cell_label}` then `c1_{procname}_{cell_label}` for each
  process.
- In `compute(param, full=False)` lay out parameters per process as
  `[c0_cell0..c0_cell(n-1), c1_cell0..c1_cell(n-1)]`. Reshape each block to
  the cell axes (broadcasting over the shape axis) and return `rnorm` of shape
  `[nbins, nproc]` where the per-bin scaling is `c₀·(1−x_m) + c₁·x_m`.
  Non-negativity of `c₀` and `c₁` is enforced via `softplus` applied inside
  `compute()`.
- Set `allowNegativeParam = True` (raw optimizer parameters are unconstrained
  reals; `softplus` is applied inside `compute()` to recover non-negative
  coefficients) and `is_linear = False`.
- Set `xparamdefault` such that `softplus(xparamdefault) = 1`, i.e.
  each raw parameter initialised to `softplus⁻¹(1) ≈ 0.5413`, giving
  `c₀=c₁=1` at initialization (flat unit background).

#### Scenario: Construction and POI count
- **WHEN** `AxisBernsteinModel` is constructed with channel `ch`, proc_spec
  `bkgBernstein`, shape_axis `mass`, cell_axes `eta,pt,charge` for a channel
  with axes `[mass(10), eta(7), pt(8), charge(2)]`
- **THEN** `model.npoi == 224` (1 process × 2 × 7 × 8 × 2), `model.npou == 0`,
  and parameter names begin with `c0_bkgBernstein_...` followed by
  `c1_bkgBernstein_...`

#### Scenario: compute() at default produces flat unit rnorm
- **WHEN** all parameters are at their default (raw value `≈ 0.5413`,
  giving `c₀=c₁=1` via softplus, i.e. `xparamdefault`)
- **THEN** the `bkgBernstein` column of `rnorm` equals 1.0 in every bin and
  all other process columns equal 1.0

#### Scenario: compute() with c₀ > c₁ gives falling background
- **WHEN** `c₀ = 2.0` and `c₁ = 0.5` for a given cell
- **THEN** the `bkgBernstein` column of `rnorm` in that cell decreases
  monotonically from mass bin 0 to the last mass bin

#### Scenario: compute() with c₀ < c₁ gives rising background
- **WHEN** `c₀ = 0.5` and `c₁ = 2.0` for a given cell
- **THEN** the `bkgBernstein` column of `rnorm` in that cell increases
  monotonically from mass bin 0 to the last mass bin

#### Scenario: compute() broadcasts across non-cell non-shape axes
- **WHEN** the channel has an axis not in cell_axes and not the shape_axis
- **THEN** the same `(c₀, c₁)` values are applied uniformly across all bins
  of that axis

#### Scenario: shape_axis in cell_axes raises ValueError
- **WHEN** `shape_axis` is also listed in `cell_axes_csv`
- **THEN** a `ValueError` is raised at construction time

#### Scenario: Wrong number of args raises ValueError
- **WHEN** `parse_args` is called with fewer than 4 or more than 4 positional
  args
- **THEN** a `ValueError` is raised

#### Scenario: Unknown channel, axis, or process raises ValueError
- **WHEN** any of channel, shape_axis, a cell axis name, or a process name
  is not found in the tensor
- **THEN** a `ValueError` is raised at construction time naming the missing item

### Requirement: AxisBernsteinModel registration in built-in loader
The system SHALL register `AxisBernsteinModel` in
`rabbit/rabbit/param_models/helpers.py` so that `--paramModel
AxisBernsteinModel ...` resolves without a full dotted module path.

#### Scenario: Short-name resolution
- **WHEN** `load_model("AxisBernsteinModel", indata, "ch", "bkg", "mass",
  "eta,pt")` is called
- **THEN** it returns an `AxisBernsteinModel` instance without raising
  `ImportError` or `KeyError`
