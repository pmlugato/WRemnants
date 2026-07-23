# Change: Add grouped-slope support to AxisExpModel

## Why

The current AxisExpModel assigns two free POIs per (process, cell): `lnAmpl`
and `slope`. With 448 cells (28η × 8pT × 2charge), 10 mass bins as constraints,
and a per-cell signal norm also free, this produces a near-degenerate Hessian:
the 2×2 sub-block for (norm_signal, lnAmpl_bkgExp) in each cell has ~97%
correlation, causing Cholesky decomposition to fail. The root cause is that 3
free parameters per cell (norm, lnAmpl, slope) are too loosely constrained by
10 mass bins when the signal and background shapes are nearly proportional.

The exponential slope is physically driven by detector acceptance and
combinatorial background kinematics, which are primarily η-dependent. Sharing
one slope per η bin (28 slopes instead of 448) is well-motivated by physics and
breaks the degeneracy without sacrificing per-cell amplitude freedom.

## What Changes

- `AxisExpModel` gains an optional 5th positional CLI argument `slope_axes_csv`:
  a comma-separated subset of `cell_axes` naming the axes over which the slope
  varies independently. Cells that share the same slope-axis coordinates share
  one slope POI.
- When `slope_axes_csv` is omitted (4 args), behavior is identical to today:
  one slope per cell (fully backward-compatible).
- `npoi` changes from `n_procs * 2 * n_cell` to
  `n_procs * (n_cell + n_slope_groups)` where
  `n_slope_groups = prod(slope axis sizes)`.
- Slope POI names are keyed by slope-axis indices only; lnAmpl names are
  unchanged (per cell).
- In `compute()`, the slope tensor is reshaped with 1s for non-slope-axis
  dimensions so TF broadcasting tiles it across the remaining cell axes — no
  explicit gather/scatter required.
- The parameterization remains `exp(lnAmpl + slope·x)` with
  `allowNegativePOI = True`; this change does not alter the functional form.
- This also corrects the spec: the existing `add-axis-exp-poi-model` delta
  describes the old `A·exp(-B·x)` form; the MODIFIED requirement here reflects
  the current `exp(lnAmpl + slope·x)` implementation.

## Impact

- Affected specs: `axis-exp-poi-model`
- Affected code: `rabbit/rabbit/poi_models/poi_model.py` only
- No changes to tensor writer, fitter, CompositePOIModel, or CLI parsing
  infrastructure
