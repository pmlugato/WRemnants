# Design: AxisExpModel grouped-slope support

## Context

AxisExpModel currently stores slope as a flat vector of length `n_cell`, reshaped
to `cell_reshape` in `compute()`. To support a coarser slope granularity we need
a different reshape for slope vs. lnAmpl.

## Goals / Non-Goals

- Goals: shared slope per any caller-specified axis subset; fully backward-
  compatible default; single-file change; no new dependencies
- Non-Goals: per-process slope grouping; runtime-changeable grouping;
  cross-channel slope sharing

## Decisions

**Decision: TF broadcasting instead of gather/scatter**

The slope tensor is given a `slope_cell_reshape` that has the actual axis size
for slope axes and 1 for all other cell axes (plus 1 for the shape axis). When
added to the per-cell `lnAmpl` (shaped with all cell axis sizes), TF
automatically tiles the slope across the non-slope dimensions. No explicit index
mapping is needed.

Example: channel axes `[mass(10), eta(28), pt(8), charge(2)]`, slope axis `eta`:
- `cell_reshape`        = `[1, 28, 8, 2]`  (lnAmpl)
- `slope_cell_reshape`  = `[1, 28, 1, 1]`  (slope)
- `x_reshaped`          = `[10,  1, 1, 1]`  (shape axis)
- `exp(lnAmpl + slope * x)` has shape `[10, 28, 8, 2]` via broadcasting ✓

Alternatives considered:
- `tf.gather`: explicit cell→group index array; correct but more code and a
  non-differentiable indexing step (TF gather IS differentiable, but broadcasting
  is simpler and equally efficient here).
- Passing a precomputed `slope_groups` array from the caller: pushes
  analysis-specific logic out of the model; rejected in favour of axis-name
  resolution inside the model using `indata.channel_info`.

**Decision: optional 5th positional arg, not a keyword flag**

Consistent with how AxisNormModel and AxisExpModel already accept all
configuration as positional args to `parse_args`. A keyword flag would require
changes to parsing infrastructure.

**Decision: slope axes must be a subset of cell axes**

Ensures the slope reshape dimensions are a subset of the lnAmpl reshape
dimensions, guaranteeing correct broadcasting. Validated at construction time
with a descriptive ValueError.

## Parameter layout (per process)

```
poi = [ lnAmpl_cell0, ..., lnAmpl_cell(n_cell-1),
        slope_group0, ..., slope_group(n_slope_groups-1) ]
```

Length = `n_cell + n_slope_groups` per process.

## Risks / Trade-offs

- Slope sharing assumes the background slope is truly homogeneous across the
  non-slope axes within each slope group. If this is violated (e.g. the slope
  varies with pT), the fit absorbs the residual into lnAmpl, potentially biasing
  per-cell amplitudes slightly. This is acceptable given the physics motivation.
- Changing `npoi` is a breaking change for any stored fit results that index
  POIs by position. Results from the old per-cell-slope fits are not compatible
  with the new layout.
