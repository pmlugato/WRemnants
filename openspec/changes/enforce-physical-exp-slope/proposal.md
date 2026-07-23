# Change: Enforce physically motivated non-positive slope in AxisExpModel

## Why

`AxisExpModel` currently parameterizes the background as
`exp(lnAmpl + slope · x_m)` with `slope` an unconstrained real. This allows
rising exponentials (`slope > 0`), which are unphysical for combinatorial
background in the B meson mass distribution: higher invariant mass combinations
are kinematically suppressed, so the background must fall (or at most be flat)
with increasing mB. Postfit plots confirm the optimizer can wander into
`slope > 0` in low-statistics cells without meaningful likelihood penalty.

## What Changes

- Re-parameterize the effective slope as `−softplus(raw_slope)` where
  `softplus(z) = log(1 + exp(z))`. This guarantees the effective slope is
  always strictly negative, enforcing a falling background, while keeping the
  optimizer variable (`raw_slope`) unconstrained and the function fully
  differentiable.
- **BREAKING**: Rename slope POIs from `slope_{proc}_{label}` to
  `raw_slope_{proc}_{label}` throughout, since the exposed parameter is now the
  raw pre-softplus value rather than the physical slope.
- Default `xpoidefault` remains `zeros`: `raw_slope = 0` maps to effective slope
  `= −log(2) ≈ −0.693`, a moderately falling background that is a better
  physical prior than the previous flat default.

## Impact

- Affected specs: `axis-exp-poi-model`
- Affected code: `rabbit/rabbit/poi_models/poi_model.py` only
- **BREAKING**: Fit-result HDF5 files and scripts referencing `slope_*` POI
  names will need updating; old fit results are incompatible with new naming.
- No changes to tensor writer, fitter, CLI parsing, `CompositePOIModel`, or
  any other model class.
