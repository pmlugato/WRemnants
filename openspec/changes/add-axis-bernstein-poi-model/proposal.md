# Change: Add AxisBernsteinModel POI model

## Why

The exponential background model (`AxisExpModel`) can suffer from near-degenerate
Hessians in cells where signal and background shapes are similar. A first-order
Bernstein polynomial is a structurally different parameterization: it describes
the background as a linear interpolation between endpoint values `c₀` (rate at
the low mass edge) and `c₁` (rate at the high mass edge). Non-negativity is
guaranteed by construction (both coefficients constrained ≥ 0 via x²
reparameterization), and the flat-background default (`c₀=c₁=1`) is an
interior point in parameter space with non-zero gradient, avoiding the boundary
degeneracy of the original exponential form.

## What Changes

- Add `AxisBernsteinModel` to `rabbit/rabbit/poi_models/poi_model.py`. The
  model assigns two POIs `(c₀, c₁)` per `(process, cell)` and produces
  `rnorm(x_m) = c₀·(1−x_m) + c₁·x_m` where `x_m ∈ [0,1]` is the normalized
  center of shape-axis bin `m`. Both coefficients must be ≥ 0; the x²
  reparameterization (`allowNegativePOI=False`) enforces this. Default
  `c₀=c₁=1` gives a flat unit background.
- Register `AxisBernsteinModel` in `rabbit/rabbit/poi_models/helpers.py` for
  short-name resolution via `--poiModel`.
- No changes to tensor writer, fitter, parsing, or `CompositePOIModel`.

## Impact

- Affected specs: `axis-bernstein-poi-model` (new capability)
- Affected code: `rabbit/rabbit/poi_models/poi_model.py`,
  `rabbit/rabbit/poi_models/helpers.py`
