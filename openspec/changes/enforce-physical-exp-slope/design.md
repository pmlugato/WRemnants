## Context

`AxisExpModel` currently uses `exp(lnAmpl + slope · x_m)` with an
unconstrained `slope` POI. The fit is free to choose `slope > 0` (rising
background), which is unphysical for combinatorial background in the B meson
mass. We need to enforce `slope_effective ≤ 0` while keeping the optimizer
variable unconstrained and the computation fully differentiable.

## Goals / Non-Goals

- Goals: guarantee effective slope ≤ 0; preserve differentiability; keep
  `raw_slope` unconstrained (no bounds needed in the optimizer); flat
  background approachable (not hard-walled at zero).
- Non-Goals: changing the amplitude parameterization; modifying other POI
  models; adding a constraint on the magnitude of the slope.

## Decisions

- **Decision: use `slope_eff = −softplus(raw_slope)`**
  - `softplus(z) = log(1 + exp(z)) > 0` for all real `z`, so
    `slope_eff = −softplus(raw_slope) < 0` always.
  - Fully differentiable everywhere; gradient `∂slope_eff/∂raw_slope =
    −sigmoid(raw_slope)` is non-zero everywhere, so the Hessian is
    non-degenerate at any finite `raw_slope`.
  - Flat background (`slope_eff → 0`) is reachable in the limit
    `raw_slope → −∞` but is never an interior point. This is an
    acceptable trade-off given the physics: the background should genuinely
    fall, so a prior away from perfectly flat is appropriate.
  - Alternatives considered:
    - `slope_eff = −exp(raw_slope)`: always negative, but the flat background
      limit (`exp → 0`) requires `raw_slope → −∞`, same trade-off, and the
      gradient vanishes exponentially fast for large negative `raw_slope`,
      creating a slow-convergence region.
    - `slope_eff = −raw_slope²`: always non-positive but zero is an interior
      point and the gradient is zero there, causing a degenerate Hessian at
      the flat case.
    - Hard bound in optimizer (`bounds=(−∞, 0]`): prevents the optimizer from
      exploring positive territory but hits active-bound pathologies in
      BFGS/TNC; rejected.
  - `softplus` is already used in `AxisBernsteinModel` for analogous reasons
    (enforcing non-negativity of Bernstein coefficients), so no new pattern
    is introduced.

- **Decision: rename POIs to `raw_slope_*`**
  - The exposed POI is no longer the physical slope; calling it `slope_*`
    would be misleading in fit-result files and impact summaries.
  - This is a breaking change for scripts and HDF5 files that reference
    `slope_*` by name; the proposal notes this explicitly.

- **Decision: keep `xpoidefault = zeros`**
  - `raw_slope = 0` → `slope_eff = −log(2) ≈ −0.693`, a moderate fall.
    This is a better physical prior than the previous flat default and places
    the initialization well inside the differentiable interior of the model.

## Risks / Trade-offs

- Flat background unreachable exactly → mitigated by physics: the background
  should fall, and the optimizer can reach `slope_eff` arbitrarily close to 0
  by driving `raw_slope` sufficiently negative.
- Breaking POI name change → mitigated by being on a feature branch; no
  external consumers of fit results yet.
- Hessian at `raw_slope → −∞` becomes ill-conditioned (gradient of softplus
  → 0) → mitigated by the default initialization being at `raw_slope = 0`,
  far from this regime.

## Open Questions

- None. Physics motivation is clear; implementation is self-contained.
