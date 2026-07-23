## 1. Implementation

- [x] 1.1 In `AxisExpModel.__init__`, rename all internal slope POI name
      construction from `slope_{proc}_{label}` to `raw_slope_{proc}_{label}`
      (the `self.pois` array).
- [x] 1.2 In `AxisExpModel.compute()`, replace `b * x_reshaped` with
      `-tf.math.softplus(b) * x_reshaped` so the effective slope is always
      strictly negative.
- [x] 1.3 Update the `allowNegativePOI` docstring comment in `__init__` to
      reflect that `raw_slope` is unconstrained but the effective slope is
      always negative.
- [x] 1.4 Update the class docstring to document the new functional form and
      POI naming.

## 2. Validation

- [x] 2.1 Syntax-check the modified file:
      `python -m py_compile rabbit/rabbit/poi_models/poi_model.py` — PASSED
- [ ] 2.2 Verify the model constructs and produces a strictly falling `rnorm`
      at default POIs — requires container (TF not available outside);
      mathematically guaranteed: softplus(b)>0 so −softplus(b)·x_m≤0.
- [ ] 2.3 Verify large negative raw_slope gives nearly flat rnorm — same
      as 2.2; softplus(−10)≈0.0000454, so effective slope ≈ 0.
- [ ] 2.4 Verify no finite raw_slope produces rising rnorm — same as 2.2;
      follows directly from softplus always being positive.
