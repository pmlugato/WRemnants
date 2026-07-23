## 1. Axis model attribute renames in `param_model.py`

- [x] 1.1 `AxisNormModel.__init__`: rename `self.pois` → `self.params`,
      `allowNegativePOI` → `allowNegativeParam`, `set_poi_default` →
      `set_param_default`; add `self.npou = 0` after `self.npoi`.
- [x] 1.2 `AxisNormModel.compute`: rename argument `poi` → `param` and all
      internal uses.
- [x] 1.3 `AxisExpModel.__init__`: rename `self.pois` → `self.params`,
      `self.allowNegativePOI` → `self.allowNegativeParam`,
      `self.xpoidefault` → `self.xparamdefault`; add `self.npou = 0`.
- [x] 1.4 `AxisExpModel.compute`: rename argument `poi` → `param` and all
      internal uses.
- [x] 1.5 `AxisBernsteinModel.__init__`: rename `self.pois` → `self.params`,
      `self.allowNegativePOI` → `self.allowNegativeParam`,
      `self.xpoidefault` → `self.xparamdefault`; add `self.npou = 0`.
      Update the inline comment referencing `allowNegativePOI`.
- [x] 1.6 `AxisBernsteinModel.compute`: rename argument `poi` → `param` and
      all internal uses.

## 2. Resolve `CompositeParamModel` conflicts

- [x] 2.1 In `CompositeParamModel.__init__`, take the upstream version:
      `self.allowNegativeParam = allowNegativeParam` and
      `self.is_linear = self.nparams == 0 or self.allowNegativeParam`.
- [x] 2.2 In `CompositeParamModel.compute`, take the upstream version:
      iterate over `self.param_models`, slice by `m.nparams`, pass raw
      `param` slice without per-sub-model squaring.

## 3. Fix `load_models` in `param_models/helpers.py`

- [x] 3.1 Add `AxisNormModel`, `AxisExpModel`, `AxisBernsteinModel` to
      `baseline_models` dict (value `"param_model"`).
- [x] 3.2 In `load_models`, derive the composite's `allowNegativeParam` as
      `any(m.allowNegativeParam for m in models)` rather than from
      `kwargs.get("allowNegativeParam", False)`, so axis models with
      `allowNegativeParam=True` are not broken when composited.

## 4. Delete `poi_models/helpers.py`

- [x] 4.1 Accept the upstream deletion of
      `rabbit/rabbit/poi_models/helpers.py` (`git rm`).

## 5. Resolve `parsing.py` conflict

- [x] 5.1 Take upstream: replace `--poiModel`/`default=[]` with
      `--paramModel`/`default=None` and update the help text.

## 6. Resolve `rabbit_fit.py` conflict

- [x] 6.1 Take upstream: use `args.paramModel or [["Mu"]]` and
      `ph.load_models(model_specs, indata, **vars(args))`.

## 7. Update `debug_poi_gradient.py`

- [x] 7.1 Replace `--poiModel` → `--paramModel` in `add_argument`.
- [x] 7.2 Replace `from rabbit.poi_models import helpers as ph` →
      `from rabbit.param_models import helpers as ph`.
- [x] 7.3 Replace `args.poiModel` → `args.paramModel`.
- [x] 7.4 Replace `poi_model.allowNegativePOI` → `poi_model.allowNegativeParam`
      and `poi_model.xpoidefault` → `poi_model.xparamdefault` in print statements.

## 8. Validation

- [x] 8.1 Syntax-check all modified Python files — PASSED (all 4 files OK).
- [x] 8.2 Confirm no remaining `allowNegativePOI`, `xpoidefault`, `set_poi_default`,
      `self\.pois`, `poiModel`, or conflict markers in rabbit source — PASSED (clean).
