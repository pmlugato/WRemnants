# Change: Migrate axis POI models to upstream ParamModel framework

## Why

Upstream refactored `rabbit/poi_models/poi_model.py` →
`rabbit/param_models/param_model.py`, renaming the base class (`POIModel` →
`ParamModel`) and all its attributes (`allowNegativePOI` → `allowNegativeParam`,
`xpoidefault` → `xparamdefault`, `pois` → `params`, `set_poi_default` →
`set_param_default`, and adding the `npou` concept for model-internal nuisances).
The three axis models (`AxisNormModel`, `AxisExpModel`, `AxisBernsteinModel`) and
the local helper (`poi_models/helpers.py`) still use the old names, leaving four
merge conflicts and one modify/delete conflict. This change resolves all five
without altering any fitted behaviour.

## What Changes

- **`param_model.py`** — rename old attributes in the three axis model classes:
  `self.pois` → `self.params`, `allowNegativePOI` → `allowNegativeParam`,
  `xpoidefault` → `xparamdefault`, `set_poi_default` → `set_param_default`,
  `compute(self, poi, ...)` → `compute(self, param, ...)`, add `self.npou = 0`.
- **`CompositeParamModel`** — resolve the two conflict regions by taking the
  upstream version (no per-sub-model squaring in `compute()`; composite
  `allowNegativeParam` set from constructor; slices by `m.nparams`).
- **`param_models/helpers.py`** — add `AxisNormModel`, `AxisExpModel`,
  `AxisBernsteinModel` entries pointing at `"param_model"`. Fix `load_models` to
  derive the composite's `allowNegativeParam` from sub-models (True if any
  sub-model has `allowNegativeParam=True`), not solely from CLI `--allowNegativeParam`.
  Without this fix, compositing axis models whose `allowNegativeParam=True` with
  the upstream default of `False` would cause the fitter to incorrectly square
  their raw parameters before passing them to `compute()`.
- **`poi_models/helpers.py`** — accept the upstream deletion (resolve
  modify/delete conflict by removing this file).
- **`parsing.py`** — take upstream: rename `--poiModel` → `--paramModel`,
  change `default=[]` → `default=None`.
- **`rabbit_fit.py`** — take upstream: use `ph.load_models` and
  `args.paramModel or [["Mu"]]`.
- **`scripts/rabbit/debug_poi_gradient.py`** — update all old-API references:
  `--poiModel` → `--paramModel`, `poi_models` import → `param_models`,
  `args.poiModel` → `args.paramModel`, `allowNegativePOI` → `allowNegativeParam`,
  `xpoidefault` → `xparamdefault`.

## Impact

- Affected specs: `axis-exp-poi-model`, `axis-norm-poi-model`,
  `axis-bernstein-poi-model`
- Affected code: `rabbit/rabbit/param_models/param_model.py`,
  `rabbit/rabbit/param_models/helpers.py`, `rabbit/rabbit/parsing.py`,
  `rabbit/bin/rabbit_fit.py`, `rabbit/rabbit/poi_models/helpers.py` (deleted),
  `scripts/rabbit/debug_poi_gradient.py`
- No change to fitted behaviour, POI names, or tensor format
- No changes to tensor writer, histmaker, or any other scripts
