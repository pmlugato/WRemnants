# Change: Add axis-driven normalization POI model to rabbit

## Why

`wremnants/postprocessing/btojpsik_poi_models.py` contains `SignalNormModel`,
which assigns one unconstrained normalization POI per cell of the non-mass axes
(pt × eta × charge) for the signal process. The channel, axis set, and process
are all hardcoded to the btojpsik calibration context. Any other channel or
process that needs the same per-bin-cell normalization structure must duplicate
this logic. Moving a configurable version into rabbit makes it a first-class
reusable primitive available to all analyses.

## What Changes

- Add `AxisNormModel` to `rabbit/rabbit/poi_models/poi_model.py` alongside the
  existing built-in models (`Ones`, `Mu`, `Mixture`, `SaturatedProjectModel`).
  The model assigns one normalization POI per bin-combination of a
  caller-specified list of axes, applied to a caller-specified set of processes
  within a named channel.
- CLI shape: `--poiModel AxisNormModel <channel> <proc_spec> <axes>`
  — always exactly 3 positional args. `proc_spec` is `all` or a
  comma-separated list of process names. `axes` is a comma-separated list of
  axis names.
- Register `AxisNormModel` in `rabbit/rabbit/poi_models/helpers.py` so the
  built-in loader resolves it by class name.
- Update `scripts/rabbit/btojpsik_tensor.py`: no change needed for
  `signal=True` since btojpsik will pass the process name explicitly.
- Remove `wremnants/postprocessing/btojpsik_poi_models.py::SignalNormModel`
  once callers have migrated to `AxisNormModel`.

## Impact

- Affected specs: `axis-norm-poi-model` (new)
- Affected code:
  - `rabbit/rabbit/poi_models/poi_model.py` (add AxisNormModel class)
  - `rabbit/rabbit/poi_models/helpers.py` (register new model)
  - `wremnants/postprocessing/btojpsik_poi_models.py` (remove SignalNormModel)
