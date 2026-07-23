# Change: Add axis-driven exponential background POI model and CLI composition

## Why

The existing btojpsik background model encodes `A·e^{−Bx}` shape via
template-morph NOIs, which are inherently linear perturbations around a seeded
shape. This makes `A` and `B` poorly identified, produces degenerate Hessians,
and couples shape and amplitude in a way that confounds postfit interpretation.
A proper parametric POI model with explicit `(A, B)` per (eta, pt, charge) cell
solves all three problems. Additionally, running two complementary POI models
(e.g. `AxisNormModel` for signal + flat background, `AxisExpModel` for the
exponential background) currently requires a custom wrapper class; a first-class
CLI composition mechanism makes this reproducible for all users.

## What Changes

- Add `AxisExpModel` to `rabbit/rabbit/poi_models/poi_model.py`. For each
  process in proc_spec and each bin of the cell axes, the model holds two
  independent POIs `(A, B)` and in `compute()` produces `rnorm = A·e^{−B·x_m}`
  where `x_m` is the normalized center of mass bin `m`. All other channels and
  processes are left at 1.0.
- CLI shape: `--poiModel AxisExpModel <channel> <proc_spec> <shape_axis> <cell_axes>`
  — always exactly 4 positional args.
- Register `AxisExpModel` in `rabbit/rabbit/poi_models/helpers.py`.
- Change `--poiModel` in `rabbit/rabbit/parsing.py` from `nargs="+"` to
  `action="append", nargs="+"` so the flag can be given multiple times.
  The fitter automatically wraps multiple models in `CompositePOIModel`.
- Update `scripts/rabbit/btojpsik_tensor.py`: add `--bkgModel exp_poi` option
  that writes all-ones templates for `flatBkg` and `bkgExp` processes (scaled
  by a per-cell yield estimate for initialization) and adds no NOIs for either.

## Impact

- Affected specs: `axis-exp-poi-model` (new), `poi-model-composition` (new)
- Affected code:
  - `rabbit/rabbit/poi_models/poi_model.py` (add `AxisExpModel`)
  - `rabbit/rabbit/poi_models/helpers.py` (register `AxisExpModel`)
  - `rabbit/rabbit/parsing.py` (`--poiModel` becomes repeatable)
  - `rabbit/bin/rabbit_fit.py` (compose multiple models if >1 given)
  - `scripts/rabbit/btojpsik_tensor.py` (add `--bkgModel exp_poi` path)
