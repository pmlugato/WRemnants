# Design: AxisExpModel and CLI POI Composition

## Context

The btojpsik mass fit has three physical components per (eta, pt, charge) cell:
1. Signal peak — per-cell normalization, handled by `AxisNormModel`.
2. Flat combinatorial background — per-cell normalization, handled by `AxisNormModel`.
3. Exponential falling background — mass-dependent shape `A·e^{−Bx}`, not
   representable as a uniform per-cell scaling.

`AxisNormModel` broadcasts a single scale factor across all mass bins in a cell.
The exponential background requires the scaling to vary **within** the mass axis,
requiring a different parametric form.

## Goals / Non-Goals

- **Goal**: `AxisExpModel` with explicit `(A, B)` POIs per (process, cell),
  mass-axis-dependent compute, and no tensor-level shape seeding.
- **Goal**: first-class CLI composition via repeated `--poiModel` flags.
- **Goal**: `btojpsik_tensor.py` `--bkgModel exp_poi` writes all-ones templates
  so the POI model owns the full background shape.
- **Non-goal**: shared `B` across cells (may be added later).
- **Non-goal**: higher-order shapes (polynomial, etc.) — `AxisExpModel` is
  restricted to the `A·e^{−Bx}` form.
- **Non-goal**: replacing or modifying `AxisNormModel`.

## Decisions

### AxisExpModel CLI shape

**Decision**: 4 positional args:
```
--poiModel AxisExpModel <channel> <proc_spec> <shape_axis> <cell_axes>
```

- `shape_axis`: the single axis along which `e^{−Bx}` is evaluated (e.g.
  `bkmm_jpsimc_mass`). Bin centers are read from `channel_info` and normalized
  to `[0, 1]` so `B` is dimensionless.
- `cell_axes`: comma-separated axes where `(A, B)` vary independently (e.g.
  `bkmm_kaon_pt,bkmm_kaon_eta,bkmm_kaon_charge`).

`shape_axis` and `cell_axes` must be disjoint. Together they need not cover all
channel axes — any unspecified axis is broadcast over (same POI for every bin
of that axis).

### POI layout per process

For each process in proc_spec:
- `n_cell = prod(cell_axis_sizes)` cells
- `n_poi_per_proc = 2 * n_cell` (one A and one B per cell)
- POI names: `expA_{procname}_{cell_label}`, `expB_{procname}_{cell_label}`
- Total: `npoi = n_procs * 2 * n_cell`

A and B for different processes are fully independent (same principle as
`AxisNormModel`).

### is_linear

`False`. The likelihood gradient with respect to B involves `x·A·e^{−Bx}`,
which is nonlinear. The SciPy minimizer handles this via automatic
differentiation through TensorFlow; no special treatment is needed beyond
setting `is_linear = False`.

### allowNegativePOI

Both A and B default to `allowNegativePOI = False`, parameterized as `x²`
internally. This keeps A positive (physical background amplitude) and B
positive (falling exponential). A rising exponential (`B < 0`) would be
unusual for combinatorial background and is excluded by default.

### x coordinate normalization

Bin centers of `shape_axis` are normalized to `[0, 1]`:
```
x_m = (center_m - center_0) / (center_{N-1} - center_0)
```
This makes B dimensionless and of order 1 for typical falling backgrounds,
which helps the optimizer. The normalized centers are stored as a constant
TF tensor at construction time.

### Default initialization

Default POI: `A = 1.0`, `B = 1.0` (before squaring for the `x²`
reparametrization). With an all-ones template, the initial predicted background
yield per mass bin is `e^{−x_m}`, i.e. roughly 1 event/bin at x=0 and
`e^{−1} ≈ 0.37` at x=1. The btojpsik tensor writer scales the all-ones
template by a rough per-cell yield estimate (see below) so the POI starts near
its true value.

### CLI composition mechanism

**Decision**: change `--poiModel` in `parsing.py` to `action="append",
nargs="+"`. Each occurrence of `--poiModel` appends one model spec (class name
+ positional args) as a sublist. `args.poiModel` becomes a list of lists.

In `rabbit_fit.py`:
```python
models = [ph.load_model(spec[0], indata, *spec[1:], **vars(args))
          for spec in args.poiModel]
poi_model = models[0] if len(models) == 1 else poi_model_module.CompositePOIModel(models)
```

Default changes from `["Mu"]` to `[["Mu"]]` to match the new structure.
Single-model invocations are fully backward compatible.

**Example btojpsik composition**:
```bash
--poiModel AxisNormModel btojpsik_stuff signal,flatBkg \
    bkmm_kaon_pt,bkmm_kaon_eta,bkmm_kaon_charge \
--poiModel AxisExpModel btojpsik_stuff bkgExp \
    bkmm_jpsimc_mass \
    bkmm_kaon_pt,bkmm_kaon_eta,bkmm_kaon_charge
```

### btojpsik tensor writer: `--bkgModel exp_poi`

New choice alongside existing `flat` and `exp`.

When `--bkgModel exp_poi` is set (implies `--signalNormPOI`):
- Adds `flatBkg` process: all-ones histogram scaled by per-cell flat yield
  estimate `max(data_yield - signal_yield, 0) / n_mass_bins`.
- Adds `bkgExp` process: all-ones histogram scaled by per-cell exp yield
  estimate (same formula, or a fixed fraction of the flat estimate).
- No NOIs added for either process.
- Signal is added with `signal=True`.

The yield scaling means both A and B POIs start near 1.0, so the optimizer
does not need to travel far from the default initialization.

The old `--bkgModel exp` (NOI-based) path is retained for backward
compatibility but not recommended.

## Risks / Trade-offs

- `is_linear = False` means the fitter cannot use the linear approximation for
  `CompositePOIModel.is_linear`. The combined model will be nonlinear even if
  `AxisNormModel` alone would be linear. This is unavoidable given the physics.
- `B > 0` restriction (via `x²`) prevents a rising background. If future use
  cases require this, `allowNegativePOI=True` can be passed.
- The yield scaling heuristic in the tensor writer is a rough estimate. If the
  signal MC is a poor description of the data before fitting, the estimated flat
  background yield could be negative (clamped to 0). In that case, the initial
  background yield is 0, which may slow convergence but does not break the fit.

## Open Questions

- None.
