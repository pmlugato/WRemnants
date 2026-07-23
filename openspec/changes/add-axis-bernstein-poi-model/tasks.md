## 1. AxisBernsteinModel class

- [x] 1.1 Add `AxisBernsteinModel` to `poi_model.py` with `parse_args` accepting
      exactly 4 positional args `(channel, proc_spec, shape_axis, cell_axes_csv)`
- [x] 1.2 In `__init__`: resolve channel, shape axis, cell axes, and processes
      from `indata.channel_info` with the same validation pattern as
      `AxisExpModel` (raise `ValueError` for unknown channel, axis, or process;
      raise `ValueError` if `shape_axis` appears in `cell_axes_csv`)
- [x] 1.3 Compute `n_cell`, `npoi = n_procs * 2 * n_cell`; generate POI names
      `c0_{proc}_{cell_label}` then `c1_{proc}_{cell_label}` for each process
- [x] 1.4 Normalize shape-axis bin centers to `[0, 1]` and store as a constant
      TF tensor `self.x_m`
- [x] 1.5 Build reshape helpers `cell_reshape`, `shape_reshape`, `full_shape`
      from channel axis ordering (same pattern as `AxisExpModel`)
- [x] 1.6 Set `allowNegativePOI=False`, `is_linear=False`, and call
      `set_poi_default(None, allowNegativePOI=False)` so `xpoidefault = ones`
      (x²=1 → c₀=c₁=1 → flat unit background at initialization)
- [x] 1.7 Implement `compute(poi, full=False)`: slice `c0_poi` and `c1_poi` per
      process, reshape each to `cell_reshape`, compute
      `scaling = c0*(1−x_reshaped) + c1*x_reshaped`, broadcast to `full_shape`,
      and accumulate into `irnorm` via the same scatter pattern as `AxisExpModel`

## 2. Registration

- [x] 2.1 Add `"AxisBernsteinModel": "poi_model"` to `helpers.py`

## 3. Validation

- [x] 3.1 `python -m py_compile rabbit/rabbit/poi_models/poi_model.py`
- [ ] 3.2 Confirm default POIs (c₀=c₁=1) produce flat unit rnorm in a smoke test
- [ ] 3.3 Confirm fit runs with
      `--poiModel AxisBernsteinModel <channel> bkgBernstein <mass_axis> <cell_axes>`
