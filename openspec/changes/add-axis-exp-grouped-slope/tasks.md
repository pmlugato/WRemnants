## 1. AxisExpModel constructor

- [x] 1.1 Add `slope_axes_csv` parameter (default `None`) to `AxisExpModel.__init__`; when `None` default to `cell_axes_csv`
- [x] 1.2 Parse and validate slope axis names: must be non-empty subset of cell axis names; raise `ValueError` naming any axis not in cell axes
- [x] 1.3 Compute `self.slope_axes`, `self.slope_axis_names`, `self.n_slope_groups`
- [x] 1.4 Compute `self.slope_cell_reshape`: same as `cell_reshape` but with 1 for every non-slope-axis dimension
- [x] 1.5 Update `self.npoi = len(self.proc_idxs) * (self.n_cell + self.n_slope_groups)`
- [x] 1.6 Update POI name generation: lnAmpl names keyed by cell axes (unchanged); slope names keyed by slope axes only
- [x] 1.7 Update `self.xpoidefault = tf.zeros([self.npoi], ...)` (size already correct via updated npoi)

## 2. AxisExpModel parse_args

- [x] 2.1 Update `parse_args` to accept 4 or 5 positional args; pass 5th as `slope_axes_csv` when present; raise `ValueError` otherwise

## 3. AxisExpModel compute()

- [x] 3.1 Update POI slicing: lnAmpl slice stays `[0..n_cell)` per process; slope slice is `[n_cell..n_cell+n_slope_groups)` per process
- [x] 3.2 Reshape `b_poi` with `self.slope_cell_reshape` instead of `self.cell_reshape`
- [x] 3.3 Verify `tf.exp(a + b * x_reshaped)` broadcasts correctly to `self.full_shape`

## 4. Validation

- [x] 4.1 `python -m py_compile rabbit/rabbit/poi_models/poi_model.py`
- [x] 4.2 Confirm backward compatibility: 4-arg invocation produces `n_slope_groups == n_cell` and identical layout to current code
- [ ] 4.3 Confirm btojpsik fit command with `bkmm_kaon_eta` as 5th arg runs without error and produces positive-definite Hessian
