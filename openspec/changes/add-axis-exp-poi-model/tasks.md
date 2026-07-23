## 1. AxisExpModel

- [ ] 1.1 Add `AxisExpModel` to `rabbit/rabbit/poi_models/poi_model.py`
       (constructor, `parse_args` requiring exactly 4 args, `compute`)
- [ ] 1.2 Register `AxisExpModel` in `rabbit/rabbit/poi_models/helpers.py`

## 2. CLI composition

- [ ] 2.1 Change `--poiModel` in `rabbit/rabbit/parsing.py` to
       `action="append", nargs="+"` and update default to `[["Mu"]]`
- [ ] 2.2 Update `rabbit/bin/rabbit_fit.py` to instantiate one model per
       sublist and wrap in `CompositePOIModel` when `len > 1`

## 3. btojpsik tensor writer

- [ ] 3.1 Add `exp_poi` to `--bkgModel` choices in
       `scripts/rabbit/btojpsik_tensor.py`
- [ ] 3.2 Implement the `exp_poi` path: all-ones templates scaled by per-cell
       yield estimates for `flatBkg` and `bkgExp`, no NOIs added

## 4. Validation

- [ ] 4.1 `python -m py_compile rabbit/rabbit/poi_models/poi_model.py`
- [ ] 4.2 Verify `AxisExpModel` POI count: 1 proc × 2 × pt(8) × eta(7) ×
       charge(2) = 224
- [ ] 4.3 Verify `compute()`: A=1, B=0 → rnorm=1.0 everywhere for target
       process; A=1, B>0 → rnorm decreases monotonically across mass bins
- [ ] 4.4 Verify two-model composition: `--poiModel AxisNormModel ... --poiModel
       AxisExpModel ...` produces a `CompositePOIModel` with correct total npoi
- [ ] 4.5 Verify single-model invocation is unchanged (no wrapping, no
       regression on existing fits)
- [ ] 4.6 Verify `--bkgModel exp_poi` tensor has no `norm_flatBkg_*` or
       `norm_bkgExp_*` NOIs
