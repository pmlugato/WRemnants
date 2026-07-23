## 1. Implementation

- [ ] 1.1 Add `AxisNormModel` to `rabbit/rabbit/poi_models/poi_model.py`
       (constructor, `parse_args` requiring exactly 3 args, `compute`)
- [ ] 1.2 Register `AxisNormModel` in `rabbit/rabbit/poi_models/helpers.py`
- [ ] 1.3 Remove `SignalNormModel` from
       `wremnants/postprocessing/btojpsik_poi_models.py`; update any callers
       that reference it via `--poiModel`

## 2. Validation

- [ ] 2.1 `python -m py_compile rabbit/rabbit/poi_models/poi_model.py`
- [ ] 2.2 Verify POI count/names: channel with axes `[mass(10), pt(4), eta(7), charge(2)]`,
       proc `BuToJpsiK_2018`, axes `pt,charge` → `npoi == 8`
- [ ] 2.3 Verify `compute()` shapes and per-process scaling: target process
       column equals POI value; other columns equal 1.0
- [ ] 2.4 Verify `proc_spec=all` scales all process columns
- [ ] 2.5 Verify wrong arg count, unknown channel, unknown axis, and unknown
       process each raise `ValueError` at construction time
- [ ] 2.6 `grep -r SignalNormModel` returns no live references after removal
