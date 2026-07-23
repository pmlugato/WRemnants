# Change: Add --histToFit m32 observable (m(B) − m(J/ψ))

## Why

The nomc B mass peak is broader than the jpsimc peak because the J/ψ mass constraint is
absent. The Q-value observable Q = m(B, nomc) − m(J/ψ, nomc) recovers peak sharpness:
both masses share the same muon 4-momenta, so the leading muon momentum-scale dependence
partially cancels in the difference, leaving a peak dominated by kaon momentum resolution.
This produces a more numerically stable fit observable while preserving kaon A/e/M
sensitivity. Muon A/e/M sensitivity is intentionally reduced by the cancellation — this
trade-off is acceptable when fit stability is the priority.

## What Changes

- Add `--histToFit` argument to `btojpsik.py` with `choices=["mB", "m32"]` and
  `default="mB"`. `m32` selects Q = m(B, nomc) − m(J/ψ, nomc) as the fit observable.
- When `--histToFit m32`: Define `bkmm_{vertex_prefix}_m32_scalar` in `build_graph` as
  `bkmm_{vertex_prefix}_mass_scalar − mm_kin_mass`. After `select_only_passing_bkmm_candidates`,
  `mm_*` columns are reindexed to scalars and `bkmm_*` columns to length-1 RVecs, so
  `mm_kin_mass` is already the scalar dimuon mass for the selected candidate.
- Update `configure_fit_histogram` to accept `hist_to_fit` and select the appropriate
  column name and axis.
- Add `bkmm_nomc_m32` and `bkmm_jpsimc_m32` axes to `btojpsik_axes.py`:
  `Regular(100, 1.8, 2.6, name=..., underflow=False, overflow=False)`.
- No changes to selection functions, tensor writer, or fit scripts. The tensor writer
  receives `--massAxis bkmm_nomc_m32` on the command line.

## Axis Range

Signal Q ≈ m(B) − m(J/ψ) ≈ 5.279 − 3.097 = 2.182 GeV. With the B mass window [5.2, 5.4]
and the dimuon kinematic-fit mass peaked near m(J/ψ), the signal peak sits comfortably
within [1.8, 2.6] GeV (100 bins, 8 MeV/bin).

## Backward Compatibility

`--histToFit mB` (the default) reproduces current behaviour exactly. The m32 column
Define and axis lookup are gated on `args.histToFit == "m32"`.

| Modified interface | New parameter | Default | Effect of default |
|-|-|-|-|
| `btojpsik.py` | `--histToFit` | `"mB"` | unchanged column/axis |
| `configure_fit_histogram` | `hist_to_fit` | `"mB"` | unchanged |

## Physical Notes

- `mm_kin_mass` is the unconstrained dimuon kinematic-fit mass — genuinely "nomc" for
  the dimuon subsystem. Combining it with `bkmm_nomc_mass` is self-consistent.
- Combining `bkmm_jpsimc_mass` with `mm_kin_mass` would be a physically odd hybrid
  (one mass with constraint, one without). The `m32` mode is meaningful in any
  `--histToFit` / `--nomc` combination but is primarily intended for `--nomc`.
- Kaon A/e/M NOIs remain fully sensitive (kaon enters only m(B)). Muon A/e/M
  sensitivity is reduced due to the partial cancellation in the difference.

## Impact

- Affected specs: histmaker-btojpsik (delta)
- Affected code: `scripts/histmakers/btojpsik.py`,
  `wremnants/production/btojpsik_axes.py`
