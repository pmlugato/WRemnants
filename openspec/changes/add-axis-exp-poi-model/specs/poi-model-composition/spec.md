## ADDED Requirements

### Requirement: Repeatable --poiModel CLI flag
The system SHALL allow `--poiModel` to be specified multiple times on the
`rabbit_fit.py` command line. Each occurrence provides one complete model
specification (class name followed by its positional arguments). When multiple
models are given, the fitter automatically wraps them in `CompositePOIModel`.
A single `--poiModel` invocation behaves identically to the current behavior.

The implementation MUST:
- Change `--poiModel` in `rabbit/rabbit/parsing.py` from `nargs="+"` to
  `action="append", nargs="+"` so `args.poiModel` becomes a list of lists.
- Update the default from `["Mu"]` to `[["Mu"]]`.
- In `rabbit/bin/rabbit_fit.py`, instantiate one model per sublist and compose
  them with `CompositePOIModel` when more than one is present.

#### Scenario: Single --poiModel unchanged
- **WHEN** the command line contains exactly one `--poiModel AxisNormModel ch p axes`
- **THEN** `args.poiModel == [["AxisNormModel", "ch", "p", "axes"]]` and the
  fit uses a plain `AxisNormModel` (no `CompositePOIModel` wrapping)

#### Scenario: Two --poiModel flags compose automatically
- **WHEN** the command line contains
  `--poiModel AxisNormModel ch signal,flatBkg axes1`
  `--poiModel AxisExpModel ch bkgExp mass axes1`
- **THEN** the fit uses a `CompositePOIModel` wrapping both models, with
  `npoi == AxisNormModel.npoi + AxisExpModel.npoi`

#### Scenario: Default Mu model unchanged
- **WHEN** `--poiModel` is not provided
- **THEN** the fit uses the default `Mu` model as before

#### Scenario: Three or more models compose correctly
- **WHEN** three `--poiModel` flags are given
- **THEN** all three are composed into a single `CompositePOIModel` and the
  combined `pois` array is the concatenation of all three models' POI names
