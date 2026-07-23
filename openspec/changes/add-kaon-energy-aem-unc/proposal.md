# Change: Use kaon energy in the e-term of A/e/M uncertainty weights

## Why

The `e` parameter in the A/e/M calibration is an energy-loss-in-material correction.
The physically correct curvature shift from an energy-loss variation scales as
dκ ∝ (E/p) · k², where E is the particle's total energy and p its total momentum.
For muons (m_μ ≪ p_μ) the ratio E/p ≈ 1 and the massless approximation used today is
accurate. For kaons with m_K ≈ 494 MeV and pT ~ 1–8 GeV, E_K/p_K can differ from 1 by
~10–30%, so the current `eUnc * k²` form under-estimates the kaon energy-loss uncertainty.
The A and M terms are unaffected (A is a scale offset, M is magnetic-field related), so
the change is strictly scoped to the e term.

## What Changes

- **C++ (`muon_calibration.hpp`)**: Add an optional `particle_mass` parameter (default
  `0.0`) to `calculateQopUnc(pt, eta, charge, AUnc, eUnc, MUnc)`. When `particle_mass > 0`,
  compute `p_total = pt * cosh(eta)`, `E_over_p = sqrt(1 + mass²/p_total²)`, and multiply
  the `eUnc * k` sub-expression by `E_over_p` before forming `kUnc`. A and M terms are
  unchanged. The massless limit (`mass = 0`) reproduces the current result exactly.
- **C++ (`muon_calibration.hpp`)**: Propagate `particle_mass` through
  `JpsiCorrectionsUncHelperSplines` (stored at construction, forwarded to every
  `calculateQopUnc` call inside the operator).
- **Python (`wremnants/production/muon_calibration.py`)**: Add a `particle_mass` kwarg
  (default `0.0`) to `make_jpsi_crctn_unc_helper`. Thread it through to the helper
  construction so it is stored and forwarded to the C++ class.
- **Histmaker (`scripts/histmakers/btojpsik.py`)**: Pass `particle_mass=KAON_MASS_GEV`
  (≈ 0.49368) when building `jpsi_crctn_data_unc_helper` and `jpsi_crctn_MC_unc_helper`
  for the kaon path. The muon uncertainty helpers (`muon_diff_weights_helper` and its
  associated unc helper) keep the default `particle_mass=0`.

## Impact

- Affected specs: `btojpsik-uncertainties` (new capability spec, no existing spec)
- Affected code:
  - `wremnants/production/include/muon_calibration.hpp` — `calculateQopUnc` and
    `JpsiCorrectionsUncHelperSplines`
  - `wremnants/production/muon_calibration.py` — `make_jpsi_crctn_unc_helper`
  - `scripts/histmakers/btojpsik.py` — `make_jpsi_crctn_helpers` call site
- Not breaking: default `particle_mass = 0` preserves the massless result for all
  existing callers (muon scale variation, W/Z histmakers, etc.)
