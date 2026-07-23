## 1. C++ — mass-aware e-term in calculateQopUnc

- [ ] 1.1 Add optional `double particle_mass = 0.0` argument to the
      `calculateQopUnc(pt, eta, charge, AUnc, eUnc, MUnc)` overload in `muon_calibration.hpp`
- [ ] 1.2 Inside that overload, compute:
      `float p_total = pt * std::cosh(eta);`
      `double E_over_p = (particle_mass > 0.0) ? std::sqrt(1.0 + particle_mass*particle_mass / (p_total*p_total)) : 1.0;`
      and replace `eUnc * k` with `eUnc * k * E_over_p` in the `kUnc` expression
- [ ] 1.3 Verify: `particle_mass = 0` → `E_over_p = 1` → identical output to pre-change

## 2. C++ — propagate mass through JpsiCorrectionsUncHelperSplines

- [ ] 2.1 Add a `double particle_mass_` member to `JpsiCorrectionsUncHelperSplines`
- [ ] 2.2 Add a constructor argument `double particle_mass = 0.0` and store it
- [ ] 2.3 Forward `particle_mass_` to every `calculateQopUnc(recPt, recEta, recCharge, AUnc, eUnc, MUnc)` call inside the operator

## 3. Python — thread mass through make_jpsi_crctn_unc_helper

- [ ] 3.1 Add `particle_mass: float = 0.0` kwarg to `make_jpsi_crctn_unc_helper`
- [ ] 3.2 Pass `particle_mass` to the `JpsiCorrectionsUncHelperSplines` constructor at
      the call site where the C++ helper is instantiated
- [ ] 3.3 Propagate the kwarg through `make_jpsi_crctn_helpers` → `make_jpsi_crctn_unc_helper`

## 4. Histmaker wiring

- [ ] 4.1 Define `KAON_MASS_GEV = 0.49368` near the top of `btojpsik.py`
- [ ] 4.2 Pass `particle_mass=KAON_MASS_GEV` in the `make_jpsi_crctn_helpers` call for
      the kaon uncertainty helper (`jpsi_crctn_data_unc_helper` / `jpsi_crctn_MC_unc_helper`)
- [ ] 4.3 Confirm the muon uncertainty helper path (nomc) keeps `particle_mass=0` (default)

## 5. Validation

- [ ] 5.1 `python -m py_compile scripts/histmakers/btojpsik.py wremnants/production/muon_calibration.py`
- [ ] 5.2 Run with `--maxFiles 1 --includeKaonScaleVariations` and confirm
      `nominal_muonScaleSyst_responseWeights` is present in the output HDF5
- [ ] 5.3 Compare e-variation alternative weights before and after: expect a nonzero
      shift, largest in low-pT kaon bins where m_K/p_K is largest
- [ ] 5.4 Confirm A and M alternative weights are unchanged (regression check)
- [ ] 5.5 Confirm nomc muon-only path (muon unc helpers with mass=0) produces identical
      output to pre-change
