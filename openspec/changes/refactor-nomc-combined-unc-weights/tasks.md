## 1. btojpsik.py — define concatenated combined kinematics
- [ ] 1.1 When `--nomc`, `--includeKaonScaleVariations`, and `diff_weights_helper` are
       all active: after the per-muon response-weight columns are defined, concatenate
       the six kinematic RVecs (recoPt, recoEta, recoCharge, genPt, genEta, genCharge)
       and the response-weight RVec from kaon + mu1 + mu2 into single 3-element RVecs
       under the prefix `btojpsik_combined_sel`
- [ ] 1.2 Call `add_jpsi_crctn_stats_unc_hists` once with `reco_sel_GF="btojpsik_combined_sel"`
       and `hist_suffix=""` — writes `nominal_muonScaleSyst_responseWeights` with the
       combined product tensor

## 2. btojpsik.py — gate existing per-particle calls
- [ ] 2.1 Remove the existing per-particle `add_jpsi_crctn_stats_unc_hists` calls
       (kaon and mu1/mu2) from the nomc + includeKaonScaleVariations code path
- [ ] 2.2 When `--validationHists`: write per-particle HistoBoosts by calling
       `add_jpsi_crctn_stats_unc_hists` separately for kaon (`hist_suffix="_kaon"`),
       mu1 (`hist_suffix="_mu1"`), and mu2 (`hist_suffix="_mu2"`)

## 3. btojpsik.py — non-nomc kaon path
- [ ] 3.1 Confirm existing kaon call in the non-nomc `includeKaonScaleVariations`
       branch is unchanged — no code change required

## 4. btojpsik_tensor.py
- [ ] 4.1 Remove the `variation_hist_mu1` / `variation_hist_mu2` loading block in the
       `if args.nomc:` branch
- [ ] 4.2 Remove the two `_add_aem_systematics(... name_suffix="_mu1")` and
       `_add_aem_systematics(... name_suffix="_mu2")` calls

## 5. muon_calibration.py
- [ ] 5.1 No interface changes required — verify by inspection
