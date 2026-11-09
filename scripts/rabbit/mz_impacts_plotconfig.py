"""Plot config for rabbit_plot_pulls_and_impacts.py: the seven-row grouped impact
plot on mZ used to compare channel combinations (Z + J/psi, Z + B+, Z + J/psi + B+).

The Z and J/psi tensors call the scale groups A/e/M, the B+ tensor calls them
scaleA/scalee/scaleM; the full fit carries both names with identical content.
Pick the grouping whose names exist in BOTH files of an overlay:
  mzscale_aem   for full vs Z + J/psi (no B+)
  mzscale_scale for full vs Z + B+ (no J/psi)
"""

_rest = ["stat", "binByBinStat", "theory"]
_rest_mb = ["resolution", "stat", "binByBinStat"]  # J/psi + B+ has no Z theory group
nuisance_grouping = {
    "mzscale_aem": ["Total", "A", "e", "M"] + _rest,
    "mzscale_scale": ["Total", "scaleA", "scalee", "scaleM"] + _rest,
    "mb_aem": ["Total", "A", "e", "M"] + _rest_mb,
}
impact_labels = {
    "resolution": "B⁺ resolution (24 free widths)",
    "A": "Magnetic Field (A)",
    "e": "Energy Loss (e)",
    "M": "Misalignment (M)",
    "scaleA": "Magnetic Field (A)",
    "scalee": "Energy Loss (e)",
    "scaleM": "Misalignment (M)",
    "stat": "Data stat.",
    "binByBinStat": "Simulation stat.",
    "theory": "Theory",
}
# rabbit_plot_cov reads this one for the axis tick labels
systematics_labels = {
    "massShiftZ100MeV": "mZ",
    "massShiftB5MeV": "mB",
    "massShiftB10MeV": "mB",
}
