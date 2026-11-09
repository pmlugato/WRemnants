"""Rename map that correlates the Z and J/psi muon-scale nuisances.

Consumed by ``rabbit/bin/merge_inputs.py --rename``, which loads this file by
path and looks for ``rename_systs``, ``rename_groups``, ``rename_processes``,
and ``rename_channels``.

Correlation in ``merge_inputs.py`` happens by name and by nothing else: two
systematics that carry the same name in different channels become one fit
parameter. So this file is the entire mechanism by which the J/psi channels
constrain the same A/e/M that shift the Z lineshape.

The index convention is shared because both histmakers build their variations
from ``muon_calibration.make_jpsi_crctn_helpers`` with the same arguments. Under
``--fitMuonScaleAndResolution`` that helper takes its no-covariance branch,
where the packed ``[eta x scale_params x unc]`` histogram is filled at

    ivar = ieta * n_scale_params + iparm,   iparm in (A, e, M) order

``lowmass_tensor.py`` reads the same axis with the same arithmetic and names the
result ``A{i}`` / ``e{i}`` / ``M{i}``; ``setupRabbit.py`` keeps the raw index and
emits ``Scale_correction_unc{N}``. Both names were read off real tensors rather
than guessed.

Verified on the probe tensors:

* the calibration file has 48 eta bins, so there are 3 x 48 = 144 nuisances
* the Z tensor carries ``Scale_correction_unc0`` .. ``Scale_correction_unc143``
* every ``unc`` index touches exactly one ``(eta, param)`` cell, checked for all
  144

The joint fits group the basis into 24 eta bins of 0.2 (both histmakers run
with ``--etaBins 24``; the B+ histmaker with ``--scaleVariationEtaBins 24``),
because at 48 adjacent nuisance pairs come out nearly degenerate. The map does
not need to know that number; see ``N_ETA_MAX``.
"""

# Size of the calibration file's native eta basis. The rename is written for
# all of it, because the index arithmetic below does not depend on how coarsely
# the histmakers grouped the basis (their --etaBins / --scaleVariationEtaBins):
# a Z tensor built on N eta groups carries Scale_correction_unc0 .. unc{3N-1},
# exactly the first 3N keys of this map, and the remaining keys match nothing.
# Which N the merge actually used is verified by check_joint_tensor.py --nEta,
# against what the merged tensor contains.
N_ETA_MAX = 48

# Order matters: it is the order the packed histogram is filled in.
PARAMS = ["A", "e", "M"]

# Base name setupRabbit gives the Z-side scale nuisances, from
# Datagroups.systHists building `baseName + axisLabel + entry`.
Z_SCALE_BASE = "Scale_correction_unc"

rename_systs = {
    f"{Z_SCALE_BASE}{ieta * len(PARAMS) + iparm}": f"{param}{ieta}"
    for ieta in range(N_ETA_MAX)
    for iparm, param in enumerate(PARAMS)
}

# Deliberately EMPTY, matching havyn's rename.py.
#
# An earlier version here mapped `scaleCrctn` (the Z side's single group) and
# `A`/`e`/`M` (the J/psi side's three) all onto one `scale_AeM`. That produced a
# grouped-impact table with a single calibration row, which is why our impacts
# chart could not be compared against havyn's -- theirs splits Magnetic Field
# (A), Energy Loss (e) and Misalignment (M) as separate rows. `merge_inputs.py`
# unions groups per systematic name, so leaving this empty gives both: the three
# J/psi groups AND the Z-side `scaleCrctn` covering the same parameters.
rename_groups = {}

# Cosmetic, and matching havyn's rename.py exactly.
rename_channels = {
    "ch0": "Z",
}
