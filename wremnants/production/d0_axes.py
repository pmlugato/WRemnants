"""Template-axis definitions shared by the D0 histmaker and its tests."""

import hist

DEFAULT_MRK_EDGES = (
    0.00,
    0.06,
    0.10,
    0.14,
    0.18,
    0.22,
    0.26,
    0.30,
    0.34,
    0.38,
    0.42,
    0.46,
    0.50,
    0.60,
    0.70,
    0.85,
    1.05,
    1.30,
    1.55,
    1.80,
)


def make_d0_template_axes(eta_bins=24, mrk_bins=None):
    """Return the etaK, mRK, and D0-mass axes used by calibration templates."""
    if eta_bins <= 0:
        raise ValueError("eta_bins must be positive")
    if mrk_bins is not None and mrk_bins <= 0:
        raise ValueError("mrk_bins must be positive when specified")

    eta_axis = hist.axis.Regular(eta_bins, -2.4, 2.4, name="etaK")
    if mrk_bins is None:
        mrk_axis = hist.axis.Variable(DEFAULT_MRK_EDGES, name="mRK")
    else:
        mrk_axis = hist.axis.Regular(
            mrk_bins, DEFAULT_MRK_EDGES[0], DEFAULT_MRK_EDGES[-1], name="mRK"
        )
    mass_axis = hist.axis.Regular(25, 1.8, 1.93, name="D0mass")
    return eta_axis, mrk_axis, mass_axis
