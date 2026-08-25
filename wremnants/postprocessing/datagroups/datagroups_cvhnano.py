"""Datagroups for the CVH NanoAOD B+ -> J/psi K analysis.

There is one simulated sample -- an inclusive b -> J/psi X production -- so the
stack cannot be built by splitting datasets the way `datagroups_btojpsik` does
for the old exclusive BuToJpsiK / BuToJpsiPi trio. It is built instead by
slicing the `genCategory` axis that `btojpsik_cvhnano.py` puts on every
histogram, with the same dataset appearing as a member of several groups under
several different `memberOp` slices.

`genCategory` encoding is defined once, in
`wremnants.production.btojpsik_cvhnano_axes`; this module imports it rather than
repeating the integers, so the two cannot drift apart.

Every group here is a **real b-hadron decay**, identified from the per-leg
generator parentage. The combinatorial bin is deliberately not among them; see
`NoBAncestor` below.
"""

from wremnants.production.btojpsik_cvhnano_axes import (
    GEN_BD,
    GEN_BS,
    GEN_BU_KX,
    GEN_BU_PIX,
    GEN_DATA,
    GEN_NO_B_ANCESTOR,
    GEN_OTHER_B,
    GEN_SIGNAL,
)
from wums import logging

logger = logging.child_logger(__name__)

GEN_CATEGORY_AXIS = "genCategory"

# Group name -> genCategory bin, in stack order (signal last, so it sits on top
# of the small backgrounds rather than under them).
MC_GROUPS = [
    ("OtherB", GEN_OTHER_B),
    ("Bs", GEN_BS),
    ("Bd", GEN_BD),
    ("BuJpsiPiX", GEN_BU_PIX),
    ("BuJpsiKX", GEN_BU_KX),
    ("BuJpsiK", GEN_SIGNAL),
]

# `NoBAncestor` is combinatorial: candidates whose three legs have no b-hadron
# common ancestor. It is a **diagnostic group, not a template**, and callers are
# expected to exclude it -- `btojpsik_cvhnano_plots.sh` passes
# `--excludeProcs NoBAncestor` for every stack except the one plot whose purpose
# is to show what the residual looks like.
#
# Why it must not be fitted: the simulation is an inclusive b -> J/psi X
# production, so it contains **none** of the prompt-charmonium or fake-muon
# combinatorial that data also has (gap G9 in `add-jpsix-histmaker`). A template
# built from one third of the combinatorial would look like a measurement of all
# of it. Combinatorial is modelled in data by the analytic background shape.
#
# It is a group rather than a bare bin so that the central plotting tool can draw
# it on request -- the composition is a result (77.8% of the gated peak window
# against 21.4% signal) and looking at it should not require a bespoke script.
DIAGNOSTIC_GROUPS = [("NoBAncestor", GEN_NO_B_ANCESTOR)]

# The groups that may enter a fit. Consumers that mean "the simulated templates"
# must use this rather than "every group that is not Data", or they will silently
# fold the combinatorial diagnostic into the signal.
TEMPLATE_GROUP_NAMES = [name for name, _ in MC_GROUPS]


def _slice_gen_category(category):
    """Return an op selecting one genCategory bin and dropping the axis.

    Histograms that carry no `genCategory` axis are passed through untouched, so
    a group definition does not have to know which histograms were booked with
    it.
    """

    def op(h, category=category):
        if GEN_CATEGORY_AXIS not in [ax.name for ax in h.axes]:
            return h
        return h[{GEN_CATEGORY_AXIS: complex(0, category)}]

    return op


def make_datagroups_cvhnano(
    dg, combine=False, pseudodata_pdfset=None, excludeGroups=None, filterGroups=None
):
    dg.groups = {}

    data_members = dg.get_members_from_results(is_data=True)
    mc_members = dg.get_members_from_results()

    dg.addGroup(
        "Data",
        members=data_members,
        memberOp=[_slice_gen_category(GEN_DATA)] * len(data_members),
    )

    for name, category in MC_GROUPS + DIAGNOSTIC_GROUPS:
        dg.addGroup(
            name,
            members=mc_members,
            memberOp=[_slice_gen_category(category)] * len(mc_members),
        )

    dg.filterGroups(filterGroups)
    dg.excludeGroups(excludeGroups)
