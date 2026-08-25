"""Axes and fit-arm description for the CVH NanoAOD B+ -> J/psi K analysis.

The CVH NanoAOD is produced from our TkAlJpsiX AlCaReco by
``Analysis/HitAnalyzer/test/runCvhBplusJpsiK.py``. It shares no physics branch
name with the BMM-Tools NanoAOD that ``btojpsik_axes.py`` describes, so the two
modules are deliberately disjoint.

Every candidate carries three independent fit arms with identical column
layouts, which is the repetition worth parameterising -- see ``ArmSpec``. The
channel itself is a constant: a descriptor over the one registered channel
would only obscure the column names.
"""

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import hist

# Candidate collection prefix in the NanoAOD. One channel is registered; a
# second three-leg channel (Bc -> J/psi pi) would change this and BACHELOR_MASS
# and nothing else.
CAND = "BuJpsiK"
BACHELOR = "kaon"
BACHELOR_MASS = 0.49368  # charged kaon, GeV

# Failed fits write -99 rather than NaN, so a one-sided kinematic bound would
# happily accept them. Anything at or below this is a sentinel, not a value.
SENTINEL = -98.0

# genCategory axis encoding.
#
# Every bin below GEN_NO_B_ANCESTOR is a *real b-hadron decay*, identified from
# the per-leg generator parentage the NanoAOD now carries. The previous encoding
# was three bins deduced from a candidate-level direction match; it is retired
# because that match cannot tell B+ -> J/psi K+ from B+ -> J/psi K*+, and the
# two do not belong in the same template (the latter sits ~200 MeV low with no
# entries at all in the peak window).
#
# GEN_NO_B_ANCESTOR is filled and is deliberately NOT a template: it is
# combinatorial, and combinatorial is modelled in data by the analytic
# background, together with the prompt and fake components this simulation does
# not contain (gap G9). Keeping the bin costs nothing and keeps the composition
# denominators available; `datagroups_cvhnano` is what excludes it from plots.
#
# The last bin is reserved for data so that data and simulation histograms are
# conformable without data sharing a bin with a simulation category.
GEN_SIGNAL = 0  # B+ -> J/psi K+, every leg a direct daughter
GEN_BU_KX = 1  # B+ with a kaon bachelor, but not that decay
GEN_BU_PIX = 2  # B+ with a pion bachelor
GEN_BD = 3  # B0
GEN_BS = 4  # Bs
GEN_OTHER_B = 5  # any other b hadron
GEN_NO_B_ANCESTOR = 6  # no b-hadron common ancestor -- not a template
GEN_DATA = 7
N_GEN_CATEGORIES = 8
GEN_CATEGORY_LABELS = {
    GEN_SIGNAL: "BuJpsiK",
    GEN_BU_KX: "BuJpsiKX",
    GEN_BU_PIX: "BuJpsiPiX",
    GEN_BD: "Bd",
    GEN_BS: "Bs",
    GEN_OTHER_B: "OtherB",
    GEN_NO_B_ANCESTOR: "NoBAncestor",
    GEN_DATA: "Data",
}

# Simulated categories that are real b-hadron decays, in stack order. The
# combinatorial bin is absent on purpose; see above.
GEN_TEMPLATE_CATEGORIES = [
    GEN_SIGNAL,
    GEN_BU_KX,
    GEN_BU_PIX,
    GEN_BD,
    GEN_BS,
    GEN_OTHER_B,
]

# All simulated bins, and all simulated bins that are not the signal. Consumers
# that need "everything under the peak that is not signal" must use the second
# of these rather than listing bins, so that adding a category cannot silently
# drop it from a background sum.
GEN_MC_CATEGORIES = [
    GEN_SIGNAL,
    GEN_BU_KX,
    GEN_BU_PIX,
    GEN_BD,
    GEN_BS,
    GEN_OTHER_B,
    GEN_NO_B_ANCESTOR,
]
GEN_NON_SIGNAL_CATEGORIES = [c for c in GEN_MC_CATEGORIES if c != GEN_SIGNAL]

# Ground-state |pdgId| of the b hadron each simulated category selects. The
# bachelor's own species then splits the channel's own mother species.
CATEGORY_MOTHER = {
    GEN_SIGNAL: 521,
    GEN_BU_KX: 521,
    GEN_BU_PIX: 521,
    GEN_BD: 511,
    GEN_BS: 531,
}

# |pdgId| of the bachelor hypothesis and of the intermediate resonance. The
# muons must come from a J/psi, and that J/psi must itself be a direct daughter
# of the b hadron -- charmonium feed-down (chi_c1 20443, psi(2S) 100443,
# chi_c0 10441) is 1.9% of B+-ancestor candidates and puts ZERO entries in the
# peak window, so it is not signal.
KAON_PDGID = 321
PION_PDGID = 211
MUON_PDGID = 13
JPSI_PDGID = 443

# Mother species for the signal truth category, |pdgId|.
#
# The producer's `genBPdgId` is the nearest *last-copy* b hadron in dR
# (JpsiXKinematicFitProducer.cc). B*+ -> B+ gamma is a forced decay, so the B*+
# is itself a last copy sitting essentially collinear with the B+ that actually
# decayed to J/psi K+ -- and it wins the dR comparison about a quarter of the
# time. Requiring |pdgId| == 521 alone therefore throws away 25% of genuine
# signal. Excited states must be folded onto their ground state before the
# species test; the same holds for B*0/B0 and Bs*/Bs, and for the Sigma_b
# states, which are not in the fragment's forced-decay list and simply decay
# strongly to Lambda_b pi.
MOTHER_PDGID = 521  # B+
GROUND_STATE = {
    523: 521,  # B*+  -> B+
    513: 511,  # B*0  -> B0
    533: 531,  # Bs*  -> Bs
    543: 541,  # Bc*  -> Bc
    5112: 5122,  # Sigma_b-   -> Lambda_b
    5114: 5122,  # Sigma_b*-  -> Lambda_b
    5212: 5122,  # Sigma_b0   -> Lambda_b
    5214: 5122,  # Sigma_b*0  -> Lambda_b
    5222: 5122,  # Sigma_b+   -> Lambda_b
    5224: 5122,  # Sigma_b*+  -> Lambda_b
}

# Forced-decay branching fractions from the generator fragment
# (condor/mc_inclusive_btojpsix_2016postvfp/fragments/BPH-RunIISummer20UL16GEN-00017-fragment.py).
#
# EvtGen was told to decay every b hadron in `list_forced_decays` into a J/psi
# mode with total probability 1. The "original total forced BR" comment on each
# Decay block is the physical probability that was divided out. Simulated events
# therefore over-represent each species by 1/BR, and the factor differs by three
# orders of magnitude across species (Bc 0.911, Lambda_b 0.00085), so the
# *composition* of the stacked simulation is wrong until these are re-applied.
#
# Keyed by ground-state |pdgId|. J/psi -> mu mu was likewise forced from 0.0593
# to 1, but that factor is common to every species and so cancels in any
# composition; it is applied too, so the weight is the true decay probability.
FORCED_BR = {
    521: 0.01756960,  # B+
    511: 0.01778030,  # B0
    531: 0.02298000,  # Bs
    541: 0.91097820,  # Bc+
    5122: 0.00085000,  # Lambda_b0
    5132: 0.00085000,  # Xi_b-
    5232: 0.00047000,  # Xi_b0
    5332: 0.00085000,  # Omega_b-
}
BR_JPSI_MUMU = 0.05930000

# Charmonium feed-down (psi(2S) 0.8238, chi_c1 0.344, ...) was forced too, and
# is a second per-event factor that FORCED_BR does not capture. It is now
# *identifiable* -- the decay chain can be walked, and feed-down is 1.9% of
# B+-ancestor candidates -- but the per-state factors are not applied, because
# feed-down candidates are excluded from the signal category rather than
# reweighted into it. See the gap record.


@dataclass
class ArmSpec:
    """One fit arm of the candidate.

    ``kvfRaw``/``kvfCvh`` are Kalman vertex fits on raw and CVH-refit tracks;
    ``jointCvh`` is the joint N-body CVH fit that the calibration exists to
    characterise.

    ``diagnostics`` is None for arms that carry no convergence information.
    Those arms are gated on ``ok_col`` and the sentinel check alone -- there is
    nothing else to gate on, and inventing a gate would be worse than none.
    """

    name: str
    label: str
    # Columns, without the CAND prefix.
    ok: str
    mass: str
    masserr: str
    pt: str
    geom: Dict[str, str] = field(default_factory=dict)
    # Vertex probability used to rank candidates. The KVF arms store one; the
    # joint arm stores chi2/ndof instead, from which it is computed.
    vtxprob: Optional[str] = None
    vtxchi2: Optional[str] = None
    # Dimuon vertex block. Only the KVF arms have one: the joint arm imposes a
    # single common vertex on all legs, so no distinct dimuon vertex exists.
    dimuon_vtx: Dict[str, str] = field(default_factory=dict)
    # Fitted dimuon kinematics. Only the joint arm has these.
    dimuon_kin: Dict[str, str] = field(default_factory=dict)
    # Convergence diagnostics, or None.
    diagnostics: Optional[Dict[str, str]] = None
    # Per-leg fitted kinematics, joint arm only.
    legs: Dict[str, str] = field(default_factory=dict)

    def col(self, short: str) -> str:
        """Full NanoAOD column name for a short name belonging to this arm."""
        return f"{CAND}_{short}"

    @property
    def ok_col(self) -> str:
        return self.col(self.ok)

    @property
    def mass_col(self) -> str:
        return self.col(self.mass)

    @property
    def sentinel_cols(self) -> List[str]:
        """Every column consumed from this arm, which must not be a sentinel.

        Deliberately includes the diagnostics: a fit that wrote -99 into its
        chi-square has not merely failed a bound, it has no chi-square.
        """
        shorts = [self.mass, self.masserr, self.pt]
        shorts += list(self.geom.values())
        shorts += list(self.dimuon_vtx.values())
        shorts += list(self.dimuon_kin.values())
        shorts += list(self.legs.values())
        if self.vtxprob:
            shorts.append(self.vtxprob)
        if self.vtxchi2:
            shorts.append(self.vtxchi2)
        if self.diagnostics:
            shorts += list(self.diagnostics.values())
        # dict preserves insertion order and de-duplicates
        return [self.col(s) for s in dict.fromkeys(shorts)]


_GEOM = ["alphaBS", "lxy", "sxy", "l3d", "sl3d", "alpha3dPV"]


def _geom(prefix: str) -> Dict[str, str]:
    return {g: f"{prefix}{g[0].upper()}{g[1:]}" for g in _GEOM}


# Unconstrained dimuon mass per arm (gap G11, closed). The joint arm stores only
# the *constrained* J/psi mass, pinned to the PDG value for every candidate, so
# this is the analysis' only fitted dimuon mass -- and the muon momentum-scale
# handle that ties this channel to the J/psi calibration.
def _dimuon_mass(prefix: str) -> Dict[str, str]:
    return {"mass": f"{prefix}DimuonMass", "massErr": f"{prefix}DimuonMassErr"}


def _dimuon_vtx(prefix: str) -> Dict[str, str]:
    d = {
        "vtxProb": f"{prefix}DimuonVtxProb",
        "alphaBS": f"{prefix}DimuonAlphaBS",
        "sxy": f"{prefix}DimuonSxy",
        "sl3d": f"{prefix}DimuonSl3d",
    }
    d.update(_dimuon_mass(prefix))
    return d


ARMS: Dict[str, ArmSpec] = {
    "kvfRaw": ArmSpec(
        name="kvfRaw",
        label="KVF, raw tracks",
        ok="kvfRawOk",
        mass="kvfRawMass",
        masserr="kvfRawMassErr",
        pt="kvfRawPt",
        geom=_geom("kvfRaw"),
        vtxprob="kvfRawVtxProb",
        vtxchi2="kvfRawVtxChi2",
        dimuon_vtx=_dimuon_vtx("kvfRaw"),
        diagnostics=None,
    ),
    "kvfCvh": ArmSpec(
        name="kvfCvh",
        label="KVF, CVH-refit tracks",
        ok="kvfCvhOk",
        mass="kvfCvhMass",
        masserr="kvfCvhMassErr",
        pt="kvfCvhPt",
        geom=_geom("kvfCvh"),
        vtxprob="kvfCvhVtxProb",
        vtxchi2="kvfCvhVtxChi2",
        dimuon_vtx=_dimuon_vtx("kvfCvh"),
        diagnostics=None,
    ),
    "jointCvh": ArmSpec(
        name="jointCvh",
        label="joint N-body CVH",
        ok="jointCvhOk",
        mass="jointCvhMass",
        masserr="jointCvhMassErr",
        pt="jointCvhPt",
        geom=_geom("jointCvh"),
        # No stored vertex probability: computed from chisq/ndof.
        vtxprob=None,
        vtxchi2=None,
        dimuon_kin={
            "mass": "jointCvhJpsiMass",
            "massErr": "jointCvhJpsiMassErr",
            "pt": "jointCvhJpsiPt",
            "eta": "jointCvhJpsiEta",
            "phi": "jointCvhJpsiPhi",
        },
        diagnostics={
            "edmRef": "jointCvhEdmRef",
            "chisq": "jointCvhChisq",
            "ndof": "jointCvhNdof",
            "massErr": "jointCvhMassErr",
        },
        legs={
            f"leg{i}{q}": f"jointCvhLeg{i}{q}"
            for i in (0, 1, 2)
            for q in ("Pt", "Eta", "Phi")
        },
    ),
}

DEFAULT_ARM = "jointCvh"

# The joint arm has no dimuon vertex of its own, so dimuon-vertex selections
# and histograms borrow this arm's block. Same track source (CVH refit), so the
# borrowed quantities are the closest available match.
DIMUON_VTX_FALLBACK = "kvfCvh"


def dimuon_vtx_arm(arm: ArmSpec) -> ArmSpec:
    """Arm supplying the dimuon vertex block for selections on `arm`."""
    return arm if arm.dimuon_vtx else ARMS[DIMUON_VTX_FALLBACK]


# ---------------------------------------------------------------------------
# Axes
# ---------------------------------------------------------------------------

axis_gen_category = hist.axis.Integer(
    0, N_GEN_CATEGORIES, name="genCategory", underflow=False, overflow=False
)


def _reg(n, lo, hi, name, **kw):
    kw.setdefault("underflow", False)
    kw.setdefault("overflow", True)
    return hist.axis.Regular(n, lo, hi, name=name, **kw)


def _arm_axes(arm: ArmSpec) -> Dict[str, hist.axis.AxesMixin]:
    """Axes for one fit arm. Ranges cover the observed distribution.

    The mother mass range is the AlCaReco window [5.0, 5.5]; the fitted arms
    can land outside it, so these axes overflow.
    """
    a = {}
    p = arm.name
    a[arm.mass_col] = _reg(100, 5.0, 5.5, f"{p}_mass")
    a[arm.col(arm.masserr)] = _reg(60, 0.0, 0.12, f"{p}_massErr")
    a[arm.col(arm.pt)] = _reg(60, 0.0, 60.0, f"{p}_pt")
    a[arm.col(arm.geom["alphaBS"])] = _reg(64, 0.0, math.pi, f"{p}_alphaBS")
    a[arm.col(arm.geom["lxy"])] = _reg(60, 0.0, 6.0, f"{p}_lxy")
    a[arm.col(arm.geom["sxy"])] = _reg(60, 0.0, 120.0, f"{p}_sxy")
    a[arm.col(arm.geom["l3d"])] = _reg(60, 0.0, 6.0, f"{p}_l3d")
    a[arm.col(arm.geom["sl3d"])] = _reg(60, 0.0, 120.0, f"{p}_sl3d")
    # 3D pointing angle to the PV (gap G8, closed). The transverse alphaBS has
    # a long tail because the missing daughter tilts the candidate out of the
    # plane; the 3D angle is the one the displacement cut should have used.
    a[arm.col(arm.geom["alpha3dPV"])] = _reg(64, 0.0, math.pi, f"{p}_alpha3dPV")
    if arm.vtxprob:
        a[arm.col(arm.vtxprob)] = _reg(50, 0.0, 1.0, f"{p}_vtxProb")
    if arm.vtxchi2:
        a[arm.col(arm.vtxchi2)] = _reg(60, 0.0, 60.0, f"{p}_vtxChi2")
    for key, short in arm.dimuon_vtx.items():
        rng = {
            "vtxProb": (50, 0.0, 1.0),
            "alphaBS": (64, 0.0, math.pi),
            "sxy": (60, 0.0, 120.0),
            "sl3d": (60, 0.0, 120.0),
            "mass": (80, 2.9, 3.3),
            "massErr": (60, 0.0, 0.06),
        }[key]
        a[arm.col(short)] = _reg(*rng, f"{p}_dimuon_{key}")
    for key, short in arm.dimuon_kin.items():
        rng = {
            "mass": (80, 2.9, 3.3),
            "massErr": (60, 0.0, 0.06),
            "pt": (60, 0.0, 60.0),
            "eta": (48, -2.4, 2.4),
            "phi": (32, -math.pi, math.pi),
        }[key]
        a[arm.col(short)] = _reg(*rng, f"{p}_jpsi_{key}")
    for key, short in arm.legs.items():
        q = key[-2:]
        rng = {
            "Pt": (60, 0.0, 30.0),
            "Et": (48, -2.4, 2.4),  # "Eta"[-2:] == "ta"; handled below
            "Ph": (32, -math.pi, math.pi),
        }
        if short.endswith("Pt"):
            a[arm.col(short)] = _reg(60, 0.0, 30.0, f"{p}_{key}")
        elif short.endswith("Eta"):
            a[arm.col(short)] = _reg(48, -2.4, 2.4, f"{p}_{key}")
        else:
            a[arm.col(short)] = _reg(32, -math.pi, math.pi, f"{p}_{key}")
    if arm.diagnostics:
        # edmRef spans ~1e-10 to 1e2 and the gate sits at 1e-2, so it is only
        # readable on a log scale. narf cannot convert transformed axes, so the
        # edges are log-spaced explicitly instead.
        a[arm.col(arm.diagnostics["edmRef"])] = hist.axis.Variable(
            [10 ** (-11 + 0.25 * i) for i in range(53)],
            name=f"{p}_edmRef",
            underflow=True,
            overflow=True,
        )
        a[arm.col(arm.diagnostics["chisq"])] = _reg(60, 0.0, 600.0, f"{p}_chisq")
        a[arm.col(arm.diagnostics["ndof"])] = hist.axis.Integer(
            0, 90, name=f"{p}_ndof", underflow=False, overflow=True
        )
    return a


def _shared_axes() -> Dict[str, hist.axis.AxesMixin]:
    """Axes that do not belong to a fit arm."""
    c = CAND
    a = {
        # Raw (pre-fit) candidate
        f"{c}_mass": _reg(100, 5.0, 5.5, "cand_mass"),
        f"{c}_pt": _reg(60, 0.0, 60.0, "cand_pt"),
        f"{c}_eta": _reg(48, -2.4, 2.4, "cand_eta"),
        f"{c}_phi": _reg(32, -math.pi, math.pi, "cand_phi"),
        f"{c}_charge": hist.axis.Integer(
            -2, 3, name="cand_charge", underflow=False, overflow=False
        ),
        f"{c}_vertexChi2": _reg(60, 0.0, 60.0, "cand_vertexChi2"),
        f"{c}_nDau": hist.axis.Integer(
            0, 5, name="cand_nDau", underflow=False, overflow=False
        ),
        # Bachelor, raw. Median pT is 0.52 GeV and p95 1.9 GeV -- a regular axis
        # over the full range, NOT the [1,2,3,8] fitting binning of the old
        # analysis, which would put most candidates in its underflow.
        f"{c}_{BACHELOR}Pt": _reg(60, 0.0, 12.0, "bachelor_pt"),
        f"{c}_{BACHELOR}Eta": _reg(48, -2.4, 2.4, "bachelor_eta"),
        f"{c}_{BACHELOR}Phi": _reg(32, -math.pi, math.pi, "bachelor_phi"),
        f"{c}_jpsiPt": _reg(60, 0.0, 60.0, "jpsi_raw_pt"),
        f"{c}_nLegsRefit": hist.axis.Integer(
            0, 5, name="nLegsRefit", underflow=False, overflow=False
        ),
        f"{c}_jointCvhNiter": hist.axis.Integer(
            0, 12, name="jointCvh_nIter", underflow=False, overflow=True
        ),
        # Event level. Keys are the dataframe column that fills them.
        f"n{c}": hist.axis.Integer(0, 15, name="nCand", underflow=False, overflow=True),
        "nPV": hist.axis.Integer(0, 60, name="nPV", underflow=False, overflow=True),
        # True pileup (gap G4, closed). nTrueInt is the Poisson mean, so it is
        # a float; the axis is integer-edged because the reweighting is, and
        # because the official data profiles are binned that way.
        "Pileup_nTrueInt": _reg(60, 0.0, 60.0, "nTrueInt"),
        "Pileup_nPU": hist.axis.Integer(
            0, 60, name="nPU", underflow=False, overflow=True
        ),
        "nTrack": hist.axis.Integer(
            0, 150, name="nTrack", underflow=False, overflow=True
        ),
        "nMuon": hist.axis.Integer(0, 15, name="nMuon", underflow=False, overflow=True),
        "PV_z": _reg(60, -15.0, 15.0, "PV_z"),
        # Per-candidate legs, pulled through the Track cross-links by
        # selections.define_leg_kinematics. Keys match those column names
        # exactly, so no key-to-column translation is needed anywhere.
        "cand_mu0Pt": _reg(60, 0.0, 30.0, "mu0_pt"),
        "cand_mu0Eta": _reg(48, -2.4, 2.4, "mu0_eta"),
        "cand_mu1Pt": _reg(60, 0.0, 30.0, "mu1_pt"),
        "cand_mu1Eta": _reg(48, -2.4, 2.4, "mu1_eta"),
        # Raw, origin-referenced (gap G13). Kept so the defect stays visible.
        "cand_bachelorDxy": _reg(60, -0.5, 0.5, "bachelor_dxy"),
        "cand_bachelorDz": _reg(60, -2.0, 2.0, "bachelor_dz"),
        # PV-referenced, which is what an impact parameter is supposed to be.
        # Ranges are tighter because the origin-referenced versions are
        # dominated by the beam offset rather than by the physics.
        "cand_bachelorD0": _reg(80, -0.2, 0.2, "bachelor_d0"),
        "cand_bachelorDzPV": _reg(80, -1.0, 1.0, "bachelor_dzPV"),
        "cand_bachelorNormChi2": _reg(60, 0.0, 6.0, "bachelor_normChi2"),
        # Impact-parameter SIGNIFICANCE, the standard combinatorial
        # discriminant, which had no input before (gap G6). Signed, so the
        # symmetry of the distribution is visible rather than folded away.
        "cand_bachelorD0Sig": _reg(60, -30.0, 30.0, "bachelor_d0Sig"),
        "cand_bachelorDzPVSig": _reg(60, -30.0, 30.0, "bachelor_dzPVSig"),
        # Relative pT uncertainty rather than the absolute one: the absolute
        # error is a proxy for pT itself and says nothing on its own.
        "cand_bachelorRelPtErr": _reg(60, 0.0, 0.12, "bachelor_relPtErr"),
        "cand_bachelorNValidPixelHits": hist.axis.Integer(
            0, 10, name="bachelor_nValidPixelHits", underflow=False, overflow=True
        ),
        "cand_bachelorTrackerLayers": hist.axis.Integer(
            0, 20, name="bachelor_trackerLayers", underflow=False, overflow=True
        ),
        "cand_bachelorHighPurity": hist.axis.Integer(
            0, 2, name="bachelor_highPurity", underflow=False, overflow=False
        ),
        "cand_bachelorNValidHits": _reg(40, 0.0, 40.0, "bachelor_nValidHits"),
        "cand_bachelorDedxHarmonic2": _reg(60, 0.0, 12.0, "bachelor_dedx"),
        # Dimuon rebuilt from the raw muon tracks. The NanoAOD's only J/psi
        # mass is the joint fit's CONSTRAINED parameter, pinned to the PDG
        # value for every candidate, so it carries no information (gap G11).
        # This is the only unconstrained dimuon mass available.
        "cand_dimuonRawMass": _reg(80, 2.9, 3.3, "dimuon_raw_mass"),
        "cand_dimuonRawPt": _reg(60, 0.0, 60.0, "dimuon_raw_pt"),
        "cand_dimuonRawEta": _reg(48, -2.4, 2.4, "dimuon_raw_eta"),
        # Truth. `genB*` is the old candidate-level direction match, kept for
        # one cycle so the parentage categorisation can be compared against it.
        f"{c}_genBPt": _reg(60, 0.0, 60.0, "genB_pt"),
        f"{c}_genBEta": _reg(48, -2.4, 2.4, "genB_eta"),
        f"{c}_genBDR": _reg(60, 0.0, 0.5, "genB_dR"),
        # Per-leg parentage (gap G2, closed). Legs 0 and 1 are the muons, leg 2
        # the bachelor -- confirmed empirically, not assumed: leg 2 is the
        # bachelor for 139 185 of 139 185 candidates over 30 files.
        f"{c}_nLegsGenMatched": hist.axis.Integer(
            0, 4, name="nLegsGenMatched", underflow=False, overflow=False
        ),
        f"{c}_leg2GenPt": _reg(60, 0.0, 12.0, "bachelor_genPt"),
        f"{c}_leg2GenDR": _reg(60, 0.0, 0.02, "bachelor_genDR"),
        # Kaon momentum response, the direct handle on the bachelor momentum
        # scale. Sentinel-guarded and only defined where the leg matched, so
        # unmatched candidates do not pile up at zero.
        "cand_bachelorPtResponse": _reg(60, 0.90, 1.10, "bachelor_ptResponse"),
    }
    return a


def build_axes() -> Dict[str, hist.axis.AxesMixin]:
    """All available axes, keyed by the dataframe column that fills them."""
    axes = _shared_axes()
    for arm in ARMS.values():
        axes.update(_arm_axes(arm))
    return axes


all_axes = build_axes()

# Columns that only exist in simulation, so they are skipped for data.
MC_ONLY_AXES = {
    f"{CAND}_genBPt",
    f"{CAND}_genBEta",
    f"{CAND}_genBDR",
    f"{CAND}_nLegsGenMatched",
    f"{CAND}_leg2GenPt",
    f"{CAND}_leg2GenDR",
    "cand_bachelorPtResponse",
    "Pileup_nTrueInt",
    "Pileup_nPU",
}

# Event-level scalars, which need no per-candidate reduction.
EVENT_SCALAR_AXES = {
    f"n{CAND}",
    "nPV",
    "nTrack",
    "nMuon",
    "Pileup_nTrueInt",
    "Pileup_nPU",
}


# --- axes for the differential fit binning ---------------------------------
# The eventual fit bins mass in cells of bachelor (pT, eta, charge). eta is on
# its intended final binning: 48 bins of 0.1 over [-2.4, 2.4]. pT is left fine
# and regular here because the quantile edges have to be derived *from* this
# histogram -- committing to quantile edges before measuring the occupancy
# would be circular.
axis_fit_bachelor_pt = hist.axis.Regular(
    100, 0.0, 10.0, name="fit_bachelor_pt", underflow=False, overflow=True
)
axis_fit_bachelor_eta = hist.axis.Regular(
    48, -2.4, 2.4, name="fit_bachelor_eta", underflow=False, overflow=False
)
# Two bins, +1 and -1, matching the old analysis' `bkmm_kaon_charge`. NOT the
# five-bin `cand_charge` axis above: that one is a diagnostic, kept wide on
# purpose so a pathological candidate charge shows in a plot rather than being
# clipped into an edge bin. Using it here put three permanently empty bins into
# the tensor and made every consumer carry a rule about which bins are real.
#
# There is no kaon-charge branch in this NanoAOD. `BuJpsiK_charge` is the
# composite candidate charge, equal to the kaon charge because the muons are
# required to be opposite-sign, so it is what fills this axis. `kaonPdgId` has
# the same sign but reading a charge off a particle-identity column would
# suggest truth is involved, and it is not -- the pdgId is the +/-321
# hypothesis.
# --- axes for the kaon-ID study -------------------------------------------
# The discriminants are studied against bachelor pT, not inclusively, because
# every one of them depends on it: displacement significance improves with
# momentum, and dE/dx separates K from pi only below ~0.7 GeV/c. A 1D histogram
# of any of these answers the wrong question.
#
# Coarse pT bins on purpose -- the point is overlays a reader can tell apart, and
# the boundaries are placed where the physics changes: 0.7 GeV is where
# ionisation stops separating kaons from pions, and 1 GeV is the 2018 production
# floor this is meant to replace.
axis_id_bachelor_pt = hist.axis.Variable(
    [0.0, 0.4, 0.7, 1.0, 1.5, 2.5, 4.0, 20.0],
    name="id_bachelor_pt",
    underflow=False,
    overflow=False,
)

# Quality variables for the kaon-ID study, keyed by the column that fills them.
# Signed where the sign carries information: a symmetric d0 distribution is the
# signature of a prompt track, and folding it away would hide that.
ID_AXES = {
    "cand_bachelorD0": _reg(80, -0.1, 0.1, "id_d0"),
    "cand_bachelorD0Sig": _reg(80, -20.0, 20.0, "id_d0Sig"),
    "cand_bachelorDzPV": _reg(80, -0.5, 0.5, "id_dzPV"),
    "cand_bachelorDzPVSig": _reg(80, -20.0, 20.0, "id_dzPVSig"),
    "cand_bachelorDedxHarmonic2": _reg(70, 0.0, 14.0, "id_dedx"),
    "cand_bachelorNValidHits": hist.axis.Integer(
        0, 36, name="id_nValidHits", underflow=False, overflow=True
    ),
    "cand_bachelorNormChi2": _reg(50, 0.0, 5.0, "id_normChi2"),
    "cand_bachelorRelPtErr": _reg(60, 0.0, 0.12, "id_relPtErr"),
}

axis_fit_bachelor_charge = hist.axis.Regular(
    2, -2.0, 2.0, name="fit_bachelor_charge", underflow=False, overflow=False
)
