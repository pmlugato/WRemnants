"""Selections for the CVH NanoAOD B+ -> J/psi K analysis.

Candidates live in one flat per-event collection (``BuJpsiK_*``). Every
selection narrows a per-candidate boolean mask ``cand_pass`` and then requires
at least one surviving candidate, so the cutflow counts events with at least
one candidate still alive.

Two properties of this NanoAOD shape the module and are worth stating once:

* Failed fits write ``-99``, not NaN. A one-sided kinematic bound accepts them,
  so sentinel rejection is a separate selection applied *before* any kinematic
  cut on that arm.
* ``jointCvhOk == 1`` is a convergence flag, not a quality verdict: it admits
  fits with chi-square up to 8e29. The quality gate is therefore mandatory and
  each of its bounds is a separate, individually counted selection.
"""

import json
import os

from wremnants.production.btojpsik_cvhnano_axes import (
    CAND,
    CATEGORY_MOTHER,
    FORCED_BR,
    GEN_BD,
    GEN_BS,
    GEN_BU_KX,
    GEN_BU_PIX,
    GEN_DATA,
    GEN_NO_B_ANCESTOR,
    GEN_OTHER_B,
    GEN_SIGNAL,
    GROUND_STATE,
    JPSI_PDGID,
    KAON_PDGID,
    MOTHER_PDGID,
    SENTINEL,
    dimuon_vtx_arm,
)
from wums import logging

logger = logging.child_logger(__name__)


# ---------------------------------------------------------------------------
# Input file hygiene
# ---------------------------------------------------------------------------


def select_good_files(filepaths):
    """Drop NanoAOD files whose producing job did not report success.

    The production directories contain truncated stubs from killed jobs: a file
    that exists is not a file that is complete. Each NanoAOD has a JSON sidecar
    written by the condor wrapper; the completion marker is
    ``exit_reason == "ok"``, and where the wrapper recorded ``out_bytes`` the
    file must still have exactly that many bytes.

    Files produced before ``out_bytes`` was added carry no size to check, so
    the sidecar's exit reason is all there is; that is still strictly better
    than trusting existence.
    """
    good, dropped = [], []
    for path in filepaths:
        sidecar = path[:-5] + ".json" if path.endswith(".root") else path + ".json"
        if not os.path.exists(sidecar):
            dropped.append((path, "no sidecar"))
            continue
        try:
            with open(sidecar, "r", encoding="utf-8") as handle:
                record = json.load(handle)
        except (OSError, ValueError) as exc:
            dropped.append((path, f"unreadable sidecar ({exc})"))
            continue
        if record.get("exit_reason") not in ("ok", "already_done"):
            dropped.append((path, f"exit_reason={record.get('exit_reason')!r}"))
            continue
        expected = record.get("out_bytes")
        if expected:
            try:
                actual = os.path.getsize(path)
            except OSError as exc:
                dropped.append((path, f"unstattable ({exc})"))
                continue
            if int(expected) != actual:
                dropped.append((path, f"size {actual} != recorded {expected}"))
                continue
        good.append(path)

    if dropped:
        logger.warning(
            "Dropping %d of %d input files that did not report a clean production.",
            len(dropped),
            len(filepaths),
        )
        for path, reason in dropped[:10]:
            logger.warning("  %s: %s", os.path.basename(path), reason)
        if len(dropped) > 10:
            logger.warning("  ... and %d more", len(dropped) - 10)
    return good


# ---------------------------------------------------------------------------
# Trigger
# ---------------------------------------------------------------------------


def validate_trigger_available(
    filepaths, path, treename="Events", nprocs=16, cache=None
):
    """Fail before the event loop if `path` is missing from any input file.

    The HLT menu is not uniform across files in one directory: 28 of 709
    columns are absent from some 2016H files. RDataFrame only discovers that
    mid-loop, after the run has been going for hours, so the check is done up
    front. Every file is checked, not a sample -- one missing column anywhere is
    enough to kill the run.

    **This used a ThreadPoolExecutor and did not parallelise.** The work looks
    like pure I/O wait, but `TFile.Open` and `GetBranch` go through PyROOT,
    which holds the GIL for the whole call, so sixteen threads ran one at a
    time. Measured at 125 ms per file single-threaded, 22 327 files cost over an
    hour of wall clock before the event loop began -- against the "minutes" this
    docstring used to claim.

    Two changes:

    * processes instead of threads, so the GIL is not shared;
    * an on-disk cache keyed by the file list, because the answer for a given
      production does not change and re-deriving it on every run is the real
      waste. Pass `cache=<path>` to enable it.
    """
    import hashlib
    import json
    import os
    from concurrent.futures import ProcessPoolExecutor

    import ROOT

    column = path if path.startswith("HLT_") else f"HLT_{path}"

    key = hashlib.sha1(
        ("\n".join(sorted(filepaths)) + "|" + column + "|" + treename).encode()
    ).hexdigest()
    if cache and os.path.exists(cache):
        try:
            with open(cache) as handle:
                seen = json.load(handle)
        except (OSError, ValueError):
            seen = {}
        if seen.get(key) == "ok":
            logger.info(
                "%s already verified present in these %d files (cache %s).",
                column,
                len(filepaths),
                cache,
            )
            return column

    logger.info(
        "Checking %s is present in all %d input files (%d processes).",
        column,
        len(filepaths),
        nprocs,
    )
    ROOT.ROOT.EnableThreadSafety()

    with ProcessPoolExecutor(max_workers=nprocs) as pool:
        for filepath, reason in pool.map(
            _check_one_file, ((f, treename, column) for f in filepaths), chunksize=8
        ):
            if filepath is not None:
                raise RuntimeError(
                    f"Input file {filepath} {reason}. The HLT menu varies "
                    "between files in this production, so a path must be "
                    "present in every input file before it can be used. Pick a "
                    "path common to all files, or restrict the input set."
                )
    logger.info("%s present in all input files.", column)

    if cache:
        try:
            seen = {}
            if os.path.exists(cache):
                with open(cache) as handle:
                    seen = json.load(handle)
            seen[key] = "ok"
            with open(cache, "w") as handle:
                json.dump(seen, handle)
        except (OSError, ValueError) as exc:
            # A cache that cannot be written is not a reason to fail a run that
            # has already answered the question correctly.
            logger.warning("Could not write the trigger cache %s: %s", cache, exc)

    return column


def _check_one_file(args):
    """One file's verdict. Module level so a process pool can pickle it."""
    import ROOT

    filepath, treename, column = args
    rootfile = ROOT.TFile.Open(filepath)
    if not rootfile or rootfile.IsZombie():
        return filepath, "cannot be opened"
    tree = rootfile.Get(treename)
    present = bool(tree) and tree.GetBranch(column) is not None
    rootfile.Close()
    return (None, None) if present else (filepath, f"has no {column}")


def apply_trigger(df, column):
    """Apply the trigger identically to data and simulation."""
    logger.info("HLT selection: %s", column)
    return df.Filter(column, column)


# ---------------------------------------------------------------------------
# Candidate mask machinery
# ---------------------------------------------------------------------------


def init_candidate_mask(df):
    """Start every candidate alive."""
    return df.Define(
        "cand_pass", f"ROOT::VecOps::RVec<bool>(static_cast<size_t>(n{CAND}), true)"
    )


def _narrow(df, condition):
    """Narrow `cand_pass` by a per-candidate C++ `condition` written in terms of `i`."""
    body = f"""
    ROOT::VecOps::RVec<bool> out(cand_pass.size(), false);
    for (size_t i = 0; i < out.size(); ++i) {{
        if (!cand_pass[i]) continue;
        out[i] = ({condition});
    }}
    return out;
    """
    return df.Redefine("cand_pass", body)


def count_candidates(df, tag):
    """Book a sum of surviving candidates, so cuts can be read per candidate."""
    column = f"cand_n_pass_{tag}"
    df = df.Define(column, "ROOT::VecOps::Sum(cand_pass)")
    return df, df.Sum(column)


def apply_selections(df, selections, cutflow, cand_cutflow=None):
    """Run `selections` in order, recording yields after each.

    Every entry is `(name, callable)`; the callable narrows `cand_pass`. After
    each one the event must still hold at least one live candidate.

    Two yields are recorded because they answer different questions. The event
    count is what the analysis keeps. The candidate count is what a cut
    actually removes -- an event with several candidates survives a cut that
    kills most of them, so the event count understates every per-candidate cut.
    The quality-gate fractions quoted in the proposal are candidate-level.
    """
    dfs_per_cut = []
    for index, (name, action) in enumerate(selections):
        logger.info("selection: %s", name)
        df = action(df)
        if cand_cutflow is not None:
            df, cand_sum = count_candidates(df, f"sel{index}")
            cand_cutflow[name] = cand_sum
        df = df.Filter("ROOT::VecOps::Any(cand_pass)", name)
        cutflow[name] = df.SumAndCount("weight")[0]
        dfs_per_cut.append(df)
    return df, cutflow, dfs_per_cut


# ---------------------------------------------------------------------------
# Derived per-candidate columns
# ---------------------------------------------------------------------------


def define_leg_kinematics(df):
    """Per-candidate leg quantities pulled through the Track cross-links.

    The candidate stores a row index into ``Track`` for each leg. Indices were
    measured to be always valid, but an out-of-range index would read past the
    end of the array, so every lookup is guarded and yields a sentinel instead.
    """
    legs = [
        ("mu0", f"{CAND}_mu0TrackIdx"),
        ("mu1", f"{CAND}_mu1TrackIdx"),
        ("bachelor", f"{CAND}_kaonTrackIdx"),
    ]
    track_vars = [
        ("Pt", "Track_pt", "float"),
        ("Eta", "Track_eta", "float"),
        ("Phi", "Track_phi", "float"),
        # Referenced to the ORIGIN, not to the primary vertex. Kept for one
        # cycle so the fix below can be demonstrated rather than asserted; the
        # phi modulation of the median is 0.187 cm here and 0.002 cm in D0.
        ("Dxy", "Track_dxy", "float"),
        ("Dz", "Track_dz", "float"),
        # Needed by the momentum-scale variations: the A/e/M shift on 1/pt is
        # charge dependent through the M term.
        ("Charge", "Track_charge", "int"),
        # PV-referenced at source (gap G13, closed). These used to be computed
        # downstream from PV_x/PV_y and the track phi, which reached a phi
        # modulation of 0.0003 cm; the producer's version does better and is
        # the same quantity, so the arithmetic is gone.
        ("D0", "Track_d0", "float"),
        ("D0Err", "Track_d0Err", "float"),
        ("DzPV", "Track_dzPV", "float"),
        ("DzPVErr", "Track_dzPVErr", "float"),
        ("PtErr", "Track_ptErr", "float"),
        ("NormChi2", "Track_normChi2", "float"),
        ("NValidHits", "Track_nValidHits", "float"),
        ("NValidPixelHits", "Track_nValidPixelHits", "float"),
        ("TrackerLayers", "Track_trackerLayers", "float"),
        ("HighPurity", "Track_highPurity", "float"),
        ("DedxHarmonic2", "Track_dedxHarmonic2", "float"),
    ]
    for tag, idx_col in legs:
        for suffix, track_col, ctype in track_vars:
            df = df.Define(
                f"cand_{tag}{suffix}",
                f"""
                ROOT::VecOps::RVec<{ctype}> out({idx_col}.size(), {SENTINEL}f);
                for (size_t i = 0; i < out.size(); ++i) {{
                    int idx = {idx_col}[i];
                    if (idx >= 0 && idx < static_cast<int>({track_col}.size()))
                        out[i] = static_cast<{ctype}>({track_col}[idx]);
                }}
                return out;
                """,
            )
        # Ratios of two looked-up columns, so both sentinels have to be tested
        # and neither denominator may be zero. A ratio of a sentinel to a
        # sentinel is 1.0, which no bound would catch.
        for name, num, den in [
            (f"cand_{tag}D0Sig", f"cand_{tag}D0", f"cand_{tag}D0Err"),
            (f"cand_{tag}DzPVSig", f"cand_{tag}DzPV", f"cand_{tag}DzPVErr"),
            (f"cand_{tag}RelPtErr", f"cand_{tag}PtErr", f"cand_{tag}Pt"),
        ]:
            df = df.Define(
                name,
                f"""
                ROOT::VecOps::RVec<float> out({num}.size(), {SENTINEL}f);
                for (size_t i = 0; i < out.size(); ++i) {{
                    if ({num}[i] <= {SENTINEL}f || {den}[i] <= {SENTINEL}f) continue;
                    if ({den}[i] <= 0.f) continue;
                    out[i] = {num}[i] / {den}[i];
                }}
                return out;
                """,
            )
        # Muon-level flags, reached through Track_muonIdx. A track with no
        # associated muon yields false rather than a sentinel: "not a global
        # muon" is the honest reading of "not a muon at all".
        for suffix, muon_col in [
            ("IsGlobal", "Muon_isGlobal"),
            ("IsTracker", "Muon_isTracker"),
        ]:
            df = df.Define(
                f"cand_{tag}{suffix}",
                f"""
                ROOT::VecOps::RVec<bool> out({idx_col}.size(), false);
                for (size_t i = 0; i < out.size(); ++i) {{
                    int idx = {idx_col}[i];
                    if (idx < 0 || idx >= static_cast<int>(Track_muonIdx.size())) continue;
                    int midx = Track_muonIdx[idx];
                    if (midx >= 0 && midx < static_cast<int>({muon_col}.size()))
                        out[i] = {muon_col}[midx];
                }}
                return out;
                """,
            )
    return df


MUON_MASS = 0.1056583745


def define_raw_dimuon(df):
    """Reconstruct the dimuon from the two muon legs' raw track momenta.

    No longer the only unconstrained dimuon mass available: gap G11 is closed and
    the NanoAOD now carries ``kvfRawDimuonMass`` / ``kvfCvhDimuonMass`` with
    their errors, from a genuine vertex fit. Those are the mass-scale
    observables; this stays for two reasons that are not aesthetic. It is the
    only dimuon **pT** there is, which `select_dimuon_pt` needs, and being built
    from the raw track momenta with no vertex constraint at all it is the
    zeroth-order reference the fitted masses are compared against.

    ``jointCvhJpsiMass`` remains useless for this purpose: it is the joint fit's
    *constrained* parameter, pinned to the PDG value for every candidate.
    """
    df = df.Define(
        "cand_dimuonRaw",
        f"""
        ROOT::VecOps::RVec<ROOT::Math::PtEtaPhiMVector> out;
        out.reserve(cand_mu0Pt.size());
        for (size_t i = 0; i < cand_mu0Pt.size(); ++i) {{
            ROOT::Math::PtEtaPhiMVector a(cand_mu0Pt[i], cand_mu0Eta[i],
                                          cand_mu0Phi[i], {MUON_MASS});
            ROOT::Math::PtEtaPhiMVector b(cand_mu1Pt[i], cand_mu1Eta[i],
                                          cand_mu1Phi[i], {MUON_MASS});
            out.push_back(a + b);
        }}
        return out;
        """,
    )
    for name, accessor in [
        ("cand_dimuonRawMass", "M()"),
        ("cand_dimuonRawPt", "Pt()"),
        ("cand_dimuonRawEta", "Eta()"),
    ]:
        df = df.Define(
            name,
            f"""
            ROOT::VecOps::RVec<double> out(cand_dimuonRaw.size(), 0.);
            for (size_t i = 0; i < out.size(); ++i) out[i] = cand_dimuonRaw[i].{accessor};
            return out;
            """,
        )
    return df


def define_arm_vertex_probability(df, arm):
    """Per-candidate vertex probability for ranking, uniform across arms.

    The KVF arms store one. The joint arm does not -- all legs share one vertex
    and it reports chi-square and degrees of freedom instead -- so the
    equivalent upper-tail probability is computed from those. Candidates whose
    chi-square or ndof is a sentinel get probability 0, which sorts them last
    rather than accidentally winning.
    """
    name = f"cand_{arm.name}_vtxProb"
    if arm.vtxprob:
        df = df.Define(name, f"ROOT::VecOps::RVec<double>({arm.col(arm.vtxprob)})")
        return df, name
    chisq = arm.col(arm.diagnostics["chisq"])
    ndof = arm.col(arm.diagnostics["ndof"])
    df = df.Define(
        name,
        f"""
        ROOT::VecOps::RVec<double> out({chisq}.size(), 0.);
        for (size_t i = 0; i < out.size(); ++i) {{
            double c = {chisq}[i];
            double n = {ndof}[i];
            if (std::isfinite(c) && std::isfinite(n) && c > 0. && n > 0.)
                out[i] = ROOT::Math::chisquared_cdf_c(c, n);
        }}
        return out;
        """,
    )
    return df, name


def define_chisq_ndof(df, arm):
    """Joint-fit chi-square per degree of freedom, as `cand_chisqNdof`.

    A non-positive or non-finite denominator yields a large sentinel rather than a
    division: an infinity would pass an upper bound on some platforms and a NaN
    passes every comparison as false, which silently changes the meaning of the
    cut rather than rejecting the candidate.
    """
    if not arm.diagnostics:
        return df
    chisq = arm.col(arm.diagnostics["chisq"])
    ndof = arm.col(arm.diagnostics["ndof"])
    return df.Define(
        "cand_chisqNdof",
        f"""
        ROOT::VecOps::RVec<double> out({chisq}.size(), 1.e9);
        for (size_t i = 0; i < out.size(); ++i) {{
            double c = {chisq}[i];
            double n = {ndof}[i];
            if (std::isfinite(c) && std::isfinite(n) && n > 0. && c >= 0.)
                out[i] = c / n;
        }}
        return out;
        """,
    )


# ---------------------------------------------------------------------------
# Selections
# ---------------------------------------------------------------------------


def select_arm_finite(df, arm):
    """Reject candidates carrying a sentinel or non-finite value in this arm.

    Applied before any kinematic cut on the arm, because a sentinel passes a
    one-sided bound: `-99 < 5.5` is true.
    """
    checks = " && ".join(
        f"(std::isfinite({col}[i]) && {col}[i] > {SENTINEL})"
        for col in arm.sentinel_cols
    )
    return _narrow(df, checks)


def select_arm_ok(df, arm):
    """Require the arm's own success flag."""
    return _narrow(df, f"{arm.ok_col}[i] == 1")


def quality_gate_selections(arm, bounds):
    """Named selections for the arm's convergence diagnostics.

    Returns an empty list for arms that carry no diagnostics -- there is
    nothing to gate on, and a fabricated gate would be worse than none.

    Each bound is separate so the cutflow shows which one bites. `edmRef` is
    the dominant term (15.7% of flagged-good candidates in data); the rest
    guard the pathological tail.
    """
    if not arm.diagnostics:
        logger.info(
            "Arm %s carries no convergence diagnostics; gated on its success "
            "flag and sentinel check alone.",
            arm.name,
        )
        return []

    d = arm.diagnostics
    edm = arm.col(d["edmRef"])
    chisq = arm.col(d["chisq"])
    ndof = arm.col(d["ndof"])
    masserr = arm.col(d["massErr"])
    return [
        (
            f"{arm.name} edmRef < {bounds['edmRef']:g}",
            lambda df: _narrow(df, f"{edm}[i] >= 0. && {edm}[i] < {bounds['edmRef']}"),
        ),
        (
            f"{arm.name} chisq < {bounds['chisq']:g}",
            lambda df: _narrow(
                df, f"{chisq}[i] > 0. && {chisq}[i] < {bounds['chisq']}"
            ),
        ),
        (
            f"{arm.name} ndof < {bounds['ndof']:g}",
            lambda df: _narrow(df, f"{ndof}[i] > 0. && {ndof}[i] < {bounds['ndof']}"),
        ),
        (
            f"{arm.name} massErr < {bounds['massErr']:g}",
            lambda df: _narrow(
                df, f"{masserr}[i] > 0. && {masserr}[i] < {bounds['massErr']}"
            ),
        ),
    ]


def select_opposite_sign_muons(df):
    """The two muon legs must have opposite charge."""
    return _narrow(
        df,
        f"""
        [&]{{
            int i0 = {CAND}_mu0TrackIdx[i], i1 = {CAND}_mu1TrackIdx[i];
            if (i0 < 0 || i1 < 0) return false;
            if (i0 >= (int)Track_charge.size() || i1 >= (int)Track_charge.size()) return false;
            return Track_charge[i0] * Track_charge[i1] < 0;
        }}()
        """,
    )


def select_muon_pt(df, min_pt):
    return _narrow(df, f"cand_mu0Pt[i] > {min_pt} && cand_mu1Pt[i] > {min_pt}")


def select_muon_eta(df, max_abs_eta):
    return _narrow(
        df,
        f"std::abs(cand_mu0Eta[i]) < {max_abs_eta} "
        f"&& std::abs(cand_mu1Eta[i]) < {max_abs_eta}",
    )


def select_muon_global(df):
    """Both muon legs must be reconstructed as global muons.

    This is the only muon-quality handle in the NanoAOD: there is no soft-muon
    MVA (gap G7), so the old analysis's `softMVA > 0.45` has no counterpart.
    """
    return _narrow(df, "cand_mu0IsGlobal[i] && cand_mu1IsGlobal[i]")


def select_bachelor_pt(df, min_pt):
    return _narrow(df, f"{CAND}_{'kaon'}Pt[i] > {min_pt}")


def select_bachelor_eta(df, max_abs_eta):
    return _narrow(df, f"std::abs({CAND}_kaonEta[i]) < {max_abs_eta}")


def select_mass_window(df, arm, low, high):
    return _narrow(df, f"{arm.mass_col}[i] > {low} && {arm.mass_col}[i] < {high}")


def select_vertex_probability(df, prob_col, min_prob):
    return _narrow(df, f"{prob_col}[i] > {min_prob}")


def select_sl3d(df, arm, min_sl3d):
    col = arm.col(arm.geom["sl3d"])
    return _narrow(df, f"{col}[i] > {min_sl3d}")


def select_alpha_bs(df, arm, max_alpha):
    col = arm.col(arm.geom["alphaBS"])
    return _narrow(df, f"{col}[i] < {max_alpha}")


def select_dimuon_vertex_probability(df, arm, min_prob):
    """Dimuon vertex probability, from the arm's own block or the fallback arm."""
    source = dimuon_vtx_arm(arm)
    col = source.col(source.dimuon_vtx["vtxProb"])
    return _narrow(df, f"{col}[i] > {min_prob} && {col}[i] > {SENTINEL}")


def select_dimuon_sl3d(df, arm, min_sl3d):
    source = dimuon_vtx_arm(arm)
    col = source.col(source.dimuon_vtx["sl3d"])
    return _narrow(df, f"{col}[i] > {min_sl3d}")


def select_dimuon_alpha_bs(df, arm, max_alpha):
    """Dimuon pointing angle to the beamspot, from the arm's block or the fallback.

    Distinct from `select_alpha_bs`, which is the *candidate's* pointing angle. The
    previous analysis cut the dimuon's, and for a three-body candidate the two are
    not the same quantity: the bachelor track moves the candidate direction away
    from the J/psi direction by an amount that grows as the kaon hardens.
    """
    source = dimuon_vtx_arm(arm)
    col = source.col(source.dimuon_vtx["alphaBS"])
    return _narrow(df, f"{col}[i] < {max_alpha} && {col}[i] > {SENTINEL}")


def select_dimuon_pt(df, min_pt):
    """Dimuon pT, from the dimuon rebuilt out of the raw muon tracks.

    Uses the raw dimuon rather than a fitted one because no arm stores a fitted
    dimuon *pT*. The unconstrained dimuon *mass* is now stored per arm (gap G11,
    closed), which is what the mass-scale work needs; the momentum is not, so a
    pT cut still has to be taken on the rebuilt quantity. It exists for both
    datasets and every arm, which a per-arm column would not.
    """
    return _narrow(df, f"cand_dimuonRawPt[i] > {min_pt}")


def select_bachelor_pt_max(df, max_pt):
    """Upper bound on bachelor pT.

    The previous analysis capped the kaon at 8 GeV. Worth keeping separate from
    the floor: the floor removes tracks that are barely reconstructed, the ceiling
    removes a region where the kaon is hard enough that the candidate kinematics
    stop looking like a B decay.
    """
    return _narrow(df, f"{CAND}_kaonPt[i] < {max_pt}")


def select_bachelor_n_valid_hits(df, min_hits):
    """Bachelor valid-hit floor.

    The column comes from `define_leg_kinematics`, which pulls it through the
    Track cross-link and names it `cand_bachelor...`, not `BuJpsiK_bachelor...`.
    A failed lookup writes the sentinel, which fails this one-sided cut anyway,
    but is rejected explicitly so the cut does not depend on that coincidence.
    """
    col = "cand_bachelorNValidHits"
    return _narrow(df, f"{col}[i] > {SENTINEL} && {col}[i] > {min_hits}")


def select_bachelor_abs_d0(df, min_abs_d0):
    """Displacement of the bachelor track from the primary vertex.

    PV-referenced, unlike the origin-referenced `bachelorDxy` kept only to keep
    gap G13 visible. A kaon from a B decay is displaced; a prompt track is not.

    The sentinel check is not optional and is not decoration. A failed Track
    lookup writes -98, and this is the one cut in the set taken on an absolute
    value, so |-98| > 0.02 would *pass* -- a broken candidate would be selected
    as maximally displaced.
    """
    col = "cand_bachelorD0"
    return _narrow(df, f"{col}[i] > {SENTINEL} && std::abs({col}[i]) > {min_abs_d0}")


def select_bachelor_abs_d0_sig(df, min_abs_sig):
    """Bachelor impact-parameter significance, |d0| / sigma(d0), w.r.t. the PV.

    The kaon ID the background composition asks for. The dominant background is a
    real J/psi plus a **soft prompt track from the primary vertex** -- two thirds
    pions, essentially none from a b decay -- while the signal kaon comes from the
    displaced B vertex. Significance rather than d0 itself, because the resolution
    varies strongly with momentum and a fixed distance cut is therefore a
    different cut at every pT: measured, significance rejects 0.513 of the
    background at 90% signal efficiency in 0.4 < pT < 0.7 GeV against 0.352 for
    raw d0.

    Sentinel-guarded on both inputs, and this is not optional: the ratio is taken
    on an absolute value, so an unguarded |-98| would pass as maximally
    displaced -- and `cand_bachelorD0Sig` is already -98 when either the distance
    or its error is a sentinel.
    """
    col = "cand_bachelorD0Sig"
    return _narrow(df, f"{col}[i] > {SENTINEL} && std::abs({col}[i]) > {min_abs_sig}")


def select_bachelor_dedx(df, min_dedx):
    """Harmonic-2 dE/dx floor on the bachelor track.

    Mostly removes candidates whose estimator is unset rather than performing any
    real particle identification, which at a median bachelor pT of 0.52 GeV this
    estimator does not deliver.
    """
    col = "cand_bachelorDedxHarmonic2"
    return _narrow(df, f"{col}[i] > {SENTINEL} && {col}[i] > {min_dedx}")


def select_joint_n_iter(df, arm, max_iter):
    """Upper bound on the joint fit's iteration count.

    Warning, and the reason this is not applied by default: the iteration count is
    a diagnostic of the very fit that produces the mass. Selecting on how hard that
    fit worked selects, through the same fit, on the mass and its resolution, so a
    peak that narrows under this cut has to be shown not to have been sculpted by
    it.
    """
    return _narrow(df, f"{CAND}_{arm.name}Niter[i] <= {max_iter}")


def select_chisq_ndof(df, max_chisq_ndof):
    """Upper bound on the joint fit's chi-square per degree of freedom.

    The normalised variable, not the raw chi-square. `jointCvhNdof` is not
    constant -- it spans 32 to 66 between the 5th and 99th percentiles -- so a raw
    chi-square bound is a cut whose tightness varies by a factor of two across
    candidates, and it varies with the number of hits on the three tracks, hence
    with pT and eta. Requires `define_chisq_ndof`.
    """
    return _narrow(df, f"cand_chisqNdof[i] < {max_chisq_ndof}")


def select_n_legs_refit(df, min_legs):
    """Require this many legs to have used a CVH refit track."""
    return _narrow(df, f"{CAND}_nLegsRefit[i] >= {min_legs}")


# ---------------------------------------------------------------------------
# Candidate reduction
# ---------------------------------------------------------------------------


def select_best_candidate(df, prob_col):
    """Keep the single surviving candidate with the highest vertex probability.

    Candidate columns are reduced to length-1 vectors rather than scalars, so
    downstream code sees a uniform shape.
    """
    df = df.Define(
        "cand_best_idx",
        f"""
        int best = -1;
        double best_prob = -1.;
        for (size_t i = 0; i < cand_pass.size(); ++i) {{
            if (!cand_pass[i]) continue;
            if ({prob_col}[i] > best_prob) {{ best_prob = {prob_col}[i]; best = (int)i; }}
        }}
        return best;
        """,
    )
    df = df.Filter("cand_best_idx >= 0", "one candidate per event")

    for name in df.GetColumnNames():
        col = str(name)
        if not (col.startswith(f"{CAND}_") or col.startswith("cand_")):
            continue
        if col in ("cand_pass", "cand_best_idx"):
            continue
        coltype = df.GetColumnType(col)
        if "RVec" in coltype or "vector<" in coltype:
            df = df.Redefine(
                col,
                f"ROOT::VecOps::Take({col}, ROOT::VecOps::RVec<int>{{cand_best_idx}})",
            )
    return df


def define_scalar(df, column, alias, ctype="double"):
    """Expose the first (only) element of a reduced candidate column as a scalar."""
    return df.Define(alias, f"static_cast<{ctype}>({column}[0])")


# ---------------------------------------------------------------------------
# Truth categorisation
# ---------------------------------------------------------------------------


def _ground_state_expr(var):
    """C++ expression folding an excited b hadron |pdgId| onto its ground state."""
    arms = " ".join(
        f"if (p == {excited}) return {ground};"
        for excited, ground in GROUND_STATE.items()
    )
    return f"[](int p) {{ {arms} return p; }}(std::abs({var}))"


# A pdgId carrying a b (anti)quark: for mesons the 5 sits in the hundreds digit
# (B+ 521, B0 511, Bs 531, Bc 541), for baryons in the thousands (Lambda_b 5122).
#
# The diquark exclusion is not decoration. The producer's ancestor search
# rejects quarks, gluons and protons but not diquarks, which are four-digit
# (1103, 2101, 2203, ...) and therefore pass -- and a diquark common ancestor is
# exactly the vacuous shared-beam-remnant match the exclusion exists to
# prevent. Measured: 71 of 9 748 ancestor matches over 60 files. Requiring a b
# *hadron* here rejects all of them, so the leak is defused downstream of a
# producer that has not been fixed yet.
_IS_B_HADRON = """
inline bool cvhnano_isBHadron(int pdg) {
  const int p = std::abs(pdg);
  if (p < 100) return false;
  const bool diquark = (p >= 1000 && p < 6000 && (p % 100) < 10);
  if (diquark) return false;
  return (p / 100) % 10 == 5 || (p / 1000) % 10 == 5;
}
"""

_declared = False


def _declare_helpers():
    """Declare the C++ helpers once. Redeclaring them is a hard JIT error."""
    global _declared
    if _declared:
        return
    import ROOT

    ROOT.gInterpreter.Declare(_IS_B_HADRON)
    _declared = True


def define_gen_category(df, is_data, dr_max=None):
    """Fill the genCategory axis from the generator parentage of each leg.

    Every simulated bin except the last is a *real b-hadron decay*, and the
    signal bin is an exclusive decay rather than a shared ancestor. Both points
    were measured, not assumed.

    **A shared ancestor is not the decay.** Requiring only that the three legs
    descend from one B+ admits two components, and neither belongs in a signal
    template: the kaon arriving via an excited kaon (K*+ 323, K*+(1410) 10323,
    K2*+ 325), 13.7% of B+-ancestor candidates, and the J/psi arriving via
    charmonium feed-down (chi_c1 20443, psi(2S) 100443, chi_c0 10441), 1.9%.
    Both sit ~200 MeV low because a daughter is missing, and between them they
    put **zero** candidates in 5.24 < m < 5.32. Folding them into the signal
    would add a 15.6% low-mass tail that is not signal, so the chain is walked:
    the bachelor must be a direct daughter of the ancestor, and the muons'
    J/psi must itself be a direct daughter.

    **`dr_max` is gone.** The old definition matched the *candidate direction*
    to the nearest last-copy b hadron and needed a tight angular cut to keep
    coincidental matches out, which cost signal and never separated
    B+ -> J/psi K+ from B+ -> J/psi K*+. Per-leg matching already requires
    dR < 0.02 **and** |dpt|/pt < 0.1 per leg in the producer, so there is no
    candidate-level radius left to choose. The argument is accepted and ignored
    so that a stale command line fails loudly at the parser rather than
    silently changing meaning here.

    The last simulated bin is candidates with no b-hadron common ancestor. It is
    filled and is not a template -- see the module comment in
    `btojpsik_cvhnano_axes` and `datagroups_cvhnano`.
    """
    if is_data:
        return df.Define("cand_genCategory", f"static_cast<double>({GEN_DATA})")

    _declare_helpers()
    return df.Define(
        "cand_genCategory",
        f"""
        const int ap = {CAND}_genAncestorPdgId[0];
        if (ap == 0) return static_cast<double>({GEN_NO_B_ANCESTOR});
        const int g = {_ground_state_expr("ap")};
        if (!cvhnano_isBHadron(g)) return static_cast<double>({GEN_NO_B_ANCESTOR});

        const int nG = static_cast<int>(Gen_pdgId.size());
        const int ai = {CAND}_genAncestorIdx[0];
        const int i0 = {CAND}_leg0GenIdx[0];
        const int i2 = {CAND}_leg2GenIdx[0];

        // Row -1 means unmatched; the producer never forces a nearest neighbour.
        auto mother = [&](int row) {{
            if (row < 0 || row >= nG) return -1;
            const int m = Gen_genPartIdxMother[row];
            return (m >= 0 && m < nG) ? m : -1;
        }};

        // The bachelor is a direct daughter of the ancestor.
        const bool kaonDirect = (i2 >= 0 && ai >= 0 && mother(i2) == ai);
        // The muons come from a J/psi that is itself a direct daughter.
        const int jrow = mother(i0);
        const bool jpsiDirect =
            (jrow >= 0 && ai >= 0 && std::abs(Gen_pdgId[jrow]) == {JPSI_PDGID}
             && mother(jrow) == ai);

        const int l2 = std::abs({CAND}_leg2GenPdgId[0]);
        if (g == {CATEGORY_MOTHER[GEN_SIGNAL]}) {{
            if (l2 != {KAON_PDGID}) return static_cast<double>({GEN_BU_PIX});
            return (kaonDirect && jpsiDirect) ? static_cast<double>({GEN_SIGNAL})
                                              : static_cast<double>({GEN_BU_KX});
        }}
        if (g == {CATEGORY_MOTHER[GEN_BD]}) return static_cast<double>({GEN_BD});
        if (g == {CATEGORY_MOTHER[GEN_BS]}) return static_cast<double>({GEN_BS});
        return static_cast<double>({GEN_OTHER_B});
        """,
    )


def define_bachelor_pt_response(df, is_data):
    """Reconstructed over generated bachelor pT, where the leg matched.

    The direct handle on the bachelor momentum scale, and the one quantity the
    old NanoAOD could not supply at all (gap G2). Unmatched legs write a
    sentinel rather than zero: a spike at zero in a response histogram reads as
    a catastrophic mismeasurement, which is the opposite of "no truth here".
    """
    if is_data:
        return df
    return df.Define(
        "cand_bachelorPtResponse",
        f"""
        ROOT::VecOps::RVec<float> out({CAND}_leg2GenPt.size(), {SENTINEL}f);
        for (size_t i = 0; i < out.size(); ++i) {{
            if ({CAND}_leg2GenPt[i] <= 0.f) continue;
            out[i] = {CAND}_kaonPt[i] / {CAND}_leg2GenPt[i];
        }}
        return out;
        """,
    )


# ---------------------------------------------------------------------------
# Momentum-scale variation inputs
# ---------------------------------------------------------------------------

# Leg order for the scale variations: bachelor first, then the two muons. This
# matches the concatenation order the 2018 histmaker uses, so the per-leg mass
# list [m_K, 0, 0] means the same thing in both channels.
#
# The NanoAOD's own leg indices run the other way round -- leg0/leg1 are the
# muons and leg2 is the bachelor -- so the mapping is stated once here rather
# than being re-derived at each use.
SCALE_VAR_LEGS = [("bachelor", 2), ("mu0", 0), ("mu1", 1)]

# The network's per-muon class tag. 443 is the J/psi-leg class; there is no
# hadron class, so the bachelor is tagged 443 as well, following d0_mass.py.
# This is the assumption the channel rests on and it is measured, not assumed:
# see task 6.5.
MUON_SOURCE_JPSI = 443

# PDG charged-kaon mass [GeV]. Enters the scale variations through the e term,
# which shifts total energy rather than pt, so pt is recomputed through
# 1/cosh(eta) and the shift is asymmetric in down/up.
KAON_MASS_GEV = 0.493677


def define_scale_variation_inputs(
    df, is_data, prefix="scaleVar", legs=None, muon_source=MUON_SOURCE_JPSI
):
    """Three-leg reco and gen kinematics for the A/e/M variation helpers.

    Returns ``(df, prefix)``. Every column is a 3-element RVec in
    ``SCALE_VAR_LEGS`` order, which is what both the ONNX-backed
    ``JpsiCorrectionsUncReweightHelper`` and the analytic
    ``JpsiCorrectionsUncHelperSplines`` expect: they loop over the legs and
    multiply the per-leg alternative weights, so one call covers the candidate.

    Gen kinematics come from indexing ``Gen_*`` with ``leg{N}GenIdx`` rather
    than from new branches -- the NanoAOD stores per-leg gen *pt* but not gen
    eta, phi or charge, and the index is all that is needed.

    Unmatched legs (``GenIdx < 0``) are written with gen == reco, so the
    response residual is exactly 1 and the leg contributes a defined weight.
    That is a formality rather than a policy: measured at kaon pT > 1 GeV,
    every candidate in the exclusive-signal and other-real-b categories has all
    three legs matched, and **no** not-fully-matched candidate is signal. The
    chain walk that assigns the truth category needs three matched legs to reach
    a common ancestor, so a partially matched signal candidate cannot exist.
    ``{prefix}_allMatched`` carries the flag so the claim is counted on the full
    sample rather than asserted from 12 files.

    Data gets no gen columns at all, so the caller must not build variations for
    it -- the variations are a property of the simulated template.
    """
    if is_data:
        return df, None

    legs = SCALE_VAR_LEGS if legs is None else legs

    for var, ctype in [("Pt", "float"), ("Eta", "float"), ("Phi", "float")]:
        df = df.Define(
            f"{prefix}_reco{var}",
            "ROOT::VecOps::RVec<%s>{%s}"
            % (
                ctype,
                ", ".join(
                    f"static_cast<{ctype}>(cand_{tag}{var}[0])" for tag, _ in legs
                ),
            ),
        )
    df = df.Define(
        f"{prefix}_recoCharge",
        "ROOT::VecOps::RVec<int>{%s}"
        % ", ".join(f"static_cast<int>(cand_{tag}Charge[0])" for tag, _ in legs),
    )

    # Gen index per leg, in the same order, guarded against an index past the
    # end of the Gen collection. -1 means "no match".
    df = df.Define(
        f"{prefix}_genIdx",
        """
        ROOT::VecOps::RVec<int> out{%s};
        for (auto &idx : out)
            if (idx >= static_cast<int>(Gen_pt.size())) idx = -1;
        return out;
        """ % ", ".join(f"{CAND}_leg{n}GenIdx[0]" for _, n in legs),
    )
    df = df.Define(
        f"{prefix}_allMatched",
        f"ROOT::VecOps::Min({prefix}_genIdx) >= 0",
    )

    for var, gen_col, ctype in [
        ("Pt", "Gen_pt", "float"),
        ("Eta", "Gen_eta", "float"),
        ("Phi", "Gen_phi", "float"),
        ("Charge", "Gen_charge", "int"),
    ]:
        df = df.Define(
            f"{prefix}_gen{var}",
            f"""
            ROOT::VecOps::RVec<{ctype}> out({prefix}_reco{var}.size());
            for (size_t i = 0; i < out.size(); ++i) {{
                const int idx = {prefix}_genIdx[i];
                out[i] = (idx >= 0) ? static_cast<{ctype}>({gen_col}[idx])
                                    : {prefix}_reco{var}[i];
            }}
            return out;
            """,
        )

    df = df.Define(
        f"{prefix}_muon_source",
        "ROOT::VecOps::RVec<int>{%s}" % ", ".join([str(int(muon_source))] * len(legs)),
    )
    return df, prefix


def define_forced_decay_weight(df, is_data, enable):
    """Undo the generator's forced b-hadron decays.

    Every b hadron in the fragment's `list_forced_decays` was decayed into a
    J/psi mode with total probability 1, and J/psi -> mu mu likewise. The
    simulation therefore over-represents each species by 1/BR, and BR spans
    0.00047 (Xi_b0) to 0.911 (Bc+) -- a factor of 1900 in relative composition
    between the extremes, and 20 between B+ and Lambda_b. Without this weight a
    stacked plot's *shape* is fine per component but its component *fractions*
    are not physical.

    **The species now comes from the muons' parentage, not from the candidate
    direction.** The weight belongs to whichever b hadron produced the J/psi,
    and that is answered exactly by walking mother links up from a matched muon
    leg to the first b hadron. The previous version read the species off
    `genBPdgId`, the nearest last-copy b hadron to the *candidate*, which is the
    right answer for the wrong reason -- it works only because the J/psi carries
    almost all the candidate momentum. That fallback is retained for candidates
    whose muon legs did not match, which is where it is genuinely the only
    information available.

    Deliberately *not* gated on the truth category. The event exists at all only
    because its J/psi passed the generator filter, so a combinatorial candidate
    in the same event needs the same weight. Gating on the category would leave
    the combinatorial bin unweighted and the composition incommensurate, which
    is wrong by orders of magnitude, not a detail.

    Normalised to the B+ so the signal category keeps weight 1: what matters is
    the ratio between species, and an absolute 1e-3 scale only makes yields
    unreadable. Lambda_b then carries 0.048 and Bc+ 51.8.
    """
    if is_data or not enable:
        return df.Define("cand_forcedDecayWeight", "1.0")

    _declare_helpers()
    ref = FORCED_BR[MOTHER_PDGID]
    arms = " ".join(
        f"if (g == {pid}) return {br / ref};" for pid, br in FORCED_BR.items()
    )
    return df.Define(
        "cand_forcedDecayWeight",
        f"""
        const int nG = static_cast<int>(Gen_pdgId.size());
        int species = 0;
        int row = {CAND}_leg0GenIdx[0];
        for (int step = 0; step < 12 && row >= 0 && row < nG; ++step) {{
            const int m = Gen_genPartIdxMother[row];
            if (m < 0 || m >= nG) break;
            if (cvhnano_isBHadron(Gen_pdgId[m])) {{ species = Gen_pdgId[m]; break; }}
            row = m;
        }}
        if (species == 0) species = {CAND}_genBPdgId[0];
        if (species == 0) return 1.0;
        const int g = {_ground_state_expr("species")};
        return [](int g) {{ {arms} return 1.0; }}(g);
        """,
    )


def define_pileup_weight(df, is_data, weights):
    """Reweight the simulation to the data pileup profile.

    Keyed on ``int(Pileup_nTrueInt)`` -- the generated Poisson mean, which is
    what a pileup profile is defined on (gap G4, closed). The previous version
    reweighted the *reconstructed* vertex count, because the NanoAOD carried no
    truth pileup at all. That stopgap removed the observable `nPV` disagreement
    (chi2/bin 201.4 -> 1.34) without claiming to reproduce the true profile,
    and it could not: `nPV` folds in the vertex-reconstruction efficiency, so
    matching it absorbs a reconstruction difference into a pileup weight.

    Now that the true profile is reweightable, the `nPV` agreement becomes a
    *test* rather than something imposed by construction, which is the whole
    point of having asked for the branch.

    Weights outside the tabulated range fall back to 1.0 rather than to the
    nearest tabulated value: extrapolating a ratio measured where there is data
    into a region where there is none invents a number.
    """
    if is_data or not weights:
        return df.Define("cand_pileupWeight", "1.0")

    # Written as a block with an explicit return, not as an immediately-invoked
    # lambda: RDataFrame treats a Define expression *beginning* with `[` as a
    # callable rather than an expression, and the JIT then fails with the
    # unhelpful "cannot form a reference to 'void'".
    arms = "\n".join(
        f"        if (n == {k}) return {v!r};" for k, v in sorted(weights.items())
    )
    return df.Define(
        "cand_pileupWeight",
        f"""
        const int n = static_cast<int>(Pileup_nTrueInt);
{arms}
        return 1.0;
        """,
    )


# ---------------------------------------------------------------------------
# Generator-truth diagnostics
#
# `define_gen_category` can only ever classify a candidate whose legs matched,
# so the population it calls signal is the population that survived a per-leg
# dR < 0.02 AND |dpt|/pt < 0.1 cut in the producer. That second cut is a
# *resolution* cut, and this analysis measures resolution -- so the signal
# template is defined by a requirement on the very quantity being measured.
#
# The two definitions below are what make that testable rather than arguable:
# `define_gen_truth_signal` finds the decay in the Gen table with no reco
# requirement at all, and `define_gen_match_diagnostics` says, for candidates
# that did reconstruct, how many legs the matcher accepted.
# ---------------------------------------------------------------------------


def define_gen_truth_signal(df, is_data):
    """Count B+ -> J/psi(-> mu mu) K+ decays by walking the Gen table.

    The one distribution in the analysis that the gen *matching* cannot bias,
    because no reco object enters it. Must be defined before the trigger and
    candidate filters or it stops being a denominator.

    A decaying B+ is identified by its daughters rather than by status code: a
    documentation copy of the B+ has a B+ daughter and no J/psi, so requiring a
    direct J/psi and a direct charged kaon picks the copy that actually decays
    without needing to know the generator's status conventions.

    Emits, per event:
      ``genTruth_n``        number of such decays
      ``genTruth_kaonPt``   gen pT of each decay's bachelor
      ``genTruth_kaonEta``  gen eta of each decay's bachelor
      ``genTruth_inAcc``    all three final-state legs inside the analysis
                            acceptance, so that acceptance and efficiency can be
                            separated instead of multiplied together
    """
    if is_data:
        return df

    df = df.Define(
        "genTruth_rows",
        f"""
        const int nG = static_cast<int>(Gen_pdgId.size());
        ROOT::VecOps::RVec<int> out;
        for (int i = 0; i < nG; ++i) {{
            if (std::abs(Gen_pdgId[i]) != {MOTHER_PDGID}) continue;
            int kaonRow = -1, jpsiRow = -1;
            for (int j = 0; j < nG; ++j) {{
                if (Gen_genPartIdxMother[j] != i) continue;
                const int p = std::abs(Gen_pdgId[j]);
                if (p == {KAON_PDGID}) kaonRow = j;
                else if (p == {JPSI_PDGID}) jpsiRow = j;
            }}
            if (kaonRow < 0 || jpsiRow < 0) continue;
            // The J/psi must go to two muons: the reco candidate needs both,
            // so a J/psi -> e e decay is not in this denominator.
            int nMu = 0;
            for (int j = 0; j < nG; ++j) {{
                if (Gen_genPartIdxMother[j] == jpsiRow &&
                    std::abs(Gen_pdgId[j]) == 13) ++nMu;
            }}
            if (nMu < 2) continue;
            out.push_back(kaonRow);
            out.push_back(jpsiRow);
        }}
        return out;
        """,
    )
    df = df.Define("genTruth_n", "static_cast<int>(genTruth_rows.size() / 2)")
    df = df.Define(
        "genTruth_kaonPt",
        """
        ROOT::VecOps::RVec<float> out;
        for (size_t i = 0; i + 1 < genTruth_rows.size(); i += 2)
            out.push_back(Gen_pt[genTruth_rows[i]]);
        return out;
        """,
    )
    df = df.Define(
        "genTruth_kaonEta",
        """
        ROOT::VecOps::RVec<float> out;
        for (size_t i = 0; i + 1 < genTruth_rows.size(); i += 2)
            out.push_back(Gen_eta[genTruth_rows[i]]);
        return out;
        """,
    )
    # Acceptance uses the analysis' own thresholds. They are hard-coded rather
    # than threaded through from the parser because this is a *reference*
    # denominator: if it moved with the selection it could not be compared
    # across selections, which is the only thing it is for.
    return df.Define(
        "genTruth_inAcc",
        """
        ROOT::VecOps::RVec<int> out;
        const int nG = static_cast<int>(Gen_pdgId.size());
        for (size_t i = 0; i + 1 < genTruth_rows.size(); i += 2) {
            const int kr = genTruth_rows[i], jr = genTruth_rows[i + 1];
            bool ok = (Gen_pt[kr] > 1.0f && std::abs(Gen_eta[kr]) < 2.4f);
            int nMuOk = 0;
            for (int j = 0; j < nG; ++j) {
                if (Gen_genPartIdxMother[j] != jr) continue;
                if (std::abs(Gen_pdgId[j]) != 13) continue;
                if (Gen_pt[j] > 3.0f && std::abs(Gen_eta[j]) < 2.4f) ++nMuOk;
            }
            out.push_back((ok && nMuOk >= 2) ? 1 : 0);
        }
        return out;
        """,
    )


def define_gen_match_diagnostics(df, is_data):
    """Per-candidate matcher outcome, so the categories can be audited.

    ``cand_nLegsMatched`` is the producer's own count. It is the number that
    decides whether a candidate can have an ancestor at all -- the ancestor
    search runs only when all three legs matched -- and therefore the number
    that decides whether a real signal decay is classified as signal or falls
    into the no-ancestor bin.

    ``cand_bachelorRelDPt`` is the quantity the matcher cut on, kept so that the
    truncation it imposes can be seen rather than inferred.
    """
    if is_data:
        return df

    _declare_helpers()
    df = df.Define(
        "cand_nLegsMatched",
        f"static_cast<double>({CAND}_nLegsGenMatched[0])",
    )
    df = df.Define(
        "cand_bachelorRelDPt",
        f"""
        const float gp = {CAND}_leg2GenPt[0];
        if (gp <= 0.f) return {SENTINEL};
        return static_cast<double>(({CAND}_kaonPt[0] - gp) / gp);
        """,
    )
    df = df.Define(
        "cand_bachelorGenPt",
        f"static_cast<double>({CAND}_leg2GenPt[0])",
    )
    # Why a candidate has no b ancestor, which the category axis cannot say.
    # `define_gen_category` sends three different failures to the same
    # no-ancestor bin, and they call for three different responses:
    #
    #   0  a leg was not matched -- the ancestor search never ran
    #   1  every leg matched but the search found nothing within its depth
    #      budget, which is a property of the SEARCH, not of the event
    #   2  an ancestor was found and it is not a b hadron -- the only one of the
    #      three that is genuinely "not a b decay"
    #   3  a b ancestor was found (so the candidate is not in the no-ancestor
    #      bin at all; kept so the axis sums to the full sample)
    #
    # State 1 is the one worth separating: this NanoAOD matches against the RECO
    # `genParticles` collection, which keeps every documentation copy, while the
    # depth budget was ported from Bmm5, which searches the *pruned* collection
    # with the copies already removed.
    return df.Define(
        "cand_ancestorState",
        f"""
        if ({CAND}_nLegsGenMatched[0] < 3) return 0.0;
        const int ap = {CAND}_genAncestorPdgId[0];
        if (ap == 0) return 1.0;
        const int g = {_ground_state_expr("ap")};
        return cvhnano_isBHadron(g) ? 3.0 : 2.0;
        """,
    )
