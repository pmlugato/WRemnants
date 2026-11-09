"""Deck: trustworthy truth in the regenerated CVH NanoAOD, and the prefits it buys.

Companion to make_cvhnano_tensor_deck.py, which is the record of the tensor and
the first fits. This deck covers what changed when the NanoAOD was regenerated
with generator parentage, and what the truth categorisation says about the
background -- the question the earlier fits could not answer.

Figures come from the central plotting chain:
  scripts/plotting/btojpsik_cvhnano_plots.sh <histmaker.hdf5>
which runs makeDataMCStackPlot.py --noData for the simulation-only stacks and
btojpsik_cvhnano_prefit.py for the data-side prefits.
"""

import json
import os
import sys

sys.path.insert(0, "/home/submit/pmlugato/.claude/skills/slides")
from deck import Deck

ASSETS = "/home/submit/pmlugato/public_html/slides/assets_cvhnano_truth"
HIST = "/ceph/submit/data/user/p/pmlugato/mz/histograms"
PLOTS = "/home/submit/pmlugato/public_html/mz"

# Measured numbers are read from the JSON the measurement wrote, not retyped.
# A missing key is a build failure rather than a blank on a slide -- a deck that
# silently prints nothing where a number should be is worse than one that fails.
with open(f"{HIST}/truth_validation_truthv2.json") as fh:
    V = json.load(fh)
with open(f"{PLOTS}/btojpsik_cvhnano_truthv2_prefit/prefit_summary.json") as fh:
    P = json.load(fh)
with open(f"{PLOTS}/cvhnano_kpt1_prefit/prefit_summary.json") as fh:
    PK = json.load(fh)
with open(f"{ASSETS}/widths.json") as fh:
    W = json.load(fh)
# Kaon pT > 1 GeV, the selection the 2018 NanoAOD applied at production and ours
# does not. Same 2 520-file inputs, so every comparison against the baseline is
# like for like.
with open(f"{HIST}/truth_validation_kpt1.json") as fh:
    VK = json.load(fh)
with open(f"{ASSETS}/kpt1/widths.json") as fh:
    WK = json.load(fh)
with open(f"{ASSETS}/id_summary.json") as fh:
    ID = json.load(fh)
with open(f"{HIST}/dimuon_cut_scan.json") as fh:
    SCAN = {r["name"]: r for r in json.load(fh)}
with open(f"{HIST}/truth_validation_kaonid15.json") as fh:
    VID = json.load(fh)


def sc(name, field, fmt="{:.3f}"):
    return fmt.format(SCAN[name][field])


def rej(var, sl, at=90):
    return f"{ID['rejection'][var][sl][f'rej_at_{at}']:.3f}"
# The 2018 cross-check. Optional only because the two 2018 histmaker passes take
# about an hour each; when the file is absent the slide is skipped and a warning
# is printed, rather than a slide being built with holes in it.
W2018 = None
_w2018_path = f"{HIST}/width_comparison_2018.json"
if os.path.exists(_w2018_path):
    with open(_w2018_path) as fh:
        W2018 = json.load(fh)
else:
    print(f"WARNING: {_w2018_path} absent -- the 2018 comparison slide is SKIPPED")


def comp(name, field):
    return V["composition"][name][field]


def pct(name):
    return f"{comp(name, 'pct_of_peak'):.2f}%"


def n(name):
    return f"{comp(name, 'total'):,.0f}".replace(",", " ")


def shape(label):
    return f"{V['shape_tests'][label]['chi2_per_bin']:.2f}"


def kcomp(name, field):
    return VK["composition"][name][field]


def kpct(name):
    return f"{kcomp(name, 'pct_of_peak'):.2f}%"


def kbkg(model, field):
    v = WK["background_models"][model][field]
    return f"{v:,.0f}".replace(",", " ") if field == "signal" else f"{v:.2f}"


# Every real b-hadron decay except the signal, as a share of the peak window.
SUM_REAL_B_PCT = "{:.2f}%".format(
    sum(comp(k, "pct_of_peak")
        for k in ("BuJpsiKX", "BuJpsiPiX", "Bd", "Bs", "OtherB"))
)

def MEV(arm, side):
    """side: 'data', 'mc_total' or 'mc_signal_only'."""
    return f"{W['widths'][arm]['lin'][side + '_mev']:.2f} MeV"


def RATIO(arm):
    return f"{W['widths'][arm]['lin']['ratio_to_total']:.3f}"


def BKG(model, field):
    v = W["background_models"][model][field]
    return f"{v:,.0f}".replace(",", " ") if field == "signal" else f"{v:.2f}"


_PF = P["BuJpsiK_jointCvhMass"]["prefit"]
PREFIT_SIGNAL = f"{_PF['signal']:,.0f}".replace(",", " ")
PREFIT_CHI2 = f"{_PF['chi2ndf']:.2f}"
_PFK = PK["BuJpsiK_jointCvhMass"]["prefit"]
PREFIT_KPT1_SIGNAL = f"{_PFK['signal']:,.0f}".replace(",", " ")
PREFIT_KPT1_CHI2 = f"{_PFK['chi2ndf']:.2f}"

d = Deck(assets=ASSETS)

d.title_slide(
    "Truth we can trust, and the prefits it buys",
    "The background model and the kaon pT floor are one problem — and 2018 had the cut we lack",
    "Pietro Lugato  ·  18 August 2026  ·  mZ / CVH NanoAOD",
)

d.bullets(
    "Where this goes",
    [
        "What the regenerated NanoAOD added, and which gaps it closed",
        "How the truth match works now, and why a shared ancestor is not enough",
        "The categories: every simulated component is a real b-hadron decay",
        "The composition, and the answer to “is there a peaking background?”",
        "The prefit, and what its residual actually says — including a claim I retracted",
    ],
)

# --------------------------------------------------------------------------
# what changed
# --------------------------------------------------------------------------

d.table(
    "Eight gaps closed by re-running the NanoAOD step",
    [
        ["gap", "what arrived", "evidence"],
        ["G1", "Gen_genPartIdxMother", "0.135% at −1, i.e. record roots only"],
        ["G2", "per-leg gen match", "kaon reco/gen pT 1.00001, rms 0.0183 (1.00017 with pT > 1)"],
        ["G4", "Pileup_nTrueInt", "mean 21.82 sim vs 24.31 data — a real offset"],
        ["G6", "track errors, hit detail", "finite and positive for 100% of tracks"],
        ["G8", "3D pointing angle", "alphaBS and lxy checked: already clean"],
        ["G11", "unconstrained m(μμ)", "median 3.09385 vs PDG 3.09690"],
        ["G13", "PV-referenced d0", "φ modulation 0.187 → 0.002 cm, factor 107"],
    ],
    lead="No AlCaReco reproduction: everything came from files already on disk.",
    widths=[0.10, 0.32, 0.58],
    note="G7 closes as NOT A GAP: candidate muons are already 100% global-and-tracker, and "
    "softMuonMva cannot be evaluated at all — innerTrack() is unavailable for 100% of "
    "persisted muons. That retires the hypothesis that a missing muon variable is why the "
    "2018 selection could not be reproduced.",
)

d.bullets(
    "How the truth match works now",
    [
        "Each leg is matched to a final-state gen particle: dR < 0.02 AND |Δpt|/pt < 0.1",
        "Charge must agree; an unmatched leg stays −1, never a forced nearest neighbour",
        "The three matched legs' mother links are walked to a common ancestor",
        "The search permutes generation depth per leg — muons sit two below the B, the kaon one",
        ("Quarks, gluons and protons are skipped, or every pair shares a beam proton", 1),
    ],
)

d.bullets(
    "Why the old definition could not do this",
    [
        "It matched the candidate DIRECTION to the nearest last-copy b hadron",
        "That needed a tight dR cut to keep coincidental matches out, and still let them in",
        "It cannot distinguish B⁺ → J/ψ K⁺ from B⁺ → J/ψ K*⁺ — the direction is the same",
        "Graded against parentage: 78.3% pure at 98.7% efficiency",
        ("Usable as a fallback; not good enough for a template", 1),
    ],
)

d.table(
    "A shared ancestor is not the decay",
    [
        ["B⁺-ancestor candidates, kaon bachelor", "share", "median m", "in 5.24–5.32"],
        ["kaon and J/ψ both DIRECT daughters", "84.4%", "5.2807", "87.9%"],
        ["kaon via K*⁺, K*⁺(1410), K₂*⁺", "13.7%", "5.0751", "0.0%"],
        ["J/ψ via χ_c1, ψ(2S), χ_c0", "1.9%", "5.0104", "0.0%"],
    ],
    lead="This is why the chain is walked rather than just the ancestor read off.",
    widths=[0.44, 0.14, 0.18, 0.24],
    note="Both rejected components sit ~200 MeV low, because a daughter is missing, and neither "
    "puts a single candidate in the peak window. Folding them into the signal would add a 15.6% "
    "low-mass tail that is not signal. The muons' mother is the J/ψ for 1535 of 1535, so the "
    "chain walk is not fragile.",
)

d.table(
    "Why they sit low — and why one does not",
    [
        ["category", "dominant decay", "missing", "endpoint", "median"],
        ["BuJpsiKX", "B⁺ → J/ψ K*⁺", "π⁰", "5.144", "5.068"],
        ["Bd", "B⁰ → J/ψ K*⁰", "π⁻", "5.140", "5.094"],
        ["BuJpsiPiX", "B⁺ → J/ψ π⁺", "nothing", "—", "5.389"],
    ],
    lead="The B mass is the same for all of them. Two mechanisms move it, in opposite directions.",
    widths=[0.20, 0.28, 0.16, 0.18, 0.18],
    note="A missing daughter X gives m(J/ψK) ≤ m_B − m_X, a hard endpoint 140 MeV low. "
    "B⁺ → J/ψ π⁺ misses nothing: the pion gets the KAON mass, so the mass comes out "
    "+84 MeV HIGH — calculated 5.3634, measured 5.3890, on the upper shoulder.",
)

d.table(
    "The categories: every one is a real b-hadron decay",
    [
        ["bin", "group", "definition"],
        ["0", "BuJpsiK", "ancestor B⁺; J/ψ and kaon both direct; bachelor a kaon"],
        ["1", "BuJpsiKX", "ancestor B⁺, kaon bachelor, not the direct chain"],
        ["2", "BuJpsiPiX", "ancestor B⁺, non-kaon bachelor (98% pion)"],
        ["3", "Bd", "ancestor B⁰"],
        ["4", "Bs", "ancestor B⁰ₛ"],
        ["5", "OtherB", "any other b hadron (Λ_b, Ξ_b, B_c, …)"],
        ["6", "NoBAncestor", "no b-hadron ancestor — filled, never plotted or fitted"],
    ],
    lead="Excited b hadrons are folded onto their ground state before the species test.",
    widths=[0.08, 0.20, 0.72],
    note="Bin 6 is combinatorial, and combinatorial is modelled in data. This simulation omits "
    "the prompt and fake components entirely (G9), so a template for one third of the "
    "background would read as a measurement of all of it.",
)

# --------------------------------------------------------------------------
# what the composition says
# --------------------------------------------------------------------------

d.table(
    "There is no peaking background",
    [
        ["category", "candidates", "% of peak window"],
        ["B⁺ → J/ψ K⁺", n("BuJpsiK"), pct("BuJpsiK")],
        ["B⁺ → J/ψ K⁺X", n("BuJpsiKX"), pct("BuJpsiKX")],
        ["B⁺ → J/ψ π⁺X", n("BuJpsiPiX"), pct("BuJpsiPiX")],
        ["B⁰", n("Bd"), pct("Bd")],
        ["B⁰ₛ", n("Bs"), pct("Bs")],
        ["other b hadron", n("OtherB"), pct("OtherB")],
        ["no b ancestor — not a template", n("NoBAncestor"), pct("NoBAncestor")],
    ],
    lead="Joint-arm mass, gated, window 5.24 < m < 5.32. Full simulation.",
    widths=[0.44, 0.28, 0.28],
    note=f"Every real-b background together is {SUM_REAL_B_PCT} of the window, i.e. "
    f"{V['real_b_background_over_signal'] * 100:.1f}% of the signal — and each peaks 100–200 MeV "
    "BELOW the B mass, not under it. One smooth analytic shape plus the signal template is the "
    "right model. The earlier fits did not fail for want of a peaking component.",
)

d.figure(
    "The simulated stack, at its own yields",
    "stack_jointCvh_mass.png",
    lead="Simulation only. Not normalised to data — it is not a model of the data.",
)

d.figure(
    "The same stack, with kaon pT > 1 GeV",
    "stack_jointCvh_mass_kpt1.png",
    lead="Same components, same treatment — only the bachelor pT floor moves from 0.15 to 1 GeV.",
)

d.figure(
    "The same stack, with everything else it reconstructs",
    "stack_residual.png",
    lead="Log scale, and the combinatorial included: the residual is what the real-b decays are not.",
)

d.figure(
    "The residual, with kaon pT > 1 GeV",
    "stack_residual_kpt1.png",
    lead="Same log scale. The grey maximum under the peak is what the pT floor removes.",
)

d.figure_row(
    "The same two, combinatorial at the bottom of the stack",
    ["stack_residual_combofirst_truthv2.png", "stack_residual_combofirst_kpt1.png"],
    lead="Left: baseline. Right: kaon pT > 1 GeV. Everything real-b is the thin skin on top.",
)

d.figure(
    "The prefit with an exponential background",
    "prefit_jointCvh_mass.png",
    lead=f"χ²/ndf {PREFIT_CHI2}. Too high at the peak, too low on both shoulders — read on.",
)

d.figure(
    "Swap only the background, and the fit closes",
    "fig_bkgmodel.png",
    lead="Same data, same template. Green is the true simulated shape; red is an exponential.",
)

d.table(
    "With our selection, the background model is the problem",
    [
        ["background", "χ²/ndf", "fitted signal"],
        ["true simulated non-signal shape", BKG("true_shape", "chi2ndf"), BKG("true_shape", "signal")],
        ["exponential", BKG("exponential", "chi2ndf"), BKG("exponential", "signal")],
        ["quadratic", BKG("quadratic", "chi2ndf"), BKG("quadratic", "signal")],
    ],
    lead="54 bins over 5.15 < m < 5.42, 1 452 352 data candidates.",
    widths=[0.48, 0.20, 0.32],
    note="The exponential inflates the signal yield 2.8×. A cubic gives 85.22, i.e. no better "
    "than the quadratic — so this is not one more polynomial order away. rabbit's available "
    "models are monotonic, and the background has a maximum: that is the gap.",
)

d.figure(
    "The same fit, drawn as a subtraction",
    "fig_width.png",
    lead="A picture of the χ²/ndf 5.62 fit — NOT independent evidence. See the caveat below.",
)

d.bullets(
    "Why that last plot proves less than it looks like",
    [
        "Data was fitted as a·S + b·B with BOTH shapes fixed from simulation",
        "Only the two normalisations floated — the bump was not free",
        "But if the fit describes the data, the residual must resemble S by construction",
        "So the honest statement is the χ²: 5.62 with fixed shapes against 745.79",
        ("And 5.62 is not 1 — there is residual shape disagreement, just 130× less", 1),
    ],
)

d.bullets(
    "A claim I retracted, and the control that caught it",
    [
        "First version subtracted a straight LINE from the data and got a ratio of 1.71",
        "The simulated background is a broad maximum on the B mass — a line cannot remove it",
        "So what was measured was the peak plus a residual bump, on the data side only",
        "Control: treat the simulated TOTAL identically. It gives 24.62 vs data 26.50, ratio 1.08",
        ("The inflation was in the method, and it inflates both sides equally", 1),
    ],
)

d.table(
    "Widths, data against the simulated total, one treatment",
    [
        ["arm", "data σ", "simulated total σ", "ratio"],
        ["joint N-body CVH", MEV("jointCvh", "data"), MEV("jointCvh", "mc_total"), RATIO("jointCvh")],
        ["KVF on CVH tracks", MEV("kvfCvh", "data"), MEV("kvfCvh", "mc_total"), RATIO("kvfCvh")],
        ["KVF on raw tracks", MEV("kvfRaw", "data"), MEV("kvfRaw", "mc_total"), RATIO("kvfRaw")],
        ["raw, pre-fit", MEV("raw", "data"), MEV("raw", "mc_total"), RATIO("raw")],
        ["dimuon, CVH tracks", MEV("dimuonCvh", "data"), MEV("dimuonCvh", "mc_total"), RATIO("dimuonCvh")],
        ["dimuon, raw tracks", MEV("dimuonRaw", "data"), MEV("dimuonRaw", "mc_total"), RATIO("dimuonRaw")],
    ],
    lead="Gaussian on a line, 5.19 < m < 5.37, identical on both sides.",
    widths=[0.34, 0.20, 0.26, 0.20],
    note="These are still line-subtracted, so the absolute values are inflated — the point is "
    "that the inflation is common. The dimuon, which shares no bachelor track, agrees to 1.02.",
)

# --------------------------------------------------------------------------
# the other fixes, each with the disagreement it was supposed to move
# --------------------------------------------------------------------------

d.table(
    "Data/simulation shape agreement, χ² per bin",
    [
        ["observable", "χ²/bin", "reads as"],
        ["bachelor dxy (origin-referenced)", shape("bachelor dxy (origin-referenced, the defect)"), "the defect, kept visible"],
        ["bachelor d0 (PV-referenced)", shape("bachelor d0 (PV-referenced, the fix)"), "the fix, at source"],
        ["bachelor d0 significance", shape("bachelor d0 significance"), "new discriminant (G6)"],
        ["bachelor nValidHits", shape("bachelor nValidHits"), "the standing disagreement"],
        ["nPV", shape("nPV"), "now a test, not a fit"],
        ["unconstrained m(μμ)", shape("unconstrained dimuon mass"), "the muon scale handle"],
    ],
    lead="Shapes only, and comparable WITHIN this run — not against earlier runs.",
    widths=[0.42, 0.16, 0.42],
    note="χ²/bin grows with statistics, so Phase 2 numbers are not commensurate. Within-run: "
    "the fraction within 0.05 cm of the vertex goes 0.286/0.175 (data/sim) → 0.893/0.884. "
    "True-pileup reweighting does NOT fix nPV: that gap is vertex reconstruction.",
)

if W2018:
    rows = [["input", "data σ", "sim total σ", "ratio"]]
    for r in W2018:
        rows.append([
            r["name"],
            f"{r['data_sigma_mev']:.2f} MeV",
            f"{r['mc_total_sigma_mev']:.2f} MeV",
            f"{r['ratio_to_total']:.3f}",
        ])
    d.table(
        "2018 against us, one window and one method",
        rows,
        lead="Common window |m − 5.279| < 74 MeV, the widest the 2018 mass axis allows.",
        widths=[0.40, 0.20, 0.22, 0.18],
        note="The 2018 simulation is EXCLUSIVE signal with no background, so its ratio is an "
        "UPPER bound on a fair comparison. Ours compares data to the simulated total, "
        "with background on both sides.",
    )

    d.figure_row(
        "The 2018 plot, with and without those two cuts",
        ["plot2018_full.png", "plot2018_nodisc.png"],
        lead="Left: full selection. Right: softMVA and the BDT dropped. Kaon |η| < 2.4 in both.",
    )

    d.bullets(
        "What those two cuts do, and what they do not",
        [
            "Dropping both takes the 2018 data from 1 868 087 to 5 179 229 candidates",
            "So they remove 64% of what survives everything else — they work",
            "But the peak SHAPE does not move: core 10.07 → 10.06 MeV, extra smear 3.5 → 3.5",
            "What grows is the pedestal, and the 2018 simulation has no background at all",
            ("So they buy purity. Data IS wider than simulation — just not because of them", 1),
        ],
    )

    d.table(
        "The selection, 2018 against now — and what proxies what",
        [
            ["stage", "2018", "now", "handle / proxy"],
            ["muon pT, |η|", "4 GeV, 1.4", "3 GeV, 2.4 — applied", "Track via mu{0,1}TrackIdx"],
            ["muon soft ID, BMM BDT", "0.45, 0.10", "no counterpart", "innerTrack absent for 100%"],
            ["dimuon pT, αBS, vtxProb, sl3d", "7, 0.4, 0.1, 4", "all available, all OFF",
             "borrowed from the kvfCvh arm"],
            ["B vertex probability", "0.3", "available, OFF", "χ²/ndof upper tail — the proxy"],
            ["B fit quality", "—", "gate only, ratio OFF", "edmRef, χ², ndof, σ(m)"],
            ["kaon pT", "1–8 GeV", "> 0.1, no ceiling", "1 GeV was a PRODUCTION cut"],
            ["kaon |η|", "1.4", "2.4", "candidate columns"],
        ],
        lead="Every 2018 handle exists as an argument. Every optional one defaults to OFF.",
        widths=[0.28, 0.16, 0.24, 0.32],
        note="The joint arm has one common vertex for all legs, so its B probability is the "
        "χ²/ndof upper tail and its dimuon block is borrowed. The J/ψ subsystem χ²/ndof — the "
        "noted proxy for route-1 vertex probabilities — is NOT in the nano.",
    )

    d.bullets(
        "A caution on every width number in this deck",
        [
            "Four treatments of the same data gave 1.71, 1.076, 1.101 and 0.957",
            "Each differs only in how the background under the peak was modelled",
            "The smearing fit is better posed but still reaches χ²/ndf 72.8 on 2018",
            "So the data/simulation width difference is real and NOT yet reliably measured",
            ("Settle it inside the fit, with a floated resolution scale, not with subtractions", 1),
        ],
    )

d.table(
    "Kaon pT > 1 GeV: a cut 2018 had and we do not",
    [
        ["", "baseline", "kaon pT > 1 GeV"],
        ["signal, % of peak window", pct("BuJpsiK"), kpct("BuJpsiK")],
        ["combinatorial, % of window", pct("NoBAncestor"), kpct("NoBAncestor")],
        ["signal candidates kept", f"{comp('BuJpsiK','in_peak'):,.0f}".replace(",", " "),
         f"{kcomp('BuJpsiK','in_peak'):,.0f}".replace(",", " ")],
        ["real-b background / signal",
         f"{V['real_b_background_over_signal']:.4f}", f"{VK['real_b_background_over_signal']:.4f}"],
    ],
    lead="Applied at production in 2018; our nano keeps kaons down to 0.15 GeV.",
    widths=[0.44, 0.28, 0.28],
    note="96.8% of the signal kept for 2.7× the purity. A production-level cut, not an "
    "analysis one, so every selection comparison so far was blind to it.",
)

d.table(
    "And it largely fixes the background model too",
    [
        ["background", "χ²/ndf, baseline", "χ²/ndf, kaon pT > 1"],
        ["true simulated shape", BKG("true_shape", "chi2ndf"), kbkg("true_shape", "chi2ndf")],
        ["exponential", BKG("exponential", "chi2ndf"), kbkg("exponential", "chi2ndf")],
        ["quadratic", BKG("quadratic", "chi2ndf"), kbkg("quadratic", "chi2ndf")],
    ],
    lead="Same fit, same window, same template treatment — only the kaon pT floor changes.",
    widths=[0.40, 0.30, 0.30],
    note="The exponential's yield bias collapses from 2.8× to 1.11×, and the quadratic now fits "
    "BETTER than the simulated shape, whose own statistical noise becomes the limit. "
    "A monotonic background is adequate once the soft-kaon combinatorial is gone.",
)

d.figure(
    "The same test, with the soft kaons gone",
    "fig_bkgmodel_kpt1.png",
    lead="Both models now track the data. The exponential's residual is ±15%, not ±30%.",
)

d.figure(
    "The prefit that an exponential can almost carry",
    "prefit_kpt1.png",
    lead=f"Kaon pT > 1 GeV: {PREFIT_KPT1_SIGNAL} signal, χ²/ndf {PREFIT_KPT1_CHI2} "
    f"against {PREFIT_CHI2} before.",
)

d.table(
    "Soft kaons give the BEST B mass resolution",
    [
        ["kaon pT [GeV]", "signal n", "core width", "w95"],
        ["0.15 – 0.50", "1 182", "7.64 MeV", "42.6"],
        ["0.50 – 1.00", "5 183", "7.84 MeV", "44.0"],
        ["1.00 – 1.50", "6 495", "8.03 MeV", "43.6"],
        ["2.00 – 3.00", "11 648", "10.41 MeV", "50.7"],
        ["5.00 – 10.0", "16 411", "12.44 MeV", "57.6"],
    ],
    lead="Simulated signal only. Half-width containing 50% of the peak, monotonic in pT.",
    widths=[0.28, 0.22, 0.26, 0.24],
    note="A factor 1.6 across the range: σ(p)/p is roughly flat, so a soft kaon carries a "
    "smaller ABSOLUTE momentum error into the B mass. Cutting kaon pT > 1 GeV therefore "
    "discards the best-measured signal — the template got WIDER, 10.93 → 11.55 MeV.",
)

d.table(
    "What the low-kaon-pT background actually is",
    [
        ["bachelor truth", "kaon pT < 1", "kaon pT > 1"],
        ["a real pion, not from a b hadron", "66.5%", "61.5%"],
        ["unmatched — fake, or >10% mismeasured", "23.9%", "12.5%"],
        ["a real kaon, not from a b hadron", "6.9%", "15.7%"],
        ["a real proton", "1.8%", "8.2%"],
        ["a real track straight from a b hadron", "0.3%", "1.4%"],
    ],
    lead="No-ancestor candidates, per-leg truth. n = 44 314 and 8 033.",
    widths=[0.48, 0.26, 0.26],
    note="It is a real J/ψ from a b hadron plus a SOFT PROMPT TRACK from the primary vertex — "
    "two thirds pions, essentially none from a b decay. That names the discriminant: "
    "displacement, not particle identity.",
)

d.bullets(
    "The gen match biases the quantity we are measuring",
    [
        "Matching requires dR < 0.02 AND |Δpt|/pt_gen < 0.1 on every leg",
        "So a badly measured kaon fails to match and its candidate leaves the signal category",
        "Measured: 3.6% of B+-like candidates lose their bachelor leg that way",
        "Those have p16–p84 of 5.126–5.411 against 5.250–5.298 for the matched ones",
        ("The template's TAILS are truncated by the truth definition, not by the detector", 1),
    ],
)

d.bullets(
    "Where we actually stand",
    [
        "No peaking template is needed — real-b backgrounds are 0.5% of the window",
        "The soft-kaon combinatorial is what an analytic shape cannot describe",
        "But cutting it away costs the best-measured signal, so a pT cut is the wrong tool",
        "Data is wider than simulation in BOTH samples — 2018 needs ~3.5 MeV extra smearing too",
        ("And part of that may be the gen match truncating the template's tails", 1),
    ],
)

d.figure(
    "The background is prompt. The signal is not.",
    "fig_id_d0sig.png",
    lead="Combinatorial spikes at zero significance; signal is broad. Every pT slice.",
)

d.table(
    "Background rejected at fixed signal efficiency",
    [
        ["variable", "0–0.4 GeV", "0.4–0.7", "0.7–1", "1.5–2.5"],
        ["d0 significance", rej("id_d0Sig", "0.0-0.4"), rej("id_d0Sig", "0.4-0.7"),
         rej("id_d0Sig", "0.7-1.0"), rej("id_d0Sig", "1.5-2.5")],
        ["d0 itself", rej("id_d0", "0.0-0.4"), rej("id_d0", "0.4-0.7"),
         rej("id_d0", "0.7-1.0"), rej("id_d0", "1.5-2.5")],
        ["dz significance", rej("id_dzPVSig", "0.0-0.4"), rej("id_dzPVSig", "0.4-0.7"),
         rej("id_dzPVSig", "0.7-1.0"), rej("id_dzPVSig", "1.5-2.5")],
        ["dE/dx", rej("id_dedx", "0.0-0.4"), rej("id_dedx", "0.4-0.7"),
         rej("id_dedx", "0.7-1.0"), rej("id_dedx", "1.5-2.5")],
    ],
    lead="Rejection at 90% signal efficiency, per bachelor pT slice.",
    widths=[0.28, 0.18, 0.18, 0.18, 0.18],
    note="Significance beats raw d0 everywhere, because the resolution varies with "
    "momentum and a fixed distance is a different cut at every pT. Transverse beats "
    "longitudinal. dE/dx is the weakest by a factor of four.",
)

d.table(
    "One cut, and it beats the pT floor",
    [
        ["|d0/σ| >", "signal eff", "background kept", "implied peak purity"],
        ["1.0", "0.863", "0.209", "53.3%"],
        ["1.5", "0.800", "0.089", "71.3%"],
        ["2.0", "0.740", "0.042", "83.0%"],
        ["3.0", "0.630", "0.015", "92.2%"],
        ["kaon pT > 1 GeV", "—", "—", "57.9%"],
    ],
    lead="Integrated over all bachelor pT, so the soft kaons are kept.",
    widths=[0.24, 0.22, 0.26, 0.28],
    note="Purity is IMPLIED: efficiencies are measured over the whole stored mass range and "
    "applied to the peak-window composition, so it assumes the cut does not sculpt the mass. "
    "A run with |d0/σ| > 1.5 applied is going now to test exactly that.",
)

d.figure(
    "dE/dx: real separation, too little of it",
    "fig_id_dedx.png",
    lead="Soft kaons do ionise more — the solid blue tail above 5. It is not enough.",
)

d.bullets(
    "Why dE/dx does not work, even where it should",
    [
        "Separation exists only below 0.4 GeV and only in the high tail",
        "Rejection there is 0.108 at 90% signal efficiency, against 0.370 for d0 significance",
        "And the estimator is MISSING for 26.9% of signal tracks at 0–0.4 GeV",
        "11% missing at 0.4–0.7, 5% at 0.7–1 — absent exactly where it was wanted",
        ("Worth having tested; not worth building a selection on", 1),
    ],
)

d.bullets(
    "Four questions this does not yet answer",
    [
        "How much of the data/simulation width gap is detector and how much is the truth match?",
        "What functional form describes a J/ψ plus a soft prompt track? It has a maximum",
        "Can a displacement-based kaon ID reach the purity a pT cut reaches, at full efficiency?",
        "Does dE/dx separate K from π where it matters — below ~0.7 GeV/c?",
        ("Each is answerable with what is now on disk, except the first", 1),
    ],
)

d.bullets(
    "Proposed order of work",
    [
        "1. Re-match on dR alone — configured, dR tightened 0.02 → 0.01, validating now",
        "2. Kaon ID on |d0/σ|: measured, and it beats the pT floor at every working point",
        "3. dE/dx: tested and dropped, on the numbers above",
        "4. Check the d0 cut does not sculpt the mass — the run is going",
        ("5. Then the differential fit, with the soft kaons still in it", 1),
    ],
)

# --------------------------------------------------------------------------
# the fair comparison with 2018: does the dimuon selection do anything?
# --------------------------------------------------------------------------

d.table(
    "The 2018 dimuon and muon cuts, each alone and all together",
    [
        ["configuration", "sig eff", "purity", "χ² exp", "bias"],
        ["baseline", "1.000", sc("baseline", "purity_peak"),
         sc("baseline", "chi2_exp", "{:.0f}"), sc("baseline", "bias_exp", "{:.2f}×")],
        ["dimuon pT > 7", sc("dimuon pT > 7", "eff"), sc("dimuon pT > 7", "purity_peak"),
         sc("dimuon pT > 7", "chi2_exp", "{:.0f}"), sc("dimuon pT > 7", "bias_exp", "{:.2f}×")],
        ["dimuon vtxProb > 0.1", sc("dimuon vtxProb > 0.1", "eff"),
         sc("dimuon vtxProb > 0.1", "purity_peak"),
         sc("dimuon vtxProb > 0.1", "chi2_exp", "{:.0f}"),
         sc("dimuon vtxProb > 0.1", "bias_exp", "{:.2f}×")],
        ["dimuon αBS < 0.4", sc("dimuon alphaBS < 0.4", "eff"),
         sc("dimuon alphaBS < 0.4", "purity_peak"),
         sc("dimuon alphaBS < 0.4", "chi2_exp", "{:.0f}"),
         sc("dimuon alphaBS < 0.4", "bias_exp", "{:.2f}×")],
        ["B vtxProb > 0.1", sc("B vtxProb > 0.1", "eff"), sc("B vtxProb > 0.1", "purity_peak"),
         sc("B vtxProb > 0.1", "chi2_exp", "{:.0f}"), sc("B vtxProb > 0.1", "bias_exp", "{:.2f}×")],
        ["muon pT > 4", sc("muon pT > 4", "eff"), sc("muon pT > 4", "purity_peak"),
         sc("muon pT > 4", "chi2_exp", "{:.0f}"), sc("muon pT > 4", "bias_exp", "{:.2f}×")],
        ["all five together", sc("all five together", "eff"),
         sc("all five together", "purity_peak"),
         sc("all five together", "chi2_exp", "{:.0f}"),
         sc("all five together", "bias_exp", "{:.2f}×")],
    ],
    lead="χ² exp is the exponential-background fit to data; bias is its yield against the true shape.",
    widths=[0.28, 0.16, 0.16, 0.18, 0.22],
    note="Three are no-ops. B vtxProb makes the background WORSE — max/min 2.67 → 13.95, the "
    "fit-quality family sculpting toward the B mass again. αBS is the only one that helps and "
    "it costs 65% of the signal. All five together: 29% efficiency for purity 0.376.",
)

d.figure(
    "The only dimuon variable with any power",
    "fig_dimuon_kvfCvh_dimuon_alphaBS.png",
    lead="The combinatorial contains a REAL J/ψ from a real b hadron, so it looks like signal.",
)

d.table(
    "Signal efficiency against kaon pT — flat means unbiased",
    [
        ["configuration", "0–0.5", "0.5–1", "1–1.5", "1.5–2.5", "4–20"],
        ["dimuon vtxProb > 0.1", "0.979", "0.978", "0.978", "0.977", "0.978"],
        ["dimuon αBS < 0.4", "0.259", "0.269", "0.286", "0.301", "0.372"],
        ["B vtxProb > 0.1", "0.868", "0.863", "0.857", "0.857", "0.853"],
        ["all five together", "0.226", "0.225", "0.238", "0.253", "0.311"],
        ["|d0/σ| > 1.5", "0.385*", "0.633*", "0.843", "0.884", "0.823"],
    ],
    lead="The dimuon cuts are unbiased in kaon pT. The kaon ID is not.",
    widths=[0.28, 0.14, 0.14, 0.14, 0.15, 0.15],
    note="*|d0/σ| slices are 0–0.4 and 0.4–0.7. σ(d0) grows at low pT through scattering, "
    "so a FIXED significance threshold is tighter there.",
)

d.table(
    "Where the signal and the background actually are",
    [
        ["kaon pT [GeV]", "% of signal", "% of background", "S/B"],
        ["0.0–0.4", "0.51", "33.08", "0.0012"],
        ["0.4–0.7", "2.58", "36.15", "0.0056"],
        ["0.7–1.0", "4.02", "15.04", "0.0210"],
        ["1.5–2.5", "16.16", "4.78", "0.2658"],
        ["4.0–20", "48.42", "0.46", "8.2046"],
    ],
    lead="Below 1 GeV: 7.1% of the signal and 84.3% of the background.",
    widths=[0.26, 0.24, 0.26, 0.24],
    note="This is why the kaon pT floor worked and the dimuon cuts cannot. It also sets the "
    "real question: is the A/e leverage of 7% of the signal, sitting where S/B is 0.001–0.02, "
    "worth the background it brings? That depends on ∂m/∂A and ∂m/∂e, not on yield.",
)

d.figure_row(
    "By eye: all five 2018 cuts, against the kaon ID",
    ["stack_residual_s18all.png", "stack_residual_kaonid.png"],
    lead="Left: the 2018 dimuon block. Right: |d0/σ| > 1.5. Combinatorial at the bottom, log scale.",
)

d.table(
    "The kaon ID, measured — and it beats everything else tried",
    [
        ["", "baseline", "kaon pT > 1", "|d0/σ| > 1.5", "2018 dimuon block"],
        ["peak purity", "0.214", "0.579", "0.677", "0.376"],
        ["signal efficiency", "1.000", "0.968†", "0.858", "0.290"],
        ["exponential χ²/ndf", "745.8", "41.2", "25.4", "273.6"],
        ["quadratic χ²/ndf", "85.5", "8.5", "6.8", "58.1"],
        ["yield bias", "2.8×", "1.11×", "1.073×", "2.15×"],
    ],
    lead="†of candidates above 1 GeV only — it discards the soft kaons entirely.",
    widths=[0.24, 0.18, 0.18, 0.18, 0.22],
    note="The kaon ID does not sculpt the signal: template width 15.37 → 15.50 MeV. Its cost is "
    "the pT bias above, so the working point is not final — a per-pT threshold at constant "
    "efficiency removes the bias but reaches only 45% implied purity.",
)

d.bullets(
    "What this settles, and the one question left",
    [
        "2018's advantage was NOT its dimuon selection — three of those cuts are no-ops here",
        "The B vertex probability actively hurts: it sculpts the background toward the B mass",
        "Every useful handle is on the bachelor track, because the background has a real J/ψ",
        "Displacement is the best of them: purity 0.677, exponential χ²/ndf 25.4, bias 1.073×",
        ("Left open: keep the soft kaons and model their background, or cut and lose the lever", 1),
    ],
)

d.section("Backup")

d.bullets(
    "Three things worth knowing",
    [
        "Data is at 21 806 of 22 328 files; the rest were still being produced",
        "The producer's ancestor search does not exclude diquarks — four-digit, so they pass",
        ("71 of 9 748 matches; defused here by requiring a b hadron, fix it upstream", 1),
        "Simulation is scaled by a PLACEHOLDER cross section, so its y axis is arbitrary",
        ("Composition and shape are unaffected; only the axis numbers are meaningless", 1),
    ],
)

d.bullets(
    "The commands",
    [
        "scripts/plotting/btojpsik_cvhnano_pileup.py — nTrueInt weights from the standard profile",
        "scripts/histmakers/btojpsik_cvhnano.py --allaxes --fitBinsHist --pileupWeights …",
        "scripts/plotting/btojpsik_cvhnano_plots.sh — stacks (--noData) then the data prefit",
        "scripts/btojpsik/cvhnano_truth_validation.py — every number on the tables above",
        ("One histmaker command, both datasets, one file: --filterProcs is variadic", 1),
    ],
)

d.save("260818_cvhnano_truth_and_prefits.pptx")
