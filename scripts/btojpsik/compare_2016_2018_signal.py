#!/usr/bin/env python3
"""2016postVFP (CVH NanoAOD) against 2018 (BMM NanoAOD), signal only, with ratios.

Both signal definitions are gen-level: 2016 is the `genCategory == BuJpsiK` slice
of the inclusive B -> J/psi X sample; 2018 is the exclusive BuToJpsiK production.

The two samples do NOT share a selection, and the mismatch is not cosmetic:
2018 caps the kaon at |eta| < 1.4 and 8 GeV where the CVH nano runs to 2.4 and
10. Every distribution is therefore shown twice -- once over the CVH selection
and once restricted to the 2018 acceptance -- because resolution degrades with
eta and an unmatched comparison would attribute that to the fit.

    python scripts/btojpsik/compare_2016_2018_signal.py
"""

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from wremnants.postprocessing.rabbit_cvhnano_helpers import open_results

H = "/ceph/submit/data/user/p/pmlugato/mz/histograms/"
OUT = os.path.expanduser("~/public_html/mz/btojpsik_cvhnano")
GEN_SIGNAL = 0
ETA_2018 = 1.4  # the 2018 kaon eta axis stops here
PT_2018 = 8.0  # and the 2018 analysis capped the kaon at 8 GeV


def vals(h):
    a = np.asarray(h.view())
    if a.dtype.names:
        return a["value"].astype(float), a["variance"].astype(float)
    return a.astype(float), a.astype(float)


def rebin(y, v, n):
    k = (len(y) // n) * n
    return y[:k].reshape(-1, n).sum(1), v[:k].reshape(-1, n).sum(1)


def load():
    """(edges, y, var) per distribution per era, on matched binning."""
    out = {}
    with open_results(H + "btojpsik_cvhnano_fullstats.hdf5") as r:
        h = r["InclusiveBToJpsiX_2016PostVFP"]["output"]["fitbins"].get()
        pt_e = np.asarray(h.axes["fit_bachelor_pt"].edges)
        eta_e = np.asarray(h.axes["fit_bachelor_eta"].edges)
        m_e = np.asarray(h.axes["jointCvh_mass"].edges)
        y, v = vals(h)
        y, v = y[..., GEN_SIGNAL], v[..., GEN_SIGNAL]  # pt, eta, charge, mass
        out["16"] = dict(pt_e=pt_e, eta_e=eta_e, m_e=m_e, y=y, v=v)
    with open_results(H + "btojpsik_sel2018full.hdf5") as r:
        o = r["BuToJpsiK_2018"]["output"]
        for key, name in (
            ("m", "nominal_bkmm_jpsimc_mass"),
            ("pt", "nominal_bkmm_jpsimc_kaon1pt"),
            ("eta", "nominal_bkmm_jpsimc_kaon1eta"),
        ):
            hh = o[name].get()
            yy, vv = vals(hh)
            out.setdefault("18", {})[key] = (np.asarray(hh.axes[0].edges), yy, vv)
    return out


def panel(ax_top, ax_rat, edges, a, b, la, lb, xlabel, logy=False):
    """Area-normalised overlay with a ratio panel. Poisson errors, propagated."""
    c = 0.5 * (edges[:-1] + edges[1:])
    w = np.diff(edges)
    out = []
    for (y, v), col, lab in ((a, "C0", la), (b, "C3", lb)):
        n = y.sum()
        d, e = y / n / w, np.sqrt(v) / n / w
        ax_top.errorbar(c, d, yerr=e, fmt="o", ms=3, lw=1, color=col, label=lab)
        out.append((d, e))
    (d1, e1), (d2, e2) = out
    ok = (d1 > 0) & (d2 > 0)
    r = np.where(ok, d2 / np.where(d1 > 0, d1, 1), np.nan)
    re = np.where(ok, np.abs(r) * np.sqrt((e1 / np.where(d1 > 0, d1, 1)) ** 2
                                          + (e2 / np.where(d2 > 0, d2, 1)) ** 2), np.nan)
    ax_rat.axhline(1.0, color="grey", lw=0.8)
    ax_rat.errorbar(c, r, yerr=re, fmt="k.", ms=4, lw=0.9)
    ax_top.set_ylabel("normalised")
    if logy:
        ax_top.set_yscale("log")
    ax_top.legend(frameon=False, fontsize=9)
    ax_rat.set_xlabel(xlabel)
    ax_rat.set_ylabel("2018 / 2016")
    ax_rat.set_xlim(edges[0], edges[-1])
    return r, re


def gauss_sigma(edges, y):
    """Core width from a Gaussian + linear fit, scanned. Returns MeV."""
    c = 0.5 * (edges[:-1] + edges[1:])
    m0 = c[np.argmax(y)]
    best = None
    for s in np.linspace(0.004, 0.040, 361):
        g = np.exp(-0.5 * ((c - m0) / s) ** 2)
        A = np.vstack([g, np.ones_like(c), c - m0]).T
        p, *_ = np.linalg.lstsq(A, y, rcond=None)
        res = float(np.sum((A @ p - y) ** 2))
        if best is None or res < best[0]:
            best = (res, s)
    return 1000 * best[1]


def fig(name, edges, a, b, la, lb, xlabel, title_note="", logy=False):
    f, (t, r) = plt.subplots(
        2, 1, figsize=(10.4, 5.0), sharex=True,
        gridspec_kw=dict(height_ratios=[2.3, 1], hspace=0.06),
    )
    panel(t, r, edges, a, b, la, lb, xlabel, logy=logy)
    p = os.path.join(OUT, name)
    f.savefig(p, dpi=150, bbox_inches="tight")
    plt.close(f)
    print("MARKER wrote", p, title_note)


def main():
    d = load()
    s16 = d["16"]
    pt_e, eta_e, m_e = s16["pt_e"], s16["eta_e"], s16["m_e"]
    y16, v16 = s16["y"], s16["v"]

    pt_c = 0.5 * (pt_e[:-1] + pt_e[1:])
    eta_c = 0.5 * (eta_e[:-1] + eta_e[1:])
    m_c = 0.5 * (m_e[:-1] + m_e[1:])

    accept = {
        "cvh": np.ones_like(eta_c, bool),
        "2018acc": np.abs(eta_c) < ETA_2018,
    }
    ptcut = {"cvh": pt_c > 1.0, "2018acc": (pt_c > 1.0) & (pt_c < PT_2018)}

    for tag in ("cvh", "2018acc"):
        sel = np.ix_(ptcut[tag], accept[tag])
        ysel = y16[sel[0], sel[1]]  # pt, eta, charge, mass
        vsel = v16[sel[0], sel[1]]

        # --- mass, on the common 5.2-5.4 window at 10 MeV
        m16 = ysel.sum(axis=(0, 1, 2))
        mv16 = vsel.sum(axis=(0, 1, 2))
        keep = (m_c > 5.2) & (m_c < 5.4)
        m16, mv16 = m16[keep], mv16[keep]
        e16 = m_e[np.r_[np.where(keep)[0], np.where(keep)[0][-1] + 1]]
        m16, mv16 = rebin(m16, mv16, 2)  # 5 -> 10 MeV
        e16r = e16[::2]

        e18, y18, v18 = d["18"]["m"]
        y18r, v18r = rebin(y18, v18, 5)  # 2 -> 10 MeV
        e18r = e18[::5]
        assert np.allclose(e16r, e18r), (e16r[:3], e18r[:3])

        s_16 = gauss_sigma(e16r, m16)
        s_18 = gauss_sigma(e18r, y18r)
        print(f"MARKER [{tag}] core sigma  2016 {s_16:.2f} MeV   2018 {s_18:.2f} MeV"
              f"   ratio {s_16/s_18:.3f}   (2016 N={m16.sum():,.0f}, 2018 N={y18r.sum():,.0f})")
        fig(f"cmp1618_mass_{tag}.png", e16r, (m16, mv16), (y18r, v18r),
            "2016postVFP, CVH joint 3-body", "2018, BMM",
            "m(J/psi K)  [GeV]", logy=True)

        # --- kaon pT, common 0.1 GeV bins over 1-8
        p16 = ysel.sum(axis=(1, 2, 3))
        pv16 = vsel.sum(axis=(1, 2, 3))
        pe18, py18, pv18 = d["18"]["pt"]
        pc18 = 0.5 * (pe18[:-1] + pe18[1:])
        k18 = (pc18 > 1.0) & (pc18 < PT_2018)
        k16 = (pt_c[ptcut[tag]] > 1.0) & (pt_c[ptcut[tag]] < PT_2018)
        edges = np.arange(1.0, PT_2018 + 1e-9, 0.1)
        fig(f"cmp1618_kaonpt_{tag}.png", edges, (p16[k16], pv16[k16]),
            (py18[k18], pv18[k18]),
            "2016postVFP, CVH joint 3-body", "2018, BMM",
            "kaon $p_T$  [GeV]")

        # --- kaon |eta|, common 0.1 bins over 0-1.4
        g16 = ysel.sum(axis=(0, 2, 3))
        gv16 = vsel.sum(axis=(0, 2, 3))
        ec = eta_c[accept[tag]]
        ge18, gy18, gv18 = d["18"]["eta"]
        gc18 = 0.5 * (ge18[:-1] + ge18[1:])
        edges = np.arange(0.0, ETA_2018 + 1e-9, 0.1)
        f16 = np.histogram(np.abs(ec), bins=edges, weights=g16)[0]
        fv16 = np.histogram(np.abs(ec), bins=edges, weights=gv16)[0]
        f18 = np.histogram(np.abs(gc18), bins=edges, weights=gy18)[0]
        fv18 = np.histogram(np.abs(gc18), bins=edges, weights=gv18)[0]
        fig(f"cmp1618_kaoneta_{tag}.png", edges, (f16, fv16), (f18, fv18),
            "2016postVFP, CVH joint 3-body", "2018, BMM",
            "kaon $|\\eta|$")


if __name__ == "__main__":
    main()
