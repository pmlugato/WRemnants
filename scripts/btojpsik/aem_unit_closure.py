#!/usr/bin/env python
"""Closure of the A/e/M unit fix in the B+ histmaker (notes/aem_unit_mismatch.md, step 2).

Rerunning ``btojpsik_cvhnano.py`` with ``--fitMuonScaleAndResolution`` changes
only the width of one unit of each scale nuisance. Three things must hold:

1. **Nominal.** Every histogram that is not the A/e/M variation is unchanged,
   to floating-point reassociation (the RDataFrame merges per-slot partial sums
   in a nondeterministic order, so bit-exactness is not available at -j > 1).
2. **Unit ratio.** The induced shift of the mean B mass per unit of each
   nuisance scales, fine kaon-eta bin by fine kaon-eta bin, by the ratio of
   the new width to the old one. The shift is measured on truth-matched signal
   with the kaon in that eta bin, over the fit window.
3. **Linearity.** The reweight is ``exp(log_r)`` with clipping, not a linear
   derivative, and one unit of e at scale 5 is 5e-2 GeV, about 5% of a 1 GeV
   kaon's momentum. So the 5/5/1 histogram against the 1/1/1 one must give 5,
   5 and 1; the deviation, with the up/down asymmetry of the shift, is the
   measured non-linearity per pT band.

Only test 1 is pass/fail. Tests 2 and 3 are measurements and are reported as
such: a reweight is not expected to be exactly linear at low kaon pT.
"""

import argparse
import importlib.util
import json
import os
import sys

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from wums import ioutils

SIGNAL_SAMPLE = "InclusiveBToJpsiX_2016PostVFP"
VAR_HIST = "nominal_muonScaleSyst_responseWeights"
GEN_SIGNAL = 0
PARAMS = ["A", "e", "M"]
# Kaon pT bands of aem_recompute_shift.py, so the slopes compare directly.
PT_BANDS = [(1.0, 2.0), (2.0, 2.8), (2.8, 5.7), (5.7, None)]
MASS_WINDOW = (5.05, 5.45)
N_GROUPS = 24


def load_guard():
    path = os.path.join(
        os.path.dirname(__file__), "..", "rabbit", "check_joint_tensor.py"
    )
    spec = importlib.util.spec_from_file_location("check_joint_tensor", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_OPEN = []  # the histograms are lazy proxies, so their files stay open


def read(path):
    f = h5py.File(path, "r")
    _OPEN.append(f)
    meta = ioutils.pickle_load_h5py(f["meta_info"])
    out = {
        s: ioutils.pickle_load_h5py(f[s])["output"]
        for s in f.keys()
        if s != "meta_info"
    }
    return meta["args"], out


def widths(guard, args):
    """[48, 3] width of one unit, as the joint-fit guard computes it."""
    c = guard.Checker()
    w = guard.implied_widths(c, "", "bplus", args)
    if c.failures:
        raise RuntimeError(f"widths inconsistent with flags: {c.failures}")
    return np.asarray(w)


def nominal_closure(ref, test):
    """Max relative difference of every non-variation histogram."""
    rows = []
    for sample in sorted(set(ref) | set(test)):
        r, t = ref.get(sample, {}), test.get(sample, {})
        for name in sorted(set(r) | set(t)):
            if name == VAR_HIST:
                continue
            if name not in r or name not in t:
                rows.append((sample, name, None, "missing in one input"))
                continue
            try:
                a = np.asarray(r[name].get().values(flow=True))
                b = np.asarray(t[name].get().values(flow=True))
            except AttributeError:
                continue
            if a.shape != b.shape:
                rows.append((sample, name, None, f"shape {a.shape} vs {b.shape}"))
                continue
            both0 = (a == 0) & (b == 0)
            with np.errstate(divide="ignore", invalid="ignore"):
                rel = np.where(
                    both0, 0.0, np.abs(b - a) / np.maximum(np.abs(a), 1e-300)
                )
            rows.append((sample, name, float(np.nanmax(rel)) if rel.size else 0.0, ""))
    return rows


def induced_shifts(h, hnom):
    """Mean-mass shift per unit, [param, fine eta, pT band], with curvature and errors.

    Signal only, kaon in the given fine eta bin, mass in the fit window; the
    nuisance is the eta group containing that fine bin. ``hnom`` is the
    unvaried ``fitbins`` histogram (same pt/eta/charge/mass axes), which gives
    the nominal mean and so the curvature
    (mu_up + mu_dn - 2 mu_0) / (mu_up - mu_dn): zero for a linear response.
    """
    hs = h[{"genCategory": GEN_SIGNAL}]
    hn = hnom[{"genCategory": GEN_SIGNAL}]
    pt_edges = hs.axes["fit_bachelor_pt"].edges
    m_axis = hs.axes["jointCvh_mass"]
    mc = m_axis.centers
    msel = (mc > MASS_WINDOW[0]) & (mc < MASS_WINDOW[1])
    mc = mc[msel]
    # values: pt(+overflow), eta, charge, mass, unc, downUp
    v = np.asarray(hs.values(flow=True))
    var = np.asarray(hs.variances(flow=True))
    # drop flow on all axes except pt overflow (pt has no underflow)
    ax = hs.axes
    def strip(x):
        sl = []
        for a in ax:
            lo = 1 if a.traits.underflow else 0
            hi = (1 if a.traits.overflow else 0)
            if a.name == "fit_bachelor_pt":
                sl.append(slice(lo, None))  # keep overflow in the last band
            else:
                sl.append(slice(lo, x.shape[len(sl)] - hi))
        return x[tuple(sl)]
    v, var = strip(v), strip(var)
    vn = np.asarray(hn.values(flow=True))
    vn = vn[1 if hn.axes[0].traits.underflow else 0 :, :hn.axes[1].size, :hn.axes[2].size, :][..., : hn.axes[3].size]
    vn = vn[:, :, :, msel].sum(axis=2)  # pt, eta, mass
    v = v[:, :, :, msel].sum(axis=2)  # sum charge -> pt, eta, mass, unc, du
    var = var[:, :, :, msel].sum(axis=2)
    npt = v.shape[0]
    pt_lo = np.append(pt_edges[:-1], pt_edges[-1])  # overflow bin starts at the last edge

    neta = v.shape[1]
    nb = len(PT_BANDS)
    shift = np.full((3, neta, nb), np.nan)
    asym = np.full((3, neta, nb), np.nan)
    err = np.full((3, neta, nb), np.nan)
    nsig = np.zeros((neta, nb))
    rebin = neta // N_GROUPS
    for ib, (lo, hi) in enumerate(PT_BANDS):
        psel = (pt_lo[:npt] >= lo - 1e-9) & (
            np.ones(npt, bool) if hi is None else (pt_lo[:npt] < hi - 1e-9)
        )
        vb = v[psel].sum(axis=0)  # eta, mass, unc, du
        wb = var[psel].sum(axis=0)
        nb0 = vn[psel].sum(axis=0)  # eta, mass
        for ieta in range(neta):
            g = ieta // rebin
            for ip in range(3):
                iu = g * 3 + ip
                dn, up = vb[ieta, :, iu, 0], vb[ieta, :, iu, 1]
                n0 = nb0[ieta]
                if up.sum() <= 0 or dn.sum() <= 0 or n0.sum() <= 0:
                    continue
                mu_up = (up * mc).sum() / up.sum()
                mu_dn = (dn * mc).sum() / dn.sum()
                mu_0 = (n0 * mc).sum() / n0.sum()
                shift[ip, ieta, ib] = 1e3 * 0.5 * (mu_up - mu_dn)
                if mu_up != mu_dn:
                    asym[ip, ieta, ib] = (mu_up + mu_dn - 2 * mu_0) / (mu_up - mu_dn)
                # conservative: up and down treated as independent samples
                n_eff = up.sum() ** 2 / max(wb[ieta, :, iu, 1].sum(), 1e-300)
                rms = np.sqrt(((up * (mc - mu_up) ** 2).sum()) / up.sum())
                err[ip, ieta, ib] = 1e3 * 0.5 * np.sqrt(2.0) * rms / np.sqrt(max(n_eff, 1))
            nsig[ieta, ib] = nb0[ieta].sum()
    return shift, asym, err, nsig


def band_label(ib):
    lo, hi = PT_BANDS[ib]
    return f"{lo:g}-{hi:g} GeV" if hi else f">{lo:g} GeV"


def plot_ratio(eta_c, s_ref, e_ref, s_test, e_test, expect, labels, title, path):
    fig, (ax, rx) = plt.subplots(
        2, 1, figsize=(8, 6), sharex=True, gridspec_kw={"height_ratios": [2.2, 1]}
    )
    ax.errorbar(eta_c, s_ref, e_ref, fmt="o", ms=3, label=labels[0])
    ax.errorbar(eta_c, s_test, e_test, fmt="s", ms=3, label=labels[1])
    ax.axhline(0, color="gray", lw=0.5)
    ax.set_ylabel("induced mean m(B) shift\nper unit [MeV]")
    ax.set_yscale("symlog", linthresh=max(1e-3, np.nanmin(np.abs(s_ref)) if np.any(np.isfinite(s_ref)) else 1e-3))
    ax.legend(fontsize=8)
    ax.set_title(title, fontsize=10)
    with np.errstate(divide="ignore", invalid="ignore"):
        r = s_test / s_ref / expect
        re = np.abs(r) * np.sqrt((e_test / s_test) ** 2 + (e_ref / s_ref) ** 2)
    rx.errorbar(eta_c, r, re, fmt="o", ms=3, color="k")
    rx.axhline(1, color="r", lw=0.8)
    rx.set_ylim(0.5, 1.5)
    rx.set_ylabel("measured /\nexpected ratio")
    rx.set_xlabel("kaon eta (fine bin centre)")
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    fig.savefig(path.replace(".png", ".pdf"))
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ref", required=True, help="old histogram (massfit widths)")
    p.add_argument("--test", required=True, help="new histogram (fitted-scale widths)")
    p.add_argument("--test1", help="new histogram at scale 1/1/1, for the linearity test")
    p.add_argument("--outdir", required=True)
    p.add_argument("--nominalTol", type=float, default=1e-12)
    a = p.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    guard = load_guard()

    ref_args, ref = read(a.ref)
    test_args, test = read(a.test)
    w_ref, w_test = widths(guard, ref_args), widths(guard, test_args)
    summary = {"ref": a.ref, "test": a.test}

    # ---- 1. nominal
    rows = nominal_closure(ref, test)
    bad = [r for r in rows if r[2] is None or r[2] > a.nominalTol]
    print(f"\n[1] nominal closure: {len(rows)} histograms compared, tolerance {a.nominalTol:g}")
    for s, n, m, why in rows:
        flag = "ok " if (m is not None and m <= a.nominalTol) else "BAD"
        print(f"   {flag} {s:<34s} {n:<36s} max rel diff {m if m is not None else why}")
    summary["nominal_pass"] = not bad
    summary["nominal_rows"] = [(s, n, m, w) for s, n, m, w in rows]
    print(f"   -> {'PASS' if not bad else 'FAIL'} ({len(bad)} differing)")

    # ---- 2. unit ratio
    s_ref, as_ref, e_ref, n_ref = induced_shifts(ref[SIGNAL_SAMPLE][VAR_HIST].get(), ref[SIGNAL_SAMPLE]["fitbins"].get())
    s_test, as_test, e_test, _ = induced_shifts(test[SIGNAL_SAMPLE][VAR_HIST].get(), test[SIGNAL_SAMPLE]["fitbins"].get())
    expect = (w_test / w_ref).T  # [param, eta]
    eta_edges = np.linspace(-2.4, 2.4, w_ref.shape[0] + 1)
    eta_c = 0.5 * (eta_edges[1:] + eta_edges[:-1])
    print("\n[2] unit ratio: measured shift ratio / expected width ratio, signal-weighted over eta")
    summary["unit_ratio"] = {}
    for ip, par in enumerate(PARAMS):
        for ib in range(len(PT_BANDS)):
            r = s_test[ip, :, ib] / s_ref[ip, :, ib] / expect[ip]
            wgt = n_ref[:, ib] * np.isfinite(r)
            mean = np.nansum(r * wgt) / wgt.sum() if wgt.sum() else np.nan
            med = np.nanmedian(r)
            summary["unit_ratio"][f"{par}_{ib}"] = dict(
                weighted_mean=float(mean), median=float(med),
                shift_ref_median=float(np.nanmedian(s_ref[ip, :, ib])),
                shift_test_median=float(np.nanmedian(s_test[ip, :, ib])),
            )
            print(
                f"   {par} pT {band_label(ib):<12s} ratio/expected: weighted mean {mean:.4f}, "
                f"median {med:.4f}   |   shift per unit, eta median: old {np.nanmedian(s_ref[ip,:,ib]):+.4f} "
                f"new {np.nanmedian(s_test[ip,:,ib]):+.3f} MeV"
            )
            plot_ratio(
                eta_c, s_ref[ip, :, ib], e_ref[ip, :, ib], s_test[ip, :, ib],
                e_test[ip, :, ib], expect[ip],
                ["old unit (massfit widths)", "new unit (fitted-scale widths)"],
                f"{par}: kaon pT {band_label(ib)}; ratio panel = (new/old) / (width ratio)",
                os.path.join(a.outdir, f"unit_ratio_{par}_pt{ib}.png"),
            )

    # ---- 3. linearity
    if a.test1:
        t1_args, t1 = read(a.test1)
        w1 = widths(guard, t1_args)
        s1, as1, e1, n1 = induced_shifts(t1[SIGNAL_SAMPLE][VAR_HIST].get(), t1[SIGNAL_SAMPLE]["fitbins"].get())
        expect5 = (w_test / w1).T
        print("\n[3] linearity: shift(scale 5/5/1) / shift(scale 1/1/1) / width ratio; curvature (mu_up+mu_dn-2mu_0)/(mu_up-mu_dn)")
        summary["linearity"] = {}
        for ip, par in enumerate(PARAMS):
            for ib in range(len(PT_BANDS)):
                r = s_test[ip, :, ib] / s1[ip, :, ib] / expect5[ip]
                wgt = n1[:, ib] * np.isfinite(r)
                mean = np.nansum(r * wgt) / wgt.sum() if wgt.sum() else np.nan
                summary["linearity"][f"{par}_{ib}"] = dict(
                    weighted_mean=float(mean),
                    width_ratio=float(np.nanmedian(expect5[ip])),
                    curvature_scale5=float(np.nanmedian(as_test[ip, :, ib])),
                    curvature_scale1=float(np.nanmedian(as1[ip, :, ib])),
                )
                print(
                    f"   {par} pT {band_label(ib):<12s} expected x{np.nanmedian(expect5[ip]):.0f}: "
                    f"measured/expected {mean:.4f}   curvature scale5 {np.nanmedian(as_test[ip,:,ib]):+.2e} "
                    f"scale1 {np.nanmedian(as1[ip,:,ib]):+.2e}"
                )
                plot_ratio(
                    eta_c, s1[ip, :, ib], e1[ip, :, ib], s_test[ip, :, ib],
                    e_test[ip, :, ib], expect5[ip],
                    ["scale 1/1/1", "scale 5/5/1"],
                    f"{par}: kaon pT {band_label(ib)}; ratio panel = (5/5/1 / 1/1/1) / width ratio",
                    os.path.join(a.outdir, f"linearity_{par}_pt{ib}.png"),
                )

    with open(os.path.join(a.outdir, "closure_summary.json"), "w") as f:
        json.dump(summary, f, indent=1, default=str)
    print(f"\nwrote {a.outdir}")
    return 0 if summary["nominal_pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
