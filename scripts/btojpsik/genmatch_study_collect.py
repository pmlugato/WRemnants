#!/usr/bin/env python3
"""Gen-matching study, collection step (worklist 3.1-3.4).

Starts from the GENERATOR side, so the result does not depend on the producer's
own verdict: every exclusive B+ -> J/psi(-> mu mu) K+ decay in the simulation is
found from the full `Gen` record (with `genPartIdxMother`), and each of its three
final-state legs is confronted with the reconstructed objects the way
JpsiXKinematicFitProducer does it -- same charge, dR and |dpt|/pt_gen -- but
with the thresholds left free, so the pass/fail decision becomes an analysis
choice and can be plotted against gen pT.

Per leg (kaon, mu+, mu-) the nearest same-charge reco object is recorded with
its dR and (pt_reco - pt_gen)/pt_gen, together with whether ANY same-charge
object passes the producer's cuts (0.02, 0.1), which separates "the nearest one
fails" from "a different one would have passed" (mis-association).

Candidate level: the "truth candidate" is the reconstructed BuJpsiK whose three
track indices are the nearest tracks to the three gen legs (within dR < 0.1, no
pt requirement). Its producer verdict, gate variables and joint-CVH mass are
stored so the fraction of true signal that the producer's thresholds discard,
and its mass distribution, can be shown per gen kaon pT (3.4).

Counts are unweighted: the forced-decay weight is a per-species constant, so
within the exclusive B+ -> J/psi K sample it cancels in every ratio shown.

Output: one .npz with flat per-leg and per-candidate arrays.
"""
import argparse
import glob
import json
import os
import sys
from multiprocessing import Pool

import numpy as np
import uproot

PDG_B = 521
PDG_JPSI = 443
PDG_K = 321
PDG_MU = 13
PDG_GAMMA = 22
SENT = -99.0

# Producer thresholds, for the "would the producer have matched it" flags.
PROD_DR = 0.02
PROD_DPT = 0.1
# Loose association used to define the truth candidate: nearest same-charge
# track within this cone, no momentum requirement.
LOOSE_DR = 0.1

GEN_COLS = ["Gen_pdgId", "Gen_genPartIdxMother", "Gen_status", "Gen_pt", "Gen_eta",
            "Gen_phi", "Gen_charge"]
TRK_COLS = ["Track_pt", "Track_eta", "Track_phi", "Track_charge", "Track_muonIdx",
            "Track_nValidHits", "Track_d0", "Track_d0Err", "Track_ptErr", "Track_normChi2",
            "Track_dxy", "Track_dz"]
MU_COLS = ["Muon_pt", "Muon_eta", "Muon_phi", "Muon_charge", "Muon_isGlobal", "Muon_isTracker", "Muon_nMatches"]
CAND_COLS = ["BuJpsiK_kaonTrackIdx", "BuJpsiK_mu0TrackIdx", "BuJpsiK_mu1TrackIdx",
             "BuJpsiK_leg0GenIdx", "BuJpsiK_leg1GenIdx", "BuJpsiK_leg2GenIdx",
             "BuJpsiK_leg0GenMotherPdgId", "BuJpsiK_leg1GenMotherPdgId",
             "BuJpsiK_leg2GenMotherPdgId", "BuJpsiK_nLegsGenMatched",
             "BuJpsiK_genAncestorPdgId", "BuJpsiK_kaonPt", "BuJpsiK_kaonEta",
             "BuJpsiK_jointCvhOk", "BuJpsiK_jointCvhEdmRef", "BuJpsiK_jointCvhChisq",
             "BuJpsiK_jointCvhNdof", "BuJpsiK_jointCvhMassErr", "BuJpsiK_jointCvhMass",
             "BuJpsiK_kvfRawMass", "BuJpsiK_leg2GenDR", "BuJpsiK_leg2GenPt"]
COLS = GEN_COLS + TRK_COLS + MU_COLS + CAND_COLS + ["event"]


def dphi(a, b):
    d = a - b
    return (d + np.pi) % (2 * np.pi) - np.pi


def delta_r(eta1, phi1, eta2, phi2):
    return np.hypot(eta1 - eta2, dphi(phi1, phi2))


def follow_copies(pid, mo, idx, code):
    """Follow idx down its own-flavour copy chain (e.g. J/psi -> J/psi radiative)."""
    for _ in range(10):
        dau = np.where(mo == idx)[0]
        same = dau[pid[dau] == code]
        if len(same) == 0:
            return idx, dau
        idx = same[0]
    return idx, np.where(mo == idx)[0]


def exclusive_decays(pid, mo, st):
    """Yield (kaon_idx, muplus_idx, muminus_idx, b_idx) for exclusive B+ -> J/psi K."""
    out = []
    for b in np.where(np.abs(pid) == PDG_B)[0]:
        dau = np.where(mo == b)[0]
        dpid = pid[dau]
        if np.any(np.abs(dpid) == PDG_B):
            continue  # not the last copy
        non_gamma = dau[dpid != PDG_GAMMA]
        if len(non_gamma) != 2:
            continue
        s = set(np.abs(pid[non_gamma]).tolist())
        if s != {PDG_JPSI, PDG_K}:
            continue
        k = non_gamma[np.abs(pid[non_gamma]) == PDG_K][0]
        j = non_gamma[pid[non_gamma] == PDG_JPSI][0]
        if st[k] != 1:
            continue
        j, jd = follow_copies(pid, mo, j, PDG_JPSI)
        mus = jd[np.abs(pid[jd]) == PDG_MU]
        if len(mus) != 2 or not np.all(st[mus] == 1):
            continue
        mup = mus[pid[mus] == -PDG_MU]  # mu+ has pdgId -13
        mum = mus[pid[mus] == PDG_MU]
        if len(mup) != 1 or len(mum) != 1:
            continue
        out.append((k, mup[0], mum[0], b))
    return out


def match_leg(gpt, geta, gphi, gq, rpt, reta, rphi, rq):
    """Confront one gen leg with a reco collection. Returns a dict of scalars."""
    sel = rq == gq
    res = dict(n_same_charge=int(sel.sum()), near_idx=-1, near_dr=SENT, near_dpt=SENT,
               near_hits=-1, any_pass=0, n_pass=0, n_in_loose=0, pass_dpt_of_nearest_in_dr=SENT,
               best_dr_passing_dpt=SENT, best_dpt_within_dr=SENT)
    if not sel.any():
        return res
    idx = np.where(sel)[0]
    dr = delta_r(geta, gphi, reta[idx], rphi[idx])
    dpt = (rpt[idx] - gpt) / gpt
    i = np.argmin(dr)
    res["near_idx"] = int(idx[i])
    res["near_dr"] = float(dr[i])
    res["near_dpt"] = float(dpt[i])
    passing = (dr < PROD_DR) & (np.abs(dpt) < PROD_DPT)
    res["any_pass"] = int(passing.any())
    res["n_pass"] = int(passing.sum())
    res["n_in_loose"] = int((dr < LOOSE_DR).sum())
    # Best dR among objects that satisfy the pt cut alone, and best |dpt| among
    # objects that satisfy the dR cut alone: the two one-sided views of a failure.
    okpt = np.abs(dpt) < PROD_DPT
    if okpt.any():
        res["best_dr_passing_dpt"] = float(dr[okpt].min())
    okdr = dr < PROD_DR
    if okdr.any():
        res["best_dpt_within_dr"] = float(dpt[okdr][np.argmin(np.abs(dpt[okdr]))])
    return res


def straight_line_doca(pt1, eta1, phi1, dxy1, dz1, pt2, eta2, phi2, dxy2, dz2):
    """The producer's trackTrackDCA, rebuilt from nano track parameters.

    reco::Track reference point = PCA to the beam line, so (vx, vy, vz) =
    (-dxy sin(phi), dxy cos(phi), dz) with dxy referenced to the origin, which is
    what Track_dxy/Track_dz store (gap G13). Directions from (pt, eta, phi)."""
    def line(pt, eta, phi, dxy, dz):
        p = np.array([-dxy * np.sin(phi), dxy * np.cos(phi), dz])
        u = np.array([np.cos(phi), np.sin(phi), np.sinh(eta)]); u /= np.linalg.norm(u)
        return p, u
    p1, u1 = line(pt1, eta1, phi1, dxy1, dz1); p2, u2 = line(pt2, eta2, phi2, dxy2, dz2)
    c = np.cross(u1, u2); cn = np.linalg.norm(c); d = p2 - p1
    if cn < 1e-12:
        return float(np.linalg.norm(d - np.dot(d, u1) * u1))
    return float(abs(np.dot(d, c)) / cn)


def process_file(path):
    sidecar = os.path.splitext(path)[0] + ".json"
    try:
        with open(sidecar) as fh:
            meta = json.load(fh)
        if meta.get("exit_reason") != "ok":
            return None
    except OSError:
        pass
    try:
        t = uproot.open(path)["Events"]
        a = t.arrays(COLS, library="np")
    except Exception as exc:  # truncated stubs: open the tree, do not count the file
        sys.stderr.write(f"skip {path}: {exc}\n")
        return None

    legs = {"K": [], "mup": [], "mum": []}
    cands = []
    nev = len(a["event"])
    for i in range(nev):
        pid = a["Gen_pdgId"][i]; mo = a["Gen_genPartIdxMother"][i]; st = a["Gen_status"][i]
        gpt = a["Gen_pt"][i]; geta = a["Gen_eta"][i]; gphi = a["Gen_phi"][i]; gq = a["Gen_charge"][i]
        decays = exclusive_decays(pid, mo, st)
        if not decays:
            continue
        tpt = a["Track_pt"][i]; teta = a["Track_eta"][i]; tphi = a["Track_phi"][i]
        tq = a["Track_charge"][i]; tmu = a["Track_muonIdx"][i]; thits = a["Track_nValidHits"][i]
        td0 = a["Track_d0"][i]; td0e = a["Track_d0Err"][i]; tpte = a["Track_ptErr"][i]; tnchi = a["Track_normChi2"][i]
        tdxy = a["Track_dxy"][i]; tdz = a["Track_dz"][i]
        mpt = a["Muon_pt"][i]; meta_ = a["Muon_eta"][i]; mphi = a["Muon_phi"][i]; mq = a["Muon_charge"][i]
        mglob = a["Muon_isGlobal"][i]; mtrk = a["Muon_isTracker"][i]; mnm = a["Muon_nMatches"][i]
        ktrk = a["BuJpsiK_kaonTrackIdx"][i]; m0trk = a["BuJpsiK_mu0TrackIdx"][i]; m1trk = a["BuJpsiK_mu1TrackIdx"][i]
        for (k, mup, mum, b) in decays:
            kin = {}
            near_trk = {}
            for name, g in (("K", k), ("mup", mup), ("mum", mum)):
                r = match_leg(gpt[g], geta[g], gphi[g], gq[g], tpt, teta, tphi, tq)
                r.update(near_d0=SENT, near_d0err=SENT, near_pterr=SENT, near_normchi2=SENT)
                if r["near_idx"] >= 0:
                    ni = r["near_idx"]
                    r["near_hits"] = int(thits[ni]); r["near_d0"] = float(td0[ni]); r["near_d0err"] = float(td0e[ni])
                    r["near_pterr"] = float(tpte[ni]); r["near_normchi2"] = float(tnchi[ni])
                r.update(gen_pt=float(gpt[g]), gen_eta=float(geta[g]), gen_phi=float(gphi[g]),
                         b_pt=float(gpt[b]), kaon_gen_pt=float(gpt[k]), event=int(a["event"][i]))
                # muons additionally against the Muon table (what the J/psi legs are built from)
                if name != "K":
                    rm = match_leg(gpt[g], geta[g], gphi[g], gq[g], mpt, meta_, mphi, mq)
                    r.update({f"muon_{kk}": vv for kk, vv in rm.items()})
                    # tight-selector proxy for the nearest Muon (the B+ channel J/psi uses GoodMuons:
                    # global & tracker & nMatches > 1, plus global-track terms not in the nano)
                    mi = rm["near_idx"]
                    r["muon_tight"] = int(mi >= 0 and rm["near_dr"] < 0.1 and bool(mglob[mi]) and bool(mtrk[mi]) and mnm[mi] > 1)
                legs[name].append(r)
                near_trk[name] = r["near_idx"] if (r["near_idx"] >= 0 and r["near_dr"] < LOOSE_DR) else -1
            # truth candidate: all three loose-nearest tracks present and one candidate uses exactly them
            c = dict(kaon_gen_pt=float(gpt[k]), kaon_gen_eta=float(geta[k]), b_pt=float(gpt[b]),
                     have_tracks=int(all(v >= 0 for v in near_trk.values())), found=0,
                     prod_leg2=SENT, prod_leg0=SENT, prod_leg1=SENT, prod_nmatched=SENT, prod_anc=SENT,
                     prod_signal=SENT, reco_kpt=SENT, reco_keta=SENT, gate=SENT, mass=SENT,
                     raw_mass=SENT, leg2_dr=SENT, leg2_dpt=SENT, ok=SENT, edm=SENT, chisq=SENT, ndof=SENT, masserr=SENT,
                     k_hits=legs["K"][-1]["near_hits"], k_d0=legs["K"][-1]["near_d0"], k_d0err=legs["K"][-1]["near_d0err"],
                     k_pterr=legs["K"][-1]["near_pterr"], k_normchi2=legs["K"][-1]["near_normchi2"],
                     mup_near_dr=legs["mup"][-1]["near_dr"], mum_near_dr=legs["mum"][-1]["near_dr"],
                     mup_tight=legs["mup"][-1]["muon_tight"], mum_tight=legs["mum"][-1]["muon_tight"],
                     doca_mu0=SENT, doca_mu1=SENT, doca_max=SENT,
                     k_near_dr=legs["K"][-1]["near_dr"],
                     k_near_dpt=legs["K"][-1]["near_dpt"], k_any_pass=legs["K"][-1]["any_pass"],
                     mup_any_pass=legs["mup"][-1]["any_pass"], mum_any_pass=legs["mum"][-1]["any_pass"],
                     event=int(a["event"][i]))
            if c["have_tracks"]:
                mus = {near_trk["mup"], near_trk["mum"]}
                same_mus = np.array([{int(x), int(y)} == mus for x, y in zip(m0trk, m1trk)], dtype=bool)
                hit = np.where((np.asarray(ktrk) == near_trk["K"]) & same_mus)[0]
                if len(hit) == 1:
                    j = int(hit[0])
                    c["found"] = 1
                    g2 = a["BuJpsiK_leg2GenIdx"][i][j]
                    c["prod_leg2"] = int(g2 >= 0)
                    c["prod_leg0"] = int(a["BuJpsiK_leg0GenIdx"][i][j] >= 0)
                    c["prod_leg1"] = int(a["BuJpsiK_leg1GenIdx"][i][j] >= 0)
                    c["prod_nmatched"] = int(a["BuJpsiK_nLegsGenMatched"][i][j])
                    anc = int(a["BuJpsiK_genAncestorPdgId"][i][j])
                    c["prod_anc"] = anc
                    # the histmaker's exclusive-signal logic, in its essentials: B+ common
                    # ancestor, kaon a direct daughter of the B, muons direct daughters of a J/psi
                    c["prod_signal"] = int(abs(anc) == PDG_B and
                                           abs(int(a["BuJpsiK_leg2GenMotherPdgId"][i][j])) == PDG_B and
                                           abs(int(a["BuJpsiK_leg0GenMotherPdgId"][i][j])) == PDG_JPSI and
                                           abs(int(a["BuJpsiK_leg1GenMotherPdgId"][i][j])) == PDG_JPSI)
                    kt, m0, m1 = int(ktrk[j]), int(m0trk[j]), int(m1trk[j])
                    dd = [straight_line_doca(tpt[kt], teta[kt], tphi[kt], tdxy[kt], tdz[kt],
                                             tpt[m], teta[m], tphi[m], tdxy[m], tdz[m]) for m in (m0, m1)]
                    c["doca_mu0"], c["doca_mu1"], c["doca_max"] = dd[0], dd[1], max(dd)
                    c["reco_kpt"] = float(a["BuJpsiK_kaonPt"][i][j])
                    c["reco_keta"] = float(a["BuJpsiK_kaonEta"][i][j])
                    ok = a["BuJpsiK_jointCvhOk"][i][j]; edm = a["BuJpsiK_jointCvhEdmRef"][i][j]
                    chi = a["BuJpsiK_jointCvhChisq"][i][j]; ndf = a["BuJpsiK_jointCvhNdof"][i][j]
                    me = a["BuJpsiK_jointCvhMassErr"][i][j]; m = a["BuJpsiK_jointCvhMass"][i][j]
                    c["gate"] = int(ok == 1 and 0 <= edm < 1e-2 and 0 < chi < 1e6 and 0 < ndf < 200
                                    and 0 < me < 1 and m > -98)
                    c.update(ok=int(ok), edm=float(edm), chisq=float(chi), ndof=float(ndf), masserr=float(me))
                    c["mass"] = float(m)
                    c["raw_mass"] = float(a["BuJpsiK_kvfRawMass"][i][j])
                    c["leg2_dr"] = float(a["BuJpsiK_leg2GenDR"][i][j])
                    lp = a["BuJpsiK_leg2GenPt"][i][j]
                    c["leg2_dpt"] = float((c["reco_kpt"] - lp) / lp) if lp > 0 else SENT
                elif len(hit) > 1:
                    c["found"] = 2  # ambiguous, should not happen
            cands.append(c)
    return dict(legs=legs, cands=cands, nev=nev)


def stack(dicts):
    if not dicts:
        return {}
    keys = dicts[0].keys()
    return {k: np.array([d[k] for d in dicts]) for k in keys}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inputDir", default="/ceph/submit/data/group/cms/store/mc/inclusive_btojpsix_2016postvfp_v3_NanoAOD_v2")
    ap.add_argument("--maxFiles", type=int, default=600)
    ap.add_argument("--nproc", type=int, default=32)
    ap.add_argument("-o", "--output", default="/ceph/submit/data/user/p/pmlugato/mz/studies/genmatch_study.npz")
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.inputDir, "*.root")))[: args.maxFiles]
    print(f"{len(files)} files", flush=True)
    with Pool(args.nproc) as pool:
        results = [r for r in pool.imap_unordered(process_file, files, chunksize=4) if r is not None]
    print(f"{len(results)} files read, {sum(r['nev'] for r in results)} events", flush=True)

    out = {}
    for leg in ("K", "mup", "mum"):
        arr = stack([d for r in results for d in r["legs"][leg]])
        for k, v in arr.items():
            out[f"{leg}__{k}"] = v
    for k, v in stack([c for r in results for c in r["cands"]]).items():
        out[f"cand__{k}"] = v
    out["meta__nfiles"] = np.array(len(results))
    out["meta__nevents"] = np.array(sum(r["nev"] for r in results))
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    np.savez_compressed(args.output, **out)
    nK = len(out["K__gen_pt"])
    print(f"exclusive B+ -> J/psi(mumu) K decays: {nK}")
    print(f"  kaon: producer-style match available for {out['K__any_pass'].mean():.4f}")
    print(f"  mu+ : {out['mup__any_pass'].mean():.4f}   mu-: {out['mum__any_pass'].mean():.4f}")
    c = out
    print(f"  truth candidates found: {c['cand__found'].sum()} of {nK} "
          f"({(c['cand__have_tracks'] == 1).sum()} with all three loose tracks)")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
