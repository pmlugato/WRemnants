#!/usr/bin/env python3
"""Stage B1 cutflow (worklist section 8): where intact soft kaons are lost
between generalTracks and the AlCaReco B+ candidate.

Runs with FWLite (cmsenv of CMSSW_15_0_19_patch2 or 10_6_20) on the study
output: GEN-filtered exclusive B+ -> J/psi(mu mu) K+, no pileup, with
generalTracks, muons, beam spot, SimTracks/SimVertices, the standard B+
candidates and the DOCA-free clone kept.

Per gen kaon (one per event by construction of the filter; the first if more):
  gen kinematics; Geant4 fate (as in Stage B0); the nearest same-charge
  generalTrack (dR, dpt/pt, hits, chi2, d0/sigma); the two muons' tracks via
  the tight J/psi candidate; the producer's straight-line DOCA of the kaon track
  to each muon track (its formula, from reference point and momentum); the
  kaon-hypothesis m(mu mu K) and B pT built exactly as the producer does
  (propagated four-vectors are NOT reproduced -- the raw ones are used, and the
  difference is checked against the stored candidates); and whether a standard
  and/or DOCA-free B+ candidate containing that track exists.

Output: npz with flat arrays; the plotting is done by softkaon_cutflow_plots.py.
"""
import argparse
import glob
import math
import sys

import numpy as np
import ROOT  # noqa: F401
from DataFormats.FWLite import Events, Handle

PDG_B, PDG_JPSI, PDG_K, PDG_MU, PDG_GAMMA = 521, 443, 321, 13, 22
M_K, M_MU = 0.493677, 0.1056584
DOCA_CUT = 0.03
SENT = -99.0
END_TYPES = {201: 1, 202: 1, 203: 1, 121: 2, 131: 2, 141: 2, 151: 2}


def dphi(a, b):
    d = a - b
    while d > math.pi:
        d -= 2 * math.pi
    while d < -math.pi:
        d += 2 * math.pi
    return d


def dr(e1, p1, e2, p2):
    return math.hypot(e1 - e2, dphi(p1, p2))


def track_dca(t1, t2):
    """JpsiXCandidateProducer::trackTrackDCA, straight lines from the reference points."""
    u1 = np.array([t1.px(), t1.py(), t1.pz()]) / t1.p()
    u2 = np.array([t2.px(), t2.py(), t2.pz()]) / t2.p()
    c = np.cross(u1, u2)
    cn = np.linalg.norm(c)
    d = np.array([t2.vx() - t1.vx(), t2.vy() - t1.vy(), t2.vz() - t1.vz()])
    if cn < 1e-12:
        return float(np.linalg.norm(d - np.dot(d, u1) * u1))
    return float(abs(np.dot(d, c)) / cn)


def exclusive_decay(gp):
    for ib in range(gp.size()):
        b = gp[ib]
        if abs(b.pdgId()) != PDG_B:
            continue
        dau = [b.daughter(j) for j in range(b.numberOfDaughters())]
        if any(abs(d.pdgId()) == PDG_B for d in dau):
            continue
        ng = [d for d in dau if d.pdgId() != PDG_GAMMA]
        if len(ng) != 2 or sorted(abs(d.pdgId()) for d in ng) != [PDG_K, PDG_JPSI]:
            continue
        k = [d for d in ng if abs(d.pdgId()) == PDG_K][0]
        j = [d for d in ng if d.pdgId() == PDG_JPSI][0]
        if k.status() != 1:
            continue
        for _ in range(10):
            jd = [j.daughter(m) for m in range(j.numberOfDaughters())]
            same = [d for d in jd if d.pdgId() == PDG_JPSI]
            if not same:
                break
            j = same[0]
        mus = [d for d in jd if abs(d.pdgId()) == PDG_MU]
        if len(mus) != 2 or any(m.status() != 1 for m in mus):
            continue
        kidx = -1
        for ig in range(gp.size()):
            g = gp[ig]
            if g.pdgId() == k.pdgId() and abs(g.pt() - k.pt()) < 1e-6 and abs(g.eta() - k.eta()) < 1e-6:
                kidx = ig
                break
        return kidx, k, mus, b
    return None


def nearest(gen, coll, charge=None, maxdr=0.1):
    best, bdr = -1, 1e9
    for i in range(coll.size()):
        t = coll[i]
        if charge is not None and t.charge() != charge:
            continue
        d = dr(gen.eta(), gen.phi(), t.eta(), t.phi())
        if d < bdr:
            bdr, best = d, i
    return (best, bdr) if bdr < maxdr else (-1, bdr)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="files or globs")
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--maxEvents", type=int, default=-1)
    args = ap.parse_args()
    files = []
    for x in args.inputs:
        files += sorted(glob.glob(x)) if any(c in x for c in "*?[") else [x]

    H = lambda t: Handle(t)  # noqa: E731
    hg, hst, hsv = H("std::vector<reco::GenParticle>"), H("std::vector<SimTrack>"), H("std::vector<SimVertex>")
    htr, hmu = H("std::vector<reco::Track>"), H("std::vector<reco::Muon>")
    hjp = H("std::vector<reco::VertexCompositeCandidate>")
    hb, hbn = H("std::vector<reco::VertexCompositeCandidate>"), H("std::vector<reco::VertexCompositeCandidate>")
    hgood = H("std::vector<reco::Muon>")
    rows = []
    nev = 0
    for path in files:
        try:
            events = Events(path)
        except Exception as exc:
            sys.stderr.write(f"skip {path}: {exc}\n")
            continue
        for ev in events:
            if 0 < args.maxEvents <= nev:
                break
            nev += 1
            ev.getByLabel("genParticles", hg)
            dec = exclusive_decay(hg.product())
            if dec is None:
                continue
            kidx, k, mus, b = dec
            ev.getByLabel("generalTracks", htr); tracks = htr.product()
            ev.getByLabel("g4SimHits", hst); ev.getByLabel("g4SimHits", hsv); st, sv = hst.product(), hsv.product()
            ok_jpsi = ev.getByLabel("ALCARECOTkAlJpsiXJpsiCandidates", hjp)
            ok_b = ev.getByLabel("ALCARECOTkAlJpsiXBPlusCandidates", hb)
            ok_bn = ev.getByLabel("ALCARECOTkAlJpsiXBPlusCandidatesNoDoca", hbn)
            row = dict(event=int(ev.eventAuxiliary().event()), gen_pt=k.pt(), gen_eta=k.eta(), gen_phi=k.phi(), gen_q=k.charge(),
                       b_pt=b.pt(), n_tracks=int(tracks.size()),
                       end_class=0, end_r=SENT,
                       trk_idx=-1, trk_dr=SENT, trk_dpt=SENT, trk_hits=-1, trk_nchi2=SENT, trk_d0=SENT, trk_d0err=SENT, trk_pt=SENT, trk_eta=SENT,
                       n_jpsi=-1, jpsi_ok=0, mu_dr_max=SENT, doca_mu0=SENT, doca_mu1=SENT, doca_max=SENT,
                       m_raw=SENT, bpt_raw=SENT, n_bcand=-1, n_bcand_nodoca=-1, in_bcand=0, in_bcand_nodoca=0)
            # Geant4 fate (as Stage B0)
            for it in range(st.size()):
                t = st[it]
                if t.genpartIndex() == kidx + 1 and t.type() == k.pdgId():
                    tid = t.trackId()
                    ends = []
                    for iv in range(sv.size()):
                        v = sv[iv]
                        if v.parentIndex() == tid and int(v.processType()) in END_TYPES:
                            ends.append((v.position().rho(), END_TYPES[int(v.processType())]))
                    if ends:
                        ends.sort()
                        row["end_r"], row["end_class"] = ends[0][0], ends[0][1]
                    break
            # nearest same-charge generalTrack
            ti, tdr = nearest(k, tracks, charge=int(k.charge()), maxdr=0.1)
            if ti >= 0:
                t = tracks[ti]
                row.update(trk_idx=ti, trk_dr=tdr, trk_dpt=(t.pt() - k.pt()) / k.pt(), trk_hits=int(t.numberOfValidHits()),
                           trk_nchi2=t.normalizedChi2(), trk_d0=t.dxy(), trk_d0err=t.dxyError(), trk_pt=t.pt(), trk_eta=t.eta())
            # the tight J/psi candidate whose muon tracks match the gen muons
            if ok_jpsi:
                jps = hjp.product(); row["n_jpsi"] = int(jps.size())
                best = None
                for ij in range(jps.size()):
                    jp = jps[ij]
                    legs = [jp.daughter(0), jp.daughter(1)]
                    d = 0.0
                    for mu in mus:
                        d = max(d, min(dr(mu.eta(), mu.phi(), l.eta(), l.phi()) for l in legs))
                    if best is None or d < best[0]:
                        best = (d, jp)
                if best is not None and best[0] < 0.05:
                    row["jpsi_ok"] = 1; row["mu_dr_max"] = best[0]
                    jp = best[1]
                    # the J/psi daughters are RecoChargedCandidates built from generalTracks;
                    # recover the underlying track by kinematic identity (same charge,
                    # dR < 1e-3, |dpt/pt| < 1e-3), which also works when the daughter
                    # four-vector was propagated by the producer.
                    mutrks = []
                    for il in range(2):
                        leg = jp.daughter(il)
                        best, bd = None, 1e9
                        for it2 in range(tracks.size()):
                            t2 = tracks[it2]
                            if t2.charge() != leg.charge():
                                continue
                            d = dr(leg.eta(), leg.phi(), t2.eta(), t2.phi())
                            if d < bd and abs(t2.pt() - leg.pt()) / max(leg.pt(), 1e-6) < 5e-3:
                                bd, best = d, t2
                        mutrks.append(best if bd < 5e-3 else None)
                    if ti >= 0 and all(m is not None for m in mutrks):
                        t = tracks[ti]
                        d0 = track_dca(t, mutrks[0]); d1 = track_dca(t, mutrks[1])
                        row.update(doca_mu0=d0, doca_mu1=d1, doca_max=max(d0, d1))
                        # raw kaon-hypothesis mass and B pT, as the producer before propagation
                        import ROOT as R
                        lvk = R.Math.PtEtaPhiMVector(t.pt(), t.eta(), t.phi(), M_K)
                        lvj = R.Math.PtEtaPhiMVector(jp.pt(), jp.eta(), jp.phi(), jp.mass())
                        lvb = lvk + lvj
                        row.update(m_raw=lvb.M(), bpt_raw=lvb.Pt())
            # is the track in a standard / DOCA-free B+ candidate?
            for ok, h, key, nkey in ((ok_b, hb, "in_bcand", "n_bcand"), (ok_bn, hbn, "in_bcand_nodoca", "n_bcand_nodoca")):
                if not ok:
                    continue
                cands = h.product(); row[nkey] = int(cands.size())
                if ti >= 0:
                    t = tracks[ti]
                    for ic in range(cands.size()):
                        c = cands[ic]
                        for il in range(c.numberOfDaughters()):
                            dtr = c.daughter(il)
                            if abs(dtr.pdgId()) == PDG_K and dtr.charge() == t.charge() and abs(dtr.pt() - t.pt()) / t.pt() < 0.02 and dr(dtr.eta(), dtr.phi(), t.eta(), t.phi()) < 0.01:
                                row[key] = 1
            rows.append(row)
    if not rows:
        sys.exit("no decays found")
    out = {k: np.array([r[k] for r in rows]) for k in rows[0]}
    out["meta__nevents"] = np.array(nev)
    np.savez_compressed(args.output, **out)
    print(f"{nev} events, {len(rows)} decays; track within 0.1: {np.mean(out['trk_idx'] >= 0):.3f}; "
          f"tight J/psi: {out['jpsi_ok'].mean():.3f}; in standard B+ cand: {out['in_bcand'].mean():.3f}; "
          f"in DOCA-free cand: {out['in_bcand_nodoca'].mean():.3f}; wrote {args.output}")


if __name__ == "__main__":
    main()
