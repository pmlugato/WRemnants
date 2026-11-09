#!/usr/bin/env python3
"""Stage B0 of worklist section 8: what happens to the kaon in Geant4.

Runs on the MC AlCaReco (not the nano) with FWLite, because the AlCaReco keeps
`SimTracks_g4SimHits` and `SimVertexs_g4SimHits`. For every exclusive
B+ -> J/psi(mu mu) K+ decay found in `genParticles`, the kaon's SimTrack is
located through `genpartIndex`, and its fate is read off the SimVertices whose
`parentIndex` is that track: a vertex with a decay or a hadronic-inelastic
process type is where the kaon stopped being a kaon. Everything else parented
by the kaon (delta rays, brems) leaves it alive.

The Geant4 process-type codes stored in SimVertex follow CMS's
G4ProcessTypeEnumerator (201 = decay, 121 = hadronic inelastic, 111 = hadronic
elastic, 2 = ionisation ...). They are cross-checked here against the daughter
SimTrack species at the vertex, and both are stored so the classification can
be redone offline.

Per kaon the output carries the gen kinematics, the event number (to join with
the nano-level study on (event, gen pt)), production and end vertex (r, z),
process type, daughter pdgIds, and the SimTrack momentum at production.

Usage (inside `cmsenv` of CMSSW_15_0_19_patch2):
    genmatch_simtrack_fate.py -o out.npz file1.root [file2.root ...]
"""
import argparse
import sys

import numpy as np
import ROOT  # noqa: F401  (sets up the FWLite environment)
from DataFormats.FWLite import Events, Handle

PDG_B, PDG_JPSI, PDG_K, PDG_MU, PDG_GAMMA = 521, 443, 321, 13, 22
# Process types treated as "the kaon is gone"; verified against daughters.
END_TYPES = {201: "decay", 202: "decay", 203: "decay", 121: "had_inelastic",
             131: "had_capture", 141: "had_fission", 151: "had_atrest"}
SENT = -99.0


def exclusive_decays(gp):
    """Indices (k, mu+, mu-, b) of exclusive B+ -> J/psi(mu mu) K+ in a GenParticleCollection."""
    out = []
    n = gp.size()
    for ib in range(n):
        b = gp[ib]
        if abs(b.pdgId()) != PDG_B:
            continue
        nd = b.numberOfDaughters()
        dau = [b.daughter(j) for j in range(nd)]
        if any(abs(d.pdgId()) == PDG_B for d in dau):
            continue  # not the last copy
        non_gamma = [d for d in dau if d.pdgId() != PDG_GAMMA]
        if len(non_gamma) != 2:
            continue
        pids = sorted(abs(d.pdgId()) for d in non_gamma)
        if pids != [PDG_K, PDG_JPSI]:
            continue
        k = [d for d in non_gamma if abs(d.pdgId()) == PDG_K][0]
        j = [d for d in non_gamma if d.pdgId() == PDG_JPSI][0]
        if k.status() != 1:
            continue
        # follow J/psi copies
        for _ in range(10):
            jd = [j.daughter(m) for m in range(j.numberOfDaughters())]
            same = [d for d in jd if d.pdgId() == PDG_JPSI]
            if not same:
                break
            j = same[0]
        mus = [d for d in jd if abs(d.pdgId()) == PDG_MU]
        if len(mus) != 2 or any(m.status() != 1 for m in mus):
            continue
        # collection index of the kaon: genParticles is unpruned and 1:1, so search by identity
        kidx = -1
        for ig in range(n):
            if gp[ig].pdgId() == k.pdgId() and abs(gp[ig].pt() - k.pt()) < 1e-6 and abs(gp[ig].eta() - k.eta()) < 1e-6:
                kidx = ig
                break
        out.append((kidx, k, mus, b))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+")
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--maxEvents", type=int, default=-1)
    args = ap.parse_args()

    hst, hsv, hg = Handle("std::vector<SimTrack>"), Handle("std::vector<SimVertex>"), Handle("std::vector<reco::GenParticle>")
    rows = []
    offset_votes = {0: 0, 1: 0}
    nev = 0
    for path in args.files:
        try:
            events = Events(path)
        except Exception as exc:
            sys.stderr.write(f"skip {path}: {exc}\n")
            continue
        for ev in events:
            if args.maxEvents > 0 and nev >= args.maxEvents:
                break
            nev += 1
            ev.getByLabel("genParticles", hg)
            gp = hg.product()
            decays = exclusive_decays(gp)
            if not decays:
                continue
            ev.getByLabel("g4SimHits", hst); ev.getByLabel("g4SimHits", hsv)
            st, sv = hst.product(), hsv.product()
            # index SimTracks by trackId and by genpartIndex; vertices by parent
            by_gen = {}
            by_id = {}
            for it in range(st.size()):
                t = st[it]
                by_id[t.trackId()] = it
                g = t.genpartIndex()
                if g >= 0:
                    by_gen.setdefault(g, []).append(it)
            children = {}  # vertex index -> list of SimTrack indices born there
            for it in range(st.size()):
                vi = st[it].vertIndex()
                if vi >= 0:
                    children.setdefault(vi, []).append(it)
            evid = ev.eventAuxiliary().event()
            for (kidx, k, mus, b) in decays:
                row = dict(event=int(evid), gen_pt=k.pt(), gen_eta=k.eta(), gen_phi=k.phi(), gen_q=k.charge(),
                           b_pt=b.pt(), b_eta=b.eta(), kidx=kidx, has_simtrack=0, sim_pt=SENT,
                           prod_r=SENT, prod_z=SENT, end_r=SENT, end_z=SENT, end_type=-1, end_class=0,
                           n_daughters=0, dau_pdg0=0, dau_pdg1=0, dau_pdg2=0, n_secondary_vtx=0)
                # genpartIndex convention: barcode = collection index + 1 (vote on which)
                cand = []
                for off in (0, 1):
                    for it in by_gen.get(kidx + off, []):
                        t = st[it]
                        if t.type() == k.pdgId() and abs(t.momentum().pt() - k.pt()) / max(k.pt(), 1e-6) < 0.05:
                            cand.append((off, it))
                if cand:
                    off, it = cand[0]
                    offset_votes[off] += 1
                    t = st[it]
                    row["has_simtrack"] = 1
                    row["sim_pt"] = t.momentum().pt()
                    if t.vertIndex() >= 0:
                        pv = sv[t.vertIndex()].position()
                        row["prod_r"], row["prod_z"] = pv.rho(), pv.z()
                    # vertices parented by the kaon
                    tid = t.trackId()
                    ends = []
                    nsec = 0
                    for iv in range(sv.size()):
                        v = sv[iv]
                        if v.parentIndex() != tid:
                            continue
                        nsec += 1
                        pt_ = v.processType()
                        dau = [st[c].type() for c in children.get(iv, [])]
                        # a vertex where a same-sign kaon continues is not the end
                        if pt_ in END_TYPES or (dau and not any(d == t.type() for d in dau) and pt_ not in (1, 2, 3, 4, 111)):
                            ends.append((v.position().rho(), v.position().z(), pt_, dau))
                    row["n_secondary_vtx"] = nsec
                    if ends:
                        ends.sort(key=lambda e: e[0])
                        r, z, pt_, dau = ends[0]
                        row["end_r"], row["end_z"], row["end_type"] = r, z, int(pt_)
                        row["n_daughters"] = len(dau)
                        for m, d in enumerate(dau[:3]):
                            row[f"dau_pdg{m}"] = int(d)
                        cls = END_TYPES.get(int(pt_), "other")
                        row["end_class"] = {"decay": 1, "had_inelastic": 2, "had_capture": 2, "had_fission": 2, "had_atrest": 2, "other": 3}[cls]
                        # inelastic but a same-type kaon continues: a kinked track, not a lost one
                        if row["end_class"] == 2 and any(d == t.type() for d in dau):
                            row["end_class"] = 4
                rows.append(row)
    if not rows:
        sys.stderr.write("no decays found\n")
        sys.exit(1)
    keys = rows[0].keys()
    out = {k: np.array([r[k] for r in rows]) for k in keys}
    out["meta__nevents"] = np.array(nev)
    out["meta__offset_votes"] = np.array([offset_votes[0], offset_votes[1]])
    np.savez_compressed(args.output, **out)
    print(f"{nev} events, {len(rows)} exclusive decays, simtrack found for {out['has_simtrack'].mean():.4f}; "
          f"genpartIndex offset votes (0:{offset_votes[0]}, 1:{offset_votes[1]}); "
          f"end classes: {np.unique(out['end_class'], return_counts=True)}; wrote {args.output}")


if __name__ == "__main__":
    main()
