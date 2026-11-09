"""Induced mean B+ mass shift per unit of each A/e/M nuisance, read from B+ tensors.

    python template_induced_shift.py [--tensors LABEL=tensor.hdf5 ...] [--params A e M]

Relative paths are taken under the ceph tensor directory. Defaults reproduce the
2026-09-21 all-legs vs kaon-only comparison.
"""
import argparse, os
import hdf5plugin, h5py, numpy as np, re, sys
T="/ceph/submit/data/user/p/pmlugato/mz/tensors/"
ap=argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--tensors", nargs="+", default=["all legs=btojpsik_cvhnano_tensor_joint24_res5_260914.hdf5",
                                                 "kaon only=btojpsik_cvhnano_tensor_joint24_kaononly_mz.hdf5"])
ap.add_argument("--params", nargs="+", default=["A","e"])
ap.add_argument("--theta", type=float, default=1.0,
                help="evaluate rabbit's interpolation exp(+-theta*logk) at this theta and "
                "divide by it: 1 reproduces the up/down templates, small theta gives the slope the fit uses")
cli=ap.parse_args()
for spec in cli.tensors:
    tag,f=spec.split("=",1)
    h=h5py.File(os.path.join(T,f)); systs=[s.decode() for s in h["hsysts"][...]]
    procs=[s.decode() for s in h["hprocs"][...]]; isig=procs.index("BuJpsiK")
    norm=np.asarray(h["hnorm"][...]).reshape(-1,len(procs))      # (nbins, nproc)
    logk=np.asarray(h["hlogk"][...]).reshape(norm.shape[0],len(procs),len(systs))
    # bins are (pt, eta, charge, mass) C-order with 3 x 24 x 2 x 20
    npt,neta,nq,nm=3,24,2,20
    s=norm[:,isig].reshape(npt,neta,nq,nm)
    m=np.linspace(5.05,5.45,nm+1); mc=0.5*(m[1:]+m[:-1])
    print(f"\n=== {tag}")
    for p in cli.params:
        idx=[systs.index(f"{p}{i}") for i in range(neta) if f"{p}{i}" in systs]
        sh=[]
        for i,j in enumerate(idx):
            lk=logk[:,isig,j].reshape(npt,neta,nq,nm)
            up=s*np.exp(cli.theta*lk); dn=s*np.exp(-cli.theta*lk)
            # mean mass in this eta group only (the slab the nuisance acts on)
            for ipt in range(npt):
                w_up=up[ipt,i].sum(axis=0); w_dn=dn[ipt,i].sum(axis=0); w0=s[ipt,i].sum(axis=0)
                if w0.sum()<=0: continue
                mu_up=(w_up*mc).sum()/w_up.sum(); mu_dn=(w_dn*mc).sum()/w_dn.sum()
                sh.append((ipt,1e3*0.5*(mu_up-mu_dn)/cli.theta))
        for ipt in range(npt):
            v=[x for j,x in sh if j==ipt]
            if v: print(f"  {p}: pT bin {ipt}: induced mean shift per unit = {np.mean(v):+.4f} MeV (rms over eta {np.std(v):.4f})")
