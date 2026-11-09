import os, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from wremnants.postprocessing.rabbit_cvhnano_helpers import open_results
H="/ceph/submit/data/user/p/pmlugato/mz/histograms/"; GEN=0
def vals(h):
    a=np.asarray(h.view()); return a["value"].astype(float) if a.dtype.names else a.astype(float)
def load(fn):
    with open_results(H+fn) as r:
        o=r["InclusiveBToJpsiX_2016PostVFP"]["output"]
        hv=o["nominal_muonScaleSyst_responseWeights"].get()
        names=[a.name for a in hv.axes]
        v=np.take(vals(hv),GEN,axis=names.index("genCategory"))
        pe=np.asarray(hv.axes["fit_bachelor_pt"].edges)
        nb=v.shape[3]; me=np.linspace(5.0,5.5,nb+1)
    return pe,0.5*(me[:-1]+me[1:]),v
pe,mc,vf=load("btojpsik_cvhnano_fullstats.hdf5")
_,_,vn=load("btojpsik_cvhnano_full_nofloor.hdf5")
pc=0.5*(pe[:-1]+pe[1:])
def shift(w):
    s=w.sum(); return float((w*mc).sum()/s) if s>0 else np.nan
bands=[(1.0,1.25),(1.25,1.5),(1.5,1.75),(1.75,2.0),(2.0,2.5),(2.5,3.0),(3.0,4.0),(4.0,5.5),(5.5,8.0)]
x=[0.5*(a+b) for a,b in bands]; xe=[0.5*(b-a) for a,b in bands]
fig,ax=plt.subplots(figsize=(9.6,4.6))
for u,lab,col in ((1,"$e$, $\\eta$ group 0","C0"),(4,"$e$, $\\eta$ group 1","C1"),(0,"$A$, $\\eta$ group 0","C2")):
    for V,ls,tag in ((vf,"-","clamped at 2 GeV"),(vn,"--","network extrapolated")):
        ys=[]
        for a,b in bands:
            m=(pc>=a)&(pc<b)
            ys.append(0.5*(shift(V[m][...,u,1].sum(axis=(0,1,2)))-shift(V[m][...,u,0].sum(axis=(0,1,2))))*1e3)
        ax.errorbar(x,ys,xerr=xe,fmt="o" if ls=="-" else "s",ls=ls,ms=4,lw=1.3,color=col,
                    mfc=col if ls=="-" else "none",
                    label=f"{lab} — {tag}")
ax.axvline(2.0,color="k",lw=1.0,ls=":")
ax.text(2.05,ax.get_ylim()[1]*0.92,"training floor",fontsize=9)
ax.axhline(0,color="grey",lw=0.7)
ax.set_xlabel("bachelor (kaon) reco $p_T$  [GeV]")
ax.set_ylabel("induced peak shift  [MeV]")
ax.set_xscale("log")
import matplotlib.ticker as mt
ax.xaxis.set_major_locator(mt.FixedLocator([1,1.5,2,3,4,5,6,8]))
ax.xaxis.set_minor_locator(mt.NullLocator())
ax.xaxis.set_major_formatter(mt.FixedFormatter(["1","1.5","2","3","4","5","6","8"]))
ax.legend(frameon=False,fontsize=8.5,ncol=2)
o=os.path.expanduser("~/public_html/mz/btojpsik_cvhnano/pt_floor_variations.png")
fig.savefig(o,dpi=150,bbox_inches="tight"); print("MARKER wrote",o)
