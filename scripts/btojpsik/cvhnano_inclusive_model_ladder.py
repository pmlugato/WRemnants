import numpy as np, hist, itertools, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from wremnants.postprocessing.rabbit_cvhnano_helpers import open_results
FN="/ceph/submit/data/user/p/pmlugato/mz/histograms/btojpsik_cvhnano_fullstats.hdf5"
LO,HI=5.05,5.45
with open_results(FN) as r:
    sl={"jointCvh_mass":slice(hist.loc(LO),hist.loc(HI))}
    hd=r["Charmonium_2016PostVFP_cvhNano"]["output"]["nominal"].get()[sl]
    hm=r["InclusiveBToJpsiX_2016PostVFP"]["output"]["nominal"].get()[sl]
    centers=hd.axes["jointCvh_mass"].centers
    data=hd.view()[:,7].astype(float); mcv=hm.view().astype(float)
nrm=lambda v: v/v.sum()
sig=nrm(mcv[:,0]); kx=nrm(mcv[:,1]); pix=nrm(mcv[:,2]); bd=nrm(mcv[:,3])
t=2*(centers-LO)/(HI-LO)-1.0
Bb=[]
for k in range(1,7):
    v=t**float(k); v-=v.mean()
    for u in Bb: v-=(v@u)/(u@u)*u
    Bb.append(v/np.sqrt(np.mean(v*v)))
Bb=np.array(Bb)
def sm_g(x,s):
    d=centers[:,None]-centers[None,:]; K=np.exp(-0.5*(d/max(s,1e-6))**2); K/=K.sum(axis=0,keepdims=True); return K@x
def sm_a(x,s,r):
    d=centers[:,None]-centers[None,:]
    w=np.where(d>=0,max(s,1e-6),max(s,1e-6)*np.clip(r,0.2,8.0))
    K=np.exp(-0.5*(d/w)**2); K/=K.sum(axis=0,keepdims=True); return K@x

def fit(name,order,asym,tmpls,nseed=12):
    nt=len(tmpls); npar=1+(2 if asym else 1)+nt+1+order
    def build(p):
        i=0; N=np.exp(p[i]); i+=1
        if asym: sh=sm_a(sig,abs(p[i]),abs(p[i+1])); i+=2
        else:    sh=sm_g(sig,abs(p[i])); i+=1
        m=N*nrm(sh)
        for T in tmpls: m=m+np.exp(p[i])*T; i+=1
        a=p[i]; i+=1
        return m+np.exp(a+(Bb[:order].T@p[i:i+order] if order else 0.0))
    def nll(p):
        mu=build(p)
        if not np.all(np.isfinite(mu)) or np.any(mu<=0): return 1e12
        return float(np.sum(mu-data*np.log(mu)))
    rng=np.random.RandomState(11); bestf,bx=1e30,None
    for k in range(nseed):
        p=[np.log(1.36e6)+rng.randn()*0.15]
        p.append([0.003,0.005,0.007,0.009][k%4])
        if asym: p.append([1.0,1.5,2.0,2.5,3.0][k%5])
        for _ in tmpls: p.append(np.log(1e5)+rng.randn()*0.7)
        p.append(np.log(data.mean()))
        p+= list(rng.randn(order)*0.05)
        x=np.array(p,float)
        for _ in range(2):
            rr=minimize(nll,x,method="L-BFGS-B",options=dict(maxiter=3000,maxfun=30000,ftol=1e-15))
            rr=minimize(nll,rr.x,method="Nelder-Mead",options=dict(maxiter=6000,maxfev=6000,xatol=1e-11,fatol=1e-9))
            x=rr.x
        if rr.fun<bestf: bestf,bx=rr.fun,x
    mu=build(bx); pull=(data-mu)/np.sqrt(mu); c2=float(np.sum(pull**2)); ndf=len(data)-npar
    sh=lambda lo,hi:100*float(np.sum(pull[(centers>=lo)&(centers<hi)]**2))/c2
    i=1; s=abs(bx[i]); r=abs(bx[i+1]) if asym else 1.0
    print(f"MARKER {name:<38s} {npar:>3d} {c2:>8.0f} {c2/ndf:>7.2f}  Nsig {np.exp(bx[0]):>10,.0f}"
          f"  sig {1000*s:5.2f} r {r:4.2f}   shoulder {sh(5.05,5.22):5.1f}% peak {sh(5.22,5.32):5.1f}% upper {sh(5.32,5.45):5.1f}%")
    return build,bx,c2,ndf,mu,pull

print(f"MARKER {'model':<38s} {'np':>3s} {'chi2':>8s} {'/ndf':>7s}")
fit("poly2 + Gaussian res (reference)",2,False,[])
fit("poly6, no res, no templates",6,False,[])
fit("poly2 + ASYM res",2,True,[])
fit("poly2 + Gauss res + Bd,PiX,KX",2,False,[bd,pix,kx])
b,x,c2,ndf,mu,pull=fit("poly2 + ASYM res + Bd,PiX,KX",2,True,[bd,pix,kx])
fit("poly3 + ASYM res + Bd,PiX,KX",3,True,[bd,pix,kx])
fit("poly4 + ASYM res + Bd,PiX,KX",4,True,[bd,pix,kx])

fig,(a1,a2)=plt.subplots(2,1,figsize=(11.5,5.6),sharex=True,gridspec_kw=dict(height_ratios=[2.3,1],hspace=0.06))
a1.errorbar(centers,data,yerr=np.sqrt(data),fmt="k.",ms=4,lw=1,label="Data 16.9 fb$^{-1}$")
a1.plot(centers,mu,"r-",lw=1.7,label="total fit")
i=3
_sh=sm_a(sig,abs(x[1]),abs(x[2])); _sh=_sh/_sh.sum()
a1.plot(centers,np.exp(x[0])*_sh,color="0.25",lw=1.4,label="signal, asym. resolution")
a1.plot(centers,np.exp(x[i+3]+Bb[:2].T@x[i+4:i+6]),"b:",lw=1.4,label="combinatorial (analytic)")
a1.plot(centers,np.exp(x[i])*bd,lw=1.3,ls="--",label="$B^0$ partially reconstructed")
a1.plot(centers,np.exp(x[i+2])*kx,lw=1.3,ls="--",label="$B^+\\to J/\\psi K X$")
a1.plot(centers,np.exp(x[i+1])*pix,lw=1.3,ls="--",label="$B^+\\to J/\\psi\\pi^+$")
a1.set_yscale("log"); a1.set_ylim(3e2,6e5)
a1.set_ylabel("candidates / 5 MeV")
a1.legend(frameon=False,fontsize=9,ncol=3,loc="upper left")
a2.axhline(0,color="grey",lw=0.8); a2.plot(centers,pull,"k.-",ms=4,lw=0.9)
a2.set_xlabel("m($J/\\psi K$)  [GeV]"); a2.set_ylabel("pull"); a2.set_ylim(-12,12)
a2.set_xlim(centers[0]-0.003,centers[-1]+0.003)
o=os.path.expanduser("~/public_html/mz/btojpsik_cvhnano/inclusive_datafit_physics.png")
fig.savefig(o,dpi=150,bbox_inches="tight"); print("MARKER wrote",o)
