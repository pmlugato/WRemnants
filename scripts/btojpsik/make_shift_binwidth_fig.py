import os, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
OUT=os.path.expanduser("~/public_html/slides/assets_btojpsik_overview")
INK,ACC,GREY="#1F1F1F","#C0392B","#7F7F7F"
# left: what one free shift parameter buys, 2018 vs ours
# right: pure averaging -- rebinning data AND model together, no refit
f,(a1,a2)=plt.subplots(1,2,figsize=(11.2,4.3))
lab=["Gaussian smear\n+ exp·poly2","+ one free\nmass shift"]
x=np.arange(2)
a1.plot(x,[66.10,1.44],"o-",color=INK,ms=9,lw=2.2,label="2018 (BMM nano)")
a1.plot(x,[53.59,43.94],"s-",color=ACC,ms=9,lw=2.2,label="2016postVFP (CVH)")
a1.set_yscale("log"); a1.set_xticks(x); a1.set_xticklabels(lab,fontsize=9)
a1.set_ylabel(r"$\chi^2$/ndf"); a1.set_xlim(-0.35,1.35)
a1.axhline(1.0,color=GREY,ls=":",lw=1.0)
a1.annotate(r"$\Delta\chi^2$ = 6145",(0.5,9),color=INK,fontsize=10,ha="center")
a1.annotate(r"$\Delta\chi^2$ = 778",(0.52,33),color=ACC,fontsize=10,ha="center")
a1.legend(frameon=False,fontsize=9.5,loc="lower left")
a1.set_title("2018's residual was a $-1.44$ MeV scale offset",fontsize=10)
a1.grid(alpha=0.15)

w18=[2,10,20]; w16=[5,10,20]
a2.plot(w18,[1.44,2.66,3.64],"o-",color=INK,ms=7,lw=2.0,label="2018, free shift")
a2.plot(w18,[66.10,396.0,1049.0],"o--",color=INK,ms=5,lw=1.3,alpha=0.55,label="2018, no shift")
a2.plot(w16,[43.94,93.02,207.07],"s-",color=ACC,ms=7,lw=2.0,label="2016, free shift")
a2.axvline(20,color=GREY,ls=":",lw=1.2)
a2.text(19.4,1.6,"2018 tensors\nused 20 MeV",fontsize=9,color=GREY,ha="right")
a2.set_yscale("log"); a2.set_xlabel("mass bin width  [MeV]")
a2.set_ylabel(r"$\chi^2$/ndf"); a2.axhline(1.0,color=GREY,ls=":",lw=1.0)
a2.legend(frameon=False,fontsize=9)
a2.set_title("Coarsening makes it worse, not invisible",fontsize=10)
a2.grid(alpha=0.15)
p=os.path.join(OUT,"shift_and_binwidth.png")
f.savefig(p,dpi=150,bbox_inches="tight"); print("MARKER wrote",p)
