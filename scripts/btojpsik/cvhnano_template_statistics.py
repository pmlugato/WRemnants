"""Is each MC template precise enough to describe the data it must describe?

Per bin: template fractional error 1/sqrt(N_MC,bin), propagated to the total
prediction through that component's data-level share, versus the data's own
Poisson precision in that bin.
"""
import numpy as np, hist
from wremnants.postprocessing.rabbit_cvhnano_helpers import open_results
FN="/ceph/submit/data/user/p/pmlugato/mz/histograms/btojpsik_cvhnano_fullstats.hdf5"
LO,HI=5.05,5.45
with open_results(FN) as r:
    sl={"jointCvh_mass":slice(hist.loc(LO),hist.loc(HI))}
    hd=r["Charmonium_2016PostVFP_cvhNano"]["output"]["nominal"].get()[sl]
    hm=r["InclusiveBToJpsiX_2016PostVFP"]["output"]["nominal"].get()[sl]
    data=hd.view()[:,7].astype(float); mcv=hm.view().astype(float)
SF=1_363_434/72_800.0   # signal-matched data/MC scale factor
names={0:"signal B+->J/psiK",1:"B+->J/psiKX",2:"B+->J/psiPiX",3:"B0",4:"Bs",6:"combinatorial(NoB)"}
NB=len(data)
print(f"MARKER {'component':<20s}{'MC evts':>10s}{'-> data':>12s}{'share':>8s}"
      f"{'MC/bin':>9s}{'tmpl err/bin':>13s}{'as frac of total':>18s}{'vs data 0.36%':>15s}")
dpb=1/np.sqrt(data.mean())
for i,nm in names.items():
    n=mcv[:,i].sum(); nd=n*SF; share=nd/data.sum()
    mpb=n/NB; terr=1/np.sqrt(max(mpb,1e-9))
    contrib=terr*share
    print(f"MARKER {nm:<20s}{n:>10,.0f}{nd:>12,.0f}{100*share:>7.1f}%{mpb:>9,.0f}"
          f"{100*terr:>12.1f}%{100*contrib:>17.2f}%{contrib/dpb:>14.1f}x")
print(f"MARKER {'':<20s}data precision per bin = {100*dpb:.2f}%   (mean {data.mean():,.0f} events/bin)")
print("MARKER")
print("MARKER per-cell (288 cells = 3 pT x 48 eta x 2 charge), MC events per template per cell:")
for i,nm in names.items():
    print(f"MARKER   {nm:<22s}{mcv[:,i].sum()/288:>10,.0f} events, spread over {NB} mass bins")
