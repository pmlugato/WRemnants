import numpy as np
L      = 16.893549842      # /fb, narf, Run2016 F+G+H Charmonium
NSIG_D = 1_363_434         # fitted, 5.05-5.45
NDAT   = 6_246_614
NSIG_M = 72_800
NMC_SEL= 282_753
NANO_M = 6_770_705
NANO_D = 87_341_047
GEN_SIG_IN_NANO = 714_190  # gen B+ -> J/psi K+ decays present in the 6.77M nano events
GEN_SIG_INACC   = 572_082
EPS_FILTER = 0.016802      # generator filter efficiency, measured on the full run
NGEN = NANO_M/EPS_FILTER

sig_bb   = 500e-6 * 1e15   # 500 ub -> fb  (sigma(pp->bbbar X), 13 TeV, 4pi)
f_u      = 0.408
BR_K     = 1.020e-3
BR_MUMU  = 5.961e-2
Nprod = sig_bb*L*2*f_u*BR_K*BR_MUMU

P=print
P(f"MARKER  DATA   lumi                          {L:.2f} /fb")
P(f"MARKER         nano events                   {NANO_D:>14,}")
P(f"MARKER         selected candidates           {NDAT:>14,}   (5.05-5.45)")
P(f"MARKER         fitted SIGNAL                 {NSIG_D:>14,}   purity {NSIG_D/NDAT:.3f}")
P(f"MARKER")
P(f"MARKER  MC     Pythia events generated       {NGEN:>14,.0f}   (nano/eps_filter, eps={EPS_FILTER})")
P(f"MARKER         nano events (filter+AlCaReco) {NANO_M:>14,}")
P(f"MARKER          of which contain gen signal  {GEN_SIG_IN_NANO:>14,}   {100*GEN_SIG_IN_NANO/NANO_M:.2f}% of nano")
P(f"MARKER          gen signal in acceptance     {GEN_SIG_INACC:>14,}")
P(f"MARKER         selected candidates           {NMC_SEL:>14,.0f}   (5.05-5.45)")
P(f"MARKER         of which gen SIGNAL           {NSIG_M:>14,}   purity {NSIG_M/NMC_SEL:.3f}")
P(f"MARKER")
P(f"MARKER  RATIO  data/MC, all selected         {NDAT/NMC_SEL:>14.1f}")
P(f"MARKER  RATIO  data/MC, SIGNAL               {NSIG_D/NSIG_M:>14.1f}   <-- the number that matters")
P(f"MARKER         MC equivalent luminosity      {L*NSIG_M/NSIG_D:>14.2f} /fb   vs {L:.2f} /fb of data")
P(f"MARKER")
P(f"MARKER  CROSS-CHECK of the data yield from first principles")
P(f"MARKER         sigma(bbbar) x L              {sig_bb*L:>14.3g} bbbar pairs")
P(f"MARKER         x 2 f_u x BR(JpsiK) x BR(mumu){Nprod:>14.3g} signal decays produced")
P(f"MARKER         observed / produced           {NSIG_D/Nprod:>14.2e}  = total efficiency")
P(f"MARKER         (acceptance x trigger x reco x selection; ~0.3% is the expected size")
P(f"MARKER          for a displaced double-muon trigger on a soft b spectrum)")
P(f"MARKER")
P(f"MARKER  EFFICIENCY per nano event, signal candidates")
P(f"MARKER         current inclusive sample      {NSIG_M/NANO_M:>14.4f}")
P(f"MARKER         exclusive B+ -> J/psi K+      {NSIG_M/GEN_SIG_IN_NANO:>14.4f}   (every event has the decay)")
P(f"MARKER         gain                          {(NSIG_M/GEN_SIG_IN_NANO)/(NSIG_M/NANO_M):>14.2f} x")
P(f"MARKER")
P(f"MARKER  PRODUCTION NEEDED, in full-simulation (nano) events")
P(f"MARKER  {'target':<22s}{'signal evts':>13s}{'exclusive':>14s}{'inclusive':>15s}{'vs sim already done':>21s}")
for lab,N in (("2x",2),("4x",4),("6x  (recommended)",6),("12x",12),("18.7x = data parity",NSIG_D/NSIG_M)):
    need=N*NSIG_M
    exc=need/(NSIG_M/GEN_SIG_IN_NANO); inc=need/(NSIG_M/NANO_M)
    P(f"MARKER  {lab:<22s}{need:>13,.0f}{exc:>14,.0f}{inc:>15,.0f}{exc/NANO_M:>20.2f}x")
