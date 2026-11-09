"""Derive true-pileup weights for the CVH NanoAOD simulation.

`Pileup_nTrueInt` is the generated Poisson mean of the pileup distribution, and
it is what a pileup profile is defined on. The weight is therefore

    w(n) = P_data(n) / P_MC(n)

with the data profile taken from the standard WRemnants 2016 postVFP file --
`pileupProfileData_2016Legacy_RunpostVFP_04June2021.root`, the same one
`wremnants.production.pileup.make_pileup_helper` uses -- and the simulation's
profile measured directly from the NanoAOD files.

Two choices here differ from `make_pileup_helper`, and both are forced:

* **The MC profile is measured, not read from a file.** The tabulated MC
  profiles in `wremnants-data` are for the central W/Z productions. This is a
  private inclusive b -> J/psi X sample with its own pileup mixing, so its
  profile has to come from its own files.
* **It is measured over all events in the files, before any selection.** A
  pileup weight corrects the generated profile. Deriving it from selected events
  would fold the selection efficiency's pileup dependence into the weight and
  then apply it again, double-counting.

Bins where the simulation has fewer than `--minEntries` events keep weight 1.0:
a ratio measured on a handful of events is noise, and propagating it would
inject that noise into every distribution rather than only into the pileup.
Weights are capped at `--cropHighWeight`, as the central helper does.

    python scripts/plotting/btojpsik_cvhnano_pileup.py -o pu.json
"""

import argparse
import glob
import json
import os

import numpy as np
import uproot

DEFAULT_MC = (
    "/ceph/submit/data/group/cms/store/mc/"
    "inclusive_btojpsix_2016postvfp_v3_NanoAOD_v2/*.root"
)
DEFAULT_PROFILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "../../wremnants-data/data/pileupProfiles",
    "pileupProfileData_2016Legacy_RunpostVFP_04June2021.root",
)

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-o", "--output", default="pileup_weights.json")
parser.add_argument("--mcFiles", default=DEFAULT_MC, help="glob of simulation NanoAODs")
parser.add_argument("--maxFiles", type=int, default=400)
parser.add_argument("--dataProfile", default=DEFAULT_PROFILE)
parser.add_argument("--profileHist", default="pileup")
parser.add_argument("--nbins", type=int, default=100)
parser.add_argument("--minEntries", type=float, default=200.0)
parser.add_argument("--cropHighWeight", type=float, default=5.0)
args = parser.parse_args()

# --- data profile ---------------------------------------------------------
with uproot.open(args.dataProfile) as f:
    prof = f[args.profileHist].to_numpy()
dv = np.asarray(prof[0], dtype=float)
dv = np.pad(dv, (0, max(0, args.nbins - len(dv))))[: args.nbins]

# --- simulated profile, from the files ------------------------------------
files = sorted(glob.glob(args.mcFiles))[: args.maxFiles]
if not files:
    raise SystemExit(f"no simulation files matched {args.mcFiles}")
mv = np.zeros(args.nbins, dtype=float)
nread = 0
for path in files:
    try:
        with uproot.open(path) as f:
            n = f["Events"]["Pileup_nTrueInt"].array(library="np")
    except Exception as exc:  # truncated file from a killed job
        print(f"  unreadable, skipping: {path} ({exc})")
        continue
    mv += np.histogram(n, bins=args.nbins, range=(0.0, float(args.nbins)))[0]
    nread += 1

print(f"read {nread} of {len(files)} simulation files, {mv.sum():.0f} events")

dn = dv / dv.sum()
mn = mv / mv.sum()

weights = {}
for i in range(args.nbins):
    if mv[i] < args.minEntries or dv[i] <= 0:
        continue
    weights[i] = float(min(dn[i] / mn[i], args.cropHighWeight))

with open(args.output, "w") as fh:
    json.dump(weights, fh, indent=1, sort_keys=True)

covered = sum(mv[i] for i in weights)
mean_d = float((np.arange(args.nbins) * dn).sum())
mean_m = float((np.arange(args.nbins) * mn).sum())
print(f"wrote {len(weights)} weights to {args.output}")
print(f"  range {min(weights.values()):.3f} - {max(weights.values()):.3f}")
print(f"  mean nTrueInt: data {mean_d:.2f}, simulation {mean_m:.2f}")
print(f"  covers {covered/mv.sum():.4f} of simulated events; the rest keep weight 1")
