#!/usr/bin/env python3
"""Task 1.0d gate: merged-channel-list output vs single-channel output.

The merged path (`channels=BPlus`) concatenates exactly one collection, so it
must reproduce the single-channel path (`srcCandidates=...BPlusResonances`)
column for column. Anything else is a bug in the merger or in the wiring.

Two differences are EXPECTED and are not failures:
  * the table name (`JpsiX` vs `BuJpsiK`), because the merged form cannot
    know it is carrying only one channel;
  * one extra column, `channelId`, which must be all zeros for a one-channel
    list.

Both files must come from the SAME BUILD. Comparing against the August
production nano would conflate this change with six weeks of upstream
Geant4ePropagator work (openspec design.md D9a-ter).

  compare_merged_vs_single.py --merged mrg_bponly.root --single tbl_BPlus.root
"""

import argparse
import sys

import numpy as np
import uproot

ap = argparse.ArgumentParser()
ap.add_argument("--merged", required=True)
ap.add_argument("--single", required=True)
ap.add_argument("--merged-table", default="JpsiX")
ap.add_argument("--single-table", default="BuJpsiK")
ap.add_argument(
    "--rtol", type=float, default=0.0, help="0 (default) demands exact equality"
)
ap.add_argument(
    "--no-channel-id",
    action="store_true",
    help="skip the channelId check. Use when comparing two "
    "single-channel files (e.g. across a rebuild), where "
    "neither side has the column.",
)
args = ap.parse_args()


def flat(tree, branch):
    a = tree[branch].array(library="np")
    if a.dtype == object:
        return np.concatenate([np.asarray(x) for x in a]) if len(a) else np.array([])
    return a


tm = uproot.open(args.merged)["Events"]
ts = uproot.open(args.single)["Events"]

cols_m = {
    k[len(args.merged_table) + 1 :]
    for k in tm.keys()
    if k.startswith(args.merged_table + "_")
}
cols_s = {
    k[len(args.single_table) + 1 :]
    for k in ts.keys()
    if k.startswith(args.single_table + "_")
}

only_m = sorted(cols_m - cols_s)
only_s = sorted(cols_s - cols_m)
shared = sorted(cols_m & cols_s)

print(f"merged table  '{args.merged_table}' : {len(cols_m)} columns")
print(f"single table  '{args.single_table}' : {len(cols_s)} columns")
print(f"shared                              : {len(shared)}")
if only_m:
    print(f"  only in merged : {only_m}")
if only_s:
    print(f"  only in single : {only_s}")

fails, checked = [], 0
for c in shared:
    a = flat(tm, f"{args.merged_table}_{c}")
    b = flat(ts, f"{args.single_table}_{c}")
    if a.shape != b.shape:
        fails.append((c, f"shape {a.shape} vs {b.shape}"))
        continue
    if a.dtype.kind == "f":
        same = (
            np.array_equal(a, b)
            if args.rtol == 0
            else np.allclose(a, b, rtol=args.rtol, equal_nan=True)
        )
        if not same:
            # NaN-aware exact comparison
            both_nan = np.isnan(a) & np.isnan(b)
            if not np.array_equal(a[~both_nan], b[~both_nan]):
                d = np.abs(a[~both_nan] - b[~both_nan])
                fails.append((c, f"max |diff| = {d.max():.6g}"))
                continue
    elif not np.array_equal(a, b):
        fails.append((c, f"{int((a != b).sum())} of {a.size} differ"))
        continue
    checked += 1

print(f"\n{checked} of {len(shared)} shared columns identical")

# channelId must exist and be all zeros for a one-channel list
cid_ok = None
if args.no_channel_id:
    print("channelId check skipped (--no-channel-id)")
elif "channelId" in cols_m:
    cid = flat(tm, f"{args.merged_table}_channelId")
    # bool(...) matters: np.all returns np.bool_, and `is True` is
    # False for np.True_, which silently turned a pass into a FAIL.
    cid_ok = bool(cid.size > 0 and np.all(cid == 0))
    print(
        f"channelId present, all zero: {cid_ok}  "
        f"(n={cid.size}, unique={np.unique(cid).tolist()})"
    )
else:
    print("channelId MISSING from the merged table")

if fails:
    print(f"\nFAIL: {len(fails)} columns differ")
    for c, why in fails[:20]:
        print(f"  {c:34s} {why}")
    sys.exit(1)
if only_s:
    print(f"\nFAIL: single-channel columns absent from merged: {only_s}")
    sys.exit(1)
if not args.no_channel_id and not cid_ok:
    print("\nFAIL: channelId missing or non-zero on a one-channel list")
    sys.exit(1)

if args.no_channel_id:
    print(f"\nPASSED: the two files agree on all {checked} shared columns")
else:
    print(
        "\nGATE PASSED: merged(BPlus) reproduces single-channel B+ exactly, "
        "plus channelId == 0"
    )
