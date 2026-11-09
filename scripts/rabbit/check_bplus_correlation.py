#!/usr/bin/env python3
"""Does the merge actually correlate B+ with the Z and J/psi channels?

merge_inputs.py correlates systematics by NAME and by nothing else, so a merge
that correlates nothing looks exactly like a merge that correlates everything:
it writes a tensor, the fit converges, and the answer is quietly wrong. This
gates the fit on the one thing that cannot be seen afterwards.

The check reads `hsysts` straight out of each input tensor and intersects the
name sets. That is exactly what merge_inputs does, so it cannot disagree with
it, and it needs nothing from rabbit's internals.

    python scripts/rabbit/check_bplus_correlation.py --nEta 24 \
        --bplus B.hdf5 --others Z.hdf5 jpsi_prompt.hdf5 jpsi_displaced.hdf5 \
        [--merged merged.hdf5] [--rename scripts/rabbit/lowmass_z_rename.py]
"""

import argparse
import importlib.util
import re
import sys

import h5py

AEM = re.compile(r"^([AeM])(\d+)$")


def names(path):
    with h5py.File(path, "r") as f:
        return [
            s.decode() if isinstance(s, bytes) else str(s) for s in f["hsysts"][...]
        ]


def load_rename(path):
    spec = importlib.util.spec_from_file_location("_rn", path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return getattr(m, "rename_systs", {})


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--bplus", required=True)
    p.add_argument("--others", nargs="+", required=True)
    p.add_argument("--merged")
    p.add_argument("--rename")
    p.add_argument("--nEta", type=int, default=24)
    a = p.parse_args()

    ren = load_rename(a.rename) if a.rename else {}
    bn = names(a.bplus)
    print(f"B+ tensor: {len(bn)} systematics")

    other = {}
    for f in a.others:
        n = [ren.get(x, x) for x in names(f)]
        other[f.split("/")[-1]] = set(n)
        print(f"  {f.split('/')[-1]:<48s} {len(n)} systematics (after rename)")

    union_other = set().union(*other.values()) if other else set()
    shared = sorted(set(bn) & union_other)
    aem_b = sorted(
        [s for s in bn if AEM.match(s)],
        key=lambda s: (AEM.match(s).group(1), int(AEM.match(s).group(2))),
    )

    print(
        f"\nB+ names matching the A/e/M basis pattern: {len(aem_b)} "
        f"(expected {3 * a.nEta})"
    )
    fail = False
    for param in ("A", "e", "M"):
        idx = sorted(
            int(AEM.match(s).group(2)) for s in aem_b if AEM.match(s).group(1) == param
        )
        missing = [i for i in range(a.nEta) if i not in idx]
        extra = [i for i in idx if i >= a.nEta]
        ok = not missing and not extra
        fail |= not ok
        print(
            f"  {'OK  ' if ok else 'FAIL'} {param}: {len(idx)} present"
            + (
                f", missing {missing[:6]}{'...' if len(missing) > 6 else ''}"
                if missing
                else ""
            )
            + (f", out of range {extra[:6]}" if extra else "")
        )

    print(f"\nnames shared between B+ and the other channels: {len(shared)}")
    per = {k: len(set(bn) & v) for k, v in other.items()}
    for k, v in per.items():
        print(f"  with {k:<48s} {v}")
    if not shared:
        print(
            "\nFAIL: B+ shares NO systematic with any other channel. Its tensor "
            "was written with the standalone naming (scaleA0..); rebuild it "
            "with --scaleParamNaming joint."
        )
        return 1
    aem_shared = [s for s in shared if AEM.match(s)]
    print(f"  ...of which are A/e/M basis parameters: {len(aem_shared)}")

    if a.merged:
        mn = names(a.merged)
        expected = len(set(bn) | union_other)
        print(
            f"\nmerged tensor: {len(mn)} systematics, "
            f"union of the inputs is {expected}"
        )
        if len(mn) != expected:
            print(
                "  WARN: the merged count is not the union of the inputs. "
                "Either something else was renamed, or a name collided that "
                "was not meant to."
            )

    ok = (len(aem_shared) == 3 * a.nEta) and not fail
    print(
        "\n"
        + ("PASS" if ok else "FAIL")
        + f": {len(aem_shared)} of {3 * a.nEta} basis parameters correlate B+ "
        "with the Z / J/psi channels"
    )
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
