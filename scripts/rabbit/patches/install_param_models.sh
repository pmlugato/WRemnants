#!/usr/bin/env bash
# Re-apply our param models to the rabbit submodule.
#
# rabbit is a submodule and gets reset. On 2026-08-27 a `git reset` to
# havyn/upstream-main-havyn silently removed AxisExpPolyModel and its registry
# entry mid-session, AFTER six fits had already used it -- the next fit failed
# with "Class AxisExpPolyModel not found". Run this after any submodule update,
# or before a campaign, and it is a no-op if the class is already present.
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
RABBIT=${RABBIT:-$(cd "$HERE/../../.." && pwd)/rabbit}
PM="$RABBIT/rabbit/param_models/param_model.py"
HP="$RABBIT/rabbit/param_models/helpers.py"
[[ -f $PM && -f $HP ]] || { echo "rabbit not found at $RABBIT" >&2; exit 1; }

python3 - "$PM" "$HP" "$HERE/axis_exp_poly_model.py" <<'PY'
import sys, importlib.util
pm, hp, src = sys.argv[1:4]
spec = importlib.util.spec_from_file_location("m", src)
mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)

s = open(pm).read()
if "class AxisExpPolyModel" in s:
    print("param_model.py: already present")
else:
    anchor = "class AxisBernsteinModel(ParamModel):"
    if anchor not in s:
        sys.exit("could not find AxisBernsteinModel to insert before")
    s = s.replace(anchor, mod.AXIS_EXP_POLY_MODEL.lstrip("\n") + "\n" + anchor, 1)
    open(pm, "w").write(s)
    print("param_model.py: AxisExpPolyModel inserted")

t = open(hp).read()
if '"AxisExpPolyModel"' in t:
    print("helpers.py: already registered")
else:
    key = '    "AxisExpModel": "param_model",'
    if key not in t:
        sys.exit("could not find AxisExpModel registry line")
    t = t.replace(key, key + '\n    "AxisExpPolyModel": "param_model",', 1)
    open(hp, "w").write(t)
    print("helpers.py: AxisExpPolyModel registered")
PY
# The rabbit_plot_cov fix: a missing --params entry raised TypeError (boolean
# mask on a list) instead of naming the parameter. Applied only if absent.
PATCH="$HERE/rabbit_plot_cov_missing_param.patch"
if [[ -s $PATCH ]]; then
  if grep -q "requested = np.asarray(args.params)" "$RABBIT/bin/rabbit_plot_cov.py"; then
    echo "rabbit_plot_cov.py: already patched"
  elif git -C "$RABBIT" apply --check "$PATCH" 2>/dev/null; then
    git -C "$RABBIT" apply "$PATCH" && echo "rabbit_plot_cov.py: patched"
  else
    echo "rabbit_plot_cov.py: patch does not apply cleanly -- reconcile by hand" >&2
  fi
fi

find "$RABBIT/rabbit" -name "__pycache__" -prune -exec rm -rf {} + 2>/dev/null || true
echo "done; stale __pycache__ cleared"
