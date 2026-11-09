"""Eigen-analysis of the fit Hessian at the prefit point (or at a saved
parameter vector) for a rabbit tensor, to name the direction of negative
curvature that makes trust-exact/BFGS die at iteration 0 and trust-krylov stall
on the observation.

Same command line as rabbit_fit.py (tensor, --paramModel ..., -t 0), plus
    --topk N     how many eigen-directions to print (default 5)

Example:
    python scripts/btojpsik/cvhnano_prefit_hessian.py <tensor.hdf5> \\
        --paramModel AxisNormModel btojpsik_cvhnano BuJpsiK fit_bachelor_eta ... -t 0
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "rabbit", "bin"))
import rabbit_fit  # noqa: E402
from rabbit import fitter as rfitter  # noqa: E402
from rabbit import inputdata  # noqa: E402
from rabbit.param_models import helpers as ph  # noqa: E402


def main():
    argv = sys.argv[1:]
    topk = 5
    if "--topk" in argv:
        i = argv.index("--topk")
        topk = int(argv[i + 1])
        del argv[i:i + 2]
    parser = rabbit_fit.make_parser()
    args = parser.parse_args(argv)

    indata = inputdata.FitInputData(args.filename, args.pseudoData)
    param_model = ph.load_models(args.paramModel or [["Mu"]], indata, **vars(args))
    fit = rfitter.Fitter(indata, param_model, args, do_blinding=False)
    fit.defaultassign()
    fit.set_nobs(indata.data_obs)

    val, grad, hess = fit.loss_val_grad_hess()
    grad = np.asarray(grad).reshape(-1)
    hess = np.asarray(hess)
    names = [n.decode() if isinstance(n, bytes) else str(n) for n in fit.parms]
    print(f"loss at prefit: {float(val):.6g}   nparams: {len(names)}")
    diag = np.diag(hess)
    neg = np.where(diag < 0)[0]
    print(f"negative diagonal entries: {len(neg)}")
    for i in neg[:20]:
        print(f"   {names[i]:<45s} H_ii = {diag[i]:+.4g}   grad = {grad[i]:+.4g}")
    w, v = np.linalg.eigh(hess)
    print(f"eigenvalues: min {w[0]:.4g}  max {w[-1]:.4g}  n negative {int((w < 0).sum())}")
    for k in range(min(topk, len(w))):
        if w[k] >= 0:
            break
        comp = v[:, k]
        order = np.argsort(-np.abs(comp))[:6]
        print(f"  lambda = {w[k]:+.4g}:")
        for i in order:
            print(f"      {comp[i]:+.3f}  {names[i]}  (grad {grad[i]:+.3g}, H_ii {diag[i]:+.3g})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
