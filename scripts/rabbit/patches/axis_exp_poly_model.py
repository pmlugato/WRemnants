"""AxisExpPolyModel -- kept OUTSIDE the rabbit submodule on purpose.

The rabbit checkout is a submodule and gets reset (2026-08-27: `git reset` to
`havyn/upstream-main-havyn` wiped this class and its registry entry mid-session,
after six fits had already used it). So the source of truth lives here and
`install_param_models.sh` re-applies it. Run that after any submodule update.
"""

AXIS_EXP_POLY_MODEL = '''

class AxisExpPolyModel(ParamModel):
    """
    Per-(process, cell) exponential-polynomial background param model, general
    order.

        rnorm = exp(lnAmpl_c + sum_{k=1..N} c_k,g_k * P_k(t))

    with t the shape-axis bin centre mapped affinely onto [-1, 1] and P_k an
    orthogonal polynomial basis built ON THE ACTUAL BIN GRID.  Order N = 1 is
    AxisExpModel (a slope), N = 2 adds curvature, and each higher order adds a
    genuinely new direction.

    **Why the basis is orthogonalised on the grid rather than nominally.**  Any
    polynomial basis of the same degree spans the same functions, so the fitted
    optimum cannot tell them apart -- but the Hessian can, and a fit that stalls
    on a near-flat direction is the failure this basis exists to avoid.  Nominal
    Legendre polynomials are orthogonal as integrals over [-1, 1], which is not
    orthogonality over 80 bin centres: anchoring t at the first and last bin
    CENTRE puts both at exactly +/-1 and leaves <P_2, 1> = +1.3e-2 with mutual
    overlaps to 8e-2.  The overlap with the CONSTANT is the one that matters,
    because lnAmpl is a free constant, one per cell.  Mapping with the axis EDGES
    and running modified Gram-Schmidt against the constant and each lower order
    gives overlaps of 5e-17 to 5e-16 for orders 1-4, for any binning.

    Grouping: lnAmpl over cell_axes, and each polynomial order over its own
    subset of the cell axes.  Fewer groupings than the order means the remaining
    orders reuse the last one.  A coarser grouping for shape than for amplitude
    is the normal case: every cell needs its own rate, but a shape that is a
    property of the background's kinematics is shared, and sharing it is what
    keeps the parameter count from tracking the cell count.

    All parameters are unconstrained reals -- exp() of anything real is positive,
    so there is no boundary, no softplus, and no coefficient-positivity
    degeneracy.  All zeros is a flat unit background, and every order nests the
    one below it at an interior point, so a chi-square difference between
    adjacent orders is a likelihood ratio for one parameter.

    Usage::

        --paramModel AxisExpPolyModel <channel> <proc_spec> <shape_axis> \\\\
            <cell_axes> <order> [<axes_for_order_1>[ <axes_for_order_2> ...]]
    """

    @classmethod
    def parse_args(cls, indata, *args, **kwargs):
        if len(args) < 5:
            raise ValueError(
                f"AxisExpPolyModel requires at least 5 positional arguments "
                f"(channel, proc_spec, shape_axis, cell_axes, order) but got "
                f"{len(args)}: {args}"
            )
        channel, proc_spec, shape_axis, cell_axes_csv, order = args[:5]
        return cls(
            indata,
            channel,
            proc_spec,
            shape_axis,
            cell_axes_csv,
            int(order),
            order_axes_csv=list(args[5:]) or None,
            **kwargs,
        )

    def __init__(
        self,
        indata,
        channel,
        proc_spec,
        shape_axis,
        cell_axes_csv,
        order,
        order_axes_csv=None,
        expectSignal=None,
        allowNegativeParam=False,
        **kwargs,
    ):
        self.indata = indata
        if order < 1:
            raise ValueError(f"AxisExpPolyModel order must be >= 1, got {order}")
        self.order = int(order)

        if channel not in indata.channel_info:
            raise ValueError(
                f"Channel '{channel}' not found in tensor. "
                f"Available: {list(indata.channel_info.keys())}"
            )
        self.channel = channel
        channel_axes = indata.channel_info[channel]["axes"]
        axis_by_name = {a.name: a for a in channel_axes}

        if shape_axis not in axis_by_name:
            raise ValueError(
                f"Shape axis '{shape_axis}' not found in channel '{channel}'. "
                f"Available: {list(axis_by_name.keys())}"
            )

        cell_names = [n.strip() for n in cell_axes_csv.split(",")]
        for name in cell_names:
            if name not in axis_by_name:
                raise ValueError(
                    f"Cell axis '{name}' not found in channel '{channel}'. "
                    f"Available: {list(axis_by_name.keys())}"
                )
            if name == shape_axis:
                raise ValueError(
                    f"Axis '{name}' appears in both shape_axis and cell_axes."
                )
        self.cell_axis_names = set(cell_names)
        self.cell_axes = [axis_by_name[n] for n in cell_names]
        self.shape_axis = shape_axis

        if order_axes_csv is None:
            order_axes_csv = [cell_axes_csv]
        groupings = []
        for k in range(self.order):
            csv = order_axes_csv[min(k, len(order_axes_csv) - 1)]
            names = [n.strip() for n in csv.split(",")]
            bad = [n for n in names if n not in self.cell_axis_names]
            if bad:
                raise ValueError(
                    f"Order-{k+1} axes {bad} are not in cell_axes "
                    f"'{cell_axes_csv}'. They must be a subset of cell_axes."
                )
            groupings.append(names)
        self.order_axis_names = [set(g) for g in groupings]
        self.order_axes = [[axis_by_name[n] for n in g] for g in groupings]
        self.n_order_groups = [
            int(np.prod([a.size for a in axes])) for axes in self.order_axes
        ]

        if proc_spec == "all":
            target_encoded = list(indata.procs)
        else:
            target_encoded = []
            for name in [p.strip() for p in proc_spec.split(",")]:
                encoded = name.encode() if isinstance(name, str) else name
                if encoded not in indata.procs:
                    raise ValueError(
                        f"Process '{name}' not found in tensor. "
                        f"Available: {[p.decode() if isinstance(p, bytes) else p for p in indata.procs]}"
                    )
                target_encoded.append(encoded)
        self.proc_idxs = [
            int(np.where(indata.procs == p)[0][0]) for p in target_encoded
        ]

        cell_shape = [a.size for a in self.cell_axes]
        self.n_cell = int(np.prod(cell_shape))
        self.stride = self.n_cell + sum(self.n_order_groups)
        self.npoi = len(self.proc_idxs) * self.stride

        names = []
        for proc_encoded in target_encoded:
            proc_name = (
                proc_encoded.decode()
                if isinstance(proc_encoded, bytes)
                else str(proc_encoded)
            )
            for idxs in itertools.product(*[range(s) for s in cell_shape]):
                label = "_".join(f"{a.name}{i}" for a, i in zip(self.cell_axes, idxs))
                names.append(f"lnAmpl_{proc_name}_{label}".encode())
            for k in range(self.order):
                axes = self.order_axes[k]
                shp = [a.size for a in axes]
                for idxs in itertools.product(*[range(s) for s in shp]):
                    label = "_".join(f"{a.name}{i}" for a, i in zip(axes, idxs))
                    names.append(f"poly{k+1}_{proc_name}_{label}".encode())
        self.params = np.array(names)

        centers = np.asarray(axis_by_name[shape_axis].centers, dtype=np.float64)
        edges = np.asarray(axis_by_name[shape_axis].edges, dtype=np.float64)
        lo, hi = float(edges[0]), float(edges[-1])
        span = max(hi - lo, 1e-6)
        t = 2.0 * (centers - lo) / span - 1.0

        basis = []
        for k in range(1, self.order + 1):
            v = t.astype(np.float64) ** k
            v = v - v.mean()
            for u in basis:
                v = v - (v @ u) / (u @ u) * u
            nrm = float(np.sqrt(np.mean(v * v)))
            if nrm < 1e-12:
                raise ValueError(
                    f"Order {k} is degenerate on this shape axis "
                    f"({len(centers)} bins): the basis collapsed. Lower the order."
                )
            basis.append(v / nrm)
        self.p_k = [tf.constant(v, dtype=indata.dtype) for v in basis]

        self.basis_max_overlap = max(
            [abs(float(v.mean())) for v in basis]
            + [
                abs(float(basis[i] @ basis[j]) / len(t))
                for i in range(len(basis))
                for j in range(i + 1, len(basis))
            ]
            or [0.0]
        )

        full_shape = [a.size for a in channel_axes]
        self.full_shape = full_shape
        self.cell_reshape = [
            a.size if a.name in self.cell_axis_names else 1 for a in channel_axes
        ]
        self.order_reshape = [
            [a.size if a.name in self.order_axis_names[k] else 1 for a in channel_axes]
            for k in range(self.order)
        ]
        self.shape_reshape = [
            a.size if a.name == shape_axis else 1 for a in channel_axes
        ]

        self.npou = 0
        self.allowNegativeParam = True
        self.is_linear = False
        self.xparamdefault = tf.zeros([self.npoi], dtype=indata.dtype)

    def compute(self, param, full=False):
        p_reshaped = [
            tf.reshape(self.p_k[k], self.shape_reshape) for k in range(self.order)
        ]
        rnorms = []
        for key, v in self.indata.channel_info.items():
            if v["masked"] and not full:
                continue
            nbins_channel = int(np.prod([a.size for a in v["axes"]]))
            irnorm = tf.ones(
                [nbins_channel, self.indata.nproc], dtype=self.indata.dtype
            )
            if key == self.channel:
                for i, proc_idx in enumerate(self.proc_idxs):
                    off = i * self.stride
                    a = tf.reshape(param[off : off + self.n_cell], self.cell_reshape)
                    expo = a
                    pos = off + self.n_cell
                    for k in range(self.order):
                        n = self.n_order_groups[k]
                        ck = tf.reshape(param[pos : pos + n], self.order_reshape[k])
                        expo = expo + ck * p_reshaped[k]
                        pos += n
                    scaling = tf.reshape(
                        tf.broadcast_to(tf.exp(expo), self.full_shape), [-1, 1]
                    )
                    proc_col = tf.one_hot(
                        proc_idx, self.indata.nproc, dtype=self.indata.dtype
                    )
                    irnorm = irnorm + (scaling - 1.0) * tf.reshape(proc_col, [1, -1])
            rnorms.append(irnorm)
        return tf.concat(rnorms, axis=0)

'''
