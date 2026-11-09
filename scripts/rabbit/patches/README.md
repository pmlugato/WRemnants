# rabbit-side changes, kept outside the submodule

`rabbit/` is a git submodule. Its working tree gets stashed and reset when
upstream branches are pulled in — on 2026-08-27 a reset to
`havyn/upstream-main-havyn` removed `AxisExpPolyModel` mid-campaign, *after* six
fits had used it, and three queued fits then reported "done" in two minutes
having done nothing.

So every rabbit-side change this analysis depends on is stored here and
re-applied by an idempotent script, rather than living only in the submodule's
working tree.

    scripts/rabbit/patches/install_param_models.sh     # run after any submodule update

Contents:

| file | what it restores |
|---|---|
| `axis_exp_poly_model.py` | `AxisExpPolyModel` + its `helpers.py` registry entry |
| `rabbit_plot_cov_missing_param.patch` | `rabbit_plot_cov.py` naming a missing `--params` entry instead of raising `TypeError` |

The installer is a no-op when the changes are already present, so it is safe to
run before every campaign.

## Reconciling with `stash@{0}`

`stash@{0}` ("local axis models + plot_cov fix") holds the pre-reset copy of the
same three files, so `git stash pop` **will conflict** with the re-applied
versions. The stash additionally contains `AxisExpQuadModel`, which the current
version deliberately drops: `BKGMODEL=expquad` is now an alias for
`AxisExpPolyModel` order 2, verified equivalent (span residual 3.2e-07; Asimov
edmval 2.2338972927e-29 against 2.2338976236e-29). Nothing in the stash is
needed that this directory does not already reproduce.
