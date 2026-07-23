# Design: AxisNormModel

## Context

`SignalNormModel` works by:
1. Reading `signal_process` from tensor metadata.
2. Assuming the channel has axes `[mass, pt, eta, charge]` and treating
   everything after index 0 as "cell axes".
3. Iterating over the Cartesian product of cell-axis bins to build one POI name
   per cell (`signorm_pt0_eta0_charge0`, …).
4. In `compute()`, reshaping the POI vector to `[1, n_pt, n_eta, n_charge]`,
   broadcasting it across the mass axis, then applying it only to the signal
   process column via a one-hot mask.

The goal is to generalize all four steps so that the channel, the target
processes, and the axes are provided via CLI rather than hardcoded.

## Goals / Non-Goals

- **Goal**: configurable channel, process selection, and axis list.
- **Goal**: always exactly 3 required positional args — no optional args, no
  ambiguity.
- **Goal**: add the class to `poi_model.py` alongside existing built-ins,
  consistent with the file's established role.
- **Goal**: register in `helpers.py` so the class is loadable by short name.
- **Non-goal**: a `signal` keyword default for proc_spec. Callers always name
  the processes or pass `all`.
- **Non-goal**: multiple-channel simultaneous operation. The model is scoped to
  one named channel per instance; two channels can be handled by
  `CompositePOIModel`.
- **Non-goal**: replacing `Mu`, `Ones`, or `SaturatedProjectModel`.

## Decisions

### CLI argument shape

**Decision**: always 3 positional args:
```
--poiModel AxisNormModel <channel> <proc_spec> <axes>
```
Examples:
```
# btojpsik: explicit signal process name
--poiModel AxisNormModel btojpsik_stuff BuToJpsiK_2018 bkmm_jpsimc_kaon1pt,bkmm_jpsimc_kaon1eta,bkmm_kaon_charge

# scale all processes
--poiModel AxisNormModel btojpsik_stuff all bkmm_jpsimc_kaon1pt,bkmm_jpsimc_kaon1eta,bkmm_kaon_charge

# multiple named processes
--poiModel AxisNormModel mychannel proc1,proc2 ax1,ax2
```

`parse_args` always expects `args = (channel, proc_spec, axes_csv)` and raises
`ValueError` if not exactly 3.

**Rationale**: removes all ambiguity. No 2-arg shorthand, no `signal` keyword,
no length-inspection in `parse_args`.

### Process resolution

| `proc_spec` | Resolved targets |
|---|---|
| `"all"` | all processes in the channel |
| `"proc1,proc2,..."` | the named processes; `ValueError` if any is absent |

Each matched process receives its **own independent** set of per-cell POIs.
`npoi = n_procs * n_cell`. POI names are prefixed with the process name:
`norm_{procname}_{ax1}{i}_{ax2}{j}...`. Two processes passed together are
equivalent to composing two single-process models — no POI is shared.

### Axis lookup

Axes are matched by `.name` against `channel_info[channel]["axes"]`. A clear
`ValueError` is raised if a requested axis name is absent.

### POI naming

Names follow `norm_{ax1name}{i}_{ax2name}{j}...` (bytes). The `signorm_`
prefix is dropped because the model is no longer restricted to signal processes.

### Broadcast in `compute()`

The POI tensor is reshaped to the full axis layout of the channel (size 1 on
axes not in the requested set, full bin size on axes that are), then
`broadcast_to` expands it to the complete channel shape before flattening.
This is identical to `SignalNormModel`'s approach.

Return shape is `[nbins, nproc]`:
- target process columns carry the per-cell POI value
- all other process columns are 1.0

## Risks / Trade-offs

- If the user requests the mass axis (or another axis that should not vary
  per-bin), nothing in the model prevents it. Convention remains user
  responsibility.
- Passing `"all"` when the tensor has many processes creates a large POI
  vector. This is intentional — each process needs its own normalization
  freedom — but the user should be aware of the scale.

## Open Questions

- None.
