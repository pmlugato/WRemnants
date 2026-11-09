<!-- OPENSPEC:START -->
# OpenSpec Instructions

These instructions are for AI assistants working in this project.

Always open `@/openspec/AGENTS.md` when the request:
- Mentions planning or proposals (words like proposal, spec, change, plan)
- Introduces new capabilities, breaking changes, architecture shifts, or big performance/security work
- Sounds ambiguous and you need the authoritative spec before coding

Use `@/openspec/AGENTS.md` to learn:
- How to create and apply change proposals
- Spec format and conventions
- Project structure and guidelines

Keep this managed block so 'openspec update' can refresh the instructions.

<!-- OPENSPEC:END -->



# Environment Setup

Commands in this repository must be run inside the apptainer/singularity container at:

```
/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/bendavid/cmswmassdocker/wmassdevrolling
```

or the corresponding docker/podman container at

```
gitlab-registry.cern.ch/bendavid/cmswmassdocker/wmassdevrolling
```

Check whether you are already inside the container by checking if the `SINGULARITY_CONTAINER` or $container environment variable is set.

If not inside the container, prefix commands with `scripts/ci/run_with_singularity.sh`:

```bash
scripts/ci/run_with_singularity.sh bash -c "source setup.sh && <command>"
```

If already inside the container, then commands can be run directly (but still sourcing setup.sh).

## Verification

To verify the environment is correctly set up, run:

```bash
python scripts/tests/testenv.py
```

Expected output ends with:
```
Successfully imported relevant packages
```
