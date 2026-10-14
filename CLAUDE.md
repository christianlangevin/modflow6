# CLAUDE.md

Guidance for working in the MODFLOW 6 repository.

## What this is

MODFLOW 6 is the U.S. Geological Survey's modular hydrologic simulation
program. The simulation engine is written in **Fortran** (`src/`); the build
system is **Meson**; the test suite and developer tooling are **Python**,
managed with **Pixi**. Current version: `version.txt` (`6.8.0.dev0`).

Key reference docs in the repo root:
- `DEVELOPER.md` — authoritative build/test/contribute guide (long; consult for specifics)
- `IDM.md` — Input Data Model / DFN-driven input processing
- `CONTRIBUTING.md`, `EXTENDED.md` (PETSc/MPI/NetCDF builds)

## Environment (Pixi)

All Python tooling and the compiler toolchain come from the Pixi environment.
Run developer tasks via `pixi run <task>` (tasks defined in `pixi.toml`). The
first-time setup also needs:

```shell
pixi run install      # installs flopy, pymake, modflowapi from git
pixi run get-exes      # builds latest release in develop mode + downloads related exes into bin/
pixi run update-flopy  # regenerate flopy.mf6 classes from DFN files
```

## Building (Meson)

Two steps: configure a build dir, then install. Binaries land in `bin/`.

```shell
pixi run setup builddir              # configure (release); add -Ddebug=true for debug
pixi run build builddir              # meson install -> bin/
```

Raw meson equivalent: `meson setup --prefix=$(pwd) --libdir=bin builddir` then
`meson install -C builddir`. Extra meson flags (`-Ddebug=true`, `-Dextended=true`,
`--wipe`) pass through the pixi command. Build options are in `meson.options`.

After building, the executable is `bin/mf6` (`mf6 -v` prints version).

**Build compiler (this machine):** the pixi *default* environment does **not**
include a Fortran compiler, so a fresh `meson setup`/`--wipe` in it fails with
"Compiler gfortran cannot compile programs". Use the **`gcc-extended-build`**
pixi environment, which bundles a pixi-managed gfortran (and handles the macOS
SDK/sysroot automatically — no `ld: library not found` issues):

```shell
pixi run -e gcc-extended-build setup _builddir_pixi_gcc
pixi run -e gcc-extended-build build _builddir_pixi_gcc   # installs bin/mf6, bin/libmf6.dylib, bin/zbud6
```

This is the preferred, no-conda build path. Notes:
- First use installs the env, which pulls PETSc/OpenMPI/HDF5/NetCDF (heavy, ~GB,
  separate solve group) — overkill for the serial build but it's the only pixi
  env that ships a compiler. It gives **gfortran 13.4.0** (was 15.1.0 until
  develop's #3016 re-solved the environment for the conda-forge test-drive).
- **After a develop sync that changes `pixi.toml`/`pixi.lock`, wipe the build
  dir first**: a compiler-version change leaves stale `.mod` files and the
  build dies with "Cannot read module file ... created by a different version
  of GNU Fortran". Fix: `pixi run -e gcc-extended-build setup _builddir_pixi_gcc
  --wipe`, then build.
- A plain `build`/`meson install` (no source changes) only regenerates
  `version.f90` and reinstalls; a clean recompile needs `--wipe` at setup.
- Only *building* needs this env. The pixi `test`/`autotest` tasks run fine in
  the default env once `bin/mf6` exists.

Alternative compiler (gfortran 14.2.0): the conda `modflow6` env also works, but
its gfortran needs the env activated for the SDK sysroot:
`source ~/miniforge3/etc/profile.d/conda.sh && conda activate modflow6` before
`meson setup --prefix=$(pwd) --libdir=bin _builddir_Darwin_gfortran_release`.

## Testing

Two kinds of tests:

**Fortran unit tests** — `test-drive` framework, files named `autotest/Test*.f90`,
driven by Meson. Build first, then:

```shell
pixi run test builddir                    # all unit tests
pixi run test builddir --verbose ArrayHandlers   # one module (name from autotest/tester.f90)
```

**Python integration tests** — `pytest`, files named `autotest/test_*.py` (~390
scripts). Must run from `autotest/` if calling pytest directly; the pixi task
runs from root.

```shell
pixi run autotest                  # parallel (-n auto), full suite
pixi run autotest -S               # smoke: fast subset (-m "not slow and not regression")
```

Pytest markers (`pytest.ini`): `slow`, `external` (needs external model repos),
`large`, `regression`. Combine with `-m "not slow and not regression"`.

Integration tests use the `TestFramework` class (`autotest/framework.py`) with
`build`/`check`/`plot` hooks; programs under test are registered in the
`targets` dict (`autotest/conftest.py`). See `DEVELOPER.md` "Writing tests".

**Before running integration tests**, ensure the test prerequisites are in place
(see "Environment" above): `bin/mf6` is freshly built, `pixi run get-exes` has
populated `bin/rebuilt` + `bin/downloaded`, and flopy is synced.

**Common failure — stale flopy.** If tests die during model construction with
`flopy.mf6.mfbase.FlopyException: Extraneous kwargs "..."`, the installed flopy
predates a DFN change (a new input option the binary supports but flopy doesn't
know about). Fix with `pixi run update-flopy`. Always run this after editing any
`.dfn` and before running the suite. (A real smoke run: ~700 pass / ~480 skip in
~35s; the skips are the `external`-model cases, which `-S` does not pull.)

## Source layout (`src/`)

- `Model/` — the model types: `GroundWaterFlow` (GWF), `GroundWaterTransport`
  (GWT), `GroundWaterEnergy` (GWE), `ParticleTracking` (PRT),
  `ChannelFlow` (CHF), `OverlandFlow` (OLF), `SurfaceWaterFlow` (SWF)
- `Solution/`, `Exchange/`, `Distributed/`, `Timing/` — solver, model coupling, parallel, time stepping
- `Utilities/` — memory manager (`Memory/`), block parser, I/O, observations, time series, etc.
- `Idm/` — autogenerated `*idm.f90` input-definition modules (do not hand-edit)
- `mf6.f90` — program entry point
- `srcbmi/` — BMI/XMI shared-library interface

## Input Data Model (IDM) & DFN workflow — IMPORTANT

Input parameters are described in **definition files** (`*.dfn`) under
`doc/mf6io/mf6ivar/dfn/` (~147 files). These drive autogenerated Fortran modules
and flopy classes. **When you change a `.dfn` file you must regenerate** and
commit the generated artifacts alongside it:

```shell
pixi run update-fortran-definitions   # regenerate src/Idm/*idm.f90 (dfn2f90.py over dfns.txt)
pixi run run-mf6ivar                  # regenerate mf6io LaTeX input docs
pixi run update-flopy                 # regenerate flopy classes
```

A package is registered with IDM by adding its dfn filename to
`utils/idmloader/dfns.txt`. Generated Fortran definition modules are checked
into source control. See `IDM.md` for the full design.

## Formatting & linting (run before a PR)

```shell
pixi run prepare-pull-request    # fix-style + python-style + makefiles + mf6ivar + fortran defs
```

Individual pieces:
- Fortran: `fprettify` (config `.fprettify.yaml`; 82-col, 2-space indent) — `pixi run check-format` / `fix-style`
- Python: `ruff` (config `ruff.toml`; 88-col, double quotes) — `pixi run check-python-lint`, `fix-python-style`
- Spelling: `codespell` (`pixi run check-spelling`)

If a `utils/` program (e.g. `mf5to6`, `zbud6`) needs a new `src/` module, add it
to that util's `pymake/extrafiles.txt`.

## Branching / contributing

git-flow: develop on **`develop`** (the default branch and PR target); `master`
tracks the latest release. Development PRs are squash-merged. The full suite must
pass in CI before merge. Conventional-commit style messages (`feat(gwf):`,
`refactor(input):`, `chore(deps):`, etc.) — see recent git history.

Deprecating input options is done via `deprecated x.y.z` / `removed x.y.z`
attributes in DFN files (kept in perpetuity); see `DEVELOPER.md` "Deprecation policy".
