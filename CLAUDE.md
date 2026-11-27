# CLAUDE.md

fortransitions is a Fortran (fpm) port of the Python package
[CosmoTransitions](https://github.com/clwainwright/CosmoTransitions).
It computes bounce actions for vacuum tunnelling. So far `fullTunneling`
(path deformation) and `SingleFieldInstanton` (1d overshoot/undershoot)
are ported, plus everything they depend on. See `README.md` for the
user-facing description.

## Build and test

```sh
fpm build --profile release
fpm test  --profile release              # all tests
fpm test  --profile release check        # a single test program
fpm run   --profile release              # app/full_tunneling_example.f90, writes CSVs
fpm test  --profile debug                # bounds checks, FPE traps (gfortran)
```

- **Always pass `--profile debug` or `--profile release`.** Without one,
  fpm uses no flags from `fpm.toml`. ifx then defaults to
  `-fp-model=fast`, which breaks the thin-wall shooting.
- gfortran, ifx (`--compiler ifx`) and LFortran (`--compiler lfortran`)
  are supported. The flags for each are under `[features]` in `fpm.toml`.
  Run the tests with all three after changes in `src/` when possible.
- LFortran: the profiles need `--cpp --legacy-array-sections
  --realloc-lhs-arrays`; release uses `--fast --no-fast-math` (fast math
  breaks `ieee_is_finite` and NaN comparisons). Debug builds need
  `llvm-dwarfdump` on the PATH (for `-g`).
- LFortran and the Python interface (`make wrapper FC=lfortran`):
  `python/fpm.toml` adds `--separate-compilation` (otherwise `lfortran -c`
  writes empty object files and the shared libraries export nothing), and
  the Makefile escapes `$ORIGIN` twice, because LFortran runs its link
  command through a shell that would expand it.
- LFortran bug (see `LFORTRAN_ISSUES.md`): a value assigned through an
  argument is lost if the call also passes an internal procedure and the
  variable is used by another internal procedure of the same host.
  Don't pass such variables as output arguments in calls that pass
  internal procedures (e.g. to `rkqs_vec`, `brentq`); use a separate
  local variable and copy it, as in `sfi_integrate_profile`.
- The debug profile traps invalid/zero/overflow FPEs under gfortran. Code
  must not rely on producing NaN/Inf silently.
- `*.csv` is gitignored, so example and test output is never committed.

## Layout

| File                    | Module                            | Python counterpart                                   |
| ----------------------- | --------------------------------- | ---------------------------------------------------- |
| `src/fortransitions.f90`| `fortransitions`                  | top-level module, re-exports the public API only     |
| `src/config.f90`        | `fortransitions__config`          | `wp`, `pi`, status codes (`status_ok`, `err_*`)       |
| `src/potentials.f90`    | `fortransitions__potentials`      | abstract `potential_1d`/`potential_nd`, `pot1d_func`  |
| `src/tunneling1d.f90`   | `fortransitions__tunneling1d`     | `tunneling1D.SingleFieldInstanton` → `single_field_instanton`, `profile1d` |
| `src/splinepath.f90`    | `fortransitions__splinepath`      | `pathDeformation.SplinePath` → `spline_path`         |
| `src/deformation.f90`   | `fortransitions__deformation`     | `pathDeformation.Deformation_Spline` → `deformation_spline` |
| `src/pathdeformation.f90`| `fortransitions__pathdeformation`| `pathDeformation.fullTunneling` → `full_tunneling`, `full_tunneling_result` |
| `src/helpers.f90`       | `fortransitions__helpers`         | `helper_functions` (Nbspld2, deriv14, cumtrapz, simpson, lstsq, ...) |
| `src/optimize.f90`      | `fortransitions__optimize`        | scipy `brentq`, `fminbound`, `fmin` replacements     |
| `src/special.f90`       | `fortransitions__special`         | scipy fractional-order Bessel `iv`, `jv`             |
| `src/c_potentials.f90`  | `fortransitions__c_potentials`    | potentials given by C function pointers (Python callbacks) |
| `src/c_tunneling1d.f90` | `fortransitions__c_tunneling1d`   | `bind(C)` interface to `single_field_instanton` (opaque handles) |
| `src/c_pathdeformation.f90`| `fortransitions__c_pathdeformation` | `bind(C)` interface to `full_tunneling` (opaque handles) |

- `test/`: `check.f90`, `double_well.f90` and `many_fields.f90` each
  assert against reference values from the Python package. They count
  failures and end with `error stop 1` if any test fails.
  `c_interface.f90` checks that the `bind(C)` layer gives results
  identical to the native API.
- `reference/`: Python scripts (`ref_*.py`) that generate those
  reference values with the original cosmoTransitions. See
  `reference/README.md` for the package versions and the reasoning
  behind each test's tolerance. The scripts need `cosmoTransitions`,
  numpy and scipy.
- `app/full_tunneling_example.f90`: the two-field example from
  `examples/fullTunneling.py`.

- `python/`: the Python package `fortransitions` (see
  `python/README.md`). `python/fpm.toml` builds `../src` as a shared
  library. The `Makefile` holds all of the build logic (`make lib`,
  `make wrapper`, `make clean`), and the `pixi.toml` tasks (`install`,
  `test`, `clean`) only call it, so the pixi and manual builds stay in
  sync. The rpath (`$ORIGIN`, plus the env's `lib/` under pixi) is set at
  link time, so no patchelf step is needed. The ctypes wrappers are in
  `python/fortransitions/` (`_lib.py` declares the C signatures), and
  the tests are in `python/test/run_tests.py`.
  `python/test/compare_cosmoTransitions.py` compares the results directly
  against cosmoTransitions; run it in the pixi environment `compare`
  (`pixi run -e compare compare`). Only Linux is supported.

External fpm dependencies: `odeint` (`rkqs_vec` Cash-Karp RK stepper),
`bspline-fortran` (`bspline_1d`, used for all splines), and `gradmin`
(`gauss_jordan` linear solve). The example program writes its CSV files
with plain `write` statements, so it needs no CSV library.

## Conventions

- Module names are `fortransitions__<file>` (`module-naming = false`).
  Every module is `implicit none` + `private` with explicit `public`
  lists. Each imported symbol gets its own line:
  `use mod, only : name`.
- Any new public API must also be re-exported from `src/fortransitions.f90`.
- Use `real(wp)` throughout. `wp = real64` and must match
  bspline-fortran's precision. Write literals as `1.0_wp`.
- Indent with 2 spaces. Document with FORD-style `!!` comments placed
  after the declaration they describe.
- Classes are derived types with type-bound procedures. Implementations
  carry a short prefix: `sfi_*` (single_field_instanton), `dfs_*`
  (deformation_spline), `spline_path_*`. The object argument is `self`.
- Optional arguments mirror the Python keyword arguments and their
  defaults. They are resolved into local `name_` variables.
- **No exceptions and no `stop` in library code.** Errors go out through
  an `integer, intent(out) :: status` argument using the codes in
  `fortransitions__config`, which mirror `PotentialError`,
  `IntegrationError` and `DeformationError`. Callers check `status` and
  propagate it. When you add a new failure mode, reuse an existing code
  or add one to `config.f90` and re-export it.
- Users supply potentials by extending `potential_nd` (deferred `v`,
  `grad`) or `potential_1d` (deferred `v`; `dv`/`d2v` default to finite
  differences). Their methods take `class(...), intent(inout) :: self`.
  `potential_nd%grad_many(x, dv)` evaluates the gradient at all points of
  an `(npts, ndim)` array. It loops over `grad` by default and can be
  overridden with a vectorized version. Code that needs the gradient at
  many points at once (the path deformation: `dfs_forces`, the `v2min`
  loop in `dfs_init`) must call `grad_many`, not `grad` in a loop.

## Python interface

- Build and test with `cd python && pixi run test` (or `--quick` via
  `pixi run python test/run_tests.py --quick`). The manual route is
  `make wrapper` in an environment with numpy. Rebuild after every
  change in `src/`, because the package loads
  `python/fortransitions/lib/libfortransitions.so`.
- The pixi environments and a manual install in another Python (e.g. a
  pyenv env) share `python/fortransitions/lib/`: the last `make lib`
  wins. The Makefile therefore uses one fpm build directory per toolchain
  (`build/<FC>-env-<pixi env>`, `build/<FC>-system`). fpm does not relink
  when only the link flags (rpath) change. A gfortran build runs in all
  of the environments (it needs only libgfortran.so.5 with GFORTRAN_10
  symbols), so no rebuild is needed when switching. Use the build of the
  toolchain under test when that matters.
- The pixi environment `ifx` (conda-forge `ifx_linux-64`) provides ifx
  and sets `FC=ifx`: `pixi run -e ifx test` tests the ifx build, and
  `pixi run -e compare-ifx compare` (features compare + ifx) compares it
  with cosmoTransitions. That
  library links Intel's runtime from `python/.pixi/envs/ifx/lib`, so it
  works in other environments only while the `ifx` environment exists.
- `make wrapper` installs in setuptools' strict editable mode, i.e. a
  link tree in `python/build/__editable__*` on a plain `sys.path` entry.
  The default import-hook mode would let the repository folder (also
  named `fortransitions`) shadow the package as an empty namespace
  package when Python runs in the repository's parent directory. New
  files in the package need a reinstall (`make wrapper`).
- The C interface takes no optional arguments. Python passes every
  value explicitly, with defaults identical to the Fortran ones, and
  "find automatically" options get a `has_*` flag. Any change to the
  arguments of `full_tunneling`, `single_field_instanton%init` or
  `%find_profile` must be carried through to `src/c_*.f90`,
  `python/fortransitions/_lib.py` and the Python wrappers.
- `c_potential_nd` overrides `grad_many` with one batched callback, so
  that Python evaluates `dV` on all path points in one numpy call. This
  is what makes the path deformation fast from Python. The Python side
  (`_Vectorization` in `_callbacks.py`) checks on the first call whether
  the user's `dV` accepts `(npts, ndim)` arrays, and falls back to
  point-by-point evaluation otherwise (`fullTunneling(vectorized=...)`).
- Python exceptions in callbacks are stored, and the callbacks return
  NaN from then on. The Fortran code must therefore end promptly on NaN
  input: write loop guards NaN-safe (`.not. (x <= limit)`), and return
  `err_numerical` on non-finite values instead of looping.
- `error stop` in library code would also kill the Python process.
  The Python wrappers validate the arguments that would reach the
  existing `error stop`s (e.g. `maxiter`, `V_spline_samples`, `nb`/`kb`).

## Porting guidelines

- Stay faithful to the Python algorithm, including its defaults and
  tolerances, so that results can be compared against the Python
  package. Where the port deliberately differs, document it in the doc
  comment and in the README's "Ported so far" section. Examples: the
  `v_mode = "direct"` path, the per-dimension 1d splines that remove
  FITPACK's 10-field limit, and the missing `save_all_steps`.
- When porting a new piece: add a test program in `test/` that asserts
  against values from a new `reference/ref_<name>.py`, list the script
  in `reference/README.md`, and update the README.
- Use relative tolerances in tests. The default
  `xtol = phitol = 1e-4` leaves errors of a few 1e-3 in the action
  (hence `rtol = 5e-3` in `check.f90`). Tighter tolerances (1e-9) give
  agreement to about 1e-7. Path-deformation trajectories can still
  differ slightly between compilers and optimization levels.
- `res%action_pot` and `res%action_kin` (Derrick-rescaled) and
  `res%fratio` are built-in quality diagnostics. Use them when checking
  whether a change has degraded accuracy.
- Before calling a change done, run the tests under both profiles. Also
  run them with ifx if it is available.
