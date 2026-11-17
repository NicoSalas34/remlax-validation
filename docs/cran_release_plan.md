# remlax to CRAN: release branch and freeze conditions

Note of 2026-10-29, on main at commit 777435a. It extends the audit
`cran_policy_audit.md` of 2026-09-29 (commit 48a75d2), made before the
arrival of the MME engine, AI-REML, the Python server, the Bayesian route and
the spatial lag.

## State of the package

`R CMD check --as-cran` of 2026-10-29, R 4.5.3, with the Python solver (RX_PY
set in `~/.Renviron`):

- no error;
- one warning, `qpdf` missing from the machine, unrelated to the package;
- one clock note ("future file timestamps"), specific to the machine;
- one note on the `rx_bayes` example, 5,1 s, and 32 s with `--run-donttest`;
- tests: 890 passed, 1 skipped, 197 s elapsed;
- vignettes: rebuilt without error, their code is not run.

The check without Python, the one the CRAN machines will see, remains to be
redone. RX_PY must be neutralised in the environment of the check itself, for
instance with `R_ENVIRON_USER=/dev/null`.

## Blockers since the audit

### Writes to the home directory

By default the solver writes to `~/.cache/remlax`:

- the XLA compilation cache (`device.py`, `activer_cache_compilation`);
- the compiled library of the C kernel of the selected inversion
  (`selinv.py`, `libselinv_<hash>.so`).

Today's check wrote there. The CRAN policy forbids any write outside the
temporary directory of the session, except in `tools::R_user_dir()`, provided
the content stays small and managed. The R side passes no `REMLAX_CACHE` to
the solver today.

Proposed fix:

- R passes `REMLAX_CACHE = tools::R_user_dir("remlax", "cache")` at each
  launch of the solver;
- during a check (`_R_CHECK_LIMIT_CORES_` set, or a non-interactive session
  without consent), the cache goes to `tempdir()`;
- a function `rx_cache_clear()` empties this folder and reports its size.

### Compiling C code at run time

`selinv.py` compiles `selinv_c.c` on first use, with the system compiler,
then loads the library through `ctypes`. A CRAN reviewer will most likely
refuse code compiled outside the installation.

Two routes are possible:

- move the kernel to `src/`. R compiles it at installation into `libs/`, with
  its own toolchain. R then passes the path of the library to the solver,
  which Python loads through `ctypes` as today. The BLAS and LAPACK pointers
  remain those of SciPy, passed at the call;
- disable the kernel in the CRAN version. The NumPy loop takes over, with the
  same result, but 5 to 15 times slower on the selected inversion (0,3 s
  against 3,9 s on a pedigree of 20 031 individuals).

I recommend the first route. It keeps the speed, and compilation goes through
the normal mechanism of an R package.

### Example too long

The `rx_bayes` example exceeds 5 s. It moves to `\donttest`.

## Points to check, not blocking

- Persistent Python server. It is disabled during a check and stopped by a
  finaliser at the end of the session. It remains to check on Windows that no
  process outlives R.
- Compute threads. `fils.py` restricts the affinity of the process. During a
  check, a single thread is used, as the policy requires. This must be
  rechecked on win-builder.
- Platforms. The package has been checked under Linux only. win-builder
  (R-devel and R-release), mac-builder and R-hub remain to be done. The audit
  had fixed `system2(env =)`, ignored under Windows, without testing it.
- DOI of Salas et al. (2026). The bioRxiv prefix 10.64898 resolves correctly
  (checked on 2026-10-29).
- Maintainer. The institutional address must remain valid after the PhD. The
  ORCID and the `cph` role of a possible employer remain to be decided.

## Decisions to make before freezing the API

Once on CRAN, renaming a function or an argument requires a deprecation
cycle. The following choices are therefore made before the freeze.

### Language of the interface

The interface mixes English and French:

- arguments in French: `marqueurs`, `carte`, `codage`, `construire`,
  `K_de`, `bloc`;
- one function in French: `rx_scan_seuil()`;
- output fields in French: `moteur`, `secondes`, `composantes_degenerees`,
  `signe_non_identifie`, `second_ordre_omis`;
- solver messages and printouts in French ("Ajustement REML",
  "parametres").

CRAN does not forbid it. A public package will still be read mostly in
English. I propose to switch the interface to English before the freeze, and
to keep the old names as deprecated aliases for one version. The internal
solver messages can stay in French at first.

### Scope of the first version

Some functions are recent and have seen little use outside our data: the
Bayesian route, the spatial lag, the block engine, the hybrid MME engine,
`rx_scan_loco()` and `rx_scan_exact()`. Two options are open:

- publish them, marked experimental in their help, with no promise of a
  stable interface;
- keep them on main, out of the CRAN version, until the next one.

I recommend the first option, with the experimental label. Removing code from
a release branch creates a divergence that is costly to maintain.

## Branch strategy

The `cran-layout` branch has already been merged into main, which is now
eight commits ahead of it. Maintaining a layout branch parallel to main cost
this merge. CRAN compliance must live on main, and `.Rbuildignore` already
excludes what does not go into the archive (`manual/`, `tests-python/`, the
benchmark archive).

Proposal:

- main stays the development branch, at version `0.1.0.9000` after the first
  submission;
- at the freeze, a `release/0.1` branch starts from a tagged commit. It
  receives only fixes requested by CRAN;
- the `v0.1.0` tag marks the submitted archive. Each fix made on
  `release/0.1` is carried over to main;
- the `cran-layout` and `structural-lag` branches, fully contained in main,
  can be deleted.

## Continuous integration

The repository has no `.github/` today, although `.Rbuildignore` anticipates
one. Two families of GitHub Actions jobs are needed:

- an `--as-cran` check without Python, under Linux, Windows and macOS, with
  R-release and R-devel (`r-lib/actions/check-r-package`). It reproduces what
  CRAN will see;
- a Linux job with Python and JAX on CPU, which runs the testthat tests with
  the solver and the pytest suite.

## Freeze conditions

The package is ready for a release branch when:

- the blockers above are fixed on main;
- the interface is frozen, language included;
- the tests added on 2026-10-29 (`test-scan-plus.R`, `test_face_cone.py`)
  are committed, and the full suite passes with and without Python;
- `verif_simulations` (remlax-validation) passes without a new failure;
- the `--as-cran` check gives 0 errors, 0 warnings and only the "New
  submission" note on win-builder, mac-builder and continuous integration;
- `NEWS.md` carries a short summary of the first version, and the current
  development log is moved elsewhere;
- `cran-comments.md` is up to date.

## Proposed order

The blockers on main come first: cache, C kernel, example too long. The
language and the scope are decided next, then the renaming is done.
Continuous integration is installed. The freeze follows: tag, `release/0.1`
branch, external checks. The submission comes last.
