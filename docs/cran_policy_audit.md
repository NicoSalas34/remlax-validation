# remlax 0.1.0 : audit against the CRAN Repository Policy (revision 6875) and the submission checklist

Branch `cran-layout`, commit 48a75d2. Local checks: R 4.5.3, `R CMD check --as-cran --timings`, with and without the Python solver.

Check results: 0 errors, 0 warnings. Notes: "New submission"; three GitHub URLs in 404 (private repositories); "unable to verify current time" (local clock check only).

| Rule | Status | Evidence or action |
|---|---|---|
| Submission form: DESCRIPTION present, Maintainer "NAME <EMAIL>", tarball from R CMD build | met | remlax_0.1.0.tar.gz, 500 KB |
| Check with --as-cran on R-devel | to do | only R 4.5.3 here; run win-builder R-devel before submitting |
| No warnings or significant notes | met | only the three notes above |
| Single designated maintainer, a person, reachable address | to confirm | nicolas.salas@supagro.fr: is it durable? |
| Copyright and authorship clear; cph role if an institution holds it | to confirm | only aut/cre for N. Salas; employer as cph? |
| ORCID in Authors@R (checklist) | to do | ORCID not provided |
| Licence in license.db, perpetual distribution | met | GPL (>= 3) |
| No binary executable code | met | archive holds .R, .Rd, .py, .csv, .json, .Rmd only |
| Package name unused on CRAN (current, archive) and Bioconductor | met | none of the three lists holds remlax |
| Strong dependencies on CRAN; no restricted software required | met | Imports: methods, Matrix, jsonlite, stats, tools, utils; asreml is not a dependency, its results are precomputed |
| Suggests used conditionally | met | RTMB/TMB via rx_tmb_available(), lme4/nlme/sommer via requireNamespace() |
| External software declared | met | SystemRequirements: Python >= 3.10 with jax, numpy, scipy |
| Informative Description, author-year citations with <doi:...>, software names in single quotes, functions as foo() | met | Patterson and Thompson (1971), Salas et al. (2026), 'JAX', 'Python', 'ASReml-R', 'RTMB' |
| Cross-platform (two major platforms) | improved, to verify | system2(env=) is ignored on Windows: variables now set by Sys.setenv() and restored; to test on win-builder |
| Size: tarball < 10 MB, data and docs < 5 MB | met | 500 KB |
| Check CPU time as small as possible | met | examples 15 s, tests 98 s CPU with Python; without Python the solver tests are skipped |
| At most two cores during checks | fixed | one computation thread when _R_CHECK_LIMIT_CORES_ is set (measured 1.2 cores; tests CPU/elapsed 1.5) |
| Examples run a few seconds each | fixed | every export has an example; rx_scan in \donttest (5.4 s); only rx_install_python stays in \dontrun (large download) |
| Internet resources fail gracefully | met | rx_install_python() stops with a message if venv or pip fail |
| No writing outside tempdir | fixed | Python wrote __pycache__ into the R library: now PYTHONDONTWRITEBYTECODE=1; bundles go to tempfile() |
| User data only in R_user_dir(), kept small and managed | fixed | rx_install_python() asks in interactive sessions; new rx_remove_python() |
| No modification of the global environment | met | rx_validate() now restores .Random.seed |
| No q(), no ::: on base packages, no tampering with loaded code | met | none found |
| No external software left running | met | the solver subprocess ends with each call |
| Messages suppressible (reviewer practice) | fixed | cat() only in print methods; solver output moved into stop(); rx_validate() uses message() |
| Secure downloads | met | pip uses PyPI over https |
