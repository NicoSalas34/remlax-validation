# Correctness

Does remlax compute the right restricted likelihood, and the right estimates?
Each folder answers the question against a different kind of reference.

| Folder | Reference | Latest run |
|---|---|---|
| [reference_software/](reference_software/) | asreml, lme4, sommer, nlme, pbkrtest and closed forms, on the same model and data | 2026-11-04, remlax 3bfae58 |
| [simulations/](simulations/) | the simulated truth, and a dense REML arbiter written in R | 2026-11-04, remlax 3bfae58 |
| [published_data/](published_data/) | published data sets refitted with the programs of their papers | 2026-11-01, remlax 0.1.0 to 0.2.0.9000 |
| [golden/](golden/) | remlax itself before a refactoring | references of 2026-11-01 |

The reader's summary is [docs/validation-summary.md](../docs/validation-summary.md).
