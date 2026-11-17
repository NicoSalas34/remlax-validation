# Building the manuscript

    pdflatex remlax && bibtex remlax && pdflatex remlax && pdflatex remlax

`remlax.tex` uses the Journal of Statistical Software class, `jss.cls`, with
`jss.bst`, in its `nojss` (preprint) form. Neither file is vendored here; both
come from the journal's author package or from CTAN (`tlmgr install jss`).

The built PDF is not versioned, so that the manuscript read from a clone is
always the one in `remlax.tex`.

## Figures

- `figures/fig_validation.png`: figure of the validation campaign of
  2026-09-28. That campaign is kept in the git tag `archive/2026-11-17`. The
  figure of the latest campaign is
  `correctness/reference_software/latest/2026-11-04_remlax-3bfae58/fig_validation.png`.
- `figures/fig_temps_n.pdf`, `fig_temps_p.pdf`, `fig_temps_C.pdf`,
  `fig_cpu_gpu.pdf`, `fig_remlax_asreml.pdf`, `fig_regle.pdf`: the figures of
  `speed/absolute/latest/2026-11-04_remlax-v14/`, regenerated without figure titles by
  `tools/synthese_abs_en.py`. This script is `speed/absolute/scripts/summarise.py`
  with its figure titles removed and the figure file names above kept, since
  `remlax.tex` includes them. The data logic is unchanged. Both scripts read the
  stored results, whose file and column names are French, and write
  `measurements.csv` with English columns (`case`, `hardware`, `route`,
  `seconds_median`, `n_failures`, `diff_best`, ...). Its rows are those of
  `speed/absolute/latest/2026-11-04_remlax-v14/measurements.csv` once columns and route names are translated
  (`remlax_froid` to `remlax_cold`, `mme_creux` to `mme_sparse`, `mme_hybride`
  to `mme_hybrid`).

      python3 tools/synthese_abs_en.py ../speed/absolute/latest/2026-11-04_remlax-v14/inputs <output dir>

The numbers in the text come from `speed/absolute/latest/2026-11-04_remlax-v14/measurements.csv`
and from the remlax repository (README, NEWS, vignettes, `dev/notes.md`).
Other numbers were read from runs now kept only in the git tag
`archive/2026-11-17`: `validation/results/` (campaigns up to 2026-10-30),
`verif_simulations/resultats/` and `docs/rapport_nuit_2026-10-27.md`. Their
latest counterparts are `correctness/reference_software/latest/`,
`correctness/simulations/latest/` and `speed/real_models/latest/`. Open
points are marked `%% TODO` in `remlax.tex`.
