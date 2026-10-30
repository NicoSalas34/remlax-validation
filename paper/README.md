# Building the manuscript

    pdflatex remlax && bibtex remlax && pdflatex remlax && pdflatex remlax

`remlax.tex` uses the Journal of Statistical Software class, `jss.cls`, with
`jss.bst`, in its `nojss` (preprint) form. Neither file is vendored here; both
come from the journal's author package or from CTAN (`tlmgr install jss`).

The built PDF is not versioned, so that the manuscript read from a clone is
always the one in `remlax.tex`.

## Figures

- `figures/fig_validation.png`: copy of `validation/results/fig_validation.png`
  (validation campaign of 2026-09-28).
- `figures/fig_temps_n.pdf`, `fig_temps_p.pdf`, `fig_temps_C.pdf`,
  `fig_cpu_gpu.pdf`, `fig_remlax_asreml.pdf`, `fig_regle.pdf`: the figures of
  `banc_kmmr/absolu/sortie_v14/`, regenerated without figure titles by
  `tools/synthese_abs_en.py`. This script is `banc_kmmr/absolu/synthese_abs.py`
  with its figure titles removed and the figure file names above kept, since
  `remlax.tex` includes them. The data logic is unchanged. Both scripts read the
  stored results, whose file and column names are French, and write
  `measurements.csv` with English columns (`case`, `hardware`, `route`,
  `seconds_median`, `n_failures`, `diff_best`, ...). Its rows are those of
  `sortie_v14/mesures.csv` once columns and route names are translated
  (`remlax_froid` to `remlax_cold`, `mme_creux` to `mme_sparse`, `mme_hybride`
  to `mme_hybrid`).

      python3 tools/synthese_abs_en.py ../banc_kmmr/absolu/resultats <output dir>

Every number in the text comes from `banc_kmmr/absolu/sortie_v14/mesures.csv`,
`validation/results/`, `verif_simulations/resultats/`,
`docs/rapport_nuit_2026-10-27.md`, or the remlax repository (README, NEWS,
vignettes, `dev/notes.md`). Open points are marked `%% TODO` in `remlax.tex`.
