# Published data sets, state of 2026-11-01

All 2414 counted checks passed on the 15 data sets. The report is
[docs/published-data.md](../../../../docs/published-data.md).

The data sets were not all refitted with the same remlax. This folder name
therefore carries no version. Refitting the older ones with the current
remlax would bring the whole table to one version.

| remlax | Data sets |
|---|---|
| 0.2.0.9000 | maize_drops_multitrait, maize_technow2014, milk_pedigreemm, wheat_crossa2010, yeast_bloom2015 (4 traits) |
| 0.2.0 | pig_cleveland2012, tilapia_joshi2018, twins_openmx, warcolak_nadiv |
| 0.1.0 | cpdata_simulated, gryphon_wilson2010, maize_drops_gwas, mice_valdar2006, potato_tetraploid, rice44k_zhao2011 |

The version of each fit is printed at the top of its `.log`.

## Files

| File | Content |
|---|---|
| `<data set>.csv` | the checks of one data set |
| `<data set>.log` | the log of its script |
| `*_timing.csv`, `*_manhattan.png` | times and Manhattan plots of the genome scans |
| `fig_published_data.png` | every counted check as discrepancy over tolerance |
