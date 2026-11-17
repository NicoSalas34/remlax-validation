# Studies

Method studies built on remlax. They use the package; they do not test it.

| Folder | Study |
|---|---|
| [cogwas/](cogwas/) | simulated wheat and alfalfa mixtures for co-GWAS; power of one-step and two-step GWAS for additive DGE and IGE QTL ([report of 2026-11-02](cogwas/report_additive_qtl_2026-11-02.md), figures in `figures_additive_qtl/`) |
| [demos/](demos/) | demonstrations on simulated field trials: 2D P-spline surface with direct and indirect effects (`spatial.R`), one-step GWAS with LOCO (`gwas_dge_ige.R`), neighbourhood kernels and identifiability of the dilution (`identifiability_d.R`, `ring_kernel.R`, `two_species_kernel.R`, `weighted_kernel.R`, `replicates_d.R`, `no_radius.R`); stored outputs in `output/` |
| [genome_scan/](genome_scan/) | the genome scan of remlax against rrBLUP and statgenGWAS, and the families of tests |
| [neighbor_gwas/](neighbor_gwas/) | multi-trait Neighbor GWAS on simulated data |
| `pathosystem/` | work in progress, not versioned |

The stored outputs of `demos/output/` predate the English pass of the scripts
and keep French file names.
