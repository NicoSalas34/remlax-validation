# Candidate public data sets: animals, humans and model organisms

Survey made on 2026-10-31 for the validation of remlax beyond additive models. Every fact below was checked on the cited primary source (article, supplementary file, repository page, CRAN documentation or package source). "Not verified" marks what could not be checked. "Computed locally" marks a value refitted during the survey; such values are not references. Small files opened during the survey sit in `/tmp/claude-1000/-home-salas-Documents-Git-IGE-analysis-2024-2025/162e877b-436f-40fc-bda3-3809b2eeed6b/scratchpad/datasets/` and will not persist.

## Ranking for non-additive models

| rank | data set | species | non-additive or structured terms | printed reference | exact target? |
|---|---|---|---|---|---|
| 1 | Bloom et al. 2015 BYxRM | yeast, 4390 segregants | A + A∘A (Hadamard) + repeatability | AI-REML fractions with SE, 20 traits | same data, same model |
| 2 | Cleveland et al. 2012 PIC | pig, 3534 animals | genomic A, D, A×A; 5 traits | 4 papers, 3 programs | plausibility band (filters and D coding differ) |
| 3 | Joshi et al. 2018 | Nile tilapia, 2524 fish | pedigree A + D + maternal; sire/dam diallel | ASReml-4 VCs and logLik, every model | yes (ASReml-4) |
| 4 | nadiv `warcolak` | simulated, 5400 | pedigree A + D + X-linked, sex-specific | simulated truth and realised effects | truth, refit against nadiv + asreml or gremlin |
| 5 | Bloom et al. 2013 BYxRM | yeast, 1008 segregants | 16 per-chromosome kernels; A + A∘A in Young & Durbin 2014 | Supp Table 1 and 2; GCTA fit | yes for the 16-kernel model |
| 6 | OpenMx `twinData`, mets `twinbmi` | human twins | ACE, ADE | ML estimates and -2lnL | yes, in ML mode |
| 7 | solarius `dat30` | simulated families | pedigree A, bivariate | SOLAR ML h2, rhog, logLik | yes, in ML mode |
| 8 | coxme GAW pedigrees | human families | kinship + IBD (two kernels), family effect | lmekin ML logLik | yes, in ML mode |
| 9 | Ober et al. 2015, 2012 | Drosophila DGRP | G + iid line + replicate; G∘G (2012) | SAS and ASReml tables | 2015 yes; 2012 only on line means |
| 10 | MCMCglmm `BTdata` | blue tit, 828 | dam + foster nest, multi-trait | sommer REML VCs in enhancer help | yes |
| 11 | pedigreemm `milk` | Holstein cows | repeated lactations (A + PE) | none verified | refit only |
| 12 | Robledo et al. 2018 | Atlantic salmon, 1480 | A, bivariate with rG near 1 | ASReml 3 h2 | yes |
| 13 | qtl2data DO, Gough F2 | mouse | LOCO GWAS; repeated weights | one h2 (Gatti 2014) | partial |

The tilapia, twin, solarius and coxme references give exact numbers on public data. The yeast and PIC references give the closest public tests of A∘A and genomic D. No public data set was found with a printed cross-software comparison (asreml, GCTA, sommer, BLUPF90) of A, D or A×A variances. The HIBLUP paper compares BLUPF90, ASReml, GCTA, MTG2, LDAK and sommer on time and memory only, on undistributed simulated data ([PMC10164590](https://pmc.ncbi.nlm.nih.gov/articles/PMC10164590/)).

## 1. Yeast BYxRM, Bloom et al. 2015

- Article: Bloom et al. 2015, Nat Commun 6:8712, [doi:10.1038/ncomms9712](https://doi.org/10.1038/ncomms9712), CC BY 4.0.
- Data: 4390 haploid segregants, 20 end-point growth traits, at least 2 replicates per segregant. Genotypes coded -1 (BY) / +1 (RM). The marker count is inconsistent in the text (28,820 and 28,220) and in the Supp Data 5 description (4390 × 28,880); not verified on the RData.
- Download: Supplementary Data 5, `genotypes.RData` + `phenotypes.RData`, 13,479,206 bytes, [MOESM729 zip](https://static-content.springer.com/esm/art%3A10.1038%2Fncomms9712/MediaObjects/41467_2015_BFncomms9712_MOESM729_ESM.zip). Replicates are stored as pairs in `pheno_raw`, so the fit can use individual records.
- Model printed: y = Za + Zi + Zp + e, a ~ N(0, A σ²A) with A from `rrBLUP::A.mat`, i ~ N(0, A∘A σ²AA), p iid strain repeatability. Custom AI-REML in R.
- Reference values: per-trait fractions of phenotypic variance with SE in [MOESM726 xls](https://static-content.springer.com/esm/art%3A10.1038%2Fncomms9712/MediaObjects/41467_2015_BFncomms9712_MOESM726_ESM.xls). Examples for the A + A∘A + P + E model:

| trait | A | A∘A | P | E |
|---|---|---|---|---|
| CobaltChloride | 0.426 (0.042) | 0.071 (0.014) | 0.108 (0.017) | 0.395 (0.011) |
| CopperSulfate | 0.461 (0.044) | 0.213 (0.018) | 0.164 (0.016) | 0.162 (0.004) |
| Hydroxyurea | 0.086 (0.013) | 0.023 (0.014) | 0.007 (0.021) | 0.884 (0.022) |
| Zeocin | 0.704 (0.055) | 0.047 (0.009) | 0.152 (0.009) | 0.097 (0.002) |

- The text quotes medians of 43.3 % (A), 9.2 % (A∘A) and 8.8 % (repeatability). Dominance is not identifiable in haploids.
- Refit with: asreml-R (`vm()` on two user matrices), sommer, gaston `lmm.aireml` with a kernel list, GCTA `--mgrm`.

## 2. PIC pig, Cleveland, Hickey & Forni 2012

- Article: G3 2:429, [doi:10.1534/g3.111.001453](https://doi.org/10.1534/g3.111.001453), CC BY 3.0. The article states that PIC made the data available for method comparison. No separate data licence is given.
- Download: the G3 supplement link returns 403, and the PMC link sits behind a JavaScript challenge. The Europe PMC bundle works: [supplementaryFiles](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3337471/supplementaryFiles), 67.4 MB, containing `supp_2.4.429_FileS1.zip` (67.2 MB).
- Content of FileS1 (CSV, CRLF): `genotypes.txt` (3534 animals × 52,843 SNPs coded 0/1/2, imputed values non-integer, 934 MB uncompressed), `pedigree.txt` (6473 IDs, 1247 founders, 0 = unknown), `phenotypes.txt` (t1 to t5, missing "."), `ebvs.txt`, `PIC Data Description.docx`.
- Non-missing records: t1 2804, t2 2715, t3 3141, t4 3152, t5 3184, matching Table 1 of the article. Phenotypes are already corrected for fixed effects; t1 and t2 are weighted progeny means. SNP order is randomised and there is no map.
- Pedigree h² from the full PIC data (Table 1): 0.07, 0.16, 0.38, 0.58, 0.62 for t1 to t5. The docx mislabels t5 as "t4=0.62".
- Models: genomic A, D, A×A, pedigree A, 5-trait multi-trait, GWAS without map. No maternal or PE term.

Printed non-additive values on these data (filters, subsets and D coding differ between papers):

| source | program | t1 | t2 | t3 | t4 | t5 |
|---|---|---|---|---|---|---|
| Da et al. 2014, Table 1, h²A / h²D | GVCBLUP | 0.03 / 7.2e-7 | 0.27 / 0.02 | 0.22 / 0.07 | 0.35 / 0.01 | 0.38 / 0.05 |
| Liu et al. 2022, Table 2, ADM model, VA/VP / VD/VP | DMU AI-REML | 0.033 / 0.000 | 0.263 / 0.015 | 0.213 / 0.066 | | |
| Da et al. 2022, Table 1, A / D / A×A | EPIHAP family, exact EGERM | 0.023 / 0.000 / 0.046 | 0.217 / 0.013 / 0.186 | 0.131 / 0.000 / 0.278 | 0.336 / 0.000 / 0.017 | 0.366 / 0.052 / 0.054 |

Sources: Da et al. 2014 PLoS ONE [doi:10.1371/journal.pone.0087666](https://doi.org/10.1371/journal.pone.0087666); Liu et al. 2022 Heredity [doi:10.1038/s41437-022-00504-6](https://doi.org/10.1038/s41437-022-00504-6); Da, Liang & Prakapenka 2022 Front Genet [doi:10.3389/fgene.2022.922369](https://doi.org/10.3389/fgene.2022.922369). Nishio & Satoh 2014 PLoS ONE ([doi:10.1371/journal.pone.0085792](https://doi.org/10.1371/journal.pone.0085792), Table 6) fitted t4 and t5 only, on 1800 randomly drawn animals. Their GBLUP-D values (t4: VA 1.735, VD 0.537, VE 3.331; t5: VA 1239.0, VD 220.8, VE 2049.3) cannot be reproduced exactly because the subset is not given.

Vitezica et al. 2013 (Genetics 195:1223) used heterogeneous-stock mice and simulation. Vitezica et al. 2017 (Genetics 206:1297) is simulation only. Neither uses the PIC data.

- Refit with: asreml-R with user G, D and G∘G via `vm()`, sommer (`A.mat`, `D.mat`, `E.mat`), GCTA (`--make-grm-d`; A×A as a user GRM), gaston multi-kernel, BLUPF90 with user G.

## 3. Nile tilapia, Joshi et al. 2018

- Article: Joshi, Woolliams, Meuwissen & Gjøen 2018, Heredity 120:452, [doi:10.1038/s41437-017-0046-x](https://doi.org/10.1038/s41437-017-0046-x), CC BY 4.0.
- Download: phenotypes [MOESM5 csv](https://static-content.springer.com/esm/art%3A10.1038%2Fs41437-017-0046-x/MediaObjects/41437_2017_46_MOESM5_ESM.csv) (181,614 bytes; 2524 fish; ID, Sire, Dam, batch, six traits BWH, BL, BD, BT, FW, FY, filleter, two dates, full-sib family, reciprocal-cross type). Pedigree [MOESM6 csv](https://static-content.springer.com/esm/art%3A10.1038%2Fs41437-017-0046-x/MediaObjects/41437_2017_46_MOESM6_ESM.csv) (61,233 bytes, 4050 IDs). Full variance tables in [MOESM3 pdf](https://static-content.springer.com/esm/art%3A10.1038%2Fs41437-017-0046-x/MediaObjects/41437_2017_46_MOESM3_ESM.pdf).
- Design: 155 full-sib families from a factorial mating with reciprocal crosses, pedigree from microsatellites.
- Model: fixed reciprocal-cross type, batch, collection day (and filleter for FW, FY). Random animal A, animal D (coefficients of fraternity, D⁻¹ from nadiv) and a maternal effect with one level per full-sib family. ASReml-4. Models marked with an asterisk re-base the pedigree at generation 20. Sire-dam diallel models are also reported.
- Reference values: Table 2 gives the ASReml logLik for every model and trait (BWH: A -28.34, AD -26.66, AM -24.15, ADM -22.19). Table S3.1 gives VA, VD, VM, VE, VP with SE; checked in the PDF for BWH:

| model | VA | VD | VM | VE | VP |
|---|---|---|---|---|---|
| A | 3216 [989] | | | 4312 [496] | 7528 [546] |
| AD | 2085 [1219] | 1688 [1090] | | 3673 [636] | 7447 [543] |
| ADM | 0 [0] | 1548 [614] | 637 [313] | 4739 [450] | 6925 [381] |
| AM | 1008 [492] | | 635 [316] | 5342 [287] | 6986 [405] |

- The ADM fit puts VA on the boundary, which tests constrained estimation. No genotypes. The genomic follow-up (Joshi et al. 2020 GSE, [doi:10.1186/s12711-019-0522-2](https://doi.org/10.1186/s12711-019-0522-2)) has genotypes under a data transfer agreement and is unusable.
- Refit with: asreml-R with `nadiv::makeD`, sommer, gremlin, MCMCglmm, WOMBAT, BLUPF90.
- Refit on 2026-10-31 (`published_data/scripts/12_tilapia_joshi2018.R`): the printed values are reproduced only with a maternal effect per dam (18 dams), not per full-sib family, and with the day as a numeric covariate beside batch (Day1 for BWH, BL, BD; Day2 for BT; Day2 + filleter for FW, FY). The Table 2 log-likelihoods equal the REML value (asreml convention) plus an exact multiple of 100 (200 for BD, 12200 for BWH). D⁻¹ is 96 % dense on this pedigree.

## 4. nadiv `warcolak` and related R package data

- nadiv 2.18.0 on [CRAN](https://cran.r-project.org/package=nadiv), GPL. Paper: Wolak 2012, Methods Ecol Evol 3:792, [doi:10.1111/j.2041-210X.2012.00213.x](https://doi.org/10.1111/j.2041-210X.2012.00213.x). Functions `makeD` (coefficients of fraternity, diagonal 1-f), `makeAA` (A∘A from the pedigree), `makeDomEpi` (AD, DD), `makeS`/`makeSd` (sex-linked) ([reference manual](https://search.r-project.org/CRAN/refmans/nadiv/html/makeD.html)).
- `warcolak` ([doc](https://search.r-project.org/CRAN/refmans/nadiv/html/warcolak.html)): 5400 individuals from 75 sets of the Fairbairn & Roff three-generation design, built by `simPedDFC`. Two uncorrelated traits. Simulated VA 0.4, VD 0.3, VR 0.3; trait2 splits VA into autosomal 0.3 and X-linked 0.1. The realised effects (`t1_a`, `t1_d`, `t2_s`, ...) are columns of the data. The help page holds the generating code. No published fitted estimates were found.
- gremlin 1.1.0 ([CRAN](https://cran.r-project.org/package=gremlin), GPL-3, Wolak) is an AI-REML program for animal models with generalised inverses. Its [README](https://github.com/matthewwolak/gremlin) fits `trait1 ~ sex - 1` with random `ID + IDD` and `ginverse = list(ID = Ainv, IDD = Dinv)` on `warcolak`, without printed output. gremlin is a further reference program for pedigree D.
- Mrode textbook data in R packages are pedigrees only, apart from small weaning-gain columns: nadiv `Mrode2`, `Mrode3` (with WWG and genetic groups), `Mrode9` (12 pigs, adapted from Example 9.1, dominance), AGHmatrix `ped.mrode`, gremlin `Mrode11` (5 records). The book solutions use known variance ratios. They test the MME solve at fixed variances, not REML. The Example 9.1 phenotypes are not shipped.

## 5. Yeast BYxRM, Bloom et al. 2013

- Article: Nature 494:234, [doi:10.1038/nature11867](https://doi.org/10.1038/nature11867).
- Data: 1008 segregants, 30,594 SNPs (11,623 unique markers), 46 traits, replicate measurements. [Princeton data page](http://genomics-pubs.princeton.edu/YeastCross_BYxRM/data.shtml), online on 2026-10-31, no licence stated: `BYxRM_GenoData.txt` (23.7 MB), `BYxRM_PhenoData.txt` (0.8 MB), `pheno_raw.Rdata` (0.9 MB), `cross.Rdata` (22.2 MB).
- Printed values: Supp Table 1 ([MOESM86 xls](https://static-content.springer.com/esm/art%3A10.1038%2Fnature11867/MediaObjects/41586_2013_BFnature11867_MOESM86_ESM.xls)) gives per-trait H² and h² with SE (Cadmium_Chloride H² 0.960, h² 0.849). h² was fitted with rrBLUP, H² with lme4. Supp Table 2 (MOESM87) splits additive variance over 16 per-chromosome kernels, with some estimates pinned at 1e-9, which makes a 16-kernel constrained REML test.
- A + A∘A on these data: Young & Durbin 2014, Genetics 198:1405, [doi:10.1534/genetics.114.170795](https://doi.org/10.1534/genetics.114.170795), GCTA AI-REML, mean pairwise share 0.10. The per-trait Table S1 / File S5 was not opened.
- Bloom et al. 2019 eLife 8:e49212 ([article](https://elifesciences.org/articles/49212)) extends this to 16 crosses, 13,950 segregants and 38 traits with A, A∘A and repeatability via the R package regress. The per-trait supplementary file was not opened; data links are in [joshsbloom/yeast-16-parents](https://github.com/joshsbloom/yeast-16-parents).

## 6. Twin data: OpenMx `twinData` and mets `twinbmi`

All printed twin estimates found are ML. No printed REML twin estimate on open data was found. A twin ACE model is a linear mixed model with A = 1 within MZ pairs and 0.5 within DZ pairs, C a pair effect, D = 1 (MZ) and 0.25 (DZ). remlax must run in ML mode for these comparisons.

- `twinData` ([doc](https://rdrr.io/cran/OpenMx/man/twinData.html); OpenMx on [CRAN](https://cran.r-project.org/package=OpenMx), Apache-2.0; [file](https://raw.githubusercontent.com/OpenMx/OpenMx/master/data/twinData.rda), 44 KB). 3808 Australian twin pairs, wide format, with wt, ht, htwt, bmi (7·log BMI), age and cohort. Pair counts computed locally: younger cohort MZFF 569, MZMM 273, DZFF 351, DZMM 206, DZOS 510.
- Printed ML reference, hard-coded Mx values in [demo/UnivariateTwinAnalysis_MatrixRaw.R](https://github.com/OpenMx/OpenMx/blob/master/demo/UnivariateTwinAnalysis_MatrixRaw.R): bmi, MZ = zyg 1, DZ = zyg 3, one common mean, FIML with incomplete pairs. VA 0.6173023, VC 5.6e-14 (boundary), VE 0.1730462, mean 21.39293, -2lnL 4067.663. Reproduced locally.
- mets `twinbmi` (mets on [CRAN](https://cran.r-project.org/package=mets), Apache-2.0; [quantitative-twin vignette](https://cran.r-project.org/web/packages/mets/vignettes/quantitative-twin.html)): 11,188 Danish records (doc says 11,411 subjects, which does not match the data). Printed ACE with `bmi ~ age + gender`: sd(A) 2.7270, sd(C) 1.7e-6, sd(E) 2.0276, intercept 18.599, age 0.11892, male 1.3846, logLik -22019.66 (df 6). Reproduced locally. ADE is not printed; computed locally logLik -22017.51.
- umx ([refman](https://cran.r-project.org/web/packages/umx/refman/umx.html)) ships GFF, iqdat, docData and us_skinfold_data, without printed variance estimates (not verified).
- Rabe-Hesketh, Skrondal & Gjessing 2008 (Biometrics 64:280, [doi:10.1111/j.1541-0420.2007.00803.x](https://doi.org/10.1111/j.1541-0420.2007.00803.x)) fit twin and family models as mixed models with printed ML numbers ([Stata 2010 slides](https://www.stata.com/meeting/germany10/germany10_rabe-hesketh.pdf)). The data files could not be obtained (Wiley supplement 403). Unusable until found.
- Refit with: OpenMx, umx, mets, asreml or sommer with the pair-level A, C and D matrices, gaston multi-kernel.

## 7. solarius `dat30` and coxme GAW pedigrees

- solarius (archived from CRAN on 2020-03-07, GPL ≥ 3, [tarball](https://cran.r-project.org/src/contrib/Archive/solarius/solarius_0.3.0.2.tar.gz) 486 KB, [repo](https://github.com/ugcd/solarius)). `dat30`: 174 individuals in 29 nuclear families, two simulated, highly correlated traits, age and sex. Printed SOLAR ML output in the [tutorial](http://ugcd.github.io/solarius/vignettes/tutorial.html): `trait1 ~ age + sex` h2r 0.8061621 (SE 0.1100465), loglik -210.8689; bivariate `trait1 + trait2 ~ 1` h2r 0.8218823 and 0.6270026, rhog 0.9728759, rhoe 0.4120487. SOLAR omits the (n/2)·log(2π) constant (159.895 here). The univariate h2 was reproduced locally. `dat50` failed in SOLAR and is useless.
- coxme ([CRAN](https://cran.r-project.org/package=coxme), LGPL-2), vignette "The lmekin function", data `vignettes/gaw.rda` in the source tarball. 1497 individuals in 23 pedigrees, plus a SOLAR IBD matrix at one locus. The GAW year is not stated by the author. Printed ML fits with `2·kinship`: gfit1 `age ~ q1 + (1|id)` logLik -4114.317; gfit2 kinship + IBD (two kernels) logLik -4252.44; gfit3 kinship + `(1|famid)` logLik -4114.079. gfit1 logLik reproduced locally. The likelihood surface is flat, so compare logLik rather than variances.
- Refit with: SOLAR, coxme::lmekin, gaston multi-kernel, asreml, sommer.

## 8. Drosophila DGRP

- The official site dgrp2.gnets.ncsu.edu did not answer on 2026-10-31; a [Wayback snapshot](https://web.archive.org/web/20250125054114/http://dgrp2.gnets.ncsu.edu/data.html) lists VCF, plink files (bed 230.8 MB), `.tgeno`, a FastLMM GRM and line-mean CSVs. Genotype mirror: [Zenodo 10.5281/zenodo.155396](https://doi.org/10.5281/zenodo.155396), CC BY 4.0, `dgrp2.vcf.gz` 193.8 MB, 4,438,427 variants. Freeze 2: 205 lines (Huang et al. 2014 Genome Res, [PMC4079974](https://pmc.ncbi.nlm.nih.gov/articles/PMC4079974/)).
- Line means for 838 phenotypes: [DGRPool](https://dgrpool.epfl.ch) (Gardeux et al. 2024 eLife, [doi:10.7554/eLife.88981](https://doi.org/10.7554/eLife.88981)), no licence stated.
- Ober et al. 2015 PLoS ONE 10:e0126880 ([doi:10.1371/journal.pone.0126880](https://doi.org/10.1371/journal.pone.0126880), CC BY), chill coma recovery: S1 Table holds 35,479 individual records for 176 lines. S2 Table (SAS Proc Mixed, pooled sexes): line 27.33 (3.17), sex×line 2.00 (0.51), rep(sex×line) 4.02 (0.37), error 49.92 (0.38). S4 Table, females: G 14.36 (15.24), iid line 18.50 (14.91), rep 3.03, error 55.38. Supports a genomic + iid line + replicate model on raw records.
- Ober et al. 2012 PLoS Genet 8:e1002685 ([doi:10.1371/journal.pgen.1002685](https://doi.org/10.1371/journal.pgen.1002685)): ASReml 2.0 with line, sex×line, replicate, G and G∘G on individual records (Tables S3, S4, with ln L; startle Model 3 ln L -27317.41). The individual records are not verified as public; public versions are line means. Exact refit not possible on what was found.
- Dominance is not identifiable in inbred lines.

## 9. Animal pedigree data with maternal, rearing and repeated-measure terms

- MCMCglmm `BTdata` / `BTped` ([CRAN](https://cran.r-project.org/package=MCMCglmm), GPL-2; Hadfield et al. 2007 J Evol Biol 20:549, [doi:10.1111/j.1420-9101.2006.01262.x](https://doi.org/10.1111/j.1420-9101.2006.01262.x)). 828 blue tit chicks from 106 full-sib families, cross-fostered; tarsus, back colour, hatch date; genetic dam and foster nest; pedigree of 1040 with 212 unrelated founders. Copy in enhancer 1.1.2 as `DT_btdata`, whose help prints a sommer REML fit of `tarsus ~ sex` with random `dam + fosternest`: dam 0.22021839 (SE 0.03228929), fosternest 0.06924283 (0.01962657), units 0.56789911 (0.02163284). The help also gives a 3-trait `usm` model without printed output. Supports animal + dam + foster nest, multi-trait.
- pedigreemm `milk` ([CRAN](https://cran.r-project.org/package=pedigreemm), GPL ≥ 2; Vazquez et al. 2010 J Anim Sci 88:497, [doi:10.2527/jas.2009-1952](https://doi.org/10.2527/jas.2009-1952), paywalled). 3397 lactations of 1359 cows, 38 sires, 57 herds, lactations 1 to 5 (1314, 1006, 640, 334, 103); milk, fat, protein, SCS; pedigree `pedCows` (6547). Repeated records support A + permanent environment and multi-trait. The help states kg with range 4065 to 19,345, but the data range is 8953 to 42,614, consistent with pounds. No printed variance components were verified. `mastitis` (1675 cows, count data) suits GLMMs only.
- Gryphon (Wilson et al. 2010) is already covered. enhancer `DT_legendre` is simulated random-regression data.

## 10. Other usable data

- Atlantic salmon AGD, Robledo et al. 2018 G3 8:1195 ([doi:10.1534/g3.118.200075](https://doi.org/10.1534/g3.118.200075), CC BY 4.0). [Europe PMC bundle](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5873910/supplementaryFiles), 11 MB: 1480 fish with gill score and amoebic load, sire/dam, 17,156 SNPs with positions. ASReml 3 h² (Table 2): gill score 0.25 (pedigree) / 0.24 (G), amoebic load 0.36 / 0.25. Bivariate genetic correlation "∼1", a boundary case.
- DO mice, Gatti et al. 2014 G3 4:1623 ([doi:10.1534/g3.114.013748](https://doi.org/10.1534/g3.114.013748)), data in [qtl2data](https://github.com/rqtl/qtl2data) `DO_Gatti2014/do.zip` (1.6 MB; the repository has no licence file). 742 mice, 7854 SNPs, WBC and NEUT. Printed neutrophil h² 0.132 (R package regress). Kinship comes from founder probabilities, so the refit depends on haplotype reconstruction.
- Gough × WSB F2 in qtl2data: 1212 mice, weekly body weight over 16 weeks, parents recorded. Suits multi-trait and random regression; no printed variance components found.
- qtl2 user guide ([vignette](https://kbroman.org/qtl2/assets/vignettes/user_guide.html)) shows `calc_kinship(pr, "loco")` and LMM scans but prints no heritability or LMM LOD values.
- GCTA test files (`test.bed/.bim/.fam/.phen`) inside the [GCTA zip](https://yanglab.westlake.edu.cn/software/gcta/) (99 MB, mostly the binary, MIT): 3925 individuals, 1000 SNPs on chr1. The GCTA documentation prints example `.hsq` outputs (GREML, `--mgrm`, bivariate) on data that are not distributed, and nothing for GREML-d. Smoke test only.
- German Holstein, Zhang et al. 2015 G3 5:615 ([doi:10.1534/g3.114.016261](https://doi.org/10.1534/g3.114.016261)): 5024 bulls, 42,551 SNPs, phenotypes are standardised EBVs, no pedigree. Weak for variance components.
- QTLMAS 2010 (simulated, epistatic QTL; [PMC3103202](https://pmc.ncbi.nlm.nih.gov/articles/PMC3103202/)): the download host did not respond; not verified.

## 11. Open human genotypes without phenotypes

- 1000 Genomes phase 3: 2504 samples, no embargo ([IGSR disclaimer](https://www.internationalgenome.org/IGSR_disclaimer)); files and pedigree at [EBI FTP](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/); PLINK 2 versions at [cog-genomics](https://www.cog-genomics.org/plink/2.0/resources). Genotypes only, usable for simulated phenotypes including a dominance GRM.
- Geuvadis E-GEUV-1 (Lappalainen et al. 2013 Nature, [doi:10.1038/nature12531](https://doi.org/10.1038/nature12531)): RNA-seq on 465 1000G LCLs, open at [EBI](https://ftp.ebi.ac.uk/biostudies/fire/E-GEUV-/001/E-GEUV-1/Files/E-GEUV-1/analysis_results/) (expression matrix 87 MB, genotypes 0.6 to 3.9 GB per chromosome). Real phenotypes on open genotypes, but no printed per-gene GREML table was found.
- HapMap: retired; files remain on [NCBI FTP](https://ftp.ncbi.nlm.nih.gov/hapmap/) (hapmap3 PLINK ped 1.2 GB). No phenotypes.
- OpenSNP: shut down with deletion of all data on 2025-04-30 ([announcement](https://tzovar.as/sunsetting-opensnp/)). Unusable.
- Personal Genome Project: CC0 waiver encouraged ([PGP](https://www.personalgenomes.org/us)); portal unreachable on 2026-10-31, status not verified.

## 12. Restricted or unavailable (unusable for public validation)

| resource | access | source |
|---|---|---|
| dbGaP studies (Framingham phs000007, GAW16, GAW19) | eRA Commons account and DAC approval | [dbGaP](https://dbgap.ncbi.nlm.nih.gov/) |
| UK Biobank | application and fee | [UKB access](https://www.ukbiobank.ac.uk/enable-your-research/apply-for-access) |
| ALSPAC | proposal and fee | [ALSPAC](https://www.bristol.ac.uk/alspac/researchers/access/) |
| TwinsUK | application with institutional email, fee | [TwinsUK](https://twinsuk.ac.uk/researchers/access-data-and-samples/request-access/) |
| Vitezica et al. 2018 Genus pigs (litter size, D and AA) | on request | [doi:10.1186/s12711-018-0437-3](https://doi.org/10.1186/s12711-018-0437-3) |
| Joshi et al. 2020 tilapia genotypes | data transfer agreement | [doi:10.1186/s12711-019-0522-2](https://doi.org/10.1186/s12711-019-0522-2) |
| Aliloo et al. 2016 Holstein, Bolormaa et al. 2015 beef cattle | no data statement found | [PMC4736671](https://pmc.ncbi.nlm.nih.gov/articles/PMC4736671/), [PMC4382858](https://pmc.ncbi.nlm.nih.gov/articles/PMC4382858/) |
| HS rats genotypes (UCSD library) | page blocked automated access; open status not verified | [UCSD](https://library.ucsd.edu/dc/object/bb29129987) |

No public salmon, trout, shrimp or oyster data set with a published genomic dominance or epistatic estimate was found.

## Suggested order of adoption

The tilapia data come first for pedigree D and maternal effects, because ASReml-4 printed every model on files under 250 KB. Bloom 2015 follows for A + A∘A + repeatability on raw replicates. The PIC data then test genomic D and A×A at n = 3534 and 52,843 SNPs, with a 934 MB genotype file to convert once. `warcolak` adds a pedigree D test against simulated truth with nadiv + asreml and gremlin. The twin, solarius and coxme references test ML mode, boundary estimates (VC = 0) and multiple kernels on small files.
