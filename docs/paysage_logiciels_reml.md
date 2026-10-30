# The landscape of REML solvers and the position of remlax

Survey made on 26 October 2026 from primary sources: official documentation, manuals, source code, and articles with a verified DOI. Appendices A, B and C describe each program and each mechanism, with numbered sources. This summary keeps only the findings that bear on the decision. The remlax measurements cited here are internal and unpublished (mvD7: n = 12 760, 51 terms, 198 parameters).

Question asked: can remlax replace ASReml for heavy computations, or become a free alternative of the sommer kind, but more efficient?

## 1. Changes in 2026

The niche "free REML, asreml-style structures, GPU" was empty a year ago. It is no longer empty.

- pyreml (Marchal, Trottier, Bry and Raimondi, IGMM and IMAG, Montpellier) appeared on PyPI in June 2026. Its bioRxiv preprint is dated 24 October 2026 (doi:10.64898/2026.10.20.757556). It performs exact REML on GPU in PyTorch, with a dense V or a dense Woodbury form, and an analytical gradient. It has no sparse matrices, no AI, and no standard errors for the variances. Its benchmark compares it with WOMBAT, not with ASReml (appendix B).
- kmmr (Melsen, van Eeuwijk and Peeters, Wageningen) appeared on arXiv on 12 October 2026 (arXiv:2610.16671). It performs multivariate AI-REML on the mixed model equations. The sparse part uses CHOLMOD on CPU with Takahashi inversion, and the dense Schur complement runs on GPU. On a wheat GBLUP with 10 375 genotypes, it reports 49 s on an A100 against 989 s for ASReml-R and 512 s for sommer. It cites the V form as future work (appendix B).
- sommer 4.5.12 (CRAN, 25 October 2026) now solves the sparse mixed model equations by default. It uses CHOLMOD, Takahashi inversion, optional METIS ordering, and optional CUDA kernels for the dense part. The direct formulation on V remains available but is heading towards deprecation. The package changes very fast: version 4.4.90 broke the control interface (appendix A).

All three projects are free software (GPL-3). Two of them explicitly target plant breeding.

## 2. Direct comparison

| | remlax | ASReml-R 4.2 | sommer 4.5.12 | kmmr | pyreml 0.5.3 |
|---|---|---|---|---|---|
| Licence | GPL ≥ 3 | commercial, from £375 (academic) | GPL-3 | GPL-3, experimental | GPL-3 |
| Computational forms | dense V, sparse MME, block V (Woodbury) | sparse MME, dense/sparse partition | sparse MME; dense V being deprecated | sparse MME + dense Schur | dense V or dense Woodbury |
| Optimiser | L-BFGS-B, AI-REML, trust-region finisher | AI-REML | AI-REML, EM/PX-EM fallback | AI-REML | L-BFGS (PyTorch), Adam |
| GPU | yes (JAX): V, blocks, dense MME | no | optional CUDA, dense kernels only | Schur complement | yes (PyTorch) |
| Sparse ordering | natural order (METIS measured, not wired in) | undocumented | AMD, optional METIS | CHOLMOD | not applicable |
| asreml-style structures | us, fa, rr, ar1⊗ar1, separable, K or K⁻¹ | complete | us, fa, rr, arbitrary Kronecker products, K or K⁻¹ | us, extended FA, AR1×AR1, kernels | iid, diag, full, fa, ar, distances |
| Standard errors of variances | yes (AI, Hessian) | yes | yes | yes (AI) | no |
| Automatic structure detection | yes (global terms and blocks) | no: partition set by the user | no | no | switches between V and Woodbury depending on q < n |
| Kernel IGE neighbourhood | yes, estimated range | by hand | by hand | no | no |
| Published speed | none | reference | 2 times slower than kmmr CPU on 10 375 genotypes | 20 times ASReml (GBLUP, A100) | more than 12 times WOMBAT (A100) |

The remlax cells rest on our internal measurements. The other columns come from appendices A and B.

## 3. Features unique to remlax

According to the sources consulted, no other tool combines the following features.

- The V and MME forms in a single solver, with the same likelihood and AI-REML. pyreml has no MME form. kmmr has no V form.
- Automatic detection of global terms and blocks, followed by V = A + UU' handled by Woodbury with the blocks batched (vmap). Appendix C finds no prior art. On mvD7, one evaluation takes 0.90 s on a 10 GB MIG slice, against 4.85 s for our low-memory mode.
- A heavy multivariate model fitted on a 10 GB MIG slice. No tool reports this use.
- The active set on the positive semi-definite cone on the variance scale, and the Moré and Sorensen finisher with an SR1 update. Appendix C finds no use of these in REML software.
- REML estimation of the continuous range of an IGE kernel together with the variances. WOMBAT only profiles the dilution. Sato et al. (2021) sweep the scale.

Several building blocks, on the other hand, have clear prior art (appendix C). They include the gradient P − Pyy'P and the dense V form (Harville 1977), and differentiation of the Cholesky factor (Smith 1995, Meyer and Smith 1996). They also include the supernodal Cholesky with selected inversion (YAMS, WOMBAT, ASReml) and nested dissection on C (WOMBAT). The likelihood taken from a bordered matrix (Smith and Graser 1986, Graser et al. 1987), AI-REML (Gilmour et al. 1995) and its damping (WOMBAT, SAS) complete the list. An article will have to cite them as foundations and claim no contribution from them.

## 4. Gaps in remlax

- Computation on several GPUs or several nodes. DISSECT, REACTA, SAIGE-GPU, ExaGeoStat and DALIA do it. In JAX, jaxmg provides a distributed Cholesky with a log-determinant (appendix B).
- The matrix-free route, with preconditioned conjugate gradients and stochastic traces. BOLT-REML, SAIGE, SLEMM, the Monte Carlo REML of BLUPF90 and GPyTorch use it. Only this route goes beyond a few hundred thousand individuals.
- The mixed FP32 then FP64 precision of pyreml. It helps on GPUs without fast FP64 arithmetic (L40S).
- Genotypes compressed to 2 bits and streamed to the GPU (miraculix, Quickdraws).
- METIS ordering in our MME engine. sommer and WOMBAT already have it. Our measurement on mvD7 shows a factor of C four times smaller (0.69 GB against 2.72 GB in natural order).
- Published benchmarks against current competitors.

## 5. Answer to the question

Replacing ASReml for heavy computations is a realistic goal. ASReml has no GPU support, is capped at 96 GB and 16 threads, and its dense/sparse partition is manual. On mvD7, one AI iteration of asreml takes 50 min with the default partition and 31 min with the genetic terms set to dense. The remlax block engine evaluates the same likelihood in 0.9 s. Yet kmmr already claims a factor of 20 over ASReml, and the comparison will have to be made against it.

Becoming "a more efficient sommer" has become harder. sommer 4.5 has adopted the same sparse techniques as remlax, together with partial GPU support. The difference will rest on measurements, on interface stability, and on features sommer lacks. These are the V form and the blocks on GPU, automatic detection, and neighbourhood models. The list of functions will no longer set the two apart.

## 6. Recommendations

1. Direct benchmarks, on the same models and at the same likelihood fixed point, against sommer 4.5.12 (mmes with METIS, with and without CUDA), kmmr, pyreml and ASReml-R (default and dense partition). Start with the families of the engine benchmark and with mvD7.
2. Wire METIS (or the "blocks then separator" order of decoupage.py) into the MME engine. The gain has been measured, and the competitors have it.
3. Stabilise the R interface before any release. Interface instability is a documented weakness of sommer.
4. Aim for multi-GPU or the matrix-free route only if the target exceeds 10⁵ genotyped individuals. For the plant breeding designs in view, the block form seems sufficient.
5. For an article, centre the contribution on the union of the V and MME forms, automatic structure detection and neighbourhood models. Cite pyreml and kmmr as contemporary work. Claim no priority on GPU REML.
6. Contact the pyreml team, which is in Montpellier like the user. Their choices (dense Woodbury, mixed precision, index grids) overlap with ours.

## 7. Unverified points

- The WOMBAT, DMU and MTG2 websites were unreachable. Their entries rest on old manuals or on excerpts (appendix A).
- The exact scale of the kmmr benchmarks (number of observations not given) needs checking in the full text. So does the log-likelihood gap against ASReml-R on an FA1 model (−3 948.66 against −3 937.65).
- The WCGALP 2022 proceedings were blocked. Possible GPU work in animal breeding has not been read.
- ASReml is closed source. Undocumented heuristics may exist.
- Priority for "the first general REML on GPU" remains disputed between pyreml and kmmr. REACTA (2014) had done it for a specialised case.

The appendices follow: A, reference software and alternatives; B, accelerated computing; C, prior art of the mechanisms.

## Appendix A. Reference software and alternatives for REML in quantitative genetics

Survey made on 26 October 2026. The sources are official pages, manuals, code repositories and articles with a DOI verified on CrossRef or PubMed. Numbers in square brackets refer to the list of sources at the end of the document. The mention "unverified" flags a point I could not confirm on a primary source.

Three official servers were unreachable during the survey: the WOMBAT site (didgeridoo.une.edu.au), the DMU site (dmu.ghpc.au.dk) and the MTG2 page (Google Sites, login required). The information on these three programs therefore rests on older manuals or on search engine excerpts, flagged as such.

### 1. Comparison table

Abbreviations. MME: Henderson's mixed model equations. Dense V: formulation on the variance matrix of the observations, of size n × n. AI: average information algorithm. K: relationship matrix supplied as is. K⁻¹: its inverse, supplied. Closed: binary without source code.

| Software | Licence | Language | Algorithm | Sparse or dense | Multivariate us / fa | ar1⊗ar1 | Relationship matrix | GPU | Latest version |
|---|---|---|---|---|---|---|---|---|---|
| ASReml-R 4.2 | Commercial, closed; academic licence from £375 [1] | Fortran engine (unverified), R interface | AI-REML [2][3] | Sparse MME, dense/sparse partition [2] | Yes: us, fa, facv, xfa, rr [4] | Yes [4] | K or K⁻¹ (ainverse, sparse matrices) [2] | No | 4.2.0.480-2, 6 May 2026 [5] |
| ASReml-SA 4.2 | Commercial, same licence as ASReml-R [4] | Unverified | AI-REML [4] | Sparse MME; workspace up to 96 GB; 16 threads by default [4] | Yes [4] | Yes [4] | K or K⁻¹ [4] | No | Release 4.2 (October 2022 specification); current build unverified [4] |
| Echidna | Free except for commercial use, closed, "alpha" version [6] | Unverified | AI-REML [6] | Unverified | Unverified | Unverified | H⁻¹ and GRM (release note) [7] | Unverified | Build 2.15, 22 October 2026 [7] |
| WOMBAT | Free for the scientific community, executable only [8] | Unverified | PX-EM then AI; EM; derivative-free Powell and Simplex [8] | Supernodal sparse MME; MMD, AMD, METIS ordering; --dense option [8] | us; reduced rank by partial Cholesky; random regression [8] | Not documented in the 2017 manual | K⁻¹ through a .gin file [8] | No | Manual of 31 May 2017; current version unverified (site unreachable) |
| BLUPF90+ (REMLF90, AIREMLF90), GIBBSF90+ | Free for research; free version limited to 30 000 genotyped animals; commercial use by agreement; binaries, old "educational" sources [9][10][11] | Fortran 90/95 [12] | EM-REML (remlf90), AI-REML (airemlf90), Gibbs [11] | Sparse MME; FSPAK then supernodal YAMS [13] | Multi-trait yes [11]; fa unverified | Not native; possible through breedR | A⁻¹ from pedigree, single-step H⁻¹ [11][14] | No | Programs page updated on 5 June 2026 [11] |
| DMU | Free for research, executables; commercial use by agreement (search excerpt) [15] | Unverified | AI-REML of Jensen et al. 1997, EM fallback (DMUAI) (excerpt) [15] | Sparse MME, FSPAK [16] | Multivariate (excerpt) [15] | Unverified | Unverified | No | Guide v6 release 5.6, June 2024 (excerpt) [15] |
| MiX99 | Commercial; Standard, Core and University licences; animal breeding only [17] | Unverified | MC-EM-REML; AI-REML cited in 2024 [17] | Iterative solution; MPI and multithreading [17] | Multi-trait [17] | No | Unverified | No | Release IX/2025 (excerpt) [17] |
| MiXBLUP 3 | Commercial: €2 400 to €24 790 per year; free one-month trial [18] | Unverified | Not stated; core derived from MiX99 [19] | Unverified | Unverified | No | Unverified | No | MiXBLUP 3 (2026 prices) [18] |
| MTG2 | Free; v2.15 sources deposited on Dryad; software licence not stated [20] | Fortran [20] | Direct AI-REML + eigendecomposition of the GRM [21] | Dense V [21] | Multi-trait, random regression, several GRMs [21] | No | Dense K [21] | Unverified | v2.16, 17 August 2020 [20] |
| GCTA (GREML) | GPL-3, parts under LGPL [22] | C++ [22] | Direct AI-REML, HE [23][24] | Dense V; fastGWA-REML on a sparse GRM [25] | Bivariate only [16] | No | K (GRM) | No | v1.94.1 (2022 tag); last commit July 2025 [22] |
| LDAK | GPL-3 [26] | C (unverified) | REML on multiple kernels [27]; algorithm details unverified | Dense V | Unverified | No | K | No | ldak7.01; repository active in October 2026 [26] |
| sommer | GPL-3 [28] | R + C++ (Armadillo, Eigen, CHOLMOD) [28] | AI-REML on sparse MME (mmes, default); direct NR/AI inversion (mmer, direct formulation); EM/PX-EM recovery [28] | Both; CHOLMOD, Takahashi selected inverse, optional METIS, stochastic PCG [28] | us, diag, fa, rr [28] | Yes, vsm(ar1m(), ar1m()) [28] | K or K⁻¹ (inverse attribute) [28] | Optional CUDA, dense kernels only [28] | 4.5.12, 25 October 2026 [28] |
| lme4 | GPL ≥ 2 [29] | R + C++ (Eigen) [29] | PLS + derivative-free optimiser [30] | Sparse Cholesky [30] | Univariate; us, diag, cs, ar1 per term since 2.0 [31] | No | Not native | No | 2.1-0, 23 October 2026 [29] |
| lme4qtl | GPL-3 [32] | R | lme4 with Z·L, L a factor of K [33] | Sparse if K is sparse [33] | Through lme4 | No | K [33] | No | Last commit January 2021 [32] |
| pedigreemm | GPL ≥ 2 [34] | R | lme4 with a factor of A [35] | Sparse | No | No | Pedigree only [36] | No | 0.3-5, 14 October 2024 [34] |
| lme4breeding | GPL ≥ 2 [37] | R | lme4 with Z·L; rotation by eigendecomposition [38] | Sparse | us documented in the preprint [38] | Unverified | K [38] | No | 1.1.5, 5 October 2026 [37] |
| breedR | GPL-3 for the R interface; BLUPF90 binaries [39] | R + Fortran binaries [39] | REMLF90 (EM), AIREMLF90 (AI), Gibbs [39] | Sparse MME | Unverified | Yes, AR1⊗AR1 and splines [39] | Pedigree [39] | No | 0.12-7 on GitHub; no longer maintained [39] |
| glmmTMB | AGPL-3 [40] | R + C++ (TMB) [40] | ML or REML by Laplace, gradient by automatic differentiation [40] | Sparse (TMB) | us, rr, diag, cs, toep, etc. [41] | ar1 per term; Kronecker product undocumented | Dense K through propto or equalto [41] | No | 1.1.15.2, 29 September 2026 [40] |
| spaMM | CeCILL-2 [42] | R + C++ [42] | ML or REML by Laplace, h-likelihood [42] | Sparse for Markov fields and AR [42] | Multivariate (fitmv) [42]; us unverified | AR1, Matérn; Kronecker unverified | K (corrMatrix) or precision [42] | No | 4.7.17, 9 October 2026 [42] |
| MCMCglmm | GPL ≥ 2 [43] | R + C++ | Bayesian MCMC [44] | Sparse A⁻¹ [45] | us, idh, multivariate [44] | No | K⁻¹ through ginverse [45] | No | 2.36, 6 May 2024 [43] |
| brms | GPL-2 [46] | R + Stan | HMC (NUTS) [47] | Dense | Multivariate [46] | Not native (unverified) | K through gr(cov = ) | OpenCL through cmdstanr [48] | 2.23.0, 9 September 2025 [46] |
| R-INLA | MIT (repository) [49] | C (GMRFLib) + R [49] | INLA, no REML [50] | Sparse precision [50] | Unverified | ar1 and product by group [51] | K⁻¹ through generic0 [51] | Unverified | 26.10.26, 26 October 2026 [49] |
| MixedModels.jl | MIT [52] | Julia [52] | PLS + derivative-free optimiser (NEWUOA by default in v5); REML as an option [52] | Block sparse | No | No | Not native | Undocumented | v5.9.1, 2 October 2026 [52] |
| statsmodels MixedLM | BSD-3 [53] | Python [53] | ML or REML after Lindstrom and Bates [54][55] | Dense per group [54] | No | No | No | No | 0.15.0, 27 August 2026 [53] |

### 2. Software entries

#### 2.1 ASReml (ASReml-R 4.2 and ASReml-SA 4.2)

ASReml is commercial software from VSN International. The academic licence starts at £375. Commercial prices are on request [1]. ASReml-R runs under the same licence as the stand-alone program ASReml-SA [4]. The licensing system relies on Reprise, with cloud management [4]. I found no free licence for low-income countries on the official page.

The engine combines the AI algorithm of Gilmour, Thompson and Cullis (1995) with sparse matrix methods [2][3]. The ASReml-R manual states that only non-zero elements are stored [2]. The terms are split between a dense set and a sparse set [2]. The inverse of the coefficient matrix is formed in full for the dense set, and only in part for the sparse set [2]. Random effects are sparse by default. A dense genomic relationship matrix can be moved into the dense set with the dense option [2]. The sparse equations are reordered to preserve sparsity [2].

The scale limits are documented. Version 4.2 raises accessible memory to 96 GB and parallelism to 16 threads [2][4]. The default workspace of ASReml-SA is 2 GB [4]. The ASReml-R manual states that some very large problems require ASReml-SA because of the overheads of R [2]. Version 4.2.0.393 of 29 August 2025 raised the predict workspace to 32 GB [5]. The variance structures cover us, fa, facv, xfa, rr, ar1 and their direct products [4]. The latest version of ASReml-R is 4.2.0.480-2 of 6 May 2026, built for R 4.6 [5]. No GPU is mentioned.

#### 2.2 Echidna

Echidna is developed by Arthur Gilmour, lead author of ASReml [6]. It is free except for commercial use. Commercial use requires a current ASReml licence [6]. The download requires registration and acceptance of a licence [6]. The source code is not distributed. The site presents the program as an alpha version updated at least monthly [6]. It runs on Mac, Linux and 64-bit Windows [6].

Variance parameter estimation relies on the AI algorithm [6]. The note for build 2.15 (22 October 2026) reports forming an H⁻¹ matrix for 201 150 pedigree individuals and 40 000 GRM individuals in 10 minutes instead of 9 hours. Each iteration with this matrix takes 24 minutes [7]. I could not read the manual or the WCGALP 2018 presentation article (dead link). The available structures and the language therefore remain unverified.

#### 2.3 WOMBAT

WOMBAT is developed by Karin Meyer (AGBU, University of New England) [8][56]. It is free for the scientific community, with an obligation to cite it [8]. Only the executable is distributed, with the manual and examples, for Linux and Windows [8][56].

By default, WOMBAT runs a few PX-EM iterations and then switches to AI-REML [8]. It also offers classical EM and two derivative-free maximisations, Powell and Simplex [8]. Computations are sparse by default, with a supernodal approach that calls BLAS and LAPACK on the dense sub-blocks [8]. The ordering can be MMD, AMD or METIS [8]. A --dense option exists for dense GRMs. Compact storage then limits the system to 65 535 equations. The --dense2 option removes this limit [8]. Matrices are unstructured by default, with optional reduced-rank estimation through a partial Cholesky factor [8]. A known covariance matrix is entered in inverse form, in a .gin file [8]. The 2017 manual documents no spatial autoregressive structure. A 2020 abstract mentions penalised REML, GRMs, association scans and single-step [57]. The current version could not be verified.

#### 2.4 The BLUPF90 family

The programs are developed at the University of Georgia by I. Misztal, D. Lourenco and co-workers [12]. They are free for research, with mandatory citation. Commercial use requires contacting the authors [9]. The research version limits the number of genotyped animals to 30 000 [11]. Binaries are supplied for Linux, Windows and Mac. Only old versions of the source code are available "for educational purposes" [10]. The current programs are therefore not free software. Some programs (BLUP90IOD3, ACCF90GS2) are distributed only under agreement [11].

BLUPF90+ groups blupf90, remlf90 (EM-REML) and airemlf90 (AI-REML). GIBBSF90+ groups the Gibbs samplers [11]. Moving from FSPAK to the supernodal package YAMS sped up AI-REML by an average factor above 10, over 23 models [13]. Trace computation gained up to a factor of 67 thanks to hash table storage [13]. Single-step (H⁻¹) is the speciality of the family [14]. The wiki reports a gradual move from the ifort compiler to ifx in 2025-2026 [10]. AR1⊗AR1 spatial structures are not native. breedR builds them in its R layer before calling the binaries [39].

#### 2.5 DMU

DMU is developed by Per Madsen and Just Jensen (QGG, Aarhus University). The official server refused connections during the survey. According to search engine excerpts of the official guide, the package is distributed as executables, free for research, and commercial use requires an agreement with QGG [15]. The DMUAI module estimates the components by AI-REML after Jensen et al. (1997), with a possible EM fallback [15]. The most recent guide found is version 6, release 5.6, of June 2024 [15]. The HIBLUP article states that DMU uses FSPAK and that the distributed version offers no parallel computation [16]. These points remain to be confirmed in the guide itself.

#### 2.6 MiX99 and MiXBLUP

MiX99 is software from Luke (Natural Resources Institute Finland) [17]. Three licences exist: Standard, Core for emerging countries, and University for academic use [17]. The licences are reserved for animal breeding applications [17]. The mix99s solver is multithreaded. mix99p and apax99p are MPI programs [17]. The cited references deal with preconditioned conjugate gradients and with a Monte Carlo EM REML. A 2024 reference deals with AI-REML [17]. The page gives no price.

MiXBLUP is developed by Wageningen University & Research together with Luke. Its core derives from MiX99 [19]. Licences cost from €2 400 per year (pedigree, fewer than one million equations) to €24 790 per year ("Large EC" licence, 3 to 10 workstations) [18]. A one-month trial is free. No academic licence is mentioned [18]. These pages do not state the REML method.

#### 2.7 MTG2

MTG2 is developed by S. Hong Lee [21]. It combines the "direct" AI-REML of Lee and van der Werf (2006) with an eigendecomposition of the GRM [21][58]. The direct algorithm works on the matrix V of the observations and not on the MME [21]. It pays off when the genetic covariance is dense [21]. The eigendecomposition requires a single full-rank GRM and a balanced design, with no missing phenotype across traits [21]. The Dryad deposit of 17 August 2020 contains version 2.16, the Fortran sources of 2.15 and the manual [20]. The software licence is not stated.

#### 2.8 GCTA and LDAK

GCTA is under GPL-3, with parts under LGPL [22]. The latest tag in the repository is v1.94.1 (2022), and the last commit dates from July 2025 [22]. GREML relies on direct AI-REML on a dense V [23][24]. The HIBLUP article states that GCTA handles only two traits at a time and has no function for random environmental effects [16]. fastGWA-REML exploits a sparse GRM for large human cohorts [25]. LDAK is under GPL-3 (version 3 or later) [26]. It estimates heritability by REML on weighted kernels [27]. Both programs target human genetics. They offer neither general us/fa structures nor spatial structures.

#### 2.9 sommer

sommer is the free R package of G. Covarrubias-Pazaran, under GPL-3 [28]. Version 4.5.12 appeared on CRAN on 25 October 2026 [28]. The package has changed a great deal since 2025.

The history of the functions can be read in the CHANGELOG [28]. mmer is the historical function, with direct inversion on the r × r matrix of the observations. mmec, added in 4.1.7 (July 2022), implemented Henderson's AI after Jensen, Madsen and Thompson (1997). In 4.4.1 (April 2025), the two engines were unified in mmes. In 4.4.81 (September 2026), mmer came back alongside mmes. In 4.5.12, mmes and mmer are exported. mmec is no longer exported. The code of mmer is gathered in a file DEPRECATED_mmer.R, "pending deprecation" [28].

The mmes engine solves the MME by default with a sparse AI-REML [28]. The sparse LDLᵀ factorisation is reused for the solutions and the derivatives. The required elements of the inverse come from the Takahashi recursion [28]. CHOLMOD is available, with optional METIS ordering or AMD by default [28]. The direct formulation (formulation = "direct") works on the dense covariance of the observations [28]. Both engines have an EM, GEM or PX-GEM recovery when AI fails [28]. Convergence requires a small projected score. A rejected step returns convergence = FALSE [28]. A PCG solver with a stochastic trace estimator (Hutchinson) exists [28]. The solveOnly mode solves the MME at known variances without forming C. The CHANGELOG reports about 10 s for one million animals with one trait and 40 s with three traits [28].

Structures combine into arbitrary Kronecker products through vsm() [28]. The spatial vignette shows an AR1⊗AR1 field as a random effect or in the residual [28]. K or K⁻¹ matrices are entered through Gu, the inverse being flagged by an attribute [28]. Since 2026, rrm() fits a true reduced rank with zero specific variances, like rr() in ASReml [28]. A build compiled with CUDA speeds up some dense kernels of the Henderson engine. The sparse work stays on CPU, and so does the direct formulation [28]. The pace of change is very fast: version 4.4.90 broke the numerical control interface [28].

#### 2.10 lme4 and its genetic extensions

lme4 is under GPL ≥ 2 [29]. Version 2.1-0 appeared on 23 October 2026 [29]. The engine solves a penalised least squares problem by sparse Cholesky, and optimises the profiled likelihood with a derivative-free optimiser [30]. Since version 2.0-1 (March 2026), lme4 accepts us, diag, cs and ar1 structures on the covariance of a random term [31]. lme4 does not accept a relationship matrix natively.

The genetic extensions replace Z by Z·L, where L is a factor of the relationship matrix [33][38]. lme4qtl follows this principle for arbitrary K matrices [33]. Its last commit dates from January 2021 [32]. pedigreemm does it for the pedigree [35][36]. Its latest version dates from October 2024 [34]. lme4breeding, by the author of sommer, adds a rotation by eigendecomposition restricted to a single relationship matrix [38]. Its version 1.1.5 dates from 5 October 2026 [37]. None of these packages offers a native multivariate us model in the ASReml sense.

#### 2.11 breedR

breedR is an R interface under GPL-3 that calls the BLUPF90 binaries (REMLF90, AIREMLF90, Gibbs) [39]. It targets forestry trials. It offers autoregressive spatial effects and splines, as well as a competition effect between neighbours [39]. The README states that the package is no longer actively maintained [39]. It is not on CRAN.

#### 2.12 glmmTMB and spaMM

glmmTMB is under AGPL-3 [40]. It fits by maximum likelihood with Laplace integration and gradients by automatic differentiation (TMB) [40]. The option REML = TRUE exists [41]. The structures include us, toep, cs, diag, ar1, hetar1, ou, exp, gau, mat, rr, propto and equalto [41]. propto fits a covariance proportional to a known matrix, hence a dense relationship matrix K [41]. No ar1⊗ar1 Kronecker product is documented in the vignette.

spaMM is under CeCILL-2 [42]. It fits models with spatially correlated effects, multivariate responses and non-Gaussian effects [42]. The likelihood or restricted likelihood uses Laplace approximations [42]. Markov fields and autoregressive models exploit sparse precision matrices [42].

#### 2.13 Bayesian routes: MCMCglmm, brms, R-INLA

MCMCglmm is under GPL ≥ 2 [43]. It fits multivariate GLMMs by MCMC [44]. The ginverse argument takes sparse inverses proportional to the covariance of a random effect, typically A⁻¹ [45]. Its latest version dates from May 2024 [43].

brms is under GPL-2 and relies on Stan (HMC, NUTS) [46][47]. It fits multivariate models. A known covariance is entered through gr(cov = ). OpenCL computation is available with the cmdstanr backend, and not with rstan [48].

R-INLA is under the MIT licence in its public repository [49]. It implements the integrated nested Laplace approximation and not REML [50]. The latent model generic0 takes a supplied precision matrix, hence K⁻¹ [51]. The ar1 model and the group construction give separable structures [51]. The version of 26 October 2026 is 26.10.26 [49].

#### 2.14 MixedModels.jl and statsmodels

MixedModels.jl is under the MIT licence [52]. Version 5.9.1 dates from 2 October 2026 [52]. Since v5, optimisation is unconstrained and the default optimiser is NEWUOA from NLopt [52]. The package follows the PLS approach of lme4. It supports neither a relationship matrix nor a native multivariate us structure.

statsmodels is under BSD-3 [53]. MixedLM is built on independent groups [54]. Crossed effects require treating the whole data set as a single group [54]. The likelihood, gradient and Hessian computations follow Lindstrom and Bates (1988) [54][55]. This module does not suit genetic models with a relationship matrix.

### 3. Published comparisons of speed and accuracy

All the comparisons found were carried out by the authors of one of the compared programs. No recent independent comparison includes ASReml, BLUPF90, WOMBAT and DMU together on large genomic data.

| Study | Software | Data | Figures | Caveat |
|---|---|---|---|---|
| Covarrubias-Pazaran 2016 [59] | sommer, ASReml, rrBLUP, regress, EMMREML, BGLR, MCMCglmm | GBLUP with 1 component, N = 5 000; 3 dense components, N = 10 585 hybrids | 1 component: sommer 240 s, ASReml 443 s, EMMREML 611 s, rrBLUP 766 s, BGLR 1 292 s, regress 1 684 s, MCMCglmm 533 556 s. 3 components: sommer 426 s, regress 1 868 s, EMMREML 1 137 s, ASReml 6 874 s, BGLR 17 364 s, MCMCglmm more than 12 days | ASReml reported an abnormal termination in both cases. The inversion of the 10 585 × 10 585 matrix required by ginverse dominates the time of ASReml and MCMCglmm. i7 3.4 GHz machine, 16 GB |
| Lee and van der Werf 2016 [21] | MTG2, GEMMA, ASReml, WOMBAT | Mice, n = 1 908, up to 5 traits; ARIC humans, n = 7 263, 8 traits | Mice: MTG2 a few seconds, "a few thousand times" faster than ASReml and WOMBAT. ARIC: ASReml and WOMBAT too slow to run. Estimates "essentially identical" | A single GRM, design balanced by imputation. Reported gain above 1 000 times |
| Masuda et al. 2015 [13] | AIREMLF90 with FSPAK against YAMS | 23 models (chicken, pig, beef and dairy cattle), including single-step | YAMS on average more than 10 times faster than FSPAK; trace computation up to 67 times faster | Comparison within the same program |
| Masuda et al. 2014, WCGALP [60] | AIREMLF90 with FSPAK against YAMS | More than 15 000 genotyped animals | 1 trait: first 5 iterations in more than 4 h (FSPAK) against 20 min (YAMS). 4 traits: FSPAK fails, 5 h with YAMS | Figures read in a search excerpt; PDF unreachable, unverified |
| Yin et al. 2023 (HIBLUP) [16] | HIBLUP, BLUPF90, ASReml-SA, GCTA, MTG2, LDAK, sommer | Simulation of 10 000 individuals; 500 000 individuals for GCTA | Near-identical estimates. ASReml converges in 9 iterations against 15 for BLUPF90, for a similar total time. sommer requires much time and memory. BLUPF90 has the edge when few phenotyped individuals are genotyped. GRM of 500 000 individuals: 21 h for HIBLUP against 56 h for GCTA | Size limited by the BLUPF90 licence. DMU excluded for lack of parallelism. Detailed figures in tables not reread |
| Covarrubias-Pazaran 2024 (lme4breeding, preprint) [38] | lme4breeding, lme4qtl 0.2.2, sommer 4.3.4, ASReml 4.1.0.182 | Wheat, 599 lines, 4 environments | Components "almost identical". Gain of 4 to 160 times for the diagonal model. Ratios of 1:1 350 against ASReml and 1:1 261 against lme4qtl for the unstructured model | Mac M1, 8 GB. Unpublished preprint. Ratios cited without the underlying times |
| Caamal-Pat et al. 2021 (lme4GS) [61] | lme4GS, BGLR 0.8, sommer 4.1.3 | Genomic kernel models | Model 6: sommer 46 s, lme4GS 81.5 s, BGLR 143 s. Model 7: lme4GS 1.3 s, sommer 2.7 s, BGLR 20.2 s | No comparison with ASReml |
| Covarrubias-Pazaran 2018 (preprint) [62] | sommer, ASReml-R, MCMCglmm | BTdata and CPdata data sets | According to an excerpt: ASReml-R fastest on BTdata, sommer slowest; sommer faster on CPdata; sommer slow when p < n | Figures not reread in the full text |

Two developer announcements complete this table without being comparisons. Echidna reports 24 minutes per iteration with an H⁻¹ of 201 150 individuals and 40 000 genotyped animals [7]. sommer reports about 10 s for one million animals in solveOnly mode, without variance estimation [28].

A poster by Rogers and Taylor (IBS Australasia, 2019) reportedly compares R mixed model packages for plant breeding. I did not find its content. The figures attributed to it online remain unverified.

### 4. Sources

1. VSN International, ASReml-R page. https://vsni.co.uk/software/asreml-r (accessed 26/10/2026).
2. Butler D.G., Cullis B.R., Gilmour A.R., Gogel B.J., Thompson R. (2023). ASReml-R Reference Manual Version 4.2. VSN International. https://asreml.kb.vsni.co.uk/wp-content/uploads/sites/3/ASReml-R-Reference-Manual-4.2.pdf (sections 1.1, 2.2, A.1, A.2).
3. Gilmour A.R., Thompson R., Cullis B.R. (1995). Average information REML: an efficient algorithm for variance parameter estimation in linear mixed models. Biometrics 51:1440–1450. https://doi.org/10.2307/2533274
4. Gilmour A.R., Gogel B.J., Cullis B.R., Welham S.J., Thompson R. (2022). ASReml User Guide Release 4.2 Functional Specification. https://asreml.kb.vsni.co.uk/wp-content/uploads/sites/3/ASReml-4.2-Functional-Specification.pdf (sections 7.11, 11.3.6, option !MP).
5. VSN International, Release notes ASReml-R 4.2. https://asreml.kb.vsni.co.uk/?p=2374
6. Echidna Mixed Models Software, home page. https://www.echidnamms.org/
7. Gilmour A.R., note "Echidna Build 2.15 October 22, 2026". https://www.echidnamms.org/echidna-build-2-15-october-22-2026/
8. Meyer K. (2017). WOMBAT User Notes, version of 31 May 2017. INRAE copy: https://forge-dga.jouy.inra.fr/attachments/download/8063/WombatManual-V20170531.pdf (chapters 1, 2, 5.3). Official site http://didgeridoo.une.edu.au/km/wombat.php unreachable on 26/10/2026.
9. BLUPF90 wiki, "Conditions of use". https://nce.ads.uga.edu/wiki/doku.php?id=condition_of_use
10. BLUPF90 wiki, "Distribution". https://nce.ads.uga.edu/wiki/doku.php?id=distribution
11. BLUPF90 wiki, "Application programs". https://nce.ads.uga.edu/wiki/doku.php?id=application_programs
12. BLUPF90 wiki, home page. https://nce.ads.uga.edu/wiki/doku.php?id=start
13. Masuda Y., Aguilar I., Tsuruta S., Misztal I. (2015). Technical note: Acceleration of sparse operations for average-information REML analyses with supernodal methods and sparse-storage refinements. J. Anim. Sci. 93:4670–4674. https://doi.org/10.2527/jas.2015-9395
14. Lourenco D., Legarra A., Tsuruta S., Masuda Y., Aguilar I., Misztal I. (2020). Single-step genomic evaluations from theory to practice: using SNP chips and sequence data in BLUPF90. Genes 11:790. https://doi.org/10.3390/genes11070790
15. Madsen P., Jensen J. DMU, a package for analyzing multivariate mixed models, guide v6 release 5.6 (2024). https://dmu.ghpc.au.dk/dmu/DMU/Doc/Current/dmuv6_guide.5.6.pdf. Server unreachable on 26/10/2026; content known from search engine excerpts only.
16. Yin L. et al. (2023). HIBLUP: an integration of statistical models on the BLUP framework for efficient genetic evaluation using big genomic data. Nucleic Acids Res. 51:3501–3512. https://doi.org/10.1093/nar/gkad074
17. Luke, MiX99 page. https://www.luke.fi/en/node/8911
18. MiXBLUP, "Licenses" page. https://www.mixblup.eu/licenses
19. MiXBLUP, home page. https://www.mixblup.eu/
20. Lee S.H., Zhou X. (2020). CORE GREML: MTG2 code and instruction. Dryad. https://doi.org/10.5061/dryad.bk3j9kd8c
21. Lee S.H., van der Werf J.H.J. (2016). MTG2: an efficient algorithm for multivariate linear mixed model analysis based on genomic information. Bioinformatics 32:1420–1422. https://doi.org/10.1093/bioinformatics/btw012
22. GCTA repository. https://github.com/jianyangqt/gcta (README, tags, commits).
23. Yang J., Lee S.H., Goddard M.E., Visscher P.M. (2011). GCTA: a tool for genome-wide complex trait analysis. Am. J. Hum. Genet. 88:76–82. https://doi.org/10.1016/j.ajhg.2010.11.011
24. See [21], introduction: direct AI is implemented in GCTA-GREML.
25. Jiang L. et al. (2019). A resource-efficient tool for mixed model association analysis of large-scale data. Nat. Genet. 51:1749–1755. https://doi.org/10.1038/s41588-019-0530-8
26. LDAK repository, file LICENCE.txt. https://github.com/dougspeed/LDAK
27. Speed D., Hemani G., Johnson M.R., Balding D.J. (2012). Improved heritability estimation from genome-wide SNPs. Am. J. Hum. Genet. 91:1011–1021. https://doi.org/10.1016/j.ajhg.2012.10.010
28. sommer 4.5.12, CRAN. https://cran.r-project.org/package=sommer. Sources read: DESCRIPTION, NAMESPACE, CHANGELOG, README, vignettes sommer.qg and sommer.spatial, R/FUN_mme_gpu.R, R/FUN_mme_control.R (archive sommer_4.5.12.tar.gz).
29. lme4 2.1-0, CRAN. https://cran.r-project.org/package=lme4
30. Bates D., Mächler M., Bolker B., Walker S. (2015). Fitting linear mixed-effects models using lme4. J. Stat. Softw. 67(1). https://doi.org/10.18637/jss.v067.i01
31. lme4, NEWS file, versions 2.0-1 and 2.1-0. https://github.com/lme4/lme4/blob/master/inst/NEWS.Rd
32. lme4qtl repository. https://github.com/variani/lme4qtl
33. Ziyatdinov A. et al. (2018). lme4qtl: linear mixed models with flexible covariance structure for genetic studies of related individuals. BMC Bioinformatics 19:68. https://doi.org/10.1186/s12859-018-2057-x
34. pedigreemm 0.3-5, CRAN. https://cran.r-project.org/package=pedigreemm
35. Vazquez A.I., Bates D.M., Rosa G.J.M., Gianola D., Weigel K.A. (2010). Technical note: An R package for fitting generalized linear mixed models in animal breeding. J. Anim. Sci. 88:497–504. https://doi.org/10.2527/jas.2009-1952
36. See [59], table 2: pedigreemm accepts a pedigree but not a covariance matrix supplied directly.
37. lme4breeding 1.1.5, CRAN. https://cran.r-project.org/package=lme4breeding
38. Covarrubias-Pazaran G. (2024). lme4breeding: enabling genetic evaluation in the era of genomic data. bioRxiv. https://doi.org/10.1101/2024.05.17.594795
39. breedR repository (README, DESCRIPTION, NAMESPACE, R/competition.R, R/spatial.R). https://github.com/famuvie/breedR
40. glmmTMB 1.1.15.2, CRAN. https://cran.r-project.org/package=glmmTMB ; Brooks M.E. et al. (2017). glmmTMB balances speed and flexibility among packages for zero-inflated generalized linear mixed modeling. R J. 9:378. https://doi.org/10.32614/RJ-2017-066
41. glmmTMB, vignette "Covariance structures with glmmTMB". https://github.com/glmmTMB/glmmTMB/blob/master/glmmTMB/vignettes/covstruct.rmd
42. spaMM 4.7.17, CRAN. https://cran.r-project.org/package=spaMM ; Rousset F., Ferdy J.-B. (2014). Testing environmental and genetic effects in the presence of spatial autocorrelation. Ecography 37:781–790. https://doi.org/10.1111/ecog.00566
43. MCMCglmm 2.36, CRAN. https://cran.r-project.org/package=MCMCglmm
44. Hadfield J.D. (2010). MCMC methods for multi-response generalized linear mixed models: the MCMCglmm R package. J. Stat. Softw. 33(2). https://doi.org/10.18637/jss.v033.i02
45. MCMCglmm, help page MCMCglmm.Rd, argument ginverse. https://github.com/jarrodhadfield/MCMCglmm
46. brms 2.23.0, CRAN. https://cran.r-project.org/package=brms
47. Bürkner P.-C. (2017). brms: an R package for Bayesian multilevel models using Stan. J. Stat. Softw. 80(1). https://doi.org/10.18637/jss.v080.i01
48. brms, file R/backends.R (handling of OpenCL by backend). https://github.com/paul-buerkner/brms
49. R-INLA repository (LICENSE, version 26.10.26). https://github.com/hrue/r-inla
50. Rue H., Martino S., Chopin N. (2009). Approximate Bayesian inference for latent Gaussian models by using integrated nested Laplace approximations. J. R. Stat. Soc. B 71:319–392. https://doi.org/10.1111/j.1467-9868.2008.00700.x
51. R-INLA documentation of the latent models generic0 and ar1. https://github.com/hrue/r-inla/tree/devel/r-inla.org/doc/latent
52. MixedModels.jl repository (README, version 5.9.1). https://github.com/JuliaStats/MixedModels.jl
53. statsmodels repository (licence, version 0.15.0). https://github.com/statsmodels/statsmodels
54. statsmodels, documentation "Linear Mixed Effects Models". https://www.statsmodels.org/stable/mixed_linear.html
55. Lindstrom M.J., Bates D.M. (1988). Newton-Raphson and EM algorithms for linear mixed-effects models for repeated-measures data. J. Am. Stat. Assoc. 83:1014–1022. https://doi.org/10.1080/01621459.1988.10478693
56. Meyer K. (2007). WOMBAT: a tool for mixed model analyses in quantitative genetics by restricted maximum likelihood (REML). J. Zhejiang Univ. Sci. B 8:815–821. https://doi.org/10.1631/jzus.2007.B0815
57. Meyer K. (2020). WOMBAT: a tool for estimation of genetic parameters, highlights and updates. ICQG 6, abstracts, p. 199. https://hdl.handle.net/1959.11/30741
58. Lee S.H., van der Werf J.H.J. (2006). An efficient variance component approach implementing an average information REML suitable for combined LD and linkage mapping with a general complex pedigree. Genet. Sel. Evol. 38:25. https://doi.org/10.1186/1297-9686-38-1-25
59. Covarrubias-Pazaran G. (2016). Genome-assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11:e0156744. https://doi.org/10.1371/journal.pone.0156744
60. Masuda Y. et al. (2014). Communication at the 10th WCGALP on YAMS in AIREMLF90. INIA copy: https://inia.uy/sites/default/files/publications/2024-10/Masuda-et-al2014WCGALP.pdf (access denied, unverified).
61. Caamal-Pat D. et al. (2021). lme4GS: an R-package for genomic selection. Front. Genet. 12:680569. https://doi.org/10.3389/fgene.2021.680569
62. Covarrubias-Pazaran G. (2018). Software update: moving the R package sommer to multivariate mixed models for genome-assisted prediction. bioRxiv. https://doi.org/10.1101/354639

## Appendix B. Accelerated computing for REML and linear mixed models, 2010-2026

Consulted on 26 October 2026. DOIs were verified on Crossref, or on DataCite for arXiv. Repository activity comes from the GitHub API. For pyreml, the 0.5.3 source code was read in full. The word "unverified" flags a point I could not confirm on a primary source.

### Comparison table

Column key. "Computation" states whether the tool estimates the variances by exact REML, by an approximation, or only does BLUP or GWAS. "Model" summarises generality. "Scale" gives the largest published case and its hardware.

| Tool | Computation | Model | Computational and memory technique | Scale and hardware | Licence, maintenance |
|---|---|---|---|---|---|
| pyreml 0.5.3 (Marchal et al.) | Exact REML, analytical gradient, L-BFGS then Adam as a fallback. No AI. | Multiple terms, multivariate, Σ⊗K, fa, us, diag, blocks, heteroscedasticity, spatial kernels, random regression. | Dense PyTorch. V form or Woodbury (dense q×q C = G⁻¹ + Z'R⁻¹Z). Index grids instead of Z. FP32 then FP64. | n = 10⁶, q = 25 000, 41 parameters, 262 to 279 s on an 80 GB A100. | GPL-3.0. Created in June 2026, 19 versions, 0.5.3 on 21/10/2026. |
| kmmr (Melsen et al.) | Exact AI-REML. | Multivariate, us, extended FA, AR1×AR1, Gaussian kernels, Kronecker, stacking. | Sparse part by supernodal CHOLMOD on CPU, Takahashi partial inverse. Dense Schur complement on GPU (cuSOLVER). | GBLUP with 10 375 genotypes: 49 s on A100 against 989 s for ASReml-R. | GPL-3.0, R. Active (push on 22/10/2026). |
| SAIGE-GPU | AI-REML on penalised quasi-likelihood (biased according to the SAIGE authors), then tests. | A single genetic term (full or sparse GRM). | PCG without forming the GRM, genotype × vector products by cuBLAS, multi-node MPI, random traces. | 635 969 individuals, 2 068 phenotypes. Summit (V100), Frontier (MI250X), A100. | GPL-3.0. Commit in April 2026. |
| BOLT-REML | Approximate REML: Monte Carlo gradient, AI, conjugate gradient. | Several components and several traits. | Memory ≈ MN/4 bytes. CPU only. | ≈ 50 000 individuals in the article. UK Biobank in a few days on 8 threads (manual). | GPLv3. Version 2.5 in June 2025. |
| DISSECT | Exact REML on a distributed dense V. | Univariate and bivariate LMM. Multiple components announced as in testing in 2015. | MPI, BLACS, ScaLAPACK. | 470 000 individuals, 590 004 SNPs, under 4 h on 8 400 cores (ARCHER). | GPLv3. Maintenance unverified. |
| REACTA | Exact REML, regional heritability. | Derived from GCTA. | MAGMA and CUDA, regions distributed by MPI. | 1 024 GPUs of Titan, 15 min. | Licence unverified. |
| GPU_REML (repository Asiandier/SPLIT) | Approximate REML: block PCG, Hutchinson, SLQ, AI or Fisher. | Several GRMs, one variance per GRM. | GRM never formed, compressed genotypes streamed to the GPU. | Unpublished. | No licence. Created in May 2026. No article. |
| TorchLIMIX | LIMIX multivariate LMM for GWAS. Estimation mode unverified. | Two or more traits, one relationship matrix. | Batched algebra on GPU (PyTorch). | Gain up to 100× over LIMIX. | Announced repository returns 404. |
| Quickdraws | Spike-and-slab by stochastic variational inference. Heritability by RHE, not by REML. | One trait at a time. | Mini-batches on GPU, genotypes on 2 bits. | 405 088 individuals on a 16 GB A10G. | GPL-3.0. Commit in August 2025. |
| Faster-LMM-D | LMM with one relationship matrix for GWAS, AI-REML taken from GEMMA. | One component. | GPU only for the dense matrix product. | Unverified. | GPL-3.0. Last commit in 2018. |
| FaST-LMM, pylmm, GEMMA | Exact REML by spectral decomposition. | One relationship matrix. mvLMM in GEMMA. | CPU. | 120 000 individuals (FaST-LMM). | Apache-2.0, MIT, GPL-3.0. No GPU. |
| Hail LinearMixedModel | Exact REML. | One K matrix. | Spectral rotation. Only the tests are distributed (Spark). | Unverified. | MIT. Removed since Hail 0.2.94. |
| REGENIE, GloWGR | Stacked ridge, no REML. | GWAS. | Segments in memory. Spark for GloWGR. | 407 746 individuals. | MIT, Apache-2.0. No GPU. |
| tensorQTL, jaxQTL | Linear models or GLMs, no random effect. | QTL. | PyTorch GPU; JAX. | Gain > 200× (tensorQTL). | BSD-3, MIT. |
| GPyTorch (BBMM) | GP marginal likelihood, approximated by CG and stochastic Lanczos beyond 800 points. No REML found. | Gaussian processes, structured kernels. | Batched matrix-matrix products, preconditioner, multi-GPU. | 10⁶ exact points in under 2 h (Wang et al. 2019). | MIT. Active. |
| KeOps | Kernel reductions without storing the n×n matrix. | Computational building block. | Generated CUDA code, autodiff. | Millions of points in seconds (abstract). | MIT. Active. |
| GPflow, GPJax, tinygp | GP likelihood, variational. | Gaussian processes. | TensorFlow or JAX. | Unverified. | Apache-2.0, MIT. |
| TensorFlow Probability | ML by Monte Carlo EM (HMC + Adam) in the LMM example. No REML. | 3 iid random factors (InstEval). | XLA on GPU. | 58 737 observations. | Apache-2.0. Active. |
| NumPyro, BlackJAX | Bayesian (NUTS, HMC, VI). | Any model written by the user. | JAX, vectorised chains, GPU and TPU. | Not applicable. | Apache-2.0. Active. |
| Stan with OpenCL | Bayesian. | Any Stan model. | Cholesky, its derivative and a few likelihoods on GPU. | Unverified. | BSD-3. Active. |
| ExaGeoStat | Exact Gaussian likelihood (Matérn), plus TLR approximations. | Spatial. | Tiled dense Cholesky, GPU and distributed memory. | Unverified in the text. | Unverified. |
| INLA-DIST, DALIA | Approximate Bayesian (INLA). | Spatio-temporal fields, multivariate for DALIA. | Block sparse factorisation on GPU, selected inversion (Serinv). | 496 GH200 (Alps) for DALIA. | Unverified. |
| jaxmg | Building block: distributed Cholesky, solve and log-determinant. | Dense matrix. | cuSOLVERMp, multi-GPU and multi-node, callable under jit. | Up to the combined memory of the GPUs. | Apache-2.0. Active. |
| miraculix in MiXBLUP 3.1 | BLUP by PCG (single-step), no REML. | Multi-trait. | Genotype × vector product on GPU, 2 bits, CUTLASS. | 2.61 M genotyped, 26.46 M in the pedigree, A100. | Apache 2.0 announced. MiXBLUP commercial. |
| BOLT (Theta Solutions) | BLUP by PCG and single-step Gibbs. No documented REML. | Bayes A, B, C, LASSO. | Several GPUs. | 13 M animals. | Proprietary. |
| MC-REML (Matilainen; Bermann) | EM-REML and AI-REML with Monte Carlo or Hutchinson traces. | Multi-trait, single-step. | CPU, iterative solvers. | 7 M animals, 330 000 genotyped, 121 h, 53 GB. | In BLUPF90 (2025). |
| SLEMM, SL_REML | REML by stochastic Lanczos. | Two components (SL_REML). SLEMM with SNP weights. | CPU. | 3 M simulated individuals (SLEMM). | Public code. |
| RHE-mc, REHE | Method of moments (Haseman-Elston), not REML. | Up to 100 components (RHE-mc). | Random traces. | 10⁶ individuals, 10⁶ SNPs in a few hours. | Unverified. |
| HIBLUP, MTG2 | Exact AI-REML (MME or V form). | Multi-trait. | Multithreaded CPU. | 500 000 simulated individuals (HIBLUP). | Free executables. No GPU. |

### General REML solvers on GPU

#### pyreml

- Article. Marchal, Trottier, Bry and Raimondi, bioRxiv preprint of 24 October 2026, doi:10.64898/2026.10.20.757556. Software archived under doi:10.5281/zenodo.21394319. The full text was read.
- Package. PyPI, from 0.1.0 (24 June 2026) to 0.5.3 (21 October 2026), 19 versions (https://pypi.org/pypi/pyreml/json). GPL-3.0-or-later. Repository with no star and no open issue on 26 October (GitHub API). Authors at IGMM and IMAG, Montpellier.
- Optimisation. Gradient descent on −2ℓ_REML, following Simchoni and Rosset (JMLR 2023, 24(156)), with the REML term added. PyTorch L-BFGS with a strong Wolfe search, and Adam as a fallback when the line search fails (pyreml/Optimizer.py). β is optimised as a parameter together with the variances (MixedModel.py, l. 297 and 304). The stopping criterion counts the remaining significant digits of the likelihood (doc/training.qmd).
- Gradient. Autodiff through the Cholesky factor is replaced by the closed-form REML gradient, injected through a "ghost loss" (MixedModel.py, REML_loss, l. 709-720). The overlap with the custom_vjp of remlax is direct.
- V or Woodbury form. The SMW route factorises C = G⁻¹ + Z'R⁻¹Z, dense of size q×q, and never forms V. It is used when q < n and R is diagonal (MixedModel.py, l. 207-225; Variance.py, class Capacitance). Otherwise a Cholesky factorisation of the dense V is used (class DirectSolve).
- Memory. The Kronecker blocks and ZGZ' are never formed. Z and R are replaced by index grids. A test checks that the largest allocation stays below 4 × max(q², np) (tests/test_biggplup.py, n = 23 600).
- Precision. Three modes: FP64, FP32, or mixed (FP32 then FP64), mixed by default. In FP32, 241 functional tests fail (preprint, section 2.3).
- Absences verified in the 0.5.3 code. No call to torch.sparse: everything is dense. No AI algorithm and no information matrix. No standard error for the variance components in the documentation. Only the fixed effects (SE) and the BLUPs (PEV, CD) get them, through the dense inverse of Henderson's equations (MixedModel.py, l. 817-830). No GLMM (preprint, discussion). No multi-GPU. No built-in likelihood ratio test.
- Structures. Eight left structures (iid, diag, full, fa, bl_resp, bl_form, kr_resp, kr_form) and seven right structures (iid, str, het, dist, eucl, ar_iso, ar_ani), GaussianComponents.py, l. 10-11. A low-level constructor accepts arbitrary varmeth functions.
- Published benchmarks. GBLUP compared with WOMBAT (EM, AI, PX-AI on 64 EPYC 9124 cores) and not with ASReml. More than 12× faster on A100, more than 2× on L40S. WOMBAT AI did not converge for n ≥ 10 000. Extreme case: 25 000 genotypes × 40 environments, 10⁶ observations, 279 s in FP64 on A100 (preprint). The documentation shows 261.65 s for a model of the same size with 41 parameters (doc/building_models.qmd). FA2 on 8 traits and 942 individuals: 429 s for pyreml against 281 s for WOMBAT, or 66 s at equal tolerance. Validation against lme4, rrBLUP, BLUPF90 and spaMM (505 tests).
- The authors describe their package as "the first library performing REML estimations on GPU". REACTA (2014) and SAIGE-GPU contradict this priority in the strict sense.

#### kmmr

- Melsen, van Eeuwijk and Peeters, arXiv:2610.16671 (October 2026), doi:10.48550/arXiv.2610.16671. Repository github.com/KillianMelsen/kmmr, GPL-3.0, R, archive doi:10.5281/zenodo.23398636. The package is labelled experimental, Linux or WSL2 with CUDA.
- Multivariate AI-REML on the mixed model equations. The sparse equations are eliminated by supernodal Cholesky (CHOLMOD). The dense Schur complement runs on GPU: potrf, potrs and potri from cuSOLVER, SpMM from cuSPARSE. The partial inverse of the sparse part comes from a supernodal Takahashi recursion (cholmod-extra).
- Benchmarks on an 80 GB A100 and an i7. Wheat GBLUP with 10 375 genotypes: ASReml-R 988.6 s, sommer 511.5 s, kmmr 129.3 s on 8 threads and 49.2 s on GPU. FA + stacked kernel model: ASReml-R does not converge in 50 iterations, kmmr 6.8 s on GPU.
- Limits stated by the authors. The memory of the Schur complement grows as the square of the number of dense equations. The V form on GPU is cited as a natural extension. On single-stage models, kmmr converges to lower log-likelihoods than ASReml-R (−3948.66 against −3937.65 for FA1).

#### GPU_REML (repository Asiandier/SPLIT)

- Repository created on 30 May 2026, with no licence and no article (https://github.com/Asiandier/SPLIT).
- REML for V = Σθ_g K_g + θ_e I, with K_g = Z_g Z_g'/m applied without forming the GRM. Compressed genotypes streamed to the GPU. Block PCG, Hutchinson traces, log-determinants by SLQ, constrained AI or Fisher updates (README).
- A research tool for heritability partitioning. No published benchmark.

#### TorchLIMIX

- Horn, Nikoloski and Lippert, bioRxiv 2026, doi:10.64898/2026.04.29.721342. A PyTorch port of the LIMIX multivariate LMM, with a gain of up to two orders of magnitude.
- The repository github.com/bi-horn/torchLIMIX returns 404. The covariance estimation mode is unverified.

### Human genetics and GWAS

#### SAIGE-GPU

- Rodriguez et al. 2026, Bioinformatics 42(3):btag032, doi:10.1093/bioinformatics/btag032.
- The GPU speeds up step 1, the fit of the null GLMM. The genotype × vector products of the PCG go through cuBLAS without forming the GRM. The genotype columns are distributed across nodes by MPI. Step 2 stays on CPU.
- Step 1 of SAIGE is an AI-REML on penalised quasi-likelihood, with random traces. Zhou et al. 2018 (doi:10.1038/s41588-018-0184-y) advise against using its τ to estimate heritability.
- 635 969 MVP participants and 2 068 phenotypes. About 72× per PCG iteration, 3 to 5× on the full fit. Summit, Frontier, A100 in the cloud.
- Code: branch SAIGE-GPU-1.3.3 of github.com/saigegit/SAIGE, GPL-3.0, archive doi:10.5281/zenodo.17642591.

#### Quickdraws

- Loya et al. 2025, Nat Genet, doi:10.1038/s41588-024-02044-7. Spike-and-slab regression by stochastic variational inference, in mini-batches on GPU (PyTorch).
- Heritability comes from randomised Haseman-Elston, not from REML. One trait at a time. 405 088 individuals on a 16 GB A10G. github.com/PalamaraLab/quickdraws, GPL-3.0.

#### Faster-LMM-D, pylmm, FaST-LMM, GEMMA

- Faster-LMM-D (github.com/genetics-statistics/faster_lmm_d, GPL-3.0, last commit in December 2018). The GPU (cuBLAS or ArrayFire) replaces only the dense matrix product (source/faster_lmm_d/optmatrix.d). No associated article.
- pylmm (github.com/nickFurlotte/pylmm). One component, a grid on h then Brent, in REML. Rewritten in October 2026 (v1.0.0, MIT). No GPU.
- FaST-LMM (Lippert et al. 2011, doi:10.1038/nmeth.1681) and GEMMA (Zhou and Stephens 2012, doi:10.1038/ng.2310; mvLMM 2014, doi:10.1038/nmeth.2848). Spectral decomposition of a single relationship matrix. No GPU.

#### BOLT-LMM and BOLT-REML

- BOLT-REML: Loh et al. 2015, Nat Genet 47:1385, doi:10.1038/ng.3431 (PMC4666835). The gradient is estimated by Monte Carlo, the Hessian by the AI matrix, and the systems by conjugate gradients. 100 random vectors in general, 15 for regional analyses.
- Accuracy. Estimates almost identical to GCTA with 100 trials on 50 000 simulated traits. With 15 trials, the standard errors increase by a factor of 1.03.
- Cost ≈ O(MN^1.5) and memory ≈ MN/4 bytes, almost independent of the number of components. 22 chromosomal components and bivariate pairs in the article. CPU only (manual v2.5, https://alkesgroup.broadinstitute.org/BOLT-LMM/BOLT-LMM_manual.html).

#### GCTA, fastGWA, REGENIE, Hail, tensorQTL, jaxQTL

- GCTA (Yang et al. 2011, doi:10.1016/j.ajhg.2010.11.011) and fastGWA (Jiang et al. 2019, doi:10.1038/s41588-019-0530-8): CPU, MKL or OpenBLAS. No GPU route found.
- REGENIE (Mbatchou et al. 2021, doi:10.1038/s41588-021-00870-7): stacked ridge, no REML, no GPU. GloWGR distributes it on Spark (github.com/projectglow/glow).
- Hail LinearMixedModel: a single K matrix, exact REML on the main node, distributed tests. Removed since version 0.2.94 (Hail documentation).
- tensorQTL (Taylor-Weiner et al. 2019, doi:10.1186/s13059-019-1836-7) and jaxQTL (Zhang et al. 2025, medRxiv, doi:10.1101/2025.01.18.25320755): no random effect.

### Distributed REML

#### DISSECT

- Canela-Xandri et al. 2015, Nat Commun 6:10162, doi:10.1038/ncomms10162. C++ with MPI, BLACS and ScaLAPACK.
- Exact REML in O(N³) on a distributed dense V. Univariate and bivariate LMM. 470 000 individuals and 590 004 SNPs in under 4 h on 8 400 cores and about 16 TB of RAM (ARCHER), REML and BLUP included.
- GPLv3 on www.dissect.ed.ac.uk. The certificate had expired at the time of consultation. Maintenance unverified.

#### REACTA

- Cebamanos et al. 2014, Bioinformatics, doi:10.1093/bioinformatics/btt754. Derived from GCTA for regional heritability.
- GRM and REML inversion by MAGMA, CUDA kernels, regions distributed by MPI. 19 746 regions on 1 024 GPUs of Titan (ORNL) in 15 min 26 s. About 3× on a K20.

#### ExaGeoStat

- Abdulah et al. 2018, IEEE TPDS 29(12), doi:10.1109/TPDS.2018.2850749 (arXiv:1708.02835). Exact Gaussian likelihood of a Matérn field by tiled dense Cholesky, on NVIDIA GPUs, KNL and distributed memory. Approximate tile low-rank variants (unverified in this text).
- Spatial domain, with no fixed effects and no generic REML.

#### Other distributed building blocks

- jaxmg (github.com/flatironinstitute/jaxmg, Apache-2.0, active). Connects JAX to cuSOLVERMp. Provides potrs with a log-determinant, on matrices distributed over several GPUs and nodes, under jit.
- Lewis et al. 2022, PNAS, doi:10.1073/pnas.2122762119. Distributed dense linear algebra on TPU pods, in float32. Product of two matrices of size 2²⁰ in about 2 min on 2 048 TPU v3 cores.
- Wang et al. 2019, "Exact Gaussian processes on a million data points", NeurIPS (arXiv:1903.08114). Partitioned kernel products distributed over several GPUs, conjugate gradients. More than one million points in under 2 h. Number and type of GPUs unverified in the text.

### Machine learning and probabilistic libraries

#### GPyTorch and BBMM

- Gardner et al. 2018, NeurIPS, arXiv:1809.11165. A modified batched conjugate gradient (mBCG) that provides in one call the solves, the log-determinant by stochastic Lanczos, and the traces. The exact cost drops from O(n³) to O(n²). Specialised preconditioner.
- Default settings (https://docs.gpytorch.ai/en/stable/settings.html): Cholesky below 800 points, 10 vectors for the traces, 20 Lanczos iterations. Beyond that, the likelihood is therefore approximate.
- Objective: marginal likelihood of Gaussian processes. No REML term (log|X'V⁻¹X|) found in the documentation (unverified in the code). github.com/cornellius-gp/gpytorch, MIT, active. The algebra is isolated in linear_operator.

#### KeOps

- Charlier et al. 2021, JMLR 22(74):1-6 (https://www.jmlr.org/papers/v22/20-275.html). Symbolic tensors defined by a formula, reduced on GPU without storing the n×n matrix, with autodiff. Python, Matlab and R bindings.
- A computational building block for large dense kernels. github.com/getkeops/keops, MIT.

#### GPflow, GPJax, tinygp

- GPflow (Matthews et al. 2017, JMLR 18(40)) on TensorFlow, Apache-2.0. GPJax (Pinder and Dodd 2022, JOSS, doi:10.21105/joss.04455) and tinygp (github.com/dfm/tinygp, MIT) on JAX.
- Gaussian processes and variational inference. No asreml-type mixed model.

#### TensorFlow Probability

- Dillon et al. 2017, arXiv:1711.10604. The official mixed model example (https://www.tensorflow.org/probability/examples/Linear_Mixed_Effects_Models) fits InstEval with three iid random factors.
- Fit by Monte Carlo EM: HMC for the E step, Adam for the M step. Maximum likelihood, no REML. GPU recommended.

#### NumPyro, BlackJAX, Stan

- NumPyro (Phan et al. 2019, arXiv:1912.11554): iterative NUTS compiled end to end by JAX. BlackJAX (Cabezas et al. 2024, arXiv:2402.10797): composable samplers on CPU, GPU and TPU. Apache-2.0, active.
- Stan (Češnovar et al. 2019, arXiv:1907.01063): OpenCL routines for the Cholesky factorisation, its derivative and a few likelihoods. Enabled by STAN_OPENCL (https://mc-stan.org/docs/cmdstan-guide/parallelization.html).
- These tools provide no ready-made mixed model. The Bayesian route of remlax overlaps with their generic use.

#### Other JAX or PyTorch packages

- jaxgam (github.com/scatt67/jaxgam): a JAX port of mgcv, Laplace REML for additive models. A personal project declared not ready for production.
- GitHub search of 26 October 2026 ("mixed model REML jax", "REML gpu"): no other JAX library for general linear mixed models found.

### Stochastic and iterative REML

#### Trace estimators

- Hutchinson 1989, Commun. Stat. Simul. Comput., doi:10.1080/03610918908812806. With m Rademacher vectors, tr(A) ≈ the mean of z'Az. The variance equals 2(‖A‖²_F − Σ a_ii²)/m.
- Stochastic Lanczos quadrature (SLQ) for tr(f(A)), including log|A|: Ubaru, Chen and Saad 2017, SIAM J. Matrix Anal. Appl., doi:10.1137/16M1104974, with bounds on the number of vectors.

#### SL_REML (Border and Becker 2019)

- BMC Bioinformatics 20:411, doi:10.1186/s12859-019-2978-z (PMC6668092). Two algorithms, SLDF_REML and L_FOMC_REML, based on the shift invariance of Krylov subspaces.
- A single iterative pass, after which the objective is re-evaluated by vector operations. Limited to the two-component model.
- Accuracy. Heritabilities consistent with GCTA. For n = 250 000, a mean squared error of 1.24 × 10⁻⁷ against 1.74 × 10⁻⁶ for BOLT-LMM v2.3.2. Python code on github.com/rborder/SL_REML.

#### SLEMM

- Cheng et al. 2023, Bioinformatics, doi:10.1093/bioinformatics/btad127. REML by stochastic Lanczos with SNP weighting. Up to 3 M individuals and 1 M simulated SNPs. CPU. github.com/jiang18/slemm.

#### Monte Carlo REML in animal breeding

- Matilainen et al. 2012, JABG, doi:10.1111/j.1439-0388.2012.01000.x (multi-trait MC EM-REML) and 2013, PLoS One, doi:10.1371/journal.pone.0080821 (MC Newton, MC AI, MC Broyden).
- Bermann et al. 2025, GSE 57:63, doi:10.1186/s12711-025-01006-9. MC-ss-GREML with EM and AI. Same estimates as the exact method, in 14 % of the time and 1 % of the memory (100 000 animals, 14 parameters). 7 M animals in the pedigree, 330 000 of them genotyped, in 121 h and 53 GB, on CPU.
- Bermann et al. 2025, GSE 57:44, doi:10.1186/s12711-025-00991-1. Approximate log-determinant, mean relative error of about 10⁻³, 2 to 500 times faster. The maximum of the approximate likelihood coincides with the exact one.
- Strandén et al. 2024, GSE 56:73, doi:10.1186/s12711-024-00939-x. Augmented AI-REML, a single solve of the augmented equations per iteration, time per iteration reduced by 75 to 86 %.

#### Methods of moments and MM

- RHE-reg (Wu and Sankararaman 2018, doi:10.1093/bioinformatics/bty253) and RHE-mc (Pazokitoroudi et al. 2020, doi:10.1038/s41467-020-17576-9). Haseman-Elston with random traces. A hundred components on 10⁶ individuals and 10⁶ SNPs in a few hours. No REML.
- REHE (Yue et al. 2021, Genet Epidemiol 45(8), doi:10.1002/gepi.22432). Haseman-Elston constrained to positive variances, with a resampling variant. Reported accuracy comparable to REML.
- MM algorithms (Zhou, Hu, Zhou and Lange 2019, JCGS 28(2), doi:10.1080/10618600.2018.1529601). Exact ML and REML by minorisation-maximisation. Global convergence to a KKT point, faster than EM beyond two components. Extends to the multivariate case with missing data. No GPU mentioned.

### Animal and plant breeding

#### miraculix and MiXBLUP

- Freudenberg et al. 2023, Front. Genet. 14:1220408, doi:10.3389/fgene.2023.1220408. Products Z·M and Z'·M with Z stored on 2 bits, on GPU through an extension of CUTLASS 2.10 (V100, A100, H100).
- ICBF data: 2.61 M genotyped animals, 26.46 M animals in the pedigree, 16.59 M records on 6 traits. ssSNPBLUP drops from 6.47 h to 2.49 h on A100. This is BLUP by PCG, without REML.
- Integrated into the hpblup solver of MiXBLUP 3.1 (Vandenplas et al., WCGALP 2022, doi:10.3920/978-90-8686-940-4_358). MiXBLUP is commercial. Repository github.com/alexfreudenberg/miraculix, Apache 2.0 announced, last push in August 2023.

#### BOLT from Theta Solutions and Gibbs samplers

- Fernando et al. 2016, GSE 48:96, doi:10.1186/s12711-016-0273-2. Bayesian single-step and PCG on four Titan X cards. 6.18 M animals in the pedigree, about 3 h per single-trait analysis.
- WCGALP 2018 presentation (https://www.icar.org/Documents/Auckland-2018/WCGALP/1415%20Dorian%20Garrick.pdf): 13 M animals, 1 000 PCG iterations in 1.5 h. Proprietary software. No documented REML.
- Heterogeneous CPU/GPU GenSel (DeVore et al. 2014, IEEE BIBM, doi:10.1109/bibm.2014.6999186): a gain of only 1.84×, since the chain remains sequential.

#### Pioneering GPU work

- Coffey et al. 2011, Interbull Bull. 44:114-116 (https://journal.interbull.org/index.php/ib/article/view/1694): product AA' on a GTX 590, from 507 s to 70 s.
- Meyer and Tier 2013, Proc. AAABG 20:270-273 (https://hdl.handle.net/1959.11/14323): construction and inversion of G on GPU, a gain of 4 to 5×.

#### Observed absences

- No documented GPU route in BLUPF90, DMU, MiX99, JWAS (Project.toml of v2.4.0 without CUDA) or HIBLUP (Yin et al. 2023, doi:10.1093/nar/gkad074, full text without "GPU"). The lack of a result does not prove an absence.
- MTG2 (Lee and van der Werf 2016, doi:10.1093/bioinformatics/btw012): direct AI-REML on V with an eigendecomposition. No GPU version found.
- ChromaX (Younis et al. 2023, doi:10.1093/bioinformatics/btad691) simulates breeding programmes on GPU in JAX, without inference.

### Computational building blocks for the MME form on GPU

- GPU-accelerated supernodal sparse Cholesky in CHOLMOD: Rennich, Stosic and Davis 2016, Parallel Computing, doi:10.1016/j.parco.2016.06.004.
- Distributed selected inversion: PSelInv (Jacquelin, Lin and Yang 2016, ACM TOMS, doi:10.1145/2786977). Selected inversion on GPU of block tridiagonal arrowhead matrices: Serinv (Maillou et al. 2025, IEEE CLUSTER, doi:10.1109/cluster59342.2025.11186484). DALIA uses it (Gaedke-Merzhäuser et al. 2025, arXiv:2507.06938, 496 GH200), and it stems from INLA-DIST (SIAM J. Sci. Comput. 2024, doi:10.1137/23m1561531).
- NVIDIA cuDSS (https://docs.nvidia.com/cuda/cudss/latest/index.html): sparse direct solver on GPU, positive definite matrices, multi-GPU and multi-node. The documentation mentions neither a selected inverse nor a log-determinant.
- jax.experimental.sparse.linalg.spsolve: sparse solver by QR, CUDA only, without vmap (https://docs.jax.dev/en/latest/_autosummary/jax.experimental.sparse.linalg.spsolve.html).

### Position of remlax

Overlaps.
- pyreml shares most of the dense route of remlax: exact REML on GPU, an analytical gradient instead of autodiff through the Cholesky factor, Woodbury, a Kronecker product never formed, fa and multivariate models.
- kmmr shares the MME route: AI-REML, supernodal Cholesky, Takahashi partial inverse, GPU for the dense part.
- The Bayesian route overlaps with the generic use of NumPyro, BlackJAX or Stan.

Features unique to remlax in the sources consulted.
- A single engine combines the dense V form, the sparse MME form and AI-REML. pyreml has neither sparse matrices nor AI. kmmr has no V form, which it cites as future work.
- Automatic detection of global terms and blocks, with V = A + UU' by Woodbury and blocks processed in batches (vmap), has no equivalent in any tool found. The Woodbury form of pyreml factorises a dense C of size q.
- Fitting heavy multivariate models on a 10 GB MIG slice is reported nowhere else.
- The AI matrix provides standard errors and serves the covariance tests. pyreml gives none for the variances.

Capabilities present elsewhere and absent from remlax, according to the description received.
- Multi-GPU and multi-node computation: DISSECT, REACTA, SAIGE-GPU, ExaGeoStat, DALIA, Wang et al. 2019. jaxmg offers a distributed Cholesky with a log-determinant in JAX.
- The matrix-free route, with PCG and stochastic traces: BOLT-REML, SAIGE, SLEMM, MC-ss-GREML, GPyTorch, GPU_REML. Only this route goes beyond a few hundred thousand individuals.
- The mixed FP32 then FP64 precision of pyreml, useful on consumer-grade GPUs (L40S).
- Genotypes compressed to 2 bits and streamed to the GPU: miraculix, SAIGE-GPU, Quickdraws, GPU_REML.

### Unverified points

- The full text of kmmr beyond 100 000 characters was not read. The authors do not give the number of observations in their benchmarks.
- The gap between 261.65 s (pyreml documentation, version 0.5.0) and 279 s (preprint) for the case with 10⁶ observations. The two configurations may differ.
- The estimation mode and the code of TorchLIMIX.
- A step 2 of SAIGE-GPU on AMD GPUs, suggested only by a branch name.
- The maintenance of DISSECT, the licences of REACTA and ExaGeoStat, and the hardware of Wang et al. 2019.
- The absence of REML in the GPyTorch code (verified only in the documentation).
- Any GPU content in the WCGALP 2022 proceedings on BLUPF90, MiX99 and JWAS: PDFs blocked by Cloudflare. EAAP, ASAS-ADSA and Interbull were not covered in full text.
- The priority claimed by pyreml ("first library performing REML estimations on GPU") is contradicted by REACTA (2014) for a specialised REML. For a general-purpose REML, priority between pyreml (June 2026) and kmmr (repository of March 2026, preprint of October 2026) remains open.

### Sources

Articles and preprints (verified DOIs)
- Abdulah S et al. 2018. ExaGeoStat. IEEE TPDS 29(12). doi:10.1109/TPDS.2018.2850749
- Bermann M et al. 2025. GSE 57:44. doi:10.1186/s12711-025-00991-1
- Bermann M et al. 2025. GSE 57:63. doi:10.1186/s12711-025-01006-9
- Border R, Becker S. 2019. BMC Bioinformatics 20:411. doi:10.1186/s12859-019-2978-z
- Cabezas A et al. 2024. BlackJAX. arXiv:2402.10797
- Canela-Xandri O et al. 2015. Nat Commun 6:10162. doi:10.1038/ncomms10162
- Cebamanos L et al. 2014. Bioinformatics. doi:10.1093/bioinformatics/btt754
- Češnovar R et al. 2019. GPU-based parallel computation support for Stan. arXiv:1907.01063
- Charlier B et al. 2021. JMLR 22(74):1-6. https://www.jmlr.org/papers/v22/20-275.html
- Cheng H et al. 2023. SLEMM. Bioinformatics. doi:10.1093/bioinformatics/btad127
- Coffey M et al. 2011. Interbull Bull. 44:114-116. https://journal.interbull.org/index.php/ib/article/view/1694
- DeVore et al. 2014. IEEE BIBM. doi:10.1109/bibm.2014.6999186
- Dillon JV et al. 2017. TensorFlow Distributions. arXiv:1711.10604
- Fernando RL et al. 2016. GSE 48:96. doi:10.1186/s12711-016-0273-2
- Freudenberg A et al. 2023. Front. Genet. 14:1220408. doi:10.3389/fgene.2023.1220408
- Gaedke-Merzhäuser L et al. 2022. Stat Comput. doi:10.1007/s11222-022-10192-1
- Gaedke-Merzhäuser L et al. 2024. SIAM J. Sci. Comput. doi:10.1137/23m1561531
- Gaedke-Merzhäuser L et al. 2025. DALIA. arXiv:2507.06938, doi:10.48550/arXiv.2507.06938
- Gardner JR et al. 2018. GPyTorch. NeurIPS. arXiv:1809.11165
- Geoga CJ, Anitescu M, Stein ML. 2020. JCGS. doi:10.1080/10618600.2019.1652616
- Gilmour AR. 2019. JABG 136:262. doi:10.1111/jbg.12398
- Horn B, Nikoloski Z, Lippert C. 2026. TorchLIMIX. bioRxiv. doi:10.64898/2026.04.29.721342
- Hutchinson MF. 1989. Commun Stat Simul Comput. doi:10.1080/03610918908812806
- Jacquelin M, Lin L, Yang C. 2016. PSelInv. ACM TOMS. doi:10.1145/2786977
- Jiang L et al. 2019. fastGWA. Nat Genet. doi:10.1038/s41588-019-0530-8
- Lee SH, van der Werf JHJ. 2016. MTG2. Bioinformatics. doi:10.1093/bioinformatics/btw012
- Lewis AGM et al. 2022. PNAS. doi:10.1073/pnas.2122762119
- Lippert C et al. 2011. FaST-LMM. Nat Methods. doi:10.1038/nmeth.1681
- Loh PR et al. 2015. BOLT-LMM. Nat Genet. doi:10.1038/ng.3190
- Loh PR et al. 2015. BOLT-REML. Nat Genet 47:1385. doi:10.1038/ng.3431
- Loya H et al. 2025. Quickdraws. Nat Genet. doi:10.1038/s41588-024-02044-7
- Maillou V et al. 2025. Serinv. IEEE CLUSTER. doi:10.1109/cluster59342.2025.11186484
- Marchal A, Trottier C, Bry X, Raimondi D. 2026. pyreml. bioRxiv. doi:10.64898/2026.10.20.757556
- Matilainen K et al. 2012. JABG. doi:10.1111/j.1439-0388.2012.01000.x
- Matilainen K et al. 2013. PLoS One 8:e80821. doi:10.1371/journal.pone.0080821
- Matthews AGG et al. 2017. GPflow. JMLR 18(40)
- Mbatchou J et al. 2021. REGENIE. Nat Genet. doi:10.1038/s41588-021-00870-7
- Melsen KAC, van Eeuwijk FA, Peeters CFW. 2026. kmmr. arXiv:2610.16671, doi:10.48550/arXiv.2610.16671
- Meyer K, Tier B. 2013. Proc. AAABG 20:270-273. https://hdl.handle.net/1959.11/14323
- Pazokitoroudi A et al. 2020. RHE-mc. Nat Commun. doi:10.1038/s41467-020-17576-9
- Phan D, Pradhan N, Jankowiak M. 2019. NumPyro. arXiv:1912.11554
- Pinder T, Dodd D. 2022. GPJax. JOSS. doi:10.21105/joss.04455
- Rennich SC, Stosic D, Davis TA. 2016. Parallel Computing. doi:10.1016/j.parco.2016.06.004
- Rodriguez A et al. 2026. SAIGE-GPU. Bioinformatics 42(3):btag032. doi:10.1093/bioinformatics/btag032
- Simchoni G, Rosset S. 2023. JMLR 24(156):1-57. https://www.jmlr.org/papers/v24/22-0501.html
- Strandén I et al. 2024. GSE 56:73. doi:10.1186/s12711-024-00939-x
- Taylor-Weiner A et al. 2019. tensorQTL. Genome Biol. doi:10.1186/s13059-019-1836-7
- Ubaru S, Chen J, Saad Y. 2017. SIAM J Matrix Anal Appl. doi:10.1137/16M1104974
- Vandenplas J et al. 2022. MiXBLUP 3.0. Proc. WCGALP. doi:10.3920/978-90-8686-940-4_358
- Wang KA et al. 2019. Exact GPs on a million data points. NeurIPS. arXiv:1903.08114
- Wu Y, Sankararaman S. 2018. RHE-reg. Bioinformatics. doi:10.1093/bioinformatics/bty253
- Yang J et al. 2011. GCTA. AJHG. doi:10.1016/j.ajhg.2010.11.011
- Yin L et al. 2023. HIBLUP. NAR 51:3501. doi:10.1093/nar/gkad074
- Younis OG et al. 2023. ChromaX. Bioinformatics. doi:10.1093/bioinformatics/btad691
- Yue K et al. 2021. REHE. Genet Epidemiol 45(8). doi:10.1002/gepi.22432
- Zhang Z et al. 2025. jaxQTL. medRxiv. doi:10.1101/2025.01.18.25320755
- Zhou H, Hu L, Zhou J, Lange K. 2019. MM algorithms. JCGS 28(2). doi:10.1080/10618600.2018.1529601
- Zhou W et al. 2018. SAIGE. Nat Genet. doi:10.1038/s41588-018-0184-y
- Zhou X, Stephens M. 2012. GEMMA. Nat Genet. doi:10.1038/ng.2310
- Zhou X, Stephens M. 2014. mvLMM. Nat Methods. doi:10.1038/nmeth.2848

Repositories, packages and documentation
- pyreml: https://github.com/ai-for-genome-interpretation-lab/pyreml ; https://pypi.org/project/pyreml/ ; https://ai-for-genome-interpretation-lab.github.io/pyreml/ ; Zenodo doi:10.5281/zenodo.21394319
- kmmr: https://github.com/KillianMelsen/kmmr ; Zenodo doi:10.5281/zenodo.23398636
- SAIGE-GPU: https://github.com/saigegit/SAIGE (branch SAIGE-GPU-1.3.3) ; Zenodo doi:10.5281/zenodo.17642591
- GPU_REML: https://github.com/Asiandier/SPLIT
- Faster-LMM-D: https://github.com/genetics-statistics/faster_lmm_d
- pylmm: https://github.com/nickFurlotte/pylmm
- FaST-LMM: https://github.com/fastlmm/FaST-LMM
- GEMMA: https://github.com/genetics-statistics/GEMMA
- BOLT-LMM, manual v2.5: https://alkesgroup.broadinstitute.org/BOLT-LMM/BOLT-LMM_manual.html
- GCTA: https://github.com/jianyangqt/gcta
- REGENIE: https://github.com/rgcgithub/regenie ; Glow: https://github.com/projectglow/glow
- Quickdraws: https://github.com/PalamaraLab/quickdraws
- tensorQTL: https://github.com/broadinstitute/tensorqtl ; jaxQTL: https://github.com/mancusolab/jaxqtl
- Hail: https://github.com/hail-is/hail
- SL_REML: https://github.com/rborder/SL_REML ; SLEMM: https://github.com/jiang18/slemm
- GPyTorch: https://github.com/cornellius-gp/gpytorch ; settings: https://docs.gpytorch.ai/en/stable/settings.html
- KeOps: https://github.com/getkeops/keops
- GPflow: https://github.com/GPflow/GPflow ; tinygp: https://github.com/dfm/tinygp
- TensorFlow Probability, LMM example: https://www.tensorflow.org/probability/examples/Linear_Mixed_Effects_Models
- NumPyro: https://github.com/pyro-ppl/numpyro ; BlackJAX: https://github.com/blackjax-devs/blackjax
- Stan OpenCL: https://mc-stan.org/docs/cmdstan-guide/parallelization.html
- jaxmg: https://github.com/flatironinstitute/jaxmg
- jaxgam: https://github.com/scatt67/jaxgam
- cuDSS: https://docs.nvidia.com/cuda/cudss/latest/index.html
- JAX spsolve: https://docs.jax.dev/en/latest/_autosummary/jax.experimental.sparse.linalg.spsolve.html
- miraculix: https://github.com/alexfreudenberg/miraculix
- BOLT (Theta Solutions), WCGALP 2018: https://www.icar.org/Documents/Auckland-2018/WCGALP/1415%20Dorian%20Garrick.pdf
- JWAS: https://github.com/reworkhow/JWAS.jl

## Appendix C. Prior art of the remlax mechanisms

Most remlax building blocks have clear prior art in quantitative genetics. The REML gradient with respect to V, AI-REML, sparse factorisation of the mixed model equations and selected inversion were published between 1977 and 2015. Nested dissection and the likelihood taken from the Cholesky factor of a bordered matrix date from the same period. They are also implemented in ASReml, WOMBAT, DFREML/MTDFREML and the BLUPF90/YAMS family. Automatic detection of the terms to handle at low rank and of independent blocks of observations, through a cost model and a greedy search, has no prior art found. The active set on the positive semi-definite cone on the variance scale, followed by a trust-region finisher, has none either in the REML software examined. Joint REML estimation of a continuous IGE kernel range with the variance components also has no antecedent found. GPU execution with automatic differentiation has contemporary prior art (pyreml, 2026), and older prior art for automatic differentiation itself (Meyer and Smith 1996).

### Summary table

| remlax mechanism | Verdict | Main references |
|---|---|---|
| 1. REML -2logL from a dense V, gradient P − Pyy'P | Clear prior art for the formula and the dense V form | [Harville 1977](https://doi.org/10.1080/01621459.1977.10480998), [Searle et al. 1992](https://doi.org/10.1002/9780470316856), [Covarrubias-Pazaran 2016](https://doi.org/10.1371/journal.pone.0156744) |
| 1b. Gradient with respect to V coded as an adjoint rule (custom_vjp), the rest by automatic differentiation | Partial prior art | [Smith 1995](https://doi.org/10.1080/10618600.1995.10474671), [Giles 2008](https://doi.org/10.1007/978-3-540-68942-3_4), [Kristensen et al. 2016](https://doi.org/10.18637/jss.v070.i05), pyreml 0.5.3 ([Marchal and Raimondi 2026](https://doi.org/10.5281/zenodo.21394319)) |
| 2a. Sparse mixed model equations, supernodal Cholesky | Clear prior art | [Masuda et al. 2014](https://doi.org/10.1111/jbg.12058), [Masuda et al. 2015](https://doi.org/10.2527/jas.2015-9395), WOMBAT manual 2017, ASReml 4.2 manual |
| 2b. Selected inversion (Takahashi) for the gradient and AI | Clear prior art | [Erisman and Tinney 1975](https://doi.org/10.1145/360680.360704), [Misztal and Perez-Enciso 1993](https://doi.org/10.3168/jds.S0022-0302%2893%2977478-0), [Masuda et al. 2014](https://doi.org/10.1111/jbg.12058), [Lin et al. 2011](https://doi.org/10.1145/1916461.1916464) |
| 2c. Failure of AMD on a locally dense C, gain from nested dissection or a "blocks then separator" order | Partial prior art (mechanism and use known, IGE finding unpublished) | [George 1973](https://doi.org/10.1137/0710032), [Karypis and Kumar 1998](https://doi.org/10.1137/S1064827595287997), [Amestoy et al. 2004](https://doi.org/10.1145/1024074.1024081), [Masuda et al. 2014](https://doi.org/10.1111/jbg.12058), WOMBAT and ASReml manuals |
| 3a. V = A + UU', Woodbury, block-diagonal A | Clear prior art for each element | [Lindstrom and Bates 1988](https://doi.org/10.1080/01621459.1988.10478693), [Mishchenko et al. 2007](https://doi.org/10.48550/arXiv.0709.0625), [Lippert et al. 2011](https://doi.org/10.1038/nmeth.1681), [Strandén and Garrick 2009](https://doi.org/10.3168/jds.2008-1929) |
| 3b. Cholesky factor of a bordered matrix [U X y]'A⁻¹[U X y] + diag(I,0,0) giving log|V|, log|X'V⁻¹X|, y'Py | Clear prior art for the algebra; partial when A contains absorbed random terms | [Smith and Graser 1986](https://doi.org/10.3168/jds.S0022-0302%2886%2980516-1), [Graser et al. 1987](https://doi.org/10.2527/jas1987.6451362x), [Meyer 1989](https://doi.org/10.1186/1297-9686-21-3-317), [Bates and DebRoy 2004](https://doi.org/10.1016/j.jmva.2004.04.013), [Bates et al. 2025](https://doi.org/10.48550/arXiv.2505.11674) |
| 3c. Automatic detection of global terms and independent blocks (cost model, greedy search) | No prior art found; only related heuristics | ASReml 4.2 manual (!DENSE, !GDENSE), WOMBAT manual (--dense), [Kristensen et al. 2016](https://doi.org/10.18637/jss.v070.i05), [Bates et al. 2025](https://doi.org/10.48550/arXiv.2505.11674), pyreml 0.5.3 |
| 4a. AI-REML on the variance scale | Clear prior art | [Gilmour et al. 1995](https://doi.org/10.2307/2533274), [Johnson and Thompson 1995](https://doi.org/10.3168/jds.S0022-0302%2895%2976654-1) |
| 4b. Active set, variance floor | Partial prior art (scalar bounds) | [Mishchenko et al. 2008](https://doi.org/10.3233/JCM-2008-81-203), ASReml 4.2 manual (codes P and B), WOMBAT manual (--pivot) |
| 4c. Active set on the positive semi-definite cone | No prior art found in the REML software examined | ASReml replaces the AI update by an EM update; WOMBAT reparameterises by Cholesky ([Meyer and Smith 1996](https://doi.org/10.1186/1297-9686-28-1-23)) |
| 4d. Levenberg-Marquardt damping | Clear prior art | WOMBAT manual (--modaim2), [Wolfinger et al. 1994](https://doi.org/10.1137/0915079) |
| 4e. L-BFGS-B bursts, Moré-Sorensen finisher with SR1 | Classical components; no use in REML found (except BFGS in Mishchenko) | [Byrd et al. 1995](https://doi.org/10.1137/0916069), [Moré and Sorensen 1983](https://doi.org/10.1137/0904038), [Conn et al. 1991](https://doi.org/10.1007/BF01594934), [Mishchenko et al. 2008](https://doi.org/10.3233/JCM-2008-81-203) |
| 4f. EM or PX-EM start | Clear prior art | [Liu et al. 1998](https://doi.org/10.1093/biomet/85.4.755), [Foulley and van Dyk 2000](https://doi.org/10.1186/1297-9686-32-2-143), WOMBAT and ASReml manuals |
| 5. REML on GPU with automatic differentiation | Partial and contemporary prior art | [Meyer and Smith 1996](https://doi.org/10.1186/1297-9686-28-1-23), [Gardner et al. 2018](https://doi.org/10.48550/arXiv.1809.11165), pyreml 0.5.3, [Horn et al. 2026](https://doi.org/10.64898/2026.04.29.721342) |
| 6a. Neighbourhood IGE model, dilution | Clear prior art | [Griffing 1967](https://doi.org/10.1071/BI9670127), [Muir 2005](https://doi.org/10.1534/genetics.104.035956), [Bijma et al. 2007](https://doi.org/10.1534/genetics.106.062711), [Bijma 2010](https://doi.org/10.1534/genetics.110.120485), WOMBAT manual (SOCIAL) |
| 6b. Kernel with a discrete order or weighted by distance | Clear prior art | [Cappa and Cantet 2008](https://doi.org/10.1515/sg-2008-0008), [Costa e Silva et al. 2017](https://doi.org/10.1038/hdy.2017.25), [Sato et al. 2021](https://doi.org/10.1038/s41437-020-00401-w) |
| 6c. Continuous range estimated by REML with the variances | No prior art found | likelihood profile of the dilution in WOMBAT; scale sweep by PVE in [Sato et al. 2021](https://doi.org/10.1038/s41437-020-00401-w) |

### Dense V form and adjoint gradient

The derivative of the REML log-likelihood is -½ tr(P ∂V/∂θ) + ½ y'P (∂V/∂θ) Py. This form appears in the review by [Harville 1977](https://doi.org/10.1080/01621459.1977.10480998) and in the textbook by [Searle et al. 1992](https://doi.org/10.1002/9780470316856). The derivative of -2logL with respect to V is therefore P − Pyy'P. remlax brings nothing new to the formula. Direct REML evaluation on a dense V is just as old. The sommer package uses it in its historical function, with Newton-Raphson and AI ([Covarrubias-Pazaran 2016](https://doi.org/10.1371/journal.pone.0156744)). Human genomics tools combine it with an eigendecomposition ([Kang et al. 2008](https://doi.org/10.1534/genetics.107.080101); [Yang et al. 2011](https://doi.org/10.1016/j.ajhg.2010.11.011)).

Automatic differentiation of the REML likelihood is also old. [Smith 1995](https://doi.org/10.1080/10618600.1995.10474671) differentiated the Cholesky algorithm in reverse mode. Murray describes this method as the manual application of automatic differentiation ideas to the Cholesky algorithm ([Murray 2016](https://doi.org/10.48550/arXiv.1602.07527)). [Meyer and Smith 1996](https://doi.org/10.1186/1297-9686-28-1-23) applied this differentiation to animal models. Their abstract reports first and second derivatives without inverting large matrices. The adjoint rules of the log-determinant and of the inverse are collected by [Giles 2008](https://doi.org/10.1007/978-3-540-68942-3_4). TMB computes the gradient of the log-determinant of a sparse Hessian by a reverse pass and the inverse subset algorithm ([Kristensen et al. 2016](https://doi.org/10.18637/jss.v070.i05), section 4.1). The pyreml 0.5.3 package exposes an "analytic_backward" option in a PyTorch implementation (source code examined locally, [Marchal and Raimondi 2026](https://doi.org/10.5281/zenodo.21394319)). We did not verify the exact content of this rule.

Verdict. The formula and the dense V form have clear prior art. An adjoint rule inserted into an automatic differentiation graph has partial prior art. The contribution of remlax lies in the JAX engineering and in CPU/GPU portability.

### Sparse mixed model equations

#### Absorption and the sparse/dense partition

The partition of the terms into a dense block and a sparse block is documented in ASReml. The 4.2 manual states that the partition is made at the !r point, unless set by !DENSE (section 6.10.1). By default, the dense block receives the fixed effects listed before the random effects, up to 800 equations by default. The !GDENSE option places the first GRM term in the dense equations. The user sets this partition. ASReml first absorbs the sparse equations, in an order it chooses, then the dense equations in the declared order (section 6.10.2). WOMBAT offers a --dense option when a relationship matrix is dense, for example a GRM (2017 manual, section 5.3.4). Its --s1step option groups the genotyped individuals into a dense block handled separately. All these partitions live in C. None of them moves a term to a V form.

#### Ordering

Nested dissection goes back to [George 1973](https://doi.org/10.1137/0710032). METIS provides a multilevel version of it ([Karypis and Kumar 1998](https://doi.org/10.1137/S1064827595287997)). The approximate minimum degree is described by [Amestoy et al. 1996](https://doi.org/10.1137/S0895479894278952). [Liu 1985](https://doi.org/10.1145/214392.214398) proposed the multiple minimum degree.

In animal breeding, these orderings have already been compared. YAMS uses AMD or multilevel nested dissection ([Masuda et al. 2014](https://doi.org/10.1111/jbg.12058)). Its authors report a factorisation up to 13 times faster than FSPAK "when an appropriate ordering strategy was applied". WOMBAT chooses the ordering according to the size of the problem (2017 manual, appendix A.1). Multiple minimum degree serves up to 25 000 equations, AMD up to 50 000, and METIS beyond. This manual recommends trying several orderings and offers --metisall, which runs through 48 combinations of METIS settings. It cites a study by Meyer (2005, Proc. AAABG 16:282–285, no DOI, content unverified). ASReml exposes !EQORDER (4.2 manual, table 5.5). Its manual describes a case close to the remlax finding. For a model with a dense IBD matrix and a sparse relationship matrix, the order declared by the user (!EQORDER -1) runs faster than the automatic order.

The mechanism of the AMD failure is documented. The SuiteSparse AMD removes rows with more than 10√n entries and places them at the end of the ordering (amd.h header; [Amestoy et al. 2004](https://doi.org/10.1145/1024074.1024081)). A locally dense C, with many rows just below this threshold, therefore degrades the ordering. WOMBAT also passes METIS a factor that flags the dense rows.

Verdict. The use of METIS for C in quantitative genetics has clear prior art. The "blocks then separator" order corresponds to the first level of a nested dissection. The quantified finding on a multivariate IGE model was not found in the literature. It deserves to be reported as an empirical result, without a methodological claim.

#### Selected inversion

The Takahashi recurrence is attributed to Takahashi, Fagan and Chin (1973, 8th PICA Conference, Minneapolis, pp. 63–69). This communication has no DOI and was not consulted. [Erisman and Tinney 1975](https://doi.org/10.1145/360680.360704) give the version published in a journal. [Misztal and Perez-Enciso 1993](https://doi.org/10.3168/jds.S0022-0302%2893%2977478-0) introduced it into EM-REML by sparse inversion. The exact algorithm they used was not checked in the text. [Johnson and Thompson 1995](https://doi.org/10.3168/jds.S0022-0302%2895%2976654-1) combine it with AI-REML for univariate animal models, according to their title. YAMS uses a left-looking supernodal factorisation and an inverse multifrontal inversion ([Masuda et al. 2014](https://doi.org/10.1111/jbg.12058)). Handling the traces with a hash table speeds up AI-REML by up to 67 times ([Masuda et al. 2015](https://doi.org/10.2527/jas.2015-9395)). WOMBAT uses a supernodal approach by default and cites Erisman and Tinney (2017 manual, section 5.3.5). ASReml 4.2 added a preliminary symbolic absorption and a "nodal absorption" (4.2 manual, section 11.5). Outside genetics, SelInv is the reference supernodal algorithm ([Lin et al. 2011](https://doi.org/10.1145/1916461.1916464)).

Verdict. Clear prior art. remlax reimplements a standard technique.

### Global terms, Woodbury and the bordered matrix

#### Woodbury and low rank

Henderson's equations are themselves a Woodbury form. [Strandén and Garrick 2009](https://doi.org/10.3168/jds.2008-1929) derive the equivalence between SNP-BLUP and GBLUP through this identity. [Mishchenko et al. 2007](https://doi.org/10.48550/arXiv.0709.0625) apply it to AI-REML with a low-rank representation, exact or truncated, of the IBD matrix. FaST-LMM uses it when the relationship matrix has rank lower than n ([Lippert et al. 2011](https://doi.org/10.1038/nmeth.1681)). Handling V as block diagonal by subject is the basic assumption of [Lindstrom and Bates 1988](https://doi.org/10.1080/01621459.1988.10478693). SAS PROC MIXED declares these blocks through SUBJECT= and applies the algorithms of [Wolfinger et al. 1994](https://doi.org/10.1137/0915079). The "block diagonal plus low rank" case also appears in single-step genomics ([Liu et al. 2014](https://doi.org/10.3168/jds.2014-7924); [Taskinen et al. 2017](https://doi.org/10.1186/s12711-017-0310-9)).

#### Eigendecomposition

Rotation by the eigenvectors of the relationship matrix makes V diagonal when a single random term is added to the residual ([Kang et al. 2008](https://doi.org/10.1534/genetics.107.080101); [Lippert et al. 2011](https://doi.org/10.1038/nmeth.1681); [Zhou and Stephens 2012](https://doi.org/10.1038/ng.2310)). GCTA also uses it ([Yang et al. 2011](https://doi.org/10.1016/j.ajhg.2010.11.011)). The multivariate version relies on a joint canonical transformation ([Meyer 1985](https://doi.org/10.2307/2530651); [Zhou and Stephens 2014](https://doi.org/10.1038/nmeth.2848); [Lee and van der Werf 2016](https://doi.org/10.1093/bioinformatics/btw012)). This trick does not extend to several terms that cannot be diagonalised simultaneously. IGE models contain several such terms. remlax differs on this point.

#### Bordered likelihood

[Smith and Graser 1986](https://doi.org/10.3168/jds.S0022-0302%2886%2980516-1) compute log|C| and y'Py by Gaussian elimination on the matrix of the equations bordered by the row [r' y'R⁻¹y]. The last pivot gives y'Py = y'R⁻¹y − s'r. This description comes from the MTDFREML workshop notes (Boldman et al., 2002, USDA-ARS, no DOI). The original text was not consulted. [Graser et al. 1987](https://doi.org/10.2527/jas1987.6451362x) and [Meyer 1989](https://doi.org/10.1186/1297-9686-21-3-317) made it the basis of derivative-free REML. [Boldman and Van Vleck 1991](https://doi.org/10.3168/jds.S0022-0302%2891%2978629-3) sped it up with a sparse factorisation (SPARSPAK). In statistics, [Bates and DebRoy 2004](https://doi.org/10.1016/j.jmva.2004.04.013) and then lme4 ([Bates et al. 2015](https://doi.org/10.18637/jss.v067.i01)) factorise the augmented matrix of [ZΛ X y]. [Bates et al. 2025](https://doi.org/10.48550/arXiv.2505.11674) show that 2log|R_ZZ| and r²_yy can be read directly on the diagonal of the factor, without solving for β or for the modes. This construction assumes R = σ²I.

The bordered matrix of remlax, [U X y]'A⁻¹[U X y] + diag(I,0,0), is algebraically the augmented Henderson equation, with a matrix A in place of R. The case of a block-diagonal R already exists in the multivariate analyses of DFREML. Any novelty lies in the content of A. remlax absorbs local random terms into it in V form and leaves the global terms in U. This V/C hybrid was not found in the sources consulted.

#### Automatic structure detection

No software examined chooses on its own, through a cost model, the terms to handle at low rank and the independent blocks of observations. The related mechanisms are fixed heuristics or choices left to the user. ASReml partitions at the !r point and leaves !DENSE and !GDENSE to the user (4.2 manual). WOMBAT chooses an ordering according to the number of equations and leaves --dense to the user (2017 manual). MixedModels.jl puts first the term with the largest number of levels, because the diagonal structure of block [1,1] is preserved in the factor ([Bates et al. 2025](https://doi.org/10.48550/arXiv.2505.11674), section 3.3). TMB detects the sparsity pattern of the Hessian automatically ([Kristensen et al. 2016](https://doi.org/10.18637/jss.v070.i05), section 4.2). This pattern serves the sparse factorisation. The choice between the V and C forms stays with the modeller. pyreml 0.5.3 uses the Woodbury route when the number of random effects is below n and the residual allows it (MixedModel.py, code examined locally). This rule bears on a single dimension and detects no blocks.

Verdict. The algebraic building blocks have clear prior art. Automatic detection by a cost model and a greedy search seems new. This conclusion rests on a non-exhaustive search. The internal engines of ASReml are closed and may contain undocumented heuristics.

### AI-REML optimisation and robustness

AI-REML was published by [Gilmour et al. 1995](https://doi.org/10.2307/2533274) and, for sparse animal models, by [Johnson and Thompson 1995](https://doi.org/10.3168/jds.S0022-0302%2895%2976654-1). ASReml applies it on the component scale. Newton-Raphson and EM for mixed models are described by [Lindstrom and Bates 1988](https://doi.org/10.1080/01621459.1988.10478693) and [Jennrich and Schluchter 1986](https://doi.org/10.2307/2530695). [Thompson 2008](https://doi.org/10.1098/rspb.2007.1417) traces this history.

Reparameterisation by the Cholesky factor, with the log of the diagonal elements, removes the constraints ([Groeneveld 1994](https://doi.org/10.1186/1297-9686-26-6-537); [Pinheiro and Bates 1996](https://doi.org/10.1007/BF00140873); [Meyer and Smith 1996](https://doi.org/10.1186/1297-9686-28-1-23)). WOMBAT uses it by default and bounds the pivots at 0.001 (option --pivot, 2017 manual, section 5.2.5). ASReml works on the variance scale. Its code P replaces a negative variance by a small positive value, then fixes it after repeated attempts (code B, 4.2 manual, section 7.7.4). For a US matrix, !GP triggers an EM update in place of the AI update when the latter leaves the cone. ASReml also fixes the whole matrix when it becomes nearly singular (4.2 manual, !EMFLAG). [Mishchenko et al. 2008](https://doi.org/10.3233/JCM-2008-81-203) show that standard AI-REML fails when the optimum lies on a bound. They propose an active-set method and an inverse BFGS update, for scalar non-negativity constraints.

Levenberg-Marquardt-type damping is documented. WOMBAT adds τI to the AI matrix, with τ = max(0, δ − λmin), as the default modification (--modaim2, appendix A.5). It also offers eigenvalue clamping and two modified Cholesky factorisations. SAS PROC MIXED minimises -2logL by a ridge-stabilised Newton-Raphson ([Wolfinger et al. 1994](https://doi.org/10.1137/0915079); SAS documentation). Starting with a few PX-EM iterations before AI is the WOMBAT default for more than 18 components (--pxai). ASReml offers local EM or PX-EM iterations. PX-EM comes from [Liu et al. 1998](https://doi.org/10.1093/biomet/85.4.755) and was adapted to Henderson's model by [Foulley and van Dyk 2000](https://doi.org/10.1186/1297-9686-32-2-143). WOMBAT keeps Powell and the simplex to check convergence in difficult cases.

The other components are classical in numerical optimisation. L-BFGS-B comes from [Byrd et al. 1995](https://doi.org/10.1137/0916069). The trust-region step by eigendecomposition follows [Moré and Sorensen 1983](https://doi.org/10.1137/0904038). The convergence of SR1 is established by [Conn et al. 1991](https://doi.org/10.1007/BF01594934). [Nocedal and Wright 2006](https://doi.org/10.1007/978-0-387-40065-5) collect these methods. Projection onto the positive semi-definite cone follows [Higham 1988](https://doi.org/10.1016/0024-3795%2888%2990223-6). No REML software examined combines an active set on the semi-definite cone, L-BFGS-B bursts and a Moré-Sorensen finisher with SR1.

Verdict. AI, damping, the floor and the EM start have clear prior art. The active set has partial prior art, limited to scalar bounds. Its extension to the semi-definite cone on the variance scale, and the complete assembly, have no prior art found.

### REML on GPU and automatic differentiation

Automatic differentiation applied to REML dates from [Meyer and Smith 1996](https://doi.org/10.1186/1297-9686-28-1-23), on CPU. TMB generalises automatic differentiation to latent effect models, with CHOLMOD and the inverse subset ([Kristensen et al. 2016](https://doi.org/10.18637/jss.v070.i05)). The parallelism described in this article goes through OpenMP and BLAS, without GPU. Outside genetics, GPyTorch computes the marginal likelihood of a Gaussian process and its derivatives on GPU, by batched conjugate gradients ([Gardner et al. 2018](https://doi.org/10.48550/arXiv.1809.11165)). GPJax does the same in JAX ([Pinder and Dodd 2022](https://doi.org/10.21105/joss.04455)). GPyTorch and GPJax maximise the full likelihood and do not handle REML.

In genetics, direct competitors appeared in 2026. pyreml fits mixed models by direct differentiation of REML in PyTorch, with GPU acceleration ([Marchal and Raimondi 2026](https://doi.org/10.5281/zenodo.21394319)). The version examined, 0.5.3, uses L-BFGS with Adam as a fallback (Optimizer.py). TorchLIMIX ports the LIMIX multivariate pipeline to GPU, with gains of up to two orders of magnitude ([Horn et al. 2026](https://doi.org/10.64898/2026.04.29.721342), preprint). SAIGE-GPU speeds up the SAIGE GLMM fit by a factor of 5 ([Rodriguez et al. 2026](https://doi.org/10.1093/bioinformatics/btag032)). The torch-openreml package announces an AI-REML engine on CPU and GPU. Its PyPI page could not be read, so it remains unverified. No reference found combines sparse mixed model equations, supernodal Cholesky and selected inversion in an automatic differentiation framework on GPU. Recent work by the Misztal team stays on CPU, with Monte Carlo REML and preconditioned conjugate gradients ([Bermann et al. 2025](https://doi.org/10.1186/s12711-025-01006-9)).

Verdict. Automatic differentiation of REML has clear prior art. REML on GPU in dense V form has partial and contemporary prior art (pyreml). The sparse C form on GPU, with validated CPU/GPU parity, has no prior art found.

### Neighbourhood IGE models

The IGE framework is established. [Griffing 1967](https://doi.org/10.1071/BI9670127) sets out selection at the group level. [Moore et al. 1997](https://doi.org/10.1111/j.1558-5646.1997.tb01458.x) and [Wolf et al. 1998](https://doi.org/10.1016/S0169-5347%2897%2901233-0) formalise indirect genetic effects. [Muir 2005](https://doi.org/10.1534/genetics.104.035956) integrates them into BLUP for trees and animals. [Bijma et al. 2007a](https://doi.org/10.1534/genetics.106.062711) and [Bijma et al. 2007b](https://doi.org/10.1534/genetics.106.062729) define the total breeding value and its estimation. The dependence on group size is treated by [Hadfield and Wilson 2007](https://doi.org/10.1534/genetics.107.075622). [Bijma 2010a](https://doi.org/10.1534/genetics.110.120485) introduces the dilution factor d. Optimal designs are studied by [Bijma 2010b](https://doi.org/10.1534/genetics.110.120493). [Arango et al. 2005](https://doi.org/10.2527/2005.8361241x) estimate competitive effects in pigs.

The spatial neighbourhood of plants has long been modelled. Neighbour adjustments in field trials go back to [Draper and Guttman 1980](https://doi.org/10.2307/2986297), [Kempton 1982](https://doi.org/10.1017/S0021859600054381) and [Besag and Kempton 1986](https://doi.org/10.2307/2531047). [Stringer et al. 2011](https://doi.org/10.1007/s13253-010-0051-5) jointly model spatial variability and genetic competition between plots. In trees, [Cappa and Cantet 2008](https://doi.org/10.1515/sg-2008-0008) weight the neighbours by "competition intensities", inverses of the distance and of the number of neighbours. They estimate the model by Gibbs sampling and by EM-REML. [Cappa et al. 2015](https://doi.org/10.1007/s11295-015-0917-3) add environmental heterogeneity. [Costa e Silva et al. 2013](https://doi.org/10.1111/nph.12035) contrast disease IGEs and competition IGEs in eucalyptus. In the same species, [Costa e Silva et al. 2017](https://doi.org/10.1038/hdy.2017.25) identify the most influential neighbour positions through a factor-analytic structure, then build a combined effect weighted by distance. The review by [Costa e Silva and Kerr 2013](https://doi.org/10.1007/s11295-012-0521-8) covers these models. [Sato et al. 2021](https://doi.org/10.1038/s41437-020-00401-w) estimate the effective range of neighbour effects at the peak of explained variance, over a grid of spatial scales. A more recent preprint mentions a decay with distance without fitting it ([Sato and Hamazaki 2026](https://doi.org/10.7554/eLife.111650.1)).

For crop mixtures, the literature mainly deals with mixing ability and genotypic interactions. [Knott and Mundt 1990](https://doi.org/10.1007/bf00210065) analyse mixing ability in wheat. [Forst et al. 2019](https://doi.org/10.1016/j.fcr.2019.107571) propose a framework for incomplete mixing designs. [Montazeaud et al. 2022](https://doi.org/10.1111/nph.17915) test allelic richness between genotypes in durum wheat. [Wuest and Niklaus 2018](https://doi.org/10.1038/s41559-018-0708-y) and [Wuest et al. 2021](https://doi.org/10.1038/s41559-021-01497-x) relate genetic diversity to the yield of mixtures. [Pélissier et al. 2023](https://doi.org/10.1371/journal.pbio.3002287) show that the neighbouring genotype changes disease susceptibility in wheat and rice.

The software used is varied. WOMBAT fits "social" effects with a fixed dilution factor (2017 manual, section 4.10.6.8). It writes the log-likelihood for each d and offers a quadratic approximation of the profile (--quapp). ASReml allows overlays of incidence matrices through and() (4.2 manual). Its use in the cited forestry IGE studies was not checked in the texts. Cappa and Cantet use their own Gibbs sampler. Sato and colleagues use rNeighborGWAS and RAINBOWR.

Verdict. The IGE model, the dilution, and the discrete or distance-weighted kernel have clear prior art. Joint REML estimation of a continuous range with the variance components has no prior art found. Existing approaches sweep a grid (Sato) or profile the likelihood (WOMBAT for d). The search on interspecific mixtures with within-species and between-species IGEs was not exhaustive.

### Uncertainties

Several primary texts were not read in full. The description of the bordered matrix of Smith and Graser 1986 comes from MTDFREML workshop notes. The results of Meyer 2005 on ordering are known only through their citation in the WOMBAT manual. The choice of ordering between AMD and nested dissection by matrix type, in Masuda et al. 2014, is known only from the abstract. The internal algorithms of ASReml are closed. An undocumented partition or constraint heuristic cannot be ruled out. The pyreml version examined is 0.5.3 from the local repository. Its exact behaviour (adjoint rule, Woodbury route) deserves a targeted reading before any published comparison. The torch-openreml page could not be loaded.

### Bibliography

Articles and chapters (DOIs verified on CrossRef or doi.org)

- Amestoy PR, Davis TA, Duff IS (1996). An approximate minimum degree ordering algorithm. SIAM J Matrix Anal Appl 17:886–905. https://doi.org/10.1137/S0895479894278952
- Amestoy PR, Davis TA, Duff IS (2004). Algorithm 837: AMD, an approximate minimum degree ordering algorithm. ACM TOMS 30:381–388. https://doi.org/10.1145/1024074.1024081
- Arango J, Misztal I, Tsuruta S, et al. (2005). Estimation of variance components including competitive effects of Large White growing gilts. J Anim Sci 83:1241–1246. https://doi.org/10.2527/2005.8361241x
- Bates D, DebRoy S (2004). Linear mixed models and penalized least squares. J Multivar Anal 91:1–17. https://doi.org/10.1016/j.jmva.2004.04.013
- Bates D, Mächler M, Bolker B, Walker S (2015). Fitting linear mixed-effects models using lme4. J Stat Softw 67(1). https://doi.org/10.18637/jss.v067.i01
- Bates D, Alday PM, Kokandakar AH (2025). Mixed-model log-likelihood evaluation via a blocked Cholesky factorization. arXiv:2505.11674 (preprint). https://doi.org/10.48550/arXiv.2505.11674
- Bermann M, Legarra A, Aguilar I, et al. (2025). Estimation of (co)variance components for very large datasets and complex single-step genomic models. Genet Sel Evol 57:63. https://doi.org/10.1186/s12711-025-01006-9
- Besag J, Kempton R (1986). Statistical analysis of field experiments using neighbouring plots. Biometrics 42:231–251. https://doi.org/10.2307/2531047
- Bijma P, Muir WM, Van Arendonk JAM (2007a). Multilevel selection 1. Genetics 175:277–288. https://doi.org/10.1534/genetics.106.062711
- Bijma P, Muir WM, Ellen ED, Wolf JB, Van Arendonk JAM (2007b). Multilevel selection 2. Genetics 175:289–299. https://doi.org/10.1534/genetics.106.062729
- Bijma P (2010a). Multilevel selection 4: modeling the relationship of indirect genetic effects and group size. Genetics 186:1029–1031. https://doi.org/10.1534/genetics.110.120485
- Bijma P (2010b). Estimating indirect genetic effects: precision of estimates and optimum designs. Genetics 186:1013–1028. https://doi.org/10.1534/genetics.110.120493
- Boldman KG, Van Vleck LD (1991). Derivative-free restricted maximum likelihood estimation in animal models with a sparse matrix solver. J Dairy Sci 74:4337–4343. https://doi.org/10.3168/jds.S0022-0302%2891%2978629-3
- Byrd RH, Lu P, Nocedal J, Zhu C (1995). A limited memory algorithm for bound constrained optimization. SIAM J Sci Comput 16:1190–1208. https://doi.org/10.1137/0916069
- Cappa EP, Cantet RJC (2008). Direct and competition additive effects in tree breeding. Silvae Genet 57:45–56. https://doi.org/10.1515/sg-2008-0008
- Cappa EP, Muñoz F, Sanchez L, Cantet RJC (2015). A novel individual-tree mixed model to account for competition and environmental heterogeneity. Tree Genet Genomes 11:120. https://doi.org/10.1007/s11295-015-0917-3
- Conn AR, Gould NIM, Toint PL (1991). Convergence of quasi-Newton matrices generated by the symmetric rank one update. Math Program 50:177–195. https://doi.org/10.1007/BF01594934
- Costa e Silva J, Potts BM, Bijma P, Kerr RJ, Pilbeam DJ (2013). Genetic control of interactions among individuals. New Phytol 197:631–641. https://doi.org/10.1111/nph.12035
- Costa e Silva J, Kerr RJ (2013). Accounting for competition in genetic analysis, with particular emphasis on forest genetic trials. Tree Genet Genomes 9:1–17. https://doi.org/10.1007/s11295-012-0521-8
- Costa e Silva J, Potts BM, Gilmour AR, Kerr RJ (2017). Genetic-based interactions among tree neighbors. Heredity 119:125–135. https://doi.org/10.1038/hdy.2017.25
- Covarrubias-Pazaran G (2016). Genome-assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11:e0156744. https://doi.org/10.1371/journal.pone.0156744
- Draper NR, Guttman I (1980). Incorporating overlap effects from neighbouring units into response surface models. Appl Stat 29:128–134. https://doi.org/10.2307/2986297
- Erisman AM, Tinney WF (1975). On computing certain elements of the inverse of a sparse matrix. Commun ACM 18:177–179. https://doi.org/10.1145/360680.360704
- Forst E, Enjalbert J, Allard V, et al. (2019). A generalized statistical framework to assess mixing ability from incomplete mixing designs. Field Crops Res 242:107571. https://doi.org/10.1016/j.fcr.2019.107571
- Foulley JL, van Dyk DA (2000). The PX-EM algorithm for fast stable fitting of Henderson's mixed model. Genet Sel Evol 32:143–163. https://doi.org/10.1186/1297-9686-32-2-143
- Gardner JR, Pleiss G, Bindel D, Weinberger KQ, Wilson AG (2018). GPyTorch: blackbox matrix-matrix Gaussian process inference with GPU acceleration. arXiv:1809.11165 (NeurIPS 2018). https://doi.org/10.48550/arXiv.1809.11165
- George A (1973). Nested dissection of a regular finite element mesh. SIAM J Numer Anal 10:345–363. https://doi.org/10.1137/0710032
- Giles MB (2008). Collected matrix derivative results for forward and reverse mode algorithmic differentiation. Lect Notes Comput Sci Eng 64:35–44. https://doi.org/10.1007/978-3-540-68942-3_4
- Gilmour AR, Thompson R, Cullis BR (1995). Average information REML. Biometrics 51:1440–1450. https://doi.org/10.2307/2533274
- Graser HU, Smith SP, Tier B (1987). A derivative-free approach for estimating variance components in animal models by restricted maximum likelihood. J Anim Sci 64:1362–1370. https://doi.org/10.2527/jas1987.6451362x
- Griffing B (1967). Selection in reference to biological groups I. Aust J Biol Sci 20:127–140. https://doi.org/10.1071/BI9670127
- Groeneveld E (1994). A reparameterization to improve numerical optimization in multivariate REML (co)variance component estimation. Genet Sel Evol 26:537–545. https://doi.org/10.1186/1297-9686-26-6-537
- Hadfield JD, Wilson AJ (2007). Multilevel selection 3. Genetics 177:667–668. https://doi.org/10.1534/genetics.107.075622
- Harville DA (1977). Maximum likelihood approaches to variance component estimation and to related problems. J Am Stat Assoc 72:320–338. https://doi.org/10.1080/01621459.1977.10480998
- Higham NJ (1988). Computing a nearest symmetric positive semidefinite matrix. Linear Algebra Appl 103:103–118. https://doi.org/10.1016/0024-3795%2888%2990223-6
- Horn B, Nikoloski Z, Lippert C (2026). TorchLIMIX: GPU-accelerated multivariate genome-wide association studies. bioRxiv (preprint). https://doi.org/10.64898/2026.04.29.721342
- Jennrich RI, Schluchter MD (1986). Unbalanced repeated-measures models with structured covariance matrices. Biometrics 42:805–820. https://doi.org/10.2307/2530695
- Johnson DL, Thompson R (1995). Restricted maximum likelihood estimation of variance components for univariate animal models using sparse matrix techniques and average information. J Dairy Sci 78:449–456. https://doi.org/10.3168/jds.S0022-0302%2895%2976654-1
- Kang HM, Zaitlen NA, Wade CM, et al. (2008). Efficient control of population structure in model organism association mapping. Genetics 178:1709–1723. https://doi.org/10.1534/genetics.107.080101
- Karypis G, Kumar V (1998). A fast and high quality multilevel scheme for partitioning irregular graphs. SIAM J Sci Comput 20:359–392. https://doi.org/10.1137/S1064827595287997
- Kempton RA (1982). Adjustment for competition between varieties in plant breeding trials. J Agric Sci 98:599–611. https://doi.org/10.1017/S0021859600054381
- Knott EA, Mundt CC (1990). Mixing ability analysis of wheat cultivar mixtures under diseased and nondiseased conditions. Theor Appl Genet 80:313–320. https://doi.org/10.1007/BF00210065
- Kristensen K, Nielsen A, Berg CW, Skaug H, Bell BM (2016). TMB: automatic differentiation and Laplace approximation. J Stat Softw 70(5). https://doi.org/10.18637/jss.v070.i05
- Lee SH, van der Werf JHJ (2016). MTG2: an efficient algorithm for multivariate linear mixed model analysis based on genomic information. Bioinformatics 32:1420–1422. https://doi.org/10.1093/bioinformatics/btw012
- Lin L, Yang C, Meza JC, et al. (2011). SelInv: an algorithm for selected inversion of a sparse symmetric matrix. ACM TOMS 37:40. https://doi.org/10.1145/1916461.1916464
- Lindstrom MJ, Bates DM (1988). Newton-Raphson and EM algorithms for linear mixed-effects models for repeated-measures data. J Am Stat Assoc 83:1014–1022. https://doi.org/10.1080/01621459.1988.10478693
- Lippert C, Listgarten J, Liu Y, et al. (2011). FaST linear mixed models for genome-wide association studies. Nat Methods 8:833–835. https://doi.org/10.1038/nmeth.1681
- Liu C, Rubin DB, Wu YN (1998). Parameter expansion to accelerate EM: the PX-EM algorithm. Biometrika 85:755–770. https://doi.org/10.1093/biomet/85.4.755
- Liu JWH (1985). Modification of the minimum-degree algorithm by multiple elimination. ACM TOMS 11:141–153. https://doi.org/10.1145/214392.214398
- Liu Z, Goddard ME, Reinhardt F, Reents R (2014). A single-step genomic model with direct estimation of marker effects. J Dairy Sci 97:5833–5850. https://doi.org/10.3168/jds.2014-7924
- Masuda Y, Baba T, Suzuki M (2014). Application of supernodal sparse factorization and inversion to the estimation of (co)variance components by residual maximum likelihood. J Anim Breed Genet 131:227–236. https://doi.org/10.1111/jbg.12058
- Masuda Y, Aguilar I, Tsuruta S, Misztal I (2015). Technical note: acceleration of sparse operations for average-information REML analyses with supernodal methods and sparse-storage refinements. J Anim Sci 93:4670–4674. https://doi.org/10.2527/jas.2015-9395
- Meyer K (1985). Maximum likelihood estimation of variance components for a multivariate mixed model with equal design matrices. Biometrics 41:153–165. https://doi.org/10.2307/2530651
- Meyer K (1989). Restricted maximum likelihood to estimate variance components for animal models with several random effects using a derivative-free algorithm. Genet Sel Evol 21:317–340. https://doi.org/10.1186/1297-9686-21-3-317
- Meyer K, Smith SP (1996). Restricted maximum likelihood estimation for animal models using derivatives of the likelihood. Genet Sel Evol 28:23–49. https://doi.org/10.1186/1297-9686-28-1-23
- Meyer K (2007). WOMBAT: a tool for mixed model analyses in quantitative genetics by REML. J Zhejiang Univ Sci B 8:815–821. https://doi.org/10.1631/jzus.2007.B0815
- Mishchenko K, Holmgren S, Rönnegård L (2007). Efficient implementation of the AI-REML iteration for variance component QTL analysis. arXiv:0709.0625 (research report). https://doi.org/10.48550/arXiv.0709.0625
- Mishchenko K, Holmgren S, Rönnegård L (2008). Newton-type methods for REML estimation in genetic analysis of quantitative traits. J Comput Methods Sci Eng 8:53–67. https://doi.org/10.3233/JCM-2008-81-203
- Misztal I, Perez-Enciso M (1993). Sparse matrix inversion for restricted maximum likelihood estimation of variance components by expectation-maximization. J Dairy Sci 76:1479–1483. https://doi.org/10.3168/jds.S0022-0302%2893%2977478-0
- Montazeaud G, Flutre T, Ballini E, et al. (2022). From cultivar mixtures to allelic mixtures. New Phytol 233:2573–2584. https://doi.org/10.1111/nph.17915
- Moore AJ, Brodie ED, Wolf JB (1997). Interacting phenotypes and the evolutionary process I. Evolution 51:1352–1362. https://doi.org/10.1111/j.1558-5646.1997.tb01458.x
- Moré JJ, Sorensen DC (1983). Computing a trust region step. SIAM J Sci Stat Comput 4:553–572. https://doi.org/10.1137/0904038
- Muir WM (2005). Incorporation of competitive effects in forest tree or animal breeding programs. Genetics 170:1247–1259. https://doi.org/10.1534/genetics.104.035956
- Murray I (2016). Differentiation of the Cholesky decomposition. arXiv:1602.07527. https://doi.org/10.48550/arXiv.1602.07527
- Nocedal J, Wright SJ (2006). Numerical Optimization, 2nd ed. Springer. https://doi.org/10.1007/978-0-387-40065-5
- Pélissier R, Ballini E, Temple C, et al. (2023). The genetic identity of neighboring plants in intraspecific mixtures modulates disease susceptibility of both wheat and rice. PLoS Biol 21:e3002287. https://doi.org/10.1371/journal.pbio.3002287
- Pinder T, Dodd D (2022). GPJax: a Gaussian process framework in JAX. J Open Source Softw 7:4455. https://doi.org/10.21105/joss.04455
- Pinheiro JC, Bates DM (1996). Unconstrained parametrizations for variance-covariance matrices. Stat Comput 6:289–296. https://doi.org/10.1007/BF00140873
- Rodriguez A, Kim Y, Nandi TN, et al. (2026). SAIGE-GPU: accelerating genome- and phenome-wide association studies using GPUs. Bioinformatics. https://doi.org/10.1093/bioinformatics/btag032
- Sato Y, Yamamoto E, Shimizu KK, Nagano AJ (2021). Neighbor GWAS. Heredity 126:597–614. https://doi.org/10.1038/s41437-020-00401-w
- Sato Y, Hamazaki K (2026). Dissecting oligogenic and polygenic indirect genetic effects through the lens of neighbor genotypic identity. eLife, reviewed preprint. https://doi.org/10.7554/eLife.111650.1
- Searle SR, Casella G, McCulloch CE (1992). Variance Components. Wiley. https://doi.org/10.1002/9780470316856
- Smith SP, Graser HU (1986). Estimating variance components in a class of mixed models by restricted maximum likelihood. J Dairy Sci 69:1156–1165. https://doi.org/10.3168/jds.S0022-0302%2886%2980516-1
- Smith SP (1995). Differentiation of the Cholesky algorithm. J Comput Graph Stat 4:134–147. https://doi.org/10.1080/10618600.1995.10474671
- Strandén I, Garrick DJ (2009). Technical note: derivation of equivalent computing algorithms for genomic predictions and reliabilities of animal merit. J Dairy Sci 92:2971–2975. https://doi.org/10.3168/jds.2008-1929
- Stringer JK, Cullis BR, Thompson R (2011). Joint modeling of spatial variability and within-row interplot competition to increase the efficiency of plant improvement. JABES 16:269–281. https://doi.org/10.1007/s13253-010-0051-5
- Taskinen M, Mäntysaari EA, Strandén I (2017). Single-step SNP-BLUP with on-the-fly imputed genotypes and residual polygenic effects. Genet Sel Evol 49:36. https://doi.org/10.1186/s12711-017-0310-9
- Thompson R (2008). Estimation of quantitative genetic parameters. Proc R Soc B 275:679–686. https://doi.org/10.1098/rspb.2007.1417
- Wolf JB, Brodie ED, Cheverud JM, Moore AJ, Wade MJ (1998). Evolutionary consequences of indirect genetic effects. Trends Ecol Evol 13:64–69. https://doi.org/10.1016/S0169-5347%2897%2901233-0
- Wolfinger R, Tobias R, Sall J (1994). Computing Gaussian likelihoods and their derivatives for general linear mixed models. SIAM J Sci Comput 15:1294–1310. https://doi.org/10.1137/0915079
- Wuest SE, Niklaus PA (2018). A plant biodiversity effect resolved to a single chromosomal region. Nat Ecol Evol 2:1933–1939. https://doi.org/10.1038/s41559-018-0708-y
- Wuest SE, Peter R, Niklaus PA (2021). Ecological and evolutionary approaches to improving crop variety mixtures. Nat Ecol Evol 5:1068–1077. https://doi.org/10.1038/s41559-021-01497-x
- Yang J, Lee SH, Goddard ME, Visscher PM (2011). GCTA: a tool for genome-wide complex trait analysis. Am J Hum Genet 88:76–82. https://doi.org/10.1016/j.ajhg.2010.11.011
- Zhou H, Stephens M (2012). Genome-wide efficient mixed-model analysis for association studies. Nat Genet 44:821–824. https://doi.org/10.1038/ng.2310
- Zhou X, Stephens M (2014). Efficient multivariate linear mixed model algorithms for genome-wide association studies. Nat Methods 11:407–409. https://doi.org/10.1038/nmeth.2848

Manuals, software and sources without a DOI

- Gilmour AR, Gogel BJ, Cullis BR, Welham SJ, Thompson R (2022). ASReml User Guide, Functional Specification, Release 4.2. VSN International. https://asreml.kb.vsni.co.uk/wp-content/uploads/sites/3/ASReml-4.2-Functional-Specification.pdf (sections 5.8, 6.10, 7.7.4, 11.5 consulted)
- Meyer K (2017). WOMBAT, a program for mixed model analyses by restricted maximum likelihood. User notes, version of 31 May 2017. https://forge-dga.jouy.inra.fr/attachments/download/8063/WombatManual-V20170531.pdf (sections 1, 4.10.6.8, 5.2–5.3, appendices A.1 and A.5 consulted)
- Boldman K, Kriese L, Van Tassell C, Kachman S, Dodenhoff J (2002). MTDFREML Workshop Notes. USDA-ARS. https://ars.usda.gov/ARSUserFiles/80420530/MTDFREML/MTDFREML_Workshop_Notes.pdf
- Marchal A, Raimondi D (2026). pyreml, general-purpose linear mixed model solver, version 0.5.3 (source code examined locally). https://doi.org/10.5281/zenodo.21394319
- SuiteSparse AMD, header amd.h (dense row threshold). https://github.com/DrTimothyAldenDavis/SuiteSparse/blob/stable/AMD/Include/amd.h
- Takahashi K, Fagan J, Chin MS (1973). Formation of a sparse bus impedance matrix and its application to short circuit study. 8th PICA Conference, Minneapolis, pp. 63–69. No DOI, not consulted.
- Meyer K (2005). Ordering strategies to reduce computational requirements in variance component estimation. Proc Assoc Advmt Anim Breed Genet 16:282–285. No DOI, not consulted.
