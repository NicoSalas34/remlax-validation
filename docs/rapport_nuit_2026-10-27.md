remlax · nuit du 26 au 27 octobre 2026

# remlax face à kmmr et asreml-R

Optimisation du solveur, validation complète et banc sur sept modèles de génétique quantitative, CPU et GPU, version figée v10.

## En bref

remlax atteint sur chacun des sept cas un optimum au moins aussi haut que kmmr et asreml-R. Sur quatre cas, il est strictement plus haut. Les deux modèles DROPS en une étape font décrocher kmmr, qui s’arrête après deux itérations. asreml plafonne leurs corrélations ar1 à 0,98, colle des variances spécifiques à zéro sur drops2_fa3 et s’arrête 3,8 unités sous remlax sur krause_fa2.

Sur A100, remlax est le plus rapide sur tous les modèles Krause : 22 fois kmmr sur krause_fa2, près de 4 fois sur krause_us. Sur CPU, il passe devant kmmr sur krause_us (316 s contre 359 s) et drops2_fa3, il est au niveau de kmmr sur la GBLUP et drops2_fa1, et il reste plus lent sur krause_fa2 (995 s contre 813 s), où il monte 20 unités plus haut. Sur les petits modèles, l’essentiel de son temps est fixe : environ 2,5 s pour démarrer le solveur Python, une fois par session R.

×22 (krause_fa2 sur A100 : 44,4 s contre 961 s pour kmmr, optimum plus haut).

×3,7 (krause_us sur A100 : 22,1 s contre 82,6 s).

357 → 29 s (drops1_fa1 sur CPU, de la v1 à la v10, mêmes nœuds).

36 / 37 (cas de la vérification par simulation ; le dernier est un aléa connu).

## Résultats du banc

kmmr et remlax tournent sur les mêmes nœuds de genotoul, AMD EPYC 7513 avec 8 cœurs ou A100 80 Go, dans la même image. asreml-R tourne sur meso, où sa licence est disponible : ses temps ne valent qu’en ordre de grandeur. Le temps de kmmr comprend l’inversion de K faite dans R avant l’appel. Celui de remlax comprend tout `rx_fit()`, démarrage du solveur compris, au second ajustement de la session (programmes compilés relus sur disque). ▼ marque un optimum plus bas que le meilleur des trois.

| cas | log-vraisemblance REML à l’optimum |  |  | CPU, 8 cœurs (s) |  |  | A100 (s) |  |
|----|----|----|----|----|----|----|----|----|
|  | kmmr | asreml-R | remlax | kmmr | remlax | asreml | kmmr | remlax |
| DROPS deux étapes, fa1 | -492.86 | -492.86 | -492.86 | 6.2 | 7.2 | 0.5 | 8.3 | 8.2 |
| DROPS deux étapes, fa3 | ▼ -368.06 | ▼ -378.25 | -367.68 | 16.0 | 9.7 | 2.4 | 19.6 | 10.0 |
| DROPS une étape, fa1 | ▼ -2154.65 | ▼ -1161.86 | -1161.21 | 3.9 | 29.1 | 15.3 | 5.0 | 42.7 |
| DROPS une étape, us(5) | ▼ -2169.30 | ▼ -1146.72 | -1146.19 | 5.8 | 18.0 | 10.5 | 5.6 | 28.3 |
| Krause, GBLUP | 389.57 | 389.57 | 389.57 | 9.1 | 9.0 | 60.2 | 7.4 | 9.2 |
| Krause, fa2 × 5 traitements | ▼ 7109.92 | ▼ 7126.36 | 7130.18 | 813 | 995 | 7429 | 961 | 44.4 |
| Krause, us(5) × 5 traitements | 7138.93 | 7138.93 | 7138.93 | 359 | 316 | 7530 | 82.6 | 22.1 |

kmmr à la tolérance 1e-6 sur la log-vraisemblance. asreml-R avec ses réglages par défaut, workspace de 30 Go pour Krause. Sur drops1, asreml retrouve l’optimum de remlax à 0,02 près quand on fige comme lui les corrélations à 0,98.

### Itérations et coût d’une évaluation

| cas | it. kmmr | it. asreml | pas AI essayés remlax | remlax s/éval CPU | remlax s/éval A100 | version |
|----|----|----|----|----|----|----|
| DROPS deux étapes, fa1 | 23 | 13 | 24 | 0.071 | 0.030 | v10 |
| DROPS deux étapes, fa3 | 38 | 30 | 138 | 0.070 | 0.033 | v10 |
| DROPS une étape, fa1 | 2 | 24 | 76 | 0.435 | 0.607 | v10 |
| DROPS une étape, us(5) | 2 | 23 | 30 | 0.416 | 0.608 | v10 |
| Krause, GBLUP | 11 | 7 | 4 | 0.676 | 0.056 | v10 |
| Krause, fa2 × 5 traitements | 14 | 15 | 56 | 28.122 | 1.252 | v10 |
| Krause, us(5) × 5 traitements | 12 | 9 | 9 | 30.670 | 1.242 | v10 |

Les pas AI essayés de remlax comptent aussi les pas rejetés, évalués sur la seule vraisemblance. Sur les modèles fa, la convergence de l’AI reste linéaire en fin de course : c’est la limite principale de remlax sur CPU.

### Progression de remlax au fil des versions figées

| cas, matériel | v1 | v2 | v3 | v6 | v10 |
|----|----|----|----|----|----|
| DROPS deux étapes, fa3, CPU | 186 s -367.68 | 12.9 s -367.68 | 10.4 s -367.68 | 10.4 s -367.68 | 9.7 s -367.68 |
| DROPS une étape, fa1, CPU | 357 s -1161.21 | 25.4 s -1161.21 | 415 s -1161.21 | 34.5 s -1161.21 | 29.1 s -1161.21 |
| Krause, GBLUP, CPU | 11.7 s 389.57 | 14.4 s 389.57 | 11.1 s 389.57 | 10.0 s 389.57 | 9.0 s 389.57 |
| Krause, us(5) × 5 traitements, CPU |  | 1286 s 7138.93 |  | 444 s 7138.93 | 316 s 7138.93 |
| Krause, fa2 × 5 traitements, A100 | 480 s 7130.18 |  |  | 66.6 s 7130.18 | 44.4 s 7130.18 |
| Krause, us(5) × 5 traitements, A100 | 266 s 7138.93 |  |  | 28.0 s 7138.93 | 22.1 s 7138.93 |

v1 à v3 : premières optimisations (départs, essais sur la valeur seule, fils BLAS, blocs pleins). v6 : brides de l’AI et voie GPU sans transfert. v9 : cache de la dernière évaluation et règle de Heywood. v10 : un pas AI accepté reprend la factorisation de son essai. La v3 trouvait drops1_fa1 en 415 s sur CPU et le manquait sur GPU ; la v10 le trouve partout.

### Où passe le temps de remlax (v10, à chaud)

#### CPU, 8 cœurs

| cas | démarrage | préparation | 1re évaluation | optimisation | finition | sorties | total solveur |
|----|----|----|----|----|----|----|----|
| DROPS deux étapes, fa1 | 2.6 | 0.3 | 2.2 | 1.9 | 0.0 | 0.1 | 7.2 |
| DROPS deux étapes, fa3 | 2.6 | 0.3 | 0.3 | 6.3 | 0.0 | 0.1 | 9.6 |
| DROPS une étape, fa1 | 2.5 | 2.3 | 0.6 | 23.0 | 0.0 | 0.2 | 28.9 |
| DROPS une étape, us(5) | 2.4 | 2.3 | 0.6 | 12.1 | 0.0 | 0.2 | 17.9 |
| Krause, GBLUP | 3.0 | 1.7 | 0.8 | 3.0 | 0.0 | 0.1 | 8.8 |
| Krause, fa2 × 5 traitements | 2.6 | 2.4 | 28.6 | 960 | 0.0 | 0.6 | 995 |
| Krause, us(5) × 5 traitements | 2.5 | 2.8 | 29.3 | 280 | 0.0 | 0.6 | 316 |

#### A100

| cas | démarrage | préparation | 1re évaluation | optimisation | finition | sorties | total solveur |
|----|----|----|----|----|----|----|----|
| DROPS deux étapes, fa1 | 2.8 | 0.3 | 0.5 | 0.8 | 0.0 | 0.1 | 8.1 |
| DROPS deux étapes, fa3 | 2.8 | 0.3 | 0.5 | 2.7 | 0.0 | 0.1 | 10.0 |
| DROPS une étape, fa1 | 2.7 | 2.9 | 0.8 | 32.2 | 0.0 | 0.3 | 42.6 |
| DROPS une étape, us(5) | 2.8 | 3.0 | 0.8 | 17.7 | 0.0 | 0.3 | 28.1 |
| Krause, GBLUP | 2.8 | 1.5 | 0.5 | 0.4 | 0.0 | 0.1 | 9.0 |
| Krause, fa2 × 5 traitements | 2.7 | 2.8 | 1.8 | 33.1 | 0.0 | 0.1 | 44.3 |
| Krause, us(5) × 5 traitements | 2.8 | 2.8 | 1.8 | 10.7 | 0.0 | 0.1 | 22.0 |

## Ce qui a changé dans remlax

| changement | défaut corrigé | effet mesuré |
|----|----|----|
| Chrono par étapes (`rx_fit()$chrono`) | aucune mesure de la préparation, de la compilation ni des sorties | chaque étape, du lancement par R à la relecture du résultat |
| Départ sur les résidus des effets fixes, par composante et par section | la variance de y brut contient les écarts de moyenne entre environnements | krause_us (600 génotypes) : 17 → 7 itérations |
| Exposition divisée par les couples (ligne, composante) chargés | division par n × t : départs t fois trop grands en format long | drops2 : 13 itérations amorties en moins |
| Essais de pas AI sur la vraisemblance seule | chaque essai rejeté coûtait gradient et matrice AI | drops2 : 42 → 12 s |
| OpenBLAS de NumPy à un fil pendant le MME | les OpenBLAS de NumPy et de SciPy se disputaient les cœurs | évaluation 1,46 → 0,25 s (m = 1428, 8 cœurs) |
| Greffon CUDA non chargé sur CPU | 2,0 s d’initialisation du backend à chaque lancement | démarrage 4,3 → 2,2 s |
| Cache de compilation persistant | recompilation à chaque lancement | GPU : premier appel 3,4 → 0,6 s |
| Table des composantes vectorisée, sonde Python à la demande | une opération JAX par entrée de Σ ; un Python lancé pour rien | 1,9 → 0,18 s ; 1,5 s |
| Blocs pleins dans le MME dense, différés hors du motif creux | 178 millions d’entrées regroupées par bincount à chaque évaluation | préparation Krause 18 → 3 s |
| Sur GPU, C et C⁻¹ restent sur la carte, blocs groupés par terme | 2,8 Go aller-retour par évaluation ; 406 opérations à compiler | évaluation Krause 6,9 → 1,25 s ; drops2 à froid 39 → 12 s |
| Dernière évaluation du MME en cache | finition et sorties refactorisaient C au point final | krause_us sur CPU : 444 → 384 s (v6 → v9) |
| Un pas AI accepté reprend la factorisation de son essai | l’évaluation complète refactorisait C juste après l’essai sur la valeur | environ 10 s de Cholesky en moins par pas accepté sur Krause en CPU |
| AI bridée : variance ÷10, paramètre de niveau 0,5 au plus par pas | sauts entre bassins d’une vraisemblance multimodale | drops1 : bon optimum quel que soit le nombre de fils |
| Variance bridée trois pas de suite libre de rejoindre le plancher | cas de Heywood d’un fa : dix paliers de plus | krause_fa2 : 101 → 56 pas essayés |
| Sortie du plancher après L-BFGS-B et après l’AI | une variance au plancher a un gradient nul en log | drops2 en forme V : −2734,67 → −1753,64 |
| METIS et CHOLMOD, ce dernier forcé dans notre ordre | SuperLU manquait de mémoire ; CHOLMOD pouvait rendre un facteur dans son propre ordre | mvD7 en 26 Go sur CPU ; 7 erreurs de la vérification → 0 |
| `engine = "auto"`, `optimizer = "auto"`, `theta_init` avec NA | choix du moteur laissé à l’utilisateur ; appel préalable pour fixer un paramètre | MME et AI quand le modèle s’y prête |

Essayés puis écartés : une correction quasi-Newton (SR1) de la matrice AI, qui augmentait les itérations partout, et un test de rapport à la manière des régions de confiance, sans effet sur les optimums et jusqu’à 60 % plus lent (laissé en option, `REMLAX_AI_RAPPORT`). Une bride ÷1000 sur les variances accélérait les fa mais manquait l’optimum de drops1_us.

## Validation de la version finale

| contrôle | résultat |
|----|----|
| tests Python, CPU | 257 passés, 1 sauté |
| tests Python, GPU (L40S) | 253 passés, 1 sauté |
| testthat | tout passe ; 1 test sauté faute du dépôt IGE |
| `R CMD check` | Status: OK, sans avertissement ni note |
| vérification par simulation (37 cas, 4 voies, arbitre REML dense en R) | 36 cas passés ; mt_g_fa2 échoue sur sa graine 1 sur le seul critère de récupération (une composante à 3,06 erreurs-types, même vraisemblance exacte sur les quatre voies) et passe sur 4 graines sur 5 |
| parité CPU/GPU (L40S, 7 paquets du banc) | évaluations à 1e-13 près ou mieux, ajustements complets à 2e-10 près |

Quatre tests ont changé, chacun justifié dans son commentaire. Sur un jeu à vraisemblance plate et multimodale, le moteur par blocs et la finition à région de confiance trouvent désormais un meilleur point que le moteur dense : les tests vérifient qu’ils ne font pas pire. La finition à région de confiance ne fait plus gagner d’évaluations depuis la correction des départs : son test ne compare plus que l’optimum. Un seuil analytique de 1e-8 était trop strict pour les sommes du GPU : il passe à 1e-7.

## Limites et décisions à prendre

- Le défaut reste `engine = "dense"` et `optimizer = "lbfgsb"`. Passer à `"auto"` donnerait la voie la plus rapide sans réglage, mais changerait le chemin d’optimisation de tes scripts, dont ceux du pipeline IGE. Je ne l’ai pas fait sans ton accord.
- Les nouveaux départs changent le chemin de tous les ajustements. Les optimums sont les mêmes sur la vérification et sur le banc, mais les ajustements du pipeline IGE devront être relancés avant toute comparaison avec les anciens.
- Les vraisemblances multimodales restent un risque. Les brides de l’AI ont réglé tous les cas observés, sans garantir l’optimum global : `n_restarts` reste la parade en cas de doute.
- Sur CPU, krause_fa2 reste plus lent que kmmr : l’AI y fait 19 itérations contre 14, plus des essais de pas réduits en début de course. C’est la prochaine cible.
- Le démarrage du solveur Python coûte environ 2,5 s une fois par session R ; c’est l’essentiel du temps sur les petits modèles.
- asreml-R : la licence de meso n’accepte que quelques exécutions simultanées ; krause_us et krause_fa2 y ont demandé deux heures chacun. L’ajustement asreml de mvD7 plante après l’itération (« double free or corruption »).
- Le dossier `_archive_bancs_20261021` est exclu de la construction du paquet, mais pas de git : à garder ou non.

## Prêt à pousser

Rien n’est commité. Dans le dépôt remlax, 45 fichiers suivis sont modifiés, et ceux-ci sont nouveaux et à ajouter :

    R/remlax_bayes.R  R/remlax_lag.R  R/remlax_mme.R  R/remlax_serveur.R
    inst/python/remlax/{ai,bayes,blocs,chrono,decoupage,econome,finition,lag,mme,mme_cli,mme_hybride,selinv,serveur}.py
    inst/python/remlax/selinv_c.c
    man/rx_bayes.Rd  man/rx_border.Rd  man/rx_fit_mme.Rd  man/rx_lag.Rd  man/rx_prior.Rd  man/rx_server_stop.Rd
    manual/ (sources ; les sorties LaTeX sont ignorées)
    tests-python/test_{bayes,blocs,decalage,econome,finition,mme_hybride,plancher,poids}.py
    tests/testthat/test-decalage.R  vignettes/structural.Rmd

Dans remlax-validation, à ajouter : `banc_kmmr/` (cas, pilotes kmmr, remlax et asreml, synthèse), `verif_simulations/`, `docs/paysage_logiciels_reml.md` et ce rapport (`docs/rapport_nuit_2026-10-27.md`). L’image `container/ige_reml.def` du dépôt IGE gagne pymetis et cvxopt ; l’image reconstruite est sur genotoul (`/work/user/nsalas/kmmr_banc/ige_reml_v2.sif`).

Versions figées, sources et empreintes : `/work/user/nsalas/kmmr_banc/` sur genotoul (`rlib_remlax` pour la v1, puis `rlib_remlax_v2` à `v10`, `VERSION_FIGEE_v*.txt`, `remlax_fige_v*_src.tgz`).
