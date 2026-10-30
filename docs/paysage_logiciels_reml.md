# Paysage des solveurs REML et positionnement de remlax

Recherche faite le 26 octobre 2026 sur sources primaires : documentation officielle, manuels, code source, articles avec DOI vérifié. Les annexes A, B et C détaillent chaque logiciel et chaque mécanisme, avec leurs sources numérotées. Cette synthèse ne reprend que les constats utiles à la décision. Les mesures de remlax citées ici sont internes et non publiées (mvD7 : n = 12 760, 51 termes, 198 paramètres).

Question posée : remlax peut-il remplacer ASReml sur les calculs lourds, ou devenir une alternative libre du type de sommer, en plus efficace ?

## 1. Ce qui a changé en 2026

Le créneau « REML libre, structures façon asreml, GPU » était vide il y a un an. Il ne l'est plus.

- pyreml (Marchal, Trottier, Bry et Raimondi, IGMM et IMAG, Montpellier) est paru sur PyPI en juin 2026. Sa prépublication bioRxiv date du 24 octobre 2026 (doi:10.64898/2026.10.20.757556). C'est un REML exact sur GPU en PyTorch, à V dense ou Woodbury dense, avec gradient analytique. Il n'a ni matrices creuses, ni AI, ni erreurs-types des variances. Son banc l'oppose à WOMBAT, pas à ASReml (annexe B).
- kmmr (Melsen, van Eeuwijk et Peeters, Wageningen) est paru sur arXiv le 12 octobre 2026 (arXiv:2610.16671). C'est un AI-REML multivarié sur équations du modèle mixte : partie creuse par CHOLMOD sur CPU, inversion de Takahashi, complément de Schur dense sur GPU. Sur un GBLUP de blé à 10 375 génotypes, il annonce 49 s sur A100 contre 989 s pour ASReml-R et 512 s pour sommer. Il cite la forme V comme travail futur (annexe B).
- sommer 4.5.12 (CRAN, 25 octobre 2026) résout désormais par défaut les équations du modèle mixte creuses : CHOLMOD, inversion de Takahashi, ordonnancement METIS optionnel, noyaux CUDA optionnels pour la partie dense. La formulation directe sur V reste disponible mais part vers la dépréciation. Le paquet change très vite : la 4.4.90 a cassé l'interface de contrôle (annexe A).

Les trois projets sont libres (GPL-3). Deux visent explicitement l'amélioration des plantes.

## 2. Comparaison directe

| | remlax | ASReml-R 4.2 | sommer 4.5.12 | kmmr | pyreml 0.5.3 |
|---|---|---|---|---|---|
| Licence | GPL ≥ 3 | commerciale, dès 375 £ (académique) | GPL-3 | GPL-3, expérimental | GPL-3 |
| Formes de calcul | V dense, MME creuses, V par blocs (Woodbury) | MME creuses, partition dense/creux | MME creuses ; V dense en voie de dépréciation | MME creuses + Schur dense | V dense ou Woodbury dense |
| Optimiseur | L-BFGS-B, AI-REML, finition à région de confiance | AI-REML | AI-REML, repli EM/PX-EM | AI-REML | L-BFGS (PyTorch), Adam |
| GPU | oui (JAX) : V, blocs, MME dense | non | CUDA optionnel, noyaux denses seulement | complément de Schur | oui (PyTorch) |
| Ordonnancement du creux | ordre naturel (METIS mesuré, non branché) | non documenté | AMD, METIS optionnel | CHOLMOD | sans objet |
| Structures façon asreml | us, fa, rr, ar1⊗ar1, séparables, K ou K⁻¹ | complètes | us, fa, rr, Kronecker arbitraires, K ou K⁻¹ | us, FA étendu, AR1×AR1, noyaux | iid, diag, full, fa, ar, distances |
| Erreurs-types des variances | oui (AI, Hessien) | oui | oui | oui (AI) | non |
| Détection automatique de structure | oui (termes globaux et blocs) | non : partition fixée par l'utilisateur | non | non | bascule V ou Woodbury selon q < n |
| Voisinage IGE à noyau | oui, portée estimée | à la main | à la main | non | non |
| Vitesse publiée | aucune | référence | 2 fois plus lent que kmmr CPU sur 10 375 génotypes | 20 fois ASReml (GBLUP, A100) | plus de 12 fois WOMBAT (A100) |

Les cases de remlax reposent sur nos mesures internes. Les autres colonnes viennent des annexes A et B.

## 3. Ce qui reste propre à remlax

D'après les sources consultées, aucun autre outil ne réunit les éléments suivants.

- Les formes V et MME dans un même solveur, avec la même vraisemblance et l'AI-REML. pyreml n'a pas de forme MME. kmmr n'a pas de forme V.
- La détection automatique des termes globaux et des blocs, suivie de V = A + UU' traitée par Woodbury avec les blocs en lots (vmap). L'annexe C ne trouve aucune antériorité. Sur mvD7, une évaluation prend 0,90 s sur une tranche MIG de 10 Go, contre 4,85 s pour notre mode économe.
- Un modèle multivarié lourd ajusté sur une tranche MIG de 10 Go. Aucun outil ne rapporte cet usage.
- L'ensemble actif sur le cône semi-défini positif en échelle des variances, et le finisseur de Moré et Sorensen avec mise à jour SR1. L'annexe C n'en trouve pas d'usage dans les logiciels REML.
- L'estimation par REML de la portée continue d'un noyau IGE avec les variances. WOMBAT profile seulement la dilution ; Sato et al. (2021) balaient l'échelle.

Plusieurs briques ont en revanche une antériorité nette (annexe C) : le gradient P − Pyy'P et la forme V dense (Harville 1977), la dérivation de la Cholesky (Smith 1995, Meyer et Smith 1996), la Cholesky supernodale avec inversion sélective (YAMS, WOMBAT, ASReml), la dissection emboîtée sur C (WOMBAT), la vraisemblance tirée d'une matrice bordée (Smith et Graser 1986, Graser et al. 1987), l'AI-REML (Gilmour et al. 1995) et son amortissement (WOMBAT, SAS). Un article devra les citer comme fondations, pas comme apports.

## 4. Ce qui manque à remlax

- Le calcul sur plusieurs GPU ou plusieurs nœuds : DISSECT, REACTA, SAIGE-GPU, ExaGeoStat et DALIA le font. En JAX, jaxmg fournit une Cholesky distribuée avec log-déterminant (annexe B).
- La voie sans matrice, par gradient conjugué préconditionné et traces stochastiques : BOLT-REML, SAIGE, SLEMM, le REML Monte Carlo de BLUPF90 et GPyTorch. Elle seule dépasse quelques centaines de milliers d'individus.
- La précision mixte FP32 puis FP64 de pyreml, utile sur les GPU sans calcul FP64 rapide (L40S).
- Les génotypes compressés sur 2 bits lus en flux vers le GPU (miraculix, Quickdraws).
- L'ordonnancement METIS dans notre moteur MME. sommer et WOMBAT l'ont déjà. Notre mesure sur mvD7 montre un facteur de C divisé par 4 (0,69 Go contre 2,72 Go en ordre naturel).
- Des bancs publiés contre les concurrents actuels.

## 5. Réponse à la question

Remplacer ASReml sur les calculs lourds est un objectif réaliste. ASReml n'a pas de GPU, plafonne à 96 Go et 16 threads, et sa partition dense/creux est manuelle. Sur mvD7, une itération AI d'asreml prend 50 min en partition par défaut et 31 min avec les termes génétiques en dense ; le moteur par blocs de remlax évalue la même vraisemblance en 0,9 s. Mais kmmr revendique déjà un facteur 20 face à ASReml, et la comparaison devra se faire contre lui.

Être « un sommer plus efficace » est devenu plus difficile. sommer 4.5 a adopté les mêmes techniques creuses que remlax et un GPU partiel. La différence ne se fera plus sur la liste des fonctions mais sur des mesures, sur la stabilité de l'interface, et sur ce que sommer n'a pas : la forme V et les blocs sur GPU, la détection automatique, les modèles de voisinage.

## 6. Recommandations

1. Bancs directs, aux mêmes modèles et au même point fixe de vraisemblance, contre sommer 4.5.12 (mmes avec METIS, avec et sans CUDA), kmmr, pyreml et ASReml-R (partition par défaut et dense). Commencer par les familles du banc moteurs et par mvD7.
2. Brancher METIS (ou l'ordre « blocs puis séparateur » de decoupage.py) dans le moteur MME. Le gain est mesuré, et les concurrents l'ont.
3. Stabiliser l'interface R avant toute diffusion : c'est un point faible documenté de sommer.
4. Ne viser le multi-GPU ou la voie sans matrice que si la cible dépasse 10⁵ individus génotypés. Pour les dispositifs d'amélioration des plantes visés, la forme par blocs semble suffire.
5. Pour un article : centrer l'apport sur la réunion des formes V et MME, la détection automatique de structure et les modèles de voisinage ; citer pyreml et kmmr comme travaux contemporains ; ne pas revendiquer de priorité sur le REML GPU.
6. Prendre contact avec l'équipe de pyreml, à Montpellier comme l'utilisateur. Leurs choix (Woodbury dense, précision mixte, grilles d'index) recoupent les nôtres.

## 7. Points non vérifiés

- Sites de WOMBAT, DMU et MTG2 inaccessibles : leurs fiches reposent sur des manuels anciens ou des extraits (annexe A).
- Échelle exacte des bancs de kmmr (nombre d'observations non donné) et écart de log-vraisemblance face à ASReml-R sur un modèle FA1 (−3 948,66 contre −3 937,65), à relire dans le texte complet.
- Actes WCGALP 2022 bloqués : d'éventuels travaux GPU en sélection animale n'ont pas été lus.
- ASReml est fermé : des heuristiques non documentées peuvent exister.
- La priorité « premier REML général sur GPU » reste disputée entre pyreml et kmmr ; REACTA (2014) l'avait fait pour un cas spécialisé.

Les annexes suivent : A, logiciels de référence et alternatives ; B, calcul accéléré ; C, antériorité des mécanismes.

## Annexe A. Logiciels de référence et alternatives pour le REML en génétique quantitative

Recherche faite le 26 octobre 2026. Les sources sont des pages officielles, des manuels, des dépôts de code et des articles avec DOI vérifié sur CrossRef ou PubMed. Les numéros entre crochets renvoient à la liste des sources en fin de document. La mention « non vérifié » signale un point que je n'ai pas pu confirmer sur une source primaire.

Trois serveurs officiels étaient inaccessibles pendant la recherche : le site de WOMBAT (didgeridoo.une.edu.au), le site de DMU (dmu.ghpc.au.dk) et la page MTG2 (Google Sites, connexion exigée). Les informations sur ces trois logiciels reposent donc sur des manuels plus anciens ou sur des extraits de moteur de recherche, signalés comme tels.

### 1. Tableau comparatif

Abréviations. MME : équations du modèle mixte de Henderson. V dense : formulation sur la matrice de variance des observations, de taille n × n. AI : algorithme d'information moyenne. K : matrice de parenté fournie telle quelle. K⁻¹ : son inverse fournie. Fermé : binaire sans code source.

| Logiciel | Licence | Langage | Algorithme | Creux ou dense | us / fa multivarié | ar1⊗ar1 | Parenté | GPU | Dernière version |
|---|---|---|---|---|---|---|---|---|---|
| ASReml-R 4.2 | Commerciale, fermée ; licence académique dès 375 £ [1] | Moteur Fortran (non vérifié), interface R | AI-REML [2][3] | MME creuses, partition dense/creux [2] | Oui : us, fa, facv, xfa, rr [4] | Oui [4] | K ou K⁻¹ (ainverse, matrices creuses) [2] | Non | 4.2.0.480-2, 6 mai 2026 [5] |
| ASReml-SA 4.2 | Commerciale, même licence que ASReml-R [4] | Non vérifié | AI-REML [4] | MME creuses ; espace de travail jusqu'à 96 Go ; 16 threads par défaut [4] | Oui [4] | Oui [4] | K ou K⁻¹ [4] | Non | Release 4.2 (spécification d'octobre 2022) ; build courant non vérifié [4] |
| Echidna | Gratuit hors usage commercial, fermé, version « alpha » [6] | Non vérifié | AI-REML [6] | Non vérifié | Non vérifié | Non vérifié | H⁻¹ et GRM (note de version) [7] | Non vérifié | Build 2.15, 22 octobre 2026 [7] |
| WOMBAT | Gratuit pour la communauté scientifique, exécutable seul [8] | Non vérifié | PX-EM puis AI ; EM ; Powell et Simplex sans dérivées [8] | MME creuses supernodales ; ordonnancement MMD, AMD, METIS ; option --dense [8] | us ; rang réduit par Cholesky partiel ; régression aléatoire [8] | Non documenté dans le manuel 2017 | K⁻¹ via fichier .gin [8] | Non | Manuel du 31 mai 2017 ; version courante non vérifiée (site inaccessible) |
| BLUPF90+ (REMLF90, AIREMLF90), GIBBSF90+ | Gratuit pour la recherche ; version libre limitée à 30 000 génotypés ; commercial sur accord ; binaires, anciennes sources « éducatives » [9][10][11] | Fortran 90/95 [12] | EM-REML (remlf90), AI-REML (airemlf90), Gibbs [11] | MME creuses ; FSPAK puis YAMS supernodal [13] | Multicaractère oui [11] ; fa non vérifié | Non natif ; possible via breedR | A⁻¹ par pedigree, H⁻¹ single-step [11][14] | Non | Page programmes mise à jour le 5 juin 2026 [11] |
| DMU | Gratuit pour la recherche, exécutables ; commercial sur accord (extrait de recherche) [15] | Non vérifié | AI-REML de Jensen et al. 1997, repli EM (DMUAI) (extrait) [15] | MME creuses, FSPAK [16] | Multivarié (extrait) [15] | Non vérifié | Non vérifié | Non | Guide v6 release 5.6, juin 2024 (extrait) [15] |
| MiX99 | Commerciale ; licences Standard, Core, University ; élevage seulement [17] | Non vérifié | MC-EM-REML ; AI-REML cité en 2024 [17] | Résolution itérative ; MPI et multithread [17] | Multicaractère [17] | Non | Non vérifié | Non | Release IX/2025 (extrait) [17] |
| MiXBLUP 3 | Commerciale : 2 400 à 24 790 €/an ; essai gratuit d'un mois [18] | Non vérifié | Non précisé ; cœur dérivé de MiX99 [19] | Non vérifié | Non vérifié | Non | Non vérifié | Non | MiXBLUP 3 (tarifs 2026) [18] |
| MTG2 | Gratuit ; sources v2.15 déposées sur Dryad ; licence logicielle non précisée [20] | Fortran [20] | AI-REML direct + décomposition propre de la GRM [21] | V dense [21] | Multicaractère, régression aléatoire, plusieurs GRM [21] | Non | K dense [21] | Non vérifié | v2.16, 17 août 2020 [20] |
| GCTA (GREML) | GPL-3, parties en LGPL [22] | C++ [22] | AI-REML direct, HE [23][24] | V dense ; fastGWA-REML sur GRM creuse [25] | Bivarié seulement [16] | Non | K (GRM) | Non | v1.94.1 (étiquette 2022) ; dernier commit juillet 2025 [22] |
| LDAK | GPL-3 [26] | C (non vérifié) | REML sur noyaux multiples [27] ; détail de l'algorithme non vérifié | V dense | Non vérifié | Non | K | Non | ldak7.01 ; dépôt actif en octobre 2026 [26] |
| sommer | GPL-3 [28] | R + C++ (Armadillo, Eigen, CHOLMOD) [28] | AI-REML sur MME creuses (mmes, défaut) ; inversion directe NR/AI (mmer, formulation directe) ; récupération EM/PX-EM [28] | Les deux ; CHOLMOD, inverse sélective de Takahashi, METIS optionnel, PCG stochastique [28] | us, diag, fa, rr [28] | Oui, vsm(ar1m(), ar1m()) [28] | K ou K⁻¹ (attribut inverse) [28] | CUDA optionnel, noyaux denses seulement [28] | 4.5.12, 25 octobre 2026 [28] |
| lme4 | GPL ≥ 2 [29] | R + C++ (Eigen) [29] | PLS + optimiseur sans dérivées [30] | Cholesky creux [30] | Univarié ; us, diag, cs, ar1 par terme depuis 2.0 [31] | Non | Non natif | Non | 2.1-0, 23 octobre 2026 [29] |
| lme4qtl | GPL-3 [32] | R | lme4 avec Z·L, L facteur de K [33] | Creux si K creuse [33] | Via lme4 | Non | K [33] | Non | Dernier commit janvier 2021 [32] |
| pedigreemm | GPL ≥ 2 [34] | R | lme4 avec facteur de A [35] | Creux | Non | Non | Pedigree seulement [36] | Non | 0.3-5, 14 octobre 2024 [34] |
| lme4breeding | GPL ≥ 2 [37] | R | lme4 avec Z·L ; rotation par décomposition propre [38] | Creux | us documenté dans le préprint [38] | Non vérifié | K [38] | Non | 1.1.5, 5 octobre 2026 [37] |
| breedR | GPL-3 pour l'interface R ; binaires BLUPF90 [39] | R + binaires Fortran [39] | REMLF90 (EM), AIREMLF90 (AI), Gibbs [39] | MME creuses | Non vérifié | Oui, AR1⊗AR1 et splines [39] | Pedigree [39] | Non | 0.12-7 sur GitHub ; plus maintenu [39] |
| glmmTMB | AGPL-3 [40] | R + C++ (TMB) [40] | ML ou REML par Laplace, gradient par différentiation automatique [40] | Creux (TMB) | us, rr, diag, cs, toep, etc. [41] | ar1 par terme ; produit de Kronecker non documenté | K dense via propto ou equalto [41] | Non | 1.1.15.2, 29 septembre 2026 [40] |
| spaMM | CeCILL-2 [42] | R + C++ [42] | ML ou REML par Laplace, h-vraisemblance [42] | Creux pour les champs markoviens et AR [42] | Multivarié (fitmv) [42] ; us non vérifié | AR1, Matérn ; Kronecker non vérifié | K (corrMatrix) ou précision [42] | Non | 4.7.17, 9 octobre 2026 [42] |
| MCMCglmm | GPL ≥ 2 [43] | R + C++ | Bayésien MCMC [44] | A⁻¹ creuse [45] | us, idh, multivarié [44] | Non | K⁻¹ via ginverse [45] | Non | 2.36, 6 mai 2024 [43] |
| brms | GPL-2 [46] | R + Stan | HMC (NUTS) [47] | Dense | Multivarié [46] | Non natif (non vérifié) | K via gr(cov = ) | OpenCL via cmdstanr [48] | 2.23.0, 9 septembre 2025 [46] |
| R-INLA | MIT (dépôt) [49] | C (GMRFLib) + R [49] | INLA, pas de REML [50] | Précision creuse [50] | Non vérifié | ar1 et produit par groupe [51] | K⁻¹ via generic0 [51] | Non vérifié | 26.10.26, 26 octobre 2026 [49] |
| MixedModels.jl | MIT [52] | Julia [52] | PLS + optimiseur sans dérivées (NEWUOA par défaut en v5) ; REML en option [52] | Creux par blocs | Non | Non | Non natif | Non documenté | v5.9.1, 2 octobre 2026 [52] |
| statsmodels MixedLM | BSD-3 [53] | Python [53] | ML ou REML selon Lindstrom et Bates [54][55] | Dense par groupe [54] | Non | Non | Non | Non | 0.15.0, 27 août 2026 [53] |

### 2. Fiches par logiciel

#### 2.1 ASReml (ASReml-R 4.2 et ASReml-SA 4.2)

ASReml est un logiciel commercial de VSN International. La licence académique commence à 375 £ ; les tarifs commerciaux sont sur devis [1]. ASReml-R tourne sous la même licence que le programme autonome ASReml-SA [4]. Le système de licence repose sur Reprise, avec gestion en nuage [4]. Je n'ai pas trouvé de licence gratuite pour les pays à faible revenu sur la page officielle.

Le moteur combine l'algorithme AI de Gilmour, Thompson et Cullis (1995) et des méthodes matricielles creuses [2][3]. Le manuel de ASReml-R précise que seuls les éléments non nuls sont stockés [2]. Les termes sont partagés entre un ensemble dense et un ensemble creux [2]. L'inverse de la matrice des coefficients est formée en entier pour l'ensemble dense, et seulement en partie pour l'ensemble creux [2]. Les effets aléatoires sont creux par défaut. Une matrice de parenté génomique dense peut être basculée dans l'ensemble dense avec l'option dense [2]. Les équations creuses sont réordonnées pour préserver la creusité [2].

Les limites d'échelle sont documentées. La version 4.2 porte la mémoire accessible à 96 Go et le parallélisme à 16 threads [2][4]. L'espace de travail par défaut de ASReml-SA est de 2 Go [4]. Le manuel de ASReml-R indique que certains très gros problèmes exigent ASReml-SA à cause des surcoûts de R [2]. La version 4.2.0.393 du 29 août 2025 a porté l'espace de travail de predict à 32 Go [5]. Les structures de variance couvrent us, fa, facv, xfa, rr, ar1 et leurs produits directs [4]. La dernière version de ASReml-R est la 4.2.0.480-2 du 6 mai 2026, compilée pour R 4.6 [5]. Aucun GPU n'est mentionné.

#### 2.2 Echidna

Echidna est développé par Arthur Gilmour, auteur principal de ASReml [6]. Il est gratuit hors usage commercial. Un usage commercial exige une licence ASReml en cours [6]. Le téléchargement demande une inscription et l'acceptation d'une licence [6]. Le code source n'est pas distribué. Le site présente le logiciel comme une version alpha mise à jour au moins tous les mois [6]. Il tourne sous Mac, Linux et Windows 64 bits [6].

L'estimation des paramètres de variance repose sur l'algorithme AI [6]. La note du build 2.15 (22 octobre 2026) annonce la formation d'une matrice H⁻¹ pour 201 150 individus en pedigree et 40 000 en GRM en 10 minutes au lieu de 9 heures. Chaque itération avec cette matrice prend 24 minutes [7]. Je n'ai pas pu lire le manuel ni l'article de présentation au WCGALP 2018 (lien mort). Les structures disponibles et le langage restent donc non vérifiés.

#### 2.3 WOMBAT

WOMBAT est développé par Karin Meyer (AGBU, University of New England) [8][56]. Il est gratuit pour la communauté scientifique, avec obligation de citation [8]. Seul l'exécutable est distribué, avec le manuel et des exemples, pour Linux et Windows [8][56].

Par défaut, WOMBAT fait quelques itérations PX-EM puis passe à l'AI-REML [8]. Il offre aussi l'EM classique et deux maximisations sans dérivées, Powell et Simplex [8]. Les calculs sont creux par défaut, avec une approche supernodale qui appelle BLAS et LAPACK sur les sous-blocs denses [8]. L'ordonnancement se choisit entre MMD, AMD et METIS [8]. Une option --dense existe pour les GRM denses. Le stockage compact limite alors le système à 65 535 équations ; l'option --dense2 lève cette limite [8]. Les matrices sont non structurées par défaut, avec estimation possible en rang réduit par un facteur de Cholesky partiel [8]. Une matrice de covariance connue s'entre sous forme inverse, dans un fichier .gin [8]. Le manuel de 2017 ne documente pas de structure autorégressive spatiale. Un résumé de 2020 mentionne le REML pénalisé, les GRM, les scans d'association et le single-step [57]. La version courante n'a pas pu être vérifiée.

#### 2.4 Famille BLUPF90

Les programmes sont développés à l'Université de Géorgie par I. Misztal, D. Lourenco et collaborateurs [12]. Ils sont gratuits pour la recherche, avec citation obligatoire. Un usage commercial passe par un contact avec les auteurs [9]. La version de recherche limite le nombre d'animaux génotypés à 30 000 [11]. Les binaires sont fournis pour Linux, Windows et Mac. Seules d'anciennes versions du code source sont disponibles « à but éducatif » [10]. Les programmes courants ne sont donc pas des logiciels libres. Certains programmes (BLUP90IOD3, ACCF90GS2) ne sont distribués que sous accord [11].

BLUPF90+ regroupe blupf90, remlf90 (EM-REML) et airemlf90 (AI-REML). GIBBSF90+ regroupe les échantillonneurs de Gibbs [11]. Le passage de FSPAK au paquet supernodal YAMS a accéléré l'AI-REML d'un facteur moyen supérieur à 10, sur 23 modèles [13]. Le calcul des traces a gagné jusqu'à un facteur 67 grâce au stockage en table de hachage [13]. Le single-step (H⁻¹) est la spécialité de la famille [14]. Le wiki signale le passage progressif du compilateur ifort à ifx en 2025-2026 [10]. Les structures spatiales AR1⊗AR1 ne sont pas natives ; breedR les construit dans sa couche R avant d'appeler les binaires [39].

#### 2.5 DMU

DMU est développé par Per Madsen et Just Jensen (QGG, Aarhus University). Le serveur officiel refusait les connexions pendant la recherche. D'après les extraits de moteur de recherche du guide officiel, le paquet est distribué sous forme d'exécutables, gratuit pour la recherche, et l'usage commercial demande un accord avec QGG [15]. Le module DMUAI estime les composantes par AI-REML selon Jensen et al. (1997), avec repli possible sur l'EM [15]. Le guide le plus récent repéré est la version 6, release 5.6, de juin 2024 [15]. L'article de HIBLUP indique que DMU utilise FSPAK et que la version distribuée n'offre pas de calcul parallèle [16]. Ces points restent à confirmer sur le guide lui-même.

#### 2.6 MiX99 et MiXBLUP

MiX99 est un logiciel de Luke (Institut finlandais des ressources naturelles) [17]. Trois licences existent : Standard, Core pour les pays émergents et University pour l'usage académique [17]. Les licences sont réservées aux applications d'élevage animal [17]. Le solveur mix99s est multithread ; mix99p et apax99p sont des programmes MPI [17]. Les références citées portent sur le gradient conjugué préconditionné et sur un REML EM par Monte Carlo. Une référence de 2024 porte sur l'AI-REML [17]. La page ne donne pas de tarif.

MiXBLUP est développé par Wageningen University & Research avec Luke. Son cœur dérive de MiX99 [19]. Les licences coûtent de 2 400 €/an (pedigree, moins d'un million d'équations) à 24 790 €/an (licence « Large EC », 3 à 10 postes) [18]. Un essai d'un mois est gratuit. Aucune licence académique n'est mentionnée [18]. La méthode REML n'est pas précisée sur ces pages.

#### 2.7 MTG2

MTG2 est développé par S. Hong Lee [21]. Il combine l'AI-REML « direct » de Lee et van der Werf (2006) avec une décomposition propre de la GRM [21][58]. L'algorithme direct travaille sur la matrice V des observations et non sur les MME [21]. Il est avantageux quand la covariance génétique est dense [21]. La décomposition propre exige une seule GRM de rang plein et un plan équilibré, sans phénotype manquant entre caractères [21]. Le dépôt Dryad du 17 août 2020 contient la version 2.16, les sources de la 2.15 en Fortran et le manuel [20]. La licence logicielle n'est pas précisée.

#### 2.8 GCTA et LDAK

GCTA est sous GPL-3, avec des parties en LGPL [22]. La dernière étiquette du dépôt est la v1.94.1 (2022), et le dernier commit date de juillet 2025 [22]. GREML repose sur l'AI-REML direct sur V dense [23][24]. L'article de HIBLUP indique que GCTA ne traite que deux caractères à la fois et n'a pas de fonction pour les effets d'environnement aléatoires [16]. fastGWA-REML exploite une GRM creuse pour les grandes cohortes humaines [25]. LDAK est sous GPL-3 (version 3 ou ultérieure) [26]. Il estime l'héritabilité par REML sur des noyaux pondérés [27]. Ces deux logiciels visent la génétique humaine. Ils n'offrent ni structures us/fa générales ni structures spatiales.

#### 2.9 sommer

sommer est le paquet R libre de G. Covarrubias-Pazaran, sous GPL-3 [28]. La version 4.5.12 est parue sur CRAN le 25 octobre 2026 [28]. Le paquet a beaucoup changé depuis 2025.

L'historique des fonctions se lit dans le CHANGELOG [28]. mmer est la fonction historique à inversion directe sur la matrice r × r des observations. mmec, ajoutée en 4.1.7 (juillet 2022), implantait l'AI de Henderson selon Jensen, Madsen et Thompson (1997). En 4.4.1 (avril 2025), les deux moteurs ont été unifiés dans mmes. En 4.4.81 (septembre 2026), mmer est revenue à côté de mmes. Dans la 4.5.12, mmes et mmer sont exportées. mmec ne l'est plus. Le code de mmer est regroupé dans un fichier DEPRECATED_mmer.R, « en attente de dépréciation » [28].

Le moteur mmes résout par défaut les MME avec un AI-REML creux [28]. La factorisation LDLᵀ creuse est réutilisée pour les solutions et les dérivées. Les éléments utiles de l'inverse sont obtenus par la récursion de Takahashi [28]. CHOLMOD est disponible, avec ordonnancement METIS optionnel ou AMD par défaut [28]. La formulation directe (formulation = "direct") travaille sur la covariance dense des observations [28]. Les deux moteurs disposent d'une récupération EM, GEM ou PX-GEM quand l'AI échoue [28]. La convergence exige un score projeté petit ; un pas refusé renvoie convergence = FALSE [28]. Un solveur PCG avec estimateur stochastique des traces (Hutchinson) existe [28]. Le mode solveOnly résout les MME à variances connues sans former C. Le CHANGELOG annonce environ 10 s pour un million d'animaux en un caractère et 40 s en trois caractères [28].

Les structures se combinent en produits de Kronecker arbitraires par vsm() [28]. La vignette spatiale montre un champ AR1⊗AR1 en effet aléatoire ou en résidu [28]. Les matrices K ou K⁻¹ s'entrent par Gu, l'inverse étant marquée par un attribut [28]. rrm() ajuste depuis 2026 un vrai rang réduit à variances spécifiques nulles, comme rr() de ASReml [28]. Une version compilée avec CUDA accélère certains noyaux denses du moteur Henderson. Le travail creux reste sur CPU, et la formulation directe reste sur CPU [28]. Le rythme de changement est très rapide : la 4.4.90 a cassé l'interface de contrôle numérique [28].

#### 2.10 lme4 et ses extensions génétiques

lme4 est sous GPL ≥ 2 [29]. La version 2.1-0 est parue le 23 octobre 2026 [29]. Le moteur résout un problème de moindres carrés pénalisés par Cholesky creux, et optimise la vraisemblance profilée avec un optimiseur sans dérivées [30]. Depuis la version 2.0-1 (mars 2026), lme4 accepte des structures us, diag, cs et ar1 sur la covariance d'un terme aléatoire [31]. lme4 n'accepte pas de matrice de parenté native.

Les extensions génétiques remplacent Z par Z·L, où L est un facteur de la matrice de parenté [33][38]. lme4qtl suit ce principe pour des matrices K quelconques [33]. Son dernier commit date de janvier 2021 [32]. pedigreemm le fait pour le pedigree [35][36] ; sa dernière version date d'octobre 2024 [34]. lme4breeding, du même auteur que sommer, ajoute une rotation par décomposition propre limitée à une seule matrice de parenté [38]. Sa version 1.1.5 date du 5 octobre 2026 [37]. Aucun de ces paquets ne propose de modèle multivarié us natif au sens de ASReml.

#### 2.11 breedR

breedR est une interface R sous GPL-3 qui appelle les binaires BLUPF90 (REMLF90, AIREMLF90, Gibbs) [39]. Il vise les essais forestiers. Il propose des effets spatiaux autorégressifs et des splines, ainsi qu'un effet de compétition entre voisins [39]. Le README indique que le paquet n'est plus maintenu activement [39]. Il n'est pas sur CRAN.

#### 2.12 glmmTMB et spaMM

glmmTMB est sous AGPL-3 [40]. Il ajuste par maximum de vraisemblance avec intégration de Laplace et gradients par différentiation automatique (TMB) [40]. L'option REML = TRUE existe [41]. Les structures incluent us, toep, cs, diag, ar1, hetar1, ou, exp, gau, mat, rr, propto et equalto [41]. propto ajuste une covariance proportionnelle à une matrice connue, donc une parenté K dense [41]. Aucun produit de Kronecker ar1⊗ar1 n'est documenté dans la vignette.

spaMM est sous CeCILL-2 [42]. Il ajuste des modèles avec effets spatialement corrélés, réponses multivariées et effets non gaussiens [42]. La vraisemblance ou la vraisemblance restreinte utilise des approximations de Laplace [42]. Les champs markoviens et les modèles autorégressifs exploitent des matrices de précision creuses [42].

#### 2.13 Voies bayésiennes : MCMCglmm, brms, R-INLA

MCMCglmm est sous GPL ≥ 2 [43]. Il ajuste des GLMM multivariés par MCMC [44]. L'argument ginverse reçoit des inverses creuses proportionnelles à la covariance d'un effet aléatoire, typiquement A⁻¹ [45]. Sa dernière version date de mai 2024 [43].

brms est sous GPL-2 et s'appuie sur Stan (HMC, NUTS) [46][47]. Il ajuste des modèles multivariés. Une covariance connue s'entre par gr(cov = ). Le calcul OpenCL est disponible avec le moteur cmdstanr, pas avec rstan [48].

R-INLA est sous licence MIT dans son dépôt public [49]. Il implémente l'approximation de Laplace imbriquée intégrée et non le REML [50]. Le modèle latent generic0 reçoit une matrice de précision fournie, donc K⁻¹ [51]. Le modèle ar1 et la construction par groupe donnent des structures séparables [51]. La version du 26 octobre 2026 est la 26.10.26 [49].

#### 2.14 MixedModels.jl et statsmodels

MixedModels.jl est sous licence MIT [52]. La version 5.9.1 date du 2 octobre 2026 [52]. Depuis la v5, l'optimisation est non contrainte et l'optimiseur par défaut est NEWUOA de NLopt [52]. Le paquet suit l'approche PLS de lme4. Il ne prend ni parenté ni structure us multivariée native.

statsmodels est sous BSD-3 [53]. MixedLM est fondé sur des groupes indépendants [54]. Les effets croisés exigent de traiter tout le jeu de données comme un seul groupe [54]. Les calculs de vraisemblance, de gradient et de hessien suivent Lindstrom et Bates (1988) [54][55]. Ce module ne convient pas aux modèles génétiques avec parenté.

### 3. Comparaisons publiées de vitesse et de justesse

Toutes les comparaisons repérées ont été conduites par les auteurs de l'un des logiciels comparés. Aucune comparaison indépendante récente n'inclut à la fois ASReml, BLUPF90, WOMBAT et DMU sur de grandes données génomiques.

| Étude | Logiciels | Données | Résultat chiffré | Réserve |
|---|---|---|---|---|
| Covarrubias-Pazaran 2016 [59] | sommer, ASReml, rrBLUP, regress, EMMREML, BGLR, MCMCglmm | GBLUP 1 composante, N = 5 000 ; 3 composantes denses, N = 10 585 hybrides | 1 composante : sommer 240 s, ASReml 443 s, EMMREML 611 s, rrBLUP 766 s, BGLR 1 292 s, regress 1 684 s, MCMCglmm 533 556 s. 3 composantes : sommer 426 s, regress 1 868 s, EMMREML 1 137 s, ASReml 6 874 s, BGLR 17 364 s, MCMCglmm plus de 12 jours | ASReml a signalé une terminaison anormale dans les deux cas. L'inversion de la matrice 10 585 × 10 585 exigée par ginverse domine le temps de ASReml et MCMCglmm. Machine i7 3,4 GHz, 16 Go |
| Lee et van der Werf 2016 [21] | MTG2, GEMMA, ASReml, WOMBAT | Souris, n = 1 908, jusqu'à 5 caractères ; humains ARIC, n = 7 263, 8 caractères | Souris : MTG2 quelques secondes, « quelques milliers de fois » plus rapide que ASReml et WOMBAT. ARIC : ASReml et WOMBAT trop lents pour tourner. Estimations « essentiellement identiques » | Une seule GRM, plan équilibré par imputation. Gain annoncé supérieur à 1 000 fois |
| Masuda et al. 2015 [13] | AIREMLF90 avec FSPAK contre YAMS | 23 modèles (poulet, porc, bovins viande et lait), dont single-step | YAMS en moyenne plus de 10 fois plus rapide que FSPAK ; calcul des traces jusqu'à 67 fois plus rapide | Comparaison interne au même logiciel |
| Masuda et al. 2014, WCGALP [60] | AIREMLF90 avec FSPAK contre YAMS | Plus de 15 000 génotypés | 1 caractère : 5 premières itérations en plus de 4 h (FSPAK) contre 20 min (YAMS). 4 caractères : échec de FSPAK, 5 h avec YAMS | Chiffres lus dans un extrait de recherche ; PDF inaccessible, non vérifié |
| Yin et al. 2023 (HIBLUP) [16] | HIBLUP, BLUPF90, ASReml-SA, GCTA, MTG2, LDAK, sommer | Simulation de 10 000 individus ; 500 000 individus pour GCTA | Estimations quasi identiques. ASReml converge en 9 itérations contre 15 pour BLUPF90, pour un temps total voisin. sommer exige beaucoup de temps et de mémoire. BLUPF90 est avantagé quand peu d'individus phénotypés sont génotypés. GRM de 500 000 individus : 21 h pour HIBLUP contre 56 h pour GCTA | Taille limitée par la licence de BLUPF90. DMU exclu faute de parallélisme. Chiffres détaillés dans des tableaux non relus |
| Covarrubias-Pazaran 2024 (lme4breeding, préprint) [38] | lme4breeding, lme4qtl 0.2.2, sommer 4.3.4, ASReml 4.1.0.182 | Blé, 599 lignées, 4 environnements | Composantes « presque identiques ». Gain de 4 à 160 fois pour le modèle diagonal. Rapports de 1:1 350 contre ASReml et 1:1 261 contre lme4qtl pour le modèle non structuré | Mac M1, 8 Go. Préprint non publié. Rapports cités sans les temps sous-jacents |
| Caamal-Pat et al. 2021 (lme4GS) [61] | lme4GS, BGLR 0.8, sommer 4.1.3 | Modèles génomiques à noyaux | Modèle 6 : sommer 46 s, lme4GS 81,5 s, BGLR 143 s. Modèle 7 : lme4GS 1,3 s, sommer 2,7 s, BGLR 20,2 s | Pas de comparaison avec ASReml |
| Covarrubias-Pazaran 2018 (préprint) [62] | sommer, ASReml-R, MCMCglmm | Jeux BTdata et CPdata | D'après un extrait : ASReml-R le plus rapide sur BTdata, sommer le plus lent ; sommer plus rapide sur CPdata ; sommer lent quand p < n | Chiffres non relus sur le texte complet |

Deux annonces de développeurs complètent ce tableau sans constituer des comparaisons. Echidna annonce 24 minutes par itération avec une H⁻¹ de 201 150 individus et 40 000 génotypés [7]. sommer annonce environ 10 s pour un million d'animaux en mode solveOnly, sans estimation des variances [28].

Un poster de Rogers et Taylor (IBS Australasie, 2019) comparerait des paquets R de modèles mixtes pour l'amélioration des plantes. Je n'ai pas trouvé son contenu. Les chiffres qui lui sont attribués en ligne restent non vérifiés.

### 4. Sources

1. VSN International, page ASReml-R. https://vsni.co.uk/software/asreml-r (consultée le 26/10/2026).
2. Butler D.G., Cullis B.R., Gilmour A.R., Gogel B.J., Thompson R. (2023). ASReml-R Reference Manual Version 4.2. VSN International. https://asreml.kb.vsni.co.uk/wp-content/uploads/sites/3/ASReml-R-Reference-Manual-4.2.pdf (sections 1.1, 2.2, A.1, A.2).
3. Gilmour A.R., Thompson R., Cullis B.R. (1995). Average information REML: an efficient algorithm for variance parameter estimation in linear mixed models. Biometrics 51:1440–1450. https://doi.org/10.2307/2533274
4. Gilmour A.R., Gogel B.J., Cullis B.R., Welham S.J., Thompson R. (2022). ASReml User Guide Release 4.2 Functional Specification. https://asreml.kb.vsni.co.uk/wp-content/uploads/sites/3/ASReml-4.2-Functional-Specification.pdf (sections 7.11, 11.3.6, option !MP).
5. VSN International, Release notes ASReml-R 4.2. https://asreml.kb.vsni.co.uk/?p=2374
6. Echidna Mixed Models Software, page d'accueil. https://www.echidnamms.org/
7. Gilmour A.R., note « Echidna Build 2.15 October 22, 2026 ». https://www.echidnamms.org/echidna-build-2-15-october-22-2026/
8. Meyer K. (2017). WOMBAT User Notes, version du 31 mai 2017. Copie INRAE : https://forge-dga.jouy.inra.fr/attachments/download/8063/WombatManual-V20170531.pdf (chapitres 1, 2, 5.3). Site officiel http://didgeridoo.une.edu.au/km/wombat.php inaccessible le 26/10/2026.
9. BLUPF90 wiki, « Conditions of use ». https://nce.ads.uga.edu/wiki/doku.php?id=condition_of_use
10. BLUPF90 wiki, « Distribution ». https://nce.ads.uga.edu/wiki/doku.php?id=distribution
11. BLUPF90 wiki, « Application programs ». https://nce.ads.uga.edu/wiki/doku.php?id=application_programs
12. BLUPF90 wiki, page d'accueil. https://nce.ads.uga.edu/wiki/doku.php?id=start
13. Masuda Y., Aguilar I., Tsuruta S., Misztal I. (2015). Technical note: Acceleration of sparse operations for average-information REML analyses with supernodal methods and sparse-storage refinements. J. Anim. Sci. 93:4670–4674. https://doi.org/10.2527/jas.2015-9395
14. Lourenco D., Legarra A., Tsuruta S., Masuda Y., Aguilar I., Misztal I. (2020). Single-step genomic evaluations from theory to practice: using SNP chips and sequence data in BLUPF90. Genes 11:790. https://doi.org/10.3390/genes11070790
15. Madsen P., Jensen J. DMU, a package for analyzing multivariate mixed models, guide v6 release 5.6 (2024). https://dmu.ghpc.au.dk/dmu/DMU/Doc/Current/dmuv6_guide.5.6.pdf. Serveur inaccessible le 26/10/2026 ; contenu connu par extraits de moteur de recherche seulement.
16. Yin L. et al. (2023). HIBLUP: an integration of statistical models on the BLUP framework for efficient genetic evaluation using big genomic data. Nucleic Acids Res. 51:3501–3512. https://doi.org/10.1093/nar/gkad074
17. Luke, page MiX99. https://www.luke.fi/en/node/8911
18. MiXBLUP, page « Licenses ». https://www.mixblup.eu/licenses
19. MiXBLUP, page d'accueil. https://www.mixblup.eu/
20. Lee S.H., Zhou X. (2020). CORE GREML: MTG2 code and instruction. Dryad. https://doi.org/10.5061/dryad.bk3j9kd8c
21. Lee S.H., van der Werf J.H.J. (2016). MTG2: an efficient algorithm for multivariate linear mixed model analysis based on genomic information. Bioinformatics 32:1420–1422. https://doi.org/10.1093/bioinformatics/btw012
22. Dépôt GCTA. https://github.com/jianyangqt/gcta (README, étiquettes, commits).
23. Yang J., Lee S.H., Goddard M.E., Visscher P.M. (2011). GCTA: a tool for genome-wide complex trait analysis. Am. J. Hum. Genet. 88:76–82. https://doi.org/10.1016/j.ajhg.2010.11.011
24. Voir [21], introduction : l'AI direct est implanté dans GCTA-GREML.
25. Jiang L. et al. (2019). A resource-efficient tool for mixed model association analysis of large-scale data. Nat. Genet. 51:1749–1755. https://doi.org/10.1038/s41588-019-0530-8
26. Dépôt LDAK, fichier LICENCE.txt. https://github.com/dougspeed/LDAK
27. Speed D., Hemani G., Johnson M.R., Balding D.J. (2012). Improved heritability estimation from genome-wide SNPs. Am. J. Hum. Genet. 91:1011–1021. https://doi.org/10.1016/j.ajhg.2012.10.010
28. sommer 4.5.12, CRAN. https://cran.r-project.org/package=sommer. Sources lues : DESCRIPTION, NAMESPACE, CHANGELOG, README, vignettes sommer.qg et sommer.spatial, R/FUN_mme_gpu.R, R/FUN_mme_control.R (archive sommer_4.5.12.tar.gz).
29. lme4 2.1-0, CRAN. https://cran.r-project.org/package=lme4
30. Bates D., Mächler M., Bolker B., Walker S. (2015). Fitting linear mixed-effects models using lme4. J. Stat. Softw. 67(1). https://doi.org/10.18637/jss.v067.i01
31. lme4, fichier NEWS, versions 2.0-1 et 2.1-0. https://github.com/lme4/lme4/blob/master/inst/NEWS.Rd
32. Dépôt lme4qtl. https://github.com/variani/lme4qtl
33. Ziyatdinov A. et al. (2018). lme4qtl: linear mixed models with flexible covariance structure for genetic studies of related individuals. BMC Bioinformatics 19:68. https://doi.org/10.1186/s12859-018-2057-x
34. pedigreemm 0.3-5, CRAN. https://cran.r-project.org/package=pedigreemm
35. Vazquez A.I., Bates D.M., Rosa G.J.M., Gianola D., Weigel K.A. (2010). Technical note: An R package for fitting generalized linear mixed models in animal breeding. J. Anim. Sci. 88:497–504. https://doi.org/10.2527/jas.2009-1952
36. Voir [59], tableau 2 : pedigreemm accepte un pedigree mais pas une matrice de covariance fournie directement.
37. lme4breeding 1.1.5, CRAN. https://cran.r-project.org/package=lme4breeding
38. Covarrubias-Pazaran G. (2024). lme4breeding: enabling genetic evaluation in the era of genomic data. bioRxiv. https://doi.org/10.1101/2024.05.17.594795
39. Dépôt breedR (README, DESCRIPTION, NAMESPACE, R/competition.R, R/spatial.R). https://github.com/famuvie/breedR
40. glmmTMB 1.1.15.2, CRAN. https://cran.r-project.org/package=glmmTMB ; Brooks M.E. et al. (2017). glmmTMB balances speed and flexibility among packages for zero-inflated generalized linear mixed modeling. R J. 9:378. https://doi.org/10.32614/RJ-2017-066
41. glmmTMB, vignette « Covariance structures with glmmTMB ». https://github.com/glmmTMB/glmmTMB/blob/master/glmmTMB/vignettes/covstruct.rmd
42. spaMM 4.7.17, CRAN. https://cran.r-project.org/package=spaMM ; Rousset F., Ferdy J.-B. (2014). Testing environmental and genetic effects in the presence of spatial autocorrelation. Ecography 37:781–790. https://doi.org/10.1111/ecog.00566
43. MCMCglmm 2.36, CRAN. https://cran.r-project.org/package=MCMCglmm
44. Hadfield J.D. (2010). MCMC methods for multi-response generalized linear mixed models: the MCMCglmm R package. J. Stat. Softw. 33(2). https://doi.org/10.18637/jss.v033.i02
45. MCMCglmm, page d'aide MCMCglmm.Rd, argument ginverse. https://github.com/jarrodhadfield/MCMCglmm
46. brms 2.23.0, CRAN. https://cran.r-project.org/package=brms
47. Bürkner P.-C. (2017). brms: an R package for Bayesian multilevel models using Stan. J. Stat. Softw. 80(1). https://doi.org/10.18637/jss.v080.i01
48. brms, fichier R/backends.R (gestion d'OpenCL selon le moteur). https://github.com/paul-buerkner/brms
49. Dépôt R-INLA (LICENSE, version 26.10.26). https://github.com/hrue/r-inla
50. Rue H., Martino S., Chopin N. (2009). Approximate Bayesian inference for latent Gaussian models by using integrated nested Laplace approximations. J. R. Stat. Soc. B 71:319–392. https://doi.org/10.1111/j.1467-9868.2008.00700.x
51. Documentation R-INLA des modèles latents generic0 et ar1. https://github.com/hrue/r-inla/tree/devel/r-inla.org/doc/latent
52. Dépôt MixedModels.jl (README, version 5.9.1). https://github.com/JuliaStats/MixedModels.jl
53. Dépôt statsmodels (licence, version 0.15.0). https://github.com/statsmodels/statsmodels
54. statsmodels, documentation « Linear Mixed Effects Models ». https://www.statsmodels.org/stable/mixed_linear.html
55. Lindstrom M.J., Bates D.M. (1988). Newton-Raphson and EM algorithms for linear mixed-effects models for repeated-measures data. J. Am. Stat. Assoc. 83:1014–1022. https://doi.org/10.1080/01621459.1988.10478693
56. Meyer K. (2007). WOMBAT: a tool for mixed model analyses in quantitative genetics by restricted maximum likelihood (REML). J. Zhejiang Univ. Sci. B 8:815–821. https://doi.org/10.1631/jzus.2007.B0815
57. Meyer K. (2020). WOMBAT: a tool for estimation of genetic parameters, highlights and updates. ICQG 6, résumés, p. 199. https://hdl.handle.net/1959.11/30741
58. Lee S.H., van der Werf J.H.J. (2006). An efficient variance component approach implementing an average information REML suitable for combined LD and linkage mapping with a general complex pedigree. Genet. Sel. Evol. 38:25. https://doi.org/10.1186/1297-9686-38-1-25
59. Covarrubias-Pazaran G. (2016). Genome-assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11:e0156744. https://doi.org/10.1371/journal.pone.0156744
60. Masuda Y. et al. (2014). Communication au 10e WCGALP sur YAMS dans AIREMLF90. Copie INIA : https://inia.uy/sites/default/files/publications/2024-10/Masuda-et-al2014WCGALP.pdf (accès refusé, non vérifié).
61. Caamal-Pat D. et al. (2021). lme4GS: an R-package for genomic selection. Front. Genet. 12:680569. https://doi.org/10.3389/fgene.2021.680569
62. Covarrubias-Pazaran G. (2018). Software update: moving the R package sommer to multivariate mixed models for genome-assisted prediction. bioRxiv. https://doi.org/10.1101/354639

## Annexe B. Calcul accéléré pour le REML et les modèles mixtes linéaires, 2010-2026

Consultation le 26 octobre 2026. Les DOI ont été vérifiés sur Crossref, ou sur DataCite pour arXiv. L'activité des dépôts vient de l'API GitHub. Pour pyreml, le code source 0.5.3 a été lu en entier. Le mot « non vérifié » signale un point que je n'ai pas pu confirmer sur une source primaire.

### Tableau comparatif

Légende des colonnes. « Calcul » dit si l'outil estime les variances par REML exact, par une approximation, ou ne fait que du BLUP ou de la GWAS. « Modèle » résume la généralité. « Échelle » reprend le plus grand cas publié et son matériel.

| Outil | Calcul | Modèle | Technique calcul et mémoire | Échelle et matériel | Licence, maintenance |
|---|---|---|---|---|---|
| pyreml 0.5.3 (Marchal et al.) | REML exact, gradient analytique, L-BFGS puis Adam en secours. Pas d'AI. | Multi-termes, multivarié, Σ⊗K, fa, us, diag, blocs, hétéroscédasticité, noyaux spatiaux, régression aléatoire. | PyTorch dense. Forme V ou Woodbury (C = G⁻¹ + Z'R⁻¹Z dense q×q). Grilles d'index au lieu de Z. FP32 puis FP64. | n = 10⁶, q = 25 000, 41 paramètres, 262 à 279 s sur A100 80 Go. | GPL-3.0. Créé en juin 2026, 19 versions, 0.5.3 le 21/10/2026. |
| kmmr (Melsen et al.) | AI-REML exact. | Multivarié, us, FA étendu, AR1×AR1, noyaux gaussiens, Kronecker, empilement. | Partie creuse par CHOLMOD supernodal sur CPU, inverse partielle de Takahashi. Complément de Schur dense sur GPU (cuSOLVER). | GBLUP 10 375 génotypes : 49 s sur A100 contre 989 s pour ASReml-R. | GPL-3.0, R. Actif (push le 22/10/2026). |
| SAIGE-GPU | AI-REML sur quasi-vraisemblance pénalisée (biaisé selon les auteurs de SAIGE), puis tests. | Un seul terme génétique (GRM pleine ou creuse). | PCG sans former la GRM, produits génotype × vecteur par cuBLAS, MPI multi-nœuds, traces aléatoires. | 635 969 individus, 2 068 phénotypes. Summit (V100), Frontier (MI250X), A100. | GPL-3.0. Commit en avril 2026. |
| BOLT-REML | REML approché : gradient Monte Carlo, AI, gradient conjugué. | Plusieurs composantes et plusieurs traits. | Mémoire ≈ MN/4 octets. CPU seulement. | ≈ 50 000 individus dans l'article. UK Biobank en quelques jours sur 8 fils (manuel). | GPLv3. Version 2.5 en juin 2025. |
| DISSECT | REML exact sur V dense distribuée. | LMM univarié et bivarié. Multi-composantes annoncé en test en 2015. | MPI, BLACS, ScaLAPACK. | 470 000 individus, 590 004 SNP, moins de 4 h sur 8 400 cœurs (ARCHER). | GPLv3. Maintenance non vérifiée. |
| REACTA | REML exact, héritabilité régionale. | Dérivé de GCTA. | MAGMA et CUDA, régions réparties par MPI. | 1 024 GPU de Titan, 15 min. | Licence non vérifiée. |
| GPU_REML (dépôt Asiandier/SPLIT) | REML approché : PCG par blocs, Hutchinson, SLQ, AI ou Fisher. | Plusieurs GRM, une variance par GRM. | GRM jamais formée, génotypes compressés lus en flux vers le GPU. | Non publié. | Aucune licence. Créé en mai 2026. Aucun article. |
| TorchLIMIX | LMM multivarié de LIMIX pour la GWAS. Mode d'estimation non vérifié. | Deux traits ou plus, une parenté. | Algèbre par lots sur GPU (PyTorch). | Gain jusqu'à 100× face à LIMIX. | Dépôt annoncé en 404. |
| Quickdraws | Spike-and-slab par inférence variationnelle stochastique. Héritabilité par RHE, pas par REML. | Un trait à la fois. | Mini-lots sur GPU, génotypes sur 2 bits. | 405 088 individus sur un A10G de 16 Go. | GPL-3.0. Commit en août 2025. |
| Faster-LMM-D | LMM à une parenté pour la GWAS, AI-REML repris de GEMMA. | Une composante. | GPU seulement pour le produit matriciel dense. | Non vérifié. | GPL-3.0. Dernier commit en 2018. |
| FaST-LMM, pylmm, GEMMA | REML exact par décomposition spectrale. | Une parenté. mvLMM dans GEMMA. | CPU. | 120 000 individus (FaST-LMM). | Apache-2.0, MIT, GPL-3.0. Pas de GPU. |
| Hail LinearMixedModel | REML exact. | Une matrice K. | Rotation spectrale. Seuls les tests sont distribués (Spark). | Non vérifié. | MIT. Retiré depuis Hail 0.2.94. |
| REGENIE, GloWGR | Ridge empilé, pas de REML. | GWAS. | Segments en mémoire. Spark pour GloWGR. | 407 746 individus. | MIT, Apache-2.0. Pas de GPU. |
| tensorQTL, jaxQTL | Modèles linéaires ou GLM, pas d'effet aléatoire. | QTL. | PyTorch GPU ; JAX. | Gain > 200× (tensorQTL). | BSD-3, MIT. |
| GPyTorch (BBMM) | Vraisemblance marginale de PG, approchée par CG et Lanczos stochastique au-delà de 800 points. Pas de REML trouvé. | Processus gaussiens, noyaux structurés. | Produits matrice-matrice par lots, préconditionneur, multi-GPU. | 10⁶ points exacts en moins de 2 h (Wang et al. 2019). | MIT. Actif. |
| KeOps | Réductions de noyaux sans stocker la matrice n×n. | Brique de calcul. | Code CUDA généré, autodiff. | Millions de points en secondes (résumé). | MIT. Actif. |
| GPflow, GPJax, tinygp | Vraisemblance de PG, variationnel. | Processus gaussiens. | TensorFlow ou JAX. | Non vérifié. | Apache-2.0, MIT. |
| TensorFlow Probability | ML par EM Monte Carlo (HMC + Adam) dans l'exemple LMM. Pas de REML. | 3 facteurs aléatoires iid (InstEval). | XLA sur GPU. | 58 737 observations. | Apache-2.0. Actif. |
| NumPyro, BlackJAX | Bayésien (NUTS, HMC, VI). | Tout modèle écrit par l'utilisateur. | JAX, chaînes vectorisées, GPU et TPU. | Non applicable. | Apache-2.0. Actifs. |
| Stan avec OpenCL | Bayésien. | Tout modèle Stan. | Cholesky, sa dérivée et quelques vraisemblances sur GPU. | Non vérifié. | BSD-3. Actif. |
| ExaGeoStat | Vraisemblance gaussienne exacte (Matérn), plus approximations TLR. | Spatial. | Cholesky dense en tuiles, GPU et mémoire distribuée. | Non vérifié dans le texte. | Non vérifié. |
| INLA-DIST, DALIA | Bayésien approché (INLA). | Champs spatio-temporels, multivariés pour DALIA. | Factorisation creuse par blocs sur GPU, inversion sélective (Serinv). | 496 GH200 (Alps) pour DALIA. | Non vérifié. |
| jaxmg | Brique : Cholesky, solve et log-déterminant distribués. | Matrice dense. | cuSOLVERMp, multi-GPU et multi-nœuds, appelable sous jit. | Jusqu'à la mémoire cumulée des GPU. | Apache-2.0. Actif. |
| miraculix dans MiXBLUP 3.1 | BLUP par PCG (single-step), pas de REML. | Multicaractère. | Produit génotype × vecteur sur GPU, 2 bits, CUTLASS. | 2,61 M génotypés, 26,46 M au pedigree, A100. | Apache 2.0 annoncé. MiXBLUP commercial. |
| BOLT (Theta Solutions) | BLUP par PCG et Gibbs single-step. Pas de REML documenté. | Bayes A, B, C, LASSO. | Plusieurs GPU. | 13 M animaux. | Propriétaire. |
| MC-REML (Matilainen ; Bermann) | EM-REML et AI-REML à traces Monte Carlo ou Hutchinson. | Multicaractère, single-step. | CPU, solveurs itératifs. | 7 M animaux, 330 000 génotypés, 121 h, 53 Go. | Dans BLUPF90 (2025). |
| SLEMM, SL_REML | REML par Lanczos stochastique. | Deux composantes (SL_REML). SLEMM avec poids de SNP. | CPU. | 3 M individus simulés (SLEMM). | Code public. |
| RHE-mc, REHE | Méthode des moments (Haseman-Elston), pas REML. | Jusqu'à 100 composantes (RHE-mc). | Traces aléatoires. | 10⁶ individus, 10⁶ SNP en quelques heures. | Non vérifié. |
| HIBLUP, MTG2 | AI-REML exact (forme MME ou V). | Multicaractère. | CPU multithread. | 500 000 individus simulés (HIBLUP). | Exécutables gratuits. Pas de GPU. |

### Solveurs REML généraux sur GPU

#### pyreml

- Article. Marchal, Trottier, Bry et Raimondi, prépublication bioRxiv du 24 octobre 2026, doi:10.64898/2026.10.20.757556. Logiciel archivé sous doi:10.5281/zenodo.21394319. Le texte complet a été lu.
- Paquet. PyPI, de la 0.1.0 (24 juin 2026) à la 0.5.3 (21 octobre 2026), 19 versions (https://pypi.org/pypi/pyreml/json). GPL-3.0-or-later. Dépôt sans étoile ni ticket ouvert au 26 octobre (API GitHub). Auteurs à l'IGMM et à l'IMAG, Montpellier.
- Optimisation. Descente de gradient sur −2ℓ_REML, dans la lignée de Simchoni et Rosset (JMLR 2023, 24(156)), avec le terme REML ajouté. L-BFGS de PyTorch avec recherche de Wolfe forte, Adam en secours quand la recherche linéaire échoue (pyreml/Optimizer.py). β est optimisé comme paramètre en même temps que les variances (MixedModel.py, l. 297 et 304). Le critère d'arrêt compte les chiffres significatifs restants de la vraisemblance (doc/training.qmd).
- Gradient. L'autodiff à travers le Cholesky est remplacée par le gradient REML en forme close, injecté par une « ghost loss » (MixedModel.py, REML_loss, l. 709-720). Le recouvrement avec le custom_vjp de remlax est direct.
- Forme V ou Woodbury. La voie SMW factorise C = G⁻¹ + Z'R⁻¹Z, dense de taille q×q, et ne forme jamais V. Elle s'active quand q < n et que R est diagonale (MixedModel.py, l. 207-225 ; Variance.py, classe Capacitance). Sinon un Cholesky de V dense (classe DirectSolve).
- Mémoire. Les blocs de Kronecker et ZGZ' ne sont jamais formés. Z et R sont remplacés par des grilles d'index. Un test vérifie que la plus grande allocation reste sous 4 × max(q², np) (tests/test_biggplup.py, n = 23 600).
- Précision. Trois modes : FP64, FP32, ou mixte (FP32 puis FP64), mixte par défaut. En FP32, 241 tests fonctionnels échouent (prépublication, section 2.3).
- Absences vérifiées dans le code 0.5.3. Aucun appel à torch.sparse : tout est dense. Aucun algorithme AI ni matrice d'information. Aucune erreur-type des composantes de variance dans la documentation : seuls les effets fixes (SE) et les BLUP (PEV, CD) en reçoivent, par l'inverse dense des équations de Henderson (MixedModel.py, l. 817-830). Pas de GLMM (prépublication, discussion). Pas de multi-GPU. Pas de test du rapport de vraisemblance intégré.
- Structures. Huit structures à gauche (iid, diag, full, fa, bl_resp, bl_form, kr_resp, kr_form) et sept à droite (iid, str, het, dist, eucl, ar_iso, ar_ani), GaussianComponents.py, l. 10-11. Un constructeur de bas niveau accepte des fonctions varmeth arbitraires.
- Bancs publiés. GBLUP comparé à WOMBAT (EM, AI, PX-AI sur 64 cœurs EPYC 9124) et non à ASReml. Plus de 12× plus rapide sur A100, plus de 2× sur L40S. WOMBAT AI n'a pas convergé pour n ≥ 10 000. Cas extrême : 25 000 génotypes × 40 environnements, 10⁶ observations, 279 s en FP64 sur A100 (prépublication) ; la documentation affiche 261,65 s pour un modèle de même taille avec 41 paramètres (doc/building_models.qmd). FA2 sur 8 traits et 942 individus : 429 s pour pyreml contre 281 s pour WOMBAT, ou 66 s à tolérance égale. Validation contre lme4, rrBLUP, BLUPF90 et spaMM (505 tests).
- Les auteurs se disent « the first library performing REML estimations on GPU ». REACTA (2014) et SAIGE-GPU contredisent cette priorité au sens strict.

#### kmmr

- Melsen, van Eeuwijk et Peeters, arXiv:2610.16671 (octobre 2026), doi:10.48550/arXiv.2610.16671. Dépôt github.com/KillianMelsen/kmmr, GPL-3.0, R, archive doi:10.5281/zenodo.23398636. Paquet qualifié d'expérimental, Linux ou WSL2 avec CUDA.
- AI-REML multivarié sur les équations du modèle mixte. Les équations creuses sont éliminées par Cholesky supernodal (CHOLMOD). Le complément de Schur dense passe sur GPU : potrf, potrs, potri de cuSOLVER, SpMM de cuSPARSE. L'inverse partielle de la partie creuse vient d'une récursion de Takahashi supernodale (cholmod-extra).
- Bancs sur A100 80 Go et i7. GBLUP blé à 10 375 génotypes : ASReml-R 988,6 s, sommer 511,5 s, kmmr 129,3 s sur 8 fils et 49,2 s sur GPU. Modèle FA + noyau empilé : ASReml-R ne converge pas en 50 itérations, kmmr 6,8 s sur GPU.
- Limites écrites par les auteurs. La mémoire du complément de Schur croît comme le carré du nombre d'équations denses. La forme V sur GPU est citée comme extension naturelle. Sur des modèles en une étape, kmmr converge vers des log-vraisemblances plus basses qu'ASReml-R (−3948,66 contre −3937,65 en FA1).

#### GPU_REML (dépôt Asiandier/SPLIT)

- Dépôt créé le 30 mai 2026, sans licence ni article (https://github.com/Asiandier/SPLIT).
- REML pour V = Σθ_g K_g + θ_e I, avec K_g = Z_g Z_g'/m appliqué sans former la GRM. Génotypes compressés lus en flux vers le GPU. PCG par blocs, traces de Hutchinson, log-déterminants par SLQ, mises à jour AI ou Fisher contraintes (README).
- Outil de recherche sur la partition de l'héritabilité. Aucun banc publié.

#### TorchLIMIX

- Horn, Nikoloski et Lippert, bioRxiv 2026, doi:10.64898/2026.04.29.721342. Portage PyTorch du LMM multivarié de LIMIX, gain jusqu'à deux ordres de grandeur.
- Le dépôt github.com/bi-horn/torchLIMIX renvoie 404. Le mode d'estimation des covariances n'est pas vérifié.

### Génétique humaine et GWAS

#### SAIGE-GPU

- Rodriguez et al. 2026, Bioinformatics 42(3):btag032, doi:10.1093/bioinformatics/btag032.
- Le GPU accélère l'étape 1, l'ajustement du GLMM nul. Les produits génotype × vecteur du PCG passent par cuBLAS sans former la GRM. Les colonnes de génotypes sont réparties entre nœuds par MPI. L'étape 2 reste sur CPU.
- L'étape 1 de SAIGE est un AI-REML sur quasi-vraisemblance pénalisée, avec traces aléatoires. Zhou et al. 2018 (doi:10.1038/s41588-018-0184-y) déconseillent son τ pour estimer l'héritabilité.
- 635 969 participants du MVP et 2 068 phénotypes. Environ 72× par itération PCG, 3 à 5× sur l'ajustement complet. Summit, Frontier, A100 en nuage.
- Code : branche SAIGE-GPU-1.3.3 de github.com/saigegit/SAIGE, GPL-3.0, archive doi:10.5281/zenodo.17642591.

#### Quickdraws

- Loya et al. 2025, Nat Genet, doi:10.1038/s41588-024-02044-7. Régression spike-and-slab par inférence variationnelle stochastique, en mini-lots sur GPU (PyTorch).
- L'héritabilité vient du Haseman-Elston aléatoire, pas du REML. Un trait à la fois. 405 088 individus sur un A10G de 16 Go. github.com/PalamaraLab/quickdraws, GPL-3.0.

#### Faster-LMM-D, pylmm, FaST-LMM, GEMMA

- Faster-LMM-D (github.com/genetics-statistics/faster_lmm_d, GPL-3.0, dernier commit en décembre 2018). Le GPU (cuBLAS ou ArrayFire) ne remplace que le produit matriciel dense (source/faster_lmm_d/optmatrix.d). Aucun article associé.
- pylmm (github.com/nickFurlotte/pylmm). Une composante, grille sur h puis Brent, en REML. Réécrit en octobre 2026 (v1.0.0, MIT). Pas de GPU.
- FaST-LMM (Lippert et al. 2011, doi:10.1038/nmeth.1681) et GEMMA (Zhou et Stephens 2012, doi:10.1038/ng.2310 ; mvLMM 2014, doi:10.1038/nmeth.2848). Décomposition spectrale d'une seule parenté. Pas de GPU.

#### BOLT-LMM et BOLT-REML

- BOLT-REML : Loh et al. 2015, Nat Genet 47:1385, doi:10.1038/ng.3431 (PMC4666835). Le gradient est estimé par Monte Carlo, la hessienne par la matrice AI, les systèmes par gradient conjugué. 100 vecteurs aléatoires en général, 15 pour les analyses régionales.
- Précision. Estimations presque identiques à GCTA avec 100 essais sur 50 000 traits simulés. Avec 15 essais, les erreurs-types augmentent d'un facteur 1,03.
- Coût ≈ O(MN^1,5) et mémoire ≈ MN/4 octets, presque indépendante du nombre de composantes. 22 composantes chromosomiques et paires bivariées dans l'article. CPU seulement (manuel v2.5, https://alkesgroup.broadinstitute.org/BOLT-LMM/BOLT-LMM_manual.html).

#### GCTA, fastGWA, REGENIE, Hail, tensorQTL, jaxQTL

- GCTA (Yang et al. 2011, doi:10.1016/j.ajhg.2010.11.011) et fastGWA (Jiang et al. 2019, doi:10.1038/s41588-019-0530-8) : CPU, MKL ou OpenBLAS. Aucune voie GPU trouvée.
- REGENIE (Mbatchou et al. 2021, doi:10.1038/s41588-021-00870-7) : ridge empilé, pas de REML, pas de GPU. GloWGR le distribue sur Spark (github.com/projectglow/glow).
- Hail LinearMixedModel : une seule matrice K, REML exact sur le nœud principal, tests distribués. Retiré depuis la version 0.2.94 (documentation Hail).
- tensorQTL (Taylor-Weiner et al. 2019, doi:10.1186/s13059-019-1836-7) et jaxQTL (Zhang et al. 2025, medRxiv, doi:10.1101/2025.01.18.25320755) : pas d'effet aléatoire.

### REML distribué

#### DISSECT

- Canela-Xandri et al. 2015, Nat Commun 6:10162, doi:10.1038/ncomms10162. C++ avec MPI, BLACS et ScaLAPACK.
- REML exact en O(N³) sur V dense distribuée. LMM univarié et bivarié. 470 000 individus et 590 004 SNP en moins de 4 h sur 8 400 cœurs et environ 16 To de RAM (ARCHER), REML et BLUP compris.
- GPLv3 sur www.dissect.ed.ac.uk. Certificat expiré lors de la consultation. Maintenance non vérifiée.

#### REACTA

- Cebamanos et al. 2014, Bioinformatics, doi:10.1093/bioinformatics/btt754. Dérivé de GCTA pour l'héritabilité régionale.
- GRM et inversion REML par MAGMA, noyaux CUDA, régions réparties par MPI. 19 746 régions sur 1 024 GPU de Titan (ORNL) en 15 min 26 s. Environ 3× sur une K20.

#### ExaGeoStat

- Abdulah et al. 2018, IEEE TPDS 29(12), doi:10.1109/TPDS.2018.2850749 (arXiv:1708.02835). Vraisemblance gaussienne exacte d'un champ de Matérn par Cholesky dense en tuiles, sur GPU NVIDIA, KNL et mémoire distribuée. Variantes approchées en faible rang par tuiles (non vérifié dans ce texte).
- Domaine spatial, sans effets fixes ni REML générique.

#### Autres briques distribuées

- jaxmg (github.com/flatironinstitute/jaxmg, Apache-2.0, actif). Relie JAX à cuSOLVERMp. Fournit potrs avec log-déterminant, sur matrices réparties sur plusieurs GPU et nœuds, sous jit.
- Lewis et al. 2022, PNAS, doi:10.1073/pnas.2122762119. Algèbre linéaire dense distribuée sur pods TPU, en float32. Produit de deux matrices de taille 2²⁰ en environ 2 min sur 2 048 cœurs TPU v3.
- Wang et al. 2019, « Exact Gaussian processes on a million data points », NeurIPS (arXiv:1903.08114). Produits de noyau partitionnés et distribués sur plusieurs GPU, gradient conjugué. Plus d'un million de points en moins de 2 h. Nombre et type de GPU non vérifiés dans le texte.

### Bibliothèques d'apprentissage automatique et probabilistes

#### GPyTorch et BBMM

- Gardner et al. 2018, NeurIPS, arXiv:1809.11165. Gradient conjugué par lots modifié (mBCG) qui fournit en un appel les solves, le log-déterminant par Lanczos stochastique et les traces. Coût exact ramené de O(n³) à O(n²). Préconditionneur spécialisé.
- Réglages par défaut (https://docs.gpytorch.ai/en/stable/settings.html) : Cholesky sous 800 points, 10 vecteurs pour les traces, 20 itérations de Lanczos. Au-delà, la vraisemblance est donc approchée.
- Objectif : vraisemblance marginale de processus gaussiens. Aucun terme REML (log|X'V⁻¹X|) trouvé dans la documentation (non vérifié dans le code). github.com/cornellius-gp/gpytorch, MIT, actif ; l'algèbre est isolée dans linear_operator.

#### KeOps

- Charlier et al. 2021, JMLR 22(74):1-6 (https://www.jmlr.org/papers/v22/20-275.html). Tenseurs symboliques définis par formule, réduits sur GPU sans stocker la matrice n×n, avec autodiff. Liaisons Python, Matlab et R.
- Brique de calcul pour les noyaux denses de grande taille. github.com/getkeops/keops, MIT.

#### GPflow, GPJax, tinygp

- GPflow (Matthews et al. 2017, JMLR 18(40)) sur TensorFlow, Apache-2.0. GPJax (Pinder et Dodd 2022, JOSS, doi:10.21105/joss.04455) et tinygp (github.com/dfm/tinygp, MIT) sur JAX.
- Processus gaussiens et inférence variationnelle. Pas de modèle mixte de type asreml.

#### TensorFlow Probability

- Dillon et al. 2017, arXiv:1711.10604. L'exemple officiel de modèle mixte (https://www.tensorflow.org/probability/examples/Linear_Mixed_Effects_Models) ajuste InstEval avec trois facteurs aléatoires iid.
- Ajustement par EM Monte Carlo : HMC pour l'étape E, Adam pour l'étape M. Maximum de vraisemblance, pas de REML. GPU recommandé.

#### NumPyro, BlackJAX, Stan

- NumPyro (Phan et al. 2019, arXiv:1912.11554) : NUTS itératif compilé de bout en bout par JAX. BlackJAX (Cabezas et al. 2024, arXiv:2402.10797) : échantillonneurs composables sur CPU, GPU et TPU. Apache-2.0, actifs.
- Stan (Češnovar et al. 2019, arXiv:1907.01063) : routines OpenCL pour le Cholesky, sa dérivée et quelques vraisemblances. Activé par STAN_OPENCL (https://mc-stan.org/docs/cmdstan-guide/parallelization.html).
- Ces outils ne fournissent pas de modèle mixte prêt à l'emploi. La voie bayésienne de remlax recoupe leur usage générique.

#### Autres paquets JAX ou PyTorch

- jaxgam (github.com/scatt67/jaxgam) : portage JAX de mgcv, REML de Laplace pour les modèles additifs. Projet personnel déclaré non prêt pour la production.
- Recherche GitHub du 26 octobre 2026 (« mixed model REML jax », « REML gpu ») : aucune autre bibliothèque JAX de modèles mixtes linéaires généraux trouvée.

### REML stochastique et itératif

#### Estimateurs de trace

- Hutchinson 1989, Commun. Stat. Simul. Comput., doi:10.1080/03610918908812806. Avec m vecteurs de Rademacher, tr(A) ≈ moyenne de z'Az. La variance vaut 2(‖A‖²_F − Σ a_ii²)/m.
- Lanczos stochastique (SLQ) pour tr(f(A)), dont log|A| : Ubaru, Chen et Saad 2017, SIAM J. Matrix Anal. Appl., doi:10.1137/16M1104974, avec bornes sur le nombre de vecteurs.

#### SL_REML (Border et Becker 2019)

- BMC Bioinformatics 20:411, doi:10.1186/s12859-019-2978-z (PMC6668092). Deux algorithmes, SLDF_REML et L_FOMC_REML, fondés sur l'invariance par décalage des sous-espaces de Krylov.
- Une seule passe itérative, puis l'objectif se réévalue par opérations vectorielles. Limité au modèle à deux composantes.
- Précision. Héritabilités cohérentes avec GCTA. Pour n = 250 000, erreur quadratique moyenne de 1,24 × 10⁻⁷ contre 1,74 × 10⁻⁶ pour BOLT-LMM v2.3.2. Code Python sur github.com/rborder/SL_REML.

#### SLEMM

- Cheng et al. 2023, Bioinformatics, doi:10.1093/bioinformatics/btad127. REML par Lanczos stochastique avec pondération des SNP. Jusqu'à 3 M individus et 1 M SNP simulés. CPU. github.com/jiang18/slemm.

#### REML Monte Carlo en sélection animale

- Matilainen et al. 2012, JABG, doi:10.1111/j.1439-0388.2012.01000.x (MC EM-REML multicaractère) et 2013, PLoS One, doi:10.1371/journal.pone.0080821 (MC Newton, MC AI, MC Broyden).
- Bermann et al. 2025, GSE 57:63, doi:10.1186/s12711-025-01006-9. MC-ss-GREML avec EM et AI. Mêmes estimations que l'exact, en 14 % du temps et 1 % de la mémoire (100 000 animaux, 14 paramètres). 7 M animaux au pedigree dont 330 000 génotypés en 121 h et 53 Go, sur CPU.
- Bermann et al. 2025, GSE 57:44, doi:10.1186/s12711-025-00991-1. Log-déterminant approché, écart relatif moyen d'environ 10⁻³, 2 à 500 fois plus rapide. Le maximum de la vraisemblance approchée coïncide avec l'exact.
- Strandén et al. 2024, GSE 56:73, doi:10.1186/s12711-024-00939-x. AI-REML augmenté, une seule résolution des équations augmentées par itération, temps par itération réduit de 75 à 86 %.

#### Méthodes des moments et MM

- RHE-reg (Wu et Sankararaman 2018, doi:10.1093/bioinformatics/bty253) et RHE-mc (Pazokitoroudi et al. 2020, doi:10.1038/s41467-020-17576-9). Haseman-Elston à traces aléatoires. Cent composantes sur 10⁶ individus et 10⁶ SNP en quelques heures. Pas de REML.
- REHE (Yue et al. 2021, Genet Epidemiol 45(8), doi:10.1002/gepi.22432). Haseman-Elston contraint à des variances positives, avec une variante par rééchantillonnage. Précision annoncée comparable au REML.
- Algorithmes MM (Zhou, Hu, Zhou et Lange 2019, JCGS 28(2), doi:10.1080/10618600.2018.1529601). ML et REML exacts par minorisation-maximisation. Convergence globale vers un point KKT, plus rapide que l'EM au-delà de deux composantes. S'étend au multivarié avec données manquantes. Pas de GPU mentionné.

### Sélection animale et végétale

#### miraculix et MiXBLUP

- Freudenberg et al. 2023, Front. Genet. 14:1220408, doi:10.3389/fgene.2023.1220408. Produits Z·M et Z'·M avec Z stockée sur 2 bits, sur GPU par extension de CUTLASS 2.10 (V100, A100, H100).
- Données ICBF : 2,61 M génotypés, 26,46 M animaux au pedigree, 16,59 M enregistrements sur 6 caractères. ssSNPBLUP passe de 6,47 h à 2,49 h sur A100. C'est du BLUP par PCG, sans REML.
- Intégré au solveur hpblup de MiXBLUP 3.1 (Vandenplas et al., WCGALP 2022, doi:10.3920/978-90-8686-940-4_358). MiXBLUP est commercial. Dépôt github.com/alexfreudenberg/miraculix, Apache 2.0 annoncé, dernier push en août 2023.

#### BOLT de Theta Solutions et échantillonneurs de Gibbs

- Fernando et al. 2016, GSE 48:96, doi:10.1186/s12711-016-0273-2. Single-step bayésien et PCG sur quatre Titan X. 6,18 M animaux au pedigree, environ 3 h par analyse unicaractère.
- Présentation WCGALP 2018 (https://www.icar.org/Documents/Auckland-2018/WCGALP/1415%20Dorian%20Garrick.pdf) : 13 M animaux, 1 000 itérations PCG en 1,5 h. Logiciel propriétaire. Aucun REML documenté.
- GenSel hétérogène CPU/GPU (DeVore et al. 2014, IEEE BIBM, doi:10.1109/bibm.2014.6999186) : gain de 1,84× seulement, la chaîne restant séquentielle.

#### Travaux GPU pionniers

- Coffey et al. 2011, Interbull Bull. 44:114-116 (https://journal.interbull.org/index.php/ib/article/view/1694) : produit AA' sur GTX 590, de 507 s à 70 s.
- Meyer et Tier 2013, Proc. AAABG 20:270-273 (https://hdl.handle.net/1959.11/14323) : construction et inversion de G sur GPU, gain de 4 à 5×.

#### Absences constatées

- Aucune voie GPU documentée dans BLUPF90, DMU, MiX99, JWAS (Project.toml de la v2.4.0 sans CUDA) ni HIBLUP (Yin et al. 2023, doi:10.1093/nar/gkad074, texte complet sans « GPU »). Une absence de résultat ne prouve pas une absence.
- MTG2 (Lee et van der Werf 2016, doi:10.1093/bioinformatics/btw012) : AI-REML direct sur V avec décomposition propre. Aucune version GPU trouvée.
- ChromaX (Younis et al. 2023, doi:10.1093/bioinformatics/btad691) simule des programmes de sélection sur GPU en JAX, sans inférence.

### Briques de calcul pour la forme MME sur GPU

- Cholesky creux supernodal de CHOLMOD accéléré sur GPU : Rennich, Stosic et Davis 2016, Parallel Computing, doi:10.1016/j.parco.2016.06.004.
- Inversion sélective distribuée : PSelInv (Jacquelin, Lin et Yang 2016, ACM TOMS, doi:10.1145/2786977). Inversion sélective sur GPU de matrices tridiagonales par blocs avec flèche : Serinv (Maillou et al. 2025, IEEE CLUSTER, doi:10.1109/cluster59342.2025.11186484), utilisé par DALIA (Gaedke-Merzhäuser et al. 2025, arXiv:2507.06938, 496 GH200) et issu de INLA-DIST (SIAM J. Sci. Comput. 2024, doi:10.1137/23m1561531).
- cuDSS de NVIDIA (https://docs.nvidia.com/cuda/cudss/latest/index.html) : solveur direct creux sur GPU, matrices définies positives, multi-GPU et multi-nœuds. La documentation ne mentionne ni inverse sélective ni log-déterminant.
- jax.experimental.sparse.linalg.spsolve : solveur creux par QR, CUDA seulement, sans vmap (https://docs.jax.dev/en/latest/_autosummary/jax.experimental.sparse.linalg.spsolve.html).

### Positionnement de remlax

Recoupements.
- pyreml partage l'essentiel de la voie dense de remlax : REML exact sur GPU, gradient analytique au lieu de l'autodiff à travers Cholesky, Woodbury, Kronecker jamais formé, fa et multivarié.
- kmmr partage la voie MME : AI-REML, Cholesky supernodal, inverse partielle de Takahashi, GPU pour la partie dense.
- La voie bayésienne recoupe l'usage générique de NumPyro, BlackJAX ou Stan.

Points propres à remlax dans les sources consultées.
- Un même moteur réunit la forme V dense, la forme MME creuse et l'AI-REML. pyreml n'a ni creux ni AI. kmmr n'a pas de forme V, qu'il cite comme travail futur.
- La détection automatique des termes globaux et des blocs, avec V = A + UU' par Woodbury et blocs traités en lots (vmap), n'a d'équivalent dans aucun outil trouvé. Le Woodbury de pyreml factorise un C dense de taille q.
- L'ajustement de modèles multivariés lourds sur une tranche MIG de 10 Go n'est rapporté nulle part ailleurs.
- La matrice AI fournit des erreurs-types et sert aux tests de covariance. pyreml n'en donne aucune pour les variances.

Capacités présentes ailleurs et absentes de remlax, d'après la description reçue.
- Le multi-GPU et le multi-nœud : DISSECT, REACTA, SAIGE-GPU, ExaGeoStat, DALIA, Wang et al. 2019. jaxmg offre en JAX un Cholesky distribué avec log-déterminant.
- La voie sans matrice, par PCG et traces stochastiques : BOLT-REML, SAIGE, SLEMM, MC-ss-GREML, GPyTorch, GPU_REML. Elle seule dépasse quelques centaines de milliers d'individus.
- La précision mixte FP32 puis FP64 de pyreml, utile sur les GPU grand public (L40S).
- Les génotypes compressés sur 2 bits et lus en flux vers le GPU : miraculix, SAIGE-GPU, Quickdraws, GPU_REML.

### Points non vérifiés

- Le texte complet de kmmr au-delà de 100 000 caractères n'a pas été lu. Le nombre d'observations de ses bancs n'est pas donné par les auteurs.
- L'écart entre 261,65 s (documentation pyreml, version 0.5.0) et 279 s (prépublication) pour le cas à 10⁶ observations. Les deux configurations ne sont peut-être pas identiques.
- Le mode d'estimation et le code de TorchLIMIX.
- Une étape 2 de SAIGE-GPU sur GPU AMD, que seul un nom de branche suggère.
- La maintenance de DISSECT, la licence de REACTA et d'ExaGeoStat, le matériel de Wang et al. 2019.
- L'absence de REML dans le code de GPyTorch (vérifiée seulement dans la documentation).
- Le contenu GPU éventuel des actes WCGALP 2022 sur BLUPF90, MiX99 et JWAS : PDF bloqués par Cloudflare. EAAP, ASAS-ADSA et Interbull n'ont pas été couverts en texte intégral.
- La priorité revendiquée par pyreml (« first library performing REML estimations on GPU ») est contredite par REACTA (2014) pour un REML spécialisé. Pour un REML généraliste, la priorité entre pyreml (juin 2026) et kmmr (dépôt de mars 2026, prépublication d'octobre 2026) reste ouverte.

### Sources

Articles et prépublications (DOI vérifiés)
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

Dépôts, paquets et documentation
- pyreml : https://github.com/ai-for-genome-interpretation-lab/pyreml ; https://pypi.org/project/pyreml/ ; https://ai-for-genome-interpretation-lab.github.io/pyreml/ ; Zenodo doi:10.5281/zenodo.21394319
- kmmr : https://github.com/KillianMelsen/kmmr ; Zenodo doi:10.5281/zenodo.23398636
- SAIGE-GPU : https://github.com/saigegit/SAIGE (branche SAIGE-GPU-1.3.3) ; Zenodo doi:10.5281/zenodo.17642591
- GPU_REML : https://github.com/Asiandier/SPLIT
- Faster-LMM-D : https://github.com/genetics-statistics/faster_lmm_d
- pylmm : https://github.com/nickFurlotte/pylmm
- FaST-LMM : https://github.com/fastlmm/FaST-LMM
- GEMMA : https://github.com/genetics-statistics/GEMMA
- BOLT-LMM, manuel v2.5 : https://alkesgroup.broadinstitute.org/BOLT-LMM/BOLT-LMM_manual.html
- GCTA : https://github.com/jianyangqt/gcta
- REGENIE : https://github.com/rgcgithub/regenie ; Glow : https://github.com/projectglow/glow
- Quickdraws : https://github.com/PalamaraLab/quickdraws
- tensorQTL : https://github.com/broadinstitute/tensorqtl ; jaxQTL : https://github.com/mancusolab/jaxqtl
- Hail : https://github.com/hail-is/hail
- SL_REML : https://github.com/rborder/SL_REML ; SLEMM : https://github.com/jiang18/slemm
- GPyTorch : https://github.com/cornellius-gp/gpytorch ; réglages : https://docs.gpytorch.ai/en/stable/settings.html
- KeOps : https://github.com/getkeops/keops
- GPflow : https://github.com/GPflow/GPflow ; tinygp : https://github.com/dfm/tinygp
- TensorFlow Probability, exemple LMM : https://www.tensorflow.org/probability/examples/Linear_Mixed_Effects_Models
- NumPyro : https://github.com/pyro-ppl/numpyro ; BlackJAX : https://github.com/blackjax-devs/blackjax
- Stan OpenCL : https://mc-stan.org/docs/cmdstan-guide/parallelization.html
- jaxmg : https://github.com/flatironinstitute/jaxmg
- jaxgam : https://github.com/scatt67/jaxgam
- cuDSS : https://docs.nvidia.com/cuda/cudss/latest/index.html
- JAX spsolve : https://docs.jax.dev/en/latest/_autosummary/jax.experimental.sparse.linalg.spsolve.html
- miraculix : https://github.com/alexfreudenberg/miraculix
- BOLT (Theta Solutions), WCGALP 2018 : https://www.icar.org/Documents/Auckland-2018/WCGALP/1415%20Dorian%20Garrick.pdf
- JWAS : https://github.com/reworkhow/JWAS.jl

## Annexe C. Antériorité des mécanismes de remlax

La plupart des briques de remlax ont une antériorité nette en génétique quantitative. Le gradient REML par rapport à V, l'AI-REML, la factorisation creuse des équations du modèle mixte et l'inversion sélective ont été publiés entre 1977 et 2015. La dissection emboîtée et la vraisemblance tirée de la Cholesky d'une matrice bordée datent de la même période. Ils sont aussi implémentés dans ASReml, WOMBAT, DFREML/MTDFREML et la famille BLUPF90/YAMS. La détection automatique des termes à traiter en rang faible et des blocs d'observations indépendants, par un modèle de coût et une recherche gloutonne, n'a pas d'antériorité trouvée. L'ensemble actif sur le cône semi-défini positif en échelle des variances, suivi d'un finisseur à région de confiance, n'en a pas non plus dans les logiciels REML examinés. L'estimation conjointe par REML d'une portée continue de noyau IGE avec les composantes de variance reste elle aussi sans antécédent trouvé. L'exécution sur GPU par dérivation automatique a une antériorité contemporaine (pyreml, 2026) et une antériorité ancienne pour la dérivation automatique elle-même (Meyer et Smith 1996).

### Tableau de synthèse

| Mécanisme de remlax | Verdict | Références principales |
|---|---|---|
| 1. -2logL REML depuis V dense, gradient P − Pyy'P | Antériorité claire pour la formule et la forme V dense | [Harville 1977](https://doi.org/10.1080/01621459.1977.10480998), [Searle et al. 1992](https://doi.org/10.1002/9780470316856), [Covarrubias-Pazaran 2016](https://doi.org/10.1371/journal.pone.0156744) |
| 1b. Gradient par rapport à V codé en règle adjointe (custom_vjp), reste par dérivation automatique | Antériorité partielle | [Smith 1995](https://doi.org/10.1080/10618600.1995.10474671), [Giles 2008](https://doi.org/10.1007/978-3-540-68942-3_4), [Kristensen et al. 2016](https://doi.org/10.18637/jss.v070.i05), pyreml 0.5.3 ([Marchal et Raimondi 2026](https://doi.org/10.5281/zenodo.21394319)) |
| 2a. Équations du modèle mixte creuses, Cholesky supernodale | Antériorité claire | [Masuda et al. 2014](https://doi.org/10.1111/jbg.12058), [Masuda et al. 2015](https://doi.org/10.2527/jas.2015-9395), manuel WOMBAT 2017, manuel ASReml 4.2 |
| 2b. Inversion sélective (Takahashi) pour le gradient et l'AI | Antériorité claire | [Erisman et Tinney 1975](https://doi.org/10.1145/360680.360704), [Misztal et Perez-Enciso 1993](https://doi.org/10.3168/jds.S0022-0302%2893%2977478-0), [Masuda et al. 2014](https://doi.org/10.1111/jbg.12058), [Lin et al. 2011](https://doi.org/10.1145/1916461.1916464) |
| 2c. Échec d'AMD sur un C localement dense, gain de la dissection emboîtée ou d'un ordre « blocs puis séparateur » | Antériorité partielle (mécanisme et usage connus, constat IGE non publié) | [George 1973](https://doi.org/10.1137/0710032), [Karypis et Kumar 1998](https://doi.org/10.1137/S1064827595287997), [Amestoy et al. 2004](https://doi.org/10.1145/1024074.1024081), [Masuda et al. 2014](https://doi.org/10.1111/jbg.12058), manuels WOMBAT et ASReml |
| 3a. V = A + UU', Woodbury, A diagonale par blocs | Antériorité claire pour chaque élément | [Lindstrom et Bates 1988](https://doi.org/10.1080/01621459.1988.10478693), [Mishchenko et al. 2007](https://doi.org/10.48550/arXiv.0709.0625), [Lippert et al. 2011](https://doi.org/10.1038/nmeth.1681), [Strandén et Garrick 2009](https://doi.org/10.3168/jds.2008-1929) |
| 3b. Cholesky d'une matrice bordée [U X y]'A⁻¹[U X y] + diag(I,0,0) donnant log|V|, log|X'V⁻¹X|, y'Py | Antériorité claire pour l'algèbre ; partielle pour A contenant des termes aléatoires absorbés | [Smith et Graser 1986](https://doi.org/10.3168/jds.S0022-0302%2886%2980516-1), [Graser et al. 1987](https://doi.org/10.2527/jas1987.6451362x), [Meyer 1989](https://doi.org/10.1186/1297-9686-21-3-317), [Bates et DebRoy 2004](https://doi.org/10.1016/j.jmva.2004.04.013), [Bates et al. 2025](https://doi.org/10.48550/arXiv.2505.11674) |
| 3c. Détection automatique des termes globaux et des blocs indépendants (modèle de coût, recherche gloutonne) | Pas d'antériorité trouvée ; heuristiques voisines seulement | manuel ASReml 4.2 (!DENSE, !GDENSE), manuel WOMBAT (--dense), [Kristensen et al. 2016](https://doi.org/10.18637/jss.v070.i05), [Bates et al. 2025](https://doi.org/10.48550/arXiv.2505.11674), pyreml 0.5.3 |
| 4a. AI-REML en échelle des variances | Antériorité claire | [Gilmour et al. 1995](https://doi.org/10.2307/2533274), [Johnson et Thompson 1995](https://doi.org/10.3168/jds.S0022-0302%2895%2976654-1) |
| 4b. Ensemble actif, plancher des variances | Antériorité partielle (bornes scalaires) | [Mishchenko et al. 2008](https://doi.org/10.3233/JCM-2008-81-203), manuel ASReml 4.2 (codes P et B), manuel WOMBAT (--pivot) |
| 4c. Ensemble actif sur le cône semi-défini positif | Pas d'antériorité trouvée dans les logiciels REML examinés | ASReml remplace la mise à jour AI par une mise à jour EM ; WOMBAT reparamètre par Cholesky ([Meyer et Smith 1996](https://doi.org/10.1186/1297-9686-28-1-23)) |
| 4d. Amortissement de Levenberg-Marquardt | Antériorité claire | manuel WOMBAT (--modaim2), [Wolfinger et al. 1994](https://doi.org/10.1137/0915079) |
| 4e. Salves de L-BFGS-B, finisseur Moré-Sorensen avec SR1 | Composants classiques ; usage en REML non trouvé (sauf BFGS chez Mishchenko) | [Byrd et al. 1995](https://doi.org/10.1137/0916069), [Moré et Sorensen 1983](https://doi.org/10.1137/0904038), [Conn et al. 1991](https://doi.org/10.1007/BF01594934), [Mishchenko et al. 2008](https://doi.org/10.3233/JCM-2008-81-203) |
| 4f. Démarrage EM ou PX-EM | Antériorité claire | [Liu et al. 1998](https://doi.org/10.1093/biomet/85.4.755), [Foulley et van Dyk 2000](https://doi.org/10.1186/1297-9686-32-2-143), manuels WOMBAT et ASReml |
| 5. REML sur GPU par dérivation automatique | Antériorité partielle et contemporaine | [Meyer et Smith 1996](https://doi.org/10.1186/1297-9686-28-1-23), [Gardner et al. 2018](https://doi.org/10.48550/arXiv.1809.11165), pyreml 0.5.3, [Horn et al. 2026](https://doi.org/10.64898/2026.04.29.721342) |
| 6a. Modèle IGE de voisinage, dilution | Antériorité claire | [Griffing 1967](https://doi.org/10.1071/BI9670127), [Muir 2005](https://doi.org/10.1534/genetics.104.035956), [Bijma et al. 2007](https://doi.org/10.1534/genetics.106.062711), [Bijma 2010](https://doi.org/10.1534/genetics.110.120485), manuel WOMBAT (SOCIAL) |
| 6b. Noyau à ordre discret ou pondéré par la distance | Antériorité claire | [Cappa et Cantet 2008](https://doi.org/10.1515/sg-2008-0008), [Costa e Silva et al. 2017](https://doi.org/10.1038/hdy.2017.25), [Sato et al. 2021](https://doi.org/10.1038/s41437-020-00401-w) |
| 6c. Portée continue estimée par REML avec les variances | Pas d'antériorité trouvée | profil de vraisemblance de la dilution dans WOMBAT ; balayage d'échelle par PVE chez [Sato et al. 2021](https://doi.org/10.1038/s41437-020-00401-w) |

### Forme V dense et gradient adjoint

La dérivée de la log-vraisemblance REML s'écrit -½ tr(P ∂V/∂θ) + ½ y'P (∂V/∂θ) Py. Cette forme figure dans la revue de [Harville 1977](https://doi.org/10.1080/01621459.1977.10480998) et dans le manuel de [Searle et al. 1992](https://doi.org/10.1002/9780470316856). La dérivée de -2logL par rapport à V vaut donc P − Pyy'P. remlax n'apporte rien sur la formule. L'évaluation REML directe sur V dense est aussi ancienne. Le paquet sommer l'utilise dans sa fonction historique, avec Newton-Raphson et AI ([Covarrubias-Pazaran 2016](https://doi.org/10.1371/journal.pone.0156744)). Les outils de génomique humaine la combinent avec une décomposition propre ([Kang et al. 2008](https://doi.org/10.1534/genetics.107.080101); [Yang et al. 2011](https://doi.org/10.1016/j.ajhg.2010.11.011)).

La dérivation automatique de la vraisemblance REML est elle aussi ancienne. [Smith 1995](https://doi.org/10.1080/10618600.1995.10474671) a dérivé l'algorithme de Cholesky en mode inverse. Murray décrit cette méthode comme l'application manuelle des idées de la dérivation automatique à l'algorithme de Cholesky ([Murray 2016](https://doi.org/10.48550/arXiv.1602.07527)). [Meyer et Smith 1996](https://doi.org/10.1186/1297-9686-28-1-23) ont appliqué cette dérivation aux modèles animaux. Leur résumé indique des dérivées premières et secondes sans inversion de grandes matrices. Les règles adjointes du log-déterminant et de l'inverse sont rassemblées par [Giles 2008](https://doi.org/10.1007/978-3-540-68942-3_4). TMB calcule le gradient du log-déterminant d'un hessien creux par une passe inverse et l'algorithme du sous-ensemble inverse ([Kristensen et al. 2016](https://doi.org/10.18637/jss.v070.i05), section 4.1). Le paquet pyreml 0.5.3 expose une option « analytic_backward » dans une implémentation PyTorch (code source examiné localement, [Marchal et Raimondi 2026](https://doi.org/10.5281/zenodo.21394319)). Nous n'avons pas vérifié le contenu exact de cette règle.

Verdict. La formule et la forme V dense ont une antériorité claire. Une règle adjointe insérée dans un graphe de dérivation automatique a une antériorité partielle. L'apport de remlax tient à l'ingénierie en JAX et à la portabilité CPU/GPU.

### Équations du modèle mixte creuses

#### Absorption et partition creux/dense

La partition des termes en un bloc dense et un bloc creux est documentée dans ASReml. Le manuel 4.2 indique que la partition se fait au point !r, sauf réglage par !DENSE (section 6.10.1). Par défaut, le bloc dense reçoit les effets fixes listés avant les effets aléatoires, jusqu'à 800 équations par défaut. L'option !GDENSE place le premier terme GRM dans les équations denses. L'utilisateur fixe cette partition. ASReml absorbe d'abord les équations creuses, dans un ordre qu'il choisit, puis les équations denses dans l'ordre déclaré (section 6.10.2). WOMBAT propose une option --dense quand une matrice de parenté est dense, par exemple une GRM (manuel 2017, section 5.3.4). Son option --s1step regroupe les individus génotypés en un bloc dense traité à part. Ces partitions vivent toutes dans C. Aucune ne bascule un terme vers une forme V.

#### Ordonnancement

La dissection emboîtée remonte à [George 1973](https://doi.org/10.1137/0710032). METIS en donne une version multiniveau ([Karypis et Kumar 1998](https://doi.org/10.1137/S1064827595287997)). Le degré minimum approché est décrit par [Amestoy et al. 1996](https://doi.org/10.1137/S0895479894278952). [Liu 1985](https://doi.org/10.1145/214392.214398) a proposé le degré minimum multiple.

En sélection animale, ces ordres sont déjà comparés. YAMS emploie l'AMD ou la dissection emboîtée multiniveau ([Masuda et al. 2014](https://doi.org/10.1111/jbg.12058)). Ses auteurs rapportent une factorisation jusqu'à 13 fois plus rapide que FSPAK « when an appropriate ordering strategy was applied ». WOMBAT choisit l'ordre selon la taille du problème (manuel 2017, annexe A.1). Le degré minimum multiple sert jusqu'à 25 000 équations, l'AMD jusqu'à 50 000, METIS au-delà. Ce manuel recommande d'essayer plusieurs ordres et propose --metisall, qui parcourt 48 combinaisons de réglages de METIS. Il cite une étude de Meyer (2005, Proc. AAABG 16:282–285, sans DOI, contenu non vérifié). ASReml expose !EQORDER (manuel 4.2, tableau 5.5). Son manuel décrit un cas proche du constat de remlax. Pour un modèle avec une matrice IBD dense et une parenté creuse, l'ordre déclaré par l'utilisateur (!EQORDER -1) s'exécute plus vite que l'ordre automatique.

Le mécanisme de l'échec d'AMD est documenté. L'AMD de SuiteSparse retire les lignes de plus de 10√n entrées et les place en fin d'ordre (en-tête amd.h ; [Amestoy et al. 2004](https://doi.org/10.1145/1024074.1024081)). Un C localement dense, avec beaucoup de lignes juste sous ce seuil, dégrade donc l'ordre. WOMBAT transmet aussi à METIS un facteur qui désigne les lignes denses.

Verdict. L'usage de METIS pour C en génétique quantitative a une antériorité claire. L'ordre « blocs puis séparateur » correspond au premier niveau d'une dissection emboîtée. Le constat chiffré sur un modèle IGE multivarié n'a pas été trouvé dans la littérature. Il mérite d'être rapporté comme résultat empirique, sans revendication de méthode.

#### Inversion sélective

La récurrence de Takahashi est attribuée à Takahashi, Fagan et Chin (1973, 8th PICA Conference, Minneapolis, pp. 63–69). Cette communication n'a pas de DOI et n'a pas été consultée. [Erisman et Tinney 1975](https://doi.org/10.1145/360680.360704) en donnent la version publiée en revue. [Misztal et Perez-Enciso 1993](https://doi.org/10.3168/jds.S0022-0302%2893%2977478-0) l'ont introduite dans l'EM-REML par inversion creuse ; l'algorithme exact employé n'a pas été vérifié sur le texte. [Johnson et Thompson 1995](https://doi.org/10.3168/jds.S0022-0302%2895%2976654-1) l'associent à l'AI-REML pour les modèles animaux univariés, d'après leur titre. YAMS emploie une factorisation supernodale « left-looking » et une inversion multifrontale inverse ([Masuda et al. 2014](https://doi.org/10.1111/jbg.12058)). Un traitement des traces par table de hachage accélère l'AI-REML jusqu'à 67 fois ([Masuda et al. 2015](https://doi.org/10.2527/jas.2015-9395)). WOMBAT emploie par défaut une approche supernodale et cite Erisman et Tinney (manuel 2017, section 5.3.5). ASReml 4.2 a ajouté une absorption symbolique préalable et une « absorption nodale » (manuel 4.2, section 11.5). Hors génétique, SelInv est l'algorithme supernodal de référence ([Lin et al. 2011](https://doi.org/10.1145/1916461.1916464)).

Verdict. Antériorité claire. remlax réimplémente une technique standard.

### Termes globaux, Woodbury et matrice bordée

#### Woodbury et rang faible

Les équations de Henderson sont elles-mêmes une forme de Woodbury. [Strandén et Garrick 2009](https://doi.org/10.3168/jds.2008-1929) dérivent l'équivalence entre SNP-BLUP et GBLUP par cette identité. [Mishchenko et al. 2007](https://doi.org/10.48550/arXiv.0709.0625) l'appliquent à l'AI-REML avec une représentation de rang faible, exacte ou tronquée, de la matrice IBD. FaST-LMM l'utilise quand la parenté est de rang inférieur à n ([Lippert et al. 2011](https://doi.org/10.1038/nmeth.1681)). Le traitement de V diagonale par blocs de sujets est l'hypothèse de base de [Lindstrom et Bates 1988](https://doi.org/10.1080/01621459.1988.10478693). SAS PROC MIXED déclare ces blocs par SUBJECT= et applique les algorithmes de [Wolfinger et al. 1994](https://doi.org/10.1137/0915079). Le cas « bloc diagonal plus rang faible » apparaît aussi en génomique en une étape ([Liu et al. 2014](https://doi.org/10.3168/jds.2014-7924); [Taskinen et al. 2017](https://doi.org/10.1186/s12711-017-0310-9)).

#### Décomposition propre

La rotation par les vecteurs propres de la parenté rend V diagonale quand un seul terme aléatoire s'ajoute au résidu ([Kang et al. 2008](https://doi.org/10.1534/genetics.107.080101); [Lippert et al. 2011](https://doi.org/10.1038/nmeth.1681); [Zhou et Stephens 2012](https://doi.org/10.1038/ng.2310)). GCTA l'emploie aussi ([Yang et al. 2011](https://doi.org/10.1016/j.ajhg.2010.11.011)). La version multivariée repose sur une transformation canonique conjointe ([Meyer 1985](https://doi.org/10.2307/2530651); [Zhou et Stephens 2014](https://doi.org/10.1038/nmeth.2848); [Lee et van der Werf 2016](https://doi.org/10.1093/bioinformatics/btw012)). Cette astuce ne s'étend pas à plusieurs termes non simultanément diagonalisables. Les modèles IGE en contiennent plusieurs. remlax s'en distingue sur ce point.

#### Vraisemblance bordée

[Smith et Graser 1986](https://doi.org/10.3168/jds.S0022-0302%2886%2980516-1) calculent log|C| et y'Py par élimination de Gauss sur la matrice des équations bordée par la ligne [r' y'R⁻¹y]. Le dernier pivot donne y'Py = y'R⁻¹y − s'r. Cette description provient des notes d'atelier MTDFREML (Boldman et al., 2002, USDA-ARS, sans DOI). Le texte original n'a pas été consulté. [Graser et al. 1987](https://doi.org/10.2527/jas1987.6451362x) et [Meyer 1989](https://doi.org/10.1186/1297-9686-21-3-317) en ont fait la base du REML sans dérivées. [Boldman et Van Vleck 1991](https://doi.org/10.3168/jds.S0022-0302%2891%2978629-3) l'ont accélérée par une factorisation creuse (SPARSPAK). En statistique, [Bates et DebRoy 2004](https://doi.org/10.1016/j.jmva.2004.04.013) puis lme4 ([Bates et al. 2015](https://doi.org/10.18637/jss.v067.i01)) factorisent la matrice augmentée de [ZΛ X y]. [Bates et al. 2025](https://doi.org/10.48550/arXiv.2505.11674) montrent que 2log|R_ZZ| et r²_yy se lisent directement sur la diagonale du facteur, sans résoudre pour β ni pour les modes. Cette construction suppose R = σ²I.

La matrice bordée de remlax, [U X y]'A⁻¹[U X y] + diag(I,0,0), est algébriquement l'équation de Henderson augmentée, avec une matrice A à la place de R. Le cas R diagonale par blocs existe déjà dans les analyses multivariées de DFREML. La nouveauté éventuelle tient au contenu de A. remlax y absorbe des termes aléatoires locaux en forme V et laisse les termes globaux dans U. Cette hybridation V/C n'a pas été trouvée dans les sources consultées.

#### Détection automatique de structure

Aucun logiciel examiné ne choisit seul, par un modèle de coût, les termes à traiter en rang faible et les blocs d'observations indépendants. Les mécanismes voisins sont des heuristiques fixes ou des choix laissés à l'utilisateur. ASReml partitionne au point !r et laisse !DENSE et !GDENSE à l'utilisateur (manuel 4.2). WOMBAT choisit un ordre selon le nombre d'équations et laisse --dense à l'utilisateur (manuel 2017). MixedModels.jl place en tête le terme au plus grand nombre de niveaux, parce que la structure diagonale du bloc [1,1] se conserve dans le facteur ([Bates et al. 2025](https://doi.org/10.48550/arXiv.2505.11674), section 3.3). TMB détecte automatiquement le motif de creux du hessien ([Kristensen et al. 2016](https://doi.org/10.18637/jss.v070.i05), section 4.2). Ce motif sert à la factorisation creuse. Le choix entre formes V et C reste au modélisateur. pyreml 0.5.3 active la voie Woodbury quand le nombre d'effets aléatoires est inférieur à n et que le résidu s'y prête (MixedModel.py, code examiné localement). Cette règle porte sur une seule dimension et ne détecte pas de blocs.

Verdict. Les briques algébriques ont une antériorité claire. La détection automatique par modèle de coût et recherche gloutonne semble nouvelle. Cette conclusion repose sur une recherche non exhaustive. Les moteurs internes d'ASReml sont fermés et peuvent contenir des heuristiques non documentées.

### Optimisation AI-REML et robustesse

L'AI-REML est publié par [Gilmour et al. 1995](https://doi.org/10.2307/2533274) et, pour les modèles animaux creux, par [Johnson et Thompson 1995](https://doi.org/10.3168/jds.S0022-0302%2895%2976654-1). ASReml l'applique en échelle des composantes. Le Newton-Raphson et l'EM pour les modèles mixtes sont décrits par [Lindstrom et Bates 1988](https://doi.org/10.1080/01621459.1988.10478693) et [Jennrich et Schluchter 1986](https://doi.org/10.2307/2530695). [Thompson 2008](https://doi.org/10.1098/rspb.2007.1417) retrace cette histoire.

La reparamétrisation par le facteur de Cholesky, avec log des éléments diagonaux, supprime les contraintes ([Groeneveld 1994](https://doi.org/10.1186/1297-9686-26-6-537); [Pinheiro et Bates 1996](https://doi.org/10.1007/BF00140873); [Meyer et Smith 1996](https://doi.org/10.1186/1297-9686-28-1-23)). WOMBAT l'emploie par défaut et borne les pivots à 0,001 (option --pivot, manuel 2017, section 5.2.5). ASReml agit en échelle des variances. Son code P remplace une variance négative par une petite valeur positive, puis la fixe après des tentatives répétées (code B, manuel 4.2, section 7.7.4). Pour une matrice US, !GP déclenche une mise à jour EM à la place de la mise à jour AI quand celle-ci sort du cône. ASReml fixe aussi la matrice entière quand elle devient quasi singulière (manuel 4.2, !EMFLAG). [Mishchenko et al. 2008](https://doi.org/10.3233/JCM-2008-81-203) montrent que l'AI-REML standard échoue quand l'optimum est sur une borne. Ils proposent une méthode d'ensemble actif et une mise à jour BFGS inverse, pour des contraintes de non-négativité scalaires.

L'amortissement du type Levenberg-Marquardt est documenté. WOMBAT ajoute τI à la matrice AI, avec τ = max(0, δ − λmin), comme modification par défaut (--modaim2, annexe A.5). Il propose aussi le plaquage des valeurs propres et deux Cholesky modifiées. SAS PROC MIXED minimise -2logL par un Newton-Raphson stabilisé par crête ([Wolfinger et al. 1994](https://doi.org/10.1137/0915079) ; documentation SAS). Le démarrage par quelques itérations PX-EM avant l'AI est le réglage par défaut de WOMBAT pour plus de 18 composantes (--pxai). ASReml propose des itérations EM ou PX-EM locales. Le PX-EM vient de [Liu et al. 1998](https://doi.org/10.1093/biomet/85.4.755) et a été adapté au modèle de Henderson par [Foulley et van Dyk 2000](https://doi.org/10.1186/1297-9686-32-2-143). WOMBAT garde Powell et le simplexe pour vérifier la convergence des cas difficiles.

Les autres composants sont classiques en optimisation numérique. L-BFGS-B vient de [Byrd et al. 1995](https://doi.org/10.1137/0916069). Le pas de région de confiance par décomposition propre suit [Moré et Sorensen 1983](https://doi.org/10.1137/0904038). La convergence de SR1 est établie par [Conn et al. 1991](https://doi.org/10.1007/BF01594934). [Nocedal et Wright 2006](https://doi.org/10.1007/978-0-387-40065-5) rassemblent ces méthodes. La projection sur le cône semi-défini positif suit [Higham 1988](https://doi.org/10.1016/0024-3795%2888%2990223-6). Aucun logiciel REML examiné ne combine un ensemble actif sur le cône semi-défini, des salves de L-BFGS-B et un finisseur Moré-Sorensen avec SR1.

Verdict. L'AI, l'amortissement, le plancher et le démarrage EM ont une antériorité claire. L'ensemble actif a une antériorité partielle, limitée aux bornes scalaires. Son extension au cône semi-défini en échelle des variances, et l'assemblage complet, n'ont pas d'antériorité trouvée.

### REML sur GPU et dérivation automatique

La dérivation automatique appliquée au REML date de [Meyer et Smith 1996](https://doi.org/10.1186/1297-9686-28-1-23), sur CPU. TMB généralise la dérivation automatique aux modèles à effets latents, avec CHOLMOD et le sous-ensemble inverse ([Kristensen et al. 2016](https://doi.org/10.18637/jss.v070.i05)). Le parallélisme décrit dans cet article passe par OpenMP et BLAS, sans GPU. Hors génétique, GPyTorch calcule la vraisemblance marginale d'un processus gaussien et ses dérivées sur GPU, par gradient conjugué par lots ([Gardner et al. 2018](https://doi.org/10.48550/arXiv.1809.11165)). GPJax fait de même en JAX ([Pinder et Dodd 2022](https://doi.org/10.21105/joss.04455)). GPyTorch et GPJax maximisent la vraisemblance complète et ne traitent pas le REML.

En génétique, des concurrents directs sont apparus en 2026. pyreml ajuste des modèles mixtes par différenciation directe du REML en PyTorch, avec accélération GPU ([Marchal et Raimondi 2026](https://doi.org/10.5281/zenodo.21394319)). La version 0.5.3 examinée utilise L-BFGS, avec Adam en repli (Optimizer.py). TorchLIMIX porte le pipeline multivarié LIMIX sur GPU, avec des gains jusqu'à deux ordres de grandeur ([Horn et al. 2026](https://doi.org/10.64898/2026.04.29.721342), prépublication). SAIGE-GPU accélère l'ajustement du GLMM de SAIGE d'un facteur 5 ([Rodriguez et al. 2026](https://doi.org/10.1093/bioinformatics/btag032)). Le paquet torch-openreml annonce un moteur AI-REML sur CPU et GPU. Sa page PyPI n'a pas pu être lue ; il reste non vérifié. Aucune référence trouvée ne combine équations du modèle mixte creuses, Cholesky supernodale et inversion sélective dans un cadre de dérivation automatique sur GPU. Les travaux récents de l'équipe Misztal restent sur CPU, avec REML Monte Carlo et gradient conjugué préconditionné ([Bermann et al. 2025](https://doi.org/10.1186/s12711-025-01006-9)).

Verdict. La dérivation automatique du REML a une antériorité claire. Le REML sur GPU en forme V dense a une antériorité partielle et contemporaine (pyreml). La forme C creuse sur GPU, avec parité CPU/GPU validée, n'a pas d'antériorité trouvée.

### Modèles IGE de voisinage

Le cadre IGE est établi. [Griffing 1967](https://doi.org/10.1071/BI9670127) pose la sélection au niveau des groupes. [Moore et al. 1997](https://doi.org/10.1111/j.1558-5646.1997.tb01458.x) et [Wolf et al. 1998](https://doi.org/10.1016/S0169-5347%2897%2901233-0) formalisent les effets génétiques indirects. [Muir 2005](https://doi.org/10.1534/genetics.104.035956) les intègre au BLUP pour les arbres et les animaux. [Bijma et al. 2007a](https://doi.org/10.1534/genetics.106.062711) et [Bijma et al. 2007b](https://doi.org/10.1534/genetics.106.062729) définissent la valeur génétique totale et son estimation. La dépendance à la taille du groupe est traitée par [Hadfield et Wilson 2007](https://doi.org/10.1534/genetics.107.075622). [Bijma 2010a](https://doi.org/10.1534/genetics.110.120485) introduit le facteur de dilution d. Les dispositifs optimaux sont étudiés par [Bijma 2010b](https://doi.org/10.1534/genetics.110.120493). [Arango et al. 2005](https://doi.org/10.2527/2005.8361241x) estiment les effets compétitifs chez le porc.

Le voisinage spatial des plantes est modélisé depuis longtemps. Les ajustements de voisinage en essais au champ remontent à [Draper et Guttman 1980](https://doi.org/10.2307/2986297), [Kempton 1982](https://doi.org/10.1017/S0021859600054381) et [Besag et Kempton 1986](https://doi.org/10.2307/2531047). [Stringer et al. 2011](https://doi.org/10.1007/s13253-010-0051-5) modélisent conjointement la variabilité spatiale et la compétition génétique entre parcelles. Chez les arbres, [Cappa et Cantet 2008](https://doi.org/10.1515/sg-2008-0008) pondèrent les voisins par des « intensités de compétition », inverses de la distance et du nombre de voisins. Ils estiment le modèle par Gibbs et par EM-REML. [Cappa et al. 2015](https://doi.org/10.1007/s11295-015-0917-3) y ajoutent l'hétérogénéité environnementale. [Costa e Silva et al. 2013](https://doi.org/10.1111/nph.12035) opposent les IGE de maladie et de compétition chez l'eucalyptus. Chez la même espèce, [Costa e Silva et al. 2017](https://doi.org/10.1038/hdy.2017.25) identifient les positions de voisins les plus influentes par une structure factorielle, puis construisent un effet combiné pondéré par la distance. La revue de [Costa e Silva et Kerr 2013](https://doi.org/10.1007/s11295-012-0521-8) couvre ces modèles. [Sato et al. 2021](https://doi.org/10.1038/s41437-020-00401-w) estiment la portée effective des effets de voisinage au pic de variance expliquée, sur une grille d'échelles spatiales. Une prépublication plus récente évoque une décroissance avec la distance sans l'ajuster ([Sato et Hamazaki 2026](https://doi.org/10.7554/eLife.111650.1)).

Pour les mélanges de cultures, la littérature porte surtout sur l'aptitude au mélange et les interactions génotypiques. [Knott et Mundt 1990](https://doi.org/10.1007/bf00210065) analysent l'aptitude au mélange chez le blé. [Forst et al. 2019](https://doi.org/10.1016/j.fcr.2019.107571) proposent un cadre pour les dispositifs de mélange incomplets. [Montazeaud et al. 2022](https://doi.org/10.1111/nph.17915) testent la richesse allélique entre génotypes chez le blé dur. [Wuest et Niklaus 2018](https://doi.org/10.1038/s41559-018-0708-y) et [Wuest et al. 2021](https://doi.org/10.1038/s41559-021-01497-x) relient la diversité génétique au rendement des mélanges. [Pélissier et al. 2023](https://doi.org/10.1371/journal.pbio.3002287) montrent que le génotype voisin modifie la sensibilité aux maladies du blé et du riz.

Les logiciels utilisés sont variés. WOMBAT ajuste des effets « sociaux » avec un facteur de dilution fixé (manuel 2017, section 4.10.6.8). Il écrit la log-vraisemblance pour chaque d et propose une approximation quadratique du profil (--quapp). ASReml permet des superpositions de matrices d'incidence par and() (manuel 4.2). Son usage dans les études IGE forestières citées n'a pas été vérifié sur les textes. Cappa et Cantet emploient leur propre échantillonneur de Gibbs. Sato et collègues emploient rNeighborGWAS et RAINBOWR.

Verdict. Le modèle IGE, la dilution et le noyau discret ou pondéré par la distance ont une antériorité claire. L'estimation conjointe d'une portée continue par REML, avec les composantes de variance, n'a pas d'antériorité trouvée. Les approches existantes balaient une grille (Sato) ou profilent la vraisemblance (WOMBAT pour d). La recherche sur les mélanges interspécifiques avec IGE intra et inter espèces n'a pas été exhaustive.

### Incertitudes

Plusieurs textes primaires n'ont pas été lus en entier. La description de la matrice bordée de Smith et Graser 1986 provient de notes d'atelier MTDFREML. Les résultats de Meyer 2005 sur l'ordonnancement ne sont connus que par leur citation dans le manuel WOMBAT. Le choix d'ordre entre AMD et dissection emboîtée selon le type de matrice, chez Masuda et al. 2014, n'est connu que par le résumé. Les algorithmes internes d'ASReml sont fermés ; une heuristique non documentée de partition ou de contrainte ne peut pas être exclue. La version de pyreml examinée est la 0.5.3 du dépôt local. Son comportement exact (règle adjointe, voie Woodbury) mérite une lecture ciblée avant toute comparaison publiée. La page de torch-openreml n'a pas pu être chargée.

### Bibliographie

Articles et chapitres (DOI vérifiés sur CrossRef ou doi.org)

- Amestoy PR, Davis TA, Duff IS (1996). An approximate minimum degree ordering algorithm. SIAM J Matrix Anal Appl 17:886–905. https://doi.org/10.1137/S0895479894278952
- Amestoy PR, Davis TA, Duff IS (2004). Algorithm 837: AMD, an approximate minimum degree ordering algorithm. ACM TOMS 30:381–388. https://doi.org/10.1145/1024074.1024081
- Arango J, Misztal I, Tsuruta S, et al. (2005). Estimation of variance components including competitive effects of Large White growing gilts. J Anim Sci 83:1241–1246. https://doi.org/10.2527/2005.8361241x
- Bates D, DebRoy S (2004). Linear mixed models and penalized least squares. J Multivar Anal 91:1–17. https://doi.org/10.1016/j.jmva.2004.04.013
- Bates D, Mächler M, Bolker B, Walker S (2015). Fitting linear mixed-effects models using lme4. J Stat Softw 67(1). https://doi.org/10.18637/jss.v067.i01
- Bates D, Alday PM, Kokandakar AH (2025). Mixed-model log-likelihood evaluation via a blocked Cholesky factorization. arXiv:2505.11674 (prépublication). https://doi.org/10.48550/arXiv.2505.11674
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
- Horn B, Nikoloski Z, Lippert C (2026). TorchLIMIX: GPU-accelerated multivariate genome-wide association studies. bioRxiv (prépublication). https://doi.org/10.64898/2026.04.29.721342
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
- Mishchenko K, Holmgren S, Rönnegård L (2007). Efficient implementation of the AI-REML iteration for variance component QTL analysis. arXiv:0709.0625 (rapport de recherche). https://doi.org/10.48550/arXiv.0709.0625
- Mishchenko K, Holmgren S, Rönnegård L (2008). Newton-type methods for REML estimation in genetic analysis of quantitative traits. J Comput Methods Sci Eng 8:53–67. https://doi.org/10.3233/JCM-2008-81-203
- Misztal I, Perez-Enciso M (1993). Sparse matrix inversion for restricted maximum likelihood estimation of variance components by expectation-maximization. J Dairy Sci 76:1479–1483. https://doi.org/10.3168/jds.S0022-0302%2893%2977478-0
- Montazeaud G, Flutre T, Ballini E, et al. (2022). From cultivar mixtures to allelic mixtures. New Phytol 233:2573–2584. https://doi.org/10.1111/nph.17915
- Moore AJ, Brodie ED, Wolf JB (1997). Interacting phenotypes and the evolutionary process I. Evolution 51:1352–1362. https://doi.org/10.1111/j.1558-5646.1997.tb01458.x
- Moré JJ, Sorensen DC (1983). Computing a trust region step. SIAM J Sci Stat Comput 4:553–572. https://doi.org/10.1137/0904038
- Muir WM (2005). Incorporation of competitive effects in forest tree or animal breeding programs. Genetics 170:1247–1259. https://doi.org/10.1534/genetics.104.035956
- Murray I (2016). Differentiation of the Cholesky decomposition. arXiv:1602.07527. https://doi.org/10.48550/arXiv.1602.07527
- Nocedal J, Wright SJ (2006). Numerical Optimization, 2e éd. Springer. https://doi.org/10.1007/978-0-387-40065-5
- Pélissier R, Ballini E, Temple C, et al. (2023). The genetic identity of neighboring plants in intraspecific mixtures modulates disease susceptibility of both wheat and rice. PLoS Biol 21:e3002287. https://doi.org/10.1371/journal.pbio.3002287
- Pinder T, Dodd D (2022). GPJax: a Gaussian process framework in JAX. J Open Source Softw 7:4455. https://doi.org/10.21105/joss.04455
- Pinheiro JC, Bates DM (1996). Unconstrained parametrizations for variance-covariance matrices. Stat Comput 6:289–296. https://doi.org/10.1007/BF00140873
- Rodriguez A, Kim Y, Nandi TN, et al. (2026). SAIGE-GPU: accelerating genome- and phenome-wide association studies using GPUs. Bioinformatics. https://doi.org/10.1093/bioinformatics/btag032
- Sato Y, Yamamoto E, Shimizu KK, Nagano AJ (2021). Neighbor GWAS. Heredity 126:597–614. https://doi.org/10.1038/s41437-020-00401-w
- Sato Y, Hamazaki K (2026). Dissecting oligogenic and polygenic indirect genetic effects through the lens of neighbor genotypic identity. eLife, prépublication évaluée. https://doi.org/10.7554/eLife.111650.1
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

Manuels, logiciels et sources sans DOI

- Gilmour AR, Gogel BJ, Cullis BR, Welham SJ, Thompson R (2022). ASReml User Guide, Functional Specification, Release 4.2. VSN International. https://asreml.kb.vsni.co.uk/wp-content/uploads/sites/3/ASReml-4.2-Functional-Specification.pdf (sections 5.8, 6.10, 7.7.4, 11.5 consultées)
- Meyer K (2017). WOMBAT, a program for mixed model analyses by restricted maximum likelihood. User notes, version du 31 mai 2017. https://forge-dga.jouy.inra.fr/attachments/download/8063/WombatManual-V20170531.pdf (sections 1, 4.10.6.8, 5.2–5.3, annexes A.1 et A.5 consultées)
- Boldman K, Kriese L, Van Tassell C, Kachman S, Dodenhoff J (2002). MTDFREML Workshop Notes. USDA-ARS. https://ars.usda.gov/ARSUserFiles/80420530/MTDFREML/MTDFREML_Workshop_Notes.pdf
- Marchal A, Raimondi D (2026). pyreml, general-purpose linear mixed model solver, version 0.5.3 (code source examiné localement). https://doi.org/10.5281/zenodo.21394319
- SuiteSparse AMD, en-tête amd.h (seuil des lignes denses). https://github.com/DrTimothyAldenDavis/SuiteSparse/blob/stable/AMD/Include/amd.h
- Takahashi K, Fagan J, Chin MS (1973). Formation of a sparse bus impedance matrix and its application to short circuit study. 8th PICA Conference, Minneapolis, pp. 63–69. Sans DOI, non consulté.
- Meyer K (2005). Ordering strategies to reduce computational requirements in variance component estimation. Proc Assoc Advmt Anim Breed Genet 16:282–285. Sans DOI, non consulté.
