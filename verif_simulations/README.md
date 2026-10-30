# Vérification de remlax par simulation

Cette suite vérifie, pour chaque typologie de modèle que remlax sait ajuster, deux choses. La logvraisemblance REML rendue doit être exacte. Les composantes de variance estimées doivent retrouver les valeurs simulées. Ce n'est pas un banc de vitesse.

## Lancer

Depuis la racine du dépôt :

```sh
Rscript verif_simulations/verif_simulations.R --coeurs=10
```

Le script relance lui-même R avec `OPENBLAS_NUM_THREADS=1`. Chaque processus JAX est limité à deux fils (`XLA_FLAGS` et `options(remlax.threads = 2)`). Sur 20 cœurs, la suite complète prend une dizaine de minutes.

Options :

- `--pkg=` chemin des sources de remlax, chargées par `pkgload::load_all` ; sans cette option, le script charge le paquet installé par `R/charger_remlax.R` ;
- `--cas=` expression régulière sur les noms de cas (par exemple `--cas='^mt_'`) ; `--liste` affiche le catalogue ;
- `--coeurs=` nombre de processus parallèles (10 par défaut) ;
- `--rep=` nombre de réplicats par cas, sur des graines successives (1 par défaut) ;
- `--moteurs=` expression régulière sur les voies (`lbfgsb`, `ai`, `mme`, `mme_creuse`, `mme_lbfgsb`) ;
- `--graine=` graine de base ; `--sortie=` dossier des résultats (par défaut `resultats/`).

## Ce que fait chaque cas

Le cas construit un dispositif et simule y selon la V vraie du modèle. Il ajuste ensuite le modèle par chaque voie applicable. Les voies sont `rx_fit(optimizer = "lbfgsb")` et `rx_fit(optimizer = "ai")`. Le moteur des équations du modèle mixte, `rx_fit(engine = "mme")`, s'ajoute quand le modèle est dans son périmètre, par AI-REML en inversion automatique (`mme`) ou creuse (`mme_creuse`). Pour une parenté donnée par sa précision et pour le champ ar1 × ar1, une voie de plus ajuste le même moteur par L-BFGS-B en inversion creuse (`mme_lbfgsb`). Elle remplace `rx_fit_sparse()`, le moteur RTMB retiré de remlax en octobre 2026. Le backend est toujours le CPU.

L'arbitre recompose V en R à partir des Sigma rendues par l'ajustement (`fit$sigmas`, `fit$sigmas_res`) et des paramètres de niveaux de `fit$varcomp`. Il utilise la géométrie vraie du dispositif : coordonnées, cellules du champ, incidences construites par le cas. Il évalue ensuite la REML par la forme en V, avec la constante 2π, qui est la convention de `fit$logLik` :

    -2 logL = (n - p) log 2π + log|V| + log|X'V⁻¹X| + y'Py

Pour le décalage spatial (`rx_lag`), l'arbitre applique le critère `reml_y` : REML de y* = (I - M) y, jacobien -2 log|det A| et terme de moyenne. La simulation passe par la même fonction de V évaluée aux paramètres vrais. Le code de l'arbitre est dans `R/outils.R`, les cas dans `R/cas.R`.

## Critères

- arbitre : |logLik - arbitre| / |arbitre| < 1e-8 ;
- accord des voies : la logLik de chaque voie est à moins de 1e-4 de la meilleure du cas ;
- convergence : `conv_decrement` ou `conv_grad_rel` pour toutes les voies ; `scipy_success` est rapporté dans le détail sans trancher ;
- récupération : |estimé - vrai| < 3 erreurs-types. L'erreur-type vient de `fit$varcomp` (méthode delta sur le Hessien). Une voie sans erreur-type reprend celle de L-BFGS-B sur le même modèle.

Un ajustement est PASSE s'il satisfait les quatre critères. REFUS_ATTENDU signale une voie hors périmètre qui refuse le modèle avec le message prévu (AI-REML sur un modèle à décalage).

## Lire les résultats

Le dossier de sortie contient quatre fichiers :

- `resume.txt` : par cas, une ligne par voie (statut, logLik, arbitre, écarts, max|z|, raisons d'échec), puis les composantes de L-BFGS-B ;
- `ajustements.csv` : une ligne par cas et par voie, avec n, nombre de paramètres, logLik, logLik_arbitre, ecart_rel, ecart_entre_moteurs, convergence et son détail, durée ;
- `composantes.csv` : une ligne par composante, avec vrai, estime, se, la source de l'erreur-type et z ;
- `cas.csv` : famille, description et durée de chaque cas.

Les noms de composantes suivent `fit$varcomp$nom`. Une entrée de Sigma s'écrit `terme[i,j]`, un paramètre de niveaux `terme!famille` (suffixé 1, 2 quand la famille en a plusieurs), une section résiduelle porte le nom de son niveau dans `dsum`.

Avec un seul réplicat et environ 150 composantes, une ou deux valeurs |z| ≥ 3 peuvent apparaître par hasard. Un `--rep=5` sur le cas concerné tranche.

Deux familles de composantes donnent des |z| élevés sans défaut de calcul. Les variances des termes `rx_spl2d` reposent sur 7 à 9 niveaux : leur distribution est très asymétrique et l'erreur-type asymptotique ne la décrit pas. Sur 24 réplicats, la moyenne de chaque variance reste à moins de deux erreurs-types de la moyenne de la vérité. La portée de `lvr` entre dans une corrélation à pointe : la logvraisemblance a un point anguleux chaque fois que la portée croise une distance entre deux positions, et le Hessien par différences finies donne une erreur-type sans signification (0,01 pour un écart réel de 0,1 à 1,4).

Le dossier `resultats_replicats/` contient quatre réplicats des cas `mt_g_fa2`, `spline_spl2d`, `ige_melange_deux_especes` et `metrique_lvr`.
