# Nettoyage de remlax : diagnostic et plan

Note du 2026-10-29. Elle porte sur l'arbre de travail de remlax : commit
777435a, plus le retrait de RTMB, le README retravaillé et les pages d'aide
régénérées par roxygen2 8.1.1. Ces trois changements ne sont pas encore
commités.

Le diagnostic vient d'une lecture complète du code, en quatre parties : le
côté R, la couche d'orchestration Python (`fit.py`, `cli.py`, `inference.py`,
`serveur.py`), le modèle et les moteurs denses (`model.py`, `levels.py`,
`econome.py`, `blocs.py`, `lag.py`), puis le moteur MME, l'AI-REML, la voie
bayésienne et le scan.

## Objectif et contraintes

Le code doit pouvoir être lu, compris et repris par quelqu'un d'autre. Chaque
fonction d'entrée doit se lire comme une suite d'étapes nommées, et chaque
étape doit être une fonction courte.

Le nettoyage ne change aucun résultat. Il garde intacts :

- les noms et arguments des fonctions R exportées ;
- les options de la ligne de commande du solveur ;
- le format du paquet de données échangé entre R et Python ;
- les clés de `result.json` et les fichiers `out_*.bin` ;
- l'ordre des évaluations de l'objectif. Le mémo de `fit.py` ne garde qu'une
  valeur par clé, et l'AI-REML en dépend : déplacer un appel change le
  nombre d'évaluations et parfois le dernier bit.

## Taille et forme du code

Le paquet compte 6 800 lignes de R et 12 300 lignes de Python. Quatre
fichiers concentrent l'essentiel de la difficulté :

| fichier | lignes | problème principal |
|---|---|---|
| `R/remlax.R` | 2 368 | tout le côté R dans un seul fichier |
| `fit.py` | 1 890 | `_fit_reml_coeur` fait 830 lignes et seize étapes |
| `mme.py` | 1 301 | `_evaluer` et `_derivees` font chacune 190 à 200 lignes |
| `ai.py` | 997 | `ai_reml_var` fait 390 lignes, avec des fermetures qui modifient l'état |

Plusieurs défauts reviennent dans tout le code :

- Des commentaires de journal de bord. Ils rapportent des mesures, des dates
  ou des bugs passés au lieu d'expliquer le code. Ils renvoient parfois à des
  documents absents du dépôt (`theorie_v2.md`, `note_remlax.md`, `mvD7`).
- Le français et l'anglais mêlés dans les noms internes.
- Des conventions recalculées à plusieurs endroits :
  - la nature des niveaux d'un terme, `tm.get("lvl") or ("fixed" if LK else "id")`, réécrite une vingtaine de fois ;
  - le découpage de theta, réécrit quatre fois en R et cinq fois en Python ;
  - l'écriture binaire du paquet, faite deux fois en R ;
  - quatre grammaires pour désigner une colonne de terme (`"terme:a"`).
- Un bloc d'import de secours (`try` relatif, `except` absolu) répété 57 fois,
  pour que les fichiers tournent aussi comme scripts.

## Code mort

Ces éléments n'ont aucun appelant. Les recherches ont été faites sur R/,
tests/, tests-python/ et inst/.

- R : `rx_here`, `rx_pkg_dir`, et quelques variables locales jamais lues.
- `mme_cli.py` : le fichier entier. Il importe des symboles qui n'existent
  plus et planterait au chargement.
- `ai.ai_reml` : remplacée par `ai_reml_var`.
- `decoupage.ordre_niveaux_mme` et `_NiveauxCreux.appliquer_K`.
- `fit.py`, lignes 641 et 644 : un calcul jamais utilisé, masqué par un
  `noqa`.
- `fit.py` : un `fill_diagonal` sans effet et un terme `* 0`.
- Des imports inutilisés et des attributs `compile_args` jamais lus.
- `core.py` : la passe arrière par morceaux, qui ne s'active qu'avec une
  variable d'environnement que rien ne pose.

Trois éléments demandent une décision plutôt qu'une suppression directe :
- `optimizer = "ai_theta"`, ancienne variante de l'AI gardée pour
  comparaison ;
- `sparsity.py`, diagnostic hors du chemin d'ajustement ;
- la structure Sigma `"fixed"`, qui n'est accessible que par un test.

## Défauts trouvés pendant la lecture

Le nettoyage ne les corrige pas, puisqu'il garantit des résultats identiques.
Chacun fera l'objet d'une correction séparée, avec son propre test, après le
nettoyage.

Côté R :

- Les analyseurs de formules évaluent `rank =`, `struct =` et `order =` dans
  leur propre environnement. `fa(trait, rank = rk)` lit alors la variable
  locale `rk` de l'analyseur, et une variable de l'utilisateur n'est pas
  trouvée.
- `rx_reml` avec `trait` donné par un nom de colonne et des lignes
  incomplètes transforme le nom en vecteur.
- `model.matrix` dépend de l'option globale `na.action`.
- `rx_scan` nomme son dossier à la seconde près : deux scans dans la même
  seconde se marchent dessus.

Côté Python :

- Avec des poids, Wald, Kenward-Roger et les prédictions reconstruisent V sans
  les poids. C'est probablement faux, à vérifier par un test.
- `grad_rel` a deux définitions. La bascule finale du MME vers le moteur dense
  dépend donc de `hessian = TRUE` ou `FALSE`.
- Le plancher relevé par le MME reste en place si l'ajustement bascule vers le
  moteur dense, mais le tableau des composantes utilise le plancher d'origine.
- Le contexte `jax.default_device` n'est jamais refermé sur les chemins
  d'erreur. Le serveur persistant le garde ouvert.
- Un terme donné seulement par `Kinv` est traité comme K = I par les moteurs
  denses, sans message.
- Le compteur `nit` de l'AI compte les points essayés, pas les itérations.
- Le noyau C de l'inversion sélective bascule sans message vers la boucle
  NumPy, 5 à 15 fois plus lente, quand la compilation échoue.
- La structure `sep` sur une section de résiduelle ne peut pas fonctionner.
  Son chemin n'est pas atteignable aujourd'hui.

## Filet de sécurité

Le filet a trois niveaux.

1. Résultats de référence de bout en bout
   (`remlax-validation/refacto/golden.R`). Le script ajuste 36 cas qui
   parcourent les chemins du solveur :
   - les moteurs dense, économe, MME dense et creux, et blocs ;
   - les optimiseurs L-BFGS-B, AI et automatique, avec et sans redémarrages ;
   - les structures us, fa, diag, ar1 × ar1, spl2d, termes explicites,
     décalage spatial et résiduelle `"none"` ;
   - les sorties PEV, Wald, Kenward-Roger, vpredict et prédictions ;
   - le scan, LOCO et le réajustement exact ;
   - la voie bayésienne, l'évaluation à theta fixé, les paramètres fixés et
     les aides de voisinage.

   Chaque objet rendu est enregistré, avec son `print()` et son `summary()`.
   Deux captures du même code sont identiques au bit près, durées imprimées
   exceptées. La référence est rangée dans
   `refacto/reference_avant_nettoyage/`. Après chaque étape, la comparaison
   doit répondre « aucun écart ».
2. Les suites existantes : 266 tests pytest et environ 890 vérifications
   testthat.
3. Pour les étapes qui touchent au MME et à l'AI, des valeurs de référence
   de (f, g, AI) à theta fixé, pour chaque voie d'inversion (dense CPU, blocs,
   CHOLMOD, SuperLU, hybride) et chaque forme de terme. On y ajoute la
   trajectoire de l'AI pas à pas, parce que l'optimum seul ne détecte pas un
   changement de chemin. Ce niveau reste à écrire.

Toutes les comparaisons se font sur CPU, avec un seul fil de calcul. Le
nombre de fils BLAS déplace la vraisemblance au sixième chiffre, et le GPU
n'est pas déterministe. `verif_simulations` sera relancé à la fin.

## Architecture cible

### Côté R

`R/remlax.R` est découpé par étape de la chaîne. Chaque fichier a une
fonction d'entrée qui enchaîne des étapes nommées.

| fichier | contenu |
|---|---|
| `counts.R` | comptes de paramètres et une source unique pour le découpage de theta |
| `term.R`, `residual.R`, `model.R` | `rx_term` (incidences, structure, parenté, niveaux), `rx_residual`, `rx_model` |
| `formula_random.R`, `formula_residual.R`, `spl2d.R` | les analyseurs de formules, avec leurs aides communes (cellules ar1, métriques 2D, options) |
| `reml.R` | `rx_reml` : effets fixes, lignes incomplètes, spl2d, analyse des formules, modèle, ajustement |
| `bundle_write.R`, `bundle_read.R` | un seul écrivain binaire, partagé par `rx_fit`, `rx_scan` et `rx_predict` ; la lecture du résultat |
| `python.R`, `launch.R`, `server.R` | trouver et vérifier Python, lancer le solveur (serveur ou sous-processus), le serveur persistant |
| `fit.R` | `rx_fit`, avec une table unique de correspondance entre options R et options du solveur, partagée avec `rx_fit_mme` |
| `predict.R`, `scan.R`, `design.R`, `ratios.R`, `bayes.R`, `lag.R`, `spatial.R`, `summary.R`, `validate.R` | une fonction d'entrée et ses étapes |

Les fonctions internes prennent un nom de la forme `.rx_<verbe>_<objet>`.
`rx_level_of` et `rx_scan_solver_args` portent aujourd'hui le préfixe public
sans être exportées. Les tests qui appellent des fonctions internes sont
mis à jour en même temps.

### Côté Python

- `cli.py`. `main` (375 lignes) devient une suite d'étapes : lire les
  arguments, préparer le processus, lire le paquet, ajuster, calculer
  l'inférence, écrire. Une table unique des tableaux de sortie remplace les
  quinze écritures et la liste noire de clés. Le contexte du périphérique
  passe par un bloc `with`.
- `fit.py` devient un sous-paquet :
  - départ : valeurs initiales, balayage des portées, bornes des loadings ;
  - moteurs : `construire_moteur`, avec un contrat unique pour dense,
    économe, blocs et MME ;
  - optimisation : le compteur d'évaluations, L-BFGS-B, AI, sortie du
    plancher, contrôle des points selles ;
  - polissage, diagnostic et sorties ;
  - `fit_reml`, écrit comme la suite de ses seize étapes.
- `inference.py` est scindé entre les composantes (noms, tableau, vpredict)
  et les effets fixes (Wald, Kenward-Roger, prédictions), avec une seule
  inversion de V et une seule covariance de theta.
- Modèle : un module `terms.py` qui porte la nature des niveaux et une
  disposition de theta unique, utilisée par tous les moteurs. `levels.py`
  est scindé en familles stationnaires, métriques et personnalisées.
- `mme.py` devient un sous-paquet. Chaque voie d'inversion (dense CPU,
  dense GPU, blocs, creuse, hybride) devient une classe avec la même
  interface, ce qui supprime cinq chaînes de `if`. Le cache du dernier point,
  aujourd'hui implicite, devient explicite. `remlax.mme` reste importable.
- `ai.py` : `ai_reml_var` est découpée en fonctions nommées (projection sur
  le cône, ensemble actif, faces, pas de Newton, recherche du pas), avec un
  petit objet d'état à la place des fermetures qui modifient leurs
  variables englobantes.
- `selinv.py` est découpé par sections : ordre, CHOLMOD, SuperLU,
  supernœuds, noyau C.

### Commentaires

Chaque commentaire garde la raison d'être du code en une ou deux phrases.
Les mesures, les dates et l'histoire des bugs vont dans
`docs/notes_solveur.md`, un document pour les développeurs qui part dans le
dépôt mais pas dans l'archive du paquet.

## Ordre des travaux

Chaque étape est un commit qui passe le filet complet.

1. Commiter l'état actuel comme base : retrait de RTMB, README et pages
   d'aide. Ouvrir une branche `nettoyage`.
2. Supprimer le code mort, les imports inutilisés et les variables jamais
   lues.
3. Déplacer les commentaires de journal vers `docs/notes_solveur.md`.
4. Côté R : découpage en fichiers, puis fonctions d'étapes, puis
   déduplications (écrivain binaire, table des options, grammaire des
   colonnes, disposition de theta).
5. Python, orchestration : `cli.py`, puis `fit.py`, puis `inference.py`.
6. Python, modèle et moteurs denses.
7. Écrire le troisième niveau du filet, puis restructurer le MME, l'AI et
   `selinv`.
8. Adapter la documentation, en particulier la vignette `engine`, qui
   décrit les modules, et ajouter un guide d'architecture pour les
   développeurs.
9. Corriger les défauts listés plus haut, un commit et un test par défaut.
   Les résultats changent alors, et chaque changement est expliqué.

C'est un travail de plusieurs séances. Les étapes 2 et 3 sont rapides et sans
risque. Les étapes 5 et 7 demandent le plus de soin.

## Décisions à prendre

1. La langue des noms internes et des commentaires. L'API publique est en
   anglais, les commentaires sont en français et les noms internes mêlent
   les deux. Je propose des noms internes en anglais, cohérents avec l'API,
   et des commentaires en anglais pour un paquet public. Garder les
   commentaires en français est possible, mais le mélange actuel ne doit pas
   rester.
2. Le sort de `optimizer = "ai_theta"`, de `sparsity.py` et de la structure
   Sigma `"fixed"`.
3. Le bloc d'import de secours. Le supprimer demande que R lance toujours le
   solveur par `python -m remlax.cli`. Aujourd'hui, R peut le lancer par le
   chemin du fichier, ce qui charge deux fois certains modules.
4. L'ordre proposé : nettoyage à résultats identiques d'abord, corrections
   ensuite.
