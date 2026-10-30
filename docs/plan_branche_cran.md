# remlax vers CRAN : branche de publication et conditions du gel

Note du 2026-10-29, sur main au commit 777435a. Elle prolonge l'audit
`cran_policy_audit.md` du 2026-09-29 (commit 48a75d2), fait avant l'arrivée
du moteur MME, de l'AI-REML, du serveur Python, de la voie bayésienne et du
décalage spatial.

## Où en est le paquet

`R CMD check --as-cran` du 2026-10-29, R 4.5.3, avec le solveur Python (RX_PY
défini dans `~/.Renviron`) :

- aucune erreur ;
- un avertissement, `qpdf` absent de la machine, sans rapport avec le paquet ;
- une note d'horloge (« future file timestamps »), propre à la machine ;
- une note sur l'exemple de `rx_bayes`, 5,1 s, et 32 s avec `--run-donttest` ;
- tests : 890 réussis, 1 ignoré, 197 s de temps écoulé ;
- vignettes : reconstruites sans erreur, leur code n'est pas exécuté.

Le check sans Python, celui que verront les machines du CRAN, reste à refaire.
Il faut neutraliser RX_PY dans l'environnement du check lui-même, par exemple
avec `R_ENVIRON_USER=/dev/null`.

## Ce qui bloque depuis l'audit

### Écritures dans le répertoire personnel

Le solveur écrit par défaut dans `~/.cache/remlax` :

- le cache de compilation XLA (`device.py`, `activer_cache_compilation`) ;
- la bibliothèque compilée du noyau C de l'inversion sélective
  (`selinv.py`, `libselinv_<hash>.so`).

Le check de ce jour y a écrit. La politique du CRAN interdit toute écriture
hors du répertoire temporaire de la session, sauf dans
`tools::R_user_dir()`, à condition que le contenu reste petit et géré. Le
côté R ne transmet aujourd'hui aucun `REMLAX_CACHE` au solveur.

Correction proposée :

- R transmet `REMLAX_CACHE = tools::R_user_dir("remlax", "cache")` à chaque
  lancement du solveur ;
- pendant un check (`_R_CHECK_LIMIT_CORES_` défini, ou session non
  interactive sans accord), le cache va dans `tempdir()` ;
- une fonction `rx_cache_clear()` vide ce dossier et en donne la taille.

### Compilation de code C à l'exécution

`selinv.py` compile `selinv_c.c` au premier usage, avec le compilateur du
système, puis charge la bibliothèque par `ctypes`. Un relecteur du CRAN
refusera très probablement du code compilé hors de l'installation.

Deux voies :

- déplacer le noyau dans `src/`. R le compile à l'installation dans `libs/`,
  avec sa propre chaîne. R transmet ensuite au solveur le chemin de la
  bibliothèque, que Python charge par `ctypes` comme aujourd'hui. Les
  pointeurs BLAS et LAPACK restent ceux de SciPy, passés à l'appel ;
- désactiver le noyau dans la version CRAN. La boucle NumPy prend le relais,
  avec le même résultat, mais 5 à 15 fois plus lentement sur l'inversion
  sélective (0,3 s contre 3,9 s sur un pedigree de 20 031 individus).

Je recommande la première voie. Elle garde la vitesse, et la compilation
passe par le mécanisme normal d'un paquet R.

### Exemple trop long

L'exemple de `rx_bayes` dépasse 5 s. Il passe en `\donttest`.

## Points à vérifier sans être bloquants

- Serveur Python persistant. Il est désactivé pendant un check et arrêté par
  un finaliseur à la fin de la session. Il reste à vérifier sous Windows
  qu'aucun processus ne survit à R.
- Fils de calcul. `fils.py` restreint l'affinité du processus. Pendant un
  check, un seul fil est utilisé, comme l'exige la politique. Il faut le
  revérifier sur win-builder.
- Plateformes. Le paquet n'a été vérifié que sous Linux. win-builder (R-devel
  et R-release), mac-builder et R-hub restent à faire. L'audit avait corrigé
  `system2(env =)`, ignoré sous Windows, sans l'avoir testé.
- DOI de Salas et al. (2026). Le préfixe 10.64898 de bioRxiv se résout bien
  (vérifié le 2026-10-29).
- Mainteneur. L'adresse institutionnelle doit rester valide après la thèse.
  L'ORCID et le rôle `cph` d'un employeur éventuel restent à décider.

## Décisions à prendre avant de geler l'API

Une fois sur le CRAN, renommer une fonction ou un argument demande un cycle
de dépréciation. Les choix suivants se font donc avant le gel.

### Langue de l'interface

L'interface mélange l'anglais et le français :

- arguments en français : `marqueurs`, `carte`, `codage`, `construire`,
  `K_de`, `bloc` ;
- une fonction en français : `rx_scan_seuil()` ;
- champs de sortie en français : `moteur`, `secondes`, `composantes_degenerees`,
  `signe_non_identifie`, `second_ordre_omis` ;
- messages et affichages du solveur en français (« Ajustement REML »,
  « parametres »).

Le CRAN ne l'interdit pas. Un paquet public sera pourtant lu surtout en
anglais. Je propose de passer l'interface en anglais avant le gel, en gardant
les anciens noms comme alias dépréciés pendant une version. Les messages
internes du solveur peuvent rester en français dans un premier temps.

### Périmètre de la première version

Certaines fonctions sont récentes et ont peu servi hors de nos données :
la voie bayésienne, le décalage spatial, le moteur par blocs, le moteur MME
hybride, `rx_scan_loco()` et `rx_scan_exact()`. Deux options :

- les publier, marquées expérimentales dans leur aide, sans promesse de
  stabilité de leur interface ;
- les garder sur main, hors de la version CRAN, jusqu'à la suivante.

Je recommande la première option, avec la mention expérimentale. Retirer du
code d'une branche de publication crée une divergence coûteuse à entretenir.

## Stratégie de branches

La branche `cran-layout` a déjà été fusionnée dans main, qui l'a dépassée de
huit commits. Entretenir une branche de mise en forme parallèle à main a
coûté cette fusion. La conformité au CRAN doit vivre sur main, et
`.Rbuildignore` exclut déjà ce qui ne part pas dans l'archive (`manual/`,
`tests-python/`, l'archive de bancs).

Proposition :

- main reste la branche de développement, en version `0.1.0.9000` après la
  première soumission ;
- au gel, une branche `release/0.1` part d'un commit étiqueté. Elle ne reçoit
  que des corrections demandées par le CRAN ;
- l'étiquette `v0.1.0` marque l'archive soumise. Chaque correction faite sur
  `release/0.1` est reportée sur main ;
- les branches `cran-layout` et `structural-lag`, entièrement contenues dans
  main, peuvent être supprimées.

## Intégration continue

Le dépôt n'a pas de `.github/` aujourd'hui, bien que `.Rbuildignore` le
prévoie. Il faudrait deux familles de jobs GitHub Actions :

- un check `--as-cran` sans Python, sous Linux, Windows et macOS, avec
  R-release et R-devel (`r-lib/actions/check-r-package`). Il reproduit ce que
  verra le CRAN ;
- un job Linux avec Python et JAX sur CPU, qui lance les tests testthat avec
  le solveur et la suite pytest.

## Conditions du gel

Le paquet est prêt pour une branche de publication quand :

- les blocages ci-dessus sont corrigés sur main ;
- l'interface est figée, langue comprise ;
- les tests ajoutés le 2026-10-29 (`test-scan-plus.R`, `test_face_cone.py`)
  sont commités, et la suite complète passe avec et sans Python ;
- `verif_simulations` (remlax-validation) passe sans nouvel échec ;
- le check `--as-cran` donne 0 erreur, 0 avertissement et seulement la note
  « New submission » sur win-builder, mac-builder et l'intégration continue ;
- `NEWS.md` porte un résumé court de la première version, et le journal de
  développement actuel est déplacé ailleurs ;
- `cran-comments.md` est à jour.

## Ordre proposé

On corrige d'abord les blocages sur main : cache, noyau C, exemple trop long.
On tranche ensuite la langue et le périmètre, puis on fait les renommages. On
installe l'intégration continue. On gèle alors : étiquette, branche
`release/0.1`, checks externes. La soumission vient en dernier.
