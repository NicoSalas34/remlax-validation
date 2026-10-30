# Plan du banc : moteurs remlax, CPU et GPU, contre asreml, sommer et lme4

Rédigé le 2026-10-22, révisé le même jour. Ce banc doit fonder la règle `engine = "auto"` de `rx_fit`, la règle d'inversion du moteur MME et le choix du périphérique. Il compare aussi remlax à asreml, sommer et lme4, sur les temps, sur la logLik atteinte et sur les paramètres estimés, dont la vérité est connue par simulation.

## Questions

- Quand le moteur MME est-il plus rapide que le moteur dense, et quand le dense est-il nécessaire ? Le dense est nécessaire hors du périmètre MME (résiduelle none, rr, rx_lag).
- Dans le moteur MME, à partir de quel remplissage du facteur l'inversion creuse (SuperLU et inversion sélective supernodale) devient-elle plus lente que l'inversion dense de C ?
- Pour chaque moteur, à partir de quelle taille le GPU l'emporte-t-il sur le CPU, compilation et transferts compris ?
- remlax retrouve-t-il le même optimum qu'asreml, sommer et lme4, et en combien de temps ?

## Voies comparées

| Voie | Moteur | Inversion | Périphérique |
|---|---|---|---|
| D-CPU | dense (V formée) | Cholesky de V | CPU |
| D-GPU | dense | Cholesky de V | GPU |
| M-dense-CPU | MME | C dense (LAPACK) | CPU |
| M-dense-GPU | MME | C dense (JAX) | GPU |
| M-creux-CPU | MME | C creuse (SuperLU et inversion sélective) | CPU |
| asreml | asreml-R 4.2 | AI-REML creux | CPU |
| sommer | sommer 4.4 (`mmer`) | formulation directe en V | CPU |
| lme4 | lme4 2.0 (`lmer`, pedigreemm pour le pédigrée) | Cholesky creux | CPU |

L'optimiseur de remlax est l'AI en échelle des variances (`optimizer = "ai"`) pour toutes les voies, avec les mêmes bornes, le même point de départ et le même critère d'arrêt. Le hessien de diagnostic est désactivé (`hessian = FALSE`).

## Familles de modèles

Chaque famille couvre un régime de rapport entre n (observations) et m (taille de C) et de remplissage de C. Un logiciel externe n'est lancé que s'il sait ajuster exactement le même modèle.

| Famille | Modèle | Externes |
|---|---|---|
| iid | génotype iid, bloc iid | asreml, sommer, lme4 |
| grm | génotype avec GRM dense, bloc iid | asreml, sommer, lme4 (Z L, K = L L') |
| us | 3 caractères, us génétique, résiduelle us par unité | asreml, sommer |
| ped | modèle animal, A^-1 creuse, bloc iid | asreml, sommer, lme4 (pedigreemm) |
| ige | DGE et IGE us(2) sur incidences de voisinage, IEE iid | asreml (n ≤ 3000) |
| champ | génotype iid, champ ar1 x ar1 en terme aléatoire | asreml |
| resid | génotype iid, résiduelle ar1(row):ar1(col) | asreml |
| exp | corrélation exponentielle sur coordonnées 1D | asreml |
| fa | 5 caractères, fa de rang 2, résiduelle diag | asreml |

sommer et lme4 n'estiment ni une corrélation ar1 ou exponentielle, ni une structure fa ; ils sont absents de ces familles.

## Tailles et répétitions

n vaut 500, 1000, 1500, 2000, 3000, 4000, 6000 et 8000. n compte les observations, tous caractères confondus. n croît en ajoutant des génotypes, et le nombre d'observations par génotype reste à peu près fixe : m, la taille de C, croît donc avec n. Un effet de n ne doit pas se lire comme un effet du nombre de répétitions par génotype.

| Famille | Niveaux du terme principal | Observations par niveau et par caractère | m à n = 500 | m à n = 8000 |
|---|---|---|---|---|
| iid | génotypes = n/5 | 5 | 126 | 1615 |
| grm | génotypes = n/4 | 4 | 152 | 2027 |
| us | génotypes = n/12 (3 caractères) | 4 | 126 | 2001 |
| ped | animaux = n | 1 | 522 | 8022 |
| ige | génotypes = n/5 | 5 | 701 | 11201 |
| champ | génotypes = n/4 | 4 | 624 | 9967 |
| resid | génotypes = n/4 | 4 | 125 | 1963 |
| exp | positions = n/3 | 3 | 167 | 2667 |
| fa | génotypes = n/20 (5 caractères) | 4 (2,5 à n = 500, plancher de 40 génotypes) | 205 | 2005 |

Les génotypes sont tirés avec remise : le dispositif n'est pas équilibré et le nombre d'observations par génotype fluctue autour de sa moyenne. `banc_moteurs/R/dimensions.R` reconstruit chaque cas à partir de sa graine et écrit `dimensions.csv`. L'analyse en tire une figure des dimensions contre n, la règle d'échelle de chaque famille dans le titre de ses panneaux, et une figure des temps contre m.

Chaque cellule (famille, n) est simulée cinq fois. La graine dépend de la famille, de n et de la répétition, et toutes les voies et tous les logiciels d'une répétition lisent les mêmes données. Les 360 cas (9 familles, 8 tailles, 5 répétitions) forment un tableau SLURM.

## Machine

Tout tourne sur le nœud io-gpu-a100-01, tranche MIG 1g.10gb (1/7 de l'A100), 6 cœurs et 32 Go par tâche, conteneur `ige_reml.sif`. Les voies CPU et les logiciels externes utilisent les cœurs de la tranche, avec OMP, OpenBLAS et MKL réglés sur 6 fils. CPU et GPU partagent donc la même machine. Chaque ajustement tourne dans un processus R à part (callr), avec un délai maximal de 1800 s. Le temps est mesuré dans ce processus, autour du seul appel d'ajustement.

## Mesures enregistrées

Une ligne par ajustement : famille, n, répétition, graine, logiciel, voie, statut (OK, ERREUR, DELAI), temps total, `compile_s`, `eval_s`, nombre d'itérations et d'évaluations, moteur et inversion effectifs, convergence déclarée, logLik native et logLik arbitre. La logLik arbitre est la REML dense commune, évaluée aux estimations de chaque logiciel : les constantes des logLik natives diffèrent d'un logiciel à l'autre. Une seconde table donne, pour chaque ajustement, la valeur vraie et la valeur estimée de chaque paramètre.

## Analyse

Dans chaque répétition, la référence est la meilleure logLik arbitre parmi les ajustements réussis. Un ajustement est juste à moins de 1e-3 de cette référence. Il est retenu pour les temps et les paramètres à moins de 0,1. Les rapports de temps se calculent dans la répétition, puis se résument par la médiane, le minimum et le maximum des cinq répétitions.

Figures produites par `banc_moteurs/R/analyser.R` :

- dimensions des cas (génotypes, observations par génotype, m) contre n ;
- état des ajustements par famille et par voie ;
- temps contre n, toutes voies, en temps total puis en temps du solveur remlax seul ;
- temps contre m ;
- gain du GPU (CPU / GPU) pour le dense et le MME ;
- MME contre dense, et inversion creuse contre dense ;
- remlax contre asreml ;
- voie la plus rapide par cellule, toutes voies puis remlax seul ;
- écart à la meilleure logLik ;
- récupération des paramètres (biais, RMSE relative) ;
- accord des estimations des logiciels externes avec remlax ;
- part de la compilation JAX dans le temps total.

La règle `auto` se construit sur la carte de la voie remlax la plus rapide, à partir de prédicteurs calculables avant l'ajustement (n, m, nnz du motif de C, remplissage symbolique, présence d'une résiduelle structurée ou d'une K dense paramétrée).

## Lancement

```
W=/scratch/users/salasn/banc_moteurs_20261022
sbatch --array=0-359%7 $W/remlax-validation/banc_moteurs/banc_mig.sbatch
Rscript banc_moteurs/R/analyser.R --res=<copie locale de $W/res>
```

Une tâche dont le CSV existe déjà s'arrête aussitôt : un tableau relancé ne refait que les cas manquants.

## Limites connues

- Une tranche MIG sous-estime le GPU pour les grandes tailles. L'A100 entière n'est pas utilisée.
- asreml `expv` en terme aléatoire ne retrouve pas l'optimum REML du modèle simulé (ρ écrasé vers 0, logLik arbitre 20 unités sous celle de remlax à n = 500, même sur un cas à une observation par position). L'écart est signalé, pas corrigé.
- asreml n'ajuste le modèle ige qu'à n ≤ 3000 : les incidences de voisinage passent par des colonnes denses du tableau de données.
