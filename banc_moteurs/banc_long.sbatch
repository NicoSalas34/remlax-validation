#!/bin/bash
# ==============================================================================
# Bancs multi-caracteres en temps longs (2026-10-23, dispositif corrige).
#
# NOMBRE DE PARAMETRES a nombre d'UNITES fixe : n = 1000 t, donc 1000 plantes et
# 200 genotypes (250 pour mtus) quel que soit t. A n fixe, le nombre de
# genotypes baissait avec t (100 a t = 8 pour 136 parametres dans G) : les
# petites valeurs propres de G (0,05) n'etaient plus estimables, l'optimum REML
# tombait au bord (rang 12 sur 16) et le banc mesurait la mauvaise
# identification autant que le nombre de parametres.
#   cas 0-35  : mtigee, mtige, mtus x t = 8,6,4,3,2,1 x 2 repetitions, n = 1000 t
#   cas 36-50 : taille, mtigee3 a n = 8000,4000,2000,1000,500, 3 repetitions
# UNE TACHE PAR (CAS, VOIE), delai de 48 h par ajustement. Voies denses CPU
# retirees pour t >= 6 (plus de 9 h par ajustement sur une tranche MIG), et voies
# MME CPU retirees pour mtigee6 et mtigee8 (plus de 8 a 12 h).
# Sorties : res_long/<voie>/.
#
#   sbatch --array=0-407%7 $W/remlax-validation/banc_moteurs/banc_long.sbatch
# ==============================================================================
#SBATCH --job-name=banc_long
#SBATCH --partition=gpu-cirad-dedicated
#SBATCH --account=dedicated-gpu@cirad-long
#SBATCH --qos=dedicated
#SBATCH --gres=gpu:nvidia_a100_80gb_pcie_1g.10gb:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --time=3-00:00:00
#SBATCH --output=/scratch/users/salasn/banc_moteurs_20261022/logs_long/%A_%a.out
#SBATCH --error=/scratch/users/salasn/banc_moteurs_20261022/logs_long/%A_%a.out
set -uo pipefail
W=${W:-/scratch/users/salasn/banc_moteurs_20261022}
SIF=${SIF:-/scratch/users/salasn/ige_reml.sif}
VOIES=(D-CPU D-GPU D-CPU-LBFGS D-GPU-LBFGS M-dense-CPU M-dense-GPU M-creux-CPU asreml)
i=${SLURM_ARRAY_TASK_ID:?tableau SLURM requis}
c=$(( i / ${#VOIES[@]} )); v=${VOIES[$(( i % ${#VOIES[@]} ))]}
if [ $c -lt 36 ]; then
  TS=(8 6 4 3 2 1); FS=(mtigee mtige mtus)
  r=$(( c % 2 + 1 )); k=$(( c / 2 )); f=${FS[$(( k % 3 ))]}; t=${TS[$(( k / 3 ))]}
  f="$f$t"; n=$(( 1000 * t ))
else
  j=$(( c - 36 )); TAILLES=(8000 4000 2000 1000 500)
  r=$(( j % 3 + 1 )); n=${TAILLES[$(( j / 3 ))]}; f=mtigee3; t=3
fi
case "$v" in
  D-CPU|D-CPU-LBFGS) [ $t -ge 6 ] && { echo "voie dense CPU retiree pour $f (t >= 6)"; exit 0; } ;;
  M-dense-CPU|M-creux-CPU)
    # mtigee (IEE au niveau de la plante) : C presque pleine, m ~ n ; plus de
    # 8 a 12 h par ajustement sur CPU pour t >= 6 (2026-10-24), retire
    case "$f" in mtigee6|mtigee8) echo "voie MME CPU retiree pour $f"; exit 0 ;; esac ;;
esac
if [ "$v" = "asreml" ]; then ARGS="--voies=none --externes=asreml"; else ARGS="--voies=$v --externes=none"; fi
sortie="$W/res_long/$v"
cle="${f}_${n}_${r}"
fcsv="$sortie/ajustements_${cle}.csv"
if [ -f "$fcsv" ] && ! grep -q "License error\|licenses in use\|workspace" "$fcsv"; then echo "deja fait : $cle $v"; exit 0; fi
NT=${SLURM_CPUS_PER_TASK:-6}
mkdir -p "$sortie" "$W/tmp" "$W/cache"
echo "=== $(hostname) | ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} | $cle | $v | $(date) ==="
apptainer exec --nv \
  --env RENV_CONFIG_AUTOLOADER_ENABLED=FALSE --env R_PROFILE_USER=/dev/null \
  --env R_LIBS="$W/rlib_mme2" --env TMPDIR="$W/tmp" --env REMLAX_CACHE="$W/cache" \
  --env OMP_NUM_THREADS=$NT --env OPENBLAS_NUM_THREADS=$NT --env MKL_NUM_THREADS=$NT \
  --env XLA_PYTHON_CLIENT_PREALLOCATE=false --env JAX_ENABLE_X64=1 \
  "$SIF" Rscript "$W/remlax-validation/banc_moteurs/R/executer.R" \
    --famille=$f --n=$n --rep=$r --sortie="$sortie" --delai=${DELAI:-172800} $ARGS
echo "=== fin $(date) | code $? ==="
