#!/bin/bash
# ==============================================================================
# Bancs IEE en temps longs : famille mtigee<t> (DGE/IGE us(2t) + IEE au niveau
# de la plante us(t) + residuelle us(t)). UNE TACHE PAR (CAS, VOIE) : les voies
# d'un cas tournent en parallele au lieu de s'additionner. Delai de 48 h par
# ajustement.
#   cas 0-11  : nombre de parametres, t = 8,6,4,3,2,1 a n = 4000, 2 repetitions
#   cas 12-26 : taille, mtigee3 a n = 8000,4000,2000,1000,500, 3 repetitions
#   voies     : D-CPU D-GPU D-CPU-LBFGS D-GPU-LBFGS M-dense-CPU M-dense-GPU
#               M-creux-CPU asreml (asreml jusqu'a n = 4000)
# Sorties : res_iee_long/<voie>/ajustements_<cas>.csv (une ligne par fichier).
#
# Version 2 (2026-10-23) : voies MME sur rlib_mme2 (reprise par le moteur dense
# des que le gradient relatif depasse 1e-5, au lieu de 1e-3) ; voies denses CPU
# retirees pour t = 6 et 8 (plus de 9 h par ajustement sur une tranche MIG).
#
#   sbatch --array=0-215%7 $W/remlax-validation/banc_moteurs/banc_iee_long2.sbatch
# ==============================================================================
#SBATCH --job-name=iee_long2
#SBATCH --partition=gpu-cirad-dedicated
#SBATCH --account=dedicated-gpu@cirad-long
#SBATCH --qos=dedicated
#SBATCH --gres=gpu:nvidia_a100_80gb_pcie_1g.10gb:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --time=3-00:00:00
#SBATCH --output=/scratch/users/salasn/banc_moteurs_20261022/logs_iee_long/%A_%a.out
#SBATCH --error=/scratch/users/salasn/banc_moteurs_20261022/logs_iee_long/%A_%a.out
set -uo pipefail
W=${W:-/scratch/users/salasn/banc_moteurs_20261022}
SIF=${SIF:-/scratch/users/salasn/ige_reml.sif}
VOIES=(D-CPU D-GPU D-CPU-LBFGS D-GPU-LBFGS M-dense-CPU M-dense-GPU M-creux-CPU asreml)
i=${SLURM_ARRAY_TASK_ID:?tableau SLURM requis}
c=$(( i / ${#VOIES[@]} )); v=${VOIES[$(( i % ${#VOIES[@]} ))]}
if [ $c -lt 12 ]; then
  FAMS=(mtigee8 mtigee6 mtigee4 mtigee3 mtigee2 mtigee1)
  r=$(( c % 2 + 1 )); f=${FAMS[$(( c / 2 ))]}; n=4000
else
  j=$(( c - 12 )); TAILLES=(8000 4000 2000 1000 500)
  r=$(( j % 3 + 1 )); n=${TAILLES[$(( j / 3 ))]}; f=mtigee3
fi
case "$v:$f" in
  D-CPU:mtigee8|D-CPU:mtigee6|D-CPU-LBFGS:mtigee8|D-CPU-LBFGS:mtigee6)
    echo "voie dense CPU retiree pour $f (plus de 9 h sur une tranche MIG)"; exit 0 ;;
esac
RLIB="$W/rlib"; case "$v" in M-*) RLIB="$W/rlib_mme2" ;; esac
if [ "$v" = "asreml" ]; then
  [ $n -gt 4000 ] && { echo "asreml non lance au-dela de n = 4000 : $f $n"; exit 0; }
  ARGS="--voies=none --externes=asreml"
else
  ARGS="--voies=$v --externes=none"
fi
sortie="$W/res_iee_long/$v"
cle="${f}_${n}_${r}"
fcsv="$sortie/ajustements_${cle}.csv"
if [ -f "$fcsv" ] && ! grep -q "License error\|licenses in use\|workspace" "$fcsv"; then echo "deja fait : $cle $v"; exit 0; fi
NT=${SLURM_CPUS_PER_TASK:-6}
mkdir -p "$sortie" "$W/tmp" "$W/cache"
echo "=== $(hostname) | ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} | $cle | $v | $(date) ==="
apptainer exec --nv \
  --env RENV_CONFIG_AUTOLOADER_ENABLED=FALSE --env R_PROFILE_USER=/dev/null \
  --env R_LIBS="$RLIB" --env TMPDIR="$W/tmp" --env REMLAX_CACHE="$W/cache" \
  --env OMP_NUM_THREADS=$NT --env OPENBLAS_NUM_THREADS=$NT --env MKL_NUM_THREADS=$NT \
  --env XLA_PYTHON_CLIENT_PREALLOCATE=false --env JAX_ENABLE_X64=1 \
  "$SIF" Rscript "$W/remlax-validation/banc_moteurs/R/executer.R" \
    --famille=$f --n=$n --rep=$r --sortie="$sortie" --delai=${DELAI:-172800} $ARGS
echo "=== fin $(date) | code $? ==="
