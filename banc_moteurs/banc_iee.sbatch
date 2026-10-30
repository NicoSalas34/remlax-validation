#!/bin/bash
# ==============================================================================
# Bancs IEE : famille mtigee<t> (DGE/IGE us(2t) + IEE au niveau de la plante
# us(t) + residuelle us(t) par plante ; R/familles_params.R). m est de l'ordre
# de n et C presque pleine, comme dans mvD7.
#   taches 0-17  : nombre de parametres, t = 8,6,4,3,2,1 a n = 4000, 3 repetitions
#   taches 18-57 : taille, mtigee3 a n = 8000..500, 5 repetitions
# Toutes les voies remlax (dense AI et L-BFGS-B, CPU et GPU ; MME dense CPU/GPU ;
# MME creux) et asreml (n <= 4000), tranche MIG 1g.10gb.
#
#   sbatch --array=0-57%6 $W/remlax-validation/banc_moteurs/banc_iee.sbatch
# ==============================================================================
#SBATCH --job-name=banc_iee
#SBATCH --partition=gpu-cirad-dedicated
#SBATCH --account=dedicated-gpu@cirad-long
#SBATCH --qos=dedicated
#SBATCH --gres=gpu:nvidia_a100_80gb_pcie_1g.10gb:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --time=16:00:00
#SBATCH --output=/scratch/users/salasn/banc_moteurs_20261022/logs_iee/%A_%a.out
#SBATCH --error=/scratch/users/salasn/banc_moteurs_20261022/logs_iee/%A_%a.out
set -uo pipefail
W=${W:-/scratch/users/salasn/banc_moteurs_20261022}
SIF=${SIF:-/scratch/users/salasn/ige_reml.sif}
i=${SLURM_ARRAY_TASK_ID:?tableau SLURM requis}
if [ $i -lt 18 ]; then
  FAMS=(mtigee8 mtigee6 mtigee4 mtigee3 mtigee2 mtigee1)
  r=$(( i % 3 + 1 )); f=${FAMS[$(( i / 3 ))]}; n=4000; sortie="$W/res_iee_params"
else
  j=$(( i - 18 )); TAILLES=(8000 6000 4000 3000 2000 1500 1000 500)
  r=$(( j % 5 + 1 )); n=${TAILLES[$(( j / 5 ))]}; f=mtigee3; sortie="$W/res_iee_taille"
fi
EXT=asreml; [ $n -gt 4000 ] && EXT=none
cle="${f}_${n}_${r}"
fcsv="$sortie/ajustements_${cle}.csv"
if [ -f "$fcsv" ] && ! grep -q "License error\|licenses in use\|workspace" "$fcsv"; then echo "deja fait : $cle"; exit 0; fi
NT=${SLURM_CPUS_PER_TASK:-6}
mkdir -p "$sortie" "$W/tmp" "$W/cache"
echo "=== $(hostname) | ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} | $cle | $(date) ==="
apptainer exec --nv \
  --env RENV_CONFIG_AUTOLOADER_ENABLED=FALSE --env R_PROFILE_USER=/dev/null \
  --env R_LIBS="$W/rlib" --env TMPDIR="$W/tmp" --env REMLAX_CACHE="$W/cache" \
  --env OMP_NUM_THREADS=$NT --env OPENBLAS_NUM_THREADS=$NT --env MKL_NUM_THREADS=$NT \
  --env XLA_PYTHON_CLIENT_PREALLOCATE=false --env JAX_ENABLE_X64=1 \
  "$SIF" Rscript "$W/remlax-validation/banc_moteurs/R/executer.R" \
    --famille=$f --n=$n --rep=$r --sortie="$sortie" --delai=${DELAI:-1800} \
    --voies=D-CPU,D-GPU,D-CPU-LBFGS,D-GPU-LBFGS,M-dense-CPU,M-dense-GPU,M-creux-CPU --externes=$EXT
echo "=== fin $(date) | code $? ==="
