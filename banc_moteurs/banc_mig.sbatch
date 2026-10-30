#!/bin/bash
# ==============================================================================
# Banc des moteurs remlax contre asreml, sommer et lme4, sur tranche MIG 1g.10gb.
#
# Une tache par cas (famille, n, repetition) : 9 familles x 8 tailles x 5
# repetitions = 360 taches. Chaque tache ajuste les cinq voies remlax (dense et
# MME, CPU et GPU) puis les trois logiciels externes, sur les MEMES donnees
# (graine fonction de famille, n et repetition). Toutes les voies CPU tournent
# sur les coeurs alloues a la tranche : CPU et GPU partagent donc la machine.
#
#   W=/scratch/users/salasn/banc_moteurs_20261022
#   sbatch --array=0-359%7 $W/remlax-validation/banc_moteurs/banc_mig.sbatch
#
# Les plus grandes tailles passent en premier (indices bas), pour que la queue
# du tableau soit faite de cas courts.
# ==============================================================================
#SBATCH --job-name=banc_mot
#SBATCH --partition=gpu-cirad-dedicated
#SBATCH --account=dedicated-gpu@cirad-long
#SBATCH --qos=dedicated
#SBATCH --gres=gpu:nvidia_a100_80gb_pcie_1g.10gb:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/users/salasn/banc_moteurs_20261022/logs/%A_%a.out
#SBATCH --error=/scratch/users/salasn/banc_moteurs_20261022/logs/%A_%a.out
set -uo pipefail

W=${W:-/scratch/users/salasn/banc_moteurs_20261022}
SIF=${SIF:-/scratch/users/salasn/ige_reml.sif}
FAMILLES=(iid grm us ped ige champ resid exp fa)
TAILLES=(8000 6000 4000 3000 2000 1500 1000 500)
NREP=5

i=${SLURM_ARRAY_TASK_ID:?tableau SLURM requis}
nf=${#FAMILLES[@]}
r=$(( i % NREP + 1 )); i=$(( i / NREP ))
f=${FAMILLES[$(( i % nf ))]}; i=$(( i / nf ))
n=${TAILLES[$i]}
cle="${f}_${n}_${r}"
# cas complet et sans echec de licence : rien a faire ; sinon executer.R reprend
# ce qui manque (fichiers .partiel, temoin en_cours_<cas>)
fcsv="$W/res/ajustements_${cle}.csv"
if [ -f "$fcsv" ] && ! grep -q "License error\|licenses in use\|workspace" "$fcsv"; then echo "deja fait : $cle"; exit 0; fi

NT=${SLURM_CPUS_PER_TASK:-6}
mkdir -p "$W/res" "$W/tmp" "$W/cache"
echo "=== $(hostname) | ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} | $cle | $NT coeurs | $(date) ==="
nvidia-smi -L

# renv desarme (le .Rprofile du projet repointerait .libPaths vers l'hote) ;
# un seul nombre de fils pour tous les logiciels ; allocateur JAX sans
# preallocation (la tranche ne fait que 10 Go).
apptainer exec --nv \
  --env RENV_CONFIG_AUTOLOADER_ENABLED=FALSE --env R_PROFILE_USER=/dev/null \
  --env R_LIBS="$W/rlib" --env TMPDIR="$W/tmp" --env REMLAX_CACHE="$W/cache" \
  --env OMP_NUM_THREADS=$NT --env OPENBLAS_NUM_THREADS=$NT --env MKL_NUM_THREADS=$NT \
  --env XLA_PYTHON_CLIENT_PREALLOCATE=false --env JAX_ENABLE_X64=1 \
  "$SIF" Rscript "$W/remlax-validation/banc_moteurs/R/executer.R" \
    --famille=$f --n=$n --rep=$r --sortie="$W/res" --delai=${DELAI:-1800}
rc=$?
echo "=== fin $(date) | code $rc ==="
exit $rc
