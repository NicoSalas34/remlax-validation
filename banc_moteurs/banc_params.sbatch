#!/bin/bash
# ==============================================================================
# "Number of parameters" benchmark: multi-trait families mtus<t> and mtige<t>
# (R/familles_params.R), t = 1..8 traits, n = 4000 observations, every remlax
# route (dense CPU/GPU with AI and L-BFGS-B, MME dense CPU/GPU, MME sparse) and
# asreml, on a MIG slice 1g.10gb. 12 families x 3 replicates = 36.
#
#   sbatch --array=0-35%6 $W/remlax-validation/banc_moteurs/banc_params.sbatch
# ==============================================================================
#SBATCH --job-name=bench_par
#SBATCH --partition=gpu-cirad-dedicated
#SBATCH --account=dedicated-gpu@cirad-long
#SBATCH --qos=dedicated
#SBATCH --gres=gpu:nvidia_a100_80gb_pcie_1g.10gb:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/users/salasn/banc_moteurs_20261022/logs_params/%A_%a.out
#SBATCH --error=/scratch/users/salasn/banc_moteurs_20261022/logs_params/%A_%a.out
set -uo pipefail
W=${W:-/scratch/users/salasn/banc_moteurs_20261022}
SIF=${SIF:-/scratch/users/salasn/ige_reml.sif}
FAMS=(mtus8 mtige8 mtus6 mtige6 mtus4 mtige4 mtus3 mtige3 mtus2 mtige2 mtus1 mtige1)
NREP=3
i=${SLURM_ARRAY_TASK_ID:?SLURM array required}
r=$(( i % NREP + 1 )); f=${FAMS[$(( i / NREP ))]}
key="${f}_4000_${r}"
fcsv="$W/res_params/fits_${key}.csv"
if [ -f "$fcsv" ] && ! grep -q "License error\|licenses in use\|workspace" "$fcsv"; then echo "already done: $key"; exit 0; fi
NT=${SLURM_CPUS_PER_TASK:-6}
mkdir -p "$W/res_params" "$W/tmp" "$W/cache"
echo "=== $(hostname) | ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} | $key | $(date) ==="
apptainer exec --nv \
  --env RENV_CONFIG_AUTOLOADER_ENABLED=FALSE --env R_PROFILE_USER=/dev/null \
  --env R_LIBS="$W/rlib" --env TMPDIR="$W/tmp" --env REMLAX_CACHE="$W/cache" \
  --env OMP_NUM_THREADS=$NT --env OPENBLAS_NUM_THREADS=$NT --env MKL_NUM_THREADS=$NT \
  --env XLA_PYTHON_CLIENT_PREALLOCATE=false --env JAX_ENABLE_X64=1 \
  "$SIF" Rscript "$W/remlax-validation/banc_moteurs/R/executer.R" \
    --family=$f --n=4000 --rep=$r --out="$W/res_params" --timeout=${TIMEOUT:-1800} \
    --routes=D-CPU,D-GPU,D-CPU-LBFGS,D-GPU-LBFGS,M-dense-CPU,M-dense-GPU,M-sparse-CPU --external=asreml
echo "=== end $(date) | code $? ==="
