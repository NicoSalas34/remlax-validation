#!/bin/bash
# ==============================================================================
# Benchmark of the remlax engines against asreml, sommer and lme4, on a MIG
# slice 1g.10gb.
#
# One task per case (family, n, replicate): 9 families x 8 sizes x 5
# replicates = 360 tasks. Each task fits the five remlax routes (dense and
# MME, CPU and GPU) then the three external programs, on the SAME data (seed
# a function of family, n and replicate). Every CPU route runs on the cores
# allocated to the slice, so CPU and GPU share the machine.
#
#   W=/scratch/users/salasn/banc_moteurs_20261022
#   sbatch --array=0-359%7 $W/remlax-validation/banc_moteurs/banc_mig.sbatch
#
# The largest sizes go first (low indices), so that the tail of the array is
# made of short cases.
# ==============================================================================
#SBATCH --job-name=bench_eng
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
FAMILIES=(iid grm us ped ige field resid exp fa)
SIZES=(8000 6000 4000 3000 2000 1500 1000 500)
NREP=5

i=${SLURM_ARRAY_TASK_ID:?SLURM array required}
nf=${#FAMILIES[@]}
r=$(( i % NREP + 1 )); i=$(( i / NREP ))
f=${FAMILIES[$(( i % nf ))]}; i=$(( i / nf ))
n=${SIZES[$i]}
key="${f}_${n}_${r}"
# case complete and without licence failure: nothing to do. Otherwise
# executer.R resumes what is missing (.partial files, marker running_<case>)
fcsv="$W/res/fits_${key}.csv"
if [ -f "$fcsv" ] && ! grep -q "License error\|licenses in use\|workspace" "$fcsv"; then echo "already done: $key"; exit 0; fi

NT=${SLURM_CPUS_PER_TASK:-6}
mkdir -p "$W/res" "$W/tmp" "$W/cache"
echo "=== $(hostname) | ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} | $key | $NT cores | $(date) ==="
nvidia-smi -L

# renv disabled (the project .Rprofile would point .libPaths back to the host);
# one thread count for all software; JAX allocator without preallocation (the
# slice has only 10 GB).
apptainer exec --nv \
  --env RENV_CONFIG_AUTOLOADER_ENABLED=FALSE --env R_PROFILE_USER=/dev/null \
  --env R_LIBS="$W/rlib" --env TMPDIR="$W/tmp" --env REMLAX_CACHE="$W/cache" \
  --env OMP_NUM_THREADS=$NT --env OPENBLAS_NUM_THREADS=$NT --env MKL_NUM_THREADS=$NT \
  --env XLA_PYTHON_CLIENT_PREALLOCATE=false --env JAX_ENABLE_X64=1 \
  "$SIF" Rscript "$W/remlax-validation/banc_moteurs/R/executer.R" \
    --family=$f --n=$n --rep=$r --out="$W/res" --timeout=${TIMEOUT:-1800}
rc=$?
echo "=== end $(date) | code $rc ==="
exit $rc
