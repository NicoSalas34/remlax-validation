#!/bin/bash
# ==============================================================================
# IEE benchmarks: family mtigee<t> (DGE/IGE us(2t) + IEE at plant level us(t)
# + residual us(t) per plant; R/familles_params.R). m is of the order of n and
# C nearly full, as in mvD7.
#   tasks 0-17  : number of parameters, t = 8,6,4,3,2,1 at n = 4000, 3 replicates
#   tasks 18-57 : size, mtigee3 at n = 8000..500, 5 replicates
# Every remlax route (dense AI and L-BFGS-B, CPU and GPU; MME dense CPU/GPU;
# MME sparse) and asreml (n <= 4000), MIG slice 1g.10gb.
#
#   sbatch --array=0-57%6 $W/remlax-validation/banc_moteurs/banc_iee.sbatch
# ==============================================================================
#SBATCH --job-name=bench_iee
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
i=${SLURM_ARRAY_TASK_ID:?SLURM array required}
if [ $i -lt 18 ]; then
  FAMS=(mtigee8 mtigee6 mtigee4 mtigee3 mtigee2 mtigee1)
  r=$(( i % 3 + 1 )); f=${FAMS[$(( i / 3 ))]}; n=4000; out="$W/res_iee_params"
else
  j=$(( i - 18 )); SIZES=(8000 6000 4000 3000 2000 1500 1000 500)
  r=$(( j % 5 + 1 )); n=${SIZES[$(( j / 5 ))]}; f=mtigee3; out="$W/res_iee_size"
fi
EXT=asreml; [ $n -gt 4000 ] && EXT=none
key="${f}_${n}_${r}"
fcsv="$out/fits_${key}.csv"
if [ -f "$fcsv" ] && ! grep -q "License error\|licenses in use\|workspace" "$fcsv"; then echo "already done: $key"; exit 0; fi
NT=${SLURM_CPUS_PER_TASK:-6}
mkdir -p "$out" "$W/tmp" "$W/cache"
echo "=== $(hostname) | ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} | $key | $(date) ==="
apptainer exec --nv \
  --env RENV_CONFIG_AUTOLOADER_ENABLED=FALSE --env R_PROFILE_USER=/dev/null \
  --env R_LIBS="$W/rlib" --env TMPDIR="$W/tmp" --env REMLAX_CACHE="$W/cache" \
  --env OMP_NUM_THREADS=$NT --env OPENBLAS_NUM_THREADS=$NT --env MKL_NUM_THREADS=$NT \
  --env XLA_PYTHON_CLIENT_PREALLOCATE=false --env JAX_ENABLE_X64=1 \
  "$SIF" Rscript "$W/remlax-validation/banc_moteurs/R/executer.R" \
    --family=$f --n=$n --rep=$r --out="$out" --timeout=${TIMEOUT:-1800} \
    --routes=D-CPU,D-GPU,D-CPU-LBFGS,D-GPU-LBFGS,M-dense-CPU,M-dense-GPU,M-sparse-CPU --external=$EXT
echo "=== end $(date) | code $? ==="
