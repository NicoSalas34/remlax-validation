#!/bin/bash
# ==============================================================================
# Long-running IEE benchmarks: family mtigee<t> (DGE/IGE us(2t) + IEE at
# plant level us(t) + residual us(t)). ONE TASK PER (CASE, ROUTE): the routes
# of a case run in parallel instead of adding up. Time limit of 48 h per fit.
#   cases 0-11  : number of parameters, t = 8,6,4,3,2,1 at n = 4000, 2 replicates
#   cases 12-26 : size, mtigee3 at n = 8000,4000,2000,1000,500, 3 replicates
#   routes      : D-CPU D-GPU D-CPU-LBFGS D-GPU-LBFGS M-dense-CPU M-dense-GPU
#                 M-sparse-CPU asreml (asreml up to n = 4000)
# Outputs: res_iee_long/<route>/fits_<case>.csv (one row per file).
#
#   sbatch --array=0-215%7 $W/remlax-validation/banc_moteurs/banc_iee_long.sbatch
# ==============================================================================
#SBATCH --job-name=iee_long
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
ROUTES=(D-CPU D-GPU D-CPU-LBFGS D-GPU-LBFGS M-dense-CPU M-dense-GPU M-sparse-CPU asreml)
i=${SLURM_ARRAY_TASK_ID:?SLURM array required}
c=$(( i / ${#ROUTES[@]} )); v=${ROUTES[$(( i % ${#ROUTES[@]} ))]}
if [ $c -lt 12 ]; then
  FAMS=(mtigee8 mtigee6 mtigee4 mtigee3 mtigee2 mtigee1)
  r=$(( c % 2 + 1 )); f=${FAMS[$(( c / 2 ))]}; n=4000
else
  j=$(( c - 12 )); SIZES=(8000 4000 2000 1000 500)
  r=$(( j % 3 + 1 )); n=${SIZES[$(( j / 3 ))]}; f=mtigee3
fi
if [ "$v" = "asreml" ]; then
  [ $n -gt 4000 ] && { echo "asreml not run beyond n = 4000: $f $n"; exit 0; }
  ARGS="--routes=none --external=asreml"
else
  ARGS="--routes=$v --external=none"
fi
out="$W/res_iee_long/$v"
key="${f}_${n}_${r}"
fcsv="$out/fits_${key}.csv"
if [ -f "$fcsv" ] && ! grep -q "License error\|licenses in use\|workspace" "$fcsv"; then echo "already done: $key $v"; exit 0; fi
NT=${SLURM_CPUS_PER_TASK:-6}
mkdir -p "$out" "$W/tmp" "$W/cache"
echo "=== $(hostname) | ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} | $key | $v | $(date) ==="
apptainer exec --nv \
  --env RENV_CONFIG_AUTOLOADER_ENABLED=FALSE --env R_PROFILE_USER=/dev/null \
  --env R_LIBS="$W/rlib" --env TMPDIR="$W/tmp" --env REMLAX_CACHE="$W/cache" \
  --env OMP_NUM_THREADS=$NT --env OPENBLAS_NUM_THREADS=$NT --env MKL_NUM_THREADS=$NT \
  --env XLA_PYTHON_CLIENT_PREALLOCATE=false --env JAX_ENABLE_X64=1 \
  "$SIF" Rscript "$W/remlax-validation/banc_moteurs/R/executer.R" \
    --family=$f --n=$n --rep=$r --out="$out" --timeout=${TIMEOUT:-172800} $ARGS
echo "=== end $(date) | code $? ==="
