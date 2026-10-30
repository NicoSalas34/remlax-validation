#!/bin/bash
# ==============================================================================
# Long-running multi-trait benchmarks (2026-10-23, corrected design).
#
# NUMBER OF PARAMETERS at a fixed number of UNITS: n = 1000 t, hence 1000 plants
# and 200 genotypes (250 for mtus) whatever t. At fixed n, the number of
# genotypes fell with t (100 at t = 8 for 136 parameters in G). The small
# eigenvalues of G (0.05) were then no longer estimable, the REML optimum fell
# on the boundary (rank 12 of 16) and the benchmark measured poor
# identification as much as the number of parameters.
#   cases 0-35  : mtigee, mtige, mtus x t = 8,6,4,3,2,1 x 2 replicates, n = 1000 t
#   cases 36-50 : size, mtigee3 at n = 8000,4000,2000,1000,500, 3 replicates
# ONE TASK PER (CASE, ROUTE), time limit of 48 h per fit. Dense CPU routes
# removed for t >= 6 (more than 9 h per fit on a MIG slice), and MME CPU routes
# removed for mtigee6 and mtigee8 (more than 8 to 12 h).
# Outputs: res_long/<route>/.
#
#   sbatch --array=0-407%7 $W/remlax-validation/banc_moteurs/banc_long.sbatch
# ==============================================================================
#SBATCH --job-name=bench_long
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
ROUTES=(D-CPU D-GPU D-CPU-LBFGS D-GPU-LBFGS M-dense-CPU M-dense-GPU M-sparse-CPU asreml)
i=${SLURM_ARRAY_TASK_ID:?SLURM array required}
c=$(( i / ${#ROUTES[@]} )); v=${ROUTES[$(( i % ${#ROUTES[@]} ))]}
if [ $c -lt 36 ]; then
  TS=(8 6 4 3 2 1); FS=(mtigee mtige mtus)
  r=$(( c % 2 + 1 )); k=$(( c / 2 )); f=${FS[$(( k % 3 ))]}; t=${TS[$(( k / 3 ))]}
  f="$f$t"; n=$(( 1000 * t ))
else
  j=$(( c - 36 )); SIZES=(8000 4000 2000 1000 500)
  r=$(( j % 3 + 1 )); n=${SIZES[$(( j / 3 ))]}; f=mtigee3; t=3
fi
case "$v" in
  D-CPU|D-CPU-LBFGS) [ $t -ge 6 ] && { echo "dense CPU route removed for $f (t >= 6)"; exit 0; } ;;
  M-dense-CPU|M-sparse-CPU)
    # mtigee (IEE at plant level): C nearly full, m ~ n; more than 8 to 12 h
    # per fit on CPU for t >= 6 (2026-10-24), removed
    case "$f" in mtigee6|mtigee8) echo "MME CPU route removed for $f"; exit 0 ;; esac ;;
esac
if [ "$v" = "asreml" ]; then ARGS="--routes=none --external=asreml"; else ARGS="--routes=$v --external=none"; fi
out="$W/res_long/$v"
key="${f}_${n}_${r}"
fcsv="$out/fits_${key}.csv"
if [ -f "$fcsv" ] && ! grep -q "License error\|licenses in use\|workspace" "$fcsv"; then echo "already done: $key $v"; exit 0; fi
NT=${SLURM_CPUS_PER_TASK:-6}
mkdir -p "$out" "$W/tmp" "$W/cache"
echo "=== $(hostname) | ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} | $key | $v | $(date) ==="
apptainer exec --nv \
  --env RENV_CONFIG_AUTOLOADER_ENABLED=FALSE --env R_PROFILE_USER=/dev/null \
  --env R_LIBS="$W/rlib_mme2" --env TMPDIR="$W/tmp" --env REMLAX_CACHE="$W/cache" \
  --env OMP_NUM_THREADS=$NT --env OPENBLAS_NUM_THREADS=$NT --env MKL_NUM_THREADS=$NT \
  --env XLA_PYTHON_CLIENT_PREALLOCATE=false --env JAX_ENABLE_X64=1 \
  "$SIF" Rscript "$W/remlax-validation/banc_moteurs/R/executer.R" \
    --family=$f --n=$n --rep=$r --out="$out" --timeout=${TIMEOUT:-172800} $ARGS
echo "=== end $(date) | code $? ==="
