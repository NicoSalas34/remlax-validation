#!/bin/bash -l
# Absolute benchmark (simulated cases): one line of tasks<LIST>_<TARGET>.txt per
# task (the stored lists are named taches<LIST>_<TARGET>.txt and are read when
# the English one is missing).
#   CPU : sbatch --array=1-N%12 --partition=workq --cpus-per-task=8 --mem=96G --export=ALL,TARGET=cpu banc_abs.sbatch
#   GPU : sbatch --array=1-N%2 --partition=gpuq --gres=gpu:nvidia_a100:1 --cpus-per-task=8 --mem=96G --export=ALL,TARGET=gpu banc_abs.sbatch
# remlax mme runs twice: "cold" (empty compilation cache), then "warm".
# Each fit is bounded by TIMEOUT seconds. A timeout is written to
# res/abs_<TARGET>_timeouts.csv.
#SBATCH --job-name=babs
#SBATCH --time=14:00:00
#SBATCH --output=/work/user/nsalas/kmmr_banc/logs/%x-%A_%a.out
module load containers/Apptainer/1.4.1
K=/work/user/nsalas/kmmr_banc
TIMEOUT=${TIMEOUT:-21600}
TASKS=$K/absolu/tasks${LIST}_${TARGET}.txt; [ -f $TASKS ] || TASKS=$K/absolu/taches${LIST}_${TARGET}.txt
read c o v t rep < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $TASKS)
P=${rep:+${rep}_}   # replicate: prefix of the label (empty for the first pass)
NC=${SLURM_CPUS_PER_TASK:-8}
ID=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
CACHE=$K/tmp/jax_$ID; RC=$K/tmp/rc_$ID
rm -rf $CACHE $RC
NVL=/usr/local/lib/python3.12/dist-packages/nvidia
LDP="$NVL/cuda_runtime/lib:$NVL/cublas/lib:$NVL/cusolver/lib:$NVL/cusparse/lib:$NVL/nvjitlink/lib:$NVL/cuda_nvrtc/lib:/.singularity.d/libs"
E="--env LD_LIBRARY_PATH=$LDP --env RENV_CONFIG_AUTOLOADER_ENABLED=FALSE --env R_PROFILE_USER=/dev/null --env R_LIBS=$K/${RLIB:-rlib_remlax_v10}:$K/rlib --env OMP_NUM_THREADS=$NC --env OPENBLAS_NUM_THREADS=$NC --env TMPDIR=$K/tmp --env REMLAX_JAX_CACHE=$CACHE --env REMLAX_CACHE=$RC"
if [ "$TARGET" = gpu ]; then NV=--nv; G=1; E="$E --env XLA_PYTHON_CLIENT_MEM_FRACTION=0.95"; else NV=; G=0; E="$E --env JAX_PLATFORMS=cpu"; fi
OUT=$K/res/abs_${TARGET}${SUFFIX}.csv
A="--case $c --tool $o --gpu $G --kmmr $K/src/kmmr --sim $K/absolu --out $OUT"
[ "$v" != "-" ] && A="$A --route $v"
[ "$t" != "-" ] && A="$A --tol $t"
echo "== task $SLURM_ARRAY_TASK_ID: $c $o $v $t ($TARGET, $NC threads, $(hostname)) $(date -Is)"
[ -n "$NV" ] && nvidia-smi -L
run() {
  timeout $TIMEOUT apptainer exec $NV $E $K/ige_reml_v2.sif Rscript $K/banc/R/executer.R $A --label ${P}$1
  rc=$?
  if [ $rc -eq 124 ]; then
    F=$K/res/abs_${TARGET}${SUFFIX}_timeouts.csv
    [ -f $F ] || echo "label,case,tool,route,gpu,node,timeout,date" > $F
    echo "${P}$1,$c,$o,${v/-/ai_tol$t},$G,$(hostname),$TIMEOUT,$(date '+%Y-%m-%d %H:%M:%S')" >> $F
    echo "== TIMEOUT EXCEEDED ($TIMEOUT s)"
  fi
  return $rc
}
if [ "$o" = remlax ] && [ "$v" = mme ]; then
  run cold && run warm
else
  run single
fi
rm -rf $CACHE $RC
echo "== end $(date -Is)"
