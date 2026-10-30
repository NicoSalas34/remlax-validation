#!/bin/bash -l
# remlax benchmark (frozen version) against kmmr on genotoul: one line of
# tasks.txt (or of the stored taches.txt) per array task.
#   CPU : sbatch --array=1-N%10 --partition=workq --cpus-per-task=8 --mem=64G --export=ALL,TARGET=cpu banc_g.sbatch
#   GPU : sbatch --array=1-N%2 --partition=gpuq --gres=gpu:nvidia_a100:1 --cpus-per-task=8 --mem=64G --export=ALL,TARGET=gpu banc_g.sbatch
# A remlax "mme" task runs TWICE in a row: "cold" (empty compilation cache),
# then "warm" (programs read back from disk).
#SBATCH --job-name=bkg
#SBATCH --time=12:00:00
#SBATCH --output=/work/user/nsalas/kmmr_banc/logs/%x-%A_%a.out
module load containers/Apptainer/1.4.1
K=/work/user/nsalas/kmmr_banc
# SIF, RLIB_REMLAX, TASKS and OUT_TAG select the image, the frozen remlax, the task
# list and the suffix of the result files (defaults: the runs of 2026-10).
SIF=${SIF:-$K/ige_reml_v2.sif}
RLIB_REMLAX=${RLIB_REMLAX:-$K/rlib_remlax}
TASKS=${TASKS:-$K/tasks.txt}; [ -f $TASKS ] || TASKS=$K/taches.txt
read c o v t < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $TASKS)
NC=${SLURM_CPUS_PER_TASK:-8}
ID=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
CACHE=$K/tmp/jax_$ID; RC=$K/tmp/rc_$ID
rm -rf $CACHE $RC
NVL=/usr/local/lib/python3.12/dist-packages/nvidia
LDP="$NVL/cuda_runtime/lib:$NVL/cublas/lib:$NVL/cusolver/lib:$NVL/cusparse/lib:$NVL/nvjitlink/lib:$NVL/cuda_nvrtc/lib:/.singularity.d/libs"
E="--env LD_LIBRARY_PATH=$LDP --env RENV_CONFIG_AUTOLOADER_ENABLED=FALSE --env R_PROFILE_USER=/dev/null --env R_LIBS=$RLIB_REMLAX:$K/rlib --env OMP_NUM_THREADS=$NC --env OPENBLAS_NUM_THREADS=$NC --env TMPDIR=$K/tmp --env REMLAX_JAX_CACHE=$CACHE --env REMLAX_CACHE=$RC"
if [ "$TARGET" = gpu ]; then NV=--nv; G=1; E="$E --env XLA_PYTHON_CLIENT_MEM_FRACTION=0.95"; else NV=; G=0; E="$E --env JAX_PLATFORMS=cpu"; fi
A="--case $c --tool $o --gpu $G --kmmr $K/src/kmmr --out $K/res/bench_${TARGET}${OUT_TAG:+_$OUT_TAG}.csv"
[ "$v" != "-" ] && A="$A --route $v"
[ "$t" != "-" ] && A="$A --tol $t"
echo "== task $SLURM_ARRAY_TASK_ID: $c $o $v $t ($TARGET, $NC threads, $(hostname)) $(date -Is)"
[ -n "$NV" ] && nvidia-smi -L
run() { apptainer exec $NV $E $SIF Rscript $K/banc/R/executer.R $A --label $1; }
if [ "$o" = remlax ] && [ "$v" = mme ]; then
  run cold
  run warm
else
  run single
fi
rm -rf $CACHE $RC
echo "== end $(date -Is)"
