#!/bin/bash
# Absolute benchmark on meso, where the asreml licence lives. Same harness, same data sets.
#   asreml (licence: in series) : sbatch --array=1-30%1 --export=ALL,LIST=asreml banc_abs_meso.sbatch
#   kmmr + remlax CPU           : sbatch --array=1-N%8 --export=ALL,LIST=meso banc_abs_meso.sbatch
#SBATCH --job-name=babs_m
#SBATCH --partition=cpu-dedicated
#SBATCH --account=dedicated-cpu@cirad-long
#SBATCH --qos=dedicated
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=1-00:00:00
#SBATCH --output=/scratch/users/salasn/kmmr_banc/logs/%x-%A_%a.out
K=/scratch/users/salasn/kmmr_banc
TIMEOUT=${TIMEOUT:-21600}
TASKS=$K/absolu/tasks_${LIST}.txt; [ -f $TASKS ] || TASKS=$K/absolu/taches_${LIST}.txt
read c o v t < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $TASKS)
ID=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
CACHE=$K/tmp/jax_$ID; RC=$K/tmp/rc_$ID
mkdir -p $K/tmp; rm -rf $CACHE $RC
E="--env RENV_CONFIG_AUTOLOADER_ENABLED=FALSE --env R_PROFILE_USER=/dev/null --env R_LIBS=$K/rlib_remlax_v10:$K/rlib --env OMP_NUM_THREADS=8 --env OPENBLAS_NUM_THREADS=8 --env JAX_PLATFORMS=cpu --env TMPDIR=$K/tmp --env REMLAX_JAX_CACHE=$CACHE --env REMLAX_CACHE=$RC"
OUT=$K/res/abs_meso.csv
A="--case $c --tool $o --gpu 0 --kmmr $K/src/kmmr --sim $K/absolu --out $OUT"
[ "$v" != "-" ] && A="$A --route $v"
[ "$t" != "-" ] && A="$A --tol $t"
echo "== task $SLURM_ARRAY_TASK_ID: $c $o $v $t (meso, $(hostname)) $(date -Is)"
run() {
  timeout $TIMEOUT apptainer exec $E /scratch/users/salasn/ige_reml_v2.sif Rscript $K/banc/R/executer.R $A --label $1
  rc=$?
  if [ $rc -eq 124 ]; then
    F=$K/res/abs_meso_timeouts.csv
    [ -f $F ] || echo "label,case,tool,route,gpu,node,timeout,date" > $F
    echo "$1,$c,$o,${v/-/ai},0,$(hostname),$TIMEOUT,$(date '+%Y-%m-%d %H:%M:%S')" >> $F
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
