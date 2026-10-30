#!/bin/bash
# Absolute benchmark on meso, on a WHOLE RESERVED NODE. The shared meso nodes
# are overloaded by unconfined jobs (load 300 for 192 cores), and the same fit
# varied there by a factor of 1 to 40 (A/B trial of 2026-10-27).
# One case per task: asreml, kmmr, remlax in turn, three replicates each, all
# pinned to the same 8 cores of one socket (taskset).
#   sbatch --array=1-30%1 banc_abs_meso_excl.sbatch      (asreml licence: in series)
#SBATCH --job-name=babs_x
#SBATCH --partition=cpu-dedicated
#SBATCH --account=dedicated-cpu@cirad-long
#SBATCH --qos=dedicated
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=1-00:00:00
#SBATCH --output=/scratch/users/salasn/kmmr_banc/logs/%x-%A_%a.out
K=/scratch/users/salasn/kmmr_banc
TIMEOUT=${TIMEOUT:-21600}
CASES=$K/absolu/cases_excl.txt; [ -f $CASES ] || CASES=$K/absolu/cas_excl.txt
c=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $CASES)
CORES=0-7
OUT=$K/res/abs_meso_x.csv
echo "== $c on $(hostname), reserved, cores $CORES, load $(cut -d' ' -f1 /proc/loadavg) $(date -Is)"
run() {  # tool route tol label
  CACHE=$K/tmp/jaxx_${SLURM_JOB_ID}_$4; RC=$K/tmp/rcx_${SLURM_JOB_ID}_$4; rm -rf $CACHE $RC
  E="--env RENV_CONFIG_AUTOLOADER_ENABLED=FALSE --env R_PROFILE_USER=/dev/null --env R_LIBS=$K/rlib_remlax_v11:$K/rlib --env OMP_NUM_THREADS=8 --env OPENBLAS_NUM_THREADS=8 --env JAX_PLATFORMS=cpu --env TMPDIR=$K/tmp --env REMLAX_JAX_CACHE=$CACHE --env REMLAX_CACHE=$RC"
  A="--case $c --tool $1 --gpu 0 --kmmr $K/src/kmmr --sim $K/absolu --out $OUT --label $4"
  [ "$2" != "-" ] && A="$A --route $2"; [ "$3" != "-" ] && A="$A --tol $3"
  timeout $TIMEOUT taskset -c $CORES apptainer exec $E /scratch/users/salasn/ige_reml_v2.sif Rscript $K/banc/R/executer.R $A
  rc=$?
  if [ $rc -eq 124 ]; then
    F=$K/res/abs_meso_x_timeouts.csv; [ -f $F ] || echo "label,case,tool,route,gpu,node,timeout,date" > $F
    echo "$4,$c,$1,${2/-/ai},0,$(hostname),$TIMEOUT,$(date '+%Y-%m-%d %H:%M:%S')" >> $F
  fi
  echo "   load after $4: $(cut -d' ' -f1 /proc/loadavg)"
  return $rc
}
for r in r1 r2 r3; do
  run asreml - - ${r}_single || [ $r != r1 ] || true
  run kmmr - 0.000001 ${r}_single
  run remlax mme - ${r}_cold && run remlax mme - ${r}_warm
done
echo "== end $(date -Is)"
