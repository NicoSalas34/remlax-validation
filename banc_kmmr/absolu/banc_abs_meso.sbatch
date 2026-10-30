#!/bin/bash
# Banc absolu sur meso, ou se trouve la licence asreml. Meme harnais, memes jeux.
#   asreml (licence : en serie) : sbatch --array=1-30%1 --export=ALL,LISTE=asreml banc_abs_meso.sbatch
#   kmmr + remlax CPU           : sbatch --array=1-N%8 --export=ALL,LISTE=meso banc_abs_meso.sbatch
#SBATCH --job-name=babs_m
#SBATCH --partition=cpu-dedicated
#SBATCH --account=dedicated-cpu@cirad-long
#SBATCH --qos=dedicated
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=1-00:00:00
#SBATCH --output=/scratch/users/salasn/kmmr_banc/logs/%x-%A_%a.out
K=/scratch/users/salasn/kmmr_banc
DELAI=${DELAI:-21600}
read c o v t < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $K/absolu/taches_${LISTE}.txt)
ID=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
CACHE=$K/tmp/jax_$ID; RC=$K/tmp/rc_$ID
mkdir -p $K/tmp; rm -rf $CACHE $RC
E="--env RENV_CONFIG_AUTOLOADER_ENABLED=FALSE --env R_PROFILE_USER=/dev/null --env R_LIBS=$K/rlib_remlax_v10:$K/rlib --env OMP_NUM_THREADS=8 --env OPENBLAS_NUM_THREADS=8 --env JAX_PLATFORMS=cpu --env TMPDIR=$K/tmp --env REMLAX_JAX_CACHE=$CACHE --env REMLAX_CACHE=$RC"
SORTIE=$K/res/abs_meso.csv
A="--cas $c --outil $o --gpu 0 --kmmr $K/src/kmmr --sim $K/absolu --sortie $SORTIE"
[ "$v" != "-" ] && A="$A --voie $v"
[ "$t" != "-" ] && A="$A --tol $t"
echo "== tache $SLURM_ARRAY_TASK_ID : $c $o $v $t (meso, $(hostname)) $(date -Is)"
lancer() {
  timeout $DELAI apptainer exec $E /scratch/users/salasn/ige_reml_v2.sif Rscript $K/banc/R/executer.R $A --etiquette $1
  rc=$?
  if [ $rc -eq 124 ]; then
    F=$K/res/abs_meso_delais.csv
    [ -f $F ] || echo "etiquette,cas,outil,voie,gpu,noeud,delai,date" > $F
    echo "$1,$c,$o,${v/-/ai},0,$(hostname),$DELAI,$(date '+%Y-%m-%d %H:%M:%S')" >> $F
    echo "== DELAI DEPASSE ($DELAI s)"
  fi
  return $rc
}
if [ "$o" = remlax ] && [ "$v" = mme ]; then
  lancer froid && lancer chaud
else
  lancer unique
fi
rm -rf $CACHE $RC
echo "== fin $(date -Is)"
