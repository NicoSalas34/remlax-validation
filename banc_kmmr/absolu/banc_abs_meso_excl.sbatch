#!/bin/bash
# Banc absolu sur meso, NOEUD RESERVE EN ENTIER. Les noeuds partages de meso
# sont surcharges par des jobs non confines (charge 300 pour 192 coeurs) : le
# meme ajustement y variait d'un facteur 1 a 40 (essai A/B du 2026-10-27).
# Un cas par tache : asreml, kmmr, remlax tour a tour, trois repetitions
# chacun, tous epingles sur les memes 8 coeurs d'un socket (taskset).
#   sbatch --array=1-30%1 banc_abs_meso_excl.sbatch      (licence asreml : en serie)
#SBATCH --job-name=babs_x
#SBATCH --partition=cpu-dedicated
#SBATCH --account=dedicated-cpu@cirad-long
#SBATCH --qos=dedicated
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=1-00:00:00
#SBATCH --output=/scratch/users/salasn/kmmr_banc/logs/%x-%A_%a.out
K=/scratch/users/salasn/kmmr_banc
DELAI=${DELAI:-21600}
c=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $K/absolu/cas_excl.txt)
COEURS=0-7
SORTIE=$K/res/abs_meso_x.csv
echo "== $c sur $(hostname), reserve, coeurs $COEURS, charge $(cut -d' ' -f1 /proc/loadavg) $(date -Is)"
lancer() {  # outil voie tol etiquette
  CACHE=$K/tmp/jaxx_${SLURM_JOB_ID}_$4; RC=$K/tmp/rcx_${SLURM_JOB_ID}_$4; rm -rf $CACHE $RC
  E="--env RENV_CONFIG_AUTOLOADER_ENABLED=FALSE --env R_PROFILE_USER=/dev/null --env R_LIBS=$K/rlib_remlax_v11:$K/rlib --env OMP_NUM_THREADS=8 --env OPENBLAS_NUM_THREADS=8 --env JAX_PLATFORMS=cpu --env TMPDIR=$K/tmp --env REMLAX_JAX_CACHE=$CACHE --env REMLAX_CACHE=$RC"
  A="--cas $c --outil $1 --gpu 0 --kmmr $K/src/kmmr --sim $K/absolu --sortie $SORTIE --etiquette $4"
  [ "$2" != "-" ] && A="$A --voie $2"; [ "$3" != "-" ] && A="$A --tol $3"
  timeout $DELAI taskset -c $COEURS apptainer exec $E /scratch/users/salasn/ige_reml_v2.sif Rscript $K/banc/R/executer.R $A
  rc=$?
  if [ $rc -eq 124 ]; then
    F=$K/res/abs_meso_x_delais.csv; [ -f $F ] || echo "etiquette,cas,outil,voie,gpu,noeud,delai,date" > $F
    echo "$4,$c,$1,${2/-/ai},0,$(hostname),$DELAI,$(date '+%Y-%m-%d %H:%M:%S')" >> $F
  fi
  echo "   charge apres $4 : $(cut -d' ' -f1 /proc/loadavg)"
  return $rc
}
for r in r1 r2 r3; do
  lancer asreml - - ${r}_unique || [ $r != r1 ] || true
  lancer kmmr - 0.000001 ${r}_unique
  lancer remlax mme - ${r}_froid && lancer remlax mme - ${r}_chaud
done
echo "== fin $(date -Is)"
