#!/bin/bash -l
# Banc absolu (cas simules) : une ligne de taches_<CIBLE>.txt par tache.
#   CPU : sbatch --array=1-N%12 --partition=workq --cpus-per-task=8 --mem=96G --export=ALL,CIBLE=cpu banc_abs.sbatch
#   GPU : sbatch --array=1-N%2 --partition=gpuq --gres=gpu:nvidia_a100:1 --cpus-per-task=8 --mem=96G --export=ALL,CIBLE=gpu banc_abs.sbatch
# remlax mme tourne deux fois : "froid" (cache de compilation vide) puis "chaud".
# Chaque ajustement est borne par DELAI secondes ; un depassement est ecrit
# dans res/abs_<CIBLE>_delais.csv.
#SBATCH --job-name=babs
#SBATCH --time=14:00:00
#SBATCH --output=/work/user/nsalas/kmmr_banc/logs/%x-%A_%a.out
module load containers/Apptainer/1.4.1
K=/work/user/nsalas/kmmr_banc
DELAI=${DELAI:-21600}
read c o v t rep < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $K/absolu/taches${LISTE}_${CIBLE}.txt)
P=${rep:+${rep}_}   # repetition : prefixe de l etiquette (vide pour la premiere passe)
NC=${SLURM_CPUS_PER_TASK:-8}
ID=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
CACHE=$K/tmp/jax_$ID; RC=$K/tmp/rc_$ID
rm -rf $CACHE $RC
NVL=/usr/local/lib/python3.12/dist-packages/nvidia
LDP="$NVL/cuda_runtime/lib:$NVL/cublas/lib:$NVL/cusolver/lib:$NVL/cusparse/lib:$NVL/nvjitlink/lib:$NVL/cuda_nvrtc/lib:/.singularity.d/libs"
E="--env LD_LIBRARY_PATH=$LDP --env RENV_CONFIG_AUTOLOADER_ENABLED=FALSE --env R_PROFILE_USER=/dev/null --env R_LIBS=$K/${RLIB:-rlib_remlax_v10}:$K/rlib --env OMP_NUM_THREADS=$NC --env OPENBLAS_NUM_THREADS=$NC --env TMPDIR=$K/tmp --env REMLAX_JAX_CACHE=$CACHE --env REMLAX_CACHE=$RC"
if [ "$CIBLE" = gpu ]; then NV=--nv; G=1; E="$E --env XLA_PYTHON_CLIENT_MEM_FRACTION=0.95"; else NV=; G=0; E="$E --env JAX_PLATFORMS=cpu"; fi
SORTIE=$K/res/abs_${CIBLE}${SUFFIXE}.csv
A="--cas $c --outil $o --gpu $G --kmmr $K/src/kmmr --sim $K/absolu --sortie $SORTIE"
[ "$v" != "-" ] && A="$A --voie $v"
[ "$t" != "-" ] && A="$A --tol $t"
echo "== tache $SLURM_ARRAY_TASK_ID : $c $o $v $t ($CIBLE, $NC fils, $(hostname)) $(date -Is)"
[ -n "$NV" ] && nvidia-smi -L
lancer() {
  timeout $DELAI apptainer exec $NV $E $K/ige_reml_v2.sif Rscript $K/banc/R/executer.R $A --etiquette ${P}$1
  rc=$?
  if [ $rc -eq 124 ]; then
    F=$K/res/abs_${CIBLE}${SUFFIXE}_delais.csv
    [ -f $F ] || echo "etiquette,cas,outil,voie,gpu,noeud,delai,date" > $F
    echo "${P}$1,$c,$o,${v/-/ai_tol$t},$G,$(hostname),$DELAI,$(date '+%Y-%m-%d %H:%M:%S')" >> $F
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
