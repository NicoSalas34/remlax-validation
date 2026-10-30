#!/bin/bash -l
# Banc remlax (version figee) contre kmmr sur genotoul : une ligne de taches.txt
# par tache du tableau.
#   CPU : sbatch --array=1-N%10 --partition=workq --cpus-per-task=8 --mem=64G --export=ALL,CIBLE=cpu banc_g.sbatch
#   GPU : sbatch --array=1-N%2 --partition=gpuq --gres=gpu:nvidia_a100:1 --cpus-per-task=8 --mem=64G --export=ALL,CIBLE=gpu banc_g.sbatch
# Une tache remlax "mme" tourne DEUX fois de suite : "froid" (cache de
# compilation vide) puis "chaud" (programmes relus sur disque).
#SBATCH --job-name=bkg
#SBATCH --time=12:00:00
#SBATCH --output=/work/user/nsalas/kmmr_banc/logs/%x-%A_%a.out
module load containers/Apptainer/1.4.1
K=/work/user/nsalas/kmmr_banc
read c o v t < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $K/taches.txt)
NC=${SLURM_CPUS_PER_TASK:-8}
ID=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
CACHE=$K/tmp/jax_$ID; RC=$K/tmp/rc_$ID
rm -rf $CACHE $RC
NVL=/usr/local/lib/python3.12/dist-packages/nvidia
LDP="$NVL/cuda_runtime/lib:$NVL/cublas/lib:$NVL/cusolver/lib:$NVL/cusparse/lib:$NVL/nvjitlink/lib:$NVL/cuda_nvrtc/lib:/.singularity.d/libs"
E="--env LD_LIBRARY_PATH=$LDP --env RENV_CONFIG_AUTOLOADER_ENABLED=FALSE --env R_PROFILE_USER=/dev/null --env R_LIBS=$K/${RLIB:-rlib_remlax}:$K/rlib --env OMP_NUM_THREADS=$NC --env OPENBLAS_NUM_THREADS=$NC --env TMPDIR=$K/tmp --env REMLAX_JAX_CACHE=$CACHE --env REMLAX_CACHE=$RC"
if [ "$CIBLE" = gpu ]; then NV=--nv; G=1; E="$E --env XLA_PYTHON_CLIENT_MEM_FRACTION=0.95"; else NV=; G=0; E="$E --env JAX_PLATFORMS=cpu"; fi
A="--cas $c --outil $o --gpu $G --kmmr $K/src/kmmr --sortie $K/res/banc_${CIBLE}${SUFFIXE}.csv"
[ "$v" != "-" ] && A="$A --voie $v"
[ "$t" != "-" ] && A="$A --tol $t"
echo "== tache $SLURM_ARRAY_TASK_ID : $c $o $v $t ($CIBLE, $NC fils, $(hostname)) $(date -Is)"
[ -n "$NV" ] && nvidia-smi -L
lancer() { apptainer exec $NV $E $K/ige_reml_v2.sif Rscript $K/banc/R/executer.R $A --etiquette $1; }
if [ "$o" = remlax ] && [ "$v" = mme ]; then
  lancer froid
  lancer chaud
else
  lancer unique
fi
rm -rf $CACHE $RC
echo "== fin $(date -Is)"
