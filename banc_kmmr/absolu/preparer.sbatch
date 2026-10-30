#!/bin/bash -l
# Ecrit les jeux simules du banc absolu dans /work/user/nsalas/kmmr_banc/absolu/sim
#SBATCH --job-name=babs_prep
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=workq
#SBATCH --output=/work/user/nsalas/kmmr_banc/logs/%x-%j.out
module load containers/Apptainer/1.4.1
K=/work/user/nsalas/kmmr_banc
apptainer exec --env RENV_CONFIG_AUTOLOADER_ENABLED=FALSE --env R_PROFILE_USER=/dev/null --env OMP_NUM_THREADS=8 --env OPENBLAS_NUM_THREADS=8 \
  $K/ige_reml_v2.sif Rscript -e "source('$K/banc/R/sim.R'); preparer_sim(readLines('$K/absolu/cas.txt'), '$K/absolu')"
