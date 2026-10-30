#!/bin/bash -l
# Writes the simulated data sets of the absolute benchmark to /work/user/nsalas/kmmr_banc/absolu/sim
#SBATCH --job-name=babs_prep
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=workq
#SBATCH --output=/work/user/nsalas/kmmr_banc/logs/%x-%j.out
module load containers/Apptainer/1.4.1
K=/work/user/nsalas/kmmr_banc
CASES=$K/absolu/cases_g45.txt; [ -f $CASES ] || CASES=$K/absolu/cas_g45.txt
mkdir -p $K/absolu/sim/copies_g1; for f in $K/absolu/sim/sim_*_t*_g[23].rds; do case $f in *_t1_*) ;; *) mv $f $K/absolu/sim/copies_g1/ ;; esac; done
apptainer exec --env RENV_CONFIG_AUTOLOADER_ENABLED=FALSE --env R_PROFILE_USER=/dev/null --env OMP_NUM_THREADS=8 --env OPENBLAS_NUM_THREADS=8 \
  $K/ige_reml_v2.sif Rscript -e "source('$K/banc/R/sim.R'); prepare_sim(readLines('$CASES'), '$K/absolu')"
