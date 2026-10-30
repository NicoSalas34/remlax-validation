#!/bin/bash
# =============================================================================
# FULL VALIDATION OF remlax ON GPU (A100 MIG slice)
#
# Everything runs in ONE image, ige_reml.sif, which holds both R (asreml, lme4,
# sommer, pbkrtest) and JAX CUDA. This makes the job possible. remlax has R
# call the solver, and two images would need apptainer nested in apptainer.
# Nesting is impossible here: the host binary does not find its libraries
# inside the container.
#
# WHAT IS CHECKED, in this order:
#   1. Python core             algebra, analytic REML, degenerate cases
#   2. cross-validation        same -2logL as lme4 and sommer, ON GPU
#   3. asreml (catalogue)      ar1, ar1xar1, 2D splines, str, vpredict, Wald
#   4. asreml (additions)      lvr, anisotropic mtrn, own, dsum, predict, K-R
#   5. full IGE model          the real case of the project
#   6. CPU/GPU parity          35 designs, same bundle, two machines
#
# Step 6 alone tests the MACHINE. The others test the code and simply run on
# the card.
#
#   sbatch scripts/tests/slurm_valid_gpu.slurm
# =============================================================================
#SBATCH --job-name="remlax_valid"
#SBATCH --output="logs/remlax_valid_%j.out"
#SBATCH --error="logs/remlax_valid_%j.out"
#SBATCH --partition=gpu-cirad-dedicated
#SBATCH --account=dedicated-gpu@cirad-long
#SBATCH --qos=dedicated
#SBATCH --gres=gpu:nvidia_a100_80gb_pcie_1g.10gb:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00

set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
mkdir -p logs output/tests

SIF=${SIF:-/scratch/users/salasn/ige_reml.sif}

# renv MUST BE DISABLED inside the container. The project has a .Rprofile that
# loads renv, and renv points .libPaths() to the HOST library, built for another
# version of R. The container then loses its own library and stops on
# "there is no package called 'jsonlite'", although everything is there.
export APPTAINERENV_RENV_CONFIG_AUTOLOADER_ENABLED=FALSE
RUN() { apptainer exec --nv "$SIF" "$@"; }

echo "=== $(hostname) | job ${SLURM_JOB_ID:-local} | $(date) ==="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
RUN python3 -c "import jax; print('[jax]', jax.__version__, jax.devices())"

status=0
step() {                        # step <number> <title> <command...>
  local n="$1" title="$2"; shift 2
  echo; echo "########## $n. $title ##########"
  if "$@"; then echo "--> $n OK"; else echo "--> $n FAIL (code $?)"; status=1; fi
}

# 1-5: the suites as they are, but with the GPU backend.
step 1 "core (Python)"              RUN python3 scripts/tests/test_remlax_core.py
step 2 "lme4/sommer cross (GPU)"  env APPTAINERENV_RK_BACKEND=gpu \
        apptainer exec --nv "$SIF" Rscript scripts/tests/test_remlax.R
step 3 "asreml: catalogue"          RUN Rscript scripts/tests/test_remlax_asreml.R
step 4 "asreml: 09/2026 additions"  env APPTAINERENV_RK_BACKEND=gpu \
        apptainer exec --nv "$SIF" Rscript scripts/tests/test_remlax_asreml2.R
step 5 "full IGE model"             RUN Rscript scripts/tests/test_remlax_ige.R

# 6: parity. R exports and Python fits. Both backends read EXACTLY the same
# bundle, which is the whole point of the on-file format.
#
# COUNT COMPILATIONS, NOT COMPUTATIONS. 34 models x 2 backends make 68 XLA
# compilations. On a MIG SLICE (7 SMs out of 108) the compiler is so slow that
# XLA itself prints "Very slow compile?". Compilation dominates the sweep, not
# the algebra: the models have 160 to 480 observations. Expect several tens of
# minutes here, and do not mistake it for a hang.
echo; echo "########## 6. CPU/GPU parity ##########"
RUN Rscript scripts/tests/export_bundles.R --out=output/tests/bundles || status=1
RUN python3 scripts/tests/parite_gpu.py output/tests/bundles || status=1

echo; echo "=== end $(date) | overall status: $([ $status -eq 0 ] && echo OK || echo FAIL) ==="
exit $status
