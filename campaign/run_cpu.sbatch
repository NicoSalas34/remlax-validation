#!/bin/bash
# =============================================================================
# Validation campaign of remlax, CPU part: every suite against asreml, lme4,
# nlme, sommer and pbkrtest, plus the Python suites of remlax.
#
# Submitted by campaign/prepare.sh, which copies the sources to $W and sets W
# and SIF. Each suite writes $W/logs/<name>.log under the name that
# validation/parse_logs.py expects. A failing suite does not stop the others.
# =============================================================================
#SBATCH --job-name=rx_valid_cpu
#SBATCH --partition=cpu-dedicated
#SBATCH --account=dedicated-cpu@cirad-long
#SBATCH --qos=dedicated
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=2-00:00:00

set -uo pipefail
: "${W:?W (work folder) is not set}" "${SIF:?SIF (image) is not set}"
L=$W/logs; mkdir -p "$L" "$W/rlib"

# renv off and an empty user profile: the image keeps its own library. R_LIBS
# puts the remlax under test first. One BLAS thread per process keeps the
# likelihood reproducible to the last digits.
export APPTAINERENV_RENV_CONFIG_AUTOLOADER_ENABLED=FALSE
export APPTAINERENV_R_PROFILE_USER=/dev/null
export APPTAINERENV_R_LIBS="$W/rlib"
export APPTAINERENV_RX_BACKEND=cpu
export APPTAINERENV_JAX_PLATFORMS=cpu
export APPTAINERENV_OMP_NUM_THREADS=1
export APPTAINERENV_OPENBLAS_NUM_THREADS=1
export APPTAINERENV_REMLAX_CACHE="$W/cache"
RUN() { apptainer exec "$SIF" "$@"; }

echo "=== $(hostname) | job ${SLURM_JOB_ID:-local} | $(date) | $SIF ==="
RUN R CMD INSTALL -l "$W/rlib" "$W/remlax" > "$L/install.log" 2>&1 \
  || { echo "remlax installation failed, see $L/install.log"; exit 1; }

RUN Rscript -e 'v <- function(p) as.character(utils::packageVersion(p));
  cat(sprintf("[R] %s | asreml %s | lme4 %s | sommer %s | nlme %s | pbkrtest %s | remlax %s\n",
      R.version.string, v("asreml"), v("lme4"), v("sommer"), v("nlme"), v("pbkrtest"), v("remlax")))' \
  > "$L/00_versions.log" 2>&1
RUN python3 -c "import jax; print('[jax]', jax.__version__)" >> "$L/00_versions.log" 2>&1
cat "$W/REVISIONS" >> "$L/00_versions.log"
cat "$L/00_versions.log"

status=0
step() {                       # step <log> <command...>
  local log=$1; shift
  local t0=$(date +%s)
  "$@" > "$L/$log" 2>&1
  local rc=$?
  printf "%-24s exit %3d  %6d s\n" "$log" "$rc" "$(( $(date +%s) - t0 ))"
  [ $rc -eq 0 ] || status=1
}

cd "$W/remlax"
step 00_pytest.log        RUN python3 -m pytest tests-python -q -p no:cacheprovider
step 01_core.log          RUN python3 tests-python/test_remlax_core.py
step 02_stress_s0.log     RUN python3 tests-python/stress_remlax.py --seed 0
step 02_stress_s50.log    RUN python3 tests-python/stress_remlax.py --seed 50
step 02_stress_s100.log   RUN python3 tests-python/stress_remlax.py --seed 100

cd "$W/validation"
step 03_lme4_sommer.log   RUN Rscript tests/R/test_remlax.R
step 07_nlme.log          RUN Rscript tests/R/test_remlax_nlme.R
step 08_lme4_extra.log    RUN Rscript tests/R/test_remlax_lme4.R
step 09_sommer_extra.log  RUN Rscript tests/R/test_remlax_sommer.R
step 04_asreml.log        RUN Rscript tests/R/test_remlax_asreml.R
step 05_asreml2.log       RUN Rscript tests/R/test_remlax_asreml2.R
step 15_asreml3.log       RUN Rscript tests/R/test_remlax_asreml3.R
step 06b_asreml_blup.log  RUN Rscript tests/R/test_remlax_asreml_blup.R
step 06_ige.log           RUN Rscript tests/R/test_remlax_ige.R
step 16_sep_ar1.log       RUN Rscript tests/R/test_sep_ar1.R
step diag_asreml3_ecarts.log RUN Rscript tests/R/diag_asreml3_ecarts.R

echo "=== end $(date) | overall: $([ $status -eq 0 ] && echo OK || echo FAIL) ==="
exit $status
