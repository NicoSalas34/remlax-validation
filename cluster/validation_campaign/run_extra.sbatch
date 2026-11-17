#!/bin/bash
# =============================================================================
# Validation campaign of remlax, extra part: the testthat suite of the package
# and the simulation checks (verif_simulations, dense REML arbiter and
# simulated truth).
#
# Submitted after run_cpu.sbatch, which installs remlax into $W/rlib. W and SIF
# are set at submission, as for the other two jobs.
# =============================================================================
#SBATCH --job-name=rx_valid_extra
#SBATCH --partition=cpu-dedicated
#SBATCH --account=dedicated-cpu@cirad-long
#SBATCH --qos=dedicated
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --time=12:00:00

set -uo pipefail
: "${W:?W (work folder) is not set}" "${SIF:?SIF (image) is not set}"
L=$W/logs; mkdir -p "$L"
export APPTAINERENV_RENV_CONFIG_AUTOLOADER_ENABLED=FALSE
export APPTAINERENV_R_PROFILE_USER=/dev/null
export APPTAINERENV_R_LIBS="$W/rlib"
export APPTAINERENV_RX_BACKEND=cpu
export APPTAINERENV_JAX_PLATFORMS=cpu
export APPTAINERENV_REMLAX_CACHE="$W/cache_extra"
RUN() { apptainer exec "$SIF" "$@"; }

echo "=== $(hostname) | job ${SLURM_JOB_ID:-local} | $(date) | $SIF ==="
for i in $(seq 1 120); do [ -f "$W/rlib/remlax/DESCRIPTION" ] && break; sleep 30; done
[ -f "$W/rlib/remlax/DESCRIPTION" ] || { echo "remlax not installed in $W/rlib"; exit 1; }

status=0
step() {                       # step <log> <command...>
  local log=$1; shift
  local t0=$(date +%s)
  "$@" > "$L/$log" 2>&1
  local rc=$?
  printf "%-24s exit %3d  %6d s\n" "$log" "$rc" "$(( $(date +%s) - t0 ))"
  [ $rc -eq 0 ] || status=1
}

cd "$W/remlax/tests"
step 30_testthat.log RUN env OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 Rscript -e '
  library(testthat); library(remlax)
  cat("remlax", as.character(packageVersion("remlax")), system.file(package = "remlax"), "\n")
  r <- as.data.frame(test_dir("testthat", package = "remlax", load_package = "installed",
                              reporter = "summary", stop_on_failure = FALSE))
  print(r[, c("file", "test", "nb", "failed", "skipped", "error")], right = FALSE)
  cat(sprintf("TOTAL expectations %d | failed %d | errors %d | skipped tests %d\n",
              sum(r$nb), sum(r$failed), sum(r$error), sum(r$skipped)))
  quit(status = as.integer(sum(r$failed) + sum(r$error) > 0))'

cd "$W/validation"
step 31_verif_simulations.log RUN Rscript correctness/simulations/scripts/verif_simulations.R \
  --cores=10 --out="$W/verif_simulations"

echo "=== end $(date) | overall: $([ $status -eq 0 ] && echo OK || echo FAIL) ==="
exit $status
