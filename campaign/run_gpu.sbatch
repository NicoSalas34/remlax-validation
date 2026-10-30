#!/bin/bash
# =============================================================================
# Validation campaign of remlax, GPU part: CPU/GPU parity on a battery of
# serialised designs. R writes the bundles; Python fits each one on both
# backends of the same node. Only the device changes, so any gap beyond the
# reduction noise is a defect.
#
# Submitted by campaign/prepare.sh after run_cpu.sbatch has installed remlax.
# The partition and the card are set at submission (see prepare.sh).
# =============================================================================
#SBATCH --job-name=rx_valid_gpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00

set -uo pipefail
: "${W:?W (work folder) is not set}" "${SIF:?SIF (image) is not set}"
L=$W/logs; mkdir -p "$L"
export APPTAINERENV_RENV_CONFIG_AUTOLOADER_ENABLED=FALSE
export APPTAINERENV_R_PROFILE_USER=/dev/null
export APPTAINERENV_R_LIBS="$W/rlib"
export APPTAINERENV_XLA_PYTHON_CLIENT_PREALLOCATE=false
export APPTAINERENV_REMLAX_CACHE="$W/cache"
RUN() { apptainer exec --nv "$SIF" "$@"; }

echo "=== $(hostname) | job ${SLURM_JOB_ID:-local} | $(date) | $SIF ==="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
# The CPU job installs remlax into $W/rlib when it starts.
for i in $(seq 1 120); do [ -f "$W/rlib/remlax/DESCRIPTION" ] && break; sleep 30; done
[ -f "$W/rlib/remlax/DESCRIPTION" ] || { echo "remlax not installed in $W/rlib"; exit 1; }
cd "$W/validation"
RUN Rscript tests/R/export_bundles.R --out="$W/bundles" > "$L/20_export_bundles.log" 2>&1 \
  || { echo "bundle export failed, see $L/20_export_bundles.log"; exit 1; }
cd "$W/remlax"
# 34 models x 2 backends: compilation dominates on a small card. Expect tens
# of minutes, not a hang.
RUN python3 tests-python/parite_gpu.py "$W/bundles" > "$L/21_parite.log" 2>&1
rc=$?
tail -5 "$L/21_parite.log"
exit $rc
