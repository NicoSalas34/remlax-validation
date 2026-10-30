#!/bin/bash
# =============================================================================
# Copies remlax and remlax-validation at their current commits to the cluster
# and submits the two jobs of the campaign.
#
#   campaign/prepare.sh [--remlax=<repo>] [--host=cluster] [--dry-run]
#
# Environment:
#   SIF   image on the cluster (default /scratch/users/salasn/ige_reml_v3.sif)
#   BASE  parent of the work folder (default /scratch/users/salasn)
#   GPU_OPTS  sbatch options of the GPU job (default: MIG slice of an A100)
#
# The work folder is $BASE/validation_<date>_<remlax commit>. Fetch the logs
# afterwards with
#   rsync -a <host>:<work folder>/logs/ validation/logs_<date>/
# then run validation/parse_logs.py on that folder.
# =============================================================================
set -euo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd)
REMLAX=$HOME/Documents/Git/remlax
HOST=cluster
DRY=0
for a in "$@"; do
  case $a in
    --remlax=*) REMLAX=${a#*=} ;;
    --host=*) HOST=${a#*=} ;;
    --dry-run) DRY=1 ;;
    *) echo "unknown argument $a"; exit 2 ;;
  esac
done
SIF=${SIF:-/scratch/users/salasn/ige_reml_v3.sif}
BASE=${BASE:-/scratch/users/salasn}
GPU_OPTS=${GPU_OPTS:---partition=gpu-cirad-dedicated --account=dedicated-gpu@cirad-long --qos=dedicated --gres=gpu:nvidia_a100_80gb_pcie_1g.10gb:1}

for r in "$REMLAX" "$HERE"; do
  if [ -n "$(git -C "$r" status --porcelain --untracked-files=no)" ]; then
    echo "$r has uncommitted changes: commit them first, the campaign records commits"; exit 1
  fi
done
RX=$(git -C "$REMLAX" rev-parse --short HEAD)
VA=$(git -C "$HERE" rev-parse --short HEAD)
W=$BASE/validation_$(date +%Y-%m-%d)_$RX
TMP=$(mktemp -d)
git -C "$REMLAX" archive --format=tar.gz --prefix=remlax/ HEAD > "$TMP/remlax.tgz"
git -C "$HERE" archive --format=tar.gz --prefix=validation/ HEAD \
  -- R tests campaign validation > "$TMP/validation.tgz"
printf "remlax commit %s\nremlax-validation commit %s\nimage %s\n" "$RX" "$VA" "$SIF" > "$TMP/REVISIONS"
echo "work folder: $HOST:$W"
[ $DRY -eq 1 ] && { ls -la "$TMP"; exit 0; }

ssh "$HOST" "mkdir -p $W && test -f $SIF" || { echo "missing $SIF on $HOST"; exit 1; }
scp -q "$TMP"/remlax.tgz "$TMP"/validation.tgz "$TMP"/REVISIONS "$HOST:$W/"
ssh "$HOST" "cd $W && tar xzf remlax.tgz && tar xzf validation.tgz && touch validation/.here && \
  c=\$(sbatch --parsable --output=$W/cpu_%j.out --export=ALL,W=$W,SIF=$SIF validation/campaign/run_cpu.sbatch) && \
  g=\$(sbatch --parsable $GPU_OPTS --dependency=after:\$c --output=$W/gpu_%j.out --export=ALL,W=$W,SIF=$SIF validation/campaign/run_gpu.sbatch) && \
  echo \"submitted: cpu \$c, gpu \$g\""
rm -rf "$TMP"
