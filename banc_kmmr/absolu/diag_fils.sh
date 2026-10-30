#!/bin/bash
# Diagnostic: the same remlax fit three times in a row, recording threads, CPU
# time and load. Arguments: root K, image, R_LIBS, output.
K=$1; SIF=$2; RL=$3; OUT=$4; CASE=${CASE:-sim_q1000_t8_r1_sdiag_x0_o100}
E="--env RENV_CONFIG_AUTOLOADER_ENABLED=FALSE --env R_PROFILE_USER=/dev/null --env R_LIBS=$RL --env OMP_NUM_THREADS=8 --env OPENBLAS_NUM_THREADS=8 --env JAX_PLATFORMS=cpu --env TMPDIR=$K/tmp --env REMLAX_JAX_CACHE=$K/tmp/jaxdiag_$SLURM_JOB_ID"
echo "node $(hostname) | visible cores $(nproc) | cpuset $(cat /sys/fs/cgroup/cpuset.cpus.effective 2>/dev/null || grep Cpus_allowed_list /proc/self/status) | $(lscpu | grep -E 'Thread|Socket|NUMA node\(s\)' | tr -s ' ' | tr '\n' ';')"
for i in 1 2 3; do
  /usr/bin/time -f "rep$i real %e s | user %U s | sys %S s | CPU %P" apptainer exec $E $SIF Rscript $K/banc/R/executer.R --case $CASE --tool remlax --route mme --gpu 0 --kmmr $K/src/kmmr --sim $K/absolu --out $OUT --label rep$i > $OUT.rep$i.log 2>&1 &
  P=$!
  while kill -0 $P 2>/dev/null; do
    echo "  rep$i $(date +%T) load $(cut -d' ' -f1 /proc/loadavg) | python threads $(ps -C python3 -o nlwp= | paste -sd+ | bc 2>/dev/null) | python CPU $(ps -C python3 -o %cpu= | paste -sd+ | bc 2>/dev/null)% | servers $(pgrep -fc 'serveur.py|server.py')"
    sleep 5
  done
  wait $P; grep "rep$i real" $OUT.rep$i.log
  # AI trace and compilation lines: "[compilation]" in older remlax logs, "[compile]" in newer ones
  grep -hE "aiv|\[compil(ation|e)\]" $OUT.rep$i.log | tail -3
done
