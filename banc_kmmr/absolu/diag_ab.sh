#!/bin/bash
# Essai A/B alterne sur un meme noeud : sans srun (A) puis avec srun --cpu-bind=cores (B).
K=/scratch/users/salasn/kmmr_banc; SIF=/scratch/users/salasn/ige_reml_v2.sif
CAS=${CAS:-sim_q1000_t8_r1_sdiag_x0_o100}
E="--env RENV_CONFIG_AUTOLOADER_ENABLED=FALSE --env R_PROFILE_USER=/dev/null --env R_LIBS=$K/rlib_remlax_v11:$K/rlib --env OMP_NUM_THREADS=8 --env OPENBLAS_NUM_THREADS=8 --env JAX_PLATFORMS=cpu --env TMPDIR=$K/tmp --env REMLAX_JAX_CACHE=$K/tmp/jaxab_$SLURM_JOB_ID"
echo "noeud $(hostname) | affinite du lot $(python3 -c 'import os;print(sorted(os.sched_getaffinity(0)))')"
for essai in A1 B1 A2 B2; do
  case $essai in A*) PRE="";; B*) PRE="srun --cpu-bind=cores --cpus-per-task=8";; esac
  t0=$(date +%s.%N)
  $PRE apptainer exec $E $SIF Rscript $K/banc/R/executer.R --cas $CAS --outil remlax --voie mme --gpu 0 --kmmr $K/src/kmmr --sim $K/absolu --sortie $K/res/diag_ab.csv --etiquette $essai > $K/res/diag_ab.$essai.log 2>&1 &
  P=$!
  pic=0; aff=""
  while kill -0 $P 2>/dev/null; do
    S=$(pgrep -u $USER -f "serveur.py" | head -1)
    if [ -n "$S" ]; then
      f=$(awk '/Threads/{print $2}' /proc/$S/status 2>/dev/null); [ -n "$f" ] && [ "$f" -gt "$pic" ] && pic=$f
      aff=$(awk '/Cpus_allowed_list/{print $2}' /proc/$S/status 2>/dev/null)
      cs=$(awk '/nonvoluntary_ctxt_switches/{print $2}' /proc/$S/status 2>/dev/null)
    fi
    sleep 2
  done
  t1=$(date +%s.%N)
  printf "%s | %6.1f s | fils du serveur (pic) %s | affinite du serveur %s | changements forces %s | charge %s\n" $essai $(echo "$t1 - $t0" | bc) "$pic" "$aff" "$cs" "$(cut -d' ' -f1 /proc/loadavg)"
  grep -h "aiv     5" $K/res/diag_ab.$essai.log
done
