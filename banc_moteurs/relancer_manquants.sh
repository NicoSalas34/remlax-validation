#!/bin/bash
# Resubmits the benchmark cases without a final CSV, or whose CSV contains an
# asreml licence failure. executer.R resumes each case where it stopped: a fit
# that had killed the task (memory) is recorded as MEMORY and skipped.
#
#   bash relancer_manquants.sh            # prints the indices
#   bash relancer_manquants.sh --submit
W=${W:-/scratch/users/salasn/banc_moteurs_20261022}
FAMILIES=(iid grm us ped ige field resid exp fa)
SIZES=(8000 6000 4000 3000 2000 1500 1000 500)
NREP=5
idx=()
for ((i = 0; i < ${#FAMILIES[@]} * ${#SIZES[@]} * NREP; i++)); do
  j=$i; r=$(( j % NREP + 1 )); j=$(( j / NREP ))
  f=${FAMILIES[$(( j % ${#FAMILIES[@]} ))]}; j=$(( j / ${#FAMILIES[@]} ))
  n=${SIZES[$j]}
  fcsv="$W/res/fits_${f}_${n}_${r}.csv"
  if [ ! -f "$fcsv" ] || grep -q "License error\|licenses in use\|workspace" "$fcsv"; then idx+=("$i"); fi
done
list=$(IFS=,; echo "${idx[*]}")
echo "${#idx[@]} cases to resubmit: ${list:-none}"
if [ "${1:-}" = "--submit" ] && [ ${#idx[@]} -gt 0 ]; then
  sbatch --array="${list}%6" "$W/remlax-validation/banc_moteurs/banc_mig.sbatch"
fi
