#!/bin/bash
# Relance les cas du banc sans CSV final, ou dont le CSV contient un echec de
# licence asreml. executer.R reprend chaque cas la ou il s'est arrete : un
# ajustement qui avait tue la tache (memoire) est note MEMOIRE et saute.
#
#   bash relancer_manquants.sh            # affiche les indices
#   bash relancer_manquants.sh --soumettre
W=${W:-/scratch/users/salasn/banc_moteurs_20261022}
FAMILLES=(iid grm us ped ige champ resid exp fa)
TAILLES=(8000 6000 4000 3000 2000 1500 1000 500)
NREP=5
idx=()
for ((i = 0; i < ${#FAMILLES[@]} * ${#TAILLES[@]} * NREP; i++)); do
  j=$i; r=$(( j % NREP + 1 )); j=$(( j / NREP ))
  f=${FAMILLES[$(( j % ${#FAMILLES[@]} ))]}; j=$(( j / ${#FAMILLES[@]} ))
  n=${TAILLES[$j]}
  fcsv="$W/res/ajustements_${f}_${n}_${r}.csv"
  if [ ! -f "$fcsv" ] || grep -q "License error\|licenses in use\|workspace" "$fcsv"; then idx+=("$i"); fi
done
liste=$(IFS=,; echo "${idx[*]}")
echo "${#idx[@]} cas a relancer : ${liste:-aucun}"
if [ "${1:-}" = "--soumettre" ] && [ ${#idx[@]} -gt 0 ]; then
  sbatch --array="${liste}%6" "$W/remlax-validation/banc_moteurs/banc_mig.sbatch"
fi
