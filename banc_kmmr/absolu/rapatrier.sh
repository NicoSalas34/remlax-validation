#!/bin/bash
# Fetches the CSV files of the absolute benchmark from genotoul and meso into absolu/resultats
D=$(dirname $(readlink -f $0))/resultats
scp -q "genotoul:/work/user/nsalas/kmmr_banc/res/abs_cpu*.csv" "genotoul:/work/user/nsalas/kmmr_banc/res/abs_gpu*.csv" $D/ 2>/dev/null
scp -q "cluster:/scratch/users/salasn/kmmr_banc/res/abs_meso_x13*.csv" $D/ 2>/dev/null
cp $(dirname $D)/axes.csv $D/
ls $D
