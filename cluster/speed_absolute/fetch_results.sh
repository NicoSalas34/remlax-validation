#!/bin/bash
# Fetches the CSV files of the absolute benchmark from genotoul and meso into
# the inputs/ folder of a new run of speed/absolute/latest/.
#   cluster/speed_absolute/fetch_results.sh speed/absolute/latest/<date>_remlax-<version>/inputs
D=${1:?usage: fetch_results.sh <run folder>/inputs}
mkdir -p "$D"
scp -q "genotoul:/work/user/nsalas/kmmr_banc/res/abs_cpu*.csv" "genotoul:/work/user/nsalas/kmmr_banc/res/abs_gpu*.csv" "$D"/ 2>/dev/null
scp -q "cluster:/scratch/users/salasn/kmmr_banc/res/abs_meso_x13*.csv" "$D"/ 2>/dev/null
cp "$(dirname "$(readlink -f "$0")")/axes.csv" "$D"/
ls "$D"
