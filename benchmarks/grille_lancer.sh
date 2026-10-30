#!/usr/bin/env bash
# Grid driver: ONE cell per process, from the smallest to the largest, within
# a time budget. The large cells are the most informative but the most
# expensive. They are reached only if the budget allows, instead of losing the
# whole sweep to an exceeded walltime.
set -u
RS="${RS:-Rscript}"; REPO="${REPO:-.}"; OUT="${OUT:-grid.csv}"
BK="${BK:-cpu}"; TAG="${TAG:-}"; BUDGET="${BUDGET:-18000}"
CELLS="${CELLS:-sparse-dense dense-dense}"
NS="${NS:-2000 4000 8000 16000}"; QS="${QS:-500 2000}"; TS="${TS:-1 2 4 8}"
CAP="${CAP:-1800}"   # time cap per cell, in seconds
start=$(date +%s)
for t in $TS; do for n in $NS; do for q in $QS; do
  [ "$q" -ge "$n" ] && continue          # q >= n: degenerate case, not identified
  for cl in $CELLS; do
    elapsed=$(( $(date +%s) - start ))
    if [ "$elapsed" -gt "$BUDGET" ]; then
      echo "[budget] $elapsed s elapsed, clean stop before $cl n=$n q=$q t=$t"
      exit 0
    fi
    case_="${cl%%-*}"; eng="${cl##*-}"
    t0=$(date +%s)
    timeout "$CAP" $RS "$REPO/benchmarks/grille.R" --case "$case_" --engine "$eng" \
      --backend "$BK" --n "$n" --q "$q" --t "$t" --out "$OUT" --tag "$TAG" \
      --repo "$REPO" 2>&1 | grep -viE "codetools|adfun|remplie a 100|100% full"
    rc=$?
    [ "$rc" -eq 124 ] && echo "[cap] $cl n=$n q=$q t=$t exceeded $CAP s, skipped"
    echo "    ($(( $(date +%s) - t0 )) s, total $(( $(date +%s) - start )) s)"
  done
done; done; done
echo "[end] $(( $(date +%s) - start )) s in total"
