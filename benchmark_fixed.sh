#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

SEQ="./nbody_sequential.exe"
PAR="./nbody_parallel.exe"
SCEN="benchmark_fixed"
GIF_DIR="benchmark_gifs/$SCEN"
threads=(1 2 4 8)
runs=3

mkdir -p "$GIF_DIR"

# CSV header
printf "scenario,seq_s"
for T in "${threads[@]}"; do
  printf ",par${T}_s,speedup${T}"
done
echo

# Sequential
echo "[SEQ] ${SCEN} sequential (${runs}×)..."
sum_seq=0
for i in $(seq 1 $runs); do
  t0=$(date +%s%N)
  $SEQ "$SCEN"   > /dev/null
  t1=$(date +%s%N)
  dt=$(awk "BEGIN{print ($t1 - $t0)/1e9}")
  sum_seq=$(awk "BEGIN{print $sum_seq + $dt}")
  for f in anime_${SCEN}_*.gif; do mv "$f" "$GIF_DIR/"; done
  echo "      run $i: ${dt}s"
done
avg_seq=$(awk "BEGIN{print $sum_seq / $runs}")
echo "[SEQ] avg: ${avg_seq}s"
line="${SCEN},${avg_seq}"

# Parallel
for T in "${threads[@]}"; do
  echo "[PAR] ${SCEN}, T=${T} (${runs}×)..."
  sum_par=0
  for i in $(seq 1 $runs); do
    t0=$(date +%s%N)
    $PAR "$SCEN" $T > /dev/null
    t1=$(date +%s%N)
    dt=$(awk "BEGIN{print ($t1 - $t0)/1e9}")
    sum_par=$(awk "BEGIN{print $sum_par + $dt}")
    for f in gif_parallel_${SCEN}_*.gif; do mv "$f" "$GIF_DIR/"; done
    echo "      run $i: ${dt}s"
  done
  avg_par=$(awk "BEGIN{print $sum_par / $runs}")
  speedup=$(awk "BEGIN{print $avg_seq / $avg_par}")
  echo "[PAR] T=${T} avg: ${avg_par}s → speedup ${speedup}"
  line+=",${avg_par},${speedup}"
done

# Emit CSV row
echo "$line"
