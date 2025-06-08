#!/usr/bin/env bash
set -euo pipefail

SEQ="./nbody_sequential_noviz.exe"
PAR="./nbody_parallel_noviz.exe"
SCEN="benchmark_fixed"
threads=(1 2 4 8)
runs=5

# CSV header
printf "scenario,seq_s"
for T in "${threads[@]}"; do
  printf ",par${T}_s,speedup${T}"
done
echo

# Sequential timing
sum_seq=0
for i in $(seq 1 $runs); do
  t0=$(date +%s%N)
  $SEQ "$SCEN"
  t1=$(date +%s%N)
  dt=$(awk "BEGIN{print ($t1 - $t0)/1e9}")
  sum_seq=$(awk "BEGIN{print $sum_seq + $dt}")
done
avg_seq=$(awk "BEGIN{print $sum_seq / $runs}")
line="$SCEN,$avg_seq"

# Parallel timing
for T in "${threads[@]}"; do
  sum_par=0
  for i in $(seq 1 $runs); do
    t0=$(date +%s%N)
    $PAR "$SCEN" $T
    t1=$(date +%s%N)
    dt=$(awk "BEGIN{print ($t1 - $t0)/1e9}")
    sum_par=$(awk "BEGIN{print $sum_par + $dt}")
  done
  avg_par=$(awk "BEGIN{print $sum_par / $runs}")
  speedup=$(awk "BEGIN{print $avg_seq / $avg_par}")
  line+=",${avg_par},${speedup}"
done

echo "$line"
