#!/usr/bin/env bash
set -euo pipefail

SCENARIO=benchmark_fixed
SEQ_EXE=./nbody_sequential.exe
PAR_EXE=./nbody_parallel.exe
THREADS=(1 2 4 8)

declare -A time_s

# run sequential
echo "Running sequential..."
t0=$(date +%s.%N)
"$SEQ_EXE" "$SCENARIO"
t1=$(date +%s.%N)
time_s[seq]=$(awk -v a=$t0 -v b=$t1 'BEGIN{printf "%.4f", b-a}')

# run parallel for each thread count
for T in "${THREADS[@]}"; do
  echo "Running parallel (T=$T)..."
  t0=$(date +%s.%N)
  "$PAR_EXE" "$SCENARIO" "$T"
  t1=$(date +%s.%N)
  time_s[$T]=$(awk -v a=$t0 -v b=$t1 'BEGIN{printf "%.4f", b-a}')
done

# output CSV
echo
echo "threads,time_s,speedup"
echo "seq,${time_s[seq]},1.0000"
for T in "${THREADS[@]}"; do
  speedup=$(awk -v s=${time_s[seq]} -v p=${time_s[$T]} 'BEGIN{printf "%.4f", s/p}')
  echo "$T,${time_s[$T]},${speedup}"
done
