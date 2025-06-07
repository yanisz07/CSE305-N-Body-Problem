#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# ─── User‐configurable paths ────────────────────────────────────────────────
SRC_SEQ="barnes_hut_sequential_NOGIF.cpp"
SRC_PAR="barnes_hut_parallel_NOGIF.cpp"
BIN_SEQ="barnes_hut_seq"
BIN_PAR="barnes_hut_par"
SCEN="benchmark_fixed"
runs=3
threads=(1 2 4 8)
NO_GIF_FLAG="-DNO_GIF"
# ─────────────────────────────────────────────────────────────────────────────

# Build (θ defaults to 0.5 in code)
echo "BUILD: compiling binaries..." >&2
g++ -std=c++17 -O2 $NO_GIF_FLAG "$SRC_SEQ" -o "$BIN_SEQ" \
    $(Magick++-config --cppflags --cxxflags --ldflags --libs) >&2
g++ -std=c++17 -O2 $NO_GIF_FLAG "$SRC_PAR" -o "$BIN_PAR" \
    $(Magick++-config --cppflags --cxxflags --ldflags --libs) >&2

# CSV header (stdout)
{
  printf "scenario,seq_s"
  for T in "${threads[@]}"; do
    printf ",par${T}_s,speedup${T}"
  done
  echo
} 

# Sequential baseline
echo "SEQ: ${SCEN} (${runs} runs)..." >&2
sum_seq=0
for i in $(seq 1 $runs); do
  t0=$(date +%s%N)
  ./"$BIN_SEQ" "$SCEN" > /dev/null
  t1=$(date +%s%N)
  dt=$(awk "BEGIN{print ($t1 - $t0)/1e9}")
  sum_seq=$(awk "BEGIN{print $sum_seq + $dt}")
  echo "  run $i: ${dt}s" >&2
done
avg_seq=$(awk "BEGIN{print $sum_seq / $runs}")
echo "SEQ avg: ${avg_seq}s" >&2

# Start building CSV row
csv_line="${SCEN},${avg_seq}"

# Parallel sweeps
for T in "${threads[@]}"; do
  echo "PAR: ${SCEN}, T=${T} (${runs} runs)..." >&2
  sum_par=0
  for i in $(seq 1 $runs); do
    t0=$(date +%s%N)
    ./"$BIN_PAR" "$SCEN" $T > /dev/null
    t1=$(date +%s%N)
    dt=$(awk "BEGIN{print ($t1 - $t0)/1e9}")
    sum_par=$(awk "BEGIN{print $sum_par + $dt}")
    echo "  run $i: ${dt}s" >&2
  done
  avg_par=$(awk "BEGIN{print $sum_par / $runs}")
  speedup=$(awk "BEGIN{print $avg_seq / $avg_par}")
  echo "PAR T=${T} avg: ${avg_par}s → speedup ${speedup}" >&2

  csv_line+=",${avg_par},${speedup}"
done

# Emit the CSV line to stdout
echo "$csv_line"
