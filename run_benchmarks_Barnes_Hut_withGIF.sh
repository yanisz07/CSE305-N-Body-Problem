#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# ─── Adjust these to your project ─────────────────────────────────────────────
SRC_SEQ="Barnes_Hut_Sequential.cpp"       # your sequential source
SRC_PAR="barnes_hut_parallel.cpp"    # your parallel source
BIN_SEQ="barnes_hut_seq"             # sequential binary name
BIN_PAR="barnes_hut_par"             # parallel binary name
SCEN="benchmark_fixed"               # config to benchmark
GIF_DIR="benchmark_gifs/$SCEN"       # where to stash GIFs (if any)
runs=3                               # repetitions per setting
threads=(1 2 4 8)                    # thread‐counts to sweep over
# ───────────────────────────────────────────────────────────────────────────────

mkdir -p "$GIF_DIR"

# CSV header
printf "scenario,seq_s"
for T in "${threads[@]}"; do
  printf ",par${T}_s,speedup${T}"
done
echo

# 1) Build both binaries (θ defaults to 0.5 in your code)
echo "[BUILD] compiling sequential → $BIN_SEQ, parallel → $BIN_PAR"
g++ -std=c++17 -O2 "$SRC_SEQ" -o "$BIN_SEQ" \
    $(Magick++-config --cppflags --cxxflags --ldflags --libs)
g++ -std=c++17 -O2 "$SRC_PAR" -o "$BIN_PAR" \
    $(Magick++-config --cppflags --cxxflags --ldflags --libs)

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
