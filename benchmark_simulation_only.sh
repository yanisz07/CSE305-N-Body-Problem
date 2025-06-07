#!/usr/bin/env bash

# benchmark_simulation_only.sh
# Runs only the integration (no GIF generation) for sequential and parallel versions,
# measures wall-clock time (averaged over multiple runs), and prints a table of results.

# Adjust these if you want different thread counts or number of iterations per measurement:
THREAD_COUNTS=(1 2 4 8)
ITER=3

# List of valid scenarios:
SCENARIOS=(
  "two_body_test"
  "earth_moon"
  "jupiter_moons"
  "solar_system"
  "milky_way"
  "large_random_simulation"
)

# Names of the “simulation-only” executables (must be compiled without any GIF code):
SEQ_SIM_BIN="./nbody_sequential_sim.exe"
PAR_SIM_BIN="./nbody_parallel_sim.exe"

# Check that both executables exist:
if [[ ! -x "$SEQ_SIM_BIN" ]]; then
  echo "Error: $SEQ_SIM_BIN not found or not executable."
  exit 1
fi
if [[ ! -x "$PAR_SIM_BIN" ]]; then
  echo "Error: $PAR_SIM_BIN not found or not executable."
  exit 1
fi

# Function to measure average wall-clock time (seconds) over $ITER runs:
measure_avg_time() {
  local cmd=("$@")
  local total=0

  for ((i=1; i<=ITER; i++)); do
    start=$(date +%s.%N)
    "${cmd[@]}" >/dev/null
    end=$(date +%s.%N)
    dt=$(echo "$end - $start" | bc -l)
    total=$(echo "$total + $dt" | bc -l)
  done

  # Compute average = total / ITER
  avg=$(echo "scale=6; $total / $ITER" | bc -l)
  printf "%.6f" "$avg"
}

# Print header for the full benchmark:
echo
echo "=================================================================================="
echo "    N-Body Simulation-Only Benchmark (averaged over $ITER runs each)             "
echo "    Sequential sim: $SEQ_SIM_BIN                                                 "
echo "    Parallel sim  : $PAR_SIM_BIN                                                 "
echo "=================================================================================="
printf "%-20s | %10s" "Scenario" "SeqTime(s)"
for th in "${THREAD_COUNTS[@]}"; do
  printf " | %10s  %10s" "Par($th)" "Speedup"
done
echo
sep="----------------------|------------"
for _ in "${THREAD_COUNTS[@]}"; do
  sep+="|----------------------"
done
echo "$sep"

# Loop over each scenario:
for scenario in "${SCENARIOS[@]}"; do
  # Measure sequential simulation-only time:
  seq_time=$(measure_avg_time "$SEQ_SIM_BIN" "$scenario")

  # Prepare arrays to store parallel times and speedups:
  declare -a par_times
  declare -a speedups

  # Measure parallel simulation-only times for each thread count:
  for th in "${THREAD_COUNTS[@]}"; do
    par_time=$(measure_avg_time "$PAR_SIM_BIN" "$scenario" "$th")
    par_times+=( "$par_time" )

    # speedup = seq_time / par_time
    speedup=$(echo "scale=3; $seq_time / $par_time" | bc -l)
    speedups+=( "$speedup" )
  done

  # Print one line of results for this scenario:
  printf "%-20s | %10.6f" "$scenario" "$seq_time"
  for idx in "${!THREAD_COUNTS[@]}"; do
    printf " | %10.6f  %10.3f" "${par_times[$idx]}" "${speedups[$idx]}"
  done
  echo
done

echo "----------------------------------------------------------------------------------"
echo "Note: Times are in seconds. Speedup = (Sequential sim time) / (Parallel sim time)."
echo "=================================================================================="
