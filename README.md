CSE305-N-Body-Problem
N-Body Simulation
Prerequisite: Magick++ must be installed.

1) Open a terminal and cd into the project folder (where body.hpp, config.hpp, nbody_sequential.cpp, and nbody_parallel.cpp live).

2) To compile, run the commands:

g++ -std=c++17 -O2 nbody_sequential.cpp $(Magick++-config --cxxflags --libs) -o nbody_sequential.exe

g++ -std=c++17 -O2 nbody_parallel.cpp   $(Magick++-config --cxxflags --libs) -o nbody_parallel.exe

g++ -std=c++17 -O2 Barnes_Hut_Sequential.cpp $(Magick++-config --cxxflags --libs) -o barnes_hut_sequential

g++ -std=c++17 -O2 Barnes_Hut_Parallel.cpp $(Magick++-config --cxxflags --libs) -o barnes_hut_parallel



3) To run the sequential simulation:

NAIVE/STRAIGHTFORWARD ALGO:
./nbody_sequential.exe <scenario_name>

BARNES HUT:
./barnes_hut_sequential.exe <scenario_name>

(Note: Ignore the NOGIF versions. They were used only for benchmarks)

4) To run the parallel simulation:

NAIVE/STRAIGHTFORWARD ALGO:
./nbody_parallel.exe <scenario_name> [num_threads]

BARNES HUT:
./barnes_hut_parallel.exe <scenario_name> [num_threads]

Note that if you omit [num_threads], the program uses all available hardware threads.

------------------------------------------------------------------------------------------

Valid <scenario_name> options:

two_body_test
earth_moon
jupiter_moons
solar_system
milky_way
large_random_simulation
benchmark_fixed


After running, the program writes an animated GIF named:

animation_<sequential/parallel>_YYYYMMDD_HHMMSS.gif


