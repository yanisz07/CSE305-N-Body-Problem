# CSE305-N-Body-Problem

# N-Body Simulation

**Prerequisite**: Magick++ must be installed.

1) Open a terminal and cd into the project folder (where n_body_basic.cpp lives).

2) To compile run the command:

g++ -std=c++17 -O2 n_body_basic.cpp $(Magick++-config --cxxflags --libs) -o n_body_basic.exe

3) To run the simulation:
./n_body_basic.exe <scenario_name>

Valid scenarios:
  two_body_test
  earth_moon
  jupiter_moons
  solar_system
  milky_way
  large_random_simulation