// config.hpp
#pragma once

#include <string>
#include <iostream>
#include <cstdlib>

struct Config {
    double distanceScale;  // maps physical [m] to normalized [0→1]
    double timestep;       // seconds per step
    int    width, height;  
    int    steps;          // number of steps

    Config(const std::string& name) {
        width  = 800;
        height = 800;

        if (name == "earth_moon") {
            distanceScale = 0.5  / 3.84e8;   // ±3.84×10⁸ m → ±0.5
            timestep      = 3600;           // 1 h
            steps         = 200;
        }
        else if (name == "jupiter_moons") {
            distanceScale = 0.5  / 4.22e8;
            timestep      = 3600;
            steps         = 200;
        }
        else if (name == "solar_system") {
            distanceScale = 0.5  / 1.496e11;
            timestep      = 86400;          // 1 day
            steps         = 365;
        }
        else if (name == "milky_way") {
            distanceScale = 0.5  / 5e20;
            timestep      = 5e6 * 3.1536e7; // 5 Myr
            steps         = 100;
        }
        else if (name == "large_random_simulation") {
            distanceScale = 0.5  / 1e10;
            timestep      = 86400;
            steps         = 365;
        }
        else if (name == "two_body_test") {
            distanceScale = 0.5  / 2.0;
            timestep      = 0.05;
            steps         = 200;
        }
        else if (name == "benchmark_fixed") {
            distanceScale = 0.5  / 1e9;
            timestep      = 3600;
            steps         = 100;
        }
        else {
            std::cerr << "Invalid config: " << name << "\n";
            std::exit(1);
        }
    }
};
