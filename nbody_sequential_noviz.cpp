// nbody_sequential_noviz.cpp
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include "config.hpp"
#include "body.hpp"

using std::vector;
using std::string;

static constexpr double G_CONST = 6.67430e-11;
static constexpr double PI      = 3.14159265358979323846;

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <scenario>\n";
        return 1;
    }
    string configName = argv[1];
    Config cfg(configName);

    // 1) Build bodies
    vector<Body> bodies;
    if (configName == "benchmark_fixed") {
        const int NX = 100, NY = 100; // N = NX*NY
        const double span = 2e9;
        const double dx = span/(NX-1), dy = span/(NY-1);
        bodies.reserve(NX*NY);
        for (int ix = 0; ix < NX; ++ix)
            for (int iy = 0; iy < NY; ++iy)
                bodies.emplace_back(
                  1e25,
                  Vec{-1e9 + ix*dx, -1e9 + iy*dy},
                  Vec{0,0}
                );
    }
    else if (configName == "earth_moon") {
        bodies = {
            {5.972e24, {0,0},        {0,0}},
            {7.34767309e22, {3.84e8,0}, {0,1022}}
        };
    }
    else if (configName == "jupiter_moons") {
        bodies = {
            {1.898e27, {0,0},       {0,0}},
            {8.93e22,  {4.22e8,0},  {0,17320}},
            {4.8e22,   {6.71e8,0},  {0,13740}},
            {1.48e23,  {1.07e9,0},  {0,10870}},
            {1.08e23,  {1.88e9,0},  {0,8200}}
        };
    }
    else if (configName == "solar_system") {
        bodies = {
            {1.989e30, {0,0},         {0,0}},
            {3.285e23, {5.79e10,0},    {0,47400}},
            {4.867e24, {1.082e11,0},   {0,35020}},
            {5.972e24, {1.496e11,0},   {0,29780}},
            {6.39e23,  {2.279e11,0},   {0,24130}},
            {1.898e27, {7.785e11,0},   {0,13070}}
        };
    }
    else if (configName == "milky_way") {
        std::mt19937_64 rng(123);
        std::uniform_real_distribution<double> rD(0,5e20),
                                              aD(0,2*PI),
                                              mD(1e30,1e32);
        for (int i = 0; i < 100; ++i) {
            double r  = rD(rng);
            double th = aD(rng);
            bodies.emplace_back(
              mD(rng),
              Vec{r*std::cos(th), r*std::sin(th)},
              Vec{0,0}
            );
        }
    }
    else if (configName == "large_random_simulation") {
        std::mt19937_64 rng(456);
        std::uniform_real_distribution<double> pD(-1e10,1e10),
                                              vD(-1e4,1e4),
                                              mD(1e20,1e25);
        for (int i = 0; i < 50; ++i) {
            bodies.emplace_back(
              mD(rng),
              Vec{pD(rng),pD(rng)},
              Vec{vD(rng),vD(rng)}
            );
        }
    }
    else if (configName == "two_body_test") {
        bodies = {
            {1e7, {-1,0}, {0,0}},
            {1e7, { 1,0}, {0,0}}
        };
    }
    else {
        std::cerr << "Unknown scenario: " << configName << "\n";
        return 1;
    }
    int N = bodies.size();

    // 2) Run simulation (no visualization)
    for (int step = 0; step < cfg.steps; ++step) {
        // zero forces
        for (auto &b : bodies) b.force = {0,0};

        // compute gravity
        for (int i = 0; i < N; ++i) {
            for (int j = i+1; j < N; ++j) {
                double dx = bodies[j].pos.x - bodies[i].pos.x;
                double dy = bodies[j].pos.y - bodies[i].pos.y;
                double d2 = dx*dx + dy*dy + 1e-12;
                double inv = 1.0/std::sqrt(d2);
                double F = G_CONST * bodies[i].m * bodies[j].m * inv*inv;
                double fx = F * dx * inv;
                double fy = F * dy * inv;
                bodies[i].force.x += fx; bodies[i].force.y += fy;
                bodies[j].force.x -= fx; bodies[j].force.y -= fy;
            }
        }

        // update
        for (auto &b : bodies) {
            b.vel.x += (b.force.x/b.m) * cfg.timestep;
            b.vel.y += (b.force.y/b.m) * cfg.timestep;
            b.pos.x += b.vel.x * cfg.timestep;
            b.pos.y += b.vel.y * cfg.timestep;
        }
    }

    return 0;
}
