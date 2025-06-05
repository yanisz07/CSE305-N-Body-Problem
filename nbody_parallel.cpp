// nbody_parallel.cpp

#include "body.hpp"
#include "config.hpp"

#include <Magick++.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <chrono>
#include <ctime>
#include <random>
#include <thread>     
#include <algorithm>

using namespace Magick;
using std::vector;
using std::string;

static constexpr double G = 6.67430e-11;

int main(int argc, char** argv) {
    if (argc < 2 || argc > 3) {
        std::cerr << "Usage: " << argv[0]
                  << " <config_name> [num_threads]\n"
                  << "Valid configs: two_body_test, earth_moon, jupiter_moons,\n"
                  << "               solar_system, milky_way, large_random_simulation\n";
        return 1;
    }
    string configName = argv[1];
    int T = (argc == 3)
                ? std::stoi(argv[2])
                : std::thread::hardware_concurrency();
    if (T < 1) T = 1;

    Config cfg(configName);
    InitializeMagick(*argv);

    vector<Body> bodies;
    if (configName == "earth_moon") {
        bodies = {
            {5.972e24,      {0.0, 0.0},       {0.0, 0.0}},  
            {7.34767309e22, {3.84e8, 0.0},    {0.0, 1022.0}}
        };
    }
    else if (configName == "jupiter_moons") {
        bodies = {
            {1.898e27, {    0.0,   0.0 }, {   0.0,    0.0   }},
            {8.93e22,  {4.22e8,   0.0 }, {   0.0, 17320.0   }},
            {4.8e22,   {6.71e8,   0.0 }, {   0.0, 13740.0   }},
            {1.48e23,  {1.07e9,   0.0 }, {   0.0, 10870.0   }},
            {1.08e23,  {1.88e9,   0.0 }, {   0.0,  8200.0   }}
        };
    }
    else if (configName == "solar_system") {
        bodies = {
            {1.989e30, {      0.0,      0.0 }, {   0.0,     0.0   }},
            {3.285e23, { 5.79e10,      0.0 }, {   0.0,  47400.0   }},
            {4.867e24, {1.082e11,      0.0 }, {   0.0,  35020.0   }},
            {5.972e24, {1.496e11,      0.0 }, {   0.0,  29780.0   }},
            {6.39e23,  {2.279e11,      0.0 }, {   0.0,  24130.0   }},
            {1.898e27, {7.785e11,      0.0 }, {   0.0,  13070.0   }}
        };
    }
    else if (configName == "milky_way") {
        std::mt19937_64 rng(std::random_device{}());
        std::uniform_real_distribution<double> rDist(0, 5e20);
        std::uniform_real_distribution<double> aDist(0, 2 * std::acos(-1));
        std::uniform_real_distribution<double> mDist(1e30, 1e32);
        const int Nstars = 100;
        for (int i = 0; i < Nstars; ++i) {
            double r = rDist(rng);
            double th = aDist(rng);
            double x = r * std::cos(th), y = r * std::sin(th);
            double m = mDist(rng);
            bodies.emplace_back(m, Vec{x, y}, Vec{0.0, 0.0});
        }
    }
    else if (configName == "large_random_simulation") {
        std::mt19937_64 rng(std::random_device{}());
        std::uniform_real_distribution<double> posD(-1e10, 1e10);
        std::uniform_real_distribution<double> velD(-1e4, 1e4);
        std::uniform_real_distribution<double> massD(1e20, 1e25);
        const int Nrand = 50;
        for (int i = 0; i < Nrand; ++i) {
            double x = posD(rng), y = posD(rng);
            double vx = velD(rng), vy = velD(rng);
            double m = massD(rng);
            bodies.emplace_back(m, Vec{x, y}, Vec{vx, vy});
        }
    }
    else if (configName == "two_body_test") {
        bodies = {
            {1e7, {-1.0, 0.0}, {0.0, 0.0}},
            {1e7, { 1.0, 0.0}, {0.0, 0.0}}
        };
    }
    else {
        std::cerr << "Unexpected config: " << configName << "\n";
        return 1;
    }

    int N = bodies.size();

    vector<vector<Vec>> forces_thread(T, vector<Vec>(N, {0.0, 0.0}));
    vector<vector<Vec>> history;
    history.reserve(cfg.steps + 1);
    {
        vector<Vec> frame;
        frame.reserve(N);
        for (auto &b : bodies) {
            frame.push_back(b.pos);
        }
        history.push_back(frame);
    }

    std::cout << "Starting parallel integration with " << T
              << " threads, " << cfg.steps << " steps.\n";

    for (int step = 0; step < cfg.steps; ++step) {
        for (int t = 0; t < T; ++t) {
            for (int i = 0; i < N; ++i) {
                forces_thread[t][i].x = 0.0;
                forces_thread[t][i].y = 0.0;
            }
        }

        auto compute_forces = [&](int t_id, int i_start, int i_end) {
            for (int i = i_start; i < i_end; ++i) {
                for (int j = i + 1; j < N; ++j) {
                    double dx = bodies[j].pos.x - bodies[i].pos.x;
                    double dy = bodies[j].pos.y - bodies[i].pos.y;
                    double d2 = dx*dx + dy*dy + 1e-12;
                    double inv = 1.0 / std::sqrt(d2);
                    double F   = G * bodies[i].m * bodies[j].m * inv * inv;
                    double fx  = F * dx * inv;
                    double fy  = F * dy * inv;

                    forces_thread[t_id][i].x +=  fx;
                    forces_thread[t_id][i].y +=  fy;
                    forces_thread[t_id][j].x += -fx;
                    forces_thread[t_id][j].y += -fy;
                }
            }
        };

        vector<std::thread> workers;
        workers.reserve(T);
        int chunk = (N + T - 1) / T;
        for (int t = 0; t < T; ++t) {
            int i_start = t * chunk;
            int i_end   = std::min(i_start + chunk, N);
            if (i_start < i_end) {
                workers.emplace_back(compute_forces, t, i_start, i_end);
            }
        }
        for (auto &th : workers) {
            th.join();
        }

        for (int i = 0; i < N; ++i) {
            bodies[i].force.x = 0.0;
            bodies[i].force.y = 0.0;
            for (int t = 0; t < T; ++t) {
                bodies[i].force.x += forces_thread[t][i].x;
                bodies[i].force.y += forces_thread[t][i].y;
            }
        }

        auto update_positions = [&](int i_start, int i_end) {
            for (int i = i_start; i < i_end; ++i) {
                bodies[i].vel.x += (bodies[i].force.x / bodies[i].m) * cfg.timestep;
                bodies[i].vel.y += (bodies[i].force.y / bodies[i].m) * cfg.timestep;
                bodies[i].pos.x += bodies[i].vel.x * cfg.timestep;
                bodies[i].pos.y += bodies[i].vel.y * cfg.timestep;
            }
        };

        workers.clear();
        for (int t = 0; t < T; ++t) {
            int i_start = t * chunk;
            int i_end   = std::min(i_start + chunk, N);
            if (i_start < i_end) {
                workers.emplace_back(update_positions, i_start, i_end);
            }
        }
        for (auto &th : workers) {
            th.join();
        }

        vector<Vec> frame;
        frame.reserve(N);
        for (auto &b : bodies) {
            frame.push_back(b.pos);
        }
        history.push_back(frame);
    }

    std::cout << "Finished integration, building " << history.size() << " frames.\n";

    vector<Image> frames;
    frames.reserve(history.size());
    vector<string> colors = {"red","blue","green","orange","purple"};

    for (size_t t = 0; t < history.size(); ++t) {
        Image img(Geometry(cfg.width, cfg.height), "white");
        img.strokeColor("black");
        img.strokeWidth(1);

        for (int i = 0; i < N; ++i) {
            double nx = history[t][i].x * cfg.distanceScale + 0.5;
            double ny = history[t][i].y * cfg.distanceScale + 0.5;
            int px = int(nx * (cfg.width - 1));
            int py = cfg.height - 1 - int(ny * (cfg.height - 1));

            img.fillColor(colors[i % colors.size()]);
            img.draw(DrawableCircle(px, py, px + 5, py));
        }

        img.fillColor("black");
        std::ostringstream ss;
        ss << "t=" << std::fixed << std::setprecision(1)
           << (t * cfg.timestep) << "s";
        img.annotate(ss.str(), NorthWestGravity);

        img.animationDelay(10);
        img.animationIterations(0);

        frames.push_back(img);
    }

    std::cout << "Building GIF, writing to disk…\n";

    auto now = std::chrono::system_clock::now();
    std::time_t ti = std::chrono::system_clock::to_time_t(now);
    std::tm* tm    = std::localtime(&ti);
    std::ostringstream fn;
    fn << "animation_parallel_" << configName << "_"
       << std::put_time(tm, "%Y%m%d_%H%M%S") << ".gif";

    writeImages(frames.begin(), frames.end(), fn.str());
    std::cout << "Done writing GIF (" << fn.str() << "). Exiting.\n";

    return 0;
}
