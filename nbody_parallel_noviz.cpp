// nbody_parallel_noviz.cpp

#include <iostream>
#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <cmath>
#include <random>
#include "config.hpp"
#include "body.hpp"

using std::vector;
using std::string;

static constexpr double G_CONST = 6.67430e-11;
static constexpr double PI      = 3.14159265358979323846;

// a reusable, sense-reversing barrier for C++11/14/17
class Barrier {
    std::mutex              mtx;
    std::condition_variable cv;
    int                     count;
    const int               threshold;
    int                     generation = 0;

public:
    explicit Barrier(int n)
      : count(n), threshold(n) {}

    void wait() {
        std::unique_lock<std::mutex> lk(mtx);
        int gen = generation;
        if (--count == 0) {
            // last thread arrives
            generation++;
            count = threshold;
            cv.notify_all();
        } else {
            // wait until generation changes
            cv.wait(lk, [&]{ return gen != generation; });
        }
    }
};

int main(int argc, char** argv) {
    if (argc < 2 || argc > 3) {
        std::cerr << "Usage: " << argv[0]
                  << " <scenario> [num_threads]\n";
        return 1;
    }
    string configName = argv[1];
    int T = (argc == 3
             ? std::stoi(argv[2])
             : std::thread::hardware_concurrency());
    if (T < 1) T = 1;

    Config cfg(configName);

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
        for (int i=0; i<100; ++i) {
            double r  = rD(rng), th = aD(rng);
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
        for (int i=0; i<50; ++i) {
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

    int N = (int)bodies.size();
    int chunk = (N + T - 1) / T;

    // per-thread force buffers
    vector<vector<Vec>> forces_thread(T, vector<Vec>(N, {0,0}));

    // three barriers: after compute, after reduce, after update
    Barrier b_compute(T), b_reduce(T), b_update(T);

    // worker lambda
    auto worker = [&](int tid) {
        int start = tid * chunk;
        int end   = std::min(start + chunk, N);

        for (int step = 0; step < cfg.steps; ++step) {
            // (1) zero local buffer
            for (int i = start; i < end; ++i)
                forces_thread[tid][i] = {0,0};

            // (2) compute forces on [start..end)
            for (int i = start; i < end; ++i) {
                for (int j = i + 1; j < N; ++j) {
                    double dx = bodies[j].pos.x - bodies[i].pos.x;
                    double dy = bodies[j].pos.y - bodies[i].pos.y;
                    double d2 = dx*dx + dy*dy + 1e-12;
                    double inv = 1.0/std::sqrt(d2);
                    double F   = G_CONST * bodies[i].m * bodies[j].m * inv*inv;
                    double fx  = F * dx * inv;
                    double fy  = F * dy * inv;
                    forces_thread[tid][i].x += fx;
                    forces_thread[tid][i].y += fy;
                    forces_thread[tid][j].x -= fx;
                    forces_thread[tid][j].y -= fy;
                }
            }

            // (3) wait for all threads to finish compute
            b_compute.wait();

            // (4) thread 0 does the reduction into bodies[].force
            if (tid == 0) {
                for (int i = 0; i < N; ++i) {
                    bodies[i].force = {0,0};
                    for (int t = 0; t < T; ++t) {
                        bodies[i].force.x += forces_thread[t][i].x;
                        bodies[i].force.y += forces_thread[t][i].y;
                    }
                }
            }

            // (5) wait for reduction
            b_reduce.wait();

            // (6) update velocities & positions on [start..end)
            for (int i = start; i < end; ++i) {
                bodies[i].vel.x += (bodies[i].force.x / bodies[i].m) * cfg.timestep;
                bodies[i].vel.y += (bodies[i].force.y / bodies[i].m) * cfg.timestep;
                bodies[i].pos.x += bodies[i].vel.x * cfg.timestep;
                bodies[i].pos.y += bodies[i].vel.y * cfg.timestep;
            }

            // (7) wait for all threads to finish update
            b_update.wait();
        }
    };

    // spawn T workers
    vector<std::thread> threads;
    threads.reserve(T);
    for (int t = 0; t < T; ++t)
        threads.emplace_back(worker, t);

    // join once at the end
    for (auto &th : threads)
        th.join();

    return 0;
}
