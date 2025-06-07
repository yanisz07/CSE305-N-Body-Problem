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
#include <queue>
#include <mutex>
#include <condition_variable>


using namespace Magick;
using std::vector;
using std::string;

static constexpr double G = 6.67430e-11;

#ifndef THETA
static constexpr double theta = 0.5;
#else
static constexpr double theta = THETA;

#endif
static constexpr double softening = 1e-12;
static constexpr double PI = std::acos(-1.0);


struct QuadNode {
    Vec center;
    double halfDim;
    double mass;
    Vec com;
    int bodyIndex;
    bool isLeaf;
    QuadNode* children[4];

    QuadNode(const Vec& c, double h)
      : center(c), halfDim(h), mass(0.0), com{0.0,0.0}, bodyIndex(-1), isLeaf(true) {
        for (int i = 0; i < 4; ++i) children[i] = nullptr;
    }
    ~QuadNode() {
        for (int i = 0; i < 4; ++i) delete children[i];
    }

    bool contains(const Vec& p) const {
        return (p.x >= center.x - halfDim && p.x <= center.x + halfDim &&
                p.y >= center.y - halfDim && p.y <= center.y + halfDim);
    }

    int getQuadrant(const Vec& p) const {
        if (p.x < center.x) return (p.y >= center.y ? 0 : 2);
        else                return (p.y >= center.y ? 1 : 3);
    }

    Vec childCenter(int q) const {
        double off = halfDim / 2.0;
        switch(q) {
            case 0: return Vec{center.x - off, center.y + off};
            case 1: return Vec{center.x + off, center.y + off};
            case 2: return Vec{center.x - off, center.y - off};
            case 3: return Vec{center.x + off, center.y - off};
        }
        return Vec{0.0,0.0};
    }

    void insert(int idx, vector<Body>& bodies) {
        const Vec& p = bodies[idx].pos;
        if (!contains(p)) return;
        if (isLeaf && bodyIndex < 0) {
            bodyIndex = idx;
            mass = bodies[idx].m;
            com = bodies[idx].pos;
        } else {
            if (isLeaf) {
                int old = bodyIndex;
                bodyIndex = -1;
                isLeaf = false;
                insertIntoChild(old, bodies);
            }
            insertIntoChild(idx, bodies);
            recomputeCOM();
        }
    }

private:
    void insertIntoChild(int idx, vector<Body>& bodies) {
        int q = getQuadrant(bodies[idx].pos);
        if (!children[q]) children[q] = new QuadNode(childCenter(q), halfDim/2.0);
        children[q]->insert(idx, bodies);
    }
    void recomputeCOM() {
        mass = 0.0;
        com = Vec{0.0,0.0};
        for (int i = 0; i < 4; ++i) {
            if (children[i] && children[i]->mass > 0.0) {
                mass += children[i]->mass;
                com.x += children[i]->com.x * children[i]->mass;
                com.y += children[i]->com.y * children[i]->mass;
            }
        }
        if (mass > 0.0) {
            com.x /= mass;
            com.y /= mass;
        }
    }
};

void computeForceOnBody(int i, QuadNode* node, vector<Body>& bodies) {
    if (!node || node->mass <= 0.0) return;
    if (node->isLeaf && node->bodyIndex >= 0 && node->bodyIndex != i) {
        int j = node->bodyIndex;
        double dx = bodies[j].pos.x - bodies[i].pos.x;
        double dy = bodies[j].pos.y - bodies[i].pos.y;
        double dist2 = dx*dx + dy*dy + softening;
        double dist = std::sqrt(dist2);
        if (dist > 0.0) {
            double F = G * bodies[i].m * bodies[j].m / dist2;
            bodies[i].force.x += F * dx / dist;
            bodies[i].force.y += F * dy / dist;
        }
        return;
    }
    double dx = node->com.x - bodies[i].pos.x;
    double dy = node->com.y - bodies[i].pos.y;
    double dist2 = dx*dx + dy*dy + softening;
    double dist = std::sqrt(dist2);
    double s = node->halfDim * 2.0;
    if (s / dist < theta) {
        double F = G * bodies[i].m * node->mass / dist2;
        bodies[i].force.x += F * dx / dist;
        bodies[i].force.y += F * dy / dist;
    } else {
        for (int c = 0; c < 4; ++c) {
            if (node->children[c])
                computeForceOnBody(i, node->children[c], bodies);
        }
    }
}
int main(int argc, char** argv) {
    if (argc < 2 || argc > 3) {
        std::cerr << "Usage: " << argv[0]
                  << " <config_name> [num_threads]\n";
        return 1;
    }
    string configName = argv[1];
    int T = (argc == 3
             ? std::max(1, std::stoi(argv[2]))
             : std::max(1u, std::thread::hardware_concurrency()));

    Config cfg(configName);
    InitializeMagick(*argv);

    vector<Body> bodies;
    
    if (configName == "benchmark_fixed") {
        // 500 bodies on a 20×25 grid in ±1e9 m
        const int NX = 20, NY = 25;
        const double span = 2e9;
        const double dx = span/(NX-1), dy = span/(NY-1);
        bodies.reserve(NX*NY);
        for (int ix = 0; ix < NX; ++ix) {
            for (int iy = 0; iy < NY; ++iy) {
                bodies.emplace_back(
                  1e25,
                  Vec{-1e9 + ix*dx, -1e9 + iy*dy},
                  Vec{0,0}
                );
            }
        }
    }
    else if (configName == "earth_moon") {
        bodies = {
            {5.972e24,      { 0.0, 0.0 }, { 0.0, 0.0 }},
            {7.34767309e22, {3.84e8, 0.0 }, { 0.0, 1022.0 }}
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

    //thread-safe queue of position snapshots:
    std::queue<vector<Vec>>      posQueue;
    std::mutex                   mtxQueue;
    std::condition_variable      cvQueue;
    bool                         doneProducing = false;

    // The frame buffer (only I/O thread pushes here, main thread reads at the end)
    vector<Image> frames;
    std::mutex     mtxFrames;

    // I/O thread: consumes posQueue, and then produces frames[]
    std::thread ioThread([&]() {
        for (;;) {
            vector<Vec> snapshot;
            {
                std::unique_lock<std::mutex> lk(mtxQueue);
                cvQueue.wait(lk, [&]() {
                    return !posQueue.empty() || doneProducing;
                });
                if (posQueue.empty() && doneProducing)
                    break;
                snapshot = std::move(posQueue.front());
                posQueue.pop();
            }
            // Build one Image from snapshot:
            Image img(Geometry(cfg.width, cfg.height), "white");
            img.strokeColor("black");
            img.strokeWidth(1);
            static const vector<string> colors = {
                "red","blue","green","orange","purple"
            };
            for (int i = 0; i < N; ++i) {
                double nx = snapshot[i].x * cfg.distanceScale + 0.5;
                double ny = snapshot[i].y * cfg.distanceScale + 0.5;
                int px = int(nx * (cfg.width - 1));
                int py = cfg.height - 1 - int(ny * (cfg.height - 1));

                img.fillColor(colors[i % colors.size()]);
                img.draw(DrawableCircle(px, py, px + 5, py));
            }
            img.animationDelay(10);
            img.animationIterations(0);

            {
                std::lock_guard<std::mutex> lk(mtxFrames);
                frames.push_back(std::move(img));
            }
        }
    });

    std::cout << "Starting parallel Barnes-Hut with "
              << T << " threads, " << cfg.steps << " steps.\n";

    for (int step = 0; step < cfg.steps; ++step) {
        for (auto &b : bodies) {
            b.force.x = b.force.y = 0.0;
        }

        //build root quadtree
        double minX = bodies[0].pos.x, maxX = minX,
               minY = bodies[0].pos.y, maxY = minY;
        for (int i = 1; i < N; ++i) {
            minX = std::min(minX, bodies[i].pos.x);
            maxX = std::max(maxX, bodies[i].pos.x);
            minY = std::min(minY, bodies[i].pos.y);
            maxY = std::max(maxY, bodies[i].pos.y);
        }
        double width = std::max(maxX - minX, maxY - minY);
        Vec    center{ (minX + maxX) / 2.0, (minY + maxY) / 2.0 };
        double halfDim = 0.5 * width + 1e-5;

        QuadNode* root = new QuadNode(center, halfDim);
        for (int i = 0; i < N; ++i)
            root->insert(i, bodies);

        //parallel force computation
        {
            vector<std::thread> pool;
            pool.reserve(T);
            int chunk = (N + T - 1) / T;
            for (int t = 0; t < T; ++t) {
                int start = t * chunk;
                int end   = std::min(start + chunk, N);
                if (start >= end) continue;
                pool.emplace_back(
                    [start,end,root,&bodies]() {
                        for (int i = start; i < end; ++i)
                            computeForceOnBody(i, root, bodies);
                    });
            }
            for (auto &th : pool) th.join();
        }

        delete root;

        // parallel integration
        {
            vector<std::thread> pool;
            pool.reserve(T);
            int chunk = (N + T - 1) / T;
            for (int t = 0; t < T; ++t) {
                int start = t * chunk;
                int end   = std::min(start + chunk, N);
                if (start >= end) continue;
                pool.emplace_back(
                    [start,end,&bodies,&cfg]() {
                        for (int i = start; i < end; ++i) {
                            auto &b = bodies[i];
                            b.vel.x += (b.force.x / b.m) * cfg.timestep;
                            b.vel.y += (b.force.y / b.m) * cfg.timestep;
                            b.pos.x += b.vel.x * cfg.timestep;
                            b.pos.y += b.vel.y * cfg.timestep;
                        }
                    });
            }
            for (auto &th : pool) th.join();
        }

        //snapshot positions and notify I/O thread
        {
            vector<Vec> snapshot;
            snapshot.reserve(N);
            for (auto &b : bodies)
                snapshot.push_back(b.pos);

            {
                std::lock_guard<std::mutex> lk(mtxQueue);
                posQueue.push(std::move(snapshot));
            }
            cvQueue.notify_one();
        }
    }

    //signal I/O thread to finish
    {
        std::lock_guard<std::mutex> lk(mtxQueue);
        doneProducing = true;
    }
    cvQueue.notify_one();
    ioThread.join();

    //write the GIF
    std::cout << "Building and writing GIF (" << frames.size() << " frames)…\n";
    auto now = std::chrono::system_clock::now();
    std::time_t ti = std::chrono::system_clock::to_time_t(now);
    std::tm* tm    = std::localtime(&ti);
    std::ostringstream fn;
    fn << "animation_barnes_hut_" << configName << "_"
       << std::put_time(tm, "%Y%m%d_%H%M%S") << ".gif";

    writeImages(frames.begin(), frames.end(), fn.str());
    std::cout << "Done: " << fn.str() << "\n";
    return 0;
}
