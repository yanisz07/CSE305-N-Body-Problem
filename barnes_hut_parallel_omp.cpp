// barnes_hut_parallel_omp.cpp

#include "body.hpp"
#include "config.hpp"

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <omp.h>

using std::vector;
using std::string;

// Constants
static constexpr double G         = 6.67430e-11;
static constexpr double softening = 1e-12;

#ifndef THETA
static constexpr double theta = 0.5;
#else
static constexpr double theta = THETA;
#endif

// 2D quadtree node for Barnes–Hut
struct QuadNode {
    Vec center;       // center of this node’s region
    double halfDim;   // half width of region
    double mass;      // total mass in node
    Vec    com;       // center of mass
    int    bodyIndex; // if leaf with one body, its index; else -1
    bool   isLeaf;
    QuadNode* children[4];

    QuadNode(const Vec& c, double h)
      : center(c), halfDim(h), mass(0.0), com{0.0,0.0}, bodyIndex(-1), isLeaf(true)
    {
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
            case 0: return Vec{center.x - off, center.y + off}; // NW
            case 1: return Vec{center.x + off, center.y + off}; // NE
            case 2: return Vec{center.x - off, center.y - off}; // SW
            case 3: return Vec{center.x + off, center.y - off}; // SE
        }
        return Vec{0.0,0.0};
    }

    void insert(int idx, vector<Body>& bodies) {
        const Vec& p = bodies[idx].pos;
        if (!contains(p)) return;

        if (isLeaf && bodyIndex < 0) {
            bodyIndex = idx;
            mass      = bodies[idx].m;
            com       = bodies[idx].pos;
            return;
        }

        if (isLeaf) {
            // subdivide
            int old = bodyIndex;
            bodyIndex = -1;
            isLeaf = false;
            insertIntoChild(old, bodies);
        }
        insertIntoChild(idx, bodies);
        recomputeCOM();
    }

  private:
    void insertIntoChild(int idx, vector<Body>& bodies) {
        int q = getQuadrant(bodies[idx].pos);
        if (!children[q])
            children[q] = new QuadNode(childCenter(q), halfDim / 2.0);
        children[q]->insert(idx, bodies);
    }

    void recomputeCOM() {
        mass = 0.0;
        com  = Vec{0.0,0.0};
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

// Recursive force accumulation using Barnes–Hut criterion
void computeForceOnBody(int i, QuadNode* node, vector<Body>& bodies) {
    if (!node || node->mass <= 0.0) return;

    // leaf with single body
    if (node->isLeaf && node->bodyIndex >= 0) {
        int j = node->bodyIndex;
        if (j == i) return;
        double dx = bodies[j].pos.x - bodies[i].pos.x;
        double dy = bodies[j].pos.y - bodies[i].pos.y;
        double d2 = dx*dx + dy*dy + softening;
        double d  = std::sqrt(d2);
        if (d > 0.0) {
            double F = G * bodies[i].m * bodies[j].m / d2;
            bodies[i].force.x += F * dx / d;
            bodies[i].force.y += F * dy / d;
        }
        return;
    }

    // internal node or group
    double dx = node->com.x - bodies[i].pos.x;
    double dy = node->com.y - bodies[i].pos.y;
    double d2 = dx*dx + dy*dy + softening;
    double d  = std::sqrt(d2);
    double s  = node->halfDim * 2.0;

    if (s / d < theta) {
        // approximate as single body
        double F = G * bodies[i].m * node->mass / d2;
        bodies[i].force.x += F * dx / d;
        bodies[i].force.y += F * dy / d;
    } else {
        // recurse into children
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

    // Set OpenMP thread count
    if (argc == 3) {
        int T = std::max(1, std::stoi(argv[2]));
        omp_set_num_threads(T);
    }

    Config cfg(configName);

    // Initialize bodies according to config
    vector<Body> bodies;
    if (configName == "benchmark_fixed") {
        const int NX = 20, NY = 25;
        const double span = 2e9;
        const double dx = span/(NX-1), dy = span/(NY-1);
        bodies.reserve(NX*NY);
        for (int ix = 0; ix < NX; ++ix)
            for (int iy = 0; iy < NY; ++iy)
                bodies.emplace_back(
                    1e25,
                    Vec{-1e9 + ix*dx, -1e9 + iy*dy},
                    Vec{0.0,0.0}
                );
    }
    // ... other configs (earth_moon, solar_system, etc.) same as before ...
    else {
        std::cerr << "Unknown config: " << configName << "\n";
        return 1;
    }

    int N = bodies.size();
    std::cout << "Starting OpenMP Barnes-Hut with "
              << omp_get_max_threads() << " threads, "
              << cfg.steps << " steps.\n";

    for (int step = 0; step < cfg.steps; ++step) {
        // zero forces
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < N; ++i)
            bodies[i].force = Vec{0.0, 0.0};

        // compute bounding box
        double minX = bodies[0].pos.x, maxX = minX;
        double minY = bodies[0].pos.y, maxY = minY;
        for (int i = 1; i < N; ++i) {
            minX = std::min(minX, bodies[i].pos.x);
            maxX = std::max(maxX, bodies[i].pos.x);
            minY = std::min(minY, bodies[i].pos.y);
            maxY = std::max(maxY, bodies[i].pos.y);
        }
        double width = std::max(maxX - minX, maxY - minY);
        Vec center{ (minX + maxX)/2.0, (minY + maxY)/2.0 };
        double halfDim = 0.5 * width + 1e-5;

        // build quadtree
        QuadNode* root = new QuadNode(center, halfDim);
        for (int i = 0; i < N; ++i)
            root->insert(i, bodies);

        // parallel force computation
        #pragma omp parallel for schedule(dynamic, 16)
        for (int i = 0; i < N; ++i)
            computeForceOnBody(i, root, bodies);

        delete root;

        // parallel integrate
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < N; ++i) {
            auto &b = bodies[i];
            b.vel.x += (b.force.x / b.m) * cfg.timestep;
            b.vel.y += (b.force.y / b.m) * cfg.timestep;
            b.pos.x += b.vel.x * cfg.timestep;
            b.pos.y += b.vel.y * cfg.timestep;
        }
    }

    std::cout << "Finished integration.\n";
    return 0;
}
