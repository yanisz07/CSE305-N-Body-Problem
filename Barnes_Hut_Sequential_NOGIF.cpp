// Barnes_Hut_Sequential_NOGIF.cpp

#include "body.hpp"
#include "config.hpp"

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <chrono>
#include <ctime>

using std::vector;
using std::string;

// Gravitational constant and softening
static constexpr double G         = 6.67430e-11;
static constexpr double softening = 1e-12;

// Theta threshold (overrideable via -DTHETA)
#ifndef THETA
static constexpr double theta = 0.5;
#else
static constexpr double theta = THETA;
#endif

struct QuadNode {
    Vec center;
    double halfDim;
    double mass;
    Vec    com;
    int    bodyIndex;
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
            mass      = bodies[idx].m;
            com       = bodies[idx].pos;
            return;
        }
        if (isLeaf) {
            // split leaf
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
            children[q] = new QuadNode(childCenter(q), halfDim/2.0);
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
    // leaf with single body
    if (node->isLeaf && node->bodyIndex >= 0) {
        int j = node->bodyIndex;
        if (j != i) {
            double dx = bodies[j].pos.x - bodies[i].pos.x;
            double dy = bodies[j].pos.y - bodies[i].pos.y;
            double dist2 = dx*dx + dy*dy + softening;
            double dist  = std::sqrt(dist2);
            if (dist > 0.0) {
                double F = G * bodies[i].m * bodies[j].m / dist2;
                bodies[i].force.x += F * dx / dist;
                bodies[i].force.y += F * dy / dist;
            }
        }
        return;
    }
    // internal node or group
    double dx = node->com.x - bodies[i].pos.x;
    double dy = node->com.y - bodies[i].pos.y;
    double dist2 = dx*dx + dy*dy + softening;
    double dist  = std::sqrt(dist2);
    double s     = node->halfDim * 2.0;
    if (s / dist < theta) {
        double F = G * bodies[i].m * node->mass / dist2;
        bodies[i].force.x += F * dx / dist;
        bodies[i].force.y += F * dy / dist;
    } else {
        for (int c = 0; c < 4; ++c)
            if (node->children[c])
                computeForceOnBody(i, node->children[c], bodies);
    }
}

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <config_name>\n";
        return 1;
    }
    string configName = argv[1];
    Config cfg(configName);

    // Initialize bodies based on config
    vector<Body> bodies;
    // ... same config blocks as before ...
    // (omitted for brevity)

    int N = bodies.size();

    std::cout << "Starting Barnes-Hut integration (“" << configName
              << "”) with " << cfg.steps << " steps.\n";

    for (int step = 0; step < cfg.steps; ++step) {
        // zero forces
        for (auto &b : bodies) b.force = Vec{0.0,0.0};

        // build quadtree
        double minX = bodies[0].pos.x, maxX = minX;
        double minY = bodies[0].pos.y, maxY = minY;
        for (int i = 1; i < N; ++i) {
            minX = std::min(minX, bodies[i].pos.x);
            maxX = std::max(maxX, bodies[i].pos.x);
            minY = std::min(minY, bodies[i].pos.y);
            maxY = std::max(maxY, bodies[i].pos.y);
        }
        double width = std::max(maxX-minX, maxY-minY);
        Vec center{ (minX+maxX)/2.0, (minY+maxY)/2.0 };
        double halfDim = 0.5*width + 1e-5;

        QuadNode* root = new QuadNode(center, halfDim);
        for (int i = 0; i < N; ++i) root->insert(i, bodies);

        // compute forces
        for (int i = 0; i < N; ++i)
            computeForceOnBody(i, root, bodies);

        delete root;

        // integrate
        for (auto &b : bodies) {
            b.vel.x += (b.force.x / b.m) * cfg.timestep;
            b.vel.y += (b.force.y / b.m) * cfg.timestep;
            b.pos.x += b.vel.x * cfg.timestep;
            b.pos.y += b.vel.y * cfg.timestep;
        }
    }

    std::cout << "Finished Barnes-Hut integration.\n";
    return 0;
}
