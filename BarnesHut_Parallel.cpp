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
static constexpr double theta = 0.5;
static constexpr double softening = 1e-12;

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
    int T = (argc==3 ? std::stoi(argv[2]) : std::thread::hardware_concurrency());
    if (T<1) T=1;

    Config cfg(configName);
    InitializeMagick(*argv);

    vector<Body> bodies;
    if (configName == "earth_moon") {
        bodies = {
            {5.972e24, {0,0}, {0,0}},
            {7.34767309e22, {3.84e8,0}, {0,1022}}
        };
    }
    else if (configName == "jupiter_moons") {
        bodies = {
            {1.898e27, {0,0}, {0,0}},
            {8.93e22, {4.22e8,0}, {0,17320}},
            {4.8e22, {6.71e8,0}, {0,13740}},
            {1.48e23, {1.07e9,0}, {0,10870}},
            {1.08e23, {1.88e9,0}, {0,8200}}
        };
    }
    else if (configName == "solar_system") {
        bodies = {
            {1.989e30, {0,0}, {0,0}},
            {3.285e23, {5.79e10,0}, {0,47400}},
            {4.867e24, {1.082e11,0}, {0,35020}},
            {5.972e24, {1.496e11,0}, {0,29780}},
            {6.39e23, {2.279e11,0}, {0,24130}},
            {1.898e27, {7.785e11,0}, {0,13070}}
        };
    }
    else if (configName == "milky_way") {
        std::mt19937_64 rng(std::random_device{}());
        std::uniform_real_distribution<double> rd(0,5e20), ad(0,2*M_PI), md(1e30,1e32);
        const int Nstars=100;
        for (int i=0;i<Nstars;++i) {
            double r=rd(rng), th=ad(rng);
            bodies.emplace_back(md(rng), Vec{r*std::cos(th),r*std::sin(th)}, Vec{0,0});
        }
    }
    else if (configName == "large_random_simulation") {
        std::mt19937_64 rng(std::random_device{}());
        std::uniform_real_distribution<double> pD(-1e10,1e10), vD(-1e4,1e4), mD(1e20,1e25);
        const int Nrand=50;
        for(int i=0;i<Nrand;++i)
            bodies.emplace_back(mD(rng), Vec{pD(rng),pD(rng)}, Vec{vD(rng),vD(rng)});
    }
    else if (configName == "two_body_test") {
        bodies = {{1e7,{-1,0},{0,0}}, {1e7,{1,0},{0,0}}};
    } else {
        std::cerr<<"Unexpected config"<<"\n";
        return 1;
    }

    int N = bodies.size();
    vector<vector<Vec>> history;
    history.reserve(cfg.steps+1);
    {
        vector<Vec> frame; frame.reserve(N);
        for (auto &b: bodies) frame.push_back(b.pos);
        history.push_back(frame);
    }

    std::cout<<"Starting parallel Barnes-Hut with "<<T<<" threads, "<<cfg.steps<<" steps.\n";

    for (int step=0; step<cfg.steps; ++step) {
        for (auto &b: bodies) { b.force = Vec{0,0}; }

        double minX=bodies[0].pos.x, maxX=minX, minY=bodies[0].pos.y, maxY=minY;
        for (int i=1;i<N;++i) {
            minX = std::min(minX, bodies[i].pos.x);
            maxX = std::max(maxX, bodies[i].pos.x);
            minY = std::min(minY, bodies[i].pos.y);
            maxY = std::max(maxY, bodies[i].pos.y);
        }
        double width = std::max(maxX-minX, maxY-minY);
        Vec center{(minX+maxX)/2.0, (minY+maxY)/2.0};
        double halfDim = 0.5*width + 1e-5;

        QuadNode* root = new QuadNode(center, halfDim);
        for (int i=0; i<N; ++i) root->insert(i, bodies);

        vector<std::thread> workers;
        workers.reserve(T);
        int chunk = (N + T - 1)/T;
        for (int t=0; t<T; ++t) {
            int start = t*chunk;
            int end   = std::min(start+chunk, N);
            if (start<end) {
                workers.emplace_back([start,end,root,&bodies]{
                    for (int i=start; i<end; ++i)
                        computeForceOnBody(i, root, bodies);
                });
            }
        }
        for (auto &th: workers) th.join();

        delete root;

        workers.clear();
        for (int t=0; t<T; ++t) {
            int start = t*chunk;
            int end   = std::min(start+chunk, N);
            if (start<end) {
                workers.emplace_back([start,end,&bodies,&cfg]{
                    for (int i=start; i<end; ++i) {
                        auto &b = bodies[i];
                        b.vel.x += (b.force.x / b.m) * cfg.timestep;
                        b.vel.y += (b.force.y / b.m) * cfg.timestep;
                        b.pos.x += b.vel.x * cfg.timestep;
                        b.pos.y += b.vel.y * cfg.timestep;
                    }});
            }
        }
        for (auto &th: workers) th.join();

        vector<Vec> frame; frame.reserve(N);
        for (auto &b: bodies) frame.push_back(b.pos);
        history.push_back(frame);
    }

    vector<Image> frames;
    frames.reserve(history.size());
    vector<string> colors = {"red","blue","green","orange","purple"};
    for (size_t t=0; t<history.size(); ++t) {
        Image img(Geometry(cfg.width, cfg.height), "white");
        img.strokeColor("black"); img.strokeWidth(1);
        for (int i=0; i<N; ++i) {
            double nx = history[t][i].x * cfg.distanceScale + 0.5;
            double ny = history[t][i].y * cfg.distanceScale + 0.5;
            int px = int(nx*(cfg.width-1));
            int py = cfg.height-1 - int(ny*(cfg.height-1));
            img.fillColor(colors[i % colors.size()]);
            img.draw(DrawableCircle(px,py, px+5,py));
        }
        std::ostringstream ss;
        ss<<"t="<<std::fixed<<std::setprecision(1)
          <<(t*cfg.timestep)<<"s";
        img.annotate(ss.str(), NorthWestGravity);
        img.animationDelay(10); img.animationIterations(0);
        frames.push_back(img);
    }

    std::ostringstream fn;
    auto now = std::chrono::system_clock::now();
    std::time_t ti = std::chrono::system_clock::to_time_t(now);
    fn<<"animation_barnes_hut_"<<configName<<"_"<<std::put_time(std::localtime(&ti), "%Y%m%d_%H%M%S")<<".gif";
    writeImages(frames.begin(), frames.end(), fn.str());
    std::cout<<"Done writing GIF ("<<fn.str()<<").\n";
    return 0;
}
