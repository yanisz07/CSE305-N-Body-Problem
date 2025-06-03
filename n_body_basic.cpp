// n_body_basic.cpp
#include <Magick++.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <chrono>
#include <ctime>
#include <random>

using namespace Magick;
using std::vector;
using std::string;

// —— Global constants ——
const double G      = 6.67430e-11;    // gravitational constant

// —— Config class to choose zoom & timestep presets ——
struct Config {
    double distanceScale;  // maps physical positions [m] to normalized [0,1]
    double timestep;       // seconds per simulation step
    int    width, height;  // image dimensions
    int    steps;          // number of integration steps

    Config(const string& name) {
        width  = 800;
        height = 800;

        if (name == "earth_moon") {
            // Earth–Moon distance ≈ 3.84×10^8 m.
            // We map that to span from 0→1 in normalized coords (±0.5 around origin).
            distanceScale = 0.5 / 3.84e8;  // so moon at 3.84e8 m → normalized = 1.0
            timestep      = 3600;          // 1 hour per step
            steps         = 200;           // total steps → 200 h (~8 days)
        }
        else if (name == "jupiter_moons") {
            // Io orbital radius ≈ 4.22×10^8 m → map to ±0.5
            distanceScale = 0.5 / 4.22e8;
            timestep      = 3600;      // 1 h per step
            steps         = 200;       // ~200 h
        }
        else if (name == "solar_system") {
            // Earth–Sun ≈1.496×10^11 m → map to ±0.5
            distanceScale = 0.5 / 1.496e11;
            timestep      = 86400;     // 1 day per step
            steps         = 365;       // ~1 year
        }
        else if (name == "milky_way") {
            // approximate disk radius ~5×10^20 m → map to ±0.5
            distanceScale = 0.5 / 5e20;
            timestep      = 5e6 * 3.1536e7; // 5 Myr per step
            steps         = 100;            // 500 Myr
        }
        else if (name == "large_random_simulation") {
            // random bodies in ±1e10 m → map that to ±0.5
            distanceScale = 0.5 / 1e10;
            timestep      = 86400;     // 1 day per step
            steps         = 365;       // ~1 year
        }
        else if (name == "two_body_test") {
            // two bodies 2 m apart → map to ±0.5
            distanceScale = 0.5 / 2.0;  
            timestep      = 0.05;      // 0.05 s per step
            steps         = 200;       // 10 s total
        }
        else {
            std::cerr << "Invalid config name: " << name << "\n";
            std::exit(1);
        }
    }
};

// —— 2D vector & Body definitions ——
struct Vec { double x, y; };
struct Body {
    double m;
    Vec    pos, vel, force;
    Body(double m_, Vec p_, Vec v_)
      : m(m_), pos(p_), vel(v_), force{0,0} {}
};

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0]
                  << " <config_name>\n"
                  << "Valid configs: two_body_test, earth_moon, jupiter_moons,\n"
                  << "               solar_system, milky_way, large_random_simulation\n";
        return 1;
    }
    string configName = argv[1];
    Config cfg(configName);

    InitializeMagick(*argv);

    // —— Define bodies for each scenario ——
    vector<Body> bodies;
    if (configName == "earth_moon") {
        bodies = {
            // Earth at origin
            {5.972e24,        {     0.0,      0.0},   {  0.0,     0.0    }},
            // Moon 384×10^6 m away, initial orbital velocity ~1022 m/s
            {7.34767309e22,   {3.84e8,       0.0},   {  0.0,   1022.0   }}
        };
    }
    else if (configName == "jupiter_moons") {
        bodies = {
            {1.898e27,    {    0.0,       0.0},    {  0.0,      0.0   }},  
            {8.93e22,     {4.22e8,       0.0},    {  0.0,   17320.0   }},  
            {4.8e22,      {6.71e8,       0.0},    {  0.0,   13740.0   }},  
            {1.48e23,     {1.07e9,       0.0},    {  0.0,   10870.0   }},  
            {1.08e23,     {1.88e9,       0.0},    {  0.0,    8200.0   }}
        };
    }
    else if (configName == "solar_system") {
        bodies = {
            {1.989e30,    {      0.0,       0.0},   {  0.0,     0.0   }},
            {3.285e23,    { 5.79e10,       0.0},   {  0.0,  47400.0   }},
            {4.867e24,    {1.082e11,       0.0},   {  0.0,  35020.0   }},
            {5.972e24,    {1.496e11,       0.0},   {  0.0,  29780.0   }},
            {6.39e23,     {2.279e11,       0.0},   {  0.0,  24130.0   }},
            {1.898e27,    {7.785e11,       0.0},   {  0.0,  13070.0   }}
        };
    }
    else if (configName == "milky_way") {
        std::mt19937_64 rng(std::random_device{}());
        std::uniform_real_distribution<double> rDist(0, 5e20);
        std::uniform_real_distribution<double> aDist(0, 2*M_PI);
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
            {1e7,  {-1.0, 0.0}, {0.0, 0.0}},
            {1e7,  { 1.0, 0.0}, {0.0, 0.0}}
        };
    }
    else {
        std::cerr << "Unexpected config: " << configName << "\n";
        return 1;
    }

    int N = bodies.size();

    // —— Record trajectory ——
    vector<vector<Vec>> history;
    history.reserve(cfg.steps + 1);
    {
        vector<Vec> frame;
        frame.reserve(N);
        for (auto &b : bodies) frame.push_back(b.pos);
        history.push_back(frame);
    }

    // —— Integrate with simple Euler ——
    for (int step = 0; step < cfg.steps; ++step) {
        for (auto &b : bodies) b.force = {0, 0};

        for (int i = 0; i < N; ++i) {
            for (int j = i + 1; j < N; ++j) {
                double dx  = bodies[j].pos.x - bodies[i].pos.x;
                double dy  = bodies[j].pos.y - bodies[i].pos.y;
                double d2  = dx*dx + dy*dy + 1e-12;
                double inv = 1.0 / std::sqrt(d2);
                double F   = G * bodies[i].m * bodies[j].m * inv * inv;
                double fx  = F * dx * inv;
                double fy  = F * dy * inv;
                bodies[i].force.x += fx;
                bodies[i].force.y += fy;
                bodies[j].force.x -= fx;
                bodies[j].force.y -= fy;
            }
        }

        for (auto &b : bodies) {
            b.vel.x += (b.force.x / b.m) * cfg.timestep;
            b.vel.y += (b.force.y / b.m) * cfg.timestep;
            b.pos.x += b.vel.x * cfg.timestep;
            b.pos.y += b.vel.y * cfg.timestep;
        }

        vector<Vec> frame;
        frame.reserve(N);
        for (auto &b : bodies) frame.push_back(b.pos);
        history.push_back(frame);
    }

    std::cout << "Finished integration (" << cfg.steps << " steps), building frames\n";

    // —— Build frames (scale by cfg.distanceScale) ——
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
            img.draw(DrawableCircle(px, py, px+5, py));
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

    std::cout << "Finished building frames, now writing GIF…\n";

    // —— Timestamped filename ——
    auto now = std::chrono::system_clock::now();
    std::time_t ti = std::chrono::system_clock::to_time_t(now);
    std::tm* tm    = std::localtime(&ti);
    std::ostringstream fn;
    fn << "anime_" << configName << "_"
       << std::put_time(tm, "%Y%m%d_%H%M%S") << ".gif";

    writeImages(frames.begin(), frames.end(), fn.str());
    std::cout << "Done writing GIF (" << fn.str() << "). Exiting.\n";
    return 0;
}