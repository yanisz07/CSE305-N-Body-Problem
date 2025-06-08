// nbody_parallel.cpp

#include <Magick++.h>
#include <iostream>
#include <vector>
#include <deque>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <cmath>
#include <random>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <ctime>

#include "config.hpp"
#include "body.hpp"

using std::vector;
using std::string;

// physical constants
static constexpr double G_CONST = 6.67430e-11;
static constexpr double PI      = 3.14159265358979323846;


class Barrier {
    std::mutex              mtx;
    std::condition_variable cv;
    int                     count, threshold, generation = 0;
public:
    explicit Barrier(int n): count(n), threshold(n) {}
    void wait() {
        std::unique_lock<std::mutex> lk(mtx);
        int gen = generation;
        if (--count == 0) {
            generation++;
            count = threshold;
            cv.notify_all();
        } else {
            cv.wait(lk, [&]{ return gen != generation; });
        }
    }
};

// per‐frame normalized positions for IO thread
struct FrameData {
    vector<Vec> normalized_positions;
};

// shared queue for frames
static std::mutex              queue_mutex;
static std::condition_variable queue_cv;
static std::deque<FrameData>   render_queue;
static bool                    physics_done = false;

int main(int argc, char** argv) {
    if (argc < 2 || argc > 3) {
        std::cerr << "Usage: " << argv[0]
                  << " <scenario> [num_threads]\n";
        return 1;
    }
    // read scenario & thread count
    string configName = argv[1];
    int T = (argc == 3
             ? std::stoi(argv[2])
             : std::thread::hardware_concurrency());
    if (T < 1) T = 1;

    // load config
    Config cfg(configName);

    // init Magick++
    Magick::InitializeMagick(*argv);

    // build timestamped GIF filename
    auto now = std::chrono::system_clock::now();
    std::time_t ti = std::chrono::system_clock::to_time_t(now);
    std::tm* tm = std::localtime(&ti);
    std::ostringstream fn;
    fn << "gif_parallel_" << configName << "_"
       << std::put_time(tm, "%Y%m%d_%H%M%S") << ".gif";
    string gif_filename = fn.str();

    // feedback
    std::cout << "Scenario:       " << configName << "\n"
              << "Threads:        " << T << "\n"
              << "Steps:          " << cfg.steps
              << ", Timestep:    " << cfg.timestep << " s\n"
              << "DistanceScale:  " << cfg.distanceScale << "\n"
              << "Output GIF:     " << gif_filename << "\n\n";

    // 1) define bodies
    vector<Body> bodies;
    if (configName == "benchmark_fixed") {
        const int NX = 25, NY = 25;
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

    int N     = bodies.size();
    int chunk = (N + T - 1) / T;

    // per-thread force buffers
    vector<vector<Vec>> forces_thread(T, vector<Vec>(N, {0,0}));

    // three barriers: after compute, after reduce, after update
    Barrier b_compute(T), b_reduce(T), b_update(T);

    // — IO thread (consumer) —
    std::thread io_thread([&](){
        std::cout << "[IO] thread started\n";
        vector<Magick::Image> frames;
        frames.reserve(cfg.steps + 1);
        while (true) {
            std::unique_lock<std::mutex> lk(queue_mutex);
            queue_cv.wait(lk, [&]{
                return !render_queue.empty() || physics_done;
            });
            if (!render_queue.empty()) {
                FrameData F = std::move(render_queue.front());
                render_queue.pop_front();
                lk.unlock();

                // draw frame
                Magick::Image img(
                  Magick::Geometry(cfg.width, cfg.height),
                  "white"
                );
                img.strokeColor("black");
                img.strokeWidth(1);
                for (int i = 0; i < N; ++i) {
                    int px = int(F.normalized_positions[i].x * (cfg.width-1));
                    int py = cfg.height-1 
                           - int(F.normalized_positions[i].y * (cfg.height-1));
                    img.fillColor("red");
                    img.draw(Magick::DrawableCircle(px,py, px+5,py));
                }
                img.fillColor("black");
                std::ostringstream ss;
                ss << "t=" << std::fixed << std::setprecision(1)
                   << (frames.size() * cfg.timestep) << "s";
                img.annotate(ss.str(), Magick::NorthWestGravity);
                img.animationDelay(10);
                img.animationIterations(0);

                frames.push_back(std::move(img));
                continue;
            }
            if (physics_done && render_queue.empty()) {
                break;
            }
        }
        std::cout << "[IO] writing GIF...\n";
        Magick::writeImages(frames.begin(), frames.end(), gif_filename);
        std::cout << "[IO] GIF done\n";
    });

    // — Physics worker threads —
    auto worker = [&](int tid){
        int start = tid * chunk;
        int end   = std::min(start + chunk, N);
        for (int step = 0; step < cfg.steps; ++step) {
            // zero local buffer
            for (int i = start; i < end; ++i)
                forces_thread[tid][i] = {0,0};

            // compute forces
            for (int i = start; i < end; ++i) {
                for (int j = i+1; j < N; ++j) {
                    double dx = bodies[j].pos.x - bodies[i].pos.x;
                    double dy = bodies[j].pos.y - bodies[i].pos.y;
                    double d2 = dx*dx + dy*dy + 1e-12;
                    double inv = 1.0/std::sqrt(d2);
                    double F   = G_CONST*bodies[i].m*bodies[j].m*inv*inv;
                    double fx  = F*dx*inv;
                    double fy  = F*dy*inv;
                    forces_thread[tid][i].x += fx;
                    forces_thread[tid][i].y += fy;
                    forces_thread[tid][j].x -= fx;
                    forces_thread[tid][j].y -= fy;
                }
            }
            b_compute.wait();

            // reduction by thread 0
            if (tid == 0) {
                for (int i = 0; i < N; ++i) {
                    bodies[i].force = {0,0};
                    for (int t = 0; t < T; ++t) {
                        bodies[i].force.x += forces_thread[t][i].x;
                        bodies[i].force.y += forces_thread[t][i].y;
                    }
                }
            }
            b_reduce.wait();

            // update positions
            for (int i = start; i < end; ++i) {
                bodies[i].vel.x += (bodies[i].force.x / bodies[i].m) * cfg.timestep;
                bodies[i].vel.y += (bodies[i].force.y / bodies[i].m) * cfg.timestep;
                bodies[i].pos.x += bodies[i].vel.x * cfg.timestep;
                bodies[i].pos.y += bodies[i].vel.y * cfg.timestep;
            }
            b_update.wait();

            // push frame data
            FrameData F;
            F.normalized_positions.resize(N);
            for (int i = 0; i < N; ++i) {
                F.normalized_positions[i].x 
                  = bodies[i].pos.x * cfg.distanceScale + 0.5;
                F.normalized_positions[i].y 
                  = bodies[i].pos.y * cfg.distanceScale + 0.5;
            }
            {
                std::lock_guard<std::mutex> lk(queue_mutex);
                render_queue.push_back(std::move(F));
            }
            queue_cv.notify_one();
        }
    };

    // launch physics threads
    vector<std::thread> phys;
    phys.reserve(T);
    for (int t = 0; t < T; ++t)
        phys.emplace_back(worker, t);

    // join physics
    for (auto &th : phys) th.join();

    // signal IO & finish
    {
        std::lock_guard<std::mutex> lk(queue_mutex);
        physics_done = true;
    }
    queue_cv.notify_one();
    io_thread.join();

    std::cout << "Done. GIF: " << gif_filename << "\n";
    return 0;
}
