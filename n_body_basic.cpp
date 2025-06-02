// n_body_basic.cpp
#include <Magick++.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <chrono>
#include <ctime>

using namespace Magick;
using std::vector;
using std::string;

// ——— Simulation parameters ———
const double G      = 6.67430e-11;    // gravitational constant
const double dt     = 0.05;           // time step (in seconds)
const int    STEPS = 200;             // number of steps (gif time = dt*STEPS seconds)
const int    WIDTH = 800;
const int    HEIGHT= 800;

// ——— 2D vector and Body definitions ———
struct Vec { double x,y; };
struct Body {
    double m;
    Vec    pos, vel, force;
    Body(double m_, Vec p_, Vec v_)
      : m(m_), pos(p_), vel(v_), force{0,0} {}
};

int main(int argc, char** argv) {
    InitializeMagick(*argv);

    // Define bodies here (hard-coded for now) ——
    // Each entry is { mass [kg], {x [m], y [m]}, {vx [m/s], vy [m/s]} }.
    vector<Body> bodies = {
        {1.0, {-1.0, 0.0}, { 0.0,  0.5}},
        {1.0, { 1.0, 0.0}, { 0.0, -0.5}}
    };
    int N = bodies.size();

    // —— Record the trajectory —
    vector<vector<Vec>> history;
    history.reserve(STEPS+1);
    {
        vector<Vec> frame;
        for(auto &b : bodies) frame.push_back(b.pos);
        history.push_back(frame);
    }

    // —— Integrate with simple Euler —
    for(int step=0; step<STEPS; ++step) {
        // zero forces
        for(auto &b : bodies) b.force = {0,0};

        // compute pairwise gravity
        for(int i=0;i<N;++i) {
            for(int j=i+1;j<N;++j) {
                double dx = bodies[j].pos.x - bodies[i].pos.x;
                double dy = bodies[j].pos.y - bodies[i].pos.y;
                double d2 = dx*dx + dy*dy + 1e-12;
                double inv = 1.0/std::sqrt(d2);
                double F   = G * bodies[i].m * bodies[j].m * inv*inv;
                double fx  = F * dx * inv;
                double fy  = F * dy * inv;
                bodies[i].force.x += fx;
                bodies[i].force.y += fy;
                bodies[j].force.x -= fx;
                bodies[j].force.y -= fy;
            }
        }

        // update velocities & positions
        for(auto &b : bodies){
            b.vel.x += (b.force.x/b.m) * dt;
            b.vel.y += (b.force.y/b.m) * dt;
            b.pos.x += b.vel.x * dt;
            b.pos.y += b.vel.y * dt;
        }

        // record
        vector<Vec> frame;
        for(auto &b : bodies) frame.push_back(b.pos);
        history.push_back(frame);
    }

    // DEBUG: report when integration is done
    std::cout << "▷ Finished integration, building " 
              << history.size() << " frames\n";

    // —— Build frames with a fixed window [-2,2] —
    vector<Image> frames;
    frames.reserve(history.size());
    vector<string> colors = {"red","blue"};

    for(size_t t=0; t<history.size(); ++t) {
        Image img(Geometry(WIDTH,HEIGHT),"white");
        img.strokeColor("black");
        img.strokeWidth(1);

        // draw each body
        for(int i=0;i<N;++i){
            double fx = (history[t][i].x + 2.0) / 4.0; // map [-2,2]→[0,1]
            double fy = (history[t][i].y + 2.0) / 4.0;
            int px = int(fx * (WIDTH-1));
            int py = HEIGHT-1 - int(fy * (HEIGHT-1));
            img.fillColor(colors[i%colors.size()]);
            img.draw(DrawableCircle(px,py, px+5,py));
        }

        // annotate time
        img.fillColor("black");
        std::ostringstream ss;
        ss<<"t="<<std::fixed<<std::setprecision(2)<<t*dt;
        img.annotate(ss.str(), NorthWestGravity);

        // slow & loop
        img.animationDelay(10);      // 0.1s per frame
        img.animationIterations(0);  // infinite

        frames.push_back(img);
    }

    // DEBUG: report before writing the GIF
    std::cout << "Finished building frames, now writing GIF…\n";

    // —— Generate a timestamped filename ——
    auto now = std::chrono::system_clock::now();
    std::time_t ti = std::chrono::system_clock::to_time_t(now);
    std::tm* tm = std::localtime(&ti);
    std::ostringstream fn;
    fn << "anime_"
       << std::put_time(tm, "%Y%m%d_%H%M%S")
       << ".gif";

    writeImages(frames.begin(), frames.end(), fn.str());

    std::cout << "Done writing GIF (" << fn.str() << "). Exiting.\n";
    return 0;
}
