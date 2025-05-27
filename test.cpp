#include <iostream>
#include <cmath>
#include <vector>
#include <thread>
#include <mutex>
#include <Magick++.h>
using namespace Magick;

std::mutex force_mutex;

void gather_input(int &n, std::vector<double> &masses, std::vector<Vector2D> &positions, std::vector<Vector2D> &velocities, double &time_step, double &total_time) {
    std::cout << "Enter the number of bodies to simulate: ";
    std::cin >> n;

    masses.resize(n);
    positions.resize(n);
    velocities.resize(n);

    for (int i = 0; i < n; ++i) {
        std::cout << "Enter data for body " << (i + 1) << ":\n";
        std::cout << "  Mass (kg): ";
        std::cin >> masses[i];
        while (masses[i] <= 0) {
            std::cout << "  Incorrect Mass. Enter again Mass (kg): ";
            std::cin >> masses[i];
        }
        std::cout << "  Initial x position (meters): ";
        std::cin >> positions[i].x;
        std::cout << "  Initial y position (meters): ";
        std::cin >> positions[i].y;
        std::cout << "  Initial x velocity (m/s): ";
        std::cin >> velocities[i].x;
        std::cout << "  Initial y velocity (m/s): ";
        std::cin >> velocities[i].y;
    }

    std::cout << "Enter the time step for the simulation (seconds): ";
    std::cin >> time_step;
    std::cout << "Enter the total simulation time (seconds): ";
    std::cin >> total_time;
}

void compute_forces_segment(const int n, const std::vector<double>& masses, const std::vector<Vector2D>& positions, std::vector<Vector2D>& forces, int start, int end, double G) {
    for (int i = start; i < end; ++i) {
        Vector2D force = {0, 0};
        for (int j = i + 1; j < n; ++j) {
            Vector2D delta = {positions[j].x - positions[i].x, positions[j].y - positions[i].y};
            double dist_squared = delta.x * delta.x + delta.y * delta.y;
            double dist = std::sqrt(dist_squared);
            double force_magnitude = G * masses[i] * masses[j] / dist_squared;
            Vector2D force_ij = {force_magnitude * delta.x / dist, force_magnitude * delta.y / dist};

            std::lock_guard<std::mutex> lock(force_mutex);
            forces[i].x += force_ij.x;
            forces[i].y += force_ij.y;
            forces[j].x -= force_ij.x; 
            forces[j].y -= force_ij.y; 
        }
    }
}

void compute_forces(const int n, const std::vector<double>& masses, const std::vector<Vector2D>& positions, std::vector<Vector2D>& forces, double G) {
    forces.assign(n, Vector2D{0, 0});
    int num_threads = std::thread::hardware_concurrency();
    int chunk_size = (n + num_threads - 1) / num_threads;
    std::vector<std::thread> threads;

    for (int i = 0; i < num_threads; ++i) {
        int start = i * chunk_size;
        int end = std::min(start + chunk_size, n);
        if (start < end) {
            threads.emplace_back(compute_forces_segment, n, std::ref(masses), std::ref(positions), std::ref(forces), start, end, G);
        }
    }

    for (auto& t : threads) {
        if (t.joinable()) {
            t.join();
        }
    }
}

void update_bodies_segment(int n, std::vector<double>& masses, std::vector<Vector2D>& positions, std::vector<Vector2D>& velocities, std::vector<Vector2D>& forces, double time_step, int start, int end) {
    for (int i = start; i < end; ++i) {
        velocities[i].x += forces[i].x / masses[i] * time_step;
        velocities[i].y += forces[i].y / masses[i] * time_step;
        positions[i].x += velocities[i].x * time_step;
        positions[i].y += velocities[i].y * time_step;
    }
}

void update_bodies(int n, std::vector<double>& masses, std::vector<Vector2D>& positions, std::vector<Vector2D>& velocities, std::vector<Vector2D>& forces, double time_step) {
    int num_threads = std::thread::hardware_concurrency();
    int chunk_size = (n + num_threads - 1) / num_threads;
    std::vector<std::thread> threads;

    for (int i = 0; i < num_threads; ++i) {
        int start = i * chunk_size;
        int end = std::min(start + chunk_size, n);
        if (start < end) {
            threads.emplace_back(update_bodies_segment, n, std::ref(masses), std::ref(positions), std::ref(velocities), std::ref(forces), time_step, start, end);
        }
    }

    for (auto& t : threads) {
        if (t.joinable()) {
            t.join();
        }
    }
}

void draw_arrow(Magick::Image &frame, int x1, int y1, double dx, double dy, const std::string &color) {
    double angle = std::atan2(dy, dx);
    const double arrow_length = 15;
    const double arrow_angle = M_PI / 6; // 30 degrees for each arrowhead wing

    int x2 = x1 + static_cast<int>(arrow_length * std::cos(angle));
    int y2 = y1 + static_cast<int>(arrow_length * std::sin(angle));

    int x3 = x2 - static_cast<int>(arrow_length/3 * std::cos(angle - arrow_angle));
    int y3 = y2 - static_cast<int>(arrow_length/3 * std::sin(angle - arrow_angle));
    int x4 = x2 - static_cast<int>(arrow_length/3 * std::cos(angle + arrow_angle));
    int y4 = y2 - static_cast<int>(arrow_length/3 * std::sin(angle + arrow_angle));

    frame.strokeColor(color);

    frame.draw(Magick::DrawableLine(x1, y1, x2, y2));
    frame.draw(Magick::DrawableLine(x2, y2, x3, y3));
    frame.draw(Magick::DrawableLine(x2, y2, x4, y4));
}

void save_frame(const std::vector<Vector2D>& positions, const std::vector<Vector2D>& velocities, const std::vector<Vector2D>& forces, int n, int t, std::vector<Magick::Image>& frames, double min_x, double max_x, double min_y, double max_y) {
    int width = 800;
    int height = 800;
    Magick::Image frame(Magick::Geometry(width, height), "white");
    frame.strokeWidth(2);
    std::vector<std::string> colors = {"red", "green", "blue", "yellow", "cyan", "magenta", "orange", "purple", "brown", "pink"};

    double range_x = max_x - min_x;
    double range_y = max_y - min_y;

    for (int i = 0; i < n; ++i) {
        int x = ((positions[i].x - min_x) / range_x) * width;
        int y = height - ((positions[i].y - min_y) / range_y) * height;

        frame.fillColor(colors[i % colors.size()]);
        frame.strokeColor("black");
        frame.draw(Magick::DrawableCircle(x, y, x + 3, y + 3));

        double vx = (velocities[i].x / range_x) * width;
        double vy = -(velocities[i].y / range_y) * height; // Negative because the y-axis is inverted in the image
        double Fx = (forces[i].x / range_x) * width;
        double Fy = -(forces[i].y / range_y) * height; // Negative because the y-axis is inverted in the image

        frame.fillColor("none"); // Reset fill color before drawing arrows
        draw_arrow(frame, x, y, vx, vy, colors[i % colors.size()]); // velocities are in the color of the object
        if (n > 1) {
            for (const auto& force : forces) {
                draw_arrow(frame, x, y, Fx, Fy, "black"); // forces are in color black
            }
        }
    }
    frame.strokeColor("black");
    std::string border_info = "Time: " + std::to_string(t) + 
                              "\nRange: [" + std::to_string(min_x) + ", " + std::to_string(max_x) + "] x " +
                              "[" + std::to_string(min_y) + ", " + std::to_string(max_y) + "]";
    frame.annotate(border_info, Magick::NorthWestGravity);
    frames.push_back(frame);
}

void visualize(const std::vector<std::vector<Vector2D>>& all_positions, const std::vector<std::vector<Vector2D>>& all_velocities, const std::vector<std::vector<Vector2D>>& all_forces, int n, double time_step, double total_time) {
    Magick::InitializeMagick(nullptr);
    std::vector<Magick::Image> frames;
    double min_x = all_positions[0][0].x;
    double max_x = all_positions[0][0].x;
    double min_y = all_positions[0][0].y;
    double max_y = all_positions[0][0].y;
    for (const auto& positions_at_time : all_positions) {
        for (const auto& pos : positions_at_time) {
            if (pos.x < min_x) min_x = pos.x;
            if (pos.x > max_x) max_x = pos.x;
            if (pos.y < min_y) min_y = pos.y;
            if (pos.y > max_y) max_y = pos.y;
        }
    }

    double margin_x = (max_x - min_x) * 0.1;
    double margin_y = (max_y - min_y) * 0.1;
    min_x -= margin_x;
    max_x += margin_x;
    min_y -= margin_y;
    max_y += margin_y;

    for (double t = 0; t < total_time; t += time_step) {
        save_frame(all_positions[t], all_velocities[t], all_forces[t], n, t*time_step, frames, min_x, max_x, min_y, max_y);
    }

    Magick::writeImages(frames.begin(), frames.end(), "nbody_simulation.gif");
}

int main() {
    int n;
    std::vector<double> masses;
    std::vector<Vector2D> positions, velocities;
    double time_step, total_time;

    gather_input(n, masses, positions, velocities, time_step, total_time);

    std::vector<Vector2D> forces(n);
    
    std::vector<std::vector<Vector2D>> all_positions, all_velocities, all_forces;
    all_positions.push_back(positions);
    all_velocities.push_back(velocities);

    for (double t = 0; t < total_time; t += time_step) {
        compute_forces(n, masses, positions, forces);
        update_bodies(n, masses, positions, velocities, forces, time_step);

        all_positions.push_back(positions);
        all_velocities.push_back(velocities);
        all_forces.push_back(forces);
                std::cout << "Time: " << t + time_step << std::endl;
        for (int i = 0; i < n; ++i) {
            std::cout << "Body " << i + 1 << ": Position (" << positions[i].x << ", " << positions[i].y << ")" << std::endl;
        }
    }    

    visualize(all_positions, all_velocities, all_forces, n, time_step, total_time);

    return 0;
}
