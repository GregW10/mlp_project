#include "datsup.hpp"

/* Acronyms:
 * - COM: Centre of Mass
 * - COMV: Centre of Mass velocity
 * - EV: Escape Velocity
 * */

// generate 3 random velocities
// make sure no body exceeds escape velocity (EV) - tricky because others might move away, maybe have threshold under EV
// subtract COM velocity from all velocities

typedef unsigned long long ull_t;

int main(int argc, char **argv) {
    gtd::parser parser{argc, argv}; // I create an instance of my `gtd::parser` class to allow easy extraction of args
    ull_t evols = parser.get_arg("--evolutions", 100'000ull); // number of evolutions to simulate
    ull_t epochs = parser.get_arg("--epochs", 10'000ull); // number of simulation epochs to write per evolution
    ull_t iters = parser.get_arg("--iterations", 1'000'000'000ull); // number of iterations per epoch
    ull_t num_m = parser.get_arg("--num_mass_means", evols/100); // number of different mean masses to use
    long double dt = parser.get_arg("--timestep", 1.0l/pow(2, 28)); // time-step used to evolve simulations
    long double box_x = parser.get_arg("--box_width", 2.0l); // length of bounding box along x
    long double box_y = parser.get_arg("--box_length", 2.0l); // length of bounding box along y
    long double box_z = parser.get_arg("--box_height", 2.0l); // length of bounding box along z
    // Note: the box is always centered around the COM, so centre coordinates are never supplied as an argument
    // This next parameter is important, as it will determine whether the evolution of the 3-body system should stop if
    // one of the 3 bodies ends up being flung far out:
    bool e_stop = parser.get_arg("-e", false);
    // Note that it could only ever be one body being flung far out, as all bodies getting flung away from each other
    // would violate the conservation of energy (assuming they all start out with velocities below their mutual EVs).
    bool verbose = parser.get_arg("-v", false);
    if (!parser.empty()) {
        std::cerr << "Unrecognised arguments:\n";
        for (const auto &arg : parser)
            std::cout << arg << '\n';
        return 1;
    }
    if (!evols) {
        // error
    }
    if (!epochs) {
        // error
    }
    if (!iters) {
        // error
    }
    if (!num_m) {
        // error
    }
    if (num_m > evols) {
        // error
    }
    if (dt <= 0) {
        // error
    }
    if (box_x <= 0) {
        // error
    }
    if (box_y <= 0) {
        // error
    }
    if (box_z <= 0) {
        // error
    }
    return 0;
}
