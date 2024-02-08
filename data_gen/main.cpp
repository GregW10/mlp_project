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

std::mutex mutex;
ull_t evols;

unsigned int running = 0;

void run_sim(gtd::sys sys, ull_t epochs) {

}

int main(int argc, char **argv) {
    gtd::parser parser{argc, argv}; // I create an instance of my `gtd::parser` class to allow easy extraction of args
    evols = parser.get_arg("--evolutions", 100'000ull); // number of evolutions to simulate
    ull_t epochs = parser.get_arg("--epochs", 10'000ull); // number of simulation epochs to write per evolution
    ull_t iters = parser.get_arg("--iterations", 1'000'000'000ull); // number of iterations per epoch
    ull_t num_m = parser.get_arg("--num_mass_means", 100ull); // number of different mean masses to use
    long double start_mass = parser.get_arg("--starting_mass", 1'000.0l); // starting mean mass value
    long double mass_step = parser.get_arg("--mass_step", TRILLION/num_m); // step in mean mass value
    // Despite not having collisions activated for the training data, I will still give the particles a non-zero radius,
    // in case I wish to render certain evolutions with my ray-tracer later:
    long double radius = parser.get_arg("--radius", 0.0625l);
    long double dt = parser.get_arg("--timestep", 1.0l/pow(2, 28)); // time-step used to evolve simulations
    long double box_x = parser.get_arg("--box_width", 2.0l); // length of bounding box along x
    long double box_y = parser.get_arg("--box_length", 2.0l); // length of bounding box along y
    long double box_z = parser.get_arg("--box_height", 2.0l); // length of bounding box along z
    // Note: the box is always centered around the COM, so centre coordinates are never supplied as an argument
    long double eps = parser.get_arg("--softening", 1.0l/65'536.0l); // softening length
    // The "minimum separation" is the minimum starting separation the bodies have to have:
    long double min_sep_def = eps*10; // the default is 10 softening lengths - if the bodies are closer, pos. resampled
    long double min_sep = parser.get_arg("--min_sep", min_sep_def);
    // This next parameter is important, as it will determine whether the evolution of the 3-body system should stop if
    // one of the 3 bodies ends up being flung far out:
    bool e_stop = parser.get_arg("-e", false) | parser.get_arg("--early_stop", false); // bitwise given short-circuit
    // Note that it could only ever be one body being flung far out, as all bodies getting flung away from each other
    // would violate the conservation of energy (assuming they all start out with velocities below their mutual EVs).
    bool verbose = parser.get_arg("-v", false) | parser.get_arg("--verbose", false);
    if (!parser.empty()) {
        std::cerr << BOLD_TXT(RED_TXT("Error!")) YELLOW_TXT(" Unrecognised arguments:") "\n" << BLUE_TXT_START;
        for (const auto &arg : parser)
            std::cout << arg << '\n';
        std::cout << RESET_TXT_FLAGS;
        return 1;
    }
    if (!evols) {
        std::cerr << "Error: cannot have zero evolutions.\n";
        return 1;
    }
    if (!epochs) {
        std::cerr << "Error: cannot have zero epochs per evolution.\n";
        return 1;
    }
    if (!iters) {
        std::cerr << "Error: cannot have zero iterations per epoch.\n";
        return 1;
    }
    if (!num_m) {
        std::cerr << "Error: cannot have zero mass means.\n";
        return 1;
    }
    if (num_m > evols) {
        std::cerr << "Error: cannot have more mass means than number of evolutions.\n";
        return 1;
    }
    if (start_mass <= 0) {
        std::cerr << "Error: starting mass must be positive.\n";
        return 1;
    }
    if (mass_step <= 0) {
        std::cerr << "Error: mass step must be positive.\n";
        return 1;
    }
    if (radius < 0) {
        std::cerr << "Error: radius cannot be negative.\n";
        return 1;
    }
    if (dt <= 0) {
        std::cerr << "Error: time-step must be positive.\n";
        return 1;
    }
    if (box_x <= 0) {
        std::cerr << "Error: box width must be positive.\n";
        return 1;
    }
    if (box_y <= 0) {
        std::cerr << "Error: box length must be positive.\n";
        return 1;
    }
    if (box_z <= 0) {
        std::cerr << "Error: box height must be positive.\n";
        return 1;
    }
    if (eps < 0) {
        std::cerr << "Error: softening-length must be non-negative.\n";
        return 1;
    }
    if (min_sep > std::min({box_x, box_y, box_z})/2.0l) {
        std::cerr << "Error: the minimum separation between the bodies cannot be greater than half the minimum side "
                     "length of the bounding box.\n";
        return 1;
    }
    char *path = new char[41];
    gtd::strcpy_c(path, "experiment_");
    gtd::strcat_c(path, gtd::get_date_and_time());
    if (verbose)
        std::cout << YELLOW_TXT("Attempting to create \"") BLUE_TXT_START << path <<
        YELLOW_TXT("\" directory...") << std::endl;
    if (mkdir(path, S_IRWXU | S_IRWXO | S_IRWXG) == -1) {
        std::cerr << BOLD_TXT(RED_TXT("Error: ")) YELLOW_TXT("could not create \"") BLUE_TXT_START << path <<
        YELLOW_TXT("\" directory.") << std::endl;
        delete [] path;
        return 1;
    }
    if (verbose)
        std::cout << YELLOW_TXT("Created \"") BLUE_TXT_START << path << YELLOW_TXT("\" directory.") << std::endl;
    if (chdir(path) == -1) {
        std::cerr << BOLD_TXT(RED_TXT("Error: ")) YELLOW_TXT("could change working directory to \"") BLUE_TXT_START
        << path << YELLOW_TXT("\".") << std::endl;
        delete [] path;
        return 1;
    }
    gtd::strcat_c(path, ".txt");
    if (verbose)
        std::cout << YELLOW_TXT("Attempting to create \"") BLUE_TXT_START << path <<
              YELLOW_TXT("\" file...") << std::endl;
    std::ofstream logf{path, std::ios_base::out | std::ios_base::trunc};
    if (!logf.good()) {
        std::cerr << BOLD_TXT(RED_TXT("Error: ")) YELLOW_TXT("could not create \"") BLUE_TXT_START << path <<
                  YELLOW_TXT("\" file.") << std::endl;
        delete [] path;
        return 1;
    }
    if (verbose)
        std::cout << YELLOW_TXT("Writing configurations to \"") BLUE_TXT_START << path <<
                  YELLOW_TXT("\" file...") << std::endl;
    logf << "Distinct evolutions = " << evols << "\nEpochs per evolution = " << epochs << "\nIterations per epoch = "
         << iters << "\nNumber of mass means = " << num_m << "\nStarting mean mass = " << start_mass
         << " kg\nStep in mean mass = " << mass_step << " kg\nRadius of bodies = " << radius << " m\nTime-step = " << dt
         << " s\nBox length along x = " << box_x << " m\nBox length along y = " << box_y << " m\nBox length along z = "
         << box_z << " m\nSoftening length = " << eps << " m\nDefault min. sep. = " << min_sep_def
         << " m\nActual min. sep. = " << min_sep << " m\nEarly stopping = " << std::boolalpha << e_stop
         << "\nVerbose output = " << verbose << '\n';
    if (verbose)
        std::cout << YELLOW_TXT("Configurations written to \"") BLUE_TXT_START << path <<
                  YELLOW_TXT("\" file.") << std::endl;
    logf.close();
    delete [] path;
    unsigned int numt = std::thread::hardware_concurrency();
    numt = numt ? numt : 1; // single-thread execution if number could not be determined
    long double max_x = box_x/2.0l;
    long double min_x = -max_x;
    long double max_y = box_y/2.0l;
    long double min_y = -max_y;
    long double max_z = box_z/2.0l;
    long double min_z = -max_z;
    gtd::sys sys;
    sys.set_timestep(dt);
    sys.set_iterations(iters);
    sys.emplace_body(false); sys.emplace_body(false); sys.emplace_body(false);
    std::mt19937_64 rng{std::random_device{}()}; // Mersenne-Twister engine
    std::uniform_real_distribution<long double> pdist_x{min_x, max_x}; // used for position value x-coordinate
    std::uniform_real_distribution<long double> pdist_y{min_y, max_y}; // used for position value y-coordinate
    std::uniform_real_distribution<long double> pdist_z{min_z, max_z}; // used for position value z-coordinate
    std::uniform_real_distribution<long double> vdist; // used for velocity values
    std::lognormal_distribution<long double> mdist; // used for mass values
    while (evols > 0) {
        do {
            sys.front().pos()[0] = pdist_x(rng); // this is more efficient than using a loop, despite more clutter
            sys.front().pos()[1] = pdist_y(rng);
            sys.front().pos()[2] = pdist_z(rng);
            sys[1].pos()[0] = pdist_x(rng);
            sys[1].pos()[1] = pdist_y(rng);
            sys[1].pos()[2] = pdist_z(rng);
            if (gtd::vec_ops::distance(sys.front().pos(), sys[1].pos()) < min_sep)//check here to avoid redundant below
                continue;
            sys.back().pos()[0] = pdist_x(rng);
            sys.back().pos()[1] = pdist_y(rng);
            sys.back().pos()[2] = pdist_z(rng);
        } while (gtd::vec_ops::distance(sys.front().pos(), sys[2].pos()) < min_sep ||
                 gtd::vec_ops::distance(sys[1].pos(), sys.back().pos()) < min_sep);
    }
    return 0;
}
