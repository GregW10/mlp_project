#include "datsup.hpp"

/* Acronyms:
 * - COM: Centre of Mass
 * - COMV: Centre of Mass Velocity
 * - EV: Escape Velocity
 * - KE: Kinetic Energy
 * - PE: Potential Energy
 * - TE: Total Energy
 * - SD: Standard Deviation
 * - LN: Log-Normal
 * */

using ptype_ln = std::lognormal_distribution<long double>::param_type;

std::mutex main_mutex;
std::mutex worker_mutex;
std::condition_variable cv;
unsigned int running = 0; // to log the rumber of worker threads running at one time

uint64_t evols;
uint64_t epochs;

bool verbose;
bool nsys;

void log_ejection(const gtd::sys &sys, char *path, uint64_t zero_index, uint64_t counter) {
    std::lock_guard<std::mutex> wguard{worker_mutex};
    *(path + zero_index) = 0;
    std::cerr << "----------------\nBody ejected in " << path << " simulation at epoch " << counter
              << '/' << epochs << ".\nBodies when ejection confirmed:\n";
    counter = 1;
    for (const gtd::bod_0f &b : sys)
        std::cerr << "------\nBody " << counter++ << "\n------\nPosition: " << b.pos()
                  << " m\nVelocity: " << b.vel() << " m/s\nAcceleration: " << b.acceleration()
                  << " m/s^2\n";
    std::cerr << "----------------\n";
    --running;
    cv.notify_one();
}

template <bool check_ejection>
void run_sim(gtd::sys sys, uint64_t batch_num, uint64_t num_in_batch) {
    if (sys.num_bodies() != 3)
        return;
    gtd::String path{true};
    path.append_back(batch_num - 1); // to make it start at 0
    path.push_back('_');
    path.append_back(num_in_batch);
    if (mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IRWXO) == -1) {
        std::lock_guard<std::mutex> guard{worker_mutex};
        std::cerr << "Error: could not create directory \"" << path << "\".\n";
        --running;
        return;
        // exit(1); // not the best idea... will find a more graceful way to deal with this (reduce evols and return?)
    }
    path.append_back("/init_params.csv");
    std::ofstream params{path.c_str(), std::ios_base::out | std::ios_base::trunc};
    if (!params.good())
        goto rest; // if for some reason the file cannot be created, it's not the end of the world
    params << "bod1_pos_x,bod1_pos_y,bod1_pos_z,bod1_vel_x,bod1_vel_y,bod1_vel_z,bod1_mass,"
              "bod2_pos_x,bod2_pos_y,bod2_pos_z,bod2_vel_x,bod2_vel_y,bod2_vel_z,bod2_mass,"
              "bod3_pos_x,bod3_pos_y,bod3_pos_z,bod3_vel_x,bod3_vel_y,bod3_vel_z,bod3_mass\r\n";
    for (const auto &b : sys)
        params << b.pos()[0] << ',' << b.pos()[1] << ',' << b.pos()[2] << ','
               << b.vel()[0] << ',' << b.vel()[1] << ',' << b.vel()[2] << ','
               << b.mass() << ',';
    params.seekp(-1, std::ios_base::cur);
    params << "\r\n";
    params.close();
    rest:
    path.erase_chars(path.get_length() - 15);
    path.append_back("log.3bod");
    gtd::f3bod<long double> file{path.c_str(), &sys, epochs + 1};
    file.add_entry();
    uint64_t counter = 1;
    if (nsys) {
        path.erase_chars(path.get_length() - 8);
        path.append_back("epoch");
        unsigned long long max_len = llroundl(floorl(log10l(epochs)) + 1);
        uint64_t epoch_len = path.get_length();
        char *npath = new char[epoch_len + max_len + 6]; // "epochXX...XX.nsys"
        gtd::strcpy_c(npath, path.c_str());
        path.clear();
        char *last_digit = npath + epoch_len + max_len - 1;
        char *ptr = npath + epoch_len;
        while (counter++ <= max_len)
            *ptr++ = '0';
        gtd::strcpy_c(ptr, ".nsys");
        sys.to_nsys(npath);
        counter = 1;
        uint64_t val;
        while (counter <= epochs) {
            sys.evolve();
            file.add_entry();
            val = counter++;
            ptr = last_digit;
            while (val > 0) {
                *ptr-- = (val % 10) + 48;
                val /= 10;
            }
            sys.to_nsys(npath);
            if constexpr (check_ejection) { // the constexpr statement avoids constant re-checking within the loop
                // check if a body's velocity is greater than its EV w.r.t. the other two
                if (gtd::will_eject(sys)) {
                    log_ejection(sys, npath, epoch_len - 6, counter);
                    delete [] npath;
                    return;
                }
            }
        }
        delete [] npath;
    }
    else {
        while (counter++ <= epochs) {
            sys.evolve();
            file.add_entry();
            if constexpr (check_ejection) {
                if (gtd::will_eject(sys)) {
                    // using const_cast<> ain't pretty, change this:
                    log_ejection(sys, const_cast<char*>(path.c_str()), path.get_length() - 9, counter);
                    return;
                }
            }
        }
    }
    std::lock_guard<std::mutex> guard{worker_mutex}; // acquire worker mutex to stop other workers from red. `running`
    --running; // reduce the counter for the number of threads running
    cv.notify_one(); // notify `main` that the worker thread finished, `main` either creates another thread or exits
}

int main(int argc, char **argv) {
    gtd::parser parser{argc, argv}; // I create an instance of my `gtd::parser` class to allow easy extraction of args
    evols = parser.get_arg("--evolutions", 100'000ull); // number of evolutions to simulate
    bool tty = isatty(STDOUT_FILENO); // if connected to a terminal, output colours, if not, don't
    if (!evols) {
        std::cerr << "Error: cannot have zero evolutions.\n";
        return 1;
    }
    epochs = parser.get_arg("--epochs", 10'000ull); // number of simulation epochs to write per evolution
    uint64_t iters = parser.get_arg("--iterations", 1'000'000ull); // number of iterations per epoch
    // Number of different mean masses to use:
    uint64_t num_m = parser.get_arg("--num_mass_means", evols /* (uint64_t) ceill(((long double) evols)/100.0l) */);
    // `num_m` has to be checked up here, as the default value for `mass_step` uses it
    if (!num_m) {
        std::cerr << "Error: cannot have zero mass means.\n";
        return 1;
    }
    if (num_m > evols) {
        std::cerr << "Error: cannot have more mass means than number of evolutions.\n";
        return 1;
    }
    const long double start_mass = parser.get_arg("--starting_mass", 1'000.0l); // starting mean mass value
    const long double mass_step = parser.get_arg("--mass_step", TRILLION/(num_m - 1)); // step in mean mass value
    const long double mass_sd_scaling = parser.get_arg("--mass_sd_scaling", 0.125l); // mass_SD = scaling*mean_mass
    // I give the particles a radius of zero, as I will be training the NN to predict the collisionless 3-body problem,
    // and so the particles are treated as point particles:
    const long double radius = parser.get_arg("--radius", 0.0l);
    const long double dt = parser.get_arg("--timestep", 1.0l/pow(2, 16)); // time-step used to evolve simulations
    const long double vel_scale = parser.get_arg("--velocity_scaling", 0.5l); // scaling factor for velocities
    const long double box_x = parser.get_arg("--box_width", 2.0l); // length of bounding box along x
    const long double box_y = parser.get_arg("--box_length", 2.0l); // length of bounding box along y
    const long double box_z = parser.get_arg("--box_height", 2.0l); // length of bounding box along z
    // Note: the box is always centered around the COM, so centre coordinates are never supplied as an argument
    const long double eps = parser.get_arg("--softening", 1.0l/65'536.0l); // softening length
    // The "minimum separation" is the minimum starting separation the bodies have to have:
    long double min_sep_def = eps*10; // the default is 10 softening lengths - if the bodies are closer, pos. resampled
    long double min_sep = parser.get_arg("--min_sep", min_sep_def);
    if (2*radius >= min_sep)
        std::cerr << BOLD_TXT_START MAGENTA_TXT("Warning:")
        YELLOW_TXT(" there is a non-zero probability that the particles generated with the given radius of ")
        BLUE_TXT_START
                  << radius << RESET_TXT_FLAGS BOLD_TXT_START GREEN_TXT(" m ") YELLOW_TXT_START
                  "will overlap when positions are randomly generated, which would cause termination of the "
                               "simulation.\nConsider choosing a smaller radius.\n" << RESET_TXT_FLAGS;
    // This next parameter is important, as it will determine whether the evolution of the 3-body system should stop if
    // one of the 3 bodies ends up being flung far out:
    bool e_stop = parser.get_arg("-e", false) | parser.get_arg("--early_stop", false); // bitwise given short-circuit
    // Note that it could only ever be one body being flung far out, as all bodies getting flung away from each other
    // would violate the conservation of energy (assuming they all start out with tot_KE < tot_PE).
    verbose = parser.get_arg("-v", false) | parser.get_arg("--verbose", false); // for verbose output
    nsys = parser.get_arg("-n", false) | parser.get_arg("--output_nsys", false); // for nsys files too
    if (!parser.empty()) {
        std::cerr << BOLD_TXT(RED_TXT("Error!")) YELLOW_TXT(" Unrecognised arguments:") "\n" << BLUE_TXT_START;
        for (const std::string &arg : parser)
            std::cout << arg << '\n';
        std::cout << RESET_TXT_FLAGS;
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
    if (start_mass <= 0) {
        std::cerr << "Error: starting mass must be positive.\n";
        return 1;
    }
    if (mass_step <= 0) {
        std::cerr << "Error: mass step must be positive.\n";
        return 1;
    }
    if (mass_sd_scaling <= 0) {
        std::cerr << "Error: the scaling factor applied to the mass to produce the S.D. must be positive.\n";
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
    if (vel_scale <= 0 || vel_scale > 1) {
        std::cerr << "Error: the velocity scaling factor must be within the (0,1] interval.\n";
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
    unsigned int numt = std::thread::hardware_concurrency();
    numt = numt ? numt : 1; // single-thread execution if number could not be determined
    logf << "Distinct evolutions = " << evols << "\nEpochs per evolution = " << epochs << "\nIterations per epoch = "
         << iters << "\nNumber of mass means = " << num_m << "\nStarting mean mass = " << start_mass
         << " kg\nStep in mean mass = " << mass_step << " kg\nSD fraction of mean = " << mass_sd_scaling
         << "\nRadius of bodies = " << radius << " m\nTime-step = " << dt
         << " s\nVelocity scaling factor = " << vel_scale
         << "\nBox length along x = " << box_x << " m\nBox length along y = " << box_y << " m\nBox length along z = "
         << box_z << " m\nSoftening length = " << eps << " m\nDefault min. sep. = " << min_sep_def
         << " m\nActual min. sep. = " << min_sep << " m\nEarly stopping = " << std::boolalpha << e_stop
         << "\nVerbose output = " << verbose << "\nOutput nsys = " << nsys << "\nNumber of parallel simulations = "
         << numt << '\n';
    if (verbose)
        std::cout << YELLOW_TXT("Configurations written to \"") BLUE_TXT_START << path <<
                  YELLOW_TXT("\" file.") << std::endl;
    logf.close();
    delete [] path;
    long double max_x = box_x/2.0l;
    long double min_x = -max_x;
    long double max_y = box_y/2.0l;
    long double min_y = -max_y;
    long double max_z = box_z/2.0l;
    long double min_z = -max_z;
    gtd::sys sys;
    sys.timestep(dt);
    sys.iters(iters);
    sys.softening(eps);
    sys.emplace_body(false).emplace_body(false).emplace_body(false);
    for (gtd::bod_0f &b : sys)
        b.set_radius(radius);
    gtd::bod_0f *bod1 = &sys.front();
    gtd::bod_0f *bod2 = &sys[1];
    gtd::bod_0f *bod3 = &sys.back();
    std::mt19937_64 rng{std::random_device{}()}; // Mersenne-Twister engine
    // Isotropic position and velocity distributions:
    std::uniform_real_distribution<long double> pdist_x{min_x, max_x}; // used for position value x-coordinate
    std::uniform_real_distribution<long double> pdist_y{min_y, max_y}; // used for position value y-coordinate
    std::uniform_real_distribution<long double> pdist_z{min_z, max_z}; // used for position value z-coordinate
    std::uniform_real_distribution<long double> vdist{-1, 1}; // used for velocity values
    std::lognormal_distribution<long double> mdist; // used for mass values
    long double final_scaling;
    gtd::vec3 com_pos;
    gtd::vec3 com_vel;
    const long double mass_batch_size = ((long double) evols)/num_m;
    uint64_t mass_counter = 1; // number of distinct mass values used
    uint64_t next_mass_switch = llroundl(mass_batch_size);
    uint64_t num_in_batch = 0;
    long double mass_mean = start_mass;
    long double mass_sd = mass_sd_scaling*mass_mean;
    if (verbose)
        std::cout << "Starting mass batch with mean mass " << mass_mean << " kg and SD " << mass_sd << " kg"
                  << std::endl;
    std::pair<long double, long double> mu_sd_gauss;
    uint64_t counter = 0;
    std::unique_lock<std::mutex> mlock{main_mutex};//, std::defer_lock};
    std::unique_lock<std::mutex> wlock{worker_mutex, std::defer_lock};
    while (/* evols > 0 */ counter < evols) {
        // First I set the masses (as they are required for COM calculations):
        if (counter == next_mass_switch) {
            mass_mean = start_mass + mass_counter*mass_step; // this approach keeps higher precision than rep. add.
            mass_sd = mass_sd_scaling*mass_mean;
            next_mass_switch = llroundl(++mass_counter*mass_batch_size);
            num_in_batch = 0;
            if (verbose)
                std::cout << "Starting mass batch with mean mass " << mass_mean << " kg and SD " << mass_sd << " kg"
                          << std::endl;
        }
        mu_sd_gauss = gtd::lognormal_to_gauss(mass_mean, mass_sd); // convert LN mean and SD into Gaussian mean and SD
        // The Gaussian mean and SD are the parameters required for `std::lognormal_distribution`:
        mdist.param(ptype_ln{mu_sd_gauss.first, mu_sd_gauss.second});
        for (gtd::bod_0f &bod : sys)
            bod.set_mass(mdist(rng));
        // Now I generate initial random positions with components in [-1,1]:
        do {
            bod1->pos()[0] = pdist_x(rng); // this is more efficient than using a loop, despite more clutter
            bod1->pos()[1] = pdist_y(rng);
            bod1->pos()[2] = pdist_z(rng);
            bod2->pos()[0] = pdist_x(rng);
            bod2->pos()[1] = pdist_y(rng);
            bod2->pos()[2] = pdist_z(rng);
            if (gtd::vec_ops::distance(bod1->pos(), bod2->pos()) < min_sep)//check here to avoid redundancy below
                continue;
            bod3->pos()[0] = pdist_x(rng);
            bod3->pos()[1] = pdist_y(rng);
            bod3->pos()[2] = pdist_z(rng);
        } while (gtd::vec_ops::distance(bod1->pos(), bod3->pos()) < min_sep ||
                 gtd::vec_ops::distance(bod2->pos(), bod3->pos()) < min_sep);
        // Now I centre the origin of the coordinate system at the COM:
        com_pos = sys.com_pos();
        for (gtd::bod_0f &b : sys)
            b.pos() -= com_pos;
        // Finally, I set the velocities. First I generate their x, y and z components in [-1,1]:
        for (gtd::bod_0f &bod : sys) {
            bod.vel()[0] = vdist(rng);
            bod.vel()[1] = vdist(rng);
            bod.vel()[2] = vdist(rng);
        }
        // Next, I subtract the COM velocity from all bodies to ensure we are in the COM frame:
        com_vel = sys.com_vel();
        for (gtd::bod_0f &bod : sys)
            bod.vel() -= com_vel; // velocities have to be in the COM frame!!!
        // I calculate the final scaling factor for velocities and multiply this by all velocities:
        final_scaling = vel_scale/sqrtl(gtd::beta_factor(sys)); // have decided NOT to take sqrt of velocity scaling
        for (gtd::bod_0f &bod : sys)
            bod.vel() *= final_scaling;
        if (running == numt)
            cv.wait(mlock, [&numt](){return running < numt;});
        wlock.lock();
        ++running;
        wlock.unlock();
        if (e_stop)
            // create worker thread and immediately detach:
            std::thread{run_sim<true>, sys, mass_counter, num_in_batch}.detach();
        else
            std::thread{run_sim<false>, sys, mass_counter, num_in_batch}.detach();
        ++counter;
        ++num_in_batch;
        if (verbose)
            std::cout << "Evolution " << counter << '/' << evols << " dispatched." << std::endl;
    }
    cv.wait(mlock, [](){return !running;}); // wait for all detached worker threads to have finished running
    std::cout << "------------------------\nAll evolutions finished.\n------------------------" << std::endl;
    return 0;
}
