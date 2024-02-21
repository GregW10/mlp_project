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
long double ejf;
long double dts;
long double dte;

// bool verbose;
bool nsys;

uint64_t ejected = 0; // number of bodies ejected - only used if ejection detection is enabled

void log_ejection(const gtd::sys &sys, char *path, uint64_t zero_index, uint64_t counter) {
    std::lock_guard<std::mutex> wguard{worker_mutex};
    *(path + zero_index) = 0;
    std::cerr << "----------------\nBody ejected in simulation " << path << " at epoch " << counter
              << '/' << epochs << ".\nBodies when ejection confirmed:\n";
    counter = 1;
    for (const gtd::bod_0f &b : sys)
        std::cerr << "------\nBody " << counter++ << "\n------\nPosition: " << b.pos()
                  << " m\nVelocity: " << b.vel() << " m/s\nAcceleration: " << b.acceleration()
                  << " m/s^2\n";
    std::cerr << "----------------\n";
    --running;
    ++ejected;
    cv.notify_one();
}

template <bool check_ejection, bool tscale>
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
    params << std::setprecision(20);
    if constexpr (tscale)
        params << "bod1_pos_x,bod1_pos_y,bod1_pos_z,bod1_vel_x,bod1_vel_y,bod1_vel_z,bod1_mass,"
                  "bod2_pos_x,bod2_pos_y,bod2_pos_z,bod2_vel_x,bod2_vel_y,bod2_vel_z,bod2_mass,"
                  "bod3_pos_x,bod3_pos_y,bod3_pos_z,bod3_vel_x,bod3_vel_y,bod3_vel_z,bod3_mass,"
                  "timestep,iterations\r\n";
    else
        params << "bod1_pos_x,bod1_pos_y,bod1_pos_z,bod1_vel_x,bod1_vel_y,bod1_vel_z,bod1_mass,"
                  "bod2_pos_x,bod2_pos_y,bod2_pos_z,bod2_vel_x,bod2_vel_y,bod2_vel_z,bod2_mass,"
                  "bod3_pos_x,bod3_pos_y,bod3_pos_z,bod3_vel_x,bod3_vel_y,bod3_vel_z,bod3_mass\r\n";
    for (const auto &b : sys)
        params << b.pos()[0] << ',' << b.pos()[1] << ',' << b.pos()[2] << ','
               << b.vel()[0] << ',' << b.vel()[1] << ',' << b.vel()[2] << ','
               << b.mass() << ',';
    if constexpr (tscale)
        params << sys.timestep() << ',' << sys.iters() << "\r\n";
    else {
        params.seekp(-1, std::ios_base::cur);
        params << "\r\n";
    }
    params.close();
    rest:
    path.erase_chars(path.get_length() - 15);
    path.append_back("log.3bod");
    gtd::f3bodw<long double> file{path.c_str(), &sys, epochs + 1};
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
                if (gtd::will_eject(sys, ejf)) {
                    log_ejection(sys, npath, epoch_len - 6, counter);
                    delete [] npath;
                    file.remaining(0); // have to change the number of entries in the .3bod file, as evol. is cut short
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
                if (gtd::will_eject(sys, ejf)) {
                    // using const_cast<> ain't pretty, change this:
                    log_ejection(sys, const_cast<char*>(path.c_str()), path.get_length() - 9, counter);
                    file.remaining(0);
                    return;
                }
            }
        }
    }
    std::lock_guard<std::mutex> guard{worker_mutex}; // acquire worker mutex to stop other workers from red. `running`
    --running; // reduce the counter for the number of threads running
    cv.notify_one(); // notify `main` that the worker thread finished, `main` either creates another thread or exits
}

template <bool e_stop, bool verbose, bool tty, bool tscale>
void main_loop(long double box_x, long double box_y, long double box_z, long double dt, uint64_t iters, long double eps,
               long double radius, uint64_t num_m, long double start_mass, long double mass_sd_scaling,
               long double mass_step, long double min_sep, long double vel_scale, unsigned int numt) {
    long double max_x = box_x/2.0l;
    long double min_x = -max_x;
    long double max_y = box_y/2.0l;
    long double min_y = -max_y;
    long double max_z = box_z/2.0l;
    long double min_z = -max_z;
    gtd::sys sys;
    if constexpr (!tscale) { // no scaling, so fixed time-step and number of iterations across all mass batches
        sys.timestep(dt);
        sys.iters(iters);
    }
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
    long double triple_mass = 3*start_mass; // only needed if `tscale == true`
    if constexpr (verbose) {
        if constexpr (tty)
            std::cout << "\033[38;5;206mStarting mass batch with mean mass \033[1m\033[38;5;214m" << mass_mean <<
            " \033[1m\033[38;5;11mkg \033[0m\033[38;5;206mand SD \033[1m\033[38;5;214m" << mass_sd <<
            " \033[1m\033[38;5;11mkg.\033[0m" << std::endl;
        else
            std::cout << "Starting mass batch with mean mass " << mass_mean << " kg and SD " << mass_sd << " kg"
                      << std::endl;
    }
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
            if constexpr (verbose) {
                if constexpr (tty)
                    std::cout << "\033[38;5;206mStarting mass batch with mean mass \033[1m\033[38;5;214m" << mass_mean
                              << " \033[1m\033[38;5;11mkg \033[0m\033[38;5;206mand SD \033[1m\033[38;5;214m" << mass_sd
                              << " \033[1m\033[38;5;11mkg.\033[0m" << std::endl;
                else
                    std::cout << "Starting mass batch with mean mass " << mass_mean << " kg and SD " << mass_sd << " kg"
                              << std::endl;
            }
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
        if constexpr (tscale) {
            long double factor = dts*powl(triple_mass/(sys[0].mass() + sys[1].mass() + sys[2].mass()), dte);
            sys.timestep(factor*dt);
            sys.iters((uint64_t) ceill(((long double) iters)/factor));
        }
        if (running == numt)
            cv.wait(mlock, [&numt](){return running < numt;});
        wlock.lock();
        if constexpr (verbose) {
            if constexpr (tscale) {
                if constexpr (tty) {
                    std::cout << "\033[1m\033[38;5;123m----------------------------\n"
                                 "\033[0m\033[38;5;10mStarting new evolution with:\n"
                                 "\033[38;5;13mTime-step \033[38;5;11m= \033[1m\033[38;5;14m"
                                 << sys.timestep() << " \033[0m\033[38;5;15ms\n"
                                 "\033[38;5;13mIterations \033[38;5;11m= \033[1m\033[38;5;14m"
                                 << sys.iters() << "\n"
                                 "\033[38;5;13mTotal sim. time \033[38;5;11m= \033[1m\033[38;5;14m"
                                 << sys.timestep()*sys.iters()*epochs
                                 << " \033[0m\033[38;5;15ms\n\033[1m\033[38;5;123m----------------------------\033[0m"
                                 << std::endl;
                } else {
                    std::cout << "----------------------------\nStarting new evolution with:\nTime-step = "
                              << sys.timestep() << " s\nIterations = " << sys.iters() << "\nTotal sim. time = "
                              << sys.timestep()*sys.iters()*epochs << " s\n----------------------------" << std::endl;
                }
            }
        }
        ++running;
        wlock.unlock();
        // create worker thread and immediately detach:
        std::thread{run_sim<e_stop, tscale>, sys, mass_counter, num_in_batch}.detach();
        ++counter;
        ++num_in_batch;
        if constexpr (verbose) {
            if constexpr (tty)
                std::cout << "\033[38;5;117mEvolution \033[1m\033[38;5;10m" << counter << "\033[38;5;11m/\033[38;5;10m"
                          << evols << " \033[0m\033[38;5;117mdispatched.\033[0m" << std::endl;
            else
                std::cout << "Evolution " << counter << '/' << evols << " dispatched." << std::endl;
        }
    }
    cv.wait(mlock, [](){return !running;}); // wait for all detached worker threads to have finished running
}

[[noreturn]] void log_error(const char *str, bool istty) {
    if (!str)
        exit(1);
    if (istty)
        std::cerr << BOLD_TXT_START RED_TXT_START "Error: " RESET_TXT_FLAGS << YELLOW_TXT_START << str << '\n';
    else
        std::cerr << "Error: " << str << '\n';
    exit(1);
}

int main(int argc, char **argv) {
    gtd::parser parser{argc, argv}; // I create an instance of my `gtd::parser` class to allow easy extraction of args
    // This next parameter is important, as it will determine whether the evolution of the 3-body system should stop if
    // one of the 3 bodies ends up being flung far out:
    bool e_stop = parser.get_arg("-e", false) | parser.get_arg("--early_stop", false); // bitwise given short-circuit
    // Note that it could only ever be one body being flung far out, as all bodies getting flung away from each other
    // would violate the conservation of energy (assuming they all start out with tot_KE < tot_PE).
    bool verbose = parser.get_arg("-v", false) | parser.get_arg("--verbose", false); // for verbose output
    nsys = parser.get_arg("-n", false) | parser.get_arg("--output_nsys", false); // for nsys files too
    evols = parser.get_arg("--evolutions", 1'000ull); // number of evolutions to simulate
    bool tty = isatty(STDOUT_FILENO); // if connected to a terminal, output colours, if not, don't
    if (!evols)
        log_error("cannot have zero evolutions.", tty);
    epochs = parser.get_arg("--epochs", 10'000ull); // number of simulation epochs to write per evolution
    if (!epochs)
        log_error("cannot have zero epochs per evolution.", tty);
    uint64_t iters = parser.get_arg("--iterations", 100'000ull); // number of iterations per epoch
    if (!iters)
        log_error("cannot have zero iterations per epoch.", tty);
    // Number of different mean masses to use:
    uint64_t num_m = parser.get_arg("--num_mass_means", evols /* (uint64_t) ceill(((long double) evols)/100.0l) */);
    // `num_m` has to be checked up here, as the default value for `mass_step` uses it
    if (!num_m)
        log_error("cannot have zero mass means.", tty);
    if (num_m > evols)
        log_error("cannot have more mass means than number of evolutions.", tty);
    const long double start_mass = parser.get_arg("--starting_mass", 1'000.0l); // starting mean mass value
    if (start_mass <= 0)
        log_error("starting mass must be positive.", tty);
    const long double mass_step = parser.get_arg("--mass_step", start_mass/10); // step in mean mass value
    if (mass_step <= 0)
        log_error("mass step must be positive.", tty);
    const long double mass_sd_scaling = parser.get_arg("--mass_sd_scaling", 0.125l); // mass_SD = scaling*mean_mass
    if (mass_sd_scaling <= 0)
        log_error("the scaling factor applied to the mass to produce the S.D. must be positive.", tty);
    // I give the particles a radius of zero, as I will be training the NN to predict the collisionless 3-body problem,
    // and so the particles are treated as point particles:
    const long double radius = parser.get_arg("--radius", 0.0l);
    if (radius < 0)
        log_error("radius cannot be negative.", tty);
    const long double dt = parser.get_arg("--timestep", 1.0l/pow(2, 24)); // time-step used to evolve simulations
    if (dt <= 0)
        log_error("time-step must be positive.", tty);
    bool tscale = !parser.get_arg("--no_timestep_scaling", false);
    auto rem = parser.remaining();
    dts = parser.get_arg("--timestep_scaling", 1.0l);
    dte = parser.get_arg("--timestep_exp", 1.0l);
    if (!tscale && parser.remaining() < rem) {
        if (verbose) {
            if (tty)
                std::cout << "\033[1m\033[38;5;9mWarning: \033[0m\033[38;5;10m\"\033[1m\033[38;5;13m--timestep_scaling"
                             "\033[0m\033[38;5;10m\"\033[38;5;11m and/or "
                             "\033[0m\033[38;5;10m\"\033[1m\033[38;5;13m--timestep_exp\033[0m\033[38;5;10m\""
                             "\033[38;5;11m were/was specified, but "
                             "\033[38;5;10m\"\033[1m\033[38;5;13m--no_timestep_scaling\033[38;5;10m\"\033[38;5;11m "
                             "was also passed. Timestep scaling is activated. If this was not your intention, remove "
                             "\033[38;5;10m\"\033[1m\033[38;5;13m--no_timestep_scaling\033[38;5;10m\"\033[38;5;11m as a"
                             " command-line argument." << std::endl;
            else
                std::cout << "Warning: \"--timestep_scaling\" and/or \"--timestep_exp\" were/was specified, but "
                             "\"--no_timestep_scaling\" was also passed. Timestep scaling is activated. If this was not"
                             " your intention, remove \"--no_timestep_scaling\" as a command-line argument."
                             << std::endl;
        }
        tscale = true;
    }
    if (dts <= 0)
        log_error("time-step scaling must be positive", tty);
    const long double vel_scale = parser.get_arg("--velocity_scaling", 0.5l); // scaling factor for velocities
    if (vel_scale <= 0 || vel_scale > 1)
        log_error("the velocity scaling factor must be within the interval (0,1].", tty);
    const long double box_x = parser.get_arg("--box_width", 2.0l); // length of bounding box along x
    if (box_x <= 0)
        log_error("box width must be positive.", tty);
    const long double box_y = parser.get_arg("--box_length", 2.0l); // length of bounding box along y
    if (box_y <= 0)
        log_error("box length must be positive.", tty);
    const long double box_z = parser.get_arg("--box_height", 2.0l); // length of bounding box along z
    // Note: the box is always centered around the COM, so centre coordinates are never supplied as an argument
    if (box_z <= 0)
        log_error("box height must be positive.", tty);
    const long double eps = parser.get_arg("--softening", 1.0l/powl(2, 25)); // softening length
    if (eps < 0)
        log_error("softening-length must be non-negative.", tty);
    // The "minimum separation" is the minimum starting separation the bodies have to have:
    long double min_sep_def = eps*10; // the default is 10 softening lengths - if the bodies are closer, pos. resampled
    long double min_sep = parser.get_arg("--min_sep", min_sep_def);
    if (min_sep > std::min({box_x, box_y, box_z})/2.0l)
        log_error("the minimum separation between the bodies cannot be greater than half the minimum side "
                  "length of the bounding box.", tty);
    ejf = parser.get_arg("--ejection_factor", std::numeric_limits<long double>::quiet_NaN());
    if (ejf != ejf) {
        if (e_stop)
            ejf = EJ_FACTOR;
    } else {
        if (ejf < 1.0l)
            log_error("the ejection factor must be at least 1.0.", tty);
        e_stop = true; // assume ejection checking is desired if the ejection factor was included
        if (verbose) {
            if (tty)
                std::cout << "\033[1m\033[38;5;9mWarning: \033[0m\033[38;5;10m\"\033[1m\033[38;5;13m--early_stop"
                             "\033[0m\033[38;5;10m\" \033[38;5;11mwas not specified, but ejection factor was passed. "
                             "Early stopping is activated. If this was not your intention, remove the value passed "
                             "through \033[38;5;10m\"\033[1m\033[38;5;13m--ejection_factor\033[38;5;10m\"\033[38;5;11m."
                             << std::endl;
            else
                std::cout << "Warning: \"--early_stop\" was not specified, but ejection factor was passed. Early "
                             "stopping is activated. If this was not your intention, remove the value passed through "
                             "\"--ejection_factor\"." << std::endl;
        }
    }
    if (verbose && 2*radius >= min_sep) {
        if (tty)
            std::cerr << BOLD_TXT_START MAGENTA_TXT("Warning:")
            YELLOW_TXT(" there is a non-zero probability that the particles generated with the given radius of ")
            BLUE_TXT_START << radius << RESET_TXT_FLAGS BOLD_TXT_START GREEN_TXT(" m ") YELLOW_TXT_START
            "will overlap when positions are randomly generated, which would cause termination of the "
            "simulation.\nConsider choosing a smaller radius.\n" << RESET_TXT_FLAGS;
        else
            std::cerr << "Warning: there is a non-zero probability that the particles generated with the given "
                         "radius of " << radius << " m will overlap when positions are randomly generated, which would "
                         "cause termination of the ""simulation.\nConsider choosing a smaller radius.\n";
    }
    const char *apath = parser.get_arg("-o");
    unsigned int max_num_threads = std::thread::hardware_concurrency();
    max_num_threads = max_num_threads ? max_num_threads : 1; // single-thread execution if value is not computable
    unsigned int numt = parser.get_arg("--threads", max_num_threads);
    if (!numt || numt > max_num_threads) // user-defined value cannot be zero or over the max. num. of conc. threads
        numt = max_num_threads;
    if (!parser.empty()) {
        if (tty) {
            std::cerr << BOLD_TXT(RED_TXT("Error!")) YELLOW_TXT(" Unrecognised arguments:") "\n";
            for (const auto &[_, arg] : parser)
                std::cerr << GREEN_TXT("\"") << BLUE_TXT_START << arg << GREEN_TXT("\"") << '\n';
            std::cerr << RESET_TXT_FLAGS;
        }
        else {
            std::cerr << "Error! Unrecognised arguments:\n";
            for (const auto &[_, arg] : parser)
                std::cerr << arg << '\n';
        }
        return 1;
    }
    if (verbose) {
        if (tty)
            std::cout << "\033[38;5;123mNumber of concurrent simulation threads to run: \033[1m\033[38;5;11m" << numt
                      << "\033[0m" << std::endl;
        else
            std::cout << "Number of concurrent simulation threads to run: " << numt << std::endl;
    }
    char *path;
    if (!apath) {
        path = new char[41];
        gtd::strcpy_c(path, "experiment_");
        gtd::strcat_c(path, gtd::get_date_and_time());
    }
    else {
        path = new char[gtd::strlen_c(apath) + 5];
        gtd::strcpy_c(path, apath);
    }
    if (verbose) {
        if (tty)
            std::cout << YELLOW_TXT("Attempting to create \"") BLUE_TXT_START << path <<
                      YELLOW_TXT("\" directory...") << std::endl;
        else
            std::cout << "Attempting to create \"" << path << "\" directory..." << std::endl;
    }
    if (mkdir(path, S_IRWXU | S_IRWXO | S_IRWXG) == -1) {
        if (tty)
            std::cerr << BOLD_TXT(RED_TXT("Error: ")) YELLOW_TXT("could not create \"") BLUE_TXT_START << path <<
            YELLOW_TXT("\" directory.") << std::endl;
        else
            std::cerr << "Error: could not create \"" << path << "\" directory." << std::endl;
        delete [] path;
        return 1;
    }
    if (verbose) {
        if (tty)
            std::cout << YELLOW_TXT("Created ") GREEN_TXT("\"") BLUE_TXT_START << path << GREEN_TXT("\"")
            YELLOW_TXT("\" directory.") << std::endl;
        else
            std::cout << "Created \"" << path << "\" directory." << std::endl;
    }
    if (chdir(path) == -1) {
        if (tty)
            std::cerr << BOLD_TXT(RED_TXT("Error: ")) YELLOW_TXT("could change working directory to \"") BLUE_TXT_START
            << path << YELLOW_TXT("\".") << std::endl;
        else
            std::cerr << "Error: could change working directory to \"" << path << "\"." << std::endl;
        delete [] path;
        return 1;
    }
    gtd::strcat_c(path, ".txt");
    if (verbose) {
        if (tty)
            std::cout << YELLOW_TXT("Attempting to create \"") BLUE_TXT_START << path <<
                      YELLOW_TXT("\" file...") << std::endl;
        else
            std::cout << "Attempting to create \"" << path << "\" file..." << std::endl;
    }
    std::ofstream logf{path, std::ios_base::out | std::ios_base::trunc};
    if (!logf.good()) {
        if (tty)
            std::cerr << BOLD_TXT(RED_TXT("Error: ")) YELLOW_TXT("could not create \"") BLUE_TXT_START << path <<
                  YELLOW_TXT("\" file.") << std::endl;
        else
            std::cerr << "Error: could not create \"" << path << "\" file." << std::endl;
        delete [] path;
        return 1;
    }
    if (verbose) {
        if (tty)
            std::cout << YELLOW_TXT("Writing configurations to \"") BLUE_TXT_START << path <<
                      YELLOW_TXT("\" file...") << std::endl;
        else
            std::cout << "Writing configurations to \"" << path << "\" file..." << std::endl;
    }
    logf << std::boolalpha
         << "Distinct evolutions = " << evols << "\nEpochs per evolution = " << epochs << "\nIterations per epoch = "
         << iters << "\nNumber of mass means = " << num_m << "\nStarting mean mass = " << start_mass
         << " kg\nStep in mean mass = " << mass_step << " kg\nSD fraction of mean = " << mass_sd_scaling
         << "\nRadius of bodies = " << radius << " m\nTime-step = " << dt
         << " s\nTime-step scaling = ";
    if (tscale)
        logf << dts << "\nTime-step mass exponent = " << dte;
    else
        logf << "(time-step variation not activated)\nTime-step mass exponent = (time-step variation not activated)";
    logf << "\nVelocity scaling factor = " << vel_scale
         << "\nBox length along x = " << box_x << " m\nBox length along y = " << box_y << " m\nBox length along z = "
         << box_z << " m\nSoftening length = " << eps << " m\nDefault min. sep. = " << min_sep_def
         << " m\nActual min. sep. = " << min_sep << " m\nEarly stopping = " << e_stop
         << "\nVerbose output = " << verbose << "\nOutput nsys = " << nsys << "\nNumber of parallel simulations = "
         << numt << '\n';
    if (verbose) {
        if (tty)
            std::cout << YELLOW_TXT("Configurations written to \"") BLUE_TXT_START << path <<
                      YELLOW_TXT("\" file.") << std::endl;
        else
            std::cout << "Configurations written to \"" << path << "\" file." << std::endl;
    }
    logf.close();
    delete [] path;
    if (e_stop) {
        if (verbose) {
            if (tty) {
                if (tscale)
                    main_loop<true, true, true, true>(box_x, box_y, box_z, dt, iters, eps, radius, num_m, start_mass,
                                                      mass_sd_scaling, mass_step, min_sep, vel_scale, numt);
                else
                    main_loop<true, true, true, false>(box_x, box_y, box_z, dt, iters, eps, radius, num_m, start_mass,
                                                       mass_sd_scaling, mass_step, min_sep, vel_scale, numt);
            }
            else {
                if (tscale)
                    main_loop<true, true, false, true>(box_x, box_y, box_z, dt, iters, eps, radius, num_m, start_mass,
                                                       mass_sd_scaling, mass_step, min_sep, vel_scale, numt);
                else
                    main_loop<true, true, false, false>(box_x, box_y, box_z, dt, iters, eps, radius, num_m, start_mass,
                                                        mass_sd_scaling, mass_step, min_sep, vel_scale, numt);
            }
        } else {
            if (tty) {
                if (tscale)
                    main_loop<true, false, true, true>(box_x, box_y, box_z, dt, iters, eps, radius, num_m, start_mass,
                                                       mass_sd_scaling, mass_step, min_sep, vel_scale, numt);
                else
                    main_loop<true, false, true, false>(box_x, box_y, box_z, dt, iters, eps, radius, num_m, start_mass,
                                                        mass_sd_scaling, mass_step, min_sep, vel_scale, numt);
            }
            else {
                if (tscale)
                    main_loop<true, false, false, true>(box_x, box_y, box_z, dt, iters, eps, radius, num_m, start_mass,
                                                        mass_sd_scaling, mass_step, min_sep, vel_scale, numt);
                else
                    main_loop<true, false, false, false>(box_x, box_y, box_z, dt, iters, eps, radius, num_m, start_mass,
                                                         mass_sd_scaling, mass_step, min_sep, vel_scale, numt);
            }
        }
    } else {
        if (verbose) {
            if (tty) {
                if (tscale)
                    main_loop<false, true, true, true>(box_x, box_y, box_z, dt, iters, eps, radius, num_m, start_mass,
                                                       mass_sd_scaling, mass_step, min_sep, vel_scale, numt);
                else
                    main_loop<false, true, true, false>(box_x, box_y, box_z, dt, iters, eps, radius, num_m, start_mass,
                                                        mass_sd_scaling, mass_step, min_sep, vel_scale, numt);
            }
            else {
                if (tscale)
                    main_loop<false, true, false, true>(box_x, box_y, box_z, dt, iters, eps, radius, num_m, start_mass,
                                                        mass_sd_scaling, mass_step, min_sep, vel_scale, numt);
                else
                    main_loop<false, true, false, false>(box_x, box_y, box_z, dt, iters, eps, radius, num_m, start_mass,
                                                         mass_sd_scaling, mass_step, min_sep, vel_scale, numt);
            }
        } else {
            if (tty) {
                if (tscale)
                    main_loop<false, false, true, true>(box_x, box_y, box_z, dt, iters, eps, radius, num_m, start_mass,
                                                        mass_sd_scaling, mass_step, min_sep, vel_scale, numt);
                else
                    main_loop<false, false, true, false>(box_x, box_y, box_z, dt, iters, eps, radius, num_m, start_mass,
                                                         mass_sd_scaling, mass_step, min_sep, vel_scale, numt);
            }
            else {
                if (tscale)
                    main_loop<false, false, false, true>(box_x, box_y, box_z, dt, iters, eps, radius, num_m, start_mass,
                                                         mass_sd_scaling, mass_step, min_sep, vel_scale, numt);
                else
                    main_loop<false, false, false, false>(box_x, box_y, box_z, dt, iters, eps, radius, num_m,start_mass,
                                                          mass_sd_scaling, mass_step, min_sep, vel_scale, numt);
            }
        }
    }
    if (tty)
        std::cout << "\033[1m\033[38;5;21m------------------------\nAll evolutions finished.\n------------------------"
                     "\n" BOLD_TXT_START MAGENTA_TXT_START << ejected << RESET_TXT_FLAGS BLUE_TXT_START " out of "
                     BOLD_TXT_START MAGENTA_TXT_START << evols
                  << RESET_TXT_FLAGS BLUE_TXT_START " evolutions resulted in an ejection.\n";
    else
        std::cout << "------------------------\nAll evolutions finished.\n------------------------\n" << ejected
                  << " out of " << evols << " evolutions resulted in an ejection.\n";
    return 0;
}
