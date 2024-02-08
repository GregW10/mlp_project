#ifndef DATSUP_HPP
#define DATSUP_HPP

#if !defined(__APPLE__) && !defined(__linux__)
#error "Can only compile on UNIX systems."
#endif

#define GREGSYS_SOFTENING

#include "../glib/nbod/gregsys.hpp"
#include "../glib/misc/gregparse.hpp"
#include <thread>
#include <mutex>

namespace gtd {
    // Here I define a function used to sample values from a log-normal distribution with the mean and standard
    // deviation being those of the original Gaussian (hence the corresponding mean and standard deviation of the
    // log-normal distribution itself have to be transformed into those of the original Gaussian before passing them to
    // the function):
    // long double log_normal_sample(long double mu_gaussian, long double sigma_gaussian) {
    //
    // }
    class not_3body_system_error : public std::logic_error {
    public:
        not_3body_system_error() : std::logic_error{"Error: the system object does not contain 3 bodies.\n"} {}
        explicit not_3body_system_error(const char *msg) : std::logic_error{msg} {}
    };
    class f3bod {
        std::ofstream *file; // pointer to `std::ofstream` object attached to same .3bod file for instance's lifespan
        const gtd::sys *sys; // pointer to the `gtd::sys` object which contains the evolving 3-body system
        void check_3bods() const {
            if (this->sys->num_bodies() != 3)
                throw not_3body_system_error{};
        }
        void write_header() {

        }
    public:
        f3bod(const char *fpath, const gtd::sys *_sys) :
        file{new std::ofstream{fpath, std::ios_base::out | std::ios_base::trunc}}, sys{_sys} {
            if (!file->good()) {
                delete file;
                throw std::ios_base::failure{"Error: could not open file.\n"};
            }
            if (!this->sys) {
                delete file;
                throw std::invalid_argument{"Error: nullptr passed as pointer to system object.\n"};
            }
            this->check_3bods();
            // write header
        }
        void add_entry() { // add an entry for the 3 bodies at the current iteration in the system's evolution

        }
        ~f3bod() {
            if (file->is_open())
                file->close();
            delete file;
        }
    };
}

#endif
