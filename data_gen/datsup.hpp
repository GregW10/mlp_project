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
    std::pair<long double, long double> esc_vel_mags_com(long double m1, long double m2, long double sep) {
        long double root = sqrtl((2*gtd::sys::G_SI)/(sep*(m1 + m2)));
        return {m2*root, m1*root}; // escape velocities for body 1 and body 2, respectively, in the COM frame
    }
    // This function returns the "beta factor" - the of the kinetic energy to the magnitude of the potential energy of a
    // system. I use this as a kind of normalisation of the velocities of the 3 bodies to ensure they never have enough
    // energy to escape the system (as this would result in poor training data, with one body flung off in a straight
    // line and only two bodies remaining):
    long double beta_factor(const gtd::sys &sys) {
        if (!sys.num_bodies())
            throw std::invalid_argument{"Error: the system cannot be empty.\n"};
        return sys.kinetic_energy()/(-sys.potential_energy()); // PE always negative, so take negative to get abs. val.
    }
    // Below I define a function which takes the mean and standard deviation of a log-normal distribution and returns
    // the mean and standard deviation of the corresponding Gaussian distribution. I derived these formulae from first
    // principles, calculating E[x] and E[x^2] via integration.
    template <typename T>
    std::pair<T, T> lognormal_to_gauss(const T &mu_ln, const T &sigma_ln) {
        T sig_over_mu = sigma_ln/mu_ln;
        T squared_p1 = 1 + sig_over_mu*sig_over_mu;
        return {logl(mu_ln/sqrtl(squared_p1)), sqrtl(logl(squared_p1))};
    }
    class not_3body_system_error : public std::logic_error {
    public:
        not_3body_system_error() : std::logic_error{"Error: the system object does not contain 3 bodies.\n"} {}
        explicit not_3body_system_error(const char *msg) : std::logic_error{msg} {}
    };
    class missing_entries_error : public std::logic_error {
    public:
        missing_entries_error() : std::logic_error{"Error: not all entries have been written to the .3bod file.\n"} {}
        explicit missing_entries_error(const char *msg) : std::logic_error{msg} {}
    };
    class entry_limit_error : public std::logic_error {
    public:
        entry_limit_error() : std::logic_error{"Error: entry limit reached. To add more entries, change the total "
                                               "number of entries in the .3bod file.\n"} {}
        explicit entry_limit_error(const char *msg) : std::logic_error{msg} {}
    };
    template <typename T> requires (std::is_floating_point_v<T> && sizeof(T) <= 255 && sizeof(long double) < 127)
    class f3bod {
        std::ofstream *file; // pointer to `std::ofstream` object attached to same .3bod file for instance's lifespan
        const gtd::sys *sys; // pointer to the `gtd::sys` object which contains the evolving 3-body system
        uint64_t _epochs;
        uint64_t rem;
        void check_all() const {
            if (!file->good()) {
                delete file;
                throw std::ios_base::failure{"Error: could not open file.\n"};
            }
            if (!this->sys) {
                delete file;
                throw std::invalid_argument{"Error: nullptr passed as pointer to system object.\n"};
            }
            if (!this->_epochs) {
                delete file;
                throw std::invalid_argument{"Error: zero epochs.\n"};
            }
            if (this->sys->num_bodies() != 3)
                throw not_3body_system_error{};
        }
        void write_header() {
            header hdr{};
            hdr.dt = sys->timestep();
            hdr.N = this->_epochs;
            T *ptr = hdr.masses;
            for (const auto &b : *this->sys)
                *ptr++ = b.mass(); // this code will break if I make `gtd::body<>::mass` return a value and not a ref.
            for (const auto &b : *this->sys) // permissible to continue with `ptr` here, as the `struct` is contig.
                *ptr++ = b.rad(); // this code will break if I make `gtd::body<>::rad` return a value and not a ref.
            if constexpr (std::same_as<T, long double>) {
                for (const auto &b : *this->sys) // if `T` is `long double`, then can continue with `ptr`
                    *ptr++ = b.restitution();
            } else {
                long double *lptr = hdr.nCORs;
                for (const auto &b : *this->sys) // if `T` is NOT `long double`, then must use new `long double` pointer
                    *lptr++ = b.restitution();
            }
            file->write((char *) &hdr, sizeof(header));
        }
    public:
#pragma pack(push, 1)
        struct header {
            char hd[2] = {'3', 'B'};
            unsigned char flt_size{sizeof(T)};
            unsigned char lds_coll{sizeof(long double) << 1}; // should be 32: `sizeof(long double)` in 7 MSBs, bool b.
            long double dt{};
            uint64_t N{};
            T masses[3]{};
            T radii[3]{};
            long double nCORs[3]{};
        };
        struct entry {
            T positions[9]{}; // x1, y1, z1, x2, y2, z2, x3, y3, z3
            T velocities[9]{}; // vx1, vy1, vz1, vx2, vy2, vz2, vx3, vy3, vz3
        };
#pragma pack(pop)
        using sys_t = gtd::system<T, T, T>;
        f3bod(const char *fpath, const sys_t *_sys, uint64_t num_epochs) :
        file{new std::ofstream{fpath, std::ios_base::out | std::ios_base::trunc}}, sys{_sys}, _epochs{num_epochs},
        rem{num_epochs} {
            this->check_all();
            this->write_header();
        }
        uint64_t epochs() const noexcept {
            return this->_epochs;
        }
        uint64_t epochs(uint64_t new_num) {
            if (!new_num)
                return this->_epochs;
            if (new_num == this->_epochs)
                return new_num;
            if (new_num > this->_epochs) {
                this->rem += new_num - this->_epochs;
                this->_epochs = new_num;
                return new_num;
            }
            if (new_num < this->_epochs) {
                uint64_t added;
                if (new_num < (added = (this->_epochs - rem))) // if `new_num` is less than the entries already added
                    throw entry_limit_error{"Error: new total number of epochs would be less than the number of entries"
                                            " already added to the .3bod file.\n"};
                this->_epochs = new_num;
                this->rem = new_num - added;
                this->file->seekp(4 + sizeof(long double)); // seek to point in the file where I wrote the num. epochs
                this->file->write((char *) &new_num, sizeof(uint64_t));
                this->file->seekp(0, std::ios_base::end); // seek back to the end for further writing
            }
            return new_num;
        }
        uint64_t remaining() const noexcept {
            return this->rem;
        }
        uint64_t remaining(uint64_t new_rem) { // convenience method for adding on a number of entries at the end
            this->epochs(this->_epochs + new_rem - this->rem);
            return this->rem;
        }
        std::ofstream::pos_type add_entry() {
            /* Adds an entry of the positions and velocities of the three bodies to the .3bod file. The method returns
             * the position in the output stream after the addition of the entry.
             * Note: when this is called for the last entry to be added, the file stream is closed and the
             * `std::ofstream` object is deleted, thus rendering the class unusable, though still in a valid state. */
            if (!rem)
                throw entry_limit_error{};
            entry _e;
            T *pptr = _e.positions;
            T *vptr = _e.velocities;
            for (const auto &b : *sys) {
                *pptr++ = b.pos()[0]; // I'd rather not have two nested loops
                *pptr++ = b.pos()[1];
                *pptr++ = b.pos()[2];
                *vptr++ = b.vel()[0];
                *vptr++ = b.vel()[1];
                *vptr++ = b.vel()[2];
            }
            file->write((char *) &_e, sizeof(entry));
            if (!--rem) {
                file->close();
                delete file;
                file = nullptr;
                return 0; // only returns zero if last entry
            }
            return this->file->tellp(); // returns the output position indicator
        }
        bool done() const noexcept {
            return !this->rem;
        }
        ~f3bod() {
            if (this->rem)
                std::cerr << BOLD_TXT(RED_TXT("Warning: "))
                YELLOW_TXT("an `f3bod` object is being destroyed, but it did not finish writing to its associated .3bod"
                           " file.\n");
            if (!file)
                return;
            if (file->is_open())
                file->close();
            delete file;
        }
    };
}

#endif
