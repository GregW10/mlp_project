#ifndef DATSUP_HPP
#define DATSUP_HPP

/* Header file containing important function and class definitions for data generation and visualisation support. */

#if !defined(__APPLE__) && !defined(__linux__)
#error "Can only compile on UNIX systems."
#endif

#define GREGSYS_SOFTENING // this will activate softening in my `gtd::system<>` class

#include "../glib/nbod/gregsys.hpp"
#include "../glib/misc/gregparse.hpp"
#include <thread>
#include <mutex>
#include <condition_variable>

#define EJ_FACTOR 10.0l // minimum ratio of largest distance between bodies to second largest required for ej. check

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
    bool will_eject(const gtd::sys &sys) {
        if (sys.num_bodies() != 3)
            throw std::invalid_argument{"Error: system object passed must have 3 bodies.\n"};
        const gtd::bod_0f &b1 = sys[0];
        const gtd::bod_0f &b2 = sys[1];
        const gtd::bod_0f &b3 = sys[2];
        const gtd::vec3 p1 = b1.pos();
        const gtd::vec3 p2 = b2.pos();
        const gtd::vec3 p3 = b2.pos();
        long double d1 = p1.magnitude();
        long double d2 = p2.magnitude();
        long double d3 = p3.magnitude();
        const gtd::bod_0f *ejbod;
        const gtd::bod_0f *rbod1;
        const gtd::bod_0f *rbod2;
        if (d1 > EJ_FACTOR*std::max(d2, d3)) {
            ejbod = &b1;
            rbod1 = &b2;
            rbod2 = &b3;
        } else if (d2 > EJ_FACTOR*std::max(d1, d3)) {
            ejbod = &b2;
            rbod1 = &b1;
            rbod2 = &b3;
        } else if (d3 > EJ_FACTOR*std::max(d1, d2)) {
            ejbod = &b3;
            rbod1 = &b1;
            rbod2 = &b2;
        }
        else
            return false;
        auto [ev1, ev2] = esc_vel_mags_com(ejbod->mass(), // mass of possibly ejected body
                                           (rbod1->mass() + rbod2->mass()), // total mass of other two bodies
                                           // distance between the pos. ej. body and the COM of the other two:
                                           (ejbod->pos() - gtd::com(*rbod1, *rbod2)).magnitude());
        // gtd::vec3 other_two_com_pos = gtd::com(rbod1, rbod2);
        gtd::vec3 other_two_com_vel = gtd::vel_com(*rbod1, *rbod2);
        if (ejbod->vel().magnitude() >= ev1 || other_two_com_vel.magnitude() >= ev2) // test both in case of round. err.
            return true;
        return false;
        /*
        long double d1 = (b1.pos() - b2.pos()).magnitude();
        long double d2 = (b1.pos() - b3.pos()).magnitude();
        long double d3 = (b2.pos() - b3.pos()).magnitude();
        if (d1 > 10*d3) {

        }
        else if (d2 > 10*d3) {

        }
        else if (d3 > 10*d1) {

        }
        long double max = std::max({d1, d2, d3}); */

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
    class invalid_3bod_format : public std::logic_error {
    public:
        invalid_3bod_format() : std::logic_error{"Error: invalid .3bod format.\n"} {}
        explicit invalid_3bod_format(const char *msg) : std::logic_error{msg} {}
    };
    class size_mismatch_error : public invalid_3bod_format {
    public:
        size_mismatch_error() : std::logic_error{"Error: floating point data type size does not match \"sizeof(T)\".\n"}
        {}
        explicit size_mismatch_error(const char *msg) : std::logic_error{msg} {}
    };
    template <typename T = long double>
    requires (std::is_floating_point_v<T> && sizeof(T) <= 255 && sizeof(long double) < 127)
    class f3bodw {
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
        f3bodw(const char *fpath, const sys_t *_sys, uint64_t num_epochs) :
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
        ~f3bodw() {
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
    template <typename T, bool alloc = true> // `alloc` dictates whether to store the entire contents in memory
    requires (std::is_floating_point_v<T> && sizeof(T) <= 255 && sizeof(long double) < 127)
    class f3bodr {
        std::ifstream in;
        typename f3bodw<T>::header hdr; // might consider removing this, as it involves unnecessary initialisation here
        class entry_it { // iterator over entries
            uint64_t tot;
            typename f3bodw<T>::entry *_e{};
        public:
            entry_it(std::ifstream &_in, uint64_t tot_num) : tot{tot_num} {}
        };
    public:
        f3bodr(const char *path) : in{path, std::ios_base::in | std::ios_base::binary} {
            struct stat buff{};
            if (stat(path, &buff) == -1)
                throw std::invalid_argument{"Error: .3bod file provided does not exist.\n"};
            if (!S_ISREG(buff.st_mode))
                throw std::invalid_argument{"Error: .3bod file provided is not a regular file.\n"};

        }
    };
    // Here I define a trivial class used to compare two `std::thread` objects, such that `std::thread` objects can be
    // added to an `std::set`, where the `Compare` type is set to `thread_comparator`:
    class thread_comparator {
    public:
        struct is_transparent {};
        bool operator()(const std::thread &t1, const std::thread &t2) const noexcept {
            return t1.get_id() < t2.get_id();
        }
        bool operator()(const std::thread &_t, const std::thread::id &_i) const noexcept {
            return _t.get_id() < _i;
        }
        bool operator()(const std::thread::id &_i, const std::thread &_t) const noexcept {
            return _i < _t.get_id();
        }
    };
}

#endif
