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
#include <iomanip>

#define EJ_FACTOR 8.0l // minimum ratio of largest distance between bodies to second largest required for ej. check

#define F3BOD_SIZE(fltSize, numEpochs) (12 + 4*sizeof(long double) + 6*fltSize*(1 + 3*numEpochs))

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
    bool will_eject(const gtd::sys &sys, long double ej_factor) {
        if (sys.num_bodies() != 3)
            throw std::invalid_argument{"Error: system object passed must have 3 bodies.\n"};
        const gtd::bod_0f &b1 = sys[0];
        const gtd::bod_0f &b2 = sys[1];
        const gtd::bod_0f &b3 = sys[2];
        const gtd::vec3 p1 = b1.pos();
        const gtd::vec3 p2 = b2.pos();
        const gtd::vec3 p3 = b3.pos();
        long double d12 = (p1 - p2).magnitude();
        long double d13 = (p1 - p3).magnitude();
        long double d23 = (p2 - p3).magnitude();
        long double ej_d12 = ej_factor*d12;
        long double ej_d13 = ej_factor*d13;
        long double ej_d23 = ej_factor*d23;
        const gtd::bod_0f *ejbod;
        const gtd::bod_0f *rbod1;
        const gtd::bod_0f *rbod2;
        if (d12 > ej_d23 && d13 > ej_d23) {
            ejbod = &b1;
            rbod1 = &b2;
            rbod2 = &b3;
        } else if (d13 > ej_d12 && d23 > ej_d12) {
            ejbod = &b3;
            rbod1 = &b1;
            rbod2 = &b2;
        } else if (d12 > ej_d13 && d23 > ej_d13) {
            ejbod = &b2;
            rbod1 = &b1;
            rbod2 = &b3;
        } else
            return false; /*
        long double d1 = p1.magnitude();
        long double d2 = p2.magnitude();
        long double d3 = p3.magnitude();
        if (d1 > ej_factor*std::max(d2, d3)) {
            std::cout << "d1: " << d1 << ", max: " << std::max(d2, d3) << std::endl;
            ejbod = &b1;
            rbod1 = &b2;
            rbod2 = &b3;
        } else if (d2 > ej_factor*std::max(d1, d3)) {
            std::cout << "d2: " << d2 << ", max: " << std::max(d1, d3) << std::endl;
            ejbod = &b2;
            rbod1 = &b1;
            rbod2 = &b3;
        } else if (d3 > ej_factor*std::max(d1, d2)) {
            std::cout << "d3: " << d3 << ", max: " << std::max(d1, d2) << std::endl;
            ejbod = &b3;
            rbod1 = &b1;
            rbod2 = &b2;
        } else
            return false; */
        auto [ev1, ev2] = esc_vel_mags_com(ejbod->mass(), // mass of possibly ejected body
                                           (rbod1->mass() + rbod2->mass()), // total mass of other two bodies
                                           // distance between the pos. ej. body and the COM of the other two:
                                           (ejbod->pos() - gtd::com(*rbod1, *rbod2)).magnitude());
        // gtd::vec3 other_two_com_pos = gtd::com(rbod1, rbod2);
        gtd::vec3 other_two_com_vel = gtd::vel_com(*rbod1, *rbod2);
        if (ejbod->vel().magnitude() >= ev1 || other_two_com_vel.magnitude() >= ev2) // test both in case of round. err.
            return true;
        return false;
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
        size_mismatch_error() : invalid_3bod_format{"Error: floating point data type size does not match \"sizeof(T)\".\n"}
        {}
        explicit size_mismatch_error(const char *msg) : invalid_3bod_format{msg} {}
    };
    class invalid_3bod_ld_size : public invalid_3bod_format {
    public:
        invalid_3bod_ld_size() : invalid_3bod_format{"Error: reported \"sizeof(long double)\" in file does not match "
                                                     "\"sizeof(long double)\" on this platform.\n"} {}
        explicit invalid_3bod_ld_size(const char *msg) : invalid_3bod_format{msg} {}
    };
    class invalid_3bod_T_size : public invalid_3bod_format {
    public:
        invalid_3bod_T_size() : invalid_3bod_format{"Error: reported \"sizeof(T)\" in file does not match "
                                                     "\"sizeof(T)\" on this platform.\n"} {}
        explicit invalid_3bod_T_size(const char *msg) : invalid_3bod_format{msg} {}
    };
    class iterator_index_error : public std::invalid_argument {
    public:
        iterator_index_error() :
        std::invalid_argument{"Error: the required offset for indexing the iterator is out-of-bounds.\n"} {}
        explicit iterator_index_error(const char *msg) : std::invalid_argument{msg} {}
    };
    template <isNumWrapper T, bool, bool>
    requires (std::is_floating_point_v<T> && sizeof(T) <= 255 && sizeof(long double) < 127)
    class f3bodr;
    template <isNumWrapper T = long double>
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
            hdr.dt = sys->timestep()*sys->iters();
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
                if (!rem) { // close file if there are zero entries remaining to be added after tot. num. epochs changed
                    file->close();
                    delete file;
                    file = nullptr;
                }
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
        template <isNumWrapper U, bool, bool>
        requires (std::is_floating_point_v<U> && sizeof(U) <= 255 && sizeof(long double) < 127)
        friend class f3bodr;
        template <bool, bool>
        friend class entry_it;
    };
    // `alloc` dictates whether to store the entire contents in memory and `chk` whether to check iterator indexing
    template <isNumWrapper T, bool alloc = true, bool bnd_chk = true>
    requires (std::is_floating_point_v<T> && sizeof(T) <= 255 && sizeof(long double) < 127)
    class f3bodr {
    public:
        using entry_type = typename f3bodw<T>::entry;
    private:
        template <bool _alloc = true, bool _chk = true> //not pretty, but I avoid declaring unnecessary member variables
        class entry_it { // iterator over entries
            uint64_t tot; // total number of iterations in .3bod file
            entry_type *_arr; // array allocated to hold all positions and velocities at each iteration
            entry_type *_ptr; // pointer to current entry within array
            int64_t pos{}; // position of current entry within array
            bool cpy = false; // boolean indicating whether `*this` has been copy constructed
            // private ctors only for internal use:
            entry_it() : tot{0}, _arr{nullptr}, _ptr{nullptr}, cpy{true} {}
            entry_it(uint64_t _tot, entry_type *_array, entry_type *_entry, int64_t _pos, bool _copy) noexcept :
            tot{_tot}, _arr{_array}, _ptr{_entry}, pos{_pos}, cpy{_copy} {}
            entry_it begin_it() const noexcept {
                return {tot, _arr, _arr, 0, true};
            }
            entry_it end_it() const noexcept {
                return {tot, _arr, _arr + tot, tot, true};
            }
        public:
            entry_it(std::ifstream &_in, uint64_t tot_num) : tot{tot_num} {
                if (!_in.good())
                    throw std::ios_base::failure{"Error: could not read from file.\n"};
                this->_arr = new entry_type[tot_num];
                _in.read((char *) this->_arr, sizeof(entry_type)*tot_num);
                _ptr = _arr; // + `pos` is already zero
            }
            entry_it(entry_it &&other) noexcept : tot{other.tot}, _arr{other._arr}, _ptr{other._ptr}, pos{other.pos} {
                other.tot = 0;
                other._arr = nullptr;
                other._ptr = nullptr;
                other.pos = -1;
                other.cpy = true;
                // memory moved here, so `this->cpy` must remain false
            }
            entry_it(const entry_it &other) :
                tot{other.tot}, _arr{other._arr}, _ptr{other._ptr}, pos{other.pos}, cpy{true} {}
            int64_t position() const noexcept {
                return this->pos;
            }
            const entry_type &operator*() const noexcept(!_chk) {
                if constexpr (_chk)
                    if (pos < 0 || pos >= tot) // cannot dereference `end()` iterator
                        throw iterator_index_error{"Error: cannot dereference an out-of-bounds iterator.\n"};
                return *_ptr;
            }
            const entry_type *operator->() const noexcept(!_chk) {
                if constexpr (_chk)
                    if (pos < 0 || pos >= tot)
                        throw iterator_index_error{"Error: cannot dereference an out-of-bounds iterator.\n"};
                return _ptr;
            }
            const entry_type &at(uint64_t index) const noexcept(!_chk) {
                if constexpr (_chk)
                    if (index >= tot)
                        throw iterator_index_error{"Error: index out of bounds.\n"};
                return *(_arr + index);
            }
            entry_it operator++(int) noexcept(!_chk) {
                if constexpr (_chk)
                    if (pos >= tot) // not if equal to `tot - 1`: need an `end()` iterator
                        throw iterator_index_error{};
                entry_it _copy = *this;
                ++_ptr;
                ++pos;
                return _copy;
            }
            entry_it &operator++() noexcept(!_chk) {
                if constexpr (_chk)
                    if (pos >= tot)
                        throw iterator_index_error{};
                ++_ptr;
                ++pos;
                return *this;
            }
            entry_it operator--(int) noexcept(!_chk) {
                if constexpr (_chk)
                    if (pos <= 0)
                        throw iterator_index_error{};
                entry_it _copy = *this;
                --_ptr;
                --pos;
                return _copy;
            }
            entry_it &operator--() noexcept(!_chk) {
                if constexpr (_chk)
                    if (pos <= 0)
                        throw iterator_index_error{};
                --_ptr;
                --pos;
                return *this;
            }
            entry_it &operator=(const entry_it &other) {
                if (&other == this)
                    return *this;
                this->tot = other.tot;
                this->_arr = other._arr;
                this->_ptr = other._ptr;
                this->pos = other.pos;
                this->cpy = true;
                return *this;
            }
            entry_it &operator=(entry_it &&other) noexcept {
                if (&other == this) // would be possible if `std::move` were used
                    return *this;
                this->tot = other.tot;
                this->_arr = other._arr;
                this->_ptr = other._ptr;
                this->pos = other.pos;
                this->cpy = false;
                other.tot = 0;
                other._arr = nullptr;
                other._ptr = nullptr;
                other.pos = -1;
                other.cpy = true;
                return *this;
            }
            void stop() {
                if (!cpy)
                    delete [] this->_arr;
                tot = 0;
                _arr = nullptr;
                _ptr = nullptr;
                pos = -1;
                cpy = true;
            }
            ~entry_it() {
                if (!cpy) // avoids attempting to deallocate the same memory twice, which is undefined behaviour
                    delete [] this->_arr;
            }
            friend bool operator==(const entry_it &it1, const entry_it &it2) {
                // not checking other members as only two iterators pointing to same array could point to same address:
                return it1._ptr == it2._ptr;
            }
            friend bool operator!=(const entry_it &it1, const entry_it &it2) {
                return it1._ptr != it2._ptr;
            }
            friend bool operator>(const entry_it &it1, const entry_it &it2) {
                return it1._ptr > it2._ptr;
            }
            friend bool operator<(const entry_it &it1, const entry_it &it2) {
                return it1._ptr < it2._ptr;
            }
            friend bool operator>=(const entry_it &it1, const entry_it &it2) {
                return it1._ptr >= it2._ptr;
            }
            friend bool operator<=(const entry_it &it1, const entry_it &it2) {
                return it1._ptr <= it2._ptr;
            }
            friend entry_it operator+(const entry_it &it1, int64_t offset) {
                if constexpr (_chk) {
                    int64_t new_pos = it1.pos + offset;
                    if (new_pos < 0 || new_pos > it1.tot)
                        throw iterator_index_error{};
                    return {it1.tot, it1._arr, it1._ptr + offset, new_pos, true};
                }
                return {it1.tot, it1._arr, it1._ptr + offset, it1.pos + offset, true};
            }
            friend entry_it operator+(int64_t offset, const entry_it &it1) {
                if constexpr (_chk) {
                    int64_t new_pos = it1.pos + offset;
                    if (new_pos < 0 || new_pos > it1.tot)
                        throw iterator_index_error{};
                    return {it1.tot, it1._arr, it1._ptr + offset, new_pos, true};
                }
                return {it1.tot, it1._arr, it1._ptr + offset, it1.pos + offset, true};
            }
            friend entry_it operator-(const entry_it &it1, int64_t offset) {
                if constexpr (_chk) {
                    int64_t new_pos = it1.pos - offset;
                    if (new_pos < 0 || new_pos > it1.tot)
                        throw iterator_index_error{};
                    return {it1.tot, it1._arr, it1._ptr - offset, new_pos, true};
                }
                return {it1.tot, it1._arr, it1._ptr - offset, it1.pos - offset, true};
            }
            friend entry_it operator-(int64_t offset, const entry_it &it1) {
                if constexpr (_chk) {
                    int64_t new_pos = it1.pos - offset;
                    if (new_pos < 0 || new_pos > it1.tot)
                        throw iterator_index_error{};
                    return {it1.tot, it1._arr, it1._ptr - offset, new_pos, true};
                }
                return {it1.tot, it1._arr, it1._ptr - offset, it1.pos - offset, true};
            }
            template <isNumWrapper U, bool, bool>
            requires (std::is_floating_point_v<U> && sizeof(U) <= 255 && sizeof(long double) < 127)
            friend class f3bodr;
            template <typename U> requires (std::is_floating_point_v<U>)
            friend class normaliser;
            template <typename U> requires (std::is_floating_point_v<U>)
            friend class mass_normaliser;
            template <typename U> requires (std::is_floating_point_v<U>)
            friend class log_mass_normaliser;
            template <typename U> requires (std::is_floating_point_v<U>)
            friend class time_normaliser;
            template <typename U> requires (std::is_floating_point_v<U>)
            friend class pos_normaliser;
            template <typename U> requires (std::is_floating_point_v<U>)
            friend class vel_normaliser;
            template <typename U> requires (std::is_floating_point_v<U>)
            friend class log_vel_normaliser;
        };
        template <bool _chk>
        class entry_it<false, _chk> {
            /* Decrementing the iterator to below `begin()` when `_chk == false` is UB. */
            uint64_t tot;
            mutable std::ifstream in; //the position of the pointer within the stream is always one entry ahead of `pos`
            entry_type _e;
            entry_type *_ptr = &_e; // is always pointer to `_e`, except when iterator is `end()` (then `nullptr`)
            int64_t pos{}; // `pos` gives the index of the entry `_e` within the .3bod file
            std::string fname; // need filename to make class copyable (`std::ifstream` objects cannot be copied)
            entry_it() : tot{0}, _ptr{nullptr} /*, fname{nullptr} */ {}
            entry_it(uint64_t _tot, const char *path, int64_t _pos) :
            tot{_tot}, in{path, std::ios_base::in | std::ios_base::binary}, _e{}, _ptr{&_e}, pos{_pos},
            fname{path} {
                if (!in.good())
                    throw std::ios_base::failure{"Error: failure to read from file provided.\n"};
                if (pos < tot) {
                    in.seekg(sizeof(hdr_t) + pos*sizeof(entry_type));
                    in.read((char *) &_e, sizeof(entry_type));
                }
                else { /*
                    if constexpr (_chk)
                        if (pos > tot)
                            throw iterator_index_error{}; */
                    in.seekg(0, std::ios_base::end).setstate(std::ios_base::eofbit);
                    _ptr = nullptr;
                }
            }
            entry_it begin_it() const noexcept {
                return {fname, tot};
            }
            entry_it end_it() const noexcept {
                return {tot, fname, tot};
            }
        public:
            // using entry_type = typename f3bodw<T>::entry;
            entry_it(const char *path, uint64_t _tot) :
                tot{_tot}, in{path, std::ios_base::in | std::ios_base::binary}, fname{path} {
                if (!in.good())
                    throw std::ios_base::failure{"Error reading from .3bod file.\n"};
                in.seekg(sizeof(hdr_t)); // seek to offset of where array of pos./vel. starts
                in.read((char *) &_e, sizeof(entry_type)); // read in first entry
            }
            entry_it(const entry_it &other) : tot{other.tot}, _e{other._e}, pos{other.pos},
            in{other.fname, std::ios_base::in | std::ios_base::binary}, fname{other.fname} {
                if (!in.good())
                    throw std::ios_base::failure{"Error: could not open file associated with other iterator.\n"};
                if (other.pos < tot) {
                    this->in.seekg(sizeof(hdr_t) + (other.pos + 1)*sizeof(entry_type));
                }
                else {
                    this->in.seekg(0, std::ios_base::end).setstate(std::ios_base::eofbit);
                    _ptr = nullptr;
                }
            }
            entry_it(entry_it &&other) : tot{other.tot}, in{std::move(other.in)}, _e{other._e},_ptr{&_e},pos{other.pos},
            fname{std::move(other.fname)} { // not marked `noexcept` because `std::ifstream` move ctor isn't either
                other._ptr = nullptr;
                other.pos = -1;
                // other.fname = nullptr;
            }
            int64_t position() const noexcept {
                return this->pos;
            }
            const entry_type &operator*() const& noexcept(!_chk) {
                if constexpr (_chk)
                    if (pos < 0 || pos >= tot) // cannot dereference `end()` iterator
                        throw iterator_index_error{"Error: cannot dereference an out-of-bounds iterator.\n"};
                return *_ptr; // causes seg. fault if `_chk` is `false` and iterator is `end()` iterator
            }
            const entry_type &operator*() const&& noexcept(!_chk) { // overload for when `*this` is an r-value
                /* NOT THREAD-SAFE: consider return-by-value */
                static entry_type entry;
                if constexpr (_chk)
                    if (pos < 0 || pos >= tot)
                        throw iterator_index_error{"Error: cannot dereference an out-of-bounds iterator.\n"};
                return (entry = *_ptr); // maybe change to `_e`, but almost want to keep UB for when `*this` is `end()`
            }
            const entry_type *operator->() const& noexcept(!_chk) {
                if constexpr (_chk)
                    if (pos < 0 || pos >= tot)
                        throw iterator_index_error{"Error: cannot dereference an out-of-bounds iterator.\n"};
                return _ptr;
            }
            const entry_type *operator->() const&& noexcept(!_chk) {
                static entry_type entry;
                if constexpr (_chk)
                    if (pos < 0 || pos >= tot)
                        throw iterator_index_error{"Error: cannot dereference an out-of-bounds iterator.\n"};
                return &(entry = *_ptr);
            }
            const entry_type &at(uint64_t index) const { // consider changing this to return by value
                /* CURRENTLY NOT THREAD-SAFE */
                if constexpr (_chk)
                    if (index >= tot)
                        throw iterator_index_error{"Error: index out of bounds.\n"};
                static entry_type entry;
                std::ios_base::iostate state = in.rdstate();
                std::ifstream::pos_type _pos = in.tellg();
                in.clear();
                in.seekg(sizeof(hdr_t) + index*sizeof(entry_type));
                in.read((char *) &entry, sizeof(entry_type));
                in.seekg(_pos);
                in.setstate(state);
                return entry;
            }
            entry_it operator++(int) noexcept(!_chk) {
                entry_it _copy = *this;
                this->operator++();
                return _copy;
            }
            entry_it &operator++() noexcept(!_chk) {
                if constexpr (_chk)
                    if (pos >= tot) // not if equal to `tot - 1`: need an `end()` iterator
                        throw iterator_index_error{};
                if (pos == tot - 1) { // CASE FOR BECOMING `end()` ITERATOR
                    _ptr = nullptr; // no further element, so must be `nullptr`
                    in.setstate(std::ios_base::eofbit);
                }
                else {
                    _ptr = &_e;
                    in.read((char *) &_e, sizeof(entry_type));
                }
                ++pos;
                return *this;
            }
            entry_it operator--(int) noexcept(!_chk) {
                entry_it _copy = *this;
                this->operator--();
                return _copy;
            }
            entry_it &operator--() noexcept(!_chk) {
                if constexpr (_chk) {
                    if (pos <= 0) // not if equal to `tot - 1`: need an `end()` iterator
                        throw iterator_index_error{};
                }
                else {
                    if (pos > tot) { // case for `pos` being above the `end()` iterator - nothing else should be done
                        --pos;
                        return *this;
                    }
                }
                if (pos == tot) {
                    _ptr = &_e;
                    in.clear(); // clear the `std::ios_base::eofbit` bit
                    in.seekg(-sizeof(entry_type), std::ios_base::end);
                }
                else
                    in.seekg(-2*sizeof(entry_type), std::ios_base::cur);
                in.read((char *) &_e, sizeof(entry_type));
                --pos;
                return *this;
            }
            entry_it &operator=(const entry_it &other) {
                if (&other == this)
                    return *this;
                this->in.close();
                this->in.open(other.fname, std::ios_base::in | std::ios_base::binary);
                if (!in.good())
                    throw std::ios_base::failure{"Error: could not open file associated with other iterator.\n"};
                this->tot = other.tot;
                this->fname = other.fname;
                if (other.pos < tot) {
                    this->in.seekg(sizeof(hdr_t) + (other.pos + 1)*sizeof(entry_type));
                    _ptr = &_e;
                }
                else {
                    this->in.seekg(0, std::ios_base::end).setstate(std::ios_base::eofbit);
                    _ptr = nullptr;
                }
                this->_e = other._e;
                this->pos = other.pos;
                return *this;
            }
            entry_it &operator=(entry_it &&other) { // cannot be marked `noexcept` as `std::ifstream::operator=` isn't
                if (&other == this)
                    return *this;
                this->tot = other.tot;
                this->_e = other._e;
                if (other.pos >= other.tot)
                    this->_ptr = nullptr;
                else
                    this->_ptr = &_e;
                this->in = std::move(other.in);
                this->pos = other.pos;
                this->fname = std::move(other.fname);
                // other.fname = nullptr;
                other.tot = 0;
                other._ptr = nullptr;
                other.pos = -1;
                return *this;
            }
            void stop() {
                if (in.is_open())
                    in.close();
                this->tot = 0;
                // this->fname = nullptr;
                this->fname.clear();
                this->fname.shrink_to_fit();
                this->_ptr = nullptr;
                this->pos = -1;
            }
            ~entry_it() {
                if (in.is_open()) // in case `*this` was moved
                    in.close();
                // this->fname.clear();
                // this->fname.shrink_to_fit();
            }
            // All comparison functions return `false` (apart from `operator!=`) if the iterators point to diff. files
            friend bool operator==(const entry_it &it1, const entry_it &it2) {
                // not checking other members as only two iterators pointing to same array could point to same address:
                return it1.pos == it2.pos && it1.fname == it2.fname; // IMPROVE THIS, USE IDS OR HASHING
            }
            friend bool operator!=(const entry_it &it1, const entry_it &it2) {
                return it1.pos != it2.pos || it1.fname != it2.fname; // thank god for short-circuit evaluation
            }
            friend bool operator>(const entry_it &it1, const entry_it &it2) {
                return it1.pos > it2.pos && it1.fname == it2.fname;
            }
            friend bool operator<(const entry_it &it1, const entry_it &it2) {
                return it1.pos < it2.pos && it1.fname == it2.fname;
            }
            friend bool operator>=(const entry_it &it1, const entry_it &it2) {
                return it1.pos >= it2.pos && it1.fname == it2.fname;
            }
            friend bool operator<=(const entry_it &it1, const entry_it &it2) {
                return it1.pos <= it2.pos && it1.fname == it2.fname;
            }
            friend entry_it operator+(const entry_it &it1, int64_t offset) {
                if constexpr (_chk) {
                    int64_t new_pos = it1.pos + offset;
                    if (new_pos < 0 || new_pos > it1.tot)
                        throw iterator_index_error{};
                    return {it1.tot, it1.fname, new_pos};
                }
                return {it1.tot, it1.fname, it1.pos + offset};
            }
            friend entry_it operator+(int64_t offset, const entry_it &it1) {
                if constexpr (_chk) {
                    int64_t new_pos = it1.pos + offset;
                    if (new_pos < 0 || new_pos > it1.tot)
                        throw iterator_index_error{};
                    return {it1.tot, it1.fname, new_pos};
                }
                return {it1.tot, it1.fname, it1.pos + offset};
            }
            friend entry_it operator-(const entry_it &it1, int64_t offset) {
                if constexpr (_chk) {
                    int64_t new_pos = it1.pos - offset;
                    if (new_pos < 0 || new_pos > it1.tot)
                        throw iterator_index_error{};
                    return {it1.tot, it1.fname, new_pos};
                }
                return {it1.tot, it1.fname, it1.pos - offset};
            }
            friend entry_it operator-(int64_t offset, const entry_it &it1) {
                if constexpr (_chk) {
                    int64_t new_pos = it1.pos - offset;
                    if (new_pos < 0 || new_pos > it1.tot)
                        throw iterator_index_error{};
                    return {it1.tot, it1.fname, new_pos};
                }
                return {it1.tot, it1.fname, it1.pos - offset};
            }
            template <isNumWrapper U, bool, bool>
            requires (std::is_floating_point_v<U> && sizeof(U) <= 255 && sizeof(long double) < 127)
            friend class f3bodr;
        };
        typename f3bodw<T>::header hdr; // might consider removing this, as it involves unnecessary initialisation here
        entry_it<alloc, bnd_chk> _it;
        uint64_t file_size;
    public:
        using hdr_t = typename f3bodw<T>::header;
        using iterator_type = entry_it<alloc, bnd_chk>;
        explicit f3bodr(const char *path) {
            struct stat buff{};
            if (stat(path, &buff) == -1)
                throw std::invalid_argument{"Error: .3bod/.3bodpp file provided does not exist.\n"};
            if (!S_ISREG(buff.st_mode))
                throw std::invalid_argument{"Error: .3bod/.3bodpp file provided is not a regular file.\n"};
            if (buff.st_size < sizeof(hdr_t))
                throw invalid_3bod_format{"Error: .3bod/.3bodpp format invalid (at least for current type 'T'.\n"};
            std::ifstream in{path, std::ios_base::in | std::ios_base::binary};
            in.read((char *) &hdr, sizeof(hdr_t));
            if (hdr.hd[0] != '3' || (hdr.hd[1] != 'B' && hdr.hd[1] != 'P'))
                throw invalid_3bod_format{"Error: invalid .3bod/.3bodpp header.\n"};
            if (hdr.flt_size != sizeof(T))
                throw invalid_3bod_T_size{"Error: reported size of f.p. data type does not match \"sizeof(T)\".\n"};
            if ((hdr.lds_coll >> 1) != sizeof(long double))
                throw invalid_3bod_ld_size{"Error: reported size of a \"long double\" does not match "
                                           "\"sizeof(long double)\" on this platform.\n"};
            if (buff.st_size != (file_size = F3BOD_SIZE(sizeof(T), hdr.N)))
                throw invalid_3bod_format{"Error: invalid total file size.\n"};
            if constexpr (alloc) {
                // construct `_it` from reference to `in`, as `_it` will only use it once to read in entire array
                _it = entry_it<alloc, bnd_chk>{in, hdr.N};
                in.close();
            }
            else {
                // construct `_it` from `path`, since `_it` will construct its own `std::ifstream` and keep that
                // shut in before constructing new `entry_it` instance
                in.close();
                _it = entry_it<alloc, bnd_chk>{path, hdr.N};
            }
        }
        const hdr_t &header() const noexcept {
            return this->hdr;
        }
        const entry_type &entry_at(uint64_t index) const noexcept {
            // return *(_it + index); // not efficient
            return _it.at(index); // much more efficient!!
        }
        iterator_type begin() const {
            return this->_it.begin_it(); // could just return `_it` itself
        }
        iterator_type end() const {
            return this->_it.end_it();
        }
        void close() {
            this->_it.stop();
        }
        template <typename U> requires (std::is_floating_point_v<U>)
        friend class normaliser;
        template <typename U> requires (std::is_floating_point_v<U>)
        friend class mass_normaliser;
        template <typename U> requires (std::is_floating_point_v<U>)
        friend class log_mass_normaliser;
        template <typename U> requires (std::is_floating_point_v<U>)
        friend class time_normaliser;
        template <typename U> requires (std::is_floating_point_v<U>)
        friend class pos_normaliser;
        template <typename U> requires (std::is_floating_point_v<U>)
        friend class vel_normaliser;
        template <typename U> requires (std::is_floating_point_v<U>)
        friend class log_vel_normaliser;
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
