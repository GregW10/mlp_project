#ifndef NNSUP_HPP
#define NNSUP_HPP

#include "../data_management/datsup.hpp"
#include <dirent.h>

#define ABS(val) (val < 0 ? -val : val)

namespace gtd {
    template <typename T> requires (std::is_floating_point_v<T>)
    class normaliser {
    protected:
        using h_type = typename f3bodr<T>::hdr_t;
        using e_type = typename f3bodr<T>::entry_type;
        T _mmass{}; // max. mass
        T _nmass{}; // min. mass
        long double _mtime{}; // longest sim. time - `long double` since time-step is always stored as such
        T _mpos{}; // maximum absolute position component
        T _mvel{}; // maximum absolute velocity component
        virtual h_type *normalise_header(h_type*) const = 0;
        virtual e_type *normalise_entry(e_type*) const = 0;
        virtual consteval uint64_t get_id() const noexcept = 0;
    public:
#pragma pack(push, 1)
        struct normv_file {
            const char header[4] = {'N', 'O', 'R', 'M'};
            const uint32_t _sizeof_T = sizeof(T);
            const uint32_t _sizeof_ld = sizeof(long double);
            const uint32_t _sizeof_u64 = sizeof(uint64_t); // always 8
            uint64_t _id;
            T _maxm;
            T _minm;
            long double _maxt;
            T _maxp;
            T _maxv;
        };
#pragma pack(pop)
        normaliser() = default;
        normaliser(T max_mass, T min_mass, long double max_time, T max_abs_poscomp, T max_abs_velcomp) :
        _mmass{max_mass}, _nmass{min_mass}, _mtime{max_time}, _mpos{max_abs_poscomp}, _mvel{max_abs_velcomp} {
            if (_mmass < 0 || _nmass < 0 || _mtime < 0 || _mpos < 0 || _mvel < 0)
                throw std::invalid_argument{"Error: all maximum values passed must be non-negative.\n"};
            if (_nmass > _mmass)
                throw std::invalid_argument{"Error: minimum mass cannot be greater than maximum mass.\n"};
        }
        normaliser(const char *normv_path) {
            this->load_normv(normv_path);
        }
        T max_mass() const noexcept {
            return _mmass;
        }
        T min_mass() const noexcept {
            return _nmass;
        }
        long double max_time() const noexcept {
            return _mtime;
        }
        T max_pos() const noexcept {
            return _mpos;
        }
        T max_vel() const noexcept {
            return _mvel;
        }
        virtual void max_mass(T _new_val) {
            if (_new_val < 0)
                throw std::invalid_argument{"Error: max. mass val. must be non-negative.\n"};
            _mmass = _new_val;
        }
        virtual void min_mass(T _new_val) {
            if (_new_val < 0)
                throw std::invalid_argument{"Error: min. mass val. must be non-negative.\n"};
            _nmass = _new_val;
        }
        virtual void max_time(long double _new_val) {
            if (_new_val < 0)
                throw std::invalid_argument{"Error: max. time val. must be non-negative.\n"};
            _mtime = _new_val;
        }
        virtual void max_pos(T _new_val) {
            if (_new_val < 0)
                throw std::invalid_argument{"Error: max. pos. val. must be non-negative.\n"};
            _mpos = _new_val;
        }
        virtual void max_vel(T _new_val) {
            if (_new_val < 0)
                throw std::invalid_argument{"Error: max. vel. val. must be non-negative.\n"};
            _mvel = _new_val;
        }
        virtual void operator()(const std::string &input, const std::string &output) {
            f3bodr<T> reader{input.c_str()};
            std::ofstream out{output, std::ios_base::out | std::ios_base::binary};
            if (!out.good()) {
                std::cerr << "Error: could not open file \"" << output << "\".\n";
                exit(1);
            }
            reader.hdr.hd[1] = 'P';
            out.write((char *) normalise_header(&reader.hdr), sizeof(h_type));
            e_type *_entries = reader._it._arr;
            e_type *_entry = _entries;
            uint64_t counter = reader.hdr.N;
            while (counter --> 0)
                normalise_entry(_entry++);
            out.write((char *) _entries, reader.hdr.N*sizeof(e_type));
            out.close();
        }
        std::ofstream::pos_type write_normv(const char *path) const {
            /* Writes the maximum values to a .normv file. The file format of a .normv file is as follows:
             *     - 4 bytes: Characters 'N', 'O', 'R' and 'M'
             *     - 4 bytes: `sizeof(T)`
             *     - 4 bytes: `sizeof(long double)`
             *     - 4 bytes: `sizeof(uint64_t)`
             *     - `sizeof(uint64_t)` bytes: ID of normaliser
             *     - `3*sizeof(T) + sizeof(long double)` bytes: Maximum values, in same order of class */
            if (!path)
                throw std::invalid_argument{"Error: path for .normv file cannot be nullptr.\n"};
            std::ofstream out{path, std::ios_base::out | std::ios_base::binary};
            if (!out.good())
                throw std::ios_base::failure{"Error: could not open .normv file for writing.\n"};
            /* out.write(_normv_hdr, 4*sizeof(char));
            out.write((char *) &_sizeof_T, sizeof(uint32_t));
            out.write((char *) &_sizeof_ld, sizeof(uint32_t));
            out.write((char *) &_sizeof_u64, sizeof(uint32_t));
            out.write((char *) &_id, sizeof(uint64_t));
            out.write((char *) &_mmass, sizeof(T));
            out.write((char *) &_mtime, sizeof(long double));
            out.write((char *) &_mpos, sizeof(T));
            out.write((char *) &_mvel, sizeof(T)); */
            normv_file _f;
            _f._id = get_id();
            _f._maxm = _mmass;
            _f._minm = _nmass;
            _f._maxt = _mtime;
            _f._maxp = _mpos;
            _f._maxv = _mvel;
            out.write((char *) &_f, sizeof(normv_file));
            std::ofstream::pos_type _pos = out.tellp(); // will just be `sizeof(normv_file)`
            out.close();
            return _pos;
        }
        virtual normv_file load_normv(const char *path, bool adopt_values = true) {
            if (!path)
                throw std::invalid_argument{"Error: path to .normv file passed cannot be nullptr.\n"};
            struct stat buff{};
            if (stat(path, &buff) == -1)
                throw std::ios_base::failure{"Error: could not obtain file information.\n"};
            if (!S_ISREG(buff.st_mode))
                throw std::ios_base::failure{"Error: specified file is not a regular file.\n"};
            if (buff.st_size != sizeof(normv_file))
                throw std::ios_base::failure{"Error: invalid .normv format (at least for given templated type T).\n"};
            std::ifstream in{path, std::ios_base::in | std::ios_base::binary};
            if (!in.good())
                throw std::ios_base::failure{"Error: could not open file.\n"};
            normv_file _f;
            in.read((char *) &_f, sizeof(normv_file));
            in.close();
            if (adopt_values) {
                if (_f._maxm < 0 || _f._minm < 0 || _f._maxt < 0 || _f._maxp < 0 || _f._maxv < 0)
                    throw std::invalid_argument{"Error: all maximum values passed must be non-negative.\n"};
                if (_f._minm > _f._maxm)
                    throw std::invalid_argument{"Error: minimum mass cannot be greater than maximum mass.\n"};
                this->_mmass = _f._maxm;
                this->_nmass = _f._minm;
                this->_mtime = _f._maxt;
                this->_mpos = _f._maxp;
                this->_mvel = _f._maxv;
            }
            return _f;
        }
    };
    template <typename T> requires (std::is_floating_point_v<T>)
    class mass_normaliser : virtual public normaliser<T> {
        consteval uint64_t get_id() const noexcept override {
            return 0b00000001; // 1
        }
        using typename normaliser<T>::h_type;
        using typename normaliser<T>::e_type;
        h_type *normalise_header(h_type *_h) const override {
            /* Divides masses by the maximum mass. This normalises masses to be within [0,1].  */
            // _h->hd[1] = 'P'; // to differentiate between the .3bod and .3bodpp formats
            _h->masses[0] /= normaliser<T>::_mmass;
            _h->masses[1] /= normaliser<T>::_mmass;
            _h->masses[2] /= normaliser<T>::_mmass;
            return _h;
        }
        e_type *normalise_entry(e_type *_e) const override {
            return _e; // mass_normaliser performs no action on positions or velocities
        }
    public:
        using normaliser<T>::normaliser; // bring all parent ctors into scope
    };
    template <typename T> requires (std::is_floating_point_v<T>)
    class log_mass_normaliser : public virtual normaliser<T> {
        T _lmax;
        consteval uint64_t get_id() const noexcept override {
            return 0b00000010; // 2
        }
        using typename normaliser<T>::h_type;
        using typename normaliser<T>::e_type;
        h_type *normalise_header(h_type *_h) const override {
            /* Takes the logarithms of masses and divides them by the maximum absolute value of all logarithms of mass
             * values. This normalises masses to be within [-1,1]. */
            // _h->hd[1] = 'P';
            _h->masses[0] = std::log10(_h->masses[0])/_lmax;
            _h->masses[1] = std::log10(_h->masses[1])/_lmax;
            _h->masses[2] = std::log10(_h->masses[2])/_lmax;
            return _h;
        }
        e_type *normalise_entry(e_type *_e) const override {
            return _e; // log_mass_normaliser performs no action on positions or velocities
        }
    public:
        log_mass_normaliser(T max_mass, T min_mass, T max_time, T max_pcomp, T max_vcomp) :
        normaliser<T>{max_mass, min_mass, max_time, max_pcomp, max_vcomp} {
            T _logmax = std::log10(max_mass);
            T _logmin = std::log10(min_mass);
            _logmax = ABS(_logmax);
            _logmin = ABS(_logmin);
            _lmax = _logmax > _logmin ? _logmax : _logmin;
        }
        void max_mass(T _new_val) override {
            if (_new_val < 0)
                throw std::invalid_argument{"Error: max. mass val. must be non-negative.\n"};
            normaliser<T>::_mmass = _new_val;
            T _logmax = std::log10(_new_val);
            _logmax = ABS(_logmax);
            if (_logmax > _lmax)
                _lmax = _logmax;
        }
        void min_mass(T _new_val) override {
            if (_new_val < 0)
                throw std::invalid_argument{"Error: min. mass val. must be non-negative.\n"};
            normaliser<T>::_nmass = _new_val;
            T _logmax = std::log10(_new_val);
            _logmax = ABS(_logmax);
            if (_logmax > _lmax)
                _lmax = _logmax;
        }
        typename normaliser<T>::normv_file load_normv(const char *path, bool adopt_values = true) {
            if (!adopt_values)
                return normaliser<T>::load_normv(path, adopt_values);
            typename normaliser<T>::normv_file _f = normaliser<T>::load_normv(path, adopt_values);
            T _logmax = std::log10(normaliser<T>::_mmass);
            T _logmin = std::log10(normaliser<T>::_nmass);
            _logmax = ABS(_logmax);
            _logmin = ABS(_logmin);
            _lmax = _logmax > _logmin ? _logmax : _logmin;
            return _f;
        }
    };
    template <typename T> requires (std::is_floating_point_v<T>)
    class time_normaliser : public virtual normaliser<T> {
        consteval uint64_t get_id() const noexcept override {
            return 0b00000100; // 4
        }
        using typename normaliser<T>::h_type;
        using typename normaliser<T>::e_type;
        h_type *normalise_header(h_type *_h) const override {
            /* Divides the time per epoch by the largest total simulation time, such that any duration of time between
             * a pair of epochs in any .3bodpp file picked would fall within the range [0,1]. */
            // _h->hd[1] = 'P';
            _h->dt /= normaliser<T>::_mtime;
            return _h;
        }
        e_type *normalise_entry(e_type *_e) const override {
            return _e; // time_normaliser performs no action on positions or velocities
        }
    public:
        using normaliser<T>::normaliser;
    };
    template <typename T> requires (std::is_floating_point_v<T>)
    class pos_normaliser : public virtual normaliser<T> {
        consteval uint64_t get_id() const noexcept override {
            return 0b00001000; // 8
        }
        using typename normaliser<T>::h_type;
        using typename normaliser<T>::e_type;
        h_type *normalise_header(h_type *_h) const override {
            // _h->hd[1] = 'P';
            return _h;
        }
        e_type *normalise_entry(e_type *_e) const override {
            /* Divides all position components by the largest absolute value position component such that all normalised
             * position components fall within [-1,1]. */
            T *ptr = (T*) _e->positions;
            unsigned int counter = 9;
            while (counter --> 0)
                *ptr++ /= normaliser<T>::_mpos;
            return _e;
        }
    public:
        using normaliser<T>::normaliser;
    };
    template <typename T> requires (std::is_floating_point_v<T>)
    class vel_normaliser : public virtual normaliser<T> {
        consteval uint64_t get_id() const noexcept override {
            return 0b00010000; // 16
        }
        using typename normaliser<T>::h_type;
        using typename normaliser<T>::e_type;
        h_type *normalise_header(h_type *_h) const override {
            // _h->hd[1] = 'P';
            return _h;
        }
        e_type *normalise_entry(e_type *_e) const override {
            /* Divides all position components by the largest absolute value position component such that all normalised
             * position components fall within [-1,1]. */
            T *ptr = (T*) _e->velocities;
            unsigned int counter = 9;
            while (counter --> 0)
                *ptr++ /= normaliser<T>::_mvel;
            return _e;
        }
    public:
        using normaliser<T>::normaliser;
    };
    /* template <typename T> requires (std::is_floating_point_v<T>)
    class log_vel_normaliser : public virtual normaliser<T> {

    }; */
    template <typename T, typename ...Args> requires (std::is_floating_point_v<T> && sizeof...(Args) > 0 &&
                                                     (std::is_base_of_v<normaliser<T>, Args> && ...))
    class gen_normaliser : public Args... {
        consteval uint64_t get_id() const noexcept override {
            return (Args::get_id() | ...); // take bitwise OR of all IDs of parent classes
        }
        using typename normaliser<T>::h_type;
        using typename normaliser<T>::e_type;
        h_type *normalise_header(h_type *_h) const override {
            return (Args::normalise_header(_h), ...);
        }
        e_type *normalise_entry(e_type *_e) const override {
            return (Args::normalise_entry(_e), ...);
        }
    public:
        using normaliser<T>::normaliser;
    };
    template <typename T> requires (std::is_floating_point_v<T>)
    std::pair<std::string, std::vector<std::string>> preprocess(const char *dpath, const normaliser<T> &_norm) {
        /* Preprocesses the entire dataset. Goes through all .3bod files in given directory, finds max. values, uses
         * these to normalise the data in each .3bod file, placing the normalised data into each corresponding .3bodpp
         * file in a newly created directory. */
        if (!dpath)
            throw std::invalid_argument{"Error: path to directory cannot be nullptr.\n"};
        DIR *dir;
        struct dirent *entry;
        if (!(dir = opendir(dpath))) {
            std::cerr << "Error: could not open path to raw data directory.\n";
            exit(1);
        }
        std::vector<std::string> files;
        T max_mass = 0;
        T min_mass = std::numeric_limits<T>::infinity();
        long double max_sim_time = 0;
        long double sim_time;
        T max_abs_pcomp = 0;
        T max_abs_vcomp = 0;
        T pv_val;
        std::string fpath = dpath;
        if (!fpath.ends_with('/'))
            fpath.push_back('/');
        uint64_t dpath_len = fpath.size();
        typename f3bodr<T>::hdr_t *header;
        typename f3bodr<T>::iterator_type _it;
        typename f3bodr<T>::iterator_type _end;
        uint64_t counter;
        T *ptr;
        while ((entry = readdir(dir))) {
            if (gtd::endswith(entry->d_name, ".3bod")) {
                fpath += entry->d_name;
                f3bodr<T> reader{fpath.c_str()}; // might throw due to `T` mismatch, catch outside
                header = &reader.header();
                if (!header->N)
                    continue; // empty .3bod files do not get considered
                counter = 0;
                while (counter < 3) {
                    if (header->masses[counter] > max_mass)
                        max_mass = header->masses[counter];
                    if (header->masses[counter] < min_mass)
                        min_mass = header->masses[counter];
                    ++counter;
                }
                if ((sim_time = header->dt*(header->N - 1)) > max_sim_time)
                    max_sim_time = sim_time;
                _it = reader.begin();
                _end = reader.end();
                while (_it != _end) {
                    ptr = (T*) _it->positions;
                    counter = 0;
                    while (counter++ < 9) {
                        if ((pv_val = ABS(*ptr)) > max_abs_pcomp)
                            max_abs_pcomp = pv_val;
                        ++ptr;
                    }
                    counter = 0;
                    while (counter++ < 9) {
                        if ((pv_val = ABS(*ptr)) > max_abs_vcomp)
                            max_abs_vcomp = pv_val;
                        ++ptr;
                    }
                    ++_it;
                }
                // files.emplace_back(entry->d_name);
                files.push_back(fpath); // I prefer the slightly increased memory use in storing the dir. name prepended
                fpath.erase(dpath_len); // before each filename, rather than having to re-prepend it later
            }
        }
        closedir(dir);
        if (files.empty()) {
            std::cerr << "Error: no .3bod files found within directory \"" << dpath << "\"\n";
            exit(1);
        }
        fpath.pop_back(); // remove the trailing slash
        fpath += "_preprocessed_";
        const char *date_and_time = gtd::get_date_and_time(); // not thread-safe
        fpath += date_and_time; // not thread-safe
        if (mkdir(fpath.c_str(), S_IRWXU | S_IRWXG | S_IRWXO) == -1) {
            std::cerr << "Error: could not create \"" << fpath << "\" directory.\n";
            exit(1);
        }
        fpath.push_back('/');
        uint64_t _new_dpath_len = dpath_len + 14 + 25 + 1;
        std::vector<std::string> _new_files;
        _new_files.reserve(files.size());
        _norm.max_mass(max_mass);
        _norm.min_mass(min_mass);
        _norm.max_time(max_sim_time);
        _norm.max_pos(max_abs_pcomp);
        _norm.max_vel(max_abs_vcomp);
        for (const std::string &_old_path : files) {
            fpath.append(_old_path.begin() + dpath_len, _old_path.end());
            _norm(_old_path, fpath);
            _new_files.push_back(fpath);
            fpath.erase(_new_dpath_len);
        }
        // Write max. values to binary data file
        fpath += "norm_vals_";
        fpath += date_and_time;
        _norm.write_normv((fpath += ".normv").c_str());
        return {fpath, _new_files};
    }
}
#endif
