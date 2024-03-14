#ifndef NNSUP_HPP
#define NNSUP_HPP

#include "../data_management/datsup.hpp"
#include <dirent.h>
#include <memory>

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
        constexpr virtual uint64_t get_id() const noexcept {return 0;}
    public:
#pragma pack(push, 1)
        struct normv_file {
            const char header[4] = {'N', 'O', 'R', 'M'};
            const uint32_t _sizeof_T = sizeof(T);
            const uint32_t _sizeof_ld = sizeof(long double);
            // const uint32_t _sizeof_u64 = sizeof(uint64_t); // always 8
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
                throw std::invalid_argument{"Error: all maximum/minimum values passed must be non-negative.\n"};
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
            if (_new_val < this->_nmass)
                throw std::invalid_argument{"Error: max. mass cannot be less than current min. mass. If you wish to "
                                            "reduce both, please set min. mass first.\n"};
            _mmass = _new_val;
        }
        virtual void min_mass(T _new_val) {
            if (_new_val < 0)
                throw std::invalid_argument{"Error: min. mass val. must be non-negative.\n"};
            if (_new_val > this->_mmass)
                throw std::invalid_argument{"Error: min. mass cannot be greater than current max. mass. If you wish to "
                                            "increase both, please set max. mass first.\n"};
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
        virtual std::ofstream::pos_type operator()(const std::string &input, const std::string &output) {
            f3bodr<T> reader{input.c_str()};
            std::ofstream out{output, std::ios_base::out | std::ios_base::binary};
            if (!out.good())
                throw std::ios_base::failure{"Error: could not open output .3bodpp file.\n"};
            reader.hdr.hd[1] = 'P';
            out.write((char *) normalise_header(&reader.hdr), sizeof(h_type));
            e_type *_entries = reader._it._arr;
            e_type *_entry = _entries;
            uint64_t counter = reader.hdr.N;
            while (counter --> 0)
                normalise_entry(_entry++);
            out.write((char *) _entries, reader.hdr.N*sizeof(e_type));
            std::ofstream::pos_type _pos = out.tellp();
            out.close();
            return _pos;
        }
        consteval size_t normv_fsize() const noexcept { // returns the size of a .normv file for the given `T` type
            return sizeof(normv_file);
        }
        std::ofstream::pos_type write_normv(const char *path) const {
            /* Writes the maximum values to a .normv file. The file format of a .normv file is as follows:
             *     - 4 bytes: Characters 'N', 'O', 'R' and 'M'
             *     - 4 bytes: `sizeof(T)`
             *     - 4 bytes: `sizeof(long double)`
             *     - `sizeof(uint64_t)` bytes: ID of normaliser
             *     - `4*sizeof(T) + sizeof(long double)` bytes: Maximum/minimum values, in same order of class */
            if (!path)
                throw std::invalid_argument{"Error: path for .normv file cannot be nullptr.\n"};
            std::ofstream out{path, std::ios_base::out | std::ios_base::binary};
            if (!out.good())
                throw std::ios_base::failure{"Error: could not open .normv file for writing.\n"};
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
                    throw std::invalid_argument{"Error: all max./min. values must be non-negative.\n"};
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
        virtual ~normaliser() = default; // required for deleting derived class object through base class pointer
    };
    template <typename T> requires (std::is_floating_point_v<T>)
    class mass_normaliser : virtual public normaliser<T> {
    protected:
        constexpr uint64_t get_id() const noexcept override {
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
    protected:
        T _lmax;
        constexpr uint64_t get_id() const noexcept override {
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
        log_mass_normaliser() = default;
        log_mass_normaliser(T max_mass, T min_mass, T max_time, T max_pcomp, T max_vcomp) :
        normaliser<T>{max_mass, min_mass, max_time, max_pcomp, max_vcomp} {
            T _logmax = std::log10(max_mass);
            T _logmin = std::log10(min_mass);
            _logmax = ABS(_logmax);
            _logmin = ABS(_logmin);
            _lmax = _logmax > _logmin ? _logmax : _logmin;
        }
        void max_mass(T _new_val) override {
            normaliser<T>::max_mass(_new_val); // takes care of exceptions
            T _logmax = std::log10(_new_val);
            _logmax = ABS(_logmax);
            if (_logmax > _lmax)
                _lmax = _logmax;
        }
        void min_mass(T _new_val) override {
            normaliser<T>::min_mass(_new_val);
            T _logmin = std::log10(_new_val);
            _logmin = ABS(_logmin);
            if (_logmin > _lmax)
                _lmax = _logmin;
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
    protected:
        constexpr uint64_t get_id() const noexcept override {
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
    protected:
        constexpr uint64_t get_id() const noexcept override {
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
    protected:
        constexpr uint64_t get_id() const noexcept override {
            return 0b00010000; // 16
        }
        using typename normaliser<T>::h_type;
        using typename normaliser<T>::e_type;
        h_type *normalise_header(h_type *_h) const override {
            // _h->hd[1] = 'P';
            return _h;
        }
        e_type *normalise_entry(e_type *_e) const override {
            /* Divides all position components by the largest absolute value velocity component such that all normalised
             * velocity components fall within [-1,1]. */
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
    protected:
        constexpr uint64_t get_id() const noexcept override {
            return (Args::get_id() | ...); // take bitwise OR of all IDs of parent classes
        }
        using h_type = typename f3bodr<T>::hdr_t;
        using e_type = typename f3bodr<T>::entry_type;
        h_type *normalise_header(h_type *_h) const override {
            return (Args::normalise_header(_h), ...);
        }
        e_type *normalise_entry(e_type *_e) const override {
            return (Args::normalise_entry(_e), ...);
        }
    public:
        // using normaliser<T>::normaliser;
        gen_normaliser() = default;
    };
    template <typename T> requires (std::is_floating_point_v<T>)
    void preprocess(const char *dpath, // path to training data directory
                    const char *vdir, // path to validation data directory (can be `nullptr`)
                    // char **normv_path_d, // made to point to preprocessed training data directory .normv file
                    // char **normv_path_v, // made to point to preprocessed validation data directory .normv file
                    std::unique_ptr<std::vector<std::string>> &dppf, // preprocessed training data file paths
                    std::unique_ptr<std::vector<std::string>> &vppf, // preprocessed validation data file paths
                    normaliser<T> &_norm, // normaliser
                    bool jpp /* whether to just preprocess */) {
        /* Preprocesses the entire dataset. Goes through all .3bod files in given directory, finds max. values, uses
         * these to normalise the data in each .3bod file, placing the normalised data into each corresponding .3bodpp
         * file in a newly created directory. */
        if (!dpath)
            throw std::invalid_argument{"Error: path to training data directory cannot be nullptr.\n"};
        // if (!jpp) {
        //     if (!dppf)
        //         throw std::invalid_argument{"Error: pointer to pointer to std::vector<std::string> for holding "
        //                                     "preprocessed training data file paths cannot be nullptr unless data is "
        //                                     "only being preprocessed.\n"};
        //     if (vdir && !vppf)
        //         throw std::invalid_argument{"Error: pointer to pointer to std::vector<std::string> for holding "
        //                                     "preprocessed validation data file paths cannot be nullptr if vdir is not "
        //                                     "nullptr, unless data is only being preprocessed.\n"};
        // }
        DIR *dir;
        struct dirent *entry;
        std::string fpath;
        if (!(dir = opendir(dpath))) {
            fpath = "Error: could not open path to training data directory \"";
            fpath += dpath;
            fpath += "\".\n";
            throw std::ios_base::failure{fpath};
        }
        fpath = dpath;
        std::vector<std::string> files;
        T max_mass = 0;
        T min_mass = std::numeric_limits<T>::infinity();
        long double max_sim_time = 0;
        long double sim_time;
        T max_abs_pcomp = 0;
        T max_abs_vcomp = 0;
        T pv_val;
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
                header = &reader.hdr;
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
            fpath = "Error: no .3bod files found within directory \"";
            fpath += dpath;
            fpath += "\".\n";
            throw std::logic_error{fpath};
        }
        if (!str_eq(dpath, ".")) {
            fpath.pop_back(); // remove the trailing slash
            fpath += "_preprocessed_";
        } else {
            fpath += "preprocessed_";
        }
        const char *date_and_time = gtd::get_date_and_time(); // not thread-safe
        fpath += date_and_time; // not thread-safe
        if (mkdir(fpath.c_str(), S_IRWXU | S_IRWXG | S_IRWXO) == -1) {
            files[0] = "Error: could not create \"";
            files[0] += fpath;
            files[0] += "\" directory.\n";
            throw std::ios_base::failure{files[0]};
        }
        _norm.max_mass(max_mass);
        _norm.min_mass(min_mass);
        _norm.max_time(max_sim_time);
        _norm.max_pos(max_abs_pcomp);
        _norm.max_vel(max_abs_vcomp);
        uint64_t _prev_len = fpath.size();
        if (jpp) { // case for just preprocessing - so no need to return preprocessed filenames
            if (vdir) {
                std::vector<std::string> *vfiles = find_files(vdir, ".3bod", true);
                if (!vfiles) {
                    fpath = "Error: no .3bod files found in validation data directory \"";
                    fpath += vdir;
                    fpath += "\".\n";
                    throw std::logic_error{fpath};
                }
                fpath += "_validation/";
                uint64_t _vdir_len = strlen_c(vdir);
                if (*(vdir + _vdir_len - 1) != '/')
                    ++_vdir_len;
                uint64_t _new_len = _prev_len + 12;
                for (const std::string &vfile : *vfiles) {
                    fpath.append(vfile.begin() + _vdir_len, vfile.end());
                    fpath += "pp";
                    _norm(vfile, fpath);
                    fpath.erase(_new_len);
                }
                delete [] vfiles;
                std::cout << "Preprocessed validation data files written to \"" << fpath << "\".\n";
                fpath += "norm_vals_for_valset_";
                fpath += date_and_time;
                _norm.write_normv((fpath += ".normv").c_str());
                std::cout << "Normalisation values within validation set written to \"" << fpath << "\".\n";
                // if (normv_path_v)
                //     *normv_path_v = strcpy_c(new char[fpath.size() + 1], fpath.c_str());
                fpath.erase(_prev_len);
            }
            fpath.push_back('/');
            uint64_t _new_dpath_len = _prev_len + 1; /* dpath_len + 14 + 25 + 2; */
            for (const std::string &_old_path : files) {
                fpath.append(_old_path.begin() + dpath_len, _old_path.end());
                fpath += "pp"; // to make the extension .3bodpp ("pp" for "preprocessed")
                _norm(_old_path, fpath);
                fpath.erase(_new_dpath_len);
            }
            std::cout << "Preprocessed training data files written to \"" << fpath << "\".\n";
            // Write max. values to binary data file
            fpath += "norm_vals_";
            fpath += date_and_time;
            _norm.write_normv((fpath += ".normv").c_str());
            std::cout << "Normalisation values within training set written to \"" << fpath << "\".\n";
            // if (normv_path_d)
            //     *normv_path_d = strcpy_c(new char[fpath.size() + 1], fpath.c_str());
            // *dppf = nullptr;
            // if (vppf)
            //     *vppf = nullptr;
            // return {fpath, nullptr}; // nullptr for `std::vector` as the names of the preprocessed files are not needed
            return;
        }
        if (vdir) {
            std::vector<std::string> *vfiles = find_files(vdir, ".3bod", true);
            if (!vfiles) {
                fpath = "Error: no .3bod files found in validation data directory \"";
                fpath += vdir;
                fpath += "\".\n";
                throw std::logic_error{fpath};
            }
            fpath += "_validation/";
            uint64_t _vdir_len = strlen_c(vdir);
            if (*(vdir + _vdir_len - 1) != '/')
                ++_vdir_len;
            uint64_t _new_len = _prev_len + 12;
            std::vector<std::string> *_new_vfiles = new std::vector<std::string>{};
            _new_vfiles->reserve(vfiles->size());
            for (const std::string &vfile : *vfiles) {
                fpath.append(vfile.begin() + _vdir_len, vfile.end());
                fpath += "pp";
                _norm(vfile, fpath);
                _new_vfiles->push_back(fpath);
                fpath.erase(_new_len);
            }
            delete [] vfiles;
            std::cout << "Preprocessed validation data files written to \"" << fpath << "\".\n";
            fpath += "norm_vals_for_valset_";
            fpath += date_and_time;
            _norm.write_normv((fpath += ".normv").c_str());
            std::cout << "Normalisation values within validation set written to \"" << fpath << "\".\n";
            // if (normv_path_v)
            //     *normv_path_v = strcpy_c(new char[fpath.size() + 1], fpath.c_str());
            fpath.erase(_prev_len);
            // *vppf = _new_vfiles;
            vppf.reset(_new_vfiles);
        }
        fpath.push_back('/');
        uint64_t _new_dpath_len = fpath.size(); /* dpath_len + 14 + 25 + 2; */
        std::vector<std::string> *_new_files = new std::vector<std::string>{};
        _new_files->reserve(files.size());
        for (const std::string &_old_path : files) {
            fpath.append(_old_path.begin() + dpath_len, _old_path.end());
            fpath += "pp";
            _norm(_old_path, fpath);
            _new_files->push_back(fpath);
            fpath.erase(_new_dpath_len);
        }
        std::cout << "Preprocessed training data files written to \"" << fpath << "\".\n";
        fpath += "norm_vals_";
        fpath += date_and_time;
        _norm.write_normv((fpath += ".normv").c_str());
        std::cout << "Normalisation values within training set written to \"" << fpath << "\".\n";
        // if (normv_path_d)
        //     *normv_path_d = strcpy_c(new char[fpath.size() + 1], fpath.c_str());
        // *dppf = _new_files;
        dppf.reset(_new_files);
        // return {fpath, _new_files};
    }
}
#endif
