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
        T _mmass{};
        T _mtime{};
        T _mpos{};
        T _mvel{};
        virtual h_type normalise_header(const h_type&) = 0;
        virtual e_type normalise_entry(const e_type&) = 0;
    public:
        normaliser() = default;
        normaliser(T max_mass, T max_time, T max_abs_poscomp, T max_abs_velcomp) : _mmass{max_mass}, _mtime{max_time},
        _mpos{max_abs_poscomp}, _mvel{max_abs_velcomp} {
            if (_mmass < 0 || _mtime < 0 || _mpos < 0 || _mvel < 0)
                throw std::invalid_argument{"Error: all maximum values passed must be non-negative.\n"};
        }
        T max_mass() const noexcept {
            return _mmass;
        }
        T max_time() const noexcept {
            return _mtime;
        }
        T max_pos() const noexcept {
            return _mpos;
        }
        T max_vel() const noexcept {
            return _mvel;
        }
        void max_mass(T _new_val) {
            if (_new_val < 0)
                throw std::invalid_argument{"Error: max. mass val. must be non-negative.\n"};
            _mmass = _new_val;
        }
        void max_time(T _new_val) {
            if (_new_val < 0)
                throw std::invalid_argument{"Error: max. time val. must be non-negative.\n"};
            _mtime = _new_val;
        }
        void max_pos(T _new_val) {
            if (_new_val < 0)
                throw std::invalid_argument{"Error: max. pos. val. must be non-negative.\n"};
            _mpos = _new_val;
        }
        void max_vel(T _new_val) {
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
            // const h_type &_current_header = reader.header();
            h_type _new_header = normalise_header(reader.header());
            out.write((char *) _new_header, sizeof(h_type));
            e_type _new_entry;
            // e_type _new_entries = new e_type[_current_header.N];
            // e_type *ptr = _new_entries;
            for (const e_type &entry : reader) {
                // *ptr++ = normalise_entry(entry);
                _new_entry = normalise_entry(entry);
                out.write((char *) _new_entry, sizeof(e_type));
            }
            // out.write((char *) _new_entries, _current_header.N*sizeof(e_type));
            // delete [] _new_entries;
            out.close();
        }
        void write_max_vals(const char *path) const {
            // write all max. vals followed by a hash of some sort indicating the norm. method (abstract method)
        }
    };
    template <typename T> requires (std::is_floating_point_v<T>)
    class mass_log_normaliser : virtual public normaliser<T> {

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
        fpath += gtd::get_date_and_time(); // not thread-safe
        if (mkdir(fpath.c_str(), S_IRWXU | S_IRWXG | S_IRWXO) == -1) {
            std::cerr << "Error: could not create \"" << fpath << "\" directory.\n";
            exit(1);
        }
        fpath.push_back('/');
        uint64_t _new_dpath_len = dpath_len + 14 + 25 + 1;
        std::vector<std::string> _new_files;
        _new_files.reserve(files.size());
        _norm.max_mass(max_mass);
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
    }
}
#endif
