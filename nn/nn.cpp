#include "../glib/ml/gregffnn.hpp"
#include "../glib/misc/gregparse.hpp"
#include "nnsup.hpp"
#include <algorithm>

#define DEF_NUM_PASSES ((uint64_t) 1'000)

[[noreturn]] void pp_error() {
    std::cerr << "Error: only one of \"--just_preprocess\", \"--preprocess\" or \"--not_preprocessed\" can be "
                 "specified.\n";
    exit(1);
}

template <typename dtype, typename wtype> requires (std::is_floating_point_v<dtype> && std::is_floating_point_v<wtype>)
void run_sim() {
    /* Runs the entire simulation. `dtype` is the data type of the data and `wtype` is the type to use in the NN. */
}

template <typename T> requires (std::is_floating_point_v<T>)
gtd::normaliser<T> *get_normaliser(const char *norm_arg) {
    // logic...
    if (!norm_arg)
        return nullptr;
    norm_arg += 7;
    if (!*norm_arg)
        return nullptr;
    bool _m = false;
    bool _l = false;
    bool _t = false;
    bool _p = false;
    bool _v = false;
    do {
        if (*norm_arg == 'm')
            _m = true;
        if (*norm_arg == 'l')
            _l = true;
        if (*norm_arg == 't')
            _t = true;
        if (*norm_arg == 'p')
            _p = true;
        if (*norm_arg == 'v')
            _v = true;
    } while (*++norm_arg);
    if (_m) {
        if (_t) {
            if (_p) {
                if (_v)
                    return new gtd::gen_normaliser<T, gtd::mass_normaliser<T>, gtd::time_normaliser<T>, gtd::pos_normaliser<T>, gtd::vel_normaliser<T>>{};
                else
                    return new gtd::gen_normaliser<T, gtd::mass_normaliser<T>, gtd::time_normaliser<T>, gtd::pos_normaliser<T>>{};
            }
            else {
                if (_v)
                    return new gtd::gen_normaliser<T, gtd::mass_normaliser<T>, gtd::time_normaliser<T>, gtd::vel_normaliser<T>>{};
                else
                    return new gtd::gen_normaliser<T, gtd::mass_normaliser<T>, gtd::time_normaliser<T>>{};
            }
        }
        else {
            if (_p) {
                if (_v)
                    return new gtd::gen_normaliser<T, gtd::mass_normaliser<T>, gtd::pos_normaliser<T>, gtd::vel_normaliser<T>>{};
                else
                    return new gtd::gen_normaliser<T, gtd::mass_normaliser<T>, gtd::pos_normaliser<T>>{};
            }
            else {
                if (_v)
                    return new gtd::gen_normaliser<T, gtd::mass_normaliser<T>, gtd::vel_normaliser<T>>{};
                else
                    return new gtd::gen_normaliser<T, gtd::mass_normaliser<T>>{};
            }
        }
    } else {
        if (_l) { // `_m` and `_l` are mutually exclusive (this was already checked by the regex)
            if (_t) {
                if (_p) {
                    if (_v)
                        return new gtd::gen_normaliser<T, gtd::log_mass_normaliser<T>, gtd::time_normaliser<T>, gtd::pos_normaliser<T>, gtd::vel_normaliser<T>>{};
                    else
                        return new gtd::gen_normaliser<T, gtd::log_mass_normaliser<T>, gtd::time_normaliser<T>, gtd::pos_normaliser<T>>{};
                }
                else {
                    if (_v)
                        return new gtd::gen_normaliser<T, gtd::log_mass_normaliser<T>, gtd::time_normaliser<T>, gtd::vel_normaliser<T>>{};
                    else
                        return new gtd::gen_normaliser<T, gtd::log_mass_normaliser<T>, gtd::time_normaliser<T>>{};
                }
            }
            else {
                if (_p) {
                    if (_v)
                        return new gtd::gen_normaliser<T, gtd::log_mass_normaliser<T>, gtd::pos_normaliser<T>, gtd::vel_normaliser<T>>{};
                    else
                        return new gtd::gen_normaliser<T, gtd::log_mass_normaliser<T>, gtd::pos_normaliser<T>>{};
                }
                else {
                    if (_v)
                        return new gtd::gen_normaliser<T, gtd::log_mass_normaliser<T>, gtd::vel_normaliser<T>>{};
                    else
                        return new gtd::gen_normaliser<T, gtd::log_mass_normaliser<T>>{};
                }
            }
        }
        else {
            if (_t) {
                if (_p) {
                    if (_v)
                        return new gtd::gen_normaliser<T, gtd::time_normaliser<T>, gtd::pos_normaliser<T>, gtd::vel_normaliser<T>>{};
                    else
                        return new gtd::gen_normaliser<T, gtd::time_normaliser<T>, gtd::pos_normaliser<T>>{};
                }
                else {
                    if (_v)
                        return new gtd::gen_normaliser<T, gtd::time_normaliser<T>, gtd::vel_normaliser<T>>{};
                    else
                        return new gtd::gen_normaliser<T, gtd::time_normaliser<T>>{};
                }
            }
            else {
                if (_p) {
                    if (_v)
                        return new gtd::gen_normaliser<T, gtd::pos_normaliser<T>, gtd::vel_normaliser<T>>{};
                    else
                        return new gtd::gen_normaliser<T, gtd::pos_normaliser<T>>{};
                }
                else {
                    // if (_v)
                        return new gtd::gen_normaliser<T, gtd::vel_normaliser<T>>{};
                    // else
                    //     return new gtd::gen_normaliser<T>{}; // this is not an option, regex anyway assures this
                }
            }
        }
    }
    // return new gtd::mass_normaliser<T>{}; // placeholder until logic above is implemented
}

int main(int argc, char **argv) {
    gtd::parser parser{argc, argv};
    const char *data_dir = parser.get_arg("--data_dir");
    if (!data_dir)
        data_dir = ".";
    /* There are 4 possibilities for running this program, selection of which is controlled by the 3 booleans below.
     * If:
     *     - All 3 booleans are false, then the program interprets the data in `data_dir` as being preprocessed and runs
     *       the neural net directly on these (thus looks for .3bodpp files)
     *     - `jpp` is true, then the program just performs preprocessing on the allegedly unprocessed data in `data_dir`
     *       and saves it all to a new folder
     *     - `pp` is true, then the allegedly unprocessed data in `data_dir` is preprocessed and written all to a new
     *       folder, and the neural network is run on this freshly preprocessed data
     *     - `npp` is true, then the neural network runs directly on the allegedly unprocessed data in `data_dir` (and
     *       thus looks for .3bod files)
     * Any other combination is invalid. */
    bool jpp = parser.get_arg("--just_preprocess", false); // only preprocess data in `data_dir`, don't run neural net
    bool pp = parser.get_arg("--preprocess", false); // preprocess data in `data_dir` and run neural net
    bool npp = parser.get_arg("--not_preprocessed", false); // data in `data_dir` isn't preprocessed and shouldn't be
    const char *norm_str = parser.get_arg(std::regex{R"(^--norm=(m|l)?t?p?v?$)"}); // improve this regex
    if (jpp) {
        if (pp || npp)
            pp_error();
        try {
            gtd::normaliser<long double> *nptr = get_normaliser<long double>(norm_str);
            gtd::preprocess(data_dir, *nptr);
            delete nptr;
        } catch (const std::exception&) {
            try {
                gtd::normaliser<double> *nptr = get_normaliser<double>(norm_str);
                gtd::preprocess(data_dir, *nptr);
                delete nptr;
            } catch (const std::exception&) {
                try {
                    gtd::normaliser<float> *nptr = get_normaliser<float>(norm_str);
                    gtd::preprocess(data_dir, *nptr);
                    delete nptr;
                } catch (const std::exception&) {
                    std::cerr << "Error: unrecognised floating point data type in data files.\n";
                    return 1;
                }
            }
        }
        return 0;
    }
    if (pp) {
        if (npp)
            pp_error();
    }
    if (npp) {

    }
    const char *layers = parser.get_arg(std::regex{R"(^--layers=\d+-\d+:\d+(,\d+-\d+:\d+)*$)"});
    if (!layers)
        layers = "--layers=22-128:1,128-64:1,64-18:0";
    const char *ftype = parser.get_arg("--fp_type");
    if (!ftype)
        ftype = "long double";
    /* uint64_t num_passes = parser.get_arg("--num_passes", DEF_NUM_PASSES);
    if (!num_passes) {
        std::cerr << "Error: number of passes over entire dataset cannot be zero.\n";
        return 1;
    }
    std::vector<std::string> files;
    struct dirent *entry;
    DIR *dir;
    if (!(dir = opendir(data_dir))) {
        std::cerr << "Error: could not open directory \"" << data_dir << "\".\n";
        return 1;
    }
    while ((entry = readdir(dir)))
        if (gtd::endswith(entry->d_name, ".3bod"))
            files.emplace_back(entry->d_name);
    closedir(dir);
    if (files.empty()) {
        std::cerr << "Error: no .3bod files found in specified directory \"" << data_dir << "\".\n";
        return 1;
    }
    uint64_t pass = 0;
    std::mt19937_64 rng{std::random_device{}()};
    while (pass++ < num_passes) {
        std::shuffle(files.begin(), files.end(), rng); // shuffle .3bod files so pass over them isn't in the same order
    } */
    return 0;
}
